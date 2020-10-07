#' Title
#'
#' @param xs0 xs0
#' @param y0 y0
#' @param lam lam
#' @param Rs Rs
#' @param NCORES number of cores
#'
#' @return
#' @export
#'
#' @useDynLib ddsPLS
do_loo <- function(xs0,y0,lam=0,Rs=1,NCORES=1,method=2){
  n <- nrow(xs0[[1]])
  q <- ncol(y0)
  p <- ncol(xs0[[1]])
  xs0_cbind <- list(do.call(cbind,xs0))
  PRESS_y <- matrix(NA,n,q)
  B <- list()
  K <- length(xs0)
  for(i in 1:n){
    X_train <- xs0_cbind[[1]][-i,,drop=F]
    X_test <- xs0_cbind[[1]][i,,drop=F]
    Y_train <- y0[-i,,drop=FALSE]
    Y_test <- y0[i,,drop=FALSE]
    m_1 <- ddsPLS2(Xs = list(X_train),Y = Y_train,lam = lam,Rs = Rs,method=method)
    B_i <- m_1$B[[1]]
    mu_k <- matrix(unlist(m_1$parameters$x$mu),nrow = 1)
    sd_k <- matrix(unlist(m_1$parameters$x$sd),nrow = 1)
    mu_y <- matrix(unlist(m_1$parameters$y$mu),nrow = 1)
    sd_y <- matrix(unlist(m_1$parameters$y$sd),nrow = 1)
    y_pred <- mu_y + (((X_test+mu_k)/sd_k)%*%B_i)*sd_y
    PRESS_y[i,] <- (y_pred-Y_test)^2
    B[[i]] <- B_i
  }
  list(id=1:n,PRESS_y=PRESS_y,B=B)
}


#' Title
#'
#' @param Xs Xs
#' @param Y Y
#' @param lambdas lambdas
#' @param tau rate learning
#' @param NCORES number of cores
#' @param NZV near zero var
#'
#' @return
#' @export
#'
#' @useDynLib ddsPLS
Q2_local_ddsPLS <- function(Xs,Y,lambdas = 0.5,
                            tau=0.0975,NCORES=1,
                            NZV=1e-3,method=2){
  n <- nrow(Xs[[1]])
  q <- ncol(Y)
  id_na <- lapply(Xs,function(xx){which(is.na(xx[,1]))})
  K <- length(Xs)
  Ls <- length(lambdas)
  ncomps = Q2_TOTAL <- rep(NA,Ls)
  q <- ncol(Y)
  Y_init <- scale(Y)
  Xs_init <- lapply(Xs,scale)
  RSS0 = RSS_h_moins_1 <- colSums(Y_init^2)
  PRESS_y = RSS_y = Q2_y <- matrix(NA,n,q)
  lambda <- rep(NA,n)
  ps <- unlist(lapply(Xs_init,ncol))
  xs0_cbind <- list(do.call(cbind,Xs_init))
  u <- matrix(NA,sum(ps),n)
  P <- matrix(0,sum(ps),n)
  B <- matrix(0,sum(ps),q)
  Us = Bs <- list()
  B_tot_LOO <- matrix(0,n,sum(ps)*q)

  ### For each 'h' component, look for the best lambda
  test <- T
  y0 <- Y_init
  h <- 0
  V <- NULL
  K_h <- 1:K
  while(test){
    h <- h + 1
    Bis <- list()
    # for(i_l in 1:length(lambdas)){
    NCORES_w <- min(NCORES,Ls)
    `%my_do%` <- ifelse(NCORES_w!=1,{
      out<-`%dopar%`
      cl <- makeCluster(NCORES_w)#cl <- parallel::makeCluster(NCORES_w)
      registerDoParallel(cl)#doParallel::registerDoParallel(cl)
      out},{
        out <- `%do%`
        out})
    pos_decoupe <- 1
    PRESS_il_y_and_i_l <- foreach(i_l=1:Ls,.combine = "rbind", .packages = c("Rcpp","ddsPLS"),.noexport ="getCOV_COV") %my_do% {
      lam <- lambdas[i_l]
      out <- do_loo(xs0_cbind,y0,lam=lam)
      c(i_l,colSums(out$PRESS_y))
    }
    if(NCORES_w!=1){
      stopCluster(cl)
    }
    if(Ls>1){
      PRESS_il_y <- PRESS_il_y_and_i_l[order(PRESS_il_y_and_i_l[,1,drop=F]),-1]
    }else{
      PRESS_il_y <- matrix(PRESS_il_y_and_i_l[-1],nrow=1)
    }
    RSS_il_y <- matrix(NA,Ls,q)
    for(i_l in 1:Ls){
      lam <- lambdas[i_l]
      m_1 <- ddsPLS2(Xs = xs0_cbind,Y = y0,lam = lam,Rs = 1,method=method)
      y_est_train <- m_1$y_est
      RSS_il_y[i_l,] <- colSums((y_est_train-y0)^2)
    }
    Q2_h_y <- 1-PRESS_il_y/matrix(rep(RSS_h_moins_1,Ls),nrow = Ls,byrow = T)
    Q2_h_sum <- 1-rowSums(PRESS_il_y)/sum(RSS_h_moins_1)
    max_Q2_y <- apply(Q2_h_y,MARGIN = 1,max)
    best_id_h <- which.max(max_Q2_y)
    if(length(best_id_h)>0){
      test_h <- max_Q2_y[best_id_h]>tau
      best_lam_h <- lambdas[best_id_h]
      lambda[h] <- best_lam_h
      RSS_y[h,] <- RSS_il_y[best_id_h,]
      RSS_h_moins_1 <- RSS_y[h,]
      PRESS_y[h,] <- PRESS_il_y[best_id_h,]
      Q2_y[h,] <- Q2_h_y[best_id_h,]
      if(!test_h){
        test <- F
      }else{
        # Get the regression matrix of the optimal model
        m_plus <- model_PLS(xs0_cbind[[1]],y0,best_lam_h,R = 1, tau=tau,method=method,NZV=NZV)
        if(norm(m_plus$U_out)>NZV){
          u[,h] <- m_plus$U_out
          t_r <- xs0_cbind[[1]]%*%u[,h,drop=F]
          P[,h] <- crossprod(xs0_cbind[[1]],t_r)/sum(t_r^2)
          if(h!=1){
            for(s_r in (h-1):1){
              u[,h] <- u[,h]-u[,s_r]*sum(P[,s_r]*u[,h])
            }
          }
          if(sum(u[,h]^2)>NZV){
            B_r <- tcrossprod(u[,h,drop=F],m_plus$V_out)
            V_r <- m_plus$V_out
            y0_plus <- xs0_cbind[[1]]%*%B_r
            Var_rel <- sum(y0_plus^2)/sum(RSS0)
            print(Var_rel)
            if(Var_rel<NZV){
              test <- F
            }else{
              B <- B + B_r
              V <- cbind(V,V_r)
              y0 <- y0 - y0_plus
            }
          }else{
            test <- F
          }
        }else{
          test <- F
        }
      }
    }
  }
  B <- B * matrix(rep(apply(Y,2,sd),sum(ps)),ncol = q,byrow = T)
  h_opt <- h - 1
  if(h_opt>0){
    for(h in 1:h_opt){
      if(h==1){
        Q2_cum <- sum(PRESS_y[h,])/sum(RSS0)
        Q2_cum_y <- PRESS_y[h,]/RSS0
      }else{
        Q2_cum <- Q2_cum*sum(PRESS_y[h,])/sum(RSS_y[h-1,])
        Q2_cum_y <- Q2_cum_y*PRESS_y[h,]/RSS_y[h-1,]
      }
    }
    Q2_cum <- 1-Q2_cum
    Q2_cum_y <- 1-Q2_cum_y
    i_0 <- 0
    for(k in 1:K){
      Us[[k]] <- u[i_0+1:ps[k],1:h_opt,drop=F]
      Bs[[k]] <- B[i_0+1:ps[k],,drop=F]
      i_0 <- ps[k]
    }
    mu_k <- colMeans(xs0_cbind[[1]])
    sd_k <- apply(xs0_cbind[[1]],2,sd)
    mu_y <- matrix(rep(colMeans(Y),n) ,nrow = n,byrow = T)
    sd_y <- matrix(rep(apply(Y,2,sd),n) ,nrow = n,byrow = T)
    y_est <- mu_y + (xs0_cbind[[1]]%*%B)*sd_y

    # Compute LOO error to get Q2_reg
    ERRORS_LOO <- rep(0,n)
    lambda_sol <- lambda[1:h_opt]
    Y_pred_all <- matrix(0,n,q)
    p <- ncol(xs0_cbind[[1]])
    for(ii in 1:n){
      X_train <- xs0_cbind[[1]][-ii,,drop=F]
      X_test <- xs0_cbind[[1]][ii,,drop=F]
      Y_train <- Y_init[-ii,,drop=FALSE]
      Y_test <- Y_init[ii,,drop=FALSE]
      mu_k <- colMeans(X_train)
      sd_k <- apply(X_train,2,sd)
      mu_y <- colMeans(Y_train)
      sd_y <- apply(Y_train,2,sd)
      ## Get LOO error
      model_ii <- get_model_ddsPLS2(X_train,Y_train,ncomp=h_opt,lambda[1:h_opt],method=method)
      B_ii <- model_ii$B
      Y_pred_all[ii,] <- (((X_test-mu_k)/sd_k)%*%B_ii)*sd_y+mu_y
    }
    ERRORS_LOO <- colSums((Y_init-Y_pred_all)^2)
    RSS0 = RSS_h_moins_1 <- colSums(Y_init^2)
    Q2_reg <- 1 - sum(ERRORS_LOO)/sum(RSS0)
    Q2_reg_y <- 1 - ERRORS_LOO/RSS0
    # Prepare outputs
    optimal_parameters <- list(lambda=lambda[1:h_opt],R=h_opt,
                               Q2_cum_y=Q2_cum_y,Q2_cum=Q2_cum,Q2_reg=Q2_reg,Q2_reg_y=Q2_reg_y,
                               ERRORS_LOO=ERRORS_LOO,Y_pred_LOO=Y_pred_all)
    parameters <- list(RSS0=RSS0,RSS_y=RSS_y[1:(h_opt+1),],PRESS_y=PRESS_y[1:(h_opt+1),],Q2_y=Q2_y[1:(h_opt+1),])
    out <- list(Us=Us,V=V,Bs=Bs,B_cbind=B,y_est=y_est,optimal_parameters=optimal_parameters,parameters=parameters,
                mu_k=mu_k,sd_k=sd_k,mu_y=colMeans(Y),sd_y=apply(Y,2,sd))
  }else{
    out <- NULL
  }
  out
}


#' Title
#'
#' @param Xs Xs
#' @param Y Y
#' @param lambdas lambdas
#' @param tau rate learning
#' @param NCORES number of cores
#' @param NZV near zero var
#'
#' @return
#' @export
#'
#' @useDynLib ddsPLS
Q2_global_ddsPLS <- function(Xs,Y,lambdas = 0.5,
                             tau=0.0975,NCORES=1,
                             NZV=1e-3,method=2){#,
  # threshs_to_test=NULL){
  n <- nrow(Xs[[1]])
  q <- ncol(Y)
  id_na <- lapply(Xs,function(xx){which(is.na(xx[,1]))})
  K <- length(Xs)
  PRESS = RSS = Q2_cum <- list()
  ncomps = Q2_cum <- rep(NA,length(lambdas))
  Q2_cum_y <- matrix(NA,length(lambdas),q)
  q <- ncol(Y)
  Y_init <- scale(Y)
  Xs_init <- lapply(Xs,scale)
  PRESS_y = RSS_y = Q2_h_y = Q2_h_sum <- list()
  xs0_cbind <- list(do.call(cbind,Xs_init))
  p <- ncol(xs0_cbind[[1]])
  for(i_l in 1:length(lambdas)){
    lam <- lambdas[i_l]
    RSS0_y = RSS_h_moins_1 <- colSums(Y_init^2)
    PRESS_y[[i_l]] = RSS_y[[i_l]] = Q2_h_y[[i_l]] <- matrix(NA,n,q)
    Q2_h_sum[[i_l]] <- rep(NA,n)
    test <- T
    xs0 <- Xs_init ; y0 <- Y_init
    h <- 0
    K_h <- 1:K
    while(test){
      h <- h + 1
      m_1 <- ddsPLS2(Xs = xs0_cbind,Y = y0,lam = lam,Rs = 1,method=method)
      y_est_train <- m_1$y_est
      RSS_y[[i_l]][h,] <- colSums((y_est_train-y0)^2)
      ERRORS_and_Bis <- do_loo(xs0_cbind,y0,NCORES=NCORES,lam=lam)
      PRESS_y[[i_l]][h,] <- colSums(ERRORS_and_Bis$PRESS_y)
      Q2_h_y[[i_l]][h,] <- 1-PRESS_y[[i_l]][h,]/RSS_h_moins_1
      Q2_h_sum[[i_l]][h] <- 1-sum(PRESS_y[[i_l]][h,])/sum(RSS_h_moins_1)
      # Update RESS h-1
      test_h <- max(Q2_h_y[[i_l]][h,])>tau
      if(!test_h){
        test <- F
        ncomps[i_l] <- h-1
      }else{
        m_plus <- model_PLS(xs0_cbind[[1]],y0,lam,R = 1, tau=tau,method=method,NZV=NZV)
        if(norm(m_plus$U_out)<NZV){
          test <- F
          ncomps[i_l] <- h-1
        }else{
          y0_plus <- y0-m_plus$e_y
          Var_rel <- sum(y0_plus^2)/sum(RSS0_y)
          if(Var_rel<NZV){
            test <- F
            ncomps[i_l] <- h-1
          }else{
            y0 <- m_plus$e_y
          }
        }
      }
      RSS_h_moins_1 <- RSS_y[[i_l]][h,]
    }
    if(ncomps[i_l]!=0){
      for(h in 1:ncomps[i_l]){
        if(h==1){
          Q2_cum[i_l] <- sum(PRESS_y[[i_l]][h,])/sum(RSS0_y)
          Q2_cum_y[i_l,] <- PRESS_y[[i_l]][h,]/RSS0_y
        }else{
          Q2_cum[i_l] <- Q2_cum[i_l]*sum(PRESS_y[[i_l]][h,])/sum(RSS_h_moins_1)
          Q2_cum_y[i_l,] <- Q2_cum_y[i_l,]*PRESS_y[[i_l]][h,]/RSS_h_moins_1
        }
      }
      Q2_cum[i_l] <- 1-Q2_cum[i_l]
      Q2_cum_y[i_l,] <- 1-Q2_cum_y[i_l,]
    }
  }
  Q2_y_max <- max(na.omit(Q2_cum))
  if(max(Q2_y_max)>tau){
    best_id <- which(Q2_cum==Q2_y_max)[1]
    if(length(best_id)>0){
      Q2_cum_out <- Q2_cum[best_id]
      Q2_cum_y_out <- Q2_cum_y[best_id,]
      best_lambda <- lambdas[best_id]
      best_ncomps <- ncomps[best_id]
      xs0 <- Xs_init ; y0 <- Y_init
      ps <- unlist(lapply(xs0,ncol))
      xs0_cbind <- list(do.call(cbind,xs0))
      model <- ddsPLS2(Xs = xs0_cbind,Y = y0,lam = best_lambda,Rs = best_ncomps,method=method)

      u <- model$u[[1]]
      t <- xs0_cbind[[1]]%*%u
      v <- model$V_super
      s <- y0%*%v
      B <- model$B[[1]]
      Us = Bs <- list()
      i_0 <- 0
      for(k in 1:(K)){
        Us[[k]] <- u[i_0+1:ps[k],,drop=F]
        Bs[[k]] <- B[i_0+1:ps[k],,drop=F]
        i_0 <- ps[k]
      }
      y_est <- model$y_est

      ERRORS_LOO <- rep(0,n)
      lambda_sol <- rep(best_lambda,best_ncomps)
      Y_pred_all <- matrix(0,n,q)
      p <- ncol(xs0_cbind[[1]])
      for(ii in 1:n){
        X_train <- xs0_cbind[[1]][-ii,,drop=F]
        X_test <- xs0_cbind[[1]][ii,,drop=F]
        Y_train <- Y_init[-ii,,drop=FALSE]
        Y_test <- Y_init[ii,,drop=FALSE]
        mu_k <- colMeans(X_train)
        sd_k <- apply(X_train,2,sd)
        mu_y <- colMeans(Y_train)
        sd_y <- apply(Y_train,2,sd)
        ## Get LOO error
        model_ii <- get_model_ddsPLS2(X_train,Y_train,ncomp=best_ncomps,lambda_sol,method=method)
        B_ii <- model_ii$B
        Y_pred_all[ii,] <- (((X_test-mu_k)/sd_k)%*%B_ii)*sd_y+mu_y
      }
      ERRORS_LOO <- colSums((Y_init-Y_pred_all)^2)
      mu_y <- matrix(rep(colMeans(Y_init),n) ,nrow = n,byrow = T)
      Q2_reg <- 1 - sum(ERRORS_LOO)/sum(RSS0_y)
      Q2_reg_y <- 1 - ERRORS_LOO/RSS0_y
      # Prepare outputs
      optimal_parameters <- list(lambda=best_lambda,R=best_ncomps,
                                 Q2_reg_y=Q2_reg_y,Q2_reg=Q2_reg,
                                 Q2_cum=Q2_cum_out,Q2_cum_y=Q2_cum_y_out,
                                 ERRORS_LOO=ERRORS_LOO,Y_pred_LOO=Y_pred_all)
      parameters <- list(RSS0_y=RSS0_y,
                         RSS_y=lapply(RSS_y,function(ri){ri[1:(h+1),]}),
                         PRESS_y=lapply(PRESS_y,function(ri){ri[1:(h+1),]}),
                         Q2_h_y = lapply(Q2_h_y,function(ri){ri[1:(h+1),]}),
                         Q2_h_sum = lapply(Q2_h_sum,function(ri){ri[1:(h+1)]}),
                         Q2_cum=Q2_cum,Q2_cum_y=Q2_cum_y)
      out <- list(model=model,Us=Us,Bs=Bs,B_cbind=do.call(rbind,Bs),y_est=y_est,
                  optimal_parameters=optimal_parameters,parameters=parameters)
    }else{
      out <- NULL
    }
  }else{
    out <- NULL
  }
  out
}


#' Title
#'
#' @param Xs Xs
#' @param Y Y
#' @param lambdas lambdas
#' @param tau rate learning
#' @param NCORES number of cores
#' @param NZV near zero var
#'
#' @return
#' @export
#'
#' @useDynLib ddsPLS
Q2_ddsPLS <- function(Xs,Y,lambdas = 0.5,
                             tau=0.0975,NCORES=1,
                             NZV=1e-3,type="local"){
  if(type=="local"){
    out <- Q2_local_ddsPLS(Xs,Y,lambdas = lambdas,
                           tau=0.0975,NCORES=NCORES,
                           NZV=NZV)
  }else{
    out <- Q2_global_ddsPLS(Xs,Y,lambdas = lambdas,
                           tau=0.0975,NCORES=NCORES,
                           NZV=NZV)
  }
  out
}


Q2_OLS <- function(X,Y){
  n <- nrow(X)
  SE <- rep(NA,n)
  for(pos_test in 1:n){
    X_train <- X[-pos_test,,drop=F]
    X_test <- X[pos_test,,drop=F]
    Y_train <- Y[-pos_test,,drop=FALSE]
    Y_test <- Y[pos_test,,drop=FALSE]
    mu_k <- colMeans(X_train)
    sd_k <- apply(X_train,2,sd)
    mu_y <- colMeans(Y_train)
    sd_y <- apply(Y_train,2,sd)
    B_th_all_ii <- MASS::ginv(scale(X_train))%*%scale(Y_train)
    y_pred <- (((X_test-mu_k)/sd_k)%*%B_th_all_ii)*sd_y+mu_y
    SE[pos_test] <- sum((Y_test-y_pred)^2)
  }
  PRESS <- mean(na.omit(SE))
  RSS <- mean(rowSums(Y-matrix(rep(colMeans(Y),n),nrow = n,byrow = T))^2)
  1-PRESS/RSS
}

Q2_OLS_th <- function(X,Y,B_th){
  n <- nrow(X)
  SE <- rep(NA,n)
  for(pos_test in 1:n){
    X_train <- X[-pos_test,,drop=F]
    X_test <- X[pos_test,,drop=F]
    Y_train <- Y[-pos_test,,drop=FALSE]
    Y_test <- Y[pos_test,,drop=FALSE]
    mu_k <- colMeans(X_train)
    sd_k <- apply(X_train,2,sd)
    mu_y <- colMeans(Y_train)
    sd_y <- apply(Y_train,2,sd)
    y_pred <- (((X_test-mu_k)/sd_k)%*%B_th)*sd_y+mu_y
    SE[pos_test] <- sum((Y_test-y_pred)^2)
  }
  PRESS <- mean(na.omit(SE))
  RSS <- mean(rowSums(Y-matrix(rep(colMeans(Y),n),nrow = n,byrow = T))^2)
  1-PRESS/RSS
}

