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
do_loo <- function(xs0,y0,lam=0,NCORES=1,method=2,deflatX=T){
  n <- nrow(xs0[[1]])
  q <- ncol(y0)
  p <- ncol(xs0[[1]])
  x0 <- do.call(cbind,xs0)
  PRESS_y <- matrix(NA,n,q)
  B <- list()
  K <- length(xs0)
  var_selected_y <- matrix(0,n,q)
  for(i in 1:n){
    X_train <- x0[-i,,drop=F]
    X_test <- x0[i,,drop=F]
    Y_train <- y0[-i,,drop=FALSE]
    Y_test <- y0[i,,drop=FALSE]
    mu_k <- colMeans(X_train)
    sd_k <- apply(X_train,2,sd)
    mu_y <- colMeans(Y_train)
    sd_y <- apply(Y_train,2,sd)
    X_train <- do.call(cbind,lapply(1:p,function(j,xx,mu_k,sd_k){if(sd_k[j]>1e-9){out <- (xx[,j]-mu_k[j])/sd_k[j]}else{out = xx[,j]};out},X_train,mu_k,sd_k))
    Y_train <- do.call(cbind,lapply(1:q,function(j,xx,mu_k,sd_k){if(sd_k[j]>1e-9){out <- (xx[,j]-mu_k[j])/sd_k[j]}else{out = xx[,j]};out},Y_train,mu_y,sd_y))
    m_1 <-  model_PLS(x = X_train,y=Y_train,lam=lam,deflatX = deflatX)# ddsPLS2(Xs = list(scale(X_train)),Y = scale(Y_train),lam = lam)
    B_i <- tcrossprod(m_1$U_out,m_1$V_out) # m_1$B[[1]]
    var_selected_y[i,which(colSums(abs(B_i))>1e-9)] <- 1
    # y_pred <- mu_y + (((X_test+mu_k)/sd_k)%*%B_i)*sd_y
    X_test_normalize <- do.call(cbind,lapply(1:p,function(j,xx,mu_k,sd_k){if(sd_k[j]>1e-9){out <- (xx[,j]-mu_k[j])/sd_k[j]}else{out = xx[,j]};out},X_test,mu_k,sd_k))
    y_pred <- mu_y + (X_test_normalize%*%B_i)*sd_y
    PRESS_y[i,] <- (y_pred-Y_test)^2
    B[[i]] <- B_i
  }
  list(id=1:n,PRESS_y=PRESS_y,B=B,var_selected_y=var_selected_y)
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
Q2_local_ddsPLS <- function(Xs,Y,lambdas = 0.5,deflatX=T,
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
  ps <- unlist(lapply(Xs_init,ncol))
  sd_x_inv_mat <- matrix(rep(unlist(lapply(apply(do.call(cbind,Xs),2,sd),function(ss){if(abs(ss)>1e-9){out <- 1/ss}else{out <- 0};out})),q),ncol = q,byrow = T)
  sd_y_mat <- matrix(rep(apply(Y,2,sd),sum(ps)),ncol = q,byrow = T)
  sd_y_x_inv <- sd_y_mat * sd_x_inv_mat
  RSS0 = RSS_h_moins_1 = RSS_h_moins_1_star <- colSums(Y_init^2)
  PRESS_y = RSS_y = Q2_y <- matrix(NA,n,q)
  PRESS_y_star = RSS_y_star = Q2_y_star <- matrix(NA,n,q)
  lambda = Q2_sum_star <- rep(NA,n)
  ps <- unlist(lapply(Xs_init,ncol))
  x0 <- do.call(cbind,Xs_init)
  u <- matrix(NA,sum(ps),n)
  P <- matrix(0,sum(ps),n)
  B <- matrix(0,sum(ps),q)
  Us = Bs = B_r_out <- list()
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
    PRESS_il_y_and_i_l <- foreach(i_l=1:Ls,.packages = "ddsPLS") %my_do% {
      out <- do_loo(list(x0),y0,lam=lambdas[i_l])
      out$PRESS_y_star <- colSums(out$PRESS_y*out$var_selected_y)
      out$PRESS_y <- colSums(out$PRESS_y)
      out$B <- NULL
      out
    }
    if(NCORES_w!=1){
      stopCluster(cl)
    }
    PRESS_il_y <- do.call(rbind,lapply(PRESS_il_y_and_i_l,function(ll){ll$PRESS_y}))
    # STAR
    PRESS_il_y_star <- do.call(rbind,lapply(PRESS_il_y_and_i_l,function(ll){ll$PRESS_y_star}))
    PRESS_il_y_star[which(rowSums(PRESS_il_y_star)<1e-9),] <- NA
    RSS_il_y = var_sel <- matrix(NA,Ls,q)
    RSS_il_y_star <- matrix(0,Ls,q)
    for(i_l in 1:Ls){
      lam <- lambdas[i_l]
      m_1 <- model_PLS(x0,y0,lam,deflatX=deflatX)
      B_i_l <- m_1$B
      id_B <- which(colSums(abs(B_i_l))>1e-9)
      y_est_train <- m_1$y_est
      var_sel[i_l,id_B] <- 1
      RSS_il_y[i_l,] <- colSums((y_est_train-y0)^2)
      if(length(id_B)>0){
        RSS_il_y_star[i_l,id_B] <- colSums(((y_est_train-y0)^2)[,id_B,drop=F])
      }
    }
    Q2_h_y <- 1-PRESS_il_y/matrix(rep(RSS_h_moins_1,Ls),nrow = Ls,byrow = T)
    Q2_h_sum <- 1-rowSums(PRESS_il_y)/sum(RSS_h_moins_1)
    # STAR
    RSS_h_moins_1_star <- t(apply(var_sel,1,function(vv){vv*RSS_h_moins_1}))
    if(ncol(RSS_h_moins_1_star)!=q) RSS_h_moins_1_star <- t(RSS_h_moins_1_star)
    RSS_h_moins_1_star[which(is.na(RSS_h_moins_1_star))] <- 0
    Q2_h_sum_star<- 1-rowSums(PRESS_il_y_star)/rowSums(RSS_h_moins_1_star)
    # max_Q2_y <- apply(Q2_h_y,MARGIN = 1,max)
    best_id_h <- which(Q2_h_sum_star==max(na.omit(Q2_h_sum_star)))[1] # which.max( max_Q2_y)#  Q2_h_sum)#
    if(length(best_id_h)>0){
      test_h <- Q2_h_sum_star[best_id_h]>tau # max_Q2_y[best_id_h]>tau # Q2_h_sum[best_id_h]>tau #
      test_h <- test_h & Q2_h_sum[best_id_h]>tau
      if(!test_h){
        test <- F
      }else{
        best_lam_h <- lambdas[best_id_h]
        lambda[h] <- best_lam_h
        RSS_y[h,] <- RSS_il_y[best_id_h,]
        if(h==1){
          RSS0_star <- RSS_h_moins_1_star[best_id_h,]
        }
        RSS_h_moins_1 <- RSS_y[h,]
        PRESS_y[h,] <- PRESS_il_y[best_id_h,]
        Q2_y[h,] <- Q2_h_y[best_id_h,]
        # STAR
        RSS_y_star[h,] <- RSS_il_y_star[best_id_h,]
        RSS_h_moins_1_star <- RSS_y_star[h,]
        PRESS_y_star[h,] <- PRESS_il_y_star[best_id_h,]
        Q2_sum_star[h] <- Q2_h_sum_star[best_id_h]
        # Get the regression matrix of the optimal model
        m_plus <- model_PLS(x0,y0,best_lam_h,R = 1,NZV=NZV,deflatX=deflatX)
        if(norm(m_plus$U_out)>NZV){
          u[,h] <- m_plus$U_out
          t_r <- x0%*%u[,h,drop=F]
          P[,h] <- crossprod(x0,t_r)/sum(t_r^2)
          if(h!=1){
            for(s_r in (h-1):1){
              u[,h] <- u[,h]-u[,s_r]*sum(P[,s_r]*u[,h])
            }
          }
          if(sum(u[,h]^2)>NZV){
            B_r <- (tcrossprod(u[,h,drop=F],m_plus$V_out))#*sd_y_x_inv
            V_r <- m_plus$V_out
            x0_plus <- tcrossprod(t_r,P[,h])
            y0_plus <- x0%*%B_r
            id_y_sel <- which(abs(V_r)>NZV)
            Var_rel <- sum(y0_plus^2)/sum(RSS0)
            if(Var_rel<NZV){
              test <- F
            }else{
              B <- B + B_r*sd_y_x_inv
              B_r_out[[h]] <- B_r*sd_y_x_inv
              V <- cbind(V,V_r)
              y0 <- y0 - y0_plus
              if(deflatX){
                x0 <- x0 - x0_plus
              }
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
  h_opt <- h - 1
  x0_center <- scale(do.call(cbind,Xs_init),scale = F)
  if(h_opt>0){
    Q2_cum_star <- 1-prod(1-Q2_sum_star[1:h_opt])
    for(h in 1:h_opt){
      if(h==1){
        Q2_cum <- 1- sum(PRESS_y[h,])/sum(RSS0)
        Q2_cum_y <- Q2_y[h,] # PRESS_y[h,]/RSS0
        # STAR
        # Q2_cum_star <- 1-sum(PRESS_y_star[h,])/sum(RSS0_star)
      }else{
        Q2_cum <- 1-(1-Q2_cum)*sum(PRESS_y[h,])/sum(RSS_y[h-1,])
        Q2_cum_y <- 1-(1-Q2_cum_y)*PRESS_y[h,]/RSS_y[h-1,]
        # STAR
        # Q2_cum_star <- 1-(1-Q2_cum_star)*sum(PRESS_y_star[h,])/sum(RSS_y_star[h-1,])
      }
    }
    i_0 <- 0
    for(k in 1:K){
      Us[[k]] <- u[i_0+1:ps[k],1:h_opt,drop=F]
      Bs[[k]] <- B[i_0+1:ps[k],,drop=F]
      i_0 <- ps[k]
    }
    mu_y <- matrix(rep(colMeans(Y),n) ,nrow = n,byrow = T)
    y_est <- mu_y + x0_center%*%B

    # Compute LOO error to get Q2_reg
    ERRORS_LOO <- rep(0,n)
    lambda_sol <- lambda[1:h_opt]
    Y_pred_all <- matrix(0,n,q)
    p <- ncol(x0_center)
    var_sel <- matrix(0,n,q)
    ## Get LOO error
    for(ii in 1:n){
      X_train <- x0_center[-ii,,drop=F]
      X_test <- x0_center[ii,,drop=F]
      Y_train <- Y[-ii,,drop=FALSE]
      Y_test <- Y[ii,,drop=FALSE]
      mu_k <- colMeans(X_train)
      mu_y <- colMeans(Y_train)
      sd_x_inv <- unlist(lapply(apply(X_train,2,sd),function(ss){if(abs(ss)>1e-9){out <- 1/ss}else{out <- 0};out}))
      m_plus <- model_PLS(X_train,Y_train,lambda_sol,R = h_opt,
                          NZV=NZV,deflatX=deflatX)
      B_ii <- m_plus$B
      var_sel[ii,which(colSums(abs(B_ii))>1e-9)] <- 1
      Y_pred_all[ii,] <- mu_y + ((X_test-mu_k)*sd_x_inv)%*%B_ii
    }
    ERRORS_LOO <- colSums((Y_init-Y_pred_all)^2)
    Q2_reg <- 1 - sum(ERRORS_LOO)/sum(RSS0)
    Q2_reg_y <- 1 - ERRORS_LOO/RSS0
    # STAR
    m_star <- model_PLS(x0,Y_init,lambda_sol,R = h_opt,
                        NZV=NZV,deflatX=deflatX)
    ERRORS_LOO_star <- colSums(((Y_init-Y_pred_all)^2)*var_sel)
    RSS0_star_model <- colSums(Y_init^2)[which(rowSums(abs(m_star$V_out))>1e-9)]
    Q2_reg_star <- 1 - sum(ERRORS_LOO_star)/sum(RSS0_star_model)
    explained_variance <- unlist(lapply(B_r_out,function(b,xx){
      sum((xx%*%b)^2)/sum(RSS0)*100
    },x0))
    # Prepare outputs
    optimal_parameters <- list(lambda=lambda_sol,R=h_opt,
                               Q2_cum_y=Q2_cum_y,Q2_cum=Q2_cum,
                               Q2_cum_star=Q2_cum_star,Q2_sum_star=Q2_sum_star,
                               Q2_reg=Q2_reg,Q2_reg_y=Q2_reg_y,Q2_reg_star=Q2_reg_star,
                               ERRORS_LOO=ERRORS_LOO,ERRORS_LOO_star=ERRORS_LOO_star,
                               Y_pred_LOO=Y_pred_all)
    parameters <- list(RSS0=RSS0,RSS0_star=RSS0_star,RSS_y=RSS_y[1:(h_opt+1),],
                       PRESS_y=PRESS_y[1:(h_opt+1),],Q2_y=Q2_y[1:(h_opt+1),])
    out <- list(optimal_parameters=optimal_parameters,Us=Us,V=V,
                explained_variance=explained_variance,
                Bs=Bs,B_cbind=B,B_r=B_r_out,
                y_est=y_est,parameters=parameters,
                mu_k=mu_k,mu_y=colMeans(Y),sd_y=apply(Y,2,sd))
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
Q2_global_ddsPLS <- function(Xs,Y,lambdas = 0.5,deflatX=T,
                             tau=0.0975,NCORES=1,
                             NZV=1e-3,method=2){#,
  # threshs_to_test=NULL){
  n <- nrow(Xs[[1]])
  q <- ncol(Y)
  id_na <- lapply(Xs,function(xx){which(is.na(xx[,1]))})
  K <- length(Xs)
  PRESS = RSS <- list()
  ncomps = Q2_cum = Q2_cum_star <- rep(NA,length(lambdas))
  Q2_cum_y <- matrix(NA,length(lambdas),q)
  q <- ncol(Y)
  Y_init <- scale(Y)
  Xs_init <- lapply(Xs,scale)
  ps <- unlist(lapply(Xs,ncol))
  sd_x_inv_mat <- matrix(rep(unlist(lapply(apply(do.call(cbind,Xs),2,sd),function(ss){if(abs(ss)>1e-9){out <- 1/ss}else{out <- 0};out})),q),ncol = q,byrow = T)
  sd_y_mat <- matrix(rep(apply(Y,2,sd),sum(ps)),ncol = q,byrow = T)
  sd_y_x_inv <- sd_y_mat * sd_x_inv_mat
  PRESS_y = RSS_y = Q2_h_y = Q2_h_sum = Q2_h_sum_star = RSS_h_star <- list()
  p <- sum(unlist(lapply(Xs,ncol)))
  for(i_l in 1:length(lambdas)){
    lam <- lambdas[i_l]
    RSS0_y = RSS_h_moins_1 <- colSums(Y_init^2)
    PRESS_y[[i_l]] = RSS_y[[i_l]] = RSS_h_star[[i_l]] = Q2_h_y[[i_l]] <- matrix(NA,n,q)
    Q2_h_sum[[i_l]] = Q2_h_sum_star[[i_l]] <- rep(NA,n)
    test <- T
    x0 <- do.call(cbind,Xs_init) ; y0 <- Y_init
    h <- 0
    K_h <- 1:K
    while(test){
      h <- h + 1
      m_1 <- model_PLS(x0,y0,lam = lam,deflatX=deflatX)
      y_est_train <- m_1$y_est
      # STAR
      var_sel <- which(rowSums(abs(m_1$V_out))>1e-9)
      if(h==1){
        RSS_h_moins_1 <- RSS0_y
      }
      RSS_y[[i_l]][h,] <- colSums((y_est_train-y0)^2)
      ERRORS_and_Bis <- do_loo(list(x0),y0,NCORES=NCORES,lam=lam)
      PRESS_y[[i_l]][h,] <- colSums(ERRORS_and_Bis$PRESS_y)
      Q2_h_y[[i_l]][h,] <- 1-PRESS_y[[i_l]][h,]/RSS_h_moins_1
      Q2_h_sum[[i_l]][h] <- 1-sum(PRESS_y[[i_l]][h,])/sum(RSS_h_moins_1)
      # STAR
      var_sel <- ERRORS_and_Bis$var_selected_y
      var_sel_model <- c(m_1$V_out)
      var_sel_model[which(abs(var_sel_model)>1e-9)] <- 1
      var_sel_model[which(abs(var_sel_model)<1e-9)] <- 0
      RSS_h_moins_1_star <- RSS_h_moins_1*var_sel_model
      PRESS_il_y_star <- colSums(ERRORS_and_Bis$PRESS_y*var_sel)
      if(sum(RSS_h_moins_1_star)==0){
        Q2_h_sum_star[[i_l]][h] <- NA
      }else{
        Q2_h_sum_star[[i_l]][h] <- 1-sum(PRESS_il_y_star)/sum(RSS_h_moins_1_star)
      }
      # Update RESS h-1
      test_h <- Q2_h_sum_star[[i_l]][h]>tau#max(Q2_h_y[[i_l]][h,])>tau
      if(is.na(test_h)) test_h <- F
      if(!test_h){
        test <- F
        ncomps[i_l] <- h-1
      }else{
        m_plus <- model_PLS(x0,y0,lam,R = 1, NZV=NZV,deflatX=deflatX)
        if(norm(m_plus$U_out)<NZV){
          test <- F
          ncomps[i_l] <- h-1
        }else{
          y0_plus <- y0-m_plus$e_y
          x0_plus <- x0-m_plus$e_x
          Var_rel <- sum(y0_plus^2)/sum(RSS0_y)
          if(Var_rel<NZV){
            test <- F
            ncomps[i_l] <- h-1
          }else{
            y0 <- m_plus$e_y
            x0 <- m_plus$e_x
          }
        }
      }
      RSS_h_moins_1 <- RSS_y[[i_l]][h,]
    }
    if(ncomps[i_l]!=0){
      for(h in 1:ncomps[i_l]){
        if(h==1){
          Q2_cum[i_l] <- 1-sum(PRESS_y[[i_l]][h,])/sum(RSS0_y)
          Q2_cum_y[i_l,] <- 1-PRESS_y[[i_l]][h,]/RSS0_y
          # STAR
          Q2_cum_star[i_l] <- Q2_h_sum_star[[i_l]][h]
        }else{
          Q2_cum[i_l] <- 1-(1-Q2_cum[i_l])*(1-Q2_h_sum[[i_l]][h])#sum(PRESS_y[[i_l]][h,])/sum(RSS_h_moins_1)
          Q2_cum_y[i_l,] <- 1-(1-Q2_cum_y[i_l,])*(1-Q2_h_y[[i_l]][h,])#PRESS_y[[i_l]][h,]/RSS_h_moins_1
          # STAR
          Q2_cum_star[i_l] <- 1-(1-Q2_cum_star[i_l])*(1-Q2_h_sum[[i_l]][h])
        }
      }
    }
  }
  Q2_cum_star_max <- max(na.omit(Q2_cum_star))
  B_r <- list()
  if(Q2_cum_star_max>tau){
    best_id <- which(Q2_cum_star==Q2_cum_star_max)[1]
    if(length(best_id)>0){
      Q2_cum_out <- Q2_cum[best_id]
      Q2_cum_y_out <- Q2_cum_y[best_id,]
      best_lambda <- lambdas[best_id]
      best_ncomps <- ncomps[best_id]
      xs0 <- Xs_init
      y0 <- Y
      ps <- unlist(lapply(xs0,ncol))
      x0 <-  do.call(cbind,Xs)# do.call(cbind,xs0)
      # model <- ddsPLS2(Xs = list(x0),Y = y0,lam = best_lambda,Rs = best_ncomps)
      model <- model_PLS(x0,y0,lam = best_lambda,R = best_ncomps)
      B_r <- model$B_r
      u <- model$U_out #model$u[[1]]
      B <- model$B
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
      Y_pred_all = var_sel <- matrix(0,n,q)
      p <- ncol(x0)
      for(ii in 1:n){
        X_train <- x0[-ii,,drop=F]
        X_test <- x0[ii,,drop=F]
        Y_train <- y0[-ii,,drop=FALSE]
        Y_test <- y0[ii,,drop=FALSE]
        mu_k <- colMeans(X_train)
        mu_y <- colMeans(Y_train)
        ## Get LOO error
        model_ii <- model_PLS(X_train,Y_train,lam = lambda_sol,R = best_ncomps)
        # get_model_ddsPLS2(X_train,Y_train,ncomp=best_ncomps,lambda_sol)
        B_ii <- model_ii$B
        var_sel_ii <- which(colSums(abs(B_ii))>1e-9)
        if(length(var_sel_ii)>0){
          var_sel[ii,var_sel_ii] <- 1
        }
        Y_pred_all[ii,] <- mu_y + (X_test-mu_k)%*%B_ii
      }
      ERRORS_LOO <- colSums((Y_init-Y_pred_all)^2)
      ERRORS_LOO_star <- colSums((Y_init-Y_pred_all)^2*var_sel)
      mu_y <- matrix(rep(colMeans(Y_init),n) ,nrow = n,byrow = T)
      Q2_reg <- 1 - sum(ERRORS_LOO)/sum(RSS0_y)
      Q2_reg_y <- 1 - ERRORS_LOO/RSS0_y
      # STAR
      var_sel_model <- which(rowSums(abs(model$V_out))>1e-9)
      Q2_reg_star <- 1 - sum(ERRORS_LOO_star)/sum(RSS0_y[var_sel_model])
      explained_variance <- unlist(lapply(B_r,function(b,xx){
        sum((xx%*%b)^2)/sum(RSS0_y)*100
      },x0))
      # Prepare outputs
      optimal_parameters <- list(lambda=best_lambda,R=best_ncomps,
                                 Q2_reg_y=Q2_reg_y,Q2_reg=Q2_reg,
                                 Q2_reg_star=Q2_reg_star,
                                 Q2_cum=Q2_cum_out,Q2_cum_y=Q2_cum_y_out,
                                 Q2_cum_star=Q2_cum_star_max,
                                 ERRORS_LOO=ERRORS_LOO,
                                 ERRORS_LOO_star=ERRORS_LOO_star,
                                 Y_pred_LOO=Y_pred_all)
      parameters <- list(RSS0_y=RSS0_y,
                         RSS_y=lapply(RSS_y,function(ri){ri[1:(h+1),]}),
                         PRESS_y=lapply(PRESS_y,function(ri){ri[1:(h+1),]}),
                         Q2_h_y = lapply(Q2_h_y,function(ri){ri[1:(h+1),]}),
                         Q2_h_sum = lapply(Q2_h_sum,function(ri){ri[1:(h+1)]}),
                         Q2_cum=Q2_cum,Q2_cum_y=Q2_cum_y)
      out <- list(optimal_parameters=optimal_parameters,model=model,Us=Us,explained_variance=explained_variance,
                  Bs=Bs,B_cbind=B,B_r=B,
                  y_est=y_est,
                  parameters=parameters)
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
Q2_ddsPLS <- function(Xs,Y,lambdas = 0.5,deflat=T,
                      tau=0.0975,NCORES=1,
                      NZV=1e-3,type="local"){
  if(type=="local"){
    out <- Q2_local_ddsPLS(Xs,Y,lambdas = lambdas,deflat=deflat,
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

