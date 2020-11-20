#' Title
#'
#' @param x0 x0
#' @param y0 y0
#' @param n n
#' @param p p
#' @param q q
#' @param COV COV
#' @param abs_COV abs_COV
#' @param max_COV max_COV
#' @param lam lam
#' @param NZV NZV
#'
#' @return
#' @export
#'
#' @useDynLib ddsPLS
do_one_component <- function(x0,y0,n,p,q,COV,abs_COV,max_COV,lam,remove_COV=NULL,NZV=1e-3){
  max_cov_y <- apply(abs_COV,1,max)
  max_cov_x <- apply(abs_COV,2,max)
  id_y_high <- which(max_cov_y>lam)
  id_x_high <- which(max_cov_x>lam)
  if(length(id_x_high)>0 & length(id_y_high)>0){
    COV_high <- COV[id_y_high,id_x_high,drop=F]
    abs_COV_high <- abs(COV_high)
    COV_COV_high <- abs_COV_high - lam
    if(!is.null(remove_COV)){
      if(norm(remove_COV)>1e-9){
        COV_COV_high <- COV_COV_high - abs(remove_COV[id_y_high,id_x_high,drop=F])
      }
    }
    COV_COV_high[which(COV_COV_high<0)] <- 0
    COV_COV_high <- COV_COV_high*sign(COV_high)
    # Do svd
    U0 <- matrix(0,p,1)
    V0 <- matrix(0,q,1)
    ## X part
    if(norm(COV_COV_high)>1e-9){
      # model_NULL <- svd(COV_high,nu = 1,nv = 1)
      # svd_XY <- svd(COV_COV_high,nv = 1,nu=1)#svd(COV_COV_high,nv = 1,nu=0)#

      svd_mix_Y <- svd(tcrossprod(COV_COV_high,COV_high),nu = 0,nv = 1)
      svd_mix_X <- svd(crossprod(COV_COV_high,COV_high),nu = 0,nv = 1)
      # u_x_no_std <- t(COV_COV_high)%*%model_NULL$u
      U0[id_x_high,] <- svd_mix_X$v#u_x_no_std/sqrt(sum(u_x_no_std^2))#svd_XY$v#
      ## Y part
      # u_y_no_std <- COV_COV_high%*%model_NULL$v
      V0[id_y_high,] <- svd_mix_Y$v#u_y_no_std/sqrt(sum(u_y_no_std^2))#svd_XY$u#
    }
  }else{
    U0 <- matrix(0,p,1)
    V0 <- matrix(0,q,1)
  }
  t <- x0%*%U0
  V_svd0 <- V0%*%crossprod(V0,crossprod(y0,t)) #crossprod(y0,t) #
  norm_t_0 <- sum(t^2)
  if(norm_t_0>NZV){
    V_svd <- V_svd0/norm_t_0
  }else{
    U0 <- matrix(0,p,1)
    V_svd <- matrix(0,q,1)
  }
  ##
  list(t=t,U0=U0,V_svd=V_svd,V0=V0)
}
########

#' Title
#'
#' @param x x
#' @param y y
#' @param deflatX deflatX
#' @param R R
#' @param NZV NZV
#' @param to.scale to.scale
#'
#' @return
#' @export
#'
#' @useDynLib ddsPLS
model_PLS <- function(x,y,lam,deflatX=T,R=1,remove_COV=NULL,RSS0=NULL,NZV=1e-3,to.scale=T,COV_init=NULL){
  p <- ncol(x)
  q <- ncol(y)
  n <- nrow(y)
  if(to.scale){
    mu_x <- matrix(rep(colMeans(x),n),ncol = p,byrow = T)
    mu_y <- matrix(rep(colMeans(y),n),ncol = q,byrow = T)
    sd_x_inv <- unlist(lapply(apply(x,2,sd),function(ss){if(abs(ss)>1e-9){out <- 1/ss}else{out <- 0};out}))
    sd_x_inv_mat <- matrix(rep(sd_x_inv,q),ncol = q,byrow = T)
    sd_y_mat <- matrix(rep(apply(y,2,sd),p),ncol = q,byrow = T)
    sd_y_x_inv <- sd_y_mat * sd_x_inv_mat
    x0 <- scale(x)
    y0 <- scale(y)
  }else{
    x0 <- x
    y0 <- y
    sd_y_mat <- matrix(rep(apply(y,2,sd),p),ncol = q,byrow = T)
  }
  if(is.null(RSS0)){
    RSS0 <- sum(sd_y_mat^2)
  }else{
    RSS0 <- sum(RSS0)
  }
  var_y_init <- sum(y^2)
  U_out = U_star <- matrix(0,p,R)
  score_x <- matrix(0,n,R)
  V_out <- matrix(0,q,R)
  bXr <- matrix(0,R,p)
  bYr <- matrix(0,R,q)
  B <- matrix(0,p,q)
  B_r <- list()
  if(to.scale){
    y_est <- mu_y
  }else{
    y_est <- matrix(0,n,q)
  }
  stop <- F
  no_model <- F
  var_expl <- rep(NA,R)
  # covs <- list()
  for(r in 1:R){
    if(r==1 & !is.null(COV_init)){
      COV <- COV_init
    }else{
      COV <- crossprod(y0,x0)/(n-1)
    }
    # covs[[r]] <- COV
    abs_COV <- abs(COV)
    max_COV <- max(na.omit(c(abs_COV)))
    lam_r <- lam
    if(length(lam)>1) lam_r <- lam[r]
    if(lam_r<max_COV){

      c_h <- do_one_component(x0 = x0,y0 = y0,n = n,p = p,q = q,COV = COV,abs_COV = abs_COV,
                              max_COV=max_COV,lam = lam_r,NZV=NZV,remove_COV=remove_COV)
      t_r <- c_h$t ; U0  <- c_h$U0 ; V_svd  <- c_h$V_svd ; V0  <- c_h$V0
      ## DEFLAT ##
      if(sum(U0^2)>NZV){
        score_x[,r] <- t_r
        bt <- crossprod(t_r,x0)/sum(t_r^2)
        U_out[,r] <- U0
        U_star[,r] <- U0
        V_out[,r] <- V_svd
        B_r_0 <- tcrossprod(U0,V_svd)
        y_est <- y_est + x0%*%B_r_0
        bXr[r,] <- bt
        bYr[r,] <- t(V_svd)
        if(r!=1){
          for(s_r in (r-1):1){
            U_star[,r] <- U_star[,r]-U_star[,s_r]*sum(bXr[s_r,]*U_star[,r])
          }
        }
        if(to.scale){
          B_r[[r]] <- tcrossprod(U_star[,r],V_svd)*sd_y_x_inv
        }else{
          B_r[[r]] <- tcrossprod(U_star[,r],V_svd)
        }
        B <- B + B_r[[r]]
        y_plus_un <- tcrossprod(t_r,V_out[,r,drop=F])
        x0_plus <- t_r%*%bt
        var_y_plus_un <- sum(y_plus_un^2)
        var_expl[r] <- var_y_plus_un/RSS0
        y0 <- y0 - y_plus_un
        if(deflatX){
          x0 <- x0 - x0_plus
        }
      }
    }else{
      V0 <- matrix(0,q,1)
    }
  }
  if(to.scale){
    y_est <- y_est - mu_x%*%B
  }
  list(no_model=no_model,U_out=U_out,U_star=U_star,
       V_out=V_out,V_optim=V0,
       B=B,B_r=B_r,var_expl=var_expl,#covs=covs,
       score_x=score_x,y_est=y_est,bXr=bXr,bYr=bYr,e_x=x0,e_y=y0)
}
########

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
do_loo <- function(xs0,y0,lam=0,NCORES=1,deflatX=T){
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
    m_1 <-  model_PLS(x = X_train,y=Y_train,lam=lam,remove_COV=remove_COV,deflatX = deflatX,to.scale = F)
    B_i <- tcrossprod(m_1$U_out,m_1$V_out)
    var_selected_y[i,which(colSums(abs(B_i))>1e-9)] <- 1
    X_test_normalize <- do.call(cbind,lapply(1:p,function(j,xx,mu_k,sd_k){if(sd_k[j]>1e-9){out <- (xx[,j]-mu_k[j])/sd_k[j]}else{out = xx[,j]};out},X_test,mu_k,sd_k))
    y_pred <- mu_y + (X_test_normalize%*%B_i)*sd_y
    PRESS_y[i,] <- (y_pred-Y_test)^2
    B[[i]] <- B_i
  }
  list(id=1:n,PRESS_y=PRESS_y,B=B,var_selected_y=var_selected_y)
}

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
do_loo <- function(xs0,y0,lam=0,remove_COV=NULL,NCORES=1,deflatX=T){
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
    m_1 <-  model_PLS(x = X_train,y=Y_train,lam=lam,remove_COV=remove_COV,deflatX = deflatX,to.scale = F)
    B_i <- tcrossprod(m_1$U_out,m_1$V_out)
    var_selected_y[i,which(colSums(abs(B_i))>1e-9)] <- 1
    X_test_normalize <- do.call(cbind,lapply(1:p,function(j,xx,mu_k,sd_k){if(sd_k[j]>1e-9){out <- (xx[,j]-mu_k[j])/sd_k[j]}else{out = xx[,j]};out},X_test,mu_k,sd_k))
    y_pred <- mu_y + (X_test_normalize%*%B_i)*sd_y
    PRESS_y[i,] <- (y_pred-Y_test)^2
    B[[i]] <- B_i
  }
  list(id=1:n,PRESS_y=PRESS_y,B=B,var_selected_y=var_selected_y)
}

bootstrap_pls <- function(X_init,Y_init,B_previous,u,v,h=1,lambdas=0){
  N_lambdas <- length(lambdas)
  n <- nrow(X_init)
  q <- ncol(Y_init)
  p <- ncol(X_init)
  id_IB <- sample(1:n,size = n,replace = T)
  id_OOB <- (1:n)[-sort(unique(id_IB))]
  while(length(id_OOB)==0){
    id_IB <- sample(1:n,size = n,replace = T)
    id_OOB <- (1:n)[-sort(unique(id_IB))]
  }
  X_train <- X_init[id_IB,,drop=F]
  X_test <- X_init[id_OOB,,drop=F]
  Y_train <- Y_init[id_IB,,drop=F]
  Y_test <- Y_init[id_OOB,,drop=F]
  mu_k <- colMeans(X_train)
  sd_k <- apply(X_train,2,sd)#rep(1,length(mu_k))#
  mu_y <- colMeans(Y_train)
  sd_y <- apply(Y_train,2,sd)#apply(Y_train,2,sd)
  # COV_init_no_sd <- crossprod(Y_train,X_train)/(n-1)
  X_train <- do.call(cbind,lapply(1:p,function(j,xx,mu_k,sd_k){if(sd_k[j]>1e-9){out <- (xx[,j]-mu_k[j])/sd_k[j]}else{out = xx[,j]};out},X_train,mu_k,sd_k))
  X_test_normalize <- do.call(cbind,lapply(1:p,function(j,xx,mu_k,sd_k){if(sd_k[j]>1e-9){out <- (xx[,j]-mu_k[j])/sd_k[j]}else{out = xx[,j]};out},X_test,mu_k,sd_k))
  Y_train <- do.call(cbind,lapply(1:q,function(j,xx,mu_k,sd_k){if(sd_k[j]>1e-9){out <- (xx[,j]-mu_k[j])/sd_k[j]}else{out = xx[,j]};out},Y_train,mu_y,sd_y))
  Y_test_normalize <- do.call(cbind,lapply(1:q,function(j,xx,mu_k,sd_k){if(sd_k[j]>1e-9){out <- (xx[,j]-mu_k[j])/sd_k[j]}else{out = xx[,j]};out},Y_test,mu_y,sd_y))
  U_reconstruct <- matrix(0,p,h)
  V_reconstruct <- matrix(0,q,h)
  t_all <- matrix(0,n,h)
  X_defla = Y_defla <- list()
  X_defla[[1]] <- X_train
  Y_defla[[1]] <- Y_train
  B_youyou <- matrix(0,p,q)
  X_r <- X_train;Y_r <- Y_train
  if(h>1){
    U_reconstruct[,1:(h-1)] <- u[,1:(h-1)]
    for(r in 1:(h-1)){
      t_r <- X_r%*%u[,r]
      t_all[,r] <- t_r
      bt <- crossprod(t_r,X_r)/sum(t_r^2)
      V_svd0 <- v[,r]%*%crossprod(v[,r],crossprod(Y_r,t_r))
      norm_t_0 <- sum(t_r^2)
      if(norm_t_0>1e-9){
        V_svd <- V_svd0/norm_t_0
      }else{
        V_svd <- V_svd0*0
      }
      y_r <- tcrossprod(t_r,V_svd)
      V_reconstruct[,r] <- V_svd
      # Do deflation
      Y_r <- Y_r - y_r
      X_r <- X_r - t_r%*%bt
      X_defla[[r+1]] <- X_r
      Y_defla[[r+1]] <- Y_r
    }
    # Create weights
    if(h>2){
      for(r in 2:(h-1)){
        for(s_r in (r-1):1){
          x_s_r <- X_defla[[s_r]]
          u_s_r <- U_reconstruct[,s_r,drop=F]
          u_r <- U_reconstruct[,r]
          t_s_r <- t_all[,s_r]
          bt <- c(crossprod(t_s_r,x_s_r)/sum(t_s_r^2))
          U_reconstruct[,r] <- u_r-u_s_r*sum(bt*u_r)
        }
      }
    }
    # Build regression matrix
    if(h>1){
      for(r in 1:(h-1)){
        B_youyou <- B_youyou + tcrossprod(U_reconstruct[,r],V_reconstruct[,r])
      }
    }
  }
  vars_expl_star = vars_expl = Q2_star = Q2 = model_exists <- rep(0,N_lambdas)
  B_out <- matrix(0,q*p,N_lambdas)

  u_out <- matrix(0,p,N_lambdas)
  V_optim_phi = V_model <- matrix(0,q,N_lambdas)
  for(i_l in 1:N_lambdas){
    m_1 <-  model_PLS(x = X_r,y=Y_r,lam=lambdas[i_l],to.scale = F)
    # B <- tcrossprod(m_1$U_out,m_1$V_out)
    # B_out[,i_l] <- c(B)
    # u_out[,i_l] <- m_1$U_out
    # V_out[,i_l] <- m_1$V_optim
    u_il <- m_1$U_out
    V_il <- m_1$V_optim
    u_out[,i_l] <- u_il
    V_optim_phi[,i_l] <- V_il
    t_r <- X_r%*%u_il
    bt <- crossprod(t_r,X_r)/sum(t_r^2)
    V_svd0 <- V_il%*%crossprod(V_il,crossprod(Y_r,t_r))
    norm_t_0 <- sum(t_r^2)
    if(norm_t_0>1e-9){
      V_svd <- V_svd0/norm_t_0
    }else{
      V_svd <- V_svd0*0
    }
    # Create weights
    U_reconstruct[,h] <- u_il
    if(h>1){
      for(s_r in (h-1):1){
        x_s_r <- X_defla[[s_r]]
        u_s_r <- U_reconstruct[,s_r,drop=F]
        u_r <- U_reconstruct[,h]
        t_s_r <- t_all[,s_r]
        bt <- crossprod(t_s_r,x_s_r)/sum(t_s_r^2)
        U_reconstruct[,h] <- u_r-u_s_r*sum(bt*u_r)
      }
    }
    # Modify regression matrix
    B_next <- tcrossprod(U_reconstruct[,h],V_svd)
    # B_youyou <- B_youyou + B_next
    # Predict values
    # # Component h
    # y_train_pred_h <- X_train%*%B_next
    # y_test_pred_h <- X_test_normalize%*%B_next
    # Y_train_pred_h <- t(apply(y_train_pred_h,1,function(yi){ mu_y + (yi)*sd_y }))
    # Y_test_pred_h <- t(apply(y_test_pred_h,1,function(yi){ mu_y + (yi)*sd_y }))
    # All components
    B_all <- B_next + B_youyou
    y_train_pred <- X_train%*%B_all
    y_test_pred <- X_test_normalize%*%B_all
    Y_train_pred <- t(apply(y_train_pred,1,function(yi){ mu_y + (yi)*sd_y }))
    Y_test_pred <- t(apply(y_test_pred,1,function(yi){ mu_y + (yi)*sd_y }))
    # Previous components
    y_test_pred_RSS <- X_test_normalize%*%B_previous
    Y_test_pred_RSS <- t(apply(y_test_pred_RSS,1,function(yi){ mu_y + (yi)*sd_y }))
    # Compute criterions
    vars_expl[i_l] <- sum( (Y_train-Y_train_pred)^2 ) / sum( Y_train^2 )
    Q2[i_l] <- 1 - sum( (Y_test-Y_test_pred)^2 ) / sum((Y_test-Y_test_pred_RSS)^2)

    if(length(which(abs(V_svd)>1e-9))>0) model_exists[i_l] <- 1

    # RSS <- sum((Y_train_pred)^2)
    # V_model[,i_l] <- m_1$V_out
    # var_selected_y <- matrix(0,length(id_OOB),q)
    # id_y_selected <- which(colSums(abs(B))>1e-9)
    # if(length(id_y_selected)>0) model_exists[i_l] <- 1
    # var_selected_y[,id_y_selected] <- 1
    # y_pred <- t(apply(X_test_normalize,1,function(xi){
    #   mu_y + (xi%*%B)*sd_y
    # }))
    # var_selected_y_train <- matrix(0,length(id_IB),q)
    # var_selected_y_train[,id_y_selected] <- 1
    # y_train_est <- m_1$y_est
    # # vars_expl_star[i_l] <- sum(((y_train_est-Y_train)*var_selected_y_train)^2)/
    # #   sum(((Y_init[id_IB,,drop=F])*var_selected_y_train)^2)
    # vars_expl[i_l] <- sum(((y_train_est-Y_train))^2)/
    #   sum(((Y_init[id_IB,,drop=F]))^2)
    # y_0 <- t(apply(X_test_normalize,1,function(xi){
    #   mu_y
    # }))
    # PRESS_y <- (y_pred-Y_test)^2
    # RSS_y <- (y_0-Y_test)^2
    # # RSS_y_star <- RSS_y*var_selected_y
    # # Q2_star[i_l] <- 1-sum(var_selected_y*PRESS_y)/sum(RSS_y_star)
    # Q2[i_l] <- 1-sum(PRESS_y)/sum(RSS_y)
  }
  #vars_expl_star=vars_expl_star,Q2_star=Q2_star
  list(u_out=u_out,V_optim_phi=V_optim_phi,id_OOB=id_OOB,model_exists=model_exists,
       vars_expl=vars_expl,Q2=Q2)#,B=B_out,cov=COV_init)#PRESS=sum(PRESS_y_star),RSS=sum(RSS_y_star),

}


Q2_OLS <- function(x,y){
  X <- scale(x)
  Y <- scale(y)
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
  RSS <- mean(rowSums(Y^2))
  1-PRESS/RSS
}

Q2_OLS_th <- function(x,y,B_th){
  X <- scale(x)
  Y <- scale(y)
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
  RSS <- mean(rowSums(Y^2))
  1-PRESS/RSS
}


#' Title
#'
#' @param Xs Xs
#' @param Y Y
#' @param lambdas lambdas
#' @param tau rate learning
#' @param NCORES number of cores
#' @param NZV near zero var
#' @param verbose verbose
#'
#' @return
#' @export
#'
#' @useDynLib ddsPLS
Q2_local_ddsPLS <- function(Xs,Y,N_lambdas = 100,
                            lambda_max=1,alpha=2/3,
                            n_B=20,
                            deflatX=T,
                            tau=0.0975,NCORES=1,center=T,
                            NZV=1e-3,verbose=F){
  ###### INTERNAL FUNCTION
  do_soft_thresh_lowest <- function(x0,y0,lambda,cov_init=NULL){
    n <- nrow(y0)
    if(is.null(cov_init)){
      cov_init <- crossprod(y0,x0)/(n-1)
    }
    abs_COV <- abs(cov_init)
    id_next <- which(abs_COV<lambda)
    if(length(id_next)==0){
      lambda_out <- lambda
    }else{
      lambda_out <- (max(abs_COV[id_next])+lambda)/2
    }
    COV_sth <- abs_COV-lambda_out#lambda#
    COV_sth[which(COV_sth<0)] <- 0
    list(cov_th=COV_sth*sign(cov_init),lambda=lambda_out)#lambda)#
  }
  get_errors_stuff <- function(x0,y0,N_lambdas,lambdas_h,remove_COV=NULL,
                               deflatX,RSS0,NCORES){
    q <- ncol(y0)
    NCORES_w <- min(NCORES,N_lambdas)
    `%my_do%` <- ifelse(NCORES_w!=1,{
      out<-`%dopar%`
      cl <- makeCluster(NCORES_w)#cl <- parallel::makeCluster(NCORES_w)
      registerDoParallel(cl)#doParallel::registerDoParallel(cl)
      out},{
        out <- `%do%`
        out})
    PRESS_il_y_and_i_l <- foreach(i_l=1:N_lambdas,.packages = "ddsPLS") %my_do% {
      out <- do_loo(list(x0),y0,lam=lambdas_h[i_l])
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
    RSS_il_y = var_sel = RSS_il_y_star <- matrix(0,N_lambdas,q)
    vars <- rep(0,length(lambdas_h))
    for(i_l in 1:N_lambdas){
      lam <- lambdas_h[i_l]
      m_1 <- model_PLS(x0,y0,lam,remove_COV=remove_COV,deflatX=deflatX,to.scale = F)
      vars[i_l] <- sum((m_1$y_est)^2)/sum(RSS0)
      B_i_l <- m_1$B
      id_B <- which(colSums(abs(B_i_l))>1e-9)
      y_est_train <- m_1$y_est
      var_sel[i_l,id_B] <- 1
      RSS_il_y[i_l,] <- colSums((y_est_train-y0)^2)
      if(length(id_B)>0){
        RSS_il_y_star[i_l,id_B] <- colSums(((y_est_train-y0)^2)[,id_B,drop=F])
      }
    }
    Q2_h_y <- 1-PRESS_il_y/matrix(rep(RSS_h_moins_1,N_lambdas),nrow = N_lambdas,byrow = T)
    Q2_h_sum <- 1-rowSums(PRESS_il_y)/sum(RSS_h_moins_1)
    # STAR
    RSS_h_moins_1_star <- t(apply(var_sel,1,function(vv){vv*RSS_h_moins_1}))
    if(ncol(RSS_h_moins_1_star)!=q) RSS_h_moins_1_star <- t(RSS_h_moins_1_star)
    RSS_h_moins_1_star[which(is.na(RSS_h_moins_1_star))] <- 0
    # Q2_h_sum_star[[h]] <- 1-rowSums(PRESS_il_y_star)/rowSums(RSS_h_moins_1_star)
    list(Q2_h_star = 1-rowSums(PRESS_il_y_star)/rowSums(RSS_h_moins_1_star),
         vars=vars,RSS_il_y=RSS_il_y,var_sel=var_sel,Q2_h_y=Q2_h_y,
         Q2_h_sum=Q2_h_sum,RSS_il_y_star=RSS_il_y_star,
         RSS_h_moins_1_star=RSS_h_moins_1_star,PRESS_il_y=PRESS_il_y,
         PRESS_il_y_star=PRESS_il_y_star)
  }
  #######
  n <- nrow(Xs[[1]])
  q <- ncol(Y)
  id_na <- lapply(Xs,function(xx){which(is.na(xx[,1]))})
  K <- length(Xs)
  Y_init <- scale(Y,center = center)
  Xs_init <- lapply(Xs,scale,center=center)
  ps <- unlist(lapply(Xs_init,ncol))
  p <- sum(ps)
  mu_x_s <- lapply(Xs,colMeans)
  mu_y <- colMeans(Y)
  if(!center){
    mu_x_s <- lapply(mu_x_s,function(mu){mu*0})
    mu_y <- colMeans(Y)*0
  }
  MU_X <- matrix(rep(unlist(mu_x_s),n),nrow = n,byrow = T)
  MU_Y <- matrix(rep(mu_y,n),nrow = n,byrow = T)
  SD_Y <- matrix(rep(apply(Y,2,sd),n),ncol = q,byrow = T)
  SD_X <- matrix(rep(apply(do.call(cbind,Xs),2,sd),n),ncol = p,byrow = T)
  sd_x_inv_mat <- matrix(rep(unlist(lapply(apply(do.call(cbind,Xs),2,sd),function(ss){if(abs(ss)>1e-9){out <- 1/ss}else{out <- 0};out})),q),ncol = q,byrow = T)
  sd_y_mat <- matrix(rep(apply(Y,2,sd),sum(ps)),ncol = q,byrow = T)
  sd_y_inv_mat <- matrix(rep(unlist(lapply(apply(Y,2,sd),function(ss){if(abs(ss)>1e-9){out <- 1/ss}else{out <- 0};out})),sum(ps)),ncol = q,byrow = T)
  sd_x_mat <- matrix(rep(apply(do.call(cbind,Xs),2,sd),q),ncol = q,byrow = T)
  sd_y_x_inv <- sd_y_mat * sd_x_inv_mat
  sd_x_y_inv <- sd_x_mat * sd_y_inv_mat
  RSS0 = RSS_h_moins_1 = RSS_h_moins_1_star <- colSums(Y_init^2)
  PRESS_y = RSS_y = Q2_y <- matrix(NA,n,q)
  PRESS_y_star = RSS_y_star = Q2_y_star <- matrix(NA,n,q)
  lambda = Q2_sum_star <- rep(NA,n)
  ps <- unlist(lapply(Xs_init,ncol))
  X_init <- do.call(cbind,Xs_init)
  x0 <- X_init
  u <- matrix(NA,sum(ps),n)
  V_phi = V <- matrix(NA,q,n)
  P <- matrix(0,sum(ps),n)
  B <- matrix(0,sum(ps),q)
  t_h <- matrix(NA,n,n)
  Us = Bs = B_r_out =
    Q2_h_sum_star = Q2_h_sum_star_sd_moins = Q2_h_sum_star_sd_plus =
    lambdas_out = vars_h_boot = vars_h_boot_sd_moins = vars_h_boot_sd_plus <- list()
  q2_max_h <- NULL
  B_tot_LOO <- matrix(0,n,sum(ps)*q)

  ### For each 'h' component, look for the best lambda
  test <- T
  y0 <- Y_init
  h <- 0
  VAR_h_s <- NULL
  K_h <- 1:K
  # Check lambdas and stuff
  lambdas_h <- seq(0,lambda_max,length.out = N_lambdas)
  ncomps = Q2_TOTAL <- rep(NA,N_lambdas)
  x0_deflated = y0_deflated = cov_mean_h = prop_models_ok <- list()
  remove_COV <- matrix(0,q,p)
  B_previous <- matrix(0,p,q)
  while(test){
    h <- h + 1
    Bis <- list()
    pos_decoupe <- 1
    Q2_Boot = vars = Q2_h_sum <- rep(0,N_lambdas)
    RSS_il_y = var_sel = Q2_h_y = RSS_il_y_star = RSS_h_moins_1_star = PRESS_il_y <- matrix(0,N_lambdas,q)
    Q2_B = vars_in = Q2_h_clas <- matrix(NA,n_B,N_lambdas)
    ## LOO
    # result_b <- get_errors_stuff(x0=x0,y0=y0,N_lambdas=N_lambdas,lambdas_h=lambdas_h,
    #                              remove_COV=remove_COV,deflatX=deflatX,RSS0=RSS0,
    #                              NCORES=NCORES)
    # RSS_il_y_star  <- result_b$RSS_il_y_star
    # PRESS_il_y_star <- result_b$PRESS_il_y_star
    # RSS_h_moins_1_star  <- result_b$RSS_h_moins_1_star
    # # Q2_Boot <- result_b$Q2_h_star
    # vars <-  result_b$vars
    # RSS_il_y <- result_b$RSS_il_y
    # var_sel <- result_b$var_sel
    # Q2_h_y <- result_b$Q2_h_y
    # Q2_h_sum <- result_b$Q2_h_sum
    # PRESS_il_y <- result_b$PRESS_il_y
    ## BOOTSTRAP
    NCORES_w <- min(NCORES,n_B)
    `%my_do%` <- ifelse(NCORES_w!=1,{
      out<-`%dopar%`;cl <- makeCluster(NCORES_w)
      registerDoParallel(cl);out},{out <- `%do%`;out})
    Q2_star_bootstrat <- foreach(i_B=1:n_B,.packages = "ddsPLS",.combine='c',.multicombine=TRUE) %my_do% {
      out <- bootstrap_pls(X_init=X_init,Y_init=Y_init,B_previous=B_previous,u=u,v=V_phi,h=h,lambdas=lambdas_h)
      res_measure <- cbind(1:length(lambdas_h),out$Q2,#out$Q2_star,
                           out$vars_expl,# out$vars_expl_star,
                           length(out$id_OOB),
                           i_B,out$model_exists)
      list(res=res_measure,V_optim_phi=out$V_optim_phi,u=out$u_out)#list(cov=out$cov,B=out$B,res=res_measure,V_model=out$V_model,u=out$u_out,V=out$V_out)
    }
    if(NCORES_w!=1)stopCluster(cl)
    KK <- length(Q2_star_bootstrat)/n_B ; id_pos_boot <- (1:n_B)*KK
    # temp <- array(unlist(Q2_star_bootstrat[id_pos_boot-5]), c(q, p, n_B))
    # cov_mean <- apply(temp, 1:2, median)
    # B_boot <- Q2_star_bootstrat[id_pos_boot-4]
    res_measure <- do.call(rbind,Q2_star_bootstrat[id_pos_boot-2])
    V_boot <- Q2_star_bootstrat[id_pos_boot-1]
    u_boot <- Q2_star_bootstrat[id_pos_boot]
    # V_boot <- Q2_star_bootstrat[id_pos_boot]
    vars_boot = vars_boot_sd_plus = vars_boot_sd_moins =
      q2_boot = q2_boot_sd_plus = q2_boot_sd_moins =
      prop_models_ok[[h]] <- rep(0,N_lambdas)
    mod_exists <- rep(0,N_lambdas)
    test_batchas <- TRUE

    i_l <- 1
    while(test_batchas){
      pos <- which(res_measure[,1]==i_l)
      toto <- na.omit(res_measure[pos,c(2,3,6)])
      if(length(toto)>0){#!any(is.na(toto[,1])) & !any(is.nan(toto[,1]))){
        Q2_qtles <- quantile(toto[,1])
        q2_boot[i_l] <- Q2_qtles[3]
        q2_boot_sd_plus[i_l] <- Q2_qtles[4]
        q2_boot_sd_moins[i_l] <- Q2_qtles[2]
        VARS_qtles <- quantile(toto[,2])
        vars_boot[i_l] <- VARS_qtles[3]
        vars_boot_sd_plus[i_l] <- VARS_qtles[4]
        vars_boot_sd_moins[i_l] <- VARS_qtles[2]
        prop_models_ok[[h]][i_l] <- sum(na.omit(toto[,3]))/length(pos)
        if(prop_models_ok[[h]][i_l]>=alpha){
          mod_exists[i_l] <- 1#floor(sum(toto[,3])/length(toto[,3]))
        }
      }
      if(is.na(q2_boot[i_l])|is.nan(q2_boot[i_l])|i_l==N_lambdas){#|mod_exists[i_l]!=1){
        test_batchas <- F
      }
      i_l <- i_l + 1
    }
    vars <- vars_boot
    id_vars_loo_ok <- which(vars_boot_sd_moins>NZV & mod_exists==1)

    Q2_h_sum_star[[h]] <- q2_boot#result_b$Q2_h_star
    Q2_h_sum_star_sd_moins[[h]] <- q2_boot_sd_moins
    Q2_h_sum_star_sd_plus[[h]] <- q2_boot_sd_plus
    lambdas_out[[h]] <- lambdas_h
    vars_h_boot[[h]] <- vars_boot
    # vars_h_boot_sd[[h]] <- vars_boot_sd
    vars_h_boot_sd_moins[[h]] <- vars_boot_sd_moins
    vars_h_boot_sd_plus[[h]] <- vars_boot_sd_plus
    # cov_mean_h[[h]] <- cov_mean
    if(length(id_vars_loo_ok)>0){
      q2_max_h <- c(q2_max_h,max(na.omit(q2_boot[id_vars_loo_ok])))
      id_s_cool <- which(q2_boot==q2_max_h[h] & q2_boot>tau &#q2_boot_sd_plus>=q2_max_h[h] & q2_boot>tau &#>=q2_max_h[h] & q2_boot>tau &#
                           vars_boot>NZV & mod_exists==1)
      if(length(id_s_cool)>0){
        best_id_h <- max(1,min(id_s_cool)-1)
        # coco_bobo <- apply(abs(cov_mean),2,max)
        # mean_coco_bobo_best <- mean(coco_bobo[order(abs(coco_bobo-lambdas_h[best_id_h]))[1:2]])
        # best_id_h <- which.min(abs(mean_coco_bobo_best-lambdas_h))
        ####
        # I <- which(q2_boot==q2_max_h[h])
        # m_I <- model_PLS(x0,y0,lambdas_h[I],R = 1,remove_COV=remove_COV,COV_init = cov_mean,
        #                  NZV=NZV,deflatX=deflatX,to.scale = F)
        # u_I <- m_I$U_star
        # t_I <- x0%*%u_I
        # phi_I <- t_I/sqrt(sum(t_I^2))
        # cor_x_tI <- abs(crossprod(phi_I,x0))/sqrt(n-1)
        # id_I <- which(abs(u_I)>1e-9);cor_min <- min(cor_x_tI[id_I])
        # id_all <- which(cor_x_tI>=cor_min)
        # cor_x_y <- abs(crossprod(y0,x0)/(n-1))
        # lambda_gogo <- min(apply(cor_x_y[,id_all,drop=F],2,max))
        # lambda_sol <- max(lambdas_h[which(lambdas_h<lambda_gogo)])
        # best_id_h <- which(lambdas_h==lambda_sol)
        ####
        if(length(best_id_h)>0){
          test_h <- Q2_h_sum_star[[h]][best_id_h]>tau
          best_lam_h <- lambdas_h[best_id_h]
          ### Find best solution
          # m_plus <- model_PLS(x0,y0,best_lam_h,R = 1,remove_COV=remove_COV,
          #                     NZV=NZV,RSS0 = RSS0,
          #                     deflatX=deflatX,to.scale = F)
          ## The regression matrix is the result of projection of x0 on bootstrap u
          u_boot_best <- lapply(u_boot,function(u_l){u_l[,best_id_h]})
          u_sol <- apply(do.call(rbind,u_boot_best),2,median)
          u_sol_boot_h <- u_sol/sqrt(sum(u_sol^2))
          t_boot_h <- x0%*%u_sol_boot_h
          V_boot_best <- lapply(V_boot,function(v_l){v_l[,best_id_h]})
          V_optim_boot_h <- apply(do.call(rbind,V_boot_best),2,median)
          V_optim_boot_h <- V_optim_boot_h/sqrt(sum(V_optim_boot_h^2))
          V_phi[,h] <- V_optim_boot_h
          V_sol_boot_h <- V_optim_boot_h%*%crossprod(V_optim_boot_h,crossprod(y0,t_boot_h))/sum(t_boot_h^2)
          B_boot_h <- tcrossprod(u_sol_boot_h,V_sol_boot_h)
          ## The regression matrix is the result of projection of x0 on bootstrap u and v
          # V_model_boot_best <- lapply(V_model_boot,function(v_l){v_l[,best_id_h]})
          # V_model_boot_h <- t(t(apply(do.call(rbind,V_model_boot_best),2,median)))
          # B_boot_second_h <- tcrossprod(u_sol_boot_h,V_model_boot_h)
          # ## The regression matrix is the result of bootstrapping --> the most general
          # B_boot_h_0 <- lapply(B_boot,function(v_l){v_l[,best_id_h]})
          # B_boot_h <- matrix(apply(do.call(rbind,B_boot_h_0),2,median),ncol=q,byrow = F)
          B_previous <- B_previous + B_boot_h
          ### End
          y_est_boot <- x0%*%B_boot_h
          varia_expl <- sum(y_est_boot^2)/sum(RSS0)
          variance_h <- varia_expl#m_plus$var_expl
          if(variance_h>NZV){
            tested_bad_comp <- F
            if(T){
              if(h>1){
                if(F){#round(min(VAR_h_s),2)+NZV<round(variance_h,2) ){
                  vars_init <- VAR_h_s
                  h <- min(which(round(VAR_h_s,2)<round(variance_h,2)))-1
                  if(h>0){
                    x0 <- x0_deflated[[h]]
                    y0 <- y0_deflated[[h]]
                    VAR_h_s <- VAR_h_s[1:h]
                  }else{
                    x0 <- do.call(cbind,Xs_init)
                    y0 <- Y_init
                    VAR_h_s <- NULL
                  }
                  res_remove <- do_soft_thresh_lowest(x0,y0,lambda=lambda[h+1])#,cov_init=cov_mean)
                  # u_x_current <- u_sol_boot_h#m_plus$U_star
                  # v_y_current <- V_model_boot_h#m_plus$V_out
                  # v_y_current <- v_y_current/sqrt(sum(v_y_current^2))
                  cov_remove_current <- res_remove$cov_th
                  nono <- sqrt(sum(B_boot_h^2))
                  gaga <- t(B_boot_h/nono)
                  alpha_current <- sum(cov_remove_current*gaga)#as.numeric(crossprod(v_y_current,cov_remove_current%*%u_x_current))
                  # b_current <- tcrossprod(v_y_current,u_x_current)
                  # b_current <- tcrossprod(m_plus$V_out,m_plus$U_star)
                  # cov_remove_current[which(abs(b_current)>1e-9)] <- 0
                  remove_COV <- cov_remove_current - alpha_current*gaga#remove_COV + cov_remove_current - alpha_current*b_current
                  if(verbose){
                    cat(paste("\n                 Conflict (",paste(round(vars_init,2),collapse = ","),
                              ")   against   ",round(variance_h,2),sep=""))
                    cat(paste("\n                 Norm of remove_cov ",
                              round(norm(remove_COV),5),sep="" ))
                    cat(paste("\n                 Next limit of lambda ",
                              round(res_remove$lambda,5),"\n",sep="" ))
                  }
                  lambdas_h <- seq(min(lambdas_h),res_remove$lambda,length.out = N_lambdas)#
                  tested_bad_comp <- T
                }else{
                  remove_COV <- 0*remove_COV
                  VAR_h_s <- c(VAR_h_s,variance_h)
                }
              }else{
                VAR_h_s <- c(VAR_h_s,variance_h)
              }
            }
            if(!tested_bad_comp){
              #Q2_h_sum_star[[h]][best_id_h]>tau
              if(!test_h){
                test <- F
              }else{
                lambda[h] <- best_lam_h
                # RSS_y[h,] <- RSS_il_y[best_id_h,]
                # if(h==1){
                #   RSS0_star <- RSS_h_moins_1_star[best_id_h,]
                # }
                # RSS_h_moins_1 <- RSS_y[h,]
                # PRESS_y[h,] <- PRESS_il_y[best_id_h,]
                # Q2_y[h,] <- Q2_h_y[best_id_h,]
                # # STAR
                # RSS_y_star[h,] <- RSS_il_y_star[best_id_h,]
                # RSS_h_moins_1_star <- RSS_y_star[h,]
                # PRESS_y_star[h,] <- PRESS_il_y_star[best_id_h,]
                # Q2_sum_star[h] <- Q2_h_sum_star[[h]][best_id_h]
                # Get the regression matrix of the optimal model
                u[,h] <- u_sol_boot_h#m_plus$U_star
                if(h!=1){
                  for(s_r in (h-1):1){
                    if(s_r!=1){
                      x_s_r <- x0_deflated[[s_r-1]]
                    }else{
                      x_s_r <- do.call(cbind,Xs_init)
                    }
                    u_s_r <- u[,s_r,drop=F]
                    t_s_r <- x_s_r%*%u_s_r
                    bt <- crossprod(t_s_r,x_s_r)/sum(t_s_r^2)
                    u[,h] <- u[,h]-u_s_r*sum(bt*u[,h])
                  }
                }
                t_r <- x0%*%u_sol_boot_h#m_plus$U_out
                t_h[,h] <- t_r
              }
              B_r <- B_boot_h#m_plus$B
              V_r <- V_sol_boot_h#m_plus$V_out
              # y0_plus <- y0 - m_plus$e_y
              # id_y_sel <- which(abs(V_r)>1e-9)
              if(verbose){
                cat(paste("\nComponent ",h,"   ",
                          "   lambda=",round(best_lam_h,3),
                          "   var.expl.=",round(variance_h*100),"% \n",sep=""))
              }
              B_r_out[[h]] <- B_r*sd_y_x_inv
              V[,h] <- V_r
              y0 <- y0-x0%*%B_boot_h#m_plus$e_y # y0 - y0_plus
              y0_deflated[[h]] <- y0
              if(deflatX){
                if(is.null(t_r)){
                  t_r <- x0%*%u_sol_boot_h
                }
                bt <- crossprod(t_r,x0)/sum(t_r^2)
                x0 <- x0 - t_r%*%bt
                x0_deflated[[h]] <- x0
              }
              lambdas_h <- seq(min(lambdas_h),best_lam_h,length.out = length(lambdas_h))
              lambdas_h <- unique(lambdas_h)
              N_lambdas <- length(lambdas_h)
            }
          }else{
            test <- F
          }
        }
      }else{
        test <- F
      }
    }else{
      test <- F
    }
  }
  h_opt <- h - 1
  x0_center <- scale(do.call(cbind,Xs_init),scale = F,center=center)
  if(h_opt>0){
    Q2_cum_star <- 1-prod(1-Q2_sum_star[1:h_opt])
    for(h in 1:h_opt){
      B <- B + B_r_out[[h]]
      # if(h==1){
      #   Q2_cum <- 1- sum(PRESS_y[h,])/sum(RSS0)
      #   Q2_cum_y <- Q2_y[h,]
      # }else{
      #   Q2_cum <- 1-(1-Q2_cum)*sum(PRESS_y[h,])/sum(RSS_y[h-1,])
      #   Q2_cum_y <- 1-(1-Q2_cum_y)*PRESS_y[h,]/RSS_y[h-1,]
      # }
    }
    i_0 <- 0
    for(k in 1:K){
      Us[[k]] <- u[i_0+1:ps[k],1:h_opt,drop=F]
      Bs[[k]] <- B[i_0+1:ps[k],,drop=F]
      i_0 <- i_0+ps[k]
    }
    lambda_sol <- lambda[1:h_opt]
    mu_y <- matrix(rep(colMeans(Y),n) ,nrow = n,byrow = T)
    y_est <- mu_y + x0_center%*%(B*sd_x_mat)
    # # Compute LOO error to get Q2_reg
    # Y_pred_all <- matrix(0,n,q)
    # p <- ncol(x0_center)
    # var_sel <- matrix(0,n,q)
    # Xs_cbind <- as.matrix(do.call(cbind,Xs_init))
    # ## Get LOO error
    # for(ii in 1:n){
    #   X_train <- Xs_cbind[-ii,,drop=F]
    #   X_test <- Xs_cbind[ii,,drop=F]
    #   Y_train <- Y_init[-ii,,drop=FALSE]
    #   Y_test <- Y_init[ii,,drop=FALSE]
    #   mu_k_ii <- colMeans(X_train)
    #   mu_y_ii <- colMeans(Y_train)
    #   sd_x_inv_ii <- unlist(lapply(apply(X_train,2,sd),function(ss){if(abs(ss)>1e-9){out <- 1/ss}else{out <- 0};out}))
    #   # sd_x_mat_ii <- apply(X_train,2,sd)
    #   sd_y_mat_ii <- apply(Y,2,sd)
    #   # sd_y_inv_ii <- unlist(lapply(apply(Y,2,sd),function(ss){if(abs(ss)>1e-9){out <- 1/ss}else{out <- 0};out}))
    #   m_plus <- model_PLS(X_train,Y_train,lambda_sol,remove_COV=remove_COV,
    #                       R = h_opt,NZV=NZV,deflatX=deflatX)
    #   var_sel[ii,which(rowSums(abs(m_plus$V_out))>1e-9)] <- 1
    #   Y_pred_all[ii,] <- mu_y_ii + ((X_test-mu_k_ii))%*%m_plus$B
    # }
    # err_LOO <- (Y_init-Y_pred_all)
    # ERRORS_LOO <- colSums(err_LOO^2)
    # RSSO_loo <- colSums(scale(Y_init,center=center)^2)
    # Q2_reg <- 1 - sum(ERRORS_LOO)/sum(RSSO_loo)
    # Q2_reg_y <- 1 - ERRORS_LOO/RSSO_loo
    # # STAR
    # X_binded <- do.call(cbind,Xs_init)
    # m_star <- model_PLS(X_binded,Y,lambda_sol,R = h_opt,remove_COV=remove_COV,
    #                     NZV=NZV,deflatX=deflatX)
    # ERRORS_LOO_star <- colSums((err_LOO^2)*var_sel)
    # RSS0_star_model <- RSSO_loo[which(rowSums(abs(m_star$V_out))>1e-9)]
    # Q2_reg_star <- 1 - sum(ERRORS_LOO_star)/sum(RSS0_star_model)
    # Prepare outputs
    #
    optimal_parameters <- list(lambda=lambda_sol,R=h_opt)
    parameters <- list(RSS0=RSS0)#,RSS0_star=RSS0_star,RSS_y=RSS_y[1:(h_opt+1),],PRESS_y=PRESS_y[1:(h_opt+1),],Q2_y=Q2_y[1:(h_opt+1),])
    bootstrap <- list(lambdas_h=lambdas_out,
                      q2_max_h=q2_max_h,
                      Q2_h_star=Q2_h_sum_star,
                      Q2_h_star_sd_plus=Q2_h_sum_star_sd_plus,
                      Q2_h_star_sd_moins=Q2_h_sum_star_sd_moins,
                      vars_h_boot=vars_h_boot,
                      vars_h_boot_sd_plus=vars_h_boot_sd_plus,
                      vars_h_boot_sd_moins=vars_h_boot_sd_moins,
                      cov_mean_h=cov_mean_h,prop_models_ok=prop_models_ok)
    res <- list(optimal_parameters=optimal_parameters,
                bootstrap=bootstrap,
                Us=Us,V=V[,1:h_opt,drop=F],B_cbind=B,
                x0_deflated=x0_deflated,y0_deflated=y0_deflated,
                t_h=t_h[,1:h_opt,drop=F],
                explained_variance=VAR_h_s*100,
                NZV=NZV,tau=tau,#score_X = m_plus$score_x,
                Bs=Bs,B_r=B_r_out,
                y_est=y_est,parameters=parameters,
                mu_k=c(mu_x_s),mu_y=mu_y,sd_y=apply(Y,2,sd))
    if(verbose){
      plot_res(res)
    }
  }else{
    res <- NULL
  }
  res
}

plot_res <- function(res){
  h_opt <- res$optimal_parameters$R
  lambdas_out <- res$bootstrap$lambdas_h
  Q2_h_sum_star <- res$bootstrap$Q2_h_star
  Q2_h_star_sd_moins <- res$bootstrap$Q2_h_star_sd_moins
  Q2_h_star_sd_plus <- res$bootstrap$Q2_h_star_sd_plus
  vars_h_boot <- res$bootstrap$vars_h_boot
  vars_h_boot_sd_moins <- res$bootstrap$vars_h_boot_sd_moins
  vars_h_boot_sd_plus <- res$bootstrap$vars_h_boot_sd_plus
  par(mfrow=c(1,h_opt))
  for(h in 1:h_opt){
    # coco_mean <- apply(abs(res$bootstrap$cov_mean_h[[h]]),2,max)
    plot(lambdas_out[[h]],Q2_h_sum_star[[h]],type="l",ylim=c(0,1.1),col="red",lwd=2,
         xlab=expression(lambda),ylab="",main=paste("Component ",h," (",
                                                    round(res$explained_variance[h]),"%)",sep="",collapse = ""))
    points(lambdas_out[[h]],Q2_h_star_sd_plus[[h]],lwd=1.5,
           lty=2,type="l",col="red")
    points(lambdas_out[[h]],Q2_h_star_sd_moins[[h]],lwd=1.5,
           lty=2,type="l",col="red")
    points(lambdas_out[[h]],vars_h_boot[[h]],type="l",col="blue")
    points(lambdas_out[[h]],vars_h_boot_sd_plus[[h]],
           lty=2,type="l",col="blue")
    points(lambdas_out[[h]],vars_h_boot_sd_moins[[h]],
           lty=2,type="l",col="blue")
    popo <- mean(range(lambdas_out[[h]]))
    max_h <- res$bootstrap$q2_max_h[h]
    points(c(res$optimal_parameters$lambda[h],res$optimal_parameters$lambda[h]),c(max_h,-2),
           type="l",col="black")
    points(res$optimal_parameters$lambda[h],max_h,pch=16,col="red",cex=1.5)
    text(x=popo,y=max_h,pos=3,labels = bquote("||Q"[B]^"2"~"||"[infinity]==.(round(max_h,2))),col='red',cex=1)
    abline(h=res$tau,lty=1,col="red",lwd=0.5)
    text(x=popo,y=res$tau,pos=3,labels = bquote(tau==.(res$tau)),col='red',cex=1)
    # abline(v=,lwd=1.5)
    # abline(v=coco_mean,col='gray',lwd=0.5,lty=1)
    legend("topleft",c(expression("Q"[B]^"2"),"Expl.var.Y"),
           fill=c("red","blue"),bty="n")
    legend("top",c("Median","25% / 75%"),lty=c(1,2),bty="n")
    points(lambdas_out[[h]],res$bootstrap$prop_models_ok[[h]],col="brown",type="l",lty=1)
  }
}
