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

      svd_mix_Y <- svd(COV_COV_high,nu = 1,nv = 1)#tcrossprod(COV_COV_high,COV_high),nu = 0,nv = 1)
      # svd_mix_X <- svd(COV_COV_high,nu = 0,nv = 1)#crossprod(COV_COV_high,COV_high),nu = 0,nv = 1)
      # u_x_no_std <- t(COV_COV_high)%*%model_NULL$u
      U0[id_x_high,] <- svd_mix_Y$v#svd_mix_X$v#u_x_no_std/sqrt(sum(u_x_no_std^2))#svd_XY$v#
      ## Y part
      # u_y_no_std <- COV_COV_high%*%model_NULL$v
      V0[id_y_high,] <- svd_mix_Y$u#u_y_no_std/sqrt(sum(u_y_no_std^2))#svd_XY$u#
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
    x_init <- scale(x)
    y_init <- scale(y)
  }else{
    x_init <- scale(x,scale = F)
    y_init <- scale(y,scale = F)
    mu_y <- matrix(0,n,q)
    sd_y_mat <- matrix(rep(apply(y,2,sd),p),ncol = q,byrow = T)
  }
  x0 <- x_init
  y0 <- y_init
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
        # var_y_plus_un <-
        var_expl[r] <- 1-sum((y_init-mu_y-y_plus_un)^2)/sum((y_init-mu_y)^2)
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
  vars_expl_star = vars_expl = vars_expl_h = Q2_star = Q2 = Q2_all = model_exists <- rep(0,N_lambdas)
  B_out <- matrix(0,q*p,N_lambdas)

  u_out <- matrix(0,p,N_lambdas)
  V_optim_phi = V_model <- matrix(0,q,N_lambdas)
  B_next_out <- list()
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
    B_all <- B_youyou + B_next
    y_train_pred <- X_train%*%B_all
    y_test_pred <- X_test_normalize%*%B_all
    # Y_train_pred <- t(apply(y_train_pred,1,function(yi){ mu_y + (yi)*sd_y }))
    # Y_test_pred <- t(apply(y_test_pred,1,function(yi){ mu_y + (yi)*sd_y }))
    # Previous components
    y_train_pred_next <- X_train%*%B_next#B_previous
    y_test_pred_RSS <- X_test_normalize%*%B_youyou#B_previous
    # Y_test_pred_RSS <- t(apply(y_test_pred_RSS,1,function(yi){ mu_y + (yi)*sd_y }))
    # Compute criterions
    vars_expl[i_l] <- 1 - sum( (Y_train-y_train_pred)^2 ) / sum( (Y_train)^2 )
    vars_expl_h[i_l] <- 1 - sum( (Y_train-y_train_pred_next)^2 ) / sum((Y_train)^2)
    Q2[i_l] <- 1 - sum( (Y_test-y_test_pred)^2 ) / sum((Y_test-y_test_pred_RSS)^2)
    MU_test <- matrix(rep(mu_y,length(id_OOB)),nrow = length(id_OOB),byrow = T)
    Q2_all[i_l] <- 1 - sum( (Y_test-y_test_pred)^2 ) / sum((Y_test)^2)#-MU_test)^2)
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
       vars_expl=vars_expl,vars_expl_h=vars_expl_h,Q2=Q2,Q2_all=Q2_all)#,B=B_out,cov=COV_init)#PRESS=sum(PRESS_y_star),RSS=sum(RSS_y_star),

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
#' @param NCORES number of cores
#' @param NZV near zero var
#' @param verbose verbose
#'
#' @return
#' @export
#'
#' @useDynLib ddsPLS
Q2_local_ddsPLS <- function(Xs,Y,N_lambdas = 100,
                            lambda_max=1,alpha=1/3,
                            n_B=20,
                            deflatX=T,NCORES=1,center=T,
                            NZV=1e-2,verbose=F){
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
    Q2_all_sum_star = Q2_all_sum_star_sd_moins = Q2_all_sum_star_sd_plus =
    lambdas_out = vars_h_boot = vars_h_boot_sd_moins = vars_h_boot_sd_plus =
    vars_h_boot_single = vars_h_boot_single_sd_moins = vars_h_boot_single_sd_plus <- list()
  q2_max_h = VAR_h_s = CUM_VAR_h_s = Q2_tot_s <- rep(NA,n)
  B_tot_LOO <- matrix(0,n,sum(ps)*q)

  ### For each 'h' component, look for the best lambda
  test <- T
  y0 <- Y_init
  h <- 0
  id_ALL_TEST_h = Q2_all_sum_star_boxplot = R2_all_sum_star_boxplot <- list()
  K_h <- 1:K
  # Check lambdas and stuff
  lambdas_h <- seq(0,lambda_max,length.out = N_lambdas)
  ncomps = Q2_TOTAL <- rep(NA,N_lambdas)
  x0_deflated = y0_deflated = cov_mean_h = prop_models_ok <- list()
  remove_COV <- matrix(0,q,p)
  B_previous <- matrix(0,p,q)
  V_boot = u_boot = res_measure <- list()
  while(test){
    h <- h + 1
    Bis <- list()
    pos_decoupe <- 1
    Q2_Boot = vars = Q2_h_sum <- rep(0,N_lambdas)
    RSS_il_y = var_sel = Q2_h_y = RSS_il_y_star = RSS_h_moins_1_star = PRESS_il_y <- matrix(0,N_lambdas,q)
    Q2_B = vars_in = Q2_h_clas <- matrix(NA,n_B,N_lambdas)
    NCORES_w <- min(NCORES,n_B)
    `%my_do%` <- ifelse(NCORES_w!=1,{
      out<-`%dopar%`;cl <- makeCluster(NCORES_w)
      registerDoParallel(cl);out},{out <- `%do%`;out})
    Q2_star_bootstrat <- foreach(i_B=1:n_B,.packages = "ddsPLS",.combine='c',.multicombine=TRUE) %my_do% {
      out <- bootstrap_pls(X_init=X_init,Y_init=Y_init,B_previous=B_previous,u=u,v=V_phi,h=h,lambdas=lambdas_h)
      res_measure <- cbind(1:length(lambdas_h),
                           out$Q2,
                           out$Q2_all,
                           out$vars_expl,
                           out$vars_expl_h,
                           length(out$id_OOB),
                           i_B,
                           out$model_exists)
      list(res=res_measure,V_optim_phi=out$V_optim_phi,u=out$u_out)#list(cov=out$cov,,res=res_measure,V_model=out$V_model,u=out$u_out,V=out$V_out)
    }
    if(NCORES_w!=1)stopCluster(cl)
    KK <- length(Q2_star_bootstrat)/n_B ; id_pos_boot <- (1:n_B)*KK
    # temp <- array(unlist(Q2_star_bootstrat[id_pos_boot-5]), c(q, p, n_B))
    # cov_mean <- apply(temp, 1:2, median)
    # B_boot <- Q2_star_bootstrat[id_pos_boot-3]
    res_measure[[h]] <- do.call(rbind,Q2_star_bootstrat[id_pos_boot-2])
    V_boot[[h]] <- Q2_star_bootstrat[id_pos_boot-1]
    u_boot[[h]] <- Q2_star_bootstrat[id_pos_boot]
    # V_boot <- Q2_star_bootstrat[id_pos_boot]
    vars_boot = vars_boot_sd_plus = vars_boot_sd_moins =
      vars_boot_h = vars_boot_h_sd_plus = vars_boot_h_sd_moins =
      q2_boot = q2_boot_sd_plus = q2_boot_sd_moins =
      q2_all_boot = q2_all_boot_sd_plus = q2_all_boot_sd_moins =
      prop_models_ok[[h]] <- rep(0,N_lambdas)
    mod_exists <- rep(0,N_lambdas)
    test_batchas <- TRUE

    i_l <- 1
    aa <- min(alpha,1-alpha)
    probabilities <- c(aa,0.5,1-aa)
    while(test_batchas){
      pos <- which(res_measure[[h]][,1]==i_l)
      toto <- na.omit(res_measure[[h]][pos,c(2,3,4,5,8)])
      if(length(toto)>0){
        Q2_qtles <- quantile(toto[,1],probs = probabilities)
        q2_boot[i_l] <- Q2_qtles[2]
        q2_boot_sd_plus[i_l] <- Q2_qtles[3]
        q2_boot_sd_moins[i_l] <- Q2_qtles[1]
        Q2_all_qtles <- quantile(toto[,2],probs = probabilities)
        q2_all_boot[i_l] <- Q2_all_qtles[2]
        q2_all_boot_sd_plus[i_l] <- Q2_all_qtles[3]
        q2_all_boot_sd_moins[i_l] <- Q2_all_qtles[1]
        VARS_qtles <- quantile(toto[,3],probs = probabilities)
        vars_boot[i_l] <- VARS_qtles[2]
        vars_boot_sd_plus[i_l] <- VARS_qtles[3]
        vars_boot_sd_moins[i_l] <- VARS_qtles[1]
        VARS_qtles <- quantile(toto[,4],probs = probabilities)
        vars_boot_h[i_l] <- VARS_qtles[2]
        vars_boot_h_sd_plus[i_l] <- VARS_qtles[3]
        vars_boot_h_sd_moins[i_l] <- VARS_qtles[1]
        prop_models_ok[[h]][i_l] <- sum(na.omit(toto[,5]))/length(pos)
        # if(prop_models_ok[[h]][i_l]>=alpha){
        #   mod_exists[i_l] <- 1#floor(sum(toto[,3])/length(toto[,3]))
        # }
      }
      if(is.na(q2_boot[i_l])|is.nan(q2_boot[i_l])|i_l==N_lambdas){#|mod_exists[i_l]!=1){
        test_batchas <- F
      }
      i_l <- i_l + 1
    }
    vars <- vars_boot
    # id_vars_loo_ok <- which(vars_boot_sd_moins>0)# & mod_exists==1)
    id_ALL_TEST <- which(q2_boot_sd_moins > 0)
    if(h!=1){
      # q2_all_boot_sd_moins > Q2_all_sum_star_sd_plus[[h-1]]
      # l_s_before <- lambdas_out[[h-1]]
      # l_s_now <- lambdas_h
      # l_s_unif <- sort(unique(c(l_s_now,l_s_before[which(l_s_before<=max(l_s_now))])))
      # N_unif <- length(l_s_unif)
      # Vec_before = Vec_now <- rep(NA,N_unif)
      # for(i_l in 1:N_unif){
      #   pos_before <- which(l_s_before==l_s_unif[i_l])
      #   pos_now <- which(l_s_now==l_s_unif[i_l])
      #   if(length(pos_now)>0){
      #     Vec_now[i_l] <- q2_all_boot_sd_moins[pos_now]
      #   }else{
      #     i_last <- max(which(!is.na(Vec_now[1:(i_l-1)])))
      #     Vec_now[i_l] <- Vec_now[i_last]
      #   }
      #   if(length(pos_before)>0){
      #     Vec_before[i_l] <- Q2_all_sum_star_sd_plus[[h-1]][pos_before]
      #   }else{
      #     i_last <- max(which(!is.na(Vec_before[1:(i_l-1)])))
      #     Vec_before[i_l] <- Vec_before[i_last]
      #   }
      # }
      # pos_ok <- which(Vec_before<Vec_now)
      # lambdas_ok <- l_s_unif[pos_ok]
      # id_test_h <- which(lambdas_h %in% lambdas_ok)
      # id_ALL_TEST <- intersect(id_ALL_TEST,id_test_h)
      id_test_h <- which(q2_all_boot_sd_moins>Q2_all_sum_star_sd_plus[[h-1]][
        which(lambdas_out[[h-1]]==lambda[h-1])
      ])
      id_ALL_TEST <- intersect(id_ALL_TEST,id_test_h)
    }
    id_ALL_TEST_h[[h]] <- id_ALL_TEST
    Q2_h_sum_star[[h]] <- q2_boot#result_b$Q2_h_star
    Q2_h_sum_star_sd_moins[[h]] <- q2_boot_sd_moins
    Q2_h_sum_star_sd_plus[[h]] <- q2_boot_sd_plus
    Q2_all_sum_star[[h]] <- q2_all_boot#result_b$Q2_h_star
    Q2_all_sum_star_sd_moins[[h]] <- q2_all_boot_sd_moins
    Q2_all_sum_star_sd_plus[[h]] <- q2_all_boot_sd_plus
    lambdas_out[[h]] <- lambdas_h
    vars_h_boot[[h]] <- vars_boot
    vars_h_boot_sd_moins[[h]] <- vars_boot_sd_moins
    vars_h_boot_sd_plus[[h]] <- vars_boot_sd_plus
    vars_h_boot_single[[h]] <- vars_boot_h
    vars_h_boot_single_sd_moins[[h]] <- vars_boot_h_sd_moins
    vars_h_boot_single_sd_plus[[h]] <- vars_boot_h_sd_plus
    # cov_mean_h[[h]] <- cov_mean
    if(length(id_ALL_TEST)>0){#id_vars_loo_ok)>0){
      q2_max_h[h] <- max(na.omit(q2_boot[id_ALL_TEST]))#id_vars_loo_ok])))
      id_s_cool <- which(q2_boot==q2_max_h[h])# & q2_boot_sd_moins>0 &#q2_boot>tau &#q2_boot_sd_plus>=q2_max_h[h] & q2_boot>tau &#>=q2_max_h[h] & q2_boot>tau &#vars_boot_sd_moins>0)#vars_boot>NZV)# & mod_exists==1)

      if(length(id_s_cool)>0){
        best_id_h <- max(1,min(id_s_cool))
        if(length(best_id_h)>0){
          test_h <- Q2_h_sum_star_sd_moins[[h]][best_id_h]>0
          if(h>1){
            best_id_h_before <- which(lambdas_out[[h-1]]==lambda[h-1])
            test2 <- Q2_all_sum_star_sd_plus[[h-1]][best_id_h_before] <
              Q2_all_sum_star_sd_moins[[h]][best_id_h]
            test_h <- test_h & test2
          }
          best_lam_h <- lambdas_h[best_id_h]
          # Maybe best lambda right of previous best lambda
          if(h<0){#>1){
            if(best_lam_h>lambda[h-1]){
              lambda_min_h <- best_lam_h
              h <- h-1
              u[,h] <- NA
              V_phi[,h] <- NA
              id_ALL_TEST_h[[h]] <- intersect(id_ALL_TEST_h[[h]],which(lambdas_h>=lambda_min_h))
              q2_max_h[h] <- max(na.omit(Q2_h_sum_star[[h]][id_ALL_TEST_h[[h]]]))
              id_s_cool <- which(Q2_h_sum_star[[h]]==q2_max_h[h])
              best_id_h <- max(1,min(id_s_cool))
              best_lam_h <- lambdas_h[best_id_h]
              if(h>1){
                x0 <- x0_deflated[[h]]
                y0 <- y0_deflated[[h]]
                B_previous <- matrix(0,p,q)
                for(s in 1:(h-1)){
                  B_previous <- B_previous + B_r_out[[s]]
                }
              }else{
                x0 <- X_init
                y0 <- Y_init
                B_previous <- matrix(0,p,q)
              }
            }else{
              lambda_min_h <- 0
            }
          }
          # coco <- crossprod(x0,y0)*0
          # for(ii in 1:n_B){
          #   i_b <- sample(1:n,size = n,replace = T)
          #   coco <- coco + crossprod(x0[i_b,,drop=F],y0[i_b,,drop=F])/((n-1)*n_B)
          # }
          # coco <- sort(abs(coco))
          # if(best_lam_h>=max(coco)){
          #   best_lam_h_hypoth <- (best_lam_h + max(coco))/2
          # }else if(best_lam_h<=min(coco)){
          #   best_lam_h_hypoth <- (best_lam_h + min(coco))/2
          # }else{
          #   po_left <- which(coco>best_lam_h)
          #   po_right <- which(coco<best_lam_h)
          #   popo_left <- max(min(po_left)-1,0)
          #   popo_right <- min(max(po_right)+1,length(coco))
          #   best_lam_h_hypoth <- (coco[popo_left]+coco[popo_right])/2
          # }
          # best_id_h <- which.min(abs(best_lam_h_hypoth-lambdas_h))
          # best_lam_h <- lambdas_h[best_id_h]
          # if(!test_h){
          #   test <- F
          # }else{
          ### Find best solution
          # temp <- array(unlist(B_boot), c(p, q, n_B))
          # B_boot_h_0 <- apply(temp, 1:2, median)
          # svdB <- svd(B_boot_h_0,nu=1,nv=1)
          # u_sol_boot_h <- svdB$u
          # t_boot_h <- x0%*%u_sol_boot_h
          #V_sol_boot_h <- svdB$v*svdB$d[1]
          # V_phi[,h] = V_optim_boot_h <- svdB$v
          # norm_t_2 <- sum(t_boot_h^2)
          # if(norm_t_2>1e-9){
          #   V_sol_boot_h <- V_optim_boot_h%*%crossprod(V_optim_boot_h,crossprod(y0,t_boot_h))/norm_t_2
          # }else{
          #   V_sol_boot_h <- matrix(0,q,1)
          # }
          # B_boot_h <- tcrossprod(u_sol_boot_h,V_sol_boot_h)
          if(T){
            m_gogo <-  model_PLS(x = x0,y=y0,lam=lambdas_h[best_id_h],to.scale = F)
            u_sol_boot_h <- m_gogo$U_out
            V_optim_boot_h <- m_gogo$V_optim
            V_phi[,h] <- V_optim_boot_h
            t_boot_h <- x0%*%u_sol_boot_h
            norm_t_0 <- sum(t_boot_h^2)
            if(norm_t_0>1e-9){
              bt <- crossprod(t_boot_h,x0)/norm_t_0
              V_sol_boot_h <- V_optim_boot_h%*%crossprod(V_optim_boot_h,crossprod(y0,t_boot_h))
              V_sol_boot_h <- V_sol_boot_h/norm_t_0
            }else{
              V_sol_boot_h <- V_sol_boot_h*0
            }
            B_boot_h <- tcrossprod(u_sol_boot_h,V_sol_boot_h)
          }
          if(F){
            ## The regression matrix is the result of projection of x0 on bootstrap u
            u_boot_best <- lapply(u_boot[[h]],function(u_l){u_l[,best_id_h]})
            uuuu_0 <- do.call(rbind,u_boot_best);sisi <- sign(uuuu_0%*%t(t(uuuu_0[1,])))
            uuuu <- uuuu_0*matrix(rep(sisi,ncol(uuuu_0)),ncol = ncol(uuuu_0),byrow = F)
            u_sol <- apply(uuuu,2,quantile,probs=0.5)#(apply(uuuu,2,quantile,probs=alpha)+apply(uuuu,2,quantile,probs=1-alpha))/2#
            norm_u <- sqrt(sum(u_sol^2))
            if(norm_u>1e-9){
              u_sol_boot_h <- u_sol/norm_u
            }else{
              u_sol_boot_h <- u_sol
              cat(paste("\n  --  Component ",h," not built, gathered model for u null.\n",sep=""))
            }
            t_boot_h <- x0%*%u_sol_boot_h
            V_boot_best <- lapply(V_boot[[h]],function(v_l){v_l[,best_id_h]})
            # V_optim_boot_h <- apply(do.call(rbind,V_boot_best),2,median)
            # vvvv <- do.call(rbind,V_boot_best)
            vvvv_0 <- do.call(rbind,V_boot_best)
            # sisi <- sign(vvvv_0%*%t(t(vvvv_0[1,])))
            vvvv <- vvvv_0*matrix(rep(sisi,ncol(vvvv_0)),ncol = ncol(vvvv_0),byrow = F)
            V_optim_boot_h <- apply(vvvv,2,quantile,probs=0.5)#(apply(vvvv,2,quantile,probs=alpha)+apply(vvvv,2,quantile,probs=1-alpha))/2
            norm_v <- sqrt(sum(V_optim_boot_h^2))
            if(norm_v>1e-9){
              V_optim_boot_h <- V_optim_boot_h/norm_v
            }else{
              cat(paste("\n  --  Component ",h," not built, gathered model for V null.\n",sep=""))
            }
            V_phi[,h] <- V_optim_boot_h
            norm_t_2 <- sum(t_boot_h^2)
            if(norm_t_2>1e-9){
              V_sol_boot_h <- V_optim_boot_h%*%crossprod(V_optim_boot_h,crossprod(y0,t_boot_h))/norm_t_2
            }else{
              V_sol_boot_h <- matrix(0,q,1)
            }
            B_boot_h <- tcrossprod(u_sol_boot_h,V_sol_boot_h)
          }
          ## The regression matrix is the result of projection of x0 on bootstrap u and v
          # V_model_boot_best <- lapply(V_model_boot,function(v_l){v_l[,best_id_h]})
          # V_model_boot_h <- t(t(apply(do.call(rbind,V_model_boot_best),2,median)))
          # B_boot_second_h <- tcrossprod(u_sol_boot_h,V_model_boot_h)
          # ## The regression matrix is the result of bootstrapping --> the most general
          # B_boot_h_0 <- lapply(B_boot,function(v_l){v_l[,best_id_h]})
          # B_boot_h <- matrix(apply(do.call(rbind,B_boot_h_0),2,median),ncol=q,byrow = F)

          ### End
          sd_y <- apply(Y,2,sd)
          ## Variance per component
          y_est_boot <- x0%*%B_boot_h
          # y_est_boot_ori <- t(apply(y_est_boot,1,function(ll){ll*sd_y+mu_y}))
          varia_expl <- 1 - sum((Y_init-y_est_boot)^2)/sum(RSS0)
          ##
          if(varia_expl>1e-9){#} varia_expl>NZV){#cum_vari_h>vars_boot_sd_moins[best_id_h] &
            VAR_h_s[h] <- varia_expl#c(VAR_h_s,varia_expl)
            Q2_tot_s[h] <- Q2_all_sum_star[[h]][best_id_h]#c(Q2_tot_s,Q2_all_sum_star[[h]][best_id_h])
            Q2_all_sum_star_boxplot[[h]] <- res_measure[[h]][which(res_measure[[h]][,1]==best_id_h),3]
            R2_all_sum_star_boxplot[[h]] <- res_measure[[h]][which(res_measure[[h]][,1]==best_id_h),4]
            lambda[h] <- best_lam_h
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

            V_r <- V_sol_boot_h#m_plus$V_out
            # y0_plus <- y0 - m_plus$e_y
            # id_y_sel <- which(abs(V_r)>1e-9)
            if(verbose){
              cat(paste("\nComponent ",h,
                        "   lambda=",round(best_lam_h,3),
                        "   var.expl._h=",round(varia_expl*100),"%",
                        "   Q2_h=",round(q2_max_h[h]*100)/100,
                        sep=""))
            }
            B_r_out[[h]] <- B_boot_h
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
            B_boot_h <- tcrossprod(u[,h],V_sol_boot_h)
            B_previous <- B_previous + tcrossprod(u[,h],V_sol_boot_h)
            ## Variance total
            y_est_boot_tot <- X_init%*%B_previous
            # y_est_boot_tot_origin <- t(apply(y_est_boot_tot,1,function(ll){ll*sd_y+mu_y}))
            cum_vari_h <- 1-sum((Y_init-y_est_boot_tot)^2)/sum(RSS0)
            CUM_VAR_h_s[h] <- cum_vari_h#c(CUM_VAR_h_s,cum_vari_h)
            lambdas_h <- seq(0,lambda_max,length.out = N_lambdas)#seq(min(lambdas_h),best_lam_h,length.out = length(lambdas_h))
            lambdas_h <- unique(lambdas_h)
            N_lambdas <- length(lambdas_h)
          }else{
            test <- F
          }
        }else{
          test <- F
        }
      }else{
        test <- F
        cat(paste("\n  --  Component ",h," not built, no accessible parameter.\n",sep=""))
      }
    }else{
      test <- F
      cat(paste("\n  --  Component ",h," not built, no accessible parameter.\n",sep=""))
    }
  }
  h_opt <- h - 1
  x0_center <- scale(do.call(cbind,Xs_init),scale = F,center=center)
  if(h_opt>0){
    Q2_cum_star <- 1-prod(1-Q2_sum_star[1:h_opt])
    for(h in 1:h_opt){
      B_r_out[[h]] <- B_r_out[[h]]*sd_y_x_inv
      B <- B + B_r_out[[h]]
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

    if(verbose){
      cat(paste("\nComplete model   ",
                "   Q2=",round(Q2_tot_s[h_opt],2)," \n\n",sep=""))
    }
    # Prepare outputs
    optimal_parameters <- list(lambda=lambda_sol,R=h_opt,
                               Q2=Q2_tot_s[h_opt])
    parameters <- list(RSS0=RSS0)#,RSS0_star=RSS0_star,RSS_y=RSS_y[1:(h_opt+1),],PRESS_y=PRESS_y[1:(h_opt+1),],Q2_y=Q2_y[1:(h_opt+1),])
    bootstrap <- list(lambdas_h=lambdas_out,
                      Q2_tot_s=na.omit(Q2_tot_s),
                      q2_max_h=q2_max_h[1:h_opt],
                      Q2_h_star=Q2_h_sum_star,
                      Q2_h_star_sd_plus=Q2_h_sum_star_sd_plus,
                      Q2_h_star_sd_moins=Q2_h_sum_star_sd_moins,
                      Q2_all_sum_star=Q2_all_sum_star,
                      Q2_all_sum_star_sd_plus=Q2_all_sum_star_sd_plus,
                      Q2_all_sum_star_sd_moins=Q2_all_sum_star_sd_moins,
                      Q2_all_sum_star_boxplot=Q2_all_sum_star_boxplot,
                      R2_all_sum_star_boxplot=R2_all_sum_star_boxplot,
                      vars_h_boot=vars_h_boot,
                      vars_h_boot_sd_plus=vars_h_boot_sd_plus,
                      vars_h_boot_sd_moins=vars_h_boot_sd_moins,
                      vars_h_boot_single=vars_h_boot_single,
                      vars_h_boot_single_sd_plus=vars_h_boot_single_sd_plus,
                      vars_h_boot_single_sd_moins=vars_h_boot_single_sd_moins,
                      cov_mean_h=cov_mean_h,prop_models_ok=prop_models_ok)
    res <- list(optimal_parameters=optimal_parameters,
                bootstrap=bootstrap,
                Us=Us,V=V[,1:h_opt,drop=F],B_cbind=B,
                id_ALL_TEST_h=id_ALL_TEST_h,
                x0_deflated=x0_deflated,y0_deflated=y0_deflated,
                alpha=alpha,
                t_h=t_h[,1:h_opt,drop=F],
                explained_variance=na.omit(VAR_h_s)*100,
                explained_variance_cum=na.omit(CUM_VAR_h_s)*100,
                NZV=NZV,#score_X = m_plus$score_x,
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
  cols <- c(RColorBrewer::brewer.pal(max(h_opt,3),"Set1")[1:h_opt],"gray")
  layout(matrix(c(1,1,2,2,3,4,4,5,5,6), 2, 5, byrow = TRUE))
  par(mar=c(3,3,2,1),mgp=c(2,1,0))
  aa <- min(res$alpha,1-res$alpha)
  ## R2 stuff
  if(T){
    ## Plot R2_h
    for(h in 1:(h_opt+1)){
      id_h <- res$id_ALL_TEST_h[[h]]
      add <- T;if(h==1) add<-F
      matplot(res$bootstrap$lambdas_h[[h]],
              cbind(res$bootstrap$vars_h_boot_single[[h]],
                    res$bootstrap$vars_h_boot_single_sd_plus[[h]],
                    res$bootstrap$vars_h_boot_single_sd_moins[[h]])
              ,lty=c(1,2,2),type="l",col=cols[h],add=add,ylim=c(-0.1,1.15),
              main=expression("R"["B,h"]^"2"),ylab="",xlab=expression(lambda))
      if(length(id_h)>0){
        matplot(res$bootstrap$lambdas_h[[h]][id_h],
                cbind(res$bootstrap$vars_h_boot_single[[h]][id_h],
                      res$bootstrap$vars_h_boot_single_sd_plus[[h]][id_h],
                      res$bootstrap$vars_h_boot_single_sd_moins[[h]][id_h])
                ,lty=c(1,2,2),type="l",col=cols[h],add=T,lwd=2)
      }
      if(h<=h_opt){
        points(res$optimal_parameters$lambda[h],
               res$bootstrap$vars_h_boot_single[[h]][which(
                 res$bootstrap$lambdas_h[[h]]==res$optimal_parameters$lambda[h])],
               col=cols[h],pch=16)
      }
      legend(bty="n","topleft",fill=cols,legend = c(unlist(lapply(1:h_opt,function(h){
        paste("Component ",h," (",round(res$explained_variance[h]),"%)",sep="",collapse = "")})),
        "Not selected component"))
      legend("top",ncol=2,c("Median",paste(round(c(aa,1-aa),2),collapse=" / ")),lty=c(1,2),bty="n")
    }
    abline(v=res$optimal_parameters$lambda,col=cols,lty=3)
    ## Plot R2
    for(h in 1:(h_opt)){
      id_h <- res$id_ALL_TEST_h[[h]]
      add <- T;if(h==1) add<-F
      matplot(res$bootstrap$lambdas_h[[h]],
              cbind(res$bootstrap$vars_h_boot[[h]],
                    res$bootstrap$vars_h_boot_sd_plus[[h]],
                    res$bootstrap$vars_h_boot_sd_moins[[h]])
              ,lty=c(1,2,2),type="l",col=cols[h],add=add,ylim=c(-0.1,1.15),
              main=expression("R"["B"]^"2"),ylab="",xlab=expression(lambda))
      if(length(id_h)>0){
        matplot(res$bootstrap$lambdas_h[[h]][id_h],
                cbind(res$bootstrap$vars_h_boot[[h]][id_h],
                      res$bootstrap$vars_h_boot_sd_plus[[h]][id_h],
                      res$bootstrap$vars_h_boot_sd_moins[[h]][id_h])
                ,lty=c(1,2,2),type="l",col=cols[h],add=T,lwd=2)
      }
      if(h<=h_opt){
        points(res$optimal_parameters$lambda[h],
               res$bootstrap$vars_h_boot[[h]][which(
                 res$bootstrap$lambdas_h[[h]]==res$optimal_parameters$lambda[h])],
               col=cols[h],pch=16)
        points(res$optimal_parameters$lambda[h],
               res$explained_variance_cum[h],
               col=cols[h],pch=4)
      }
      legend(bty="n","topleft",fill=cols,legend = c(unlist(lapply(1:h_opt,function(h){
        paste("Component ",h," (",round(res$explained_variance[h]),"%)",sep="",collapse = "")})),
        "Not selected component"))
      legend("top",ncol=2,c("Median",paste(round(c(aa,1-aa),2),collapse=" / ")),lty=c(1,2),bty="n")
    }
    abline(v=res$optimal_parameters$lambda,col=cols,lty=3)
    ## Boxplot R2
    ddff <- data.frame(do.call(rbind,lapply(1:h_opt,function(ii,bb){cbind(ii,bb[[ii]])},
                                            res$bootstrap$R2_all_sum_star_boxplot)))
    names(ddff) <- c("h","R2")
    bobo <- boxplot(R2~h,ddff,ylim=c(-0.1,1.15),border=cols,main=expression("R"[B]^"2"),xlab="Component",ylab="")
    # bobo <- boxplot(R2~h,ddff,plot=F)
    # bobo$out <- NULL
    # bobo$group <- NULL
    # for(h in 1:h_opt){
    #   qt <- quantile(res$bootstrap$R2_all_sum_star_boxplot[[h]],probs = c(1/4,res$alpha,0.5,1-res$alpha,1-1/4))
    #   bobo$stats[,h] <- matrix(qt,ncol=1)
    #   out <- res$bootstrap$R2_all_sum_star_boxplot[[h]][-which(
    #     res$bootstrap$R2_all_sum_star_boxplot[[h]]<qt[2] |
    #       res$bootstrap$R2_all_sum_star_boxplot[[h]]>qt[4]
    #   )]
    #   bobo$out <- c(bobo$out,out)
    #   bobo$group <- c(bobo$group,rep(h,length(out)))
    # }
    # bxp(bobo,ylim=c(-0.1,1.15),border=cols,main=expression("R"[B]^"2"),xlab="Component",ylab="")
  }
  ## Q2 stuff
  if(T){
    ## Plot Q2_h
    for(h in 1:(h_opt+1)){
      id_h <- res$id_ALL_TEST_h[[h]]
      add <- T;if(h==1) add<-F
      matplot(res$bootstrap$lambdas_h[[h]],
              cbind(res$bootstrap$Q2_h_star[[h]],
                    res$bootstrap$Q2_h_star_sd_plus[[h]],
                    res$bootstrap$Q2_h_star_sd_moins[[h]])
              ,lty=c(1,2,2),type="l",col=cols[h],add=add,ylim=c(-0.1,1.15),
              main=expression("Q"["B,h"]^"2"),ylab="",xlab=expression(lambda))
      if(length(id_h)>0){
        matplot(res$bootstrap$lambdas_h[[h]][id_h],
                cbind(res$bootstrap$Q2_h_star[[h]][id_h],
                      res$bootstrap$Q2_h_star_sd_plus[[h]][id_h],
                      res$bootstrap$Q2_h_star_sd_moins[[h]][id_h])
                ,lty=c(1,2,2),type="l",col=cols[h],add=T,lwd=2)
      }
      if(h<=h_opt){
        points(res$optimal_parameters$lambda[h],
               res$bootstrap$Q2_h_star[[h]][which(
                 res$bootstrap$lambdas_h[[h]]==res$optimal_parameters$lambda[h])],
               col=cols[h],pch=16)
      }
      legend(bty="n","topleft",fill=cols,legend = c(unlist(lapply(1:h_opt,function(h){
        paste("Component ",h," (",round(res$explained_variance[h]),"%)",sep="",collapse = "")})),
        "Not selected component"))
      legend("top",ncol=2,c("Median",paste(round(c(aa,1-aa),2),collapse=" / ")),lty=c(1,2),bty="n")
    }
    abline(v=res$optimal_parameters$lambda,col=cols,lty=3)
    ## Plot Q2
    for(h in 1:(h_opt)){
      id_h <- res$id_ALL_TEST_h[[h]]
      add <- T;if(h==1) add<-F
      matplot(res$bootstrap$lambdas_h[[h]],
              cbind(res$bootstrap$Q2_all_sum_star[[h]],
                    res$bootstrap$Q2_all_sum_star_sd_plus[[h]],
                    res$bootstrap$Q2_all_sum_star_sd_moins[[h]])
              ,lty=c(1,2,2),type="l",col=cols[h],add=add,ylim=c(-0.1,1.15),
              main=expression("Q"[B]^"2"),ylab="",xlab=expression(lambda))
      if(length(id_h)>0){
        matplot(res$bootstrap$lambdas_h[[h]][id_h],
                cbind(res$bootstrap$Q2_all_sum_star[[h]][id_h],
                      res$bootstrap$Q2_all_sum_star_sd_plus[[h]][id_h],
                      res$bootstrap$Q2_all_sum_star_sd_moins[[h]][id_h])
                ,lty=c(1,2,2),type="l",col=cols[h],add=T,lwd=2)
      }
      if(h<=h_opt){
        points(res$optimal_parameters$lambda[h],
               res$bootstrap$Q2_all_sum_star[[h]][which(
                 res$bootstrap$lambdas_h[[h]]==res$optimal_parameters$lambda[h])],
               col=cols[h],pch=16)
      }
      legend(bty="n","topleft",fill=cols,legend = c(unlist(lapply(1:h_opt,function(h){
        paste("Component ",h," (",round(res$explained_variance[h]),"%)",sep="",collapse = "")})),
        "Not selected component"))
      legend("top",ncol=2,c("Median",paste(round(c(aa,1-aa),2),collapse=" / ")),lty=c(1,2),bty="n")
    }
    abline(v=res$optimal_parameters$lambda,col=cols,lty=3)
    ## Boxplot Q2
    ddff <- data.frame(do.call(rbind,lapply(1:h_opt,function(ii,bb){cbind(ii,bb[[ii]])},
                                            res$bootstrap$Q2_all_sum_star_boxplot)))
    names(ddff) <- c("h","Q2")
    bobo <- boxplot(Q2~h,ddff,ylim=c(-0.1,1.15),border=cols,main=expression("Q"[B]^"2"),xlab="Component",ylab="")
    # for(h in 1:h_opt){
    #   bobo$stats[,h] <- matrix(quantile(res$bootstrap$Q2_all_sum_star_boxplot[[h]],
    #                                     probs = c(1/4,res$alpha,0.5,1-res$alpha,1-1/4)),ncol=1)
    # }
    # bxp(bobo)
  }

}
