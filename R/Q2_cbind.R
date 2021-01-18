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
#'
#' @return
#' @export
#'
#' @useDynLib ddsPLS
do_one_component <- function(x0,y0,n,p,q,COV,abs_COV,max_COV,lam,remove_COV=NULL){
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
      svd_mix_Y <- svd(COV_COV_high,nu = 1,nv = 1)
      U0[id_x_high,] <- svd_mix_Y$v
      V0[id_y_high,] <- svd_mix_Y$u
    }
  }else{
    U0 <- matrix(0,p,1)
    V0 <- matrix(0,q,1)
  }
  t <- x0%*%U0
  y0_mask <- y0;y0_mask[which(abs(V0)<1e-9),] <- 0;y0_mask[which(abs(V0)>=1e-9),] <- 1
  V_svd0 <- crossprod(y0_mask,t)#
  # V_svd0 <- V0%*%crossprod(V0,crossprod(y0,t)) #crossprod(y0,t) #
  norm_t_0 <- sum(t^2)
  if(norm_t_0>1e-9){
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
#' @param to.scale to.scale
#'
#' @return
#' @export
#'
#' @useDynLib ddsPLS
model_PLS <- function(x,y,lam,deflatX=T,R=1,remove_COV=NULL,RSS0=NULL,
                      to.scale=T,COV_init=NULL){
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
    if(length(which(!(is.na(c(abs_COV)))))==0) browser()
    max_COV <- max(na.omit(c(abs_COV)))
    lam_r <- lam
    if(length(lam)>1) lam_r <- lam[r]
    if(lam_r<max_COV){
      c_h <- do_one_component(x0 = x0,y0 = y0,n = n,p = p,q = q,COV = COV,abs_COV = abs_COV,
                              max_COV=max_COV,lam = lam_r,remove_COV=remove_COV)
      t_r <- c_h$t ; U0  <- c_h$U0 ; V_svd  <- c_h$V_svd ; V0  <- c_h$V0
      ## DEFLAT ##
      if(sum(U0^2)>1e-9){
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

bootstrap_pls <- function(X_init,Y_init,h=1,lambdas=0,B_previous,
                          lambda_prev=NA,u,v){
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
  if(length(na.omit(lambda_prev))==0){
    u <- matrix(NA,p,n)
    v <- matrix(NA,q,n)
  }
  if(h>1){
    U_reconstruct[,1:(h-1)] <- u[,1:(h-1)]
    for(r in 1:(h-1)){
      if(length(na.omit(lambda_prev))==0){
        m_gogo <-  model_PLS(x = X_r,y=Y_r,lam=lambda_prev[r],
                             to.scale = F)
        u[,r] <- m_gogo$U_out
        v[,r] <- m_gogo$V_optim
      }
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

sPLS_lasso <- function(x,y,R=1,kx,ky,NZV=1e-9,to.scale=T,deflatX=T){
  ########################################################################
  ########################################################################
  do_lasso_pls <- function(a,b,w,keep,NZV=1e-9){
    k <- ncol(b)
    Score <- a%*%w
    theta_0 <- crossprod(b,Score)
    if(k>keep){
      l <- sort(abs(theta_0))[k-keep]
    }else{
      l <- 0
    }
    theta <- abs(theta_0)-l
    theta[which(theta<0)] <- 0
    theta <- theta*sign(theta_0)
    no_theta <-sqrt(sum(theta^2))
    if(no_theta>NZV){
      theta <- theta/no_theta
    }else{
      theta <- theta*0
    }
    theta
  }
  build_one_comp_lasso <- function(xr,yr,kx,ky,NZV=1e-9,maxIter=20){
    M <- crossprod(xr,yr);ssvvdd <- svd(M,nu=1,nv=1)
    v <- ssvvdd$v;u <- ssvvdd$u;err <- 1
    iter <- 1
    while(err>NZV & iter<maxIter){
      v_init <- v
      u_init <- u
      u <- do_lasso_pls(yr,xr,v,kx,NZV)
      v <- do_lasso_pls(xr,yr,u,ky,NZV)
      err <- max(sqrt(sum((v-v_init)^2)),sqrt(sum((u-u_init)^2)))
      iter <- iter + 1
    }
    list(u=u,v=v)
  }
  ########################################################################
  ########################################################################
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
  var_y_init <- sum(y^2)
  U_out = U_star <- matrix(0,p,R)
  score_x <- matrix(0,n,R)
  C <- matrix(0,p,R)
  D <- matrix(0,q,R)
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
  for(r in 1:R){
    c_h <- build_one_comp_lasso(x0,y0,kx=kx[r],ky=ky[r],NZV=NZV)
    U0  <- c_h$u ; V_svd  <- c_h$v# ; V0  <- c_h$V0
    t_r <- x0%*%U0
    ## DEFLAT ##
    if(sum(U0^2)>1e-9){
      score_x[,r] <- t_r
      bt <- crossprod(t_r,x0)/sum(t_r^2)
      C[,r] <- t(bt)
      D[,r] <- crossprod(y0,t_r)/sum(t_r^2)
      U_out[,r] <- U0
      U_star[,r] <- U0
      V_out[,r] <- V_svd
      B_r_0 <- tcrossprod(U0,V_svd)
      # y_est <- y_est + x0%*%B_r_0
      bXr[r,] <- bt
      bYr[r,] <- t(V_svd)
      # if(r!=1){
      #   for(s_r in (r-1):1){
      #     U_star[,r] <- U_star[,r]-U_star[,s_r]*sum(bXr[s_r,]*U_star[,r])
      #   }
      # }
      # if(to.scale){
      #   B_r[[r]] <- tcrossprod(U_star[,r],V_svd)*sd_y_x_inv
      # }else{
      #   B_r[[r]] <- tcrossprod(U_star[,r],V_svd)
      # }
      # B <- B + B_r[[r]]
      U_star_cl <- U_out[,1:r]%*%solve(crossprod(C[,1:r],U_out[,1:r]))
      B <- tcrossprod(U_star_cl,D)
      y_plus_un <- tcrossprod(t_r,D[,r])#tcrossprod(t_r,V_out[,r,drop=F])
      x0_plus <- tcrossprod(t_r,C[,r])#t_r%*%bt
      var_expl[r] <- 1-sum((y_init-mu_y-y_plus_un)^2)/sum((y_init-mu_y)^2)
      y0 <- y0 - y_plus_un
      if(deflatX){
        x0 <- x0 - x0_plus
      }
    }
  }
  if(to.scale){
    y_est <- x_init%*%B - mu_x%*%B
  }
  classic <- list(C=C,D=D,U=U_out,
                  U_star=U_star_cl,B=B)
  list(no_model=no_model,U_out=U_out,U_star=U_star,
       V_out=V_out,B=B,B_r=B_r,var_expl=var_expl,classic=classic,
       score_x=score_x,y_est=y_est,bXr=bXr,bYr=bYr,e_x=x0,e_y=y0)

}

bootstrap_sPLS_lasso <- function(X_init,Y_init,h=1,paras_keep=paras_keep,
                                 paras_prev = paras_prev,u,v,sparse=T,NZV=1e-9){
  N_paras <- nrow(paras_keep)
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
  no_prev <- nrow(na.omit(paras_prev))==0
  C <- matrix(0,p,h)
  D <- matrix(0,q,h)
  if(T){
    u <- matrix(NA,p,n)
    v <- matrix(NA,q,n)
  }
  if(h>1){
    for(r in 1:(h-1)){
      if(T){
        kX_r <- p;kY_r <- q
        if(sparse){
          kX_r <- paras_prev[r,1];kY_r <- paras_prev[r,2]
        }
        m_gogo <- sPLS_lasso(x=X_r,y=Y_r,kx=kX_r,ky=kY_r,
                             to.scale=F,R=1)
        u[,r] <- m_gogo$U_out
        v[,r] <- m_gogo$V_out
      }
      t_r <- X_r%*%u[,r]
      t_all[,r] <- t_r
      bt <- crossprod(t_r,X_r)/sum(t_r^2)
      C[,r] <- t(bt)
      D[,r] <- crossprod(Y_r,t_r)/sum(t_r^2)
      U_star_cl <- u[,1:r,drop=F]%*%solve(crossprod(C[,1:r,drop=F],u[,1:r,drop=F]))
      B_youyou <- tcrossprod(U_star_cl,D[,1:r,drop=F])
      y_r <- tcrossprod(t_r,D[,r,drop=F])#tcrossprod(t_r,V_out[,r,drop=F])
      V_svd <- v[,r]
      # Do deflation
      Y_r <- Y_r - y_r
      X_r <- X_r - tcrossprod(t_r,C[,r,drop=F])
      X_defla[[r+1]] <- X_r
      Y_defla[[r+1]] <- Y_r
    }
    # Create weights
    # Build regression matrix
  }
  vars_expl_star = vars_expl = vars_expl_h = Q2_star = Q2 = Q2_all = model_exists <- rep(0,N_paras)
  B_out <- matrix(0,q*p,N_paras)
  u_out <- matrix(0,p,N_paras)
  V_optim_phi = V_model <- matrix(0,q,N_paras)
  B_next_out <- list()
  for(i_l in 1:N_paras){
    kX_r <- p;kY_r <- q
    if(sparse){
      kX_r <- paras_keep[i_l,1];kY_r <- paras_keep[i_l,2]
    }
    m_1 <- sPLS_lasso(x=X_r,y=Y_r,kx=kX_r,ky=kY_r,
                      to.scale=F,R=1)
    u_il <- m_1$U_out
    V_il <- m_1$V_out

    u_out[,i_l] <- u_il
    V_optim_phi[,i_l] <- V_il
    t_r <- X_r%*%u_il
    bt <- crossprod(t_r,X_r)/sum(t_r^2)
    C[,h] <- t(bt)
    D[,h] <- crossprod(Y_r,t_r)/sum(t_r^2)
    u[,h] <- u_il
    U_star_cl <- u[,1:h,drop=F]%*%solve(crossprod(C[,1:h,drop=F],u[,1:h,drop=F]))
    B_all <- tcrossprod(U_star_cl,D)
    y_r <- tcrossprod(t_r,D[,h,drop=F])#tcrossprod(t_r,V_out[,r,drop=F])
    V_svd <- V_il
    # Create weights

    # Modify regression matrix
    # B_next <- tcrossprod(U_reconstruct[,h],V_svd)
    # All components
    # B_all <- B_youyou + B_next
    y_train_pred <- X_train%*%B_all
    y_test_pred <- X_test_normalize%*%B_all
    # Previous components
    y_train_pred_next <- X_train%*%(B_all-B_youyou)
    y_test_pred_RSS <- X_test_normalize%*%B_youyou

    # y_train_pred_next <- X_train%*%B_next
    # y_test_pred_RSS <- X_test_normalize%*%B_youyou
    # Compute criterions
    vars_expl[i_l] <- 1 - sum( (Y_train-y_train_pred)^2 ) / sum( (Y_train)^2 )
    vars_expl_h[i_l] <- 1 - sum( (Y_train-y_train_pred_next)^2 ) / sum((Y_train)^2)
    Q2[i_l] <- 1 - sum( (Y_test-y_test_pred)^2 ) / sum((Y_test-y_test_pred_RSS)^2)
    MU_test <- matrix(rep(mu_y,length(id_OOB)),nrow = length(id_OOB),byrow = T)
    Q2_all[i_l] <- 1 - sum( (Y_test-y_test_pred)^2 ) / sum((Y_test)^2)#-MU_test)^2)

  }
  list(u_out=u_out,V_optim_phi=V_optim_phi,id_OOB=id_OOB,
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
#' @param keepXs keepXs
#' @param keepXs keepXs
#' @param n_B n_B
#' @param deflatX deflatX
#' @param center center
#' @param NCORES number of cores
#' @param NZV NZV
#' @param verbose verbose
#' @param sparse sparse
#'
#' @return
#' @export
#'
#' @useDynLib ddsPLS
Q2_boot_sPLS <- function(Xs,Y,keepXs = 1,
                         keepYs=1,
                         n_B=20,whichCriterion="Q2",
                         deflatX=T,NCORES=1,center=T,
                         NZV=1e-9,verbose=F,sparse=F){
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
  Q2_sum_star = id_paras_out <- rep(NA,n)
  paras_prev <- matrix(NA,n,2)
  ps <- unlist(lapply(Xs_init,ncol))
  X_init <- do.call(cbind,Xs_init)
  x0 <- X_init
  u <- matrix(NA,sum(ps),n)
  V_phi = V <- matrix(NA,q,n)
  P <- matrix(0,sum(ps),n)
  B <- matrix(0,sum(ps),q)
  C <- matrix(0,p,n)
  D <- matrix(0,q,n)
  t_h <- matrix(NA,n,n)
  Us = Bs = B_r_out = Q2_mean = Q2_all_mean =
    Q2_h_sum_star = Q2_h_sum_star_sd_moins = Q2_h_sum_star_sd_plus =
    Q2_all_sum_star = Q2_all_sum_star_sd_moins = Q2_all_sum_star_sd_plus =
    vars_h_boot = vars_h_boot_sd_moins = vars_h_boot_sd_plus =
    vars_h_boot_single = vars_h_boot_single_sd_moins = vars_h_boot_single_sd_plus <- list()
  q2_max_h = VAR_h_s = CUM_VAR_h_s = Q2_tot_s <- rep(NA,n)
  B_tot_LOO <- matrix(0,n,sum(ps)*q)

  ### For each 'h' component, look for the best keepX/keepY
  test <- T
  y0 <- Y_init
  h <- 0
  id_ALL_TEST_h = Q2_all_sum_star_boxplot = R2_all_sum_star_boxplot <- list()
  K_h <- 1:K
  # Check lambdas and stuff
  paras_keep <- expand.grid(keepXs,keepYs)
  N_paras <- nrow(paras_keep)
  ncomps = Q2_TOTAL <- rep(NA,N_paras)
  x0_deflated = y0_deflated = cov_mean_h = prop_models_ok <- list()
  remove_COV <- matrix(0,q,p)
  B_previous <- matrix(0,p,q)
  V_boot = u_boot = res_measure <- list()
  while(test){
    h <- h + 1
    Bis <- list()
    pos_decoupe <- 1
    Q2_Boot = vars = Q2_h_sum <- rep(0,N_paras)
    RSS_il_y = var_sel = Q2_h_y = RSS_il_y_star = RSS_h_moins_1_star = PRESS_il_y <- matrix(0,N_paras,q)
    Q2_B = vars_in = Q2_h_clas <- matrix(NA,n_B,N_paras)
    NCORES_w <- min(NCORES,n_B)
    `%my_do%` <- ifelse(NCORES_w!=1,{
      out<-`%dopar%`;cl <- makeCluster(NCORES_w)
      registerDoParallel(cl);out},{out <- `%do%`;out})
    Q2_star_bootstrat <- foreach(i_B=1:n_B,.packages = "ddsPLS",.combine='c',.multicombine=TRUE) %my_do% {
      source('~/Documents/GitHub/ddsPLS/R/Q2_cbind.R')
      out <- bootstrap_sPLS_lasso(X_init=X_init,Y_init=Y_init,
                                  u=u,v=V_phi,h=h,paras_keep=paras_keep,
                                  paras_prev = paras_prev,sparse=sparse)
      res_measure_in <- cbind(1:N_paras,
                              out$Q2,
                              out$Q2_all,
                              out$vars_expl,
                              out$vars_expl_h,
                              length(out$id_OOB),
                              i_B)
      list(res=res_measure_in,V_optim_phi=out$V_optim_phi,u=out$u_out)
    }
    if(NCORES_w!=1)stopCluster(cl)
    KK <- length(Q2_star_bootstrat)/n_B ; id_pos_boot <- (1:n_B)*KK
    res_measure[[h]] <- do.call(rbind,Q2_star_bootstrat[id_pos_boot-2])
    V_boot[[h]] <- Q2_star_bootstrat[id_pos_boot-1]
    u_boot[[h]] <- Q2_star_bootstrat[id_pos_boot]
    vars_boot = vars_boot_sd_plus = vars_boot_sd_moins =
      vars_boot_h = vars_boot_h_sd_plus = vars_boot_h_sd_moins =
      q2_boot = q2_boot_sd_plus = q2_boot_sd_moins =
      q2_all_boot = q2_all_boot_sd_plus = q2_all_boot_sd_moins =
      q2_boot_mean = q2_all_boot_mean <- rep(0,N_paras)
    test_batchas <- TRUE

    i_l <- 1
    while(test_batchas){
      pos <- which(res_measure[[h]][,1]==i_l)
      toto <- na.omit(res_measure[[h]][pos,c(2,3,4,5)])
      if(length(toto)>0){
        m <- 1
        q2_boot[i_l] <- mean(toto[,m])#Q2_qtles[2]
        # Q2_qtles <- mean(toto[,m])+c(-sd(toto[,m]),0,sd(toto[,m]))#quantile(toto[,1],probs = probabilities)#quantile(toto[,1],probs = probabilities)#quantile(toto[,m],c(1/10,1/2,9/10))#
        q2_boot_sd_plus[i_l] <- mean(toto[,m])+sd(toto[,m])#Q2_qtles[3]
        q2_boot_sd_moins[i_l] <- mean(toto[,m])-sd(toto[,m])#Q2_qtles[1]
        q2_boot_mean[i_l] <- mean(toto[,m])/(sd(toto[,m]))#(1+sd(toto[,m])/n)
        #Q2_qtles[2]/(Q2_qtles[3]-Q2_qtles[1])#mean(toto[,m])/sd(toto[,m])
        m <- 2
        # Q2_all_qtles <- mean(toto[,m])+c(-sd(toto[,m]),0,sd(toto[,m]))#quantile(toto[,2],probs = probabilities)#quantile(toto[,2],probs = probabilities)
        q2_all_boot[i_l] <- mean(toto[,m])#Q2_all_qtles[2]
        q2_all_boot_sd_plus[i_l] <- mean(toto[,m])+sd(toto[,m])#Q2_all_qtles[3]
        q2_all_boot_sd_moins[i_l] <- mean(toto[,m])-sd(toto[,m])#Q2_all_qtles[1]
        q2_all_boot_mean[i_l] <- mean(toto[,m])/(sd(toto[,m]))#q2_boot[i_l]/(q2_boot_sd_plus[i_l]-q2_boot_sd_moins[i_l])
        m <- 3
        # VARS_qtles <- mean(toto[,m])+c(-sd(toto[,m]),0,sd(toto[,m]))#quantile(toto[,3],probs = probabilities)##quantile(toto[,3],probs = probabilities)
        vars_boot[i_l] <- mean(toto[,m])#VARS_qtles[2]
        vars_boot_sd_plus[i_l] <- mean(toto[,m])+sd(toto[,m])#VARS_qtles[3]
        vars_boot_sd_moins[i_l] <- mean(toto[,m])-sd(toto[,m])#VARS_qtles[1]
        m <- 4
        # VARS_qtles <- mean(toto[,m])+c(-sd(toto[,m]),0,sd(toto[,m]))#quantile(toto[,4],probs = probabilities)##quantile(toto[,4],probs = probabilities)
        vars_boot_h[i_l] <- mean(toto[,m])#VARS_qtles[2]
        vars_boot_h_sd_plus[i_l] <- mean(toto[,m])+sd(toto[,m])#VARS_qtles[3]
        vars_boot_h_sd_moins[i_l] <- mean(toto[,m])-sd(toto[,m])#VARS_qtles[1]
      }
      if(is.na(q2_boot[i_l])|is.nan(q2_boot[i_l])|i_l==N_paras){
        test_batchas <- F
      }
      i_l <- i_l + 1
    }
    vars <- vars_boot
    id_ALL_TEST <- which(q2_boot > 0)
    if(h!=1){
      id_test_h <- which(
        q2_all_boot>Q2_all_sum_star[[h-1]][id_paras_out[h-1]])
      id_ALL_TEST <- intersect(id_ALL_TEST,id_test_h)
    }
    id_ALL_TEST_h[[h]] <- id_ALL_TEST
    Q2_h_sum_star[[h]] <- q2_boot
    Q2_h_sum_star_sd_moins[[h]] <- q2_boot_sd_moins
    Q2_h_sum_star_sd_plus[[h]] <- q2_boot_sd_plus
    Q2_mean[[h]] <- q2_boot_mean
    Q2_all_sum_star[[h]] <- q2_all_boot
    Q2_all_sum_star_sd_moins[[h]] <- q2_all_boot_sd_moins
    Q2_all_sum_star_sd_plus[[h]] <- q2_all_boot_sd_plus
    Q2_all_mean[[h]] <- q2_all_boot_mean
    vars_h_boot[[h]] <- vars_boot
    vars_h_boot_sd_moins[[h]] <- vars_boot_sd_moins
    vars_h_boot_sd_plus[[h]] <- vars_boot_sd_plus
    vars_h_boot_single[[h]] <- vars_boot_h
    vars_h_boot_single_sd_moins[[h]] <- vars_boot_h_sd_moins
    vars_h_boot_single_sd_plus[[h]] <- vars_boot_h_sd_plus
    if(length(id_ALL_TEST)>0){
      if(whichCriterion=="Q2"){
        diff_R2_Q2 <- abs(q2_boot)
        q2_max_h[h] <- max(na.omit(diff_R2_Q2[id_ALL_TEST]))#max(na.omit(q2_boot_mean[id_ALL_TEST]))
        id_s_cool <- which(diff_R2_Q2==q2_max_h[h])
      }else{
        diff_R2_Q2 <- abs(q2_boot-vars_boot)
        q2_max_h[h] <- min(na.omit(diff_R2_Q2[id_ALL_TEST]))#max(na.omit(q2_boot_mean[id_ALL_TEST]))
        id_s_cool <- which(diff_R2_Q2==q2_max_h[h])
      }
      if(length(id_s_cool)>0){
        best_id_h <- id_s_cool
        if(length(best_id_h)>0){
          test_h <- q2_boot[best_id_h]>0
          if(h>1){
            best_id_h_before <- id_paras_out[h-1]
            test2 <-
              Q2_all_sum_star[[h-1]][best_id_h_before]<Q2_all_sum_star[[h]][best_id_h]
            test_h <- test_h & test2
          }
          # best_lam_h <- lambdas_h[best_id_h]
          # if(sparse){
          #   m_gogo <- mixOmics::spls(X = x0,Y=y0,ncomp = 1,scale = F,
          #                            keepX = paras_keep[best_id_h,1],
          #                            keepY = paras_keep[best_id_h,2],tol = 1e-9)
          # }else{
          #   m_gogo <-  mixOmics::pls(X = x0,Y=y0,ncomp = 1,scale = F,tol = 1e-9)
          # }
          # model_PLS(x = x0,y=y0,lam=lambdas_h[best_id_h],to.scale = F)
          # u_sol_boot_h <- m_gogo$loadings$X
          # V_optim_boot_h <- m_gogo$loadings$Y
          kX_r <- p;kY_r <- q
          if(sparse){
            kX_r <- paras_keep[best_id_h,1];kY_r <- paras_keep[best_id_h,2]
          }
          m_gogo <- sPLS_lasso(x=x0,y=y0,kx=kX_r,ky=kY_r,
                               to.scale=F,R=1)
          u_sol_boot_h <- m_gogo$U_out
          V_optim_boot_h <- m_gogo$V_out
          V_phi[,h] <- V_optim_boot_h
          t_r <- x0%*%u_sol_boot_h
          bt <- crossprod(t_r,x0)/sum(t_r^2)
          C[,h] <- t(bt)
          D[,h] <- crossprod(y0,t_r)/sum(t_r^2)
          u[,h] <- u_sol_boot_h
          U_star_cl <- u[,1:h,drop=F]%*%solve(crossprod(C[,1:h,drop=F],u[,1:h,drop=F]))
          B_all <- tcrossprod(U_star_cl,D[,1:h,drop=F])
          x_r <- tcrossprod(t_r,C[,h,drop=F])#tcrossprod(t_r,V_out[,r,drop=F])
          y_r <- tcrossprod(t_r,D[,h,drop=F])#tcrossprod(t_r,V_out[,r,drop=F])

          V_sol_boot_h <- V_optim_boot_h
          B_boot_h <- B_all-B_previous
          ### End
          ## Variance per component
          y_est_boot <- X_init%*%B_boot_h
          varia_expl <- 1 - sum((Y_init-y_est_boot)^2)/sum(RSS0)
          ##
          if(varia_expl>1e-9){
            VAR_h_s[h] <- varia_expl
            Q2_tot_s[h] <- Q2_all_sum_star[[h]][best_id_h[1]]
            Q2_all_sum_star_boxplot[[h]] <- res_measure[[h]][which(res_measure[[h]][,1]==best_id_h[1]),3]
            R2_all_sum_star_boxplot[[h]] <- res_measure[[h]][which(res_measure[[h]][,1]==best_id_h[1]),4]
            id_paras_out[h] <- best_id_h[1]
            paras_prev[h,1] <- paras_keep[best_id_h[1],1]
            paras_prev[h,2] <- paras_keep[best_id_h[1],2]
            # Get the regression matrix of the optimal model
            t_h[,h] <- t_r

            V_r <- V_sol_boot_h
            if(verbose){
              cat(paste("\nComponent ",h,
                        "   keepX=",paras_keep[best_id_h[1],1],
                        "   keepY=",paras_keep[best_id_h[1],2],
                        "   var.expl._h=",round(varia_expl*100),"%",
                        "   Q2_h=",round(q2_max_h[h]*100)/100,
                        sep=""))
            }
            B_r_out[[h]] <- B_all
            V[,h] <- V_r
            y0 <- Y_init-y_r
            y0_deflated[[h]] <- y0
            if(deflatX){
              x0 <- X_init - x_r
              x0_deflated[[h]] <- x0
            }
            B_previous <- B_all
            ## Variance total
            y_est_boot_tot <- X_init%*%B_previous
            cum_vari_h <- 1-sum((Y_init-y_est_boot_tot)^2)/sum(RSS0)
            CUM_VAR_h_s[h] <- cum_vari_h
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
    }
    B <- B_r_out[[h_opt]]
    i_0 <- 0
    for(k in 1:K){
      Us[[k]] <- u[i_0+1:ps[k],1:h_opt,drop=F]
      Bs[[k]] <- B[i_0+1:ps[k],,drop=F]
      i_0 <- i_0+ps[k]
    }
    keepX_sol <- paras_keep[na.omit(id_paras_out),1]
    keepY_sol <- paras_keep[na.omit(id_paras_out),2]
    mu_y <- matrix(rep(colMeans(Y),n) ,nrow = n,byrow = T)
    y_est <- mu_y + x0_center%*%(B*sd_x_mat)
    if(verbose){
      cat(paste("\nComplete model   ",
                "   Q2=",round(Q2_tot_s[h_opt],2)," \n\n",sep=""))
    }
    # Prepare outputs
    optimal_parameters <- list(keepX_sol=keepX_sol,keepY_sol=keepY_sol,R=h_opt,
                               Q2=Q2_tot_s[h_opt])
    parameters <- list(RSS0=RSS0)
    bootstrap <- list(paras_keep=paras_keep,
                      Q2_mean=Q2_mean,
                      Q2_all_mean=Q2_all_mean,
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
                t_h=t_h[,1:h_opt,drop=F],
                explained_variance=na.omit(VAR_h_s)*100,
                explained_variance_cum=na.omit(CUM_VAR_h_s)*100,
                NZV=NZV,
                Bs=Bs,B_r=B_r_out,
                y_est=y_est,parameters=parameters,
                mu_k=c(mu_x_s),mu_y=mu_y,sd_y=apply(Y,2,sd))
    # if(verbose){
    #   plot_res(res)
    # }
  }else{
    res <- NULL
  }
  res
}

#' Title
#'
#' @param Xs Xs
#' @param Y Y
#' @param lambdas lambdas
#' @param NCORES number of cores
#' @param verbose verbose
#'
#' @return
#' @export
#'
#' @useDynLib ddsPLS
Q2_local_ddsPLS <- function(Xs,Y,lambdas=NULL,
                            n_B=20,
                            lowExplainedVariance=0,
                            deflatX=T,NCORES=1,center=T,verbose=F){
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
  Us = Bs = B_r_out = Q2_mean = Q2_all_mean =
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
  # lambdas_h <- seq(0,lambda_max,length.out = N_lambdas)
  if(is.null(lambdas)){
    lambdas_h <- seq(0,1,length.out = 100)
    N_lambdas <- 100
  }else{
    lambdas_h <- lambdas
    N_lambdas <- length(lambdas)
  }
  ncomps = Q2_TOTAL <- rep(NA,N_lambdas)
  x0_deflated = y0_deflated = cov_mean_h = prop_models_ok <- list()
  remove_COV <- matrix(0,q,p)
  B_previous <- matrix(0,p,q)
  V_boot = u_boot = res_measure <- list()
  vars_boot_50 = vars_boot_25 = vars_boot_75 =
    vars_boot_h_50 = vars_boot_h_25 = vars_boot_h_75 =
    q2_boot_50 = q2_boot_25 = q2_boot_75 =
    q2_all_boot_50 = q2_all_boot_25 = q2_all_boot_75 <- list()
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
      out <- bootstrap_pls(X_init=X_init,Y_init=Y_init,B_previous=B_previous,
                           u=u,v=V_phi,h=h,lambdas=lambdas_h,
                           lambda_prev = lambda)
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
      vars_boot_50[[h]] = vars_boot_25[[h]] = vars_boot_75[[h]] =
      vars_boot_h = vars_boot_h_sd_plus = vars_boot_h_sd_moins =
      vars_boot_h_50[[h]] = vars_boot_h_25[[h]] = vars_boot_h_75[[h]] =
      q2_boot = q2_boot_sd_plus = q2_boot_sd_moins =
      q2_boot_50[[h]] = q2_boot_25[[h]] = q2_boot_75[[h]] =
      q2_all_boot = q2_all_boot_sd_plus = q2_all_boot_sd_moins =
      q2_all_boot_50[[h]] = q2_all_boot_25[[h]] = q2_all_boot_75[[h]] =
      q2_boot_mean = q2_all_boot_mean =
      prop_models_ok[[h]] <- rep(0,N_lambdas)
    mod_exists <- rep(0,N_lambdas)
    test_batchas <- TRUE

    i_l <- 1
    probabilities <- c(1/4,0.5,1-1/4)#c(aa,0.5,1-aa)
    while(test_batchas){
      pos <- which(res_measure[[h]][,1]==i_l)
      toto <- na.omit(res_measure[[h]][pos,c(2,3,4,5,8)])
      if(length(toto)>0){
        m <- 1
        # Q2_qtles <- mean(toto[,m])+c(-sd(toto[,m]),0,sd(toto[,m]))#quantile(toto[,1],probs = probabilities)#quantile(toto[,1],probs = probabilities)#quantile(toto[,m],c(1/10,1/2,9/10))#
        q2_boot[i_l] <- mean(toto[,m])#Q2_qtles[2]
        q2_boot_sd_plus[i_l] <- mean(toto[,m])+sd(toto[,m])#Q2_qtles[3]
        q2_boot_sd_moins[i_l] <- mean(toto[,m])-sd(toto[,m])#Q2_qtles[1]
        q2_boot_50[[h]][i_l] <- quantile(toto[,m],probs = 0.5)
        q2_boot_75[[h]][i_l] <- quantile(toto[,m],probs = 0.75)
        q2_boot_25[[h]][i_l] <- quantile(toto[,m],probs = 0.25)
        q2_boot_mean[i_l] <- mean(toto[,m])/(sd(toto[,m]))#(1+sd(toto[,m])/n)
        #Q2_qtles[2]/(Q2_qtles[3]-Q2_qtles[1])#mean(toto[,m])/sd(toto[,m])
        m <- 2
        # Q2_all_qtles <- mean(toto[,m])+c(-sd(toto[,m]),0,sd(toto[,m]))#quantile(toto[,2],probs = probabilities)#quantile(toto[,2],probs = probabilities)
        q2_all_boot[i_l] <- mean(toto[,m])#Q2_all_qtles[2]
        q2_all_boot_sd_plus[i_l] <- mean(toto[,m])+sd(toto[,m])#Q2_all_qtles[3]
        q2_all_boot_sd_moins[i_l] <- mean(toto[,m])-sd(toto[,m])#Q2_all_qtles[1]
        q2_all_boot_50[[h]][i_l] <- quantile(toto[,m],probs = 0.5)
        q2_all_boot_75[[h]][i_l] <- quantile(toto[,m],probs = 0.75)
        q2_all_boot_25[[h]][i_l] <- quantile(toto[,m],probs = 0.25)
        q2_all_boot_mean[i_l] <- mean(toto[,m])/(sd(toto[,m]))#q2_boot[i_l]/(q2_boot_sd_plus[i_l]-q2_boot_sd_moins[i_l])
        m <- 3
        # VARS_qtles <- mean(toto[,m])+c(-sd(toto[,m]),0,sd(toto[,m]))#quantile(toto[,3],probs = probabilities)##quantile(toto[,3],probs = probabilities)
        vars_boot[i_l] <- mean(toto[,m])#VARS_qtles[2]
        vars_boot_sd_plus[i_l] <- mean(toto[,m])+sd(toto[,m])#VARS_qtles[3]
        vars_boot_sd_moins[i_l] <- mean(toto[,m])-sd(toto[,m])#VARS_qtles[1]
        vars_boot_50[[h]][i_l] <- quantile(toto[,m],probs = 0.5)
        vars_boot_75[[h]][i_l] <- quantile(toto[,m],probs = 0.75)
        vars_boot_25[[h]][i_l] <- quantile(toto[,m],probs = 0.25)
        m <- 4
        # VARS_qtles <- mean(toto[,m])+c(-sd(toto[,m]),0,sd(toto[,m]))#quantile(toto[,4],probs = probabilities)##quantile(toto[,4],probs = probabilities)
        vars_boot_h[i_l] <- mean(toto[,m])#VARS_qtles[2]
        vars_boot_h_sd_plus[i_l] <- mean(toto[,m])+sd(toto[,m])#VARS_qtles[3]
        vars_boot_h_sd_moins[i_l] <- mean(toto[,m])-sd(toto[,m])#VARS_qtles[1]
        vars_boot_h_50[[h]][i_l] <- quantile(toto[,m],probs = 0.5)
        vars_boot_h_75[[h]][i_l] <- quantile(toto[,m],probs = 0.75)
        vars_boot_h_25[[h]][i_l] <- quantile(toto[,m],probs = 0.25)
        prop_models_ok[[h]][i_l] <- sum(na.omit(toto[,5]))/length(pos)
      }
      if(is.na(q2_boot[i_l])|is.nan(q2_boot[i_l])|i_l==N_lambdas){#|mod_exists[i_l]!=1){
        test_batchas <- F
      }
      i_l <- i_l + 1
    }
    vars <- vars_boot
    id_ALL_TEST <- which(q2_boot > 0 & vars_boot_h > lowExplainedVariance)
    if(h!=1){
      id_test_h <- which(q2_all_boot>Q2_all_sum_star[[h-1]][which(lambdas_out[[h-1]]==lambda[h-1])        ])
      id_ALL_TEST <- intersect(id_ALL_TEST,id_test_h)
    }
    id_ALL_TEST_h[[h]] <- id_ALL_TEST
    Q2_h_sum_star[[h]] <- q2_boot
    Q2_h_sum_star_sd_moins[[h]] <- q2_boot_sd_moins
    Q2_h_sum_star_sd_plus[[h]] <- q2_boot_sd_plus
    Q2_mean[[h]] <- q2_boot_mean
    Q2_all_sum_star[[h]] <- q2_all_boot
    Q2_all_sum_star_sd_moins[[h]] <- q2_all_boot_sd_moins
    Q2_all_sum_star_sd_plus[[h]] <- q2_all_boot_sd_plus
    Q2_all_mean[[h]] <- q2_all_boot_mean
    lambdas_out[[h]] <- lambdas_h
    vars_h_boot[[h]] <- vars_boot
    vars_h_boot_sd_moins[[h]] <- vars_boot_sd_moins
    vars_h_boot_sd_plus[[h]] <- vars_boot_sd_plus
    vars_h_boot_single[[h]] <- vars_boot_h
    vars_h_boot_single_sd_moins[[h]] <- vars_boot_h_sd_moins
    vars_h_boot_single_sd_plus[[h]] <- vars_boot_h_sd_plus
    if(length(id_ALL_TEST)>0){
      q2_max_h[h] <- max(na.omit(q2_boot_mean[id_ALL_TEST]))
      id_s_cool <- which(q2_boot_mean==q2_max_h[h])
      if(length(id_s_cool)>0){
        best_id_h <- max(1,min(id_s_cool))
        if(length(best_id_h)>0){
          test_h <- q2_boot_mean[best_id_h]>0
          if(h>1){
            best_id_h_before <- which(lambdas_out[[h-1]]==lambda[h-1])
            test2 <-
              Q2_all_sum_star[[h-1]][best_id_h_before]<Q2_all_sum_star[[h]][best_id_h]
            test_h <- test_h & test2
          }
          best_lam_h <- lambdas_h[best_id_h]
          # Maybe best lambda right of previous best lambda
          m_gogo <-  model_PLS(x = x0,y=y0,lam=lambdas_h[best_id_h],to.scale = F)
          u_sol_boot_h <- m_gogo$U_out
          V_optim_boot_h <- m_gogo$V_optim
          V_phi[,h] <- V_optim_boot_h
          t_boot_h <- x0%*%u_sol_boot_h
          norm_t_0 <- sum(t_boot_h^2)
          if(norm_t_0>1e-9){
            bt <- crossprod(t_boot_h,x0)/norm_t_0
            y0_mask <- y0;y0_mask[,which(abs(V_optim_boot_h)<1e-9)] <- 0
            V_sol_boot_h <- crossprod(y0_mask,t_boot_h)#
            # V_sol_boot_h <- V_optim_boot_h%*%crossprod(V_optim_boot_h,crossprod(y0,t_boot_h))
            V_sol_boot_h <- V_sol_boot_h/norm_t_0
          }else{
            V_sol_boot_h <- V_sol_boot_h*0
          }
          B_boot_h <- tcrossprod(u_sol_boot_h,V_sol_boot_h)
          ### End
          sd_y <- apply(Y,2,sd)
          ## Variance per component
          y_est_boot <- x0%*%B_boot_h
          varia_expl <- 1 - sum((Y_init-y_est_boot)^2)/sum(RSS0)
          ##
          if(varia_expl>1e-9){
            VAR_h_s[h] <- varia_expl
            Q2_tot_s[h] <- Q2_all_sum_star[[h]][best_id_h]#c(Q2_tot_s,Q2_all_sum_star[[h]][best_id_h])
            Q2_all_sum_star_boxplot[[h]] <- res_measure[[h]][which(res_measure[[h]][,1]==best_id_h),3]
            R2_all_sum_star_boxplot[[h]] <- res_measure[[h]][which(res_measure[[h]][,1]==best_id_h),5]
            lambda[h] <- best_lam_h
            # Get the regression matrix of the optimal model
            u[,h] <- u_sol_boot_h
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
            t_r <- x0%*%u_sol_boot_h
            t_h[,h] <- t_r
            V_r <- V_sol_boot_h
            if(verbose){
              if(h==1){
                cat("\n ==============================")
              }
              cat(paste("\n Component ",h,
                        "   lambda=",round(best_lam_h,3),
                        "   var.expl._h=",round(varia_expl*100),"%",
                        "   Q2_h=",round(q2_all_boot[best_id_h]*100)/100,
                        sep=""))
            }
            B_r_out[[h]] <- B_boot_h
            V[,h] <- V_r
            y_r <- tcrossprod(t_r,D[,h,drop=F])
            x_r <- tcrossprod(t_r,C[,h,drop=F])
            y0 <- y0 - y_r
            if(deflatX){
              x0 <- x0 - x_r
            }
            B_previous <- B_all
            ## Variance total
            y_est_boot_tot <- X_init%*%B_previous
            cum_vari_h <- 1-sum((Y_init-y_est_boot_tot)^2)/sum(RSS0)
            CUM_VAR_h_s[h] <- cum_vari_h
          }else{
            test <- F
          }
        }else{
          test <- F
        }
      }else{
        test <- F
        cat(paste("\n  --  Component ",h," not built, no accessible parameter.",sep=""))
        cat("\n ==============================")
      }
    }else{
      test <- F
      cat(paste("\n  --  Component ",h," not built, no accessible parameter.",sep=""))
      cat("\n ==============================")
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
      cat(paste("\n Complete model   ",
                "   Q2=",round(Q2_tot_s[h_opt],2),sep=""))
      cat("\n ==============================")
    }
    # Prepare outputs
    optimal_parameters <- list(lambda=lambda_sol,R=h_opt,
                               Q2=Q2_tot_s[h_opt])
    parameters <- list(RSS0=RSS0)
    quartiles <- list(vars_boot_50=vars_boot_50,vars_boot_25=vars_boot_25,vars_boot_75=vars_boot_75,
                      vars_boot_h_50=vars_boot_h_50,vars_boot_h_25=vars_boot_h_25,vars_boot_h_75=vars_boot_h_75,
                      q2_boot_50=q2_boot_50,q2_boot_25=q2_boot_25,q2_boot_75=q2_boot_75,
                      q2_all_boot_50=q2_all_boot_50,q2_all_boot_25=q2_all_boot_25,q2_all_boot_75=q2_all_boot_75)
    bootstrap <- list(lambdas_h=lambdas_out,
                      Q2_mean=Q2_mean,
                      Q2_all_mean=Q2_all_mean,
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
                      cov_mean_h=cov_mean_h,prop_models_ok=prop_models_ok,
                      quartiles = quartiles)
    res <- list(optimal_parameters=optimal_parameters,
                bootstrap=bootstrap,
                Us=Us,V=V[,1:h_opt,drop=F],B_cbind=B,
                id_ALL_TEST_h=id_ALL_TEST_h,
                x0_deflated=x0_deflated,y0_deflated=y0_deflated,
                t_h=t_h[,1:h_opt,drop=F],
                explained_variance=na.omit(VAR_h_s)*100,
                explained_variance_cum=na.omit(CUM_VAR_h_s)*100,
                Bs=Bs,B_r=B_r_out,
                y_est=y_est,parameters=parameters,
                mu_k=c(mu_x_s),mu_y=mu_y,sd_y=apply(Y,2,sd))
    if(verbose){
      if(h_opt>0){
        plot_res(res)
      }
    }
  }else{
    res <- NULL
  }
  res
}


#' Title
#'
#' @param Xs Xs
#' @param Y Y
#' @param lambdas lambdas
#' @param NCORES number of cores
#' @param verbose verbose
#'
#' @return
#' @export
#'
#' @useDynLib ddsPLS
Q2_UNIK_ddsPLS <- function(Xs,Y,lambdas=NULL,
                           n_B=20,
                           lowExplainedVariance=0,
                           deflatX=T,NCORES=1,center=T,verbose=F){
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
  Us = Bs = B_r_out = Q2_mean = Q2_all_mean =
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
  # lambdas_h <- seq(0,lambda_max,length.out = N_lambdas)
  if(is.null(lambdas)){
    lambdas_h <- seq(0,1,length.out = 100)
    N_lambdas <- 100
  }else{
    lambdas_h <- lambdas
    N_lambdas <- length(lambdas)
  }
  ncomps = Q2_TOTAL <- rep(NA,N_lambdas)
  x0_deflated = y0_deflated = cov_mean_h = prop_models_ok <- list()
  remove_COV <- matrix(0,q,p)
  B_previous <- matrix(0,p,q)
  V_boot = u_boot = res_measure <- list()
  vars_boot_50 = vars_boot_25 = vars_boot_75 =
    vars_boot_h_50 = vars_boot_h_25 = vars_boot_h_75 =
    q2_boot_50 = q2_boot_25 = q2_boot_75 =
    q2_all_boot_50 = q2_all_boot_25 = q2_all_boot_75 <- list()
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
      out <- bootstrap_pls(X_init=X_init,Y_init=Y_init,B_previous=B_previous,
                           u=u,v=V_phi,h=h,lambdas=lambdas_h,
                           lambda_prev = lambda)
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
      vars_boot_50[[h]] = vars_boot_25[[h]] = vars_boot_75[[h]] =
      vars_boot_h = vars_boot_h_sd_plus = vars_boot_h_sd_moins =
      vars_boot_h_50[[h]] = vars_boot_h_25[[h]] = vars_boot_h_75[[h]] =
      q2_boot = q2_boot_sd_plus = q2_boot_sd_moins =
      q2_boot_50[[h]] = q2_boot_25[[h]] = q2_boot_75[[h]] =
      q2_all_boot = q2_all_boot_sd_plus = q2_all_boot_sd_moins =
      q2_all_boot_50[[h]] = q2_all_boot_25[[h]] = q2_all_boot_75[[h]] =
      q2_boot_mean = q2_all_boot_mean =
      prop_models_ok[[h]] <- rep(0,N_lambdas)
    mod_exists <- rep(0,N_lambdas)
    test_batchas <- TRUE

    i_l <- 1
    probabilities <- c(1/4,0.5,1-1/4)#c(aa,0.5,1-aa)
    while(test_batchas){
      pos <- which(res_measure[[h]][,1]==i_l)
      toto <- na.omit(res_measure[[h]][pos,c(2,3,4,5,8)])
      if(length(toto)>0){
        m <- 1
        # Q2_qtles <- mean(toto[,m])+c(-sd(toto[,m]),0,sd(toto[,m]))#quantile(toto[,1],probs = probabilities)#quantile(toto[,1],probs = probabilities)#quantile(toto[,m],c(1/10,1/2,9/10))#
        q2_boot[i_l] <- mean(toto[,m])#Q2_qtles[2]
        q2_boot_sd_plus[i_l] <- mean(toto[,m])+sd(toto[,m])#Q2_qtles[3]
        q2_boot_sd_moins[i_l] <- mean(toto[,m])-sd(toto[,m])#Q2_qtles[1]
        q2_boot_50[[h]][i_l] <- quantile(toto[,m],probs = 0.5)
        q2_boot_75[[h]][i_l] <- quantile(toto[,m],probs = 0.75)
        q2_boot_25[[h]][i_l] <- quantile(toto[,m],probs = 0.25)
        q2_boot_mean[i_l] <- mean(toto[,m])/(sd(toto[,m]))#(1+sd(toto[,m])/n)
        #Q2_qtles[2]/(Q2_qtles[3]-Q2_qtles[1])#mean(toto[,m])/sd(toto[,m])
        m <- 2
        # Q2_all_qtles <- mean(toto[,m])+c(-sd(toto[,m]),0,sd(toto[,m]))#quantile(toto[,2],probs = probabilities)#quantile(toto[,2],probs = probabilities)
        q2_all_boot[i_l] <- mean(toto[,m])#Q2_all_qtles[2]
        q2_all_boot_sd_plus[i_l] <- mean(toto[,m])+sd(toto[,m])#Q2_all_qtles[3]
        q2_all_boot_sd_moins[i_l] <- mean(toto[,m])-sd(toto[,m])#Q2_all_qtles[1]
        q2_all_boot_50[[h]][i_l] <- quantile(toto[,m],probs = 0.5)
        q2_all_boot_75[[h]][i_l] <- quantile(toto[,m],probs = 0.75)
        q2_all_boot_25[[h]][i_l] <- quantile(toto[,m],probs = 0.25)
        q2_all_boot_mean[i_l] <- mean(toto[,m])/(sd(toto[,m]))#q2_boot[i_l]/(q2_boot_sd_plus[i_l]-q2_boot_sd_moins[i_l])
        m <- 3
        # VARS_qtles <- mean(toto[,m])+c(-sd(toto[,m]),0,sd(toto[,m]))#quantile(toto[,3],probs = probabilities)##quantile(toto[,3],probs = probabilities)
        vars_boot[i_l] <- mean(toto[,m])#VARS_qtles[2]
        vars_boot_sd_plus[i_l] <- mean(toto[,m])+sd(toto[,m])#VARS_qtles[3]
        vars_boot_sd_moins[i_l] <- mean(toto[,m])-sd(toto[,m])#VARS_qtles[1]
        vars_boot_50[[h]][i_l] <- quantile(toto[,m],probs = 0.5)
        vars_boot_75[[h]][i_l] <- quantile(toto[,m],probs = 0.75)
        vars_boot_25[[h]][i_l] <- quantile(toto[,m],probs = 0.25)
        m <- 4
        # VARS_qtles <- mean(toto[,m])+c(-sd(toto[,m]),0,sd(toto[,m]))#quantile(toto[,4],probs = probabilities)##quantile(toto[,4],probs = probabilities)
        vars_boot_h[i_l] <- mean(toto[,m])#VARS_qtles[2]
        vars_boot_h_sd_plus[i_l] <- mean(toto[,m])+sd(toto[,m])#VARS_qtles[3]
        vars_boot_h_sd_moins[i_l] <- mean(toto[,m])-sd(toto[,m])#VARS_qtles[1]
        vars_boot_h_50[[h]][i_l] <- quantile(toto[,m],probs = 0.5)
        vars_boot_h_75[[h]][i_l] <- quantile(toto[,m],probs = 0.75)
        vars_boot_h_25[[h]][i_l] <- quantile(toto[,m],probs = 0.25)
        prop_models_ok[[h]][i_l] <- sum(na.omit(toto[,5]))/length(pos)
      }
      if(is.na(q2_boot[i_l])|is.nan(q2_boot[i_l])|i_l==N_lambdas){#|mod_exists[i_l]!=1){
        test_batchas <- F
      }
      i_l <- i_l + 1
    }
    vars <- vars_boot
    id_ALL_TEST <- which(q2_boot > 0 & vars_boot_h > lowExplainedVariance)
    if(h!=1){
      id_test_h <- which(q2_all_boot>Q2_all_sum_star[[h-1]][which(lambdas_out[[h-1]]==lambda[h-1])        ])
      id_ALL_TEST <- intersect(id_ALL_TEST,id_test_h)
    }
    id_ALL_TEST_h[[h]] <- id_ALL_TEST
    Q2_h_sum_star[[h]] <- q2_boot
    Q2_h_sum_star_sd_moins[[h]] <- q2_boot_sd_moins
    Q2_h_sum_star_sd_plus[[h]] <- q2_boot_sd_plus
    Q2_mean[[h]] <- q2_boot_mean
    Q2_all_sum_star[[h]] <- q2_all_boot
    Q2_all_sum_star_sd_moins[[h]] <- q2_all_boot_sd_moins
    Q2_all_sum_star_sd_plus[[h]] <- q2_all_boot_sd_plus
    Q2_all_mean[[h]] <- q2_all_boot_mean
    lambdas_out[[h]] <- lambdas_h
    vars_h_boot[[h]] <- vars_boot
    vars_h_boot_sd_moins[[h]] <- vars_boot_sd_moins
    vars_h_boot_sd_plus[[h]] <- vars_boot_sd_plus
    vars_h_boot_single[[h]] <- vars_boot_h
    vars_h_boot_single_sd_moins[[h]] <- vars_boot_h_sd_moins
    vars_h_boot_single_sd_plus[[h]] <- vars_boot_h_sd_plus
    if(length(id_ALL_TEST)>0){
      q2_max_h[h] <- max(na.omit(q2_boot_mean[id_ALL_TEST]))
      id_s_cool <- which(q2_boot_mean==q2_max_h[h])
      if(length(id_s_cool)>0){
        best_id_h <- max(1,min(id_s_cool))
        if(length(best_id_h)>0){
          test_h <- q2_boot_mean[best_id_h]>0
          if(h>1){
            best_id_h_before <- which(lambdas_out[[h-1]]==lambda[h-1])
            test2 <-
              Q2_all_sum_star[[h-1]][best_id_h_before]<Q2_all_sum_star[[h]][best_id_h]
            test_h <- test_h & test2
          }
          best_lam_h <- lambdas_h[best_id_h]
          # Maybe best lambda right of previous best lambda
          m_gogo <-  model_PLS(x = x0,y=y0,lam=lambdas_h[best_id_h],to.scale = F)
          u_sol_boot_h <- m_gogo$U_out
          V_optim_boot_h <- m_gogo$V_optim
          V_phi[,h] <- V_optim_boot_h
          t_boot_h <- x0%*%u_sol_boot_h
          norm_t_0 <- sum(t_boot_h^2)
          if(norm_t_0>1e-9){
            bt <- crossprod(t_boot_h,x0)/norm_t_0
            y0_mask <- y0;y0_mask[which(abs(V_optim_boot_h)<1e-9),] <- 0;y0_mask[which(abs(V_optim_boot_h)>=1e-9),] <- 1
            V_sol_boot_h <- crossprod(y0_mask,t_boot_h)#
            # V_sol_boot_h <- V_optim_boot_h%*%crossprod(V_optim_boot_h,crossprod(y0,t_boot_h))
            V_sol_boot_h <- V_sol_boot_h/norm_t_0
          }else{
            V_sol_boot_h <- V_sol_boot_h*0
          }
          B_boot_h <- tcrossprod(u_sol_boot_h,V_sol_boot_h)
          ### End
          sd_y <- apply(Y,2,sd)
          ## Variance per component
          y_est_boot <- x0%*%B_boot_h
          varia_expl <- 1 - sum((Y_init-y_est_boot)^2)/sum(RSS0)
          ##
          if(varia_expl>1e-9){
            VAR_h_s[h] <- varia_expl
            Q2_tot_s[h] <- Q2_all_sum_star[[h]][best_id_h]#c(Q2_tot_s,Q2_all_sum_star[[h]][best_id_h])
            Q2_all_sum_star_boxplot[[h]] <- res_measure[[h]][which(res_measure[[h]][,1]==best_id_h),3]
            R2_all_sum_star_boxplot[[h]] <- res_measure[[h]][which(res_measure[[h]][,1]==best_id_h),5]
            lambda[h] <- best_lam_h
            lambdas_h <- best_lam_h
            N_lambdas <- 1
            # Get the regression matrix of the optimal model
            u[,h] <- u_sol_boot_h
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
            t_r <- x0%*%u_sol_boot_h
            t_h[,h] <- t_r
            V_r <- V_sol_boot_h
            if(verbose){
              if(h==1){
                cat("\n ==============================")
              }
              cat(paste("\n Component ",h,
                        "   lambda=",round(best_lam_h,3),
                        "   var.expl._h=",round(varia_expl*100),"%",
                        "   Q2_h=",round(q2_all_boot[best_id_h]*100)/100,
                        sep=""))
            }
            B_r_out[[h]] <- B_boot_h
            V[,h] <- V_r
            y0 <- y0-x0%*%B_boot_h
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
            cum_vari_h <- 1-sum((Y_init-y_est_boot_tot)^2)/sum(RSS0)
            CUM_VAR_h_s[h] <- cum_vari_h
          }else{
            test <- F
          }
        }else{
          test <- F
        }
      }else{
        test <- F
        cat(paste("\n  --  Component ",h," not built, no accessible parameter.",sep=""))
        cat("\n ==============================")
      }
    }else{
      test <- F
      cat(paste("\n  --  Component ",h," not built, no accessible parameter.",sep=""))
      cat("\n ==============================")
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
      cat(paste("\n Complete model   ",
                "   Q2=",round(Q2_tot_s[h_opt],2),sep=""))
      cat("\n ==============================")
    }
    # Prepare outputs
    optimal_parameters <- list(lambda=lambda_sol,R=h_opt,
                               Q2=Q2_tot_s[h_opt])
    parameters <- list(RSS0=RSS0)
    quartiles <- list(vars_boot_50=vars_boot_50,vars_boot_25=vars_boot_25,vars_boot_75=vars_boot_75,
                      vars_boot_h_50=vars_boot_h_50,vars_boot_h_25=vars_boot_h_25,vars_boot_h_75=vars_boot_h_75,
                      q2_boot_50=q2_boot_50,q2_boot_25=q2_boot_25,q2_boot_75=q2_boot_75,
                      q2_all_boot_50=q2_all_boot_50,q2_all_boot_25=q2_all_boot_25,q2_all_boot_75=q2_all_boot_75)
    bootstrap <- list(lambdas_h=lambdas_out,
                      Q2_mean=Q2_mean,
                      Q2_all_mean=Q2_all_mean,
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
                      cov_mean_h=cov_mean_h,prop_models_ok=prop_models_ok,
                      quartiles = quartiles)
    res <- list(optimal_parameters=optimal_parameters,
                bootstrap=bootstrap,
                Us=Us,V=V[,1:h_opt,drop=F],B_cbind=B,
                id_ALL_TEST_h=id_ALL_TEST_h,
                x0_deflated=x0_deflated,y0_deflated=y0_deflated,
                t_h=t_h[,1:h_opt,drop=F],
                explained_variance=na.omit(VAR_h_s)*100,
                explained_variance_cum=na.omit(CUM_VAR_h_s)*100,
                Bs=Bs,B_r=B_r_out,
                y_est=y_est,parameters=parameters,
                mu_k=c(mu_x_s),mu_y=mu_y,sd_y=apply(Y,2,sd))
  }else{
    res <- NULL
  }
  res
}

plot_res_init <- function(res){
  h_opt <- res$optimal_parameters$R
  lambdas_out <- res$bootstrap$lambdas_h
  cols <- c(RColorBrewer::brewer.pal(max(h_opt,3),"Set1")[1:h_opt],"gray80")
  layout(matrix(c(1,1,2,2,3,7,7,4,4,5,5,6,7,7), 2, 7, byrow = TRUE))
  par(mar=c(3,3,2,1),mgp=c(2,1,0))
  aa <- min(1/4,1-1/4)
  quart <- res$bootstrap$quartiles
  ls <- res$bootstrap$lambdas_h[[1]]
  plots <- list(
    R2_h=list(q50=quart$vars_boot_h_50,q75=quart$vars_boot_h_75,q25=quart$vars_boot_h_25,
              main=expression("R"["B,h"]^"2"),
              lam_opt=res$optimal_parameters$lambda,
              vars_single=res$bootstrap$vars_h_boot_single),
    R2=list(q50=quart$vars_boot_50,q75=quart$vars_boot_75,q25=quart$vars_boot_25,
            main=expression("R"["B"]^"2"),
            lam_opt=res$optimal_parameters$lambda,
            vars_single=res$bootstrap$vars_boot_single),
    Q2_h=list(q50=quart$q2_boot_50,q75=quart$q2_boot_75,q25=quart$q2_boot_25,
              main=expression("Q"["B,h"]^"2"),
              lam_opt=res$optimal_parameters$lambda,
              vars_single=res$bootstrap$Q2_h_star),
    Q2=list(q50=quart$q2_all_boot_50,q75=quart$q2_all_boot_75,q25=quart$q2_all_boot_25,
            main=expression("Q"["B"]^"2"),
            lam_opt=res$optimal_parameters$lambda,
            vars_single=res$bootstrap$Q2_all_sum_star)
  )
  widths <- (ls[-1]+ls[-length(ls)])/2;widths <- c(ls[1],widths,ls[length(ls)])
  for(i in 1:4){
    plot(0,0,xlim=range(ls),ylim=c(-0.1,1.15),
         main=plots[[i]]$main,ylab="",xlab=expression(lambda),col="white")
    abline(h=0,lty=5,col=1)
    add <- T
    oo<-lapply(
      1:length(ls),
      function(ii){
        points(rep(ls[ii],2),c(plots[[i]]$q25[[h_opt+1]][ii],plots[[i]]$q75[[h_opt+1]][ii]),
               col="gray80",type="l",lwd=1)
        points(c(widths[max(1,ii)],widths[min(length(ls),ii+1)]),
               c(plots[[i]]$q25[[h_opt+1]][ii],plots[[i]]$q25[[h_opt+1]][ii]),
               col="gray80",type="l")
        points(c(widths[max(1,ii)],widths[min(length(ls),ii+1)]),
               c(plots[[i]]$q75[[h_opt+1]][ii],plots[[i]]$q75[[h_opt+1]][ii]),
               col="gray80",type="l")
      })
    points(ls,plots[[i]]$q50[[h_opt+1]],col="gray80",pch=1)
    for(h in (h_opt):1){
      id_h <- res$id_ALL_TEST_h[[h]]
      add <- T
      oo<-lapply(
        1:length(ls),
        function(ii){
          points(rep(ls[ii],2),c(plots[[i]]$q25[[h]][ii],plots[[i]]$q75[[h]][ii]),
                 col="gray50",type="l")
          points(c(widths[max(1,ii)],widths[min(length(ls),ii+1)]),
                 c(plots[[i]]$q25[[h]][ii],plots[[i]]$q25[[h]][ii]),
                 col="gray50",type="l")
          points(c(widths[max(1,ii)],widths[min(length(ls),ii+1)]),
                 c(plots[[i]]$q75[[h]][ii],plots[[i]]$q75[[h]][ii]),
                 col="gray50",type="l")
        })
      points(ls,plots[[i]]$q50[[h]],col="gray50",pch=1)
      if(length(id_h)>0){
        points(ls[id_h],plots[[i]]$q50[[h]][id_h],col=cols[h],pch=16)
      }
      legend(bty="n","topleft",fill=cols,legend = c(unlist(lapply(1:h_opt,function(h){
        paste("Component ",h," (",round(res$explained_variance[h]),"%)",sep="",collapse = "")})),
        "Not selected component"))
    }
    abline(v=res$optimal_parameters$lambda,col=cols,lty=3)
    if(i==2){
      abline(v=res$optimal_parameters$lambda,col=cols,lty=3)
      ## Boxplot R2
      ddff <- data.frame(do.call(rbind,lapply(1:h_opt,function(ii,bb){cbind(ii,bb[[ii]])},
                                              res$bootstrap$R2_all_sum_star_boxplot)))
      names(ddff) <- c("h","R2")
      bobo <- boxplot(R2~h,ddff,ylim=c(-0.1,1.15),border=cols,main=expression("R"["B,h"]^"2"),xlab="Component",ylab="")
      points(1:h_opt,res$explained_variance/100,col=1,pch=18,cex=2)
      abline(h=0,lty=5,col=1)
    }else if(i==4){
      ddff <- data.frame(do.call(rbind,lapply(1:h_opt,function(ii,bb){cbind(ii,bb[[ii]])},
                                              res$bootstrap$Q2_all_sum_star_boxplot)))
      names(ddff) <- c("h","Q2")
      bobo <- boxplot(Q2~h,ddff,ylim=c(-0.1,1.15),border=cols,main=expression("Q"[B]^"2"),xlab="Component",ylab="")
      abline(h=0,lty=5,col=1)
      for(h in 1:(h_opt+1)){
        id_h <- res$id_ALL_TEST_h[[h]]
        if(h!=1){
          points(res$bootstrap$lambdas_h[[h]],res$bootstrap$Q2_mean[[h]],
                 type="l",lty=1,col=cols[h])
          points(res$bootstrap$lambdas_h[[h]][id_h ],res$bootstrap$Q2_mean[[h]][id_h ],
                 type="l",lty=1,col=cols[h],lwd=2)
        }else{
          plot(res$bootstrap$lambdas_h[[h]],res$bootstrap$Q2_mean[[h]],
               type="l",lty=1,col=cols[h],ylim=range(na.omit(unlist(res$bootstrap$Q2_mean))),
               main=expression(mu[r]/sigma[r]),ylab="",xlab=expression(lambda),lwd=1)
          abline(h=0,lty=5,col=1)
          points(res$bootstrap$lambdas_h[[h]][id_h ],res$bootstrap$Q2_mean[[h]][id_h ],
                 type="l",lty=1,col=cols[h],lwd=2)
        }
        abline(v=res$optimal_parameters$lambda,col=cols,lty=3)
      }
    }
  }
}
