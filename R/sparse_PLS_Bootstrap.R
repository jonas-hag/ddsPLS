#' Title
#'
#' @param x x
#' @param y y
#' @param deflatX deflatX
#' @param lam Lambda coefficient
#' @param R R
#' @param to.scale to.scale
#'
#' @return Internal object
#' @export
#'
#' @useDynLib ddsPLS
model_PLS <- function(x,y,lam,deflatX=T,R=1,#remove_COV=NULL,
                      to.scale=T){
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
    x_init <- scaleRcpp(x)
    y_init <- scaleRcpp(y)
  }else{
    II <- matrix(1,ncol = n)
    mu_y <- crossprod(II,II%*%y)
    mu_x <- crossprod(II,II%*%x)
    # x_init <- x-mu_x
    # y_init <- y-mu_y
    x_init <- scale(x,scale = F)
    y_init <- scale(y,scale = F)
    mu_y <- matrix(0,n,q)
    # sd_y_mat <- matrix(rep(apply(y,2,sd),p),ncol = q,byrow = T)
    sd_y_mat <- crossprod(II,sqrt(colSums((y-crossprod(II,II%*%y)/n )^2/n)))
    # sd_y_mat <- matrix(1,nrow = n)%*%apply(y,2,sd)
    # sd_y_mat <- get_sd_matrixRcpp(y)
  }
  x0 <- x_init
  y0 <- y_init
  # if(is.null(RSS0)){
  #   RSS0 <- sum(sd_y_mat^2)
  # }else{
  #   RSS0 <- sum(RSS0)
  # }
  # RSS0 <- sum(sd_y_mat^2)
  # var_y_init <- sum(y^2)
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
    COV <- crossprod(y0,x0)/(n-1)
    abs_COV <- abs(COV)
    max_COV <- max(na.omit(c(abs_COV)))
    lam_r <- lam
    if(length(lam)>1) lam_r <- lam[r]
    if(lam_r<max_COV){
      c_h <- do_one_componentCpp(x0 = x0,y0 = y0,COV = COV,lam = lam_r)
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
  list(U_out=U_out,U_star=U_star,
       V_out=V_out,V_optim=V0,
       B=B,B_r=B_r,var_expl=var_expl,#covs=covs,
       score_x=score_x,y_est=y_est,bXr=bXr,bYr=bYr,e_x=x0,e_y=y0,
       x_init=x_init,y_init=y_init)
}

#' Title
#'
#' @param X_init matrix of covariates
#' @param Y_init matrix of response
#' @param h number of components to build
#' @param lambdas vector
#' @param lambda_prev lambdas used to build the previous components
#' @param u X weights of the previous components
#' @param v Y weights of the previous components
#'
#' @return Internal object
#' @export
#'
#' @useDynLib ddsPLS
bootstrap_pls_CT <- function(X_init,Y_init,h=1,lambdas=0,
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
  C <- matrix(0,p,h)
  D <- matrix(0,q,h)
  B_youyou <- matrix(0,p,q)
  X_r <- X_train;Y_r <- Y_train
  if(length(na.omit(lambda_prev))==0){
    u <- matrix(NA,p,n)
    v <- matrix(NA,q,n)
  }
  if(h>1){
    U_reconstruct[,1:(h-1)] <- u[,1:(h-1)]
    r <- 0
    test_no_null <- T
    while(test_no_null){
      r <- r + 1
      m_gogo <-  model_PLS(x = X_r,y=Y_r,lam=lambda_prev[r],to.scale = F)
      # m_gogo <- modelddsPLSCpp(x=X_r,y=Y_r,lam = lambda_prev[r])
      u[,r] <- m_gogo$U_out
      v[,r] <- m_gogo$V_optim
      t_r <- X_r%*%u[,r]
      t_all[,r] <- t_r
      norm_2 <- sum(t_r^2)
      if(norm_2>1e-9){
        bt <- crossprod(t_r,X_r)/norm_2
        C[,r] <- t(bt)
        Y_r_mask <- Y_r;Y_r_mask[,which(abs(v[,r])<1e-9)] <- 0
        D[,r] <- crossprod(Y_r_mask,t_r)/norm_2
        # U_star_cl  <- u[,1:r,drop=F]%*%solve(crossprod(C[,1:r,drop=F],u[,1:r,drop=F]))
        U_star_cl <- u[,1:r,drop=F]%*%solve(crossprod(C[,1:r,drop=F],u[,1:r,drop=F]))
        B_youyou <- tcrossprod(U_star_cl,D[,1:r,drop=F])
        y_r <- tcrossprod(t_r,D[,r,drop=F])
        x_r <- tcrossprod(t_r,C[,r,drop=F])
        # Do deflation
        Y_r <- Y_r - y_r
        X_r <- X_r - x_r
        X_defla[[r+1]] <- X_r
        Y_defla[[r+1]] <- Y_r
      }else{
        test_no_null <- F
      }
      if(r==h-1) test_no_null <- F
    }
  }else{
    norm_2 <- 1
  }
  test_previous_ok <- norm_2>1e-9
  vars_expl_star = vars_expl = vars_expl_h = Q2_star = Q2 = Q2_all = model_exists <- rep(0,N_lambdas)
  u_out <- matrix(0,p,N_lambdas)
  V_il <- matrix(0,q,1)
  V_optim_phi = V_model <- matrix(0,q,N_lambdas)
  B_next_out <- list()
  B_all <- B_youyou
  for(i_l in 1:N_lambdas){
    if(test_previous_ok){
      m_1 <- model_PLS(x = X_r,y=Y_r,lam=lambdas[i_l],to.scale = F)
      # m_1 <- modelddsPLSCpp(x=X_r,y=Y_r,lam = lambdas[i_l])
      u_il <- m_1$U_out
      V_il <- m_1$V_optim
      u_out[,i_l] <- u_il
      V_optim_phi[,i_l] <- V_il
      t_r <- X_r%*%u_il
      norm_2 <- sum(t_r^2)
      if(norm_2>1e-9){
        bt <- crossprod(t_r,X_r)/norm_2
        C[,h] <- t(bt)
        Y_r_mask <- Y_r;Y_r_mask[,which(abs(V_il)<1e-9)] <- 0
        D[,h] <- crossprod(Y_r_mask,t_r)/sum(t_r^2)
        if(h>1){
          u_cur_il <- cbind(u[,1:(h-1)],u_il)
        }else{
          u_cur_il <- u_il
        }
        # U_star_cl  <- u_cur_il%*%solve(crossprod(C,u_cur_il))
        U_star_cl <- u_cur_il%*%solve(crossprod(C,u_cur_il))
      }else{
        U_star_cl <- matrix(0,p,h)
      }
      B_all <- tcrossprod(U_star_cl,D)
      y_r <- tcrossprod(t_r,D[,h])
      x_r <- tcrossprod(t_r,C[,h])
    }

    y_train_pred <- X_train%*%B_all
    y_test_pred <- X_test_normalize%*%B_all
    # Previous components
    y_train_pred_next <- X_train%*%(B_all-B_youyou)
    y_test_pred_RSS <- X_test_normalize%*%B_youyou
    # Compute criterions
    vars_expl[i_l] <- 1 - sum( (Y_train-y_train_pred)^2 ) / sum( (Y_train)^2 )
    vars_expl_h[i_l] <- 1 - sum( (Y_train-y_train_pred_next)^2 ) / sum((Y_train)^2)
    Q2[i_l] <- 1 - sum( (Y_test-y_test_pred)^2 ) / sum((Y_test-y_test_pred_RSS)^2)
    MU_test <- matrix(rep(mu_y,length(id_OOB)),nrow = length(id_OOB),byrow = T)
    Q2_all[i_l] <- 1 - sum( (Y_test-y_test_pred)^2 ) / sum((Y_test)^2)
    if(length(which(abs(V_il)>1e-9))>0) model_exists[i_l] <- 1
  }
  list(u_out=u_out,V_optim_phi=V_optim_phi,id_OOB=id_OOB,model_exists=model_exists,
       vars_expl=vars_expl,vars_expl_h=vars_expl_h,Q2=Q2,Q2_all=Q2_all)
}

#' Title
#'
#' @param Xs Xs
#' @param Y Y
#' @param lambdas list
#' @param n_B NUmber of Bootstrap samples
#' @param lowExplainedVariance Q_2 value above which component is not build. Default to "0".
#' @param deflatX wether or not to deflat X. Default to TRUE. Do not change unless good reasons.
#' @param center wether or not to center X and Y. Default to TRUE. Do not change unless good reasons.
#' @param NCORES number of cores
#' @param verbose verbose
#'
#' @return A "sparse_PLS_Bootstrap" object with results of bootstrap simulations
#'
#' @importFrom stats quantile
#'
#' @export
#'
#' @useDynLib ddsPLS
sparse_PLS_Bootstrap <- function(Xs,Y,
                                 lambdas=NULL,
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
  RSS0 = RSS_h_moins_1 <- colSums(Y_init^2)
  PRESS_y = RSS_y = Q2_y <- matrix(NA,n,q)
  PRESS_y_star = RSS_y_star = Q2_y_star <- matrix(NA,n,q)
  Q2_sum_star <- rep(NA,n)
  ps <- unlist(lapply(Xs_init,ncol))
  X_init <- do.call(cbind,Xs_init)
  x0 <- X_init
  u <- matrix(NA,sum(ps),n)
  V_phi = V <- matrix(NA,q,n)
  P <- matrix(0,sum(ps),n)
  C <- matrix(0,sum(ps),n)
  D <- matrix(0,q,n)
  B <- matrix(0,sum(ps),q)
  t_h <- matrix(NA,n,n)
  Us = Bs = B_r_out = Q2_mean = Q2_all_mean =
    Q2_h_sum_star = Q2_h_sum_star_sd_moins = Q2_h_sum_star_sd_plus =
    Q2_all_sum_star = Q2_all_sum_star_sd_moins = Q2_all_sum_star_sd_plus =
    lambdas_out_list = vars_h_boot = vars_h_boot_sd_moins = vars_h_boot_sd_plus =
    vars_h_boot_single = vars_h_boot_single_sd_moins = vars_h_boot_single_sd_plus <- list()
  q2_max_h = VAR_h_s = CUM_VAR_h_s = Q2_tot_s <- rep(NA,n)
  B_tot_LOO <- matrix(0,n,sum(ps)*q)

  ### For each 'h' component, look for the best lambda
  test <- T
  y0 <- Y_init
  h <- 0
  id_ALL_TEST_h = Q2_all_boxplot = Q2_h_boxplot = R2_h_boxplot = R2_all_boxplot <- list()
  K_h <- 1:K
  # Check lambdas and stuff
  lambdas_in <- matrix(lambdas,ncol=1) ; lambdas_out <- matrix(NA,n,1)
  N_lambdas <- nrow(lambdas_in)
  # x0_deflated = y0_deflated <- list()
  prop_models_ok <- list()
  # remove_COV <- matrix(0,q,p)
  V_boot = u_boot = res_measure <- list()
  vars_boot_50 = vars_boot_25 = vars_boot_75 =
    vars_boot_h_50 = vars_boot_h_25 = vars_boot_h_75 =
    q2_boot_50 = q2_boot_25 = q2_boot_75 =
    q2_all_boot_50 = q2_all_boot_25 = q2_all_boot_75 =
    Q2_Boot_y <- list()
  while(test){
    h <- h + 1
    Bis <- list()
    pos_decoupe <- 1
    Q2_Boot_y[[h]] <- matrix(0,N_lambdas,q)
    i_B <- 1
    NCORES_w <- min(NCORES,n_B)
    `%my_do%` <- ifelse(NCORES_w!=1,{
      out<-`%dopar%`;cl <- makeCluster(NCORES_w)
      registerDoParallel(cl);out},{out <- `%do%`;out})
    Q2_star_bootstrat <- foreach(i_B=1:n_B,.packages = "ddsPLS",.combine='c',.multicombine=TRUE) %my_do% {
      out <- bootstrap_pls_CT(X_init=X_init,Y_init=Y_init,
                              u=u,v=V_phi,h=h,lambdas=lambdas_in,
                              lambda_prev = lambdas_out)
      res_measure_boot <- cbind(1:N_lambdas,
                                out$Q2,
                                out$Q2_all,
                                out$vars_expl,
                                out$vars_expl_h,
                                length(out$id_OOB),
                                i_B,
                                out$model_exists)
      list(res=res_measure_boot,V_optim_phi=out$V_optim_phi,u=out$u_out)#list(cov=out$cov,,res=res_measure,V_model=out$V_model,u=out$u_out,V=out$V_out)
    }
    if(NCORES_w!=1)stopCluster(cl)
    KK <- length(Q2_star_bootstrat)/n_B ; id_pos_boot <- (1:n_B)*KK
    res_measure[[h]] <- do.call(rbind,Q2_star_bootstrat[id_pos_boot-2])
    V_boot[[h]] <- Q2_star_bootstrat[id_pos_boot-1]
    u_boot[[h]] <- Q2_star_bootstrat[id_pos_boot]
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
        q2_boot[i_l] <- mean(toto[,m])
        q2_boot_sd_plus[i_l] <- mean(toto[,m])+sd(toto[,m])
        q2_boot_sd_moins[i_l] <- mean(toto[,m])-sd(toto[,m])
        q2_boot_50[[h]][i_l] <- quantile(toto[,m],probs = 0.5)
        q2_boot_75[[h]][i_l] <- quantile(toto[,m],probs = 0.75)
        q2_boot_25[[h]][i_l] <- quantile(toto[,m],probs = 0.25)
        q2_boot_mean[i_l] <- mean(toto[,m])/(sd(toto[,m]))
        m <- 2
        q2_all_boot[i_l] <- mean(toto[,m])#
        q2_all_boot_sd_plus[i_l] <- mean(toto[,m])+sd(toto[,m])
        q2_all_boot_sd_moins[i_l] <- mean(toto[,m])-sd(toto[,m])
        q2_all_boot_50[[h]][i_l] <- quantile(toto[,m],probs = 0.5)
        q2_all_boot_75[[h]][i_l] <- quantile(toto[,m],probs = 0.75)
        q2_all_boot_25[[h]][i_l] <- quantile(toto[,m],probs = 0.25)
        q2_all_boot_mean[i_l] <- mean(toto[,m])/(sd(toto[,m]))
        m <- 3
        vars_boot[i_l] <- mean(toto[,m])
        vars_boot_sd_plus[i_l] <- mean(toto[,m])+sd(toto[,m])
        vars_boot_sd_moins[i_l] <- mean(toto[,m])-sd(toto[,m])
        vars_boot_50[[h]][i_l] <- quantile(toto[,m],probs = 0.5)
        vars_boot_75[[h]][i_l] <- quantile(toto[,m],probs = 0.75)
        vars_boot_25[[h]][i_l] <- quantile(toto[,m],probs = 0.25)
        m <- 4
        vars_boot_h[i_l] <- mean(toto[,m])
        vars_boot_h_sd_plus[i_l] <- mean(toto[,m])+sd(toto[,m])
        vars_boot_h_sd_moins[i_l] <- mean(toto[,m])-sd(toto[,m])
        vars_boot_h_50[[h]][i_l] <- quantile(toto[,m],probs = 0.5)
        vars_boot_h_75[[h]][i_l] <- quantile(toto[,m],probs = 0.75)
        vars_boot_h_25[[h]][i_l] <- quantile(toto[,m],probs = 0.25)
        prop_models_ok[[h]][i_l] <- sum(na.omit(toto[,5]))/length(pos)
      }
      if(is.na(q2_boot[i_l])|is.nan(q2_boot[i_l])|i_l==N_lambdas){
        test_batchas <- F
      }
      i_l <- i_l + 1
    }
    vars <- vars_boot
    id_ALL_TEST <- which(q2_boot > lowExplainedVariance)# > 0 & vars_boot_h > lowExplainedVariance)
    if(h!=1){
      id_test_h <- which(q2_all_boot>Q2_all_sum_star[[h-1]][
        which(rowSums(
          (lambdas_out_list[[h-1]]
           -
             matrix(rep(lambdas_out[h-1,],N_lambdas),nrow = N_lambdas,byrow = T))^2
        )==0)
      ])
      id_ALL_TEST <- intersect(id_ALL_TEST,id_test_h)
    }
    max_coco <- max(abs(crossprod(x0,y0)/(n-1)))
    id_ALL_TEST <- intersect(id_ALL_TEST,which(lambdas_in<max_coco))
    id_ALL_TEST_h[[h]] <- id_ALL_TEST
    Q2_h_sum_star[[h]] <- q2_boot
    Q2_h_sum_star_sd_moins[[h]] <- q2_boot_sd_moins
    Q2_h_sum_star_sd_plus[[h]] <- q2_boot_sd_plus
    Q2_mean[[h]] <- q2_boot_mean
    Q2_all_sum_star[[h]] <- q2_all_boot
    Q2_all_sum_star_sd_moins[[h]] <- q2_all_boot_sd_moins
    Q2_all_sum_star_sd_plus[[h]] <- q2_all_boot_sd_plus
    Q2_all_mean[[h]] <- q2_all_boot_mean
    lambdas_out_list[[h]] <- lambdas_in
    vars_h_boot[[h]] <- vars_boot
    vars_h_boot_sd_moins[[h]] <- vars_boot_sd_moins
    vars_h_boot_sd_plus[[h]] <- vars_boot_sd_plus
    vars_h_boot_single[[h]] <- vars_boot_h
    vars_h_boot_single_sd_moins[[h]] <- vars_boot_h_sd_moins
    vars_h_boot_single_sd_plus[[h]] <- vars_boot_h_sd_plus
    if(length(id_ALL_TEST)>0){
      # q2_max_h[h] <- max(na.omit(q2_boot[id_ALL_TEST]))#max(na.omit(q2_boot_mean[id_ALL_TEST]))#
      # id_s_cool <- which(q2_boot==q2_max_h[h])#which(q2_boot_mean==q2_max_h[h])#
      diff_R2_Q2 <- vars_boot-q2_all_boot#abs(vars_boot-q2_all_boot)#abs(q2_boot-vars_boot_h)#vars_boot)
      q2_max_h[h] <- min(na.omit(diff_R2_Q2[id_ALL_TEST]))#max(na.omit(q2_boot_mean[id_ALL_TEST]))
      id_s_cool <- which(diff_R2_Q2==q2_max_h[h])#which(q2_boot_mean==q2_max_h[h])
      # oo <- vars_boot_h-q2_boot
      # dd <- sign(oo[-1])-sign(oo[-length(oo)])
      # id_s_cool <- which(dd!=0)
      if(length(id_s_cool)>0){
        best_id_h <- max(1,min(id_s_cool))#id_s_cool#
        if(length(best_id_h)>0){
          test_h <- q2_boot_mean[best_id_h]>0
          if(h>1){
            best_id_h_before <- which(
              rowSums(
                (lambdas_out_list[[h-1]]
                 -
                   matrix(rep(lambdas_out[h-1,],N_lambdas),nrow = N_lambdas,byrow = T))^2
              )==0)
            test2 <-
              Q2_all_sum_star[[h-1]][best_id_h_before]<Q2_all_sum_star[[h]][best_id_h]
            test_h <- test_h & test2
            best_id_h <- max(best_id_h[which(test_h)])
          }
          best_lambdas_h <- lambdas_in[best_id_h,]
          m_gogo <-  model_PLS(x = x0,y=y0,lam=best_lambdas_h,to.scale = F)
          # m_gogo <- modelddsPLSCpp(x=x0,y=y0,lam = best_lambdas_h)
          u_sol_boot_h <- m_gogo$U_out
          V_optim_boot_h <- m_gogo$V_optim
          u[,h] <- u_sol_boot_h
          if(sum(u_sol_boot_h^2)>1e-9){
            t_r <- x0%*%u_sol_boot_h
            bt <- crossprod(t_r,x0)/sum(t_r^2)
            C[,h] <- t(bt)
            Y_r_mask <- y0;Y_r_mask[,which(abs(V_optim_boot_h)<1e-9)] <- 0
            D[,h] <- crossprod(Y_r_mask,t_r)/sum(t_r^2)
            # D[,h] <- crossprod(y0,t_r)/sum(t_r^2)
            U_star_cl <- u[,1:h,drop=F]%*%solve(crossprod(C[,1:h,drop=F],u[,1:h,drop=F]))
            B_current <- tcrossprod(U_star_cl,D[,1:h,drop=F])
            y_r <- tcrossprod(t_r,D[,h,drop=F])
            x_r <- tcrossprod(t_r,C[,h,drop=F])
            ### End
            sd_y <- apply(Y,2,sd)
            ## Variance per component
            varia_expl <- 1 - sum((Y_init-y_r)^2)/sum(RSS0)
            ##
            VAR_h_s[h] <- varia_expl
            Q2_tot_s[h] <- Q2_all_sum_star[[h]][best_id_h]#c(Q2_tot_s,Q2_all_sum_star[[h]][best_id_h])
            Q2_all_boxplot[[h]] <- res_measure[[h]][which(res_measure[[h]][,1]==best_id_h),3]
            Q2_h_boxplot[[h]] <- res_measure[[h]][which(res_measure[[h]][,1]==best_id_h),2]
            R2_h_boxplot[[h]] <- res_measure[[h]][which(res_measure[[h]][,1]==best_id_h),4]
            R2_all_boxplot[[h]] <- res_measure[[h]][which(res_measure[[h]][,1]==best_id_h),5]
            lambdas_out[h,] <- best_lambdas_h
            if(verbose){
              if(h==1){
                cat("\n ==============================")
              }
              cat(paste("\n Component ",h,
                        "   lambdas=",paste(round(best_lambdas_h,3)),
                        "   var.expl._h=",round(varia_expl*100),"%",
                        "   Q2_h=",round(q2_boot[best_id_h]*100)/100,
                        sep=""))
            }
            B_r_out[[h]] <- B_current
            # Do deflation
            y0 <- y0 - y_r
            x0 <- x0 - x_r
            V[,h] <- V_optim_boot_h
            t_h[,h] <- t_r
            # y0 <- y0-x0%*%B_boot_h
            # y0_deflated[[h]] <- y0
            # x0_deflated[[h]] <- x0
            ## Variance total
            y_est_boot_tot <- X_init%*%B_current
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
        if(verbose){
          cat(paste("\n  --  Component ",h," not built, no accessible parameter.",sep=""))
          cat("\n ==============================")
        }
      }
    }else{
      test <- F
      if(verbose){
        cat(paste("\n  --  Component ",h," not built, no accessible parameter.",sep=""))
        cat("\n ==============================")
      }
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
    lambdas_sol <- lambdas_out[1:h_opt,]
    mu_y <- matrix(rep(colMeans(Y),n) ,nrow = n,byrow = T)
    y_est <- mu_y + x0_center%*%(B*sd_x_mat)

    if(verbose){
      cat(paste("\n Complete model   ",
                "   Q2=",round(Q2_tot_s[h_opt],2),sep=""))
      cat("\n ==============================")
    }
    # Prepare outputs
    optimal_parameters <- list(lambdas=lambdas_sol,R=h_opt,
                               Q2=Q2_tot_s[h_opt])
    parameters <- list(RSS0=RSS0)
    quartiles <- list(vars_boot_50=vars_boot_50,vars_boot_25=vars_boot_25,vars_boot_75=vars_boot_75,
                      vars_boot_h_50=vars_boot_h_50,vars_boot_h_25=vars_boot_h_25,vars_boot_h_75=vars_boot_h_75,
                      q2_boot_50=q2_boot_50,q2_boot_25=q2_boot_25,q2_boot_75=q2_boot_75,
                      q2_all_boot_50=q2_all_boot_50,q2_all_boot_25=q2_all_boot_25,q2_all_boot_75=q2_all_boot_75)
    bootstrap <- list(lambdas_out=lambdas_out_list,
                      Q2_mean=Q2_mean,
                      Q2_all_mean=Q2_all_mean,
                      Q2_tot_s=Q2_tot_s[1:h_opt],
                      q2_max_h=q2_max_h[1:h_opt],
                      Q2_h_star=Q2_h_sum_star,
                      Q2_h_star_sd_plus=Q2_h_sum_star_sd_plus,
                      Q2_h_star_sd_moins=Q2_h_sum_star_sd_moins,
                      Q2_all_sum_star=Q2_all_sum_star,
                      Q2_all_sum_star_sd_plus=Q2_all_sum_star_sd_plus,
                      Q2_all_sum_star_sd_moins=Q2_all_sum_star_sd_moins,
                      R2_all_boxplot=R2_all_boxplot,
                      R2_h_boxplot=R2_h_boxplot,
                      Q2_all_boxplot=Q2_all_boxplot,
                      Q2_h_boxplot=Q2_h_boxplot,
                      vars_h_boot=vars_h_boot,
                      vars_h_boot_sd_plus=vars_h_boot_sd_plus,
                      vars_h_boot_sd_moins=vars_h_boot_sd_moins,
                      vars_h_boot_single=vars_h_boot_single,
                      vars_h_boot_single_sd_plus=vars_h_boot_single_sd_plus,
                      vars_h_boot_single_sd_moins=vars_h_boot_single_sd_moins,
                      prop_models_ok=prop_models_ok,
                      quartiles = quartiles)
    res <- list(optimal_parameters=optimal_parameters,
                bootstrap=bootstrap,
                Us=Us,V=V[,1:h_opt,drop=F],B_cbind=B,
                id_ALL_TEST_h=id_ALL_TEST_h,
                # x0_deflated=x0_deflated,y0_deflated=y0_deflated,
                t_h=t_h[,1:h_opt,drop=F],
                explained_variance=VAR_h_s[1:h_opt]*100,
                explained_variance_cum=CUM_VAR_h_s[1:h_opt]*100,
                Bs=Bs,B_r=B_r_out,
                y_est=y_est,parameters=parameters,
                mu_k=c(mu_x_s),mu_y=mu_y,sd_y=apply(Y,2,sd))
    class(res) <- "sparse_PLS_Bootstrap"
    if(verbose){
      if(h_opt>0){
        plot(res)
      }
    }
  }else{
    res <- NULL
    class(res) <- "sparse_PLS_Bootstrap"
  }
  res
}
