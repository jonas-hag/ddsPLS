#' Title
#'
#' @param Xs Xs
#' @param Y Y
#' @param type "CT" and "l1
#' @param paras list
#' @param NCORES number of cores
#' @param verbose verbose
#'
#' @return
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
    for(r in 1:(h-1)){
      if(T){#length(na.omit(lambda_prev))==0){
        m_gogo <-  model_PLS(x = X_r,y=Y_r,lam=lambda_prev[r],
                             to.scale = F)
        u[,r] <- m_gogo$U_out
        v[,r] <- m_gogo$V_optim
      }
      t_r <- X_r%*%u[,r]
      t_all[,r] <- t_r
      bt <- crossprod(t_r,X_r)/sum(t_r^2)
      C[,r] <- t(bt)
      Y_r_mask <- Y_r;Y_r_mask[,which(abs(v[,r])<1e-9)] <- 0#;Y_r_mask[which(abs(v[,r])>=1e-9),] <- 1
      D[,r] <- crossprod(Y_r_mask,t_r)/sum(t_r^2)
      U_star_cl <- u[,1:r,drop=F]%*%solve(crossprod(C[,1:r,drop=F],u[,1:r,drop=F]))
      B_youyou <- tcrossprod(U_star_cl,D[,1:r,drop=F])
      y_r <- tcrossprod(t_r,D[,r,drop=F])
      x_r <- tcrossprod(t_r,C[,r,drop=F])
      # Do deflation
      Y_r <- Y_r - y_r
      X_r <- X_r - x_r
      X_defla[[r+1]] <- X_r
      Y_defla[[r+1]] <- Y_r
    }
  }
  vars_expl_star = vars_expl = vars_expl_h = Q2_star = Q2 = Q2_all = model_exists <- rep(0,N_lambdas)

  u_out <- matrix(0,p,N_lambdas)
  V_optim_phi = V_model <- matrix(0,q,N_lambdas)
  B_next_out <- list()
  for(i_l in 1:N_lambdas){
    m_1 <-  model_PLS(x = X_r,y=Y_r,lam=lambdas[i_l],to.scale = F)
    u_il <- m_1$U_out
    V_il <- m_1$V_optim
    u_out[,i_l] <- u_il
    V_optim_phi[,i_l] <- V_il
    t_r <- X_r%*%u_il
    bt <- crossprod(t_r,X_r)/sum(t_r^2)
    C[,h] <- t(bt)
    Y_r_mask <- Y_r;Y_r_mask[,which(abs(V_il)<1e-9)] <- 0#;Y_r_mask[which(abs(V_il)>=1e-9),] <- 1
    D[,h] <- crossprod(Y_r_mask,t_r)/sum(t_r^2)
    if(h>1){
      u_cur_il <- cbind(u[,1:(h-1)],u_il)
    }else{
      u_cur_il <- u_il
    }
    U_star_cl <- u_cur_il%*%solve(crossprod(C,u_cur_il))
    B_all <- tcrossprod(U_star_cl,D)
    y_r <- tcrossprod(t_r,D[,h])
    x_r <- tcrossprod(t_r,C[,h])

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
#' @param type "CT" and "l1
#' @param paras list
#' @param NCORES number of cores
#' @param verbose verbose
#'
#' @return
#' @export
#'
#' @useDynLib ddsPLS
sparse_PLS_Bootstrap <- function(Xs,Y,type="CT",
                                 paras=NULL,
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
    paras_out_list = vars_h_boot = vars_h_boot_sd_moins = vars_h_boot_sd_plus =
    vars_h_boot_single = vars_h_boot_single_sd_moins = vars_h_boot_single_sd_plus <- list()
  q2_max_h = VAR_h_s = CUM_VAR_h_s = Q2_tot_s <- rep(NA,n)
  B_tot_LOO <- matrix(0,n,sum(ps)*q)

  ### For each 'h' component, look for the best lambda
  test <- T
  y0 <- Y_init
  h <- 0
  id_ALL_TEST_h = Q2_all_sum_star_boxplot = R2_all_sum_star_boxplot <- list()
  K_h <- 1:K
  # Check paras and stuff
  paras_in <- matrix(paras,ncol=1) ; paras_out <- matrix(NA,n,1)
  if(type=="l1"){
    paras_in <- expand.grid(paras)
    paras_out <- matrix(NA,n,2)
  }
  N_paras <- nrow(paras_in)
  x0_deflated = y0_deflated = prop_models_ok <- list()
  remove_COV <- matrix(0,q,p)
  V_boot = u_boot = res_measure <- list()
  vars_boot_50 = vars_boot_25 = vars_boot_75 =
    vars_boot_h_50 = vars_boot_h_25 = vars_boot_h_75 =
    q2_boot_50 = q2_boot_25 = q2_boot_75 =
    q2_all_boot_50 = q2_all_boot_25 = q2_all_boot_75 <- list()
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
      if(type=="l1"){
        out <- bootstrap_sparsePLS(X_init=X_init,Y_init=Y_init,
                                   u=u,v=V_phi,h=h,paras_in=paras_in,
                                   paras_prev = paras_out)
      }else if(type=="CT"){
        out <- bootstrap_pls_CT(X_init=X_init,Y_init=Y_init,
                                u=u,v=V_phi,h=h,lambdas=paras_in,
                                lambda_prev = paras_out)
      }
      res_measure_boot <- cbind(1:N_paras,
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
      prop_models_ok[[h]] <- rep(0,N_paras)
    mod_exists <- rep(0,N_paras)
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
      if(is.na(q2_boot[i_l])|is.nan(q2_boot[i_l])|i_l==N_paras){
        test_batchas <- F
      }
      i_l <- i_l + 1
    }
    vars <- vars_boot
    id_ALL_TEST <- which(q2_boot > lowExplainedVariance)# > 0 & vars_boot_h > lowExplainedVariance)
    if(h!=1){
      id_test_h <- which(q2_all_boot>Q2_all_sum_star[[h-1]][
        which(rowSums(
          (paras_out_list[[h-1]]
           -
             matrix(rep(paras_out[h-1,],N_paras),nrow = N_paras,byrow = T))^2
        )==0)
      ])
      id_ALL_TEST <- intersect(id_ALL_TEST,id_test_h)
    }
    if(type=="CT"){
      max_coco <- max(abs(crossprod(x0,y0)/(n-1)))
      id_ALL_TEST <- intersect(id_ALL_TEST,which(paras_in<max_coco))
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
    paras_out_list[[h]] <- paras_in
    vars_h_boot[[h]] <- vars_boot
    vars_h_boot_sd_moins[[h]] <- vars_boot_sd_moins
    vars_h_boot_sd_plus[[h]] <- vars_boot_sd_plus
    vars_h_boot_single[[h]] <- vars_boot_h
    vars_h_boot_single_sd_moins[[h]] <- vars_boot_h_sd_moins
    vars_h_boot_single_sd_plus[[h]] <- vars_boot_h_sd_plus
    if(length(id_ALL_TEST)>0){
      # q2_max_h[h] <- max(na.omit(q2_boot[id_ALL_TEST]))#max(na.omit(q2_boot_mean[id_ALL_TEST]))#
      # id_s_cool <- which(q2_boot==q2_max_h[h])#which(q2_boot_mean==q2_max_h[h])#
      diff_R2_Q2 <- abs(q2_boot-vars_boot)
      q2_max_h[h] <- min(na.omit(diff_R2_Q2[id_ALL_TEST]))#max(na.omit(q2_boot_mean[id_ALL_TEST]))
      id_s_cool <- which(diff_R2_Q2==q2_max_h[h])#which(q2_boot_mean==q2_max_h[h])
      if(length(id_s_cool)>0){
        best_id_h <- max(1,min(id_s_cool))
        if(length(best_id_h)>0){
          test_h <- q2_boot_mean[best_id_h]>0
          if(h>1){
            best_id_h_before <- which(
              rowSums(
                (paras_out_list[[h-1]]
                 -
                   matrix(rep(paras_out[h-1,],N_paras),nrow = N_paras,byrow = T))^2
              )==0)
            test2 <-
              Q2_all_sum_star[[h-1]][best_id_h_before]<Q2_all_sum_star[[h]][best_id_h]
            test_h <- test_h & test2
          }
          best_paras_h <- paras_in[best_id_h,]
          if(type=="CT"){
            m_gogo <-  model_PLS(x = x0,y=y0,lam=best_paras_h,to.scale = F)
          }else if(type=="l1"){
            kX_r <- p;kY_r <- q
            if(sparse){
              kX_r <- best_paras_h[1];kY_r <- best_paras_h[2]
            }
            m_gogo <- sPLS_lasso(x=x0,y=y0,kx=kX_r,ky=kY_r,
                                 to.scale=F,R=1)
          }
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
            Q2_all_sum_star_boxplot[[h]] <- res_measure[[h]][which(res_measure[[h]][,1]==best_id_h),3]
            R2_all_sum_star_boxplot[[h]] <- res_measure[[h]][which(res_measure[[h]][,1]==best_id_h),5]
            paras_out[h,] <- best_paras_h
            if(verbose){
              if(h==1){
                cat("\n ==============================")
              }
              cat(paste("\n Component ",h,
                        "   paras=",paste(round(best_paras_h,3)),
                        "   var.expl._h=",round(varia_expl*100),"%",
                        "   Q2_h=",round(q2_all_boot[best_id_h]*100)/100,
                        sep=""))
            }
            B_r_out[[h]] <- B_current
            # Do deflation
            y0 <- y0 - y_r
            x0 <- x0 - x_r
            V[,h] <- V_optim_boot_h
            t_h[,h] <- t_r
            # y0 <- y0-x0%*%B_boot_h
            y0_deflated[[h]] <- y0
            x0_deflated[[h]] <- x0
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
    }
    B <- B_r_out[[h_opt]]
    i_0 <- 0
    for(k in 1:K){
      Us[[k]] <- u[i_0+1:ps[k],1:h_opt,drop=F]
      Bs[[k]] <- B[i_0+1:ps[k],,drop=F]
      i_0 <- i_0+ps[k]
    }
    paras_sol <- paras_out[1:h_opt,]
    mu_y <- matrix(rep(colMeans(Y),n) ,nrow = n,byrow = T)
    y_est <- mu_y + x0_center%*%(B*sd_x_mat)

    if(verbose){
      cat(paste("\n Complete model   ",
                "   Q2=",round(Q2_tot_s[h_opt],2),sep=""))
      cat("\n ==============================")
    }
    # Prepare outputs
    optimal_parameters <- list(paras=paras_sol,R=h_opt,
                               Q2=Q2_tot_s[h_opt])
    parameters <- list(RSS0=RSS0)
    quartiles <- list(vars_boot_50=vars_boot_50,vars_boot_25=vars_boot_25,vars_boot_75=vars_boot_75,
                      vars_boot_h_50=vars_boot_h_50,vars_boot_h_25=vars_boot_h_25,vars_boot_h_75=vars_boot_h_75,
                      q2_boot_50=q2_boot_50,q2_boot_25=q2_boot_25,q2_boot_75=q2_boot_75,
                      q2_all_boot_50=q2_all_boot_50,q2_all_boot_25=q2_all_boot_25,q2_all_boot_75=q2_all_boot_75)
    bootstrap <- list(paras_out=paras_out_list,
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
                      Q2_all_sum_star_boxplot=Q2_all_sum_star_boxplot,
                      R2_all_sum_star_boxplot=R2_all_sum_star_boxplot,
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
                x0_deflated=x0_deflated,y0_deflated=y0_deflated,
                t_h=t_h[,1:h_opt,drop=F],
                explained_variance=VAR_h_s[1:h_opt]*100,
                explained_variance_cum=CUM_VAR_h_s[1:h_opt]*100,
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
#' @param res res
#' @param h_opt h_opt
#'
#' @return
#' @export
#'
#' @useDynLib ddsPLS
plot_res <- function(res,h_opt=NULL){
  if(is.null(h_opt)){
    h_opt <- res$optimal_parameters$R
  }
  lambdas_out <- res$bootstrap$paras_out
  cols <- c(RColorBrewer::brewer.pal(max(h_opt,3),"Set1")[1:h_opt],"gray80")
  layout(matrix(c(1,1,3,2,2,8,8,8,
                  4,4,5,5,6,8,8,8,
                  rep(7,5),8,8,8), nrow=3, byrow = TRUE))
  # layout(matrix(c(1,1,3,2,2,4,4,5,5,6), 2, 5, byrow = TRUE))
  # layout(matrix(c(1,1,2,2,3,7,7,4,4,5,5,6,7,7), 2, 7, byrow = TRUE))
  par(mar=c(3,3,2,1),mgp=c(2,1,0))
  aa <- min(1/4,1-1/4)
  quart <- res$bootstrap$quartiles
  ls <- res$bootstrap$paras_out[[1]]
  plots <- list(
    R2_h=list(q50=quart$vars_boot_h_50,q75=quart$vars_boot_h_75,q25=quart$vars_boot_h_25,
              main=expression("R"["B,h"]^"2"),
              lam_opt=res$optimal_parameters$paras,
              vars_single=res$bootstrap$vars_h_boot_single),
    R2=list(q50=quart$vars_boot_50,q75=quart$vars_boot_75,q25=quart$vars_boot_25,
            main=expression("R"["B"]^"2"),
            lam_opt=res$optimal_parameters$paras,
            vars_single=res$bootstrap$vars_boot_single),
    Q2_h=list(q50=quart$q2_boot_50,q75=quart$q2_boot_75,q25=quart$q2_boot_25,
              main=expression("Q"["B,h"]^"2"),
              lam_opt=res$optimal_parameters$paras,
              vars_single=res$bootstrap$Q2_h_star),
    Q2=list(q50=quart$q2_all_boot_50,q75=quart$q2_all_boot_75,q25=quart$q2_all_boot_25,
            main=expression("Q"["B"]^"2"),
            lam_opt=res$optimal_parameters$paras,
            vars_single=res$bootstrap$Q2_all_sum_star)
  )
  widths <- (ls[-1]+ls[-length(ls)])/2;widths <- c(ls[1],widths,ls[length(ls)])
  for(i in 1:4){
    plot(0,0,xlim=range(ls),ylim=c(-0.1,1.15),
         main=plots[[i]]$main,ylab="",xlab="Parameter",col="white")
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
    # points(ls,plots[[i]]$q50[[h_opt+1]],col="gray80",pch=1)
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
    abline(v=res$optimal_parameters$paras,col=cols,lty=3)
    if(i==2){
      abline(v=res$optimal_parameters$paras,col=cols,lty=3)
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

      id_rev <- rev(1:h_opt)
      uu<-res$Us[[1]][,id_rev];uu[which(uu==0)] <- NA
      matplot(uu,pch=id_rev,col="white",xlab="Index",ylab="Weight",main="Weights")
      abline(h=0,col="gray80")
      matplot(uu,pch=id_rev,col=cols[id_rev],add=T)

      f <- do.call(cbind,lapply(
        1:h_opt,
        function(hh){
          id_h <- res$id_ALL_TEST_h[[hh]]
          out <- abs(res$bootstrap$Q2_h_star[[hh]]-res$bootstrap$vars_h_boot[[hh]])
          out[-id_h] <- NA
          out
        }))
      matplot(ls,f,pch=1:h_opt,col=cols,xlab="Parameter",ylab="Weight",
              main=expression("|"~bar("R"["B,h"]^"2")~"-"~bar("Q"["B,h"]^"2")~"|"))
      abline(v=res$optimal_parameters$paras,col=cols,lty=3)


      # for(h in 1:(h_opt+1)){
      #   id_h <- res$id_ALL_TEST_h[[h]]
      #   if(h!=1){
      #     points(res$bootstrap$paras_out[[h]],res$bootstrap$Q2_mean[[h]],
      #            type="l",lty=1,col=cols[h])
      #     points(res$bootstrap$paras_out[[h]][id_h ],res$bootstrap$Q2_mean[[h]][id_h ],
      #            type="l",lty=1,col=cols[h],lwd=2)
      #   }else{
      #     plot(res$bootstrap$paras_out[[h]],res$bootstrap$Q2_mean[[h]],
      #          type="l",lty=1,col=cols[h],ylim=range(na.omit(unlist(res$bootstrap$Q2_mean))),
      #          main=expression(mu[r]/sigma[r]),ylab="",xlab="Parameter",lwd=1)
      #     abline(h=0,lty=5,col=1)
      #     points(res$bootstrap$paras_out[[h]][id_h ],res$bootstrap$Q2_mean[[h]][id_h ],
      #            type="l",lty=1,col=cols[h],lwd=2)
      #   }
      #   abline(v=res$optimal_parameters$paras,col=cols,lty=3)
      # }
    }
  }
}
