get_desgin_A <- function(p0=0.1,d=20,p=10,N_env=20){
  A <- matrix(0,d,p)
  for(i in 1:p){
    TEST <- T
    while(TEST){
      sol <- NULL
      test <- T;iter <- 1
      while(test){
        evenements <- seq(-1,1,length.out = N_ev)
        sol <- c(sol,sample(evenements,size = 1)*rbinom(n = 1,size=1,prob = p0))
        iter <- iter + 1
        if(iter==d+1) test <- F
      }
      TEST <- sum(sol^2)==0
    }
    A[,i] <- sol/sqrt(sum(sol^2))
  }
  A
}

get_A_C <- function(p0=0.1,d=20,p=10,N_env=20){
  A <- get_desgin_A(p0=p0,d=d,p=p,N_env=N_env)
  C <- get_desgin_A(p0=p0,d=d,p=q,N_env=N_env)
  Phi_sel_Y <- which(rowSums(abs(C))!=0)
  Phi_sel_X <- which(rowSums(abs(A))!=0)
  Phi_common_X_Y <- intersect(Phi_sel_Y,Phi_sel_X)
  if(length(Phi_common_X_Y)>0){
    X_sel_model <- which(colSums(abs(A[Phi_common_X_Y,,drop=F]))!=0)
    Y_sel_model <- which(colSums(abs(C[Phi_common_X_Y,,drop=F]))!=0)
  }else{
    X_sel_model = Y_sel_model <- NULL
  }
  list(A=A,C=C,selected=list(x=X_sel_model,y=Y_sel_model))
}

compare_selection <- function(B_hat,B_th){
  vY_th <- which(colSums(abs(B_th))>1e-9)
  vY_no_sel_th <- which(colSums(abs(B_th))<1e-9)
  vX_th <- which(rowSums(abs(B_th))>1e-9)
  vX_no_sel_th <- which(rowSums(abs(B_th))<1e-9)
  vY_hat <- which(colSums(abs(B_hat))>1e-9)
  vX_hat <- which(rowSums(abs(B_hat))>1e-9)
  c(length(intersect(vX_th,vX_hat)),length(vX_hat),
    length(intersect(vY_th,vY_hat)),length(vY_hat))
}


plotCor <- function(){
  mat_cor <- data.frame(cor(cbind(do.call(cbind,Xs),Y)))
  namesX <- unlist(lapply(1:length(Xs),function(i,Xs){paste("X","_Block",i,"_Var",1:ncol(Xs[[i]]),sep = "")},Xs))
  namesY <- paste("Y","_Var",1:ncol(Y),sep = "")
  colnames(mat_cor) =rownames(mat_cor) <- c(namesX,namesY)
  par(mfrow=c(1,1))
  corrplot::corrplot(as.matrix(mat_cor),method = "number" )
}

do_sPLS_one_comp <- function(x,y,Ks=NULL,y_0,kYs=NULL,NCORES=1,tau=0.0975,sparse=T){
  n <- nrow(y);p <- ncol(x);q <- ncol(y)
  kfolds <- n
  fold <- 1:n
  if(!sparse) Ks <- p
  if(!sparse) kYs <- q

  paras <- expand.grid(1:n,Ks,kYs)
  colnames(paras) <- c("n","keepX","keepY")
  n_paras <- nrow(paras)
  if(NCORES>n_paras){
    decoupe <- 1:n_paras
  }else{
    decoupe <- replicate(n_paras/NCORES + 1,
                         sample(1:NCORES))[1:n_paras]
  }
  NCORES_w <- min(NCORES,nrow(paras))
  `%my_do%` <- ifelse(NCORES_w!=1,{
    out<-`%dopar%`
    cl <- makeCluster(NCORES_w)
    registerDoParallel(cl)
    out},{
      out <- `%do%`
      out})
  pos_decoupe <- 1
  PRESS_il_y_and_i_l <- foreach(pos_decoupe=1:min(NCORES,n_paras),.packages = c("mixOmics")) %my_do% {
    paras_here_pos <- which(decoupe==pos_decoupe)
    paras_here <- paras[paras_here_pos,,drop=F]
    n_i <- nrow(paras_here)
    SE <- matrix(NA,n_i,q)
    var_sel <- matrix(0,n_i,q)
    for(i in 1:n_i){
      i_fold <- paras_here[i,1]
      kX <- paras_here[i,2]
      kY <- paras_here[i,3]
      pos_train <- which(fold!=i_fold)
      X_train <- x[pos_train,,drop=FALSE]
      X_test <- x[-pos_train,,drop=FALSE]
      Y_train <- y[pos_train,,drop=FALSE]
      Y_test <- y[-pos_train,,drop=FALSE]
      if(sparse){
        model <- mixOmics::spls(X_train,Y_train,ncomp = 1,keepX = kX,keepY = kY)
      }else{
        model <- mixOmics::pls(X_train,Y_train,ncomp = 1)
      }
      B_est <- predict(model,newdata = X_train)$B.hat[,,1]
      Y_est <- predict(model,newdata = X_test)$predict[,,1]
      SE[i,] <- (Y_est-Y_test)^2
      var_sel[i,which(colSums(abs(B_est))>1e-9)] <- 1
    }
    list(paras_here=paras_here,SE=SE,var_sel=var_sel)
  }
  if(NCORES_w!=1){
    stopCluster(cl)
  }
  PARAS <- do.call(rbind,lapply(
    1:length(PRESS_il_y_and_i_l),
    function(ii){
      out <- cbind(PRESS_il_y_and_i_l[[ii]]$paras_here,ii)
      out <- cbind(out,1:nrow(out))
      out
    }))
  colnames(PARAS) <- c(colnames(paras),"id_list","id_in_list")
  ## Get RSS and PRESS
  RSS <- colSums(scale(y,scale = F)^2)
  paras_K_X_Y <- expand.grid(Ks,kYs)
  n_paras <- nrow(paras_K_X_Y)
  PRESS = Q2_h <- matrix(0,n_paras,q)
  PRESS_star <- rep(0,n_paras)
  Q2_sum = Q2_h_sum_star = Q2_sum_star <- rep(NA,n_paras)
  VAR_sel <- list()
  for(i in 1:n_paras){
    kX <- paras_K_X_Y[i,1]
    kY <- paras_K_X_Y[i,2]
    pos <- which(PARAS$keepX==kX & PARAS$keepY==kY)
    pos_list <- PARAS$id_list[pos]
    pos_in_list <- PARAS$id_in_list[pos]
    VAR_sel[[i]] <- list(kX=kX,kY=kY,select = matrix(0,length(pos),q),
                         PRESS=matrix(0,length(pos),q))
    for(j in 1:length(pos)){
      PRESS_j <- PRESS_il_y_and_i_l[[pos_list[j]]]$SE[pos_in_list[j],]
      var_sel_j <- which(PRESS_il_y_and_i_l[[pos_list[j]]]$var_sel[pos_in_list[j],]>1e-9)
      PRESS[i,] <- PRESS[i,] + PRESS_j
      # PRESS_star[i] <- PRESS_star[i] + sum(PRESS_j[var_sel_j])
      VAR_sel[[i]]$select[j,var_sel_j] <- 1
      VAR_sel[[i]]$PRESS[j,] <- PRESS_j
    }

    if(sparse){
      model <- mixOmics::spls(x,y,ncomp = 1,keepX = kX,keepY = kY)
    }else{
      model <- mixOmics::pls(x,y,ncomp = 1)
    }
    var_sel_model <- which(as.numeric(rowSums(abs(model$loadings$Y)))>1e-9)
    Q2_h[i,] <- 1-PRESS[i,]/RSS
    # Q2_reg[i] <- 1- sum((y_0-Y_pred_all)^2)/sum(y_0^2)
    Q2_sum[i] <- 1-sum(PRESS[i,])/sum(RSS)
    # Q2_h_sum_star[i] <- 1-sum(PRESS_star[i])/sum(RSS[var_sel_model])
  }
  # maxQ2_h_sum <- max(Q2_sum)#max(Q2_h_sum_star) # max(na.omit(c(Q2_h)))
  maxQ2_sum <- max(na.omit(Q2_sum))
  if(maxQ2_sum>tau){
    id_best <- which(Q2_sum==maxQ2_sum)#which(Q2_h_sum_star==maxQ2_h_sum,arr.ind = T)[1]#which(Q2_sum==maxQ2_sum)#
    kX_best <- paras_K_X_Y[id_best,1]
    kY_best <- paras_K_X_Y[id_best,2]
    model <- mixOmics::spls(x,y,ncomp = 1,keepX = kX_best,keepY = kY_best)
    prepre <- predict(model,newdata = x)
    y_est <- prepre$predict[,,1]
    y_res <- y-y_est
    t_r <- model$variates$X
    x_est <- t_r%*%crossprod(t_r,as.matrix(x))/sum(t_r^2)
    x_res <- x-x_est
    B_sol <- prepre$B.hat[,,1]
    variance_expl <- sum((scale(x,scale=F)%*%prepre$B.hat[,,1])^2)
    if(norm(B_sol)==0){# | model$explained_variance$Y<NZV){
      out <- NULL
    }else{
      out <- list(paras_K_X_Y=paras_K_X_Y,kX_best=kX_best,kY_best=kY_best,
                  B_sol=B_sol,x_res=x_res,y_res=y_res,RSS=RSS,PRESS=PRESS,
                  Q2_h=Q2_h,Q2_h_sum_star=Q2_h_sum_star,VAR_sel=VAR_sel,
                  variance_expl=variance_expl)
    }
  }else{
    out <- NULL
  }
  out
}

do_sPLS_all_comp <- function(x,y,Ks=NULL,kYs=NULL,B_th,tau=0.0975,NCORES=1,NZV = 0.01,sparse=T){
  n <- nrow(y);p <- ncol(x);q <- ncol(y)
  if(is.null(Ks)){
    Ks <- 1:p
  }
  if(is.null(kYs)){
    kYs <- 1:q
  }
  test <- TRUE
  h <- 0
  PRESS = RSS <- list()
  Q2_h <- 0
  Q2_reg <- NULL
  # B_all <- list()
  best_keepXs <- NULL
  y_res = y_0 <- scale(y)
  x_res = x_0 <- scale(x)
  Q2_cum_star_h = Q2_r_star_h <- NULL
  VAR_h_s <- NULL
  RSS0 <- colSums(y_res^2)
  while(test){
    ONE_COMP <- do_sPLS_one_comp(x_res,y_res,Ks=Ks,kYs=kYs,y_0 = y_0,
                                 NCORES=NCORES,tau=tau,sparse=sparse)
    if(is.null(ONE_COMP)){
      test <- F
      if(h==0){
        Q2_cum = Q2_reg <- 0
      }
    }else{
      variance_explained <- ONE_COMP$variance_expl/sum(RSS0)
      if(T){#variance_explained>NZV){
        if(h==0){
          Q2_h <- ONE_COMP$Q2_h
          best_keepXs <- ONE_COMP$kX_best
          best_keepYs <- ONE_COMP$kY_best
          id_kX <- which(ONE_COMP$paras_K_X_Y[,1]==best_keepXs)
          id_kY <- which(ONE_COMP$paras_K_X_Y[,2]==best_keepYs)
          id_all <- intersect(id_kX,id_kY)
          Q2_h_best <- ONE_COMP$Q2_h[id_all]
          # Q2_cum <- 1-sum(ONE_COMP$PRESS[id_all,])/sum(RSS0)#sum(ONE_COMP$RSS)
          # Q2_cum_star = Q2_cum_star_h = Q2_r_star_h <- ONE_COMP$Q2_h_sum_star[id_all]
        }else{
          id_kX <- which(ONE_COMP$paras_K_X_Y[,1]==ONE_COMP$kX_best)
          id_kY <- which(ONE_COMP$paras_K_X_Y[,2]==ONE_COMP$kY_best)
          id_all <- intersect(id_kX,id_kY)
          # Q2_r_star <- ONE_COMP$Q2_h_sum_star[id_all]
          Q2_h_best <- c(Q2_h_best,ONE_COMP$Q2_h[id_all])
          # Q2_cum_pot <- 1-(1-Q2_cum)*sum(ONE_COMP$PRESS[id_all,])/sum(RSS[[h]])#sum(ONE_COMP$RSS)
          # Q2_cum_star_pot <- 1-(1-Q2_cum_star)*(1-Q2_r_star)
          # if(Q2_cum_star_pot>Q2_cum_star){
          # Q2_cum_star <- 1-(1-Q2_cum_star)*(1-Q2_r_star)#Q2_cum_star_pot
          # Q2_cum <- 1-(1-Q2_cum)*sum(ONE_COMP$PRESS[id_all,])/sum(RSS[[h]])#sum(ONE_COMP$RSS)#Q2_cum_pot
          # Q2_cum_star_h <- c(Q2_cum_star_h,Q2_cum_star)
          # Q2_r_star_h <- c(Q2_r_star_h,Q2_r_star)
          Q2_h <- rbind(Q2_h,ONE_COMP$Q2_h)
          best_keepXs <- c(best_keepXs,ONE_COMP$kX_best)
          best_keepYs <- c(best_keepYs,ONE_COMP$kY_best)
          # }else{
          #   test <- F
          # }
        }
        if(test){
          PRESS[[h+1]] <- ONE_COMP$PRESS
          RSS[[h+1]] <- sum(ONE_COMP$y_res^2)#ONE_COMP$RSS
          # STAR
          kX_sol <- best_keepXs[length(best_keepXs)]
          kY_sol <- best_keepYs[length(best_keepYs)]
          id_var_sel <- do.call(rbind,lapply(ONE_COMP$VAR_sel,function(vv){c(vv$kX,vv$kY)}))
          id_all <- intersect(which(id_var_sel[,1]==kX_sol),which(id_var_sel[,2]==kY_sol))
          PRESS_star_sel <- ONE_COMP$VAR_sel[[id_all]]$PRESS
          Var_sel_loo <- ONE_COMP$VAR_sel[[id_all]]$select
          # PRESS_star_sel <- sum(Var_sel_loo*PRESS_star_sel)
          y_res_0 <- y_res
          y_res <- ONE_COMP$y_res
          x_res <- ONE_COMP$x_res
          Var_h <- 1-sum((y_res-y_0)^2)/sum(y_0^2)#(sum(scale(y_res_0,scale = F)^2)-sum(scale(y_res,scale = F)^2))/sum(RSS[[1]])
          #1-sum(scale(y_res,scale = F)^2)/sum(RSS[[1]])
          if(h==0){
            VAR_h_s <- Var_h
            # if(Var_h<NZV){
            #   test <- F
            # }
          }else{
            # if(Var_h<NZV){
            #   test <- F
            # }
            VAR_h_s <- c(VAR_h_s,Var_h)
          }
          h <- h + 1
        }
      }else{
        test <- F
      }
    }
  }
  if(h!=0){
    R_hat <- length(best_keepXs)
    # y_res <- scale(y)
    # x_res <- scale(x)
    if(sparse){
      model <- mixOmics::spls(x_0,y_0,ncomp = R_hat,keepX = best_keepXs,keepY = best_keepYs)
    }else{
      model <- mixOmics::pls(x_0,y_0,ncomp = R_hat)
    }
    # y_sc <- scale(y,scale = F)
    ERRORS_LOO <- rep(0,n)
    Y_pred_all <- matrix(0,n,q)
    var_sel <- matrix(0,n,q)
    for(ii in 1:n){
      X_train <- x[-ii,,drop=F]
      X_test <- x[ii,,drop=F]
      Y_train <- y[-ii,,drop=FALSE]
      Y_test <- y[ii,,drop=FALSE]
      mu_k_ii <- colMeans(X_train)
      mu_y_ii <- colMeans(Y_train)
      # sd_x_inv_ii <- unlist(lapply(apply(X_train,2,sd),function(ss){if(abs(ss)>1e-9){out <- 1/ss}else{out <- 0};out}))
      # sd_x_mat_ii <- apply(X_train,2,sd)
      # sd_y_mat_ii <- apply(Y,2,sd)
      # sd_y_inv_ii <- unlist(lapply(apply(Y,2,sd),function(ss){if(abs(ss)>1e-9){out <- 1/ss}else{out <- 0};out}))
      m_plus <- mixOmics::spls(X_train,Y_train,ncomp = R_hat,keepX = best_keepXs,keepY = best_keepYs)
      B_hat <- predict(m_plus,newdata = X_train)$B.hat[,,R_hat]
      var_sel[ii,which(colSums(abs(B_hat))>1e-9)] <- 1
      Y_pred_all[ii,] <- mu_y_ii + as.matrix(X_test-mu_k_ii)%*%B_hat
    }
    # err_LOO <- (y_res-Y_pred_all)
    # ERRORS_LOO <- colSums(err_LOO^2)
    # ERRORS_LOO_star <- colSums((err_LOO^2)*var_sel)
    # RSS0_star_model <- sum(colSums(y_res^2)[as.vector(which(abs(rowSums(model$loadings$Y))>1e-9))])
    # Q2_reg_star <- 1 - sum(ERRORS_LOO_star)/sum(RSS0_star_model)
    # Q2_reg_star <- 1-PRESS_star_sel/sum(y_res[,as.vector(which(abs(rowSums(model$loadings$Y))>1e-9))]^2)
    Q2_reg <- 1- sum((y-Y_pred_all)^2)/sum(y_0^2)
    # Q2_reg <- 1 - sum(PRESS[[h]][id_all,])/sum(y_res^2)
    B_hat <- predict(model,newdata = x_0)$B.hat[,,R_hat]
    stat_sel <- compare_selection(B_hat,B_th)
    dist_B <- sum((B_hat-B_th)^2)
    out <- list(R_hat=R_hat,B_hat=B_hat,SEL=stat_sel,
                explained_variance=VAR_h_s,
                Q2_cum=Q2_cum,Q2_reg=Q2_reg,Q2_h=Q2_h,
                Q2_h_best=Q2_h_best,
                best_keepXs=best_keepXs,best_keepYs=best_keepYs,model=model)
  }
}

do_mixomics <- function(Xs,Y,kxs,kys,ncores_i,R_max = 10){
  n <- nrow(Y)
  fold <- min(n,100)
  x <- data.frame(Xs[[1]])
  para_keep <- expand.grid(kxs,kys)
  idid <- rep((nrow(para_keep)-nrow(para_keep)%%ncores_i)/ncores_i,ncores_i)
  if(nrow(para_keep)%%ncores_i>0){
    idid[1:(nrow(para_keep)%%ncores_i)] <- idid[1:(nrow(para_keep)%%ncores_i)] + 1
  }
  id_core <- unlist(lapply(1:length(idid),function(ii){rep(ii,idid[ii])}))
  r_oo <- 1
  test_Q2 <- T
  kXi=kYi <- NULL
  Q2_final <- -20
  Q2_h <- NULL
  while(test_Q2){
    cat(r_oo);cat(" ");cat(Q2_final);cat(" ")
    cl <- makeCluster(ncores_i)
    registerDoParallel(cl)
    res_all <- foreach(ii=1:ncores_i,.packages = "mixOmics",.combine='rbind') %dopar% {
      id_here <- which(id_core==ii)
      out <- matrix(NA,length(id_here),4)
      for(id_ro in 1:length(id_here)){
        ro <- id_here[id_ro]
        kxi <- c(kXi,para_keep[ro,1]);kyi <- c(kYi,para_keep[ro,2])
        model.spls = spls(x, scale(Y), ncomp = r_oo, mode = 'regression',
                          keepX = kxi, keepY = kyi)
        if(n!=fold){
          model.spls.val <- perf(model.spls, validation = "Mfold", folds = fold)
        }else{
          model.spls.val <- perf(model.spls, validation = "loo")
        }
        Q2_tot <- model.spls.val$Q2.total[r_oo]
        Q2_true <- 1-sum(model.spls.val$PRESS[r_oo,])/sum(model.spls.val$RSS[1,])
        out[id_ro,] <- c(Q2_tot,para_keep[ro,1],para_keep[ro,2],Q2_true)
      }
      out
    }
    stopCluster(cl)
    id_ok <- which(res_all[,1]>0.0975)
    if(length(id_ok)>0){
      id_best  <- which.max(res_all[,1])
      kXi <- c(kXi,res_all[id_best,2])
      kYi <- c(kYi,res_all[id_best,3])
      r_oo <- r_oo+1
      Q2_final <- res_all[id_best,4]
      Q2_h <- c(Q2_h,res_all[id_best,1])
      if(r_oo == R_max+1){
        test_Q2 <- F
        r_oo <- r_oo-1
      }
    }else{
      test_Q2 <- F
      r_oo <- r_oo-1
    }
  }
  model.spls = mixOmics::spls(x, scale(Y), ncomp = r_oo, mode = 'regression',
                              keepX = kXi, keepY = kYi)
  Bb <- predict(model.spls,x)$B.hat[,,r_oo]
  list(optimal_parameters=list(Q2=Q2_final,Q2_h=Q2_h,R=r_oo,keepX=kXi,keepY=kYi),
       B_cbind=Bb,U=model.spls$loadings$X,V=model.spls$loadings$Y)
}

give_me_plot_select <- function(){
  par(mfrow=c(2,2))
  cols <- RColorBrewer::brewer.pal(length(method)+1,"Set1")[-6]
  col_box <- RColorBrewer::brewer.pal(9,"Pastel1")[7]

  ylim <- c(-0.25,1)
  boxplot(VP_X~method*n,df,
          main="True Positive (TP) for X selected variables (true=3)",
          border=cols[rep(1:l_m,i)],col=col_box,las=2,xlab="");abline(h=3,lty=2)
  abline(v=0.5+c(0:(length(ns)) )*length(unique(df$method)),lty=3)
  boxplot(SEL_X~method*n,df,main="Number of X selected variables (true=3)",
          border=cols[rep(1:l_m,i)],col=col_box,las=2,xlab="");abline(h=3,lty=2)
  abline(v=0.5+c(0:(length(ns)) )*length(unique(df$method)),lty=3)
  boxplot(VP_Y~method*n,df,main="True Positive (TP) for Y selected variables (true=2)",
          border=cols[rep(1:l_m,i)],col=col_box,las=2,xlab="");abline(h=2,lty=2)
  abline(v=0.5+c(0:(length(ns)) )*length(unique(df$method)),lty=3)
  boxplot(SEL_Y~method*n,df,main="Number of Y selected variables (true=2)",
          border=cols[rep(1:l_m,i)],col=col_box,las=2,xlab="");abline(h=2,lty=2)
  abline(v=0.5+c(0:(length(ns)) )*length(unique(df$method)),lty=3)
}

give_me_plot <- function(){
  par(mfrow=c(1,3),mar=c(3,3,3,2) )
  coco <- rep(1,1)
  # layout(matrix(c(coco,2*coco,6*coco,7*coco,
  #                 3*coco,4*coco,8*coco,9*coco,
  #                 5*coco,5*coco,10*coco,10*coco), 3, byrow = TRUE))
  layout(matrix(c(coco,2*coco,
                  3*coco,4*coco,
                  5*coco,6*coco), 3, byrow = TRUE))
  cols <- RColorBrewer::brewer.pal(length(method)+1,"Set1")[-6]
  col_box <- RColorBrewer::brewer.pal(9,"Pastel1")[7]

  ylim <- c(-0.25,1)

  boxplot(Q2_CUM_star~method*n,border=cols[rep(1:l_m,i)],col=col_box,df,main=expression("Q"[1:R]^"2,*"),
          xlab="",xaxt="n",ylim=ylim,ylab="")
  abline(h=0.0975)
  abline(v=0.5+c(0:(length(ns)) )*length(levels(df$method)),lty=3)
  axis(side = 1,at = length(levels(df$method))/2+0.5+c(0:(length(ns)-1) )*length(levels(df$method)),
       labels = paste("n=",ns,sep=""),cex.axis=0.9)
  legend("bottom",legend = level_legend,fill = cols,ncol = length(method),bg = "white",cex = 0.6)
  boxplot(Q2_CUM~method*n,col=col_box,border=cols[rep(1:l_m,i)],df,main=expression("Q"[1:R]^2),xlab="",xaxt="n",ylim=ylim,ylab="")
  abline(h=0.0975)
  abline(v=0.5+c(0:(length(ns)) )*length(levels(df$method)),lty=3)
  axis(side = 1,at = length(levels(df$method))/2+0.5+c(0:(length(ns)-1) )*length(levels(df$method)),
       labels = paste("n=",ns,sep=""),cex.axis=0.9)
  legend("bottom",legend = level_legend,fill = cols,ncol = length(method),bg = "white",cex = 0.6)

  boxplot(Q2_star~method*n,col=col_box,border=cols[rep(1:l_m,i)],df,main=expression("Q"^"2,*"),xlab="",xaxt="n",ylim=ylim,ylab="")
  abline(h=0.0975)
  abline(v=0.5+c(0:(length(ns)) )*length(levels(df$method)),lty=3)
  axis(side = 1,at = length(levels(df$method))/2+0.5+c(0:(length(ns)-1) )*length(levels(df$method)),
       labels = paste("n=",ns,sep=""),cex.axis=0.9)
  legend("bottom",legend = level_legend,fill = cols,ncol = length(method),bg = "white",cex = 0.6)
  boxplot(Q2~method*n,col=col_box,border=cols[rep(1:l_m,i)],df,main=expression("Q"^2),xlab="",xaxt="n",ylim=ylim,ylab="")
  abline(h=0.0975)
  abline(v=0.5+c(0:(length(ns)) )*length(levels(df$method)),lty=3)
  axis(side = 1,at = length(levels(df$method))/2+0.5+c(0:(length(ns)-1) )*length(levels(df$method)),
       labels = paste("n=",ns,sep=""),cex.axis=0.9)
  legend("bottom",legend = level_legend,fill = cols,ncol = length(method),bg = "white",cex = 0.6)

  boxplot(SE_B~method*n,df[-which(df$method=="OLS_esp"),],main=expression('||'~hat(B)-A^"+"*C~'||/||'~A^"+"*C~'||'),
          log="y",col=col_box,border=cols[rep(1:l_m,i)],xlab="",xaxt="n",ylab="");abline(h=0,lty=2)
  abline(v=0.5+c(0:(length(ns)) )*length(unique(df$method)),lty=3)
  axis(side = 1,at = length(levels(df$method))/2+0.5+c(0:(length(ns)-1) )*length(levels(df$method)),
       labels = paste("n=",ns,sep=""),cex.axis=0.9)
  legend("top",legend = level_legend,fill = cols,ncol = length(method),bg = "white",cex = 0.6)

  # boxplot(VP_X~method*n,df,
  #         main="True Positive (TP) for X selected variables (true=3)",
  #         col=cols[rep(1:l_m,i)],las=2,xlab="");abline(h=3,lty=2)
  # abline(v=0.5+c(0:(length(ns)) )*length(unique(df$method)),lty=3)
  # boxplot(SEL_X~method*n,df,main="Number of X selected variables (true=3)",
  #         col=cols[rep(1:l_m,i)],las=2,xlab="");abline(h=3,lty=2)
  # abline(v=0.5+c(0:(length(ns)) )*length(unique(df$method)),lty=3)
  # boxplot(VP_Y~method*n,df,main="True Positive (TP) for Y selected variables (true=2)",
  #         col=cols[rep(1:l_m,i)],las=2,xlab="");abline(h=2,lty=2)
  # abline(v=0.5+c(0:(length(ns)) )*length(unique(df$method)),lty=3)
  # boxplot(SEL_Y~method*n,df,main="Number of Y selected variables (true=2)",
  #         col=cols[rep(1:l_m,i)],las=2,xlab="");abline(h=2,lty=2)
  # abline(v=0.5+c(0:(length(ns)) )*length(unique(df$method)),lty=3)

  boxplot(R~method*n,df,main=expression(hat(R)),col=col_box,border=cols[rep(1:l_m,i)],xlab="",xaxt="n",ylab="")
  abline(h=length(which(svd(B_th_all)$d>1e-9)),lty=2)
  abline(v=0.5+c(0:(length(ns)-1) )*length(unique(df$method)),lty=3)
  axis(side = 1,at = length(levels(df$method))/2+0.5+c(0:(length(ns)-1) )*length(levels(df$method)),
       labels = paste("n=",ns,sep=""),cex.axis=0.9)
  legend("top",legend = level_legend,fill = cols,ncol = length(method),bg = "white",cex = 0.6)
}

give_me_plot_comm <- function(){
  par(mar=c(3,3,3,2))
  # layout(matrix(c(1,2,
  #                 1,2,
  #                 3,4,
  #                 3,5), 4, byrow = TRUE),mar=c(3,3,3,2))
  # layout(matrix(c(1,2), nrow = 1, byrow = TRUE))
  # layout(matrix(c(1,2,1,2,1,2,1,2,3,4,3,4,3,4,3,4,5,5),ncol=2,byrow = T))
  layout(matrix(c(rep(c(1,2),10),3,3),ncol=2,byrow = T))
  # par(mfrow=c(1,3),mar=c(3,3,3,2))
  cols <- RColorBrewer::brewer.pal(length(method)+1,"Set1")[-6]
  col_box <- RColorBrewer::brewer.pal(9,"Pastel1")[7]

  ylim <- c(-0.05,1)
  lwd <- 1.5
  ncols <- 3
  cex.leg <- 1
  b<-boxplot(Q2~method*n,border=cols,col=col_box,df,main=expression("Q"^"2"),#border=cols[rep(1:l_m,i)]
          xlab="",xaxt="n",ylab="",lwd=lwd,pch=1:length(level_legend));#abline(h=c(0,0.0975),lty=3,lwd=2)
  id_meth <- rep(1:nlevels(df$method),length(unique(df$n)))
  coli_s <- cols[id_meth]
  for(ii in 1:length(coli_s)){
    rect(xleft = ii-0.5+0.1,xright = ii+0.5-0.1,
         ytop = b$stats[4,ii],ybottom = b$stats[2,ii],
         density = density[id_meth[ii]],angle=angle[id_meth[ii]],col=coli_s[ii])
  }
  abline(v=0.5+c(0:(length(ns)) )*length(levels(df$method)),lty=3)
  abline(h=(2*eps^2)/3,lty=1,col="black")
  axis(side = 1,at = length(levels(df$method))/2+0.5+c(0:(length(ns)-1) )*length(levels(df$method)),
       labels = paste("n=",ns,sep=""),cex.axis=0.9)
  # legend("bottomright",legend = level_legend,fill = cols,ncol = 1,bg = "white",cex = cex.leg)
  text(x = 20,y = 0.85,labels = expression(gamma^"*"==frac(2,3)~.~epsilon^2==0.6534),pos = 3)
  axis(2,at = (2*eps^2)/3,labels = expression(gamma^"*"),las=1)

  # boxplot(time~method*n,border=cols,col=col_box,df,main="Raw time (not adjusted on number of parameters)",
  #         xlab="",xaxt="n",ylab="",lwd=lwd);#abline(h=c(0,0.0975),lty=3,lwd=2)
  # abline(v=0.5+c(0:(length(ns)) )*length(levels(df$method)),lty=3)
  # axis(side = 1,at = length(levels(df$method))/2+0.5+c(0:(length(ns)-1) )*length(levels(df$method)),
  #      labels = paste("n=",ns,sep=""),cex.axis=0.9)
  # legend("topright",legend = level_legend,fill = cols,ncol = 1,bg = "white",cex = cex.leg)

  b<-boxplot(log(1-DIST_B_hat)~method*n,df,#[-which(df$method=="OLS_esp"),],
          main=expression("||A"~hat(B)~"-D||"^2/"||D||"^2~" (logarithmic scale)"),
          col=col_box,border=cols,xlab="",yaxt="n",xaxt="n",ylab="",lwd=lwd,pch=1:length(level_legend))
  id_meth <- rep(1:nlevels(df$method),length(unique(df$n)))
  coli_s <- cols[id_meth]
  for(ii in 1:length(coli_s)){
    rect(xleft = ii-0.5+0.1,xright = ii+0.5-0.1,
         ytop = b$stats[4,ii],ybottom = b$stats[2,ii],
         density = density[id_meth[ii]],angle=angle[id_meth[ii]],col=coli_s[ii])
  }
  iii0 <- c((1e-1)^(rev(-1:4)))
  iii <- sort(c(iii0,5*iii0))
  abline(h=log(iii),lty=2,lwd=1/2,col="gray")
  axis(2,at=log(iii),labels=iii,las=2)
  abline(v=0.5+c(0:(length(ns)) )*length(unique(df$method)),lty=3)
  axis(side = 1,at = length(levels(df$method))/2+0.5+c(0:(length(ns)-1) )*length(levels(df$method)),
       labels = paste("n=",ns,sep=""),cex.axis=0.9)
  # legend("topright",legend = level_legend,fill = cols,ncol = 1,bg = "white",cex = cex.leg)

  # boxplot(NORM_B_hat~method*n,df,#[-which(df$method=="OLS_esp"),],
  #         main=expression("||"~hat(B)~"||"^2),
  #         col=col_box,border=cols,xlab="",xaxt="n",ylab="",lwd=lwd)#,ylim=c(0,1))
  # abline(h=0,lty=2)
  # abline(v=0.5+c(0:(length(ns)) )*length(unique(df$method)),lty=3)
  # axis(side = 1,at = length(levels(df$method))/2+0.5+c(0:(length(ns)-1) )*length(levels(df$method)),
  #      labels = paste("n=",ns,sep=""),cex.axis=0.9)
  # legend("topright",legend = level_legend,fill = cols,ncol = 1,bg = "white",cex = cex.leg)
  par(mar=rep(0,4))
  plot(0,type='n',axes=FALSE,ann=FALSE)
  legend("top",legend = level_legend,angle=angle,border=cols,fill=cols,density=density,col = cols,pch=1:length(level_legend),
         ncol = 6,bty = "n",cex = cex.leg,pt.cex=2)
}

plot_R <- function(){
  n_s <- sort(unique(df$n))
  n_n_s <- length(n_s)
  n_meth <- length(method)
  layout(matrix(1:(n_n_s*n_meth),nrow = n_n_s,byrow = T))
  cc <- 0.9
  par(mar=c(2,2,1,0),cex.main=cc,cex.lab=cc,cex.axis=cc)
  cols <- RColorBrewer::brewer.pal(length(method)+1,"Set1")[-6]
  breaks <- (min(c(2,na.omit(df$R)))-0.5):(max(c(2,na.omit(df$R)))+0.5)
  for(i_n in 1:n_n_s){
    for(i_mm in 1:n_meth){
      Ri <- df$R[which(df$method==method[i_mm] & df$n==ns[i_n])]
      # Ri[which(is.na(Ri))] <- 0

      hist(rep(2,100),breaks = breaks,xaxt="n",yaxt="n",col="gray",ylim=c(0,110),
           main="",xlab="",ylab="",density=20 , angle=33,border="gray80")
      axis(side = 1,unique(na.omit(df$R)),line = -1,labels = unique(na.omit(df$R)),tick = F)
      axis(side = 2,(0:10)*10,line = -1,labels = (0:10)*10,tick = F,las=2,gap.axis=1/4)
      abline(h=(0:5)*20,lty=2,col="gray")
      hist(Ri,breaks = breaks,col=cols[i_mm],add=T,border="gray20",
           density=density[i_mm] , angle=angle[i_mm])
      title(paste(level_legend[i_mm]," n=",ns[i_n],sep=""),line = 0)
      title(xlab=expression(hat(R)),ylab="Frequency",line = 1)
      xx <- unique(na.omit(df$R))
      yy <- unlist(lapply(xx,function(xxx){length(which(Ri==xxx))}))
      text(xx,yy,labels = yy,pos = 3,col=1,cex=cc*0.9)
    }
  }
}


plot_TVP_TFP_X <- function(){
  n_n_s <- length(unique(df$n))
  p <- ncol(A)
  var_star <- which(colSums(abs(A))>0)
  p_star <- length(var_star)
  sel_x <- df[,c(1,3,7,8)]
  TVP <- (sel_x$VP_X)/(p_star)
  TFP <- (sel_x$SEL_X- sel_x$VP_X)/(p-p_star)
  sel_x <- cbind(sel_x,TVP,TFP)

  layout(matrix(c(rep(c(1,2,3),15),rep(4,3)),ncol=3,byrow = T))
  par(mar=c(3.1, 3.1, 4.1, 0.2))
  cols <- RColorBrewer::brewer.pal(length(method)+1,"Set1")[-6]
  col_box <- RColorBrewer::brewer.pal(9,"Pastel1")[7]

  ylim <- c(-0.05,1.4*1)
  lwd <- 1.5
  ncols <- 3
  cex.leg <- 1

  b<-boxplot(SEL_X~method*n,border=cols,col=col_box,sel_x,main="Number of selected variables",#border=cols[rep(1:l_m,i)]
          xlab="",xaxt="n",ylab="",lwd=lwd,pch=1:length(level_legend))
  id_meth <- rep(1:nlevels(df$method),length(unique(df$n)))
  coli_s <- cols[id_meth]
  for(ii in 1:length(coli_s)){
    rect(xleft = ii-0.5+0.1,xright = ii+0.5-0.1,
         ytop = b$stats[4,ii],ybottom = b$stats[2,ii],
         density = density[id_meth[ii]],angle=angle[id_meth[ii]],col=coli_s[ii])
  }
  abline(h=100,lty=1,col=1)
  axis(side = 1,at = length(levels(sel_x$method))/2+0.5+c(0:(length(ns)-1) )*length(levels(sel_x$method)),
       labels = paste("n=",ns,sep=""),cex.axis=0.9)
  abline(v=0.5+c(0:(length(ns)) )*length(unique(sel_x$method)),lty=3)
  # legend("top",legend = level_legend,fill = cols,
  #        ncol = 2,bg = "white",cex = cex.leg,bty="n")

  b<-boxplot(TVP~method*n,border=cols,col=col_box,sel_x,main="True Positive Rate (TPR)",#border=cols[rep(1:l_m,i)]
          xlab="",xaxt="n",ylab="",lwd=lwd,pch=1:length(level_legend));#abline(h=c(0,0.0975),lty=3,lwd=2)
  id_meth <- rep(1:nlevels(df$method),length(unique(df$n)))
  coli_s <- cols[id_meth]
  for(ii in 1:length(coli_s)){
    rect(xleft = ii-0.5+0.1,xright = ii+0.5-0.1,
         ytop = b$stats[4,ii],ybottom = b$stats[2,ii],
         density = density[id_meth[ii]],angle=angle[id_meth[ii]],col=coli_s[ii])
  }
  abline(h=1,lty=1,col=1)
  axis(side = 1,at = length(levels(sel_x$method))/2+0.5+c(0:(length(ns)-1) )*length(levels(sel_x$method)),
       labels = paste("n=",ns,sep=""),cex.axis=0.9)
  abline(v=0.5+c(0:(length(ns)) )*length(unique(sel_x$method)),lty=3)
  # legend("top",legend = level_legend,fill = cols,
  #        ncol = 2,bg = "white",cex = cex.leg,bty="n")

  b<-boxplot(TFP~method*n,border=cols,col=col_box,sel_x,main="False Positive Rate (FPR)",#border=cols[rep(1:l_m,i)]
          xlab="",xaxt="n",ylab="",lwd=lwd,pch=1:length(level_legend));#abline(h=c(0,0.0975),lty=3,lwd=2)
  id_meth <- rep(1:nlevels(df$method),length(unique(df$n)))
  coli_s <- cols[id_meth]
  for(ii in 1:length(coli_s)){
    rect(xleft = ii-0.5+0.1,xright = ii+0.5-0.1,
         ytop = b$stats[4,ii],ybottom = b$stats[2,ii],
         density = density[id_meth[ii]],angle=angle[id_meth[ii]],col=coli_s[ii])
  }
  abline(h=0,lty=1,col=1)
  abline(v=0.5+c(0:(length(ns)) )*length(unique(sel_x$method)),lty=3)
  axis(side = 1,at = length(levels(sel_x$method))/2+0.5+c(0:(length(ns)-1) )*length(levels(sel_x$method)),
       labels = paste("n=",ns,sep=""),cex.axis=0.9)
  # legend("top",legend = level_legend,fill = cols,
  #        ncol = 2,bg = "white",cex = cex.leg,bty="n")
  par(mar=rep(0,4))
  plot(0,type='n',axes=FALSE,ann=FALSE)
  legend("top",legend = level_legend,angle=angle,border=cols,fill=cols,density=density,col = cols,pch=1:length(level_legend),
         ncol = 6,bty = "n",cex = cex.leg,pt.cex=2)
}

plot_TVP_TFP_Y <- function(){
  n_n_s <- length(unique(df$n))
  p <- ncol(D)
  var_star <- which(colSums(abs(D))>0)
  p_star <- length(var_star)
  sel_x <- df[,c(1,3,9,10)]
  TVP <- (sel_x$VP_Y)/(p_star)
  TFP <- (sel_x$SEL_Y- sel_x$VP_Y)/(p-p_star)
  sel_x <- cbind(sel_x,TVP,TFP)

  layout(matrix(c(rep(c(1,2,3),15),rep(4,3)),ncol=3,byrow = T))
  par(mar=c(3.1, 3.1, 4.1, 0.2))
  cols <- RColorBrewer::brewer.pal(length(method)+1,"Set1")[-6]
  col_box <- RColorBrewer::brewer.pal(9,"Pastel1")[7]

  ylim <- c(-0.05,1.4*1)
  lwd <- 1.5
  ncols <- 3
  cex.leg <- 1

  b<-boxplot(SEL_Y~method*n,border=cols,col=col_box,df,main="Number of selected variables",#border=cols[rep(1:l_m,i)]
          xlab="",xaxt="n",ylab="",lwd=lwd,pch=1:length(level_legend))
  id_meth <- rep(1:nlevels(df$method),length(unique(df$n)))
  coli_s <- cols[id_meth]
  for(ii in 1:length(coli_s)){
    rect(xleft = ii-0.5+0.1,xright = ii+0.5-0.1,
         ytop = b$stats[4,ii],ybottom = b$stats[2,ii],
         density = density[id_meth[ii]],angle=angle[id_meth[ii]],col=coli_s[ii])
  }
  abline(h=2,lty=1,col=1)
  axis(side = 1,at = length(levels(df$method))/2+0.5+c(0:(length(ns)-1) )*length(levels(df$method)),
       labels = paste("n=",ns,sep=""),cex.axis=0.9)
  abline(v=0.5+c(0:(length(ns)) )*length(unique(df$method)),lty=3)
  # legend("top",legend = level_legend,fill = cols,
  #        ncol = 2,bg = "white",cex = cex.leg,bty="n")


  b<-boxplot(TVP~method*n,border=cols,col=col_box,df,main="True Positive Rate (TPR)",#border=cols[rep(1:l_m,i)]
          xlab="",xaxt="n",ylab="",lwd=lwd,pch=1:length(level_legend));#abline(h=c(0,0.0975),lty=3,lwd=2)
  id_meth <- rep(1:nlevels(df$method),length(unique(df$n)))
  coli_s <- cols[id_meth]
  for(ii in 1:length(coli_s)){
    rect(xleft = ii-0.5+0.1,xright = ii+0.5-0.1,
         ytop = b$stats[4,ii],ybottom = b$stats[2,ii],
         density = density[id_meth[ii]],angle=angle[id_meth[ii]],col=coli_s[ii])
  }
  abline(h=1,lty=1,col=1)
  abline(v=0.5+c(0:(length(ns)) )*length(unique(df$method)),lty=3)
  axis(side = 1,at = length(levels(df$method))/2+0.5+c(0:(length(ns)-1) )*length(levels(df$method)),
       labels = paste("n=",ns,sep=""),cex.axis=0.9)
  # legend("top",legend = level_legend,fill = cols,
  #        ncol = 2,bg = "white",cex = cex.leg,bty="n")

  b<-boxplot(TFP~method*n,border=cols,col=col_box,df,main="False Positive Rate (FPR)",#border=cols[rep(1:l_m,i)]
          xlab="",xaxt="n",ylab="",lwd=lwd,pch=1:length(level_legend));#abline(h=c(0,0.0975),lty=3,lwd=2)
  id_meth <- rep(1:nlevels(df$method),length(unique(df$n)))
  coli_s <- cols[id_meth]
  for(ii in 1:length(coli_s)){
    rect(xleft = ii-0.5+0.1,xright = ii+0.5-0.1,
         ytop = b$stats[4,ii],ybottom = b$stats[2,ii],
         density = density[id_meth[ii]],angle=angle[id_meth[ii]],col=coli_s[ii])
  }
  abline(h=0,lty=1,col=1)
  abline(v=0.5+c(0:(length(ns)) )*length(unique(df$method)),lty=3)
  axis(side = 1,at = length(levels(df$method))/2+0.5+c(0:(length(ns)-1) )*length(levels(df$method)),
       labels = paste("n=",ns,sep=""),cex.axis=0.9)
  # legend("top",legend = level_legend,fill = cols,
  #        ncol = 2,bg = "white",cex = cex.leg,bty="n")

  par(mar=rep(0,4))
  plot(0,type='n',axes=FALSE,ann=FALSE)
  legend("top",legend = level_legend,angle=angle,border=cols,fill=cols,density=density,col = cols,pch=1:length(level_legend),
         ncol = 6,bty = "n",cex = cex.leg,pt.cex=2)
}

plot_sel_simu_x <- function(ncolX=8){
  CUMS <- 5
  n_n_s <- length(unique(df$n))
  lay_0 <- rep(1,n_n_s*CUMS)
  lay_1 <- do.call(cbind,lapply(1+1:n_n_s,function(ii){matrix(rep(ii,CUMS*6),ncol=CUMS)}))
  iioo <- c(n_n_s,n_n_s)+2
  lay_2 <- t(t(c(1,rev(c(iioo,iioo+1,iioo+2)))))
  layout(cbind(lay_2,rbind(lay_0,lay_1)))

  col_sel_x <- RColorBrewer::brewer.pal(ncolX,"Set2")
  col_sel_x[2:3] <- col_sel_x[c(3,2)]

  par(mar=c(0,0,0,0))
  plot(0,type='n',axes=FALSE,ann=FALSE)
  # plot(c(-1,1),c(-1,1),axes=FALSE,col="white")#,xlab="",ylab="")
  legend("center",ncol=8,legend = 1:8,fill=col_sel_x,cex=1.5,bty="n")

  sel_x <- df[,c(1,3,7)]
  nsiii <- sort(unique(sel_x$n))
  paras <- expand.grid(nsiii,unique(as.character(sel_x$method[which(sel_x$method %in% method[1:3])])))
  sel_x_good <- matrix(0,nrow(paras),max(na.omit(sel_x$SEL_X)))
  for(iii in 1:nrow(paras)){
    n_iii <- paras[iii,1];met <- as.character(paras[iii,2])
    popo <- which(df$n==n_iii & df$method==met)
    tata <- table(sel_x[popo,3])
    sel_x_good[iii,as.numeric(names(tata))] <- tata
  }
  for(i_n in 1:length(unique(df$n))){
    n_i <- nsiii[i_n]
    popo <- which(paras[,1]==n_i)
    par(mar=c(2,1,1,1))
    barplot(height = t(sel_x_good[popo,]),space=0.2,horiz=T,xlim=c(0,100),las=1,col = col_sel_x,names.arg = c("", "", ""),main=paste("n=",n_i,sep=""))
  }
  labs <- c("ddsPLS Boot", "sPLS [16] Boot", "SPLS [4] Boot")
  for(ii in 1:3){
    plot(c(-1,1),c(-1,1),axes=FALSE,col="white",xlab="",ylab="")
    text(cex=1.5,0,0, labs[ii], xpd=TRUE,srt=90)
  }
}

plot_sel_simu_y <- function(){
  CUMS <- 5
  n_n_s <- length(unique(df$n))
  lay_0 <- rep(1,n_n_s*CUMS)
  lay_1 <- do.call(cbind,lapply(1+1:n_n_s,function(ii){matrix(rep(ii,CUMS*6),ncol=CUMS)}))
  iioo <- c(n_n_s,n_n_s)+2
  lay_2 <- t(t(c(1,rev(c(iioo,iioo+1,iioo+2)))))
  layout(cbind(lay_2,rbind(lay_0,lay_1)))

  col_sel_y <- RColorBrewer::brewer.pal(3,"Set2")

  par(mar=c(0,0,0,0))
  plot(0,type='n',axes=FALSE,ann=FALSE)
  # plot(c(-1,1),c(-1,1),axes=FALSE,col="white")#,xlab="",ylab="")
  legend("center",ncol=3,legend = 1:3,fill=col_sel_y,cex=1.5,bty="n")

  sel_y <- df[,c(1,3,10)]
  nsiii <- sort(unique(sel_y$n))
  paras <- expand.grid(nsiii,unique(as.character(sel_y$method[which(sel_y$method %in% method[1:3])])))
  sel_y_good <- matrix(0,nrow(paras),max(as.numeric(names(table(sel_y$SEL_Y)))))
  for(iii in 1:nrow(paras)){
    n_iii <- paras[iii,1];met <- as.character(paras[iii,2])
    popo <- which(df$n==n_iii & df$method==met)
    tata <- table(sel_y[popo,3])
    sel_y_good[iii,as.numeric(names(tata))] <- tata
  }
  for(i_n in 1:length(unique(df$n))){
    n_i <- nsiii[i_n]
    popo <- which(paras[,1]==n_i)
    par(mar=c(2,1,1,1))
    barplot(height = t(sel_y_good[popo,]),space=0.2,horiz=T,xlim=c(0,100),las=1,
            col = col_sel_y,names.arg = c("", "", ""),main=paste("n=",n_i,sep=""))
  }
  labs <- c("ddsPLS Boot", "sPLS [16] Boot", "SPLS [4] Boot")
  for(ii in 1:3){
    plot(c(-1,1),c(-1,1),axes=FALSE,col="white",xlab="",ylab="")
    text(cex=1.5,0,0, labs[ii], xpd=TRUE,srt=90)
  }
}

test <- function(){
  if(T){
    library(doParallel)
    library(ddsPLS)

    NCORES <- 22
    NZV <- 1e-2

    eps1=eps2=eps3=epsY=eps <- 0.99#0.8#

    ff <- function(cc){out <- cc/sqrt(sum(cc^2));if(is.na(out[1])) out <- cc;out}
    A <- apply(cbind(c(1,1),c(0,1),c(0,0)),2,ff)
    d <- nrow(A)
    A <- eps1*cbind(A,matrix(0,nrow = d,ncol = 2))
    # A2 <- apply(cbind(c(2,1)),2,ff)
    # A2 <- eps1*cbind(A2,matrix(0,nrow = d,ncol = 1))
    # A3 <- eps3*matrix(0,nrow = d,ncol = 2)
    p1 = p2 <- 50; p <- 1000
    R1 <- 3; R2 <- 2; R3 <- 0;d <- R1+R2+R3
    A <- rbind(
      matrix(rep(c(rep(1,p1),rep(0,p2),rep(0,p-p1-p2)),R1),nrow = R1,byrow = T),
      matrix(rep(c(rep(0,p1),rep(2,p2),rep(0,p-p1-p2)),R2),nrow = R2,byrow = T)#,
      # matrix(rep(c(rep(0,p1),rep(1,p2),rep(0,p-p1-p2)),R3),nrow = R3,byrow = T)
    )
    A <- eps1*apply(A,2,ff)
    D <- rbind(
      matrix(rep(c(2,0,0),R1),nrow = R1,byrow = T),
      matrix(rep(c(0,1,0),R2),nrow = R2,byrow = T)#,
      # matrix(rep(c(0,1,0),R3),nrow = R3,byrow = T)
    )
    D <- eps1*apply(D,2,ff)
    # C <- epsY*apply(cbind(c(2,1),c(0,-1),rep(0,d)),2,ff)

    A_all <- A#cbind(A1,A2,A3)
    B_th_all=B_th <- MASS::ginv(A_all)%*%D
    L_total <- ncol(D)+ncol(A_all)

    LAMBDAS <- seq(0,1,length.out = 66)
    KXS <- unique(sort(c(round(seq(1,ncol(A_all),length.out = 22)),100) ))
    KYS <- unique(round(seq(1,ncol(D),length.out = 3)))

    ns <- c(25,50,100,200,400)#unique(round(seq(20,300,length.out = 8)))#unique(round(seq(20,150,length.out = 5)))
    NCORES_S <- rep(20,5)
    ALPHA <- rep(1/2,5)
    n_Bs <- c(1000,500,300,100,100)
    Ns <- 1:100
    paras <- expand.grid(ns,Ns)
    NNs <- nrow(paras)
    Q2_cum = error_B <- rep(NA,NNs)
    method <- c("ddsPLS Boot",
                "ddsPLS Unik Boot",
                "sPLS Boot",
                "PLS Boot",
                "SPLS Boot",
                "sPLS classik")
    level_legend <- c("ddsPLS","ddsPLS Unik Boot","sPLS [4] Boot",
                      "PLS Boot","SPLS [16] Boot","sPLS [4] classik")
    l_m <- length(method)
    names_df <- c("n","id","method","Q2","DIST_B_hat","NORM_B_hat",
                  "VP_X","SEL_X","VP_Y","SEL_Y","R","time")
    mat_0 <- expand.grid(ns,Ns,method)
    mat_plus <- matrix(NA,nrow(mat_0),length(names_df)-3)
    df <- data.frame(cbind(mat_0,mat_plus))
    names(df) <- names_df
    id_sel <- which(names(df) %in% c("VP_X","SEL_X","VP_Y","SEL_Y"))
    datas <- list(Xs=list(),Y=list(),phi=list())
    ALL_FUCKING_MODELS <- list()
    for(method_i in method){
      ALL_FUCKING_MODELS[[method_i]] <- list()
    }
    # load(../../Hadrien/data_last.RData")#load("../data_simu/data_signalFaible.RData")
    # i <- 6 ; i_m <- 1
    file_data <- "/Users/hlorenzo/Documents/GitHub/data_last_with_unik.RData"
    load(file_data)
    pos_old_unik <- which(df$method=="ddsPLS Unik Boot")
    if(length(pos_old_unik)>0){
      levels(df$method)[which(levels(df$method)=="ddsPLS Unik Boot")] <- "sPLS Boot Var-Q2"
      df$method[pos_old_unik] <- "sPLS Boot Var-Q2"
    }
    method[2] <- "sPLS Boot Var-Q2"
    level_legend[1] <- "ddsPLS"
    level_legend[2] <- "sPLS [7] 1"
    level_legend[3] <- "sPLS [7] 2"
    level_legend[4] <- "NIPALS-PLS"
    level_legend[5] <- "SPLS [1]"
    level_legend[6] <- "sPLS [7] 3"

    angle <- c(0,45,90,0,135,45,90)
    density <- c(50,50,50,25,25,50,25)
    # density=density ,angle=angle

  }
  # posINIT <- which(df$method==method[1] & df$R>2)
  posNA <- which(is.na(df$Q2))
  for(i in 41:500){#c(1:63)[posINIT]){#447:500){#386
    n <- paras[i,1]
    cat("\n\n______________________________________________________________________________")
    cat(paste("\n n =",n,"  ---   i =",i))
    pos <- intersect(which(df$n==n),which(df$id==paras[i,2]))
    pos_method <- unlist(lapply(method,function(mm){pos[which(df$method[pos]==mm)]}))
    # # Do data
    # psi <- MASS::mvrnorm(n,mu = rep(0,d+L_total),Sigma = diag(d+L_total))
    # phi <- psi[,1:d,drop=F];pt <- d
    # # SIs <- lapply(list(A1,A2,A3,C),function(M){do.call(cbind,lapply(sqrt(1-diag(crossprod(M))),function(sisi){rep(sisi,n)}))})
    # # X1 <- phi%*%A1 + SIs[[1]]*psi[,pt+1:ncol(A1),drop=F];pt <- pt + ncol(A1)
    # # X2 <- phi%*%A2 + SIs[[2]]*psi[,pt+1:ncol(A2),drop=F];pt <- pt + ncol(A2)
    # # X3 <- phi%*%A3 + SIs[[3]]*psi[,pt+1:ncol(A3),drop=F];pt <- pt + ncol(A3)
    # SIs <- lapply(list(A,D),function(M){do.call(cbind,lapply(sqrt(1-diag(crossprod(M))),function(sisi){rep(sisi,n)}))})
    # Xs <- list(x=phi%*%A + SIs[[1]]*psi[,pt+1:ncol(A),drop=F]);pt <- pt + ncol(A)
    # Y <- phi%*%D + SIs[[2]]*psi[,pt+1:ncol(D),drop=F];pt <- pt + ncol(Y)
    # datas$Xs[[i]] <- Xs
    # datas$Y[[i]] <- Y
    # datas$phi[[i]] <- phi
    # if(i%%250==0){
    #   save(datas,df,ALL_FUCKING_MODELS,file = file_data)
    #   #save(datas,df,file = "../data_simu/data_signalFaible.RData")
    #   # save(datas,df,varExplained,LAMBDAS_SOL,file = "../../Hadrien/data_signalFaible.RData")#save(datas,df,file = "../data_simu/data_signalFaible.RData")
    # }
    # Load data
    datas$phi[[i]] -> phi #; phi_ok <- phi[,1:2]
    datas$Xs[[i]] -> Xs
    datas$Y[[i]] -> Y
    x <- do.call(cbind,Xs)
    ##############
    for(i_m in 1:length(pos)){
      pos_i <- pos[i_m]
      method_i <- method[which(pos_method==pos_i)]
      cat(paste("\n     <- ",method_i,"... ",sep=""))
      toPlot <- pos_i %in% posNA#method_i %in% method[c(5,6)]
      if(toPlot){
        time_1 <- Sys.time()
        if(method_i %in% c("ddsPLS Boot","PLS Boot")){
          lambdas <- LAMBDAS
          pos_n <- which(ns==n)
          n_b_i <- n_Bs[pos_n]
          ncores_i <- NCORES_S[pos_n]
          if(method_i=="PLS Boot"){
            lambdas <- 0
            ncores_i <- 7
          }
          # res <- Q2_local_ddsPLS(Xs,Y,lambdas=lambdas,
          #                        n_B = n_b_i,NCORES=ncores_i,verbose = T)
          res <- sparse_PLS_Bootstrap(Xs,Y,type="CT",
                                      paras=lambdas,
                                      n_B=n_b_i,
                                      lowExplainedVariance=0,
                                      deflatX=T,NCORES=ncores_i,center=T,verbose=T)
          cat("\n")
        }else if(method_i %in% "ddsPLS Unik Boot"){
          lambdas <- LAMBDAS
          pos_n <- which(ns==n)
          n_b_i <- n_Bs[pos_n]
          ncores_i <- NCORES_S[pos_n]
          res <- Q2_UNIK_ddsPLS(Xs,Y,lambdas=lambdas,
                                n_B = n_b_i,NCORES=ncores_i,verbose = T)
        }else if(method_i %in% c("sPLS Boot","sPLS Boot Var-Q2") ){
          pos_n <- which(ns==n)
          n_b_i <- n_Bs[pos_n]
          ncores_i <- NCORES_S[pos_n]
          kxs <- KXS;kys <- KYS; sparse <- T
          # if(!(method_i %in% c("sPLS Boot")) ){
          #   kxs <- ncol(x);kys <- ncol(Y); sparse <- F;ncores_i <- 7
          # }
          # res <- Q2_boot_sPLS(Xs,Y,keepXs = kxs,keepYs=kys,
          #                     n_B=n_b_i,deflatX=T,NCORES=ncores_i,center=T,
          #                     NZV=1e-9,verbose=T,sparse = sparse)
          if(method_i == c("sPLS Boot")){
            whichCriterion <- "Q2"
          }else{
            whichCriterion <- "Q2-Var"
          }
          res <- Q2_boot_sPLS(Xs,Y,keepXs = kxs,keepYs=kys,whichCriterion = whichCriterion,
                              n_B=n_b_i,deflatX=T,NCORES=ncores_i,center=T,
                              NZV=1e-9,verbose=T,sparse = sparse)
        }else if(method_i %in% "SPLS Boot"){
          # etas <- seq(0.05,0.95,length.out = 5)
          # kappas <- seq(0.05,0.45,length.out = 5)
          # res = tryCatch({
          #   Q2_chun(x,Y,K=5,B_th=B_th,etas=etas,kappas=kappas,NCORES = length(kappas),NZV=NZV)
          # }, error = function(error_condition) {
          #   NULL
          # })
          etas <- seq(0.1,0.9,0.1)
          Ks <- 1:10
          fold <- min(n,100)
          cl <- makeCluster(length(etas))
          registerDoParallel(cl)
          res_all <- foreach(et=1:length(etas),.packages = "spls",.combine='rbind') %dopar% {
            cvmo <- cv.spls(Xs[[1]],Y,fold=fold,K = Ks,eta = etas[et],scale.y = T,plot.it = F)
            c(max(1-cvmo$mspemat),cvmo$K.opt,etas[et])
          }
          stopCluster(cl)
          id_max <- which.max(res_all[,1])
          mo <- spls::spls(Xs[[1]],Y,K = res_all[id_max,2],eta = res_all[id_max,3])
          res <- list(optimal_parameters=list(Q2=max(res_all[,1]),R=as.numeric(res_all[id_max,2])),B_cbind=mo$betahat)
        }else if(method_i %in% "sPLS classik"){
          ncores_i <- NCORES_S[pos_n]
          kxs <- KXS;kys <- KYS; sparse <- T
          res <- do_mixomics(Xs,Y,kxs,kys,ncores_i)
        }
        if(!is.null(res)){
          df[pos_i,]$Q2 <- res$optimal_parameters$Q2
          df[pos_i,]$DIST_B_hat <- 1-sum((A%*%res$B_cbind-D)^2)/sum(D^2)
          df[pos_i,]$NORM_B_hat <- sqrt(sum((res$B_cbind)^2))
          df[pos_i,]$R <- res$optimal_parameters$R
          if(method_i %in% c("sPLS classik",method[c(2,3)])){
            if(is.list(res$U)){
              df[pos_i,id_sel] <- compare_selection(tcrossprod(res$U[[1]],res$V),B_th_all)
            }else{
              df[pos_i,id_sel] <- compare_selection(tcrossprod(res$U,res$V),B_th_all)
            }
          }else{
            df[pos_i,id_sel] <- compare_selection(res$B_cbind,B_th_all)
          }
        }
        ALL_FUCKING_MODELS[[method_i]][[i]] <- res
        df$time[pos_i] <- as.numeric(difftime(Sys.time(), time_1, units = "secs"))
      }
      cat("\n\n______________________________________________________________________________")
      ##########################
      ########## PLOT ##########
      ##########################
      if(toPlot){
        # pdf(file = "/Users/hlorenzo/Documents/GitHub/Results_Last/Simulations_2.pdf",width = 15,height = 6)
        postscript("/Users/hlorenzo/Documents/GitHub/Results_Last/Simulations_2.eps",onefile=TRUE, horizontal=F,
                   width=9, height=2,pointsize=3)
        give_me_plot_comm()
        dev.off()

        # pdf(file = "/Users/hlorenzo/Documents/GitHub/Results_Last/Simulations_R_2.pdf",width = 15,height = 10)
        postscript("/Users/hlorenzo/Documents/GitHub/Results_Last/Simulations_R_2.eps", onefile=TRUE, horizontal=F,
                   width=11, height=5,pointsize=8)
        plot_R()
        dev.off()

        # pdf(file = "/Users/hlorenzo/Documents/GitHub/Results_Last/Simulations_sel_x_2.pdf",width = 15,height = 5)
        postscript("/Users/hlorenzo/Documents/GitHub/Results_Last/Simulations_sel_x_2.eps", onefile=TRUE, horizontal=F,
                   width=9, height=2,pointsize =2)
        # plot_X()
        plot_TVP_TFP_X()
        dev.off()

        # pdf(file = paste("/Users/hlorenzo/Documents/GitHub/Results_Last/Simulations_sel_y_2.pdf",sep=""),width = 15,height = 5)
        postscript("/Users/hlorenzo/Documents/GitHub/Results_Last/Simulations_sel_y_2.eps", onefile=TRUE, horizontal=F,
                   width=9, height=2,pointsize =2)
        # plot_sel_simu_y()
        plot_TVP_TFP_Y()
        dev.off()
      }
    }
  }
}
