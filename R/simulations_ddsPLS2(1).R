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

do_sPLS_one_comp <- function(x,y,Ks=NULL,kYs=NULL,NCORES=1,tau=0.0975){
  n <- nrow(y);p <- ncol(x);q <- ncol(y)
  kfolds <- n
  fold <- 1:n
  if(is.null(Ks)) Ks <- p
  if(is.null(kYs)) kYs <- q

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
      model <- mixOmics::spls(X_train,Y_train,ncomp = 1,keepX = kX,keepY = kY)
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
      PRESS_star[i] <- PRESS_star[i] + sum(PRESS_j[var_sel_j])
      VAR_sel[[i]]$select[j,var_sel_j] <- 1
      VAR_sel[[i]]$PRESS[j,] <- PRESS_j
    }
    model <- mixOmics::spls(x,y,ncomp = 1,keepX = kX,keepY = kY)
    var_sel_model <- which(as.numeric(rowSums(abs(model$loadings$Y)))>1e-9)
    Q2_h[i,] <- 1-PRESS[i,]/RSS
    Q2_sum[i] <- 1-sum(PRESS[i,])/sum(RSS)
    Q2_h_sum_star[i] <- 1-sum(PRESS_star[i])/sum(RSS[var_sel_model])
  }
  maxQ2_h_sum <- max(Q2_h_sum_star) # max(na.omit(c(Q2_h)))
  # maxQ2_sum <- max(na.omit(Q2_sum))
  if(maxQ2_h_sum>tau){
    id_best <- which(Q2_h_sum_star==maxQ2_h_sum,arr.ind = T)[1]#which(Q2_sum==maxQ2_sum)#
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

do_sPLS_all_comp <- function(x,y,Ks=NULL,kYs=NULL,B_th,tau=0.0975,NCORES=1,NZV = 0.01){
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
  y_res <- scale(y)
  x_res <- scale(x)
  Q2_cum_star_h = Q2_r_star_h <- NULL
  VAR_h_s <- NULL
  RSS0 <- colSums(y_res^2)
  while(test){
    ONE_COMP <- do_sPLS_one_comp(x_res,y_res,Ks=Ks,kYs=kYs,
                                 NCORES=NCORES,tau=tau)
    if(is.null(ONE_COMP)){
      test <- F
      if(h==0){
        Q2_cum = Q2_reg <- 0
      }
    }else{
      variance_explained <- ONE_COMP$variance_expl/sum(RSS0)
      if(variance_explained>NZV){
        if(h==0){
          Q2_h <- ONE_COMP$Q2_h
          best_keepXs <- ONE_COMP$kX_best
          best_keepYs <- ONE_COMP$kY_best
          id_kX <- which(ONE_COMP$paras_K_X_Y[,1]==best_keepXs)
          id_kY <- which(ONE_COMP$paras_K_X_Y[,2]==best_keepYs)
          id_all <- intersect(id_kX,id_kY)
          Q2_cum <- 1-sum(ONE_COMP$PRESS[id_all,])/sum(RSS0)#sum(ONE_COMP$RSS)
          Q2_cum_star = Q2_cum_star_h = Q2_r_star_h <- ONE_COMP$Q2_h_sum_star[id_all]
        }else{
          id_kX <- which(ONE_COMP$paras_K_X_Y[,1]==ONE_COMP$kX_best)
          id_kY <- which(ONE_COMP$paras_K_X_Y[,2]==ONE_COMP$kY_best)
          id_all <- intersect(id_kX,id_kY)
          Q2_r_star <- ONE_COMP$Q2_h_sum_star[id_all]
          # Q2_cum_pot <- 1-(1-Q2_cum)*sum(ONE_COMP$PRESS[id_all,])/sum(RSS[[h]])#sum(ONE_COMP$RSS)
          # Q2_cum_star_pot <- 1-(1-Q2_cum_star)*(1-Q2_r_star)
          # if(Q2_cum_star_pot>Q2_cum_star){
          Q2_cum_star <- 1-(1-Q2_cum_star)*(1-Q2_r_star)#Q2_cum_star_pot
          Q2_cum <- 1-(1-Q2_cum)*sum(ONE_COMP$PRESS[id_all,])/sum(RSS[[h]])#sum(ONE_COMP$RSS)#Q2_cum_pot
          Q2_cum_star_h <- c(Q2_cum_star_h,Q2_cum_star)
          Q2_r_star_h <- c(Q2_r_star_h,Q2_r_star)
          Q2_h <- rbind(Q2_h,ONE_COMP$Q2_h)
          best_keepXs <- c(best_keepXs,ONE_COMP$kX_best)
          best_keepYs <- c(best_keepYs,ONE_COMP$kY_best)
          # }else{
          #   test <- F
          # }
        }
        if(test){
          PRESS[[h+1]] <- ONE_COMP$PRESS
          RSS[[h+1]] <- as.numeric(colSums(ONE_COMP$y_res^2))#ONE_COMP$RSS
          # STAR
          kX_sol <- best_keepXs[length(best_keepXs)]
          kY_sol <- best_keepYs[length(best_keepYs)]
          id_var_sel <- do.call(rbind,lapply(ONE_COMP$VAR_sel,function(vv){c(vv$kX,vv$kY)}))
          id_all <- intersect(which(id_var_sel[,1]==kX_sol),which(id_var_sel[,2]==kY_sol))
          PRESS_star_sel <- ONE_COMP$VAR_sel[[id_all]]$PRESS
          Var_sel_loo <- ONE_COMP$VAR_sel[[id_all]]$select
          PRESS_star_sel <- sum(Var_sel_loo*PRESS_star_sel)
          y_res_0 <- y_res
          y_res <- ONE_COMP$y_res
          x_res <- ONE_COMP$x_res
          Var_h <- (sum(scale(y_res_0,scale = F)^2)-sum(scale(y_res,scale = F)^2))/sum(RSS[[1]])
          #1-sum(scale(y_res,scale = F)^2)/sum(RSS[[1]])
          if(h==0){
            VAR_h_s <- Var_h
            if(Var_h<NZV){
              test <- F
            }
          }else{
            if(Var_h<NZV){
              test <- F
            }
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
    y_res <- scale(y)
    x_res <- scale(x)
    model <- mixOmics::spls(x_res,y_res,ncomp = R_hat,keepX = best_keepXs,keepY = best_keepYs)
    # y_sc <- scale(y,scale = F)
    ERRORS_LOO <- rep(0,n)
    Y_pred_all <- matrix(0,n,q)
    var_sel <- matrix(0,n,q)
    for(ii in 1:n){
      X_train <- x_res[-ii,,drop=F]
      X_test <- x_res[ii,,drop=F]
      Y_train <- y_res[-ii,,drop=FALSE]
      Y_test <- y_res[ii,,drop=FALSE]
      mu_k_ii <- colMeans(X_train)
      mu_y_ii <- colMeans(Y_train)
      sd_x_inv_ii <- unlist(lapply(apply(X_train,2,sd),function(ss){if(abs(ss)>1e-9){out <- 1/ss}else{out <- 0};out}))
      # sd_x_mat_ii <- apply(X_train,2,sd)
      sd_y_mat_ii <- apply(Y,2,sd)
      # sd_y_inv_ii <- unlist(lapply(apply(Y,2,sd),function(ss){if(abs(ss)>1e-9){out <- 1/ss}else{out <- 0};out}))
      m_plus <- mixOmics::spls(x_res,y_res,ncomp = R_hat,keepX = best_keepXs,keepY = best_keepYs)
      B_hat <- predict(m_plus,newdata = X_train)$B.hat[,,R_hat]
      var_sel[ii,which(rowSums(abs(m_plus$loadings$Y))>1e-9)] <- 1
      Y_pred_all[ii,] <- mu_y_ii + ((X_test-mu_k_ii))%*%B_hat
    }
    err_LOO <- (y_res-Y_pred_all)
    ERRORS_LOO <- colSums(err_LOO^2)
    ERRORS_LOO_star <- colSums((err_LOO^2)*var_sel)
    RSS0_star_model <- sum(colSums(y_res^2)[as.vector(which(abs(rowSums(model$loadings$Y))>1e-9))])
    Q2_reg_star <- 1 - sum(ERRORS_LOO_star)/sum(RSS0_star_model)
    # Q2_reg_star <- 1-PRESS_star_sel/sum(y_res[,as.vector(which(abs(rowSums(model$loadings$Y))>1e-9))]^2)
    Q2_reg <- 1- sum(ERRORS_LOO)/sum(y_res^2)
    # Q2_reg <- 1 - sum(PRESS[[h]][id_all,])/sum(y_res^2)
    B_hat <- predict(model,newdata = x_res)$B.hat[,,R_hat]
    stat_sel <- compare_selection(B_hat,B_th)
    dist_B <- sum((B_hat-B_th)^2)
    out <- list(R_hat=R_hat,B_hat=B_hat,SEL=stat_sel,
                explained_variance=VAR_h_s,
                Q2_cum=Q2_cum,Q2_reg=Q2_reg,Q2_h=Q2_h,
                Q2_cum_star=Q2_cum_star,Q2_reg_star=Q2_reg_star,
                Q2_cum_star_h=Q2_cum_star_h,Q2_r_star_h=Q2_r_star_h,
                best_keepXs=best_keepXs,best_keepYs=best_keepYs,model=model)
  }
}

give_me_plot_select <- function(){
  par(mfrow=c(2,2))
  cols <- RColorBrewer::brewer.pal(length(method)+1,"Set1")[-6]
  col_box <- RColorBrewer::brewer.pal(9,"Pastel1")[7]

  level_legend <- c("ddsPLS","sPLS [16]","SPLS [4]", "PLS2 [26]",expression(X^"+"~Y),expression(A^"+"~C))

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

  level_legend <- c("ddsPLS","sPLS [16]","SPLS [4]", "PLS2 [26]",expression(X^"+"~Y),expression(A^"+"~C))

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
  # par(mar=c(3,3,3,2))
  # layout(matrix(c(1,2,
  #                 1,2,
  #                 3,4,
  #                 3,5), 4, byrow = TRUE))
  par(mfrow=c(2,2),mar=c(3,3,3,2))
  cols <- RColorBrewer::brewer.pal(length(method)+1,"Set1")[-6]
  col_box <- RColorBrewer::brewer.pal(9,"Pastel1")[7]

  level_legend <- c("ddsPLS","sPLS [16]","SPLS [4]", "PLS2 [26]",expression(X^"+"~Y),expression(A^"+"~C))

  ylim <- c(-0.25,0.9)
  lwd <- 1.5
  ncols <- 3
  cex.leg <- 1

  boxplot(Q2~method*n,border=cols[rep(1:l_m,i)],col=col_box,df,main=expression("Q"^"2"),
          xlab="",xaxt="n",ylim=ylim,ylab="",lwd=lwd);abline(h=0.0975,lty=3,lwd=2)
  abline(v=0.5+c(0:(length(ns)) )*length(levels(df$method)),lty=3)
  axis(side = 1,at = length(levels(df$method))/2+0.5+c(0:(length(ns)-1) )*length(levels(df$method)),
       labels = paste("n=",ns,sep=""),cex.axis=0.9)
  legend("bottomright",legend = level_legend,fill = cols,ncol = ncols,bg = "white",cex = cex.leg)

  boxplot(Q2_CUM~method*n,col=col_box,border=cols[rep(1:l_m,i)],df,main=expression("Q"[1:R]^2),
          xlab="",xaxt="n",ylim=ylim,ylab="",lwd=lwd);abline(h=0.0975,lty=3,lwd=2)
  abline(v=0.5+c(0:(length(ns)) )*length(levels(df$method)),lty=3)
  axis(side = 1,at = length(levels(df$method))/2+0.5+c(0:(length(ns)-1) )*length(levels(df$method)),
       labels = paste("n=",ns,sep=""),cex.axis=0.9)
  legend("bottomright",legend = level_legend,fill = cols,ncol = ncols,bg = "white",cex = cex.leg)

  boxplot(R~method*n,df,main=expression(hat(R)),col=col_box,border=cols[rep(1:l_m,i)],lwd=lwd,xlab="",xaxt="n",ylab="")
  abline(h=length(which(svd(B_th_all)$d>1e-9)),lty=3,lwd=2)
  abline(v=0.5+c(0:(length(ns)-1) )*length(unique(df$method)),lty=3)
  axis(side = 1,at = length(levels(df$method))/2+0.5+c(0:(length(ns)-1) )*length(levels(df$method)),
       labels = paste("n=",ns,sep=""),cex.axis=0.9)
  legend("bottomright",legend = level_legend,fill = cols,ncol = ncols,bg = "white",cex = cex.leg)

  boxplot(sqrt(SE_B)~method*n,df[-which(df$method=="OLS_esp"),],main=expression('||'~hat(B)-A^"+"*C~'||/||'~A^"+"*C~'||'),
          col=col_box,border=cols[rep(1:l_m,i)],xlab="",xaxt="n",ylab="",lwd=lwd,ylim=c(0,3.5));abline(h=0,lty=2)
  abline(v=0.5+c(0:(length(ns)) )*length(unique(df$method)),lty=3)
  axis(side = 1,at = length(levels(df$method))/2+0.5+c(0:(length(ns)-1) )*length(levels(df$method)),
       labels = paste("n=",ns,sep=""),cex.axis=0.9)
  legend("topright",legend = level_legend,fill = cols,ncol = ncols,bg = "white",cex = cex.leg)

  # boxplot(SEL_X~method*n,df,main="Number of X selected variables (true=3)",
  #         col=col_box,border=cols[rep(1:l_m,i)],las=2,xlab="",xaxt="n",ylab="",lwd=lwd);abline(h=3,lty=3,lwd=2)
  # axis(side = 1,at = length(levels(df$method))/2+0.5+c(0:(length(ns)-1) )*length(levels(df$method)),
  #      labels = paste("n=",ns,sep=""),cex.axis=0.9)
  # abline(v=0.5+c(0:(length(ns)) )*length(unique(df$method)),lty=3)
  # legend("bottomright",legend = level_legend,fill = cols,ncol = ncols,bg = "white",cex = cex.leg)
  #
  # boxplot(SEL_Y~method*n,df,main="Number of Y selected variables (true=2)",
  #         col=col_box,border=cols[rep(1:l_m,i)],las=2,xlab="",xaxt="n",ylab="",lwd=lwd);abline(h=2,lty=3,lwd=2)
  # axis(side = 1,at = length(levels(df$method))/2+0.5+c(0:(length(ns)-1) )*length(levels(df$method)),
  #      labels = paste("n=",ns,sep=""),cex.axis=0.9)
  # abline(v=0.5+c(0:(length(ns)) )*length(unique(df$method)),lty=3)
  # legend("right",legend = level_legend,fill = cols,ncol = ncols,bg = "white",cex = cex.leg)
}

plot_sel_simu_x <- function(){
  CUMS <- 5
  lay_0 <- rep(1,8*CUMS)
  lay_1 <- do.call(cbind,lapply(1+1:8,function(ii){matrix(rep(ii,CUMS*6),ncol=CUMS)}))
  lay_2 <- t(t(c(1,rev(c(10,10,11,11,12,12)))))
  layout(cbind(lay_2,rbind(lay_0,lay_1)))

  col_sel_x <- RColorBrewer::brewer.pal(8,"Set2")
  col_sel_x[2:3] <- col_sel_x[c(3,2)]

  par(mar=c(0,0,0,0))
  plot(0,type='n',axes=FALSE,ann=FALSE)
  # plot(c(-1,1),c(-1,1),axes=FALSE,col="white")#,xlab="",ylab="")
  legend("center",ncol=8,legend = 1:8,fill=col_sel_x,cex=1.5,bty="n")

  sel_x <- df[,c(1,3,8)]
  nsiii <- sort(unique(sel_x$n))
  paras <- expand.grid(nsiii,unique(as.character(sel_x$method[which(sel_x$method %in% method[1:3])])))
  sel_x_good <- matrix(0,nrow(paras),length(table(sel_x$SEL_X)))
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
    barplot(height = t(sel_x_good[popo,]),space=0.2,horiz=T,xlim=c(0,100),las=1,
            col = col_sel_x,names.arg = c("", "", ""),main=paste("n=",n_i,sep=""))
  }
  labs <- c("ddsPLS", "sPLS [16]", "SPLS [4]")
  for(ii in 1:3){
    plot(c(-1,1),c(-1,1),axes=FALSE,col="white",xlab="",ylab="")
    text(cex=1.5,0,0, labs[ii], xpd=TRUE,srt=90)
  }
}

plot_sel_simu_y <- function(){
  CUMS <- 5
  lay_0 <- rep(1,8*CUMS)
  lay_1 <- do.call(cbind,lapply(1+1:8,function(ii){matrix(rep(ii,CUMS*6),ncol=CUMS)}))
  lay_2 <- t(t(c(1,rev(c(10,10,11,11,12,12)))))
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
  labs <- c("ddsPLS", "sPLS [16]", "SPLS [4]")
  for(ii in 1:3){
    plot(c(-1,1),c(-1,1),axes=FALSE,col="white",xlab="",ylab="")
    text(cex=1.5,0,0, labs[ii], xpd=TRUE,srt=90)
  }
}

test <- function(){
  if(T){
    library(doParallel)

    NCORES <- 15
    NZV <- 1e-2

    a <- 1/sqrt(2)
    a2 <- 2/sqrt(5)
    b <- 1/sqrt(5)
    A0 <- rbind(t(matrix(c(a,0,0,0,a,1,0,0,0,0,1,1,0,0,0,0),nrow = 4)),c(0,0,0,0))
    A20 <- rbind(t(matrix(c(a2,0,0,0,sqrt(1-a2^2),1,0,0,0,0,1,-a,0,0,0,a),nrow = 4))[,c(1,4)],c(0,0))
    A30 <- cbind(c(0,0,1,0,0),c(0,0,2*b,b,0))
    C0 <- matrix(c(2*b,0,0,b,-1,0,0,0,0,0,0,0,0,0,1),nrow = 5,byrow = T)

    ## Add noise
    # eps <- 0.99;eps2 <- 0.99;eps3 <- 0.99;epsY <- 0.99#eps <- 0.95;eps2 <- 0.9;eps3 <- 0.99;epsY <- 0.9
    eps=eps2=eps3=epsY <- 0.99#0.8#
    nbVars <- c(ncol(A0),ncol(A20),ncol(A30),ncol(C0))
    totalVar <- sum(nbVars)
    I_total <- diag(rep(1,totalVar))
    A <- rbind(A0*eps,I_total[,1:nbVars[1]]*sqrt(1-eps^2))
    A2 <- rbind(A20*eps2,I_total[,nbVars[1]+1:(nbVars[2])]*sqrt(1-eps2^2))
    A3 <- rbind(A30*eps3,I_total[,sum(nbVars[1:2])+1:(nbVars[3])]*sqrt(1-eps3^2))
    C <- rbind(C0*epsY,I_total[,sum(nbVars[1:3])+1:(nbVars[4])]*sqrt(1-epsY^2))


    d <- nrow(C)
    A_all <- cbind(A,A2,A3)
    B_th_all=B_th <- MASS::ginv(A_all)%*%C

    lambdas <- seq(0,1,length.out = 100)

    ns <- unique(round(seq(20,300,length.out = 8)))#unique(round(seq(20,150,length.out = 5)))
    Ns <- 1:100
    paras <- expand.grid(ns,Ns)
    NNs <- nrow(paras)
    Q2_cum = error_B <- rep(NA,NNs)
    method <- c("ddsPLS_local",
                "sPLS (Le Cao)",
                "SPLS (Chun and Keles)",
                "PLS",
                "OLS",
                "OLS_esp")
    l_m <- length(method)
    names_df <- c("n","id","method","Q2","Q2_CUM",
                  "SE_B","VP_X","SEL_X","VP_Y","SEL_Y","R","Q2_star","Q2_CUM_star")
    mat_0 <- expand.grid(ns,Ns,method)
    mat_plus <- matrix(NA,nrow(mat_0),length(names_df)-3)
    df <- data.frame(cbind(mat_0,mat_plus))
    names(df) <- names_df
    id_sel <- which(names(df) %in% c("VP_X","SEL_X","VP_Y","SEL_Y"))
    datas <- list(Xs=list(),Y=list(),phi=list())
    LAMBDAS_SOL <- list()
    varExplained <- list()
    load("../../Hadrien/data_signalFort_no_1_2.RData")#load("../data_simu/data_signalFaible.RData")
    i <- 6 ; i_m <- 1
  }
  posINIT <- unique(which(df$method==method[1] & df$n %in% c(20,220) ))
  for(i in posINIT[which(posINIT>=305)]){#1:NNs){#386
    n <- paras[i,1]
    pos <- intersect(which(df$n==n),which(df$id==paras[i,2]))
    pos_method <- unlist(lapply(method,function(mm){pos[which(df$method[pos]==mm)]}))
    LAMBDAS_SOL[[i]]=varExplained[[i]] <- list()
    # # Do data
    # phi <- matrix(rnorm(n*d),nrow = n)
    # X <- phi%*%A;Y <- phi%*%C;X2 <- phi%*%A2;X3 <- phi%*%A3
    # Xs <- list(X,X2,X3)
    # datas$Xs[[i]] <- Xs
    # datas$Y[[i]] <- Y
    # datas$phi[[i]] <- phi
    if(i%%5==0){
      save(datas,df,varExplained,LAMBDAS_SOL,file = "../../Hadrien/data_signalFort_no_1_2.RData")#save(datas,df,file = "../data_simu/data_signalFaible.RData")
      # save(datas,df,varExplained,LAMBDAS_SOL,file = "../../Hadrien/data_signalFaible.RData")#save(datas,df,file = "../data_simu/data_signalFaible.RData")
    }
    # Load data
    datas$Xs[[i]] -> Xs
    datas$Y[[i]] -> Y
    x <- do.call(cbind,Xs)
    sel_no_sp <- c(length(which(rowSums(abs(B_th))>1e-9)),ncol(x),length(which(colSums(abs(B_th))>1e-9)),ncol(Y))
    sensib_no_sp_X <- sel_no_sp[1]/(sel_no_sp[1]+0)
    specif_no_sp_X <- 0/(0+sel_no_sp[2]-sel_no_sp[1])
    sensib_no_sp_Y <- sel_no_sp[3]/(sel_no_sp[3]+0)
    specif_no_sp_Y <- 0/(0+sel_no_sp[4]-sel_no_sp[3])
    ##############
    cat("\n   ");cat(i);cat("   ")
    for(i_m in 1:length(pos)){
      pos_i <- pos[i_m]
      method_i <- method[which(pos_method==pos_i)]
      cat(paste(method_i,"... ",sep=""))
      toPlot <- method_i %in% method[1]
      if(toPlot){#T){#
        if(method_i %in% c("ddsPLS_local","ddsPLS_global","PLS")){
          lambda_max <- 1
          N_lambdas <- 100
          if(method_i=="PLS"){
            lambda_max <- 0
            N_lambdas <- 1
            df[pos_i,id_sel] <- sel_no_sp
          }
          if(n==20) n_B <- 1000
          if(n==220) n_B <- 100
          res <- Q2_local_ddsPLS(Xs,Y,N_lambdas = N_lambdas,lambda_max = lambda_max,
                                 n_B = n_B,NCORES=20,verbose = T)
          if(F){
            R_hat <- res$optimal_parameters$R
            uu <- do.call(rbind,res$Us)
            ts <- scale(x)%*%uu
            tau_r <- apply(ts,2,sd)
            W_star <- uu%*%diag(1/tau_r)
            A_hat <- MASS::ginv(W_star)
            AA_th <- cbind(A,A2,A3)[1:2,]
            cor(AA_th,A_hat)
          }
          if(!is.null(res)){
            df[pos_i,]$Q2 <- res$optimal_parameters$Q2
            # df[pos_i,]$Q2_star <- res$optimal_parameters$Q2_reg_star
            df[pos_i,]$SE_B <- sum((res$B_cbind-B_th_all)^2)
            # df[pos_i,]$Q2_CUM <- res$optimal_parameters$Q2_cum
            # df[pos_i,]$Q2_CUM_star <- res$optimal_parameters$Q2_cum_star
            df[pos_i,]$R <- res$optimal_parameters$R
            df[pos_i,id_sel] <- compare_selection(res$B_cbind,B_th_all)
            if(method_i=="PLS"){
              df[pos_i,]$Q2_star <- df[pos_i,]$Q2
            }
            LAMBDAS_SOL[[i]][[method_i]] <- res$optimal_parameters$lambda
            varExplained[[i]][[method_i]] <- res$explained_variance
            # cat(" |-");
            # cat(i);
            # cat(", n=");cat(df[pos_i,]$n)
            # cat(", R=");cat(df[pos_i,]$R)
            # cat(", Q2_cum=");cat(round(df[pos_i,]$Q2_CUM,2))
            # cat(", Q2_reg=");cat(round(df[pos_i,]$Q2,2))
            # cat(", Q2_cum_star=");cat(round(df[pos_i,]$Q2_CUM_star,2))
            # cat(", Q2_reg_star=");cat(round(df[pos_i,]$Q2_star,2));cat("-| \n")
          }else if(method_i=="sPLS (Le Cao)"){
            res <- tryCatch({
              do_sPLS_all_comp(as.data.frame(x),Y,B_th=B_th,Ks=1:8,kYs=1:3,NCORES = 10,NZV = NZV)
            }, error = function(error_condition) {
              NULL
            })
            # if(!is.null(res)){
            df$Q2_CUM[pos_i] <- res$Q2_cum
            df$Q2[pos_i] <- res$Q2_reg
            df$Q2_CUM_star[pos_i] <- res$Q2_cum_star
            df$Q2_star[pos_i] <- res$Q2_reg_star
            df$SE_B[pos_i] <- sum((res$B_hat-B_th)^2)
            df[pos_i,id_sel] <- res$SEL
            df$R[pos_i] <- res$R_hat
            LAMBDAS_SOL[[i]][[method_i]] <- list(keepX = res$best_keepXs,keepY = res$best_keepYs)
            varExplained[[i]][[method_i]] <- res$VAR_h_s
          }
        }else if(method_i=="SPLS (Chun and Keles)"){
          etas <- seq(0.05,0.95,length.out = 5)
          kappas <- seq(0.05,0.45,length.out = 5)
          res = tryCatch({
            Q2_chun(x,Y,K=5,B_th=B_th,etas=etas,kappas=kappas,NCORES = length(kappas),NZV=NZV)
          }, error = function(error_condition) {
            NULL
          })
          if(!is.null(res)){
            df$Q2_CUM[pos_i] <- res$Q2_cum
            df$Q2[pos_i] <- res$Q2
            df$Q2_CUM_star[pos_i] <- res$Q2_cum_star
            df$Q2_star[pos_i] <- res$Q2_star
            df$SE_B[pos_i] <- res$SE_B
            df[pos_i,id_sel] <- res[,c("VP_X","SEL_X","VP_Y","SEL_Y")]
            df$R[pos_i] <- res$R
            LAMBDAS_SOL[[i]][[method_i]] <- list(eta = res$eta,R=res$R,kappa = res$kappa)
          }
        }else if(method_i=="OLS"){
          df[pos_i,]$Q2 <- Q2_OLS(x,Y)
          df$Q2_star[pos_i] <-  df[pos_i,]$Q2
          df[pos_i,]$SE_B <- sum((MASS::ginv(scale(do.call(cbind,Xs)))%*%scale(Y)-B_th_all)^2)
          df[pos_i,id_sel] <- sel_no_sp
        }else if(method_i=="OLS_esp"){
          df[pos_i,]$Q2 <- Q2_OLS_th(x,Y,B_th_all)
          df$Q2_star[pos_i] <-  df[pos_i,]$Q2
          df[pos_i,]$SE_B <- 0
          df[pos_i,id_sel] <- sel_no_sp
        }
      }
    }
    # n_ok <- df$n[1:i]
    # matou <- do.call(rbind,lapply(LAMBDAS_SOL[1:i],function(o){a<-o$ddsPLS_local;c(a,rep(NA,20-length(a)))}))
    # comp_ok <- which(unlist(lapply(apply(matou,2,table),length))!=0)
    # matou_ok <- matou[,comp_ok]
    # par(mfrow=c(length(ns),length(comp_ok)))
    # for(i_n in 1:length(ns)){
    #   n_i <- ns[i_n]
    #   popo <- which(n_ok==n_i)
    #   for(r in comp_ok){
    #     vals <- matou_ok[popo,r]
    #     if(length(na.omit(vals))>0){
    #       hist(vals,40,main=paste(n_i,r),xlim=c(0,1))
    #     }else{
    #       plot(0,0,col="white")
    #     }
    #   }
    # }
    ##########################
    ########## PLOT ##########
    ##########################
    if(T){
      # postscript("/Users/hlorenzo/Dropbox/Results/simulations_plotsFaible.eps", width=16, height=10, onefile=TRUE, horizontal=FALSE)
      # give_me_plot()
      # dev.off()
      #
      # pdf(file = paste("/Users/hlorenzo/Dropbox/Results/simu_all_DoneFaible.pdf",sep=""),width = 14,height = 12)
      # give_me_plot()
      # dev.off()
      #
      # pdf(file = paste("/Users/hlorenzo/Dropbox/Results/simu_all_Done_selectFaible.pdf",sep=""),width = 14,height = 12)
      # give_me_plot_select()
      # dev.off()



      # postscript("/Users/hlorenzo/Dropbox/Results/simulations_plots.eps", width=16, height=10, onefile=TRUE, horizontal=FALSE)
      # give_me_plot()
      # dev.off()
      # pdf(file = paste("/Users/hlorenzo/Dropbox/Results/simu_all_Done.pdf",sep=""),width = 14,height = 12)
      # give_me_plot()
      # dev.off()
      # pdf(file = paste("/Users/hlorenzo/Dropbox/Results/simu_all_Done_select.pdf",sep=""),width = 14,height = 12)
      # give_me_plot_select()
      # dev.off()

      # pdf(file = "/Users/hlorenzo/Dropbox/Results/Simulations.pdf",width = 15,height = 8)
      postscript("/Users/hlorenzo/Dropbox/Results/Simulations.eps", width=30, height=10, onefile=TRUE, horizontal=T)
      give_me_plot_comm()
      dev.off()

      pdf(file = "/Users/hlorenzo/Dropbox/Results/Simulations_sel_x.pdf",width = 14,height = 9)
      # postscript("/Users/hlorenzo/Dropbox/Results/Simulations_sel_x.eps", width=14, height=9, onefile=TRUE, horizontal=FALSE)
      plot_sel_simu_x()
      dev.off()

      pdf(file = paste("/Users/hlorenzo/Dropbox/Results/Simulations_sel_y.pdf",sep=""),width = 14,height = 9)
      # postscript("/Users/hlorenzo/Dropbox/Results/Simulations_sel_y.eps", width=14, height=9, onefile=TRUE, horizontal=F)
      plot_sel_simu_y()
      dev.off()
    }
    ##########################
    ########## STAT ##########
    ##########################
    # sel_th_X <- sel_no_sp[1]
    # ALL_X <- sel_no_sp[2]
    # sensib_X <- df$VP_X/sel_th_X
    # FN <- sel_th_X-df$VP_X
    # VN <- (ALL_X-df$SEL_X)-FN
    # specif_X <- VN/(VN+df$SEL_X-df$VP_X)
    #
    # sel_th_Y <- sel_no_sp[3]
    # ALL_Y <- sel_no_sp[4]
    # sensib_Y <- df$VP_Y/sel_th_Y
    # FN <- sel_th_Y-df$VP_Y
    # VN <- (ALL_Y-df$SEL_Y)-FN
    # specif_Y <- VN/(VN+df$SEL_Y-df$VP_Y)
    #
    # df_sen_spe <- data.frame(list(n=df$n,method=df$method,sensibility=sensib_X,specificity=specif_X))
    # boxplot(sensibility~method*n,df_sen_spe,main="Number of Y selected variables (true=2)",
    #         col=cols[rep(1:l_m,i)],las=2,xlab="");abline(h=2,lty=2)
    # abline(v=0.5+c(0:(length(ns)) )*length(unique(df_sen_spe$method)),lty=3)
    # boxplot(specificity~method*n,df_sen_spe,main="Number of Y selected variables (true=2)",
    #         col=cols[rep(1:l_m,i)],las=2,xlab="");abline(h=2,lty=2)
    # abline(v=0.5+c(0:(length(ns)) )*length(unique(df_sen_spe$method)),lty=3)

    # paras_stat <- expand.grid(ns,method)
    # n_stat <- nrow(paras_stat)
    # out_stat <- data.frame(cbind(paras_stat,matrix(NA,n_stat,ncol(df)-3-3)))
    # names(out_stat) <- names(df)[-c(2,which(names(df) %in% paste("V",14:16,sep="")))]
    # out_stat_sens_spec <- data.frame(cbind(paras_stat,matrix(NA,n_stat,4)))
    # nana <- expand.grid(c("Sensibility","Specificity"),c("X","Y"))
    # names(out_stat_sens_spec) <- c(names(df)[1:2],paste(nana[,1],nana[,2]))
    # for(i_stat in 1:n_stat){
    #   m_stat <-   paras_stat[i_stat,2]
    #   n_stat <-   paras_stat[i_stat,1]
    #   po <- which(df$n==n_stat & df$method==m_stat)
    #   for(j in 3:ncol(out_stat)){
    #     d_i_s <- na.omit(df[po,j+1])
    #     if(any(is.infinite(d_i_s ))){
    #       d_i_s <- d_i_s[-which(is.infinite(d_i_s ))]
    #     }
    #     if(length(d_i_s)>2){
    #       if(names(out_stat)[j]=="SE_B"){
    #         # norm_th <- sqrt(sum(B_th^2))
    #         # stat_SE <- sqrt(d_i_s)/norm_th
    #         out_stat[i_stat,j] <- paste("$",round(mean(stat_SE,na.rm = T),2)," \\pm  ",
    #                                     round(sd(stat_SE,na.rm = T),2),"$",sep="")
    #
    #       }else{
    #         out_stat[i_stat,j] <- paste("$",round(mean(d_i_s,na.rm = T),2)," \\pm  ",
    #                                     round(sd(d_i_s,na.rm = T),2),"$",sep="")
    #       }
    #     }else{
    #       out_stat[i_stat,j] <- ""
    #     }
    #   }
    #   if(F){
    #     ## SENSIBILITY
    #     ##### X
    #     d_i_s_sens <- na.omit(sensib_X[po])
    #     if(any(is.infinite(d_i_s_sens ))){
    #       d_i_s_sens <- d_i_s_sens[-which(is.infinite(d_i_s_sens ))]
    #     }
    #     if(length(d_i_s_sens)>2){
    #       out_stat_sens_spec$`Sensibility X`[i_stat] <- paste("$",round(mean(d_i_s_sens,na.rm = T),2)," \\pm  ",
    #                                                           round(sd(d_i_s_sens,na.rm = T),2),"$",sep="")
    #     }else{
    #       out_stat_sens_spec$`Sensibility X`[i_stat] <- ""
    #     }
    #     ##### Y
    #     d_i_s_sens <- na.omit(sensib_Y[po])
    #     if(any(is.infinite(d_i_s_sens ))){
    #       d_i_s_sens <- d_i_s_sens[-which(is.infinite(d_i_s_sens ))]
    #     }
    #     if(length(d_i_s_sens)>2){
    #       out_stat_sens_spec$`Sensibility Y`[i_stat] <- paste("$",round(mean(d_i_s_sens,na.rm = T),2)," \\pm  ",
    #                                                           round(sd(d_i_s_sens,na.rm = T),2),"$",sep="")
    #     }else{
    #       out_stat_sens_spec$`Sensibility Y`[i_stat] <- ""
    #     }
    #     ## SPECIFICITY
    #     ##### X
    #     d_i_s_sens <- na.omit(specif_X[po])
    #     if(any(is.infinite(d_i_s_sens ))){
    #       d_i_s_sens <- d_i_s_sens[-which(is.infinite(d_i_s_sens ))]
    #     }
    #     if(length(d_i_s_sens)>2){
    #       out_stat_sens_spec$`Specificity X`[i_stat] <- paste("$",round(mean(d_i_s_sens,na.rm = T),2)," \\pm  ",
    #                                                           round(sd(d_i_s_sens,na.rm = T),2),"$",sep="")
    #     }else{
    #       out_stat_sens_spec$`Specificity X`[i_stat] <- ""
    #     }
    #     ##### Y
    #     d_i_s_sens <- na.omit(specif_Y[po])
    #     if(any(is.infinite(d_i_s_sens ))){
    #       d_i_s_sens <- d_i_s_sens[-which(is.infinite(d_i_s_sens ))]
    #     }
    #     if(length(d_i_s_sens)>2){
    #       out_stat_sens_spec$`Specificity Y`[i_stat] <- paste("$",round(mean(d_i_s_sens,na.rm = T),2)," \\pm  ",
    #                                                           round(sd(d_i_s_sens,na.rm = T),2),"$",sep="")
    #     }else{
    #       out_stat_sens_spec$`Specificity Y`[i_stat] <- ""
    #     }
    #   }
    # }
    # out_stat$n <- paste(out_stat$n)
    # out_stat[,c(1,2)] <- out_stat[,c(2,1)]
    # names(out_stat)[c(1,2)] <- names(out_stat)[c(2,1)]
    # object <- xtable::xtable(x=out_stat[,-c(3,4,11,12)])
    # print(object,sanitize.text.function=function(x){x}, include.rownames=FALSE, file="/Users/hlorenzo/Dropbox/res_stat.txt")
  }
}
