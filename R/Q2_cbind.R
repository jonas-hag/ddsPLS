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
do_one_component <- function(x0,y0,n,p,q,COV,abs_COV,max_COV,lam,NZV=1e-3){
  max_cov_y <- apply(abs_COV,1,max)
  max_cov_x <- apply(abs_COV,2,max)
  id_y_high <- which(max_cov_y>lam)
  id_x_high <- which(max_cov_x>lam)
  if(length(id_x_high)>0 & length(id_y_high)>0){
    COV_high <- COV[id_y_high,id_x_high,drop=F]
    abs_COV_high <- abs(COV_high)
    COV_COV_high <- abs_COV_high - lam
    COV_COV_high[which(COV_COV_high<0)] <- 0
    COV_COV_high <- COV_COV_high*sign(COV_high)
    # Do svd
    U0 <- matrix(0,p,1)
    V0 <- matrix(0,q,1)
    ## X part
    model_NULL <- svd(COV_high,nu = 1,nv = 1)
    # svd_XY <- svd(COV_COV_high,nv = 1,nu=1)#svd(COV_COV_high,nv = 1,nu=0)#
    u_x_no_std <- t(COV_COV_high)%*%model_NULL$u
    U0[id_x_high,] <- u_x_no_std/sqrt(sum(u_x_no_std^2))#svd_XY$v#
    ## Y part
    u_y_no_std <- COV_COV_high%*%model_NULL$v
    V0[id_y_high,] <- u_y_no_std/sqrt(sum(u_y_no_std^2))#svd_XY$u#
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
model_PLS <- function(x,y,lam,kX=NULL,deflatX=T,R=1,NZV=1e-3,to.scale=T){
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
  for(r in 1:R){
    COV <- crossprod(y0,x0)/(n-1)
    abs_COV <- abs(COV)
    max_COV <- max(na.omit(c(abs_COV)))
    lam_r <- lam
    if(length(lam)>1) lam_r <- lam[r]
    if(!is.null(kX)){
      kX_r <- kX
      if(length(kX)>1) kX_r <- kX[r]
      if(kX_r==p){
        lam_r <- 0
      }else{
        lam_r <- sort(apply(abs(COV),2,max),decreasing = F)[p-kX_r]
      }
    }
    if(lam_r<max_COV){
      c_h <- do_one_component(x0 = x0,y0 = y0,n = n,p = p,q = q,COV = COV,abs_COV = abs_COV,
                              max_COV=max_COV,lam = lam_r,NZV=NZV)
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
        y0 <- y0 - y_plus_un
        if(deflatX){
          x0 <- x0 - x0_plus
        }
      }
    }
  }
  if(to.scale){
    y_est <- y_est - mu_x%*%B
  }
  list(no_model=no_model,U_out=U_out,U_star=U_star,V_out=V_out,
       B=B,B_r=B_r,
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
do_loo <- function(xs0,y0,lam=0,kX=NULL,NCORES=1,deflatX=T){
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
    m_1 <-  model_PLS(x = X_train,y=Y_train,lam=lam,kX=kX,deflatX = deflatX,to.scale = F)
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
do_loo <- function(xs0,y0,lam=0,kX=NULL,NCORES=1,deflatX=T){
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
    m_1 <-  model_PLS(x = X_train,y=Y_train,lam=lam,kX=kX,deflatX = deflatX,to.scale = F)
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
#' @param x0 x0
#' @param y0 y0
#' @param lam lam
#' @param Rs Rs
#' @param NCORES number of cores
#'
#' @return
#' @export
#'
#' @useDynLib ddsPLS
do_boot <- function(x0,y0,Y_init,lam=0,deflatX=T){
  n <- nrow(x0)
  q <- ncol(y0)
  p <- ncol(x0)
  id_IB <- sample(1:n,size = n,replace = T)
  id_OOB <- (1:n)[-sort(unique(id_IB))]
  while(length(id_OOB)==0){
    id_IB <- sample(1:n,size = n,replace = T)
    id_OOB <- (1:n)[-sort(unique(id_IB))]
  }
  X_train <- x0[id_IB,,drop=F]
  X_test <- x0[id_OOB,,drop=F]
  Y_train <- y0[id_IB,,drop=F]
  Y_test <- y0[id_OOB,,drop=F]
  mu_k <- colMeans(X_train)
  sd_k <- apply(X_train,2,sd)
  mu_y <- colMeans(Y_train)
  sd_y <- apply(Y_train,2,sd)
  X_train <- do.call(cbind,lapply(1:p,function(j,xx,mu_k,sd_k){if(sd_k[j]>1e-9){out <- (xx[,j]-mu_k[j])/sd_k[j]}else{out = xx[,j]};out},X_train,mu_k,sd_k))
  Y_train <- do.call(cbind,lapply(1:q,function(j,xx,mu_k,sd_k){if(sd_k[j]>1e-9){out <- (xx[,j]-mu_k[j])/sd_k[j]}else{out = xx[,j]};out},Y_train,mu_y,sd_y))
  m_1 <-  model_PLS(x = X_train,y=Y_train,lam=lam,deflatX = deflatX,to.scale = F)
  B <- tcrossprod(m_1$U_out,m_1$V_out)
  var_selected_y <- matrix(0,length(id_OOB),q)
  var_selected_y[,which(colSums(abs(B))>1e-9)] <- 1
  X_test_normalize <- do.call(cbind,lapply(1:p,function(j,xx,mu_k,sd_k){if(sd_k[j]>1e-9){out <- (xx[,j]-mu_k[j])/sd_k[j]}else{out = xx[,j]};out},X_test,mu_k,sd_k))
  y_pred <- t(apply(X_test_normalize,1,function(xi){
    mu_y + (xi%*%B)*sd_y
  }))
  var_selected_y_train <- matrix(0,length(id_IB),q)
  var_selected_y_train[,which(colSums(abs(B))>1e-9)] <- 1
  # vars_expl_star = sum((X_train%*%B-Y_train)^2)#/sum(Y_init^2*var_selected_y_train)
  vars_expl_star <- sum(((X_train%*%B-Y_train)*var_selected_y_train)^2)/
    sum(((Y_init[id_IB,,drop=F])*var_selected_y_train)^2)
  y_0 <- t(apply(X_test_normalize,1,function(xi){
    mu_y
  }))
  PRESS_y <- (y_pred-Y_test)^2
  PRESS_y_star <- var_selected_y*PRESS_y
  RSS_y <- (y_0-Y_test)^2
  RSS_y_star <- RSS_y*var_selected_y
  list(id_OOB=id_OOB,Q2=1-sum(PRESS_y_star)/sum(RSS_y_star),
       PRESS=sum(PRESS_y_star),RSS=sum(RSS_y_star),
       vars_expl_star=vars_expl_star)
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
                            lambda_max=1,use_lambda=T,
                            n_B=20,
                            deflatX=T,
                            tau=0.0975,NCORES=1,center=T,
                            NZV=1e-3,verbose=F){
  ###### INTERNAL FUNCTION
  get_errors_stuff <- function(x0,y0,N_lambdas,lambdas_h,keepX_h_ordered,
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
      out <- do_loo(list(x0),y0,lam=lambdas_h[i_l],kX=keepX_h_ordered[i_l])
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
      m_1 <- model_PLS(x0,y0,lam,deflatX=deflatX,to.scale = F)
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
  sd_x_inv_mat <- matrix(rep(unlist(lapply(apply(do.call(cbind,Xs),2,sd),function(ss){if(abs(ss)>1e-9){out <- 1/ss}else{out <- 0};out})),q),ncol = q,byrow = T)
  sd_y_mat <- matrix(rep(apply(Y,2,sd),sum(ps)),ncol = q,byrow = T)
  sd_y_inv_mat <- matrix(rep(unlist(lapply(apply(Y,2,sd),function(ss){if(abs(ss)>1e-9){out <- 1/ss}else{out <- 0};out})),sum(ps)),ncol = q,byrow = T)
  sd_x_mat <- matrix(rep(apply(do.call(cbind,Xs),2,sd),q),ncol = q,byrow = T)
  sd_y_x_inv <- sd_y_mat * sd_x_inv_mat
  sd_x_y_inv <- sd_x_mat * sd_y_inv_mat
  RSS0 = RSS_h_moins_1 = RSS_h_moins_1_star <- colSums(Y_init^2)
  PRESS_y = RSS_y = Q2_y <- matrix(NA,n,q)
  PRESS_y_star = RSS_y_star = Q2_y_star <- matrix(NA,n,q)
  lambda = kX = Q2_sum_star <- rep(NA,n)
  ps <- unlist(lapply(Xs_init,ncol))
  x0 <- do.call(cbind,Xs_init)
  u <- matrix(NA,sum(ps),n)
  P <- matrix(0,sum(ps),n)
  B <- matrix(0,sum(ps),q)
  Us = Bs = B_r_out = Q2_h_sum_star = lambdas_out = vars_h_boot<- list()
  B_tot_LOO <- matrix(0,n,sum(ps)*q)

  ### For each 'h' component, look for the best lambda
  test <- T
  y0 <- Y_init
  h <- 0
  V = VAR_h_s <- NULL
  K_h <- 1:K
  # Check lambdas and stuff
  # coco <- abs(crossprod(Y_init,do.call(cbind,Xs_init)))/(n-1)
  lambdas_h <- seq(0,lambda_max,length.out = N_lambdas)#seq(0,min(max(coco),lambda_max),length.out = N_lambdas+2)[2:(N_lambdas+1)]#seq(min(coco),min(max(coco),lambda_max),length.out = N_lambdas+2)[2:(N_lambdas+1)]
  if(!use_lambda){
    lambdas_h_ordered <- sort(apply(coco,2,max),decreasing = F)
    keepX_h_ordered <- rev(0:(p-1))
    if(N_lambdas>=p){
      lambdas_h <- lambdas_h_ordered
      N_lambdas <- p
    }else{
      lambdas_h <- lambdas_h_ordered[(p-N_lambdas+1):p]
      keepX_h_ordered <- keepX_h_ordered[(p-N_lambdas+1):p]
    }
  }else{
    keepX_h_ordered <- NULL
  }
  ncomps = Q2_TOTAL <- rep(NA,N_lambdas)
  x0_deflated = y0_deflated <- list()
  while(test){
    h <- h + 1
    Bis <- list()
    pos_decoupe <- 1
    sample_id_B <- cbind(1:n,
                         matrix(sample(1:n,size = n*(n_B-1),replace = T),ncol = n_B-1))
    Q2_Boot = vars = Q2_h_sum <- rep(0,N_lambdas)
    RSS_il_y = var_sel = Q2_h_y = RSS_il_y_star = RSS_h_moins_1_star = PRESS_il_y <- matrix(0,N_lambdas,q)
    Q2_B = vars_in = Q2_h_clas <- matrix(NA,n_B,N_lambdas)
    ## LOO
    result_b <- get_errors_stuff(x0,y0,N_lambdas,lambdas_h,keepX_h_ordered,deflatX,RSS0,NCORES)
    RSS_il_y_star  <- result_b$RSS_il_y_star
    PRESS_il_y_star <- result_b$PRESS_il_y_star
    RSS_h_moins_1_star  <- result_b$RSS_h_moins_1_star
    # Q2_Boot <- result_b$Q2_h_star
    vars <-  result_b$vars
    RSS_il_y <- result_b$RSS_il_y
    var_sel <- result_b$var_sel
    Q2_h_y <- result_b$Q2_h_y
    Q2_h_sum <- result_b$Q2_h_sum
    PRESS_il_y <- result_b$PRESS_il_y
    ## BOOTSTRAP
    NCORES_w <- min(NCORES,n_B)
    `%my_do%` <- ifelse(NCORES_w!=1,{
      out<-`%dopar%`;cl <- makeCluster(NCORES_w)
      registerDoParallel(cl);out},{out <- `%do%`;out})
    Q2_star_bootstrat <- foreach(i_B=1:n_B,.packages = "ddsPLS",.combine=rbind) %my_do% {
      Q2_out <- matrix(NA,N_lambdas,6)
      for(i_l in 1:N_lambdas){
        out <- do_boot(x0,y0,lam=lambdas_h[i_l],Y_init=Y_init)
        # Q2_out[i_l,] <- c(i_l,out$PRESS,out$RSS)
        m_plus <- model_PLS(x0,y0,lambdas_h[i_l],R = 1,NZV=NZV,deflatX=deflatX,to.scale = F)
        mod_ok <- 0
        if(sum(m_plus$U_star^2)>1e-9){
          mod_ok <- 1
        }
        Q2_out[i_l,] <- c(i_l,out$Q2,out$vars_expl_star,length(out$id_OOB),mod_ok,i_B)
      }
      Q2_out
    }
    if(NCORES_w!=1)stopCluster(cl)
    q2_boot = vars_boot <- rep(NA,N_lambdas)
    mod_exists <- rep(0,N_lambdas)
    test_batchas <- TRUE
    # for(i_l in 1:N_lambdas){
    i_l <- 1
    while(test_batchas){
      pos <- which(Q2_star_bootstrat[,1]==i_l)
      toto <- Q2_star_bootstrat[pos,2:5]
      q2_boot[i_l] <- mean(toto[,1])
      vars_boot[i_l] <- mean(toto[,2])
      mod_exists[i_l] <- floor(sum(toto[,4])/length(toto[,4]))
      if(is.na(q2_boot[i_l])|is.nan(q2_boot[i_l])|i_l==N_lambdas|mod_exists[i_l]!=1){
        test_batchas <- F
      }
      i_l <- i_l + 1
    }
    vars <- vars_boot
    id_vars_loo_ok <- which(vars_boot>NZV & mod_exists==1)# & !is.infinite(result_b$Q2_h_star))

    if(F){
      if(length(id_vars_loo_ok)>0){

        plot_simus <- function(rere){
          i_B_s <- sort(unique(rere[,6]))
          i_l_s <- sort(unique(rere[,1]))
          oo <- matrix(NA,length(i_B_s),length(i_l_s))
          for(i_B in i_B_s){
            for(i_l in i_l_s){
              pos <- which( rere[,6]==i_B & rere[,1]==i_l )
              oo[i_B,i_l] <- rere[pos,2]
            }
          }
          oo
        }

        best_id_h <- which(q2_boot==max(na.omit(q2_boot[id_vars_loo_ok])))

        par(mfrow=c(1,2))
        # ok_simus <- plot_simus(Q2_star_bootstrat)
        # matplot(lambdas_h,ok_simus,col="gray",lwd=0.8,lty=1,type="l",ylim=c(0,1))
        plot(lambdas_h,q2_boot,type="l")#,ylim=c(0,1))
        points(lambdas_h,result_b$Q2_h_star,col="red",type="l")
        points(lambdas_h,mod_exists,col="brown",type="l",lty=2)
        abline(h=tau,col="blue")
        abline(h=NZV,col="red")
        abline(v=lambdas_h[best_id_h],col="orange")
        points(lambdas_h,vars_boot,type="l",col="pink",lwd=2)

        plot(lambdas_h,vars_boot+q2_boot,type="l",col="green",lwd=2)


        browser()
      }
      #which(vars>NZV)
      # q2_mean <- colMeans(Q2_star_bootstrat)
      # if(length(id_vars_ok)>0){
      #   par(mfrow=c(1,2))
      #   id_good <- which(q2_mean==max(na.omit(q2_mean[id_vars_ok])))
      #   matplot(lambdas_h,t(Q2_star_bootstrat),lty=1,type='l',col="gray")
      #   points(lambdas_h,q2_mean,col="red",type="l")
      #   points(lambdas_h,Q2_Boot,col="blue",type="l")
      #   abline(v=lambdas_h[id_good])
      #   plot(lambdas_h,vars,type="l")
      #   abline(v=lambdas_h[id_good])
      #   browser()
      # }
      # browser()
    }

    Q2_h_sum_star[[h]] <- q2_boot#result_b$Q2_h_star
    lambdas_out[[h]] <- lambdas_h
    vars_h_boot[[h]] <- vars_boot
    if(length(id_vars_loo_ok)>0){
      q2_max <- max(na.omit(q2_boot[id_vars_loo_ok]))
      id_s_cool <- which(q2_boot>=q2_max-NZV & q2_boot>tau &
                           vars_boot>NZV & mod_exists==1)
      if(length(id_s_cool)>0){
        best_id_h <- min(id_s_cool)
        if(length(best_id_h)>0){
          test_h <- Q2_h_sum_star[[h]][best_id_h]>tau
          best_lam_h <- lambdas_h[best_id_h]
          m_plus <- model_PLS(x0,y0,best_lam_h,R = 1,NZV=NZV,
                              deflatX=deflatX,to.scale = F)
          if(sum(m_plus$U_star^2)>1e-9){
            tested_bad_comp <- F
            if(h>1){
              print("-----")
              print(h)
              print(VAR_h_s)
              print(result_b$vars[best_id_h])
              if(min(VAR_h_s)<result_b$vars[best_id_h]){
                h <- min(which(VAR_h_s<result_b$vars[best_id_h]))
                VAR_h_s <- VAR_h_s[1:h]
                lambdas_h <- seq(min(lambdas_h),lambdas_h[best_id_h],length.out = N_lambdas)
                x0 <- x0_deflated[[h-1]]
                y0 <- y0_deflated[[h-1]]
                # V <- cbind(V[,1:h],V_r)
                u[,h] <- m_plus$U_star
                t_r <- x0%*%m_plus$U_out
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
                V_r <- crossprod(y0,t_r)/sum(t_r^2)#m_plus$V_out
                browser()
              }else{
                VAR_h_s <- c(VAR_h_s,result_b$vars[best_id_h])
              }
            }else{
              VAR_h_s <- c(VAR_h_s,result_b$vars[best_id_h])
            }
            if(!tested_bad_comp){
              #Q2_h_sum_star[[h]][best_id_h]>tau
              if(!test_h){
                test <- F
              }else{
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
                Q2_sum_star[h] <- Q2_h_sum_star[[h]][best_id_h]
                # Get the regression matrix of the optimal model
                u[,h] <- m_plus$U_star
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
                t_r <- x0%*%m_plus$U_out
              }
              B_r <- m_plus$B
              V_r <- m_plus$V_out
              y0_plus <- y0 - m_plus$e_y#x0%*%B_r
              id_y_sel <- which(abs(V_r)>NZV)
              if(verbose){
                cat(paste("\nComponent ",h," built with",
                          "\n          lambda=",round(best_lam_h,3),
                          "\n          kX=",kX[h],
                          "\n          Q_2^star=",round(Q2_h_sum_star[[h]][best_id_h],2),
                          "\n          var.expl.=",round(result_b$vars[best_id_h]*100,2),"%",sep=""))
              }
              B <- B + B_r*sd_y_x_inv
              B_r_out[[h]] <- B_r*sd_y_x_inv
              V <- cbind(V,V_r)
              y0 <- m_plus$e_y # y0 - y0_plus
              y0_deflated[[h]] <- y0
              if(deflatX){
                x0 <- m_plus$e_x # x0 - x0_plus
                x0_deflated[[h]] <- x0
              }
              lambdas_h <- seq(min(lambdas_h),best_lam_h,length.out = length(lambdas_h))
              lambdas_h <- unique(lambdas_h)
              N_lambdas <- length(lambdas_h)
              if(!use_lambda){
                coco <- abs(crossprod(y0,x0))/(n-1)
                lambdas_h_ordered <- sort(apply(coco,2,max),decreasing = F)
                keepX_h_ordered <- rev(0:(p-1))
                if(N_lambdas>=p){
                  lambdas_h <- lambdas_h_ordered
                  N_lambdas <- p
                }else{
                  lambdas_h <- lambdas_h_ordered[(p-N_lambdas+1):p]
                  keepX_h_ordered <- keepX_h_ordered[(p-N_lambdas+1):p]
                }
              }
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
      i_0 <- i_0+ps[k]
    }
    mu_y <- matrix(rep(colMeans(Y),n) ,nrow = n,byrow = T)
    y_est <- mu_y + x0_center%*%(B*sd_x_mat)
    # Compute LOO error to get Q2_reg
    lambda_sol <- lambda[1:h_opt]
    if(!use_lambda){
      kX_sol <- kX[1:h]
    }else{
      kX_sol <- NULL
    }
    Y_pred_all <- matrix(0,n,q)
    p <- ncol(x0_center)
    var_sel <- matrix(0,n,q)
    Xs_cbind <- as.matrix(do.call(cbind,Xs_init))
    ## Get LOO error
    for(ii in 1:n){
      X_train <- Xs_cbind[-ii,,drop=F]
      X_test <- Xs_cbind[ii,,drop=F]
      Y_train <- Y_init[-ii,,drop=FALSE]
      Y_test <- Y_init[ii,,drop=FALSE]
      mu_k_ii <- colMeans(X_train)
      mu_y_ii <- colMeans(Y_train)
      sd_x_inv_ii <- unlist(lapply(apply(X_train,2,sd),function(ss){if(abs(ss)>1e-9){out <- 1/ss}else{out <- 0};out}))
      # sd_x_mat_ii <- apply(X_train,2,sd)
      sd_y_mat_ii <- apply(Y,2,sd)
      # sd_y_inv_ii <- unlist(lapply(apply(Y,2,sd),function(ss){if(abs(ss)>1e-9){out <- 1/ss}else{out <- 0};out}))
      m_plus <- model_PLS(X_train,Y_train,lambda_sol,kX = kX_sol,R = h_opt,NZV=NZV,deflatX=deflatX)
      var_sel[ii,which(rowSums(abs(m_plus$V_out))>1e-9)] <- 1
      Y_pred_all[ii,] <- mu_y_ii + ((X_test-mu_k_ii))%*%m_plus$B
    }
    err_LOO <- (Y_init-Y_pred_all)
    ERRORS_LOO <- colSums(err_LOO^2)
    RSSO_loo <- colSums(scale(Y_init,center=center)^2)
    Q2_reg <- 1 - sum(ERRORS_LOO)/sum(RSSO_loo)
    Q2_reg_y <- 1 - ERRORS_LOO/RSSO_loo
    # STAR
    X_binded <- do.call(cbind,Xs_init)
    m_star <- model_PLS(X_binded,Y,lambda_sol,kX = kX_sol,R = h_opt,
                        NZV=NZV,deflatX=deflatX)
    ERRORS_LOO_star <- colSums((err_LOO^2)*var_sel)
    RSS0_star_model <- RSSO_loo[which(rowSums(abs(m_star$V_out))>1e-9)]
    Q2_reg_star <- 1 - sum(ERRORS_LOO_star)/sum(RSS0_star_model)
    # Prepare outputs
    optimal_parameters <- list(lambda=lambda_sol,kX=kX_sol,R=h_opt,
                               Q2_cum_y=Q2_cum_y,Q2_cum=Q2_cum,
                               Q2_cum_star=Q2_cum_star,Q2_sum_star=Q2_sum_star,
                               Q2_reg=Q2_reg,Q2_reg_y=Q2_reg_y,Q2_reg_star=Q2_reg_star,
                               ERRORS_LOO=ERRORS_LOO,ERRORS_LOO_star=ERRORS_LOO_star,
                               Y_pred_LOO=Y_pred_all)
    parameters <- list(RSS0=RSS0,RSS0_star=RSS0_star,RSS_y=RSS_y[1:(h_opt+1),],
                       PRESS_y=PRESS_y[1:(h_opt+1),],Q2_y=Q2_y[1:(h_opt+1),])
    bootstrap <- list(lambdas_h=lambdas_out,
                      Q2_h_star=Q2_h_sum_star,
                      vars_h_boot=vars_h_boot)
    out <- list(optimal_parameters=optimal_parameters,
                bootstrap=bootstrap,
                Us=Us,V=V,
                explained_variance=VAR_h_s*100,
                score_X = m_plus$score_x,
                Bs=Bs,B_cbind=B,B_r=B_r_out,
                y_est=y_est,parameters=parameters,
                mu_k=c(mu_x_s),mu_y=mu_y,sd_y=apply(Y,2,sd))
  }else{
    out <- NULL
  }
  out
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

