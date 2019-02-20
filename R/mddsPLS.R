#' The core function of the Multi-Data-Driven sparse PLS function.
#'
#' This function should not be used directly by the user.
#'
#' @param Xs A matrix, if there is only one block, or a list of matrices,, if there is more than one block, of \emph{n} rows each, the number of individuals. Some rows must be missing. The different matrices can have different numbers of columns. The length of Xs is denoted by \emph{K}.
#' @param Y A matrix of n rows of a vector of length n detailing the response matrix. No missing values are allowed in that matrix.
#' @param lambda A real \eqn{[0,1]} where 1 means just perfect correlations will be used and 0 no regularization is used.
#' @param R A strictly positive integer detailing the number of components to build in the model.
#' @param L0 An integer non nul parameter giving the largest number of X variables that can be selected.
#' @param mode A character chain. Possibilities are "\emph{reg}", which implies regression problem or anything else which means clustering is considered. Default is "\emph{reg}".
#' @param id_na A list of na indices for each block. Initialized to NULL.
#' @param verbose Logical. If TRUE, the function cats specificities about the model. Default is FALSE.
#'
#' @return A list containing the following objects:
#' \describe{
#'   \item{u}{A list of length \emph{K}. Each element is a \emph{p_kXR} matrix : the
#'    weights per block per axis.}
#'   \item{u_t_super}{A list of length \emph{K}. Each element is a \emph{p_kXR} matrix : the
#'    weights per block per axis scaled on the super description of the dataset. Denoted as
#'    \emph{scaled super-weights}.}
#'   \item{v}{A \emph{qXR} matrix : the weights for the \emph{Y} part.}
#'   \item{ts}{A list of length \emph{R}. Each element is a \emph{nXK} matrix : the
#'    scores per axis per block.}
#'   \item{(t,s)}{Two \emph{nXR} matrices, super-scores of the \emph{X} and \emph{Y} parts.}
#'   \item{(t_ort,s_ort)}{Two \emph{nXR} matrices, final scores of the \emph{X} and \emph{Y} part.
#'    They correspond to \emph{PLS} scores of \emph{(t,s)} scores and so \emph{t_ort^T s_ort} is diagonal,
#'    \emph{t_ort}, respectively \emph{s_ort}, carries the same information as \emph{t}, respectively \emph{s}.}
#'   \item{B}{A list of length \emph{K}. Each element is a \emph{p_kXq} matrix : the
#'    regression matrix per block.}
#'   \item{(mu_x_s,sd_x_s)}{Two lists of length \emph{K}. Each element is a \emph{p_k} vector : the
#'    mean and standard deviation variables per block.}
#'   \item{(mu_y,sd_y)}{Two vectors of length \emph{q} : the mean and the standard deviation variables for \emph{Y} part.}
#'   \item{R}{Given as an input.}
#'   \item{q}{A non negatvie integer : the number of variables of \emph{Y} matrix. }
#'   \item{Ms}{A list of length \emph{K}. Each element is a \emph{qXp_k} matrix : the
#'    soft-thresholded empirical variance-covariance matrix \eqn{Y^TX_k/(n-1)}.}
#'   \item{lambda}{Given as an input.}
#' }
#'
#' @importFrom stats sd model.matrix
#' @importFrom MASS lda
#'
MddsPLS_core <- function(Xs,Y,lambda=0,R=1,mode="reg",L0=NULL,verbose=FALSE,id_na=NULL){
  is.multi <- is.list(Xs)&!(is.data.frame(Xs))
  if(!is.multi){
    Xs <- list(Xs)
  }
  K <- length(Xs)
  ps <- lapply(Xs,ncol)
  ## Standardize Xs
  mu_x_s <- lapply(Xs,colMeans)
  sd_x_s <- lapply(Xs,function(X){apply(X,2,sd)})
  Xs <- lapply(Xs,scale)
  pos_0 <- lapply(sd_x_s,function(sdi){which(sdi==0)})
  if(length(unlist(pos_0))>0){
    for(i_0 in which(lapply(pos_0,function(pp){length(pp)})>0)){
      Xs[[i_0]][,pos_0[[i_0]]] <- 0
    }
  }
  ## Standardize Y
  Y_0 <- Y
  if(mode=="reg"){
    if(!(is.matrix(Y)|is.data.frame(Y))){
      Y <- as.matrix(Y)
    }
    if(is.data.frame(Y)){
      Y <- as.matrix(Y)
    }
  }
  else{
    Y_df <- data.frame(Y)
    Y <- scale(model.matrix( ~ Y - 1, data=Y_df))
  }
  mu_y <- colMeans(Y)
  sd_y <- apply(Y,2,sd)
  for(q_j in 1:length(sd_y)){
    if(sd_y[q_j]!=0){
      Y[,q_j] <- scale(Y[,q_j])
    }
  }
  q <- ncol(Y)
  n <- nrow(Y)
  ## Create soft-thresholded matrices
  lambda_in <- lambda
  if(length(lambda_in)==1){
    lambda_in <- rep(lambda_in,K)
  }
  if(!is.null(L0)){
    ps <- unlist(lapply(Xs,ncol))
    sum_ps <- sum(ps);cum_ps <- cumsum(c(0,ps))
    all_maxs <- rep(NA,sum_ps)
    for(k in 1:K){
      ii <- cum_ps[k]
      if(is.null(id_na)){
        c_k <- abs(cor(Y,Xs[[k]]))
      }else{
        if(length(id_na[[k]])>0){
          c_k <- abs(cor(Y[-id_na[[k]],,drop=F],Xs[[k]][-id_na[[k]],,drop=F]))
        }else{
          c_k <- abs(cor(Y,Xs[[k]]))
        }
      }
      all_maxs[ii+1:(ps[k])] <- apply(c_k,2,max)
      # if(is.null(id_na)){
      #   all_maxs[ii+1:(ps[k])] <- apply(abs(crossprod(Y,Xs[[k]])/(n-1)),MARGIN = 2,max)
      # }else{
      #   if(length(id_na[[k]])>0){
      #     all_maxs[ii+1:(ps[k])] <- apply(abs(crossprod(Y[-id_na[[k]],],Xs[[k]][-id_na[[k]],])/(n-1)),MARGIN = 2,max)
      #   }else{
      #     all_maxs[ii+1:(ps[k])] <- apply(abs(crossprod(Y,Xs[[k]])/(n-1)),MARGIN = 2,max)
      #   }
      # }
    }
    lambda_L0 <- sort(all_maxs,decreasing = T)[min(sum_ps,1+L0)]
    lambda_in <- rep(lambda_L0,K)
  }
  Ms <- lapply(1:K,function(k,Xs,Y,l,n){
    M0 <- crossprod(Y,Xs[[k]])/(n-1)
    M <- abs(M0) - l[k]
    M[which(M<0)] <- 0
    M <- sign(M0)*M
  },Xs,Y,lambda_in,n)
  if(verbose){
    N_max <- sum(unlist(lapply(Ms,function(m){length(which(colSums(abs(m))!=0))})))
    cat(paste("At most ",N_max," variable(s) can be selected in the X part",sep=""));cat("\n")
  }
  ## Solve optimization problem
  tete <- min(R,q)#min(R,min(unlist(lapply(Ms,dim))))
  if(R<1){
    stop("Choose R superior to 1",
         call. = FALSE)
  }
  if(tete!=R){
    # warning(paste("R modified to ",tete," due to dimension consraints",sep=""),
    #         call. = FALSE)
    R <- tete
  }
  if(floor(R)!=ceiling(R)){
    R <- floor(R)
    warning(paste("R not integer and estimated to ",R,sep=""),
            call. = FALSE)
  }
  #### Inside problems
  u_t_r = u_t_r_0 <- list()
  t_r <- list()
  z_r <- list()
  z_t=t_t <- list()
  # BETA_r <- list()
  for(k in 1:K){
    if(norm(Ms[[k]])==0){
      svd_k <- list(v=matrix(0,nrow = ncol(Ms[[k]]),ncol = R),
                    d=rep(0,R))
    }
    else{
      # R_k <- min(R,min(dim(Ms[[k]])))
      # svd_k <- svd(Ms[[k]],nu = 0,nv = R)
      R_init <- min(q,sum(ps))
      svd_k_init <- svd(Ms[[k]],nu = 0,nv = R_init)
      eigen_YXk <- apply(crossprod(Y,Xs[[k]])%*%svd_k_init$v,2,function(t)sum(t^2))
      eigen_YXk[which(svd_k_init$d==0)] <- 0
      ordo_YXk <- order(eigen_YXk,decreasing = T)[1:min(R,length(eigen_YXk))]
      svd_k <- list(v=svd_k_init$v[,ordo_YXk,drop=F],
                    d=svd_k_init$d[ordo_YXk])
      ## Complete coefficients if needed
      length_val_prop <- length(svd_k$d)
      if(length_val_prop<R){
        svd_k$d <- c(svd_k$d,rep(0,R-length_val_prop))
      }
      ## Complete basis if needed
      ncol_V <- ncol(svd_k$v)
      if(ncol_V<R){
        svd_k$v <- cbind(svd_k$v,matrix(0,nrow(svd_k$v),R-ncol_V))
      }
      ## Test coefficients and put corresponding vectors to 0 if needed
      for(r in 1:R){
        if(svd_k$d[r]==0){
          svd_k$v[,r] <- 0
        }
      }
    }
    u_t_r[[k]] = u_t_r_0[[k]] <- svd_k$v
    if(k==1){
      for(r in 1:R){
        t_r[[r]] <- matrix(0,n,K)
        z_r[[r]] <- matrix(0,q,K)
      }
    }
    for(r in 1:R){
      if(svd_k$d[r]!=0){
        t_r[[r]][,k] <- Xs[[k]]%*%u_t_r[[k]][,r]
        z_r[[r]][,k] <- Ms[[k]]%*%u_t_r[[k]][,r]
      }
    }
    z_t[[k]] <- Ms[[k]]%*%u_t_r[[k]]
    t_t[[k]] <- Xs[[k]]%*%u_t_r[[k]]#crossprod(Y,Xs[[k]]%*%u_t_r[[k]])
  }
  t <- matrix(NA,n,R)
  v <- matrix(0,q,R)
  ## Big SVD solution ######################### -----------------
  U_t_super <- list()
  Z <- do.call(cbind,z_t)
  R_opt <- min(dim(Z))
  svd_Z <- svd(Z,nu = R_opt,nv = R_opt)
  ## Reorder well ## -----
  crosY_T <- do.call(cbind,t_t)
  power <- colSums((crosY_T%*%svd_Z$v)^2)
  orderGood <- order(power,decreasing = T)
  svd_Z$d <- (svd_Z$d[orderGood])[1:R]
  svd_Z$u <- (svd_Z$u[,orderGood,drop=F])[,1:R,drop=F]
  svd_Z$v <- (svd_Z$v[,orderGood,drop=F])[,1:R,drop=F]
  ## END --- Reorder well ## -----
  beta_all <- svd_Z$v
  beta_list <- list()
  V_super <- Z%*%beta_all
  for(k in 1:K){#r in 1:R){
    #beta_list[[r]] <- beta_all[K*(r-1)+1:K,,drop=F]
    beta_list[[k]] <- beta_all[R*(k-1)+1:R,,drop=F]
  }
  v = V_super <- svd_Z$u
  if(length(svd_Z$d)<R){
    svd_Z$d <- c(svd_Z$d,rep(0,R-length(svd_Z$d)))
    v <- cbind(v,matrix(0,nrow=nrow(v),ncol=R-length(svd_Z$d)))
  }
  T_super <- matrix(0,nrow=n,ncol=R)
  for(k in 1:K){
    # U_t_super[[k]] <- matrix(0,nrow=nrow(u_t_r[[k]]),ncol=R)
    # for(r in 1:R){
    #   U_t_super[[k]] <- U_t_super[[k]] + u_t_r[[k]][,r,drop=F]%*%beta_list[[r]][k,,drop=F]
    # }
    U_t_super[[k]] <- u_t_r[[k]]%*%beta_list[[k]]
    T_super <- T_super + Xs[[k]]%*%U_t_super[[k]]
  }
  ###########################################################
  vars_current <- rep(0,R)
  for(r in 1:R){
    sc_r <- T_super[,r,drop=F]#scale(x$mod$t[,r,drop=F],scale=F)
    var_t_super_r <- sum(sc_r^2)
    if(var_t_super_r!=0){
      # VAR_SUPER_COMPS_ALL_Y[r] <- sum(abs(crossprod(sc_r,y_obs)/sqrt(var_t_super_r*VAR_TOT)))
      # b <- solve(crossprod(sc_r))%*%crossprod(sc_r,Y)
      # vars_current[r] <- (norm(sc_r%*%b,"f"))^2
      deno <- norm(tcrossprod(Y),'f')*norm(tcrossprod(sc_r),'f')
      vars_current[r] <- sum(diag(tcrossprod(Y)%*%tcrossprod(sc_r)))/deno
    }
  }
  l_cur <- length(vars_current)
  if(l_cur>1){
    ord <- order(vars_current,decreasing = T)
    if(sum(abs(ord-1:R))!=0){
      T_super <- T_super[,ord,drop=F]
      v=V_super <- v[,ord,drop=F]
      for(k in 1:K){
        U_t_super[[k]] <- U_t_super[[k]][,ord,drop=F]
        # u_t_r[[k]] <- u_t_r[[k]][,ord,drop=F]
        beta_list[[k]] <- beta_list[[k]][,ord,drop=F]
      }
    }
  }
  ###########################################################
  ## -------------------------- ######################### -----------------
  S_super <- Y%*%V_super
  T_S <- crossprod(T_super,S_super)
  T_T <- crossprod(T_super)
  # svd_ort <- svd(S_T,nu = R,nv = R)
  svd_ort_T_super <- svd(T_super,nu = 0,nv = R)
  # u_ort <- svd_ort_T_super$u
  v_ort <- svd_ort_T_super$v
  Delta_ort <- svd_ort_T_super$d^2
  if(sum(Delta_ort)!=0){
    t_ort <- T_super%*%v_ort
    s_ort <- S_super%*%v_ort
    D_0_inv <- matrix(0,nrow = length(Delta_ort),ncol = length(Delta_ort))
    del_0 <- which(Delta_ort<1e-10)
    if(length(del_0)>0){
      diag(D_0_inv)[-del_0] <- 1/Delta_ort[-del_0]
    }else{
      diag(D_0_inv) <- 1/Delta_ort
    }
    B_0 <- v_ort%*%tcrossprod(D_0_inv,v_ort)%*%T_S
    A <- matrix(0,R,R)
    for(r in 1:R){
      A[r,r] <- as.numeric(crossprod(t_ort[,r],s_ort[,r]))/as.numeric(crossprod(t_ort[,r]))
    }
  }else{
    t_ort=s_ort <- matrix(0,nrow = nrow(T_super),ncol=R)
    B_0 <- matrix(0,nrow = R,ncol=R)
    A <- matrix(0,R,R)
    V_super <- matrix(0,q,R)
  }
  u <- beta_all#beta# deprecated
  s <- S_super
  t <- T_super
  if(mode=="reg"){
    B <- list()
    count <- 1
    for(k in 1:K){
      B_k <- tcrossprod(U_t_super[[k]]%*%B_0,V_super)
      if(anyNA(B_k)){
        B_k <- matrix(0,nrow(B_k),ncol(B_k))
      }
      B[[k]] <- B_k
    }
  }else{
    dataf <- data.frame(cbind(Y_0,t));colnames(dataf)[1]<-"Y"
    for( cc in 2:ncol(dataf)){
      dataf[,cc] <- as.numeric(levels(dataf[,cc])[dataf[,cc]])
    }
    sds <- apply(dataf[,-1,drop=FALSE],2,sd)
    if(any(sds==0)){
      pos_sd0 <- as.numeric(which(sds==0))
      if(length(pos_sd0)==length(sds)){
        B <- NULL
      }else{
        dataf <- dataf[,-c(1+pos_sd0)]
        B <- lda(Y ~ ., data = dataf)
        B <- list(B=B,sds=sds)
      }
    }
    else{
      B <- lda(Y ~ ., data = dataf)
    }
  }
  if(verbose){
    a<-lapply(u_t_r,function(u){apply(u,2,function(u){length(which(abs(u)>1e-9))})})
    cat("    For each block of X, are selected in order of component:");cat("\n")
    for(k in 1:K){
      cat(paste("        @ (",paste(a[[k]],collapse = ","),") variable(s)",sep=""));cat("\n")
    }
    cat("    For the Y block, are selected in order of component:");cat("\n")
    cat(paste("        @ (",paste(apply(v,2,function(u){length(which(abs(u)>1e-9))}),
                                  collapse = ","),") variable(s)",sep=""));cat("\n")
  }
  list(u=u_t_r,u_t_super=U_t_super,v=v,ts=t_r,beta_comb=u,t=t,s=s,t_ort=t_ort,s_ort=s_ort,B=B,
       mu_x_s=mu_x_s,sd_x_s=sd_x_s,mu_y=mu_y,sd_y=sd_y,R=R,q=q,Ms=Ms,lambda=lambda_in[1])
}


#' Multi-Data-Driven sparse PLS function.
#'
#' This function takes a set \eqn{X} of \eqn{K} matrices defining the same \eqn{n} individuals and a matrix \eqn{Y} defining also those individuals. According to the num-
#' ber of components \eqn{R}, the user fixes the number of components the model
#' must be built on. The coefficient lambda regularizes the quality of proximity to the data choosing to forget the least correlated bounds between
#' \eqn{X} and \eqn{Y} datasets.
#'
#' @param Xs A matrix, if there is only one block, or a list of matrices,, if there is more than one block, of \emph{n} rows each, the number of individuals. Some rows must be missing. The different matrices can have different numbers of columns. The length of Xs is denoted by \emph{K}.
#' @param Y A matrix of \emph{n} rows of a vector of length \emph{n} detailing the response matrix. No missing values are allowed in that matrix.
#' @param lambda A real \eqn{[0,1]} where 1 means just perfect correlations will be used and 0 no regularization is used.
#' @param R A strictly positive integer detailing the number of components to build in the model.
#' @param L0 An integer non nul parameter giving the largest number of X variables that can be selected.
#' @param mode A character chain. Possibilities are "\emph{reg}", which implies  regression problem or anything else which means clustering is considered.  Default is "\emph{reg}".
#' @param errMin_imput Positive real. Minimal error in the Tribe Stage of the Koh-Lanta algorithm. Default is \eqn{1e-9}.
#' @param maxIter_imput Positive integer. Maximal number of iterations in the Tribe Stage of the Koh-Lanta algorithm. If equals to \eqn{0}, mean imputation is  considered. Default is \eqn{5}.
#' @param verbose Logical. If TRUE, the function cats specificities about the model. Default is FALSE.
#'
#' @return A list containing a mddsPLS object, see \code{\link{MddsPLS_core}}.
#'
#' @export
#'
#' @importFrom stats na.omit
#'
#' @seealso \code{\link{summary.mddsPLS}}, \code{\link{plot.mddsPLS}}, \code{\link{predict.mddsPLS}}, \code{\link{perf_mddsPLS}}, \code{\link{summary.perf_mddsPLS}}, \code{\link{plot.perf_mddsPLS}}
#'
#' @examples
#' # Single-block example :
#' ## Classification example :
#' data("penicilliumYES")
#' X <- penicilliumYES$X
#' X <- scale(X[,which(apply(X,2,sd)>0)])
#' Y <- as.factor(unlist(lapply(c("Melanoconidiu","Polonicum","Venetum"),function(tt){rep(tt,12)})))
#' mddsPLS_model_class <- mddsPLS(Xs = X,Y = Y,lambda = 0.958,R = 2,mode = "clas",verbose = TRUE)
#' summary(mddsPLS_model_class)
#'
#' ## Regression example :
#' data("liver.toxicity")
#' X <- scale(liver.toxicity$gene)
#' Y <- scale(liver.toxicity$clinic)
#' mddsPLS_model_reg <- mddsPLS(Xs = X,Y = Y,lambda=0.9,R = 1, mode = "reg",verbose = TRUE)
#' summary(mddsPLS_model_reg)
#'
#' # Multi-block example :
#' ## Classification example :
#' data("penicilliumYES")
#' X <- penicilliumYES$X
#' X <- scale(X[,which(apply(X,2,sd)>0)])
#' Xs <- list(X[,1:1000],X[,-(1:1000)])
#' Xs[[1]][1:5,]=Xs[[2]][6:10,] <- NA
#' Y <- as.factor(unlist(lapply(c("Melanoconidiu","Polonicum","Venetum"),function(tt){rep(tt,12)})))
#' mddsPLS_model_class <- mddsPLS(Xs = Xs,Y = Y,lambda = 0.95,R = 2,mode = "clas",verbose = TRUE)
#' summary(mddsPLS_model_class)
#'
#' ## Regression example :
#' data("liver.toxicity")
#' X <- scale(liver.toxicity$gene)
#' Xs <- list(X[,1:1910],X[,-(1:1910)])
#' Xs[[1]][1:5,]=Xs[[2]][6:10,] <- NA
#' Y <- scale(liver.toxicity$clinic)
#' mddsPLS_model_reg <- mddsPLS(Xs = Xs,Y = Y,lambda=0.9,R = 1, mode = "reg",verbose = TRUE)
#' summary(mddsPLS_model_reg)
mddsPLS <- function(Xs,Y,lambda=0,R=1,mode="reg",L0=NULL,
                    errMin_imput=1e-9,maxIter_imput=50,
                    verbose=FALSE){
  get_variances <- function(x){
    Xs <- x$Xs
    K <- length(Xs)
    y_obs <- x$Y_0
    R <- length(x$mod$ts)
    y_pred <- predict(x,Xs)
    mode <- x$mode
    if(mode=="reg"){
      if(!is.matrix(y_obs)){
        y_obs <- matrix(y_obs,ncol=1)
      }
      if(!is.matrix(y_pred)){
        y_pred <- matrix(y_pred,ncol=1)
      }
    }else{
      y_obs <- as.matrix(model.matrix( ~ y_obs - 1, data=data.frame(y_obs,ncol=1)))
      colnames(y_obs) <- levels(as.factor(x$Y_0))
      q <- ncol(y_obs)
      y_obs <- scale(y_obs,scale = F)
    }
    q <- ncol(y_obs)
    VAR_TOT=VAR_TOT_FROB <- norm(y_obs,"f")
    VAR_COMPS=VAR_COMPS_FROB <- matrix(0,K,R)
    VAR_SUPER_COMPS=VAR_SUPER_COMPS_FROB <- matrix(0,q,R)
    VAR_SUPER_COMPS_ALL_Y=VAR_SUPER_COMPS_ALL_Y_FROB <- matrix(0,1,R)
    for(r in 1:R){
      for(k in 1:K){
        t_r <- x$mod$ts[[r]]
        t_k_r <- scale(t_r[,k,drop=F],scale=F)
        var_t_k_r_all <- sum(t_k_r^2)
        if(var_t_k_r_all!=0){
          b <- solve(crossprod(t_k_r))%*%crossprod(t_k_r,y_obs)
          VAR_COMPS[k,r] <- (norm(t_k_r%*%b,"f")/VAR_TOT)^2
          deno <- norm(tcrossprod(y_obs),'f')*norm(tcrossprod(var_t_k_r_all),'f')
          VAR_COMPS_FROB[k,r] <- sum(diag(tcrossprod(y_obs)%*%tcrossprod(t_k_r)))/deno
        }
      }
    }
    for(j in 1:q){
      Y_j <- y_obs[,j,drop=F]
      var_j <- norm(Y_j,"f")
      if(var_j!=0){
        for(r in 1:R){
          t_super_r <- scale(x$mod$t[,r,drop=F],scale=F)
          var_t_super_r <- sum(t_super_r^2)
          if(var_t_super_r!=0){
            b <- solve(crossprod(t_super_r))%*%crossprod(t_super_r,Y_j)
            VAR_SUPER_COMPS[j,r] <- (norm(t_super_r%*%b,"f")/var_j)^2
            deno <- norm(tcrossprod(t_super_r),'f')*norm(tcrossprod(Y_j),'f')
            VAR_SUPER_COMPS_FROB[j,r] <- sum(diag(tcrossprod(Y_j)%*%tcrossprod(t_super_r)))/
                                                     deno
          }
        }
      }
    }
    for(r in 1:R){
      sc_r <- scale(x$mod$t[,r,drop=F],scale=F)
      var_t_super_r <- sum(sc_r^2)
      if(var_t_super_r!=0){
        b <- solve(crossprod(sc_r))%*%crossprod(sc_r,y_obs)
        VAR_SUPER_COMPS_ALL_Y[r] <- (norm(sc_r%*%b,"f")/VAR_TOT)^2
        deno <- norm(tcrossprod(y_obs),'f')*norm(tcrossprod(sc_r),'f')
        VAR_SUPER_COMPS_ALL_Y_FROB[r] <- sum(diag(tcrossprod(sc_r)%*%tcrossprod(y_obs)))/
                                            deno
      }
    }
    if(is.null(names(Xs))){
      legend_names_in <- paste("Block",1:K,sep=" ")
    }else{
      legend_names_in <- names(Xs)
      for(k in 1:K){
        if(nchar(names(Xs)[k])==0){
          legend_names_in[k] <- paste("Block",k)
        }
      }
    }
    rownames(VAR_COMPS)=rownames(VAR_COMPS_FROB) <- legend_names_in;
    colnames(VAR_COMPS)=colnames(VAR_COMPS_FROB) <- paste("Comp.",1:R)
    colnames(VAR_SUPER_COMPS)= names(VAR_SUPER_COMPS_ALL_Y) =
      colnames(VAR_SUPER_COMPS_FROB)= names(VAR_SUPER_COMPS_ALL_Y_FROB)<- paste("Super Comp.",1:R)
    rownames(VAR_SUPER_COMPS)=rownames(VAR_SUPER_COMPS_FROB) <- colnames(y_obs)
    return(list(
      Linear=list(VAR_SUPER_COMPS_ALL_Y=VAR_SUPER_COMPS_ALL_Y,
                  VAR_SUPER_COMPS=VAR_SUPER_COMPS,VAR_COMPS=VAR_COMPS),
      Frobenius=list(VAR_SUPER_COMPS_ALL_Y=VAR_SUPER_COMPS_ALL_Y_FROB,
                  VAR_SUPER_COMPS=VAR_SUPER_COMPS_FROB,VAR_COMPS=VAR_COMPS_FROB)))
  }
  if(lambda<0|lambda>1){
    stop("Choose lambda regularization parameter between 0 and 1",
         call. = FALSE)
  }
  is.multi <- is.list(Xs)&!(is.data.frame(Xs))
  if(!is.multi){
    Xs <- list(Xs)
  }
  K <- length(Xs)
  ps <- lapply(Xs,ncol)
  Y_0 <- Y
  if(!(is.matrix(Y)|is.data.frame(Y))){
    Y <- as.matrix(Y)
  }
  n <- nrow(Y)
  q <- ncol(Y)
  has_converged <- maxIter_imput
  id_na <- lapply(Xs,function(x){which(is.na(x[,1]),arr.ind = TRUE)})
  any_na_no_all <- lapply(Xs,function(x){
    oo <- which(is.na(x),arr.ind = TRUE)[,1]
    pi <- ncol(x)
    table_o <- table(oo)
    toto <- which(table_o!=pi)
    out <- NA
    if(length(toto)!=0){
      out <- names(table_o)[toto]
    }
    as.numeric(out)
  })
  if(length(na.omit(unlist(any_na_no_all)))!=0){
    which.block <- which(unlist(lapply(any_na_no_all,function(oo){length(na.omit(oo))!=0})))
    mess1 <- "Block(s) with values missing not for all the variables:\n"
    mess2 <- paste("(",paste(which.block,collapse=","),")\n",sep="",collapse=",")
    mess3 <- "Corresponding individuals for each block:\n"
    ouou <- paste(unlist(lapply(which.block,function(i){paste(
      "(",paste(any_na_no_all[[i]],collapse=",",sep=""),
      ")",sep="")})),collapse=",")
    mess4 <- paste(ouou,"\n",sep="",collapse="")
    stop(paste(mess1,mess2,mess3,mess4),
         call. = FALSE)
  }
  Models_init <- list()
  if(length(unlist(id_na))==0){
    ## If ther is no missing sample
    mod <- MddsPLS_core(Xs,Y,lambda=lambda,R=R,mode=mode,L0=L0,verbose=verbose)
  }else{
    if(!is.null(L0)){
      ps_init <- unlist(lapply(Xs,ncol))
      sum_ps_init <- sum(ps_init);cum_ps_init <- cumsum(c(0,ps_init))
      all_maxs_init <- rep(NA,sum_ps_init)
      for(k in 1:K){
        ii <- cum_ps_init[k]
        c_k <- abs(cor(Y,Xs[[k]],use = "pairwise"))
        all_maxs_init[ii+1:(ps_init[k])] <- apply(c_k,2,max)
        # if(is.null(id_na)){
        #   all_maxs_init[ii+1:(ps_init[k])] <- apply(abs(crossprod(Y,Xs[[k]])/(n-1)),MARGIN = 2,max)
        # }else{
        #   l_k <- length(id_na[[k]])
        #   if(l_k>0){
        #     all_maxs_init[ii+1:(ps_init[k])] <-
        #       apply(abs(crossprod(Y[-id_na[[k]],],Xs[[k]][-id_na[[k]],])/(n-1)),MARGIN = 2,max)
        #   }else{
        #     all_maxs_init[ii+1:(ps_init[k])] <- apply(abs(crossprod(Y,Xs[[k]])/(n-1)),MARGIN = 2,max)
        #   }
        # }
      }
      lambda_init <- sort(all_maxs_init,decreasing = T)[min(sum_ps_init,1+L0)]
    }else{
      lambda_init <- lambda
    }
    ## If ther are some missing samples
    for(k in 1:K){## ## Types of imputation for initialization
      if(length(id_na[[k]])>0){
        ## ## ## Imputation to mean
        # mu_k <- colMeans(Xs[[k]],na.rm = TRUE)
        # for(k_ik in 1:length(id_na[[k]])){
        #   Xs[[k]][id_na[[k]][k_ik],] <- mu_k
        # }
        ## ## ## Imputation with prediction by block
        ## Be careful if this is a classification case ##############################
        y_train <- Xs[[k]][-id_na[[k]],,drop=F]
        if(mode!="reg"){
          x_df <- data.frame(y=Y)
          x <- scale(model.matrix( ~ y - 1, data=x_df))
          x_train <- x[-id_na[[k]],,drop=F]
          x_test <- x[id_na[[k]],,drop=F]
        }else{
          x_train <- Y[-id_na[[k]],,drop=F]
          x_test <- Y[id_na[[k]],,drop=F]
        }
        model_init <- mddsPLS(x_train,y_train,R=R,lambda = lambda_init)
        y_test <- predict(model_init,x_test)
        Xs[[k]][id_na[[k]],] <- y_test
      }else{
        model_init <- mddsPLS(Y,Xs[[k]],R=R)
      }
      Models_init[[k]] <- model_init
    }
    if(K>1){
      # Xs_init <- Xs
      mod_0 <- MddsPLS_core(Xs,Y,lambda=lambda,R=R,mode=mode,L0=L0,id_na=id_na)
      if(sum(abs(as.vector(mod_0$s)))!=0){
        Mat_na <- matrix(0,n,K)
        for(k in 1:K){
          Mat_na[id_na[[k]],k] <- 1
        }
        err <- 2
        iter <- 0
        while(iter<maxIter_imput&err>errMin_imput){
          iter <- iter + 1
          for(k in 1:K){
            if(length(id_na[[k]])>0){
              no_k <- (1:K)[-k]
              i_k <- id_na[[k]]
              Xs_i <- mod_0$s[-i_k,,drop=FALSE]
              newX_i <- mod_0$s[i_k,,drop=FALSE]
              ## ## ## Look for selected variables
              Var_selected_k <- which(rowSums(abs(mod_0$u[[k]]))!=0)
              if(length(Var_selected_k)>0){
                ## ## ## ## Impute on the selected variables
                Y_i_k <- Xs[[k]][-i_k,Var_selected_k,drop=FALSE]
                model_here <- MddsPLS_core(Xs_i,Y_i_k,lambda=mod_0$lambda,R=R,L0=NULL)
                mod_i_k <- list(mod=model_here,R=R,mode="reg",maxIter_imput=maxIter_imput)
                class(mod_i_k) <- "mddsPLS"
                Xs[[k]][i_k,Var_selected_k] <- predict.mddsPLS(mod_i_k,newX_i)
              }
            }
          }
          mod <- MddsPLS_core(Xs,Y,lambda=mod_0$lambda,R=R,mode=mode,L0=L0)#NULL)#######################L0)#
          if(sum(abs(mod$t_ort))*sum(abs(mod_0$t_ort))!=0){
            err <- 0
            for(r in 1:R){
              n_new <- sqrt(sum(mod$t_ort[,r]^2))
              n_0 <- sqrt(sum(mod_0$t_ort[,r]^2))
              if(n_new*n_0!=0){
                err_i <- abs(1-as.numeric(abs(diag(crossprod(mod$t_ort[,r],mod_0$t_ort[,r]))))/(n_new*n_0))
                err <- err + err_i
              }
            }
          }
          else{
            err <- 0
          }
          if(iter>=maxIter_imput){
            has_converged <- 0
          }
          if(err<errMin_imput){
            has_converged <- iter
          }
          mod_0 <- mod
        }
      }
    }
    mod <- MddsPLS_core(Xs,Y,lambda=lambda,R=R,mode=mode,verbose=verbose,L0=L0)
  }
  out <- list(mod=mod,Xs=Xs,Y_0=Y_0,lambda=lambda,mode=mode,id_na=id_na,
              maxIter_imput=maxIter_imput,has_converged=has_converged,L0=L0)
  class(out) <- "mddsPLS"
  out$Variances <- get_variances(out)
  out
}

