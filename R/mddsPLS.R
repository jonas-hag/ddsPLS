#' The core function of the Multi-Data-Driven sparse PLS function
#'
#' This function should not be used directly by the user.
#'
#' @format A list containing the following objects:
#' \describe{
#'   \item{Xs}{A data-frame of a matrix or a list of data-frames or matrices of
#' \emph{n} rows each, the number of individuals. Some rows must be missing. The
#' different matrices can have different numbers of columns. The length of Xs is
#' denoted by \emph{K}.}
#'   \item{Y}{A matrix of n rows of a vector of length n detailing the
#' response matrix. No missing values are allowed in that matrix.}
#'   \item{lambda}{A real \eqn{[0,1]} where 1 means just perfect correlations
#' will be used and 0 no regularization is used.}
#'   \item{R}{A strictly positive integer detailing the number of components to
#' build in the model.}
#'   \item{mode}{A character chain. Possibilities are "\emph{reg}", which implies
#'  regression problem or anything else which means clustering is considered.
#'  Default is "\emph{reg}".}
#'   \item{verbose}{Logical. If TRUE, the function cats specificities about the
#' model. Default is FALSE.}
#' }
#'
#'
#' @return A list containing the following objects:
#' \describe{
#'   \item{u}{A list of length \emph{K}. Each element is a \emph{p_kXR} matrix : the
#'    weights per block per axis.}
#'   \item{v}{A \emph{qXR} matrix : the weights for the \emph{Y} part.}
#'   \item{ts}{A list of length \emph{R}. Each element is a \emph{nXK} matrix : the
#'    scores per axis per block.}
#'   \item{t}{A \emph{nXR} matrix, the final score of the \emph{X} part.}
#'   \item{s}{A \emph{nXR} matrix, the score of the \emph{Y} part.}
#'   \item{B}{A list of length \emph{K}. Each element is a \emph{p_kXq} matrix : the
#'    regression matrix per block.}
#'   \item{mu_x_s}{A list of length \emph{K}. Each element is a \emph{p_k} vector : the
#'    mean variables per block.}
#'   \item{sd_x_s}{A list of length \emph{K}. Each element is a \emph{p_k} vector : the
#'    standard deviation variables per block.}
#'   \item{mu_y}{A vector of length \emph{q} : the mean variables for \emph{Y} part.}
#'   \item{sd_y}{A vector of length \emph{q} : the standard deviation variables for \emph{Y} part.}
#'   \item{R}{Given as an input.}
#'   \item{q}{A non negatvie integer : the number of variables of \emph{Y} matrix. }
#'   \item{Ms}{A list of length \emph{K}. Each element is a \emph{qXp_k} matrix : the
#'    soft-thresholded empirical variance-covariance matrix \eqn{Y^TX_k/(n-1)}.}
#'   \item{lambda}{Given as an input.}
#' }
MddsPLS_core <- function(Xs,Y,lambda=0,R=1,mode="reg",verbose=FALSE){
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
  Ms <- lapply(1:K,function(k,Xs,Y,l,n){
    M0 <- crossprod(Y,Xs[[k]])/(n-1)
    M <- abs(M0) - l[k]
    M[which(M<0)] <- 0
    M <- sign(M0)*M
  },Xs,Y,lambda_in,n)
  if(verbose){
    N_max <- sum(unlist(lapply(Ms,function(m){length(which(colSums(abs(m))!=0))})))
    cat(paste("At most ",N_max," variable(s) can be selected",sep=""));cat("\n")
  }
  ## Solve optimization problem
  #### Inside problems
  u_t_r <- list()
  t_r <- list()
  z_r <- list()
  for(k in 1:K){
    if(norm(Ms[[k]])==0){
      svd_k <- list(v=matrix(0,
                             nrow = ncol(Ms[[k]]),
                             ncol = R))
    }
    else{
      svd_k <- svd(Ms[[k]],nu = 0,nv = R)
    }
    u_t_r[[k]] <- svd_k$v
    if(k==1){
      for(r in 1:R){
        t_r[[r]] <- matrix(NA,n,K)
        z_r[[r]] <- matrix(NA,q,K)
      }
    }
    for(r in 1:R){
      t_r[[r]][,k] <- Xs[[k]]%*%u_t_r[[k]][,r]
      z_r[[r]][,k] <- Ms[[k]]%*%u_t_r[[k]][,r]
    }
  }
  t <- matrix(NA,n,R)
  v <- matrix(0,q,R)

  # t_all <- do.call(cbind,t_r)
  t_all <- do.call(cbind,t_r)
  z_all <- do.call(cbind,z_r)
  svd_all <- svd(z_all,nu=R,nv=R)#svd(t_all,nu=R,nv=R)
  u <- svd_all$v
  v0 <- svd_all$u
  # v <- crossprod(Y,v0)
  # for(r in 1:R){
  #   v[,r] <- v[,r]/sqrt(sum(v[,r]^2))
  # }
  v <- v0#
  s <- Y%*%v0#Y%*%v
  t <-  t_all%*%u
  alphas <- rep(0,R)
  for(r in 1:R){
    n_t_2<-sum(diag(crossprod(t[,r])))
    if(n_t_2!=0){
      alphas[r] <- sum(diag(crossprod(s[,r],t[,r])))/n_t_2
    }else{
      alphas[r] <- 0
    }
  }

  if(mode=="reg"){
    B <- list()
    for(k in 1:K){
      beta_k <- u[(k-1)*R+1:R,,drop=FALSE]
      B[[k]] <- u_t_r[[k]]%*%beta_k
      for(r in 1:R){
        B[[k]][,r] <- B[[k]][,r]*alphas[r]
      }
      B[[k]]  <- tcrossprod(B[[k]],v)
    }
  }
  else{
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
        B <- MASS::lda(Y ~ ., data = dataf)
        B <- list(B=B,sds=sds)
      }
    }
    else{
      B <- MASS::lda(Y ~ ., data = dataf)
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
  list(u=u_t_r,v=v,ts=t_r,beta_comb=u,t=t,s=s,B=B,mu_x_s=mu_x_s,sd_x_s=sd_x_s,mu_y=mu_y,
       sd_y=sd_y,R=R,q=q,Ms=Ms,lambda=lambda)
}


#' Multi-Data-Driven sparse PLS function
#'
#' This function takes a set \eqn{X} of \eqn{K} matrices defining the same \eqn{n} individuals and a matrix \eqn{Y} defining also those individuals. According to the num-
#' ber of components \eqn{R}, the user fixes the number of components the model
#' must be built on. The coefficient lambda regularizes the quality of proximity to the data choosing to forget the least correlated bounds between
#' \eqn{X} and \eqn{Y} datasets.
#'
#' @format A list containing the following objects:
#' \describe{
#'   \item{Xs}{A data-frame of a matrix or a list of data-frames or matrices of
#' \emph{n} rows each, the number of individuals. Some rows must be missing. The
#' different matrices can have different numbers of columns. The length of Xs is
#' denoted by \emph{K}.}
#'   \item{Y}{A matrix of \emph{n} rows of a vector of length \emph{n} detailing the
#' response matrix. No missing values are allowed in that matrix.}
#'   \item{lambda}{A real \eqn{[0,1]} where 1 means just perfect correlations
#' will be used and 0 no regularization is used.}
#'   \item{R}{A strictly positive integer detailing the number of components to
#' build in the model.}
#'   \item{mode}{A character chain. Possibilities are "\emph{reg}", which implies
#'  regression problem or anything else which means clustering is considered.
#'  Default is "\emph{reg}".}
#'   \item{errMin_imput}{Positive real. Minimal error in the Tribe Stage of the
#' Koh-Lanta algorithm. Default is \eqn{1e-9}.}
#'   \item{maxIter_imput}{Positive integer. Maximal number of iterations in the
#' Tribe Stage of the Koh-Lanta algorithm. If equals to \eqn{0}, mean imputation is
#'  considered. Default is \eqn{5}.}
#'   \item{verbose}{Logical. If TRUE, the function cats specificities about the
#' model. Default is FALSE.}
#' }
#'
#' @return A list containing the following objects:
#' \describe{
#'   \item{mod}{A mddsPLS object, see
#'    }
#' }
#'
#' @export
#'
#' @seealso \code{\link{predict.mddsPLS}}
#'
#' @examples
#' # Classification example :
#' data("penicilliumYES")
#' X <- penicilliumYES$X
#' X <- scale(X[,which(apply(X,2,sd)>0)])
#' Y <- as.factor(unlist(lapply(c("Melanoconidiu","Polonicum","Venetum"),function(tt){rep(tt,12)})))
#' mddsPLS_model_class <- mddsPLS(Xs = X,Y = Y,lambda = 0.958,R = 2,mode = "clas",verbose = TRUE)
#'
#' # Regression example :
#' data("liver.toxicity")
#' X <- scale(liver.toxicity$gene)
#' Y <- scale(liver.toxicity$clinic)
#' mddsPLS_model_reg <- mddsPLS(Xs = X,Y = Y,lambda=0.9,R = 1, mode = "reg",verbose = TRUE)
mddsPLS <- function(Xs,Y,lambda=0,R=1,mode="reg",
                    errMin_imput=1e-9,maxIter_imput=50,
                    verbose=FALSE){
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
  if(length(unlist(id_na))==0){
    ## If ther is no missing sample
    mod <- MddsPLS_core(Xs,Y,lambda=lambda,R=R,mode=mode,verbose=verbose)
  }else{
    ## If ther are some missing samples
    for(k in 1:K){## ## Imputation to mean
      if(length(id_na[[k]])>0){
        mu_k <- colMeans(Xs[[k]],na.rm = TRUE)
        for(k_ik in 1:length(id_na[[k]])){
          Xs[[k]][id_na[[k]][k_ik],] <- mu_k
        }
      }
    }
    if(K>1){
      Xs_init <- Xs
      mod_0 <- MddsPLS_core(Xs,Y,lambda=lambda,R=R,mode=mode)
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
                model_here <- MddsPLS_core(Xs_i,Y_i_k,lambda=lambda)
                mod_i_k <- list(mod=model_here,R=R,mode=mode,maxIter_imput=maxIter_imput)
                class(mod_i_k) <- "mddsPLS"
                Xs[[k]][i_k,Var_selected_k] <- predict(mod_i_k,newX_i)
              }
            }
          }
          mod <- MddsPLS_core(Xs,Y,lambda=lambda,R=R,mode=mode)
          if(sum(abs(mod$s))*sum(abs(mod_0$s))!=0){
            err_cr <- rep(0,K)
            for(k in 1:K){
              err_cr[k] <- min(abs(diag(crossprod(mod$u[[k]],mod_0$u[[k]]))))
            }
            which_0 <- which(err_cr==0)
            if(length(which_0)>0){
              if(length(which_0)!=length(err_cr)){
                err <- 1-min(err_cr[-which_0])
              }else{
                err <- 0
              }
            }
            else{
              err <- 1-min(err_cr)
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
    mod <- MddsPLS_core(Xs,Y,lambda=lambda,R=R,mode=mode,verbose=verbose)
  }
  out <- list(mod=mod,Xs=Xs,Y_0=Y_0,lambda=lambda,mode=mode,
              maxIter_imput=maxIter_imput,has_converged=has_converged)
  class(out) <- "mddsPLS"
  out
}


#' The predict function of a mdd-sPLS model
#'
#' @param mod_0 A mdd-sPLS object, output from the mddsPLS function.
#' @param newX A data-set where individuals are described by the same as for mod_0
#'
#' @return Predicted values for those individuals
#' @export
#'
#' @examples
#' mod_0 <- mddsPLS(X,Y)
#' Y_est <- predict(mod_0,X)
predict.mddsPLS  <- function(mod_0,newX){
  fill_X_test <- function(mod_0,X_test){
    lambda <- mod_0$lambda
    R <- mod_0$mod$R
    id_na_test <- unlist(lapply(X_test,function(x){anyNA(x)}))
    mod <- mod_0$mod
    if(any(id_na_test)){
      ## Create covariable matrix train
      pos_ok <- which(!id_na_test)
      t_X_here <- do.call(cbind,lapply(1:R,function(ii,ti){
        ti[[ii]][,pos_ok]
      },mod$ts))
      u_X_here <- mod$u[pos_ok]
      mu_x_here <- mod$mu_x_s[pos_ok]
      sd_x_0 <- mod$sd_x_s[pos_ok]
      ## Create to be predicted matrix train
      pos_no_ok <- (1:K)[-pos_ok]
      pos_vars_Y_here <- lapply(mod$u[pos_no_ok],function(u){which(rowSums(abs(u))!=0)})
      if(sum(unlist(pos_vars_Y_here))!=0){
        nvars_Y_here_TOTAL <- length(unlist(pos_vars_Y_here))
        vars_Y_here <- matrix(0,nrow(t_X_here),nvars_Y_here_TOTAL)
        C_pos <- 1
        for(k_id in 1:length(pos_no_ok)){
          vars_k_id <- pos_vars_Y_here[[k_id]]
          if(length(vars_k_id)>0){
            vars_Y_here[,C_pos+(0:(length(vars_k_id)-1))] <- mod_0$Xs[[pos_no_ok[k_id]]][,vars_k_id,drop=FALSE]
            C_pos <- C_pos + length(vars_k_id)
          }
        }
      }
      else{
        vars_Y_here <- matrix(0,nrow(t_X_here),1)
      }
      ## Generate model
      model_impute_test <- mddsPLS(t_X_here,vars_Y_here,lambda = lambda,R = R,maxIter_imput = mod_0$maxIter_imput)
      ## Create test dataset
      n_test <- nrow(X_test[[1]])
      t_X_test <- matrix(NA,n_test,ncol(t_X_here))
      K_h <- sum(1-id_na_test)
      for(r_j in 1:R){
        for(k_j in 1:K_h){
          kk <- pos_ok[k_j]
          pos_col <- (r_j-1)*K_h+k_j
          xx <- X_test[[kk]]
          for(id_xx in 1:n_test){
            xx[id_xx,] <- xx[id_xx,]-mu_x_here[[k_j]]
            xx[id_xx,which(sd_x_0[[k_j]]!=0)] <-
              xx[id_xx,which(sd_x_0[[k_j]]!=0)]/sd_x_0[[k_j]][which(sd_x_0[[k_j]]!=0)]
          }
          t_X_test[,pos_col] <- xx%*%u_X_here[[k_j]][,r_j]
        }
      }
      ## Estimate missing values
      res <- predict(model_impute_test,t_X_test)
      ## Put results inside Xs
      C_pos <- 1
      for(k_id in 1:length(pos_no_ok)){
        vars_k_id <- pos_vars_Y_here[[k_id]]
        X_test[[pos_no_ok[k_id]]] <- matrix(mod$mu_x_s[[pos_no_ok[k_id]]],nrow = 1)
        if(length(vars_k_id)>0){
          X_test[[pos_no_ok[k_id]]][1,vars_k_id] <- res[C_pos+(0:(length(vars_k_id)-1))]
          C_pos <- C_pos + length(vars_k_id)
        }
      }
    }
    X_test
  }

  is.multi <- is.list(newX)&!(is.data.frame(newX))
  if(!is.multi){
    newX <- list(newX)
  }
  n_new <- nrow(newX[[1]])
  mod <- mod_0$mod
  q <- mod$q
  if(n_new==1){
    K <- length(newX)
    id_na_test <- unlist(lapply(newX,function(x){anyNA(x)}))
    if(any(id_na_test)){
      if(K>1 & mod_0$maxIter_imput>0){
        newX <- fill_X_test(mod_0,newX)
      }
      else{
        for(k in 1:K){
          if(id_na_test[k]){
            newX[[k]][1,] <- mod_0$mod$mu_x_s[[k]]
          }
        }
      }
    }
    mode <- mod_0$mode
    Y_0 <- mod_0$Y_0
    mu_x_s <- mod$mu_x_s
    sd_x_s <- mod$sd_x_s
    mu_y <- mod$mu_y
    sd_y <- mod$sd_y
    R <- mod$R
    K <- length(mu_x_s)
    for(k in 1:K){
      for(i in 1:n_new){
        newX[[k]][i,]<-(newX[[k]][i,]-mu_x_s[[k]])
        ok_sd <- which(sd_x_s[[k]]!=0)
        newX[[k]][i,ok_sd] <- newX[[k]][i,ok_sd]/sd_x_s[[k]][ok_sd]
      }
    }
    if(mode=="reg"){
      newY <- matrix(0,n_new,q)
      for(k in 1:K){
        newY <- newY + newX[[k]]%*%mod$B[[k]]
      }
      for(i in 1:n_new){
        newY[i,]<-newY[i,]*sd_y+mu_y
      }
    }
    else{
      t_r_new <- list()
      for(k in 1:K){
        if(k==1){
          for(r in 1:R){
            t_r_new[[r]] <- matrix(NA,n_new,K)
          }
        }
        for(r in 1:R){
          t_r_new[[r]][,k] <- newX[[k]]%*%mod_0$mod$u[[k]][,r]
        }
      }
      df_new <- data.frame(do.call(cbind,t_r_new)%*%mod_0$mod$beta_comb)
      colnames(df_new) <- paste("X",2:(ncol(df_new)+1),sep="")
      if(is.null(mod_0$mod$B)){
        newY <- list(class=sample(levels(mod_0$Y_0),size = 1,
                                  prob = table(mod_0$Y_0)/sum(table(mod_0$Y_0))))
      }
      else if(!is.null(mod_0$mod$B$sds)){
        pos_sds_0 <- 1+which(mod_0$mod$B$sds)
        newY <- predict(mod_0$mod$B,df_new[,c(1,pos_sds_0)])
      }else{
        newY <- predict(mod_0$mod$B,df_new)
      }
    }
  }
  else{
    newY <- matrix(NA,n_new,q)
    for(i_new in 1:n_new){
      newY[i_new,] <- predict(mod_0,lapply(newX,function(nx,ix){nx[ix,,drop=FALSE]},i_new))
    }
  }
  newY
}


