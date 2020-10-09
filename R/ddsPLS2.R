#' Title
#'
#' @param Xs Xs
#' @param Y position of the legend
#' @param x.intersp legend x space
#' @param y.intersp legend y space
#'
#' @return
#' @export
#'
#' @useDynLib ddsPLS
predict_ddsPLS2 <- function(model,xs){


}


#' Title
#'
#' @param Xs Xs
#' @param Y Y
#' @param lam lam
#' @param tau tau
#' @param NZV NZV
#' @param Rs Rs
#'
#' @return
#' @export
#'
#' @useDynLib ddsPLS
ddsPLS2 <- function(Xs,Y,lam,tau=0.0975,method=2,NZV=1e-3,Rs=NA,standardize=F){
  ### Scale matrices
  MU_X_K = SD_X_K = xs0 <- list()
  K <- length(Xs)
  for(k in 1:K){
    MU_X_K[[k]] <- colMeans(Xs[[k]],na.rm = T)
    SD_X_K[[k]] <- apply(Xs[[k]],2,sd,na.rm = T)
    if(standardize){
      xs0[[k]] <- scale(Xs[[k]])
    }else{
      xs0[[k]] <- Xs[[k]]
    }
  }
  MU_Y <- colMeans(Y,na.rm = T)
  SD_Y <- apply(Y,2,sd,na.rm = T)
  if(standardize){
    y0 <- scale(Y)
  }else{
    y0 <- Y
  }
  K <- length(Xs)
  if(length(Rs)==0){
    Rs <- rep(1,K)
  }else if(length(Rs)==1){
    if(is.na(Rs)){
      Rs <- rep(1,K)
    }
  }
  Us = Ts <- list()
  n <- nrow(Y)
  for(k in 1:K){
    Us[[k]] <- matrix(0,ncol(Xs[[k]]),Rs[k])
    Ts[[k]] <- matrix(0,n,Rs[k])
  }
  h<-1
  beta_plus = t_K = t_plus = B = V_super = S_super = us <- list()
  for(k in 1:K){
    B[[k]] <- matrix(0,ncol(xs0[[k]]),ncol(y0))
  }
  B_r_out <- list()
  for(h in 1:max(Rs)){
    K_h <- which(Rs>=h)
    m_plus <- auto_ddsPLS(xs0[K_h],y0,lambdas = lam,plotVarB = F,
                          R_max = 1,tau=tau,scaleMat = F)
    for(i_k in 1:length(K_h)){
      k <- K_h[i_k]
      if(length(m_plus$model$model_per_block[[i_k]]$U_hat)!=0){
        Us[[k]][,h] <- m_plus$model$model_per_block[[i_k]]$U_hat
      }
      Ts[[k]][,h] <- xs0[[k]]%*%Us[[k]][,h]
    }
    if(length(K_h)>1){
      beta_h_plus <- m_plus$model$u_super
    }else{
      beta_h_plus <- matrix(1,1,1)
    }
    beta_plus[[h]] <- list(K_h=K_h,beta=beta_h_plus)
    l_beta <- length(beta_h_plus)
    t_K[[h]] = t_plus[[h]] <- list()
    if(length(K_h)>1){
      V_super[[h]] <- m_plus$model$model_super$RES$V
      if(!is.null(m_plus$model$model_super$RES$V)){
        S_super[[h]] <- y0%*%V_super[[h]]
      }
    }else{
      V_super[[h]] <- m_plus$model$model_per_block[[1]]$RES$V
      if(!is.null(V_super[[h]])){
        S_super[[h]] <- y0%*%V_super[[h]]
      }
    }
    for(i_k in 1:length(K_h)){
      u_i_k <- m_plus$model$model_per_block[[i_k]]$U_hat
      if(length(u_i_k)>0){
        t_K[[h]][[i_k]] <- list(K_h=K_h,t=xs0[[K_h[i_k]]]%*%u_i_k)
        t_plus[[h]][[i_k]] <- list(K_h=K_h,t=t_K[[h]][[i_k]]$t%*%beta_plus[[h]]$beta[[i_k]])
        B[[K_h[i_k]]] <- B[[K_h[i_k]]] + tcrossprod(u_i_k%*%beta_plus[[h]]$beta[[i_k]],V_super[[h]])
      }else{
        t_K[[h]][[i_k]] = t_plus[[h]][[i_k]] <- NA
      }
    }
    B_r_out[[h]] <- B[[1]]
    if(length(K_h)>1){
      y0 <- m_plus$model$model_super$RES$residuals$y
    }else{
      y0 <- m_plus$model$model_per_block[[1]]$RES$residuals$y
    }

    ### No deflation on X
    # for(i_k_h in 1:length(K_h)){
    #   k <- K_h[i_k_h]
    #   xs0[[k]] <- m_plus$model$model_per_block[[i_k_h]]$RES$residuals$x
    #   id_na_k <- id_na[[k]]
    #   if(length(id_na_k)!=0){
    #     xs0[[k]][id_na_k,] <- NA
    #   }
    # }
    ###

  }

  ## Do the prediction
  y_est <- matrix(rep(MU_Y,n),nrow = n,byrow = T)
  for(k in 1:K){
    mu_k <- matrix(rep(MU_X_K[[k]],n),nrow = n,byrow = T)
    if(standardize){
      sd_k <- matrix(rep(SD_X_K[[k]],n),nrow = n,byrow = T)
      sd_y <- matrix(rep(SD_Y,n),nrow = n,byrow = T)
      y_est <- y_est + (((xs0[[k]]+mu_k)/sd_k)%*%B[[k]])*SD_Y
    }else{
      y_est <- y_est + (((xs0[[k]]+mu_k))%*%B[[k]])
    }
  }

  ## Remove bad components
  posU0 <- lapply(Us,function(u){which(colSums(u^2)<NZV)})
  uOk = tOk <- list()
  for(k in 1:K){
    if(length(posU0[[k]])>0){
      uOk[[k]] <- Us[[k]][,-posU0[[k]],drop=F]
    }else{
      uOk[[k]] <- Us[[k]]
    }
    tOk[[k]] <- xs0[[k]]%*%uOk[[k]]
  }

  list(B=B,B_r=B_r_out,u=uOk,t=tOk,y_est=y_est,
       V_super=do.call(cbind,V_super),S_super=do.call(cbind,S_super),
       parameters=list(x=list(mu=MU_X_K,sd=SD_X_K),y=list(mu=MU_Y,sd=SD_Y)),residuals=list(y=y0))
}

get_model_ddsPLS2 <- function(x,y,ncomp,lams){
  p <- ncol(x)
  q <- ncol(y)
  n <- nrow(x)
  B <- matrix(0,p,q)
  U <- matrix(0,p,ncomp)
  P <- matrix(0,p,ncomp)
  Q <- matrix(0,q,ncomp)
  Ts <- matrix(0,n,ncomp)
  y0 <- y
  x0 <- x
  for(r in 1:ncomp){
    model <- ddsPLS2(Xs = list(x0),Y = y0,lam = lams[r],Rs = 1)
    if(length(model$u[[1]])!=0){
      U[,r] <- model$u[[1]]
      Ts[,r] <- model$t[[1]]
      P[,r] <- crossprod(x,model$t[[1]])/sum(model$t[[1]]^2)
      Q[,r] <- crossprod(y0,model$t[[1]])/sum(model$t[[1]]^2)
      Q[,r] <- tcrossprod(model$V_super)%*%Q[,r]
      if(r!=1){
        for(s_r in (r-1):1){
          U[,r] <- U[,r]-U[,s_r]*sum(P[,s_r]*U[,r])
        }
      }
      B <- B + tcrossprod(U[,r,drop=F],model$V_super)
      y0 <- y0 - tcrossprod(Ts[,r],Q[,r,drop=F])
      x0 <- x0 - tcrossprod(Ts[,r],P[,r,drop=F])
    }
  }
  list(U=U,Ts=Ts,B=B,P=P,residuals=list(e_x=x0,e_y=y0))
}
