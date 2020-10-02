# do_deflat <- function(x,y,y_init,lambda,R,NZV=1e-3){
#   n <- nrow(x)
#   p <- ncol(x)
#   q <- ncol(y)
#   v <- matrix(0,p,R)
#   u <- matrix(0,q,R)
#   t <- matrix(0,n,R)
#   s <- matrix(0,n,R)
#   B <- matrix(0,p,q)
#   bXr <- matrix(0,R,p)
#   bYr <- matrix(0,R,q)
#   # Phi_r <- diag(p)
#   # Phi_y_r <- diag(q)
#   x0 <- x/sqrt(n)
#   y0 <- y/sqrt(n)
#   COV <- crossprodC(y0,x0)
#   COV_th <- abs(COV)-lambda
#   pos_negat <- which(COV_th<0)
#   if(length(pos_negat)>0){
#     COV_th[pos_negat] <- 0
#   }
#   COV_th <- sign(COV)*COV_th
#   svd_r <- svd(crossprod(COV,COV_th),nu=1,nv=1)
#   test <- FALSE
#   if(svd_r$d[1]!=0) test <- TRUE
#   r <- 0
#   res_norm_2 <- list(X=rep(0,R),Y=rep(0,R))
#   while(test & r<R){
#     r <- r+1
#     res_norm_2$X[r] <- sum(x0^2)
#     res_norm_2$Y[r] <- sum(y0^2)
#     if(svd_r$d[1]!=0){
#       v[,r] <- svd_r$v
#       u[,r] <- COV_th%*%v[,r]  #svd_r$u
#       u[,r] <- u[,r]/sqrt(sum(u[,r]^2))
#       t[,r] <- mmultC(x0,v[,r,drop=F])
#       s[,r] <- mmultC(y0,u[,r,drop=F])
#       bXr[r,] <- (crossprod(t[,r],x0))/sum(t[,r]^2)
#       bYr[r,] <- t(u[,r]) #(crossprod(s[,r],y0))/sum(s[,r]^2)
#       defX <- t[,r,drop=F]%*%bXr[r,,drop=F]
#       defY <- t[,r,drop=F]%*%bYr[r,,drop=F]
#       x0 <- x0 - defX
#       y0 <- y0 - defY
#       COV <- crossprodC(y0,x0)
#       COV_th <- abs(COV)-lambda
#       pos_negat <- which(COV_th<0)
#       if(length(pos_negat)>0){
#         COV_th[pos_negat] <- 0
#       }
#       COV_th <- sign(COV)*COV_th
#       svd_r <- svd(crossprod(COV,COV_th),nu=1,nv=1)
#     }else{
#       test <- FALSE
#     }
#   }
#   U_out <- matrix(NA,p,r)
#   V_out <- matrix(NA,q,r)
#   U_out[,1] <- v[,1]
#   V_out[,1] <- u[,1]
#   if(ncol(u)>1){
#     for(r_s in 2:r){
#       id_m <- rev(1:(r_s-1))
#       u_cur <- v[,r_s,drop=F]
#       v_cur <- u[,r_s,drop=F]
#       for(m in id_m){
#         u_m <- v[,m]
#         bX_m <- bXr[m,]
#         u_cur <- u_cur - u_m*sum(bX_m*u_cur)
#         v_m <- u[,m]
#         bY_m <- bYr[m,]
#         v_cur <- v_cur - v_m*sum(bY_m*v_cur)
#       }
#       U_out[,r_s] <- u_cur
#       V_out[,r_s] <- v_cur
#     }
#   }
#   list(u=V_out,v=U_out,t=t,s=s,d=sqrt(colSums(t^2)),
#        res_norm_2=res_norm_2)
# }
#
#
# getOptimalRandLambda <- function(X,Y,lambdas,gamma=0,
#                                  tau=1e-2,R_max=NULL,
#                                  NZV=1e-3,
#                                  nb_na=0){
#   ## Functions --
#   get_u_and_others <- function(deflation,COV_0=NULL,
#                                lam,R_max,X=NULL,Y=NULL){
#     if(!deflation){
#       ### Construct Covariance matrix
#       COV_here <- COV_0
#       if(prod(dim(COV_0))!=1){
#         R_i_l <- min(nrow(COV_here),min(ncol(COV_here),R_max))
#       }else{
#         COV_here <- COV_here[1]
#         R_i_l <- 1
#       }
#       COV_here_0 <- abs(COV_here)-lam
#       test_cove_nul <- which(COV_here_0<0)
#       if(length(test_cove_nul)>0)  COV_here_0[test_cove_nul] <- 0
#       COV_here <- sign(COV_here)*COV_here_0
#       if(prod(dim(COV_0))==1) COV_here <- matrix(COV_here,1,1)
#       ### Construct SVD
#       svd_cov_here <- svd(crossprod(COV_0,COV_here),nu = R_i_l,nv=R_i_l)
#       deltas <- svd_cov_here$d[1:R_i_l]
#       U_here <- svd_cov_here$v
#       V_here <- svd_cov_here$u
#       t_here <- X%*%U_here
#       s_here <- Y%*%V_here
#       suppl <- list()
#     }else{
#       R_i_l <- min(nrow(X),R_max)
#       svd_cov_here <- do_deflat(X,Y,lambda = lam,R=R_i_l)
#       deltas <- svd_cov_here$d
#       U_here <- svd_cov_here$v
#       V_here <- svd_cov_here$u
#       t_here <- svd_cov_here$t
#       s_here <- svd_cov_here$s
#       suppl <- list(res_norm_2=svd_cov_here$res_norm_2)
#     }
#     comp_ok <- 1:length(deltas)
#     delta_null <- which(deltas<NZV)
#     if(length(delta_null)>0) comp_ok <- comp_ok[-delta_null]
#     B <- matrix(0,p,q)
#     U_here=V_here=B_0 <- NULL
#     if(length(delta_null)!=length(deltas)){
#       U_here <- svd_cov_here$v[,comp_ok,drop=F]
#       V_here <- svd_cov_here$u[,comp_ok,drop=F]
#       t_here <- t_here[,comp_ok,drop=F]
#       s_here <- s_here[,comp_ok,drop=F]
#       B_0 <- tcrossprod(solve(crossprod(t_here)),t_here)%*%s_here
#       B <- tcrossprod((U_here%*%B_0),V_here)
#     }
#     list(deltas=deltas[comp_ok],
#          R_i_l=min(R_i_l,max(comp_ok)),B=B,
#          u=U_here,v=V_here,t=t_here,s=s_here,suppl=suppl)
#   }
#   ## Main --
#   n <- nrow(Y)
#   COV_0 <- crossprod(Y,X)/n
#   p <- ncol(X); q <- ncol(Y)
#   lam <- lambdas[1]
#   n_lambdas <- length(lambdas)
#   R_max_in <- R_max
#   if(is.null(R_max_in)) R_max_in <- min(p,q)
#   alpha = MSE = varXU = resiX = resiY = resiF = varCA <- matrix(NA,n_lambdas,R_max_in)
#   alpha_j <- list()
#   for(j in 1:q){
#     alpha_j[[j]] <- matrix(NA,n_lambdas,R_max_in)
#   }
#   ## Compute alpha_0
#   F_t <- COV_0*n
#   norm_2_F_t <- sum(F_t^2)
#   alpha_R_0 <- norm_2_F_t/((n-nb_na)*q*p)
#
#   ## CONSTRUCT ALPHA
#   for(i_l in 1:n_lambdas){
#     lam <- lambdas[i_l]
#     if(lam<max(abs(COV_0))){
#       u_and_stuff <- get_u_and_others(
#         COV_0=COV_0,
#         lam=lam,R_max=R_max_in,X=X,Y=Y,y_init=)
#       R_i_l <- u_and_stuff$R_i_l
#       deltas <- u_and_stuff$deltas
#       ## build alpha per axis
#       for(r in 1:R_i_l){
#         U <- matrix(u_and_stuff$u[,1:r,drop=F],ncol=r)
#         V <- matrix(u_and_stuff$v[,1:r,drop=F],ncol=r)
#         t <- matrix(u_and_stuff$t[,1:r,drop=F],ncol=r)
#         s <- matrix(u_and_stuff$s[,1:r,drop=F],ncol=r)
#         pred <- tcrossprod(COV_0%*%U,U)*n
#         F_t <- COV_0*n-pred
#         MSE_i <- sum(F_t^2)/(n*p*q)
#         if(!deflation){
#           resi_Y <- sum((Y-tcrossprod(Y%*%V,V))^2)/(n*q)
#           resi_X <- sum((X-tcrossprod(X%*%U,U))^2)/(n*p)
#         }else{
#           resi_Y <- u_and_stuff$suppl$res_norm_2$Y/(n*q)
#           resi_X <- u_and_stuff$suppl$res_norm_2$X/(n*p)
#         }
#         expl_X <- sum((X%*%U)^2)/(n*p)
#         MSE[i_l,r] <- MSE_i#/(lam^2)
#         resiX[i_l,r] <- resi_X[r]
#         resiY[i_l,r] <- resi_Y[r]
#         varXU[i_l,r] <- (1-resi_X[r])#expl_X#expl_X*resi_Y
#         resiF[i_l,r] <- 1-(1-resi_Y[r])*(1-resi_X[r])
#         varCA[i_l,r] <- sum(deltas[1:r]^2)
#         norm_A_sur_p <- 1-resi_X[r]-expl_X
#         alpha[i_l,r] <- (MSE_i+
#                            gamma*(1-lam^2)*r+
#                            resiY[i_l,r]-resiX[i_l,r]
#                          -resiY[i_l,r]*resiX[i_l,r])
#
#         # (1-resiX[i_l,r])*
#         # resiY[i_l,r]+
#       }
#     }
#   }
#   ## Do not look if alpha worse than alpha_R_0
#   if(length(na.omit(c(alpha)))==0){
#     R_hat <- 0
#     lambda_hat <- NA
#     alpha_hat <- NA
#     alpha_mins <- NA
#   }else{
#     varXU[which(is.na(varXU))] <- 0
#     pos_var_nul <- which(varXU<NZV,arr.ind = T)
#     alp_work <- alpha
#     test <- TRUE
#     while(test){
#       alpha_vect <- na.omit(c(alp_work))
#       min_alphas <- min(alpha_vect)
#       alpha_range <- diff(range(alpha_vect))
#       al_ra_pot <- min_alphas+(tau)
#
#       pos_ok <- which(alp_work<=al_ra_pot,arr.ind = T)
#       R_hat <- min(pos_ok[,2])
#       lambda_pos_ok <- min(which(alp_work[,R_hat]<=al_ra_pot,arr.ind = T))
#       lambda_hat <- lambdas[lambda_pos_ok]
#       R_hat_bad <- which(pos_var_nul[,2]==R_hat)
#       lam_hat_bad <- which(pos_var_nul[,1]==lambda_pos_ok)
#       if(length(intersect(R_hat_bad,lam_hat_bad))!=0){
#         if(R_hat>1){
#           alp_work <- alpha[,1:(R_hat-1),drop=F]
#         }else{
#           R_hat = 0
#           lambda_hat = NA
#           test <- FALSE
#         }
#       }else{
#         test <- FALSE
#       }
#     }
#     alpha_hat <- alpha[lambda_pos_ok,R_hat]
#     alpha_vect <- na.omit(c(alpha))
#     alpha_mins <- lapply(1:ncol(alpha),function(r_i){
#       al_ri <- na.omit(alpha[,r_i])
#       if(length(al_ri)==0){
#         out <- max(alpha_vect)
#       }
#       else{out <- min(al_ri)}
#       out
#     })
#   }
#   if(R_hat!=0){
#     model_optim <- get_u_and_others(
#       COV_0=COV_0,
#       lam=lambda_hat,R_max=R_hat,X=X,Y=Y)
#     B <- model_optim$B
#     U <- model_optim$u
#     V <- model_optim$v
#     deltas <- model_optim$deltas
#     if(any(deltas<NZV)){
#       nul_var <- which(deltas<NZV)
#       if(min(nul_var)>1){
#
#       }
#     }
#   }else{
#     B <- NULL
#     U <- NULL
#     V <- NULL
#     deltas <- NULL
#   }
#   #### Prepare output
#   out <- list(lambda_hat=lambda_hat,R_hat=R_hat,alpha_hat=alpha_hat,
#               alpha_mins=alpha_mins,lambdas=lambdas,
#               alpha=alpha,MSE=MSE,varXU=varXU,
#               resiX=resiX,resiY=resiY,resiF=resiF,varCA=varCA,
#               alpha_R_0=alpha_R_0,U_hat=U,B_hat=B)
#   out
# }

# do_one_loop <- function(x,y,y_init,COV,COV_COV,F_mat=NULL,
#                         tau=1e-2){
#   svd_loop <- svd(COV_COV,nv = 1,nu=0)
#   d <- svd_loop$d[1]
#   all_variances <- sum(y_init^2)
#   p <- ncol(x)
#   q <- ncol(y)
#   U <- matrix(0,p,1)
#   V <- matrix(0,q,1)
#   if(d!=0){
#     U <- svd_loop$v
#     V0 <- COV%*%U
#     V <- V0/sqrt(sum(V0^2))
#     var_here <- sum((y_init%*%V0)^2)
#     if(var_here/all_variances>tau){
#     }
#   }
#   list(val_sg=d,U=U,V=V)
# }


#' Title
#'
#' @param res res
#' @param colsK colsK
#' @param ymin ymin
#' @param ymax ymax
#' @param sig_f sig_f
#'
#' @return
#' @export
#'
#' @useDynLib ddsPLS
plot_res_optim <- function(res,colsK="black",ymin=0,ymax=5,sig_f=NULL){
  cols <- RColorBrewer::brewer.pal(5,"Dark2")
  R_max_in <- ncol(res$alpha)
  if(res$R_hat!=0){
    cat(paste("Lambda_hat=",round(res$lambda_hat,3),
              " and R_hat=",res$R_hat,sep=""));cat("\n")
    cat(paste("Error empty model=",round(res$alpha_R_0,1),sep=""));cat("\n")
  }else{
    cat("No solution model found");cat("\n")
    cat(paste("Error empty model=",round(res$alpha_R_0,1),sep=""));cat("\n")
  }
  for( r in 1:R_max_in){
    plot(res$lambdas,res$alpha[,r],lty=1,col=colsK,
         type="l",lwd=2,ylim=c(ymin,ymax),
         main=paste("r=",r," over ",R_max_in,sep=""),
         ylab=expression(alpha[t]),xlab=expression(lambda))
    abline(h=c(0,1),v=0,lwd=2,col="gray")
    if(!is.null(sig_f)){
      ## plot sigma_f
      abline(v=sig_f,col=cols[2],lty=2)
      text(sig_f-0.05,4,labels = bquote(sigma[f]),col=cols[2])
    }
    if(res$R_hat!=0){
      col_R <- "blue" ; lty_R <- 1
      if(r!=res$R_hat){
        col_R <- cols[1] ; lty_R <- 2
      }
      abline(h=res$alpha_mins[r],lty=lty_R,col=col_R)
      text(0.2,y = res$alpha_mins[r]+0.2,labels = bquote(min~alpha[t]),
           col=col_R)
      if(r==res$R_hat){
        points(res$lambda_hat,res$alpha_hat,col="red",cex=1.7,pch=16)
      }
    }
  }
}


#' Title
#'
#' @param model model
#' @param names_blocks names_blocks
#' @param ymax ymax
#' @param ymax_super ymax_super
#' @param cexLegend cexLegend
#' @param pos pos
#'
#' @return
#' @export
#'
#' @useDynLib ddsPLS
plot_all_in_one <- function(model,names_blocks=NULL,
                            ymax_super=NULL,
                            cexLegend=1,pos="bottom",
                            singul=NULL){
  K <- length(model$Xs)
  R_all <- ncol(model$model$model_per_block[[1]]$alpha )
  if(is.null(names_blocks)){
    names_blocks <- paste("Block",1:K)
  }
  ymin <- min(unlist(lapply(model$model$model_per_block,
                            function(m){
                              toto <- na.omit(c(m$alpha))
                              out <- NULL
                              if(length(toto)>0){
                                out <- min(toto)
                              }
                            })))
  ymax <- max(unlist(lapply(model$model$model_per_block,
                            function(m){
                              toto <- na.omit(c(m$alpha))
                              out <- NULL
                              if(length(toto)>0){
                                out <- max(toto)
                              }
                            })))
  ymin_s <- min(na.omit(c(model$model$model_super$alpha)))
  bs_norm <- unlist(lapply(model$model$B,norm))
  cols <- RColorBrewer::brewer.pal(8,name = "Dark2")[-5]
  lty <- c(rep(1,7),rep(2,7),rep(3,7))
  # par(mfrow=c(2,round((K+1)/2)),mar=c(2.5,2.5,2.5,0.3),mgp=c(1.5,0.5,0))
  for(k in 1:K){
    mod_k <- model$model$model_per_block[[k]]
    matplot(mod_k$lambdas,mod_k$alpha,ylim=c(ymin,ymax),
            type = "l",lty = lty,lwd=2,col=cols,
            xlab=expression(lambda),
            ylab=expression(alpha[t]),
            main=paste("Super model",
                       ", R_opt=",mod_k$R_hat,
                       ", lam_opt=",round(mod_k$lambda_hat,2),sep="")
    )
    if(!is.null(singul)){
      abline(v=singul[[k]],lty=2)
    }
    pch <- 1
    if(bs_norm[k]!=0) pch <- 16
    points(mod_k$lambdas,
           mod_k$alpha_hat,col="red",pch=pch,cex=1.8)
    # legend(pos,legend = paste("R=",1:R_all,sep=""),
    #        col = cols,lty=lty,lwd = 2,bty="n",
    #        ncol=2,cex=cexLegend)
  }
  mod_super <- model$model$model_super
  if(is.null(ymax_super)){
    ymax_super <- ymax
  }
  if(length(mod_super)>0){
    matplot(mod_super$lambdas,mod_super$alpha,#ylim=c(ymin_s,ymax_super),
            type = "l",lty = lty,lwd=2,col=cols,
            xlab=expression(lambda),
            ylab=expression(alpha),
            main=paste("Super model",
                       ", R_opt=",mod_super$R_hat,
                       ", lam_opt=",mod_super$lambda_hat,sep="")
    )
    abline(h=c(0,1),lty=2,col="gray")
    pch <- 1
  }

}



#' Title
#'
#' @param model model
#' @param posLegend position of the legend
#' @param x.intersp legend x space
#' @param y.intersp legend y space
#'
#' @return
#' @export
#'
#' @useDynLib ddsPLS
plot_alpha <- function(model,posLegend="right",
                       x.intersp = 0.7,y.intersp = 0.7,
                       names_blocks = NULL){
  ## Reset personnal plot par() settings
  opar <- par(no.readonly =TRUE)
  on.exit(par(opar))
  ## -----------------------------------
  K <- length(model$model$model_per_block)
  if(is.null(names_blocks)){
    names_blocks <- paste("Block",1:K)
  }
  if(K+1<=8){
    cols=RColorBrewer::brewer.pal(K+1,name = "Dark2")
  }else{
    if(K+1<=12){
      cols=RColorBrewer::brewer.pal(K+1,name = "Paired")
    }else{
      cols=1:(K+1)
    }
  }
  dim_plot <- ceiling(sqrt(K+1))
  par(mfrow=c(dim_plot,dim_plot),mar=c(4,3,3,0.3),mgp=c(2,1,0))
  for(k in 1:K){
    alpha <- model$model$model_per_block[[k]]$alpha
    lambdas <- model$model$model_per_block[[k]]$lambdas
    if(length(na.omit(alpha))!=0){
      lambdas <- model$model$model_per_block[[k]]$lambdas
      lam <- model$model$model_per_block[[k]]$lambda_hat
      # alpha <- model$model$model_per_block[[k]]$alpha_hat
      # V12 <- model$model$model_per_block[[k]]$V12
      R <- model$model$model_per_block[[k]]$R_hat
      ylim <- range(na.omit(alpha))

      if(is.numeric(lam)){
        main <- bquote(.(names_blocks[k])~", R"==.(R)~
                         ", "~lambda==.(round(lam,2) ))
      }else{
        main <- bquote(.(names_blocks[k])~", Empty model")
      }
      plot(lambdas,
           alpha,
           ylab="Expl. Var.",ylim=c(ylim[1],ylim[2]),
           type="l",col=cols[k],lwd=2,
           main=main,
           xlab="")
      points(lambdas,alpha,lwd=1,type="l")
      mtext(expression(lambda),side=1,line=2)
      # points(lambdas,alpha,type="p",pch=k,col=cols[k])
      pch <- 1;if(length(model$model$selected_x[[k]])>0)pch <- 16
      if(is.numeric(lam)){
        points(lam,alpha[which(lambdas==lam)],col=cols[k],pch=pch,cex=2)
        points(lam,alpha[which(lambdas==lam)],pch=8)
      }
      # legend(posLegend,legend = c(expression(alpha),expression(v[1]),
      #                             expression(v[1~2])),
      #        lty=1:3,col=cols[k],x.intersp = x.intersp,y.intersp = y.intersp,
      #        bty="n",lwd=c(2,1,1))
    }else{
      plot(lambdas,
           rep(0,length(lambdas)),
           ylab="",yaxt="n",
           type="l",col="white",lwd=2,
           main=bquote(.(names_blocks[k])~", Nothing kept"),
           xlab="")
    }
  }
  alpha <- model$model$model_super$alpha
  if(length(na.omit(alpha))!=0){
    lambdas <- model$model$model_super$lambdas
    lam <- model$model$model_super$lambda_hat
    alpha <- model$model$model_super$alpha_hat
    # V12 <- model$model$model_super$V12
    R <- model$model$model_super$R_hat
    ylim <- range(na.omit(alpha))
    if(is.numeric(lam)){
      main <- bquote("Super block, R"==.(R)~
                       ", "~lambda==.(round(lam,2) ))
    }else{
      main <- bquote("Super block, Empty model")
    }
    plot(lambdas,
         alpha,
         ylab="Expl. Var.",ylim=c(ylim[1],ylim[2]),
         type="l",col=cols[k+1],lwd=2,
         main=main,
         xlab="")
    points(lambdas,alpha,lwd=1,type="l")
    mtext(expression(lambda),side=1,line=2)
    # points(lambdas,alpha,type="p",pch=k+1,col=cols[k+1])
    if(is.numeric(lam)){
      points(lam,alpha[which(lambdas==lam)],col=cols[k+1],pch=16,cex=2)
      points(lam,alpha[which(lambdas==lam)],pch=8)
    }
    # legend(posLegend,legend = c(expression(alpha),expression(v[1]),
    #                             expression(v[1~2])),
    #        x.intersp = x.intersp,y.intersp = y.intersp,
    #        lty=1:3,col=cols[k+1],bty="n",lwd=c(2,1,1))
  }else{
    if(is.numeric(lam)){
      main <- bquote("Super block, R"==.(R)~
                       ", "~lambda==.(round(lam,2) ))
    }else{
      main <- bquote("Super block, Empty model")
    }
    plot(lambdas,
         rep(0,length(lambdas)),
         ylab="Expl. Var.",
         type="l",col="white",yaxt="n",lwd=2,
         main=main,
         xlab="")
  }
}

#' Title
#'
#' @param model model
#' @param posLegend position of the legend
#' @param x.intersp legend x space
#' @param y.intersp legend y space
#'
#' @return
#' @export
#'
#' @useDynLib ddsPLS
plot_var_B <- function(model,posLegend="bottomright",
                       plotSuper=TRUE,
                       x.intersp = 0.7,y.intersp = 0.7,
                       names_blocks = NULL,low_lambda=NULL){
  ## Reset personnal plot par() settings
  opar <- par(no.readonly =TRUE)
  on.exit(par(opar))
  ## -----------------------------------
  K <- length(model$model$model_per_block)
  if(is.null(names_blocks)){
    names_blocks <- paste("Block",1:K)
  }
  if(K+1<=8){
    cols=RColorBrewer::brewer.pal(max(3,K+1),name = "Dark2")
  }else{
    if(K+1<=12){
      cols=RColorBrewer::brewer.pal(K+1,name = "Paired")
    }else{
      cols=1:(K+1)
    }
  }
  ylim <- range(na.omit(unlist(lapply(model$model$model_per_block,
                                      function(mm){
                                        oo <- na.omit(mm$alpha_hat)
                                        out <- NA
                                        if(length(oo)>0){
                                          out <- range(oo)
                                        }
                                        out
                                      }))))
  xlim <- range(na.omit(unlist(lapply(model$model$model_per_block,
                                      function(mm){
                                        oo <- na.omit(mm$lambdas)
                                        out <- NA
                                        if(length(oo)>0){
                                          out <- range(oo)
                                        }
                                        out
                                      }))))
  if(plotSuper){
    par(mfrow=c(1,2),mar=c(3.5,3.5,3,0.3),mgp=c(2,1,0))
    main <- expression("Covariate blocks")
  }else{
    par(mfrow=c(2,1),mar=c(3,3,3,0.3),mgp=c(2,1,0))
    main <- expression(Var[B](lambda)~" for missing samples removed")
  }
  plotted_first <- F
  for(k in 1:K){
    alpha <- model$model$model_per_block[[k]]$alpha
    lambdas <- model$model$model_per_block[[k]]$lambdas
    if(length(na.omit(alpha))!=0){
      lambdas <- model$model$model_per_block[[k]]$lambdas
      lam <- model$model$model_per_block[[k]]$lambda_hat
      R <- model$model$model_per_block[[k]]$R_hat
      if(!plotted_first){
        plot(lambdas,
             alpha,
             ylab=expression(Var[B](lambda)),ylim=ylim,xlim=xlim,
             type="l",col=cols[k],lwd=2,
             main=main,
             xlab="")
        plotted_first <- T
      }else{
        points(lambdas,col=cols[k],lwd=2,type="l",
               alpha)

      }
      mtext(expression(lambda),side=1,line=2)
      if(plotSuper){
        pch <- 1;if(length(model$model$selected_x[[k]])>0)pch <- 16
        if(is.numeric(lam)){
          points(lam,alpha[which(lambdas==lam)],col=cols[k],pch=pch,cex=2)
          points(lam,alpha[which(lambdas==lam)],pch=8)
        }
      }
    }
  }
  legend(posLegend,
         legend = names_blocks,
         col=cols,
         lty=1,
         x.intersp = x.intersp,y.intersp = y.intersp,
         bty="n",lwd=2)
  if(plotSuper){
    alpha <- model$model$model_super$alpha
    if(length(na.omit(alpha))!=0){
      lambdas <- model$model$model_super$lambdas
      lam <- model$model$model_super$lambda_hat
      alpha <- model$model$model_super$alpha_hat
      R <- model$model$model_super$R_hat
      ylim <- range(na.omit(alpha))
      if(is.numeric(lam)){
        main <- bquote("Super block, R"==.(R)~
                         ", "~lambda==.(round(lam,2) ))
      }else{
        main <- bquote("Super block, Empty model")
      }
      plot(lambdas,
           alpha,
           ylab=expression(Var[B](lambda)),ylim=c(ylim[1],ylim[2]),
           type="l",col=cols[k+1],lwd=2,
           main=main,
           xlab="")
      points(lambdas,alpha,lwd=1,type="l")
      mtext(expression(lambda),side=1,line=2)
      if(is.numeric(lam)){
        points(lam,alpha[which(lambdas==lam)],col=cols[k+1],pch=16,cex=2)
        points(lam,alpha[which(lambdas==lam)],pch=8)
      }
    }else{
      main <- bquote("Super block, Empty model")
      plot(lambdas,
           rep(0,length(lambdas)),
           ylab="Expl. Var.",
           type="l",col="white",yaxt="n",lwd=2,
           main=main,
           xlab="")
    }

  }else{
    ylim <- range(na.omit(unlist(lapply(model$model$model_per_block,
                                        function(mm){
                                          oo <- na.omit(mm$R)
                                          out <- NA
                                          if(length(oo)>0){
                                            out <- range(oo)
                                          }
                                          out
                                        }))))
    for(k in 1:K){
      Rs <- model$model$model_per_block[[k]]$R
      lambdas <- model$model$model_per_block[[k]]$lambdas
      if(length(na.omit(Rs))!=0){
        lambdas <- model$model$model_per_block[[k]]$lambdas
        if(k==1){
          plot(lambdas,
               Rs,
               ylab=expression(R),ylim=ylim,xlim=xlim,
               type="l",col=cols[k],lwd=2,
               main=expression("Number of components"),
               xlab="",log="y")
          mtext(expression(lambda),side=1,line=2)
        }else{
          points(lambdas,col=cols[k],lwd=2,type="l",
                 Rs)

        }
        low_k <- low_lambda[[k]]
        if(length(low_k)>0){
          points(lambdas[low_k],Rs[low_k],
                 type="l",lwd=3,col="red")
        }
      }
    }
    max_all <- lapply(low_lambda,function(aa){
      if(length(aa)>0) max(aa)})
    lambda_max <- lapply(
      1:K,
      function(k){
        ls <-  model$model$model_per_block[[k]]$lambdas
        ls[max_all[[k]]+1]
      }
    )
    abline(v=max(unlist(lambda_max)),col="red",lwd=1,lty=2)
    l_m <- round(max(unlist(lambda_max)),digits = 2)
    text(labels = bquote(lambda[min]==.(l_m)),pos = 4,
         x = l_m*1.05,y = (ylim[2]-ylim[1])/2,col="red" )
  }
}

## One component ##
do_one_component <- function(x0,y0,method=2,n,p,q,COV,abs_COV,max_COV,lam,tau=1e-2,NZV=1e-3){
  max_cov_y <- apply(abs_COV,1,max)
  max_cov_x <- apply(abs_COV,2,max)
  id_y_high <- which(max_cov_y>lam)
  id_x_high <- which(max_cov_x>lam)
  if(p-length(id_x_high)>=0 & q-length(id_y_high)>=0){
    COV_high <- COV[id_y_high,id_x_high,drop=F]
    abs_COV_high <- abs(COV_high)
    COV_COV_high <- COV_high - lam
    COV_COV_high[which(COV_COV_high<0)] <- 0
    if(method==2){
      COV_COV_high <- crossprod(COV_high,COV_COV_high*sign(COV_high))
      # COV_COV <- getCOV_COV(COV,lam)
    }else{
      COV_COV_high <- crossprod(COV_COV_high*sign(COV_high))
    }
    svd_loop <- svd(COV_COV_high,nv = 1,nu=0)
    # svd_loop <- rsvd(COV_COV,nv = 1,nu=0,k=1)# Ne pas utiliser rsvd car pas stable...
    U0 <- matrix(0,p,1)
    U0[id_x_high,] <- svd_loop$v
  }else{
    U0 <- matrix(0,p,1)
  }
  t <- x0%*%U0
  V_svd0 <- crossprod(y0,t)
  norm_t_0 <- sum(t^2)
  if(norm_t_0>NZV){
    proj_comp <- tcrossprod(t,V_svd0/norm_t_0)
    var_comp <- sum(proj_comp^2)
    #Build descriptor
    B_r <- tcrossprod(U0[id_x_high,,drop=F],V_svd0/norm_t_0)
    numer <- sum(diag(crossprod(B_r,COV_COV_high%*%B_r)))^2
    Y_here <- y0[,id_y_high,drop=F]/n
    V_svd <- V_svd0/norm_t_0
  }else{
    U0 <- matrix(0,p,1)
    V_svd <- matrix(0,q,1)
  }
  ##
  list(t=t,U0=U0,V_svd=V_svd)
}
########

## Model ##
model_PLS <- function(x,y,lam,tau=1e-2,method=2,R=1,NZV=1e-3){
  p <- ncol(x)
  q <- ncol(y)
  n <- nrow(y)
  x0 <- x
  y0 <- y
  U_out <- matrix(0,p,R)
  V_out <- matrix(0,q,R)
  bXr <- matrix(0,R,p)
  bYr <- matrix(0,R,q)
  e <- x0
  stop <- F
  no_model <- F
  for(r in 1:R){
    COV <- crossprod(y0,x0)/n
    abs_COV <- abs(COV)
    max_COV <- max(na.omit(abs_COV))
    if(lam<max_COV){
      c_h <- do_one_component(x0,y0,method=method,n,p,q,COV,abs_COV,max_COV,lam,tau=tau,NZV=NZV)
      t <- c_h$t ; U0  <- c_h$U0 ; V_svd  <- c_h$V_svd
      ## DEFLAT ##
      if(sum(U0^2)>NZV){
        bt <- crossprod(t,x0)/sum(t^2)
        e <- x0 # - t%*%bt  #  On teste sans deflation sur x : westerhuis and smilde 2001
        U_out[,r] <- U0
        V_out[,r] <- V_svd
        bXr[r,] <- bt
        bYr[r,] <- t(V_svd)
        y0 <- y0 - tcrossprod(t,V_out[,r,drop=F])
      }
    }
  }
  list(no_model=no_model,U_out=U_out,V_out=V_out,bXr=bXr,bYr=bYr,e_x=e,e_y=y0)
}
########

deflat_in_loop <- function(x,y,R_max=NA,lam,method=2,
                           lam_min=0,NZV=1e-3,tau=1e-2){
  ### ---SCRIPT--- ###
  alpha=angle=term1=term2 =norm_th=norm_est=term_croise<- NA
  all_variances <- sum(y^2)
  ## Y part
  n <- nrow(y)
  denom <- 1
  p <- ncol(x)
  q <- ncol(y)
  ## model ##

  m_h <- model_PLS(x,y,lam,tau=tau,method=method,NZV=NZV)
  no_model <- m_h$no_model ; U_out <- m_h$U_out ; V_out <- m_h$V_out ; bXr <- m_h$bXr ; bYr <- m_h$bYr
  residuals <- list(x=m_h$e_x,y=m_h$e_y)
  ########
  if(!no_model){
    pos_na <- which(is.na(V_out[1,]))
    if(length(pos_na)>0){
      R_END <- min(pos_na)-1
    }else{
      R_END <- R_max
    }
    ALPHAS=TERM_CROISE=NORM_TH <- rep(NA,R_END)
    for(R_end in 1:R_END){
      V <- V_out[,1:R_end,drop=F]
      U <- U_out[,1:R_end,drop=F]
      bXr_here <- bXr[1:R_end,,drop=F]
      bYr_here <- bYr[1:R_end,,drop=F]
      U_out_here <- U
      V_out_here <- V
      if(R_end>1){
        for(r_s in 2:R_end){
          id_m <- rev(1:(r_s-1))
          u_cur <- U[,r_s,drop=F]
          v_cur <- V[,r_s,drop=F]
          for(m in id_m){
            u_m <- U[,m]
            bX_m <- bXr_here[m,]
            u_cur <- u_cur - u_m*sum(bX_m*u_cur)
            v_m <- V[,m]
            bY_m <- bYr_here[m,]
            v_cur <- v_cur
          }
          U_out_here[,r_s] <- u_cur
          V_out_here[,r_s] <- v_cur
        }
      }
      B <- tcrossprod(U_out_here,V_out_here)
      ## YtX part
      COV <- crossprod(y,x)/n
      # COV_COV <- getCOV_COV(COV,lam)
      abs_COV <- abs(COV)
      max_cov_x <- apply(abs_COV,2,max)
      id_x_high <- which(max_cov_x>lam)
      max_cov_y <- apply(abs_COV,1,max)
      id_y_high <- which(max_cov_y>lam)
      COV_high <- COV[id_y_high,id_x_high,drop=F]
      abs_COV_high <- abs(COV_high)
      COV_COV_high <- COV_high - lam
      COV_COV_high[which(COV_COV_high<0)] <- 0
      # COV_COV_high <- crossprod(COV_high,COV_COV_high*sign(COV_high))
      if(method==2){
        COV_COV_high <- crossprod(COV_high,COV_COV_high*sign(COV_high))
        # COV_COV <- getCOV_COV(COV,lam)
      }else{
        COV_COV_high <- COV_COV_high*sign(COV_high)
      }

      varY <- all_variances

      #Build descriptor
      B_o <- B[id_x_high,,drop=F]
      Y_here <- y#y[,id_y_high,drop=F]
      y_est <- x%*%B
      y_here_y_est <- crossprod(Y_here,y_est)
      numer <- sum(y_here_y_est^2)
      denom <- sum(Y_here^2)^2
      ALPHAS[R_end] <- numer/denom
    }
    R_END <- which.max(ALPHAS)[1]
    alpha <- max(ALPHAS)[1]
    term_croise <- TERM_CROISE[R_END]
    norm_th <- NORM_TH[R_END]
    U_out <- U_out_here[,1:R_END,drop=F]
    V_out <- V_out_here[,1:R_END,drop=F]
    B <- tcrossprod(U_out,V_out)
  }else{
    R_END <- NULL
    U_out <- matrix(0,p,0)
    V_out <- matrix(0,q,0)
    term_croise <- NA
    norm_th <- NA
    alpha <- NA
    B <- NULL
  }
  list(alpha=alpha,
       norm_th=norm_th,
       term_croise=term_croise,
       U=U_out,V=V_out,B=B,R=R_END,residuals=list(x=residuals$x,y=residuals$y))
}

get_lambda_U_R <- function(x,y,lams,tau,R_max,method=2,NZV=1e-8){
  lam_opt <- lams[1]
  R <- 1
  alpha <- 1

  if(length(lam_opt)!=0){
    RES_0 <- deflat_in_loop(x,y,R_max,lam_opt,method=method,lam_min=lams[1],
                            tau=tau,NZV=1e-3)
  }else{
    RES_0 <- list(R=NULL,B=NULL,U=0,residuals=list(x=x,y=y))
  }
  list(lambda_hat=lam_opt,R_hat=RES_0$R,alpha_hat=alpha,R=R,
       lambdas=lams,U_hat=RES_0$U,B_hat=RES_0$B,
       RES=RES_0)
}

auto_loop_ddsPLS <- function(Xs,Y,
                             lambdas=NULL,lambdas_super=NULL,
                             tau=1e-2,
                             NZV=1e-3,method=2,
                             R_max=NA,
                             id_na=NULL,var_selected=NULL,
                             keep_model_per_block=T){
  if(is.na(R_max)){
    R_max = n <- nrow(Y)
  }else{
    n <- nrow(Y)
  }
  K <- length(Xs)
  q <- ncol(Y)
  ps <- unlist(lapply(Xs,ncol))
  u_k = t_k <- list()
  model_per_block = lams <- list()

  if(is.null(var_selected)){
    ## Get NA positions per block
    id_na <- lapply(Xs,function(xx){which(is.na(xx[,1]))})
  }else{
    for(k in 1:K){
      id_na_k <- id_na[[k]]
      if(length(id_na_k)!=0){
        Xs[[k]][id_na_k,] <- NA
      }
    }
  }
  for(k in 1:K){
    if(is.null(lambdas)){
      al_cors <- range(abs(cor(Xs[[k]],Y,use = "pairwise")))
      lams[[k]] <- seq(min(na.omit(lambda_min,al_cors[1])),al_cors[2],
                       length.out = n_lambdas)
    }else{
      lams[[k]] <- lambdas
    }
  }
  ## Imputation using, or not, the selected variables
  block_selected <- rep(0,K)
  for(k in 1:K){
    p_k <- ps[k]
    pos_na_k <- id_na[[k]]
    nb_na_k <- 0
    if(length(pos_na_k)>0){
      nb_na_k <- length(pos_na_k)
      y_train <- Xs[[k]][-pos_na_k,,drop=F]
      test_ok <- TRUE
      if(!is.null(var_selected)){
        var_k <- var_selected[[k]]
        if(length(var_k)>0){
          y_train <- y_train[,var_selected[[k]],drop=F]
        }else{
          test_ok <- FALSE
        }
      }
      p_k <- ncol(y_train)
      if(test_ok){
        mu_k_0 <- colMeans(y_train)
        sd_k_0 <- apply(y_train, 2, sd)
        x_train <- Y[-pos_na_k,,drop=F]
        mu_y_0 <- colMeans(x_train)
        sd_y_0 <- apply(x_train, 2, sd)
        x_test <- Y[pos_na_k,,drop=F]
        for(j in 1:q){
          if(sd_y_0[j]!=0){
            x_test[,j] <- (x_test[,j]-mu_y_0[j])/sd_y_0[j]
            x_train[,j] <- (x_train[,j]-mu_y_0[j])/sd_y_0[j]
          }else{
            x_test[,j] <- 0
            x_train[,j] <- 0
          }
        }
        for(i in 1:p_k){
          if(sd_k_0[i]!=0){
            y_train[,i] <- (y_train[,i]-mu_k_0[i])/sd_k_0[i]
          }else{
            y_train[,i] <- 0
          }
        }
        model_init <- get_lambda_U_R(
          x_train,y_train,lams=lams[[k]],tau=tau,
          R_max=min(R_max,ncol(x_train)),method=method)
        ##
        B_k <- model_init$B_hat
        y_test <- matrix(0,length(pos_na_k),p_k)
        if(!is.null(B_k)){
          y_test <- x_test%*%B_k
        }
        for(i in 1:p_k){
          y_test[,i] <- y_test[,i]*sd_k_0[i]+mu_k_0[i]
        }
      }
      MU_K_ALL <- colMeans(Xs[[k]][-pos_na_k,,drop=F])
      for(i in 1:ps[[k]]){
        Xs[[k]][pos_na_k,i] <- MU_K_ALL[i]
      }
      if(!is.null(var_selected)){
        var_k <- var_selected[[k]]
        if(length(var_k)>0){
          Xs[[k]][pos_na_k,var_k] <- y_test
        }
      }else{
        Xs[[k]][pos_na_k,] <- y_test
      }
    }
    ########
    Xs_k <- Xs[[k]]
    if(!is.null(var_selected)){
      var_k <- var_selected[[k]]
      if(length(var_k)>0){
        Xs_k <- Xs_k[,var_k,drop=F]
      }else{
        Xs_k <- matrix(0,nrow(Xs_k),1)
      }
    }
    model_k <- get_lambda_U_R(Xs_k,Y,lams=lams[[k]],method=method,tau=tau,
                              R_max=min(R_max,ncol(Xs_k)))
    if(!is.null(var_selected)){
      var_k <- var_selected[[k]]
      if(length(var_k)>0){
        R_k_hat <- model_k$R_hat
        if(R_k_hat!=0){
          U_k_hat <- matrix(0,ncol(Xs[[k]]),R_k_hat)
          U_k_hat[var_k,] <- model_k$U_hat
          model_k$U_hat =model_k$RES$U <- U_k_hat
          B_k_hat <- matrix(0,ncol(Xs[[k]]),ncol(Y))
          B_k_hat[var_k,] <- model_k$B_hat
          model_k$B_hat = model_k$RES$B <- B_k_hat
        }
      }
    }
    ########
    if(keep_model_per_block){
      model_per_block[[k]] <- model_k
    }
    if(!is.null(model_k$RES$R)){
      u_k[[k]] <- model_k$RES$U
      t_k[[k]] <- Xs[[k]]%*%u_k[[k]]
    }else{
      u_k[[k]] <- NA
      t_k[[k]] <- NULL
    }
  }
  ## Prepare output
  B_K_TOTAL <- model_per_block[[1]]$B_hat
  if(is.null(B_K_TOTAL)){
    SELECTED_VARS <- NULL
    B_K_TOTAL <- matrix(0,ncol(Xs[[1]]),q)
  }else{
    SELECTED_VARS <- which(rowSums(abs(B_K_TOTAL))>NZV)
  }
  out <- list(B=list(B_K_TOTAL),Xs=Xs,selected_x=SELECTED_VARS,
              u=u_k,model_per_block=model_per_block)
}



#' Automatic
#'
#' @param Xs list
#' @param Y matrix
#' @param lambda_min coef
#' @param n_lambdas number of coeffs
#' @param tau rate of learning
#' @param NZV near zero variance
#' @param keep_model_per_block must be put to true
#' @param plotVarB logical
#'
#' @return
#' @export
#'
#' @useDynLib ddsPLS
auto_ddsPLS <- function(Xs,Y,
                        lambdas=NULL,
                        tau=1e-2,R_max=NA,
                        lambdas_super=NULL,method=2,
                        NZV=1e-3,
                        scaleMat=T,
                        keep_model_per_block=T,
                        plotVarB=T){
  id_na <- lapply(Xs,function(xx){which(is.na(xx[,1]))})
  K <- length(id_na)
  n <- nrow(Xs[[1]])
  q <- ncol(Y)
  MU_X_K = SD_X_K <- list()
  for(k in 1:K){
    if(scaleMat){
      MU_X_K[[k]] <- colMeans(Xs[[k]],na.rm = T)
      SD_X_K[[k]] <- apply(Xs[[k]],2,sd,na.rm = T)
      Xs[[k]] <- scale(Xs[[k]])
    }else{
      p <- ncol(Xs[[k]])
      MU_X_K[[k]] <- rep(0,p)
      SD_X_K[[k]] <- rep(1,p)
    }
  }
  if(scaleMat){
    MU_Y <- colMeans(Y,na.rm = T)
    SD_Y <- apply(Y,2,sd,na.rm = T)
    Y <- scale(Y)
  }else{
    MU_Y <- rep(0,q)
    SD_Y <- rep(1,q)
  }
  Xs_init_scaled <- Xs
  model <- auto_loop_ddsPLS(
    Xs_init_scaled,Y,
    lambdas=lambdas,lambdas_super=lambdas_super,
    tau=tau,R_max=R_max,
    NZV=NZV,method=method,
    keep_model_per_block=keep_model_per_block)
  var_sel_0 <- model$selected_x
  model_0 <- model
  are_there_na <- sum(unlist(id_na))!=0
  if(length(var_sel_0)>0){
    iteration <- 1
    if(are_there_na){
      test <- TRUE # '== one_at_least_is_different'
      while(test){
        model <- auto_loop_ddsPLS(
          Xs_init_scaled,Y,
          lambdas=lambdas,lambdas_super=lambdas_super,
          tau=tau,R_max=R_max,method=method,NZV=NZV,
          id_na=id_na,var_selected=var_sel_0,
          keep_model_per_block=keep_model_per_block)
        var_sel <- model$selected_x
        if(length(var_sel)==0){
          test <- FALSE
        }else{
          for(k in 1:K){
            new_sel_k_is_empty <- FALSE
            test <- test &
              !setequal(var_sel[[k]],var_sel_0[[k]])
          }
        }
        var_sel_0 <- var_sel
        model_0 <- model
        iteration <- iteration + 1
      }
    }

  }else{
    iteration <- 0
  }
  Xs_out <- model$Xs
  for(k in 1:K){
    MU_k_rep <- matrix(rep(MU_X_K[[k]],n),nrow = n,byrow = T)
    SD_k_rep <- matrix(rep(SD_X_K[[k]],n),nrow = n,byrow = T)
    Xs_out[[k]] <- Xs_out[[k]]*SD_k_rep+MU_k_rep
  }
  model$Xs <- NULL
  parameters <- list(
    lambdas=lambdas,lambdas_super=lambdas_super,tau=tau,NZV=NZV,
    keep_model_per_block=keep_model_per_block
  )
  out <- list(model=model,nb_iterations=iteration,
              Xs=Xs_out,mu_X=MU_X_K,sd_X=SD_X_K,q=q,
              mu_Y=MU_Y,sd_Y=SD_Y,parameters=parameters)

  if(plotVarB){
    plot_var_B(out)
  }
  out
}

#' Title
#'
#' @param Xs Xs
#' @param Y Y
#' @param posLegend posLegend
#' @param x.intersp x.intersp
#' @param y.intersp y.intersp
#'
#' @return
#' @export
#'
#' @useDynLib ddsPLS
est_l_min <- function(Xs,Y,names_blocks=NULL,
                      posLegend = "topright",
                      x.intersp = 0.7,y.intersp = 0.7){
  K <- length(Xs)
  if(is.null(names_blocks)){
    names_blocks <- paste("Block",1:K)
  }
  if(K+1<=8){
    cols=RColorBrewer::brewer.pal(K+1,name = "Dark2")
  }else{
    if(K+1<=12){
      cols=RColorBrewer::brewer.pal(K+1,name = "Paired")
    }else{
      cols=1:(K+1)
    }
  }
  denx=deny <- list()
  hist_p <- list()
  for(k in 1:K){
    cor_k <- abs(cor(Xs[[k]],Y,use = "pairwise.complete"))
    hist_p[[k]] <- hist(cor_k,100,plot = F)
    # hist_p_k <- hist(cor_k,100,plot = F)
    nn <- length(hist_p[[k]]$breaks)
    Sum_count <- numeric(nn)
    for(i in 1:(nn-1) ){
      Sum_count[i] <- sum(hist_p[[k]]$counts[i:(nn-1)])
    }
    hist_p[[k]]$counts <- Sum_count
    # den <- density(cor_k)
    # denx[[k]] <- den$x
    # deny[[k]] <- den$y
  }
  # ylim <- range(c(deny))
  ylim <- range(unlist(lapply(hist_p,function(f){f$counts})))
  for(k in 1:K){
    if(k==1){
      main <- expression("Cumulative counts of the observed covariate coefficients")
      plot(hist_p[[k]]$breaks,
           hist_p[[k]]$counts,
           ylim=ylim,xlim=c(0,1),type="l",
           xlab=expression(lambda),ylab=expression(Count),
           col=cols[k],lwd=2,main=main)
      # plot(denx[[k]],deny[[k]],type="l",
      #      xlab=expression(lambda),ylab=expression(density),
      #      ylim=ylim,col=cols[k],lwd=2,main=main)
    }else{
      # points(denx[[k]],deny[[k]],type="l",col=cols[k],lwd=2)
      points(hist_p[[k]]$breaks,
             hist_p[[k]]$counts,type="l",col=cols[k],lwd=2)
    }
  }
  abline(h=0,lty=2)
  legend(posLegend,legend = names_blocks,col=cols,
         lty=1,x.intersp = x.intersp,y.intersp = y.intersp,
         bty="n",lwd=2)
  # lams <- seq(lambda_0,1,length.out = nlambda)
  # Var_B = Lambdas <- list()
  # Mo <- list()
  # Mo$model <- list()
  # Mo$model$model_per_block <- list()
  # lam_01 <- rep(NA,K)
  # low_lambda <- list()
  # for(k in 1:K){
  #   id_na_k <- which(is.na(Xs[[k]][,1]))
  #   if(length(id_na_k)>0){
  #     xk <- list(Xs[[k]][-id_na_k,])
  #     yk <- Y[-id_na_k,]
  #   }else{
  #     xk <- list(Xs[[k]])
  #     yk <- Y
  #   }
  #   Mo$model$model_per_block[[k]] <- auto_ddsPLS(
  #     xk,yk,tau=tau,lambdas = lams)$model$model_per_block[[1]]
  #   R_k <- Mo$model$model_per_block[[k]]$R
  #   lams_k <- Mo$model$model_per_block[[k]]$lambdas
  #   if(length(na.omit(R_k))>0){
  #     bad_val <- boxplot.stats(R_k)$out
  #     low_lambda[[k]] <- which(R_k %in% bad_val & R_k > median(R_k,na.rm = T))
  #   }
  # }
  # max_all <- lapply(low_lambda,function(aa){
  #   if(length(aa)>0) max(aa)})
  # lambda_max <- lapply(1:K,
  #                      function(k){
  #                        ls <-  Mo$model$model_per_block[[k]]$lambdas
  #                        ls[max_all[[k]]+1]
  #                      }
  # )
  # if(plotResults){
  #   plot_var_B(Mo,plotSuper = F,posLegend = posLegend,
  #              low_lambda=low_lambda)
  # }
  # max(unlist(lambda_max))
}

#' Title
#'
#' @param model model
#' @param newX new dataset
#'
#' @return
#' @export
#'
#' @useDynLib ddsPLS
predict_auto_ddsPLS <- function(model,newX){
  para <- model$parameters
  K <- length(newX)
  n <- nrow(model$Xs[[1]])
  n_new <- nrow(newX[[1]])
  q <- model$q
  id_na_new <- lapply(newX,function(xx){which(is.na(xx[,1]))})
  are_there_na <- sum(unlist(id_na_new))!=0
  Xs_in <- list()
  if(!are_there_na){
    newX_IN <- newX
    MU_super <- model$model$mu
    B_super <- model$model$B
    Y_est <- list()
    count <- 0
    Y_est <- matrix(0,n_new,q)
    for(k in 1:K){
      ## Center new individual
      MU_k <- matrix(rep(model$mu_X[[k]],n_new),nrow = n_new,byrow = T)
      SD_k <- matrix(rep(model$sd_X[[k]],n_new),nrow = n_new,byrow = T)
      newX_IN[[k]] <- (newX_IN[[k]]-MU_k)/SD_k
      pos_Inf <- which(newX_IN[[k]]==Inf)
      if(length(pos_Inf)>0){
        newX_IN[[k]][pos_Inf] <- 0
      }
      u_k <- model$model$u[[k]]
      if(!anyNA(u_k)){
        m_k <- count + 1:ncol(u_k)
        for(m in m_k){
          Y_est <- Y_est +
            matrix(rep(model$model$mu[m,],n_new),nrow = n_new,byrow = T)
        }
        Y_est <- Y_est + newX_IN[[k]]%*%model$model$B[[k]]
      }
    }
    Y_est <- Y_est*matrix(rep(model$sd_Y,n_new),ncol = q,byrow = T)+
      matrix(rep(model$mu_Y,n_new),ncol = q,byrow = T)
  }
  else{
    ### If there are missing values in the TEST
    newX_IN <- newX
    if(length(model$model$model_super)>0){
      if(model$model$model_super$R_hat!=0){
        id_k_m_r <- cbind(1,1:nrow(model$model$mu),1)
        u_super <- model$model$u_super
        count <- 0
        for(k in 1:K){
          if(count<nrow(id_k_m_r)){
            m <- nrow(u_super[[k]])
            if(length(m)>0){
              id_k_m_r[count+1:m,1] <- k
              id_k_m_r[count+1:m,3] <- 1:m
              count <- count + m
            }
          }
        }
        for(i in 1:n_new){
          x_i <- lapply(newX_IN,function(xx){xx[i,,drop=F]})
          pos_na <- which(unlist(lapply(x_i,function(xx){is.na(xx[1])})))
          if(length(pos_na)>0){
            for(k in pos_na){
              newX_IN[[k]][i,] <- model$mu_X[[k]]
            }
            x_i_to_comp <- x_i[pos_na]
            m_mis <- id_k_m_r[which(id_k_m_r[,1]%in%pos_na),2]
            m_not_mis <- (1:nrow(id_k_m_r))[-m_mis]
            if(length(m_mis)>0 & length(m_not_mis)!=0){
              t_i_not_mis <- matrix(0,n,nrow(id_k_m_r)-length(m_mis))
              t_i_mis <- matrix(0,1,nrow(id_k_m_r)-length(m_mis))
              for(jj in 1:length(m_not_mis)){
                m_cur <- m_not_mis[jj]
                k_cur <- id_k_m_r[m_cur,1]
                r_cur <- id_k_m_r[m_cur,3]
                u_cur <- model$model$model_per_block[[k_cur]]$U_hat[,r_cur]
                X_cur <- model$Xs[[k_cur]]
                MU_notmiss <- matrix(rep(model$mu_X[[k_cur]],n),nrow = n,byrow = T)
                SD_notmiss <- matrix(rep(model$sd_X[[k_cur]],n),nrow = n,byrow = T)
                t_i_not_mis_jj <- (X_cur-MU_notmiss)/SD_notmiss
                t_i_mis_jj <- (x_i[[k_cur]]-model$mu_X[[k_cur]])/model$sd_X[[k_cur]]
                pos_inf <- which(model$sd_X[[k_cur]]<para$NZV)
                if(length(pos_inf)>0){
                  for(jjj in pos_inf){
                    t_i_mis_jj[,jjj] <- 0
                  }
                }
                t_i_not_mis[,jj] <- t_i_not_mis_jj%*%u_cur
                t_i_mis[,jj] <- t_i_mis_jj%*%u_cur
              }
              k_mis <- unique(id_k_m_r[m_mis,1])
              X_mis <- list()
              ps_mis <- NULL
              for(id_k_mis in 1:length(k_mis)){
                k_mis_j <- k_mis[id_k_mis]
                var_k_mis_j <- model$model$selected_x[[k_mis_j]]
                if(length(var_k_mis_j)>0){
                  X_mis[[id_k_mis]] <- model$Xs[[k_mis_j]][,var_k_mis_j,drop=F]
                  ps_mis <- c(ps_mis,length(var_k_mis_j))
                }else{
                  X_mis[[id_k_mis]] <- NULL
                  ps_mis <- c(ps_mis,0)
                }
              }
              X_mis <- do.call(cbind,X_mis)
              if(!is.null(X_mis)){
                MU_X_mis <- colMeans(X_mis)
                SD_X_mis <- apply(X_mis,2,sd)
                X_mis <- scale(X_mis)
                model_imput <- auto_ddsPLS(
                  Xs = list(t_i_not_mis),Y = X_mis,lambdas=para$lambdas,
                  tau=para$tau,method=method,
                  NZV=para$NZV,
                  keep_model_per_block=para$keep_model_per_block)
                X_predicted <- predict_auto_ddsPLS(model_imput,list(t_i_mis))
                X_predicted <- X_predicted*SD_X_mis+MU_X_mis
                count <- 0
                for(id_k_mis in 1:length(k_mis)){
                  k_mis_j <- k_mis[id_k_mis]
                  if(length(model$model$selected_x)>0){
                    var_k_mis_j <- model$model$selected_x[[k_mis_j]]#model_imput$model$selected_x[[1]]#
                    p_mis_k <- ps_mis[id_k_mis]#length(var_k_mis_j)#
                    if(p_mis_k>0 & norm(model_imput$model$B[[1]])>0){
                      newX_IN[[k_mis_j]][i,var_k_mis_j] <-
                        X_predicted[,count + 1:p_mis_k]
                    }
                    count <- count + p_mis_k
                  }
                }
              }
            }
            else{
              for(k in pos_na){
                newX_IN[[k]][i,] <- model$mu_X[[k]]
              }
            }
          }
        }
        Y_est <- predict_auto_ddsPLS(model,newX_IN)

      }else{
        Y_est <- matrix(rep(model$mu_Y,n_new),ncol = q,byrow = T)
      }
    }else{
      Y_est <- matrix(rep(model$mu_Y,n_new),ncol = q,byrow = T)
    }
  }
  Y_est
}


#' Title
#'
#' @param Xs Xs
#' @param Y Y
#' @param n_rate number of rates
#' @param rate_min min rate
#' @param rate_max max rate
#' @param tau rate learning
#' @param R_max max
#' @param lambda_min min lambda
#' @param n_lambdas number of lambdas
#' @param kfolds number of folds
#' @param fold_fixed yes or no
#' @param NCORES number of cores
#' @param NZV near zero var
#'
#' @return
#' @export
#'
#' @useDynLib ddsPLS
perf_auto_ddsPLS <- function(Xs,Y,
                             lambdas = 0.5,lambdas_super=NULL,
                             tau=1e-2,
                             kfolds="loo",method=2,
                             fold_fixed=NULL,NCORES=1,
                             NZV=1e-3){#,
  # threshs_to_test=NULL){
  n <- nrow(Xs[[1]])
  q <- ncol(Y)
  ## CV design
  if(kfolds=="loo"){
    kfolds <- n
    fold <- 1:n
  }else if(!is.null(fold_fixed)){
    fold <- fold_fixed
  }else{
    fold <- replicate(n/kfolds+1,sample(1:kfolds))[1:n]
  }

  paras <- expand.grid(1:max(fold),tau)

  if(NCORES>nrow(paras)){
    decoupe <- 1:nrow(paras)
  }else{
    decoupe <- replicate(nrow(paras)/NCORES + 1,
                         sample(1:NCORES))[1:nrow(paras)]
  }
  NCORES_w <- min(NCORES,nrow(paras))
  `%my_do%` <- ifelse(NCORES_w!=1,{
    out<-`%dopar%`
    cl <- makeCluster(NCORES_w)#cl <- parallel::makeCluster(NCORES_w)
    registerDoParallel(cl)#doParallel::registerDoParallel(cl)
    out},{
      out <- `%do%`
      out})
  pos_decoupe <- 1
  ERRORS <- foreach(pos_decoupe=1:min(NCORES,nrow(paras)),
                    .combine = rbind,.packages = "Rcpp") %my_do% {
                      source('~/Documents/Hadrien/GitHub/ddsPLS/R/auto_ddsPLS.R')
                      Rcpp::sourceCpp('src/basics.cpp')
                      n <- nrow(Xs[[1]])
                      q <- ncol(Y)
                      K <- length(Xs)
                      paras_here_pos <- which(decoupe==pos_decoupe)
                      paras_here <- paras[paras_here_pos,,drop=FALSE]
                      errors <- matrix(NA,nrow(paras_here),q)
                      R_est=lambda_est <- matrix(NA,nrow(paras_here),K+1)
                      for(i in 1:nrow(paras_here)){
                        i_fold <- paras_here[i,1]
                        tau_i <- paras_here[i,2]
                        pos_train <- which(fold!=i_fold)
                        X_train <- lapply(Xs,function(xx){
                          xx[pos_train,,drop=F]})
                        X_test <- lapply(Xs,function(xx){
                          xx[-pos_train,,drop=F]})
                        Y_train <- Y[pos_train,,drop=FALSE]
                        Y_test <- Y[-pos_train,,drop=FALSE]
                        model_cv <- auto_ddsPLS(
                          X_train,Y_train,
                          tau=tau_i,
                          lambdas=lambdas,method=method,lambdas_super = lambdas_super)
                        for(k in 1:K){
                          R_k <- model_cv$model$model_per_block[[k]]$R_hat
                          lambda_k <- model_cv$model$model_per_block[[k]]$lambda_hat
                          if(!is.null(R_k)){
                            R_est[i,k] <- R_k
                            lambda_est[i,k] <- lambda_k
                          }
                        }
                        R_k <- model_cv$model$model_super$R_hat
                        lambda_k <- model_cv$model$model_super$lambda_hat
                        if(!is.null(R_k)){
                          R_est[i,K+1] <- R_k
                          lambda_est[i,K+1] <- lambda_k
                        }
                        Y_est <- predict_auto_ddsPLS(model_cv,X_test)
                        errors[i,] <- colMeans((Y_test-Y_est)^2)
                      }
                      out <- cbind(errors,R_est,lambda_est,paras_here)
                    }
  if(NCORES_w!=1){
    stopCluster(cl)
  }
  K <- length(Xs)
  MSE <- matrix(NA,length(tau),q)
  for(i_tau in 1:length(tau)){
    tau_i <- tau[i_tau]
    pos_cv <- which(ERRORS$Var2==tau_i)
    MSE[i_tau,] <- colMeans(ERRORS[pos_cv,1:q,drop=F])
  }
  if(is.null(colnames(Y))){
    colnames(MSE) <- paste("Y",1:q)
  }else{
    colnames(MSE) <- colnames(Y)
  }

  # MSE <- ERRORS[,1:q,drop=F]
  # if(is.null(colnames(Y))){
  #   colnames(MSE) <- paste("Y",1:q)
  # }else{
  #   colnames(MSE) <- colnames(Y)
  # }
  # R <- ERRORS[,q+1:(K+1),drop=F]
  # lambda <- ERRORS[,q+K+1+1:(K+1),drop=F]
  # colnames(R)=colnames(lambda) <- c(paste("Block",1:K),"Super block")
  # list(MSE=MSE,R=R,lambda=lambda)
  list(MSE=MSE,tau=tau,source=ERRORS)
}
