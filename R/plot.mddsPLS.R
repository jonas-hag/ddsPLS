#' Function to plot \emph{mddsPLS}
#'
#' That function must be applied to a \emph{mddsPLS} object. Extra parameters are
#'  avalaible to control the plot quality.
#'
#' @param x The perf_mddsPLS object.
#' @param weights logical. Whether to plot the weight values. If \emph{FALSE} then the variates are plotted. Initialized to \emph{TRUE}.
#' @param super logical. If \emph{TRUE} barplots are filled with **Super-Weights** in the case of **weights** of with général **X** and **Y** components else.
#' @param block vector of intergers indicating which components must be plotted. If equals \emph{NULL} then all the components are plotted. Initialized to \emph{NULL}.
#' @param comp vector of intergers indicating which blocks must be plotted. If equals \emph{NULL} then all the blocks are plotted. Initialized to \emph{NULL}.
#' @param addY logical. Whether or not to plot **Block Y**. Initialized to \emph{FALSE}.
#' @param mar_left positive float. Extra lines to add to the left margins, where the variable names are written.
#' @param pos_legend Initialized to "topright"
#' @param legend_names vector of character. Indicates the names of the blocks. Initialized to NULL and in this case just gets positions in the Xs list.
#' @param ... Other plotting parameters to affect the plot.
#'
#' @return The plot visualisation
#'
#' @importFrom graphics abline arrows barplot legend par
#'
#' @seealso  \code{\link{mddsPLS}}, \code{\link{summary.mddsPLS}}
#'
#' @export
#'
#' @examples
#' library(doParallel)
#' # Classification example :
#' data("penicilliumYES")
#' X <- penicilliumYES$X
#' X <- scale(X[,which(apply(X,2,sd)>0)])
#' Y <- as.factor(unlist(lapply(c("Melanoconidiu","Polonicum","Venetum"),
#' function(tt){rep(tt,12)})))
#' x <- mddsPLS(Xs = X,Y = Y,R = 3, mode = "clas",lambda=0.8)
#' plot(x)
#'
#' # Regression example :
#' data("liver.toxicity")
#' X <- scale(liver.toxicity$gene)
#' Y <- scale(liver.toxicity$clinic)
#' #res_cv_reg <- ddsPLS(Xs = X,Y = Y,lambda=0.8,R = 2)
#' #plot(res_cv_reg)
plot.mddsPLS <- function(x,weights=TRUE,super=FALSE,mar_left=2,
                         block=NULL,comp=NULL,addY=FALSE,
                         pos_legend="topright",legend_names=NULL,
                         ...){
  R <- x$mod$R
  if(is.null(comp)){
    comp_in <- 1:R
  }else{
    comp_in <- comp
  }
  R_in <- length(comp_in)
  K <- length(x$Xs)
  block_in <- block
  if(is.null(block_in)){
    block_in <- 1:K
  }
  if (any(!comp_in %in% 1:R)) {
    stop("One asked component does not exist",
         call. = FALSE)
  }else if(any(!block_in %in% 1:K)){
    stop("One asked block does not exist",
         call. = FALSE)
  }
  isReg <- x$mode=="reg"
  Y_in <- x$Y_0
  if(!isReg){
    names_Y <- levels(Y_in)
    q <- nlevels(Y_in)
  }else{
    if(!(is.matrix(Y_in)|is.data.frame(Y_in))){
      Y_in <- as.matrix(Y_in)
    }
    if(is.data.frame(Y_in)){
      Y_in <- as.matrix(Y_in)
    }
    q <- ncol(Y_in)
    names_Y <- colnames(Y_in)
    if(is.null(names_Y)){
      names_Y <- 1:q
    }
  }
  if(is.null(legend_names) & is.null(names(x$Xs))){
    legend_names_in <- paste("Block",block_in,sep=" ")
  }
  l_bl <- K+1
  if(l_bl<3){
    colors <- 1:l_bl
  }else if(l_bl>8){
    pal <- grDevices::colorRampPalette(colors)
    colors <- pal(l_bl)
  }else{
    colors <- RColorBrewer::brewer.pal(l_bl, "Dark2")
  }
  if(weights){
    viz <- x$mod$u
    viz_y <- x$mod$v
    if(super){
      viz <- x$mod$u_t_super
    }
    toplot <- list()
    ind_1 <- 1:length(block_in)
    ind_2 <- 1:R_in
    if(!super){
      if(addY){
        par(mfrow=c(R_in,length(block_in)+1),
            mar=c(5,4+mar_left,4,2)+0.1)
      }else{
        par(mfrow=c(R_in,length(block_in)),
            mar=c(5,4+mar_left,4,2)+0.1)
      }
    }else{
      if(addY){
        par(mfrow=c(R_in,1+1),
            mar=c(5,4+mar_left,4,2)+0.1)
      }else{
        par(mfrow=c(R_in,1),
            mar=c(5,4+mar_left,4,2)+0.1)
      }
    }
    for(i_r in ind_2){
      for(i_k in ind_1){
        r <- comp_in[i_r]
        k <- block_in[i_k]
        if(i_r==1){
          toplot[[k]] <- list()
        }
        viz_k <- viz[[k]]
        viz_k_r <- viz_k[,r]
        pos_no_nul <- which(abs(viz_k_r)>1e-12)
        main <- paste(legend_names_in[i_k],", component ",r,sep="")
        if(length(pos_no_nul)>0){
          toplot[[k]][[r]] <- viz_k[pos_no_nul,r]
          names(toplot[[k]][[r]]) <- colnames(x$Xs[[k]])[pos_no_nul]
          toplot[[k]][[r]] <- toplot[[k]][[r]][order(abs(toplot[[k]][[r]]),
                                                     decreasing = T)]

          if(!super){
            barplot(toplot[[k]][[r]],horiz = T,las=2,col=colors[k],xlim = c(-1,1),
                    main=main)
          }
        }else{
          toplot[[k]][[r]] <- 0
          if(!super){
            plot(1, type="n", axes=F, xlab="", ylab="",main=main)
          }
        }
      }
      if(addY & !super){
        y_como <- viz_y[,r]
        pos_no_nul <- which(abs(y_como)>1e-12)
        if(length(pos_no_nul)>0){
          y_como <- y_como[pos_no_nul]
          names(y_como) <- names_Y[pos_no_nul]
        }
        toplot_y <- y_como[order(abs(y_como),decreasing = T)]
        barplot(toplot_y,horiz = T,las=2,col=colors[K+1],xlim = c(-1,1),
                main=paste("Bloc Y, component ",r,sep=""))
        legeds <- c(legend_names_in,"Block Y")
        colOut <- colors[c(block_in,K+1)]
      }else{
        legeds <- legend_names_in
        colOut <- colors[block_in]
      }
    }
    if(super){
      for(i_r in 1:R_in){
        r <- comp_in[i_r]
        plotR <- NULL
        cols <- NULL
        for(k in block_in){
          if(toplot[[k]][[r]][1]!=0){
            plotR <- c(plotR,toplot[[k]][[r]])
            cols <- c(cols,rep(colors[k],length(toplot[[k]][[r]])))
          }
        }
        main <- paste("Bloc Xs, component ",r,sep="")
        if(is.null(plotR)){
          plot(1, type="n", axes=F, xlab="", ylab="",main=main)
        }else{
          oo <- order(abs(plotR),decreasing = T)
          barplot(plotR[oo],horiz = T,las=2,col=cols[oo],main=main)
        }
        if(addY){
          y_como <- viz_y[,r]
          pos_no_nul <- which(abs(y_como)>1e-12)
          if(length(pos_no_nul)>0){
            y_como <- y_como[pos_no_nul]
            names(y_como) <- names_Y[pos_no_nul]
          }
          toplot_y <- y_como[order(abs(y_como),decreasing = T)]
          barplot(toplot_y,horiz = T,las=2,col=colors[K+1],xlim = c(-1,1),
                  main=paste("Bloc Y, component ",r,sep=""))
          legeds <- c(legend_names_in,"Block Y")
          colOut <- colors[c(block_in,K+1)]
        }else{
          legeds <- legend_names_in
          colOut <- colors[block_in]
        }
      }
    }
    legend(pos_legend,legend = legeds,fill = colOut)
  }else{
    viz <- x$mod$ts
  }
}
