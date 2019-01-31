#' Function to plot \emph{mddsPLS}
#'
#' That function must be applied to a \emph{mddsPLS} object. Extra parameters are
#'  avalaible to control the plot quality.
#'
#' @param x The perf_mddsPLS object.
#' @param weights logical. Whether to plot the weight values. If \emph{FALSE} then the variates are plotted. Initialized to \emph{TRUE}.
#' @param super logical. If \emph{TRUE} barplots are filled with **Super-Weights** in the case of **weights** of with général **X** and **Y** components else.
#' @param block vector of intergers. If equals \emph{NULL} then all the components are plotted. Which components must be plotted. Initialized to \emph{NULL}.
#' @param mar_left positive float. Extra lines to add to the left margins, where the variable names are written.
#' @param pos_legend Initialized to "topright"
#' @param legend_names vector of character. Indicates the names of the blocks. Initialized to NULL and in this case just gets positions in the Xs list.
#' @param ... Other plotting parameters to affect the plot.
#'
#' @return The plot visualisation
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
#' x <- mddsPLS(Xs = X,Y = Y,R = 1, mode = "clas",lambda=0.8)
#' plot(x)
#'
#' # Regression example :
#' data("liver.toxicity")
#' X <- scale(liver.toxicity$gene)
#' Y <- scale(liver.toxicity$clinic)
#' #res_cv_reg <- perf_mddsPLS(Xs = X,Y = Y,lambda_min=0.8,n_lambda=2,R = 1,
#' # mode = "reg")
#' #plot(res_cv_reg)
plot.mddsPLS <- function(x,weights=TRUE,super=FALSE,block=NULL,mar_left=2,
                         pos_legend="topright",legend_names=NULL,
                         ...){
  browser()
  R <- x$mod$R
  K <- length(x$Xs)
  if(is.null(block)){
    block <- 1:K
  }
  if(is.null(legend_names) & is.null(names(x$Xs))){
    legend_names <- paste("Block",block,sep=" ")
  }
  if(is.null(block)){
    block <- 1:K
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
    if(!super){
      par(mfrow=c(length(block),R+1),
          mar=c(5,4+mar_left,4,2)+0.1)
    }else{
      par(mfrow=c(1,R+1),
          mar=c(5,4+mar_left,4,2)+0.1)
    }
    for(i_k in 1:length(block)){
      k <- block[i_k]
      toplot[[k]] <- list()
      viz_k <- viz[[k]]
      for(r in 1:R){
        viz_k_r <- viz_k[,r]
        pos_no_nul <- which(abs(viz_k_r)>1e-12)
        if(length(pos_no_nul)>0){
          toplot[[k]][[r]] <- viz_k[pos_no_nul,r]
          names(toplot[[k]][[r]]) <- colnames(x$Xs[[k]])[pos_no_nul]
          toplot[[k]][[r]] <- toplot[[k]][[r]][order(abs(toplot[[k]][[r]]),
                                                     decreasing = T)]

          if(!super){
            barplot(toplot[[k]][[r]],horiz = T,las=2,col=colors[k],xlim = c(-1,1),
                    main=paste(legend_names[i_k],", component ",R,sep=""))
          }
        }
        y_como <- viz_y[,r]
        names(y_como) <- colnames(x$mod$)[pos_no_nul]
        toplot[[k]][[r]] <- toplot[[k]][[r]][order(abs(toplot[[k]][[r]]),
                                                   decreasing = T)]
        barplot(toplot[[k]][[r]],horiz = T,las=2,col=colors[k],xlim = c(-1,1),
                main=paste(legend_names[i_k],", component ",R,sep=""))
      }
    }
    if(super){
      for(r in 1:R){
        plotR <- NULL
        cols <- NULL
        for(k in block){
          plotR <- c(plotR,toplot[[k]][[r]])
          cols <- c(cols,rep(colors[k],length(toplot[[k]][[r]])))
        }
        oo <- order(abs(plotR),decreasing = T)
        barplot(plotR[oo],horiz = T,las=2,col=cols[oo])
      }
    }
    legend(pos_legend,legend = legend_names,fill = colors[block])
  }else{
    viz <- x$mod$ts
  }


}
