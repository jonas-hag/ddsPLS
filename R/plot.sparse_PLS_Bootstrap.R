#' Function to plot bootstrap performance results of the ddsPLS algorithm.
#'
#' @param x A sparse_PLS_Bootstrap object
#' @param h_opt Optimal number of components to be plotted
#' @param type Type of the weights and regression coefficients
#' @param lty Shape of the lines, if type=="l"
#' @param no_plot Internal logical parameter
#' @param ... Other plotting parameters to affect the plot.
#'
#'
#' @importFrom graphics layout boxplot
#'
#' @export
#'
#' @useDynLib ddsPLS
plot.sparse_PLS_Bootstrap <- function(x,h_opt=NULL,type="p",
                                      lty=NA,no_plot=F,...){
  if(!no_plot){
    ## Reset personnal plot par() settings
    opar <- par(no.readonly =TRUE)
    on.exit(par(opar))
    ## -----------------------------------
  }
  if(is.null(h_opt)){
    h_opt <- x$optimal_parameters$R
  }
  if(length(h_opt)>0){
    lambdas_out <- x$bootstrap$lambdas_out
    cols <- c(RColorBrewer::brewer.pal(max(h_opt,3),"Set1")[1:h_opt],"gray80")
    col_B <- RColorBrewer::brewer.pal(nrow(x$V),"Paired")
    layout(matrix(c(1,1,2,3,3,4,rep(10,3),
                    5,5,6,7,7,8,rep(10,3),
                    rep(9,9),
                    rep(11,9)), nrow=4, byrow = TRUE))
    # layout(matrix(c(1,1,3,2,2,4,4,5,5,6), 2, 5, byrow = TRUE))
    # layout(matrix(c(1,1,2,2,3,7,7,4,4,5,5,6,7,7), 2, 7, byrow = TRUE))
    par(mar=c(3,3,2,1),mgp=c(2,1,0))
    aa <- min(1/4,1-1/4)
    quart <- x$bootstrap$quartiles
    ls <- x$bootstrap$lambdas_out[[1]]
    plots <- list(
      R2_h=list(q50=quart$vars_boot_h_50,q75=quart$vars_boot_h_75,q25=quart$vars_boot_h_25,
                main=expression("R"["B,h"]^"2"),
                lam_opt=x$optimal_parameters$lambdas,
                vars_single=x$bootstrap$vars_h_boot_single),
      R2=list(q50=quart$vars_boot_50,q75=quart$vars_boot_75,q25=quart$vars_boot_25,
              main=expression("R"["B"]^"2"),
              lam_opt=x$optimal_parameters$lambdas,
              vars_single=x$bootstrap$vars_boot_single),
      Q2_h=list(q50=quart$q2_boot_50,q75=quart$q2_boot_75,q25=quart$q2_boot_25,
                main=expression("Q"["B,h"]^"2"),
                lam_opt=x$optimal_parameters$lambdas,
                vars_single=x$bootstrap$Q2_h_star),
      Q2=list(q50=quart$q2_all_boot_50,q75=quart$q2_all_boot_75,q25=quart$q2_all_boot_25,
              main=expression("Q"["B"]^"2"),
              lam_opt=x$optimal_parameters$lambdas,
              vars_single=x$bootstrap$Q2_all_sum_star)
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
        id_h <- x$id_ALL_TEST_h[[h]]
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
          paste("Component ",h," (",round(x$explained_variance[h]),"%)",sep="",collapse = "")})),
          "Not selected component"))
      }
      abline(v=x$optimal_parameters$lambdas,col=cols,lty=3)
      abline(v=x$optimal_parameters$lambdas,col=cols,lty=3)
      if(i==1){
        ddff <- data.frame(do.call(rbind,lapply(1:h_opt,function(ii,bb){cbind(ii,bb[[ii]])},x$bootstrap$R2_h_boxplot)))
        names(ddff) <- c("h","R2")
        bobo <- boxplot(R2~h,ddff,ylim=c(-0.1,1.15),border=cols,main=expression("R"["B,h"]^"2"),xlab="Component",ylab="")
        points(1:h_opt,x$explained_variance[1:h_opt]/100,col=1,pch=18,cex=2)
        abline(h=0,lty=5,col=1)
      }else if(i==2){
        ddff <- data.frame(do.call(rbind,lapply(1:h_opt,function(ii,bb){cbind(ii,bb[[ii]])},x$bootstrap$R2_all_boxplot )))
        names(ddff) <- c("h","R2")
        bobo <- boxplot(R2~h,ddff,ylim=c(-0.1,1.15),border=cols,main=expression("R"["B"]^"2"),xlab="Component",ylab="")
        abline(h=0,lty=5,col=1)
      }else if(i==3){
        ddff <- data.frame(do.call(rbind,lapply(1:h_opt,function(ii,bb){cbind(ii,bb[[ii]])},x$bootstrap$Q2_h_boxplot )))
        names(ddff) <- c("h","Q2")
        bobo <- boxplot(Q2~h,ddff,ylim=c(-0.1,1.15),border=cols,main=expression("Q"["B,h"]^"2"),xlab="Component",ylab="")
        abline(h=0,lty=5,col=1)
      }else{
        ddff <- data.frame(do.call(rbind,lapply(1:h_opt,function(ii,bb){cbind(ii,bb[[ii]])},x$bootstrap$Q2_all_boxplot )))
        names(ddff) <- c("h","Q2")
        bobo <- boxplot(Q2~h,ddff,ylim=c(-0.1,1.15),border=cols,main=expression("Q"["B"]^"2"),xlab="Component",ylab="")
        abline(h=0,lty=5,col=1)
      }
    }

    id_rev <- rev(1:h_opt)
    uu<-lapply(x$Us,function(u){u[,id_rev,drop=F]})
    uu <- do.call(rbind,uu)
    uu[which(uu==0)] <- NA
    matplot(uu,pch=id_rev,col="white",xlab="Index",ylab="Weight",main="Weights")
    abline(h=0,col="gray80")
    if(is.na(lty)){
      matplot(uu,pch=id_rev,col=cols[id_rev],add=T,type=type)
    }else{
      matplot(uu,pch=id_rev,col=cols[id_rev],add=T,type=type,lty=lty)
    }
    ps <- c(0,unlist(lapply(x$Us,nrow )))
    abline(v=cumsum(ps),lty=5)
    ############
    f <- do.call(cbind,lapply(
      1:h_opt,
      function(hh){
        id_h <- x$id_ALL_TEST_h[[hh]]
        out <- (x$bootstrap$vars_h_boot[[hh]]-x$bootstrap$Q2_all_sum_star[[hh]])#abs(x$bootstrap$vars_h_boot[[hh]]-x$bootstrap$Q2_all_sum_star[[hh]])
        out[-id_h] <- NA
        out
      }))
    matplot(ls,f,pch=1:h_opt,col=cols,xlab="Parameter",ylab="Weight",
            main=expression(bar("R"["B"]^"2")~"-"~bar("Q"["B"]^"2")))
    abline(v=x$optimal_parameters$lambdas,col=cols,lty=3)
    #########
    bb <- x$B_cbind
    id_b <- 1:ncol(bb)
    bb[which(bb==0)] <- NA
    matplot(bb,col="white",xlab="Index",ylab="Value",main="Regression coefficients (B)",pch=id_b)
    abline(h=0,col="gray80")
    if(is.na(lty)){
      matplot(bb,col=col_B,add=T,type=type,pch=id_b)
    }else{
      matplot(bb,col=col_B,add=T,type=type,pch=id_b,lty=lty)
    }
    abline(v=cumsum(ps),lty=5)
  }
}
