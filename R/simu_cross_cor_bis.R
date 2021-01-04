give_me_plot_comm_bis <- function(){
  par(mar=c(3,3,3,2))
  # layout(matrix(c(1,2,
  #                 1,2,
  #                 3,4,
  #                 3,5), 4, byrow = TRUE),mar=c(3,3,3,2))
  # layout(matrix(c(1,2), nrow = 1, byrow = TRUE))
  layout(matrix(c(1,2,3,4),ncol=2,byrow = T))
  # par(mfrow=c(1,3),mar=c(3,3,3,2))
  cols <- RColorBrewer::brewer.pal(length(method)+1,"Set1")[-6]
  col_box <- RColorBrewer::brewer.pal(9,"Pastel1")[7]

  ylim <- c(-0.05,1.1*(2*eps^2)/3)
  lwd <- 1.5
  ncols <- 3
  cex.leg <- 1
  boxplot(Q2~method*rho,border=cols,col=col_box,df,main=expression("Q"^"2"),#border=cols[rep(1:l_m,i)]
          xlab="",xaxt="n",ylab="",lwd=lwd);#abline(h=c(0,0.0975),lty=3,lwd=2)
  abline(v=0.5+c(0:1 )*length(levels(df$method)),lty=3)
  abline(h=0,lty=2,col="black")
  legend("bottomright",legend = level_legend,fill = cols,
         ncol = 1,bg = "white",cex = cex.leg,bty="n")

  boxplot(time~method*rho,border=cols,col=col_box,df,main="Raw time (not adjusted on number of parameters)",
          xlab="",xaxt="n",ylab="",lwd=lwd);#abline(h=c(0,0.0975),lty=3,lwd=2)
  abline(v=0.5+c(0:(length(n_rho_s)+1) )*length(levels(df$method)),lty=3)
  abline(h=0,lty=2)
  legend("topright",legend = level_legend,fill = cols,
         ncol = 1,bg = "white",cex = cex.leg,bty="n")

  boxplot(1-DIST_B_hat~method*rho,df,#[-which(df$method=="OLS_esp"),],
          main=expression("||A"~hat(B)~"-C||"^2/"||C||"^2),
          col=col_box,border=cols,xlab="",xaxt="n",ylab="",lwd=lwd)#,ylim=c(0,1))
  abline(v=0.5+c(0:(length(n_rho_s)+1) )*length(unique(df$method)),lty=3)
  abline(h=0,lty=2)
  legend("topright",legend = level_legend,fill = cols,ncol = 1,bg = "white",cex = cex.leg,bty="n")

  boxplot(NORM_B_hat~method*rho,df,#[-which(df$method=="OLS_esp"),],
          main=expression("||"~hat(B)~"||"^2),
          col=col_box,border=cols,xlab="",xaxt="n",ylab="",lwd=lwd)#,ylim=c(0,1))
  abline(h=0,lty=2)
  abline(v=0.5+c(0:(length(n_rho_s)+1) )*length(unique(df$method)),lty=3)
  legend("topright",legend = level_legend,fill = cols,ncol = 1,bg = "white",cex = cex.leg,bty="n")

}

plot_R_rho <- function(){
  rho_s <- sort(unique(df$rho))
  n_rho_s <- length(rho_s)
  labs <- lapply(rho_s,function(rr){bquote(rho==.(rr))})
  n_meth <- length(method)
  layout(matrix(1:(n_rho_s*n_meth),nrow = n_rho_s,byrow = T))
  par(mar=c(2,2,2,3))
  cols <- RColorBrewer::brewer.pal(length(method)+1,"Set1")[-6]
  breaks <- (min(c(2,na.omit(df$R)))-0.5):(max(c(2,na.omit(df$R)))+0.5)
  for(i_n in 1:n_rho_s){
    for(i_mm in 1:n_meth){
      Ri <- df$R[which(df$method==method[i_mm] & df$rho==rho_s[i_n])]
      # Ri[which(is.na(Ri))] <- 0

      hist(rep(3,100),breaks = breaks,xaxt="n",yaxt="n",col="gray90",ylim=c(0,110),
           main="",xlab="",ylab="")
      axis(side = 1,unique(na.omit(df$R)),line = -1,labels = unique(na.omit(df$R)),tick = F)
      axis(side = 2,(0:5)*20,line = -1,labels = (0:5)*20,tick = F)
      abline(h=(0:5)*20,lty=2,col="gray")
      hist(Ri,breaks = breaks,col=cols[i_mm],add=T)
      title(bquote(.(method[i_mm])~"\n "~rho==.(rho_s[i_n])),line = -1)
      title(xlab=expression(hat(R)),ylab="Frequency",line = 1)
      xx <- unique(na.omit(df$R))
      yy <- unlist(lapply(xx,function(xxx){length(which(Ri==xxx))}))
      if(length(which(yy!=0))>0){
        text(xx[which(yy!=0)],yy[which(yy!=0)],labels = yy[which(yy!=0)],pos = 3,col=cols[i_mm])
      }
    }
  }
}

plot_TVP_TFP_X_rho <- function(){
  rho_s <- sort(unique(df$rho))
  n_rho_s <- length(rho_s)
  labs <- lapply(rho_s,function(rr){bquote(rho==.(rr))})
  p <- ncol(A)
  var_star <- which(colSums(abs(A))>0)
  p_star <- length(var_star)
  sel_x <- df[,c(1,3,7,8)]
  TVP <- (sel_x$VP_X)/(p_star)
  TFP <- (sel_x$SEL_X- sel_x$VP_X)/(p-p_star)
  sel_x <- cbind(sel_x,TVP,TFP)

  layout(matrix(c(1,2),ncol=2,byrow = T))
  # par(mfrow=c(1,3),mar=c(3,3,3,2))
  cols <- RColorBrewer::brewer.pal(length(method)+1,"Set1")[-6]
  col_box <- RColorBrewer::brewer.pal(9,"Pastel1")[7]

  ylim <- c(-0.05,1.1*1)
  lwd <- 1.5
  ncols <- 3
  cex.leg <- 1
  boxplot(TVP~method*rho,border=cols,col=col_box,df,main="True Positive Rate (TPR)",#border=cols[rep(1:l_m,i)]
          xlab="",xaxt="n",ylim=ylim,ylab="",lwd=lwd);#abline(h=c(0,0.0975),lty=3,lwd=2)
  abline(h=c(0,1),lty=1,col="gray")
  axis(side = 1,at = length(levels(df$method))/2+0.5+c(0:(length(rho_s)-1) )*length(levels(df$method)),
       labels = labs,cex.axis=0.9)
  abline(v=0.5+c(0:(length(rho_s)) )*length(unique(df$method)),lty=3)
  legend("top",legend = level_legend,fill = cols,
         ncol = 2,bg = "white",cex = cex.leg,bty="n")

  boxplot(TFP~method*rho,border=cols,col=col_box,df,main="False Positive Rate (FPR)",#border=cols[rep(1:l_m,i)]
          xlab="",xaxt="n",ylim=ylim,ylab="",lwd=lwd);#abline(h=c(0,0.0975),lty=3,lwd=2)
  abline(h=c(0,1),lty=1,col="gray")
  abline(v=0.5+c(0:(length(rho_s)) )*length(unique(df$method)),lty=3)
  axis(side = 1,at = length(levels(df$method))/2+0.5+c(0:(length(rho_s)-1) )*length(levels(df$method)),
       labels = labs,cex.axis=0.9)
  legend("top",legend = level_legend,fill = cols,
         ncol = 2,bg = "white",cex = cex.leg,bty="n")

}

plot_TVP_TFP_Y_rho <- function(){
  rho_s <- sort(unique(df$rho))
  n_rho_s <- length(rho_s)
  labs <- lapply(rho_s,function(rr){bquote(rho==.(rr))})
  p <- ncol(D)
  var_star <- which(colSums(abs(D))>0)
  p_star <- length(var_star)
  sel_x <- df[,c(1,3,9,10)]
  TVP <- (sel_x$VP_Y)/(p_star)
  TFP <- (sel_x$SEL_Y- sel_x$VP_Y)/(p-p_star)
  sel_x <- cbind(sel_x,TVP,TFP)

  layout(matrix(c(1,2),ncol=2,byrow = T))
  # par(mfrow=c(1,3),mar=c(3,3,3,2))
  cols <- RColorBrewer::brewer.pal(length(method)+1,"Set1")[-6]
  col_box <- RColorBrewer::brewer.pal(9,"Pastel1")[7]

  ylim <- c(-0.05,1.1*1)
  lwd <- 1.5
  ncols <- 3
  cex.leg <- 1
  boxplot(TVP~method*rho,border=cols,col=col_box,df,main="True Positive Rate (TPR)",#border=cols[rep(1:l_m,i)]
          xlab="",xaxt="n",ylim=ylim,ylab="",lwd=lwd);#abline(h=c(0,0.0975),lty=3,lwd=2)
  abline(h=c(0,1),lty=1,col="gray")
  axis(side = 1,at = length(levels(df$method))/2+0.5+c(0:(length(rho_s)-1) )*length(levels(df$method)),
       labels = labs,cex.axis=0.9)
  abline(v=0.5+c(0:(length(rho_s)) )*length(unique(df$method)),lty=3)
  legend("top",legend = level_legend,fill = cols,
         ncol = 2,bg = "white",cex = cex.leg,bty="n")

  boxplot(TFP~method*rho,border=cols,col=col_box,df,main="False Positive Rate (FPR)",#border=cols[rep(1:l_m,i)]
          xlab="",xaxt="n",ylim=ylim,ylab="",lwd=lwd);#abline(h=c(0,0.0975),lty=3,lwd=2)
  abline(h=c(0,1),lty=1,col="gray")
  abline(v=0.5+c(0:(length(rho_s)) )*length(unique(df$method)),lty=3)
  axis(side = 1,at = length(levels(df$method))/2+0.5+c(0:(length(rho_s)-1) )*length(levels(df$method)),
       labels = labs,cex.axis=0.9)
  legend("top",legend = level_legend,fill = cols,
         ncol = 2,bg = "white",cex = cex.leg,bty="n")

}

test <- function(){
  if(T){
    source('~/Documents/GitHub/ddsPLS/R/simulations_ddsPLS2(1).R')
    library(ddsPLS)
    library(doParallel)

    NCORES <- 15
    NZV <- 1e-2

    eps1=eps2=eps3=epsY=eps <- 0.99#0.8#
    n <- c(100)#unique(round(seq(20,300,length.out = 8)))#unique(round(seq(20,150,length.out = 5)))
    NCORES_S <- 6
    ALPHA <- rep(1/2,5)
    n_Bs <- 300#c(1000,500,300,100,100)
    Ns <- 1:100
    rhos <- 1

    ff <- function(cc){out <- cc/sqrt(sum(cc^2));if(is.na(out[1])) out <- cc;out}
    p1 = p2 = p3 <- 50; p <- 1000
    R1 <- 3; R2 <- 2; R3 <- 2;d <- R1+R2+R3

    paras <- expand.grid(rhos,Ns)
    NNs <- nrow(paras)
    Q2_cum = error_B <- rep(NA,NNs)
    method <- c("ddsPLS Boot",
                "ddsPLS Unik Boot",
                "sPLS Boot",
                "PLS Boot",
                "SPLS Boot")
    level_legend <- c("ddsPLS Boot","ddsPLS Unik Boot","sPLS [4] Boot",
                      "PLS Boot","SPLS [16] Boot")
    l_m <- length(method)
    names_df <- c("rho","id","method","Q2","DIST_B_hat","NORM_B_hat",
                  "VP_X","SEL_X","VP_Y","SEL_Y","R","time")
    mat_0 <- expand.grid(rhos,Ns,method)
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
    file_data <- "/Users/hlorenzo/Documents/GitHub/data_last_with_rho_bis.RData"
    # load(file_data)
  }
  # posINIT <- which(df$method==method[1] & df$R>2)
  for(i in 28:nrow(paras)){
    rho <- paras[i,1]
    A <- rbind(
      matrix(rep(c(rep(1,p1),rep(1,p2),rep(0,p-p1-p2)),R1),nrow = R1,byrow = T),
      matrix(rep(c(rep(1,p1),rep(0,p2),rep(0,p-p1-p2)),R2),nrow = R2,byrow = T),
      matrix(rep(c(rep(0,p1),rep(0,p2),rep(1,p3),rep(0,p-p1-p2-p3)),R3),nrow = R3,byrow = T)
    )
    A <- eps1*apply(A,2,ff)
    D <- rbind(
      matrix(rep(c(0,0,0,0),R1),nrow = R1,byrow = T),
      matrix(rep(c(1,0,0,0),R2),nrow = R2,byrow = T),
      matrix(rep(c(0,1,0,0),R3),nrow = R3,byrow = T)
    )
    D <- eps1*apply(D,2,ff)
    q <- ncol(D)
    LAMBDAS <- seq(0,1,length.out = 66)
    KXS <- unique(sort(c(round(seq(1,p,length.out = 22)),100,50,150) ))
    KYS <- unique(round(seq(1,q,length.out = 3)))
    L_total <- q+p
    A_all <- A
    B_th_all=B_th <- MASS::ginv(A_all)%*%D
    cat("\n\n______________________________________________________________________________")
    cat(paste("\n rho =",rho,"  ---   i =",i))
    pos <- intersect(which(df$rho==rho),which(df$id==paras[i,2]))
    pos_method <- unlist(lapply(method,function(mm){pos[which(df$method[pos]==mm)]}))
    # # Do data
    psi <- MASS::mvrnorm(n,mu = rep(0,d+L_total),Sigma = diag(d+L_total))
    phi <- psi[,1:d,drop=F];pt <- d
    # SIs <- lapply(list(A1,A2,A3,C),function(M){do.call(cbind,lapply(sqrt(1-diag(crossprod(M))),function(sisi){rep(sisi,n)}))})
    # X1 <- phi%*%A1 + SIs[[1]]*psi[,pt+1:ncol(A1),drop=F];pt <- pt + ncol(A1)
    # X2 <- phi%*%A2 + SIs[[2]]*psi[,pt+1:ncol(A2),drop=F];pt <- pt + ncol(A2)
    # X3 <- phi%*%A3 + SIs[[3]]*psi[,pt+1:ncol(A3),drop=F];pt <- pt + ncol(A3)
    SIs <- lapply(list(A,D),function(M){do.call(cbind,lapply(sqrt(1-diag(crossprod(M))),function(sisi){rep(sisi,n)}))})
    Xs <- list(x=phi%*%A + SIs[[1]]*psi[,pt+1:ncol(A),drop=F]);pt <- pt + ncol(A)
    Y <- phi%*%D + SIs[[2]]*psi[,pt+1:ncol(D),drop=F];pt <- pt + ncol(Y)
    datas$Xs[[i]] <- Xs
    datas$Y[[i]] <- Y
    datas$phi[[i]] <- phi
    if(i%%50==0){
      save(datas,df,ALL_FUCKING_MODELS,file = file_data)
      #save(datas,df,file = "../data_simu/data_signalFaible.RData")
      # save(datas,df,varExplained,LAMBDAS_SOL,file = "../../Hadrien/data_signalFaible.RData")#save(datas,df,file = "../data_simu/data_signalFaible.RData")
    }
    # Load data
    # datas$phi[[i]] -> phi #; phi_ok <- phi[,1:2]
    # datas$Xs[[i]] -> Xs
    # datas$Y[[i]] -> Y
    x <- do.call(cbind,Xs)
    ##############
    for(i_m in 1:length(pos)){
      pos_i <- pos[i_m]
      method_i <- method[which(pos_method==pos_i)]
      cat(paste("\n     <- ",method_i,"... ",sep=""))
      toPlot <- method_i %in% method[c(1,3,4)]
      if(toPlot){
        time_1 <- Sys.time()
        if(method_i %in% c("ddsPLS Boot","PLS Boot")){
          lambdas <- LAMBDAS
          n_b_i <- n_Bs
          ncores_i <- NCORES_S
          if(method_i=="PLS Boot"){
            lambdas <- 0
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
          n_b_i <- n_Bs
          ncores_i <- NCORES_S
          res <- Q2_UNIK_ddsPLS(Xs,Y,lambdas=lambdas,
                                n_B = n_b_i,NCORES=ncores_i,verbose = T)
        }else if(method_i %in% c("sPLS Boot") ){
          n_b_i <- n_Bs
          ncores_i <- NCORES_S
          kxs <- KXS;kys <- KYS; sparse <- T
          if(!(method_i %in% c("sPLS Boot")) ){
            kxs <- ncol(x);kys <- ncol(Y); sparse <- F;ncores_i <- 7
          }
          res <- Q2_boot_sPLS(Xs,Y,keepXs = kxs,keepYs=kys,
                              n_B=n_b_i,deflatX=T,NCORES=ncores_i,center=T,
                              NZV=1e-9,verbose=T,sparse = sparse)
        }
        if(!is.null(res)){
          df[pos_i,]$Q2 <- res$optimal_parameters$Q2
          df[pos_i,]$DIST_B_hat <- 1-sum((A%*%res$B_cbind-D)^2)/sum(D^2)
          df[pos_i,]$NORM_B_hat <- sqrt(sum((res$B_cbind)^2))
          df[pos_i,]$R <- res$optimal_parameters$R
          df[pos_i,id_sel] <- compare_selection(res$B_cbind,B_th_all)
        }
        ALL_FUCKING_MODELS[[method_i]][[i]] <- res
        df$time[pos_i] <- as.numeric(difftime(Sys.time(), time_1, units = "secs"))
      }
      cat("\n\n______________________________________________________________________________")
      ##########################
      ########## PLOT ##########
      ##########################
      if(toPlot & length(pos)==length(pos)){
        pdf(file = "/Users/hlorenzo/Documents/GitHub/Simulations_2_bis.pdf",width = 15,height = 8)
        # postscript("/Users/hlorenzo/Dropbox/Results_Last/Simulations.eps", width=30, height=10, onefile=TRUE, horizontal=T)
        give_me_plot_comm_bis()
        dev.off()

        pdf(file = "/Users/hlorenzo/Documents/GitHub/Simulations_R_2_bis.pdf",width = 12,height = 10)
        plot_R_rho()
        dev.off()

        pdf(file = "/Users/hlorenzo/Documents/GitHub/Simulations_sel_x_2_bis.pdf",width = 14,height = 9)
        # postscript("/Users/hlorenzo/Dropbox/Results/Simulations_sel_x.eps", width=14, height=9, onefile=TRUE, horizontal=FALSE)
        # plot_X()
        plot_TVP_TFP_X_rho()
        dev.off()

        pdf(file = paste("/Users/hlorenzo/Documents/GitHub/Simulations_sel_y_2_bis.pdf",sep=""),width = 14,height = 9)
        # postscript("/Users/hlorenzo/Dropbox/Results/Simulations_sel_y.eps", width=14, height=9, onefile=TRUE, horizontal=F)
        # plot_sel_simu_y()
        plot_TVP_TFP_Y_rho()
        dev.off()
      }
    }
  }
}
