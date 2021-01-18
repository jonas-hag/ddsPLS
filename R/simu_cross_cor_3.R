give_me_plot_comm_rho <- function(){
  rho_s <- sort(unique(df$rho))
  n_rho_s <- length(rho_s)
  labs <- lapply(rho_s,function(rr){bquote(rho==.(rr))})
  par(mar=c(3,3,3,2))
  # layout(matrix(c(1,2,
  #                 1,2,
  #                 3,4,
  #                 3,5), 4, byrow = TRUE),mar=c(3,3,3,2))
  # layout(matrix(c(1,2), nrow = 1, byrow = TRUE))
  # layout(matrix(c(1,2,1,2,1,2,1,2,3,4,3,4,3,4,3,4,rep(5,2)),ncol=2,byrow = T))
  layout(matrix(c(rep(c(1,2),10),3,3),ncol=2,byrow = T))
  # par(mfrow=c(1,3),mar=c(3,3,3,2))
  cols <- RColorBrewer::brewer.pal(length(method)+1,"Set1")[-6]
  col_box <- RColorBrewer::brewer.pal(9,"Pastel1")[7]

  ylim <- c(-0.05,1.1*(2*eps^2)/3)
  lwd <- 1.5
  ncols <- 3
  cex.leg <- 1
  b<-boxplot(Q2~method*rho,border=cols,col=col_box,df,main=expression("Q"^"2"),#border=cols[rep(1:l_m,i)]
             xlab=expression(rho),xaxt="n",ylab="",lwd=lwd,pch=1:length(level_legend));#abline(h=c(0,0.0975),lty=3,lwd=2)
  id_meth <- rep(1:nlevels(df$method),length(unique(df$rho)))
  coli_s <- cols[id_meth]
  for(ii in 1:length(coli_s)){
    rect(xleft = ii-0.5+0.1,xright = ii+0.5-0.1,
         ytop = b$stats[4,ii],ybottom = b$stats[2,ii],
         density = density[id_meth[ii]],angle=angle[id_meth[ii]],col=coli_s[ii])
  }
  abline(h=(0:20)/20,lty=2,lwd=1/2,col="gray")
  abline(v=0.5+c(0:(length(n_rho_s)+1) )*length(levels(df$method)),lty=3)
  abline(h=0,lty=2,col="black")
  axis(side = 1,
       at = length(levels(df$method))/2+0.5+c(0:(length(n_rho_s)) )*length(levels(df$method)),
       labels = rho_s,cex.axis=0.9)
  # legend("bottomright",legend = level_legend,fill = cols,
  #        ncol = 1,bg = "white",cex = cex.leg,bty="n")
  text(x = 10,y = 1.1*(3*eps^2)/4,labels = expression(gamma^"*"==frac(3,4)~.~epsilon^2==0.7351),pos = 3)
  axis(2,at = (3*eps^2)/4,labels = expression(gamma^"*"),las=1)
  abline(h=(3*eps^2)/4,lty=1,col="black")

  # boxplot(time~method*rho,border=cols,col=col_box,df,main="Raw time (not adjusted on number of parameters)",
  #         xlab=expression(rho),xaxt="n",ylab="",lwd=lwd);#abline(h=c(0,0.0975),lty=3,lwd=2)
  # abline(v=0.5+c(0:(length(n_rho_s)+1) )*length(levels(df$method)),lty=3)
  # abline(h=0,lty=2)
  # axis(side = 1,at = length(levels(df$method))/2+0.5+c(0:(length(rho_s)-1) )*length(levels(df$method)),
  #      labels = rho_s,cex.axis=0.9)
  # legend("topright",legend = level_legend,fill = cols,
  #        ncol = 1,bg = "white",cex = cex.leg,bty="n")

  b<-boxplot(log(1-DIST_B_hat)~method*rho,df,#[-which(df$method=="OLS_esp"),],
             main=expression("||A"~hat(B)~"-D||"^2/"||D||"^2~"(logarithmic scale)"),#log="y",
             col=col_box,border=cols,xlab=expression(rho),xaxt="n",yaxt="n",ylab="",
             lwd=lwd,pch=1:length(level_legend))#,ylim=c(0,1))
  id_meth <- rep(1:nlevels(df$method),length(unique(df$rho)))
  coli_s <- cols[id_meth]
  for(ii in 1:length(coli_s)){
    rect(xleft = ii-0.5+0.1,xright = ii+0.5-0.1,
         ytop = b$stats[4,ii],ybottom = b$stats[2,ii],
         density = density[id_meth[ii]],angle=angle[id_meth[ii]],col=coli_s[ii])
  }
  iii0 <- c((1e-1)^(rev(0:3)))
  iii <- sort(c(iii0,5*iii0))
  abline(h=log(iii),lty=2,lwd=1/2,col="gray")
  axis(2,at=log(iii),labels=iii,las=2)
  abline(v=0.5+c(0:(length(n_rho_s)+1) )*length(unique(df$method)),lty=3)
  # abline(h=0,lty=2)
  axis(side = 1,at = length(levels(df$method))/2+0.5+c(0:(length(rho_s)-1) )*length(levels(df$method)),
       labels =rho_s,cex.axis=0.9)
  # legend("topright",legend = level_legend,fill = cols,ncol = 1,bg = "white",cex = cex.leg,bty="n")

  # boxplot(NORM_B_hat~method*rho,df,#[-which(df$method=="OLS_esp"),],
  #         main=expression("||"~hat(B)~"||"^2),
  #         col=col_box,border=cols,xlab=expression(rho),xaxt="n",ylab="",lwd=lwd)#,ylim=c(0,1))
  # abline(h=0,lty=2)
  # abline(v=0.5+c(0:(length(n_rho_s)+1) )*length(unique(df$method)),lty=3)
  # axis(side = 1,at = length(levels(df$method))/2+0.5+c(0:(length(rho_s)-1) )*length(levels(df$method)),
  #      labels = rho_s,cex.axis=0.9)
  # legend("topright",legend = level_legend,fill = cols,ncol = 1,bg = "white",cex = cex.leg,bty="n")
  par(mar=rep(0,4))
  plot(0,type='n',axes=FALSE,ann=FALSE)
  legend("top",legend = level_legend,angle=angle,border=cols,fill=cols,density=density,col = cols,pch=1:length(level_legend),
         ncol = 6,bty = "n",cex = cex.leg,pt.cex=2)
}

plot_R_rho <- function(){
  rho_s <- sort(unique(df$rho))
  n_rho_s <- length(rho_s)
  labs <- lapply(rho_s,function(rr){bquote(rho==.(rr))})
  n_meth <- length(method)
  layout(matrix(1:(n_rho_s*n_meth),nrow = n_rho_s,byrow = T))
  cc <- 0.9
  par(mar=c(2,2,2,0),cex.main=cc,cex.lab=cc,cex.axis=cc)
  cols <- RColorBrewer::brewer.pal(length(method)+1,"Set1")[-6]
  breaks <- (min(c(2,na.omit(df$R)))-0.5):(max(c(2,na.omit(df$R)))+0.5)
  for(i_n in 1:n_rho_s){
    for(i_mm in 1:n_meth){
      Ri <- df$R[which(df$method==method[i_mm] & df$rho==rho_s[i_n])]
      # Ri[which(is.na(Ri))] <- 0

      hist(rep(3,100),breaks = breaks,xaxt="n",yaxt="n",col="gray90",ylim=c(0,110),
           main="",xlab=expression(rho),ylab="",density=20 , angle=33,border="gray80")
      axis(side = 1,unique(na.omit(df$R)),line = -1,
           labels = unique(na.omit(df$R)),tick = F,gap.axis=1/6)
      axis(side = 2,(0:10)*10,line = -1,labels = (0:10)*10,tick = F,las=2)
      abline(h=(0:5)*20,lty=2,col="gray")
      hist(Ri,breaks = breaks,col=cols[i_mm],add=T,border = "gray20",
           density=density[i_mm] , angle=angle[i_mm])
      title(bquote(.(level_legend[i_mm])~"\n "~rho==.(rho_s[i_n])),line = -1)
      title(xlab=expression(hat(R)),ylab="Frequency",line = 1)
      xx <- unique(na.omit(df$R))
      yy <- unlist(lapply(xx,function(xxx){length(which(Ri==xxx))}))
      if(length(which(yy!=0))>0){
        text(xx[which(yy!=0)],yy[which(yy!=0)],labels = yy[which(yy!=0)],pos = 3,col=1)#cols[i_mm])
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

  layout(matrix(c(rep(c(1,2,3),15),rep(4,3)),ncol=3,byrow = T))
  par(mar=c(3.1, 3.1, 4.1, 0.2))
  cols <- RColorBrewer::brewer.pal(length(method)+1,"Set1")[-6]
  col_box <- "white"#RColorBrewer::brewer.pal(9,"Pastel1")[7]

  ylim <- c(-0.05,1.4*1)
  lwd <- 1.5
  ncols <- 3
  cex.leg <- 1

  b<-boxplot(SEL_X~method*rho,border=cols,col=col_box,df,main="Number of selected variables",#border=cols[rep(1:l_m,i)]
             xlab=expression(rho),xaxt="n",ylab="",lwd=lwd,pch=1:length(level_legend))
  id_meth <- rep(1:nlevels(df$method),length(unique(df$rho)))
  coli_s <- cols[id_meth]
  for(ii in 1:length(coli_s)){
    rect(xleft = ii-0.5+0.1,xright = ii+0.5-0.1,
         ytop = b$stats[4,ii],ybottom = b$stats[2,ii],
         density = density[id_meth[ii]],angle=angle[id_meth[ii]],col=coli_s[ii])
  }
  abline(h=c(150),lty=1,col=1)
  axis(side = 1,at = length(levels(df$method))/2+0.5+c(0:(length(rho_s)-1) )*length(levels(df$method)),
       labels = rho_s,cex.axis=0.9)
  abline(v=0.5+c(0:(length(rho_s)) )*length(unique(df$method)),lty=3)
  # legend("top",legend = level_legend,fill = cols,
  #        ncol = 2,bg = "white",cex = cex.leg,bty="n")

  b<-boxplot(TVP~method*rho,border=cols,col=col_box,df,main="True Positive Rate (TPR)",#border=cols[rep(1:l_m,i)]
             xlab=expression(rho),xaxt="n",ylab="",lwd=lwd,pch=1:length(level_legend));#abline(h=c(0,0.0975),lty=3,lwd=2)
  id_meth <- rep(1:nlevels(df$method),length(unique(df$rho)))
  coli_s <- cols[id_meth]
  for(ii in 1:length(coli_s)){
    rect(xleft = ii-0.5+0.1,xright = ii+0.5-0.1,
         ytop = b$stats[4,ii],ybottom = b$stats[2,ii],
         density = density[id_meth[ii]],angle=angle[id_meth[ii]],col=coli_s[ii])
  }
  abline(h=1,lty=1,col=1)
  axis(side = 1,at = length(levels(df$method))/2+0.5+c(0:(length(rho_s)-1) )*length(levels(df$method)),
       labels = rho_s,cex.axis=0.9)
  abline(v=0.5+c(0:(length(rho_s)) )*length(unique(df$method)),lty=3)
  # legend("top",legend = level_legend,fill = cols,
  #        ncol = 2,bg = "white",cex = cex.leg,bty="n")

  b<-boxplot(TFP~method*rho,border=cols,col=col_box,df,main="False Positive Rate (FPR)",#border=cols[rep(1:l_m,i)]
             xlab=expression(rho),xaxt="n",ylab="",lwd=lwd,pch=1:length(level_legend));#abline(h=c(0,0.0975),lty=3,lwd=2)
  id_meth <- rep(1:nlevels(df$method),length(unique(df$rho)))
  coli_s <- cols[id_meth]
  for(ii in 1:length(coli_s)){
    rect(xleft = ii-0.5+0.1,xright = ii+0.5-0.1,
         ytop = b$stats[4,ii],ybottom = b$stats[2,ii],
         density = density[id_meth[ii]],angle=angle[id_meth[ii]],col=coli_s[ii])
  }
  abline(h=0,lty=1,col=1)
  abline(v=0.5+c(0:(length(rho_s)) )*length(unique(df$method)),lty=3)
  axis(side = 1,at = length(levels(df$method))/2+0.5+c(0:(length(rho_s)-1) )*length(levels(df$method)),
       labels = rho_s,cex.axis=0.9)
  # legend("top",legend = level_legend,fill = cols,
  #        ncol = 2,bg = "white",cex = cex.leg,bty="n")
  par(mar=rep(0,4))
  plot(0,type='n',axes=FALSE,ann=FALSE)
  legend("top",legend = level_legend,angle=angle,border=cols,fill=cols,density=density,col = cols,pch=1:length(level_legend),
         ncol = 6,bty = "n",cex = cex.leg,pt.cex=2)

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

  layout(matrix(c(rep(c(1,2,3),15),rep(4,3)),ncol=3,byrow = T))
  par(mar=c(3.1, 3.1, 4.1, 0.2))
  cols <- RColorBrewer::brewer.pal(length(method)+1,"Set1")[-6]
  col_box <- RColorBrewer::brewer.pal(9,"Pastel1")[7]

  ylim <- c(-0.05,1.4*1)
  lwd <- 1.5
  ncols <- 3
  cex.leg <- 1

  b<-boxplot(SEL_Y~method*rho,border=cols,col=col_box,df,main="Number of selected variables",#border=cols[rep(1:l_m,i)]
             xlab=expression(rho),xaxt="n",ylab="",lwd=lwd,pch=1:length(level_legend))
  id_meth <- rep(1:nlevels(df$method),length(unique(df$rho)))
  coli_s <- cols[id_meth]
  for(ii in 1:length(coli_s)){
    rect(xleft = ii-0.5+0.1,xright = ii+0.5-0.1,
         ytop = b$stats[4,ii],ybottom = b$stats[2,ii],
         density = density[id_meth[ii]],angle=angle[id_meth[ii]],col=coli_s[ii])
  }
  abline(h=c(3),lty=1,col=1)
  axis(side = 1,at = length(levels(df$method))/2+0.5+c(0:(length(rho_s)-1) )*length(levels(df$method)),
       labels = rho_s,cex.axis=0.9)
  abline(v=0.5+c(0:(length(rho_s)) )*length(unique(df$method)),lty=3)
  # legend("top",legend = level_legend,fill = cols,
  #        ncol = 2,bg = "white",cex = cex.leg,bty="n")

  b<-boxplot(TVP~method*rho,border=cols,col=col_box,df,main="True Positive Rate (TPR)",#border=cols[rep(1:l_m,i)]
             xlab=expression(rho),xaxt="n",ylab="",lwd=lwd,pch=1:length(level_legend));#abline(h=c(0,0.0975),lty=3,lwd=2)
  id_meth <- rep(1:nlevels(df$method),length(unique(df$rho)))
  coli_s <- cols[id_meth]
  for(ii in 1:length(coli_s)){
    rect(xleft = ii-0.5+0.1,xright = ii+0.5-0.1,
         ytop = b$stats[4,ii],ybottom = b$stats[2,ii],
         density = density[id_meth[ii]],angle=angle[id_meth[ii]],col=coli_s[ii])
  }
  abline(h=1,lty=1,col=1)
  axis(side = 1,at = length(levels(df$method))/2+0.5+c(0:(length(rho_s)-1) )*length(levels(df$method)),
       labels = rho_s,cex.axis=0.9)
  abline(v=0.5+c(0:(length(rho_s)) )*length(unique(df$method)),lty=3)
  # legend("top",legend = level_legend,fill = cols,
  #        ncol = 2,bg = "white",cex = cex.leg,bty="n")

  b<-boxplot(TFP~method*rho,border=cols,col=col_box,df,main="False Positive Rate (FPR)",#border=cols[rep(1:l_m,i)]
             xlab=expression(rho),xaxt="n",ylab="",lwd=lwd,pch=1:length(level_legend));#abline(h=c(0,0.0975),lty=3,lwd=2)
  id_meth <- rep(1:nlevels(df$method),length(unique(df$rho)))
  coli_s <- cols[id_meth]
  for(ii in 1:length(coli_s)){
    rect(xleft = ii-0.5+0.1,xright = ii+0.5-0.1,
         ytop = b$stats[4,ii],ybottom = b$stats[2,ii],
         density = density[id_meth[ii]],angle=angle[id_meth[ii]],col=coli_s[ii])
  }
  abline(h=0,lty=1,col=1)
  abline(v=0.5+c(0:(length(rho_s)) )*length(unique(df$method)),lty=3)
  axis(side = 1,at = length(levels(df$method))/2+0.5+c(0:(length(rho_s)-1) )*length(levels(df$method)),
       labels = rho_s,cex.axis=0.9)
  # legend("top",legend = level_legend,fill = cols,s
  #        ncol = 2,bg = "white",cex = cex.leg,bty="n")
  par(mar=rep(0,4))
  plot(0,type='n',axes=FALSE,ann=FALSE)
  legend("top",legend = level_legend,angle=angle,border=cols,fill=cols,density=density,col = cols,pch=1:length(level_legend),
         ncol = 6,bty = "n",cex = cex.leg,pt.cex=2)
}

test <- function(){
  if(T){
    source('~/Documents/GitHub/ddsPLS/R/simulations_ddsPLS2(1).R')
    library(ddsPLS)
    library(doParallel)

    NCORES <- 15
    NZV <- 1e-2

    eps1=eps2=eps3=epsY=eps <- 0.9#0.8#
    n <- c(100)#unique(round(seq(20,300,length.out = 8)))#unique(round(seq(20,150,length.out = 5)))
    NCORES_S <- 22
    ALPHA <- rep(1/2,5)
    n_Bs <- 300#c(1000,500,300,100,100)
    Ns <- 1:100
    rhos <- c(0.1,0.9)

    ff <- function(cc){out <- cc/sqrt(sum(cc^2));if(is.na(out[1])) out <- cc;out}
    p1 = p2 = p3 <- 50; p <- 1000
    q <- 4
    R1 <- 3; R2 <- 2; R3 <- 2;d <- R1+R2+R3

    LAMBDAS <- seq(0,1,length.out = 66)
    KXS <- unique(sort(c(round(seq(1,p,length.out = 22)),100) ))
    KYS <- unique(round(seq(1,q,length.out = 3)))


    paras <- expand.grid(rhos,Ns)
    NNs <- nrow(paras)
    Q2_cum = error_B <- rep(NA,NNs)
    method <- c("ddsPLS Boot",
                "ddsPLS Unik Boot",
                "sPLS Boot",
                "PLS Boot",
                "SPLS Boot",
                "sPLS classik")
    level_legend <- c("ddsPLS Boot","ddsPLS Unik Boot","sPLS [4] Boot",
                      "PLS Boot","SPLS [16] Boot","sPLS classik")
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
    file_data <- "/Users/hlorenzo/Documents/GitHub/data_last_with_rho3.RData"
    load(file_data)
    pos_old_unik <- which(df$method=="ddsPLS Unik Boot")
    if(length(pos_old_unik)>0){
      levels(df$method)[which(levels(df$method)=="ddsPLS Unik Boot")] <- "sPLS Boot Var-Q2"
      df$method[pos_old_unik] <- "sPLS Boot Var-Q2"
    }
    method[2] <- "sPLS Boot Var-Q2"
    level_legend[1] <- "ddsPLS"
    level_legend[2] <- "sPLS [7] 1"
    level_legend[3] <- "sPLS [7] 2"
    level_legend[4] <- "NIPALS-PLS"
    level_legend[5] <- "SPLS [1]"
    level_legend[6] <- "sPLS [7] 3"

    angle <- c(0,45,90,0,135,45,90)
    density <- c(50,50,50,25,25,50,25)
  }
  # posINIT <- which(df$method==method[1] & df$R>2)
  i = i_m <- 1
  # posNA <- which(is.na(df$Q2))
  for(i in 1:200){
    rho <- paras[i,1]
    A <- rbind(
      matrix(rep(c(rep(1,p1),rep(0,p-p1)),R1),nrow = R1,byrow = T),
      matrix(rep(c(rep(0,p1),rep(1,p2),rep(0,p-p1-p2)),R2),nrow = R2,byrow = T),
      matrix(rep(c(rep(0,p1),rep(1,p2),rep(1,p3),rep(0,p-p1-p2-p3)),R3),nrow = R3,byrow = T)
    )
    A <- eps1*apply(A,2,ff)
    al1 <- (sqrt(1+rho)+sqrt(1-rho))/2;al2 <- (sqrt(1+rho)-sqrt(1-rho))/2
    D <- rbind(
      matrix(rep(c(1,0,0,0),R1),nrow = R1,byrow = T),
      matrix(rep(c(0,al1,al2,0),R2),nrow = R2,byrow = T),
      matrix(rep(c(0,al2,al1,0),R3),nrow = R3,byrow = T)
    )
    D <- eps1*apply(D,2,ff)
    q <- ncol(D)
    L_total <- q+p
    A_all <- A
    B_th_all=B_th <- MASS::ginv(A_all)%*%D
    cat("\n\n______________________________________________________________________________")
    cat(paste("\n rho =",rho,"  ---   i =",i))
    pos <- intersect(which(df$rho==rho),which(df$id==paras[i,2]))
    pos_method <- unlist(lapply(method,function(mm){pos[which(df$method[pos]==mm)]}))
    # Do data
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
      # toPlot <- pos_i %in% posNA#method_i %in% method[1:4]
      toPlot<-T
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
        }else if(method_i %in% c("sPLS Boot","sPLS Boot Var-Q2") ){
          n_b_i <- n_Bs
          ncores_i <- NCORES_S
          kxs <- KXS;kys <- KYS; sparse <- T
          if(method_i == c("sPLS Boot")){
            whichCriterion <- "Q2"
          }else{
            whichCriterion <- "Q2-Var"
          }
          res <- Q2_boot_sPLS(Xs,Y,keepXs = kxs,keepYs=kys,whichCriterion = whichCriterion,
                              n_B=n_b_i,deflatX=T,NCORES=ncores_i,center=T,
                              NZV=1e-9,verbose=T,sparse = sparse)
        }else if(method_i %in% "SPLS Boot"){
          etas <- seq(0.1,0.9,0.1)
          Ks <- 1:10
          fold <- min(n,100)
          cl <- makeCluster(length(etas))
          registerDoParallel(cl)
          res_all <- foreach(et=1:length(etas),.packages = "spls",.combine='rbind') %dopar% {
            cvmo <- cv.spls(Xs[[1]],Y,fold=fold,K = Ks,eta = etas[et],scale.y = T,plot.it = F)
            c(max(1-cvmo$mspemat),cvmo$K.opt,etas[et])
          }
          stopCluster(cl)
          id_max <- which.max(res_all[,1])
          mo <- spls::spls(Xs[[1]],Y,K = res_all[id_max,2],eta = res_all[id_max,3])
          res <- list(optimal_parameters=list(Q2=max(res_all[,1]),R=as.numeric(res_all[id_max,2])),B_cbind=mo$betahat)
        }else if(method_i %in% "sPLS classik"){
          ncores_i <- NCORES_S
          kxs <- KXS;kys <- KYS; sparse <- T
          res <- do_mixomics(Xs,Y,kxs,kys,ncores_i)
        }
        if(!is.null(res)){
          df[pos_i,]$Q2 <- res$optimal_parameters$Q2
          df[pos_i,]$DIST_B_hat <- 1-sum((A%*%res$B_cbind-D)^2)/sum(D^2)
          df[pos_i,]$NORM_B_hat <- sqrt(sum((res$B_cbind)^2))
          df[pos_i,]$R <- res$optimal_parameters$R
          if(method_i %in% c("sPLS classik",method[c(2,3)])){
            if(is.list(res$U)){
              df[pos_i,id_sel] <- compare_selection(tcrossprod(res$U[[1]],res$V),B_th_all)
            }else{
              df[pos_i,id_sel] <- compare_selection(tcrossprod(res$U,res$V),B_th_all)
            }
          }else{
            df[pos_i,id_sel] <- compare_selection(res$B_cbind,B_th_all)
          }
        }
        ALL_FUCKING_MODELS[[method_i]][[i]] <- res
        df$time[pos_i] <- as.numeric(difftime(Sys.time(), time_1, units = "secs"))
      }
      cat("\n\n______________________________________________________________________________")
      ##########################
      ########## PLOT ##########
      ##########################
      if(toPlot & i>1){
        # pdf(file = "/Users/hlorenzo/Documents/GitHub/Simulations_rho.pdf",width = 15,height = 5)
        postscript("/Users/hlorenzo/Documents/GitHub/Simulations_rho3.eps", onefile=TRUE, horizontal=F,
                   width=6, height=2,pointsize=0.1)
        give_me_plot_comm_rho()
        dev.off()

        # pdf(file = "/Users/hlorenzo/Documents/GitHub/Simulations_R_rho.pdf",width = 15,height = 5)
        postscript("/Users/hlorenzo/Documents/GitHub/Simulations_R_rho3.eps", onefile=TRUE, horizontal=F,
                   width=11, height=3,pointsize=8)
        plot_R_rho()
        dev.off()

        # pdf(file = "/Users/hlorenzo/Documents/GitHub/Simulations_sel_x_rho.pdf",width = 15,height = 5)
        postscript("/Users/hlorenzo/Documents/GitHub/Simulations_sel_x_rho3.eps", onefile=TRUE, horizontal=F,
                   width=6, height=2,pointsize=0.1)
        # postscript("/Users/hlorenzo/Dropbox/Results/Simulations_sel_x.eps", width=14, height=9, onefile=TRUE, horizontal=FALSE)
        # plot_X()
        plot_TVP_TFP_X_rho()
        dev.off()

        # pdf(file = paste("/Users/hlorenzo/Documents/GitHub/Simulations_sel_y_rho.pdf",sep=""),width = 15,height = 5)
        postscript("/Users/hlorenzo/Documents/GitHub/Simulations_sel_y_rho3.eps", onefile=TRUE, horizontal=F,
                   width=6, height=2,pointsize=0.1)
        # postscript("/Users/hlorenzo/Dropbox/Results/Simulations_sel_y.eps", width=14, height=9, onefile=TRUE, horizontal=F)
        # plot_sel_simu_y()
        plot_TVP_TFP_Y_rho()
        dev.off()
      }
    }
  }
}
