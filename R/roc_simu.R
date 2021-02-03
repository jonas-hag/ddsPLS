# import numpy as np
# import pyPLS
# from matplotlib import pyplot as plt
# np.random.seed(5)
# X1, S1, C1 = pyPLS.simulateData(50, 4, 500, 20, signalToNoise=100, ConcVarFact=0.8, n_min_peaks=5)
# X2, S2, C2 = pyPLS.simulateData(50, 4, 5000, 30, signalToNoise=100, ConcVarFact=0.8, n_min_peaks=5)
# Ct = np.zeros((50,3))
# Ct[:,0:3] = C1[:,0:3]
# X2 = Ct @ S2[0:3,:] + 4*C2[:,3:] @ S2[3:,:] + np.random.normal(size=(50, 5000))
# np.savetxt("X1-v2.csv", X, delimiter=",")
# np.savetxt("X2-v2.csv", X2, delimiter=",")
# np.savetxt("S1.csv", S, delimiter=",")
# np.savetxt("S2.csv", S2, delimiter=",")
# np.savetxt("Y-v2.csv", Ct, delimiter=",")

simulateData <- function(n, ncp, p, sigma, sigmaNoise=0.1, ConcVarFact=0.8, n_min_peaks=3){
  meanPeakSigma <- sigma
  sigPeakSigma <- sigma / 4
  axis <- 1:p
  S <- matrix(0,ncp,p)
  C <- matrix(0,n,ncp)
  for(i in 1:ncp){
    npeaks <- 3+ceiling(10*runif(1))
    peakheights <- runif(npeaks)
    sigmas <- runif(npeaks) * sigPeakSigma + meanPeakSigma
    position <- runif(npeaks) * p
    for(j in 1:npeaks){
      S[i,] <- S[i,] + peakheights[j] * exp(-0.5 * ((axis - position[j]) / sigmas[j])^2)
    }
  }
  meanC <- sort(10^runif(ncp),decreasing = T)
  varC <- ConcVarFact * meanC * runif(ncp)
  for(i in 1:ncp){
    C[,i] <- rnorm(n = n,mean = meanC[i], sd = varC[i]/2)
  }
  X <- C%*%S;X <- X/max(abs(X))
  E <- matrix(rnorm(n*p,sd = sigmaNoise),nrow = n,ncol = p)
  X <- X*sqrt(1-sigmaNoise^2) + E
  list(X=X, C=C, S=S, E=E)
}

simulateMulti <- function(seed=1,n=50,q=5,p1=500,p2=5000,
                          sigma1=0.01,sigma2=0.01,sigmaY=0.1,
                          ncpX=10,ncpXCom=5,ncpXYCom=3,plot=F){
  set.seed(seed)
  # ncpX for each X separately
  Data_1 <- simulateData(n=n, ncp=ncpX, p=p1, sigma=20, sigmaNoise = sigma1, ConcVarFact=0.8, n_min_peaks=5)
  Data_2 <- simulateData(n=n, ncp=ncpX, p=p2, sigma=30, sigmaNoise = sigma2, ConcVarFact=0.8, n_min_peaks=5)
  S1 <- Data_1$S;C1 <- Data_1$C;X1 <- Data_1$X
  S2 <- Data_2$S;C2 <- Data_2$C
  # ncpXCom in common
  C2[,1:ncpXCom] <- C1[,1:ncpXCom,drop=F]
  X2 <- C2%*%S2;X2 <- X2/max(abs(X2))
  E2 <- matrix(rnorm(n*p2,sd = sigma2),nrow = n,ncol = p2)
  X2 <- X2*sqrt(1-sigma2^2) + E2
  # Build Y on ncpXYCom components
  Y <- scale(Data_1$C[,1:ncpXYCom,drop=F])
  # Add extra variables and noise
  E_y <- matrix(rnorm(n*q,sd = sigmaY),nrow = n,ncol = q)
  if(q>ncpXYCom){
    Y <- cbind(Y,matrix(rnorm(n*(q-ncpXYCom)),nrow = n))*sqrt(1-sigmaY^2) + E_y
  }else{
    Y <- Y*sqrt(1-sigmaY^2) + E_y
  }
  if(plot){
    layout(matrix(c(1,2,2,3,3,3),nrow = 2,byrow = T))
    matplot(t(X1),lty=1,type="l")
    matplot(t(X2),lty=1,type="l")
    corXY <- cor(cbind(X1,X2),Y)
    matplot(corXY,type="l")
  }
  list(Xs=list(X1=X1,X2=X2),Y=Y,S=list(S1=S1,S2=S2,SXY=list(S1=S1[1:ncpXYCom,],S2=S2[1:ncpXYCom,])))
}

get_tpr_fpr_auc <- function(p,U,S,verbose=T,mode="norm2"){
  Rs <- ncol(U)
  out <- list()
  for(R in 1:Rs){
    sel <- rep(0,p);sel[which(rowSums(abs(U[,1:R,drop=F]))>1e-9)] <- 1
    if(sum(sel)!=p & sum(sel)!=0){
      if(mode=="max"){
        pred <- ROCR::prediction(apply(S,2,max),sel)
      }else if(mode=="max"){
        pred <- ROCR::prediction(apply(S,2,sum),sel)
      }else if(mode=="norm2"){
        pred <- ROCR::prediction(apply(S,2,function(ss){sum(ss^2)}),sel)
      }
      perf <- ROCR::performance(pred,"tpr","fpr")
      fpr <- perf@x.values[[1]]
      tpr <- perf@y.values[[1]]
      perf <- ROCR::performance(pred,"auc")
      auc <- perf@y.values[[1]]
    }else{
      fpr <- NA
      tpr <- NA
      auc <- NA
    }
    out[[R]] <- list(fpr=fpr,tpr=tpr,auc=auc)
  }
  if(verbose){
    cat(unlist(lapply(out,function(rr){rr$auc})));cat("\n")
  }
  out
}




# X1 <- read.csv("../Sartorius1/data/newSartorius/X1-v2.csv",header = F)
# X2 <- read.csv("../Sartorius1/data/newSartorius/X2-v2.csv",header = F)
# Y <- read.csv("../Sartorius1/data/newSartorius/Y-v2.csv",header = F)
# S1 <- read.csv("../Sartorius1/data/newSartorius/S1.csv",header = F)
# S2 <- read.csv("../Sartorius1/data/newSartorius/S2.csv",header = F)

# S <- cbind(S1,S2)
# Xs <- list(X1,X2)

test <- function(){
  library(doParallel)
  library(ddsPLS)
  source('~/Documents/GitHub/ddsPLS/R/simulations_ddsPLS2(1).R')

  p <- 5500
  B <- 500;lambdas <- seq(0,1,length.out = 100)
  kxs <- unique(round(seq(1,p,length.out = 20)));kys <- 1:5; sparse <- T
  etas <- seq(0.1,0.9,0.1);Ks <- 1:5;fold <- 50
  level_legend <- c("ddsPLS","sPLS [7] 1","sPLS [7] 2","NIPALS-PLS","SPLS [1]","sPLS [7] 3")
  id_meth <- c(1,2,3,5,6)
  cols <- RColorBrewer::brewer.pal(6+1,"Set1")[-6]
  col_box <- RColorBrewer::brewer.pal(9,"Pastel1")[7]
  angle <- c(0,45,90,0,135,45,90)
  density <- c(50,50,50,25,25,50,25)
  sigma_y <- 0.2


  # RESULTS <- list()

  name_file <- "/Users/hlorenzo/Documents/GitHub/results_design_3.RData"
  name_pdf <- "/Users/hlorenzo/Documents/GitHub/results_design_3_plot"

  Q2_all=R_all <- matrix(NA,20,5)
  colnames(Q2_all)=colnames(R_all) <- level_legend[id_meth]
  postscript("/Users/hlorenzo/Documents/GitHub/simu_cloarec_eps.eps", onefile=TRUE, horizontal=F,
             width=6, height=1.5,pointsize=1)
  par(mar=c(4.5,4,3,2))
  layout(matrix(c(1,2,2,2),nrow = 1,byrow = T))
  matplot(t(datas$Xs$X1[1:5,]),lty=1,type="l",ylab="Intensities",
          main=expression(bold(X)[1]~", Spectroscopic data"),
          xlab="Variable index",lwd=0.1,col=5:1+1)
  matplot(t(datas$Xs$X2[1:5,]),lty=1,type="l",ylab="Intensities",
          main=expression(bold(X)[2]~", NMR spectroscopic data"),
          xlab="Variable index",lwd=0.1/3,col=5:1+1)
  dev.off()
  for(seed in 3:100){
    cat("\n_________________________________________________\n")
    cat(seed)
    # RESULTS[[seed]] <- list()
    datas <- simulateMulti(seed=seed,sigma1 = 0.05,sigma2 = 0.05,sigmaY = sigma_y,plot=F);range(cor(cbind(datas$Xs$X1,datas$Xs$X2),datas$Y))
    Xs <- datas$Xs;X_c <- cbind(Xs$X1,Xs$X2);Y <- datas$Y
    S_common <- cbind(datas$S$SXY$S1,datas$S$SXY$S2)
    p <- sum(unlist(lapply(Xs,ncol)))

    lambdas <- 0
    res <- sparse_PLS_Bootstrap(Xs,Y,paras=lambdas,n_B=B,lowExplainedVariance=0,
                                deflatX=T,NCORES=10,center=T,verbose=T)
    RESULTS[[seed]][[level_legend[4]]] <- res


    # res <- sparse_PLS_Bootstrap(Xs,Y,paras=lambdas,n_B=B,lowExplainedVariance=0,
    #                             deflatX=T,NCORES=22,center=T,verbose=T)
    # RESULTS[[seed]][[level_legend[1]]] <- res
    # # res0 <- sparse_PLS_Bootstrap(Xs,Y,paras=0,n_B=400,lowExplainedVariance=0,deflatX=T,NCORES=22,center=T,verbose=T)
    #
    # resmix_Q2_var <- Q2_boot_sPLS(list(x=X_c),Y,keepXs = kxs,keepYs=kys,whichCriterion = "Q2-Var",
    #                               n_B=B,deflatX=T,NCORES=22,center=T,
    #                               NZV=1e-9,verbose=F,sparse = sparse)
    # RESULTS[[seed]][[level_legend[2]]] <- resmix_Q2_var
    #
    # resmix_Q2 <- Q2_boot_sPLS(list(x=X_c),Y,keepXs = kxs,keepYs=kys,whichCriterion = "Q2",
    #                           n_B=B,deflatX=T,NCORES=22,center=T,
    #                           NZV=1e-9,verbose=F,sparse = sparse)
    # RESULTS[[seed]][[level_legend[3]]] <- resmix_Q2

    # print(Q2_all[seed,1] <- RESULTS[[seed]][[level_legend[1]]]$optimal_parameters$Q2)
    # print(R_all[seed,1] <- RESULTS[[seed]][[level_legend[1]]]$optimal_parameters$R)
    # print(Q2_all[seed,2] <- RESULTS[[seed]][[level_legend[2]]]$optimal_parameters$Q2)
    # print(R_all[seed,2] <- RESULTS[[seed]][[level_legend[2]]]$optimal_parameters$R)
    # print(Q2_all[seed,3] <- RESULTS[[seed]][[level_legend[3]]]$optimal_parameters$Q2)
    # print(R_all[seed,3] <- RESULTS[[seed]][[level_legend[3]]]$optimal_parameters$R)
    # eta_i <- RESULTS[[seed]][[level_legend[5]]]$mo$eta
    # K_i <- RESULTS[[seed]][[level_legend[5]]]$mo$K
    # rr <- RESULTS[[seed]][[level_legend[5]]]$res_all
    # id <- which(rr[,2]==K_i & rr[,3]==eta_i)
    # print(Q2_all[seed,4] <- RESULTS[[seed]][[level_legend[5]]]$res_all[id,1])
    # print(R_all[seed,4] <- K_i)
    # print(Q2_all[seed,5] <- RESULTS[[seed]][[level_legend[6]]]$optimal_parameters$Q2)
    # print(R_all[seed,5] <- RESULTS[[seed]][[level_legend[6]]]$optimal_parameters$R)

    # par(mfrow=c(1,2))
    # boxplot(Q2_all,col=cols[id_meth],lty=1,type="l");abline(h=3/5*(1-sigma_y^2))
    # boxplot(R_all,col=cols[id_meth],lty=1,type="l")


    # cl <- makeCluster(length(etas));registerDoParallel(cl)
    # res_all <- foreach(et=1:length(etas),.packages = "spls",.combine='rbind') %dopar% {
    #   cvmo <- cv.spls(X_c,Y,fold=50,K = Ks,
    #                   eta = etas[et],scale.y = T,plot.it = F)
    #   c(max(1-cvmo$mspemat),cvmo$K.opt,etas[et])
    # };stopCluster(cl)
    # id_max <- which.max(res_all[,1])
    # mo <- spls::spls(X_c,Y,K = res_all[id_max,2],eta = res_all[id_max,3])
    # RESULTS[[seed]][[level_legend[5]]] <- list(mo=mo,res_all=res_all)
    #
    # res_mix <- do_mixomics(list(x=X_c),Y,kxs,kys,22)
    # RESULTS[[seed]][[level_legend[6]]] <- res_mix
    #
    # u1 <- rbind(res$Us[[1]],res$Us[[2]])
    # if(is.null(resmix_Q2_var$Us)){
    #   u2 <- matrix(0,p,1)
    # }else{
    #   u2 <- do.call(rbind,resmix_Q2_var$Us)
    # }
    # if(is.null(resmix_Q2$Us)){
    #   u3 <- matrix(0,p,1)
    # }else{
    #   u3 <- do.call(rbind,resmix_Q2$Us)
    # }
    # Us <- list(u1,u2,u3,
    #            do.call(cbind,lapply(mo$betamat,function(b){rowSums(abs(b))})),
    #            res_mix$U)
    # RESULTS[[seed]]$Us <- Us
    #
    # all_stat <- lapply(Us,function(uu){get_tpr_fpr_auc(p,U = uu,S = S_common)})
    # RESULTS[[seed]]$all_stat <- all_stat
    #
    # # save(RESULTS,file = name_file)
    #
    # tpr <- do.call(cbind,lapply(all_stat,function(aa){aa[[length(aa)]]$tpr}))
    # fpr <- do.call(cbind,lapply(all_stat,function(aa){aa[[length(aa)]]$fpr}))
    # auc <- c(lapply(all_stat,function(aa){unlist(lapply(aa,function(bb){bb$auc}))}))
    # R_max <- max(unlist(lapply(auc,length)))
    # auc_mat <- do.call(cbind,lapply(auc,function(aa){c(aa,rep(NA,R_max-length(aa)))}))
    #
    # pdf(file = paste(name_pdf,"_",seed,".pdf",sep = ""),width = 8,height = 6)
    # layout(mat = matrix(c(1,1,1,1,1,1,1,
    #                       1,1,1,1,1,1,1,
    #                       1,1,2,2,2,2,1,
    #                       1,1,2,2,2,2,1,
    #                       1,1,2,2,2,2,1,
    #                       1,1,1,1,1,1,1),nrow = 6,byrow = T))
    # matplot(fpr,tpr,xlim=c(-0,1),ylim=c(-0,1),type="l",xlab="FPR",ylab="TPR",
    #         main="ROC curves for the built models",bty="n",col=cols[id_meth])
    # legend(x=0.8,y=0.9 ,legend = level_legend[id_meth],lty=1:length(id_meth),
    #        col=cols[id_meth],ncol = 1,pt.cex=2,pch=1:length(id_meth))
    # matplot(auc_mat,type="l",bty="n",main="AUC of the ROC curves",
    #         ylab="AUC",xlab="Component",col=cols[id_meth])
    # abline(h=(0:10)/10,col="gray90");abline(v=1:7,col="gray90")
    # matplot(auc_mat,pch=1:length(id_meth),add=T,type = "p",col=cols[id_meth])#;box("figure")
    # dev.off()
  }

  aucs=Rs=NB_VARS=Q2s <- matrix(F,100,6)
  isSparse <- matrix(1,100,6)
  fprs=tprs <- list()

  for(seed in 1:nrow(aucs)){
    datas <- simulateMulti(seed=seed,sigma1 = 0.05,sigma2 = 0.05,sigmaY = 0.2,plot=F)
    Xs <- datas$Xs;X_c <- cbind(Xs$X1,Xs$X2);Y <- datas$Y
    S_common <- cbind(datas$S$SXY$S1,datas$S$SXY$S2)
    p <- sum(unlist(lapply(Xs,ncol)))
    oo_id <- c(1,2,3,4,5,6)
    for(i_id in 1:6){
      id <- oo_id[i_id]
      res_i <- RESULTS[[seed]][[level_legend[id]]]
      NB_VARS[seed,i_id] <-
        if(id==1){
          sum((rowSums(abs(rbind(res_i$Us[[1]],res_i$Us[[2]]) ))>1e-9)*1)
        }else if(id %in% c(2,3)){
          if(!is.null(res_i$Us[[1]])){
            sum((rowSums(abs(res_i$Us[[1]]))>1e-9)*1)
          }else{
            0
          }
        }else if(id==5){
          sum((rowSums(abs(res_i$mo$betahat))>1e-9)*1)
        }else if(id==4){
          p
        }else if(id==6){
          sum((rowSums(abs(res_i$U))>1e-9)*1)
        }

      Q2s[seed,i_id] <-
        if(id %in% c(1,2,3,4)){
          if(!is.null(res_i$optimal_parameters$Q2)){
            res_i$optimal_parameters$Q2
          }else{
            0
          }
        }else if(id==5){
          max(res_i$res_all[,1])
        }else if(id==6){
          res_i$optimal_parameters$Q2
        }

      Rs[seed,i_id] <- if(id %in% c(1,2,3,4)){
        if(!is.null(res_i$Us[[1]])){
          res_i$optimal_parameters$R
        }else{
          0
        }
      }else if(id==5){
        length(res_i$mo$betamat)
      }else if(id==6){
        res_i$optimal_parameters$R
      }
    }

    RESULTS[[seed]]$Us -> Us
    all_stat <- lapply(Us,function(uu){get_tpr_fpr_auc(p,U = uu,S = S_common)})

    tpr <- do.call(cbind,lapply(all_stat,function(aa){aa[[length(aa)]]$tpr}))
    fpr <- do.call(cbind,lapply(all_stat,function(aa){aa[[length(aa)]]$fpr}))
    tprs[[seed]] <- tpr
    fprs[[seed]] <- fpr
    auc <- c(lapply(all_stat,function(aa){unlist(lapply(aa,function(bb){bb$auc}))}))
    aucs[seed,c(1,2,3,5,6)] <- unlist(lapply(auc,function(aa){aa[length(aa)]}))
    isSparse[seed,which(is.na(aucs[seed,]))] <- 0
    R_max <- max(unlist(lapply(auc,length)))
    auc_mat <- do.call(cbind,lapply(auc,function(aa){c(aa,rep(NA,R_max-length(aa)))}))
    # matplot(fpr,tpr,xlim=c(-0,1),ylim=c(-0,1),type="l",xlab="FPR",ylab="TPR",
    #         main="ROC curves for the built models",bty="n",col=cols[id_meth],add=T)
  }

  range_r <- range(Rs)
  mat_R_s <- matrix(0,max(Rs),6)
  id_meth_R <- 1:6
  colnames(mat_R_s) <- level_legend[id_meth_R]
  for(i_r in 1:max(Rs)){
    for(i_m in 1:length(id_meth_R)){
      mat_R_s[i_r,i_m] <- sum((Rs[,i_m]==i_r)*1)
    }
  }

  ###########################
  ## PLOT Number Variables ##
  ###########################
  postscript("/Users/hlorenzo/Documents/GitHub/Simulations_design3_var.eps", onefile=TRUE, horizontal=F,
             width=6, height=1.5,pointsize=1)
  par(mfrow=c(1,5))
  meth_id <- c(1,2,3,5,6)
  n_bins <- 25
  breaks=unique(round(seq(min(NB_VARS),max(NB_VARS),length.out = n_bins)))
  par(mar=c(2,3.5,2,1))
  mama <- 35
  for(ii in 1:5){
    h=hist(NB_VARS[,ii],breaks=breaks,plot=F)
    ij <-(meth_id)[ii]
    plot(-10,-1,ylim=c(0,n_bins),xlim=c(0,37),bty="n",main="",yaxt="n",xaxt="n",xlab="Count",ylab="")
    # abline(v=seq(0,mama,by = 5),lty=2,lwd=1/2,col="gray70")
    abline(v=(0:10)*10+5,lty=2,col="gray90",lwd=0.7)
    abline(v=(0:10)*10,lty=2,col="gray65",lwd=0.7)
    abline(v=c(0,25,50,75,100),lty=1,col="gray50",lwd=1.1)
    bbi <- rev(h$breaks)
    print(range(h$counts))
    df_hh <- matrix(rev(h$counts),ncol = 1);rownames(df_hh) <- (bbi[-1]+bbi[-length(bbi)])/2
    barplot(df_hh,horiz=T,col=cols[ij],lwd=0.7,border=cols[ij],angle=angle[ij]+0,density=density[ij],add=T,beside=T)
    title(main = level_legend[meth_id[ii]],line = 0)
    axis(1,at=seq(0,mama,by = 5),labels = seq(0,mama,by = 5))
    axis(2,at=1:length(bbi),labels = bbi,las=2)
  }
  dev.off()
  ############################
  ## PLOT Number Components ##
  ############################
  bb <- mat_R_s*0
  bb[3,] <- max(mat_R_s)
  rownames(bb) <- 1:max(Rs)
  cex.leg=1/2
  postscript("/Users/hlorenzo/Documents/GitHub/Simulations_design3_R.eps", onefile=TRUE, horizontal=F,
             width=6, height=1.5,pointsize=1)
  layout(matrix(c(rep(1,9),rep(1+1,1)),byrow = T,ncol = 1))
  # par(mar=c(2,3,2.5,1.5),cex=1)
  par(mar=c(2,3.5,2,1))
  barplot(t(bb),beside = TRUE,col="gray90",border = "gray90")#,
  # main="Number of components")
  title(main = "Number of components",line = 0)
  abline(h=(0:10)*10+5,lty=2,col="gray90",lwd=0.7)
  abline(h=(0:10)*10,lty=2,col="gray65",lwd=0.7)
  abline(h=c(0,25,50,75,100),lty=1,col="gray50",lwd=1.1)
  barplot(t(mat_R_s),beside = TRUE,density = density[id_meth_R],
          angle=angle[id_meth_R],col=cols[id_meth_R],
          border=cols[id_meth_R],add=T)
  abline(v=(length(id_meth_R)+1)*(0:(max(Rs)+1))+0.5,lty=1,col="gray65",lwd=0.5)
  par(mar=rep(0,4))
  plot(0,type='n',axes=FALSE,ann=FALSE)
  legend("center",legend = level_legend[id_meth_R],angle=angle[id_meth_R],
         border=cols[id_meth_R],fill=cols[id_meth_R],density=density[id_meth_R],
         col = cols[id_meth_R],ncol = length(id_meth_R),bty = "n",cex = 2*cex.leg,pt.cex=2)
  dev.off()
  ############################ END

  fprs_met <- lapply(1:ncol(fprs[[1]]),function(ii){do.call(cbind,lapply(fprs,function(xx){xx[,ii]}))})
  tprs_met <- lapply(1:ncol(tprs[[1]]),function(ii){do.call(cbind,lapply(tprs,function(xx){xx[,ii]}))})
  postscript("/Users/hlorenzo/Documents/GitHub/Simulations_design3_auc.eps", onefile=TRUE, horizontal=F,
             width=10, height=5,pointsize=1)
  layout(matrix(c(rep(1:5,4),rep(c(6,6,7,7,7),4*2),rep(8,5)),ncol = 5,byrow = T))
  par(mar=c(2,2,2,1),cex=1.1)
  tpr_mean=fpr_mean <- NULL
  for(id in 1:length(id_meth)){
    plot(2,2,xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n",bty='L',
         xlab="",ylab="",main=paste("ROC curve of",level_legend[id_meth[id]]))
    axis(1,at = c(0,0.5,1),labels = c("0","1/2","1"),line = -1,tick = F)
    axis(2,at = c(0,0.5,1),labels = c("0","1/2","1"),line = 0,tick = F,las=2)
    abline(0,1,lty=2)
    abline(h=(c(0:10)/10)[-c(1,6,11)],v=(c(0:10)/10)[-c(1,6,11)],col="gray80",lty=2,lwd=1/2)
    abline(h=c(0,1,1/2),v=c(0,1,1/2),col="gray80",lty=1,lwd=1/2)
    matplot(fprs_met[[id]],tprs_met[[id]],type="l",col=cols[id_meth[id]],lty = 1,add=T)
    mtext("FPR", side=1, line=-1,at=1)
    mtext("TPR", side=2, line=-2,las=2,at = 1)
    points(rowMeans(fprs_met[[id]],na.rm = T),rowMeans(tprs_met[[id]],na.rm = T),
           type="l",lty=1,lwd=3,col=col_box)
    points(rowMeans(fprs_met[[id]],na.rm = T),rowMeans(tprs_met[[id]],na.rm = T),
           type="l",lty=1,lwd=1.5,col=1)
    tpr_mean <- cbind(tpr_mean,rowMeans(tprs_met[[id]],na.rm = T))
    fpr_mean <- cbind(fpr_mean,rowMeans(fprs_met[[id]],na.rm = T))
  }
  plot(2,2,xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n",bty='L',
       xlab="",ylab="",main="Mean ROC curves")
  # matplot(fpr_mean,tpr_mean,type="l",col=cols[id_meth],pch=id_meth,lwd=2,
  #         xlab="FPR",ylab="TPR",main="Mean ROC curves")
  abline(0,1,lty=2)
  abline(h=(c(0:10)/10)[-c(1,6,11)],v=(c(0:10)/10)[-c(1,6,11)],col="gray80",lty=2,lwd=1/2)
  abline(h=c(0,1,1/2),v=c(0,1,1/2),col="gray80",lty=1,lwd=1/2)
  matplot(fpr_mean,tpr_mean,type="l",col=cols[id_meth],lty = 1,add=T)
  axis(1,at = c(0,0.5,1),labels = c("0","1/2","1"),line = -1,tick = F)
  mtext("FPR", side=1, line=-2,at=1)
  axis(2,at = c(0,0.5,1),labels = c("0","1/2","1"),line = 0,tick = F,las=2)
  mtext("TPR", side=2, line=-2,las=2.5,at = 1)
  id_pch <- seq(1,nrow(fpr_mean),length.out = 6)
  matplot(fpr_mean[id_pch,],tpr_mean[id_pch,],pch=1:length(level_legend[id_meth]),
          lty=1:length(level_legend[id_meth]),col = cols[id_meth],add=T,cex=2)
  abline(0,1,lty=2)
  colnames(aucs) <- level_legend[id_meth]
  b <- boxplot(aucs,border=cols[id_meth],col=col_box,
               main="AUC of variable selection",
               ylab="",pch=(1:length(level_legend))[id_meth])
  abline(h=(0:20)/20,col="gray60",lty=2)
  coli_s <- cols[id_meth]
  for(ii in 1:length(coli_s)){
    rect(xleft = ii-0.5+0.1,xright = ii+0.5-0.1,
         ytop = b$stats[4,ii],ybottom = b$stats[2,ii],
         density = density[id_meth[ii]],angle=angle[id_meth[ii]],col=coli_s[ii])
  }
  par(mar=rep(0,4))
  plot(0,type='n',axes=FALSE,ann=FALSE)
  legend("center",legend = level_legend[id_meth],lty=1,
         border = cols[id_meth],col = cols[id_meth],fill = cols[id_meth],
         density = density[id_meth],angle=angle[id_meth],
         pch=1:length(level_legend[id_meth]),
         ncol = 6,bty = "n",pt.cex=2,lwd=2)
  dev.off()

  postscript("/Users/hlorenzo/Documents/GitHub/Simulations_design3_auc_short.eps", onefile=TRUE, horizontal=F,
             width=10, height=1.5)
  # layout(matrix(c(rep(1,9),rep(2,1)),ncol = 1,byrow = T))
  par(mar=c(2,0.1,2,0.1),cex=0.7)
  colnames(aucs) <- level_legend[id_meth]
  b <- boxplot(aucs,border=cols[id_meth],col=col_box,ylim=c(0.47,1),
               main="AUC",xlab="",xaxt="n",horizontal=T,
               ylab="",yaxt="n",pch=(1:length(level_legend))[id_meth])
  abline(v=(0:10)/10+0.05,lty=2,col="gray90",lwd=0.7)
  abline(v=(0:10)/10,lty=2,col="gray65",lwd=0.7)
  abline(v=c(75,100)/100,lty=1,col="gray50",lwd=1.1)
  axis(1,at=seq(0,1,by = 0.1),labels = seq(0,1,by = 0.1))
  coli_s <- cols[id_meth]
  for(ii in 1:length(coli_s)){
    rect(ybottom = ii-0.5+0.1,ytop = ii+0.5-0.1,
         xright = b$stats[4,ii],xleft = b$stats[2,ii],
         density = density[id_meth[ii]],angle=angle[id_meth[ii]],col=coli_s[ii])
    legend(x = 0.47,y = ii+0.5,legend = level_legend[id_meth[ii]],
           border = cols[id_meth[ii]],col = cols[id_meth[ii]],fill = cols[id_meth[ii]],
           density = density[id_meth[ii]],angle=angle[id_meth[ii]],
           pch=ii,
           ncol = 1,bty = "n",cex=0.7)
  }
  dev.off()
}



plot_R2_Q2 <- function(){


  layout(matrix(c(rep(c(1:5),10),rep(6,5)),ncol=5,byrow = T))
  # par(mar=c(3.1, 3.1, 4.1, 0.2))
  cols <- RColorBrewer::brewer.pal(6+1,"Set1")[-6]
  col_box <- RColorBrewer::brewer.pal(9,"Pastel1")[7]

  ylim <- c(-0.05,1.4*1)
  lwd <- 1.5
  ncols <- 3
  cex.leg <- 1
  max_width <- max(apply(Rs,2,function(rr){max(table(rr))}))
  R_max <- 12
  tab_mean_sd <- matrix(NA,5,R_max)
  for(id in 1:5){
    id_meth_i <- c(1,2,3,5,6)[id]
    # width_i <- table(Rs[,id])/max_width
    dff <- data.frame(cbind(R=Rs[,id],Q2=Q2s[,id]))
    for(rr in 1:R_max){
      idid <- which(Rs[,id]==rr)
      if(length(idid)>0){
        # tab_mean_sd[id,rr] <- paste(round(mean(Q2s[idid,id]),3)," (\\pm",round(sd(Q2s[idid,id]),3),")",sep="")
        tab_mean_sd[id,rr] <- round(median(Q2s[idid,id]),3)
      }
    }
    # b<-boxplot(Q2~R,border=cols[id_meth_i],col="white",
    #            dff,xlim=c(1,12),ylim=c(0,1),
    #            xlab="R",ylab="",bty="l",lwd=lwd,pch=id_meth_i,width=width_i)
    # abline(h=3/5*(1-sigma_y^2 ),lty=2)
    # popo <- 1.2
    # for(ii in 1:length(width_i)){
    #   rect(xleft = ii-popo*width_i[ii]/2,xright = ii+popo*width_i[ii]/2,
    #        ytop = b$stats[4,ii],ybottom = b$stats[2,ii],
    #        density = density[id_meth_i],angle=angle[id_meth_i],col=coli_s[id_meth_i])
    # }
  }
  rownames(tab_mean_sd) <- level_legend[c(1,2,3,5,6)]

  par(mar=rep(0,4))
  plot(0,type='n',axes=FALSE,ann=FALSE)
  legend("top",legend = level_legend,angle=angle,border=cols,fill=cols,density=density,col = cols,pch=1:length(level_legend),
         ncol = 6,bty = "n",cex = cex.leg,pt.cex=2)
}

plot_tpr_fpr_roc <- function(alpha,S_com,id,RESULTS,N_simu=20,
                             level_legend){
  get_U <- function(id,seed,RESULTS,level_legend){
    res_i <- RESULTS[[seed]][[level_legend[id]]]
    if(id==1){
      (rowSums(abs(rbind(res_i$Us[[1]],res_i$Us[[2]]) ))>1e-9)*1
    }else if(id %in% c(2,3)){
      if(!is.null(res_i$Us[[1]])){
        (rowSums(abs(res_i$Us[[1]]))>1e-9)*1
      }else{
        matrix(0,5500,1)
      }
    }else if(id==5){
      (rowSums(abs(res_i$mo$betahat))>1e-9)*1
    }else if(id==6){
      (rowSums(abs(res_i$U))>1e-9)*1
    }
  }
  my_tpr_fpr_roc <- function(sel,S_com,alpha){
    nn <- length(alpha)
    vect_th <- apply(S_com,2,max)
    fpr = tpr <- rep(NA,nn)
    for(i in 1:nn){
      sel_th <- (vect_th>alpha[i])*1
      tpr[i] <- sum((sel_th==1)*(sel==1))/sum(sel==1)#tabi[2,2]/sum(tabi[,2])#
      fpr[i] <- sum((sel_th==1)*(sel==0))/sum(sel==0)#tabi[2,1]/sum(tabi[,1])#
    }
    aa <- abs((fpr[-1]-fpr[-nn])*(tpr[-1]+tpr[-nn])/2)
    auc <- sum(na.omit(aa))
    list(tpr=tpr,fpr=fpr,auc=auc)
  }
  library(ggplot2)

  nn <- length(alpha)
  tprs = fprs <- rep(NA,N_simu*nn)
  ii <- 1
  for(seed in 1:N_simu){
    sel <- get_U(id,seed,RESULTS,level_legend)
    datas <- simulateMulti(seed=seed,sigma1 = 0.05,sigma2 = 0.05,sigmaY = sigma_y,plot=F)
    S_common <- cbind(datas$S$SXY$S1,datas$S$SXY$S2)
    S_com <- S_common#/max(S_common)
    result <- my_tpr_fpr_roc(sel,S_com,alpha)
    tprs[(seed-1)*nn+1:nn] <- result$tpr
    fprs[(seed-1)*nn+1:nn] <- result$fpr
    cat(" ");cat(seed)
    if(is.na(tprs[1])) browser()
  }
  all_alphas <- rep(alpha,N_simu)
  df_FPR <- as.data.frame(cbind(all_alphas,fprs))
  names(df_FPR) <- c("alpha","FPR")
  df_TPR <- as.data.frame(cbind(all_alphas,tprs))
  names(df_TPR) <- c("alpha","TPR")
  df_ROC <- as.data.frame(cbind(all_alphas,tprs,fprs))
  names(df_ROC) <- c("alpha","TPR","FPR")


  # library(gridExtra)
  # grid.arrange(tp, fp, roc, ncol=3, nrow = 1)
  list(df_FPR=df_FPR,df_TPR=df_TPR,df_ROC=df_ROC)
}



plot_all_bins <- function(all_matrices,all_id,bins,bins_alpha=20,
                          dens_min1 = 0,dens_min2 = 0,
                          dens_max1= 3000,dens_max2= 4000,
                          cex.txt=4,
                          col_name="viridis"){
  all_plots <- list()
  # all_alphas <- all_matrices[[1]]$df_TPR[,1]
  # alphas <- unique(all_alphas)
  # means_tpr_fpr <- matrix(NA,length(alphas),2)
  all_plots <- list()
  # all_alphas <- all_matrices[[1]]$df_TPR[,1]
  # alphas <- unique(all_alphas)
  # means_tpr_fpr <- matrix(NA,length(alphas),2)
  Count <- c(1,2)
  shots <- rep('Label',2)
  xc <-c(1:2)
  yc <-c(1:2)
  df00 <-data.frame(shots,Count,xc,yc)
  p1 <- ggplot() +
    geom_point(data = df00, aes(x = xc, y = yc, fill = Count), shape=22, size=0,color = 'white', stroke=1) +
    ylim(0,1)+
    #Color and Legend
    scale_fill_continuous(type = col_name,limits = c(dens_min1, dens_max1),na.value ="red") +
    #Theme
    theme(text = element_text(size=cex.txt),legend.title=element_text(size=cex.txt),
          legend.text=element_text(size=cex.txt),panel.background = element_rect(fill = "transparent",colour = NA),plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),plot.title = element_text(size = 14, hjust = 1, vjust = 1),plot.background = element_rect(fill = "transparent", colour = NA),      axis.title=element_blank(),      axis.text = element_blank(),      axis.ticks = element_blank(),      legend.position = c(.9,.75),           legend.background = element_rect(fill = "transparent")    )+ theme(legend.position = "right")
  p2 <- ggplot() +
    geom_point(data = df00, aes(x = xc, y = yc, fill = Count), shape=22, size=0,color = 'white', stroke=1) +
    ylim(0,1)+
    #Color and Legend
    scale_fill_continuous(type = col_name,limits = c(dens_min1, dens_max1),na.value ="red") +
    #Theme
    theme(text = element_text(size=cex.txt),legend.title=element_text(size=cex.txt),
          legend.text=element_text(size=cex.txt),panel.background = element_rect(fill = "transparent",colour = NA),plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),plot.title = element_text(size = 14, hjust = 1, vjust = 1),plot.background = element_rect(fill = "transparent", colour = NA),      axis.title=element_blank(),      axis.text = element_blank(),      axis.ticks = element_blank(),      legend.position = c(.9,.75),           legend.background = element_rect(fill = "transparent")    )+ theme(legend.position = "right")

  p3 <- ggplot() +
    geom_point(data = df00, aes(x = xc, y = yc, fill = Count), shape=22, size=0,color = 'white', stroke=1) +
    ylim(0,1)+xlim(0,1)+
    #Color and Legend
    scale_fill_continuous(type = col_name,limits = c(dens_min2, dens_max2),na.value ="red") +
    #Theme
    theme(text = element_text(size=cex.txt),legend.title=element_text(size=cex.txt),
          legend.text=element_text(size=cex.txt),panel.background = element_rect(fill = "transparent",colour = NA),plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),plot.title = element_text(size = 14, hjust = 1, vjust = 1),plot.background = element_rect(fill = "transparent", colour = NA),      axis.title=element_blank(),      axis.text = element_blank(),      axis.ticks = element_blank(),      legend.position = c(.9,.75),        legend.background = element_rect(fill = "transparent")    )+ theme(legend.position = "right")


  all_plots[[1]] <- p1
  all_plots[[2]] <- p2
  all_plots[[3]] <- p3
  count <- length(all_plots)+1
  for(i_id in 1:length(all_id)){
    id <- all_id[i_id]
    # youou_tpr <- all_matrices[[id]]$df_TPR
    # youou_fpr <- all_matrices[[id]]$df_FPR
    # for(ii in 1:length(alphas)){
    #   pos_a <- which(all_alphas==alphas[ii])
    #   means_tpr_fpr[ii,] <- c(
    #     mean(youou_tpr[pos_a,2],na.rm = T),
    #     mean(youou_fpr[pos_a,2],na.rm = T)
    #   )
    # }
    # means_tpr_fpr <- data.frame(list(alpha=alphas,TPR=means_tpr_fpr[,1],FPR=means_tpr_fpr[,2]))
    legn.pos <- "none"
    all_plots[[count]] <- ggplot(all_matrices[[id]]$df_TPR, aes(x=alpha, y=TPR)) +
      geom_hex(bins = c(bins_alpha,bins[i_id,1]))+xlab(expression(alpha~"\n")) +ylab(expression("\n"~TPR)) +
      scale_fill_continuous(type = col_name,limits = c(dens_min1, dens_max1),na.value ="red")+#geom_line(data=means_tpr_fpr,aes(x=alpha, y=TPR),color="red")+
      ggtitle(bquote("\n TPR versus"~alpha~"for"~.(level_legend[id])))+
      theme(legend.position = legn.pos,text = element_text(size=cex.txt))
    # theme_bw()
    count <- count + 1
    all_plots[[count]] <- ggplot(all_matrices[[id]]$df_FPR, aes(x=alpha, y=FPR)) +
      geom_hex(bins = c(bins_alpha,bins[i_id,2]))+xlab(expression(alpha~"\n")) +ylab(expression("\n"~FPR)) +
      scale_fill_continuous(type = col_name,limits = c(dens_min1, dens_max1),na.value ="red")+
      ggtitle(bquote("\n FPR versus"~alpha~"for"~.(level_legend[id])))  + theme(legend.position = legn.pos,text = element_text(size=cex.txt))
    # theme_bw()
    count <- count + 1
    all_plots[[count]] <- ggplot(all_matrices[[id]]$df_ROC, aes(x=FPR, y=TPR)) +
      geom_hex(bins = c(bins[i_id,3],bins[i_id,3])) +xlab(expression(FPR~"\n")) +ylab(expression("\n"~FPR))+
      scale_fill_continuous(type = col_name,
                            limits = c(dens_min2, dens_max2),na.value ="red")+
      geom_abline(intercept = 0, color = "black", size = 0.5,lty=2)+
      ggtitle(paste("\n ROC Curve for",level_legend[id]))+
      theme(legend.position = legn.pos,text = element_text(size=cex.txt))
    # theme_bw()
    count <- count + 1
  }
  return(do.call("grid.arrange", c(all_plots, ncol=6,as.table=F)))
}


plot_cool <- function(N_simu=100){
  library(gridExtra)
  oo_alpha<-c(lapply(1:N_simu,function(seed){
    datas <- simulateMulti(seed=seed,sigma1 = 0.05,sigma2 = 0.05,sigmaY = sigma_y,plot=F)
    S_common <- cbind(datas$S$SXY$S1,datas$S$SXY$S2)
    S_com <- S_common#/max(S_common)
    # sel <- sample(x = 0:1,replace = T,size = ncol(S_com))
    # pred <- ROCR::prediction(apply(S_com,2,max),sel)
    # perf <- ROCR::performance(pred,"tpr","fpr")
    # alpha <- perf@alpha.values[[1]]
    apply(S_com,2,max)}))
  alpha_0 <- do.call(c,oo_alpha)
  N_bins <- 50
  hh <- hist(alpha_0,breaks = N_bins)
  N_tot <- 5500/N_bins
  coco_norm <- ceiling(hh$counts/N_tot)
  alpha <- unlist(lapply(1:length(coco_norm),function(ii){
    set.seed(101)
    runif(min = hh$breaks[ii],hh$breaks[ii+1],n = coco_norm[ii])
  }))
  id_alpha <- seq(1,length(alpha),by = 10)
  # alpha <- alpha[id_alpha]
  # alpha <- rnorm(1000)^2
  # alpha <- (exp(seq(0,1,length.out = 600))-1)/(exp(1)-1)
  # alpha2 <- (exp(runif(0,0.1,n = 100))-1)/(exp(1)-1)
  # alpha <- unique(sort(c(alpha,alpha2)))
  # alpha <- seq(0,1,length.out = 10000)
  all_matrices <- list()
  all_id <- c(1,2,3,5,6)
  for(id in all_id){
    cat("---------------------------------------")
    cat(id);cat(" ")
    all_matrices[[id]] <- plot_tpr_fpr_roc(alpha = alpha,S_com = S_com,
                                           id = id,RESULTS = RESULTS,
                                           level_legend = level_legend,
                                           N_simu = N_simu)
  }
  all_matrices_backup <- all_matrices
  library(viridis)
  bins <- matrix(25,5,3)
  bins[,1]=bins[,2] <- 100
  postscript("/Users/hlorenzo/Documents/GitHub/Simulations_design3_hexbin_01000_01000_02000.eps", onefile=TRUE, horizontal=F,
             width=8, height=4.5,pointsize=0.1)
  plot_all_bins(all_matrices,all_id,bins,col_name="viridis",
                dens_min1 = 0,dens_min2 = 0,
                dens_max1 = 4000,dens_max2 = 3000,
                bins_alpha = 30,cex.txt=3.5)
  dev.off()



  postscript("/Users/hlorenzo/Documents/GitHub/Simulations_design3_auc_short.eps", onefile=TRUE, horizontal=F,
             width=10, height=1.5)
  par(mar=c(2,0.1,2,0.1),cex=0.7)
  aucs2 <- aucs[,c(1,2,3,5,6)]
  colnames(aucs2) <- level_legend[id_meth]
  b <- boxplot(aucs2,border=cols[id_meth],col=col_box,ylim=c(0.47,1),
               main="AUC",xlab="",yaxt="n",horizontal=T,
               ylab="",pch=(1:length(level_legend))[id_meth])
  abline(v=(0:20)/20,col="gray60",lty=2)
  coli_s <- cols[id_meth]
  for(ii in 1:length(coli_s)){
    rect(ybottom = ii-0.5+0.1,ytop = ii+0.5-0.1,
         xright = b$stats[4,ii],xleft = b$stats[2,ii],
         density = density[id_meth[ii]],angle=angle[id_meth[ii]],col=coli_s[ii])
    legend(x = 0.47,y = ii+0.5,legend = level_legend[id_meth[ii]],
           border = cols[id_meth[ii]],col = cols[id_meth[ii]],fill = cols[id_meth[ii]],
           density = density[id_meth[ii]],angle=angle[id_meth[ii]],
           pch=ii,
           ncol = 1,bty = "n",cex=0.8)
  }
  dev.off()

  # postscript("/Users/hlorenzo/Documents/GitHub/Simulations_design3_hexbin_3000115000_5008500_20005000.eps", onefile=TRUE, horizontal=F,
  #            width=8, height=4.5,pointsize=0.1)
  # plot_all_bins(all_matrices,all_id,bins,col_name="viridis",
  #               dens_min1 = 3000,dens_min2 = 500,
  #               dens_max1 = 115000,dens_max2 = 8500,
  #               bins_alpha = 25,cex.txt=3.5)
  # dev.off()

  #########
  ## Histo extreme FPR
  #########

  postscript("/Users/hlorenzo/Documents/GitHub/hist_oo.eps",
             onefile=TRUE, horizontal=F,
             width=10, height=1.5)
  thre <- 0.9
  count_max <- 20000
  par(mfrow=c(1,5),cex=0.3,mar=c(2,3,2,0.5))
  for(id_m in c(1,2,3,5,6)){
    if(id_m>1){
      par(mar=c(2,0.5,2,0.5))
    }
    id_TPR_ok <- which(all_matrices[[id_m]]$df_TPR$TPR>thre)
    plot(-10,-10,xlim=c(0,1),ylim=c(0,count_max),main=level_legend[id_m],xaxt="n",yaxt="n",bty="n")
    abline(v=seq(0,1,length.out = 6),col="gray80",lty=1)
    abline(h=seq(0,count_max,length.out = 6),col="gray80",lty=1)
    hist(all_matrices[[id_m]]$df_FPR$FPR[id_TPR_ok],35,xlim = c(0,1),ylim = c(0,count_max),
         xlab="",xaxt="n",ylab="",yaxt="n",col=cols[id_m],bty="n",add=T)#,angle=angle[id_m],density=density[id_m])
    if(id_m==1){
      axis(2,at = seq(0,count_max,length.out = 6),labels = paste(seq(0,count_max,length.out = 6)/1000,"k"),line = 0,las=2)
      text(0.1,y = 19000,labels = "Frequency")
    }
    text(0.96,y = 600,labels = "FPR")
    axis(1,at = seq(0,1,length.out = 6),labels = seq(0,1,length.out = 6),line = 0)
    abline(v=mean(all_matrices[[id_m]]$df_FPR$FPR[id_TPR_ok],na.rm = T)*c(1,1),lwd=c(2,1),lty=c(1,2),col=c(1,cols[id_m]))
  }
  dev.off()

  ## HISTO et AUC

  postscript("/Users/hlorenzo/Documents/GitHub/hist_AUC.eps",
             onefile=TRUE, horizontal=F,
             width=4, height=4.5)
  thre <- 0.9
  count_max <- 20000
  layout(matrix(c(6,6,6,6,6,6,1,1,2,2,3,3,7,4,4,5,5,7),ncol = 6,byrow = T))
  par(cex=0.3,mar=c(2,3,2,0.5))
  for(id_m in c(1,2,3,5,6)){
    if(id_m<1){
      par(mar=c(2,0.5,2,0.5))
    }
    id_TPR_ok <- which(all_matrices[[id_m]]$df_TPR$TPR>thre)
    plot(-10,-10,xlim=c(0,1),ylim=c(0,count_max),
         main=bquote("FPR"["TPR>"~.(thre)]~"for"~.(level_legend[id_m])),xaxt="n",yaxt="n",bty="n")
    abline(v=seq(0,1,length.out = 6),col="gray80",lty=1)
    abline(h=seq(0,count_max,length.out = 6),col="gray80",lty=1)
    hist(all_matrices[[id_m]]$df_FPR$FPR[id_TPR_ok],lwd=0.01,35,xlim = c(0,1),ylim = c(0,count_max),
         xlab="",xaxt="n",ylab="",yaxt="n",col=cols[id_m],bty="n",add=T)#,angle=angle[id_m],density=density[id_m])
    if(T){
      axis(2,at = seq(0,count_max,length.out = 6),labels = paste(seq(0,count_max,length.out = 6)/1000,"k"),line = 0,las=2)
      text(0.1,y = 19000,labels = "Frequency")
    }
    text(0.96,y = 600,labels = "FPR")
    axis(1,at = seq(0,1,length.out = 6),labels = seq(0,1,length.out = 6),line = 0)
    abline(v=mean(all_matrices[[id_m]]$df_FPR$FPR[id_TPR_ok],na.rm = T)*c(1,1),lwd=c(2,1),lty=c(1,2),col=c(1,cols[id_m]))
  }
  par(mar=c(2,0.1,2,0.1),cex=0.7)
  aucs2 <- aucs[,-4]
  colnames(aucs2) <- level_legend[id_meth]
  b <- boxplot(aucs2,border=cols[id_meth],col=col_box,ylim=c(0.47,1),
               main="AUC",xlab="",xaxt="n",horizontal=T,
               ylab="",yaxt="n",pch=(1:length(level_legend))[id_meth])
  abline(v=(0:10)/10+0.05,lty=2,col="gray90",lwd=0.7)
  abline(v=(0:10)/10,lty=2,col="gray65",lwd=0.7)
  abline(v=c(75,100)/100,lty=1,col="gray50",lwd=1.1)
  axis(1,at=seq(0,1,by = 0.1),labels = seq(0,1,by = 0.1))
  coli_s <- cols[id_meth]
  for(ii in 1:length(coli_s)){
    rect(ybottom = ii-0.5+0.1,ytop = ii+0.5-0.1,
         xright = b$stats[4,ii],xleft = b$stats[2,ii],
         density = density[id_meth[ii]],angle=angle[id_meth[ii]],col=coli_s[ii])
    legend(x = 0.47,y = ii+0.5,legend = level_legend[id_meth[ii]],
           border = cols[id_meth[ii]],col = cols[id_meth[ii]],fill = cols[id_meth[ii]],
           density = density[id_meth[ii]],angle=angle[id_meth[ii]],
           pch=ii,
           ncol = 1,bty = "n",cex=0.7)
  }
  plot(-1,-1,xlim=c(0,1),ylim=c(0,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n")
  dev.off()
}


