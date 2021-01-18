# import numpy as np
# import pyPLS
# from matplotlib import pyplot as plt
# np.random.seed(5)
# X, S, Y = pyPLS.simulateData(50, 4, 500, 20, signalToNoise=100, ConcVarFact=0.8, n_min_peaks=5)
# X2, S2, Y2 = pyPLS.simulateData(50, 4, 5000, 30, signalToNoise=100, ConcVarFact=0.8, n_min_peaks=5)
# Yt = np.zeros((50,3))
# Yt[:,0:3] = Y[:,0:3]
# X2 = Yt @ S2[0:3,:] + 4*Y2[:,3:] @ S2[3:,:] + np.random.normal(size=(50, 5000))
# np.savetxt("X1-v2.csv", X, delimiter=",")
# np.savetxt("X2-v2.csv", X2, delimiter=",")
# np.savetxt("S1.csv", S, delimiter=",")
# np.savetxt("S2.csv", S2, delimiter=",")
# np.savetxt("Y-v2.csv", Yt, delimiter=",")

get_tpr_fpr_auc <- function(p,U,S,verbose=T){
  Rs <- ncol(U)
  out <- list()
  for(R in 1:Rs){
    sel <- rep(0,p);sel[which(rowSums(abs(U[,1:R,drop=F]))>1e-9)] <- 1
    if(sum(sel)!=p){
      pred <- prediction(apply(S,2,max),sel)
      perf <- performance(pred,"tpr","fpr")
      fpr <- perf@x.values[[1]]
      tpr <- perf@y.values[[1]]
      perf <- performance(pred,"auc")
      auc <- perf@y.values[[1]]
    }else{
      fpr <- c(0,1)
      tpr <- c(0,1)
      auc <- 1/2
    }
    out[[R]] <- list(fpr=fpr,tpr=tpr,auc=auc)
  }
  if(verbose){
    cat(unlist(lapply(out,function(rr){rr$auc})));cat("\n")
  }
  out
}


X1 <- read.csv("../Sartorius1/data/newSartorius/X1-v2.csv",header = F)
X2 <- read.csv("../Sartorius1/data/newSartorius/X2-v2.csv",header = F)
Y <- read.csv("../Sartorius1/data/newSartorius/Y-v2.csv",header = F)
S1 <- read.csv("../Sartorius1/data/newSartorius/S1.csv",header = F)
S2 <- read.csv("../Sartorius1/data/newSartorius/S2.csv",header = F)

S <- cbind(S1,S2)

Xs <- list(X1,X2)
p <- ncol(X1)+ncol(X2)

lambdas <- seq(0,1,length.out = 100)
res <- sparse_PLS_Bootstrap(Xs,Y,paras=lambdas,n_B=1000,lowExplainedVariance=0,deflatX=T,NCORES=22,center=T,verbose=T)
res0 <- sparse_PLS_Bootstrap(Xs,Y,paras=0,n_B=400,lowExplainedVariance=0,deflatX=T,NCORES=22,center=T,verbose=T)

etas <- seq(0.1,0.9,0.1);Ks <- 1:8;fold <- min(n,50)
cl <- makeCluster(length(etas));registerDoParallel(cl)
res_all <- foreach(et=1:length(etas),.packages = "spls",.combine='rbind') %dopar% {
  cvmo <- cv.spls(cbind(X1,X2),Y,fold=50,K = Ks,eta = etas[et],scale.y = T,plot.it = F)
  c(max(1-cvmo$mspemat),cvmo$K.opt,etas[et])
};stopCluster(cl)
id_max <- which.max(res_all[,1])
mo <- spls::spls(cbind(X1,X2),Y,K = res_all[id_max,2],eta = res_all[id_max,3])

kxs <- unique(round(seq(1,p,length.out = 22)));kys <- c(1,2); sparse <- T
res_mix <- do_mixomics(list(x=cbind(X1,X2)),Y,kxs,kys,22)

resmix_Q2 <- Q2_boot_sPLS(list(x=cbind(X1,X2)),Y,keepXs = kxs,keepYs=kys,whichCriterion = "Q2",
                          n_B=1000,deflatX=T,NCORES=ncores_i,center=T,
                          NZV=1e-9,verbose=T,sparse = sparse)

resmix_Q2_var <- Q2_boot_sPLS(list(x=cbind(X1,X2)),Y,keepXs = kxs,keepYs=kys,whichCriterion = "Q2-Var",
                              n_B=1000,deflatX=T,NCORES=ncores_i,center=T,
                              NZV=1e-9,verbose=T,sparse = sparse)

Us <- list(rbind(res$Us[[1]],res$Us[[2]]),do.call(cbind,lapply(mo$betamat,function(b){rowSums(abs(b))})),
           res_mix$U,do.call(rbind,resmix_Q2$Us),do.call(rbind,resmix_Q2_var$Us))
all_stat <- lapply(Us,function(uu){get_tpr_fpr_auc(p,U = uu,S = S)})

tpr <- do.call(cbind,lapply(all_stat,function(aa){aa[[length(aa)]]$tpr}))
fpr <- do.call(cbind,lapply(all_stat,function(aa){aa[[length(aa)]]$fpr}))
auc <- c(lapply(all_stat,function(aa){unlist(lapply(aa,function(bb){bb$auc}))}))
R_max <- max(unlist(lapply(auc,length)))
auc_mat <- do.call(cbind,lapply(auc,function(aa){c(aa,rep(NA,R_max-length(aa)))}))

layout(mat = matrix(c(1,1,1,1,1,1,1,1,1,1,
                      1,1,1,1,1,1,1,1,1,1,
                      1,1,1,1,1,1,1,1,1,1,
                      1,1,1,1,1,1,1,1,1,1,
                      1,1,1,1,1,2,2,2,2,1,
                      1,1,1,1,1,2,2,2,2,1,
                      1,1,1,1,1,2,2,2,2,1,
                      1,1,1,1,1,1,1,1,1,1),nrow = 8,byrow = T))
matplot(fpr,tpr,xlim=c(-0,1),ylim=c(-0,1),type="l",xlab="FPR",ylab="TPR",col=1,
        main="ROC curves for the built models")
matplot(fpr,tpr,type="l",add=T,col=1)
abline(0,1,lty=5)
matplot(auc_mat,type="l",ylim=c(0,1),bty="n",main="AUC of the ROC curves",
        ylab="AUC",xlab="Component",col=1)
abline(h=(0:10)/10,col="gray90");abline(v=1:7,col="gray90")
matplot(auc_mat,pch=16,add=T,type = "p",col=1)
matplot(auc_mat,add=T,type = "l",col=1)

