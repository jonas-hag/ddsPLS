## ----include=FALSE-------------------------------------------------------
library(htmltools)
tagList(rmarkdown::html_dependency_font_awesome())

## ---- eval=F-------------------------------------------------------------
#  # For no named matrices
#  Xs <- list(X1,X2,X3)
#  
#  # For names matrices
#  Xs <- list(X_dna=X1,X_proteo=X2,X_RNA=X3)

## ---- fig.show='hold',message=FALSE--------------------------------------
library(ddsPLS)

## ---- fig.show='hold',message=FALSE,eval=T-------------------------------
data("liver.toxicity")
X <- scale(liver.toxicity$gene)
Y <- scale(liver.toxicity$clinic)
mddsPLS_model_reg <- mddsPLS(Xs = X,Y = Y,lambda=0.85,R = 3,
                             mode = "reg",verbose = TRUE)

## ----fig.height=7,fig.width=10-------------------------------------------
plot(mddsPLS_model_reg,vizu = "weights",variance = "Linear",
     super = T,comp = c(1,2),addY = T,mar_left = 3)

## ----fig.height=7,fig.width=10-------------------------------------------
plot(mddsPLS_model_reg,vizu = "heatmap",comp = 1)

## ----fig.height=7,fig.width=10-------------------------------------------
plot(mddsPLS_model_reg,vizu = "heatmap",comp = 2)

## ---- fig.show='hold',message=FALSE,eval=T-------------------------------
summary(mddsPLS_model_reg,plot_present_indiv = F)

## ----message=FALSE,eval=F------------------------------------------------
#  n_lambda <- 50
#  NCORES <- 7
#  res_cv_reg <- perf_mddsPLS(Xs = X,Y = Y,
#                             R = 1,lambda_min=0.5,n_lambda=n_lambda,
#                             mode = "reg",NCORES = NCORES,kfolds = "loo")
#  
#  res_cv_reg_L0 <- perf_mddsPLS(Xs = X,Y = Y,
#                             R = 1,L0=1:20,
#                             mode = "reg",NCORES = NCORES,kfolds = "loo")

## ----echo=F--------------------------------------------------------------
# save(res_cv_reg,file="res_cv_reg.RData")
# save(res_cv_reg_L0,file="res_cv_reg_L0.RData")
load("res_cv_reg_noXs.RData")
load("res_cv_reg_L0_noXs.RData")
res_cv_reg$Xs <- list(X)
res_cv_reg_L0$Xs <- list(X)

## ---- fig.show='hold',message=FALSE,eval=T,fig.width=12, fig.height=12----
plot(res_cv_reg,which_sd_plot = c(5,7),ylim=c(-0.3,1.1),alpha.f = 0.4,plot_mean = T,legend_names = colnames(Y))

## ---- fig.show='hold',message=FALSE,eval=T,fig.width=12, fig.height=7----
plot(res_cv_reg_L0,which_sd_plot = c(5,7),ylim=c(-0.2,1.1),alpha.f = 0.4,
     plot_mean = T,legend_names = colnames(Y),pos_legend = "top",no_occurence = T)

## ---- fig.show='hold',message=FALSE,eval=T-------------------------------
summary(res_cv_reg,plot_res_cv =F)

## ---- fig.show='hold',fig.width=7, fig.height=5,message=FALSE,eval=T-----
data("penicilliumYES")
X <- penicilliumYES$X
X <- scale(X[,which(apply(X,2,sd)>0)])
classes <- c("Melanoconidium","Polonicum","Venetum")
Y <- as.factor(unlist(lapply(classes,
                             function(tt){rep(tt,12)})))
mddsPLS_model_class <- mddsPLS(Xs = X,Y = Y,lambda = 0.956,R = 2,
                               mode = "clas",verbose = TRUE)

## ----fig.height=7,fig.width=10-------------------------------------------
plot(mddsPLS_model_class,vizu = "weights",super=T,comp = c(1,2),addY = T,mar_left = 2)

## ----fig.height=7,fig.width=10-------------------------------------------
plot(mddsPLS_model_class,vizu = "heatmap",comp = 1)

## ----fig.height=7,fig.width=10-------------------------------------------
plot(mddsPLS_model_class,vizu = "heatmap",comp = 2)

## ---- fig.show='hold',message=FALSE,eval=T-------------------------------
summary(mddsPLS_model_class,plot_present_indiv = F)

## ---- fig.show='hold',fig.width=10, fig.height=5,message=FALSE,eval=T----
plot(mddsPLS_model_class$mod$T_super,col=Y,pch=as.numeric(Y)+15,cex=2,
     xlab="1st X component, 2 var. selected",
     ylab="2nd X component, 2 var. selected")
legend(-2,0,legend=classes,col=1:3,pch=15+(1:3),box.lty=0,y.intersp=2)

## ----fig.width=7, fig.height=6,message=FALSE,eval=F----------------------
#  n_lambda <- 50
#  NCORES <- 7
#  res_cv_class <- perf_mddsPLS(X,Y,R = 2,lambda_min=0.95,n_lambda=n_lambda,
#                               mode = "clas",NCORES = NCORES,
#                               fold_fixed = rep(1:12,3))

## ----echo=F--------------------------------------------------------------
# save(res_cv_class,file="res_cv_class.RData")
load(file="res_cv_class_noXs.RData")
res_cv_class$Xs <- list(X)

## ----fig.width=10, fig.height=6,message=FALSE,eval=T---------------------
plot(res_cv_class,legend_names = levels(Y),pos_legend="bottomleft")

## ----eval=F--------------------------------------------------------------
#  NCORES <- 7
#  res_cv_class_L0 <- perf_mddsPLS(X,Y,R = 2,L0s = 1:4,
#                               mode = "clas",NCORES = NCORES,
#                               fold_fixed = rep(1:12,3))

## ----eval=T,echo=F-------------------------------------------------------
# save(res_cv_class_L0,file="res_cv_class_L0.RData")
load(file="res_cv_class_L0_noXs.RData")
res_cv_class_L0$Xs <- list(X)

## ----eval=T,fig.height=6,fig.width=10------------------------------------
plot(res_cv_class_L0,legend_names = levels(Y),pos_legend="bottomright")

## ----eval=T--------------------------------------------------------------
model_lambda <- mddsPLS(X,Y,R = 2,lambda = 0.957,mode = "cla" )
plot(model_lambda,super = T)

## ----eval=T--------------------------------------------------------------
model_L0 <- mddsPLS(X,Y,R = 2,L0 = 3,mode = "cla")
plot(model_L0,super = T,addY = T,pos_legend = "right")

## ------------------------------------------------------------------------
data("liver.toxicity")
X <- scale(liver.toxicity$gene)
Y <- scale(liver.toxicity$clinic)
p1=p2 <- 1000
p3 <- ncol(X)-p1-p2
Xs <- list(Block1=X[,1:p1],Matrix2=X[,p1+1:p2],Other_Variables=X[,p1+p2+1:p3])
Xs$Block1[1:10,] <- NA
Xs$Matrix2[5+1:10,] <- NA
Xs$Other_Variables[10+1:20,] <- NA

model_multi_vizu <- mddsPLS(Xs,Y,lambda = 0.8,R = 3)

## ----fig.height=10,fig.width=10------------------------------------------
plot(model_multi_vizu,vizu = "weights",super = T,comp=1 ,addY = T,mar_left = 5)

## ----fig.height=10,fig.width=10------------------------------------------
plot(model_multi_vizu,vizu = "weights",super = T,comp=2 ,addY = T,mar_left = 5)

## ----fig.height=10,fig.width=10------------------------------------------
plot(model_multi_vizu,vizu = "weights",super = T,comp=3 ,addY = T,mar_left = 5)

## ------------------------------------------------------------------------
plot(model_multi_vizu,vizu = "heatmap",comp = 1)

## ------------------------------------------------------------------------
plot(model_multi_vizu,vizu = "heatmap",comp = 2)

## ---- fig.show='hold',message=FALSE,eval=T-------------------------------
summary(model_multi_vizu,plot_present_indiv = T)

