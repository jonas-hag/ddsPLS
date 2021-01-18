simu_bair_2006 <- function(){

  p <- 5000
  n <- 100
  X <- matrix(NA,n,p)
  U <- matrix(runif(n*3,0,1),n,3)
  for(i in 1:p){
    for(j in 1:n){
      if(i<=50 & j<=50){
        X[j,i] <- 3 + rnorm(1)
      }else if(i<=50 & j>50){
        X[j,i] <- 4 + rnorm(1)
      }else if(50<i & i<=100){
        X[j,i] <- 3.5+1.5*(U[j,1]<0.4)*1 +rnorm(1)
      }else if(100<i & i<=200){
        X[j,i] <- 3.5+0.5*(U[j,2]<0.7)*1 +rnorm(1)
      }else if(200<i & i<=300){
        X[j,i] <- 3.5-1.5*(U[j,3]<0.3)*1 +rnorm(1)
      }else if(300<i){
        X[j,i] <- 3.5+rnorm(1)
      }
    }
  }

  y <- matrix(NA,n,1)
  for(j in 1:n){
    y[j,1] <- 1/(25)*sum(X[j,1:50]) + rnorm(1,sd=1.5)
  }

  lambdas  <- seq(0,1,length.out = 30)

  res_bair <- sparse_PLS_Bootstrap(list(x=X),y,type="CT",
                                   paras=lambdas,
                                   n_B=100,
                                   lowExplainedVariance=0,
                                   deflatX=T,NCORES=6,center=T,verbose=T)
}
