#' The summary method of the \emph{mddsPLS} function.
#'
#' This function is easy to use and gives information about the dataset and the model.
#'
#' @param object The object of class mddsPLS
#' @param plot_present_indiv logical. If \emph{TRUE}, plots a venn diagram of the missing
#'  individuals in each \emph{X} dataset. If the dataset does not appear, it means that it
#'   does not have missing values
#' @param ... Other parameters.
#'
#' @importFrom eulerr euler
#' @importFrom graphics plot
#'
#' @seealso  \code{\link{mddsPLS}}
#'
#' @export
#'
#' @examples
#' library(ddsPLS)
#' data("liver.toxicity")
#' X <- scale(liver.toxicity$gene)
#' Y <- scale(liver.toxicity$clinic)
#' X1<-X[,1:10];X1[1,]<-NA
#' X2<-X[,11:20];X2[2:5,]<-NA
#' X3<-X[,21:30];X3[4:20,]<-NA
#' X4<-X[,31:40]
#' Xs <- list(x1=X1,x2=X2,aaaa=X3,X4)
#' # object <- mddsPLS(Xs = Xs,Y = Y[,1],lambda=0.1,R = 1, mode = "reg",verbose = TRUE)
#' # summary(object,plot_present_indiv = T)
summary.mddsPLS <- function (object,plot_present_indiv=TRUE,
                             ...)
{
  K <- length(object$Xs);    sent_K <- paste("Number of blocks:",K)
  R <- object$mod$R;    sent_R <- paste("Number of dimensions:",R)
  if(!is.null(object$L0)){
    sent_lambda <- paste("Regularization coefficient:",object$L0)
  }else{
    sent_lambda <- paste("Regularization coefficient:",object$lambda)
  }
  n <- nrow(object$Xs[[1]]);    sent_n <- paste("Number of individuals:",n)
  ps <- unlist(lapply(object$Xs,ncol))
  na_x <- unlist(lapply(object$id_na,function(oo){if(length(oo)==0){out <- 0}else{out <- length(oo)};out}))
  prop_na <- signif(na_x/n*100,3)
  names_X_block <- names(object$Xs)
  if(length(names_X_block)==0){
    names_X_block <- 1:K
  }else{
    names_X_block <- unlist(lapply(1:K,function(k){
      out <- names_X_block[k]
      if(out==""){
        out <- paste("Block",k)
      }
      out
    }))
  }
  mat_miss <- matrix(NA,3,1+K)
  colnames(mat_miss) <- c(names_X_block,"Total")
  rownames(mat_miss) <- c("Number of variables","Number of missing samples","Proportion of missing samples (%)")
  mat_miss[1,] <- c(ps,sum(ps))
  mat_miss[2,] <- c(na_x,sum(na_x))
  mat_miss[3,] <- c(prop_na,signif(sum(prop_na)/K,3))
  df_miss <- data.frame(matrix(1,n,K))
  for(k in 1:K){popo <- as.numeric(object$id_na[[k]]);if(length(popo)>0){df_miss[popo,k]<-0};
  df_miss[,k] <- factor(df_miss[,k],levels = c(0,1))}
  names(df_miss) <- paste(names_X_block," (",unlist(lapply(object$id_na,length)),")",sep="")
  q <- ncol(object$Y_0); if(is.null(q)){q <- length(object$Y_0)};    sent_q <- paste("Number of variables in Y part:",q)
  mode <- object$mode;if(mode=="reg"){mode <- "regression"}else{mode <- "classification"}
  sent_mode <- paste("Model built in mode",mode)
  maxit <- object$maxIter_imput;    sent_maxit <- paste("Maximum number of iterations in the imputation process:",maxit)
  has_con <- object$has_converged!=0;if(has_con){sent_con <- ""}else{sent_con <- " not"}
  sent_con <- paste("Algorithm of imputation has",sent_con," converged",sep="")

  df_num_var_sel <- data.frame(matrix(NA,K,R))
  rownames(df_num_var_sel) <- names_X_block
  colnames(df_num_var_sel) <- paste("Comp.",1:R)
  for(r in 1:R){
    for(k in 1:K){
      df_num_var_sel[k,r] <- length(which(abs(object$mod$u[[k]][,r])>1e-9))
    }
  }


  cat("=====================================================");cat("\n")
  cat("              ddsPLS object description    ");cat("\n")
  cat("=====================================================");cat("\n")
  cat("\n")
  cat(sent_K);cat("\n")
  cat(sent_R);cat("\n")
  cat(sent_lambda);cat("\n")
  cat(sent_n);cat("\n")
  cat(sent_q);cat("\n")
  cat(sent_mode);cat("\n")
  cat(sent_maxit);cat("\n")
  cat(sent_con);cat("\n")
  cat("\n")
  cat("\n")
  cat("     Variance explained of Y (%)    ");cat("\n")
  cat("---------------------------------");cat("\n")
  cat("By the predicted values");cat("\n")
  print(signif(object$Variances$VAR_FINAL[1],2)*100)
  cat("\n")
  cat("By the Super Components");cat("\n")
  print(signif(object$Variances$VAR_SUPER_COMPS,2)*100)
  cat("\n")
  cat("By each Component of each Block");cat("\n")
  print(signif(object$Variances$VAR_COMPS,2)*100)
  cat("\n")
  cat("\n")
  cat("    Missing value information    ");cat("\n")
  cat("---------------------------------");cat("\n")
  cat("\n")
  print(data.frame(mat_miss));cat("\n")
  cat("\n")
  cat("         mddsPLS results         ");cat("\n")
  cat("---------------------------------");cat("\n")
  cat("\n")
  N_max <- sum(unlist(lapply(object$mod$Ms,function(m){length(which(colSums(abs(m))!=0))})))
  cat(paste("At most ",N_max," variable(s) can be selected in the X part",sep=""));cat("\n")
  a<-lapply(object$mod$u,function(u){apply(u,2,function(u){length(which(abs(u)>1e-9))})})
  cat("    For each block of X, are selected");cat("\n")
  print(df_num_var_sel)
  cat("    For the Y block, are selected");cat("\n")
  cat(paste("        @ (",paste(apply(object$mod$v,2,function(u){length(which(abs(u)>1e-9))}),
                                collapse = ","),") variable(s)",sep=""));cat("\n")
  cat("\n")
  cat("\n")
  cat("                 Thank's for using me      ");cat("\n")
  cat("-----------------------------------------------------");cat("\n")
  cat("                                      Hadrien Lorenzo");cat("\n")
  cat("                       hadrien.lorenzo.2015@gmail.com");cat("\n")
  cat("=====================================================");cat("\n")
  if(plot_present_indiv){
    if (!requireNamespace("eulerr", quietly = TRUE)) {
      stop("Package \"eulerr\" needed for this function to work. Please install it.",
           call. = FALSE)
    }else if(requireNamespace("eulerr", quietly = TRUE)){
      requireNamespace("eulerr")
      model_euler <- eulerr::euler(df_miss)
      plot(model_euler, counts = T, factor_names=T,quantities=T,
           main="Missing samples in each dataset and intersections")
    }
  }
}
