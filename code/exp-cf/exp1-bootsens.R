#######设置工作路径#######
setwd("D:/ECNU/23_3summer/conformal/code/ITE-Supplementary Materials/modified_code/exp-cf")
#########包的导入###########
library("cfcausal")
library("devtools")
if(exists("cfcausal:::summary_CI")){
  rm(list = c("summary_CI"))
}
library("dplyr")
library("ggplot2")
library("bootsens")
library("rJava")
devtools::load_all(".")
library("bartMachine")
rm(list = ls())
#### Get parameters
  suppressPackageStartupMessages(library("argparse"))
  parser <- ArgumentParser()
  parser$add_argument("--n", type = "integer", default = 3000, help = "Sample size")
  parser$add_argument("--d", type = "integer", default = 20, help = "Dimension")
  parser$add_argument("--gmm_star", type = "double", default = 3, help = "SA parameter, >=1")
  parser$add_argument("--alpha", type="double", default=0.2, help="mis-coverage")
  parser$add_argument("--dtype", type = "character", default = 'homo', help = 'data type')
  parser$add_argument("--cftype", type = "integer", default = 2, help = 'confounding type')
  parser$add_argument("--seed", type = "double", default = 1, help = "random seed")
  parser$add_argument("--ntrial", type = "integer", default = 1, help = "number of trials")
  parser$add_argument("--path", type = "character", default = './out/exp1_ppl/')
  args <- parser$parse_args()
  n <- args$n
  d <- args$d
  alpha <- args$alpha
  dtype <- args$dtype
  cftype <- args$cftype
  ntrial<- args$ntrial
  seed <- args$seed
  path <- args$path
  save <- TRUE

# BART. bartMachine package needed
BARTs <- function(Y, X, Xtest,ndpost = 100, ...){
  transform <- FALSE
  if(length(class(X))>1){
    transform <- TRUE
  }else if(class(X) != "data.frame"){
    transform <- TRUE
  }
  if (transform){
    X <- as.data.frame(X)
    Xtest <- as.data.frame(Xtest)
    names(Xtest) <- names(X)
  }
  fit <- bartMachine::bartMachine(X, Y, verbose = FALSE)
  ci <- bartMachine::calc_prediction_intervals(fit, Xtest, pi_conf = 1-alpha)
  return(ci$interval)
}

#######数据生成########

if(dtype=='homo'){
  Xfun <- function(n, d){
    matrix(runif(n * d), nrow = n, ncol = d)
  }
  sdfun <- function(X){
    rep(1, nrow(X))
  }
} else{
  rho <- 0.9
  Xfun <- function(n, d){
    X <- matrix(rnorm(n * d), nrow = n, ncol = d)
    fac <- rnorm(n)
    X <- X * sqrt(1 - rho) + fac * sqrt(rho)
    pnorm(X)
  }
  sdfun <- function(X){
    runif(nrow(X),0.5,1.5)
  }
}
taufun <- function(X){
  2 / (1 + exp(-5 * (X[, 1] - 0.5))) * 2 / (1 + exp(-5 * (X[, 2] - 0.5)))
}
pscorefun <- function(X){
  (1 + pbeta(1-X[, 1], 2, 4)) / 4
}
errdist <- rnorm

get_Yobs <- function(X){
  return(taufun(X) + sdfun(X) * errdist(dim(X)[1]))
}

gmm_list <- seq(1,4,1)
record <- matrix(0,nrow=ntrial*length(gmm_list),ncol=3)
print(path)
dir.create(path, recursive=TRUE, showWarnings = FALSE)

#######训练和测试##########
for(gg in 1:length(gmm_list)){
  gmm_star <- gmm_list[gg]
  for (trial in 1:ntrial){
    ##----------------------------------------------------------------
    ##                           Training                           --
    ##----------------------------------------------------------------
    X <- Xfun(n,d)
    Y <- get_Yobs(X)
    ps <- pscorefun(X)
    T <- as.numeric(runif(n)<ps)
    Y[!T] <- 0
    # 估计因果效应的极大值和极小值
    ci_point <- extrema.os(T, X, Y, gamma = log(gmm_star)) # point estimate
    # 估计因果效应的置信区间
    ci <- bootsens.os(T, X, Y, gamma = log(gmm_star), parallel = FALSE, alpha = alpha)
    ci <- unname(ci)
    ##---------------------------------------------------------------
    ##                           Testing                           --
    ##---------------------------------------------------------------

    # 利用BART模型进行区间估计
    # 贝叶斯模型对CPU消耗极大，注意调整样本数量，否则R会崩溃
    ntest <- 10000
    Xtest <- Xfun(ntest,d)
    ci_bart <- BARTs(Y, X, Xtest)
    pstest <- pscorefun(Xtest)
    Ttest <- as.numeric(runif(ntest)<pstest)

    # test outcome & evaluation
    id1 <- which(Ttest==1)
    id0 <- which(Ttest==0)
    Ytest <- rep(NA,ntest)
    Ytest[id1] <- get_Yobs(Xtest[id1,])
    Ytest[id0] <- samplecf(Xtest[id0,],taufun, sdfun, case=cftype, gmm=gmm_star)
    ite <- Ytest-0
    coverage_boot <- mean((ite >= ci[1]) & (ite <= ci[2]),na.rm = TRUE)
    coverage_bart <- mean((ite >= ci_bart[,1]) & (ite <= ci_bart[,2]),na.rm = TRUE)
    #print(paste0("ate-coverage, ",coverage_boot))
    print(paste0("ate-coverage, ",coverage_boot, ', bart-coverage ', coverage_bart))
    #record[(gg-1)*ntrial+trial,] <- c(coverage_boot, gmm_star)
    record[(gg-1)*ntrial+trial,] <- c(coverage_boot, coverage_bart, gmm_star)
  }
  if(save){
    write.csv(record, paste0(path,'ppl-coverage-', dtype, '.csv'), row.names = FALSE)
  }
}
