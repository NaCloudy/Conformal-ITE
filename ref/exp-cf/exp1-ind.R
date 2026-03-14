########包的导入########

rm(list = ls())
library("devtools")
devtools::load_all(".")
library("cfcausal")
library("dplyr")
library("ggplot2")

####参数的设置#######
  suppressPackageStartupMessages(library("argparse"))

  parser <- ArgumentParser()
  parser$add_argument("--n", type = "integer", default = 1000, help = "Sample size")
  parser$add_argument("--d", type = "integer", default = 10, help = "Dimension")
  parser$add_argument("--U", type = "double", default = 5, help = "bound of tilting function, >=1")
  parser$add_argument("--alpha", type="double", default=0.2, help="miscoverage")
  parser$add_argument("--seed", type = "double", default = 1, help = "random seed")
  parser$add_argument("--cftype", type = "integer", default = 2, help = 'confounding type')
  parser$add_argument("--ntrial", type = "integer", default = 3, help = "number of trials")
  parser$add_argument("--path", type = "character", default = 'exp1/out/', help = "verbose or not")

  args <- parser$parse_args()
  n <- args$n
  d <- args$d
  U <- args$U
  alpha <- args$alpha
  ntrial<- args$ntrial
  cftype <- args$cftype
  seed <- args$seed
  q<- c(0.5*alpha, 1- (0.5*alpha))
  path = args$path

########数据的生成##########
Xfun <- function(n, d){
  matrix(runif(n * d), nrow = n, ncol = d)
}
taufun <- function(X){
  2 / (1 + exp(-5 * (X[, 1] - 0.5))) * 2 / (1 + exp(-5 * (X[, 2] - 0.5)))
}
sdfun <- function(X){
  rep(1, nrow(X))
}
##倾向得分
pscorefun <- function(X){
  (1 + pbeta(1-X[, 1], 2, 4)) / 4
}
errdist <- rnorm
get_Yobs <- function(X){
  return(taufun(X) + sdfun(X) * errdist(dim(X)[1]))
}

datalist = list() #储存结果
idx <- 1
for (gmm_star in seq(1,4,0.5)){
  print(gmm_star)

  X <- Xfun(n,d)
  Y <- get_Yobs(X)
  ps <- pscorefun(X)
  T <- as.numeric(runif(n)<ps)
  Y[!T] <- NA # 只有接受了处理的观测值才有结果
  
  # 训练基于均值、分位数和ite的模型
  obj_mean <- conformal_SA(X, Y, gmm_star, type = "mean",quantiles=list(), outfun='RF')
  obj_cqr <- conformal_SA(X, Y, gmm_star, type = "CQR", quantiles=q, outfun='quantRF')
  obj_ite <- conformalCf(X, Y, type = 'mean', outfun ='RF', useCV = FALSE)

  # 测试
  ntest <- 1000
  Xtest <- Xfun(ntest,d)
  pstest <- pscorefun(Xtest)
  Ttest <- as.numeric(runif(ntest)<pstest)

  id1 <- which(Ttest==1)
  id0 <- which(Ttest==0)
  Ytest <- rep(NA,ntest)
  Ytest[id1] <- taufun(Xtest[id1,]) + sdfun(Xtest[id1,]) * errdist(length(id1))
  Ytest_cf <- matrix(0,length(id0),1)
  Ytest_cf[,1] <- samplecf(Xtest[id0,],taufun, sdfun, case=cftype, gmm=gmm_star)
  ITE <- Ytest
  ITE[id0] <- c(Ytest_cf)
  # 分别进行预测
  objs <- list(obj_mean,obj_cqr,obj_ite)
  groups <- c('Mean','CQR','ITE')
  for (rs in 1:length(objs)){
    if (groups[rs]=='ITE'){
      ci <- predict(objs[[rs]], Xtest, alpha = alpha)
    }else{
      ci <- predict.conformalmsm(objs[[rs]], Xtest,alpha = alpha)
    }
    dat <- data.frame(lower = ci$lower, upper=ci$upper, ite = ITE,
                      T=Ttest, cover=as.numeric((ITE>=ci$lower)&(ITE<=ci$upper)),
                      X=Xtest)
    dat$group <- groups[rs]
    dat$gmm <- gmm_star
    datalist[[idx]] <- dat
    idx <- idx+1
  }
}

########数据储存########

# 创建新的文件路径
folder<- paste0(path, 'individual', "/")
print(folder)
dir.create(folder, recursive=TRUE)
# 写入新文件
all_data = do.call(rbind, datalist)
write.csv(all_data, paste0(folder,'alpha',10*alpha,'n',n,'.csv'), row.names = FALSE)



