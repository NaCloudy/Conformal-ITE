########设定########

#### 引入库
library("devtools")
if(exists("cfcausal:::summary_CI")){
  rm(list = c("summary_CI"))
}
devtools::load_all(".")
library("cfcausal")
library("dplyr")
library("ggplot2")
library("bannerCommenter")
options(scipen=999)

#### 参数设定
  suppressPackageStartupMessages(library("argparse"))
  parser <- ArgumentParser()
  parser$add_argument("--n", type = "integer", default = 3000, help = "Sample size")
  parser$add_argument("--d", type = "integer", default = 20, help = "Dimension")
  parser$add_argument("--gmm_star", type = "double", default = 3, help = "SA parameter, >=1")
  parser$add_argument("--alpha", type="double", default=0.2, help="miscoverage")
  parser$add_argument("--cftype", type="integer", default=2, help="confounding type")
  parser$add_argument("--dtype", type="character", default='homo', help="data type, homo or het")
  parser$add_argument("--fct", type="double", default=1, help="shrinkage, <=1")
  parser$add_argument("--save", type="logical", default=TRUE, help="save")
  parser$add_argument("--seed", type = "double", default = 1, help = "random seed")
  parser$add_argument("--ntrial", type = "integer", default =5, help = "number of trials")
  parser$add_argument("--path", type = "character", default = './loop_fct/', help = "save location")
  args <- parser$parse_args()
  n <- args$n
  d <- args$d
  alpha <- args$alpha
  ntrial<- args$ntrial
  dtype <- args$dtype
  cftype<- args$cftype
  fct <- args$fct
  seed <- args$seed
  gmm_star <- args$gmm_star
  save <- args$save
  path = args$path
  #分位数
  q <- c(alpha/2, 1-alpha/2)

#######产生数据的定义#######
##fct=用于数据集缩小的shrinkage factor；如果<0，作变换
if(fct <=0){
  fct <- exp(fct)
}

##根据同方差/异方差，定义Xfun和sdfun
#同方差，X iid~unif(0,1)，sd = 1
#异方差，X = pnorm(X'), X' = X*根号(1-rho)+fac*根号(rho)；sd~unif(0.5,1.5)
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

##定义倾向得分e(X)
pscorefun <- function(X){
  (1 + pbeta(1-X[, 1], 2, 4)) / 4
}

##定义偏置的分布：正态分布
errdist <- rnorm

###公式26，定义Y_i(1)
##定义taufun为公式26的f(X_{i1})*f(X_{i2})
taufun <- function(X){
  2 / (1 + exp(-5 * (X[, 1] - 0.5))) * 2 / (1 + exp(-5 * (X[, 2] - 0.5)))
}
##定义Y_i(1)，公式26
get_Y1obs <- function(X){
  #sdfun(X) * errdist(dim(X)[1]) ~ N(0,sd^2)
  return(taufun(X) + sdfun(X) * errdist(dim(X)[1]))
}

###公式27，定义Y_i(0)
##定义taufun0为公式27的f(X_{i1})*f(X_{i2})+...
taufun0 <- function(X){
  taufun(X) + 10*sin(X[, 3])*(1/(1+exp(-5*X[, 3])))
}
##定义结果Y_i(0)，公式27
get_Y0obs <- function(X){
  return(taufun0(X) + sdfun(X) * errdist(dim(X)[1]))
}

##根据sfct来缩小现有数据集的范围来，创建新的数据集：原始数据集和缩放因子
##用以判断某种CSA方法的sharpness（fc之后，还能不能达到覆盖率1-alpha）
#新的半区间长度 = 原来半区间长度*fc，因此新区间长度 = 原来长度*fc
shrink <- function(set,fc){
  newset <- set
  #排除输入数据集中缺失或非数字值的行，选择剩下的行index
  idx <- is.finite(set[,1])
  #计算缩小后区间的中心点（求第一和第二列的平均值）
  center <- (set[idx,2] + set[idx,1])/2
  #计算缩小后区间的半长度（求第一和第二列差值的一半）
  halflen <- (set[idx,2] - set[idx,1])/2
  #
  newset[idx,] <- cbind(center-halflen*fc, center+halflen*fc)
  return(newset)
}

##定义结果的数据结构（全空）
print_list <- list("sa_mean", "sa_cqr", "ite_nuc", "sa_naive")
record <- replicate(length(print_list),matrix(0,nrow=ntrial,ncol=3), simplify=FALSE)

##### 进行ntrial次实验 #####
for (trial in 1:ntrial){
  ### 产生观测数据observed
  #X
  X <- Xfun(n,d)
  #e(X)
  ps <- pscorefun(X)
  #treatment
  T <- as.numeric(runif(n)<ps)
  #由方法定义出的Y0、Y1
  Y0 <- get_Y0obs(X)
  Y1 <- get_Y1obs(X)
  #真实观测到的数据：T和Y下标相符的那些数据
  Y_obs <- Y1*T + Y0*(1-T)
  #不符的那些数据：mask掉
  Y1[which(T==0)] <- NA
  Y0[which(T==1)] <- NA

  ### 训练：四种方法得到预测结果
  #BART方法-Bonferroni correction
  obj1_ite <- conformal_SA(X, Y1, gmm_star, type = "mean", outfun='RF')
  obj0_ite <- conformal_SA(X, Y0, gmm_star, type = "mean", outfun='RF')

  #ITE-NUC
  #conformalCF和conformalITe和sampleCF
  CIfun_inexact <- conformalIte(X, Y_obs, T, alpha = alpha,
                                algo = "nest", exact=FALSE, type = "CQR",
                                #lofun = 'RF', upfun = 'RF', citype = "mean",
                                quantiles = c(alpha/2, 1- (alpha/2)), outfun = "quantRF",  useCV = FALSE)

  ##CSA-M和CSA-Q
  #在第一组上进行训练
  obj_mean <- nested_conformalSA(X, Y1, Y0, T, gmm_star, type = "mean",quantiles=list(), outfun='RF')
  obj_cqr <- nested_conformalSA(X, Y1, Y0, T, gmm_star, type = "CQR",quantiles=q, outfun='quantRF')
  #在第二组上获得预测区间
  obj_bands_mean <- predict.nested(obj_mean, X, Y_obs, T, alpha = alpha)
  obj_bands_cqr <- predict.nested(obj_cqr, X, Y_obs, T, alpha = alpha)

  ### 测试
  ##生成测试数据
  ntest <- 10000
  Xtest <- Xfun(ntest,d)
  pstest <- pscorefun(Xtest)
  Ttest <- as.numeric(runif(ntest)<pstest)
  id1 <- which(Ttest==1)
  id0 <- which(Ttest==0)
  #Y0
  Y0test <- rep(NA,ntest)
  Y0test[id0] <- get_Y0obs(Xtest[id0,])
  Y0test[id1] <- samplecf(Xtest[id1,], taufun0, sdfun, case=cftype, gmm=gmm_star)
  #Y1
  Y1test <- rep(NA,ntest)
  Y1test[id1] <- get_Y1obs(Xtest[id1,])
  Y1test[id0] <- samplecf(Xtest[id0,],taufun, sdfun, case=cftype, gmm=gmm_star)

  ## csa-m&csa-q 进行ite estimation
  ##得到结果如：lower、upper、y1_mean、y0_mean；共ntest行
  #ITE真值
  ite <- Y1test - Y0test
  #CSA-M
  ci_mean_copy <-fit_and_predict_band(obj_bands_mean,Xtest, 'quantRF')
  ci_mean <- shrink(ci_mean_copy, fc=fct)
  ci_mean[, 3] <- ci_mean_copy[,3]
  ci_mean[, 4] <- ci_mean_copy[, 4]
  #CSA-Q
  ci_cqr_copy <- fit_and_predict_band(obj_bands_cqr,Xtest, 'quantRF')
  ci_cqr <- shrink(ci_cqr_copy, fc=fct)
  ci_cqr[, 3] <- ci_cqr_copy[,3]
  ci_cqr[, 4] <- ci_cqr_copy[, 4]

  ## bonferroni
  ci0_ite <- predict.conformalmsm(obj0_ite, Xtest,alpha = alpha/2)
  ci1_ite <- predict.conformalmsm(obj1_ite, Xtest,alpha = alpha/2)
  ci_ite <- cbind(ci1_ite[,1] - ci0_ite[,2], ci1_ite[,2] - ci0_ite[,1])

  ## ite-nuc
  ci_inexact <- CIfun_inexact(Xtest)

  ## 最后得到四组结果：CSA-M, CSA-Q, ITE-NUC, BART
  ci_list <- list(ci_mean, ci_cqr, ci_inexact, ci_ite)
  for(i in 1:length(ci_list)){
    #保形区间
    ci <- ci_list[[i]]
    #覆盖率
    coverage <- mean((ite >= ci[, 1]) & (ite <= ci[, 2]),na.rm = TRUE)
    #区间长度
    diff <- ci[, 2] - ci[, 1]
    #平均区间长度（有限值）
    len <- mean(diff[is.finite(diff)])
    #无限长度的区间个数
    n_inf <- sum(is.infinite(diff))

    #输出
    print(paste0(print_list[i], " coverage, ",coverage, ', lens ', len))
    #第i组，第trial次实验的：覆盖率、区间长度、差值（区间长度）是否有限
    record[[i]][trial,] <- c(coverage,len,n_inf)
  }
  print("#################")
  #示例输出：
  #[1] "t = 1 interval length 5.04373544853956"
  #[1] "t = 0 interval length 4.76176775211672"
  #[1] "t = 1 interval length 6.15681244007641"
  #[1] "t = 0 interval length 5.84476704240692"
  #[1] "sa_mean coverage, 0.8728, lens 5.47052226522426"
  #[1] "sa_cqr coverage, 0.9471, lens 6.67399457432251"
  #[1] "ite_nuc coverage, 0.7293, lens 4.22504101028801"
  #[1] "sa_naive coverage, 0.9927, lens 9.03805548461452"
  #[1] "#################"
}

######存储结果#######
###创建文件路径
folder<- paste0(path, dtype)
if(fct <1){
  folder <- paste0(folder,"/", "fct_", fct,"/")
}else{
  folder<- paste0(folder, "/","gmm_tmp",gmm_star,"/")
}
print(folder)
dir.create(folder, recursive=TRUE, showWarnings = FALSE)

###覆盖率coverage
coverage <-c()
for (i in 1:length(print_list)){coverage[[i]]<- as.vector(record[[i]][,1])}

data <- data.frame(Coverage=unlist(coverage),
                   group=rep(c("CSA-M","CSA-Q","ITE-NUS", "CSA-B"),
                             each=ntrial))
#存储
if(save){
  write.csv(data, paste0(folder,'coverage', '.csv'), row.names = FALSE)
}

###区间长度length
Interval_length <-c()
for (i in 1:length(print_list)){Interval_length[[i]]<- as.vector(record[[i]][,2])}
data <- data.frame(Interval_length= unlist(Interval_length),
                   group=rep(c("CSA-M","CSA-Q","ITE-NUS", "CSA-B"),each=ntrial))
#存储
if(save){
  write.csv(data, paste0(folder,'len','.csv'), row.names = FALSE)
}

