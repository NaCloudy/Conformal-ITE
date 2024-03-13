#########设置###########
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
parser$add_argument("--path", type = "character", default = './out/', help = "save location")
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
## 寻找收缩因子
get_fct<- function(set, ite){
  fcts <- seq(1, 0.5, by=-0.01)
  for(fct in fcts){
    ci <-shrink(set, fc=fct)
    coverage <-  mean((ite >= ci[, 1]) & (ite <= ci[, 2]),na.rm = TRUE)
    if(coverage < 0.8){
      record <- fct + 0.01
      break
    }}
  return(record)}

Xfun <- function(n, d){
  matrix(runif(n * d), nrow = n, ncol = d)
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
print_list <- list("sa_mean", "sa_cqr")
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

  ### 训练：两种方法得到预测结果

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

  ci_mean <-fit_and_predict_band(obj_bands_mean,Xtest, 'quantRF')
  ci_cqr <- fit_and_predict_band(obj_bands_cqr,Xtest, 'quantRF')

  record_mean <- 1
  record_cqr <- 1
  record_mean <- get_fct(ci_mean, ite)
  print("now it's the cqr")
  record_cqr <- get_fct(ci_cqr, ite)
  ci_mean[,3] <- record_mean
  ci_cqr[,3] <- record_cqr


  ## 最后得到两组结果：CSA-M, CSA-Q
  ci_list <- list(ci_mean, ci_cqr)

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
    print(paste0(print_list[i], " coverage, ",coverage, ', lens ', len, 'fct',mean(ci[,3] )))
    #第i组，第trial次实验的：覆盖率、区间长度、差值（区间长度）是否有限、最大fct
    record[[i]][trial,] <- c(coverage,len,mean(ci[,3] ))
  }
  print("#################")
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

###收缩因子fct
fct_val <-c()
for (i in 1:length(print_list)){fct_val[[i]]<- as.vector(record[[i]][,3])}
data <- data.frame(fct_val= unlist(fct_val),
                   group=rep(c("CSA-M","CSA-Q"),each=ntrial))

#存储
if(save){
  write.csv(data, paste0(folder,'fct','.csv'), row.names = FALSE)
}
