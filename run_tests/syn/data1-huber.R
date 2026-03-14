#########设置###########
library("devtools")
if(exists("cfcausal:::summary_CI")){
  rm(list = c("summary_CI"))
}
devtools::load_all(".")
library("cfcausal")
library("dplyr")
library("ggplot2")
library("bannerCommenter")
library(MASS)
library("h2o")
h2o.init()
options(scipen=999)
#### Get parameters
suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser()
parser$add_argument("--gmm_star", type = "double", default = 1, help = "SA parameter, >=1")
parser$add_argument("--alpha", type="double", default=0.2, help="miscoverage")
parser$add_argument("--cftype", type="integer", default=2, help="confounding type")
parser$add_argument("--fct", type="double", default=1, help="shrinkage, <=1")
parser$add_argument("--save", type="logical", default=TRUE, help="save")
parser$add_argument("--seed1", type = "double", default = 123, help = "data random seed")
parser$add_argument("--seed2", type = "double", default = 456, help = "model random seed")
parser$add_argument("--ntrial", type = "integer", default = 5, help = "number of trials,50")
parser$add_argument("--path", type = "character", default = './results/synthetic/data1_huber_para/', help = "save location")
parser$add_argument("--ntrain", type = "integer", default = 5, help = "training numbers,3000")
parser$add_argument("--ntest", type = "integer", default = 5, help = "testing numbers,10000")
parser$add_argument("--errdist", type = "character", default = 'norm', help = "error distribution,norm,heavy,norm_p")
parser$add_argument("--huber_alpha", type = "integer", default = 0.9, help = "huber alpha, [0,1]")
args <- parser$parse_args()
alpha <- args$alpha
gmm_star <- args$gmm_star
cftype<- args$cftype
fct <- args$fct
ntrial<- args$ntrial
seed1 <- args$seed1
seed2 <- args$seed2
save <- args$save
errdist <- args$errdist
n1 <- args$ntrain   # 训练集个数
ntest <- args$ntest # 测试集个数
path <- args$path
huber_alpha <- args$huber_alpha
q<- c(alpha/2, 1- (alpha/2))

# 导入数据
# 读取VD数据集
vd <- read.csv("data/VD.csv")
# 筛选处理组和控制组
A <- as.numeric(vd$Group == "A")
# 定义协变量矩阵
Xdata <- vd[, c("Sex", "Age", "Height", "BW","FIB4","APRI","VD0","AST0","ALT0","Plt0","TGF0","TIMP0","MMP0","P3NP0")]

# 后面是离散数据
R <- cov(Xdata[,-c(1:2)]) # 连续变量的样本协方差
m <- colMeans(Xdata[,-c(1:2)]) # 连续变量样本均值

Xfun<-function(n1,X,R,m){#需要读入原始数据矩阵X和需要生成的样本个数n1
  mydata<-mvrnorm(n1,m,R)#多元正态分布生成模拟连续数据
  #第一列 性别
  n<-sum(X[,1]-1)#第一列中变量为2的个数
  p<-n/dim(X)[1]
  x1<-rmultinom(n=n1, size=1, prob=c(1-p,p))
  v1<-rep(NA,n1)
  v1[which(x1[1,]==1)]<-1
  v1[which(x1[2,]==1)]<-2
  #第二列 年龄
  p4<-length(which(X[,2]=='61-70'))/dim(X)[1]
  p3<-length(which(X[,2]=='51-60'))/dim(X)[1]
  p2<-length(which(X[,2]=='41-50'))/dim(X)[1]
  p1<-length(which(X[,2]=='31-40'))/dim(X)[1]
  p5<-length(which(X[,2]=='71-80'))/dim(X)[1]
  x2<-rmultinom(n=n1, size=1, prob=c(p1,p2,p3,p4,p5))
  v2<-rep(NA,n1)
  v2[which(x2[1,]==1)]<-1
  v2[which(x2[2,]==1)]<-2
  v2[which(x2[3,]==1)]<-3
  v2[which(x2[4,]==1)]<-4
  v2[which(x2[5,]==1)]<-5#1\2\3\4\5分别代表31-40、41-50.。。几个年龄群体
  #离散和连续数据合并形成最终的模拟数据,离散数据做one hot 处理
  V<-as.matrix(mydata)
  df<-cbind(V,v1,v2)
  #df<-as.data.frame(df)
  df[,14]<-factor(df[,14])
  colnames(df)[13]<-'Sex'
  colnames(df)[14]<-'Age'
  df<-as.data.frame(df)
  df$Age<-as.factor(df$Age)
  X1<-model.matrix(~ .-1,df )
  return(X1)
}

sdfun <- function(X){
    rep(1, nrow(X))
}

taufun <- function(X){
  2 / (1 + exp(-5 * (0.01*X[, 1] - 0.5))) * 2 / (1 + exp(-5 * (0.01*X[, 2] - 0.5)))
}

pscorefun <- function(X){
  x1<-X[, 1]
  x1<-(x1-min(x1))/(max(x1)-min(x1))
  (1 + pbeta(1-x1, 2, 4)) / 4
}

get_Y1obs <- function(X){
  if(errdist=='norm'){
    return(taufun(X) + sdfun(X) * rnorm(dim(X)[1]))
  }else if(errdist=='heavy'){
    return(taufun(X) + sdfun(X) * rlogis(dim(X)[1], -0.1075211, 1.2992436))
  }else if(errdist=='norm_p'){
    return(taufun(X) + sdfun(X) * rnorm(dim(X)[1], 0, 2.5))
  }
}

taufun0 <- function(X){
  taufun(X) + 10*sin(X[, 3])*(1/(1+exp(-0.05*X[, 3])))
}

get_Y0obs <- function(X){
  if(errdist=='norm'){
    return(taufun0(X) + sdfun(X) * rnorm(dim(X)[1]))
  }else if(errdist=='heavy'){
    return(taufun0(X) + sdfun(X) * rlogis(dim(X)[1], 0.07311189, 1.08353882))
  }
  else if(errdist=='norm_p'){
    y <- rnorm(dim(X)[1],0,1)
    sam <- sample(1:dim(X)[1], 0.07*dim(X)[1])
    y[sam] = y[sam] + 6
    return(taufun(X) + sdfun(X) * y)
  }
}

shrink <- function(set,fc){
  newset <- set
  idx <- is.finite(set[,1])
  center <- (set[idx,2] + set[idx,1])/2
  halflen <- (set[idx,2] - set[idx,1])/2
  newset[idx,] <- cbind(center-halflen*fc, center+halflen*fc)
  return(newset)
}


print_list <- list("sa_huber_0.1", "sa_huber_0.3", "sa_huber_0.5", "sa_huber_0.7", "sa_huber_0.9")
record <- replicate(length(print_list),matrix(0,nrow=ntrial,ncol=3), simplify=FALSE)

for (trial in 1:ntrial){
  ##---------------------------------------------------------------
  ##                        Generate observed data                         -
  ##---------------------------------------------------------------
  X<-Xfun(n1,Xdata,R,m)
  n<-nrow(X)
  ps <- pscorefun(X)
  T <- as.numeric(runif(n)<ps)

  Y0 <- get_Y0obs(X)
  Y1 <- get_Y1obs(X)

  Y_obs <- Y1*T + Y0*(1-T)

  Y1[which(T==0)] <- NA
  Y0[which(T==1)] <- NA

  ##----------------------------------------------------------------
  ##                            Train on Group1
  ##----------------------------------------------------------------
  obj_huber_1 <- nested_conformalSA(X, Y1, Y0, T_obs, gmm_star, type = "mean", quantiles=list(), outfun='huberBoosting', outparams=list(huber_alpha = 0.1))
  obj_bands_huber_1 <- predict.nested(obj_huber_1, X, Y_obs, T, alpha = alpha)

  obj_huber_3 <- nested_conformalSA(X, Y1, Y0, T_obs, gmm_star, type = "mean", quantiles=list(), outfun='huberBoosting', outparams=list(huber_alpha = 0.3))
  obj_bands_huber_3 <- predict.nested(obj_huber_3, X, Y_obs, T, alpha = alpha)

  obj_huber_5 <- nested_conformalSA(X, Y1, Y0, T_obs, gmm_star, type = "mean", quantiles=list(), outfun='huberBoosting', outparams=list(huber_alpha = 0.5))
  obj_bands_huber_5 <- predict.nested(obj_huber_5, X, Y_obs, T, alpha = alpha)

  obj_huber_7 <- nested_conformalSA(X, Y1, Y0, T_obs, gmm_star, type = "mean", quantiles=list(), outfun='huberBoosting', outparams=list(huber_alpha = 0.7))
  obj_bands_huber_7 <- predict.nested(obj_huber_7, X, Y_obs, T, alpha = alpha)

  obj_huber_9 <- nested_conformalSA(X, Y1, Y0, T_obs, gmm_star, type = "mean", quantiles=list(), outfun='huberBoosting', outparams=list(huber_alpha = 0.9))
  obj_bands_huber_9 <- predict.nested(obj_huber_9, X, Y_obs, T, alpha = alpha)

  ##Testing
  Xtest <- Xfun(ntest,Xdata,R,m)
  pstest <- pscorefun(Xtest)
  Ttest <- as.numeric(runif(ntest)<pstest)
  id1 <- which(Ttest==1)
  id0 <- which(Ttest==0)

  Y0test <- rep(NA,ntest)
  Y0test[id0] <- get_Y0obs(Xtest[id0,])
  Y0test[id1] <- samplecf(Xtest[id1,], taufun0, sdfun, case=cftype, gmm=gmm_star)
  #这里用的是case2，也就是appendix C里的方法产生反事实预测Y0
  Y1test <- rep(NA,ntest)
  Y1test[id1] <- get_Y1obs(Xtest[id1,])
  Y1test[id0] <- samplecf(Xtest[id0,],taufun, sdfun, case=cftype, gmm=gmm_star)

  ##----------------------------------------------------------------
  ##                       ITE & evaluation                    -
  ##----------------------------------------------------------------

  ite <- Y1test - Y0test

  #CSA-huber
  ci_huber_copy_1 <- fit_and_predict_band(obj_bands_huber_1, Xtest,'quantRF')
  ci_huber_1 <- shrink(ci_huber_copy_1, fc=fct)
  ci_huber_1[, 3] <- ci_huber_1[,3]
  ci_huber_1[, 4] <- ci_huber_1[, 4]

  ci_huber_copy_3 <- fit_and_predict_band(obj_bands_huber_3, Xtest,'quantRF')
  ci_huber_3 <- shrink(ci_huber_copy_3, fc=fct)
  ci_huber_3[, 3] <- ci_huber_3[,3]
  ci_huber_3[, 4] <- ci_huber_3[, 4]

  ci_huber_copy_5 <- fit_and_predict_band(obj_bands_huber_5, Xtest,'quantRF')
  ci_huber_5 <- shrink(ci_huber_copy_5, fc=fct)
  ci_huber_5[, 3] <- ci_huber_5[,3]
  ci_huber_5[, 4] <- ci_huber_5[, 4]

  ci_huber_copy_7 <- fit_and_predict_band(obj_bands_huber_7, Xtest,'quantRF')
  ci_huber_7 <- shrink(ci_huber_copy_7, fc=fct)
  ci_huber_7[, 3] <- ci_huber_7[,3]
  ci_huber_7[, 4] <- ci_huber_7[, 4]

  ci_huber_copy_9 <- fit_and_predict_band(obj_bands_huber_9, Xtest,'quantRF')
  ci_huber_9 <- shrink(ci_huber_copy_9, fc=fct)
  ci_huber_9[, 3] <- ci_huber_9[,3]
  ci_huber_9[, 4] <- ci_huber_9[, 4]

  ## 最后得到四组结果：CSA-huber, CSA-M, CSA-Q, ITE-NUC, BART
  ci_list <- list(ci_huber_1, ci_huber_3, ci_huber_5, ci_huber_7, ci_huber_9)

  for(i in 1:length(ci_list)){
    #保形区间
    ci <- ci_list[[i]]
    #区间长度
    diff <- ci[, 2] - ci[, 1]
    # 找出符合条件的索引
    index <- which(diff > 9999999) # 人为设定
    # 将符合条件的值修改为Inf
    ci[index, 2] <- Inf # 人为修改
    ci[index, 1] <- Inf # 人为修改
    diff[index] <- Inf  # 人为修改

    #覆盖率
    coverage <- mean((ite >= ci[, 1]) & (ite <= ci[, 2]),na.rm = TRUE)
    #平均区间长度（有限值）
    len <- mean(diff[is.finite(diff)])
    #无限长度的区间个数
    n_inf <- sum(is.infinite(diff))

    #输出
    print(paste0(print_list[i], " coverage, ",coverage, ', lens ', len))
    #第i组，第trial次实验的：覆盖率、区间长度、差值（区间长度）是否有限
    record[[i]][trial,] <- c(coverage,len,n_inf)
  }
  print(paste0("################# trial ",trial," #################"))
}


##----------------------------------------------------------------
##                         Save results                         --
##----------------------------------------------------------------

#create a new path for files
if(fct <1){
  folder <- paste0(path,"/", errdist, "/", "fct_", fct,"/") #忽略，我们取fct==1
}else{
  folder<- paste0(path, "/", errdist, "/", "gmm_",gmm_star,"/")
}
print(folder)
dir.create(folder, recursive=TRUE, showWarnings = FALSE)

#coverage data
coverage <-c()
for (i in 1:length(print_list)){coverage[[i]]<- as.vector(record[[i]][,1])}

data <- data.frame(Coverage=unlist(coverage),
                   group=rep(c("sa_huber_0.1", "sa_huber_0.3", "sa_huber_0.5", "sa_huber_0.7", "sa_huber_0.9"),
                             each=ntrial))

if(save){
  write.csv(data, paste0(folder,'coverage', '.csv'), row.names = FALSE)
}

#length data
Interval_length <-c()
for (i in 1:length(print_list)){Interval_length[[i]]<- as.vector(record[[i]][,2])}
data <- data.frame(Interval_length= unlist(Interval_length),
                   group=rep(c("sa_huber_0.1", "sa_huber_0.3", "sa_huber_0.5", "sa_huber_0.7", "sa_huber_0.9"),each=ntrial))

if(save){
  write.csv(data, paste0(folder,'len','.csv'), row.names = FALSE)
}
