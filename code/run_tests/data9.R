#########设置###########
options (warn = -1)
library("devtools")
if(exists("cfcausal:::summary_CI")){
  rm(list = c("summary_CI"))
}
devtools::load_all(".")
library("cfcausal")
library("dplyr")
library("ggplot2")
library("bannerCommenter")
library("readxl")
options(scipen=999)
#### Get parameters
suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser()
parser$add_argument("--gmm_star", type = "double", default = 1, help = "SA parameter, >=1")
parser$add_argument("--alpha", type="double", default=0.2, help="miscoverage")
parser$add_argument("--save", type="logical", default=TRUE, help="save")
parser$add_argument("--seed", type = "double", default = 1, help = "random seed")
parser$add_argument("--ntrial", type = "integer", default = 5, help = "number of trials")
parser$add_argument("--method", type = "character", default = 'mean', help = "mean or cqr")
parser$add_argument("--save_par", type = "character", default = './results/ITE/', help = "save parent location")
parser$add_argument("--data_name", type = "character", default = 'data9', help = "data name")
args <- parser$parse_args()
method <- args$method
alpha <- args$alpha
gmm_star <- args$gmm_star
ntrial<- args$ntrial
seed <- args$seed
save_par <- args$save_par
data_name <- args$data_name
save_path <- paste0(save_par,data_name,"-",method,"/")
q<- c(alpha/2, 1- (alpha/2))

# 导入数据
# 读取this_data数据集
this_dat <- read.csv(paste0("data/",data_name,".csv"))
# preprocess
this_dat[this_dat$gender == 2,]$gender <- 0
this_dat$bun <- as.numeric(this_dat$bun)
this_dat$rr <- as.numeric(this_dat$rr)
this_dat$hr <- as.numeric(this_dat$hr)
this_dat$hospital <- as.factor(this_dat$hospital)
# missing values
row.has.na <- apply(this_dat, 1, function(x){any(is.na(x))})
this_data <- this_dat[!row.has.na,]
rownames(this_data) <- 1:nrow(this_data)
# 筛选处理组和控制组
A <- as.numeric(this_data$steroid == "1")
# 定义协变量矩阵
X <- this_data[, c("age", "gender","hospital", "adl","wheeze", "bun",
    "rr", "ams", "hr", "hot")]
X1 <- model.matrix(~ . - 1, X)
X1 <- subset(X1, select = -hospital5)
# 定义响应变量
Y_all <- this_data$time_to_stability
record <- replicate(2,matrix(0,nrow=ntrial,ncol=3), simplify=FALSE)
# 创建文件路径
folder<- paste0(save_path,'alpha_',alpha,'_gmm_',gmm_star, '/')
dir.create(folder, recursive=TRUE, showWarnings = FALSE)
##########训练测试##########
for (iter in 1:ntrial){
  n<- length(Y_all)
  trainprop <- 0.8
  set.seed(123)
  trainid <- sample(n, floor(n * trainprop))
  set.seed(NULL)
  print(paste0("alpha is ",alpha))
  # 训练集数据的划分
  Y_obs <- Y_all[trainid]
  X <- X1[trainid,]
  T_obs <- A[trainid]
  Y1 <- Y_obs
  Y1[which(T_obs==0)] <- NA
  Y0 <- Y_obs
  Y0[which(T_obs==1)] <- NA
  id <- seq(1, n)
  testid<- id[!(id %in% trainid)]
  Xtest <- X1[testid,]

  # 生成预测区间
  colnames(X) <- c("age", "gender",
                   "hospital1", "hospital2", "hospital3", "hospital4", #"hospital5",
                   "adl","wheeze", "bun", "rr", "ams", "hr", "hot")
  colnames(Xtest) <- c("age", "gender",
                       "hospital1", "hospital2", "hospital3", "hospital4", #"hospital5",
                       "adl","wheeze", "bun", "rr", "ams", "hr", "hot")
  if(method == "mean"){
    obj_mean <- nested_conformalSA(X, Y1, Y0, T_obs, gmm_star, type = "mean",quantiles=list(), outfun='RF',psparams = list(bag.fraction = 0.8,n.minobsinnode = 5) )
    obj_bands_mean <- predict.nested(obj_mean, X, Y_obs, T_obs, alpha = alpha)
    ci_mean <- fit_and_predict_band(obj_bands_mean, Xtest,'RF')
    ci_list <- list(ci_mean)
    print_list <- list("ci_mean")
    data <- cbind(ci_mean)
    colnames(data) <- c("mean_low", "mean_high", "mean_y1_mean","mean_y0_mean")
  }else if (method == "cqr"){
    obj_cqr <- nested_conformalSA(X, Y1, Y0, T_obs, gmm_star, type = "CQR",quantiles=q, outfun='quantRF',psparams = list(bag.fraction = 0.8,n.minobsinnode = 5))
    obj_bands_cqr <- predict.nested(obj_cqr, X, Y_obs, T_obs, alpha = alpha)
    ci_cqr <- fit_and_predict_band(obj_bands_cqr, Xtest,'RF')
    ci_list <- list(ci_cqr)
    print_list <- list("ci_cqr")
    data <- cbind(ci_cqr)
    colnames(data) <- c("cqr_low", "cqr_high", "cqr_y1_mean","cqr_y0_mean")
  }
  df <- as.data.frame(t(data))
  write.csv(data, file=paste0(folder, 'ntrial_', iter, '.csv'))
}
