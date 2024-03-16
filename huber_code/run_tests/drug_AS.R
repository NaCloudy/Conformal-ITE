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
parser$add_argument("--path", type = "character", default = './results/drug_AS/ITE/', help = "save location")
args <- parser$parse_args()
alpha <- args$alpha
gmm_star <- args$gmm_star
ntrial<- args$ntrial
seed <- args$seed
save <- args$save
path = args$path
q<- c(alpha/2, 1- (alpha/2))

# 导入数据
# 读取drugged_AS.csv
islander <- read.csv("data/drugged_AS.csv")
islander$Happy_Sad_group <- ifelse(islander$Happy_Sad_group == "H", 1, 0)
library(dplyr)
data_1 <- filter(islander, T == 1 | T == 0)
data_2 <- filter(islander, T == 2 | T == 0)
data_3 <- filter(islander, T == 3 | T == 0)
# 定义一个向量，包含三个数据集的名称
datasets <- c("data_1", "data_2", "data_3")
# 对每个数据集进行循环
for (i in 1:length(datasets)){
  # 使用get函数，根据名称获取数据集
  data <- get(datasets[i])
  # 筛选处理组和控制组
  A <- as.numeric(data$T == i)
  # 定义协变量矩阵
  X <- data[, c("age", "Mem_Score_Before", "Happy_Sad_group")]
  X1 <- model.matrix(~ . - 1, X)
  # 定义响应变量
  Y_all <- data$Diff
  record <- replicate(2,matrix(0,nrow=ntrial,ncol=3), simplify=FALSE)
  # 创建文件路径，根据数据集的名称修改文件名
  folder<- paste0(path,'alpha_',alpha,'_gmm_',gmm_star, '/')
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
    obj_mean <- nested_conformalSA(X, Y1, Y0, T_obs, gmm_star, type = "mean",quantiles=list(), outfun='RF',psparams = list(bag.fraction = 0.8,n.minobsinnode = 5) )
    # obj_cqr <- nested_conformalSA(X, Y1, Y0, T_obs, gmm_star, type = "CQR",quantiles=q, outfun='quantRF',psparams = list(bag.fraction = 0.8,n.minobsinnode = 5))
    obj_bands_mean <- predict.nested(obj_mean, X, Y_obs, T_obs, alpha = alpha)
    # obj_bands_cqr <- predict.nested(obj_cqr, X, Y_obs, T_obs, alpha = alpha)
    ci_mean <- fit_and_predict_band(obj_bands_mean, Xtest,'RF')
    # ci_cqr <- fit_and_predict_band(obj_bands_cqr, Xtest,'RF')
    ci_list <- list(ci_mean)
    print_list <- list("ci_mean")
    data <- cbind(ci_mean)
    df <- as.data.frame(t(data))
    # 根据数据集的名称修改文件名
    write.csv(data, file=paste0(folder, datasets[i], '_ntrial_', iter, '.csv'))
  }
}






