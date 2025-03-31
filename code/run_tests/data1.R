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
parser$add_argument("--seed1", type = "double", default = 123, help = "data random seed")
parser$add_argument("--seed2", type = "double", default = 456, help = "model random seed")
parser$add_argument("--ntrial", type = "integer", default = 3, help = "number of trials")
parser$add_argument("--method", type = "character", default = 'mean', help = "mean or cqr")
parser$add_argument("--save_par1", type = "character", default = './results/ITE/', help = "save ITE parent location")
parser$add_argument("--save_par2", type = "character", default = './results/class/', help = "save classification parent location")
parser$add_argument("--data_name", type = "character", default = 'data1', help = "data name")
args <- parser$parse_args()
method <- args$method
alpha <- args$alpha
gmm_star <- args$gmm_star
ntrial<- args$ntrial
data.seed <- args$seed1
model.seed <- args$seed2
data_name <- args$data_name
save_par1 <- args$save_par1
save_path1 <- paste0(save_par1,data_name,"-",method,"/")
save_par2 <- args$save_par2
save_path2 <- paste0(save_par2,data_name,"-",method,"/")
q<- c(alpha/2, 1- (alpha/2))

# 导入数据
# 读取this_data数据集
this_data <- read_xlsx(paste0("data/",data_name,".xlsx"))
# 筛选处理组和控制组
A <- as.numeric(this_data$group == "1")
# 定义协变量矩阵
X <- this_data[, c("Sex", "AgeGroup","Education", "Smoker", "Income", "Living",
    "TNf0", "GPx10", "IL80", "SOD30")]
X[X$Sex == 2,]$Sex <- 0
X$AgeGroup = as.factor(X$AgeGroup)
X$Education = as.factor(X$Education)
X$Smoker = as.factor(X$Smoker)
X$Income = as.factor(X$Income)
X$Living = as.factor(X$Living)
X1 <- model.matrix(~ . - 1, X)
# 定义响应变量
Y_all <- this_data$SOD32 #GPx12#IL82#SOD32#TNf2
record <- replicate(2,matrix(0,nrow=ntrial,ncol=3), simplify=FALSE)
# 创建文件路径
folder1 <- paste0(save_path1,'alpha_',alpha,'_gmm_',gmm_star, '/')
folder2 <- paste0(save_path2,'alpha_',alpha,'_gmm_',gmm_star, '_ntrial_', ntrial, '/')
dir.create(folder1, recursive=TRUE, showWarnings = FALSE)
dir.create(folder2, recursive=TRUE, showWarnings = FALSE)
##########训练测试##########

##对同一批人
n <- length(Y_all)
trainprop <- 0.8
set.seed(data.seed)
trainid <- sample(n, floor(n * trainprop))
#set.seed(NULL)
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
train_dat <- NA

for (iter in 1:ntrial){
  # 生成预测区间
  colnames(X) <- c("Sex", "AgeGroup1", "AgeGroup2", "AgeGroup3", "AgeGroup4", "Education2", "Smoker2", "Income2", "Income3", "Living2",
                   "TNf0", "GPx10", "IL80", "SOD30")
  colnames(Xtest) <- c("Sex", "AgeGroup1", "AgeGroup2", "AgeGroup3", "AgeGroup4", "Education2", "Smoker2", "Income2", "Income3", "Living2",
                       "TNf0", "GPx10", "IL80", "SOD30")
  if(method == "mean"){
    obj_mean <- nested_conformalSA(X, Y1, Y0, T_obs, gmm_star, type = "mean",quantiles=list(), outfun='RF',psparams = list(bag.fraction = 0.8,n.minobsinnode = 5),data.seed = data.seed, model.seed = model.seed)
    obj_bands_mean <- predict.nested(obj_mean, X, Y_obs, T_obs, alpha = alpha, data.seed = data.seed, model.seed = model.seed)
    ci_mean <- fit_and_predict_band(obj_bands_mean, Xtest, testid, 'RF',data.seed = data.seed, model.seed = model.seed)
    ci_list <- list(ci_mean)

    class_dat <- cbind(Xtest, effective=ci_mean$effective)
    rownames(class_dat) <- NULL
    train_dat <- rbind(train_dat, class_dat)

    data <- cbind(ci_mean)
    colnames(data) <- c("mean_low", "mean_high", "mean_y1_mean","mean_y0_mean", "id", "effectiveness")
  }else if (method == "cqr"){
    obj_cqr <- nested_conformalSA(X, Y1, Y0, T_obs, gmm_star, type = "CQR",quantiles=q, outfun='quantRF',psparams = list(bag.fraction = 0.8,n.minobsinnode = 5),data.seed = data.seed, model.seed = model.seed)
    obj_bands_cqr <- predict.nested(obj_cqr, X, Y_obs, T_obs, alpha = alpha,data.seed = data.seed, model.seed = model.seed)
    ci_cqr <- fit_and_predict_band(obj_bands_cqr, Xtest, testid, 'RF',data.seed = data.seed, model.seed = model.seed)
    ci_list <- list(ci_cqr)

    class_dat <- cbind(Xtest, effective=ci_cqr$effective)
    rownames(class_dat) <- NULL
    train_dat <- rbind(train_dat, class_dat)

    data <- cbind(ci_cqr)
    colnames(data) <- c("cqr_low", "cqr_high", "cqr_y1_mean","cqr_y0_mean", "id", "effectiveness")
  }

  df <- as.data.frame(t(data))
  write.csv(data, file=paste0(folder1, 'ntrial_', iter, '.csv'))
}


####test####
# train_dat <- train_dat[c(-1),]
# class_model <- randomForest::randomForest(x=train_dat[,c(1:14)], y= train_dat[,c(14)])
# imp <- randomForest::importance(class_model)
# imp <- cbind(imp, imp[,1]/sum(imp[,1])*100)
# imp <- as.data.frame(imp[order(-imp[, 1]), ])
# colnames(imp) <- c("IncNodePurity", "Percentage")
# plot <- randomForest::varImpPlot(class_model)
#
# write.csv(as.data.frame(imp), file=paste0(folder2, 'rf.csv'))


