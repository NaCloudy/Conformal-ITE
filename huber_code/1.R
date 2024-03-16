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
options(scipen=999)
#### Get parameters
suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser()
parser$add_argument("--gmm_star", type = "double", default = 3, help = "SA parameter, >=1")
parser$add_argument("--alpha", type="double", default=0.2, help="miscoverage")

parser$add_argument("--save", type="logical", default=TRUE, help="save")
parser$add_argument("--seed", type = "double", default = 1, help = "random seed")
parser$add_argument("--ntrial", type = "integer", default = 5, help = "number of trials")
parser$add_argument("--path", type = "character", default = './results/fish/ITE/', help = "save location")
args <- parser$parse_args()
alpha <- args$alpha
gmm_star <- args$gmm_star
ntrial<- args$ntrial
seed <- args$seed
save <- args$save
path = args$path
q<- c(alpha/2, 1- (alpha/2))


#######数据处理#########
load('./exp-fish/data/fish.Rda')
# 定义high fish consumption
A <- as.numeric(nhanes.fish$fish.level == "high")
# 定义协变量矩阵
X <- nhanes.fish[, c("gender", "age", "income", "income.missing", "race", "education", "smoking.ever", "smoking.now")]
X$race <- factor(X$race)
X1 <- model.matrix(~ . - 1, X)
# 定义响应变量
Y_all <- log2(nhanes.fish$o.LBXTHG)
record <- replicate(2,matrix(0,nrow=ntrial,ncol=3), simplify=FALSE)
# 创建文件路径
folder<- paste0(path, 'alpha_',alpha,'/','gmm_',gmm_star, '/')
dir.create(folder, recursive=TRUE, showWarnings = FALSE)

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
  Ttest <- A[testid]
  Ytest <- Y_all[testid]
  # X的列数
  ncol(X)
  
  # 将Y列附加到X中
  X <- cbind(X, Y_obs)
  colnames(X)[ncol(X)] <- "Y"
  h2o.init()
  # 这里df是一个R的数据框
  h2o_df <- as.h2o(X)
  # h2o_df$T <- as.factor(h2o_df$T)
  h2o_df_test <- as.h2o(Xtest)
  # h2o_df_test$T <- as.factor(h2o_df_test$T)
  fit <- h2o::h2o.gbm(training_frame = h2o_df,
                      x=1:13,    ## the predictor columns, by column index
                      y=14,
                      distribution = "huber", huber_alpha = 0.1)
  res <- h2o::h2o.predict(fit, newdata = h2o_df_test)
  