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
library("h2o")
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
  obj_mean <- nested_conformalSA(X, Y1, Y0, T_obs, gmm_star, type = "mean",outfun='huberBoosting')
  object <- obj_mean
  Y <- Y_obs
  T <- T_obs
  alpha = 0.2
  wthigh = 20
  wtlow = 0.05
  PSmodel <- object$PSmodel
  wtfun <- object$wtfun
  group2id <- object$group2id
  ite <- object$ite
  Xtest <- X[group2id, , drop=FALSE]
  Ytest <- Y[group2id]
  Ttest <- T[group2id]
  #ite_test <- ite[group2id]
  
  inds1 <- which(Ttest == 1)
  inds0 <- which(Ttest == 0)
  
  Xtest_treated <- Xtest[inds1, , drop=FALSE]
  Ytest_treated <- Ytest[inds1]
  Xtest_control <- Xtest[inds0, , drop=FALSE]
  Ytest_control <- Ytest[inds0]
  ###############
  object <- object$treated
  type <- object$type
  Yhat_test <- object$Ymodel(Xtest)
  wt_test <- object$wtfun(Xtest,object$gmm)
  # mean prediction
  Yslack <- cutoff_SA(object$Yscore,object$wt,wt_test,alpha)
  Ylo <- Yhat_test - Yslack
  Yup <- Yhat_test + Yslack
  Ylo <- as.vector(Ylo)
  Yup <- as.vector(Yup)
  interval <- data.frame(lower = Ylo, upper = Yup)
  
  
  
  
  
  obj_bands_mean <- predict.nested(obj_mean, X, Y_obs, T_obs, alpha = alpha)
  ci_mean <- fit_and_predict_band(obj_bands_mean, Xtest,'RF')
  ci_cqr <- fit_and_predict_band(obj_bands_cqr, Xtest,'RF')
  
  ci_list <- list(ci_mean,ci_cqr)
  print_list <- list("ci_mean", "ci_cqr")
  data <- cbind(ci_mean, ci_cqr)
  colnames(data) <- c("mean_low", "mean_high", "mean_y1_mean","mean_y0_mean","cqr_low", "cqr_high"
                      ,"cqr_y1_mean","cqr_y0_mean")
  df <- as.data.frame(t(data))
  write.csv(data, file=paste0(folder, 'ntrial_', iter, '.csv'))
  
  for(i in 1:length(ci_list)){
    ci <- ci_list[[i]]
    # 负效应和正效应的比例
    negative <- sum(ci[, 2] <=0)
    positive <- sum(ci[, 1] >=0)
    neg_percentage <- negative/dim(ci)[1]
    pos_percentage <- positive/dim(ci)[1]
    # 差异
    diff <- ci[, 2] - ci[, 1]
    # 区间长度
    len <-mean(diff[is.finite(diff)])
    print(print_list[i])
    print(paste0( ' neg_percentage: ', neg_percentage,
                  ' pos_percentage: ', pos_percentage, "  len: ", len))
    print("###############")
    record[[i]][iter,] <- c(neg_percentage,pos_percentage,len)
  }






