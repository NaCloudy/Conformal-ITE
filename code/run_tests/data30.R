#############################
######### setting ###########
#############################
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


######### get paras ###########
suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser()
parser$add_argument("--gmm_star", type = "double", default = 1, help = "SA parameter, >=1")
parser$add_argument("--alpha", type="double", default=0.2, help="miscoverage")
parser$add_argument("--save", type="logical", default=TRUE, help="save")
parser$add_argument("--seed1", type = "double", default = 123, help = "data random seed")
parser$add_argument("--seed2", type = "double", default = 456, help = "model random seed")
parser$add_argument("--ntrial", type = "integer", default = 10, help = "number of trials")
parser$add_argument("--method", type = "character", default = 'mean', help = "mean or cqr")
parser$add_argument("--save_par", type = "character", default = './results/ITE/', help = "save parent location")
parser$add_argument("--data_name", type = "character", default = 'data30', help = "data name")

args <- parser$parse_args()
method <- args$method
alpha <- args$alpha
gmm_star <- args$gmm_star
ntrial<- args$ntrial
data.seed <- args$seed1
model.seed <- args$seed2
save_par <- args$save_par
data_name <- args$data_name
save_path <- paste0(save_par,data_name,"-",method,"/")
q<- c(alpha/2, 1- (alpha/2))

###############################
######### read data ###########
###############################
this_dat <- read_xlsx(paste0("data/",data_name,".xlsx"))

######## treatment and control ########
A <- as.numeric(this_dat$Ngroup == "1")

####### define X matrix ###########
X <- this_dat[, c("Male","Age","Height", "Weight","BMI")]
X[X$Male == 2,]$Male <- 0
X$Male <- as.factor(X$Male)
X1 <- model.matrix(~ . - 1, X)

########## define y ##########
Y_all <- this_dat$UricAcid #HU
record <- replicate(2,matrix(0,nrow=ntrial,ncol=3), simplify=FALSE)

########## save path ##########
folder<- paste0(save_path,'alpha_',alpha,'_gmm_',gmm_star, '/')
dir.create(folder, recursive=TRUE, showWarnings = FALSE)


###############################
######## train & test #########
###############################
  n<- length(Y_all)
  trainprop <- 0.8
  set.seed(data.seed)
  trainid <- sample(n, floor(n * trainprop))
  print(paste0("alpha is ",alpha))

  ####### train and test data #######
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

  colnames(X) <- c("Male0", "Male1", "Age","Height", "Weight","BMI")
  colnames(Xtest) <- c("Male0", "Male1", "Age","Height", "Weight","BMI")

for (iter in 1:ntrial){
  ####### prediction band #######

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
    #print(paste0("#effective:",sum(ci_cqr$effective)))

  }

  ####### save #######
  df <- as.data.frame(t(data))
  write.csv(data, file=paste0(folder, 'ntrial_', iter, '.csv'))
}
