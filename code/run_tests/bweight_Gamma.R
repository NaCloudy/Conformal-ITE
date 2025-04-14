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
parser$add_argument("--gmm_step", type = "double", default = 1, help = "Gmm step, start from 1")
parser$add_argument("--alpha", type="double", default=0.2, help="miscoverage")
parser$add_argument("--save", type="logical", default=TRUE, help="save")
parser$add_argument("--seed1", type = "double", default = 123, help = "data random seed")
parser$add_argument("--seed2", type = "double", default = 456, help = "model random seed")
parser$add_argument("--ntrial", type = "integer", default = 2, help = "number of trials")
parser$add_argument("--method", type = "character", default = 'mean', help = "mean or cqr")
parser$add_argument("--save_par", type = "character", default = './results/ITE/', help = "save parent location")
parser$add_argument("--data_name", type = "character", default = 'sampled_bweight_500', help = "data name")

args <- parser$parse_args()
method <- args$method
alpha <- args$alpha
gmm_step <- args$gmm_step
ntrial<- args$ntrial
data.seed <- args$seed1
model.seed <- args$seed2
save_par <- args$save_par
data_name <- args$data_name
save_path <- paste0(save_par,data_name,"-log","-",method,"/")
q<- c(alpha/2, 1- (alpha/2))
gmm_list <- seq(1,5,gmm_step)

###############################
######### read data ###########
###############################
this_dat <- read.csv(paste0("./data/",data_name,".csv"))

######## treatment and control ########
A <- as.numeric(this_dat$MomSmoke == "1")

####### define X matrix ###########
X <- this_dat[, c("Black","Married","Boy", "MomAge","MomWtGain", "Visit", "MomEdLevel")]
X$Black <- as.factor(X$Black)
X$Married <- as.factor(X$Married)
X$Boy <- as.factor(X$Boy)
X$Visit <- as.factor(X$Visit)
X$MomEdLevel <- as.factor(X$MomEdLevel)
X1 <- model.matrix(~ . - 1, X)

########## define y ##########
Y_all <- this_dat$Weight
record <- replicate(2,matrix(0,nrow=ntrial,ncol=3), simplify=FALSE)

########## save path ##########
folder<- paste0(save_path,'alpha_',alpha,'_gmmStep_',gmm_step, '/')
dir.create(folder, recursive=TRUE, showWarnings = FALSE)


###############################
######## train & test #########
###############################
n<- length(Y_all)
trainprop <- 0.8
set.seed(123)
trainid <- sample(n, floor(n * trainprop))
print(paste0("alpha is ",alpha))
print(paste0("gmm",gmm_list))

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

colnames(X) <- c("Black0", "Black1", "Married1", "Boy1", "MomAge", "MomWtGain", "Visit1", "Visit2", "Visit3", "MomEdLevel1", "MomEdLevel2", "MomEdLevel3")
colnames(Xtest) <- c("Black0", "Black1", "Married1", "Boy1", "MomAge", "MomWtGain", "Visit1", "Visit2", "Visit3", "MomEdLevel1", "MomEdLevel2", "MomEdLevel3")

####### get ITE function #######
get_ITE <- function(X, Y1, Y0, T_obs, gmm_star, data.seed, model.seed){
  if(method == "mean"){
    obj_mean <- nested_conformalSA(X, Y1, Y0, T_obs, gmm_star, type = "mean",quantiles=list(), outfun='RF',psparams = list(bag.fraction = 0.8), data.seed = data.seed, model.seed = model.seed)
    obj_bands_mean <- predict.nested(obj_mean, X, Y_obs, T_obs, alpha = alpha, data.seed = data.seed, model.seed = model.seed)
    ci_mean <- fit_and_predict_band(obj_bands_mean, Xtest, testid, 'RF', data.seed = data.seed, model.seed = model.seed)
    return(ci_mean)
  }else if (method == "cqr"){
    obj_cqr <- nested_conformalSA(X, Y1, Y0, T_obs, gmm_star, type = "CQR",quantiles=q, outfun='quantRF',psparams = list(bag.fraction = 0.8), data.seed = data.seed, model.seed = model.seed)
    obj_bands_cqr <- predict.nested(obj_cqr, X, Y_obs, T_obs, alpha = alpha, data.seed = data.seed, model.seed = model.seed)
    ci_cqr <- fit_and_predict_band(obj_bands_cqr, Xtest, testid, 'RF', data.seed = data.seed, model.seed = model.seed)
    return(ci_cqr)
  }
}

for (iter in 1:ntrial){
  ####### find Gamma-value #######
  flag = 0
  for(gmm_star in gmm_list){
    ci <- get_ITE(X, Y1, Y0, T_obs, gmm_star, data.seed, model.seed)
    colnames(ci) <- c("low", "high", "y1_mean","y0_mean", "id", "effectiveness")
    print(paste0("#effective:",sum(ci_mean$effective)))
    if(flag == 1){
      break
    }
  }

  ####### save #######
  df <- as.data.frame(t(data))
  print(paste0("-1=", sum(df$effectiveness==-1), ", 0=",sum(df$effectiveness==0), ", 1=",sum(df$effectiveness==1)))
  write.csv(data, file=paste0(folder, 'ntrial_', iter, '.csv'))
}


#temp1 <- read.csv("./results/ITE/bweight-mean/alpha_0.2_gmm_1/ntrial_5.csv")
#print(paste0("-1=", sum(temp1$effectiveness==-1), ", 0=",sum(temp1$effectiveness==0), ", 1=",sum(temp1$effectiveness==1)))

