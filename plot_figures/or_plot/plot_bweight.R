# import
options (warn = -1)
library("devtools")
if(exists("cfcausal:::summary_CI")){
  rm(list = c("summary_CI"))
}
devtools::load_all(".")
source("plot_figures/or_plot/get_ORplots.R")
library(readxl)

######vk######
# 导入数据
data_name <- "bweight"
this_dat <- read.csv(paste0("./data/",data_name,".csv"))
# 筛选处理组和控制组
A_all <- as.numeric(this_dat$MomSmoke == "1")
# 定义协变量矩阵
X <- this_dat[, c("Black","Married","Boy", "MomAge","MomWtGain", "Visit", "MomEdLevel")]
X$Black <- as.factor(X$Black)
X$Married <- as.factor(X$Married)
X$Boy <- as.factor(X$Boy)
X$Visit <- as.factor(X$Visit)
X$MomEdLevel <- as.factor(X$MomEdLevel)
X1 <- model.matrix(~ . - 1, X)
# 定义响应变量
Y_all <- this_dat$Weight


main <- data_name
#path <- paste0("figures/OR/",title)
OR_bw <- get_plots(X1, A_all, Y_all, main)

text_data <- c(summary(OR_bw),OR_bw)
write(text_data, file="figures/OR/bweight.txt")
