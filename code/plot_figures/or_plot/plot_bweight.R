library("readxl")
# 导入自定义函数
source("plot_figures/or_plot/get_hist.R")
#source("plot_figures/or_plot/get_box.R")

######vk######
# 导入数据
data_name <- "bweight"
this_dat <- read.csv(paste0("./data/",data_name,".csv"))
# 筛选处理组和控制组
A <- as.numeric(this_dat$MomSmoke == "1")
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
file_name <- data_name
OR <- get_hist(X1, A, Y_all, file_name)
get_box(X1, A, Y_all, file_name)
