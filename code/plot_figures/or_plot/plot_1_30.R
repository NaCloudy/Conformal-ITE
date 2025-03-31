library("readxl")
# 导入自定义函数
source("plot_figures/or_plot/get_plots.R")
source("plot_figures/or_plot/get_hist.R")
######drug_AS######
# 导入数据
######1######
# 读取data.xlsx文件的第一个工作表
data1 <- read_xlsx("data/data1.xlsx")
# 数据预处理
A <- as.numeric(data1$group == "1")
# 定义协变量矩阵
X <- data1[, c("Sex", "AgeGroup","Education", "Smoker", "Income", "Living","TNf0", "GPx10", "IL80", "SOD30")]
X[X$Sex == 2,]$Sex <- 0
X$Sex <- as.factor(X$Sex)
X$AgeGroup <- as.factor(X$AgeGroup)
X$Education <- as.factor(X$Education)
X$Smoker <- as.factor(X$Smoker)
X$Income <- as.factor(X$Income)
X$Living <- as.factor(X$Living)
X1 <- model.matrix(~ . - 1, X)
# 定义响应变量
Y_all <- data1$SOD32#TNf2
# # 调用boosting算法进行预测
#colnames(X1) <- c("Sex", "AgeGroup","Education", "Smoker", "Income", "Living","TNf0", "GPx10", "IL80", "SOD30")
file_name <- "Data1-0322"
R <- get_plots(X1, A, Y_all, file_name)
R <- get_hist(X1, A, Y_all, file_name)

######30######
# 导入数据
# 读取data.xlsx文件的第一个工作表
data30 <- read_xlsx("data/data30.xlsx")
# 筛选处理组和控制组
A <- as.numeric(data30$Ngroup == "1")
# 定义协变量矩阵
X <- data30[, c("Male", "Age","Height", "Weight","BMI")]
X[X$Male == 2,]$Male <- 0
X$Male <- as.factor(X$Male)
X1 <- model.matrix(~ . - 1, X)
# 定义响应变量
Y_all <- data30$UricAcid
# 调用boosting算法进行预测
file_name <- "Data30"
R <- get_plots(X1, A, Y_all, file_name)
R <- get_hist(X1, A, Y_all, file_name)




######9######
# 读取data.xlsx文件的第一个工作表
this_dat <- read.csv("data/data9.csv")
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
# # 调用boosting算法进行预测
#colnames(X1) <- c("Sex", "AgeGroup","Education", "Smoker", "Income", "Living","TNf0", "GPx10", "IL80", "SOD30")
file_name <- "Data9"
R <- get_hist(X1, A, Y_all, file_name)
R <- get_box(X1, A, Y_all, file_name)
