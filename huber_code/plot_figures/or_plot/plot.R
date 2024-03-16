library("readxl")
# 导入自定义函数
source("get_hist.R")
source("get_box.R")
######drug_AS######
# 导入数据
# 读取data.xlsx文件的第一个工作表
islander<-read.csv('/Users/niubei/Desktop/保形预测/datasets/drugged_AS.csv')
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
  # 将drug_AS和i拼接成字符串，作为文件名
  file_name <- paste("drug_AS", i, sep = "_")
  get_hist(X1, A, Y_all, file_name)
  get_box(X1, A, Y_all, file_name)

}
######drug_TS######
# 导入数据
# 读取data.xlsx文件的第一个工作表
islander<-read.csv('/Users/niubei/Desktop/保形预测/datasets/drugged_TS.csv')
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
  # 调用自定义函数
  # 将drug_TS和i拼接成字符串，作为文件名
  file_name <- paste("drug_TS", i, sep = "_")
  get_hist(X1, A, Y_all, file_name)
  get_box(X1, A, Y_all, file_name)
}
######vd######
# 读取data.xlsx文件的第一个工作表
vd <- read_excel("/Users/niubei/Desktop/保形预测/datasets/Data_set.xlsx", sheet = 1)
# 数据预处理
A <- as.numeric(vd$Group == "A")
# 定义协变量矩阵
X <- vd[, c("Sex", "Age", "Height", "BW","FIB4","APRI","VD0","AST0","ALT0","Plt0","TGF0","TIMP0","MMP0","P3NP0")]
X1 <- model.matrix(~ . - 1, X)
# 定义响应变量
Y_all <- log2(vd$TGF6)
# 调用boosting算法进行预测
colnames(X1) <- c("Sex","Age31.40","Age41.50","Age51.60","Age61.70","Age71.80","Height","BW","FIB4","APRI","VD0","AST0","ALT0","Plt0","TGF0","TIMP0","MMP0","P3NP0")
file_name <- "vd"
get_hist(X1, A, Y_all, file_name)
get_box(X1, A, Y_all, file_name)

######vk######
# 导入数据
# 读取data.xlsx文件的第一个工作表
vk <-read.csv('/Users/niubei/Desktop/保形预测/datasets/VK2.csv')
vk$Gender <- ifelse(vk$Gender  == "Male", 1, 0)
vk$Access <- ifelse(vk$Access  == "AVFistula", 1, 0)
# 使用ifelse函数，把"Yes"替换为1，把"No"替换为0
# 将需要变成0-1变量的列进行转换
# 生成模型矩阵并赋值给mm
mm <- model.matrix(~ HTN + DM + HCV + Smoking + Heartfailure + ISHD, data = vk)
# 把mm数据框中的虚拟编码替换到vk数据框中
# 创建一个向量col_names，存储mm数据框中除了截距项之外的列名
col_names <- colnames(mm)[-1]

# 使用for循环遍历col_names向量中的每个元素
for (col in col_names) {
  # 使用赋值符号<-把mm数据框中对应列的值覆盖到vk数据框中对应列上
  vk[, col] <- mm[, col]
}
# 筛选处理组和控制组
A <- as.numeric(vk$T == 1)
# 定义协变量矩阵
X <- vk[, c("Gender","HTNYes","DMYes" , "HCVYes","SmokingYes", "HeartfailureYes" ,"ISHDYes","Access","Age","Durationofdialysis", "PTH", "Ca.Pre","PHPre", "CaxPProductPre")]
X1 <- model.matrix(~ . - 1, X)
# 定义响应变量
Y_all <- vk$MGPPre
file_name <- "vk"
get_hist(X1, A, Y_all, file_name)
get_box(X1, A, Y_all, file_name)
