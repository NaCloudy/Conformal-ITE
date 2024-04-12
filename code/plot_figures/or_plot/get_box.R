# 定义获取箱线图的函数
get_box <- function(X1, A_all, Y_all, main) {
  library("devtools")
  devtools::load_all(".")
  psfun <- Boosting
  n <- length(Y_all)
  trainprop <- 0.8
  set.seed(123)
  # 随机划分训练集
  trainid <- sample(n, floor(n * trainprop))
  set.seed(NULL)
  # 划分训练集
  X <- X1[trainid,]
  Y <- Y_all[trainid]
  A <- A_all[trainid]
  Xtest <- X1[-trainid,]
  ps <- do.call(psfun, list(Xtest = Xtest, Y = A , X = X))
  record <- c()
  for(j in 1:dim(X)[2]){
    # 对除了当前协变量外其他变量进行boosting预测
    ps_j <- do.call(psfun, list(Xtest = Xtest[,-c(j)], Y = A , X = X[,-c(j)]))
    # 计算当前变量优势比，即当前自变量对响应变量的影响程度
    odds_j <- (ps/(1-ps))/(ps_j/(1-ps_j))
    record <- c(record,odds_j)
  }
  # 找出record向量中大于等于1的优势比的位置
  idx <- which(record>=1)
  R <- record
  # 保留原来大于等于1的优势比
  R[idx] <- record[idx]
  # 将原来小于1的优势比转换为大于1的优势比
  R[-idx] <- 1/record[-idx]
  quantile(R, seq(0, 1, by=.05))
  summary(R)
  # 将直方图保存为png文件
  png(paste("figures/OR/",main,"_box.png", sep = ""))
  # boxplot的标题为main
  boxplot(R, main = main)
  dev.off()
}
