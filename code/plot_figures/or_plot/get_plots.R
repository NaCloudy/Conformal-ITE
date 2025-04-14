# 定义获取箱线图的函数
get_plots <- function(X1, A_all, Y_all, main) {
  #library("devtools")
  #devtools::load_all(".")
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
  ps <- do.call(psfun, list(Xtest = Xtest, Y = A , X = X, n.minobsinnode = 1, bag.fraction = 0.7))
  record <- c()
  for(j in 1:dim(X)[2]){
    # 对除了当前协变量外其他变量进行boosting预测
    ps_j <- do.call(psfun, list(Xtest = Xtest[,-c(j)], Y = A , X = X[,-c(j)], n.minobsinnode = 1, bag.fraction = 0.7))
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
  #quantile(R, seq(0, 1, by=.05))

  #save OR summary stats
  #Summ <- summary(R)
  #print(Summ)
  #write.csv(as.data.frame(Summ), file=paste0("figures/OR/",main,'.csv'))

  # save boxplot
  png(paste("figures/OR/",main,"_box.png", sep = ""))
  boxplot(R, main = main)
  dev.off()

  # save histogram
  png(paste("figures/OR/",main,"_hist.png", sep = ""))
  hist(R, breaks = seq(1, 5, 0.5), prob=TRUE, xlab='(Adjusted) odds ratio',main=main,
       cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  dens <- density(R, from=1, width=1.5)
  lines(dens)
  dev.off()

  return(R)
}
