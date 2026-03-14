# 数据预处理
Aall <- as.numeric(nhanes.fish$fish.level == "high")
Xall <- nhanes.fish[, c("gender", "age", "income", "income.missing",
"race", "education", "smoking.ever", "smoking.now")]
Xall$race <- factor(Xall$race)
X1 <- model.matrix(~ . - 1, Xall)
Yall <- log2(nhanes.fish$o.LBXTHG)

# 初始设定
psfun <- Boosting
n<- length(Yall)
trainprop <- 0.8
set.seed(123)
trainid <- sample(n, floor(n * trainprop))
set.seed(NULL)

# 划分训练集
X <- X1[trainid,]
Y <- Yall[trainid]
A <- Aall[trainid]
Xtest <- X1[-trainid,]

# do.call用于调用可变长参数函数
# 调用boosting算法进行预测
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
boxplot(R)
# 将图像保存到pdf
pdf("./calib.pdf",width=5.5,height=6)
hist(R, breaks = 10, prob=TRUE, xlab='(Adjusted) odds ratio',main=NULL,
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
dens <- density(R, from=1, width=1.5)
lines(dens)
par(mar = c(4, 4, 0.1, 0.1))
dev.off()
