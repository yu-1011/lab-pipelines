####
# Describe: use affy packages to preprocess .CEL
# Environment: R 3.3.2
# Author: Yu Chen (SKLMG)
# Input:.CEL file
# Output:Expression matrix
# R Package requirement: affy, tcltk
# Data: 2017-08-15
####

library(affy)
library(tcltk)
dir <- tk_choose.dir(caption = "Select folder")
filters <- matrix(c("CEL file", ".[Cc][Ee][Ll]", "All", ".*"), ncol = 2, byrow = T)
cel.files <- tk_choose.files(caption = "Select CELs", multi = TRUE,
                             filters = filters, index = 1)
basename(cel.files)
data.raw <- ReadAffy(filenames = cel.files)
sampleNames(data.raw) <- paste("CHIP", 1:length(cel.files), sep = "-")
sampleNames(data.raw)
pData(data.raw)

#####背景处理
################RMA
data.rma <- bg.correct(data.raw, method="rma")
class(data.rma)
head(pm(data.raw)-pm(data.rma),2)
head(mm(data.raw)-mm(data.rma),2)#差值为全部为0，说明rma方法没有对mm数值进行处理
identical(mm(data.raw), mm(data.rma))
################MAS5
data.mas <- bg.correct(data.raw, method="mas")
class(data.mas)
head(pm(data.mas),2)
head(pm(data.raw)-pm(data.mas),2)
identical(mm(data.mas), mm(data.mas))

#####归一化处理
################分位数（quantile）方法
data.mas.qt <- normalize(data.mas, method = "quantiles")
head(pm(data.mas.qt)/pm(data.mas), 5)

data.rma.qt <- normalize(data.rma, method = "quantiles")
head(pm(data.rma.qt)/pm(data.rma), 5)


#####汇总
eset.rma.liwong <- computeExprSet(data.rma.qt, pmcorrect.method="pmonly",
                                  summary.method="liwong")
eset.mas.liwong <- computeExprSet(data.mas.qt, pmcorrect.method="pmonly",
                                  summary.method="liwong")

###############一步完成预处理
eset.rma <- rma(data.raw)
eset.mas5 <- mas5(data.raw)

###############计算基因表达量
emat.rma.log2 <- exprs(eset.rma.liwong)
results.mas<- exprs(eset.mas.liwong)
class(emat.rma.log2)
emat.rma.nologs <- 2^emat.rma.log2
results.rma <- data.frame(emat.rma.log2)#我们这套数据没有技术重复

emat.rma.log2 <- exprs(eset.rma)
head(emat.rma.log2)

################  p/m/a calls()#若无mm probe 可略过此步
#data.mas5calls <- mas5calls(data.raw)
mas5calls.AffyBatch(data.raw, ids = NULL, verbose = TRUE, tau = 0.015,
                    alpha1 = 0.04, alpha2 = 0.06,
                    ignore.saturated=TRUE) #调整参数选择p<0.05的结果
eset.mas5calls <- exprs(data.mas5calls )
head(eset.mas5calls)

##################如果有技术重复, 应取平均值
#计算平均值，并做对数转换
#results.rma <- data.frame((emat.rma.log2[,c(1,3,5,7)] + emat.rma.log2[,c(2,4,6,8)])/2)
##################计算表达量差异倍数
#results.rma$fc.1h <- results.rma[,2]-results.rma[,1] #因为是对数,直接相减即可
#results.rma$fc.24h <- results.rma[,3]-results.rma[,1]
#results.rma$fc.7d <- results.rma[,4]-results.rma[,1]
#head(results.rma, 2)

############写出
write.exprs(eset.rma , file="eset.rma.txt")
