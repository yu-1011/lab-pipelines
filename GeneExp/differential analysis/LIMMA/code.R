####
# Describe: use LIMMA packages to make differential analysis
# Environment: R 3.3.2
# Author: Yu Chen (SKLMG)
# Input: gene expression matrix
# Output: differential analysis matrix
# R Package requirement: LIMMA
# Data: 2017-04-15
####

#设置工作目录,确保该目录下有存放芯片文本文件以及注释文件
setwd(“G:/R”)
#下载安装脚本
source("http://bioconductor.org/biocLite.R") 
#安装所有 Bioconductor 核心包
biocLite() 
#在bioconductor网站下载limma包
biocLite("limma")
#加载limma包
library(limma)
#读取芯片文本文件
date<-read.table("GSE56015_series_matrix.txt", header = TRUE, sep = "",comment.char = "!",row.names="ID_REF")
#显示基因表达矩阵前六行
 head(date)
#构建实验设计矩阵
condition=factor(c(1,2,3,1,2,3,1,2,3,4,5,6,4,5,6,4,5,6))
design <- model.matrix(~-1+condition)
# 构 建 对比 模 型，比较两个实验条件下表达数据，其中condition2代表只用DMSO处理过的HSPC细胞，condition5代表用石蒜碱处理过的HSPC
contrast.matrix <- makeContrasts (contrasts = "condition2 - condition5", levels = design)
#线性模型拟合
fit <- lmFit(date, design)
#根据比对模型进行差值计算#
fit1 <- contrasts.fit(fit, contrast.matrix)
# 贝 叶 斯 检验
fit2 <- eBayes(fit1) 
# 生成所有 基因的检验结果报表
dif <- topTable(fit2, coef = "condition2 - condition5",  n = nrow(fit2), lfc = log2(1.5))
# 根据P.Value对结果进行筛选，得到全部差异表达基因
dif <- dif[dif[, "P.Value"] < 0.01, ] 
#查看结果的前6行
head(dif)
#将结果保存在dif文件中
write.table(dif,"dif.txt",row.names=T, sep="\t",quote=F)
