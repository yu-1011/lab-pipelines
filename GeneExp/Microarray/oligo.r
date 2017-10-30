####
# Describe: use affy packages to preprocess .CEL
# Environment: R 3.3.2
# Author: Yu Chen (SKLMG)
# Input:.CEL file
# Output:Expression matrix
# R Package requirement: oligo, pd.mirna.4.0#这个包是你相应的注释文件
# Data: 2017-08-15
####

library("oligo")
library("pd.mirna.4.0")
data.dir <- "D:/calendar_chen/乐于助人/戴佳成/microRNA/CEL"
(celfiles <- list.files(data.dir, "\\.CEL$"))
data.raw <- read.celfiles(filenames = file.path(data.dir, celfiles))

basename(celfiles)
sampleNames(data.raw) <- paste("CHIP",1:length(celfiles))
sampleNames(data.raw)
pData(data.raw)

par(mfrow = c(8, 2))
MAplot(rawData[, 1:4], pairs = FALSE)
data.eset <- rma(data.raw)
data.exprs <- exprs(data.eset)
head(data.exprs)

pinfo <- getProbeInfo(data.raw)
head(pinfo)
fids <- pinfo$fid
data.exprs <- data.exprs[rownames(data.exprs) %in% fids, ]
nrow(data.exprs)

write.table(data.exprs  , file="data.exprs2.txt")
