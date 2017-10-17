##bioconductor和edgeR包的安装（若已安装可略过）
>source("http://bioconductor.org/biocLite.R")
>biocLite("edgeR")

##加载edgeR
>library(edgeR)
##读入数据
>x<-read.table("C:/Users/Administrator.SKY-20140313YIU/Desktop/report.csv",sep=",",header=TRUE)
##检查读入是否正确
> head(x)
##将数值和名称存入DGEList并赋值给y
> y<-DGEList(counts=x[,-1],genes=x[,1])
##设计design matrix
> group <- factor(c(1,1,2,2))
> design <- model.matrix(~group)
##推测dispersion（离散度）
> y<-estimateGLMCommonDisp(y,design,verbose=TRUE)
Disp = 0.31901 , BCV = 0.5648 
> y<-estimateGLMTrendedDisp(y, design)
> y<-estimateGLMTagwiseDisp(y, design)
##差异表达基因，to perform quasi-likelihood F-tests（法一）:
> fit <- glmQLFit(y,design)
> qlf <- glmQLFTest(fit,coef=2)
> topTags(qlf) ##前10个差异表达基因
 
##or 差异表达基因，to perform likelihood ratio tests:（法二）
>fit<-glmFit(y, design)
>lrt<-glmLRT(fit)
>topTags(lrt)##前10个差异表达基因
Coefficient:  group2 
                              genes      logFC    logCPM         F     PValue
3876                          21928  -5.996999  5.717751 10.932235 0.01328651
5949                          65969   7.016326  5.761071 10.927485 0.01329951
10998                        110454 -11.175373 10.568768 10.866447 0.01346814
3908                          21984 -11.275652  9.374238 10.390352 0.01488354
8218                          71775   2.689848  3.821093 10.438160 0.01507383
8156                          71660  -6.737531  7.650645 10.328192 0.01508235
9184                          75396   8.051870  6.474123 10.122665 0.01576475
1811                          15430  -5.152156  4.657149  9.994247 0.01621156
24    gi|165932359|ref|NM_009263.2|  -8.516929  9.621551  9.945919 0.01638397
5902                          64661  -9.605995  8.097023  9.896044 0.01656440
           FDR
3876  0.999835
5949  0.999835
10998 0.999835
3908  0.999835
8218  0.999835
8156  0.999835
9184  0.999835
1811  0.999835
24    0.999835
5902  0.999835
##绘制火山图
> summary(de<-decideTestsDGE(qlf))  ##qlf或可改为lrt
   [,1] 
-1     0
0  16268
1      0
> detags<-rownames(y)[as.logical(de)]
> plotSmear(qlf, de.tags=detags)
> abline(h=c(-4,4),col='blue') ##蓝线为2倍差异表达基因，差异表达的数据在qlf中
