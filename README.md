# lab-pipelines
这里是所有分析汇总流程  
一级文件夹名字为分类名称，DEG/GENOME/Integration  
二级文件夹名字为方法名称  
到时这里会是一个目录的啦  
[回到顶部](#lab-pipelines)  
[目录](##目录)  
[格式说明](##格式说明)  
 

## 目录
DEG 
Cell type adjustment—
Connectivity Outlier
Covariates—VOOM
Network—WGCNA
Module conservation—Zsummary
Functional enrichment—DAVID, Gorilla

Methylation 

DMP, DMR—
Feature—ChromHMM
GWAS overlap—INRICH, 


Call Variants
Imputation
LDSR—LD regression
π1 statistics
GWAS  
CVAS  
RVAS  

QTL—matrix eQTL

Causal relationship—SMR, NEO, CIT, metaSCAN  

Evolution
HAR
Dn/Ds  

## 格式说明
每个文件夹至少包括两个部分  
1.code  
2.说明文档

请在code前加上作为代码描述部分  
\####  
\# Describe: 代码作用  
\# Environment:语言环境  
\# Author: 作者  
\# Input:输入文件  
\# Output:输出文件  
\# Package requirement:相关依赖包  
\# Data: 时间  
\####

例:    
\####  
\# Describe: use affy packages to preprocess .CEL  
\# Environment: R 3.3.2  
\# Author: Yu Chen (SKLMG)  
\# Input:.CEL file  
\# Output:Expression matrix  
\# R Package requirement: affy, tcltk  
\# Data: 2017-08-15  
\####  

-------------
对于说明文档 
请加上流程的具体作用以便学习

