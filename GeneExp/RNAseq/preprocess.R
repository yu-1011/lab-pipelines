################
# Describe: PreprocessForRNA-seqData
# Environment: R 3.3.2
# Author: Sihan Liu (SKLMG)
# Input:metData(clinicInfo),exprData(fpkm,tpm,count etc.)
# Output:Expression matrix
# R Package requirement: ggfortify, WGCNA, preprocessCore, lme4
# Data: 2017-06-05
###############

## Remove outliers based on network connectivity z-scores
par(mfrow=c(1,1), mar=c(5,4,2,2))
normadj <- (0.5+0.5*bicor(datExpr))^2 ## Calculate connectivity
netsummary <- fundamentalNetworkConcepts(normadj)
ku <- netsummary$Connectivity
z.ku <- (ku-mean(ku))/sqrt(var(ku))
plot(1:length(z.ku),z.ku,col=datMeta$Group,pch=19,main="Outliers", ylab="Connectivity Z score",xlab="")
legend("bottomright",legend = levels(datMeta$Group), col = 1:3,pch=19,cex=.7)
abline(h=-2, lty=2)
outliers = (z.ku < -2)
table(outliers)

#PCA
tpm <-read.csv(file="samplefiler6_qn.csv",header=T,sep="\t",row.names=1)
pdf("re-pca_qn.pdf")
pca <-prcomp(t(tpm),scale=T)
library(ggfortify)
group <-factor(c(rep("wuhan",30),rep("xiehe",121)))
colour_group <-cm.colors(length(unique(group)))
colour <-colour_group[as.numeric(factor(group))]
colour
group2 <-data.frame(group)
data2 <-cbind(t(tpm),group2)
autoplot(pca,data=data2,colour='group')
dev.off()

###
#Batch Clustering
library(WGCNA)
datMeta = read.csv("./metaData_RNAseqFreeze.csv",row.names=1, stringsAsFactors = T) ##corvirance
batch = model.matrix(~0+LibraryBatch+FlowcellBatch+RNAIsolationBatch+BrainBank, data=datMeta)
idx = match(rownames(batch), rownames(datMeta))
tree=hclust(dist((batch)),method="ward.D2")
colors=with(datMeta[idx,], cbind(labels2colors(LibraryBatch), 
                                 labels2colors(FlowcellBatch), 
                                 labels2colors(RNAIsolationBatch)))
                                 labels2colors(BrainBank)))
colors = cbind(colors, numbers2colors(PC[idx,]))            
plotDendroAndColors(tree,colors, groupLabels = c("LibraryBatch", "FlowcellBatch", "RNAIsolationBatch","BrainBank", paste0("PC",1:5)),cex.dendroLabels = .5,main="Batch Clustering")


cutsteps = seq(min(tree$height), 0.98*max(tree$height), by = (max(tree$height) - min(tree$height))/50)
r2out = numeric()
for(h in cutsteps) {
  c=as.factor(cutree(tree,h = h))
  s=0
  for(i in 1:5) s=s+summary(lm(PC[idx,1]~c))$adj.r.squared/5
  r2out=c(r2out,s)
}
plot(cutsteps, r2out)
datMeta$final_batch = NA
datMeta$final_batch[idx] = factor(cutree(tree, h=cutsteps[which.max(r2out)]))

dev.off()


##quantile
library(preprocessCore)
tpm <-read.csv(file="samplefiler6_qn.csv",header=T,sep="\t",row.names=1)
a=as.matrix(tpm)
b<-normalize.quantiles(a)
write.csv(b,file="tpm..qn.csv")

##corvirance select
library(lme4)
####### Edit these variables according to user defined parameters and the path to your data and data files names ##########

#myPath <- "C:/Users/WangKangli/Desktop/methylationdata"
#theGene_expression_file <- "RAW_62REPLICATES.TXT"
#theMethylation_file <- "methylation.csv"
#theExperiment_data_file <- "expinfo_tab_delimited.TXT"
#theMethylation_data_file <- "sample.csv"
pct_threshold = .5876 # Amount of variability desired to be explained by the principal components.  Set to match the results in book chapter and SAS code.  User can adjust this to a higher (>= 0.8) number but < 1.0

### In addition, be sure to modify the mixed linear model by adding the appropriate random effects terms in the model

################################################

theGEDFilePath = "adjust_151sample.csv"
theExpDataFilePath = "test2.txt"

########## Load data ##########

theDataMatrix <- read.csv(theGEDFilePath, header = TRUE, row.names = 1)
dataRowN <- nrow(theDataMatrix)
dataColN <- ncol(theDataMatrix)

########## Center the data (center rows) ##########
theDataMatrixCentered <- matrix(data = 0, nrow = dataRowN, ncol = dataColN)
theDataMatrixCentered_transposed = apply(theDataMatrix, 1, scale, center = TRUE, scale = FALSE)
theDataMatrixCentered = t(theDataMatrixCentered_transposed)

exp_design <- read.table(theExpDataFilePath, header = TRUE, row.names = 1, sep="\t")
exp_design <- exp_design[,c(1,2,3,4,5,6,7,8,9,10,11,12,13)]
expDesignRowN <- nrow(exp_design)
expDesignColN <- ncol(exp_design)
myColNames <- names(exp_design)


########## Compute correlation matrix ##########

theDataCor <- cor(theDataMatrixCentered)

########## Obtain eigenvalues ##########

eigenData <- eigen(theDataCor)
eigenValues = eigenData$values
ev_n <- length(eigenValues)
eigenVectorsMatrix = eigenData$vectors
eigenValuesSum = sum(eigenValues)
percents_PCs = eigenValues /eigenValuesSum 

########## Merge experimental file and eigenvectors for n components ##########

my_counter_2 = 0
my_sum_2 = 1
for (i in ev_n:1){
my_sum_2  = my_sum_2 - percents_PCs[i]
	if ((my_sum_2) <= pct_threshold ){
		my_counter_2 = my_counter_2 + 1
	}

}
if (my_counter_2 < 3){
	pc_n  = 3

}else {
	pc_n = my_counter_2 
}

# pc_n is the number of principal components to model

pc_data_matrix <- matrix(data = 0, nrow = (expDesignRowN*pc_n), ncol = 1)
mycounter = 0
for (i in 1:pc_n){
	for (j in 1:expDesignRowN){
	mycounter <- mycounter + 1
		pc_data_matrix[mycounter,1] = eigenVectorsMatrix[j,i]

	}
}

AAA <- exp_design[rep(1:expDesignRowN,pc_n),]

Data <- cbind(AAA,pc_data_matrix)

####### Edit these variables according to your factors #######

#Data$sex <- as.factor(Data$sex)
#Data$age <- as.factor(Data$age)
#Data$Index <- as.factor(Data$Index)

variables <- c(colnames(exp_design))
    for (i in 1:length(variables)) {
        Data$variables[i] <- as.factor(Data$variables[i])
    }

########## Mixed linear model ##########
op <- options(warn = (-1))
    effects_n = expDesignColN + choose(expDesignColN, 2) + 1
    randomEffectsMatrix <- matrix(data = 0, nrow = pc_n, ncol = effects_n)

model.func <- c()
    index <- 1
    for (i in 1:length(variables)) {
        mod = paste("(1|", variables[i], ")", sep = "")
        model.func[index] = mod
        index = index + 1
    }
    for (i in 1:(length(variables) - 1)) {
        for (j in (i + 1):length(variables)) {
            mod = paste("(1|", variables[i], ":", variables[j], 
                ")", sep = "")
            model.func[index] = mod
            index = index + 1
        }
    }
    function.mods <- paste(model.func, collapse = " + ")
    
    for (i in 1:pc_n) {
        y = (((i - 1) * expDesignRowN) + 1)
        funct <- paste("pc_data_matrix", function.mods, sep = " ~ ")
        Rm1ML <- lmer(funct, Data[y:(((i - 1) * expDesignRowN) + 
            expDesignRowN), ], REML = TRUE, control=lmerControl(optCtrl=list(maxfun=1e6),check.nobs.vs.nlev = "ignore",check.nobs.vs.rankZ = "ignore",check.nobs.vs.nRE="ignore"),verbose = FALSE, 
            na.action = na.omit)
        randomEffects <- Rm1ML
        randomEffectsMatrix[i, ] <- c(unlist(VarCorr(Rm1ML)), 
            resid = sigma(Rm1ML)^2)
    }
    effectsNames <- c(names(getME(Rm1ML, "cnms")), "resid")
########## Standardize Variance ##########

randomEffectsMatrixStdze <- matrix(data = 0, nrow = pc_n, ncol = effects_n)
for (i in 1:pc_n){
	mySum = sum(randomEffectsMatrix[i,])
	for (j in 1:effects_n){
		randomEffectsMatrixStdze[i,j] = randomEffectsMatrix[i,j]/mySum	
	}
}

########## Compute Weighted Proportions ##########

randomEffectsMatrixWtProp <- matrix(data = 0, nrow = pc_n, ncol = effects_n)
for (i in 1:pc_n){
	weight = eigenValues[i]/eigenValuesSum
	for (j in 1:effects_n){
		randomEffectsMatrixWtProp[i,j] = randomEffectsMatrixStdze[i,j]*weight
	}
}

########## Compute Weighted Ave Proportions ##########

randomEffectsSums <- matrix(data = 0, nrow = 1, ncol = effects_n)
randomEffectsSums <-colSums(randomEffectsMatrixWtProp)
totalSum = sum(randomEffectsSums)
randomEffectsMatrixWtAveProp <- matrix(data = 0, nrow = 1, ncol = effects_n)

for (j in 1:effects_n){
	randomEffectsMatrixWtAveProp[j] = randomEffectsSums[j]/totalSum 	
	
}
c(84,89,83,87,1,90,58,57,81,61,92)
pdf("pvca_new.pdf",height=3.5)
bp <- barplot(randomEffectsMatrixWtAveProp[,c(84,89,83,87,1,90,58,57,81,61,92)],  xlab = "Effects", ylab = "Weighted average proportion variance", ylim= c(0,0.5),col = c("Red"), las=2,title="pvca of BE")

axis(1, at = bp, labels = effectsNames[c(84,89,83,87,1,90,58,57,81,61,92)], xlab = "Effects", las=2,cex.axis=0.5)

## replace the above code of "axis(1, at = bp, labels = effectsNames, xlab = "Effects", cex.axis = 0.5, las=2)" if you want rotate the x axis labels. 
## text(bp, par("usr")[3]-0.02, srt = 45, adj = 1,labels = effectsNames, xpd = TRUE,cex=0.8)

values = randomEffectsMatrixWtAveProp[,c(84,89,83,87,1,90,58,57,81,61,92)]
new_values = round(values , 3)
text(bp,randomEffectsMatrixWtAveProp[,c(84,89,83,87,1,90,58,57,81,61,92)],labels = new_values, pos=3,cex=0.8) # place numbers on top of bars 
dev.off()


##corvirance correction
Y <- read.table("tpm.qn.2.txt",header=T,row.names=1,sep="\t")
datExp<-log2(0.001+Y)
data<-cbind(t(datExp),datMeta)
test<- lm(t(datExp) ~ Group + RNAseqBatch+CauseDeath:ageDeath+RNAIsolationBatch:ageDeath+RNAIsolationBatch:RIN+Group:RNAseqBatch+Group:LibraryBatch+seqPC4:seqPC5+LibraryBatch+Group:RNAIsolationBatch,data = data)
test.cor<-t(test$residuals)+rowMeans(datExp)
write.csv(test.cor,file="adjust.csv")
