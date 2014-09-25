setwd("C:/Users/Courtney/Dropbox/Ober Lab/Fertility")

##Load lumi
library(lumi)

##Create object with name of data file:
alldata = c('totaldataclean.csv')

##Extract raw data from lumi file and preparing for removing bad probes:
data.lumi.all = lumiR.batch(alldata, lib.mapping=NULL, convertNuID=F,annotationColumn=c('ACCESSION', 'SYMBOL', 'PROBE_SEQUENCE', 'PROBE_START', 'CHROMOSOME', 'PROBE_CHR_ORIENTATION', 'PROBE_COORDINATES', 'DEFINITION','PROBE_ID'))


### NORMALIZATION: log2 stabilized and quantile normalization ###
data.norm.all <- lumiExpresso(data.lumi.all, bg.correct=TRUE, bgcorrect.param=list(method='forcePositive'), variance.stabilize=TRUE, varianceStabilize.param = list(method="log2"), normalize=TRUE, normalize.param=list(method="quantile"), QC.evaluation=TRUE, QC.param=list(), verbose=TRUE)
expr_quant.all <- data.norm.all@assayData$exprs

#Create Supplementary normalization figure 
samplekey = read.csv('AllColumnNames_forsup.csv')
rem = c(99,  107,	28,	75,	37,	68,	48,	79,	40,	51,	77,	106,	72,	108,	11,	21,	39,	63,	83,	112,	27,	54,	60,	76)
data.lumi.clean = data.lumi.all[,-rem]
clean.key = samplekey[-rem,]
sampleNames(data.lumi.clean) = clean.key$SampleID
data.norm.clean = data.norm.all[,-rem]
sampleNames(data.norm.clean) = clean.key$SampleID

boxplot(data.lumi.clean, main = "Pre-normalization Microarray Intensity", ylab = "Intensity")
boxplot(data.norm.clean, main = "Post-normalization Microarray Intensity", ylab = "Intensity")
plot(data.norm.clean, what='sampleRelation')

###Find the column that is lumi_ID in feature data usually this column
head(data.norm.all@featureData[[1]])
##[1] "ILMN_1343291" "ILMN_1343295" "ILMN_1651199" "ILMN_1651209" "ILMN_1651210" "ILMN_1651221"
###Convert expr_quant rownames to
rownames(expr_quant.all)=data.norm.all@featureData[[1]]

#Filter out probes that have a p<.01 in more than 70% of the sample
pval = read.table('totaldatacleanpval2.txt', as.is=T, header=T)
pvalnorowID = pval[,-1]
detect_quant.all= rowSums(pval<0.05)
detect.ind.all <- which(detect_quant.all > 79)
expr_quant.clean = expr_quant.all[detect.ind.all,]
expr_quant.all = expr_quant.clean
rm(expr_quant.clean)

##Look at plots of array data (boxplot, density, etc) :
plot(data.lumi.all, what='boxplot')
plot(data.lumi.all, what='density')
plot(data.norm.all, what='boxplot')
plot(data.norm.all, what='density')

##Check that replicates are most related
plot(data.norm.all, what='sampleRelation')

##Add information about sample names:
samplenames = read.csv('AllColumnNames.csv')
colnames(expr_quant.all) = samplenames[,1]

###Subset expression by Darren good probes
goodprobes= read.table('HT-12v4_Probes_inhg19EnsemblGenes_NoCEUHapMapSNPs_Stranded.txt', header=T)
probes = goodprobes$probeID
## Convert from factor to character
probes = as.character(goodprobes$probeID)
expr_quant.all.clean = expr_quant.all[rownames(expr_quant.all) %in% probes, ]
# dim(expr_quant.all.clean)  
# 12,383 probes of the expressed 17,208 "good" probes are in this data set
expr_quant.all= expr_quant.all.clean
remove(expr_quant.all.clean)

##Load in Covariates
#Converted chip and batch to a factor so they are categorical 
samplenames = read.csv('AllColumnNames.csv')
batch = samplenames[,2]
chip = samplenames[,3]
age =samplenames[,5]
bmi = samplenames[,6]
prmisc = samplenames[,7]
hypoth = samplenames[,8]
race = samplenames[,9]
class = samplenames[,10]
season = samplenames[,11]
#Converted categorical covariates to a factor so they are levels 
batch.f = as.factor(batch)
chip.f = as.factor(chip)
hypoth.f = as.factor(hypoth)
race.f = as.factor(race)
season.f = as.factor(season)
class.f = as.factor(class)
#Converted numerical covariates to a numeric so they are continuous
prmisc.num = as.numeric(prmisc)


##To look for correlations between covariates get the R^2 of the Pearson correlation
cor(age, bmi)^2
cor(age, prmisc.num)^2
cor(bmi, prmisc.num)^2
cor(chip, batch)^2
cor(chip, bmi)^2
cor(chip, age)^2
cor(chip, prmisc.num)^2
cor(batch, prmisc.num)^2
cor(batch, age)^2
cor(batch, bmi)^2
cor(age, hypoth)^2
cor(bmi, hypoth)^2
cor(chip, hypoth)^2
cor(batch, hypoth)^2
cor(prmisc.num, hypoth)^2
rm(batch, chip, hypoth, season, class, race, prmisc)

#Make a heatmap with raw data
library("gplots")
cor <- cor(expr_quant.all,method="pearson", use="complete.obs")
heatmap.2(cor, main="Probe expression correlation", key=T, revC=T, density.info="none", trace="none")

#Make PC plots with raw data
library(TeachingDemos)
op <- par(mfrow = c(3,3), ## split region
          oma = c(5,0,4,0) + 0.1, ## create outer margin
          mar = c(5,4,2,2) + 0.1) ## shrink some margins
tmp1 <- cnvrt.coords( 0.5, 0, input='plt' )$tdev
cor <- cor(expr_quant.all, method="pearson", use="complete.obs")
title.PC = "PCA of Gene Exp"
sum.PC <- prcomp(na.omit(cor))
sumsum <- summary(sum.PC)
#prints out plots in c(#rows, #columns)
color = samplenames[,13]
par(mfrow = c(2,2),oma=c(0,0,2,0)) 
plot(c(1:112),sum.PC$rotation[,1],col=color, xlab="Index of Samples",pch = 20, ylab=paste("PC 1 -",(sumsum$importance[2,1]*100),"% of variance",sep=" "),main=title.PC)
text(c(1:112),sum.PC$rotation[,1], samplenames[,1], cex = .5, pos=3)   
for(i in 2:4) {
  plot(sum.PC$rotation[,1], sum.PC$rotation[,i], col=color,pch=20,main=title.PC, xlab=paste("PC 1 -", (sumsum$importance[2,1]*100),"% of variance", sep=" "), ylab=paste("PC",i,"-",(sumsum$importance[2,i]*100),"% of variance", sep=" "))
  }  
par(op)
#All plotting is clipped to the device region
par(xpd=NA)
## To check for BMI and age effects
par(mfrow = c(1,1))
sampleorder3<-read.table("sampleorder3.txt")
plot(c(1:20),sum.PC$rotation[,1],col=color,xlab="Index of Samples",pch = 20, ylab=paste("PC 1 -",(sumsum$importance[2,1]*100),"% of variance",sep=" "),main=title.PC)
text(c(1:20),sum.PC$rotation[,1], sampleorder3[,4], cex = .5, pos=3)
text(c(1:20),sum.PC$rotation[,1], sampleorder3[,2], cex = .5, pos=1)

##To see if covariates are correlated with a PC (looking at PC1-7)
sum.PC <- prcomp(na.omit(cor))
summary(lm(sum.PC$rotation[, 1:7] ~ chip.f))
summary(lm(sum.PC$rotation[, 1:7] ~ season.f))
summary(lm(sum.PC$rotation[, 1:7] ~ bmi))
summary(lm(sum.PC$rotation[, 1:7] ~ age))
summary(lm(sum.PC$rotation[, 1:7] ~ hypoth.f))
summary(lm(sum.PC$rotation[, 1:7] ~ race.f))
summary(lm(sum.PC$rotation[, 1:7] ~ prmisc.num))
summary(lm(sum.PC$rotation[, 1:7] ~ class.f))
summary(lm(sum.PC$rotation[, 1:7] ~ batch.f))

## Automated PCA analysis
#This could be cleaned up using lapply or something like that, but here is a loop and code save pca regression pvalues and adj R sqs to a file:

# put all of your covariates into a list
covars<-list(chip.f,season.f,hypoth.f, race.f, bmi, age, prmisc.num, class.f, batch.f)

# this function takes your pca object, the list of covars, and the number of PCS you want to regress, the result is a matrix, N rows is the # of covariates, N columns is 2* the number of PCS
# each PC has 2 columns in the resulting matrix, column 1 is the pvalue, column 2 is the model adj R2, this is the % of variation in PC1 explained by your covariate.
npcs = 4
sum.PC <- prcomp(na.omit(cor))
results<-c()
for (f in covars) {
  for (i in 1:npcs)
  {
    s = summary(lm(sum.PC$rotation[,i]~f));
    results<-c(results,pf(s$fstatistic[[1]],
                          s$fstatistic[[2]],s$fstatistic[[3]], lower.tail = FALSE),
               s$adj.r.squared)
  }
}
resultsM<-matrix(nrow = length(covars), ncol = 2*npcs, data =
                   results, byrow = TRUE)
row.names(resultsM) = c("chip","season","hypoth", "race", "bmi", "age", "prmisc", "class", "batch")
colnames(resultsM) =c("PC1 pval", "PC1 R2","PC2 pval", "PC2 R2","PC3 pval", "PC3 R2","PC4 pval", "PC4 R2")
resultsM

##Regress out batch and chip and add the intercept back in
batch.chip.residual.int = matrix(nrow= nrow(expr_quant.all), ncol = ncol(expr_quant.all))
rownames(batch.chip.residual.int) = rownames(expr_quant.all)
colnames(batch.chip.residual.int) = colnames(expr_quant.all)
for (i in 1:nrow(expr_quant.all)) {
  model= lm(expr_quant.all[i,]~ batch.f + chip.f)
  batch.chip.residual.int[i,] = resid(model) + model$coefficients[1]
}
cor.bc.int <- cor(batch.chip.residual.int,method="pearson", use="complete.obs")
heatmap.2(cor.bc.int, main="Probe expression after batch + chip control", key=T, revC=T, density.info="histogram", trace="none", dendrogram = "column")
sum.PC.bc.int <- prcomp(na.omit(cor.bc.int))
bc.resid.int = batch.chip.residual.int
remove(batch.chip.residual.int)

#To get the correlation between replicates vs non-replicates
cor.bc <- cor(bc.resid.int,method="spearman", use="complete.obs")
replicates=c()
non.replicates=c()
names = samplenames[,12]
for(i in 1:111){
  for(j in (i+1):112){
    if(names[i]==names[j]){
      replicates=c(replicates,cor.bc[i,j])
    }
    else{
      non.replicates=c(non.replicates,cor.bc[i,j])
    }
  }
}
boxplot.n(replicates,non.replicates, main = "Correlation of Samples", ylab = "Correlation", xlab = "Replicates vs Non-Replicates")

##Remove the after progesterone treated samples to see how that changes the sample correlation
afterprog = c(28,37,40,48,51,59,68,75,77,79,99,106,107)
clean.bc.resid.int=bc.resid.int[,- afterprog]
samplenames.clean = read.csv('CleanColNames.csv')
#To get the correlation between replicates vs non-replicates
cor.bc.clean <- cor(clean.bc.resid.int,method="spearman", use="complete.obs")
replicates=c()
non.replicates=c()
names = samplenames.clean[,1]
for(i in 1:98){
  for(j in (i+1):99){
    if(names[i]==names[j]){
      replicates=c(replicates,cor.bc.clean[i,j])
    }
    else{
      non.replicates=c(non.replicates,cor.bc.clean[i,j])
    }
  }
}
boxplot.n(replicates,non.replicates, main = "Correlation of Samples", ylab = "Correlation", xlab = "Replicates                                              Non-Replicates")

##Remove the after progesterone treated samples and fertile controls to see how thar changes the sample correlation
afterprogfc = c(28,37,40,48,51,59,68,75,77,79,99,106,107,21,27,39,54,63,83,112)
clean.bc.resid.int=bc.resid.int[,- afterprogfc]
samplenames.clean = samplenames[- afterprogfc,]
#To get the correlation between replicates vs non-replicates
cor.bc.clean <- cor(clean.bc.resid.int,method="spearman", use="complete.obs")
replicates=c()
non.replicates=c()
names = samplenames.clean[,12]
for(i in 1:91){
  for(j in (i+1):92){
    if(names[i]==names[j]){
      replicates=c(replicates,cor.bc.clean[i,j])
    }
    else{
      non.replicates=c(non.replicates,cor.bc.clean[i,j])
    }
  }
}
boxplot.n(replicates,non.replicates, main = "Correlation of Samples", ylab = "Correlation", xlab = "Replicates                                              Non-Replicates")
summary(replicates)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.9338  0.9613  0.9655  0.9686  0.9831  0.9875 
summary(non.replicates)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.8604  0.9336  0.9472  0.9440  0.9574  0.9856 

##To see if covariates are correlated with a PC (looking at PC1-7)
summary(sum.PC.bc.int)
summary(lm(sum.PC.bc.int$rotation[, 1:7] ~ chip.f))
summary(lm(sum.PC.bc.int$rotation[, 1:7] ~ season.f))
summary(lm(sum.PC.bc.int$rotation[, 1:7] ~ bmi))
summary(lm(sum.PC.bc.int$rotation[, 1:7] ~ age))
summary(lm(sum.PC.bc.int$rotation[, 1:7] ~ hypoth.f))
summary(lm(sum.PC.bc.int$rotation[, 1:7] ~ race.f))
summary(lm(sum.PC.bc.int$rotation[, 1:7] ~ prmisc.num))
summary(lm(sum.PC.bc.int$rotation[, 1:7] ~ class.f))
summary(lm(sum.PC.bc.int$rotation[, 1:7] ~ batch.f))

##Add information about sample names:
samplenames = read.csv('AllColumnNames.csv')
colnames(bc.resid.int) = samplenames[,1]
##Check sample names
head(bc.resid.int)

## Combine replicates
##The after progesterone treated samples are not included
## afterprog = c(28,37,40,48,51,59,68,75,77,79,99,106,107)
m436 = apply(bc.resid.int[, c(12,26,30)],1,mean)
m306 = apply(bc.resid.int[, c(16,95)],1,mean)
m334 = apply(bc.resid.int[, c(17,69)],1,mean)
mFCP01 = apply(bc.resid.int[, c(21,39,63)],1,mean)
m336 = apply(bc.resid.int[, c(23,45)],1,mean)
m274 = apply(bc.resid.int[, c(24,56)],1,mean)
mD0143 = apply(bc.resid.int[, c(25,34)],1,mean)
m334 = apply(bc.resid.int[, c(17,69)],1,mean)
mFCP06 = apply(bc.resid.int[, c(27,54)],1,mean)
m165 = apply(bc.resid.int[, c(49,64)],1,mean)
m329 = apply(bc.resid.int[, c(91,29,41)],1,mean)
m226 = apply(bc.resid.int[, c(84,32)],1,mean)
m309 = apply(bc.resid.int[, c(33,71)],1,mean)
m347 = apply(bc.resid.int[, c(104,43)],1,mean)
m451 = apply(bc.resid.int[, c(35,67)],1,mean)
m297 = apply(bc.resid.int[, c(88,36)],1,mean)
m235 = apply(bc.resid.int[, c(109,82)],1,mean)
m370 = apply(bc.resid.int[, c(87,38)],1,mean)
m348 = apply(bc.resid.int[, c(111,42)],1,mean)
m128 = apply(bc.resid.int[, c(44,59)],1,mean)
m325 = apply(bc.resid.int[, c(80,89)],1,mean)
m441 = apply(bc.resid.int[, c(50,86)],1,mean)
m216 = apply(bc.resid.int[, c(105,52)],1,mean)
m363 = apply(bc.resid.int[, c(53,70)],1,mean)
m305 = apply(bc.resid.int[, c(55,92)],1,mean)
m453 = apply(bc.resid.int[, c(60,76)],1,mean)
m387 = apply(bc.resid.int[, c(3,110,61)],1,mean)
m257 = apply(bc.resid.int[, c(65,73)],1,mean)
m156 = apply(bc.resid.int[, c(66,74)],1,mean)
mD0124 = apply(bc.resid.int[, c(108,72)],1,mean)
m321 = apply(bc.resid.int[, c(94,78)],1,mean)
m356 = apply(bc.resid.int[, c(47,81)],1,mean)
mFCP04 = apply(bc.resid.int[, c(112,83)],1,mean)
m385 = apply(bc.resid.int[, c(100,85)],1,mean)
m423 = apply(bc.resid.int[, c(103,90)],1,mean)
m391 = apply(bc.resid.int[, c(96,97)],1,mean)
m395 = apply(bc.resid.int[, c(98,101)],1,mean)
m408 = apply(bc.resid.int[, c(102,57)],1,mean)
##Pull all the rows together to get 59 final samples
final_probe= cbind(bc.resid.int[, c(1:2,4:11,13:15,18:20,22,31,93,46,58,62)],m306, m436, m334,mFCP01, m336, m274,mD0143,mFCP06,m165,m329,m226,m309,m451,m297,m235,m370,m348,m347,m128,m325,m441, m216,m363,m305,m453,m387,m257,m156,mD0124,m321,m356,mFCP04,m385,m423,m391,m395,m408)

#If you want to look at sample correlation before avg replicates
## Finding the unique gene names matching probes to gene names using Darren's good probe list
gene_names=c()
for(i in 1:dim(bc.resid.int)[1]){
  gene_names=c(gene_names,as.vector(goodprobes[as.vector(goodprobes[,4])==row.names(bc.resid.int)[i],8]))
}

symbolsUniq = unique(gene_names)
length(symbolsUniq)


## This loop will give the most 3' value for multiple probes within the same gene. In the end you get a simple file with all genes that are expressed with the corresponding mean intensity expression levels across its different probes.
expr_gene = matrix(NA, ncol=112, nrow=length(unique(gene_names)))
i=0
for(gene in unique(gene_names)){
  i = i+1
  
  currRows = which(gene_names == gene)
  if(length(currRows)>1){
    if(goodprobes[currRows[1],6]=="+"){
      keepRow = currRows[which.max(goodprobes[currRows,2])]
    }
    else{
      keepRow = currRows[which.min(goodprobes[currRows,2])]
    }
  }
  else{
    keepRow=currRows[1]
  }
  expr_gene[i,] = bc.resid.int[keepRow,]
  
} 
dim(expr_gene)
rownames(expr_gene) = unique(gene_names)
colnames(expr_gene) = colnames(bc.resid.int)

#To get the correlation between replicates vs non-replicates
cor.gene <- cor(expr_gene2,method="spearman", use="complete.obs")
replicates=c()
non.replicates=c()
names = samplenames[,12]
for(i in 1:111){
  for(j in (i+1):112){
    if(names[i]==names[j]){
      replicates=c(replicates,cor.gene[i,j])
    }
    else{
      non.replicates=c(non.replicates,cor.gene[i,j])
    }
  }
}
boxplot.n(replicates,non.replicates, main = "Correlation of Samples", ylab = "Correlation", xlab = "Replicates                                              Non-Replicates")
## Finding the unique gene names matching probes to gene names using Darren's good probe list
gene_names=c()
for(i in 1:dim(final_probe)[1]){
  gene_names=c(gene_names,as.vector(goodprobes[as.vector(goodprobes[,4])==row.names(final_probe)[i],8]))
}

symbolsUniq = unique(gene_names)
length(symbolsUniq)

#If you want to merge at probes use this to create final gene data otherwise avg gene data below
## This loop will give the most 3' value for multiple probes within the same gene. In the end you get a simple file with all genes that are expressed with the corresponding mean intensity expression levels across its different probes.
expr_gene = matrix(NA, ncol=59, nrow=length(unique(gene_names)))
i=0
for(gene in unique(gene_names)){
  i = i+1
  
  currRows = which(gene_names == gene)
  if(length(currRows)>1){
    if(goodprobes[currRows[1],6]=="+"){
      keepRow = currRows[which.max(goodprobes[currRows,2])]
    }
    else{
      keepRow = currRows[which.min(goodprobes[currRows,2])]
    }
  }
  else{
    keepRow=currRows[1]
  }
  expr_gene[i,] = final_probe[keepRow,]
  
} 
dim(expr_gene)
rownames(expr_gene) = unique(gene_names)
colnames(expr_gene) = colnames(final_probe)

## It doesn't matter when you combine the replicates because you take the same probe in either case
## Combine the replicates after converting probe to gene
##The after progesterone treated samples are not included
## afterprog = c(28,37,40,48,51,59,68,75,77,79,99,106,107)
m436 = apply(expr_gene[, c(12,26,30)],1,mean)
m306 = apply(expr_gene[, c(16,95)],1,mean)
m334 = apply(expr_gene[, c(17,69)],1,mean)
mFCP01 = apply(expr_gene[, c(21,39,63)],1,mean)
m336 = apply(expr_gene[, c(23,45)],1,mean)
m274 = apply(expr_gene[, c(24,56)],1,mean)
mD0143 = apply(expr_gene[, c(25,34)],1,mean)
m334 = apply(expr_gene[, c(17,69)],1,mean)
mFCP06 = apply(expr_gene[, c(27,54)],1,mean)
m165 = apply(expr_gene[, c(49,64)],1,mean)
m329 = apply(expr_gene[, c(91,29,41)],1,mean)
m226 = apply(expr_gene[, c(84,32)],1,mean)
m309 = apply(expr_gene[, c(33,71)],1,mean)
m347 = apply(expr_gene[, c(104,43)],1,mean)
m451 = apply(expr_gene[, c(35,67)],1,mean)
m297 = apply(expr_gene[, c(88,36)],1,mean)
m235 = apply(expr_gene[, c(109,82)],1,mean)
m370 = apply(expr_gene[, c(87,38)],1,mean)
m348 = apply(expr_gene[, c(111,42)],1,mean)
m128 = apply(expr_gene[, c(44,59)],1,mean)
m325 = apply(expr_gene[, c(80,89)],1,mean)
m441 = apply(expr_gene[, c(50,86)],1,mean)
m216 = apply(expr_gene[, c(105,52)],1,mean)
m363 = apply(expr_gene[, c(53,70)],1,mean)
m305 = apply(expr_gene[, c(55,92)],1,mean)
m453 = apply(expr_gene[, c(60,76)],1,mean)
m387 = apply(expr_gene[, c(3,110,61)],1,mean)
m257 = apply(expr_gene[, c(65,73)],1,mean)
m156 = apply(expr_gene[, c(66,74)],1,mean)
mD0124 = apply(expr_gene[, c(108,72)],1,mean)
m321 = apply(expr_gene[, c(94,78)],1,mean)
m356 = apply(expr_gene[, c(47,81)],1,mean)
mFCP04 = apply(expr_gene[, c(112,83)],1,mean)
m385 = apply(expr_gene[, c(100,85)],1,mean)
m423 = apply(expr_gene[, c(103,90)],1,mean)
m391 = apply(expr_gene[, c(96,97)],1,mean)
m395 = apply(expr_gene[, c(98,101)],1,mean)
m408 = apply(expr_gene[, c(102,57)],1,mean)
##Pull all the rows together to get 59 final samples
final_gene= cbind(expr_gene[, c(1:2,4:11,13:15,18:20,22,31,93,46,58,62)],m306, m436, m334,mFCP01, m336, m274,mD0143,mFCP06,m165,m329,m226,m309,m451,m297,m235,m370,m348,m347,m128,m325,m441, m216,m363,m305,m453,m387,m257,m156,mD0124,m321,m356,mFCP04,m385,m423,m391,m395,m408)

##expr_gene = final_gene!!

## To make heatmap (with key and title)
library(gplots)
cor.q <- cor(expr_quant,method="pearson", use="complete.obs")
heatmap.2(cor.q, main="Menopause Gene Expression", key=T, revC=T, density.info="none", trace="none")

## To make heatmap (with key and title) using 3' most probe of the gene
library(gplots)
cor.q2 <- cor(final_gene,method="pearson", use="complete.obs")
heatmap.2(cor.q2, main="Endometrial Gene Expression", key=T, revC=T, density.info="none", trace="none")


## To make another version of heatmap with color codes for menopause class where red = pre-menopause ##and blue = post-menopause
library(gplots)
colnames(expr_quant) = colnames(expr_gene_q)
covarmi= c("red", " blue", " red", " red", " red", " blue", " blue", " red", 
           " blue", " blue", " red", " blue", " red", " red", " red", " blue", 
           " blue", " red", " blue", " blue")
heatmap.2(cor.q, main= "Correlation", key=T, revC=T,ColSideColors=covarmi, density.info="none", trace="none")
##Using probe averages
heatmap.2(cor.q2, main= "Correlation", key=T, revC=T,ColSideColors=covarmi, density.info="none", trace="none")

##Make a heatmap using averages of samples and omitting 107932
cor.q3 <- cor(finalexpr,method="pearson", use="complete.obs")
heatmap.2(cor.q3, main= "Correlation of Individuals", key=T, revC=T,ColSideColors=colorfinal, density.info="none", trace="none")

##PCA all probes
plot_colors<-c("black","purple")
library(TeachingDemos)
op <- par(mfrow = c(3,3), ## split region
          oma = c(5,0,4,0) + 0.1, ## create outer margin
          mar = c(5,4,2,2) + 0.1) ## shrink some margins
tmp1 <- cnvrt.coords( 0.5, 0, input='plt' )$tdev
cor.m <- cor(expr_quant, method="pearson", use="complete.obs")
title.PC = "PCA of Gene Exp"
sum.PC <- prcomp(na.omit(cor.m))
sumsum <- summary(sum.PC)
#Red dots are pre-menopause and blue are post-menopause
color = covarmi
#prints out plots in c(#rows, #columns)
par(mfrow = c(2,2),oma=c(0,0,2,0)) 
plot(c(1:112),sum.PC$rotation[,1],xlab="Index of Samples",pch = 20, ylab=paste("PC 1 -",(sumsum$importance[2,1]*100),"% of variance",sep=" "),main=title.PC)
text(c(1:112),sum.PC$rotation[,1], samplenames[,2], cex = .5, pos=3)   
for(i in 2:4) {
  plot(sum.PC$rotation[,1], sum.PC$rotation[,i], pch=20,main=title.PC, xlab=paste("PC 1 -", (sumsum$importance[2,1]*100),"% of variance", sep=" "), ylab=paste("PC",i,"-",(sumsum$importance[2,i]*100),"% of variance", sep=" "))
  text(sum.PC$rotation[,1], sum.PC$rotation[,i], samplenames[,2], cex = .5, pos=3)}  
par(op)
#All plotting is clipped to the device region
par(xpd=NA)
## To check for BMI and age effects
par(mfrow = c(1,1))
sampleorder3<-read.table("sampleorder3.txt")
plot(c(1:20),sum.PC$rotation[,1],col=color,xlab="Index of Samples",pch = 20, ylab=paste("PC 1 -",(sumsum$importance[2,1]*100),"% of variance",sep=" "),main=title.PC)
text(c(1:20),sum.PC$rotation[,1], sampleorder3[,4], cex = .5, pos=3)
text(c(1:20),sum.PC$rotation[,1], sampleorder3[,2], cex = .5, pos=1)


## To see if BMI or age are clustering omitting sample 107932
par(mfrow = c(1,1))
plot(c(1:18),sum.PC3$rotation[,1],col=color2,xlab="Index of Samples",pch = 20, ylab=paste("PC 1 -",(sumsum3$importance[2,1]*100),"% of variance",sep=" "),main=title.PC3)
text(c(1:18),sum.PC3$rotation[,1], sampleorder2[,4], cex = .5, pos=3)
text(c(1:18),sum.PC3$rotation[,1], sampleorder2[,2], cex = .5, pos=1)
## For multiple PCA's
par(mfrow = c(2,2),oma=c(0,0,2,0)) 
plot(c(1:18),sum.PC3$rotation[,1],col=color2,xlab="Index of Samples",pch = 20, ylab=paste("PC 1 -",(sumsum3$importance[2,1]*100),"% of variance",sep=" "),main=title.PC3)
text(c(1:18),sum.PC3$rotation[,1], sampleorder2[,5], cex = .5, pos=3)   
for(i in 2:4) { plot(sum.PC3$rotation[,1], sum.PC3$rotation[,i], pch=20, main=title.PC3, xlab=paste("PC 1 -", (sumsum$importance[2,1]*100),"% of variance", sep=" "), ylab=paste("PC",i,"-",(sumsum$importance[2,i]*100),"% of variance", sep=" "))
                text(sum.PC3$rotation[,1], sum.PC3$rotation[,i], samplenamesclean[,2], cex = .5, pos=3)}

##Extraction PC data from PCA
PCS1 = sum.PC$rotation
View(PCS1)
PCS2 =sum.PC2$rotation
View(PCS2)
PCS3 =sum.PC3$rotation
View(PCS3)
#Check this
#PCs are columns, inds are rows, like this:
#               PC1          PC2          PC3           PC4          PC5
#87_MH  -0.08176976 -0.132745281 -0.069088998  0.1525655116 -0.106039703
#125_MH -0.18216685 -0.083451452  0.079887299 -0.1936390876 -0.204922215
#168_ML -0.01938816 -0.205549632 -0.027594568 -0.0324146775 -0.223454860
#209_FL -0.12189288 -0.129940372 -0.052735639  0.0120791405 -0.178501057


#Check technical replicates 
par(mfrow = c(2,2),oma=c(0,0,1,0))
plot(expr_quant.all[1:1000,3],expr_quant.all[1:1000,111], main = "387")
plot(expr_quant.all[1:1000,3],expr_quant.all[1:1000,12], main = "387+Random")


#Averages of technical replicates
one<-apply(expr_gene_o[,1:2],1,mean)
two<-apply(expr_gene_o[,3:4],1,mean)
three<-apply(expr_gene_o[,5:6],1,mean)
four<-apply(expr_gene_o[,7:8],1,mean)
five<-apply(expr_gene_o[,9:10],1,mean)
six<-apply(expr_gene_o[,11:12],1,mean)
seven<-apply(expr_gene_o[,13:14],1,mean)
eight<-apply(expr_gene_o[,15:16],1,mean)
nine<-apply(expr_gene_o[,17:18],1,mean)
ten<-apply(expr_gene_o[,19:20],1,mean)
averages_expr_gene<-cbind(one,two,three,four, five, six, seven, eight, nine, ten)

#Looking for corellation of PC's to covariates
#Look at p-value
summary(lm(sum.PC.final$rotation[1:9, 1:7] ~ prot.f))
summary(lm(sum.PC.final$rotation[-2, 1:7] ~ bmi))
#Plotting interesting things
plot(c(1:9),sum.PC.final$rotation[,1],col=colorfinal,xlab="Index of Samples",pch = 20, ylab=paste("PC 1 -",(sumsumfinal$importance[2,1]*100),"% of variance",sep=" "),main= "PCA 1 with Location, Cell, Age, and Protocol")
text(c(1:9),sum.PC.final$rotation[,1], loc, cex = .5, pos=3)
text(c(1:9),sum.PC.final$rotation[,1], cell, cex = .5, pos=2)
text(c(1:9),sum.PC.final$rotation[,1], prot, cex = .5, pos=1)
text(c(1:9),sum.PC.final$rotation[,1], age, cex = .5, pos=4)
plot(c(1:9),sum.PC.final$rotation[,5],col=colorfinal,xlab="Index of Samples",pch = 20, ylab=paste("PC 5 -","8.06% of variance",sep=" "),main="PCA 5 with Date")
text(c(1:9),sum.PC.final$rotation[,5], date, cex = .5, pos=3)
plot(c(1:9),sum.PC.final$rotation[,4],col=colorfinal,xlab="Index of Samples",pch = 20, ylab=paste("PC 4 -","10.04% of variance",sep=" "),main="PCA 4 with Menopause Status")
text(c(1:9),sum.PC.final$rotation[,4], fsampleorder[,5], cex = .5, pos=3)

#To look for correlations between the covariates -> get R2 of a Pearson correlation test example:
cor(cell, date)^2
cor(age, prot)^2
cor(prot[-2], bmi)^2

# T test
pvalmenp = vector()
for (i in 1:length(averages_expr_gene[,1])) {
  pvalmenp = c(pvalmenp, t.test(averages_expr_gene[i,] ~ menp.f)$p.value)}
#OR
pvalmenp = vector()
for (i in 1:length(finalexpr[,1])) {
  pvalmenp = c(pvalmenp, t.test(finalexpr[i,] ~ menp.f2)$p.value)}
ttestpval = cbind(finalexpr[,0], pvalmenp)
write.table(ttestpval,"conservativepval.txt")

#Histogram of p-values
hist(pvalmenp, main= "P-values of T-test")

#QQ Plot
e = -log10(ppoints(length(pvalmenp)))
o = -log10(pvalmenp)
qqplot(e, o, main = "QQ Plot with final expression")
abline(0,1)

#Residuals
batch.residual = matrix(nrow= nrow(expr_quant.all), ncol = ncol(expr_quant.all))
rownames(batch.residual) = rownames(expr_quant.all)
colnames(batch.residual) = colnames(expr_quant.all)
for (i in 1:nrow(expr_quant.all)) {
  model= lm(expr_quant.all[i,]~ batch)
  batch.residual[i,] = resid(model)
}
cor.r <- cor(batch.residual,method="pearson", use="complete.obs")
heatmap.2(cor.r, main="Probe expression after batch control", key=T, revC=T, density.info="none", trace="none")

chip.residual = matrix(nrow= nrow(expr_quant.all), ncol = ncol(expr_quant.all))
rownames(chip.residual) = rownames(expr_quant.all)
colnames(chip.residual) = colnames(expr_quant.all)
for (i in 1:nrow(expr_quant.all)) {
  model= lm(expr_quant.all[i,]~ chip)
  chip.residual[i,] = resid(model)
}
cor.c <- cor(chip.residual,method="pearson", use="complete.obs")
heatmap.2(cor.c, main="Probe expression after chip control", key=T, revC=T, density.info="none", trace="none")

batch.chip.residual = matrix(nrow= nrow(expr_quant.all), ncol = ncol(expr_quant.all))
rownames(batch.chip.residual) = rownames(expr_quant.all)
colnames(batch.chip.residual) = colnames(expr_quant.all)
for (i in 1:nrow(expr_quant.all)) {
  model= lm(expr_quant.all[i,]~ batch + chip)
  batch.chip.residual[i,] = resid(model)
}
cor.bc <- cor(batch.chip.residual,method="pearson", use="complete.obs")
heatmap.2(cor.bc, main="Probe expression after batch + chip control", key=T, revC=T, density.info="none", trace="none")

b.c.r.bmi.residual = matrix(nrow= nrow(expr_quant.all), ncol = ncol(expr_quant.all))
rownames(b.c.r.bmi.residual) = rownames(expr_quant.all)
colnames(b.c.r.bmi.residual) = colnames(expr_quant.all)
for (i in 1:nrow(expr_quant.all)) {
  model= lm(expr_quant.all[i,]~ batch + chip + race + bmi)
  b.c.r.bmi.residual[i,] = resid(model)
}
cor.b.c.r.bmi <- cor(b.c.r.bmi.residual,method="pearson", use="complete.obs")
sum.PC.b.c.r.bmi <- prcomp(na.omit(cor.b.c.r.bmi))
summary(lm(sum.PC.b.c.r.bmi$rotation[, 1:7] ~ chip))

heatmap.2(cor.bc, main="Probe expression after batch + chip control", key=T, revC=T, density.info="none", trace="none")

#New PC on residuals
cor.final.resid <- cor(residual, method="pearson", use="complete.obs")
sum.PC.final.resid <- prcomp(na.omit(cor.final.resid))
sumsum.resid = summary(sum.PC.final.resid)
PCSf2 = sum.PC.final.resid$rotation
write.table(summary(sum.PC.final.resid)$importance,"summaryPCA2.txt")
summary(lm(sum.PC.final.resid$rotation[1:9, 1:7] ~ prot.f))
plot(c(1:9),sum.PC.final.resid$rotation[,1],col=colorfinal,xlab="Index of Samples",pch = 20, ylab=paste("PC 1 -",(sumsum.resid$importance[2,1]*100),"% of variance",sep=" "),main=title.PC)
text(c(1:9),sum.PC.final.resid$rotation[,1], fsampleorder[,5], cex = .5, pos=3)

#To get for p-values of covariates with PC on residuals
summary(lm(sum.PC.final.resid$rotation[1:9, 1:6] ~ prot.f))
summary(lm(sum.PC.final.resid$rotation[-2, 1:6] ~ bmi))

#T-test on residuals, histogram, heatmap and QQ plot
pvalr = vector()
for (i in 1:length(residual[,1])) {
  pvalr = c(pvalr, t.test(residual[i,] ~ menp.f2)$p.value)}
ttestresid = cbind(residual[,0], pvalr)
write.table(ttestresid,"residpval.txt")
hist(pvalr, main= "P-values of T-test on Residuals")
er = -log10(ppoints(length(pvalr)))
or = -log10(pvalr)
qqplot(er, or, main = "QQ Plot", xlab = "-log10 expected", ylab = "-log10 observed")
abline(0,1)
heatmap.2(cor.final.resid, main= "Correlation", key=T, revC=T,ColSideColors=colorfinal, density.info="none", trace="none")

#To pull out expression information
#For the expression of CDKL3 gene
finalexpr["CDKL3",]
finalstatus = c("pre", "post", "pre", "pre", "post", "post", "pre", "post", "post")
#To make a boxplot
boxplot(finalexpr["CDKL3",]~finalstatus, main = "Expression of CDKL3", ylab = "expression")

#GO Analysis
cutoff = .01
hgCutoff 0.001
params= new("GOHyperGParams",geneIds = ttestresid, universeGeneIds = symbols,annotation = "hgu95av2", ontology ="BP", pvalueCutoff = .01 , conditional = FALSE, testDirection = "over")
hgOver = hyperGTest(params)
conditional(params) = TRUE
hgCondOver = hyperGTest(params)

#Fixed Effects
menp = c(1,2,1,1,1,2,2,1,2,2)
menp.f = as.factor(menp)
cell= c(2,2,2,2,2,1,1,1,2,2)
cell.f= as.factor(cell)
loc = c(1,3,1,3,1,2,2,2,1,1)
loc.f= as.factor(loc)
prot = c(2,2,2,2,2,1,1,1,1,1)
prot.f= as.factor(prot)
date= c(4,1,2,3,1,7,6,7,5,7)
date.f = as.factor(date)
age = c(51, 53,53,45,49,45,53,49,53,50)
bmi = c(30.1,NA, 24.6,24.4,27,21.3,28.3,30.4,28.5,32.4)
## Later the fourth sample was deleted and these were also fixed

##Log likelihood test and then chi-squre analysis - didn't work
full.lm.fn <- function(x){
  ckagan<- try(lm( x ~ menp.f + age + cell.f + loc.f + prot.f + date.f + bmi))
  out <- as.numeric(logLik(ckagan))
  try(out)
}

system.time(logLik.full<- apply(averages_expr_gene,1,full.lm.fn))
logLik.full <- as.matrix(logLik.full)
write.table(logLik.full,"full_averagesoftechnicalreplicates.txt")

# Intercepts:
#bmi
int.lm.fn1 <- function(x){
  ckagan<- try(lm( x ~ menp.f + age + cell.f + loc.f + prot.f + date.f))
  out <- as.numeric(logLik(ckagan))
  try(out)
}
system.time(logLik.int1<- apply(averages_expr_gene,1,int.lm.fn1))
logLik.int1 <- as.matrix(logLik.int1)
write.table(logLik.int1,"interceptbmi_averagesoftechnicalreplicates.txt")

#cell count
int.lm.fn2 <- function(x){
  ckagan<- try(lm( x ~ menp.f + age + bmi + loc.f + prot.f + date.f))
  out <- as.numeric(logLik(ckagan))
  try(out)
}
system.time(logLik.int2<- apply(averages_expr_gene,1,int.lm.fn2))
logLik.int2 <- as.matrix(logLik.int2)
write.table(logLik.int2,"interceptcell_averagesoftechnicalreplicates.txt")

#location of cell culturing
int.lm.fn3 <- function(x){
  ckagan<- try(lm( x ~ menp.f + age + bmi + cell.f + prot.f + date.f))
  out <- as.numeric(logLik(ckagan))
  try(out)
}
system.time(logLik.int3<- apply(averages_expr_gene,1,int.lm.fn3))
logLik.int3 <- as.matrix(logLik.int3)
write.table(logLik.int3,"interceptloc_averagesoftechnicalreplicates.txt")

#protocol of cell culturing
int.lm.fn4 <- function(x){
  ckagan<- try(lm( x ~ menp.f + age + bmi + cell.f + loc.f + date.f))
  out <- as.numeric(logLik(ckagan))
  try(out)
}
system.time(logLik.int4<- apply(averages_expr_gene,1,int.lm.fn4))
logLik.int4 <- as.matrix(logLik.int4)
write.table(logLik.int4,"interceptprot_averagesoftechnicalreplicates.txt")

#age
int.lm.fn5 <- function(x){
  ckagan<- try(lm( x ~ menp.f + prot.f + bmi + cell.f + loc.f + date.f))
  out <- as.numeric(logLik(ckagan))
  try(out)
}
system.time(logLik.int5<- apply(averages_expr_gene,1,int.lm.fn5))
logLik.int5 <- as.matrix(logLik.int5)
write.table(logLik.int5,"interceptage_averagesoftechnicalreplicates.txt")

#date
int.lm.fn6 <- function(x){
  ckagan<- try(lm( x ~ menp.f + prot.f + bmi + cell.f + loc.f + age))
  out <- as.numeric(logLik(ckagan))
  try(out)
}
system.time(logLik.int6<- apply(averages_expr_gene,1,int.lm.fn6))
logLik.int6 <- as.matrix(logLik.int6)
write.table(logLik.int6,"interceptdate_averagesoftechnicalreplicates.txt")

#chisquare and p-value calculation 
full<-read.table("full_averagesoftechnicalreplicates.txt")
bmilm = read.table("interceptbmi_averagesoftechnicalreplicates.txt")
celllm = read.table("interceptcell_averagesoftechnicalreplicates.txt")
loclm = read.table("interceptloc_averagesoftechnicalreplicates.txt")
protlm = read.table("interceptprot_averagesoftechnicalreplicates.txt")
datelm = read.table("interceptdate_averagesoftechnicalreplicates.txt")
agelm = read.table("interceptage_averagesoftechnicalreplicates.txt")

chisquarebmi<-(-2*((bmilm[,1])-full[,1]))
chisquarecell<-(-2*((celllm[,1])-full[,1]))
chisquareloc<-(-2*((loclm[,1])-full[,1]))
chisquareprot<-(-2*((protlm[,1])-full[,1]))
chisquaredate<-(-2*((datelm[,1])-full[,1]))
chisquareage<-(-2*((agelm[,1])-full[,1]))

pvaluesbmi<-pchisq(chisquarebmi,df=1,lower.tail=FALSE)
pvaluescell<-pchisq(chisquarecell,df=1,lower.tail=FALSE)
pvaluesloc<-pchisq(chisquareloc,df=1,lower.tail=FALSE)
pvaluesprot<-pchisq(chisquareprot,df=1,lower.tail=FALSE)
pvaluesdate<-pchisq(chisquaredate,df=1,lower.tail=FALSE)
pvaluesage<-pchisq(chisquareage,df=1,lower.tail=FALSE)

#FDR corrected
write.table((p.adjust(pvaluesbmi, method="fdr")),"pvaluesbmi_fdr adjusted_averages.txt")
write.table((p.adjust(pvaluescell, method="fdr")),"pvaluescell_fdr adjusted_averages.txt")
write.table((p.adjust(pvaluesloc, method="fdr")),"pvaluesloc_fdr adjusted_averages.txt")
write.table((p.adjust(pvaluesprot, method="fdr")),"pvaluesprot_fdr adjusted_averages.txt")
write.table((p.adjust(pvaluesdate, method="fdr")),"pvaluesdate_fdr adjusted_averages.txt")
write.table((p.adjust(pvaluesage, method="fdr")),"pvaluesage_fdr adjusted_averages.txt")

p.adjbmi<-p.adjust(pvaluesbmi, method="fdr")
p.adjcell<-p.adjust(pvaluescell, method="fdr")
p.adjloc<-p.adjust(pvaluesloc, method="fdr")
p.adjprot<-p.adjust(pvaluesprot, method="fdr")
p.adjdate<-p.adjust(pvaluesdate, method="fdr")
p.adjage<-p.adjust(pvaluesage, method="fdr")

sum(p.adjbmi<=0.01)
sum(p.adjcell<=0.01)
sum(p.adjloc<=0.01)
sum(p.adjprot<=0.01)
sum(p.adjdate<=0.01)
sum(p.adjage<=0.01)

#Bonferroni corrected
bonferroni<-0.05/11392
sum(pvaluesbmi<=bonferroni)
sum(pvaluescell<=bonferroni)
sum(pvaluesloc<=bonferroni)
sum(pvaluesprot<=bonferroni)
sum(pvaluesdate<=bonferroni)
sum(pvaluesage<=bonferroni)

#calculate fold-change for differentially expressed genes
indbmi<-which(p.adjbmi<=0.01)
indcell<-which(p.adjcell<=0.01)
indloc<-which(p.adjloc<=0.01)
indprot<-which(p.adjprot<=0.01)
inddate<-which(p.adjdate<=0.01)
indage<-which(p.adjage<=0.01)

dexpgenesbmi<-symbolsUniq [indbmi]
dexpgenescell<-symbolsUniq [indcell]
dexpgenesloc<-symbolsUniq [indloc]
dexpgenesprot<-symbolsUniq [indprot]
dexpgenesdate<-symbolsUniq [inddate]
dexpgenesage<-symbolsUniq [indage]

#Final model
full.lm.fn2 <- function(x){
  ckagan<- try(lm( x ~ menp.f + age +loc.f + date.f + bmi))
  out <- as.numeric(logLik(ckagan))
  try(out)
}

system.time(logLik.full2<- apply(averages_expr_gene,1,full.lm.fn2))
logLik.full2 <- as.matrix(logLik.full2)
write.table(logLik.full2,"finalfull_averagesoftechnicalreplicates.txt")

# Menopause Intercept
int.lm.fn7 <- function(x){
  ckagan<- try(lm( x ~ age + loc.f +date.f + bmi))
  out <- as.numeric(logLik(ckagan))
  try(out)
}
system.time(logLik.int7<- apply(averages_expr_gene,1,int.lm.fn7))
logLik.int7 <- as.matrix(logLik.int7)
write.table(logLik.int7,"interceptmenp_averagesoftechnicalreplicates.txt")

full2<-read.table("finalfull_averagesoftechnicalreplicates.txt")
menp<-read.table("interceptmenp_averagesoftechnicalreplicates.txt")
chisquaremenp<-(-2*((menp[,1])-full2[,1]))
pvaluesmenp<-pchisq(chisquaremenp,df=1,lower.tail=FALSE)
p.adjmenp<-p.adjust(pvaluesmenp, method="fdr")
write.table(p.adjmenp,"pvaluesmenp_fdr adjusted_averages.txt")
sum(p.adjmenp<=0.01)
sum(pvaluesmenp<=bonferroni)
ind<-which(p.adjmenp<=0.01)
dexpgenesmenp<-symbolsUniq [ind]

####Re-do without the fourth indiv
pre<-c(1,3,4,5,8)
post<-c(2,6,7,9,10)
pre.data<-averages_expr_gene[,pre]
post.data<-averages_expr_gene[,post]
meanpre<-apply(pre.data,1,mean)
meanpost<-apply(post.data,1,mean)
foldchange<-meanpre/meanpost
foldchange<-as.matrix(foldchange)
foldchange2<- rownames(foldchange)
foldformenp<-match(dexpgenesmenp,foldchange2)
foldchange3<-foldchange[foldformenp,]
foldchange4<-which(foldchange3>=1.5)

#Volcano plot 
x.axis<-log2(foldchange)
y.axis<-(-log10(pvaluesmenp))
logfold<-log2(foldchange)
plot(x.axis,y.axis,pch=20,cex=0.5,xlim=c(-1,1),xlab="log2(foldchange)",ylab="-log10(pvalue)",main="Genes diff expressed between pre- and post- menopausal women")
threshold<-log2(1.5)
threshold2<-(-threshold)
threshold3<--log10(0.05/10641)
abline(v=threshold,col="blue",pch=2)
abline(v=threshold2,col="blue",pch=2)
abline(h=threshold3,col="blue",pch=2)

beta <- function(x){
  ckagan<- lm( x ~ menp.f)}
output = apply(averages_expr_gene,1,beta)
summary(output)
plot(expr_gene_o["POLR3F",], col=color)

##Example of plotting expression
rv = c(2,4,6,8,10,12,14,16,18,20,22,24)
unstim = c(1,3,5,7,9,11,13,15,17,19,21,23)
cases = c(5,6,7,8,9,10,11,12,21,22,23,24)
control = c(1,2,3,4,13,14,15,16,17,18,19,20)
boxplot(expr_gene["TNFSF13B",control],expr_gene["TNFSF13B",cases], main = "Expression of BAFF by Case/Control", ylab = "Normalized Expression", xlab = "Control                                                         Cases")
boxplot(expr_gene["TNFSF13B",unstim],expr_gene["TNFSF13B",rv], main = "Expression of BAFF by Treatment Group", ylab = "Normalized Expression", xlab = "Unstimulated                                                       RV Stimulated")

