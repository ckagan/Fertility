#PCA at before avg replicates
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
covars<-list(chip.f,season.f,hypoth.f, race.f, bmi, age, prmisc.num, class.f, batch.f)

samplekey = read.csv('AllColumnNames_forsup.csv')
rem = c(99,  107,  28,	75,	37,	68,	48,	79,	40,	51,	77,	106,	72,	108,	11,	21,	39,	63,	83,	112,	27,	54,	60,	76)
clean.key = samplekey[-rem,]

#x.pca = prcomp(na.omit(expr_quant.all), scale = T, center = T)
x.pca = prcomp(na.omit(bc.resid.int), scale = T, center = T)
x.pca.sum = summary(x.pca)
plot(x.pca$rotation[-rem,1], x.pca$rotation[-rem,2], xlab = paste('PC1 (',x.pca.sum$importance[2,1], ')', sep = ''),
     ylab = paste('PC2 (',x.pca.sum$importance[2,2], ')', sep = ''), main = "PC1/2", pch = 20, col = clean.key$Batch) 
legend(x = "topleft", pch = 20, col = c(1,2), c("Batch 1", "Batch 2"))

plot(x.pca$rotation[-rem,1], x.pca$rotation[-rem,3], xlab = paste('PC1 (',x.pca.sum$importance[2,1], ')', sep = ''),
     ylab = paste('PC3 (',x.pca.sum$importance[2,3], ')', sep = ''), main = "PC1/3", pch = 20, col = clean.key$Batch)
legend(x = "topleft", pch = 20, col = c(1,2), c("Batch 1", "Batch 2"))

plot(x.pca$rotation[-rem,1], x.pca$rotation[-rem,4], xlab = paste('PC1 (',x.pca.sum$importance[2,1], ')', sep = ''),
     ylab = paste('PC4 (',x.pca.sum$importance[2,4], ')', sep = ''), main = "PC1/4", pch = 20, col = clean.key$Batch)
legend(x = "topleft", pch = 20, col = c(1,2), c("Batch 1", "Batch 2"))



lmPCA = function(pca, covars, npcs)
{
  results<-c()
  for (f in covars) {
    for (i in 1:npcs)
    {
      s = summary(lm(pca$rotation[-rem,i]~f[-rem]));
      results<-c(results,pf(s$fstatistic[[1]],
                            s$fstatistic[[2]],s$fstatistic[[3]], lower.tail = FALSE),
                 s$adj.r.squared)
    }
  }
  resultsM<-matrix(nrow = length(covars), ncol = 2*npcs, data =
                     results, byrow = TRUE)
  resultsM
  
  
}

pcaresults = lmPCA(x.pca,covars,4)
rownames(pcaresults) = c("chip","season","hypoth", "race", "bmi", "age", "prmisc", "class", "batch")
colnames(pcaresults) = c("PC1 pval", "PC1 adj R sqs","PC2 pval", "PC2 adj R sqs","PC3 pval", "PC3 adj R sqs","PC4 pval", "PC4 adj R sqs")
write.table(pcaresults, 'PCAresults_beforereg.txt', sep='\t', quote=F)
write.table(pcaresults.reg, 'PCAresults_afterreg.txt', sep='\t', quote=F)


#PCA after averaging replicates
final_probe_noreg= cbind(expr_quant.all[, c(1:2,4:10,13:15,18:20,22,31,93,46,58,62)],m306, 
                         m436, m334, m336, m274,mD0143,m165,m329,m226,m309,m451,m297,
                         m235,m370,m348,m347,m128,m325,m441, m216,m363,m305,m387,m257,
                         m156,m321,m356,m385,m423,m391,m395,m408)

list = colnames(final_probe_noreg)
list.r = sub("m", "s", list)
list.m = as.matrix(list.r)
colnames(list.m) = c("Name")
sample.key.53 = read.csv('AllColumnNames_final53.csv')
avg_sample.key.53 = merge(list.m, sample.key.53, by.x = "Name", by.y = "RLP.Only", sort=F)


#Do PCA on averaged samples
x.pca = prcomp(na.omit(final_probe_noreg), scale = T, center = T)
x.pca.sum = summary(x.pca)
plot(x.pca$rotation[,1], x.pca$rotation[,2], xlab = paste('PC1 (',x.pca.sum$importance[2,1], ')', sep = ''),
     ylab = paste('PC2 (',x.pca.sum$importance[2,2], ')', sep = ''), main = "PC1/2", pch = 20, col = avg_sample.key.53$Batch) 
legend(x = "topright", pch = 20, col = c(1,2), c("Batch 1", "Batch 2"))
plot(x.pca$rotation[,2], x.pca$rotation[,3], xlab = paste('PC2 (',x.pca.sum$importance[2,2], ')', sep = ''),
     ylab = paste('PC3 (',x.pca.sum$importance[2,3], ')', sep = ''), main = "PC2/3", pch = 20, col = avg_sample.key.53$Chip)


plot(x.pca$rotation[,3], x.pca$rotation[,4], xlab = paste('PC3 (',x.pca.sum$importance[2,3], ')', sep = ''),
     ylab = paste('PC4 (',x.pca.sum$importance[2,4], ')', sep = ''), main = "PC3/4", pch = 20, col = avg_sample.key.53$race)

batch = avg_sample.key.53$Batch
chip = avg_sample.key.53$Chip
age =avg_sample.key.53$Age
bmi = avg_sample.key.53$BMI
prmisc = avg_sample.key.53$prmisc
hypoth = avg_sample.key.53$hypothyroid
race = avg_sample.key.53$race
class = avg_sample.key.53$class
season = avg_sample.key.53$season
#Converted categorical covariates to a factor so they are levels 
batch.f = as.factor(batch)
chip.f = as.factor(chip)
hypoth.f = as.factor(hypoth)
race.f = as.factor(race)
season.f = as.factor(season)
class.f = as.factor(class)
#Converted numerical covariates to a numeric so they are continuous
prmisc.num = as.numeric(prmisc)
covars<-list(chip.f,season.f,hypoth.f, race.f, bmi, age, prmisc.num, class.f, batch.f)

lmPCA = function(pca, covars, npcs)
{
  results<-c()
  for (f in covars) {
    for (i in 1:npcs)
    {
      s = summary(lm(pca$rotation[,i]~f));
      results<-c(results,pf(s$fstatistic[[1]],
                            s$fstatistic[[2]],s$fstatistic[[3]], lower.tail = FALSE),
                 s$adj.r.squared)
    }
  }
  resultsM<-matrix(nrow = length(covars), ncol = 2*npcs, data =
                     results, byrow = TRUE)
  resultsM
  
  
}

pcaresults = lmPCA(x.pca,covars,4)
rownames(pcaresults) = c("chip","season","hypoth", "race", "bmi", "age", "prmisc", "class", "batch")
colnames(pcaresults) = c("PC1 pval", "PC1 adj R sqs","PC2 pval", "PC2 adj R sqs","PC3 pval", "PC3 adj R sqs","PC4 pval", "PC4 adj R sqs")


#After regressing out batch and chip
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

final_probe_reg= cbind(bc.resid.int[, c(1:2,4:10,13:15,18:20,22,31,93,46,58,62)],m306, 
                         m436, m334, m336, m274,mD0143,m165,m329,m226,m309,m451,m297,
                         m235,m370,m348,m347,m128,m325,m441, m216,m363,m305,m387,m257,
                         m156,m321,m356,m385,m423,m391,m395,m408)
colnames(bc.resid.int) = samplekey$SampleID
plot(hclust(dist(t(bc.resid.int[,-rem]))), xlab = "hclust/Euclidean distance", main = "Dendrogram")
plot(hclust(as.dist(1-cor(as.matrix(bc.resid.int[,-rem])))),xlab = "Pearson", main = "Dendrogram using Pearson Correlation")
colnames(bc.resid.int) = samplekey$RLP.Only

cor.bc.clean <- cor(bc.resid.int[,-rem],method="spearman", use="complete.obs")
replicates=c()
non.replicates=c()
samplekey.clean = samplekey[-rem,]
names = as.character(samplekey.clean$RLP.Only)
for(i in 1:87){
  for(j in (i+1):88){
    if(names[i]==names[j]){
      replicates=c(replicates,cor.bc.clean[i,j])
    
    } else {
      non.replicates=c(non.replicates,cor.bc.clean[i,j])
    }
  }
}
boxplot(replicates,non.replicates, main = "Correlation of Samples", ylab = "Correlation", xlab = "Replicates                                              Non-Replicates")

summary(replicates)
summary(non.replicates)


cor.expr.clean <- cor(expr_quant.all[,-rem],method="spearman", use="complete.obs")
replicates=c()
non.replicates=c()
samplekey.clean = samplekey[-rem,]
names = as.character(samplekey.clean$RLP.Only)
for(i in 1:87){
  for(j in (i+1):88){
    if(names[i]==names[j]){
      replicates=c(replicates,cor.expr.clean[i,j])
      
    } else {
      non.replicates=c(non.replicates,cor.expr.clean[i,j])
    }
  }
}
boxplot(replicates,non.replicates, main = "Correlation of Samples", ylab = "Correlation", xlab = "Replicates                                              Non-Replicates")

summary(replicates)
summary(non.replicates)