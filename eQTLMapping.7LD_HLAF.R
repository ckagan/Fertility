setwd("C:/Users/Courtney/Dropbox/Impute eqtl")

#Pull all SNPs in .7 or higher LD with rs 2523393 in excel
seventy = read.table('rs2523393.7LD.txt', header=T, as.is=T)
map = read.table('imputedEVEmap.txt', header=T, as.is=T)
seventyloc = merge(seventy, map, by.x = "BP_B", by.y = "pos", all.x=T, all.y=F)

map_geno = read.table('C:/Users/Courtney/Dropbox/Ober Lab/Fertility/Final Paper eQTL Analysis/eQTL Mapping/QCgeno_all57_map_fullSNPname.txt', header=T, as.is=T)
seventyloc = merge(seventyloc, map_geno, by.x = "BP_B", by.y = "pos", all.x=T, all.y=F)

seventyloc_geno = seventyloc[rowSums(is.na(seventyloc)) < 4,]

#Pull out SNPs in original genotyping
genoSNPs = seventyloc_geno$SNP.y
genoSNPs = na.omit(genoSNPs)

#Pull out imputed SNPs that are not in the original genotyping
impSNPs = seventyloc_geno[which(rowSums(is.na(seventyloc_geno)) >= 3), 3]
imp = as.matrix(impSNPs)
colnames(imp)= c("Marker")
impmap = merge(imp, map, by.x = "Marker", by.y = "SNP", all.x=T, all.y=F)
colnames(impmap) = c("SNP", "chr_snp", "pos")

#Filter imputed genotypes
impgeno = read.table('imputedEVEdose.txt', header=T, as.is=T)
impdose = merge(impgeno, imp, by.x = "row.names", by.y = "Marker", all.x=F, all.y=T)
impdose = na.omit(impdose)
row.names(impdose) <- NULL 

#Filter genotypes
QCgeno = read.table('C:/Users/Courtney/Dropbox/Ober Lab/Fertility/Final Paper eQTL Analysis/eQTL Mapping/QCgeno_noFC_53.txt', header=T, as.is=T)
geno = as.matrix(genoSNPs)
QCgenofilt = merge(QCgeno, geno, by.x = "row.names", by.y = "V1", all.x=F, all.y=T)

genomap = merge(map_geno, geno, by.x = "SNP", by.y = "V1", all.x=F, all.y=T)

#The first row is a SNP with the same position on chr 17 - remove it from both geno and map
eqtlgeno = QCgenofilt[-1,]
eqtlmap = genomap[-1,-2]
colnames(eqtlmap) = c("SNP", "chr_snp", "pos")
row.names(eqtlmap) <- NULL 
row.names(eqtlgeno) <- NULL 

#Now append all genotypes and map
totalmap = rbind(impmap, eqtlmap)

#I'm not sure how to sort columns so sort columns and merge genotypes in excel
write.table(impdose, 'imputeddose_HLAF.7.txt', quote=F, sep='\t')
write.table(eqtlgeno, 'QCgeno_HLAF.7.txt', quote=F, sep='\t')

write.table(totalmap, 'QCgeno+impdose_HLAF.7_map.txt', quote=F, sep='\t')