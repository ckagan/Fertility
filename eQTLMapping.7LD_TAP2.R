setwd("C:/Users/Courtney/Dropbox/Impute eqtl")

#Pull all SNPs in .7 or higher LD with rs 2071473 in excel
seventy = read.table('rs2071473.7LD.txt', header=T, as.is=T)
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
write.table(impdose, 'imputeddose_TAP2.7.txt', quote=F, sep='\t')
write.table(eqtlgeno, 'QCgeno_TAP2.7.txt', quote=F, sep='\t')

write.table(totalmap, 'QCgeno+impdose_TAP2.7_map.txt', quote=F, sep='\t')

#####eQTL Mapping####
library("MatrixEQTL")
useModel = modelLINEAR
#Example All Genotypes
SNP_file_name = "QCgeno+impdose_TAP2.7.txt" 
expression_file_name = "HLARegion_FinalGeneExpr.txt"

#In case of no covariates set the variable covariates_file_name to character() in R ([] in Matlab).

covariates_file_name = character() 
output_file_name = "eQTL_results_TAP2_imputedandgeno_trans.txt"
output_file_name.cis = "eQTL_results_TAP2_imputedandgeno_cis.txt"

##Next, choose the p-value threshold. Only associations significant at this level will be saved in the output file. Note that for larger datasets the threshold should be lower. Setting the threshold to a high value for a large dataset may cause excessively large output (many gigabytes).

pvOutputThreshold = 0
cisDist = 200000
pvOutputThreshold.cis = 1

##Finally, define the covariance matrix for the error term. This a seldom used parameter. If the matrix is a multiple of identity, set the parameter to numeric() in R ([] in Matlab).

errorCovariance = numeric()

#The next section of the sample code contains three very similar parts loading the files with genotype, gene expression, and covariates. In each part one can set the file delimiter (i.e. tab "\t", comma ",", or space " "), the string indicating missing values, the number of rows of column labels, and the columns of row labels. Finally, one can change the number of the variables in a slice for the file reading procedure (do not change if not sure).

snps = SlicedData$new()
snps$fileDelimiter = "\t"
# the TAB character
snps$fileOmitCharacters = "NA" 
# denote missing values;
snps$fileSkipRows = 1 
# one row of column labels
snps$fileSkipColumns = 1 
# one column of row labels
snps$fileSliceSize = 8000 
# read file in pieces of 10,000 rows
snps$LoadFile(SNP_file_name)

gene = SlicedData$new()
gene$fileDelimiter = '\t'
# a space
gene$fileOmitCharacters = 'NA'
# denote missing values;
gene$fileSkipRows = 1
# one row of column labels
gene$fileSkipColumns = 1
# one column of row labels
gene$fileSliceSize = 10000
# read file in pieces of 10,000 rows
gene$LoadFile(expression_file_name)

cvrt = SlicedData$new()
cvrt$fileDelimiter = '\t'
# the TAB character
cvrt$fileOmitCharacters = 'NA'
# denote missing values;
cvrt$fileSkipRows = 1
# one row of column labels
cvrt$fileSkipColumns = 1
# one column of row labels
cvrt$fileSliceSize = snps$nCols()+1
# read file in one piece
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name)
}
snpspos = read.table('QCgeno+impdose_TAP2.7_map.txt', header = TRUE, stringsAsFactors = FALSE);
genepos = read.table('HLARegion_genelocmap.txt', header = TRUE, stringsAsFactors = FALSE);


#Finally, the main matrix eQTL function is called:
me = Matrix_eQTL_main(snps,
                      gene,
                      cvrt,
                      output_file_name,
                      pvOutputThreshold,
                      useModel,
                      errorCovariance,
                      verbose = TRUE,
                      output_file_name.cis,
                      pvOutputThreshold.cis,
                      snpspos,
                      genepos,
                      cisDist,
                      pvalue.hist = 100)

#####Merging with TTP results####
#First add location to results
setwd("C:/Users/Courtney/Dropbox/Impute eqtl")
eqtl = read.table('eQTL_results_TAP2_imputedandgeno_cis.txt', as.is=T, header=T)
totalmap = read.table('QCgeno+impdose_TAP2.7_map.txt', header=T, as.is=T)
colnames(totalmap) = c("Marker", "chr_snp", "pos") 
eqtl_wloc = merge(eqtl,totalmap, by.x = "SNP", by.y = "Marker", all.x=T, all.y = F)

#Merge with TTP results
setwd("C:/Users/Courtney/Dropbox/Ober Lab/Fertility/Final Paper eQTL Analysis/HutteriteTTP")
ttptap2 = read.table('TAP2_0.5MbRegion_forR.txt', as.is=T, header=T)
supt5 = merge(ttptap2, eqtl_wloc, by.x = "stopBP", by.y = "pos", all.x=T, all.y=T)

write.table(supt5, 'TAP2_TTP+eqtlwgeno+imp.txt', quote=F, sep='\t')
