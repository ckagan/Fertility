setwd("C:/Users/Courtney/Dropbox/Ober Lab/Fertility/Final Paper eQTL Analysis/eQTL Mapping")

library("MatrixEQTL")
useModel = modelLINEAR
#Example All Genotypes
SNP_file_name = "QCgeno_noFC_53.txt" 
expression_file_name = "Final_GeneExprDetected.7_woFC.txt"

#In case of no covariates set the variable covariates_file_name to character() in R ([] in Matlab).

covariates_file_name = character() 
output_file_name = "eQTL_results_trans_allpval_QCgeno.txt"
output_file_name.cis = "eQTL_results_cis_allpval_QCgeno_trans.txt"

##Next, choose the p-value threshold. Only associations significant at this level will be saved in the output file. Note that for larger datasets the threshold should be lower. Setting the threshold to a high value for a large dataset may cause excessively large output (many gigabytes).

pvOutputThreshold = .000001

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
snpspos = read.table('QCgeno_all57_map_fullSNPname.txt', header = TRUE, stringsAsFactors = FALSE)[,2:4];
genepos = read.table('minal_gene_location.txt', header = TRUE, stringsAsFactors = FALSE);

cisDist = 200000
#Finally, the main matrix eQTL function is called:
pvOutputThreshold.cis = 0
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
