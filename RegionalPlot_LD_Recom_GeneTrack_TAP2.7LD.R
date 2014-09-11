setwd("C:/Users/Courtney/Dropbox/Ober Lab/Fertility/Final Paper eQTL Analysis/Regional Plot")

tap2 = read.table('TAP2_0.5MbRegion_eQTL_forR.txt', as.is=TRUE, header=TRUE)
gscale = gray(((1-tap2[which(tap2$R2.with.rs2071473 > .7),9 ])/5), alpha=.3)
manhattan <- function(x, chr="CHR", bp="BP", p="P", snp="SNP", r2 = "R2.with.rs2071473", 
                       col=c("gray10", "gray60"), ymax=NULL, 
                       suggestiveline=-log10(1e-3), genomewideline=-log10(1e-3), 
                       highlight=NULL, logp=TRUE, ...) {
  
  # Not sure why, but package check will warn without this.
  P=index=NULL
  
  # Check for sensible dataset
  ## Make sure you have chr, bp and p columns.
  if (!(chr %in% names(x))) stop(paste("Column", chr, "not found!"))
  if (!(bp %in% names(x))) stop(paste("Column", bp, "not found!"))
  if (!(p %in% names(x))) stop(paste("Column", p, "not found!"))
  ## warn if you don't have a snp column
  if (!(snp %in% names(x))) warning(paste("No SNP column found. OK unless you're trying to highlight."))
  ## make sure chr, bp, and p columns are numeric.
  if (!is.numeric(x[[chr]])) stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
  if (!is.numeric(x[[bp]])) stop(paste(bp, "column should be numeric."))
  if (!is.numeric(x[[p]])) stop(paste(p, "column should be numeric."))
  
  # Create a new data.frame with columns called CHR, BP, and P.
  d=data.frame(CHR=x[[chr]], BP=x[[bp]], P=x[[p]], R2=x[[r2]])
  
  # If the input data frame has a SNP column, add it to the new data frame you're creating.
  if (!is.null(x[[snp]])) d=transform(d, SNP=x[[snp]])
  
  # Set positions, ticks, and labels for plotting
  ## Sort, keep only SNPs with p-values between 0 and 1
  d=subset(d[order(d$CHR, d$BP), ], (P>0 & P<=1 & is.numeric(P)))
  #d$logp <- ifelse(logp, yes=-log10(d$P), no=d$P)
  if (logp) {
    d$logp <- -log10(d$P)
  } else {
    d$logp <- d$P
  }
  d$pos=NA
  
  # Set y maximum. If ymax is undefined, not numeric, or negative, set it
  # equal to the most significant SNP.
  if (is.null(ymax)) { # still null
    ymax = ceiling(max(d$logp))
    message("Ymax will be set automatically based on most significant SNP")
  } else if (!is.numeric(ymax)){  # not numeric
    ymax = ceiling(max(d$logp))
    warning('non-numeric ymax argument.')
  } else if (ymax < 0) { # negative
    ymax = ceiling(max(d$logp))
    warning('negative ymax argument.')
  } else if (is.null(ymax)) { # still null
    ymax = ceiling(max(d$logp))
    message(paste("Using", ymax, "as ymax"))
  }
  message(paste("Using", ymax, "as ymax"))
  
  # Fixes the bug where one chromosome is missing by adding a sequential index column.
  d$index=NA
  ind = 0
  for (i in unique(d$CHR)){
    ind = ind + 1
    d[d$CHR==i,]$index = ind
  }
  
  # This section sets up positions and ticks. Ticks should be placed in the
  # middle of a chromosome. The a new pos column is added that keeps a running
  # sum of the positions of each successive chromsome. For example:
  # chr bp pos
  # 1   1  1
  # 1   2  2
  # 2   1  3
  # 2   2  4
  # 3   1  5
  nchr = length(unique(d$CHR))
  if (nchr==1) { ## For a single chromosome
    d$pos=d$BP
    ticks=floor(length(d$pos))/2+1
    xlabel = paste('Chromosome',unique(d$CHR),'position')
    labs = ticks
  } else { ## For multiple chromosomes
    lastbase=0
    ticks=NULL
    for (i in unique(d$index)) {
      if (i==1) {
        d[d$index==i, ]$pos=d[d$index==i, ]$BP
      } else {
        lastbase=lastbase+tail(subset(d,index==i-1)$BP, 1)
        d[d$index==i, ]$pos=d[d$index==i, ]$BP+lastbase
      }
      ticks=c(ticks, d[d$index==i, ]$pos[floor(length(d[d$index==i, ]$pos)/2)+1])
    }
    xlabel = 'Chromosome'
    #labs = append(unique(d$CHR),'') ## I forgot what this was here for... if seems to work, remove.
    labs=unique(d$CHR)
  }
  
  # Initialize plot
  xmax = 32807000
  xmin = 32740000
  ymin = -4.2
  plot(NULL, xaxt='n', bty='n', xaxs='i', yaxs='i', xlim=c(xmin,xmax), ylim=c(ymin,ymax),
       xlab=xlabel, ylab=expression(-log[10](italic(p))), las=1, pch=20, ...)
  
  # Add an axis. 
  if (nchr==1) { #If single chromosome, ticks and labels automatic.
    axis(1, ...)
  } else { # if multiple chrs, use the ticks and labels you created above.
    axis(1, at=ticks, labels=labs, ...)
  }
  
  
  # Add suggestive and genomewide lines
  #if (suggestiveline) abline(h=suggestiveline, col="black")
  if (genomewideline) {
    d.highlight.snp=d[which(d$SNP == 'rs2071473_r_T'), ]
    d.highlight.snp1=d[which(d$R2>.7), ]
    d.highlight.snp2=d[which(d$R2<.7), ]
    with(d.highlight.snp2, points(pos, logp, col="skyblue", pch=16,cex=1, ...))
    with(d.highlight.snp, points(pos, logp, col="red", pch=18,cex=4, ...)) 
    with(d.highlight.snp1, points(pos, logp, col=gscale, pch=16,cex=1.5, ...))    
    
  }
  
  # Highlight snps from a character vector
  if (!is.null(highlight)) {
    if (any(!(highlight %in% d$SNP))) warning("You're trying to highlight SNPs that don't exist in your results.")
    d.highlight=d[which(d$SNP %in% highlight), ]
    d.highlight1 = d.highlight[which(-log10(d.highlight$P) < genomewideline),]
    d.highlight2 = d.highlight[which(-log10(d.highlight$P) > genomewideline),]
    with(d.highlight1, points(pos, logp, col="forestgreen", pch=20, ...)) 
    with(d.highlight2, points(pos, logp, col="blue", pch=20, ...)) 
  }
  
}
manhattan(tap2)
genelist <- read.table(paste("known_genes_buildhg19_chr6_strand.txt"), header=T)
xmax = 32807000
xmin = 32740000
genes <- subset(genelist, ( genelist$STOP > xmin & genelist$START < xmax ) )

##Remove duplicate genes (picked longest)
rem=c(2,4)
genes.in.locus = genes[-rem,]

range = 5
offset <- ( range * 4 / 3 ) - range
big.range <- range + offset 

for ( i in 1:2 ) { 
  if ( genes.in.locus[i,]$STRAND == "+" ) {
    arrows(max(genes.in.locus[i,]$START), -offset, min(genes.in.locus[i,]$STOP), -offset, length=0.10, lwd=2, code=2, lty="solid", col="darkgreen")
  } else {  	
    arrows(max(genes.in.locus[i,]$START), -offset, min(genes.in.locus[i,]$STOP), -offset, length=0.10, lwd=2, code=1, lty="solid", col="darkgreen")
  }
  if ( ! is.na(genes.in.locus[i,]$GENE) ) {
    text(genes.in.locus[i,]$START + (genes.in.locus[i,]$SIZE / 2), -offset + ( big.range / 20 ), labels=genes.in.locus[i,]$GENE, cex=0.8)
  }
}



abline(h=0, col="black", lty=2)
abline(h=-4.2, col="black", lty=1)
abline(v=xmax, col="black")
abline(v=xmin, col="black")
chr=6
max.pos = 32807000
min.pos = 32740000
recomb <- read.table(paste("genetic_map_chr", chr, ".txt", sep=""), header=T)
keep.recomb <- subset(recomb, recomb[,1] > min.pos & recomb[,1] < max.pos)
ystart.recomb <- - offset + (big.range / 8)
lines(keep.recomb[,1], ystart.recomb + ( ( keep.recomb[,2] / 60 ) * ( 6 * big.range / 8 )), type="l", col="black", lwd=1, xlim=c(min.pos, max.pos), ylim=c(-offset,range), xlab="", ylab="",  axes=F)
axis(4, at=c( ystart.recomb, ystart.recomb + (big.range / 4), ystart.recomb + ( 2 * big.range / 4)), labels=c("0","20","40"), las=1)
mtext("Recombination rate (cM/Mb)", side=4, at=(-offset+big.range/2), line=2)


#Add in functional tracks 
setwd("C:/Users/Courtney/Dropbox/Ober Lab/Fertility/Final Paper eQTL Analysis/HutteriteTTP")
fairelist <- read.table("FAIRE_region.txt", header=T)
faire <- subset(fairelist, ( fairelist$Start > xmin & fairelist$Start < xmax ))
for ( i in 1:dim(faire)[1] ) { 
  segments(faire[i,]$Start, -2.2, faire[i,]$Stop, -2.2, col="red", lty = 1, lwd =5)
}
#text(32747000, -2.2, labels="FAIRE", cex= .8, col = "red")

dnalist <- read.table("DNase1_region.txt", header=T)
dna <- subset(dnalist, ( dnalist$Start > xmin & dnalist$Start < xmax ))
for ( i in 1:dim(dna)[1] ) { 
  segments(dna[i,]$Start, -2.7, dna[i,]$Stop, -2.7, col="blue", lty = 1, lwd =5)
}
#text(32747000, -2.7, labels="DNaseI", cex= .8, col = "blue")

hklist <- read.table("H3K4me3_region.txt", header=T)
hk <- subset(hklist, ( hklist$Start > xmin & hklist$Start < xmax ))
for ( i in 1:dim(hk)[1] ) { 
  segments(hk[i,]$Start, -3.2, hk[i,]$Stop, -3.2, col="black", lty = 1, lwd =5)
}
#text(32747000, -3.2, labels="H3K4me3", cex= .8, col = "black")


nr2list <- read.table("NR2F2_region.txt", header=T)
nr2 <- subset(nr2list, ( nr2list$Start > xmin & nr2list$Start < xmax ))
for ( i in 1:dim(nr2)[1] ) { 
  segments(nr2[i,]$Start, -3.7, nr2[i,]$Stop, -3.7, col="purple", lty = 1, lwd =5)
}
#text(32747000, -3.7, labels="NR2F2", cex= .8, col = "purple")

##To make key
plot(1:10,0:-9)
text(2, -2.2, labels="FAIRE", cex= .8, col = "red")
text(2, -2.7, labels="DNaseI", cex= .8, col = "blue")
text(2, -3.2, labels="H3K4me3", cex= .8, col = "black")
text(2, -3.7, labels="NR2F2", cex= .8, col = "purple")


#Make a legend for the LD
layout(matrix(1:2,ncol=2), width = c(2,1),height = c(1,1))
plot(1:2, 1:2, pch = 19, cex=2, col = c("#0F0F0F4D" ,"black") )

legend_image <- as.raster(matrix(c("black","#0F0F0F4D" ), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'legend title')
text(x=1.5, y = seq(.5,1,l=5), labels = seq(.7,1,l=5))
rasterImage(legend_image, 0, .5, 1,1)

legend_image2 <- as.raster(matrix(c("skyblue" ), ncol=1))
text(x=1.5, y = seq(-.3,0,l=1), labels = seq(0,l=1))
text(x=1.5, y = .35, labels = "<0.7")
rasterImage(legend_image2, 0, .3, 1,.5)

#Generate a list of SNPs in FAIRE-seq sites
setwd("C:/Users/Courtney/Dropbox/Ober Lab/Fertility/Final Paper eQTL Analysis/Regional Plot")
tap2 = read.table('TAP2_0.5MbRegion_eQTL_forR.txt', as.is=TRUE, header=TRUE)
sub = tap2[which(tap2$R2.with.rs2071473 > .7),]

xmax = 32807000
xmin = 32740000
setwd("C:/Users/Courtney/Dropbox/Ober Lab/Fertility/Final Paper eQTL Analysis/HutteriteTTP")
fairelist <- read.table("FAIRE_region.txt", header=T)
faire <- subset(fairelist, ( fairelist$Start > xmin & fairelist$Start < xmax ))

overlapfaire =matrix(NA, ncol=ncol(sub))
colnames(overlapfaire) = colnames(sub)
#Created a 10bp window around FAIRE-seq site
for ( i in 1:dim(faire)[1] ) { 
  temp = sub[which(sub$startBP >= (faire[i,]$Start) & sub$startBP <= (faire[i,]$Stop)),]
  overlapfaire = rbind(overlapfaire, temp)
}

