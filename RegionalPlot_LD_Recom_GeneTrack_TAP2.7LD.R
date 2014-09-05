tap2 = read.table('TAP2_0.5MbRegion_eQTL_forR_.70LD.txt', as.is=TRUE, header=TRUE)
manhattan4 <- function(x, chr="CHR", bp="BP", p="P", snp="SNP", r2 = "R2.with.rs2071473",fdr = "FDR", 
                       col=c("gray10", "gray60"), ymax=NULL, 
                       suggestiveline=-log10(1e-3), genomewideline=-log10(5e-8), 
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
  d=data.frame(CHR=x[[chr]], BP=x[[bp]], P=x[[p]], R2=x[[r2]], FDR=x[[fdr]])
  
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
  xmax = 33032600
  xmin = 32534300
  ymin = -2
  plot(NULL, xaxt='n', bty='n', xaxs='i', yaxs='i', xlim=c(xmin,xmax), ylim=c(ymin,ymax),
       xlab=xlabel, ylab=expression(-log[10](italic(p))), las=1, pch=20, ...)
  
  # Add an axis. 
  if (nchr==1) { #If single chromosome, ticks and labels automatic.
    axis(1, ...)
  } else { # if multiple chrs, use the ticks and labels you created above.
    axis(1, at=ticks, labels=labs, ...)
  }
  
  # Create a vector of alternatiting colors
  col=rep(col, max(d$CHR))
  
  # Add points to the plot
  if (nchr==1) {
    with(d, points(pos, logp, pch=20, ...))
  } else {
    # if multiple chromosomes, need to alternate colors and increase the color index (icol) each chr.
    icol=1
    for (i in unique(d$index)) {
      with(d[d$index==unique(d$index)[i], ], points(pos, logp, col=col[icol], pch=20, ...))
      icol=icol+1
    }
  }
  
  # Add suggestive and genomewide lines
  if (suggestiveline) abline(h=suggestiveline, col="black")
  if (genomewideline) {
    d.highlight2=d[which(d$R2<1 & d$R2>.8), ]
    d.highlight3=d[which(d$R2<.8 & d$R2>.6), ]
    d.highlight4=d[which(d$R2<.6 & d$R2>.4), ]
    d.highlight5=d[which(d$R2<.4 & d$R2>.2), ]
    d.highlight6=d[which(d$R2<.2), ]
    d.highlight1=d[which(d$R2==1), ]
    d.highlight = d[which(d$SNP == "rs2071473_r_T"),]
    with(d.highlight6, points(pos, logp, col="#F2F2F2", pch=20, ...))
    with(d.highlight5, points(pos, logp, col="#C1C1C1", pch=20, ...))
    with(d.highlight4, points(pos, logp, col="#919191", pch=20, ...))
    with(d.highlight3, points(pos, logp, col="#606060", pch=20, cex = 1.25,...))
    with(d.highlight2, points(pos, logp, col="#303030", pch=20, cex=1.5,...))
    with(d.highlight1, points(pos, logp, col="#000000", pch=20, cex = 1.5,...))
    with(d.highlight, points(pos, logp, col="#000000", pch=18, cex = 2, ...))
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
manhattan4(tap2)
max.pos = 33032600
min.pos = 32534300
genelist <- read.table("known_genes_buildhg19_chr6_strand.txt", header=T)
genes <- subset(genelist, ( genelist$START > min.pos & genelist$START < max.pos ) | ( genelist$STOP > min.pos & genelist$STOP < max.pos) )

##Remove duplicate genes (picked longest)
rem=c(1,2,13,16,18,19)
all.genes.in.locus = genes[-rem,]
#only expressed genes
notexpr = c(4,6,7,8,9,10,12)
genes.in.locus = all.genes.in.locus[-notexpr,]

range = 5
offset <- ( range * 4 / 3 ) - range
big.range <- range + offset 
#Expressed genes green
for ( i in c(3,7)) { 
  if ( genes.in.locus[i,]$STRAND == "+" ) {
    arrows(max(genes.in.locus[i,]$START, min.pos), -1.8, min(genes.in.locus[i,]$STOP, max.pos), -1.8, length=0.05, lwd=2, code=2, lty="solid", col="darkgreen")
  } else {    
    arrows(max(genes.in.locus[i,]$START, min.pos), -1.8, min(genes.in.locus[i,]$STOP, max.pos), -1.8, length=0.05, lwd=2, code=1, lty="solid", col="darkgreen")
  }
  if ( ! is.na(genes.in.locus[i,]$GENE) ) {
    if ( genes.in.locus[i,]$STRAND == "+" ) {
      text(genes.in.locus[i,]$START + (genes.in.locus[i,]$SIZE / 2), -1.8 + (.3), labels=genes.in.locus[i,]$GENE, cex=0.8)
    } else {    
      text(genes.in.locus[i,]$START + (genes.in.locus[i,]$SIZE / 2), -1.8 + (.3), labels=genes.in.locus[i,]$GENE, cex=0.8)
    }}
  
}
##To shift BRD2 over
for ( i in 1) { 
  
  arrows(max(genes.in.locus[i,]$START, min.pos), -1.8, min(genes.in.locus[i,]$STOP, max.pos), -1.8, length=0.05, lwd=2, code=1, lty="solid", col="darkgreen")
  
  if ( ! is.na(genes.in.locus[i,]$GENE) ) {
    
    text((genes.in.locus[i,]$START+13000) + (genes.in.locus[i,]$SIZE / 2), -1.8 + (.30), labels=genes.in.locus[i,]$GENE, cex=0.8)
  }
  
}

for ( i in 4) { 
  
  arrows(max(genes.in.locus[i,]$START, min.pos), -1.3, min(genes.in.locus[i,]$STOP, max.pos), -1.3, length=0.05, lwd=2, code=1, lty="solid", col="darkgreen")
  
  if ( ! is.na(genes.in.locus[i,]$GENE) ) {
    
    text((genes.in.locus[i,]$START-13000) + (genes.in.locus[i,]$SIZE / 2), -1, labels=genes.in.locus[i,]$GENE, cex=0.8)
  }
  
}


for ( i in c(2,6,8 )) { 
  arrows(max(genes.in.locus[i,]$START, min.pos), -1.3, min(genes.in.locus[i,]$STOP, max.pos), -1.3, length=0.05, lwd=2, code=1, lty="solid", col="darkgreen")
  
  if ( ! is.na(genes.in.locus[i,]$GENE) ) {
    
    text(genes.in.locus[i,]$START + (genes.in.locus[i,]$SIZE / 2), -1, labels=genes.in.locus[i,]$GENE, cex=0.8)
  }
  
}


abline(h=0, col="black", lty=2)
abline(h=-2, col="black", lty=1)
abline(v=32534300, col="black", lty=1)
chr=6
recomb <- read.table(paste("genetic_map_chr", chr, ".txt", sep=""), header=T)
keep.recomb <- subset(recomb, recomb[,1] > min.pos & recomb[,1] < max.pos)
ystart.recomb <- - offset + (big.range / 8)
lines(keep.recomb[,1], ystart.recomb + ( ( keep.recomb[,2] / 60 ) * ( 6 * big.range / 8 )), type="l", col="black", lwd=1, xlim=c(min.pos, max.pos), ylim=c(-offset,range), xlab="", ylab="",  axes=F)
axis(4, at=c( ystart.recomb, ystart.recomb + (big.range / 4), ystart.recomb + ( 2 * big.range / 4),ystart.recomb + ( 3 * big.range / 4)), labels=c("0","20","40", "60"), las=1)
mtext("Recombination rate (cM/Mb)", side=4, at=(-offset+big.range/2), line=2)