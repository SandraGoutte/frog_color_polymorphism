################################################################################
###	      PLOT OUTPUT OF GWAS ON GREEN/BROWN P. ROBEENSIS
################################################################################
### What this script does: 
### 1. import the output of GWAS conducted on green and brown P. robeensis 
### 2. generate and save Manhattan plots

# install the qqman package
install.packages("qqman")
# load library
library(qqman)

# Set working directory
setwd("/path/to/directory")

################################################################################
### 1. load association results from plink analysis
################################################################################
assoc_color=read.table("robeensis_filtered.P2.assoc",header=T)
head(assoc_color)

## check the number of chromosomes (we should have 12)
class(assoc_color$X.CHR)
assoc_color$X.CHR = factor(assoc_color$X.CHR)
length(levels(assoc_color$X.CHR))

## re-order the chromosomes by length: longer to shorter (as it should be!)
{
assoc_color$X.CHR2 = sapply(strsplit(assoc_color$CHR, "_"), "[", 8)
head(assoc_color)
assoc_color$X.CHR2 = 1/(as.numeric(assoc_color$X.CHR2))
head(assoc_color)
assoc_color$X.CHR2 = as.factor(assoc_color$X.CHR2)
head(assoc_color)
levels(assoc_color$X.CHR2)

for(i in 1:length(levels(assoc_color$X.CHR2))){
  levels(assoc_color$X.CHR2)[i] = i
}
head(assoc_color)
class(assoc_color$X.CHR2)
assoc_color$CHR = as.numeric(assoc_color$X.CHR2)
head(assoc_color)
}

## check the association table
head(assoc_color)

### get a list of the 100 most significant SNPs and save it
## order association results based on p-value and extract the first 100 SNPs
best_assoc = head(assoc_color[order(assoc_color$P),] ,100)
## save it
write.table(best_assoc, "GWAS_robeensis_color_best_SNPs_trueCHR.txt", sep=" ",
            col.names= T, quote=F)

### prepare data to plot
## remove all the text in the column CHR of our assoc table
assoc_color$CHR = gsub("[a-zA-Z ]", "", assoc_color$CHR)
assoc_color$CHR=sapply(strsplit(assoc_color$CHR, "_"), "[", 2)
assoc_color$contig=sapply(strsplit(assoc_color$CHR, "_"), "[", 3)

assoc_color= assoc_color[assoc_color$CHR!=0,]
assoc_color$CHR = as.numeric(assoc_color$CHR)
assoc_ordered=assoc_color[order(assoc_color$CHR),]
head(assoc_ordered)
## remove NA 
assoc_ordered2 = assoc_ordered[!is.na(assoc_ordered$P),]

## save it
write.table(assoc_ordered, "assoc_color_older_ordered.txt", sep=" ", 
            col.names= T, quote=F)
#assoc_ordered = read.table("assoc_color_older_ordered.txt",header=T)

################################################################################
### 2. generate and save Manhattan plots
################################################################################
## rename columns
colnames(assoc_ordered)[3] = "BP"
colnames(assoc_ordered)[2] = "SNP"

## keep only the 12 largest scaffolds (= chromosomes)
assoc_ordered2_chr1_12 = assoc_ordered2[assoc_ordered2$CHR==1,]

## create and save the plot
job::job({
jpeg(filename=paste0("Manhattan_Plot_GWAS_robeensis_color_chr1-12", ".jpg"),
     width=20, height=10, units="cm", res=300)
manhattan(assoc_ordered2_chr1_12, col=c("lightgrey", "black"), 
          suggestiveline=F,genomewideline=F, ylim=c(0,25), cex=1, cex.sub=1, 
          cex.axis=1, bty="l", xaxs="i")
dev.off()
})

## create and save a plot zoomed on chr 7 
assoc_CHR7= assoc_ordered2[which(assoc_ordered2$CHR==7),]

job::job({
jpeg(filename=paste0("Manhattan_Plot_GWAS_robeensis_color_zoom_chr7_26300000-26500000", ".jpg"), 
     width=20, height=10, units="cm", res=300)
manhattan(assoc_CHR7, col=c("black"), suggestiveline=F, genomewideline=F, 
          ylim=c(0,25), cex=1, cex.sub=1, cex.axis=1, bty="l", xaxs="i", 
          xlim=c(26300000,26500000))
dev.off()
})

## zoom even more
job::job({
jpeg(filename=paste0("Manhattan_Plot_GWAS_robeensis_color_zoom_chr7_26400000-26450000", ".jpg"),
     width=20, height=10, units="cm", res=300)
manhattan(assoc_CHR7, col="black", suggestiveline=F, genomewideline=F, 
          ylim=c(0,25), cex=1, cex.sub=1, cex.axis=1, bty="l", xaxs="i",
          xlim=c(26400000,26450000))
dev.off()
})
