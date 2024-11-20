################################################################################
###	         Differential expression of FoxD3 in P. robeensis
################################################################################
### What this script does:
### 1. reads RNAseq and phenotype data
### 2. Compares FoxD3 expression levels between brown and green-genotyped  
### individuals across 4 development stages
### 3. Compares FoxD3 expression levels between inside vs outside the stripe 
### in the wide striped individuals 

### Data:
### countData.csv raw count data from RNAseq experiment
### phenotypes2.csv comma delimited phenotype data with row names as
### samples and each column as a trait or sample descriptor
### e.g. stage, color, stripe, zone

### load libraries
library(Glimma)
library(edgeR)
library(DESeq2)
library(apeglm)
library(ggplot2)
library(reshape2)

################################################################################
### 1. Read RNAseq and phenotype data
################################################################################
## list the individual raw count
countData <- as.matrix(read.csv("countData.csv", row.names=1)) 

## load phenotype data
colData = read.csv("phenotypes2.csv", row.names=1)
  
## make sure that the samples for both matrices are in the same order and with the same names
head(colData)  
nrow(colData)
head(countData)
ncol(countData)

## remove the samples SB1021BS and SB1030BS because they are hybrids
## and 15-138 and 15-139 because they are ventral skin
ncol(countData) #62
countData = countData[,-which(colnames(countData)=="sg5533_L20220921.1_SB1021_BS_S30")]
countData = countData[,-which(colnames(countData)=="sg5533_L20220921.1_SB1030_BS_S39")]
countData = countData[,-which(colnames(countData)=="Sample_15_138")]
countData = countData[,-which(colnames(countData)=="Sample_15_139")]
ncol(countData) #58
  
## check that all samples IDs are the same in both tables
all(rownames(colData) %in% colnames(countData))
rownames(colData)
colnames(countData)
  
## they are in the same order, let's make sure their spellings match
colnames(countData) = rownames(colData)
all(rownames(colData) %in% colnames(countData))
rownames(colData)
colnames(countData)
  
## check the class of the phenotypes
colnames(colData)
class(colData$genotype_color)

## change phenotypes into factors
colData$color
colData$color = as.factor(colData$color)
colData$stripe = as.factor(colData$stripe)
colData$stage = as.factor(colData$stage)
colData$zone = as.factor(colData$zone)
colData$side = as.factor(colData$side)
colData$genotype_color = as.factor(colData$genotype_color)
colData$stageskincolor = as.factor(colData$stageskincolor)
  
## create a color*zone variable
colData$colorzone=paste(colData$genotype_color,colData$zone, sep=".")
colData$colorzone = as.factor(colData$colorzone)

### remove SB755_OUT : very low counts
rownames(colData)
colData2 = colData[-c(which(row.names(colData)=="sg5533_L20220921.1_SB755_OUT_S15")),]
remove = c(which(row.names(colData)=="sg5533_L20220921.1_SB755_OUT_S15"))
countData2 = countData[,-c(remove)]

### remove ventral skin
remove = c(which(colData2$side=="belly"))
colData3 = colData2[-c(remove),]
countData3 = countData2[,-c(remove)]
nrow(colData3) # 51

################################################################################
### 2. Compare FoxD3 expression levels between brown and green-genotyped  
### individuals across 4 development stages
################################################################################
### 2.1 Prepare the data
## for some individuals, we have samples IN and OUT the vertebral stripe
## for this analysis, remove the OUT to have only 1 sample per individual
remove = c(which(colData3$zone=="OUT"))
colData4 = colData3[-c(remove),]
countData4 = countData3[,-c(remove)]
nrow(colData4) # 33

## convert into the DGE list format and normalize the count data
{
  dge <- DGEList(counts=countData4)   
  head(dge)
  
  ## plots before any filtering or normalization
  lcpm <- cpm(dge, log=TRUE)
  par(mar=c(5, 4, 4, 2))
  plotDensities(lcpm, legend = FALSE, main = "Before filtering")
  abline(v = 0, lty = 3)
  
  ## non-specific filter to remove low count features
  # only keep genes which have cpm greater than 1 in at least 1/3 of the samples 
  cpm=cpm(dge)
  keep.exprs <- rowSums(cpm > 1) >= 0.3*nrow(colData4)
  
  ## filter the dge
  dge <- dge[keep.exprs, , keep.lib.sizes = FALSE] 
  dim(dge)
  
  ## plots after filtering 
  lcpm <- cpm(dge, log=TRUE)
  par(mar=c(5, 4, 4, 2))
  plotDensities(lcpm, legend = FALSE, main = "After filtering")
  abline(v = 0, lty = 3)
  
  ## visualize the distributions for each sample after filtering and before normalization
  par(mar=c(8, 4, 4, 2))
  boxplot(lcpm, las=2, main = "Before normalization", cex.axis=0.6, col="steelblue")
  
  # TMM normalization
  dge <- calcNormFactors(dge, method="TMM")
  dim(dge)
  
  ## visualize the distribution of gene expression levels after normalization
  lcpm <- cpm(dge, log=TRUE)
  par(mar=c(8, 4, 4, 2))
  boxplot(lcpm, las=2, main = "After normalization", cex.axis=0.6, col="steelblue")
}

### 2.2 Compare FoxD3 expression levels between green and brown frogs
## calculate design matrix 
design <- model.matrix(~0+ genotype_color + stage, data=colData4)
colnames(design) <- make.names(colnames(design)) # chack validity of names
cont.matrix <- makeContrasts(
  brownvsgreen=genotype_colorbrown - genotype_colorgreen,
  levels=design)

## calculate weights to correct for Poisson count noise due to discrete nature 
## of RNA-seq.
v <- voom(dge,design,plot=TRUE)

## compute moderated t-statistics for all genes
fit <- lmFit(v,design)
fit2 <- contrasts.fit(fit, cont.matrix)
efit2 <- eBayes(fit2)

## focus on the expression of FoxD3 = MSTRG.111671
Foxd3 <- "MSTRG.111671"

## retrieve the coefficients and p-values for FoxD3
gene_fit2 <- efit2[rownames(efit2) == Foxd3, ]

## save the results
write.csv(gene_fit2, file="foxd3_DE_brownvsgreen.csv")
topTable(gene_fit2)
topTreat(gene_fit2)

### 2.3 Plot FoxD3 expression levels
## extract the values for FoxD3 = MSTRG.111671
foxd3 = lcpm[which(rownames(lcpm)=="MSTRG.111671"),]
colData4$foxd3 = foxd3

## reorder the developmental stages factor levels to be chronological 
colData4_new <- colData4                             
colData4_new$stage <- factor(colData4_new$stage,     
                             c("tadpole", "metamorph", "juvenile", "adult"))

## plot the expression level of Foxd3 as function of color and stage
ggplot(colData4_new , aes(x=stage, y=foxd3,  fill=genotype_color)) + 
 geom_violin(alpha=1)+
  geom_point(position = position_jitterdodge(), alpha=1,stroke = .5, shape = 21,
             size=6, aes(color="black", fill=genotype_color))+
  scale_fill_manual(values=c( "brown","forestgreen"))+
  scale_color_manual(values=c("black"))+
  theme_minimal()

################################################################################
### 3. Compare FoxD3 expression levels between inside vs outside the stripe 
### in the wide striped individuals 
### (inside can be green or brown, whereas outside is brown in all individuals)
################################################################################
### 3.1 Prepare the data
colData4 = colData3[-c(which(colData3$stripe=="thin"), 
                       which(colData3$stripe=="unstriped")),]
## convert into the DGE list format and normalize the count data
{
  dge <- DGEList(counts=countData4)   
  head(dge)
  
  ## plots before any filtering or normalization
  lcpm <- cpm(dge, log=TRUE)
  par(mar=c(5, 4, 4, 2))
  plotDensities(lcpm, legend = FALSE, main = "Before filtering")
  abline(v = 0, lty = 3)
  
  ## non-specific filter to remove low count features
  # only keep genes which have cpm greater than 1 in at least 1/3 of the samples 
  cpm=cpm(dge)
  keep.exprs <- rowSums(cpm > 1) >= 0.3*nrow(colData4)
  
  ## filter the dge
  dge <- dge[keep.exprs, , keep.lib.sizes = FALSE] 
  dim(dge)
  
  ## plots after filtering 
  lcpm <- cpm(dge, log=TRUE)
  par(mar=c(5, 4, 4, 2))
  plotDensities(lcpm, legend = FALSE, main = "After filtering")
  abline(v = 0, lty = 3)
  
  ## visualize the distributions for each sample after filtering and before normalization
  par(mar=c(8, 4, 4, 2))
  boxplot(lcpm, las=2, main = "Before normalization", cex.axis=0.6, col="steelblue")
  
  # TMM normalization
  dge <- calcNormFactors(dge, method="TMM")
  dim(dge)
  
  ## visualize the distribution of gene expression levels after normalization
  lcpm <- cpm(dge, log=TRUE)
  par(mar=c(8, 4, 4, 2))
  boxplot(lcpm, las=2, main = "After normalization", cex.axis=0.6, col="steelblue")
}

### 3.2 comparing inside vs outside the stripe
## accounting for individual since each individual has IN and OUT samples
## create a column with individual
colData4 <- colData4 %>%
  mutate(indiv = sapply(strsplit(as.character(X.1), "_"), `[`, 1))
head(colData4)
colData4$colorzone = as.character(as.factor(colData4$colorzone))
colData4$indiv = as.character(as.factor(colData4$indiv))

## prepare a design matrix
design <- model.matrix(~0+ colorzone + stage + indiv, data=colData4)
colnames(design) <- make.names(colnames(design)) # check validity of names
cont.matrix <- makeContrasts(
  brownOUTvsbrownIN=colorzonebrown.OUT - colorzonebrown.IN,
  greenOUTvsgreenIN=colorzonegreen.OUT - colorzonegreen.IN,
  levels=design)

## calculate weights to correct for Poisson count noise due to discrete nature 
## of RNA-seq.
v <- voom(dge,design,plot=TRUE)

## compute moderated t-statistics for all genes
fit <- lmFit(v,design)
fit2 <- contrasts.fit(fit, cont.matrix)
efit2 <- eBayes(fit2)

## focus on the expression of FoxD3 = MSTRG.111671
Foxd3 <- "MSTRG.111671"

## retrieve the coefficients and p-values for FoxD3
gene_fit2 <- efit2[rownames(efit2) == Foxd3, ]

## save the results
write.csv(gene_fit2, file="foxd3_DE_INvsOUT.csv")
topTable(gene_fit2)
topTreat(gene_fit2)

### 3.3 comparing green vs brown IN and green vs brown OUT (not same individuals)
## prepare the design matrix
design <- model.matrix(~0+ colorzone + stage, data=colData4)
colnames(design) <- make.names(colnames(design)) # chack validity of names
cont.matrix <- makeContrasts(
  brownOUTvsgreenOUT=colorzonebrown.OUT - colorzonegreen.OUT,
  brownINvsgreenIN=colorzonebrown.IN - colorzonegreen.IN,
  levels=design)

## calculate weights to correct for Poisson count noise due to discrete nature 
## of RNA-seq.
v <- voom(dge,design,plot=TRUE)

## compute moderated t-statistics for all genes
fit <- lmFit(v,design)
fit2 <- contrasts.fit(fit, cont.matrix)
efit2 <- eBayes(fit2)

## focus on the expression of FoxD3 = MSTRG.111671
Foxd3 <- "MSTRG.111671"

## retrieve the coefficients and p-values for FoxD3
gene_fit2 <- efit2[rownames(efit2) == Foxd3, ]

## save the results
write.csv(gene_fit2, file="foxd3_DE_brownINvsgreenIN_brownOUTvsgreenOUT.csv")
topTable(gene_fit2)
topTreat(gene_fit2)

### 3.4 plot Foxd3 expression levels as a function of genotype and stripe zone
ggplot(colData4_new , aes(x=genotype_color, y=foxd3, fill=colorzone)) + 
  geom_violin()+
  geom_point(position = position_jitterdodge(), alpha=1,stroke = .5, shape = 21,
             size=6, aes(color="black", fill=colorzone))+
  scale_color_manual(values=c("black"))+
  scale_fill_manual(values=c( "burlywood1","brown","forestgreen", "brown"))+
  theme_minimal()
