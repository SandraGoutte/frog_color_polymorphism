################################################################################
###     Plot  in R
################################################################################
### What this script does:
### 1. Import the output from Betascan (1 file per chromosome)
### 2. Obtain the 1st and 2nd percentiles highest values across the genome
### 3. Plot betascores

################################################################################
### 1. Import the output from Betascan (1 file per chromosome)
################################################################################
betascores_chr1 = read.table("robeensis_chr1.w3000_thetamap.betascores.txt", sep="\t", header=T)
betascores_chr2 = read.table("robeensis_chr2.w3000_thetamap.betascores.txt", sep="\t", header=T)
betascores_chr3 = read.table("robeensis_chr3.w3000_thetamap.betascores.txt", sep="\t", header=T)
betascores_chr4 = read.table("robeensis_chr4.w3000_thetamap.betascores.txt", sep="\t", header=T)
betascores_chr5 = read.table("robeensis_chr5.w3000_thetamap.betascores.txt", sep="\t", header=T)
betascores_chr6 = read.table("robeensis_chr6.w3000_thetamap.betascores.txt", sep="\t", header=T)
betascores_chr7 = read.table("robeensis_chr7.w3000_thetamap.betascores.txt", sep="\t", header=T)
betascores_chr8 = read.table("robeensis_chr8.w3000_thetamap.betascores.txt", sep="\t", header=T)
betascores_chr9 = read.table("robeensis_chr9.w3000_thetamap.betascores.txt", sep="\t", header=T)
betascores_chr10 = read.table("robeensis_chr10.w3000_thetamap.betascores.txt", sep="\t", header=T)
betascores_chr11 = read.table("robeensis_chr11.w3000_thetamap.betascores.txt", sep="\t", header=T)
betascores_chr12 = read.table("robeensis_chr12.w3000_thetamap.betascores.txt", sep="\t", header=T)

## assign chromosome number 
betascores_chr1$chr = 1
betascores_chr2$chr = 2
betascores_chr3$chr = 3
betascores_chr4$chr = 4
betascores_chr5$chr = 5
betascores_chr6$chr = 6
betascores_chr7$chr = 7
betascores_chr8$chr = 8
betascores_chr9$chr = 9
betascores_chr10$chr = 10
betascores_chr11$chr = 11
betascores_chr12$chr = 12
betascores_all_chr = rbind(betascores_chr1, betascores_chr2, betascores_chr3, betascores_chr4, betascores_chr5, betascores_chr6, betascores_chr7, betascores_chr8, betascores_chr9, betascores_chr10, betascores_chr11, betascores_chr12)

## save values for the entire genome
write.table(betascores_all_chr,"betascores_robeensis_all_chr_w3000.txt", sep="\t")

################################################################################
### 2. Obtain the 1st and 2nd percentiles highest values across the genome
################################################################################
## get the top 1% values across the entire genome
betascores_all_chr_sorted = betascores_all_chr[order(betascores_all_chr$Beta1., decreasing=T),]
top_1percent_Betascores_all_chr = betascores_all_chr_sorted[1:nrow(betascores_all_chr_sorted)/100,]
# 1% = 39.18029

## get the top 2% values across the entire genome
min(betascores_all_chr_sorted[1:(2*nrow(betascores_all_chr_sorted))/100,]$Beta1.) 
## 2% = 29.96764

################################################################################
### 3. Plot betascores
################################################################################
## plot a histogram of betascores
hist(betascores_all_chr$Beta1., breaks=100)

### plots betscores for each chromosome
par(mfrow=c(6,1), mar=c(2, 2, 2, 0))
plot(betascores_chr1$Position, betascores_chr1$Beta1.,pch=16, cex=0.5, col="steelblue",  bty="l", main ="Chr1", 
     ylim=c(min(betascores_chr1$Beta1.), max(betascores_chr1$Beta1.)), xlim=c(0, max(betascores_chr1$Position)))
plot(betascores_chr2$Position, betascores_chr2$Beta1.,pch=16, cex=0.5, col="steelblue",  bty="l", main ="Chr2", 
     ylim=c(min(betascores_chr1$Beta1.), max(betascores_chr2$Beta1.)), xlim=c(0, max(betascores_chr2$Position)))
plot(betascores_chr3$Position, betascores_chr3$Beta1.,pch=16, cex=0.5, col="steelblue",  bty="l", main ="Chr3", 
     ylim=c(min(betascores_chr1$Beta1.), max(betascores_chr3$Beta1.)), xlim=c(0, max(betascores_chr3$Position)))
plot(betascores_chr4$Position, betascores_chr4$Beta1.,pch=16, cex=0.5, col="steelblue",  bty="l", main ="Chr4", 
     ylim=c(min(betascores_chr1$Beta1.), max(betascores_chr4$Beta1.)), xlim=c(0, max(betascores_chr4$Position)))
plot(betascores_chr5$Position, betascores_chr5$Beta1.,pch=16, cex=0.5, col="steelblue",  bty="l", main ="Chr5", 
     ylim=c(min(betascores_chr1$Beta1.), max(betascores_chr5$Beta1.)), xlim=c(0, max(betascores_chr5$Position)))
plot(betascores_chr6$Position, betascores_chr6$Beta1.,pch=16, cex=0.5, col="steelblue",  bty="l", main ="Chr6", 
     ylim=c(min(betascores_chr1$Beta1.), max(betascores_chr6$Beta1.)), xlim=c(0, max(betascores_chr6$Position)))
plot(betascores_chr7$Position, betascores_chr7$Beta1.,pch=16, cex=0.5, col="steelblue",  bty="l", main ="Chr7", 
     ylim=c(min(betascores_chr1$Beta1.), max(betascores_chr7$Beta1.)), xlim=c(0, max(betascores_chr7$Position)))
plot(betascores_chr8$Position, betascores_chr8$Beta1.,pch=16, cex=0.5, col="steelblue",  bty="l", main ="Chr8", 
     ylim=c(min(betascores_chr1$Beta1.), max(betascores_chr8$Beta1.)), xlim=c(0, max(betascores_chr8$Position)))
plot(betascores_chr9$Position, betascores_chr9$Beta1.,pch=16, cex=0.5, col="steelblue",  bty="l", main ="Chr9", 
     ylim=c(min(betascores_chr1$Beta1.), max(betascores_chr9$Beta1.)), xlim=c(0, max(betascores_chr9$Position)))
plot(betascores_chr10$Position, betascores_chr10$Beta1.,pch=16, cex=0.5, col="steelblue",  bty="l", main ="Chr10", 
     ylim=c(min(betascores_chr1$Beta1.), max(betascores_chr10$Beta1.)), xlim=c(0, max(betascores_chr10$Position)))
plot(betascores_chr11$Position, betascores_chr11$Beta1.,pch=16, cex=0.5, col="steelblue",  bty="l", main ="Chr11", 
     ylim=c(min(betascores_chr1$Beta1.), max(betascores_chr11$Beta1.)), xlim=c(0, max(betascores_chr11$Position)))
plot(betascores_chr12$Position, betascores_chr12$Beta1.,pch=16, cex=0.5, col="steelblue",  bty="l", main ="Chr12", 
     ylim=c(min(betascores_chr1$Beta1.), max(betascores_chr12$Beta1.)), xlim=c(0, max(betascores_chr12$Position)))

## plot chr7 region of interest
pdf("chr7_betascores_w3000_zoom7.pdf", width=15, height=4)
plot(betascores_chr7$Position, betascores_chr7$Beta1.,pch=16, cex=0.5, col="steelblue",  bty="l", main ="Chr7", 
     ylim=c(0, max(betascores_chr7$Beta1.)),  xlim=c(26415000,26418000))
abline(h=39.18029, lty=2, col="red")
dev.off()


