################################################################################
###	                PLOT RECOMBINATION MAP IN P. ROBEENSIS
################################################################################
### What this script does: 
### 1. reads the recombination map outputted by LDhat
### 2. plots recombination values in the region of interest on chromosome 7

################################################################################
### 1. Load recombination data for chromosome 7
################################################################################
all_chr7 = read.table("chr7_highcov_all.joint.txt", sep=" ", header=F)

# name the columns of our data frame
colnames(all_chr7)=c("Loci",	"Mean_rho",	"Median",	"L95",	"U95", "pos")

################################################################################
### 1. Plot recombination values in the region of interest on chromosome 7
################################################################################
# set the plot
par(mfrow=c(2,2), mar=c(5, 4, 4, 2))
plot(all_chr7$pos, all_chr7$Mean_rho,pch=16, col="white", bty="l", main="Chr 7", xlim=c(xmin, xmax), ylim=c(0,30), xaxs="i", yaxs="i", las=1, xaxt="n", yaxt="n")
abline(h=30, lty=1,col=alpha("grey", .3))
abline(h=20, lty=1,col=alpha("grey", .3))
abline(h=10, lty=1,col=alpha("grey", .3))
abline(h=25, lty=1,col=alpha("grey", .3))
abline(h=15, lty=1,col=alpha("grey", .3))
abline(h=5, lty=1,col=alpha("grey", .3))
abline(h=0, lty=1,col=alpha("grey", .3))
# highlight the region of our association peak
rect(xleft = 26411500, xright = 26422300, ybottom = -10, ytop = 35, col=alpha("grey", .3), border=NA)
par(new=T)
plot(all_chr7$pos, all_chr7$Mean_rho,pch=16, col=alpha("#666666", .6), bty="l", main="Chr 7", xlim=c(xmin, xmax), ylim=c(0,30), xaxs="i", yaxs="i", las=1, cex.axis=2)
# add the chromosome-wide average and 90th percentile value as horizontal lines
abline(h=mean(result$mean_rho, na.rm=T), lty=2, col="black")
abline(h=quantile(all_chr7$Mean_rho, probs =0.9), lty=2, col="red")