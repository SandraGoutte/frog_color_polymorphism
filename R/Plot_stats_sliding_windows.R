################################################################################
###	      PLOT SLIDING WINDOW STATISTICS        
################################################################################
### What this script does: 
### 1. load sliding window statistics 
###   1.1 Pi (P. robeensis, P. nana, P. levenorum, and P. erlangeri)
###   1.2 Dxy (P. robeensis, P. nana, P. levenorum, and P. erlangeri)
###   1.3 Fst (P. robeensis, P. nana, P. levenorum, and P. erlangeri)
###   1.4 Tajima's D (P. robeensis)
###   1.5 Betascore (P. robeensis)
###   1.6 GWAS (P. robeensis)
### 2. create a multi-panel plot for all 6 stats in P. robeensis
### 3. create a multi-panel plot for Pi, Dxy, and Fst in the 3 other Ptychadena


################################################################################
### 1. Load data outputed by Pixy using our custom 3kb window with 1kb overlap
################################################################################
## set working directory
setwd("/path/to/directory")

################################################################################
### 1.1 Pi (P. robeensis, P. nana, P. erlangeri, P. levenorum)
################################################################################
## pi for the green and brown  separately
pi = read.table("pixy_4sp_customwindow_chr7_color_pi.txt", 
                header=T)

## get each population (species + color) separately
brobpi = pi[pi$pop=="brown_robeensis",]
brobpi = brobpi[is.na(brobpi$avg_pi)==FALSE,]
grobpi = pi[pi$pop=="green_robeensis",]
grobpi = grobpi[is.na(grobpi$avg_pi)==FALSE,]
bnanpi = pi[pi$pop=="brown_nana",]
bnanpi = bnanpi[is.na(bnanpi$avg_pi)==FALSE,]
gnanpi = pi[pi$pop=="green_nana",]
gnanpi = gnanpi[is.na(gnanpi$avg_pi)==FALSE,]
blevpi = pi[pi$pop=="brown_levenorum",]
blevpi = blevpi[is.na(blevpi$avg_pi)==FALSE,]
glevpi = pi[pi$pop=="green_levenorum",]
glevpi = glevpi[is.na(glevpi$avg_pi)==FALSE,]
berlpi = pi[pi$pop=="brown_erlangeri",]
berlpi = berlpi[is.na(berlpi$avg_pi)==FALSE,]
gerlpi = pi[pi$pop=="green_erlangeri",]
gerlpi = gerlpi[is.na(gerlpi$avg_pi)==FALSE,]

## pi for the green and brown together 
pic = read.table("pixy_4sp_customwindow_chr7_species_pi.txt", header=T)

### get each species separately
  robpi = pi[pi$pop=="robeensis",]
  robpi = robpi[is.na(robpi$avg_pi)==FALSE,]
  nanpi = pi[pi$pop=="nana",]
  nanpi = nanpi[is.na(nanpi$avg_pi)==FALSE,]
  levpi = pi[pi$pop=="levenorum",]
  levpi = levpi[is.na(levpi$avg_pi)==FALSE,]
  erlpi = pi[pi$pop=="erlangeri",]
  erlpi = erlpi[is.na(erlpi$avg_pi)==FALSE,]

## get the average Pi over the entire genome
## load the pi values per chromosome
{
pi_chr1 = read.table("pixy_4sp_customwindow_chr1_species_pi.txt", header=T)
head(pi_chr1)
pi_chr2 = read.table("pixy_4sp_customwindow_chr2_species_pi.txt", header=T)
head(pi_chr2)
pi_chr3 = read.table("pixy_4sp_customwindow_chr3_species_pi.txt", header=T)
pi_chr4 = read.table("pixy_4sp_customwindow_chr4_species_pi.txt", header=T)
pi_chr5 = read.table("pixy_4sp_customwindow_chr5_species_pi.txt", header=T)
pi_chr6 = read.table("pixy_4sp_customwindow_chr6_species_pi.txt", header=T)
pi_chr7 = read.table("pixy_4sp_customwindow_chr7_species_pi.txt", header=T)
pi_chr8 = read.table("pixy_4sp_customwindow_chr8_species_pi.txt", header=T)
pi_chr9 = read.table("pixy_4sp_customwindow_chr9_species_pi.txt", header=T)
pi_chr10 = read.table("pixy_4sp_customwindow_chr10_species_pi.txt", header=T)
pi_chr11 = read.table("pixy_4sp_customwindow_chr11_species_pi.txt", header=T)
pi_chr12 = read.table("pixy_4sp_customwindow_chr12_species_pi.txt", header=T)
}

## join all the pi values
pis = list(pi_chr1,pi_chr2,pi_chr3,pi_chr4,pi_chr5,pi_chr6,pi_chr7,pi_chr8,pi_chr9,pi_chr10,pi_chr11,pi_chr12 )
names(pis)=c("pi_chr1","pi_chr2","pi_chr3","pi_chr4",'pi_chr5','pi_chr6','pi_chr7',"pi_chr8","pi_chr9","pi_chr10","pi_chr11",'pi_chr12' )

robpi_allchr = list()
nanpi_allchr = list()
levpi_allchr = list()
erlpi_allchr = list()

for(i in 1:12){
  robpi_allchr[[i]] = pis[[i]][pis[[i]]$pop=="robeensis",]
  robpi_allchr[[i]] = robpi_allchr[[i]][is.na(robpi_allchr[[i]]$avg_pi)==FALSE,]
  nanpi_allchr[[i]] = pis[[i]][pis[[i]]$pop=="nana",]
  nanpi_allchr[[i]] = nanpi_allchr[[i]][is.na(nanpi_allchr[[i]]$avg_pi)==FALSE,]
  levpi_allchr[[i]] = pis[[i]][pis[[i]]$pop=="levenorum",]
  levpi_allchr[[i]] = levpi_allchr[[i]][is.na(levpi_allchr[[i]]$avg_pi)==FALSE,]
  erlpi_allchr[[i]] = pis[[i]][pis[[i]]$pop=="erlangeri",]
  erlpi_allchr[[i]] = erlpi_allchr[[i]][is.na(erlpi_allchr[[i]]$avg_pi)==FALSE,]
}
    
## get the average (calculation as suggested on the pixy website) for each species
## P. robeensis
head(robpi_allchr[[i]])
diffs_rob = sum(unlist(lapply(robpi_allchr, function(df) df$count_diffs)), 
                na.rm = TRUE)
comparisons_rob = sum(unlist(lapply(robpi_allchr, 
                                    function(df) df$count_comparisons)),
                      na.rm = TRUE)
mean_pi_rob = diffs_rob/comparisons_rob # 0.02429921
#mean_pi_rob = 0.02429921

## P. nana
diffs_nan = sum(unlist(lapply(nanpi_allchr, function(df) df$count_diffs)), 
                na.rm = TRUE)
comparisons_nan = sum(unlist(lapply(nanpi_allchr,
                                    function(df) df$count_comparisons)), 
                      na.rm = TRUE)
mean_pi_nan = diffs_nan/comparisons_nan # 0.08239055
#mean_pi_nan = 0.08239055

## P. levenorum
diffs_lev = sum(unlist(lapply(levpi_allchr, function(df) df$count_diffs)), 
                na.rm = TRUE)
comparisons_lev = sum(unlist(lapply(levpi_allchr,
                                    function(df) df$count_comparisons)),
                      na.rm = TRUE)
mean_pi_lev = diffs_lev/comparisons_lev # 0.07400841
#mean_pi_lev = 0.07400841

## P. erlangeri
diffs_erl = sum(unlist(lapply(erlpi_allchr, function(df) df$count_diffs)),
                na.rm = TRUE)
comparisons_erl = sum(unlist(lapply(erlpi_allchr, 
                                    function(df) df$count_comparisons)), 
                      na.rm = TRUE)
mean_pi_erl = diffs_erl/comparisons_erl # 0.08715879. dxy=0.08587214
#mean_pi_erl = 0.08715879

## get the average and the 90th and 95th percentile values for Pi
## P. robeensis
avg_pi_rob = unlist(lapply(robpi_allchr, function(df) df$avg_pi))
percentiles_pi_rob <- quantile(avg_pi_rob, probs = c(0.90, 0.95))

## P. nana
avg_pi_nan = unlist(lapply(nanpi_allchr, function(df) df$avg_pi))
percentiles_pi_nan <- quantile(avg_pi_nan, probs = c(0.90, 0.95))

## P. levenorum
avg_pi_lev = unlist(lapply(levpi_allchr, function(df) df$avg_pi))
percentiles_pi_lev <- quantile(avg_pi_lev, probs = c(0.90, 0.95))

## P. erlangeri
avg_pi_erl = unlist(lapply(erlpi_allchr, function(df) df$avg_pi))
percentiles_pi_erl <- quantile(avg_pi_erl, probs = c(0.90, 0.95))

################################################################################
### 1.2 Dxy  (P. robeensis, P. nana, P. erlangeri, P. levenorum)
################################################################################
## load the Dxy data for all the chromosomes
for(i in c(1:12)){
  file=paste0("pixy_4sp_customwindow_chr", i, "_color_dxy.txt")
  assign(paste0("dxy_chr", i),read.table(file,header=T))
}

## put everything in a list and separate per species
dxys = list(dxy_chr1,dxy_chr2,dxy_chr3,dxy_chr4,dxy_chr5,dxy_chr6,dxy_chr7,
            dxy_chr8,dxy_chr9,dxy_chr10,dxy_chr11,dxy_chr12 )
names(dxys)=c("dxy_chr1","dxy_chr2","dxy_chr3","dxy_chr4",'dxy_chr5',
              'dxy_chr6','dxy_chr7',"dxy_chr8","dxy_chr9","dxy_chr10",
              "dxy_chr11",'dxy_chr12' )

robdxy_allchr = list()
nandxy_allchr = list()
levdxy_allchr = list()
erldxy_allchr = list()

for(i in 1:12){
  robdxy_allchr[[i]] = dxys[[i]][dxys[[i]]$pop1=="green_robeensis" & 
                                   dxys[[i]]$pop2=="brown_robeensis",]
  robdxy_allchr[[i]] = robdxy_allchr[[i]][is.na(robdxy_allchr[[i]]$avg_dxy)==FALSE,]
  nandxy_allchr[[i]] = dxys[[i]][dxys[[i]]$pop1=="brown_nana" &
                                   dxys[[i]]$pop2=="green_nana",]
  nandxy_allchr[[i]] = nandxy_allchr[[i]][is.na(nandxy_allchr[[i]]$avg_dxy)==FALSE,]
  levdxy_allchr[[i]] = dxys[[i]][dxys[[i]]$pop1=="brown_levenorum" &
                                   dxys[[i]]$pop2=="green_levenorum",]
  levdxy_allchr[[i]] = levdxy_allchr[[i]][is.na(levdxy_allchr[[i]]$avg_dxy)==FALSE,]
  erldxy_allchr[[i]] = dxys[[i]][dxys[[i]]$pop1=="brown_erlangeri" &
                                   dxys[[i]]$pop2=="green_erlangeri",]
  erldxy_allchr[[i]] = erldxy_allchr[[i]][is.na(erldxy_allchr[[i]]$avg_dxy)==FALSE,]
}

## get the average following the method suggested on the Pixy website
# P. robeensis
diffs_rob_dxy = sum(unlist(lapply(robdxy_allchr, function(df) df$count_diffs)),
                    na.rm = TRUE)
comparisons_rob_dxy = sum(unlist(lapply(robdxy_allchr, 
                                        function(df) df$count_comparisons)),
                          na.rm = TRUE)
mean_dxy_rob = diffs_rob_dxy/comparisons_rob_dxy #0.02434962

# P. nana
diffs_nan_dxy = sum(unlist(lapply(nandxy_allchr, function(df) df$count_diffs)),
                    na.rm = TRUE)
comparisons_nan_dxy = sum(unlist(lapply(nandxy_allchr, 
                                        function(df) df$count_comparisons)), 
                          na.rm = TRUE)
mean_dxy_nan = diffs_nan_dxy/comparisons_nan_dxy #0.0833317

# P. levenorum
diffs_lev_dxy = sum(unlist(lapply(levdxy_allchr, function(df) df$count_diffs)),
                    na.rm = TRUE)
comparisons_lev_dxy = sum(unlist(lapply(levdxy_allchr, 
                                        function(df) df$count_comparisons)),
                          na.rm = TRUE)
mean_dxy_lev = diffs_lev_dxy/comparisons_lev_dxy #0.0756415

# P. erlangeri
diffs_erl_dxy = sum(unlist(lapply(erldxy_allchr, function(df) df$count_diffs)),
                    na.rm = TRUE)
comparisons_erl_dxy = sum(unlist(lapply(erldxy_allchr,
                                        function(df) df$count_comparisons)),
                          na.rm = TRUE)
mean_dxy_erl = diffs_erl_dxy/comparisons_erl_dxy #0.08848693

################################################################################
### 1.3 Fst  (P. robeensis, P. nana, P. erlangeri, P. levenorum)
################################################################################
## load the Fst data for all the chromosomes
for(i in c(1:12)){
 file=paste0"pixy_4sp_customwindow_chr", i, "_color_fst.txt"
 assign(paste0("fst_chr", i),read.table(file,header=T))
}

## put everything in a list and separate per species
fsts = list(fst_chr1,fst_chr2,fst_chr3,fst_chr4,fst_chr5,fst_chr6,fst_chr7,
            fst_chr8,fst_chr9,fst_chr10,fst_chr11,fst_chr12 )
names(fsts)=c("fst_chr1","fst_chr2","fst_chr3","fst_chr4",'fst_chr5',
              'fst_chr6','fst_chr7',"fst_chr8","fst_chr9","fst_chr10",
              "fst_chr11",'fst_chr12' )

robfst_allchr = list()
nanfst_allchr = list()
levfst_allchr = list()
erlfst_allchr = list()

for(i in 1:12){
  robfst_allchr[[i]] = fsts[[i]][fsts[[i]]$pop1=="green_robeensis" &
                                   fsts[[i]]$pop2=="brown_robeensis",]
  robfst_allchr[[i]] = robfst_allchr[[i]][is.na(robfst_allchr[[i]]$avg_wc_fst)==FALSE,]
  nanfst_allchr[[i]] = fsts[[i]][fsts[[i]]$pop1=="brown_nana" &
                                   fsts[[i]]$pop2=="green_nana",]
  nanfst_allchr[[i]] = nanfst_allchr[[i]][is.na(nanfst_allchr[[i]]$avg_wc_fst)==FALSE,]
  levfst_allchr[[i]] = fsts[[i]][fsts[[i]]$pop1=="brown_levenorum" &
                                   fsts[[i]]$pop2=="green_levenorum",]
  levfst_allchr[[i]] = levfst_allchr[[i]][is.na(levfst_allchr[[i]]$avg_wc_fst)==FALSE,]
  erlfst_allchr[[i]] = fsts[[i]][fsts[[i]]$pop1=="brown_erlangeri" &
                                   fsts[[i]]$pop2=="green_erlangeri",]
  erlfst_allchr[[i]] = erlfst_allchr[[i]][is.na(erlfst_allchr[[i]]$avg_wc_fst)==FALSE,]
}

## get the average values per species
head(robfst_allchr[[i]])
mean_fst_rob = mean(unlist(lapply(robfst_allchr, function(df) df$avg_wc_fst)),
                    na.rm = TRUE) #-0.007191086
mean_fst_nan = mean(unlist(lapply(nanfst_allchr, function(df) df$avg_wc_fst)),
                    na.rm = TRUE) #-0.00575761
mean_fst_lev = mean(unlist(lapply(levfst_allchr, function(df) df$avg_wc_fst)),
                    na.rm = TRUE)#-0.02642259
mean_fst_erl = mean(unlist(lapply(erlfst_allchr, function(df) df$avg_wc_fst)), 
                    na.rm = TRUE) #-0.007119631


################################################################################
### 1.4 Tajima's D (P. robeensis)
################################################################################
## load the daa
tajima = read.table("TajimasD_robeensis_highcov.txt", header=T)

## get the chr 7 values only
tajima_ch7 = tajima[tajima$CHROM=="Scaffold_4__1_contigs__length_105873206",]

################################################################################
### 1.5 Betascores (P. robeensis)
################################################################################
## load the data
betascores_all_chr =read.table("betascores_all_chr_w3000.txt", sep="\t",
                               header=T)
betascores_chr7 = read.table("robeensis_chr7.w3000.betascores.txt", sep="\t",
                             header=T)
betascores_chr7$chr = 7

## get the top values for the whole genome
betascores_all_chr_sorted = betascores_all_chr[order(betascores_all_chr$Beta1., decreasing=T),]
top_1percent_Betascores_all_chr = betascores_all_chr_sorted[1:nrow(betascores_all_chr_sorted)/1000,]
  
min(betascores_all_chr_sorted[1:(2*nrow(betascores_all_chr_sorted))/100,]$Beta1.) 
# top 2% threshold = 29.96764, top 1% threshold = 39.18029


################################################################################
### 1.6 GWAS (P. robeensis)
################################################################################
## load the data
assoc_CHR7 = read.csv("assoc_CHR7.csv", header=T)
best_snps = as.character(read.csv("best_snps.csv", header=T)[,2])


################################################################################
### 2. Plot P. robeensis stats in a single multi-panel plot
################################################################################
## load libraries
library(viridis)
library(ggsci)

## define the parameters
xmin = 26300000
xmax = 26500000
col_rob = alpha("#666666",.6)
col_green = "#8DB634"
col_brown = "#9F583A"
cex_point = 2
lwd=1.5

### plot 
pdf(file=("Probeensis_stats.pdf"), width=16, height=4)
{
par(mfrow=c(2,3), mar=c(6, 4, 5, 2))

## Pi
ymin=0
ymax=0.3
plot(NULL, axes = FALSE, xlab = "", ylab = "", xlim = c(xmin,xmax), ylim = c(ymin, ymax), yaxs="i")
abline(h=c(0,0.1, 0.2, 0.3), col="gray95")
rect(xleft = 26411500, xright = 26422300, ybottom = -10, ytop = 65, col=alpha("grey", .3), border=NA)
par(new=T)
plot(grobpi$window_pos_1+1500, grobpi$avg_pi,main="Pi",yaxp=c(0,0.3, 3), las=1, cex.lab=2,cex.axis=2, ylab="P. robeensis", type="l", bty="l", lwd=lwd, col=col_green,xlim=c(xmin,xmax), xaxs="i", yaxs="i", ylim=c(ymin, ymax))
lines(robpi$window_pos_1+1500,robpi$avg_pi,las=1, type="l", bty="l", lwd=1.5, col="#666666",xlim=c(xmin,xmax), ylim=c(ymin, ymax), main="Pi",xaxs="i")
lines(brobpi$window_pos_1+1500, brobpi$avg_pi, type="l", bty="l", lwd=lwd, col=col_brown, xaxs="i", xaxt = "n")
abline(h=mean_pi_rob, lty=2)
abline(h=percentiles_pi_rob[[2]], col="red", lty=2)

## Dxy
ymin=0
ymax=0.3
plot(NULL, axes = FALSE, xlab = "", ylab = "", xlim = c(xmin,xmax), ylim = c(ymin,ymax), yaxs="i")
abline(h=c(0,0.1, 0.2, 0.3, 0.4), col="gray95")
rect(xleft = 26411500, xright = 26422300, ybottom = -10, ytop = 65, col=alpha("grey", .3), border=NA)
par(new=T)
plot(rob_dxy$window_pos_1+1500, rob_dxy$avg_dxy, yaxp=c(0,0.3, 3), cex.lab=2,cex.axis=2,main="Dxy green vs. brown", las=1, ylab="P. robeensis", 
type="l", bty="l", lwd=lwd, col=col_rob,xlim=c(xmin,xmax),  ylim = c(ymin,ymax), xaxs="i", yaxs="i")
abline(h=mean_dxy_rob, lty=2)

## Fst
ymin=-1
ymax=1
plot(NULL, axes = FALSE, xlab = "", ylab = "", xlim = c(xmin,xmax), ylim = c(ymin,ymax), yaxs="i")
abline(h=seq(from = ymin, to = ymax, length.out = 5), col="gray95")
rect(xleft = 26411500, xright = 26422300, ybottom = -10, ytop = 65, col=alpha("grey", .3), border=NA)
par(new=T)
plot(rob_fst$window_pos_1+1500, rob_fst$avg_wc_fst,main="Fst green vs. brown", cex.lab=2,cex.axis=2, las=1, xaxs="i", type="l", ylim=c(-1,1), bty="l", lwd=lwd, col=col_rob,xlim=c(xmin,xmax), ylab="P. robeensis", xaxs="i", yaxs="i", cex=cex_point)
abline(h=mean_fst_rob, lty=2)

## Tajima's D
lwd=2
plot(NULL, axes=F, xlab = "", ylab = "", xlim = c(xmin,xmax), ylim = c(-3,3), yaxs="i")
rect(xleft = 26411500, xright = 26422300, ybottom = -10, ytop = 65, col=alpha("grey", .3), border=NA)
abline(h=3, lty=1,col=alpha("grey", .3))
abline(h=2, lty=1,col=alpha("grey", .3))
abline(h=1, lty=1,col=alpha("grey", .3))
abline(h=-1, lty=1,col=alpha("grey", .3))
abline(h=-2, lty=1,col=alpha("grey", .3))
abline(h=0, lty=1, col=alpha("grey", .3))
par(new=T)
plot(tajima_ch7$BIN_START+1500, tajima_ch7$TajimaD, 
main="", type="l", bty="l", lwd=lwd, col=col_rob,
xlim=c(xmin,xmax), ylab="Tajima's D", xaxs="i",ylim = c(-3,3), las=1, yaxs="i",yaxp=c(-3,3, 6), cex.lab=2,cex.axis=2) # xaxt = "n" to suppress the axis
# get the genome-wide average value for Tajima's D
abline(h=mean(tajima$TajimaD), lty=2, col="black")

## Betascores 
plot(NULL, axes = FALSE, xlab = "", ylab = "", xlim = c(xmin,xmax), ylim = c(0, 70), yaxs="i")
rect(xleft = 26411500, xright = 26422300, ybottom = -10, ytop = 70, col=alpha("grey", .3), border=NA)
abline(h=70, lty=1,col=alpha("grey", .3))
abline(h=60, lty=1,col=alpha("grey", .3))
abline(h=50, lty=1,col=alpha("grey", .3))
abline(h=40, lty=1,col=alpha("grey", .3))
abline(h=30, lty=1,col=alpha("grey", .3))
abline(h=20, lty=1,col=alpha("grey", .3))
abline(h=10, lty=1,col=alpha("grey", .3))
par(new=T)
plot(betascores_chr7$Position, betascores_chr7$Beta1.,pch=16,cex.lab=2,cex.axis=2, cex=cex_point-.5,
    col=col_rob,  bty="l", main ="",  xaxs="i",
    ylim=c(0, 70),  xlim=c(xmin,xmax), cex.sub=1, ylab="Beta",las=1, yaxs="i")
abline(h=39.18029, lty=2, col="black") 
  
## GWAS
plot(NULL, axes = FALSE, xlab = "", ylab = "", xlim = c(xmin,xmax), ylim = c(0, 25), xaxs="i", yaxs="i")
rect(xleft = 26411500, xright = 26422300, ybottom = -10, ytop = 35, col=alpha("grey", .3), border=NA)
abline(h=25, lty=1,col=alpha("grey", .3))
abline(h=20, lty=1,col=alpha("grey", .3))
abline(h=15, lty=1,col=alpha("grey", .3))
abline(h=10, lty=1,col=alpha("grey", .3))
abline(h=5, lty=1,col=alpha("grey", .3))
abline(h=0, lty=1,col="black")
par(new=T)
manhattan(assoc_CHR7, col=c(col_rob), suggestiveline=F, 
            genomewideline=F, ylim=c(0,25), cex=cex_point,
            cex.sub=1, cex.axis=1, xlim=c(xmin,xmax), highlight=best_snps, ylab="GWAS", xaxs="i", yaxs="i",cex.lab=2,cex.axis=2) #  xaxt = "n" to suppress the axes
  
}
dev.off()

################################################################################
### 3. Plot Pi, Dxy and Fst for the 3 other Ptychadena species
################################################################################
## set the parameters
xmin = 26300000
xmax = 26500000
colpi="#666666"

## plot
pdf(file=("Ptychadena_3sp_stats.pdf"), width=16, height=4)
par(mfrow=c(3,3), mar=c(2, 4, 1, 2))

## pi
ymin=0
ymax=0.4
# P. nana
plot(NULL, axes = FALSE, xlab = "", ylab = "", xlim = c(xmin,xmax), 
     ylim=c(ymin, ymax),xaxs="i")
abline(h=c(0,0.1, 0.2, 0.3, 0.4), col="gray95")
rect(xleft = 26411500, xright = 26422300, ybottom = -10, ytop = 65,
     col=alpha("grey", .3), border=NA)
par(new=T)
plot(gnanpi$window_pos_1+1500, gnanpi$avg_pi, las=1,  cex.lab=2,
     cex.axis=2,ylab="P. nana",type="l", bty="l", lwd=lwd, col=col_green,
     xlim=c(xmin,xmax), xaxs="i", xaxt = "n", ylim=c(ymin,ymax))
lines(bnanpi$window_pos_1+1500, bnanpi$avg_pi, type="l", bty="l", lwd=lwd,
      col=col_brown, xaxs="i", xaxt = "n")
lines(nanpi$window_pos_1+1500, nanpi$avg_pi,las=1, type="l", bty="l", lwd=1.5,
      col=colpi,xlim=c(xmin,xmax), ylim=c(ymin, ymax),xaxs="i", xaxt = "n")
abline(h=mean_pi_nan, lty=2)
abline(h=percentiles_pi_nan[[2]], col="red", lty=2)
# P. levenorum
plot(NULL, axes = FALSE, xlab = "", ylab = "", xlim = c(xmin,xmax),
     ylim=c(ymin, ymax),xaxs="i")
abline(h=c(0,0.1, 0.2, 0.3, 0.4), col="gray95")
rect(xleft = 26411500, xright = 26422300, ybottom = -10, ytop = 65,
     col=alpha("grey", .3), border=NA)
par(new=T)
plot(glevpi$window_pos_1+1500, glevpi$avg_pi, las=1, cex.lab=2,cex.axis=2,
     ylab="P. levenorum",type="l", bty="l", lwd=lwd, col=col_green,
     xlim=c(xmin,xmax), xaxs="i", xaxt = "n", ylim=c(ymin,ymax))
lines(blevpi$window_pos_1+1500, blevpi$avg_pi, type="l", bty="l", lwd=lwd,
      col=col_brown, xaxs="i", xaxt = "n")
lines(levpi$window_pos_1+1500, levpi$avg_pi,las=1, type="l", bty="l", lwd=1.5,
      col=colpi,xlim=c(xmin,xmax), ylim=c(ymin, ymax),xaxs="i", xaxt = "n")
abline(h=mean_pi_lev, lty=2)
abline(h=percentiles_pi_lev[[2]], col="red", lty=2)
# P. erlangeri
plot(NULL, axes = FALSE, xlab = "", ylab = "", xlim = c(xmin,xmax),
     ylim=c(ymin, ymax),xaxs="i")
abline(h=c(0,0.1, 0.2, 0.3, 0.4), col="gray95")
rect(xleft = 26411500, xright = 26422300, ybottom = -10, ytop = 65,
     col=alpha("grey", .3), border=NA)
par(new=T)
plot(gerlpi$window_pos_1+1500, gerlpi$avg_pi, las=1, cex.lab=2,cex.axis=2,
     ylab="P. erlangeri",type="l", bty="l", lwd=lwd, col=col_green,
     xlim=c(xmin,xmax), xaxs="i",  ylim=c(ymin,ymax),xaxt = "n")
lines(berlpi$window_pos_1+1500, berlpi$avg_pi, type="l", bty="l", lwd=lwd, 
      col=col_brown, xaxs="i", xaxt = "n")
lines(erlpi$window_pos_1+1500, erlpi$avg_pi,las=1, type="l", bty="l", lwd=1.5, 
      col=colpi,xlim=c(xmin,xmax), ylim=c(ymin, ymax), xaxs="i", xaxt = "n")
abline(h=mean_pi_erl, lty=2)
abline(h=percentiles_pi_erl[[2]], col="red", lty=2)

## Dxy
ymin=0
ymax=0.3
# P. nana
plot(NULL, axes = FALSE, xlab = "", ylab = "", xlim = c(xmin,xmax), 
     ylim = c(ymin,ymax))
abline(h=c(0,0.1, 0.2, 0.3, 0.4), col="gray95")
rect(xleft = 26411500, xright = 26422300, ybottom = -10, ytop = 65,
     col=alpha("grey", .3), border=NA)
par(new=T)
plot(nan_dxy$window_pos_1+1500, nan_dxy$avg_dxy, yaxp=c(0,0.3, 3), cex.lab=2,
     cex.axis=2, las=1, ylab="P. nana",type="l", bty="l",lwd=lwd, col=col_nan,
     xlim=c(xmin,xmax), ylim = c(ymin,ymax), xaxs="i", xaxt = "n")
abline(h=mean_dxy_nan, lty=2)
# P. levenorum
plot(NULL, axes = FALSE, xlab = "", ylab = "", xlim = c(xmin,xmax), 
     ylim = c(ymin,ymax))
abline(h=c(0,0.1, 0.2, 0.3, 0.4), col="gray95")
rect(xleft = 26411500, xright = 26422300, ybottom = -10, ytop = 65, 
     col=alpha("grey", .3), border=NA)
par(new=T)
plot(lev_dxy$window_pos_1+1500, lev_dxy$avg_dxy,yaxp=c(0,0.3, 3),cex.lab=2,
     cex.axis=2, las=1,ylab="P. levenorum",type="l", bty="l", lwd=lwd, 
     col=col_lev,xlim=c(xmin,xmax),  ylim = c(ymin,ymax), xaxs="i", xaxt = "n")
abline(h=mean_dxy_lev, lty=2)
# P. erlangeri
plot(NULL, axes = FALSE, xlab = "", ylab = "", xlim = c(xmin,xmax), 
     ylim = c(ymin,ymax))
abline(h=c(0,0.1, 0.2, 0.3, 0.4), col="gray95")
rect(xleft = 26411500, xright = 26422300, ybottom = -10, ytop = 65,
     col=alpha("grey", .3), border=NA)
par(new=T)
plot(erl_dxy$window_pos_1+1500, erl_dxy$avg_dxy,yaxp=c(0,0.3, 3), 
     cex.lab=2,cex.axis=2, las=1,ylab="P. erlangeri",type="l", bty="l",
     lwd=lwd, col=col_erl,xlim=c(xmin,xmax),  ylim = c(ymin,ymax), xaxs="i",
     xaxt = "n")
abline(h=mean_dxy_erl, lty=2)

## Fst
ymin=-1
ymax=1
# P. nana
plot(NULL, axes = FALSE, xlab = "", ylab = "", xlim = c(xmin,xmax), 
     ylim = c(ymin,ymax))
abline(h=seq(from = ymin, to = ymax, length.out = 5), col="gray95")
rect(xleft = 26411500, xright = 26422300, ybottom = -10, ytop = 65, 
     col=alpha("grey", .3), border=NA)
par(new=T)
plot(nan_fst$window_pos_1+1500, nan_fst$avg_wc_fst,cex.lab=2,cex.axis=2, 
     las=1,ylab="P. nana",type="l", bty="l", lwd=lwd, ylim=c(-1,1),col=col_nan,
     xlim=c(xmin,xmax),xaxs="i")
abline(h=mean_fst_nan, lty=2)
# P. levenorum
plot(NULL, axes = FALSE, xlab = "", ylab = "", xlim = c(xmin,xmax),
     ylim = c(ymin,ymax))
abline(h=seq(from = ymin, to = ymax, length.out = 5), col="gray95")
rect(xleft = 26411500, xright = 26422300, ybottom = -10, ytop = 65, 
     col=alpha("grey", .3), border=NA)
par(new=T)
plot(lev_fst$window_pos_1+1500, lev_fst$avg_wc_fst,cex.lab=2,cex.axis=2,  
     las=1,ylab="P. levenorum",type="l", bty="l", lwd=lwd, ylim=c(-1,1),
     col=col_lev,xlim=c(xmin,xmax),xaxs="i")
abline(h=mean_fst_lev, lty=2)
# P. erlangeri
plot(NULL, axes = FALSE, xlab = "", ylab = "", xlim = c(xmin,xmax), 
     ylim = c(ymin,ymax))
abline(h=seq(from = ymin, to = ymax, length.out = 5), col="gray95")
rect(xleft = 26411500, xright = 26422300, ybottom = -10, ytop = 65, 
     col=alpha("grey", .3), border=NA)
par(new=T)
plot(erl_fst$window_pos_1+1500, erl_fst$avg_wc_fst, cex.lab=2,cex.axis=2,
     las=1,ylab="P. erlangeri",type="l", bty="l", lwd=lwd,ylim=c(-1,1), 
     col=col_erl,xlim=c(xmin,xmax), xaxs="i")
abline(h=mean_fst_erl, lty=2)
dev.off()