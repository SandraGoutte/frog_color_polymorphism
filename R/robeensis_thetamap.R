################################################################################
###	   OBTAIN MUTATION RATE MAP FOR EACH CHROMOSOME TO USE IN BETASCAN
################################################################################
### What this script does:
### 1. import the SNP density files for each chromosome (output from vcftools)
### 2. Obtain a theta map for each chromosome and export it

################################################################################
### 1. import the SNP density files for each chromosome (output from vcftools)
################################################################################
setwd("/.../Betascan/Thetamaps/")
snpdens1 = read.table("robeensis_chr1_12kb_window.snpden", sep="\t", header=T)
snpdens2 = read.table("robeensis_chr2_12kb_window.snpden", sep="\t", header=T)
snpdens3 = read.table("robeensis_chr3_12kb_window.snpden", sep="\t", header=T)
snpdens4 = read.table("robeensis_chr4_12kb_window.snpden", sep="\t", header=T)
snpdens5 = read.table("robeensis_chr5_12kb_window.snpden", sep="\t", header=T)
snpdens6 = read.table("robeensis_chr6_12kb_window.snpden", sep="\t", header=T)
snpdens7 = read.table("robeensis_chr7_12kb_window.snpden", sep="\t", header=T)
snpdens8 = read.table("robeensis_chr8_12kb_window.snpden", sep="\t", header=T)
snpdens9 = read.table("robeensis_chr9_12kb_window.snpden", sep="\t", header=T)
snpdens10 = read.table("robeensis_chr10_12kb_window.snpden", sep="\t", header=T)
snpdens11 = read.table("robeensis_chr11_12kb_window.snpden", sep="\t", header=T)
snpdens12 = read.table("robeensis_chr12_12kb_window.snpden", sep="\t", header=T)

################################################################################
### 2. Obtain a theta map for each chromosome and export it
################################################################################
## divide by the harmonic number for the number of haploid indiv a = 20
a= vector()
for(i in 1:19){
  a[i]=1/i
}
a_n = sum(a) #3.54774

maps = list(snpdens1, snpdens2, snpdens3, snpdens4, snpdens5, snpdens6, 
            snpdens7, snpdens8, snpdens9, snpdens10, snpdens11, snpdens12)

## loop through each chromosome
for(i in 1:12){
  snpdens = maps[[i]]
  map = matrix(ncol=3, nrow=nrow(snpdens))
  map=as.data.frame(map)
  map[,1]= snpdens$BIN_START
  map[1:(nrow(snpdens)-1),2]= snpdens$BIN_START[-1]
  map[nrow(snpdens),2] = map[(nrow(snpdens)-1),2]+12000
  map[,3] = snpdens$SNP_COUNT/(12000*a_n)
  write.table(map, paste0("theta_map_chr", i, ".12kb_window.txt"), sep="\t")
  
  ## to get the mutation rate per pb per gen, divide by 4*Ne
  head(map)
  map$V3 = map$V3/(4*143601)
  
  ## export values
  write.table(map, paste0("theta_map_chr", i, ".12kb_window_per_bp_per_gen.txt"), sep="\t")
}


