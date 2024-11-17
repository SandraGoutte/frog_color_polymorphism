################################################################################
###	                 PLOT OUTPUT OF TWISST ANALYSIS 
################################################################################
### What this script does: 
### 1. import output of TWISST analysis
### 2. Identify tree topologies where green individuals group together
### 3. Sum the weights of green-grouping topologies
### 4. Plot the result to identify the limits of our region of interest 

## get plot_twisst R functions from https://github.com/simonhmartin/twisst/blob/master/plot_twisst.R
source("plot_twisst.R")

################################################################################
### 1. Import output from the TWISST analysis
################################################################################
## TWISST analysis was run on 4 species, 2 phenotypes (8 groups in total) 
## on 50 SNPs non-overlapping windows. This was done on 50 VCFs starting 
## at 1 SNP interval in order to obtain overlapping sliding windows

## list the tree files to be imported
window_data_files= vector()
window_data_files[1]='renamed_phased_ptychadena_4sp_chr7_w50SNPs_0.phyml_bionj.data.tsv'
for(i in 1:49){
  window_data_files = c(window_data_files, paste0("renamed_phased_ptychadena_4sp_chr7_w50SNPs_",i,".phyml_bionj.data.tsv"))
}

## list the weigth files to be imported
weights_files = vector()
weights_files[1] = "renamed_phased_ptychadena_4sp_chr7_w50SNPs_0.phyml_bionj.weights.tsv"
for(i in 1:49){
  weights_files = c(weights_files, paste0("renamed_phased_ptychadena_4sp_chr7_w50SNPs_",i,".phyml_bionj.weights.tsv"))
}

## list the vcfs (we have 50 vcfs in total)
names=vector()
names[1] = "vcf_0"
for(i in 1:49){
  names = c(names, paste0("vcf_",i))
}
names 

## load the data
twisst_data <- import.twisst(weights_files, window_data_files)

################################################################################
### 2. Identify topologies where green individuals form a monophyletic clade
################################################################################
## group tips of interest: the green individuals of the 4 species
greens <- c("robeensis_green", "nana_green", "levenorum_green", "erlangeri_green")

## function to check if the 4 tips are grouped together in a given tree
is_grouped <- function(tree, tips) {
  return(is.monophyletic(tree, tips))
}

## initialize lists to store grouped and non-grouped trees
grouped_trees <- list()
non_grouped_trees <- list()

# iterate through each tree and classify it
index_grouped=vector()
index_non_grouped=vector()
for (i in seq_along(twisst_data$topos)) {
  tree <- twisst_data$topos[[i]]
  if (is_grouped(tree, greens)) {
    grouped_trees[[length(grouped_trees) + 1]] <- tree
    index_grouped = c(index_grouped, i)
  } else {
    non_grouped_trees[[length(non_grouped_trees) + 1]] <- tree
    index_non_grouped = c(index_non_grouped, i)
  }
}
length(index_grouped) #225
length(index_non_grouped) #10170

## function to color tree tips
color_tips <- function(tree, tips_of_interest) {
  # default color for tips
  tip_colors <- rep("black", length(tree$tip.label))
  # assign color to tips of interest
  tip_colors[tree$tip.label %in% tips_of_interest] <- "forestgreen"
  return(tip_colors)
}

## plot a topology where greens are grouped and one where they are not
par(mfrow=c(1,2))
colo = color_tips(grouped_trees[[1]], greens)
plot(grouped_trees[[1]],tip.color= colo)
colo = color_tips(non_grouped_trees[[1]], greens)
plot(non_grouped_trees[[1]],tip.color= colo)

################################################################################
### 3. Sum the weights of green-grouping topologies and non green-grouping 
###    topologies
################################################################################
twisst_data$weights_grpall = twisst_data$weights
for(j in 1:length(twisst_data$weights_grpall)){
  for(i in 1:nrow(twisst_data$weights_grpall[[j]])){
    twisst_data$weights_grpall[[j]]$topo1201[i] = sum(twisst_data$weights_grpall[[j]][index_grouped][i,])
  }
  for(i in 1:nrow(twisst_data$weights_grpall[[j]])){
    twisst_data$weights_grpall[[j]]$topo1[i] = sum(twisst_data$weights_grpall[[j]][index_non_grouped][i,])
  }
}

for(j in 1:length(twisst_data$weights_grpall)){
  twisst_data$weights_grp2[[j]] = matrix()
  twisst_data$weights_grp2[[j]] = as.data.frame(twisst_data$weights_grp2[[j]])
  twisst_data$weights_grp2[[j]] = cbind(twisst_data$weights_grpall[[j]]$topo1, twisst_data$weights_grpall[[j]]$topo1201)
}

################################################################################
### 4. Plot the result to identify the limits of our region of interest 
################################################################################
## set up the colors for our plot
newcols = c( "#E0E0DE05","#2BCE4825")
{
par(mfrow=c(2,1), mar=c(2,4,6,4),xaxs="i")
plot.weights(weights_dataframe=twisst_data$weights_grp2[[1]], positions=twisst_data$window_data[[1]][,c("start","end")],
              fill_cols=newcols, line_cols=NULL, stacked=F, xlim=c(26300000,26500000))
for(i in 1:49){
  plot.weights(weights_dataframe=twisst_data$weights_grp2[[i]], positions=twisst_data$window_data[[i]][,c("start","end")],
                fill_cols=newcols, line_cols=NULL, stacked=F, xlim=c(26300000,26500000), add=T)
  }
## add all the snps to see the distribution of snps
depth = read.table("phased_ptychadena_4sp_chr7_26000000_27000000.ldepth", header=T)
segments(x0 = depth$POS, y0 = -1, x1 = depth$POS, y1 = 0, col = "darkgrey", lwd = 2)
## add the SNPs that are polymorphic across the 4 species
snps = read.csv("polys_all4.csv")
segments(x0 = snps$pos, y0 = -1, x1 = snps$pos, y1 = 0, col = "red", lwd = 2)
  
## zoom on the region of interest
par(mar=c(6,4,2,4))
plot.weights(weights_dataframe=twisst_data$weights_grp2[[1]], positions=twisst_data$window_data[[1]][,c("start","end")],
            fill_cols=newcols, line_cols=NULL, stacked=F, xlim=c(26412000,26420000), ylim=c(0,.4))
for(i in 1:49){
  plot.weights(weights_dataframe=twisst_data$weights_grp2[[i]], positions=twisst_data$window_data[[i]][,c("start","end")],
              fill_cols=newcols, line_cols=NULL, stacked=F, xlim=c(26412000,26420000), add=T, ylim=c(0,.4))
}
## add boundaries for our window of interest
abline(h=0.05, lty=2, col="black")
abline(v=26414430, lty=2, col="red")
abline(v=26417815, lty=2, col="red")
## add all the snps to see the distribution of snps
depth = read.table("phased_ptychadena_4sp_chr7_26000000_27000000.ldepth", header=T)
segments(x0 = depth$POS, y0 = -1, x1 = depth$POS, y1 = 0, col = "darkgrey", lwd = 2)
## add the SNPs that are polymorphic across the 4 species
snps = read.csv("polys_all4.csv")
segments(x0 = snps$pos, y0 = -1, x1 = snps$pos, y1 = 0, col = "red", lwd = 2)
}

