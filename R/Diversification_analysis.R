################################################################################
###	       COLOR-DEPENDENT DIVERSIFICATION ANALYSIS IN ANURANS
################################################################################
### What this script does:
### 1. load and manipulate color and phyloegentic data
### 2. fit and compare models color state-dpeendent and independent models of
### diversification

################################################################################
###	1. Data check and manipulation
################################################################################
setwd("/path/to/directory")

## load color data
color = read.csv("color_data.csv", header=T)

## load phylogenetic data from Portik et al 2023
tree = read.tree("TreePL-Rooted_Anura_bestTree.tre")
length(tree$tip.label) # 5326 tips
tips_portik = tree$tip.label

## add the tips of the full tree on the data matrix to have them all match
tips = tree$tip.label
tips= data.frame(tips,NA)
colnames(color)
tip.color = merge(tips, color, by.x="tips", by.y="tips_portik", all.x=TRUE)
head(tip.color)
color=tip.color

## make the row names the species names
row.names(color) = color[,1]

### remove non-anurans
## get the node number for the most recent common ancestor 
## of Leiopelma hamiltoni and Staurois parvus
getMRCA(tree,c("Leiopelma_hamiltoni", "Staurois_parvus")) #5411
## extract only anurans from our tree
anura = extract.clade(tree, 5411, root.edge = 0, collapse.singles = TRUE,
                      interactive = FALSE)
length(anura$tip.label) # 5242
tree2 = anura

## re-order data to match tip labels
color<-color[tree2$tip.label,]
row.names(color) = tree2$tip.label

## get a column for green with levels "green", "no_green" and "green+no_green"
## for the diversification analysis, we have all species so we also have NAs in this column
color$greenbin = NA
for(i in 1:nrow(color)){
  if(is.na(color$green[i])==FALSE){
    if(color$green[i]==1){
      color$greenbin[i]="green"
    }
    if(color$green[i]==0){
      color$greenbin[i]="no_green"
    }
    if(color$green[i]==0.5){
      color$greenbin[i]="green+no_green"
    }
  }
}
color$greenbin

################################################################################
###	2. Fit color-dependent and color-independent models of diversification
################################################################################
## check our data
nrow(color) #5242
length(tree2$tip.label)#5242

## create a color vector
dipoly=color$greenbin
names(dipoly) = rownames(color)
dipoly = as.factor(dipoly)

## set the fractions of sampled species in each of our groups
## here  we don't know what we didn't sample
length(dipoly) # we have 5242 species sampled
## if we consider that we sampled randomly each category and there are 
## 7651 anuran species in total (according to https://amphibiansoftheworld.amnh.org/)
5242/7651 = 0.6851392
## we have sampled 0.6851392 of the total (known!) anuran diversity
## we count the number of taxa  sampled int he FULL phylogeny, not only the species with color data
## so we set that fraction for each category
rho<-setNames(c(0.6851392,0.6851392,0.6851392), 1:3) 

## fit MuSSE models for diversification implemented in the diversitree package
library(diversitree)
## make the tree ultrametric to be compatible with diversitree
is.ultrametric(tree2)
tree4 = force.ultrametric(tree2, method=c("nnls"))
is.ultrametric(tree4)
head(dipoly)

## change the levels of dipoly from 0,1,2 to 1,2,3 to be compatible with diversitree
## corresponding to: 1 = green 2 = green+no_green 3 = no_green 
levels(dipoly) <- c("1", "2", "3") 
head(dipoly)
dipoly = as.integer(dipoly)
names(dipoly) = row.names(color)

## create a general model incorporating our 3 color states
general.musse.model<-make.musse(tree=tree4, states=dipoly ,k=3,sampling.f=rho) 

## model MuSSE: to reduce the number of variables to estimate, we will simplify the model
## lets consider that our states are ordered: 
## green <-> green+no_green <-> no_green
# the transition green <-> no_green is not allowed
ordered.musse.model<-constrain(general.musse.model, q13~0,q31~0) 
argnames(ordered.musse.model) # 10 parameters 

## model ER.MuSSE: transition rates between color states are equal
orderedER.musse.model<-constrain(ordered.musse.model, q23~q12,q32~q12,q21~q12) 
argnames(orderedER.musse.model) # 7 parameters

## NULL model: the speciation and extrinction rates do not vary as function of color state
## transition rates between color states are different
ordered.musse.null<-constrain(ordered.musse.model, lambda3~lambda1,lambda2~lambda1, mu3~mu1,mu2~mu1)
argnames(ordered.musse.null) # 6 parameters

## ER NULL model: the speciation and extrinction rates do not vary as function of color state
## transition rates between color states are equal
orderedER.musse.null<-constrain(orderedER.musse.model, lambda3~lambda1,lambda2~lambda1, mu3~mu1,mu2~mu1) 
argnames(orderedER.musse.null)# 3 parameters

### fit the models
## get a starting point
p = starting.point.musse(tree4, k=3)

## fit ordered MuSSE model 
ordered.musse.mle<-find.mle(ordered.musse.model, x.init=p[argnames(ordered.musse.model)]) 
save(ordered.musse.mle, file="ordered.musse.mle_portik_completetree.Rdata")

## fit ordered ER MuSSE model 
orderedER.musse.mle<-find.mle(orderedER.musse.model, x.init=p[argnames(orderedER.musse.model)]) 
save(orderedER.musse.mle, file="orderedER.musse.mle_portik_completetree.Rdata")

## fit ordered null model 
ordered.musse.null.mle<-find.mle(ordered.musse.null, x.init=p[argnames(ordered.musse.null)]) 
save(ordered.musse.null.mle, file="ordered.musse.null.mle_portik_completetree.Rdata")

## fit ordered ER null model 
orderedER.musse.null.mle<-find.mle(orderedER.musse.null, x.init=p[argnames(orderedER.musse.null)]) 
save(orderedER.musse.null.mle, file="orderedER.musse.null.mle_portik_completetree.Rdata")

## compile the results in a table
musseAnova<-anova(orderedER.musse.null.mle, 
                  Null=ordered.musse.null.mle,
                  ER.MuSSE=orderedER.musse.mle,
                  MuSSE=ordered.musse.mle) 
musseAnova # the FULL model is significantly a better fit than the null model
write.csv(musseAnova, "musseAnova_portik_completetree.csv")

coefs_musse = coef(ordered.musse.mle)
coefs_musse
write.csv(coefs_musse, "coefs_musse_portik_fulltree.csv")

## likelihood ratio tests
library(lmtest)
# likelihood ratio test comparing ER and SYM models
lrtest(orderedER.musse.null.mle, ordered.musse.mle)
lrtest(orderedER.musse.null.mle, orderedER.musse.mle)
lrtest(ordered.musse.null.mle, orderedER.musse.mle)
lrtest(ordered.musse.null.mle, ordered.musse.mle)

## check the number of species in each category
sum(color$greenbin=="green", na.rm=T) #  186
sum(color$greenbin=="green+no_green", na.rm=T) # 191
sum(color$greenbin=="no_green", na.rm=T) # 1952
sum(is.na(color$greenbin)) # 2913

