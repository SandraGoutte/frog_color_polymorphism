################################################################################
###     JOINT EVOLUTION OF GREEN COLOR AND HABITAT PREFERENCE IN ANURANS
################################################################################
### What this script does:
### 1. fits can compare Mk models of evolution for habitat preferences 
### in anurans
### 2. reconstruct the ancestral states of habitat preferences in anurans
### 3. Fit and compare models of habitat-green color joint-evolution 

################################################################################
### 1. Fit Mk models for habitat preferences
################################################################################
### prune tree removing the species for which we dont have habitat data
## Check how many species (= tree tips) are in the tree
length(tree2$tip.label) #2363
colnames(color)

## remove unchecked species
remove = row.names(color[which(is.na(color$Habitat)),])
length(remove) # 374 species to remove

## prune the tree accordingly
tree3<- drop.tip(tree2,remove)
length(tree3$tip.label) # 1989 species included in the analysis
  
## re-order data to match tip labels
color2<-color[tree3$tip.label,]
nrow(color2) #1989
row.names(color2)
row.names(color2) = tree3$tip.label

## create a habitat vector
color2$Habitat = as.factor(color2$Habitat)
levels(color2$Habitat) #"A"   "Aq"  "AT"  "FT"  "T"   "TAq"
habitats = color2$Habitat
habitats = setNames(color2$Habitat, rownames(color2))
length(habitats)

### fit different Mk models and choose the best fitting one
## model ER: a single transition rate between all habitats
job::job({
    habitatsER_fitMk<-fitMk(tree3,habitats,model="ER", quiet=T, pi="fitzjohn")
    save(habitatsER_fitMk,  file="habitatsER_fitMk.Rdata")
})

## model SYM: symetrical transition rates 
job::job({
    habitatsSYM <- fitDiscrete(tree3,habitats,ordered=F,model="SYM", 
                               quiet=T, pi="fitzjohn")
    save(habitatsSYM,  file="habitatsSYM.Rdata")
})
  
## model ARD: all rate different
job::job({
    habitatsARD <- fitDiscrete(tree3,habitats,ordered=F,model="ARD", 
                               quiet=T, pi="fitzjohn")
    save(habitatsARD,  file="habitatsARD.Rdata")
})

## design a custom Ordered model where only transitions between adjacent
## habitats are allowed:  FT <-> T <-> AT <-> A and T <-> TAq <-> Aq
# create a design matrix
ordered.model <- matrix(c(
    0,0,1,0,0,0,
    0,0,0,0,0,2,
    3,0,0,0,4,0,
    0,0,0,0,5,0,
    0,0,6,7,0,8,
    0,9,0,0,10,0),6,6,byrow=T, dimnames=list(levels(habitats),
                                             levels(habitats)))
ordered.model
  
## fit the custom ordred model
job::job({
    habitatsOrdered <- fitDiscrete(tree3,habitats,model=ordered.model, 
                                   quiet=T, pi="fitzjohn", supressWarnings=T)
    save(habitatsOrdered,  file="habitatsOrdered.Rdata")
})

## custom Reduced model 1: based on the ARD model: only transitions >= 0.001
## create a design matrix
reduced.model1 <- matrix(c(
    0,0,3,0,0,0,
    0,0,0,0,0,9,
    1,0,0,0,6,10,
    0,13,12,0,7,14,
    0,0,4,5,0,11,
    0,2,0,0,8,0),6,6,byrow=T, dimnames=list(levels(habitats), levels(habitats)))

## fit the custom model
job::job({
    habitatsReduced1 <- fitDiscrete(tree3,habitats,model=reduced.model1, 
                                    quiet=T, pi="fitzjohn", supressWarnings=T)
    save(habitatsReduced1,  file="habitatsReduced1.Rdata")
})

## custom Reduced model 2: based on the ARD model: only transitions > 0.001
## create a design matrix
reduced.model2 <- matrix(c(
    0,0,3,0,0,0,
    0,0,0,0,0,9,
    1,0,0,0,6,10,
    0,0,0,0,7,11,
    0,0,4,5,0,12,
    0,2,0,0,8,0),6,6,byrow=T, dimnames=list(levels(habitats), levels(habitats)))
  ## fit the custom model
  job::job({
    habitatsReduced2 <- fitDiscrete(tree3,habitats,model=reduced.model2, 
                                    quiet=T, pi="fitzjohn", supressWarnings=T)
    save(habitatsReduced2,  file="habitatsReduced2.Rdata")
  })
  
#### Compare the different models of evolution
aic <- setNames(c(AIC(habitatsER), AIC(habitatsSYM),AIC(habitatsOrdered),
                  AIC(habitatsReduced1), AIC(habitatsReduced2),  
                  AIC(habitatsARD)), c("ER","SYM", "Ordered","Reduced1",
                                       "Reduced2", "ARD"))
aic ## habitatsReduced2 has the lowest AIC

## compute Akaike weights = weight of evidence in support of the model
aic.w(aic)
  
## put everything together in a table
habitat_models_compare = round(data.frame(k=c(habitatsER$opt$k,habitatsSYM$opt$k,habitatsOrdered$opt$k,
                                              habitatsReduced1$opt$k,habitatsReduced2$opt$k,habitatsARD$opt$k),
                                            logL=c(logLik(habitatsER),logLik(habitatsSYM),logLik(habitatsOrdered),
                                                   logLik(habitatsReduced1),logLik(habitatsReduced1),logLik(habitatsARD)),
                                            AIC=aic, Akaike.w=as.vector(aic.w(aic))),3)
write.csv(habitat_models_compare, "habitat_models_compare.csv")
  
## plot the different models
par(mfrow=c(2,3), mar=c(1,1,1,1))
plot(habitatsER) 
#plot(habitatsER, width=TRUE,color=TRUE,mar=rep(0.1,4), cex=2,show.zeros=FALSE)
mtext(paste("a) unordered ER, AIC:",
              round(AIC(habitatsER),2)),adj=0, cex=.5)
plot(habitatsSYM)
mtext(paste("c) unordered SYM, AIC:",
              round(AIC(habitatsSYM),2)),adj=0,cex=.5)
plot(habitatsARD)
mtext(paste("c) unordered ARD, AIC:",
              round(AIC(habitatsARD),2)),adj=0, cex=.5)
plot(habitatsOrdered, show.zeros=FALSE)
mtext(paste("c) unordered POLY, AIC:",
              round(AIC(habitatsOrdered),2)),adj=0, cex=.5)
plot(habitatsReduced1, show.zeros=FALSE)
mtext(paste("c) unordered POLY, AIC:",
              round(AIC(habitatsReduced1),2)),adj=0, cex=.5)
plot(habitatsReduced2, show.zeros=FALSE)
mtext(paste("c) unordered POLY, AIC:",
              round(AIC(habitatsReduced2),2)),adj=0, cex=.5)

################################################################################
### 2. Reconstruct ancestral state of habitat using reduced.model2
################################################################################
## create stochastic maps
  job::job({
    mtree_habitats<-make.simmap(tree3, habitats, model=reduced.model2, nsim=100,
                                pi="fitzjohn")
    save(mtree_habitats,  file="mtree_habitats.Rdata")
  })
 
## plot
job::job({
  pdf("ASR_habitats_reduced.model2.pdf", width=25, height=25)
  cols<-setNames(c("green","cyan", "darkgreen","black", "brown", "steelblue"),
                  c("A", "Aq", "AT", "FT", "T", "TAq"))
  ## plot an outline tree:
  plotTree(mtree_habitats[[1]],ftype="off",lwd=1.1, type="fan",color="white", 
           offset=0, show.tip.label = F) 
  par(fg="transparent",lend=1)
  ## now plot our 100 stochastic map trees with 99% transparency
  for(i in 1:100) plotSimmap(mtree_habitats[[i]],
                              colors=sapply(cols,make.transparent,alpha=0.02),
                             add=TRUE,lwd=1.2,ftype="off",offset=0,  type="fan")

    par(fg="black")
    habs = c("arboreal",  "aquatic" , "arboreal/terrestrial" ,
             "fossorial/terrestrial" , "terrestrial" ,  "terrestrial/aquatic")
    dev.off()
})
  
## create a consensus from the 100 stochastic maps
sum = summary(mtree_habitats)
allStates<-rbind(habitatsReduced2$data[tree3$tip.label,],sum$ace)
map<-reorder(tree3,"cladewise")
map<-paintSubTree(map,node=Ntip(map)+1,
                    state=names(which(sum$ace[1,]==max(sum$ace[1,]))))
for(i in 1:nrow(map$edge)){
    states<-sapply(map$edge[i,],function(x,aa) 
      names(which(aa[x,]==max(aa[x,]))),
      aa=allStates)
    if(length(unique(states))!=1)
      map<-paintSubTree(map,node=map$edge[i,2],state=states[2],
                        stem=0.5)
}

#### plot and compare to our plotted 100 Simmap 
job::job({
    pdf("ASR_habitats_reduced2_consensus_nolegends.pdf", width=25, height=25)
    cols<-setNames(c("green","cyan", "darkgreen","black", "brown", "steelblue"),
                   c("A", "Aq", "AT", "FT", "T", "TAq"))
    ## plot our consensus painted tree
    plotSimmap(map,colors=cols,lwd=1.2,ftype="off",offset=0,  type="fan")
    nodelabels(pie=nodes,piecol=cols,cex=0.2)
    par(fg="black")
    habs = c("arboreal",  "aquatic" , "arboreal/terrestrial" ,
             "fossorial/terrestrial" , "terrestrial" ,  "terrestrial/aquatic")
    dev.off()
  })
  mapped.habitats = map
  save(mapped.habitats, file="mapped.habitats.Rdata")
}
## we saved the consensus of the mapped tree for 6 habitats
## we will use this mapped consensus for the downstream analyses!

################################################################################
### 3. Test models of habitat-green color joint-evolution 
################################################################################
## create a variable with habitat and color state combined
habitats = color2$Habitat
greens = color2$greenbin
names(greens)=row.names(color2)
xy<-setNames(interaction(greens,habitats),names(greens))
head(xy)
levels(xy)<-gsub(".","+",levels(xy),fixed=TRUE)
head(xy)
# we now have every possible combination of habitat+color state as states 
# (3*6=18 combinations in total)

## To create the different models of evolution, we populate a 18*18 matrix
## with 0 is the transition is not allowed, and a coefficient if it is
## create the transition matrix
MODEL<-matrix(0,18,18,dimnames=list(levels(xy),levels(xy)))
MODEL
write.csv(MODEL, file="MODEL.csv")
## edit the matrix manually

## FULL model: every possible transition rate can take a different value.
## to limit the number of parameters, direct transitions between "green" 
## and "no_green" are not allowed
## only habitat transitions allowed in the habitat reduced model 2 are allowed
## simultaneous transitions between habitats and color states are not allowed
FULL = read.csv("MODEL2_model4_strict.csv")
row.names(FULL)=FULL[,1]
FULL = FULL[,-1]
FULL=as.matrix(FULL)
colnames(FULL)=rownames(FULL)
FULL

## model 2: ARD for green and ARD for habitat transitions
## but habitat transitions are independent of color state
## and color transitions are independent of habitat
model2 = read.csv("MODEL4_model4_strict.csv")
row.names(model2)=model2[,1]
model2 = model2[,-1]
model2=as.matrix(model2)
colnames(model2)=rownames(model2)

## model 3: ARD for color state transitions in each habiats, 
## but the transitions between habiatats are independent of color state
model3 = read.csv("MODEL5_model4_strict.csv")
row.names(model3)=model3[,1]
model3 = model3[,-1]
model3=as.matrix(model3)
colnames(model3)=rownames(model3)

## model 4: ARD for color state transition but independently of habitat 
## habitat transitions are ARD and color state-dependent 
model4 = read.csv("model6_model4_strict.csv")
row.names(model4)=model4[,1]
model4 = model4[,-1]
model4=as.matrix(model4)
colnames(model4)=rownames(model4)

## fit the models to our data
job::job({
  fit_xy_FULL<-fitMk(tree3,xy,model=FULL)
  save(fit_xy_FULL, file="fit_xy_FULL.Rdata")
})

job::job({
  fit_xy_model2<-fitMk(tree3,xy,model=model2)
  save(fit_xy_model2, file="fit_xy_model2.Rdata")
})

job::job({
  fit_xy_model3<-fitMk(tree3,xy,model=model3)
  save(fit_xy_model3, file="fit_xy_model3.Rdata")
})

job::job({
  fit_xy_model4<-fitMk(tree3,xy,model=model4)
  save(fit_xy_model4, file="fit_xy_model4.Rdata")
})

## plot the fitted models
par(mfrow=c(2,2))
plot(fit_xy_FULL,width=FALSE, lwd=2,color=TRUE,show.zeros=FALSE,
     mar=rep(0,4),spacer=0, tol=1e-2)
plot(fit_xy_model2,width=FALSE, lwd=2,color=TRUE,show.zeros=FALSE,
     mar=rep(0,4),spacer=0, tol=1e-2)
plot(fit_xy_model3,width=TRUE,color=TRUE,show.zeros=FALSE,
     mar=rep(0,4),spacer=0, tol=2e-2)
plot(fit_xy_model4,width=TRUE,color=TRUE,show.zeros=FALSE,
     mar=rep(0,4),spacer=0, tol=1e-2)

## compare the models using anova
green.anova = anova(fit_xy_FULL, fit_xy_model2,fit_xy_model3,fit_xy_model4)
write.csv(green.anova, "compare_green_habitats_models.csv")

## compare the different models of evolution with likelihood ratio tests
library(lmtest)
lrtest(fit_xy_model2, fit_xy_FULL)
lrtest(fit_xy_model3, fit_xy_FULL)
lrtest(fit_xy_model4, fit_xy_FULL)
