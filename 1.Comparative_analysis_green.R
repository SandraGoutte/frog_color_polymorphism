#############################################################################
###	COMPARATIVE ANALYSES OF GREEN COLOR POLYMORPHISM IN ANURANS
################################################################################
### What this script does:
### 1. reads and filter color and phylogenetic data
### 2. fits and compares different Mk nodels of evolution for the green dorsal 
### coloration in anurans  
### 3. reconstruct the ancestral states of the green coloration in anurans
### 4. plots the number of changes between color states
### 5. plot the evolutionary time spent in each color state


### install packages
{
  install.packages("job")
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("ggtree")
  
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("ggtreeExtra")
  
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("treeio")
}
### Load the libraries necessary for the analysis
{
  library(phytools)
  packageVersion("phytools")
  library(ape)
  library(geiger)
  library(job)
  library(ggimage)
  library(ggtree)
  library(TDbook)
  library(ggtreeExtra)
  library(treeio)
  library(tidytree)
  library(ggstar)
  library(ggplot2)
  library(ggnewscale)
  library(caper)
  library(phangorn)
  library(corHMM)
  library(dplyr)
  library(lmtest)
  library(ggplot2)
  library(tidyverse)
}
### Set working directory 
setwd("/path/to/directory")

################################################################################
### 1. Data import and manipulation
################################################################################
## Load the phenotype data (.csv file)
color = read.csv("color_data.csv", header=T)

## Load the phylogenetic tree from Portik et al. 2023
tree = read.tree("TreePL-Rooted_Anura_bestTree.tre")
length(tree$tip.label) # 5326 anurans
tips_portik = tree$tip.label

### Data check and manipulation
## Check that all columns containing numbers are in class "numeric"
class<-vector()
for(i in 1:ncol(color)){
  class[i]<-class(color[,i])
}
class<-as.data.frame(class)
row.names(class)<-colnames(color)
class

## add the tips of the full tree on the data matrix to have them all match
tips = tree$tip.label
tips= data.frame(tips,NA)
colnames(color)
tip.color = merge(tips, color, by.x="tips", by.y="tips_portik", all.x=TRUE)
head(tip.color)
color=tip.color

## make the row names the species names
colnames(color)
row.names(color) = color[,1]

### remove unchecked species
## check the structure of the tree data
str(tree)

## check how many species (= tree tips) are in the tree
length(tree$tip.label) # 5326

## remove unchecked species
remove = row.names(color[which(is.na(color$nb_frogs_checked_color)),])
tree2 = anura
length(remove) # 2963 species to remove
tree2<- drop.tip(tree,remove)
length(tree2$tip.label) # 2363 species included in the analysis

## re-order data to match tip labels
color<-color[tree2$tip.label,]
row.names(color) = tree2$tip.label
length(tree2$tip.label) # 2363

## some checked species have NAs in the color columns. Let's change them to 0
colnames(color)[19:31] # are the color columns. Should be only 0, 0.5, 1
color[, 19:31][is.na(color[, 19:31])] <- 0

## create a "green polymorphism" column
color$greenbin = NA
for(i in 1:nrow(color)){
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
color$greenbin

## count number of green species:
## rows where greenbin == "green" OR greenbin == "green+no_greenown"
count_rows <- color %>%
  filter(greenbin == "green+no_green" | greenbin == "green" ) %>%
  nrow()

print(count_rows)
print(count_rows)/nrow(color)

################################################################################
###	2. fit Mk models to the green color data 
################################################################################
## create a vector with the color data
green = color$greenbin
green = setNames(color$greenbin, rownames(color))
green = as.factor(green)

## model ER: equal transition rates for all
job::job({
    greenER_Mk = fitMk(tree2, green, model="ER", quiet=T, pi="fitzjohn")
    save(greenER_Mk, file="greenER_Mk.Rdata")
})
  
## model ARD: all rates different model
job::job({
    greenARD_Mk<-fitMk(tree2,green,model="ARD", quiet=T, pi="fitzjohn")
    save(greenARD_Mk,  file="greenARD_Mk.Rdata")
})

## model SYM: symetrical
job::job({
    greenSYM_Mk <- fitMk(tree2,green,ordered=F,model="SYM", quiet=T, pi="fitzjohn")
    save(greenSYM_Mk,  file="greenSYM_Mk.Rdata")
})
  
## model POLY: must go through the "green+no_green" state to transit 
## from "green" to "no_green" and vice-versa
Q<-matrix(0,3,3)
rownames(Q)<-colnames(Q)<-c("green","green+no_green","no_green")
Q[1,]<-c(0,1,0)
Q[2,]<-c(2,0,3)
Q[3,]<-c(0,4,0)
#diag(Q)<--rowSums(Q)
Q
job::job({
  green_custom_poly = fitMk(tree2,green , model=Q, quiet=T, pi="fitzjohn", ordered=TRUE)
  save(green_custom_poly, file="green_custom_poly.Rdata")
})

## plot the results
par(mfrow=c(2,2))
plot(greenER_Mk, width=TRUE,color=TRUE,mar=rep(0.1,4), cex=2,show.zeros=FALSE)
mtext(paste("a) unordered ER, AIC:",
              round(AIC(greenER_Mk),2)),adj=0)
plot(greenARD_Mk, width=TRUE,color=TRUE,mar=rep(0.1,4), cex=2,show.zeros=FALSE)
mtext(paste("c) unordered ARD, AIC:",
              round(AIC(greenARD_Mk),2)),adj=0)
plot(greenSYM_Mk, width=TRUE,color=TRUE,mar=rep(0.1,4), cex=2,show.zeros=FALSE)
mtext(paste("c) unordered SYM, AIC:",
              round(AIC(greenSYM_Mk),2)),adj=0)
plot(green_custom_poly, width=TRUE,color=TRUE,mar=rep(0.1,4), cex=2,show.zeros=FALSE)
mtext(paste("c) unordered POLY, AIC:",
              round(AIC(green_custom_poly),2)),adj=0)

## compare the fit of the models using anova
green.anova = anova(greenER_Mk,greenSYM_Mk,greenARD_Mk,green_custom_poly)
write.csv(green.anova, "green.anova_4models_fitMk_portik.csv")

## function for likelihood ratio test
lrtest<-function(model1,model2){
    lik1<-logLik(model1)
    lik2<-logLik(model2)
    LR<--2*(lik1-lik2)
    #degf<-attr(lik2,"df")-attr(lik1,"df")
    degf<-length(model2$rates)-length(model1$rates)
    P<-pchisq(LR,df=degf,lower.tail=FALSE)
    cat(paste("Likelihood ratio = ",
              signif(LR,5),"(df=",degf,") P = ",
              signif(P,4),"\n",sep=""))
  invisible(list(likelihood.ratio=LR,p=P))
}

# likelihood ratio test comparing pairs of models
lrtest(green_custom_poly,greenARD_Mk) # Likelihood ratio = 4.2935(df=2) P = 0.1169
lrtest(greenSYM_Mk,green_custom_poly) #Likelihood ratio = 221.8(df=1) P = 3.659e-50
lrtest(greenER_Mk,green_custom_poly) #Likelihood ratio = 338.01(df=3) P = 5.894e-73

################################################################################
### 3. Ancestral state reconstruction of color states
################################################################################
## use make simmap on the fitted model
job::job({
    green.smap.trees_Mk_1000<-make.simmap(tree2,greenARD_Mk$data,
                                     Q=as.Qmatrix(greenARD_Mk),pi=greenARD_Mk$pi,
                                     nsim=1000)
    save(green.smap.trees_Mk_1000, file="green.smap.trees_Mk_1000.Rdata")
})

## plot
job::job({
  green.cols<-setNames(c("darkgreen","lightgreen", "lightgrey"),c("green", "green+no_green", "no_green"))
  pdf("ASR_green_Mk_1000.pdf", width=25, height=25)
  ##plot an outline tree:
  par(fg="transparent",lend=1)
  plotTree(tree2,ftype="off",lwd=1, type="fan",color="white", offset=0) 
  ## now plot our 1000 stochastic map trees with 99% transparency
  for(i in 1:1000) plot(green.smap.trees_Mk_1000[[i]],
                        colors=sapply(green.cols,make.transparent,alpha=0.3),
                        add=TRUE,lwd=4,ftype="off",offset=0,  type="fan")
    par(fg="black")
    habs = c("green", "green+no_green", "no_green")
    legend(x="bottomleft",habs,pch=22, pt.bg=green.cols,pt.cex=1.5,bty="n",cex=3)
    dev.off()
})
  
green.smap.trees
green.obj<-summary(green.smap.trees)
green.smap.trees_Mk
green.obj_Mk<-summary(green.smap.trees_Mk)
  
## have a look at the trees
plot(green.smap.trees,ftype="off",lwd=1,colors=green.cols,type="fan")
  
## get posterior probabilities at the nodes
par(fg="transparent")
plot(summary(green.smap.trees),colors=green.piecol,type="fan",ftype="off",
       lwd=1,cex=c(0.5,0.3))
par(fg="black")
legend(x="topleft",legend=colnames(green.asr),
         pt.cex=1.8,pch=16,cex=0.8,col=green.cols,
         bty="n")
  

################################################################################
### 4. plot the number of color state changes
################################################################################
## plot change Map function
plot.changesMap<-function(x,...){
  if(hasArg(bty)) bty<-list(...)$bty
  else bty<-"l"
  if(hasArg(alpha)) alpha<-list(...)$alpha
  else alpha<-0.3
  if(hasArg(xlim)) xlim<-list(...)$xlim
  else xlim<-NULL
  if(hasArg(ylim)) ylim<-list(...)$ylim
  else ylim<-NULL
  if(hasArg(main)) main<-list(...)$main
  else main<-NULL
  if(hasArg(colors)){ 
    colors<-list(...)$colors
    nn<-names(colors)
    colors<-setNames(make.transparent(colors,alpha),nn)
  } else { 
    colors<-if(length(x$trans)==2) 
      setNames(make.transparent(c("blue","red"),alpha),x$trans)
    else
      setNames(rep(make.transparent("blue",alpha),length(x$trans)),
               x$trans)
  }
  if(length(colors)<length(x$trans))
    colors<-rep(colors,ceiling(length(x$trans)/length(colors)))[1:length(x$trans)]
  if(is.null(names(colors))) colors<-setNames(colors,x$trans)
  if(hasArg(transition)){ 
    transition<-list(...)$transition
    if(length(transition)>1){
      cat("transition should be of length 1; truncating to first element.\n")
      transition<-transition[1]
    }
  } else transition<-NULL
  p<-x$p
  hpd<-x$hpd
  bw<-x$bw
  if(length(x$trans)==2&&is.null(transition)){
    plot(p[[1]]$mids,p[[1]]$density,xlim=if(is.null(xlim)) 
      c(min(x$mins)-1,max(x$maxs)+1) else xlim,
      ylim=if(is.null(ylim)) c(0,1.2*max(c(p[[1]]$density,
                                           p[[2]]$density))) else ylim,
      type="n",xlab="number of changes",
      ylab="relative frequency across stochastic maps",
      bty=bty)
    y2<-rep(p[[1]]$density,each=2)
    y2<-y2[-length(y2)]
    x2<-rep(p[[1]]$mids-bw/2,each=2)[-1]
    x3<-c(min(x2),x2,max(x2))
    y3<-c(0,y2,0)
    polygon(x3,y3,col=colors[x$trans[1]],border=FALSE)
    lines(p[[1]]$mids-bw/2,p[[1]]$density,type="s")
    y2<-rep(p[[2]]$density,each=2)
    y2<-y2[-length(y2)]
    x2<-rep(p[[2]]$mids-bw/2,each=2)[-1]
    x3<-c(min(x2),x2,max(x2))
    y3<-c(0,y2,0)
    polygon(x3,y3,col=colors[x$trans[2]],border=FALSE)
    lines(p[[2]]$mids-bw/2,p[[2]]$density,type="s")
    dd<-0.01*diff(par()$usr[3:4])
    lines(hpd[[1]],rep(max(p[[1]]$density)+dd,2))
    lines(rep(hpd[[1]][1],2),c(max(p[[1]]$density)+dd,
                               max(p[[1]]$density)+dd-0.005))
    lines(rep(hpd[[1]][2],2),c(max(p[[1]]$density)+dd,
                               max(p[[1]]$density)+dd-0.005))
    CHARS<-strsplit(x$trans[1],"->")[[1]]
    CHARS[1]<-paste("HPD(",CHARS[1],collapse="")
    CHARS[2]<-paste(CHARS[2],")",collapse="")
    T1<-bquote(.(CHARS[1])%->%.(CHARS[2]))
    text(mean(hpd[[1]]),max(p[[1]]$density)+dd,
         T1,pos=3)
    lines(hpd[[2]],rep(max(p[[2]]$density)+dd,2))
    lines(rep(hpd[[2]][1],2),c(max(p[[2]]$density)+dd,
                               max(p[[2]]$density)+dd-0.005))
    lines(rep(hpd[[2]][2],2),c(max(p[[2]]$density)+dd,
                               max(p[[2]]$density)+dd-0.005))
    CHARS<-strsplit(x$trans[2],"->")[[1]]
    CHARS[1]<-paste("HPD(",CHARS[1],collapse="")
    CHARS[2]<-paste(CHARS[2],")",collapse="")
    T2<-bquote(.(CHARS[1])%->%.(CHARS[2]))
    text(mean(hpd[[2]]),max(p[[2]]$density)+dd,
         T2,pos=3)
    CHARS<-strsplit(x$trans[1],"->")[[1]]
    T1<-bquote(.(CHARS[1])%->%.(CHARS[2]))
    CHARS<-strsplit(x$trans[2],"->")[[1]]
    T2<-bquote(.(CHARS[1])%->%.(CHARS[2]))
    legend("topleft",legend=c(T1,T2),pch=22,pt.cex=2.2,bty="n",
           pt.bg=colors[x$trans])
  } else {
    k<-if(is.null(transition)) length(x$states) else 1
    if(k>1) par(mfrow=c(k,k))
    ii<-if(is.null(transition)) 1 else which(x$trans==transition)
    max.d<-max(unlist(lapply(p,function(x) x$density)))
    for(i in 1:k){
      for(j in 1:k){
        if(i==j&&is.null(transition)) plot.new()
        else {
          CHARS<-strsplit(x$trans[ii],"->")[[1]]
          MAIN<-if(is.null(main)) bquote(.(CHARS[1])%->%.(CHARS[2])) else
            main
          plot(p[[ii]]$mids,p[[ii]]$density,xlim=if(is.null(xlim)) 
            c(min(x$mins)-1,max(x$maxs)+1) else xlim,
            ylim=if(is.null(ylim)) c(0,1.2*max.d) else ylim,
            type="n",xlab="number of changes",
            ylab="relative frequency",main=MAIN,font.main=1,
            bty=bty)
          y2<-rep(p[[ii]]$density,each=2)
          y2<-y2[-length(y2)]
          x2<-rep(p[[ii]]$mids-bw/2,each=2)[-1]
          x3<-c(min(x2),x2,max(x2))
          y3<-c(0,y2,0)
          polygon(x3,y3,col=colors[x$trans[ii]],border=FALSE)
          lines(p[[ii]]$mids-bw/2,p[[ii]]$density,type="s")
          dd<-0.03*diff(par()$usr[3:4])
          lines(hpd[[ii]],rep(max(p[[ii]]$density)+dd,2))
          text(mean(hpd[[ii]]),max(p[[ii]]$density)+dd,"HPD",pos=3)
          ii<-ii+1
        }
      }
    }
  }
}

job::job({
  density_consensus_Mk_1000<-density(green.smap.trees_Mk_1000)
  save(density_consensus_Mk_1000, file="density_consensus_green_1000_Mk.Rdata")
})

## get the values, including the mean time spent in each state
sum = summary(green.smap.trees_Mk_1000)

## count the number of changes 
changes = countSimmap(green.smap.trees_Mk_1000,states=names(green.cols))

## plot the number of color state transitions
par(mfrow=c(3,1))
plot.changesMap(density_consensus_Mk_1000, transition="green->green+no_green",ylim = c(0, 0.08), colors="lightgrey", cex.axes=4)
par(new=T)
plot.changesMap(density_consensus_Mk_1000, transition="green+no_green->green",ylim = c(0, 0.08), colors="darkgreen", cex=2)

plot.changesMap(density_consensus_Mk_1000, transition="green+no_green->no_green",type = "l", ylim = c(0, 0.08),colors="grey", cex=2)
par(new=T)
plot.changesMap(density_consensus_Mk_1000, transition="no_green->green+no_green", ylim = c(0, 0.08), colors="lightgreen", cex=2)

plot.changesMap(density_consensus_Mk_1000, transition="no_green->green", ylim = c(0, 0.08), colors="darkgreen", cex=2)
par(new=T)
plot.changesMap(density_consensus_Mk_1000, transition="green->no_green", ylim = c(0, .08), colors="grey", cex=2)

################################################################################
### 5. plot the distribution of time spent in each state
################################################################################
ss<-lapply(green.smap.trees_Mk_1000,summary)
head(ss,2)
str(ss[1])
ss[[1]]$times[2,1]

## time as proportion 
n.times_prop<-t(sapply(ss,function(x) c(x$times[2,1],x$times[2,2],x$times[2,3])))
dimnames(n.times_prop)<-list(1:length(green.smap.trees_Mk_1000),
                            c("green","green+no_green", "no_green"))
n.times_prop = as.data.frame(n.times_prop)
n.times_prop = n.times_prop[,-3]

## reshape the data frame into long format
n.times_prop_long <- n.times_prop %>%
  pivot_longer(cols = everything(), names_to = "group", values_to = "value")

## create the density plot
ggplot(n.times_prop_long, aes(x = value, fill = group)) +
  geom_density(alpha = 0.7,color = NA) +
  labs(title = "Time spent in each state", x = "Proportion of time", y = "Density") +
  scale_fill_manual(values = c("green" = "darkgreen", "green+no_green" = "lightgreen", "no_green" = "grey"))+
  theme_minimal()


