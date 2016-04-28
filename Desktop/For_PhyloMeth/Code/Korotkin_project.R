#Get BiSSE

#setwd("/Users/Lepiota/Desktop/For_PhyloMeth")
setwd("/Users/Hailee/Desktop/For_PhyloMeth")

library(ape)
library(diversitree)
library(geiger)
library(phylotools)
library(deSolve)
library(subplex)
library(Rcpp)
library(phangorn)
library(corHMM)
library(phytools)
library(MuMIn)
library(rotl)
library(laser)

tree<-read.tree("tree_rooted.tre")
ultra<-chronopl(tree,lambda=0)
write.nexus(ultra, file="ultrametric1.tre")
plot(ultra,cex=0.7)

traits<-read.csv("traits.csv",row.names=1)
char1<-traits[,1]
names(char1)<-row.names(traits)
name.check(ultra,char1)


charLabel<-character(length(ultra$tip.label))
names(charLabel)<-names(char1)
charLabel[char1==0]<-"red"
charLabel[char1==1]<-"green"
plot(ultra,cex=0.4, label.offset=0.008, no.margin=TRUE, show.tip.label = TRUE)
tiplabels(pch=21,cex=0.8, bg=charLabel[ultra$tip.label],lwd=0.25)

#BiSSE
#full model
lik<-make.bisse(ultra,char1)
argnames(lik) # to see the arguments of lik
p <- starting.point.bisse(ultra)
fit<-find.mle(lik,p)
st<-asr.marginal(lik,coef(fit))
round(coef(fit),5)
plot(ultra,st,type="phylogram",cex=0.6)
nodelabels(thermo=t(st),piecol=c("blue","red"),cex=.4)


trait.changes.asr(ultra,traits,3)
