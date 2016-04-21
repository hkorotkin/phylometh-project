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
#gives me the exact same results every time
	#is it because the starting point is always the same?


replicate(10,st)
?replicate
?starting.point.bisse
?find.mle




























#fit<- find.mle(lik, p)
#fit[1:2]
#round(coef(fit),5)



#lik<-make.bisse(ultra,char1)
#fit<-find.mle(lik,p,method="subplex")
#st<-asr.marginal(lik,coef(fit))
#nodelabels(thermo=t(st),piecol=c("blue","red"),cex=.5)
#new.st<-as.data.frame(st)

?plot.new
?make.simmap
?trait.plot
?sim.history
?nodelabels



























################################################################################
#Full model
#accounting for different sampling depending on the character
# old sampling - sampling.f <-as.numeric(c(226/1650,69/300))
sampling.f <-as.numeric(c(10/260,38/429))
lik.s<-make.bisse(ultra,char1, sampling.f=sampling.f)
p<-starting.point.bisse(ultra)
fit.s<-find.mle(lik.s,p)
fit.s[1:2]
round(coef(fit.s),5)

#to do all sampling together with 10000 generations
# Full model
sampling.f <-as.numeric(c(10/260,38/429))
p<-starting.point.bisse(ultra)
prior<-make.prior.exponential(1/(2*(p[1]-p[3])))

lik.s<-make.bisse(ultra,char1, sampling.f=sampling.f)
fit.s<-find.mle(lik.s,p)
tmp <- mcmc(lik.s, fit.s$par, nsteps=100, prior=prior, lower=0, w=rep(1, 6), print.every=0)
w<-diff(sapply(tmp[2:7],range))
samples<-mcmc(lik.s, fit.s$par, nsteps=10000, w=w, lower=0, prior=prior, print.every=1000)
