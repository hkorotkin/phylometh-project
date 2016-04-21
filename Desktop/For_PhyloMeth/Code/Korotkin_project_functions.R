#Counting changes in tip states for discrete trait in ancestral state reconstruction
#Randomizing character traits then sampling over x amount of iterations and then counting the number of changes

count <- 0

trait.changes.asr <- function (tree, traits) {
	for (i in 1:1000) {
		random.traits<-traits
		random.traits[,1]<-sample(random.traits[,1],replace=FALSE)
		#this samples the traits in the "traits" column randomly, without replacement
		lik<-make.bisse(tree,random.traits)
		p<-starting.point.bisse(tree)
		fit<-find.mle(lik,p)
		st<-asr.marginal(lik,coef(fit))
		round(coef(fit),5)
	}
	
}





#loop will be the shuffling and bisse functions
replicate()
?sample
#how to count changes?

#table with randomized values to use in the loop. copy original table. sample from randomized column