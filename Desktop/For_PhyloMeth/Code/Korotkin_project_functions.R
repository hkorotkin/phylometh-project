#Counting changes in tip states for discrete trait in ancestral state reconstruction
#Randomizing character traits then sampling over x amount of iterations and then counting the number of changes


trait.changes.asr <- function (tree, traits, nreps=5) {
	res.final <- c()
	for (i in 1:nreps) {
		print(i)	
		tree.tmp <- tree
		random.traits<-traits
		random.traits[,1]<-sample(random.traits[,1],replace=FALSE)
		#this samples the traits in the "traits" column randomly, without replacement
		char1<-random.traits[,1]
		names(char1) <- row.names(traits)
		lik<-make.bisse(tree.tmp,char1)
		p<-starting.point.bisse(tree.tmp)
		fit<-find.mle(lik,p)
		st<-asr.marginal(lik,coef(fit))
		tree.tmp$node.label<-apply(t(st), 1, which.max)
		total.states<-c(random.traits[,1],tree.tmp$node.label)
		anc <- unique(tree.tmp$edge[,1])
		tree.tmp <- reorder(tree.tmp, "pruningwise")
		nb.nodes <- Nnode(tree.tmp)
		change.count <- 0
		for(j in seq(from=1, length.out=nb.nodes)){
			focal <- anc[j]
			desRows <- which(tree.tmp$edge[,1]==focal)
			desNodes <- tree.tmp$edge[desRows,2]
			for(desIndex in sequence(length(desRows))){
				if(!total.states[focal] == total.states[desNodes[desIndex]]){
					change.count = change.count+1
				}
			}			
		}
		res.final <- c(res.final, change.count)
	}
	return(res.final)
}

