#Written by: Dylan DeBaun

library(ape)
library(phytools)
library(ggplot2)

#Create a phylogentic tree representing the Reticulations through time so that we can use already established ltt analysis tools.

rtt <- read.table() #create a data file where the first column is timing of reticulation in order 

namertt <- stri_rand_strings(dim(rtt)[1], 3, pattern = "[A-Z]") #assign arbitrary names for creating the newick file

newick = paste0("(",namertt[1],":",abs(rtt[1,1]),", ~:",abs(rtt[1,1]),")")
for(i in 2:dim(rtt)[1]){
  next1 = paste0("(",namertt[i],":",abs(rtt[i,1]),", ~:",abs(rtt[i,1]),")")
  newick = gsub("~", next1,newick)
}
newick = gsub("~", "Z",newick)
newick = paste0(newick,";") #creates the newick text

newtree = read.tree(text = newick) #reads that in as a tree
t <- force.ultrametric(newtree) #extend the tips to extant 
ltt(t,plot=T, log.lineages = T) #create the reticulation through time plot

#simulate RTTs under pure birth model (http://phytools.org/mexico2018/ex/10/Diversification-models.html)
trees<-pbtree(b=b,n=Ntip(t),t=h,nsim=100,method="direct",
              quiet=TRUE)
obj<-ltt(trees,plot=FALSE)
plot(obj,col="grey",main="LTT compared to simulated LTTs")
lines(c(0,h),log(c(2,Ntip(t))),lty="dashed",lwd=2,col="red")
ltt(t,add=TRUE,lwd=2)

#plot the 95th confidence envelope for those birth rtts and overlay the true rtt
ltt95(trees,log=TRUE)
title(main="RTT compared to simulated RTTs")
ltt(t,add=TRUE,log.lineages=FALSE,col="red",lwd=2)
ninetyfive <-ltt95(trees,log=TRUE)

#Find rates for the rtt under two different models (birth and birth death): https://lukejharmon.github.io/ilhabela/instruction/2015/06/02/diversification/
library(diversitree)
# first fit a Yule model
pbModel <- make.yule(t)
pbMLFit <- find.mle(pbModel, 0.1)

# next fit a Birth-death model
bdModel <- make.bd(t)
bdMLFit <- find.mle(bdModel, c(0.1, 0.05), method = "optim", lower = 0)

# compare models
anova(bdMLFit, pure.birth = pbMLFit)

#run mcmc on the better fit model to come up with a distribution of values for the parameter estimation
bSamples <- mcmc(pbModel, pbMLFit$par, nsteps = 1e+05, lower = c(0), upper = c(Inf), w = c(0.1), fail.value = -Inf, print.every = 10000)
postSamples <- bSamples[c("lambda")]
profiles.plot(postSamples, col.line = c("red"), las = 1, legend = "topright")
#get the mean parameter estimates
mean(postSamples$lambda)


