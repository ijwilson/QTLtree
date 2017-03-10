library(snptree)
data(snptreeExample)


## OK remove SNPs with a very low frequency
haps[haps==1] <- -1
haps <- haps+1
snpmean <- colMeans(haps)
haps <- haps[,snpmean>0.02 & snpmean < 0.98]


simple <- bifurcation_split(haps)
simple$root.edge=15
plot(simple, root.width=100, edge.width = c(simple$counts/100,100), root.edge=T, type="cladogram")

## OK now need to be able to directly read 

nleaves <- Ntip(simple)

TIPS <- simple$edge[simple$edge[,2] <= nleaves, 2] ## gives the ordering of TIPS
                                                   ## For my wel lbehaved trees this is 1:nleaves but it needed be


