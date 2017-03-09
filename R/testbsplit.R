library(snptree)
data(snptreeExample)


## OK remove SNPs with a very low frequency
haps[haps==1] <- -1
haps <- haps+1
snpmean <- colMeans(haps)
haps <- haps[,snpmean>0.02 & snpmean < 0.98]


simple <- bifurcation_split(haps)
plot(simple, edge.width = simple$counts/100)
