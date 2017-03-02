library(snptree)
data(snptreeExample)
qtl <- rnorm(nrow(haps))
s <- SplitQTL(haps, qtl)
class(s)
plot(s)

catData <- sample(1:5, nrow(haps), replace=T)

category <- SplitCategory(haps, catData)
