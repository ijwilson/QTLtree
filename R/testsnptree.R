library(snptree)
data(snptreeExample)
qtl <- rnorm(nrow(haps))
s <- SplitQTL(haps, qtl)

qtl <- rnorm(nrow(haps))
qtl[s$labels[['45']]+1] <- qtl[s$labels[['45']]+1] + 0.1
tst <- splitQTLTest(haps,qtl)
tst$p.value

## now add a small value equal to 1/5 of a standard deviaito nto all those in
## the last haplotype
qtl[s$labels[['45']]+1] <- qtl[s$labels[['45']]+1] + 0.1
tst <- splitQTLTest(haps,qtl)
tst$p.value

class(s)
plot(s)

## Try with a different statistic
## now add a small value equal to 1/5 of a standard deviaito nto all those in
## the last haplotype

tst <- splitQTLTest(haps,qtl, pickStat = "A")
tst$p.value


tst <- splitQTLTest(haps,qtl, pickStat = "G")
tst$p.value


## category needs to be fixed to remove extra information being printed.
catData <- sample(1:5, nrow(haps), replace=T)
category <- SplitCategory(haps, catData)
#################################################
library(snptree)
data(snptreeExample)
## We can see that there appears to be some relation by looking at a
## Table of haplotypes by Case and control
haplotype <- apply(haps, 1, paste, collapse="")
table(haplotype, sample)
chisq.test(table(haplotype, sample), simulate=TRUE, B=2000)
## Compare this to a split tree
tr <- Split(haps, which(sample=="Case"))

firstlabel <- sapply(tr$labels, function(x) x[1]+1 )  ## haplotypes at the 
table(haplotype[firstlabel] == rownames(table(haplotype, sample))) 
## tips of the tree are equvalent to the lexical splitting of the haplotypes

tip_haplotypes <- haplotype_string[firstlabel,]

## We should be able to repeat this by looking at the tips of
## out haplotype tree only, whihc we can do by setting maxk=45

s <- splitTest(haps, which(sample=="Case"), reps=1000, pickStat="G", maxk=45)
print(s$p.value)

s <- splitTest(haps, which(sample=="Case"), reps=1000, pickStat="AbsSevon")
print(s$p.value)

## try with a random sample of cases
random_cases <- sort(sample(nrow(haps), size=nrow(haps)/2, replace=FALSE))
s <- splitTest(haps , reps=1000, pickStat="G")
print(s$p.value)

s <- splitTest(haps, which(sample=="Case"), reps=10000, pickStat="Gtest")
print(s$p.value)



haplotypestring <- function(tr, )
