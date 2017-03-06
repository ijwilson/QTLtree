library(snptree)
data(snptreeExample)

simple <- splitSimple(haps)
plot(simple)
simple$tip.haplotypes
data.frame(haps = simple$tip.haplotypes, n=sapply(simple$labels, length))

## QTL test
trait1 <- rnorm(nrow(haps))                ## get a random trait
s <- split_qtrait(haps, trait1)
class(s)
class(s$tree)
plot(s)
trait2 <- trait1
trait2[s$labels[['45']]] <- trait2[s$labels[['45']]] + 0.1
trait2[s$labels[['41']]] <- trait2[s$labels[['41']]] - 0.2
trait2[s$labels[['36']]] <- trait2[s$labels[['36']]] - 0.2
sb <- split_qtrait(haps, trait2)
plot(sb)

tst <- splitQTLTest(haps, trait2)
tst$p.value


## Try with a different statistic
## now add a small value equal to 1/5 of a standard deviaito nto all those in
## the last haplotype

tst <- splitQTLTest(haps, trait2, pickStat = "N")
tst$p.value

tst <- splitQTLTest(haps, trait2, pickStat = "P")
tst$p.value

tst <- splitQTLTest(haps, trait2, pickStat = "Z")
tst$p.value

tst <- splitQTLTest(haps, trait2, pickStat = "A")
tst$p.value
#################################################
catData <- sample(1:5, nrow(haps), replace=T)
category <- split_category(haps, catData)
class(category)
plot(category)

#################################################

## We can see that there appears to be some relation by looking at a
## Table of haplotypes by Case and control
haplotype <- apply(haps, 1, paste, collapse="")
table(haplotype, sample)
chisq.test(table(haplotype, sample), simulate=TRUE, B=2000)
## Compare this to a split tree
tr <- splitCaseControl(haps, which(sample=="Case"))

firstlabel <- sapply(tr$labels, function(x) x[1])  ## haplotypes at the 
table(haplotype[firstlabel] == rownames(table(haplotype, sample))) ## Order is  the same for 0/1 trees
## tips of the tree are equvalent to the lexical splitting of the haplotypes

tip_haplotypes <- haplotype[firstlabel]

## We should be able to repeat this by looking at the tips of
## out haplotype tree only, whihc we can do by setting maxk=45

s <- splitTestCC(haps, which(sample=="Case"), reps=1000, pickStat="G", maxk=5)
print(s$p.value)

s <- splitTestCC(haps, which(sample=="Case"), reps=1000, pickStat="AbsSevon")
print(s$p.value)

## try with a random sample of cases
random_cases <- sort(sample(nrow(haps), size=nrow(haps)/2, replace=FALSE))
s <- splitTestCC(haps , reps=1000,  cases =random_cases, pickStat="G")
print(s$p.value)

s <- splitTestCC(haps, which(sample=="Case"), reps=10000, pickStat="Gtest")
print(s$p.value)
