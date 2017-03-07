library(snptree)
data(snptreeExample)

simple <- split_simple(haps)
plot(simple)                ## this just plots the numbers for each haplotype

simple$tip.haplotypes
data.frame(haps = simple$tip.haplotypes, n=sapply(simple$labels, length))

## QTL test
trait1 <- rnorm(nrow(haps))                ## get a random trait
s <- split_qtrait(haps, trait1)
class(s)
class(s$tree)
plot(s)
trait2 <- trait1
trait2[s$labels[['45']]] <- trait2[s$labels[['45']]] + 0.15
trait2[s$labels[['41']]] <- trait2[s$labels[['41']]] - 0.15
trait2[s$labels[['36']]] <- trait2[s$labels[['36']]] - 0.15
sb <- split_qtrait(haps, trait2)
plot(sb)

tst <- splitQTLTest(haps, trait2)
tst$p.value
## We would expect that the third p-value would be the the that compares best with
## Other methods, are there are 3 terminal nodes that have a systematic change

## Try with a different statistic.  The default is "A"

tst <- splitQTLTest(haps, trait2, pickStat = "N")
tst$p.value

tst <- splitQTLTest(haps, trait2, pickStat = "P")
tst$p.value

tst <- splitQTLTest(haps, trait2, pickStat = "Z")
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
tr <- split_casecontrol(haps, which(sample=="Case"))
plot(tr)

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

