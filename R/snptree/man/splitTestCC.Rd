\name{splitTestCC}
\alias{splitTestCc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Test of Association using a tree of haplotypes}
\description{
   Function to test for assosiation using Case and control status by
 building a tree of haplotypes.  
}
\usage{
splitTestCC(d, cases, positions, maxk = 4, reps = 1000)
}
\arguments{
  \item{d}{ A matrix of haplotypes with samples by row }
  \item{cases}{ A matrix giving the labels of the cases.  By default the
    last half of the sample are cases - note that the algorithm is
    symmetric between cases and controls.}
  \item{reps}{The number of permutations to test the significance with}
  \item{positions}{ The positions (and their order) to build the
    haplotype tree}
  \item{maxk}{The maximum value for k in the test statistic}
    \item{pickStat}{The Statistic used to test significance in the randomisation test}

}
\details{
  The test is currently written with test statistic which is the
  number of nodes on the tree (both internal and leaves which contain
  more than one haplotype) that include both cases and controls..  Extreme
  values of this test statistic are low.  This test is done for haplotypes, so is
 suitable for mtDNA and Y chromosome based tests.  For autosomal
 haplotypes one would want to use a randomisation procedure that
 randomised the two haplotypes at an individual together.
}
\value{
  A vector with elements:
  \item{teststat }{The test statistics for values of k from 1 to maxk}
  \item{randTestStat}{The randomised Test Statistics}
  \item{p.value}{The p value (from simulation) from 1 to maxk}
  \item{leaves}{The number of leaves - the number of multi-locus haplotypes}
}
\references{\url{http://www.staff.ncl.ac.uk/i.j.wilson} 

I. J. Wilson, R. A. J. Howey, D. T. Houniet and M. Santibanez-Koref.  BMC
Proceedings. vol. 5(supp 9):S98 2011.

Recent mitochondrial DNA mutations increase the risk of developing common late-onset human diseases
G Hudson, A Gomez-Duran, IJ Wilson, PF Chinnery
PLoS Genet 10 (5), e1004369

Sevon P, Toivonen H, Ollikainen V: TreeDT: tree pattern mining for gene mapping. IEEE/ACM Trans Comput Biol Bioinform. 2006, 3: 174-185.

}
\author{ Ian Wilson \email{Ian.Wilson@ncl.ac.uk}}
\note{
  The test statistic and p-value will depend on the order in
  which the positions are added to the tree.  In this case it is only
  sesible to add them in a set order either as distance away from a
  fixed point or sequentially along the haplotypes.  Test statistics available are
  
\describe{
\item{S}{The Sevon test statistic, described in Sevon et al}
\item{Q}{The sQuared Sevon statistic.  The Sevon statistic Squared}
\item{A}{The Absolute Sevon test statistic.  The absolute value of the Sevon statistic.}
\item{G}{a G-test statistic}
\item{P}{The exact binomial tail probability}
\item{N}{A Normalised exact binomial tail probability (G test statistic)}
\item{T}{The Tree test statistic}
\item{C}{The ’Cherries’ test statistic}
\item{H}{The ’Height’ test statistic}
}

}
\examples{
library(snptree)
data(snptreeExample)
## We can see that there appears to be some relation by looking at a
## Table of haplotypes by Case and control
haplotype <- apply(haps, 1, paste, collapse="")
table(haplotype, sample)
chisq.test(table(haplotype, sample), simulate=TRUE, B=2000)
## We should be able to repeat this by looking at the tips of
## out haplotype tree only, whihc we can do by setting maxk=45
s <- splitTestCC(haps, which(sample=="Case"), reps=1000, pickStat="G", maxk=45)
print(s$p.value)

s <- splitTestCC(haps, which(sample=="Case"), reps=1000, pickStat="AbsSevon")
print(s$p.value)

## try with a random sample of cases
random_cases <- sort(sample(nrow(haps), size=nrow(haps)/2, replace=FALSE))
s <- splitTestCC(haps , reps=1000, pickStat="G")
print(s$p.value)

s <- splitTestCC(haps, which(sample=="Case"), reps=10000, pickStat="Gtest")
print(s$p.value)

}
\keyword{ htest }
