\name{splitQTLTest}
\alias{splitQTLTest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Test of Association between haplotype and a quantitative trait}
\description{
   Function to test for association using qtl by building a tree of haplotypes.  
}
\usage{
splitQTLTest(haplotypes, qtrait, positions, maxk = 4, reps = 1000)
}
\arguments{
  \item{haplotypes}{ A matrix of haplotypes with samples by row }
  \item{qtrait}{A vector giving the continuous trait associated with each individual.}
  \item{reps}{The number of permutations to test the significance with}
  \item{positions}{ The positions (and their order) to build the
    haplotype tree}
  \item{maxk}{The maximum value for k in the test statistic}
  \item{pickStat}{The test statistic to use}
}
\details{
  The test is currently written with test statistic which is the
  number of nodes on the tree (both internal and leaves which contain
  more than one haplotype) that include both cases and controls.  Extreme
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
  \item{pickStat}{The Statistic used to test significance in the randomisation test}
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
\item{Z}{The Sevon test statistic, described in Sevon et al}
\item{A}{The Absolute Sevon test statistic.  The absolute value of the Sevon statistic.}
\item{P}{a G-test statistic}
\item{N}{A Normalised exact binomial tail probability (G test statistic)}
}

}
\examples{
library(snptree)
data(snptreeExample)
## We can see that there appears to be some relation by looking at a
## Table of haplotypes by Case and control
library(snptree)
data(snptreeExample)
s <- split_simple(haps)
## We can see that there appears to be some relation by looking at a
## Table of haplotypes by Case and control


trait1 <- rnorm(nrow(haps))
tst <- splitQTLTest(haps, trait1)
tst$p.value

## now add a small value equal to 1/5 of a standard deviation nto all those in
## the last haplotype
trait1[s$labels[['45']]] <- trait1[s$labels[['45']]] + 0.2
tst <- splitQTLTest(haps, trait1)
tst$p.value

}
\keyword{ htest }
