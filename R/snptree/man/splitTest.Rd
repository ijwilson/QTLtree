\name{splitTest}
\alias{splitTest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Test of Association using a tree of haplotypes}
\description{
   Function to test for assosiation using Case and control status by
 building a tree of haplotypes.  
}
\usage{
splitTest(d, cases, positions, maxk = 4, reps = 1000)
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


}
\author{ Ian Wilson \email{Ian.Wilson@ncl.ac.uk}}
\note{
  The test statistic and p-value will depend on the order in
  which the positions are added to the tree.  In this case it is only
  sesible to add them in a set order either as distance away from a
  fixed point or sequentially along the haplotypes.  Test statistics available are
  
S The Sevon test statistic

Q the squared Sevon statistic

A the absolute Sevon test statistic

G a G-test statistic

P the exact binomial tail probability

N a normalised exact binomial tail probability (G test statistic)

T The Tree test statistic

C The ’Cherries’ test statistic

H The ’Height’ test statistic
  
  
  
}
\examples{
library(snptree)
data(snptreeExample)
s <- splitTest(haps, which(sample=="Case"), reps=1000, pickStat="Sevon")
print(s$p.value)

s <- splitTest(haps, which(sample=="Case"), reps=1000, pickStat="AbsSevon")
print(s$p.value)

s <- splitTest(haps, which(sample=="Case"), reps=10000, pickStat="Gtest")
print(s$p.value)

print(s$p.value)
hapstr <- apply(haps, 1, paste, collapse="")
chisq.test(table(hapstr, sample), simulate=TRUE)

}
\keyword{ htest }
