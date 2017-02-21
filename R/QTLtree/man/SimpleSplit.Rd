\name{SimpleSplit}
\alias{SimpleSplit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Haplotype splitting of a set of data }
\description{
  Split a matrix of haplotypes into a haplotype tree and return the
  tree.  Mainly used as a test of the method.
}
\usage{
SimpleSplit(d)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{ A matrix of haplotypes.}
}
\details{
  The splitting is done left to right over all the SNPs.  If
  not all terminal nodes are of size one then the remainder are
  split in an arbitrary manner.
}
\value{
  An tree of class \code{phylo}.
  
}
%\references{ ~put references to the literature/web site here ~ }
\author{ Ian Wilson \email{I.J.Wilson@ncl.ac.uk}
  }

\seealso{ \code{\link{Split}} }
\examples{
library(ARG)
a <- mutate(simARG(100,1000,0.01),100)
plot(SimpleSplit(a$haplotype))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
