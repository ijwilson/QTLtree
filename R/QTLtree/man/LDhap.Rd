\name{LDhap}
\alias{LDhap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculate the multiallelic D'}
\description{
  Returns the multiallelic D' as defined by Lewontin (1964) for a pair of loci.
}
\usage{
LDhap(haps)
}

\arguments{
  \item{haps}{ A matrix with two columns }
}

\value{
  Returns D'
}
\references{ Lewontin R.C. (1964) Genetics.}
  \author{Ian Wilson }
\examples{
a <-
structure(c(5, 2, 3, 3, 3, 3, 2, 1, 3, 1, 4, 3, 5, 3, 1, 2, 1, 
5, 3, 3, 3, 2, 2, 4, 2, 3, 1, 1, 4, 1, 2, 3, 4, 2, 4, 2, 1, 2, 
2, 1, 2, 4, 4, 1, 1, 1, 4, 3, 1, 4, 3, 4, 2, 1, 4, 4, 3, 1, 3, 
4, 4, 1, 4, 3, 3, 4, 3, 4, 3, 4, 2, 3, 4, 3, 1, 3, 4, 4, 6, 3, 
3, 1, 1, 3, 2, 2, 1, 2, 3, 2, 1, 1, 2, 2, 4, 1, 2, 3, 2, 1),
.Dim = as.integer(c(20, 
5)))
LDhap(a[,c(1:2)])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ univar }
