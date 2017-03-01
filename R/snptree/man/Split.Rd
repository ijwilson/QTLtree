\name{Split}
\alias{Split}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Split a set of haplotypes}
\description{
  Splits a set of haplotypes into a tree using a lexical sorting.
}
\usage{
Split(d, cases, SplitPositions, positions)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{ A matrix of haplotypes with samples by row }
  \item{cases}{ A matrix giving the labels of the cases.  By default the
    last half of the sample are cases - note that the algorithm is
    symmetric between cases and controls.}
  \item{SplitPositions}{ The order in which to split the haplotypes }
  \item{positions}{ Positions for the sites given in \code{d} }
}

\details{
  ~~ If necessary, more details than the description above ~~
}

\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}

\references{\url{http://www.staff.ncl.ac.uk/i.j.wilson} }

\author{ Ian Wilson \email{Ian.Wilson@ncl.ac.uk}}

\seealso{ \code{\link{splitTest}}}
\examples{
data(snptreeExample)
s <- Split(haps, cases = which(sample=="Case"))
plot(s)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
