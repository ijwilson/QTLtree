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
 arguments{
  \item{d}{ A matrix of haplotypes with samples by row }
  \item{cases}{ A matrix giving the labels of the cases.  By default the
    last half of the sample are cases - note that the algorithm is
    symmetric between cases and controls.}
  \item{SplitPositions}{ ~~Describe \code{SplitPositions} here~~ }
  \item{positions}{ ~~Describe \code{positions} here~~ }
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
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ 
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
data(snptreeExample)
s <- Split(haps, cases = which(sample=="Case"))
plot(s)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
