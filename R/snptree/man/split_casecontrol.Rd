\name{split_casecontrol}
\alias{split_casecontrol}
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


\value{
A list with the following elements.
  \item{tree }{A tree of class \code{ape:phylo} which can be plotted by the \code{ape} library}
  \item{nodepos}{A vector which gives the ordering of the internal nodes of the tree.}
  \item{ccnode}{Matrix of the number of cases and controls at each internal node of the tree.}
  \item{cctip}{Matrix of the numbers of cases and controls at each leaf}
  \item{n}{Total Number of samples on the tree}
  \item{cases}{Vector of indices of cases}
  \item{labels}{a list whihc gives the indicies of all the samples at each leaf.}
}

\references{\url{http://www.staff.ncl.ac.uk/i.j.wilson} 

I. J. Wilson, R. A. J. Howey, D. T. Houniet and M. Santibanez-Koref.  BMC
Proceedings. vol. 5(supp 9):S98 2011.

}

\author{ Ian Wilson \email{Ian.Wilson@ncl.ac.uk}}

\seealso{ \code{\link{splitTest}}}
\examples{
data(snptreeExample)
s <- split_casecontrol(haps, cases = which(sample=="Case"))
class(s)
plot(s)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ tree }
