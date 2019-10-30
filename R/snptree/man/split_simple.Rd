\name{split_simple}
\alias{split_simple}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Split a set of haplotypes}

\description{
  Splits a set of haplotypes into a tree using lexical sorting.
}
\usage{
split_simple(d, SplitPositions, positions, quiet = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{ A matrix of haplotypes with samples by row }
  \item{SplitPositions}{ The order in which to split the haplotypes }
  \item{positions}{ Positions for the sites given in \code{d} }
}

\value{
A list with the following elements.
  \item{tree }{A tree of class \code{ape:phylo} which can be plotted by the \code{ape} library}
  \item{nodepos}{A vector which gives the ordering of the internal nodes of the tree.}
  \item{n}{Total Number of samples on the tree}
  \item{cases}{Vector of indices of cases}
  \item{labels}{a list which gives the indices of all the samples at each leaf.}
  \item{tip_haplotypes}{The haplotypes at the tips of the branches}
}
\references{\url{http://www.staff.ncl.ac.uk/i.j.wilson} 

I. J. Wilson, R. A. J. Howey, D. T. Houniet and M. Santibanez-Koref.  BMC
Proceedings. vol. 5(supp 9):S98 2011.

}
\author{ Ian Wilson \email{Ian.Wilson@ncl.ac.uk}}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{split_casecontrol}}, \code{\link{split_qtrait}}
}
\examples{
data(snptreeExample)
s <- split_simple(haps)
class(s)
plot(s)
data.frame(haps = simple$tip.haplotypes, n=sapply(simple$labels, length))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ tree }


