\name{plot.split}
\alias{plot.split}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot a bifucating tree produced using \code{Split} using the \code{ape} library.
}
\description{
Plots the tree
}
\usage{
plot.split(a, pie = FALSE, length = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{
 The object of class \code{snptree}}
  \item{pie}{
  Put a pie chart of relative frequencies of cases and controls at the nodes?
  }
  \item{length}{
%%     ~~Describe \code{length} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{\url{http://www.staff.ncl.ac.uk/i.j.wilson} }
\author{ Ian Wilson \email{Ian.Wilson@ncl.ac.uk}}
\note{
%%  ~~further notes~~
}


\seealso{
\code{\link{Split}}}
\examples{
data(snptreeExample)
s <- Split(haps,which(sample=="Case"), position)
plot(s)
}
\keyword{ graphics }
