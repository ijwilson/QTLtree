\name{plot.split}
\alias{plot.split}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot a bifucating tree produced using \code{Split} using the \code{ape} library.
}
\description{
Plots the tree.  This plot function generally gets the cases and control for each point on the 
tree. 
}
\usage{
plot.split(a, pie = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{
 The object of class \code{snptree}}
  \item{pie}{
  Put a pie chart of relative frequencies of cases and controls at the nodes?
  }
}
\details{
 If you just want the raw tree then using the \code{plot(object.tree, show.tip.label=FALSE)} works well.
 }
\value{
No value returned.  This function called for side effect.
}
\references{\url{http://www.staff.ncl.ac.uk/i.j.wilson} }
\author{ Ian Wilson \email{Ian.Wilson@ncl.ac.uk}}


\seealso{
\code{\link{split_casecontrol}}}, \code{\link{split_simple}}}, \\code{\link{split_category}}}, \code{\link{split_qtrait}}}
\examples{
data(snptreeExample)
s <- split_casecontrol(haps, which(sample=="Case"), position)
plot(s)
s2 <- split_casecontrol(haps, which(sample=="Case"), rev(position))
plot(s2)
plot(s, pie=TRUE)
}
\keyword{ graphics }
