\name{rowdiffmat}
\alias{coldiffmat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Get Manhattan distances between rows and columns of matrices }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rowdiffmat(m)
coldiffmat(m)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{ A matrix }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  An upper triangular matrix of differences between rows or columns
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(m,maxlen)
{
  nr <- nrow(m)
  len <- (nr)*(nr-1)/2
  
  r <- .C("coldiffmat"
     ,as.integer(m)
     ,as.integer(ncol(m))
     ,as.integer(nrow(m))
     ,res=as.integer(numeric(len))
     ,PACKAGE="ijwtools")[[4]]

  a <- matrix(ncol=nr,nrow=nr) 
  a[lower.tri(a)] <- r
  t(a)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
