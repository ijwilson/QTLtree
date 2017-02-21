\name{f2}
\alias{f2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
f2(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ ~~Describe \code{x} here~~ }
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

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(x) {
    Fv=x[-npops]
    thet=x[npops]
    if (min(F)<0.0)  return(1E10)
    else if (max(F)>=1) return(1E10)
    if (thet<=0.0)  return(1E10)
    set.seed(seed);
    return(
           -.C("likelihoodDriven"
               ,as.integer(xx)
               ,as.integer(location)
               ,as.integer(length(location))
               ,as.double(F)
               ,as.double(driveF)
               ,as.integer(nperms)
               ,as.double(theta)
               ,res=as.double(numeric(1))
               )$res
           )
    }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
