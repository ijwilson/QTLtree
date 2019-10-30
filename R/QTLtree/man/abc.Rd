\name{abc}
\alias{abc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Approximate Bayesian Computation  }
\description{
  Get estimates of parameter values for Approximate Bayesian Computation
  (ABC) using the correction of Beaumont et al. (2002).  
}
\usage{
abc(obsStats, simPars, simStats, delta = 0.1, sweep = TRUE, correct = TRUE, plot = FALSE, dis = "rootsumsq")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obsStats}{ A vector of statistics for the observed data }
  \item{simPars}{ A matrix of parameter values used in the simulations }
  \item{simStats}{A matrix of statistics for the simulated parameter
    values -- one per row. }
  \item{delta}{ The proportion of the sample closest to the observed
    statistics to use. }
  \item{sweep}{ Sweep out the means and variances for the statistics? }
  \item{correct}{ Use Beaumont et al. 2002's correction? }
  \item{plot}{ Get plots? }
  \item{dis}{ Which distance measure to use.  Default is root sum of
    squares "rootsumsq", but can be "sumsq", "absdis","max", or "min"}
  }
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  A data.fram with three rows:~Describe the value returned
  corrected:  The corrected parameter values using Beaumont et als
  (2002) correct.  Weighted parameter values using only the Kernal to
  select the proportion without the correction.  The prior values. 
}
\references{ Beaumont, M.A., Zhang, W. & Balding, D. 2002.  Approximate
  Bayesian computation in population genetics. Genetics, 2002, 162, 2025-2035 }
\author{ Ian Wilson \email{I.J.Wilson@ncl.ac.uk} }
\note{ ~~further notes~~
}

\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
data(dosSantos)
getStats <- function(dd) {
    n <- nrow(dd)
    t1 <- table(dd[,1])
    t2 <- table(dd[,2])
    ta <- table(dd[,1],dd[,2])
    c(piA=sum(t1^2)/n^2+sum(t2^2)/n^2,piAB=sum(ta^2)/n^2,nA=length(t1>0)+length(t2>0),nAB=length(ta>0),Dprime=LDhap(dd))
  }

collectStats <- function(param,n=667) {
    require(ARG)
    a <- simARG(n,2,rec=param$rho,growthmodel=paste("exponential(",param$alpha,")"))
    b <- mutateSTR(a,rep(param$theta,2),c(1,2))
    s <- getStats(b$d)
    ARG.remove(a)
    c(s,params=unlist(param))
  }

getpars <- function() {
    list(theta=rgamma(1,3,.3),rho= rgamma(1,1,1),alpha=rgamma(1,1,1))
  }
## not run -- it takes too long
if (FALSE) {
res <- t(replicate(10000,collectStats(getpars())))
res <- res[complete.cases(res),]
tb= abc(getStats(dosSantos),res[,6:8],res[,1:5],plot=T,delta=0.1,dis="sumsq")
tb
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
