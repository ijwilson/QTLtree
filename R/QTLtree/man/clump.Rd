\name{clump}
\alias{clump}
\alias{gclump}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to do a clumped Chi-squared test of association}
\description{
  \code{clump} is a port to R of the c code of David Curtis for finding the
  maximum chi-squared test statistic for a table and then testing
  significance by getting the simulated p-value.  \code{gclump} is my
  version of the code which uses the G test statistic (which is
  additive). 
}
\usage{
clump(xxx, nsims = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xxx}{ A two row matrix of counts }
  \item{nsims}{ The number of Monte-Carlo rearrangements for simulated p-values}
}
\details{
  From document distributed with original source:
CLUMP is a program designed to assess the significance of the
departure of observed values in a contingency table from the
expected values conditional on the marginal totals. The present
implementation works on 2 x N tables and was designed for use in
genetic case-control association studies, but the program should
be useful for any 2 x N contingency table, especially where N is
large and the table is sparse. The significance is assessed
using a Monte Carlo approach, by performing repeated simulations
to generate tables having the same marginal totals as the one
under consideration, and counting the number of times that a
chi-squared value associated with the real table is achieved by
the randomly simulated data. This means that the significance
levels assigned should be unbiased (with accuracy dependent on
the number of simulations performed) and that no special account
needs to be taken of continuity corrections or small expected
values. The method is described in full in: Sham PC & Curtis D.
1995. Monte Carlo tests for associations between disease and
alleles at highly polymorphic loci. Ann Hum Genet. 59: 97-105.
Please cite this reference when using the CLUMP program.

}
\value{
  A data frame with test statistics, simulated p-values and p-values
  where appropriate.
}
\references{
  Sham PC & Curtis D. 1995. Monte Carlo tests for associations between
  disease and alleles at highly polymorphic loci. Ann Hum Genet. 59:
  97-105.
  
  \url{http://www.mds.qmul.ac.uk/statgen/dcurtis/software.html}
}
\author{ Ian Wilson \email{I.J.Wilson@ncl.ac.uk} }
%\note{ ~~further notes~~ 

% ~Make other sections like Warning with \section{Warning }{....} ~
%}

\examples{
aa = matrix(c(1,2,1,1,2,2,3,2,1,1,4,0,0,3,0,0,1,0,4,0),nrow=2,byrow=TRUE)
clump(aa)
gclump(aa)
data(centenarians)
attach(centenarians)
clump(xtabs(count~status+haplogroup,subset=(area=="North"&sex=="Male")))
fisher.test(xtabs(count~status+haplogroup,subset=(area=="North"&sex=="Male")))

}
\keyword{ htest }

