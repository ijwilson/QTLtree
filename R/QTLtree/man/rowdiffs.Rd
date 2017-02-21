\name{rowdiffs}
\alias{rowdiffs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fund the absolute difference between rows for a 0/1 matrix}
\description{
  Used for various things for SNP data, this function calculates
  frequency spectrum of row
  difference for a matrix.  Can also be used for the column difference
  by using the transpose.
}
\usage{
rowdiffs(m, maxlen)
}
\arguments{
  \item{m}{ a 0/1 matrix.}
  \item{maxlen}{ The maximum value to keep track of.}
}
\value{
  A vector which gives the frequency spectrum of Manhattan differences
  between  rows of a matrix.
}
\author{ Ian Wilson }
\examples{
data(b)
rd <- rowdiffs(b[[1]]$haplotypes)
barplot(rd[1:49])

rd1 <- rowdiffs(b[[1]]$haplotypes[1:40,])
rd2 <- rowdiffs(b[[1]]$haplotypes[41:80,])
barplot(rbind(rd1[1:49],rd2[1:49]),beside=TRUE,col=c("blue","cyan"))

}
\seealso{ \link{rowdiffmat}, \link{coldiffmat}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
