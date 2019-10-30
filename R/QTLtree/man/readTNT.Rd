\name{readTNT}
\alias{readTNT2D}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Read TNT data}
\description{
  A function to read TNT style data from a file or connection.  The
  format a value for the length of the array, then the rest of the array
  in free format.
  For Array2D data the format gives the number of rows, the number of
  columns
  and then the rest of the matrix in free format.
}
\usage{
readTNT(con, skip = 0, close.con = FALSE)
readTNT2D(con, skip = 0, close.con = FALSE)
writeTNT(data,con)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{con}{ a connection of a filesname }
  \item{skip}{how many lines to skip }
  \item{close.con}{close connection on exit ?}
}
\value{
  returns a vector or matrix (for readTNT and readTNT2D respectively)
  or writes a vector or matrix in TNT style for writeTNT.
}
\author{ Ian J. Wilson}
\examples{
cat("30",rnorm(30),sep="\n",file="/tmp/testfile")
readTNT("/tmp/testfile")
cat("10 3",rnorm(30),file="/tmp/testfileb")
readTNT2D("/tmp/testfileb")
}
\seealso{ \code{\link{simplexml}}}
\keyword{IO }% at least one, from doc/KEYWORDS
\keyword{connections }% __ONLY ONE__ keyword per line


