\name{simplexml}
\alias{simplexml}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Read certains types of simple xml file }
\description{
  Reads simple xml files and returns a list with the tag names producing
  the tag names.   It reads until the end tag of \code{opening}.
}
\usage{
simplexml(conn, texttags = NULL, numbertags = NULL, vectags = NULL, mattags = NULL, opening)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{conn}{ A connection }
  \item{texttags}{ tags that should be read (if seen) as strings }
  \item{numbertags}{ tags that should be read as single numbers }
  \item{vectags}{ Tags read as TNT style arrays }
  \item{mattags}{ Tags read as TNT style Array2Ds }
  \item{opening}{ The opening tag - read until this closes }
}
\details{
  If necessary, more details than the description above 
}
\value{
  A list with entries matching the tags that have been read from the
  possible ones given to the function.
}
\references{  For the TNT C++ library \link{ http://math.nist.gov/tnt/}
  }
\author{ Ian Wilson}

\seealso{  \code{\link{readTNT}},\code{\link{readTNT2D}} }
\examples{
tmp=tempfile()
cat("<mat>2 2 1 2 3 4\n</mat>","<text>A test string</test>"
,"</start>",file=tmp,sep="")
simplexml(file(tmp,"r"),texttags=c("text"),mattags=c("mat"),opening="start")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
