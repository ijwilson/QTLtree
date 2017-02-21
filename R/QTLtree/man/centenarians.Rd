\name{centinarians}
\alias{centinarians}
\docType{data}
\title{ Longievity and mtDNA haplogroup in humans.}
\description{
  Data on Longievity and mtDNA haplogroup in Italians.  The data are divided by sex and
by whether they are from Northern of Southern Italy.
}
\usage{data(centenarians)}
\format{
  A data frame with 80 observations on the following 5 variables.
  \describe{
    \item{\code{haplogroup}}{a factor with levels \code{H} \code{I} \code{J} \code{K} \code{T} \code{U} \code{V} \code{W} \code{X} \code{Other}}
    \item{\code{count}}{The number of individuals}
    \item{\code{status}}{The disease status, a factor with levels \code{Control} \code{Case}.  Cases are centenarians, 
Controls younger individuals.}
    \item{\code{sex}}{A factor with levels \code{Female} \code{Male}}
    \item{\code{area}}{The location in Italy.  A factor with levels \code{North} \code{South}}
  }
}
\source{
 De Benedictis et al.  1999.  Mitochondrial DNA inherited variants are associated with
successful ageing and longevity in humans.  FASEB 13: 1532-1536.
}
\examples{
data(centenarians)

}
\keyword{datasets}
