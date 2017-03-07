\name{split_category}
\alias{split_category}
\alias{split_qtrait}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Produce a bifurcating split of haplotypes that each have a category.}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
split_category(d, categories, SplitPositions, positions, quiet = TRUE)
split_qtrait(d, qtrait, SplitPositions, positions, quiet = TRUE)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{
%%     ~~Describe \code{d} here~~
}
  \item{categories}{
%%     ~~Describe \code{categories} here~~
}
\item{qtrait}{
  \item{qtrait}{A vector giving the continuous trait associated with each individual.}
}
  \item{SplitPositions}{
%%     ~~Describe \code{SplitPositions} here~~
}
  \item{positions}{
%%     ~~Describe \code{positions} here~~
}
  \item{quiet}{
if \code{quiet} is \code{FALSE}, then no extra information is printed from this function.
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


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{split_simple}}, \code{\link{split_casecontrol}}
}
\examples{
data(snptreeExample)

catData <- sample(1:5, nrow(haps), replace=T)
category <- split_category(haps, catData)
class(category)
plot(category)

trait1 <- rnorm(nrow(haps))                ## get a random trait
s <- split_qtrait(haps, trait1)
class(s)
class(s$tree)
plot(s)
trait2 <- trait1
trait2[s$labels[['45']]] <- trait2[s$labels[['45']]] + 0.15
trait2[s$labels[['41']]] <- trait2[s$labels[['41']]] - 0.15
trait2[s$labels[['36']]] <- trait2[s$labels[['36']]] - 0.15
sb <- split_qtrait(haps, trait2)
plot(sb)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ tree }
