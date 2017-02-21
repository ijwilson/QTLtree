\name{SplitQTL}
\alias{SplitQTL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Split a set of haplotypes based based on a tree structure
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SplitQTL(d, qtl, SplitPositions, positions, quiet = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{
    The set of haplotypes
    
}
\item{qtl}{
  Values of the quantitative trait loci (the same length as the
  number of haplotypes.  This means values are repeated for
  diploids.
%%     ~~Describe \code{qtl} here~~
}
  \item{SplitPositions}{
%%     ~~Describe \code{SplitPositions} here~~
}
  \item{positions}{
%%     ~~Describe \code{positions} here~~
}
  \item{quiet}{
%%     ~~Describe \code{quiet} here~~
}
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
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (d, qtl, SplitPositions, positions, quiet = TRUE) 
{
    n <- nrow(d)
    if (missing(SplitPositions)) 
        SplitPositions <- 1:ncol(d)
    if (missing(positions)) 
        positions <- 1:ncol(d)
    if (length(positions) != ncol(d)) 
        stop("length of positions does not match columns")
    if (min(SplitPositions) < 1 || max(SplitPositions) > ncol(d)) 
        stop("values in Splitpositions should lie between 1 and ncol(d)")
    if (length(qtl) != n) 
        stop("We need length of qtl to be the same as the number of samples")
    maxedges <- 4 * (n - 1)
    a <- .C("GetQTLSplit", as.integer(t(d)), as.integer(n), as.integer(ncol(d)), 
        as.integer(SplitPositions - 1), as.integer(length(SplitPositions)), 
        as.double(qtl), edge = as.integer(numeric(maxedges)), 
        leaves = as.integer(numeric(1)), leafcount = as.integer(numeric(n)), 
        labels = as.integer(numeric(n)), nodepos = as.integer(numeric(2 * 
            n)), tippos = as.integer(numeric(2 * n)), qtlnode = as.double(numeric(n)), 
        qtlleaf = as.double(numeric(n)), PACKAGE = "ijwtools")
    leaves <- a$leaves
    nedge <- 4 * (leaves - 1)
    edge <- matrix(a$edge[1:(4 * (leaves - 1))], ncol = 2, byrow = FALSE)
    labs <- tapply(a$labels, rep(1:leaves, a$leafcount[1:leaves]), 
        c)
    nodepos <- matrix(a$nodepos[1:(2 * (leaves - 1))], ncol = 2, 
        byrow = T) + 1
    NodePos <- list(left = nodepos[, 1], right = nodepos[, 2])
    tippos <- matrix(a$tippos[1:(2 * leaves)], ncol = 2, byrow = T) + 
        1
    TipPos <- list(left = tippos[, 1], right = tippos[, 2])
    bb <- list(edge = edge, Nnode = leaves - 1, edge.length = rep(1, 
        2 * (leaves - 1)), tip.label = labs)
    class(bb) <- c("phylo")
    if (!quiet) 
        cat(leaves, " leaves on tree\n")
    b <- list(tree = bb, nodepos = NodePos, tippos = TipPos, 
        qtlnode = a$qtlnode[1:(leaves - 1)], qtlleaf = a$qtlleaf[1:leaves], 
        n = n, labels = labs)
    class(b) <- "splitqtl"
    b
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
