
"SplitQTL" <- function(d,qtl,SplitPositions,positions,quiet=TRUE)
{
  if (missing(SplitPositions))  SplitPositions <- 1:ncol(d)
  if (missing(positions)) positions <- 1:ncol(d)
  if (length(positions)!=ncol(d))
    stop("length of positions does not match columns")
  if (min(SplitPositions)<1||max(SplitPositions)>ncol(d))
    stop("values in Splitpositions should lie between 1 and ncol(d)")
  
  n <- nrow(d)
  if (length(qtl) != n)
    stop("We need length of qtl to be the same as the number of samples")
  
  maxedges<- 4*(n-1)
  
  a <- .C("GetQTLSplit"
          ,as.integer(t(d))
          ,as.integer(n)
          ,as.integer(ncol(d))
          ,as.integer(SplitPositions-1)
          ,as.integer(length(SplitPositions))
          ,as.double(qtl)
          ,edge=as.integer(numeric(maxedges))
          ,leaves=as.integer(numeric(1))      # number of leaves
          ,leafcount=as.integer(numeric(n))
          ,labels=as.integer(numeric(n))      # labels at the leaves
          ,nodepos=as.integer(numeric(2*n))   # stuff for ape
          ,tippos= as.integer(numeric(2*n))
          ,qtlnode = as.double(numeric(n))
          ,qtlleaf = as.double(numeric(n))
          ,PACKAGE="ijwtools")
  
  leaves <- a$leaves
  nedge <- 4*(leaves-1)
  
  edge <- matrix(a$edge[1:(4*(leaves-1))],ncol=2,byrow=FALSE)
  
  
  
  labs <- tapply(a$labels,rep(1:leaves,a$leafcount[1:leaves]),c)
  
  nodepos <- matrix(a$nodepos[1:(2*(leaves-1))],ncol=2,byrow=T)+1	
  NodePos <- list(left=nodepos[,1],right=nodepos[,2])
  
  tippos <- matrix(a$tippos[1:(2*leaves)],ncol=2,byrow=T)+1
  TipPos <- list(left=tippos[,1],right=tippos[,2])
  
  bb <- list(edge=edge,Nnode=leaves-1,edge.length=rep(1,2*(leaves-1)),tip.label=labs)
  class(bb) <- c("phylo")
  
  if (!quiet) cat(leaves," leaves on tree\n")
  b <- list(tree=bb,nodepos=NodePos,tippos=TipPos
            ,leafcount=a$leafcount[1:leaves],nodecount=a$leafcount[(leaves+1):(2*leaves)],
            qtlnode=a$qtlnode[1:(leaves-1)],qtlleaf=a$qtlleaf[1:leaves],n=n,labels=labs)
  class(b) <- "splitqtl"
  b
}