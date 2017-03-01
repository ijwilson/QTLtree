
"Split" <- function(d, cases, SplitPositions, positions, quiet=TRUE)
{
  if (missing(SplitPositions))  SplitPositions <- 1:ncol(d)
  if (missing(positions)) positions <- 1:ncol(d)
  if (length(positions)!=ncol(d))
    stop("length of positions does not match columns")
  if (min(SplitPositions)<1||max(SplitPositions)>ncol(d))
    stop("values in Splitpositions should lie between 1 and ncol(d)")
  if (min(cases)<1||max(cases)>nrow(d))
    stop("values in cases between 1 and number of samples")
  
  n <- nrow(d)
  if (missing(cases))   cases <- (n/2+1):n
  maxedges<- 4*(n-1)
  
  a <- .C("GetSplit"
          ,as.integer(t(d))
          ,as.integer(n)
          ,as.integer(ncol(d))
          ,as.integer(SplitPositions-1)
          ,as.integer(length(SplitPositions))
          ,as.integer(cases-1)
          ,as.integer(length(cases))
          ,edge=as.integer(numeric(maxedges))
          ,ccleaf=as.integer(numeric(2*n))
          ,ccnode=as.integer(numeric(2*n))
          ,leaves=as.integer(numeric(1))
          ,labels=as.integer(numeric(n))
          ,nodepos=as.integer(numeric(2*n))
          ,PACKAGE="snptree")
  
  leaves <- a$leaves
  nedge <- 4*(leaves-1)
  
  edge <- matrix(a$edge[1:(4*(leaves-1))],ncol=2,byrow=FALSE)
  ccnode <- matrix(a$ccnode[1:(2*(leaves-1))],ncol=2,byrow=FALSE) 
  ccleaf <- matrix(a$ccleaf[1:(2*leaves)],ncol=2,byrow=FALSE)
  
  nl <- rowSums(ccleaf)
  #   labs <- tapply(a$labels,rep(1:leaves,nl),paste,collapse=" ")
  labs <- tapply(a$labels,rep(1:leaves,nl),c)
    
  bb <- list(edge=edge,Nnode=leaves-1,edge.length=rep(1,2*(leaves-1)),tip.label=labs)
  class(bb) <- c("phylo")
  
  if (!quiet) cat(leaves," leaves on tree\n")
  
  b <- list(tree=bb,nodepos=a$nodepos[1:(leaves-1)]+1,
            ccnode=ccnode,cctip=ccleaf,n=n,cases=cases,labels=labs)
  
  class(b) <- "split"
  b
}

"splitTest" <- function(d, cases, positions, maxk=4, reps=1000, pickStat="Sevon")
  {
    if (missing(positions))
      positions <- 1:ncol(d)
 
    n <- nrow(d)
    if (missing(cases))   cases <- (n/2+1):n
 
    a <- .C("splitTest"
            ,as.integer(t(d))
            ,as.integer(n)
            ,as.integer(ncol(d))
            ,as.integer(positions-1)
            ,as.integer(length(positions))
            ,as.integer(cases)
            ,as.integer(length(cases))
            ,as.integer(reps)
            ,as.integer(maxk)
            ,teststat=as.double(numeric(maxk))
            ,randteststats=as.double(numeric(maxk*reps))
            ,leaves=as.integer(numeric(1))
            ,as.character(pickStat)
       ,PACKAGE="snptree")
       
    rs <- matrix(a$randteststats,ncol=maxk,byrow=TRUE)
    
    p <-  1.0 - (apply(rbind(a$teststat,rs),2,function(x) rank(x)[1]))/(reps+1.0)    

    list(testStat=a$teststat,randTestStats=rs,leaves=a$leaves,p.value =p)
  }
