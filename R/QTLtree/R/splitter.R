"splittest" <- function(d,cases,reps=1000,positions)
  {
    if (missing(positions))
      positions <- 1:ncol(d)
    # void splittest(int *data, int *samplesize, int *nSNP, int *positions, int *npos, int *cases
	#	 , int *ncases, int *reps, double *pval)
    n <- nrow(d)
    if (missing(cases))   cases <- (n/2+1):n
 

    a <- .C("splittest"
            ,as.integer(t(d))
            ,as.integer(n)
            ,as.integer(ncol(d))
            ,as.integer(positions-1)
            ,as.integer(length(positions))
            ,as.integer(cases)
            ,as.integer(length(cases))
            ,as.integer(reps)
            ,teststat=as.double(numeric(1))
            ,p=as.double(numeric(1))
            ,leaves=as.integer(numeric(1))
       ,PACKAGE="ijwtools")

    list(teststat=a$teststat,p.value=a$p,leaves=a$leaves)
  }
"Split" <- function(d,cases,SplitPositions,positions,quiet=TRUE)
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
            ,tippos= as.integer(numeric(2*n))
            ,PACKAGE="ijwtools")

    leaves <- a$leaves
    nedge <- 4*(leaves-1)

    edge <- matrix(a$edge[1:(4*(leaves-1))],ncol=2,byrow=FALSE)
    ccnode <- matrix(a$ccnode[1:(2*(leaves-1))],ncol=2,byrow=FALSE) 
    ccleaf <- matrix(a$ccleaf[1:(2*leaves)],ncol=2,byrow=FALSE)

    nl <- rowSums(ccleaf)
 #   labs <- tapply(a$labels,rep(1:leaves,nl),paste,collapse=" ")
    labs <- tapply(a$labels,rep(1:leaves,nl),c)

    nodepos <- matrix(a$nodepos[1:(2*(leaves-1))],ncol=2,byrow=T)+1	
    NodePos <- list(left=nodepos[,1],right=nodepos[,2])
 
    tippos <- matrix(a$tippos[1:(2*leaves)],ncol=2,byrow=T)+1
    TipPos <- list(left=tippos[,1],right=tippos[,2])
  
    bb <- list(edge=edge,Nnode=leaves-1,edge.length=rep(1,2*(leaves-1)),tip.label=labs)
    class(bb) <- c("phylo")
   
    if (!quiet) cat(leaves," leaves on tree\n")
    b <- list(tree=bb,nodepos=NodePos,tippos=TipPos,
              ccnode=ccnode,cctip=ccleaf,n=n,cases=cases,labels=labs)
    class(b) <- "split"
    b
  }

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

"SplitCategory" <- function(d,categories,SplitPositions,positions,quiet=TRUE)
  {
    n <- nrow(d)
    if (missing(SplitPositions))  SplitPositions <- 1:ncol(d)
    if (missing(positions)) positions <- 1:ncol(d)
    if (length(positions)!=ncol(d))
      stop("length of positions does not match columns")
    if (min(SplitPositions)<1||max(SplitPositions)>ncol(d))
      stop("values in Splitpositions should lie between 1 and ncol(d)")
    if (length(categories) != n)
      stop("We need length of qtl to be the same as the number of samples")
   
    maxedges<- 4*(n-1)
    ncat <- max(categories)

    a <- .C("GetCategorySplit"
            ,as.integer(t(d))
            ,as.integer(n)
            ,as.integer(ncol(d))
            ,as.integer(SplitPositions-1)
            ,as.integer(length(SplitPositions))
            ,as.integer(categories)
            ,edge=as.integer(numeric(maxedges))
            ,leaves=as.integer(numeric(1))      # number of leaves
            ,leafcount=as.integer(numeric(n))
            ,labels=as.integer(numeric(n))      # labels at the leaves
            ,nodepos=as.integer(numeric(2*n))   # stuff for ape
            ,tippos= as.integer(numeric(2*n))
            ,nodecat = as.integer(numeric(ncat*n))
            ,leafcat = as.integer(numeric(ncat*n))
            ,as.integer(ncat)
            ,PACKAGE="ijwtools")

    #print(a)

    leaves <- a$leaves
    nedge <- 4*(leaves-1)

    edge <- matrix(a$edge[1:(4*(leaves-1))],ncol=2,byrow=FALSE)

    print(length(a$labels))
    print(leaves)
    print(a$leafcount[1:leaves])
    labs <- tapply(a$labels,rep(1:leaves,a$leafcount[1:leaves]),c)

    nodepos <- matrix(a$nodepos[1:(2*(leaves-1))],ncol=2,byrow=T)+1	
    NodePos <- list(left=nodepos[,1],right=nodepos[,2])
 
    tippos <- matrix(a$tippos[1:(2*leaves)],ncol=2,byrow=T)+1
    TipPos <- list(left=tippos[,1],right=tippos[,2])
  
    bb <- list(edge=edge,Nnode=leaves-1,edge.length=rep(1,2*(leaves-1)),tip.label=labs)
    class(bb) <- c("phylo")
   
    if (!quiet) cat(leaves," leaves on tree\n")
    b <- list(tree=bb,nodepos=NodePos,tippos=TipPos,
              nodecat=matrix(a$nodecat,ncol=ncat,byrow=T)[1:(leaves-1),]
              ,leafcat=matrix(a$leafcat,ncol=ncat,byrow=T)[1:leaves,]
              ,n=n
              ,labels=labs
              )
    class(b) <- "splitcategory"
    b
  }


plot.split <-
  function (a, pie = FALSE, length = FALSE) 
{
  require(ape)
  plot(a$tree, show.tip.label = FALSE)
  if (pie == TRUE) {
    nodelabels(pie = a$ccnode)
    tiplabels(pie = a$cctip)
  }
  else {
      nodelabels(a$ccnode[, 2], adj = c(1.3, 1.2), frame = "n", 
                 cex = 0.6)
      nodelabels(a$ccnode[, 1], adj = c(1.3, -0.2), frame = "n", 
                 cex = 0.6)
      tiplabels(a$cctip[, 2], adj = c(1.3, 1.2), frame = "n", 
                cex = 0.6)
      tiplabels(a$cctip[, 1], adj = c(1.3, -0.2), frame = "n", 
                cex = 0.6)
    if (length==FALSE) {
      propcase <- length(a$cases)/a$n
      nn <- rowSums(a$ccnode)
      s <- (a$ccnode[, 1] - nn * propcase)/(sqrt(nn * propcase * 
                                                 (1 - propcase)))
      hl <- cut(s, c(-1e+100, -2, 2, 1e+100))
      cols <- c("green", "white", "red")[as.numeric(hl)]
      nodelabels(round(s, 2), adj = -0.2, cex = 0.6, font = 2, 
                 bg = cols)
      nn <- rowSums(a$cctip)
      s <- (a$cctip[, 1] - nn * propcase)/(sqrt(nn * propcase * 
                                                (1 - propcase)))
      hl <- cut(s, c(-1e+100, -2, 2, 1e+100))
      cols <- c("green", "white", "red")[as.numeric(hl)]
      tiplabels(round(s, 2), adj = -0.2, cex = 0.6, font = 2, 
                bg = cols)

    } else {
      l <- a$nodepos$right-a$nodepos$left
      nodelabels(l , adj = -0.2, cex = 0.6, font = 2)
      lt <- a$tippos$right-a$tippos$left
      tiplabels(lt , adj = -0.2, cex = 0.6, font = 2)
    }
  }
}

"plot.splitqtl" <-
  function (a, length = FALSE) 
{
  require(ape)
  plot(a$tree, show.tip.label = FALSE)


  
  nodelabels(round(a$qtlnode,2),cex=0.5)
  tiplabels(round(a$qtlleaf,2),cex=0.5)



  
}

sortedpos <- function(positions,centre)
  {
    order(abs(positions-centre))
  }

## get the maximum pairwise haplotype lengths for pairs of individuals 
"maxhaplength" <- function(d,SplitPosition,samps,positions)
  {
    
     if (length(samps)!=2)
        stop("Expect two samples")
      if (max(samps)>nrow(d)||min(samps)<1)
        stop("Samples should be between 1 and number of samples")
      if (samps[1]==samps[2]) {
        warning("samples should be different")
        return(0.0)
      }
      if (SplitPosition<1||SplitPosition>ncol(d)) {
        stop("position should be between 1 and sites")
      }
  
     
     sites <- ncol(d)
     ## may want to use position on the chromosome rather than
     ## SNP order!!
     ##   if (SplitPosition<positions[1]) SP <- 1
     ##   else if (SplitPositions>positions[sites]) SP <- sites
     ##  else {
     ##    for (i in 1:sites) if (SplitPosition<positions
     ##  }
     
      
      for (Left in seq(SplitPosition,1,-1)) 
        if (d[samps[1],Left]!=d[samps[2],Left])
          break;

     if (Left==1&&d[samps[1],1]==d[samps[2],1])
        LeftPos <- positions[1]
      else LeftPos <- positions[Left+1]
      
      for (Right in (SplitPosition+1):sites) 
        if (d[samps[1],Right]!=d[samps[2],Right])
          break

      if (Right==sites&&d[samps[1],sites]==d[samps[2],sites])
        RightPos <- positions[sites]
      else RightPos <- positions[Right-1]

     if (RightPos<=LeftPos) return(1)

     RightPos-LeftPos
     # c(LeftPos,RightPos)      
  }

## Split everything from left to right - splitting
## those nodes with no information at random    */
"SimpleSplit" <- function(d)
  {
    require(ape)
    n <- nrow(d)
    maxedges<- 4*(n-1)
    
    a <- .C("simplesplit"
            ,as.integer(t(d))
            ,as.integer(n)
            ,as.integer(ncol(d))
            ,edge=as.integer(numeric(maxedges))
            ,labels=as.integer(numeric(n))
            ,PACKAGE="ijwtools")

    bb <- list(edge= matrix(a$edge,ncol=2,byrow=FALSE)
               ,Nnode=n-1,edge.length=rep(1,2*(n-1)),tip.label=a$labels)
    class(bb) <- c("phylo")
    bb
  }
