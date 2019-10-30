"rowdiffs" <- function(m,maxlen)
{
  if (missing(maxlen)) 
    maxlen <- ncol(m)
  
  .C("rowdiffs"
     ,as.integer(t(m))
     ,as.integer(nrow(m))
     ,as.integer(ncol(m))
    ,res=as.integer(numeric(maxlen))
     ,as.integer(maxlen)
     ,PACKAGE="ijwtools")[[4]]

}

"coldiffmat" <- function(m)
{
  nc <- ncol(m)
  len <- (nc)*(nc-1)/2
  
  r <- .C("coldiffmat"
     ,as.integer(t(m))
     ,as.integer(nrow(m))
     ,as.integer(ncol(m))
     ,res=as.integer(numeric(len))
     ,PACKAGE="ijwtools")[[4]]

  a <- matrix(ncol=nc,nrow=nc) 
  a[lower.tri(a)] <- r
  t(a)
}

"rowdiffmat" <- function(m)
{
  nr <- nrow(m)
  len <- (nr)*(nr-1)/2
  
  r <- .C("coldiffmat"
     ,as.integer(m)
     ,as.integer(ncol(m))
     ,as.integer(nrow(m))
     ,res=as.integer(numeric(len))
     ,PACKAGE="ijwtools")[[4]]

  a <- matrix(ncol=nr,nrow=nr) 
  a[lower.tri(a)] <- r
  t(a)
}
