
## simulate data from the Ewen's Sampling formula
## and from other simple models for the allele frequency
## spectra of genetic data
rEwens <- function(theta,n)
### Simulate data from the ESF
  {
    x <- 1;
    k <- 2;
    for (i in 2:n) {
      p <- c(x,theta)/(i-1+theta);
      a <- sample(k,1,prob=p);
      if (a==k) {
        k <- k+1;
        x <- c(x,1);
      } else x[a] <- x[a]+1;
    }
    x;  
  }


## this is the Ewens sampling formula when we have
## a[1] alleles represented once, a[2] alleles represented
## twice etc..
dEwens <-  function(a,theta,log=FALSE)
  {
    xx <- 1:length(a)
    ## get the total sample size
    n <- sum(xx*a)
    ld <- lgamma(theta)-lgamma(theta+n)+lgamma(n+1)
    
    l <- a*log(theta) - a*log(xx) -lgamma(a+1)
    
    if (log) return (ld +sum(l))
    else return (exp(ld+sum(l)))
  }

## this is the Ewens sampling formula when we have
## a[1] alleles represented once, a[2] alleles represented
## twice etc..
Stirling <- function(n,k)
  {
    if (n==k) return(1)
    tmp <- matrix(0,nrow=n,ncol=k)
    tmp[1,1] <- 1

    for (i in 2:n) {
      tmp[i,1] <- -(i-1)*tmp[i-1,1]
      if (k>1) {
        for (j in 2:k) {
          tmp[i,j] = tmp[i-1,j-1] - (i-1) * tmp[i-1,j]
        }
      }
    }
    tmp[n,k]
  }
## Density for the number of alleles
dEwensk <-  function(n,k,theta,log=FALSE)
  {
    return( abs(Stirling(n,k))*(theta^k)*gamma(theta)/gamma(n+theta))
  }


