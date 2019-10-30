

## linkage disequilibrium calculations
## calculate the linkage disequilibrium for
## a pair of loci (haps is a two column matrix
## with rows representing individuals
LDhap <- function(haps)
  {
    Dprime <- function(pa,pb,Dt)
      {
        Dm <- pmin(outer(pa,pb),outer(1-pa,1-pb) )
        pos <-  pmin(outer(1-pa,pb),outer(pa,1-pb))
        Dm[Dt>0] <- pos[Dt>0]
        abs(Dt/Dm)
      }
    n <- nrow(haps)
    ## first get the allele frequencies
    pA <- table(haps[,1])/n
    pB <- table(haps[,2])/n
    pAB <- outer(pA,pB)
    D <- table(haps[,1],haps[,2])/n -pAB
    Dp <- Dprime(pA,pB,D)

    sum(pAB*Dp)

  }
