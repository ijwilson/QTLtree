
## function written for Mauro
## Call the C++ code from R
mutsim <- function(cases=500,
                   penetrance=c(0.0,0.2,1.0),
                   var=40,
                   minfreq=0.45,
                   maxfreq=0.55,
                   minmarkerfreq=0.005,
                   theta=40,
                   alpha=10,
                   rho=10,
                   highratesites=100,
                   seed) {

## returns the *Diploid" sample size required for a simulation
"DipCaseSampleSize" <-
function(n,minfreq,penetrance,pfail=0.00001)
  {
    if (minfreq>0.5)
      stop("Expect the minimum population frequency to be less than 0.5")
    if (length(penetrance) !=3)
      stop("expect three penetrance parameters in .DipCaseSampleSize")
    
    pcase <-  penetrance[3]*minfreq*minfreq+2*penetrance[2]*(1-minfreq)*minfreq
    +penetrance[1]*(1-minfreq)*(1-minfreq)

    q <- qnorm(pfail)
    
   a <- polyroot(c(n*n,-2*n*pcase-q*q*pcase*(1-pcase),pcase*pcase))
    ceiling(Re(a[2]))
    
  }

  library(ijwtools)
  tmpfile=tempfile();
  penetrances <- paste("--pi",c("00","01","11"),"=",penetrance,sep="",collapse=" ")
  if (!missing(seed)) set.seed(seed)
  SampleSizeRequired = DipCaseSampleSize(cases,minfreq,penetrance)
  repeat {
    localSeed <-  floor(runif(1)*1E7)
 
    command <- paste("mutsimrec --var=",var," --high=",highratesites," --ss=",2*SampleSizeRequired," --alpha="
                     ,alpha," --hmtheta=",theta," --seed=",localSeed," --rho=",rho
                     ," | DipCaseControlxml --edge=0 --hr ",penetrances," --ss1=",cases," --ss2="
                     ,cases," --min=",minfreq," --max=",maxfreq," > ",tmpfile,sep="")
    ret <- system(command)
    if (ret==0) break
  }

  a <- readEnrich(tmpfile)$replicate1
  unlink(tmpfile)
  a$locations <- factor(1-a$locations,labels=c("Case","Control"))
  cm <- colMeans(a$haplotypes)
  u <- cm>minmarkerfreq
  a$haplotype <- a$haplotypes[,u]
  a$haplotypes <- NULL
  a$position <- a$position[u]
  a
}



