## note that this assumes diploids so that we are using
## vectors CaseControl and regions that have length equal to
## ncol(genotypes)/2
subpoprand <- function(genotypes,CaseControl,regions,reps=1000)
  {
    n <- ncol(genotypes)
    nSNP <- nrow(genotypes)
    if (length(CaseControl)==n)
      CaseControl <- CaseControl[seq(1,n,2)]
    if (length(regions)==n)
      regions <- regions[seq(1,n,2)]

   # print(length(CaseControl))
  #  print(n)
  #  print(length(regions))

    if (length(CaseControl)!=round(n/2,0))
      stop("length of CaseControl should equal n/2 or n")

    if (length(regions)!=round(n/2,0))
      stop("length of regions should equal n/2 or n")

    if (!is.factor(CaseControl)) CaseControl <- factor(CaseControl)
    if (!is.factor(regions)) regions <- factor(regions)

    res <- .C("subpoprand"
            ,as.integer(t(genotypes))
            ,as.integer(n/2)
            ,as.integer(nSNP)
            ,as.integer(as.integer(CaseControl)-1)
            ,as.integer(as.integer(regions)-1)
            ,as.integer(reps)
            ,statistic=numeric(nSNP)
            ,p=numeric(nSNP)
            ,PACKAGE="ijwtools"
            ,NAOK=TRUE)

    cbind(statistic=res$statistic,p=(res$p+0.5)/(reps+0.5))
  }
