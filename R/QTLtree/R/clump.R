"chiStat" <- function(x) {
  ## my own calculation of the chi-squared statistic
  sr <- rowSums(x)
  sc <- colSums(x)
  E <- outer(sr, sc, "*")/sum(sr)
  ##  chisq.test(table(sample$d[,locus],sample$location))$statistic
  sum(sort((x - E)^2/E, decreasing = TRUE))
}

"clump" <- function(xxx,nsims=1000) {
  if (nrow(xxx)!=2)
    stop("should pass matrix with 2 rows to clump")
  
  res <- .C("clumpR"
           ,as.double(t(xxx))
            ,as.integer(ncol(xxx))
            ,as.integer(nsims)
            ,teststat=as.double(numeric(4))
            ,over=as.integer(integer(4))
            ,redcol=as.integer(integer(1))
            ,PACKAGE="ijwtools"
          )
  p <- c(1.-pchisq(res$teststat[1],ncol(xxx)-1),1.-pchisq(res$teststat[2],res$redcol-1),NA,NA)
  a <- cbind(testStatistic=res$teststat,simulatedp=(res$over+1)/(nsims+1),pvalue=p)
  rownames(a) <- c("all","Cells with small totals clumped","Clumped to maximise chi-squared","Maximum for each column")
  a
}


"gclump" <- function(xxx,nsims=1000) {
  if (nrow(xxx)!=2)
    stop("should pass matrix with 2 rows to clump")
  
  res <- .C("gclump"
           ,as.double(t(xxx))
            ,as.integer(ncol(xxx))
            ,as.integer(nsims)
            ,teststat=as.double(numeric(4))
            ,over=as.integer(integer(4))
            ,PACKAGE="ijwtools"
            )
  p <- c(1.-pchisq(res$teststat[1],ncol(xxx)-1),NA,1.-pchisq(res$teststat[3],ncol(xxx)-1),NA)
  a <- cbind(testStatistic=res$teststat,simulatedp=(res$over+1)/(nsims+1),pvalue=p)
  rownames(a) <- c("all G","G Clumped to maximise G","all chi-squares","chi-squared clumped to maximise G")
  a
}


"omnistat" <- function(xxx,nsims=1000,ncc) {
  if (nrow(xxx)!=2)
    stop("should pass matrix with 2 rows to clump")
  if (missing(ncc)) 
    res <- .C("omnibus1"
              ,as.double(t(xxx))
              ,as.integer(ncol(xxx))
              ,as.integer(nsims)
              ,teststat=as.double(numeric(2))
              ,over=as.integer(integer(2))
              ,PACKAGE="ijwtools"
              )
  else 
    res <- .C("omnibus2"
              ,as.double(t(xxx))
              ,as.integer(ncol(xxx))
              ,as.integer(nsims)
              ,as.integer(ncc)
              ,teststat=as.double(numeric(2))
              ,over=as.integer(integer(2))
              ,PACKAGE="ijwtools"
              )
#  p <- c(1.-pchisq(res$teststat[1],ncol(xxx)-1),NA,1.-pchisq(res$teststat[3],ncol(xxx)-1),NA)
  a <- cbind(testStatistic=res$teststat,simulatedp=(res$over+1)/(nsims+1))
  rownames(a) <- c("omnibus","Chi-squared")
  a
}
