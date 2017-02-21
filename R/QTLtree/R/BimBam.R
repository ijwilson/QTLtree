writeBimBam <- function(aa,stem="out1",rep=1)
  ## function to write BimBam input file from my own input
  {
    res <- aa$replicate1
    
    con <- file(paste(stem,"geno",sep="."),"w")
    con2 <-  file(paste(stem,"position",sep="."),"w")
    nsites <- ncol(res$haplotypes)
    ndip <- nrow(res$haplotypes)/2
    
    cat(ndip,"\n",nsites,"\n",file=con,sep="")

    snpnames <- paste("rs",1:nsites,sep="")
    
    for (i in 1:nsites) {
      aaa <- matrix(res$haplotypes[,i],ncol=2,byrow=T)
      cat(snpnames[i],",",file=con)
      cat(c("AA","AT","TT")[rowSums(aaa)+1],sep=",",file=con)
      cat("\n",file=con)

      cat(snpnames[i],res$position[i],sep=",",file=con2)
      cat("\n",file=con2)
    }

    close(con)
    close(con2)
    ## now write the phenotype file
    pheno <- gl(2,ndip/2,labels=c(0,1))
    write(paste(pheno),file=paste(stem,"pheno",sep="."))
  }

writeBimBamA <- function(res,stem="out1",rep=1)
  ## function to write BimBam input file from my own input
  {
    con <- file(paste(stem,"geno",sep="."),"w")
    con2 <-  file(paste(stem,"position",sep="."),"w")
    nsites <- ncol(res$haplotype)
    ndip <- nrow(res$haplotype)/2
    
    cat(ndip,"\n",nsites,"\n",file=con,sep="")

    snpnames <- paste("rs",1:nsites,sep="")
    
    for (i in 1:nsites) {
      aaa <- matrix(res$haplotype[,i],ncol=2,byrow=T)
      cat(snpnames[i],",",file=con)
      cat(c("AA","AT","TT","?A","?T","6","??")[rowSums(aaa)+1],sep=",",file=con)
      cat("\n",file=con)

      cat(snpnames[i],res$position[i],sep=",",file=con2)
      cat("\n",file=con2)
    }

    close(con)
    close(con2)
    ## now write the phenotype file
    pheno <- gl(2,ndip/2,labels=c(0,1))
    write(paste(pheno),file=paste(stem,"pheno",sep="."))
  }

AddMissing <- function(a,prop=0.05)
  {
    ncol=ncol(a$haplotype)
    nrow=nrow(a$haplotype)
    miss <- matrix(rbinom(ncol*nrow,1,prop),ncol=ncol)
    a$haplotype[miss==1] <- 3
    a
  }
