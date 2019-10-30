

## a function to write a data file out as a beagle phased
## file along with positions and traits files
write.beagle <- function(xx,filestem,traits=c("Case","Control"))
  {
    use <- xx$location %in% traits
    markers <- paste("SNP",xx$position,sep="_")
    yy <- cbind("M",markers, t(xx$haplotype[use,]))
    
    con <- file(paste(filestem,".phased",sep=""),"w")
    cat("#written from R function write.beagle.R\n",file=con)
    write.table(yy,row.names=FALSE,col.names=FALSE,quote=FALSE,file=con)
    close(con)
    
    con <- file(paste(filestem,".trait",sep=""),"w")
    cat("A trait ",file=con)
    cat(xx$location[use],file=con)
    close(con)
    #position file
    con <- file(paste(filestem,".position",sep=""),"w")
    cat("# written from R function write.beagle.R\n",file=con)
    write.table(cbind(markers,xx$position),row.names=FALSE,col.names=FALSE,quote=FALSE,file=con)
    close(con)
  }

pCaseControl <- function(b) {
  apply(b$haplotype,2,function(x) chisq.test(table(x,b$location))$p.value)

}
