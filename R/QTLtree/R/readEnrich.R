`readEnrich` <-
function(filename)
  {
    conn <- file(filename,"r")
    
    rep <- 1
    results=list()
    while (skipToTag(conn,"replicate")) {
      results[[paste("replicate",rep,sep="")]] <- simplexml(conn,numbertags=c("use","diseasepos"),vectags=c("position","locations","diseaseSNP","labels"),mattags=c("haplotypes"),opening="replicate")
      class(results[[paste("replicate",rep,sep="")]]) <- "GenomicHaplotype"
      rep <- rep+1
    }
    close(conn)
    results
  }

