`readTNT` <-
function(con,skip=0,close.con=FALSE,what=double(0))
  {
    if (is.character(con)) {
      con <- file(con,"r");
    }
    ## cat("skipping ... ",readLines(con,skip),"\n");
    a <- scan(con,n=1,skip=skip,quiet=TRUE);
    res <- scan(con,n=a,what=what,quiet=TRUE);
    ## print(class(res))
    if (close.con) close(con)
    res
  }


"writeTNT" <- function(data,con)
  {
    if (is.character(con)) 
      con <- file(con,"w")
    
    if (is.matrix(data)) {
      cat(nrow(data)," ",ncol(data),"\n",file=con)
      write.table(data,row.names=FALSE,col.names=FALSE,file=con)
    } else {
      cat(length(data),"\n",file=con)
      cat(data,"\n",file=con)
    }
  }

`readTNT2D` <-
function(con,skip=0,close.con=FALSE)
  {
    if (is.character(con)) {
      con <- file(con,"r");
    }
    a <- scan(con,n=2,skip=skip,quiet=TRUE);
      cl <-a[2]
    rw <- a[1];
    res <- matrix(scan(con,n=cl*rw,quiet=TRUE),ncol=cl,nrow=rw,byrow=TRUE)
    if (close.con) close(con)
    res;
  }

