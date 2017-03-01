
combineLRplot <- function(stem,wh=1,cols=c("green","blue"),cutoff=1000,addbad=0,...)
{
  Left <- get(paste(stem,"L",sep=""))
  Left <- Left[-nrow(Left),]
  Right <- get(paste(stem,"R",sep=""))[-1,]
  between <- (Right$position+Left$position)/2
  column <- 2*wh
  wL <- which(Left[,column+1]>=cutoff)
  wR <- which(Right[,column+1]>=cutoff)
  
  mx <- round(between[length(between)]/1000000,0)

  plot(between,Left[,column],type="l",ylab="log p-value",xlab="position (MB)",axes=F,col=cols[1],...)
  axis(1,at=seq(0,mx,1)*1000000,labels=seq(0,mx,1))
  axis(2)
  lines(between,Right[,column],col=cols[2])
  
  rug(between[wL],col=cols[1],side=3)
  rug(between[wR],col=cols[2],side=3,ticksize=-0.03)

  both <-  between[wR[wR %in% wL]]
  if (length(both)>0)
    rug(between[wR[wR %in% wL]],side=1,col="cyan")
  

  if (addbad>0) {
    ch <- getStatistics(addbad)
    badpos <- ch$pos[ch$good_clustering==0]
    points(badpos,rep(min(Left[,column]),length(badpos)),col="red",pch=2)
    }
  r <- range(both)
  #if (r[2]-r[1]<1E6) cat(paste("chr",addbad,":",round(r[1],0),"..",round(r[2],2),sep=""),"\n")
  cl <- getclusters(both)
  lapply(cl,function(x)  cat(paste("chr",addbad,":",round(x[1],0),"..",round(x[2],0),sep=""),"\n"))
  both
}

checkdata <- function(pattern="c*_*_58C,HT.phased")
  {
    tfile <- tempfile()
    system(paste("find . -name ",pattern," > ",tfile,sep=""))
    con <- file(tfile,"r")
    b <- readLines(con)
    #print(b)
    unlink(tfile)
    close(con)
    res <- strsplit(matrix(unlist(strsplit(b,"/")),ncol=3,byrow=T)[,2],"Chr")
    table(as.numeric(matrix(unlist(res),ncol=2,byrow=T)[,2]))
  }

testtable <- function(a,reps=1000,k=4)
  {
    wh <- seq(3,by=2,length=k)
    tb <- which(a[,wh]==reps,arr.ind=TRUE)
    table(a[tb[,1],1],tb[,2])
  }

getdata <- function(chr=1,stem="/home/dijw/IHG/WTCCC/Experiment2/",ex="",first=0,last=250,step=10,k=10,reps=1000,stat='P',seed=1,TreeLength=10,pops="58C,HT",direction="",skip=0)
{
  myreadtable <- function(x,columnnames) {
    if (file.exists(x)) {
      return(read.table(x,col.names=columnnames,skip=skip))
    } else {
      a = data.frame(matrix(ncol=length(colnames),nrow=0))
      colnames(a) <- columnnames
      return(a)
    }
  }
  colnames <- c("position",paste(c("S","p"),rep(1:k,rep(2,k)),sep=""),"indexmin","posmin","posmax","nleaves")
 
  allstem <- paste(stem,"Chr",chr,"/",ex,"c",chr,sep="")
 
  names <-  paste(allstem,seq(first,last-step,step),seq(first+step,last,step),paste(pops,"results",sep="."),sep="_")
  names <- paste(names,paste(TreeLength,stat,reps,seed,sep="."),sep=".")
  if (direction!="")
    names <- paste(names,direction,sep=".")
  print(names)
  xx <- lapply(names, myreadtable,columnnames=colnames)
  do.call("rbind",xx)
}
## a very simple estimate of the clusters
getclusters <- function(positions,minss=8E-2)
  {
    if (length(positions)==1) return(list(c(positions,positions+1)))
    if (length(positions)==2) {
      if (diff(positions)>10000)
        return(list(c(positions[1],positions[1]+1),c(positions[2],positions[2]+1)))
      else return(list(positions))
    }
    require(cluster)
    pos <- positions/diff(range(positions))
    d <- dist(pos)
    ss <- sum(pos*pos)
    if (ss < minss) return(list(range(positions)))
    maxk <- min(40,length(positions)-1)
    if (maxk>=2){
      for (k in 2:maxk) {
        km <- kmeans(d,k)
        wss <- km$withinss
        if (max(wss)<minss) {   ## finished
          return(tapply(positions,km$cluster,range)) 
        }  
      }
    }
    
    
}
      
  

plotdata <- function(xx,wh=1,add=FALSE,side=3,rugcol="red",...)
  {
    cl <- 2*wh
    w <- which(xx[,cl+1]==1000)
    if (add==FALSE)
      plot(xx[,1],xx[,cl],type="l",ylab="log p-value",xlab="position (MB)",axes=F,...)
    else lines(xx[,1],xx[,cl],...)
    axis(2)
    mx <- round(xx[nrow(xx),1]/10000000,0)*10
    print(mx)
    axis(1,at=seq(0,mx,10)*1000000,labels=seq(0,mx,10))
    rug(xx[w,1],col=rugcol,side=side)
  }

a2 <- getdata("/home/dijw/IHG/WTCCC/Experiment2/Chr2/c2")
a1 <-  getdata(last=250)


GeneLocator(1)

GeneLocator <- function(chr=2) {
  require(RSNPper)
  GetGenesb <- function(x) {
    GenesFromPos <- function(x) {
      r <- itemsInRange("genes",paste("chr",chr,sep="")
                        ,as.character(as.integer(x)-10)
                        ,as.character(as.integer(x)+10))
      if (length(r)>0) return(list(name=r[[1]]["NAME"],position=x))
    }
    b <- lapply(x,GenesFromPos)
    matrix(unlist(b),ncol=2,byrow=T) 
  }
  res <- GetGenesb(locator()$x)
}

a1 <- getdata(1)
a1T <- getdata(1,ex="aaa_",stat="T",first=210,last=230,k=4)
a1Tl <- getdata(1,ex="aaa_",stat="T",first=210,last=230,k=4,TreeLength=40)
a2  <- getdata(2)
a2c <- getdata(2,pops="58C,NBS")
a2T <- getdata(2,stat="T",k=4,TreeLength=50)
a2Tc <- getdata(2,stat="T",k=4,TreeLength=50,pops="58C,NBS")
a2tc <- getdata(2,stat="t",k=6,TreeLength=50,pops="58C,NBS")

a2t <- getdata(2,stat="t",k=6,TreeLength=50)
a2Tb <- getdata(2,stat="T",k=4,TreeLength=50,pops="NBS,HT")

a2t.100 <- getdata(2,stat="t",k=6,TreeLength=100)
a2tc.100 <- getdata(2,stat="t",k=6,TreeLength=100,pops="58C,NBS")
 

a3 <- getdata(3,last=200)
a3L <-  getdata(3,stat="P",dir="L",last=90,skip=3)
a3R <- getdata(3,stat="P",dir="R",last=90,skip=3)
a4 <- getdata(4,last=200)
a4L <-  getdata(4,stat="P",dir="L",last=140,skip=3)
a4R <- getdata(4,stat="P",dir="R",last=140,skip=3)
a5 <- getdata(5,last=190)
a6 <- getdata(6,last=180)
a7 <- getdata(7,last=160)
a8 <- getdata(8,last=150)
a8L <-  getdata(8,stat="P",dir="L",last=150)
a8R <- getdata(8,stat="P",dir="R",last=150)
a9 <- getdata(9,last=140)
a10 <- getdata(10,last=140)
a11 <- getdata(11,last=140)
a12 <- getdata(12,last=140)
a13 <- getdata(13,last=120)
a13L <-  getdata(13,stat="P",dir="L",last=120)
a13R <- getdata(13,stat="P",dir="R",last=120)
a14 <- getdata(14,last=110)
a14L <-  getdata(14,stat="P",dir="L",last=110,skip=3)
a14R <- getdata(14,stat="P",dir="R",last=110,skip=3)
a15 <- getdata(15,last=110)
a15L <-  getdata(15,stat="P",dir="L",last=110)
a15R <- getdata(15,stat="P",dir="R",last=110)
a16 <- getdata(16,last=90)
a16L <-  getdata(16,stat="P",TreeLength=10,dir="L",skip=3)
a16R <- getdata(16,stat="P",TreeLength=10,dir="R",skip=3)
a16t.50 <- getdata(16,stat="t",k=6,TreeLength=50)
a16t.50L <- getdata(16,stat="t",k=6,TreeLength=50,direction="L")
a16t.50R <- getdata(16,stat="t",k=6,TreeLength=50,direction="R")

a17 <- getdata(17,last=80)
a17L <-  getdata(17,stat="P",dir="L",last=120)
a17R <- getdata(17,stat="P",dir="R",last=120)
a18 <- getdata(18,last=80)
a18L <-  getdata(19,stat="P",dir="L",last=80,skip=3)
a18R <- getdata(19,stat="P",dir="R",last=80,skip=3)
a19 <- getdata(19,last=70)
a19L <-  getdata(19,stat="P",dir="L",last=110,skip=3)
a19R <- getdata(19,stat="P",dir="R",last=110,skip=3)
a20 <- getdata(20,last=70)
a20L <-  getdata(20,stat="P",dir="L",last=110,skip=3)
a20R <- getdata(20,stat="P",dir="R",last=110,skip=3)
a21 <- getdata(21,last=50)
a21L <-  getdata(21,stat="P",dir="L",last=110,skip=3)
a21R <- getdata(21,stat="P",dir="R",last=110)
a22 <- getdata(22,last=50)
a22L <-  getdata(22,stat="P",dir="L",last=110,skip=3)
a22R <- getdata(22,stat="P",dir="R",last=110,skip=3)


postscript(file="plot1.eps",height=7,width=11)
par(mfrow=c(3,1))
plotdata(a1,main="Chromosome 1")
plotdata(a2,main="Chromosome 2")
plotdata(a3,main="Chromosome 3")
dev.off()

postscript(file="plot2.eps",height=7,width=11)
par(mfrow=c(3,1),mar=c(4,4,1,1))
plotdata(a4,main="Chromosome 4")
plotdata(a5,main="Chromosome 5")
plotdata(a6,main="Chromosome 6")
dev.off()

postscript(file="plot3.eps",height=7,width=11)
par(mfrow=c(3,1),mar=c(4,4,1,1))
plotdata(a7,main="Chromosome 7")
plotdata(a8,main="Chromosome 8")
plotdata(a9,main="Chromosome 9")
dev.off()

postscript(file="plot4.eps",height=7,width=11)
#pdf(file="plot4.pdf",paper="a4")
par(mfrow=c(3,1),mar=c(4,4,1,1))
plotdata(a10 ,main="Chromosome 10")
plotdata(a11,main="Chromosome 11")
plotdata(a12,main="Chromosome 12")
dev.off()

postscript(file="plot5.eps",height=7,width=11)
#pdf(file="plot5.pdf")
par(mfrow=c(3,2))
plotdata(a13,main="Chromosome 13")
plotdata(a14,main="Chromosome 14")
plotdata(a15,main="Chromosome 15")
plotdata(a16,main="Chromosome 16")
plotdata(a17,main="Chromosome 17")
plotdata(a18,main="Chromosome 18")
dev.off()

postscript(file="plot6.eps",height=7,width=11)
par(mfrow=c(3,2))
plotdata(a19,main="Chromosome 19")
plotdata(a20,main="Chromosome 20")
plotdata(a21,main="Chromosome 21")
plotdata(a22,main="Chromosome 22")
dev.off()
