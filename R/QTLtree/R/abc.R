

## this really does not work very well at the moment except for
## for a single parameter - with more than a single parameter then
## things start to go very wrong - and this is expected as at the
## moment the correction is applied to all the parameters independently
## Based on Baeumont et al 2003, Genetics
"abc" <- function(obsStats,simPars,simStats,delta=0.1
                ,sweep=TRUE,correct=TRUE,plot=FALSE,dis="rootsumsq")
  {
    ## the Epanechnikov kernel
    epan <- function(x,delta) {
      ret <- (1-(x/delta)^2)/delta;
      ret[x>delta] <- 0; 
      ret;
    }
    ## weighted kernel density estimates
    wdensity <- function(x,w,kernel="normal",pts=1024)
      {
        r <- range(x);
        xval <- seq(r[1]-0.1*(r[2]-r[1]),r[2]+0.1*(r[2]-r[1]),length=pts);
        x <- x[w>0];
        w <- w[w>0];
        y <- numeric(pts);
        bw=bw.nrd(x);
        for (i in 1:pts) {
          y[i] <- sum(w*dnorm(xval[i]-x,0,bw));
        }
        y <- y/sum(w);
        data.frame(x=xval,y=y);
      }
    if (sweep) { ## sweep out the standard deviation to standardise
      simStats <- rbind(obsStats,simStats)
      simStats <- scale(simStats)
      obsStats <- simStats[1,]
      simStats <- simStats[-1,]
    }
    ## calculate the individual distances between the simulated statistics and the
    ## test statistics  
    sdis <- sweep(simStats,2,obsStats)
    ## and summarise these distances across each different statistic measured
    if (dis=="rootsumsq")
      s2 <- sqrt(rowSums(sdis*sdis))
    else if (dis=="sumsq")
       s2 <- rowSums(sdis*sdis)
    else if (dis=="absdis")
      s2 <- rowSums(abs(sdis))
    else if (dis=="max")
      s2 <- apply(abs(sdis),1,max)
    else if (dis=="min")
      s2 <- apply(abs(sdis),1,min)
    
    ## what value of the distance is equivalent to get the correct delta
    kdelta <- quantile(s2,probs=c(delta));
    ## get the weighting using the Epanechnikov kernel
    w <- epan(s2,kdelta);
    sw <- sum(w)
    mw <- matrix(w,nrow=length(w),ncol=ncol(simPars))
  
    rawres <- colSums(simPars*mw)/sw
    if (correct) {
      np <- ncol(simPars)
      mod <- lm(simPars~sdis,weights=w)
      if (plot==TRUE)   {
        par(mfrow=c(2,np),mar=c(5,4,1,1))
        for (i in 1:np) {
          plot(simPars[w>0,i]~s2[w>0],ylab=colnames(simPars)[i],xlab="Distance",axes=FALSE)
          axis(1)
          axis(2)
          abline(mod$coef[1,i],mod$coef[2,i],col="blue")
          abline(rawres[i],0,col="red")
        }
      }
      coef <- matrix(mod$coef[1,],ncol=np,byrow=T,nrow=length(w))
      corrected <- simPars+coef-mod$fitted;
      
      if (plot==TRUE) {
        for (i in 1:np) {
          dcorr <- wdensity(corrected[,i],w)
          draw <- wdensity(simPars[,i],w)
          dprior <- density(simPars[,i])
          my <- max(c(dcorr$y,draw$y,dprior$y))
          plot(dcorr,ylim=c(0,my),col="blue",type="l",xlab=colnames(simPars)[i],ylab="density",axes=FALSE)
          axis(1)
          axis(2)
          lines(draw,col="red")
          lines(dprior,col="black")
        }
        ##      for (i in 1:np)
        ##        plot(corrected[w>0,i],w[w>0])
      }
      sw <- sum(w)
      w <- matrix(w,nrow=length(w),ncol=np)
      return(rbind(corrected=colSums(corrected*w)/sw,raw=rawres,prior=colMeans(simPars)))
    }
    return(rawres)
  }
