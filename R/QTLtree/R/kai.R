testDip <- function(cases,penetrance=c(0,0.5,1.0),theta=40,seed)
  {
    require(ijwtools)
    if (missing(seed))
      a <- mutsim(cases/2,penetrance,theta=theta)
    else  a <- mutsim(cases/2,penetrance,theta=theta,seed=seed)

    haps <- apply(a$haplotype,1,paste,collapse="")
    omnistat(table(gl(2,cases),haps))[1,2]
  }


testHap <- function(cases,penetrance=c(0.0,1.0),theta=40,seed)
  {
    require(ijwtools)
 
    if (missing(seed))
      a <- mutsimhap(cases,200,penetrance[2],penetrance[1],theta=theta)
    else
      a <- mutsimhap(cases,200,penetrance[2],penetrance[1],theta=theta,seed=seed)
    a
  }


mutsimhap <- function(cases,
                   trysites=200,penetrance=0.5,
                   sporadic=0.36,
                   sites=20,
                   minfreq=0.45,
                   maxfreq=0.55,
                   minmarkerfreq=0.005,
                   theta=40,
                   alpha=10,
                      highmutsites=200,
                   seed) {
  trysites=200
  tmpfile=tempfile()
  if (missing(seed)) seed <-  floor(runif(1)*1E7)
  command <- paste("mutsim --s1=",trysites," --s2=",highmutsites," --s3=",trysites," --ss=5000 --alpha=",alpha," --t=",theta," --seed=",seed
    ," | enrich --loud --edge=",trysites," --penetrance=",penetrance," --sporadic=",sporadic," --ss1=",cases," --ss2="
    ,cases," --min=",minfreq," --max=",maxfreq," | stripmiddle ",sites," ",trysites," ",minmarkerfreq," > ",tmpfile,sep="")
  repeat {
    ret <- system(command)
    if (ret!=0) break
  }
  
 
    
  a=matrix(scan(tmpfile,quiet=TRUE),ncol=sites,byrow=TRUE)
  
  haps <- apply(a,1,paste,collapse="")

  omnistat(table(gl(2,cases),haps))[1,2]
}

