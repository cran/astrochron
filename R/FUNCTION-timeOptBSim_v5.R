### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2026 Stephen R. Meyers
###
###########################################################################
### function timeOptBSim - (SRM: Feb. 24, 2026)
###########################################################################

timeOptBSim <- function(res,Nsim=1000,genplot=TRUE,check=TRUE,verbose=TRUE)
 {

  if(verbose) 
   {
    cat("\n----- TimeOptBSim: TimeOpt Astrochronologic Testing -----\n")
    cat("     Using the Method of Malinverno and Meyers (2024) \n\n")
   }

#######################################################################################
#### (1) data preparation
#######################################################################################
  Ndv = res['Ndv']
  dx = res['dx']
  mvopt = res['mvopt']
  tGa = res['tGa']
  umin = res['umin']
  umax = res['umax']
  Nu = res['Nu']
  kmin = res['kmin']
  kmax = res['kmax']
  Nk = res['Nk']
  phi1 = res['phi1']
  phi2 = res['phi2']
  rsqdatamap = res['rsqdatamap']
  rsqdatamapecc = res['rsqdatamapecc']
  rsqdatamapobl = res['rsqdatamapobl']
  rsqdatamappre = res['rsqdatamappre']
  detrend = res['detrend']

# calc zv, starting from zero
  zv0=(0:(Ndv-1))*dx  

  prior <- .getpriorpdfastropara(tGa,mvopt)         # parameters for priors, excluding u
# g's and s's are set to their prior mean value, which may differ from the
# mu parameter listed in Hoang et al. 2021
# here we are initializing mv to the prior means
  mv = prior$snmean                    # g's, s's, k (and u, to be added)
  priormeanastro = prior$snmean 
  priorsigmaastro = prior$sigma
# add element for u (set to zero)
  mv=append(mv,0)
# note that mv for u and k will be replaced during simluations. identify their indices in mv
  ik <- .getmvindex('k',mvopt)
  iu <- .getmvindex('u',mvopt)

# determine number of frequencies used for precession, obliquity and eccentricity fit
  eop <- .compeopfreq(mv,mvopt)
  ev=eop$ev
  ov=eop$ov
  pv=eop$pv
  nev= length(ev)
  nov= length(ov)
  npv= length(pv)

  je= 1:(nev*2)
  if(any(is.na(ov)))
   {
    nov=0
   }else{ 
    jo=(nev*2)+(1:(nov*2))
   }
  if (any(is.na(pv)))
   {
    npv=0
   }else{
    jp=((nev+nov)*2)+(1:(npv*2))
   }


#######################################################################################
#### (2) compliance checks
#######################################################################################
   if(check)
    {
# check that Nsim is at least a reasonable minimum value
     if (Nsim<1000)
      { 
       cat('**** WARNING: Nsim must be at least 1000. Nsim reset to 1000\n\n')
       Nsim=1000
      }
# check that the number of data points is at least a reasonable minimum value
     if (Ndv < 100) stop("\n**** ERROR: It is recommended that the number of stratigraphic data points should be at least 100. TERMINATING NOW!\n")
    }

# below is a mandatory check
     if (length(res)!=17) stop("\n**** ERROR: res not formatted correctly. TERMINATING NOW!\n")

#######################################################################################
#### (3) Monte Carlo significance calculation
#######################################################################################
# compute Nsim AR(2) time series in columns of matrix Dsim
   cat('=== Monte Carlo significance calculation:',Nsim,'\n')
   Dsim=matrix(0,nrow=Ndv,ncol=Nsim)
   dvsim=double(Ndv+2)
   for (j in 1:Nsim)
    {
# AR(2) time series in Dsim
     randev=rnorm(Ndv+2)
     dvsim[1]=randev[1]
     dvsim[2]=phi1*dvsim[1]+randev[2]
     for (i in 3:(Ndv+2)) dvsim[i]=phi1*dvsim[i-1]+phi2*dvsim[i-2]+randev[i]
     if (detrend) 
        {
          lm.1 <- lm(dvsim[3:(Ndv+2)] ~ seq(1:Ndv))
          dvsim = dvsim[3:(Ndv+2)]- (lm.1$coeff[2]*seq(1:Ndv) + lm.1$coeff[1])
        }      
     dvsim=scale(dvsim)[,1]   # zero mean, unit variance
     Dsim[,j]=dvsim
    }
   sumsqrdvsim=Ndv-1                   # sum of squares of normalized dvsim vectors

# set up progress display
  if(verbose)
   {
    cat('Significance calculation progress:\n')
    cat("0%       25%       50%       75%       100%\n")
    progress <- utils::txtProgressBar(min = 0, max = Nu, style = 1, width=43)
   }   
  
# compute max. R^2 of AR(2) simulated data in u, k grid
# sedimentation rate grid
  uv=seq(umin,umax,length.out=Nu)
# precession frequency grid  
  kv=seq(kmin,kmax,length.out=Nk)
# set up r2 vectors, initialize elements to zero  
  rsqdatasim=double(Nsim)
  rsqdatasimecc=double(Nsim)
  rsqdatasimobl=double(Nsim)
  rsqdatasimpre=double(Nsim)

# loop over sedimentation rate grid
  for (j in 1:Nu)
   {
    if(verbose) utils::setTxtProgressBar(progress, j) 
    tv=zv0/uv[j]                                               # calibrate meters to kyr
    mv[iu]=uv[j]
# loop over precession frequency grid
    for (i in 1:Nk)
     {
      mv[ik]=kv[i]
# compute matrix G
      G <- .compGmats(mv,mvopt,tv)$Gdv
#     H=(G'*G)\G';
      H=solve(t(G)%*%G,t(G))
# loop over simulations
      for (jsim in 1:Nsim)
       { 
        dvjsim=Dsim[,jsim]
# matrix H times dv gives fitted sine/cosine amplitudes afitv
        afitv=H%*%dvjsim
        thisrsq <- .comprsq(dvjsim,sumsqrdvsim,G,afitv)
        if (thisrsq>rsqdatasim[jsim])
         {
# record maximum r2 for this simulation (over all uv and kv)
          rsqdatasim[jsim]=thisrsq                             # R^2 for all astro. signals
          rsqdatasimecc[jsim] <- .comprsq(dvjsim,sumsqrdvsim,G,afitv,je)
          if (nov==0)
           {
            rsqdatasimobl[jsim]=0
           }else{
            rsqdatasimobl[jsim] <- .comprsq(dvjsim,sumsqrdvsim,G,afitv,jo)
           }
          if (npv==0)
           {              
            rsqdatasimpre[jsim]=0
           }else{
            rsqdatasimpre[jsim] <- .comprsq(dvjsim,sumsqrdvsim,G,afitv,jp)
           }
         }  
# end Nsim loop
       }
# end Nk look
     } 
# end Nu loop
   }   
# close progress display
  close(progress)  

# compute and output p-values
  cat('\n======= p-values  =======\n')
  vals <- .comppvalue(rsqdatamap,rsqdatasim)
  pval=vals[[1]]
  pvaln=vals[[2]]
  if(pvaln>0) cat('All astronomical cycles, p =',pval,'\n')
  if(pvaln==0) cat('All astronomical cycles, p <',pval,'\n')

  vals <- .comppvalue(rsqdatamapecc,rsqdatasimecc) 
  pvalecc=vals[[1]]
  pvaleccn=vals[[2]]
  if(pvaleccn>0) cat('Eccentricity only, p =',pvalecc,'\n')
  if(pvaleccn==0) cat('Eccentricity onl, p <',pvalecc,'\n')

  if (nov>0)
   {
    vals <- .comppvalue(rsqdatamapobl,rsqdatasimobl)
    pvalobl=vals[[1]]
    pvalobln=vals[[2]]
    if(pvalobln>0)cat('Obliquity only, p =',pvalobl,'\n')
    if(pvalobln==0) cat('Obliquity only, p <',pvalobl,'\n')
   }

  if (npv>0)
   {
    vals <- .comppvalue(rsqdatamappre,rsqdatasimpre)
    pvalpre=vals[[1]]
    pvalpren=vals[[2]]
    if(pvalpren>0) cat('Climatic precession only, p =',pvalpre,'\n')
    if(pvalpren==0) cat('Climatic precession only, p <',pvalpre,'\n')
   }


# ========================================================================
# ============================      plot      ============================
# ========================================================================
  if(genplot)
   {
    dev.new(height=7,width=4.5)
# presently, this only allows for omission of obliquity    
    if(nov>0) layoutA <- layout(matrix(c(1,2,3,4), 4, 1)) 
    if(nov == 0) layoutA <- layout(matrix(c(1,2,3), 3, 1)) 
# margins: bottom, left, top, right
    par(mar = c(3.5, 3.6, 2.1, 1))
#         title, text label, axis_line
    par(mgp = c(3, 0.5, 0))

# plot results for all astronomical cycles
    minmax=c(0,max(rsqdatasim,(1.2*rsqdatamap)))
    plot(density(rsqdatasim),col="gray",bty="L",xlim=minmax,type="l",xlab="",ylab="",main="")
    polygon(density(rsqdatasim),col="gray",border=NA)
    abline(v=rsqdatamap,lty=3,col="red",lwd=2)
    mtext(expression("Pearson r"^2),side=1,line=2,cex=0.8)
    mtext(c("All astronomical cycles"),side=3,line=0,font=4)
    mtext(bquote(r^2==.(round(rsqdatamap,digits=3)) ~ " "),side=3,at=rsqdatamap,line=-1.5,adj=1,cex=0.8,col="red")
    if(pvaln>0) mtext(bquote(p==.(round(pval,digits=3)) ~ " "),side=3,at=rsqdatamap,line=-2.7,adj=1,cex=0.8,col="red")
    if(pvaln==0) mtext(bquote(p<.(pval) ~ " "),side=3,at=rsqdatamap,line=-2.7,adj=1,cex=0.8,col="red")

    plot(density(rsqdatasimecc),col="gray",bty="L",xlim=minmax,type="l",xlab="",ylab="",main="")
    polygon(density(rsqdatasimecc),col="gray",border=NA)
    abline(v=rsqdatamapecc,lty=3,col="red",lwd=2)
    mtext(expression("Pearson r"^2),side=1,line=2,cex=0.8)
    mtext(c("Eccentricity only"),side=3,line=0,font=4)
    mtext(bquote(" " ~ r^2==.(round(rsqdatamapecc,digits=3))),side=3,at=rsqdatamapecc,line=-1.5,adj=0,cex=0.8,col="red")
    if(pvaleccn>0) mtext(bquote(" " ~ p==.(round(pvalecc,digits=3))),side=3,at=rsqdatamapecc,line=-2.7,adj=0,cex=0.8,col="red")
    if(pvaleccn==0) mtext(bquote(" " ~ p<.(pvalecc)),side=3,at=rsqdatamapecc,line=-2.7,adj=0,cex=0.8,col="red")

# cross out results for eccentricity if pval for all cycles > 0.1
    if(pval>0.1)
     {
      minmaxD=range(density(rsqdatasimecc)$y)
      lines(c(minmax),c(minmaxD),col="red",lwd=10)
      lines(c(minmax),c(minmaxD[2],minmaxD[1]),col="red",lwd=10) 
     }   

    if(nov>0)
     {
      plot(density(rsqdatasimobl),col="gray",bty="L",xlim=minmax,type="l",xlab="",ylab="",main="")
      polygon(density(rsqdatasimobl),col="gray",border=NA)
      abline(v=rsqdatamapobl,lty=3,col="red",lwd=2)
      mtext(expression("Pearson r"^2),side=1,line=2,cex=0.8)
      mtext(c("Obliquity only"),side=3,line=0,font=4)
      mtext(bquote(" " ~ r^2==.(round(rsqdatamapobl,digits=3))),side=3,at=rsqdatamapobl,line=-1.5,adj=0,cex=0.8,col="red")
      if(pvalobln>0) mtext(bquote(" " ~ p==.(round(pvalobl,digits=3))),side=3,at=rsqdatamapobl,line=-2.7,adj=0,cex=0.8,col="red")
      if(pvalobln==0) mtext(bquote(" " ~ p<.(pvalobl)),side=3,at=rsqdatamapobl,line=-2.7,adj=0,cex=0.8,col="red")

# cross out results for obliquity if pval for all cycles > 0.1
      if(pval>0.1)
       {
        minmaxD=range(density(rsqdatasimobl)$y)
        lines(c(minmax),c(minmaxD),col="red",lwd=10)
        lines(c(minmax),c(minmaxD[2],minmaxD[1]),col="red",lwd=10) 
       }   
     } 

    plot(density(rsqdatasimpre),col="gray",bty="L",xlim=minmax,type="l",xlab="",ylab="",main="")
    polygon(density(rsqdatasimpre),col="gray",border=NA)
    abline(v=rsqdatamappre,lty=3,col="red",lwd=2)
    mtext(expression("Pearson r"^2),side=1,line=2,cex=0.8)
    mtext(c("Climatic precession only"),side=3,line=0,font=4)
    mtext(bquote(" " ~ r^2==.(round(rsqdatamappre,digits=3))),side=3,at=rsqdatamappre,line=-1.5,adj=0,cex=0.8,col="red")
    if(pvalpren>0) mtext(bquote(" " ~ p==.(round(pvalpre,digits=3))),side=3,at=rsqdatamappre,line=-2.7,adj=0,cex=0.8,col="red")
    if(pvalpren==0) mtext(bquote(" " ~ p<.(pvalpre)),side=3,at=rsqdatamappre,line=-2.7,adj=0,cex=0.8,col="red")

# cross out results for climatic precession if pval for all cycles > 0.1
    if(pval>0.1)
     {
      minmaxD=range(density(rsqdatasimpre)$y)
      lines(c(minmax),c(minmaxD),col="red",lwd=10)
      lines(c(minmax),c(minmaxD[2],minmaxD[1]),col="red",lwd=10) 
     }

# end genplot
   }
# end function timeOptBSim
 }
 