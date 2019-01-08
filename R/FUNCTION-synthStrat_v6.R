### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2018 Stephen R. Meyers
###
###########################################################################
### synthStrat: model stratigraphy, create a visual representation
###          (SRM: August 1-2, 2018; August 4, 2018)
###
###########################################################################


synthStrat <- function(signal=NULL,nfacies=4,clip=T,flip=F,fmax=0.1,output=F,genplot=2,verbose=T)
 {

  if(verbose) cat("\n----- Generating synthetic stratigraphy -----\n")

  if(!is.null(signal)) signal=data.frame(signal)
  if(is.null(signal)) signal=cycles(start=0,end=1000,dt=1,genplot=F,verbose=F)

# standardize  
  signal[2] = (signal[2] - mean(signal[,2]))/sd(signal[,2])
  if(flip) signal[2]=-1*signal[2]
  
# create stratigraphy using nfacies
  strat=signal
  minSig=min(signal[2])
  maxSig=max(signal[2])

# truncate at mean  
  if(clip) 
   {
     inc=maxSig/nfacies
     for(i in nfacies:1) strat[(signal[2]<=inc*i),2]=i
   }  
# do not truncate  
  if(!clip) 
   {
     inc=(maxSig-minSig)/nfacies
     for(i in nfacies:1) strat[(signal[2]<=(minSig+inc*i)),2]=i
# ensure max signal was assigned correctly (allowing for round off error)
     strat[which.max(signal[,2]),2]=nfacies    
   }  

# initialize plot
  if(genplot==1)  
   {
    dev.new(height=7,width=6,title="synthStrat")
    par(mfrow=c(1,3))
   } 
  if(genplot==2)
   {
    dev.new(height=7,width=10,title="synthStrat")
    m <- cbind(c(1,1),c(2,2),c(3,3),c(4,5),c(4,5))
    layout(m)
   }
  par(mar = c(5.1, 4.1, 2.1, 0.5))
  plot(signal[,2],signal[,1],ylim=c(max(signal[1]),min(signal[1])),xlab="",ylab="",bty="n",type="l",lwd=2)
  mtext("Forcing",side=3,font=2)
  mtext("Normalized Forcing Signal",side=1,line=2.5)
  mtext("Depth (cm) or Time (ka)",side=2,line=2.5)
  if(clip) 
   {
     abline(v=0,lty=3,col="red")
     for(i in nfacies:1) abline(v=inc*i,lty=3,col="red")
   } 
  if(!clip) 
   {
     for(i in nfacies:1) abline(v=minSig+inc*i,lty=3,col="red")
     abline(v=minSig,lty=3,col="red")
   }  
#  plot(strat[,2],strat[,1],xlim=c(0,nfacies),ylim=c(max(strat[1]),min(strat[1])),xlab="",ylab="",bty="n")
  plot(strat[,2],strat[,1],xlim=c(0,nfacies),ylim=c(max(strat[1]),min(strat[1])),xlab="",ylab="",bty="n",type="l")
  mtext("Facies",side=3,font=2)
  mtext("Sedimentary Facies",side=1,line=2.5)
  plot(-100,-100,xlim=c(0,1),ylim=c(max(strat[1]),min(strat[1])),xlab="",ylab="",bty="n",xaxt="n",cex=0)
  mtext("Stratigraphy",side=3,font=2)
# add polygons
# define colors for fill
  fill=gray.colors(n=nfacies,start=0.9,end=0.1,gamma=2.2)
#  polygon(c(-1,1,1,-1),c(min(strat[1]),min(strat[1]),max(strat[1]),max(strat[1])),col=fill[1],border=NA)
# find transition points, values != 0 identify transitions
  trans=diff(strat[,2])
  transPts=which(trans!=0)
  beds=length(transPts)+1 
# plot first bed
  polygon(c(-1,1,1,-1),c(min(strat[1]),min(strat[1]),strat[transPts[1],1],strat[transPts[1],1]),col=fill[strat[transPts[1],2]],border=fill[strat[transPts[1],2]])
  for(i in 1:(beds-1)) polygon(c(-1,1,1,-1),c(strat[transPts[i],1],strat[transPts[i],1],strat[transPts[i+1],1],strat[transPts[i+1],1]),col=fill[strat[transPts[i+1],2]],border=fill[strat[transPts[i+1],2]])
# plot last bed
  polygon(c(-1,1,1,-1),c(strat[transPts[beds-1],1],strat[transPts[beds-1],1],max(strat[1]),max(strat[1])),col=fill[strat[transPts[beds-1]+1,2]],border=fill[strat[transPts[beds-1]+1,2]])

  if(genplot==2)
   {
     forcingSpec=periodogram(signal,detrend=T,output=1,genplot=F,verbose=F)
     faciesSpec=periodogram(strat,detrend=T,output=1,genplot=F,verbose=F)
     par(mar = c(5.1, 5.1, 4.1, 1.5))
     plot(forcingSpec[,1],forcingSpec[,2],xlim=c(0,fmax),type="l",xlab="",ylab="")
     mtext("Forcing Spectrum",side=3,line=2,font=2)
     mtext("Frequency",side=1,line=2.5)
     mtext("Amplitude",side=2,line=2.5)
     plot(faciesSpec[,1],faciesSpec[,2],xlim=c(0,fmax),type="l",xlab="",ylab="")
     mtext("Facies Spectrum",side=3,line=2,font=2)
     mtext("Frequency",side=1,line=2.5)
     mtext("Amplitude",side=2,line=2.5)
   }  
   
  if(output) return(strat)
# end function synthStrat  
 }