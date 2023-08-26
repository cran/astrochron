### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2023 Stephen R. Meyers
###
###########################################################################
### function bicoherence : Calculate bispectrum and bicoherence using WOSA method,
###                       as detailed in Choudhury, Shah, and Thornhill (2008).
###                       This version uses Kim and Powers (1979) bicoherence normalization,
###                       which is bounded [0-1]. (SRM: Nov 16-30, 2012; March 19, 2013; 
###                       Feb. 13-17, 2022; Mar. 8, 2023; July 22-24, 2023; August 23, 2023)
###########################################################################


bicoherence <- function (dat,overlap=50,segments=8,CF=0,CL=95,padfac=2,demean=T,detrend=T,taper=T,maxF=Nyq,output=0,genplot=T,color=1,id=NULL,logpwr=F,logbis=F,check=T,verbose=T)
{

# options below are hardwired
# normfft : normalize fft results by segpts (T or F)
 normfft=T

 if(verbose) cat("CALCULATE BISPECTRUM AND BICOHERENCE USING THE WOSA METHOD\n")

 d <- data.frame(dat)
 npts <- length(d[,1]) 
 dt = abs(d[2,1]-d[1,1])
 
# error checking
 if(check)
   {
     dtest <- d[2:npts,1]-d[1:(npts-1),1] 
     epsm=1e-9
     if( (max(dtest)-min(dtest)) > epsm ) 
       {
         cat("\n**** ERROR: sampling interval is not uniform.\n")
         stop("**** TERMINATING NOW!")
       }
# padfac must be >= 1
     if(padfac < 1) padfac=1  
     if(overlap>50)
       {
         if(verbose) cat("* overlap cannot be greater than 50. Resetting to 50\n")
         overlap=50
       }
   }

# Nyquist
 Nyq <- 1/(2*dt)
# number of data points per segment
 segpts=floor(npts/(segments-(segments-1)*overlap/100))
# Rayleigh Frequency
 Ray <- 1/(segpts*dt)
# guarantee an even number of data points after padding, so Nyquist exists.   
 if((segpts*padfac)%%2 != 0) padpts= (segpts*padfac) + 1
 if((segpts*padfac)%%2 == 0) padpts= (segpts*padfac)
# frequency grid, post-padding
 df <- 1/(padpts*dt)
# set index vector for frequencies (to total padded length)
 i <- seq(1,padpts,by=1)
# make frequency vector. neg frequences not assigned correctly here (that's OK, will ignore).
 freq <- df*(i-1)
# number of frequencies from f(0) to Nyquist (post-padding)
 zero2Nyq = 0.5*(padpts) + 1

 if(verbose) 
  {
    cat("* Number of data points=", npts,"\n")
    cat("* Sample interval=", dt,"\n")
    cat("* Nyquist frequency=", Nyq,"\n")
    cat("\n* Number of data points per segment=",segpts,"\n")
    cat("* Rayleigh frequency=",Ray,"\n")
  }

# set up array to hold fft results, for selected number of segments
 ft<- rep(NA,padpts*segments)
 dim(ft) <- c(padpts,segments)

 if(verbose) cat("\n* Calculating WOSA Fourier Coefficients and Power Spectrum  \n")
# loop over time series segments
 ihold=1
 for (i in 1:segments)
  {
# shuffle in data
#   if(verbose) cat("  - Processing segment",i,":",ihold,",",(ihold+segpts-1),"\n")
   dd <- data.frame(cbind( d[ihold:(ihold+segpts-1),1], d[ihold:(ihold+segpts-1),2] ))
# remove mean and linear trend if desired
   if (demean) 
    { 
      dave <- colMeans(dd[2])
      dd[2] <- dd[2] - dave
#      if(verbose) cat("     Mean value removed=",dave,"\n")
    }
# use least-squares fit to remove linear trend
   if (detrend) 
    {
      lm.0 <- lm(dd[,2] ~ dd[,1])
      dd[2] <- dd[2] - (lm.0$coeff[2]*dd[1] + lm.0$coeff[1])
#      if(verbose) cat("     Linear trend removed: m=",lm.0$coeff[2],"b=",lm.0$coeff[1],"\n")
    }
# apply Hanning taper if desired
   if (taper)
    {
      w=.5*(1-cos(2*pi*(1:segpts)/(segpts+1)))
      dd[2] = dd[2] * w   
    } 
# pad and also convert from data frame to numeric
   pad <- as.numeric(dd[,2])
   if (padpts>segpts) pad <- append(pad, rep(0, (padpts - segpts)))
# take fft
   ft[,i] <- fft(pad)
   if(normfft) ft[,i] <- ft[,i]/segpts
# now set ihold for beginning of next segment
   ihold = ihold + floor(segpts*(100-overlap)/100)
# end WOSA fft segments loop
 }

# calculate power for each segment and average for WOSA
 pwr <- rowMeans(Mod(ft)^2)


##########################################################
### BISPECTRUM 
##########################################################

 if(verbose) cat("\n* Calculating Bispectrum \n")
# number of frequencies for bispectrum and bicoherence
 nfreq = floor(maxF/df)

# set up nfreq x nfreq matrix for complex results, initialize to zero   
 B = complex(nfreq*nfreq)
 dim(B) <- c(nfreq,nfreq)

 for (ii in 1:segments)
  {
    for (i in 1:nfreq) 
     {
       B[i,1:nfreq] = B[i,1:nfreq] + ft[i,ii]*ft[1:nfreq,ii]*Conj(ft[i+(1:nfreq)-1,ii])
     }
  }    
 B = B/segments
# squared magnitude
 Bmag2 = Mod(B)^2
 if(logbis) Bmag3 = log(Bmag2)
 if(!logbis) Bmag3 = Bmag2

##########################################################
### BICOHERENCE
##########################################################

 if(verbose) cat("* Calculating Bicoherence: Kim and Powers method \n")   
# see eq 3.27 of Choudhury et al. (2008)
# set up nfreq x nfreq matrix for real results, initialize to zero   
 f1f2 = double(nfreq*nfreq)
 dim(f1f2) <- c(nfreq,nfreq)
 f3 = double(nfreq*nfreq)
 dim(f3) <- c(nfreq,nfreq)
 for (ii in 1:segments)
  {
    for (i in 1:nfreq)
     {
          f1f2[i,1:nfreq] = f1f2[i,1:nfreq] + Mod( ft[i,ii]*ft[1:nfreq,ii] )^2
          f3[i,1:nfreq] = f3[i,1:nfreq] + Mod( ft[i+(1:nfreq)-1,ii] )^2
     }
  }    
 f1f2 = f1f2/segments
 f3 = f3/segments

# set up nfreq x nfreq matrix for real results, initialize to zero  
 den = double(nfreq*nfreq)
 dim(den) <- c(nfreq,nfreq)    
# calculate denomenator of bicoherence without conditioning factor
 for (i in 1:nfreq)
  {
    den[i,1:nfreq] = (f1f2[i,1:nfreq]*f3[i,1:nfreq])
  }  

 if(CF==0)
  {     
 if(verbose) cat("\n* Calculating Conditioning Factor \n")
# calculate 75th percentile for Magnitude-squared bicoherence denominator, as suggested on pg. 80 of Choudhury et al., 2008 
# again, here we are evaluating both sides of the symmetric matrix, as well as the diagonal.
     p95 = sort(den)[floor(length(den)*.95)]
     p90 = sort(den)[floor(length(den)*.9)]
     p75 = sort(den)[floor(length(den)*.75)]
     p50 = sort(den)[floor(length(den)*.5)]
     p25 = sort(den)[floor(length(den)*.25)]
     if(verbose) 
      {
        cat("* 95th percentile for Magnitude-squared bicoherence denominator is at",p95,"\n")
        cat("* 90th percentile for Magnitude-squared bicoherence denominator is at",p90,"\n")
        cat("* 75th percentile for Magnitude-squared bicoherence denominator is at",p75,"\n")
        cat("* 50th percentile for Magnitude-squared bicoherence denominator is at",p50,"\n")
        cat("* 25th percentile for Magnitude-squared bicoherence denominator is at",p25,"\n")
      }  
# automatically assign to CF
        if(verbose) cat("* Automatically assigning 75th percentile value as conditioning factor for bicoherence denominator \n\n")      
        CF = p75
# end CF==0   
  }

 if(genplot)
  {
### plot the denominator of the squared bicoherence (see Choudhury et al., 2008) 
     dev.new(title=paste("Denominator of Magnitude-squared Bicoherence"),height=8,width=7)
### plot the WOSA Power results
    par(mfrow=c(3,1))
    if(logpwr) plot(freq[1:nfreq],log(pwr[1:nfreq]), type="l",col="red", ylab="Log(Power)", xlim=c(0,maxF),xlab="Frequency", main="WOSA Power Spectrum")
    if(!logpwr) plot(freq[1:nfreq],pwr[1:nfreq], type="l",col="red", ylab="Power", xlim=c(0,maxF),xlab="Frequency", main="WOSA Power Spectrum")
    maxd=max(den)
    mind=min(den)
# here we are evaluating both sides of the symmetric matrix, as well as the diagonal
    plot(freq[1:nfreq],den[,1],type="l",xlab="Frequency (f1)", ylab="Bicoherence Denominator",ylim=c(mind,maxd),main="Denominator of Magnitude-squared Bicoherence")
    lines(rep(c(NA,freq[1:nfreq]),nfreq-1),rbind(NA,den[,2:nfreq]))
    abline(h=CF,lwd=2,col="red")
    plot(freq[1:nfreq],log(den[,1]),type="l",xlab="Frequency (f1)", ylab="Log(Bicoherence Denominator)",ylim=c(log(mind),log(maxd)),main="Denominator of Magnitude-squared Bicoherence")
    lines(rep(c(NA,freq[1:nfreq]),nfreq-1),rbind(NA,log(den[,2:nfreq])))
    abline(h=log(CF),lwd=2,col="red")
  }
        
# NOW ADD CONDITIONING FACTOR AND CALCULATE BICOHERENCE
 den = den + CF
# set up nfreq x nfreq matrix for real results, initialize to zero   
 bic2 = double(nfreq*nfreq)
 dim(bic2) <- c(nfreq,nfreq)  
 bic2 = Bmag2/den

### evaluate signficance of bicoherence magnitude at each individual bifrequency, using method outlined on pg. 83 of Choudhury et al. (2008)
 if(verbose) cat("* Calculating Magnitude-squared Bicoherence Confidence Levels \n")
 bic2CL = double(nfreq*nfreq)
 dim(bic2CL) <- c(nfreq,nfreq)   
 for (i in 1:nfreq)
  {
    for (j in 1:nfreq) bic2CL[i,j] = 100*pchisq(2*segments*bic2[i,j],df=2,lower.tail=T) 
  }
# confidence level for plot contour
 CLset=qchisq(CL/100,df=2)/(2*segments)
 if(verbose) 
  {
    cat("  - 80 percent confidence level=",qchisq(.80,df=2)/(2*segments),"\n")
    cat("  - 90 percent confidence level=",qchisq(.90,df=2)/(2*segments),"\n")
    cat("  - 95 percent confidence level=",qchisq(.95,df=2)/(2*segments),"\n")
    cat("  - 99 percent confidence level=",qchisq(.99,df=2)/(2*segments),"\n")
  }

 if(genplot)
  {            
# spectrum, bispectrum, bicoherence, confidence level             
    dev.new(title=paste("Bispectrum results"),height=6.5,width=5.6)

    if(color==0) colscale=gray((0:499)/499)
    if(color==1) colscale=hcl.colors(12, "YlOrRd", rev = TRUE)
 
    ply=""
    plx=""
    if(logpwr) 
     {
       ply="y"
       plx="x"
     }  

# fig = x1, x2, y1, y2   
    par(fig=c(0.1,0.85,0.1,0.85))
    image(freq[1:nfreq],freq[1:nfreq],Bmag3,col=colscale,useRaster=T,xlab="",ylab="")
    if(is.null(id)) abline(a=0,b=1,lwd=4,col="#FFFFFF70")  
    if(!is.null(id)) 
     {
        for(i in 1:length(id)) abline(a=id[i],b=-1,lwd=1,lty=2)      
     }
    mtext("Frequency",side=2,line=2.5) 
    mtext("Frequency",side=1,line=2.5)  
    mtext("Bispectrum",side=3,line=5,cex=1.75)
    par(fig=c(0.1,0.85,0,0.85))
    image.plot(legend.only=T,col=colscale,zlim=range(Bmag3),horizontal=T)
    par(fig=c(0.1,0.85,0.565,0.975), new=TRUE)
    plot(freq[1:nfreq],pwr[1:nfreq], type='l',axes=FALSE,xaxs="i",xlab="",ylab="",lwd=2,log=ply)           
    par(fig=c(0.625,1,0.1,0.85),new=TRUE)
    plot(pwr[1:nfreq],freq[1:nfreq], type='l',axes=FALSE,yaxs="i",xlab="",ylab="",lwd=2,log=plx)
    par(fig=c(0.1,0.85,0.8,1),new=TRUE)

    dev.new(title=paste("Bicoherence results"),height=6.5,width=5.6)
    par(fig=c(0.1,0.85,0.1,0.85))
    image(freq[1:nfreq],freq[1:nfreq],bic2,col=colscale,useRaster=T,xlab="",ylab="")
    contour(freq[1:nfreq],freq[1:nfreq],bic2,add=T,levels=CLset,drawlabels=F,lwd=2,xlab="",ylab="")   
    if(is.null(id)) abline(a=0,b=1,lwd=4,col="#FFFFFF70")  
    if(!is.null(id)) 
     {
        for(i in 1:length(id)) abline(a=id[i],b=-1,lwd=1,lty=2)      
     }
    mtext("Frequency",side=2,line=2.5) 
    mtext("Frequency",side=1,line=2.5)  
    mtext("Bicoherence",side=3,line=5,cex=1.75)
    par(fig=c(0.1,0.85,0,0.85))
    image.plot(legend.only=T,col=colscale,zlim=range(bic2),horizontal=T)
    par(fig=c(0.1,0.85,0.565,0.975), new=TRUE)
    plot(freq[1:nfreq],pwr[1:nfreq],type='l',axes=FALSE,xaxs="i",xlab="",ylab="",lwd=2,log=ply)                       
    par(fig=c(0.625,1,0.1,0.85),new=TRUE)
    plot(pwr[1:nfreq],freq[1:nfreq],type='l',axes=FALSE,yaxs="i",xlab="",ylab="",lwd=2,log=plx)

    dev.new(title=paste("Bicoherence CL results"),height=6.5,width=5.6)
    par(fig=c(0.1,0.85,0.1,0.85))
    image(freq[1:nfreq],freq[1:nfreq],bic2CL,col=colscale,useRaster=T,xlab="",ylab="")
    contour(freq[1:nfreq],freq[1:nfreq],bic2CL,add=T,levels=c(95),drawlabels=F,lwd=2,xlab="",ylab="")   
    if(is.null(id)) abline(a=0,b=1,lwd=4,col="#FFFFFF70")  
    if(!is.null(id)) 
     {
        for(i in 1:length(id)) abline(a=id[i],b=-1,lwd=1,lty=2)      
     }
    mtext("Frequency",side=2,line=2.5) 
    mtext("Frequency",side=1,line=2.5)  
    mtext("Bicoherence CL",side=3,line=5,cex=1.75)
    par(fig=c(0.1,0.85,0,0.85))
    image.plot(legend.only=T,col=colscale,zlim=range(bic2CL),horizontal=T)
    par(fig=c(0.1,0.85,0.565,0.975), new=TRUE)
    plot(freq[1:nfreq],pwr[1:nfreq], type='l',axes=FALSE,xaxs="i",xlab="",ylab="",lwd=2,log=ply)           
    par(fig=c(0.625,1,0.1,0.85),new=TRUE)
    plot(pwr[1:nfreq],freq[1:nfreq], type='l',axes=FALSE,yaxs="i",xlab="",ylab="",lwd=2,log=plx)
 }

# output
 frequency<-freq[1:nfreq]
 WOSA_Power <- cbind(frequency,pwr[1:nfreq]) 
 colnames(Bmag2) <- frequency 
 WOSA_Bispectrum <- cbind(frequency,Bmag2)
 colnames(bic2) <- frequency 
 WOSA_Bicoherence <- cbind(frequency,bic2)
 if(output==0 || output == 3 || output==4)
  {
    colnames(bic2CL) <- frequency
    WOSA_BicoherenceCL <- cbind(frequency,bic2CL)
  }
     
 if (output == 1) 
  {
    if(verbose) cat("* Returning Magnitude-squared Bispectrum  \n")
    return(WOSA_Bispectrum)   
  }   

 if (output == 2) 
  {
    if(verbose) cat("* Returning Magnitude-squared Bicoherence  \n")
    return(WOSA_Bicoherence)   
  }   

 if (output == 3) 
  {
    if(verbose) cat("* Returning Magnitude-squared Bicoherence Confidence Level  \n")
    return(WOSA_BicoherenceCL)   
  } 

 if (output == 4)
  {
    if(verbose) 
      {
        cat("* Returning the following:  \n")
        cat("1 = WOSA Power \n")
        cat("2 = Magnitude-squared Bispectrum  \n")
        cat("3 = Magnitude-squared Bicoherence  \n")
      }  
    if(verbose) 
     {
       cat("4 = Magnitude-squared Bicoherence Confidence Level  \n")
       cat("5 = Maxima in Magnitude-squared Bicoherence \n")
     }  
    return( list( WOSA_Power, WOSA_Bispectrum, WOSA_Bicoherence, WOSA_BicoherenceCL ))
 }
  
#### END function bicoherence
}