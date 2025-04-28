### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2025 Stephen R. Meyers
###
###########################################################################
### function testTilt - (SRM: December 27, 2015; January 5-6, 2016; 
###                           June 1-8, 2017; June 20, 2018; June 22, 2018;
###                           June 24, 2018; January 14, 2021; January 12, 2025)
###
### Perform astrochonologic testing using obliquity modulations
### as in Zeeden et al. (2019).
###########################################################################


testTilt <- function(dat,nsim=1000,gen=1,edge=0.025,cutoff=1/150,maxNoise=0.25,rho=NULL,detrendEnv=T,solution=NULL,output=F,genplot=T,verbose=T)
{
   if(verbose) { cat("\n----- PERFORMING ASTROCHRONOLOGIC TESTING OF OBLIQUITY AM-----\n") }
   dat <- data.frame(dat)

# sort data to ensure increasing order
   if(verbose) cat(" * Sorting data into ensure increasing order, removing empty entries\n")
   dat <- dat[order(dat[,1],na.last=NA,decreasing=F),]
   npts <- length(dat[,1]) 

# standardize dat
   dat[2] <- dat[2] - mean(dat[,2])
   dat[2] <- dat[2]/sd(dat[,2])
         
# error checking 
   maxdat=max(dat[,1])
   if(maxdat > 50000) 
      {
       cat("\n**** ERROR: this approach is only applicable when testing astrochronologies <= 50 Ma!\n")
       stop("**** TERMINATING NOW!")
     }

# error checking 
   mindat=min(dat[,1])
   if(mindat < 0) 
      {
       cat("\n**** ERROR: time (first column of input series) must be positive.\n")
       stop("**** TERMINATING NOW!")
     }

   dtest <- dat[2:npts,1]-dat[1:(npts-1),1] 
   epsm=1e-9
   if( (max(dtest)-min(dtest)) > epsm ) 
     {
       cat("\n**** ERROR: temporal sampling interval is not uniform.\n")
       stop("**** TERMINATING NOW!")
     }

   dt = dat[2,1]-dat[1,1]

   if(dt > 10)
     {
       cat("\n**** ERROR: temporal sampling interval must be <= 10 ka.\n")
       stop("**** TERMINATING NOW!")
     }

   if(verbose) { cat(" * Number of data points=", npts,"\n") }
   if(verbose) { cat(" * Sample interval=", dt,"\n") }
   
   if(!is.null(solution)) 
    {
# check to ensure that obliquity covers the analyzed series    
     min2=min(solution[,1])
     max2=max(solution[,1])
     if(min2>mindat || max2<maxdat)
      {
       cat("\n**** ERROR: obl does not span the entire temporal interval of your input series.\n")
       stop("**** TERMINATING NOW!")
      }     
    }
      
###########################################################################
### determine theoretical obliquity target signal for time interval
###########################################################################
   if(is.null(solution)) solution=getLaskar(sol="la04",verbose=F)
   xout <- seq(dat[1,1],dat[npts,1],by=dt)
# save axial obliquity 
   oblT <- data.frame(approx(solution[,1],solution[,3],xout,method="linear",n=npts))
# create demodulated obliquity
   demod <- oblT
   demod[2]= demod[2]-mean(demod[,2])
   demod[2]= demod[2]/sd(demod[,2]) 
   demod[2]=demod[2]/hilbert(demod,genplot=F,check=F,verbose=F)[2]

###########################################################################
### process obliquity target
###########################################################################
# bandpass filter to be consistent with tuned data processing
   oblT_BP=taner(oblT,padfac=2,flow=0.015,fhigh=0.038,roll=10^20,demean=T,addmean=T,detrend=F,genplot=F,check=F,verbose=F)
# determine instantaneous amplitude
   oblT_BP_Hil=hilbert(oblT_BP,padfac=2,demean=T,addmean=T,detrend=F,genplot=F,check=F,verbose=F)
# lowpass instantaneous amplitude
   oblT_amp=taner(oblT_BP_Hil,padfac=2,fhigh=cutoff,roll=10^20,demean=T,addmean=T,detrend=F,genplot=F,check=F,verbose=F) 
   if(detrendEnv) oblT_amp=detrend(oblT_amp,genplot=F,verbose=F)

###########################################################################
### determine region of edge effects, to remove
###########################################################################
  if(edge > 0) ii=floor(edge*length(oblT_amp[,1])):ceiling((1-edge)*length(oblT_amp[,1]))
  if(edge == 0) ii=1:length(oblT_amp[,1])

###########################################################################
### process demodulated obliquity to determine residual correlation with
### full obliquity signal modulation
###########################################################################
# bandpass filter to be consistent with tuned data processing
   demod_BP=taner(demod,padfac=2,flow=0.015,fhigh=0.038,roll=10^20,demean=T,addmean=T,detrend=F,genplot=F,check=F,verbose=F)
# determine instantaneous amplitude
   demod_BP_Hil=hilbert(demod_BP,padfac=2,demean=T,addmean=T,detrend=F,genplot=F,check=F,verbose=F)
# lowpass instantaneous amplitude
   demod_amp=taner(demod_BP_Hil,padfac=2,fhigh=cutoff,roll=10^20,demean=T,addmean=T,detrend=F,genplot=F,check=F,verbose=F)   
   if(detrendEnv) demod_amp=detrend(demod_amp,genplot=F,verbose=F)
# calculate spearman correlation. remove edge effects using 'ii' indices
   demodCor=cor(demod_amp[ii,2],oblT_amp[ii,2],method=c("spearman"))
   if(verbose) 
    {
      cat("\nINITIAL SCREENING USING THEORY: \n")
      cat(" * Theoretical obliquity AM vs. demodulated theoretical\n")
      cat("      obliquity AM, r=",demodCor,"\n")
    }
 
###########################################################################
### determine appropriate ar1 noise level to add, if needed to supress 
### residual imposed modulation, then add it to original stratigraphic
### series
###########################################################################
# switch to indentify if noise addition simulations were performed
   inoise=F
   dat2=dat

 # calculate rho for noise addition simulations and/or Monte Carlo if selected
   if(demodCor>0 || gen == 2)
    {
# npts is the length of data vector 'dat', mean of dat is already zero
      if(is.null(rho)) 
       { 
         rho=sum(dat[1:(npts-1),2] * dat[2:npts,2]) / sum(dat[,2]^2)  
         if(verbose) cat("\n * Estimated AR1 coefficient=",rho,"\n")
       }
    }

# we will take a coservative approach, and add noise until a correlation of <= 0 is achieved 
   if(demodCor>0)
    {
      if (verbose) 
       {
         cat(" * PLEASE WAIT: Estimating noise level required to suppress\n")
         cat("      residual imposed obliquity modulations. \n")
       }  
# standardize the obliquity signal
      demodStd=demod
      demodStd[2]= demodStd[2]-mean(demodStd[,2])
      demodStd[2]= demodStd[2]/sd(demodStd[,2]) 
      for (i in 1:1000)
       {
         demodNow=demodStd
         ar1Now=ar1(npts=npts,dt=dt,rho=rho,mean=0,sdev=maxNoise*i/1000,nsim=1,genplot=F,verbose=F)
         demodNow[2]=demodNow[2]+ar1Now[2]
# bandpass filter to be consistent with tuned data processing
         demodNow_BP=taner(demodNow,padfac=2,flow=0.015,fhigh=0.038,roll=10^20,demean=T,addmean=T,detrend=F,genplot=F,check=F,verbose=F)
# determine instantaneous amplitude
         demodNow_BP_Hil=hilbert(demodNow_BP,padfac=2,demean=T,addmean=T,detrend=F,genplot=F,check=F,verbose=F)
# lowpass instantaneous amplitude
         demodNow_amp=taner(demodNow_BP_Hil,padfac=2,fhigh=cutoff,roll=10^20,demean=T,addmean=T,detrend=F,genplot=F,check=F,verbose=F)   
         if(detrendEnv) demodNow_amp=detrend(demodNow_amp,genplot=F,verbose=F)
# calculate spearman correlation. remove edge effects using 'ii' indices
         demodCorNoise=cor(demodNow_amp[ii,2],oblT_amp[ii,2],method=c("spearman"))
         if(demodCorNoise <= 0) break
         if(demodCorNoise >0 && i == 1000) 
           {
             cat("\n**** ERROR: Noise simulations were not able to removed imposed obliquity modulations.\n")
             cat("****        Try increasing parameter maxNoise.\n")
             stop("**** TERMINATING NOW!")
           }              
      }
      if(verbose) 
       {
         cat(" * Strength of noise added (standardized) to demodulated theoretical\n")
         cat("      obliquity =",maxNoise*i/1000,", iteration =",i,"\n")
         cat(" * Theoretical obliquity AM vs. 'demodulated theoretical obliquity + \n")
         cat("      noise' AM, r=",demodCorNoise,"\n")
       }  
# add noise to standardized data
      dat2[2]=dat[2]+ar1Now[2]
      inoise=T
    } 
    
###########################################################################
### process tuned stratigraphic series WITHOUT noise addition
###########################################################################    
# bandpass
   obl_BP_noNoise=taner(dat,padfac=2,flow=0.015,fhigh=0.038,roll=10^20,demean=T,addmean=T,detrend=F,genplot=F,check=F,verbose=F)
# determine instantaneous amplitude
   obl_BP_Hil_noNoise=hilbert(obl_BP_noNoise,padfac=2,demean=T,addmean=T,detrend=F,genplot=F,check=F,verbose=F)
# lowpass instantaneous amplitude
   obl_BP_Hil_Lowpass_noNoise=taner(obl_BP_Hil_noNoise,padfac=2,fhigh=cutoff,roll=10^20,demean=T,addmean=T,detrend=F,genplot=F,check=F,verbose=F)
# save a copy for plotting later...
   plotSeries=obl_BP_Hil_Lowpass_noNoise
   if(detrendEnv) obl_BP_Hil_Lowpass_noNoise=detrend(obl_BP_Hil_Lowpass_noNoise,genplot=F,verbose=F)
# calculate spearman correlation. remove edge effects using 'ii' indices 
   datCor_noNoise=cor(obl_BP_Hil_Lowpass_noNoise[ii,2],oblT_amp[ii,2],method=c("spearman"))
   if(verbose) 
    {
      cat("\nEVALUATION OF YOUR DATA: \n")
      cat(" * Theoretical obliquity AM vs. observed obliquity AM, r=",datCor_noNoise,"\n")
    }

###########################################################################
### process tuned stratigraphic series WITH noise addition
###########################################################################    
# bandpass
   obl_BP=taner(dat2,padfac=2,flow=0.015,fhigh=0.038,roll=10^20,demean=T,addmean=T,detrend=F,genplot=F,check=F,verbose=F)
# determine instantaneous amplitude
   obl_BP_Hil=hilbert(obl_BP,padfac=2,demean=T,addmean=T,detrend=F,genplot=F,check=F,verbose=F)
# lowpass instantaneous amplitude
   obl_BP_Hil_Lowpass=taner(obl_BP_Hil,padfac=2,fhigh=cutoff,roll=10^20,demean=T,addmean=T,detrend=F,genplot=F,check=F,verbose=F)
   if(detrendEnv) obl_BP_Hil_Lowpass=detrend(obl_BP_Hil_Lowpass,genplot=F,verbose=F)
# calculate spearman correlation. remove edge effects using 'ii' indices 
   datCor=cor(obl_BP_Hil_Lowpass[ii,2],oblT_amp[ii,2],method=c("spearman"))
   if(verbose && inoise) cat(" * Theoretical obliquity AM vs. 'observed obliquity + noise' AM, r=",datCor,"\n\n")

###########################################################################
### summary plot
###########################################################################    
# generate plot comparing data results to target
   if(genplot)
    {
     dev.new(title = "testObliquity Amplitude Modulation Assessment", height=7,width=8)
     par(mfrow=c(3,1))
     ymaxPlot=max(obl_BP_noNoise[,2],obl_BP_Hil_noNoise[,2],plotSeries[,2])
     yminPlot=min(obl_BP_noNoise[,2],obl_BP_Hil_noNoise[,2],plotSeries[,2])
# plot 1: data without noise addition
     plot(obl_BP_noNoise,type="l",xlab="Time (kiloyears)",ylab="Climate Proxy Value",main="Filtered record (black), instantaneous AM (red) and lowpassed AM (blue)",ylim=c(yminPlot,ymaxPlot))
     lines(obl_BP_Hil_noNoise,col="red")
     lines(plotSeries,col="blue")
# plot 2: data lowpassed AM vs theory
     s_obl_BP_Hil_Lowpass_noNoise=s(obl_BP_Hil_Lowpass_noNoise,verbose=F)
     s_oblT_amp=s(oblT_amp,verbose=F)
     ymaxPlot=max(s_obl_BP_Hil_Lowpass_noNoise[,2],s_oblT_amp[,2])
     yminPlot=min(s_obl_BP_Hil_Lowpass_noNoise[,2],s_oblT_amp[,2])
     plot(s_obl_BP_Hil_Lowpass_noNoise,type="l",col="blue",xlab="Time (kiloyears)",ylab="Standardized Variables",main="Comparison of theoretical (black) and observed (blue) modulation",ylim=c(yminPlot,ymaxPlot),lwd=2)
     lines(s_oblT_amp,lwd=2)
     if(edge > 0) 
      { 
        abline(v=c(dat[ii[1],1]),lty=2,lwd=2,col="gray")
        abline(v=c(dat[ii[length(ii)],1]),lty=2,lwd=2,col="gray")
      }  
# now with noise addition
     if(inoise)
      { 
        s_obl_BP_Hil_Lowpass=s(obl_BP_Hil_Lowpass,verbose=F)
        plot(s_obl_BP_Hil_Lowpass,type="l",col="blue",xlab="Time (kiloyears)",ylab="Standardized Variables",main="Comparison of theoretical (black) and observed + noise (blue) modulation",ylim=c(yminPlot,ymaxPlot),lwd=2)
        lines(s_oblT_amp,lwd=2)
        if(edge > 0) 
         { 
           abline(v=c(dat[ii[1],1]),lty=2,lwd=2,col="gray")
           abline(v=c(dat[ii[length(ii)],1]),lty=2,lwd=2,col="gray")
         }  
       }  
    }

    if(nsim==0)
      { 
        out = data.frame(cbind(obl_BP[1],obl_BP[2],obl_BP_Hil[2],obl_BP_Hil_Lowpass[2],oblT_amp[2],obl_BP_noNoise[2],obl_BP_Hil_noNoise[2],obl_BP_Hil_Lowpass_noNoise[2]))
        colnames(out)[1]<-"Time"
        colnames(out)[2]<-"Obliquity_BP"
        colnames(out)[3]<-"Obliquity_BP_AM"
        colnames(out)[4]<-"Obliquity_BP_AM_LP"
        colnames(out)[5]<-"Obliquity_Target_AM"
        colnames(out)[6]<-"Obliquity+Noise_BP"
        colnames(out)[7]<-"Obliquity+Noise_BP_AM"
        colnames(out)[8]<-"Obliquity+Noise_BP_AM_LP"
     }
     
# check for correlation <= 0   
   if(datCor<=0)
    {
       cat("\n**** ERROR: correlation must be positive to proceed.\n")
       if(nsim==0 && output) return(out) 
       stop("**** TERMINATING NOW!")
    }
       
###########################################################################
### perform Monte Carlo simulations using phase-randomized or AR1 surrogates
###########################################################################   
if(nsim>0)
{   
# note that we are generating surrogates from stratigraphic series + noise
   if(gen==1)sur=surrogates(dat2,nsim=nsim,preserveMean=T,std=T,genplot=F,verbose=T)
   if(gen==2)sur=ar1(npts=npts,dt=dt,mean=0,sdev=1,rho=rho,shuffle=F,nsim=nsim,genplot=F,verbose=T)

   if(verbose) 
    {  
      cat("\n * PLEASE WAIT: Performing", nsim,"simulations\n")
      cat("\n0%       25%       50%       75%       100%\n")
# create a progress bar
      progress = utils::txtProgressBar(min = 0, max = nsim, style = 1, width=43)
    }
    
   simcor=double(nsim)
   for (i in 1:nsim)
    {
      if(verbose) utils::setTxtProgressBar(progress, i)
      surdat=data.frame(cbind(dat2[,1],sur[,i]))
      sur_BP=taner(surdat,padfac=2,flow=0.015,fhigh=0.038,roll=10^20,demean=T,addmean=T,detrend=F,genplot=F,check=F,verbose=F)
# determine instantaneous amplitude
      sur_BP_Hil=hilbert(sur_BP,padfac=2,demean=T,addmean=T,detrend=F,genplot=F,check=F,verbose=F)
# lowpass instantaneous amplitude
      sur_BP_Hil_Lowpass=taner(sur_BP_Hil,padfac=2,fhigh=cutoff,roll=10^20,demean=T,addmean=T,detrend=F,genplot=F,check=F,verbose=F)
      if(detrendEnv) sur_BP_Hil_Lowpass=detrend(sur_BP_Hil_Lowpass,genplot=F,verbose=F)
# calculate spearman rank correlation coefficient. remove edge effects using 'ii' indices
      simcor[i]=cor(sur_BP_Hil_Lowpass[ii,2],oblT_amp[ii,2],method=c("spearman"))
     }
       
    if(genplot)
      {
        dev.new(title = "testTilt Monte Carlo Results", height = 5, width = 6)
        par(mfrow=c(1,1))
        plot(density(simcor, bw="nrd0"),xlim=c(-1,1),type="l",col="black",xlab="Correlation Coefficient",main="testTilt Monte Carlo Results",cex.lab=1.1,lwd=2)
        polygon(density(simcor),col="red",border=NA)
        abline(v=datCor,col="blue",lwd=2,lty=3)
        mtext(round(datCor,digits=5),side=3,line=0,at=datCor,cex=1,font=4,col="blue")
    }  

# now sort simulation correlation results, determine how many have r > your result
# this follows Ebisuzaki (1997)
    pnumgt = sum(abs(simcor)>abs(datCor))
    cat("\n * Number of simulations with |r<sim>| > |r<dat>| = ", pnumgt,"\n")
    ppvalue=pnumgt/nsim
    if(ppvalue < (10/nsim) && (10/nsim) <=1 ) 
     {
       cat("\n * P-value <", 10/nsim,"\n")
       if(genplot) mtext(paste("p <",10/nsim),side=3,line=0,at=-0.75,cex=1,font=4,col="red")
     }  
    if(ppvalue >= (10/nsim) && (10/nsim) <=1 ) 
     {
       cat("\n * P-value =", ppvalue,"\n")    
       if(genplot) mtext(paste("p =",round(ppvalue,digits=5)),side=3,line=0,at=-0.75,cex=1,font=4,col="red")
     }  
    if((10/nsim) > 1 ) 
     {
       cat("\n * P-value = 1 \n") 
       if(genplot) mtext("p = 1",side=3,line=0,at=-0.75,cex=1,font=4,col="red") 
     }  

    if(verbose) close(progress) 
# return simlulation results for plotting
    if(output) return(simcor)

# end nsim>0 section
} 
 
if(nsim==0 && output) return(out)
     
# end function testTilt 
}     