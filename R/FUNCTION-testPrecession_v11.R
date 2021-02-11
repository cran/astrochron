### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2021 Stephen R. Meyers
###
###########################################################################
### function testPrecession - (SRM: July 10, 2013; August 9, 2013; August 14, 2013;
###                         Sept. 16-17, 2014; Jan. 8, 2015; Jan. 20-22, 2015;
###                         Jan. 29, 2015; February 4, 2015; February 23, 2015;
###                         March 11, 2015; September 10, 2015; July 22, 2016;
###                         October 26, 2016; June 7-8, 2017; June 20, 2018;
###                         June 24, 2018; January 14, 2021)
###
### Perform astrochonologic testing as in Zeeden et al. (2015).
###########################################################################

testPrecession <- function(dat,nsim=1000,gen=1,edge=0.025,maxNoise=1,rho=NULL,detrendEnv=T,solution=NULL,output=F,genplot=T,verbose=T)
{

   if(verbose) { cat("\n----- PERFORMING ZEEDEN ET AL. (2015) ASTROCHRONOLOGIC TESTING -----\n") }
   dat <- data.frame(dat)

# sort data to ensure increasing order
   if(verbose) cat(" * Sorting data into ensure increasing order, removing empty entries\n")
   dat <- dat[order(dat[,1],na.last=NA,decreasing=F),]
   npts <- length(dat[,1]) 

# standardize
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

   if(dt > 5)
     {
       cat("\n**** ERROR: temporal sampling interval must be <= 5 ka.\n")
       stop("**** TERMINATING NOW!")
     }

   if(verbose) { cat(" * Number of data points=", npts,"\n") }
   if(verbose) { cat(" * Sample interval=", dt,"\n") }
   
   if(!is.null(solution)) 
    {
# check to ensure that solution covers the analyzed series    
     min2=min(solution[,1])
     max2=max(solution[,1])
     if(min2>mindat || max2<maxdat)
      {
       cat("\n**** ERROR: solution does not span the entire temporal interval of your input series.\n")
       stop("**** TERMINATING NOW!")
      }     
    }
    
###########################################################################
### get theoretical solutions if needed
###########################################################################
   if(is.null(solution)) solution=getLaskar(sol="la04",verbose=F)

###########################################################################
### extract sinw and esinw
###########################################################################
   xout <- seq(dat[1,1],dat[npts,1],by=dt)
# save axial precession
   axial <- data.frame(approx(solution[,1],solution[,2],xout,method="linear",n=npts))
   axial[2] = sin(axial[2])
   eccIn <- data.frame(approx(solution[,1],solution[,4],xout,method="linear",n=npts))
# determine esinw
   precT = axial
   precT[2] = eccIn[2]*precT[2]
# standardize axial precession and esinw
   axial[2]= axial[2]-mean(axial[,2])
   axial[2]= axial[2]/sd(axial[,2])     
   precT[2]= precT[2]-mean(precT[,2])
   precT[2]= precT[2]/sd(precT[,2])

###########################################################################
### determine eccentricity signal for time interval
###########################################################################
# bandpass filter to be consistent with tuned data processing
   precT_BP=taner(precT,padfac=2,flow=0.029,fhigh=0.12,roll=10^3,demean=T,detrend=F,genplot=F,verbose=F)
# determine instantaneous amplitude
   precT_BP_Hil=hilbert(precT_BP,padfac=2,demean=T,detrend=F,genplot=F,verbose=F)
# lowpass instantaneous amplitude
   ecc=taner(precT_BP_Hil,padfac=2,fhigh=0.013,flow=-0.013,roll=10^4,demean=T,detrend=F,genplot=F,verbose=F)
   if(detrendEnv) ecc=detrend(ecc,genplot=F,verbose=F)

###########################################################################
### determine region of edge effects, to remove
###########################################################################
  if(edge > 0) ii=floor(edge*length(ecc[,1])):ceiling((1-edge)*length(ecc[,1]))
  if(edge == 0) ii=1:length(ecc[,1])

###########################################################################
### process axial prececession (sinw) to determine residual correlation
### with full eccentricity signal modulation
###########################################################################
# bandpass filter to be consistent with tuned data processing
   axial_BP=taner(axial,padfac=2,flow=0.029,fhigh=0.12,roll=10^3,demean=T,addmean=T,detrend=F,genplot=F,verbose=F)
# determine instantaneous amplitude
   axial_BP_Hil=hilbert(axial_BP,padfac=2,demean=T,addmean=T,detrend=F,genplot=F,verbose=F)
# lowpass instantaneous amplitude
   axial_amp=taner(axial_BP_Hil,padfac=2,fhigh=0.013,flow=-0.013,roll=10^4,demean=T,addmean=T,detrend=F,genplot=F,verbose=F)   
   if(detrendEnv) axial_amp=detrend(axial_amp,genplot=F,verbose=F)
# calculate spearman correlation. remove edge effects using 'ii' indices
   axialCor=cor(axial_amp[ii,2],ecc[ii,2],method=c("spearman"))
   if(verbose) 
    {
      cat("\nINITIAL SCREENING USING THEORY: \n")
      cat(" * Theoretical esinw AM vs. theoretical sinw AM, r=",axialCor,"\n")
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
   if(axialCor>0 || gen == 2)
    {
      if(is.null(rho)) 
       {
        lag0 <- dat[1:(npts-1),2]
        lag1 <- dat[2:npts,2]
        rho = cor(lag0,lag1)
       }
    }

# we will take a coservative approach, and add noise until a correlation of <= 0 is achieved 
   if(axialCor>0)
    {
      if(verbose) 
        {
          cat(" * PLEASE WAIT: Estimating noise level required to suppress\n")
          cat("      residual imposed eccentricity modulations. \n")
        }
      for (i in 1:1000)
       {
         axialNow=axial
         ar1Now=ar1(npts=npts,dt=dt,rho=rho,mean=0,sdev=maxNoise*i/1000,nsim=1,genplot=F,verbose=F)
         axialNow[2]=axialNow[2]+ar1Now[2]
# bandpass filter to be consistent with tuned data processing
         axialNow_BP=taner(axialNow,padfac=2,flow=0.029,fhigh=0.12,roll=10^3,demean=T,addmean=T,detrend=F,genplot=F,verbose=F)
# determine instantaneous amplitude
         axialNow_BP_Hil=hilbert(axialNow_BP,padfac=2,demean=T,addmean=T,detrend=F,genplot=F,verbose=F)
# lowpass instantaneous amplitude
         axialNow_amp=taner(axialNow_BP_Hil,padfac=2,fhigh=0.013,flow=-0.013,roll=10^4,demean=T,addmean=T,detrend=F,genplot=F,verbose=F)   
         if(detrendEnv) axialNow_amp=detrend(axialNow_amp,genplot=F,verbose=F)
# calculate spearman correlation. remove edge effects using 'ii' indices
         axialCorNoise=cor(axialNow_amp[ii,2],ecc[ii,2],method=c("spearman"))
         if(axialCorNoise <= 0) break
         if(axialCorNoise >0 && i == 1000) 
           {
             cat("\n**** ERROR: Noise simulations were not able to remove imposed eccentricity modulations.\n")
             cat("****        Try increasing parameter maxNoise.\n")
             stop("**** TERMINATING NOW!")
           }              
      }
      if(verbose) 
       {
         cat(" * Strength of noise added (standardized) to sinw =",maxNoise*i/1000,", iteration =",i,"\n")
         cat(" * Theoretical esinw AM vs. 'sinw + noise' AM, r=",axialCorNoise,"\n")
       }  
# add noise to standardized data
      dat2[2]=dat[2]+ar1Now[2]
      inoise=T
    } 
    
###########################################################################
### process tuned data series WITHOUT noise addition
###########################################################################    
# bandpass
   prec_BP_noNoise=taner(dat,padfac=2,flow=0.029,fhigh=0.12,roll=10^3,demean=T,detrend=F,genplot=F,verbose=F)
# determine instantaneous amplitude
   prec_BP_Hil_noNoise=hilbert(prec_BP_noNoise,padfac=2,demean=T,detrend=F,genplot=F,verbose=F)
# lowpass instantaneous amplitude
   prec_BP_Hil_Lowpass_noNoise=taner(prec_BP_Hil_noNoise,padfac=2,fhigh=0.013,flow=-0.013,roll=10^4,demean=T,detrend=F,genplot=F,verbose=F)
# save a copy for plotting later...
   plotSeries=prec_BP_Hil_Lowpass_noNoise
   if(detrendEnv) prec_BP_Hil_Lowpass_noNoise=detrend(prec_BP_Hil_Lowpass_noNoise,genplot=F,verbose=F)
# calculate spearman correlation. remove edge effects using 'ii' indices 
   datCor_noNoise=cor(prec_BP_Hil_Lowpass_noNoise[ii,2],ecc[ii,2],method=c("spearman"))
   if(verbose) 
    {
      cat("\nEVALUATION OF YOUR DATA: \n")
      cat(" * Theoretical esinw AM vs. observed precession AM, r=",datCor_noNoise,"\n")
    }

###########################################################################
### process tuned stratigraphic series WITH noise addition
###########################################################################    
# bandpass
   prec_BP=taner(dat2,padfac=2,flow=0.029,fhigh=0.12,roll=10^3,demean=T,addmean=T,detrend=F,genplot=F,verbose=F)
# determine instantaneous amplitude
   prec_BP_Hil=hilbert(prec_BP,padfac=2,demean=T,addmean=T,detrend=F,genplot=F,verbose=F)
# lowpass instantaneous amplitude
   prec_BP_Hil_Lowpass=taner(prec_BP_Hil,padfac=2,fhigh=0.013,flow=-0.013,roll=10^4,demean=T,addmean=T,detrend=F,genplot=F,verbose=F)
   if(detrendEnv) prec_BP_Hil_Lowpass=detrend(prec_BP_Hil_Lowpass,genplot=F,verbose=F)
# calculate spearman correlation. remove edge effects using 'ii' indices 
   datCor=cor(prec_BP_Hil_Lowpass[ii,2],ecc[ii,2],method=c("spearman"))
   if(verbose && inoise) cat(" * Theoretical esinw AM vs. 'observed precession + noise' AM, r=",datCor,"\n\n")

###########################################################################
### summary plot
###########################################################################    
# generate plot comparing data results to target
   if(genplot)
    {
     dev.new(title = "testPrecession Amplitude Modulation Assessment", height=7,width=8)
     if(inoise) par(mfrow=c(3,1))
     if(!inoise) par(mfrow=c(2,1))
     ymaxPlot=max(prec_BP_noNoise[,2],prec_BP_Hil_noNoise[,2],plotSeries[,2])
     yminPlot=min(prec_BP_noNoise[,2],prec_BP_Hil_noNoise[,2],plotSeries[,2])
# plot 1: data without noise addition
     plot(prec_BP_noNoise,type="l",xlab="Time (kiloyears)",ylab="Climate Proxy Value",main="Filtered record (black), instantaneous AM (red) and lowpassed AM (blue)",ylim=c(yminPlot,ymaxPlot))
     lines(prec_BP_Hil_noNoise,col="red")
     lines(plotSeries,col="blue")
# plot 2: data lowpassed AM vs theory
     s_prec_BP_Hil_Lowpass_noNoise=s(prec_BP_Hil_Lowpass_noNoise,verbose=F)
     s_ecc=s(ecc,verbose=F)
     ymaxPlot=max(s_prec_BP_Hil_Lowpass_noNoise[,2],s_ecc[,2])
     yminPlot=min(s_prec_BP_Hil_Lowpass_noNoise[,2],s_ecc[,2])
     plot(s_prec_BP_Hil_Lowpass_noNoise,type="l",col="blue",xlab="Time (kiloyears)",ylab="Standardized Variables",main="Comparison of theoretical (black) and observed (blue) modulation",ylim=c(yminPlot,ymaxPlot),lwd=2)
     lines(s_ecc,lwd=2)
     if(edge > 0) 
      { 
        abline(v=c(dat[ii[1],1]),lty=2,lwd=2,col="gray")
        abline(v=c(dat[ii[length(ii)],1]),lty=2,lwd=2,col="gray")
      }  
# now with noise addition
     if(inoise)
      { 
        s_prec_BP_Hil_Lowpass=s(prec_BP_Hil_Lowpass,verbose=F)
        plot(s_prec_BP_Hil_Lowpass,type="l",col="blue",xlab="Time (kiloyears)",ylab="Standardized Variables",main="Comparison of theoretical (black) and observed + noise (blue) modulation",ylim=c(yminPlot,ymaxPlot),lwd=2)
        lines(s_ecc,lwd=2)
        if(edge > 0) 
         { 
           abline(v=c(dat[ii[1],1]),lty=2,lwd=2,col="gray")
           abline(v=c(dat[ii[length(ii)],1]),lty=2,lwd=2,col="gray")
         }
       }    
    }

    if(nsim==0)
     { 
       out = data.frame(cbind(prec_BP_noNoise[1],prec_BP_noNoise[2],prec_BP_Hil_noNoise[2],prec_BP_Hil_Lowpass_noNoise[2],ecc[2],prec_BP[2],prec_BP_Hil[2],prec_BP_Hil_Lowpass[2]))
       colnames(out)[1]<-"Time"
       colnames(out)[2]<-"Precession_BP"
       colnames(out)[3]<-"Precession_BP_Amp"
       colnames(out)[4]<-"Precession_BP_Amp_LP"
       colnames(out)[5]<-"Eccentricity"
       colnames(out)[6]<-"Precession+Noise_BP"
       colnames(out)[7]<-"Precession+Noise_BP_Amp"
       colnames(out)[8]<-"Precession+Noise_BP_Amp_LP"
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
# note that we are generating surrogates from stratigraphic series + noise (if added)
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
      surdat=data.frame(cbind(dat[,1],sur[,i]))
      sur_BP=taner(surdat,padfac=2,flow=0.029,fhigh=0.12,roll=10^3,demean=T,detrend=F,genplot=F,verbose=F)
# determine instantaneous amplitude
      sur_BP_Hil=hilbert(sur_BP,padfac=2,demean=T,detrend=F,genplot=F,verbose=F)
# lowpass instantaneous amplitude
      sur_BP_Hil_Lowpass=taner(sur_BP_Hil,padfac=2,fhigh=0.013,flow=-0.013,xmax=.1,roll=10^4,demean=T,detrend=F,genplot=F,verbose=F)
      if(detrendEnv) sur_BP_Hil_Lowpass=detrend(sur_BP_Hil_Lowpass,genplot=F,verbose=F)
# calculate spearman rank correlation coefficient 
      simcor[i]=cor(sur_BP_Hil_Lowpass[,2],ecc[,2],method=c("spearman"))
     }

     if(genplot)
      {
        dev.new(title = "testPrecession Monte Carlo Results", height = 5, width = 6)
        par(mfrow=c(1,1))
        plot(density(simcor, bw="nrd0"),xlim=c(-1,1),type="l",col="black",xlab="Correlation Coefficient",main="testPrecession Monte Carlo Results",cex.lab=1.1,lwd=2)
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
     
# end function testPrecession     
}     