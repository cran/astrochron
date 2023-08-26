### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2019 Stephen R. Meyers
###
###########################################################################
### confAdjust function - (SRM: November 14-30, 2017; December 6, 2017; 
###                             April 13, 2018; April 17, 2018; 
###                             January 8, 2019; June 18, 2022)
###
### apply multiple comparisons correction to spectrum confidence levels
###########################################################################


confAdjust <- function (spec,npts,dt,tbw=3,ntap=5,flow=NULL,fhigh=NULL,output=T,xmin=df,xmax=NULL,pl=1,genplot=T,verbose=T)
{

if(verbose) cat("\n----- PERFORMING multiple comparisons correction to spectral confidence levels -----\n")

spec <- data.frame(spec)
ipts <- length(spec[,1])
df <- spec[2,1]-spec[1,1]

# if there are 3 columns, they must be in the order: freq, power, background
if(length(spec)==3)
 {
   freq <- spec[,1]
   pwrRaw <- spec[,2]
   RawAR <- spec[,3]
 }
 
# if there are more than 3 columns, we assume the results come from lowspec, mtm, mtmPL, 
#  mtmML96, or periodogram. Note that periodgram has a different order!
if(length(spec)>3)
 {
   freq <- spec[,1]
   if(colnames(spec)[2]=="Power" || colnames(spec)[2]== "Prewhite_power") 
    {
      pwrRaw <- spec[,2]
      RawAR <- spec[,5]
    }
# for results from periodogram
   if(colnames(spec)[2]=="Amplitude") 
    {
      pwrRaw <- spec[,3]
      RawAR <- spec[,6]
    }  
 }
 

if(length(spec)<3 || length(spec)>9)
     {
       cat("\n**** ERROR: spec does not have the correct column structure.\n")
       stop("**** TERMINATING NOW!")
     }

if(length(spec)>3 && length(spec)<8)
 {
    cat("\n**** ERROR: spec does not have the correct column structure.\n")
    stop("**** TERMINATING NOW!")
 }

if(is.null(flow)) flow=freq[1]
if(is.null(fhigh)) fhigh=freq[ipts]
if(is.null(xmax)) xmax=freq[ipts]
plotBounds=T
if(min(flow)==freq[1] && max(fhigh)==freq[ipts]) plotBounds=F

# determine degrees of freedom
dof = (2*ntap)
### calculate Rayleigh frequency
Ray <- 1/(dt*npts)

# number of 'independent' frequencies in periodogram (positive frequencies only)
nf1=(freq[ipts]-freq[1])/Ray
ifreqs=0
for (i in 1:length(flow)) ifreqs=append(ifreqs,which( (freq >= flow[i]) & (freq <= fhigh[i]) ))
# remove first entry, which is 0
ifreqs=ifreqs[2:length(ifreqs)]
# remove duplicate values
ifreqs=unique(ifreqs)
# determine number of 'independent' points for multiple correction
nf2=(length(ifreqs)/length(freq))*nf1

if (verbose) 
 {
   cat(" * Number of data points in stratigraphic series:",npts,"\n")
   cat(" * Stratigraphic series length (space or time):",(npts-1)*dt,"\n")
   cat(" * Sampling interval (space or time):",dt,"\n")
   cat(" * Number of frequencies in spectrum:",ipts,"\n")
   cat(" * Frequency spacing:",df,"\n")
   cat(" * Smallest frequency in spectrum:",freq[1],"\n")
   cat(" * Largest frequency in spectrum:",freq[ipts],"\n")
   cat(" * Rayleigh frequency:",Ray,"\n")
   cat(" * MTM Power spectrum halfwidth resolution:",tbw/(npts*dt),"\n")
   cat(" * Number of independent frequencies in full spectrum (for a periodogram):",nf1,"\n")
   cat(" * Number of independent frequencies in region of interest (for a periodogram):",nf2,"\n")
 }

# now apply the correction
AR1_90_local <- RawAR*qchisq(1.0-((1-0.9)/nf2), df=dof)/dof
AR1_95_local <- RawAR*qchisq(1.0-((1-0.95)/nf2), df=dof)/dof
AR1_99_local <- RawAR*qchisq(1.0-((1-0.99)/nf2), df=dof)/dof

AR1_90_global <- RawAR*qchisq(1.0-((1-0.9)/nf1), df=dof)/dof
AR1_95_global <- RawAR*qchisq(1.0-((1-0.95)/nf1), df=dof)/dof
AR1_99_global <- RawAR*qchisq(1.0-((1-0.99)/nf1), df=dof)/dof

AR1_90 <- AR1_90_global
AR1_95 <- AR1_95_global
AR1_99 <- AR1_99_global

AR1_90[ifreqs] <- AR1_90_local[ifreqs]
AR1_95[ifreqs] <- AR1_95_local[ifreqs]
AR1_99[ifreqs] <- AR1_99_local[ifreqs]

### generate plots
if(genplot)
 {
# first plot power spectrum, with red noise model and confidence levels if requested
   dev.new(height=5.5,width=7)
   par(mfrow=c(1,1))
   mtitle=c("Power (black); Noise Model fit (red); 90%CL, 95%CL, 99%CL (dotted)")

# linear frequency, log power
   if(pl==1) plLog="y"
 # log frequency, log power   
   if(pl==2) plLog="xy"
 # linear frequency, linear power    
   if(pl==3) plLog=""
 # log frequency, linear power    
   if(pl==4) plLog="x"
   xlim=c(xmin,xmax)
   ylim=c(min(pwrRaw),max(pwrRaw,AR1_90))
   plot(freq,pwrRaw,xlim=xlim,ylim=ylim,type="l", col="black", xlab="Frequency",ylab="Power",main=mtitle,cex.axis=1.1,cex.lab=1.1,lwd=0,log=plLog)
   if(plotBounds && pl<=2) rect(flow,min(AR1_99,pwrRaw)*10^-4,fhigh,max(AR1_99,pwrRaw),col="#BEBEBE50",border="NA")
   if(plotBounds && pl>2) rect(flow,-1,fhigh,max(AR1_99,pwrRaw),col="#BEBEBE50",border="NA") 
   lines(freq,pwrRaw,col="black",lwd=2)
   lines(freq,RawAR,col="red",lwd=2)
   lines(freq,AR1_90,col="red",lwd=1,lty=3)
   lines(freq,AR1_95,col="red",lwd=1,lty=3)
   lines(freq,AR1_99,col="red",lwd=1,lty=3)
# end genplot section
 }
  
if (output) 
 {
       spectrum <- data.frame(cbind(freq,pwrRaw,RawAR,AR1_90,AR1_95,AR1_99))
       colnames(spectrum)[1] <- 'Frequency'
       colnames(spectrum)[2] <- 'Power'
       colnames(spectrum)[3] <- 'Background_fit'
       colnames(spectrum)[4] <- 'CL90_power'
       colnames(spectrum)[5] <- 'CL95_power'
       colnames(spectrum)[6] <- 'CL99_power'
       return(spectrum)
 }

#### END function confAdjust
}
