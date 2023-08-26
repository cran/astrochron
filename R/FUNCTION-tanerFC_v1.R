### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2023 Stephen R. Meyers
###
###########################################################################
### function tanerFC - taner bandpass filter. this version allows one to 
###                     pass Fourier coefficients into the function 
###                     (SRM: July 30, 2022)
###
###########################################################################
# modified from FUNCTION-taner_v9.R

tanerFC <- function (fc,npts,flow=NULL,fhigh=NULL,roll=10^3,output=1,genplot=T,verbose=T)
{

   if(verbose)   { cat("\n----- APPLY TANER BANDPASS FILTER TO FOURIER COEFFICIENTS-----\n") } 
 
# fc contains Fourier coefficients calculated with 'periodogram' (using padfac >=1 and nrm=0).
# when padfac >=1, periodogram guarantees an even number of points so Nyquist exists, as
# required below.
   fc <- data.frame(fc)
   freq <- fc[,1]
   ft <- complex(real=fc[,2],imaginary=fc[,3])
   nf <- length(fc[,1])      
   nyqfreq = 0.5*nf + 1
   negfreq = 0.5*nf + 2
   dt=1/(2*freq[nyqfreq])
   if(is.null(npts) && verbose) 
     {
       cat("\n**** ERROR: npts (number of points prior to padding) must be specified.\n")
       stop("**** TERMINATING NOW!")
     }  

# when flow is not specified, automatically perform lowpass filtering
   if(is.null(flow)) 
    {
      if(verbose)  { cat("\n **** NOTE: flow not specified. Will perform lowpass filtering.\n")}    
      flow = -1*fhigh
    }   
       
   if(flow > fhigh)  {fh <- fhigh; fhigh <- flow; flow <- fh}
   
### assign parameters for taner filter, then generate filter. this portion modified from 
### L. Hinnov FORTRAN code TANER.FOR
### also see: http://www.rocksolidimages.com/pdf/attrib_revisited.htm#_Toc328470897
### center of filter
    fcent=(flow+fhigh)/2
### convert to angular frequency
    wl = 2*pi*flow
    wc = 2*pi*fcent
    wh = 2*pi*fhigh
    bw = wh - wl
### amp2 is the maximum value of the taper
### set to be similar to Butterworth filter
    amp2 = 1/sqrt(2)
### roll is the roll-off octave
    arg1 = 1 - (roll*log(10)) / (20*log(amp2))
    arg1 = log(arg1)
    arg2 = ((bw+2)/bw)^2
    arg2 = log(arg2)
    twod = 2*arg1/arg2
### generate filter for positive frequencies
    dw = 2*pi/(nf*dt)
    w = ((1:nyqfreq)-1)*dw
    arg = (2*abs(w-wc)/bw)^twod
    darg=-1*arg
    filter_pos =  (amp2 * exp(darg))
### generate filter for negative frequencies
    w = ((negfreq:nf)-nf-1)*dw
    aw = abs(w)
    arg = (2*abs(aw-wc)/bw)^twod
    filter_neg = (amp2 * exp(-arg))

### normalize window maximum to 1 (rather than 1/sqrt(2)
    filter_pos=filter_pos/max(filter_pos)
    filter_neg=filter_neg/max(filter_neg)

    taper=c(filter_pos,filter_neg)   

### apply filter
    ft = ft*taper

# inverse FFT
###  convert to real number and normalize iFFT output by length of record, supress warning
###   about discarding imaginary component (it is zero!).
    ifft <- suppressWarnings(as.double(fft(ft,inverse=TRUE)/nf))
### isolate prepadded portion
    bp <- ifft[1:npts]
    d<- data.frame(cbind(((1:npts)-1)*dt, bp))
    colnames(d) <- c("Location","Value")
       
    if(genplot) plot(d, type="l",col="red", ylab="Value", xlab="Location", main="Bandpassed Signal",bty="n")
   
    if(output == 1) {return(d)}     
    if(output == 2) return(data.frame( cbind(freq[1:nyqfreq],filter_pos) ) )
       
#### END function taner
}