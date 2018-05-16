### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2018 Stephen R. Meyers
###
###########################################################################
### makeNoise: make noise given specified power spectrum background, 
###             following the algorithm of Timmer and Konig  
###             (SRM: December 2-4, 2017; April 13, 2018)
### this function is called by pwrLaw and BPL
###########################################################################

makeNoise <- function (S,dt=1,mean=0,sdev=1,addPt=F,nsim=1,genplot=T,verbose=T)
{

  S=as.numeric(S)
  nfreq=length(S)
# number of points for time series
  npts=(nfreq-1)*2

  if(verbose) cat("\n ----- GENERATING SURROGATE SERIES FROM SPECTRUM-----\n")

# Nyquist
   Nyq <- 1/(2*dt)
# Rayleigh Frequency
   Ray <- 1/(npts*dt)
   df <- Ray
# FFT: Locate real components for zero, nyquist, first neg freq., (zero-df)
   izero = 1
   nyqfreq = 0.5*npts + 1
   negfreq = 0.5*npts + 2
   minusdf = npts

# initialize matrix for simulation results
   if(!addPt) 
     {
       noise2 <- double(npts*nsim)
       dim(noise2) <- c(npts,nsim)
     }
   if(addPt) 
     {
       noise2 <- double((npts-1)*nsim)
       dim(noise2) <- c(npts-1,nsim)
     }

# start simulation loop
for (i in 1:nsim)
 { 
# Intialize arrays
   f <-double(npts)
   a <-double(npts)
   b <-double(npts)
# f[1],  a[1] and b[1] will stay as zero
   
### For each postive frequency (f), draw two Gaussian distributed random numbers,
#    and multiply them by the square root of the power spectrum.  
#    The first random number will represent the real and the second the imaginary
#    Fourier coefficient. Assign f(0) = 0, only generate real component for f(Nyq).
#    Make sure the frequencies are appropriately ordered:
#    f(0) -> f(Nyq); f(-Nyq+df) -> f(0-df)

   f[2:nyqfreq] = (1:(nyqfreq-1))*df
   a[2:(nyqfreq-1)] = rnorm(nyqfreq-2) * sqrt(S[2:(nyqfreq-1)])
   b[2:(nyqfreq-1)] = rnorm(nyqfreq-2) * sqrt(S[2:(nyqfreq-1)])
#  Now f(Nyq)
   f[nyqfreq] = Nyq
   a[nyqfreq] = rnorm(1) * sqrt(S[nyqfreq])
   b[nyqfreq] = 0

### Assign Fourier coefficients to each negative frequency (f).
#    These should be assigned as the complex conjugate of the positive frequencies.
   f[negfreq:minusdf] = ( (npts/2) - (1:length(negfreq:minusdf)) ) * df * -1
   a[negfreq:minusdf] = a[(nyqfreq-1):2]
   b[negfreq:minusdf] = -1 * b[(nyqfreq-1):2]
 
### Calculate power spectrum (positive frequencies only) for plotting
if(genplot && nsim==1)
 {
   pwr <- double(nyqfreq)
   pwr[1:nyqfreq] <- a[1:nyqfreq]^2 +b[1:nyqfreq]^2
 }
    
###  Now perform the inverse FFT to obtain time series simulation
###  put Fourier coefficients into complex variable for iFFT
   x <- complex(npts)
   x <- a + 1i * b
###  convert to real number and normalize iFFT output by length of record
   noise1 <- suppressWarnings(as.double( fft(x,inverse=TRUE)/npts ))

### if addPt=T, remove the extra data point
   if(addPt) 
     {
       npts <- npts-1
       noise1 <- noise1[1:npts]
     }     
 
### standardize to specified mean and variance
###  note, first mean is function call, second instance is desired mean value
   noise1 = noise1 - mean(noise1)
   noise1 = noise1 * sdev/sd(noise1)
   noise1 = noise1 + mean
   noise2[,i] = noise1
}

### now make time axis
   time <-seq(0,(npts-1)*dt,by=dt)
   noise <- data.frame(cbind(time,noise2))

   if(genplot && nsim==1)
    {
      par(mfrow=c(3,2))
      suppressWarnings(plot(f[1:nyqfreq],pwr,type="l",xlab="Log10(Frequency)", ylab="Log10(Power)", main="Log-Log Power Spectrum of Surrogate Series",log="xy"))
      plot(f[1:nyqfreq],pwr,type="l",xlab="Frequency", ylab="Power", main="Power Spectrum of Surrogate Series")
      plot(noise, cex=.5, xlab="Location", ylab="Noise Value", main="Surrogate Noise Series"); lines(noise)
### plot the denisty and the histogram together
      hist(noise[,2],freq=F,xlab="Noise Value",main="Histogram of Noise Values"); lines(density(noise[,2], bw="nrd"),col="red"); grid()
### boxplot
      boxplot(noise[,2],ylab="Noise Value",main="Boxplot of Noise Values")
### Normal probabilty plot (Normal Q-Q Plot)
      qqnorm(noise[,2],xlab="Noise Value"); qqline(noise[,2], col="red"); grid()
     }
      
   return(noise)

#### END function makeNoise
}