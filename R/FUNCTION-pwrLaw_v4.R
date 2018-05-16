### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2018 Stephen R. Meyers
###
###########################################################################
### pwrLaw : make power law noise - (SRM: January 20-24, 2012; Nov. 17, 2012;
###                                       October 30, 2015; Nov. 1, 2015;
###                                       December 3-4, 2017; April 13, 2018)
### This function generates a power law time series following the algorithm
###  of Timmer and Konig (1995), On Generating Power Law Noise, Astronomy
###  and Astrophysics.  The power law is according to log(freq).
###  It allows a flat plateau at low frequencies
###########################################################################


pwrLaw <- function (npts=1024,dt=1,mean=0,sdev=1,beta=2,fcut=0,nsim=1,genplot=T,verbose=T)
{

   if(verbose) cat("\n ----- GENERATING POWER LAW SURROGATES -----\n")

   addPt=F
   if(npts%%2 != 0) 
     {
# Increasing npts by 1, so Nyquist exists
       npts=npts+1
       addPt=T
     }

#### Calculate Nyquist
   Nyq <- 1/(2*dt)
#### Calculate Rayleigh Frequency
   Ray <- 1/(npts*dt)
   df <- Ray
   nyqfreq = 0.5*npts + 1

# Intialize arrays
   f <-double(nyqfreq)
   S <-double(nyqfreq)
# f[1] and S[1] will stay as zero

# Calculate the theoretical BPL power spectrum. See eq B5 of Vaughan et al., 2011     
   f[2:nyqfreq] = (1:(nyqfreq-1))*df

# if plateau desired
   if(fcut>0)
    {
# identify frequencies below fcut, if plateau desired
      ifreq=which(f<fcut)   
      S[ifreq[-1]] = fcut^-beta
# generate power law section of spectrum      
      S[(max(ifreq)+1):nyqfreq] = f[(max(ifreq)+1):nyqfreq]^-beta
    }
# if no plateau desired
   if(fcut==0) S[2:nyqfreq] = f[2:nyqfreq]^-beta

   noise=makeNoise(S,dt=dt,mean=mean,sdev=sdev,addPt=addPt,nsim=nsim,genplot=genplot,verbose=F)
      
   return(noise)

#### END function pwrLaw
}