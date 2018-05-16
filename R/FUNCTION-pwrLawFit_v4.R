### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2018 Stephen R. Meyers
###
###########################################################################
### Estimate power law (1/f) fit - (SRM: November 1-2, 2015; 
###                              November 29, 2017; December 4, 2017)
###
### This follows a recipe adapted from Vaughan (2005):
###  (1) use least squares to fit a line to log10 power/ log10 frequency
###  (2) estimate beta and N from the fit
###  (3) perform a bias correction of N
###  (4) estimate confidence levels from bias corrected background
###
###  One can define the frequency range for fitting. 
###  For MTM, avoid the lowest frequences in the fit, which are biased.
###########################################################################

pwrLawFit <- function (spec,dof=2,flow=NULL,fhigh=NULL,output=1,genplot=T,verbose=T)
{

# spec contains linear frequency (column 1) and linear power (column 2).
# dof is the degrees of freedom for power, at each frequency. default is 2 (periodogram)
# flow : values >= flow are investigated in fit
# fhigh : values <= fhigh are investigated in fit
# this function assumes that f(0) and Nyquist are not included in input spectrum.

if(verbose) cat("\n----- ESTIMATING Power Law (1/f) Fit -----\n")

spec <- data.frame(spec)
nfreq <- length(spec[,1])
df <- spec[2,1]-spec[1,1]

# error checking 
   if(df<0)
     { 
       if (verbose) cat("\n * Sorting spectrum into increasing frequency, removing empty entries\n")
       spec <- spec[order(spec[1], na.last = NA, decreasing = F), ]
       df <- spec[2,1]-spec[1,1]
       nfreq <- length(spec[,1])
     }

if (verbose) 
 {
   cat(" * Number of frequencies in spectrum:",nfreq,"\n")
   cat(" * Minimum frequency:",spec[1,1],"\n")
   cat(" * Maximum frequency:",spec[nfreq,1],"\n")
   cat(" * Frequency spacing:",df,"\n")
 }

# power law fit and confidence level estimation
if (is.null(flow)) flow=spec[1,1]
if (is.null(fhigh)) fhigh=spec[nfreq,1]

# isolate frequencies for fitting
specfit <- subset(spec, (spec[1] >= flow) & (spec[1] <= fhigh) )

# fit line to log10(power) and log10(freq), which is the equivalent
#  of an exponential continuum model
#   predict power given frequency
lm.0 <- lm(log10(specfit[,2]) ~ log10(specfit[,1]))

if (verbose) 
 {
   cat(" * Slope from log(power) vs. log(frequency) fit (m):",lm.0$coefficients[2],"\n")
   cat(" * Y-intercept from log(power) vs. log(frequency) fit (b):",lm.0$coefficients[1],"\n")
 }

# save linear fit, which is in log power
lm.0.fit= as.vector(lm.0$fit)

# see Vaughan (2005, pg. 3) and Abramowitz & Stegun (1964): Here we want the 
#  expectation value of the log of the spectrum, which is not the expectation 
#  of log(spectrum). We must apply a bias correction for y = mx + b, 
#  and the power law relationship P=Nf^-beta
#  beta = - m 
beta =-1*lm.0$coefficients[2]
# log(N) = b - bias
bias = (digamma(dof/2)-log(dof/2))/log(10)
logN = lm.0$coefficients[1] - bias

if (verbose) 
 {
   cat(" * beta =",beta,"\n")
   cat(" * log(N) =",logN,"\n")
   cat(" * estimated bias =",bias,"\n")
 }

# this is the bias corrected power law fit, as log10(power)
PLfit = logN - beta*log10(specfit[,1])

# convert PLfit from log10 power to power
PLfit=10^PLfit

chiPL <- (specfit[,2]/PLfit) * dof
chiCL <- pchisq(chiPL, df=dof)
### 90, 95 and 99% confidence levels
CL_90 <- PLfit*qchisq(0.9, df=dof)/dof
CL_95 <- PLfit*qchisq(0.95, df=dof)/dof
CL_99 <- PLfit*qchisq(0.99, df=dof)/dof

### generate plots
if(genplot)
 {
   par(mfrow=c(2,1))
   mtitle=c("1/f fit (blue), unbiased 1/f fit (red), 90%CL, 95%CL, 99%CL (dotted)")
     plot(spec[,1],spec[,2],type="l", col="black", main=mtitle,xlab="Log Frequency",ylab="Log Power",cex.axis=1.1,cex.lab=1.1,lwd=2,bty="n",log="xy")
      lines(specfit[,1],10^lm.0.fit,col="blue",lwd=2)
      lines(specfit[,1],PLfit,col="red",lwd=2)
      lines(specfit[,1],CL_90,col="red",lwd=1,lty=3)
      lines(specfit[,1],CL_95,col="red",lwd=1,lty=3)
      lines(specfit[,1],CL_99,col="red",lwd=1,lty=3)
     plot(spec[,1],spec[,2],type="l", col="black", main=mtitle,xlab="Linear Frequency",ylab="Linear Power",cex.axis=1.1,cex.lab=1.1,lwd=2,bty="n")
      lines(specfit[,1],10^lm.0.fit,col="blue",lwd=2)
      lines(specfit[,1],PLfit,col="red",lwd=2)
      lines(specfit[,1],CL_90,col="red",lwd=1,lty=3)
      lines(specfit[,1],CL_95,col="red",lwd=1,lty=3)
      lines(specfit[,1],CL_99,col="red",lwd=1,lty=3)
 }

if (output==1) 
 {
   spectrum <- data.frame(cbind(specfit[,1],specfit[,2],chiCL*100,PLfit,CL_90,CL_95,CL_99))
   colnames(spectrum)[1] <- 'Frequency'
   colnames(spectrum)[2] <- 'Power'
   colnames(spectrum)[3] <- 'PowerLaw_CL'
   colnames(spectrum)[4] <- 'PowerLaw_fit'
   colnames(spectrum)[5] <- 'CL_90'
   colnames(spectrum)[6] <- 'CL_95'
   colnames(spectrum)[7] <- 'CL_99'

   return(spectrum)
 }

if (output==2) 
 {
   out=data.frame(cbind(beta,logN,lm.0$coefficients[1]))
   colnames(out) <- c("beta","unbiased_logN","biased_logN")
   rownames(out) <- NULL
   return(out)
 }
   
#### END function pwrLawFit
}
