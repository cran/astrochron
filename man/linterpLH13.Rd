\name{linterpLH13}
\alias{linterpLH13}
\title{Piecewise linear interpolation of stratigraphic series using the approach of Laepple and Huybers (2013)}
\description{
Linearly interpolate stratigraphic series using the approach of Laepple and Huybers (2013), including anti-alias filtering.
}
\usage{
linterpLH13(dat,dt=NULL,start=NULL,dtMin=NULL,thresh=0.7,beta=1,tbw=20,smooth=0.05,
            logF=F, antialias=T,roll=10^10,nsim=1,ncores=2,output=T,maxF=NULL,pl=1,
            genplot=T,check=T,verbose=1)
}
\arguments{
  \item{dat}{Stratigraphic series for piecewise linear interpolation. First column should be location (e.g., depth), second column should be data value.}
  \item{dt}{New optimal sampling interval. If dt is not specified, it will be estimated using the approach of Laepple and Huybers (2013).}
  \item{start}{Start interpolating at what time/depth/height value? By default, the first value of the stratigraphic series will be used.}
  \item{dtMin}{Very fine sampling interval for noise generation. If dt is not specified, the finest sample spacing of dat is used.}
  \item{thresh}{Threshold value for R2 (see details). Default value is 0.7 as in Laepple and Huybers (2013)}
  \item{beta}{Power law coefficient for stochastic surrogates, >=0.}
  \item{tbw}{MTM time-bandwidth product.}
  \item{smooth}{Smoothing parameter for spectral ratio R (see details).}
  \item{logF}{Smooth ratio R using log10 transformed frequency? (T or F)}  
  \item{antialias}{Apply anti-alias filter? (T or F)}  
  \item{roll}{Taner filter roll-off rate for anti-alias filter, in dB/octave.}  
  \item{nsim}{Number of Monte Carlo simulations.}  
  \item{ncores}{Number of cores to use for parallel processing.}  
  \item{maxF}{Maximum frequency for spectrum plotting.}  
  \item{pl}{Plot logarithm of spectral power (1) or linear spectral power (2)?}  
  \item{output}{Return interpolated series as new data frame? (T or F)}  
  \item{genplot}{Generate summary plots? (T or F)}
  \item{check}{Conduct compliance checks before processing? (T or F). In general this should be activated; the option is included for Monte Carlo simulation.}
  \item{verbose}{Verbose output? (0=none, 1=some reporting, 2=extensive reporting)}
}

\details{
This function seeks to estimate the optimal sampling interval for piecewise linear interpolation of an unevenly sampled stratigraphic series ("D"), so as to reduce bias when
reconstructing the spectral background.  The algorithm is slightly modified from Laepple and Huybers (2013), and consists of the following steps:

(1) Generate a 1/f stochastic noise signal, N. A sampling interval equivalent to the finest step of data set D is used for the noise model if dtMin is NULL.

(2) Subsample the 1/f noise (N) on the actual data set D sampling grid, yielding N2.

(3) Interpolate the subsampled 1/f noise (N2) to dtMin, yielding N3.

(4) Calculate the MTM power spectra for N (Spec_N) and N3 (Spec_N3), using a low resolution spectrum to provide a relatively smooth estimate of the continuum.

(5) Divide the resampled power spectrum (Spec_N3) by unresampled spectrum (Spec_N), yielding the ratio R.

(6) Smooth the ratio R using a Gaussian kernel, yielding R2. 

(7) Estimate the highest reliable frequency (F) as the lowest frequency at which the ratio R2 < 0.7 (or specified thresh value).

(8) Calculate the optimal interpolation resolution, dt, as 0.5/F.

If anti-alias filtering is selected:

(9) Interpolate the data set D to 0.1*dt, yielding D2.

(10) Filter D2 using a Taner lowpass (anti-alias) filter with a cutoff frequency of 1.2/dt, yielding D3.

(11) Resample D3 at the optimal resolution dt, yielding D4.

If anti-alias filtering is not selected, instead of steps 9-11, resample D at the optimal resolution dt.

When anti-alias filtering is applied, the first and last value of the interpolated series (D4) 
have a tendency to be strongly biased if they closely align with the original sampling grid of D. 
To address this issue, the first and last interpolated points are removed if they are 
within 0.5*dt (the optimal interpolation interval) of the first or last sampling location of D.

NOTE: A finer (smaller increment) interpolation interval than the 'optimal' one proposed by linterpLH13 
may be appropriate for astronomical signal detection, although the magnitude of the observed astronomical 
variability may be substantially underestimated. 
}

\references{
T. Laepple and P. Huybers, 2013,
\emph{Reconciling discrepancies between Uk37 and Mg/Ca reconstructions of Holocene marine temperature variability}: Earth and Planetary Science Letters, v. 375, p. 418-429.
}

\examples{
\donttest{
# generate example series from ar1 noise, 5 kyr sampling interval
ex = ar1(npts=501,dt=10)

# jitter sampling times
ex[1]=ex[1]+rnorm(501,sd=1)
# sort
ex = ex[order(ex[,1],na.last=NA,decreasing=FALSE),]

# view sampling statistics
strats(ex)

# run linterpLH13
linterpLH13(ex,output=FALSE)

# run linterpLH13 with multiple simulations
linterpLH13(ex,nsim=500,output=FALSE)
}
}