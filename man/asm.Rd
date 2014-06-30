\name{asm}
\alias{asm}
\title{Average Spectral Misfit}
\description{
Calculate Average Spectral Misfit with Monte Carlo spectra simulations, as updated in Meyers et al. (2012).
}
\usage{
asm(freq,target,fper=NULL,rayleigh,nyquist,sedmin=1,sedmax=5,numsed=50,
    linLog=1,iter=100000,output=F,genplot=T)
}
\arguments{
  \item{freq}{A vector of candidate astronomical cycles observed in your data spectrum (cycles/m). Maximum allowed is 500.}
  \item{target}{A vector of astronomical frequencies to evaluate (1/ka). These must be in order of increasing frequency (e.g., e1,e2,e3,o1,o2,p1,p2). Maximum allowed is 50 frequencies.}
  \item{fper}{A vector of uncertainties on each target frequency (1/ka).  Values should be from 0-1, representing uncertainty as a percent of each target frequency.  The order of the uncertainties must follow that of the target vector. By default, no uncertainty is assigned.}
  \item{rayleigh}{Rayleigh frequency (cycles/m).}
  \item{nyquist}{Nyquist frequency (cycles/m).}
  \item{sedmin}{Minimum sedimentation rate for investigation (cm/ka).}
  \item{sedmax}{Maximum sedimentation rate for investigation (cm/ka).}
  \item{numsed}{Number of sedimentation rates to investigate in ASM optimization grid. Maximum allowed is 500.}
  \item{linLog}{Use linear or logarithmic scaling for sedimentation rate grid spacing? (0=linear, 1=log)}
  \item{iter}{Number of Monte Carlo simulations for significance testing. Maximum allowed is 100,000.}
  \item{output}{Return output as a new data frame? (T or F)}
  \item{genplot}{Generate summary plots? (T or F)}
}

\details{
This function will caculate the Average Spectral Misfit between a data spectrum and astronomical target spectrum, following Meyers and Sageman (2007), 
and the improvements of Meyers et al. (2012).  
}

\value{
A data frame containing:
Sedimentation rate (cm/ka), ASM (cycles/ka), Null hypothesis significance level (0-100 percent), Number of astronomical terms fit.
}

\seealso{ \code{\link{eAsm}} }

\references{
S.R. Meyers and  B.B. Sageman, 2007,
\emph{Quantification of Deep-Time Orbital Forcing by Average Spectral Misfit}: American Journal of Science, v. 307, p. 773-792.

S.R. Meyers, B.B. Sageman and M.A. Arthur, 2012,
\emph{Obliquity forcing of organic matter accumulation during Oceanic Anoxic Event 2}: Paleoceanography, 27, PA3212, doi:10.1029/2012PA002286.
}

\examples{
## these frequencies are from modelA (type '?astrochron' for more information). Units are cycles/m  
freq <- c(0.1599833,0.5332776,1.5998329,2.6797201,3.2796575,3.8795948,5.5194235,6.5459830)
freq <- data.frame(freq)

## Rayleigh frequency in cycles/m
rayleigh <- 0.1245274

## Nyquist frequency in cycles/m
nyquist <- 6.66597

## orbital target in 1/ky. Predicted periods for 94 Ma (see Meyers et al., 2012)
target <- c(1/405.47,1/126.98,1/96.91,1/37.66,1/22.42,1/18.33)

## percent uncertainty in orbital target
fper=c(0.023,0.046,0.042,0.008,0.035,0.004)

asm(freq=freq,target=target,fper=fper,rayleigh=rayleigh,nyquist=nyquist,sedmin=0.5,sedmax=3,
    numsed=100,linLog=1,iter=100000,output=FALSE)

}