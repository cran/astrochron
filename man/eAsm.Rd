\name{eAsm}
\alias{eAsm}
\title{Evolutive Average Spectral Misfit}
\description{
Calculate Evolutive Average Spectral Misfit with Monte Carlo spectra simulations, as updated in Meyers et al. (2012).
}
\usage{
eAsm(spec,siglevel=0.9,target,fper=NULL,rayleigh,nyquist,sedmin=1,sedmax=5,
      numsed=50,linLog=1,iter=100000,ydir=1,output=4,genplot=F)
}
\arguments{
  \item{spec}{Time-frequency spectral results to evaluate. Must have the following format: column 
  1=frequency; remaining columns (2 to n)=probability; titles for columns 2 to 
  n must be the location (depth or height). Note that this format is ouput by function eha.}
  \item{siglevel}{Threshold level for filtering peaks.}
  \item{target}{A vector of astronomical frequencies to evaluate (1/ka). These must be in order of increasing frequency (e.g., e1,e2,e3,o1,o2,p1,p2). Maximum allowed is 50 frequencies.}
  \item{fper}{A vector of uncertainties on each target frequency (1/ka).  Values should be from 0-1, representing uncertainty as a percent of each target frequency.  The order of the uncertainties must follow that of the target vector. By default, no uncertainty is assigned.}
  \item{rayleigh}{Rayleigh frequency (cycles/m).}
  \item{nyquist}{Nyquist frequency (cycles/m).}
  \item{sedmin}{Minimum sedimentation rate for investigation (cm/ka).}
  \item{sedmax}{Maximum sedimentation rate for investigation (cm/ka).}
  \item{numsed}{Number of sedimentation rates to investigate in ASM optimization grid. Maximum allowed is 500.}
  \item{linLog}{Use linear or logarithmic scaling for sedimentation rate grid spacing? (0=linear, 1=log)}
  \item{iter}{Number of Monte Carlo simulations for significance testing. Maximum allowed is 100,000.}
  \item{ydir}{Direction for y-axis in plots (depth or height). -1 = values increase downwards (slower plotting!), 1 = values increase upwards.}
  \item{output}{Return output as a new data frame? (0 = nothing, 1 = Ho-SL, 2 = ASM, 3 = # astronomical terms, 4 = everything)}
  \item{genplot}{Generate summary plots? (T or F)}
}

\details{
Please see function asm for details.
}

\seealso{ \code{\link{asm}}, \code{\link{eAsmTrack}}, \code{\link{eha}}, \code{\link{testPrecession}}, \code{\link{timeOpt}}, and \code{\link{timeOptSim}} }

\references{
S.R. Meyers and  B.B. Sageman, 2007,
\emph{Quantification of Deep-Time Orbital Forcing by Average Spectral Misfit}: American Journal of Science, v. 307, p. 773-792.

S.R. Meyers, 2012,
\emph{Seeing Red in Cyclic Stratigraphy: Spectral Noise Estimation for Astrochronology}: Paleoceanography, 27, PA3228, doi:10.1029/2012PA002307.

S.R. Meyers, B.B. Sageman and M.A. Arthur, 2012,
\emph{Obliquity forcing of organic matter accumulation during Oceanic Anoxic Event 2}: Paleoceanography, 27, PA3212, doi:10.1029/2012PA002286.
}

\examples{
\donttest{
# use modelA as an example
data(modelA)

# interpolate to even sampling interval
modelAInterp=linterp(modelA)

# perform EHA analysis, save harmonic F-test confidence level results to 'spec'
spec=eha(modelAInterp,win=8,step=2,pad=1000,output=4)

# perform Evolutive Average Spectral Misfit analysis, save results to 'res'
res=eAsm(spec,target=c(1/405.47,1/126.98,1/96.91,1/37.66,1/22.42,1/18.33),rayleigh=0.1245274,
         nyquist=6.66597,sedmin=0.5,sedmax=3,numsed=100,siglevel=0.8,iter=10000,output=4)

# identify minimum Ho-SL in each record and plot
pl(1)
eAsmTrack(res[1],threshold=0.05)

# extract Ho-SL result at 18.23 m
HoSL18.23=extract(res[1],get=18.23,pl=1)

# extract ASM result at 18.23 m
asm18.23=extract(res[2],get=18.23,pl=0)
}
}