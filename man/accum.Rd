\name{accum}
\alias{accum}
\title{Calculate geochemical proxy accumulation rates}
\description{
Calculate geochemical proxy accumulation rates using defined sedimentation rate and dry bulk density. Sedimentation rate and dry bulk density data are piecewise linearly interpolated to data locations, to estimate bulk accumulation rate.
}
\usage{
accum(dat,sedrate=NULL,density=NULL,genplot=T,check=T,verbose=T)
}
\arguments{
  \item{dat}{Stratigraphic series containing geochemical proxy data. First column should be location (e.g., depth), second column should be geochemical concentration in weight percent (0-100).}
  \item{sedrate}{Sedimentation rate. First column should be location (e.g., depth), second column should be sedimentation rate (cm/ka). Alternatively, if a single value is given, a constant sedimentation rate is applied.}
  \item{density}{Dry bulk density. First column should be location (e.g., depth), second column should be dry bulk density (g/cm3/ka). Alternatively, if a single value is given, a constant bulk density is applied.}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{check}{Conduct compliance checks before processing? (T or F). In general this should be activated; the option is included for Monte Carlo simulation.}
  \item{verbose}{Verbose output? (T or F)}
}

\examples{
# Generate an example geochemical data set, using a stochastic model
    ex=ar1(npts=1001)
# the first column in 'ex' is depth (m), and we will take the second column to be CaCO3
#  rescale second column to range from 0-100% CaCO3
    m=100/(max(ex[2])-min(ex[2]))
    b=100-(m*max(ex[2]))
    ex[2]= m*ex[2] +b
    autoPlot(ex)
# generate example sedimentation rate history, with a stepwise tripling 
#  from 1m/kyr to 3m/kyr at 500 m
    sr=ex
    sr[1:500,2]=1
    sr[501:1001,2]=3
    autoPlot(sr)

# calculate accumulation rates
    res=accum(dat=ex,sedrate=sr,density=2.65)   
}