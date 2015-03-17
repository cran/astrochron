\name{astrochron-package}
\alias{astrochron-package}
\alias{astrochron}
\docType{package}
\title{astrochron: A Computational Tool for Astrochronology}
\description{
'astrochron' is a computational tool for astrochronology. It includes routines for astrochronologic testing, astronomical time scale construction, and time series analysis. Also included are a range of statistical analysis and modeling routines that are relevant to time scale development and paleoclimate analysis.
}
\details{
\tabular{ll}{
Package: \tab astrochron\cr
Type: \tab Package\cr
Version: \tab 0.4.1\cr
Date: \tab 2015-03-17\cr
License: \tab GPL-3\cr
}
%%~~ An overview of how to use the package, including the most important functions ~~
}
\author{
Stephen Meyers

Maintainer: Stephen Meyers <smeyers@geology.wisc.edu>
}

\note{
Please note that this version of 'astrochron' is undergoing BETA TESTING.

Development of 'astrochron' is partially supported by the U.S. National Science Foundation:

\emph{CAREER: Deciphering the Beat of a Timeless Rhythm - The Future of Astrochronology (EAR 1151438 to S. Meyers)}

\emph{Collaborative Research: Evolution of the Climate Continuum- Late Paleogene to Present (OCE 1003603 to S. Meyers and L. Hinnov)}

TO CITE PACKAGE 'astrochron' IN PUBLICATIONS, PLEASE USE:

Meyers, S.R. (2014). \emph{astrochron: An R Package for Astrochronology}. http://cran.r-project.org/package=astrochron 

Also cite the original research papers that document the relevant algorithms, as referenced on the help pages for specific functions.

}

%%~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation directory ~~
\keyword{ package }
\examples{
%%~~ simple examples of the most important functions ~~
### EXAMPLES OF SOME FUNCTIONS AVAILABLE IN 'astrochron':

### This demo will use a model (series are usually read using the function 'read').
data(modelA)

### Interpolate the model stratigraphic series to its median sampling interval
modelAInterp=linterp(modelA)

### Calculate MTM spectrum using 2pi Slepian tapers, include AR1 condfidence level estimates,
###  plot power with linear scale
mtm(modelAInterp,tbw=2,ar=TRUE,pl=2)

### Perform Evolutive Harmonic Analysis using 2pi Slepian tapers, a window of 8 meters,
###  pad to 1000 points, and output Harmonic F-test confidence level results
fCL=eha(modelAInterp,win=8,pad=1000,output=4)

### Extract Harmonic F-test spectrum at approximately 22 meters height
spec=extract(fCL,22)
### In this extracted spectrum, identify F-test peak maxima exceeding 90% confidence level
freqs=peak(spec,level=0.9)[2]
### Conduct ASM testing on these peaks
# set Rayleigh frequency in cycles/m
rayleigh=0.1245274
# set Nyquist frequency in cycles/m
nyquist=6.66597
# set orbital target in 1/ky
target=c(1/405.47,1/126.98,1/96.91,1/37.66,1/22.42,1/18.33)
# execute ASM
asm(freq=freqs,target=target,rayleigh=rayleigh,nyquist=nyquist,sedmin=0.5,sedmax=3,numsed=100,
    linLog=1,iter=100000,output=FALSE)

# Check to see if this is an interactive R session, for compliance with CRAN standards.
# YOU CAN SKIP THE FOLLOWING LINE IF YOU ARE USING AN INTERACTIVE SESSION.
if(interactive()) {

### Interactively track obliquity term in EHA harmonic F-test confidence level results
freqs=trackFreq(fCL,fmin=1.2,fmax=2.4,threshold=0.8)

### Convert the spatial frequencies to sedimentation rates
sedrate=freq2sedrate(freqs,period=37.66)

### Convert the sedimentation rate curve to a time-space map
time=sedrate2time(sedrate)

### Tune the stratigraphic series using the time-space map
modelATuned=tune(modelAInterp,time)

### Interpolate the tuned series
modelATunedInterp=linterp(modelATuned)

### Perform Evolutive Harmonic Analysis on the tuned series
eha(modelATunedInterp)

}

}
