\name{astrochron-package}
\alias{astrochron-package}
\alias{astrochron}
\docType{package}
\title{astrochron: A Computational Tool for Astrochronology}
\description{
This software provides routines for astrochronologic testing, astronomical time scale 
construction, and time series analysis <doi:10.1016/j.earscirev.2018.11.015>. Also 
included are a range of statistical analysis and modeling routines that are relevant 
to time scale development and paleoclimate analysis.
}
\details{
\tabular{ll}{
Package: \tab astrochron\cr
Type: \tab Package\cr
Version: \tab 1.2\cr
Date: \tab 2023-08-25\cr
License: \tab GPL-3\cr
}
%%~~ An overview of how to use the package ~~
}
\author{
Stephen Meyers

Maintainer: Stephen Meyers <smeyers@geology.wisc.edu>
}

\note{

Development of the 'astrochron' package is partially supported by the U.S. National Science Foundation and the Heising-Simons Foundation:

\emph{Leveraging the Geologic Record to Constrain Solar System Evolution, Earth-Moon Dynamics, Paleoclimate Change and Geological Time (Heising-Simons Foundation Award 2021-2797)}

\emph{Collaborative Research: Improving the Late Cretaceous-Eocene geomagnetic polarity time scale by integrating the global magnetic anomaly record and astrochronology (U.S. National Science Foundation Award OCE 2051616)}

\emph{CAREER: Deciphering the Beat of a Timeless Rhythm - The Future of Astrochronology (U.S. National Science Foundation Award EAR 1151438)}

\emph{Collaborative Research: Evolution of the Climate Continuum - Late Paleogene to Present (U.S. National Science Foundation Award OCE 1003603)}




TO CITE THIS PACKAGE IN PUBLICATIONS, PLEASE USE:

Meyers, S.R. (2014). \emph{Astrochron: An R Package for Astrochronology}. https://cran.r-project.org/package=astrochron 

Also cite the original research papers that document the relevant algorithms, as referenced on the help pages for specific functions.

}

%%~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation directory ~~
\keyword{ package }

\references{
S.R. Meyers, 2019, 
\emph{Cyclostratigraphy and the problem of astrochronologic testing}: Earth-Science Reviews
v.190, 190-223, doi:10.1016/j.earscirev.2018.11.015.
}

\examples{
%%~~ simple examples of the most important functions ~~
### EXAMPLES OF SOME FUNCTIONS AVAILABLE IN THIS SOFTWARE:

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

### Perform Average Spectral Misfit analysis
### Extract Harmonic F-test spectrum at approximately 22 meters height
spec=extract(fCL,22)
### In this extracted spectrum, identify F-test peak maxima exceeding 90% confidence level
freqs=peak(spec,level=0.9)[2]
### Conduct ASM testing on these peaks
# set Rayleigh frequency in cycles/m
rayleigh=0.1245274
# set Nyquist frequency in cycles/m
nyquist=6.66597
# set astronomical target in 1/ky
target=c(1/405.47,1/126.98,1/96.91,1/37.66,1/22.42,1/18.33)
# execute ASM
asm(freq=freqs,target=target,rayleigh=rayleigh,nyquist=nyquist,sedmin=0.5,sedmax=3,numsed=100,
    linLog=1,iter=100000,output=FALSE)

# Check to see if this is an interactive R session (for compliance with CRAN standards).
# YOU SHOULD SKIP THE FOLLOWING LINE IF YOU ARE USING AN INTERACTIVE SESSION.
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
