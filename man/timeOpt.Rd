\name{timeOpt}
\alias{timeOpt}
\title{TimeOpt: Evaluation of eccentricity-related amplitude modulation and bundling in paleoclimate data}
\description{
TimeOpt: Evaluation of eccentricity-related amplitude modulation and bundling in paleoclimate data, as in Meyers (2015).
}
\usage{
timeOpt(dat,sedmin=0.5,sedmax=5,numsed=100,linLog=1,limit=T,fit=1,r2max=1,
        fitModPwr=T,flow=NULL,fhigh=NULL,roll=NULL,targetE=NULL,targetP=NULL,
        detrend=T,output=0,title=NULL,genplot=T,check=T,verbose=T)
}
\arguments{
  \item{dat}{Stratigraphic series for astrochronologic assessment. First column should be depth or height (in meters), second column should be data value.}
  \item{sedmin}{Minimum sedimentation rate for investigation (cm/ka).}
  \item{sedmax}{Maximum sedimentation rate for investigation (cm/ka).}
  \item{numsed}{Number of sedimentation rates to investigate in optimization grid.}
  \item{linLog}{Use linear or logarithmic scaling for sedimentation rate grid spacing? (0=linear, 1=log; default value  is 1)}
  \item{limit}{Limit evaluated sedimentation rates to region in which full target signal can be recovered? (T or F)}
  \item{fit}{Test for (1) precession amplitude modulation or (2) short eccentricity amplitude modulation?}
  \item{r2max}{Which maximum in r2 do you want to use for the time model calibration? (1) r^2_opt, (2) r^2_spectral, (3) r^2_envelope}
  \item{fitModPwr}{Include the modulation periods in the spectral fit? (T or F)}
  \item{flow}{Low frequency cut-off for Taner bandpass (half power point; in cycles/ka)}
  \item{fhigh}{High frequency cut-off for Taner bandpass (half power point; in cycles/ka)}
  \item{roll}{Taner filter roll-off rate, in dB/octave.}  
  \item{targetE}{A vector of eccentricity periods to evaluate (in ka). These must be in order of decreasing period, with a first value of 405 ka.}
  \item{targetP}{A vector of precession periods to evaluate (in ka). These must be in order of decreasing period.}
  \item{detrend}{Remove linear trend from data series? (T or F)}
  \item{output}{Which results would you like to return to the console? (0) no output; (1) return sedimentation rate grid, r^2_envelope, r^2_power, r^2_opt; (2) return optimal time series, filtered precession, precession envelope, TimeOpt-reconstructed eccentricity model, full TimeOpt-regression model}
  \item{title}{A character string (in quotes) specifying the title for the graphics window (optional)}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{check}{Conduct compliance checks before processing? (T or F) In general this should be activated; the option is included for Monte Carlo simulation.}
  \item{verbose}{Verbose output? (T or F)}
}

\details{
TimeOpt is an astronomical testing algorithm for untuned (spatial) stratigraphic data. 
The algorithm identifies the sedimentation rate(s) that simultaneously optimizes: 
(1) eccentricity amplitude modulations within the precession band, and (2) the 
concentration of spectral power at specified target astronomical periods.

For each temporal calibration investigated (i.e., sedimentation rate), the observed precession 
band amplitude envelope is extracted using bandpass filtering and the Hilbert transform. 
The fit of the extracted precession envelope to the eccentricity periods is evaluated 
using a linear regression onto sine and cosine terms that reflect the five dominant 
eccentricity periods (~405.7, 130.7, 123.8, 98.9 and 94.9 kyr); amplitude and phase of the 
eccentricity terms are not assigned, but are determined during the linear model optimization. 
This approach is advantageous, as (1) the transfer functions associated with the climate 
and depositional systems can alter the amplitude and phase of the theoretical eccentricity 
terms (e.g, Laurin et al., 2005), and (2) the amplitude and phase of the eccentricity terms 
are unconstrained for deep-time investigations (>50 Ma). The quality of the "fit" is 
estimated by calculation of the correlation of the fitted eccentricity model time series 
to the observed precession band envelope (r^2_envelope), indicating the fraction of 
variance shared between the model and envelope.

The concentration of power at the target astronomical periods is evaluated using
a linear regression of the temporally-calibrated series onto sine and cosine terms 
that reflect the dominant eccentricity and precession periods. As above, the amplitude
and phase of each term is determined during the linear model optimization, and the
quality of the "fit" is estimated by calculation of the correlation of the fitted
astronomical model series to the temporally-calibrated series (r^2_spectral).

The final measure of fit (r^2_opt) is determined as: 

r^2_opt = r^2_envelope * r^2_spectral

which is simply the product of the fraction of variance shared between "model and envelope" 
and "model and time-calibrated data". This optimization approach identifies the sedimentation 
rate at which the precession envelope strongly expresses expected eccentricity modulation, 
while simultaneously, spectral power is concentrated at the target astronomical periods. 
r^2_opt can take on values ranging from 0 to 1 (a perfect fit to the astronomical model),
and provides a measure of overall quality of the astronomically calibrated time series.  
A similar approach is applicable to evaluate short eccentricity amplitude modulations. 
The statistical significance of the r^2_opt is determined via Monte Carlo simulation 
(see timeOptSim).

}

\value{
if output = 1, a data frame containing the following will be returned:
Sedimentation rate (cm/ka), r^2_envelope, r^2_spectral, r^2_opt

if output = 2, a data frame containing the following will be returned:
Time (ka), tuned time series, filtered precession, precession envelope, TimeOpt-reconstructed eccentricity model, full TimeOpt-regression model
}

\seealso{ \code{\link{asm}}, \code{\link{eAsmTrack}}, \code{\link{testPrecession}}, \code{\link{timeOptPlot}}, and \code{\link{timeOptSim}} }

\references{
S.R. Meyers, 2015,
\emph{The evaluation of eccentricity-related amplitude modulation and bundling in paleoclimate data: An inverse approach for astrochronologic testing and time scale optimization}: Paleoceanography, 30, doi:10.1002/2015PA002850.
}

\examples{
 \dontrun{
# generate a test signal with precession and eccentricity
ex=etp(tmin=1,tmax=1000,dt=5,pWt=1,oWt=0,eWt=1,esinw=TRUE,genplot=FALSE,verbose=FALSE)
# convert to meters with sedimentation rate of 2 cm/kyr
ex[1]<-ex[1]*0.02
timeOpt(ex,sedmin=0.5,sedmax=5,numsed=100,fit=1,output=0)
# evaluate short eccentricity modulations
timeOpt(ex,sedmin=0.5,sedmax=5,numsed=100,fit=2,output=0)
 }
}