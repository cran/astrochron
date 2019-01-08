\name{timeOptTemplatePlot}
\alias{timeOptTemplatePlot}
\title{TimeOptTemplatePlot: Generate summary figure for TimeOptTemplate analyses}
\description{
TimeOptTemplatePlot: Generate summary figure for TimeOptTemplate analyses.
}
\usage{
timeOptTemplatePlot(dat=NULL,template=NULL,detrend=T,detrendTemplate=F,flipTemplate=F,
   srMin=NULL,srMax=NULL,res1=NULL,simres=NULL,fit=1,flow=NULL,fhigh=NULL,roll=NULL,
   targetE=NULL,targetP=NULL,xlab="Depth (m)",ylab="Proxy Value",fitR=NULL,output=0,
   verbose=T)
}
\arguments{
  \item{dat}{Stratigraphic series used for astrochronologic assessment. First column should be depth or height (in meters), second column should be data value.}
  \item{template}{Instantaneous sedimentation rate template to fit. This represents a unitless proportional sedimentation rate history. Default model is a copy of dat, which will be scaled for instantaneous accumulation optimization.}
  \item{detrend}{Remove linear trend from data series? (T or F)}
  \item{detrendTemplate}{Remove linear trend from sedimentation rate template? (T or F)}
  \item{flipTemplate}{Flip direction of sedimentation rate template? (T or F)}
  \item{srMin}{Minimum sedimentation rate for template}
  \item{srMax}{Maximum sedimentation rate for template}
  \item{res1}{Data frame containing TimeOpt results: sedimentation rate grid, r^2_envelope, r^2_power, r^2_opt.}
  \item{simres}{Data frame containing the r^2_opt value for each Monte Carlo simulation.}
  \item{fit}{Test for (1) precession amplitude modulation or (2) short eccentricity amplitude modulation?}
  \item{flow}{Low frequency cut-off for Taner bandpass (half power point; in cycles/ka).}
  \item{fhigh}{High frequency cut-off for Taner bandpass (half power point; in cycles/ka).}
  \item{roll}{Taner filter roll-off rate, in dB/octave.}
  \item{targetE}{A vector of eccentricity periods to evaluate (in ka). These must be in order of decreasing period, with a first value of 405 ka.}
  \item{targetP}{A vector of precession periods to evaluate (in ka). These must be in order of decreasing period.}
  \item{xlab}{Label for the depth/height axis.}
  \item{ylab}{Label for proxy variable evaluated.}
  \item{fitR}{The r2 value at the optimal sedimentation rate.}
  \item{output}{Which results you like to return to console? (0) no output; (1) return sedimentation rate grid, r2; (2) return optimal time series, bandpassed series, Hilbert and fitted periods}
  \item{verbose}{Verbose output? (T or F)}
}

\seealso{ \code{\link{asm}}, \code{\link{eAsmTrack}}, \code{\link{testPrecession}}, \code{\link{timeOpt}}, \code{\link{timeOptSim}}, and \code{\link{timeOptTemplate}}}

\references{
S.R. Meyers, 2015,
\emph{The evaluation of eccentricity-related amplitude modulations and bundling in 
paleoclimate data: An inverse approach for astrochronologic testing and time scale 
optimization}: Paleoceanography, v.30, 1625-1640.

S.R. Meyers, 2019, 
\emph{Cyclostratigraphy and the problem of astrochronologic testing}: Earth-Science Reviews
v. 190, 190-223.
}