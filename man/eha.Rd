\name{eha}
\alias{eha}
\title{Evolutive Harmonic Analysis & Evolutive Power Spectral Analysis}
\description{
Evolutive Harmonic Analysis & Evolutive Power Spectral Analysis using the Thomson multitaper method (Thomson, 1982)
}
\usage{
eha(dat,tbw=2,pad,fmin,fmax,step,win,demean=T,detrend=T,siglevel=0.90,
    sigID=F,ydir=1,output=0,pl=1,palette=6,centerZero=T,ncolors=100,xlab,ylab,
    genplot=2,verbose=T)
}
\arguments{
  \item{dat}{Stratigraphic series to analyze. First column should be location (e.g., depth), second column should be data value.}
  \item{tbw}{MTM time-bandwidth product (<=10)}
  \item{pad}{Pad with zeros to how many points? Must not factor into a prime number >23. Maximum number of points is 200,000.}
  \item{fmin}{Smallest frequency for analysis and plotting.}
  \item{fmax}{Largest frequency for analysis and plotting.}
  \item{step}{Step size for EHA window, in units of space or time.}
  \item{win}{Window size for EHA, in units of space or time.}
  \item{demean}{Remove mean from data series? (T or F)}
  \item{detrend}{Remove linear trend from data series? (T or F)}
  \item{siglevel}{Significance level for peak identification/filtering (0-1)}
  \item{sigID}{Identify signficant frequencies on power, amplitude, and probabilty plots. Only applies when one spectrum is calculated. (T or F)}
  \item{ydir}{Direction for y-axis in EHA plots (depth,height,time). -1 = values increase downwards (slower plotting), 1 = values increase upwards}
  \item{output}{Return output as new data frame? 0=no; 1=all results; 2=power; 3=amplitude; 4=probability; 5=significant frequencies (only for one spectrum); 6=significant frequencies and their probabilities (only for one spectrum)}
  \item{pl}{Plot logarithm of spectral power (1) or linear spectral power (2)?}
  \item{palette}{What color palette would you like to use? (1) rainbow, (2) grayscale, (3) blue, (4) red, (5) blue-white-red (if values are negative and positive, white is centered on zero), (6) viridis}
  \item{centerZero}{Center color scale on zero (use an equal number of postive and negative color divisions)? (T or F)}
  \item{ncolors}{Number of colors steps to use in palette.}
  \item{xlab}{Label for x-axis. Default = "Frequency"}
  \item{ylab}{Label for y-axis. Default = "Location"}
  \item{genplot}{Plotting options. 0= no plots; 1= power, amplitude, f-test, probability; 2=data series, power, amplitude, probability; 3= data series, power, normalized amplitude (maximum in each window normalized to unity), normalized amplitude filtered at specified siglevel; 4= data series, normalized power (maximum in each window normalized to unity), normalized amplitude (maximum in each window normalized to unity), normalized amplitude filtered at specified siglevel}
  \item{verbose}{Verbose output? (T or F)}
}

\seealso{ \code{\link{extract}}, \code{\link{lowspec}}, \code{\link{mtmAR}}, \code{\link{mtmML96}}, \code{\link{periodogram}}, \code{\link{trackFreq}} and \code{\link{traceFreq}} }

\references{
Thomson, D. J., 1982, \emph{Spectrum estimation and harmonic analysis}, Proc. IEEE, 70, 1055-1096, doi:10.1109/PROC.1982.12433.
}

\examples{
## as an example, evaluate the modelA
data(modelA)

## interpolate to even sampling interval of 0.075 m
ex1=linterp(modelA, dt=0.075)
  
## perform EHA with a time-bandwidth parameter of 2, using an 7.95 meter window, 0.15 m step, 
## and pad to 1000 points
## set labels for plots (optional)
eha(ex1,tbw=2,win=7.95,step=0.15,pad=1000,xlab="Frequency (cycles/m)",ylab="Height (m)")

## for comparison generate spectrum for entire record, using time-bandwidth parameter of 3, and 
## pad to 5000 points
## start by making a new plot
pl(1)
eha(ex1,tbw=3,win=38,pad=5000,xlab="Frequency (cycles/m)")
}