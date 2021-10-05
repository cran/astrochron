\name{integratePower}
\alias{integratePower}
\title{Determine the total power within a given bandwidth}
\description{
Determine the total power within a given bandwidth, and also the ratio of this power to the total power in the spectrum (or up to a specified frequency). If bandwidth is not specified, generate interactive plots for bandwidth selection. For use with the function eha, integratePower can process spectrograms (time-frequency) or single spectra.
}
\usage{
integratePower(spec,flow=NULL,fhigh=NULL,fmax=NULL,unity=F,f0=T,xmin=NULL,
               xmax=NULL,ymin=NULL,ymax=NULL,npts=NULL,pad=NULL,ydir=1,
               palette=6,ncolors=100,h=6,w=9,ln=F,genplot=T,verbose=T)
}
\arguments{
  \item{spec}{Spectral results to evaluate. If the data frame contains time-frequency results,
  it must have the following format: column 1=frequency; remaining columns (2 to n)=power; titles for columns 2 to 
  n must be the location (depth or height). Note that this format is ouput by function eha. If the data frame contains
  one spectrum, it must have the following format: column 1=frequency, column 2=power.}
  \item{flow}{Low frequency cutoff for integration. If flow or fhigh are not specified, interactive plotting is activated.}
  \item{fhigh}{High frequency cutoff for integration. If flow or fhigh are not specified, interactive plotting is activated.}
  \item{fmax}{Integrate total power up to this frequency.}
  \item{unity}{Normalize spectra such that total variance (up to fmax) is unity. (T of F)}
  \item{f0}{Is f(0) included in the spectra? (T or F)}
  \item{xmin}{Minimum frequency for PLOTTING.}
  \item{xmax}{Maximum frequency for PLOTTING.}
  \item{ymin}{Minimum depth/height/time for PLOTTING. Only used if processing time-frequency results.}
  \item{ymax}{Maximum depth/height/time for PLOTTING. Only used if processing time-frequency results.} 
  \item{npts}{The number of points in the processed time series window. This is needed for proper spectrum normalization.} 
  \item{pad}{The total padded length of the processed time series window. This is needed for proper spectrum normalization.} 
  \item{ydir}{Direction for y-axis in plots (depth or height). -1 = values increase downwards (slower plotting!), 1 = values increase upwards. Only used if processing time-frequency results.} 
  \item{palette}{What color palette would you like to use? (1) rainbow, (2) grayscale, (3) blue, (4) red, (5) blue-white-red, (6) viridis}
  \item{ncolors}{Number of colors to use in plot. Only used if processing time-frequency results.}
  \item{h}{Height of plot in inches.}
  \item{w}{Width of plot in inches.}  
  \item{ln}{Plot natural log of spectral results? (T or F)}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}

\details{
Depending on the normalization used, you may want to preprocess the power spectra prior to integration.
}

\seealso{ \code{\link{eha}} }

\examples{
# generate etp signal over past 10 Ma
ex=etp(tmax=10000)

# evolutive power
pwr=eha(ex,win=500,fmax=.1,pad=2000,output=2,pl=2)

# integrate power from main obliquity term
integratePower(pwr,flow=0.02,fhigh=0.029,npts=501,pad=2000)
}