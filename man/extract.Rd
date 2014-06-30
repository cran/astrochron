\name{extract}
\alias{extract}
\title{Extract record from EHA time-frequency output or eAsm output}
\description{
Extract record from EHA time-frequency output or eAsm output: Use interactive graphical interface to identify record.
}
\usage{
extract(spec,get=NULL,xmin=NULL,xmax=NULL,ymin=NULL,ymax=NULL,h=6,w=4,ydir=-1,pl=0,
        ncolors=100,genplot=T,verbose=T)
}
\arguments{
  \item{spec}{Time-frequency spectral results to evaluate, or alternatively, eAsm results to evaluate. For time-frequency results, must have the following format: column 
  1=frequency; remaining columns (2 to n)=power, amplitude or probability; titles for columns 2 to 
  n must be the location (depth or height). Note that this format is ouput by function eha. For eAsm results, must have the following format: column 
  1=sedimentation rate; remaining columns (2 to n)=Ho-SL or ASM; titles for columns 2 to n must be the location (depth or height). }
  \item{get}{Record to extract (height/depth/time). If no value given, graphical interface is activated.}
  \item{xmin}{Minimum frequency or sedimentation rate for PLOTTING.}
  \item{xmax}{Maximum frequency or sedimentation rate for PLOTTING.}
  \item{ymin}{Minimum depth/height for PLOTTING.}
  \item{ymax}{Maximum depth/height for PLOTTING.} 
  \item{h}{Height of plot in inches.}
  \item{w}{Width of plot in inches.}
  \item{ydir}{Direction for y-axis in plots (depth or height). -1 = values increase downwards (slower plotting!), 1 = values increase upwards.} 
  \item{pl}{An option for the color plots (0=do nothing; 1=plot log of value [useful for plotting power], 2=normalize to maximum value [useful for plotting amplitude]).} 
  \item{ncolors}{Number of colors to use in plot.}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}

\seealso{ \code{\link{eha}} }