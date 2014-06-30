\name{plotEha}
\alias{plotEha}
\title{Create color time-frequency plots from eha results}
\description{
Create color time-frequency plots from eha results.
}
\usage{
plotEha(spec,xmin=NULL,xmax=NULL,ymin=NULL,ymax=NULL,h=6,w=4,ydir=1,pl=0,norm=NULL,
        xaxis=c("Frequency (cycles/ka)"),yaxis=c("Time (ka)"),ncolors=100,colorscale=F,
        filetype=0,output=T,verbose=T)
}
\arguments{
  \item{spec}{Time-frequency spectral results to evaluate. Must have the following format: column 
  1=frequency; remaining columns (2 to n)=power, amplitude or probability; titles for columns 2 to 
  n must be the location (depth or height). Note that this format is ouput by function eha.}
  \item{xmin}{Minimum frequency for PLOTTING.}
  \item{xmax}{Maximum frequency for PLOTTING.}
  \item{ymin}{Minimum depth/height for PLOTTING.}
  \item{ymax}{Maximum depth/height for PLOTTING.} 
  \item{h}{Height of plot in inches.}
  \item{w}{Width of plot in inches.}
  \item{ydir}{Direction for y-axis in plots (depth or height). -1 = values increase downwards (slower plotting!), 1 = values increase upwards.} 
  \item{pl}{An option for the color plots (0=do nothing; 1=plot log of value [useful for plotting power], 2=normalize to maximum value [useful for plotting amplitude], 3=use normalization provided in norm.} 
  \item{norm}{Optional amplitude normalization divisor, consisting of a single column dataframe. This option is provided in case you'd like to normalize a set of EHA results using the same scheme (e.g., before and after removal of spectral lines).}
  \item{xaxis}{Label for x-axis.}
  \item{yaxis}{Label for y-axis.}
  \item{ncolors}{Number of colors to use in plot.}
  \item{colorscale}{Include a color scale in the plot? (T or F)}
  \item{filetype}{Generate .pdf, .jpeg or .png file? (0=no; 1=pdf; 2=jpeg; 3=png)}
  \item{output}{If amplitude is normalized (pl = 2), output normalization used? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}