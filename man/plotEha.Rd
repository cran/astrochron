\name{plotEha}
\alias{plotEha}
\title{Create color time-frequency plots from eha results.}
\description{
Create color time-frequency plots from eha results.
}
\usage{
plotEha(spec,xmin,xmax,ymin,ymax,h=6,w=4,ydir=1,pl=0,norm,palette=6,
        centerZero=T,ncolors=100,colorscale=F,xlab,ylab,filetype=0,output=T,verbose=T)
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
  \item{pl}{An option for the color plots (0=linear scale; 1=plot log of value [useful for plotting power], 2=normalize to maximum value in each window [useful for plotting amplitude], 3=use normalization provided in norm.} 
  \item{norm}{Optional amplitude normalization divisor, consisting of a single column dataframe. This option is provided in case you'd like to normalize a set of EHA results using the same scheme (e.g., before and after removal of spectral lines).}
  \item{palette}{What color palette would you like to use? (1) rainbow, (2) grayscale, (3) blue, (4) red, (5) blue-white-red (if values are negative and positive, white is centered on zero), (6) viridis}
  \item{centerZero}{Center color scale on zero (use an equal number of postive and negative color divisions)? (T or F)}
  \item{ncolors}{Number of colors steps to use in palette.}
  \item{colorscale}{Include a color scale in the plot? (T or F)}
  \item{xlab}{Label for x-axis. Default = "Frequency"}
  \item{ylab}{Label for y-axis. Default = "Location"}
  \item{filetype}{Generate .pdf, .jpeg, .png or tiff file? (0=no; 1=pdf; 2=jpeg; 3=png; 4=tiff)}
  \item{output}{If amplitude is normalized (pl = 2), output normalization used? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}
\examples{
## as an example, evaluate the modelA
data(modelA)

## interpolate to even sampling interval of 0.075 m
ex1=linterp(modelA, dt=0.075)
  
## perform EHA with a time-bandwidth parameter of 2, using an 7.95 meter window, 0.15 m step, 
## and pad to 1000 points, output amplitude
res=eha(ex1,tbw=2,win=7.95,step=0.15,pad=1000,genplot=0,output=3)

# plot EHA amplitude, normalized to maximum value in each window
plotEha(res,xlab="Frequency (cycles/m)",ylab="Height (m)",pl=2)
}