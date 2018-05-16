\name{trackPeak}
\alias{trackPeak}
\title{A tool to interactively select points to track peak trajectories on plots}
\description{
A tool to interactively select points to track peak trajectories on plots, for results from functions such as eTimeOpt, EHA, eAsm.
}
\usage{
trackPeak(dat,threshold=NULL,pick=T,minVal=NULL,maxVal=NULL,dmin=NULL,dmax=NULL,
          xmin=NULL,xmax=NULL,ymin=NULL,ymax=NULL,ydir=-1,h=6,w=4,ncolors=100,
          genplot=T,verbose=T)
}
\arguments{
  \item{dat}{Data frame with results to evaluate. It must have the following format: column 
  1=parameter to track (e.g., frequency, sedimentation rate, etc.; x-axis of plot); remaining columns 
  (2 to n)=parameter to evaluate for peak identification (color on plot); titles for columns 2 to n 
  must be the location (depth/height/time; y-axis of plot). Note that this format is ouput by functions 
  eha, eTimeOpt, eAsm.}
  \item{threshold}{Threshold level for filtering peaks. By default all peak maxima reported.}
  \item{pick}{Pick the peaks of interest using a graphical interface? (T or F). Only activated if genplot=T.}
  \item{minVal}{Minimum parameter value for analysis (e.g., frequency, sedimentation rate, etc.).}
  \item{maxVal}{Maximum parameter value for analysis  (e.g., frequency, sedimentation rate, etc.).}
  \item{dmin}{Minimum depth/height/time for analysis. NOT ACTIVATED YET!}
  \item{dmax}{Maximum depth/height/time for analysis. NOT ACTIVATED YET!}  
  \item{xmin}{Minimum parameter value for PLOTTING.}
  \item{xmax}{Maximum parameter value for PLOTTING.}
  \item{ymin}{Minimum depth/height/time for PLOTTING.}
  \item{ymax}{Maximum depth/height/time for PLOTTING.} 
  \item{ydir}{Direction for y-axis in plots (depth or height). -1 = values increase downwards (slower plotting!), 1 = values increase upwards.} 
  \item{h}{Height of plot in inches.}
  \item{w}{Width of plot in inches.}
  \item{ncolors}{Number of colors to use in plot.}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}

\seealso{ \code{\link{eha}} and \code{\link{eTimeOpt}} }