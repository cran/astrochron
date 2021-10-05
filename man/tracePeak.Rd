\name{tracePeak}
\alias{tracePeak}
\title{A tool to interactively trace peak trajectories on plots}
\description{
A tool to interactively trace peak trajectories on plots, for results from such functions as eTimeOpt, EHA, eAsm.
}
\usage{
tracePeak(dat,color=2,h=6,w=4,xmin=NULL,xmax=NULL,ymin=NULL,ymax=NULL,
          ydir=-1,palette=6,ncolors=100,path=1)
}
\arguments{
  \item{dat}{Data frame with results to evaluate. It must have the following format: column 
  1=parameter to track (e.g., frequency, sedimentation rate, etc.; x-axis of plot); remaining columns 
  (2 to n)=parameter to evaluate for peak identification (color on plot); titles for columns 2 to n 
  must be the location (depth/height/time; y-axis of plot). Note that this format is ouput by functions 
  eha, eTimeOpt, eAsm.}
  \item{color}{Line color for tracing. 1 = transparent black; 2 = transparent white; 3 = transparent yellow}
  \item{h}{Height of plot in inches.}
  \item{w}{Width of plot in inches.}
  \item{xmin}{Minimum parameter value to plot.}
  \item{xmax}{Maximum parameter value to plot.}
  \item{ymin}{Minimum depth/height/time to plot.}
  \item{ymax}{Maximum depth/height/time to plot.}  
  \item{ydir}{Direction for y-axis in plots (depth/height/time). -1 = values increase downwards, 1 = values increase upwards.}
  \item{palette}{What color palette would you like to use? (1) rainbow, (2) grayscale, (3) blue, (4) red, (5) blue-white-red, (6) viridis}
  \item{ncolors}{Number of colors to use in plot.}
  \item{path}{How do you want to represent the path?: 1=lines and points; 2=lines; 3=points}
}

\seealso{ \code{\link{eha}} and \code{\link{eTimeOpt}} }
