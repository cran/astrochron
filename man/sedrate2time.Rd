\name{sedrate2time}
\alias{sedrate2time}
\title{Integrate sedimentation rate curve to obtain time-space map}
\description{
Integrate sedimentation rate curve to obtain time-space map.
}
\usage{
sedrate2time(sedrates,timedir=1,genplot=T,check=T,verbose=T)
}
\arguments{
  \item{sedrates}{Data frame containing depth/height in first column (meters) and sedimentation rates in second column (cm/ka).}
  \item{timedir}{Floating time scale direction: 1= time increases with depth/height; 2= time decreases with depth/height.)}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{check}{Conduct compliance checks before processing? (T or F) In general this should be activated; the option is included for Monte Carlo simulation.}
  \item{verbose}{Verbose output? (T or F)}
}