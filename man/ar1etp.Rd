\name{ar1etp}
\alias{ar1etp}
\title{AR(1) + ETP simulation Routine}
\description{
Simulate a combined AR(1) + ETP signal, plot spectrum and confidence levels
}
\usage{
ar1etp(etpdat=NULL,iter=100,rho=0.9,ARvar=1,sig=90,tbw=2,padfac=5,ftest=F,fmax=0.1,
       speed=0.5,pl=2,output=F,graphfile=0)
}
\arguments{
  \item{etpdat}{Eccentricity, tilt, precession astronmical series. First column = time, second column = ETP. If not entered, will use default series.}
  \item{iter}{Number of simulations.}
  \item{rho}{AR(1) coefficient.}
  \item{ARvar}{Multiplicative factor for AR1 noise (1= eqivalent to ETP variance).}
  \item{sig}{Demarcate what confidence level (percent) on plots?}
  \item{tbw}{MTM time-bandwidth product.}
  \item{padfac}{Padding factor.}
  \item{ftest}{Include MTM harmonic f-test results? (T or F)}
  \item{fmax}{Maximum frequency for plotting.}
  \item{speed}{Set the amount of time to pause before plotting new graph, in seconds.}
  \item{pl}{Plot log power (1) or linear power (2)?}
  \item{output}{Output modeled time series? (T or F). If selected, only one simulation performed.}
  \item{graphfile}{Output a pdf or jpg image of each plot? 0 = no, 1 = pdf, 2 = jpeg.  If yes, there will be no output to screen. Individual graphic files will be produced for each simluation, for assembling into a movie.}
}