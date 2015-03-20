\name{ar1etp}
\alias{ar1etp}
\title{AR(1) + ETP simulation Routine}
\description{
Simulate a combined AR(1) + ETP signal, plot spectrum and confidence levels
}
\usage{
ar1etp(etpdat=NULL,nsim=100,rho=0.9,wtAR=1,sig=90,tbw=2,padfac=5,ftest=F,fmax=0.1,
       speed=0.5,pl=2,graphfile=0)
}
\arguments{
  \item{etpdat}{Eccentricity, tilt, precession astronmical series. First column = time, second column = ETP. If not entered will use default series from Laskar et al. (2004), spanning 0-1000 kyr.}
  \item{nsim}{Number of simulations.}
  \item{rho}{AR(1) coefficient for noise modeling.}
  \item{wtAR}{Multiplicative factor for AR1 noise (1= eqivalent to ETP variance).}
  \item{sig}{Demarcate what confidence level (percent) on plots?}
  \item{tbw}{MTM time-bandwidth product.}
  \item{padfac}{Pad with zeros to (padfac*npts) points, where npts is the number of data points.}
  \item{ftest}{Include MTM harmonic f-test results? (T or F)}
  \item{fmax}{Maximum frequency for plotting.}
  \item{speed}{Set the amount of time to pause before plotting new graph, in seconds.}
  \item{pl}{Plot log power (1) or linear power (2)?}
  \item{graphfile}{Output a pdf or jpg image of each plot? 0 = no, 1 = pdf, 2 = jpeg.  If yes, there will be no output to screen. Individual graphic files will be produced for each simluation, for assembling into a movie.}
}

\details{
Note: Setting wtAR=1 will provide equal variance contributions from the etp model and the ar1 model.
More generally, set wtAR to the square root of the desired variance contribution (wtAR=0.5 will generate an AR1 model with variance that is 25\% of the etp model).

Note: You may use the function etp to generate eccentricity-tilt-precession models.
}

\seealso{ \code{\link{getLaskar}}, and \code{\link{etp}}}

\references{ 
Laskar, J., Robutel, P., Joutel, F., Gastineau, M., Correia, A.C.M., Levrard, B., 2004, 
\emph{A long term numerical solution for the insolation quantities of the Earth}: Astron. Astrophys., Volume 428, 261-285.
}

\examples{
\donttest{

# run simulations using the default settings
ar1etp()

# compare with a second model:
# generate etp model spanning 0-2000 ka, with sampling interval of 5 ka.
ex1=etp(tmin=0,tmax=2000,dt=5)
# run simulations, with rho=-.7, and scaling noise to have 50% of the etp model variance
ar1etp(etpdat=ex1,rho=0.7,wtAR=sqrt(0.5))

}
}