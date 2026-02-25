\name{timeOptBSim}
\alias{timeOptBSim}
\title{TimeOptBSim: Significance testing via Monte Carlo simulation for TimeOptB}
\description{
TimeOptBSim: Perform Monte Carlo AR(2) simulations to evaluate significance of TimeOptB results, as in Malinverno & Meyers (2024).
}
\usage{
timeOptBSim(res,Nsim=1000,genplot=TRUE,check=TRUE,verbose=TRUE)           
}
\arguments{
  \item{res}{Vector returned from timeOptB, which contains the parameters for the Monte Carlo simulation.}
  \item{Nsim}{Number of Monte Carlo simulations.}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{check}{Conduct compliance checks before processing? (T or F) In general this should be activated.}
  \item{verbose}{Verbose output? (T or F)}
}

\details{
The TimeOptB algorithm (Malinverno & Meyers, 2024) is a Bayesian version of TimeOpt (Meyers, 2015),
which evaluates stratigraphic data sets to assess (1) the concentration of spectral power at specified 
target astronomical periods (spectral fit), and (2) eccentricity amplitude modulations within the precession band 
(envelope fit).  The envelope fit can optionally be omitted.
  
For a given cyclostratigraphic data set, TimeOptB calculates the posterior probability density function (PDF) 
of the axial precession frequency (k) and of the sedimentation rate (u), using a 2D grid of u and k. The 
Solar System frequencies g_i and s_i are set to fixed values based on the age of the strata (following 
Table 2 in Hoang et al., 2021). The Bayesian prior for the precession frequency is age-specific, based 
on a polynomial fit to the model results of Farhat et al. (2022), with uncertainties from Waltham (2015).
The Bayesian prior for sedimentation rate is defined as a uniform distribution (between 'sedmin' and 'sedmax'). 

TimeOptBSim evaluates the significance of the TimeOptB results -- for the maximum a posteriori (MAP) values of k and u -- 
using Monte Carlo simulation with AR(2) stochastic series.  The reported p-values for 'Eccentricity only',
'Obliquity only', and 'Climatic precession only' should only be used when the results for 'All astronomical cycles'
are considered significant (here taken as p-value less than or equal to 0.1).


For additional guidance on the application of TimeOptB and TimeOptBSim, please see Malinverno & Meyers (2024).
}

\seealso{ \code{\link{timeOpt}}, \code{\link{timeOptB}}, \code{\link{timeOptMCMC}}, \code{\link{timeOptBMCMC}} }

\references{
Meyers, S.R., 2015,
\emph{The evaluation of eccentricity-related amplitude modulation and bundling in paleoclimate data: An inverse approach for astrochronologic testing and time scale optimization}: Paleoceanography, 30, doi:10.1002/2015PA002850.

Malinverno, A. and Meyers, S.R., 2024,
\emph{Bayesian estimation of past astronomical frequencies, lunar distance, and length of day from sediment cycles}: Geochemistry, Geophysics, Geosystems, 25, e2023GC011176.

Hoang, N.H., Mogavero, F., Laskar, J., 2021,
\emph{Chaotic diffusion of the fundamental frequencies in the Solar System}: Astronomy & Astrophysics, 654, A156. 

Waltham, 2015,
\emph{Milankovitch period uncertainties and their impact on cyclostratigraphy}: Journal of Sedimentary Research, 85, 990-998. 
}

\examples{
 \dontrun{
# Obtain the Xiamaling Cu/Al dataset from the Astrochron server
CuAl=getData("Xiamaling-CuAl")

# Isolate interval of interest and interpolate the data to the median sampling interval of 0.012 m.
CuAl_0.012=linterp(iso(CuAl,xmin=263.42,xmax=265.496,genplot=FALSE))

# Run timeOptB
res = timeOptB(CuAl_0.012, sedmin=0.3, sedmax=0.4, age=1400, kmin=78, kmax=94)

# Run timeOptBSim
timeOptBSim(res,Nsim=1000)
 }
}