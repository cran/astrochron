\name{timeOptBMCMCplot}
\alias{timeOptBMCMCplot}
\title{TimeOptBMCMCplot: Generate summary figures and statistics for TimeOptBMCMC results}
\description{
TTimeOptBMCMCplot: Generate summary figures and statistics for TimeOptBMCMC results.
}
\usage{
timeOptBMCMCplot(dat,mcmcres,pdfpara,nburnin=NULL,fc_dat=NULL)
}
\arguments{
  \item{dat}{Stratigraphic series for astrochronologic assessment. First column should be depth or height (in meters), second column should be data value.}
  \item{mcmcres}{Results from timeOptBMCMC}
  \item{pdfpara}{Parameter file from timeOptBMCMC}
  \item{nburnin}{Number or MCMC samples to discard as burn-in interval. If set to NULL, a burn-in detection algorithm is implemented, using the median value from second half of the MCMC chain as the threshold value.}
  \item{fc_dat}{Fourier Coefficients for dat (optional)}
}

\details{
The TimeOptBMCMC algorithm (Malinverno & Meyers, 2024) is a Bayesian version of TimeOpt (Meyers, 2015),
which evaluates stratigraphic data sets to assess (1) the concentration of spectral power at specified 
target astronomical periods (spectral fit), and (2) eccentricity amplitude modulations within the precession band 
(envelope fit).  The envelope fit can optionally be omitted ('env' set to FALSE).
  
For a given cyclostratigraphic data set, TimeOptBMCMC calculates the posterior probability density function (PDF) 
of the axial precession frequency (k), sedimentation rate (u), Solar System frequencies g_i, and Solar System frequencies s_i, 
using an adaptive Markov chain Monte Carlo approach. The Bayesian priors for Solar System frequencies g_i and s_i 
are age-specific, following Table 2 in Hoang et al. (2021). The Bayesian prior for the precession frequency is 
age-specific, based on a polynomial fit to the model results of Farhat et al. (2022), with uncertainties from Waltham (2015).
The Bayesian prior for sedimentation rate is defined as a uniform distribution (between 'sedmin' and 'sedmax'). 

Function timeOptBMCMCplot generates summary figures and stistics for a given analysis, and is automatically called from 
function timeOptBMCMC. It can also be called independently; when doing so, be sure to use the prepared data that is returned
from timeOptBMCMC (see examples below).

For additional guidance on the application of TimeOptBMCMC, please see Malinverno & Meyers (2024).
}

\references{
Meyers, S.R., 2015,
\emph{The evaluation of eccentricity-related amplitude modulation and bundling in paleoclimate data: An inverse approach for astrochronologic testing and time scale optimization}: Paleoceanography, 30, doi:10.1002/2015PA002850.

Meyers, S.R. and Malinverno, A, 2018,
\emph{Proterozoic Milankovitch cycles and the history of the solar system}: Proceedings of the National Academy of Sciences, www.pnas.org/cgi/doi/10.1073/pnas.1717689115.

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

# Run timeOptBMCMC and plot results
res = timeOptBMCMC(CuAl_0.012, sedmin=0.3, sedmax=0.4, age=1400, nsamples=50000)

# If you want to re-plot the results from timeOptBMCMC, note that the Cu/Al data
# was detrended (if selected) and standardized in timeOptBMCMC. Use res$dat, returned
# from timeOptBMCMC
timeOptBMCMCplot(res$dat,res$mcmcres,res$pdfpara)
 }
}