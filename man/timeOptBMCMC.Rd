\name{timeOptBMCMC}
\alias{timeOptBMCMC}
\title{TimeOptBMCMC: A Bayesian implementation of TimeOpt that calculates posterior PDFs for precession frequency, sedimentation rate, and Solar System g- and s-frequencies, using adaptive MCMC}
\description{
TimeOptBMCMC: A Bayesian implementation of TimeOpt that calculates posterior PDFs for precession 
frequency, sedimentation rate, and solar system g- and s-frequencies, using adaptive Markov chain Monte Carlo.
This function follows the approach introduced in Malinverno and Meyers (2024), with some improvements, 
and is an update of the method introduced in Meyers and Malinverno (2018).
}
\usage{
timeOptBMCMC(dat,age=NULL,env=TRUE,mvopt=0,sedmin=0.5,sedmax=5,ftol=0.005,roll=10^7,
             nsamples=50000,test=FALSE,detrend=TRUE,savefile=FALSE,genplot=TRUE,
             check=TRUE,verbose=TRUE)
}
\arguments{
  \item{dat}{Stratigraphic series for astrochronologic assessment. First column should be depth or height (in meters), second column should be data value.}
  \item{age}{Age of stratagraphic series, in Ma.}
  \item{env}{Include envelope fit in the posterior estimates? (T or F)}
  \item{mvopt}{This option specifies the s-terms to include (for obliquity). (0) Include s1 through s6; (1) Include s3 and s4; (2) No s-terms included.}
  \item{sedmin}{Minimum sedimentation rate for investigation (cm/ka).}
  \item{sedmax}{Maximum sedimentation rate for investigation (cm/ka).}
  \item{ftol}{Tolerance in cycles/ka used to define the precession bandpass. This value is added to the highest precession frequency, and subtracted from the lowest precession frequency, to define the half power points for the Taner bandpass filter.}
  \item{roll}{Taner filter roll-off rate, in dB/octave.}  
  \item{nsamples}{Number of MCMC samples (length of MCMC chain).}
  \item{test}{Activate epsilon testing mode? This option will assign all MCMC samples a log-likelihood of unity. This provides a diagnostic check to ensure that the applied epsilon values are sampling the entire range of parameter values. (T or F)}
  \item{detrend}{Remove linear trend from stratigraphic series? (T or F)}
  \item{savefile}{Save MCMC samples to file MCMCsamples.csv? (T or F)}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{check}{Conduct compliance checks before processing? (T or F) In general this should be activated.}
  \item{verbose}{Verbose output? (T or F)}
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

For additional guidance on the application of TimeOptBMCMC, please see Malinverno & Meyers (2024).
}

\seealso{ \code{\link{timeOpt}}, \code{\link{timeOptB}}, \code{\link{timeOptMCMC}}, \code{\link{timeOptBSim}} }

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