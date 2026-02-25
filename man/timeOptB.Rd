\name{timeOptB}
\alias{timeOptB}
\title{TimeOptB: A Bayesian implementation of TimeOpt that calculates posterior PDFs for precession frequency and sedimentation rate}
\description{
TimeOptB: A Bayesian implementation of TimeOpt that calculates posterior PDFs for precession frequency and sedimentation rate. This function follows the approach introduced in Malinverno and Meyers (2024), with some improvements.
}
\usage{
timeOptB(dat,age=NULL,env=TRUE,mvopt=0,sedmin=0.5,sedmax=5,numsed=120,kmin=NULL,kmax=NULL,
         numk=120,ftol=0.005,roll=10^7,detrend=TRUE,palette=6,genplot=TRUE,
         check=TRUE,verbose=TRUE)           
}
\arguments{
  \item{dat}{Stratigraphic series for astrochronologic assessment. First column should be depth or height (in meters), second column should be data value.}
  \item{age}{Age of stratagraphic series, in Ma.}
  \item{env}{Include envelope fit in the posterior estimates? (T or F)}
  \item{mvopt}{This option specifies the s-terms to include (for obliquity). (0) Include s1 through s6; (1) Include s3 and s4; (2) No s-terms included.}
  \item{sedmin}{Minimum sedimentation rate for investigation (cm/ka).}
  \item{sedmax}{Maximum sedimentation rate for investigation (cm/ka).}
  \item{numsed}{Number of sedimentation rates to investigate in grid.}
  \item{kmin}{Minimum precession frequency for investigation (arcsec/yr).}
  \item{kmax}{Maximum precession frequency for investigation (arcsec/yr).}
  \item{numk}{Number of precession frequencies to investigate in grid.}
  \item{ftol}{Tolerance in cycles/ka used to define the precession bandpass. This value is added to the highest precession frequency, and subtracted from the lowest precession frequency, to define the half power points for the Taner bandpass filter.}
  \item{roll}{Taner filter roll-off rate, in dB/octave.}  
  \item{detrend}{Remove linear trend from stratigraphic series? (T or F)}
  \item{palette}{What color palette would you like to use? (1) rainbow, (2) grayscale, (3) blue, (4) red, (5) blue-white-red (if values are negative and positive, white is centered on zero), (6) viridis}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{check}{Conduct compliance checks before processing? (T or F) In general this should be activated.}
  \item{verbose}{Verbose output? (T or F)}
}

\details{
The TimeOptB algorithm (Malinverno & Meyers, 2024) is a Bayesian version of TimeOpt (Meyers, 2015),
which evaluates stratigraphic data sets to assess (1) the concentration of spectral power at specified 
target astronomical periods (spectral fit), and (2) eccentricity amplitude modulations within the precession band 
(envelope fit).  The envelope fit can optionally be omitted (when 'env' set to FALSE).
  
For a given cyclostratigraphic data set, TimeOptB calculates the posterior probability density function (PDF) 
of the axial precession frequency (k) and of the sedimentation rate (u), using a 2-dimensional grid of u and k. The 
Solar System frequencies g_i and s_i are set to fixed values based on the age of the strata (following 
Table 2 in Hoang et al., 2021). The Bayesian prior for the precession frequency is age-specific, based 
on a polynomial fit to the model results of Farhat et al. (2022), with uncertainties from Waltham (2015).
The Bayesian prior for sedimentation rate is defined as a uniform distribution (between 'sedmin' and 'sedmax'). 

For additional guidance on the application of TimeOptB, please see Malinverno & Meyers (2024).
}

\seealso{ \code{\link{timeOpt}}, \code{\link{timeOptBSim}}, \code{\link{timeOptMCMC}}, \code{\link{timeOptBMCMC}} }

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