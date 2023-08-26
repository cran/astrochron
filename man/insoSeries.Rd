\name{insoSeries}
\alias{insoSeries}
\title{Laskar et al. (2004) insolation calculations}
\description{
Laskar et al. (2004) insolation calculations.  This is a wrapper function to calculate 
insolation using the palinsol package, following astrochron syntax.}
\usage{
insoSeries(tmin=0,tmax=1000,dt=1,opt=1,long=90,lat=65,threshold=400,l1=0,l2=70,S0=1365,
         genplot=TRUE,verbose=TRUE) 
}
\arguments{
  \item{tmin}{Minimum time for evaluation, in kyr (negative for future, positive for past; valid from -21000 to 51000)}
  \item{tmax}{Maximum time for evaluation, in kyr (negative for future, positive for past; valid from -21000 to 51000)}
  \item{dt}{Sampling interval in kyr. Must be >= 1}
  \item{opt}{1=Calculate daily mean insolation for given day and latitude; 2=calculate caloric insolation for given latitude; 3=calculate integrated insolation for given latitude; 4=calculating integrated insolation for given latitude, with threshold}
  \item{long}{True solar longitude, in degrees, 0-360}
  \item{lat}{Latitude on the Earth for evaluation, in degrees, -90 to 90}
  \item{threshold}{Threshold insolation, in W/m^2}
  \item{l1}{Lower true solar longitude bound of the time-integral, in degrees, 0-360}
  \item{l2}{Upper true solar longitude bound of the time-integral, in degrees, 0-360}
  \item{S0}{Solar constant, W/m^2}    
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}

\details{
This wrapper function performs four different types of insolation calculations based on 
the Laskar et al. (2004) astronomical parameters, using the 'palinsol' package (Crucifix, 2016).

The following 'palinsol' functions are used:

opt 1: Insol - Calculate Mean Insolation for Given Day and Latitude

opt 2: calins - Calculate Caloric Insolation for Given Latitude (Milankovitch; Berger, 1978)

opt 3: Insol_l1l2 - Calculate Integrated Insolation for Given Latitude (Berger et al., 2010)

opt 4: thrins - Calculate Integrated Insolation for Given Latitude, with Threshold (Huybers and Tziperman, 2008)

la04 -  Calculate astronomical parameters (Laskar et al., 2004)



Please see those functions for further details.

If opt 1 is selected: specify 'long' and 'lat'

If opt 2 is selected: specify 'lat'

If opt 3 is selected: specify 'lat', 'l1' and 'l2'

If opt 4 is selected: specify 'lat' and 'threshold'



APPROXIMATE CONVENTIONAL DATES FOR A GIVEN MEAN LONGITUDE (long):

0: 21 March (equinox)

30: 21 April

60: 21 May

90: 21 June (solstice)

120: 21 July

150: 21 August

180: 21 September (equinox)

210: 21 October

240: 21 November

270: 21 December (solstice)

300: 21 January

330: 21 February
}


\references{

Berger, A., 1978, 
\emph{Long-term variations of caloric insolation resulting from the earth's orbital elements}:
Quaternary Research, 9, 139 - 167.

Berger, A., Loutre, M.F., Yin Q., 2010, 
\emph{Total irradiation during any time interval of the year using elliptic integrals}:
Quaternary Science Reviews, 29, 1968 - 1982, doi:10.1016/j.quascirev.2010.05.007

Crucifix, M., 2016,
\emph{palinsol: Insolation for Palaeoclimate Studies. R package version 0.93}: https://CRAN.R-project.org/package=palinsol

Huybers, P., Tziperman, E., 2008, 
\emph{Integrated summer insolation forcing and 40,000-year glacial cycles: The perspective from an ice-sheet/energy-balance model}:
Paleoceanography, 23.

J. Laskar, P. Robutel, F. Joutel, M. Gastineau, A.C.M. Correia, and B. Levrard, B., 2004, 
\emph{A long term numerical solution for the insolation quantities of the Earth}: 
Astron. Astrophys., Volume 428, 261-285.
}