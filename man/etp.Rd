\name{etp}
\alias{etp}
\title{Generate Eccentricity-Tilt-Precession Models}
\description{
Calculate eccentricity-tilt-precession time series using solutions from Laskar et al. (2004, 2011)
}
\usage{
etp(laskar,tmin=0,tmax=1000,dt=1,eWt=1,oWt=1,pWt=1,sol=5,esinw=F,standardize=T,genplot=T,
    verbose=T)
}
\arguments{
  \item{laskar}{Laskar et al. (2004, 2011) astronomical solutions.}
  \item{tmin}{Start time (ka before present, J2000) for ETP.}
  \item{tmax}{End time (ka before present, J2000) for ETP.}
  \item{dt}{Sample interval for ETP (ka). Minimum = 1 ka.}
  \item{eWt}{Relative weight applied to eccentricity solution.}
  \item{oWt}{Relative weight applied to obliquity solution.}
  \item{pWt}{Relative weight applied to precession solution.}
  \item{sol}{Eccentricity solution to use: 1= LA10a, 2= LA10b, 3= LA10c, 4= LA10d, 5= LA04.}
  \item{esinw}{Use esinw in ETP calculation? (T or F).}
  \item{standardize}{Standardize (subtract mean, divide by standard deviation) precession, obliquity and eccentricity series before applying weight and combining? (T or F)}
  \item{genplot}{Generate summary plots? (T or F).}
  \item{verbose}{Verbose output? (T or F).}
}

\value{
Eccentricity + tilt + precession.
}

\references{ 
Laskar, J., Robutel, P., Joutel, F., Gastineau, M., Correia, A.C.M., Levrard, B., 2004, 
\emph{A long term numerical solution for the insolation quantities of the Earth}: Astron. Astrophys., Volume 428, 261-285.

Laskar, J., Fienga, A., Gastineau, M., Manche, H., 2011, 
\emph{La2010: A new orbital solution for the long-term motion of the Earth}: Astron. Astrophys., Volume 532, A89.
}