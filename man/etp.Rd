\name{etp}
\alias{etp}
\title{Generate eccentricity-tilt-precession models}
\description{
Calculate eccentricity-tilt-precession time series using the theoretical astronomical solutions.
By default, the Laskar et al. (2004) solutions will be downloaded. Alternatively, one can
specify the astronomical solution.
}
\usage{
etp(tmin=NULL,tmax=NULL,dt=1,eWt=1,oWt=1,pWt=1,esinw=T,solution=NULL,standardize=T,
     genplot=T,verbose=T)
}
\arguments{
  \item{tmin}{Start time (ka before present, J2000) for ETP. Default value is 0 ka, unless the data frame 'solution' is specified, in which case the first time datum is used.}
  \item{tmax}{End time (ka before present, J2000) for ETP. Default value is 1000 ka, unless the data frame 'solution' is specified, in which case the last time datum is used.}
  \item{dt}{Sample interval for ETP (ka). Minimum = 1 ka.}
  \item{eWt}{Relative weight applied to eccentricity solution.}
  \item{oWt}{Relative weight applied to obliquity solution.}
  \item{pWt}{Relative weight applied to precession solution.}
  \item{esinw}{Use e*sinw in ETP calculation? (T or F). If set to false, sinw is used.}
  \item{solution}{A data frame containing the astronomical solution to use. The data frame must have four columns: Time (ka, positive and increasing), Precession Angle, Obliquity, Eccentricity.}
  \item{standardize}{Standardize (subtract mean, divide by standard deviation) precession, obliquity and eccentricity series before applying weight and combining? (T or F)}
  \item{genplot}{Generate summary plots? (T or F).}
  \item{verbose}{Verbose output? (T or F).}
}

\details{
Note: If you plan to repeatedly execute the etp function, it is advisable to download the 
astronomical solution once using the function getLaskar.

Note: It is common practice to construct ETP models that have specified variance ratios
(e.g., 1:1:1 or 1:0.5:0.5) for eccentricity, obliquity and precession. In order to 
construct such models, it is necessary to choose 'standardize=T', and to set the individual 
weights (eWt, oWt, pWt) to the square root of the desired variance contribution.
}

\value{
Eccentricity + tilt + precession.
}

\seealso{ \code{\link{getLaskar}}}

\references{ 
Laskar, J., Robutel, P., Joutel, F., Gastineau, M., Correia, A.C.M., Levrard, B., 2004, 
\emph{A long term numerical solution for the insolation quantities of the Earth}: Astron. Astrophys., Volume 428, 261-285.

Laskar, J., Fienga, A., Gastineau, M., Manche, H., 2011, 
\emph{La2010: A new orbital solution for the long-term motion of the Earth}: Astron. Astrophys., Volume 532, A89.

Laskar, J., Gastineau, M., Delisle, J.-B., Farres, A., Fienga, A.: 2011,
\emph{Strong chaos induced by close encounters with Ceres and Vesta}: Astron. Astrophys., Volume 532, L4.
}

\examples{
# create an ETP model from 10000 ka to 20000 ka, with a 5 ka sampling interval
# this will automatically download the astronomical solution
ex=etp(tmin=10000,tmax=20000,dt=5)

# alternatively, download the astronomical solution first
ex2=getLaskar()
ex=etp(tmin=10000,tmax=20000,dt=5,solution=ex2)
}