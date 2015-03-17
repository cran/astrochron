\name{testPrecession}
\alias{testPrecession}
\title{Astrochronologic testing via the amplitude modulation approach of Zeeden et al. (2015).}
\description{
Astrochronologic testing via the amplitude modulation approach of Zeeden et al. (2015).
}
\usage{
testPrecession(dat,nsim=1000,gen=1,rho=NULL,esinw=NULL,output=T,genplot=T,verbose=T)
}
\arguments{
  \item{dat}{Stratigraphic series to analyze. First column should be location (time in kyr, a positive value), second column should be data value.}
  \item{nsim}{Number of Monte Carlo simulations (phase-randomized surrogates or AR1 surrogates).}
  \item{gen}{Monte Carlo simulation generator: (1) use phase-randomized surrogates, (2) use AR1 surrogates.}
  \item{rho}{Specified lag-1 correlation coefficient (rho). This value is only used if gen=2. If rho is not specified, it will be calculated within the function.}  
  \item{esinw}{Theoretical target 'eccentricity * sin(omega)' used for astrochronologic testing.  By default this is automatically determined within the function, using the solution of Laskar et al. (2004).}  
  \item{output}{Return results as a new data frame? (T or F)}  
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}

\details{
This astrochronologic testing method compares observed precession-scale amplitude modulations to those expected from the theoretical eccentricity solutions. It is applicable for testing astrochronologies spanning 0-50 Ma. 
The technique implements a series of filters to guard against artificial introduction of eccentricity modulations during tuning and data processing, and evaluates the statistical significance of the results using Monte Carlo simulation (Zeeden et al., 2015). 

The astronomically-tuned data series under evaluation should consist of two columns: time in kiloyears & data value. 
Note that time must be positive. The default astronomical solutions used for the astrochronologic testing come from Laskar et al. (2004).

When reporting a p-value for your result, it is important to consider the number of simulations used. 
A factor of 10 is appropriate, such that for 1000 simulations one would report a minimum p-value of "p<0.01", and for 10000 simulations one would report a minimum p-value of "p<0.001".

Please be aware that the kernel density estimate plots, which summarize the simulations, represent 'smoothed' models. Due to the smoothing bandwidth, they can sometimes give the impression of simulation values that are larger or smaller than actually present.  However, the reported p-value does not suffer from these issues.
}


\value{
When nsim is set to zero, the function will output a data frame with five columns:
 
1=time, 2=precession bandpass filter output, 3=amplitude envelope of (2), 4=lowpass filter output of (3), 5=theoretical eccentricity (as extracted from precession modulations using the filtering algorithm)

When nsim is > 0, the function will output the correlation coefficients for each simulation.
}


\references{
C. Zeeden, S.R. Meyers, L.J. Lourens, and F.J. Hilgen, 2015 (accepted),
\emph{Testing astronomically tuned age models}: Paleoceanography.

J. Laskar, P. Robutel, F. Joutel, M. Gastineau, A.C.M. Correia, and B. Levrard, B., 2004, 
\emph{A long term numerical solution for the insolation quantities of the Earth}: 
Astron. Astrophys., Volume 428, 261-285.
}

\examples{
\donttest{
### as a test series, use the three dominant precession terms from Berger et al. (1992)
ex<-cycles(start=0,end=1000,dt=2)

### now conduct astrochronologic testing
res1=testPrecession(ex)


### if you plan to run testPrecession repeatedly, it is advisable to download the astronomical
### solution and construct esinw first
ex2<-getLaskar()
ex3<-etp(tmin=0,tmax=1000,dt=2,eWt=0,oWt=0,pWt=1,esinw=TRUE,solution=ex2,standardize=FALSE)

### now conduct astrochronologic testing
res2<-testPrecession(ex,esinw=ex3)
}
}