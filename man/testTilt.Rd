\name{testTilt}
\alias{testTilt}
\title{Astrochronologic testing via the obliquity amplitude modulation approach of Zeeden et al. (2019 submitted).}
\description{
Astrochronologic testing via the obliquity amplitude modulation approach of Zeeden et al. (2019 submitted).
}
\usage{
testTilt(dat,nsim=1000,gen=1,edge=0.025,cutoff=1/150,maxNoise=0.25,rho=NULL,detrendEnv=T,
               solution=NULL,output=F,genplot=T,verbose=T)
}
\arguments{
  \item{dat}{Stratigraphic series to analyze. First column should be location (time in ka, a positive value), second column should be data value.}
  \item{nsim}{Number of Monte Carlo simulations (phase-randomized surrogates or AR1 surrogates).}
  \item{gen}{Monte Carlo simulation generator: (1) use phase-randomized surrogates, (2) use AR1 surrogates.}
  \item{edge}{Percentage of record to exclude from beginning and end of data series, to remove edge effects. (0-1)}  
  \item{cutoff}{Cutoff frequency for lowpass filtering.}  
  \item{maxNoise}{Maximum noise level to add in simulations. A value of 1 will apply maximum noise that is equivalent to 1 sd of data.}
  \item{rho}{Specified lag-1 correlation coefficient (rho). This value is only used if gen=2. If rho is not specified, it will be calculated within the function.}  
  \item{detrendEnv}{Linearly detrend envelope? (T or F)}  
  \item{solution}{Theoretical solution used for astrochronologic testing. Solution should be in the format: time (ka), precession angle, obliquity, eccentricity (the output from function 'getLaskar'). By default this is automatically determined within the function, using the solution of Laskar et al. (2004).}  
  \item{output}{Return results as a new data frame? (T or F)}  
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}

\details{
This astrochronologic testing method compares observed obliquity-scale amplitude modulations to those expected from the theoretical solutions. It is applicable for testing astrochronologies spanning 0-50 Ma. 
The technique implements a series of filters to guard against artificial introduction of modulations during tuning and data processing, and evaluates the statistical significance of the results using Monte Carlo simulation.
The algorithm includes an adaptive noise addition step to improvement the significance testing approach. See Zeeden et al. (2019 submitted) for additional information. 

The astronomically-tuned data series under evaluation should consist of two columns: time in kiloyears & data value. 
Note that time must be positive. The default obliquity solution used for the astrochronologic testing comes from Laskar et al. (2004).

When reporting a p-value for your result, it is important to consider the number of simulations used. 
A factor of 10 is appropriate, such that for 1000 simulations one would report a minimum p-value of "p<0.01", and for 10000 simulations one would report a minimum p-value of "p<0.001".

Please be aware that the kernel density estimate plots, which summarize the simulations, represent 'smoothed' models. Due to the smoothing bandwidth, they can sometimes give the impression of simulation values that are larger or smaller than actually present.  However, the reported p-value does not suffer from these issues.
}


\value{
When nsim is set to zero, the function will output a data frame with five columns:
 
1=time, 2=obliquity bandpass filter output, 3=amplitude envelope of (2), 4=lowpass filter output of (3), 5=theoretical obliquity (as extracted from modulations using the filtering algorithm), 6=(2) + noise, 7=amplitude envelope of (6), 8=lowpass filter output of (7)

When nsim is > 0, the function will output the correlation coefficients for each simulation.
}


\seealso{ \code{\link{asm}}, \code{\link{eAsmTrack}}, \code{\link{timeOpt}}, and \code{\link{timeOptSim}} }


\references{
C. Zeeden, S.R. Meyers, F.J. Hilgen, L.J. Lourens, and J. Laskar, 2019 submitted,
\emph{Time scale evaluation and the quantification of obliquity forcing}: Quaternary Science Reviews.

C. Zeeden, S.R. Meyers, L.J. Lourens, and F.J. Hilgen, 2015,
\emph{Testing astronomically tuned age models}: Paleoceanography, 30, doi:10.1002/2014PA002762.

J. Laskar, P. Robutel, F. Joutel, M. Gastineau, A.C.M. Correia, and B. Levrard, B., 2004, 
\emph{A long term numerical solution for the insolation quantities of the Earth}: 
Astron. Astrophys., Volume 428, 261-285.
}

\examples{
 \dontrun{
### as a test series, use the obliquity series from Laskar et al. (2004), spanning 
### the past 4 million years
ex<-etp(tmin=0,tmax=4000,dt=2,eWt=0,oWt=1,pWt=0,solution=solution,standardize=FALSE)

### now conduct astrochronologic testing
res1=testTilt(ex)

### if you plan to run testTilt repeatedly, it is advisable to download the astronomical
### solution
solution<-getLaskar()

### now conduct astrochronologic testing
res<-testTilt(ex,solution=solution)
 }
}