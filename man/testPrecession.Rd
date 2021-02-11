\name{testPrecession}
\alias{testPrecession}
\title{Astrochronologic testing via the precession amplitude modulation approach of Zeeden et al. (2015).}
\description{
Astrochronologic testing via the precession amplitude modulation approach of Zeeden et al. (2015), as updated in Zeeden et al. (2018 submitted).
}
\usage{
testPrecession(dat,nsim=1000,gen=1,edge=0.025,maxNoise=1,rho=NULL,detrendEnv=T,
               solution=NULL,output=F,genplot=T,verbose=T)
}
\arguments{
  \item{dat}{Stratigraphic series to analyze. First column should be location (time in ka, a positive value), second column should be data value.}
  \item{nsim}{Number of Monte Carlo simulations (phase-randomized surrogates or AR1 surrogates).}
  \item{gen}{Monte Carlo simulation generator: (1) use phase-randomized surrogates, (2) use AR1 surrogates.}
  \item{edge}{Percentage of record to exclude from beginning and end of data series, to remove edge effects. (a value from 0-1)}  
  \item{maxNoise}{Maximum noise level to add in simulations. A value of 1 will apply maximum noise that is equivalent to 1 standard deviation of the data.}
  \item{rho}{Specified lag-1 correlation coefficient (rho). If rho is not specified, it will be calculated within the function.}  
  \item{detrendEnv}{Linearly detrend envelope? (T or F)}  
  \item{solution}{Theoretical solution used for astrochronologic testing. Solution should be in the format: time (ka), precession angle, obliquity, eccentricity (the output from function 'getLaskar'). By default this is automatically determined within the function, using the solution of Laskar et al. (2004).}  
  \item{output}{Return results as a new data frame? (T or F)}  
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}

\details{
This astrochronologic testing method compares observed precession-scale amplitude modulations to those expected from the theoretical eccentricity solutions. It is applicable for testing astrochronologies spanning 0-50 Ma. 
The technique implements a series of filters to guard against artificial introduction of eccentricity modulations during tuning and data processing, and evaluates the statistical significance of the results using Monte Carlo simulation (Zeeden et al., 2015). 

The algorithm includes an improvement in the significance testing approach. Specifically, as a safeguard against artificially imposed modulations, an adaptive noise addition step is implemented (as outlined in Zeeden et al., submitted).

The astronomically-tuned data series under evaluation should consist of two columns: time in kiloyears & data value. 
Note that time must be positive. The default astronomical solutions used for the astrochronologic testing come from Laskar et al. (2004).

When reporting a p-value for your result, it is important to consider the number of simulations used. 
A factor of 10 is appropriate, such that for 1000 simulations one would report a minimum p-value of "p<0.01", and for 10000 simulations one would report a minimum p-value of "p<0.001".

Please be aware that the kernel density estimate plots, which summarize the simulations, represent 'smoothed' models. Due to the smoothing bandwidth, they can sometimes give the impression of simulation values that are larger or smaller than actually present.  However, the reported p-value does not suffer from these issues.

IMPORTANT CHANGES (June 20, 2018): Note that this version has been updated to use 'solution' instead of 'esinw', for consistency with the function 'testTilt'. If you are invoking the default option, you do not need to make any changes to your script.
Also note that the new option 'edge' has been added, which by default will truncate your data series by 5 percent (2.5 percent on each end of the record), to guard against edge effects that can be present in the amplitude envelope. Set edge to 0 to reconstruct the
original (now legacy) 'testPrecession' approach.
}


\value{
When nsim is set to zero, the function will output a data frame with five columns:
 
1=time, 2=precession bandpass filter output, 3=amplitude envelope of (2), 4=lowpass filter output of (3), 5=theoretical eccentricity (as extracted from precession modulations using the filtering algorithm), 6=(2) + noise, 7=amplitude envelope of (6), 8=lowpass filter output of (7)

When nsim is > 0, the function will output the correlation coefficients for each simulation.
}


\seealso{ \code{\link{asm}}, \code{\link{eAsmTrack}}, \code{\link{timeOpt}}, and \code{\link{timeOptSim}} }


\references{
C. Zeeden, S.R. Meyers, L.J. Lourens, and F.J. Hilgen, 2015,
\emph{Testing astronomically tuned age models}: Paleoceanography, 30, doi:10.1002/2014PA002762.

C. Zeeden, S.R. Meyers, F.J. Hilgen, L.J. Lourens, and J. Laskar, submitted,
\emph{Time scale evaluation and the quantification of obliquity forcing}: Quaternary Science Reviews.

J. Laskar, P. Robutel, F. Joutel, M. Gastineau, A.C.M. Correia, and B. Levrard, B., 2004, 
\emph{A long term numerical solution for the insolation quantities of the Earth}: 
Astron. Astrophys., Volume 428, 261-285.
}

\examples{
\dontrun{
### as a test series, use the three dominant precession terms from Berger et al. (1992)
ex<-cycles(start=0,end=1000,dt=2)

### now conduct astrochronologic testing
res1=testPrecession(ex)


### if you plan to run testPrecession repeatedly, it is advisable to download the astronomical
### solution first
solution<-getLaskar()

### now conduct astrochronologic testing
res2<-testPrecession(ex,solution=solution)
}
}