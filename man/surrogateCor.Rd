\name{surrogateCor}
\alias{surrogateCor}
\title{Estimate correlation coefficient and significance for serially correlated data}
\description{
Estimate correlation coefficient and significance for serially correlated data. This 
algorithm permits the analysis of data sets with different sampling grids, as 
discussed in Baddouh et al. (2016). The sampling grid from the data set with fewer 
points (in the common interval) is used for resampling. Resampling is conducted 
using piecewise-linear interpolation. 

If either dat1 or dat2 have only one column, the resampling is skipped. 

The significance of the correlation is determined using the method of Ebisuzaki W. (1997).
}
\usage{
surrogateCor(dat1,dat2,firstDiff=F,cormethod=1,nsim=1000,output=2,genplot=T,verbose=T)
}
\arguments{
  \item{dat1}{Data series with one or two columns. If two columns, first should be location (e.g., depth), second column should be data value.}
  \item{dat2}{Data series with one or two columns. If two columns, first should be location (e.g., depth), second column should be data value.}
  \item{firstDiff}{Calculate correlation using first differences? (T or F)}
  \item{cormethod}{Method used for calculation of correlation coefficient (1=Pearson, 2=Spearman rank, 3=Kendall)}
  \item{nsim}{Number of phase-randomized surrogate series to generate. If nsim <=1, simulation is deactivated.}
  \item{output}{Return which of the following?: 1= correlation coefficients for each simulation;  2= correlation coefficient for data series; 3= data values used in correlation estimate (resampled)}  
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}

\details{
Paraphrased from Baddouh et al. (2016):
To provide a quantitative evaluation of the correlation between two data sets that do
not share a common sampling grid, we introduce a statistical approach that employs sample 
interpolation, and significance testing with phase-randomized surrogate data (Ebisuzaki, 1997). 
The sparser sampling grid is used to avoid over-interpolation. Correlation is evaluated 
using Pearson, Spearman Rank, or Kendall rank coefficients. The statistical significance of the 
resulting correlation coefficients are estimated via Monte Carlo simulations 
using phase-randomized surrogates; the surrogates are subject to the same interpolation 
process, and compensate for potential serial correlation of data (Ebisuzaki, 1997).  

The first-difference series of each variable can also evaluated, to assess correlation in 
the magnitude of change between sequential stratigraphic samples rather than absolute 
magnitude.  
}

\seealso{ \code{\link{surrogates}} }

\references{
M. Baddouh, S.R. Meyers, A.R. Carroll, B.L. Beard, C.M. Johnson , 2016,
\emph{Lacustrine 87-Sr/86-Sr as a tracer to reconstruct Milankovitch forcing of the Eocene hydrologic cycle}: Earth and Planetary Science Letters.

W. Ebisuzaki, 1997,
\emph{A Method to Estimate the Statistical Significance of a Correlation When the Data Are Serially Correlated}: Journal of Climate, v. 10, p. 2147-2153.
}

\examples{
# generate two stochastic AR1 series
ex1 <- ar1(npts=100,dt=5)
ex2 <- ar1(npts=100,dt=6)

# calculate pearson correlation coefficient and p-value 
surrogateCor(ex1,ex2)
}