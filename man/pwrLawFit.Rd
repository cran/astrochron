\name{pwrLawFit}
\alias{pwrLawFit}
\title{Estimate power law (1/f) fit to power spectrum}
\description{
Estimate power law (1/f) fit to power spectrum, following the algorithm of Vaughan (2005).
}
\usage{
pwrLawFit(spec,dof=2,flow=NULL,fhigh=NULL,output=1,genplot=T,verbose=T)

}
\arguments{
  \item{spec}{Power spectrum. First column is frequency, second column is raw power (linear). Do not include the zero frequency and Nyquist.}
  \item{dof}{Degrees of freedom for power spectral estimate. Default is 2, for a simple periodogram.}
  \item{flow}{Lowest frequency to include in 1/f fit}
  \item{fhigh}{Highest frequency to include in 1/f fit}
  \item{output}{Output results of 1/f fit? (0=none; 1=Frequency,Power,Power Law CL,Unbiased Power Law fit,CL_90,CL_95,CL_99; 2=beta, unbiased log10N, biased log10N)}
  \item{genplot}{generate summary plots (T or F)}
  \item{verbose}{verbose output (T or F)}
}

\references{
Vaughan, S. (2005), \emph{A simple test for periodic signals in red noise}, Astronomy & Astrophysics.
}

\examples{
# generate example series with periods of 400 ka, 100 ka, 40 ka and 20 ka
ex = cycles(freqs=c(1/400,1/100,1/40,1/20),start=1,end=1000,dt=5)

# add AR1 noise
noise = ar1(npts=200,dt=5,sd=.5)
ex[2] = ex[2] + noise[2]

# calculate periodogram
res=periodogram(ex,output=1,padfac=1)

# extract power and remove the Nyquist frequency
resPwr=cb(res,c(1,3))
resPwr=resPwr[-length(resPwr[,1]),]

pwrLawFit(resPwr)
}