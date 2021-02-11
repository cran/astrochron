\name{synthStrat}
\alias{synthStrat}
\title{Synthesize stratigraphy from forcing function}
\description{
Synthesize stratigraphy from forcing function.
}
\usage{
synthStrat(signal=NULL,nfacies=4,clip=T,flip=F,fmax=0.1,output=F,genplot=2,verbose=T)
}
\arguments{
  \item{signal}{Forcing signal. First column should be time (in ka), second column should be forcing.}
  \item{nfacies}{Number of sedimentary facies to model.}
  \item{clip}{Clip forcing signal at mean value? (T or F)}
  \item{flip}{Reverse the sign of the forcing? (T or F)}
  \item{fmax}{Maximum frequency for spectra (if genplot=2).}
  \item{output}{Output facies series? (T or F)}
  \item{genplot}{Generate summary plots? (1) plot stratigraphy, (2) plot statigraphy and spectra.}
  \item{verbose}{Verbose output? (T or F)}
}

\value{
modeled stratigraphic series.
}

\examples{
 \dontrun{
# EX.1: precession, unclipped
signal=etp(tmin=8400,tmax=8900,pWt=1,oWt=0,eWt=0)
synthStrat(signal,nfacies=4,clip=FALSE,genplot=2)

# EX.2: more finely resolved facies
#synthStrat(signal,nfacies=15,clip=FALSE,genplot=2)

# EX.3: couplets
#synthStrat(signal,nfacies=2,clip=FALSE,genplot=2)

# EX.4: precession, clipped
#synthStrat(signal,nfacies=4,genplot=2)

# EX.5: noise
noise=ar1(npts=501,rho=0.8)
#synthStrat(noise,nfacies=4,genplot=2)

# EX.6: precession + noise 
#signal2=signal
#signal2[2]=signal2[2]+0.75*noise[2]
#synthStrat(signal2,nfacies=4,genplot=2)

# EX.7: p-0.5t, clipped (demonstrates interference pattern; compare with EX.4
#signal3=etp(tmin=8400,tmax=8900,pWt=1,oWt=-0.5,eWt=0)
#synthStrat(signal3,nfacies=4,genplot=2)

# EX.8: ice sheet model, using p-0.5t
#ice=imbrie()
#synthStrat(ice,nfacies=5,clip=FALSE,genplot=2)

# EX.9: precession, clipped, ramping sedimentation rate
#synthStrat(linterp(sedRamp(signal,genplot=FALSE),genplot=FALSE),nfacies=6,
# clip=TRUE,genplot=2,fmax=10)

 }
}