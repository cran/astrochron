\name{diffAccum}
\alias{diffAccum}
\title{Model differential accumulation}
\description{
Model differential accumulation. The input variable (e.g., insolation, proxy value) is rescaled to sedimentation rate curve varying from sedmin to sedmax. Input series must be evenly sampled in time.
}
\usage{
diffAccum(dat,sedmin=0.01,sedmax=0.02,dir=1,genplot=T,verbose=T)
}
\arguments{
  \item{dat}{Model input series with two columns. First column must be time in ka, second column should be data value. Data series must be evenly sampled in time.}
  \item{sedmin}{Minimum sedimentation rate (m/ka)}
  \item{sedmax}{Maximum sedimentation rate (m/ka)}
  \item{dir}{1=peaks have higher accumulation rate, -1=troughs have higher accumulation rate}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}

\examples{
# generate model with one 20 ka cycle
ex <- cycles(1/20)

diffAccum(ex)
}