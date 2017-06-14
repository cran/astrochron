\name{prewhiteAR}
\alias{prewhiteAR}
\title{Prewhiten stratigraphic series with autoregressive filter, order selected by Akaike Information Criterion}
\description{
Prewhiten stratigraphic series using autoregressive (AR) filter. Appropriate AR order can be automatically determined using the Akaike Information Criterion, or alternatively, the order may be predefined.
}
\usage{
prewhiteAR(dat,order=0,method="mle",aic=T,genplot=T,verbose=T)
}
\arguments{
  \item{dat}{Stratigraphic series for prewhitening. First column should be location (e.g., depth), second column should be data value for prewhitening. Series must have uniform sampling interval.}
  \item{order}{AR order for prewhitening (if aic=F), or alternatively, the maximum AR order to investigate (if aic=T). If order is set to <=0, will evaluate up to maximum default order (this varies based on method).}
  \item{method}{Method for AR parameter estimation: ("yule-walker", "burg", "ols", "mle", "yw")}
  \item{aic}{Select model using AIC? if F, will use order. AIC is only strictly valid if method is "mle".}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}

\references{
Akaike, H. (1974), \emph{A new look at the statistical model identification}, IEEE Trans. Autom. Control, 19, 716-723, doi:10.1109/TAC.1974.1100705.
}

\seealso{ \code{\link{ar}}, \code{\link{arcsinT}}, \code{\link{bandpass}}, \code{\link{demean}}, \code{\link{detrend}}, \code{\link{divTrend}}, \code{\link{logT}}, \code{\link{lowpass}}, \code{\link{noKernel}}, and \code{\link{prewhiteAR1}} }