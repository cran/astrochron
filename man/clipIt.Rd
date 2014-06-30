\name{clipIt}
\alias{clipIt}
\title{Create non-linear response by clipping stratigraphic series}
\description{
Create non-linear response by clipping stratigraphic series below a threshold value. Alternatively, mute response below a threshold value using a contant divisor.  Both approaches will enhance power in modulator (e.g., eccentricity) and diminish power the carrier (e.g., precession).
}
\usage{
clipIt(dat,thresh=NULL,clipval=NULL,clipdiv=NULL,genplot=T,verbose=T)
}
\arguments{
  \item{dat}{Stratigraphic series. First column should be location (e.g., depth), second column should be data value.}
  \item{thresh}{Clip below what theshold value? By default will clip at mean value.}
  \item{clipval}{What number should be assigned to the clipped values? By default, the value of thresh is used.}
  \item{clipdiv}{Clip using what divisor? A typical value is 2. By default, clipdiv is unity.}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}
