\name{sortNave}
\alias{sortNave}
\title{Remove missing entries, sort data, average duplicates}
\description{
Sort and average duplicates in stratigraphic series, as performed in 'read' function.
}
\usage{
sortNave(dat,sortDecr=F,ave=T,xmin=NULL,xmax=NULL,genplot=1,verbose=T)
}
\arguments{
  \item{dat}{Stratigraphic series for processing. First column should be location (e.g., depth), second column should be data value.}
  \item{sortDecr}{Sorting direction? (F=increasing, T=decreasing)}
  \item{ave}{Average duplicate values? (T or F)}
  \item{xmin}{Minimum x-axis value for plotting}
  \item{xmax}{Maximum x-axis value for plotting}
  \item{genplot}{Generate summary plots? 0=none, 1=stratigraphic series, distribution, box plot, Q-Q, 2=stratigraphic series}
  \item{verbose}{Verbose output? (T or F)}
}