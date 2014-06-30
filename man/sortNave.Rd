\name{sortNave}
\alias{sortNave}
\title{Remove Missing Entries, Sort Data, Average Duplicates}
\description{
Sort and average duplicates in stratigraphic series, as performed in 'read' function.
}
\usage{
sortNave(dat,sortDecr=F,ave=T,genplot=T,verbose=T)
}
\arguments{
  \item{dat}{Stratigraphic series for processing. First column should be location (e.g., depth), second column should be data value.}
  \item{sortDecr}{Sorting direction? (F=increasing, T=decreasing)}
  \item{ave}{Average duplicate values? (T or F)}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}