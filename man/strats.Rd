\name{strats}
\alias{strats}
\title{Summary Statistics for Stratigraphic Series}
\description{
Summary statistics for stratigraphic series: sampling interval and proxy values.
}
\usage{
strats(dat,output=0,genplot=1)
}
\arguments{
  \item{dat}{Stratigraphic series to evaluate. The first column should contain location (e.g., depth), and the second column should contain data value. This function also accepts non-stratigraphic (single column) input, in which case the sampling interval assessment is skipped.}
  \item{output}{Output: (0) nothing, (1) cumulative dt as percent of data points, (2) cumulative dt as percent of total interval duration, (3) dt by location}
  \item{genplot}{Generate summary plots? (0) none, (1) include plot of cumulative dt, (2) include dt histogram/density plot}
}

\details{
This function will generate a range of summary statistics for time series, including 
sampling interval information and the statistical distribution of proxy values.
}