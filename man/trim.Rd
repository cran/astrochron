\name{trim}
\alias{trim}
\title{Remove Outliers from Stratigraphic Series}
\description{
Automatically remove outliers from stratigraphic series, using 'boxplot' algorithm.
}
\usage{
trim(dat,c=1.5,genplot=T,verbose=T)
}
\arguments{
  \item{dat}{Stratigraphic series for outlier removal. First column should be location (e.g., depth), second column should be data value.}
  \item{c}{'c' defines the 'coef' variable for boxplot.stats. For more information: ?boxplot.stats}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}

\seealso{ \code{\link{delPts}}, \code{\link{idPts}}, \code{\link{iso}} and \code{\link{trimAT}} }