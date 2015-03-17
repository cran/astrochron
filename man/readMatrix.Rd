\name{readMatrix}
\alias{readMatrix}
\title{Read Data Matrix from File}
\description{
Read data matrix from a file, either tab-delimited, CSV, or semicolon-delimited.
}
\usage{
readMatrix(d=1,h="auto",output=1,genplot=F)
}
\arguments{
  \item{d}{What column delimiter is used? (0 = tab/.txt, 1 = comma/.csv, 2 = semicolon). CSV is the default option, which interfaces well with EXCEL.}
  \item{h}{Does the data file have column titles/headers? ("yes", "no", "auto"). "auto" will auto detect column titles/headers, which must be single strings and start with a character.}
  \item{output}{Return data as: 1= matrix, 2=data frame}
  \item{genplot}{generate summary plots (T or F).}
}
\details{
Missing values (in the file that you are reading from) should be indicated by 'NA'.
If you have included characters in the column titles that are not permitted by R, they will be modified!
}