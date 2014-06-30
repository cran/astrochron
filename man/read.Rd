\name{read}
\alias{read}
\title{Read Data from File}
\description{
Read stratigraphic data series from a file, either tab-delimited or CSV. First column must contain location data (depth, height, time). The function will remove missing entries, sort by location, average duplicate values, and generate summary plots.
}
\usage{
read(d=1,h="auto",srt=T,ave=T,genplot=T)
}
\arguments{
  \item{d}{What column delimiter is used? (0 = tab/.txt, 1 = comma/.csv). CSV is the default option, which interfaces well with EXCEL.}
  \item{h}{Does the data file have column titles/headers? ("yes", "no", "auto"). "auto" will auto detect column titles/headers, which must be single strings and start with a character.}
  \item{srt}{Sort data values by first column? (T or F)}
  \item{ave}{Average duplicate values? (T or F). Only applies if input file has 2 columns}
  \item{genplot}{generate summary plots (T or F).}
}
\details{
Missing values (in the file that you are reading from) should be indicated by 'NA'.
If you have included characters in the column titles that are not permitted by R, they will be modified!
}