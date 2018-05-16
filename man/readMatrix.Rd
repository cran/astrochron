\name{readMatrix}
\alias{readMatrix}
\title{Read data matrix from file}
\description{
Read data matrix from a file, either tab-delimited, CSV, or semicolon-delimited.
}
\usage{
readMatrix(file=NULL,d=1,h="auto",skip=0,output=1,check=T,genplot=F)
}
\arguments{
  \item{file}{An optional file name, which must be in quotes (use the full directory path if the file is not in your present working directory). When a file name is not specified (the default), the file will be selected using a graphical user interface.}
  \item{d}{What column delimiter is used? (0 = tab/.txt, 1 = comma/.csv, 2 = semicolon). CSV is the default option, which interfaces well with EXCEL.}
  \item{h}{Does the data file have column titles/headers? ("yes", "no", "auto"). "auto" will auto detect column titles/headers, which must be single strings and start with a character.}
  \item{skip}{Number of lines to skip before beginning to read file}
  \item{output}{Return data as: 1= matrix, 2=data frame}
  \item{check}{Check for empty entries in the matrix? (T or F).}
  \item{genplot}{generate summary plots (T or F).}
}
\details{
Missing values (in the file that you are reading from) should be indicated by 'NA'.
If you have included characters in the column titles that are not permitted by R, they will be modified!
}