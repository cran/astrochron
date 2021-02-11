\name{pl}
\alias{pl}
\title{Set up plots}
\description{
Open new device and set up for multiple plots, output to screen or PDF if desired.
}
\usage{
pl(n,r,c,h,w,mar,file,title)
}
\arguments{
  \item{n}{Number of plots per page (1-25). When specified, this parameter takes precedence and the default settings for r and c are used (the r and c options below are ignored).}
  \item{r}{Number of rows of plots.}
  \item{c}{Number of columns of plots.}
  \item{h}{Height of new page (a.k.a. "device").}
  \item{w}{Width of new page (a.k.a. "device").}
  \item{mar}{A numerical vector of the form c(bottom, left, top, right) which gives the margin size specified in inches.}
  \item{file}{File name, in quotes. Accepted file formats include .pdf, .jpg, .png, .tiff, .bmp; the format must be indicated using the appropriate filename extension at the end of the file name. If a file name is not designated, the plot is output to the screen instead.}
  \item{title}{Plot title (must be in quotes)}
}
