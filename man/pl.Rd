\name{pl}
\alias{pl}
\title{Set Up Plots}
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
  \item{file}{PDF file name, in quotes. If a file name is not designated, then the plot is output to the screen instead.}
  \item{title}{Plot title (must be in quotes)}
}
