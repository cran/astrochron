\name{xplot}
\alias{xplot}
\title{Generate cross-plot with kernel density estimates on axes}
\description{
Generate a cross-plot with kernel density estimates on axes. If multiple data points are superposed
in cross-plot, transparency of points reflects data density. Custom axes titles optional.
}
\usage{
xplot(x,y,xlab=NULL,ylab=NULL,main=NULL,fill=T)
}
\arguments{
  \item{x}{Variable 1}
  \item{y}{Variable 2}
  \item{xlab}{Label for the x-axis, in quotes}
  \item{ylab}{Label for the y-axis, in quotes}
  \item{main}{Label for the plot, in quotes}
  \item{fill}{Use gray fill for density plots? (T or F)}
}

\examples{
# random numbers from a normal distribution
ex1<-rnorm(1000)
# random numbers from an exponential distribution
ex2<-rexp(1000)

xplot(ex1,ex2)
}