\name{cb}
\alias{cb}
\title{Combine Multiple Vectors}
\description{
Take any number of input vectors, bind the columns together, and return as a data.frame
}
\usage{
cb(a,b)
}
\arguments{
  \item{a}{first input vector OR a data frame with >1 column.}
  \item{b}{second input vector OR if a is a data frame with > 1 column, a list of columns to bind.}
}
\examples{
# example dataset
x<-rnorm(100)
dim(x)<-c(10,10)
x<-data.frame(x)

# bind two columns
cb(x[1],x[5])

# bind five columns
cb(x,c(1,2,4,7,9))
}