\name{idPts}
\alias{idPts}
\title{Interactively identify points in plot}
\description{
Interactively identify points in x,y plot.
}
\usage{
idPts(dat1,dat2=NULL,ptsize=1,xmin=NULL,xmax=NULL,ymin=NULL,ymax=NULL,
       logx=F,logy=F,plotype=1,annotate=1,output=1,verbose=T)
}
\arguments{
  \item{dat1}{Data frame with one or two columns. If one column, dat2 must also be specified.}
  \item{dat2}{Data frame with one column.}
  \item{ptsize}{Size of plotted points.}
  \item{xmin}{Minimum x-value (column 1) to plot}
  \item{xmax}{Maximum x-value (column 1) to plot}  
  \item{ymin}{Minimum y-value (column 2) to plot}
  \item{ymax}{Maximum y-value (column 2) to plot}  
  \item{logx}{Plot x-axis using logarithmic scaling? (T or F)}
  \item{logy}{Plot y-axis using logarithmic scaling? (T or F)}
  \item{plotype}{Type of plot to generate: 1= points and lines, 2 = points, 3 = lines}
  \item{annotate}{Annotate plot with text indicating coordinates?: 0=none, 1=annotate above point, 2=annotate below point}
  \item{output}{Return identified points as a data frame? (0) no, (1) return x and y, (2) return index, x and y}
  \item{verbose}{Verbose output? (T or F)}
}

\seealso{ \code{\link{delPts}}, \code{\link{iso}}, \code{\link{trim}} and \code{\link{trimAT}} }