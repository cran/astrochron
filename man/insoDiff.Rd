\name{insoDiff}
\alias{insoDiff}
\title{insolation difference map using Laskar et al. (2004) insolation calculations}
\description{
Generate map that show difference in insolation between two times
using Laskar et al. (2004) insolation calculations.  
This is a wrapper function to calculate insolation using the palinsol package, 
following astrochron syntax.}
\usage{
insoDiff(t1,t2,S0=1365,verbose=T)
}
\arguments{
  \item{t1}{Time 1, in kyr (negative for future, positive for past; valid from -21000 to 51000)}
  \item{t2}{Time 2, in kyr (negative for future, positive for past; valid from -21000 to 51000)}
  \item{S0}{Solar constant, W/m^2}    
  \item{verbose}{Verbose output? (T or F)}
}

\details{
This wrapper function generates maps with insolation calculations based on 
the Laskar et al. (2004) astronomical parameters, using the 'palinsol' package (Crucifix, 2016).
}


\references{

Crucifix, M., 2016,
\emph{palinsol: Insolation for Palaeoclimate Studies. R package version 0.93}: https://CRAN.R-project.org/package=palinsol

J. Laskar, P. Robutel, F. Joutel, M. Gastineau, A.C.M. Correia, and B. Levrard, B., 2004, 
\emph{A long term numerical solution for the insolation quantities of the Earth}: 
Astron. Astrophys., Volume 428, 261-285.
}