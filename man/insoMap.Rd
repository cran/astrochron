\name{insoMap}
\alias{insoMap}
\title{generate an insolation map using Laskar et al. (2004) insolation calculations}
\description{
Generate insolation maps using Laskar et al. (2004) insolation calculations.  
This is a wrapper function to calculate insolation using the palinsol package, 
following astrochron syntax.}
\usage{
insoMap(t=0,pick=0,S0=1365,output=F,verbose=T) 
}
\arguments{
  \item{t}{Time for evaluation, in kyr (negative for future, positive for past; valid from -21000 to 51000)}
  \item{pick}{Interactively pick (1) latitude to examine seasonal insolation at that latitude, or (2) day to examine insolation across latitude}
  \item{S0}{Solar constant, W/m^2}    
  \item{output}{Output insolation matrix (if pick=0) or selected insolation transect (if pick = 1 or 2)? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}

\details{
This wrapper function generates maps with insolation calculations based on 
the Laskar et al. (2004) astronomical parameters, using the 'palinsol' package (Crucifix, 2016).

Option 'pick' will allow you to extract insolation by (1) latitude or (2) season.
}


\references{

Crucifix, M., 2016,
\emph{palinsol: Insolation for Palaeoclimate Studies. R package version 0.93}: https://CRAN.R-project.org/package=palinsol

J. Laskar, P. Robutel, F. Joutel, M. Gastineau, A.C.M. Correia, and B. Levrard, B., 2004, 
\emph{A long term numerical solution for the insolation quantities of the Earth}: 
Astron. Astrophys., Volume 428, 261-285.
}