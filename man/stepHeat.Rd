\name{stepHeat}
\alias{stepHeat}
\title{Ar/Ar Geochronology: Generate an Ar/Ar age spectrum and calculate step-heating plateau age.}
\description{
The stepHeat function will evaluate data from stepwise heating experiments, producing an Ar/Ar age spectrum, 
a weighted mean age with uncertainty, and other helpful statistics/plots (with interactive graphics for data culling). 
The function includes the option to generate results using the approach of IsoPlot 3.70 (Ludwig, 2008) or ArArCALC (Koppers, 2002).}
\usage{
stepHeat(dat,unc=1,lambda=5.463e-10,J=NULL,Jsd=NULL,CI=2,cull=T,del=NULL,output=F,
         idPts=T,size=NULL,unit=1,setAr=95,color="black",genplot=T,verbose=T)
}
\arguments{
  \item{dat}{dat must be a data frame with six columns, as follows: (1) \%Ar39 released, (2) date, (3) date uncertainty (one or two sigma), (4) K/Ca, (5) \%Ar40*, (6) F, and (7) F uncertainty (one or two sigma). NOTE: F is the ratio Ar40*/Ar39K (see Koppers, 2002).}
  \item{unc}{What is the uncertainty on your input dates? (1) one sigma, or (2) 2 sigma. DEFAULT is one sigma. This also applies to the F uncertainty, and the J-value uncertainty (if specified) }
  \item{lambda}{Total decay constant of K40, in units of 1/year. The default value is 5.463e-10/year (Min et al., 2000).}
  \item{J}{Neutron fluence parameter}
  \item{Jsd}{Uncertainty for J-value (neutron fluence parameter; as one or two sigma)}
  \item{CI}{Which convention would you like to use for the 95\% confidence intervals? (1) ISOPLOT (Ludwig, 2008), (2) ArArCALC (Koppers, 2002)}
  \item{cull}{Would you like to remove one or more dates using the graphical interface? (T or F)}
  \item{del}{A vector of indices indicating dates to remove from weighted mean calculation. If specified, this takes precedence over cull.}
  \item{output}{Return weighted mean results as new data frame? (T or F)}
  \item{idPts}{Identify datum number on each point? (T or F)}
  \item{size}{Multiplicative factor to increase or decrease size of symbols and fonts. The default is 1.4}
  \item{unit}{The time unit for your results. (1) = Ma, (2) = Ka}
  \item{setAr}{Set the \%Ar40* level to be illustrated on the plot. The default is 95\%.}
  \item{color}{Color to use for symbols. Default is black.}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}

\details{
This function performs weighted mean age calculations, including estimation of age uncertainties, 
mean square weighted deviation, and probability of fit. 

The following plots are produced:

(1) \%Ar40* versus \%Ar39 released
(2) K/Ca versus \%Ar39 released
(3) Ar/Ar age spectrum, with 2 sigma uncertainties for each step, and weighted mean with 95\% confidence interval (in red)

If the J-value and its uncertainty are input, stepHeat will calculate and include the uncertainty associated with J. The uncertainty is calculated and propagated following equation 18 of Koppers (2002).

}

\references{
A.A.P. Koppers, 2002,
\emph{ArArCALC- software for 40Ar/39Ar age calculations}: Computers & Geosciences, v. 28, p. 605-619.

K.R. Ludwig, 2008,
\emph{User's Manual for Isoplot 3.70: A Geochronological Toolkit for Microsoft Excel}: Berkeley Geochronology Center Special Publication No. 4, Berkeley, 77 p.

I. McDougall and T.M. Harrison, 1991,
\emph{Geochronology and Thermochronology by the 40Ar/39Ar Method}: Oxford University Press, New York, 269 pp.

K. Min, R. Mundil, P.R. Renne, and K. Ludwig, 2000,
\emph{A test for systematic errors in 40Ar/39Ar geochronology through comparison with U/Pb analysis of a 1.1-Ga rhyolite}: Geochimica et Cosmochimica Acta, v. 64, p. 73-98.

I. Wendt and C. Carl, 1991,
\emph{The statistical distribution of the mean squared weighted deviation}: Chemical Geology, v. 86, p. 275-285.
}

\seealso{ \code{\link{wtMean}} }

\examples{
\donttest{
# Check to see if this is an interactive R session, for compliance with CRAN standards.
# YOU CAN SKIP THE FOLLOWING LINE IF YOU ARE USING AN INTERACTIVE SESSION.
if(interactive()) {

# Sample MT-09-09 incremental heating Ar/Ar data from Sageman et al. (2014).
perAr39 <- c(4.96,27.58,19.68,39.9,6.25,1.02,0.42,0.19)
age <- c(90.08,89.77,89.92,89.95,89.89,89.55,87.71,86.13)
sd <- c(0.18,0.11,0.08,0.06,0.14,0.64,1.5,3.22)
KCa <- c(113,138,101,195,307,27,17,24)
perAr40 <- c(93.42,99.42,99.64,99.79,99.61,97.99,94.64,90.35)
Fval <- c(2.148234,2.140643,2.144197,2.145006,2.143627,2.135163,2.090196,2.051682)
Fsd <- c(0.00439,0.00270,0.00192,0.00149,0.00331,0.01557,0.03664,0.07846)
ex <- data.frame(cbind(perAr39,age,sd,KCa,perAr40,Fval,Fsd))

stepHeat(ex)

# plot without points identified
stepHeat(ex,size=0,idPts=F,cull=F)

}
}
}