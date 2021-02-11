\name{wtMean}
\alias{wtMean}
\title{Ar/Ar Geochronology: calculate weighted mean age, age uncertainty, and other associated statistics/plots (with interactive graphics for data culling).}
\description{
The wtMean function is designed for Ar/Ar Geochronology, but is also useful as a general purpose weighted mean estimator. 
It will calculate weighted mean age, age uncertainty, and other helpful statistics/plots (with interactive graphics for data culling).
The function includes the option to generate results using the approach of IsoPlot 3.70 (Ludwig, 2008) or ArArCALC (Koppers, 2002).
}
\usage{
wtMean(dat,sd=NULL,unc=1,lambda=5.463e-10,J=NULL,Jsd=NULL,CI=2,cull=-1,del=NULL,
        sort=1,output=F,idPts=T,size=NULL,unit=1,setAr=95,color="black",
        genplot=T,verbose=T)
}
\arguments{
  \item{dat}{dat must contain one of the following: (1) a vector of dates/values for weighted mean calculation, (2) a matrix with two columns: date or value and uncertainty (one or two sigma), or (3) a matrix with six columns, as follows: date, date uncertainty (one or two sigma), K/Ca, \%Ar40*, F, and F uncertainty (one or two sigma). NOTE: F is the ratio Ar40*/Ar39K (see Koppers, 2002). See "details" for more information.}
  \item{sd}{Vector of uncertainties associated with each date or value in 'dat', as one or two sigma. This option is ignored if dat has more than one column}
  \item{unc}{What is the uncertainty on your input dates/values? (1) one sigma, or (2) two sigma. DEFAULT is one sigma. This also applies to the F uncertainty, and the J-value uncertainty (if specified) }
  \item{lambda}{Relevant for Ar/Ar only- Total decay constant of K40, in units of 1/year. The default value is 5.463e-10/year (Min et al., 2000).}
  \item{J}{Relevant for Ar/Ar only- Neutron fluence parameter}
  \item{Jsd}{Relevant for Ar/Ar only- Uncertainty for J-value (neutron fluence parameter; as one or two sigma)}
  \item{CI}{Which convention would you like to use for the 95\% confidence intervals? (1) ISOPLOT (Ludwig, 2008), (2) ArArCALC (Koppers, 2002) (see below for details)}
  \item{cull}{Would you like select dates/data with a graphical interface? (0=no, 1=select points to retain, -1=select points to remove)}
  \item{del}{A vector of indices indicating data points to remove from weighted mean calculation. If specified, this takes precedence over cull.}
  \item{sort}{Sort by date/values? (0=no; 1=sort into increasing order; 2=sort into decreasing order)}
  \item{output}{Return weighted mean results as new data frame? (T or F)}
  \item{idPts}{Identify datum number on each point? (T or F)}
  \item{size}{Multiplicative factor to increase or decrease size of symbols and fonts for plot.}
  \item{unit}{Relevant for geochronology only- The time unit for your results. (1) = Ma, (2) = Ka}
  \item{setAr}{Relevant for Ar/Ar only- Set the \%Ar40* level to be illustrated on the plot. The default is 95\%.}
  \item{color}{Color to use for symbols. Default is black.}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{verbose}{Verbose output? (T or F)}
}

\details{
This function performs weighted mean age calculations, including estimation of age uncertainties, 
mean square weighted deviation, and probability of fit, following the approaches used in 
IsoPlot 3.70 (Ludwig, 2008) and ArArCALC (Koppers, 2002). It is also useful as a general purpose
weighted mean estimator.

The function accepts input in three formats:

(1) each date/value and its uncertainty can be entered as individual vectors ('dat' and 'sd').

(2) a two column matrix can be input as 'dat', with each date or value (first column) and its uncertainty (second column).

(3) a six column matrix can be input as 'dat', with each date, its uncertainty, the associated K/Ca value, \%Ar40*, F, and F uncertainty (one or two sigma).
This option must be used if you wish to calculate and include the uncertainty associated with J. The uncertainty is calculated and propagated following equation 18 of Koppers (2002).

The following plots are produced:

(1) A normal Q-Q plot for the dates/values (in essence this is the same as IsoPlot's linearized probability plot).

(2) A cumulative Gaussian plot for the dates/values (a.k.a. cumulative probability plot). This is derived by summing the individual normal distributions for each date/value.

(3) A plot of each date/value with its 2-sigma uncertainties.

In addition, K/Ca and Ar40* data are plotted if provided.

A NOTE regarding confidence intervals: There are two conventions that can be used to calculate the confidence intervals, selected with the option 'CI':

(1) ISOPLOT convention (Ludwig, 2008). When the probability of fit is >= 0.15, the confidence interval is based on 1.96*sigma.  When the probability of fit is < 0.15, the confidence interval is based on t*sigma*sqrt(MSWD). 

(2) ArArCALC convention (Koppers, 2002). When MSWD <=1, the confidence interval is based on 1.96*sigma. When MSWD > 1, the confidence interval is based on 1.96*sigma*sqrt(MSWD).

ADDITIONAL ADVICE: Use the function readMatrix to load your data in R (rather than the function read).
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

\seealso{ \code{\link{stepHeat}} }

\examples{
\dontrun{
# Check to see if this is an interactive R session, for compliance with CRAN standards.
# YOU CAN SKIP THE FOLLOWING LINE IF YOU ARE USING AN INTERACTIVE SESSION.
if(interactive()) {

# Sample NE-08-01 Ar/Ar data from Meyers et al. (2012) supplement
age <- c(93.66,94.75,94.6,94.22,86.87,94.64,94.34,94.03,93.56,93.85,88.55,93.45,93.84,
          94.39,94.11,94.48,93.82,93.81,94.18,93.78,94.41,93.49,95.07,94.19)
sd2<- c(5.83,4.10,8.78,2.5,8.86,3.37,4.63,3.18,8.35,5.73,4.23,2.56,2.3,1.7,3.1,2.78,
         1.62,.92,.98,1.41,1.21,1.38,1.48,0.93)
sd <- sd2/2
wtMean(age,sd)

# to calculate the weighted mean without interactive plots and data culling
wtMean(age,sd,cull=0,output=TRUE,genplot=FALSE,verbose=FALSE)
}
}
}