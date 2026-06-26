\name{hoangPriors}
\alias{hoangPriors}
\title{Priors for the fundamental frequencies of the Solar System, using the formulas in Table 2 of Hoang et al. (2021)}
\description{
Priors for the fundamental frequencies of the Solar System, using the formulas in Table 2 of Hoang et al. (2021). 
This function will return the value of the prior pdf of the i-th fundamental frequency 
(g_i or s_i) at a given age, from the formulas in Table 2 of Hoang et al. 2021.
}
\usage{
hoangPriors(ssftype=NULL,i=NULL,age=NULL,ssfi=NULL)
}
\arguments{
  \item{ssftype}{The type of fundamental frequency to calculate ('g' or 's').}
  \item{i}{Index indicating the fundamental frequency of interest. For 'g', integer values from 1-5 are permitted. For 's', integer values from 1-4 and 6 are permitted.}
  \item{age}{Age in billions of years before present, Ga.}
  \item{ssfi}{Optional values of the fundamental frequencies for which you'd like to determine the prior. If not specified, the time-dependent parameters of the skew normal pdf are reported.}
}

\references{
Hoang, N.H., Mogavero, F., Laskar, J., 2021,
\emph{Chaotic diffusion of the fundamental frequencies in the Solar System}: Astronomy & Astrophysics, 654, A156. 
}

\examples{

# View the priors for the g and s frequencies at 1 Ga, for comparison 
#  with Figure 9 of Hoang et al., 2021

# set up figure
dev.new(height = 4.5, width = 9.5, title = "Fundamental Frequencies")
par(mfrow = c(2, 5))

# age in billions of years before present
  age=1
  
# g1 through g5
  gmin=c(5, 7.40, 16.9, 17.5, 4.2573)
  gmax=c(6, 7.54, 17.7, 18.5, 4.2576)
# loop over g frequencies
  for (i in 1:5)
   {
# compute range of g_i
        x=seq(gmin[i],gmax[i],(gmax[i]-gmin[i])/500)     
# compute pdf and its integral
        pdfRes=hoangPriors(ssftype="g",i,age,x)
# plot prior pdf of g_i
        plot(x,pdfRes$ssfipdf,type="l",lwd=2,col="red",xlab="arcsec",ylab="",main=paste("g",i))
    }
   
# s1, s2, s3, s4, s6 (s5 set as NaN)
  xmin=c(-6.5, -8.5, -19.5, -18.3, NaN ,-26.3489) 
  xmax=c(-4.8, -6, -18, -17.4, NaN, -26.3469)
# loop over s frequencies
  for (i in 1:6)
   {
      if(i != 5) {
# compute range of s_i
        x=seq(xmin[i],xmax[i],(xmax[i]-xmin[i])/500)     
# compute pdf and its integral
        pdfRes=hoangPriors(ssftype="s",i,age,x)
# plot prior pdf of s_i
        plot(x,pdfRes$ssfipdf,type="l",lwd=2,col="red",xlab="arcsec",ylab="",main=paste("s",i))
      }
   } 
   
}