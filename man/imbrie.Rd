\name{imbrie}
\alias{imbrie}
\title{Imbrie and Imbrie (1980) ice sheet model}
\description{
An implementation of the Imbrie and Imbrie (1980) ice sheet model
}
\usage{
imbrie(insolation=NULL,Tm=17,b=0.6,times=NULL,initial=0,burnin=100,standardize=T,
       output=T,genplot=1,check=T,verbose=T)
}
\arguments{
  \item{insolation}{Insolation, in ka (negative for future, positive for past). Default is insolation over the past 1000 ka from 65 deg. North, 21 June.}
  \item{Tm}{Vector of mean time constants in ka. Default is 17 ka.  The order of the Tm values should match vectors b and times.}
  \item{b}{Vector of nonlinearity coefficient (a value ranging from 0 to 1). Default is 0.6. The order of the b values should match vectors Tm and times.}
  \item{times}{Vector of start times for each Tm and b listed above. Leave as NULL if you only need to model one Tm and b value.}
  \item{initial}{Initial value for ice volume, relative to centered record. Default is 0.}
  \item{burnin}{Number of points for model burn-in. This is required to achieve stable model results. Default is 100 points.}
  \item{standardize}{Standardize model output to maximum value of one and minimum value of zero? (T or F)}  
  \item{output}{Output model results? (T or F)}
  \item{genplot}{Generate summary plots? (1) plot insolation and ice volume series, (2) plot animated insolation, ice volume and phase portrait.)}
  \item{check}{Conduct compliance checks before processing? (T or F). In general this should be activated; the option is included for Monte Carlo simulation.}
  \item{verbose}{Verbose output? (T or F)}
}

\details{
This function will implement the ice volume model of Imbrie and Imbrie (1980), following the conventions of Paillard et al. (1996).

When using the 'times' vector, consider the following example:

times= c(500,1000)

Tm=c(15,5)

b=c(0.6,0.3)

In this case, a Tm of 15 (b of 0.6) will be applied to model from 0-500 ka,
and a Tm of 5 (b of 0.3) will be applied to model 500-1000 ka.

}

\references{
Imbrie, J., and Imbrie, J.Z., (1980),
\emph{Modeling the Climatic Response to Orbital Variations}: Science, v. 207, p. 943-953.

Lisiecki, L. E., and M. E. Raymo, 2005, 
\emph{A Pliocene-Pleistocene stack of 57 globally distributed benthic d18O records}, Paleoceanography, 20, PA1003, doi:10.1029/2004PA001071.

Paillard, D., L. Labeyrie and P. Yiou, (1996),
\emph{Macintosh program performs time-series analysis}: Eos Trans. AGU, v. 77, p. 379.
}

\examples{
\dontrun{
# make a very simple forcing (on/off)
forcing=cycles(0,end=300)
forcing[50:150,2]=1
plot(forcing,type="l")

# use this forcing to drive the imbrie ice model
# set b=0, Tm = 1
imbrie(forcing,b=0,Tm=1,output=F)

# let's view the evolution of the ice sheet
imbrie(forcing,b=0,Tm=1,output=F,genplot=2)

# now increase the response time
imbrie(forcing,b=0,Tm=10,output=F,genplot=2)

# now model slow growth, fast decay
imbrie(forcing,b=0.5,Tm=10,output=F,genplot=2)

# now make a 100 ka cyclic forcing
forcing=cycles(1/100,end=300)
imbrie(forcing,b=0,Tm=1,output=F,genplot=2)
imbrie(forcing,b=0,Tm=10,output=F,genplot=2)
imbrie(forcing,b=0.5,Tm=10,output=F,genplot=2)
# show burn-in
imbrie(forcing,b=0.5,Tm=10,output=F,genplot=2,burnin=0)

# now examine Malutin Milankovitch's hypothesis: 65 deg N, summer solstice
imbrie(b=0.5,Tm=10,output=F,genplot=2,burnin=900)

# use the ice model output to make a synthetic stratigraphic section
res=imbrie(b=0.5,T=10,output=T,genplot=1,burnin=100)
synthStrat(res,clip=F)

# generate ice model for last 5300 ka, using 65 deg. N insolation, 21 June
# allow b and Tm values to change as in Lisiecki and Raymo (2005):
insolation=getLaskar("insolation")
insolation=iso(insolation,xmin=0,xmax=5300)
#  b is 0.3 from 5300 to 3000 ka, then linearly increases to 0.6 between 3000 and 1500 ka.
#  b is 0.6 from 1500 ka to present.
set_b=linterp(cb(c(1500,3000),c(0.6,0.3)),dt=1)
set_b=rbind(set_b,c(5400,0.3))
#  Tm is 5 ka from 5300 to 3000 ka, then linearly increases to 15 ka between 3000 and 1500 ka.
#  Tm is 15 ka from 1500 ka to present.
set_Tm=linterp(cb(c(1500,3000),c(15,5)),dt=1)
set_Tm=rbind(set_Tm,c(5400,5))
# now run model
ex=imbrie(insolation=insolation,Tm=set_Tm[,2],b=set_b[,2],times=set_b[,1])
# time-frequency analysis of model result
eha(ex,fmax=0.1,win=500,step=10,pad=5000,genplot=4,pl=2)

 }
}