\name{timeOptTemplate}
\alias{timeOptTemplate}
\title{TimeOpt analysis using variable sedimentation models}
\description{
Evaluation of eccentricity-related amplitude modulations and bundling in paleoclimate data, as in Meyers (2015)
and Meyers (2019), adapted to allow the evaluation of a wide range of variable sedimentation models, 
including: differential accumulation across bedding couplets, linear accumulation rate change, step 
changes in sedimentation rate, etc.
}
\usage{
timeOptTemplate(dat,template=NULL,sedmin=0.5,sedmax=5,difmin=NULL,difmax=NULL,fac=NULL,
         numsed=50,linLog=1,limit=T,fit=1,fitModPwr=T,iopt=3,flow=NULL,fhigh=NULL,
         roll=NULL,targetE=NULL,targetP=NULL,cormethod=1,detrend=T,detrendTemplate=F,
         flipTemplate=F,ncores=1,output=0,genplot=1,check=T,verbose=1)
 }
\arguments{
  \item{dat}{Stratigraphic series for astrochronologic assessment. First column should be depth or height (in meters), second column should be data value.}
  \item{template}{Instantaneous sedimentation rate template to fit. This represents a unitless proportional sedimentation rate history. Default model is a copy of dat, which will be scaled for instantaneous accumulation optimization.}
  \item{sedmin}{Minimum AVERAGE sedimentation rate for investigation (cm/ka).}
  \item{sedmax}{Maximum AVERAGE sedimentation rate for investigation (cm/ka).}
  \item{difmin}{Minimum instantaneous sedimentation rate to investigate (cm/ka).}  
  \item{difmax}{Maximum instantaneous sedimentation rate to investigate (cm/ka). By default, this is ignored, and fac is used.} 
  \item{fac}{Maximum instantaneous accumulation factor. Maximum rate is scaled to each investigated sedrate as fac*sedrate. Default value of 5 is based on experimentation. If larger than this, risk getting into local minimum during fit.}
  \item{numsed}{Number of sedimentation rates to investigate in optimization grid.}
  \item{linLog}{Use linear or logarithmic scaling for sedimentation rate grid spacing? (0=linear, 1=log)}
  \item{limit}{Limit evaluated sedimentation rates to region in which full target signal can be recovered? (T or F)}
  \item{fit}{Test for (1) precession amplitude modulations or (2) short eccentricity amplitude modulations?}
  \item{fitModPwr}{Include the modulation periods in the spectral fit? (T or F)}
  \item{iopt}{Optimize on (1) modulations, (2) spectral power, (3) modulations*spectral power}
  \item{flow}{Low frequency cut-off for Taner bandpass (half power point; in cycles/ka)}
  \item{fhigh}{High frequency cut-off for Taner bandpass (half power point; in cycles/ka)}
  \item{roll}{Taner filter roll-off rate, in dB/octave.  Default value is 10^3.}  
  \item{targetE}{A vector of eccentricity periods to evaluate (in ka). These must be in order of decreasing period, with first value of 405 ka.}
  \item{targetP}{A vector of precession periods to evaluate (in ka). These must be in order of decreasing period.}
  \item{cormethod}{Method used for calculation of correlation coefficient (1=Pearson, 2=Spearman)}
  \item{detrend}{Remove linear trend from data series? (T or F)}
  \item{detrendTemplate}{Remove linear trend from sedimentation rate template? (T or F)}
  \item{flipTemplate}{Flip direction of sedimentation rate template? (T or F)}
  \item{ncores}{Number of cores to use for parallel processing}
  \item{output}{Which results you like to return to console? (0) no output; (1) return sedimentation rate grid, r2; (2) return optimal time series, bandpassed series, Hilbert and fitted periods}
  \item{genplot}{Generate summary plots? (0 = nothing, 1=summary plot, 2=progress + summary plots)}
  \item{check}{Conduct compliance checks before processing? (T or F) In general this should be activated; the option is included for Monte Carlo simulation.}
  \item{verbose}{Verbose output? (0 = nothing, 1=minimal, 2=everything)}
}

\details{
TimeOpt employs a probabilistic linear regression model framework to investigate amplitude
modulation and frequency ratios (bundling) in stratigraphic data, while simultaneously
determining the optimal time scale. This function further develops the method to optimize
upon complex sedimentation templates. The approach is demonstrated below with a series of 
examples. 

The statistical significance of the r^2_opt is determined via Monte Carlo simulation 
(see timeOptSim). See timeOpt for more information on the basic methodology.
}

\value{
if output = 1, a data frame containing the following will be returned:
Sedimentation rate (cm/ka), r-squared value (r^2_envelope, r^2_spectra, or r^2_opt)

if output = 2, a data frame containing the following will be returned:
Time (ka), tuned time series, bandpassed series, envelope, reconstructed model
}

\seealso{ \code{\link{asm}}, \code{\link{eAsmTrack}}, \code{\link{testPrecession}}, \code{\link{timeOpt}}, \code{\link{timeOptSim}}, and \code{\link{timeOptTemplatePlot}}}

\references{
S.R. Meyers, 2015,
\emph{The evaluation of eccentricity-related amplitude modulations and bundling in 
paleoclimate data: An inverse approach for astrochronologic testing and time scale 
optimization}: Paleoceanography, v.30, 1625-1640.

S.R. Meyers, 2019, 
\emph{Cyclostratigraphy and the problem of astrochronologic testing}: Earth-Science Reviews
v.190, 190-223.
}

\examples{
\dontrun{
# EXAMPLE (1): Differential accumulation across bedding couplets
ex=etp(tmin=0,tmax=1000,dt=5,pWt=1,oWt=0,eWt=1,esinw=TRUE)
ex2=diffAccum(ex,0.01,.05)
ex2=linterp(ex2)
# first with the nominal timeOpt approach
timeOpt(ex2,sedmin=1,sedmax=4,numsed=100)
timeOptSim(ex2,sedmin=1,sedmax=4,numsed=100,numsim=2000)
# then with the timeOptTemplate approach
timeOptTemplate(ex2,sedmin=1,sedmax=4,difmin=.5,difmax=6,numsed=100,ncores=2)
timeOptTemplateSim(ex2,sedmin=1,sedmax=4,difmin=.5,difmax=6,numsed=100,numsim=1000,
 ncores=2)


# EXAMPLE (2): Linear sedimentation rate increase
ex=etp(tmin=0,tmax=1000,dt=5,pWt=1,oWt=0,eWt=1,esinw=TRUE)
ex2=sedRamp(ex,srstart=0.01,srend=0.05)
ex2=linterp(ex2)
# first with the nominal timeOpt approach
timeOpt(ex2,sedmin=1,sedmax=4,numsed=100)
# then with the timeOptTemplate approach
# create linear model for input. the magnitude does not matter, it will be rescaled. 
# (it just needs to be a line)
template=ex2; template[2]=ex2[1]
timeOptTemplate(ex2,template=template,sedmin=1,sedmax=4,difmin=.5,difmax=6,numsed=100,ncores=2)
# view optimization procedure (must set ncores=1)
timeOptTemplate(ex2,template=template,sedmin=2.75,sedmax=3.25,difmin=.5,difmax=6,numsed=20,
ncores=1,genplot=2)

# EXAMPLE (3): Step increase in sedimentation rate, from 1 cm/kyr to 2 cm/kyr at 7 meters depth
ex=etp(tmin=0,tmax=1000,dt=5,pWt=1,oWt=0,eWt=1,esinw=TRUE)
ex2=ex
ex2[1]=ex[1]*.01
ex2[141:201,1]=ex2[141:201,1]*2-7
ex2=linterp(ex2)
# first with the nominal timeOpt approach
timeOpt(ex2,sedmin=1,sedmax=4,numsed=100)
# then with the timeOptTemplate approach
# create step model for input. the magnitude does not matter, it will be rescaled. 
template=ex2; template[1:140,2]=1; template[141:261,2]=2
timeOptTemplate(ex2,template=template,sedmin=1,sedmax=4,numsed=100,ncores=2)
# view optimization procedure (must set ncores=1)
timeOptTemplate(ex2,template=template,sedmin=1,sedmax=2,numsed=20,ncores=1,genplot=2)


# EXAMPLE (4): A record with a 100 kyr hiatus at 10 meters depth
ex=etp(tmin=0,tmax=1000,dt=5,pWt=1,oWt=0,eWt=1,esinw=TRUE)
ex2=delPts(ex,del=101:121)
# use a background sedimentation rate of 2 cm/kyr
ex2[1]=0:179*5*0.02
# first evaluate the distorted record with the nominal timeOpt approach
timeOpt(ex2,sedmin=1,sedmax=4,numsed=100)
# then with the timeOptTemplate approach
#  create a constant sedimentation rate template with possible hiatus of unknown
#  duration at 10 m
template=ex2; template[2]=10; template[101,2]=1
timeOptTemplate(ex2,template=template,sedmin=1,sedmax=3,difmax=3,numsed=100,ncores=2)
# now perform a finer grid search near the maximum, using power only
#  notice the oscillatory nature of the power fit.
res=timeOptTemplate(ex2,template=template,sedmin=1.5,sedmax=2,difmax=3,numsed=100,
 ncores=2,iopt=2,output=2)
# compare true eccentricity to TimeOpt-derived eccentricity
pl(2)
plot(ex,type="l",main="True Eccentricity Series",xlab="True Time (kyr)",ylab="")
plot(res[,1],res[,4],type="l",main="Black=TimeOpt precession AM;  Red=TimeOpt eccentricity model",
xlab="TimeOpt derived time (kyr)",ylab="")
lines(res[,1],res[,5],col="red",lwd=2)
}
}