\name{periodogram}
\alias{periodogram}
\title{Simple periodogram}
\description{
Calculate periodogram for stratigraphic series
}
\usage{
periodogram(dat,padfac=2,demean=T,detrend=F,nrm=1,background=0,medsmooth=0.2,
             bc=F,output=0,f0=F,fNyq=T,xmin=0,xmax=Nyq,pl=1,genplot=T,check=T,
            verbose=T)
}
\arguments{
  \item{dat}{Stratigraphic series to analyze. First column should be location (e.g., depth), second column should be data value.}
  \item{padfac}{Pad with zeros to (padfac*npts) points, where npts is the original number of data points. padfac will automatically promote the total padded series length to an even number, to ensure the Nyquist frequency is calculated. padfac must be >= 1.}
  \item{demean}{Remove mean from data series? (T or F)}
  \item{detrend}{Remove linear trend from data series? (T or F)}
  \item{nrm}{Power normalization: 0 = no normalization; 1 = divide Fourier coefficients by npts.}
  \item{background}{Estimate noise model background spectrum and confidence levels? (0= No, 1= AR1, 2= Power Law, 3= Robust AR1)}
  \item{medsmooth}{Robust AR1 median smoothing parameter (1 = use 100\% of spectrum; 0.20 = use 20\%)}
  \item{bc}{Apply Bonferroni correction? (T or F)}
  \item{output}{Return output as new data frame? 0= no; 1= frequency, amplitude, power, phase (+ background fit and confidence levels if background selected); 2= frequency, real coeff., imag. coeff. If option 2 is selected, all frequencies from negative to positive Nyquist are returned, regardless of the f0 and fNyq settings; 3= noise model coefficients (if background selected); 4= option 1+3}
  \item{f0}{Return results for the zero frequency? (T or F)}
  \item{fNyq}{Return results for the Nyquist frequency? (T or F)}
  \item{xmin}{Smallest frequency for plotting.}
  \item{xmax}{Largest frequency for plotting.}
  \item{pl}{Power spectrum plotting: 1 = log power, 2 = linear power}
  \item{genplot}{Generate summary plots? (T or F)}
  \item{check}{Conduct compliance checks before processing? (T or F) In general this should be activated; the option is included for Monte Carlo simulation.}
  \item{verbose}{Verbose output? (T or F)}
}

\details{
The power spectrum normalization approach applied here (nrm=1) divides the Fourier coefficients 
by the number of points (npts) in the stratigraphic series, which is equivalent to dividing the 
power by (npts*npts). The (npts*npts) normalization has the convenient property whereby
-- for an unpadded series -- the sum of the power in the positive frequencies is equivalent 
to half of variance; the other half of the variance is in the negative frequencies.
}

\seealso{ \code{\link{mtm}} and \code{\link{lowspec}}}


\examples{
# ***** PART 1: Demonstrate the impact of tapering
# generate example series with 10 periods: 100, 40, 29, 21, 19, 14, 10, 5, 4 and 3 ka.
ex=cycles(c(1/100,1/40,1/29,1/21,1/19,1/14,1/10,1/5,1/4,1/3),amp=c(1,.75,0.01,.5,.25,
             0.01,0.1,0.05,0.001,0.01))

# set zero padding amount for spectral analyses 
# (pad= 1 results in no zero padding, pad = 2 will pad the series to two times its original length)
# start with pad = 1, then afterwards evaluate pad=2
pad=1

# calculate the periodogram with no tapering applied (a "rectangular window")
res=periodogram(ex,output=1,padfac=pad)

# save the frequency grid and the power for plotting
freq=res[1]
pwr_rect=res[3]

# now compare with results obtained after applying four different tapers: 
#  Hann, 30% cosine taper, DPSS with a time-bandwidth product of 1, and DPSS 
#  with a time-bandwidth product of 3
pwr_hann=periodogram(hannTaper(ex,demean=FALSE),output=1,padfac=pad)[3]
pwr_cos=periodogram(cosTaper(ex,p=.3,demean=FALSE),output=1,padfac=pad)[3]
pwr_dpss1=periodogram(dpssTaper(ex,tbw=1,demean=FALSE),output=1,padfac=pad)[3]
pwr_dpss3=periodogram(dpssTaper(ex,tbw=3,demean=FALSE),output=1,padfac=pad)[3]

# now plot the results
ymin=min(rbind (log(pwr_rect[,1]),log(pwr_hann[,1]),log(pwr_cos[,1]),log(pwr_dpss1[,1]),
          log(pwr_dpss3[,1]) ))
ymax=max(rbind (log(pwr_rect[,1]),log(pwr_hann[,1]),log(pwr_cos[,1]),log(pwr_dpss1[,1]),
          log(pwr_dpss3[,1]) ))

pl(2)
plot(freq[,1],log(pwr_rect[,1]),type="l",ylim=c(ymin,ymax),lwd=2,ylab="log(Power)",
      xlab="Frequency (cycles/ka)",
      main="Comparison of rectangle (black), cosine (blue) and Hann (orange) taper",
      cex.main=1)
lines(freq[,1],log(pwr_hann[,1]),col="orange",lwd=2)
lines(freq[,1],log(pwr_cos[,1]),col="blue")
points(c(1/100,1/40,1/29,1/21,1/19,1/14,1/10,1/5,1/4,1/3),rep(ymax,10),cex=.5,
       col="purple")

plot(freq[,1],log(pwr_rect[,1]),type="l",ylim=c(ymin,ymax),lwd=2,ylab="log(Power)",
      xlab="Frequency (cycles/ka)",
      main="Comparison of rectangle (black), 1pi DPSS (green) and 3pi DPSS (red) taper",
      cex.main=1)
lines(freq[,1],log(pwr_dpss1[,1]),col="green")
lines(freq[,1],log(pwr_dpss3[,1]),col="red",lwd=2)
points(c(1/100,1/40,1/29,1/21,1/19,1/14,1/10,1/5,1/4,1/3),rep(ymax,10),cex=.5,
       col="purple")


# ***** PART 2: Now add a very small amount of red noise to the series 
#               (with lag-1 correlation = 0.5)
ex2=ex
ex2[2]=ex2[2]+ar1(rho=.5,dt=1,npts=500,sd=.005,genplot=FALSE)[2]

# compare the original series with the series+noise
pl(2)
plot(ex,type="l",lwd=2,lty=3,col="black",xlab="time (ka)",ylab="signal",
      main="signal (black dotted) and signal+noise (red)"); lines(ex2,col="red")
plot(ex[,1],ex2[,2]-ex[,2],xlab="time (ka)",ylab="difference",
      main="Difference between the two time series (very small!)")

# calculate the periodogram with no tapering applied (a "rectangular window")
res.2=periodogram(ex2,output=1,padfac=pad)

# save the frequency grid and the power for plotting
freq.2=res.2[1]
pwr_rect.2=res.2[3]

# now compare with results obtained after applying four different tapers: 
#  Hann, 30% cosine taper, DPSS with a time-bandwidth product of 1, and DPSS 
#  with a time-bandwidth product of 3
pwr_hann.2=periodogram(hannTaper(ex2,demean=FALSE),output=1,padfac=pad)[3]
pwr_cos.2=periodogram(cosTaper(ex2,p=.3,demean=FALSE),output=1,padfac=pad)[3]
pwr_dpss1.2=periodogram(dpssTaper(ex2,tbw=1,demean=FALSE),output=1,padfac=pad)[3]
pwr_dpss3.2=periodogram(dpssTaper(ex2,tbw=3,demean=FALSE),output=1,padfac=pad)[3]

# now plot the results
ymin=min(rbind (log(pwr_rect.2[,1]),log(pwr_hann.2[,1]),log(pwr_cos.2[,1]),
         log(pwr_dpss1.2[,1]),log(pwr_dpss3.2[,1]) ))
ymax=max(rbind (log(pwr_rect.2[,1]),log(pwr_hann.2[,1]),log(pwr_cos.2[,1]),
         log(pwr_dpss1.2[,1]),log(pwr_dpss3.2[,1]) ))

pl(2)
plot(freq.2[,1],log(pwr_rect.2[,1]),type="l",ylim=c(ymin,ymax),lwd=2,ylab="log(Power)",
     xlab="Frequency (cycles/ka)",
     main="Comparison of rectangle (black), cosine (blue) and Hann (orange) taper",
     cex.main=1)
lines(freq.2[,1],log(pwr_hann.2[,1]),col="orange",lwd=2)
lines(freq.2[,1],log(pwr_cos.2[,1]),col="blue")
points(c(1/100,1/40,1/29,1/21,1/19,1/14,1/10,1/5,1/4,1/3),rep(ymax,10),cex=.5,
        col="purple")

plot(freq.2[,1],log(pwr_rect.2[,1]),type="l",ylim=c(ymin,ymax),lwd=2,ylab="log(Power)",
      xlab="Frequency (cycles/ka)",
      main="Comparison of rectangle (black), 1pi DPSS (green) and 3pi DPSS (red) taper",
      cex.main=1)
lines(freq.2[,1],log(pwr_dpss1.2[,1]),col="green")
lines(freq.2[,1],log(pwr_dpss3.2[,1]),col="red",lwd=2)
points(c(1/100,1/40,1/29,1/21,1/19,1/14,1/10,1/5,1/4,1/3),rep(ymax,10),cex=.5,
      col="purple")

# ***** PART 3: Return to PART 1, but this time increase the zero padding to 2 (pad=2)
}