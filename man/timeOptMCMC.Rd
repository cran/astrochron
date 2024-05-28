\name{timeOptMCMC}
\alias{timeOptMCMC}
\title{TimeOptMCMC: Evaluation of eccentricity-related amplitude modulation and bundling 
in paleoclimate data ("TimeOpt"), with uncertainties via Markov-Chain Monte Carlo}
\description{
TimeOptMCMC: Evaluation of eccentricity-related amplitude modulation and bundling in 
paleoclimate data ("TimeOpt"; Meyers, 2015), with uncertainties on all fitting parameters 
via Markov-Chain Monte Carlo (MCMC). This function follows the approach of 
Meyers and Malinverno (2018).  MCMC is implemented using the Metropolis-Hastings algorithm.  
Optimization is conducted upon the sedimentation rate (constant within the study interval), 
the fundamental frequencies g1-g5, the precession constant k, and four hyperparameters 
associated with the residuals from the spectral and envelope fit. The priors for the 
k and g's are Gaussian, while other parameters (sedrate, hyperparameters) are 
uniform (uninformative).
}
\usage{
timeOptMCMC(dat,iopt=1,sedmin=0.5,sedmax=5,sedstart=NULL,gAve=NULL,
        gSd=NULL,gstart=NULL,kAve=NULL,kSd=NULL,kstart=NULL,
        rhomin=0,rhomax=0.9999,rhostart=NULL,sigmamin=NULL,sigmamax=NULL,sigmastart=NULL,
        ran=T,fit=1,ftol=0.01,roll=10^3,nsamples=1000,epsilon=NULL,test=F,burnin=-1,
        detrend=T,output=1,savefile=F,genplot=1,verbose=T)
}
\arguments{
  \item{dat}{Stratigraphic series for astrochronologic assessment. First column should be depth or height (in meters), second column should be data value.}
  \item{iopt}{(1) fit power and envelope, (2) fit power only.}
  \item{sedmin}{Minimum sedimentation rate for investigation (cm/ka).}
  \item{sedmax}{Maximum sedimentation rate for investigation (cm/ka).}
  \item{sedstart}{Initial sedimentation rate for MCMC search (cm/ka). Default is 0.5*(sedmin+sedmax). Alternatively, if set to negative number, a random value is selected from the prior distribution.}
  \item{gAve}{Vector which contains the average values for the g1 through g5 fundamental frequencies (arcsec/year). Must be in the following order: g1,g2,g3,g4,g5.}
  \item{gSd}{Vector which contains the standard deviation for the g1 through g5 fundamental frequencies (arcsec/year). Must be in the following order: g1,g2,g3,g4,g5.}
  \item{gstart}{Vector which contains the initial values for the g1 through g5 fundamental frequencies (arcsec/year). Must be in the following order: g1,g2,g3,g4,g5. Default is 0.5*(gmin+gmax). Alternatively, if set to negative number, a random value is selected from the prior distribution.}
  \item{kAve}{Average value for the precession constant (arcsec/year).}
  \item{kSd}{Standard deviation for the precession constant (arcsec/year).}
  \item{kstart}{Initial value for the precession constant (arcsec/year). Default is 0.5*(kmin+kmax). Alternatively, if set to negative number, a random value is selected from the prior distribution.}
  \item{rhomin}{Minimum value for residual lag-1 autocorrelation (for both spectral and envelope fit). Default is 0.}
  \item{rhomax}{Maximum value for residual lag-1 autocorrelation (for both spectral and envelope fit). Default is 0.9999}
  \item{rhostart}{Initial value for residual lag-1 autocorrelation (for both spectral and envelope fit). Default 0.5. Alternatively, if set to negative number, a random value is selected from the prior distribution.}
  \item{sigmamin}{Minimum value for residual sigma (for both spectral and envelope fit).}
  \item{sigmamax}{Maximum value for residual sigma (for both spectral and envelope fit).}
  \item{sigmastart}{Initial value for residual sigma (for both spectral and envelope fit). Default 0.5*(data standard deviation). Alternatively, if set to negative number, a random value is selected from the prior distribution.}
  \item{ran}{Would you like to randomly select the parameter for updating (T), or simultaneously update all the parameters (F)?}
  \item{fit}{Test for (1) precession amplitude modulation or (2) short eccentricity amplitude modulation? Option 2 is not yet functional!}
  \item{ftol}{Tolerance in cycles/ka used to define the precession bandpass. It is added to the highest precession frequency, and subtracted from the lowest precession frequency, to define the half power points for the Taner bandpass filter.}
  \item{roll}{Taner filter roll-off rate, in dB/octave.}  
  \item{nsamples}{Number of candidate MCMC simluations to perform.}
  \item{epsilon}{Vector of dimension 11, which controls how large the jump is between each candidate value, e.g. sedimentation rate. For example, a value of 0.2 will yield maximum jump +/- 10 percent of sedimentation rate range. The vector must be arranged in the the following order: sedrate,k,g1,g2,g3,g4,g5,spec_rho,spec_sigma,env_rho,env_sigma. If NULL, all epsilon values will be assigned 0.2}
  \item{test}{Activate epsilon testing mode? This option will assign all MCMC samples a log-likelihood of unity. This provides a diagnostic check to ensure that the applied epsilon values are sampling the entire range of parameter values. (T or F)}
  \item{burnin}{Threshold for detection of MCMC stability.}
  \item{detrend}{Remove linear trend from data series? (T or F)}
  \item{output}{Which results would you like to return to the console? (0) no output; (1) return all MCMC candidates}
  \item{savefile}{Save MCMC samples to file MCMCsamples.csv? (T or F). If true, results are output after every 1000 iterations (last iterations will not be reported if you do not end on an even thousand!)}
  \item{genplot}{Generate summary plots? (0= none; 1=display all summary plots; 2=also include progress plot during iterations; 3=be quiet, but save all plots as pngs to the working directory)}
  \item{verbose}{Verbose output? (T or F)}
}

\details{
TimeOpt is an astronomical testing algorithm for untuned (spatial) stratigraphic data. 
The algorithm identifies the sedimentation rate(s) that simultaneously optimizes: 
(1) eccentricity amplitude modulations within the precession band, and (2) the 
concentration of spectral power at specified target astronomical periods.

This version of TimeOpt uses MCMC via Metropolis-Hastings to estimate the parameters 
and their uncertainties. The priors for the k and g's are Gaussian, while the other 
parameters (sedrate, hyperparameters) are uniform (uninformative). 

When ran=T, the following approach is used to select the parameter to modify:

0.25 probability of changing sedimentation rate

0.25 probability of changing k

0.30 probability of changing g1,g2,g3,g4,g5 (simultaneously)

0.10 probability of changing sigma_spec,rho_spec (simultaneously)

0.10 probability of changing sigma_env,rho_env (simultaneously)
  
This is motivated by sensitivity tests, and the fact that we are most interested 
in g, k and s;  moving each group of parameters (sedrate, k or g's) has specific
consequences we can isolate. 


For additional guidance on the application of TimeOptMCMC, please see the guide "Bayesian Inversion with TimeOptMCMC" (Meyers, 2024), and Meyers and Malinverno (2018).
}

\references{
S.R. Meyers, 2024,
\emph{Bayesian Inversion with TimeOptMCMC}: EarthArXiv, https://doi.org/10.31223/X5DQ3H

S.R. Meyers, 2015,
\emph{The evaluation of eccentricity-related amplitude modulation and bundling in paleoclimate data: An inverse approach for astrochronologic testing and time scale optimization}: Paleoceanography, 30, doi:10.1002/2015PA002850.

S.R. Meyers and A. Malinverno, 2018,
\emph{Proterozoic Milankovitch cycles and the history of the solar system}: Proceedings of the National Academy of Sciences, www.pnas.org/cgi/doi/10.1073/pnas.1717689115.
}