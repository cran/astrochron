### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2025 Stephen R. Meyers
###
###########################################################################
### function hoangPriors- (SRM: August 2, 2022; June 10, 2024; May 23, 2025)
#      return value of prior pdf of the i-th fundamental Solar system frequency
#      (g_i or s_i), from the formulas in Table 2 of Hoang et al. 2021. 
#
###########################################################################

# This is an R-translation of Alberto Malinverno's function ssfipriorpdf.m
# Input:
# - ssftype: 'g' or 's', a single character
# - i: index of Solar system frequency
# - age: in Ga
# - ssfi: value(s) of Solar system frequency for calculation of the prior
#   pdf value ssfipdf (optional) 
#
# Output:
# - ssfipdf, prior pdf at value(s) ssfi (NaN if ssfi is not in the 
#   argument list)
# - mu, sigma, alpha: parameters of skew normal pdf
#
# Note: if alpha==0, the mean of the skew normal pdf equals mu; otherwise,
#
# delta=alpha/sqrt(1+alpha^2);
# ssfimean=mu+sigma*delta*sqrt(2/pi);

hoangPriors <- function(ssftype=NULL,i=NULL,age=NULL,ssfi=NULL)
{

if(ssftype=="g")
 {
# check input i
    if (i%%1!=0 || i<1 || i>5)  stop('gipriorpdf: i must be an integer between 1 and 5')
# Values in Table 2 of Hoang et al. 2021 for g_1 to g_5
   amu=c(5.759, 7.448, 17.269, 17.896, 4.257454)      # mu_0 in Table 2
   bmu=c(0.006, -0.004, 0.002, 0.005, -2.1E-6)        # time-dependent change of mu_0
   asigmasq=c(3.37E-2, 4.17E-4, 6.63E-3, 6.88E-3, 4.63E-10)  # a in Table 2
   bsigmasq=c(0.52, 0.70, 0.43, 0.41, 0.88)           # exponent b of time-dependent sigmasq
   aalpha=c(-2.25, 1.38, 0, 0, 0)                     # skewness parameter alpha_0 in Table 2
   balpha=c(-0.5, 0.21, 0, 0, 0)                      # time-dependent change of alpha_0
 }

if(ssftype=="s")
 {
# NOTE: mu_0 for s_6 is erroneously set to -2.634787 in Table 2
   if (i%%1!=0 || i<1 || i>6)  stop('gipriorpdf: i must be an integer between 1 and 6')
   if (i==5) stop('ssfipriorpdf: s_5 is fixed to zero and has no prior pdf')
   amu=c(-5.652, -6.709, -18.773, -17.707, NaN, -26.34787)
   bmu=c(-0.032, 0.030, 0.009, 0.013, NaN, 1.5E-5)
   asigmasq=c(2.68E-2, 1.20E-1, 2.86E-2, 1.19E-2, NaN, 1.21E-8)
   bsigmasq=c(0.83, 0.76, 0.56, 0.68, NaN, 0.85)
   aalpha=c(1.12, -2.94, -3.40, -1.73, NaN, 0)
   balpha=c(0.16, -1.23, -0.08, -0.28, NaN, 0)
 }

if(ssftype != "s" && ssftype != "g") stop('ssfipriorpdf: ssftype must be g or s')


# This is an R-translation of Alberto Malinverno's function skewnormalpdf.m
# compute a skew-normal pdf for parameters mu, sigma, 
# alpha at the value(s) in x
skewnormalpdf <- function(mu,sigma,alpha,x)
{
  xnorm=((x-mu)/sigma)
  smallphi=(1/sqrt(2*pi))*exp(-0.5*xnorm^2)
# set up error function
  erf <- function(y) 2 * pnorm(y * sqrt(2)) - 1
  bigphi=0.5*(1+erf((alpha/sqrt(2))*xnorm))
  snpdf=(2/sigma)*smallphi*bigphi

  return(snpdf)
}


# parameters mu, sigma, alpha of prior pdf of i-th ssf
mu=amu[i]+bmu[i]*age
sigmasq=asigmasq[i]*age^(bsigmasq[i])
alpha=aalpha[i]+balpha[i]*age                        # skewness parameter of skew normal pdf
sigma=sqrt(sigmasq)

# compute prior value of pdf for value(s) in ssfi, if requested
if(is.null(ssfi)) ssfipdf=NaN
if(!is.null(ssfi))
 { 
   ssfipdf=skewnormalpdf(mu,sigma,alpha,ssfi)
   if (ssftype=='g' && i==4)                         # add secondary mode to g_4
    {
      mu1=17.6755
      sigmasq1=0.0034
      A1=0.11-0.012*age
      secpdf=skewnormalpdf(mu1,sqrt(sigmasq1),0,ssfi) # normal pdf, alpha=0
      ssfipdf=(1-A1)*ssfipdf+A1*secpdf
    }

   if (ssftype=='s' && i==3)                         # add secondary mode to s_3
    {
      mu1=-18.5256
      sigmasq1=0.0028
      A1=0.023
      secpdf=skewnormalpdf(mu1,sqrt(sigmasq1),0,ssfi)# normal pdf, alpha=0
      ssfipdf=(1-A1)*ssfipdf+A1*secpdf
    }
 }
if(is.null(ssfi)) return(data.frame(mu,sigma,alpha))
if(!is.null(ssfi)) return(data.frame(ssfipdf))
}