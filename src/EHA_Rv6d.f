c This code is a component of astrochron: An R Package for Astrochronology
c Copyright (C) 2016 Stephen R. Meyers
c
c Contact Stephen Meyers (smeyers@geology.wisc.edu) for information on
c updates. 
c
c Legal Notice: EHA is a program to perform evolutive harmonic analysis,
c and evolutive power spectral analysis, developed by Stephen R. Meyers. 
c The program includes subroutines from a range of sources (see subroutines for 
c further information).
c
c This program is distributed in the hope that it will be useful, 
c but WITHOUT ANY WARRANTY; without even the implied warranty of 
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

      subroutine EHA_Rv6(npts,inc,nspec,nfreq,dt,imean,itrend,
     &  ifstart,ifend,bufx,tbw,kmany,newpts,tapers,evals,f,amp,phase,
     &  ftest,power,height,ier)

c INPUT:
c npts = number of points in moving window
c inc = number of points to increment each window
c nspec = number of windows
c nfreq = number of frequencies
c dt = sample spacing
c imean = demean window (1=yes, 2=no)
c itrend = remove linear trend from window (1=yes, 2=no)
c ifstart = index fo minimum frequency for analysis
c ifend = index of maximum frequency for analysis
c bufx = data value
c tbw = mtm time-bandwidth product
c kmany = number of tapers
c newpts = padded length
c tapers = a vector (npts x kmany long), contining the dpss tapers (not normalized to RMS=1)
c evals = eigenvalues for each dpss in tapers

c OUTPUT:
c f = frequency output
c amp = amplitude output
c phase = phase output
c ftest = ftest output
c power = power output
c height = height output
c ier = returns 0 if FFT conducted correctly, 1 if error

c====================================================================== 
c     STEP 1: Define variable types, constants, and set up arrays.
c====================================================================== 
c     Set variable types.
      implicit real(8) (a-h,o-z)                                         
c     Define maximum number of data points in data series.
      parameter (mxdata=200000)
c     Define maximum number of DPSS tapers (max number of fft's).
      parameter (mxdpss=20)
c     Dimension generic work arrays.
      dimension xx(mxdata),yy(mxdata)
c     Dimension dpss arrays. NOTE: dpss has reverse array dimensions
      dimension u10(mxdpss),u20(mxdpss),dpss(mxdata,mxdpss),
     &          evals(mxdpss),bias(mxdpss), 
     &          tapers(mxdata*mxdpss)                         ! added for R
c     Dimension input arrays.
      dimension bufx(mxdata)                                 
c     Dimension adaptive weight arrays.
      dimension dg(mxdata),spw(mxdpss),dcf(mxdpss,mxdata)
c     Dimension amplitude and f-test work arrays.      
      dimension x(mxdata),y(mxdata),x1(mxdpss,mxdata),y1(mxdpss,mxdata)
      dimension f(nfreq),ftest(nspec,nfreq)                   ! added for R
      dimension power(nspec,nfreq),amp(nspec,nfreq)           ! added for R
      dimension phase(nspec,nfreq),height(nspec)              ! added for R
c     Define constants
      pi = 4.d0 * datan(1.d0)
      rad = 180.d0 / pi                                                 


c====================================================================== 
c     STEP 2: Define moving window MTM options.
c====================================================================== 

      if(kmany.gt.mxdpss) then
c        write(6,*) '* Maximum number of tapers set to ',mxdpss
c        stop
         goto 1000            ! return; R -June 21, 2013
      endif
c     First taper to use.
      k1 = 1
c     Last taper to use.
      k2 = kmany          

c     Periodogram frequency resolution (Rayleigh resolution). Thomson
c     (1990) suggests that F-testing works well down to this resolution.
c     We will employ it as the f-test bandwidth.
      fdf = 1.d0 / (dfloat(npts) * dt)
c     Frequency grid due to padding.
      df = 1.d0 / (dfloat(newpts) * dt)
c     The bandwidth resolution (halfwidth) of the MTM power spectra, fw, 
c     is determined by the "time-bandwidth product", tbw = tot * fw, 
c     where tot is the total length of the data window.
c     Calculate length of data window.
      tot = dfloat(npts) * dt
c     Frequency resolution (halfwidth) of MTM power spectra.      
      fw = tbw / tot
      tbwp = tbw/dfloat(npts)
c     Calculate number of points in output spectra.
c     If we have real data, only want positive frequencies (fft symmetry).
c      nfreq = nint((f2-f1)/df) + 1    !now this is determined in R, May 2013
c     Calculate starting point for frequency output.
c      ifstart=nint(f1/df) + 1         !now this is determined in R, May 2013
c     Calculate ending point for frequency output.
c      ifend=nint(f2/df) + 1           !now this is determined in R, May 2013
      
c====================================================================== 
c     STEP 3: Calculate dpss, normalize and determine bias.
c             The dpss are calculated in R using the library 'multitaper',
c             and these tapers are not normalized to RMS=1.
c====================================================================== 
c place tapers in array dpss
        itick=1
        do k = k1,k2
          do n = 1,npts
            dpss(n,k)=tapers(itick)
            itick=itick+1
          end do
        end do
          
c     Calculate bias for DPSS window.
        do k = k1,k2
          bias(k)=1.d0-evals(k)
c e is the limiting broadband bias weight
          e = 1.d0 / dsqrt(evals(k))
c normalize DPSS to rms=1 to preserve power for white process
          RMSnorm=0.d0
          do n=1,npts
             RMSnorm=RMSnorm+(dpss(n,k)*dpss(n,k))
          end do
          RMSnorm=dsqrt(RMSnorm/npts)
          do n=1,npts
             dpss(n,k)=dpss(n,k)/RMSnorm
          end do
        end do

c====================================================================== 
c     Step 4: Compute raw DPSS FFT's.           
c====================================================================== 
c     Use xx and yy as real and imaginary work buffers to hold 
c     DPSS for FFT processing. 
      do k = k1,k2
        do n = 1,npts                                                 
          xx(n) = dpss(n,k)                                       
          yy(n) = 0.d0 
        end do 
c     Pad DPSS with zeros.  Only reason to do this is to ensure
c     we obey FFT length rule. Don't need to pad to get better 
c     resolution at f0.
        call PAD(xx,npts,newpts)
        call PAD(yy,npts,newpts)
c     Compute the FFT using Singleton's FFT.
        call SINGLETON(xx,yy,newpts,newpts,newpts,1,ier)
c     Check for FFT error.
c        if(ier.ne.0) write(6,*) '*  FFT ERROR'
c     Save values of the Fourier coefficients for f = 0.0    
c     in the u10 and u20 buffers. 
        u10(k) = xx(1)   
c     Don't need to save yy, since it should be zero (Thomson,1982)
c     but will store anyway... still used in some calculations below
c     purely for ease of code interpretation.
        u20(k) = yy(1)
      end do 
      
c====================================================================   
c     STEP 5: Calculate (amplitude)**2 response of tapers at f0.               
c====================================================================   
c     Compute the sum of the squared kmany FFTed DPSS at f = 0,  
c     and store in u. This will be used in the least squares regression
c     for determination of the amplitudes and in the F-test.  
      u = 0.d0                                                          
      do k = k1,k2                                                  
c     u is the squared modulus at f(0) for each dpss window. 
c     u20(K) is always zero, so don't need to include in calculation.
c       u = u + u10(k)**2 + u20(k)**2 
        u = u + u10(k)**2
      end do                                 

c===================================================================    
c     STEP 6: Begin evolutive loop, preprocess time series.                    
c===================================================================    
c     Start through loop for evolutive spectra.         
      l1 = -inc + 1                                                     
      do i = 1,nspec                                               
c     Calculate depth / height / time indicator.
c        height(i)=first+(dfloat(i-1)*dt*dfloat(inc))+(dfloat(npts-1)*    !removed 'first' (addressed in R script), May 2013
        height(i)=(dfloat(i-1)*dt*dfloat(inc))+(dfloat(npts-1)*        
     $         dt/2.d0)
        l1 = l1 + inc  
        l2 = l1 + npts - 1  
c     Put window of interest into work buffers x (real) and y (imaginary).
        do l = l1,l2                                                  
          l0 = l - l1 + 1  
          x(l0) = bufx(l) 
          y(l0) = 0.d0 
        end do
c     remove mean value from real components.
c     We do this to remove the static offset (large power at zero frequency).
        if (imean.eq.1) call DEMEAN(npts,x)
c     Remove linear trend from real components.
        if (itrend.eq.1) call DETREND(npts,x)
c     Find the sum of the squares (sig2) of the demeaned data window.
        sig2 = 0.d0  
        do n = 1,npts                                                
          sig2 = sig2 + x(n)**2 
        end do
c     Choose your scaling factor, which is based on the sum of 
c     squares (sig2).  This scaling factor will be used in the 
c     adaptive amplitude routine to estimate the spectral energy 
c     at a given frequency that leaks in from outside the frequency 
c     band of interest (a.k.a. bias).
        sbias = sig2/dfloat(npts**2)
c     Now divide sig2 by npts-1 to get the variance of the data window.
c     This can be used to normalize the power results.
        sig2 = sig2 / dfloat(npts-1) 
 
c====================================================================== 
c     STEP 7: Apply tapers.                                   
c====================================================================== 
c     Begin taper/pad/fft loop.
        do k = 1,kmany
          do n = 1,npts
c     Taper real component.
            xx(n) = x(n) * dpss(n,k)  !SRM
c     Don't need to taper imaginary component (it's zero).
c           yy(n) = y(n) * dpss(n,k) 
            yy(n) = y(n) 
          end do   

c====================================================================== 
c     STEP 8: Compute FFTs of tapered series.
c====================================================================== 
c     Pad time series with zeros to satisfy length requirement of 
c     FFT, and for frequency domain interpolation.
          call PAD(xx,npts,newpts) 
          call PAD(yy,npts,newpts) 
c     Enter FFT to compute the eigencoefficients for tapered series k.
c     Use Singleton's FFT.
          call SINGLETON(xx,yy,newpts,newpts,newpts,1,ier)
c     Check for FFT error.
c          if(ier.ne.0) write(6,*) '*  FFT ERROR'

c     Store these eigencoefficients in work buffers x1 (real)
c     and y1 (imaginary).
          do n = 1,newpts                                               
c normalization
            x1(k,n) = xx(n)/dble(npts) 
            y1(k,n) = yy(n)/dble(npts)
          end do  
        end do  
c     End taper/pad/fft loop.                                                        

c====================================================================== 
c     STEP 9: Adaptive weighting iterative procedure.                          
c====================================================================== 
c     Compute an "adaptive weighted" sum of the kmany eigenspectra      
c     by estimating (via iteration) the weights of the individual 
c     eigenspectra at each frequency.  
c
c     Define the maximum number of iterations to perform when
c     calculating the adaptive amplitude weights.  
        jj=40
c     Define the tolerance level for convergence.
        tol = 0.0003d0 
c     Compute the squared moduli, normalized by the variance. 
c     Restrict analysis to frequencies of interest.                         
        do n = ifstart,ifend 
          ff = dfloat(n - 1) * df ! Modified for R, May 2013, SRM
          do k = k1,k2                                                  
            spw(k) = (x1(k,n)**2 + y1(k,n)**2) / sbias  
          end do  
c     Take the average of the first two estimates. This is the first guess
c     to initialize the iteration procedure.                       
          as = (spw(k1) + spw(k1+1)) / 2.d0  
c     Set the error flag to 1. Will stay at 1 if convergence not reached.
          iflag = 1   
c     Start the iterative procedure.
c     Find coefficients to the selected tolerance level tol for each    
c     of the kmany eigenspectra at the nth frequency.                   
          do j = 1,jj                                                   
            fn = 0.d0    
            fx = 0.d0    
            do k = k1,k2 
c     This is equation 11 of Park et al. (1987). Note that we can determine
c     the bias from energy leakage outside of the frequency 
c     band using the approximation that the bias is a constant fraction
c     of the total variance. The bias array represents this fraction.
c     The variance term was determined during computation of spw.
c     Evals=eigenvalue, bias=1-eigenvalue 
              a1 = dsqrt(evals(k)) * as / (evals(k) * as + bias(k)) 
c     Square the weight estimate, since this is what we'll want to 
c     normalize our spectral estimate to.       
              a1 = a1**2
c     Add up the (weights*spectral_estimates) at each frequency.
              fn = fn + a1 * spw(k)
c     Add up the weights at each frequency.                                
              fx = fx + a1 
            end do 
c     Divide the (weights*spectral_estimates) by the weights at each frequency.
c     This is equation 10 of Park et al. (1987), which determines the estimate 
c     of the spectrum for this iteration.                                             
            ax = fn / fx   
            das = dabs(ax - as)  
c     Check for convergence: if achieved, set error flag to 0.            
            if(das.lt.tol) iflag = 0
            as = ax 
          end do  
          if(iflag.eq.1) then
c            write(6,*) 'No convergence at frequency=',ff ! Modified for R, May 2013, SRM
          endif
                                                                        
c     Now save the final weights to the vector dcf.
          degs = 0.d0  
          do k = k1,k2   
            dcf(k,n) = dsqrt(evals(k)) * as / (evals(k) * as + bias(k))
c     degs is the stability of the estimate, aka the 
c     degrees of freedom at each frequency, associated with the 
c     adaptive amplitude methodology. 
            degs=degs+(dcf(k,n)**2)/(dcf(1,n)**2)  
          end do
          dg(n) = 2.d0 * degs  
        end do

c====================================================================   
c     STEP 10: Calculate the real and imaginary components for 
c              amplitude, which will be used in f-testing.
c====================================================================   
c     Calculate Thomson's (1982) regression Equation (13.5) to get the
c     expected values mu(f) of the eigencoefficients, and store the real
c     and imaginary values temporarily in x and y (respectively).
c     Restict analysis to frequencies of interest.
        do n = ifstart,ifend                   
c          ff = dfloat(n - 1) * df ! Modified for R, May 2013, SRM
c     x and y are the real and imaginary amplitude buffers for the final
c     amplitude estimiate.
          x(n) = 0.d0  
          y(n) = 0.d0  
          do k = k1,k2                                                  
c     Add up all the amplitude estimates from each dpss window.
c     Note: we are performing a complex conjugation here.
c     These results will be the least squares fit of the real and
c     imaginary Fourier coefficients (eigencoefficients). 
c     RECALL that U20(k) is zero!  Just showing for ease of code interp.
              x(n)= x(n) + (x1(k,n)*u10(k) + y1(k,n)*u20(k))          
              y(n)= y(n) + (y1(k,n)*u10(k) - x1(k,n)*u20(k))
          end do  
          x(n) = x(n) / u  
          y(n) = y(n) / u  
        end do
                                                                        
c====================================================================   
c      STEP 11: Perform F-testing
c====================================================================   
c     Compute F-statistics and store the results in the xx buffer for   
c     future examination for significance.                              
c     Restict analysis to frequencies of interest.
c      begin frequency loop
        ick=1
        do n = ifstart,ifend 
c     Find frequency location.                                                
          f(ick) = dfloat(n - 1) * df ! Modified for R, May 2013, SRM    
          backgrnd = 0.d0 
          do k = k1,k2                                                  
c     Estimate the continuous spectrum (background) by subtracting the
c     estimated "amplitude" (harmonic line) means from each 
c     of the eigencoefficients. This is the denomenator of the f-test, 
c     which represents the variance associated with the residuals 
c     (Sum of Squared Error), aka, the variance not explained by the model.
c     RECALL that U20(k) is zero!  Just showing for ease of code interp.
               sum1 = x1(k,n) - (x(n)*U10(k) - y(n)*U20(k))            
               sum2 = y1(k,n) - (x(n)*U20(k) + y(n)*U10(k))            
               backgrnd = backgrnd + sum1**2 + sum2**2
          end do  
c     F-test calculation. 
c     The numerator of the F-test is an estimate of the power in the 
c     line components (aka the variance explained by the model).  
c     The F-test compares this power (model variance) to the 
c     background spectra (residual variance).   
          ftest(i,ick)=(dfloat(kmany)-1.D0)*(x(n)**2+y(n)**2)*u/backgrnd 
  
c====================================================================   
c     STEP 12: Calculate amplitude and phase.
c====================================================================   
c     Amplitude calculation.
          amp(i,ick) = dsqrt(x(n)**2 + y(n)**2) 
          amp(i,ick)= 2.d0*dble(npts)*amp(i,ick)
c     Phase calculation.
          px = x(n) 
          py = y(n)
          phase(i,ick) = datan2(py,px) * rad 

c====================================================================   
c     STEP 13: Calculate power.
c====================================================================  
c     Using adaptive weights.
           pwr1=0.d0
           pwr2=0.d0
           do k= k1,k2
             pwr1 = pwr1 + (dcf(k,n)**2) * (x1(k,n)**2 + y1(k,n)**2)
             pwr2 = pwr2 + dcf(k,n)**2
           end do
           power(i,ick) = pwr1 / pwr2

c increment frequency           
           ick = ick + 1
c end frequency loop               
        end do
c end moving window loop        
      end do  
                                                        
1000  return                                                             
      end   
