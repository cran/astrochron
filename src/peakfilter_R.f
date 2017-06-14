c This code is a component of astrochron: An R Package for Astrochronology
c Copyright (C) 2016 Stephen R. Meyers

      SUBROUTINE peakfilter_R(numpeak,nfreq,tbwRay,siglevel,
     &  freqloc,probmax,freq,background,pwr,cl,loc,nout)

      IMPLICIT REAL(8) (A-H,O-Z)
      INTEGER numpeak,nfreq,freqloc,loc,nout
      DIMENSION freqloc(numpeak),probmax(numpeak)
      DIMENSION freq(nfreq),background(nfreq),pwr(nfreq)
      DIMENSION cl(nfreq),loc(numpeak)


c numpeak = number of significant F-test peaks
c nfreq = number of frequencies in spectrum
c tbwRay = ralyeigh frequency * time-bandwidth product
c siglevel = signifiance level for filtering
c freqloc = index for each F-test peak
c probmax = probabilities associated with freqloc
c freq = frequencies associated with background and cl
c background = spectral background estimate
c pwr = power
c cl = confidence level from noise model
c loc = index for filter significant peaks
c nout = number of peaks returned

c loop over F-test peaks (these already filtered at siglevel)
      ii=1
      do i = 1, numpeak
        ifreq = freqloc(i)
        itest=0
c examine power spectrum +/- tbw*Ray to see if there is power at required noise CL;   
c also require that F-test peak is on a power peak relative the local background
        do j = 1, nfreq
           bandlow = freq(ifreq) - tbwRay
           bandhigh = freq(ifreq) + tbwRay
           if(freq(j).ge.bandlow.and.freq(j).le.bandhigh) then
             if(cl(j).ge.siglevel) then
               if(pwr(ifreq).ge.background(ifreq)) itest=1 
             endif
           endif  
        end do  
        if(itest.eq.1) then
           loc(ii) = ifreq
           ii=ii+1
        endif  
      end do
      nout=ii-1
      
      return
      end
