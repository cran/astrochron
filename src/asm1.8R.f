c This code is a component of astrochron: An R Package for Astrochronology
c Copyright (C) 2014 Stephen R. Meyers

c Contact Stephen Meyers (smeyers@geology.wisc.edu) for information on
c updates. 
c
c Legal Notice: ASM is a program to calculate average spectral misfit 
c and Monte Carlo spectra simulations, Copyright (C) 2014 Stephen R. Meyers.
c
c This program is distributed in the hope that it will be useful, 
c but WITHOUT ANY WARRANTY; without even the implied warranty of 
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


      subroutine asm18_R(freq,numfreq,fBerger,iterms,fper,Rayleigh,
     $ Quist,sedmin,sedmax,numsed,linLog,ispecgen,irepl,isetfreq,
     $ iterations,newNumsed,asmSedrate,asm,asmTerms,asmCL)

c initialize variable types
      implicit real*8 (a-h,o-z)

c set maximum number of iterations
      parameter (mxdata=100000)
c set maximum number of sedimentation rates
      parameter (mxsr=500)
c set maximum number of frequencies
      parameter (mxfreq=500)
c dimension arrays
      dimension fBerger(50), freq(mxfreq), freq2(mxfreq)              !SRM: June 6, 2013
      dimension asm(mxsr),asmTerms(mxsr),asmSedrate(mxsr),asmCL(mxsr)
      dimension rAsm(mxsr),saveAsm(mxsr,mxdata)
      dimension rAsmTerms(mxsr),rAsmSedrate(mxsr)
      dimension store(mxdata),sedrate(mxsr),temp(mxsr)
      dimension fper(50)                                              !SRM: June 6, 2013
      dimension freqgen(mxfreq)                                       !SRM: June 6, 2013
      
c set machine precision constant x 10^3
      epsm=1.11022302D-13

c default option for incorporting uncertainty in measured frequency (1=use, 2=ignore)
       ires = 1

c check to see if any of the input frequencies fall below the Rayleigh freq.
       iremove = 0
       do i=1,numfreq
          freqTest=Rayleigh-freq(i)
          if(freqTest.gt.epsm) then
             iremove = iremove + 1
          endif
       end do
c if frequency falls below Rayleigh freq., exclude it from the analysis
       if (iremove.gt.0) then
          numfreq=numfreq-iremove
          do i = 1,numfreq
             freq(i)=freq(i+iremove)
          end do
c          write(0,*) ' '
c          write(0,*) '  Warning: Some frequencies < Rayleigh frequency'
c          write(0,*) '           These frequences will be excluded'
c          write(0,*) '           New number of frequencies=',numfreq
       endif

c generate sedrate array
c  linear scaling option
      if(linLog.eq.0) then
        sedinc = (sedmax-sedmin)/dble(numsed-1)
        do i=1,numsed
           sedrate(i) = sedmin + (dble(i-1)*sedinc)
        end do
      endif
c  logarithmic scaling option
      if(linLog.eq.1) then
        sedinc = (DLOG10(sedmax)-DLOG10(sedmin))/dble(numsed-1)
        do i=1,numsed
          sedrate(i) = DLOG10(sedmin) + (dble(i-1)*sedinc)
          sedrate(i) = 10.d0**sedrate(i) 
        end do
      endif

c check to see if any sedrates exceed Milankovitch detection range
      iStoreMin=1
      iStoreMax=numsed
      berSedmin=(100.d0/Quist)*fBerger(1)
      berSedmax=(100.d0/Rayleigh)*fBerger(iterms)
      isedmin=0
      newNumsed=numsed
      if(sedrate(1).lt.berSedmin.or.sedrate(numsed).gt.
     $   berSedmax) then
           do i=1,numsed
             if(sedrate(i).lt.berSedmin)then
                newNumsed=newNumsed-1
             elseif(isedmin.eq.0.and.sedrate(i).ge.berSedmin)
     $       then
                sedminNew=sedrate(i)
                iStoreMin=i
                isedmin=1
             elseif(sedrate(i).le.berSedmax)then
                sedmaxNew=sedrate(i)
                iStoreMax=i
             elseif(sedrate(i).gt.berSedmax) then
                newNumsed=newNumsed-1
             endif
           end do

         sedmin=sedminNew
         sedmax=sedmaxNew

c         write(0,*) '   - Some sedrates exceeded Milankovitch ',
c     $            'detection range'
c         write(0,*) '     New range of sedimentation rates:',
c     $             sedmin,sedmax
c         write(0,*) ' '
       endif


c calculate average spectral misfit
      call specmisfit(epsm,freq,Rayleigh,sedrate,iStoreMin,
     $          iStoreMax,numfreq,mxfreq,mxsr,fBerger,asm,
     $          asmTerms,asmSedrate,iterms,Quist,fper,ires)

c ************ DETERMINE SIGNIFICANCE LEVELS **************

c initialize random number.
c it calls the fortran 95 routine in file init_random_seed.f95            !SRM: Jan. 13, 2014
       CALL init_random_seed()                                            !SRM: Jan. 13, 2014
c initialize tickers for random frequency generation routine
       iexceedNyq=0
       iaddRay=0

c if using the bootstrap approach
       if(ispecgen.eq.2) then                                             !This whole section SRM: June 6, 2013
c calculate frequency spacing to use in bootstrap sampling
         do i=2,numfreq
          freqgen(i-1)=freq(i)-freq(i-1)
         end do
       endif
       
c   Start model iteration loop
       do iboot=1,iterations

c using the approach in Meyers and Sageman (2007)
        if(ispecgen.eq.1) then                                             !SRM: June 6, 2013
c generate frequencies
c these frequencies were originally determined using a flat                !SRM: Jan. 13, 2014
c white noise ran number generator (Num. Rec. ran2).                       !SRM: Jan. 13, 2014
c now replaced with intrinsic subroutine random_number                     !SRM: Jan. 13, 2014
c generator returns value between 0 and 1, exclusive of 1,                 !SRM: Jan. 13, 2014
c using the KISS algorithm                                                 !SRM: Jan. 13, 2014
35       do i=1,numfreq
           CALL RANDOM_NUMBER(e)
c linearly rescale random value to frequency within the range of interest
c new values are beteween Rayleigh and Nyquist.
          freq2(i)=( (Quist-Rayleigh)*e ) + Rayleigh
         end do

c order frequencies from lowest to highest.
c if error, ierr1 < 0                                                     !SRM: Jan. 13, 2014
         CALL  dlasrtastro('I',numfreq,freq2,ierr1)                       !SRM: June 26, 2014

c make sure the frequencies are seperated by Rayleigh resolution.
c if not add Rayleigh resolution
         do i=2,numfreq
           diff=freq2(i)-freq2(i-1)

50         if(diff.lt.Rayleigh) then
              freq2(i)=freq2(i)+Rayleigh
              iaddRay=iaddRay+1
           endif
  
c check to ensure that Rayleigh spacing is satistified
c this may not be so if previous frequencies were adjusted
           diff=freq2(i)-freq2(i-1)
           if(diff.lt.Rayleigh) goto 50
                  
           if(freq2(i).gt.Quist) then
              iexceedNyq=iexceedNyq+1
c if maximum frequency exceeds Nyquist then start this iteration over
              goto 35
           endif
         end do
c end Meyers and Sageman (2007) section        
        endif
      

c if using the bootstrap approach
        if(ispecgen.eq.2) then                                             !This whole section SRM: June 6-7, 2013
c generate frequencies
c OPTION 1: Pick one random number from flat white noise, scaled to be between Rayleigh and Nyquist
70        if(isetfreq.eq.1) then
            CALL RANDOM_NUMBER(e)
c rescale random value to frequency within the range of interest
c new values are between Rayleigh and Nyquist.
            freq2(1)=( (Quist-Rayleigh)*e ) + Rayleigh
          endif  

c OPTION 2: Pick 'numfreq' random numbers scaled between Rayliegh and Nyquist (as in Monte Carlo 
c simulations; ispecgen.eq.1), then sort and use the lowest one.
          if(isetfreq.eq.2) then
            do i=1,numfreq
             CALL RANDOM_NUMBER(e)
c linearly rescale random value to frequency within the range of interest
c new values are beteween Rayleigh and Nyquist.
             freq2(i)=( (Quist-Rayleigh)*e ) + Rayleigh
            end do
c order frequencies from lowest to highest
c if error, ierr2 < 0                                                     !SRM: Jan. 13, 2014
            CALL  dlasrtastro('I',numfreq,freq2,ierr2)                    !SRM: June 26, 2014
c NOTE: we are not checking to ensure Rayleigh spacing here (as in the Monte Carlo simulations of
c       specgen.eq.1).  The goal here is simply to ensure that we have a fair chance of
c       "anchoring" the bootstrap spectra in the lower portion of the spectrum.
c NOTE: freq2(1) will be preserved, but all other freq2 entries will be written over below.
          endif

c OPTION 3: The lowest frequency is identical to the data spectrum (freq(1))
          if(isetfreq.eq.3) then
            freq2(1)=freq(1)
          endif

c determine remaining frequencies by bootstrap sampling (with or without replacement)
c from the frequency spacing of the vector 'freq' (that is, the input frequencies
c that were obtained from the data).
c if sampling without replacement ('irepl' = 1), freqgen will be rewritten, 
c  so most recalculate each time
          if(irepl.eq.1) then
            do i=2,numfreq
              freqgen(i-1)=freq(i)-freq(i-1)
            end do        
          endif

          do i=2,numfreq
85         CALL RANDOM_NUMBER(e)
c the random number is from 0-1 excluding the end values
c linearly rescale random value to be an integer between 1 and numfreq-1
c note that DINT truncates toward zero, that is, the fractional portion of its magnitude 
c is truncated and its sign is preserved. 
c we will rescale values from 1 to numfreq, and afterwards use DINT to truncate to numfreq-1
           fe=( (numfreq - 1)*e ) + 1
           ife = DINT(fe)
c if bootstrapping without replacement
           if(freqgen(ife).lt.0.d0) goto 85
           freq2(i)=freq2(i-1) + freqgen(ife) 
c if bootstrapping without replacement
           if(irepl.eq.1) freqgen(ife)=-1.d0
           if(freq2(i).gt.Quist) then
              iexceedNyq=iexceedNyq+1
c if maximum frequency exceeds Nyquist then start this iteration over
              goto 70
           endif
          end do
c end boostrap section
        endif
      

c calculate ASM for this iteration
        call specmisfit(epsm,freq2,Rayleigh,sedrate,iStoreMin,
     $          iStoreMax,numfreq,mxfreq,mxsr,fBerger,rAsm,
     $          rAsmTerms,rAsmSedrate,iterms,Quist,fper,ires)

c save misfit values for sorting
        do i=iStoreMin,iStoreMax
          saveAsm(i,iboot)= rAsm(i)
        end do

c end iterations loop
       end do

c now sort ASM model results,
c loop over sedrates 
       do i=iStoreMin,iStoreMax
         do j=1,iterations
c convert 2-dim array to 1-dim array
           store(j)=saveAsm(i,j)
         end do
c sort misfit results
c if error, ierr3 < 0                                                  !SRM: Jan. 13, 2014
         CALL  dlasrtastro('I',iterations,store,ierr3)                 !SRM: June 26, 2014
c the array 'store' contains the probability distribution function
c for the given sedimentation rate
c identify 90%, 95%, 99% confidence levels
         locate=iterations-INT(iterations*0.9)
         CL90=store(locate)
         locate=iterations-INT(iterations*0.95)
         CL95=store(locate)
         locate=iterations-INT(iterations*0.99)
         CL99=store(locate)
         locate=iterations-(INT(iterations*0.5))
         CL50=store(locate)

c calculate probabilities from ASM values
         test11=asm(i)
c initialize iprob for the case when the observed ASM value  !SRM: Feb 2012
c  exceeds all simulated ASM values                          !SRM: Feb 2012
         iprob = iterations                                  !SRM: Feb 2012
c search through the probability distribution
c  identify all values lower than to the measured ASM
          do j=iterations,1,-1
              test22=store(j)
              testit=test22-test11
              if(testit.gt.epsm) iprob=j
          end do
          prob=dble(iterations-iprob)
          prob=prob/dble(iterations)
          prob=(1.d0-prob)*100.d0
c save probability
          asmCL(i)=prob
c output median ASM from simulation
c          write(2,*) sedrate1,CL50          
c end sedimentation rate loop
       end do

c shuffle evaluated sedrates into first 'newNumsed' vector elements, for pass to R
      newNumsed=iStoreMax-iStoreMin+1
      ii=1
      do i=iStoreMin,iStoreMax
         temp(ii)=asmSedrate(i)
         ii=ii+1
      end do
      do i=1,ii-1
         asmSedrate(i)=temp(i)
      end do
      
      ii=1
      do i=iStoreMin,iStoreMax
         temp(ii)=asm(i)
         ii=ii+1
      end do
      do i=1,ii-1
         asm(i)=temp(i)
      end do      

      ii=1
      do i=iStoreMin,iStoreMax
         temp(ii)=asmTerms(i)
         ii=ii+1
      end do
      do i=1,ii-1
         asmTerms(i)=temp(i)
      end do      

      ii=1
      do i=iStoreMin,iStoreMax
         temp(ii)=asmCL(i)
         ii=ii+1
      end do
      do i=1,ii-1
         asmCL(i)=temp(i)
      end do      

      return
      end



c*********************************************************************
      subroutine specmisfit(epsm,freq,Rayleigh,sedrate,
     $          iStoreMin,iStoreMax,numfreq,mxfreq,mxsr,
     $          fBerger,asm,asmTerms,asmSedrate,iterms,
     $          Quist,ffper,ires)

c initialize variable types
      implicit real*8 (a-h,o-z)
      real*8 misfit
c dimension arrays
      dimension freq(mxfreq),freqT(mxfreq)
      dimension misfit(mxfreq)
      dimension fBerger(50)         !SRM: June 6, 2013
      dimension asm(mxsr),asmTerms(mxsr),asmSedrate(mxsr)
      dimension sedrate(mxsr)
      dimension ffper(50)           !SRM: June 6, 2013

c now loop over sedrates
      do i=iStoreMin, iStoreMax

c calculate effective resolution
         resol=rayleigh*sedrate(i)/100.d0
c if we want to incorporate uncertainty in measured frequency
         if(ires.eq.1) resol2 = resol
c this option is only used if we want to ignore uncertainty in measured frequency
         if(ires.eq.2) resol2 = 0.d0
c calculate maximum frequency
         QuistN=Quist*sedrate(i)/100.d0
c determine minimum and maximum period we can assess
         MinBerger=1
         MaxBerger=iterms

        do j=1,iterms
           fmin=fBerger(j) - ( fBerger(j) * ffper(j) )
           fmax=fBerger(j) + ( fBerger(j) * ffper(j) )
           if(resol-fmax.gt.epsm) MinBerger=MinBerger+1
           if(QuistN-fmin.lt.epsm) MaxBerger=MaxBerger-1
        end do

c calculate freqs (cycles/kyr) in measured spectrum
        do j=1,numfreq
c calibrate spatial frequencies
           freqT(j)=freq(j)*sedrate(i)/100.d0
        end do

c loop over predicted frequencies of Berger
        do j=MinBerger,MaxBerger
c find the measured peak that is closest to the predicted value for each
c component
           fit=100000.d0
           do k=1,numfreq
             test=DABS(freqT(k)-fBerger(j))
             if(test.lt.fit) then
               fit=test
               locfit=k
             endif
           end do

c determine spectral misfit (cycles/kyr) for this component. allow for
c uncertainties in measured (calibrated) frequency and orbital target
c frequency

c put uncertainty bars on the orbital target frequency
           fBerMin = fBerger(j) - ( fBerger(j)*ffper(j) )
           fBerMax = fBerger(j) + ( fBerger(j)*ffper(j) )
c put uncertainty bars on the measured (calibrated) frequency
           freqTmin = freqT(locfit) - (0.5d0*resol2)
           freqTmax = freqT(locfit) + (0.5d0*resol2)
c evaluate spectral misfit
           if(fBerger(j).ge.freqt(locfit)) then
                misfit(j) = fBerMin-freqTmax
                if(misfit(j).lt.0.d0) misfit(j) = 0.d0
           endif
           if(fBerger(j).lt.freqt(locfit)) then
                misfit(j) = freqTmin-fBerMax
                if(misfit(j).lt.0.d0) misfit(j) = 0.d0
           endif
        end do

c now calculate misfit parameters
        cumulMisfit=0.d0
        do j=MinBerger,MaxBerger
          cumulMisfit=cumulMisfit+misfit(j)
        end do
        aveMisfit=cumulMisfit/dble(MaxBerger-MinBerger+1)
c save results
        asmSedrate(i)=sedrate(i)
        asm(i)=aveMisfit
        asmTerms(i)=MaxBerger-MinBerger+1

      end do

      return
      end
