c This code is a component of astrochron: An R Package for Astrochronology
c Copyright (C) 2016 Stephen R. Meyers

      SUBROUTINE peak_R(npts,y,loc,iplat,numpeak,numplat)

      IMPLICIT REAL(8) (A-H,O-Z)
      INTEGER loc,iplat,numpeak,numplat
      DIMENSION y(npts),loc(npts),iplat(npts)
      
c Search for local maxima
c Initialize peak and plateau counters
      j=1
      k=1
      do i=1,npts-2
      if (y(i).lt.y(i+1).and.y(i+1).gt.y(i+2)) then
c Write peaks
         loc(j) = i+1
         j=j+1
c Check for plateaus
       elseif(y(i).eq.y(i+1).or.y(i+1).eq.y(i+2)) then
c Write message to standard error
         iplat(k) = i+1
         k=k+1
       endif
      end do
      numpeak =  j-1
      numplat = k-1
      return
      end
