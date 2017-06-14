c This code is a component of astrochron: An R Package for Astrochronology
c Copyright (C) 2016 Stephen R. Meyers

      SUBROUTINE trough_R(npts,y,loc,iplat,numtrough,numplat)

      IMPLICIT REAL(8) (A-H,O-Z)
      INTEGER loc,iplat,numtrough,numplat
      DIMENSION y(npts),loc(npts),iplat(npts)
      
c Search for local minima
c Initialize peak and plateau counters
      j=1
      k=1
      do i=1,npts-2
      if (y(i).gt.y(i+1).and.y(i+1).lt.y(i+2)) then
c Write trough
         loc(j) = i+1
         j=j+1
c Check for plateaus
       elseif(y(i).eq.y(i+1).or.y(i+1).eq.y(i+2)) then
c Write message to standard error
         iplat(k) = i+1
         k=k+1
       endif
      end do
      numtrough =  j-1
      numplat = k-1
      return
      end
