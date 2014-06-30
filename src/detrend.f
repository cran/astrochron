c This code is a component of astrochron: An R Package for Astrochronology
c Copyright (C) 2014 Stephen R. Meyers

      SUBROUTINE DETREND (ipts,y)
c determine least squares fit.
      real*8   y(*)                                         
      real*8   m,mn,md,b,aline
      real*8   xx,xy,sumx,sumy
c intialized least squares variables
      xx = 0.0d0 
      sumx = 0.0d0 
      xy = 0.0d0 
      sumy= 0.d0
c calculate running totals
      do i = 1,ipts                                                  
c sum of x^2
        xx = xx + (dfloat(i)*dfloat(i))                          
c sum of x
        sumx = sumx + dfloat(i) 
c sum of y
        sumy = sumy + y(i)
c sum of xy
        xy = xy + ( dfloat(i) * y(i) ) 
      end do
c calculate the slope (Mendenhall and Sincich)
      mn=xy - ((sumx*sumy)/dfloat(ipts))
      md= xx - (dfloat(ipts)*(sumx/dfloat(ipts))**2)
      m=mn/md
c calculate the y-intercept (Mendenhall and Sincich)
      b=(sumy/dfloat(ipts)) - (m * sumx/dfloat(ipts))
c now removed the linear trend
      do i = 1,ipts
        aline = (m*dfloat(i)) + b 
        y(i)= y(i) - aline
      end do
      return                                                           
      end                     
