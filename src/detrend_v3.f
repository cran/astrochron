c This code is a component of astrochron: An R Package for Astrochronology
c Copyright (C) 2023 Stephen R. Meyers

c changed dfloat to dble for CRAN compliance. SRM: June 24, 2023

      SUBROUTINE DETREND (ipts,y)
c determine least squares fit.
      real(8)   y(*)                   ! this modified for R compliance                            
      real(8)   m,mn,md,b,aline        ! this modified for R compliance 
      real(8)   xx,xy,sumx,sumy        ! this modified for R compliance 
c intialized least squares variables
      xx = 0.0d0 
      sumx = 0.0d0 
      xy = 0.0d0 
      sumy= 0.d0
c calculate running totals
      do i = 1,ipts                                                  
c sum of x^2
        xx = xx + (dble(i)*dble(i))                          
c sum of x
        sumx = sumx + dble(i) 
c sum of y
        sumy = sumy + y(i)
c sum of xy
        xy = xy + ( dble(i) * y(i) ) 
      end do
c calculate the slope (Mendenhall and Sincich)
      mn=xy - ((sumx*sumy)/dble(ipts))
      md= xx - (dble(ipts)*(sumx/dble(ipts))**2)
      m=mn/md
c calculate the y-intercept (Mendenhall and Sincich)
      b=(sumy/dble(ipts)) - (m * sumx/dble(ipts))
c now removed the linear trend
      do i = 1,ipts
        aline = (m*dble(i)) + b 
        y(i)= y(i) - aline
      end do
      return                                                           
      end                     
