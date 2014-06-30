c This code is a component of astrochron: An R Package for Astrochronology
c Copyright (C) 2014 Stephen R. Meyers

      SUBROUTINE DEMEAN (ipts,y)
      real*8   y(*)                                         
      real*8   mean
c calculate mean
      mean=0.d0
      do i = 1,ipts                                          
        mean = mean + y(i)
      end do
      mean = mean / dfloat(ipts)
c demean series 
      do i = 1,ipts                                                 
        y(i) = y(i) - mean 
      end do
      return
      end
