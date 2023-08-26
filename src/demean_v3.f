c This code is a component of astrochron: An R Package for Astrochronology
c Copyright (C) 2023 Stephen R. Meyers

c changed dfloat to dble for CRAN compliance. SRM: June 24, 2023

      SUBROUTINE DEMEAN (ipts,y)
      real(8)   y(*)             ! this modified for R compliance                                
      real(8)   mean             ! this modified for R compliance
c calculate mean
      mean=0.d0
      do i = 1,ipts                                          
        mean = mean + y(i)
      end do
      mean = mean / dble(ipts)
c demean series 
      do i = 1,ipts                                                 
        y(i) = y(i) - mean 
      end do
      return
      end
