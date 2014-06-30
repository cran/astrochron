c This code is a component of astrochron: An R Package for Astrochronology
c Copyright (C) 2014 Stephen R. Meyers

      subroutine PAD(x,npts,new)                                      
      real*8 x(*)                                        
      do j= (npts+1), new
        x(j)=0.d0
      end do
      return                                                            
      end
