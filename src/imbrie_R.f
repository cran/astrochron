c This code is a component of astrochron: An R Package for Astrochronology
c Copyright (C) 2024 Stephen R. Meyers

      SUBROUTINE imbrie_R(npts,x,ipts,Tm2,b2,yinit,y)

      IMPLICIT REAL(8) (A-H,O-Z)
      dimension x(npts),y(npts),dydt(npts)
      dimension Tm2(ipts),b2(ipts)

c  x: the centered (zero mean) insolation, multiplied by -1 (Analyseries convention)
c  y: the resultant ice sheet history
c  yinit: the initial value for y
c  npts: number of elements in x, y, dydt
c  Tm2: time constants
c  b2: non-linearity coefficients
c  ipts: number of elements in Tm2, b2 (must be either 1, or equivalent to npts)

c calculate first dydt
      y(1)=yinit
c following convention of Analyseries
c grow  
      if(x(1).gt.y(1)) dydt(1)= ((1.d0 - b2(1))/Tm2(1))*(x(1)-y(1))
c decay
      if(x(1).le.y(1)) dydt(1)= ((1.d0 + b2(1))/Tm2(1))*(x(1)-y(1))
      if(ipts.eq.1) then
c now calculate for remainder of record
        do i=2,npts
          y(i)= y(i-1) + dydt(i-1) 
          if(x(i).gt.y(i)) dydt(i)= ((1.d0 - b2(1))/Tm2(1))*(x(i)-y(i))
          if(x(i).le.y(i)) dydt(i)= ((1.d0 + b2(1))/Tm2(1))*(x(i)-y(i))   
        end do
      endif  
            
      if(ipts.eq.npts) then
c now calculate for remainder of record
        do i=2,npts
          y(i)= y(i-1) + dydt(i-1) 
          if(x(i).gt.y(i)) dydt(i)= ((1.d0 - b2(i))/Tm2(i))*(x(i)-y(i))
          if(x(i).le.y(i)) dydt(i)= ((1.d0 + b2(i))/Tm2(i))*(x(i)-y(i))   
        end do
      endif  

      return
      end
