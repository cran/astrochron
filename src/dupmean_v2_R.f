c This code is a component of astrochron: An R Package for Astrochronology
c Copyright (C) 2016 Stephen R. Meyers

      SUBROUTINE dupmean_R(ipts,x,y,npts,xx,yy)     ! this added for R
c define variables
      implicit real(8) (a-h,o-z)                    ! this modified for R compliance
      dimension x(ipts),y(ipts),xx(ipts),yy(ipts)

c search data for dups
      irmpts=0
      do i=1,ipts
        k=1
        if (i+irmpts.ge.ipts) then
          yy(i)=y(i+irmpts)
          xx(i)=x(i+irmpts)
          goto 100
        endif
        yy(i)=y(i+irmpts)
        xx(i)=x(i+irmpts)
           do j=i+irmpts+1,ipts
             if(xx(i).eq.x(j)) then
       	        yy(i)=yy(i)+y(j)
                k=k+1
                irmpts=irmpts+1
              endif
           end do
           if (k.gt.1) yy(i)=yy(i)/dfloat(k)  ! this changed for R
       end do

100     npts=ipts-irmpts            ! added for R
        return                      ! added for R
        end
