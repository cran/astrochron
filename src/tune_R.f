c This code is a component of astrochron: An R Package for Astrochronology
c Copyright (C) 2014 Stephen R. Meyers

      SUBROUTINE tune_R(ipts,x,ictrl,ctrl,t,tuned)     ! this added for R
c set up variables, constants and arrays
      implicit real*8 (a-h,o-z)
      dimension x(ipts),ctrl(ictrl),t(ictrl)
      dimension tuned(ipts),sedrate(ictrl)

c calculate sedimentation rates based on time-space map.
      do i=2,ictrl
        sedrate(i-1)=(ctrl(i)-ctrl(i-1))/(t(i)-t(i-1))
      end do

c loop over data series and convert height location to 
c new temporal location (vector 'tuned')
c number of sedimentation rates created
      irates=ictrl-1
c initialize tickers to record # points outside time-space map
      itop=0
      ibtm=0
      do j=1,irates
        do i=1,ipts
c first deal with all data that is inside the time-space map
          if(x(i).ge.ctrl(j).and.x(i).le.ctrl(j+1)) 
     $        tuned(i)=(x(i)-ctrl(j))/sedrate(j)+t(j)
        end do
      end do
      do i=1,ipts
c if data series starts before first time control point,
c use sedrate from first conversion block
          if(x(i).lt.ctrl(1)) then
            tuned(i)=(x(i)-ctrl(1))/sedrate(1)+t(1)
            ibtm=ibtm+1
c if data series ends after the last time control point, 
c use sedrate from last conversion block
          elseif(x(i).gt.ctrl(ictrl)) then
            tuned(i)=(x(i)-ctrl(ictrl))/sedrate(irates)+t(ictrl)
            itop=itop+1
          endif
      end do

200   return
      end
