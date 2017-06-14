c This code is a component of astrochron: An R Package for Astrochronology
c Copyright (C) 2017 Stephen R. Meyers

c mwin_R.f: code to determine start and end points for a moving window of fixed
c           duration (e.g. 500 kiloyears). for use with unevenly sampled data.
c
c           input is a single column (vector) with time/depth/height locations.  
c           must be sorted into increasing order. no duplicate values permitted.
c           routine ticks up one point at a time.
c
c           output is vector with start points and vector with end points for each
c           'window'.  to be executed from R. (SRM: March 15, 2013).

c removed 'stop' statement, replaced with return. SRM: June 21, 2013
c additional modification for R compliance. SRM: May 24, 2017; June 13, 2017

       subroutine mwin_R(ipts,winsize,x,npts,n1,n2,avex,midx1,midx2)

c ipts = number of data points
c winsize = size of moving window
c x = input vector of time/height/depth locations (sorted to increasing order)
c n1 = vector with [data point number] for start of each window
c n2 = vector with [data point number] for end of each window
c npts = number of elements in n1 and n2
c avex = vector with mean location for each window
c midx1 = vector with middle location for each window
c midx2 = vector with midpoint for each window (in case window was not actually winsize)

c define variables and constants, set up arrays
       implicit real(8) (a-h,o-z)                      ! this modified for R compliance   
       real(8) midx1,midx2                             ! this modified for R compliance   
       DIMENSION x(ipts),midx1(ipts),midx2(ipts),avex(ipts)
       DIMENSION n1(ipts),n2(ipts)

c set machine precision constant x 10^3
      epsm=1.11022302D-13

c determine when to stop moving window analysis
       stopwin=x(ipts)-winsize
c determine where to start and end first window
       winstart=x(1)
       winend=x(1)+winsize

c npts will be the number of points in the output vectors n1 and n2
       npts=0
c loop over data, calculate moving average
       do i=1,ipts
c intialize variables, arrays
         avex(i)=0.d0
         n1(i) = i
c identify data in current window 
         do j=i,ipts
           if(x(j)-winend.le.epsm) iout=j
         end do
         n2(i)=iout
         npts=npts+1
         do j=i,iout
           avex(i)=avex(i)+x(j)
c keep track of the last depth
           xhold=x(j)
         end do
c kpts stores the number of points in the moving window
         kpts=iout-i+1
c also calculate the average hieght/depth
         avex(i)=avex(i)/dfloat(kpts)

c determine depth/height location
c in 1st scenario, assign data to the middle of the the 'winsize' window
         midx1(i)=(winend+winstart)/2.d0
c in 2nd scenario, adjust in case window was not actually winsize
         midx2(i)= xhold - ( (xhold-x(i))/2.d0 )
c in 3rd scenario, assign data to average x (we've already calculated this)
c set start and end values of next window
         winstart=x(i+1)
         winend=x(i+1) + winsize
c check to see if we are at the end of the data series
         if(x(i+1)-stopwin.ge.epsm) then
           goto 40
         endif
c end big loop
       end do
40     continue

       return 
       end
