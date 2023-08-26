c This code is a component of astrochron: An R Package for Astrochronology
c Copyright (C) 2023 Stephen R. Meyers
     
c mwincenter_R.f: code to determine start and end points for a moving window of fixed
c           duration (e.g. 500 kiloyears). for use with unevenly sampled data.
c
c           input is a single column (vector) with time/depth/height locations.  
c           must be sorted into increasing order. no duplicate values permitted.
c           routine ticks up one point at a time.
c
c           output is vector with start points and vector with end points for each
c           'window'.  to be executed from R. (SRM: March 15, 2013).
c

c removed 'stop' statement, replaced with return. SRM: June 21, 2013
c additional modification for R compliance. SRM: May 24, 2017
c this version searches +/- 0.5 * winsize around each input datum (x). 
c (SRM: June 13-14, 2017)
c changed dfloat to dble for CRAN compliance. SRM: June 24, 2023

      subroutine mwincenter_R(ipts,winsize,x,npts,n1,n2,avex,midx1,
     $midx2)

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

c determine where to start first window
      do i=ipts,1,-1
         if( (x(i)-winsize/2.d0) - x(1) .ge. -1*epsm) istart=i
      end do
      
c determine where to end last window
      do i=1,ipts
        if(x(i) - (x(ipts)-winsize/2.d0) .le. epsm) iend=i
      end do

c npts is the number of points in the output vectors n1 and n2
      npts=iend-istart+1
            
c loop over data, calculate moving average
      k=1
      do i=istart,iend
c initialize variables, arrays
         avex(k)=0.d0
c identify data in current window 
         winstart=x(i)-winsize/2.d0
         do j=ipts,1,-1
           if(x(j)-winstart. ge. -1*epsm) iin=j
         end do
         n1(k) = iin
         winend=x(i)+winsize/2.d0
         do j=1,ipts
           if(x(j)-winend .le. epsm) iout=j
         end do
         n2(k)=iout
         do j=iin,iout
           avex(k)=avex(k)+x(j)
         end do
c kpts stores the number of points in the moving window
         kpts=iout-iin+1
c also calculate the average hieght/depth
         avex(k)=avex(k)/dble(kpts)

c determine depth/height location
c in 1st scenario, assign data to the middle of the the 'winsize' window
         midx1(k)=x(i)
c in 2nd scenario, adjust in case window was not actually winsize
         midx2(k)= x(iin) + ( x(iout)-x(iin) ) / 2.d0 
c in 3rd scenario, assign data to average x (we've already calculated this)
c end big loop
         k=k+1
      end do

      return 
      end
