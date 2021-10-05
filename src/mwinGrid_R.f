c This code is a component of astrochron: An R Package for Astrochronology
c Copyright (C) 2020 Stephen R. Meyers

c mwinGrid_R.f: code to determine start and end indices for a moving window of fixed
c           duration (e.g. 500 kiloyears), using a fixed sampling grid (space or time)
c           that is evenly spaced.
c
c           input is a single column (vector) with time/depth/height locations.  
c           must be sorted into increasing order. no duplicate values permitted.
c
c           output is vector with start points and vector with end points for each
c           'window'.  to be executed from R. (SRM: August 24-30, 2020).

c this modified from mwin_R_v3.f

       subroutine mwinGrid_R(ipts,x,start,iGridPts,step,winsize,n1,n2)

c ipts = number of data points in x
c x = input vector of time/height/depth locations (sorted to increasing order)
c begin = start moving window at this time/height/depth
c iGridPts = number of points in the output grid
c step = desired grid spacing (space or time) 
c winsize = size of moving window (space or time)
c n1 = vector with [data point number] for start of each window
c n2 = vector with [data point number] for end of each window


c define variables and constants, set up arrays
       implicit real(8) (a-h,o-z)                      ! this modified for R compliance   
       DIMENSION x(ipts)
       DIMENSION n1(iGridPts),n2(iGridPts)

c set machine precision constant x 10^3
c       epsm=1.11022302D-13

c move through the grid      
       do k=1,iGridPts
c determine where to start and end window
         winstart = start + step*(k-1)
         winend = winstart + winsize
c loop over data, identify data points in current window
         do i=ipts,1,-1
c epsm option removed to capture first point
c           if(x(i)-winstart.ge.epsm) iin=i
            if(x(i).ge.winstart) iin=i
         end do
         do i=1,ipts
c           if(x(i)-winend.le.epsm) iout=i
            if(x(i).le.winend) iout=i
         end do  
         n1(k)=iin
         n2(k)=iout
       end do

       return 
       end
