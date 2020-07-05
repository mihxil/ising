c--------------------------------------------------------------------------
c correlation
c This function calculates the correlation between two part of the `sample'
c It searches its information int he common arrays zsbnd and iz0s and with
c site1/2 offset1/2 and size1/2 one  indicate which 2 parts it should compare.

      function correlation(site1, offset1, nint1, position1,
     *                     site2, offset2, nint2, position2,
     *                     sgn, length)
      include 'parameters.f'
      integer*2 iz0s(mxxy), nin(mxxy)
      real*8    zsbnd(mxxy, maxbnd)     
      common/mats/ zsbnd, iz0s, nin
      real*8   position1, position2
      real  lastpos1, lastpos2
      integer site1, offset1, site2, offset2, sgn
      integer nint1, nint2
      integer*2 int1, int2, isign
      real *8    correlation, correlation2, dist1, dist2, dist, length
         
c     start values
      int1 = offset1 + 1
      int2 = offset2 + 1
      isign = sgn

      if(nint1 .ge. 1) then
         dist1 = zsbnd(site1, int1) - position1 ! distance of not yet dealt with part until next interface
         lastpos1 = zsbnd(site1, int1)                ! position of 
      else
         dist1 = length
         lastpos1 = length + position1
      endif
      if(nint2 .ge. 1) then
         dist2 = zsbnd(site2, int2) - position2
         lastpos2 = zsbnd(site2, int2)         
      else
         dist2 = length
         lastpos2 = length + position2
      endif

      dist = dist1 ! for the case that we even don't jump into the do-while loop

      do while (      int1 .le. (offset1 + nint1) 
     *          .or.  int2 .le. (offset2 + nint2) )  ! check this, could be wrong
         if(dist1 .le. dist2) then  ! what happens if dist1 .eq. dist2 ? 
            dist = dist1
            dist2 = dist2 - dist1
            ! should something happen with lastintpos2??
            int1 = int1 + 1
            if(int1 .gt. offset1 + nint1) then
               dist1 = position1 + length - lastpos1
               lastpos1 = position1 + length ! superfluous (we shouldn't need it anymore)
            else
               dist1 = zsbnd(site1, int1) - lastpos1
               lastpos1 = zsbnd(site1, int1)
            endif
         else ! everything the same but 1 and 2 interchanged
            dist = dist2
            dist1 = dist1 - dist2
            int2 = int2 + 1
            if(int2 .gt. offset2 + nint2) then 
               dist2 = position2 + length - lastpos2
               lastpos2 = position2 + length ! superfluous
            else
               dist2 = zsbnd(site2, int2) - lastpos2
               lastpos2 = zsbnd(site2, int2)
            endif
         endif
         correlation2 = correlation2 + isign * dist
         isign = -isign                  
      enddo   ! while      

      ! deal with the last resting part:
      correlation2 = correlation2 + isign * dist ! something like this

      correlation = correlation2
      return 
      end ! correlation

