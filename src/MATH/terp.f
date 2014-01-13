c     interpolation and extrapolation by m-th order polynomial
c     maximum m = 3. Change nmax if needed.
c     Input x and y arrays, returns y value y0 at requested x value x0.
c     Dies on error.

      subroutine terp (x, y, n, m, x0, y0)
      implicit double precision (a-h, o-z)

      dimension x(n), y(n)

c     Find out between which x points x0 lies
      i = locat (x0, n, x)
      k = min( max(i-m/2,1) , n-m )
      call polint( x(k), y(k), m+1, x0, y0, dy)

      return
      end

      function locat (x, n, xx)
      integer  u, m, n
      double precision x, xx(n)

c     Binary search for index of grid point immediately below x.
c     Array xx required to be monotonic increasing.
c     Returns
c     0            x <  xx(1)
c     1            x =  xx(1)
c     i            x =  xx(i)
c     n            x >= xx(n)

      locat = 0
      u = n+1

   10 if (u-locat .gt. 1)  then
         m = (u + locat) / 2
         if (x .lt. xx(m))  then
            u = m
         else
            locat = m
         endif
         goto 10
      endif

      return
      end


c     These routines, terp1 and locat1, are special versions to
c     be used with ff2chi, which uses some single and some double
c     precision.  They are the same as the routines in terp.f.

      subroutine terp1 (x, y, n, x0, y0)
      implicit double precision (a-h, o-z)

      real x(n), y(n)

c     Find out between which x points x0 lies
      i = locat1 (x0, n, x)
c     if i < 1, set i=1, if i > n-1, set i=n-1
      i = max (i, 1)
      i = min (i, n-1)

      if (x(i+1) - x(i) .eq. 0)  stop 'TERP-1'

      y0 = y(i) +  (x0 - x(i)) * (y(i+1) - y(i)) / (x(i+1) - x(i))

      return
      end

      function locat1 (x, n, xx)
      integer  u, m, n
      double precision x
      real xx(n)

c     Binary search for index of grid point immediately below x.
c     Array xx required to be monotonic increasing.
c     Returns
c     0            x <  xx(1)
c     1            x =  xx(1)
c     i            x =  xx(i)
c     n            x >= xx(n)

      locat1 = 0
      u = n+1

   10 if (u-locat1 .gt. 1)  then
         m = (u + locat1) / 2
         if (x .lt. xx(m))  then
            u = m
         else
            locat1 = m
         endif
         goto 10
      endif

      return
      end
