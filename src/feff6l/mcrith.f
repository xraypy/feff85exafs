      subroutine mcrith (npat, ipat, ri, indbet,
     1                   ipot, nncrit, fbetac, ckspc, xheap)

      include 'const.h'
      include 'dim.h'
      dimension ipat(npatx)
      dimension ri(npatx+1), indbet(npatx+1)
      dimension ipot(0:natx)
      parameter (necrit=9, nbeta=40)
      dimension fbetac(-nbeta:nbeta,0:npotx,necrit), ckspc(necrit)

c     Decide if we want the path added to the heap.

      if (ipat(npat) .eq. 0 .or. npat.le.2)  then
c        Partial path is used for xheap, not defined for ss and
c        triangles.  Special case: central atom added to end of path 
c        necessary for complete tree, but not a real path, again,
c        xheap not defined.  Return -1 as not-defined flag.
         xheap = -1
      else
c        Calculate xheap and see if we want to add path to heap.
c        Factor for comparison is sum over nncrit of
c        f(beta1)*f(beta2)*..*f(beta npat-2)/(rho1*rho2*..*rho npat-1).
c        Compare this to sum(1/p), multiply by 100 so we can think 
c        in percent.  Allow for degeneracy when setting crit.
         xheap = 0
         spinv = 0
         do 340  icrit = 1, nncrit
            x = ckspc(icrit) ** (-(npat-1)) * ri(npat-1)
            do 320  i = 1, npat-2
               ipot0 = ipot(ipat(i))
               x = x * fbetac(indbet(i),ipot0,icrit) / ri(i)
  320       continue
            spinv = spinv + 1/ckspc(icrit)
            xheap = xheap + x
  340    continue
         xheap = 100 * xheap / spinv

c        Factor for comparison is sum over nncrit of
c        New xheap:
c        Full chi is
c f(beta1)*f(beta2)*..*f(beta npat)cos(beta0)/(rho1*rho2*..*rho nleg).
c Some of this stuff may change when the path is modified --
c we can't use rho nleg or nleg-1, beta0, beta(npat) or beta(npat-1).
c We DO want to normalize wrt first ss path, f(pi)/(rho nn)**2.
c
c So save f(pi)/(rho nn)**2, 
c calculate 
c f(beta1)*f(beta2)*..*f(beta npat-2)/(rho1*rho2*..*rho npat-1).
c divide nn ss term by stuff we left out -- beta(npat), beta(npat-1),
c cos(beta0), rho nleg, rho nleg-1.
c
c Sum this over nncrit and try it out.
*
c        Sum over nncrit of
c        1/(rho1+rho2+..+rho npat-1).
*        reff = 0
*        do 350  i = 1, npat-1
*           reff = reff + ri(i)
* 350    continue
*        xss = 0
*        do 360  icrit = 1, nncrit
*           rho = ckspc(icrit) * reff
*           xss = xss + 1/rho
* 360    continue
*        xheap = 100 * xheap / xss
      endif

      return
      end
