      subroutine outcrt (npat, ipat, ckspc,
     1    nncrit, fbetac, xlamc, ne, ik0, cksp, 
     1    fbeta, xlam, ipot,
     1    xport, xheap, xheapr,
     1    xout, xcalcx)

c     This make pw importance factor for pathsd, also recalculates
c     pathfinder criteria for output.  Pathfinder recalculation
c     is hacked from ccrit, so be sure to update this if ccrit
c     is changed.

      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'
      dimension ipat(npatx)
      dimension ipot(0:natx)
      parameter (necrit=9, nbeta=40)
      dimension fbetac(-nbeta:nbeta,0:nphx,necrit), ckspc(necrit)
      dimension fbeta(-nbeta:nbeta,0:nphx,nex), cksp(nex)
      dimension xlamc(necrit), xlam(nex)

c     local variables
      dimension ri(npatx+1), beta(npatx+1), indbet(npatx+1)
      dimension xporti(nex)
      parameter (eps = 1.0e-6)

c     Space for variables for time reversed path (used in xheapr
c     calculation below)
      dimension ipat0(npatx)
      dimension ri0(npatx+1), indbe0(npatx+1)

c     mrb is 'efficient' way to get only ri and beta
c     note that beta is cos(beta)
      call mrb (npat, ipat, ri, beta)

c     Make index into fbeta array (this is nearest cos(beta) grid point,
c     code is a bit cute [sorry!], see prcrit for grid).
      do 290  i = 1, npat+1
         tmp = abs(beta(i))
         n = int(tmp / 0.025)
         del = tmp - n*0.025
         if (del .gt. 0.0125)  n = n+1
         if (beta(i) .lt. 0)  n = -n
         indbet(i) = n
  290 continue

c     Make pw importance factor by integrating over all points
c     above the edge
c     Path importance factor is integral d|p| of
c        (product of f(beta)/rho for the scatterers) * cos(beta0)/rho0
c     Include mean free path factor, exp(-rtot/xlam)
      rtot = 0
      do 510  i = 1, npat+1
         rtot = rtot + ri(i)
  510 continue
      do 560  ie = ik0, ne
         rho = ri(npat+1) * cksp(ie)
         crit = max (abs(beta(npat+1)), 0.3) / rho
         do 520  iat = 1, npat
            rho = ri(iat) * cksp(ie)
            ipot0 = ipot(ipat(iat))
            crit = crit * fbeta(indbet(iat),ipot0,ie) / rho
  520    continue
         crit = crit * exp (-rtot/xlam(ie))
         xporti(ie) =  abs(crit)
  560 continue


c     integrate from ik0 to ne
      nmax = ne - ik0 + 1
      call strap (cksp(ik0), xporti(ik0), nmax, xport)

c     Stuff for  output.
c     Heap crit thing (see ccrit and mcrith for comments)
c     If a path got time reversed, its xheap may be smaller than
c     it was before it got time-reversed.  So calculate it both
c     ways.
c     xheap for path, xheapr for time-reversed path

      xheap  = -1
      xheapr = -1
      call mcrith (npat, ipat, ri, indbet,
     1             ipot, nncrit, fbetac, ckspc, xheap)

c     Prepare arrays for time reversed path and make xheapr
c     See timrev.f for details on indexing here.

      nleg = npat+1
c     ri
      do 200  i = 1, nleg
         ri0(i) = ri(nleg+1-i)
  200 continue
c     indbet  and ipat
      indbe0(nleg) = indbet(nleg)
      do 210  i = 1, nleg-1
         indbe0(i) = indbet(nleg-i)
         ipat0(i) = ipat(nleg-i)
  210 continue

      call mcrith (npat, ipat0, ri0, indbe0,
     1             ipot, nncrit, fbetac, ckspc, xheapr)

c     Keep crit thing (see mcritk for comments)
      call mcritk (npat, ipat, ri, beta, indbet,
     1             ipot, nncrit, fbetac, xlamc, ckspc, xout, xcalcx)

      return
      end
