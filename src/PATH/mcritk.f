      subroutine mcritk (npat, ipat, ri, beta, indbet,
     1      ipot, nncrit, fbetac, xlamc, ckspc, xout, xcalcx)

      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'
      dimension ipat(npatx)
      dimension ri(npatx+1), beta(npatx+1), indbet(npatx+1)
      dimension ipot(0:natx)
      parameter (necrit=9, nbeta=40)
      dimension fbetac(-nbeta:nbeta,0:nphx,necrit), ckspc(necrit)
      dimension xlamc(necrit)

cc    xcalcx is max xcalc encountered so far.  Set to -1 to reset it --
cc    otherwise it gets passed in and out as mcritk gets called.
c     calculation of xcalcx changed by ala. It is calculated only
c     a first call, i.e. for the NN SS path, and is not recalculated

c     We may want path in heap so that other paths built from this
c     path will be considered, but do not want this path to be
c     written out for itself.  Decide that now and save the flag
c     in the heap, so we won't have to re-calculate the mpprm
c     path parameters later.

c     Do not want it for output if last atom is central atom,
c     use xout = -1 as flag for undefined, don't keep it.
      if (ipat(npat) .eq. 0)  then
         xout = -1
         return
      endif

c     Make xout, output inportance factor.  This is sum over p of
c     (product of f(beta)/rho for the scatterers) * 
c                                 (cos(beta0)/rho(npat+1).
c     Compare this to xoutx, max xout encountered so far.
c     Use mean free path factor, exp(-rtot/xlam)
c     Multiply by 100 so we can think in percent.

      xcalc = 0
      rtot = 0
      do 410  i = 1, npat+1
         rtot = rtot + ri(i)
  410 continue
      do 460  icrit = 1, nncrit
         rho = ri(npat+1) * ckspc(icrit)
c        when beta(0)=90 degrees, get zero, so fudge with cos=.2
         x = max (abs(beta(npat+1)), 0.3) / rho
         do 420  iat = 1, npat
            rho = ri(iat) * ckspc(icrit)
            ipot0 = ipot(ipat(iat))
            x = x * fbetac(indbet(iat),ipot0,icrit) / rho
  420    continue
         x = x * exp (-rtot/xlamc(icrit))
         xcalc = xcalc + x
  460 continue
      if (xcalcx.le.0)  xcalcx = xcalc
      xout = 100 * xcalc / xcalcx
      return
      end
