      subroutine eu2au
c     generic routine to
c     transform all code variables from exp. units (Angstrom, eV, etc.)
c     to the atomic units (bohrs, hartrees, etc.)

      implicit double precision (a-h, o-z)
      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'

c     sig2g is not transformed
      rmax  = rmax  / bohr
      rfms1 = rfms1 / bohr
      rfms2 = rfms2 / bohr
      rdirec = rdirec / bohr
      vr0   = vr0 / hart
      vi0   = vi0 / hart
      vrcorr = vrcorr / hart 
      vicorr = vicorr / hart 
      gamach = gamach / hart
      ecv   = ecv   / hart
      emin  = emin  / hart
      emax  = emax  / hart
      eimag = eimag / hart
      vixan = vixan / hart
      xkstep = xkstep * bohr
      xkmax  = xkmax  * bohr
      totvol = totvol / bohr**3
      do 10 iat = 1, nat
      do 10 i = 1,3
        rat(i,iat) = rat (i, iat) / bohr
  10  continue
      do 20 iph = 0, nph
      do 20 iovr = 1, novr(iph)
         rovr(iovr,iph) = rovr(iovr,iph) / bohr
  20  continue

      return
      end
