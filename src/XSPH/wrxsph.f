      subroutine wrxsph (phpad,
     1       nsp, ne, ne1, ne3, nph, ihole, rnrmav,xmu,edge,
     2       ik0, ixc, rs, vint,
     3       em, eref, lmax, iz, potlbl, ph, rkk)
      implicit double precision (a-h, o-z)
c     writes down file 'phase.pad' to be read by rphpad
c  Energy grid information
c     em   - complex energy grid
c     eref - V_int + i*gamach/2 + self-energy correction
c     ne   - total number of points in complex energy grid
c     ne1  - number of points on main horizontal axis
c     ne2  - number of points on vertical vertical axis ne2=ne-ne1-ne3
c     ne3  - number of points on auxilary horizontal axis (need for f')
c     xmu  - Fermi energy
c     edge - x-ray frequency for final state at Fermi level
c     ik0  - grid point index at Fermi level
c  Potential type information
c     ixc    - potential model (this, rs, vint added to phase.pad for sake of onepath.f)
c     rs     - r_s estimate from rhoint (4/3 r_s**3 * rhoint = 1)
c     vint   - muffin-tin zero energy (interstitial potential) 
c     nph    - number of potential types
c     iz     - charge of nuclei (atomic number)
c     potlbl - label for each potential type
c     lmax   - max orb momentum for each potential type
c     ihole  - index of core-hole orbital for absorber (iph=0)
c     rnrmav - average Norman radius (used in headers only)
c  Main output of xsect and phases module (except that in xsect.bin)
c     ph  - complex scattering phase shifts
c     rkk - complex multipole matrix elements

      include '../HEADERS/dim.h'

      character*256 phpad

      character*6  potlbl
      dimension  potlbl(0:nphx)

      complex*16 ph(nex,-ltot:ltot,nspx,0:nphx), eref(nex,nspx), em(nex)
      complex*16 rkk(nex, 8, nspx)
      dimension lmax(0:nphx)
      dimension iz(0:nphx)

c     Local staff
c     npadx control padlib precision (see padlib package)
      parameter (npadx=8)
c     use temp to write ph, rkk, since ne < nex
      complex*16 temp(nex*(2*ltot+1))
      dimension dum(3)

      
c     intialize temp to all 0+i0
      do 5 i = 1, nex*(2*ltot+1)
         temp(i) = dcmplx(0.0,0.0)
 5    continue

c      print *, '>', phpad(1:istrln(phpad)), '<'
      open (unit=1, file=phpad, status='unknown', iostat=ios)
      call chopen (ios, 'phase.pad', 'wrxsph')

      write(1,10) nsp, ne, ne1, ne3, nph, ihole, ik0, npadx, ixc,
     &     rs, vint
 10   format (9(1x,i4), 2(1x,f10.5))

      dum(1) = rnrmav
      dum(2) = xmu 
      dum(3) = edge
      call wrpadd(1, npadx, dum(1), 3)

      call wrpadx(1, npadx, em(1), ne)
      ii = 0
      do 60 isp = 1, nsp
      do 60 ie=1, ne
        ii = ii + 1
        temp(ii) = eref (ie, isp)
  60  continue
      call wrpadx (1, npadx, temp(1), ii)

      do 80  iph = 0, nph
         write(1, 20) lmax(iph), iz(iph), potlbl(iph)
  20     format(2(1x,i3), 1x, a6)
         do 75  isp = 1, nsp
            ii = 0
            do 70  ie = 1, ne
            do 70  ll = -lmax(iph), lmax(iph)
               ii = ii+ 1
               temp(ii) = ph(ie, ll, isp, iph)
   70       continue
            call wrpadx (1, npadx, temp(1), ii )
   75    continue
   80 continue

      ii = 0
      do 90 isp = 1, nsp
      do 90 kdif = 1, 8
      do 90 ie=1, ne
        ii = ii + 1
        temp(ii) = rkk (ie, kdif, isp)
  90  continue
      call wrpadx (1, npadx, temp(1), ii)

      close (unit=1)

      return
      end
