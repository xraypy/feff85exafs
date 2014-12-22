      subroutine rdxsph ( phpad, 
     1     ne, ne1, ne3, nph, ihole, rnrmav, xmu, edge,
     2     ik0, ixc, rs, vint,
     3     em, eref, iz, potlbl, ph, rkk, lmax, lmaxp1)
      implicit double precision (a-h, o-z)
c     reads file 'phase.pad' 
c
c     phpad - specify path to phase.pad     (character*256)
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
c     ixc    - potential model
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
      complex*16 rkk(nex,8,nspx)
      dimension lmax0(0:nphx), lmax(nex,0:nphx)
      dimension iz(0:nphx)
c     kinit, linit, ilinit,  - initial state kappa and ang. mom.
c     lmaxp1  -largest lmax in problem + 1

c     phmin is min value to use for |phase shift|
      parameter (phmin = 1.0d-7)

c     Local staff
c     use temp to write ph, rkk, since ne < nex
      complex*16 temp(nex*(2*ltot+1))
      dimension dum(3)

      call triml(phpad)
c      print *, istrln(phpad), '--', phpad(1:istrln(phpad)), '--'
      open (unit=1, file=phpad, status='old', iostat=ios, err=3)
      goto 6
 3    continue
      open (unit=1, file='phase.pad', status='old', iostat=ios, err=4)
      goto 5

 4    continue
      stop 'cannot find phase.pad in rdxsph'

 5    continue
      phpad = 'phase.pad'

 6    continue
      call chopen (ios, 'phase.pad', 'rdxsph')

      ixc = 0
      rs = 0.
      vint = 0.
      read(1,10) nsp, ne, ne1, ne3, nph, ihole, ik0, npadx, ixc,
     &     rs, vint
  10  format (9(1x,i4), 2(1x,f10.5))

      call rdpadd(1, npadx, dum(1), 3)
      rnrmav = dum(1)
      xmu    = dum(2)
      edge   = dum(3)

      call rdpadx(1, npadx, em(1), ne)
c     call rdpadx(1, npadx, eref(1), ne)
      call rdpadx (1, npadx, temp(1), ne*nsp)
      ii = 0
      do 60 isp = 1, nsp
      do 60 ie=1, ne
        ii = ii + 1
        eref (ie, isp) = temp(ii)
  60  continue

      do 80  iph = 0, nph
         read(1, 20)  lmax0(iph), iz(iph), potlbl(iph)
  20     format(2(1x,i3), 1x, a6)

         do 75 isp = 1,nsp 
            ii = ne * (2*lmax0(iph)+1)
            call rdpadx (1, npadx, temp(1), ii )
            ii = 0
            do 70  ie = 1, ne
            do 70  ll = -lmax0(iph), lmax0(iph)
               ii = ii+ 1
               ph(ie,ll,isp,iph) = temp(ii)
   70       continue
   75    continue
   80 continue

      call rdpadx (1, npadx, temp(1), ne*8*nsp)
      ii = 0
      do 90 isp = 1,nsp 
      do 90 kdif = 1, 8
      do 90 ie=1, ne
        ii = ii + 1
        rkk (ie, kdif, isp) = temp(ii)
  90  continue

      close (unit=1)

c     make additional data for output
      lmaxp1 = 0
      do 180  iph = 0, nph
      do 180  ie = 1, ne
c        Set lmax to include only non-zero phases
         do 160  il =  lmax0(iph), 0, -1
            lmax(ie,iph) = il
            if (abs(sin(ph(ie, il, 1, iph))) .gt. phmin .or.
     3          abs(sin(ph(ie, il,nsp,iph))) .gt. phmin)  goto 161
  160    continue
  161    continue
         if (lmax(ie,iph)+1 .gt. lmaxp1)  lmaxp1 = lmax(ie,iph)+1
  180 continue

      return
      end
