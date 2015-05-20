      subroutine rdpot ( ntitle, title, rnrmav, xmu, vint, rhoint,
     1                  emu, s02, erelax, wp, ecv,rs,xf, qtotel, 
     2                  imt, rmt, inrm, rnrm, folp, folpx, xnatph,
     3                  dgc0, dpc0, dgc, dpc, adgc, adpc,
     3                  edens, vclap, vtot, edenvl, vvalgs, dmag, xnval,
     4                  eorb, kappa, iorb, qnrm, xnmues, nohole, ihole,
     5                  inters, totvol, iafolp, xion, iunf, iz, jumprm)
c  opens pot.pad file and reads following information
c  General:
c     ntitle - number of title lines
c     title  - title itself
c     emu    - edge position (x-ray energy for final state at Fermi level)
c  Muffin-tin geometry
c     rmt    - muffin-tin radii
c     imt    - index of radial grid just below muffin-tin radii
c     rnrm   - Norman radii
c     inrm   - index of radial grid just below Norman radii
c     rnrmav - average Norman radius
c     folp   - overlap parameter for rmt
c     folpx  - maximum value for folp
c     xnatph - number of atoms of each potential type
c  Atomic orbitals info (Dirac spinors)
c     dgc0   - upper component for initial orbital
c     dpc0   - lower component for initial orbital
c     dgc    - upper components for all atomic orbitals
c     dpc    - lower components for all atomic orbitals
c     adgc   - development coefficient for upper components
c     adpc   - development coefficient for lower components
c     xnval  - number of valence electrons for each atomic orbital
c              used for core-valence separation and non-local exchange
c     eorb  - atomic enrgy of each orbital for the absorber
c  Electron density information 
c     rhoint - interstitial density
c     rs     - r_s estimate from rhoint (4/3 r_s**3 * rhoint = 1)
c     xf     - estimate of momentum at Fermi level from rhoint
c     edens  - total electron density
c     edenvl - density from valence electrons
c     dmag   - density for spin-up minus density for spin-down
c     qtotel - total charge of a cluster
c     qnrm   - charge accumulated inside Norman sphere as result of SCF
c     xnmues - occupation numbers of valence orbitals from SCF procedure
c  Potential information
c     xmu    - Fermi level position
c     ecv    - core-valence separation energy
c     vint   - muffin-tin zero energy (interstitial potential)
c     vclap  - Coulomb potential
c     vtot   - vclap + xc potential from edens
c     vvalgs - vclap + xc potential from edenvl (EXCHANGE 5 model)
c  Specific data for convolution with excitation spectrum (see mbconv)
c     s02    - many body reduction factor S_0^2 
c     erelax - estimate of relaxation energy = efrozen - emu, where
c              efrozen is edge position estimate from Koopmans theorem
c     wp     - estimate of plasmon frequency from rhoint

      implicit double precision (a-h, o-z)
      include '../HEADERS/dim.h'
      dimension adgc(10, 30, 0:nphx)
      dimension adpc(10, 30, 0:nphx)
      dimension dgc(251, 30, 0:nphx)
      dimension dgc0(251)
      dimension dmag(251, 0:nphx)
      dimension dpc(251, 30, 0:nphx)
      dimension dpc0(251)
      dimension edens(251, 0:nphx)
      dimension edenvl(251, 0:nphx)
      dimension eorb(30)
      dimension folp(0:nphx)
      dimension folpx(0:nphx)   
      dimension imt(0:nphx)
      dimension inrm(0:nphx)
      dimension iorb(-4:3,0:nphx)
      dimension iz(0:nphx)
      dimension kappa(30)
      dimension qnrm(0:nphx)
      dimension rmt(0:nphx)
      dimension rnrm(0:nphx)
      dimension vclap(251, 0:nphx)
      dimension vtot(251, 0:nphx)
      dimension vvalgs(251, 0:nphx)
      dimension xion(0:nphx)
      dimension xnatph(0:nphx)
      dimension xnmues(0:lx,0:nphx)
      dimension xnval(30,0:nphx)

      character*80 title(nheadx)

      dimension dum(13)

  10  format(a)

      open (unit=3, file='pot.pad', status='old')
      read(3,30) ntitle, nph, npadx, nohole, ihole, inters, iafolp,
     1            jumprm, iunf
  30  format(9(1x,i4))
c     nph and npadx are not passed to calling subroutine
      do 133  i  = 1, ntitle
         read(3,10) title(i)
         call triml(title(i))
  133 continue
c     Misc double precision stuff from pot.pad
      call rdpadd(3, npadx, dum(1), 13)
      rnrmav = dum(1)
      xmu    = dum(2)
      vint   = dum(3)
      rhoint = dum(4)
      emu    = dum(5)
      s02    = dum(6)
      erelax = dum(7)
      wp     = dum(8)
      ecv    = dum(9)
      rs     = dum(10)
      xf     = dum(11)
      qtotel = dum(12)
      totvol = dum(13)

c     read imt
      read (3, 40) (imt(i),i=0,nph)
  40  format(20(1x,i4))
      call rdpadd(3, npadx, rmt(0), nph+1)
c     read inrm
      read (3, 40) (inrm(i),i=0,nph)
      read (3, 40) (iz(i),i=0,nph)
      read (3, 40) (kappa(i),i=1,30)
      call rdpadd(3, npadx, rnrm(0), nph+1)
      call rdpadd(3, npadx, folp(0), nph+1)
      call rdpadd(3, npadx, folpx(0), nph+1)
      call rdpadd(3, npadx, xnatph(0), nph+1)
      call rdpadd(3, npadx, xion(0), nph+1)
      call rdpadd(3, npadx, dgc0(1), 251)
      call rdpadd(3, npadx, dpc0(1), 251)
      call rdpadd(3, npadx, dgc(1,1,0), 251*30*(nph+1) )
      call rdpadd(3, npadx, dpc(1,1,0), 251*30*(nph+1) )
      call rdpadd(3, npadx, adgc(1,1,0), 10*30*(nph+1) )
      call rdpadd(3, npadx, adpc(1,1,0), 10*30*(nph+1) )
      call rdpadd(3, npadx, edens(1,0), 251*(nph+1) )
      call rdpadd(3, npadx, vclap(1,0), 251*(nph+1) )
      call rdpadd(3, npadx, vtot(1,0), 251*(nph+1) )
      call rdpadd(3, npadx, edenvl(1,0), 251*(nph+1) )
      call rdpadd(3, npadx, vvalgs(1,0), 251*(nph+1) )
      call rdpadd(3, npadx, dmag(1,0), 251*(nph+1) )
      call rdpadd(3, npadx, xnval(1,0), 30*(nph+1) )
      call rdpadd(3, npadx, eorb(1), 30)
      do 50 iph=0,nph
 50   read (3, 60) (iorb(i,iph),i=-4,3)
 60   format(8(1x,i2))
      call rdpadd(3, npadx, qnrm(0), nph+1 )
      nn = (lx+1)*(nph+1)
      call rdpadd(3, npadx, xnmues(0,0), nn )
      close (unit=3)

      return
      end
