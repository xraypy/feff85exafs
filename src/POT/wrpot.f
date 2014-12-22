      subroutine wrpot ( nph, ntitle, title, rnrmav, xmu, vint, rhoint,
     1                  emu, s02, erelax, wp, ecv,rs,xf, qtotel, 
     2                  imt, rmt, inrm, rnrm, folp, folpx, xnatph,
     3                  dgc0, dpc0, dgc, dpc, adgc, adpc,
     3                  edens, vclap, vtot, edenvl, vvalgs, dmag, xnval,
     4                  eorb, kappa, iorb, qnrm, xnmues, nohole, ihole,
     5                  inters, totvol, iafolp, xion, iunf, iz, jumprm)
c  opens pot.pad file and writes following information
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
c  Atomic orbitals info (Dirac spinors)
c     dgc0   - upper component for initial orbital
c     dpc0   - lower component for initial orbital
c     dgc    - upper components for all atomic orbitals
c     dpc    - lower components for all atomic orbitals
c     adgc   - development coefficient for upper components
c     adpc   - development coefficient for lower components
c     xnval  - number of valence electrons for each atomic orbital
c     eorb   - atomic energy of occupied orbitals
c     kappa  - quantum number kappa of occupied orbitals
c     iorb   - index of last occupied orbital for each kappa
c              used for core-valence separation and non-local exchange
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
      parameter (npadx=8)
      dimension imt(0:nphx), rmt(0:nphx), inrm(0:nphx),  rnrm(0:nphx)
      dimension folp(0:nphx), folpx(0:nphx), dgc0(251), dpc0(251)
      dimension dgc(251, 30, 0:nphx), dpc(251, 30, 0:nphx)
      dimension adgc(10, 30, 0:nphx), adpc(10, 30, 0:nphx)
      dimension edens(251, 0:nphx), vclap(251, 0:nphx)
      dimension vtot(251, 0:nphx), edenvl(251, 0:nphx)
      dimension vvalgs(251, 0:nphx), dmag(251, 0:nphx)
      dimension xnval(30,0:nphx), qnrm(0:nphx), xnmues(0:lx,0:nphx)
      dimension eorb(30), kappa(30)
      dimension iorb(-4:3,0:nphx), iz(0:nphx), xion(0:nphx)
      dimension xnatph(0:nphx)
      character*80 title(nheadx)

      dimension dum(13)
      character*75 wfmt

  10  format(a)
      write (wfmt, 35) nph + 1
  35  format( '(', i3, '(1x,i4))' )

      open (unit=3, file='pot.pad', status='unknown', iostat=ios)
      call chopen (ios, 'pot.pad', 'pot')
      write(3,20) ntitle, nph, npadx, nohole, ihole, inters, iafolp,
     1            jumprm, iunf
  20  format (9(1x,i4))
      do 133  i  = 1, ntitle
         ll = istrln(title(i))
         write(3,10) title(i)(1:ll)
  133 continue
c     Misc stuff from pot.pad
      dum(1)  = rnrmav
      dum(2)  = xmu
      dum(3)  = vint
      dum(4)  = rhoint
      dum(5)  = emu
      dum(6)  = s02
      dum(7)  = erelax
      dum(8)  = wp
      dum(9)  = ecv
      dum(10)  = rs
      dum(11)  = xf
      dum(12)  = qtotel
      dum(13)  = totvol
      call wrpadd(3, npadx, dum(1), 13)

      write (3, 40) (imt(i),i=0,nph)
  40  format(20(1x,i4))
      
      call wrpadd(3, npadx, rmt(0), nph+1)

      write (3, 40) (inrm(i),i=0,nph)
      write (3, 40) (iz(i),i=0,nph)
      write (3, 40) (kappa(i),i=1,30)
      
      call wrpadd(3, npadx, rnrm(0), nph+1)
      call wrpadd(3, npadx, folp(0), nph+1)
      call wrpadd(3, npadx, folpx(0), nph+1)
      call wrpadd(3, npadx, xnatph(0), nph+1)
      call wrpadd(3, npadx, xion(0), nph+1)
      call wrpadd(3, npadx, dgc0(1), 251)
      call wrpadd(3, npadx, dpc0(1), 251)
      call wrpadd(3, npadx, dgc(1,1,0), 251*30*(nph+1) )
      call wrpadd(3, npadx, dpc(1,1,0), 251*30*(nph+1) )
      call wrpadd(3, npadx, adgc(1,1,0), 10*30*(nph+1) )
      call wrpadd(3, npadx, adpc(1,1,0), 10*30*(nph+1) )
      call wrpadd(3, npadx, edens(1,0), 251*(nph+1) )
      call wrpadd(3, npadx, vclap(1,0), 251*(nph+1) )
      call wrpadd(3, npadx, vtot(1,0), 251*(nph+1) )
      call wrpadd(3, npadx, edenvl(1,0), 251*(nph+1) )
      call wrpadd(3, npadx, vvalgs(1,0), 251*(nph+1) )
      call wrpadd(3, npadx, dmag(1,0), 251*(nph+1) )
      call wrpadd(3, npadx, xnval(1,0), 30*(nph+1) )
      call wrpadd(3, npadx, eorb(1), 30)
      
      do 50 iph=0,nph
        write (3, 60) (iorb(i,iph),i=-4,3)
 50   continue
 60   format(8(1x,i2))
      call wrpadd(3, npadx, qnrm(0), nph+1 )
      call wrpadd(3, npadx, xnmues(0,0), (lx+1)*(nph+1) )
      close (unit=3)

      return
      end
