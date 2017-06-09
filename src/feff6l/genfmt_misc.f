      subroutine genfmt_writefeffdat(ipath, ntext, text,
     $     nleg, deg, reff, rnrmav, edge, rat, iz, ipot,
     $     potlbl, l0, il0, ne, xk, ck, ph, cchi)
c
c  write feffNNNN.dat file for NNNN=ipath
      implicit none

      include 'const.h'
      include 'dim.h'
      include 'vers.h'
      double precision eps
      parameter (eps = 1.0d-16)

      integer ipath, ntext, nleg, l0, il0, ne, it, istrln
      integer ipot(0:legtot), iz(0:npotx), jpot, ileg, ie, ios

      double precision deg, reff, rnrmav, edge, getxk
      double precision rat(3, 0:legtot+1)
      double precision xlam, redfac, phff, phffo, cdelt, cdelto
      double precision xk(nex), ckmag(nex)

      complex*16 ph(nex, ltot+1, 0:npotx)
      complex*16 cchi(nex), ck(nex), cfms

      character fname*16
      character*80 text(40)
      character*6 potlbl(0:npotx)

      external istrln, getxk

 50   format (a)
 60   format (1x, a)
 70   format (1x, 79('-'))

      write(fname, 241)  ipath
 241  format ('feff', i4.4, '.dat')
      open (unit=3, file=fname, status='unknown', iostat=ios)
      call chopen (ios, fname, 'genfmt')
c     put header on feff.dat
      do 245  it = 1, ntext
         write(3,60)  text(it)(1:istrln(text(it)))
 245  continue
      write(3,250) ipath, 2, vfeff, vgenfm
 250  format (' Path', i5, '      icalc ', i7, t57, 2a12)
      write(3,70)
      write(3,290)  nleg, deg, reff*bohr, rnrmav, edge*ryd
 290  format (1x, i3, f8.3, f9.4, f10.4, f11.5,
     1     ' nleg, deg, reff, rnrmav(bohr), edge')

c
      write(3,300)
 300  format ('        x         y         z   pot at#')
      jpot = ipot(nleg)
      write(3,310)  rat(1,nleg)*bohr, rat(2,nleg)*bohr,
     $     rat(3,nleg)*bohr, jpot, iz(jpot), potlbl(jpot)
 310  format (1x, 3f10.4, i3, i4, 1x, a6, '   absorbing atom')
 311  format (1x, 3f10.4, i3, i4, 1x, a6)
      do 330  ileg = 1, nleg-1
         jpot = ipot(ileg)
         write(3,311)  rat(1,ileg)*bohr, rat(2,ileg)*bohr,
     $        rat(3,ileg)*bohr, jpot, iz(jpot), potlbl(jpot)
 330  continue


      write(3,340)
 340  format    ('    k   real[2*phc]   mag[feff]  phase[feff]',
     1     ' red factor   lambda      real[p]@#')

c     Make the feff.dat stuff and write it to feff.dat
      do 900  ie = 1, ne
c     Consider chi in the standard XAFS form.  Use R = rtot/2.
         xlam = 1.0e10
         if (abs(dimag(ck(ie))) .gt. eps) xlam = 1/dimag(ck(ie))
         redfac = exp (-2 * dimag (ph(ie,il0,ipot(nleg))))
         cdelt = 2*dble(ph(ie,il0,ipot(nleg)))
         cfms = cchi(ie) * xk(ie) * reff**2 *
     1        exp(2*reff/xlam) / redfac
         if (abs(cchi(ie)) .lt. eps)  then
            phff = 0
         else
            phff = atan2 (dimag(cchi(ie)), dble(cchi(ie)))
         endif
c     remove 2 pi jumps in phases
         if (ie .gt. 1)  then
            call pijump (phff, phffo)
            call pijump (cdelt, cdelto)
         endif
         phffo = phff
         cdelto = cdelt

c     write 1 k, momentum wrt fermi level k=sqrt(p**2-kf**2)
c           2 central atom phase shift (real part),
c           3 magnitude of feff,
c           4 phase of feff,
c           5 absorbing atom reduction factor,
c           6 mean free path = 1/(Im (p))
c           7 real part of local momentum p
         write(3,640)
     1        xk(ie)/bohr,
     2        cdelt + l0*pi,
     3        abs(cfms) * bohr,
     4        phff - cdelt - l0*pi,
     5        redfac,
     6        xlam * bohr,
     7        dble(ck(ie))/bohr
 640     format (1x, f6.3, 1x, 3(1pe11.4,1x),0pe11.4,1x,
     1        2(1pe11.4,1x))
 900  continue

      close(3)
      return
      end

      subroutine calc_zabinsky(ne, ik0, deg, ck, cchi, xportx, crit)
c
c calculate Zabinsky Curved Wave Importance Factor for a Path
c
c  input:  ne    number of energy points
c          ik0
c          deg   degeneracy
c          ck    complex k
c          cchi  complex chi(k)
c          xport max importance factor so far
c  output: xport max importance factor so far
c          crit  importance / max importance

      implicit none
      include 'dim.h'
      complex*16 cchi(nex), ck(nex)
      double precision ckmag(nex), ffmag(nex)
      double precision xport, xportx, deg, crit
      integer ie, ne, ik0, nemax

      do 10  ie = 1, ne
         ckmag(ie) = abs(ck(ie))
         ffmag(ie) = abs(cchi(ie))
 10   continue

c     integrate from edge (ik0) to ne
      nemax = ne - ik0 + 1
      call trap(ckmag(ik0), ffmag(ik0), nemax, xport)
      xport = abs(deg*xport)
      if (xport .gt. xportx)  xportx = xport
      crit = 100 * xport / xportx

      return
      end
