      subroutine genfmt (ipr3, critcw, sig2g, iorder)
      implicit double precision (a-h, o-z)

      include 'const.h'
      include 'dim.h'
      include 'clmz.h'
      include 'fmatrx.h'
      include 'lambda.h'
      include 'pdata.h'
      include 'nlm.h'
      include 'rotmat.h'

      include 'vers.h'
      include 'pola.h'

      complex*16  rho(legtot), pmati(lamtot,lamtot,2)
      complex*16  pllp, ptrac, srho, prho, cdel1, cfac
      complex*16  cchi(nex), cfms, mmati
      dimension   mmati(-mtot:mtot,-mtot:mtot)
      dimension   t3j(-mtot-1:mtot+1,-1:1)
      dimension   xk(nex), ckmag(nex)
      complex*16  ck(nex)
      dimension   ffmag(nex)

      character*12 fname, messag*128

      logical done

c     Input flags:
c     iorder, order of approx in f-matrix expansion (see setlam)
c             (normal use, 2.  Do ss exactly regardless of iorder)

c     used for divide-by-zero and trig tests
      parameter (eps = 1.0e-16)

c     Read phase calculation input, data returned via commons
      open (unit=1, file='phase.bin', status='old',
     1      access='sequential', form='unformatted', iostat=ios)
      call chopen (ios, 'phase.bin', 'genfmt')
      call rphbin (1)
      close (unit=1)

c     Open path input file (unit in) and read title.  Use unit 1.
      ntitle = 5
      open (unit=1, file='paths.dat', status='old', iostat=ios)
      call chopen (ios, 'paths.dat', 'genfmt')
      call rdhead (1, ntitle, title, ltitle)
      if (ntitle .le. 0)  then
         title(1) = ' '
         ltitle(1) = 1
      endif

c     cgam = gamma in mean free path calc (eV).  Set to zero in this
c     version.  Set it to whatever you want if you need it.
c     cgam = 0
c     cgam = cgam / ryd
c     add cnst imag part to eref
c     do 20  ie = 1, ne
c        eref(ie) = eref(ie) - coni*cgam/2
c  20 continue

   50 format (a)
   60 format (1x, a)
   70 format (1x, 79('-'))

c     Save filenames of feff.dat files for use by ff2chi
      open (unit=2, file='files.dat', status='unknown', iostat=ios)
      call chopen (ios, 'files.dat', 'genfmt')
c     Put phase header on top of files.dat
      do 100  itext = 1, ntext
         write(2,60)  text(itext)(1:ltext(itext))
  100 continue
      write(2,70)
      write(2,120)
  120 format ('    file        sig2   amp ratio    ',
     1        'deg    nlegs  r effective')

c     Set crit0 for keeping feff.dat's
      if (ipr3 .le. 0)  crit0 = 2*critcw/3
c     Make a header for the running messages.
       write(messag, 130) critcw
       call echo(messag)
  130 format ('    Curved wave chi amplitude ratio', f7.2, '%')
      if (ipr3 .le. 0)  then
         write(messag, 131) crit0
         call echo(messag)
       endif
  131 format ('    Discard feff.dat for paths with cw ratio <',
     1         f7.2, '%')
       write(messag, 132)
       call echo(messag)
  132 format ('    path  cw ratio     deg    nleg  reff')

c     Set nlm factors in common /nlm/ for use later
      call snlm (ltot+1, mtot+1)

      if (pola) then
c        Make 3j factors in t3j  (multiplied by sqrt(3*(2l0+1)) for
c        further convinience - the same expression for chi)
c        l0 - final momentum, initial momentum = l0-1.
         do 140  m0 = -l0+1,l0-1
            t3j(m0, 1) = (-1)**(l0+1+m0)*sqrt(3.0d0*(l0+m0)*(l0+m0+1)
     1                /(2*l0)/(2*l0-1))
            t3j(m0, 0) = (-1)**(l0+m0)*sqrt(3.0d0*(l0*l0-m0*m0)/
     1                l0/(2*l0-1))
  140    continue
         do 145  m0 = -l0+1,l0-1
            t3j(m0,-1) = t3j(-m0,1)
  145    continue
      endif

c     While not done, read path, find feff.
      open (unit=4,file='nstar.dat', status='unknown', iostat=ios)
      write(4,198, iostat=ios) evec
  198 format('polarization  ',3f8.4)
      write(4,199, iostat=ios)
  199 format('npath  nstar')
      npath = 0
      ntotal = 0
      nused = 0
      xportx = -1
  200 continue

c        Read current path
         call rdpath (1, pola, done,xstar)
         icalc = iorder
         if (done)  goto  1000
         npath = npath + 1
         ntotal = ntotal + 1

         write (4,201,iostat=ios) npath, xstar
  201    format (i5, f8.4)

c        Need reff
         reff = 0
         do 220  i = 1, nleg
            reff = reff + ri(i)
  220    continue
         reff = reff/2

c        Set lambda for low k
         call setlam (icalc, 1)

c        Calculate and store rotation matrix elements
c        Only need to go to (il0, il0, ...) for isc=nleg and
c        nleg+1 (these are the paths that involve the 'z' atom
         call rot3i (il0, il0, nleg)
         do 400  isc = 1, nsc
            call rot3i (lmaxp1, mmaxp1, isc)
  400    continue
         if (pola) then
c           one more rotation in polarization case
            call rot3i (il0, il0, nleg+1)
            call mmtr(t3j,mmati)
         endif 


c        Big energy loop
         do 800  ie = 1, ne

c           real momentum (k)
            xk(ie) = getxk (em(ie) - edge)

c           complex momentum (p)
            ck(ie) = sqrt (em(ie) - eref(ie))
            ckmag(ie) = abs(ck(ie))
c           complex rho
            do 420  ileg = 1, nleg
               rho(ileg) = ck(ie) * ri(ileg)
  420       continue

c           if ck is zero, xafs is undefined.  Make it zero and jump
c           to end of calc part of loop.
            if (abs(ck(ie)) .le. eps)  then
               cchi(ie) = 0
               goto 620
            endif

c           Calculate and store spherical wave factors c_l^(m)z^m/m!
c           in a matrix clmi(il,im,ileg), ileg=1...nleg.
c           Result is that common /clmz/ is updated for use by fmtrxi.

c           zero clmi arrays
            do 440  ileg = 1, legtot
            do 440  il = 1, ltot+1
            do 440  im = 1, mtot+ntot+1
               clmi(il,im,ileg) = 0
  440       continue

            mnmxp1 = mmaxp1 + nmax

            lxp1 = max (lmax(ie,ipot(1))+1, l0+1)
            mnp1 = min (lxp1, mnmxp1)
            call sclmz (rho, lxp1, mnp1, 1)

            lxp1 = max (lmax(ie,ipot(nsc))+1, l0+1)
            mnp1 = min (lxp1, mnmxp1)
            call sclmz (rho, lxp1, mnp1, nleg)

            do 460  ileg = 2, nleg-1
               isc0 = ileg-1
               isc1 = ileg
               lxp1 = max (lmax(ie,ipot(isc0))+1, lmax(ie,ipot(isc1))+1)
               mnp1 = min (lxp1, mnmxp1)
               call sclmz (rho, lxp1, mnp1, ileg)
  460       continue

c           Calculate and store scattering matrices fmati.

            if (pola) then
c              Polarization version, make new m matrix
c              this will fill fmati(...,nleg) in common /fmtrxi/
               call mmtrxi (laml0x, mmati, ie, 1, nleg)
            else 
c              Termination matrix, fmati(...,nleg)
               iterm = 1
               call fmtrxi (laml0x, laml0x, ie, iterm, 1, nleg)
            endif

            iterm = -1
c           First matrix
            call fmtrxi (lamx, laml0x, ie, iterm, 2, 1)
c           Last matrix if needed
            if (nleg .gt. 2)  then
               call fmtrxi (laml0x, lamx, ie, iterm, nleg, nleg-1)
            endif
c           Intermediate scattering matrices
            do 480  ilegp = 2, nsc-1
               ileg = ilegp + 1
               call fmtrxi (lamx, lamx, ie, iterm, ileg, ilegp)
  480       continue

c           Big matrix multiplication loops.
c           Calculates trace of matrix product
c           M(1,N) * f(N,N-1) * ... * f(3,2) * f(2,1), as in reference.
c           We will (equivalently) calculate the trace over lambda_N of
c           f(N,N-1) * ... * f(3,2) * f(2,1) * M(1,N), working from
c           right to left.
c           Use only 2 pmati arrays, alternating indp (index p)
c           1 and 2.

c           f(2,1) * M(1,N) -> pmat(1)
            indp = 1
            do 520  lmp = 1, laml0x
            do 520  lm = 1, lamx
               pllp = 0
               do 500  lmi = 1, laml0x
                  pllp = pllp + fmati(lm,lmi,1) * fmati(lmi,lmp,nleg)
  500          continue
               pmati(lm,lmp,indp)=pllp
  520       continue

c           f(N,N-1) * ... * f(3,2) * [f(2,1) * M(1,N)]
c           Term in [] is pmat(1)
            do 560 isc = 2, nleg-1
c              indp is current p matrix, indp0 is previous p matrix
               indp = 2 - mod(isc,2)
               indp0 = 1 + mod(indp,2)
               do 550  lmp = 1, laml0x
               do 550  lm = 1, lamx
                  pllp=0
                  do 540 lmi = 1, lamx
                     pllp = pllp +
     1                      fmati(lm,lmi,isc)*pmati(lmi,lmp,indp0)
  540             continue
  550          pmati(lm,lmp,indp) = pllp
  560       continue

c           Final trace over matrix
            ptrac=0
            do 580  lm = 1, laml0x
               ptrac = ptrac + pmati(lm,lm,indp)
  580       continue

c           Calculate xafs
c           srho=sum pr(i), prho = prod pr(i)
            srho=0
            prho=1
            do 600  ileg = 1, nleg
               srho = srho + rho(ileg)
               prho = prho * rho(ileg)
  600       continue
c           Complex chi (without 2kr term)
c           ipot(nleg) is central atom
            cdel1 = exp(2*coni*ph(ie,il0,ipot(nleg)))
            cfac = cdel1 * exp(coni*(srho-2*xk(ie)*reff)) / prho

            cchi(ie) = ptrac * cfac/(2*l0+1)

c           When ck(ie)=0, xafs is set to zero.  Calc above undefined.
c           Jump to here from ck(ie)=0 test above.
  620       continue

c        end of energy loop
  800    continue

c        Make importance factor, deg*(integral (|chi|*d|p|))
c        make ffmag (|chi|)
c        xport   importance factor
         do 810  ie = 1, ne
               ffmag(ie) = abs(cchi(ie))
  810    continue

c        integrate from edge (ik0) to ne
         nemax = ne - ik0 + 1
         call trap (ckmag(ik0), ffmag(ik0), nemax, xport)
         xport = abs(deg*xport)
         if (xport .gt. xportx)  xportx = xport
         crit = 100 * xport / xportx

c        Write output if path is important enough (ie, path is

c        Write feff.dat if we need it.
         if (ipr3 .ge. 1  .or.  crit .ge. crit0)  then
c           Prepare output file feffnnnn.dat (unit 3)
            write(fname,241)  ipath
  241       format ('feff', i4.4, '.dat')
            open (unit=3, file=fname, status='unknown', iostat=ios)
            call chopen (ios, fname, 'genfmt')
c           put header on feff.dat
            do 245  itext = 1, ntext
               write(3,60)  text(itext)(1:ltext(itext))
  245       continue
            write(3,250) ipath, icalc, vfeff, vgenfm
  250       format (' Path', i5, '      icalc ', i7, t57, 2a12)
            write(3,70)
            write(3,290)  nleg, deg, reff*bohr, rnrmav, edge*ryd
  290       format (1x, i3, f8.3, f9.4, f10.4, f11.5, 
     1              ' nleg, deg, reff, rnrmav(bohr), edge')
            write(3,300)
  300       format ('        x         y         z   pot at#')
            write(3,310)  (rat(j,nleg)*bohr,j=1,3), ipot(nleg),
     1                    iz(ipot(nleg)), potlbl(ipot(nleg))
  310       format (1x, 3f10.4, i3, i4, 1x, a6, '   absorbing atom')
            do 330  ileg = 1, nleg-1
               write(3,320)  (rat(j,ileg)*bohr,j=1,3), ipot(ileg),
     1                       iz(ipot(ileg)), potlbl(ipot(ileg))
  320          format (1x, 3f10.4, i3, i4, 1x, a6)
  330       continue

            write(3,340)
  340       format    ('    k   real[2*phc]   mag[feff]  phase[feff]',
     1                 ' red factor   lambda      real[p]@#')

c           Make the feff.dat stuff and write it to feff.dat
            do 900  ie = 1, ne
c              Consider chi in the standard XAFS form.  Use R = rtot/2.
               xlam = 1.0e10
               if (abs(dimag(ck(ie))) .gt. eps) xlam = 1/dimag(ck(ie))
               redfac = exp (-2 * dimag (ph(ie,il0,ipot(nleg))))
               cdelt = 2*dble(ph(ie,il0,ipot(nleg)))
               cfms = cchi(ie) * xk(ie) * reff**2 *
     1              exp(2*reff/xlam) / redfac
               if (abs(cchi(ie)) .lt. eps)  then
                  phff = 0
               else
                  phff = atan2 (dimag(cchi(ie)), dble(cchi(ie)))
               endif
c              remove 2 pi jumps in phases
               if (ie .gt. 1)  then
                  call pijump (phff, phffo)
                  call pijump (cdelt, cdelto)
               endif
               phffo = phff
               cdelto = cdelt

c              write 1 k, momentum wrt fermi level k=sqrt(p**2-kf**2)
c                    2 central atom phase shift (real part),
c                    3 magnitude of feff,
c                    4 phase of feff,
c                    5 absorbing atom reduction factor,
c                    6 mean free path = 1/(Im (p))
c                    7 real part of local momentum p

               write(3,640)
     1            xk(ie)/bohr,
     2            cdelt + l0*pi,
     3            abs(cfms) * bohr,
     4            phff - cdelt - l0*pi,
     5            redfac,
     6            xlam * bohr,
     7            dble(ck(ie))/bohr
  640          format (1x, f6.3, 1x, 3(1pe11.4,1x),0pe11.4,1x,
     1                               2(1pe11.4,1x))
  900       continue

c           Done with feff.dat
            close (unit=3)

c           Put feff.dat and stuff into files.dat
            write(2,820) fname, sig2g, crit, deg,
     1                   nleg, reff*bohr
  820       format(1x, a, f8.5, 2f10.3, i6, f9.4)

c           Tell user about the path we just did
            write(messag, 210) ipath, crit, deg, nleg, reff*bohr
            call echo(messag)            
  210       format (3x, i4, 2f10.3, i6, f9.4)
            nused = nused+1

         else
c           path unimportant, tell user
            write(messag, 211) ipath, crit, deg, nleg, reff*bohr
            call echo(messag)
  211       format (3x, i4, 2f10.3, i6, f9.4, ' neglected')
         endif

c        Do next path
         goto 200

c     Done with loop over paths
 1000 continue
c     close paths.dat, files.dat
      close (unit=1)
      close (unit=2)
      close (unit=4)
       write(messag,1010) nused, ntotal
       call echo(messag)
 1010 format (1x, i4, ' paths kept, ', i4, ' examined.')

      return
      end
