c///////////////////////////////////////////////////////////////////////
c FEFF PROGRAMS (referred below as a System)
c Copyright (c) 1986-2002, University of Washington.
c 
c END-USER LICENSE 
c 
c A signed End-user License Agreement from the University of Washington
c Office of Technology Transfer is required to use these programs and
c subroutines.
c 
c See the URL: http://leonardo.phys.washington.edu/feff/
c 
c USE RESTRICTIONS:
c 
c 1. The End-user agrees that neither the System, nor any of its
c components shall be used as the basis of a commercial product, and
c that the System shall not be rewritten or otherwise adapted to
c circumvent the need for obtaining additional license rights.
c Components of the System subject to other license agreements are
c excluded from this restriction.
c
c 2. Modification of the System is permitted, e.g., to facilitate
c its performance by the End-user. Use of the System or any of its
c components for any purpose other than that specified in this Agreement
c requires prior approval in writing from the University of Washington.
c
c 3. The license granted hereunder and the licensed System may not be
c assigned, sublicensed, or otherwise transferred by the End-user.  
c
c 4. The End-user shall take reasonable precautions to ensure that
c neither the System nor its components are copied, or transferred out
c side of his/her current academic or government affiliated laboratory
c or disclosed to parties other than the End-user.
c 
c 5. In no event shall the End-user install or provide this System
c on any computer system on which the End-user purchases or sells
c computer-related services.
c 
c 6. Nothing in this agreement shall be construed as conferring rights
c to use in advertising, publicity, or otherwise any trademark or the
c names of the System or the UW.   In published accounts of the use or
c application of FEFF the System should be referred to  by this name,
c with an appropriate literature reference:
c 
c FEFF8: A.L. Ankudinov, B. Ravel, J.J. Rehr, and S.D. Conradson,
c        Phys. Rev. B 58, pp. 7565-7576 (1998).
c
c LIMITATION OF LIABILITY:
c
c 1.   THE UW MAKES NO WARRANTIES , EITHER EXPRESSED OR IMPLIED, AS TO
c THE CONDITION OF THE SYSTEM, ITS MERCHANTABILITY, OR ITS FITNESS FOR
c ANY PARTICULAR PURPOSE.  THE END-USER AGREES TO ACCEPT THE SYSTEM
c 'AS IS' AND IT IS UNDERSTOOD THAT THE UW IS NOT OBLIGATED TO PROVIDE
c MAINTENANCE, IMPROVEMENTS, DEBUGGING OR SUPPORT OF ANY KIND.
c
c 2. THE UW SHALL NOT BE LIABLE FOR ANY DIRECT, INDIRECT, SPECIAL,
c INCIDENTAL OR CONSEQUENTIAL DAMAGES SUFFERED BY THE END-USER OR ANY
c OTHER PARTIES FROM THE USE OF THE SYSTEM.
c
c 3.  The End-user agrees to indemnify the UW for liability resulting
c from the use of the System by End-user. The End-user and the UW each
c agree to hold the other harmless for their own negligence.
c
c TITLE:
c
c 1.  Title patent, copyright and trademark rights to the System are
c retained by the UW. The End-user shall take all reasonable precautions
c to preserve these rights.
c 
c 2.  The UW reserves the right to license or grant any other rights to
c the System to other persons or entities.
c///////////////////////////////////////////////////////////////////////
c License is applicable for routines below, until otherwise specified.
c
c     sub-program exchange
      program ffmod6
c     subroutine ffmod6 (iabs)

c     final calculations for various spectroscopies
c     (EXAFS, XANES, FPRIME, DANES, XES)
c     written by a.ankudinov 2000
c     modified by a.ankudinov 2001 for new I/O structure

c     INPUT: mod6.inp global.dat xsect.bin fms.bin list.dat and feff.bin
c     OUTPUT: xmu.dat (chi.dat for EXAFS)

      implicit double precision (a-h, o-z)

c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}
c={../HEADERS/parallel.h
      integer par_type, this_process, numprocs, my_rank
      logical master, worker, parallel_run
      real*8 wall_comm, time_comm
      common /timing/ wall_comm, time_comm
      common /parallel/ numprocs, my_rank, this_process, 
     .          master, worker, parallel_run, par_type
c= ../HEADERS/parallel.h}
      integer elnes,ipmin,ipmax,ipstep  !KJ my variables 1-06 
      integer absolu !KJ 3-06     


      call par_begin
      if (worker) go to 400

c     sub-program exchange
      iabs = 1
c     

c     open the log file, unit 11.  See subroutine wlog.
      open (unit=11, file='log6.dat', status='unknown', iostat=ios)
      call chopen (ios, 'log6.dat', 'feff')

c     read  input files
      call reff2x(mchi, ispec, ipr6, idwopt, critcw, s02, sig2g,
     1                   tk, thetad, mbconv, absolu,  !KJ added absolu 3-06
     1                   vrcorr, vicorr, alphat, thetae, iabs, nabs,
     4            elnes,ipmin,ipmax,ipstep)   !KJ added this line 1-06     

      if (mchi .eq. 1)  then
         call wlog(' Calculating chi...')
         if (ispec.gt.0 .and. ispec.lt.3) then 
c           using FMS+Paths method
            call ff2xmu (ispec, ipr6, idwopt, critcw, s02, sig2g, 
     1                   tk, thetad, mbconv, absolu,  !KJ added absolu 3-06
     1                   vrcorr, vicorr, alphat, thetae, iabs, nabs,
     4            elnes,ipmin,ipmax,ipstep)   !KJ added this line  1-06     
         elseif (ispec.eq.3 .or. ispec.eq.4) then 
c           using FMS+Paths method
            call ff2afs ( ipr6, idwopt, critcw, s02, sig2g, 
     1                   tk, thetad, mbconv, absolu,  !KJ added absolu 4-06
     1                   vrcorr, vicorr, alphat, thetae, iabs, nabs,
     4            elnes,ipmin,ipmax,ipstep)   !KJ added this line 4-06     
         else
c           using MS Paths expansion
            call ff2chi (ispec, ipr6, idwopt, critcw, s02, sig2g, 
     1                   tk, thetad, mbconv, absolu,  !KJ added absolu 3-06
     1                   vrcorr, vicorr, alphat, thetae, iabs, nabs,
     4            elnes,ipmin,ipmax,ipstep)   !KJ added this line 1-06     
         endif
         call wlog(' Done with module 6: DW + final sum over paths.')
      endif
      close (unit=11)

  400 call par_barrier
      call par_end
      
c     sub-program exchange
      stop  
c     return

      end
      subroutine exconv (omega, nk, efermi, s02, erelax, wp, xmu)
      parameter (nfinex = 601)
c     convolution of xmu(e) with excitation spectrum, which is modeled
c     by: f(e) = s02*delta(e) + theta(e)*exp(-e/ed)*x1 + fp(e)
c     plasmon contribution modeled by fp(e) = theta(e-wp)*exp(-e/wp)*x2
c     normalization factors x1, x2 and distribution parameter ed are
c     found from conditions: 1) integral d(e)*f(e) = 1 
c     2) integral d(e)*fp(e) = wwp  0<=wwp<1  s02+wwp<=1
c     3) integral d(e)*e*f(e) = erelax
c     Input:
c       omega - enrgy grid (e)
c       nk    - number of points in energy grid
c       efermi- fermi level
c       s02   - overlap with fully relaxed orbitals
c       erelax- relaxation energy 
c       wp    - plasmon frequency
c       xmu   - original absorption coefficient
c     Output
c       xmu  - result of convolution, rewritten at the end. 
c     This subroutine uses the fact, that if convolution is made for
c     e(i), then for e(i+1), the convolution integral with exp(-e/ed)
c     for e<e(i) is simply scaled by exp((e(i)-e(i+1)) / ed). This makes
c     this convolution fast.
c     written by ala. december 1995
c     
      implicit double precision (a-h,o-z)
      dimension  omega(nk), xmu(nk)
c     work space
      dimension  slope(nfinex), dmu(nfinex), xmup(nfinex)

      if (s02 .ge. 0.999) return
      if (wp .le. 0.0) wp = 0.00001
      if (nk .gt. nfinex) then
         call par_stop('check nfinex in subroutine exconv')
      endif
c     change weight for plasmon excitation here
      wwp = 0.00
c     sm1 - weight for shake up (off) processes
      sm1 = 1.0 - s02 - wwp
      edp = wp
      ed = (erelax - wwp * (wp + edp)) / sm1
      i0 = locat (efermi, nk, omega)
      do 10 i = 1, i0
         slope(i) = 0.0
         dmu(i) = 0.0
 10   continue
      do 20 i = i0, nk - 1
 20     slope(i) = ed * (xmu(i+1) - xmu(i)) / (omega(i+1) - omega(i))
      call terp (omega, xmu, nk, 1, efermi, xmuf)

c     start induction
      xmult = exp ((efermi - omega(i0+1)) / ed)
      dmu(i0+1) = xmu(i0+1) - slope(i0) - xmult * (xmuf - slope(i0))
      do 50 i = i0 + 1, nk - 1
         xmult = exp ((omega(i) - omega(i+1)) / ed)
         dmu(i+1) = xmu(i+1) - slope(i) + xmult*(dmu(i)-xmu(i)+slope(i))
 50   continue
      do 55 i = 1, nk
 55   xmup(i) = s02 * xmu(i) + sm1 * dmu(i)

c     do convolution with plasmon pole
      do 60 i = i0, nk - 1
 60     slope(i) = slope(i) / ed * edp
      xmult = exp ((efermi - omega(i0+1)) / edp)
      dmu(i0+1) = xmu(i0+1) - slope(i0) - xmult * (xmuf - slope(i0))
      do 70 i = i0 + 1, nk - 1
         xmult = exp ((omega(i) - omega(i+1)) / edp)
         dmu(i+1) = xmu(i+1) - slope(i) + xmult*(dmu(i)-xmu(i)+slope(i))
 70   continue

      do 90 i = 1, nk
        en = omega(i) - wp
        j0 = locat(en, nk, omega)
        if (en .gt. efermi) then
           xmult = exp ((omega(j0) - en) / edp)
           dif = xmu(j0) - slope(j0)
           xmup(i) = xmup(i) + wwp * (xmult * (dmu(j0) - dif) + dif +
     1                            slope(j0) * (en - omega(j0)) / edp)
        endif
 90   continue

      do 200 i = 1, nk
 200  xmu(i) = xmup(i)

      return
      end
       subroutine feffdt(ntotal,iplst,nptot,ntext,text,ne,npot,
     $      ihole, iorder, l0, rnrmav, xmu, edge, potlbl,
     $      iz,phc,ck,xk,index,
     $      nleg,deg,nepts,reff,crit,ipot,rat,achi,phchi)
c
c     writes feffnnnn.dat files and files.dat 
c     for compatibility with the old feff
c
      implicit double precision (a-h, o-z)

c={../HEADERS/const.h
      parameter (pi = 3.14159 26535 89793 23846 26433d0)
      parameter (one = 1, zero = 0)
      parameter (third = one/3)
      parameter (raddeg = 180 / pi)
      complex*16 coni
      parameter (coni = (0,1))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919 158 292 677 512 811d0)

      parameter (bohr = 0.529 177 249d0, ryd  = 13.605 698d0)
      parameter (hart = 2 * ryd)
      parameter (alpinv = 137.035 989 56d0)
c     fine structure alpha
      parameter (alphfs = 1 / alpinv)
c= ../HEADERS/const.h}
c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}
c={../HEADERS/vers.h
      character*12 vfeff
c                       123456789012  
      parameter (vfeff='Feff 8.50   ')
c= ../HEADERS/vers.h}
      parameter (eps4 = 1.0e-4)
      parameter (eps = 1.0e-16)

      parameter (npx=15000)
      character*12 fname(npx)
      character*512 slog
      dimension iplst(npx)

c     Stuff from feff.bin, note that floating point numbers are
c     single precision
cc      character*78 string
      real rnrmav, xmu, edge
cc      dimension ltext(nheadx)
      character*80 text(nheadx)
      character*6  potlbl(0:nphx)
      dimension iz(0:nphx)
c     central atom phase shift at l0
      complex phc(nex)
      complex ck(nex)
      real xk(nex)
      dimension index(npx)
      dimension nleg(npx)
      real deg(npx), reff(npx), crit(npx)
      dimension ipot(legtot,npx)
      real rat(3,legtot,npx)
cc      real beta(legtot,npx)
cc      real eta(legtot,npx)
cc      real ri(legtot,npx)
      real achi(nex,npx), phchi(nex,npx)
       integer istrln
       complex*16 cchi, cfms
       external istrln

       call wlog (' feffdt, feff.bin to feff.dat conversion ' // vfeff)

c     read feff.bin
c     Use single precision for all fp numbers in feff.bin
      do 20  itext = 1, ntext
         ltxt = istrln(text(itext))
c        text(itext) does not have carriage control
         call wlog (' ' // text(itext)(1:ltxt))
   20 continue

      write(slog,60)  nptot
   60 format (i8, ' paths to process')
      call wlog (slog)

c     make files.dat
  150 format (a)
  160 format (1x, a)
  170 format (1x, 71('-'))

c     Save filenames of feff.dat files
      open (unit=2, file='files.dat', status='unknown', iostat=ios)
      call chopen (ios, 'files.dat', 'genfmt')
c     Put phase header on top of files.dat
      do 200  itext = 1, ntext
         ltxt = istrln( text(itext))
         write(2,160)  text(itext)(1:ltxt)
  200 continue
      write(2,170)
      write(2,210)
  210 format ('    file        sig2   amp ratio    ',
     1        'deg    nlegs  r effective')
c     do each path
      call wlog ('    path     filename')

      do 500  ilist = 1, ntotal
c        find index of path
         do 410  j = 1, nptot
            if (iplst(ilist) .eq. index(j))  then
               i = j
               goto 430
            endif
  410    continue
         write(slog,420)  ilist, iplst(ilist)
  420    format (' did not find path i, iplst(i) ', 2i10)
         call wlog(slog)
  430    continue
c        Path i is the path from feff.bin that corresponds to
c        the path ilist in list.dat.  The index of the path is
c        iplst(ilist) and index(i).

c        Prepare output file feffnnnn.dat
         write(fname(i),220)  index(i)
  220    format ('feff', i4.4, '.dat')
         write(slog,230)  i, fname(i)
  230    format (i8, 5x, a)
         call wlog(slog)
c        zero is debye-waller factor column
         write(2,240) fname(i), zero, crit(i), deg(i),
     1                   nleg(i), reff(i)*bohr
  240    format(1x, a, f8.5, 2f10.3, i6, f9.4)

         ip = i
c     Write feff.dat's
         open (unit=3, file=fname(ip), status='unknown', iostat=ios)
         call chopen (ios, fname(ip), 'feffdt')
c        put header on feff.dat
         do 300  itext = 1, ntext
            ltxt = istrln(text(itext))
            write(3,160)  text(itext)(1:ltxt)
  300    continue
         write(3,310) ip, iorder
  310    format (' Path', i5, '      icalc ', i7)
         write(3,170)
         write(3,320)  nleg(ip), deg(ip), reff(ip)*bohr, rnrmav, 
     1                 edge*hart
  320    format (1x, i3, f8.3, f9.4, f10.4, f11.5, 
     1           ' nleg, deg, reff, rnrmav(bohr), edge')
         write(3,330)
  330    format ('        x         y         z   pot at#')
         write(3,340)  (rat(j,nleg(ip),ip)*bohr,j=1,3), 
     1                 ipot(nleg(ip),ip),
     1                 iz(ipot(nleg(ip),ip)), potlbl(ipot(nleg(ip),ip))
  340    format (1x, 3f10.4, i3, i4, 1x, a6, '   absorbing atom')
         do 360  ileg = 1, nleg(ip)-1
            write(3,350)  (rat(j,ileg,ip)*bohr,j=1,3), ipot(ileg,ip),
     1                    iz(ipot(ileg,ip)), potlbl(ipot(ileg,ip))
  350       format (1x, 3f10.4, i3, i4, 1x, a6)
  360    continue

         write(3,370)
  370    format    ('    k   real[2*phc]   mag[feff]  phase[feff]',
     1              ' red factor   lambda     real[p]@#')

c        Make the feff.dat stuff and write it to feff.dat
c        Also write out for inspection to fort.66
c        note that dimag takes complex*16 argument, aimag takes
c        single precision complex argument.  Stuff from feff.bin
c        is single precision, cchi is complex*16
         do 450  ie = 1, ne
c           Consider chi in the standard XAFS form.  Use R = rtot/2.
            cchi = achi(ie,ip) * exp (coni*phchi(ie,ip))
            xlam = 1.0e10
            if (abs(aimag(ck(ie))) .gt. eps) xlam = 1/aimag(ck(ie))
            redfac = exp (-2 * aimag (phc(ie)))
            cdelt = 2*dble(phc(ie))
            cfms = cchi * xk(ie) * reff(ip)**2 *
     1           exp(2*reff(ip)/xlam) / redfac
            if (abs(cchi) .lt. eps)  then
               phff = 0
            else
               phff = atan2 (dimag(cchi), dble(cchi))
            endif
c           remove 2 pi jumps in phases
            if (ie .gt. 1)  then
               call pijump (phff, phffo)
               call pijump (cdelt, cdelto)
            endif
            phffo = phff
            cdelto = cdelt

c           write 1 k, momentum wrt fermi level k=sqrt(p**2-kf**2)
c                 2 central atom phase shift (real part),
c                 3 magnitude of feff,
c                 4 phase of feff,
c                 5 absorbing atom reduction factor,
c                 6 mean free path = 1/(Im (p))
c                 7 real part of local momentum p

            write(3,400)
     1         xk(ie)/bohr,
     2         cdelt + l0*pi,
     3         abs(cfms) * bohr,
     4         phff - cdelt - l0*pi,
     5         redfac,
     6         xlam * bohr,
     7         dble(ck(ie))/bohr
  400       format (1x, f6.3, 1x, 3(1pe11.4,1x),1pe10.3,1x,
     1                            2(1pe11.4,1x))

  450    continue

c        Done with feff.dat
         close (unit=3)
  500 continue
      close (unit=2)

      return
      end
      subroutine ff2afs (ipr4, idwopt, critcw, s02, sig2g,
     1                   tk, thetad, mbconv, absolu,  !KJ added absolu 4-06
     1                   vrcorr, vicorr, alphat, thetae, iabs, nabs,
     4            elnes,ipmin,ipmax,ipstep)   !KJ added this line  4-06
c     calculate anomalous scattering amplitude for a given edge
c     Writes down main output: chi.dat and xmu.dat
      implicit double precision (a-h, o-z)

c={../HEADERS/const.h
      parameter (pi = 3.14159 26535 89793 23846 26433d0)
      parameter (one = 1, zero = 0)
      parameter (third = one/3)
      parameter (raddeg = 180 / pi)
      complex*16 coni
      parameter (coni = (0,1))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919 158 292 677 512 811d0)

      parameter (bohr = 0.529 177 249d0, ryd  = 13.605 698d0)
      parameter (hart = 2 * ryd)
      parameter (alpinv = 137.035 989 56d0)
c     fine structure alpha
      parameter (alphfs = 1 / alpinv)
c= ../HEADERS/const.h}
c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}
      parameter (eps4 = 1.0e-4)
      integer ipmin,ipmax,ipstep,elnes !KJ my variables  4-06    
      integer absolu !KJ 4-06  


c     header from list.dat
      dimension lhead(nheadx)
      character*80  head(nheadx)
      complex*16 gtr(nex),gtrful(ipmin:ipmax,nex) !KJ added gtrful 4-06

      parameter (npx=15000)
!KJ      parameter (npx = 1200)
c     indices of paths to do, read from list.dat
      dimension ip(npx)
      real sig2u(npx)

      complex*16 cchi(nex), ckp
c     to keep Im part of cchi 11.18.97 ala
      dimension rchtot(nex), xkp(nex)
      complex*16 chia(nex)

      logical dwcorr
      character*512 slog
      character*2 coment
      parameter (coment='# ')

c     Stuff from feff.bin, note that floating point numbers are
c     single precision.  Be careful throughout this routine, especially
c     when passing things to subroutines or intrinsic functions.
      real rnrmav, xmu, edge
      character*80 title(nheadx), titfms
      character*6  potlbl(0:nphx)
      dimension iz(0:nphx)
c     central atom phase shift at l0
      complex phc(nex)
      complex ck(nex)
      real xk(nex)
      dimension index(npx)
      dimension nleg(npx)
      real deg(npx), reff(npx), crit(npx)
      dimension ipot(legtot,npx)
      real rat(3,legtot,npx), beta(legtot,npx), eta(legtot,npx)
      real ri(legtot,npx)
      real achi(nex,npx), phchi(nex,npx)

c     stuff from xsect.bin
      complex*16 emxs(nex), xsec(nex)
      dimension omega(nex), xkxs(nex), xsnorm(nex), fpp(nex)
      
c !KJ locals  1-06
      integer iip,nip
      logical cross 
      character*9 f1,f2
      character*10 f0,f3
      complex*16 gtrtemp(nex*(1+ipmax-ipmin))
      complex*16 kxsec(nex)     
c !KJ end my variables




c     get gtr - result of FMS
      do 112 ie =1,nex
      gtr(ie)= 0
      do 112 iip =ipmin,ipmax  !KJ I added this variable ip 1-06
  112 gtrful(iip,ie) = 0
      ntfms = 0
      nip=ipmax-ipmin+1 !KJ 1-06
      
      open (unit=1, file='fms.bin', status='old', iostat=ios)
      if (ios.le.0) then
         ntfms = 1
         read(1, 113) titfms
  113    format(a)
         read(1, 115) ne, ne1, ne3, nph, npadx
  115    format(5(1x,i3))
         call rdpadx(1, npadx, gtrtemp, ne*nip)  !KJ I added *nip, changed gtr to gtrtemp  1-06
      endif
      close (unit=1)
!KJ Next lines my addition to read several spectra at once. 1-06
	i=1
	do iip=ipmin,ipmax
	do j=1,ne
        gtrful(iip,j)=gtrtemp(i+j-1)
	enddo
	  i=i+ne
	enddo
!KJ Now we don't need gtrtemp anymore.  End my changes



c !KJ loop over iip added to process several spectra at once  4-06
c !KJ reading of feff.bin and list.dat moved inside the loop (used to be before reading
c !KJ xsect.bin      
      do iip=ipmin,ipmax,ipstep
        cross=(.not.(iip.eq.1.or.iip.eq.10.or.iip.eq.5.or.iip.eq.9))
      
c !KJ choose different filename for each spectrum.
        if(iip.eq.1) then
	  f1(1:9)='chi.dat  '
	  f2(1:9)='xmu.dat  '
	  f0(1:10)='feff.bin  '
	  f3(1:10)='list.dat  ' 	  
	elseif(iip.eq.10) then
	  f1(1:9)='chi10.dat'
	  f2(1:9)='xmu10.dat'
	  f0(1:10)='feff10.bin'
	  f3(1:10)='list10.dat'	  
	elseif(iip.gt.1.and.iip.lt.10) then
	  f1(1:4)='chi0'
	  f1(5:5)= char(48+iip)
	  f1(6:9)='.dat'
	  f2(1:4)='xmu0'
	  f2(5:5)= char(48+iip)
	  f2(6:9)='.dat'
	  f0(1:5)='feff0'
	  f0(6:6)= char(48+iip)
	  f0(7:10)='.bin'	
	  f3(1:5)='list0'
	  f3(6:6)= char(48+iip)
	  f3(7:10)='.dat'
	else
	  stop 'crazy iip in ff2xmu'
	endif
	do i=1,nex
	gtr(i)=gtrful(iip,i)
	enddo
!KJ


c     open list.dat and read list of paths we want
      open (unit=1, file= f3, status='old', iostat=ios)!KJ changed 'list.dat' to f3 1-06
      ntotal = 0
      if (ios.le.0) then
        call chopen (ios, f3, 'ff2afs')  !KJ id.
        nhead = nheadx
        call rdhead (1, nhead, head, lhead)
c       skip a label line
        read(1,*)
c       ip is index of path, sig2u is debye-waller from user
        do 100  i = 1, npx
           read(1,*,end=110)  ip(i), sig2u(i)
           ntotal = i
  100   continue
  110   continue
      endif
      close (unit=1)


c     lines below allow to skip FMS module for DANES 
c     after XANES calculations
      open (unit=1, file='phase.bin', status='old', iostat=ios)
      if (ios.le.0) then
        read(1,*) ne3, ne3, ne3, ne3
      endif
      close (unit=1)

       call rdfbin (f0, nphx, nex, npx, legtot,   !KJ changed 'feff.bin' to f0  4-06
     $     nptot, ne, npot, ihole, iorder, ilinit, 
     $     rnrmav, xmu, edge, potlbl, iz, phc, ck, xk,
     $     index, nleg, deg, reff,
     $     crit, ipot, rat, beta, eta, ri, achi, phchi)

c     read xsect.bin file
      call  rdxbin (s02p, erelax, wp, edgep, s02, gamach, ne1, ik0,
     2  emxs, omega, xkxs, xsnorm, xsec, nxsec, mbconv, title, ntitle)
cc !!KJ these comments copied from ff2chi and not necessarily relevant  4-06.
c !KJ I have put rdxbin inside the loop since omega is 'recycled' below,
c !KJ which is a problem if the loop executes more than once.
c !KJ Simply reading the file again and again is the lazy solution,
c !KJ but it avoids confusing changes to the code (eg., new variables).



c     make combined title
      if (ntfms.eq.1) then
        ntitle = ntitle + 1
        title(ntitle) = titfms
      endif
      do 120 ihead = 1, nhead
  120 title(ntitle+ihead) = head(ihead)
      ntitle = ntitle + nhead

c     write feffnnnn.dat
      if (ipr4.eq.3) then
         call feffdt(ntotal,ip,nptot,ntitle,title,ne,npot,
     $        ihole, iorder, ilinit, rnrmav, xmu, edge, potlbl,
     $        iz,phc,ck,xk,index,
     $        nleg,deg,nepts,reff,crit,ipot,rat,achi,phchi)
       end if

c     If there is a vicorr, will need a mean free path factor xlam0.
c     Use it as  chi(ie) * exp (2 * reff * xlam0)
c     ckp is ck' = ck prime.
      if (abs(vicorr) .ge. eps4) then
         do 180  ie = 1, ne
            ckp = sqrt (ck(ie)**2 + 2*coni*vicorr)
            xlam0 = aimag(ck(ie)) - dimag(ckp)
            do 170  ipath = 1, nptot
               achi(ie,ipath) = achi(ie,ipath) * 
     1               exp (2 * reff(ipath) * xlam0)
  170       continue
  180    continue
      endif

c     k'**2 = k**2 + vr. If there is no real correction
c     (vrcorr = 0), these two grids will be the same.
c           k' is value for output,  k is  value used for
c           interpolations with original grid.

c     vrcorr shifts the edge and the k grid
      if (abs(vrcorr) .gt. eps4)  then
         edge = edge - vrcorr
      endif

c     ik0 is index at fermi level
      do 250  i = 1, ne
         temp = xk(i)*abs(xk(i)) + 2*vrcorr
         if (temp.ge. 0) then
           xkp(i) = sqrt(temp)
         else
           xkp(i) = - sqrt(-temp)
         endif
  250 continue
     

      dwcorr = .false.
      if (tk .gt. 1.0e-3)  dwcorr = .true.

c     Open chi.dat and xmu.dat (output) and start headers
      if (iabs.eq.nabs) then
         open (unit=3, file=f1, status='unknown', iostat=ios) !KJ changed chi.dat to f1
         call chopen (ios, f1, 'ff2afs') !KJ id.  4-06
         open (unit=8, file=f2, status='unknown', iostat=ios) !KJ changed xmu.dat to f2
         call chopen (ios, f2, 'ff2afs') !KJ id.

c        write miscellaneous staff into headers  !KJ corrected typo
         call wrhead (8, ntitle, title, dwcorr, s02,
     1     tk, thetad, sig2g, alphat, vrcorr, vicorr, critcw)

c        also write information on the screen
         if (alphat .gt. zero)  then
            write(slog,322) alphat
  322       format ('    1st and 3rd cumulants, alphat = ', 1pe20.4)
            call wlog(slog)
         endif
         if (abs(vrcorr).ge.eps4 .or. abs(vicorr).ge.eps4)  then
            write(slog,343) vrcorr*hart, vicorr*hart
  343       format ('    Energy zero shift, vr, vi ', 1p, 2e14.5)
            call wlog(slog)
         endif

         write(slog,370) critcw
         call wlog(slog)
  370    format ('    Use all paths with cw amplitude ratio', f7.2, '%')
         if (dwcorr)  then
            write(slog,380) s02, tk, thetad, sig2g
            call wlog(slog)
         else
            write(slog,381) s02, sig2g
            call wlog(slog)
         endif
  380    format('    S02', f7.3, '  Temp', f8.2, '  Debye temp', f8.2,
     1           '  Global sig2', f9.5)
  381    format('    S02', f7.3, '  Global sig2', f9.5)
      endif


c     make chi and sum it
      do 400  i = 1, nex
         cchi(i) = 0
  400 continue
      do 402  ik = 1, ne
         cchi(ik)= s02 * gtr(ik)
  402 continue
c     add correction due to vicorr
      if (vicorr.gt.eps4) then
         call conv(omega,cchi,ne1,vicorr)
c        call conv(omega,xsec,ne1,vicorr)
      endif


c     add Debye-Waller factors
      ispec = 3
      call dwadd (ntotal, nptot, idwopt, ip, index, crit, critcw, sig2g,
     1  sig2u, dwcorr, rnrmav, nleg, deg, reff, iz, ipot, rat,tk,thetad,
     2  alphat, thetae, mbconv, s02, ne1, ck, achi, phchi, ne, xk, xkp,
     3  xkp, cchi, iabs, nabs, ispec, ipr4, ntitle,
     4  title, vrcorr, vicorr,  nused)

c     read or initialize chia - result of configuration average
      if (iabs.eq.1) then
         do 635 ie =1, nex
            chia(ie) = 0
  635    continue
      else
         open (unit=1, file='chia.bin', status='old',
     1   access='sequential', form='unformatted', iostat=ios)
         do 640 ie = 1,ne
  640    read(1) chia(ie)
         close (unit=1, status='delete')
      endif

      if(iabs.eq.1) then
c        compare grids in xsect.bin and feff.bin
         do 680 i = 1, nxsec
           del = xk(i)**2 - xkxs(i)**2
           if (abs(del) .gt.  10*eps4)  then
             call wlog(' Emesh in feff.bin and xsect.bin different.')
c            stop 
           endif
  680    continue
      endif

c     add contribution from an absorber iabs 
c     present scheme assumes that xsec is the same for all iabs.
      do 701 ik = 1, ne
         chia(ik)   = chia(ik)   + cchi(ik)/ nabs
  701 continue
      if (iabs.lt.nabs) then
c        save chia in chia.bin for averaging
         open (unit=1, file='chia.bin', status='unknown',
     1   access='sequential', form='unformatted', iostat=ios)
         do 760 ie=1,ne
  760    write(1) chia(ie)
         close(unit=1)
      endif

      if (iabs.eq.nabs) then
c        The loop over absorbers is finished. Write out the results.
         write(8,600) coment, nused, ntotal
  600    format (a2,1x, i4, '/', i4, ' paths used')
  610    format (a2,1x, 71('-'))

         do 702 ik = 1, ne
            rchtot(ik) = dimag (chia(ik))
  702    continue
c        prepare the output grid omega
         efermi = edge + omega(1) - dble(emxs(1))

c        do convolution with excitation spectrum
         if (mbconv .gt. 0) then
            wp = wp / 2.
            call  exconv
     1      (omega, ne1, efermi, s02p, erelax, wp, xsnorm)
            call  exconv
     1      (omega, ne1, efermi, s02p, erelax, wp, rchtot)
         endif

c        normalize to xsec at 50 eV above edge
c        and prepare the output energy grid omega
         edg50 = efermi + 50 / hart
         call terp (omega, xsnorm,  ne1, 1, edg50, xsedge)
         if (absolu.eq.1) xsedge=dble(1)  !KJ 4-06 don't normalize	 
         write(8,660)  coment, xsedge 
  660    format (a2,' xsedge+ 50, used to normalize mu ', 1pe20.4)
         write(8,610) coment
         write(8,665) coment
  665    format (a2,' omega    e    k    mu    mu0     chi     @#')

c        transform from cross section in Angstrom**2 to f"/m*c**2
         do 670 ie = 1,ne
            energy = dble(emxs(ie)) + efermi 
            prefac = 4 * pi * alpinv / energy * bohr**2
c           add alpha**2 to convert to units for f'
            xsec(ie) = xsec(ie) / prefac * alpinv**2
            xsnorm(ie) = xsnorm(ie) / prefac * alpinv**2
  670    continue

         do i=1,nex
	 if (.not.cross) then !KJ I added this block 4-06
	   kxsec(i)=xsec(i)
	 else
	   kxsec(i)=dcmplx(0,0)
	 endif !KJ end my code
	 enddo


c        do correction using brouder method
         ne2 = ne - ne1 - ne3
         call fprime(efermi, emxs, ne1, ne3, ne, ik0, kxsec,xsnorm,chia, !KJ changed xsec to kxsec 4-06
     1       vrcorr, vicorr, cchi)
         do 850 ie=1,ne1
           fpp(ie)=xsnorm(ie) + dimag(xsnorm(ie)*chia(ie))
           rchtot(ie)=dble(xsnorm(ie)*chia(ie)+cchi(ie))
  850    continue
         do 855 ie=1,ne
           chia(ie) = 0
  855    continue
         call fprime(efermi, emxs, ne1, ne3, ne, ik0, kxsec,xsnorm,chia, !KJ id.
     1       vrcorr, vicorr, cchi)
         do 860 ie = 1, ne1
            em0 = dble(emxs(ie))
            xsec0 = dble( cchi(ie))
            chi0  = (rchtot(ie) - xsec0)
            if (ne2.gt.0) then
c             DANES
c             - signs to comply with Cromer-Liberman notation for f', f"
              write(8,700)  omega(ie)*hart, em0*hart, xkp(ie)/bohr,
     1             -rchtot(ie), -xsec0, -chi0
  700         format (1x, 2f11.3, f8.3, 1p, 3e13.5)
            else
c             FPRIME
              write(8,710)  omega(ie)*hart, em0*hart, 
     1              -rchtot(ie), -xsec0, fpp(ie), xsnorm(ie)
  710         format (1x, 2f11.3, 4e13.5)
            endif
  860    continue

         close (unit=8)
         close (unit=3, status='delete')
      endif
c     for if (iabs=abs); or the last absorber

      enddo !KJ of my iip=ipmin,ipmax,ipstep loop  1-06


      return
      end
      subroutine ff2chi (ispec, ipr4, idwopt, critcw, s02, sig2g,
     1                   tk, thetad, mbconv, absolu,  !KJ added absolu 3-06
     1                   vrcorr, vicorr, alphat, thetae, iabs, nabs,
     4            elnes,ipmin,ipmax,ipstep)   !KJ added this line  1-06     
c     adds the contributions from each path and absorber, including
c     Debye-Waller factors. Writes down main output: chi.dat and xmu.dat
      implicit double precision (a-h, o-z)

c={../HEADERS/const.h
      parameter (pi = 3.14159 26535 89793 23846 26433d0)
      parameter (one = 1, zero = 0)
      parameter (third = one/3)
      parameter (raddeg = 180 / pi)
      complex*16 coni
      parameter (coni = (0,1))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919 158 292 677 512 811d0)

      parameter (bohr = 0.529 177 249d0, ryd  = 13.605 698d0)
      parameter (hart = 2 * ryd)
      parameter (alpinv = 137.035 989 56d0)
c     fine structure alpha
      parameter (alphfs = 1 / alpinv)
c= ../HEADERS/const.h}
c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}
      parameter (eps4 = 1.0e-4)
      integer ipmin,ipmax,ipstep,elnes !KJ my variables  1-06    
      integer absolu !KJ 3-06  

c     header from list.dat
      dimension lhead(nheadx)
      character*80  head(nheadx)

      parameter (npx=15000)
!KJ      parameter (npx = 1200)
c     indices of paths to do, read from list.dat
      dimension ip(npx)
      real sig2u(npx)

      parameter (nfinex = 601)
      complex*16 cchi(nfinex), ckck(nfinex), ccc, ckp
c     to keep Im part of cchi 11.18.97 ala
      dimension rchtot(nfinex)
      complex*16 chia(nfinex)
      dimension xkp(nfinex), xk0(nfinex)

      logical dwcorr
      character*512 slog
      character*2 coment
      parameter (coment='# ')

c     Stuff from feff.bin, note that floating point numbers are
c     single precision.  Be careful throughout this routine, especially
c     when passing things to subroutines or intrinsic functions.
      real rnrmav, xmu, edge
      character*80 title(nheadx)
      character*6  potlbl(0:nphx)
      dimension iz(0:nphx)
c     central atom phase shift at l0
      complex phc(nex)
      complex ck(nex)
      real xk(nex)
      dimension index(npx)
      dimension nleg(npx)
      real deg(npx), reff(npx), crit(npx)
      dimension ipot(legtot,npx)
      real rat(3,legtot,npx), beta(legtot,npx), eta(legtot,npx)
      real ri(legtot,npx)
      real achi(nex,npx), phchi(nex,npx)

c     stuff from xsect.bin
      complex*16 emxs(nex), xsec(nex)
      dimension omega(nex), xkxs(nex), xsnorm(nex), fpp(nex)
      dimension omegax(nfinex)
c#mn
      external getxk

c !KJ locals  1-06
      integer iip,nip
      logical cross 
      character*9 f1,f2
      character*10 f0,f3
      complex*16 kxsec(nex)     
!KJ end my variables      
      

c     lines below allow to skip FMS module for DANES
c     after XANES calculations
      open (unit=1, file='phase.bin', status='old', iostat=ios)
      if (ios.le.0 .and. abs(ispec).eq.3) then
        read(1,*) ne3, ne3, ne3, ne3
      endif
      close (unit=1)


c !KJ loop over iip added to process several spectra at once  1-06
c !KJ reading of feff.bin and list.dat moved inside the loop (used to be before reading
c !KJ xsect.bin      
      do iip=ipmin,ipmax,ipstep
        cross=(.not.(iip.eq.1.or.iip.eq.10.or.iip.eq.5.or.iip.eq.9))
      
c !KJ choose different filename for each spectrum.
        if(iip.eq.1) then
	  f1(1:9)='chi.dat  '
	  f2(1:9)='xmu.dat  '
	  f0(1:10)='feff.bin  '
	  f3(1:10)='list.dat  ' 	  
	elseif(iip.eq.10) then
	  f1(1:9)='chi10.dat'
	  f2(1:9)='xmu10.dat'
	  f0(1:10)='feff10.bin'
	  f3(1:10)='list10.dat'	  
	elseif(iip.gt.1.and.iip.lt.10) then
	  f1(1:4)='chi0'
	  f1(5:5)= char(48+iip)
	  f1(6:9)='.dat'
	  f2(1:4)='xmu0'
	  f2(5:5)= char(48+iip)
	  f2(6:9)='.dat'
	  f0(1:5)='feff0'
	  f0(6:6)= char(48+iip)
	  f0(7:10)='.bin'	
	  f3(1:5)='list0'
	  f3(6:6)= char(48+iip)
	  f3(7:10)='.dat'
	else
	  stop 'crazy iip in ff2xmu'
	endif

c     open list.dat and read list of paths we want
      open (unit=1, file = f3, status='old', iostat=ios) !KJ changed 'list.dat' to f3 1-06
      call chopen (ios, f3, 'ff2chi') !KJ id.
      nhead = nheadx
      call rdhead (1, nhead, head, lhead)
c     skip a label line
      read(1,*)
      ntotal = 0
c     ip is index of path, sig2u is debye-waller from user
      do 100  i = 1, npx
         read(1,*,end=110)  ip(i), sig2u(i)
         ntotal = i
  100 continue
  110 continue
      close (unit=1)
      
      
c     read 'xsect.bin'
      call  rdxbin (s02p, erelax, wp, edgep, s02, gamach, ne1, ik0,
     2  emxs, omega, xkxs, xsnorm, xsec, nxsec, mbconv, title, ntitle)
c !KJ I have put rdxbin inside the loop since omega is 'recycled' below,
c !KJ which is a problem if the loop executes more than once.
c !KJ Simply reading the file again and again is the lazy solution,
c !KJ but it avoids confusing changes to the code (eg., new variables).

       call rdfbin (f0, nphx, nex, npx, legtot,  !KJ changed 'feff.bin' to f0  1-06
     $      nptot, ne, npot, ihole, iorder, ilinit, 
     $      rnrmav, xmu, edge, potlbl, iz, phc, ck, xk, index, 
     $      nleg, deg, reff, crit, ipot, 
     $      rat, beta, eta, ri, achi, phchi)
     
c     make combined title
      do 120 ihead = 1, nhead
  120 title(ntitle+ihead) = head(ihead)
      ntitle = ntitle + nhead

c     write feffnnnn.dat
      if (ipr4.eq.3) then
         call feffdt(ntotal,ip,nptot,ntitle,title,ne1,npot,
     $        ihole, iorder, ilinit, rnrmav, xmu, edge, potlbl,
     $        iz,phc,ck,xk,index,
     $        nleg,deg,nepts,reff,crit,ipot,rat,achi,phchi)
       end if

      if (iabs.eq.1) then
c        compare grids in xsect.bin and feff.bin
         do 680 i = 1, nxsec
           del = xk(i)**2 - xkxs(i)**2
           if (abs(ispec).ne.3 .and. abs(del) .gt. 10*eps4)  then
             call wlog(' Emesh in feff.bin and xsect.bin different.')
             call wlog
     1       (' Results may be meaningless, check input files.')
             call wlog
     1       (' Either use XANES card or remove xsect.bin file.')
             write(slog,670)  i, xk(i)/bohr, xkxs(i)/bohr, del
             call wlog(slog)
  670        format(i7, 1p, 3e13.5)
             call par_stop('FF2CHI-1') 
           endif
  680    continue
      endif

c     If there is a vicorr, will need a mean free path factor xlam0.
c     Use it as  chi(ie) * exp (2 * reff * xlam0)
c     ckp is ck' = ck prime.
      if (abs(vicorr) .ge. eps4) then
         do 170  ipath = 1, nptot
            do 180  ie = 1, ne
               ckp = sqrt (ck(ie)**2 + coni*2*vicorr)
               xlam0 = aimag(ck(ie)) - dimag(ckp)
               achi(ie,ipath) = achi(ie,ipath) * 
     1              exp (2 * reff(ipath) * xlam0)
 180        continue
 170     continue
      endif

c     Decide on fine grid.  We need two, k' evenly spaced by 
c     delk (0.05 invA) and k0 being the place in the original k 
c     grid corresponding to each k'.  k0 will not in general be on 
c     an original grid point.  Define k' by k'**2 = k**2 + vr.
c     If there is no real correction (vrcorr = 0), these two grids
c     will be the same.
c           k' is value for output, k0 is k value used for
c           interpolations with original grid.

c     vrcorr shifts the edge and the k grid
      if (abs(vrcorr) .gt. eps4)  then
         edge = edge - vrcorr
      endif

c     Find xkmin, beginning of k' grid
      delk = 0.05 * bohr
      tmp = sign (real(one), xk(1))
      e = tmp * xk(1)**2 / 2 + vrcorr
      xkpmin = getxk (e)
      n = xkpmin / delk
c     need 1st int ABOVE xkpmin/delk
      if (xkpmin .gt. 0)  n = n + 1
c     First k grid point moved by vrcorr
      xkmin = n * delk

c     Make xkp (k') and xk0 (k0) fine grids
c     ik0 is index at fermi level
      if (abs(ispec).ne.3) ik0 = 1
      ik0p = 1
      do 250  i = 1, nfinex
         xkp(i) = xkmin + delk * (i - 1)
         tmp = sign (one, xkp(i))
         e = tmp * xkp(i)**2 /2 - vrcorr
         xk0(i) = getxk(e)
         if (xk0(i).lt.eps4)  ik0p = i
         if (xk0(i) .gt. xk(ne1)+eps4)  goto 260
         nkx = i
  250 continue
  260 continue

      dwcorr = .false.
      if (tk .gt. 1.0e-3)  dwcorr = .true.

c     Open chi.dat and xmu.dat (output) and start headers
      if (iabs.eq.nabs) then
         open (unit=3, file=f1, status='unknown', iostat=ios) !KJ changed chi.dat to f1 1-06
         call chopen (ios, f1, 'ff2chi') !KJ id.
         open (unit=8, file=f2, status='unknown', iostat=ios) !KJ changed xmu.dat to f2  1-06
         call chopen (ios, f2, 'ff2chi') !KJ id.

c        write miscellaneous staff into headers  !KJ corrected typo
         call wrhead (3, ntitle, title, dwcorr, s02,
     1        tk, thetad, sig2g, alphat, vrcorr, vicorr, critcw)

         call wrhead (8, ntitle, title, dwcorr, s02,
     1        tk, thetad, sig2g, alphat, vrcorr, vicorr, critcw)

c        also write information on the screen
         if (alphat .gt. zero)  then
            write(slog,322) alphat
  322       format ('    1st and 3rd cumulants, alphat = ', 1pe20.4)
            call wlog(slog)
         endif
         if (abs(vrcorr).ge.eps4 .or. abs(vicorr).ge.eps4)  then
            write(slog,343) vrcorr*hart, vicorr*hart
  343       format ('    Energy zero shift, vr, vi ', 1p, 2e14.5)
            call wlog(slog)
         endif

         write(slog,370) critcw
         call wlog(slog)
  370    format ('    Use all paths with cw amplitude ratio', f7.2, '%')
         if (dwcorr)  then
            write(slog,380) s02, tk, thetad, sig2g
  380       format('    S02', f7.3, '  Temp', f8.2, '  Debye temp',f8.2,
     1           '  Global sig2', f9.5)
            call wlog(slog)
         else
            write(slog,381) s02, sig2g
  381       format('    S02', f7.3, '  Global sig2', f9.5)
            call wlog(slog)
         endif
      endif


c     make chi and sum it
      do 400  i = 1, nfinex
         cchi(i) = 0
  400 continue

c     add Debye-Waller factors
      call dwadd (ntotal, nptot, idwopt, ip, index, crit, critcw, sig2g,
     1  sig2u, dwcorr, rnrmav, nleg, deg, reff, iz, ipot,rat, tk,thetad,
     2  alphat, thetae, mbconv, s02, ne1, ck, achi, phchi, nkx, xk, xk0,
     3  xkp, cchi, iabs, nabs, ispec, ipr4, ntitle,
     4  title, vrcorr, vicorr, nused)

c     read or initialize chia - result of configuration average
      if (iabs.eq.1) then
         do 635 ie =1, nfinex
            chia(ie) = 0
  635    continue
      else
         open (unit=1, file='chia.bin', status='old',
     1   access='sequential', form='unformatted', iostat=ios)
         do 640 ie = 1,nkx
  640    read(1) chia(ie)
         close (unit=1, status='delete')
      endif

c     add contribution from an absorber iabs 
c     present scheme assumes that xsec is the same for all iabs.
      do 701 ik = 1, nkx
         chia(ik)   = chia(ik)   + cchi(ik)/ nabs
  701 continue
      if (iabs.lt.nabs) then
c        save chia in chia.bin for averaging
         open (unit=1, file='chia.bin', status='unknown',
     1   access='sequential', form='unformatted', iostat=ios)
         do 760 ie=1,nkx
  760    write(1) chia(ie)
         close(unit=1)
      endif

      if (iabs.eq.nabs) then
c        the loop over absorbers finished, ready to report results

c        Write it out
         write(3,600)  coment, nused, ntotal
         write(8,600)  coment, nused, ntotal
  600    format (a2, 1x, i4, '/', i4, ' paths used')
         write(3,610) coment
  610    format (a2, 1x, 71('-'))
         write(3,620) coment
  620    format(a2,
     1         '      k          chi          mag           phase @#')

         do 702 ik = 1, nkx
           if (abs(ispec).ne.3) then
            rchtot(ik) = dimag (chia(ik))
           else
            rchtot(ik) = dble (chia(ik))
           endif
  702    continue
c        prepare the output grid omegax
         efermi = edge + omega(1) - dble(emxs(1))
         do 590  ik = 1, nkx
            if (xkp(ik) .lt. 0.0) then
               omegax(ik) = - xkp(ik) * xkp(ik) / 2  + efermi
            else
               omegax(ik) = xkp(ik) * xkp(ik) / 2  + efermi
            endif
  590    continue

c        do convolution with excitation spectrum
c        it is currently screwed up since xsnorm is rewritten
c        fix later
         if (mbconv .gt. 0) then
            wp = wp / 2.
            call  exconv
     1      (omega, ne1, efermi, s02p, erelax, wp, xsnorm)
            call  exconv
     1      (omegax, nkx, efermi, s02p, erelax, wp, rchtot)
         endif


c        write to 'chi.dat'
         do 660 ik = 1, nkx
            ccc = chia(ik)
            phase = 0
            if (abs(ccc) .gt. 0)  then
               phase = atan2 (dimag(ccc), dble(ccc))
            endif
            if (ik .gt. 1)  call pijump (phase, phase0)
            phase0 = phase
            if (ipr4.ne.4) then
              write(3,630)  xkp(ik)/bohr, rchtot(ik), abs(ccc), phase0
  630         format (1x, f10.4, 3x, 3(1pe13.6,1x))
            else
c             need to report ck into chi.dat for Conradson's program
c             complex*16 should be used in terpc
              do 625 i=1,ne
  625         ckck(i) = dble(real(ck(i))) +coni*dble(aimag(ck(i)))
              call terpc (xkxs, ckck, ne, 3, xk0(ik), ckp)
              write(3,650)  xkp(ik)/bohr, rchtot(ik), abs(ccc), phase0,
     1        dble(ckp)/bohr, dimag(ckp)/bohr
  650         format (1x, f10.4, 3x, 5(1pe13.6,1x))
            endif
  660    continue
         close (unit=3)
   
c        write to 'xmu.dat'
c        normalize to xsec at 50 eV above edge
c        and prepare the output energy grid omegax
         edg50 = efermi + 50 / hart
         call terp (omega, xsnorm,  ne1, 1, edg50, xsedge)
	 if (absolu.eq.1) xsedge=dble(1) !KJ 1-06 don't normalize
         write(8,690)  coment, xsedge 
  690    format (a2, ' xsedge+50, used to normalize mu ', 1pe20.4)
         write(8,610) coment
         write(8,695) coment
  695    format (a2,' omega    e    k    mu    mu0     chi     @#')

         do i=1,nex
	 if (.not.cross) then !KJ I added this block 1-06
	   kxsec(i)=xsec(i)
	 else
	   kxsec(i)=dcmplx(0,0)
	 endif !KJ end my code
	 enddo
	 

c        do edge correction and write down results to xmu.dat, chi.dat
         do 710 ie = 1, ne
  710    chia(ie) = 0
         if (abs(ispec).eq.3) then
c          transform from cross section in Angstrom**2 to f"/m*c**2
           do 697 ie = 1,ne
             energy = dble(emxs(ie)) + efermi
             prefac = 4 * pi * alpinv / energy * bohr**2
c            add alpha**2 to convert to units for f'
             kxsec(ie) = kxsec(ie) / prefac * alpinv**2   !KJ changed xsec to kxsec  1-06
             xsnorm(ie) = xsnorm(ie) / prefac * alpinv**2
  697      continue
           ne2 = ne - ne1 - ne3
           call fprime(efermi, emxs, ne1, ne3,ne,ik0, kxsec,xsnorm,chia,
     1       vrcorr, vicorr, cchi)  !KJ changed xsec to kxsec 1-06
           do 850 ie=1,ne1
             omega(ie) = dble(cchi(ie))
  850      continue
         else
           call xscorr (ispec, emxs, ne1, ne, ik0, kxsec, xsnorm, chia,
     1       vrcorr, vicorr, cchi) !KJ xsec to kxsec 7/06
c          omega is not used as energy array, but as xsec array below
           do 711 ie = 1, ne1
  711      omega(ie) = dimag(kxsec(ie)+cchi(ie))  !KJ xsec to kxsec 7/06
         endif

         do 750  ik = 1, nkx
            em0 = omegax(ik) - efermi + edge
            call terp (xkxs, omega,  ne1, 1, xk0(ik), xsec0)
            call terp (xkxs, xsnorm,  ne1, 1, xk0(ik), xsnor0)
            if (omegax(ik).ge.efermi) then
              chi0 = xsnor0 * rchtot(ik)
            else
              chi0 = xsnor0 * rchtot(ik0p)
            endif
            if (abs(ispec).ne.3) then
              write(8,700)  omegax(ik)*hart, em0*hart, xkp(ik)/bohr,
     1              ( chi0 + dble(xsec0) )/xsedge,
     1              xsec0 /xsedge, rchtot(ik)
            else
c             signs to comply with Cromer-Liberman notation for f', f"
              write(8,700)  omegax(ik)*hart, em0*hart, xkp(ik)/bohr,
     1             -(xsec0+chi0), -xsec0, -chi0
            endif
  700       format (1x, 2f11.3, f8.3, 1p, 3e13.5)
  750    continue
         close (unit=8)
      endif
c     for if (iabs=abs); or the last absorber


      enddo !KJ of my iip=ipmin,ipmax,ipstep loop   1-06

      return
      end
      subroutine rdxbin (s02p, erelax, wp, edgep, s02, gamach, ne1, ik0,
     1                   emxs, omega, xkxs, xsnorm, xsec, nxsec, mbconv,
     2                   title, ntitle)

      implicit double precision (a-h, o-z)
c={../HEADERS/const.h
      parameter (pi = 3.14159 26535 89793 23846 26433d0)
      parameter (one = 1, zero = 0)
      parameter (third = one/3)
      parameter (raddeg = 180 / pi)
      complex*16 coni
      parameter (coni = (0,1))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919 158 292 677 512 811d0)

      parameter (bohr = 0.529 177 249d0, ryd  = 13.605 698d0)
      parameter (hart = 2 * ryd)
      parameter (alpinv = 137.035 989 56d0)
c     fine structure alpha
      parameter (alphfs = 1 / alpinv)
c= ../HEADERS/const.h}
c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}

c     header from xsect.bin
      dimension ltitle(nheadx)
      character*80 title(nheadx)
      complex*16 emxs(nex), xsec(nex)
      dimension omega(nex), xkxs(nex), xsnorm(nex)

      open (unit=8, file='xsect.bin', status='old', iostat=ios)
c     read xsect.bin
      ntitle = nheadx
      call rdhead (8, ntitle, title, ltitle)
c     read method for xsec calculation
      read(8,*)  s02p, erelax, wp, edgep, emu
      if (mbconv .gt.0 .or. s02.le.0.1) s02=s02p
c     read gamach (in eV) for use in atan at absorption edge
c     and convert to code units
      read(8,*)  gamach, ne1, ik0
      gamach = gamach / hart
c     skip label and read after it
      read(8,*)
      i = 1
  300    read(8,*,end=310)  ereal, eimag, xsnorm(i), dum1, dum2
         xsec(i) = dum1 + coni*dum2
c        xsect.bin is in eV and invA, convert to code units here
         emxs(i) = (ereal + coni*eimag) / hart
         xkxs(i) = getxk (dble(emxs(i)) - edgep)
         omega(i) = dble(emxs(i)) - edgep + emu
         nxsec = i
         i = i + 1
         if (i.le.nex) goto 300
  310 continue
      close(unit=8)
 
      return
      end

      subroutine wrhead (iunit, nhead, head, dwcorr, s02,
     1  tk, thetad, sig2g, alphat, vrcorr, vicorr, critcw)

      implicit double precision (a-h, o-z)

c={../HEADERS/const.h
      parameter (pi = 3.14159 26535 89793 23846 26433d0)
      parameter (one = 1, zero = 0)
      parameter (third = one/3)
      parameter (raddeg = 180 / pi)
      complex*16 coni
      parameter (coni = (0,1))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919 158 292 677 512 811d0)

      parameter (bohr = 0.529 177 249d0, ryd  = 13.605 698d0)
      parameter (hart = 2 * ryd)
      parameter (alpinv = 137.035 989 56d0)
c     fine structure alpha
      parameter (alphfs = 1 / alpinv)
c= ../HEADERS/const.h}
c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}
c={../HEADERS/vers.h
      character*12 vfeff
c                       123456789012  
      parameter (vfeff='Feff 8.50   ')
c= ../HEADERS/vers.h}
      parameter (eps4 = 1.0e-4)
      character*80  head(nheadx)
      logical dwcorr
      character*2 coment
      parameter (coment='# ')

c     write miscellanious staff into headers
c     add feff verdion to the first line
      ll = istrln(head(1))
      if (ll .lt. 0)  then
        head(1)= 'Untitled'
        ll = istrln(head(1))
      endif
      write(iunit,310) coment, head(1)(1:), vfeff
  310 format (a2, a55, t66, a12)

c     the rest of the title
      do 330  ihead = 2, nhead
         ll = istrln(head(ihead))
         if (ll .gt. 0)  then
            write(iunit,320) coment, head(ihead)(1:ll)
         endif
  320    format (a2, a)
  330 continue
      if (dwcorr)  then
         write(iunit,340)  coment, s02, tk, thetad, sig2g
  340    format (a2,' S02=', f5.3, '  Temp=', f7.2,'  Debye_temp=',f7.2,
     1        '  Global_sig2=', f8.5)
      else
         write(iunit,341)  coment, s02, sig2g
  341    format (a2, ' S02=', f5.3,
     1   '                                        Global_sig2=', f8.5)
      endif
      if (alphat .gt. zero)  then
         write(iunit,321)  coment, alphat
  321    format (a2, ' 1st and 3rd cumulants, alphat = ', 1pe20.4)
      endif

      if (abs(vrcorr).ge.eps4 .or. abs(vicorr).ge.eps4)  then
         write(iunit,342) coment, vrcorr*hart, vicorr*hart
      endif
  342 format (a2, ' Energy zero shift, vr, vi ', 1p, 2e14.5)

      if (critcw .gt. 0)  write(iunit,350) coment, critcw
  350 format (a2, ' Curved wave amplitude ratio filter ', f7.3, '%')
      write(iunit,360) coment
  360 format (a2, '    file         sig2 tot  cw amp ratio   deg',
     1        '  nlegs   reff  inp sig2')
c     stop writing misc. staff to files

      return
      end


      subroutine dwadd (ntotal,nptot,idwopt,ip,index,crit,critcw,sig2g,
     1  sig2u, dwcorr, rnrmav, nleg, deg, reff, iz, ipot, rat,tk,thetad,
     2  alphat, thetae, mbconv, s02, ne1, ck, achi, phchi, nkx, xk, xk0,
     3  xkp, cchi, iabs, nabs, ispec, ipr4, nhead,
     4  head, vrcorr, vicorr, nused)
      implicit double precision (a-h, o-z)

c={../HEADERS/const.h
      parameter (pi = 3.14159 26535 89793 23846 26433d0)
      parameter (one = 1, zero = 0)
      parameter (third = one/3)
      parameter (raddeg = 180 / pi)
      complex*16 coni
      parameter (coni = (0,1))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919 158 292 677 512 811d0)

      parameter (bohr = 0.529 177 249d0, ryd  = 13.605 698d0)
      parameter (hart = 2 * ryd)
      parameter (alpinv = 137.035 989 56d0)
c     fine structure alpha
      parameter (alphfs = 1 / alpinv)
c= ../HEADERS/const.h}
c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}
      parameter (eps4 = 1.0e-4)
      character*80  head(nheadx)
      parameter (npx=15000)
c     indices of paths to do, read from list.dat
      dimension ip(npx)
      real sig2u(npx)

      parameter (nfinex = 601)
      complex*16 cchi(nfinex), ccpath(nfinex), ccc, ckp
c     to keep Im part of cchi 11.18.97 ala
      complex*16 dw, dw1, dw3
      dimension xkp(nfinex), xk0(nfinex)

      logical dwcorr
      dimension rattmp(3,0:legtot)
      dimension iztmp(0:legtot)
      character*512 slog
      character*12 fname
      real rnrmav
      dimension iz(0:nphx)
c     central atom phase shift at l0
      complex ck(nex)
      real xk(nex)
      dimension index(npx)
      dimension nleg(npx)
      real deg(npx), reff(npx), crit(npx)
      dimension ipot(legtot,npx)
      real rat(3,legtot,npx)
      real achi(nex,npx), phchi(nex,npx)
      dimension sig2x(0:nphx, 0:nphx)
      character*2 coment
      parameter (coment='# ')

c     Keep stats on paths used to make chi
      nused = 0
      xkref = dble(ck(1)**2) - xk(1)*abs(xk(1)) 

c     open the files for sigrm and sigem
      if (idwopt.eq.1) then
         iem = 111
         open(unit=iem,file='s2_em.dat',status='unknown', iostat=ios)
         call chopen (ios, 's2_em.dat', 'sigem')
      elseif (idwopt.eq.2) then
         irm1 =111
         open(unit=irm1,file='s2_rm2.dat',status='unknown', iostat=ios)
         call chopen (ios, 's2_rm2.dat', 'sigrm')
         irm2 = 112
         open(unit=irm2,file='s2_rm1.dat',status='unknown', iostat=ios)
         call chopen (ios, 's2_rm1.dat', 'sigrm')
      endif
      if (alphat .gt. 0) then
        icum = 113
        open(unit=icum, file='cum.dat', status='unknown', iostat=ios)
        call chopen (ios, 'cum.dat', 'sig3')
        Write(icum, 363)
  363  format('# first and third icumulant for single scattering paths')
        write(icum,364) thetae, alphat
        write(icum,365)
  364   format ('# Einstein-Temp. =',f9.2 ,'   ', 'alpha=',f9.5)
  365   format ('#       file   sig1    sig2    sig3 ')
      endif

      if (idwopt.ge.1) then
c        initialize statistics for max DW for sigrm
         sig2mx=0
         do 400 iph1=0,nphx
         do 400 iph2=0,nphx
  400    sig2x(iph1, iph2) = 0
      endif
     

c     cycle over all paths in the list
      do 560  ilist = 1, ntotal
c        find index of path
         do 410  j = 1, nptot
            if (ip(ilist) .eq. index(j))  then
               ipath = j
               goto 430
            endif
  410    continue
         write(slog,420)  ilist, ip(ilist)
  420    format (' did not find path i, ip(i) ', 2i10)
         call wlog(slog)
  430    continue
c        Path ipath is the path from feff.bin that corresponds to
c        the path ilist in list.dat.  The index of the path is
c        ip(ilist) and index(ipath).

c        Use this path if it passes critcw filter
         if (crit(ipath) .lt. critcw)  goto 550

c        do debye-waller factors, get sig2d from correlated debye 
c        model if required
c        A note about units:  sig2g, sig2u() and sig2d are all in
c        Angs**2.  Convert to code units after we've summed them.
         sig2 = sig2g + sig2u(ilist)
         if (dwcorr .and. idwopt.ge.0)  then
c           note that stuff from feff.bin is single precision and
c           mostly in multidim. arrays.  sigms is double precision
c           and its arrays are dimensioned for a single path, so
c           use tmp variables to call it.  tk, thetad and sig2d are 
c           all dp, and therefore OK.  Also note that sigms takes
c           inputs in angstroms, except for rs which is in bohr.
            rs = rnrmav
            do 460  ileg = 1, nleg(ipath)
               iztmp(ileg) = iz(ipot(ileg,ipath))
               do 450  j = 1, 3
                  rattmp(j,ileg) = rat(j,ileg,ipath) * bohr
  450          continue
  460       continue
            iztmp(0) = iztmp(nleg(ipath))
            do 470  j = 1,3
               rattmp(j,0) = rattmp(j,nleg(ipath))
  470       continue
            if (idwopt.eq.0) then 
c             use CD model
              call sigms (tk, thetad, rs, legtot, nleg(ipath), 
     1                  rattmp, iztmp, sig2d)
            elseif (idwopt.eq.1) then 
c             use EM method
              call sigem
     1        (sig2mx,sig2x,iem,tk,ipath,nleg(ipath),rattmp,sig2d)
            elseif (idwopt.eq.3) then  !KJ 7/06 added this section
c             use CL model
              call sigcl (tk, thetad, rs, legtot, nleg(ipath), 
     1                  rattmp, iztmp, sig2d)   
            else 
c             use RM
              call sigrm
     1        (sig2mx,sig2x,irm1,irm2,tk,ipath,nleg(ipath),rattmp,sig2d)
            endif
            sig2 = sig2 + sig2d
         endif
         sig2 = sig2 / (bohr**2)

c        Do first and third cumulants
         sig1 = 0
         sig3 = 0
         if (alphat .gt. zero  .and. nleg(ipath) .eq. 2)  then
           if (thetae.le.0.d0) then
c            call sig3  to get sig1 and sig3 for single scattering paths
c           use reff(ipath) for r, note that reff is single precision
             iz1 = iztmp(nleg(ipath))
             iz2 = iztmp(1)
             call sigte3(iz1, iz2, sig2, alphat, thetad, reff(ipath),
     1                   sig1, sig3)
           else
c            this gets sig1 and sig3 for single scattering paths
c            using Morse potential
             call sigm3(sig1, sig2, sig3, tk, alphat, thetae)
           endif
           write(icum,475) index(ipath),  sig1 * bohr,
     1                 sig2*(bohr**2), sig3*(bohr**3)
  475      format( i10,f9.5,f9.5,' ',f9.7)
         endif

c        put the debye-waller factor and other cumulants into 
c        achi and phchi
         if (mbconv .gt. 0) s02 = 1.0
         do 480  i = 1, ne1
            dw = exp(-2 * sig2 * ck(i)**2)
            dw1 = exp (2 * coni * ck(i) * sig1)
            dw3 = exp ((-4 * coni * ck(i)**3 * sig3) / 3)
            dw = dw * dw1 * dw3
            phdw = 0.0
            if (abs(dw).gt.0) phdw = atan2 (dimag(dw), dble(dw))
            achi(i,ipath) = achi(i,ipath) * abs(dw) * s02 * deg(ipath)
            phchi(i,ipath) = phchi(i,ipath) + phdw
  480    continue
c        make sure no 2pi jumps in phase
         do 490  i = 2, ne1
c           phchi is single precision, so use tmp variables
            curr = phchi (i, ipath)
            old = phchi (i-1, ipath)
            call pijump (curr, old)
            phchi (i, ipath) = curr
  490    continue

         do 500  ik = 1, nkx
            call terp1 (xk, achi(1,ipath),  ne1, xk0(ik), achi0)
            call terp1 (xk, phchi(1,ipath), ne1, xk0(ik), phchi0)
            ccpath(ik) = 
     1        achi0 * exp (coni * (2 * xk0(ik) * reff(ipath) + phchi0))
c           note that this already includes s02, deg, sig2, etc.
c           sum total complex chi
            cchi(ik) = cchi(ik) + ccpath(ik)
  500    continue
         nused = nused + 1

         if (iabs.eq.nabs) then
c           Put path into chi.dat, xmu.dat as required
            if (abs(sig2u(ilist)) .gt. 0.000001)  then
              write(3,515)  coment, index(ipath), sig2*(bohr**2),
     1          crit(ipath), deg(ipath), nleg(ipath), reff(ipath)*bohr, 
     2          sig2u(ilist)
              write(8,515)  coment, index(ipath), sig2*(bohr**2),
     1          crit(ipath), deg(ipath), nleg(ipath), reff(ipath)*bohr, 
     2            sig2u(ilist)
            else
              write(3,515) coment, index(ipath), sig2*(bohr**2),
     1          crit(ipath), deg(ipath), nleg(ipath), reff(ipath)*bohr
              write(8,515) coment, index(ipath), sig2*(bohr**2),
     1          crit(ipath), deg(ipath), nleg(ipath), reff(ipath)*bohr
            endif
  515       format(a2, 1x, i10, 5x, f9.5, 2f10.2, i6, f9.4, f9.5)
         endif

c        write out a chinnnn.dat for this path, if necessary.
         if (ipr4 .eq. 2 .and. iabs.eq.nabs .and. ispec.eq.0)  then
c           make filename chipnnnn.dat
            write(fname,520)  index(ipath)
  520       format('chip', i4.4, '.dat')
            open (unit=9, file=fname, status='unknown',iostat=ios)
            call chopen (ios, fname, 'ff2chi')
            do 530  ihead = 1, nhead
               lhead = istrln(head(ihead))
               if (lhead .gt. 0)  then
                  write(9,320) head(ihead)(1:lhead)
  320             format (a)
               endif
  530       continue
            if (dwcorr)  then
               write(9,340)  s02, tk, thetad, sig2g
  340          format (' S02', f7.3, '  Temp', f8.2,'  Debye temp',f8.2,
     1        '  Global sig2', f9.5)
            else
               write(9,341)  s02, sig2g
  341          format (' S02', f7.3,
     1      '                                        Global sig2', f9.5)
            endif
            if (alphat .gt. zero)  then
               write(9,321)  alphat
  321          format (' 1st and 3rd cumulants, alphat = ', 1pe20.4)
            endif

            if (abs(vrcorr).ge.eps4 .or. abs(vicorr).ge.eps4)  then
               write(9,342)  vrcorr, vicorr
  342          format (' Energy zero shift, vr, vi ', 1p, 2e14.5)
            endif
            write(9,*) 'Debye-waller factor ', sig2, sig3

            write(9,610)
  610       format (1x, 71('-'))
            write(9,535)
  535       format ('       k         chi           mag          ',
     1              'phase        phase-2kr  @#')
            do 540  i = 1, nkx
               ckp = sqrt (xkp(i)*abs(xkp(i)) + xkref)
c              it would be better to use interpolation for ckp
c              fix later if complaints about chipnnn.dat files, ala
               xlam0 =  - dimag(ckp)
               ccc = ccpath(i) * exp(2 * reff(ipath) * xlam0)
               phase = 0
               if (abs(ccc) .gt. 0)  phase=atan2(dimag(ccc), dble(ccc))
               if (i .gt. 1)  call pijump (phase, phase0)
               phase0 = phase
               write(9,630)  xkp(i)/bohr, dimag(ccc), abs(ccc), phase,
     1                       phase-2*xk0(i)*reff(ipath)
  630          format (1x, f10.4, 3x, 4(1pe13.6,1x))
  540       continue
            close (unit=9)
         endif

  550    continue
  560 continue

c     close files opened for sigem and sigrem
      if (idwopt.eq.1) then
        close (unit=iem)
      elseif (idwopt.eq.2) then
        close (unit=irm1)
        close (unit=irm2)
      endif
      if (alphat .gt. 0) then
        close (unit=icum)
      endif

      return
      end
      subroutine ff2xmu (ispec, ipr4, idwopt, critcw, s02, sig2g,
     1                   tk, thetad, mbconv, absolu,  !KJ added absolu 3-06
     1                   vrcorr, vicorr, alphat, thetae, iabs, nabs,
     4            elnes,ipmin,ipmax,ipstep)   !KJ added this line 1-06     
c     adds the contributions from each path and absorber, including
c     Debye-Waller factors. Writes down main output: chi.dat and xmu.dat
      implicit double precision (a-h, o-z)

c={../HEADERS/const.h
      parameter (pi = 3.14159 26535 89793 23846 26433d0)
      parameter (one = 1, zero = 0)
      parameter (third = one/3)
      parameter (raddeg = 180 / pi)
      complex*16 coni
      parameter (coni = (0,1))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919 158 292 677 512 811d0)

      parameter (bohr = 0.529 177 249d0, ryd  = 13.605 698d0)
      parameter (hart = 2 * ryd)
      parameter (alpinv = 137.035 989 56d0)
c     fine structure alpha
      parameter (alphfs = 1 / alpinv)
c= ../HEADERS/const.h}
c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}
      parameter (eps4 = 1.0e-4)
      integer ipmin,ipmax,ipstep,elnes !KJ my variables 1-06 
      integer absolu !KJ 3-06     

c     header from list.dat
      dimension lhead(nheadx)
      character*80  head(nheadx)
      complex*16 gtr(nex),gtrful(ipmin:ipmax,nex) !KJ added gtrful 1-06

      parameter (npx=15000)
!KJ      parameter (npx = 1200)
c     indices of paths to do, read from list.dat
      dimension ip(npx)
      real sig2u(npx)

      complex*16 cchi(nex), ckp
c     to keep Im part of cchi 11.18.97 ala
      dimension rchtot(nex), xkp(nex)
      complex*16 chia(nex)

      logical dwcorr
      character*512 slog
      character*2 coment
      parameter (coment='# ')

c     Stuff from feff.bin, note that floating point numbers are
c     single precision.  Be careful throughout this routine, especially
c     when passing things to subroutines or intrinsic functions.
      real rnrmav, xmu, edge
      character*80 title(nheadx), titfms
      character*6  potlbl(0:nphx)
      dimension iz(0:nphx)
c     central atom phase shift at l0
      complex phc(nex)
      complex ck(nex)
      real xk(nex)
      dimension index(npx)
      dimension nleg(npx)
      real deg(npx), reff(npx), crit(npx)
      dimension ipot(legtot,npx)
      real rat(3,legtot,npx), beta(legtot,npx), eta(legtot,npx)
      real ri(legtot,npx)
      real achi(nex,npx), phchi(nex,npx)

c     stuff from xsect.bin
      complex*16 emxs(nex), xsec(nex)
      dimension omega(nex), xkxs(nex), xsnorm(nex)
      
c !KJ locals  1-06
      integer iip,nip
      logical cross 
      character*9 f1,f2
      character*10 f0,f3
      complex*16 gtrtemp(nex*(1+ipmax-ipmin))
      complex*16 kxsec(nex)     
c !KJ end my variables


c     get gtr - result of FMS
      do 112 ie =1,nex
      gtr(ie)= 0
      do 112 iip =ipmin,ipmax  !KJ I added this variable ip 1-06
  112 gtrful(iip,ie) = 0
      ntfms = 0
      nip=ipmax-ipmin+1 !KJ 1-06
      
      open (unit=1, file='fms.bin', status='old', iostat=ios)
      if (ios.le.0) then
         ntfms = 1
         read(1, 113) titfms
  113    format(a)
         read(1, 115) ne, ne1, ne3, nph, npadx
  115    format(5(1x,i3))
         call rdpadx(1, npadx, gtrtemp, ne*nip)  !KJ I added *nip, changed gtr to gtrtemp  1-06
      endif
      close (unit=1)
!KJ Next lines my addition to read several spectra at once. 1-06
	i=1
	do iip=ipmin,ipmax
	do j=1,ne
        gtrful(iip,j)=gtrtemp(i+j-1)
	enddo
	  i=i+ne
	enddo
!KJ Now we don't need gtrtemp anymore.  End my changes



c     read xsect.bin file
      call  rdxbin (s02p, erelax, wp, edgep, s02, gamach, ne1, ik0,
     2  emxs, omega, xkxs, xsnorm, xsec, nxsec, mbconv, title, ntitle)


c !KJ loop over iip added to process several spectra at once  1-06
c !KJ reading of feff.bin moved inside the loop (used to be before reading
c !KJ xsect.bin      
      do iip=ipmin,ipmax,ipstep
        cross=(.not.(iip.eq.1.or.iip.eq.10.or.iip.eq.5.or.iip.eq.9))
      
c !KJ choose different filename for each spectrum.
        if(iip.eq.1) then
	  f1(1:9)='chi.dat  '
	  f2(1:9)='xmu.dat  '
	  f0(1:10)='feff.bin  '
	  f3(1:10)='list.dat  ' 	  
	elseif(iip.eq.10) then
	  f1(1:9)='chi10.dat'
	  f2(1:9)='xmu10.dat'
	  f0(1:10)='feff10.bin'
	  f3(1:10)='list10.dat'	  
	elseif(iip.gt.1.and.iip.lt.10) then
	  f1(1:4)='chi0'
	  f1(5:5)= char(48+iip)
	  f1(6:9)='.dat'
	  f2(1:4)='xmu0'
	  f2(5:5)= char(48+iip)
	  f2(6:9)='.dat'
	  f0(1:5)='feff0'
	  f0(6:6)= char(48+iip)
	  f0(7:10)='.bin'	
	  f3(1:5)='list0'
	  f3(6:6)= char(48+iip)
	  f3(7:10)='.dat'
	else
	  stop 'crazy iip in ff2xmu'
	endif
	do i=1,nex
	gtr(i)=gtrful(iip,i)
	enddo
	


c     open list.dat and read list of paths we want
      open (unit=1, file= f3, status='old', iostat=ios)!KJ changed 'list.dat' to f3 4-06
      ntotal = 0
      if (ios.le.0) then
        call chopen (ios, f3, 'ff2afs')  !KJ id.
        nhead = nheadx
        call rdhead (1, nhead, head, lhead)
c       skip a label line
        read(1,*)
c       ip is index of path, sig2u is debye-waller from user
        do 100  i = 1, npx
           read(1,*,end=110)  ip(i), sig2u(i)
           ntotal = i
  100   continue
  110   continue
      endif
      close (unit=1)


       call rdfbin (f0, nphx, nex, npx, legtot, !KJ changed 'feff.bin' to f0  1-06
     $     nptot, ne, npot, ihole, iorder, ilinit, 
     $     rnrmav, xmu, edge, potlbl, iz, phc, ck, xk,
     $     index, nleg, deg, reff,
     $     crit, ipot, rat, beta, eta, ri, achi, phchi)

c     make combined title
      if (ntfms.eq.1) then
        ntitle = ntitle + 1
        title(ntitle) = titfms
      endif
      do 120 ihead = 1, nhead
 120  title(ntitle+ihead) = head(ihead)
      ntitle = ntitle + nhead

c     write feffnnnn.dat
      if (ipr4.eq.3) then
         call feffdt(ntotal,ip,nptot,ntitle,title,ne,npot,
     $        ihole, iorder, ilinit, rnrmav, xmu, edge, potlbl,
     $        iz,phc,ck,xk,index,
     $        nleg,deg,nepts,reff,crit,ipot,rat,achi,phchi)
       end if

c     If there is a vicorr, will need a mean free path factor xlam0.
c     Use it as  chi(ie) * exp (2 * reff * xlam0)
c     ckp is ck' = ck prime.
      if (abs(vicorr) .ge. eps4) then
         do 180  ie = 1, ne
            ckp = sqrt (ck(ie)**2 + coni*2*vicorr)
            xlam0 = aimag(ck(ie)) - dimag(ckp)
            do 170  ipath = 1, nptot
               achi(ie,ipath) = achi(ie,ipath) * 
     1               exp (2 * reff(ipath) * xlam0)
  170       continue
  180    continue
      endif

c     k'**2 = k**2 + vr. If there is no real correction
c     (vrcorr = 0), these two grids will be the same.
c           k' is value for output,  k is  value used for
c           interpolations with original grid.

c     vrcorr shifts the edge and the k grid
      if (abs(vrcorr) .gt. eps4)  then
         edge = edge - vrcorr
      endif

c     ik0 is index at fermi level
      do 250  i = 1, ne
         temp = xk(i)*abs(xk(i)) + 2*vrcorr
         if (temp.ge. 0) then
           xkp(i) = sqrt(temp)
         else
           xkp(i) = - sqrt(-temp)
         endif
  250 continue
     

      dwcorr = .false.
      if (tk .gt. 1.0e-3)  dwcorr = .true.

c     Open chi.dat and xmu.dat (output) and start headers
      if (iabs.eq.nabs) then
         open (unit=3, file=f1, status='unknown', iostat=ios) !KJ changed chi.dat to f1
         call chopen (ios, f1, 'ff2chi') !KJ id.  1-06
         open (unit=8, file=f2, status='unknown', iostat=ios) !KJ changed xmu.dat to f2
         call chopen (ios, f2, 'ff2chi') !KJ id.

c        write miscellaneous staff into headers  !KJ corrected typo
         call wrhead (8, ntitle, title, dwcorr, s02,
     1     tk, thetad, sig2g, alphat, vrcorr, vicorr, critcw)
     
c        also write information on the screen
         if (alphat .gt. zero)  then
            write(slog,322) alphat
  322       format ('    1st and 3rd cumulants, alphat = ', 1pe20.4)
            call wlog(slog)
         endif
         if (abs(vrcorr).ge.eps4 .or. abs(vicorr).ge.eps4)  then
            write(slog,343) vrcorr*hart, vicorr*hart
  343       format ('    Energy zero shift, vr, vi ', 1p, 2e14.5)
            call wlog(slog)
         endif

         write(slog,370) critcw
         call wlog(slog)
  370    format ('    Use all paths with cw amplitude ratio', f7.2, '%')
         if (dwcorr)  then
            write(slog,380) s02, tk, thetad, sig2g
            call wlog(slog)
         else
            write(slog,381) s02, sig2g
            call wlog(slog)
         endif
  380    format('    S02', f7.3, '  Temp', f8.2, '  Debye temp', f8.2,
     1           '  Global sig2', f9.5)
  381    format('    S02', f7.3, '  Global sig2', f9.5)
      endif


c     make chi and sum it
      do 400  i = 1, nex
         cchi(i) = 0
  400 continue
      do 402  ik = 1, ne
         cchi(ik)= s02 * gtr(ik)
  402 continue


c     add Debye-Waller factors
      call dwadd (ntotal, nptot, idwopt, ip, index, crit, critcw, sig2g,
     1  sig2u, dwcorr, rnrmav, nleg, deg, reff, iz, ipot, rat,tk,thetad,
     2  alphat, thetae, mbconv, s02, ne1, ck, achi, phchi, ne, xk, xkp,
     3  xkp, cchi, iabs, nabs, ispec, ipr4, ntitle,
     4  title, vrcorr, vicorr,  nused)

c     read or initialize chia - result of configuration average
      if (iabs.eq.1) then
         do 635 ie =1, nex
            chia(ie) = 0
  635    continue
      else
         open (unit=1, file='chia.bin', status='old',
     1   access='sequential', form='unformatted', iostat=ios)
         do 640 ie = 1,ne
  640    read(1) chia(ie)
         close (unit=1, status='delete')
      endif

      if(iabs.eq.1) then
c        compare grids in xsect.bin and feff.bin
         do 680 i = 1, nxsec
           del = xk(i)**2 - xkxs(i)**2
           if (abs(del) .gt.  10*eps4)  then
             call wlog(' Emesh in feff.bin and xsect.bin different.')
             call par_stop('FF2XMU-1') 
           endif
  680    continue
      endif

c     add contribution from an absorber iabs 
c     present scheme assumes that xsec is the same for all iabs.
      do 701 ik = 1, ne
         chia(ik)   = chia(ik)   + cchi(ik)/ nabs
  701 continue
      if (iabs.lt.nabs) then
c        save chia in chia.bin for averaging
         open (unit=1, file='chia.bin', status='unknown',
     1   access='sequential', form='unformatted', iostat=ios)
         do 760 ie=1,ne
  760    write(1) chia(ie)
         close(unit=1)
      endif

      if (iabs.eq.nabs) then
c        The loop over absorbers is finished. Write out the results.
         write(8,600)  coment, nused, ntotal
  600    format ( a2, 1x, i4, '/', i4, ' paths used')
  610    format ( a2, 1x, 71('-'))

         do 702 ik = 1, ne
            rchtot(ik) = dimag (chia(ik))
  702    continue
c        prepare the output grid omega
         efermi = edge + omega(1) - dble(emxs(1))

c        do convolution with excitation spectrum
         if (mbconv .gt. 0) then
            wp = wp / 2.
            call  exconv
     1      (omega, ne1, efermi, s02p, erelax, wp, xsnorm)
            call  exconv
     1      (omega, ne1, efermi, s02p, erelax, wp, rchtot)
         endif

c        normalize to xsec at 50 eV above edge
c        and prepare the output energy grid omega
         edg50 = efermi + 50 / hart
         if (ispec.eq.2) edg50 = efermi
         call terp (omega, xsnorm,  ne1, 1, edg50, xsedge)
         if (absolu.eq.1) xsedge=dble(1)  !KJ 3-06 don't normalize	 
         write(8,660)  coment, xsedge 
  660    format (a2, ' xsedge+ 50, used to normalize mu ', 1pe20.4)
         write(8,610) coment 
         write(8,670) coment
  670    format (a2,' omega    e    k    mu    mu0     chi     @#')


         do i=1,nex
	 if (.not.cross) then !KJ I added this block 1-06
	   kxsec(i)=xsec(i)
	 else
	   kxsec(i)=dcmplx(0,0)
	 endif !KJ end my code
	 enddo
	 
	
c        do correction using brouder method
         vi0 = 0
         call xscorr(ispec,emxs, ne1, ne, ik0, kxsec,xsnorm,chia,
     1       vrcorr, vi0, cchi) !KJ changed xsec to kxsec  1-06
         do 850 ie=1,ne1
           rchtot(ie)=dimag( kxsec(ie)+xsnorm(ie)*chia(ie)+cchi(ie)) !KJ id.
  850    continue
  
  
         do 855 ie=1,ne
           chia(ie) = 0
  855    continue
         call xscorr(ispec, emxs, ne1, ne, ik0, kxsec,xsnorm,chia,
     1       vrcorr, vi0, cchi) !KJ changed xsec to kxsec  1-06
         do 856 ie = 1, ne1
 856     cchi(ie) = dimag(kxsec(ie)+cchi(ie)) * coni+rchtot(ie) !KJ id.

         if (vicorr.gt.eps4 .and. ntotal.eq.0) then
c           add correction due to vicorr
            call conv(omega,cchi,ne1,vicorr)
c           call conv(omega,xsec,ne1,vicorr)
         endif


         do 860 ie = 1, ne1
            em0 = dble(emxs(ie))
            xsec0 = dimag(cchi(ie))
            rchtot(ie) = dble (cchi(ie))
            chi0  = (rchtot(ie) - xsec0)/xsedge
            write(8,700)  omega(ie)*hart, em0*hart, xkp(ie)/bohr,
     1              rchtot(ie)/xsedge, xsec0/xsedge, chi0

c           if you want f'' at the output in el. units use next line
c    1          rchtot(ie)*omega(ie)*prefac, xsec0*omega(ie)*prefac, chi0
c   with        prefac = alpinv / 4 / pi /bohr**2

  700       format (1x, 2f11.3, f8.3, 1p, 3e13.5)
  860    continue

         close (unit=8)
         close (unit=3, status='delete')
      endif
c     for if (iabs=abs); or the last absorber

      enddo !KJ of my iip=ipmin,ipmax,ipstep loop  1-06
      
      
      return
      end
      subroutine fprime( ei, emxs ,ne1, ne3, ne, ik0, xsec, xsnorm,chia,
     1                  vrcorr, vicorr, cchi)
c     calculate f' including solid state and lifetime effects.
c     using algorithm in Ankudinov, Rehr danes paper.
c     the output correction is returned via cchi. The rest is input
c      mu(omega) = xsec + xsnorm*chia  + (cchi)

      implicit double precision (a-h, o-z)
c={../HEADERS/const.h
      parameter (pi = 3.14159 26535 89793 23846 26433d0)
      parameter (one = 1, zero = 0)
      parameter (third = one/3)
      parameter (raddeg = 180 / pi)
      complex*16 coni
      parameter (coni = (0,1))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919 158 292 677 512 811d0)

      parameter (bohr = 0.529 177 249d0, ryd  = 13.605 698d0)
      parameter (hart = 2 * ryd)
      parameter (alpinv = 137.035 989 56d0)
c     fine structure alpha
      parameter (alphfs = 1 / alpinv)
c= ../HEADERS/const.h}
c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}

      dimension  xsnorm(nex), omega(nex)
      complex*16 emxs(nex), xsec(nex), chia(nex), cchi(nex) 
      complex*16 xmu(nex), aa, bb, c1, x1, x2, ec, temp
      complex*16 xmup(nex)
      dimension emp(nex)
      parameter (eps4 = 1.0d-4)
      complex*16 lorenz, funlog, value
      external lorenz, funlog
      dimension dout(7,nex)
      character*72 string
      dimension oscstr(14), enosc(14)
      integer ient
      data ient /0/

c     read data from fpf0.dat
      open (unit=16, file='fpf0.dat', status='old', iostat=ios)
      read  (16,*)  string
      read  (16,*)  eatom
      read  (16,*)  nosc
      do 5 i=1, nosc
        read (16,*) oscstr(i), enosc(i)
   5  continue
c     the rest is f0(Q) and is not currently needed
      close (unit=16)

      ient = ient+1
      ifp = 1
      efermi = dble(emxs(ne1+1)) 
      xloss = dimag(emxs(1))
      ne2 = ne-ne1-ne3
      if (ne2.gt.0) then
c        DANES
         do 10 ie = 1,ne1
   10    xmu(ie) = coni*xsnorm(ie) +  xsnorm(ie)*chia(ie)
         do 11 ie = ne1+1,ne1+ne2
   11    xmu (ie) = xsnorm(ie)*chia(ie)
         do 12 ie = ne-ne3+1, ne
   12    xmu (ie) =  coni*xsnorm(ie)
      else
c        FPRIME
         do 13 ie = 1,ne
   13    xmu (ie) = xsec(ie) + xsnorm(ie)*chia(ie)
      endif

      if (abs(vrcorr).gt.eps4) then
         bb = xmu(ik0)
         efermi = efermi - vrcorr
         do 20 ie = 1,ne1
   20    omega(ie) = dble(emxs(ie))
         call terpc(omega, xmu ,ne1, 1, efermi, bb)
         do 30 ie = 1, ne2
   30    emxs(ne1+ie) = emxs(ne1+ie) - vrcorr
         if (abs(xmu(ik0)).gt. eps4) bb = bb/xmu(ik0)
c        rescale values on vertical axis
         do 60 ie = ne1+1, ne-ne3
   60    xmu(ie) = xmu (ie) * bb 
      endif
              

      if (vicorr.gt.eps4) then
         xloss = xloss + vicorr
         do 40 ie=1,ne2
   40    omega(ie) = dimag(emxs(ne1+ie))
         call terpc(omega, xmu(ne1+1) ,ne2, 1, xloss, aa)
         do 50 ie = 1, ne1
            xx = vicorr**2 /(vicorr**2 + (dble(emxs(ie))-efermi)**2)
            xmu(ie) = xmu(ie)*(1.0d0 - xx) + aa * xx
            emxs(ie) = emxs(ie) + coni*vicorr
   50    continue
      endif

      do 200 ie = 1, ne1
c        cycle over energy points on horizontal grid

         dout(1,ie) = dble(emxs(ie)) * hart
         dele = dble(emxs(ie)) - efermi
c        delp correspond to pole with negative frequency
c        see Sakurai for details

         delp = -dele - 2*ei
c        delp = dele
c        dele = delp

         cchi(ie) = 0
         if (ne2.gt.0) then
            if (abs(dele).lt.eps4) dele = 0.0d0
            w1 = dimag(emxs(ne1+1))
            w2 = dimag(emxs(ne1+2))
            w3 = dimag(emxs(ne1+3))

c           matsubara pole
            temp = lorenz(ifp,xloss,w1,dele)*xmu(ne1+1)*2*coni*w1
            temp = temp + lorenz(ifp,xloss,w1,delp)*xmu(ne1+1)*2*coni*w1
            dout(2,ie)=dble(temp)
c           sommerfeld correction
            temp = coni*w1**2/ 6*(lorenz(ifp,xloss,w3,dele)*xmu(ne1+3)-
     2      lorenz(ifp,xloss,w2,dele)*xmu(ne1+2)) / (w3-w2) 
            dout(3,ie)=dble(temp)

            cchi(ie) = lorenz(ifp,xloss,w1,dele)*xmu(ne1+1) *2*coni*w1
     1      + coni * w1**2 / 6 * (lorenz(ifp,xloss,w3,dele)*xmu(ne1+3)-
     2      lorenz(ifp,xloss,w2,dele)*xmu(ne1+2)) / (w3-w2) 
c           from negative pole has additional minus sign
            cchi(ie) = cchi(ie) + 
     1      lorenz(ifp,xloss,w1,delp)*xmu(ne1+1) *2*coni*w1
     1      + coni * w1**2 / 6 * (lorenz(ifp,xloss,w3,delp)*xmu(ne1+3)-
     2      lorenz(ifp,xloss,w2,delp)*xmu(ne1+2)) / (w3-w2) 

c           theta funcion contribution only for positive pole
            if (dele .lt. eps4)    cchi(ie) = cchi(ie) - xmu(ie)
            if (abs(dele).lt.eps4) cchi(ie) = cchi(ie) + xmu(ie)/2

c           anomalous contribution
            temp = 0
            wp = 2*ei
            if (dele.ge.eps4) temp = xmu(ie)
            if (abs(dele).lt.eps4) temp = xmu(ie)/2
            temp = temp + xmu(ik0)*  funlog(1,xloss,wp,dele)
c               xmu(iko) + xsec(ik0)  if n3 >0
            dout(4,ie)=dble(temp) 

c           integration over vertical axis to final point
            n1 = ne1+2
            n2 = ne-ne3
            call fpint (emxs, xmu, n1, n2, dele, xloss, eps4, efermi,
     1                  value)
            cchi(ie) = cchi(ie) + value
c           add contribution from other pole
            call fpint (emxs, xmu, n1, n2, delp, xloss, eps4, efermi,
     1                  value)
            cchi(ie) = cchi(ie) + value
         endif 

c        integration over horizontal axis to final point
         temp = 0
         if (ne2.gt.0) then
c           DANES
            n1 = ne1-ik0 + 1
            do 120 i = ik0, ne1
              emp(i-ik0+1) = dble(emxs(i))
              xmup(i-ik0+1) = coni*xsnorm(i)
  120       continue
            do 130 i = 1, ne3
              emp(i+n1) = dble(emxs(i+ne-ne3))
              xmup(i+n1) = xmu(i+ne-ne3)
  130       continue
            n2 = n1 + ne3
         else
c           FPRIME
            n1 = 0
            do 140 i = 1, ne1
              if (n1.eq.0 .and. dble(emxs(i)).gt. dble(emxs(ne1+1))) 
     1            n1 = i
  140       continue
            do 150 i = 1, ne3
               emp(i) =  dble(emxs(ne1+i))
               xmup(i) =  xmu(ne1+i)
  150       continue
            n2 = ne3
         endif
         call fpintp (emp, xmup , n2, dele, xloss, efermi, value)
         temp  = temp + value
c        add contribution from other pole
         call fpintp (emp, xmup , n2, delp, xloss, efermi, value)
         temp  = temp + value

c         was used before
cc          contribution to fp from poles of the core states
c           temp=0
c           do 110  i=2, nosc
cc             eif = E_f- E_i  in hartrees
cc             eif = enosc(i)-enosc(1) 
cc             deltaf = deltaf - oscstr(i)*2*alpinv**2/eif
c              temp = temp + alpinv**2 * oscstr(i)* (dele -
c    1      enosc(i)+efermi-1)/ ((dele-enosc(i)+efermi-1)**2+xloss**2)
c              temp = temp + alpinv**2 * oscstr(i)* (delp -
c    1      enosc(i)+efermi-1)/ ((delp-enosc(i)+efermi-1)**2+xloss**2)
c 110       continue

         dout(5,ie) = dble(temp)
         cchi(ie) = cchi(ie) + temp

c        total contribution (not normalized)
         temp = xmu(ie) + cchi(ie)
         dout(6,ie) = dble(temp)
c        (integral w2 to wmax) minus (cusp formula)
         dout (7,ie) = dout(6,ie)-dout(4,ie)
  200 continue

c     restore the input energy mesh
      if (vicorr.gt.eps4) then
         do 250 ie = 1, ne1
  250    emxs(ie) = emxs(ie) - coni*vicorr
      endif
      if (abs(vrcorr).gt.eps4) then
         do 260 ie = 1, ne2
  260    emxs(ne1+ie) = emxs(ne1+ie) + vrcorr
      endif

c     if (ient.eq.1) then
      open(unit=3,file='danes.dat', status='unknown', iostat=ios)
      write(3,310) '# E  matsub. sommerf. anomal. tale, total, differ.'
  310 format (a)
      do 300 ie = 1, ne1
         write(3,320) (dout(i,ie), i=1,7)
  320    format ( 7(1x,1pe11.4))
  300 continue
      close(unit=3)
c     endif

      return
      end

      complex*16 function funlog (icase, xloss, w, dele)
c     anomalous fp should have all main features of total fp
c     except smooth difference 
c     analytic expression for anomalous fp (without integral)
c     is obtained by adding and subtracting G(Ef + i*Gamma) / E-w
c     and performing integral for Im axis analytically
c     icase = 1 simplified expression (compared to 2) 
c     icase=2  use real w 
c     icase=3  pure imaginary w (absolute value is input)
      implicit double precision (a-h, o-z)
c={../HEADERS/const.h
      parameter (pi = 3.14159 26535 89793 23846 26433d0)
      parameter (one = 1, zero = 0)
      parameter (third = one/3)
      parameter (raddeg = 180 / pi)
      complex*16 coni
      parameter (coni = (0,1))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919 158 292 677 512 811d0)

      parameter (bohr = 0.529 177 249d0, ryd  = 13.605 698d0)
      parameter (hart = 2 * ryd)
      parameter (alpinv = 137.035 989 56d0)
c     fine structure alpha
      parameter (alphfs = 1 / alpinv)
c= ../HEADERS/const.h}
      parameter (eps4 = 1.0d-4)

      if (icase.eq.1) then 
         if (abs(dele).ge.eps4) then 
            funlog= coni/2/pi*
     1      (log((-xloss+coni*dele)/w)+ log((xloss+coni*dele)/w))

         else
            funlog= coni/pi*log(abs(xloss/w))
         endif

      elseif (icase.eq.2) then
        if (abs(dele).ge.eps4) then
          funlog= coni/2/pi* (w+coni*xloss) * (
     1    ( log((-xloss+coni*dele)/w)) / (w+dele+coni*xloss) +
     2    ( log(( xloss+coni*dele)/w)) / (w+dele-coni*xloss))
        else
          funlog= coni/pi*(log(abs(xloss/w)))*
     1    (1 + coni*xloss/(w-coni*xloss))
        endif

      elseif (icase.eq.3) then
        if (abs(dele).ge.eps4) then
          funlog= -(w+xloss)/2/pi* (
     1    log((-xloss+coni*dele)/w) / (dele+coni*(w+xloss)) +
     2    log(( xloss+coni*dele)/w) / (dele+coni*(w-xloss)) )
        else
          funlog= coni/pi* log(abs(xloss/w))*
     1    (1 + xloss/(w-xloss))
        endif
      
      endif

      return
      end

      subroutine fpint (emxs, xmu, n1, n2, dele, xloss, eps4, efermi,
     1                  value)
c     performs integral for fp calculations between points n1 and n2.
      implicit double precision (a-h, o-z)
c={../HEADERS/const.h
      parameter (pi = 3.14159 26535 89793 23846 26433d0)
      parameter (one = 1, zero = 0)
      parameter (third = one/3)
      parameter (raddeg = 180 / pi)
      complex*16 coni
      parameter (coni = (0,1))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919 158 292 677 512 811d0)

      parameter (bohr = 0.529 177 249d0, ryd  = 13.605 698d0)
      parameter (hart = 2 * ryd)
      parameter (alpinv = 137.035 989 56d0)
c     fine structure alpha
      parameter (alphfs = 1 / alpinv)
c= ../HEADERS/const.h}
c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}

      complex*16 emxs(nex), xmu(nex), value
      complex*16  z1, z2, aa, bb, c1

c     last interval - similar to Matsubara pole ( shift and - sign)
c     notice that this also works for horizontal axis if last value
c     is small
      z1 = emxs(n2)-efermi
      z2 = emxs(n2-1)-efermi
      value =  - coni/pi * (z1-dele) / (xloss**2+(z1-dele)**2)
     1          *xmu(n2) * (2 * (z1-z2))
c     all other intervals
      do  300 i = n1, n2-2
         z1 = emxs(i) - efermi
         z2 = emxs(i+1) - efermi
         bb=(xmu(i+1)*(z2-dele) - xmu(i)*(z1-dele)) / xloss / (z2-z1)
         aa = xmu(i)*(z1-dele)/xloss - bb * z1
         c1 = (aa+bb*(dele+coni*xloss )) / 2 /coni
         if (abs(dele-dble(z1)).lt.eps4 .and.
     1       abs(dele-dble(z2)).lt.eps4) then
            value = value  -  coni/pi *c1*
     1      log( abs((z2-dele-coni*xloss)/(z1-dele-coni*xloss)) )
         else
            value    = value   -  coni/pi *c1*
     1      log((z2-dele-coni*xloss)/(z1-dele-coni*xloss))
         endif
         c1 = -(aa+bb*(dele-coni*xloss )) / 2 /coni
         value    = value    -  coni/pi *c1*
     1   log((z2-dele+coni*xloss)/(z1-dele+coni*xloss))
  300  continue

      return
      end

      subroutine fpintp (em, xmu, n2, dele, xloss, efermi, value)
c     performs integral for fp calculations between points 1 and n2.
c     and adds tail to infinity
      implicit double precision (a-h, o-z)
c={../HEADERS/const.h
      parameter (pi = 3.14159 26535 89793 23846 26433d0)
      parameter (one = 1, zero = 0)
      parameter (third = one/3)
      parameter (raddeg = 180 / pi)
      complex*16 coni
      parameter (coni = (0,1))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919 158 292 677 512 811d0)

      parameter (bohr = 0.529 177 249d0, ryd  = 13.605 698d0)
      parameter (hart = 2 * ryd)
      parameter (alpinv = 137.035 989 56d0)
c     fine structure alpha
      parameter (alphfs = 1 / alpinv)
c= ../HEADERS/const.h}
c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}

      dimension em(nex)
      complex*16 xmu(nex), value
      complex*16  z1, z2, aa, bb, cc

      value = 0
c     all intervals 
      do  300 i = 1, n2-1
         x1 = em(i) - efermi
         x2 = em(i+1) - efermi
         de = (x2-x1)/2
         x0 = (em(i) + em(i+1)) / 2
         call terpc(em, xmu, n2, 3, x0, aa)
         bb=(xmu(i+1) - xmu(i)) / (x2-x1)
         cc = (xmu(i+1) - aa - bb * de) / de**2
         z1 =  dele - x0 + efermi - coni*xloss
         z2 =  dele - x0 + efermi + coni*xloss
         value    = value  + 2*de*bb + 2*z1*de*cc +
     1    log((de-z1)/(-de-z1)) * (aa+bb*z1+cc*z1**2)
         value    = value  + 2*de*bb + 2*z2*de*cc +
     1    log((de-z2)/(-de-z2)) * (aa+bb*z2+cc*z2**2)
  300 continue

c     tail of xmu to infinity approximated by aa/(w-bb)**2
      x1 = em(n2-1)
      x2 = em(n2)
      a = sqrt ( dble(xmu(n2-1)/xmu(n2)) )
      b = ( a*x1 - x2) / (a-1)
      if (b.gt. x1) b = 0
      aa = xmu(n2) * (x2-b)**2
      z1 = dele -coni*xloss - b
      z2 = dele +coni*xloss - b
      x0 = x2 - b
      value = value + log( x0/(x0-z1) ) *aa/z1**2 - aa/z1/x0
      value = value + log( x0/(x0-z2) ) *aa/z2**2 - aa/z2/x0

c     multiply by constant factor
      value = - coni /2 /pi *value

      return
      end
       subroutine rdfbin(fbfile, nphx, nex, npathx, nlegx, 
     $      npaths, ne, npot, ihole, iorder, ilinit, 
     $      rnrmav, xmu, edge,  potlbl, iz, phc, ck, xk, index, 
     $      nleg, deg, reff, crit, ipot, 
     $      rat, beta, eta, ri, achi, phchi)
c
c read path information from PAD-format feff.bin
c  arguments:
c   fbile   name of feff.bin file                             [in]
c   nphx    dimension of  potlbl, iz (both (0:nphx)           [in]
c             max # of potentials
c   nex     dimension of many energy arrays                   [in]
c             max # of energy points
c   npathx     dimension of index,nleg,ipot,deg,reff,crit,    [in]
c           rat,beta, eta, ri, achi, phchi
c             max # of paths
c   nlegx  dimension of  ipot, rat, beta,eta,ri              [in]
c            max # of legs in a path
c   npaths  number of paths read                             [out]
c   ne      maximum number of energy points read             [out]
c   npot    number of potentials read                        [out]
c   rnrmav  average norman radius                            [out]
c   edge    shift in edge energy (?)                         [out]
c   iorder  order of genfmt matrix used                      [out]
c   potlbl  array of potential labels                        [out]
c   iz      array of atomic numbers for potentials           [out]
c   phc     array of central atom phase shift (complex)      [out]
c   ck      array of wavenumbers/momentum (complex)          [out]
c   xk      array of wavenumbers/momentum (real)             [out]
c   index   array of path indices                            [out]
c   nleg    array:  number of legs in path                   [out]
c   deg     array:  path degeneracy                          [out]
c   reff    array:  half path length of path                 [out]
c   crit    array:  importance factor for path               [out]
c   ipot    array:  pots, in order, that make up the path    [out]
c   rat     array:  atomic positions of atoms in path        [out]
c   beta    array:  euler angle for path                     [out]
c   eta     array:  second euler angle for path              [out]
c   ri      array:  path leg distances for path              [out]
c   achi    array:  amplitude of chi for path                [out]
c   phchi   array:  phase of chi for path                    [out]
c
c notes:
c   the data in feff.bin is written completely in printable
c   ascii characters.  The file is however, highly formatted
c   and kept fairly small.  all text is stored as is, but most
c   numerical data in arrays (both real and complex) is stored
c   in a special Packed Ascii Data (PAD) format which uses 6
c   printable characters to represent a real number.
c
c   special markers in the first 1 or 2 characters of each line
c   give hints about the contents of the line:
c      #_    top 2 lines.  The first line must begin "#_feff.bin"
c      #"    title lines / plain text
c      #&    misc info about potentials, calc method
c      #@    potential labels and iz
c      ##    path index,  deg, reff, crit, ipots involved
c      !     PAD characters to be read as a real array
c      $     PAD characters to be read as a complex array
c
c copyright (c) 1999  matt newville:  jan 1999
c modified by alex ankudinov: feb 2000; few fixes for feff8.2
c------------------------------------------------------------------
       integer nphx, nex, npathx, nlegx, npaths
       integer i, j, ivers, nexmax
       character*(*) fbfile
       character*(*) potlbl(0:nphx), filnam*128, str*128, msg*256
       integer iorder,  index(npathx), nleg(npathx)
       integer ne, npot, ipot(nlegx,npathx), iz(0:nphx)
       integer istrln, ier1, ier2, ier3, nwords, npadx, nwordx
       real    deg(npathx), reff(npathx), crit(npathx)
       real    rnrmav, edge, xk(nex)
       double precision bohr, tmpdp
       parameter (bohr = 0.529 177 249d0, nwordx = 20)
       parameter(nexmax = 256)
       character*20 words(nwordx)
       real     rat(3,nlegx,npathx), beta(nlegx,npathx)
       real     eta(nlegx,npathx),  ri(nlegx,npathx)
       real     achi(nex,npathx), phchi(nex,npathx), tmpr(nexmax)
       complex  phc(nex), ck(nex), tmpc(nexmax)
       external  istrln

c open feff.bin
       filnam = ' '
       filnam = fbfile
       call triml(filnam)
       il     = istrln(filnam)
       if (filnam.eq.' ') filnam = 'feff.bin'
cc       print*, ' RDFBIN!  ', filnam(1:il),':',il
       open (unit=3, file=filnam, status='old', err=450)
 10    format(a)
c first line, must match  "#_feff.bin"
       read(3,10,err=920) str
       call triml(str)
       if ((str(1:10).ne.'#_feff.bin')) go to 900
c check version of feff.bin : only support version 3 here!!
       ivers = 0
       if ((str(1:14).eq.'#_feff.bin fil')) ivers = 1
       if ((str(1:14).eq.'#_feff.bin v02')) ivers = 2
       if ((str(1:14).eq.'#_feff.bin v03')) ivers = 3
       if (ivers.ne.3) go to 930
c second line:   npot, ne
       read(3,10,err=920) str
       call triml(str)
       if ((str(1:2).ne.'#_')) go to 900
       nwords = 3
       str    = str(3:)
       call bwords(str,nwords,words)
       if (nwords.ne.3) go to 905
       call str2in(words(1), npot,  ier1)
       call str2in(words(2), ne,    ier2)
       call str2in(words(3), npadx, ier3)
       if ((ier1.ne.0).or.(ier2.ne.0).or.(ier3.ne.0)) go to 910

c  read in misc stuff:  (rnrmav, edge, iorder )
       read(3,10,err=920) str
       call triml(str)
       if (str(1:2).ne.'#&') go to 900
       nwords = 6
       str    = str(3:)
       call bwords(str,nwords,words)
       if (nwords.ne.6) go to 905
       call str2in(words(1), ihole,  ier1)
       call str2in(words(2), iorder, ier2)
       call str2in(words(3), ilinit, ier3)
       if ((ier1.ne.0).or.(ier2.ne.0).or.(ier3.ne.0)) go to 910
       call str2re(words(4), rnrmav, ier1)
       call str2re(words(5), xmu   , ier2)
       call str2re(words(6), edge  , ier3)
       if ((ier1.ne.0).or.(ier2.ne.0).or.(ier3.ne.0)) go to 910
c  read pot label and iz line
       read(3,10,err=920) str
       call triml(str)
       if (str(1:2).ne.'#@') go to 900
       nwords = 2 * npot + 2
c  note: potlbl cannot be blank!!
       str    = str(3:)
       call bwords(str, nwords, words)
       if (nwords.ne.(2 + 2*npot)) go to 905
       do 200 i = 0, npot
          potlbl(i) = words(i+1)
          iz(i) = -1
          call str2in(words(2+npot+i),iz(i),ier1)
          if (ier1.ne.0)  go to 910
 200   continue

c read  numerical data that are the same for all paths
       call rdpadc(3,npadx, phc, ne)
       call rdpadc(3,npadx, ck,ne)
       call rdpadr(3,npadx, xk,ne)
       npaths = 0
c now, for each path:
       do 300  i = 1, npathx
          index(i) = 0
c  read path  info "##" line  and retrieve all the stuff from it
          read(3,10,end=450,err=920) str
          call triml(str)
          if (str(1:2).ne.'##')   go to 900
          nwords = nwordx
          str    = str(3:)
          call bwords(str,nwords,words)
          call str2in(words(1),  index(i), ier1)
          call str2in(words(2),  nleg(i), ier2)
          call str2re(words(3),  deg(i), ier3)
          if ((ier1.ne.0).or.(ier2.ne.0).or.(ier3.ne.0)) go to 910
          call str2dp(words(4),  tmpdp, ier2)
          reff(i) = tmpdp / bohr
          call str2dp(words(5),  tmpdp, ier3)
          crit(i) = tmpdp
          if ((ier1.ne.0).or.(ier2.ne.0).or.(ier3.ne.0)) go to 910
          npaths = npaths + 1
          do 230 j = 1, nleg(i)
             call str2in(words(5+j),ipot(j,i),ier1)
             if (ier1.ne.0) go to 910
 230      continue
c
c  next, read padded arrays for rat,beta, ..., achi, phchi
          call rdpadr(3,npadx, rat(1,1,i),3*nleg(i))
          call rdpadr(3,npadx, beta(1,i),   nleg(i))
          call rdpadr(3,npadx, eta(1,i),    nleg(i))
          call rdpadr(3,npadx, ri(1,i),     nleg(i))
          call rdpadr(3,npadx, achi(1, i),  ne)
          call rdpadr(3,npadx, phchi(1, i), ne)
c  fill in rest of achi and phchi with zeros
          do 270 j = ne+1, nex
             achi(j,i)  = 0
             phchi(j,i) = 0
 270      continue
 300    continue
 450    continue
       close(3)
cc       print*, ' RDFBIN done!'
       return
 900   call wlog (' -- rdfbin error: wrong format : at line')
       go to 990
 905   call wlog (' -- rdfbin error: missing data : at line')
       go to 990
 910   call wlog (' -- rdfbin error:   bad data   : at line')
       go to 990
 920   call wlog (' -- rdfbin error: unknown error: at line')
       go to 990
 930   call wlog (' -- rdfbin error: unknown version of feff.bin')
       go to 990

 990   call wlog (str)
       call par_stop(' -- fatal error reading feff.bin -- ')
       end

      subroutine reff2x(mchi, ispec, ipr6, idwopt, critcw, s02, sig2g,
     1                   tk, thetad, mbconv, absolu,  !KJ added absolu 3-06
     1                   vrcorr, vicorr, alphat, thetae, iabs, nabs,
     4            elnes,ipmin,ipmax,ipstep)   !KJ added this line  1-06     

      implicit double precision (a-h, o-z)

c={../HEADERS/const.h
      parameter (pi = 3.14159 26535 89793 23846 26433d0)
      parameter (one = 1, zero = 0)
      parameter (third = one/3)
      parameter (raddeg = 180 / pi)
      complex*16 coni
      parameter (coni = (0,1))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919 158 292 677 512 811d0)

      parameter (bohr = 0.529 177 249d0, ryd  = 13.605 698d0)
      parameter (hart = 2 * ryd)
      parameter (alpinv = 137.035 989 56d0)
c     fine structure alpha
      parameter (alphfs = 1 / alpinv)
c= ../HEADERS/const.h}
c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}

cc    global.dat 
c       configuration average
        integer nabs, iphabs
c       global polarization data
        integer  ipol, ispin, le2
        double precision evec(3), xivec(3), spvec(3), elpty,angks,rclabs
        complex*16 ptz(-1:1, -1:1)
cc    mod6.inp
        integer  mchi, idwopt, ipr6, mbconv, absolu  !KJ added absolu 3-06
        double precision  vrcorr, vicorr, s02, tk, thetad
        double precision  alphat, thetae, sig2g

        integer elnes,ipmin,ipmax,ipstep  !KJ my variables  1-06
c     Local stuff
      character*512 slog

c     standard formats for string, integers and real numbers
  10  format(a)
  20  format (20i4)
  30  format (6f13.5)

cc    read global.inp
      open (file='global.dat', unit=3, status='unknown',iostat=ios)
c       configuration average data
        read  (3, 10) slog
        read  (3, 45) nabs, iphabs, rclabs
  45    format ( 2i8, f13.5)
      close(3)
c     read mod6.inp
      open (file='mod6.inp', unit=3, status='old',iostat=ios)
        read (3,10)  slog
        read (3,20)  mchi, ispec, idwopt, ipr6, mbconv, absolu !KJ added absolu 3-06
        read (3,10)  slog
        read (3,30)  vrcorr, vicorr, s02, critcw
        read (3,10)  slog
        read (3,30)  tk, thetad, alphat, thetae, sig2g
      close(3)

      
c  !KJ Next section added to read ELNES variables       1-06
c     read eels.inp
      elnes=0
      open(file='eels.inp',unit=3,status='old',err=900)
        read(3,*,err=900,end=900) 
	read(3,20,err=900,end=900) elnes
        read(3,*,err=900,end=900)
        read(3,*,err=900,end=900)
	read(3,*,err=900,end=900)
	read(3,20,err=900,end=900) ipmin,ipstep,ipmax
      close(3)
      goto 901
900   continue
      elnes=0
901   continue
      if(elnes.eq.0) then
        ipstep=1
	ipmax=1
	ipmin=1
      endif
               
c  !KJ end my changes      

c     transform energies to atomic units
      vrcorr = vrcorr / hart
      vicorr = vicorr / hart

      return
      end
      subroutine xscorr(ispec, emxs ,ne1, ne, ik0, xsec, xsnorm, chia,
     1                  vrcorr, vicorr, cchi)
c     convolute xmu(E)=xsec+xsnorm*chia with lorentzian using 
c     calculations in the complex energy plane

c     Input: ispec - type of spectroscopy
c       emxs - complex energy grid
c       ne1 - number of points on horizonatal axis 
c       ne - total number of points (ne-ne1) points on vertical axis
c       ik0 - Fermi level index on horizontal axis
c       xsec, xsnorm, chia - give function f in complex energy plain
c           xmu(ie) = xsec + xsnorm*chia
c       vrcorr = correction for the shift of the Fermi level
c       vicorr = 0 (disabled)
c     Output: cchi(w) - result of convolution for w = dble(emxs)
c       cchi(w) = \int_C dE xmu(E)*xloss/pi/((E-w)**2+xloss**2) = 
c       xmu(w+i*xloss)* [1/2+atan(w-efermi/xloss)/pi] +
c       \int_C dE ff(E)*xloss/pi/((E-w)**2+xloss**2)
c       where ff(E)=xmu(E)-xmu(w+i*xloss) for w<efermi we use 
c       xmu(efermi+i*xloss) instead of xmu(w+i*xloss);
c       contour C starts at efermi, goes vertically to efermi+i*xloss 
c       and then goes horizontally to infinity + i*xloss

      implicit double precision (a-h, o-z)
c={../HEADERS/const.h
      parameter (pi = 3.14159 26535 89793 23846 26433d0)
      parameter (one = 1, zero = 0)
      parameter (third = one/3)
      parameter (raddeg = 180 / pi)
      complex*16 coni
      parameter (coni = (0,1))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919 158 292 677 512 811d0)

      parameter (bohr = 0.529 177 249d0, ryd  = 13.605 698d0)
      parameter (hart = 2 * ryd)
      parameter (alpinv = 137.035 989 56d0)
c     fine structure alpha
      parameter (alphfs = 1 / alpinv)
c= ../HEADERS/const.h}
c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}

      dimension  xsnorm(nex), omega(nex)
      complex*16 emxs(nex), xsec(nex), chia(nex), cchi(nex) 
      complex*16 xmu(nex), aa, bb, c1, f1, f2, ff(nex), xmu0
      parameter (eps4 = 1.0d-4)
      complex*16 ec(nex), fc(nex), e1,e2, z1,z2, corr
      complex*16 lorenz
      external lorenz, astep

      ne2 = ne-ne1
      efermi = dble(emxs(ne)) 
      xloss = dimag(emxs(1))

c     xmu - analytic function in complex energy plain
      do  ie = 1,ne
        xmu (ie) = xsec(ie) + xsnorm(ie)*chia(ie)
      enddo
c     real frequencies
      do ie = 1, ne1
        omega(ie) = dble(emxs(ie))
      enddo

      if (abs(vrcorr).gt.eps4) then
c       account for the fermi level shift
        bb = xmu(ik0)
        efermi = efermi - vrcorr
        call terpc(omega, xmu ,ne1, 1, efermi, bb)

c       shift the vertical axis
        do ie = 1, ne2
          emxs(ne1+ie) = emxs(ne1+ie) - vrcorr
        enddo

c       rescale values on vertical axis
        bb = bb/xmu(ik0)
        do ie = ne1+1, ne
          xmu(ie) = xmu (ie) * bb 
        enddo
      else
        bb = 1
      endif

c     construct the integration countur C
      nc = 0
c     start with points on vertical axis below xloss
      do ie = 1,ne2
        if (dimag(emxs(ne1+ie)).lt.xloss) then
          nc = nc+1
          ec(nc) = emxs(ne1+ie)
          fc(nc) = xmu(ne1+ie)
        endif
      enddo
c     add corner at efermi + xloss*i
      nc = nc+1
      ic0 = nc
      if (abs(vrcorr).gt.eps4) then
        ec(nc) = efermi + coni*xloss
        fc(nc) = bb * xmu(ik0)
      else
        ec(nc) = emxs(ik0)
        fc(nc) = xmu(ik0)
      endif
c     add points on horizontal axis above efermi
      if (ispec.ne.2) then
        do ie = 1,ne1
          if (dble(emxs(ie))-efermi.gt.eps4) then
            nc = nc+1
            ec(nc) = emxs(ie)
            fc(nc) = xmu(ie)
          endif
        enddo
      else
c       ispec=2 - emission calculations- need points below E_fermi
        do ie = ne1,1,-1
          if (efermi-dble(emxs(ie)).gt.eps4) then
            nc = nc+1
            ec(nc) = emxs(ie)
            fc(nc) = xmu(ie)
          endif
        enddo
      endif
c     endo of countour construction
              
c     cycle over frequency points 
      do ie = 1, ne1
        if (omega(ie).ge.efermi) then
          xmu0 = xmu(ie)
          if (ispec.eq.2) xmu0 = xmu(ik0)*bb
        else
          xmu0 = xmu(ik0)*bb
          if (ispec.eq.2) xmu0 = xmu(ie)
        endif
        e1 = omega(ie) + coni*xloss
        e2 = omega(ie) - coni*xloss
        do ic = 1, nc 
          ff(ic) = fc(ic) - xmu0
        enddo
        dele = omega(ie) - efermi
        cchi(ie) = xmu0 * astep( xloss, dele)
        if (ispec.eq.2) cchi(ie) = xmu0 - cchi(ie)
        corr = 0

        if (abs(dele).lt.eps4) dele = 0.0d0
        w1 = dimag(ec(1))
        w2 = dimag(ec(2))
        w3 = dimag(ec(3))
        ip =0

c       add half matsubara pole contribution
c       equivalent to integral from efermi to efermi+i*w1
        corr = corr + lorenz(ip,xloss,w1,dele)*ff(1) *coni*w1
        if (nc0.gt.3) then
c       add sommerfeld correction (correction for derivative)
c         corr = corr + coni * w1**2 / 6   / (w3-w2) *
c    2   (lorenz(ip,xloss,w3,dele)*ff(3)-lorenz(ip,xloss,w2,dele)*ff(2))
        endif


c       cycle over contour points 
        do ic = 1,nc-1
c         perform integration over contour from efermi+i*2*w1 to 
c         efermi+i*xloss; linear interpolation of ff between  z1 and z2
          z1 = ec(ic)
          z2 = ec(ic+1)
c         if (ic.eq.1) z1 = efermi+coni*2*w1
          f1 = ff(ic)
          f2 = ff(ic+1)
c         if (ic.eq.1) f1 = (f1*(z2-z1) + f2*(z1-ec(ic))) / (z2-ec(ic))
c         add correction from pole above real axis
          aa = 0
          if (abs(z1-e1).gt.eps4 .and. abs(z2-e1).gt.eps4) then
            aa = log((z2-e1)/(z1-e1)) *(f1*(z2-e1)+f2*(e1-z1))
c           z1 or z2 equal to e1; in this case corr is exactly zero
          endif
c         second pole 
          aa = aa - log((z2-e2)/(z1-e2)) *(f1*(z2-e2)+f2*(e2-z1))
          corr = corr + aa/ (z2-z1) /2/pi/coni
        enddo
c       end of cycle over contour points
        if (ispec.eq.2) corr = -corr
c       if (ispec.eq.2) corr = 0

        cchi(ie) = cchi(ie) +  corr
c       return the result of convolution minus bare value
        cchi(ie) = cchi(ie) - xmu(ie)
      enddo
c     end of cycle over frequency points

c     restore the input energy mesh
      if (abs(vrcorr).gt.eps4) then
        do  ie = ne1+1, ne
          emxs(ie) = emxs(ie) + vrcorr
        enddo
      endif

      return
      end

      complex*16 function lorenz (ifp, xloss, w, dele)
      implicit double precision (a-h, o-z)
c={../HEADERS/const.h
      parameter (pi = 3.14159 26535 89793 23846 26433d0)
      parameter (one = 1, zero = 0)
      parameter (third = one/3)
      parameter (raddeg = 180 / pi)
      complex*16 coni
      parameter (coni = (0,1))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919 158 292 677 512 811d0)

      parameter (bohr = 0.529 177 249d0, ryd  = 13.605 698d0)
      parameter (hart = 2 * ryd)
      parameter (alpinv = 137.035 989 56d0)
c     fine structure alpha
      parameter (alphfs = 1 / alpinv)
c= ../HEADERS/const.h}
c     ifp is dummy now. correspond to ifp=0 in old code
c     can remove it and change calls to lorenz in other routines

      lorenz = xloss /pi / (xloss**2+(coni*w-dele)**2)

      return
      end

      double precision function astep ( xloss, dele)
      implicit double precision (a-h, o-z)
c={../HEADERS/const.h
      parameter (pi = 3.14159 26535 89793 23846 26433d0)
      parameter (one = 1, zero = 0)
      parameter (third = one/3)
      parameter (raddeg = 180 / pi)
      complex*16 coni
      parameter (coni = (0,1))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919 158 292 677 512 811d0)

      parameter (bohr = 0.529 177 249d0, ryd  = 13.605 698d0)
      parameter (hart = 2 * ryd)
      parameter (alpinv = 137.035 989 56d0)
c     fine structure alpha
      parameter (alphfs = 1 / alpinv)
c= ../HEADERS/const.h}

      astep = 0.5d0 + atan(dele/xloss) /pi
      if (astep.lt.0.d0) astep = 0.d0
      if (astep.gt.1.d0) astep = 1.d0

      return
      end
c///////////////////////////////////////////////////////////////////////
c FEFF PROGRAMS (referred below as a System)
c Copyright (c) 1986-2002, University of Washington.
c 
c END-USER LICENSE 
c 
c A signed End-user License Agreement from the University of Washington
c Office of Technology Transfer is required to use these programs and
c subroutines.
c 
c See the URL: http://leonardo.phys.washington.edu/feff/
c 
c USE RESTRICTIONS:
c 
c 1. The End-user agrees that neither the System, nor any of its
c components shall be used as the basis of a commercial product, and
c that the System shall not be rewritten or otherwise adapted to
c circumvent the need for obtaining additional license rights.
c Components of the System subject to other license agreements are
c excluded from this restriction.
c
c 2. Modification of the System is permitted, e.g., to facilitate
c its performance by the End-user. Use of the System or any of its
c components for any purpose other than that specified in this Agreement
c requires prior approval in writing from the University of Washington.
c
c 3. The license granted hereunder and the licensed System may not be
c assigned, sublicensed, or otherwise transferred by the End-user.  
c
c 4. The End-user shall take reasonable precautions to ensure that
c neither the System nor its components are copied, or transferred out
c side of his/her current academic or government affiliated laboratory
c or disclosed to parties other than the End-user.
c 
c 5. In no event shall the End-user install or provide this System
c on any computer system on which the End-user purchases or sells
c computer-related services.
c 
c 6. Nothing in this agreement shall be construed as conferring rights
c to use in advertising, publicity, or otherwise any trademark or the
c names of the System or the UW.   In published accounts of the use or
c application of FEFF the System should be referred to  by this name,
c with an appropriate literature reference:
c 
c FEFF8: A.L. Ankudinov, B. Ravel, J.J. Rehr, and S.D. Conradson,
c        Phys. Rev. B 58, pp. 7565-7576 (1998).
c
c LIMITATION OF LIABILITY:
c
c 1.   THE UW MAKES NO WARRANTIES , EITHER EXPRESSED OR IMPLIED, AS TO
c THE CONDITION OF THE SYSTEM, ITS MERCHANTABILITY, OR ITS FITNESS FOR
c ANY PARTICULAR PURPOSE.  THE END-USER AGREES TO ACCEPT THE SYSTEM
c 'AS IS' AND IT IS UNDERSTOOD THAT THE UW IS NOT OBLIGATED TO PROVIDE
c MAINTENANCE, IMPROVEMENTS, DEBUGGING OR SUPPORT OF ANY KIND.
c
c 2. THE UW SHALL NOT BE LIABLE FOR ANY DIRECT, INDIRECT, SPECIAL,
c INCIDENTAL OR CONSEQUENTIAL DAMAGES SUFFERED BY THE END-USER OR ANY
c OTHER PARTIES FROM THE USE OF THE SYSTEM.
c
c 3.  The End-user agrees to indemnify the UW for liability resulting
c from the use of the System by End-user. The End-user and the UW each
c agree to hold the other harmless for their own negligence.
c
c TITLE:
c
c 1.  Title patent, copyright and trademark rights to the System are
c retained by the UW. The End-user shall take all reasonable precautions
c to preserve these rights.
c 
c 2.  The UW reserves the right to license or grant any other rights to
c the System to other persons or entities.
c///////////////////////////////////////////////////////////////////////
c License is applicable for routines below, until otherwise specified.
c
c---------------------------------------------------------------------
c     program sigms.f
c
c     calculates debye-waller factors for each multiple
c     scattering path using Debye-Model correlations
c
c     files:  input  pathd_all.dat  multiple scattering path data
c             output fort.3  sig**2 vs path
c                    fort.2  long output
c
c     version 1  (29 july 91)
c
c     coded by j. rehr
c     path data from s. zabinsky
c
c     modified to use pdata.inp, Dec 1991, siz
c     Subroutine version, Dec 1991, siz
c
c---------------------------------------------------------------------

      subroutine sigms (tk, thetad, rs, nlegx, nleg, rat, iz, sig2)
c               tk temperature in degrees K
c               thetad debye temp in degrees K
c               rs=wigner seitz or norman radius in bohr, averaged
c                  over entire problem
c                  (4pi/3)*rs**3 = sum( (4pi/3)rnrm**3 ) / N
c                  (sum is over all atoms in the problem)
c               nlegx used in dimensions of rat and iz
c               nleg nlegs in path
c               rat positions of each atom in path
c               iz atomic number of each atom in path
c               NB Units of distance in this routine
c                  are angstroms, including sig**2.  rs is in bohr.
c               sig2 is output debye waller factor

      implicit double precision (a-h,o-z)

c={../HEADERS/const.h
      parameter (pi = 3.14159 26535 89793 23846 26433d0)
      parameter (one = 1, zero = 0)
      parameter (third = one/3)
      parameter (raddeg = 180 / pi)
      complex*16 coni
      parameter (coni = (0,1))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919 158 292 677 512 811d0)

      parameter (bohr = 0.529 177 249d0, ryd  = 13.605 698d0)
      parameter (hart = 2 * ryd)
      parameter (alpinv = 137.035 989 56d0)
c     fine structure alpha
      parameter (alphfs = 1 / alpinv)
c= ../HEADERS/const.h}

c     nlegx is max number of atoms in any one path
      dimension rat(3,0:nlegx)
      dimension iz(0:nlegx)
c#mn
       external dist

c      parameters
c               x = k_d*R   (distance parameter)
c               R distance in angstroms
c               y = hbar omegad/kT = thetad/t
c               thetad debye temp in degrees K
c               tk temperature in degrees K
c               k_d = (6*pi**2 N/V) = debye wave number
c               N/V=1/(4pi/3rs**3)
c               rs=wigner seitz or norman radius in bohr
c               ami, amj masses at sites i and j in amu
c               I = int_0^1 (y/x) dw sin(wx)coth(wy/2)

c     Note:  There are nleg atoms including the central atom
c            index 0 and index nleg both refer to central atom,
c            which makes special code unnecessary later.

      sigtot=0
      do 800 il=1,nleg
      do 800 jl=il,nleg

c        calculate r_i-r_i-1 and r_j-r_j-1

         rij = dist (rat(1,il), rat(1,jl))
         call corrfn (rij, cij, thetad, tk, iz(il), iz(jl), rs)
         sig2ij=cij

         rimjm = dist (rat(1,il-1), rat(1,jl-1))
         call corrfn (rimjm, cimjm, thetad, tk, iz(il-1), iz(jl-1), rs)
         sig2ij=sig2ij+cimjm

         rijm = dist (rat(1,il), rat(1,jl-1))
         call corrfn (rijm, cijm, thetad, tk, iz(il), iz(jl-1), rs)
         sig2ij=sig2ij-cijm

         rimj = dist (rat(1,il-1), rat(1,jl))
         call corrfn (rimj, cimj, thetad, tk, iz(il-1), iz(jl), rs)
         sig2ij=sig2ij-cimj

         riim = dist (rat(1,il), rat(1,il-1))
         rjjm = dist (rat(1,jl), rat(1,jl-1))

         ridotj=(rat(1,il)-rat(1,il-1))*(rat(1,jl)-rat(1,jl-1))+
     1          (rat(2,il)-rat(2,il-1))*(rat(2,jl)-rat(2,jl-1))+
     2          (rat(3,il)-rat(3,il-1))*(rat(3,jl)-rat(3,jl-1))
         ridotj=ridotj/(riim*rjjm)

c        double count i .ne. j  terms
         if(jl.ne.il) sig2ij=2*sig2ij
         sig2ij=sig2ij*ridotj
         sigtot=sigtot+sig2ij

  800 continue
      sig2=sigtot/4

c     sig2 is in bohr**2, just as we wanted for ff2chi
      return
      end



      subroutine corrfn(rij,cij,thetad,tk,iz1,iz2,rsavg)
c     subroutine calculates correlation function
c     c(ri,rj)=<xi xj> in the Debye approximation
c
c             =(1/N)sum_k exp(ik.(Ri-Rj))(1/sqrt(mi*mj))*
c              (hbar/2w_k)*coth(beta hbar w_k/2)
c             = (3kT/mu w_d**2)*sqrt(mu**2/mi*mj)*I
c
c      parameters
c               x = k_d*R   (distance parameter)
c               R distance in angstroms
c               y = hbar omegad/kT = thetad/t
c               thetad debye temp in degrees K
c               tk temperature in degrees K
c               k_d = (6*pi**2 N/V) = debye wave number
c               N/V=1/(4pi/3rs**3)
c               rs=wigner seitz or norman radius in bohr
c               ami, amj masses at sites i and j in amu
c               I = int_0^1 (y/x) dw sin(wx)coth(wy/2)
c
c      solution by numerical integration
c
      implicit double precision (a-h, o-z)
      common /xy/ x, yinv

c={../HEADERS/const.h
      parameter (pi = 3.14159 26535 89793 23846 26433d0)
      parameter (one = 1, zero = 0)
      parameter (third = one/3)
      parameter (raddeg = 180 / pi)
      complex*16 coni
      parameter (coni = (0,1))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919 158 292 677 512 811d0)

      parameter (bohr = 0.529 177 249d0, ryd  = 13.605 698d0)
      parameter (hart = 2 * ryd)
      parameter (alpinv = 137.035 989 56d0)
c     fine structure alpha
      parameter (alphfs = 1 / alpinv)
c= ../HEADERS/const.h}

c     con=hbar**2/kB*amu)*10**20   in ang**2 units
c     hbar = 1.054 572 666 e-34, amu = 1.660 540 e-27, 
c     kB = 1.380 6581 d-23
      parameter (con = 48.508 459 393 094)
c#mn
       external atwtd

c     external fn
c     rij=2.55
c     tk=295
c     thetad=315
c     ami=amj=63.55 at wt for Cu
c     rs=2.7

      ami=atwtd(iz1)
      amj=atwtd(iz2)
      rs=rsavg
c     thetad in degrees K, t temperature in degrees K
c     y=thetad/tk
      yinv=tk/thetad
      xkd=(9*pi/2)**(third)/(rs*bohr)
      fac=(3/2.)*con/(thetad*sqrt(ami*amj))
      rj=rij
      x=xkd*rj
c     call numerical integration
      call bingrt (grater, eps, nx)
      cij=fac*grater
      return
      end
      double precision function fn(w)
      implicit double precision (a-h,o-z)
      common/xy/x,yinv
c     fn=(sin(wx)/x)*coth(wy/2)
c     change code to allow t=0 without bombing
c     fn=2/y
      fn=2*yinv
      if(w.lt.1.e-20) return
      fac=w
      if(x.gt.0.) fac=sin(w*x)/x
      emwy=0.
      if(yinv.gt.0.0125) emwy=exp(-w/yinv)
      emwy=exp(-w/yinv)
      fn=fac*(1+emwy)/(1-emwy)
      return
      end
c-----------------------------------------------
      subroutine bingrt (b, eps, n)
c     subroutine calculates integrals between [0,1]
c      b = int_0^1 f(z) dz
c     by trapezoidal rule and binary refinement
c     (romberg integration)
c     coded by j rehr (10 Feb 92)
c     see, e.g., numerical recipes for discussion
c     and a much fancier version
c-----------------------------------------------
c     del=dz  itn=2**n tol=1.e-5
c     starting values
      implicit double precision (a-h,o-z)
      common /xy/x,yinv
      character*512 slog
c     external fn
c     error is approximately 2**(-2n) ~ 10**(-.6n)
c     so nmax=10 implies an error of 1.e-6
      parameter(nmax = 10, tol = 1.e-5)
      parameter(zero=0, one=1)
      n=0
      itn=1
      del=1.
      bn=(fn(zero)+fn(one))/2
      bo=bn
 10   continue
c     nth iteration
c     b_n+1=(b_n)/2+deln*sum_0^2**n f([2n-1]deln)
      n=n+1
      if(n.gt.nmax) go to 40
      del=del/2
      sum=0.
      do 20 i=1, itn
      zi=(2*i-1)*del
 20   sum=sum+fn(zi)
c     bnp1=b_n+1 is current value of integral
      bnp1=bn/2+del*sum
c     cancel leading error terms b=[4b-bn]/3
c     note: this is the first term in the
c     neville table - remaining errors were
c     found too small to justify the added code
      b=(4*bnp1-bn)/3
      eps=abs((b-bo)/b)
      if(eps.lt.tol) goto 60
      bn=bnp1
      bo=b
      itn=itn*2
      goto 10
 40   write(slog,50) n,itn, b,eps
      call wlog(slog)
 50   format(' not converged, n,itn,b,eps=',
     1  2i4,2e14.6)
      return
 60   continue
      return
      end
c---------------------------------------------------------------------
c     program sigem
c
c     calculate the Debye-Waller factors for each MS path
c     using the equation-of-motion methods
c
c     input files:  feff.inp and spring.inp
c
c     version 2  ( January 99)
c
c     coded by  A. Poiarkova
c
c---------------------------------------------------------------------
c  References:  
c             for the EM method: Phys. Rev. B , 59, p.948, 1999
c     also see dissertation 
c        "X-ray Absorption Fine Structure Debye-Waller Factors"
c         by Anna V. Poiarkova
c
c---------------------------------------------------------------------
c         tk temperature in degrees K
c         nleg  nlegs in path
c         rat   positions of each atom in path
c         NB Units of distance in this routine
c            are angstroms, including sig2. 
c         sig2 is output DW factor
c
c---------------------------------------------------------------------
      subroutine sigem (sig2mx, sig2x, iem, tk, ipath, nleg, rat, sig2)
      implicit double precision (a-h, o-z)

c={dwpar.h
c-*-fortran-*-
c nlegx1 MUST be the same as legtot, the maximum number of scattering
c       legs in a path
c nphx1 MUST be the same as nphx, the maximum number of atomic species

      parameter (natxdw = 200)
      parameter (nlegx1 = 9)
      parameter (nphx1=7)

c= dwpar.h}
c={../HEADERS/parallel.h
      integer par_type, this_process, numprocs, my_rank
      logical master, worker, parallel_run
      real*8 wall_comm, time_comm
      common /timing/ wall_comm, time_comm
      common /parallel/ numprocs, my_rank, this_process, 
     .          master, worker, parallel_run, par_type
c= ../HEADERS/parallel.h}

c feff parameters (from dim.h):
c     parameter (legtot=9) 
c     parameter (nphx = 7)

      parameter (nphx = nphx1)
      parameter (natx = natxdw)

c local parameters:
      parameter (amu0  = 1.660 54)
      parameter (pi = 3.14159 26535 89793 23846 26433d0)
      parameter (nwx = 700)

      double precision sig2mx, sig2x(0:nphx,0:nphx)
      dimension iphat(natx), izph(0:nphx)

c variables shared with rdspr.f:
      dimension rat1(3,natx), iz(natx)
      dimension dm(3,3,natx,natx)
      dimension rnn(3,natx,natx)
      dimension nnl(natx,natx)

c local variables:
      dimension rat(3,0:nlegx1)
      dimension nconv(0:nlegx1)
      dimension q0(3,natx)
      dimension gr(nwx), w(nwx)
      dimension nq0(0:nlegx1) 
      dimension uu(3,natx), up(3,natx), ff(3,natx)

      character*30  fname
      parameter (ntitx1 = 10)
      character*71  title(ntitx1)
      dimension ltit(ntitx1)
c     character*80  titlep(ntitx1)

      character*512 slog
      logical iem_open

      save 
      data nsigc /0/
c-------------------------------------------------------------

      inquire(unit=iem,opened=iem_open)
      if (nsigc.eq.0) then

c Read coordinates and potentials from feff.inp
      call dwrdin (rat1, iphat, izph, natom,
     1            ntitle, title, ltit)

      if (natom.gt.natx) natom=natx
      do 5 iat=1, natom
         iz(iat) = izph(iphat(iat))
         if (iphat(iat).eq.0) i0=iat
  5   continue

c Read spring.inp and build dynamical matrix
      call rdspr(rat1, iz, natom, i0, 
     1           dm, rnn, 
     1           acut, res, wmax, dosfit, zshell, w0,
     1           rintr, iprdos, nnl)

            write(slog,7)
   7        format(2x,'Calculating Debye-Waller factors via EMM...')
            call wlog(slog)
            write(slog,9)
   9        format(2x,'This might take a while.')
            call wlog(slog)

      if (ipath.ne.0.and.iem_open) then
c           Echo title cards to s2_em.dat
            do 10  i = 1, ntitle
               write(iem,12)  title(i)(1:ltit(i))
  10        continue
  12        format (1x, a)
            write(iem,17) tk, natom
  17        format(1x,'temperature =',f7.2,2x,'N_at =',i4)
            write(iem,19)
  19        format (1x, 71('-'))
            write(iem,25)
            write(slog,25)
            call wlog(slog)
  25        format(3x,'ipath',4x,'nleg',3x,'sig2',5x,
     1            'mu_ipath',2x,'check0(%)')
      endif

c Integration parameters:
      wmaxx=sqrt(zshell)
      dt=2.*pi/wmaxx/15.
c top limit in t integration:
      cutoff=2.*sqrt(2.*acut)/res/wmaxx 
      nstep=cutoff/dt
      xlam=acut/(cutoff)**2
      wl=0.0000001
c top limit in w integration:
      wm=wmax*wmaxx 
      dw=0.01 
      nw=(wm-wl)/dw + 1
      if (nw .gt. nwx) then
          nw = nwx
          dw = (wm-wl)/(nw -1)
      endif
      nfit = dosfit*nw/20.

      endif
c------------------------------------

cc    Open path input file (unit in) and read title.  Use unit 2.
c     ntitle2 = 5
c     open(unit=2,file='paths.dat',status='old', iostat=ios)
c     call chopen (ios, 'paths.dat', 'sigem')
c     call rdhead (2, ntitle2, titlep, ltit)
cc    if (ntitle2 .le. 0)  then
cc       titlep(1) = ' '
cc    endif

c 84  continue
c     read(2,*,end=1010) ipath, nleg
c     skip label (x y z ipot rleg beta eta) and read the path
c     read(2,*)
      do 78 ileg=0,nleg
c        read(2,*,end=1010) (rat(j,ileg),j=1,3)
         nconv(ileg)=0
  78  continue

      do 88 n=1,3
         aa = rat(n,nleg)
         do 87 i=0,(nleg-2)
            j=nleg-i
            rat(n,j)=rat(n,j-1)
  87     continue
         rat(n,1)=aa
  88  continue
      do 89 i=1,nleg
         nq0(i)=0.
  89  continue

c nconv converts # of an atom in the nleg list of coordinates (rat) to
c its # in the full list of all atomic coordinates (rat1)
      do 94 i=1,natom
         do 91 n=1,3
   91    q0(n,i)=0.
         do 95 jl=1,nleg
            m=0
            do 93 n=1,3
               l=nint(100.*rat(n,jl))
               l1=nint(100.*rat1(n,i))
               if (abs(l-l1).le.1) m=m+1
   93       continue
            if (m.eq.3) then
              nconv(jl)=i
            endif
   95    continue
   94 continue
c     check that all path atoms are found
      do 96 jl=1,nleg
        if (nconv(jl).eq.0) then
           print*,' did not find atom jl=', jl
           print*, rat(1,jl),rat(2,jl),rat(3,jl)
           call par_stop('SIGREM-1')
        endif
  96  continue

      atmu=0.
      iq0=0
      nconv(0)=nconv(nleg)
      do 100 il=1,nleg
         l=nconv(il)
         do 101 jq=1,iq0
  101    if(nq0(jq).eq.l) go to 102
         iq0=iq0+1
         nq0(iq0)=l
  102    continue
         nq0x=iq0
         i=nconv(il)
         im=nconv(il-1)
         ip=nconv(il+1)
c        if (il.eq.1) im=nconv(nleg)
         if (il.eq.nleg) ip=nconv(1)
         atmass=atwtd(iz(i))
      do 100 n=1,3
         atmu=atmu + 0.25*( rnn(n,i,im)+rnn(n,i,ip) )**2 /atmass
  100 continue
      atmu=1./atmu
      icount=1
  108 continue
      icount= icount+1
      if (icount.gt.10) call par_stop('SIGREM-2')

      do 115 i=1,natom
      do 115 n=1,3
  115 q0(n,i)=0.

c Build initial state vector |Q_j(0)> for the current path
      do 116 n=1,3
      do 116 il=1,nleg
         i=nconv(il)
         im=nconv(il-1)
         ip=nconv(il+1)
         if (il.eq.1) im=nconv(nleg)
         if (il.eq.nleg) ip=nconv(1)
         atmass=atwtd(iz(i))
         q0(n,i)=q0(n,i)+sqrt(atmu/atmass)*(rnn(n,im,i)-rnn(n,i,ip))/2.
  116 continue

c make sure it's normalized <Q_j(0)|Q_j(0)>=1
      q0q0=0.
      do 120 iq0=1,nq0x
      i=nq0(iq0)
      do 120 n=1,3
         q0q0=q0q0+q0(n,i)*q0(n,i)
 120  continue
      p00=nint(q0q0*1000.d0)/1000.d0
      if (abs(p00-1.d0).gt.5.d-4) then
         atmu=atmu/q0q0
         go to 108
      endif

c     to get THz units:
      wnorm=100.*w0/sqrt(amu0*10.) 
c*** moments
      a0=0.
      do 132 il=1,nq0x
      do 132 im=1,nq0x
         l=nq0(il)
         m=nq0(im)
      do 132 n1=1,3
      do 132 n2=1,3
         a0 = a0 + q0(n1,l)*dm(n1,n2,l,m)*q0(n2,m)/w0/w0
  132 continue
      a0=wnorm*sqrt(a0)

      do 125 kw=1, nwx
         gr(kw) = 0.
  125 w(kw) = (wl+(kw-1)*dw)

c  make file prdennnnn.dat
      if (master.and.ipath.ne.0.and.ipath.le.iprdos) then 
         write(fname,130)  ipath
  130    format('prden', i4.4, '.dat')
         open (unit=25, file=fname, status='unknown',iostat=ios)
         call chopen (ios, fname, 'sigem')
         do 134  i = 1, ntitle
            write(25,136)  title(i)(1:ltit(i))
  134    continue
         write(25,135) natom
  135    format('#',1x,'N_at =', i4)
  136    format ('#',1x, a)
         write(25,138)
  138    format ('#',1x, 71('-'))
      endif


c  set initial conditions
      do 150 i=1, natom
      do 150 n=1,3
         uu(n,i)=q0(n,i)
         up(n,i)=uu(n,i)
  150 continue

c Solve 3*natom equations of motion and find projected VDOS (gr)
      dt2=dt*dt
      t=dt/2.
      do 200 kstep = 1, nstep
c        damping factor:
         e1=exp(-xlam*t*t) 
         xat=0.
         do 167 i=1, natom
         do 167 n=1,3
  167    xat = xat + uu(n,i)*q0(n,i)
         xat=xat*e1
         do 170 kw=1, nw
  170    gr(kw) = gr(kw) + xat*cos(w(kw)*t)*dt
         if(kstep.eq.nstep) go to 200

         do 175 i=1,natom
         do 175 n=1,3
  175    ff(n,i)=0.

         do 180 i=1,natom
            jn=1
  185       if (nnl(i,jn).ne.0) then
               j=nnl(i,jn)
               am=w0*w0
               do 187 n1=1,3
               do 187 n2=1,3
                  ff(n1,i)=ff(n1,i)-dm(n1,n2,i,j)*uu(n2,j)/am
                  if(i.ne.j) ff(n1,j)=ff(n1,j)-dm(n1,n2,j,i)*uu(n2,i)/am
  187          continue
               jn = jn + 1
               go to 185
            endif
  180    continue

         do 199 i=1,natom
         do 199 n=1,3
            put=2.*uu(n,i)-up(n,i)+dt2*ff(n,i)
            up(n,i)=uu(n,i)
            uu(n,i)=put
  199    continue

  200 t=t+dt

      afit = 0.
      if (nfit.ne.0) then
         if (w(nfit).ne.0.) afit=gr(nfit)/(w(nfit)**4)
      endif

c fit vibr.density to A*w^4, for low w
      do 225 kw=1, nfit
         gr(kw)=afit*w(kw)**4
  225 continue

c Normalization of the pr.density of modes 
c (it's the 2/pi factor which was left out till now with,
c perhaps, a small diffrence due to the fit) 
      gr(nw)=0.
      if (gr(1).lt.0.) gr(1)=0.
      xx=(gr(1)+gr(nw))*dw/2.
      do 247 kw=2, (nw-1)
         if (gr(kw).lt.0.) gr(kw)=0.
  247 xx = xx + gr(kw)*dw
      cn1=1./xx

      if (master.and.ipath.ne.0.and.ipath.le.iprdos) then
c to get THz units:
         wnorm=100.*w0/sqrt(amu0*10.) 
         write(25,349) ipath, nleg
  349    format('#',2x,'ipath =',i3,2x,'nleg =',i2)
         write(25,350)
c 350    format(1h#,6x,5hw,THz,18x,6hrho(w))
  350    format(1h#,6x,'cm^-1',18x,6hrho(w))
         do 370 kw=1,nw
            write(25,360) w(kw)*wnorm*100./6./pi, gr(kw)*cn1/wnorm
c           write(25,360) w(kw)*wnorm, gr(kw)*cn1/wnorm
  360       format(2x,f10.3,15x,f10.7)
  370    continue
         close (unit=25)
      endif

      wt=tk/187.64/w0
      ccc=cn1
      check0 = abs((2./pi - cn1)/(2./pi))
      check0=check0*100.
      coef = ccc*0.5*0.2587926/atmu/w0
c integrate over w to get sig2
      cth=0.
      s2=0.
c     gr(1)=0.
      do 400 kw=2, (nw-1)
         cth = 1./tanh( w(kw)/(2.*wt) )
         s2 = s2 + coef*gr(kw)*cth*dw/w(kw)
  400 continue
      sig2 = s2

      if (ipath.ne.0.and.iem_open) then
         write(iem,473) ipath, nleg, sig2,atmu,check0
         write(slog,473) ipath, nleg, sig2,atmu,check0
         call wlog(slog)
  473    format(4x,i3,4x,i3,4x,f7.5,3x,f7.3,4x,f5.2)
      endif

      nsigc = nsigc + 1
c1000 go to 84
c1010 continue
c     close (unit=2)
      if (sig2.gt.1.0) then
        sig2 = 1.0d0
        call wlog (' WARNING: Found sig**2>1. Set sig2=1. ')
        write (slog,1011) nconv(1), nconv(2)
 1011   format('          Possible zero ferquency modes with atoms', i4,
     1   ' or', i4)
         call wlog(slog)
         call wlog('          Check springs.inp')
      endif
      if (check0.gt.5.0) then
        write (slog,*) ' WARNING: Failed check0 test:missing VDOS.',
     1  ' Reduce dosfit and/or increase wmax.'
        call wlog(slog)
      endif

c     update maximum DW factors
      if (sig2.gt.sig2mx) sig2mx=sig2
      if (sig2.gt.sig2x( iphat(nconv(1)),  iphat(nconv(2)) )) then
         sig2x( iphat(nconv(1)),  iphat(nconv(2)) ) = sig2
         sig2x( iphat(nconv(2)),  iphat(nconv(1)) ) = sig2
      endif

      return
      end
c----------------------------------------------------
      subroutine dwrdin (rat, iphat, izph, nat,
     1            ntitle, title, ltit)

c     Read feff.inp for sigem.f
c     (here we need only coordinates and potentials)

      implicit double precision (a-h, o-z)

c={dwpar.h
c-*-fortran-*-
c nlegx1 MUST be the same as legtot, the maximum number of scattering
c       legs in a path
c nphx1 MUST be the same as nphx, the maximum number of atomic species

      parameter (natxdw = 200)
      parameter (nlegx1 = 9)
      parameter (nphx1=7)

c= dwpar.h}

c feff parameters:
c     parameter (nphx = 7)

      parameter (nphx = nphx1)
      parameter (natx = natxdw)

      dimension iphat(natx)
      dimension rat(3,natx)
      dimension iatph(0:nphx)

      character*6  potlbl(0:nphx)

c     Local stuff
      character*150  line
      parameter (nwordx = 20)
      character*20 words(nwordx)

      parameter (ntitx = 10)
      character*71  title(ntitx)
      dimension ltit(ntitx)
      dimension izph(0:nphx)
      logical iscomm
      parameter (nssx = 16)

      parameter (big = 1.0e5)
      character*512 slog

   10 format (a)
   20 format (bn, i15)
   30 format (bn, f15.0)

cc    initialize things

      ntitle = 0

      nat = 0
      do 100  iat = 1, natx
         iphat(iat) = -1
  100 continue

      nph = 0
      do 110  iph = 0, nphx
         iatph(iph) = 0
         izph(iph) = 0
         potlbl(iph) = ' '
  110 continue

c     Open feff.inp, the input file we're going to read
      open (unit=1, file='feff.inp', status='old', iostat=ios)
      call chopen (ios, 'feff.inp', 'rdinp')

c     tokens  0 if not a token
c             1 if ATOM (ATOMS)
c             7 if TITL (TITLE)
c            10 if DEBY (DEBYE)
c            13 if PRIN (PRINT)
c            14 if POTE (POTENTIALS)
c            -1 if END  (end)
c     mode flag  0 ready to read a keyword card
c                1 reading atom positions
c                2 reading overlap instructions for unique pot
c                3 reading unique potential definitions

      mode = 0
  200 read(1,10,iostat=ios)  line
         if (ios .lt. 0)  line='END'
         call triml (line)
         if (iscomm(line))  goto 200
         nwords = nwordx
         call bwords (line, nwords, words)
         itok = itoken (words(1),'feff.inp')

c        process the card using current mode
  210    continue

         if (mode .eq. 0)  then
            if (itok .eq. 1)  then
c              ATOM
c              Following lines are atom postions, one per line
               mode = 1

            elseif (itok .eq. 7)  then
c              TITLE title...
               ntitle = ntitle + 1
               if (ntitle .le. ntitx)  then
                  title(ntitle) = line(6:)
                  call triml (title(ntitle))
               else
                  call wlog(' Too many title lines, title ignored')
                  call wlog(' ' // line(1:71))
               endif
               mode = 0

c           elseif (itok .eq. 10)  then
cc             DEBYE  temp debye-temp
cc                  temps in kelvin
cc                  These add to any sig2 from SIG2 card or files.dat
c              read(words(2),30,err=900)  tk
c              read(words(3),30,err=900)  thetad
c              idwopt=0
c              read(words(4),20,err=900)  idwopt
c              mode = 0
            elseif (itok .eq. 14)  then
c              POTENTIALS
c              Following lines are unique potential defs, 1 per line
               mode = 3

            elseif (itok .eq. -1)  then
cc             END
               goto 220
            else
               mode = 0
c *            write(slog,'(1x,a)') line(1:70)
c *            call wlog(slog)
c *            write(slog,'(1x,a)') words(1)
c *            call wlog(slog)
c *            write(slog,'(a,i8)') ' Token ', itok
c *            call wlog(slog)
c *            call wlog(' Keyword unrecognized.')
c *            call wlog(' See FEFF document -- some old features')
c *            call wlog(' are no longer available.')
c *            call par_stop('DWRDIN-1')
            endif
         elseif (mode .eq. 1)  then
            if (itok .ne. 0)  then
cc             We're done reading atoms.
cc             Change mode and process current card.
               mode = 0
               goto 210
            endif
            nat = nat+1
            if (nat .gt. natx)  then
               write(slog,'(1x,a,i5)') 'Too many atoms, max is ', natx
               call wlog(slog)
               write(slog,'(1x,a,i5,a)') 'Only', natx,
     1       ' atoms will be considered in the DW factor calculations.'
               call wlog(slog)
               nat = nat-1
               mode = 0
               goto 210
c              call par_stop('DWRDIN-2')
            endif
            if (nat.le.natx) then
               read(words(1),30,err=900)  rat(1,nat)
               read(words(2),30,err=900)  rat(2,nat)
               read(words(3),30,err=900)  rat(3,nat)
               read(words(4),20,err=900)  iphat(nat)
               if (iphat(nat).eq.0) iat0 = nat
            else
               mode = 0
               goto 210
            endif
         elseif (mode .eq. 3)  then
            if (itok .ne. 0)  then
cc             We're done reading unique potential definitions
cc             Change mode and process current card.
               mode = 0
               goto 210
            endif
            read(words(1),20,err=900)  iph
            if (iph .lt. 0  .or.  iph .gt. nphx)  then
               write(slog,'(a,i8)') 
     1             'Unique potentials must be between 0 and ',
     1             nphx
               call wlog(slog)
               write(slog,'(i8,a)') iph, ' not allowed'
               call wlog(slog)
               write(slog,'(1x,a)') line(1:71)
               call wlog(slog)
               call par_stop('DWRDIN-3')
            endif
            read(words(2),20,err=900)  izph(iph)
cc          No potential label if user didn't give us one
cc          Default set above is potlbl=' '
            if (nwords .ge. 3)  potlbl(iph) = words(3)
         else
            write(slog,'(a,i8)') 
     .        'DWRDIN-4: Mode unrecognized, mode ', mode
c           call wlog(slog)
            call par_stop(slog)
         endif
      goto 200
  220 continue

cc    We're done reading the input file, close it.
      close (unit=1)
            if (nat .gt. natx)  then
               write(slog,'(a,i8)') 'Too many atoms for DW calculations,
     1         max is ', natx
               call wlog(slog)
               write(slog,'(a,i8,a)') 'Only atoms up to #',natx,
     1         '  will be considered'
               call wlog(slog)
            endif

      do 250 iat = 1, nat
      do 250 i = 1,3
        if (iat.ne. iat0) rat(i,iat) = rat(i,iat) - rat(i,iat0)
 250  continue
      do 251 i = 1,3
 251  rat(i,iat0) = 0.d0

cc    Find out how many unique potentials we have
      nph = 0
      do 300  iph = nphx, 0, -1
         if (izph(iph) .gt. 0)  then
            nph = iph
            goto 301
         endif
  300 continue
  301 continue
cc    Must have central atom
      if (izph(0) .le. 0)  then
         call wlog(' No absorbing atom (unique pot 0) was defined.')
         call par_stop('DWRDIN-5')
      endif
cc    Find central atom (only 1 permitted)
      iatabs = -1
      do 400  iat = 1, nat
         if (iphat(iat) .eq. 0)  then
            if (iatabs .lt. 0)  then
               iatabs = iat
            else
               call wlog(' More than one absorbing atom (potential 0)')
               call wlog(' Only one absorbing atom allowed')
               call par_stop('DWRDIN-6')
            endif
         endif
  400 continue

cc    Then find model atoms for unique pots that have them
cc    Use atom closest to absorber for model
      do 330  iph = 0, nphx
         rabs = big
         do 320  iat = 1, nat
            if (iph .eq. iphat(iat))  then
               tmp = dist (rat(1,iat), rat(1,iatabs))
               if (tmp .lt. rabs)  then
cc                this is the closest so far
                  rabs = tmp
                  iatph(iph) = iat
               endif
            endif
  320    continue
  330 continue
cc    if iatph > 0, a model atom has been found.

      if (ntitle .le. 0)  then
         ntitle = 1
         title(1) = 'Null title'
      endif
      do 490  i = 1, ntitle
         ltit(i) = istrln (title(i))
  490 continue

      return

  900 continue
      call wlog(' Error reading input, bad line follows:')
      write(slog,'(1x,a)') line(1:71)
      call wlog(slog)
      call par_stop('DWRDIN-7 fatal error.')

c      return
      end

c----------------------------------------------------------
      subroutine rdspr(rat1, iz, natom, i0, 
     1           dm, rnn, 
     1           acut, res, wmax, dosfit, zshell, w0, 
     1           rintr, iprdos, nnl)

c     Read spring.inp for multiple scattering feff and
c     build dynamical matrix.

      implicit double precision (a-h, o-z)

c={dwpar.h
c-*-fortran-*-
c nlegx1 MUST be the same as legtot, the maximum number of scattering
c       legs in a path
c nphx1 MUST be the same as nphx, the maximum number of atomic species

      parameter (natxdw = 200)
      parameter (nlegx1 = 9)
      parameter (nphx1=7)

c= dwpar.h}
c={../HEADERS/parallel.h
      integer par_type, this_process, numprocs, my_rank
      logical master, worker, parallel_run
      real*8 wall_comm, time_comm
      common /timing/ wall_comm, time_comm
      common /parallel/ numprocs, my_rank, this_process, 
     .          master, worker, parallel_run, par_type
c= ../HEADERS/parallel.h}

c feff parameters:

c     parameter (nphx = nphx1)
      parameter (natx = natxdw)

c new local parameters:
      parameter (nangx = 7*natx)
      parameter (nsprx = 40)
      parameter (nshx = 100)

c variables shared with sigem.f:
      dimension rat1(3,natx), iz(natx)
      dimension dm(3,3,natx,natx)
      dimension rnn(3,natx,natx)
      dimension nnl(natx,natx)

c local variables:
      dimension rshell(natx,0:nshx)
      dimension nspr(2,nsprx), drij(natx,natx)
      dimension str(natx,natx)
      dimension ang(nangx), dang(nangx)
      dimension nang(3,nangx)
      dimension dmstr(3,3,natx,natx),dma(3,3,natx,natx)
      dimension si(3), sj(3), sk(3)

      character*150  line
      parameter (nwordx = 20)
      character*20  words(nwordx)
      character*512 slog

      logical iscomm

   10 format (a)
   20 format (bn, i15)
   30 format (bn, f15.0)

c initialize things

      do 40 i=1, natom
      do 40 j=1, natom
         str(i,j)=0.
         drij(i,j)=0.02
   40 continue
      do 47 na=1,nangx
         ang(na)=0.
         dang(na)=0.
      do 47 m=1,3
         nang(m,na)=0
   47 continue

      do 50 ispr=1, nsprx
      do 50 n=1,2
   50 nspr(n,ispr)=0

      acut=3.
      res=0.05
      dosfit=0.
      wmax=1.
      na=1
      nintr=0
      strx=10000.
      ispr=1
      iprdos = 0
      ddrij=0.02
      ddang=0.02

      open(unit=1,file='spring.inp',status='old', iostat=ios)
      call chopen (ios, 'spring.inp', 'rdspr')

c     tokens  0 if not a token
c             1 if STRE (STRETCHES)
c             2 if ANGL (ANGLES)
c             3 if VDOS
c             4 if PRDOS 
c            -1 if END  (end)
c     mode flag  0 ready to read a keyword card
c                1 reading stretches
c                2 reading angle-bends
 
      mode = 0
  200 read(1,10,iostat=ios)  line
         if (ios .lt. 0)  line='END'
         call triml (line)
         if (iscomm(line))  goto 200
         nwords = nwordx
         call bwords (line, nwords, words)
         itok = itoken (words(1),'spring.inp')

c        process the card using current mode
  210    continue

         if (mode .eq. 0)  then
            if (itok .eq. 1)  then
c              STRE
c              Following lines are stretches, one per line
c              read(words(2),20,err=900)  nintr
               mode = 1
            elseif (itok .eq. 2)  then
c              ANGL
c              Following are ...
               mode = 2
            elseif (itok .eq. 3)  then
c              VDOS
c              VDOS  resolution, a_cut, wmax, dosfit
c               0 - do not run modules, 1 - run module
               read(words(2),30,err=900)  res
               read(words(3),30,err=900)  wmax
               read(words(4),30,err=900)  dosfit
               if (nwords.gt.4) then
                   read(words(5),30,err=900)  acut
               endif
               mode = 0
            elseif (itok .eq. 4)  then
c              PRINT  iprdos
c              to print or not to print prdennnnn.dat files;
c              if the card is present, these files will be
c              printed for paths 1 through iprdos
               iprdos = 1 
               read(words(2),20,err=900)  iprdos
               mode = 0
            elseif (itok .eq. -1)  then
c              END
               goto 220
            else
               write(slog,'(1x,a)') line(1:71)
               call wlog(slog)
               write(slog,'(1x,a)') words(1)
               call wlog(slog)
               write(slog,'(a,i8)') ' Token ', itok
             call wlog(slog)
               call wlog(' Keyword unrecognized.')
               call wlog(' See FEFF document -- some old features')
               call wlog(' are no longer available.')
             call par_stop('RDSPR-1')
            endif
         elseif (mode .eq. 1)  then
            if (itok .ne. 0)  then
c              We're done reading stretches
c              Change mode and process current card.
               mode = 0
               goto 210
            endif
            read(words(1),20,err=900) ii
            i=ii+1
            call chekin (ii, natom, line)
            read(words(2),20,err=900) jj
            j=jj+1
            call chekin (jj, natom, line)
            read(words(3),30,err=900) str(i,j)
            if (str(i,j).lt.strx) then
               strx=str(i,j)
               ix=i
               jx=j
            endif
            nspr(1,ispr)=i
            nspr(2,ispr)=j
            ispr=ispr+1
            read(words(4),30,err=900) ddrij
            drij(i,j) = abs(ddrij)/100.
            drij(j,i) = abs(ddrij)/100.
         elseif (mode .eq. 2)  then
            if (itok .ne. 0)  then
c              We're done reading angle-bends
c              Change mode and process current card.
               mode = 0
               goto 210
            endif
            read(words(1),20,err=900) ii
            i=ii+1
            call chekin (ii, natom, line)
            read(words(2),20,err=900) jj
            j=jj+1
            call chekin (jj, natom, line)
            read(words(3),20,err=900) kk
            k=kk+1
            call chekin (kk, natom, line)
            read(words(4),30,err=900) ang(na)
            nang(1,na)=i
            nang(2,na)=j
            nang(3,na)=k
            read(words(5),30,err=900) ddang
            dang(na) = abs(ddang)/100.
            na=na+1
         else
            write(slog,'(a,i8)') 'Mode unrecognized, mode ', mode
            call wlog(slog)
            call par_stop('RDSPR-2')
         endif
      goto 200
  220 continue

c     We're done reading the input file, close it.
      close (unit=1)
      nax=na-1

c     write statistics on found bonds and angles into spring.dat 
      if (master) then
        open (unit=2, file='spring.dat', status='unknown',iostat=ios)
        call chopen (ios, 'spring.dat', 'spring')
        write(2,*) ' Statistics on spring constants in spring.inp.'
        write(2,*) '   STRETCHES  i  j  aa   found_number'  
      endif

c find all stretching bonds
      do 321 jspr=1, (ispr-1)
         icnt=0
         i=nspr(1,jspr)
         j=nspr(2,jspr)
         aa = str(i,j)
         ddrij=drij(i,j)
         ip=iz(i)
         jp=iz(j)
         rij = dist (rat1(1,i), rat1(1,j))
         if (aa.eq.0.) go to 321
         do 320 k=1, natom
            do 320 l=k+1, natom
               kp=iz(k)
               lp=iz(l)
               rkl = dist (rat1(1,k), rat1(1,l))
               comp = abs(rij/rkl - 1.)
               if (comp.gt.ddrij) go to 320
               if (ip.ne.kp.or.jp.ne.lp) then
                  if (ip.ne.lp.or.jp.ne.kp) go to 320
               endif
               str(k,l) = aa
               str(l,k) = aa
calex       to check the bonds, that were found
calex        print*, k,l,aa
                icnt = icnt+1
  320    continue
         str(j,i) = aa
         if (master) write (2,*) i-1, j-1, aa, icnt
  321 continue
      if (master) write(2,*) '   BENDS   i  j  k   aa   found_number'  

c find all bending angles
      naxx=nax
      do 323 na=1,nax
         icnt=1
         i=nang(1,na)
         j=nang(2,na)
         k=nang(3,na)
         ddrij=drij(i,j)
         ddrkj=drij(k,j)
         ip=iz(i)
         jp=iz(j)
         kp=iz(k)
         call coss(rat1(1,i),rat1(1,j),rat1(1,k),cosijk)
         rij = dist (rat1(1,i), rat1(1,j))
         rkj = dist (rat1(1,k), rat1(1,j))
         aa=ang(na)
c        print*, na, i,j,k, aa
         do 326 ii=1, natom
         do 326 jj=1, natom
            if (ii.eq.jj) go to 326
            rrij=dist (rat1(1,ii), rat1(1,jj))
            do 322 kk=ii+1, natom
               if (kk.eq.jj) go to 322
               rrkj=dist (rat1(1,kk), rat1(1,jj))
               comp1 = abs(rrij/rij - 1.)
               comp2 = abs(rrkj/rkj - 1.)
               if (comp1.gt.ddrij.or.comp2.gt.ddrkj) then
                  comp1 = abs(rrkj/rij - 1.)
                  comp2 = abs(rrij/rkj - 1.)
                  if (comp1.gt.ddrij.or.comp2.gt.ddrkj) go to 322
               endif
               iip=iz(ii)
               jjp=iz(jj)
               kkp=iz(kk)
            if (iip.ne.ip.or.jjp.ne.jp.or.kkp.ne.kp) then
               if (kkp.ne.ip.or.jjp.ne.jp.or.iip.ne.kp) go to 322
            endif
               call coss(rat1(1,ii),rat1(1,jj),rat1(1,kk),cssijk)
               if (dacos(cosijk).eq.0.) go to 322
                  comp = abs( dacos(cssijk)/dacos(cosijk) -1.)
               if (comp.ge.dang(na)) go to 322
               do 324 na1=1,naxx
                  ii1=nang(1,na1)
                  jj1=nang(2,na1)
                  kk1=nang(3,na1)
               if (ii.eq.ii1.and.jj.eq.jj1.and.kk.eq.kk1) go to 322
               if (kk.eq.ii1.and.jj.eq.jj1.and.ii.eq.kk1) go to 322
 324           continue
               naxx=naxx+1
               ang(naxx)=aa
               nang(1,naxx)=ii
               nang(2,naxx)=jj
               nang(3,naxx)=kk
calex          to check the bends, that were found
c              print*, naxx, ii,jj,kk,aa
               icnt = icnt + 1
               if (naxx.eq.nangx) goto 333
 322        continue
 326     continue
         if (master) write (2,*) i-1, j-1, k-1, aa, icnt
 323  continue
 333  continue

      if (master) close (unit=2)

      do 325 i=1, natom
      do 325 nshell=0, nshx
         rshell(i,nshell)=0.
 325  continue

c find shells
      rintr=0.
      nintr=1
      do 330 i=1, natom
      nshell=0
      do 335 j=1, natom
         if (j.eq.i) go to 332
         if (nshell.gt.nshx) go to 332
         rij = dist (rat1(1,i), rat1(1,j))
         ddrij=drij(i,j)
         ncount=0
         do 331 ish=0, nshell
            b = real(rshell(i,ish))
            dif=1.
            if (b.ne.0.) dif = abs(rij -b)/b
            if (dif.le.ddrij) ncount=ncount+1 
 331     continue
         if (ncount.eq.0) then
            nshell = nshell + 1
            if (str(i,j).ne.0.and.rij.gt.rintr) rintr=rij
            rshell(i,nshell) = rij
         endif
 332     do 335 n=1,3
         rnn(n,i,j)=0.
         do 335 m=1,3
            dmstr(n,m,i,j)=0.
            dma(n,m,i,j)=0.
            dm(n,m,i,j)=0.
 335  continue
c sort rshell into ascending numerical order
c and find maximum order of interacting neighbor nintr
      do 342 jsh=2,nshell 
         aa = rshell(i,jsh)
         do 341 ish=jsh-1,1,-1
            if(rshell(i,ish).le.aa) go to 340
            rshell(i,ish+1)=rshell(i,ish)
 341     continue
         ish=0
 340     rshell(i,ish+1) = aa
         if (aa.le.rintr.and.(ish+1).gt.nintr) nintr = ish+1
 342  continue
 330  continue

      zshell=0.
      i1=0
      do 350 i=1,natom
         do 352 in=1,natx
 352     nnl(i,in)=0
         do 350 j=i+1,natom
            dx=rat1(1,j)-rat1(1,i)
            dy=rat1(2,j)-rat1(2,i)
            dz=rat1(3,j)-rat1(3,i)
            dr=sqrt(dx*dx+dy*dy+dz*dz)
            rnn(1,i,j)=dx/dr
            rnn(2,i,j)=dy/dr
            rnn(3,i,j)=dz/dr
            rnn(1,j,i)=-rnn(1,i,j)
            rnn(2,j,i)=-rnn(2,i,j)
            rnn(3,j,i)=-rnn(3,i,j)
            rrij = abs( dr/rshell(1,1) -1.)
c           if (i.eq.1.and.rrij.le.drij(i,j)) zshell=zshell+1
            if (i.eq.i0.and.rrij.le.drij(i,j)) then
               zshell=zshell+1
               if (i1.eq.0.and.str(i,j).ne.0.) i1 = j
            endif
  350 continue

c Build dynm. matrix for angle bends 
      nan=0
      do 355 na=1,naxx
c        print*,na, naxx, nangx
         i=nang(1,na)
         j=nang(2,na)
         k=nang(3,na)
         if (i.eq.j.or.j.eq.k) go to 355
         if(ang(na).eq.0.) go to 355
         nan=nan+1
         rij = dist(rat1(1,i),rat1(1,j))
         rkj = dist(rat1(1,k),rat1(1,j))
         if (rij.gt.rintr.or.rkj.gt.rintr) go to 355
         call sang (i, j, k, rat1, si, sj, sk)
         do 357 n1=1,3
         do 357 n2=1,3
            dma(n1,n2,i,j)=dma(n1,n2,i,j)+ang(na)*si(n1)*sj(n2)
            dma(n1,n2,j,k)=dma(n1,n2,j,k)+ang(na)*sj(n1)*sk(n2)
            dma(n1,n2,i,k)=dma(n1,n2,i,k)+ang(na)*si(n1)*sk(n2)

            dma(n1,n2,i,i)=dma(n1,n2,i,i)+ang(na)*si(n1)*si(n2)
            dma(n1,n2,j,j)=dma(n1,n2,j,j)+ang(na)*sj(n1)*sj(n2)
            dma(n1,n2,k,k)=dma(n1,n2,k,k)+ang(na)*sk(n1)*sk(n2)

            dma(n2,n1,j,i)=dma(n1,n2,i,j)
            dma(n2,n1,k,j)=dma(n1,n2,j,k)
            dma(n2,n1,k,i)=dma(n1,n2,i,k)
  357    continue
  355 continue

c Build dynm. matrix for stretches 
      do 375 l=1,natom
      do 375 m=l,natom
         do 373 n1=1,3
            x2=str(l,m)*rnn(n1,l,m)
         do 373 n2=1,3
         dmi=0.
         if (l.eq.m) then
            do 377 i=1,natom
               if (real(str(i,m)).eq.0.) go to 377
               dmi = dmi + str(i,m)*rnn(n1,i,m)*rnn(n2,i,m)
 377        continue
         endif
         dmstr(n1,n2,l,m) = dmi - x2*rnn(n2,l,m)
         dmstr(n2,n1,m,l)=dmstr(n1,n2,l,m)
 373  continue
 375  continue

c Add two dynm. matrices D_str+D_ang
      lnx=0
      do 380 i=1,natom
      ami=sqrt(atwtd(iz(i)))
      in=0
      do 380 j=1,natom
      amj=sqrt(atwtd(iz(j)))
      sumdm=0.
      do 381 n1=1,3
      do 381 n2=1,3
         dmdm=dma(n1,n2,i,j)+dmstr(n1,n2,i,j)
         sumdm=sumdm+abs(dmdm)
         dm(n1,n2,i,j)=dmdm/ami/amj
 381  continue
      if (real(sumdm).ne.0.and.i.le.j) then
         in=in+1
         nnl(i,in)=j
      endif
      if (in.ge.lnx) lnx=in
 380  continue

      atmu = 1./(1./atwtd(iz(i0)) + 1./atwtd(iz(i1)))
      a0=0.
      do 450 i=1,2
      do 450 j=1,2
      if (i.eq.1) l=i0
      if (i.eq.2) l=i1
      if (j.eq.1) m=i0
      if (j.eq.2) m=i1
      do 450 n1=1,3
      do 450 n2=1,3
         fact = (-1)**(i+j)*atmu
         atmass = 1./atwtd(iz(l))/atwtd(iz(m))
         a0 = a0 + fact*sqrt(atmass)*rnn(n1,i0,i1)*
     1             dm(n1,n2,l,m)*rnn(n2,i0,i1)
  450 continue
c     effective freq. for the 1st shell:
      w0=sqrt(a0) 
      if (w0.eq.0.) then
         atmux = 1./(1./atwtd(iz(ix)) + 1./atwtd(iz(jx)))
         w0=sqrt(strx/atmux)
      endif

      return

  900 continue
      call wlog(' Error reading input, bad line follows:')
      write(slog,'(1x,a)') line(1:71)
      call wlog(slog)
      call par_stop('RDSPR fatal error.')

cc      return
      end

c---------------------------------------------------

      subroutine sang (i, j, k, rat1, si, sj, sk)

c*calculates coeficients s  (Sint=sum {si*ui}) connecting internal
c*coordinate delta_phi(ijk)=(valence ijk angle bend) with ui (atomic
c*displacements)

      implicit double precision (a-h, o-z)

c={dwpar.h
c-*-fortran-*-
c nlegx1 MUST be the same as legtot, the maximum number of scattering
c       legs in a path
c nphx1 MUST be the same as nphx, the maximum number of atomic species

      parameter (natxdw = 200)
      parameter (nlegx1 = 9)
      parameter (nphx1=7)

c= dwpar.h}

      parameter (natx = natxdw)

      dimension rat1(3,natx), rji(3), rjk(3)
      dimension eji(3), ejk(3), ej(3)
      dimension si(3), sk(3), sj(3)

      dji=0.
      djk=0.
      dik=0.
      do 905 m = 1, 3
         rji(m) = rat1(m,i) - rat1(m,j)
         rjk(m) = rat1(m,k) - rat1(m,j)
         dji = dji + rji(m)**2
         djk = djk + rjk(m)**2
         dik = dik + ( rat1(m,k) - rat1(m,i) )**2
         si(m) = 0.
         sj(m) = 0.
         sk(m) = 0.
  905 continue
      dji = sqrt(dji)
      djk = sqrt(djk)
      dik = sqrt(dik)

      dotj=0.
      do 910 m = 1, 3
         eji(m) = rji(m)/dji
         ejk(m) = rjk(m)/djk
         dotj = dotj + eji(m) * ejk(m)
  910 continue
c     ri = dji
c     rk = djk
c     rj = sqrt(dji*djk)
c     rj = dik
      rj=1.
      call vect (eji, ejk, ej, sinj)
      do 920 m = 1, 3
         si(m) = rj*(dotj * eji(m) - ejk(m))/dji/sinj
         sk(m) = rj*(dotj * ejk(m) - eji(m))/djk/sinj
         sj(m) = rj*((dji - djk * dotj)*eji(m) +
     1              (djk - dji * dotj)*ejk(m))/dji/djk/sinj
  920 continue
      return
      end

c-----------------------------------------------------------
      subroutine vect (v1, v2, v3, sin12)

c*calculates vector product v3 = [v1 x v2] and sin of the angle

      implicit double precision (a-h, o-z)

      dimension v1(3), v2(3), v3(3)

      v3(1) = v1(2)*v2(3) - v1(3)*v2(2)
      v3(2) = v1(3)*v2(1) - v1(1)*v2(3)
      v3(3) = v1(1)*v2(2) - v1(2)*v2(1)
      d1 = 0.
      d2 = 0.
      d3 = 0.
      do 990 m = 1, 3
         d1 = d1 + v1(m)**2
         d2 = d2 + v2(m)**2
         d3 = d3 + v3(m)**2
  990 continue
      sin12 = sqrt(d3/d1/d2)
      return
      end

c-----------------------------------------------------------
      subroutine coss (v1,v2,v3,cos12)
c* calculates cos between two vectors v1-v2 and v3-v2
      implicit double precision (a-h, o-z)
      dimension v1(3), v2(3), v3(3)

      vv1=0.
      vv2=0.
      scal=0.
      do 995 m=1,3
         vv1=vv1+(v1(m)-v2(m))**2
         vv2=vv1+(v3(m)-v2(m))**2
         scal=scal+(v1(m)-v2(m))*(v3(m)-v2(m))
  995 continue
      cos12=scal/vv1/vv2
      return
      end

c-----------------------------------------------------------
      subroutine chekin (i, natom, line)
      character*150  line
      character*512  slog

            if (i .gt. (natom-1) .or. i .lt. 0) then
               write(slog,'(a,i8)')
     1             'the atomic indexes must be between 0 and ',
     1             (natom - 1)
               call wlog(slog)
               write(slog,'(i8,a)') i, ' not allowed'
               call wlog(slog)
               write(slog,'(1x,a)') line(1:71)
               call wlog(slog)
               call par_stop('RDSPR')
            endif
        return
        end

c---------------------------------------------------------------------
c     program sigrm
c
c     calculate the Debye-Waller factors for each MS path
c     using the recursion method
c
c     input files:  feff.inp and spring.inp
c
c     version 2  ( January 99)
c
c     coded by  A. Poiarkova
c
c---------------------------------------------------------------------
c  References:
c             for the RM: J. Synchrotron Rad., 1999 (to bu published)
c     also see dissertation
c        "X-ray Absorption Fine Structure Debye-Waller Factors"
c         by Anna V. Poiarkova
c
c---------------------------------------------------------------------
c         tk temperature in degrees K
c         nleg  nlegs in path
c         rat   positions of each atom in path
c         NB Units of distance in this routine
c            are angstroms, including sig2.
c         sig2 is output DW factor
c
c---------------------------------------------------------------------
      subroutine sigrm (sig2mx, sig2x,ir1, ir2, tk,ipath,nleg,rat,sig2)
      implicit double precision (a-h, o-z)

c={dwpar.h
c-*-fortran-*-
c nlegx1 MUST be the same as legtot, the maximum number of scattering
c       legs in a path
c nphx1 MUST be the same as nphx, the maximum number of atomic species

      parameter (natxdw = 200)
      parameter (nlegx1 = 9)
      parameter (nphx1=7)

c= dwpar.h}

c feff parameters (from dim.h):
c     parameter (legtot=9) 
c     parameter (nphx = 7)

      parameter (nphx = nphx1)
      parameter (natx = natxdw)

c local parameters:
      parameter (amu0  = 1.660 54)
      double precision sig2mx, sig2x(0:nphx,0:nphx)
      dimension iphat(natx), izph(0:nphx)

c variables shared with rdspr.f:
      dimension rat1(3,natx), iz(natx)
      dimension dm(3,3,natx,natx)
      dimension rnn(3,natx,natx)
      dimension nnl(natx,natx)

c local variables:
      dimension rat(3,0:nlegx1)
      dimension nconv(0:nlegx1)
      dimension q0(3,natx)
c   list of atoms |0>=|Q>:
      dimension nq0(0:nlegx1)  
c   state |1>=D|Q>:
      dimension q1(3,natx)  
c   list of atoms in |1>:
      dimension nq1(natx)   

c     character*30  fname
      parameter (ntitx1 = 10)
      character*71  title(ntitx1)
      dimension ltit(ntitx1)
c     character*80  titlep(ntitx1)

      character*512 slog

      logical ir1_open, ir2_open

      save 
      data nsigc /0/
c-------------------------------------------------------------

      inquire(unit=ir1,opened=ir1_open)
      inquire(unit=ir2,opened=ir2_open)

      if (nsigc.eq.0) then
c        Read coordinates and potentials from feff.inp
         call dwrdin (rat1, iphat, izph, natom,
     1            ntitle, title, ltit)

         if (natom.gt.natx) natom=natx
         do 5 iat=1, natom
            iz(iat) = izph(iphat(iat))
            if (iphat(iat).eq.0) i0=iat
  5      continue

         write(slog,7)
   7     format(2x,'Calculating Debye-Waller factors via RM...')
         call wlog(slog)
         write(slog,9)
   9     format(2x,'This might take a while.')
         call wlog(slog)

c        Read spring.inp and build dynamical matrix
         call rdspr(rat1, iz, natom, i0, 
     1           dm, rnn, 
     1           acut, res, wmax, dosfit, zshell, w0,
     1           rintr, iprdos, nnl)

         if (ipath.ne.0) then
	    if(ir1_open) then
c             Echo title cards to s2_rm2.dat
              do 10  i = 1, ntitle
                write(ir1,12)  title(i)(1:ltit(i))
  10          continue
  12          format (1x, a)
              write(ir1,17) tk, natom
  17          format(1x,'temperature =',f7.2,2x,'N_at =',i4)
              write(ir1,19)
  19          format (1x, 71('-'))
              write(ir1,25)
              write(slog,25)
  25          format(1x,'ipath',2x,'nleg',4x,'sig2',3x,'mu_ipath',4x,
     1          'w_1',6x,'w_2',7x,'A1',5x,'A2')
              call wlog(slog)
            endif
            if (iprdos.ne.0.and.ir2_open) then
c              Echo title cards to s2_rm1.dat
               do 30  i = 1, ntitle
                  write(ir2,12)  title(i)(1:ltit(i))
  30           continue
               write(ir2,17) tk, natom
               write(ir2,19)
               write(ir2,35)
  35           format(1x,'ipath',2x,'nleg',4x,'sig2',3x,'mu_ipath',
     1              4x,'w_e')
            endif
         endif
      endif
      nsigc = nsigc + 1
c---- end of first time reading -------

cc    Open path input file (unit in) and read title.  Use unit 2.
c     ntitle2 = 5
c     open(unit=2,file='paths.dat',status='old', iostat=ios)
c     call chopen (ios, 'paths.dat', 'sigrm')
c     call rdhead (2, ntitle2, titlep, ltit)
cc    if (ntitle2 .le. 0)  then
cc       titlep(1) = ' '
cc    endif

c 84  continue
c     read(2,*,end=1010) ipath, nleg
c     skip label (x y z ipot rleg beta eta) and read the path
c     read(2,*)
      do 78 ileg=0,nleg
c        read(2,*,end=1010) (rat(j,ileg),j=1,3)
         nconv(ileg)=0
  78  continue

      do 88 n=1,3
         aa = rat(n,nleg)
         do 87 i=0,(nleg-2)
            j=nleg-i
            rat(n,j)=rat(n,j-1)
  87     continue
         rat(n,1)=aa
  88  continue
      do 89 i=1,nleg
         nq0(i)=0.
  89  continue

c nconv converts # of an atom in the nleg list of coordinates (rat) to
c its # in the full list of all atomic coordinates (rat1)
      do 94 i=1,natom
         do 91 n=1,3
            q1(n,i)=0.
   91    q0(n,i)=0.
            do 95 jl=1,nleg
            m=0
            do 93 n=1,3
               l=nint(100.*rat(n,jl))
               l1=nint(100.*rat1(n,i))
               if (abs(l-l1).le.1) m=m+1
   93       continue
            if (m.eq.3) then
              nconv(jl)=i
              go to 95
            endif
   95    continue
   94 continue

      atmu=0.
      inn=1
      nq1(inn)=1
      iq0=0
      nconv(0)=nconv(nleg)
      do 100 il=1,nleg
         l=nconv(il)
         do 101 jq=1,iq0
  101    if(nq0(jq).eq.l) go to 102
         iq0=iq0+1
         nq0(iq0)=l
  102    continue
         do 105 ii=1,natom
            a = dist(rat1(1,ii),rat1(1,l))
            if (a.le.rintr) then
               do 103 jn=1,inn
  103          if(nq1(jn).eq.ii) go to 105
               inn=inn+1
               nq1(inn)=ii
            endif
  105    continue
         nq1x=inn
         nq0x=iq0
         i=nconv(il)
         im=nconv(il-1)
         ip=nconv(il+1)
c        if (il.eq.1) im=nconv(nleg)
         if (il.eq.nleg) ip=nconv(1)
         atmass=atwtd(iz(i))
      do 100 n=1,3
         atmu=atmu + 0.25*( rnn(n,i,im)+rnn(n,i,ip) )**2 /atmass
  100 continue
      atmu=1./atmu
  108 continue

      do 115 i=1,natom
      do 115 n=1,3
  115 q0(n,i)=0.

c Build initial state vector |Q_j(0)> for the current path
      do 116 n=1,3
      do 116 il=1,nleg
         i=nconv(il)
         im=nconv(il-1)
         ip=nconv(il+1)
         if (il.eq.1) im=nconv(nleg)
         if (il.eq.nleg) ip=nconv(1)
         atmass=atwtd(iz(i))
         q0(n,i)=q0(n,i)+sqrt(atmu/atmass)*(rnn(n,im,i)-rnn(n,i,ip))/2.
  116 continue

c make sure it's normalized <Q_j(0)|Q_j(0)>=1
      q0q0=0.
      do 120 iq0=1,nq0x
      i=nq0(iq0)
      do 120 n=1,3
         q0q0=q0q0+q0(n,i)*q0(n,i)
 120  continue
      p00=nint(q0q0*1000.)/1000.
      if (abs(p00-1.d0).gt.5.d-4) then
         atmu=atmu/q0q0
         go to 108
      endif

c     to get THz units:
      wnorm=100.*w0/sqrt(amu0*10.) 
c*** moments
      a0=0.
      do 132 il=1,nq0x
      do 132 im=1,nq0x
         l=nq0(il)
         m=nq0(im)
      do 132 n1=1,3
      do 132 n2=1,3
         a0 = a0 + q0(n1,l)*dm(n1,n2,l,m)*q0(n2,m)/w0/w0
  132 continue
      we=wnorm*sqrt(a0)
      if (we.lt.1) then
c        recursion method is inapplicable, use statistics to set sig2
         sig2 = sig2x ( iphat(nconv(1)),  iphat(nconv(1)) )
         if (sig2.lt.1.d-6) sig2 = sig2mx
         return
      endif

      do 137 iset=1,nq1x
         i=nq1(iset)
      do 137 n1=1,3
            q1i=0.
            do 138 im=1,nq0x
               m=nq0(im)
            do 138 n2=1,3
               q1i=q1i+dm(n1,n2,i,m)*q0(n2,m)/w0/w0
  138       continue
         q1(n1,i) = q1i - a0*q0(n1,i)
  137 continue

      b0=0.
      do 139 i=1,natom
      do 139 n1=1,3
         b0=b0+q1(n1,i)*q1(n1,i)
  139 continue

      a1=0.
      do 150 iset=1,nq1x
         i=nq1(iset)
         do 150 n1=1,3
         q2=0.
         do 151 jset=1,nq1x
         j=nq1(jset)
            do 151 n2=1,3
               q2 = q2 + dm(n1,n2,i,j)*q1(n2,j)/w0/w0
  151       continue
            a1 = a1 + q1(n1,i)*q2
  150 continue

      a0=a0*wnorm**2
      a1=a1/b0
      a1=a1*wnorm**2
      b0=b0*wnorm**4

c** recursion sigma^2
      dd = (a0+a1)**2 - 4.*(a0*a1-b0)
      x1 = (a0+a1+sqrt(dd))/2.
      x2 = (a0+a1-sqrt(dd))/2.
      aa2 = (a1-x2)/(x1-x2)
c     aa2 = (a1-x2)/(x1-x2)*9./8.
      aa1 = (x1-a1)/(x1-x2)
      w1 = sqrt(x1)
      w2 = sqrt(x2)
      s1 = 3.1746/(atmu*w1*tanh(w1*7.6383/2./tk))
      s2 = 3.1746/(atmu*w2*tanh(w2*7.6383/2./tk))
      sigma2 = aa1*s1+aa2*s2
      sig2e = 3.1746/(atmu*we*tanh(we*7.6383/2./tk))

      if (ipath.ne.0) then
         write(slog,250) ipath,nleg,sigma2,atmu,w1,w2,aa1,aa2
         call wlog(slog)
         if (ir1_open) then
           write(ir1,250) ipath,nleg,sigma2,atmu,w1,w2,aa1,aa2
  250      format(1x,i3,4x,i1,3x,f9.5,2x,f7.3,2x,f7.2,2x,f7.2,
     1       4x,f5.3,2x,f5.3)
         endif
         if (iprdos.ne.0.and.ir2_open) then
            write(ir2,260) ipath,nleg,sig2e,atmu,we
  260       format(1x,i3,4x,i3,3x,f7.5,2x,f7.3,2x,f7.2)
         endif
      endif
      sig2 = sigma2

c     update maximum DW factors
      if (sig2.gt.sig2mx) sig2mx=sig2
      if (sig2.gt.sig2x( iphat(nconv(1)),  iphat(nconv(2)) )) then
         sig2x( iphat(nconv(1)),  iphat(nconv(2)) ) = sig2
         sig2x( iphat(nconv(2)),  iphat(nconv(1)) ) = sig2
      endif

      return
      end
c----------------------------------------------------
      subroutine sigte3 (iz1,iz2, sig2, alphat, thetad, reff, sig1,sig3)
c     single scattering only.

c     input: sig2
c     iz1, iz2 are iz at central atom and neighbor
c     alphat coeef of thermal expansion at high T
c     reff

c     output: sig1 sig3
      implicit double precision (a-h, o-z)
      real reff

c     con=hbar**2/kB*amu)*10**20   in ang**2 units
c     hbar = 1.054 572 666 e-34, amu = 1.660 540 e-27, 
c     kB = 1.380 6581 d-23
      parameter (con = 48.508 459 393 094)
      parameter (hbar = 1.054 572 666 e-34)
      parameter (amu = 1.660 540 e-27)
      parameter (xkb = 1.380 6581 e-23)

      ami=atwtd(iz1)*amu
      amj=atwtd(iz2)*amu

c     reduced mass
      xmu = 1 / (1/ami + 1/amj)
c     Einstein frequency
      omega = (2 * xkb * thetad) / (3 * hbar)
      xks = xmu * omega**2
      xk3 = xks**2 * reff * alphat / (3 * xkb)
      sig02 = hbar * omega / xks
      sig1 = -3 * (xk3 / xks) * sig2
      sig3 = 2 - (4.0/3.0) * (sig02 / sig2)**2
      sig3 = sig3 * (sig1 * sig2)

      return
      end
      subroutine sigm3(sig1, sig2, sig3, tk, alphat, thetae)
c     using correlated Einstein-model with a morse potential
c     Nguyen Van Hung & J.J.Rehr Phys. Rev. B 56 , 43

      implicit double precision (a-h, o-z)
      real sig02,  sig01, z
c     dimension alphat=[1/anstroems]
      parameter (bohr = 0.529 177 249d0)
      parameter(three= 3)
      parameter(four= 4 )
      parameter(fourthird= four/three)
      parameter(threequater= three/four)
              
      alphat= alphat * bohr  
      z=exp(- thetae/tk)
      sig02= (1-z)/ (1+z) * sig2
      sig01 = alphat * sig02 * threequater
      sig1 = sig01 * sig2 / sig02
      sig3 = (2- fourthird * (sig02/sig2) **2)* sig1 * sig2

      return
      end
c---------------------------------------------------------------------
c     program sigcl.f
c
c     calculates debye-waller factors for each multiple
c     scattering path using Debye-Model correlations
c
c     files:  input  pathd_all.dat  multiple scattering path data
c             output fort.3  sig**2 vs path
c                    fort.2  long output
c
c     version 1  (29 july 91)
c
c     coded by j. rehr
c     path data from s. zabinsky
c
c     modified to use pdata.inp, Dec 1991, siz
c     Subroutine version, Dec 1991, siz
c     Hacked for classical Debye Model, 7/2006, Kevin Jorissen  !KJ
c
c---------------------------------------------------------------------

      subroutine sigcl (tk, thetad, rs, nlegx, nleg, rat, iz, sig2)
c               tk temperature in degrees K
c               thetad debye temp in degrees K
c               rs=wigner seitz or norman radius in bohr, averaged
c                  over entire problem
c                  (4pi/3)*rs**3 = sum( (4pi/3)rnrm**3 ) / N
c                  (sum is over all atoms in the problem)
c               nlegx used in dimensions of rat and iz
c               nleg nlegs in path
c               rat positions of each atom in path
c               iz atomic number of each atom in path
c               NB Units of distance in this routine
c                  are angstroms, including sig**2.  rs is in bohr.
c               sig2 is output debye waller factor

      implicit double precision (a-h,o-z)

c={../HEADERS/const.h
      parameter (pi = 3.14159 26535 89793 23846 26433d0)
      parameter (one = 1, zero = 0)
      parameter (third = one/3)
      parameter (raddeg = 180 / pi)
      complex*16 coni
      parameter (coni = (0,1))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919 158 292 677 512 811d0)

      parameter (bohr = 0.529 177 249d0, ryd  = 13.605 698d0)
      parameter (hart = 2 * ryd)
      parameter (alpinv = 137.035 989 56d0)
c     fine structure alpha
      parameter (alphfs = 1 / alpinv)
c= ../HEADERS/const.h}

c     nlegx is max number of atoms in any one path
      dimension rat(3,0:nlegx)
      dimension iz(0:nlegx)
c#mn
       external dist

c      parameters
c               x = k_d*R   (distance parameter)
c               R distance in angstroms
c               y = hbar omegad/kT = thetad/t
c               thetad debye temp in degrees K
c               tk temperature in degrees K
c               k_d = (6*pi**2 N/V) = debye wave number
c               N/V=1/(4pi/3rs**3)
c               rs=wigner seitz or norman radius in bohr
c               ami, amj masses at sites i and j in amu
c               I = int_0^1 (y/x) dw sin(wx)coth(wy/2)

c     Note:  There are nleg atoms including the central atom
c            index 0 and index nleg both refer to central atom,
c            which makes special code unnecessary later.

      sigtot=0
      do 800 il=1,nleg
      do 800 jl=il,nleg

c        calculate r_i-r_i-1 and r_j-r_j-1

         rij = dist (rat(1,il), rat(1,jl))
         call corrfn2 (rij, cij, thetad, tk, iz(il), iz(jl), rs)
         sig2ij=cij

         rimjm = dist (rat(1,il-1), rat(1,jl-1))
         call corrfn2 (rimjm, cimjm, thetad, tk, iz(il-1), iz(jl-1), rs)
         sig2ij=sig2ij+cimjm

         rijm = dist (rat(1,il), rat(1,jl-1))
         call corrfn2 (rijm, cijm, thetad, tk, iz(il), iz(jl-1), rs)
         sig2ij=sig2ij-cijm

         rimj = dist (rat(1,il-1), rat(1,jl))
         call corrfn2 (rimj, cimj, thetad, tk, iz(il-1), iz(jl), rs)
         sig2ij=sig2ij-cimj

         riim = dist (rat(1,il), rat(1,il-1))
         rjjm = dist (rat(1,jl), rat(1,jl-1))

         ridotj=(rat(1,il)-rat(1,il-1))*(rat(1,jl)-rat(1,jl-1))+
     1          (rat(2,il)-rat(2,il-1))*(rat(2,jl)-rat(2,jl-1))+
     2          (rat(3,il)-rat(3,il-1))*(rat(3,jl)-rat(3,jl-1))
         ridotj=ridotj/(riim*rjjm)

c        double count i .ne. j  terms
         if(jl.ne.il) sig2ij=2*sig2ij
         sig2ij=sig2ij*ridotj
         sigtot=sigtot+sig2ij

  800 continue
      sig2=sigtot/4

c     sig2 is in bohr**2, just as we wanted for ff2chi
      return
      end



      subroutine corrfn2(rij,cij,thetad,tk,iz1,iz2,rsavg)
c     subroutine calculates correlation function
c     c(ri,rj)=<xi xj> in the Debye approximation
c
c             =(1/N)sum_k exp(ik.(Ri-Rj))(1/sqrt(mi*mj))*
c              (hbar/2w_k)*coth(beta hbar w_k/2)
c             = (3kT/mu w_d**2)*sqrt(mu**2/mi*mj)*I
c
c      parameters
c               x = k_d*R   (distance parameter)
c               R distance in angstroms
c               y = hbar omegad/kT = thetad/t
c               thetad debye temp in degrees K
c               tk temperature in degrees K
c               k_d = (6*pi**2 N/V) = debye wave number
c               N/V=1/(4pi/3rs**3)
c               rs=wigner seitz or norman radius in bohr
c               ami, amj masses at sites i and j in amu
c               I = int_0^1 (y/x) dw sin(wx)coth(wy/2)
c
c      solution by numerical integration
c
      implicit double precision (a-h, o-z)
      common /xy/ x, yinv

c={../HEADERS/const.h
      parameter (pi = 3.14159 26535 89793 23846 26433d0)
      parameter (one = 1, zero = 0)
      parameter (third = one/3)
      parameter (raddeg = 180 / pi)
      complex*16 coni
      parameter (coni = (0,1))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919 158 292 677 512 811d0)

      parameter (bohr = 0.529 177 249d0, ryd  = 13.605 698d0)
      parameter (hart = 2 * ryd)
      parameter (alpinv = 137.035 989 56d0)
c     fine structure alpha
      parameter (alphfs = 1 / alpinv)
c= ../HEADERS/const.h}

c     con=hbar**2/kB*amu)*10**20   in ang**2 units
c     hbar = 1.054 572 666 e-34, amu = 1.660 540 e-27, 
c     kB = 1.380 6581 d-23
      parameter (con = 48.508 459 393 094)
c#mn
       external atwtd

c     external fn
c     rij=2.55
c     tk=295
c     thetad=315
c     ami=amj=63.55 at wt for Cu
c     rs=2.7

      ami=atwtd(iz1)
      amj=atwtd(iz2)
      rs=rsavg
c     thetad in degrees K, t temperature in degrees K
c     y=thetad/tk
      yinv=tk/thetad
      xkd=(9*pi/2)**(third)/(rs*bohr)
      fac=(3/2.)*con/(thetad*sqrt(ami*amj))
      rj=rij
      x=xkd*rj
c     call numerical integration
      call bingrt2 (grater, eps, nx)
      cij=fac*grater
      return
      end
      double precision function fn2(w)
      implicit double precision (a-h,o-z)
      common/xy/x,yinv
c     fn=(sin(wx)/x)*coth(wy/2)
c     change code to allow t=0 without bombing
c     fn=2/y
      fn2=2*yinv
      if(w.lt.1.e-20) return
      fac=w
      if(x.gt.0.) fac=sin(w*x)/x
!      emwy=0.
!      if(yinv.gt.0.0125) emwy=exp(-w/yinv)
!      emwy=exp(-w/yinv)
!      fn2=fac*(1+emwy)/(1-emwy)
      fn2=fac*2*yinv/w  !KJ coup de grace
      return
      end
c-----------------------------------------------
      subroutine bingrt2 (b, eps, n)
c     subroutine calculates integrals between [0,1]
c      b = int_0^1 f(z) dz
c     by trapezoidal rule and binary refinement
c     (romberg integration)
c     coded by j rehr (10 Feb 92)
c     see, e.g., numerical recipes for discussion
c     and a much fancier version
c-----------------------------------------------
c     del=dz  itn=2**n tol=1.e-5
c     starting values
      implicit double precision (a-h,o-z)
      common /xy/x,yinv
      character*512 slog
c     external fn
c     error is approximately 2**(-2n) ~ 10**(-.6n)
c     so nmax=10 implies an error of 1.e-6
      parameter(nmax = 10, tol = 1.e-5)
      parameter(zero=0, one=1)
      n=0
      itn=1
      del=1.
      bn=(fn2(zero)+fn2(one))/2
      bo=bn
 10   continue
c     nth iteration
c     b_n+1=(b_n)/2+deln*sum_0^2**n f([2n-1]deln)
      n=n+1
      if(n.gt.nmax) go to 40
      del=del/2
      sum=0.
      do 20 i=1, itn
      zi=(2*i-1)*del
 20   sum=sum+fn2(zi)
c     bnp1=b_n+1 is current value of integral
      bnp1=bn/2+del*sum
c     cancel leading error terms b=[4b-bn]/3
c     note: this is the first term in the
c     neville table - remaining errors were
c     found too small to justify the added code
      b=(4*bnp1-bn)/3
      eps=abs((b-bo)/b)
      if(eps.lt.tol) goto 60
      bn=bnp1
      bo=b
      itn=itn*2
      goto 10
 40   write(slog,50) n,itn, b,eps
      call wlog(slog)
 50   format(' not converged, n,itn,b,eps=',
     1  2i4,2e14.6)
      return
 60   continue
      return
      end
c///////////////////////////////////////////////////////////////////////
c Distribution:  COMMON 1.0
c Copyright (c) [2002] University of Washington
c 
c This software was prepared in part with US Government Funding under
c DOE contract DE-FG03-97ER45623.

c Redistribution and use of this Distribution in source and binary
c formats, with or without modification is permitted, provided the 
c following conditions are met:
c 
c Redistributions must retain the above notices and the following list
c of conditions and disclaimer;
c 
c Modified formats carry the marking
c     "Based on or developed using Distribution: COMMON 1.0
c      COMMON 1.0 Copyright (c) [2002] University of Washington"
c 
c Recipient acknowledges the right of the University of Washington to
c prepare uses of this Distribution and its modifications that may be
c substantially similar or functionally equivalent to
c Recipient-prepared modifications.
c
c Recipient and anyone obtaining access to the Distribution through
c recipient's actions accept all risk associated with possession and
c use of the Distribution.
c
c THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESS OR IMPLIED
c WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
c MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
c IN NO EVENT SHALL THE UNIVERSITY OF WASHINGTON OR CONTRIBUTORS TO THE
c DISTRIBUTION BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
c EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
c PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
c REVENUE; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
c LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
c NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
c SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
c///////////////////////////////////////////////////////////////////////
c License is applicable for routines below, until otherwise specified.
c
      subroutine chopen (ios, fname, mod)
c     Writes error msg and stops if error in ios flag from open
c     statement.  fname is filename, mod is module with failed open.
      character*(*) fname, mod
      character*512 slog

c     open successful
      if (ios .le. 0)  return

c     error opening file, tell user and die.
      i = istrln(fname)
      j = istrln(mod)
      write(slog,100)  fname(1:i), mod(1:j)
      call wlog(slog)

  100 format (' Error opening file, ', a, 
     2        ' in module ', a)

      call wlog(' Fatal error')
      call par_stop('CHOPEN')
      end
      subroutine fixdsp (dxorg, dxnew, dgc0, dpc0, dgcx, dpcx, jnew)

c     This fixes up the dirac spinor components (dgc and dpc) from ATOM
c     for the xsect code.

      implicit double precision (a-h, o-z)

c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}
c={../HEADERS/const.h
      parameter (pi = 3.14159 26535 89793 23846 26433d0)
      parameter (one = 1, zero = 0)
      parameter (third = one/3)
      parameter (raddeg = 180 / pi)
      complex*16 coni
      parameter (coni = (0,1))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919 158 292 677 512 811d0)

      parameter (bohr = 0.529 177 249d0, ryd  = 13.605 698d0)
      parameter (hart = 2 * ryd)
      parameter (alpinv = 137.035 989 56d0)
c     fine structure alpha
      parameter (alphfs = 1 / alpinv)
c= ../HEADERS/const.h}

      dimension dgc0(251), dpc0(251)
      dimension dgcx(nrptx), dpcx(nrptx)

      dimension xorg(nrptx), xnew(nrptx)

      parameter (xx00 = 8.8)

c     statement functions to do indexing.  delta is 'dx' for current
c     grid.  jjj is index of grid point immediately before 'r'
      xxx(j) = -xx00 + (j-1)*delta
      rrr(j) = exp (-xx00 + (j-1)*delta)
      jjj(r) = (log(r) + xx00) / delta + 1

c     Use linear interpolation in x whether necessary or not.  If
c     new grid is same as old, it shouldn't make any difference.

c     relation between x, r, and j.  xx00 = 8.8 for all grids
c     in this version, change it if more flexibility is necessary.
c     xx = -xx00 + (j-1)*delta
c     rr = exp (xx)
c     jj = (log(r) + xx00) / delta + 1; this is j immediately BELOW r

c     The dgc and dpc arrays are zero beyond a certain point, usually
c     inside the muffin tin radius.  Find this distance.
      do 100  i = 251, 1, -1
         if ( abs(dgc0(i)) .ge. 1.0d-11 .or. 
     1        abs(dpc0(i)) .ge. 1.0d-11 )  then
            imax = i
            goto 16
         endif
  100 continue
      call wlog(' Should never see this line from sub fixdsp')
   16 continue
c     jmax is the first point where both dpc and dgc are zero in
c     the original grid
      jmax = imax + 1
      if (jmax.gt.251) jmax = 251

      delta = dxorg
      do 10  j = 1, jmax
         xorg(j) = xxx(j)
   10 continue
      rmax = rrr(jmax)

c     How far out do we go in the new grid?  To the last new grid
c     point before jmax.  Everything will be zero beyond jmax.
      delta = dxnew
      jnew = jjj(rmax)
      do 20  j = 1, jnew
         xnew(j) = xxx(j)
   20 continue

c     interpolate to new grid using x, only inside of rmax
      do 30  j = 1, jnew
         call terp (xorg, dgc0,  jmax, 3, xnew(j), dgcx(j))
         call terp (xorg, dpc0,  jmax, 3, xnew(j), dpcx(j))
   30 continue

c     and zero the arrays past rmax
      do 32  j = jnew+1, nrptx
         dgcx(j) = 0
         dpcx(j) = 0
   32 continue

      return
      end
      subroutine fixdsx (iph, dxorg, dxnew, dgc, dpc, dgcn, dpcn)

c     This fixes up the dirac spinor components (dgc and dpc) from ATOM
c     for the xsect and phase codes.

      implicit double precision (a-h, o-z)

c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}
c={../HEADERS/const.h
      parameter (pi = 3.14159 26535 89793 23846 26433d0)
      parameter (one = 1, zero = 0)
      parameter (third = one/3)
      parameter (raddeg = 180 / pi)
      complex*16 coni
      parameter (coni = (0,1))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919 158 292 677 512 811d0)

      parameter (bohr = 0.529 177 249d0, ryd  = 13.605 698d0)
      parameter (hart = 2 * ryd)
      parameter (alpinv = 137.035 989 56d0)
c     fine structure alpha
      parameter (alphfs = 1 / alpinv)
c= ../HEADERS/const.h}

      dimension dgc(251,30,0:nphx+1), dpc(251,30,0:nphx+1)
      dimension dgcn(nrptx,30), dpcn(nrptx,30)

      dimension xorg(nrptx), xnew(nrptx)

      parameter (xx00 = 8.8)

c     statement functions to do indexing.  delta is 'dx' for current
c     grid.  jjj is index of grid point immediately before 'r'
      xxx(j) = -xx00 + (j-1)*delta
      rrr(j) = exp (-xx00 + (j-1)*delta)
      jjj(r) = (log(r) + xx00) / delta + 1

c     Use linear interpolation in x whether necessary or not.  If
c     new grid is same as old, it shouldn't make any difference.

c     relation between x, r, and j.  xx00 = 8.8 for all grids
c     in this version, change it if more flexibility is necessary.
c     xx = -xx00 + (j-1)*delta
c     rr = exp (xx)
c     jj = (log(r) + xx00) / delta + 1; this is j immediately BELOW r

c     The dgc and dpc arrays are zero beyond a certain point, usually
c     inside the muffin tin radius.  Find this distance.

      delta = dxorg
      do 10  j = 1, 251
         xorg(j) = xxx(j)
   10 continue

      delta = dxnew
      do 20  j = 1, nrptx
         xnew(j) = xxx(j)
   20 continue

      do 200 iorb = 1, 30
         imax = 0
         do 100  i = 251, 1, -1
            if ( abs(dgc(i,iorb,iph)) .ge. 1.0d-11 .or. 
     1           abs(dpc(i,iorb,iph)) .ge. 1.0d-11 )  then
               imax = i
               goto 16
            endif
  100    continue
   16    continue
         if (imax .eq. 0) then
            jnew = 0
            goto 35
         endif
c        jmax is the first point where both dpc and dgc are zero in
c        the original grid
         jmax = imax + 1
         if (jmax .gt. 251) jmax = 251

         delta = dxorg
         rmax = rrr(jmax)

c        How far out do we go in the new grid?  To the last new grid
c        point before jmax.  Everything will be zero beyond jmax.
         delta = dxnew
         jnew = jjj(rmax)

c        interpolate to new grid using x, only inside of rmax
         do 30  j = 1, jnew
            call terp(xorg,dgc(1,iorb,iph),jmax,3, xnew(j),dgcn(j,iorb))
            call terp(xorg,dpc(1,iorb,iph),jmax,3, xnew(j),dpcn(j,iorb))
   30    continue

c        and zero the arrays past rmax
   35    do 40  j = jnew+1, nrptx
            dgcn(j,iorb) = 0
            dpcn(j,iorb) = 0
   40    continue
  200 continue

      return
      end
      subroutine fixvar (rmt, edens, vtot, dmag,
     1                   vint, rhoint, dxorg, dxnew, jumprm,
     2                   vjump, ri, vtotph, rhoph, dmagx)

      implicit double precision (a-h, o-z)

c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}
c={../HEADERS/const.h
      parameter (pi = 3.14159 26535 89793 23846 26433d0)
      parameter (one = 1, zero = 0)
      parameter (third = one/3)
      parameter (raddeg = 180 / pi)
      complex*16 coni
      parameter (coni = (0,1))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919 158 292 677 512 811d0)

      parameter (bohr = 0.529 177 249d0, ryd  = 13.605 698d0)
      parameter (hart = 2 * ryd)
      parameter (alpinv = 137.035 989 56d0)
c     fine structure alpha
      parameter (alphfs = 1 / alpinv)
c= ../HEADERS/const.h}


      dimension edens(251), vtot (251), dmag(251)
      dimension vtotph(nrptx), rhoph(nrptx), dmagx(nrptx)
      dimension ri(nrptx)
      dimension xorg(nrptx), xnew(nrptx)

      parameter (xx00 = 8.8)

c     statement functions to do indexing.  delta is 'dx' for current
c     grid.  jjj is index of grid point immediately before 'r'
      xxx(j) = -xx00 + (j-1)*delta
      rrr(j) = exp (-xx00 + (j-1)*delta)
      jjj(r) = (log(r) + xx00) / delta + 1

c     PHASE needs
c     vtot = total potential including gs xcorr, no r**2
c     edens = rho, charge density, no factor of 4*pi, no r**2
c     From overlapping, vtot = potential only, ok as is
c                       edens = density*4*pi, so fix this here.
c     ri = r grid through imt+1

c     Only values inside the muffin tin are used, except that XCPOT
c     (in PHASE) uses values at imt+1 and requires these to be the
c     interstitial values.  So set the last part of the arrays to
c     interstitial values...

c     Use linear interpolation in x whether necessary or not.  If
c     new grid is same as old, it shouldn't make any difference.

c     relation between x, r, and j.  xx00 = 8.8 for all grids
c     in this version, change it if more flexibility is necessary.
      
c     xx = -xx00 + (j-1)*delta
c     rr = exp (xx)
c     jj = (log(r) + xx00) / delta + 1; this is j immediately BELOW r

      delta = dxorg
      jmtorg = jjj(rmt)
      jriorg = jmtorg + 1
      jrior1 = jriorg + 1
      do 10  j = 1, jrior1
         xorg(j) = xxx(j)
   10 continue

      delta = dxnew
      jmtnew = jjj(rmt)
      jrinew = jmtnew + 1
      jrine1 = jrinew + 1
      do 20  j = 1, jrine1
         xnew(j) = xxx(j)
   20 continue

c     interpolate to new grid using x, only inside of muffintin
c     jri (first interstitial point) must be set to interstitial value
      do 30  j = 1, jrinew
         call terp (xorg, vtot,  jriorg, 3, xnew(j), vtotph(j))
         call terp (xorg, edens, jrior1, 3, xnew(j), rhoph(j))
         call terp (xorg, dmag,  jrior1, 3, xnew(j), dmagx(j))
   30 continue

      if (jumprm .eq. 1) then
         xmt = log(rmt)
         call terp (xorg, vtot,  jriorg, 3, xmt, vmt)
         vjump = vint - vmt
      endif
      if (jumprm .gt. 0) then
         do 90  j = 1, jrinew
            vtotph(j) = vtotph(j) + vjump
   90    continue
      endif

      delta = dxnew
      do 180  j = 1, nrptx
         ri(j) = rrr(j)
  180 continue
      do 190  j = 1, jrinew
         rhoph(j) = rhoph(j)/(4*pi)
  190 continue
      do 200  j = jrinew+1, nrptx
         vtotph(j) = vint
         rhoph(j) = rhoint/(4*pi)
c fix later : need to calculate interstitial dmint
c      want interpolation beyond mt also
         dmagx(j) = 0.0d0
  200 continue

      return
      end
      subroutine getorb (iz, ihole, xion, iunf, norb, norbco, iorb,
     1                  iholep, nqn, nk, xnel, xnval, xmag)
c     Gets orbital data for chosen element.  Input is:
c       iz - atomic number of desired element,
c       ihole - index of core-hole orbital
c       xion  - ionicity (usually zero)
c     other arguments are output.
c       norb - total number of orbitals
c       norbco - number of core orbitals
c       iorb - index of orbital for making projections (last occupied)
c       iholep - index of core hole orbital in compacted list
c       nqn - principal quantum number for each orbital
c       nk - quantum number kappa for each orbital
c       xnel - occupation for each orbital
c       xnval - valence occupation for each orbital
c       xmag - spin magnetization for each orbital
c     Feel free to change occupation numbers for element of interest.
c     ival(i) is necessary only for partly nonlocal exchange model.
c     iocc(i) and ival(i) can be fractional
c     But you have to keep the sum of iocc(i) equal to nuclear charge.
c     Also ival(i) should be equal to iocc(i) or zero. 
c     Otherwise you have to change this subroutine or contact authors 
c     for help.

      implicit double precision (a-h, o-z)

c     Written by Steven Zabinsky, July 1989
c     modified (20 aug 1989)  table increased to at no 99
c     Recipe for final state configuration is changed. Valence
c     electron occupations are added. ala 17.1.1996

c     Table for each element has occupation of the various levels.
c     The order of the levels in each array is:

c     element  level     principal qn (nqn), kappa qn (nk)
c           1  1s        1  -1
c           2  2s        2  -1
c           3  2p1/2     2   1
c           4  2p3/2     2  -2
c           5  3s        3  -1
c           6  3p1/2     3   1
c           7  3p3/2     3  -2
c           8  3d3/2     3   2
c           9  3d5/2     3  -3
c          10  4s        4  -1
c          11  4p1/2     4   1
c          12  4p3/2     4  -2
c          13  4d3/2     4   2
c          14  4d5/2     4  -3
c          15  4f5/2     4   3
c          16  4f7/2     4  -4
c          17  5s        5  -1
c          18  5p1/2     5   1
c          19  5p3/2     5  -2
c          20  5d3/2     5   2
c          21  5d5/2     5  -3
c          22  5f5/2     5   3
c          23  5f7/2     5  -4
c          24  6s        6  -1
c          25  6p1/2     6   1
c          26  6p3/2     6  -2
c          27  6d3/2     6   2
c          28  6d5/2     6  -3
c          29  7s        7  -1

      dimension nqn(30), nk(30), xnel(30), xnval(30), xmag(30)
      dimension kappa (29)
      real iocc, ival, ispn
      dimension iocc (100, 29), ival (100, 29), ispn (100, 29)
      dimension nnum (29), iorb(-4:3)
      character*512 slog

c     kappa quantum number for each orbital
c     k = - (j + 1/2)  if l = j - 1/2
c     k = + (j + 1/2)  if l = j + 1/2
      data kappa /-1,-1, 1,-2,-1,   1,-2, 2,-3,-1,   1,-2, 2,-3, 3,
     1            -4,-1, 1,-2, 2,  -3, 3,-4,-1, 1,  -2, 2,-3,-1/

c     principal quantum number (energy eigenvalue)
      data nnum  /1,2,2,2,3,  3,3,3,3,4,  4,4,4,4,4,
     1            4,5,5,5,5,  5,5,5,6,6,  6,6,6,7/

c     occupation of each level for z = 1, 99
      data (iocc( 1,i),i=1,29)  /1,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival( 1,i),i=1,29)  /1,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn( 1,i),i=1,29)  /1,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc( 2,i),i=1,29)  /2,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival( 2,i),i=1,29)  /2,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn( 2,i),i=1,29)  /1,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc( 3,i),i=1,29)  /2,1,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival( 3,i),i=1,29)  /0,1,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn( 3,i),i=1,29)  /0,1,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc( 4,i),i=1,29)  /2,2,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival( 4,i),i=1,29)  /0,2,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn( 4,i),i=1,29)  /0,1,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc( 5,i),i=1,29)  /2,2,1,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival( 5,i),i=1,29)  /0,2,1,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn( 5,i),i=1,29)  /0,0,1,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

c     data (iocc( 6,i),i=1,29)  /2,2,2,0,0,  0,0,0,0,0,  0,0,0,0,0,
      data (iocc( 6,i),i=1,29)  /2,1,2,1,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
c     data (ival( 6,i),i=1,29)  /0,2,2,0,0,  0,0,0,0,0,  0,0,0,0,0,
      data (ival( 6,i),i=1,29)  /0,1,2,1,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn( 6,i),i=1,29)  /0,0,1,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc( 7,i),i=1,29)  /2,2,2,1,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival( 7,i),i=1,29)  /0,2,2,1,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn( 7,i),i=1,29)  /0,0,0,1,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc( 8,i),i=1,29)  /2,2,2,2,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival( 8,i),i=1,29)  /0,2,2,2,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn( 8,i),i=1,29)  /0,0,0,1,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc( 9,i),i=1,29)  /2,2,2,3,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival( 9,i),i=1,29)  /0,2,2,3,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn( 9,i),i=1,29)  /0,0,0,1,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(10,i),i=1,29)  /2,2,2,4,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(10,i),i=1,29)  /0,0,2,4,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(10,i),i=1,29)  /0,0,0,1,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(11,i),i=1,29)  /2,2,2,4,1,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(11,i),i=1,29)  /0,0,2,4,1,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(11,i),i=1,29)  /0,0,0,0,1,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(12,i),i=1,29)  /2,2,2,4,1,  1,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(12,i),i=1,29)  /0,0,0,0,1,  1,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(12,i),i=1,29)  /0,0,0,0,1,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(13,i),i=1,29)  /2,2,2,4,2,  1,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(13,i),i=1,29)  /0,0,0,0,2,  1,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(13,i),i=1,29)  /0,0,0,0,0,  1,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(14,i),i=1,29)  /2,2,2,4,2,  2,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(14,i),i=1,29)  /0,0,0,0,2,  2,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(14,i),i=1,29)  /0,0,0,0,0,  1,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(15,i),i=1,29)  /2,2,2,4,2,  2,1,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(15,i),i=1,29)  /0,0,0,0,2,  2,1,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(15,i),i=1,29)  /0,0,0,0,0,  0,1,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(16,i),i=1,29)  /2,2,2,4,2,  2,2,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(16,i),i=1,29)  /0,0,0,0,2,  2,2,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(16,i),i=1,29)  /0,0,0,0,0,  0,1,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(17,i),i=1,29)  /2,2,2,4,2,  2,3,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(17,i),i=1,29)  /0,0,0,0,2,  2,3,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(17,i),i=1,29)  /0,0,0,0,0,  0,1,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(18,i),i=1,29)  /2,2,2,4,2,  2,4,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(18,i),i=1,29)  /0,0,0,0,2,  2,4,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(18,i),i=1,29)  /0,0,0,0,0,  0,1,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(19,i),i=1,29)  /2,2,2,4,2,  2,4,0,0,1,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(19,i),i=1,29)  /0,0,0,0,2,  2,4,0,0,1,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(19,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,1,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(20,i),i=1,29)  /2,2,2,4,2,  2,4,0,0,1,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(20,i),i=1,29)  /0,0,0,0,0,  2,4,0,0,1,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(20,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,1,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(21,i),i=1,29)  /2,2,2,4,2,  2,4,1,0,1,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(21,i),i=1,29)  /0,0,0,0,0,  2,4,1,0,1,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(21,i),i=1,29)  /0,0,0,0,0,  0,0,1,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(22,i),i=1,29)  /2,2,2,4,2,  2,4,2,0,1,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(22,i),i=1,29)  /0,0,0,0,0,  2,4,2,0,1,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(22,i),i=1,29)  /0,0,0,0,0,  0,0,2,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(23,i),i=1,29)  /2,2,2,4,2,  2,4,3,0,1,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(23,i),i=1,29)  /0,0,0,0,0,  2,4,3,0,1,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(23,i),i=1,29)  /0,0,0,0,0,  0,0,3,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(24,i),i=1,29)  /2,2,2,4,2,  2,4,4,0,1,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(24,i),i=1,29)  /0,0,0,0,0,  2,4,4,0,1,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(24,i),i=1,29)  /0,0,0,0,0,  0,0,4,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(25,i),i=1,29)  /2,2,2,4,2,  2,4,4,1,1,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(25,i),i=1,29)  /0,0,0,0,0,  0,0,4,1,1,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(25,i),i=1,29)  /0,0,0,0,0,  0,0,4,1,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(26,i),i=1,29)  /2,2,2,4,2,  2,4,4,2,1,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(26,i),i=1,29)  /0,0,0,0,0,  0,0,4,2,1,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(26,i),i=1,29)  /0,0,0,0,0,  0,0,2,2,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(27,i),i=1,29)  /2,2,2,4,2,  2,4,4,3,1,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(27,i),i=1,29)  /0,0,0,0,0,  0,0,4,3,1,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(27,i),i=1,29)  /0,0,0,0,0,  0,0,0,3,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(28,i),i=1,29)  /2,2,2,4,2,  2,4,4,4,1,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(28,i),i=1,29)  /0,0,0,0,0,  0,0,4,4,1,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(28,i),i=1,29)  /0,0,0,0,0,  0,0,0,1,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(29,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,1,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(29,i),i=1,29)  /0,0,0,0,0,  0,0,4,6,1,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(29,i),i=1,29)  /0,0,0,0,0,  0,0,0,1,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(30,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,1,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(30,i),i=1,29)  /0,0,0,0,0,  0,0,4,6,1,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(30,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,1,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(31,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(31,i),i=1,29)  /0,0,0,0,0,  0,0,4,6,2,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(31,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(32,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(32,i),i=1,29)  /0,0,0,0,0,  0,0,4,6,2,  2,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(32,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(33,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,1,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(33,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,2,  2,1,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(33,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,1,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(34,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,2,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(34,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,2,  2,2,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(34,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,1,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(35,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,3,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(35,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,2,  2,3,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(35,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,1,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(36,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(36,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,2,  2,4,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(36,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,1,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(37,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,0,0,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(37,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,2,  2,4,0,0,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(37,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(38,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,0,0,0,
     1                           0,2,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(38,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  2,4,0,0,0,
     1                           0,2,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(38,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(39,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,1,0,0,
     1                           0,2,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(39,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  2,4,1,0,0,
     1                           0,2,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(39,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,1,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(40,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,2,0,0,
     1                           0,2,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(40,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  2,4,2,0,0,
     1                           0,2,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(40,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,2,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(41,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,0,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(41,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  2,4,4,0,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(41,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,3,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(42,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,1,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(42,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  2,4,4,1,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(42,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,4,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(43,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,1,0,
     1                           0,2,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(43,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,4,1,0,
     1                           0,2,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(43,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,4,1,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(44,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,3,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(44,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,4,3,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(44,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,2,2,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(45,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,4,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(45,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,4,4,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(45,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,2,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(46,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(46,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,4,6,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(46,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,1,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(47,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(47,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,4,6,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(47,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(48,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,2,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(48,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,4,6,0,
     1                           0,2,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(48,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(49,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,2,1,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(49,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,4,6,0,
     1                           0,2,1,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(49,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,1,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(50,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,2,2,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(50,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,4,6,0,
     1                           0,2,2,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(50,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,1,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(51,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,2,2,1,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(51,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,4,6,0,
     1                           0,2,2,1,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(51,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,1,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(52,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,2,2,2,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(52,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,2,2,2,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(52,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,1,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(53,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,2,2,3,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(53,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,2,2,3,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(53,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,1,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(54,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,2,2,4,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(54,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,2,2,4,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(54,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,1,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(55,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,2,2,4,0,  0,0,0,1,0,  0,0,0,0/
      data (ival(55,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,2,2,4,0,  0,0,0,1,0,  0,0,0,0/
      data (ispn(55,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,1,0,  0,0,0,0/

      data (iocc(56,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,2,2,4,0,  0,0,0,2,0,  0,0,0,0/
      data (ival(56,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,2,2,4,0,  0,0,0,2,0,  0,0,0,0/
      data (ispn(56,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,1,0,  0,0,0,0/

      data (iocc(57,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,2,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ival(57,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,2,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ispn(57,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,1,  0,0,0,0,0,  0,0,0,0/

      data (iocc(58,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,1,
     1                           0,2,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ival(58,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,1,
     1                           0,0,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ispn(58,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,1,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(59,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,2,
     1                           0,2,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ival(59,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,2,
     1                           0,0,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ispn(59,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,2,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(60,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,3,
     1                           0,2,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ival(60,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,3,
     1                           0,0,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ispn(60,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,3,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(61,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,4,
     1                           0,2,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ival(61,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,4,
     1                           0,0,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ispn(61,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,4,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(62,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,5,
     1                           0,2,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ival(62,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,5,
     1                           0,0,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ispn(62,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,5,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(63,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           0,2,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ival(63,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,6,
     1                           0,0,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ispn(63,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,6,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(64,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           1,2,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ival(64,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,6,
     1                           1,0,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ispn(64,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,6,
     1                           1,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(65,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           2,2,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ival(65,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,6,
     1                           2,0,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ispn(65,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,4,
     1                           2,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(66,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           3,2,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ival(66,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,6,
     1                           3,0,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ispn(66,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,2,
     1                           3,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(67,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           4,2,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ival(67,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,6,
     1                           4,0,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ispn(67,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           4,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(68,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           5,2,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ival(68,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,6,
     1                           5,0,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ispn(68,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           3,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(69,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           6,2,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ival(69,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,6,
     1                           6,0,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ispn(69,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           2,0,0,0,0,  0,0,0,2,0,  0,0,0,0/

      data (iocc(70,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           7,2,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ival(70,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,6,
     1                           7,0,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ispn(70,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           1,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(71,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ival(71,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,6,
     1                           8,0,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ispn(71,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,1,  0,0,0,0,0,  0,0,0,0/

      data (iocc(72,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,2,  0,0,0,2,0,  0,0,0,0/
      data (ival(72,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,6,
     1                           8,0,2,4,2,  0,0,0,2,0,  0,0,0,0/
      data (ispn(72,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,2,  0,0,0,0,0,  0,0,0,0/

      data (iocc(73,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,3,  0,0,0,2,0,  0,0,0,0/
      data (ival(73,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,6,
     1                           8,0,0,0,3,  0,0,0,2,0,  0,0,0,0/
      data (ispn(73,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,3,  0,0,0,0,0,  0,0,0,0/

      data (iocc(74,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,3,  1,0,0,2,0,  0,0,0,0/
c    1                           8,2,2,4,4,  0,0,0,2,0,  0,0,0,0/
      data (ival(74,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,6,
     1                           8,0,0,0,3,  1,0,0,2,0,  0,0,0,0/
c    1                           8,0,0,0,4,  0,0,0,2,0,  0,0,0,0/
      data (ispn(74,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,4,  0,0,0,0,0,  0,0,0,0/

      data (iocc(75,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  1,0,0,2,0,  0,0,0,0/
      data (ival(75,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,4,  1,0,0,2,0,  0,0,0,0/
      data (ispn(75,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,4,  1,0,0,0,0,  0,0,0,0/

      data (iocc(76,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  2,0,0,2,0,  0,0,0,0/
      data (ival(76,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,4,  2,0,0,2,0,  0,0,0,0/
      data (ispn(76,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,2,  2,0,0,0,0,  0,0,0,0/

      data (iocc(77,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  3,0,0,2,0,  0,0,0,0/
      data (ival(77,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,4,  3,0,0,2,0,  0,0,0,0/
      data (ispn(77,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  3,0,0,0,0,  0,0,0,0/

      data (iocc(78,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  5,0,0,1,0,  0,0,0,0/
      data (ival(78,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,4,  5,0,0,1,0,  0,0,0,0/
      data (ispn(78,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  2,0,0,0,0,  0,0,0,0/

      data (iocc(79,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,1,0,  0,0,0,0/
      data (ival(79,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,4,  6,0,0,1,0,  0,0,0,0/
      data (ispn(79,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  1,0,0,0,0,  0,0,0,0/

      data (iocc(80,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,0,  0,0,0,0/
      data (ival(80,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,4,  6,0,0,2,0,  0,0,0,0/
      data (ispn(80,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,1,0,  0,0,0,0/

      data (iocc(81,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,1,  0,0,0,0/
      data (ival(81,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,4,  6,0,0,2,1,  0,0,0,0/
      data (ispn(81,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,1,  0,0,0,0/

      data (iocc(82,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,2,  0,0,0,0/
      data (ival(82,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,4,  6,0,0,2,2,  0,0,0,0/
      data (ispn(82,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,1,  0,0,0,0/

      data (iocc(83,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,2,  1,0,0,0/
      data (ival(83,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,4,  6,0,0,2,2,  1,0,0,0/
      data (ispn(83,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  1,0,0,0/

      data (iocc(84,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,2,  2,0,0,0/
      data (ival(84,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,4,  6,0,0,2,2,  2,0,0,0/
      data (ispn(84,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  1,0,0,0/

      data (iocc(85,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,2,  3,0,0,0/
      data (ival(85,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,2,2,  3,0,0,0/
      data (ispn(85,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  1,0,0,0/

      data (iocc(86,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,2,  4,0,0,0/
      data (ival(86,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,2,2,  4,0,0,0/
      data (ispn(86,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  1,0,0,0/

      data (iocc(87,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,2,  4,0,0,1/
      data (ival(87,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,2,2,  4,0,0,1/
      data (ispn(87,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,1/

      data (iocc(88,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,2,  4,0,0,2/
      data (ival(88,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,2,  4,0,0,2/
      data (ispn(88,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,1/

      data (iocc(89,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,2,  4,1,0,2/
      data (ival(89,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,2,  4,1,0,2/
      data (ispn(89,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,1,0,0/

      data (iocc(90,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,2,  4,2,0,2/
      data (ival(90,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,2,  4,2,0,2/
      data (ispn(90,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,2,0,0/

      data (iocc(91,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,2,0,2,2,  4,1,0,2/
      data (ival(91,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,2,0,0,2,  4,1,0,2/
      data (ispn(91,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,2,0,0,0,  0,0,0,0/

      data (iocc(92,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,3,0,2,2,  4,1,0,2/
      data (ival(92,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,3,0,0,2,  4,1,0,2/
      data (ispn(92,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,1.5,0,0,0,  0,0,0,0/

      data (iocc(93,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,4,0,2,2,  4,1,0,2/
      data (ival(93,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,4,0,0,2,  4,1,0,2/
      data (ispn(93,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,4,0,0,0,  0,0,0,0/

      data (iocc(94,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,6,0,2,2,  4,0,0,2/
      data (ival(94,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,6,0,0,2,  4,0,0,2/
      data (ispn(94,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,5,0,0,0,  0,0,0,0/

      data (iocc(95,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,6,1,2,2,  4,0,0,2/
      data (ival(95,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,6,1,0,2,  4,0,0,2/
      data (ispn(95,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,5,1,0,0,  0,0,0,0/

      data (iocc(96,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,6,2,2,2,  4,0,0,2/
      data (ival(96,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,6,2,0,2,  4,0,0,2/
      data (ispn(96,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,5,2,0,0,  0,0,0,0/

      data (iocc(97,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,6,3,2,2,  4,0,0,2/
      data (ival(97,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,6,3,0,2,  4,0,0,2/
      data (ispn(97,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,3,3,0,0,  0,0,0,0/

      data (iocc(98,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,6,4,2,2,  4,0,0,2/
      data (ival(98,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,6,4,0,2,  4,0,0,2/
      data (ispn(98,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,1,4,0,0,  0,0,0,0/

      data (iocc(99,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,6,5,2,2,  4,0,0,2/
      data (ival(99,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,6,5,0,2,  4,0,0,2/
      data (ispn(99,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,4,0,0,  0,0,0,0/

      data (iocc(100,i),i=1,29) /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,6,6,2,2,  4,0,0,2/
      data (ival(100,i),i=1,29) /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,6,6,0,2,  4,0,0,2/
      data (ispn(100,i),i=1,29) /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,3,0,0,  0,0,0,0/

      if (iz .lt. 1  .or.  iz .gt. 99)  then
    8    format(' Atomic number ', i5, ' not available.')
         write(slog,8)  iz
         call wlog(slog)
         call par_stop('GETORB-0')
      endif

      ion = nint(xion)
      delion=xion-ion
      index = iz - ion
      ilast = 0
      iscr = 0
      iion = 0
      iholep = ihole

c     find last occupied orbital (ilast) and iion for delion.ge.0
      do 30 i=29,1,-1
         if (iion.eq.0 .and. iocc(index,i).gt.delion) iion=i
         if (ilast.eq.0 .and. iocc(index,i).gt.0) ilast=i
 30   continue

      if (ihole.gt.0) then
         if ( iocc(index,ihole) .lt. 1 ) then
           call wlog(' Cannot remove an electron from this level')
           call par_stop('GETORB-1')
         endif
      endif
      if (ihole.eq.ilast) then 
         if ( iocc(index,ihole)-delion.lt.1) then
           call wlog(' Cannot remove an electron from this level')
           call par_stop('GETORB-1')
        endif
      endif

c        the recipe for final state atomic configuration is changed
c        from iz+1 prescription, since sometimes it changed occupation
c        numbers in more than two orbitals. This could be consistent
c        only with s02=0.0. New recipe remedy this deficiency.

c     find where to put screening electron
      index1 = index + 1
      do 10  i = 1, 29
 10   if (iscr.eq.0 .and. (iocc(index1,i)-iocc(index,i)).gt.0.5) iscr=i
c     special case of hydrogen like ion
c     if (index.eq.1) iscr=2

c     find where to add or subtract charge delion (iion).
c     if (delion .ge. 0) then
c        removal of electron charge
c        iion is already found
      if (delion .lt. 0) then
c        addition of electron charge
         iion = iscr
c        except special cases
         if (ihole.ne.0 .and.
     1       iocc(index,iscr)+1-delion.gt.2*abs(kappa(iscr))) then
             iion = ilast
             if (ilast.eq.iscr .or. iocc(index,ilast)-delion.gt.
     1                          2*abs(kappa(ilast)) ) iion = ilast + 1
         endif
      endif

      norb = 0
      do 19 i=-4, 3
 19   iorb(i) = 0
      do 20  i = 1, 29
         if (iocc(index,i).gt.0 .or. (i.eq.iscr .and. ihole.gt.0)
     1       .or. (i.eq.iion .and. iocc(index,i)-delion.gt.0))  then
            if (i.ne.ihole .or. iocc(index,i).ge.1) then
               norb = norb + 1
               nqn(norb) = nnum(i)
               nk(norb)  = kappa(i)
               xnel(norb) = iocc(index,i)
               if (i.eq.ihole) then
                  xnel(norb) = xnel(norb) - 1
                  iholep = norb
               endif
               if (i.eq.iscr .and. ihole.gt.0)  xnel(norb)=xnel(norb)+1
               xnval(norb)= ival(index,i)
               if ((kappa(i).eq.-4 .or. kappa(i).eq.3) .and. iunf.eq.0)
     1           xnval(norb) = 0
               xmag(norb) = ispn(index,i)
               iorb(nk(norb)) = norb
               if (i.eq.ihole .and. xnval(norb).ge.1)
     1                         xnval(norb) = xnval(norb) - 1
               if (i.eq.iscr .and. ihole.gt.0) 
     1                         xnval(norb) = xnval(norb) + 1
               if (i.eq.iion)  xnel(norb) = xnel(norb) - delion
               if (i.eq.iion)  xnval(norb) = xnval(norb) - delion
            endif
         endif
   20 continue
      norbco = norb

c     check that all occupation numbers are within limits
      do 50 i = 1, norb
         if ( xnel(i).lt.0 .or.  xnel(i).gt.2*abs(nk(i)) .or.
     1       xnval(i).lt.0 .or. xnval(i).gt.2*abs(nk(i)) ) then
            write (slog,55) i
   55       format(' error in getorb.f. Check occupation number for ',
     1      i3, '-th orbital. May be a problem with ionicity.')
            call wlog(slog)
            call par_stop('GETORB-99')
         endif
  50  continue
c      do 60 i=1,norb
c60    xnval(i) = 0.0d0
c60    xnval(i) = xnel(i)
            
      return
      end
      double precision function getxk (e)
      implicit double precision (a-h, o-z)

c     Make xk from energy(in Hartrees) as
c          k =  sqrt(2*e)  for e > 0  (above the edge)
c          k = -sqrt(-2*e)  for e < 0  (below the edge)
      getxk = sqrt(abs(2*e))
      if (e .lt. 0)  getxk = - getxk
      return
      end
      subroutine sthead (ntitle, title, nph, iz, rmt, rnrm,
     1                  xion, ihole, ixc,
     2                  vr0, vi0, gamach, xmu, xf, vint, rs,
     2                  nohole, lreal,  rgrd)

c     SeT HEAD
c     This routine makes the file header, returned in head array.
c     header lines do not include a leading blank.
c     Last header line is not --------- end-of-header line

c     title lines coming into sthead include carriage control, since
c     they were read from potph.bin

      implicit double precision (a-h, o-z)

c={../HEADERS/const.h
      parameter (pi = 3.14159 26535 89793 23846 26433d0)
      parameter (one = 1, zero = 0)
      parameter (third = one/3)
      parameter (raddeg = 180 / pi)
      complex*16 coni
      parameter (coni = (0,1))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919 158 292 677 512 811d0)

      parameter (bohr = 0.529 177 249d0, ryd  = 13.605 698d0)
      parameter (hart = 2 * ryd)
      parameter (alpinv = 137.035 989 56d0)
c     fine structure alpha
      parameter (alphfs = 1 / alpinv)
c= ../HEADERS/const.h}
c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}
c={../HEADERS/vers.h
      character*12 vfeff
c                       123456789012  
      parameter (vfeff='Feff 8.50   ')
c= ../HEADERS/vers.h}

      dimension xion(0:nphx)
      dimension iz(0:nphx)
      dimension rmt(0:nphx)
      dimension rnrm(0:nphx)

      character*80 title(nheadx), store
      character*16 s1, s2

      character*10 shole(0:29)
      character*8  sout(0:7)
      data shole /'no hole',   'K  shell',  'L1 shell',  'L2 shell',
     2            'L3 shell',  'M1 shell',  'M2 shell',  'M3 shell',
     3            'M4 shell',  'M5 shell',  'N1 shell',  'N2 shell',
     4            'N3 shell',  'N4 shell',  'N5 shell',  'N6 shell',
     5            'N7 shell',  'O1 shell',  'O2 shell',  'O3 shell',
     6            'O4 shell',  'O5 shell',  'O6 shell',  'O7 shell',
     7            'P1 shell',  'P2 shell',  'P3 shell',  'P4 shell',
     8            'P5 shell',  'R1 shell'/
      data sout /'H-L exch', 'D-H exch', 'Gd state', 'DH - HL ',
     1           'DH + HL ', 'val=s+d ', 'sigmd(r)', 'sigmd=c '/


c     Fills head arrray, n = number of lines used.
c     Does not include line of dashes at the end.

      if (ntitle .ge. 1 ) then
         ii = istrln(title(1)) 
         if (ii.gt.1)  then
            write(store,100)  title(1)(1:), vfeff
         else
            write(store,102)  vfeff
         endif
      else
         write(store,102)   vfeff
      endif
  100 format( a55, t66, a12)
  102 format( t66, a12)
      title(1) = store
      nstor = 1

c     remove empty title lines
      do 120  ititle = 2, ntitle
         ii = istrln ( title (ititle) ) 
         if (ii.le.1)  goto 120
         nstor = nstor+1
         title(nstor) = title (ititle)
  120 continue
      ntitle = nstor

c     add more title lines
      if (xion(0) .ne. 0)  then
         ntitle = ntitle + 1
         write(title(ntitle),130)  iz(0), rmt(0)*bohr,
     1                    rnrm(0)*bohr, xion(0), shole(ihole)
      else
         ntitle = ntitle + 1
         write(title(ntitle),140)  iz(0), rmt(0)*bohr,
     1                    rnrm(0)*bohr, shole(ihole)
      endif
  130 format('Abs   Z=',i2, ' Rmt=',f6.3, ' Rnm=',f6.3,
     1       ' Ion=',f5.2,  1x,a10)
  140 format('Abs   Z=',i2,' Rmt=',f6.3,' Rnm=',f6.3, 1x,a10)
c     if (nohole.ge.0)  then
c        ntitle = ntitle + 1
c        write(title(ntitle),142)
c 142    format ('Calculations done with no core hole.')
c     endif
      if (lreal.ge.1 .or. (abs(rgrd - 0.05) .gt. 1.0e-5)) then
        ntitle = ntitle + 1
        s1 = ' '
        if (lreal.gt.1)  then
c        write(title(ntitle),144)
c 144    format ('Calculations done using only real phase shifts.')
         s1 = 'RPHASES'
        elseif (lreal.eq.1) then
c        ntitle = ntitle + 1
c        write(title(ntitle),145)
c 145    format ('Calculations done using only real self energy.')
         s1 = 'RSIGMA'
        endif
        s2 = '  '
        if (abs(rgrd - 0.05) .gt. 1.0e-5)  then
         write(s2,146)  rgrd
  146    format ('  RGRID', f7.4)
        endif
        ilen = istrln(s1)
        title(ntitle) = s1(1:ilen) // s2
      endif

      do 150  iph = 1, nph
         if (xion(iph) .ne. 0)  then
            ntitle = ntitle + 1
            write(title(ntitle),160)  iph, iz(iph),  rmt(iph)*bohr,
     1           rnrm(iph)*bohr, xion(iph)
         else
            ntitle = ntitle + 1
            write(title(ntitle),170)  iph, iz(iph),  rmt(iph)*bohr,
     1           rnrm(iph)*bohr
         endif
  150 continue
  160 format('Pot',i2,' Z=',i2,' Rmt=',f6.3,' Rnm=',f6.3,' Ion=',f5.2)
  170 format('Pot',i2,' Z=',i2,' Rmt=',f6.3,' Rnm=',f6.3)
      if (abs(vi0) .gt. 1.0e-8 .or. abs(vr0) .gt. 1.0e-8)  then
         ntitle = ntitle + 1
         write(title(ntitle),180)  gamach*hart, sout(ixc), vi0*hart,
     1                           vr0*hart
      else
         ntitle = ntitle + 1
         write(title(ntitle),190)  gamach*hart, sout(ixc)
      endif
      ntitle = ntitle + 1
  180 format('Gam_ch=',1pe9.3, 1x,a8, ' Vi=',1pe10.3, ' Vr=',1pe10.3)
  190 format('Gam_ch=',1pe9.3, 1x,a8)
  200 format('Mu=',1pe10.3, ' kf=',1pe9.3, ' Vint=',1pe10.3,
     x        ' Rs_int=',0pf6.3)
      write(title(ntitle),200)  xmu*hart, xf/bohr, vint*hart, rs

      return
      end

      subroutine wthead (io, ntitle, title)
c     Dump title lines to unit io, which must be open. 
      integer io, i, ll
      character*80 title(ntitle)

c     nice for UNIX to use with gnuplot etc.,
      do 310 i = 1, ntitle
         ll = istrln(title(i))
         write(io,300)  title(i)(1:ll)
  300    format (a)
  310 continue

      return
      end
      function itoken (word,flname)
c     chars in word assumed upper case, left justified
c     returns 0 if not a token, otherwise returns token

      character*(*) word
      character*4   w
      character*20 flname
      integer itoken

      w = word(1:4)
      call upper(w)
      
c     Tokens for feff.inp
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (flname(1:8).eq.'feff.inp') then
         if     (w .eq. 'ATOM')  then
            itoken = 1
         elseif (w .eq. 'HOLE')  then
            itoken = 2
         elseif (w .eq. 'OVER')  then
            itoken = 3
         elseif (w .eq. 'CONT')  then
            itoken = 4
         elseif (w .eq. 'EXCH')  then
            itoken = 5
         elseif (w .eq. 'ION ')  then
            itoken = 6
         elseif (w .eq. 'TITL')  then
            itoken = 7
         elseif (w .eq. 'FOLP')  then
            itoken = 8
         elseif (w .eq. 'RPAT' .or. w .eq. 'RMAX')  then
            itoken = 9
         elseif (w .eq. 'DEBY')  then
            itoken = 10
         elseif (w .eq. 'RMUL')  then
            itoken = 11
         elseif (w .eq. 'SS  ')  then
            itoken = 12
         elseif (w .eq. 'PRIN')  then
            itoken = 13
         elseif (w .eq. 'POTE')  then
            itoken = 14
         elseif (w .eq. 'NLEG')  then
            itoken = 15
         elseif (w .eq. 'CRIT')  then
            itoken = 16
         elseif (w .eq. 'NOGE')  then
            itoken = 17
         elseif (w .eq. 'IORD')  then
            itoken = 18
         elseif (w .eq. 'PCRI')  then
            itoken = 19
         elseif (w .eq. 'SIG2')  then
            itoken = 20
         elseif (w .eq. 'XANE')  then
            itoken = 21
         elseif (w .eq. 'CORR')  then
            itoken = 22
         elseif (w .eq. 'AFOL')  then
            itoken = 23
         elseif (w .eq. 'EXAF')  then
            itoken = 24
         elseif (w .eq. 'POLA')  then
            itoken = 25
         elseif (w .eq. 'ELLI')  then
            itoken = 26
         elseif (w .eq. 'RGRI')  then
            itoken = 27
         elseif (w .eq. 'RPHA')  then
            itoken = 28
         elseif (w .eq. 'NSTA')  then
            itoken = 29
         elseif (w .eq. 'NOHO')  then
            itoken = 30
         elseif (w .eq. 'SIG3')  then
            itoken = 31
         elseif (w .eq. 'JUMP')  then
            itoken = 32
         elseif (w .eq. 'MBCO')  then
            itoken = 33
         elseif (w .eq. 'SPIN')  then
            itoken = 34
         elseif (w .eq. 'EDGE')  then
            itoken = 35
         elseif (w .eq. 'SCF ')  then
            itoken = 36
         elseif (w .eq. 'FMS ')  then
            itoken = 37
         elseif (w .eq. 'LDOS')  then
            itoken = 38
         elseif (w .eq. 'INTE')  then
            itoken = 39
         elseif (w .eq. 'CFAV')  then
            itoken = 40
         elseif (w .eq. 'S02 ')  then
            itoken = 41
         elseif (w .eq. 'XES ')  then
            itoken = 42
         elseif (w .eq. 'DANE')  then
            itoken = 43
         elseif (w .eq. 'FPRI')  then
            itoken = 44
         elseif (w .eq. 'RSIG')  then
            itoken = 45
         elseif (w .eq. 'XNCD')  then
            itoken = 46
         elseif (w .eq. 'XMCD')  then
            itoken = 46
         elseif (w .eq. 'MULT')  then
            itoken = 47
         elseif (w .eq. 'UNFR')  then
            itoken = 48
         elseif (w .eq. 'TDLD')  then
            itoken = 49
         elseif (w .eq. 'PMBS')  then
            itoken = 50
         elseif (w .eq. 'PLAS')  then
            itoken = 51
         elseif (w .eq. 'SO2C')  then
            itoken = 52
         elseif (w .eq. 'SELF')  then
            itoken = 53
         elseif (w .eq. 'SFSE')  then
            itoken = 54
         elseif (w .eq. 'RCONV') then
            itoken = 55
         elseif (w .eq. 'ELNE') then !KJ new card for EELS 1-06
            itoken = 56
         elseif (w .eq. 'EXEL') then !KJ new card for EELS 1-06
            itoken = 57
         elseif (w .eq. 'MAGI') then !KJ new card for EELS 1-06
            itoken = 58
         elseif (w .eq. 'ABSO') then !KJ new card 3-06
            itoken = 59    
         elseif (w .eq. 'EGRI')  then !Josh Kas
            itoken = 60
         elseif (w .eq. 'END ')  then
            itoken = -1            
         else
            itoken = 0
         endif
      elseif (flname(1:10).eq.'spring.inp') then
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     These tokens are for spring.inp (input for eq of motion method)
         if (w .eq. 'STRE')  then
            itoken = 1
         elseif (w .eq. 'ANGL')  then
            itoken = 2
         elseif (w .eq. 'VDOS')  then
            itoken = 3
         elseif (w .eq. 'PRDOS') then
            itoken = 4
         elseif (w .eq. 'END ')  then
            itoken = -1            
         else
            itoken = 0
         endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      endif
      return
      end


c====================================================================
      integer function nxtunt(iunit)

c  this function returns the value of the next unopened unit number
c  equal to or larger than iunit.  it will return neither unit numbers
c  0, 5, or 6 nor a negative unit number
c $Id: nxtunt.f,v 1.1.1.1 2006/01/12 06:37:42 hebhop Exp $
c $Log: nxtunt.f,v $
c Revision 1.1.1.1  2006/01/12 06:37:42  hebhop
c New version of feff. feff8.5 (Extension of feff8.4)
c Includes:
c 	1) All feff8.4 capabilities.
c 	2) Screened core hole (calculation of W).
c 	3) Multiple pole self energy calculation.
c 	4) Convolution with spectral function.
c New cards and options:
c 	1) NOHOLE 2      (screened hole)
c 	2) PLASMON ipl   (multiple pole self energy)
c 	3) SO2CONV       (convolve output with spectral function)
c 	4) SELF          (print on shell self energy as calculated by Luke)
c 	5) SFSE k0        (print off shell self energy Sigma(k0,e) )
c
c Revision 1.1.1.1  2000/02/11 02:23:58  alex
c Initialize feff82
c
c Revision 1.10  1999/04/02 21:32:47  newville
c cleaned up nxtunt (matt)
c
c Revision 1.9  1999/02/11 20:08:08  alex
c x39 version: dim.h + misc. small changes
c
c Revision 1.8  1998/12/29 23:59:07  alex
c feff8x35 version
c
c Revision 1.7  1998/11/19 03:23:11  alex
c feff8x32 version
c
c Revision 1.6  1998/10/26 14:11:16  ravel
c no comments beyond column 71
c
c Revision 1.5  1998/10/18 21:47:51  alex
c feff8x30 version implements Broyden algorithm for self-consistency
c
c Revision 1.4  1998/02/24 18:31:37  ravel
c I should really be more careful.  This is the last commitment done
c      cright.
c
c Revision 1.1.1.1  1997/04/27 20:18:03  ravel
c Initial import of xanes sources, version 0.37
c
c Revision 1.1  1996/06/23 16:05:02  bruce
c Initial revision
c

       integer iunit
       logical open

       nxtunt = max(1, iunit) - 1
 10    continue
       nxtunt = nxtunt + 1
       if ((nxtunt.eq.5).or.(nxtunt.eq.6)) nxtunt = 7
       inquire (unit=nxtunt, opened=open)
       if (open) go to 10
       return
c  end integer function nxtunt
       end

c====================================================================
c     Periodic table of the elements
c     Written by Steven Zabinsky, Feb 1992.  Deo Soli Gloria

c     atwts(iz)  single precision fn, returns atomic weight
c     atwtd(iz)  double precision fn, returns atomic weight
c     atsym(iz)  character*2 fn, returns atomic symbol

      double precision function atwtd (iz)
      double precision weight
      common /atwtco/ weight(103)
      atwtd = weight(iz)
      return
      end

      real function atwts (iz)
      double precision weight
      common /atwtco/ weight(103)
      atwts = weight(iz)
      return
      end

      character*2 function atsym (iz)
      character*2 sym
      common /atsyco/ sym(103)
      atsym = sym(iz)
      return
      end

      block data prtbbd
c     PeRiodic TaBle Block Data

c     Atomic weights from inside front cover of Ashcroft and Mermin.

      double precision weight
      common /atwtco/ weight(103)

      character*2 sym
      common /atsyco/ sym(103)

      data weight /
     1   1.0079, 4.0026, 6.941,  9.0122, 10.81,   12.01,
     2   14.007, 15.999, 18.998, 20.18,  22.9898, 24.305,
     3   26.982, 28.086, 30.974, 32.064, 35.453,  39.948,
     4   39.09,  40.08,  44.956, 47.90,  50.942,  52.00,
     5   54.938, 55.85,  58.93,  58.71,  63.55,   65.38,
     6   69.72,  72.59,  74.922, 78.96,  79.91,   83.80,
     7   85.47,  87.62,  88.91,  91.22,  92.91,   95.94,
     8   98.91,  101.07, 102.90, 106.40, 107.87,  112.40,
     9   114.82, 118.69, 121.75, 127.60, 126.90,  131.30,
     x   132.91, 137.34, 138.91, 140.12, 140.91,  144.24,
     1   145,    150.35, 151.96, 157.25, 158.92,  162.50,
     2   164.93, 167.26, 168.93, 173.04, 174.97,  178.49,
     3   180.95, 183.85, 186.2,  190.20, 192.22,  195.09,
     4   196.97, 200.59, 204.37, 207.19, 208.98,  210,
     5   210,    222,    223,    226,    227,     232.04,
     6   231,    238.03, 237.05, 244,    243,     247,
     7   247,    251,    254,    257,    256,     254,
     8   257/

      data sym /  'H', 'He','Li','Be','B', 'C', 'N', 'O', 'F', 'Ne',
     1            'Na','Mg','Al','Si','P', 'S', 'Cl','Ar','K', 'Ca',
     2            'Sc','Ti','V', 'Cr','Mn','Fe','Co','Ni','Cu','Zn',
     3            'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y', 'Zr',
     4            'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',
     5            'Sb','Te','I', 'Xe','Cs','Ba','La','Ce','Pr','Nd',
     6            'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
     7            'Lu','Hf','Ta','W', 'Te','Os','Ir','Pt','Au','Hg',
     8            'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',
     9            'Pa','U', 'Np','Pu','Am','Cm','Bk','Cf','Es','Fm',
     x            'Md','No','Lw'/

      end
      subroutine pijump (ph, old)
      implicit double precision (a-h, o-z)

c     removes jumps of 2*pi in phases

c     ph = current value of phase (may be modified on output, but
c          only by multiples of 2*pi)
c     old = previous value of phase

c={../HEADERS/const.h
      parameter (pi = 3.14159 26535 89793 23846 26433d0)
      parameter (one = 1, zero = 0)
      parameter (third = one/3)
      parameter (raddeg = 180 / pi)
      complex*16 coni
      parameter (coni = (0,1))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919 158 292 677 512 811d0)

      parameter (bohr = 0.529 177 249d0, ryd  = 13.605 698d0)
      parameter (hart = 2 * ryd)
      parameter (alpinv = 137.035 989 56d0)
c     fine structure alpha
      parameter (alphfs = 1 / alpinv)
c= ../HEADERS/const.h}
      parameter (twopi = 2 * pi)
      dimension xph(3)

      xph(1) = ph - old
      jump =  (abs(xph(1))+ pi) / twopi
      xph(2) = xph(1) - jump*twopi
      xph(3) = xph(1) + jump*twopi


      xphmin = min (abs(xph(1)), abs(xph(2)), abs(xph(3)))
      isave = 0
      do 10  i = 1, 3
         if (abs (xphmin - abs(xph(i))) .le. 0.01)  isave = i
   10 continue
      if (isave .eq. 0)  then
         call par_stop('pijump')
      endif

      ph = old + xph(isave)

      return
      end
      subroutine rdhead (io, nhead, head, lhead)
      implicit double precision (a-h, o-z)

c     Reads title line(s) from unit io.  Returns number of lines
c     read.  If more than nheadx lines, skips over them.  End-of-header
c     marker is a line of 1 blank, 71 '-'s.
c     lhead is length of each line w/o trailing blanks.
c     header lines returned will have 1st space on line blank for
c     carriage control

      character*80 head(nhead)
      dimension lhead(nhead)
      character*80  line

      n = 0
      nheadx = nhead
      nhead = 0
   10 read(io,20)  line
   20    format(a)
         if (line(4:11) .eq. '--------')  goto 100
         n = n+1
         if (n .le. nheadx)  then
            head(n) = line
            lhead(n) = istrln(head(n))
            nhead = n
         endif
      goto 10
  100 continue
      return
      end
      subroutine rdpot ( ntitle, title, rnrmav, xmu, vint, rhoint,
     1                  emu, s02, erelax, wp, ecv,rs,xf, qtotel, 
     2                  imt, rmt, inrm, rnrm, folp, folpx, xnatph,
     3                  dgc0, dpc0, dgc, dpc, adgc, adpc,
     3                  edens, vclap, vtot, edenvl, vvalgs, dmag, xnval,
     4                  eorb, kappa, iorb, qnrm, xnmues, nohole, ihole,
     5                  inters, totvol, iafolp, xion, iunf, iz, jumprm)
c  opens pot.bin file and reads following information
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
c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}
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

  10  format(a)
   20 format (bn, i15)

      open (unit=3, file='pot.bin', status='old')
      read(3,30) ntitle, nph, npadx, nohole, ihole, inters, iafolp,
     1            jumprm, iunf
  30  format(9(1x,i4))
c     nph and npadx are not passed to calling subroutine
      do 133  i  = 1, ntitle
         read(3,10) title(i)
         call triml(title(i))
  133 continue
c     Misc double precision stuff from pot.bin
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
      subroutine rdxsph ( ne, ne1, ne3, nph, ihole, rnrmav,xmu,edge,
     1               ik0, em, eref, iz, potlbl, ph, rkk, lmax, lmaxp1)
      implicit double precision (a-h, o-z)
c     reads file 'phase.bin' 
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
c     nph - number of potential types
c     iz  - charge of nuclei (atomic number)
c     potlbl - label for each potential type
c     lmax - max orb momentum for each potential type
c     ihole - index of core-hole orbital for absorber (iph=0)
c     rnrmav - average Norman radius (used in headers only)
c  Main output of xsect and phases module (except that in xsect.bin)
c     ph  - complex scattering phase shifts
c     rkk - complex multipole matrix elements

c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}

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

      open (unit=1, file='phase.bin', status='old', iostat=ios)
      call chopen (ios, 'phase.bin', 'rdxsph')

      read(1,10) nsp, ne, ne1, ne3, nph, ihole, ik0, npadx
  10  format (8(1x,i4))

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
      subroutine setkap(ihole, kinit, linit)
      implicit double precision (a-h, o-z)

c     Set initial state ang mom and quantum number kappa
c     ihole  initial state from ihole    
c     1      K    1s      L=0 -> linit=0 
c     2      LI   2s      L=0 -> linit=0
c     3      LII  2p 1/2  L=1 -> linit=1
c     4      LIII 2p 3/2  L=1 -> linit=1
c     5+     etc.
      if (ihole.le. 2 .or. ihole.eq. 5 .or. ihole.eq.10 .or.
     1    ihole.eq.17 .or. ihole.eq.24 .or. ihole.eq.27)  then
c        hole in s state
         linit = 0
         kinit = -1
      elseif (ihole.eq. 3 .or. ihole.eq. 6 .or. ihole.eq.11 .or.
     1        ihole.eq.18 .or. ihole.eq.25 .or. ihole.eq.30)  then
c        hole in p 1/2 state
         linit = 1
         kinit = 1
      elseif (ihole.eq. 4 .or. ihole.eq. 7 .or. ihole.eq.12 .or.
     1        ihole.eq.19 .or. ihole.eq.26)  then
c        hole in p 3/2 state
         linit = 1
         kinit = -2
      elseif (ihole.eq. 8 .or. ihole.eq.13 .or.
     1        ihole.eq.20 .or. ihole.eq.27)  then
c        hole in d 3/2 state
         linit = 2
         kinit = 2
      elseif (ihole.eq. 9 .or. ihole.eq.14 .or.
     1        ihole.eq.21 .or. ihole.eq.28)  then
c        hole in d 5/2 state
         linit = 2
         kinit = -3
      elseif (ihole.eq.15 .or. ihole.eq.22)  then
c        hole in  f 5/2 state
         linit = 3
         kinit = 3
      elseif (ihole.eq.16 .or. ihole.eq.23)  then
c        hole in  f 7/2 state
         linit = 3
         kinit = -4
      else
c        some unknown hole
         call par_stop('invalid hole number in setkap')
      endif

      return
      end
C FUNCTION ISTRLN (STRING)  Returns index of last non-blank
C                           character.  Returns zero if string is
C                           null or all blank.

      FUNCTION ISTRLN (STRING)
      CHARACTER*(*)  STRING
      CHARACTER BLANK, TAB
      PARAMETER (BLANK = ' ', TAB = '	')
C     there is a tab character here  ^

C  -- If null string or blank string, return length zero.
      ISTRLN = 0
      IF (STRING (1:1) .EQ. CHAR(0))  RETURN
      IF (STRING .EQ. ' ')  RETURN

C  -- Find rightmost non-blank character.
      ILEN = LEN (STRING)
      DO 20  I = ILEN, 1, -1
         IF (STRING(I:I).NE.BLANK .AND. STRING(I:I).NE.TAB)  GOTO 30
   20 CONTINUE
   30 ISTRLN = I

      RETURN
      END
C SUBROUTINE TRIML (STRING)  Removes leading blanks.

      SUBROUTINE TRIML (STRING)
      CHARACTER*(*)  STRING
      CHARACTER*200  TMP
      CHARACTER BLANK, TAB
      PARAMETER (BLANK = ' ', TAB = '	')
C     there is a tab character here  ^

      JLEN = ISTRLN (STRING)

C  -- All blank and null strings are special cases.
      IF (JLEN .EQ. 0)  RETURN

C  -- FInd first non-blank char
      DO 10  I = 1, JLEN
         IF (STRING(I:I).NE.BLANK .AND. STRING(I:I).NE.TAB)  GOTO 20
   10 CONTINUE
   20 CONTINUE

C  -- If I is greater than JLEN, no non-blanks were found.
      IF (I .GT. JLEN)  RETURN

C  -- Remove the leading blanks.
      TMP = STRING (I:)
      STRING = TMP
      RETURN
      END
C SUBROUTINE UPPER (STRING)  Changes a-z to upper case.

      SUBROUTINE UPPER (STRING)
      CHARACTER*(*)  STRING

      JLEN = ISTRLN (STRING)

      DO 10  I = 1, JLEN
         IC = ICHAR (STRING (I:I))
         IF ((IC .LT. 97)  .OR.  (IC .GT. 122))  GOTO 10
         STRING (I:I) = CHAR (IC - 32)
   10 CONTINUE

      RETURN
      END
C SUBROUTINE LOWER (STRING)  Changes A-Z to lower case.

      SUBROUTINE LOWER (STRING)
      CHARACTER*(*)  STRING

      JLEN = ISTRLN (STRING)

      DO 10  I = 1, JLEN
         IC = ICHAR (STRING (I:I))
         IF ((IC .LT. 65) .OR.  (IC .GT. 90))  GOTO 10
         STRING (I:I) = CHAR (IC + 32)
   10 CONTINUE

      RETURN
      END
C***********************************************************************
C
      SUBROUTINE BWORDS (S, NWORDS, WORDS)
C
C     Breaks string into words.  Words are seperated by one or more
C     blanks or tabs, or a comma and zero or more blanks.
C
C     ARGS        I/O      DESCRIPTION
C     ----        ---      -----------
C     S            I       CHAR*(*)  String to be broken up
C     NWORDS      I/O      Input:  Maximum number of words to get
C                          Output: Number of words found
C     WORDS(NWORDS) O      CHAR*(*) WORDS(NWORDS)
C                          Contains words found.  WORDS(J), where J is
C                          greater then NWORDS found, are undefined on
C                          output.
C
C      Written by:  Steven Zabinsky, September 1984
C      Tab char added July 1994.
C
C**************************  Deo Soli Gloria  **************************

C  -- No floating point numbers in this routine.
      IMPLICIT INTEGER (A-Z)

      CHARACTER*(*) S, WORDS(NWORDS)

      CHARACTER BLANK, COMMA, TAB
      PARAMETER (BLANK = ' ', COMMA = ',', TAB = '	')
C     there is a tab character here               ^.

C  -- BETW    .TRUE. if between words
C     COMFND  .TRUE. if between words and a comma has already been found
      LOGICAL BETW, COMFND

C  -- Maximum number of words allowed
      WORDSX = NWORDS

C  -- SLEN is last non-blank character in string
      SLEN = ISTRLN (S)

C  -- All blank string is special case
      IF (SLEN .EQ. 0)  THEN
         NWORDS = 0
         RETURN
      ENDIF

C  -- BEGC is beginning character of a word
      BEGC = 1
      NWORDS = 0

      BETW   = .TRUE.
      COMFND = .TRUE.

      DO 10  I = 1, SLEN
         IF (S(I:I) .EQ. BLANK .OR. S(I:I) .EQ. TAB)  THEN
            IF (.NOT. BETW)  THEN
               NWORDS = NWORDS + 1
               WORDS (NWORDS) = S (BEGC : I-1)
               BETW = .TRUE.
               COMFND = .FALSE.
            ENDIF
         ELSEIF (S(I:I) .EQ. COMMA)  THEN
            IF (.NOT. BETW)  THEN
               NWORDS = NWORDS + 1
               WORDS (NWORDS) = S(BEGC : I-1)
               BETW = .TRUE.
            ELSEIF (COMFND)  THEN
               NWORDS = NWORDS + 1
               WORDS (NWORDS) = BLANK
            ENDIF
            COMFND = .TRUE.
         ELSE
            IF (BETW)  THEN
               BETW = .FALSE.
               BEGC = I
            ENDIF
         ENDIF

         IF (NWORDS .GE. WORDSX)  RETURN

   10 CONTINUE

      IF (.NOT. BETW  .AND.  NWORDS .LT. WORDSX)  THEN
         NWORDS = NWORDS + 1
         WORDS (NWORDS) = S (BEGC :SLEN)
      ENDIF

      RETURN
      END
      SUBROUTINE UNTAB (STRING)
C REPLACE TABS WITH BLANKS :    TAB IS ASCII DEPENDENT
      INTEGER        ITAB , I, ILEN, ISTRLN
      PARAMETER      (ITAB = 9)
      CHARACTER*(*)  STRING, TAB*1
      EXTERNAL ISTRLN
      TAB  = CHAR(ITAB)
      ILEN = MAX(1, ISTRLN(STRING))
 10   CONTINUE 
        I = INDEX(STRING(:ILEN), TAB ) 
        IF (I .NE. 0) THEN
            STRING(I:I) = ' '
            GO TO 10
        END IF
      RETURN
C END SUBROUTINE UNTAB
      END

      logical function iscomm (line)
c     returns true if line is a comment or blank line, false otherwise
c#mn{ rewritten to allow ";*%#" as comment characters
       character*(*) line
       iscomm = ((line.eq.' ').or.(index(';*%#',line(1:1)).ne.0))
c#mn}
      return
      end
      subroutine str2dp(str,dpval,ierr)
c  return dp number "dpval" from character string "str"
c  if str cannot be a number, ierr < 0 is returned.
      character*(*) str, fmt*15 
      double precision dpval
      integer  ierr , lenmax
      parameter ( lenmax = 40)
      logical  isnum
      external isnum
      ierr = -99
      if (isnum(str)) then
         ierr = 0
         write(fmt, 10) min(lenmax, len(str))
 10      format('(bn,f',i3,'.0)')
         read(str, fmt, err = 20, iostat=ierr) dpval
      end if    
      if (ierr.gt.0) ierr = -ierr
      return
 20   continue
      ierr = -98
      return
c end subroutine str2dp
      end
      subroutine str2re(str,val,ierr)
c  return real from character string "str"
      character*(*) str 
      double precision dpval
      real     val
      integer  ierr
      call str2dp(str,dpval,ierr)
      if (ierr.eq.0) val = dpval
      return
c end subroutine str2re
      end
      subroutine str2in(str,intg,ierr)
c  return integer from character string "str"
c  returns ierr = 1 if value was clearly non-integer
      character*(*) str 
      double precision val, tenth
      parameter (tenth = 1.d-1)
      integer  ierr, intg
      call str2dp(str,val,ierr)
      if (ierr.eq.0) then
         intg = int(val)
         if ((abs(intg - val) .gt. tenth))  ierr = 1
       end if
      return
c end subroutine str2in
      end
       logical function isnum (string)
c  tests whether a string can be a number. not foolproof!
c  to return true, string must contain:
c    - only characters in  'deDE.+-, 1234567890' (case is checked)
c    - no more than one 'd' or 'e' 
c    - no more than one '.'
c  matt newville
       character*(*)  string,  number*20
c note:  layout and case of *number* is important: do not change!
       parameter (number = 'deDE.,+- 1234567890')
       integer   iexp, idec, i, j, istrln
       external  istrln
       iexp  = 0
       idec  = 0
       isnum = .false. 
       do 100  i = 1, max(1, istrln(string))
          j = index(number,string(i:i))
          if (j.le.0)               go to 200
          if((j.ge.1).and.(j.le.4)) iexp = iexp + 1
          if (j.eq.5)               idec = idec + 1
 100   continue
c  every character in "string" is also in "number".  so, if there are 
c  not more than one exponential and decimal markers, it's a number
       if ((iexp.le.1).and.(idec.le.1)) isnum = .true.
 200   continue
       return
c  end logical function isnum
       end
      subroutine wlog (string)
      character*(*) string

c={../HEADERS/parallel.h
      integer par_type, this_process, numprocs, my_rank
      logical master, worker, parallel_run
      real*8 wall_comm, time_comm
      common /timing/ wall_comm, time_comm
      common /parallel/ numprocs, my_rank, this_process, 
     .          master, worker, parallel_run, par_type
c= ../HEADERS/parallel.h}

c     This output routine is used to replace the PRINT statement
c     for output that "goes to the terminal", or to the log file.
c     If you use a window based system, you can modify this routine
c     to handle the running output elegantly.
c     Handle carriage control in the string you pass to wlog.
c
c     The log file is also written here, hard coded here.

c     The log file is unit 11.  The log file is opened in the
c     main program, program feff.

c     make sure not to write trailing blanks

   10 format (a)

c     Suppress output in sequential loops
      if (par_type .eq. 2) return

      il = istrln (string)
      if (il .eq. 0)  then
         print10
         if (par_type .ne. 3) write(11,10)
      else
         print10, string(1:il)
         if (par_type .ne. 3) write(11,10) string(1:il)
      endif
      return
      end
      subroutine lblank (string)
      character*(*) string
c     add a leading blank, useful for carriage control
      string = ' ' // string
      return
      end
      double precision function xx (j)
      implicit double precision (a-h, o-z)
c     x grid point at index j, x = log(r), r=exp(x)
      parameter (delta = 0.050 000 000 000 000)
      parameter (c88   = 8.800 000 000 000 000)
c     xx = -8.8 + (j-1)*0.05
      xx = -c88 + (j-1)*delta
      return
      end

      double precision function rr(j)
      implicit double precision (a-h, o-z)
c     r grid point at index j
      rr = exp (xx(j))
      return
      end

      function ii(r)
      implicit double precision (a-h, o-z)
c     index of grid point immediately below postion r
      parameter (delta = 0.050 000 000 000 000)
      parameter (c88   = 8.800 000 000 000 000)
c     ii = (log(r) + 8.8) / 0.05 + 1
      ii = (log(r) + c88) / delta + 1
      return
      end
c
c PAD library:   Packed Ascii Data 
c   these routines contain code for handling packed-ascii-data  
c   (pad) arrays for writing printable character strings that 
c   represent real or complex scalars and arrays to a file.
c
c routines included in padlib are (dp==double precision):
c   wrpadd     write a dp array as pad character strings
c   wrpadx     write a dp complex array as pad character strings
c   rdpadr     read a pad character array as a real array
c   rdpadd     read a pad character array as a dp  array
c   rdpadc     read a pad character array as a complex array
c   rdpadx     read a pad character array as a dp complex array
c   pad        internal routine to convert dp number to pad string
c   unpad      internal routine to pad string to dp number
c
c routines not included, but required by padlib:
c     triml, istrln, wlog
c
c//////////////////////////////////////////////////////////////////////
c Copyright (c) 1997--2001 Matthew Newville, The University of Chicago
c Copyright (c) 1992--1996 Matthew Newville, University of Washington
c
c Permission to use and redistribute the source code or binary forms of
c this software and its documentation, with or without modification is
c hereby granted provided that the above notice of copyright, these
c terms of use, and the disclaimer of warranty below appear in the
c source code and documentation, and that none of the names of The
c University of Chicago, The University of Washington, or the authors
c appear in advertising or endorsement of works derived from this
c software without specific prior written permission from all parties.
c
c THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
c EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
c MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
c IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
c CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
c TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
c SOFTWARE OR THE USE OR OTHER DEALINGS IN THIS SOFTWARE.
c//////////////////////////////////////////////////////////////////////
c License is applicable for routines below, until otherwise specified.
c
       subroutine wrpadd(iout,npack,array,npts)
c
c write a dp array to a file in packed-ascii-data format
c
c inputs:  [ no outputs / no side effects ]
c   iout   unit to write to (assumed open)
c   npack  number of characters to use (determines precision)
c   array  real array 
c   npts   number of array elements to read
c notes:
c   real number converted to packed-ascii-data string using pad
c={padlib.h
c padlib.h -*-fortran-*-
c  header of parameters for packed-ascii-data (pad) routines
       implicit none
       character cpadr, cpadi, cpadc
       integer   maxlen, ibase, ioff, ihuge, ibas2
       double precision  ten, tenlog, huge, tiny, one, zero, base
       parameter(cpadr = '!', cpadc = '$', cpadi = '%')
       parameter(ibase = 90, ioff = 37, ihuge = 38, maxlen = 82)
       parameter(ibas2 = ibase/2, zero=0d0, one=1.d0, ten = 10.d0)
       parameter(tenlog= 2.302585092994045684d0)
       parameter(huge = ten**ihuge, tiny = one/huge)
       parameter(base = ibase*one)
c
c= padlib.h}
       integer    iout, npack, npts, mxl, js, i
       character  str*128
       double precision array(*), xr
       js  = 0
       str = ' '
       mxl = maxlen - npack + 1
       do 20 i = 1, npts
          js = js+npack
          xr = array(i)
          call pad(xr, npack, str(js-npack+1:js))
          if ((js.ge.mxl).or.(i.eq.npts)) then
             write(iout,100) cpadr, str(1:js)
             js = 0
          end if
 20    continue
       return
 100   format(a1,a)
       end
c --padlib--
       subroutine wrpadx(iout,npack,array,npts)
c write complex*16 array as pad string
c={padlib.h
c padlib.h -*-fortran-*-
c  header of parameters for packed-ascii-data (pad) routines
       implicit none
       character cpadr, cpadi, cpadc
       integer   maxlen, ibase, ioff, ihuge, ibas2
       double precision  ten, tenlog, huge, tiny, one, zero, base
       parameter(cpadr = '!', cpadc = '$', cpadi = '%')
       parameter(ibase = 90, ioff = 37, ihuge = 38, maxlen = 82)
       parameter(ibas2 = ibase/2, zero=0d0, one=1.d0, ten = 10.d0)
       parameter(tenlog= 2.302585092994045684d0)
       parameter(huge = ten**ihuge, tiny = one/huge)
       parameter(base = ibase*one)
c
c= padlib.h}
       integer    iout, npack, npts, mxl, js, i
       complex*16 array(*)
       character  str*128
       double precision xr, xi
       js = 0
       str  = ' '
       mxl  = maxlen - 2 * npack + 1
       do 20 i = 1, npts
          js = js  + 2 * npack
          xr = dble(array(i))
          xi = dimag(array(i))
          call pad(xr, npack, str(js-2*npack+1:js-npack))
          call pad(xi, npack, str(js-npack+1:js))
          if ((js.ge.mxl).or.(i.eq.npts)) then
             write(iout,100) cpadc, str(1:js)
             js = 0
          end if
 20    continue
       return
 100   format(a1,a)
       end
c --padlib--
       subroutine wrpadr(iout,npack,array,npts)
c
c write a real array to a file in packed-ascii-data format
c
c inputs:  [ no outputs / no side effects ]
c   iout   unit to write to (assumed open)
c   npack  number of characters to use (determines precision)
c   array  real array 
c   npts   number of array elements to read
c notes:
c   real number converted to packed-ascii-data string using pad
c={padlib.h
c padlib.h -*-fortran-*-
c  header of parameters for packed-ascii-data (pad) routines
       implicit none
       character cpadr, cpadi, cpadc
       integer   maxlen, ibase, ioff, ihuge, ibas2
       double precision  ten, tenlog, huge, tiny, one, zero, base
       parameter(cpadr = '!', cpadc = '$', cpadi = '%')
       parameter(ibase = 90, ioff = 37, ihuge = 38, maxlen = 82)
       parameter(ibas2 = ibase/2, zero=0d0, one=1.d0, ten = 10.d0)
       parameter(tenlog= 2.302585092994045684d0)
       parameter(huge = ten**ihuge, tiny = one/huge)
       parameter(base = ibase*one)
c
c= padlib.h}
       integer    iout, npack, npts, mxl, js, i
       character  str*128
       real    array(*)
       double precision xr
       js  = 0
       str = ' '
       mxl = maxlen - npack + 1
       do 20 i = 1, npts
          js = js+npack
          xr = dble(array(i))
          call pad(xr, npack, str(js-npack+1:js))
          if ((js.ge.mxl).or.(i.eq.npts)) then
             write(iout,100) cpadr, str(1:js)
             js = 0
          end if
 20    continue
       return
 100   format(a1,a)
       end
c --padlib--
       subroutine wrpadc(iout,npack,array,npts)
c write complex (*8) array as pad string
c={padlib.h
c padlib.h -*-fortran-*-
c  header of parameters for packed-ascii-data (pad) routines
       implicit none
       character cpadr, cpadi, cpadc
       integer   maxlen, ibase, ioff, ihuge, ibas2
       double precision  ten, tenlog, huge, tiny, one, zero, base
       parameter(cpadr = '!', cpadc = '$', cpadi = '%')
       parameter(ibase = 90, ioff = 37, ihuge = 38, maxlen = 82)
       parameter(ibas2 = ibase/2, zero=0d0, one=1.d0, ten = 10.d0)
       parameter(tenlog= 2.302585092994045684d0)
       parameter(huge = ten**ihuge, tiny = one/huge)
       parameter(base = ibase*one)
c
c= padlib.h}
       integer    iout, npack, npts, mxl, js, i
       complex    array(*)
       character  str*128
       double precision xr, xi
       js = 0
       str  = ' '
       mxl  = maxlen - 2 * npack + 1
       do 20 i = 1, npts
          js = js  + 2 * npack
          xr = dble(array(i))
          xi = aimag(array(i))
          call pad(xr, npack, str(js-2*npack+1:js-npack))
          call pad(xi, npack, str(js-npack+1:js))
          if ((js.ge.mxl).or.(i.eq.npts)) then
             write(iout,100) cpadc, str(1:js)
             js = 0
          end if
 20    continue
       return
 100   format(a1,a)
       end
c --padlib--
       subroutine rdpadd(iou,npack,array,npts)
c read dparray from packed-ascii-data file
c arguments:
c   iou    unit to read from (assumed open)                   (in)
c   npack  number of characters to use (determines precision) (in)
c   array  real array                                         (out)
c   npts   number of array elements to read / number read     (in/out)
c notes:
c   packed-ascii-data string converted to real array using  unpad
c={padlib.h
c padlib.h -*-fortran-*-
c  header of parameters for packed-ascii-data (pad) routines
       implicit none
       character cpadr, cpadi, cpadc
       integer   maxlen, ibase, ioff, ihuge, ibas2
       double precision  ten, tenlog, huge, tiny, one, zero, base
       parameter(cpadr = '!', cpadc = '$', cpadi = '%')
       parameter(ibase = 90, ioff = 37, ihuge = 38, maxlen = 82)
       parameter(ibas2 = ibase/2, zero=0d0, one=1.d0, ten = 10.d0)
       parameter(tenlog= 2.302585092994045684d0)
       parameter(huge = ten**ihuge, tiny = one/huge)
       parameter(base = ibase*one)
c
c= padlib.h}
       integer iou, npack, npts, ndline, i, istrln, ipts, iread
       double precision    array(*), unpad , tmp
       character  ctest, ccomp
       character  str*128
       external  unpad, istrln, iread
       ccomp = cpadr
       ipts = 0
 10    continue 
          i = iread(iou, str)
          if (i.lt.0) go to 50
          call triml(str)
          ctest  = str(1:1)
          str    = str(2:)
          ndline = i/npack
          if ((ctest.ne.ccomp).or.(ndline.le.0)) go to 200
          do 30 i = 1, ndline
             ipts  = ipts + 1
             tmp   = unpad(str(1-npack+i*npack:i*npack),npack)
             array(ipts) = tmp
             if (ipts.ge.npts) go to 50
 30       continue 
          go to 10
 50    continue 
       return
 200   continue
       call wlog (' -- Read_PAD error:  bad data at line:')
       i = istrln(str)
       call wlog (str(:i))
       stop ' -- fatal error in reading PAD data file -- '
       end
c --padlib--
       subroutine rdpadr(iou,npack,array,npts)
c read real array from packed-ascii-data file
c arguments:
c   iou    unit to read from (assumed open)                   (in)
c   npack  number of characters to use (determines precision) (in)
c   array  real array                                         (out)
c   npts   number of array elements to read / number read     (in/out)
c notes:
c   packed-ascii-data string converted to real array using  unpad
c={padlib.h
c padlib.h -*-fortran-*-
c  header of parameters for packed-ascii-data (pad) routines
       implicit none
       character cpadr, cpadi, cpadc
       integer   maxlen, ibase, ioff, ihuge, ibas2
       double precision  ten, tenlog, huge, tiny, one, zero, base
       parameter(cpadr = '!', cpadc = '$', cpadi = '%')
       parameter(ibase = 90, ioff = 37, ihuge = 38, maxlen = 82)
       parameter(ibas2 = ibase/2, zero=0d0, one=1.d0, ten = 10.d0)
       parameter(tenlog= 2.302585092994045684d0)
       parameter(huge = ten**ihuge, tiny = one/huge)
       parameter(base = ibase*one)
c
c= padlib.h}
       integer iou, npack, npts, ndline, i, istrln, ipts, iread
       real    array(*)
       double precision unpad , tmp
       character  ctest, ccomp
       character  str*128
       external  unpad, istrln, iread
       ccomp = cpadr
       ipts = 0
 10    continue 
          i = iread(iou, str)
          if (i.lt.0) go to 50
          call triml(str)
          ctest  = str(1:1)
          str    = str(2:)
          ndline = i/npack
          if ((ctest.ne.ccomp).or.(ndline.le.0)) go to 200
          do 30 i = 1, ndline
             ipts  = ipts + 1
             tmp   = unpad(str(1-npack+i*npack:i*npack),npack)
             array(ipts) = real(tmp)
             if (ipts.ge.npts) go to 50
 30       continue 
          go to 10
 50    continue 
       return
 200   continue
       call wlog (' -- Read_PAD error:  bad data at line:')
       i = istrln(str)
       call wlog (str(:i))
       stop ' -- fatal error in reading PAD data file -- '
       end
c --padlib--
       subroutine rdpadc(iou,npack,array,npts)
c read complex array from packed-ascii-data file
c arguments:
c   iou    unit to read from (assumed open)                  (in)
c   npack  number of characters to use (determines precision)(in)
c   array  complex array                                     (out)
c   npts   number of array elements to read / number read    (in/out)
c notes:
c   packed-ascii-data string converted to real array using  unpad
c={padlib.h
c padlib.h -*-fortran-*-
c  header of parameters for packed-ascii-data (pad) routines
       implicit none
       character cpadr, cpadi, cpadc
       integer   maxlen, ibase, ioff, ihuge, ibas2
       double precision  ten, tenlog, huge, tiny, one, zero, base
       parameter(cpadr = '!', cpadc = '$', cpadi = '%')
       parameter(ibase = 90, ioff = 37, ihuge = 38, maxlen = 82)
       parameter(ibas2 = ibase/2, zero=0d0, one=1.d0, ten = 10.d0)
       parameter(tenlog= 2.302585092994045684d0)
       parameter(huge = ten**ihuge, tiny = one/huge)
       parameter(base = ibase*one)
c
c= padlib.h}
       integer iou, npack,npts, ndline, i, istrln, ipts, np, iread
       double precision  unpad, tmpr, tmpi
       complex  array(*)
       character  ctest, ccomp
       character  str*128
       external  unpad, istrln, iread
       ccomp = cpadc
       ipts = 0
       np   = 2 * npack
 10    continue 
          i = iread(iou, str)
          if (i.lt.0) go to 50
          call triml(str)
          ctest  = str(1:1)
          str    = str(2:)
          ndline = i / np
          if ((ctest.ne.ccomp).or.(ndline.le.0)) go to 200
          do 30 i = 1, ndline
             ipts = ipts + 1
             tmpr = unpad(str(1-np+i*np:-npack+i*np),npack)
             tmpi = unpad(str(1-npack+i*np:i*np),npack)
             array(ipts) = cmplx(tmpr, tmpi)
             if (ipts.ge.npts) go to 50
 30       continue 
          go to 10
 50    continue 
       return
 200   continue
       call wlog (' -- Read_PAD error:  bad data at line:')
       i = istrln(str)
       call wlog (str(:i))
       stop ' -- fatal error in reading PAD data file -- '
       end
       subroutine rdpadx(iou,npack,array,npts)
c read complex*16 array from packed-ascii-data file
c arguments:
c   iou    unit to read from (assumed open)                  (in)
c   npack  number of characters to use (determines precision)(in)
c   array  complex array                                     (out)
c   npts   number of array elements to read / number read    (in/out)
c notes:
c   packed-ascii-data string converted to real array using  unpad
c={padlib.h
c padlib.h -*-fortran-*-
c  header of parameters for packed-ascii-data (pad) routines
       implicit none
       character cpadr, cpadi, cpadc
       integer   maxlen, ibase, ioff, ihuge, ibas2
       double precision  ten, tenlog, huge, tiny, one, zero, base
       parameter(cpadr = '!', cpadc = '$', cpadi = '%')
       parameter(ibase = 90, ioff = 37, ihuge = 38, maxlen = 82)
       parameter(ibas2 = ibase/2, zero=0d0, one=1.d0, ten = 10.d0)
       parameter(tenlog= 2.302585092994045684d0)
       parameter(huge = ten**ihuge, tiny = one/huge)
       parameter(base = ibase*one)
c
c= padlib.h}
       integer iou, npack,npts, ndline, i, istrln, ipts, np, iread
       double precision  unpad, tmpr, tmpi
       complex*16  array(*)
       character  ctest, ccomp
       character  str*128
       external  unpad, istrln, iread
       ccomp = cpadc
       ipts = 0
       np   = 2 * npack
 10    continue 
          i = iread(iou, str)
          if (i.lt.0) go to 50
          call triml(str)
          ctest  = str(1:1)
          str    = str(2:)
          ndline = i / np
          if ((ctest.ne.ccomp).or.(ndline.le.0)) go to 200
          do 30 i = 1, ndline
             ipts = ipts + 1
             tmpr = unpad(str(1-np+i*np:-npack+i*np),npack)
             tmpi = unpad(str(1-npack+i*np:i*np),npack)
             array(ipts) = cmplx(tmpr, tmpi)
             if (ipts.ge.npts) go to 50
 30       continue 
          go to 10
 50    continue 
       return
 200   continue
       call wlog (' -- Read_PAD error:  bad data at line:')
       i = istrln(str)
       call wlog (str(:i))
       stop ' -- fatal error in reading PAD data file -- '
       end

c --padlib--
       subroutine pad(xreal,npack,str)
c  convert dp number *xreal* to packed-ascii-data string *str*
c={padlib.h
c padlib.h -*-fortran-*-
c  header of parameters for packed-ascii-data (pad) routines
       implicit none
       character cpadr, cpadi, cpadc
       integer   maxlen, ibase, ioff, ihuge, ibas2
       double precision  ten, tenlog, huge, tiny, one, zero, base
       parameter(cpadr = '!', cpadc = '$', cpadi = '%')
       parameter(ibase = 90, ioff = 37, ihuge = 38, maxlen = 82)
       parameter(ibas2 = ibase/2, zero=0d0, one=1.d0, ten = 10.d0)
       parameter(tenlog= 2.302585092994045684d0)
       parameter(huge = ten**ihuge, tiny = one/huge)
       parameter(base = ibase*one)
c
c= padlib.h}
       integer  iexp, itmp, isgn, i, npack, j
       double precision xreal, xwork, xsave,onem, tenth
       parameter (onem  =  0.99999999997d0)
       parameter (tenth =  0.099999999994d0)
       character str*(*)
c
       str      = ' '
       xsave    = min(huge, max(-huge, xreal))
       isgn     = 1
       if (xsave.le.0) isgn = 0
c
       xwork    = dabs( xsave )
       iexp     = 0
       if ((xwork.lt.huge).and.(xwork.gt.tiny))  then
          iexp  =   1 + int(log(xwork) / tenlog  )
       else if (xwork.ge.huge) then
          iexp  = ihuge
          xwork = one
       else if (xwork.le.tiny)  then
          xwork = zero
       end if
c force xwork between ~0.1 and ~1
c note: this causes a loss of precision, but 
c allows backward compatibility
       xwork    = xwork / (ten ** iexp)
 20    continue
       if (xwork.ge.one) then
          xwork = xwork * 0.100000000000000d0
          iexp  = iexp + 1
       else if (xwork.le.tenth) then
          xwork = xwork * ten
          iexp  = iexp - 1
       endif
       if (xwork.ge.one) go to 20

       itmp     = int ( ibas2 * xwork ) 
       str(1:1) = char(iexp  + ioff + ibas2 )
       str(2:2) = char( 2 * itmp + isgn + ioff)
       xwork    = xwork * ibas2 - itmp
       if (npack.gt.2) then
          do 100 i = 3, npack
             itmp     = int( base * xwork + 1.d-9)
             str(i:i) = char(itmp + ioff)
             xwork    = xwork * base - itmp
 100      continue
       end if
       if (xwork.ge.0.5d0) then
          i = itmp + ioff + 1
          if (i.le.126) then
             str(npack:npack)= char(i)
          else 
             j = ichar(str(npack-1:npack-1))
             if (j.lt.126) then
                str(npack-1:npack-1) = char(j+1)
                str(npack:npack)     = char(37)
             endif 
          endif
       endif
       return
       end
c --padlib--
       double precision function unpad(str,npack)
c
c  convert packed-ascii-data string *str* to dp number *unpad*
c={padlib.h
c padlib.h -*-fortran-*-
c  header of parameters for packed-ascii-data (pad) routines
       implicit none
       character cpadr, cpadi, cpadc
       integer   maxlen, ibase, ioff, ihuge, ibas2
       double precision  ten, tenlog, huge, tiny, one, zero, base
       parameter(cpadr = '!', cpadc = '$', cpadi = '%')
       parameter(ibase = 90, ioff = 37, ihuge = 38, maxlen = 82)
       parameter(ibas2 = ibase/2, zero=0d0, one=1.d0, ten = 10.d0)
       parameter(tenlog= 2.302585092994045684d0)
       parameter(huge = ten**ihuge, tiny = one/huge)
       parameter(base = ibase*one)
c
c= padlib.h}
       double precision sum
       integer   iexp, itmp, isgn, i, npack
       character str*(*)
       unpad = zero
       if (npack.le.2) return
       iexp  =     (ichar(str(1:1)) - ioff   ) - ibas2
       isgn  = mod (ichar(str(2:2)) - ioff, 2) * 2 - 1
       itmp  =     (ichar(str(2:2)) - ioff   ) / 2
       sum   = dble(itmp/(base*base))
c       do 100 i = 3, npack
c          sum = sum + dble(ichar(str(i:i)) - ioff) / base**i
c 100   continue
       do 100 i = npack, 3, -1
          sum = sum + dble(ichar(str(i:i)) - ioff) / base**i
 100   continue
       unpad = 2 * isgn * ibase * sum * (ten ** iexp)
cc       print*, sum, iexp,unpad
       return
       end
c --padlib--
c end of pad library
c ----------
       integer function iread(lun,string)
c
c generalized internal read:
c    read a string the next line of an opened file 
c    unit, returning the real length of string
c 
c inputs:   
c   lun     opened file unit number
c outputs:
c   string  string read from file
c returns:
c   iread   useful length of string, as found from 
c                  sending string to 'sclean' to 
c                  remove non-printable characters
c                   and then istrln  
c           or
c              -1   on 'end-of-file'
c              -2   on 'error'
c
c copyright (c) 1999  Matthew Newville
       implicit none
       character*(*) string
       integer    lun, istrln
       external   istrln
       string = ' '
 10    format(a)
       read(lun, 10, end = 40, err = 50) string
       call sclean(string)
       iread = istrln(string)
       return
 40    continue 
       string = ' '
       iread = -1
       return
 50    continue 
       string = ' '
       iread = -2
       return
       end
       subroutine sclean(str) 
c
c  clean a string, especially for strings passed between 
c  different file systems, or from C functions:
c
c   1. characters in the range char(0), or char(10)...char(15) 
c      are interpreted as end-of-line characters, so that all
c      remaining characters are explicitly blanked.
c   2. all other characters below char(31) (including tab) are
c      replaced by a single blank
c
c  this is mostly useful when getting a string generated by a C 
c  function and for handling dos/unix/max line-endings.
c
c copyright (c) 1999  Matthew Newville
       character*(*) str, blank*1
       parameter (blank = ' ')
       integer i,j,is
       do 20 i = 1, len(str)
          is = ichar(str(i:i))
          if ((is.eq.0) .or. ((is.ge.10) .and. (is.le.15))) then
             do 10 j= i, len(str)
                str(j:j) = blank
 10          continue
             return
          elseif (is.le.31)  then
             str(i:i)  = blank
          end if
 20    continue 
       return
c end subroutine sclean
       end

      SUBROUTINE rdcmt(iUnt,Cmt)
      INTEGER iUnt, i1
      CHARACTER(300) line
      CHARACTER(4) Cmt
      CHARACTER TmpCmt(4), ch
      LOGICAL CmtLin

      CmtLin = .true.
      DO i1 = 1, 4
         TmpCmt(i1) = Cmt(i1:i1)
      END DO
 5    CONTINUE
      READ(iUnt,*,END=10) ch
      DO i1 = 1, 4
         IF(ch.eq.TmpCmt(i1)) goto 5
      END DO
      
      BACKSPACE(iUnt)
      
 10   CONTINUE
      
      RETURN
      END
      subroutine setgam (iz, ihole, gamach)

c     Sets gamach, core hole lifetime.  Data comes from graphs in
c     K. Rahkonen and K. Krause,
c     Atomic Data and Nuclear Data Tables, Vol 14, Number 2, 1974.
c     output gamach is in eV

      implicit double precision (a-h, o-z)

      dimension zh(8,16), gamh(8,16)

      dimension zk(8), gamkp(8)
      parameter (ryd  = 13.605 698d0)
      parameter (hart = 2*ryd)
      character*512 slog


c     Note that 0.99 replaces 1.0, 95.1 replaces 95.0 to avoid roundoff
c     trouble.
c     Gam arrays contain the gamma values.
c     We will take log10 of the gamma values so we can do linear
c     interpolation from a log plot.

      data  zh   / 0.99,  10.0, 20.0, 40.0, 50.0, 60.0, 80.0, 95.1,
     2              0.99, 18.0, 22.0, 35.0, 50.0, 52.0, 75.0,  95.1,
     3              0.99,  17.0, 28.0, 31.0, 45.0, 60.0,  80.0, 95.1,
     4              0.99,  17.0, 28.0, 31.0, 45.0, 60.0,  80.0, 95.1,
     5              0.99,  20.0, 28.0, 30.0, 36.0, 53.0,  80.0, 95.1,
     6              0.99,  20.0, 22.0, 30.0, 40.0, 68.0,  80.0, 95.1,
     7              0.99,  20.0, 22.0, 30.0, 40.0, 68.0,  80.0, 95.1,
     8              0.99,  36.0, 40.0, 48.0, 58.0, 76.0,  79.0, 95.1,
     9              0.99,  36.0, 40.0, 48.0, 58.0, 76.0,  79.0, 95.1,
     *              0.99,  30.0, 40.0, 47.0, 50.0, 63.0,  80.0, 95.1,
     1              0.99,  40.0, 42.0, 49.0, 54.0, 70.0,  87.0, 95.1,
     2              0.99,  40.0, 42.0, 49.0, 54.0, 70.0,  87.0, 95.1,
     3              0.99,  40.0, 50.0, 55.0, 60.0, 70.0,  81.0, 95.1,
     4              0.99,  40.0, 50.0, 55.0, 60.0, 70.0,  81.0, 95.1,
     5              0.99,  71.0, 73.0, 79.0, 86.0, 90.0,  95.0,100.0,
     6              0.99,  71.0, 73.0, 79.0, 86.0, 90.0,  95.0,100.0/

      data  gamh / 0.02,  0.28, 0.75,  4.8, 10.5, 21.0, 60.0, 105.0,
     2              0.07,  3.9,  3.8,  7.0,  6.0,  3.7,  8.0,  19.0,
     3              0.001, 0.12,  1.4,  0.8,  2.6,  4.1,   6.3, 10.5,
     4              0.001, 0.12, 0.55,  0.7,  2.1,  3.5,   5.4,  9.0,
     5              0.001,  1.0,  2.9,  2.2,  5.5, 10.0,  22.0, 22.0,
     6              0.001,0.001,  0.5,  2.0,  2.6, 11.0,  15.0, 16.0,
     7              0.001,0.001,  0.5,  2.0,  2.6, 11.0,  10.0, 10.0,
     8              0.0006,0.09, 0.07, 0.48,  1.0,  4.0,   2.7,  4.7,
     9              0.0006,0.09, 0.07, 0.48, 0.87,  2.2,   2.5,  4.3,
     *              0.001,0.001,  6.2,  7.0,  3.2, 12.0,  16.0, 13.0,
     1              0.001,0.001,  1.9, 16.0,  2.7, 13.0,  13.0,  8.0,
     2              0.001,0.001,  1.9, 16.0,  2.7, 13.0,  13.0,  8.0,
     3              0.001,0.001, 0.15,  0.1,  0.8,  8.0,   8.0,  5.0,
     4              0.001,0.001, 0.15,  0.1,  0.8,  8.0,   8.0,  5.0,
     5              0.001,0.001, 0.05, 0.22,  0.1, 0.16,   0.5,  0.9,
     6              0.001,0.001, 0.05, 0.22,  0.1, 0.16,   0.5,  0.9/

c     Since feff8 can be called any number of times . ALA

      if (ihole .le. 0)  then
         gamach = 0
         write(slog,'(a,1pe13.5)') ' No hole in SETGAM, gamach = ', 
     1                             gamach
         call wlog(slog)
         return
      endif
      if (ihole .gt. 16)  then
         call wlog(' This version of FEFF will set gamach = 0.1 eV ' //
     1             ' for O1 and higher hole')
         call wlog(' You can use CORRECTIONS card  to set ' //
     1   ' gamach = 0.1 + 2*vicorr ')
c        stop 'SETGAM-2'
      endif

      zz = iz
      if (ihole .le. 16)  then
         do 10  i = 1, 8
            gamkp(i) = log10 (gamh(i,ihole))
            zk(i) = zh(i,ihole)
   10    continue
         call terp (zk, gamkp, 8, 1, zz, gamach)
      else
c     include data from the tables later.
c     Now gamach=0.1eV for any O-hole for any element.
         gamach = -1.0
      endif

c     Change from log10 (gamma) to gamma
      gamach = 10.0 ** gamach


      return
      end
      subroutine iniptz(ptz,iptz,modus)
        !KJ This routine rewrites the ptz-matrix.

      implicit none
c     which polarization tensor to create      
      integer iptz
c     two ways of working
      integer modus
c     the polarization tensor
      complex*16 ptz(-1:1,-1:1)

      complex*16 zero,one,coni
      parameter (zero=(0,0),one=(1,0),coni=(0,1))
      integer i,j
      real*8 unity(3,3)



      if (iptz.lt.1.or.iptz.gt.10) then
          call wlog('Inieln sees weird iptz - returns without 
     1      changing ptz - danger of calculating nonsense !!')
      endif


      do i=1,3
      do j=1,3
	unity(i,j)=dble(0)
      enddo
      unity(i,i)=dble(1)/dble(3)
      enddo
      do i=-1,1
      do j=-1,1
        ptz(i,j)=zero
      enddo
      enddo
      

      if (modus.eq.1) then
! work in spherical coordinates

         if(iptz.eq.10) then
        do i=-1,1
	do j=-1,1
	   ptz(i,j)=unity(i+2,j+2)
        enddo
	enddo
         else
            i=(iptz-1)/3+1  ! row index
            j=iptz-3*(i-1)  ! column index
            i=i-2 !shift from 1..3 to -1..1
            j=j-2
            ptz(i,j)=one
         endif  


      elseif(modus.eq.2) then
! work in carthesian coordinates      

      if (iptz.eq.10) then ! orientation averaged spectrum
        do i=-1,1
	do j=-1,1
	   ptz(i,j)=unity(i+2,j+2)
        enddo
	enddo
	
        elseif (iptz.eq.1) then   ! x x*
          ptz(1,1)=one/dble(2)
        ptz(-1,-1)=one/dble(2)
        ptz(-1,1)=-one/dble(2)
        ptz(1,-1)=-one/dble(2)
        elseif (iptz.eq.5) then ! y y*
          ptz(1,1)=one/dble(2)
        ptz(-1,-1)=one/dble(2)
        ptz(-1,1)=one/dble(2)
        ptz(1,-1)=one/dble(2)
         elseif (iptz.eq.9) then ! z z*
        ptz(0,0)=one
        elseif (iptz.eq.2) then ! x y*
          ptz(1,1)=one*coni/dble(2)
        ptz(-1,-1)=-one*coni/dble(2)
        ptz(-1,1)=-one*coni/dble(2)
        ptz(1,-1)=one*coni/dble(2)
        elseif (iptz.eq.4) then ! x* y
          ptz(1,1)=-one*coni/dble(2)
        ptz(-1,-1)=one*coni/dble(2)
        ptz(-1,1)=-one*coni/dble(2)
        ptz(1,-1)=one*coni/dble(2)
        elseif (iptz.eq.3) then ! x z*
          ptz(-1,0)=one/dsqrt(dble(2))
        ptz(1,0)=-one/dsqrt(dble(2))
        elseif (iptz.eq.7) then ! x* z
          ptz(0,-1)=one/dsqrt(dble(2))
        ptz(0,1)=-one/dsqrt(dble(2))
        elseif (iptz.eq.6) then ! y z*
          ptz(-1,0)=-one*coni/dsqrt(dble(2))
        ptz(1,0)=-one*coni/dsqrt(dble(2))
        elseif (iptz.eq.8) then ! y* z
          ptz(0,-1)=one*coni/dsqrt(dble(2))
        ptz(0,1)=one*coni/dsqrt(dble(2))
      endif
      
      
        else
          stop 'alien modus in inieln'
        endif


      return
      end


C From HDK@psuvm.psu.edu Thu Dec  8 15:27:16 MST 1994
C 
C The following was converted from Algol recursive to Fortran iterative
C by a colleague at Penn State (a long time ago - Fortran 66, please
C excuse the GoTo's). The following code also corrects a bug in the
C Quicksort algorithm published in the ACM (see Algorithm 402, CACM,
C Sept. 1970, pp 563-567; also you younger folks who weren't born at
C that time might find interesting the history of the Quicksort
C algorithm beginning with the original published in CACM, July 1961,
C pp 321-322, Algorithm 64). Note that the following algorithm sorts
C integer data; actual data is not moved but sort is affected by sorting
C a companion index array (see leading comments). The data type being
C sorted can be changed by changing one line; see comments after
C declarations and subsequent one regarding comparisons(Fortran
C 77 takes care of character comparisons of course, so that comment
C is merely historical from the days when we had to write character
C compare subprograms, usually in assembler language for a specific
C mainframe platform at that time). But the following algorithm is
C good, still one of the best available.


      SUBROUTINE QSORTI (ORD,N,A)
C
C==============SORTS THE ARRAY A(I),I=1,2,...,N BY PUTTING THE
C   ASCENDING ORDER VECTOR IN ORD.  THAT IS ASCENDING ORDERED A
C   IS A(ORD(I)),I=1,2,...,N; DESCENDING ORDER A IS A(ORD(N-I+1)),
C   I=1,2,...,N .  THIS SORT RUNS IN TIME PROPORTIONAL TO N LOG N .
C
C
C     ACM QUICKSORT - ALGORITHM #402 - IMPLEMENTED IN FORTRAN 66 BY
C                                 WILLIAM H. VERITY, WHV@PSUVM.PSU.EDU
C                                 CENTER FOR ACADEMIC COMPUTING
C                                 THE PENNSYLVANIA STATE UNIVERSITY
C                                 UNIVERSITY PARK, PA.  16802
C
      IMPLICIT INTEGER (A-Z)
C
      DIMENSION ORD(N),POPLST(2,20)
      DOUBLE PRECISION X,XX,Z,ZZ,Y
C
C     TO SORT DIFFERENT INPUT TYPES, CHANGE THE FOLLOWING
C     SPECIFICATION STATEMENTS; FOR EXAMPLE, FOR FORTRAN CHARACTER
C     USE THE FOLLOWING:  CHARACTER *(*) A(N)
C
      DOUBLE PRECISION A(N)
C
      NDEEP=0
      U1=N
      L1=1
      DO 1  I=1,N
    1 ORD(I)=I
    2 IF (U1.LE.L1) THEN
         RETURN
      END IF
C
    3 L=L1
      U=U1
C
C PART
C
    4 P=L
      Q=U
C     FOR CHARACTER SORTS, THE FOLLOWING 3 STATEMENTS WOULD BECOME
C     X = ORD(P)
C     Z = ORD(Q)
C     IF (A(X) .LE. A(Z)) GO TO 2
C
C     WHERE "CLE" IS A LOGICAL FUNCTION WHICH RETURNS "TRUE" IF THE
C     FIRST ARGUMENT IS LESS THAN OR EQUAL TO THE SECOND, BASED ON "LEN"
C     CHARACTERS.
C
      X=A(ORD(P))
      Z=A(ORD(Q))
      IF (X.LE.Z) GO TO 5
      Y=X
      X=Z
      Z=Y
      YP=ORD(P)
      ORD(P)=ORD(Q)
      ORD(Q)=YP
    5 IF (U-L.LE.1) GO TO 15
      XX=X
      IX=P
      ZZ=Z
      IZ=Q
C
C LEFT
C
    6 P=P+1
      IF (P.GE.Q) GO TO 7
      X=A(ORD(P))
      IF (X.GE.XX) GO TO 8
      GO TO 6
    7 P=Q-1
      GO TO 13
C
C RIGHT
C
    8 Q=Q-1
      IF (Q.LE.P) GO TO 9
      Z=A(ORD(Q))
      IF (Z.LE.ZZ) GO TO 10
      GO TO 8
    9 Q=P
      P=P-1
      Z=X
      X=A(ORD(P))
C
C DIST
C
   10 IF (X.LE.Z) GO TO 11
      Y=X
      X=Z
      Z=Y
      IP=ORD(P)
      ORD(P)=ORD(Q)
      ORD(Q)=IP
   11 IF (X.LE.XX) GO TO 12
      XX=X
      IX=P
   12 IF (Z.GE.ZZ) GO TO 6
      ZZ=Z
      IZ=Q
      GO TO 6
C
C OUT
C
   13 CONTINUE
      IF (.NOT.(P.NE.IX.AND.X.NE.XX)) GO TO 14
      IP=ORD(P)
      ORD(P)=ORD(IX)
      ORD(IX)=IP
   14 CONTINUE
      IF (.NOT.(Q.NE.IZ.AND.Z.NE.ZZ)) GO TO 15
      IQ=ORD(Q)
      ORD(Q)=ORD(IZ)
      ORD(IZ)=IQ
   15 CONTINUE
      IF (U-Q.LE.P-L) GO TO 16
      L1=L
      U1=P-1
      L=Q+1
      GO TO 17
   16 U1=U
      L1=Q+1
      U=P-1
   17 CONTINUE
      IF (U1.LE.L1) GO TO 18
C
C START RECURSIVE CALL
C
      NDEEP=NDEEP+1
      POPLST(1,NDEEP)=U
      POPLST(2,NDEEP)=L
      GO TO 3
   18 IF (U.GT.L) GO TO 4
C
C POP BACK UP IN THE RECURSION LIST
C
      IF (NDEEP.EQ.0) GO TO 2
      U=POPLST(1,NDEEP)
      L=POPLST(2,NDEEP)
      NDEEP=NDEEP-1
      GO TO 18
C
C END SORT
C END QSORT
C
      END
c///////////////////////////////////////////////////////////////////////
c Distribution:  FEFF_MATH 1.0
c Copyright (c) [2002] University of Washington
c 
c This software was prepared in part with US Government Funding under
c DOE contract DE-FG03-97ER45623.

c Redistribution and use of this Distribution in source and binary
c formats, with or without modification is permitted, provided the 
c following conditions are met:
c 
c Redistributions must retain the above notices and the following list
c of conditions and disclaimer;
c 
c Modified formats carry the marking
c     "Based on or developed using Distribution: FEFF_MATH 1.0
c      FEFF_MATH 1.0 Copyright (c) [2002] University of Washington"
c 
c Recipient acknowledges the right of the University of Washington to
c prepare uses of this Distribution and its modifications that may be
c substantially similar or functionally equivalent to
c Recipient-prepared modifications.
c
c Recipient and anyone obtaining access to the Distribution through
c recipient's actions accept all risk associated with possession and
c use of the Distribution.
c
c THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESS OR IMPLIED
c WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
c MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
c IN NO EVENT SHALL THE UNIVERSITY OF WASHINGTON OR CONTRIBUTORS TO THE
c DISTRIBUTION BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
c EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
c PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
c REVENUE; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
c LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
c NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
c SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
c///////////////////////////////////////////////////////////////////////
c License is applicable for routines below, until otherwise specified.
c
      subroutine bcoef(kinit, ipol, ptz, le2, ltrace, ispin, angks,
     1                 kind, lind, bmat)
c     written by alexei ankudinov; march 2000
c     calculate bmat: the energy independent sum over polarization and
c     angular momenta indices
c     bmat = \sum_{p,p', all m_j} <LS|J><J|R|J1><J1|\alpha_p exp(i kz)|I>
c                    ptz(p,p') 
c            <I|\alpha_p'^* exp(-i kz) J2><J2|R'|J'><J'|L'S'>
c     where R is rotation from spin vector to x-ray k-vector
c     and R' is rotation back
c     see Eq.10 and 11 in Ankudinov,Rehr, Phys.Rev.B (accepted),
c     Theory of solid state contribution to the x-ray elastic scattering
c     aditional rotation matrices are needed when x-ray k-vector
c     is not along the spin-axis (see rotations in rdinp)

c     more precisely it is
c     bmat(l1 l1' j l ml ms; l2 l2' j' l' ml' ms') =
c        (-)**(j-j'+l2'+1) i**(l'-l) \sum_{p,p',mi,m1,mj,m2,mj'}
c        <LS|J>   r^j_{m1,mj}(angks)   3j( j l1 i ; -m1 p mi)
c        (-p)**(l1+l1'+1) ptz(p,p') (-p')**(l2+l2'+1) 
c        3j( j' l2 i ; -m2  p' mi)   r^j'_{m2,mj'}(angks)   <J'|L'S'>
c     where l1 l1' are set by the multipole moment(E1-l1=1,l1'=0;
c     E2-l1=2,l1'=1; M1-l1=1,l1'=1; etc.;
c     j and l define quantum number kappa and for each multipole moment
c     Only few final kappa's are allowed and  it is convinient
c     to denote (l1 l1' j l) by one index 'k'
c     thus  k=1-8 to include both E1 and E2 transitions;
c     ml and ms are projections of orbital and spin moments.

c     bmat  is used to calculate absorption fine structure (chi) via
c       chi = \sum_{k ms,k' ms'}  rkk(k,ms)  rkk(k',ms')
c       \sum_{ml,ml'}  bmat(k ml ms; k' ml' ms')  G_(l' ml' ms';l ml ms)
c     where sum over spins can be moved from first sum to second for
c     spin independent systems. The above expression is suitable for FMS
c     and for MS expansion on can use Eq.15 in RA paper to obtain
c     expression for the termination   matrix
c     T_{lam1 ms,lamN ms'} = \sum_{k k'} rkk(k,ms) rkk(k',ms')
c       \sum_{ml,ml'}  bmat(k ml ms; k' ml' ms') gam(l,lam1,rho1,ms)
c        gamtl(l',lamN,rhoN,ms')
c     Notice that for spin-dependent systems the scattering F matrices
c     in RA paper also should have additional spin indices. In genfmt
c     we currently neglect spin-flip processes which simplifies
c     calculations with MS expansion. (T and F are diagonal in ms,ms')
       
c     This subroutine is written for general spin-dependent asymmetric
c     system and arbitrary polarization tenzor. The symmetry of the 
c     system and polarization tenzor can be used
c     to speed up FMS or MS calculations in appropriate subroutines.
c     (see comments in subroutines mpprmp, fmstot)

c     input:
c       kinit - kappa for initial orbital
c       ipol - polarization type measurement
c       ptz  - polarization tensor (needed only for ipol=1 case)
c       le2  - sets which multipole moments to include (see mkptz)
c       ltrace- .true. for xsect.f, where need to perform trace over ml
c       angks - angle between k-vector and spin-vector 

c     output
c       lind  - orb.mom.(kappa)  needed in fmstot only (for indexing)
c       bmat  - energy independent matrix to calculate absorption 
c       in many cases bmat is diagonal due to the choice of xyz frame,
c       but for general case full 16*(2*lx+1)*16*(2*lx+1) matrix is kept

      implicit double precision (a-h,o-z)
c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}
      complex*16 coni
      parameter (coni = (0,1))

c     need only parameter lx to set max orb momentum
      complex*16 ptz, bmat, pmat, tmat
      dimension ptz(-1:1,-1:1),  bmat(-lx:lx,0:1,8, -lx:lx,0:1,8)
c       to include all possible dipole and quadrupole transitions 
c       final kp, and kpp have 8 possibilities
      logical ltrace

c     local staff
      dimension  t3j( 8, 0:1, -lx:lx+1), x3j(8, -1:1, -lx:lx+1)
c     qmat = <J2|R'|J'><J'|L'S'> - diagonal in kappa index
      dimension qmat( -lx:lx+1, -lx:lx, 0:1, 8)
c     pmat = <J1|\alpha_j exp(i kz)|I> ptz <I|\alpha_k^* exp(-i kz)|J2>
      dimension pmat( -lx:lx+1, 8, -lx:lx+1, 8)
c     tmat = pmat*qmat ; bmat = qmat^T*tmat
      dimension tmat( -lx:lx+1, 8, -lx:lx, 0:1, 8)
c     total and orbital momenta for 8 possible final kappa
      dimension jind(8), lind(8), kind(8)

      external cwig3j

      do 10 i6 = 1, 8
      do 10 i5 = 0 ,1
      do 10 i4 = -lx,lx
      do 10 i3 = 1, 8
      do 10 i2 = 0 ,1
      do 10 i1 = -lx,lx
         bmat( i1, i2, i3, i4, i5, i6) = 0
  10  continue

c     3 dipole transitions
      do 20 k=-1,1
         kap=kinit+k
         if (k.eq.0) kap=-kap
         jkap = abs(kap)
         lkap = kap
         if (kap.le.0) lkap = abs(kap) -1
c        check that orbital momentum does not exceed max allowed
         if (lkap .gt. lx) then
c          set final j and l to unphysical values
           jkap = 0
           lkap = -1 
           kap = 0
         endif
         jind(k+2) = jkap
         lind(k+2) = lkap
         kind(k+2) = kap
  20  continue

c     include 5 quadrupole or 3 mag.dipole  transitions
      do 120 k=-2,2
         jkap = abs(kinit) + k
         if (jkap.le.0) jkap = 0
         kap= jkap
         if (kinit.lt.0 .and. abs(k).ne.1) kap=-jkap
         if (kinit.gt.0 .and. abs(k).eq.1) kap=-jkap
         lkap = kap
         if(kap.le.0) lkap = - kap - 1
         if (lkap.gt.lx .or. le2.eq.0
     1                  .or. (le2.eq.1 .and. abs(k).eq.2)) then
c           set unphysical jkap and lkap to make shorter calculations
            jkap = 0
            lkap = -1
            kap = 0
         endif
         jind(k+6) = jkap
         lind(k+6) = lkap
         kind(k+6) = kap
 120  continue

      if (ipol.eq.0) then
c       polarization average case; bmat is diagonal and simple
        do 100 k = 1, 8
        do 100 ms = 0 ,1
        do 100 ml = -lind(k), lind(k)
c         i2 = (2*l1+1) , where l1 is defined by multipole moment
          i2 = 3
          if (le2.eq.2 .and. k.gt.3) i2 = 5
          bmat(ml,ms,k, ml,ms,k) = 0.5d0 / (2*lind(k)+1.d0) / i2
          if (k.le.3) bmat(ml,ms,k, ml,ms,k) = - bmat(ml,ms,k, ml,ms,k)
 100    continue
      else
c       more complicated bmat for linear(ipol=1) and circular(ipol=2)
c       polarizations
c       Put 3j factors in x3j and t3j. t3j are multiplied by
c       sqrt(2*j'+1) for  further convinience.
        do 30  mp=-lx,lx+1
        do 30  ms=0,1
        do 30  k1=1,8
  30    t3j(k1,ms,mp) = 0.0d0
        do 40  mp=-lx,lx+1
        do 40  ms=-1,1
        do 40  k1=1,8
  40      x3j(k1,ms,mp) = 0.0d0

        do 70  k1 = 1,8
        do 70  mp = -jind(k1)+1,jind(k1)
          do 50 ms=0,1
            j1 = 2 * lind(k1)
            j2 = 1
            j3 = 2 * jind(k1) - 1
            m1 = 2*(mp-ms)
            m2 = 2*ms - 1
            t3j(k1,ms,mp)=sqrt(j3+1.0d0) * cwig3j(j1,j2,j3,m1,m2,2)
            if (mod( (j2-j1-m1-m2)/2 , 2) .ne.0) 
     1          t3j(k1,ms,mp) = - t3j(k1,ms,mp)
c           t3j(m0,i)    are Clebsch-Gordon coefficients
  50      continue
          do 60 i=-1,1
            j1 = 2 * jind(k1) - 1
            j2 = 2
            if (k1.gt.3 .and. le2.eq.2) j2 = 4
            j3 = 2 * abs(kinit) - 1
            m1 = -2*mp + 1
            m2 = 2*i
            x3j(k1,i,mp)= cwig3j(j1,j2,j3,m1,m2,2)
  60      continue
  70    continue

c       calculate qmat
        do 220 i=1,8
        do 220 ms=0,1
        do 220 ml= -lind(i), lind(i)
        do 220 mj= -jind(i)+1, jind(i)
          mp = ml+ms
          jj = 2*jind(i) - 1
          mmj = 2*mj - 1
          mmp = 2*mp - 1
          value = rotwig(angks, jj, mmj, mmp, 2)
          qmat(mj,ml,ms,i) = value * t3j(i,ms,mp)
 220    continue

c       calculate pmat
        do 240 i2 = 1,8
        do 240 m2 = -jind(i2)+1, jind(i2)
        do 240 i1 = 1,8
        do 240 m1 = -jind(i1)+1, jind(i1)
          pmat(m1,i1,m2,i2) = 0
          if (abs(m2-m1).le.2) then
            do 230 j=-1,1
            do 230 i=-1,1
c             check that initial moment is the same
              if (m1-i.eq.m2-j) then
                is = 1
c               (-p) factors for M1 transitions
                if (le2.eq.1 .and. i.gt.0 .and. i1.gt.3) is = -is
                if (le2.eq.1 .and. j.gt.0 .and. i2.gt.3) is = -is
                pmat(m1,i1,m2,i2) = pmat(m1,i1,m2,i2) +
     1          is * x3j(i1,i,m1) * ptz(i,j) * x3j(i2,j,m2)
              endif
 230        continue
c           multiply by (-)^(j-j'+l2'+1) i**(l'-l) factor
c           additional (-) is from Eq.10 (-2*ck)
            is = 1
            if (mod(jind(i1)-jind(i2), 2) .ne.0) is = -is
            if (i2.le.3) is = -is
            pmat(m1,i1,m2,i2) = pmat(m1,i1,m2,i2) * is
     1           * coni**(lind(i2)-lind(i1))
          endif
 240    continue

c       calculate tmat = pmat*qmat
        do 270 i1=1,8
        do 270 ms=0,1
        do 270 ml=-lind(i1), lind(i1)
        do 270 i2=1,8
        do 270 mj=-jind(i2)+1, jind(i2)
          tmat(mj,i2, ml,ms,i1) = 0
          do 260 mp = -jind(i1)+1, jind(i1)
            tmat(mj,i2, ml,ms,i1) = tmat(mj,i2, ml,ms,i1)+
     1           pmat(mj,i2,mp,i1) * qmat(mp,ml,ms,i1)
 260      continue
 270    continue
         
c       calculate bmat = qmat^T * tmat
        do 300 i1=1,8
        do 300 ms1=0,1
        do 300 ml1=-lind(i1), lind(i1)
        do 300 i2=1,8
        do 300 ms2=0,1
        do 300 ml2=-lind(i2), lind(i2)
          bmat(ml2,ms2,i2, ml1,ms1,i1) = 0
          do 280 mj=-jind(i2)+1, jind(i2)
            bmat(ml2,ms2,i2, ml1,ms1,i1) = bmat(ml2,ms2,i2, ml1,ms1,i1)+
     1      qmat(mj,ml2,ms2,i2) * tmat(mj,i2,ml1,ms1,i1) 
 280      continue
 300    continue
c       end of ipol=1,2 cases
      endif 

      if (ltrace) then
c       need to trace bmat over ml for xsect.f
        do 390 i1 = 1, 8
        do 390 ms1 = 0,1
        do 390 i2 = 1, 8
        do 390 ms2 = 0,1
          if (lind(i1).ne.lind(i2) .or. ms1.ne.ms2) then
               bmat(0,ms2,i2, 0,ms1,i1) = 0
          else
             do 360 ml = 1, lind(i1)
               bmat(0,ms1,i2, 0,ms1,i1) =  bmat(0,ms1,i2, 0,ms1,i1) +
     1         bmat(-ml,ms1,i2, -ml,ms1,i1) + bmat(ml,ms1,i2, ml,ms1,i1)
 360         continue
          endif
 390    continue
      endif

      if (ispin .eq. 0) then
c       G(Ls,L's') is spin diagonal; trace over spin
        do 480 i1 = 1, 8
        do 480 i2 = 1, 8
        do 480 ml1 = -lind(i1), lind(i1)
        do 480 ml2 = -lind(i2), lind(i2)
           bmat(ml2,0,i2, ml1,0,i1) =   bmat(ml2,0,i2, ml1,0,i1) +
     1                                  bmat(ml2,1,i2, ml1,1,i1)
 480    continue
      elseif (ispin.eq.2 .or. (ispin.eq.1 .and. nspx.eq.1)) then
c       move spin up part into the position of spin-down
        do 490 i1 = 1, 8
        do 490 i2 = 1, 8
        do 490 ml1 = -lind(i1), lind(i1)
        do 490 ml2 = -lind(i2), lind(i2)
           bmat(ml2,0,i2, ml1,0,i1) =   bmat(ml2,1,i2, ml1,1,i1)
 490    continue

      endif

      return
      end
      subroutine besjn (x, jl, nl)

c-----------------------------------------------------------------------
c
c     purpose:  to calculate the spherical bessel functions jl and nl
c               for l = 0 to 30 (no offset)
c
c     arguments:
c       x = argument of jl and nl
c       jl = jl bessel function (abramowitz conventions)
c       nl = nl bessel function (abramowitz yl conventions)
c            Note that this array nl = abramowitz yl.
c       jl and nl must be dimensioned 
c            complex*16 jl(ltot+2), nl(ltot+2), with ltot defined in 
c            dim.h.
c
c     notes:  jl and nl should be calculated at least to 10 place
c             accuracy for the range 0<x<100 according to spot
c             checks with tables
c
c     error messages written with PRINT statement.
c
c     first coded by r. c. albers on 14 dec 82
c
c     version 3
c
c     last modified: 27 jan 83 by r. c. albers
c     dimension of jl,nl changed from 31 to 26  (10 aug 89) j. rehr
c     modified again, siz, June 1992
c
c-----------------------------------------------------------------------

      implicit double precision (a-h, o-z)
c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}

      complex*16 x
      complex*16 jl(ltot+2), nl(ltot+2)
      complex*16 cjl(ltot+2), sjl(ltot+2), cnl(ltot+2), snl(ltot+2)

      complex*16 xjl,xnl,asx,acx
      complex*16 xi,xi2,xi3,xi4,xi5,xi6,xi7,xi8,xi9,xi10,xi11

      parameter (xcut = 1.d0, xcut1 = 7.51d0, xcut2 = 5.01d0)

      if (dble(x) .le. 0)  stop 'Re(x) is .le. zero in besjn'

      lmaxp1 = ltot+2

      if (dble(x) .lt. xcut .and. abs(dimag(x)) .lt. xcut)  then
c        case Re(x) < 1, just use series expansion
         do 10 il = 1,lmaxp1
            l = il-1
            ifl = 0
            call bjnser (x,l,xjl,xnl,ifl)
            jl(il) = xjl
            nl(il) = xnl
   10    continue

      elseif (dble(x) .lt. xcut1 .and. abs(dimag(x)) .lt. xcut1)  then

c        case 1 <= Re(x) < 7.5

         call bjnser (x,lmaxp1-1,xjl,xnl,1)
         jl(lmaxp1) = xjl

         call bjnser (x,lmaxp1-2,xjl,xnl,1)
         jl(lmaxp1-1) = xjl

         if (dble(x) .lt. xcut2 .and. abs(dimag(x)) .lt. xcut2)  then
c           Re(x) < 5
            call bjnser (x,0,xjl,xnl,2)
            nl(1) = xnl
            call bjnser (x,1,xjl,xnl,2)
            nl(2) = xnl
         else
c           Re(x) >= 5
            asx = sin(x)
            acx = cos(x)
            xi = 1 / x
            xi2 = xi**2
            nl(1) = -acx*xi
            nl(2) = -acx*xi2 - asx*xi
         endif

c        Use recursion relation 10.1.19 to get nl and jl
         do 50 lp1 = 3, lmaxp1
            l = lp1 - 2
            tlxp1 = 2*l + 1
            nl(lp1) = tlxp1 * nl(lp1-1) / x  -  nl(lp1-2)
   50    continue

         do 60 lxx = 3,lmaxp1
            lp1 = lmaxp1+1-lxx
            l = lp1-1
            tlxp3 = 2*l + 3
            jl(lp1) = tlxp3 * jl(lp1+1) / x  -  jl(lp1+2)
   60    continue

      else
c        case Re(x) > 7.5
c        Use AS 10.1.8 and 10.1.9, sjl=P, qjl=Q, note that AS formulae
c        use cos (z - n*pi/2), etc., so cos and sin terms get a bit
c        scrambled (mod 4) here, since n is integer.  These are hard-
c        coded into the terms below.
         xi = 1 / x
         xi2  = xi*xi
         xi3  = xi*xi2
         xi4  = xi*xi3
         xi5  = xi*xi4
         xi6  = xi*xi5
         xi7  = xi*xi6
         xi8  = xi*xi7
         xi9  = xi*xi8
         xi10 = xi*xi9
         xi11 = xi*xi10

         sjl(1) = xi
         sjl(2) = xi2
         sjl(3) = 3*xi3 - xi
         sjl(4) = 15*xi4 - 6*xi2
         sjl(5) = 105*xi5 - 45*xi3 + xi
         sjl(6) = 945*xi6 - 420*xi4 + 15*xi2
         sjl(7) = 10395*xi7 - 4725*xi5 + 210*xi3 - xi
         sjl(8) = 135135*xi8 - 62370*xi6 + 3150*xi4 - 28*xi2
         sjl(9) = 2027025*xi9 - 945945*xi7 + 51975*xi5 
     1            - 630*xi3 + xi
         sjl(10) = 34459425*xi10 - 16216200*xi8 + 945945*xi6 
     1            - 13860*xi4 + 45*xi2
         sjl(11) = 654729075*xi11 - 310134825*xi9 + 18918900*xi7 
     1            - 315315*xi5 + 1485*xi3 - xi
         cjl(1) = 0
         cjl(2) = -xi
         cjl(3) = -3*xi2
         cjl(4) = -15*xi3 + xi
         cjl(5) = -105*xi4 + 10*xi2
         cjl(6) = -945*xi5 + 105*xi3 - xi
         cjl(7) = -10395*xi6 + 1260*xi4 - 21*xi2
         cjl(8) = -135135*xi7 + 17325*xi5 - 378*xi3 + xi
         cjl(9) = -2027025*xi8 + 270270*xi6 - 6930*xi4 + 36*xi2
         cjl(10) = -34459425*xi9 + 4729725*xi7 - 135135*xi5 
     1             + 990*xi3 - xi
         cjl(11) = -654729075*xi10 + 91891800*xi8 - 2837835*xi6 
     1             + 25740*xi4 - 55*xi2
         do 80 ie = 1,11
            snl(ie) = cjl(ie)
            cnl(ie) = -sjl(ie)
   80    continue
         do 90 lp1 = 12,lmaxp1
            l = lp1-2
            tlxp1 = float(2*l+1)
            sjl(lp1) = tlxp1*xi*sjl(lp1-1)-sjl(lp1-2)
            cjl(lp1) = tlxp1*xi*cjl(lp1-1)-cjl(lp1-2)
            snl(lp1) = tlxp1*xi*snl(lp1-1)-snl(lp1-2)
            cnl(lp1) = tlxp1*xi*cnl(lp1-1)-cnl(lp1-2)
   90    continue
         asx = sin(x)
         acx = cos(x)
         do 110 lp1 = 1,lmaxp1
            jl(lp1) = asx*sjl(lp1)+acx*cjl(lp1)
            nl(lp1) = asx*snl(lp1)+acx*cnl(lp1)
  110    continue
      endif

      return
      end
      subroutine besjh (x, lbmax, jl, hl)

c-----------------------------------------------------------------------
c
c     purpose:  to calculate the spherical bessel functions jl and hl
c               for l = 0 to lbmax (no offset)
c
c     arguments:
c       x = argument of jl and nl
c       lbmax
c       jl = jl bessel function (abramowitz conventions)
c       hl = hl^+ bessel function (messiah conventions) for Im x >=0
c       hl = hl^- bessel function (messiah conventions) for Im x < 0
c       jl and hl must be dimensioned 
c            complex*16 jl(0:lbmax), hl(0:lbmax), 
c
c     notes:  jl and hl should be calculated at least to 10 place
c             accuracy for the range 0<x<100 according to spot
c             checks with tables
c
c     error messages written with PRINT statement.
c
c     first coded by r. c. albers on 14 dec 82
c
c     version 3
c
c     last modified: 27 jan 83 by r. c. albers
c     dimension of jl,nl changed from 31 to 26  (10 aug 89) j. rehr
c     modified again, siz, June 1992
c     rewritten for jl and hl by a.l. ankudinov feb 2000
c
c-----------------------------------------------------------------------

      implicit double precision (a-h, o-z)
c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}

      complex*16 x
      complex*16 jl(0:lbmax), nl(ltot+2)
      complex*16 hl(0:lbmax)
      complex*16 cjl(ltot+2), sjl(ltot+2)

      complex*16 xjl,xnl,asx,acx, epx
      complex*16 xi,xi2,xi3,xi4,xi5,xi6,xi7,xi8,xi9,xi10,xi11

      parameter (xcut = 1.d0, xcut1 = 7.51d0, xcut2 = 5.01d0)
      complex*16 coni
      parameter (coni=(0,1))

      if (dble(x) .lt. 0)  stop 'Re(x) is .lt. zero in besjh'

      lmax = min(lbmax, ltot+1)
      lmaxp1 = lmax + 1

      if (dble(x) .lt. xcut .and. abs(dimag(x)) .lt. xcut)  then
c        case Re(x) < 1, just use series expansion
         do 10 ll = 0,lmax
            ifl = 0
            call bjnser (x,ll,xjl,xnl,ifl)
            jl(ll) = xjl
            hl(ll) = -xnl + coni*xjl
   10    continue

      elseif (dble(x) .lt. xcut1 .and. abs(dimag(x)) .lt. xcut1)  then

c        case 1 <= Re(x) < 7.5

         call bjnser (x,lmax,xjl,xnl,1)
         jl(lmax) = xjl

         call bjnser (x,lmax-1,xjl,xnl,1)
         jl(lmax-1) = xjl

         if (dble(x) .lt. xcut2 .and. abs(dimag(x)) .lt. xcut2)  then
c           Re(x) < 5
            call bjnser (x,0,xjl,xnl,2)
            nl(1) = xnl
            call bjnser (x,1,xjl,xnl,2)
            nl(2) = xnl
         else
c           Re(x) >= 5
            asx = sin(x)
            acx = cos(x)
            xi = 1 / x
            xi2 = xi**2
            nl(1) = -acx*xi
            nl(2) = -acx*xi2 - asx*xi
         endif

c        Use recursion relation 10.1.19 to get nl and jl
         do 50 lp1 = 3, lmaxp1
            l = lp1 - 2
            tlxp1 = 2*l + 1
            nl(lp1) = tlxp1 * nl(lp1-1) / x  -  nl(lp1-2)
   50    continue

         do 60 lxx = 3,lmaxp1
            lp1 = lmaxp1+1-lxx
            l = lp1-1
            tlxp3 = 2*l + 3
            jl(l) = tlxp3 * jl(l+1) / x  -  jl(l+2)
   60    continue

         do 65 il = 1, lmaxp1
            l = il - 1
            hl(l) = -nl(il) + coni*jl(l)
   65    continue

      else
c        case Re(x) > 7.5
c        Use AS 10.1.8 and 10.1.9, sjl=P, qjl=Q, note that AS formulae
c        use cos (z - n*pi/2), etc., so cos and sin terms get a bit
c        scrambled (mod 4) here, since n is integer.  These are hard-
c        coded into the terms below.
         xi = 1 / x
         xi2  = xi*xi
         xi3  = xi*xi2
         xi4  = xi*xi3
         xi5  = xi*xi4
         xi6  = xi*xi5
         xi7  = xi*xi6
         xi8  = xi*xi7
         xi9  = xi*xi8
         xi10 = xi*xi9
         xi11 = xi*xi10

         sjl(1) = xi
         sjl(2) = xi2
         sjl(3) = 3*xi3 - xi
         sjl(4) = 15*xi4 - 6*xi2
         sjl(5) = 105*xi5 - 45*xi3 + xi
         sjl(6) = 945*xi6 - 420*xi4 + 15*xi2
         sjl(7) = 10395*xi7 - 4725*xi5 + 210*xi3 - xi
         sjl(8) = 135135*xi8 - 62370*xi6 + 3150*xi4 - 28*xi2
         sjl(9) = 2027025*xi9 - 945945*xi7 + 51975*xi5 
     1            - 630*xi3 + xi
         sjl(10) = 34459425*xi10 - 16216200*xi8 + 945945*xi6 
     1            - 13860*xi4 + 45*xi2
         sjl(11) = 654729075*xi11 - 310134825*xi9 + 18918900*xi7 
     1            - 315315*xi5 + 1485*xi3 - xi
         cjl(1) = 0
         cjl(2) = -xi
         cjl(3) = -3*xi2
         cjl(4) = -15*xi3 + xi
         cjl(5) = -105*xi4 + 10*xi2
         cjl(6) = -945*xi5 + 105*xi3 - xi
         cjl(7) = -10395*xi6 + 1260*xi4 - 21*xi2
         cjl(8) = -135135*xi7 + 17325*xi5 - 378*xi3 + xi
         cjl(9) = -2027025*xi8 + 270270*xi6 - 6930*xi4 + 36*xi2
         cjl(10) = -34459425*xi9 + 4729725*xi7 - 135135*xi5 
     1             + 990*xi3 - xi
         cjl(11) = -654729075*xi10 + 91891800*xi8 - 2837835*xi6 
     1             + 25740*xi4 - 55*xi2
         do 90 lp1 = 12,lmaxp1
            l = lp1-2
            tlxp1 = float(2*l+1)
            sjl(lp1) = tlxp1*xi*sjl(lp1-1)-sjl(lp1-2)
            cjl(lp1) = tlxp1*xi*cjl(lp1-1)-cjl(lp1-2)
   90    continue
         asx = sin(x)
         acx = cos(x)
         if (dimag(x).ge. 0.d0) then
           epx = exp(coni*x)
         else 
           epx = exp(-coni*x)
         endif
         do 110 ll = 0,lmax
            lp1 = ll + 1
            jl(ll) = asx*sjl(lp1)+acx*cjl(lp1)
            if (dimag(x).ge. 0.d0) then
              hl(ll) = (sjl(lp1)+coni*cjl(lp1)) * epx
            else
              hl(ll) = (sjl(lp1)-coni*cjl(lp1)) * epx
            endif
  110    continue
      endif

      return
      end
      subroutine bjnser (x, l, jl, nl, ifl)

c-----------------------------------------------------------------------
c
c     subroutine: bjnser (x,l,jl,nl,ifl)
c
c     purpose:  to calculate the spherical bessel functions jl and nl
c
c     arguments:
c       x = argument of jl and nl
c       l = l value calculated (no offset)
c       jl = jl bessel function (abramowitz conventions)
c       nl = nl bessel function (abramowitz yl conventions)
c       ifl = 0 return both jl and nl
c             1 return jl only
c             2 return nl only
c
c     notes:  jl and nl are calculated by a series
c             expansion according to 10.1.2 and 10.1.3
c             in abramowitz and stegun (ninth printing),
c             page 437
c
c             error msgs written with PRINT statements.
c
c     first coded by r. c. albers on 26 jan 83
c
c     version 2
c
c     last modified: 27 jan 83 by r. c. albers
c
c-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

      complex*16 x,u,ux,del,pj,pn
      complex*16 jl,nl

      character*512 slog

      parameter (niter = 160, tol = 1.e-15)

      if (l .lt. 0) then
         call wlog(' l .lt. 0 in bjnser')
         stop 'bjnser 1'
      endif
      if (dble(x).lt. 0.) then
         write(slog,30) x
         call wlog(slog)
   30    format (' x = ', 1p, 2e14.6, ' is .le. 0 in bjnser')
         stop 'bjnser 2'
      endif

      lp1 = l+1
      u = x**2 / 2

c     make djl = 1 * 3 * 5 * ... * (2*l+1),
c          dnl = 1 * 3 * 5 * ... * (2*l-1)
      djl = 1
      fac = -1
      do 50 il = 1, lp1
         fac = fac + 2
         djl = fac * djl
   50 continue
      dnl = djl / (2*l+1)


      if (ifl .eq. 2)   goto 90
c     make jl
c     pj is term in { } in 10.1.2, del is last factor in the series
c     convergence test is (last factor)/(total term) <= tol
      pj = 1
      nf = 1
      nfac = 2*l + 3
      den = nfac
      sgn = -1
      ux = u
      do 60 il = 1, niter
         del = sgn*ux / den
         pj = pj + del
         trel = abs (del / pj)
         if (trel .le. tol)  goto 80
         sgn = -sgn
         ux = u*ux
         nf = nf+1
         nfac = nfac+2
         den = nf * nfac * den
   60 continue
      stop  'jl does not converge in bjnser'
   80 jl = pj * (x**l) / djl

   90 if (ifl.eq.1) return
c     make nl
c     pn is term in { } in 10.1.3, del is last factor in the series
c     convergence test is (last factor)/(total term) <= tol
      pn = 1
      nf = 1
      nfac = 1 - 2*l
      den = nfac
      sgn = -1
      ux = u
      do 100  il = 1, niter
         del = sgn * ux / den
         pn = pn + del
         trel = abs (del / pn)
         if (trel .le. tol) goto 120
         sgn = -sgn
         ux = u*ux
         nf = nf+1
         nfac = nfac+2
         den = nf * nfac * den
  100 continue
      stop  'nl does not converge in bjnser'
  120 nl = -pn * dnl / (x**lp1)

      return
      end
      subroutine conv(omega,xsec,ne1,vicorr)
c     multiply xsec by theta(omega-efermi) and
c     convolute xsec(omega) with  xloss/((omega-omega0)**2+xloss**2)/pi
c     the result is xsec0(omega0)

      implicit double precision (a-h, o-z)
c={../HEADERS/const.h
      parameter (pi = 3.14159 26535 89793 23846 26433d0)
      parameter (one = 1, zero = 0)
      parameter (third = one/3)
      parameter (raddeg = 180 / pi)
      complex*16 coni
      parameter (coni = (0,1))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919 158 292 677 512 811d0)

      parameter (bohr = 0.529 177 249d0, ryd  = 13.605 698d0)
      parameter (hart = 2 * ryd)
      parameter (alpinv = 137.035 989 56d0)
c     fine structure alpha
      parameter (alphfs = 1 / alpinv)
c= ../HEADERS/const.h}
c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}
      dimension  omega(nex)
      complex*16 xsec(nex), xsec0(nex), xsecdx

      complex*16 conv1
      external conv1

      do 100 ie = 1,ne1
         xsec0(ie) = 0.0d0
         omega0 = omega(ie)
c        Add one more point to correct for the finite grid
c        at large energies. Use linear interpolation.
         dx = max( omega(ne1) - omega(ne1-1), 50*vicorr)
         xlast = omega(ne1)+dx
         dx = dx / ( omega(ne1) - omega(ne1-1))
         xsecdx = xsec(ne1)+ (xsec(ne1)-xsec(ne1-1)) * dx

c        first interval
         do 50  i = 1, ne1-1
            xsec0(ie) = xsec0(ie) +
     1      conv1(omega(i),omega(i+1),xsec(i),xsec(i+1),omega0,vicorr)
  50     continue
c        last interval
         xsec0(ie) = xsec0(ie) +
     1   conv1(omega(ne1),xlast,xsec(ne1),xsecdx,omega0,vicorr)
         xsec0(ie) = xsec0(ie) /real(pi)
  100 continue
      do 200 ie = 1, ne1
  200 xsec(ie) = xsec0(ie)

      return
      end

      complex*16 function conv1(x1,x2,y1,y2,x0,xloss)
c     convolution of function 1/(omega-omega0-i*xloss)/pi
c     makes linear interpolation for function between x1,x2 and
c     takes advantage that the integral can be taken analytically.
      implicit double precision (a-h, o-z)
      complex*16  y1, y2, t, coni,dum, a, b
      parameter (coni = (0.0,1.0))

      d = (x2-x1) / 2.0
      a = dble(y2-y1) / 2.0
      b = dble(y2+y1) / 2.0
      t = d / ( (x1+x2)/2 - x0 - coni*xloss )
      if (abs(t) .ge. 0.1) then
         dum = 2.0*a + (b - a/t) * log((1+t)/(1-t))
      else
         dum = 2.0*b*(t+t**3 / 3.0) - 2.0/3.0 * a*t**2
      endif
      conv1 = dimag (dum)

      d = (x2-x1) / 2.0
      a = dimag(y2-y1) / 2.0
      b = dimag(y2+y1) / 2.0
      t = d / ( (x1+x2)/2 - x0 - coni*xloss )
      if (abs(t) .ge. 0.1) then
         dum = 2.0*a + (b - a/t) * log((1+t)/(1-t))
      else
         dum = 2.0*b*(t+t**3 / 3.0) - 2.0/3.0 * a*t**2
      endif
      conv1 = conv1 + coni* dimag( dum)

      return
      end
      subroutine cpl0 (x, pl0, lmaxp1)
      implicit double precision (a-h, o-z)

c-----------------------------------------------------------------------
c
c     cpl0:  Calculate associated legendre polynomials p_l0(x)
c            by recursion.
c            Adapted from aslgndr.
c
c     first written: (25 june 86) by j. j. rehr
c
c     version 1 (25 june 86) (aslgndr)
c     version 2 (March, 1992) siz
c
c-----------------------------------------------------------------------

      dimension pl0 (lmaxp1)

      lmax = lmaxp1-1

c     calculate legendre polynomials p_l0(x) up to l=lmax
      pl0(1) = 1
      pl0(2) = x
      do 10  il = 2, lmax
         l = il-1
         pl0(il+1) = ( (2*l+1)*x*pl0(il) - l*pl0(l) ) / il
   10 continue

      return
      end
      subroutine csomm (dr,dp,dq,dpas,da,m,np)
c Modified to use complex p and q.  SIZ 4/91
c integration by the method of simpson of (dp+dq)*dr**m from 
c 0 to r=dr(np)
c dpas=exponential step;
c for r in the neighborhood of zero (dp+dq)=cte*r**da
c **********************************************************************
      implicit double precision (a-h,o-z)
      dimension dr(*)
      complex*16  dp(*),dq(*),da,dc
      mm=m+1
      d1=da+mm
      da=0.0
      db=0.0
      do 70 i=1,np
      dl=dr(i)**mm
      if (i.eq.1.or.i.eq.np) go to 10
      dl=dl+dl
      if ((i-2*(i/2)).eq.0) dl=dl+dl
   10 dc=dp(i)*dl
      da=da+dc
      dc=dq(i)*dl
      da=da+dc
   70 continue
      da=dpas*da/3
      dd=exp(dpas)-1.0
      db=d1*(d1+1.0)*dd*exp((d1-1.0)*dpas)
      db=dr(1)*(dr(2)**m)/db
      dd=(dr(1)**mm)*(1.0+1.0/(dd*(d1+1.0)))/d1
      da=da+dd*(dp(1)+dq(1))-db*(dp(2)+dq(2))
      return
      end
      subroutine csomm2 (dr,dp,dpas,da,rnrm,np)
c Modified to use complex p and q.  SIZ 4/91
c Modified to use double simpson integration ALA 3/97
c integration by the method of simpson of dp*dr from 
c 0 to r=rnrm  with proper end corrections
c dpas=exponential step;
c for r in the neighborhood of zero dp=cte*r**da
c **********************************************************************
      implicit double precision (a-h,o-z)
      dimension dr(*)
      complex*16  dp(*),da,dc

      d1=dble(da)+1
      da=0.0
      db=0.0
c      np-2=inrm -point of grid just below rnrm
      a1=log(rnrm/dr(np-2)) / dpas
      a2=a1**2/8.0d0
      a3=a1**3/12.0d0
      do 70 i=1,np
         if (i.eq.1) then
            dc=dp(i) *dr(i)*9.0d0/24.0d0
         elseif (i.eq.2) then
            dc=dp(i) *dr(i)*28.0d0/24.0d0
         elseif (i.eq.3) then
            dc=dp(i)*dr(i)*23.0d0/24.0d0
         elseif (i.eq.np-3) then
            dc=dp(i)*dr(i)*(25.0d0/24.0d0-a2+a3)
         elseif (i.eq.np-2) then
            dc=dp(i)*dr(i)*(0.5d0+a1-3*a2-a3)
         elseif (i.eq.np-1) then
            dc=dp(i)*dr(i)*(-1.0d0/24.0d0+5*a2-a3)
         elseif (i.eq.np) then
            dc=dp(i)*dr(i)*(-a2+a3)
         else
c           like trapesoidal rule
            dc=dp(i)*dr(i)
         endif
         da=da+dc
   70 continue
      da=dpas*da

c     add initial point (r=0) correction
      dd=exp(dpas)-1.0
      db=d1*(d1+1.0)*dd*exp((d1-1.0)*dpas)
      db=dr(1)/db
      dd=(dr(1))*(1.0+1.0/(dd*(d1+1.0)))/d1
      da=da+dd*dp(1)-db*dp(2)
      return
      end
      double precision function cwig3j (j1,j2,j3,m1,m2,ient)
c     wigner 3j coefficient for integers  (ient=1)
c                         or semiintegers (ient=2)
c     other arguments should be multiplied by ient
 
      implicit double precision (a-h,o-z)
      parameter (idim = 58)
      character*512 slog
c     dimensions  modified for larger arguments by ala 12.12.94
      dimension al(idim+1),m(12)
      save ini, al
      data ini/1/
c     idim-1 is the largest argument of factorial to calculate

      m3=-m1-m2
      if (ini) 1,21,1
c        initialisation of the log's of the factorials
 1    ini=0
      al(1)=0.0d 00
      do 11 i=1,idim
         b=i
 11      al(i+1)=al(i)+ log(b)
 21   cwig3j=0.0d 00
      if (((ient-1)*(ient-2)).ne.0) go to 101
      ii=ient+ient
c        test triangular inequalities, parity and maximum values of m
      if (( abs(m1)+ abs(m2)).eq.0.and.mod(j1+j2+j3,ii).ne.0) go to 99
      m(1)=j1+j2-j3
      m(2)=j2+j3-j1
      m(3)=j3+j1-j2
      m(4)=j1+m1
      m(5)=j1-m1
      m(6)=j2+m2
      m(7)=j2-m2
      m(8)=j3+m3
      m(9)=j3-m3
      m(10)=j1+j2+j3+ient
      m(11)=j2-j3-m1
      m(12)=j1-j3+m2
      do 41 i=1,12
         if (i.gt.10) go to 31
         if (m(i).lt.0) go to 99
 31      if (mod(m(i),ient).ne.0) go to 101
         m(i)=m(i)/ient
         if (m(i).gt.idim) go to 101
 41   continue

c        calculate 3j coefficient
      max0= max(m(11),m(12),0)+1
      min0= min(m(1),m(5),m(6))+1
      isig=1
      if (mod(max0-1,2).ne.0) isig=-isig
      c=-al(m(10)+1)
      do 61 i=1,9
 61   c=c+al(m(i)+1)
      c=c/2.0d 00
      do 71 i=max0,min0
      j=2-i
      b=al(i)+al(j+m(1))+al(j+m(5))+al(j+m(6))+al(i-m(11))+al(i-m(12))
      cwig3j=cwig3j+isig* exp(c-b)
 71   isig=-isig
      if (mod(j1-j2-m3,ii).ne.0) cwig3j=-cwig3j
 99   return
 101     write(slog,'(a,6i5)') 'error in cwig3j ',j1,j2,j3,m1,m2,ient
         call wlog(slog)
      stop
      end
      double precision function determ(array,nord,nrows)
c
c     calculate determinate of a square matrix
c        (from bevington "data reduction and error analysis
c         for the physical sciences" pg 294)
c     array: matrix to be analyzed
c     nord: order of matrix
c     nrows:  first dimension of matrix in calling routine
c
      double precision array(nrows,nrows)
      determ = 1.
      do 150 k=1,nord
c
c
        if (array(k,k).ne.0) go to 130
        do 100 j=k,nord
          if (array(k,j).ne.0) go to 110
  100   continue
        determ = 0.
        go to 160
c
  110   do 120 i=k,nord
          saved = array(i,j)
          array(i,j) = array(i,k)
  120   array(i,k) = saved
        determ = -determ
c
  130   determ = determ*array(k,k)
        if (k.ge.nord) go to 150
        k1 = k+1
        do 140 i=k1,nord
          do 140 j=k1,nord
  140   array(i,j) = array(i,j)-array(i,k)*array(k,j)/array(k,k)
  150 continue
  160 return
c end double precision function determ
      end
      double precision function dist (r0, r1)
c     find distance between cartesian points r0 and r1
      implicit double precision (a-h, o-z)
      dimension r0(3), r1(3)
      dist = 0
      do 10  i = 1, 3
         dist = dist + (r0(i) - r1(i))**2
   10 continue
      dist = sqrt (dist)
      return
      end
      double precision function rotwig (beta, jj, m1, m2, ient)
c     uses Wigner formula (Messiah eq.C.72) to calculate rotation matrix
c     for integers  (ient=1)  or semiintegers (ient=2)
c     other arguments (except beta) should be multiplied by ient
 
      implicit double precision (a-h,o-z)
      parameter (idim = 58)
c     dimensions  modified for larger arguments by ala 12.12.94
      dimension al(idim+1),m(12)
      save ini, al
      data ini/1/
c     idim-1 is the largest argument of factorial to calculate

      if (((ient-1)*(ient-2)).ne.0) stop ' Illegal ient in rotwig.'

      if (ini.eq.1) then
c       initialisation of the log's of the factorials
        ini=0
        al(1)=0.0d 00
        do 11 i=1,idim
           b=i
 11        al(i+1)=al(i)+ log(b)
      endif
      rotwig = 0.d0

      if ( m1.ge.0 .and. abs(m1).ge.abs(m2)) then
         m1p = m1 
         m2p = m2
         betap = beta
         isign = 1
      elseif (m2.ge.0 .and. abs(m2).ge.abs(m1)) then
         m1p = m2
         m2p = m1
         betap = - beta
         isign = 1
      elseif (m1.le.0 .and. abs(m1).ge.abs(m2)) then
         m1p = - m1
         m2p = - m2
         betap = beta
         isign = (-1)**( (m1-m2)/ient ) 
      else
         m1p = - m2
         m2p = - m1
         betap = - beta
         isign = (-1)**( (m2-m1)/ient ) 
      endif

      temp = 0.d0
      zeta = cos ( betap / 2.d0 )
      eta  = sin ( betap / 2.d0 )
      do 100 it = m1p - m2p, jj - m2p, ient
        m(1) = 1 + (jj+m1p) / ient
        m(2) = 1 + (jj-m1p) / ient
        m(3) = 1 + (jj+m2p) / ient
        m(4) = 1 + (jj-m2p) / ient
        m(5) = 1 + (jj+m1p-it) / ient
        m(6) = 1 + (jj-m2p-it) / ient
        m(7) = 1 + it / ient
        m(8) = 1 + (m2p-m1p+it) / ient
        m(9)  = (2*jj+m1p-m2p-2*it) / ient 
        m(10) = (2*it-m1p+m2p) / ient 
        factor = 0.d0
        do 110 i = 1,4
  110     factor = factor + al(m(i))/2.d0 - al(m(i+4))
c       special cases to resolve 0.d0**0 problem (a.ankudinov, may 2001)
        if (m(10).eq.0 .and. m(9).eq.0) then
          temp = temp + (-1)**(it/ient)*exp(factor)
        elseif (m(10).eq.0) then
          temp = temp + (-1)**(it/ient)*zeta**m(9)*exp(factor)
        elseif (m(9).eq.0) then
          temp = temp + (-1)**(it/ient)*eta**m(10)*exp(factor)
        else
c         general expression
          temp = temp+ (-1)**(it/ient)*zeta**m(9)*eta**m(10)*exp(factor)
        endif
  100 continue

      rotwig = isign * temp
     
      return
      end
      subroutine phamp (rmt, pu, qu, ck, jl, nl, jlp, nlp, ikap,
     1                  ph, amp)
c     calculate phase shift at mt radius
c     needs to calculate atan of complex variable (coded below)
      implicit double precision (a-h, o-z)
c={../HEADERS/const.h
      parameter (pi = 3.14159 26535 89793 23846 26433d0)
      parameter (one = 1, zero = 0)
      parameter (third = one/3)
      parameter (raddeg = 180 / pi)
      complex*16 coni
      parameter (coni = (0,1))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919 158 292 677 512 811d0)

      parameter (bohr = 0.529 177 249d0, ryd  = 13.605 698d0)
      parameter (hart = 2 * ryd)
      parameter (alpinv = 137.035 989 56d0)
c     fine structure alpha
      parameter (alphfs = 1 / alpinv)
c= ../HEADERS/const.h}
c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}
      external besjn, atan2c

      complex*16 pu, qu, ck,  jl, nl, jlp, nlp, ph, amp
      complex*16 xkr, a, b, factor

c     initialize staff
      xkr = ck*rmt
      isign=1
      if (ikap.lt.0) isign = -1
      a = ck*alphfs
      factor = isign*a/(1+sqrt(1+a**2))

c     find a and b that pu = rmt*(a*jl+b*nl), qu=factor*rmt*(a*jlp+b*nlp)
      a = isign*ck*xkr* (pu*nlp - qu*nl/factor)
      b = isign*ck*xkr* (qu*jl/factor - pu*jlp)

c     pu =  amp * rmt * (jl*cos(ph) - nl*sin(ph))
c     qu =  amp * rmt * (jlp*cos(ph) - nlp*sin(ph)) * factor
c     tan(ph) = - b/a
      b = -b
      call atan2c ( a, b, amp, ph)

      return
      end
      subroutine atancc(temp, phx)
c     phx=atan(temp), for complex numbers
      implicit double precision (a-h, o-z)
      complex*16 temp, phx

      xx = dble (temp)
      yy = dimag(temp)
      if (xx .ne. 0)  then
         alph = (1 - xx**2 - yy**2)
         alph = sqrt(alph**2 + 4*xx**2) - alph
         alph = alph / (2 * xx)
         alph = atan (alph)
      else
         alph = 0
      endif
      beta = (xx**2 + (yy+1)**2) / (xx**2 + (yy-1)**2)
      beta = log(beta) / 4
      phx = dcmplx (alph, beta)

      return
      end

      subroutine atan2c(a, b, ampl, phx)
c     for complex a, b find complex ampl, phx such that:
c     a= ampl*cos(phx)  and  b= ampl*sin(phx)
c     phx=atan(b/a)
      implicit double precision (a-h, o-z)
      parameter (pi = 3.14159 26535 89793 23846 26433d0)
      complex*16 a, b, ampl, phx, temp

      aa = abs(a)
      bb = abs(b)
      if (aa+bb.eq. 0) then
         ampl=0.d0
         phx =0.d0
      elseif ( aa.gt.bb) then
         temp = b/a
         call atancc ( temp, phx)
         ampl = a / cos(phx)
      else
         temp = a/b
         call atancc ( temp, phx)
         phx = pi / 2 - phx
         ampl = b/sin(phx)
      endif

      if (dble(ampl).lt. 0.d0) then
         ampl = -ampl
         phx = phx + pi
      endif

      return
      end
      subroutine exjlnl (z, l, jl, nl)

c     purpose:  to calculate the spherical bessel functions jl and nl
c               for l = 0 to 6  using exact analytic expression
c
c     arguments:
c       z = argument of jl and nl
c       l = integer order of spherical bessel function
c       jl = jl bessel function (abramowitz conventions)
c       nl = nl bessel function (abramowitz yl conventions)
c            Note that this nl = abramowitz yl.
c
c       analytic expressions from abramowitz 10.1.11 and 10.1.12
c       recurrence relation to get analytic j4,n4  eqns 10.1.19-22 ala

      implicit double precision (a-h, o-z)

      complex*16 z, jl, nl

      complex*16 cosz, sinz

c     Exact formulae unstable for very small z, so use series
c     expansion there.  Limit of .3 chosen for 9 digit agreement.
      if (abs(z) .lt. 0.3)  then
         call bjnser (z, l, jl, nl, 0)
      else
c        use analytic formulae
         cosz = cos(z)
         sinz = sin(z)

         if (l .eq. 0)  then
            jl =  sinz / z
            nl = -cosz / z

         elseif (l .eq. 1)  then
            jl =  sinz/z**2 - cosz/z
            nl = -cosz/z**2 - sinz/z

         elseif (l .eq. 2)  then
            jl = ( 3/z**3 - 1/z)*sinz - 3*cosz/z**2
            nl = (-3/z**3 + 1/z)*cosz - 3*sinz/z**2

         elseif (l .eq. 3)  then
            jl = ( 15/z**4 - 6/z**2)*sinz + (-15/z**3 + 1/z)*cosz
            nl = (-15/z**4 + 6/z**2)*cosz + (-15/z**3 + 1/z)*sinz

         elseif (l .eq. 4)  then
            jl = ( 105/z**5 - 45/z**3 + 1/z )*sinz + 
     1                ( -105/z**4 + 10/z**2 )*cosz
            nl = (-105/z**5 + 45/z**3 - 1/z )*cosz + 
     1                ( -105/z**4 + 10/z**2 )*sinz

         elseif (l .eq. 5)  then
            jl = ( 945/z**6 - 420/z**4 + 15/z**2 )*sinz + 
     1              ( -945/z**5 + 105/z**3 - 1/z )*cosz
            nl = (-945/z**6 + 420/z**4 - 15/z**2 )*cosz + 
     1              ( -945/z**5 + 105/z**3 - 1/z )*sinz

         elseif (l .eq. 6)  then
            jl = ( 10395/z**7 - 4725/z**5 + 210/z**3 - 1/z )*sinz + 
     1              ( -10395/z**6 + 1155/z**4 - 21/z**2 )*cosz
            nl = (-10395/z**7 + 4725/z**5 - 210/z**3 + 1/z )*cosz + 
     1              ( -10395/z**6 + 1155/z**4 - 21/z**2 )*sinz

         else
            stop 'exjlnl, l out of range'
         endif
      endif

      return
      end
      subroutine polint( xa, ya, n, x, y, dy)
c     draws a polynimial P(x) of order (n-1) through n points.
c     returns y = P(x) and dy - estimate of the error
c     adapted  from numerical recipies in fortran by Press et al.

      implicit double precision (a-h,o-z)
      integer n, nmax
      parameter (nmax=4)
      dimension xa(nmax), ya(nmax), c(nmax), d (nmax)

      ns = 1
      dif = abs (x-xa(1))
      do 10 i=1,n
         dift = abs(x-xa(i))
         if (dift.lt.dif) then
            ns = i
            dif = dift
         endif
         c(i) = ya(i)
         d(i) = ya(i)
  10  continue
      y = ya(ns)
      ns = ns-1
      do 30 m=1,n-1
         do 20 i=1,n-m
            ho = xa(i)-x
            hp = xa(i+m)-x
            w = c(i+1) - d(i)
            den = ho-hp
            if (den.eq.0) pause 'failure in polint'
            den = w/den
            d(i) = hp*den
            c(i) = ho*den
  20     continue
         if (2*ns .lt. n-m) then
            dy = c(ns+1)
         else
            dy = d(ns)
            ns = ns-1
         endif
         y = y + dy
  30  continue

      return
      end
      function sdist (r0, r1)
c     find distance squared between cartesian points r0 and r1
c     single precision
      dimension r0(3), r1(3)
      sdist = 0
      do 10  i = 1, 3
         sdist = sdist + (r0(i) - r1(i))**2
   10 continue
      sdist = sqrt(sdist)
      return
      end
      subroutine somm (dr,dp,dq,dpas,da,m,np)
c
c integration by the method of simpson of (dp+dq)*dr**m from
c 0 to r=dr(np)
c dpas=exponential step;
c for r in the neighborhood of zero (dp+dq)=cte*r**da
c **********************************************************************
      implicit double precision (a-h,o-z)
      dimension dr(np), dp(np), dq(np)
      mm=m+1
      d1=da+mm
      da=0.0
      db=0.0
      do 70 i=1,np
      dl=dr(i)**mm
      if (i.eq.1.or.i.eq.np) go to 10
      dl=dl+dl
      if ((i-2*(i/2)).eq.0) dl=dl+dl
   10 dc=dp(i)*dl
      if (dc) 20,40,30
   20 db=db+dc
      go to 40
   30 da=da+dc
   40 dc=dq(i)*dl
      if (dc) 50,70,60
   50 db=db+dc
      go to 70
   60 da=da+dc
   70 continue
      da = dpas * (da + db) / 3.0
      dc=exp(dpas)-1.0
      db=d1*(d1+1.0)*dc*exp((d1-1.0)*dpas)
      db=dr(1)*(dr(2)**m)/db
      dc=(dr(1)**mm)*(1.0+1.0/(dc*(d1+1.0)))/d1
      da=da+dc*(dp(1)+dq(1))-db*(dp(2)+dq(2))
      return
      end
      subroutine somm2 (dr,dp,dpas,da,rnrm,m,np)
c Modified to use complex p and q.  SIZ 4/91
c Modified to use double simpson integration ALA 3/97
c integration by the method of simpson of dp*dr from 
c 0 to r=rnrm  with proper end corrections
c dpas=exponential step;
c for r in the neighborhood of zero dp=cte*r**da
c **********************************************************************
      implicit double precision (a-h,o-z)
      dimension dr(*)
      dimension  dp(*)

      mm = m + 1
      d1=dble(da)+mm
      da=0.0
      db=0.0
c      np-2=inrm -point of grid just below rnrm
      a1=log(rnrm/dr(np-2)) / dpas
      a2=a1**2/8.0d0
      a3=a1**3/12.0d0
      do 70 i=1,np
         if (i.eq.1) then
            dc=dp(i) *dr(i)**mm*9.0d0/24.0d0
         elseif (i.eq.2) then
            dc=dp(i) *dr(i)**mm*28.0d0/24.0d0
         elseif (i.eq.3) then
            dc=dp(i)*dr(i)**mm*23.0d0/24.0d0
         elseif (i.eq.np-3) then
            dc=dp(i)*dr(i)**mm*(25.0d0/24.0d0-a2+a3)
         elseif (i.eq.np-2) then
            dc=dp(i)*dr(i)**mm*(0.5d0+a1-3*a2-a3)
         elseif (i.eq.np-1) then
            dc=dp(i)*dr(i)**mm*(-1.0d0/24.0d0+5*a2-a3)
         elseif (i.eq.np) then
            dc=dp(i)*dr(i)**mm*(-a2+a3)
         else
c           like trapesoidal rule
            dc=dp(i)*dr(i)**mm
         endif
         da=da+dc
   70 continue
      da=dpas*da

c     add initial point (r=0) correction
      dd=exp(dpas)-1.0
      db=d1*(d1+1.0)*dd*exp((d1-1.0)*dpas)
      db=dr(1)*(dr(2)**m)/db
      dd=(dr(1)**mm)*(1.0+1.0/(dd*(d1+1.0)))/d1
      da=da+dd*dp(1)-db*dp(2)
      return
      end
      subroutine strap (x, y, n, sum)

c     Trapeziodal integration of y(x), result in sum
c     SINGLE PRECISION
c     modified by ala to handle cases for E<Efermi
c     sum only positive numbers

      dimension x(n), y(n)

      sum = y(1) * abs(x(2) - x(1))
      do 10  i = 2, n-1
         sum = sum + y(i) * abs(x(i+1) - x(i-1))
   10 continue
      sum = sum + y(n) * abs(x(n) - x(n-1))
      sum = sum/2

      return
      end
c     interpolation and extrapolation by m-th order polynomial
c     maximum m = 3. Change nmax if needed.
c     Input x and y arrays, returns y value y0 at requested x value x0.
c     Dies on error.

      subroutine terp (x, y, n, m, x0, y0)
      implicit double precision (a-h, o-z)

      dimension x(n), y(n)

c     Find out between which x points x0 lies
      i = locat (x0, n, x)
      k = min( max(i-m/2,1) , n-m )
      call polint( x(k), y(k), m+1, x0, y0, dy)

      return
      end

      function locat (x, n, xx)
      integer  u, m, n
      double precision x, xx(n)

c     Binary search for index of grid point immediately below x.
c     Array xx required to be monotonic increasing.
c     Returns
c     0            x <  xx(1)
c     1            x =  xx(1)
c     i            x =  xx(i)
c     n            x >= xx(n)

      locat = 0
      u = n+1

   10 if (u-locat .gt. 1)  then
         m = (u + locat) / 2
         if (x .lt. xx(m))  then
            u = m
         else
            locat = m
         endif
         goto 10
      endif

      return
      end


c     These routines, terp1 and locat1, are special versions to
c     be used with ff2chi, which uses some single and some double
c     precision.  They are the same as the routines in terp.f.

      subroutine terp1 (x, y, n, x0, y0)
      implicit double precision (a-h, o-z)

      real x(n), y(n)

c     Find out between which x points x0 lies
      i = locat1 (x0, n, x)
c     if i < 1, set i=1, if i > n-1, set i=n-1
      i = max (i, 1)
      i = min (i, n-1)

      if (x(i+1) - x(i) .eq. 0)  stop 'TERP-1'

      y0 = y(i) +  (x0 - x(i)) * (y(i+1) - y(i)) / (x(i+1) - x(i))

      return
      end

      function locat1 (x, n, xx)
      integer  u, m, n
      double precision x
      real xx(n)

c     Binary search for index of grid point immediately below x.
c     Array xx required to be monotonic increasing.
c     Returns
c     0            x <  xx(1)
c     1            x =  xx(1)
c     i            x =  xx(i)
c     n            x >= xx(n)

      locat1 = 0
      u = n+1

   10 if (u-locat1 .gt. 1)  then
         m = (u + locat1) / 2
         if (x .lt. xx(m))  then
            u = m
         else
            locat1 = m
         endif
         goto 10
      endif

      return
      end
c     interpolation and extrapolation by m-th order polynomial
c     maximum m = 3. Change nmax if needed.
c     Input x and y arrays, returns y value y0 at requested x value x0.
c     Dies on error.

      subroutine terpc (x, y, n, m, x0, y0)
      implicit double precision (a-h, o-z)

      complex*16 y, y0, dy
      dimension x(n), y(n)

c     Find out between which x points x0 lies
      i = locat (x0, n, x)
      k = min( max(i-m/2,1) , n-m )
      call polinc( x(k), y(k), m+1, x0, y0, dy)

      return
      end

      subroutine polinc( xa, ya, n, x, y, dy)
c     draws a polynimial P(x) of order (n-1) through n points.
c     returns y = P(x) and dy - estimate of the error
c     adapted  from numerical recipies in fortran by Press et al.

      implicit double precision (a-h,o-z)
      complex*16 ya,y,dy,c,d,w,den
      integer n, nmax
      parameter (nmax=4)
      dimension xa(nmax), ya(nmax), c(nmax), d (nmax)

      ns = 1
      dif = abs (x-xa(1))
      do 10 i=1,n
         dift = abs(x-xa(i))
         if (dift.lt.dif) then
            ns = i
            dif = dift
         endif
         c(i) = ya(i)
         d(i) = ya(i)
  10  continue
      y = ya(ns)
      ns = ns-1
      do 30 m=1,n-1
         do 20 i=1,n-m
            ho = xa(i)-x
            hp = xa(i+m)-x
            w = c(i+1) - d(i)
            den = ho-hp
            if (den.eq.0) stop 'failure in polint'
            den = w/den
            d(i) = hp*den
            c(i) = ho*den
  20     continue
         if (2*ns .lt. n-m) then
            dy = c(ns+1)
         else
            dy = d(ns)
            ns = ns-1
         endif
         y = y + dy
  30  continue

      return
      end
      subroutine trap (x, y, n, sum)
      implicit double precision (a-h, o-z)

c     Trapeziodal integration of y(x), result in sum

      dimension x(n), y(n)

      sum = y(1) * (x(2) - x(1))
      do 10  i = 2, n-1
         sum = sum + y(i) * (x(i+1) - x(i-1))
   10 continue
      sum = sum + y(n) * (x(n) - x(n-1))
      sum = sum/2

      return
      end
      SUBROUTINE CQdrtc(Coef,Sol,NSol)
c     Combutes the zeros of a quadratic polynomial
ccccccccccccccccccccccccccccccccccccccccccccccccc            
c     Input
c     Coef - array of coefficients
      COMPLEX*16 Coef(3)
ccccccccccccccccccccccccccccccccccccccccccccccccc
c     Output
c     Sol  - Array of solutions
c     NSol - # of solutions (only one if Coef(1) = 0 etc.)
c     NSol = -1 means a and b are zero
      COMPLEX*16 Sol(2)
      INTEGER NSol
ccccccccccccccccccccccccccccccccccccccccccccccccc
c     Local Variables
      COMPLEX*16 q, Sqrt
      DOUBLE PRECISION Sgn

      IF(Coef(1).eq.0.d0) THEN
         IF(Coef(2).eq.0.d0) THEN
            NSol = -1
            RETURN
         ELSE
            NSol = 1
            Sol(1) = -Coef(3)/Coef(2)
         END IF
      ELSE
         NSol = 2
         Root = Sqrt(Coef(2)**2-4.d0*Coef(1)*Coef(3))
         Sgn  = SIGN(DBLE(CONJG(Coef(2))*Root),1.d0)
         q    = -0.5d0*(Coef(2) + Sgn*Root)
         
         Sol(1) = q/Coef(1)
         Sol(2) = Coef(3)/q
      END IF

      RETURN
      END


      SUBROUTINE CCubic(Coef,Sol,NSol)
c     Combutes the zeros of a cubic polynomial
ccccccccccccccccccccccccccccccccccccccccccccccccc            
c     Input
c     Coef - array of coefficients
      COMPLEX*16 Coef(4)
ccccccccccccccccccccccccccccccccccccccccccccccccc
c     Output
c     Sol  - Array of solutions
c     NSol - # of solutions (only one if Coef(1) = 0 etc.)
c     NSol = -1 means a, b, and c are zero
      COMPLEX*16 Sol(4)
      INTEGER NSol
ccccccccccccccccccccccccccccccccccccccccccccccccc
c     Local Variables
      COMPLEX*16 P1, P2, Q, R, Coef2(3), a, b, c
      DOUBLE PRECISION Sgn, Theta
c     PARAMETERS
      COMPLEX*16 I
      PARAMETER(I = (0.d0, 1.d0))
      DOUBLE PRECISION Pi
      PARAMETER(Pi = 3.141592653589793238462643d0)

      IF(Coef(1).eq.0.d0) THEN
         Coef2(1) = Coef(2)
         Coef2(2) = Coef(3)
         Coef2(3) = Coef(4)         
         CALL CQdrtc(Coef2,Sol,NSol)
      ELSE
         a = Coef(2)/Coef(1)
         b = Coef(3)/Coef(1)
         c = Coef(4)/Coef(1)
         NSol = 3
         Q = (a**2 - 3.d0*b)/9.d0
         R = (2.d0*a**3 - 9.d0*a*b + 27.d0*c)/54.d0

         IF(((DIMAG(Q).eq.0.d0).and.(DIMAG(R).eq.0.d0)).and.
     &        (DIMAG(R**2).lt.DIMAG(Q**3))) THEN
            Theta = ACOS (DBLE(R/SQRT(Q**3)))
            Sol(1) = -2*SQRT(Q)*Cos(Theta/3.d0) - a/3.d0
            Sol(2) = -2*SQRT(Q)*Cos((Theta+2.d0*Pi)/3.d0) - a/3.d0
            Sol(3) = -2*SQRT(Q)*Cos((Theta-2.d0*Pi)/3.d0) - a/3.d0
         ELSE
            Sgn = SIGN(1.d0, DBLE(CONJG(R)*SQRT(R**2-Q**3)))
            P1 = -(R + Sgn*SQRT(R**2-Q**3))**(1.d0/3.d0)
            IF(P1.eq.0.d0) THEN
               P2 = 0.d0
            ELSE
               P2 = Q/P1
            END IF
            Sol(1) = (P1 + P2) - a/3.d0
            Sol(2) = -0.5d0*(P1 + P2) - a/3.d0 +
     &           I*SQRT(3.d0)/2.d0*(P1-P2)
            Sol(3) = -0.5d0*(P1 + P2) - a/3.d0 -
     &           I*SQRT(3.d0)/2.d0*(P1-P2)
         END IF
      END IF

      RETURN
      END
      
c///////////////////////////////////////////////////////////////////////
c PAR Subroutines
c Written by J. Sims, NIST, 2001

c This software was developed at the National Institute of Standards
c and Technology by employees of the Federal Government in the course
c of their official duties. Pursuant to title 17 Section 105 of
c the United States Code this software is not subject to copyright
c protection and is in the public domain. PAR is an experimental
c system. NIST assumes no responsibility whatsoever for its use by
c other parties, and makes no guarantees, expressed or implied, about 
c its quality, reliability, or any other characteristic. We would
c appreciate acknowledgement if the software is used.

c This software can be redistributed and/or modified freely provided
c that any derivative works bear some notice that they are derived from
c it, and any modified versions bear some notice that they have been
c modified.
c///////////////////////////////////////////////////////////////////////
c License is applicable for routines below, until otherwise specified.
c
c  **************************************************
c  Parallel feff8 routines
c  Jim Sims
c  **************************************************

      subroutine par_begin
c  **************************************************
c  Initializations for parallel version(s)
c  **************************************************

c={../HEADERS/parallel.h
      integer par_type, this_process, numprocs, my_rank
      logical master, worker, parallel_run
      real*8 wall_comm, time_comm
      common /timing/ wall_comm, time_comm
      common /parallel/ numprocs, my_rank, this_process, 
     .          master, worker, parallel_run, par_type
c= ../HEADERS/parallel.h}
c={mpif.h
c={mpif.h
      include 'mpif.h'
c= not found}
c= not found}

c-- So cvd or dbx can attach to a running process
c     call sleep(30) 

      call MPI_INIT(ierrorflag)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierrorflag)
      call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierrorflag)
      this_process = my_rank

      par_type = 0
      parallel_run = .true.
c-- The following variable will be used for IO that should only be
c-- done in one process.
      master = (my_rank .eq. 0)

      worker = (.not. master)
      if (worker) par_type = 1

c     write(6,*) 'this process = ',this_process, ' worker = ',worker

      if (master) write(6,*) 'Number of processors = ',numprocs

      return
      end

      subroutine par_stop (string)
c  **************************************************
c  Abnormal termination of the parallel session
c  **************************************************
c={../HEADERS/parallel.h
      integer par_type, this_process, numprocs, my_rank
      logical master, worker, parallel_run
      real*8 wall_comm, time_comm
      common /timing/ wall_comm, time_comm
      common /parallel/ numprocs, my_rank, this_process, 
     .          master, worker, parallel_run, par_type
c= ../HEADERS/parallel.h}
c={mpif.h
c={mpif.h
      include 'mpif.h'
c= not found}
c= not found}
c     For abnormal exits 
c     If open, close unit = 11
c     Go to the barrier that workers are sitting at
c     Then everyone will call par_end and stop
      logical is_open
      character*(*) string

      inquire(unit=11,opened=is_open)
      if (is_open) then
        call wlog(string)
        close(unit=11)
      else if (string .ne. ' ') then
	print *,string
	print *,'Abnormal termination on processor ',this_process
      endif
      call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierrorflag)

      stop ' '
      end

      subroutine par_end
c  **************************************************
c  Terminate the parallel session
c  **************************************************
      call MPI_FINALIZE(ierrorflag)
      return
      end

      subroutine par_barrier
c  **************************************************
c  Calls mpi_barrier
c  **************************************************
c={mpif.h
c={mpif.h
      include 'mpif.h'
c= not found}
c= not found}
      call MPI_BARRIER(MPI_COMM_WORLD,ierrorflag)
      return
      end

      subroutine par_send_int(buf,count,dest,tag)
c  **************************************************
c  Call mpi_send for integer arrays
c  **************************************************
      integer count,dest,tag
      integer buf(*)
c={mpif.h
c={mpif.h
      include 'mpif.h'
c= not found}
c= not found}
      call MPI_SEND(buf,count,MPI_INTEGER,dest,tag,
     .              MPI_COMM_WORLD,ierrorflag)
      return
      end

      subroutine par_send_cmplx(buf,count,dest,tag)
c  **************************************************
c  Call mpi_send for complex arrays
c  **************************************************
      integer count,dest,tag
      complex buf(*)
c={mpif.h
c={mpif.h
      include 'mpif.h'
c= not found}
c= not found}
      call MPI_SEND(buf,count,MPI_COMPLEX,dest,tag,
     .              MPI_COMM_WORLD,ierrorflag)
      return
      end

      subroutine par_send_dc(buf,count,dest,tag)
c  **************************************************
c  Call mpi_send for double_complex arrays
c  **************************************************
      integer count,dest,tag
      complex*16 buf(*)
c={mpif.h
c={mpif.h
      include 'mpif.h'
c= not found}
c= not found}
      call MPI_SEND(buf,count,MPI_DOUBLE_COMPLEX,dest,tag,
     .              MPI_COMM_WORLD,ierrorflag)
      return
      end

      subroutine par_recv_int(buf,count,source,tag)
c  **************************************************
c  Call mpi_recv for integer arrays
c  **************************************************
      integer count,source,tag
      integer buf(*)
c={mpif.h
c={mpif.h
      include 'mpif.h'
c= not found}
c= not found}
      integer istat(mpi_status_size)
      call MPI_RECV(buf,count,MPI_INTEGER,source,tag,
     .              MPI_COMM_WORLD,istat,ierrorflag)
      return
      end

      subroutine par_recv_cmplx(buf,count,source,tag)
c  **************************************************
c  Call mpi_recv for complex arrays
c  **************************************************
      integer count,source,tag
      complex buf(*)
c={mpif.h
c={mpif.h
      include 'mpif.h'
c= not found}
c= not found}
      integer istat(mpi_status_size)
      call MPI_RECV(buf,count,MPI_COMPLEX,source,tag,
     .              MPI_COMM_WORLD,istat,ierrorflag)
      return
      end

      subroutine par_recv_dc(buf,count,source,tag)
c  **************************************************
c  Call mpi_recv for double complex arrays
c  **************************************************
      integer count,source,tag
      complex*16 buf(*)
c={mpif.h
c={mpif.h
      include 'mpif.h'
c= not found}
c= not found}
      integer istat(mpi_status_size)
      call MPI_RECV(buf,count,MPI_DOUBLE_COMPLEX,source,tag,
     .              MPI_COMM_WORLD,istat,ierrorflag)
      return
      end

      subroutine par_bcast_int(buf,count,source)
c  **************************************************
c  Call mpi_bcast for integer arrays
c  **************************************************
      integer count,source
      integer buf(*)
c={mpif.h
c={mpif.h
      include 'mpif.h'
c= not found}
c= not found}
      call MPI_BCAST(buf,count,MPI_INTEGER,source,
     .              MPI_COMM_WORLD,ierrorflag)
      return
      end

      subroutine par_bcast_cmplx(buf,count,source)
c  **************************************************
c  Call mpi_bcast for complex arrays
c  **************************************************
      integer count,source
      complex buf(*)
c={mpif.h
c={mpif.h
      include 'mpif.h'
c= not found}
c= not found}
      call MPI_BCAST(buf,count,MPI_COMPLEX,source,
     .              MPI_COMM_WORLD,ierrorflag)
      return
      end

      subroutine par_bcast_dc(buf,count,source)
c  **************************************************
c  Call mpi_bcast for double_complex arrays
c  **************************************************
      integer count,source
      complex*16 buf(*)
c={mpif.h
c={mpif.h
      include 'mpif.h'
c= not found}
c= not found}
      call MPI_BCAST(buf,count,MPI_DOUBLE_COMPLEX,source,
     .              MPI_COMM_WORLD,ierrorflag)
      return
      end

      subroutine MPE_DECOMP1D( n, num_procs, myid, s, e )
c  ******************************************************
c  A routine for producing a decomposition of a 1-d 
c  array when given a number of processors.  It may 
c  be used in "direct" product decomposition.  The 
c  values returned assume a "global" domain in [1:n]
c  ******************************************************
c  MPE_Decomp1d - Compute a balanced decomposition of
c  a 1-D array
c  ******************************************************
c  Input Parameters:
c  n  - Length of the array
c  num_procs - Number of processors in decomposition
c  myid  - Rank of this processor in the decomposition 
c  (0 <= rank < size)
c  ******************************************************
c  Output Parameters:
c  s,e - Array my_particles are s:e, with the original 
c  array considered as 1:n.  
c  ******************************************************

      integer n, num_procs, myid, s, e
      integer nloc, deficit
 
      nloc  = n / num_procs
      s       = myid * nloc + 1
      deficit = mod(n,num_procs)
      s       = s + min(myid,deficit)
      if (myid .lt. deficit) then
        nloc = nloc + 1
      endif
      e = s + nloc - 1
      if (e .gt. n .or. myid .eq. num_procs-1) e = n

      return
      end

      SUBROUTINE SECONDS( W)
c  ***************************************************
c  SECONDS returns the wall clock times for a process
c  in seconds.
c  ***************************************************

      real*8 W, MPI_Wtime

      W = MPI_Wtime()

      RETURN
      END
