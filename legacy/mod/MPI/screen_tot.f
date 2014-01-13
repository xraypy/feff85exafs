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
c=======================================================================
c     FFMOD8
c=======================================================================

c     sub-program screen
      program  ffmod8
c     subroutine ffmod8
c     Written by Y. Takimoto

      implicit none
c={../HEADERS/declare.h
c     const.h
      double precision pi
      integer one, zero
      double precision third
      double precision raddeg
      double precision fa
      double precision bohr, ryd
      double precision hart
      double precision alpinv
      double precision alphfs

c     dim.h
      integer nclusx
      integer nclxtd
      integer nspx
      integer natx
      integer nattx
      integer lx
      integer nphx
      integer ltot
      integer nrptx
      integer nex, lamtot
      integer mtot, ntot
      integer npatx
      integer legtot
      integer novrx
      integer nheadx
      integer MxPole
c= ../HEADERS/declare.h}
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
c={../HEADERS/parallel.h
      integer par_type, this_process, numprocs, my_rank
      logical master, worker, parallel_run
      real*8 wall_comm, time_comm
      common /timing/ wall_comm, time_comm
      common /parallel/ numprocs, my_rank, this_process, 
     .          master, worker, parallel_run, par_type
c= ../HEADERS/parallel.h}
      
c     geom.inp -------------------------------------------------------
c     number of atoms in problem
      integer nat
c     number of unique potentials      
      integer nph
c     Unique potential input data
c     given unique pot, which atom is model?
c     (0 if none specified for this unique pot)
      integer iatph(0:nphx) 
                             
      integer ibounc(natx)
c     Specific atom input data
c     given specific atom, which unique pot?
      integer iphat(natx)
c     cartesian coords of specific atom
      double precision rat(3,natx)
c     ----------------------------------------------------------------
      integer ixc0
      double precision dx, x0, vr0
c     ldos.inp -------------------------------------------------------
      integer mldos, lfms2, ixc, ispin, minv
      double precision emin, emax, eimag, rgrd
      real rfms2, rdirec, toler1, toler2
      integer  lmaxph(0:nphx)
c     -----------------------------------------------------------------
      double precision ri(nrptx)
c     W, the result of calculation
      double precision vch(nrptx), wscrn(nrptx)
      integer i, ilast, nrx
      parameter ( nrx = 251 )
      
      double precision wall_start, wall_end
      integer ios
c     mod1.inp ---------------------------------------------------------
      integer nohole
      character dummy, mpot
c      call par_begin
c
c      ! Initialize clock
c      call seconds(wall_start)
c      wall_comm = 0.0
c
c      ! open the log file, unit 11.  See subroutine wlog.
c      if (master) then
c        open (unit=11, file='logscrn.dat', status='unknown', iostat=ios)
c        call chopen (ios, 'logscrn.dat', 'feff')
c      else
c        par_type = 2
c      endif

      call par_begin
      if (worker) go to 400

c     Josh - if nohole eq 2 run SCREEN, else skip.
      open (unit=12, file='mod1.inp', status='unknown', iostat=ios)
      call chopen (ios, 'mod1.inp', 'feff')
      read(12,'(a)') dummy
      read(12,*) mpot
      read(12,'(a)') dummy
      read(12,*) nohole, nohole
      
      if (nohole.ne.2) goto 5
c     Josh END

      open (unit=11, file='logscrn.dat', status='unknown', iostat=ios)
      call chopen (ios, 'logscrn.dat', 'feff')
      
c=================================================================
c     read  INPUT data files: geom.dat, mod2.inp, and ldos.inp.
c=================================================================
      call rdgeom(nat, nph, iatph, ibounc, iphat, rat) ! read geom.dat
      call rdldos(mldos, lfms2, ixc, ispin, minv,emin, emax, eimag,rgrd,&
     &            rfms2, rdirec, toler1, toler2, lmaxph)
            
      call wlog(' Calculating SCREEN ...')
c=================================================================
c     Prepare, then start the calculation
c=================================================================
c     make radial grid
      dx   = 0.05d0
      x0   = 8.8d0
      call setri(dx, x0, nrx, ri)
      
      vr0  = 0.0d0
      ixc0 = 0

c geom.dat => first line, ldos.inp => second line      
      call prep( nat, nph, iatph, ibounc, iphat, rat,                   &
     &           ixc, toler1, toler2, lmaxph, lfms2,                    &
     &           vr0, ixc0, nrx, ri, dx, x0, ilast, vch, wscrn)
     
      call wlog(' Done with SCREEN.') 
c=================================================================
c     open files for output
c=================================================================
      open (unit=24, file='wscrn.dat', status='unknown', iostat=ios)
      
      do i = 1, ilast
        write(24, '(33e20.10)') ri(i), wscrn(i), vch(i)
      end do
     
c     close opened files
      close (unit=24)
      close (unit=11)

c      if (master) close (unit=11)
c
c      ! Time at end of run
c      call seconds(wall_end)
c      if (master .and. parallel_run) then
c        write (6,*) 'total time    ', wall_end - wall_start
c        write (6,*) 'communicate time', wall_comm
c      endif
c      call par_end
  5   continue

  400 call par_barrier
      call par_end
c     sub-pro exchange point
      stop
c     return

      end
c=======================================================================
c     END FFMOD8
c=======================================================================

c={./rdgeom.f
c=======================================================================
c     Read GEOM.INP
c=======================================================================
      subroutine rdgeom(nat, nph, iatph, ibounc, iphat, rat)
      implicit none
c={../HEADERS/declare.h
c     const.h
      double precision pi
      integer one, zero
      double precision third
      double precision raddeg
      double precision fa
      double precision bohr, ryd
      double precision hart
      double precision alpinv
      double precision alphfs

c     dim.h
      integer nclusx
      integer nclxtd
      integer nspx
      integer natx
      integer nattx
      integer lx
      integer nphx
      integer ltot
      integer nrptx
      integer nex, lamtot
      integer mtot, ntot
      integer npatx
      integer legtot
      integer novrx
      integer nheadx
      integer MxPole
c= ../HEADERS/declare.h}
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

c     number of atoms in problem
      integer nat
c     number of unique potentials      
      integer nph
c     Unique potential input data
c     given unique pot, which atom is model?
c     (0 if none specified for this unique pot)
      integer iatph(0:nphx)
      
      integer ibounc(natx)
c     Specific atom input data
c     given specific atom, which unique pot?
      integer iphat(natx)
c     cartesian coords of specific atom
      double precision rat(3,natx)

c      Local stuff
      character*512 slog
      character*80 head(10)
      integer lhead(10)
c     character(len=80), dimension(nheadx) :: head
c     integer(I4), dimension(nheadx)       :: lhead
      integer nhead, idum, ios, i, j, i1b

c     Read  geom.dat file
      open (file='geom.dat', unit=3, status='old',iostat=ios)
c       read header from geom.dat
        nhead = nheadx
        call rdhead (3, nhead, head, lhead)
        nat = 0
        nph = 0
        do i = 0, nphx
          iatph(i) = 0
        end do
        
        do i = 1, natx
          nat = nat + 1
          
          if (nat .gt. natx)  then
            write(slog, '(a, 2i10)') ' nat, natx ', nat, natx
            call wlog(slog)
            stop 'Bad input'
          endif
          read(3,*,end=60)  idum, (rat(j,nat),j=1,3), iphat(nat), i1b
          if ( iphat(nat) .gt. nph) nph = iphat(nat)
          if ( iatph(iphat(nat)) .eq. 0 ) then
            iatph(iphat(nat)) = nat
          end if
        end do
        
  60    continue
  
        nat = nat - 1
      close(3)

      do i = 1, nat
        do j = 1, 3
          rat(j,i) = rat(j,i) / bohr
        end do
      end do

      return
c      end subroutine rdgeom
      end
c= ./rdgeom.f}
c={./rdldos.f
c=======================================================================
c     Read LDOS.INP
c=======================================================================
      subroutine rdldos(mldos, lfms2, ixc, ispin, minv,                 &
     &                   emin, emax, eimag, rgrd,                       &
     &                   rfms2, rdirec, toler1, toler2, lmaxph)
      implicit none
c={../HEADERS/declare.h
c     const.h
      double precision pi
      integer one, zero
      double precision third
      double precision raddeg
      double precision fa
      double precision bohr, ryd
      double precision hart
      double precision alpinv
      double precision alphfs

c     dim.h
      integer nclusx
      integer nclxtd
      integer nspx
      integer natx
      integer nattx
      integer lx
      integer nphx
      integer ltot
      integer nrptx
      integer nex, lamtot
      integer mtot, ntot
      integer npatx
      integer legtot
      integer novrx
      integer nheadx
      integer MxPole
c= ../HEADERS/declare.h}
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
      
      integer mldos, lfms2, ixc, ispin, minv
      double precision emin, emax, eimag, rgrd
      real rfms2, rdirec, toler1, toler2
      integer  lmaxph(0:nphx)
      character*26 string(18)

c     Local stuff
      character*512 slog
      integer ios, iph

c     standard formats for string, integers and real numbers
  10  format(a)
  20  format(20i4)
  30  format(6f13.5)

c     read ldos.inp
      open (file='ldos.inp', unit=3, status='old',iostat=ios)
        read (3,10)  slog
        read (3,20)  mldos, lfms2, ixc, ispin, minv
        read (3,10)  slog
        read (3,30)  rfms2, emin, emax, eimag, rgrd
        read (3,10)  slog
        read (3,30)  rdirec, toler1, toler2
        read (3,10)  slog
        read (3,20)  lmaxph
      close(3)

c     transform to code units (bohrs and hartrees - atomic unuts)
      rfms2  = rfms2  / bohr
      rdirec = rdirec / bohr
      emin   = emin   / hart
      emax   = emax   / hart
      eimag  = eimag  / hart

      return
c      end subroutine rdldos
      end
c= ./rdldos.f}
c={./rdscrn.f
c=======================================================================
c     Read SCRN.INP
c=======================================================================
      subroutine rdscrn(ner, nei, maxl, irrh, iend, lfxc, nrptx0,       &
     &                   emin, emax, eimax, ermin, rfms)
      implicit none
c={../HEADERS/declare.h
c     const.h
      double precision pi
      integer one, zero
      double precision third
      double precision raddeg
      double precision fa
      double precision bohr, ryd
      double precision hart
      double precision alpinv
      double precision alphfs

c     dim.h
      integer nclusx
      integer nclxtd
      integer nspx
      integer natx
      integer nattx
      integer lx
      integer nphx
      integer ltot
      integer nrptx
      integer nex, lamtot
      integer mtot, ntot
      integer npatx
      integer legtot
      integer novrx
      integer nheadx
      integer MxPole
c= ../HEADERS/declare.h}
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
    
      integer ner, nei, maxl, irrh, iend, lfxc, nrptx0
      double precision emin, emax, eimax, ermin, rfms
        
c     === Local stuff ================================================
      integer ios, i, j
      character*512 slog
      character*8 strs(50)
      character*3 str
      double precision vars(5)
c     standard formats for string, integers and real numbers
      character*15 fmt_s, fmt_i, fmt_r
      parameter (fmt_s  = '(a)')
      parameter (fmt_i  = '(a8, i4)')
      parameter (fmt_r  = '(a8, f13.5)')
c     =================================================================
c     set default values
      ner   = 100
      nei   = 20
      maxl  = 4
      irrh  = 1
      iend  = 0
      emin  = -40.0d0
      emax  = 0.0d0
      eimax = 4.0d0
      ermin = 0.001d0
      lfxc  = 0
      rfms  = 4.0d0
      nrptx0 = 251
      
c     read scrn.inp
      open (file='scrn.inp', unit=3, status='old', iostat=ios)
        if (ios .ne. 0) goto 60
        read (3,fmt_s) slog
        
        do i = 1, 12
        
          read(3,fmt_r,end=60)  strs(i), vars(i)
          str = strs(i)

          if (str .eq. 'ner') ner   = vars(i)
          if (str .eq. 'nei') nei   = vars(i)
          if (str .eq. 'max') maxl  = vars(i)
          if (str .eq. 'irr') irrh  = vars(i)
          if (str .eq. 'ien') iend  = vars(i)
          if (str .eq. 'lfx') lfxc  = vars(i)
          if (str .eq. 'emi') emin  = vars(i)
          if (str .eq. 'ema') emax  = vars(i)
          if (str .eq. 'eim') eimax = vars(i)
          if (str .eq. 'erm') ermin = vars(i)
          if (str .eq. 'rfm') rfms  = vars(i)
          if (str .eq. 'nrp') nrptx0  = vars(i)
      
        end do

      close(3)

60    continue

c     transform to code units (bohrs and hartrees - atomic unuts)
      rfms   = rfms   / bohr
      emin   = emin   / hart
      emax   = emax   / hart
      eimax  = eimax  / hart
      ermin  = ermin  / hart
      
      return
c      end subroutine rdscrn
      end
c= ./rdscrn.f}
c={./frgrid.f
c=======================================================================
c     RADIAL GRID ROUTINES
c=======================================================================
      subroutine setri(dx, x0, ilast, ri)
        implicit none
c={../HEADERS/declare.h
c     const.h
      double precision pi
      integer one, zero
      double precision third
      double precision raddeg
      double precision fa
      double precision bohr, ryd
      double precision hart
      double precision alpinv
      double precision alphfs

c     dim.h
      integer nclusx
      integer nclxtd
      integer nspx
      integer natx
      integer nattx
      integer lx
      integer nphx
      integer ltot
      integer nrptx
      integer nex, lamtot
      integer mtot, ntot
      integer npatx
      integer legtot
      integer novrx
      integer nheadx
      integer MxPole
c= ../HEADERS/declare.h}
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
        
        double precision ri(nrptx), dx, x0
        integer ilast, i

        
        do i = 1, ilast
          ri(i) = exp(-x0+(i-1)*dx)
        end do

c      end subroutine setri
      end

      function getiat(x0, dx, rref)
        implicit none
        integer getiat
        double precision rref, x0, dx
        
        getiat = (log(rref) + x0) / dx + 1
      
      return
c      end function getiat
      end
c=======================================================================
c     END RADIAL GRID ROUTINES
c=======================================================================
c= ./frgrid.f}
c={./fegrid.f
c=======================================================================
c     ENERGY GRID ROUTINES
c=======================================================================            
      subroutine setegi(emin, emax, eimax, ermin, ner, nei, em, ne)
      implicit none
c={../HEADERS/declare.h
c     const.h
      double precision pi
      integer one, zero
      double precision third
      double precision raddeg
      double precision fa
      double precision bohr, ryd
      double precision hart
      double precision alpinv
      double precision alphfs

c     dim.h
      integer nclusx
      integer nclxtd
      integer nspx
      integer natx
      integer nattx
      integer lx
      integer nphx
      integer ltot
      integer nrptx
      integer nex, lamtot
      integer mtot, ntot
      integer npatx
      integer legtot
      integer novrx
      integer nheadx
      integer MxPole
c= ../HEADERS/declare.h}
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
c     makes a grid in complex e-plane like below
c     en_grid is comlex energy in hartrees
c
c     ...................
c     .                 .
c     1                 n
c     ---------------------------->Re (Im = 0)
      double precision emin, emax, eimax, ermin
      complex*16 em(nex)
      integer ner, nei, ne

c     Local
      double precision s_real, csreal, csimag
      complex*16 de, dei, delta, em2(nex)
      integer i, neg
      
c     never do calculations on real axis.
      if (ermin .le. 0.0) ermin = 0.05d0
      s_real = emax - emin

c     differential energy step along real axis
      de  = (emax - emin) / (ner - 1)
c     differential energy step along imaginary axis
      dei = coni*(eimax - ermin) / (nei - 1)
      
      neg = 1

c=================================================================
c     Prepare Grids
c=================================================================      
      em2(neg) = emax + coni*ermin
      csreal  = 0.0d0
      csimag  = ermin
      delta   = dei
            
      do neg = 2, nex**2
c        neg = neg + 1
        
        if (dreal(em2(neg-1)) .lt. emin) then
c         going down along the imag axis
          delta = -dei
          if (dimag(em2(neg-1)) .le. ermin) then
            goto 60
          end if
        else
          if (abs(csimag) .ge. eimax) then
c           change direction, steps along the real axis
            delta  = -de
            csimag = 0.0d0
          end if
        end if
        
        csimag = csimag + abs(dimag(delta))
        csreal = csreal + abs(dreal(delta))
        
        em2(neg) = em2(neg-1) + delta
        
      end do

60    continue

      if (dimag(em2(neg-1)) .le. 0.0) then
        neg = neg - 2
      else
        neg = neg - 1
      end if
      
c     reverse the grid      
      do i = 1, neg
        em(i) = em2(neg+1-i)
      end do

      if (neg .gt. nex) then
        print *, 'Error: too many energy points:', neg
        stop
      end if
            
      ne = neg

      return
      
c      end subroutine setegi
      end
c=======================================================================
c     END ENERGY GRID ROUTINES
c=======================================================================
c= ./fegrid.f}
c={./fxc.f
c=======================================================================
c     LDA EXCHANGE
c=======================================================================
      subroutine ldafxc(ilast, ri, edens, exchange, fxc)
      implicit none
c={../HEADERS/declare.h
c     const.h
      double precision pi
      integer one, zero
      double precision third
      double precision raddeg
      double precision fa
      double precision bohr, ryd
      double precision hart
      double precision alpinv
      double precision alphfs

c     dim.h
      integer nclusx
      integer nclxtd
      integer nspx
      integer natx
      integer nattx
      integer lx
      integer nphx
      integer ltot
      integer nrptx
      integer nex, lamtot
      integer mtot, ntot
      integer npatx
      integer legtot
      integer novrx
      integer nheadx
      integer MxPole
c= ../HEADERS/declare.h}
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

      integer ilast, exchange
      double precision ri(nrptx), edens(251)
      double precision fxc(nrptx)
      double precision rs
      integer i
    
      do i= 1, ilast
c       also calculate the local exchange term fxc = vxc*r_s**3/r**2
        if (edens(i) .le. 0) then
          rs = 100
          fxc(i) = 0
        else
          rs = (edens(i)/3)**(-third)
          fxc(i) = rs**3 / ri(i)**2 / 6 * (-1.222/rs -0.75924/(11.4+rs))
          if (exchange .eq. 2) then
            fxc(i) = rs**3 / ri(i)**2 / 6 * (-1.222/rs)
          end if
        endif
c     vvbh from Von Barth Hedin paper, 1971
c     see eq.60-61 of Gross&Kohn, Adv. Quant. Chem. 21, p.255(1990).
c     eq.60 fxc = d V_xc / d rho * (2\ell + 1) /(4*pi*r**2)
c     where second factor comes from \delta(r-r')
c     for dipole field \ell=1 
c     fxc(i) = rs**3 / ri(i)**2 / 6 * (-1.22177412/rs -1.512/(30+rs))
c     c below are coefficients in Zangwill/Soven paper
      end do

      return
c      end subroutine ldafxc
      end
c= ./fxc.f}
c={./prep.f
c=======================================================================
c     PREP
c=======================================================================
c geom.dat => first line, ldos.inp => second line 
      subroutine prep (nat, nph, iatph, ibounc, iphat, rat,             &
     &                  ixc, toler1,toler2,lmaxph, lfms2,               &
     &                  vr0, ixc0, nrx, ri, dx, x0,                     &
     &                  ilast, vch, wscrn)

      implicit none
c={../HEADERS/declare.h
c     const.h
      double precision pi
      integer one, zero
      double precision third
      double precision raddeg
      double precision fa
      double precision bohr, ryd
      double precision hart
      double precision alpinv
      double precision alphfs

c     dim.h
      integer nclusx
      integer nclxtd
      integer nspx
      integer natx
      integer nattx
      integer lx
      integer nphx
      integer ltot
      integer nrptx
      integer nex, lamtot
      integer mtot, ntot
      integer npatx
      integer legtot
      integer novrx
      integer nheadx
      integer MxPole
c= ../HEADERS/declare.h}
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
      
c     =================================================================
c     function getiat
      integer getiat
c     geom.inp -------------------------------------------------------
c     number of atoms in problem
      integer nat
c     number of unique potentials
      integer nph
c     Unique potential input data
c     given unique pot, which atom is model?
c     (0 if none specified for this unique pot)
      integer iatph(0:nphx)
       
      integer ibounc(natx)
c     Specific atom input data
c     given specific atom, which unique pot? 
      integer iphat(natx)
c     cartesian coords of specific atom
      double precision rat(3,natx)
c     -----------------------------------------------------------------
      real srat(3,natx)
c     mod2.inp --------------------------------------------------------
      double precision vr0
      integer ixc0
c     ldos.inp
      integer ixc, lfms2
      real toler1, toler2, rdirec, rfms2
c     scrn.inp
      integer ner, nei, maxl, irrh, iend, lfxc, nrptx0
      double precision emin, emax, eimax, ermin, rfms
c     pot.bin --------------------------------------------------------
      integer ntitle
      character*80 title(nheadx)
      integer nohole, ihole, inters, iafolp, jumprm, iunf
      double precision rnrmav, xmu, vint, rhoint
      double precision emu, s02, erelax, wp, ecv, rs, xf, qtotel, totvol
      integer imt(0:nphx), inrm(0:nphx)
      double precision rmt(0:nphx), rnrm(0:nphx)
      double precision folp(0:nphx),folpx(0:nphx)
      double precision dgc0(251), dpc0(251)
      double precision dgc(251, 30, 0:nphx), dpc(251, 30, 0:nphx)
      double precision adgc(10, 30, 0:nphx), adpc(10, 30, 0:nphx)
      double precision edens(251, 0:nphx), vclap(251, 0:nphx)
      double precision vtot(251, 0:nphx), edenvl(251, 0:nphx)
      double precision vvalgs(251, 0:nphx), dmag(251, 0:nphx)
      double precision xnval(30,0:nphx)
      double precision eorb(30)
      integer kappa(30)
      double precision qnrm(0:nphx)
      double precision xnmues(0:lx,0:nphx)
      integer iorb(-4:3,0:nphx)
      integer iz(0:nphx)
      double precision xion(0:nphx)
      double precision xnatph(0:nphx)
c     -----------------------------------------------------------------
      double precision vtotph(nrptx,0:nphx), vvalph(nrptx,0:nphx)
      double precision rhoph(nrptx), rhphvl(nrptx), dmagx(nrptx)
      double precision dgcx(nrptx), dpcx(nrptx)
      double precision dgcn(nrptx,30,0:nphx), dpcn(nrptx,30,0:nphx)
      integer iph, jnew
      double precision dx, x0, edge, vjump
      double precision ri(nrptx)
      
      integer jri, jri1, ne, ne1, ne3, i, j
      integer inclus(0:nphx), lmax(0:nphx)
      complex*16 ph(nex, ltot+1, 0:nphx)
      complex gtr(0:lx, 0:nphx, nex)
      double precision de, enext
      complex*16 eref(nex), em(nex)
      integer lmaxph(0:nphx)
      integer nrx, ilast
c     W, the result of calculation
      double precision vch(nrx), wscrn(nrx)
c     =================================================================      
c     read scrn.inp
      call rdscrn(ner, nei, maxl, irrh, iend, lfxc, nrptx0,             &
     &                  emin, emax, eimax, ermin, rfms)
c     read pot.bin
      call rdpot ( ntitle, title, rnrmav, xmu, vint, rhoint,            &
     &                  emu, s02, erelax, wp, ecv,rs,xf, qtotel,        &
     &                  imt, rmt, inrm, rnrm, folp, folpx, xnatph,      &
     &                  dgc0, dpc0, dgc, dpc, adgc, adpc,               &
     &                  edens, vclap, vtot, edenvl, vvalgs, dmag, xnval,&
     &                  eorb, kappa, iorb, qnrm, xnmues, nohole, ihole, &
     &                  inters, totvol, iafolp, xion, iunf, iz, jumprm)
          
      rfms2  = rfms
      rdirec = rfms*2
      
      edge = xmu - vr0
      emu  = emu - vr0
      
c     make energy grid
      call setegi(xmu+emin, xmu+emax, eimax, ermin, ner, nei, em, ne)
      
c     compute phase shifts
      do iph = 0, nph
        print *, '     get phase shift for potential type ', iph
        lmax(iph) = lx
      
        call fixvar (rmt(iph),edens(1,iph),vtot(1,iph),dmag(1,iph),vint,&
     &           rhoint,dx,dx,jumprm,vjump,ri,vtotph(1,iph),rhoph,dmagx)
        call fixdsx (iph, dx, dx, dgc, dpc, dgcn(1,1,iph),dpcn(1,1,iph))

        jri     = getiat(x0, dx, rmt(iph)) + 1
        jri1    = jri + 1
        eref(1) = vtotph(jri1,iph)
        do i = 1, jri1
          vtotph(i,iph) = vtotph(i,iph) - eref(1)
        end do
        
        if (ixc .ge. 5) then
          do i = 1, jri1
            vvalph(i,iph) = vvalph(i,iph) - eref(1)
          end do
        else
          do i = 1, jri1
            vvalph(i,iph) = vtotph(i,iph)
          end do
        endif
         
        call getph( ri, dx, x0, rmt(iph), rnrm(iph), ne, em, ixc,       &
     &             vtotph(1,iph), vvalph(1,iph),                        &
     &             dgcn(1,1,iph), dpcn(1,1,iph), eref(1),               &
     &             adgc(1,1,iph), adpc(1,1,iph),                        &
     &             iz(iph), xion(iph), iunf,                            &
     &             ihole, lx, xnval(1,iph), ph(1,1,iph))
      end do

c=================================================================
c     Main part of calculation
c=================================================================
      iph = 0
c     transform to single precision
      do i = 1, nat
        do j = 1, 3
          srat(j,i) = real(rat(j,i))
        end do
      end do
      
c     only the potential type 0 (central atom) is evaluated
c geom.dat => 1st, ldos.inp => 2nd, scrn.inp => 3rd, pot.bin => 4,5th
      call screen( nat, nph, iphat, srat,                               &
     &            toler1, toler2, lmaxph, rdirec, rfms2, lfms2,         &
     &            maxl, irrh, iend, lfxc, nrptx0,                       &
     &            ihole, xmu, adgc, adpc, iz,                           &
     &            xion, iunf, xnval, edens, dgc0, dpc0,                 &
     &            ri, dx, x0, rmt, rnrm, em, ne,                        &
     &            ixc0, iph, vtotph, vvalph, eref, dgcn, dpcn, ph,      &
     &            ilast, vch, wscrn)
  
      return
c      end subroutine prep
      end
c=======================================================================
c     END PREP
c=======================================================================
c= ./prep.f}
c={./getph.f
      subroutine getph ( ri, dx, x0, rmt, rnrm, ne, em, ixc,            &
     &                  vtot, vvalgs,                                   &
     &                  dgcn, dpcn, eref, adgc, adpc, iz, xion, iunf,   &
     &                  ihole, lmaxsc, xnval, ph)
      
      implicit none
c={../HEADERS/declare.h
c     const.h
      double precision pi
      integer one, zero
      double precision third
      double precision raddeg
      double precision fa
      double precision bohr, ryd
      double precision hart
      double precision alpinv
      double precision alphfs

c     dim.h
      integer nclusx
      integer nclxtd
      integer nspx
      integer natx
      integer nattx
      integer lx
      integer nphx
      integer ltot
      integer nrptx
      integer nex, lamtot
      integer mtot, ntot
      integer npatx
      integer legtot
      integer novrx
      integer nheadx
      integer MxPole
c= ../HEADERS/declare.h}
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

c     INPUT
c     dx, x0, ri(nr)
c                  Loucks r-grid, ri=exp((i-1)*dx-x0)
c     ne, em(ne)   number of energy points,  complex energy grid
c     ixc          0  Hedin-Lunqist + const real & imag part
c                  1  Dirac-Hara + const real & imag part
c                  2  ground state + const real & imag part
c                  3  Dirac-Hara + HL imag part + const real & imag part
c                  5  Dirac-Fock exchange with core electrons +
c                     ixc=0 for valence electron density
c     rmt          r muffin tin
c     rnrm         r norman
c     vtot(nr)     total potential, including gsxc, final state
c     dgcn(dpcn)   large (small) dirac components for central atom
c     adgc(adpc)   their development coefficients
c
c     OUTPUT
c     ph

c     function getiat
      integer getiat
      
      integer ne, ixc, iunf, ihole, lmaxsc, iz
      double precision xion
c     output
      complex*16 ph(nex, ltot+1)

      double precision ri(nrptx)
      double precision vtot(nrptx), vvalgs(nrptx)
      complex*16 vtotc(nrptx), vvalc(nrptx)
      double precision xnval(30)
      double precision dgcn(nrptx,30), dpcn(nrptx,30)
      double precision adgc(10,30), adpc(10,30)

c     energy grid in complex e-plane
      complex*16 em(nex), eref

c     work space for dfovrg: regular and irregular solutions
      complex*16 pn(nrptx), qn(nrptx)
      complex*16 p2, xkmt, ck
      complex*16 pu, qu
      complex*16 temp, phx
      
      complex*16 jl, jlp1, nl, nlp1
      double precision dx, rmt, rnrm, x0
      integer i, ie, ilast, irr, ikap, lll, ic3
      integer ncycle, jri, jnrm, lmax

      lmax = lmaxsc
      if (lmax .gt. lx) lmax = lx
      if (iz .le. 4)   lmax = 2
      if (iz .le. 2)   lmax = 1
      
      do i = 1, nrptx
        vtotc(i) = vtot(i)
        vvalc(i) = vvalgs(i)
      end do
      
c     set imt and jri (use general Loucks grid)
c     rmt is between imt and jri (see function ii(r) in file xx.f)
      jri  = getiat(x0, dx, rmt) + 1
      if (jri .gt. nrptx)  call par_stop('jri .gt. nrptx in phase')
      jnrm = getiat(x0, dx, rnrm) + 1
c     ilast is the last integration point
c     it is larger than jnrm for better interpolations
      ilast = jnrm + 6
      if (ilast .gt. nrptx) ilast = nrptx

      do ie = 1, ne
c       p2 is (complex momentum)**2 referenced to energy dep xc
        p2 = em(ie) - eref
        if (mod(ixc,10) .lt. 5) then
          ncycle = 0
        else
          ncycle = 3
        endif
        ck = sqrt(2*p2+ (p2*alphfs)**2)
        xkmt = rmt * ck

        do lll = 0, lx
          if (lll .gt. lmax) then
             ph(ie, lll+1)  = 0
             goto 10
          endif

c         may want to use ihole=0 for new screening. 
c         don't want ro use it now
c         ihole = 0
          ikap = -1-lll
          irr  = -1
          ic3  = 1
          if (lll .eq. 0) ic3 = 0
          call dfovrg( ncycle, ikap, rmt, ilast, jri, p2, dx,           &
     &                  ri, vtotc, vvalc, dgcn, dpcn, adgc, adpc, xnval,&
     &                  pu, qu, pn, qn, iz, ihole, xion, iunf, irr, ic3)
          call exjlnl(xkmt, lll, jl, nl)
          call exjlnl(xkmt, lll+1, jlp1, nlp1)
          call phamp(rmt, pu, qu, ck, jl,nl,jlp1,nlp1, ikap, phx, temp)
          ph(ie, lll+1) = phx
          
        end do
   10   continue
      end do

      return
      end
c= ./getph.f}
c={./screen.f
c=======================================================================
c     wscrn
c=======================================================================
c inputs: geom.dat => 1st line, ldos.inp => 2nd, scrn.inp => 3rd
c         pot.bin => 4,5th, grids => 6th, atomic and potentials => 6,7th
c output: wscrn
      subroutine screen( nat, nph, iphat, rat,                          &
     &                   toler1, toler2, lmaxph, rdirec, rfms2, lfms2,  &
     &                   maxl, irrh, iend, lfxc, nrptx0,                &
     &                   ihole, xmu, adgc, adpc, iz,                    &
     &                   xion, iunf, xnval,edens, dgc0, dpc0,           &
     &                   ri, dx, x0, rmt, rnrm, em, ne,                 &
     &                   ixc0, iph, vtot, vvalgs, eref, dgcn, dpcn, ph, &
     &                   ilast, vch, wscrn)
c     linear response calculation
c     calculate static screening
      
      implicit none
c={../HEADERS/declare.h
c     const.h
      double precision pi
      integer one, zero
      double precision third
      double precision raddeg
      double precision fa
      double precision bohr, ryd
      double precision hart
      double precision alpinv
      double precision alphfs

c     dim.h
      integer nclusx
      integer nclxtd
      integer nspx
      integer natx
      integer nattx
      integer lx
      integer nphx
      integer ltot
      integer nrptx
      integer nex, lamtot
      integer mtot, ntot
      integer npatx
      integer legtot
      integer novrx
      integer nheadx
      integer MxPole
c= ../HEADERS/declare.h}
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
      external dgetrf, dgetrs
c     function getiat
      integer getiat
c     pot.bin --------------------------------------------------------
      integer ihole, iz, iunf
      double precision xmu, adgc(10,30), adpc(10,30)
      double precision xion, xnval(30), edens(251)
      double precision dgc0(251), dpc0(251)
c     -----------------------------------------------------------------
      double precision vtot(nrptx), vvalgs(nrptx)
      double precision dgcn(nrptx,30), dpcn(nrptx,30)
      complex*16 eref
c
      complex*16 gtrl(0:lx, nex)
      integer ixc0, ispin
c     -----------------------------------------------------------------
c     radial grid
      double precision rmt, rnrm, x0, dx, ri(nrptx)
      integer kinit, linit, jri, jri1, jnrm, ilast, iph, ic3,irr, ikap
      integer i, j, k, m, n, ir1, ir2, ll, lfin, ilp, ie, ncycle, ios
c     work space for xcpot
      complex*16 v(nrptx), vval(nrptx)
c     work space for fovrg
      complex*16 pr(nrptx), qr(nrptx), pn(nrptx), qn(nrptx)
      complex*16 p2, ck, xkmt, xkmtp, xck
      complex*16 pu, qu, dum1, factor, de
      complex*16 xfnorm, xirf, temp, ph0
      complex*16 bessj(8), bessh(8)
      complex*16 jl, jlp1, nl, nlp1

c     fms routine (variables are in single precision!)
      integer nat
      real rfms2, rdirec, toler1, toler2, rpart, aipart, rat(3,natx)
      integer il, ix, im, nsp, minv, iverb, nph, ipp, ill, lfms2
      complex gg(nspx*(lx+1)**2,nspx*(lx+1)**2,0:nphx)
      complex gtr(0:lx, 0:nphx, nex)
      complex xphase(nspx, -lx:lx, 0:nphx)
      complex cks(nspx)
      complex conis
      parameter ( conis = (0.0,1.0) )
      logical lcalc(0:lx)
      integer lmaxph(0:nphx)
      integer inclus, iphat(natx)
c     phase shift      
      complex*16 ph(nex, ltot+1, 0:nphx)
      
c     Linear Response
c     scrn.inp -------------------------------------------------------
      integer maxl, irrh, iend, lfxc, nrptx0
c     -----------------------------------------------------------------
      integer nrx
      parameter ( nrx = 251 )
      double precision percent
c     energy mesh
      complex*16 em(nex)
      integer ne
c     matrix inversion
      character*1 trans
      parameter ( trans = 'N' )
      integer info, nrhs, lda, ldb
      integer ipiv(nrx)
c     work space
c     core hole potential
      double precision vch(nrx)
c     Coulomb Exchange      
      double precision Kmat(nrx, nrx)
c     LDA fxc      
      double precision fxc(nrx)
c     chi_0(r,r')
      complex*16 chi0r(nrx, nrx)
c     chi_0(r,r',e)
      complex*16 chi0re(nrx, nrx)
c     (1 - K Chi0)
      double precision Amat(nrx, nrx)
c     W, the result of calculation
      double precision wscrn(nrx)
c     =================================================================      
      call setkap(ihole, kinit, linit)

c     set imt and jri (use general Loucks grid)
c     rmt is between imt and jri (see function ii(r) in file xx.f)
      rmt   = rmt
      rnrm  = rnrm
      jri   = getiat(x0, dx, rmt) + 1
      jri1  = jri+1
      if ( jri1 .gt. nrptx )  call par_stop('jri .gt. nrptx in phase')
c     nesvi: define jnrm
      jnrm  = getiat(x0, dx, rnrm) + 1
c      ilast is the last integration point
c     it is larger than jnrm for better interpolations
      ilast = jnrm + 6 + iend
      if (ilast .gt. nrx) ilast = nrx
      
c     initialize values to zero
      do i = 1, nrx
        do j = 1, nrx
          chi0r (i,j) = 0.0d0
          chi0re(i,j) = 0.0d0
        end do
      end do
      percent  = 0

      call yprep(iph, nat, inclus, nph, iphat, rfms2, rat, iz, rdirec)
      
      if (inclus .gt. 1) then
        print *, ' Doing FMS for a cluster of ', inclus,                &
     &         ' atoms around iph = ', iph
      end if

c===============================================================
c     start energy cycle (big loop!)
c===============================================================
      do i = 1, nrptx
        v(i)    = vtot(i)
        vval(i) = vvalgs(i)
      end do

      do ie = 1, ne
c       set the method to calculate atomic cross section
c       p2 is (complex momentum)**2 referenced to energy dep xc
        p2   = em(ie) - eref
        ck   = sqrt(2*p2 + (p2*alphfs)**2)
        rpart  = real( dble(ck))
        aipart = real(dimag(ck))
        cks(1) = cmplx(rpart, aipart)
        xkmt = rmt * ck
        !write(99,*) dreal(em(ie))*hart, dimag(em(ie))*hart, ie

        if ( mod(ixc0,10) .lt. 5 ) then
          ncycle = 0
        else
c         fix later . may be ncycle can be less
          ncycle = 3
        endif
        
c       the maximum angular momentum        
        do k = 1, maxl

          ll = k - 1
          ikap = -k

          ic3 = 1
          if (ll .eq. 0) ic3 = 0
          
c=========================================================
c         get ''regular'' solution
c=========================================================
          irr = -1
          
          call dfovrg( ncycle, ikap, rmt, ilast, jri, p2, dx,           &
     &                  ri , v, vval, dgcn, dpcn,  adgc, adpc, xnval,   &
     &                  pu, qu, pr, qr, iz, ihole, xion, iunf, irr, ic3)
          
          call exjlnl(xkmt, ll, jl, nl)
          call exjlnl(xkmt, ll+1, jlp1, nlp1)
          call phamp(rmt,pu,qu,ck,jl,nl,jlp1,nlp1,ikap,ph0,temp)
          
          factor = ck*alphfs
          factor = - factor/(1+sqrt(1+factor**2))
          dum1   = 1 / sqrt(1+factor**2)
          xfnorm = 1 / temp * dum1
          
c         normalization factor
c         dum1 is relativistic correction to normalization
c         normalize regular solution
          do i = 1, ilast
            pr(i) = pr(i)*xfnorm
            qr(i) = qr(i)*xfnorm
          end do
c===========================================================
c         get ''irregular'' solution
c===========================================================
          irr = 1
c         pu, qu are upper and lower components at muffin tin
c         set pu, qu - initial condition for irregular solution
          
          pu =   (nl*cos(ph0) +   jl*sin(ph0)) * rmt * dum1
          qu = (nlp1*cos(ph0) + jlp1*sin(ph0)) * rmt * dum1 * factor
          
          if (irrh .eq. 1) then
            call besjh(xkmt, maxl+1, bessj, bessh)
            pu = bessh(ll+1)* exp(coni*ph0) *rmt * dum1
            qu = bessh(ll+2)* exp(coni*ph0) *rmt * dum1* factor
          end if
          
          call dfovrg( ncycle, ikap, rmt, ilast, jri, p2, dx,           &
     &                  ri, v, vval, dgcn, dpcn, adgc, adpc, xnval,     &
     &                  pu, qu, pn, qn, iz, ihole, xion, iunf, irr, ic3)

c         set N- irregular solution , which is outside
c         N=(nlp1*cos(ph0)+jlp1*sin(ph0))*factor *rmt * dum1
c         N = i*R - H*exp(i*ph0)
          temp = exp(coni*ph0)
c         calculate wronskian (qu)
          qu = 2*alpinv *temp*(pn(jri)*qr(jri)-pr(jri)*qn(jri))
          qu = 1 /qu / ck
c         qu should be close to 1
          do i = 1, ilast
            pn(i) = temp * pn(i)*qu
            qn(i) = temp * qn(i)*qu
          end do
          
c         Use exact solution to continue solutions beyond rmt
          do j = jri, ilast
            xck  = ck * ri(j)
            call exjlnl(xck, ll, jl, nl)
            call exjlnl(xck, ll+1, jlp1, nlp1)
            call besjh(xck, maxl+1, bessj, bessh)
            
            pr(j) =  (jl*cos(ph0) -   nl*sin(ph0))*ri(j)*dum1
            qr(j) =(jlp1*cos(ph0) - nlp1*sin(ph0))*ri(j)*dum1*factor
            pn(j) = bessh(ll+1)*exp(coni*ph0) *ri(j)*dum1
            qn(j) = bessh(ll+2)*exp(coni*ph0) *ri(j)*dum1*factor
          end do
c=============================================================
c         setup prefactors
c=============================================================
          ll = k - 1
          
c         the factor includes averaged over spins
          factor = -1.0d0/(2*(pi**2))*(2*ll+1.0d0)*(2.0d0*ck)**2*dx**2
c=============================================================
c         FMS
c=============================================================
          if (inclus .gt. 1) then
            do ipp = 0, nph
              do ill = -lmaxph(ipp), lmaxph(ipp)
                rpart  = dble( ph(ie, abs(ill)+1, ipp))
                aipart = dimag(ph(ie, abs(ill)+1, ipp))
                xphase(1, ill, ipp) = cmplx(rpart, aipart)
              end do
            end do
          
            iverb = 0
            nsp = 1
            ispin = 0
            do i = 0, lx
              lcalc(i) = .true.
            end do
            minv = 0

            do i = 0, lx
              do j = 0, nphx
                gtr(i,j,ie) = 0.0d0
              end do
            end do

            call fms(lfms2, nsp, ispin, inclus, nph, cks,lmaxph,        &
     &              xphase,ie,iverb,minv,rdirec,toler1,toler2,lcalc, gg)
            
            do il = 0, lmaxph(iph)
                ix = il**2
                do im = 1, 2*il+1
                  gtr(il,iph,ie) = gtr(il,iph,ie) + gg(ix+im,ix+im,iph)
                end do
     
                gtrl(il,ie) = (        dble(real( gtr(il,iph,ie)))      &
     &                          + coni*dble(aimag(gtr(il,iph,ie)))     )&
     &                        * exp(2*coni*ph(ie,il+1,iph))/(2*il+1.0d0)
            end do
          end if
c=============================================================
c         construct response function
c=============================================================         
c         Construct G_l(r,r')G_l(r',r) for W0 calculation, G(r,r') = R(r<)*H(r>)
          do m = 1, ilast
            do n = m, ilast
              chi0re(m,n) = chi0re(m,n) + factor * ri(m)*ri(n)          &
     &                             * pr(m)*pr(m) * pn(n)*pn(n)
            end do

          end do

          if (inclus .gt. 1) then
            do m = 1, jnrm
              do n = m, jnrm
                chi0re(m,n) = chi0re(m,n)+factor * ri(m)*ri(n) * (      &
     &              2 * gtrl(ll,ie)* pr(m)*pr(m) * pr(n)*pn(n)          &
     &            + gtrl(ll,ie)**2 * pr(m)*pr(m) * pr(n)*pr(n)   )
              end do

            end do
          end if
       
        end do

c===============================================================
c       integration over energy
c===============================================================
        if (ie .eq. 1) then
          de = (em(ie+1) - em( ie ))/2.0d0
        else if (ie .eq. ne) then
          de = (em( ie ) - em(ie-1))/2.0d0
        else
          de = (em(ie+1) - em(ie-1))/2.0d0
        end if

        do ir1 = 1, ilast
          do i = ir1, ilast
            chi0r( ir1,i) = chi0r(ir1,i) + chi0re(ir1,i)*de
            chi0re(ir1,i) = (0.0d0, 0.0d0)
          end do
        end do
        
c       print the current progress to the wscrn
        if ( 100.0*ie/ne .gt. percent) then
          print '(a,f3.0,a)', '  ', percent, '% of energy integration'
          percent = percent + 10.0
        end if
      end do
      
      print *, '100.% -- end of energy integration --'
c=================================================================
c     end of energy cycle
c=================================================================
      
c=================================================================
c     setupt K-matrix used for calculating screening effect
c=================================================================
      do m = 1, ilast
        do n = m, ilast
          Kmat(m,n) =  4.0d0*pi * 1 * 1/ri(n)
        end do

      end do
c=================================================================
c     LDA Exchange Correlation (actually this probably doesn't work)
c=================================================================  
      if (lfxc .gt. 0) then
        call ldafxc(ilast, ri, edens, lfxc, fxc)
        print *, 'Local exchange on'
        do i = 1, ilast
          Kmat(i,i) = Kmat(i,i) + 4.0d0*pi * fxc(i)
        end do
      end if
c=================================================================
c     get response function
c=================================================================
      print *, 'Prepare response function'
      
      do ir1 = 2, ilast
        do i = 1, ir1-1
          chi0r(ir1,i) = chi0r(i,ir1)
          Kmat(ir1,i) = Kmat(i,ir1)
        end do
      end do

c     here performs v_ch(r) = \int (dgc0'^2 + dpc0'^2)*1/r> dr'
      do i = 1, ilast
        vch(i)  = (dpc0(i)**2 + dgc0(i)**2) * dx * ri(i)
      end do
c     print *, 'Check normalization', sum(vch)
      
      do i = 1, ilast
        wscrn(i) = 0.0d0
        do j = 1, i
          wscrn(i) = wscrn(i) + vch(j)
        end do
        wscrn(i) = wscrn(i)/ri(i)
        do j = i + 1, ilast
          wscrn(i) = wscrn(i) + vch(j)/ri(j)
        end do     
      end do
      
      do i = 1, ilast
        vch(i) = wscrn(i)
      end do

c=================================================================
c     start matrix inversion
c=================================================================
c     the matrix to be inverted ( .I. N = N x N Identity Matrix )
      print *, 'Start matrix inversion'
c     Amat = Amat - matmul(Kmat, Chi0)
      do i = 1, ilast
        do j = 1, ilast
          if (i .eq. j) then
            Amat(i,j) = 1.0d0
          else
            Amat(i,j) = 0.0d0
          end if
          do k = 1, nrx
            Amat(i,j) = Amat(i,j)-Kmat(i,k)*dimag(chi0r(k,j))
          end do
        end do
      end do
      
      print *, 'Compute (1 - K Chi0)^-1 v_ch'
      
      nrhs = 1
c     size(Amat,1)
      lda  = nrx
c     size(wscrn)
      ldb  = nrx
c     actual size of nxn matrix
      n    = ilast
      
      call dgetrf(n, n, Amat, lda, ipiv, info)
      
      if (info .eq. 0) then
        call dgetrs(trans, n, nrhs, Amat, lda, ipiv, wscrn, ldb, info)
        
        if (info .ne. 0) print *, info
      else
        print *, 'The factor U is singular'
      end if
c=================================================================
c     end of matrix inversion
c=================================================================
      
      return
c      end subroutine screen
      end
c=======================================================================
c     END wscrn
c=======================================================================
c= ./screen.f}
      subroutine fmsie( iph0, nph, lipotx, ie, em, eref, ph, iz,
     1                 rfms, lfms, nat, iphat, rath, gtr)

c     full multiple scattering code for single energy point
c     written by a.ankudinov 06.1997 using earlier written subroutines
c     coded by b.ravel
c     modified by a.ankudinov 2001 for new matrix inversion algorithms
c     Feb. 2002, a.ankudinov: fixed logic for MPI calculations
c       lfms=0  - extended system calculataions (e.g. crystal)
c       lfms=1  - small system calculations (e.g. molecule)
c       lfms=2  - same as 1 for MPI run (forces call yprep)

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

c     input
      dimension iphat(natx), rath(3,natx)
      real rat(3,natx), rfms, rdirec, toler1, toler2
      real rpart,aipart
      integer nph
      dimension iz(0:nphx)
      complex*16 ph(lx+1, 0:nphx)

c     work space
      integer iph0
      complex*16 em, eref
      character*512 slog
      logical lcalc
      dimension lcalc(0:lx)
c     fms staff
      integer lipotx(0:nphx)
      complex gg(nspx*(lx+1)**2, nspx*(lx+1)**2, 0:nphx)
      complex gtr(0:lx, 0:nphx)
      complex xphase(nspx, -lx:lx, 0:nphx), ck(nspx)
      complex*16 dck
      complex conis
      parameter (conis = (0,1))
      real  temper, thetax, sig2
      save

      if (rfms .le. 0.0) goto 900

c     set default (LU) inv method
      minv = 0
      rdirec = 2*rfms
      toler1 = 0.e0
      toler2 = 0.e0

      do 30 iat=1,nat
      do 30 j=1,3
   30 rat(j,iat) = real (rath(j,iat))

c     transform to single precision
      temper =0.0e0
      thetax =0.0e0
      sig2  = 0.0e0

c      it will be nice to call yprep once for all energy points,
c      fix later, and now call it every time
      if (ie.eq.1 .or. lfms.eq.0 .or. lfms.eq.2) 
     1  call yprep(iph0, nat, inclus, nph, iphat, rfms, rat,
     2     iz, rdirec )

      if (inclus.gt.1) then

cc     call fms for a cluster around central atom
       if (ie.eq.1) then
          write (slog,35) inclus, iph0
  35      format ('        Doing FMS for a cluster of ',i3,
     1    ' atoms around iph = ',i2)
          call wlog (slog)
       endif

       dck=sqrt(2*(em-eref))
       rpart  = dble(dck)
       aipart = real(dimag(dck))
       ck(1) = cmplx(rpart, aipart)
       do 1020 ipp = 0,nph
         do 1010 ill = -lipotx(ipp), lipotx(ipp)
           rpart  = dble( ph( 1+abs(ill), ipp))
           aipart = dimag(ph( 1+abs(ill), ipp)) 
           xphase(1, ill, ipp) = cmplx(rpart, aipart)
 1010    continue
 1020  continue
       iverb=0
       if (ie.eq.1) iverb = 1
       nsp = 1
       ispin = 0
       do 1011 ill = 0,lx
 1011  lcalc(ill) = .true.
       call fms(lfms, nsp, ispin, inclus, nph, ck, lipotx, xphase, ie,
     1  iverb, minv, rdirec, toler1, toler2, lcalc, gg)

c      make ck= i, since coni is c*16
       do 1030 ip=0,nph
         if (lfms.ne.0 .or. ip.eq.iph0) then
           do 1040 il=0,lipotx(ip)
             ix = il**2
             do 1050 im=1,2*il+1
               gtr(il, ip) = gtr(il, ip) + gg(ix+im,ix+im,ip)
 1050        continue
             gtr(il,ip)= gtr(il,ip)*
     1            exp(2*conis*xphase(1,il,ip))/(2*il+1)
 1040      continue
         endif
 1030  continue
      endif

 900  continue
      return
      end
      subroutine fms(lfms, nsp, ispin, inclus, npot, ck, lipotx, xphase,
     1   ik, iverb, minv, rdirec, toler1, toler2, lcalc, gg)

      implicit real (a-h,o-z)
      implicit integer (i-n)
c--------------------------------------------------------------------
c  compute full multiple scattering within some cluster at some
c  energy 
c  This uses the LU decomposition package from LAPACK.  Driver
c  routines: cgetrf (decomposition), cgetrs (backsubstitution)
c  coded by b.ravel
c  modified by a.l.ankudinov to include spin and SO interactions
c  feb 2000
c
c  dim.h and xparam.h must be included in the calling routine
c
c  most of the information needed by this package is set into common
c  blocks the companion package xprep.  In that package, the lists of
c  atomic coordinates and potential indeces are organized so that the
c  first npot+1 entries are examples of each of the unique potentials.
c  Consequently, only the upper left hand corner of the FMS matrix
c  need be recomposed to get the set of submatrices necessary to
c  compute chi for each type of atom in the cluster.
c
c  See subroutine fmstot.f for an example of decoding the output of this
c  subroutine. The third index of gg refers to the unique potential with
c  element 0 being the absorbing atom.  
c  The first two indeces are related to the |lms> state by the
c  formula:
c       nsp=1, no spin indeces
c       lm  = ( l**2 + 1 ) + ( l + m )
c            thus {1..(lx+1)^2} ==>
c            {0,0 1,-1 1,0 1,1 2,-2 2,-1 2,0 2,1 2,2 ...
c                   lx,lx-1 lx,lx}
c       nsp=2, with spin indeces
c       lms  = 2*( l**2 + 1 ) + 2*( l + m ) + (s-1/2)
c            thus {1...2*(lx+1)^2} ==>
c            {0, 0,-1/2  0. 0,1/2
c             1,-1,-1/2  1,-1,1/2  1,0,-1/2  1,0,1/2  1,1,-1/2 1,1,1/2
c             2,-2,-1/2  2,-2,1/2  2,-1,-1/2 2,-1,1/2 ...    lx,lx,1/2}
c
c  The calling protocol for xpreppack and fmspack is;
c          include 'dim.h'
c          include 'xparam.h'
c          ...
c          call xprep(nat, inclus, npot, iphat, rmax, rat,
c     $            xnrm, izx, temper, thetad)
c          energy loop {
c             ...
c             call fms(nsp, inclus, npot, ck, lipotx, xphase,
c                      ik, iverb, gg)
c             ... }
c
c  fmspack contains the following routines:
c    fms.f:     main routine of fmspack
c    kets.f:    compute all state kets for current energy
c    xclmz.f:   compute hankle-like polynomials for current energy
c    xgllm.f:   compute z-axis propagators for current energy
c    cgetrf.f:  LU decomposition of matrix
c    cgetrs.f:  backsubstitution of LU decomposed matrix
c    lu_misc.f: various routines called by LU package
c
c---------------------------------------------------------------------
c  input
c    nsp:    1) no spin indeces 2) with spin indeces
c    inclus: number of atoms in cluster
c    npot:   number of unique potentials in cluster
c    ck:     complex momentum of current energy point
c    lipotx: (0:nphasx) max l for each unique potential
c    xphase: (0:lx, 0:nphasx) single complex array of partial wave
c            phase shifts for each unique potential
c    ik:     current energy grid index, used for run-time messages
c    iverb:  do nothing when iverb <= 0
c            1  => write a message about grid point and matrix size
c
c  passed in common from xprep package (xstruc.h)
c    xrat:   (3,nclusx) array of coordinates with first npot+1
c            elements each a unique potential
c    xphi:   (nclusx, nclusx) angles between z axis and vectors
c            connecting the atoms in the cluster
c    iphx:   (nclusx) potential index of each atom in the cluster
c    drix:   huge matrix containing all rotation matrix elements
c            needed for computation of free electron propagators
c    xnlm:   matrix of legendre polynomial normalization factors
c    xpsile: matrix containing wave functions for hybridization
c            calculation
c    sigsqr: (nclusx,nclusx) matrix of pair-wise mean square
c            displacements about interatomic distances.  Currently only
c            calculated by the correlated debye model.
c
c  output
c    gg:  (nsp*lx**2, nsp*lx**2, 0:nphasx) submatrix spanning the entire
c          angular momentum basis for each unique potential

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
c={xparam.h
c nphasx MUST be the same as nphx, the maximum number of unique
c        potentials
c natxx MUST be the same as natx, the maximum number of atoms in the
c       extendeed cluster
c nexx MUST be the same as nex, the maximum number of energy points
      parameter (nphasx=nphx)
      parameter (natxx=natx)
      parameter (nexx=nex)
      parameter (istatx=(lx+1)**2*nclusx*nspx)
      parameter (nkmin=1)
c     parameter (nkmin=-9)
c= xparam.h}
c====================================================================
c  This common block contains the structural information about the
c  cluster to be used for the full multiple scattering calculation
c  xphi:  matrix of angles between the z axis and pairs of atoms
c  xrat:  xyz coordinates of the atoms in the cluster, the first
c         npot+1 entries are examples of each unique potential
c  iphx:  potential indeces of each atom in the cluster, ordered like
c         xrat
      common /xstruc/ xphi(nclusx,nclusx), xrat(3,nclusx),
     $            iphx(nclusx)
      save /xstruc/
c********************************************************************
c**** rotation matrices for entire cluster
c
      complex drix
      common /rotx/ drix(-lx:lx, -lx:lx, 0:lx, 0:1, nclusx, nclusx)
      save /rotx/
c********************************************************************
c common blocks for saving rotation matrices between xanes and rotxan
       integer    jsavx, jsav, jbmagk
       parameter (jsavx = 150, roteps = 1.e-12,jbmagk=-9999)
       dimension drisav(-lx:lx,-lx:lx,0:lx,jsavx), betsav(jsavx)
       integer   ldsav(jsavx), mdsav(jsavx)
       common /rotsav/  drisav, betsav, ldsav, mdsav, jsav
       save  /rotsav/

c********************************************************************
c**** legendre polynomial normalization constants
c
      common /lnlm/ xnlm(0:lx,0:lx)
      save   /lnlm/

c********************************************************************
c**** legendre polynomial normalization constants
c
      common /xdwf/ sigsqr(nclusx,nclusx)
      save   /xdwf/
c********************************************************************
c**** save Clebsch-Gordon coefficients: <LS|J>
      dimension t3jp(0:lx, -lx:lx, 2), t3jm(0:lx, -lx:lx, 2)
      common /t3j/ t3jp, t3jm
      save   /t3j/

c********************************************************************
      parameter (pi = 3.14159 26535 89793 23846 26433e0)
      parameter (bohr = 0.529 177 249e0)
      parameter (one = 1, zero = 0)
      complex coni
      parameter (coni = (0,1))
      complex   term, prefac, gllmz, ck(nspx)
      complex   clm(lx+2, 2*lx+3), xclm(0:lx, 0:lx, nclusx, nclusx,nspx)
      complex   xrho( nclusx, nclusx, nspx)
      integer   lipotx(0:nphasx)

c********************************************************************
c**** array of state kets at current energy
      common /stkets/ lrstat(4, istatx), istate
      save   /stkets/
      complex   xphase(nspx, -lx:lx, 0:nphasx)
      complex   tmatrx(nspx, istatx)
c     big work matrices
      complex   g0(istatx,istatx), g0t(istatx,istatx)
      logical lcalc
      dimension lcalc(0:lx)
c     return matrix containing info about each unique potential
      complex   gg(nspx*(lx+1)**2, nspx*(lx+1)**2, 0:nphasx)

      integer i0(0:nphx)
      character*3  cerr, dec
      character*13 trans
      character*75 messg

 400  format(i4)

      do 10 i=0,nphx
        if (lipotx(i).le.0)  lipotx(i) = lx
        if (lipotx(i).gt.lx) lipotx(i) = lx
        i0(i) = -1
 10   continue
c     initialize gg to zero
      do 20 i = 0, nphasx
        do 18 j = 1, nspx*(lx+1)**2
          do 16 k = 1, nspx*(lx+1)**2
            gg( k, j, i) = cmplx( zero, zero)
 16       continue
 18     continue
 20   continue

      if (lfms.eq.0) then
        ipi = iphx(1)
        ipf = iphx(1)
      else
        ipi = 0
        ipf = npot
      endif
c --- get basis kets; output array 'lrstat' passed through common
      call getkts(nsp, inclus, npot, iphx, lipotx, i0)

c --- sanity check for i0(ip)
      do 30 ip = ipi, ipf
        if (i0(ip) .lt. 0) then
          call wlog (' Cannot find all representative atoms')
          call wlog (' Increase FMS radius and rerun.')
          call par_stop(' In subroutine FMS')
        endif
  30  continue

c --- runtime message if requested
      if (iverb.gt.0 .and. minv.eq.0) then
         dec = 'LUD'
         write(messg, 4010)this_process,dec, ik, istate
 4010    format('  ',i3,'   FMS matrix (', a, ') at point ', i3,
     $               ', number of state kets =', i4)
         call wlog(messg)
      endif

c --- get all c_lm(z) values for this energy, i,j sum over all atom
c     pairs xrho and xclm are symmetric in ij
c+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
c  nota bene, in the code for setting the clmz, the indexing starts
c  at 1 rather than 0.  To my mind, that is confusing, so here I
c  reindex when I copy from clm to xclm.  See the note about this in
c  subroutine xclmz
c+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
      lplus1 = lx+1
      mplus1 = lx+1
      do 140  i=1,inclus
        do 130 j=1,i

c ------- get and store rho for this pair of atoms   bohr units on
c         r and ck
          r   = zero
          do 100 ix=1,3
            r = r + (xrat(ix,i) - xrat(ix,j))**2
 100      continue
          r   = sqrt(r)

          do 125 isp = 1,nsp
             xrho(i,j,isp) = ck(isp) * r
             xrho(j,i,isp) = xrho(i,j,isp)

c ------- store the c_lm(z) for all the rhos at this energy
c            xclm(i,j) = xclm(j,i) by symmetry
             if (i.ne.j) call xclmz(lplus1,mplus1,xrho(i,j,isp),clm)
             do 120 ll = 0,lx
               do 110 mm = 0,lx
                 if (i.eq.j) then
                     xclm(mm,ll,j,i,isp) = cmplx(zero,zero)
                 else
                     xclm(mm,ll,j,i,isp) = clm(ll+1,mm+1)
                     xclm(mm,ll,i,j,isp) = clm(ll+1,mm+1)
                 endif
 110           continue
 120         continue
 125      continue

 130    continue
 140  continue

c --- fill the G0 and T matrices for this energy
      rdir2 = rdirec**2
      do 220 ist1=1,istate
        iat1 = lrstat(1, ist1)
        l1   = lrstat(2, ist1)
        m1   = lrstat(3, ist1)
        isp1 = lrstat(4, ist1)

        do 210 ist2=1,istate
          iat2 = lrstat(1, ist2)
          l2   = lrstat(2, ist2)
          m2   = lrstat(3, ist2)
          isp2 = lrstat(4, ist2)

          rr = (xrat(1,iat1)-xrat(1,iat2))**2 +
     1    (xrat(2,iat1)-xrat(2,iat2))**2 +(xrat(3,iat1)-xrat(3,iat2))**2

c                               equation 9 in Rehr, Albers
c                               <LR| G |L'R'>

          if (iat1.eq.iat2) then
c             same atom: G=0, calculate T-matrix 
              g0(ist1,ist2)     = cmplx(zero,zero)
c             notice that T is tri-diagonal, due to conservation of
c             total momentum.(will be broken by nonspherical potential)
c             --- potential index for this atom
              iph = iphx(iat1)
            if (nsp.eq.1.and.ispin.eq.0) then
              if (ist1.eq.ist2) tmatrx(1, ist1) =
     $                    ( exp(2*coni*xphase(isp1,l1,iph)) - one )
     $                    / (2*coni) 
            else
              if (ist1.eq.ist2) then
c                set spin index for t3jm and t3jp
                 is = isp1
                 if (nsp.eq.1) then
c                  special case
                   is = 1
                   if (ispin.gt.0) is = 2
                 endif

c                diagonal matrix element
                 tmatrx(1, ist1) =
     $                    ( exp(2*coni*xphase(isp1,l1,iph)) - one )
     $                    / (2*coni) * t3jm (l1, m1, is)**2  + 
     $                    ( exp(2*coni*xphase(isp1,-l1,iph)) - one )
     $                    / (2*coni) * t3jp (l1, m1, is)**2 
              elseif (nsp.eq.2.and.l1.eq.l2.and.m1+isp1.eq.m2+isp2) then
c                same orb. mom. and total momentum projections conserved
c                calculate off-diagonal T-matrix element
c                tmatrx(2, ist1) = here only if nspx equal to 2
                 tmatrx(nsp, ist1) =
     $             ( exp(2*coni*xphase(isp1, l1,iph)) - one +
     $               exp(2*coni*xphase(isp2, l1,iph)) - one ) / (4*coni) 
     1             * t3jm (l1, m1, isp1) * t3jm (l1, m2, isp2)  + 
     $             ( exp(2*coni*xphase(isp1,-l1,iph)) - one +
     $               exp(2*coni*xphase(isp2,-l1,iph)) - one ) / (4*coni) 
     1             * t3jp (l1, m1, isp1) * t3jp (l1, m2, isp2)
              endif
            endif
          elseif (isp1.eq.isp2 .and. rr.le.rdir2) then
c           different atoms, same spin: T=0, calculate G
            g0(ist1,ist2) = cmplx(zero,zero)
            do 200 mu=-l1,l1
c             --- third arg in drix: 0==>beta, 1==>-beta
              muabs = abs(mu)
              call xgllm(muabs, ist1, ist2, lrstat,
     1                   xclm(0,0,1,1,isp1), gllmz )
              g0(ist1,ist2) = g0(ist1,ist2) +
     2             drix(mu,m1,l1,1,iat2,iat1) *  gllmz *
     3             drix(m2,mu,l2,0,iat2,iat1)
 200        continue
            prefac = exp(coni*xrho(iat1,iat2,isp1)) /
     $                  xrho(iat1,iat2,isp1)
c           use correlated debye model, sigsqr is in AA^2
            prefac = prefac * exp(-1 * sigsqr(iat1,iat2) *
     $                  ck(isp1)**2 / bohr**2)
            g0(ist1,ist2) = prefac * g0(ist1,ist2)
          else
c           different atoms, different spins:T=G=0
            g0(ist1,ist2) = cmplx(zero,zero)
          endif

c -----   end of loops over states
 210    continue
 220  continue

      if (minv.eq.0) then
         call gglu ( nsp, i0, ipi, ipf, lipotx, g0, tmatrx, g0t, gg)
      elseif (minv.eq.1) then
         dec = 'VdV'
         call ggbi ( nsp, i0, ipi, ipf, lipotx, g0, tmatrx, g0t, gg,
     1              toler1, toler2, lcalc, msord)
      elseif (minv.eq.2) then
         dec = 'LLU'
         call ggrm ( nsp, i0, ipi, ipf, lipotx, g0, tmatrx, g0t, gg,
     1              toler1, toler2, lcalc, msord)
      elseif (minv.eq.3) then
         dec = 'GMS'
         call gggm ( nsp, i0, ipi, ipf, lipotx, g0, tmatrx, g0t, gg,
     1              toler1, toler2, lcalc, msord)
      else
         dec = 'TF'
         call ggtf ( nsp, i0, ipi, ipf, lipotx, g0, tmatrx, g0t, gg,
     1              toler1, toler2, lcalc, msord)
      endif
      if (minv.ne.0) then
         write(messg, 410)this_process,dec, ik, istate, msord
 410     format('  ',i3,'. Iterative FMS (', a, ') at point ', i3,
     $               '; matrix size =', i4,'; MS order =',i5)
         call wlog(messg)
      endif

      return
      end
c--------------------------------------------------------------------
      subroutine getkts(nsp, nat, npot, iphx, lipotx, i0)

      implicit real (a-h,o-z)
      implicit integer (i-n)
c--------------------------------------------------------------------
c  construct state kets |iat,l,m> at this energy
c--------------------------------------------------------------------
c  input
c    nat:    number of atoms in cluster
c    npot:   number of unique potentials
c    iphx:   (nclusx) potential index of each atom in the cluster
c    lipotx: (nphasx) maximum angular momentum to consider for each
c            ipot
c  output
c   (istate: number of states  ---  passed in kets.h)
c    i0:     index shift for each potential representative
c   (lrstat: (4, istatx) state kets |iat,l,m> --- passed in kets.h)
c--------------------------------------------------------------------
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
c={xparam.h
c nphasx MUST be the same as nphx, the maximum number of unique
c        potentials
c natxx MUST be the same as natx, the maximum number of atoms in the
c       extendeed cluster
c nexx MUST be the same as nex, the maximum number of energy points
      parameter (nphasx=nphx)
      parameter (natxx=natx)
      parameter (nexx=nex)
      parameter (istatx=(lx+1)**2*nclusx*nspx)
      parameter (nkmin=1)
c     parameter (nkmin=-9)
c= xparam.h}
c********************************************************************
c**** array of state kets at current energy
      common /stkets/ lrstat(4, istatx), istate
      save   /stkets/
      integer   lipotx(0:nphasx), iphx(nclusx), i0(0:nphx)

      istate = 0
      do 120 iat=1,nat
        ip = iphx(iat)
c       i0(ip) - index for the ip-representative atom
c       need for simple find of states for ip-representative.
        if (i0(ip).lt.0) i0(ip) = istate
        lim = min(lx, lipotx(ip))
        do 110 l=0,lim
          do 100 m = -l, l
          do 100 isp = 1, nsp
            istate = istate + 1
            if (istate.gt.istatx) then
                call wlog('Exceeded maximum number of LR states.'//
     $                      '  Stopping')
                call par_stop('GETKTS-1')
            endif
            lrstat(1,istate) = iat
            lrstat(2,istate) = l
            lrstat(3,istate) = m
            lrstat(4,istate) = isp
 100      continue
 110    continue
 120  continue

      return
c end subroutine kets
      end
c    ----------------------------------------------------------------
      subroutine xclmz(lmaxp1,mmaxp1,rho,clm)
      implicit real(a-h,o-z)

c     calculates energy dependent factors needed in subroutine gllm
c     c(il,im) = c_l^(m)z**m/m!=c_lm             by recursion
c     c_l+1,m  = c_l-1,m-(2l+1)z(c_l,m-c_l,m-1)  l ne m
c     c_m,m    = (-z)**m (2m)!/(2**m m!)         with z=1/i rho
c
c  input:
c    lmaxp1, mmaxp1:  largest angular momentum under consideration + 1
c    rho:  distance between atoms * complex momentum at this energy
c          point
c  output:
c    clm(lx+1,lx+1):  Hankle-like polynomials from RA

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
c={xparam.h
c nphasx MUST be the same as nphx, the maximum number of unique
c        potentials
c natxx MUST be the same as natx, the maximum number of atoms in the
c       extendeed cluster
c nexx MUST be the same as nex, the maximum number of energy points
      parameter (nphasx=nphx)
      parameter (natxx=natx)
      parameter (nexx=nex)
      parameter (istatx=(lx+1)**2*nclusx*nspx)
      parameter (nkmin=1)
c     parameter (nkmin=-9)
c= xparam.h}
      parameter (one = 1, zero = 0)
      complex coni
      parameter (coni = (0,1))
      parameter (ltotb=lx+1,mtotb=ltotb,ntotb=ltotb,mntot=mtotb+ntotb)
      complex z, cmm, clm(ltotb+1,mntot+1), rho

      cmm  = cmplx(one, zero)
      z    = (-coni)/rho

      clm(1,1) = cmplx(one,zero)
      clm(2,1) = clm(1,1) - z

      lmax = lmaxp1-1
      do 20 il=2,lmax
        clm(il+1,1) = clm(il-1,1) - (z * (2*il-1) * clm(il,1))
 20   continue
c+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
c  nota bene:  the 2l-1 factor above is correct, even though in Rehr,
c  Albers equation 4 appears with a 2l+1.  The reason has to do with
c  the indexing.  in RA the subscripts on the c's start at 0.  In this
c  piece of code, the subscripts start at 1.  If you sub l-1 for
c      l, 2l+1 --> 2l-1
c+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+


      mmxp1 = min(lmaxp1, mmaxp1)
      do 40 im=2,mmxp1
        m    = im-1
        imp1 = im+1
        cmm  = (-cmm) * (2*m-1) * z
        clm(im,im)   = cmm
        clm(imp1,im) = cmm * (2*m+1) * (1-im*z)
        do 30 il=imp1,lmax
          clm(il+1,im) = clm(il-1,im) - (2*il-1) * z *
     $                               ( clm(il,im)+clm(il,m) )
c           l = il-1
c           clm(il+1,im) = clm(l,im) - (2*l+1) * z *
c      $                               ( clm(il,im)+clm(il,m) )
 30     continue
 40   continue

      return
c  end subroutine xclmz
      end
      subroutine xgllm(mu, ist1, ist2, lrstat, xclm, gllmz)
c--------------------------------------------------------------------
c  this calculates equations 11,12 from Rehr, Albers PRB v.41,#12,
c  p.8139,  the output is the G term in equation 9 from that paper
c
c  input:
c    mu:         abs val of magnetic state in sum in eqn 11 RA, mu>=0
c    ist1, ist2: state indices of mat. elem., first index of lrstat
c    lrstat:     (4,istatx,nkmin:nex) array of LR states
c    xclm:       (0:lx,0:lx,nclusx,nclusx) array of c_lm(z) for
c                present energy value
c  output:
c    gllmz:      g_ll'^|m|(z), for present state & energy, eqn 11 RA
c--------------------------------------------------------------------
c  this requires that N_lm normalization factors and c_lm(z)
c  polynomials have already been calculated.
c--------------------------------------------------------------------
      implicit real(a-h,o-z)
      implicit integer (i-n)

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
c={xparam.h
c nphasx MUST be the same as nphx, the maximum number of unique
c        potentials
c natxx MUST be the same as natx, the maximum number of atoms in the
c       extendeed cluster
c nexx MUST be the same as nex, the maximum number of energy points
      parameter (nphasx=nphx)
      parameter (natxx=natx)
      parameter (nexx=nex)
      parameter (istatx=(lx+1)**2*nclusx*nspx)
      parameter (nkmin=1)
c     parameter (nkmin=-9)
c= xparam.h}
c====================================================================
c  This header file contains the structural information about the
c  cluster to be used for the full multiple scattering calculation

      common /xstruc/ xphi(nclusx,nclusx), xrat(3,nclusx),
     $            iphx(nclusx)
      save /xstruc/

c  xphi:  matrix of angles between the z axis and pairs of atoms
c  xrat:  xyz coordinates of the atoms in the cluster, the first
c         npot+1 entries are examples of each unique potential
c  iphx:  potential indeces of each atom in the cluster, ordered like
c         xrat
c********************************************************************
c**** rotation matrices for entire cluster
c
      complex drix
      common /rotx/ drix(-lx:lx, -lx:lx, 0:lx, 0:1, nclusx, nclusx)
      save /rotx/
c********************************************************************
c common blocks for saving rotation matrices between xanes and rotxan
       integer    jsavx, jsav, jbmagk
       parameter (jsavx = 150, roteps = 1.e-12,jbmagk=-9999)
       dimension drisav(-lx:lx,-lx:lx,0:lx,jsavx), betsav(jsavx)
       integer   ldsav(jsavx), mdsav(jsavx)
       common /rotsav/  drisav, betsav, ldsav, mdsav, jsav
       save  /rotsav/

c********************************************************************
c**** legendre polynomial normalization constants
c
      common /lnlm/ xnlm(0:lx,0:lx)
      save   /lnlm/

c********************************************************************
c**** legendre polynomial normalization constants
c
      common /xdwf/ sigsqr(nclusx,nclusx)
      save   /xdwf/

c  end of xstruc.h
c********************************************************************

      parameter (zero=0.e0)
      integer    lrstat(4, istatx)
      complex xclm(0:lx, 0:lx, nclusx, nclusx), sum, gllmz
      complex gam, gamtl

      iat1     = lrstat(1, ist1)
      l1       = lrstat(2, ist1)
      iat2     = lrstat(1, ist2)
      l2       = lrstat(2, ist2)
      numax    = min(l1, l2-mu)

      sum      = cmplx(zero, zero)
      do 100 nu=0,numax
        mn    = mu+nu

c       bug for xnlm with nspx=2
        gamtl = (2*l1+1) * xclm(nu, l1, iat2, iat1) / xnlm(mu, l1)
        gam   = (-1)**mu * xclm(mn, l2, iat2, iat1) * xnlm(mu, l2)

        sum   = sum + gamtl * gam
 100  continue

      gllmz = sum

      return
c  end subroutine gllm
      end

c====================================================================
      subroutine gglu( nsp, i0, ipi, ipf, lipotx, g0, tmatrx, g0t, gg)

      implicit real (a-h,o-z)
      implicit integer (i-n)
c  output
c    gg:  (nsp*lx**2, nsp*lx**2, 0:nphasx) submatrix spanning the entire
c          angular momentum basis for each unique potential

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
c={xparam.h
c nphasx MUST be the same as nphx, the maximum number of unique
c        potentials
c natxx MUST be the same as natx, the maximum number of atoms in the
c       extendeed cluster
c nexx MUST be the same as nex, the maximum number of energy points
      parameter (nphasx=nphx)
      parameter (natxx=natx)
      parameter (nexx=nex)
      parameter (istatx=(lx+1)**2*nclusx*nspx)
      parameter (nkmin=1)
c     parameter (nkmin=-9)
c= xparam.h}
      integer  i0 (0:nphx),  lipotx(0:nphx)

      parameter (one = 1, zero = 0)
      integer   ipiv(istatx)

c**** array of state kets at current energy
      common /stkets/ lrstat(4, istatx), istate

      complex   tmatrx(nspx, istatx)
c     big work matrices
      complex   g0( istatx, istatx), g0t( istatx, istatx)
      complex   g0s( istatx, nspx*(lx+1)**2)
c     return matrix containing info about each unique potential
      complex   gg(nspx*(lx+1)**2, nspx*(lx+1)**2, 0:nphasx)

      character*3  cerr
      character*13 trans

 400  format(i4)


c -------------------- LU gg
c     multiply T and G0 matrices together, construct g0t = 1 - G0*T
c     notice that the signs below for g0t ARE correct since 1 is the
c     unit matrix
c     since t is tri-diagonal, this product can be computed in n^2 time
c     also fill up some work matrices for use in eigenvalue and
c     determinant calculations and elsewhere
      do 320 icol = 1,istate
        do 310 irow = 1,istate
c         T diagonal contribution
          g0t(irow, icol) = - g0(irow, icol) * tmatrx(1, icol)
c         T off-diagonal contribution
          l1   = lrstat(2,icol)
          m1   = lrstat(3,icol)
          isp1 = lrstat(4,icol)
          m2 = m1+isp1
          if (nsp.eq.2 .and. m2.gt.-l1+1 .and. m2.lt.l1+2) then
             ist2 = icol + (-1)**isp1
             g0t(irow, icol) = g0t(irow, icol)
     1                   - g0(irow, ist2) * tmatrx(nsp, icol)
          endif
 310    continue

        g0t(icol, icol) = g0t(icol, icol) + one

 320  continue

c --- invert matrix by LU decomposition
c     call cgetrf from lapack.  this performs an LU decomposition on
c     the matrix g0t = 1 - g0*T
      call cgetrf( istate, istate, g0t, istatx, ipiv, info )
      if (info.lt.0) then
          call wlog('    *** Error in cgetrf when computing G')
          write(cerr,400)abs(info)
          call wlog('        Argument #'//cerr//
     $                ' had an illegal value.')
      elseif (info.gt.0) then
          call wlog('    *** Error in cgetrf when computing G')
          write(cerr,400)info
          call wlog('        g0t('//cerr// ','//cerr//
     $                ') is exactly 0 -- '//
     $                'this matrix cannot be decomposed.')
      endif

c     now we want g_c = (g0t)^-1 * g0.  Rather than calculating
c     the inverse of g0t from the LU decomposition, we can compute
c     g0t^-1 * g0 directly by backsubstituting the columns of G0.
c     See sect. 2.3 in Numerical Recipes or LAPACK Users' Guide
c     sect. 2.3

c     third arg in number of output columns, istate for full
c     matrix, ipart(ik) for just the parts of the matrix needed
c     to contruct fine structure + DOS functions

      do 620 ip=ipi, ipf
        ipart = nsp*(lipotx(ip)+1)**2
        do 590 is1 = 1, istate
        do 590 is2 = 1, ipart
          g0s(is1,is2) = g0(is1, is2 + i0(ip))
  590   continue

        trans = 'NotTransposed'
        call cgetrs(trans, istate, ipart, g0t, istatx,
     $                ipiv, g0s, istatx, info)
        if (info.lt.0) then
            call wlog('    *** Error in cgetrf')
            write(cerr,400) abs(info)
            call wlog('        Argument #'//cerr//
     $              ' had an invalid value.')
        endif

c **** at this point g0s contains the full MS ****

c  pack FMS matrix into an nsp*(lx+1)^2 x nsp*(lx+1)^2 matrix for each
c  ipot

        do 600 is2=1,ipart
        do 600 is1=1,ipart
          gg( is1, is2, ip) = g0s( is1+i0(ip), is2)
 600    continue
 620  continue

      return
      end
      subroutine ggbi( nsp, i0, ipi, ipf, lipotx, g0, tmatrx, g0t, gg,
     1                 toler1, toler2, lcalc, msord)

      implicit real (a-h,o-z)
      implicit integer (i-n)
c  output
c    gg:  (nsp*lx**2, nsp*lx**2, 0:nphasx) submatrix spanning the entire
c          angular momentum basis for each unique potential
c     BiCGStab algorithm: Saad, Iterative methods for ..., p. 220 (1996)

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
c={xparam.h
c nphasx MUST be the same as nphx, the maximum number of unique
c        potentials
c natxx MUST be the same as natx, the maximum number of atoms in the
c       extendeed cluster
c nexx MUST be the same as nex, the maximum number of energy points
      parameter (nphasx=nphx)
      parameter (natxx=natx)
      parameter (nexx=nex)
      parameter (istatx=(lx+1)**2*nclusx*nspx)
      parameter (nkmin=1)
c     parameter (nkmin=-9)
c= xparam.h}
      integer  i0 (0:nphx),  lipotx(0:nphx)

      parameter (one = 1, zero = 0)
      complex coni
      parameter (coni = (0,1))

c**** array of state kets at current energy
      common /stkets/ lrstat(4, istatx), istate

      complex   tmatrx(nspx, istatx)
c     big work matrices
      complex   g0( istatx, istatx), g0t( istatx, istatx)
      logical lcalc
      dimension lcalc(0:lx)

c     Lanczos method variables
      complex xvec( istatx), yvec(istatx), avec(istatx), asve(istatx)
      complex rvec(istatx), pvec(istatx), svec( istatx)
      complex aa, dd, aw, wa, ww
      complex del, delp, omega, chi, psi
c     return matrix containing info about each unique potential
      complex   gg(nspx*(lx+1)**2, nspx*(lx+1)**2, 0:nphasx)

c      notice that in gglu we invert (1-Gt), but here (1-tG).
c     multiply T and G0 matrices together, construct g0t = 1 - T*G0
c     notice that the signs below for g0t ARE correct since 1 is the
c     unit matrix
c     since t is tri-diagonal, this product can be computed in n^2 time
c     also fill up some work matrices for use in eigenvalue and
c     determinant calculations and elsewhere
c     cycle over dimensions of matrix g0t
      do 10 icol = 1,istatx
      do 10 irow = 1,istatx
 10   g0t(irow,icol) = 0

      do 30 icol = 1,istate
        do 20 irow = 1,istate
c         T diagonal contribution T(irow, irow)
          if ( abs( g0(irow, icol)) .gt. toler2 )
     1    g0t(irow,icol)=g0t(irow,icol) - tmatrx(1,irow) * g0(irow,icol) 

c         T off-diagonal contribution T(ist2, irow) in tmatr(2,irow)
c         T off-diagonal contribution T(irow, ist2) in tmatr(2,ist2)
          l1   = lrstat(2,irow)
          m1   = lrstat(3,irow)
          isp1 = lrstat(4,irow)
          m2 = m1+isp1
          if (nsp.eq.2 .and. m2.gt.-l1+1 .and. m2.lt.l1+2) then
c            spin-flip contribution
             ist2 = irow + (-1)**isp1
             if ( abs( g0(ist2, icol)) .gt. toler2)
     1       g0t(irow, icol) = g0t(irow, icol)
     2                   - tmatrx(nsp, ist2) * g0(ist2, icol) 
          endif
 20     continue

        g0t(icol, icol) = g0t(icol, icol) + one
 30   continue

      do 920 ip=ipi, ipf
        ipart = nsp*(lipotx(ip)+1)**2
        do 910 is1 = 1, ipart
          is2 = is1+i0(ip)
          l1   = lrstat(2,is2)
          if (.not.lcalc(l1)) goto 910

c         start first tier with xvec=0
          istart = -1
          msord = 0
          do 40 is = 1, istate
          avec(is) = 0
  40      xvec(is) = 0

c         RESTART here if necessary
  50      continue
          istart = istart+1

          if (istart.gt.0) call matvec( istatx,istate,g0t,xvec,avec,1)
          do 90 is = 1,istate
 90       rvec(is) = - avec(is)
c         rvec = bvec - g0t*xvec , in our case bvec(is) = delta_{is,is2}
          rvec(is2) = rvec(is2) + 1
cc          Check convergence criteria: |r_n+1| < tol
            ipass = 1
            do 92 is = 1, istate
              if ( abs(real(rvec(is))).gt.toler1) goto 93
              if ( abs(aimag(rvec(is))).gt.toler1) goto 93
 92         continue
            ipass = 0
 93         continue
            if (ipass.eq.0) goto 700

          do 95 is = 1, istate
 95       pvec(is) = rvec(is)
          call matvec( istatx,istate,g0t,pvec,avec,1)
          msord = msord + 1

c         choose yvec that del and delp close to one
          call cdot( istatx, istate, avec, avec, aa)
          call cdot( istatx, istate, rvec, avec, wa)
          aw = real(wa) - coni* aimag(wa)
          call cdot( istatx, istate, rvec, rvec, ww)
          dd = aa*ww - aw*wa
          if (abs(dd/aa/ww) .lt.1.e-8) then
            do 96 is = 1,istate
  96        yvec(is) = rvec(is) / ww
          else
            ww = ( ww - aw ) / dd
            aa = ( wa - aa) / dd
            do 97 is = 1,istate
  97        yvec(is) = rvec(is) * aa + avec(is) * ww
          endif
          call cdot( istatx, istate, yvec, rvec, del)

c         it seems ran out of precision for nit>150
          nitx = 30
          do 500 nit = 0, nitx
            call cdot( istatx, istate, yvec, avec, delp)
            omega = del / delp

            do 120 is = 1, istate
 120        svec(is) = rvec(is) - omega * avec(is)
cc          Check convergence criteria: |s_n+1| < tol
            ipass = 1
            do 122 is = 1, istate
              if ( abs(real(svec(is))).gt.toler1) goto 123
              if ( abs(aimag(svec(is))).gt.toler1) goto 123
 122        continue
            ipass = 0
 123        continue
            if (ipass.eq.0)  then
              do 124 is = 1, istate
 124          xvec(is) = xvec(is) + omega*pvec(is)
              goto 700
            endif

            call matvec( istatx,istate,g0t,svec,asve,1)
            msord = msord + 1
            call cdot( istatx, istate, asve, asve, aa)
            call cdot( istatx, istate, asve, svec, wa)
            chi = wa / aa
            do 125 is = 1, istate
 125        xvec(is) = xvec(is) + omega*pvec(is) + chi*svec(is)

            do 130 is = 1, istate
 130        rvec(is) = svec(is) - chi* asve(is)

cc          Check convergence criteria: |r_n+1| < tol
            ipass = 1
            do 370 is = 1, istate
              if ( abs(real(rvec(is))).gt.toler1) goto 380
              if ( abs(aimag(rvec(is))).gt.toler1) goto 380
 370        continue
            ipass = 0
 380        continue
            if (ipass.eq.0) goto 700

c           prepare for next iteration
            call cdot( istatx, istate, yvec, rvec, del)
            psi = del / (delp * chi)

            do 135 is = 1, istate
 135        pvec(is) = rvec(is) + psi * (pvec(is)-chi*avec(is))
            call matvec( istatx,istate,g0t,pvec,avec,1)
            msord = msord + 1

 500      continue
c         restart since ran out of iterations
          goto 50

c         exit if tolerance has been achieved
 700      continue
c         print*, ' BI iterations:', nit + istart*nitx
c         end of BI iterations

c         at this point xvec = (1-tG)**-1 * bvec  with chosen tolerance
c         pack FMS matrix into an nsp*(lx+1)^2 x nsp*(lx+1)^2 matrix 
c         for each ipot
          do 800 is2=1,ipart
            gg( is2, is1, ip) = zero
            do 790 is = 1,istate
              gg( is2, is1, ip) = gg( is2, is1, ip) +
     1        g0( is2+i0(ip), is) * xvec(is)
 790        continue
 800      continue

 910    continue
 920  continue

      return
      end
      subroutine ggrm( nsp, i0, ipi, ipf, lipotx, g0, tmatrx, g0t, gg,
     1                 toler1, toler2, lcalc, msord)

      implicit real (a-h,o-z)
      implicit integer (i-n)
c  output
c    gg:  (nsp*lx**2, nsp*lx**2, 0:nphasx) submatrix spanning the entire
c          angular momentum basis for each unique potential

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
c={xparam.h
c nphasx MUST be the same as nphx, the maximum number of unique
c        potentials
c natxx MUST be the same as natx, the maximum number of atoms in the
c       extendeed cluster
c nexx MUST be the same as nex, the maximum number of energy points
      parameter (nphasx=nphx)
      parameter (natxx=natx)
      parameter (nexx=nex)
      parameter (istatx=(lx+1)**2*nclusx*nspx)
      parameter (nkmin=1)
c     parameter (nkmin=-9)
c= xparam.h}
      integer  i0 (0:nphx),  lipotx(0:nphx)

      parameter (one = 1, zero = 0)
      complex coni
      parameter (coni = (0,1))

c**** array of state kets at current energy
      common /stkets/ lrstat(4, istatx), istate

      complex   tmatrx(nspx, istatx)
c     big work matrices
      complex   g0( istatx, istatx), g0t( istatx, istatx)
      logical lcalc
      dimension lcalc(0:lx)

      complex xvec( istatx), xket( istatx), xbra( istatx)
      complex xketp(istatx), xbrap(istatx)
      complex zvec(istatx), rvec(istatx), svec(istatx)
      complex tket(istatx), tbra(istatx)
      double precision  dum1, dum2
      complex alphac, betac, aa, bb, yy, aac, bbc, gamma
      real alpha, beta
c     return matrix containing info about each unique potential
      complex   gg(nspx*(lx+1)**2, nspx*(lx+1)**2, 0:nphasx)

c      notice that in gglu we invert (1-Gt), but here (1-tG).
c     multiply T and G0 matrices together, construct g0t = 1 - T*G0
c     notice that the signs below for g0t ARE correct since 1 is the
c     unit matrix
c     since t is tri-diagonal, this product can be computed in n^2 time
c     also fill up some work matrices for use in eigenvalue and
c     determinant calculations and elsewhere
c     cycle over dimensions of matrix g0t
      do 10 icol = 1,istatx
      do 10 irow = 1,istatx
 10   g0t(irow,icol) = 0

      do 30 icol = 1,istate
        do 20 irow = 1,istate
c         T diagonal contribution T(irow, irow)
          if ( abs( g0(irow, icol)) .gt. toler2 )
     1    g0t(irow,icol)=g0t(irow,icol) - tmatrx(1,irow) * g0(irow,icol) 

c         T off-diagonal contribution T(ist2, irow) in tmatr(2,irow)
c         T off-diagonal contribution T(irow, ist2) in tmatr(2,ist2)
          l1   = lrstat(2,irow)
          m1   = lrstat(3,irow)
          isp1 = lrstat(4,irow)
          m2 = m1+isp1
          if (nsp.eq.2 .and. m2.gt.-l1+1 .and. m2.lt.l1+2) then
c            spin-flip contribution
             ist2 = irow + (-1)**isp1
             if ( abs( g0(ist2, icol)) .gt. toler2)
     1       g0t(irow, icol) = g0t(irow, icol)
     2                   - tmatrx(nsp, ist2) * g0(ist2, icol) 
          endif
 20     continue

        g0t(icol, icol) = g0t(icol, icol) + one
 30   continue

      do 920 ip=ipi, ipf
        ipart = nsp*(lipotx(ip)+1)**2
        do 910 is1 = 1, ipart
          is2 = is1+i0(ip)
          l1   = lrstat(2,is2)
          if (.not.lcalc(l1)) goto 910

c         start first tier with xvec=0
          istart = -1
          msord=0
          do 40 is = 1, istate
          rvec(is) = 0
  40      xvec(is) = 0

c         RESTART here if necessary
  50      continue
          istart = istart+1

          if (istart.gt.0) call matvec( istatx,istate,g0t,xvec,rvec,1)
c         rvec = g0t*xvec - bvec, in our case bvec(is) = delta_{is,is2}
          rvec(is2) = rvec(is2) - 1
          do 90 is = 1,istate
 90       xket(is) = - rvec(is)
          call cdot( istatx, istate, xket, xket, bb)
          if (abs(bb).eq.0) goto 700

          xfnorm = 1.e0 / real(dble(bb))
          do 91 is = 1, istate
 91       xbra(is) = xket(is) * xfnorm
c         |t> = A |n> ; |n> - xket, |n-1> - xketp
          call matvec ( istatx, istate, g0t, xket, tket, 1)
          msord = msord + 1
          call cdot( istatx, istate, xbra, tket, aa)
          aac = real(aa) - coni*aimag(aa)
          bb = 0
          bbc= 0
          betac = aa
          yy = 1
c         initialize vectors
          do 110 is = 1,istate
            xketp(is) = 0
            xbrap(is) = 0
            zvec(is) = xket(is)
            xvec(is) = xvec(is) + zvec(is)/betac
 110      continue

          do 120 is = 1, istate
 120      svec(is) = tket(is)
          do 130 is = 1, istate
 130      rvec(is) = rvec(is) + svec(is) / betac

c         it seems ran out of precision for nit>150
          nitx = 100
          do 300 nit = 1, nitx
c           use recursion method to calculate a_n+1, b_n, |n+1>, <n+1|
            do 140 is = 1, istate
 140        tket(is) = tket(is) - aa*xket(is) - bb*xketp(is)
            call matvec ( istatx, istate, g0t, xbra, tbra, 2)
            do 150 is = 1, istate
 150        tbra(is) = tbra(is) - aac*xbra(is) - bbc*xbrap(is)
            call cdot( istatx, istate, tbra, tket, bb)
            if (abs(bb).eq.0) goto 700

            bb = sqrt (bb)
            bbc = real(bb) - coni*aimag(bb)
            do 160 is = 1, istate
              xketp(is) = xket(is)
              xbrap(is) = xbra(is)
 160        continue
            do 170 is = 1, istate
              xket(is) = tket(is) / bb
              xbra(is) = tbra(is) / bbc
 170        continue
            call matvec ( istatx, istate, g0t, xket, tket, 1)
            msord = msord + 1
            call cdot( istatx, istate, xbra, tket, aa)
            aac = real(aa) - coni*aimag(aa)
            
c           update iterative solution xvec, 
c           and residual rvec = g0t*xvec - |1>
            alphac = bb / betac
            do 210 is = 1, istate
 210        zvec(is) = xket(is) - alphac * zvec(is)
            do 220 is = 1, istate
 220        svec(is) = tket(is) - alphac * svec(is)

            betac = aa - alphac*bb
            yy = - alphac * yy
            gamma = yy / betac
            do 230 is = 1, istate
 230        xvec(is) = xvec(is) + gamma * zvec(is)
            do 240 is = 1, istate
 240        rvec(is) = rvec(is) + gamma * svec(is)

cc          Check convergence criteria: | rvec | < tol
c           call vecvec( istatx, istate, rvec, rvec, dum2)
c           if (dum2.le.tol) goto 700
cc          Check convergence criteria: | rvec | < tol
            ipass = 1
            do 250 is = 1, istate
              if ( abs(real(rvec(is))).gt.toler1) goto 260
              if ( abs(aimag(rvec(is))).gt.toler1) goto 260
 250        continue
            ipass = 0
 260        continue
            if (ipass.eq.0) goto 700

 300      continue
c         restart since ran out of iterations
          goto 50

c         exit if tolerance has been achieved
 700      continue
c         end of RM iterations

c         at this point xvec = (1-tG)**-1 * bvec  with chosen tolerance
c         pack FMS matrix into an nsp*(lx+1)^2 x nsp*(lx+1)^2 matrix 
c         for each ipot
          do 800 is2=1,ipart
            gg( is2, is1, ip) = zero
            do 790 is = 1,istate
              gg( is2, is1, ip) = gg( is2, is1, ip) +
     1        g0( is2+i0(ip), is) * xvec(is)
 790        continue
 800      continue

 910    continue
 920  continue

      return
      end

      subroutine cdot ( istatx, istate, abra, aket, cc)
c     dot product of two vectors
c     notice that we keep bra  vector as it's complex conjugate
c     thus need to conjugate abra here
      implicit real (a-h,o-z)
      implicit integer (i-n)
      complex coni
      parameter (coni = (0,1))
      complex abra, aket, cc, aa
      dimension abra(istatx), aket(istatx)

      cc = 0
      do 10 is = 1,istate
        aa = real(abra(is)) - coni*aimag(abra(is))
        cc = cc + aa * aket(is)
 10   continue
      return
      end

      subroutine vecvec ( istatx, istate, avec, bvec, cc)
c     dot product of two vectors
      implicit real (a-h,o-z)
      implicit integer (i-n)
      complex avec, bvec
      double precision cc, aa, bb
      dimension avec(istatx), bvec(istatx)

      cc = 0
      do 10 is = 1,istate
        aa = dble(real(avec(is))) * dble(real(bvec(is)))
        bb = dble(aimag(avec(is))) * dble(aimag(bvec(is)))
        cc = cc + aa + bb
 10   continue
      return
      end

      subroutine matvec (istatx, istate, amat, bvec, cvec, itrans)
c     itrans = 1  cvec = amat * bvec
c     itrans = 2  cvec = amat^+ * bvec
c     itrans = 3  cvec = amat^T * bvec
      implicit real (a-h,o-z)
      implicit integer (i-n)
      complex coni, aa
      parameter (coni = (0,1))
      complex amat, bvec, cvec
      dimension amat(istatx, istatx), bvec(istatx), cvec(istatx)

c     initialize cvec
      do 10 is = 1,istatx
 10   cvec(is) = 0

c     cycle over dimensions of amat
      do 20 icol = 1,istate
      do 20 irow = 1,istate
        if (itrans.eq.1) then
          cvec(irow) = cvec(irow) + amat(irow, icol) * bvec(icol)
        elseif(itrans.eq.2) then
          aa = real(amat(irow, icol)) - coni*aimag(amat(irow, icol))
          cvec(icol) = cvec(icol) + aa * bvec(irow)
        else
          cvec(icol) = cvec(icol) + amat(irow, icol) * bvec(irow)
        endif
 20   continue

      return
      end
      subroutine gggm( nsp, i0, ipi, ipf, lipotx, g0, tmatrx, g0t, gg,
     1                 toler1, toler2, lcalc, msord)

      implicit real (a-h,o-z)
      implicit integer (i-n)
c  output
c    gg:  (nsp*lx**2, nsp*lx**2, 0:nphasx) submatrix spanning the entire
c          angular momentum basis for each unique potential
c     Lanczos algorithm: Graves-Morris,Salam, Num.Algor.21,p.213(1999)

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
c={xparam.h
c nphasx MUST be the same as nphx, the maximum number of unique
c        potentials
c natxx MUST be the same as natx, the maximum number of atoms in the
c       extendeed cluster
c nexx MUST be the same as nex, the maximum number of energy points
      parameter (nphasx=nphx)
      parameter (natxx=natx)
      parameter (nexx=nex)
      parameter (istatx=(lx+1)**2*nclusx*nspx)
      parameter (nkmin=1)
c     parameter (nkmin=-9)
c= xparam.h}
      integer  i0 (0:nphx),  lipotx(0:nphx)

      parameter (one = 1, zero = 0)
      complex coni
      parameter (coni = (0,1))

c**** array of state kets at current energy
      common /stkets/ lrstat(4, istatx), istate

      complex   tmatrx(nspx, istatx)
c     big work matrices
      complex   g0( istatx, istatx), g0t( istatx, istatx)
      logical lcalc
      dimension lcalc(0:lx)

c     Lanczos method variables
      complex xvec( istatx), wvec(istatx), x0(istatx), x1(istatx)
      complex avec(istatx), bvec(istatx)
      complex r0(istatx), r1(istatx), t0(istatx), t1(istatx)
      complex aa, dd, aw, wa, ww, e0, e1, alpha, beta, theta, q0, q1
c     return matrix containing info about each unique potential
      complex   gg(nspx*(lx+1)**2, nspx*(lx+1)**2, 0:nphasx)

c      notice that in gglu we invert (1-Gt), but here (1-tG).
c     multiply T and G0 matrices together, construct g0t = 1 - T*G0
c     notice that the signs below for g0t ARE correct since 1 is the
c     unit matrix
c     since t is tri-diagonal, this product can be computed in n^2 time
c     also fill up some work matrices for use in eigenvalue and
c     determinant calculations and elsewhere
c     cycle over dimensions of matrix g0t
      do 10 icol = 1,istatx
      do 10 irow = 1,istatx
 10   g0t(irow,icol) = 0

      do 30 icol = 1,istate
        do 20 irow = 1,istate
c         T diagonal contribution T(irow, irow)
          if ( abs( g0(irow, icol)) .gt. toler2 )
     1    g0t(irow,icol)=g0t(irow,icol) + tmatrx(1,irow) * g0(irow,icol) 

c         T off-diagonal contribution T(ist2, irow) in tmatr(2,irow)
c         T off-diagonal contribution T(irow, ist2) in tmatr(2,ist2)
          l1   = lrstat(2,irow)
          m1   = lrstat(3,irow)
          isp1 = lrstat(4,irow)
          m2 = m1+isp1
          if (nsp.eq.2 .and. m2.gt.-l1+1 .and. m2.lt.l1+2) then
c            spin-flip contribution
             ist2 = irow + (-1)**isp1
             if ( abs( g0(ist2, icol)) .gt. toler2)
     1       g0t(irow, icol) = g0t(irow, icol)
     2                   + tmatrx(nsp, ist2) * g0(ist2, icol) 
          endif
 20     continue

c       g0t(icol, icol) = g0t(icol, icol) + one
 30   continue

      do 920 ip=ipi, ipf
        ipart = nsp*(lipotx(ip)+1)**2
        do 910 is1 = 1, ipart
          is2 = is1+i0(ip)
          l1   = lrstat(2,is2)
          if (.not.lcalc(l1)) goto 910

c         start first tier with xvec=0
          istart = -1
          msord = 0
          do 40 is = 1, istate
          bvec(is) = 0
  40      xvec(is) = 0
c         rvec = bvec - A*xvec , in our case bvec(is) = delta_{is,is2}
          bvec(is2) = 1

c         RESTART here if necessary
  50      continue
          istart = istart+1

          if (istart.gt.0) then
            do 60 is = 1, istate
  60        xvec(is) = xvec(is) + x0(is) / q0
            call matvec( istatx,istate,g0t,xvec,avec,1)
            do 70 is = 1, istate
  70        bvec(is) = avec(is) - xvec(is)
            bvec(is2) = bvec(is2) + 1
          endif
          do 80 is = 1,istate
 80       r0(is) = bvec(is)
          do 90 is = 1,istate
 90       x0(is) = 0
          do 95 is = 1, istate
 95       x1(is) = bvec(is)
          call matvec( istatx,istate,g0t,bvec,r1,1)
          msord = msord + 1

c         choose wvec that del and delp close to one
          call cdot( istatx, istate, r0, r0, ww)
          call cdot( istatx, istate, r1, r1, aa)
          call cdot( istatx, istate, r0, r1, wa)
          aw = real(wa) - coni* aimag(wa)
          dd = aa*ww - aw*wa
          if (abs(dd/aa/ww) .lt.1.e-8) then
            do 96 is = 1,istate
  96        wvec(is) = r0(is) / ww
          else
            ww = ( ww - aw ) / dd
            aa = ( wa - aa) / dd
            do 97 is = 1,istate
  97        wvec(is) = r0(is) * aa + r1(is) * ww
          endif
c         update dot products to avoid round off errors
          call cdot( istatx, istate, wvec, r0, e0)
          call cdot( istatx, istate, wvec, r1, e1)
          q0 = 1
          q1 = 1

c         it seems ran out of precision for nit>150
          nitx = 10
          do 500 nit = 1, nitx
            tol = toler1 * abs(q1) /10
cc          Check convergence criteria: |r1| < tol / 10
cc          so mostly code will not exit here
            ipass = 1
            do 98 is = 1, istate
              if ( abs(real(r1(is))).gt.tol) goto 99
              if ( abs(aimag(r1(is))).gt.tol) goto 99
  98        continue
            ipass = 0
  99        continue
            if (ipass.eq.0) then
              do 100 is = 1, istate
 100          xvec(is) = xvec(is) + x1(is) / q1
              goto 700
            endif

            alpha = e1 / e0
            do 130 is = 1, istate
 130        t0(is) = r1(is) - alpha* r0(is)
            call matvec( istatx,istate,g0t,t0,t1,1)
            msord = msord + 1

            call cdot( istatx, istate, t0, t1, wa)
            call cdot( istatx, istate, t0, t0, ww)
            call cdot( istatx, istate, t1, t1, aa)
            aw = real(wa) - coni* aimag(wa)
            theta = (wa - aa) / (ww - aw)

            do 145 is = 1, istate
 145        r0(is) = t1(is) - theta * t0(is)
            dd = 1- theta
            do 150 is = 1, istate
 150        x0(is) = t0(is) + dd * (x1(is) - alpha*x0(is))
            q0 = dd * (q1 - alpha*q0)
            tol = toler1 * abs(q0)

cc          Check convergence criteria: |r0| < tol
            ipass = 1
            do 370 is = 1, istate
              if ( abs(real(r0(is))).gt.tol) goto 380
              if ( abs(aimag(r0(is))).gt.tol) goto 380
 370        continue
            ipass = 0
 380        continue
            if (ipass.eq.0) then 
              do 390 is = 1, istate
 390          xvec(is) = xvec(is) + x0(is) / q0
              goto 700
            endif

c           prepare for next iteration
            call cdot( istatx, istate, wvec, r0, e0)
            beta = e0 / e1
            do 255 is = 1, istate
 255        t0(is) = r0(is) - beta * r1(is)
            call matvec( istatx,istate,g0t,t0,avec,1)
            msord = msord + 1
            dd = beta * theta
            do 260 is = 1, istate
 260        r1(is) = avec(is) + dd * r1(is)
            call cdot( istatx, istate, wvec, r1, e1)

            dd = beta * (1-theta)
            do 270 is = 1, istate
 270        x1(is) = x0(is) - dd * x1(is) + t0(is)
            q1 = q0 - (1-theta) * beta * q1
 500      continue
c         restart since ran out of iterations
          goto 50

c         exit if tolerance has been achieved
 700      continue
c         end of GM iterations

c         at this point xvec = (1-tG)**-1 * bvec  with chosen tolerance
c         pack FMS matrix into an nsp*(lx+1)^2 x nsp*(lx+1)^2 matrix 
c         for each ipot
          do 800 is2=1,ipart
            gg( is2, is1, ip) = zero
            do 790 is = 1,istate
              gg( is2, is1, ip) = gg( is2, is1, ip) +
     1        g0( is2+i0(ip), is) * xvec(is)
 790        continue
 800      continue

 910    continue
 920  continue

      return
      end
      subroutine ggtf( nsp, i0, ipi, ipf, lipotx, g0, tmatrx, g0t, gg,
     1                 toler1, toler2, lcalc, msord)

      implicit real (a-h,o-z)
      implicit integer (i-n)
c  output
c    gg:  (nsp*lx**2, nsp*lx**2, 0:nphasx) submatrix spanning the entire
c          angular momentum basis for each unique potential
c     TFQMR: Saad, Iterative Methods for Sparse Matrices, p.225 (1996).

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
c={xparam.h
c nphasx MUST be the same as nphx, the maximum number of unique
c        potentials
c natxx MUST be the same as natx, the maximum number of atoms in the
c       extendeed cluster
c nexx MUST be the same as nex, the maximum number of energy points
      parameter (nphasx=nphx)
      parameter (natxx=natx)
      parameter (nexx=nex)
      parameter (istatx=(lx+1)**2*nclusx*nspx)
      parameter (nkmin=1)
c     parameter (nkmin=-9)
c= xparam.h}
      integer  i0 (0:nphx),  lipotx(0:nphx)

      parameter (one = 1, zero = 0)
      complex coni
      parameter (coni = (0,1))

c**** array of state kets at current energy
      common /stkets/ lrstat(4, istatx), istate

      complex   tmatrx(nspx, istatx)
c     big work matrices
      complex   g0( istatx, istatx), g0t( istatx, istatx)
      logical lcalc
      dimension lcalc(0:lx)

      complex xvec(istatx), uvec(istatx), avec(istatx), wvec(istatx)
      complex dvec(istatx), rvec(istatx), vvec(istatx)
      complex alpha, beta, aa, rho, eta
      real tau, nu, cm, err
c     return matrix containing info about each unique potential
      complex   gg(nspx*(lx+1)**2, nspx*(lx+1)**2, 0:nphasx)

c      notice that in gglu we invert (1-Gt), but here (1-tG).
c     multiply T and G0 matrices together, construct g0t = 1 - T*G0
c     notice that the signs below for g0t ARE correct since 1 is the
c     unit matrix
c     since t is tri-diagonal, this product can be computed in n^2 time
c     also fill up some work matrices for use in eigenvalue and
c     determinant calculations and elsewhere
c     cycle over dimensions of matrix g0t
      do 10 icol = 1,istatx
      do 10 irow = 1,istatx
 10   g0t(irow,icol) = 0

      do 30 icol = 1,istate
        do 20 irow = 1,istate
c         T diagonal contribution T(irow, irow)
          if ( abs( g0(irow, icol)) .gt. toler2 )
     1    g0t(irow,icol)=g0t(irow,icol) - tmatrx(1,irow) * g0(irow,icol) 

c         T off-diagonal contribution T(ist2, irow) in tmatr(2,irow)
c         T off-diagonal contribution T(irow, ist2) in tmatr(2,ist2)
          l1   = lrstat(2,irow)
          m1   = lrstat(3,irow)
          isp1 = lrstat(4,irow)
          m2 = m1+isp1
          if (nsp.eq.2 .and. m2.gt.-l1+1 .and. m2.lt.l1+2) then
c            spin-flip contribution
             ist2 = irow + (-1)**isp1
             if ( abs( g0(ist2, icol)) .gt. toler2)
     1       g0t(irow, icol) = g0t(irow, icol)
     2                   - tmatrx(nsp, ist2) * g0(ist2, icol) 
          endif
 20     continue

        g0t(icol, icol) = g0t(icol, icol) + one
 30   continue

      do 920 ip=ipi, ipf
        ipart = nsp*(lipotx(ip)+1)**2
        do 910 is1 = 1, ipart
          is2 = is1+i0(ip)
          l1   = lrstat(2,is2)
          if (.not.lcalc(l1)) goto 910

c         start first tier with xvec=0
          istart = -1
          msord = 0
          do 40 is = 1, istate
          rvec(is) = 0
          avec(is) = 0
  40      xvec(is) = 0

c         RESTART here if necessary
  50      continue
          istart = istart+1

          if (istart.gt.0) call matvec( istatx,istate,g0t,xvec,avec,1)
          do 90 is = 1,istate
 90       uvec(is) = - avec(is)
c         uvec = bvec - g0t*xvec , in our case bvec(is) = delta_{is,is2}
          uvec(is2) = uvec(is2) + 1
          call matvec( istatx,istate,g0t,uvec,avec,1)
          msord = msord + 1
          do 95 is = 1, istate
 95       wvec(is) = uvec(is)
          do 96 is = 1, istate
 96       vvec(is) = avec(is)
          do 97 is = 1, istate
 97       dvec(is) = 0
          call cdot( istatx, istate, uvec, uvec, aa)
          tau = sqrt(real(aa))
          nu = 0
          eta = 0
c         choose rvec = uvec /aa so that dot products are about 1
          do 98 is = 1, istate
 98       rvec(is) = uvec(is) / aa
          rho = 1

c         it seems ran out of precision for nit>150
          nitx = 20
          do 300 nit = 0, nitx
            if (mod(nit,2).eq.0) then
              call cdot( istatx, istate, rvec, vvec, aa)
              alpha = rho / aa
            else
              call matvec( istatx,istate,g0t,uvec,avec,1)
              msord = msord + 1
            endif

            do 115 is = 1, istate
 115        wvec(is) = wvec(is) - alpha * avec(is)

            aa = nu**2 * eta / alpha
            do 120 is = 1, istate
 120        dvec(is) = uvec(is) + aa * dvec(is)

            call cdot( istatx, istate, wvec, wvec, aa)
            nu = sqrt(real(aa)) / tau
            cm = 1 / sqrt(1+nu**2)
            tau = tau * nu * cm
            eta = cm**2 * alpha
            do 140 is = 1, istate
 140        xvec(is) = xvec(is) + eta * dvec(is)

cc          Check convergence criteria: | rvec | < tol
            err = (1.e0 + nit) / istate
            err = tau * sqrt(err) * 10
            if ( abs(err).lt.toler1) goto 700

            if (mod(nit,2) .ne.0) then
              aa = rho
              call cdot( istatx, istate, rvec, wvec, rho)
              beta = rho / aa
              do 210 is = 1, istate
 210          uvec(is) = wvec(is) + beta * uvec(is)

              do 215 is = 1, istate
 215          vvec(is) = beta * ( avec(is) + beta * vvec(is))
              call matvec( istatx,istate,g0t,uvec,avec,1)
              msord = msord + 1
              do 220 is = 1, istate
 220          vvec(is) = avec(is) + vvec(is)
            else
              do 230 is = 1, istate
 230          uvec(is) = uvec(is) - alpha * vvec(is)
            endif
 300      continue
c         restart since ran out of iterations
          goto 50

c         exit if tolerance has been achieved
 700      continue
c         end of TFQMR iterations

c         at this point xvec = (1-tG)**-1 * bvec  with chosen tolerance
c         pack FMS matrix into an nsp*(lx+1)^2 x nsp*(lx+1)^2 matrix 
c         for each ipot
          do 800 is2=1,ipart
            gg( is2, is1, ip) = zero
            do 790 is = 1,istate
              gg( is2, is1, ip) = gg( is2, is1, ip) +
     1        g0( is2+i0(ip), is) * xvec(is)
 790        continue
 800      continue

 910    continue
 920  continue

      return
      end
      subroutine yprep(iph0, nat, inclus, npot, iphat, rmax, rat,
     $            izx, rdirec)
c    yprep is the same as xprep for negative idwopt
c    simlifies calls in SCF and LDOS where DW factors should not enter

      implicit real (a-h,o-z)
      implicit integer (i-n)

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
c={xparam.h
c nphasx MUST be the same as nphx, the maximum number of unique
c        potentials
c natxx MUST be the same as natx, the maximum number of atoms in the
c       extendeed cluster
c nexx MUST be the same as nex, the maximum number of energy points
      parameter (nphasx=nphx)
      parameter (natxx=natx)
      parameter (nexx=nex)
      parameter (istatx=(lx+1)**2*nclusx*nspx)
      parameter (nkmin=1)
c     parameter (nkmin=-9)
c= xparam.h}
c====================================================================
c  This header file contains the structural information about the
c  cluster to be used for the full multiple scattering calculation

      common /xstruc/ xphi(nclusx,nclusx), xrat(3,nclusx),
     $            iphx(nclusx)
      save /xstruc/

c  xphi:  matrix of angles between the z axis and pairs of atoms
c  xrat:  xyz coordinates of the atoms in the cluster, the first
c         npot+1 entries are examples of each unique potential
c  iphx:  potential indeces of each atom in the cluster, ordered like
c         xrat
c********************************************************************
c**** rotation matrices for entire cluster
c
      complex drix
      common /rotx/ drix(-lx:lx, -lx:lx, 0:lx, 0:1, nclusx, nclusx)
      save /rotx/
c********************************************************************
c common blocks for saving rotation matrices between xanes and rotxan
       integer    jsavx, jsav, jbmagk
       parameter (jsavx = 150, roteps = 1.e-12,jbmagk=-9999)
       dimension drisav(-lx:lx,-lx:lx,0:lx,jsavx), betsav(jsavx)
       integer   ldsav(jsavx), mdsav(jsavx)
       common /rotsav/  drisav, betsav, ldsav, mdsav, jsav
       save  /rotsav/
c********************************************************************
c**** legendre polynomial normalization constants
c
      common /lnlm/ xnlm(0:lx,0:lx)
      save   /lnlm/

c********************************************************************
c**** legendre polynomial normalization constants
c
      common /xdwf/ sigsqr(nclusx,nclusx)
      save   /xdwf/

c  end of xstruc.h
c********************************************************************
      parameter(zero=0.e0)
      parameter (bohr = 0.529 177 249e0)
      integer   iphat(natxx), iphat2(natxx), izx(0:nphasx), izpair(0:2)
      dimension rat(3,natxx), rat2(3,natxx)
      double precision ra(natxx)
      character*78 line
c     sigms is written in double precision.  these are the variables
c     that it uses
      double precision dtemp, dthet, drs, dsigsq, pair(3,0:2)
      double precision sig2mx, sig2x(0:nphx,0:nphx)
c     iwarn - needed to wrtite waqrning just one time
      integer iwarn
      save iwarn
      data  iwarn /0/

c  initialize geometrical arrays
      do 30 i=1,nclusx
        do 10 j=1,nclusx
          xphi(j,i) = zero
 10     continue
        do 20 j=1,3
          xrat(j,i) = zero
 20     continue
        iphx(i) = 0
 30   continue
      inclus = 0

c --- find the central atom, ipot=iph0 (iph0=0 for the absorbing atom)
      icen = 0
      do 40 i=1,nat
        iphat2(i) = iphat(i)
        if (iphat(i).eq.iph0) then
            if (icen.eq.0) then
                icen = i
            elseif (iph0.eq.0) then
                call wlog('* * * ERROR!  More than one atom '//
     $                      'in the extended cluster have ipot=0')
                call wlog('      You may only have one central atom.')
                call wlog('      Stopping in xprep.')
                call par_stop('YPREP-1')
            endif
        endif
 40   continue
c --- make sure central atom is at (0,0,0)
      do 45 i=1,nat
        rat2(1,i) = rat(1,i)-rat(1,icen)
        rat2(2,i) = rat(2,i)-rat(2,icen)
        rat2(3,i) = rat(3,i)-rat(3,icen)
 45   continue

c --- sort the atoms from extended cluster by distance from central
c     atom.
      call atheap(nat, rat2, iphat2, ra)

c --- define cluster from extended cluster by as those closer than
c     rmax to central atom
      inclus=0
      rmax2 = rmax**2
      do 50 i=1,nat
        rr = (rat2(1,i)**2 + rat2(2,i)**2 + rat2(3,i)**2)
        if (rr.gt.rmax2) then
            inclus = i-1
            goto 60
        endif
 50   continue
 60   continue
      if (inclus.eq.0) inclus=nat

c --- sanity check size of cluster
      if (inclus.gt.nclusx) then
        if (iwarn.eq.0) then
          call wlog('* * * WARNING preparing cluster for '//
     $                'FMS calculation.')
          write(line,400) inclus
 400      format('      You specified a cluster of ', i3,
     $                ' atoms for the FMS calculation.')
          call wlog(line)
          write(line,410)nclusx
          call wlog(line)
 410      format('      This exceeds the hard wired limit of ', i3,
     $                ' atoms.')
          write(line,420)nclusx
          call wlog(line)
 420      format('      The cluster size was reset to ', i3,
     $                ' and the calculation will continue.')
          iwarn = 1
        endif
        inclus = nclusx
      endif

c --- make the first few entries in xrat represent each of the
c     unique potentials, sorting around the chosen center
c     (iph0=0 for the absorbing atom)
c     call sortat(iph0, inclus, npot, iphat2, iphx, rat2, xrat)
      do 430 iat = 1, inclus
          iphx(iat) = iphat2(iat)
          xrat(1,iat) = real (rat2(1,iat))
          xrat(2,iat) = real (rat2(2,iat))
          xrat(3,iat) = real (rat2(3,iat))
 430  continue


c --- Calculate and store rotation matrix elements and phi angles
c     the k loop calculates the forward then the backward rotation
c     for an atom pair (ij). k = 0-->forward, 1-->backward
      call rotint
      lplus1 = lx+1
      mplus1 = lx+1
      do 150  i=1,inclus
        do 140 j=1,inclus
          rr = (xrat(1,i)-xrat(1,j))**2 + (xrat(2,i)-xrat(2,j))**2
     1       + (xrat(3,i)-xrat(3,j))**2
c         if (rr.gt.rdirec**2) goto 140

          call getang(nclusx, xrat, i, j, xbeta, xphi(i,j))
          if (i.eq.j) goto 140
          do 130 k=0,1
            if (k.eq.1) xbeta = (-1) * xbeta
            call rotxan(lplus1, mplus1, xbeta, i, j, k)
 130      continue
 140    continue
 150  continue

c --- calculate spherical harmonic normalization factors
      call xanlm(lplus1,mplus1)

      do 200 iat2=1,nclusx
      do 200 iat1=1,nclusx
        sigsqr(iat1,iat2) = zero
 200  continue

      return
      end
      subroutine atheap(nat, rat, iphat, ra)

c--------------------------------------------------------------
c  copyright 1993 university of washington         bruce ravel
c  modified by alexei ankudinov in march 1999
c--------------------------------------------------------------
      implicit real (a-h,o-z)
      implicit integer (i-n)
c      implicit double precision (a-h,o-z)
c-------------------------------------------------------------------
c  heapsort adapted from numerical recipes.  sort atoms by distance.
c  all the pesky little do loops are for transferring rows
c  of temp into toss.
c-------------------------------------------------------------------
c  alexei ankudinov: needed to avoid unnecessary permutations when atoms
c  are at the same distance from the central atom, in order to comply 
c  feff document: the sample atom should be the nearest to absorber or
c  first in the list among equidistant
c  Add small contribution 10**-8 * number to the sorting variable ra
c  in order to achieve this.
c-------------------------------------------------------------------
c  natx:   dimension parameter from calling program
c-------------------------------------------------------------------
      dimension rat(3, nat), toss(3), iphat(nat)
      double precision ra(nat), dum 

      if (nat.lt.2) return

      l=0
      do 10 i=1,nat
         ra(i) = dble( rat(1,i)**2 + rat(2,i)**2 + rat(3,i)**2 ) +
     1           i*1.d-8
c        small addition at to prefer the old ordering
         if (l.eq.0 .and.i.gt.1) then
             if (ra(i).lt.ra(i-1)) l=1
         endif
  10  continue
c     check if array is already in order
      if (l.eq.0) return

      l  = nat/2+1
      ir = nat
 110  continue
         if (l.gt.1) then
            l = l-1
            do 120 index=1,3
               toss(index)=rat(index,l)
 120        continue
            itoss = iphat(l)
            dum = ra(l)
         else
            do 130 index=1,3
               toss(index) = rat(index,ir)
 130        continue
            itoss = iphat(ir)
            dum = ra(ir)
            do 140 index=1,3
               rat(index,ir) = rat(index,1)
 140        continue
            iphat(ir) = iphat(1)
            ra(ir) = ra(1)
            ir=ir-1
            if (ir.eq.1) then
               do 150 index=1,3
                  rat(index,1)=toss(index)
 150           continue
               iphat(1) = itoss
               ra(1) = dum
c              sort is finished
               goto 300
            endif
         endif
         i=l
         j=l+l

 160     if (j.le.ir) then
            if (j.lt.ir) then
               if ( ra(j) .lt. ra(j+1) ) then
                  j  = j + 1
               endif
            endif

            if ( dum .lt. ra(j) ) then
               do 170 index=1,3
                  rat(index,i) = rat(index,j)
 170           continue
               iphat(i) = iphat(j)
               ra(i) = ra(j)
               i=j
               j=j+j
            else
               j=ir+1
            endif
            goto 160
         endif

         do 180 index=1,3
            rat(index,i) = toss(index)
 180     continue
         iphat(i) = itoss
         ra(i) = dum

      goto 110
 300  continue

      return
c end subroutine atheap
      end
      subroutine getang(nclusx, rat, i, j, theta, phi)

c------------------------------------------------------------------
c  determine theta and phi polar angles of the vector between two
c  atom positions
c
c  inputs
c    rat:   (3,nclusx) x,y,z of all atoms in cluster
c    i, j:  indices of atoms at ends of vector Ri-Rj
c
c  outputs
c    theta: polar angle theta of vector Ri-Rj
c    phi:   polar angle phi of vector Ri-Rj
c------------------------------------------------------------------

      implicit real (a-h,o-z)
      implicit integer (i-n)

c       include 'dim.h'
c       include 'xparam.h'
      dimension rat(3,nclusx)
      parameter(tiny=1.e-7, zero=0.e0, pi=3.141592654)

      x = rat(1,i) - rat(1,j)
      y = rat(2,i) - rat(2,j)
      z = rat(3,i) - rat(3,j)
      r = sqrt(x**2 + y**2 + z**2)

c  this fails to calculate phi correctly for, as an example,
c  x=0.5e-7 and y=2e-7.  However, those numbers are below the
c  precision of the numbers stored in potph.bin.

      phi = zero
      theta  = zero
      if (i.ne.j) then
c           phi = atan2(y,x)
c        all of these conditionals will do the work for a machine that
c        cannot correctly handle a zero value for the second argument
c        of atan2
          if (abs(x).lt.tiny) then
              if (abs(y).lt.tiny) then
                  phi = zero
              elseif (y.gt.tiny) then
                  phi = pi/2
              else
                  phi = -pi/2
              endif
          else
              phi = atan2(y,x)
          endif
          if (r.gt.tiny) then
            if (z.le.-r) then
             theta = pi
            elseif ( z.lt.r) then
             theta = acos(z/r)
            endif
          endif
      endif

      return
c  end subroutine getang
      end


c====================================================================
      subroutine rotxan (lxp1, mxp1, betax, i, j, k)
      implicit real (a-h,o-z)

c     input:  lxp1, mxp1: lmax+1 & mmax+1, largest L states in matrix
c             betax is the rotation angle
c             i and j are the indeces of the atoms, thus denote
c                 which pair of atoms this is the rotation matrix for
c             k=0 for forward rotation, k=1 for backward rotation
c     output: drix(L,k,j,i) in common /rotx/
c+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+
c     adapted by BR from rot3i, version for genfmt by SIZ
c        new data structure for rotation matrices to accomodate
c        xanes calculation
c+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+
c     subroutine rot3 calculates rotation matrices for l = 0,lxp1-1

c     subroutine rot3 calculates the beta dependence of rotation
c     matrix elements using recursion of an iterated version of
c     formula (4.4.1) in edmonds.
c
c     first written:(september 17,1986) by j. mustre
c     version 2  (17 sep 86)
c     version 3  (22 feb 87) modified by j. rehr
c     version for genfmt, modified by s. zabinsky, Sept 1991
c     Initialized dri0.  Some elements may be used before being
c        initialized elsewhere -- rot3i needs to be carefully
c        checked.  S. Zabinsky, April 1993
c
c******************** warning****************************************
c     lxx must be at least lxp1 or overwriting will occur
c     nmax must be at least nm or overwriting will occur
c--------------------------------------------------------------------
c     notation dri0(l,m,n) =  drot_i(l'm'n')
c     l = l'+1, n' = n-l, m' = m-l, primes denoting subscripts
c     thus dri0(1,1,1) corresponds to the rotation matrix with
c     l' = 0, and n' and m' = 0; dri0(3,5,5) : l' = 2,n' = 2,m' = 2.
c--------------------------------------------------------------------

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
c={xparam.h
c nphasx MUST be the same as nphx, the maximum number of unique
c        potentials
c natxx MUST be the same as natx, the maximum number of atoms in the
c       extendeed cluster
c nexx MUST be the same as nex, the maximum number of energy points
      parameter (nphasx=nphx)
      parameter (natxx=natx)
      parameter (nexx=nex)
      parameter (istatx=(lx+1)**2*nclusx*nspx)
      parameter (nkmin=1)
c     parameter (nkmin=-9)
c= xparam.h}
c       include 'xstruc.h'
c====================================================================
c  This header file contains the structural information about the
c  cluster to be used for the full multiple scattering calculation

      common /xstruc/ xphi(nclusx,nclusx), xrat(3,nclusx),
     $            iphx(nclusx)
      save /xstruc/

c  xphi:  matrix of angles between the z axis and pairs of atoms
c  xrat:  xyz coordinates of the atoms in the cluster, the first
c         npot+1 entries are examples of each unique potential
c  iphx:  potential indeces of each atom in the cluster, ordered like
c         xrat
c********************************************************************
c**** rotation matrices for entire cluster
c
      complex drix
      common /rotx/ drix(-lx:lx, -lx:lx, 0:lx, 0:1, nclusx, nclusx)
      save /rotx/
c********************************************************************
c#mn{
c common blocks for saving rotation matrices between xanes and rotxan
       integer    jsavx, jsav, jbmagk
       parameter (jsavx = 150, roteps = 1.e-12,jbmagk=-9999)
       dimension drisav(-lx:lx,-lx:lx,0:lx,jsavx), betsav(jsavx)
       integer   ldsav(jsavx), mdsav(jsavx)
       common /rotsav/  drisav, betsav, ldsav, mdsav, jsav
       save  /rotsav/
c#mn}

c********************************************************************
c**** legendre polynomial normalization constants
c
      common /lnlm/ xnlm(0:lx,0:lx)
      save   /lnlm/

c********************************************************************
c**** legendre polynomial normalization constants
c
      common /xdwf/ sigsqr(nclusx,nclusx)
      save   /xdwf/

c  end of xstruc.h
c********************************************************************
      parameter (one = 1, zero = 0)
c      needed for commented out diagnostic file
c      logical open
      parameter(lxx=24)
      parameter (pi = 3.14159 26535 89793 23846 26433e0)
      complex coni, dum 
      parameter (coni = (0,1))

c     dri0 is larger than needed for genfmt, but necessary for
c     this calculation algorithm.  Copy result into smaller
c     dri arrays (in common) at end of this routine.
      dimension  dri0 (lxx+1, 2*lxx+1, 2*lxx+1)

c#mn{
c  check whether a rotation matrix for this {beta(ileg),lxp1,mxp1} has
c  been calculated and saved.  If so, just use the saved value
       do 90 isav = 1, jsav
          if (betsav(isav).eq.jbmagk) go to 95
          if ((lxp1.eq.ldsav(isav)).and.(mxp1.eq.mdsav(isav)).and.
     $         (abs(betax-betsav(isav)).le.roteps) ) then
cc             print*, 'using drisav for ', isav, betax, lxp1, mxp1
             do 85 il = 0, lx
             do 85 m1 = -il, il
             do 85 m2 = -il, il
               drix(m2,m1,il,k,j,i)=cmplx(drisav(m2,m1,il,isav),zero)
 85          continue
             go to 770
          end if
 90    continue
 95    continue
c#mn}


c     initialize dri0
      do 150 in = 1, 2*lxx+1
        do 150 im = 1, 2*lxx+1
          do 150 il = 1, lxx+1
            dri0(il,im,in) = zero
 150  continue

      nm  = mxp1
      ndm = lxp1+nm-1
      xc  = cos(betax/2)
      xs  = sin(betax/2)
      s   = sin(betax)
      dri0(1,1,1) =  1
      dri0(2,1,1) =  xc**2
      dri0(2,1,2) =  s/sqrt(2*one)
      dri0(2,1,3) =  xs**2
      dri0(2,2,1) = -dri0(2,1,2)
      dri0(2,2,2) =  cos(betax)
      dri0(2,2,3) =  dri0(2,1,2)
      dri0(2,3,1) =  dri0(2,1,3)
      dri0(2,3,2) = -dri0(2,2,3)
      dri0(2,3,3) =  dri0(2,1,1)
      do 230  l = 3, lxp1
        ln = 2*l - 1
        lm = 2*l - 3
        if (ln .gt. ndm)  ln = ndm
        if (lm .gt. ndm)  lm = ndm
        do 220  n = 1, ln
          do 210  m = 1, lm
            t1   = (2*l-1-n) * (2*l-2-n)
            t    = (2*l-1-m) * (2*l-2-m)
            f1   = sqrt(t1/t)
            f2   = sqrt( (2*l-1-n) * (n-1) / t )
            t3   = (n-2) * (n-1)
            f3   = sqrt(t3/t)
            dlnm = f1 * xc**2 * dri0(l-1,n,m)
            if (n-1 .gt. 0) dlnm = dlnm - f2*s*dri0(l-1,n-1,m)
            if (n-2 .gt. 0) dlnm = dlnm + f3*xs**2*dri0(l-1,n-2,m)
            dri0(l,n,m) = dlnm
            if (n .gt. (2*l-3))
     $                  dri0(l,m,n) = (-1)**(n-m) * dri0(l,n,m)
 210      continue
          if (n .gt. (2*l-3)) then
              dri0(l,2*l-2,2*l-2) =  dri0(l,2,2)
              dri0(l,2*l-1,2*l-2) = -dri0(l,1,2)
              dri0(l,2*l-2,2*l-1) = -dri0(l,2,1)
              dri0(l,2*l-1,2*l-1) =  dri0(l,1,1)
          endif
 220    continue
 230  continue


c     initialize drix
      do 310 il = 0, lx
      do 310 m1 = -lx, lx
      do 310 m2 = -lx, lx
        drix(m2,m1,il,k,j,i) = cmplx(zero,zero)
        drix(m2,m1,il,k,i,i) = cmplx(zero,zero)
 310  continue

c     Copy result into drix(...,k,j,i) in /rotx/
      do 390  il = 1, lxp1
        mmx = min (il-1, mxp1-1)
        do 380  m1 = -mmx, mmx
        do 380  m2 = -mmx, mmx
          drix(m2, m1, il-1, k, j, i)=cmplx(dri0(il,m1+il,m2+il),zero)
 380    continue
 390  continue
c#mn{
c      save dri if there's room
       if (jsav.lt.jsavx) then
          jsav = jsav + 1
cc          print*, 'saving dri to ',  jsav, betax, lxp1, mxp1
          betsav(jsav) = betax
          ldsav(jsav)  = lxp1
          mdsav(jsav)  = mxp1
          do 720 il = 0, lx
          do 720 m1 = -il, il
          do 720 m2 = -il, il
            drisav(m2,m1,il,jsav) = real(drix(m2,m1,il,k,j,i))
 720      continue
       else
cc          print*, 'not saving dri to ',  betax, lxp1, mxp1
       end if
 770   continue
c#mn}

c-----test sum rule on d
c       if (idbg(1).eq.1) then
c           inquire(file='rotmat.dat', opened=open)
c           if (.not.open) then
c               iun = nxtunt(25)
c               open (iun,file='rotmat.dat',status='unknown')
c           endif
c           write(iun,*)'  '
c           write(iun,*)'atom #s : ',i,j
c           write(iun,*)  ' il, im, sum, beta'
c           write(iun,*) ' (drix(il,im,in,k,j,i),in = -il,il)'
c           do 880 il = 0,lxp1-1
c             do 870 im = -il,il
c               sum = 0
c               do 850 in = -il,il
c                 term = drix(in,im,il,k,j,i)
c                 sum = sum+term**2
c  850           continue
c               write(iun,860) il,im,sum,betax
c               write(iun,862) (drix(in,im,il,k,j,i),in = -il,il)
c  860          format(2i3,1x,f16.12,1x,f8.4)
c  862          format(5f14.6)
c  870         continue
c  880       continue
c c          close(iun)
c       endif
c-----end test------------------------

        do 920 il = 0, lx
        do 920 m1 = -il, il
          dum = coni * m1 * (xphi(i,j)-pi)
          if (k.eq.1) dum = -dum
          dum = exp( dum )
          do 910 m2 = -il, il
            if (k.eq.1) then
              drix(m2,m1,il,k,j,i) = drix(m2,m1,il,k,j,i) * dum
            else
              drix(m1,m2,il,k,j,i) = drix(m1,m2,il,k,j,i) * dum
            endif
 910       continue
 920     continue

      return
c  end subroutine rotxan
      end
c====================================================================
c#mn{
       subroutine rotint
       implicit real (a-h,o-z)
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
c={xparam.h
c nphasx MUST be the same as nphx, the maximum number of unique
c        potentials
c natxx MUST be the same as natx, the maximum number of atoms in the
c       extendeed cluster
c nexx MUST be the same as nex, the maximum number of energy points
      parameter (nphasx=nphx)
      parameter (natxx=natx)
      parameter (nexx=nex)
      parameter (istatx=(lx+1)**2*nclusx*nspx)
      parameter (nkmin=1)
c     parameter (nkmin=-9)
c= xparam.h}
c        include 'xstruc.h'
c====================================================================
c  This header file contains the structural information about the
c  cluster to be used for the full multiple scattering calculation

      common /xstruc/ xphi(nclusx,nclusx), xrat(3,nclusx),
     $            iphx(nclusx)
      save /xstruc/

c  xphi:  matrix of angles between the z axis and pairs of atoms
c  xrat:  xyz coordinates of the atoms in the cluster, the first
c         npot+1 entries are examples of each unique potential
c  iphx:  potential indeces of each atom in the cluster, ordered like
c         xrat
c********************************************************************
c**** rotation matrices for entire cluster
c
      complex drix
      common /rotx/ drix(-lx:lx, -lx:lx, 0:lx, 0:1, nclusx, nclusx)
      save /rotx/
c********************************************************************
c#mn{
c common blocks for saving rotation matrices between xanes and rotxan
       integer    jsavx, jsav, jbmagk
       parameter (jsavx = 150, roteps = 1.e-12,jbmagk=-9999)
       dimension drisav(-lx:lx,-lx:lx,0:lx,jsavx), betsav(jsavx)
       integer   ldsav(jsavx), mdsav(jsavx)
       common /rotsav/  drisav, betsav, ldsav, mdsav, jsav
       save  /rotsav/
c#mn}

c********************************************************************
c**** legendre polynomial normalization constants
c
      common /lnlm/ xnlm(0:lx,0:lx)
      save   /lnlm/

c********************************************************************
c**** legendre polynomial normalization constants
c
      common /xdwf/ sigsqr(nclusx,nclusx)
      save   /xdwf/

c  end of xstruc.h
c********************************************************************
c initialize /rotsav/
       jsav = 0
       do 100 js = 1, jsavx
          betsav(js) = jbmagk
          ldsav(js)  = 0
          mdsav(js)  = 0
          do 90 il  = 0, lx
             do 80 m1 = -lx, lx
                do 70 m2 = -lx, lx
                   drisav(m2,m1,il,js) = 0
 70             continue
 80          continue
 90       continue
 100   continue
       return
c#mn}
       end
      subroutine sortat(iph0, nat, npot, iphat, iphx, rat, xrat)

      implicit real (a-h,o-z)
      implicit integer (i-n)
c--------------------------------------------------------------------
c  this subroutine sorts the atoms in xrat such that the first npot
c  entries are each a representative atom of a unique potential.  This
c  will mean that the upper left corner of the full MS matrix will
c  contain all of the information needed to compute the fine structure
c  and all of the electron densities.
c  NOTA BENE:  the atoms *must* have already been sorted by radial
c    distance!
c--------------------------------------------------------------------
c  input:
c    iph0:    potential index for central atom in LDOS (added by ala)
c                (iph0=0 for absorbing atom as the central atom)
c     nat:    number of atoms in cluster
c    npot:    number of unique potentials in cluster
c    iphat:   (nclusx) potential index of each atom in cluster as read
c             from geometry file
c    rat:     (3, nclusx) coordinates of each atom in cluster as read
c             from geometry file
c  output:
c    iphx:    (nclusx) potential index of each atom in cluster sorted
c             so that the first npot+1 entries are examples of each
c             ipot
c    xrat:    (3, nclusx) coordinates of each atom in cluster sorted
c             so that the first npot+1 entries are examples of each
c             ipot
c--------------------------------------------------------------------

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
c={xparam.h
c nphasx MUST be the same as nphx, the maximum number of unique
c        potentials
c natxx MUST be the same as natx, the maximum number of atoms in the
c       extendeed cluster
c nexx MUST be the same as nex, the maximum number of energy points
      parameter (nphasx=nphx)
      parameter (natxx=natx)
      parameter (nexx=nex)
      parameter (istatx=(lx+1)**2*nclusx*nspx)
      parameter (nkmin=1)
c     parameter (nkmin=-9)
c= xparam.h}
      dimension rat(3,natxx), xrat(3,nclusx)
      integer   iphat(natxx), iphx(nclusx), ipoint
      dimension ipoint(0:nphasx)
      integer iph0, ip, ilast

      do 10 i=0,nphasx
        ipoint(i) = 0
 10   continue
      do 30 ic=1,nat
        iphx(ic) = iphat(ic)
        do 20 ix=1,3
          xrat(ix,ic) = rat(ix,ic)
 20     continue
 30   continue

c     (iph0=0 for absorbing atom as the central atom)
      if (iphx(1).ne.iph0) then
          call wlog('* * * ERROR in sortat * * *')
          call wlog('            The first atom in xrat is not '//
     $                'the central atom.')
          call wlog('            Complain to Bruce immediately!')
          call par_stop('SORTAT-1')
      endif

c       if (idbg(4).eq.1) print*,'SORTAT: nat,npot: ',nat,npot
c       if (idbg(4).eq.1) print*,'SORTAT: xcen,ycen,zcen: ',
c      $            xcen,ycen,zcen

c --- find the example of each unique potential that is closest to the
c     central atom.  This will presumably be well within the cluster
c     that was used to compute the overlapped potentials
      ipoint(iph0) = 1
      do 150 ip=0,npot
        if (ip .ne. iph0) then
          do 130 iat=2,nat
            if (iphx(iat).eq.ip .and. ipoint(ip).eq.0) then
                ipoint(ip) = iat
c                print*,'>>>>> ip, ipoint(ip)', ip, ipoint(ip)
            endif
 130      continue
        endif
 150  continue

c --- now swap the first few atoms with the atoms found above
      do 200 ip=0,npot

c ----- some potentials might not be in the xanes cluster
        if (ipoint(ip).eq.0) goto 200
c ----- don't swap two potentials if examples live in the first npot
c       entries
        if (ipoint(ip).le.ip+1) goto 200

        xx  = xrat(1,1+ip)
        yy  = xrat(2,1+ip)
        zz  = xrat(3,1+ip)
        iph = iphx(1+ip)

        xrat(1,1+ip) = xrat(1,ipoint(ip))
        xrat(2,1+ip) = xrat(2,ipoint(ip))
        xrat(3,1+ip) = xrat(3,ipoint(ip))
        iphx(1+ip)  = iphx(ipoint(ip))

        xrat(1,ipoint(ip)) = xx
        xrat(2,ipoint(ip)) = yy
        xrat(3,ipoint(ip)) = zz
        iphx(ipoint(ip))  = iph

c       added by ala
c       check that substituted atom was not some ip example
c          ???BR Jan 16 1998???
        do 190 ipp = ip+1, npot
          if (ipoint(ipp).eq.ip+1) ipoint(ipp) = ipoint(ip)
  190   continue
c       set the correct pointer to ip example
        ipoint(ip) = ip+1

 200  continue

c     added by ala
c     Notice that fms will take the last atom of given type ip
c     from first npot atoms in the list as an example for ip.
c     Make more permutaions if necesary.
      ilast = -1
      nmin = min (npot+1, nat)
      do 210 ip = 0, npot
        if (ipoint(ip).ne.0) then
          do 205 iat = 1,nmin
  205     if (iphx(iat).eq.ip) ilast = iat

          if (ilast.ne.ipoint(ip)) then
            xx  = xrat(1,ilast)
            yy  = xrat(2,ilast)
            zz  = xrat(3,ilast)

            xrat(1,ilast)= xrat(1,ipoint(ip))
            xrat(2,ilast)= xrat(2,ipoint(ip))
            xrat(3,ilast)= xrat(3,ipoint(ip))

            xrat(1,ipoint(ip)) = xx
            xrat(2,ipoint(ip)) = yy
            xrat(3,ipoint(ip)) = zz
c           now ipoint(ip) = ilast, but don't need ipoint anymore
          endif
        endif
  210 continue

c       if (idbg(4).eq.1) then
c           do 220 i=1,npot+1
c             print *,i,xrat(1,i),xrat(2,i),xrat(3,i),iphx(i)
c  220      continue
c       endif
      return
c  end subroutine sortat
      end
      subroutine xanlm(lmaxp1,mmaxp1)

c------------------------------------------------------------------
c  calculate and store all of the legendre polynomial normalization
c  factors needed in the problem
c     xnlm= sqrt ((2l+1)(l-m)!/(l+m)!)
c  see, for instance, Arfken section 12.6.  Note that this lacks the
c  factor of sqrt(4*pi)
c
c  inputs:
c     lmaxp1, nmaxp1:  maximun l and m considered in the problem +1
c                      i.e. lmaxp1 = l_max+1
c
c  outputs:
c     all normalization factors passed in common /xnlm/
c------------------------------------------------------------------

      implicit real(a-h,o-z)
c       parameter(ltot=6,mtot=3)

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
c={xparam.h
c nphasx MUST be the same as nphx, the maximum number of unique
c        potentials
c natxx MUST be the same as natx, the maximum number of atoms in the
c       extendeed cluster
c nexx MUST be the same as nex, the maximum number of energy points
      parameter (nphasx=nphx)
      parameter (natxx=natx)
      parameter (nexx=nex)
      parameter (istatx=(lx+1)**2*nclusx*nspx)
      parameter (nkmin=1)
c     parameter (nkmin=-9)
c= xparam.h}
c       include 'xstruc.h'
c====================================================================
c  This header file contains the structural information about the
c  cluster to be used for the full multiple scattering calculation

      common /xstruc/ xphi(nclusx,nclusx), xrat(3,nclusx),
     $            iphx(nclusx)
      save /xstruc/

c  xphi:  matrix of angles between the z axis and pairs of atoms
c  xrat:  xyz coordinates of the atoms in the cluster, the first
c         npot+1 entries are examples of each unique potential
c  iphx:  potential indeces of each atom in the cluster, ordered like
c         xrat
c********************************************************************
c**** rotation matrices for entire cluster
c
      complex drix
      common /rotx/ drix(-lx:lx, -lx:lx, 0:lx, 0:1, nclusx, nclusx)
      save /rotx/
c********************************************************************
c common blocks for saving rotation matrices between xanes and rotxan
       integer    jsavx, jsav, jbmagk
       parameter (jsavx = 150, roteps = 1.e-12,jbmagk=-9999)
       dimension drisav(-lx:lx,-lx:lx,0:lx,jsavx), betsav(jsavx)
       integer   ldsav(jsavx), mdsav(jsavx)
       common /rotsav/  drisav, betsav, ldsav, mdsav, jsav
       save  /rotsav/

c********************************************************************
c**** legendre polynomial normalization constants
c
      common /lnlm/ xnlm(0:lx,0:lx)
      save   /lnlm/

c********************************************************************
c**** legendre polynomial normalization constants
c
      common /xdwf/ sigsqr(nclusx,nclusx)
      save   /xdwf/

c  end of xstruc.h
c********************************************************************

      common/afctr/afac,flzero,flg(0:50)
c      common/afctr/afac,flzero,flg(0:210)
c      common/afctr/afac,flzero,flg(0:110) vax change

      call xfctst
      do 50 il=1,lmaxp1
        mmxp1 = min(mmaxp1,il)
        do 40 im=1,mmxp1
          l    = il-1
          m    = im-1
          cnlm = (2*l+1) * flg(l-m) / flg(l+m)
          cnlm = sqrt(cnlm) * afac**m
          xnlm(m,l) = cnlm

 40     continue
 50   continue
      return
c  end subroutine xlm
      end


      subroutine xfctst
c  same as feff's factst, but with a different name
      implicit real (a-h,o-z)
c     program for s3j and s6j symbols obtained from
c     steve younger of n.b.s.   modified by j.b.mann
c--------------------------------------------------------------------
c     a set to 1/64 to prevent overflow on vax
c     range on  flg set to 0:210, rather than flg(210)
c--------------------------------------------------------------------
cBR   This allows calculation of a large factorial (~100) without
cBR   overflow problems -- factor in a power of a small number then
cBR   factor it out
c--------------------------------------------------------------------
      common /afctr/ a, flzero, flg(0:50)
c      common /afctr/ a, flzero, flg(0:210)
      a=0.03125
c     a=0.015625
      flzero = 1.0
      flg(0) = 1.0
      flg(1) = a
      do 10 i=2,50
        flg(i) = flg(i-1) * i * a
 10   continue
      return
      end



c====================================================================
      subroutine istprm ( nph, nat, iphat, rat, iatph, xnatph,
     1                novr, iphovr, nnovr, rovr, folp, folpx, iafolp,
     1                edens, edenvl, idmag,
     2                dmag, vclap, vtot, vvalgs, imt, inrm, rmt, rnrm, 
     2                ixc, rhoint, vint, rs, xf, xmu, xmunew,
     3                rnrmav, qtotel, inters, totvol)

c     Finds interstitial parameters, rmt, vint, etc.
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

      dimension iphat(natx)
      dimension rat(3,natx)
      dimension iatph(0:nphx)
      dimension xnatph(0:nphx)
      dimension novr(0:nphx)
      dimension iphovr(novrx,0:nphx)
      dimension nnovr(novrx,0:nphx)
      dimension rovr(novrx,0:nphx)
      dimension folp(0:nphx), folpx(0:nphx)
      dimension edens(251,0:nphx), edenvl(251,0:nphx)
      dimension dmag(251,0:nphx+1)
      dimension vclap(251,0:nphx)
      dimension vtot (251,0:nphx), vvalgs (251,0:nphx)
      dimension imt(0:nphx)
      dimension inrm(0:nphx)
      dimension rmt(0:nphx)
      dimension rnrm(0:nphx)
      parameter (big = 5000)
      character*512 slog
      logical lnear
      dimension lnear(0:nphx), inn(0:nphx), rnnmin(0:nphx)
c#mn
       external dist

c     work space for linear algebra
      dimension ri(251)
      parameter (novp=40)
      complex cmovp(novp*(nphx+1)+1,novp*(nphx+1)+1)
      integer ipiv(novp*(nphx+1)+1)
      save lnear

c Find muffin tin radii.  We'll find rmt based on norman prescription,
c ie, rmt(i) = R * folp * rnrm(i) / (rnrm(i) + rnrm(j)),
c a simple average
c based on atoms i and j.  We average the rmt's from each pair of
c atoms, weighting by the volume of the lense shape formed by the
c overlap of the norman spheres.
c NB, if folp=1, muffin tins touch without overlap, folp>1 gives
c overlapping muffin tins.
c
c rnn is distance between sphere centers
c rnrm is the radius of the norman sphere
c xl_i is the distance to the plane containing the circle of the
c    intersection
c h_i  = rnrm_i - xl_i is the height of the ith atom's part of
c    the lense
c vol_i = (pi/3)*(h_i**2 * (3*rnrm_i - h_i))
c
c xl_i = (rnrm_i**2 - rnrm_j**2 + rnn**2) / (2*rnn)

c     find rmt from rnrm only on first call of istprm (rmt(0)=-1)
      if (rmt(0).le.0.0) then
      do 10 iph=0,nph
  10  lnear(iph)=.false.
      do 140  iph = 0, nph
         voltot = 0
         rmtavg = 0
         inrm(iph) = ii(rnrm(iph))
         if (novr(iph) .gt. 0)  then
c           Overlap explicitly defined by overlap card
            rnear = big
            inters = mod(inters,6)
c           use Norman prescription only in this case

            do 124  iovr = 1, novr(iph)
               rnn  = rovr(iovr,iph)
               inph = iphovr(iovr,iph)
               if (rnn .le. rnear) then
                  rnear = rnn
                  rnnmin(iph) = rnn
                  inn(iph) = inph
               endif
c              Don't avg if norman spheres don't overlap
               if (rnrm(iph)+rnrm(inph) .le. rnn)  goto 124
               voltmp = calcvl (rnrm(iph), rnrm(inph), rnn)
               voltmp = voltmp + calcvl (rnrm(inph), rnrm(iph), rnn)
               rmttmp = rnn * folp(iph) * rnrm(iph) /
     1                  (rnrm(iph) + rnrm(inph))
               ntmp = nnovr(iovr,iph)
               rmtavg = rmtavg + rmttmp*voltmp*ntmp
               voltot = voltot + voltmp*ntmp
  124       continue
         else
            iat = iatph(iph)
            rnear = big
            rmt(iph) = big
            do 130  inat = 1, nat
               if (inat .eq. iat)  goto 130
               rnn = dist (rat(1,inat), rat(1,iat))
               inph = iphat(inat)
               if (rnn .le. rnear) then
                  rnear = rnn
                  rnnmin(iph) = rnn
                  inn(iph) = inph
               endif
c              Don't avg if norman spheres don't overlap
               if (rnrm(iph)+rnrm(inph) .lt. rnn)  goto 130

               if (inters.lt.6) then
c                Norman prescription
                 voltmp = calcvl (rnrm(iph), rnrm(inph), rnn)
                 voltmp = voltmp + calcvl (rnrm(inph), rnrm(iph), rnn)
                 rmttmp = rnn * folp(iph) * rnrm(iph) /
     1                  (rnrm(iph) + rnrm(inph))
                 rmtavg = rmtavg + rmttmp*voltmp
                 voltot = voltot + voltmp
               else
c                Matching point prescription
                 do 125 i=inrm(iph),1,-1
                   j=ii(rnn-rnrm(iph))
                   if (vclap(i,iph).le.vclap(j,inph)) then
                     d1 = (vclap(i+1,iph)-vclap(i,iph))/(rr(i+1)-rr(i))
                     d2 =(vclap(j,inph)-vclap(j-1,inph))/(rr(j)-rr(j-1))
                     rmtavg = rr(i) + 
     1               (vclap(j,inph)+d2*(rnn-rr(i)-rr(j))-vclap(i,iph))
     2               /(d1+d2)
                     goto 127
c                    exit from the loop
                   endif
  125            continue
  127            continue
                 if (rmtavg.lt.rmt(iph)) rmt(iph) = rmtavg
               endif
  130       continue
         endif

c        special situation if rnrm is too close or larger than
c        the nearest neighbor distance
         if (rnrm(iph).ge.rnear) lnear(iph) = .true.

         if (rmtavg .le. 0)  then
            write(slog,132) iat, iph
            call wlog(slog)
  132       format (' WARNING: NO ATOMS CLOSE ENOUGH TO OVERLAP ATOM',
     1              i5, ',  UNIQUE POT', i5, '!!  ', 
     2              'Rmt set to Rnorman.  May be error in ',
     3              'input file.')
            rmt(iph) = rnrm(iph)
         elseif(inters.lt.6) then
c           Norman prescription
            rmt(iph) = rmtavg / voltot
            if (rmt(iph) .ge. rnear)  then
c              print*,iph, rmt(iph), rnear
               call wlog(' Rmt >= distance to nearest neighbor.  ' //
     1            'Not physically, meaningful.')
               call wlog(' FEFF may crash.  Look for error in ATOM '//
     1            'list or OVERLAP cards.')
            endif
            if (rnrm(iph) .ge. rnear) then
              imax = ii(rnear) - 1
c             begin until loop
 133            if (vclap(imax,iph).lt.vclap(imax+1,iph)) goto 134
                imax = imax-1
                goto 133
c             end of until loop
 134          continue
              rmt(iph) = exp(xx(imax)) - 0.0001
            endif
         endif

  140 continue

c     set maximum value for folp(iph) if AFOLP is in use
c     LMTO lore says no more than 15% overlap
c     do 144 iph = 0, nph
c 144 folpx(iph) = 1.15
c     already done in pot.f

      do 145 iph = 0, nph
         if (iafolp. gt. 0 ) then
            temp = 0.2 + 0.8 * rnrm(iph) / rmt(iph)
         else
            temp = 0.3 + 0.7 * rnrm(iph) / rmt(iph)
         endif
         if (temp.lt.folpx(iph)) folpx(iph) = temp
         temp = rnnmin(iph)/rmt(iph)/1.06d0
         if (temp.lt.folpx(iph)) folpx(iph) = temp
         temp = exp( -(novp-3)*0.05d0)
c      make sure that with given folpx(iph) the construction
c      of the overlapping matrix in movrlp will not fail
         if (lnear(iph)) then
c           lnear=.true. only when hydrogens are present in the system.
c           want to scale both rmt for iph and inn, so that overlapping
c           matrix calculations will not fail
            temp = rnnmin(iph) / (rmt(iph)*1.05d0 + temp*rmt(inn(iph)))
            if (temp.lt.folpx(iph)) folpx(iph) = temp
            if (temp.lt.folpx(inn(iph))) folpx(inn(iph)) = temp
         else
            temp = (rnnmin(iph) - rnrm(iph))/ (temp*rmt(inn(iph)))
            if (temp.lt.folpx(inn(iph))) folpx(inn(iph)) = temp
         endif
  145 continue

      endif
c     end of finding rmt from rnrm on first call of istprm.

c     Need potential with ground state xc, put it into vtot
      do 160  iph = 0, nph
         call sidx (edens(1,iph), 250, rmt(iph), rnrm(iph),
     1              imax, imt(iph), inrm(iph))
         do 150  i = 1, imax
            if (edens(i,iph).le.0) then
             if(mod(i,10).eq.0) then
               write(slog, 149) 'negative dens ', i,iph
  149          format (a, 2i3)
               call wlog(slog)
             endif
             rs = 100
             xmag=1.0
            else
              rs = (edens(i,iph)/3)**(-third)
c     spin dependent xc potential for ground state from Von Barth, Hedin
c     J.Phys.C:Solid State Phys., 5, 1629 (1972).
c     xmag/2 -fraction of spin up or down, depending on sign in renorm.f
c     put xmag = 1.0 to calculate cmd with external potential difference
              xmag = 1.0 + idmag*dmag(i,iph)
            endif
c           wrong for ferromagnets, need to overlap dmag(i)

c           vvbh from Von Barth Hedin paper, 1971
            call vbh(rs,xmag,vvbh)
            vtot(i,iph) = vclap(i,iph) + vvbh

            if (mod(ixc,10).eq.5) then
              rsval = 10.0
              if (edenvl(i,iph) .gt. 0.00001) 
     1           rsval = (edenvl(i,iph)/3)**(-third)
              if (rsval.gt.10.0) rsval = 10.0
              xmagvl = 1.0 + idmag * dmag(i,iph) 
     1                      * edens(i,iph) / edenvl(i,iph)
              call vbh(rsval,xmagvl,vvbhvl)
              vvalgs(i,iph) = vclap(i,iph) + vvbhvl
            elseif (mod(ixc,10) .ge. 6) then
              if (edens(i,iph).le.edenvl(i,iph)) then
                 rscore =101.0
              else
                 rscore = ((edens(i,iph)-edenvl(i,iph)) / 3)**(-third)
              endif
              rsmag = (edens(i,iph)*(1+idmag*dmag(i,iph)) / 3)**(-third)
              xfmag = fa/rsmag
              call edp(rscore,xfmag,vrdh)
              vvalgs(i,iph) = vclap(i,iph) + vvbh - vrdh
            else
              vvalgs(i,iph) = 0.d0
            endif
  150    continue
  160 continue

c     What to do about interstitial values?
c     Calculate'em for all atoms, print'em out for all unique pots along
c     with derivative quantities, like fermi energy, etc.
c     Interstitial values will be average over all atoms in problem.

c     rnrmav is averge norman radius,
c     (4pi/3)rnrmav**3 = (sum((4pi/3)rnrm(i)**3)/n, sum over all atoms
c     in problem
      rnrmav = 0
      xn = 0
c     volint is total interstitial volume
      volint = 0
      do 180  iph = 0, nph
         rnrmav = rnrmav + xnatph(iph) * rnrm(iph)**3
         volint=volint-xnatph(iph) * rmt(iph)**3
         xn = xn + xnatph(iph)
  180 continue
      if (totvol.le.0.0d0) then
         volint=4*pi/3 *(volint+rnrmav)
      else
         volint=4*pi/3 *volint + totvol
      endif
c     volume of lenses from overlapping mt spheres is added in movrlp.
      rnrmav = (rnrmav/xn) ** third

      rs = 0
      vint   = 0
      rhoint = 0
      rsval = 0

      call movrlp(nph, nat, iphat, rat, iatph, xnatph,
     1            novr, iphovr, nnovr, rovr,
     2            imt, rmt, rnrm, ri, lnear,
     3            cmovp,ipiv, volint,inters)

c     If no contribution to interstitial from any atom, die.
      if (volint .le. 0)  then
         call wlog(' No interstitial density.  Check input file.')
         call par_stop('ISTPRM')
      endif

c     find interstitial density

      call ovp2mt(nph, edens, 0, qtotel, ri, xnatph, lnear,
     1            inrm, imt, rnrm, rmt, cmovp,ipiv, rhoint,inters)
      rhoint = 4*pi * rhoint / volint

      if (ixc.ge.5) then
c        find valence potential inside mt sphere (vintvl -dummy)
         call ovp2mt(nph, vvalgs, 1, qtotel, ri, xnatph, lnear,
     1           inrm, imt, rnrm, rmt, cmovp, ipiv, vintvl,inters)
      endif

c     find potential inside mt sphere and vint
      call ovp2mt(nph, vtot, 1, qtotel, ri, xnatph, lnear,
     1            inrm, imt, rnrm, rmt, cmovp, ipiv, vint,inters)

      if (vint.ge.xmu) then
        write(slog,'(a)')
     1  ' WARNING:interstitial level found above Fermi level'
        call wlog(slog)
        write(slog,'(a)')
     1  '  Results may be unreliable. See manual for details'
        call wlog(slog)
        vint = xmu - 0.05d0
        call ovp2mt(nph, vtot, 2, qtotel, ri, xnatph, lnear,
     1            inrm, imt, rnrm, rmt, cmovp, ipiv, vint,inters)
      endif
      call fermi (rhoint, vint, xmunew, rs, xf)

      return
      end

      double precision function calcvl (r1, r2, r)
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
      xl = (r1**2 - r2**2 + r**2) / (2*r)
      h = r1 - xl
      calcvl = (pi/3) * h**2 * (3*r1 - h)
      return
      end
      subroutine movrlp ( nph, nat, iphat, rat, iatph, xnatph,
     1                novr, iphovr, nnovr, rovr,
     2                imt, rmt, rnrm, ri, lnear,
     3                cmovp, ipiv, volint, inters)

c     Constructs overlap matrix based on geometry of overlapped
c     muffin-tin spheres. Uses LU decomposition for inversion of matrix
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

      dimension iphat(natx)
      dimension rat(3,natx)
      dimension iatph(0:nphx)
      dimension xnatph(0:nphx)
      dimension novr(0:nphx)
      dimension iphovr(novrx,0:nphx)
      dimension nnovr(novrx,0:nphx)
      dimension rovr(novrx,0:nphx)
      dimension imt(0:nphx)
      dimension rmt(0:nphx)
      dimension rnrm(0:nphx)
      logical lnear
      dimension lnear(0:nphx)
c     local
      character*512 slog
c     work space for linear algebra
      dimension ri(251)
      parameter (novp=40)
      complex cmovp(novp*(nphx+1)+1,novp*(nphx+1)+1)
      real bmat(nphx+1,novp*(nphx+1))
      integer ipiv(novp*(nphx+1)+1)
c#mn
       external dist, ii

c     get ipot and irav from inters
      ipot = mod(inters,2)
      irav = (inters-ipot) / 2
      do 20 i=1,251
  20  ri(i)=exp(-8.8d0+(i-1)*0.05d0)
      exphx=exp(0.025d0)

c     initiallly cmovp is a unit matrix up to ncp
      ncp = novp*(nph+1)+1
      do 30 i2=1,ncp
      do 30 i1=1,ncp
        cmovp(i1,i2) = 0.d0
        if ( i1.eq.i2 ) cmovp(i1,i2) = 1.d0
        if (i2.eq.ncp) cmovp(i1,i2) = 0.01d0
  30  continue
      do 40 i2=1,ncp-1
      do 40 i1=1,nph+1
        bmat (i1,i2) = 0.d0
  40  continue
      xn = 0.d0

      do 200 ip1=0,nph
        if (novr(ip1) .gt. 0 ) then
           nlast = novr(ip1)
        else
           iat0 = iatph(ip1)
           ntmp = 1
           nlast = nat
        endif
        if (irav.eq.1) then
          rav = (rmt(ip1) + rnrm(ip1)) / 2
        elseif (irav.eq.0) then
          rav =  rnrm(ip1)
        else
          rav=ri(imt(ip1)+1)
        endif
        if (lnear(ip1)) rav=ri(imt(ip1)+1)

        do 190 iat = 1,nlast
          if (novr(ip1) .gt. 0 ) then
             ntmp = nnovr(iat,ip1)
             ip2 = iphovr(iat,ip1)
             rnn = rovr(iat,ip1)
          else
            if (iat.eq.iat0) goto 190
            ip2 = iphat(iat)
            rnn = dist (rat(1,iat0), rat(1,iat))
          endif

c         correct for double counting volume and area
          if (rnn .lt. rmt(ip1)+rmt(ip2)) then
c            correct interstitial volume
             volint = volint + xnatph(ip1) * ntmp *
     1       (calcvl( rmt(ip1), rmt(ip2), rnn) +
     2       calcvl(rmt(ip1), rmt(ip2), rnn)) / 2.d0
          endif

c         using expression for vtot(jri) ,(jri=i1)
c         first fill  matrix bmat
          ix1 = ip1+1

          if (rav+rmt(ip2) .le. rnn) goto 100
          imin2 = ii( rnn-rav )
          if (imt(ip2)-imin2 .ge. novp-1) then
             write(slog,132) ip1
  132        format(' FOLP for POTENTIAL type ',i3,' is too big.')
             call wlog (slog)
             write(slog,'(a)') ' Reduce overlap using FOLP and rerun'
             call wlog (slog)
             call par_stop('MOVRLP-1')
          endif
          imin2=imt(ip2)-novp+1

          do 80 i2 = imin2,imt(ip2)
             r1=ri(i2)/exphx
             r2=ri(i2)*exphx
             if (i2.eq.imt(ip2)) r2=rmt(ip2)
             if (i2.eq.imt(ip2))   r1=(r1+2*ri(imt(ip2))-rmt(ip2))/2.d0
             if (i2.eq.imt(ip2)-1) r2=(r2+2*ri(imt(ip2))-rmt(ip2))/2.d0
             if (r2+rav .lt. rnn) goto 80
             if (r1+rav .lt. rnn) then
c               use linear interpolation between cases xr=0, xr=1
                xr = (rnn-rav-r1)/ (r2-r1)
                r1 = rnn-rav   
                temp =  (r2**2 - r1**2) / (4*rnn*rav) * ntmp
                ind2=i2+1
                if (i2.eq.imt(ip2))  ind2=i2-1
                xr = xr * (r2-ri(i2)) / (ri(ind2)-ri(i2))
                ix2 = ip2*novp + i2 - imin2 + 1
                bmat (ix1,ix2) = bmat (ix1,ix2) + real(temp*(1-xr))
                ix2 = ip2*novp + ind2 - imin2 + 1
                bmat (ix1,ix2)=bmat (ix1,ix2) + real(temp*xr)
             else
                temp = (r2**2 - r1**2) / (4*rnn*rav   ) * ntmp
                ix2 = ip2*novp + i2 - imin2 + 1
                bmat (ix1,ix2) = bmat (ix1,ix2) + real( temp)
             endif
  80      continue

c         using expression for vtot(i) ,(i<jri)
c         construct matrix  cmovp
 100      if (rmt(ip1)+rmt(ip2) .le. rnn) goto 190

          imin1=ii(rnn-rmt(ip2))
          imin2=ii(rnn-rmt(ip1))
          if (imt(ip1)-imin1.ge.novp-1 .or. imt(ip2)-imin2.ge.novp-1) 
     1               call par_stop('tell authors to INCREASE NOVP')
          imin1=imt(ip1)-novp+1
          imin2=imt(ip2)-novp+1

          do 180 i1 = imin1,imt(ip1)
            ri1=ri(i1)/exphx
            ri2=ri(i1)*exphx
            if (i1.eq.imt(ip1)) ri2=rmt(ip1)
            if (i1.eq.imt(ip1)) ri1=(ri1+2*ri(imt(ip1))-rmt(ip1))/2.d0
            if (i1.eq.imt(ip1)-1)
     1                         ri2=(ri2+2*ri(imt(ip1))-rmt(ip1))/2.d0
            ix1 = i1-imin1+1  + ip1*novp
            do 170 i2 = imin2,imt(ip2)
              r1=ri(i2)/exphx
              r2=ri(i2)*exphx
              if (i2.eq.imt(ip2)) r2=rmt(ip2)
              if (i2.eq.imt(ip2))   r1=(r1+2*ri(imt(ip2))-rmt(ip2))/2.d0
              if (i2.eq.imt(ip2)-1) r2=(r2+2*ri(imt(ip2))-rmt(ip2))/2.d0
              if (r2+ri2.lt.rnn) goto 170

c             calculate volume of intersection
              temp = calcvl(ri2,r2,rnn) + calcvl(r2,ri2,rnn)
              if (ri1+r2.gt.rnn)
     1          temp = temp - calcvl(ri1,r2,rnn) - calcvl(r2,ri1,rnn)
              if (ri2+r1.gt.rnn)
     1          temp = temp - calcvl(ri2,r1,rnn) - calcvl(r1,ri2,rnn)
              if (ri1+r1.gt.rnn)
     1          temp = temp + calcvl(ri1,r1,rnn) + calcvl(r1,ri1,rnn)
c             volume of intersection (temp) should be devided by volume
c             volume between spheres ri1 and ri2
              temp=temp / ( 4.d0/3.d0*pi * (ri2**3-ri1**3) ) * ntmp

              if (r1+ri2.lt.rnn) then
c               use linear interpolation between cases xr=0, xr=1
                xr = (rnn-ri(i1)-r1)/ (r2-r1)

                ind2=i2+1
                if (i2.eq.imt(ip2))  ind2=i2-1
                xr = xr * (r2-ri(i2)) / (ri(ind2)-ri(i2))
                ix2 = i2-imin2+1 + ip2*novp
                cmovp(ix1,ix2)=cmovp(ix1,ix2) 
     1                              +cmplx (temp*(1-xr))
                ix2 = ind2-imin2+1 + ip2*novp
                cmovp(ix1,ix2)=cmovp(ix1,ix2) 
     1                               +cmplx (temp*xr)
                r1=rnn-ri2
              else
                ix1 = i1-imin1+1 + ip1*novp
                ix2 = i2-imin2+1 + ip2*novp
                cmovp(ix1,ix2)=cmovp(ix1,ix2)  +cmplx (temp)
              endif
 170        continue
 180      continue

 190     continue
         xn = xn + xnatph(ip1)
  200 continue

c     using matrix bmat fill in the last row of matrix cmvovp
c     this is additional equation to find Vint.
c     switch to local equation from average over all atoms
      if (ipot .eq. 0) then
         do 260 iph=0, nph
c          xn may differ from nat, if atom list have more natx atoms
c          see rdinp.f
           aa = xnatph(iph)/xn
           do 250 ix1 = 1, ncp-1
  250      cmovp(ncp,ix1) = cmovp(ncp, ix1) + aa*bmat(iph+1,ix1)
  260    continue
      else  
         iph=0
         do 270 ix1 = 1, ncp-1
  270    cmovp(ncp,ix1) = cmovp(ncp, ix1) + bmat(iph+1,ix1)
      endif

c --- invert matrices by LU decomposition
c     call cgetrf from lapack.  this performs an LU decomposition on
c     the matrix 
      istatx=novp*(nphx+1) + 1
      call cgetrf( ncp, ncp, cmovp, istatx, ipiv, info )
      if (info.ne.0) then
          call wlog('    *** Error in cgetrf when computing cmovp')
      endif

c     have to check that the last was not permuted, otherwise
c     the density calculation will be wrong
c     this is also why we put 0.01 in last column and not 1.0
      if (ipiv(ncp).ne.ncp) 
     .  call par_stop('illegal permutation in ipiv ')

      return
      end
      subroutine ovp2mt( nph, vtot, lrewr, qtot,ri,xnatph,lnear,
     1             inrm, imt, rnrm, rmt, cmovp, ipiv, vint, inters)
c  INPUT: nph - number of diferent potentials
c   vtot(i,iph) - potential OR density at point i for potential iph
c   lrewr       - if lrewr .gt. 0 potential will be overwritten
c                   density is never overwritten (lcoul.lt.0)
c                  lrewr=0 density calculation
c                  lrewr=1 potential calculation, vint estimated
c                  lrewr=2 potential calculation, vint is fixed
c   lcoul       -  .gt.0  (potential only) calculate charge for each iph
c                  .eq.0  (potential only) flat interstitial potential 
c                  .lt.0  (density only) calc charge inside MT spheres
c   qtot       -  for density only, total electron charge of cluster
c   ri         -  loucks radial grid
c   xnatph     -  number of atoms of type iph in the cluster
c   cmovp      -  LU decomposed overlapped matrix from movrlp.f
c   ipiv       -  pivoting indices for matrix cmovp
c  OUTPUT
c    vtot    if lrewr.gt.0  decomposed overlapped potential
c            if lrewr.le.0  old prescription for potential inside MT
c              spheres or don't want to overwrite densities
c    vint    mt zero level for potentials; charge outside mt spheres for
c            densities

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

      dimension vtot(251,0:nphx), xnatph(0:nphx)
      dimension inrm(0:nphx), imt(0:nphx), rmt(0:nphx), rnrm(0:nphx)
      dimension vtotav(0:nphx)
c     work space for linear algebra
      parameter (novp=40)
      complex cmovp(novp*(nphx+1)+1,novp*(nphx+1)+1)
      complex cvovp(novp*(nphx+1)+1)
      integer ipiv(novp*(nphx+1)+1)
      dimension  ri(251)
      character*13 trans
      dimension  crho(251)
      logical lnear
      dimension lnear(0:nphx)
cpot      character*30 fname

c      get ipot and irav from inters
      ipot = mod(inters,2)
      irav = (inters-ipot) / 2
c     prepare cvovp and bvec from vtot
      ncp=0
      do 25 ip1=0,nph
      do 25 i=1,novp
        ncp = ncp + 1
        ix1 = imt(ip1)-novp + i
        cvovp(ncp)= real( vtot(ix1,ip1) )
       if (lrewr.eq.2) cvovp(ncp) = cvovp(ncp) - vint
  25  continue
      do 27 ip1=0,nph
         if (irav .eq. 1) then
           rav = (rmt(ip1) + rnrm(ip1)) / 2
         elseif(irav.eq.0) then
           rav =  rnrm(ip1)
         else
           rav = ri(imt(ip1)+1)
         endif
         if (lnear(ip1)) rav = ri(imt(ip1)+1)
         call terp(ri,vtot(1,ip1),inrm(ip1)+2,3,rav,vtotav(ip1))
  27  continue
      istatx=novp*(nphx+1)+1
      trans = 'NotTransposed'
      nrhs = 1

c     find parameters for interstitial potential
      if (lrewr.gt.0) then
c        dealing with potentials
         if (lrewr.eq.1) then
c           additional equation to find vint
            ncp = ncp + 1
            cvovp(ncp) = 0
            bsum = 0
c           switch from average equation for vint to the local one
            nphlst = 0
            if (ipot .eq. 0) nphlst = nph
            do 430 iph=0,nphlst
               cvovp(ncp) = cvovp(ncp) + vtotav(iph)*xnatph(iph)
               bsum = bsum + xnatph(iph)
  430       continue
            cvovp(ncp) = cvovp(ncp) / bsum
         endif

         call cgetrs(trans, ncp, nrhs, cmovp, istatx,
     $               ipiv, cvovp, istatx, info)
         if (info.lt.0) then
             call par_stop('    *** Error in cgetrf')
c            stop
         endif

         if (lrewr.eq.1) vint = dble(real(cvovp(ncp))) /100.0

c        rewrite vtot
         do 550 iph=0,nph
 
cpot  to write out ovp tot pot and it's mt approxim, comment out cpot
cpot         write(fname,172)  iph
cpot  172    format('potp', i2.2, '.dat')
cpot         open (unit=1, file=fname, status='unknown', iostat=ios)
cpot         call chopen (ios, fname, 'wpot')

            do 500 i=1,novp
              index1=imt(iph)-novp + i
              index2=i+novp*iph

cpot            write(1,176) i, ri(index1), 
cpot     1             vtot(index1,iph),  dble(real(cvovp(index2)))+vint
cpot  176       format (1x, i4, 1p, 3e12.4)

              vtot(index1,iph) = dble(real(cvovp(index2)))+vint
  500       continue

cpot         close (unit=1)

c           use second order extrapolation
            j=imt(iph)+1
            call terp (ri,vtot(1,iph),imt(iph),2,ri(j),vtot(j,iph))
            do 505 j=imt(iph)+2, 251
  505       vtot(j,iph) = vint
  550    continue
      else
c        dealing with  density calculations. vint  is the total
c        charge inside mt spheres.
c        Divided by interstitial volume in istprm

         call cgetrs(trans, ncp, nrhs, cmovp, istatx,
     $            ipiv, cvovp, istatx, info)
         if (info.lt.0) then
             call par_stop('    *** Error in cgetrf')
c            stop
         endif

         vint = 0
         do 450 iph=0,nph
            do 440 i=1,imt(iph)+2
               if (i.lt.imt(iph)-novp+1) then
                 crho(i) =  vtot(i,iph)*ri(i)**2
               elseif (i.le. imt(iph)) then
                 ix1 = novp*iph +i-imt(iph)+novp
                 crho(i) = real(cvovp(ix1)) * ri(i)**2
c                crho(i) =  vtot(i,iph)*ri(i)**2
               else
                 call terp(ri,crho,imt(iph),2,ri(i), crho(i) )
               endif
  440       continue
            np = imt(iph) + 2
            cdum = 0
            dpas = 0.05d0
            call somm2 (ri,crho,dpas,cdum,rmt(iph),0,np)
            vint = vint + xnatph(iph) * cdum
  450    continue
         vint=qtot-vint
      endif

      return
      end
      subroutine fermi (rhoint, vint, xmu, rs, xf)

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

c     calculate fermi level of the system (mu) according to formula
c     mu=vcoulomb(interstitial)+vxc(interstitial)+kf(interstitial)^2
c     formula  2.13 in lee and beni, phys. rev. b15,2862(1977)

c     note that vint includes both coulomb and ground state
c     exchange-correlation potentials

c     den is the interstitial density
c     rs is the density parameter
c     xf is the interstital fermi momentum
c     xmu is the fermi level in hartrees

      den = rhoint / (4*pi)
      rs = (3 / (4*pi*den)) ** third
      xf = fa / rs
      xmu = vint + xf**2 / 2

      return
      end
      subroutine sidx (rholap, npts, rmt, rnrm, imax, imt, inrm)

      implicit double precision (a-h, o-z)
      dimension rholap (npts)
      character*512 slog
c#mn
      external ii, rr

      imt = ii (rmt)
      inrm = ii (rnrm)

c     Set imax (last non-zero rholap data)
      do 220  i = imt, npts
         if (rholap(i) .le. 1.0e-5)  goto 230
         imax = i
  220 continue
  230 continue

c     We need data up to the norman radius, so move norman
c     radius if density is zero inside rnrm.
      if (inrm .gt. imax)  then
         inrm = imax
         rnrm = rr (inrm)
  232    format(a,1pe13.5)
         write(slog,232) ' Moved rnrm.  New rnrm (au) ', rnrm
         call wlog(slog)
      endif
      if (imt .gt. imax)  then
         imt = imax
         rmt = rr (imt)
         write(slog,232) ' Moved rmt.  New rmt (au) ', rmt
         call wlog(slog)
      endif
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
      subroutine xcpot (iph, ie, index, lreal, ifirst, jri,
     1                  em, xmu,
     2                 vtot, vvalgs, densty, dmag, denval,
     3                  eref, v, vval, iPl, WpCorr, Gamma, AmpFac,
     4                  vxcrmu, vxcimu, gsrel, vvxcrm, vvxcim,rnrm)

      implicit double precision (a-h, o-z)
c     calculate self-energy correction
c     first coded j. mustre de leon
c     last modified a.ankudinov 1996 for non-local self-energies
c     Ankudinov, Rehr, J. Physique IV, vol. 7, C2-121, (1997).

c     INPUT
c     iph, ie used only for debug and labels.
c     index       0  Hedin-Lunqvist + const real & imag part
c                 1  Dirac-Hara + const real & imag part
c                 2  ground state + const real & imag part
c                 3  Dirac-Hara + HL imag part + const real & imag part
c                 4  See rdinp for comment
c     lreal       not equal zero for real self energy
c     ifirst      first entry flag, set to zero before first call for
c                 each unique potential, see vxcrmu and vxcimu below
c     jri         index of first interstitial point in current
c                 Loucks r grid
c     em          current energy grid point
c     xmu         fermi level
c     vi0         const imag part to subtract from potential
c     gamach      core hole lifetime
c     vtot(nr)    total potential (coulomb and gs exchange corr)
c     vvalgs(nr)  total coulomb + gs xc potential from valence electrons
c     densty(nr)  electron density
c     dmag(nr)    density magnetization
c     denval(nr)  valence electron density
c     iPl         Control for many pole self energy (Josh)
c
c     OUTPUT
c     eref        complex energy reference for current energy
c     v(nr)       complex potential including energy dep xc
c     vval(nr)    as above,but xc from valence electrons only
c     em          current energy
c
c     WORKSPACE
c     vxcrmu and vxcimu are calculated only on first entry for a
c     particular unique potential, re-used on subsequent entries.
c     vxcrmu(nr)  real part of xc at fermi level
c     vxcimu(nr)  imag part of xc at fermi level
c     gsrel(nr) ratio of gs xc potentials with and without magnetization
c     vvxcrm(nr)  real part of xc at fermi level from valence electrons
c     vvxcim(nr)  imag part of xc at fermi level from valence electrons


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

      dimension   vtot(nrptx), vvalgs(nrptx), densty(nrptx)
      dimension   dmag(nrptx), denval(nrptx)
      complex*16  em, eref, v(nrptx), vval(nrptx)
      dimension   vxcrmu(nrptx), vxcimu(nrptx)
      dimension   vvxcrm(nrptx), vvxcim(nrptx), gsrel(nrptx)

c     Josh added variables:
c     ZRnrm      - renormalization constant
c     ZTmp       - Temp var for ZRnrm
c     csig       - control for using many pole self energy
c     NRPts      - Number of points to use inside atom.
c                  Other points are linearly interpolated.
c     WpCorr     - Array of frequencies for many pole self energy.
c     rsTmp      - Temp var for rs
c     WpTmp      - Temp var for Wp
c     AmpFac     - g_i (pole strengths)
c     Gamma      - pole broadening
c     RsInt      - Rs in the intersitial
c     DRs        - RsCore - RsInt
c     delrHL     - Re[delta sigma] for many pole self energy
c     deliHL     - Im[delta sigma] for many pole self energy
c     Rs1(NRPts) - Array of Rs points for interpolation
      complex*16  delta, deltav, ZRnrm, ZTemp
      character*512 slog
      logical csig, rdmpse
      integer NRPts, iexist, lastPl
      parameter (tol=0.0004)
      parameter (NRPts=10)
      double precision WpCorr(MxPole), rsTmp, WpTmp, 
     &     AmpFac(MxPole), Gamma(MxPole)
      double precision RsInt, DRs, delrHL(NRPts), deliHL(NRPts),
     &     Rs1(NRPts), dx, x0, volume, totvol, rnrm, ri, riInt, omp,
     &     ompmax 
      complex*16 delavg
c    Josh END

c     First calculate vxc to correct the local momentum dispersion
c     relation, delta = vxc(e,k) - vxc(mu,k), and
c               p^2 = k^2 -mu + kf^2 - delta.
c     In jr theory, v(e,r) = vcoul(r) + vxc(e,r) =
c                          = vcoul(r) + vxcgs(r) + delta(e,r).

c    at jri potential is smooth continuation of potential to r(jri)
c    at this point potential jumps to interstitial value at jri+1
c     Atom r grid
      dx = 0.05d0
      x0 = 8.8d0
      totvol = 4.d0/3.d0*pi*rnrm**3
      delavg = 0.d0 
c      totvol = 0.d0
      ZRnrm = 0.d0
      csig=.false.
      ompm1 = 0.d0
      jri1 = jri + 1
      rsTmp = RsCorr
      nmax=1
      nul=0
      ibp = index / 10
      ixc = mod(index,10)
      ixcTmp=ixc
      DO i = 1, MxPole
         IF(WpCorr(i).le.0.d0) then
            lastPl = i-1
            GOTO 5
         END IF
      END DO
   5  CONTINUE        
      if((ixc.eq.0).and.(iPl.gt.0)) then
         csig=.true.
      end if
      if (ixc .eq. 2 .or. dble(em).le.xmu)  then
         do 10  i = 1, jri1
            v(i) = vtot(i)
            vval(i) = vvalgs(i)
   10    continue
c        Ground state exchange, no self energy calculation
         goto 888
      endif

c     Josh - Added CSigma to calculate HL Sigma with broadening and such.
c     Calculate Rs at the core and interstitial densities.
      if (densty(jri1).le.0) then
         RsInt =10
      else
         RsInt = (3 / (4*pi*densty(jri1))) ** third
      endif
      if (densty(1).le.0) then
         rscore =101.d0
      else
         rscore = (3 / (4*pi*densty(1))) ** third
      endif
      DRs = (RsInt-rscore)/(NRPts-1)
      omp = SQRT(3.d0/RsInt**3)*hart
      ompmax=omp*WpCorr(lastPl)
      PRINT*, omp, ompmax
c     Now calculate delta sigma as a function of Rs and energy
      if (csig) then
         do i= NRPts, 1, -1
            rdmpse = .false.
c            if((ifirst.eq.0).and.(i.eq.NRPts)) then
c               open(unit=23,file='mpse.bin',status='old',iostat=iexist)
c               if(iexist.eq.0) then
c                  rdmpse = .true.
c               else
c                  open(unit=23,file='mpse.bin',status='replace',
c     &                 iostat=iexist)
c               end if
c            end if
c            if(rdmpse) then
c               read(23,*) RsTmp, Rs1(i), delrHL(i), deliHL(i), ZTemp
c               if(abs(1-abs(RsTmp/DBLE(em))).gt.0.001) then
c                  goto 16
c               else
c                  goto 17
c               end if
c            end if
c 16         continue
            delrHL(i) = 0.d0
            deliHL(i) = 0.d0
            Rs1(i)=rscore+DBLE(i-1)*DRs
            
c           If iPl > 1, use renormalization, else not
c           If iPl = 2, use Sigma(r) = Sigma[Wp(r)*Wp/Wp(RsInt)]
c           If iPl = 3, use Sigma(r) = 0 outside of intersitial region
c                       actually linearly interpolates to zero at the 
c                       first RPt.
c           If iPl > 3, use Sigma(r) = Sigma(RsInt) (Sigma as a bulk property)
            if(iPl.gt.1) then
               if((iPl.eq.2).or.(i.eq.NRPts)) then                  
                  call CSigZ(em,xmu,Rs1(i),delrHL(i),deliHL(i),ZTemp,
     &                 WpCorr,Gamma,AmpFac)
               elseif(iPl.eq.3) then
                  delrHL(i) = 0.d0
                  deliHL(i) = 0.d0
               else
                  delrHL(i) = delrHL(NRPts)
               end if
               if(i.eq.NRPts) ZRnrm = ZTemp
            else
               call CSigma(em,xmu,Rs1(i),delrHL(i),deliHL(i),WpCorr,
     &              Gamma,AmpFac)
            end if
c            Josh Kas - Write self energy to mpse.bin for fast processing later
c            write(23,'(6f30.10)') DBLE(em), Rs1(i), delrHL(i),
c     &           deliHL(i), Ztemp
c 17         continue
c           debugging output of deltaSigma(em, rs)
c           write(44,'(6f30.10)') DBLE(em), Rs1(i), delrHL(i), deliHL(i),
c     &           dble(ZRnrm), dimag(ZRnrm)
         end do
c        write(44,*)
      end if
c     END Josh
      
c     Add the self energy correction
      do 20  i =  jri1,1,-1
         ri = exp((i-1)*dx - x0)
         if(i.eq.jri1) then
            riInt = ri
         end if
         niter = 0
         if (densty(i).le.0) then
            rs =10
         else
            rs = (3 / (4*pi*densty(i))) ** third
         endif
c         write(22,*) 1.d0*exp(dble(i)*0.01), densty(i)         
c        Josh - If csigma is turned on, interpolate onto rs.
c        Then skip to 15 (skip other calculations and self
c        consistency)
         if(csig) then          
            omp = SQRT(3.d0/rs**3)*hart
            if(iPl.ge.4) then
               delr = delrHL(NRPts)
               deli = deliHL(NRPts)
            else
               call terp (Rs1, delrHL, NRPts, 1, rs, delr)
               call terp (Rs1, deliHL, NRPts, 1, rs, deli)
            end if
            if((iPl.ne.5).or.(omp.lt.ompmax)) then
               goto 15
            end if
         end if
c        END Josh
         
c        xf = 1.9191.../rs
         xf = fa / rs
         rsm = rs / (1+dmag(i))**third
         xfm = fa / rsm

         if (ixc.eq.5) then
            if ( denval(i) .gt. 0.00001) then
               rsval = (3 / (4*pi*denval(i))) ** third
               if (rsval.gt.10.0) rsval=10.0
            else
               rsval = 10.0
            endif
            xfval = fa / rsval
         elseif (ixc.ge.6) then
            if (densty(i) .le. denval(i) ) then
               rscore = 101.0
            else
               rscore = (3 / (4*pi*(densty(i)-denval(i)))) ** third
            endif
         endif

         if (ifirst .eq. 0)  then
c           vxc_mu indep of energy, calc only once
c           Calculate vxc at fermi level e = mu, j.m. 1/12/89
            xk = xf * 1.00001
            gsrel(i) = 1.0d0
            if (ixc .lt. 5) then
              call sigma(ixc, ibp,rs,rscore,xk,vxcrmu(i),vxcimu(i))
              if (index .eq. 0) then
c  do not need 4 following lines for gs difference in potential
c                xmag = 1.0d0+ dmag(i)
c                call vbh(rs,xmag,v1)
c                call vbh(rs, 1.0d0,v0)
c                if (v0 .ne. 0) gsrel(i) = v1/v0
              endif
            else
              call sigma(nul,ibp, rs, rscore,xk,vxcrmu(i),vxcimu(i))
            endif
            if (ixc.eq.5 ) then
               xkpp = xfval * 1.00001
               call sigma 
     1         (ixc, ibp, rsval, rscore, xkpp, vvxcrm(i),vvxcim(i))
               if (ixc.eq.5 .and. i.eq.jri1) then
                  vvxcrm(jri1) =  vxcrmu(jri1)
                  vvxcim(jri1) =  vxcimu(jri1)
               endif
            elseif (ixc .ge. 6) then
               call sigma 
     1         (ixc, ibp, rs, rscore, xk, vvxcrm(i), vvxcim(i))
               if (ixc.eq.6 .and. i.eq.jri1) then
                  vvxcrm(jri1) =  vxcrmu(jri1)
                  vvxcim(jri1) =  vxcimu(jri1)
               endif
            else
               vvxcrm(i) = 0.0d0
               vvxcim(i) = 0.0d0
            endif
         endif

c        xk2 is the local momentum squared, p^2 = k^2 - 2*mu + kf^2,
c        k^2 represents energy measured from vacuum.
c        See formula 2.15 in Lee and Beni's paper with the last 2
c        terms neglected.  (complete reference?)
         xk2 = 2 * (dble(em) - xmu) + xf**2
         xk = sqrt(xk2)
         xkm2 = 2 * (dble(em) - xmu) + xfm**2
c        quick fix
         if (xkm2.lt.0) xkm2=xk2
         xkm = sqrt(xkm2)

c        find \delta_1
         if (ixc .lt. 5) then
            call sigma (ixc, ibp, rs, rscore, xk, vxcr, vxci)
         else
            call sigma (nul, ibp, rs, rscore, xk, vxcr, vxci)
         endif
         del1r = gsrel(i) * (vxcr - vxcrmu(i))

c        Correct local momentum according to the formula
c        p^2 = k^2 - 2*mu + kf^2 - 2*delta.  Note that imag part
c        of delta is ignored, since xk2 is a real quantity.

c        find xk(em) by iterative solution of dyson equation
  50     continue
         xk2 = 2*(dble(em) - xmu - del1r) + xf**2
         if (xk2 .lt. 0)  then
            write(slog,'(1pe13.5, 3i8, a)')
     1         xk2, i, ie, iph, ' xk2, i, ie, iph'
            call wlog(slog)
            call wlog(' em, xf**2, xmu, delta')
            write(slog,'(1p, 5e13.5)') dble(em), xf**2, xmu, del1r
            call wlog(slog)
            call par_stop('XCPOT-2')
         endif
         xk = sqrt (xk2)

c        calculate \delta_2 and \delta_v,2 with the corrected
c        local momentum
         call sigma (ixc, ibp, rs, rscore, xk, vxcr, vxci)
c        delta corrected calculated with new local momentum
         delr = gsrel(i) * (vxcr - vxcrmu(i))
         deli = vxci-vxcimu(i)

         if (ixc.ge.5 .and. i.eq.jri1 .and. xk.gt.xf) then
            if (ixc.eq.5 .or. ixc.eq.6) then
               delvr = delr
               delvi = deli
            endif
         endif

         if (niter.lt.nmax) then
            del1r=delr
            niter=niter+1
            go to 50
         endif

         if (ixc .ge. 5 .and. i.lt.jri1 .and. xk.gt.xf) then
            if (ixc.eq.5) then
               xkpp=sqrt(xk**2-xf**2+xfval**2)
               call sigma (ixc, ibp, rsval,rscore,xkpp,vxcvr,vxcvi)
            else
               call sigma (ixc, ibp, rs, rscore, xk, vxcvr, vxcvi)
            endif
            delvr = vxcvr-vvxcrm(i)
            delvi = vxcvi-vvxcim(i)
         endif

c        Josh - Skip SC loop if CSigma is called. CSigma calculates self consistently.
 15      continue
         
         delta = dcmplx(delr,deli)

c	 Josh - write out delta sigma at interstitial level to sigma.dat.
         if(i.eq.jri1) then
            write(45,'(X,20e14.6)') (DBLE(em) - xmu)*hart, delr*hart, 
     &                        deli*hart, DBLE(ZRnrm), DIMAG(ZRnrm), 
     &                        SQRT(DBLE(ZRnrm)**2+DIMAG(ZRnrm)**2), 
     &                        ATAN2(DIMAG(ZRnrm),DBLE(ZRnrm)),
     &                        SQRT(DBLE(em-xmu)/2.d0)/ABS(deli)*bohr
         end if
c	 Josh END

         if (ixc .eq. 5) delta = dcmplx(delr,delvi)
         v(i) = vtot(i) + delta
         if (ixc .ge. 5) then
            deltav = dcmplx(delvr,delvi)
            vval(i) = vvalgs(i) + deltav
         endif
         if(i.eq.jri1) then
            volume = 0.d0
         elseif(i.eq.jri) then
            volume = 4.d0*pi/3.d0*(rnrm**3 - exp(3.d0*((i-1)*dx - x0)))
         else
            volume = 4.d0*pi/3.d0*exp(3.d0*((i-1)*dx-x0))*(exp(3.d0*dx)
     &           - 1.d0)
         end if
         if(volume.lt.0.d0) volume = 0.d0
         omp = SQRT(3.d0/rs**3)
         write(39,'(I5,20f30.10)') i, dble(em-xmu), volume/totvol,
     &        exp((i-1)*dx - x0)*bohr,
     &        rnrm*bohr, densty(i), denval(i), dble(volume*delta),
     &        dimag(volume*delta), omp*hart, omp-ompm1
         ompm1 = omp
         delavg = delavg + volume*delta
 20   continue
 25   continue
      write(39,*)
      
      ifirst = 1
      delavg = delavg/totvol
      write(38,'(X,20e14.6)') (DBLE(em) - xmu)*hart, dble(delavg)*hart, 
     &     dimag(delavg)*hart,
     &     SQRT(DBLE(em-xmu)/2.d0)/ABS(dimag(delavg))*bohr,
     &     totvol
c     Reference the potential with respect to mt potential, ie,
c     first interstitial point.  v(jri1) = 0

c     Note that the reference does not contain the core hole lifetime
c     since the total atomic potential should have it. However in the
c     perturbation  deltav = v - vmt it cancels out.
c     ( deltav = vat - igamma - (vatmt-igamma) ).

 888  eref = v(jri1)
      do 910 i = 1, jri1
  910 v(i) = v(i) - eref
      if (ixc.ge.5) then
         do 920 i = 1, jri1
  920    vval(i) = vval(i) - eref
      else
         do 930 i = 1, jri1
  930    vval(i) = v(i)
      endif

c     Real self energy, zero imag part
      if (lreal.gt.0)  then
         do 950  i = 1, jri1
            v(i) = dble(v(i))
            if (ixc.gt.4)  vval(i) = dble(vval(i))
  950    continue
         eref = dble(eref)
      endif

      return
      end

      subroutine sigma (ixc, ibp, rs, rscore, xk, vr, vi)
      implicit double precision (a-h, o-z)

      if ((ixc.eq.0 .or. ixc.ge.5) .and. ibp .eq. 0) then
         call rhl (rs, xk, vr, vi)
      elseif ((ixc.eq.0.or. ixc.ge.5) .and. ibp .eq. 1) then
         call rhlbp (rs, xk, vr, vi)
      elseif (ixc .eq. 1) then
         vi = 0
         call edp(rs,xk,vr)
      elseif (ixc .eq. 3) then
         call edp(rs,xk,vr)
         call imhl (rs,xk,vi,icusp)
      endif

      if (ixc .ge. 6) then
         call edp(rscore,xk,vrp)
         vr = vr - vrp
      endif

      return
      end

      subroutine cubic (xk0, wp, alph, rad, qplus, qminus)

c     input:  xk0, wp, alph
c     output: rad, qplus, qminus

      implicit double precision (a-h, o-z)
      complex*16 s1,s13
      parameter (three = 3)
      parameter (third = 1/three)

c     this subroutine finds the roots of the equation
c     4xk0 * q^3  +  (alph-4xk0^2) * q^2  +  wp^2 = 0
c     see abramowitz and stegun pg 17 for formulae.

      a2 = (alph / (4*xk0**2)  -  1) * xk0
      a0 = wp**2 / (4*xk0)
      a1 = 0
      q = a1/3 - a2**2/9
      r = (a1*a2 - 3*a0)/6  -  a2**3/27
      rad = q**3 + r**2
      if (rad .gt. 0) then
         qplus = 0
         qminus = 0
         return
      endif

      s13 = dcmplx (r, sqrt(-rad))
      s1 = s13 ** third
      qz1 = 2*s1 - a2/3
c     qz2 = -(s1 + sqrt(three)*dimag(s1) + a2/3)
      qz3 = -(s1 - sqrt(three)*dimag(s1) + a2/3)
      qplus = qz1
      qminus = qz3

      return
      end
c***********************************************************************
c
c     this subroutine calculates the ' energy dependent
c     exchange-correlation potential' (or 'dirac- hara potential')
c     ref.: paper by s.h.chou, j.j.rehr, e.a.stern, e.r.davidson (1986)
c
c     inputs:    rs in a.u.
c                xk momentum in a.u.
c     outputs:   vr --- dirac potential (Hartrees)
c     written by j. mustre 8/31/87
c**********************************************************************

      subroutine edp (rs, xk, vr)
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

      vr = 0.0d0
      if (rs .le. 100.0) then
c       p = sqrt (k^2 + kf^2) is the local momentum, and x = p / kf
c       Reference formula 23 in Role of Inelastic effects in EXAFS
c       by Rehr and Chou. EXAFS1 conference editted by Bianconi.
c       x is local momentum in units of fermi momentum

        xf = fa / rs
        x = xk / xf
        x = x + 1.0e-5
c       set to fermi level if below fermi level
        if (x .lt. 1.00001) x = 1.00001
        c = abs( (1+x) / (1-x) )
        c = log(c)
        vr = - (xf/pi) * (1 + c * (1-x**2) / (2*x))
      endif

      return
      end
      double precision function ffq (q, ef, xk, wp, alph)
      implicit double precision (a-h,o-z)

c     input:  q, wp, alph, ef, xk
c             q is dimensionless, normalized to fermi momentum
c             xk is momentum in invBohrs
c     output: ffq only

      wq = sqrt (wp**2 + alph*q**2 + q**4)
      ffq = (wp+wq)/(q**2) + alph/(2*wp)
      ffq = ((ef*wp) / (4*xk))  * log(ffq)

      return
      end
      subroutine imhl (rs, xk, eim, icusp)
      implicit double precision (a-h,o-z)

c     what is xk?  k**2 - mu + kf**2?

c written by j. mustre (march 1988)
c code is based on analytical expression derived by john rehr.
c it leaves the real part, calculated in rhl unchanged.
c
c modified by j. rehr  (oct 1991) - adds quinn approximation for
c losses due to electron-hole pairs below the plasmon turn on
c see new subroutine quinn.f, which incorporates r. albers coding of
c j.j. quinn's approximations for details.

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
c     alph is Hedin-Lundquist parameter
      parameter (alph = 4.0 / 3.0)
      external ffq

      icusp=0
      xf = fa / rs
      ef = xf**2 / 2

c     xk0 is xk normalized by k fermi.
      xk0 = xk/xf
c     set to fermi level if below fermi level
      if (xk0 .lt. 1.00001) then
         xk0 = 1.00001
      endif

c     wp is given in units of the fermi energy in the formula below.
      wp = sqrt (3 / rs**3) / ef
      xs = wp**2 - (xk0**2 - 1)**2

      eim = 0
      if (xs .lt. 0.)  then
         q2 = sqrt ( (sqrt(alph**2-4*xs) - alph) / 2 )
         qu = min (q2, (1+xk0))
         d1 = qu - (xk0 - 1)
         if (d1 .gt. 0)  then
            eim = ffq (qu,ef,xk,wp,alph) - ffq (xk0-1,ef,xk,wp,alph)
         endif
      endif
      call cubic (xk0, wp, alph, rad, qplus, qminus)

      if (rad .le. 0) then
         d2 = qplus - (xk0 + 1)
         if (d2 .gt. 0)  then
            eim = eim + ffq (qplus,ef,xk,wp,alph) - 
     1                  ffq (xk0+1,ef,xk,wp,alph)
         endif
         d3 = (xk0-1) - qminus
         if (d3 .gt. 0)  then
            eim = eim + ffq (xk0-1,ef,xk,wp,alph) - 
     1                  ffq (qminus,ef,xk,wp,alph)
c           beginning of the imaginary part and position of the cusp x0
            icusp = 1
         endif
      endif

      call quinn (xk0, rs, wp, ef, ei)
      if (eim .ge. ei)  eim = ei

      return
      end
      subroutine quinn (x, rs, wp, ef, ei)
      implicit double precision (a-h, o-z)

c     input  x, rs, wp, ef
c     output ei

c***********************************************************************
c
c     quinn: calculates low energy gamma (approx. proportional to e**2)
c             formula taken from john j. quinn, phys. rev. 126,
c             1453 (1962); equation (7).
c             a cut-off is set up at quinn's cutoff + ef = ekc; it is a
c             rounded inverted step function (a fermi function)
c             theta = 1/( 1 + exp((e-ekc)/gam)) )
c             where the rounding factor gam is set to be about 0.3 ekc.
c     modified by j. rehr (oct 1991) based on coding of r. albers
c     subroutines quinn.f and quinnc.f
c
c     variables:
c        x  = p/pf
c        rs = ws density parameter
c        ei = imaginary self energy
c        pfqryd = quinn's prefactor in atomic-rydberg units
c        wkc = quinn's plasmon threshold
c
c***********************************************************************

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

      parameter (alphaq = 1/ fa)

c     calculate quinn prefactor in atomin Hartree units
      pisqrt = sqrt(pi)
      pfq = pisqrt / (32 * (alphaq*rs)**1.5)
      temp1 = atan (sqrt (pi / (alphaq*rs)))
      temp2 = sqrt(alphaq*rs/pi) / (1 + alphaq*rs/pi)
      pfq = pfq * (temp1 + temp2)

c     calculate quinn cutoff
c     wkc = quinn's plasmon threshold
c     wkc is cut-off of quinn, pr126, 1453, 1962, eq. (11)
c     in formulae below wp=omegap/ef
      wkc = (sqrt(1+wp) - 1)**2
      wkc = (1 + (6./5.) * wkc / wp**2) * wp * ef

c     we add fermi energy to get correct energy for
c     plasma excitations to turn on
      ekc = wkc + ef

c     calculate gamma
c     gamryd = 2 * (pfqryd/x) * (x**2-1)**2
      gam = (pfq/x) * (x**2-1)**2

c     put in fermi function cutoff
      eabs = ef * x**2
      arg = (eabs-ekc) / (0.3*ekc)
      f = 0
      if (arg .lt. 80)  f = 1 / (1 + exp(arg))

      ei = -gam * f / 2

      return
      end
      subroutine rhl (rs, xk, erl, eim)
      implicit double precision (a-h, o-z)

c     input:  rs, xk
c     output: erl, eim

c     This is a new hl subroutine, using interpolation for the
c     real part while the imaginary part is calculated analytically.
c     It uses hl to calculate values at the mesh points for the inter-
c     polation of the real part. The imaginary part is calculated
c     using subroutine imhl.
c
c     written by jose mustre
c     polynomial in rs has a 3/2 power term. j.m.


c     for the right branch the interpolation has the form:
c     hl(rs,x) = e/x + f/x**2 + g/x**3
c     where e is known and
c        f = sum (i=1,3) ff(i) rs**(i+1)/2
c        g = sum (i=1,3) gg(i) rs**(i+1)/2
c
c
c     lrs=number of rs panels, in this case one has 4 panels
c     nrs=number of standard rs values, also order of rs expansion
c     if you change nrs you need to change the expansion of hl
c     in powers of rs that only has 3 terms!
c     nleft=number of coefficients for x<x0
c     nright=number of coefficients for x>x0

      parameter (lrs=4, nrs=3, nleft=4, nright=2)
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

      dimension cleft(nleft), cright(nright)

      dimension rcfl(lrs,nrs,nleft), rcfr(lrs,nrs,nright)
      data rcfr/-0.173963d+00,-0.173678d+00,-0.142040d+00,-0.101030d+00,
     1     -0.838843d-01,-0.807046d-01,-0.135577d+00,-0.177556d+00,
     2     -0.645803d-01,-0.731172d-01,-0.498823d-01,-0.393108d-01,
     3     -0.116431d+00,-0.909300d-01,-0.886979d-01,-0.702319d-01,
     4      0.791051d-01,-0.359401d-01,-0.379584d-01,-0.419807d-01,
     5     -0.628162d-01, 0.669257d-01, 0.667119d-01, 0.648175d-01/
      data rcfl/ 0.590195d+02, 0.478860d+01, 0.812813d+00, 0.191145d+00,
     1     -0.291180d+03,-0.926539d+01,-0.858348d+00,-0.246947d+00,
     2      0.363830d+03, 0.460433d+01, 0.173067d+00, 0.239738d-01,
     3     -0.181726d+03,-0.169709d+02,-0.409425d+01,-0.173077d+01,
     4      0.886023d+03, 0.301808d+02, 0.305836d+01, 0.743167d+00,
     5     -0.110486d+04,-0.149086d+02,-0.662794d+00,-0.100106d+00,
     6      0.184417d+03, 0.180204d+02, 0.450425d+01, 0.184349d+01,
     7     -0.895807d+03,-0.318696d+02,-0.345827d+01,-0.855367d+00,
     8      0.111549d+04, 0.156448d+02, 0.749582d+00, 0.117680d+00,
     9     -0.620411d+02,-0.616427d+01,-0.153874d+01,-0.609114d+00,
     1      0.300946d+03, 0.109158d+02, 0.120028d+01, 0.290985d+00,
     2      -0.374494d+03,-0.535127d+01,-0.261260d+00,-0.405337d-01/

c
c     calculate hl using interpolation coefficients
      rkf = fa/rs
      ef  = rkf**2/2
      wp  = sqrt (3/rs**3)
c    quick fix to remove jump at wp in rhl. ala 08.01.95
c    use smooth transition between 2 curves in energy range dwp
      dwp = wp/3.0

      call imhl (rs, xk, eim, icusp)

c     eim already has a factor of ef in it j.m.
c     eim also gives the position of the cusp

      xx = xk / rkf
c     set to fermi level if below fermi level
      if (xx .lt. 1.00001) then
          xx = 1.00001
      endif
c    quick fix to remove jump at wp in rhl. ala 08.01.95
      deltae = ((xx**2-1.0)*ef - wp-dwp)/dwp

c     calculate right hand side coefficients
      if (rs .lt. 0.2) then
         mrs=1
      elseif (rs .lt. 1.0) then
         mrs=2
      elseif (rs .lt. 5.0) then
         mrs=3
      else
         mrs=4
      endif

      do 210 j=1,nright
         cright(j) = rcfr(mrs,1,j)*rs + rcfr(mrs,2,j)*rs*sqrt(rs)
     1               + rcfr(mrs,3,j)*rs**2
  210 continue
      eee=-pi*wp/(4*rkf*ef)

c     if (icusp .ne. 1) then
c    quick fix to remove jump at wp in rhl. ala 08.01.95
      if (icusp .ne. 1 .or. abs(deltae).lt.1.0) then

         do 230 j=1,nleft
            cleft(j) = rcfl(mrs,1,j)*rs + rcfl(mrs,2,j)*rs**1.5
     1                 + rcfl(mrs,3,j)*rs**2
  230    continue
         erl=cleft(1)
         do 250 j=2,nleft
            erl=erl+cleft(j)*xx**(j-1)
  250    continue

c     else
c    quick fix to remove jump at wp in rhl. ala 08.01.95
      endif
      if(icusp .eq. 1 .or. abs(deltae).lt.1.0) then
c        right branch
         erlr=eee/xx
         do 280 j=1,nright
            erlr=erlr+cright(j)/xx**(j+1)
  280    continue
         if (abs(deltae).lt.1.0) then
            if (deltae.lt.0) then
               wr = (1.0 + deltae)**2/2.0
            else
               wr = 1.0 - (1.0-deltae)**2/2.0
            endif
            erl=wr*erlr + (1.0-wr)*erl
         else
            erl= erlr
         endif
      endif

      erl = erl * ef

      return
      end
      subroutine rhlbp (rs, xk, erl, eim)
c     This is a new broadened plasmon hl subroutine, 
c     using interpolation for the real and imaginary part.
c     test of multi-pole pole model
c     input:  
c        rs - r_s 
c        xk - k in a.u.
c     output: 
c        erl, eim - Re and Im part of self energy normalized to k_f**2/2

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
      parameter (nrs=21, nx=51 )
      dimension rsmesh(nrs), xmesh(nx), sigma(nrs,nx,2)
      save ifirst, rsmesh, xmesh, sigma
      data  ifirst /0/

      xf = fa / rs
      ef = xf *xf / 2.  
      wp = sqrt (3 / rs**3) / ef
      xk0 = xk / xf
      xx = (xk0 ** 2  - 1) / sqrt(rs)

      if (ifirst .eq. 0) then
c        read self energy for grid points from bphl.dat
         open (unit=2, file='bphl.dat', status='old', iostat=ios)
         call chopen (ios, 'bphl.dat', 'rhlbp')
         xmesh(1) = 0.0
         do 200 irs = 1, nrs
            sigma (irs, 1, 1) = 0.0
            sigma (irs, 1, 2) = 0.0
c           irs correspond to grid in r_s: rs = 10.0**(0.1 * irs)
            do  100 ik = 2, nx
c              xmesh define grid in k-space as follows:
c              xmesh = ((ik-1) / 20.0) * (1.0 + ((ik-1) / 20.0)**4.0)
c              xmesh = (xk**2 - 1) / sqrt (rs)
c              xk = sqrt (xmesh * sqrt(rs) + 1.0)
c              xk = k / k_f
               read(2, *) rsmesh(irs), xmesh(ik), 
     1              sigma(irs, ik, 1), sigma(irs, ik, 2)
 100        continue
 200     continue
         ifirst = 1
         close (unit=2)
      endif

c     delev = xdel * ef * hart * rs
      call terp2d (rsmesh, xmesh, sigma(1, 1, 1), nrs, nx, rs, xx, erl)
      call terp2d (rsmesh, xmesh, sigma(1, 1, 2), nrs, nx, rs, xx, eim)
c     transfer to atomic units
      erl = erl / rs / hart
      eim = eim / rs / hart

      call quinn (xk0, rs, wp, ef, ei)
      if (eim .ge. ei)  eim = ei

      return
      end

      subroutine terp2d (x, y, z, nx, ny, x0, y0, z0)
c     Linear interpolation and extrapolation.
c     2d analog of terp.f
      implicit double precision (a-h, o-z)

      dimension x(nx), y(ny), z(nx,ny)

c     Find out between which x points x0 lies
      ix = locat (x0, nx, x)
c     if i < 1, set i=1, if i > n-1, set i=n-1
      ix = max (ix, 1)
      ix = min (ix, nx-1)
      if (x(ix+1) - x(ix) .eq. 0)  call par_stop('TERP-1')
c     Find out between which y points y0 lies
      iy = locat (y0, ny, y)
c     if i < 1, set i=1, if i > n-1, set i=n-1
      iy = max (iy, 1)
      iy = min (iy, ny-1)
      if (y(iy+1) - y(iy) .eq. 0)  call par_stop('TERP-1')

      dx = (x0 - x(ix)) / (x(ix+1) - x(ix))
      dy = (y0 - y(iy)) / (y(iy+1) - y(iy))
      z1 = z(ix,iy) +  dx * (z(ix+1,iy) - z(ix,iy))
      z2 = z(ix,iy) +  dx * (z(ix+1,iy) - z(ix,iy))
      z0 = z1 + dy * (z2 - z1)

      return
      end
      subroutine vbh(rs,xmag,vxc)

      implicit double precision (a-h, o-z)

c   INPUT: density parameter rs, 2* fraction of given spin orientation.
c   OUTPUT: xc potential for given spin orientation.
c   Reference: Von Barth, Hedin, J.Phys.C, 5, 1629, (1972). eq.6.2
c   xmag is twice larger than 'x' in their paper
c   effect of tau was also found to be small. thus tau is not used

c     parameter (asm = 2.0**(-1.0/3.0) )
c     parameter (gamma = 4.0/3.0*asm/(1-asm) )
c APS parameter (gamma = 5.129762496709890 ) changed
      parameter (gamma = 5.129762802484097 )

      vxc = 0.0
      if (rs.gt.1000) goto 999
      epc = -0.0504 * flarge(rs/30)
      efc = -0.0254 * flarge(rs/75)
      xmup = -0.0504*log(1.0+30.0/rs)
c     xmuf = -0.0254*log(1.0+75.0/rs)
      vu = gamma*(efc - epc)
c     tau = xmuf-xmup-(efc-epc)*4.0/3.0
     
      alg = -1.22177412/rs + vu
      blg = xmup - vu
      vxc = alg*xmag**(1.0/3.0) + blg
c     vxc = alg*xmag**(1.0/3.0) + blg +tau*fsmall(xmag/2.0)

 999  continue
c     transform to code units (Hartrees) from Rydbergs
      vxc = vxc / 2.d0

      return
      end

      double precision function flarge(x)
      implicit double precision (a-h, o-z)
        flarge = (1+x**3)*log(1+1/x) + x/2 - x**2 - 1.0/3.0
      return
      end

c     double precision function fsmall(x)
c     implicit double precision (a-h, o-z)
c     parameter (a = 2.0**(-1.0/3.0) )
c       fsmall = ( x**(4/3) + (1.0-x)**(4/3) - a ) / (1.0-a)
c     return
c     end
      SUBROUTINE CSigma(Energy, Mu, Rs, ReSig, ImSig, WpScl, Gamma,
     &     AmpFac)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Written by Josh Kas
c     This subroutine calculates the self energy Sigma(k(Energy),Energy)
c     based on an electron gas model of epsilon^-1.
c      
c     Solve: k0**2 = 2*Energy - 2*Mu -
c                    2*(Sigma(k0,Energy)-Sigma(kFermi,EFermi))
c            
c     Steps:
c
c            1. k0**2  = 2*(Energy-Mu) + SigmaF (SigmaF is self energy at fermi level).
c            2. Sigma0 = Sigma(k0,Energy)
c            3. k1**2  = 
c                  k0**2 - 2*(Sigma0-SigmaF)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
c     Parameters:
c     MxPole - Maximum number of poles      
      INTEGER MxPole
      PARAMETER(MxPole=1000)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Input:
c     Energy - Energy at which to evaluate Sigma
c     Mu     - Fermi energy as calculated by FEFF
c     Rs     - R sub s (sphere of radius Rs holds charge e)
c     WpScl  - Scale Wp in interstitial by WpScl
c     Gamma  - Use broadening Gamma when calculating Sigma
c     AmpFac - Use amplitude AmpFac for plasmon pole.
c     Note: Atomic units are used.
      DOUBLE PRECISION Rs, WpScl(MxPole), Gamma(MxPole), AmpFac(MxPole),
     &     Mu
      COMPLEX*16 Energy
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Output:
c     ReSig  - Re[Sigma(Energy,k(Energy))]
c     ImSig  - Im[Sigma(Energy,k(Energy))]
      DOUBLE PRECISION ReSig, ImSig
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Local Variables:
c     kFermi - Fermi momentum calculated from Rs via electron gas
c              approximation.
c     EFermi - Fermi energy = kFermi^2/2
c     Wp     - Electron gas plasmon frequency.
c     Gam    - broadening for broadened pole.
c     ckF    - complex variable to store kFermi
c     ck0    - complex momentum
c     SigmaF - Self energy at the fermi energy and fermi momentum
c              Does not include the Hartree Fock part.
c     Sigma0 - Single pole self energy evaluated at
c              Wp = Wp(ipole)/Wp(R_Interstitial)*Wp(Rs)
c     dSgdE  - Derivative of sigma w.r.t. Energy
c     ZTot   - Renormalization factor Z = 1/(1-dSgdE)
c     RelEn  - Energy relative to the fermi energy from FEFF      
c              Relen = Energy - Mu + EFermi
c     SigTot - Total many pole deltaSigma = Sigma(E,k(E))-Sigma(EF,kF)
c     DelHF  - Sigma_HartreeFock(k) - Sigma_HartreeFock(kF)      
      DOUBLE PRECISION kFermi, EFermi, Wp, Gam
      COMPLEX*16 ckF, ck0, Sig, SigmaF, Sigma0,
     &     RelEn, SigTot, DelHF
      INTEGER i1, i2

c     Parameters:
      DOUBLE PRECISION DPZero, h
      PARAMETER(DPZero = 0.d0, h = 1.d-3)
      INTEGER MxIter
      LOGICAL MPole
      PARAMETER(MxIter = 1)

c     Externals:
      COMPLEX*16 Sigma1, dSigma, HFExc
      EXTERNAL Sigma1, dSigma, HFExc
      
c     Initialization
      kFermi = fa/Rs
      EFermi = kFermi*kFermi/2.d0
      SigTot=0.d0
      SigmaF = 0.d0 
      Gam = 0.d0

c     Loop1: Start self consistency loop.
      DO i2 = 1, MxIter
c        Loop2: Loop over poles to find SigmaF
         DO i1 = 1, MxPole
            IF(WpScl(i1).lt.-1000.d0) GOTO 5
            
c           Wp is in Hartrees
            Wp = SQRT(3.d0/rs**3)*WpScl(i1)
            
c           find Sigma_Fermi (SigmaF)
            ckF = kFermi*1.00001d0
            RelEn = EFermi
            SigmaF = SigmaF + Sigma1(ckF,RelEn,Wp,Gam,AmpFac(i1),
     &           kFermi,EFermi)
         END DO
 5       CONTINUE
         
c        Loop3: Loop over poles
         DO i1 = 1, MxPole
            IF(WpScl(i1).lt.-1000.d0) GOTO 10
c           Wp is in Hartrees
            Wp = SQRT(3.d0/rs**3)*WpScl(i1)

c           Start with ck0=Sqrt[Re(Energy)-Mu+EFermi]
            RelEn = DBLE(Energy) - Mu + EFermi
            ck0 = SQRT(2.d0*DBLE(RelEn))
            
c           Find Sigma0
            Sigma0 = Sigma1(ck0,RelEn,Wp,Gam,AmpFac(i1),kFermi,
     &           EFermi)
            
            SigTot = SigTot + Sigma0
            
c        End loop over poles.            
         END DO
 10      CONTINUE
         
c     End self-consistency loop
      END DO

c     Form delta sigma and retur.n
      SigTot = SigTot - SigmaF
      DelHF = HFExc(ck0,EFermi,kFermi) - HFExc(ckF,EFermi,kFermi)
      SigTot = SigTot + DelHF
      
      ReSig = DBLE(SigTot)
      ImSig = DIMAG(SigTot)
      
      RETURN
      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      FUNCTION Sigma1(ck,Energy,Wi,Gamma,Amp,kFermi,EFermi)
c     Written by Josh Kas
c     Function Sigma calculates the energy dependent part
c     of Sigma(ck,Energy) from H.L. electron gas model.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
c     Input:
c     ck     - complex momentum
c     Energy - Energy
c     Wi     - Plasmon pole energy
c     Gamma  - Broadening of plasmon pole.
c     Amp    - Amplitude of plasmon pole.
c     kFermi - Fermi momentum
c     EFermi - Fermi energy
c              This is used when calculating dSigma/dE
      DOUBLE PRECISION Wi, Gamma, Amp, kFermi, EFermi
      COMPLEX*16 ck, Energy
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Output:
c     Sigma(ck,Energy)
      COMPLEX*16 Sigma1
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Local Variables:
c     NSing  - Number of singularities of integrand (used in cgratr)
c     NCalls - Number of fcn calls used in cgratr
c     MaxNR  - Max number of regions used in cgratr
c     DPPar  - Array of double precision parameters passed to cgratr
c              to be used in the functions r1, r2, and r3
c     CPar   - Array of complex parameters passed to cgratr
c              to be used in the functions r1, r2, and r3
c     Limit1 - Lower limit of integration
c     Limit2 - Upper limit of integration
c     HLInt1 - Integral of r2 (first integral in eq. 13 of H.L.)
c     HLInt2 - Integral of r1 (second integral in eq. 13 of H.L.)
c     HLInt3 - Integral of r1 or r3 (3rd or 4th integral)
c     XSing  - Array of singularities of integrand (used by cgratr)      
      INTEGER NSing, NCalls, MaxNR
      DOUBLE PRECISION DPPar(10), Wp, Beta
      COMPLEX*16 CPar(10), Limit1, Limit2, HLInt1, HLInt2, HLInt3,
     &     XSing(20)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Loop variables:
      INTEGER i1, i2
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Parameters:
c     ZeroPl - lower bound Limit1
c     Inf    - Upper bound of Limit2
c     AbsErr - absolute error used by cgratr
c     RelErr - Relative error used by cgratr
c     Error  - used for error codes by cgratr      
      DOUBLE PRECISION ZeroPl, Inf, AbsErr, RelErr, Error
      PARAMETER(ZeroPl = 1.d-5, Inf = 1.d2)
      PARAMETER(AbsErr = 1.d-5, RelErr = 1.d-4)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Externals:
c     Externals:
c     cgratr      - integration routine
c     dr1,dr2,dr3 - functions to integrate
c     HFExc       - Calculates Hartree Fock exchange      
      COMPLEX*16 cgratr, r1, r2, r3, HFExc, Intgrl
      EXTERNAL cgratr, r1, r2, r3, HFExc, Intgrl
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c     Initialization:
      NSing  = 0
      NCalls = 0
      MaxNR  = 0
c     DPPar is array of dp parameters to evaluate functions in cgratr.
c     Everything is in dimensionless units.
c     1. xwg
      DPPar(1) = (Wi/EFermi)
c     2. xgam
      DPPar(2) = gamma/EFermi
c     3. xe
      DPPar(3) = DBLE(Energy)/EFermi
c     4. xeg (gap energy)
      DPPar(4) = 0.d0
c     CPar is array of complex parameters to evaluate functions in cgratr.
c     ck in dimensionless units.
      CPar(1) = ck/kFermi
c     2. ce (complex energy)
      CPar(2) = Energy/EFermi
      
c     Josh - This is a possible fix for functions that overlap zero by a large
c     amount so that Wp does not equal Wi. 
*      Beta = 1.d0*Gamma*SQRT(2.d0/(Wi**2+Gamma**2))
c      Beta  = 0.9d0*Gamma*Gamma/(Wi**2+Gamma**2)
c      Wp =  2.d0*Gamma*LOG(Gamma/Beta) + 2.d0*ATAN2(Wi,Gamma) -
c     &     Gamma*LOG(Gamma**2 + Wi**2)
c      Wp = SQRT(Wi*Wp)
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      DO i1 = 1, 1
         IF(i1.eq.2) THEN
            DPPar(1) = 0.d0
            DPPar(2) = 0.9d0*Gamma/EFermi
         END IF         
c     Calculate integrals in eq. 13 of H.L.
c     1)
c     Integral { ln[((ck+q)**2-E+Wq)/((ck-q)**2-E+Wq)] }
c     from ck + kFermi to Inf.
         Limit1 = ck/kFermi+1.d0
         Limit2 = Inf
      
c     Find singularities in r2
         iFcn = 2
         CALL FndSng(Limit1, Limit2, NSing, XSing, DPPar, CPar, iFcn)
      
c     Calculate integral
         HLInt1 = cgratr(r2,DPPar,CPar,Limit1,Limit2,AbsErr,RelErr,
     &        NSing, XSing,Error,NCalls,MaxNR)
         DO i2 = 1, NSing
            XSing(i2) = (1.d0,0)*XSing(i2)
         END DO

c     2)
c     Integral { ln[(kFermi**2-E-Wq)/(kFermi**2-E+Wq)*
c                     ((ck+q)**2-E+Wq)/((ck-q)**2-E-Wq)] }
c     From ck - kFermi to ck + kFermi
         Limit1 = MAX(ABS(DBLE(ck)/kFermi-1.d0), ZeroPl)
         Limit2 = ck/kFermi+1.d0

c     Find singularities in r1
         IFcn = 1
         CALL FndSng(Limit1, Limit2, NSing, XSing, DPPar, CPar, iFcn)

c     Calculate integral
         HLInt2 = cgratr(r1,DPPar,CPar,Limit1,Limit2,AbsErr,RelErr,
     &        NSing, XSing,Error,NCalls,MaxNR)
         DO i2 = 1, NSing
            XSing(i2) = (1.d0,0)*XSing(i2)
         END DO      
            
c     3)
c     Theta(kFermi-Re(ck)) *
c     Integral { ln[((ck+q)**2-E-Wq)/((ck-q)**2-E-Wq)] } +
c     Theta(Re(ck)-kFermi) *
c     Integral { ln[((ck+q)**2-E+Wq)/((ck-q)**2-E+Wq)] }
c     (Integrals from 0 to kFermi - k and 0 to k - kFermi)
         Limit1 = ZeroPl
         Limit2 = ABS(DBLE(ck)/kfermi-1.d0)

c     If ck = kFermi, HLInt3 = 0
         IF((ABS(DBLE(ck)-kFermi).lt.ZeroPl).or.
     &        (DBLE(Limit2).le.DBLE(Limit1))) THEN
            HLInt3 = 0.d0
c     If ck < kFermi, HLint3 = Integral { ln[((ck+q)**2-E-Wq)/((ck-q)**2-E-Wq)] }
         ELSEIF(DBLE(ck).lt.kFermi) THEN
            Limit2 = 1.d0 - DBLE(ck)/kFermi
            
c     Find singularities in r3
            iFcn = 3         
            CALL FndSng(Limit1, Limit2, NSing, XSing, DPPar, CPar, iFcn)
c     Calculate integral
            HLInt3 = cgratr(r3,DPPar,CPar,Limit1,Limit2,AbsErr,RelErr,
     &           NSing,XSing,Error,NCalls,MaxNR)
            DO i2 = 1, NSing
               XSing(i2) = (1.d0,0)*XSing(i2)
            END DO
         
c     Else ck > kFermi, HLint3 = Integral { ln[((ck+q)**2-E+Wq)/((ck-q)**2-E+Wq)] }
         ELSE
            Limit2 = DBLE(ck)/kFermi - 1.d0
            
c     Find singularities in r2
            iFcn = 2
            CALL FndSng(Limit1, Limit2, NSing, XSing, DPPar, CPar, iFcn)
            
c     Calculate integral
            HLInt3 = cgratr(r2,DPPar,CPar,Limit1,Limit2,AbsErr,RelErr,
     &           NSing,XSing,Error,NCalls,MaxNR)
            DO i2 = 1, NSing
               XSing(i2) = (1.d0,0)*XSing(i2)
            END DO
         END IF
                           
         IF(i1.eq.1) THEN
            Sigma1 = - Amp*Wi*(Wi-coni*Gamma)/(2.d0*pi*EFermi*ck)*
     &           (HLInt1 + HLInt2 + HLInt3)
         ELSE
            Sigma1 = Sigma1 - Amp*Beta*
     &           Wi*(coni*Gamma)/(pi*EFermi*ck)*
     &           (HLInt1 + HLInt2 + HLInt3)
         END IF
      END DO
      
      RETURN
      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      FUNCTION dSigma(ck,Energy,Wi,Gamma,Amp,kFermi,EFermi)
c     Written by Josh Kas
c     Function dSigma calculates dSigma(ck,Energy)/dE from H.L.
c     electron gas model.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
c     Input:
c     ck     - complex momentum
c     Energy - Energy
c     Wi     - Plasmon pole energy
c     Gamma  - Broadening of plasmon pole.
c     Amp    - Amplitude of plasmon pole.
c     kFermi - Fermi momentum
c     EFermi - Fermi energy
c              This is used when calculating dSigma/dE
      DOUBLE PRECISION Wi, Gamma, Amp, kFermi, EFermi
      COMPLEX*16 ck, Energy
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Output:
c     Sigma(ck,Energy)
      COMPLEX*16 dSigma
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Local Variables:
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Local Variables:
c     Wp     - Plasmon frequency, for broadened poles, Wp is not the
c              same as Wi.
c     Beta   - g_i for the negative weight pole at zero. 
c              Adding the negative weight pole at zero corrects the
c              diverging sum rule for epsilon^-1. This is irrelevant
c              for unbroadened poles.
c     NSing  - Number of singularities of integrand (used in cgratr)
c     NCalls - Number of fcn calls used in cgratr
c     MaxNR  - Max number of regions used in cgratr
c     DPPar  - Array of double precision parameters passed to cgratr
c              to be used in the functions r1, r2, and r3
c     CPar   - Array of complex parameters passed to cgratr
c              to be used in the functions r1, r2, and r3
c     Limit1 - Lower limit of integration
c     Limit2 - Upper limit of integration
c     HLInt1 - Integral of dr2 (derivative of first integral in eq. 13 of H.L.)
c     HLInt2 - Integral of dr1 (derivative second integral in eq. 13 of H.L.)
c     HLInt3 - Integral of dr1 or dr3 (derivative of 3rd or 4th integral)
c     XSing  - Array of singularities of integrand (used by cgratr)      
      INTEGER NSing, NCalls, MaxNR
      DOUBLE PRECISION DPPar(10), Wp, Beta
      COMPLEX*16 CPar(10), Limit1, Limit2, HLInt1, HLInt2, HLInt3,
     &     XSing(20)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c     Loop Variables:
      INTEGER i1, i2
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c     Parameters:
c     ZeroPl - lower bound Limit1
c     Inf    - Upper bound of Limit2
c     AbsErr - absolute error used by cgratr
c     RelErr - Relative error used by cgratr
c     Error  - used for error codes by cgratr
      DOUBLE PRECISION ZeroPl, Inf, AbsErr, RelErr, Error
      PARAMETER(ZeroPl = 1.d-5, Inf = 1.d2)
      PARAMETER(AbsErr = 1.d-5, RelErr = 1.d-4)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Externals:
c     cgratr      - integration routine
c     dr1,dr2,dr3 - functions to integrate
c     HFExc       - Calculates Hartree Fock exchange
      COMPLEX*16 cgratr, dr1, dr2, dr3, HFExc
      EXTERNAL cgratr, dr1, dr2, dr3, HFExc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc            
c     Initialization:
      NSing  = 0
      NCalls = 0
      MaxNR  = 0
c     DPPar is array of dp parameters to evaluate functions in cgratr.
c     Everything is in dimensionless units.
c     1. xwg
      DPPar(1) = Wi/EFermi
c     2. xgam
      DPPar(2) = Gamma/EFermi
c     3. xe
      DPPar(3) = Energy/EFermi
c     4. xeg
      DPPar(4) = 0.d0
c     CPar is array of complex parameters to evaluate functions in cgratr.
c     ck in dimensionless units.
      CPar(1) = ck/kFermi
c     2. ce (complex energy)
      CPar(2) = Energy/EFermi + coni*DPPar(2)
      
c     Josh - This is a possible fix for functions that overlap zero by a large
c     amount so that Wp does not equal Wi. 
c     Wp= pi*Wi/2 + Wi*ArcTan[Wi/Gamma] - Gamma*Log[Beta] + Gamma*Log[Gamma] - 
c               1/2*Gamma*Log[Wi**2 + Gamma**2]
c      Beta = 1.d0*Gamma*SQRT(2.d0/(Wi**2+Gamma**2))
c      Beta  = 0.9d0*Gamma*Gamma/(Wi**2+Gamma**2)
c      Wp =  2*Gamma*LOG(Gamma/Beta) + 2*ATAN2(Wi,Gamma) -
c     &     Gamma*LOG(Gamma**2 + Wi**2)
c      Wp = SQRT(Wi*Wp)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      DO i1 = 1, 1
         IF(i1.eq.2) THEN
            DPPar(1) = 0.d0
            DPPar(2) = 0.9d0*Gamma/EFermi
         END IF         
         
c     Calculate derivatives of integrals in eq. 13 of H.L.      
c     1)
c     Integral d/dE{ ln[((ck+q)**2-E+Wq)/((ck-q)**2-E+Wq)] }
c     from ck + kFermi to Inf.      
         Limit1 = ck/kFermi+1.d0
         Limit2 = Inf
c     Find singularities in dr2
         iFcn = 2
         CALL FndSng(Limit1, Limit2, NSing, XSing, DPPar, CPar, iFcn)
         
c     Calculate integral
         HLInt1 = cgratr(dr2,DPPar,CPar,Limit1,Limit2,AbsErr,RelErr,
     &        NSing, XSing,Error,NCalls,MaxNR)
         DO i2 = 1, NSing
            XSing(i2) = (1.d0,0)*XSing(i2)
         END DO
         
c     2)
c     Integral d/dE{ ln[(kFermi**2-E-Wq)/(kFermi**2-E+Wq)*
c     ((ck+q)**2-E+Wq)/((ck-q)**2-E-Wq)] }
c     From ck - kFermi to ck + kFermi
         Limit1 = MAX(ABS(DBLE(ck)/kFermi-1.d0), ZeroPl)
         Limit2 = ck/kFermi+1.d0
         
c     Find singularities in dr1
         iFcn = 1
         CALL FndSng(Limit1, Limit2, NSing, XSing, DPPar, CPar, iFcn)
         HLInt2 = cgratr(dr1,DPPar,CPar,Limit1,Limit2,AbsErr,RelErr,
     &        NSing, XSing,Error,NCalls,MaxNR)
         DO i2 = 1, NSing
            XSing(i2) = (1.d0,0)*XSing(i2)
         END DO
         
c     3)
c     Theta(kFermi-Re(ck)) *
c     Integral { ln[((ck+q)**2-E-Wq)/((ck-q)**2-E-Wq)] +
c     Theta(Re(ck)-kFermi) *
c     Integral { ln[((ck+q)**2-E+Wq)/((ck-q)**2-E+Wq)]
c     (Integrals from 0 to kFermi - k and 0 to k - kFermi)
         Limit1 = ZeroPl
         Limit2 = ABS(DBLE(ck)/kfermi-1.d0)
         
c     If ck = kFermi, HLInt3 = 0      
         IF((ABS(DBLE(ck)-kFermi).lt.ZeroPl).or.
     &        (DBLE(Limit2).le.DBLE(Limit1))) THEN         
            HLInt3 = 0.d0
c     If ck < kFermi, HLint3 = Integral d/dE{ ln[((ck+q)**2-E-Wq)/((ck-q)**2-E-Wq)] }          
         ELSEIF(DBLE(ck).lt.kFermi) THEN
            Limit2 = 1.d0 - DBLE(ck)/kFermi
            
c     Find singularities in r3        
            iFcn = 3
            CALL FndSng(Limit1, Limit2, NSing, XSing, DPPar, CPar, iFcn)
            
c     Calculate integral         
            HLInt3 = cgratr(dr3,DPPar,CPar,Limit1,Limit2,AbsErr,RelErr,
     &           NSing,XSing,Error,NCalls,MaxNR)
            DO i2 = 1, NSing
               XSing(i2) = (1.d0,0)*XSing(i2)
            END DO
c     Else ck > kFermi, HLint3 = Integral d/dE{ ln[((ck+q)**2-E+Wq)/((ck-q)**2-E+Wq)] }         
         ELSE
            Limit2 = DBLE(ck)/kFermi - 1.d0
            
c     Find singularities in r2         
            iFcn = 2
            CALL FndSng(Limit1, Limit2, NSing, XSing, DPPar, CPar, iFcn)
            
c     Calculate integral         
            HLInt3 = cgratr(dr2,DPPar,CPar,Limit1,Limit2,AbsErr,RelErr,
     &           NSing,XSing,Error,NCalls,MaxNR)
            DO i2 = 1, NSing
               XSing(i2) = (1.d0,0)*XSing(i2)
            END DO
         END IF

         IF(i1.eq.1) THEN
            dSigma = - Amp*Wi*(Wi-coni*Gamma)/(2.d0*pi*EFermi*ck)*
     &           (HLInt1 + HLInt2 + HLInt3)
         ELSE
            dSigma = dSigma - Amp*Beta*
     &           Wi*(coni*Gamma)/(pi*EFermi*ck)*
     &           (HLInt1 + HLInt2 + HLInt3)
         END IF
      END DO
      
      RETURN
      END      

      FUNCTION HFExc(ckIn, EFermi, kFermi)
c     returns dirac-hara hartree-fock exchange
c     ck - complex momentum in units of kFermi
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
      COMPLEX*16 ckIn, ck, HFExc, c
      DOUBLE PRECISION EFermi, kFermi
      ck = ckIn/kFermi
      c=-2.d0*EFermi/(pi*kFermi)
      IF(ABS(ck-1.d0).le.0.00001d0) THEN
         HFExc = c
      ELSE
         HFExc = c*(1.d0+(1.d0/ck-ck)* log( (1.d0+ck)/(ck-1.d0) )/2.d0)
      END IF
      RETURN
      END
      
c****************************************************************************
c     the following function routines are used for evaluating integals and
c     their derivatives.
c****************************************************************************
      complex*16 function r1(q,dppar,cpar)
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
      
c     Input:
      double precision dppar(4)
c     dppar contains:
c     xe, xeg, xwg, xgam
      complex*16 cpar(2)

c     Local Variables:
      complex*16 fq,fqq,fiq,a1,a2,a3,a4,t1,t2,q,ck,xe
      external fq

      ck=CPar(1)
c     3. xe
      xe = CPar(2)
c     4. xeg
      xeg = DPPar(4)
c     print*, 'call fq(q),ck', q, ck
      fqq=fq(q,dppar)
c     print*,'fqq=', fqq
      fiq=1./(q*fqq)
      a1=1.d0-xeg-xe-fqq - coni*1.d-10
      a2=(ck+q)**2-xe+fqq - coni*1.d-10
      a3=(ck-q)**2-xe-fqq - coni*1.d-10
      a4=1.d0+xeg-xe+fqq - coni*1.d-10
      t1=(a1*a2)
      t2=(a3*a4)
c     print*,'a1,a2,a3,a4,t1,t2',a1,a2,a3,a4,t1,t2
c      t1=t1/t2
      r1=fiq*(log(a1)+log(a2)-log(a3)-log(a4))
c     Test with r=E
c      r1=xe
c     print*,'r1 return to cgratr', r1
      return
      end
c****************************************************************************
      complex*16 function dr1(q,dppar,cpar)
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
      
c     Input:
      double precision dppar(4)
c     dppar contains:
c     xe, xeg, xwg, xgam
      complex*16 cpar(2)

c     Local Variables:
      complex*16 fq,fqq,fiq,a1,a2,a3,a4,t1,t2,q,ck,xe
      external fq

      ck=CPar(1)
c     3. xe
      xe = CPar(2)
c     4. xeg
      xeg = DPPar(4)
c     print*, 'call fq(q),ck', q, ck
      fqq=fq(q,dppar)
c     print*,'fqq=', fqq
      fiq=1./(q*fqq)
      a1=1.d0-xeg-xe-fqq - coni*1.d-10
      a2=(ck+q)**2-xe+fqq - coni*1.d-10
      a3=(ck-q)**2-xe-fqq - coni*1.d-10
      a4=1.d0+xeg-xe+fqq - coni*1.d-10
c     print*,'a1,a2,a3,a4,t1,t2',a1,a2,a3,a4,t1,t2
c      t1=t1/t2
      dr1 = -fiq*(1.d0/a1+1.d0/a2-1.d0/a3-1.d0/a4)
c     Test with r=E
c      dr1=1.d0
c      write(51,*) dble(q), dble(dr1)
c     print*,'r1 return to cgratr', r1
      return
      end
c**********************************************************************
      complex*16 function r2(q,dppar,cpar)
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
      
c     Input:
      double precision dppar(4)
c     dppar contains:
c     xe, xeg, xwg, xgam
      complex*16 cpar(2)

c     Local Variables:
      complex*16 fqq,fq,fiq,a1,a2,q,ck,xe
      ck = CPar(1)
c     3. xe
      xe = CPar(2)
c     4. xeg
c      xeg = DPPar(4)
      
      fqq=fq(q,dppar)
      fiq=1.d0/(q*fqq)
      a1=((ck+q)**2-xe+fqq) - coni*1.d-10
      a2=((ck-q)**2-xe+fqq) - coni*1.d-10
      r2=fiq*(log(a1)-log(a2))
c     Test with r=E
c      r2=xe      
30    return
      end
c**********************************************************************
      complex*16 function dr2(q,dppar,cpar)
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
      
c     Input:
      double precision dppar(4)
c     dppar contains:
c     xe, xeg, xwg, xgam
      complex*16 cpar(2)

c     Local Variables:
      complex*16 fqq,fq,fiq,a1,a2,q,ck,xe
      ck = CPar(1)
c     3. xe
      xe = CPar(2)
c     4. xeg
c      xeg = DPPar(4)
      
      fqq=fq(q,dppar)
      fiq=1.d0/(q*fqq)
      a1=((ck+q)**2-xe+fqq) - coni*1.d-10
      a2=((ck-q)**2-xe+fqq) - coni*1.d-10
      dr2=-fiq*(1.d0/a1-1.d0/a2)
c     Test with r=E
c      dr2=1.d0      
c      write(52,*) dble(q), dble(dr2)
30    return
      end  
c**********************************************************************
      complex*16 function r3(q,dppar,cpar)
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

c     Input:
      double precision dppar(4)
c     dppar contains:
c     xe, xeg, xwg, xgam
      complex*16 cpar(2)

c     Local Variables:
      complex*16 fqq,fq,fiq,a1,a2,q,ck,xe
      ck = CPar(1)
c     3. xe
      xe = CPar(2)
c     4. xeg
c      xeg = DPPar(4)
c     valid only for k<kf, q<kf-k ?
      fqq=fq(q,dppar)
      fiq=1.d0/(q*fqq)
      a1=( (ck+q)**2-xe-fqq) - coni*1.d-10
      a2=( (ck-q)**2-xe-fqq) - coni*1.d-10
      r3=fiq*(log(a1) - log(a2))
c     Test with r=E
c      r3=xe      
30    return
      end
c**********************************************************************
      complex*16 function dr3(q,dppar,cpar)
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

c     Input:
      double precision dppar(4)
c     dppar contains:
c     xe, xeg, xwg, xgam
      complex*16 cpar(2)

c     Local Variables:
      complex*16 fqq,fq,fiq,a1,a2,q,ck,xe
      ck = CPar(1)
c     3. xe
      xe = CPar(2)
c     4. xeg
c      xeg = DPPar(4)
c     valid only for k<kf, q<kf-k ?
      fqq=fq(q,dppar)
      fiq=1.d0/(q*fqq)
      a1=( (ck+q)**2-xe-fqq) - coni*1.d-10
      a2=( (ck-q)**2-xe-fqq) - coni*1.d-10
      dr3=-fiq*(1.d0/a1-1.d0/a2)
c     Test with r=E
c      dr3=1.d0
c      write(53,*) dble(q), dble(dr2)
30    return
      end     
c**********************************************************************
      complex*16 function fq(q,dppar)
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
      complex*16 q
      double precision dppar(4)

c     1. xwg
      xwg = DPPar(1)
c     2. xgam
      xgam = DPPar(2)
      
c     Here I am going to change the dispersion relation to 
c     wq = wp + 1/2 * q**2 
c     This makes calculation of broadened poles easier. 
c      fq = xwg-coni*xgam + q**2
      
c     fq(q)=w1(q)=sqrt(w1**2+((omega(q)-omega_p)/omega_f)**2)
c     omega(q)**2=omega_p**2+omega_g**2(q)
c     fq(q)=xwg+a2*q**2+a4*q**4    xwg=(w1/ef)**2
c     electron gas parameters xwg=wp**2 a2=4/3, a4=1
c     uncomment the following 4 lines to use the old dispersion relation.
      a2=4.d0/3.d0
      a4=1.d0
      fq=(xwg-coni*xgam)**2 + a2*q*q + a4*q**4
      fq=sqrt(fq)
      return
      end

      FUNCTION Intgrl(func, a, b, NPts, DpPar, CPar)
c     Function Integ integrates func by trapezoidal rule from a to b
c     using NPts steps size.
      COMPLEX*16 Intgrl, a, b, CPar(10)
      DOUBLE PRECISION DpPar(10)
      INTEGER NPts
      COMPLEX*16 func
      EXTERNAL func

      INTEGER i1
      COMPLEX*16 dx, sum, y1, y2, x1, x2


      sum = 0.d0
      dx = (b-a)/DBLE(NPts-1)
      
      x1 = a
      y1 = func(x1, DpPar, CPar)
      
      DO i1 = 1, NPts
         x2 = a + i*dx
         y2 = func(x2, DpPar, CPar)
         sum = sum + (y1+y2)
         y1 = y2
         x1 = x2
      END DO

      Intgrl = sum*dx/2.d0

      RETURN
      END
               
c*********************************************************************
c   This is Steve White's rewrite of Mike Teter's integration routine.  
c   Modified by J. Rehr for complex integration.
c   The following is a listing of the arguments in the initial function 
c   statement:
c     fn    -- routine requires external function statement in MAIN
c     xmin  -- lower limit
c     xmax  -- upper limit
c     abr   -- absolute tolerable error
c     rlr   -- relative tolerable error
c     nsing -- number of singularities or regions requiring 
c     special attention
c     xsing -- array of locations of singularities or endpoints
c     of special regions
c     error -- output for routine error messages
c     numcal-- the number of times fn was called
c     maxns -- the maximum number of regions being considered simultaneously
c     function cgratr(fn,xmin,xmax,abr,rlr,nsing,xsing,error,numcal,maxns)
c     fn declared double precision
c     double precision function cgratr(fn,xmin,xmax,abr,rlr,
c     fn declared complex*16
      
      complex*16 function cgratr(fn,dppar,cpar,xmin,xmax,abr,rlr,
     &     nsing,xsing,error,numcal,maxns)
      implicit double precision (a-h,o-z)
      parameter (mx=1500)
      integer nsing
      complex*16 fn,value,valu,fval(3,mx),xmax,xmin,del,del1
      complex*16 xleft(mx), xsing(20), cpar(10)
      double precision dppar(10)
      external fn
c     dimension xleft(mx),fval(3,mx),dx(3),wt(3)
      dimension wt9(9),dx(3),wt(3)
c     dimension xsing(20)
      logical atsing
      save dx,wt,wt9
      data dx/0.1127016653792583  ,0.5  ,0.8872983346207417  /
      data wt/0.277777777777777778  ,0.4444444444444444444  ,
     1     0.2777777777777777778  /
      data wt9/0.0616938806304841571  ,0.108384229110206161  ,
     1     0.0398463603260281088  ,0.175209035316976464  ,
     2     0.229732989232610220  ,0.175209035316976464  ,
     3     0.0398463603260281088  ,0.108384229110206161  ,
     4     0.0616938806304841571  /
c     nstack is the number of different intervals into which the 
c     integration region is currently divided. The number of regions can
c     grow if more accuracy is needed by dividing the right-most region
c     into three regions. The number of regions shrinks when the integral
c     over the right-most region is accurate enough, in which case that
c     integral is added to the total (stored in cgratr) and the region
c     is removed from consideration (and a new region is the right-most).
      nstack=nsing+1
      maxns = nstack
      error=0.  
      cgratr=0.  
c     The array xleft stores the boundary points of the regions.
c     The singular points just govern the initial placement of the regions.
      xleft(1)=xmin
      xleft(nsing+2)=xmax
      if(nsing.gt.0) then
         do 9 j=1,nsing
            xleft(j+1)=xsing(j)
 9       continue
      endif
c     For each region, calculate the function and store at three selected points.
      do 1 jj=1,nstack
         del=xleft(jj+1)-xleft(jj)
c     print*, 'fn call j= ,'
         do 1 j=1,3
c     print*, 'fn call in cgratr j= ',j
            fval(j,jj)=fn(xleft(jj)+del*dx(j),dppar,cpar)
 1    continue
c     print*, 'output of fn call, fval(j,jj)',fval(j,jj)
      numcal = nstack * 3
 6    continue
      if(nstack+3.ge.mx) then
         write(*,*) ' TOO MANY REGIONS'
         stop 0006
      endif
c     Divide the rightmost region into three subregions.  
      del=xleft(nstack+1)-xleft(nstack)
      xleft(nstack+3)=xleft(nstack+1)
      xleft(nstack+1)=xleft(nstack)+del*dx(1)*2.
      xleft(nstack+2)=xleft(nstack+3)-del*dx(1)*2.
c     The three data points already found for the region become the 
c     middle data points (number 2 in first index of fval) for each region.
      fval(2,nstack+2)=fval(3,nstack)
      fval(2,nstack+1)=fval(2,nstack)
      fval(2,nstack)=fval(1,nstack)
c     Now do the integral over the right-most region in two different ways-
c     a three point integral (valu) over each of the three subregions 
c     and a more accurate nine-point integral (value) over whole region.
c     valu is used only for the error estimate.
      icount=0
      value=0.  
      valu=0.  
      do 3 j=nstack,nstack+2
         del1=xleft(j+1)-xleft(j)
c     print*, 'fn call 2'
         fval(1,j)=fn(xleft(j)+dx(1)*del1,dppar,cpar)
         fval(3,j)=fn(xleft(j)+dx(3)*del1,dppar,cpar)
c     print*, 'fn call 2'
         numcal = numcal + 2
         do 5 k=1,3
            icount=icount+1
            value=value+wt9(icount)*fval(k,j)*del
            valu=valu+fval(k,j)*wt(k)*del1
 5       continue
 3    continue
      dif=abs(value-valu)
c     If the following condition is true, add in this integral to the total,
c     and reduce the number of regions under consideration.
      frac = del / (xmax - xmin)
      atsing = .false.
      if(frac .le. 1.0e-8) atsing = .true.
      if(dif .le. abr*frac .or. dif.le.rlr*abs(value) .or. 
     1     (atsing .and. 
     2     (frac .le. 1.0e-15 .or. dif .le. abr*0.1  ))) then
c     The following commented out line is Teeter's old error criterion.
c     if(dif.le.abr.or.dif.le.rlr*abs(value))then
         cgratr=cgratr+value
         error=error+abs(dif)
         nstack=nstack-1
c     If no more regions, we are done.
         if(nstack.le.0) return
      else
c     If the integration is insufficiently accurate, make each of the 
c     three subregions of the right-most region into regions.
c     On next pass the right-most of these is the new current region.
         nstack=nstack+2
         maxns = max(maxns,nstack)
      endif
      go to 6
      end
      SUBROUTINE CSigZ(Energy, Mu, Rs, ReSig, ImSig, ZTot, WpScl, Gamma,
     &     AmpFac)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Written by Josh Kas
c     This subroutine calculates the self energy Sigma(k(Energy),Energy)
c     based on an electron gas model of epsilon^-1.
c      
c     Solve: k0**2 = 2*Energy - 2*Mu -
c                    2*(Sigma(k0,Energy)-Sigma(kFermi,EFermi))
c            
c     Steps:
c
c            1. k0**2  = 2*(Energy-Mu) + SigmaF (SigmaF is self energy at fermi level).
c            2. Sigma0 = Sigma(k0,Energy)
c            3. Find derivative w.r.t. E dSgdE
c            4. k1**2  = 
c                  k0**2 - 2*(Sigma0-SigmaF)/(1-dSgdE)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
c     Parameters:
c     MxPole - Maximum number of poles
      INTEGER MxPole
      PARAMETER(MxPole=1000)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Input:
c     Energy - Energy at which to evaluate Sigma
c     Mu     - Fermi energy as calculated by FEFF.
c     Rs     - R sub s (sphere of radius Rs holds charge e)
c     WpScl  - Scale Wp in interstitial by WpScl
c     Gamma  - Use broadening Gamma when calculating Sigma
c     AmpFac - Use amplitude AmpFac for plasmon pole.
      DOUBLE PRECISION Rs, WpScl(MxPole), Gamma(MxPole),
     &     AmpFac(MxPole), Mu
      COMPLEX*16 Energy
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Output:
c     ReSig  - Re[Sigma(k,e)]
c     ImSig  - Im[Sigma(k,e)]
c     ZTot   - Renormalization factor Z = 1/(1-dSgdE)
c     Note: Atomic units are used.
      DOUBLE PRECISION ReSig, ImSig
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Local Variables:
c     kFermi - Fermi momentum calculated from Rs via electron gas
c              approximation.
c     EFermi - Fermi energy = kFermi^2/2
c     Wp     - Electron gas plasmon frequency.
c     Gam    - broadening for broadened pole.
c     ckF    - complex variable to store kFermi
c     ck0    - complex momentum
c     SigmaF - Self energy at the fermi energy and fermi momentum
c              Does not include the Hartree Fock part.
c     Sigma0 - Single pole self energy evaluated at
c              Wp = Wp(ipole)/Wp(R_Interstitial)*Wp(Rs)
c     dSgdE  - Derivative of sigma w.r.t. Energy
c     ZTot   - Renormalization factor Z = 1/(1-dSgdE)
c     RelEn  - Energy relative to the fermi energy from FEFF      
c              Relen = Energy - Mu + EFermi
c     SigTot - Total many pole deltaSigma = Sigma(E,k(E))-Sigma(EF,kF)
c     DelHF  - Sigma_HartreeFock(k) - Sigma_HartreeFock(kF)
      DOUBLE PRECISION kFermi, EFermi, Wp, Gam
      COMPLEX*16 ckF, ck0, SigmaF, Sigma0, dSgdE, ZTot, RelEn, RelEnP,
     &     SigTot, DelHF, SigmaP
c     Loop variables
      INTEGER i1, i2

c     Parameters:
      DOUBLE PRECISION DPZero
      PARAMETER(DPZero = 0.d0)
      INTEGER MxIter
      PARAMETER(MxIter = 1)

c     Externals:
c     Sigma1 - calculates the energy dependent part of self energy
c              for a single pole.
c     dSigma - calculates derivative of self energy w.r.t energy.
c     HFExt  - calculates Hartree Fock exchange
      COMPLEX*16 Sigma1, dSigma, HFExc
      EXTERNAL Sigma1, dSigma, HFExc
      
c     Initialization
      ZTot = 0.d0
      kFermi = fa/Rs
      EFermi = kFermi*kFermi/2.d0
      SigTot=0.d0
      dSgdE = 0.d0
      SigmaF = 0.d0
      Gam = 0.d0

c     Loop1: Start self consistency loop.
c     This does not seem to work, so MxIter = 1
      DO i2 = 1, MxIter

c        Loop2: Loop over poles to get SigmaF
         DO i1 = 1, MxPole
            IF(WpScl(i1).lt.0) GOTO 5
            
c           Wp is in Hartrees            
            Wp = SQRT(3.d0/Rs**3)*WpScl(i1)
c            Gam = Gamma(i1)
            
c           find Sigma_Fermi (SigmaF)
            ckF = kFermi*1.00001d0
            RelEn = EFermi
            SigmaF = SigmaF + Sigma1(ckF,RelEn,Wp,Gam,AmpFac(i1),
     &           kFermi, EFermi)
         END DO
c        End Loop2
 5       CONTINUE
         
         dsgdE = 0.d0
c        Loop3: Loop over poles         
         DO i1 = 1, MxPole
            IF(WpScl(i1).lt.0) GOTO 10
c           Wp is in Hartrees
            Wp = SQRT(3.d0/Rs**3)*WpScl(i1)
c            Gam = Gamma(i1)
c           Start with ck0=Sqrt[Re(Energy)-Mu+EFermi]
            RelEn = DBLE(Energy) - Mu + EFermi
            ck0 = SQRT(2.d0*DBLE(RelEn))
            
c           Find Sigma0 = Sigma(ck0,E); ck0=SQRT(2*(Energy-Mu))
            Sigma0 = Sigma1(ck0,RelEn,Wp,Gam,AmpFac(i1),kFermi,
     &           EFermi)
                  
            RelEnP = RelEn*0.001d0
            SigmaP = Sigma1(ck0,RelEnP,Wp,Gam,AmpFac(i1),kFermi,
     &           EFermi)

c            write(71,*) DBLE(RelEn-EFermi), DBLE(dSgdE), DIMAG(dSgDE)
            dSgdE = dSgdE + (SigmaP - Sigma0)/(RelEnP-RelEn)
c            dSgdE = dSgdE + 
c     &          dSigma(ck0,RelEn,Wp,Gam,AmpFac(i1),kFermi,EFermi)
c           Uncomment the following line to print derivative            
c            write(72,*) DBLE(RelEn-EFermi), DBLE(dSgdE2), DIMAG(dSgDE2)         

c           SigTot is sum of poles
            SigTot = SigTot + Sigma0
            
c        End Loop3: loop over poles.            
         END DO         
 10      CONTINUE
         
c     End Loop1: self-consistency loop      
      END DO


c     Add Hartree Fock part of delta Sigma.
      DelHF = HFExc(ck0,EFermi,kFermi) - HFExc(ckF,EFermi,kFermi)
      SigTot = SigTot - SigmaF + DelHF

c     Form ZTot and return Re and Im parts of Sigma.
      ZTot = 1.d0/(1.d0-dSgdE)
c      ZTot = 1.d0
      SigTot = ZTot*(SigTot)
      
      ReSig = DBLE(SigTot)
      ImSig = DIMAG(SigTot)
      
      RETURN
      END

      SUBROUTINE FndSng(Limit1, Limit2, NSing, XSing, DPPar, CPar, iFcn)
c     Josh Kas
c     This subroutine finds the singularities in the integrands of eq. 13
c     in
c     Single-particle Spectrum of the Degenerate Electron Gas
c     II. Numerical Results for Electrons Coupled to Plasmons
c     Phys. kondens. Materie, Bd. 6 (1967)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c     In practice this routine only solves for the singularities of one
c     of the three integrands, then checks to see that the singularity
c     is within the limits of integration, and throws out singularities
c     that are not.
c     In units of the Fermi energy the equations to solve are:
c     1) +/- k*q**3 + 2*(3*k**2 - E - 2/3)*q**2 +/- 4*k*(k**2 - E)*q +
c                     [(k**2 - E)**2 - Wp**2] = 0 
c     2) q**4 + 4/3*q**2 + Wp**2 - (1 - E)**2 = 0
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c     Input:
c     Limit1 - Lower limit of integration
c     Limit2 - Upper limit of integration
c     CPar   - Array of complex parameters passed to function
c              CPar(1) = ck/kFermi
c              CPar(2) = Energy/EFermi + i*Gamma/EFermi
c     DPPar  - Array of double precision parameters passed to function
c              DPPar(1) = Wp/EFermi
c              DPPar(2) = Gamma/EFermi
c              DPPar(3) = Energy/EFermi
c              DPPar(4) = xeg (gap energy)
c     iFcn   - Integer denoting which function is the integrand
c              iFcn = 1: solve eqs 1 and 2 for q
c              iFcn = 2: solve eq 1 for q      
      COMPLEX*16 Limit1, Limit2, CPar(10)
      DOUBLE PRECISION DPPar(10)
      INTEGER iFcn
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Output:
c     XSing  - Array of singularities
c     NSing  - Number of singularities
      COMPLEX*16 XSing(20)
      INTEGER NSing
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Local Variables:
c     Coef   - Coefficients of q**n of eq. to solve.
c              eq = Coef(1)*q**n + Coef(2)*q**(n-1)...
c     Sol(4) - Array of solutions to the equation.
c     XSing2 - Temp XSing      
c     Test   - Used to test solution of equation.
c     Zero   - Tolerance for testing solution to eqs.
c     NSol   - Number of solutions to eq.
c     Order  - Used to order singulaties from smallest to largest
      COMPLEX*16 Coef(4), Sol(4), XSing2(4)
      DOUBLE PRECISION Test, Zero
      INTEGER NSol, Order(100)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Loop variables
      INTEGER i1, i2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Initialization
      NSing = 0
      Zero=1.d-4

c     Solve eq 1 for q with + sign
      Coef(1) = 4.d0*CPar(1)
      Coef(2) = 2.d0*(3.d0*CPar(1)**2 - DPPar(3) - 2.d0/3.d0)
      Coef(3) = 4.d0*CPar(1)*(CPar(1)**2 - DPPar(3))
      Coef(4) = (CPar(1)**2 - DPPar(3))**2 - DPPar(1)**2
      
      CALL CCubic(Coef, Sol, NSol)
      
c     Test solutions. If Sol is a solution and it is real
c     and it lies between Limit1 and Limit2, add it to list
c     of singularities.         
      DO i1 = 1, NSol
         Test = ABS((CPar(1)+Sol(i1))**2 - DPPar(3) +
     &        SQRT(Sol(i1)**4 + 4.d0/3.d0*Sol(i1)**2 + DPPar(1)**2))
         IF(Test.lt.Zero) THEN               
            IF((DBLE(Sol(i1)).ge.DBLE(Limit1)).and.
     &           (DBLE(Sol(i1)).le.DBLE(Limit2)).and.
     &           (ABS(DIMAG(Sol(i1))).le.Zero)) THEN               
               NSing = NSing + 1
               XSing(NSing) = DBLE(Sol(i1))
            END IF
         END IF
      END DO
      
c     Now solve eq. 1 for q with - sign
      Coef(1) = -Coef(1)
      Coef(3) = -Coef(3)
      
      CALL CCubic(Coef, Sol, NSol)
      
c     Test solutions as before.
      DO i1 = 1, NSol
         Test = ABS((CPar(1)-Sol(i1))**2 - DPPar(3) -
     &        SQRT(Sol(i1)**4 + 4.d0/3.d0*Sol(i1)**2 + DPPar(1)**2))
         IF(Test.lt.Zero) THEN
            IF((DBLE(Sol(i1)).ge.DBLE(Limit1)).and.
     &           (DBLE(Sol(i1)).le.DBLE(Limit2)).and.
     &           (ABS(DIMAG(Sol(i1))).le.Zero)) THEN
               NSing = NSing + 1
               XSing(NSing) = DBLE(Sol(i1))
            END IF
         END IF
      END DO

c     If iFcn = 1 (Solving for singularities of r1(q))
      IF(iFcn.eq.1) THEN         
c        Solve eq. 2 for q
         Coef(1) = 1.d0
         Coef(2) = 4.d0/3.d0
         Coef(3) = DPPar(1)**2

         CALL CQdrtc(Coef,Sol,NSol)
         DO i1 = 1, NSol
            XSing2(2*i1-1) =  DBLE(SQRT(Sol(i1)))
            XSing2(2*i1)   = -DBLE(SQRT(Sol(i1)))
         END DO

c        Test Solutions
         DO i1 = 1, 2*NSol
            IF((DBLE(XSing2(i1)).ge.DBLE(Limit1)).and.
     &           (DBLE(XSing2(i1)).le.DBLE(Limit2)).and.
     &              (ABS(DIMAG(Sol(i1))).le.Zero)) THEN
               NSing = NSing + 1
               XSing(NSing) = XSing2(i1)
            END IF
         END DO
      END IF
      
c     Sort singularities
      CALL QSORTI(Order,NSing,Xsing)
      DO i1 = 1, NSing
         XSing2(i1) = XSing(i1)
      END DO
      DO i1 = 1, NSing
         XSing(i1) = XSing2(Order(i1))
      END DO

      RETURN
      END
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
      complex*16  function aprdec(ala,bla,lla)
c     the result of this function is the coefficient for the term of
c     power (l-1) for the product of two polynomes, whose coefficients
c     are in rows a and b
 
      implicit double precision (a-h, o-z)
      complex*16 ala (10)
      integer lla
      dimension bla(10)
 
      aprdec = (0.0d0, 0.0d0)
      do 11 m = 1, lla
 11      aprdec = aprdec + ala(m) * bla(lla+1-m)
      return
      end
      double precision function aprdep (a,b,l)
c     need to be in library for ATOM and PHASE; renamed aprdev
c     the result of this function is the coefficient for the term of 
c     power (l-1) for the product of two polynomes, whose coefficients
c     are in rows a and b 
 
      implicit double precision (a-h,o-z)
      dimension a(10),b(10)
 
      aprdep=0.0d 00
      do 11 m=1,l
 11      aprdep=aprdep+a(m)*b(l+1-m)
      return
      end
      subroutine dfovrg (ncycle, ikap, rmt, jlast, jri, p2, dx,
     1                  ri, vxc, vxcval, dgcn, dpcn, adgc, adpc,
     2                  xnval, pu, qu, ps, qs,
     2                  iz, ihole, xion, iunf, irr, ic3)
c     Dirac equation solver for complex energy
c     coded by a.ankudinov 1996
c     modified by a.ankudinov 1997 to get irregular solution 

c     fully relativistic version of subroutine fovrg.f
c     input:
c        ncycle  times to calculate photoelectron wave function
c                with nonlocal exchange
c        ikap    quantum number kappa for photoelectron
c        rmt     muffin-tin radius
c        jri     first interstitial grid point (imt + 1)
c        jlast   last point for integration of Dirac eq.
c        p2      current complex energy
c        dx      dx in loucks' grid (usually .05)
c        ri(nr)  loucks' position grid, r = exp ((i-1)*dx - 8.8)
c        vxc(nr) coulomb+xc potential for total density
c        vxcval  coulomb+xc potential for valence density
c        both vxc and vxcval include coulomb and nuclear potential
c        dgcn(dpcn) large(small) dirac components for 'iph' atom
c        adgc(adpc) their development coefficients
c     work space:
c        must be dimensioned in calling program.  coded like this
c        to make using different r-grids with different nrmax easy.
c
c     output:
c        pu, qu  upper and lower components at muffin tin
c        ps and qs are  upper and lower components for photoelectron

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

      complex*16 vxc(nrptx), vxcval(nrptx), p2
      dimension ri(nrptx)
      complex*16 ph0, amp, pu, qu, vu, vm(nrptx)
      complex*16 ps(nrptx), qs(nrptx), aps(10),aqs(10)

c     all atoms' dirac components and their development coefficients
      dimension dgcn(nrptx,30), dpcn(nrptx,30)
      dimension adgc(10,30), adpc(10,30)
 
c     iph atom's dirac components and their development coefficients
      common/dff/cg(nrptx,30), cp(nrptx,30), bg(10,30), bp(10,30),
     1             fl(30), fix(30), ibgp
c     fl power of the first term of development limits.
c     ibgp first dimension of the arrays bg and bp (=10)

      complex*16 gg,gp,ag,ap,dv,av,bid
      common/comdic/cl,dz,gg(nrptx),ag(10),gp(nrptx),ap(10),
     1              dv(nrptx),av(10),bid(2*nrptx+20)
c      gg,gp are the input and output for solout
      common/itescf/testy,rap(2),teste,nz,norb,norbsc
      common/mulabc/afgkc
      dimension afgkc(-ltot-1:ltot,30,0:3)
      common/messag/dlabpr,numerr
      character*8 dlabpr
c      xnel here - number of core electrons only
      common/ratom1/xnel(30),en(30),scc(30),scw(30),sce(30),
     1nq(30),kap(30),nmax(30)
      common/scrhf1/eps(435),nre(30),ipl
      common/snoyac/dvn(nrptx),anoy(10),nuc
      common/tabtec/hx,dr(nrptx),test1,test2,ndor,np,nes,method,idim
      dimension xnval(30)

c     initialize the data and test parameters
      ndor = 3
      cl = alpinv
      if (irr.gt.0) then
c        for irregular solution
         ndor=2
         aps(1) =  pu
         aqs(1) =  qu
         do 5 i=1, jri
           gg(i) = ps(i)
           gp(i) = qs(i)
 5       continue
      endif
      do 9 i = jri+1,nrptx
         vxc(i)=vxc(jri+1)
 9       vxcval(i)=vxc(jri+1)
      ibgp=10
      numerr = 0
      nz = iz
      hx = dx
      idim= 1 + nint(250*0.05/dx)
      if (idim .gt. nrptx) idim = nrptx
      if (mod(idim,2) .eq. 0) idim=idim-1
      
c     numerical integration of Dirac eq. works if you have 6 grid points
c     for one period of oscillations, switch to analytical expression
c     for a steplike potential  at large distances
      aa = 0.5
c     if (irr.gt.0) aa = 0.05
      rwkb = aa / dx / sqrt(abs(2*p2+(p2/cl)**2))
      x0 = 8.8
      iwkb= (log(rwkb) + x0) / dx  +  2
      if (iwkb.gt.idim) iwkb = idim
      if (iwkb.lt. 10) iwkb = 10
      
c     copy information into common's of atomic code
      do 13 j=1,30
      do 13 i=1,10
         bg(i,j)=adgc(i,j) 
 13      bp(i,j)=adpc(i,j) 
      do 15 j=1,30
      do 15 i=1,idim
         cg(i,j)=dgcn(i,j) 
 15      cp(i,j)=dpcn(i,j) 

      call inmuac (ihole, xion, iunf, ikap)
      nmax(norb)=jlast
      if (iwkb.ge. jlast-1) iwkb = idim
c     note that here norb correspond to photoelectron

c     calculate initial photoelectron orbital using lda
      call diff (vxc,ri,ikap,cl,hx,jri,vm)
      do 18 i = jri, nrptx
  18  vm(i)=0.0d0
      call wfirdc (p2,kap,nmax,vxc,ps,qs,aps,aqs,irr,ic3,vm,
     1             rmt,jri, iwkb)

      if (numerr .ne. 0) call par_stop('error in wfirdc')
      if (ncycle .eq. 0) go to 999

c     to get orthogonalized photo e w.f., use alternative exit below
c     in general it should not be orthogonolized. Use for testing only 
c     ala

c     further need only core electrons for exchange term
      do 40 i=1, norb-1
  40  xnel(i) = xnel(i) - xnval(i)
c     take vxcval at the origin as vxcval=vcoul +const1 + i*const2
      av(2)=av(2)+(vxcval(1)-vxc(1))/cl
      do 50 i=1,iwkb
  50  dv(i)=vxcval(i)/cl
c     keep dv=vxc/cl above iwkb

      nter=0
 
c     angular coefficients 
      call muatcc(xnval)

c     no orthogonalization needed. Looking for g.f., not w.f.
c     if (ipl.ne.0) call ortdac (ikap,ps,qs,aps,aqs)
c     ortdac orthogonalizes photoelectron orbital to core orbitals
c     have to use exchange 5 card to exit here; also want vxc=vxcval
c     if (ncycle .eq. 0) go to 999

c     iteration over the number of cycles
 101  continue
         nter=nter+1
c        calculate exchange potential
         jriwkb = min (jri, iwkb)
         call potex( ps, qs, aps, aqs, jriwkb, p2)

c        resolution of the dirac equation
         if (irr.lt.0) then
            call solout (p2, fl(norb), aps(1), aqs(1), ikap, rmt,
     1        jri, nmax(norb), ic3, vm, iwkb)
         else
            call solin (p2, fl(norb), pu, qu, ikap, rmt,
     1        jri, nmax(norb), ic3, vm, iwkb)
         endif

c     no orthogonalization needed. Looking for g.f., not w.f.
c        if (ipl.ne.0) call ortdac (ikap,gg,gp,ag,ap)

c        acceleration of the convergence 
         scc(norb)=1.0d0
         do 151 i=1,idim
            ps(i)=gg(i)
 151        qs(i)=gp(i)
         do 155 i=1,ndor
            aps(i) =ag(i) 
 155        aqs(i) =ap(i) 

      if (nter.le.ncycle) go to 101

 999  if (numerr .eq. 0) then
        if (irr.lt.0 ) then
cc        need pu, qu for regular solution
cc        want to have vxc(jri)-smooth and vxc(jri+1)=v_mt
cc        assume no exchange beyond jri 
           vu=vxc(jri+1)
           call flatv 
     1     (ri(jri), rmt, ps(jri), qs(jri), p2, vu, ikap, pu, qu)
           jlast = nmax(norb)
c          jlast might change on very rare occasion
        endif

      else
        call par_stop('error in dfovrg.f')
      endif

      return
      end

      subroutine flatv (r1, r2, p1, q1, en, vav, ikap, p2, q2)
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
c     solution of Dirac equation for flat potential for ikap is known
c     exactly (see e.g. in Loucks T.L. eq. 4-19)
c     given p1 and q1 at point r1 this subrotuine finds p2, q2 at r2
c     for given energy(en) and average potential (vav)
c     en and vav in hartrees
      external besjn, atancc

      complex*16 p1, q1, en, vav, p2, q2
      complex*16 ck, xkr, jl(ltot+2), nl(ltot+2), a,b, factor 

c     initialize staff
      ck = sqrt(2*(en-vav) + (alphfs*(en-vav))**2)
      xkr = ck*r1
      if (ikap.lt.0) then
        isign = -1
        lp = -ikap - 1
        lq = lp + 1
      else
        isign = 1
        lp = ikap
        lq = lp - 1
      endif
      a = ck * alphfs
      factor = isign*a/(1+sqrt(1+a**2))

c     find a and b that p1 = r1*(a*jl+b*nl), q1=factor*r1*(a*jl'+b*nl')
      call besjn (xkr, jl, nl)
      a = isign*ck*xkr* (p1*nl(lq+1) - q1*nl(lp+1)/factor)
      b = isign*ck*xkr* (q1*jl(lp+1)/factor - p1*jl(lq+1))

c     get values at r2
      xkr = ck * r2
      call besjn (xkr, jl, nl)
      p2 =  r2 * (jl(lp+1)*a + nl(lp+1)*b)
      q2 =  r2* factor * (jl(lq+1)*a + nl(lq+1)*b)

      return
      end

      subroutine diff (v, dr, kap, cl, dx, n, vm)
c     calculate  vm(i)=(dV/dx)*r(i)*(kap+1)/cl
c     needed for c3 term to calculate j-average phase shift
c     ref. koelling,harmon j.phys.c,3107(1977). eq.14
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
      
      complex*16 v(n), vm(n), vt(nrptx)
      dimension dr(n)
      do 5 i = 1,n
 5    vt(i) = v(i) * dr(i)**2

      vm(1)=((6.0*vt(2)+6.66666666667*vt(4)+1.2*vt(6))-(2.45*vt(1)+7.
     1 5*vt(3)+3.75*vt(5)+.166666666667*vt(7)))/dx
      vm(2)=((6.0*vt(3)+6.66666666667*vt(5)+1.2*vt(7))-(2.45*vt(2)+7.
     1 5*vt(4)+3.75*vt(6)+.166666666667*vt(8)))/dx
      nm2=n-2
      do 10 i=3,nm2
   10 vm(i)=((vt(i-2)+8.0*vt(i+1))-(8.0*vt(i-1)+vt(i+2)))/12.0/dx
      vm(n-1)=(vt(n)-vt(n-2))/(2.0*dx)
      vm(n)=(vt(n-2)*.5-2.0*vt(n-1)+1.5*vt(n))/dx

      do 20 i = 1,n
 20   vm(i) = (vm(i)-2*vt(i))/dr(i) *(kap+1.0)/cl
      return
      end
      complex*16 function dsordc(j,a,dg,dp,ag,ap)
c              * calculation of overlap integrals*
c        integration by simpson method of the   hg*(r**0)
c        hg(l)=dg(l)*cg(l,j)+dp(l)*cp(l,j)
c                cg,cp(l,j)  orbital j
c        a is such that dg,dp or hg following the case
c        behave at the origin as cte*r**a
c        the development limits at the origin (used for calculation
c        of integral form 0 to dr(1) ) of functions dg,dp and hg are
c        supposed to be in blocks ag,ap and chg respectively
c        this program uses   aprdec
c
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
      complex*16 aprdec
      common/dff/ cg(nrptx,30), cp(nrptx,30), bg(10,30), bp(10,30),
     1              fl(30), fix(30), ibgp
      complex*16 dg(nrptx),ag(10),dp(nrptx),ap(10)
      complex*16 hg(nrptx),chg(10)
c     common/ratom1/xnel(30),en(30),scc(30),scw(30),sce(30),
c    1   nq(30),kap(30),nmax(30)
      common/tabtec/hx,dr(nrptx),test1,test2,ndor,np,nes,method,idim
      dimension bgj(10),bpj(10)

c        construction of the array hg
      do  15 l= 1,ibgp
        bgj(l) = bg(l,j)
 15     bpj(l) = bp(l,j)

      do 221 l=1,idim
 221  hg(l)=dg(l)*cg(l,j)+dp(l)*cp(l,j)
      b=a+fl(j)
      do 241 l=1,ndor
 241     chg(l) = aprdec(ag,bgj,l) + aprdec(ap,bpj,l)
 
c        integration of the hg
      dsordc = (0.0d0, 0.0d0)
      do 305 l=1,idim
 305     hg(l)=hg(l)*dr(l)
      do 311 l=2,idim,2
 311     dsordc=dsordc+hg(l)+hg(l)+hg(l+1)
      dsordc=hx*(dsordc+dsordc+hg(1)-hg(idim))/3.0d0
c        integral from 0 to dr(1)
      do 331 l=1,ndor
         b=b+1.0d 00
 331     dsordc=dsordc+chg(l)*(dr(1)**b)/b
      return
      end
      subroutine inmuac (ihole, xionin, iunf, ikap)
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
      common/dff/cg(nrptx,30),cp(nrptx,30),bg(10,30),bp(10,30),fl(30),
     1    fix(30), ibgp
      common/itescf/testy,rap(2),teste,nz,norb,norbsc
c the meaning of common variables is described below
      common/ratom1/xnel(30),en(30),scc(30),scw(30),sce(30),
     1nq(30),kap(30),nmax(30)
c en one-electron energies
c scc factors for acceleration of convergence
c scw precisions of wave functions
c sce precisions of one-electron energies
c nmax number of tabulation points for orbitals
      common/scrhf1/eps(435),nre(30),ipl
c eps non diagonal lagrange parameters
c nre distingue: - the shell is closed (nre <0)
c                  the shell is open (nre>0)
c                - the orbitals in the integral rk if abs(nre) > or =2
c ipl define the existence of lagrange parameters (ipl>0)
      common/snoyac/dvn(nrptx),anoy(10),nuc
c dvn nuclear potential
c anoy development coefficients at the origin of nuclear potential
c this development is supposed to be written anoy(i)*r**(i-1)
c nuc index of nuclear radius (nuc=1 for point charge)
      common/tabtec/hx,dr(nrptx),test1,test2,ndor,np,nes,method,idim
      dimension xnval(30), iorb(-4:3)
      data nucm/11/

      testy=10.**(-5)
c testy precision for the wave functions

      call getorb (nz, ihole, xionin, iunf, norb, norbsc, iorb,
     1            iholep, nq, kap, xnel, xnval, en)
c     don't need xmag here, so use en as a dummy

      ipl=0
      do 40 i=1,norb
         en(i) = 0.d0
         nre(i)=-1
         llq= abs(kap(i))
         l=llq+llq
c       find last tabulation point
         nmax(i)=0
         do 100  j = idim, 1, -1
            if ( abs(cg(j,i)) .ge. 1.0d-11 .or.
     1           abs(cp(j,i)) .ge. 1.0d-11 )  then
               nmax(i) = j
               goto 16
            endif
  100    continue
   16    continue

         scc(i)=0.3
         if (xnel(i) .lt. l)  nre(i)=1
         if (ikap.eq.kap(i)) ipl=ipl+1
  40  continue
      norbsc=norb
      norb = norb+1
      xnel(norb)=1
      kap(norb)=ikap
      nq(norb) =9
c nz atomic number     noi ionicity (nz-number of electrons)
c norb number of orbitals
c xnel(i) number of electrons on orbital i.
      nuc=nucm
c nuc number of points inside nucleus (11 by default)

      return
      end
      subroutine intout (en,i0, kap,max0,ic3,vm)
c                  resolution of the dirac equation
c                   p' - kap*p/r = - ( en/cl-v )*g - eg/r
c                   g' + kap*g/r = ( 2*cl+en/cl-v )*p + ep/r
c at the origin v approximately is -z/(r*cl) due to the point nucleus
c en one-electron energy in atomic units and negative
c at the origin of the large(small)component
c kap quantum number kappa
c max0 the last point of tabulation of the wave function
 
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
      parameter (npi=6, test=1.0d+5)
      complex*16 en,c3,vmh
      complex*16 gg,ag,gp,ap,dv,av,eg,ceg,ep,cep, vm(nrptx)
      common/comdic/cl,dz,gg(nrptx),ag(10),gp(nrptx),ap(10),dv(nrptx),
     1   av(10),eg(nrptx),ceg(10),ep(nrptx),cep(10)

      complex*16 ec,eph,egh,f,g,ac,bc,acp,bcp,dg,dp, dv1,dv2,vh
      complex*16 dg2, dp2, dg3, dp3, dg4, dp4
      dimension dg(npi), dp(npi)

c gg,gp -output, dv,eg,ep - input
c
c cl speed of light (approximately 137.037 in atomic units)
c dz nuclear charge
c gg (gp) large (small) component
c dv direct potential (v)     eg and ep exchange potentials
c ag,ap,av,ceg and cep are respectively the
c development coefficients for gg,gp,dv,eg and ep
c
      common/tabtec/hx,dr(nrptx),test1,test2,ndor,np,nes,method,idim
c hx exponential step
c dr radial mesh
c test1,test2,nes,method are dummy.
c  ndor number of terms for the developments at the origin
c np maximum number of the tabulation points
c idim dimension of the block dr


      ccl=cl+cl
      exphx = exp (hx/2)
      ihard = 0
      ec=en/cl
 
c            solution of the inhomogenios dirac equation
c gg gp initially exch. terms, at the time of return are wave functions
c ag and ap development coefficients of  gg and gp
c en one-electron energy
c fl power of the first development term at the origin

c     runge-kutta for first npi points
      i = i0
      j=1
      f = (ec - dv(i))*dr(i)
      g = f + ccl * dr(i)
      c3 = ic3*vm(i)/g**2
      dg(j) = hx * (g*gp(i) - kap*gg(i) + ep(i))
      dp(j) = hx * (kap*gp(i) - (f-c3)*gg(i) - eg(i))

 44   continue
      if (i.ge.max0) goto 999
      ac = gg(i) + 0.5d0 * dg(j)
      bc = gp(i) + 0.5d0 * dp(j)
      rh = dr(i) *exphx
c     find potential and exchange terms between 2 points
c     use linear interpolation with imp. nonlinearity correction
      xm1 = (dr(i+1)-rh) / (dr(i+1)-dr(i))
      xm2 = (rh - dr(i)) / (dr(i+1)-dr(i))
      if (dble(av(1)) .lt. 0.0 .and. i0.eq.1) then
c        point nucleus case
c        important nonlinearity from z/r term
         dv1 = dv(i) - av(1)/dr(i)
         dv2 = dv(i+1) - av(1)/dr(i+1)
         vh = dv1*xm1 + dv2*xm2
         vh = vh + av(1)/rh
         vmh = (xm1*vm(i)*dr(i) +xm2*vm(i+1)*dr(i+1))/rh
      elseif (i0.eq.1) then
c        finite nucleus
c        important nonlinearity from z*r**2 term
         dv1 = dv(i) - av(4)*dr(i)**2
         dv2 = dv(i+1) - av(4)*dr(i+1)**2
         vh = (dv1*(dr(i+1)-rh)+dv2*(rh-dr(i))) / (dr(i+1)-dr(i))
         vh = vh + av(4)*rh**2
         vmh = (xm1*vm(i)/dr(i)**2 +xm2*vm(i+1)/dr(i+1)**2)*rh**2
      else
c        outward integration of irregular solution near jri
         vh = dv(i)*xm1 + dv(i+1)*xm2
         vmh = xm1*vm(i) +xm2*vm(i+1)
      endif
      eph = ep(i) * xm1 + ep(i+1) * xm2
      egh = eg(i) * xm1 + eg(i+1) * xm2

      f = (ec - vh)*rh
      g = f + ccl * rh
      c3 = ic3*vmh/g**2
      dg2 = hx * (g*bc - kap*ac + eph)
      dp2 = hx * (kap*bc - (f-c3)*ac - egh)
      ac = ac + 0.50*(dg2-dg(j))
      bc = bc + 0.50*(dp2-dp(j))
      dg3 = hx * (g*bc - kap*ac + eph)
      dp3 = hx * (kap*bc - (f-c3)*ac - egh)
      ac = ac + dg3 - 0.50*dg2
      bc = bc + dp3 - 0.50*dp2

      i=i+1
      j=j+1
      f = (ec - dv(i))*dr(i)
      g = f + ccl * dr(i)
      c3 = ic3*vm(i)/g**2
      dg4 = hx * (g*bc - kap*ac + ep(i))
      dp4 = hx * (kap*bc - (f-c3)*ac - eg(i))
      gg(i) = gg(i-1)+(dg(j-1) + 2.0*(dg2+dg3)+dg4)/6.0
      gp(i) = gp(i-1)+(dp(j-1) + 2.0*(dp2+dp3)+dp4)/6.0
      dg(j) = hx * (g*gp(i) - kap*gg(i) + ep(i))
      dp(j) = hx * (kap*gp(i) - (f-c3)*gg(i) - eg(i))
      if (j.lt.npi) goto 44

c     scale derivatives for milne method
      do 51 i = 1,npi
        dg(i) = dg(i)/hx
 51     dp(i) = dp(i)/hx

c     integration of the inhomogenious system
      a1 = hx * 3.3
      a2 = -hx * 4.2
      a3 = hx * 7.8
      a4 = hx * 14.0/45.0
      a5 = hx * 64.0/45.0
      a6 = hx * 24.0/45.0
      do 55 i = npi+i0-1,max0-1
         nit = 0
c        predictor
         acp=gg(i-5)+a1*(dg(npi)+dg(npi-4))+a2*(dg(npi-1)+dg(npi-3))
     1       +a3*dg(npi-2)
         bcp=gp(i-5)+a1*(dp(npi)+dp(npi-4))+a2*(dp(npi-1)+dp(npi-3))
     1       +a3*dp(npi-2)
c        ac,bc -corrector w/o contribution from derivatives at i+1
         ac=gg(i-3)+a4*dg(npi-3)+a5*(dg(npi)+dg(npi-2))+a6*dg(npi-1)
         bc=gp(i-3)+a4*dp(npi-3)+a5*(dp(npi)+dp(npi-2))+a6*dp(npi-1)
         do 61 j=1,npi-1
            dg(j)=dg(j+1)
 61         dp(j)=dp(j+1)
         f=(ec-dv(i+1))*dr(i+1)
         g=f+ccl*dr(i+1)
         c3 = ic3*vm(i+1)/g**2
 64      dg(npi)=g*bcp-kap*acp+ep(i+1)
         dp(npi)=kap*bcp-(f-c3)*acp-eg(i+1)
c        corrected values
         gg(i+1)=ac+a4*dg(npi)
         gp(i+1)=bc+a4*dp(npi)
         if ( abs(test*(gg(i+1)-acp)) .gt. abs(gg(i+1)) .or.
     1        abs(test*(gp(i+1)-bcp)) .gt. abs(gp(i+1)) ) then
c           test failed
            if (nit.lt.40) then
               acp = gg(i+1)
               bcp = gp(i+1)
               nit = nit + 1
               goto 64
            else
               ihard = ihard+1
            endif
         endif
 55   continue

 999  do 741 i=max0+1,np
         gg(i)=0.0d 00
 741     gp(i)=0.0d 00

      return
      end
      subroutine muatcc(xnval) 
c               * angular coefficients *
c        sous programmes utilises  cwig3j
c
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
      dimension xnval(30)
      common/itescf/testy,rap(2),teste,nz,norb,norbsc
      common/mulabc/afgkc
      dimension afgkc(-ltot-1:ltot,30,0:3)
      common/ratom1/xnel(30),en(30),scc(30),scw(30),sce(30),
     1nq(30),kap(30),nmax(30)
c#mn
       external cwig3j

      do 511 i=-ltot-1,ltot
      do 511 j=1,30
      do 511 k=0,3
 511  afgkc(i,j,k)=0.0d 00
 601  do 701 ikap=-ltot-1,ltot
         if (ikap .eq. 0) go to 701
         li= abs(ikap)*2-1
         do 700 j=1,norb-1
            lj= abs(kap(j))*2-1
            kmax=(li+lj)/2
            kmin= abs(li-lj)/2
            if ((ikap*kap(j)).lt.0) kmin=kmin+1
            if (xnval(j) .gt. 0.0d0) goto 700
c calculate b_k(i,j)
            do 675 k = kmin, kmax,2
               index=(k-kmin)/2
               afgkc(ikap,j,index)=xnel(j)*(cwig3j(li,k*2,lj,1,0,2)**2)
 675        continue
 700     continue
 701  continue
      return
      end
      subroutine nucdec (av,dr,dv,dz,hx,nuc,np,ndor,dr1)
c        * construction of nuclear potential *
c av coefficients of the development at the origin of nuclear potential
c dr  tabulation points
c dv  nuclear potential 
c dz  nuclear charge 
c hx  exponential step
c nuc index of the nuclear radius
c np  number of tabulation points
c ndor number of the coefficients for development at the origin
c the declared below arguments are saved, dr1 is the first
 
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
      dimension av(10),dr(nrptx),dv(nrptx),at(nrptx)

c    specify atomic mass and thickness of nuclear shell
c a atomic mass (negative or null for the point charge)
c epai parameter of the fermi density distribution
c (negative or null for uniform distribution), which is
c       cte / (1. + exp((r-rn)/epai) )
c with nuclear radius rn= 2.2677e-05 * (a**(1/3))

c calculate radial mesh
      a = 0.0
      epai = 0.0

      if (a.le.1.0d-01) then
         nuc=1
      else
         a=dz*(a**(1./3.))*2.2677d-05
         b=a/ exp(hx*(nuc-1))
         if (b.le.dr1) then
            dr1=b
         else
            b=log(a/dr1)/hx
            nuc=3+2*int(b/2.0)
            if (nuc.ge.np) call par_stop('dr1 too small')
c           index of atomic radius larger than dimension of dr
            dr1=a*exp(-(nuc-1)*hx)
         endif
      endif

      dr(1)=dr1/dz
      do 181 l=2,np
 181  dr(l)=dr(1)* exp(hx*(l-1))

      if (ndor.lt.5) then
c       * it should be at least 5 development coefficients
         call wlog('stopped in programm nucdec, ndor should be > 4.')
         call par_stop('NUCDEC-1')
      endif
c  calculate nuclear potential on calculated radial mesh
      do 11 i=1,ndor
 11      av(i)=0.0d 00
      if (epai.le.0.0) then
         do 15 i=1,np
 15         dv(i)=-dz/dr(i)
         if (nuc.le.1) then
            av(1)=-dz
         else
            av(2)=-3.0d 00*dz/(dr(nuc)+dr(nuc))
            av(4)=-av(2)/(3.0d 00*dr(nuc)*dr(nuc))
            l=nuc-1
            do 25 i=1,l
 25            dv(i)=av(2)+av(4)*dr(i)*dr(i)
         endif
      else
         b= exp(-dr(nuc)/epai)
         b=1.0d 00/(1.0d 00+b)
         av(4)=b
         av(5)=epai*b*(b-1.0d 00)
         if (ndor.le.5) go to 45
         at(1)=1.0d 00
         at(2)=1.0d 00
         nf=1
         do 41 i=6,ndor
            n=i-4
            nf=n*nf
            dv(1)=n*at(1)
            n1=n+1
            dv(n1)=1.0d 00
            do 35 j=2,n
 35         dv(j)=(n-j+2)*at(j-1)+(n-j+1)*at(j)
            do 37 j=1,n1
               m=n+1-j
               l=1
               if (mod(j,2).eq.0) l=-l
               av(i)=av(i)+l*dv(j)*(b**m)
 37            at(j)=dv(j)
 41         av(i)=b*av(i)*(epai**n)/nf
 45      do 47 i=1,np
            b=1.0d 00+ exp((dr(i)-dr(nuc))/epai)
            if ((b*av(4)).gt.1.0d+15) go to 51
            dv(i)=dr(i)*dr(i)*dr(i)/b
 47         l=i
 51      if (l.ge.(np-1)) l=np-2
         k=l+1
         do 55 i=k,np
 55         dv(i)=0.0d 00
         at(1)=0.0d 00
         at(2)=0.0d 00
         k=2
         do 61 i=4,ndor
            k=k+1
            do 58 j=1,2
 58         at(j)=at(j)+av(i)*(dr(j)**k)/k
            av(i)=av(i)/(k*(k-1))
 61         av(2)=av(2)+av(i)*(dr(1)**k)
         a=hx/2.4d+01
         b=a*1.3d+01
         k=l+1
         do 71 i=3,k
 71      at(i)=at(i-1)+b*(dv(i-1)+dv(i))-a*(dv(i-2)+dv(i+1))
         dv(l)=at(l)
         do 75 i=k,np
 75      dv(i)=dv(l)
         e= exp(hx)
         c=1.0d 00/(e*e)
         i=l-1
 83      dv(i)=dv(i+1)/e+b*(at(i+1)/e+at(i))-a*(at(i+2)*c+at(i-1)*e)
         i=i-1
         if (i-1) 85,85,83
 85      dv(1)=dv(3)*c+hx*(at(1)+4.0d 00*at(2)/e+at(3)*c)/3.0d 00
         av(2)=(av(2)+dv(1))/dr(1)
         a=-dz/dv(l)
         do 95 i=4,ndor
 95      av(i)=-a*av(i)
         av(2)=a*av(2)
         do 97 i=1,np
 97      dv(i)=a*dv(i)/dr(i)
      endif

      return
      end
      subroutine ortdac(ikap,ps,qs,aps,aqs)
c        * orthogonalization by the schmidt procedure*
c the ia orbital is orthogonalized toa all orbitals of the same
c symmetry if ia is positive, otherwise all orbitals of the same
c symmetry are orthogonalized
c        this program uses dsordc
 
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
      complex*16 dsordc
      complex*16 ps(nrptx), qs(nrptx), aps(10),aqs(10)
      common/dff/ cg(nrptx,30), cp(nrptx,30), bg(10,30), bp(10,30),
     1             fl(30), fix(30), ibgp
      common/itescf/testy,rap(2),teste,nz,norb,norbsc
      common/ratom1/xnel(30),en(30),scc(30),scw(30),sce(30),
     1   nq(30),kap(30),nmax(30)
      common/tabtec/hx,dr(nrptx),test1,test2,ndor,np,nes,method,idim
      complex*16 a
 
      do 51 j=1,norb-1
         if (kap(j).ne.ikap .or. xnel(j).le.0) go to 51
         a = dsordc(j,fl(norb),ps,qs,aps,aqs)
         do 41 i=1,idim
            ps(i)=ps(i)-a*cg(i,j)
 41         qs(i)=qs(i)-a*cp(i,j)
         do 42 i=1,ndor
            aps(i)=aps(i)-a*bg(i,j)
 42         aqs(i)=aqs(i)-a*bp(i,j)
 51   continue
      return
      end
      subroutine potdvp
c     this programm uses aprdep,multrk,yzkrdf
c     to calculate potential development coefficients
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
      common/dff/ cg(nrptx,30), cp(nrptx,30), bg(10,30), bp(10,30),
     1              fl(30), fix(30), ibgp
      complex*16 dg,ag,dp,ap,dv,av,eg,ceg,ep,cep
      common/comdic/cl,dz,dg(nrptx),ag(10),dp(nrptx),ap(10),dv(nrptx),
     2         av(10),eg(nrptx),ceg(10),ep(nrptx),cep(10)
c     dg,dp to get data from yzkrdf, dv,eg,ep -output for soldir
      common/itescf/testy,rap(2),teste,nz,norb,norbsc
      common/ratom1/xnel(30),en(30),scc(30),scw(30),sce(30),
     1nq(30),kap(30),nmax(30)
      common/snoyac/dvn(nrptx),anoy(10),nuc
      common/tabtec/hx,dr(nrptx),test1,test2,ndor,np,nes,method,idim
      dimension bgj(10),bpj(10)
c#mn
       external aprdep

      do 9 i=1,10
 9       av(i)=anoy(i)
 
c     calculate density development coefficients
      do 31 i=1,ndor
 31   ag(i)=0.0d 00
      do 51 j=1,norb-1
         do 33 i = 1,10
            bgj(i) = bg(i,j)
 33         bpj(i) = bp(i,j)
         n=2* abs(kap(j))
         l=ndor+2-n
         if (l.le.0) go to 51
         do 41 i=1,l
            m=n-2+i
 41         ag(m)=ag(m)+xnel(j)*(aprdep(bgj,bgj,i)+
     1            aprdep(bpj,bpj,i))*fix(j)**2
 51   continue

c     transform density coefficients into ones for potential
      ap(1)=0.0d 00 
      do 15 i=1,ndor
         ag(i)=ag(i)/(i+2)/(i+1)
         ap(1)=ap(1)+ag(i)*dr(1)**(i+1)
 15   continue

      do 61 i=1,ndor
         l=i+3
         if (l.gt.ndor) go to 61
         av(l)=av(l)-ag(i)
 61   continue
c     av(2)=avoy(2) + ap(1)+(vxcvzl(1)-dvn(1)) in order 
c     to have sum av(i)*dr(1)**(i-2)=vxcval(1)
      av(2)=av(2)+ap(1)
 
c addition of nuclear potential and division of potentials and
c       their development limits by speed of light
      do 527 i=1,10
 527     av(i)=av(i)/cl
      return
      end
      subroutine potex( ps, qs, aps, aqs, jri, p2)
c        this programm uses bkeato,aprdec,multrk,yzkrdc
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
      complex*16 aprdec, p2
      complex*16 ps(nrptx),qs(nrptx),aps(10),aqs(10)
      common/dff/cg(nrptx,30),cp(nrptx,30),bg(10,30),bp(10,30),
     1             fl(30), fix(30), ibgp
      complex*16 dg,ag,dp,ap,dv,av,eg,ceg,ep,cep
      common/comdic/cl,dz,dg(nrptx),ag(10),dp(nrptx),ap(10),dv(nrptx),
     2          av(10),eg(nrptx),ceg(10),ep(nrptx),cep(10)
c     dg,dp to get data from yzkrdc, dv,eg,ep -output for soldir
      common/itescf/testy,rap(2),teste,nz,norb,norbsc
      common/ratom1/xnel(30),en(30),scc(30),scw(30),sce(30),
     1nq(30),kap(30),nmax(30)
      common/tabtec/hx,dr(nrptx),test1,test2,ndor,np,nes,method,idim
      common/mulabc/afgkc
      dimension afgkc(-ltot-1:ltot,30,0:3)
      dimension bgj(10),bpj(10)
c#mn
       external aprdec
 
c     ia=norb
      jia=2* abs(kap(norb))-1
      do 9 i=1,10
         cep(i)=0.0d 00
 9       ceg(i)=0.0d 00
      do 11 i=1,idim
         ep(i)=0.0d 00
 11      eg(i)=0.0d 00
 
c     exchange terms
      do 201 j=1,norb-1
 105     jj=2* abs(kap(j))-1
         kma=(jj+jia)/2
         k= abs(jj-kma)
         if ((kap(j)*kap(norb)).lt.0) k=k+1
         kmin = k
c        kma=min(kma,15)
c        if (k.lt.kma) goto 201

c111     a=bkeato(j,ia,k)/xnel(ia)
 111     a=afgkc(kap(norb),j,(k-kmin)/2)
         if (a.eq.0.0d 00) go to 151
         call yzkrdc (j,k,fl(norb),ps,qs,aps,aqs, p2, norb)
         do 121 i=1,idim
            eg(i)=eg(i)+a*dg(i)*cg(i,j)
 121        ep(i)=ep(i)+a*dg(i)*cp(i,j)
         n=k+1+ abs(kap(j))- abs(kap(norb))
c         differrent for irregular solution
         if (fl(norb) .lt.0.0) n=k+1+ abs(kap(j)) + abs(kap(norb))
         if (n.gt.ndor) go to 141
         do 135 i=n,ndor
            ceg(i)=ceg(i)+bg(i+1-n,j)*a*ap(1)*fix(j)/fix(norb)
 135        cep(i)=cep(i)+bp(i+1-n,j)*a*ap(1)*fix(j)/fix(norb)
 141     i=2* abs(kap(j))+1
         if (i.gt.ndor) go to 151
         do 143 ix = 1,10
            bgj(ix) = bg(ix,j)
 143        bpj(ix) = bp(ix,j)
         do 145 n=i,ndor
            nx = n + 1 - i
            ceg(n) = ceg(n) - a * aprdec(ag,bgj,nx)*fix(j)**2
 145        cep(n) = cep(n) - a * aprdec(ag,bpj,nx)*fix(j)**2
 151     k=k+2
         if (k.le.kma) go to 111
 201  continue
 
c    division of potentials and
c    their development limits by speed of light
      do 527 i=1,ndor
         cep(i)=cep(i)/cl
 527     ceg(i)=ceg(i)/cl
      do 531 i=1,jri
         ep(i)=ep(i)/cl
 531     eg(i)=eg(i)/cl
      do 532 i=jri+1,nrptx
         ep(i)=0.0d0
 532     eg(i)=0.0d0

      return
      end
      subroutine solin (en,fl,agi,api,kap,rmt,jri,imax,ic3,vm, iwkb)
c                  resolution of the dirac equation
c                   p' - kap*p/r = - ( en/cl-v )*g - eg/r
c                   g' + kap*g/r = ( 2*cl+en/cl-v )*p + ep/r
c at the origin v approximately is -z/(r*cl) due to the point nucleus
c en one-electron energy in atomic units and negative
c fl power of the first term in development at the origin
c agi (api) initial values of the first development coefficient
c at the origin of the large(small)component
c kap quantum number kappa
c imax the last point of tabulation of the wave function

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
      parameter (npi=6, test=1.0d+5)
      complex*16 en,agi,api,c3,vmh
      complex*16 gg,ag,gp,ap,dv,av,eg,ceg,ep,cep, vm(nrptx)
      common/comdic/cl,dz,gg(nrptx),ag(10),gp(nrptx),ap(10),dv(nrptx),
     1   av(10),eg(nrptx),ceg(10),ep(nrptx),cep(10)

      complex*16 ec,eph,egh,f,g,ac,bc,acp,bcp,dg,dp, vh
      complex*16 dg2, dp2, dg3, dp3, dg4, dp4
      dimension dg(npi), dp(npi)

c gg,gp -output, dv,eg,ep - input
c
c cl speed of light (approximately 137.037 in atomic units)
c dz nuclear charge
c gg (gp) large (small) component
c dv direct potential (v)     eg and ep exchange potentials
c ag,ap,av,ceg and cep are respectively the
c development coefficients for gg,gp,dv,eg and ep
c
      common/tabtec/hx,dr(nrptx),test1,test2,ndor,np,nes,method,idim
c hx exponential step
c dr radial mesh
c test1,test2,nes,method are dummy.
c  ndor number of terms for the developments at the origin
c np maximum number of the tabulation points
c idim dimension of the block dr
      complex*16 jl(0:ltot+1), hl(0:ltot+1), xkmt, ck, dum1, factor
      external besjh

      ccl=cl+cl
      ihard = 0
      ec=en/cl
      do 115 i=2,ndor
         ag(i)=0.0d0
 115     ap(i)=0.0d0
c     integration of the inhomogenious system
c     no need in normalization, since we can use 
c     normalization agi=ag(1)=const
 
c            solution of the inhomogenios dirac equation
c gg gp initially exch. terms, at the time of return are wave functions
c ag and ap development coefficients of  gg and gp
c en one-electron energy
c fl power of the first development term at the origin
c agi (api) initial values of the first development coefficients
c at the origin of a large (small) component
 
c     started with h_l above jri inside dfovrg
      vmh = cl * dv(jri+1)
      ck = sqrt(2*(en-vmh) + (alphfs*(en-vmh))**2)
      il = abs(kap)
      if (kap.lt. 0) il = il - 1
      ilp = il - 1
      if (kap .lt. 0) ilp = il + 1
      ilx = il+1
      if (ilp.gt.il) ilx=ilp+1
      xsign = -1.d0
      if (kap.gt.0) xsign = 1.d0
      factor = ck*alphfs
      factor = xsign * factor/(1+sqrt(1+factor**2))
      dum1 = 1/ sqrt(1+factor**2)

      iflat = min ( jri, iwkb)
      do i = jri, imax
        j= iflat + npi - i
        xkmt = ck * dr(i)
        call besjh( xkmt, ilx, jl, hl)
        gg(i) = hl(il) * dr(i) * dum1
        gp(i) = hl(ilp) * dr(i) * dum1 * factor
        if (j.gt.0) then
          f = (ec - dv(i))*dr(i)
          g = f + ccl * dr(i)
          c3 = ic3*vm(i)/g**2
          dg(j) = -(  g*gp(i) - kap*gg(i) )
          dp(j) = -(  kap*gp(i) - (f-c3)*gg(i) )
c         neglect exchage term outside jri
c         dg(j) = -(  g*gp(i) - kap*gg(i) + ep(i) )
c         dp(j) = -(  kap*gp(i) - (f-c3)*gg(i) - eg(i) )
        endif
      enddo

c     use flatv between iwkb and jri
      do i = jri-1, iflat, -1
         j= iflat + npi - i
         if (i.eq.iwkb) then
            eph = cl* ( 3*dv(iwkb+1) - dv(iwkb+2)) /2
            if (iwkb.eq.jri-1) eph=  cl* (dv(i) + dv(i+1)) /2
         else
            eph = cl* (dv(i) + dv(i+1)) /2
         endif
         if (ic3.gt.0) then
           rav = (dr(i)+dr(i+1)) / 2
           ec = rav**3 * ( ccl+ (en - eph) / cl )**2
           eph = eph + ic3 * cl / ec * (vm(i) + vm(i+1)) / 2
         endif
         call flatv( dr(i+1), dr(i), gg(i+1), gp(i+1), en, eph, kap,
     1               gg(i), gp(i))
         if (j.gt.0) then
          f = (ec - dv(i))*dr(i)
          g = f + ccl * dr(i)
          c3 = ic3*vm(i)/g**2
          dg(j) = -(  g*gp(i) - kap*gg(i) + ep(i) )
          dp(j) = -(  kap*gp(i) - (f-c3)*gg(i) - eg(i) )
         endif
      enddo

c     integration of the inhomogenious system
      a1 = hx * 3.3
      a2 = -hx * 4.2
      a3 = hx * 7.8
      a4 = hx * 14.0/45.0
      a5 = hx * 64.0/45.0
      a6 = hx * 24.0/45.0
c     do 55 i = jri - npi + 1 , 2, -1
      do 55 i = iflat, 2, -1
         nit = 0
c        predictor
         acp=gg(i+5)+a1*(dg(npi)+dg(npi-4))+a2*(dg(npi-1)+dg(npi-3))
     1       +a3*dg(npi-2)
         bcp=gp(i+5)+a1*(dp(npi)+dp(npi-4))+a2*(dp(npi-1)+dp(npi-3))
     1       +a3*dp(npi-2)
c        ac,bc -corrector w/o contribution from derivatives at i+1
         ac=gg(i+3)+a4*dg(npi-3)+a5*(dg(npi)+dg(npi-2))+a6*dg(npi-1)
         bc=gp(i+3)+a4*dp(npi-3)+a5*(dp(npi)+dp(npi-2))+a6*dp(npi-1)
         do 61 j=1,npi-1
            dg(j)=dg(j+1)
 61         dp(j)=dp(j+1)
         f=(ec-dv(i-1))*dr(i-1)
         g=f+ccl*dr(i-1)
         c3 = ic3*vm(i-1)/g**2
 64      dg(npi)= -( g*bcp-kap*acp+ep(i-1) )
         dp(npi)= -( kap*bcp-(f-c3)*acp-eg(i-1) )
c        corrected values
         gg(i-1)=ac+a4*dg(npi)
         gp(i-1)=bc+a4*dp(npi)
         if ( abs(test*(gg(i-1)-acp)) .gt. abs(gg(i-1)) .or.
     1        abs(test*(gp(i-1)-bcp)) .gt. abs(gp(i-1)) ) then
c           test failed
            if (nit.lt.40) then
               acp = gg(i-1)
               bcp = gp(i-1)
               nit = nit + 1
               goto 64
            else
               ihard = ihard+1
            endif
         endif
 55   continue

      do 741 i=imax+1,np
         gg(i)=0.0d 00
 741     gp(i)=0.0d 00
      ag(1)=gg(1)* dr(1)**(-fl)
      ap(1)=gp(1)* dr(1)**(-fl)

      return
      end
      subroutine solout(en, fl, agi, api, kap, rmt,
     1                  jri, max0, ic3, vm, iwkb)
c                  resolution of the dirac equation
c                   p' - kap*p/r = - ( en/cl-v )*g - eg/r
c                   g' + kap*g/r = ( 2*cl+en/cl-v )*p + ep/r
c at the origin v approximately is -z/(r*cl) due to the point nucleus
c en one-electron energy in atomic units and negative
c fl power of the first term in development at the origin
c agi (api) initial values of the first development coefficient
c at the origin of the large(small)component
c kap quantum number kappa
c max0 the last point of tabulation of the wave function
 
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
      parameter (npi=6, test=1.0d+5)
      parameter (ccl=2*alpinv, csq=ccl**2 )
      complex*16 en,agi,api
      complex*16 gg,ag,gp,ap,dv,av,eg,ceg,ep,cep, vm(nrptx)
      common/comdic/cl,dz,gg(nrptx),ag(10),gp(nrptx),ap(10),dv(nrptx),
     1   av(10),eg(nrptx),ceg(10),ep(nrptx),cep(10)

      complex*16 ec,eph,f,g

c gg,gp -output, dv,eg,ep - input
c
c cl speed of light (approximately 137.037 in atomic units)
c dz nuclear charge
c gg (gp) large (small) component
c dv direct potential (v)     eg and ep exchange potentials
c ag,ap,av,ceg and cep are respectively the
c development coefficients for gg,gp,dv,eg and ep
c
      common/tabtec/hx,dr(nrptx),test1,test2,ndor,np,nes,method,idim
c hx exponential step
c dr radial mesh
c test1,test2,nes,method are dummy.
c  ndor number of terms for the developments at the origin
c np maximum number of the tabulation points
c idim dimension of the block dr


c{#mn: g77 chokes on taking real of a double-precision complex
c      if (real(av(1)).lt.0.0d 00.and.kap.gt.0) api=-agi*(kap+fl)/av(1)
c      if (real(av(1)).lt.0.0d 00.and.kap.lt.0) api=-agi*av(1)/(kap-fl)
      if (dble(av(1)).lt.0.0d 00.and.kap.gt.0) api=-agi*(kap+fl)/av(1)
      if (dble(av(1)).lt.0.0d 00.and.kap.lt.0) api=-agi*av(1)/(kap-fl)
c#mn}
      ec=en/cl
      ag(1)=agi
      ap(1)=api
      do 115 i=2,ndor
         ag(i)=ceg(i-1)
 115     ap(i)=cep(i-1)
c     integration of the inhomogenious system
c     no need in normalization, since we can use 
c     normalization agi=ag(1)=const
 
c            solution of the inhomogenios dirac equation
c gg gp initially exch. terms, at the time of return are wave functions
c ag and ap development coefficients of  gg and gp
c en one-electron energy
c fl power of the first development term at the origin
c agi (api) initial values of the first development coefficients
c at the origin of a large (small) component
 
c     initial values for the outward integration
      if (ic3.eq.0) then
c       Desclaux power expansion
         do 35 j=2,ndor
            k=j-1
            a=fl+kap+k
            b=fl-kap+k
            eph=a*b+av(1)*av(1)
            f=(ec+ccl)*ap(k)+ap(j)
            g=ec*ag(k)+ag(j)
            do 31 i=1,k
               f=f-av(i+1)*ap(j-i)
 31            g=g-av(i+1)*ag(j-i)
 
            ag(j)=(b*f+av(1)*g)/eph
 35         ap(j)=(av(1)*f-a*g)/eph

         do  41 i = 1,1
            gg(i)=0.0d 00
            gp(i)=0.0d 00
         do 41 j=1,ndor
            a=fl+j-1
            b=dr(i)**a
            gg(i)=gg(i)+b*ag(j)
 41         gp(i)=gp(i)+b*ap(j)
      else
c        see fovrg.f in feff6, be aware of different units
         twoz = -dble(av(1)) * 2.0*cl
         rat1 = twoz/ccl
         rat2 = rat1**2
         rat3 = csq/twoz
         il = -kap
         if (kap.gt.0) il = kap+1
         l0 = il-1
         ag(1) = agi
         if (twoz.le.0.0) then
            ap(1) = -ec/(2.0*il+1.0)*dr(1)*ag(1)
            ag(2) = 0.0
            ap(2) = 0.0
            ag(3) = 0.0
            ap(3) = 0.0
         else
            ap(1) = (fl-il)*rat3*ag(1)
            ag(2) = (3.0*fl-rat2)/(2.0*fl+1.0) * ag(1)
            ap(2)= rat3*( (fl -l0)*ag(2) - ag(1) ) -ap(1)
            ag(3)=( (fl+3.0*il)*ag(2) - 3.0*l0*ag(1) + 
     1      (fl+il+3.0)/rat3*ap(2) ) /(fl+1.0)/4.0
            ap(3)=( rat3*(2.0*l0*(fl+2.0-il)-l0-rat2)*ag(2)
     1      - 3.0*l0*rat3*(fl+2.0-il)*ag(1) + (fl+3.0-2.0*il-rat2)
     2      *ap(2) ) /(fl+1.0)/4.0
            ap(1) = ap(1)/ccl
            ag(2)= ag(2)*rat3
            ap(2)= ap(2)*rat3/ccl
            ag(3)= ag(3)*rat3**2
            ap(3)= ap(3)*rat3**2/ccl
         endif
         gg(1)=dr(1)**fl * (ag(1)+dr(1)*(ag(2)+dr(1)*ag(3)))
         gp(1)=dr(1)**fl * (ap(1)+dr(1)*(ap(2)+dr(1)*ap(3)))
      endif

      i0=1
      iflat = min ( jri, iwkb)
      call intout (en, i0, kap, iflat, ic3, vm)

      do 100 i = iflat, max0-1
         if (i.eq.iwkb) then
            eph = cl* ( 3*dv(iwkb+1) - dv(iwkb+2)) /2
            if (iwkb.eq.jri-1) eph=  cl* (dv(i) + dv(i+1)) /2
         else
            eph = cl* (dv(i) + dv(i+1)) /2
         endif
         if (ic3.gt.0 .and. i.lt.jri) then
           rav = (dr(i)+dr(i+1)) / 2
           ec = rav**3 * ( ccl+ (en - eph) / cl )**2
           eph = eph + ic3 * cl / ec * (vm(i) + vm(i+1)) / 2
         endif
         call flatv( dr(i), dr(i+1), gg(i), gp(i), en, eph, kap,
     1               gg(i+1), gp(i+1))
  100 continue

      return
      end
      subroutine wfirdc (eph,kap,nmax,vxc,ps,qs,aps,aqs,irr,ic3,vm,
     1                   rmt,jri, iwkb) 
c     calculate photoelectron orbital using lda in dirac equation
c     cg (cp) large (small) radial components
c     bg (bp) development coefficients at the origin of cg (cp)
c     eph one-electron energy of photoelectron
c     fl power of the first term of development at the origin
c     kap quantum number "kappa"
c     nmax number of tabulation points for the orbitals
c     vxc  is initial lda potential for photoelectron
c     ibgp first dimension of the arrays bg and bp
c        this programmes utilises nucdec,dentfa,soldir et messer
 
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
      common/dff/cg(nrptx,30),cp(nrptx,30),bg(10,30),bp(10,30),
     1             fl(30), fix(30), ibgp
      dimension kap(30),nmax(30)
c    for photoelectron potential and wavefunction will be complex
      complex*16 eph,dg,ag,dp,ap,dv,av,eg,ceg,ep,cep,vxc(nrptx)
      complex*16 ps(nrptx),qs(nrptx),aps(10),aqs(10),vm(nrptx)
      common/comdic/cl,dz,dg(nrptx),ag(10),dp(nrptx),ap(10),
     1dv(nrptx),av(10),eg(nrptx),ceg(10),ep(nrptx),cep(10)
      common/itescf/testy,rap(2),teste,nz,norb,norbsc
      common/messag/dlabpr,numerr
      character*8 dlabpr
      common/snoyac/dvn(nrptx),anoy(10),nuc
      common/tabtec/hx,dr(nrptx),test1,test2,ndor,np,nes,method,idim
 
      cl=1.370373d+02
c     speed of light in atomic units
      dz = nz
c     make r-mesh and calculate nuclear potential
c     hx exponential step
c     dr1 first tabulation point multiplied by nz
      dr1= nz*exp(-8.8)
      call nucdec (anoy,dr,dvn,dz,hx,nuc,idim,10,dr1)
c     notice that here nuc=1, 
c     unless you specify nuclear mass and thickness in nucdec.f


      a=(dz/cl)**2
      if (nuc.gt.1) a=0.0d 00
      do 11 j=1,norb
         b=kap(j)*kap(j)-a
         if (j.eq.norb) b=b+(kap(j)+1)*ic3
         fl(j)= sqrt(b)
 11      fix(j) = dr(1)**(fl(j)-abs(kap(j)))
c     if irregular solution
      if (irr.gt.0) then
         fl(norb) = -fl(norb)
         fix(norb) = 1.0/fix(norb)
      endif

c     use lda potential to calculate initial w.f.
      do 21 i=1,jri-1
 21   dv(i)= vxc(i)/cl
      do  i=jri,idim
        dv(i)= vxc(jri+1)/cl
      enddo
      if (numerr.ne.0) return
      do 51 i=1,idim
         eg(i)=0.0d 00
 51      ep(i)=0.0d 00
      do 61 i=1,ibgp
         ceg(i)=0.0d 00
 61      cep(i)=0.0d 00
      call potdvp
      av(2)=av(2)+(vxc(nuc)-dvn(nuc))/cl

c     resolution of the dirac equation to get initial orbital
      if (irr.lt.0) then
         if (a .gt. 0.0d0) then 
            aps(1) = 1.0
            if (kap(norb) .lt. 0) then
               aqs(1)=aps(1)*dz/(cl*(kap(norb)-fl(norb)))
            else
               aqs(1)=aps(1)*cl*(kap(norb)+fl(norb))/dz
            endif
         else
            if (kap(norb).lt.0)then
               aps(1)=1.0d 00
               aqs(1)=0.0d 00
            else
               aps(1)=0.0d 00
               aqs(1)=1.0d 00
            endif
         endif
      endif

 211  np=1+(8.8 + log(10.0))/hx
c     exp(-8.8+(np-1)*hx) = 10.0 bohrs - max distance
      if (idim .lt. np) np=idim
      if (nmax(norb) .gt. np) nmax(norb)=np
         
      if (irr.lt.0) then
         call solout( eph, fl(norb), aps(1), aqs(1), kap(norb), rmt,
     1              jri, nmax(norb), ic3, vm, iwkb)
      else
         call solin( eph, fl(norb), aps(1), aqs(1), kap(norb), rmt,
     1              jri, nmax(norb), ic3, vm, iwkb)
      endif
         
      do 261 i=1,10
         aps(i)=ag(i)
 261     aqs(i)=ap(i)
      do 271 i=1,idim
         ps(i)=dg(i)
 271     qs(i)=dp(i)
      return
      end
      subroutine yzkrdc (i,k,flps,ps,qs,aps,aqs,p2, norb)
c       * calculate  function yk *
c yk = r * integral of f(s)*uk(r,s)
c uk(r,s) = rinf**k/rsup**(k+1)   rinf=min(r,s)   rsup=max(r,s)
c j=norb for photoelectron
c f(s)=cg(s,i)*cg(s,j)+cp(s,i)*cp(s,j)
c f(s) is constructed by the calling programm  if i < or =0
c in the last case a function f (lies in the block dg) is supposedly
c tabulated untill point dr(j), and its' devlopment coefficients
c at the origin are in ag and the power in r of the first term is k+2

c the output functions yk and zk are in the blocks dp and dg.
c at the origin  yk = cte * r**(k+1) - developement limit,
c cte lies in ap(1) and development coefficients in ag.
c        this programm uses aprdec and yzktec
 
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
      complex*16 aprdec,p2, dyzk
c     complex*16 a1,a2,b1,b2,coni
c     complex*16 xck, temp, ck, phx
      parameter (coni=(0.d0,1.d0))
      complex*16 ps(nrptx),qs(nrptx),aps(10),aqs(10)
      common/dff/cg(nrptx,30), cp(nrptx,30), bg(10,30), bp(10,30),
     1             fl(30), fix(30), ibgp
      complex*16 dg,ag,dp,ap,bidcom, chg(10)
      common/comdic/cl,dz,dg(nrptx),ag(10),dp(nrptx),ap(10),
     1   bidcom(3*nrptx+30)
      common/ratom1/xnel(30),en(30),scc(30),scw(30),sce(30),
     1   nq(30),kap(30),nmax(30)
      common/tabtec/hx,dr(nrptx),test1,test2,ndor,np,nes,method,idim
      dimension bgi(10),bpi(10)
c#mn
       external aprdec
 
c     construction of the function f
      do  5 l= 1,ibgp
        bgi(l) = bg(l,i)
  5     bpi(l) = bp(l,i)
      id=min(nmax(i),np)
      ap(1)=fl(i)+flps
      do 11 l=1,id
 11   dg(l)=cg(l,i)*ps(l)+cp(l,i)*qs(l)
      do 12 l = id+1,idim
 12    dg(l) = 0.0d0
      do 21 l=1,ndor
 21   ag(l) = aprdec(aps,bgi,l) + aprdec(aqs,bpi,l)

      dyzk = 0
c     if (id .ge. nmax(norb)) then
c        id = nmax(norb)-1
c        ck0 = log(cg(id,i)/cg(id+1,i))  / (dr(id+1)-dr(id))
c        ck = sqrt(2*p2)
c        xck = ck/cl
c        xck = -xck/(1+sqrt(1+xck**2))
c        temp = -ps(id+1) / qs(id+1) *xck
c        xx = dble (temp)
c        yy = dimag(temp)
c        if (xx .ne. 0)  then
c            alph = (1 - xx**2 - yy**2)
c            alph = sqrt(alph**2 + 4*xx**2) - alph
c            alph = alph / (2 * xx)
c            alph = atan (alph)
c        else
c            alph = 0
c        endif
c        beta = (xx**2 + (yy+1)**2) / (xx**2 + (yy-1)**2)
c        beta = log(beta) / 4

c        phx = dcmplx (alph, beta)
c        a1 =   ps(id+1) / sin(phx)
c        a2 = - qs(id+1) / cos(phx)
c        xck=ck*dr(id+1)
c        phx = phx -xck
c        a1 = a1*cg(id+1,i)/2/coni
c        a2 = a2*cp(id+1,i)/2
c        b1=exp(coni*phx) * (a1 - a2)
c        b2=exp(-coni*phx) * (-a1 - a2)
c        xck = (ck0 - coni*ck)*dr(id+1)
c        n = k +1
c        dyzk = dyzk + b1*exp(-xck)/xck
c        dyzk = dyzk + b1*expint(n,xck)
c        xck = (ck0 + coni*ck)*dr(id+1)
c        dyzk = dyzk + b2*exp(-xck)/xck
c        dyzk = dyzk + b2*expint(n,xck)
c        dyzk = dyzk*dr(id+1)
c     endif

      call yzktec (dg,ag,dp,chg,dr,ap(1),hx,k,ndor,id,idim, dyzk)
      return
      end

c     complex*16 function expint(n,x)
c     implicit double precision (a-h,o-z)
c     integer n, maxit
c     complex*16 x, b, c, d, h, del, fact, zero
c     parameter (zero=(0.d0,0.d0))
c     parameter (maxit=100, eps=1.d-7, fpmin=1.d-30, euler=.5772156649)

c     nm1 = n - 1
c     if (n.lt.0 .or. (dble(x).lt.0.d0 .and. dimag(x).eq.0.d0) .or.
c    1     (x.eq.zero .and. (n.eq.0.or.n.eq.1))) then
c        call par_stop('Bad arguments in expint')
c     elseif (n.eq.0) then
c        expint = exp(-x) / x
c     elseif (x.eq.0) then
c        expint = 1.d0 /nm1
c     elseif (dble(x).gt.1) then
c        b = x + n
c        c = 1/fpmin
c        d = 1/b
c        h = d
c        do 10 i=1,maxit
c           a = -i*(nm1+i)
c           b = b + 2
c           d = 1 / (a*d+b)
c           c = b + a/c
c           del = c*d
c           h = h*del
c           if (abs(del-1) .lt. eps) then
c              expint = h * exp(-x)
c              return
c           endif
c 10     continue
c        call par_stop(' continued fraction failed in expint')
c     else
c        if (nm1.ne.0) then
c           expint = 1/nm1
c        else
c           expint = -log(x) - euler
c        endif
c        fact = 1
c        do 30 i=1,maxit
c           fact = - fact *x / i
c           if (i.ne.nm1) then
c              del = - fact / (i-nm1)
c           else
c              psi = - euler
c              do 20 ii=1,nm1
c                 psi = psi + 1.d0 / ii
c 20           continue
c              del = fact*(-log(x)+psi)
c           endif
c           expint = expint + del
c           if (abs(del).lt.abs(expint)*eps) return
c 30     continue
c        call par_stop('series failed in expint')
c     endif
c     return
c     end
      subroutine yzktec (f,af,g,ag,dr,ap,h,k,nd,np,idim, dyzk)
c calculation of yk(r)=zk(r)+ r**(k+1) * integral from r to 
c   infinity of  f(u) * u**(-k-1)
c zk(r) = r**(-k) * integral from 0 to r of f(u) * u**k

c at the origin f(r)=sum from i=1 to nd of af(i)*r**(ap+i-1)
c dr tabulation points   h exponential step
c np number of tabulation points for f
c idim dimension of the blocks f,g and dr

c at the origin yk=cte*r**(k+1)-developement limit
c the constant for yk lies in ap
c output functions yk and zk lie in f and g, and their
c development coefficients at the origin in af and ag.

c integration from point to point by a 4 points method.
c integral from r to r+h = h*(-f(r-h)+13*f(r)+13*f(r+h)-f(r+h+h))/24

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
      complex*16 f,af,g,ag,ap, dyzk
      dimension f(nrptx),af(10),g(nrptx),ag(10),dr(nrptx)
 
c    initialisation and development coefficients of yk
      np= min(np,idim-1)
      f(np+1)=0.0d0
      b = dble(ap)
      ap=0.0d 00
      g(1)=0.0d 00
      do 15 i=1,nd
         b=b+1.0d 00
         ag(i)=af(i)/(b+k)
         if (af(i).ne.0.0d 00) then
            c=dr(1)**b
            g(1)=g(1)+ag(i)*c
c         for irregular solution b-k-1 can become zero
            if (abs(b-k-1) .le. 0.00001) then
               af(i) = 0.0
               b = b - 1.0d0
            else
               af(i)=(k+k+1)*ag(i)/(b-k-1)
            endif
            ap=ap+af(i)*c
         endif
 15   continue
      do 21 i=1,np
 21   f(i)=f(i)*dr(i)

c     calcualation of zk
      hk=h*k
      e = exp(-h)
      ehk = e**k 

      if (k.ne.0)then
       b1 = (ehk-1.0d0 +hk) / (hk*k)
      else
       b1=h/2.0
      endif

      b0 = h-(1.0+hk)*b1
      do 51 i=1,np
 51      g(i+1)=g(i)*ehk+b0*f(i)+f(i+1)*b1
 
c     calculation of yk
      f(np+1)=g(np+1) + dyzk
      ehk=ehk*e
      i=k+k+1
      hk=hk+h
      b1 = i*(ehk-1.0d0 +hk) / (hk*(k+1))
      b0 = i*h-(1.0+hk)*b1
      do 75  i=np,1,-1
 75      f(i) = f(i+1)*ehk+b0*g(i+1)+b1*g(i)

      ap=(ap+f(1))/(dr(1)**(k+1))
      return
      end
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
      
      subroutine cgetrf( M, N, A, LDA, IPIV, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX            A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  CGETRF computes an LU factorization of a general M-by-N matrix A
!  using partial pivoting with row interchanges.
!
!  The factorization has the form
!     A = P * L * U
!  where P is a permutation matrix, L is lower triangular with unit
!  diagonal elements (lower trapezoidal if m > n), and U is upper
!  triangular (upper trapezoidal if m < n).
!
!  This is the right-looking Level 3 BLAS version of the algorithm.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) COMPLEX array, dimension (LDA,N)
!          On entry, the M-by-N matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, U(i,i) is exactly zero. The
!                factorization has been completed, but the factor U
!                is exactly singular, and division by zero will occur
!                if it is used to solve a system of equations.
!
!  ==================================================================
!
!     .. Parameters ..
      COMPLEX            ONE
      PARAMETER          ( ONE = ( 1.0E+0, 0.0E+0 ) )
!     ..
!     .. Local Scalars ..
      INTEGER            I, IINFO, J, JB, NB
!     ..
!     .. External Subroutines ..
      EXTERNAL           CGEMM, CGETF2, CLASWP, CTRSM, XERBLA
!     ..
!     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CGETRF', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 )                                          &
     &   RETURN
!
!     Determine the block size for this environment.
!
      NB = ILAENV( 1, 'CGETRF', ' ', M, N, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) THEN
!
!        Use unblocked code.
!
         CALL CGETF2( M, N, A, LDA, IPIV, INFO )
      ELSE
!
!        Use blocked code.
!
         DO 20 J = 1, MIN( M, N ), NB
            JB = MIN( MIN( M, N )-J+1, NB )
!
!           Factor diagonal and subdiagonal blocks and test for exact
!           singularity.
!
            CALL CGETF2( M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO )
!
!           Adjust INFO and the pivot indices.
!
            IF( INFO.EQ.0 .AND. IINFO.GT.0 )                            &
     &         INFO = IINFO + J - 1
            DO 10 I = J, MIN( M, J+JB-1 )
               IPIV( I ) = J - 1 + IPIV( I )
   10       CONTINUE
!
!           Apply interchanges to columns 1:J-1.
!
            CALL CLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )
!
            IF( J+JB.LE.N ) THEN
!
!              Apply interchanges to columns J+JB:N.
!
               CALL CLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1,     &
     &                      IPIV, 1 )
!
!              Compute block row of U.
!
             CALL CTRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB,   &
     &                   N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ),   &
     &                     LDA )
               IF( J+JB.LE.M ) THEN
!
!                 Update trailing submatrix.
!
                CALL CGEMM( 'No transpose', 'No transpose', M-J-JB+1,   &
     &                       N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA,     &
     &                       A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ),   &
     &                       LDA )
               END IF
            END IF
   20    CONTINUE
      END IF
      RETURN
!
!     End of CGETRF
!
      END

      SUBROUTINE CGETF2( M, N, A, LDA, IPIV, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX            A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  CGETF2 computes an LU factorization of a general m-by-n matrix A
!  using partial pivoting with row interchanges.
!
!  The factorization has the form
!     A = P * L * U
!  where P is a permutation matrix, L is lower triangular with unit
!  diagonal elements (lower trapezoidal if m > n), and U is upper
!  triangular (upper trapezoidal if m < n).
!
!  This is the right-looking Level 2 BLAS version of the algorithm.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) COMPLEX array, dimension (LDA,N)
!          On entry, the m by n matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value
!          > 0: if INFO = k, U(k,k) is exactly zero. The factorization
!               has been completed, but the factor U is exactly
!               singular, and division by zero will occur if it is
!               used to solve a system of equations.
!
!  ==================================================================
!
!     .. Parameters ..
      COMPLEX            ONE, ZERO
      PARAMETER          ( ONE = ( 1.0E+0, 0.0E+0 ),                    &
     &                   ZERO = ( 0.0E+0, 0.0E+0 ) )
!     ..
!     .. Local Scalars ..
      INTEGER            J, JP
!     ..
!     .. External Functions ..
      INTEGER            ICAMAX
      EXTERNAL           ICAMAX
!     ..
!     .. External Subroutines ..
      EXTERNAL           CGERU, CSCAL, CSWAP, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CGETF2', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 )                                          &
     &   RETURN
!
      DO 10 J = 1, MIN( M, N )
!
!        Find pivot and test for singularity.
!
         JP = J - 1 + ICAMAX( M-J+1, A( J, J ), 1 )
         IPIV( J ) = JP
         IF( A( JP, J ).NE.ZERO ) THEN
!
!           Apply the interchange to columns 1:N.
!
            IF( JP.NE.J )                                               &
     &         CALL CSWAP( N, A( J, 1 ), LDA, A( JP, 1 ), LDA )
!
!           Compute elements J+1:M of J-th column.
!
            IF( J.LT.M )                                                &
     &         CALL CSCAL( M-J, ONE / A( J, J ), A( J+1, J ), 1 )
!
         ELSE IF( INFO.EQ.0 ) THEN
!
            INFO = J
         END IF
!
         IF( J.LT.MIN( M, N ) ) THEN
!
!           Update trailing submatrix.
!
            CALL CGERU( M-J, N-J, -ONE, A( J+1, J ), 1, A( J, J+1 ),    &
     &                  LDA, A( J+1, J+1 ), LDA )
         END IF
   10 CONTINUE
      RETURN
!
!     End of CGETF2
!
      END
      SUBROUTINE cgetrs( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, LDA, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX            A( LDA, * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  CGETRS solves a system of linear equations
!     A * X = B,  A**T * X = B,  or  A**H * X = B
!  with a general N-by-N matrix A using the LU factorization computed
!  by CGETRF.
!
!  Arguments
!  =========
!
!  TRANS   (input) CHARACTER*1
!          Specifies the form of the system of equations:
!          = 'N':  A * X = B     (No transpose)
!          = 'T':  A**T * X = B  (Transpose)
!          = 'C':  A**H * X = B  (Conjugate transpose)
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  A       (input) COMPLEX array, dimension (LDA,N)
!          The factors L and U from the factorization A = P*L*U
!          as computed by CGETRF.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  IPIV    (input) INTEGER array, dimension (N)
!          The pivot indices from CGETRF; for 1<=i<=N, row i of the
!          matrix was interchanged with row IPIV(i).
!
!  B       (input/output) COMPLEX array, dimension (LDB,NRHS)
!          On entry, the right hand side matrix B.
!          On exit, the solution matrix X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  ==================================================================
!
!     .. Parameters ..
      COMPLEX            ONE
      PARAMETER          ( ONE = ( 1.0E+0, 0.0E+0 ) )
!     ..
!     .. Local Scalars ..
      LOGICAL            NOTRAN
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           CLASWP, CTRSM, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.        &
     &    LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CGETRS', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 .OR. NRHS.EQ.0 )                                       &
     &   RETURN
!
      IF( NOTRAN ) THEN
!
!        Solve A * X = B.
!
!        Apply row interchanges to the right hand sides.
!
         CALL CLASWP( NRHS, B, LDB, 1, N, IPIV, 1 )
!
!        Solve L*X = B, overwriting B with X.
!
         CALL CTRSM( 'Left', 'Lower', 'No transpose', 'Unit', N, NRHS,  &
     &               ONE, A, LDA, B, LDB )
!
!        Solve U*X = B, overwriting B with X.
!
         CALL CTRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N,    &
     &               NRHS, ONE, A, LDA, B, LDB )
      ELSE
!
!        Solve A**T * X = B  or A**H * X = B.
!
!        Solve U'*X = B, overwriting B with X.
!
         CALL CTRSM( 'Left', 'Upper', TRANS, 'Non-unit', N, NRHS, ONE,  &
     &               A, LDA, B, LDB )
!
!        Solve L'*X = B, overwriting B with X.
!
         CALL CTRSM( 'Left', 'Lower', TRANS, 'Unit', N, NRHS, ONE, A,   &
     &               LDA, B, LDB )
!
!        Apply row interchanges to the solution vectors.
!
         CALL CLASWP( NRHS, B, LDB, 1, N, IPIV, -1 )
      END IF
!
      RETURN
!
!     End of CGETRS
!
      END
      SUBROUTINE XERBLA( SRNAME, INFO )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER*6        SRNAME
      INTEGER            INFO
!     ..
!
!  Purpose
!  =======
!
!  XERBLA  is an error handler for the LAPACK routines.
!  It is called by an LAPACK routine if an input parameter has an
!  invalid value.  A message is printed and execution stops.
!
!  Installers may consider modifying the STOP statement in order to
!  call system-specific exception-handling facilities.
!
!  Arguments
!  =========
!
!  SRNAME  (input) CHARACTER*6
!          The name of the routine which called XERBLA.
!
!  INFO    (input) INTEGER
!          The position of the invalid parameter in the parameter list
!          of the calling routine.
!
! ==================================================================
!
!     .. Executable Statements ..
!
      WRITE( *, FMT = 9999 )SRNAME, INFO
!
      STOP
!
 9999 FORMAT( ' ** On entry to ', A6, ' parameter number ',I2,' had ',  &
     &      'an illegal value' )
!
!     End of XERBLA
!
      END
      subroutine  cswap (n,cx,incx,cy,incy)
!
!     interchanges two vectors.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      complex cx(*),cy(*),ctemp
      integer i,incx,incy,ix,iy,n
!
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
!
!       code for unequal increments or equal increments not equal
!         to 1
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        ctemp = cx(ix)
        cx(ix) = cy(iy)
        cy(iy) = ctemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!
!       code for both increments equal to 1
   20 do 30 i = 1,n
        ctemp = cx(i)
        cx(i) = cy(i)
        cy(i) = ctemp
   30 continue
      return
      end
      subroutine  cscal(n,ca,cx,incx)
!
!     scales a vector by a constant.
!     jack dongarra, linpack,  3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      complex ca,cx(*)
      integer i,incx,n,nincx
!
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
!
!        code for increment not equal to 1
!
      nincx = n*incx
      do 10 i = 1,nincx,incx
        cx(i) = ca*cx(i)
   10 continue
      return
!
!        code for increment equal to 1
!
   20 do 30 i = 1,n
        cx(i) = ca*cx(i)
   30 continue
      return
      end
      SUBROUTINE CGERU ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
!     .. Scalar Arguments ..
      COMPLEX            ALPHA
      INTEGER            INCX, INCY, LDA, M, N
!     .. Array Arguments ..
      COMPLEX            A( LDA, * ), X( * ), Y( * )
!     ..
!
!  Purpose
!  =======
!
!  CGERU  performs the rank 1 operation
!
!     A := alpha*x*y' + A,
!
!  where alpha is a scalar, x is an m element vector, y is an n
!  element vector and A is an m by n matrix.
!
!  Parameters
!  ==========
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of the
!           matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - COMPLEX         .
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  X      - COMPLEX          array of dimension at least
!           ( 1 + ( m - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the m
!           element vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  Y      - COMPLEX          array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCY ) ).
!           Before entry, the incremented array Y must contain the n
!           element vector y.
!           Unchanged on exit.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!  A      - COMPLEX          array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients. On exit, A is
!           overwritten by the updated matrix.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as
!           declared in the calling (sub) program. LDA must be
!           at least max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
!     .. Parameters ..
      COMPLEX            ZERO
      PARAMETER        ( ZERO = ( 0.0E+0, 0.0E+0 ) )
!     .. Local Scalars ..
      COMPLEX            TEMP
      INTEGER            I, INFO, IX, J, JY, KX
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF     ( M.LT.0 )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 7
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'CGERU ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )               &
     &   RETURN
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
      IF( INCY.GT.0 )THEN
         JY = 1
      ELSE
         JY = 1 - ( N - 1 )*INCY
      END IF
      IF( INCX.EQ.1 )THEN
         DO 20, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               DO 10, I = 1, M
                  A( I, J ) = A( I, J ) + X( I )*TEMP
   10          CONTINUE
            END IF
            JY = JY + INCY
   20    CONTINUE
      ELSE
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( M - 1 )*INCX
         END IF
         DO 40, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               IX   = KX
               DO 30, I = 1, M
                  A( I, J ) = A( I, J ) + X( IX )*TEMP
                  IX        = IX        + INCX
   30          CONTINUE
            END IF
            JY = JY + INCY
   40    CONTINUE
      END IF
!
      RETURN
!
!     End of CGERU .
!
      END
      SUBROUTINE CLASWP( N, A, LDA, K1, K2, IPIV, INCX )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      INTEGER            INCX, K1, K2, LDA, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX            A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  CLASWP performs a series of row interchanges on the matrix A.
!  One row interchange is initiated for each of rows K1 through
!  K2 of A.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.
!
!  A       (input/output) COMPLEX array, dimension (LDA,N)
!          On entry, the matrix of column dimension N to which the row
!          interchanges will be applied.
!          On exit, the permuted matrix.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.
!
!  K1      (input) INTEGER
!          The first element of IPIV for which a row interchange will
!          be done.
!
!  K2      (input) INTEGER
!          The last element of IPIV for which a row interchange will
!          be done.
!
!  IPIV    (input) INTEGER array, dimension (M*abs(INCX))
!          The vector of pivot indices.  Only the elements in
!          positions K1 through K2 of IPIV are accessed.
!          IPIV(K) = L implies rows K and L are to be interchanged.
!
!  INCX    (input) INTEGER
!          The increment between successive values of IPIV.  If IPIV
!          is negative, the pivots are applied in reverse order.
!
! ==================================================================
!
!     .. Local Scalars ..
      INTEGER            I, IP, IX
!     ..
!     .. External Subroutines ..
      EXTERNAL           CSWAP
!     ..
!     .. Executable Statements ..
!
!     Interchange row I with row IPIV(I) for each of rows K1
!     through K2.
!
      IF( INCX.EQ.0 )                                                   &
     &   RETURN
      IF( INCX.GT.0 ) THEN
         IX = K1
      ELSE
         IX = 1 + ( 1-K2 )*INCX
      END IF
      IF( INCX.EQ.1 ) THEN
         DO 10 I = K1, K2
            IP = IPIV( I )
            IF( IP.NE.I )                                               &
     &         CALL CSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
   10    CONTINUE
      ELSE IF( INCX.GT.1 ) THEN
         DO 20 I = K1, K2
            IP = IPIV( IX )
            IF( IP.NE.I )                                               &
     &         CALL CSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
            IX = IX + INCX
   20    CONTINUE
      ELSE IF( INCX.LT.0 ) THEN
         DO 30 I = K2, K1, -1
            IP = IPIV( IX )
            IF( IP.NE.I )                                               &
     &         CALL CSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
            IX = IX + INCX
   30    CONTINUE
      END IF
!
      RETURN
!
!     End of CLASWP
!
      END
      SUBROUTINE CTRSM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA,A,LDA,   &
     &                   B, LDB )
!     .. Scalar Arguments ..
      CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, LDA, LDB
      COMPLEX            ALPHA
!     .. Array Arguments ..
      COMPLEX            A( LDA, * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  CTRSM  solves one of the matrix equations
!
!     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
!
!  where alpha is a scalar, X and B are m by n matrices, A is a unit,
!  or non-unit,  upper or lower triangular matrix  and  op( A )  is
!  one of
!
!     op( A ) = A   or   op( A ) = A'   or   op( A ) = conjg( A' ).
!
!  The matrix X is overwritten on B.
!
!  Parameters
!  ==========
!
!  SIDE   - CHARACTER*1.
!           On entry, SIDE specifies whether op( A ) appears on the
!           left or right of X as follows:
!
!              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
!
!              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
!
!           Unchanged on exit.
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the matrix A is an upper
!           or lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be used
!           in the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n'   op( A ) = A.
!
!              TRANSA = 'T' or 't'   op( A ) = A'.
!
!              TRANSA = 'C' or 'c'   op( A ) = conjg( A' ).
!
!           Unchanged on exit.
!
!  DIAG   - CHARACTER*1.
!           On entry, DIAG specifies whether or not A is unit
!           triangular as follows:
!
!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of B. M must
!           be at least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of B.  N
!           must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - COMPLEX         .
!           On entry,  ALPHA specifies the scalar  alpha. When alpha
!           is zero then  A is not referenced and  B need not be set
!           before entry.
!           Unchanged on exit.
!

!  A - COMPLEX array of DIMENSION ( LDA, k ), where k is m
!           when SIDE = 'L' or 'l' and is n when SIDE = 'R' or 'r'.
!           Before entry with UPLO = 'U' or 'u', the leading k by k
!           upper triangular part of the array A must contain the
!           upper triangular matrix and the strictly lower triangular
!           part of A is not referenced.
!           Before entry  with  UPLO = 'L' or 'l',  the  leading  k
!           by k lower triangular part of the array  A must contain
!           the lower triangular matrix  and the strictly upper
!           triangular part of A is not referenced.
!           Note that when  DIAG = 'U' or 'u',  the diagonal elements
!           of A  are not referenced either,  but are assumed to be
!           unity.  Unchanged on exit.

!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as
!           declared in the calling (sub) program.  When
!           SIDE = 'L' or 'l'  then LDA  must be at least
!           max( 1, m ),  when  SIDE = 'R' or 'r'
!           then LDA must be at least max( 1, n ).
!           Unchanged on exit.
!
!  B      - COMPLEX          array of DIMENSION ( LDB, n ).
!           Before entry,  the leading  m by n part of the array
!           B must contain  the  right-hand  side  matrix  B,  and
!           on exit  is overwritten by the solution matrix  X.
!
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as
!           declared in  the  calling  (sub)  program.   LDB  must
!           be  at  least max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 3 Blas routine.
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC          CONJG, MAX
!     .. Local Scalars ..
      LOGICAL            LSIDE, NOCONJ, NOUNIT, UPPER
      INTEGER            I, INFO, J, K, NROWA
      COMPLEX            TEMP
!     .. Parameters ..
      COMPLEX            ONE
      PARAMETER        ( ONE  = ( 1.0E+0, 0.0E+0 ) )
      COMPLEX            ZERO
      PARAMETER        ( ZERO = ( 0.0E+0, 0.0E+0 ) )
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      LSIDE  = LSAME( SIDE  , 'L' )
      IF( LSIDE )THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      NOCONJ = LSAME( TRANSA, 'T' )
      NOUNIT = LSAME( DIAG  , 'N' )
      UPPER  = LSAME( UPLO  , 'U' )
!
      INFO   = 0
      IF(      ( .NOT.LSIDE                ).AND.                       &
     &         ( .NOT.LSAME( SIDE  , 'R' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.UPPER                ).AND.                       &
     &         ( .NOT.LSAME( UPLO  , 'L' ) )      )THEN
         INFO = 2
      ELSE IF( ( .NOT.LSAME( TRANSA, 'N' ) ).AND.                       &
     &         ( .NOT.LSAME( TRANSA, 'T' ) ).AND.                       &
     &         ( .NOT.LSAME( TRANSA, 'C' ) )      )THEN
         INFO = 3
      ELSE IF( ( .NOT.LSAME( DIAG  , 'U' ) ).AND.                       &
     &         ( .NOT.LSAME( DIAG  , 'N' ) )      )THEN
         INFO = 4
      ELSE IF( M  .LT.0               )THEN
         INFO = 5
      ELSE IF( N  .LT.0               )THEN
         INFO = 6
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 9
      ELSE IF( LDB.LT.MAX( 1, M     ) )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'CTRSM ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( N.EQ.0 )                                                      &
     &   RETURN
!
!     And when  alpha.eq.zero.
!
      IF( ALPHA.EQ.ZERO )THEN
         DO 20, J = 1, N
            DO 10, I = 1, M
               B( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
         RETURN
      END IF
!
!     Start the operations.
!
      IF( LSIDE )THEN
         IF( LSAME( TRANSA, 'N' ) )THEN
!
!           Form  B := alpha*inv( A )*B.
!
            IF( UPPER )THEN
               DO 60, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 30, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   30                CONTINUE
                  END IF
                  DO 50, K = M, 1, -1
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT )                                    &
     &                     B( K, J ) = B( K, J )/A( K, K )
                        DO 40, I = 1, K - 1
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   40                   CONTINUE
                     END IF
   50             CONTINUE
   60          CONTINUE
            ELSE
               DO 100, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 70, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   70                CONTINUE
                  END IF
                  DO 90 K = 1, M
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT )                                    &
     &                     B( K, J ) = B( K, J )/A( K, K )
                        DO 80, I = K + 1, M
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   80                   CONTINUE
                     END IF
   90             CONTINUE
  100          CONTINUE
            END IF
         ELSE
!
!           Form  B := alpha*inv( A' )*B
!           or    B := alpha*inv( conjg( A' ) )*B.
!
            IF( UPPER )THEN
               DO 140, J = 1, N
                  DO 130, I = 1, M
                     TEMP = ALPHA*B( I, J )
                     IF( NOCONJ )THEN
                        DO 110, K = 1, I - 1
                           TEMP = TEMP - A( K, I )*B( K, J )
  110                   CONTINUE
                        IF( NOUNIT )                                    &
     &                     TEMP = TEMP/A( I, I )
                     ELSE
                        DO 120, K = 1, I - 1
                           TEMP = TEMP - CONJG( A( K, I ) )*B( K, J )
  120                   CONTINUE
                        IF( NOUNIT )                                    &
     &                     TEMP = TEMP/CONJG( A( I, I ) )
                     END IF
                     B( I, J ) = TEMP
  130             CONTINUE
  140          CONTINUE
            ELSE
               DO 180, J = 1, N
                  DO 170, I = M, 1, -1
                     TEMP = ALPHA*B( I, J )
                     IF( NOCONJ )THEN
                        DO 150, K = I + 1, M
                           TEMP = TEMP - A( K, I )*B( K, J )
  150                   CONTINUE
                        IF( NOUNIT )                                    &
     &                     TEMP = TEMP/A( I, I )
                     ELSE
                        DO 160, K = I + 1, M
                           TEMP = TEMP - CONJG( A( K, I ) )*B( K, J )
  160                   CONTINUE
                        IF( NOUNIT )                                    &
     &                     TEMP = TEMP/CONJG( A( I, I ) )
                     END IF
                     B( I, J ) = TEMP
  170             CONTINUE
  180          CONTINUE
            END IF
         END IF
      ELSE
         IF( LSAME( TRANSA, 'N' ) )THEN
!
!           Form  B := alpha*B*inv( A ).
!
            IF( UPPER )THEN
               DO 230, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 190, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  190                CONTINUE
                  END IF
                  DO 210, K = 1, J - 1
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 200, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  200                   CONTINUE
                     END IF
  210             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 220, I = 1, M
                        B( I, J ) = TEMP*B( I, J )
  220                CONTINUE
                  END IF
  230          CONTINUE
            ELSE
               DO 280, J = N, 1, -1
                  IF( ALPHA.NE.ONE )THEN
                     DO 240, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  240                CONTINUE
                  END IF
                  DO 260, K = J + 1, N
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 250, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  250                   CONTINUE
                     END IF
  260             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 270, I = 1, M
                       B( I, J ) = TEMP*B( I, J )
  270                CONTINUE
                  END IF
  280          CONTINUE
            END IF
         ELSE
!
!           Form  B := alpha*B*inv( A' )
!           or    B := alpha*B*inv( conjg( A' ) ).
!
            IF( UPPER )THEN
               DO 330, K = N, 1, -1
                  IF( NOUNIT )THEN
                     IF( NOCONJ )THEN
                        TEMP = ONE/A( K, K )
                     ELSE
                        TEMP = ONE/CONJG( A( K, K ) )
                     END IF
                     DO 290, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  290                CONTINUE
                  END IF
                  DO 310, J = 1, K - 1
                     IF( A( J, K ).NE.ZERO )THEN
                        IF( NOCONJ )THEN
                           TEMP = A( J, K )
                        ELSE
                           TEMP = CONJG( A( J, K ) )
                        END IF
                        DO 300, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  300                   CONTINUE
                     END IF
  310             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 320, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  320                CONTINUE
                  END IF
  330          CONTINUE
            ELSE
               DO 380, K = 1, N
                  IF( NOUNIT )THEN
                     IF( NOCONJ )THEN
                        TEMP = ONE/A( K, K )
                     ELSE
                        TEMP = ONE/CONJG( A( K, K ) )
                     END IF
                     DO 340, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  340                CONTINUE
                  END IF
                  DO 360, J = K + 1, N
                     IF( A( J, K ).NE.ZERO )THEN
                        IF( NOCONJ )THEN
                           TEMP = A( J, K )
                        ELSE
                           TEMP = CONJG( A( J, K ) )
                        END IF
                        DO 350, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  350                   CONTINUE
                     END IF
  360             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 370, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  370                CONTINUE
                  END IF
  380          CONTINUE
            END IF
         END IF
      END IF
!
      RETURN
!
!     End of CTRSM .
!
      END
      SUBROUTINE CGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA,B,LDB,  &
     &                   BETA, C, LDC )
!     .. Scalar Arguments ..
      CHARACTER*1        TRANSA, TRANSB
      INTEGER            M, N, K, LDA, LDB, LDC
      COMPLEX            ALPHA, BETA
!     .. Array Arguments ..
      COMPLEX            A( LDA, * ), B( LDB, * ), C( LDC, * )
!     ..
!
!  Purpose
!  =======
!
!  CGEMM  performs one of the matrix-matrix operations
!
!     C := alpha*op( A )*op( B ) + beta*C,
!
!  where  op( X ) is one of
!
!     op( X ) = X   or   op( X ) = X'   or   op( X ) = conjg( X' ),
!
!  alpha and beta are scalars, and A, B and C are matrices,
!  with op( A ) an m by k matrix,  op( B )  a  k by n matrix and
!  C an m by n matrix.
!
!  Parameters
!  ==========
!
!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be
!           used in the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n',  op( A ) = A.
!
!              TRANSA = 'T' or 't',  op( A ) = A'.
!
!              TRANSA = 'C' or 'c',  op( A ) = conjg( A' ).
!
!           Unchanged on exit.
!
!  TRANSB - CHARACTER*1.
!           On entry, TRANSB specifies the form of op( B ) to be
!           used in the matrix multiplication as follows:
!
!              TRANSB = 'N' or 'n',  op( B ) = B.
!
!              TRANSB = 'T' or 't',  op( B ) = B'.
!
!              TRANSB = 'C' or 'c',  op( B ) = conjg( B' ).
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry,  M  specifies  the number  of rows  of the
!           matrix op( A )  and of the  matrix  C.  M  must  be at
!           least  zero. Unchanged on exit.
!
!  N      - INTEGER.
!           On entry,  N  specifies the number  of columns of the
!           matrix op( B ) and the number of columns of the matrix C.
!           N must be at least zero.
!           Unchanged on exit.
!
!  K      - INTEGER.
!           On entry,  K  specifies  the number of columns of the
!           matrix op( A ) and the number of rows of the matrix
!           op( B ). K must be at least  zero.
!           Unchanged on exit.
!
!  ALPHA  - COMPLEX         .
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A   - COMPLEX          array of DIMENSION ( LDA, ka ), where ka is
!        k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
!        Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
!        part of the array  A  must contain the matrix  A,  otherwise
!        the leading  k by m  part of the array  A  must contain  the
!        matrix A.
!        Unchanged on exit.
!
!  LDA - INTEGER.
!        On entry, LDA specifies the first dimension of A as declared
!        in the calling (sub) program. When  TRANSA = 'N' or 'n' then
!        LDA must be at least  max( 1, m ), otherwise  LDA must be at
!        least  max( 1, k ).
!        Unchanged on exit.
!
!  B   - COMPLEX          array of DIMENSION ( LDB, kb ), where kb is
!        n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
!        Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
!        part of the array  B  must contain the matrix  B,  otherwise
!        the leading  n by k  part of the array  B  must contain  the
!        matrix B.
!        Unchanged on exit.
!
!  LDB - INTEGER.
!        On entry, LDB specifies the first dimension of B as declared
!        in the calling (sub) program. When  TRANSB = 'N' or 'n' then
!        LDB must be at least  max( 1, k ), otherwise  LDB must be at
!        least  max( 1, n ).
!        Unchanged on exit.
!
!  BETA   - COMPLEX         .
!        On entry,  BETA  specifies the scalar  beta.  When  BETA  is
!        supplied as zero then C need not be set on input.
!        Unchanged on exit.
!
!  C      - COMPLEX          array of DIMENSION ( LDC, n ).
!        Before entry, the leading  m by n  part of the array  C must
!        contain the matrix  C,  except when  beta  is zero, in which
!        case C need not be set on entry.
!        On exit, the array  C  is overwritten by the  m by n  matrix
!        ( alpha*op( A )*op( B ) + beta*C ).
!
!  LDC    - INTEGER.
!        On entry, LDC specifies the first dimension of C as declared
!        in  the  calling  (sub)  program.   LDC  must  be  at  least
!        max( 1, m ).
!        Unchanged on exit.
!
!
!  Level 3 Blas routine.
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC          CONJG, MAX
!     .. Local Scalars ..
      LOGICAL            CONJA, CONJB, NOTA, NOTB
      INTEGER            I, INFO, J, L, NCOLA, NROWA, NROWB
      COMPLEX            TEMP
!     .. Parameters ..
      COMPLEX            ONE
      PARAMETER        ( ONE  = ( 1.0E+0, 0.0E+0 ) )
      COMPLEX            ZERO
      PARAMETER        ( ZERO = ( 0.0E+0, 0.0E+0 ) )
!     ..
!     .. Executable Statements ..
!
!  Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
!  conjugated or transposed, set  CONJA and CONJB  as true if  A  and
!  B  respectively are to be  transposed but  not conjugated  and set
!  NROWA, NCOLA and  NROWB  as the number of rows and  columns  of  A
!  and the number of rows of  B  respectively.
!
      NOTA  = LSAME( TRANSA, 'N' )
      NOTB  = LSAME( TRANSB, 'N' )
      CONJA = LSAME( TRANSA, 'C' )
      CONJB = LSAME( TRANSB, 'C' )
      IF( NOTA )THEN
         NROWA = M
         NCOLA = K
      ELSE
         NROWA = K
         NCOLA = M
      END IF
      IF( NOTB )THEN
         NROWB = K
      ELSE
         NROWB = N
      END IF
!
!     Test the input parameters.
!
      INFO = 0
      IF(      ( .NOT.NOTA                 ).AND.                       &
     &         ( .NOT.CONJA                ).AND.                       &
     &         ( .NOT.LSAME( TRANSA, 'T' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.NOTB                 ).AND.                       &
     &         ( .NOT.CONJB                ).AND.                       &
     &         ( .NOT.LSAME( TRANSB, 'T' ) )      )THEN
         INFO = 2
      ELSE IF( M  .LT.0               )THEN
         INFO = 3
      ELSE IF( N  .LT.0               )THEN
         INFO = 4
      ELSE IF( K  .LT.0               )THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 8
      ELSE IF( LDB.LT.MAX( 1, NROWB ) )THEN
         INFO = 10
      ELSE IF( LDC.LT.MAX( 1, M     ) )THEN
         INFO = 13
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'CGEMM ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.                                  &
     &    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE )))   &
     &   RETURN
!
!     And when  alpha.eq.zero.
!
      IF( ALPHA.EQ.ZERO )THEN
         IF( BETA.EQ.ZERO )THEN
            DO 20, J = 1, N
               DO 10, I = 1, M
                  C( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               DO 30, I = 1, M
                  C( I, J ) = BETA*C( I, J )
   30          CONTINUE
   40       CONTINUE
         END IF
         RETURN
      END IF
!
!     Start the operations.
!
      IF( NOTB )THEN
         IF( NOTA )THEN
!
!           Form  C := alpha*A*B + beta*C.
!
            DO 90, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 50, I = 1, M
                     C( I, J ) = ZERO
   50             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 60, I = 1, M
                     C( I, J ) = BETA*C( I, J )
   60             CONTINUE
               END IF
               DO 80, L = 1, K
                  IF( B( L, J ).NE.ZERO )THEN
                     TEMP = ALPHA*B( L, J )
                     DO 70, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
   70                CONTINUE
                  END IF
   80          CONTINUE
   90       CONTINUE
         ELSE IF( CONJA )THEN
!
!           Form  C := alpha*conjg( A' )*B + beta*C.
!
            DO 120, J = 1, N
               DO 110, I = 1, M
                  TEMP = ZERO
                  DO 100, L = 1, K
                     TEMP = TEMP + CONJG( A( L, I ) )*B( L, J )
  100             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  110          CONTINUE
  120       CONTINUE
         ELSE
!
!           Form  C := alpha*A'*B + beta*C
!
            DO 150, J = 1, N
               DO 140, I = 1, M
                  TEMP = ZERO
                  DO 130, L = 1, K
                     TEMP = TEMP + A( L, I )*B( L, J )
  130             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  140          CONTINUE
  150       CONTINUE
         END IF
      ELSE IF( NOTA )THEN
         IF( CONJB )THEN
!
!           Form  C := alpha*A*conjg( B' ) + beta*C.
!
            DO 200, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 160, I = 1, M
                     C( I, J ) = ZERO
  160             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 170, I = 1, M
                     C( I, J ) = BETA*C( I, J )
  170             CONTINUE
               END IF
               DO 190, L = 1, K
                  IF( B( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*CONJG( B( J, L ) )
                     DO 180, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  180                CONTINUE
                  END IF
  190          CONTINUE
  200       CONTINUE
         ELSE
!
!           Form  C := alpha*A*B'          + beta*C
!
            DO 250, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 210, I = 1, M
                     C( I, J ) = ZERO
  210             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 220, I = 1, M
                     C( I, J ) = BETA*C( I, J )
  220             CONTINUE
               END IF
               DO 240, L = 1, K
                  IF( B( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*B( J, L )
                     DO 230, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  230                CONTINUE
                  END IF
  240          CONTINUE
  250       CONTINUE
         END IF
      ELSE IF( CONJA )THEN
         IF( CONJB )THEN
!
!           Form  C := alpha*conjg( A' )*conjg( B' ) + beta*C.
!
            DO 280, J = 1, N
               DO 270, I = 1, M
                  TEMP = ZERO
                  DO 260, L = 1, K
                  TEMP = TEMP + CONJG( A( L, I ) )*CONJG( B( J, L ) )
  260             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  270          CONTINUE
  280       CONTINUE
         ELSE
!
!           Form  C := alpha*conjg( A' )*B' + beta*C
!
            DO 310, J = 1, N
               DO 300, I = 1, M
                  TEMP = ZERO
                  DO 290, L = 1, K
                     TEMP = TEMP + CONJG( A( L, I ) )*B( J, L )
  290             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  300          CONTINUE
  310       CONTINUE
         END IF
      ELSE
         IF( CONJB )THEN
!
!           Form  C := alpha*A'*conjg( B' ) + beta*C
!
            DO 340, J = 1, N
               DO 330, I = 1, M
                  TEMP = ZERO
                  DO 320, L = 1, K
                     TEMP = TEMP + A( L, I )*CONJG( B( J, L ) )
  320             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  330          CONTINUE
  340       CONTINUE
         ELSE
!
!           Form  C := alpha*A'*B' + beta*C
!
            DO 370, J = 1, N
               DO 360, I = 1, M
                  TEMP = ZERO
                  DO 350, L = 1, K
                     TEMP = TEMP + A( L, I )*B( J, L )
  350             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  360          CONTINUE
  370       CONTINUE
         END IF
      END IF
!
      RETURN
!
!     End of CGEMM .
!
      END

      INTEGER          FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3,  &
     &                 N4 )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER*( * )    NAME, OPTS
      INTEGER            ISPEC, N1, N2, N3, N4
!     ..
!
!  Purpose
!  =======
!
!  ILAENV is called from the LAPACK routines to choose problem-
!  dependent parameters for the local environment.  See ISPEC for
!  a description of the parameters.
!
!  This version provides a set of parameters which should give good,
!  but not optimal, performance on many of the currently available
!  computers.  Users are encouraged to modify this subroutine to set
!  the tuning parameters for their particular machine using the option
!  and problem size information in the arguments.
!
!  This routine will not function correctly if it is converted to all
!  lower case.  Converting it to all upper case is allowed.
!
!  Arguments
!  =========
!
!  ISPEC   (input) INTEGER
!       Specifies the parameter to be returned as the value of
!       ILAENV.
!       = 1: the optimal blocksize; if this value is 1, an unblocked
!            algorithm will give the best performance.
!       = 2: the minimum block size for which the block routine
!            should be used; if the usable block size is less than
!            this value, an unblocked routine should be used.
!       = 3: the crossover point (in a block routine, for N less
!            than this value, an unblocked routine should be used)
!       = 4: the number of shifts, used in the nonsymmetric
!            eigenvalue routines
!       = 5: the minimum column dimension for blocking to be used;
!            rectangular blocks must have dimension at least k by m,
!            where k is given by ILAENV(2,...) and m by ILAENV(5,...)
!       = 6: the crossover point for the SVD (when reducing an m by n
!            matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
!            this value, a QR factorization is used first to reduce
!            the matrix to a triangular form.)
!       = 7: the number of processors
!       = 8: the crossover point for the multishift QR and QZ methods
!            for nonsymmetric eigenvalue problems.
!
!  NAME    (input) CHARACTER*(*)
!          The name of the calling subroutine, in either upper case or
!          lower case.
!
!  OPTS    (input) CHARACTER*(*)
!          The character options to the subroutine NAME, concatenated
!          into a single character string.  For example, UPLO = 'U',
!          TRANS = 'T', and DIAG = 'N' for a triangular routine would
!          be specified as OPTS = 'UTN'.
!
!  N1      (input) INTEGER
!  N2      (input) INTEGER
!  N3      (input) INTEGER
!  N4      (input) INTEGER
!          Problem dimensions for the subroutine NAME; these may not
!          all be required.
!
! (ILAENV) (output) INTEGER
!          >= 0: the value of the parameter specified by ISPEC
!          < 0:  if ILAENV = -k, the k-th argument had an illegal
!                value.
!
!  Further Details
!  ===============
!
!  The following conventions have been used when calling ILAENV from
!  the LAPACK routines:
!  1)  OPTS is a concatenation of all of the character options to
!      subroutine NAME, in the same order that they appear in the
!      argument list for NAME, even if they are not used in
!      determining the value of the parameter specified by ISPEC.
!  2)  The problem dimensions N1, N2, N3, N4 are specified in the
!      order that they appear in the argument list for NAME.  N1 is
!      used first, N2 second, and so on, and unused problem dimensions
!      are passed a value of -1.
!  3)  The parameter value returned by ILAENV is checked for validity
!      in the calling subroutine.  For example, ILAENV is used to
!      retrieve the optimal blocksize for STRTRI as follows:
!
!      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
!      IF( NB.LE.1 ) NB = MAX( 1, N )
!
!  ==================================================================
!
!     .. Local Scalars ..
      LOGICAL            CNAME, SNAME
      CHARACTER*1        C1
      CHARACTER*2        C2, C4
      CHARACTER*3        C3
      CHARACTER*6        SUBNAM
      INTEGER            I, IC, IZ, NB, NBMIN, NX
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          CHAR, ICHAR, INT, MIN, REAL
!     ..
!     .. Executable Statements ..
!
      GO TO ( 100, 100, 100, 400, 500, 600, 700, 800 ) ISPEC
!
!     Invalid value for ISPEC
!
      ILAENV = -1
      RETURN
!
  100 CONTINUE
!
!     Convert NAME to upper case if the first character is lower case.
!
      ILAENV = 1
      SUBNAM = NAME
      IC = ICHAR( SUBNAM( 1:1 ) )
      IZ = ICHAR( 'Z' )
      IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
!
!        ASCII character set
!
         IF( IC.GE.97 .AND. IC.LE.122 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 10 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.97 .AND. IC.LE.122 )                           &
     &            SUBNAM( I:I ) = CHAR( IC-32 )
   10       CONTINUE
         END IF
!
      ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
!
!        EBCDIC character set
!
         IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.                         &
     &       ( IC.GE.145 .AND. IC.LE.153 ) .OR.                         &
     &       ( IC.GE.162 .AND. IC.LE.169 ) ) THEN
            SUBNAM( 1:1 ) = CHAR( IC+64 )
            DO 20 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.                   &
     &             ( IC.GE.145 .AND. IC.LE.153 ) .OR.                   &
     &             ( IC.GE.162 .AND. IC.LE.169 ) )                      &
     &            SUBNAM( I:I ) = CHAR( IC+64 )
   20       CONTINUE
         END IF
!
      ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
!
!        Prime machines:  ASCII+128
!
         IF( IC.GE.225 .AND. IC.LE.250 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 30 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.225 .AND. IC.LE.250 )                          &
     &            SUBNAM( I:I ) = CHAR( IC-32 )
   30       CONTINUE
         END IF
      END IF
!
      C1 = SUBNAM( 1:1 )
      SNAME = C1.EQ.'S' .OR. C1.EQ.'D'
      CNAME = C1.EQ.'C' .OR. C1.EQ.'Z'
      IF( .NOT.( CNAME .OR. SNAME ) )                                   &
     &   RETURN
      C2 = SUBNAM( 2:3 )
      C3 = SUBNAM( 4:6 )
      C4 = C3( 2:3 )
!
      GO TO ( 110, 200, 300 ) ISPEC
!
  110 CONTINUE
!
!     ISPEC = 1:  block size
!
!     In these examples, separate code is provided for setting NB for
!     real and complex.  We assume that NB will take the same value in
!     single or double precision.
!
      NB = 1
!
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.    &
     &            C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'PO' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NB = 1
         ELSE IF( SNAME .AND. C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            NB = 64
         ELSE IF( C3.EQ.'TRD' ) THEN
            NB = 1
         ELSE IF( C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.         &
     &          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.         &
     &          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.         &
     &          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.         &
     &          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.         &
     &          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.         &
     &          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.         &
     &          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.         &
     &          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( C2.EQ.'GB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'PB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'TR' ) THEN
         IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'LA' ) THEN
         IF( C3.EQ.'UUM' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'ST' ) THEN
         IF( C3.EQ.'EBZ' ) THEN
            NB = 1
         END IF
      END IF
      ILAENV = NB
      RETURN
!
  200 CONTINUE
!
!     ISPEC = 2:  minimum block size
!
      NBMIN = 2
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.         &
     &       C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 8
            ELSE
               NBMIN = 8
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.         &
     &          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.         &
     &          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.         &
     &          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.         &
     &          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.         &
     &          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.         &
     &          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.         &
     &          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.         &
     &          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      END IF
      ILAENV = NBMIN
      RETURN
!
  300 CONTINUE
!
!     ISPEC = 3:  crossover point
!
      NX = 0
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.         &
     &       C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NX = 1
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NX = 1
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.         &
     &          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.         &
     &          C4.EQ.'BR' ) THEN
               NX = 128
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.         &
     &          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.         &
     &          C4.EQ.'BR' ) THEN
               NX = 128
            END IF
         END IF
      END IF
      ILAENV = NX
      RETURN
!
  400 CONTINUE
!
!     ISPEC = 4:  number of shifts (used by xHSEQR)
!
      ILAENV = 6
      RETURN
!
  500 CONTINUE
!
!     ISPEC = 5:  minimum column dimension (not used)
!
      ILAENV = 2
      RETURN
!
  600 CONTINUE
!
!     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
!
      ILAENV = INT( REAL( MIN( N1, N2 ) )*1.6E0 )
      RETURN
!
  700 CONTINUE
!
!     ISPEC = 7:  number of processors (not used)
!
      ILAENV = 1
      RETURN
!
  800 CONTINUE
!
!     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
!
      ILAENV = 50
      RETURN
!
!     End of ILAENV
!
      END

      integer function icamax(n,cx,incx)
!
!     finds the index of element having max. absolute value.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      complex cx(*)
      real smax
      integer i,incx,ix,n
      complex zdum
      real cabs1
      cabs1(zdum) = abs(real(zdum)) + abs(aimag(zdum))
!
      icamax = 0
      if( n.lt.1 .or. incx.le.0 ) return
      icamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
!
!        code for increment not equal to 1
!
      ix = 1
      smax = cabs1(cx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(cabs1(cx(ix)).le.smax) go to 5
         icamax = i
         smax = cabs1(cx(ix))
    5    ix = ix + incx
   10 continue
      return
!
!        code for increment equal to 1
!
   20 smax = cabs1(cx(1))
      do 30 i = 2,n
         if(cabs1(cx(i)).le.smax) go to 30
         icamax = i
         smax = cabs1(cx(i))
   30 continue
      return
      end

      LOGICAL          FUNCTION LSAME( CA, CB )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     January 31, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          CA, CB
!     ..
!
!  Purpose
!  =======
!
!  LSAME returns .TRUE. if CA is the same letter as CB regardless of
!  case.
!
!  Arguments
!  =========
!
!  CA      (input) CHARACTER*1
!  CB      (input) CHARACTER*1
!          CA and CB specify the single characters to be compared.
!
! ==================================================================
!
!     .. Intrinsic Functions ..
      INTRINSIC          ICHAR
!     ..
!     .. Local Scalars ..
      INTEGER            INTA, INTB, ZCODE
!     ..
!     .. Executable Statements ..
!
!     Test if the characters are equal
!
      LSAME = CA.EQ.CB
      IF( LSAME )                                                       &
     &   RETURN
!
!     Now test for equivalence if both characters are alphabetic.
!
      ZCODE = ICHAR( 'Z' )
!
!     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
!     machines, on which ICHAR returns a value with bit 8 set.
!     ICHAR('A') on Prime machines returns 193 which is the same as
!     ICHAR('A') on an EBCDIC machine.
!
      INTA = ICHAR( CA )
      INTB = ICHAR( CB )
!
      IF( ZCODE.EQ.90 .OR. ZCODE.EQ.122 ) THEN
!
!        ASCII is assumed - ZCODE is the ASCII code of either lower or
!        upper case 'Z'.
!
         IF( INTA.GE.97 .AND. INTA.LE.122 ) INTA = INTA - 32
         IF( INTB.GE.97 .AND. INTB.LE.122 ) INTB = INTB - 32
!
      ELSE IF( ZCODE.EQ.233 .OR. ZCODE.EQ.169 ) THEN
!
!        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower
!        or upper case 'Z'.
!
         IF( INTA.GE.129 .AND. INTA.LE.137 .OR.                         &
     &       INTA.GE.145 .AND. INTA.LE.153 .OR.                         &
     &       INTA.GE.162 .AND. INTA.LE.169 ) INTA = INTA + 64
         IF( INTB.GE.129 .AND. INTB.LE.137 .OR.                         &
     &       INTB.GE.145 .AND. INTB.LE.153 .OR.                         &
     &       INTB.GE.162 .AND. INTB.LE.169 ) INTB = INTB + 64
!
      ELSE IF( ZCODE.EQ.218 .OR. ZCODE.EQ.250 ) THEN
!
!        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
!        plus 128 of either lower or upper case 'Z'.
!
         IF( INTA.GE.225 .AND. INTA.LE.250 ) INTA = INTA - 32
         IF( INTB.GE.225 .AND. INTB.LE.250 ) INTB = INTB - 32
      END IF
      LSAME = INTA.EQ.INTB
!
!     RETURN
!
!     End of LSAME
!
      END
!       SUBROUTINE CTRTRI( UPLO, DIAG, N, A, LDA, INFO )
! *
! *  -- LAPACK routine (version 2.0) --
! *     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
! *     Courant Institute, Argonne National Lab, and Rice University
! *     September 30, 1994
! *
! *     .. Scalar Arguments ..
!       CHARACTER          DIAG, UPLO
!       INTEGER            INFO, LDA, N
! *     ..
! *     .. Array Arguments ..
!       COMPLEX            A( LDA, * )
! *     ..
! *
! *  Purpose
! *  =======
! *
! *  CTRTRI computes the inverse of a complex upper or lower
! *  triangular matrix A.
! *
! *  This is the Level 3 BLAS version of the algorithm.
! *
! *  Arguments
! *  =========
! *
! *  UPLO    (input) CHARACTER*1
! *          = 'U':  A is upper triangular;
! *          = 'L':  A is lower triangular.
! *
! *  DIAG    (input) CHARACTER*1
! *          = 'N':  A is non-unit triangular;
! *          = 'U':  A is unit triangular.
! *
! *  N       (input) INTEGER
! *          The order of the matrix A.  N >= 0.
! *
! *  A       (input/output) COMPLEX array, dimension (LDA,N)
! *       On entry, the triangular matrix A.  If UPLO = 'U', the
! *       leading N-by-N upper triangular part of the array A
! *       contains the upper triangular matrix, and the strictly lower
! *       triangular part of A is not referenced.  If UPLO = 'L', the
! *       leading N-by-N lower triangular part of the array A contains
! *       the lower triangular matrix, and the strictly upper
! *       triangular part of A is not referenced.  If DIAG = 'U', the
! *       diagonal elements of A are also not referenced and are
! *       assumed to be 1.
! *       On exit, the (triangular) inverse of the original matrix, in
! *       the same storage format.
! *
! *  LDA     (input) INTEGER
! *          The leading dimension of the array A.  LDA >= max(1,N).
! *
! *  INFO    (output) INTEGER
! *       = 0: successful exit
! *       < 0: if INFO = -i, the i-th argument had an illegal value
! *       > 0: if INFO = i, A(i,i) is exactly zero.  The triangular
! *            matrix is singular and its inverse can not be computed.
! *
! *  ================================================================
! *
! *     .. Parameters ..
!       COMPLEX            ONE, ZERO
!       PARAMETER          ( ONE = ( 1.0E+0, 0.0E+0 ),
!      $                   ZERO = ( 0.0E+0, 0.0E+0 ) )
! *     ..
! *     .. Local Scalars ..
!       LOGICAL            NOUNIT, UPPER
!       INTEGER            J, JB, NB, NN
! *     ..
! *     .. External Functions ..
!       LOGICAL            LSAME
!       INTEGER            ILAENV
!       EXTERNAL           LSAME, ILAENV
! *     ..
! *     .. External Subroutines ..
!       EXTERNAL           CTRMM, CTRSM, CTRTI2, XERBLA
! *     ..
! *     .. Intrinsic Functions ..
!       INTRINSIC          MAX, MIN
! *     ..
! *     .. Executable Statements ..
! *
! *     Test the input parameters.
! *
!       INFO = 0
!       UPPER = LSAME( UPLO, 'U' )
!       NOUNIT = LSAME( DIAG, 'N' )
!       IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
!          INFO = -1
!       ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
!          INFO = -2
!       ELSE IF( N.LT.0 ) THEN
!          INFO = -3
!       ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
!          INFO = -5
!       END IF
!       IF( INFO.NE.0 ) THEN
!          CALL XERBLA( 'CTRTRI', -INFO )
!          RETURN
!       END IF
! *
! *     Quick return if possible
! *
!       IF( N.EQ.0 )
!      $   RETURN
! *
! *     Check for singularity if non-unit.
! *
!       IF( NOUNIT ) THEN
!          DO 10 INFO = 1, N
!             IF( A( INFO, INFO ).EQ.ZERO )
!      $         RETURN
!    10    CONTINUE
!          INFO = 0
!       END IF
! *
! *     Determine the block size for this environment.
! *
!       NB = ILAENV( 1, 'CTRTRI', UPLO // DIAG, N, -1, -1, -1 )
!       IF( NB.LE.1 .OR. NB.GE.N ) THEN
! *
! *        Use unblocked code
! *
!          CALL CTRTI2( UPLO, DIAG, N, A, LDA, INFO )
!       ELSE
! *
! *        Use blocked code
! *
!          IF( UPPER ) THEN
! *
! *           Compute inverse of upper triangular matrix
! *
!             DO 20 J = 1, N, NB
!                JB = MIN( NB, N-J+1 )
! *
! *              Compute rows 1:j-1 of current block column
! *
!             CALL CTRMM( 'Left', 'Upper', 'No transpose', DIAG, J-1,
!      $                  JB, ONE, A, LDA, A( 1, J ), LDA )
!             CALL CTRSM( 'Right', 'Upper', 'No transpose', DIAG, J-1,
!      $                  JB, -ONE, A( J, J ), LDA, A( 1, J ), LDA )
! *
! *              Compute inverse of current diagonal block
! *
!             CALL CTRTI2( 'Upper', DIAG, JB, A( J, J ), LDA, INFO )
!    20       CONTINUE
!          ELSE
! *
! *           Compute inverse of lower triangular matrix
! *
!             NN = ( ( N-1 ) / NB )*NB + 1
!             DO 30 J = NN, 1, -NB
!                JB = MIN( NB, N-J+1 )
!                IF( J+JB.LE.N ) THEN
! *
! *                 Compute rows j+jb:n of current block column
! *
!                CALL CTRMM( 'Left', 'Lower', 'No transpose', DIAG,
!      $                     N-J-JB+1, JB, ONE, A( J+JB, J+JB ), LDA,
!      $                     A( J+JB, J ), LDA )
!                CALL CTRSM( 'Right', 'Lower', 'No transpose', DIAG,
!      $                     N-J-JB+1, JB, -ONE, A( J, J ), LDA,
!      $                     A( J+JB, J ), LDA )
!                END IF
! *
! *              Compute inverse of current diagonal block
! *
!             CALL CTRTI2( 'Lower', DIAG, JB, A( J, J ), LDA, INFO )
!    30       CONTINUE
!          END IF
!       END IF
! *
!       RETURN
! *
! *     End of CTRTRI
! *
!       END
      SUBROUTINE CTRTI2( UPLO, DIAG, N, A, LDA, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          DIAG, UPLO
      INTEGER            INFO, LDA, N
!     ..
!     .. Array Arguments ..
      COMPLEX            A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  CTRTI2 computes the inverse of a complex upper or lower triangular
!  matrix.
!
!  This is the Level 2 BLAS version of the algorithm.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          Specifies whether the matrix A is upper or lower
!          triangular.
!          = 'U':  Upper triangular
!          = 'L':  Lower triangular
!
!  DIAG    (input) CHARACTER*1
!          Specifies whether or not the matrix A is unit triangular.
!          = 'N':  Non-unit triangular
!          = 'U':  Unit triangular
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input/output) COMPLEX array, dimension (LDA,N)
!       On entry, the triangular matrix A.  If UPLO = 'U', the
!       leading n by n upper triangular part of the array A contains
!       the upper triangular matrix, and the strictly lower
!       triangular part of A is not referenced.  If UPLO = 'L', the
!       leading n by n lower triangular part of the array A contains
!       the lower triangular matrix, and the strictly upper
!       triangular part of A is not referenced.  If DIAG = 'U', the
!       diagonal elements of A are also not referenced and are
!       assumed to be 1.
!
!       On exit, the (triangular) inverse of the original matrix, in
!       the same storage format.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value
!
!  ==================================================================
!
!     .. Parameters ..
      COMPLEX            ONE
      PARAMETER          ( ONE = ( 1.0E+0, 0.0E+0 ) )
!     ..
!     .. Local Scalars ..
      LOGICAL            NOUNIT, UPPER
      INTEGER            J
      COMPLEX            AJJ
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           CSCAL, CTRMV, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      NOUNIT = LSAME( DIAG, 'N' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CTRTI2', -INFO )
         RETURN
      END IF
!
      IF( UPPER ) THEN
!
!        Compute inverse of upper triangular matrix.
!
         DO 10 J = 1, N
            IF( NOUNIT ) THEN
               A( J, J ) = ONE / A( J, J )
               AJJ = -A( J, J )
            ELSE
               AJJ = -ONE
            END IF
!
!           Compute elements 1:j-1 of j-th column.
!
            CALL CTRMV( 'Upper', 'No transpose', DIAG, J-1, A, LDA,     &
     &                  A( 1, J ), 1 )
            CALL CSCAL( J-1, AJJ, A( 1, J ), 1 )
   10    CONTINUE
      ELSE
!
!        Compute inverse of lower triangular matrix.
!
         DO 20 J = N, 1, -1
            IF( NOUNIT ) THEN
               A( J, J ) = ONE / A( J, J )
               AJJ = -A( J, J )
            ELSE
               AJJ = -ONE
            END IF
            IF( J.LT.N ) THEN
!
!              Compute elements j+1:n of j-th column.
!
               CALL CTRMV( 'Lower', 'No transpose', DIAG, N-J,          &
     &                     A( J+1, J+1 ), LDA, A( J+1, J ), 1 )
               CALL CSCAL( N-J, AJJ, A( J+1, J ), 1 )
            END IF
   20    CONTINUE
      END IF
!
      RETURN
!
!     End of CTRTI2
!
      END
      SUBROUTINE CTRMV ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
!     .. Scalar Arguments ..
      INTEGER            INCX, LDA, N
      CHARACTER*1        DIAG, TRANS, UPLO
!     .. Array Arguments ..
      COMPLEX            A( LDA, * ), X( * )
!     ..
!
!  Purpose
!  =======
!
!  CTRMV  performs one of the matrix-vector operations
!
!     x := A*x,   or   x := A'*x,   or   x := conjg( A' )*x,
!
!  where x is an n element vector and  A is an n by n unit, or
!  non-unit, upper or lower triangular matrix.
!
!  Parameters
!  ==========
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the matrix is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANS  - CHARACTER*1.
!           On entry, TRANS specifies the operation to be performed as
!           follows:
!
!              TRANS = 'N' or 'n'   x := A*x.
!
!              TRANS = 'T' or 't'   x := A'*x.
!
!              TRANS = 'C' or 'c'   x := conjg( A' )*x.
!
!           Unchanged on exit.
!
!  DIAG   - CHARACTER*1.
!           On entry, DIAG specifies whether or not A is unit
!           triangular as follows:
!
!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.
!
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  A      - COMPLEX          array of DIMENSION ( LDA, n ).
!        Before entry with  UPLO = 'U' or 'u', the leading n by n
!        upper triangular part of the array A must contain the upper
!        triangular matrix and the strictly lower triangular part of
!        A is not referenced.
!        Before entry with UPLO = 'L' or 'l', the leading n by n
!        lower triangular part of the array A must contain the lower
!        triangular matrix and the strictly upper triangular part of
!        A is not referenced.
!        Note that when  DIAG = 'U' or 'u', the diagonal elements of
!        A are not referenced either, but are assumed to be unity.
!        Unchanged on exit.
!
!  LDA    - INTEGER.
!        On entry, LDA specifies the first dimension of A as declared
!        in the calling (sub) program. LDA must be at least
!        max( 1, n ).
!        Unchanged on exit.
!
!  X      - COMPLEX          array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the n
!           element vector x. On exit, X is overwritten with the
!           tranformed vector x.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
!     .. Parameters ..
      COMPLEX            ZERO
      PARAMETER        ( ZERO = ( 0.0E+0, 0.0E+0 ) )
!     .. Local Scalars ..
      COMPLEX            TEMP
      INTEGER            I, INFO, IX, J, JX, KX
      LOGICAL            NOCONJ, NOUNIT
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC          CONJG, MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF     ( .NOT.LSAME( UPLO , 'U' ).AND.                            &
     &         .NOT.LSAME( UPLO , 'L' )      )THEN
         INFO = 1
      ELSE IF( .NOT.LSAME( TRANS, 'N' ).AND.                            &
     &         .NOT.LSAME( TRANS, 'T' ).AND.                            &
     &         .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 2
      ELSE IF( .NOT.LSAME( DIAG , 'U' ).AND.                            &
     &         .NOT.LSAME( DIAG , 'N' )      )THEN
         INFO = 3
      ELSE IF( N.LT.0 )THEN
         INFO = 4
      ELSE IF( LDA.LT.MAX( 1, N ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'CTRMV ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( N.EQ.0 )                                                      &
     &   RETURN
!
      NOCONJ = LSAME( TRANS, 'T' )
      NOUNIT = LSAME( DIAG , 'N' )
!
!     Set up the start point in X if the increment is not unity. This
!     will be  ( N - 1 )*INCX  too small for descending loops.
!
      IF( INCX.LE.0 )THEN
         KX = 1 - ( N - 1 )*INCX
      ELSE IF( INCX.NE.1 )THEN
         KX = 1
      END IF
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
      IF( LSAME( TRANS, 'N' ) )THEN
!
!        Form  x := A*x.
!
         IF( LSAME( UPLO, 'U' ) )THEN
            IF( INCX.EQ.1 )THEN
               DO 20, J = 1, N
                  IF( X( J ).NE.ZERO )THEN
                     TEMP = X( J )
                     DO 10, I = 1, J - 1
                        X( I ) = X( I ) + TEMP*A( I, J )
   10                CONTINUE
                     IF( NOUNIT )                                       &
     &                  X( J ) = X( J )*A( J, J )
                  END IF
   20          CONTINUE
            ELSE
               JX = KX
               DO 40, J = 1, N
                  IF( X( JX ).NE.ZERO )THEN
                     TEMP = X( JX )
                     IX   = KX
                     DO 30, I = 1, J - 1
                        X( IX ) = X( IX ) + TEMP*A( I, J )
                        IX      = IX      + INCX
   30                CONTINUE
                     IF( NOUNIT )                                       &
     &                  X( JX ) = X( JX )*A( J, J )
                  END IF
                  JX = JX + INCX
   40          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 60, J = N, 1, -1
                  IF( X( J ).NE.ZERO )THEN
                     TEMP = X( J )
                     DO 50, I = N, J + 1, -1
                        X( I ) = X( I ) + TEMP*A( I, J )
   50                CONTINUE
                     IF( NOUNIT )                                       &
     &                  X( J ) = X( J )*A( J, J )
                  END IF
   60          CONTINUE
            ELSE
               KX = KX + ( N - 1 )*INCX
               JX = KX
               DO 80, J = N, 1, -1
                  IF( X( JX ).NE.ZERO )THEN
                     TEMP = X( JX )
                     IX   = KX
                     DO 70, I = N, J + 1, -1
                        X( IX ) = X( IX ) + TEMP*A( I, J )
                        IX      = IX      - INCX
   70                CONTINUE
                     IF( NOUNIT )                                       &
     &                  X( JX ) = X( JX )*A( J, J )
                  END IF
                  JX = JX - INCX
   80          CONTINUE
            END IF
         END IF
      ELSE
!
!        Form  x := A'*x  or  x := conjg( A' )*x.
!
         IF( LSAME( UPLO, 'U' ) )THEN
            IF( INCX.EQ.1 )THEN
               DO 110, J = N, 1, -1
                  TEMP = X( J )
                  IF( NOCONJ )THEN
                     IF( NOUNIT )                                       &
     &                  TEMP = TEMP*A( J, J )
                     DO 90, I = J - 1, 1, -1
                        TEMP = TEMP + A( I, J )*X( I )
   90                CONTINUE
                  ELSE
                     IF( NOUNIT )                                       &
     &                  TEMP = TEMP*CONJG( A( J, J ) )
                     DO 100, I = J - 1, 1, -1
                        TEMP = TEMP + CONJG( A( I, J ) )*X( I )
  100                CONTINUE
                  END IF
                  X( J ) = TEMP
  110          CONTINUE
            ELSE
               JX = KX + ( N - 1 )*INCX
               DO 140, J = N, 1, -1
                  TEMP = X( JX )
                  IX   = JX
                  IF( NOCONJ )THEN
                     IF( NOUNIT )                                       &
     &                  TEMP = TEMP*A( J, J )
                     DO 120, I = J - 1, 1, -1
                        IX   = IX   - INCX
                        TEMP = TEMP + A( I, J )*X( IX )
  120                CONTINUE
                  ELSE
                     IF( NOUNIT )                                       &
     &                  TEMP = TEMP*CONJG( A( J, J ) )
                     DO 130, I = J - 1, 1, -1
                        IX   = IX   - INCX
                        TEMP = TEMP + CONJG( A( I, J ) )*X( IX )
  130                CONTINUE
                  END IF
                  X( JX ) = TEMP
                  JX      = JX   - INCX
  140          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 170, J = 1, N
                  TEMP = X( J )
                  IF( NOCONJ )THEN
                     IF( NOUNIT )                                       &
     &                  TEMP = TEMP*A( J, J )
                     DO 150, I = J + 1, N
                        TEMP = TEMP + A( I, J )*X( I )
  150                CONTINUE
                  ELSE
                     IF( NOUNIT )                                       &
     &                  TEMP = TEMP*CONJG( A( J, J ) )
                     DO 160, I = J + 1, N
                        TEMP = TEMP + CONJG( A( I, J ) )*X( I )
  160                CONTINUE
                  END IF
                  X( J ) = TEMP
  170          CONTINUE
            ELSE
               JX = KX
               DO 200, J = 1, N
                  TEMP = X( JX )
                  IX   = JX
                  IF( NOCONJ )THEN
                     IF( NOUNIT )                                       &
     &                  TEMP = TEMP*A( J, J )
                     DO 180, I = J + 1, N
                        IX   = IX   + INCX
                        TEMP = TEMP + A( I, J )*X( IX )
  180                CONTINUE
                  ELSE
                     IF( NOUNIT )                                       &
     &                  TEMP = TEMP*CONJG( A( J, J ) )
                     DO 190, I = J + 1, N
                        IX   = IX   + INCX
                        TEMP = TEMP + CONJG( A( I, J ) )*X( IX )
  190                CONTINUE
                  END IF
                  X( JX ) = TEMP
                  JX      = JX   + INCX
  200          CONTINUE
            END IF
         END IF
      END IF
!
      RETURN
!
!     End of CTRMV .
!
      END
      subroutine caxpy(n,ca,cx,incx,cy,incy)
!
!     constant times a vector plus a vector.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      complex cx(*),cy(*),ca
      integer i,incx,incy,ix,iy,n
!
      if(n.le.0)return
      if (abs(real(ca)) + abs(aimag(ca)) .eq. 0.0 ) return
      if(incx.eq.1.and.incy.eq.1)go to 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        cy(iy) = cy(iy) + ca*cx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!
!        code for both increments equal to 1
!
   20 do 30 i = 1,n
        cy(i) = cy(i) + ca*cx(i)
   30 continue
      return
      end
      subroutine  ccopy(n,cx,incx,cy,incy)
!
!     copies a vector, x, to a vector, y.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      complex cx(*),cy(*)
      integer i,incx,incy,ix,iy,n
!
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        cy(iy) = cx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!
!        code for both increments equal to 1
!
   20 do 30 i = 1,n
        cy(i) = cx(i)
   30 continue
      return
      end

      SUBROUTINE CGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX,           &
     &                   BETA, Y, INCY )
!     .. Scalar Arguments ..
      COMPLEX            ALPHA, BETA
      INTEGER            INCX, INCY, LDA, M, N
      CHARACTER*1        TRANS
!     .. Array Arguments ..
      COMPLEX            A( LDA, * ), X( * ), Y( * )
!     ..
!
!  Purpose
!  =======
!
!  CGEMV  performs one of the matrix-vector operations
!
!     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,   or
!
!     y := alpha*conjg( A' )*x + beta*y,
!
!  where alpha and beta are scalars, x and y are vectors and A is an
!  m by n matrix.
!
!  Parameters
!  ==========
!
!  TRANS  - CHARACTER*1.
!           On entry, TRANS specifies the operation to be performed as
!           follows:
!
!              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
!
!              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
!
!              TRANS = 'C' or 'c'   y := alpha*conjg( A' )*x + beta*y.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - COMPLEX         .
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - COMPLEX          array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit.
!
!  X      - COMPLEX          array of DIMENSION at least
!           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
!           Before entry, the incremented array X must contain the
!           vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  BETA   - COMPLEX         .
!           On entry, BETA specifies the scalar beta. When BETA is
!           supplied as zero then Y need not be set on input.
!           Unchanged on exit.
!
!  Y      - COMPLEX          array of DIMENSION at least
!           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
!           Before entry with BETA non-zero, the incremented array Y
!           must contain the vector y. On exit, Y is overwritten by the
!           updated vector y.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
!     .. Parameters ..
      COMPLEX            ONE
      PARAMETER        ( ONE  = ( 1.0E+0, 0.0E+0 ) )
      COMPLEX            ZERO
      PARAMETER        ( ZERO = ( 0.0E+0, 0.0E+0 ) )
!     .. Local Scalars ..
      COMPLEX            TEMP
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY, LENX, LENY
      LOGICAL            NOCONJ
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC          CONJG, MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF     ( .NOT.LSAME( TRANS, 'N' ).AND.                            &
     &         .NOT.LSAME( TRANS, 'T' ).AND.                            &
     &         .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 1
      ELSE IF( M.LT.0 )THEN
         INFO = 2
      ELSE IF( N.LT.0 )THEN
         INFO = 3
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'CGEMV ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.                                  &
     &    ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )                   &
     &   RETURN
!
      NOCONJ = LSAME( TRANS, 'T' )
!
!     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
!     up the start points in  X  and  Y.
!
      IF( LSAME( TRANS, 'N' ) )THEN
         LENX = N
         LENY = M
      ELSE
         LENX = M
         LENY = N
      END IF
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( LENX - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( LENY - 1 )*INCY
      END IF
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
!     First form  y := beta*y.
!
      IF( BETA.NE.ONE )THEN
         IF( INCY.EQ.1 )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, LENY
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, LENY
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, LENY
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, LENY
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO )                                               &
     &   RETURN
      IF( LSAME( TRANS, 'N' ) )THEN
!
!        Form  y := alpha*A*x + y.
!
         JX = KX
         IF( INCY.EQ.1 )THEN
            DO 60, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  DO 50, I = 1, M
                     Y( I ) = Y( I ) + TEMP*A( I, J )
   50             CONTINUE
               END IF
               JX = JX + INCX
   60       CONTINUE
         ELSE
            DO 80, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IY   = KY
                  DO 70, I = 1, M
                     Y( IY ) = Y( IY ) + TEMP*A( I, J )
                     IY      = IY      + INCY
   70             CONTINUE
               END IF
               JX = JX + INCX
   80       CONTINUE
         END IF
      ELSE
!
!        Form  y := alpha*A'*x + y  or  y := alpha*conjg( A' )*x + y.
!
         JY = KY
         IF( INCX.EQ.1 )THEN
            DO 110, J = 1, N
               TEMP = ZERO
               IF( NOCONJ )THEN
                  DO 90, I = 1, M
                     TEMP = TEMP + A( I, J )*X( I )
   90             CONTINUE
               ELSE
                  DO 100, I = 1, M
                     TEMP = TEMP + CONJG( A( I, J ) )*X( I )
  100             CONTINUE
               END IF
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  110       CONTINUE
         ELSE
            DO 140, J = 1, N
               TEMP = ZERO
               IX   = KX
               IF( NOCONJ )THEN
                  DO 120, I = 1, M
                     TEMP = TEMP + A( I, J )*X( IX )
                     IX   = IX   + INCX
  120             CONTINUE
               ELSE
                  DO 130, I = 1, M
                     TEMP = TEMP + CONJG( A( I, J ) )*X( IX )
                     IX   = IX   + INCX
  130             CONTINUE
               END IF
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  140       CONTINUE
         END IF
      END IF
!
      RETURN
!
!     End of CGEMV .
!
      END
      double precision function dcabs1(z)
      double complex z,zz
      double precision t(2)
      equivalence (zz,t(1))
      zz = z
      dcabs1 = dabs(t(1)) + dabs(t(2))
      return
      end
      SUBROUTINE DGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB,&
     &                   BETA, C, LDC )
!     .. Scalar Arguments ..
      CHARACTER*1        TRANSA, TRANSB
      INTEGER            M, N, K, LDA, LDB, LDC
      DOUBLE PRECISION   ALPHA, BETA
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * )
!     ..
!
!  Purpose
!  =======
!
!  DGEMM  performs one of the matrix-matrix operations
!
!     C := alpha*op( A )*op( B ) + beta*C,
!
!  where  op( X ) is one of
!
!     op( X ) = X   or   op( X ) = X',
!
!  alpha and beta are scalars, and A, B and C are matrices, with op( A )
!  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
!
!  Parameters
!  ==========
!
!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n',  op( A ) = A.
!
!              TRANSA = 'T' or 't',  op( A ) = A'.
!
!              TRANSA = 'C' or 'c',  op( A ) = A'.
!
!           Unchanged on exit.
!
!  TRANSB - CHARACTER*1.
!           On entry, TRANSB specifies the form of op( B ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSB = 'N' or 'n',  op( B ) = B.
!
!              TRANSB = 'T' or 't',  op( B ) = B'.
!
!              TRANSB = 'C' or 'c',  op( B ) = B'.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry,  M  specifies  the number  of rows  of the  matrix
!           op( A )  and of the  matrix  C.  M  must  be at least  zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry,  N  specifies the number  of columns of the matrix
!           op( B ) and the number of columns of the matrix C. N must be
!           at least zero.
!           Unchanged on exit.
!
!  K      - INTEGER.
!           On entry,  K  specifies  the number of columns of the matrix
!           op( A ) and the number of rows of the matrix op( B ). K must
!           be at least  zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
!           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
!           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
!           part of the array  A  must contain the matrix  A,  otherwise
!           the leading  k by m  part of the array  A  must contain  the
!           matrix A.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
!           LDA must be at least  max( 1, m ), otherwise  LDA must be at
!           least  max( 1, k ).
!           Unchanged on exit.
!
!  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
!           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
!           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
!           part of the array  B  must contain the matrix  B,  otherwise
!           the leading  n by k  part of the array  B  must contain  the
!           matrix B.
!           Unchanged on exit.
!
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
!           LDB must be at least  max( 1, k ), otherwise  LDB must be at
!           least  max( 1, n ).
!           Unchanged on exit.
!
!  BETA   - DOUBLE PRECISION.
!           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
!           supplied as zero then C need not be set on input.
!           Unchanged on exit.
!
!  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
!           Before entry, the leading  m by n  part of the array  C must
!           contain the matrix  C,  except when  beta  is zero, in which
!           case C need not be set on entry.
!           On exit, the array  C  is overwritten by the  m by n  matrix
!           ( alpha*op( A )*op( B ) + beta*C ).
!
!  LDC    - INTEGER.
!           On entry, LDC specifies the first dimension of C as declared
!           in  the  calling  (sub)  program.   LDC  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 3 Blas routine.
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     .. Local Scalars ..
      LOGICAL            NOTA, NOTB
      INTEGER            I, INFO, J, L, NCOLA, NROWA, NROWB
      DOUBLE PRECISION   TEMP
!     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Executable Statements ..
!
!     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
!     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
!     and  columns of  A  and the  number of  rows  of  B  respectively.
!
      NOTA  = LSAME( TRANSA, 'N' )
      NOTB  = LSAME( TRANSB, 'N' )
      IF( NOTA )THEN
         NROWA = M
         NCOLA = K
      ELSE
         NROWA = K
         NCOLA = M
      END IF
      IF( NOTB )THEN
         NROWB = K
      ELSE
         NROWB = N
      END IF
!
!     Test the input parameters.
!
      INFO = 0
      IF(      ( .NOT.NOTA                 ).AND.                       &
     &         ( .NOT.LSAME( TRANSA, 'C' ) ).AND.                       &
     &         ( .NOT.LSAME( TRANSA, 'T' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.NOTB                 ).AND.                       &
     &         ( .NOT.LSAME( TRANSB, 'C' ) ).AND.                       &
     &         ( .NOT.LSAME( TRANSB, 'T' ) )      )THEN
         INFO = 2
      ELSE IF( M  .LT.0               )THEN
         INFO = 3
      ELSE IF( N  .LT.0               )THEN
         INFO = 4
      ELSE IF( K  .LT.0               )THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 8
      ELSE IF( LDB.LT.MAX( 1, NROWB ) )THEN
         INFO = 10
      ELSE IF( LDC.LT.MAX( 1, M     ) )THEN
         INFO = 13
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGEMM ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.                                  &
     &    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) ) &
     &   RETURN
!
!     And if  alpha.eq.zero.
!
      IF( ALPHA.EQ.ZERO )THEN
         IF( BETA.EQ.ZERO )THEN
            DO 20, J = 1, N
               DO 10, I = 1, M
                  C( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               DO 30, I = 1, M
                  C( I, J ) = BETA*C( I, J )
   30          CONTINUE
   40       CONTINUE
         END IF
         RETURN
      END IF
!
!     Start the operations.
!
      IF( NOTB )THEN
         IF( NOTA )THEN
!
!           Form  C := alpha*A*B + beta*C.
!
            DO 90, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 50, I = 1, M
                     C( I, J ) = ZERO
   50             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 60, I = 1, M
                     C( I, J ) = BETA*C( I, J )
   60             CONTINUE
               END IF
               DO 80, L = 1, K
                  IF( B( L, J ).NE.ZERO )THEN
                     TEMP = ALPHA*B( L, J )
                     DO 70, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
   70                CONTINUE
                  END IF
   80          CONTINUE
   90       CONTINUE
         ELSE
!
!           Form  C := alpha*A'*B + beta*C
!
            DO 120, J = 1, N
               DO 110, I = 1, M
                  TEMP = ZERO
                  DO 100, L = 1, K
                     TEMP = TEMP + A( L, I )*B( L, J )
  100             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  110          CONTINUE
  120       CONTINUE
         END IF
      ELSE
         IF( NOTA )THEN
!
!           Form  C := alpha*A*B' + beta*C
!
            DO 170, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 130, I = 1, M
                     C( I, J ) = ZERO
  130             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 140, I = 1, M
                     C( I, J ) = BETA*C( I, J )
  140             CONTINUE
               END IF
               DO 160, L = 1, K
                  IF( B( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*B( J, L )
                     DO 150, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  150                CONTINUE
                  END IF
  160          CONTINUE
  170       CONTINUE
         ELSE
!
!           Form  C := alpha*A'*B' + beta*C
!
            DO 200, J = 1, N
               DO 190, I = 1, M
                  TEMP = ZERO
                  DO 180, L = 1, K
                     TEMP = TEMP + A( L, I )*B( J, L )
  180             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  190          CONTINUE
  200       CONTINUE
         END IF
      END IF
!
      RETURN
!
!     End of DGEMM .
!
      END
      SUBROUTINE DGER  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
!     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA
      INTEGER            INCX, INCY, LDA, M, N
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
!     ..
!
!  Purpose
!  =======
!
!  DGER   performs the rank 1 operation
!
!     A := alpha*x*y' + A,
!
!  where alpha is a scalar, x is an m element vector, y is an n element
!  vector and A is an m by n matrix.
!
!  Parameters
!  ==========
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  X      - DOUBLE PRECISION array of dimension at least
!           ( 1 + ( m - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the m
!           element vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  Y      - DOUBLE PRECISION array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCY ) ).
!           Before entry, the incremented array Y must contain the n
!           element vector y.
!           Unchanged on exit.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients. On exit, A is
!           overwritten by the updated matrix.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
!     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, J, JY, KX
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF     ( M.LT.0 )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 7
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGER  ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )               &
     &   RETURN
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
      IF( INCY.GT.0 )THEN
         JY = 1
      ELSE
         JY = 1 - ( N - 1 )*INCY
      END IF
      IF( INCX.EQ.1 )THEN
         DO 20, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               DO 10, I = 1, M
                  A( I, J ) = A( I, J ) + X( I )*TEMP
   10          CONTINUE
            END IF
            JY = JY + INCY
   20    CONTINUE
      ELSE
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( M - 1 )*INCX
         END IF
         DO 40, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               IX   = KX
               DO 30, I = 1, M
                  A( I, J ) = A( I, J ) + X( IX )*TEMP
                  IX        = IX        + INCX
   30          CONTINUE
            END IF
            JY = JY + INCY
   40    CONTINUE
      END IF
!
      RETURN
!
!     End of DGER  .
!
      END
      subroutine  dscal(n,da,dx,incx)
!
!     scales a vector by a constant.
!     uses unrolled loops for increment equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      double precision da,dx(*)
      integer i,incx,m,mp1,n,nincx
!
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
!
!        code for increment not equal to 1
!
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dx(i) = da*dx(i)
   10 continue
      return
!
!        code for increment equal to 1
!
!
!        clean-up loop
!
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dx(i) = da*dx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
   50 continue
      return
      end
      subroutine  dswap (n,dx,incx,dy,incy)
!
!     interchanges two vectors.
!     uses unrolled loops for increments equal one.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
!
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
!
!       code for unequal increments or equal increments not equal
!         to 1
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dx(ix)
        dx(ix) = dy(iy)
        dy(iy) = dtemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!
!       code for both increments equal to 1
!
!
!       clean-up loop
!
   20 m = mod(n,3)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
   30 continue
      if( n .lt. 3 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,3
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
        dtemp = dx(i + 1)
        dx(i + 1) = dy(i + 1)
        dy(i + 1) = dtemp
        dtemp = dx(i + 2)
        dx(i + 2) = dy(i + 2)
        dy(i + 2) = dtemp
   50 continue
      return
      end
      SUBROUTINE DTRSM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, &
     &                   B, LDB )
!     .. Scalar Arguments ..
      CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, LDA, LDB
      DOUBLE PRECISION   ALPHA
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  DTRSM  solves one of the matrix equations
!
!     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
!
!  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
!  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
!
!     op( A ) = A   or   op( A ) = A'.
!
!  The matrix X is overwritten on B.
!
!  Parameters
!  ==========
!
!  SIDE   - CHARACTER*1.
!           On entry, SIDE specifies whether op( A ) appears on the left
!           or right of X as follows:
!
!              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
!
!              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
!
!           Unchanged on exit.
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the matrix A is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n'   op( A ) = A.
!
!              TRANSA = 'T' or 't'   op( A ) = A'.
!
!              TRANSA = 'C' or 'c'   op( A ) = A'.
!
!           Unchanged on exit.
!
!  DIAG   - CHARACTER*1.
!           On entry, DIAG specifies whether or not A is unit triangular
!           as follows:
!
!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of B. M must be at
!           least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of B.  N must be
!           at least zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
!           zero then  A is not referenced and  B need not be set before
!           entry.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m
!           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
!           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
!           upper triangular part of the array  A must contain the upper
!           triangular matrix  and the strictly lower triangular part of
!           A is not referenced.
!           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
!           lower triangular part of the array  A must contain the lower
!           triangular matrix  and the strictly upper triangular part of
!           A is not referenced.
!           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
!           A  are not referenced either,  but are assumed to be  unity.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
!           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
!           then LDA must be at least max( 1, n ).
!           Unchanged on exit.
!
!  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).
!           Before entry,  the leading  m by n part of the array  B must
!           contain  the  right-hand  side  matrix  B,  and  on exit  is
!           overwritten by the solution matrix  X.
!
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in  the  calling  (sub)  program.   LDB  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 3 Blas routine.
!
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     .. Local Scalars ..
      LOGICAL            LSIDE, NOUNIT, UPPER
      INTEGER            I, INFO, J, K, NROWA
      DOUBLE PRECISION   TEMP
!     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      LSIDE  = LSAME( SIDE  , 'L' )
      IF( LSIDE )THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      NOUNIT = LSAME( DIAG  , 'N' )
      UPPER  = LSAME( UPLO  , 'U' )
!
      INFO   = 0
      IF(      ( .NOT.LSIDE                ).AND.                       &
     &         ( .NOT.LSAME( SIDE  , 'R' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.UPPER                ).AND.                       &
     &         ( .NOT.LSAME( UPLO  , 'L' ) )      )THEN
         INFO = 2
      ELSE IF( ( .NOT.LSAME( TRANSA, 'N' ) ).AND.                       &
     &         ( .NOT.LSAME( TRANSA, 'T' ) ).AND.                       &
     &         ( .NOT.LSAME( TRANSA, 'C' ) )      )THEN
         INFO = 3
      ELSE IF( ( .NOT.LSAME( DIAG  , 'U' ) ).AND.                       &
     &         ( .NOT.LSAME( DIAG  , 'N' ) )      )THEN
         INFO = 4
      ELSE IF( M  .LT.0               )THEN
         INFO = 5
      ELSE IF( N  .LT.0               )THEN
         INFO = 6
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 9
      ELSE IF( LDB.LT.MAX( 1, M     ) )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DTRSM ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( N.EQ.0 )                                                      &
     &   RETURN
!
!     And when  alpha.eq.zero.
!
      IF( ALPHA.EQ.ZERO )THEN
         DO 20, J = 1, N
            DO 10, I = 1, M
               B( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
         RETURN
      END IF
!
!     Start the operations.
!
      IF( LSIDE )THEN
         IF( LSAME( TRANSA, 'N' ) )THEN
!
!           Form  B := alpha*inv( A )*B.
!
            IF( UPPER )THEN
               DO 60, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 30, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   30                CONTINUE
                  END IF
                  DO 50, K = M, 1, -1
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT )                                    &
     &                     B( K, J ) = B( K, J )/A( K, K )
                        DO 40, I = 1, K - 1
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   40                   CONTINUE
                     END IF
   50             CONTINUE
   60          CONTINUE
            ELSE
               DO 100, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 70, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   70                CONTINUE
                  END IF
                  DO 90 K = 1, M
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT )                                    &
     &                     B( K, J ) = B( K, J )/A( K, K )
                        DO 80, I = K + 1, M
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   80                   CONTINUE
                     END IF
   90             CONTINUE
  100          CONTINUE
            END IF
         ELSE
!
!           Form  B := alpha*inv( A' )*B.
!
            IF( UPPER )THEN
               DO 130, J = 1, N
                  DO 120, I = 1, M
                     TEMP = ALPHA*B( I, J )
                     DO 110, K = 1, I - 1
                        TEMP = TEMP - A( K, I )*B( K, J )
  110                CONTINUE
                     IF( NOUNIT )                                       &
     &                  TEMP = TEMP/A( I, I )
                     B( I, J ) = TEMP
  120             CONTINUE
  130          CONTINUE
            ELSE
               DO 160, J = 1, N
                  DO 150, I = M, 1, -1
                     TEMP = ALPHA*B( I, J )
                     DO 140, K = I + 1, M
                        TEMP = TEMP - A( K, I )*B( K, J )
  140                CONTINUE
                     IF( NOUNIT )                                       &
     &                  TEMP = TEMP/A( I, I )
                     B( I, J ) = TEMP
  150             CONTINUE
  160          CONTINUE
            END IF
         END IF
      ELSE
         IF( LSAME( TRANSA, 'N' ) )THEN
!
!           Form  B := alpha*B*inv( A ).
!
            IF( UPPER )THEN
               DO 210, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 170, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  170                CONTINUE
                  END IF
                  DO 190, K = 1, J - 1
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 180, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  180                   CONTINUE
                     END IF
  190             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 200, I = 1, M
                        B( I, J ) = TEMP*B( I, J )
  200                CONTINUE
                  END IF
  210          CONTINUE
            ELSE
               DO 260, J = N, 1, -1
                  IF( ALPHA.NE.ONE )THEN
                     DO 220, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  220                CONTINUE
                  END IF
                  DO 240, K = J + 1, N
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 230, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  230                   CONTINUE
                     END IF
  240             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 250, I = 1, M
                       B( I, J ) = TEMP*B( I, J )
  250                CONTINUE
                  END IF
  260          CONTINUE
            END IF
         ELSE
!
!           Form  B := alpha*B*inv( A' ).
!
            IF( UPPER )THEN
               DO 310, K = N, 1, -1
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( K, K )
                     DO 270, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  270                CONTINUE
                  END IF
                  DO 290, J = 1, K - 1
                     IF( A( J, K ).NE.ZERO )THEN
                        TEMP = A( J, K )
                        DO 280, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  280                   CONTINUE
                     END IF
  290             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 300, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  300                CONTINUE
                  END IF
  310          CONTINUE
            ELSE
               DO 360, K = 1, N
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( K, K )
                     DO 320, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  320                CONTINUE
                  END IF
                  DO 340, J = K + 1, N
                     IF( A( J, K ).NE.ZERO )THEN
                        TEMP = A( J, K )
                        DO 330, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  330                   CONTINUE
                     END IF
  340             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 350, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  350                CONTINUE
                  END IF
  360          CONTINUE
            END IF
         END IF
      END IF
!
      RETURN
!
!     End of DTRSM .
!
      END
      integer function idamax(n,dx,incx)
!
!     finds the index of element having max. absolute value.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      double precision dx(*),dmax
      integer i,incx,ix,n
!
      idamax = 0
      if( n.lt.1 .or. incx.le.0 ) return
      idamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
!
!        code for increment not equal to 1
!
      ix = 1
      dmax = dabs(dx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(dabs(dx(ix)).le.dmax) go to 5
         idamax = i
         dmax = dabs(dx(ix))
    5    ix = ix + incx
   10 continue
      return
!
!        code for increment equal to 1
!
   20 dmax = dabs(dx(1))
      do 30 i = 2,n
         if(dabs(dx(i)).le.dmax) go to 30
         idamax = i
         dmax = dabs(dx(i))
   30 continue
      return
      end
      integer function izamax(n,zx,incx)
!
!     finds the index of element having max. absolute value.
!     jack dongarra, 1/15/85.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      double complex zx(*)
      double precision smax
      integer i,incx,ix,n
      double precision dcabs1
!
      izamax = 0
      if( n.lt.1 .or. incx.le.0 )return
      izamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
!
!        code for increment not equal to 1
!
      ix = 1
      smax = dcabs1(zx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(dcabs1(zx(ix)).le.smax) go to 5
         izamax = i
         smax = dcabs1(zx(ix))
    5    ix = ix + incx
   10 continue
      return
!
!        code for increment equal to 1
!
   20 smax = dcabs1(zx(1))
      do 30 i = 2,n
         if(dcabs1(zx(i)).le.smax) go to 30
         izamax = i
         smax = dcabs1(zx(i))
   30 continue
      return
      end
      subroutine  zcopy(n,zx,incx,zy,incy)
!
!     copies a vector, x, to a vector, y.
!     jack dongarra, linpack, 4/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      double complex zx(*),zy(*)
      integer i,incx,incy,ix,iy,n
!
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        zy(iy) = zx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!
!        code for both increments equal to 1
!
   20 do 30 i = 1,n
        zy(i) = zx(i)
   30 continue
      return
      end
      SUBROUTINE ZGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB,&
     &                   BETA, C, LDC )
!     .. Scalar Arguments ..
      CHARACTER*1        TRANSA, TRANSB
      INTEGER            M, N, K, LDA, LDB, LDC
      COMPLEX*16         ALPHA, BETA
!     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), B( LDB, * ), C( LDC, * )
!     ..
!
!  Purpose
!  =======
!
!  ZGEMM  performs one of the matrix-matrix operations
!
!     C := alpha*op( A )*op( B ) + beta*C,
!
!  where  op( X ) is one of
!
!     op( X ) = X   or   op( X ) = X'   or   op( X ) = conjg( X' ),
!
!  alpha and beta are scalars, and A, B and C are matrices, with op( A )
!  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
!
!  Parameters
!  ==========
!
!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n',  op( A ) = A.
!
!              TRANSA = 'T' or 't',  op( A ) = A'.
!
!              TRANSA = 'C' or 'c',  op( A ) = conjg( A' ).
!
!           Unchanged on exit.
!
!  TRANSB - CHARACTER*1.
!           On entry, TRANSB specifies the form of op( B ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSB = 'N' or 'n',  op( B ) = B.
!
!              TRANSB = 'T' or 't',  op( B ) = B'.
!
!              TRANSB = 'C' or 'c',  op( B ) = conjg( B' ).
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry,  M  specifies  the number  of rows  of the  matrix
!           op( A )  and of the  matrix  C.  M  must  be at least  zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry,  N  specifies the number  of columns of the matrix
!           op( B ) and the number of columns of the matrix C. N must be
!           at least zero.
!           Unchanged on exit.
!
!  K      - INTEGER.
!           On entry,  K  specifies  the number of columns of the matrix
!           op( A ) and the number of rows of the matrix op( B ). K must
!           be at least  zero.
!           Unchanged on exit.
!
!  ALPHA  - COMPLEX*16      .
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - COMPLEX*16       array of DIMENSION ( LDA, ka ), where ka is
!           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
!           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
!           part of the array  A  must contain the matrix  A,  otherwise
!           the leading  k by m  part of the array  A  must contain  the
!           matrix A.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
!           LDA must be at least  max( 1, m ), otherwise  LDA must be at
!           least  max( 1, k ).
!           Unchanged on exit.
!
!  B      - COMPLEX*16       array of DIMENSION ( LDB, kb ), where kb is
!           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
!           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
!           part of the array  B  must contain the matrix  B,  otherwise
!           the leading  n by k  part of the array  B  must contain  the
!           matrix B.
!           Unchanged on exit.
!
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
!           LDB must be at least  max( 1, k ), otherwise  LDB must be at
!           least  max( 1, n ).
!           Unchanged on exit.
!
!  BETA   - COMPLEX*16      .
!           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
!           supplied as zero then C need not be set on input.
!           Unchanged on exit.
!
!  C      - COMPLEX*16       array of DIMENSION ( LDC, n ).
!           Before entry, the leading  m by n  part of the array  C must
!           contain the matrix  C,  except when  beta  is zero, in which
!           case C need not be set on entry.
!           On exit, the array  C  is overwritten by the  m by n  matrix
!           ( alpha*op( A )*op( B ) + beta*C ).
!
!  LDC    - INTEGER.
!           On entry, LDC specifies the first dimension of C as declared
!           in  the  calling  (sub)  program.   LDC  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 3 Blas routine.
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC          DCONJG, MAX
!     .. Local Scalars ..
      LOGICAL            CONJA, CONJB, NOTA, NOTB
      INTEGER            I, INFO, J, L, NCOLA, NROWA, NROWB
      COMPLEX*16         TEMP
!     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER        ( ONE  = ( 1.0D+0, 0.0D+0 ) )
      COMPLEX*16         ZERO
      PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) )
!     ..
!     .. Executable Statements ..
!
!     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
!     conjugated or transposed, set  CONJA and CONJB  as true if  A  and
!     B  respectively are to be  transposed but  not conjugated  and set
!     NROWA, NCOLA and  NROWB  as the number of rows and  columns  of  A
!     and the number of rows of  B  respectively.
!
      NOTA  = LSAME( TRANSA, 'N' )
      NOTB  = LSAME( TRANSB, 'N' )
      CONJA = LSAME( TRANSA, 'C' )
      CONJB = LSAME( TRANSB, 'C' )
      IF( NOTA )THEN
         NROWA = M
         NCOLA = K
      ELSE
         NROWA = K
         NCOLA = M
      END IF
      IF( NOTB )THEN
         NROWB = K
      ELSE
         NROWB = N
      END IF
!
!     Test the input parameters.
!
      INFO = 0
      IF(      ( .NOT.NOTA                 ).AND.                       &
     &         ( .NOT.CONJA                ).AND.                       &
     &         ( .NOT.LSAME( TRANSA, 'T' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.NOTB                 ).AND.                       &
     &         ( .NOT.CONJB                ).AND.                       &
     &         ( .NOT.LSAME( TRANSB, 'T' ) )      )THEN
         INFO = 2
      ELSE IF( M  .LT.0               )THEN
         INFO = 3
      ELSE IF( N  .LT.0               )THEN
         INFO = 4
      ELSE IF( K  .LT.0               )THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 8
      ELSE IF( LDB.LT.MAX( 1, NROWB ) )THEN
         INFO = 10
      ELSE IF( LDC.LT.MAX( 1, M     ) )THEN
         INFO = 13
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'ZGEMM ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.                                  &
     &    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) ) &
     &   RETURN
!
!     And when  alpha.eq.zero.
!
      IF( ALPHA.EQ.ZERO )THEN
         IF( BETA.EQ.ZERO )THEN
            DO 20, J = 1, N
               DO 10, I = 1, M
                  C( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               DO 30, I = 1, M
                  C( I, J ) = BETA*C( I, J )
   30          CONTINUE
   40       CONTINUE
         END IF
         RETURN
      END IF
!
!     Start the operations.
!
      IF( NOTB )THEN
         IF( NOTA )THEN
!
!           Form  C := alpha*A*B + beta*C.
!
            DO 90, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 50, I = 1, M
                     C( I, J ) = ZERO
   50             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 60, I = 1, M
                     C( I, J ) = BETA*C( I, J )
   60             CONTINUE
               END IF
               DO 80, L = 1, K
                  IF( B( L, J ).NE.ZERO )THEN
                     TEMP = ALPHA*B( L, J )
                     DO 70, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
   70                CONTINUE
                  END IF
   80          CONTINUE
   90       CONTINUE
         ELSE IF( CONJA )THEN
!
!           Form  C := alpha*conjg( A' )*B + beta*C.
!
            DO 120, J = 1, N
               DO 110, I = 1, M
                  TEMP = ZERO
                  DO 100, L = 1, K
                     TEMP = TEMP + DCONJG( A( L, I ) )*B( L, J )
  100             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  110          CONTINUE
  120       CONTINUE
         ELSE
!
!           Form  C := alpha*A'*B + beta*C
!
            DO 150, J = 1, N
               DO 140, I = 1, M
                  TEMP = ZERO
                  DO 130, L = 1, K
                     TEMP = TEMP + A( L, I )*B( L, J )
  130             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  140          CONTINUE
  150       CONTINUE
         END IF
      ELSE IF( NOTA )THEN
         IF( CONJB )THEN
!
!           Form  C := alpha*A*conjg( B' ) + beta*C.
!
            DO 200, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 160, I = 1, M
                     C( I, J ) = ZERO
  160             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 170, I = 1, M
                     C( I, J ) = BETA*C( I, J )
  170             CONTINUE
               END IF
               DO 190, L = 1, K
                  IF( B( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*DCONJG( B( J, L ) )
                     DO 180, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  180                CONTINUE
                  END IF
  190          CONTINUE
  200       CONTINUE
         ELSE
!
!           Form  C := alpha*A*B'          + beta*C
!
            DO 250, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 210, I = 1, M
                     C( I, J ) = ZERO
  210             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 220, I = 1, M
                     C( I, J ) = BETA*C( I, J )
  220             CONTINUE
               END IF
               DO 240, L = 1, K
                  IF( B( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*B( J, L )
                     DO 230, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  230                CONTINUE
                  END IF
  240          CONTINUE
  250       CONTINUE
         END IF
      ELSE IF( CONJA )THEN
         IF( CONJB )THEN
!
!           Form  C := alpha*conjg( A' )*conjg( B' ) + beta*C.
!
            DO 280, J = 1, N
               DO 270, I = 1, M
                  TEMP = ZERO
                  DO 260, L = 1, K
                     TEMP = TEMP +                                      &
     &                      DCONJG( A( L, I ) )*DCONJG( B( J, L ) )
  260             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  270          CONTINUE
  280       CONTINUE
         ELSE
!
!           Form  C := alpha*conjg( A' )*B' + beta*C
!
            DO 310, J = 1, N
               DO 300, I = 1, M
                  TEMP = ZERO
                  DO 290, L = 1, K
                     TEMP = TEMP + DCONJG( A( L, I ) )*B( J, L )
  290             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  300          CONTINUE
  310       CONTINUE
         END IF
      ELSE
         IF( CONJB )THEN
!
!           Form  C := alpha*A'*conjg( B' ) + beta*C
!
            DO 340, J = 1, N
               DO 330, I = 1, M
                  TEMP = ZERO
                  DO 320, L = 1, K
                     TEMP = TEMP + A( L, I )*DCONJG( B( J, L ) )
  320             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  330          CONTINUE
  340       CONTINUE
         ELSE
!
!           Form  C := alpha*A'*B' + beta*C
!
            DO 370, J = 1, N
               DO 360, I = 1, M
                  TEMP = ZERO
                  DO 350, L = 1, K
                     TEMP = TEMP + A( L, I )*B( J, L )
  350             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  360          CONTINUE
  370       CONTINUE
         END IF
      END IF
!
      RETURN
!
!     End of ZGEMM .
!
      END
      SUBROUTINE ZGERU ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
!     .. Scalar Arguments ..
      COMPLEX*16         ALPHA
      INTEGER            INCX, INCY, LDA, M, N
!     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), X( * ), Y( * )
!     ..
!
!  Purpose
!  =======
!
!  ZGERU  performs the rank 1 operation
!
!     A := alpha*x*y' + A,
!
!  where alpha is a scalar, x is an m element vector, y is an n element
!  vector and A is an m by n matrix.
!
!  Parameters
!  ==========
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - COMPLEX*16      .
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  X      - COMPLEX*16       array of dimension at least
!           ( 1 + ( m - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the m
!           element vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  Y      - COMPLEX*16       array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCY ) ).
!           Before entry, the incremented array Y must contain the n
!           element vector y.
!           Unchanged on exit.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!  A      - COMPLEX*16       array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients. On exit, A is
!           overwritten by the updated matrix.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
!     .. Parameters ..
      COMPLEX*16         ZERO
      PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) )
!     .. Local Scalars ..
      COMPLEX*16         TEMP
      INTEGER            I, INFO, IX, J, JY, KX
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF     ( M.LT.0 )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 7
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'ZGERU ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )               &
     &   RETURN
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
      IF( INCY.GT.0 )THEN
         JY = 1
      ELSE
         JY = 1 - ( N - 1 )*INCY
      END IF
      IF( INCX.EQ.1 )THEN
         DO 20, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               DO 10, I = 1, M
                  A( I, J ) = A( I, J ) + X( I )*TEMP
   10          CONTINUE
            END IF
            JY = JY + INCY
   20    CONTINUE
      ELSE
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( M - 1 )*INCX
         END IF
         DO 40, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               IX   = KX
               DO 30, I = 1, M
                  A( I, J ) = A( I, J ) + X( IX )*TEMP
                  IX        = IX        + INCX
   30          CONTINUE
            END IF
            JY = JY + INCY
   40    CONTINUE
      END IF
!
      RETURN
!
!     End of ZGERU .
!
      END
      subroutine  zscal(n,za,zx,incx)
!
!     scales a vector by a constant.
!     jack dongarra, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      double complex za,zx(*)
      integer i,incx,ix,n
!
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
!
!        code for increment not equal to 1
!
      ix = 1
      do 10 i = 1,n
        zx(ix) = za*zx(ix)
        ix = ix + incx
   10 continue
      return
!
!        code for increment equal to 1
!
   20 do 30 i = 1,n
        zx(i) = za*zx(i)
   30 continue
      return
      end
      subroutine  zswap (n,zx,incx,zy,incy)
!
!     interchanges two vectors.
!     jack dongarra, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      double complex zx(*),zy(*),ztemp
      integer i,incx,incy,ix,iy,n
!
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
!
!       code for unequal increments or equal increments not equal
!         to 1
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        ztemp = zx(ix)
        zx(ix) = zy(iy)
        zy(iy) = ztemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!
!       code for both increments equal to 1
   20 do 30 i = 1,n
        ztemp = zx(i)
        zx(i) = zy(i)
        zy(i) = ztemp
   30 continue
      return
      end
      SUBROUTINE ZTRSM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, &
     &                   B, LDB )
!     .. Scalar Arguments ..
      CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, LDA, LDB
      COMPLEX*16         ALPHA
!     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  ZTRSM  solves one of the matrix equations
!
!     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
!
!  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
!  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
!
!     op( A ) = A   or   op( A ) = A'   or   op( A ) = conjg( A' ).
!
!  The matrix X is overwritten on B.
!
!  Parameters
!  ==========
!
!  SIDE   - CHARACTER*1.
!           On entry, SIDE specifies whether op( A ) appears on the left
!           or right of X as follows:
!
!              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
!
!              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
!
!           Unchanged on exit.
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the matrix A is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n'   op( A ) = A.
!
!              TRANSA = 'T' or 't'   op( A ) = A'.
!
!              TRANSA = 'C' or 'c'   op( A ) = conjg( A' ).
!
!           Unchanged on exit.
!
!  DIAG   - CHARACTER*1.
!           On entry, DIAG specifies whether or not A is unit triangular
!           as follows:
!
!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of B. M must be at
!           least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of B.  N must be
!           at least zero.
!           Unchanged on exit.
!
!  ALPHA  - COMPLEX*16      .
!           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
!           zero then  A is not referenced and  B need not be set before
!           entry.
!           Unchanged on exit.
!
!  A      - COMPLEX*16       array of DIMENSION ( LDA, k ), where k is m
!           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
!           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
!           upper triangular part of the array  A must contain the upper
!           triangular matrix  and the strictly lower triangular part of
!           A is not referenced.
!           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
!           lower triangular part of the array  A must contain the lower
!           triangular matrix  and the strictly upper triangular part of
!           A is not referenced.
!           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
!           A  are not referenced either,  but are assumed to be  unity.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
!           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
!           then LDA must be at least max( 1, n ).
!           Unchanged on exit.
!
!  B      - COMPLEX*16       array of DIMENSION ( LDB, n ).
!           Before entry,  the leading  m by n part of the array  B must
!           contain  the  right-hand  side  matrix  B,  and  on exit  is
!           overwritten by the solution matrix  X.
!
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in  the  calling  (sub)  program.   LDB  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 3 Blas routine.
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC          DCONJG, MAX
!     .. Local Scalars ..
      LOGICAL            LSIDE, NOCONJ, NOUNIT, UPPER
      INTEGER            I, INFO, J, K, NROWA
      COMPLEX*16         TEMP
!     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER        ( ONE  = ( 1.0D+0, 0.0D+0 ) )
      COMPLEX*16         ZERO
      PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) )
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      LSIDE  = LSAME( SIDE  , 'L' )
      IF( LSIDE )THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      NOCONJ = LSAME( TRANSA, 'T' )
      NOUNIT = LSAME( DIAG  , 'N' )
      UPPER  = LSAME( UPLO  , 'U' )
!
      INFO   = 0
      IF(      ( .NOT.LSIDE                ).AND.                       &
     &         ( .NOT.LSAME( SIDE  , 'R' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.UPPER                ).AND.                       &
     &         ( .NOT.LSAME( UPLO  , 'L' ) )      )THEN
         INFO = 2
      ELSE IF( ( .NOT.LSAME( TRANSA, 'N' ) ).AND.                       &
     &         ( .NOT.LSAME( TRANSA, 'T' ) ).AND.                       &
     &         ( .NOT.LSAME( TRANSA, 'C' ) )      )THEN
         INFO = 3
      ELSE IF( ( .NOT.LSAME( DIAG  , 'U' ) ).AND.                       &
     &         ( .NOT.LSAME( DIAG  , 'N' ) )      )THEN
         INFO = 4
      ELSE IF( M  .LT.0               )THEN
         INFO = 5
      ELSE IF( N  .LT.0               )THEN
         INFO = 6
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 9
      ELSE IF( LDB.LT.MAX( 1, M     ) )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'ZTRSM ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( N.EQ.0 )                                                      &
     &   RETURN
!
!     And when  alpha.eq.zero.
!
      IF( ALPHA.EQ.ZERO )THEN
         DO 20, J = 1, N
            DO 10, I = 1, M
               B( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
         RETURN
      END IF
!
!     Start the operations.
!
      IF( LSIDE )THEN
         IF( LSAME( TRANSA, 'N' ) )THEN
!
!           Form  B := alpha*inv( A )*B.
!
            IF( UPPER )THEN
               DO 60, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 30, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   30                CONTINUE
                  END IF
                  DO 50, K = M, 1, -1
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT )                                    &
     &                     B( K, J ) = B( K, J )/A( K, K )
                        DO 40, I = 1, K - 1
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   40                   CONTINUE
                     END IF
   50             CONTINUE
   60          CONTINUE
            ELSE
               DO 100, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 70, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   70                CONTINUE
                  END IF
                  DO 90 K = 1, M
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT )                                    &
     &                     B( K, J ) = B( K, J )/A( K, K )
                        DO 80, I = K + 1, M
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   80                   CONTINUE
                     END IF
   90             CONTINUE
  100          CONTINUE
            END IF
         ELSE
!
!           Form  B := alpha*inv( A' )*B
!           or    B := alpha*inv( conjg( A' ) )*B.
!
            IF( UPPER )THEN
               DO 140, J = 1, N
                  DO 130, I = 1, M
                     TEMP = ALPHA*B( I, J )
                     IF( NOCONJ )THEN
                        DO 110, K = 1, I - 1
                           TEMP = TEMP - A( K, I )*B( K, J )
  110                   CONTINUE
                        IF( NOUNIT )                                    &
     &                     TEMP = TEMP/A( I, I )
                     ELSE
                        DO 120, K = 1, I - 1
                           TEMP = TEMP - DCONJG( A( K, I ) )*B( K, J )
  120                   CONTINUE
                        IF( NOUNIT )                                    &
     &                     TEMP = TEMP/DCONJG( A( I, I ) )
                     END IF
                     B( I, J ) = TEMP
  130             CONTINUE
  140          CONTINUE
            ELSE
               DO 180, J = 1, N
                  DO 170, I = M, 1, -1
                     TEMP = ALPHA*B( I, J )
                     IF( NOCONJ )THEN
                        DO 150, K = I + 1, M
                           TEMP = TEMP - A( K, I )*B( K, J )
  150                   CONTINUE
                        IF( NOUNIT )                                    &
     &                     TEMP = TEMP/A( I, I )
                     ELSE
                        DO 160, K = I + 1, M
                           TEMP = TEMP - DCONJG( A( K, I ) )*B( K, J )
  160                   CONTINUE
                        IF( NOUNIT )                                    &
     &                     TEMP = TEMP/DCONJG( A( I, I ) )
                     END IF
                     B( I, J ) = TEMP
  170             CONTINUE
  180          CONTINUE
            END IF
         END IF
      ELSE
         IF( LSAME( TRANSA, 'N' ) )THEN
!
!           Form  B := alpha*B*inv( A ).
!
            IF( UPPER )THEN
               DO 230, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 190, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  190                CONTINUE
                  END IF
                  DO 210, K = 1, J - 1
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 200, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  200                   CONTINUE
                     END IF
  210             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 220, I = 1, M
                        B( I, J ) = TEMP*B( I, J )
  220                CONTINUE
                  END IF
  230          CONTINUE
            ELSE
               DO 280, J = N, 1, -1
                  IF( ALPHA.NE.ONE )THEN
                     DO 240, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  240                CONTINUE
                  END IF
                  DO 260, K = J + 1, N
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 250, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  250                   CONTINUE
                     END IF
  260             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 270, I = 1, M
                       B( I, J ) = TEMP*B( I, J )
  270                CONTINUE
                  END IF
  280          CONTINUE
            END IF
         ELSE
!
!           Form  B := alpha*B*inv( A' )
!           or    B := alpha*B*inv( conjg( A' ) ).
!
            IF( UPPER )THEN
               DO 330, K = N, 1, -1
                  IF( NOUNIT )THEN
                     IF( NOCONJ )THEN
                        TEMP = ONE/A( K, K )
                     ELSE
                        TEMP = ONE/DCONJG( A( K, K ) )
                     END IF
                     DO 290, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  290                CONTINUE
                  END IF
                  DO 310, J = 1, K - 1
                     IF( A( J, K ).NE.ZERO )THEN
                        IF( NOCONJ )THEN
                           TEMP = A( J, K )
                        ELSE
                           TEMP = DCONJG( A( J, K ) )
                        END IF
                        DO 300, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  300                   CONTINUE
                     END IF
  310             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 320, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  320                CONTINUE
                  END IF
  330          CONTINUE
            ELSE
               DO 380, K = 1, N
                  IF( NOUNIT )THEN
                     IF( NOCONJ )THEN
                        TEMP = ONE/A( K, K )
                     ELSE
                        TEMP = ONE/DCONJG( A( K, K ) )
                     END IF
                     DO 340, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  340                CONTINUE
                  END IF
                  DO 360, J = K + 1, N
                     IF( A( J, K ).NE.ZERO )THEN
                        IF( NOCONJ )THEN
                           TEMP = A( J, K )
                        ELSE
                           TEMP = DCONJG( A( J, K ) )
                        END IF
                        DO 350, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  350                   CONTINUE
                     END IF
  360             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 370, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  370                CONTINUE
                  END IF
  380          CONTINUE
            END IF
         END IF
      END IF
!
      RETURN
!
!     End of ZTRSM .
!
      END
      SUBROUTINE CLACON( N, V, X, EST, KASE )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      INTEGER            KASE, N
      REAL               EST
!     ..
!     .. Array Arguments ..
      COMPLEX            V( N ), X( N )
!     ..
!
!  Purpose
!  =======
!
!  CLACON estimates the 1-norm of a square, complex matrix A.
!  Reverse communication is used for evaluating matrix-vector products.
!
!  Arguments
!  =========
!
!  N      (input) INTEGER
!         The order of the matrix.  N >= 1.
!
!  V      (workspace) COMPLEX array, dimension (N)
!         On the final return, V = A*W,  where  EST = norm(V)/norm(W)
!         (W is not returned).
!
!  X      (input/output) COMPLEX array, dimension (N)
!         On an intermediate return, X should be overwritten by
!               A * X,   if KASE=1,
!               A' * X,  if KASE=2,
!         where A' is the conjugate transpose of A, and CLACON must be
!         re-called with all the other parameters unchanged.
!
!  EST    (output) REAL
!         An estimate (a lower bound) for norm(A).
!
!  KASE   (input/output) INTEGER
!         On the initial call to CLACON, KASE should be 0.
!         On an intermediate return, KASE will be 1 or 2, indicating
!         whether X should be overwritten by A * X  or A' * X.
!         On the final return from CLACON, KASE will again be 0.
!
!  Further Details
!  ======= =======
!
!  Contributed by Nick Higham, University of Manchester.
!  Originally named CONEST, dated March 16, 1988.
!
!  Reference: N.J. Higham, "FORTRAN codes for estimating the one-norm of
!  a real or complex matrix, with applications to condition estimation",
!  ACM Trans. Math. Soft., vol. 14, no. 4, pp. 381-396, December 1988.
!
!  Last modified:  April, 1999
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER            ITMAX
      PARAMETER          ( ITMAX = 5 )
      REAL               ONE, TWO
      PARAMETER          ( ONE = 1.0E0, TWO = 2.0E0 )
      COMPLEX            CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0E0, 0.0E0 ),                    &
     &                   CONE = ( 1.0E0, 0.0E0 ) )
!     ..
!     .. Local Scalars ..
      INTEGER            I, ITER, J, JLAST, JUMP
      REAL               ABSXI, ALTSGN, ESTOLD, SAFMIN, TEMP
!     ..
!     .. External Functions ..
      INTEGER            ICMAX1
      REAL               SCSUM1, SLAMCH
      EXTERNAL           ICMAX1, SCSUM1, SLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL           CCOPY
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, AIMAG, CMPLX, REAL
!     ..
!     .. Save statement ..
      SAVE
!     ..
!     .. Executable Statements ..
!
      SAFMIN = SLAMCH( 'Safe minimum' )
      IF( KASE.EQ.0 ) THEN
         DO 10 I = 1, N
            X( I ) = CMPLX( ONE / REAL( N ) )
   10    CONTINUE
         KASE = 1
         JUMP = 1
         RETURN
      END IF
!
      GO TO ( 20, 40, 70, 90, 120 )JUMP
!
!     ................ ENTRY   (JUMP = 1)
!     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X.
!
   20 CONTINUE
      IF( N.EQ.1 ) THEN
         V( 1 ) = X( 1 )
         EST = ABS( V( 1 ) )
!        ... QUIT
         GO TO 130
      END IF
      EST = SCSUM1( N, X, 1 )
!
      DO 30 I = 1, N
         ABSXI = ABS( X( I ) )
         IF( ABSXI.GT.SAFMIN ) THEN
            X( I ) = CMPLX( REAL( X( I ) ) / ABSXI,                     &
     &               AIMAG( X( I ) ) / ABSXI )
         ELSE
            X( I ) = CONE
         END IF
   30 CONTINUE
      KASE = 2
      JUMP = 2
      RETURN
!
!     ................ ENTRY   (JUMP = 2)
!     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY CTRANS(A)*X.
!
   40 CONTINUE
      J = ICMAX1( N, X, 1 )
      ITER = 2
!
!     MAIN LOOP - ITERATIONS 2,3,...,ITMAX.
!
   50 CONTINUE
      DO 60 I = 1, N
         X( I ) = CZERO
   60 CONTINUE
      X( J ) = CONE
      KASE = 1
      JUMP = 3
      RETURN
!
!     ................ ENTRY   (JUMP = 3)
!     X HAS BEEN OVERWRITTEN BY A*X.
!
   70 CONTINUE
      CALL CCOPY( N, X, 1, V, 1 )
      ESTOLD = EST
      EST = SCSUM1( N, V, 1 )
!
!     TEST FOR CYCLING.
      IF( EST.LE.ESTOLD )                                               &
     &   GO TO 100
!
      DO 80 I = 1, N
         ABSXI = ABS( X( I ) )
         IF( ABSXI.GT.SAFMIN ) THEN
            X( I ) = CMPLX( REAL( X( I ) ) / ABSXI,                     &
     &               AIMAG( X( I ) ) / ABSXI )
         ELSE
            X( I ) = CONE
         END IF
   80 CONTINUE
      KASE = 2
      JUMP = 4
      RETURN
!
!     ................ ENTRY   (JUMP = 4)
!     X HAS BEEN OVERWRITTEN BY CTRANS(A)*X.
!
   90 CONTINUE
      JLAST = J
      J = ICMAX1( N, X, 1 )
      IF( ( ABS( X( JLAST ) ).NE.ABS( X( J ) ) ) .AND.                  &
     &    ( ITER.LT.ITMAX ) ) THEN
         ITER = ITER + 1
         GO TO 50
      END IF
!
!     ITERATION COMPLETE.  FINAL STAGE.
!
  100 CONTINUE
      ALTSGN = ONE
      DO 110 I = 1, N
         X( I ) = CMPLX( ALTSGN*( ONE+REAL( I-1 ) / REAL( N-1 ) ) )
         ALTSGN = -ALTSGN
  110 CONTINUE
      KASE = 1
      JUMP = 5
      RETURN
!
!     ................ ENTRY   (JUMP = 5)
!     X HAS BEEN OVERWRITTEN BY A*X.
!
  120 CONTINUE
      TEMP = TWO*( SCSUM1( N, X, 1 ) / REAL( 3*N ) )
      IF( TEMP.GT.EST ) THEN
         CALL CCOPY( N, X, 1, V, 1 )
         EST = TEMP
      END IF
!
  130 CONTINUE
      KASE = 0
      RETURN
!
!     End of CLACON
!
      END
      SUBROUTINE ZGBTF2( M, N, KL, KU, AB, LDAB, IPIV, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      INTEGER            INFO, KL, KU, LDAB, M, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         AB( LDAB, * )
!     ..
!
!  Purpose
!  =======
!
!  ZGBTF2 computes an LU factorization of a complex m-by-n band matrix
!  A using partial pivoting with row interchanges.
!
!  This is the unblocked version of the algorithm, calling Level 2 BLAS.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  KL      (input) INTEGER
!          The number of subdiagonals within the band of A.  KL >= 0.
!
!  KU      (input) INTEGER
!          The number of superdiagonals within the band of A.  KU >= 0.
!
!  AB      (input/output) COMPLEX*16 array, dimension (LDAB,N)
!          On entry, the matrix A in band storage, in rows KL+1 to
!          2*KL+KU+1; rows 1 to KL of the array need not be set.
!          The j-th column of A is stored in the j-th column of the
!          array AB as follows:
!          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
!
!          On exit, details of the factorization: U is stored as an
!          upper triangular band matrix with KL+KU superdiagonals in
!          rows 1 to KL+KU+1, and the multipliers used during the
!          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
!          See below for further details.
!
!  LDAB    (input) INTEGER
!          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!          > 0: if INFO = +i, U(i,i) is exactly zero. The factorization
!               has been completed, but the factor U is exactly
!               singular, and division by zero will occur if it is used
!               to solve a system of equations.
!
!  Further Details
!  ===============
!
!  The band storage scheme is illustrated by the following example, when
!  M = N = 6, KL = 2, KU = 1:
!
!  On entry:                       On exit:
!
!      *    *    *    +    +    +       *    *    *   u14  u25  u36
!      *    *    +    +    +    +       *    *   u13  u24  u35  u46
!      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
!     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
!     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
!     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *
!
!  Array elements marked * are not used by the routine; elements marked
!  + need not be set on entry, but are required by the routine to store
!  elements of U, because of fill-in resulting from the row
!  interchanges.
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ),                    &
     &                   ZERO = ( 0.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      INTEGER            I, J, JP, JU, KM, KV
!     ..
!     .. External Functions ..
      INTEGER            IZAMAX
      EXTERNAL           IZAMAX
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZGERU, ZSCAL, ZSWAP
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     KV is the number of superdiagonals in the factor U, allowing for
!     fill-in.
!
      KV = KU + KL
!
!     Test the input parameters.
!
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KL.LT.0 ) THEN
         INFO = -3
      ELSE IF( KU.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDAB.LT.KL+KV+1 ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGBTF2', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 )                                          &
     &   RETURN
!
!     Gaussian elimination with partial pivoting
!
!     Set fill-in elements in columns KU+2 to KV to zero.
!
      DO 20 J = KU + 2, MIN( KV, N )
         DO 10 I = KV - J + 2, KL
            AB( I, J ) = ZERO
   10    CONTINUE
   20 CONTINUE
!
!     JU is the index of the last column affected by the current stage
!     of the factorization.
!
      JU = 1
!
      DO 40 J = 1, MIN( M, N )
!
!        Set fill-in elements in column J+KV to zero.
!
         IF( J+KV.LE.N ) THEN
            DO 30 I = 1, KL
               AB( I, J+KV ) = ZERO
   30       CONTINUE
         END IF
!
!        Find pivot and test for singularity. KM is the number of
!        subdiagonal elements in the current column.
!
         KM = MIN( KL, M-J )
         JP = IZAMAX( KM+1, AB( KV+1, J ), 1 )
         IPIV( J ) = JP + J - 1
         IF( AB( KV+JP, J ).NE.ZERO ) THEN
            JU = MAX( JU, MIN( J+KU+JP-1, N ) )
!
!           Apply interchange to columns J to JU.
!
            IF( JP.NE.1 )                                               &
     &         CALL ZSWAP( JU-J+1, AB( KV+JP, J ), LDAB-1,              &
     &                     AB( KV+1, J ), LDAB-1 )
            IF( KM.GT.0 ) THEN
!
!              Compute multipliers.
!
               CALL ZSCAL( KM, ONE / AB( KV+1, J ), AB( KV+2, J ), 1 )
!
!              Update trailing submatrix within the band.
!
               IF( JU.GT.J )                                            &
     &            CALL ZGERU( KM, JU-J, -ONE, AB( KV+2, J ), 1,         &
     &                        AB( KV, J+1 ), LDAB-1, AB( KV+1, J+1 ),   &
     &                        LDAB-1 )
            END IF
         ELSE
!
!           If pivot is zero, set INFO to the index of the pivot
!           unless a zero pivot has already been found.
!
            IF( INFO.EQ.0 )                                             &
     &         INFO = J
         END IF
   40 CONTINUE
      RETURN
!
!     End of ZGBTF2
!
      END
      SUBROUTINE ZGBTRF( M, N, KL, KU, AB, LDAB, IPIV, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      INTEGER            INFO, KL, KU, LDAB, M, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         AB( LDAB, * )
!     ..
!
!  Purpose
!  =======
!
!  ZGBTRF computes an LU factorization of a complex m-by-n band matrix A
!  using partial pivoting with row interchanges.
!
!  This is the blocked version of the algorithm, calling Level 3 BLAS.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  KL      (input) INTEGER
!          The number of subdiagonals within the band of A.  KL >= 0.
!
!  KU      (input) INTEGER
!          The number of superdiagonals within the band of A.  KU >= 0.
!
!  AB      (input/output) COMPLEX*16 array, dimension (LDAB,N)
!          On entry, the matrix A in band storage, in rows KL+1 to
!          2*KL+KU+1; rows 1 to KL of the array need not be set.
!          The j-th column of A is stored in the j-th column of the
!          array AB as follows:
!          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
!
!          On exit, details of the factorization: U is stored as an
!          upper triangular band matrix with KL+KU superdiagonals in
!          rows 1 to KL+KU+1, and the multipliers used during the
!          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
!          See below for further details.
!
!  LDAB    (input) INTEGER
!          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!          > 0: if INFO = +i, U(i,i) is exactly zero. The factorization
!               has been completed, but the factor U is exactly
!               singular, and division by zero will occur if it is used
!               to solve a system of equations.
!
!  Further Details
!  ===============
!
!  The band storage scheme is illustrated by the following example, when
!  M = N = 6, KL = 2, KU = 1:
!
!  On entry:                       On exit:
!
!      *    *    *    +    +    +       *    *    *   u14  u25  u36
!      *    *    +    +    +    +       *    *   u13  u24  u35  u46
!      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
!     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
!     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
!     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *
!
!  Array elements marked * are not used by the routine; elements marked
!  + need not be set on entry, but are required by the routine to store
!  elements of U because of fill-in resulting from the row interchanges.
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ),                    &
     &                   ZERO = ( 0.0D+0, 0.0D+0 ) )
      INTEGER            NBMAX, LDWORK
      PARAMETER          ( NBMAX = 64, LDWORK = NBMAX+1 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, I2, I3, II, IP, J, J2, J3, JB, JJ, JM, JP,  &
     &                   JU, K2, KM, KV, NB, NW
      COMPLEX*16         TEMP
!     ..
!     .. Local Arrays ..
      COMPLEX*16         WORK13( LDWORK, NBMAX ),                       &
     &                   WORK31( LDWORK, NBMAX )
!     ..
!     .. External Functions ..
      INTEGER            ILAENV, IZAMAX
      EXTERNAL           ILAENV, IZAMAX
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZCOPY, ZGBTF2, ZGEMM, ZGERU, ZLASWP,   &
     &                   ZSCAL, ZSWAP, ZTRSM
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     KV is the number of superdiagonals in the factor U, allowing for
!     fill-in
!
      KV = KU + KL
!
!     Test the input parameters.
!
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KL.LT.0 ) THEN
         INFO = -3
      ELSE IF( KU.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDAB.LT.KL+KV+1 ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGBTRF', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 )                                          &
     &   RETURN
!
!     Determine the block size for this environment
!
      NB = ILAENV( 1, 'ZGBTRF', ' ', M, N, KL, KU )
!
!     The block size must not exceed the limit set by the size of the
!     local arrays WORK13 and WORK31.
!
      NB = MIN( NB, NBMAX )
!
      IF( NB.LE.1 .OR. NB.GT.KL ) THEN
!
!        Use unblocked code
!
         CALL ZGBTF2( M, N, KL, KU, AB, LDAB, IPIV, INFO )
      ELSE
!
!        Use blocked code
!
!        Zero the superdiagonal elements of the work array WORK13
!
         DO 20 J = 1, NB
            DO 10 I = 1, J - 1
               WORK13( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
!
!        Zero the subdiagonal elements of the work array WORK31
!
         DO 40 J = 1, NB
            DO 30 I = J + 1, NB
               WORK31( I, J ) = ZERO
   30       CONTINUE
   40    CONTINUE
!
!        Gaussian elimination with partial pivoting
!
!        Set fill-in elements in columns KU+2 to KV to zero
!
         DO 60 J = KU + 2, MIN( KV, N )
            DO 50 I = KV - J + 2, KL
               AB( I, J ) = ZERO
   50       CONTINUE
   60    CONTINUE
!
!        JU is the index of the last column affected by the current
!        stage of the factorization
!
         JU = 1
!
         DO 180 J = 1, MIN( M, N ), NB
            JB = MIN( NB, MIN( M, N )-J+1 )
!
!           The active part of the matrix is partitioned
!
!              A11   A12   A13
!              A21   A22   A23
!              A31   A32   A33
!
!           Here A11, A21 and A31 denote the current block of JB columns
!           which is about to be factorized. The number of rows in the
!           partitioning are JB, I2, I3 respectively, and the numbers
!           of columns are JB, J2, J3. The superdiagonal elements of A13
!           and the subdiagonal elements of A31 lie outside the band.
!
            I2 = MIN( KL-JB, M-J-JB+1 )
            I3 = MIN( JB, M-J-KL+1 )
!
!           J2 and J3 are computed after JU has been updated.
!
!           Factorize the current block of JB columns
!
            DO 80 JJ = J, J + JB - 1
!
!              Set fill-in elements in column JJ+KV to zero
!
               IF( JJ+KV.LE.N ) THEN
                  DO 70 I = 1, KL
                     AB( I, JJ+KV ) = ZERO
   70             CONTINUE
               END IF
!
!              Find pivot and test for singularity. KM is the number of
!              subdiagonal elements in the current column.
!
               KM = MIN( KL, M-JJ )
               JP = IZAMAX( KM+1, AB( KV+1, JJ ), 1 )
               IPIV( JJ ) = JP + JJ - J
               IF( AB( KV+JP, JJ ).NE.ZERO ) THEN
                  JU = MAX( JU, MIN( JJ+KU+JP-1, N ) )
                  IF( JP.NE.1 ) THEN
!
!                    Apply interchange to columns J to J+JB-1
!
                     IF( JP+JJ-1.LT.J+KL ) THEN
!
                        CALL ZSWAP( JB, AB( KV+1+JJ-J, J ), LDAB-1,     &
     &                              AB( KV+JP+JJ-J, J ), LDAB-1 )
                     ELSE
!
!                       The interchange affects columns J to JJ-1 of A31
!                       which are stored in the work array WORK31
!
                        CALL ZSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1,   &
     &                              WORK31( JP+JJ-J-KL, 1 ), LDWORK )
                        CALL ZSWAP( J+JB-JJ, AB( KV+1, JJ ), LDAB-1,    &
     &                              AB( KV+JP, JJ ), LDAB-1 )
                     END IF
                  END IF
!
!                 Compute multipliers
!
                  CALL ZSCAL( KM, ONE / AB( KV+1, JJ ), AB( KV+2, JJ ), &
     &                        1 )
!
!                 Update trailing submatrix within the band and within
!                 the current block. JM is the index of the last column
!                 which needs to be updated.
!
                  JM = MIN( JU, J+JB-1 )
                  IF( JM.GT.JJ )                                        &
     &               CALL ZGERU( KM, JM-JJ, -ONE, AB( KV+2, JJ ), 1,    &
     &                           AB( KV, JJ+1 ), LDAB-1,                &
     &                           AB( KV+1, JJ+1 ), LDAB-1 )
               ELSE
!
!                 If pivot is zero, set INFO to the index of the pivot
!                 unless a zero pivot has already been found.
!
                  IF( INFO.EQ.0 )                                       &
     &               INFO = JJ
               END IF
!
!              Copy current column of A31 into the work array WORK31
!
               NW = MIN( JJ-J+1, I3 )
               IF( NW.GT.0 )                                            &
     &            CALL ZCOPY( NW, AB( KV+KL+1-JJ+J, JJ ), 1,            &
     &                        WORK31( 1, JJ-J+1 ), 1 )
   80       CONTINUE
            IF( J+JB.LE.N ) THEN
!
!              Apply the row interchanges to the other blocks.
!
               J2 = MIN( JU-J+1, KV ) - JB
               J3 = MAX( 0, JU-J-KV+1 )
!
!              Use ZLASWP to apply the row interchanges to A12, A22, and
!              A32.
!
               CALL ZLASWP( J2, AB( KV+1-JB, J+JB ), LDAB-1, 1, JB,     &
     &                      IPIV( J ), 1 )
!
!              Adjust the pivot indices.
!
               DO 90 I = J, J + JB - 1
                  IPIV( I ) = IPIV( I ) + J - 1
   90          CONTINUE
!
!              Apply the row interchanges to A13, A23, and A33
!              columnwise.
!
               K2 = J - 1 + JB + J2
               DO 110 I = 1, J3
                  JJ = K2 + I
                  DO 100 II = J + I - 1, J + JB - 1
                     IP = IPIV( II )
                     IF( IP.NE.II ) THEN
                        TEMP = AB( KV+1+II-JJ, JJ )
                        AB( KV+1+II-JJ, JJ ) = AB( KV+1+IP-JJ, JJ )
                        AB( KV+1+IP-JJ, JJ ) = TEMP
                     END IF
  100             CONTINUE
  110          CONTINUE
!
!              Update the relevant part of the trailing submatrix
!
               IF( J2.GT.0 ) THEN
!
!                 Update A12
!
                  CALL ZTRSM( 'Left', 'Lower', 'No transpose', 'Unit',  &
     &                        JB, J2, ONE, AB( KV+1, J ), LDAB-1,       &
     &                        AB( KV+1-JB, J+JB ), LDAB-1 )
!
                  IF( I2.GT.0 ) THEN
!
!                    Update A22
!
                     CALL ZGEMM( 'No transpose', 'No transpose', I2, J2,&
     &                           JB, -ONE, AB( KV+1+JB, J ), LDAB-1,    &
     &                           AB( KV+1-JB, J+JB ), LDAB-1, ONE,      &
     &                           AB( KV+1, J+JB ), LDAB-1 )
                  END IF
!
                  IF( I3.GT.0 ) THEN
!
!                    Update A32
!
                     CALL ZGEMM( 'No transpose', 'No transpose', I3, J2,&
     &                           JB, -ONE, WORK31, LDWORK,              &
     &                           AB( KV+1-JB, J+JB ), LDAB-1, ONE,      &
     &                           AB( KV+KL+1-JB, J+JB ), LDAB-1 )
                  END IF
               END IF
!
               IF( J3.GT.0 ) THEN
!
!                 Copy the lower triangle of A13 into the work array
!                 WORK13
!
                  DO 130 JJ = 1, J3
                     DO 120 II = JJ, JB
                        WORK13( II, JJ ) = AB( II-JJ+1, JJ+J+KV-1 )
  120                CONTINUE
  130             CONTINUE
!
!                 Update A13 in the work array
!
                  CALL ZTRSM( 'Left', 'Lower', 'No transpose', 'Unit',  &
     &                        JB, J3, ONE, AB( KV+1, J ), LDAB-1,       &
     &                        WORK13, LDWORK )
!
                  IF( I2.GT.0 ) THEN
!
!                    Update A23
!
                     CALL ZGEMM( 'No transpose', 'No transpose', I2, J3,&
     &                           JB, -ONE, AB( KV+1+JB, J ), LDAB-1,    &
     &                           WORK13, LDWORK, ONE, AB( 1+JB, J+KV ), &
     &                           LDAB-1 )
                  END IF
!
                  IF( I3.GT.0 ) THEN
!
!                    Update A33
!
                     CALL ZGEMM( 'No transpose', 'No transpose', I3, J3,&
     &                           JB, -ONE, WORK31, LDWORK, WORK13,      &
     &                           LDWORK, ONE, AB( 1+KL, J+KV ), LDAB-1 )
                  END IF
!
!                 Copy the lower triangle of A13 back into place
!
                  DO 150 JJ = 1, J3
                     DO 140 II = JJ, JB
                        AB( II-JJ+1, JJ+J+KV-1 ) = WORK13( II, JJ )
  140                CONTINUE
  150             CONTINUE
               END IF
            ELSE
!
!              Adjust the pivot indices.
!
               DO 160 I = J, J + JB - 1
                  IPIV( I ) = IPIV( I ) + J - 1
  160          CONTINUE
            END IF
!
!           Partially undo the interchanges in the current block to
!           restore the upper triangular form of A31 and copy the upper
!           triangle of A31 back into place
!
            DO 170 JJ = J + JB - 1, J, -1
               JP = IPIV( JJ ) - JJ + 1
               IF( JP.NE.1 ) THEN
!
!                 Apply interchange to columns J to JJ-1
!
                  IF( JP+JJ-1.LT.J+KL ) THEN
!
!                    The interchange does not affect A31
!
                     CALL ZSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1,      &
     &                           AB( KV+JP+JJ-J, J ), LDAB-1 )
                  ELSE
!
!                    The interchange does affect A31
!
                     CALL ZSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1,      &
     &                           WORK31( JP+JJ-J-KL, 1 ), LDWORK )
                  END IF
               END IF
!
!              Copy the current column of A31 back into place
!
               NW = MIN( I3, JJ-J+1 )
               IF( NW.GT.0 )                                            &
     &            CALL ZCOPY( NW, WORK31( 1, JJ-J+1 ), 1,               &
     &                        AB( KV+KL+1-JJ+J, JJ ), 1 )
  170       CONTINUE
  180    CONTINUE
      END IF
!
      RETURN
!
!     End of ZGBTRF
!
      END
      SUBROUTINE ZGETF2( M, N, A, LDA, IPIV, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  ZGETF2 computes an LU factorization of a general m-by-n matrix A
!  using partial pivoting with row interchanges.
!
!  The factorization has the form
!     A = P * L * U
!  where P is a permutation matrix, L is lower triangular with unit
!  diagonal elements (lower trapezoidal if m > n), and U is upper
!  triangular (upper trapezoidal if m < n).
!
!  This is the right-looking Level 2 BLAS version of the algorithm.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
!          On entry, the m by n matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value
!          > 0: if INFO = k, U(k,k) is exactly zero. The factorization
!               has been completed, but the factor U is exactly
!               singular, and division by zero will occur if it is used
!               to solve a system of equations.
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ),                    &
     &                   ZERO = ( 0.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      INTEGER            J, JP
!     ..
!     .. External Functions ..
      INTEGER            IZAMAX
      EXTERNAL           IZAMAX
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZGERU, ZSCAL, ZSWAP
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGETF2', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 )                                          &
     &   RETURN
!
      DO 10 J = 1, MIN( M, N )
!
!        Find pivot and test for singularity.
!
         JP = J - 1 + IZAMAX( M-J+1, A( J, J ), 1 )
         IPIV( J ) = JP
         IF( A( JP, J ).NE.ZERO ) THEN
!
!           Apply the interchange to columns 1:N.
!
            IF( JP.NE.J )                                               &
     &         CALL ZSWAP( N, A( J, 1 ), LDA, A( JP, 1 ), LDA )
!
!           Compute elements J+1:M of J-th column.
!
            IF( J.LT.M )                                                &
     &         CALL ZSCAL( M-J, ONE / A( J, J ), A( J+1, J ), 1 )
!
         ELSE IF( INFO.EQ.0 ) THEN
!
            INFO = J
         END IF
!
         IF( J.LT.MIN( M, N ) ) THEN
!
!           Update trailing submatrix.
!
            CALL ZGERU( M-J, N-J, -ONE, A( J+1, J ), 1, A( J, J+1 ),    &
     &                  LDA, A( J+1, J+1 ), LDA )
         END IF
   10 CONTINUE
      RETURN
!
!     End of ZGETF2
!
      END
      SUBROUTINE ZGETRF( M, N, A, LDA, IPIV, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  ZGETRF computes an LU factorization of a general M-by-N matrix A
!  using partial pivoting with row interchanges.
!
!  The factorization has the form
!     A = P * L * U
!  where P is a permutation matrix, L is lower triangular with unit
!  diagonal elements (lower trapezoidal if m > n), and U is upper
!  triangular (upper trapezoidal if m < n).
!
!  This is the right-looking Level 3 BLAS version of the algorithm.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
!          On entry, the M-by-N matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
!                has been completed, but the factor U is exactly
!                singular, and division by zero will occur if it is used
!                to solve a system of equations.
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      INTEGER            I, IINFO, J, JB, NB
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZGEMM, ZGETF2, ZLASWP, ZTRSM
!     ..
!     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGETRF', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 )                                          &
     &   RETURN
!
!     Determine the block size for this environment.
!
      NB = ILAENV( 1, 'ZGETRF', ' ', M, N, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) THEN
!
!        Use unblocked code.
!
         CALL ZGETF2( M, N, A, LDA, IPIV, INFO )
      ELSE
!
!        Use blocked code.
!
         DO 20 J = 1, MIN( M, N ), NB
            JB = MIN( MIN( M, N )-J+1, NB )
!
!           Factor diagonal and subdiagonal blocks and test for exact
!           singularity.
!
            CALL ZGETF2( M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO )
!
!           Adjust INFO and the pivot indices.
!
            IF( INFO.EQ.0 .AND. IINFO.GT.0 )                            &
     &         INFO = IINFO + J - 1
            DO 10 I = J, MIN( M, J+JB-1 )
               IPIV( I ) = J - 1 + IPIV( I )
   10       CONTINUE
!
!           Apply interchanges to columns 1:J-1.
!
            CALL ZLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )
!
            IF( J+JB.LE.N ) THEN
!
!              Apply interchanges to columns J+JB:N.
!
               CALL ZLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1,     &
     &                      IPIV, 1 )
!
!              Compute block row of U.
!
               CALL ZTRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB, &
     &                     N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ), &
     &                     LDA )
               IF( J+JB.LE.M ) THEN
!
!                 Update trailing submatrix.
!
                  CALL ZGEMM( 'No transpose', 'No transpose', M-J-JB+1, &
     &                        N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA,    &
     &                        A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ),  &
     &                        LDA )
               END IF
            END IF
   20    CONTINUE
      END IF
      RETURN
!
!     End of ZGETRF
!
      END
      SUBROUTINE ZGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, LDA, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  ZGETRS solves a system of linear equations
!     A * X = B,  A**T * X = B,  or  A**H * X = B
!  with a general N-by-N matrix A using the LU factorization computed
!  by ZGETRF.
!
!  Arguments
!  =========
!
!  TRANS   (input) CHARACTER*1
!          Specifies the form of the system of equations:
!          = 'N':  A * X = B     (No transpose)
!          = 'T':  A**T * X = B  (Transpose)
!          = 'C':  A**H * X = B  (Conjugate transpose)
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  A       (input) COMPLEX*16 array, dimension (LDA,N)
!          The factors L and U from the factorization A = P*L*U
!          as computed by ZGETRF.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  IPIV    (input) INTEGER array, dimension (N)
!          The pivot indices from ZGETRF; for 1<=i<=N, row i of the
!          matrix was interchanged with row IPIV(i).
!
!  B       (input/output) COMPLEX*16 array, dimension (LDB,NRHS)
!          On entry, the right hand side matrix B.
!          On exit, the solution matrix X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      LOGICAL            NOTRAN
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLASWP, ZTRSM
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.        &
     &    LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGETRS', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 .OR. NRHS.EQ.0 )                                       &
     &   RETURN
!
      IF( NOTRAN ) THEN
!
!        Solve A * X = B.
!
!        Apply row interchanges to the right hand sides.
!
         CALL ZLASWP( NRHS, B, LDB, 1, N, IPIV, 1 )
!
!        Solve L*X = B, overwriting B with X.
!
         CALL ZTRSM( 'Left', 'Lower', 'No transpose', 'Unit', N, NRHS,  &
     &               ONE, A, LDA, B, LDB )
!
!        Solve U*X = B, overwriting B with X.
!
         CALL ZTRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N,    &
     &               NRHS, ONE, A, LDA, B, LDB )
      ELSE
!
!        Solve A**T * X = B  or A**H * X = B.
!
!        Solve U'*X = B, overwriting B with X.
!
         CALL ZTRSM( 'Left', 'Upper', TRANS, 'Non-unit', N, NRHS, ONE,  &
     &               A, LDA, B, LDB )
!
!        Solve L'*X = B, overwriting B with X.
!
         CALL ZTRSM( 'Left', 'Lower', TRANS, 'Unit', N, NRHS, ONE, A,   &
     &               LDA, B, LDB )
!
!        Apply row interchanges to the solution vectors.
!
         CALL ZLASWP( NRHS, B, LDB, 1, N, IPIV, -1 )
      END IF
!
      RETURN
!
!     End of ZGETRS
!
      END
      SUBROUTINE ZLASWP( N, A, LDA, K1, K2, IPIV, INCX )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      INTEGER            INCX, K1, K2, LDA, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  ZLASWP performs a series of row interchanges on the matrix A.
!  One row interchange is initiated for each of rows K1 through K2 of A.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.
!
!  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
!          On entry, the matrix of column dimension N to which the row
!          interchanges will be applied.
!          On exit, the permuted matrix.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.
!
!  K1      (input) INTEGER
!          The first element of IPIV for which a row interchange will
!          be done.
!
!  K2      (input) INTEGER
!          The last element of IPIV for which a row interchange will
!          be done.
!
!  IPIV    (input) INTEGER array, dimension (M*abs(INCX))
!          The vector of pivot indices.  Only the elements in positions
!          K1 through K2 of IPIV are accessed.
!          IPIV(K) = L implies rows K and L are to be interchanged.
!
!  INCX    (input) INTEGER
!          The increment between successive values of IPIV.  If IPIV
!          is negative, the pivots are applied in reverse order.
!
!  Further Details
!  ===============
!
!  Modified by
!   R. C. Whaley, Computer Science Dept., Univ. of Tenn., Knoxville, USA
!
! =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I, I1, I2, INC, IP, IX, IX0, J, K, N32
      COMPLEX*16         TEMP
!     ..
!     .. Executable Statements ..
!
!     Interchange row I with row IPIV(I) for each of rows K1 through K2.
!
      IF( INCX.GT.0 ) THEN
         IX0 = K1
         I1 = K1
         I2 = K2
         INC = 1
      ELSE IF( INCX.LT.0 ) THEN
         IX0 = 1 + ( 1-K2 )*INCX
         I1 = K2
         I2 = K1
         INC = -1
      ELSE
         RETURN
      END IF
!
      N32 = ( N / 32 )*32
      IF( N32.NE.0 ) THEN
         DO 30 J = 1, N32, 32
            IX = IX0
            DO 20 I = I1, I2, INC
               IP = IPIV( IX )
               IF( IP.NE.I ) THEN
                  DO 10 K = J, J + 31
                     TEMP = A( I, K )
                     A( I, K ) = A( IP, K )
                     A( IP, K ) = TEMP
   10             CONTINUE
               END IF
               IX = IX + INCX
   20       CONTINUE
   30    CONTINUE
      END IF
      IF( N32.NE.N ) THEN
         N32 = N32 + 1
         IX = IX0
         DO 50 I = I1, I2, INC
            IP = IPIV( IX )
            IF( IP.NE.I ) THEN
               DO 40 K = N32, N
                  TEMP = A( I, K )
                  A( I, K ) = A( IP, K )
                  A( IP, K ) = TEMP
   40          CONTINUE
            END IF
            IX = IX + INCX
   50    CONTINUE
      END IF
!
      RETURN
!
!     End of ZLASWP
!
      END
      SUBROUTINE DGETF2( M, N, A, LDA, IPIV, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1992
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  DGETF2 computes an LU factorization of a general m-by-n matrix A
!  using partial pivoting with row interchanges.
!
!  The factorization has the form
!     A = P * L * U
!  where P is a permutation matrix, L is lower triangular with unit
!  diagonal elements (lower trapezoidal if m > n), and U is upper
!  triangular (upper trapezoidal if m < n).
!
!  This is the right-looking Level 2 BLAS version of the algorithm.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the m by n matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value
!          > 0: if INFO = k, U(k,k) is exactly zero. The factorization
!               has been completed, but the factor U is exactly
!               singular, and division by zero will occur if it is used
!               to solve a system of equations.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            J, JP
!     ..
!     .. External Functions ..
      INTEGER            IDAMAX
      EXTERNAL           IDAMAX
!     ..
!     .. External Subroutines ..
      EXTERNAL           DGER, DSCAL, DSWAP, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETF2', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 )                                          &
     &   RETURN
!
      DO 10 J = 1, MIN( M, N )
!
!        Find pivot and test for singularity.
!
         JP = J - 1 + IDAMAX( M-J+1, A( J, J ), 1 )
         IPIV( J ) = JP
         IF( A( JP, J ).NE.ZERO ) THEN
!
!           Apply the interchange to columns 1:N.
!
            IF( JP.NE.J )                                               &
     &         CALL DSWAP( N, A( J, 1 ), LDA, A( JP, 1 ), LDA )
!
!           Compute elements J+1:M of J-th column.
!
            IF( J.LT.M )                                                &
     &         CALL DSCAL( M-J, ONE / A( J, J ), A( J+1, J ), 1 )
!
         ELSE IF( INFO.EQ.0 ) THEN
!
            INFO = J
         END IF
!
         IF( J.LT.MIN( M, N ) ) THEN
!
!           Update trailing submatrix.
!
            CALL DGER( M-J, N-J, -ONE, A( J+1, J ), 1, A( J, J+1 ), LDA,&
     &                 A( J+1, J+1 ), LDA )
         END IF
   10 CONTINUE
      RETURN
!
!     End of DGETF2
!
      END
      SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  DGETRF computes an LU factorization of a general M-by-N matrix A
!  using partial pivoting with row interchanges.
!
!  The factorization has the form
!     A = P * L * U
!  where P is a permutation matrix, L is lower triangular with unit
!  diagonal elements (lower trapezoidal if m > n), and U is upper
!  triangular (upper trapezoidal if m < n).
!
!  This is the right-looking Level 3 BLAS version of the algorithm.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the M-by-N matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
!                has been completed, but the factor U is exactly
!                singular, and division by zero will occur if it is used
!                to solve a system of equations.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, IINFO, J, JB, NB
!     ..
!     .. External Subroutines ..
      EXTERNAL           DGEMM, DGETF2, DLASWP, DTRSM, XERBLA
!     ..
!     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETRF', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 )                                          &
     &   RETURN
!
!     Determine the block size for this environment.
!
      NB = ILAENV( 1, 'DGETRF', ' ', M, N, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) THEN
!
!        Use unblocked code.
!
         CALL DGETF2( M, N, A, LDA, IPIV, INFO )
      ELSE
!
!        Use blocked code.
!
         DO 20 J = 1, MIN( M, N ), NB
            JB = MIN( MIN( M, N )-J+1, NB )
!
!           Factor diagonal and subdiagonal blocks and test for exact
!           singularity.
!
            CALL DGETF2( M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO )
!
!           Adjust INFO and the pivot indices.
!
            IF( INFO.EQ.0 .AND. IINFO.GT.0 )                            &
     &         INFO = IINFO + J - 1
            DO 10 I = J, MIN( M, J+JB-1 )
               IPIV( I ) = J - 1 + IPIV( I )
   10       CONTINUE
!
!           Apply interchanges to columns 1:J-1.
!
            CALL DLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )
!
            IF( J+JB.LE.N ) THEN
!
!              Apply interchanges to columns J+JB:N.
!
               CALL DLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1,     &
     &                      IPIV, 1 )
!
!              Compute block row of U.
!
               CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB, &
     &                     N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ), &
     &                     LDA )
               IF( J+JB.LE.M ) THEN
!
!                 Update trailing submatrix.
!
                  CALL DGEMM( 'No transpose', 'No transpose', M-J-JB+1, &
     &                        N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA,    &
     &                        A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ),  &
     &                        LDA )
               END IF
            END IF
   20    CONTINUE
      END IF
      RETURN
!
!     End of DGETRF
!
      END
      SUBROUTINE DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, LDA, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  DGETRS solves a system of linear equations
!     A * X = B  or  A' * X = B
!  with a general N-by-N matrix A using the LU factorization computed
!  by DGETRF.
!
!  Arguments
!  =========
!
!  TRANS   (input) CHARACTER*1
!          Specifies the form of the system of equations:
!          = 'N':  A * X = B  (No transpose)
!          = 'T':  A'* X = B  (Transpose)
!          = 'C':  A'* X = B  (Conjugate transpose = Transpose)
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
!          The factors L and U from the factorization A = P*L*U
!          as computed by DGETRF.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  IPIV    (input) INTEGER array, dimension (N)
!          The pivot indices from DGETRF; for 1<=i<=N, row i of the
!          matrix was interchanged with row IPIV(i).
!
!  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
!          On entry, the right hand side matrix B.
!          On exit, the solution matrix X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            NOTRAN
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLASWP, DTRSM, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.        &
     &    LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETRS', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 .OR. NRHS.EQ.0 )                                       &
     &   RETURN
!
      IF( NOTRAN ) THEN
!
!        Solve A * X = B.
!
!        Apply row interchanges to the right hand sides.
!
         CALL DLASWP( NRHS, B, LDB, 1, N, IPIV, 1 )
!
!        Solve L*X = B, overwriting B with X.
!
         CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', N, NRHS,  &
     &               ONE, A, LDA, B, LDB )
!
!        Solve U*X = B, overwriting B with X.
!
         CALL DTRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N,    &
     &               NRHS, ONE, A, LDA, B, LDB )
      ELSE
!
!        Solve A' * X = B.
!
!        Solve U'*X = B, overwriting B with X.
!
         CALL DTRSM( 'Left', 'Upper', 'Transpose', 'Non-unit', N, NRHS, &
     &               ONE, A, LDA, B, LDB )
!
!        Solve L'*X = B, overwriting B with X.
!
         CALL DTRSM( 'Left', 'Lower', 'Transpose', 'Unit', N, NRHS, ONE,&
     &               A, LDA, B, LDB )
!
!        Apply row interchanges to the solution vectors.
!
         CALL DLASWP( NRHS, B, LDB, 1, N, IPIV, -1 )
      END IF
!
      RETURN
!
!     End of DGETRS
!
      END
      SUBROUTINE DLASWP( N, A, LDA, K1, K2, IPIV, INCX )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      INTEGER            INCX, K1, K2, LDA, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  DLASWP performs a series of row interchanges on the matrix A.
!  One row interchange is initiated for each of rows K1 through K2 of A.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the matrix of column dimension N to which the row
!          interchanges will be applied.
!          On exit, the permuted matrix.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.
!
!  K1      (input) INTEGER
!          The first element of IPIV for which a row interchange will
!          be done.
!
!  K2      (input) INTEGER
!          The last element of IPIV for which a row interchange will
!          be done.
!
!  IPIV    (input) INTEGER array, dimension (M*abs(INCX))
!          The vector of pivot indices.  Only the elements in positions
!          K1 through K2 of IPIV are accessed.
!          IPIV(K) = L implies rows K and L are to be interchanged.
!
!  INCX    (input) INTEGER
!          The increment between successive values of IPIV.  If IPIV
!          is negative, the pivots are applied in reverse order.
!
!  Further Details
!  ===============
!
!  Modified by
!   R. C. Whaley, Computer Science Dept., Univ. of Tenn., Knoxville, USA
!
! =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I, I1, I2, INC, IP, IX, IX0, J, K, N32
      DOUBLE PRECISION   TEMP
!     ..
!     .. Executable Statements ..
!
!     Interchange row I with row IPIV(I) for each of rows K1 through K2.
!
      IF( INCX.GT.0 ) THEN
         IX0 = K1
         I1 = K1
         I2 = K2
         INC = 1
      ELSE IF( INCX.LT.0 ) THEN
         IX0 = 1 + ( 1-K2 )*INCX
         I1 = K2
         I2 = K1
         INC = -1
      ELSE
         RETURN
      END IF
!
      N32 = ( N / 32 )*32
      IF( N32.NE.0 ) THEN
         DO 30 J = 1, N32, 32
            IX = IX0
            DO 20 I = I1, I2, INC
               IP = IPIV( IX )
               IF( IP.NE.I ) THEN
                  DO 10 K = J, J + 31
                     TEMP = A( I, K )
                     A( I, K ) = A( IP, K )
                     A( IP, K ) = TEMP
   10             CONTINUE
               END IF
               IX = IX + INCX
   20       CONTINUE
   30    CONTINUE
      END IF
      IF( N32.NE.N ) THEN
         N32 = N32 + 1
         IX = IX0
         DO 50 I = I1, I2, INC
            IP = IPIV( IX )
            IF( IP.NE.I ) THEN
               DO 40 K = N32, N
                  TEMP = A( I, K )
                  A( I, K ) = A( IP, K )
                  A( IP, K ) = TEMP
   40          CONTINUE
            END IF
            IX = IX + INCX
   50    CONTINUE
      END IF
!
      RETURN
!
!     End of DLASWP
!
      END
      REAL             FUNCTION SCSUM1( N, CX, INCX )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      INTEGER            INCX, N
!     ..
!     .. Array Arguments ..
      COMPLEX            CX( * )
!     ..
!
!  Purpose
!  =======
!
!  SCSUM1 takes the sum of the absolute values of a complex
!  vector and returns a single precision result.
!
!  Based on SCASUM from the Level 1 BLAS.
!  The change is to use the 'genuine' absolute value.
!
!  Contributed by Nick Higham for use with CLACON.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The number of elements in the vector CX.
!
!  CX      (input) COMPLEX array, dimension (N)
!          The vector whose elements will be summed.
!
!  INCX    (input) INTEGER
!          The spacing between successive values of CX.  INCX > 0.
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I, NINCX
      REAL               STEMP
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS
!     ..
!     .. Executable Statements ..
!
      SCSUM1 = 0.0E0
      STEMP = 0.0E0
      IF( N.LE.0 )                                                      &
     &   RETURN
      IF( INCX.EQ.1 )                                                   &
     &   GO TO 20
!
!     CODE FOR INCREMENT NOT EQUAL TO 1
!
      NINCX = N*INCX
      DO 10 I = 1, NINCX, INCX
!
!        NEXT LINE MODIFIED.
!
         STEMP = STEMP + ABS( CX( I ) )
   10 CONTINUE
      SCSUM1 = STEMP
      RETURN
!
!     CODE FOR INCREMENT EQUAL TO 1
!
   20 CONTINUE
      DO 30 I = 1, N
!
!        NEXT LINE MODIFIED.
!
         STEMP = STEMP + ABS( CX( I ) )
   30 CONTINUE
      SCSUM1 = STEMP
      RETURN
!
!     End of SCSUM1
!
      END
      INTEGER          FUNCTION ICMAX1( N, CX, INCX )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      INTEGER            INCX, N
!     ..
!     .. Array Arguments ..
      COMPLEX            CX( * )
!     ..
!
!  Purpose
!  =======
!
!  ICMAX1 finds the index of the element whose real part has maximum
!  absolute value.
!
!  Based on ICAMAX from Level 1 BLAS.
!  The change is to use the 'genuine' absolute value.
!
!  Contributed by Nick Higham for use with CLACON.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The number of elements in the vector CX.
!
!  CX      (input) COMPLEX array, dimension (N)
!          The vector whose elements will be summed.
!
!  INCX    (input) INTEGER
!          The spacing between successive values of CX.  INCX >= 1.
!
! =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I, IX
      REAL               SMAX
      COMPLEX            ZDUM
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, REAL
!     ..
!     .. Statement Functions ..
      REAL               CABS1
!     ..
!     .. Statement Function definitions ..
!
!     NEXT LINE IS THE ONLY MODIFICATION.
      CABS1( ZDUM ) = ABS( REAL( ZDUM ) )
!     ..
!     .. Executable Statements ..
!
      ICMAX1 = 0
      IF( N.LT.1 )                                                      &
     &   RETURN
      ICMAX1 = 1
      IF( N.EQ.1 )                                                      &
     &   RETURN
      IF( INCX.EQ.1 )                                                   &
     &   GO TO 30
!
!     CODE FOR INCREMENT NOT EQUAL TO 1
!
      IX = 1
      SMAX = CABS1( CX( 1 ) )
      IX = IX + INCX
      DO 20 I = 2, N
         IF( CABS1( CX( IX ) ).LE.SMAX )                                &
     &      GO TO 10
         ICMAX1 = I
         SMAX = CABS1( CX( IX ) )
   10    CONTINUE
         IX = IX + INCX
   20 CONTINUE
      RETURN
!
!     CODE FOR INCREMENT EQUAL TO 1
!
   30 CONTINUE
      SMAX = CABS1( CX( 1 ) )
      DO 40 I = 2, N
         IF( CABS1( CX( I ) ).LE.SMAX )                                 &
     &      GO TO 40
         ICMAX1 = I
         SMAX = CABS1( CX( I ) )
   40 CONTINUE
      RETURN
!
!     End of ICMAX1
!
      END
      INTEGER          FUNCTION IEEECK( ISPEC, ZERO, ONE )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1998
!
!     .. Scalar Arguments ..
      INTEGER            ISPEC
      REAL               ONE, ZERO
!     ..
!
!  Purpose
!  =======
!
!  IEEECK is called from the ILAENV to verify that Infinity and
!  possibly NaN arithmetic is safe (i.e. will not trap).
!
!  Arguments
!  =========
!
!  ISPEC   (input) INTEGER
!          Specifies whether to test just for inifinity arithmetic
!          or whether to test for infinity and NaN arithmetic.
!          = 0: Verify infinity arithmetic only.
!          = 1: Verify infinity and NaN arithmetic.
!
!  ZERO    (input) REAL
!          Must contain the value 0.0
!          This is passed to prevent the compiler from optimizing
!          away this code.
!
!  ONE     (input) REAL
!          Must contain the value 1.0
!          This is passed to prevent the compiler from optimizing
!          away this code.
!
!  RETURN VALUE:  INTEGER
!          = 0:  Arithmetic failed to produce the correct answers
!          = 1:  Arithmetic produced the correct answers
!
!     .. Local Scalars ..
      REAL               NAN1, NAN2, NAN3, NAN4, NAN5, NAN6, NEGINF,    &
     &                   NEGZRO, NEWZRO, POSINF
!     ..
!     .. Executable Statements ..
      IEEECK = 1
!
      POSINF = ONE / ZERO
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      NEGINF = -ONE / ZERO
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      NEGZRO = ONE / ( NEGINF+ONE )
      IF( NEGZRO.NE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      NEGINF = ONE / NEGZRO
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      NEWZRO = NEGZRO + ZERO
      IF( NEWZRO.NE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      POSINF = ONE / NEWZRO
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      NEGINF = NEGINF*POSINF
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      POSINF = POSINF*POSINF
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
!
!
!
!
!     Return if we were only asked to check infinity arithmetic
!
      IF( ISPEC.EQ.0 )                                                  &
     &   RETURN
!
      NAN1 = POSINF + NEGINF
!
      NAN2 = POSINF / NEGINF
!
      NAN3 = POSINF / POSINF
!
      NAN4 = POSINF*ZERO
!
      NAN5 = NEGINF*NEGZRO
!
      NAN6 = NAN5*0.0
!
      IF( NAN1.EQ.NAN1 ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      IF( NAN2.EQ.NAN2 ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      IF( NAN3.EQ.NAN3 ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      IF( NAN4.EQ.NAN4 ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      IF( NAN5.EQ.NAN5 ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      IF( NAN6.EQ.NAN6 ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      RETURN
      END
      REAL             FUNCTION SLAMCH( CMACH )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992 
!
!     .. Scalar Arguments ..
      CHARACTER          CMACH
!     ..
!
!  Purpose
!  =======
!
!  SLAMCH determines single precision machine parameters.
!
!  Arguments
!  =========
!
!  CMACH   (input) CHARACTER*1
!          Specifies the value to be returned by SLAMCH:
!          = 'E' or 'e',   SLAMCH := eps
!          = 'S' or 's ,   SLAMCH := sfmin
!          = 'B' or 'b',   SLAMCH := base
!          = 'P' or 'p',   SLAMCH := eps*base
!          = 'N' or 'n',   SLAMCH := t
!          = 'R' or 'r',   SLAMCH := rnd
!          = 'M' or 'm',   SLAMCH := emin
!          = 'U' or 'u',   SLAMCH := rmin
!          = 'L' or 'l',   SLAMCH := emax
!          = 'O' or 'o',   SLAMCH := rmax
!
!          where
!
!          eps   = relative machine precision
!          sfmin = safe minimum, such that 1/sfmin does not overflow
!          base  = base of the machine
!          prec  = eps*base
!          t     = number of (base) digits in the mantissa
!          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
!          emin  = minimum exponent before (gradual) underflow
!          rmin  = underflow threshold - base**(emin-1)
!          emax  = largest exponent before overflow
!          rmax  = overflow threshold  - (base**emax)*(1-eps)
!
! =====================================================================
!
!     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            FIRST, LRND
      INTEGER            BETA, IMAX, IMIN, IT
      REAL               BASE, EMAX, EMIN, EPS, PREC, RMACH, RMAX, RMIN,&
     &                   RND, SFMIN, SMALL, T
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           SLAMC2
!     ..
!     .. Save statement ..
      SAVE               FIRST, EPS, SFMIN, BASE, T, RND, EMIN, RMIN,   &
     &                   EMAX, RMAX, PREC
!     ..
!     .. Data statements ..
      DATA               FIRST / .TRUE. /
!     ..
!     .. Executable Statements ..
!
      IF( FIRST ) THEN
         FIRST = .FALSE.
         CALL SLAMC2( BETA, IT, LRND, EPS, IMIN, RMIN, IMAX, RMAX )
         BASE = BETA
         T = IT
         IF( LRND ) THEN
            RND = ONE
            EPS = ( BASE**( 1-IT ) ) / 2
         ELSE
            RND = ZERO
            EPS = BASE**( 1-IT )
         END IF
         PREC = EPS*BASE
         EMIN = IMIN
         EMAX = IMAX
         SFMIN = RMIN
         SMALL = ONE / RMAX
         IF( SMALL.GE.SFMIN ) THEN
!
!           Use SMALL plus a bit, to avoid the possibility of rounding
!           causing overflow when computing  1/sfmin.
!
            SFMIN = SMALL*( ONE+EPS )
         END IF
      END IF
!
      IF( LSAME( CMACH, 'E' ) ) THEN
         RMACH = EPS
      ELSE IF( LSAME( CMACH, 'S' ) ) THEN
         RMACH = SFMIN
      ELSE IF( LSAME( CMACH, 'B' ) ) THEN
         RMACH = BASE
      ELSE IF( LSAME( CMACH, 'P' ) ) THEN
         RMACH = PREC
      ELSE IF( LSAME( CMACH, 'N' ) ) THEN
         RMACH = T
      ELSE IF( LSAME( CMACH, 'R' ) ) THEN
         RMACH = RND
      ELSE IF( LSAME( CMACH, 'M' ) ) THEN
         RMACH = EMIN
      ELSE IF( LSAME( CMACH, 'U' ) ) THEN
         RMACH = RMIN
      ELSE IF( LSAME( CMACH, 'L' ) ) THEN
         RMACH = EMAX
      ELSE IF( LSAME( CMACH, 'O' ) ) THEN
         RMACH = RMAX
      END IF
!
      SLAMCH = RMACH
      RETURN
!
!     End of SLAMCH
!
      END
!
!***********************************************************************
!
      SUBROUTINE SLAMC1( BETA, T, RND, IEEE1 )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      LOGICAL            IEEE1, RND
      INTEGER            BETA, T
!     ..
!
!  Purpose
!  =======
!
!  SLAMC1 determines the machine parameters given by BETA, T, RND, and
!  IEEE1.
!
!  Arguments
!  =========
!
!  BETA    (output) INTEGER
!          The base of the machine.
!
!  T       (output) INTEGER
!          The number of ( BETA ) digits in the mantissa.
!
!  RND     (output) LOGICAL
!          Specifies whether proper rounding  ( RND = .TRUE. )  or
!          chopping  ( RND = .FALSE. )  occurs in addition. This may not
!          be a reliable guide to the way in which the machine performs
!          its arithmetic.
!
!  IEEE1   (output) LOGICAL
!          Specifies whether rounding appears to be done in the IEEE
!          'round to nearest' style.
!
!  Further Details
!  ===============
!
!  The routine is based on the routine  ENVRON  by Malcolm and
!  incorporates suggestions by Gentleman and Marovich. See
!
!     Malcolm M. A. (1972) Algorithms to reveal properties of
!        floating-point arithmetic. Comms. of the ACM, 15, 949-951.
!
!     Gentleman W. M. and Marovich S. B. (1974) More on algorithms
!        that reveal properties of floating point arithmetic units.
!        Comms. of the ACM, 17, 276-277.
!
! =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            FIRST, LIEEE1, LRND
      INTEGER            LBETA, LT
      REAL               A, B, C, F, ONE, QTR, SAVEC, T1, T2
!     ..
!     .. External Functions ..
      REAL               SLAMC3
      EXTERNAL           SLAMC3
!     ..
!     .. Save statement ..
      SAVE               FIRST, LIEEE1, LBETA, LRND, LT
!     ..
!     .. Data statements ..
      DATA               FIRST / .TRUE. /
!     ..
!     .. Executable Statements ..
!
      IF( FIRST ) THEN
         FIRST = .FALSE.
         ONE = 1
!
!        LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BETA,
!        IEEE1, T and RND.
!
!        Throughout this routine  we use the function  SLAMC3  to ensure
!        that relevant values are  stored and not held in registers,  or
!        are not affected by optimizers.
!
!        Compute  a = 2.0**m  with the  smallest positive integer m such
!        that
!
!           fl( a + 1.0 ) = a.
!
         A = 1
         C = 1
!
!+       WHILE( C.EQ.ONE )LOOP
   10    CONTINUE
         IF( C.EQ.ONE ) THEN
            A = 2*A
            C = SLAMC3( A, ONE )
            C = SLAMC3( C, -A )
            GO TO 10
         END IF
!+       END WHILE
!
!        Now compute  b = 2.0**m  with the smallest positive integer m
!        such that
!
!           fl( a + b ) .gt. a.
!
         B = 1
         C = SLAMC3( A, B )
!
!+       WHILE( C.EQ.A )LOOP
   20    CONTINUE
         IF( C.EQ.A ) THEN
            B = 2*B
            C = SLAMC3( A, B )
            GO TO 20
         END IF
!+       END WHILE
!
!        Now compute the base.  a and c  are neighbouring floating point
!        numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and so
!        their difference is beta. Adding 0.25 to c is to ensure that it
!        is truncated to beta and not ( beta - 1 ).
!
         QTR = ONE / 4
         SAVEC = C
         C = SLAMC3( C, -A )
         LBETA = C + QTR
!
!        Now determine whether rounding or chopping occurs,  by adding a
!        bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a.
!
         B = LBETA
         F = SLAMC3( B / 2, -B / 100 )
         C = SLAMC3( F, A )
         IF( C.EQ.A ) THEN
            LRND = .TRUE.
         ELSE
            LRND = .FALSE.
         END IF
         F = SLAMC3( B / 2, B / 100 )
         C = SLAMC3( F, A )
         IF( ( LRND ) .AND. ( C.EQ.A ) )                                &
     &      LRND = .FALSE.
!
!        Try and decide whether rounding is done in the  IEEE  'round to
!        nearest' style. B/2 is half a unit in the last place of the two
!        numbers A and SAVEC. Furthermore, A is even, i.e. has last  bit
!        zero, and SAVEC is odd. Thus adding B/2 to A should not  change
!        A, but adding B/2 to SAVEC should change SAVEC.
!
         T1 = SLAMC3( B / 2, A )
         T2 = SLAMC3( B / 2, SAVEC )
         LIEEE1 = ( T1.EQ.A ) .AND. ( T2.GT.SAVEC ) .AND. LRND
!
!        Now find  the  mantissa, t.  It should  be the  integer part of
!        log to the base beta of a,  however it is safer to determine  t
!        by powering.  So we find t as the smallest positive integer for
!        which
!
!           fl( beta**t + 1.0 ) = 1.0.
!
         LT = 0
         A = 1
         C = 1
!
!+       WHILE( C.EQ.ONE )LOOP
   30    CONTINUE
         IF( C.EQ.ONE ) THEN
            LT = LT + 1
            A = A*LBETA
            C = SLAMC3( A, ONE )
            C = SLAMC3( C, -A )
            GO TO 30
         END IF
!+       END WHILE
!
      END IF
!
      BETA = LBETA
      T = LT
      RND = LRND
      IEEE1 = LIEEE1
      RETURN
!
!     End of SLAMC1
!
      END
!
!***********************************************************************
!
      SUBROUTINE SLAMC2( BETA, T, RND, EPS, EMIN, RMIN, EMAX, RMAX )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      LOGICAL            RND
      INTEGER            BETA, EMAX, EMIN, T
      REAL               EPS, RMAX, RMIN
!     ..
!
!  Purpose
!  =======
!
!  SLAMC2 determines the machine parameters specified in its argument
!  list.
!
!  Arguments
!  =========
!
!  BETA    (output) INTEGER
!          The base of the machine.
!
!  T       (output) INTEGER
!          The number of ( BETA ) digits in the mantissa.
!
!  RND     (output) LOGICAL
!          Specifies whether proper rounding  ( RND = .TRUE. )  or
!          chopping  ( RND = .FALSE. )  occurs in addition. This may not
!          be a reliable guide to the way in which the machine performs
!          its arithmetic.
!
!  EPS     (output) REAL
!          The smallest positive number such that
!
!             fl( 1.0 - EPS ) .LT. 1.0,
!
!          where fl denotes the computed value.
!
!  EMIN    (output) INTEGER
!          The minimum exponent before (gradual) underflow occurs.
!
!  RMIN    (output) REAL
!          The smallest normalized number for the machine, given by
!          BASE**( EMIN - 1 ), where  BASE  is the floating point value
!          of BETA.
!
!  EMAX    (output) INTEGER
!          The maximum exponent before overflow occurs.
!
!  RMAX    (output) REAL
!          The largest positive number for the machine, given by
!          BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating point
!          value of BETA.
!
!  Further Details
!  ===============
!
!  The computation of  EPS  is based on a routine PARANOIA by
!  W. Kahan of the University of California at Berkeley.
!
! =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            FIRST, IEEE, IWARN, LIEEE1, LRND
      INTEGER            GNMIN, GPMIN, I, LBETA, LEMAX, LEMIN, LT,      &
     &                   NGNMIN, NGPMIN
      REAL               A, B, C, HALF, LEPS, LRMAX, LRMIN, ONE, RBASE, &
     &                   SIXTH, SMALL, THIRD, TWO, ZERO
!     ..
!     .. External Functions ..
      REAL               SLAMC3
      EXTERNAL           SLAMC3
!     ..
!     .. External Subroutines ..
      EXTERNAL           SLAMC1, SLAMC4, SLAMC5
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
!     ..
!     .. Save statement ..
      SAVE               FIRST, IWARN, LBETA, LEMAX, LEMIN, LEPS, LRMAX,&
     &                   LRMIN, LT
!     ..
!     .. Data statements ..
      DATA               FIRST / .TRUE. / , IWARN / .FALSE. /
!     ..
!     .. Executable Statements ..
!
      IF( FIRST ) THEN
         FIRST = .FALSE.
         ZERO = 0
         ONE = 1
         TWO = 2
!
!        LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values of
!        BETA, T, RND, EPS, EMIN and RMIN.
!
!        Throughout this routine  we use the function  SLAMC3  to ensure
!        that relevant values are stored  and not held in registers,  or
!        are not affected by optimizers.
!
!        SLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1.
!
         CALL SLAMC1( LBETA, LT, LRND, LIEEE1 )
!
!        Start to find EPS.
!
         B = LBETA
         A = B**( -LT )
         LEPS = A
!
!        Try some tricks to see whether or not this is the correct  EPS.
!
         B = TWO / 3
         HALF = ONE / 2
         SIXTH = SLAMC3( B, -HALF )
         THIRD = SLAMC3( SIXTH, SIXTH )
         B = SLAMC3( THIRD, -HALF )
         B = SLAMC3( B, SIXTH )
         B = ABS( B )
         IF( B.LT.LEPS )                                                &
     &      B = LEPS
!
         LEPS = 1
!
!+       WHILE( ( LEPS.GT.B ).AND.( B.GT.ZERO ) )LOOP
   10    CONTINUE
         IF( ( LEPS.GT.B ) .AND. ( B.GT.ZERO ) ) THEN
            LEPS = B
            C = SLAMC3( HALF*LEPS, ( TWO**5 )*( LEPS**2 ) )
            C = SLAMC3( HALF, -C )
            B = SLAMC3( HALF, C )
            C = SLAMC3( HALF, -B )
            B = SLAMC3( HALF, C )
            GO TO 10
         END IF
!+       END WHILE
!
         IF( A.LT.LEPS )                                                &
     &      LEPS = A
!
!        Computation of EPS complete.
!
!        Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3)).
!        Keep dividing  A by BETA until (gradual) underflow occurs. This
!        is detected when we cannot recover the previous A.
!
         RBASE = ONE / LBETA
         SMALL = ONE
         DO 20 I = 1, 3
            SMALL = SLAMC3( SMALL*RBASE, ZERO )
   20    CONTINUE
         A = SLAMC3( ONE, SMALL )
         CALL SLAMC4( NGPMIN, ONE, LBETA )
         CALL SLAMC4( NGNMIN, -ONE, LBETA )
         CALL SLAMC4( GPMIN, A, LBETA )
         CALL SLAMC4( GNMIN, -A, LBETA )
         IEEE = .FALSE.
!
         IF( ( NGPMIN.EQ.NGNMIN ) .AND. ( GPMIN.EQ.GNMIN ) ) THEN
            IF( NGPMIN.EQ.GPMIN ) THEN
               LEMIN = NGPMIN
!            ( Non twos-complement machines, no gradual underflow;
!              e.g.,  VAX )
            ELSE IF( ( GPMIN-NGPMIN ).EQ.3 ) THEN
               LEMIN = NGPMIN - 1 + LT
               IEEE = .TRUE.
!            ( Non twos-complement machines, with gradual underflow;
!              e.g., IEEE standard followers )
            ELSE
               LEMIN = MIN( NGPMIN, GPMIN )
!            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
!
         ELSE IF( ( NGPMIN.EQ.GPMIN ) .AND. ( NGNMIN.EQ.GNMIN ) ) THEN
            IF( ABS( NGPMIN-NGNMIN ).EQ.1 ) THEN
               LEMIN = MAX( NGPMIN, NGNMIN )
!            ( Twos-complement machines, no gradual underflow;
!              e.g., CYBER 205 )
            ELSE
               LEMIN = MIN( NGPMIN, NGNMIN )
!            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
!
         ELSE IF( ( ABS( NGPMIN-NGNMIN ).EQ.1 ) .AND.                   &
     &            ( GPMIN.EQ.GNMIN ) ) THEN
            IF( ( GPMIN-MIN( NGPMIN, NGNMIN ) ).EQ.3 ) THEN
               LEMIN = MAX( NGPMIN, NGNMIN ) - 1 + LT
!            ( Twos-complement machines with gradual underflow;
!              no known machine )
            ELSE
               LEMIN = MIN( NGPMIN, NGNMIN )
!            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
!
         ELSE
            LEMIN = MIN( NGPMIN, NGNMIN, GPMIN, GNMIN )
!         ( A guess; no known machine )
            IWARN = .TRUE.
         END IF
!**
! Comment out this if block if EMIN is ok
         IF( IWARN ) THEN
            FIRST = .TRUE.
            WRITE( 6, FMT = 9999 )LEMIN
         END IF
!**
!
!        Assume IEEE arithmetic if we found denormalised  numbers above,
!        or if arithmetic seems to round in the  IEEE style,  determined
!        in routine SLAMC1. A true IEEE machine should have both  things
!        true; however, faulty machines may have one or the other.
!
         IEEE = IEEE .OR. LIEEE1
!
!        Compute  RMIN by successive division by  BETA. We could compute
!        RMIN as BASE**( EMIN - 1 ),  but some machines underflow during
!        this computation.
!
         LRMIN = 1
         DO 30 I = 1, 1 - LEMIN
            LRMIN = SLAMC3( LRMIN*RBASE, ZERO )
   30    CONTINUE
!
!        Finally, call SLAMC5 to compute EMAX and RMAX.
!
         CALL SLAMC5( LBETA, LT, LEMIN, IEEE, LEMAX, LRMAX )
      END IF
!
      BETA = LBETA
      T = LT
      RND = LRND
      EPS = LEPS
      EMIN = LEMIN
      RMIN = LRMIN
      EMAX = LEMAX
      RMAX = LRMAX
!
      RETURN
!
 9999 FORMAT( / / ' WARNING. The value EMIN may be incorrect:-',        &
     &      '  EMIN = ', I8, /                                          &
     &      ' If, after inspection, the value EMIN looks',              &
     &      ' acceptable please comment out ',                          &
     &      / ' the IF block as marked within the code of routine',     &
     &      ' SLAMC2,', / ' otherwise supply EMIN explicitly.', / )
!
!     End of SLAMC2
!
      END
!
!***********************************************************************
!
      REAL             FUNCTION SLAMC3( A, B )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      REAL               A, B
!     ..
!
!  Purpose
!  =======
!
!  SLAMC3  is intended to force  A  and  B  to be stored prior to doing
!  the addition of  A  and  B ,  for use in situations where optimizers
!  might hold one of these in a register.
!
!  Arguments
!  =========
!
!  A, B    (input) REAL
!          The values A and B.
!
! =====================================================================
!
!     .. Executable Statements ..
!
      SLAMC3 = A + B
!
      RETURN
!
!     End of SLAMC3
!
      END
!
!***********************************************************************
!
      SUBROUTINE SLAMC4( EMIN, START, BASE )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      INTEGER            BASE, EMIN
      REAL               START
!     ..
!
!  Purpose
!  =======
!
!  SLAMC4 is a service routine for SLAMC2.
!
!  Arguments
!  =========
!
!  EMIN    (output) EMIN
!          The minimum exponent before (gradual) underflow, computed by
!          setting A = START and dividing by BASE until the previous A
!          can not be recovered.
!
!  START   (input) REAL
!          The starting point for determining EMIN.
!
!  BASE    (input) INTEGER
!          The base of the machine.
!
! =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I
      REAL               A, B1, B2, C1, C2, D1, D2, ONE, RBASE, ZERO
!     ..
!     .. External Functions ..
      REAL               SLAMC3
      EXTERNAL           SLAMC3
!     ..
!     .. Executable Statements ..
!
      A = START
      ONE = 1
      RBASE = ONE / BASE
      ZERO = 0
      EMIN = 1
      B1 = SLAMC3( A*RBASE, ZERO )
      C1 = A
      C2 = A
      D1 = A
      D2 = A
!+    WHILE( ( C1.EQ.A ).AND.( C2.EQ.A ).AND.
!    $       ( D1.EQ.A ).AND.( D2.EQ.A )      )LOOP
   10 CONTINUE
      IF( ( C1.EQ.A ) .AND. ( C2.EQ.A ) .AND. ( D1.EQ.A ) .AND.         &
     &    ( D2.EQ.A ) ) THEN
         EMIN = EMIN - 1
         A = B1
         B1 = SLAMC3( A / BASE, ZERO )
         C1 = SLAMC3( B1*BASE, ZERO )
         D1 = ZERO
         DO 20 I = 1, BASE
            D1 = D1 + B1
   20    CONTINUE
         B2 = SLAMC3( A*RBASE, ZERO )
         C2 = SLAMC3( B2 / RBASE, ZERO )
         D2 = ZERO
         DO 30 I = 1, BASE
            D2 = D2 + B2
   30    CONTINUE
         GO TO 10
      END IF
!+    END WHILE
!
      RETURN
!
!     End of SLAMC4
!
      END
!
!***********************************************************************
!
      SUBROUTINE SLAMC5( BETA, P, EMIN, IEEE, EMAX, RMAX )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      LOGICAL            IEEE
      INTEGER            BETA, EMAX, EMIN, P
      REAL               RMAX
!     ..
!
!  Purpose
!  =======
!
!  SLAMC5 attempts to compute RMAX, the largest machine floating-point
!  number, without overflow.  It assumes that EMAX + abs(EMIN) sum
!  approximately to a power of 2.  It will fail on machines where this
!  assumption does not hold, for example, the Cyber 205 (EMIN = -28625,
!  EMAX = 28718).  It will also fail if the value supplied for EMIN is
!  too large (i.e. too close to zero), probably with overflow.
!
!  Arguments
!  =========
!
!  BETA    (input) INTEGER
!          The base of floating-point arithmetic.
!
!  P       (input) INTEGER
!          The number of base BETA digits in the mantissa of a
!          floating-point value.
!
!  EMIN    (input) INTEGER
!          The minimum exponent before (gradual) underflow.
!
!  IEEE    (input) LOGICAL
!          A logical flag specifying whether or not the arithmetic
!          system is thought to comply with the IEEE standard.
!
!  EMAX    (output) INTEGER
!          The largest exponent before overflow
!
!  RMAX    (output) REAL
!          The largest machine floating-point number.
!
! =====================================================================
!
!     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
!     ..
!     .. Local Scalars ..
      INTEGER            EXBITS, EXPSUM, I, LEXP, NBITS, TRY, UEXP
      REAL               OLDY, RECBAS, Y, Z
!     ..
!     .. External Functions ..
      REAL               SLAMC3
      EXTERNAL           SLAMC3
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MOD
!     ..
!     .. Executable Statements ..
!
!     First compute LEXP and UEXP, two powers of 2 that bound
!     abs(EMIN). We then assume that EMAX + abs(EMIN) will sum
!     approximately to the bound that is closest to abs(EMIN).
!     (EMAX is the exponent of the required number RMAX).
!
      LEXP = 1
      EXBITS = 1
   10 CONTINUE
      TRY = LEXP*2
      IF( TRY.LE.( -EMIN ) ) THEN
         LEXP = TRY
         EXBITS = EXBITS + 1
         GO TO 10
      END IF
      IF( LEXP.EQ.-EMIN ) THEN
         UEXP = LEXP
      ELSE
         UEXP = TRY
         EXBITS = EXBITS + 1
      END IF
!
!     Now -LEXP is less than or equal to EMIN, and -UEXP is greater
!     than or equal to EMIN. EXBITS is the number of bits needed to
!     store the exponent.
!
      IF( ( UEXP+EMIN ).GT.( -LEXP-EMIN ) ) THEN
         EXPSUM = 2*LEXP
      ELSE
         EXPSUM = 2*UEXP
      END IF
!
!     EXPSUM is the exponent range, approximately equal to
!     EMAX - EMIN + 1 .
!
      EMAX = EXPSUM + EMIN - 1
      NBITS = 1 + EXBITS + P
!
!     NBITS is the total number of bits needed to store a
!     floating-point number.
!
      IF( ( MOD( NBITS, 2 ).EQ.1 ) .AND. ( BETA.EQ.2 ) ) THEN
!
!        Either there are an odd number of bits used to store a
!        floating-point number, which is unlikely, or some bits are
!        not used in the representation of numbers, which is possible,
!        (e.g. Cray machines) or the mantissa has an implicit bit,
!        (e.g. IEEE machines, Dec Vax machines), which is perhaps the
!        most likely. We have to assume the last alternative.
!        If this is true, then we need to reduce EMAX by one because
!        there must be some way of representing zero in an implicit-bit
!        system. On machines like Cray, we are reducing EMAX by one
!        unnecessarily.
!
         EMAX = EMAX - 1
      END IF
!
      IF( IEEE ) THEN
!
!        Assume we are on an IEEE machine which reserves one exponent
!        for infinity and NaN.
!
         EMAX = EMAX - 1
      END IF
!
!     Now create RMAX, the largest machine number, which should
!     be equal to (1.0 - BETA**(-P)) * BETA**EMAX .
!
!     First compute 1.0 - BETA**(-P), being careful that the
!     result is less than 1.0 .
!
      RECBAS = ONE / BETA
      Z = BETA - ONE
      Y = ZERO
      DO 20 I = 1, P
         Z = Z*RECBAS
         IF( Y.LT.ONE )                                                 &
     &      OLDY = Y
         Y = SLAMC3( Y, Z )
   20 CONTINUE
      IF( Y.GE.ONE )                                                    &
     &   Y = OLDY
!
!     Now multiply by BETA**EMAX to get RMAX.
!
      DO 30 I = 1, EMAX
         Y = SLAMC3( Y*BETA, ZERO )
   30 CONTINUE
!
      RMAX = Y
      RETURN
!
!     End of SLAMC5
!
      END
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
