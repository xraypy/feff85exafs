c     sub-program exchange point
      program rdinp 
c     subroutine rdinp (nabs,ceels)

c    reads 'feff.inp' file and writes several files in special format
c    ready for the use by other modules: geom.dat, global.dat,
c    mod1.inp, mod2.inp, mod3.inp mod4.inp mod5.inp mod6.inp ldos.inp .
c    The subroutine output 'nabs' is needed for configurational average
c    The rest of output, passed to wrtall via common blocks (allinp.h)

c     coded s. zabinski 1994
c     last modified by a.l.ankudinov march 2001  for new i/o structure

      implicit double precision (a-h, o-z)
      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'
      include '../RDINP/allinp.h'
      include '../HEADERS/vers.h'
      include '../HEADERS/parallel.h'

c     Single scattering path to go with Overlap information
      parameter (nssx = 16)
      dimension indss(nssx), iphss(nssx)
      dimension degss(nssx), rss(nssx)

c     Local stuff
      character*150  line
      dimension ltit(nheadx)
      parameter (nwordx = 20)
      character*20 words(nwordx)
      integer iatph(0:nphx)
      integer icnt  !KJ 1-06 this is just a local index that does not need to be saved

      parameter (nbr=30)
      logical nogeom
      parameter (big = 1.0e5)
      character*512 slog
      character*38 tmpstr
      logical ceels  !KJ for monolithic version 5-6

      external dist

      icnt = 0

c   10 format (a)
   20 format (bn, i15)
   30 format (bn, f15.0)

      call par_begin
      if (worker) go to 400

c     open the log file, unit 11.  See subroutine wlog.
      open (unit=11, file='log.dat', status='unknown', iostat=ios)
      call chopen (ios, 'log.dat', 'feff')

      tmpstr = vfeff // 'release '//vf85e
      call triml (tmpstr)
      call wlog(' ' // tmpstr)

c     initialize all things to be passed
      nabs = 1
      call iniall

c     initialize local staff 
      iatom = 0
      ifolp = 0
      iovrlp = 0
      iphabs = 0
      lxnat = 0
      folpx = 1.15d0
      nogeom = .false.
      rclabs = big
      rmult = 1.0d0
      s02h = 1.0d0
      natt = 0
      nss = 0
      do 90  iss = 1, nssx
         indss(iss) = 0
         iphss(iss) = 0
         degss(iss) = 0
         rss(iss) = 0
  90  continue
      do 95 iph = 0, nphx
  95  iatph(iph) = 0

c     tokens  0 if not a token
c             1 if ATOM (ATOMS)
c             2 if HOLE
c             3 if OVER (OVERLAP)
c             4 if CONT (CONTROL)
c             5 if EXCH (EXCHANGE)
c             6 if ION
c             7 if TITL (TITLE)
c             8 if FOLP
c             9 if RPATH or RMAX
c            10 if DEBY (DEBYE)
c            11 if RMUL (RMULTIPLIER)
c            12 if SS
c            13 if PRIN (PRINT)
c            14 if POTE (POTENTIALS)
c            15 if NLEG
c            16 if CRIT (CRITERIA)
c            17 if NOGEOM
c            18 if IORDER
c            19 if PCRI (PCRITERIA)
c            20 if SIG2
c            21 if XANE (XANES)
c            22 if CORR (CORRECTIONS)
c            23 if AFOL (AFOLP)
c            24 if EXAF (EXAFS)
c            25 if POLA (POLARIZATION)
c            26 if ELLI (ELLIPTICITY) 
c            27 if RGRI (RGRID)
c            28 if RPHA (RPHASES), real phase shifts
c            29 if NSTA (NSTAR), n* for co-linear polarization
c            30 if NOHO (NOHOLE), use no hole for potentials
c            31 if SIG3 third and first cumulants for ss paths
c            32 if JUMP (JUMPRM), remove jumps of potential   
c            33 if MBCO (MBCONV), do convolution with exitation spectrum
c            34 if SPIN do calculation for spin-up(down) photoelectron  
c            35 if EDGE to specify edge by name
c            36 if SCF  do self-consistency loop
c            37 if FMS  use FMS for cluster of the size rfms
c            38 if LDOS print out l-dos for specified energy range
c            39 if INTE how to find interstitial parameters
c            40 if CFAV to do configuration average
c            41 if S02  to specify S_0^2
c            45 if RSIG (RSIGMA), real self-energy 
c            46 if XNCD natural dichroism
c            47 if MULT for quadrupolar etc. transitions
c            48 if UNFR unfreeze f-electrons
c            49 if TDLDA use TDLDA background
c            50 if PMBSE use BSE for background
c            51 if PLASMON       - Added by Josh Kas
c                                - PLASMON
c                                - With this card set, ffmod2 will read exc.dat and
c                                - use a multiple pole self energy
c            52 if S02C (S02CONV) compute S_0^2 from response function
c            53 if SELF print on shell self energy as a function of E.
c            54 if SFSE print off shell self energy and spectral function.
c            55 if RCONV print running convolution with spectral function.
c            56 if ELNE calculate ELNES  !KJ 1-06
c            57 if EXEL calculate EXELFS !KJ 1-06
c            58 if MAGI plot magic angle !KJ 1-06
c            59 if ABSO don't normalize spectrum !KJ 3-06
c            60 if EGRID (Gives user control of grid through grid.inp)
c            -1 if END  (end)
c     mode flag  0 ready to read a keyword card
c                1 reading atom positions
c                2 reading overlap instructions for unique pot
c                3 reading unique potential definitions
c                4 reading EELS input  !KJ

c#mn{
c  replaced read of feff.inp with  call to rdline, which will:
c    1. read from feff.inp if found, otherwise will stop and complain
c       (support for reading from standard input would be easy to add)
c    2. handles line processing tasks like 
c         = ignoring comment lines and blank lines
c         = tab removal
c    3. allows 'include' files in input file
c    4. for initial call, set jinit = -1, line = input_file_name
c
      mode  = 0
      jinit = -1
      line  = 'feff.inp'
  200 continue 
         call rdline(jinit,line)
         if (line .eq. 'read_line_end')    line='END'
         if (line .eq. 'read_line_error')  line='END'
c#mn}

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
               iatom  = iatom  +1
            elseif (itok .eq. 2)  then
c              HOLE     1  1.0
c                   holecode s02
               read(words(2),20,err=900)  ihole
               if (nwords.gt.2) read(words(3),30,err=900)  s02h
               mode = 0
            elseif (itok .eq. 3)  then
c              OVERLAP iph
c                  iph  n  r
               read(words(2),20,err=900)  iph
               call phstop(iph,line)
               call warnex(' OVERLAP:')
               mode = 2
               iovrlp = iovrlp +1
            elseif (itok .eq. 4)  then
c              CONTROL  mphase, mpath, mfeff, mchi
c               0 - do not run modules, 1 - run module
               if (nwords.eq.5) then
c                 feff7 input file
                  read(words(2),20,err=900)  mpot
                  mphase = mpot
                  mfms = mpot
                  read(words(3),20,err=900)  mpath
                  read(words(4),20,err=900)  mfeff
                  read(words(5),20,err=900)  mchi
               else
c                 feff8 input file
                  read(words(2),20,err=900)  mpot
                  read(words(3),20,err=900)  mphase
                  read(words(4),20,err=900)  mfms
                  read(words(5),20,err=900)  mpath
                  read(words(6),20,err=900)  mfeff
                  read(words(7),20,err=900)  mchi
               endif
               mode = 0
            elseif (itok .eq. 5)  then
c              EXCHANGE  ixc  vr0  vi0 (ixc0)
c              ixc=0  Hedin-Lunqvist + const real & imag part
c              ixc=1  Dirac-Hara + const real & imag part
c              ixc=2  ground state + const real & imag part
c              ixc=3  Dirac-Hara + HL imag part + const real & imag part
c              ixc=5  partially nonlocal: Dirac-Fock for core + HL for
c                     valence electrons, + const real & imag part
c              ixc=10 same as ixc=0 with broadened plasmon HL selfenergy
c              ixc=13 same as ixc=3 with broadened plasmon HL selfenergy
c              ixc=15 same as ixc=5 with broadened plasmon HL selfenergy
c              vr0 is const imag part of potential
c              vi0 is const imag part of potential
c              Default is HL. (ixc=0, vr0=0, vi0=0, ixc0 = 2)
               vr0=0.0
               vi0=0.0
               read(words(2),20,err=900)  ixc
!              if (nwords.ge.3) (read(words(3),30,err=900)  vr0
                read(words(3),30,err=900)  vr0
!              if (nwords.ge.4) read(words(4),30,err=900)  vi0
                read(words(4),30,err=900)  vi0
               if (nwords .gt. 4) read(words(5),20,err=900)  ixc0
               if (ixc .ge. 3)  call warnex(' EXCHANGE >= 3:')
               mode = 0
            elseif (itok .eq. 6)  then
c              ION  iph xion(iph)
               read(words(2),20,err=900)  iph
               call phstop(iph,line)
               read(words(3),30,err=900)  xion(iph)
               call warnex(' ION:')
               mode = 0
            elseif (itok .eq. 7)  then
c              TITLE title...
               ntitle = ntitle + 1
               if (ntitle .le. nheadx)  then
                  title(ntitle) = line(6:)
                  call triml (title(ntitle))
               else
                  call wlog(' Too many title lines, title ignored')
                  call wlog(' ' // line(1:71))
               endif
               mode = 0
            elseif (itok .eq. 8)  then
c              FOLP iph folp (overlap factor, default 1)
               ifolp = 1
               read(words(2),20,err=900)  iph
               call phstop(iph,line)
               read(words(3),30,err=900)  folp(iph)
               call warnex(' FOLP:')
               mode = 0
            elseif (itok .eq. 9)  then
c              RPATH rmax (max r for ss and pathfinder)
               read(words(2),30,err=900)  rmax
            elseif (itok .eq. 10)  then
c              DEBYE  temp debye-temp ( idwopt )
c                   temps in kelvin
c                   idwopt = 0 use CD model
c                   idwopt = 1 use EM method
c                   idwopt = 2 use RM method
c                   idwopt = -1,-2,... don't calculate DW factors
c                   These add to any sig2 from SIG2 card or files.dat
               read(words(2),30,err=900)  tk
               read(words(3),30,err=900)  thetad
               idwopt=0 
               if (nwords.gt.3) then
                 read(words(4),20,err=900)  idwopt
                 if (idwopt.gt.2) then
                    write(slog,'(a,i5,2x,a)')
     1                 ' Option idwopt=',idwopt,'is not available.'
                    call wlog(slog)
                    write(slog,'(a)')
     1                   '...setting idwopt=2 to use RM.' 
                    call wlog(slog)
                 endif
               endif
               mode = 0
            elseif (itok .eq. 11)  then
c              RMULTIPLIER  rmult
c              Multiples atom coord, rss, overlap and rmax distances by
c              rmult (default 1).  DOES NOT modify sig2g
               read(words(2),30,err=900)  rmult
               mode = 0
            elseif (itok .eq. 12)  then
c              SS index ipot deg rss
               nss = nss + 1
               if (nss .gt. nssx)  then
                  write(slog,'(a,i8)')
     1               ' Too many ss paths requested, max is ', nssx
                  call wlog(slog)
                  call par_stop('RDINP')
               endif
               read(words(2),20,err=900)  indss(nss)
               read(words(3),20,err=900)  iphss(nss)
               read(words(4),30,err=900)  degss(nss)
               read(words(5),30,err=900)  rss(nss)
               mode = 0
            elseif (itok .eq. 13)  then
c              PRINT  ipr1  ipr2  ipr3  ipr4 ipr5 ipr6
c              print flags for various modules
c              ipr1 potph  0 pot.pad only
c                          1 add misc.dat
c                          2 add pot.dat
c                          5 add atom.dat
c                          6 add central atom dirac stuff
c                          7 stop after doing central atom dirac stuff
c              ipr2 xsph   0 phase.pad only
c                          2 add  phase.dat
c                          3 add  emesh.dat
c              ipr3 fmstot  currently is dummy
c              ipr4 pathfinder  0 paths.dat only
c                               1 add crit.dat
c                               2 keep geom.dat
c                               3 add fbeta files
c                               5 special magic code, crit&geom only
c                                 not paths.dat.  Use for path studies
c              ipr5 genfmt 0 files.dat, feff.dats that pass 2/3 of
c                            curved wave importance ratio
c                          1 keep all feff.dats
c              ipr6 ff2chi 0 chi.dat
c                          1 add sig2.dat with debye waller factors
c                          2 add chipnnnn.dat for each path
c                          3 add feffnnnn.dat for each path, and
c                            do not add chipnnnn.dat for each path
c                          4 add both feffnnnn.dat and chipnnnn.dat
c                            for each path
               if (nwords.eq.5) then
c                 feff7 input file
                  read(words(2),20,err=900)  ipr1
                  ipr2 = ipr1
                  ipr3 = ipr1
                  read(words(3),20,err=900)  ipr4
                  read(words(4),20,err=900)  ipr5
                  read(words(5),20,err=900)  ipr6
               else
c                 feff8 input file
                  read(words(2),20,err=900)  ipr1
                  read(words(3),20,err=900)  ipr2
                  read(words(4),20,err=900)  ipr3
                  read(words(5),20,err=900)  ipr4
                  read(words(6),20,err=900)  ipr5
                  read(words(7),20,err=900)  ipr6
               endif
               mode = 0
            elseif (itok .eq. 14)  then
c              POTENTIALS
c              Following lines are unique potential defs, 1 per line
               mode = 3
            elseif (itok .eq. 15)  then
c              NLEG nlegmax (for pathfinder)
               read(words(2),20,err=900)  nlegxx
               mode = 0
            elseif (itok .eq. 16)  then
c              CRIT critcw critpw
               read(words(2),30,err=900)  critcw
               read(words(3),30,err=900)  critpw
               mode = 0
            elseif (itok .eq. 17)  then
c              NOGEOM (do not write geom.dat) (disabled)
               nogeom = .true.
               mode = 0
            elseif (itok .eq. 18)  then
c              IORDER  iorder (used in genfmt, see setlam for meaning)
               read(words(2),20,err=900)  iorder
               call warnex(' IORDER:')
               mode = 0
            elseif (itok .eq. 19)  then
c              PCRIT  pcritk pcrith
c                     (keep and heap criteria for pathfinder)
               read(words(2),30,err=900)  pcritk
               read(words(3),30,err=900)  pcrith
               mode = 0
            elseif (itok .eq. 20)  then
c              SIG2  sig2g   global sig2 used by ff2chi, summed with
c              correlated debye model if DEBYE card used, and with
c              sig2 from files.dat if non-zero.
c              Units are Ang**2
               read(words(2),30,err=900)  sig2g
               mode = 0
            elseif (itok .eq. 21)  then
c              XANES ( xkmax  xkstep vixan)
               slog='XANES is not available in this version of feff8.5'
               call wlog(slog)
               call wlog("Ignoring XANES")
            elseif (itok .eq. 22)  then
c              CORRECTIONS  e0-shift, lambda correction
c              e0 shift is in eV, edge will be edge-e0
c              lambda corr is a const imag energy in eV
c              e0 and lambda corr same as vr0 and vi0 in EXCH card
               read(words(2),30,err=900)  vrcorr
               read(words(3),30,err=900)  vicorr
               mode = 0
            elseif (itok .eq. 23)  then
c              AFOLP use generalized automatic folp
               folpx = 1.15
               if (nwords.ge.2) read(words(2),30,err=900)  folpx
               mode =0
            elseif (itok .eq. 24)  then
c              EXAFS  xkmax for energy grid
               read(words(2),30,err=900)  xkmax
               mode = 0
            elseif (itok .eq. 25)  then
c              POLARIZATION  X Y Z
               ipol = 1
c              run linear polarization code 
               read(words(2),30,err=900)  evec(1)
               read(words(3),30,err=900)  evec(2)
               read(words(4),30,err=900)  evec(3)
               mode = 0
            elseif (itok .eq. 26)  then
c              ELLIPTICITY  E incident direction
               read(words(2),30,err=900)  elpty
               read(words(3),30,err=900)  xivec(1)
               read(words(4),30,err=900)  xivec(2)
               read(words(5),30,err=900)  xivec(3)
               mode = 0
            elseif (itok .eq. 27)  then
c              RGRID  rgrd
c              rgrd will be dpas, default is 0.03 in feff7
               read(words(2),30,err=900)  rgrd
               call warnex(' RGRID:')
               write(slog,'(a,1pe13.5)') ' RGRID, rgrd; ', rgrd
               call wlog(slog)
               i = 1 + int (12.5d0 / rgrd)
               if (mod(i,2) .eq. 0) i = i + 1
               if (i.gt.nrptx) then
                 write(slog,'(a,i6)') 
     1           ' FATAL error in RGRID: increase in dim.h nrptx to', i
                 call wlog(slog)
                 call par_stop(' ')
               endif
               mode = 0
            elseif (itok .eq. 28)  then
c              RPHASES (real phase shifts only)
               call warnex(' RPHASES:')
               call wlog(' Real phase shifts only will be used.  ' //
     1                   'FEFF results will be unreliable.')
               lreal = 2
               mode = 0
            elseif (itok .eq. 29)  then
c              NSTAR, write out n* for colinear polarization
               wnstar = .true.
               mode = 0
            elseif (itok .eq. 30)  then
c              NOHOLE
               if (nohole.lt.0) then
                  nohole = 0
                  if (nwords.ge.2) read(words(2),20,err=900)  nohole
                  call warnex(' NOHOLE:')
               end if
            elseif (itok .eq. 31)  then
c              SIG3 alphat  thetae   first and third cumulants for ss paths
               read(words(2),30,err=900)  alphat
               if (nwords.ge.3) read(words(3),20,err=900)  thetae
               call warnex(' SIG3:')
               write(slog,'(a,1pe13.5)') ' SIG3, alphat ; ', alphat
               call wlog(slog)
               mode = 0
            elseif (itok .eq. 32)  then
c              JUMPRM remove potential jumps at muffin tin radii
               jumprm = 1
            elseif (itok .eq. 33)  then
c              MBCONV do many body convolution with excitation spectrum
               mbconv = 1
            elseif (itok .eq. 34)  then
c              SPIN  specifies spin direction on central atom 
               read(words(2),20,err=900)  ispin 
c              set default spin along z axis
               if (ispin.ne.0) spvec(3) = 1.d0
               if (nwords.gt.2) read(words(3),30,err=900)  spvec(1)
               if (nwords.gt.3) read(words(4),30,err=900)  spvec(2)
               if (nwords.gt.4) read(words(5),30,err=900)  spvec(3)
            elseif (itok .eq. 35)  then
c              EDGE     L3 
c                   holecode
               call setedg (words(2), ihole)
               mode = 0
            elseif (itok .eq. 36)  then
c              SCF    rfms [ lfms nscmt  ca1 nmix  ecv icoul]
c              number of cycles, mode of calculating coulomb potential,
c              convergence accelerator
               nscmt = nbr
               ca1 = 0.2d0
               read(words(2),30,err=900)  rfms1
               if (nwords.gt.2) read(words(3),20,err=900)  lfms1
               if (nwords.gt.3) read(words(4),20,err=900)  nscmt
               if (nwords.gt.4) read(words(5),30,err=900)  ca1
               if (nwords.gt.5) read(words(6),20,err=900)  nmix
               if (nwords.gt.6) read(words(7),30,err=900)  ecv
               if (nwords.gt.7) read(words(8),20,err=900)  icoul
               if (nscmt.le.0 .or. nscmt.gt.nbr) nscmt = nbr
               if (lfms1.gt.0) lfms1 = 1
c              sanity checks for ca1
               if (ca1.lt.0) ca1 =0
               if (ca1.gt.0.5) then
                 call wlog(' Reduce convergence factors in SCF ')
                 call par_stop
     .            (' Cannot run with specified ca1 in SCF card.')
               endif
               if (ecv.ge.0) ecv = -40.0
               if (nmix.le.0) nmix=1
               if (nmix.gt.30) nmix=30
            elseif (itok .eq. 37)  then
c              FMS   rfms2  (lfms2 minv toler1 toler2 rdirec)
c              radius of the cluster to do FMS
               slog="FMS is not available in this version of FEFF8.5"
               call wlog(slog)
               call wlog("Ignoring FMS")
            elseif (itok .eq. 38)  then
               slog="LDOS is not available in this version of FEFF8.5"
               call wlog(slog)
               call wlog("Ignoring LDOS")
            elseif (itok .eq. 39)  then
c              INTERSTITIAL  inters  totvol
c              inters = 1 local V_int (around central atom)
c              inters = 0 extended V_int (average over all atoms)
c              more obscure options described in manual
               read(words(2),20,err=900)  inters
               if (nwords.ge.3) read(words(3),30,err=900)  totvol
            elseif (itok .eq. 40) then
c              CFAV  iphabs nabs rclabs
               read(words(2),20,err=900)  iphabs
               read(words(3),20,err=900)  nabs
               read(words(4),30,err=900)  rclabs
               if (rclabs.lt.0.5) rclabs=big
               mode = 0
            elseif (itok .eq. 41) then
c              S02  s02
               read(words(2),30,err=900)  s02
               mode = 0
            elseif (itok .eq. 42)  then
               slog="XES is not available in this version of FEFF8.5"
               call wlog(slog)
               call wlog("Ignoring XES")
            elseif (itok .eq. 43)  then
               slog="DANES is not available in this version of FEFF8.5"
               call wlog(slog)
               call wlog("Ignoring DANES")
c              DANES ( xkmax  xkstep vixan)
            elseif (itok .eq. 44)  then
c              FPRIME  emin emax estep
               slog="FPRIME is not available in this version of FEFF8.5"
               call wlog(slog)
               call wlog("Ignoring FPRIME")
            elseif (itok .eq. 45)  then
c              RSIGMA  (real self energy only)
               call warnex(' RSIGMA :')
               call wlog(' Real self energy only will be used.  ' //
     1                   'FEFF results will be unreliable.')
               if (lreal.lt.1) lreal = 1
               mode = 0
            elseif (itok .eq. 46)  then
c              XNCD or XMCD
               slog="XMCD is not available in this version of FEFF8.5"
               call wlog(slog)
               call wlog("Ignoring XMCD")
            elseif (itok .eq. 47)  then
c              MULTIPOLES le2 (l2lp)
               read(words(2),20,err=900)  le2
               if (nwords.gt.2) read(words(3),20,err=900)  l2lp
               mode = 0
            elseif (itok .eq. 48)  then
c              UNFREEZEF   
               iunf = 1
               mode = 0
            elseif (itok .eq. 49)  then
c              TDLDA 
               slog="TDLDA is not available in this version of FEFF8.5"
               call wlog(slog)
               call wlog("Ignoring TDLDA")
            elseif (itok .eq. 50)  then
c              PMBSE 
               slog="PMBSE is not available in this version of FEFF8.5"
               call wlog(slog)
               call wlog("Ignoring PMBSE")
            elseif (itok .eq. 51)  then ! Added by Josh Kas
c              PLASMON
               if(nwords.gt.1) then
                  read(words(2),20,err=900) iPlsmn
               else
                  iPlsmn = 1
               end if
            elseif (itok .eq. 52)  then ! Added by Josh Kas
c              S02CONV
               mso2conv = 1
            elseif (itok .eq. 53)  then ! Added by Josh Kas
c              SELF (print out on shell self energy Sig(k(E),E) )
               ipse = 1
            elseif (itok .eq. 54)  then ! Added by Josh Kas
c              SFSE k0 (print out self energy Sig(k0,E) ) 
               ipsk = 1
               read(words(2),30,err=900)  wsigk
            elseif (itok .eq. 55) then ! Added by Josh Kas
c              RCONV (print running convolution with file cfname at energy cen)
c              RCONV cen cname
               read(words(2),30,err=900) cen
               cfname = words(3)(1:12)
            elseif (itok.eq.56) then  !KJ added this card 1-06
c               ELNES
               slog="ELNES is not available in this version of FEFF8.5"
               call wlog(slog)
               call wlog("Ignoring ELNES")
            elseif (itok.eq.57) then  !KJ added this card 1-06
c               EXELFS
               slog="EXELFS is not available in this version of FEFF8.5"
               call wlog(slog)
               call wlog("Ignoring EXELFS")
            elseif (itok .eq. 58) then !KJ added this card 1-06
c               MAGIC card
               slog="MAGIC is not available in this version of FEFF8.5"
               call wlog(slog)
               call wlog("Ignoring MAGIC")
            elseif (itok .eq. 59) then !KJ added this card 3-06
c               ABSOLUTE card
               slog="ABSOLUTE is not available in this" // 
     &                  "version of FEFF8.5"
               call wlog(slog)
               call wlog("Ignoring ABSOLUTE")
            elseif ( itok .eq. 60) then
c               EGRID card
               slog="EGRID is not available in this version of FEFF8.5"
               call wlog(slog)
               call wlog("Ignoring EGRID")
            elseif (itok .eq. -1)  then
c              END
               goto 220
            else
               write(slog,'(1x,a)') line(1:70)
               call wlog(slog)
               write(slog,'(1x,a)') words(1)
               call wlog(slog)
               write(slog,'(a,i8)') ' Token ', itok
               call wlog(slog)
               call wlog(' Keyword unrecognized.')
               call wlog(' See FEFF document -- some old features')
               call wlog(' are no longer available.')
               call par_stop('RDINP-2')
            endif
         elseif (mode .eq. 1)  then
            if (itok .ne. 0)  then
c              We're done reading atoms.
c              Change mode and process current card.
               mode = 0
               goto 210
            endif
            natt = natt+1
            if (natt.gt. nattx)  then
               write(slog,'(a,i8)') 'Too many atoms, max is ', nattx
               call wlog(slog)
               call par_stop('RDINP-3')
            endif
            read(words(1),30,err=900)  ratx(1,natt)
            read(words(2),30,err=900)  ratx(2,natt)
            read(words(3),30,err=900)  ratx(3,natt)
            read(words(4),20,err=900)  iphatx(natt)
            if (iatph(iphatx(natt)) .le. 0) iatph(iphatx(natt)) = natt
         elseif (mode .eq. 2)  then
            if (itok .ne. 0)  then
c              We're done reading these overlap instructions.
c              Change mode and process current card.
               mode = 0
               goto 210
            endif
            novr(iph) = novr(iph)+1
            iovr = novr(iph)
            if (iovr .gt. novrx)  then
               write(slog,'(a,i8)') 'Too many overlap shells, max is ',
     1                               novrx
               call wlog(slog)
               call par_stop('RDINP-5')
            endif
            read(words(1),20,err=900) iphovr(iovr,iph)
            read(words(2),20,err=900) nnovr(iovr,iph)
            read(words(3),30,err=900) rovr(iovr,iph)
         elseif (mode .eq. 3)  then
            if (itok .ne. 0)  then
c              We're done reading unique potential definitions
c              Change mode and process current card.
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
               call par_stop('RDINP')
            endif
            read(words(2),20,err=900)  iz(iph)
            if (iz(iph).lt. 6) then
               lmaxsc(iph) = 1
            elseif (iz(iph).lt.55) then
               lmaxsc(iph) = 2
            else
               lmaxsc(iph) = 3
            endif
c           No potential label if user didn't give us one
c           Default set above is potlbl=' '
            if (nwords .ge. 3)  potlbl(iph) = words(3)(1:6)
            if (nwords .ge. 4)  then
              read(words(4),20,err=900) ltmp
              if (ltmp.ge.1 .and. ltmp.le.lx) lmaxsc(iph) = ltmp
            endif
            lmaxph(iph) = 3
            if (iz(iph).lt.6) lmaxph(iph) = 2
            if (nwords .ge. 5)  then
              read(words(5),20,err=900) ltmp
              if (ltmp.ge.1 .and. ltmp.le.lx) lmaxph(iph) = ltmp
            endif
            if (nwords .ge. 6) then
              read(words(6),30,err=900) xnatph(iph)
              lxnat = 1
            endif
            if (nwords .ge. 7) then
              read(words(7),30,err=900) spinph(iph)
            endif
           elseif (mode.eq.4) then  !KJ 1-06 this mode added to read ELNES card
             if(icnt.eq.5) then
                 read(words(1),30,err=900) ebeam   ! read beam energy in keV
                 ebeam=ebeam * dble(1000)  ! convert to eV
                 read(words(2),20,err=1011) aver ! average over sample to beam orientation?
                 read(words(3),20,err=1011) cross ! calculate cross terms?
                 read(words(4),20,err=1011) relat ! use relativistic q-vector?
                 read(words(5),20,err=1012) iinput ! read xmu.dat or opconsKK.dat or ... ?   !KJ 5/6
                 read(words(6),20,err=1013) spcol !column that has spectrum
                 if (aver.eq.1) icnt=icnt-1 !skip the line for beam orientation
                 goto 1011
 1012            iinput=1                
 1013            spcol=4
                 if(iinput.eq.2) spcol=3
 1011          continue ! Josh - Should these 1011 be 900? !KJ No.  Optional input.
              elseif(icnt.eq.4) then
                 read(words(1),30,err=900) xivec(1)  ! read direction of incoming beam
                 read(words(2),30,err=900) xivec(2)  ! in arbitrary units
                 read(words(3),30,err=900) xivec(3)
                 xinorm=dsqrt(xivec(1)**2+xivec(2)**2+xivec(3)**2)
                 if (xinorm.gt.0.0) then
                    do i=1,3
                       xivec(i)=xivec(i)/xinorm    ! normalize this vector.
                    enddo
                 elseif(.not.(aver.eq.1)) then
                    call wlog('WARNING : beam direction unspecified
     1                  in orientation sensitive EELS calculation.
     2                  Please correct before running EELS module.')
                 endif
             elseif(icnt.eq.3) then
                 read(words(1),30,err=900) acoll  ! collection semiangle in mrad
                 read(words(2),30,err=900) aconv  ! convergence semiangle in mrad
                 acoll=acoll/dble(1000);aconv=aconv/dble(1000) ! convert from mrad to rad
             elseif(icnt.eq.2) then
                 read(words(1),20,err=900) nqr    ! specify q-mesh, radial parameter
                 read(words(2),20,err=900) nqf    ! specify q-mesh, angular parameter
                 if(nqr*nqf.eq.0) then
                    call wlog('WARNING : zero q-mesh points specified
     1               for EELS calculation.  Please correct before
     2               running EELS module.')
                 endif
           elseif(icnt.eq.1) then
                 read(words(1),30,err=900) thetax ! detector position in plane perpendicular to beam ; angle in mrad
                 read(words(2),30,err=900) thetay ! detector position in plane perpendicular to beam ; angle in mrad
                 mode=0  ! finished reading ELNES card
!! initialize evec to be nonzero and perpendicular to xivec
!               if(dabs(xivec(1)-xivec(2)).gt.0.0001.or.
!     1            dabs(xivec(2)-xivec(3)).gt.0.0001) then
!                     evec(1)=xivec(2)*xivec(1)-xivec(3)**2
!                     evec(2)=xivec(3)*xivec(2)-xivec(1)**2
!                     evec(3)=xivec(1)*xivec(3)-xivec(2)**2
!                 else
!                     evec(1)=xivec(2)
!                     evec(2)=dble(0)
!                     evec(3)=-xivec(2)
!                 endif
             endif
             icnt=icnt-1    ! now read the next line
         !KJ end my changes                            
         else
            write(slog,'(a,i8)') 'Mode unrecognized, mode ', mode
            call wlog(slog)
            call par_stop('RDINP-6')
         endif
      goto 200
  220 continue
c done reading input file, 
c#{mn
c call rdline with jinit=0 to clean up all input files
       jinit = 0
       call rdline(jinit,line)
c#mn}

c     Fix up defaults, error check limits, figure out free atoms, etc.



c !KJ added this check 1-06
      if(magic.eq.1.and.(eels.ne.1)) then
          call wlog('To use MAGIC card you must have ELNES card.')
          call wlog('Ignoring MAGIC card.')
          magic=0
        endif
c !KJ

c  !KJ another check for eels 1-06
      if((eels.eq.1).and.(aver.eq.1).and.(cross.eq.1)) then
          call wlog('WARNING : you have asked to calculate an
     1   orientation averaged spectrum, but you have also asked
     2   to calculate cross-terms.  Averaging kills the cross terms.
     3   Hence the program ignores your request and does not
     4   calculate cross terms.')
      endif
c  !KJ

c  !KJ  set up a variable needed for elnes 1-06
        if(eels.eq.1) then
          if(aver.eq.1) then
             ipstep=1
             ipmin=10
             ipmax=10
          else
            ipmin=1
            ipmax=9
            if(cross.eq.1) then
               ipstep=1
            else
               ipstep=4
            endif
          endif
        endif
c  !KJ


c     need smaller rgrid for nonlocal exchange
      if (ixc0.lt.0) ixc0 = 0
      if (mod(ixc,10).ge.5 .and. rgrd.gt.0.03) rgrd=0.03d0 
      if (mod(ixc0,10).ge.5 .and. rgrd.gt.0.03) rgrd=0.03d0 
c     must use linear polarization to use nstar
      if (wnstar)  then
         if (ipol.ne.1)  then
            call wlog(' Must have linear polarization to use NSTAR.')
            call wlog(' NSTAR will be turned off.')
            wnstar = .false.
         endif
      endif

c     Do not use ihole .le. 0
      if (ihole .le. 0)  then
         call wlog(' Use NOHOLE to calculate without core hole.')
         call wlog(' Only ihole greater than zero are allowed.')
         call par_stop('RDINP')
      endif

c     Find out how many unique potentials we have
c     in POTENTIAL card
      nph = 0
      do 300  iph = nphx, 0, -1
         if (iz(iph) .gt. 0)  then
            nph = iph
            goto 301
         endif
  300 continue
  301 continue

c     cannot use OVERLAP and ATOMS cards together
      if (iatom .gt. 0 .and. iovrlp .gt. 0)  then
        call wlog(' Cannot use ATOMS and OVERLAP in the same feff.inp.')
        call par_stop('RDINP')
      endif

c     cannot use OVERLAP and CFAVERAGE   cards together
      if (novr(0) .gt. 0) then
c        OVERLAP is used, cannot do configuration average
         iphabs = 0
         nabs = 1
         rclabs = big
      endif

c     Must have central atom
      if (iz(0) .le. 0)  then
         if (iphabs .gt. 0) then
c           central atom is of the iphabs type
            iz(0) = iz(iphabs)
            potlbl(0) = potlbl(iphabs)
            lmaxsc(0) = lmaxsc(iphabs)
            lmaxph(0) = lmaxph(iphabs)
            xion(0) = xion(iphabs)
         else
            call wlog(' No absorbing atom (unique pot 0) was defined.')
            call par_stop('RDINP')
         endif
      endif

c     No gaps allowed in unique pots.  Make sure we have enough
c     to overlap all unique pots 0 to nph.
      if (iphabs.gt.0 .and. iatph(0).le.0)   iatph(0) = iatph(iphabs)
      do 340  iph = 0, nph
         if (iatph(iph) .le. 0  .and.  novr(iph) .le. 0)  then
c           No model atom, no overlap cards, can't do this unique pot
            write(slog,'(a,i8)') 
     1       ' No atoms or overlap cards for unique pot ', iph
            call wlog(slog)
            call wlog(' Cannot calculate potentials, etc.')
            call par_stop('RDINP-')
         endif
c        by default freeze f-electrons and reset lmaxsc=2
         if (iunf.eq.0 .and. lmaxsc(iph).gt.2) lmaxsc(iph)=2
  340 continue

c     Need number of atoms of each unique pot, count them.  If none,
c     set to one. Do statistics for all atoms in feff.inp.
      do 350  iph = 0, nph
        if (lxnat.eq.0) then 
          xnatph(iph) = 0
          do 346  iat = 1, natt
              if (iphatx(iat) .eq. iph)  xnatph(iph) = xnatph(iph)+1
  346     continue
          if (iph.gt.0 .and. iph.eq.iphabs) xnatph(iph) = xnatph(iph)-1
        else
          if (xnatph(iph).le. 0.01) then
            if (iph.eq.0) then
              xnatph(iph) = 0.01d0
            else
              write (slog,'(a,i4)') ' Inconsistency in POTENTIAL card'//
     1                              ' is detected for unique pot ', iph
              call wlog (slog)
              call wlog (' Results might be meaningless.')
            endif
          endif
        endif
        if (xnatph(iph) .le. 0)  xnatph(iph) = 1
  350 continue
      if (lxnat.ne.0) then
c        normalize statistics to hav one absorber
         do 351 iph = 1, nph
  351    xnatph(iph) = xnatph(iph) /xnatph(0)
         xnatph(0) = 1
      endif
      xnat = 0
      do 352 iph = 0,nph
  352 xnat = xnat + xnatph(iph)

c     Find distance to nearest and most distant atom (use overlap card
c     if no atoms specified.)
      if (natt .lt. 2)  then
         ratmin = rovr(1,0)
         ratmax = rovr(novr(0),0)
      else
         ratmax = 0
         ratmin = 1.0e10
         iatabs = iatph(0)
         icount = 0
         if (iatabs.le.0) iatabs = iatph( iphabs)
         if (iatabs.le.0) call par_stop('RDINP fatal error: iatabs=NaN')

         do 412  iat = 1, natt
           if (iphatx(iat) .eq. iphabs .or. iphatx(iat).eq.0)
     1        icount = icount +1
           if (iat.ne.iatabs) then
c           skip absorbing atom
            tmp = dist (ratx(1,iat), ratx(1,iatabs))
            if (tmp .gt. ratmax)  ratmax = tmp
            if (tmp .lt. ratmin)  ratmin = tmp
           endif
  412    continue
         if (nabs.le.0) nabs = icount
      endif

c     Set total volume
      if (totvol.gt.0) totvol = totvol * ratmin**3 * xnat

c     Set rfms if they are too small
      if (rfms1 .lt. ratmin) rfms1 = -1.e0
      if (rfms2 .lt. ratmin) rfms2 = -1.e0
      if (rfms2 .lt. ratmin .and. ispec.lt.2) ispec = - ispec 
      if (rfms2 .lt. ratmin .and. ispec.eq.3) ispec = - ispec 
c     if ispec.le.0 MS expansion will be used, else - FMS method.
      

c     Set rmax if necessary
      if (rmax.le.0 .and. nss.le.0 .and. ispec.le.0)  then
c        set to min (2+ times ratmin, ratmax) (magic numbers to
c        avoid roundoff, note that rmax is single precision).
         rmax = min (2.2 * real(ratmin), 1.01 * real(ratmax))
      endif

c     Set core hole lifetime (central atom quantity) and s02
      iph = 0
      call setgam (iz(iph), ihole, gamach)
      if (s02.eq.1.d0) s02=s02h

c     Convert everything to code units, and use rmult factor
c     rmax is for pathfinder, so leave it in Ang.
      rmax = rmax * real(rmult)
      rfms1 = rfms1 * real(rmult)
      rfms2 = rfms2 * rmult 
      totvol = totvol * rmult**3
c     Use rmult factor.  Leave distances in Ang.
      do 430  iat = 1, natt
         do 420  i = 1, 3
            ratx(i,iat) = ratx(i,iat) * rmult
  420    continue
  430 continue
      do 460  iph = 0, nph
         do 450  iovr = 1, novr(iph)
            rovr(iovr,iph) = rovr(iovr,iph) * rmult
  450    continue
  460 continue
      do 462  iss = 1, nss
c        rss used only to make paths.dat, so leave it in Angstroms.
         rss(iss) = rss(iss) * rmult
  462 continue

c     Clean up control flags
      if (mpot .ne. 0)  mpot = 1
      if (mphase .ne. 0)  mphase = 1
      if (mfms .ne. 0)  mfms = 1
      if (mpath  .ne. 0)  mpath = 1
      if (mfeff  .ne. 0)  mfeff = 1
      if (mchi   .ne. 0)  mchi = 1
      if (nss    .le. 0)  ms = 1
      if (ifolp  .ne. 0)  iafolp = -1
      if (natt.le.0) then
c       Overalp geometry
        mfms = 0
        mpath = 0
        ms = 0
c       no SCF loop
        nscmt = 0
        do 464 iph = 0, nph
          if (novr(iph).le.0) call par_stop('Bad OVERLAP cards.')
  464   continue
      endif

      if (iafolp .ge. 0) then
         do 485 i = 0, nphx
  485    folp(i) = folpx
      endif

      if (ntitle .le. 0)  then
         ntitle = 1
         title(1) = 'Null title'
      endif
      do 490  i = 1, ntitle
         ltit(i) = istrln (title(i))
  490 continue
      nttl = ntitle

c     write atoms.dat, global.inp, modN.inp and ldos.inp
c      call wrtall (nabs)
      call wrtjsn (nabs)

c     In case of OVERLAP and SS calculateions write 'paths.dat'
c     without invoking the pathfinder. Single scattering paths only.
      if (nss .gt. 0  .and.  mpath .eq. 1)  then
         open (unit=1, file='paths.dat', status='unknown', iostat=ios)
         call chopen (ios, 'paths.dat', 'rdinp')
         do 750  i = 1, ntitle
            write(1,748)  title(i)(1:ltit(i))
  748       format (1x, a)
  750    continue
         write(1,751)
  751    format (' Single scattering paths from ss lines cards',
     1           ' in feff input')
         write(1,706)
  706    format (1x, 71('-'))
         do 760  iss = 1, nss
            if (rmax.le.0  .or.  rss(iss).le.rmax)  then
c              NB, rmax and rss are in angstroms
               write(1,752) indss(iss), 2, degss(iss),
     2              rss(iss)
  752          format ( 2i4, f8.3,
     1             '  index,nleg,degeneracy,r=', f8.4)
               write(1,766)
  766          format (' single scattering')
               write(1,754) rss(iss), zero, zero, iphss(iss),
     1                      potlbl(iphss(iss))
               write(1,753) zero, zero, zero, 0, potlbl(0)
  753          format (3f12.6, i4,  1x, '''', a6, '''', '  x,y,z,ipot')
  754          format (3f12.6, i4,  1x, '''', a6, '''')
            endif
  760    continue
         close (unit=1)
      endif

      do 120  i = 1, ntitle
         call wlog(' ' // title(i)(1:ltit(i)))
  120 continue

c     if user does not want geom.dat, don't do it
      if (nogeom)  then
c        don't delete geom.dat when done with it either...
         if (ipr4 .lt. 2)  ipr4 = 2
         if (nabs.gt.1) call 
     1     par_stop('NOGEOM and CFAVERAGE are incompatible')
      else
c       temporarily call ffsort. here
        iabs = 1
c !KJ 1-06 : If the user does EELS and does not calculate cross terms for an
c       orientation sensitive calculation, FEFF mustn't change the
c       coordinate system, as this would lead to the appearance of
c       cross terms after all.  Therefore, I added an argument to the
c       calling sequence of ffsort.
c       To be precise, giving '.false.' disables the call of ffsort to mkptz.
c       Giving '.true.' makes ffsort work exactly as it always has.
        if((eels.eq.1)) then
           call ffsort(iabs,.false.)
        else
           call ffsort(iabs,.true.)
        endif   !KJ end my changes
       endif
       ceels=(eels.eq.1) !KJ 5-6 for monolithic version

      close(unit=11)
  400 call par_barrier
      call par_end

c     sub-program exchange
      stop
c     return

c     normal end of rdinp

  900 continue
      call wlog(' Error reading input, bad line follows:')
      write(slog,'(1x,a)') line(1:71)
      call wlog(slog)
      call par_stop('RDINP fatal error.')

      end

      subroutine phstop (iph,line)
      implicit double precision (a-h, o-z)
      character*(*) line
      include '../HEADERS/dim.h'
      character*512 slog
      if (iph .lt. 0  .or.  iph .gt. nphx)  then
         write(slog,10) iph, nphx, line
         call wlog(slog)
   10    format (' Unique potential index', i5, ' out of range.', 
     1           ' Must be between 0 and', i5, '.  Input line:', 
     2           1x, a)
         call par_stop('RDINP - PHSTOP')
      endif
      return
      end

      subroutine warnex (string)
      implicit double precision (a-h, o-z)
c     This prints a warning message if the user is using an
c     expert option.
      character*(*) string

      call wlog(string)
      call wlog(' Expert option, please read documentation ' //
     1          'carefully and check your results.')
      return
      end
