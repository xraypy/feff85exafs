      subroutine json_read_libpotph(
c     TITLE
     1       ntitle, title,
c     ATOMS
     2       nat, rat, iphat,
c     POTENTIALS
     3       nph, iz, potlbl, lmaxsc, lmaxph, xnatph, spinph,
c     HOLE/EDGE
     4       ihole,
c     SCF
     5       rfms1, lfms1, nscmt, ca1, nmix, ecv, icoul,
c     POLARIZATION, ELLIPTICITY
     6       ipol, evec, elpty, xivec,
c     SPIN
     7       ispin, spvec, angks,
c     computed
     8       ptz, gamach,
c     EXCHANGE
     9       ixc, vr0, vi0, ixc0,
c     AFOLP, FOLP, ION, RGRID, UNFREEZEF
     _       iafolp, folp, xion, rgrd, iunf,
c     INTERSTITIAL, JUMPRM, NOHOLE
     1       inters, totvol, jumprm, nohole, iplsmn)




c$$$(nat, rat, iphat,
c$$$     2       le2, elpty, angks, evec, xivec, spvec, ptz,
c$$$     1       nabs, iphabs, rclabs, ipol, ispin,
c$$$     5       mpot, rgrd, ntitle, title, ipr1, ispec,
c$$$     6       nohole, ihole, gamach, nph, iz, lmaxsc, xnatph,
c$$$     7       xion, iunf, ixc, jumprm, iafolp, folp, inters, totvol,
c$$$     8       rfms1, lfms1, nscmt, ca1, nmix, ecv, icoul,
c$$$     3       mphase, ipr2, vixan, xkstep, xkmax,
c$$$     7       lmaxph, potlbl, spinph, vr0, vi0, ixc0, lreal, 
c$$$     1       rfms2, lfms2, l2lp, iPl, iGrid,
c$$$     2       izstd, ifxc, ipmbse, itdlda, nonlocal, ibasis)

      use json_module
      implicit double precision (a-h, o-z)

      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'

      logical :: found
      type(json_file) :: json   !the JSON structure read from the file:
      integer,dimension(:),allocatable :: intgs
      character*80,dimension(:),allocatable :: strings
      double precision,dimension(:),allocatable :: dbpcx, dbpcy, dbpcz
      double precision,dimension(:),allocatable :: dbpc1, dbpc2, dbpc3
      double precision,dimension(:),allocatable :: dbpc4, dbpc5, dbpc6
      double precision,dimension(:),allocatable :: dbpc7, dbpc8, dbpc9
      double precision,dimension(:),allocatable :: dbpc10

c     dimension/types of atoms.json things
c     dimension/types of global.json things
      integer nat, iphat(natx), ipol, ispin
      double precision evec(3), xivec(3), spvec(3), 
     1       rat(3,natx), elpty, angks
      complex*16 ptz(-1:1, -1:1)


c     dimension/type os mod1/pot things
      character*80 title(nheadx)
      integer nph, ntitle, ihole, iafolp, ixc, iplsmn,
     1       iunf, nmix, nohole, jumprm, inters, nscmt, icoul, lfms1
      integer iz(0:nphx), lmaxsc(0:nphx)
      real rfms1
      double precision gamach, rgrd, ca1, ecv, totvol
      double precision  xnatph(0:nphx), folp(0:nphx), xion(0:nphx)
c     for OVERLAP option
c$$ovr$$      integer novr(0:nphx), iphovr(novrx,0:nphx), nnovr(novrx,0:nphx)
c$$ovr$$      double precision  rovr(novrx,0:nphx)


      double precision vr0, vi0, spinph(0:nphx)
      integer lmaxph(0:nphx)
      character*6  potlbl(0:nphx)

c      integer iphabs, ipr1, le2, nabs
c      integer mpot, mphase, ipr2, ixc0, ispec, lreal, lfms2, l2lp, iPl, 
c     1       iGrid
c     double precision xkstep, xkmax, vixan, rclabs
c     real rfms2
c     integer izstd, ifxc, ipmbse, itdlda, nonlocal, ibasis

      character*32 s1, s2, s3

      call json%load_file('libpotph.json')
      if (json_failed()) then   !if there was an error reading the file
         print *, "failed to read libpotph.json"
         stop
      else

c        TITLE
         call json%get('ntitle', ntitle, found)
                if (.not. found) call bailout('ntitle', 'libpotph.json')
         call json%get('titles', strings, found)
                if (.not. found) call bailout('titles', 'libpotph.json')
         do 10 itit = 1, nheadx
            title(itit) = strings(itit)            
 10      continue


c        ATOMS
         call json%get('natt', nat, found)
                if (.not. found) call bailout('natt',   'libpotph.json')
         call json%get('x',    dbpcx, found)
                if (.not. found) call bailout('x',      'libpotph.json')
         call json%get('y',    dbpcy, found)
                if (.not. found) call bailout('y',      'libpotph.json')
         call json%get('z',    dbpcz, found)
                if (.not. found) call bailout('z',      'libpotph.json')
         call json%get('iphatx', intgs, found)
                if (.not. found) call bailout('iphatx', 'libpotph.json')
         do 20 i=1,nat
            iphat(i)  = intgs(i)
            rat(1,i)  = dbpcx(i)
            rat(2,i)  = dbpcy(i)
            rat(3,i)  = dbpcz(i)
 20      continue

c        POTENTIALS
         call json%get('nph',    nph, found)
                if (.not. found) call bailout('nph',    'libpotph.json')
         call json%get('iz', intgs, found)
                if (.not. found) call bailout('iz',     'libpotph.json')
         do 30 iph = 0, nphx
            iz(iph) = intgs(iph+1)            
 30      continue
         call json%get('potlbl', strings, found)
                if (.not. found) call bailout('potlbl', 'libpotph.json')
         do 40 itit = 1, nphx
c            potlbl(itit-1) = strings(itit)
            potlbl(itit-1) = strings(itit)(1:6)
 40      continue
         call json%get('lmaxsc', intgs, found)
                if (.not. found) call bailout('lmaxsc', 'libpotph.json')
         do 50 iph = 0, nphx
            lmaxsc(iph) = intgs(iph+1)            
 50      continue
         call json%get('lmaxph', intgs, found)
                if (.not. found) call bailout('lmaxph', 'libpotph.json')
         do 60 iph = 0, nphx
            lmaxph(iph) = intgs(iph+1)            
 60      continue
         call json%get('xnatph', dbpc1, found)
                if (.not. found) call bailout('xnatph', 'libpotph.json')
         do 70 iph = 0, nphx
            xnatph(iph) = dbpc1(iph+1)            
 70      continue
         call json%get('spinph', dbpc2, found)
                if (.not. found) call bailout('spinph', 'libpotph.json')
         do 80 iph = 0, nphx
            spinph(iph) = dbpc2(iph+1)            
 80      continue

c        HOLE/EDGE
         call json%get('ihole',  ihole, found)
                if (.not. found) call bailout('ihole',  'libpotph.json')

c        SCF
         call json%get('rfms1',  toss, found)   
                if (.not. found) call bailout('rfms1',  'libpotph.json')
         rfms1 = real(toss)
         call json%get('lfms1',  lfms1, found)
                if (.not. found) call bailout('lfms1',  'libpotph.json')
         call json%get('nscmt',  nscmt, found)
                if (.not. found) call bailout('nscmt',  'libpotph.json')
         call json%get('ca1',    ca1, found)
                if (.not. found) call bailout('ca1',    'libpotph.json')
         call json%get('nmix',   nmix, found)
                if (.not. found) call bailout('nmix',   'libpotph.json')
         call json%get('ecv',    ecv, found)
                if (.not. found) call bailout('ecv',    'libpotph.json')
         call json%get('icoul',  icoul, found)
                if (.not. found) call bailout('icoul',  'libpotph.json')

c        POLARIZATION, ELLIPTICITY
         call json%get('ipol', ipol, found)
                if (.not. found) call bailout('ipol',   'libpotph.json')
         call json%get('evec',  dbpc3, found)
                if (.not. found) call bailout('evec',   'libpotph.json')
         do 90 i=1,3
            evec(i) = dbpc3(i)
 90      continue
         call json%get('elpty', elpty, found)
                if (.not. found) call bailout('elpty',  'libpotph.json')
         call json%get('xivec',  dbpc4, found)
                if (.not. found) call bailout('xivec',  'libpotph.json')
         do 100 i=1,3
            xivec(i) = dbpc4(i)
 100     continue

c        SPIN
         call json%get('ispin', ispin, found)
                if (.not. found) call bailout('ispin',  'libpotph.json')
         call json%get('spvec',  dbpc5, found)
                if (.not. found) call bailout('spvec',  'libpotph.json')
         do 110 i=1,3
            spvec(i) = dbpc5(i)
 110     continue
         call json%get('angks', angks, found)
                if (.not. found) call bailout('angks',  'libpotph.json')


c        computed: ptz and gamach
         call json%get('ptz0',  dbpc6, found)
                if (.not. found) call bailout('ptz0',  'libpotph.json')
         ptz(-1,-1) = dcmplx(dbpc6(1), dbpc6(2))
         ptz( 0,-1) = dcmplx(dbpc6(3), dbpc6(4))
         ptz( 1,-1) = dcmplx(dbpc6(5), dbpc6(6))
         call json%get('ptz1',  dbpc7, found)
                if (.not. found) call bailout('ptz1',  'libpotph.json')
         ptz(-1, 0) = dcmplx(dbpc7(1), dbpc7(2))
         ptz( 0, 0) = dcmplx(dbpc7(3), dbpc7(4))
         ptz( 1, 0) = dcmplx(dbpc7(5), dbpc7(6))
         call json%get('ptz2',  dbpc8, found)
                if (.not. found) call bailout('ptz2',  'libpotph.json')
         ptz(-1, 1) = dcmplx(dbpc8(1), dbpc8(2))
         ptz( 0, 1) = dcmplx(dbpc8(3), dbpc8(4))
         ptz( 1, 1) = dcmplx(dbpc8(5), dbpc8(6))
         call json%get('gamach', gamach, found)
                if (.not. found) call bailout('gamach', 'libpotph.json')


c        EXCHANGE
         call json%get('ixc',    ixc, found)
                if (.not. found) call bailout('ixc',    'libpotph.json')
         call json%get('vro',   vr0, found)
                if (.not. found) call bailout('vr0',    'libpotph.json')
         call json%get('vio',   vi0, found)
                if (.not. found) call bailout('vi0',    'libpotph.json')
         call json%get('ixc0',   ixc0, found)
                if (.not. found) call bailout('ixc0',   'libpotph.json')

c        AFOLP, FOLP, ION, RGRID, UNFREEZEF
         call json%get('iafolp', iafolp, found)
                if (.not. found) call bailout('iafolp', 'libpotph.json')
         call json%get('folp', dbpc9, found)
                if (.not. found) call bailout('folp',   'libpotph.json')
         do 120 iph = 0, nphx
            folp(iph) = dbpc9(iph+1)            
 120    continue
         call json%get('xion', dbpc10, found)
                if (.not. found) call bailout('xion',   'libpotph.json')
         do 130 iph = 0, nphx
            xion(iph) = dbpc10(iph+1)            
 130    continue
         call json%get('rgrd',   rgrd, found)
                if (.not. found) call bailout('rgrd',   'libpotph.json')
         call json%get('iunf',   iunf, found)
                if (.not. found) call bailout('iunf',   'libpotph.json')
         

c        INTERSTITIAL, JUMPRM, NOHOLE, PLASMON
         call json%get('inters', inters, found)
                if (.not. found) call bailout('inters', 'libpotph.json')
         call json%get('totvol', totvol, found)
                if (.not. found) call bailout('totvol', 'libpotph.json')
         call json%get('nohole', nohole, found)
                if (.not. found) call bailout('nohole', 'libpotph.json')
         call json%get('jumprm', jumprm, found)
                if (.not. found) call bailout('jumprm', 'libpotph.json')
c         call json%get('iplsmn', iplsmn, found)
c                if (.not. found) call bailout('iplsmn', 'libpotph.json')


c The rest are parameters associated with:
c    CFAVERAGE, MULTIPOLE, XANES, FMS, PMBSE, TDLDA, RPHASES
c    CONTROL, PRINT 

c        content of global.json
c         call json%get('nabs', nabs, found)
c                if (.not. found) call bailout('nabs',   'libpotph.json')
c         call json%get('iphabs', iphabs, found)
c                if (.not. found) call bailout('iphabs', 'libpotph.json')
c         call json%get('rclabs', rclabs, found)
c                if (.not. found) call bailout('rclabs', 'libpotph.json')
c         call json%get('le2', le2, found)
c                if (.not. found) call bailout('le2',    'libpotph.json')



c        content of pot/mod1.json (except overlap stuff)
c         call json%get('mpot',   mpot, found)
c                if (.not. found) call bailout('mpot',   'libpotph.json')
c         call json%get('ipr1',   ipr1, found)
c                if (.not. found) call bailout('ipr1',   'libpotph.json')
c         call json%get('ispec',  ispec, found)
c                if (.not. found) call bailout('ispec',  'libpotph.json')




c        content of xsph/mod2.json
c         call json%get('mphase',   mphase, found)
c                if (.not. found) call bailout('mphase', 'libpotph.json')
c         call json%get('ipr2',   ipr2, found)
c                if (.not. found) call bailout('ipr2',   'libpotph.json')
c         call json%get('lreal',   lreal, found)
c                if (.not. found) call bailout('lreal',  'libpotph.json')
c         call json%get('lfms2',   lfms2, found)
c                if (.not. found) call bailout('lfms2',  'libpotph.json')
c         call json%get('l2lp',   l2lp, found)
c                if (.not. found) call bailout('l2lp',   'libpotph.json')
c         call json%get('iGrid',   iGrid, found)
c                if (.not. found) call bailout('iGrid',  'libpotph.json')


c         call json%get('rfms2',  toss, found)
c                if (.not. found) call bailout('rfms2',  'libpotph.json')
c         rfms2 = real(toss)
c         call json%get('xkstep',   xkstep, found)
c                if (.not. found) call bailout('xkstep', 'libpotph.json')
c         call json%get('xkmax',   xkmax, found)
c                if (.not. found) call bailout('xkmax',  'libpotph.json')
c         call json%get('vixan',   vixan, found)
c                if (.not. found) call bailout('vixan',  'libpotph.json')

c         call json%get('izstd',   izstd, found)
c                if (.not. found) call bailout('izstd')
c         call json%get('ifxc',   ifxc, found)
c                if (.not. found) call bailout('ifxc')
c         call json%get('ipmbse',   ipmbse, found)
c                if (.not. found) call bailout('ipmbse')
c         call json%get('itdlda',   itdlda, found)
c                if (.not. found) call bailout('itdlda')
c         call json%get('nonlocal',   nonlocal, found)
c                if (.not. found) call bailout('nonlocal')
c         call json%get('ibasis',   ibasis, found)
c                if (.not. found) call bailout('ibasis')


         call json%destroy()
      end if

c********************************************************************************
c     the 3000s are taken straight from reapot.f
c     transform to code units (bohrs and hartrees - atomic unuts)
      rfms1 = rfms1 / real(bohr)
      gamach = gamach / hart
      ecv   = ecv   / hart
      totvol = totvol / bohr**3
      do 3010 iat = 1, nat
         do 3000 i = 1,3
            rat(i,iat) = rat (i, iat) / bohr
 3000    continue
 3010 continue

c$$ovr$$      do 3030 iph = 0, nph
c$$ovr$$         do 3020 iovr = 1, novr(iph)
c$$ovr$$            rovr(iovr,iph) = rovr(iovr,iph) / bohr
c$$ovr$$ 3020    continue
c$$ovr$$ 3030 continue

c     add lines to the title
c      if (mpot.eq.1) then
      ntitle = ntitle + 1
      if (nat.gt.1) then
         if (rfms1.lt.0) rfms1 = 0
         if (nscmt.gt.0) then
            write(s1, 3230) nscmt, rfms1*bohr, lfms1
 3230       format(' POT  SCF', i4, f8.4, i4)
         else
            write(s1, 3235) 
 3235       format(' POT  Non-SCF' )
         endif
      else
         write(s1, 3240) 
 3240    format(' POT  used OVERLAP geometry,')
      endif
      if (nohole.eq.0) then
         write(s2, 3310) 
 3310    format(', NO core-hole,')
      elseif (nohole.eq.2) then
         write(s2, 3315) 
 3315    format(', screened core-hole,')
      else
         write(s2, 3320) 
 3320    format(', core-hole,')
      endif
      if (iafolp.lt.0) then
         write(s3, 3330) folp(0)
 3330    format(' FOLP (folp(0)=', f6.3, ')' )
      else
         write(s3, 3340) folp(0)
 3340    format(' AFOLP (folp(0)=', f6.3, ')' )
      endif
c     concatenate 3 strings into 1
      title(ntitle) = ' '
      ilen   = istrln(s1)
      istart = 1
      iend   = ilen
      title(ntitle)(istart:iend) = s1(1:ilen)
      ilen   = istrln(s2)
      istart = iend + 1
      iend   = iend + ilen
      title(ntitle)(istart:iend) = s2(1:ilen)
      ilen   = istrln(s3)
      istart = iend + 1
      iend   = iend + ilen
      title(ntitle)(istart:iend) = s3(1:ilen)
c     endif


      return
      end
