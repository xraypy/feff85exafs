      subroutine json_read_libpotph(nat, rat, iphat,
     2       le2, elpty, angks, evec, xivec, spvec, ptz,
     1       nabs, iphabs, rclabs, ipol, ispin,
     5       mpot, rgrd, ntitle, title, ipr1, ispec,
     6       nohole, ihole, gamach, nph, iz, lmaxsc, xnatph,
     7       xion, iunf, ixc, jumprm, iafolp, folp, inters, totvol,
     8       rfms1, lfms1, nscmt, ca1, nmix, ecv, icoul,
     3       mphase, ipr2, vixan, xkstep, xkmax,
     7       lmaxph, potlbl, spinph, vr0, vi0, ixc0, lreal, 
     1       rfms2, lfms2, l2lp, iPl, iGrid,
     2       izstd, ifxc, ipmbse, itdlda, nonlocal, ibasis)

      use json_module
      implicit double precision (a-h, o-z)

      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'

      logical :: found
      type(json_file) :: json   !the JSON structure read from the file:
      integer,dimension(:),allocatable :: intgs
      character*80,dimension(:),allocatable :: strings
      double precision,dimension(:),allocatable :: dbpcs, dbpcy, dbpcz

c     dimension/types of atoms.json things
c     dimension/types of global.json things
      integer nat, iphat(natx), nabs, iphabs, ipol, ispin, le2
      double precision evec(3), xivec(3), spvec(3), 
     1       rat(3,natx), elpty, angks, rclabs
      complex*16 ptz(-1:1, -1:1)


c     dimension/type os mod1/pot things
      character*80 title(nheadx)
      integer mpot, nph, ntitle, ihole, ipr1, iafolp, ixc, 
     1       iunf, nmix, nohole, jumprm, inters, nscmt, icoul, lfms1
      integer iz(0:nphx), lmaxsc(0:nphx)
      real rfms1
      double precision gamach, rgrd, ca1, ecv, totvol
      double precision  xnatph(0:nphx), folp(0:nphx), xion(0:nphx)
c     for OVERLAP option
c      integer novr(0:nphx), iphovr(novrx,0:nphx), nnovr(novrx,0:nphx)
c      double precision  rovr(novrx,0:nphx)


c     dimension/type os mod2/xpsh things
      integer mphase, ipr2, ixc0, ispec, lreal, lfms2, l2lp, iPl, 
     1       iGrid
      double precision xkstep, xkmax, vixan
      double precision vr0, vi0, spinph(0:nphx)
      real rfms2
      integer lmaxph(0:nphx)
      character*6  potlbl(0:nphx)
      integer izstd, ifxc, ipmbse, itdlda, nonlocal, ibasis


      call json%load_file('libpotph.json')
      if (json_failed()) then   !if there was an error reading the file
         print *, "failed to read libpotph.json"
         stop
      else

c        content of atoms.json
         call json%get('natt', nat, found)
                if (.not. found) call bailout('natt',   'libpotph.json')
         call json%get('x',    dbpcs, found)
                if (.not. found) call bailout('x',      'libpotph.json')
         call json%get('y',    dbpcy, found)
                if (.not. found) call bailout('y',      'libpotph.json')
         call json%get('z',    dbpcz, found)
                if (.not. found) call bailout('z',      'libpotph.json')
         call json%get('iphatx', intgs, found)
                if (.not. found) call bailout('iphatx', 'libpotph.json')
         do 10 i=1,nat
            iphat(i)  = intgs(i)
            rat(1,i)  = dbpcs(i)
            rat(2,i)  = dbpcy(i)
            rat(3,i)  = dbpcz(i)
 10      continue


c        content of global.json
         call json%get('nabs', nabs, found)
                if (.not. found) call bailout('nabs',   'libpotph.json')
         call json%get('iphabs', iphabs, found)
                if (.not. found) call bailout('iphabs', 'libpotph.json')
         call json%get('rclabs', rclabs, found)
                if (.not. found) call bailout('rclabs', 'libpotph.json')
         call json%get('ipol', ipol, found)
                if (.not. found) call bailout('ipol',   'libpotph.json')
         call json%get('ispin', ispin, found)
                if (.not. found) call bailout('ispin',  'libpotph.json')
         call json%get('le2', le2, found)
                if (.not. found) call bailout('le2',    'libpotph.json')
         call json%get('elpty', elpty, found)
                if (.not. found) call bailout('elpty',  'libpotph.json')
         call json%get('angks', angks, found)
                if (.not. found) call bailout('angks',  'libpotph.json')

         call json%get('evec',  dbpcs, found)
                if (.not. found) call bailout('evec',   'libpotph.json')
         do 110 i=1,3
            evec(i) = dbpcs(i)
 110     continue
         call json%get('xivec',  dbpcs, found)
                if (.not. found) call bailout('xivec',  'libpotph.json')
         do 120 i=1,3
            xivec(i) = dbpcs(i)
 120     continue
         call json%get('spvec',  dbpcs, found)
                if (.not. found) call bailout('spvec',  'libpotph.json')
         do 130 i=1,3
            spvec(i) = dbpcs(i)
 130     continue

         call json%get('ptz0',  dbpcs, found)
                if (.not. found) call bailout('ptz0',  'libpotph.json')
         ptz(-1,-1) = dcmplx(dbpcs(1), dbpcs(2))
         ptz( 0,-1) = dcmplx(dbpcs(3), dbpcs(4))
         ptz( 1,-1) = dcmplx(dbpcs(5), dbpcs(6))
         call json%get('ptz1',  dbpcs, found)
                if (.not. found) call bailout('ptz1',  'libpotph.json')
         ptz(-1, 0) = dcmplx(dbpcs(1), dbpcs(2))
         ptz( 0, 0) = dcmplx(dbpcs(3), dbpcs(4))
         ptz( 1, 0) = dcmplx(dbpcs(5), dbpcs(6))
         call json%get('ptz2',  dbpcs, found)
                if (.not. found) call bailout('ptz2',  'libpotph.json')
         ptz(-1, 1) = dcmplx(dbpcs(1), dbpcs(2))
         ptz( 0, 1) = dcmplx(dbpcs(3), dbpcs(4))
         ptz( 1, 1) = dcmplx(dbpcs(5), dbpcs(6))


c        content of pot/mod1.json (except overlap stuff)
         call json%get('mpot',   mpot, found)
                if (.not. found) call bailout('mpot',   'libpotph.json')
         call json%get('nph',    nph, found)
                if (.not. found) call bailout('nph',    'libpotph.json')
         call json%get('ntitle', ntitle, found)
                if (.not. found) call bailout('ntitle', 'libpotph.json')
         call json%get('ihole',  ihole, found)
                if (.not. found) call bailout('ihole',  'libpotph.json')
         call json%get('ipr1',   ipr1, found)
                if (.not. found) call bailout('ipr1',   'libpotph.json')
         call json%get('iafolp', iafolp, found)
                if (.not. found) call bailout('iafolp', 'libpotph.json')
         call json%get('ixc',    ixc, found)
                if (.not. found) call bailout('ixc',    'libpotph.json')
         call json%get('ispec',  ispec, found)
                if (.not. found) call bailout('ispec',  'libpotph.json')

         call json%get('nmix',   nmix, found)
                if (.not. found) call bailout('nmix',   'libpotph.json')
         call json%get('nohole', nohole, found)
                if (.not. found) call bailout('nohole', 'libpotph.json')
         call json%get('jumprm', jumprm, found)
                if (.not. found) call bailout('jumprm', 'libpotph.json')
         call json%get('inters', inters, found)
                if (.not. found) call bailout('inters', 'libpotph.json')
         call json%get('nscmt',  nscmt, found)
                if (.not. found) call bailout('nscmt',  'libpotph.json')
         call json%get('icoul',  icoul, found)
                if (.not. found) call bailout('icoul',  'libpotph.json')
         call json%get('lfms1',  lfms1, found)
                if (.not. found) call bailout('lfms1',  'libpotph.json')
         call json%get('iunf',   iunf, found)
                if (.not. found) call bailout('iunf',   'libpotph.json')

         call json%get('gamach', gamach, found)
                if (.not. found) call bailout('gamach', 'libpotph.json')
         call json%get('rgrd',   rgrd, found)
                if (.not. found) call bailout('rgrd',   'libpotph.json')
         call json%get('ca1',    ca1, found)
                if (.not. found) call bailout('ca1',    'libpotph.json')
         call json%get('ecv',    ecv, found)
                if (.not. found) call bailout('ecv',    'libpotph.json')
         call json%get('totvol', totvol, found)
                if (.not. found) call bailout('totvol', 'libpotph.json')
         call json%get('rfms1',  toss, found)   
                if (.not. found) call bailout('rfms1',  'libpotph.json')
         rfms1 = real(toss)

         call json%get('titles', strings, found)
                if (.not. found) call bailout('titles', 'libpotph.json')
         do 1000 itit = 1, nheadx
            title(itit) = strings(itit)            
 1000    continue
         call json%get('iz', intgs, found)
                if (.not. found) call bailout('iz',     'libpotph.json')
         do 1010 iph = 0, nphx
            iz(iph) = intgs(iph+1)            
 1010    continue
         call json%get('lmaxsc', intgs, found)
                if (.not. found) call bailout('lmaxsc', 'libpotph.json')
         do 1020 iph = 0, nphx
            lmaxsc(iph) = intgs(iph+1)            
 1020    continue
         call json%get('xnatph', dbpcs, found)
                if (.not. found) call bailout('xnatph', 'libpotph.json')
         do 1030 iph = 0, nphx
            xnatph(iph) = dbpcs(iph+1)            
 1030    continue
         call json%get('xion', dbpcs, found)
                if (.not. found) call bailout('xion',   'libpotph.json')
         do 1040 iph = 0, nphx
            xion(iph) = dbpcs(iph+1)            
 1040    continue
         call json%get('folp', dbpcs, found)
                if (.not. found) call bailout('folp',   'libpotph.json')
         do 1050 iph = 0, nphx
            folp(iph) = dbpcs(iph+1)            
 1050    continue



c        content of xsph/mod2.json
         call json%get('mphase',   mphase, found)
                if (.not. found) call bailout('mphase', 'libpotph.json')
         call json%get('ipr2',   ipr2, found)
                if (.not. found) call bailout('ipr2',   'libpotph.json')
         call json%get('ixc',   ixc, found)
                if (.not. found) call bailout('ixc',    'libpotph.json')
         call json%get('ixc0',   ixc0, found)
                if (.not. found) call bailout('ixc0',   'libpotph.json')
         call json%get('lreal',   lreal, found)
                if (.not. found) call bailout('lreal',  'libpotph.json')
         call json%get('lfms2',   lfms2, found)
                if (.not. found) call bailout('lfms2',  'libpotph.json')
         call json%get('l2lp',   l2lp, found)
                if (.not. found) call bailout('l2lp',   'libpotph.json')
         call json%get('iPlsmn', iPl, found)
                if (.not. found) call bailout('iPlsmn', 'libpotph.json')
         call json%get('iGrid',   iGrid, found)
                if (.not. found) call bailout('iGrid',  'libpotph.json')

         call json%get('vro',   vr0, found)
                if (.not. found) call bailout('vr0',    'libpotph.json')
         call json%get('vio',   vi0, found)
                if (.not. found) call bailout('vi0',    'libpotph.json')

         call json%get('rfms2',  toss, found)
                if (.not. found) call bailout('rfms2',  'libpotph.json')
         rfms2 = real(toss)
         call json%get('gamach',   gamach, found)
                if (.not. found) call bailout('gamach', 'libpotph.json')
         call json%get('xkstep',   xkstep, found)
                if (.not. found) call bailout('xkstep', 'libpotph.json')
         call json%get('xkmax',   xkmax, found)
                if (.not. found) call bailout('xkmax',  'libpotph.json')
         call json%get('vixan',   vixan, found)
                if (.not. found) call bailout('vixan',  'libpotph.json')

         call json%get('izstd',   izstd, found)
                if (.not. found) call bailout('izstd')
         call json%get('ifxc',   ifxc, found)
                if (.not. found) call bailout('ifxc')
         call json%get('ipmbse',   ipmbse, found)
                if (.not. found) call bailout('ipmbse')
         call json%get('itdlda',   itdlda, found)
                if (.not. found) call bailout('itdlda')
         call json%get('nonlocal',   nonlocal, found)
                if (.not. found) call bailout('nonlocal')
         call json%get('ibasis',   ibasis, found)
                if (.not. found) call bailout('ibasis')

         call json%get('potlbl', strings, found)
                if (.not. found) call bailout('potlbl', 'libpotph.json')
         do 2000 itit = 1, nphx
c            potlbl(itit-1) = strings(itit)
            potlbl(itit-1) = strings(itit)(1:6)
 2000    continue
         call json%get('lmaxph', intgs, found)
                if (.not. found) call bailout('lmaxph', 'libpotph.json')
         do 2010 iph = 0, nphx
            lmaxph(iph) = intgs(iph+1)            
 2010    continue
         call json%get('spinph', dbpcs, found)
                if (.not. found) call bailout('spinph', 'libpotph.json')
         do 2020 iph = 0, nphx
            spinph(iph) = dbpcs(iph+1)            
 2020    continue


         call json%destroy()
      end if

      return
      end
