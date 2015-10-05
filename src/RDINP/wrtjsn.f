      subroutine wrtjsn (nabs)
c     writes data stored in common blocks of allinp.h to 
c     all necessary input files for other modules.
c     version 1.0 written by Alexei Ankudinov, March 2001

c     Note: to add input variable one has to add it to the 
c        appropriate common block in allinp.h, properly initialize
c        it in subroutine iniall and modify subroutine wrtall
c        to write it to the appropriate input file.
c        (i.e. one has to make modifications in 3 places)

      use json_module

      integer ntoss
      ntoss = nabs

!     initialize the module:
c      call json_initialize()

      call json_initialize()
      call json_mod1()
      call json_mod2()
c      call json_mod3()
      call json_mod4()
      call json_mod5()
      call json_mod6()

c      call json_s02()
      call json_atoms()
      call json_global(nabs)

      call json_libpotph()

      return
      end


c$$$ template for a json export

c$$$      subroutine json_modN
c$$$
c$$$      use json_module
c$$$
c$$$      implicit double precision (a-h, o-z)
c$$$      include '../HEADERS/dim.h'
c$$$      include '../HEADERS/parallel.h'
c$$$      include '../RDINP/allinp.h'
c$$$
c$$$      integer  iunit
c$$$      type(json_value),pointer :: mN
c$$$
c$$$
c$$$      ! root
c$$$      call json_value_create(mN)      ! create the value and associate the pointer
c$$$      call to_object(mN,'modN.json')  ! add the file name as the name of the overall structure
c$$$
c$$$      
c$$$      call json_value_add(mN, 'somevariable', somevariable)
c$$$      
c$$$      open(newunit=iunit, file='modN.json', status='REPLACE')
c$$$      call json_print(mN,iunit)
c$$$      close(iunit)
c$$$
c$$$!     cleanup:
c$$$      call json_destroy(mN)
c$$$
c$$$
c$$$      return
c$$$      end




      subroutine json_mod1

      use json_module

      implicit double precision (a-h, o-z)
      include '../HEADERS/vers.h'
      include '../HEADERS/dim.h'
      include '../HEADERS/parallel.h'
      include '../RDINP/allinp.h'

      integer  iunit
      type(json_value),pointer :: m1
      character*7 vname
      integer itoss(novrx)
      double precision rtoss(novrx)

c     root
      call json_value_create(m1)      ! create the value and associate the pointer
      call to_object(m1,'pot.json')  ! add the file name as the name of the overall structure
      
      call json_value_add(m1, 'mod',    1)
      call json_value_add(m1, 'vfeff',  vfeff)
      call json_value_add(m1, 'vf85e',  vf85e)

      call json_value_add(m1, 'mpot',   mpot)
      call json_value_add(m1, 'nph',    nph)
      call json_value_add(m1, 'ntitle', ntitle)
      call json_value_add(m1, 'ihole',  ihole)
      call json_value_add(m1, 'ipr1',   ipr1)
      call json_value_add(m1, 'iafolp', iafolp)
      call json_value_add(m1, 'ixc',    ixc)
      call json_value_add(m1, 'ispec',  ispec)

      call json_value_add(m1, 'nmix',   nmix)
      call json_value_add(m1, 'nohole', nohole)
      call json_value_add(m1, 'jumprm', jumprm)
      call json_value_add(m1, 'inters', inters)
      call json_value_add(m1, 'nscmt',  nscmt)
      call json_value_add(m1, 'icoul',  icoul)
      call json_value_add(m1, 'lfms1',  lfms1)
      call json_value_add(m1, 'iunf',   iunf)

      call json_value_add(m1, 'gamach', gamach)
      call json_value_add(m1, 'rgrd',   rgrd)
      call json_value_add(m1, 'ca1',    ca1)
      call json_value_add(m1, 'ecv',    ecv)
      call json_value_add(m1, 'totvol', totvol)
      call json_value_add(m1, 'rfms1',  dble(rfms1))

      call json_value_add(m1, 'titles', title)

      call json_value_add(m1, 'iz',     iz)
      call json_value_add(m1, 'lmaxsc', lmaxsc)
      call json_value_add(m1, 'xnatph', xnatph)
      call json_value_add(m1, 'xion',   xion)
      call json_value_add(m1, 'folp',   folp)
      call json_value_add(m1, 'novr',   novr)

c     the following disentangles the 2D arrays used for overlap into
c     arrays that can be stored as json arrays, see POT/reapot.f line
c     188 and following for how this gets reconstructed into a 2D array
      do 10 iph = 0, nph
         write (vname, "(A6,I1)") "iphovr", iph
         do 20 iovr = 1, novrx
            itoss(iovr) = iphovr(iovr, iph)
 20      continue
         call json_value_add(m1, vname, itoss)

         write (vname, "(A5,I1)") "nnovr", iph
         do 30 iovr = 1, novrx
            itoss(iovr) = nnovr(iovr, iph)
 30      continue
         call json_value_add(m1, vname, itoss)

         write (vname, "(A4,I1)") "rovr", iph
         do 40 iovr = 1, novrx
            rtoss(iovr) = rovr(iovr, iph)
 40      continue
         call json_value_add(m1, vname, rtoss)
 10   continue

      open(newunit=iunit, file='pot.json', status='REPLACE')
      call json_print(m1,iunit)
      close(iunit)

!     cleanup:
      call json_destroy(m1)


      return
      end




      subroutine json_mod2

      use json_module

      implicit double precision (a-h, o-z)
      include '../HEADERS/vers.h'
      include '../HEADERS/dim.h'
      include '../HEADERS/parallel.h'
      include '../RDINP/allinp.h'

      integer  iunit
      type(json_value),pointer :: m2


      ! root
      call json_value_create(m2)      ! create the value and associate the pointer
      call to_object(m2,'xsph.json')  ! add the file name as the name of the overall structure

      call json_value_add(m2, 'mod',      2)
      call json_value_add(m2, 'vfeff',    vfeff)
      call json_value_add(m2, 'vf85e',    vf85e)

      call json_value_add(m2, 'mphase',   mphase)
      call json_value_add(m2, 'ipr2',     ipr2)
      call json_value_add(m2, 'ixc',      ixc)
      call json_value_add(m2, 'ixc0',     ixc0)
      call json_value_add(m2, 'ispec',    ispec)
      call json_value_add(m2, 'lreal',    lreal)
      call json_value_add(m2, 'lfms2',    lfms2)
      call json_value_add(m2, 'nph',      nph)
      call json_value_add(m2, 'l2lp',     l2lp)
      call json_value_add(m2, 'iPlsmn',   iPlsmn)
      call json_value_add(m2, 'iGrid',    iGrid)

      call json_value_add(m2, 'vro',      vr0)
      call json_value_add(m2, 'vio',      vi0)

      call json_value_add(m2, 'rgrd',     rgrd)
      call json_value_add(m2, 'rfms2',    rfms2)
      call json_value_add(m2, 'gamach',   gamach)
      call json_value_add(m2, 'xkstep',   xkstep)
      call json_value_add(m2, 'xkmax',    xkmax)
      call json_value_add(m2, 'vixan',    vixan)
      call json_value_add(m2, 'izstd',    izstd)
      call json_value_add(m2, 'ifxc',     ifxc)
      call json_value_add(m2, 'ipmbse',   ipmbse)
      call json_value_add(m2, 'itdlda',   itdlda)
      call json_value_add(m2, 'nonlocal', nonlocal)
      call json_value_add(m2, 'ibasis',   ibasis)

      call json_value_add(m2, 'lmaxph', lmaxph)
      call json_value_add(m2, 'potlbl', potlbl)
      call json_value_add(m2, 'spinph', spinph)



      open(newunit=iunit, file='xsph.json', status='REPLACE')
      call json_print(m2,iunit)
      close(iunit)

!     cleanup:
      call json_destroy(m2)


      return
      end



      subroutine json_mod3

      use json_module

      implicit double precision (a-h, o-z)
      include '../HEADERS/vers.h'
      include '../HEADERS/dim.h'
      include '../HEADERS/parallel.h'
      include '../RDINP/allinp.h'

      integer  iunit
      type(json_value),pointer :: m3


      ! root
      call json_value_create(m3)      ! create the value and associate the pointer
      call to_object(m3,'fms.json')  ! add the file name as the name of the overall structure

      call json_value_add(m3, 'mod',    3)
      call json_value_add(m3, 'vfeff',  vfeff)
      call json_value_add(m3, 'vf85e',  vf85e)

      call json_value_add(m3, 'mfms',   mfms)
      call json_value_add(m3, 'idwopt', idwopt)
      call json_value_add(m3, 'minv',   minv)
      call json_value_add(m3, 'rfms2',  rfms2)
      call json_value_add(m3, 'rdirec', dble(rdirec))
      call json_value_add(m3, 'toler1', dble(toler1))
      call json_value_add(m3, 'toler2', dble(toler2))
      call json_value_add(m3, 'tk',     tk)
      call json_value_add(m3, 'thetad', thetad)
      call json_value_add(m3, 'sig2g',  sig2g)


      call json_value_add(m3, 'lmaxph', lmaxph)

      
      open(newunit=iunit, file='fms.json', status='REPLACE')
      call json_print(m3,iunit)
      close(iunit)

!     cleanup:
      call json_destroy(m3)


      return
      end




      subroutine json_mod4

      use json_module

      implicit double precision (a-h, o-z)
      include '../HEADERS/vers.h'
      include '../HEADERS/dim.h'
      include '../HEADERS/parallel.h'
      include '../RDINP/allinp.h'

      integer  iunit
      type(json_value),pointer :: m4


      ! root
      call json_value_create(m4)      ! create the value and associate the pointer
      call to_object(m4,'path.json')  ! add the file name as the name of the overall structure

      call json_value_add(m4, 'mod',    4)
      call json_value_add(m4, 'vfeff',  vfeff)
      call json_value_add(m4, 'vf85e',  vf85e)

      call json_value_add(m4, 'mpath',  mpath)
      call json_value_add(m4, 'ms',     ms)
      call json_value_add(m4, 'nncrit', nncrit)
      call json_value_add(m4, 'nlegxx', nlegxx)
      call json_value_add(m4, 'ipr4',   ipr4)
      call json_value_add(m4, 'critpw', dble(critpw))
      call json_value_add(m4, 'pcritk', dble(pcritk))
      call json_value_add(m4, 'pcrith', dble(pcrith))
      call json_value_add(m4, 'rmax',   dble(rmax))
      call json_value_add(m4, 'rfms2',  dble(rfms2))
      
      open(newunit=iunit, file='path.json', status='REPLACE')
      call json_print(m4,iunit)
      close(iunit)

!     cleanup:
      call json_destroy(m4)


      return
      end



      subroutine json_mod5

      use json_module

      implicit double precision (a-h, o-z)
      include '../HEADERS/vers.h'
      include '../HEADERS/dim.h'
      include '../HEADERS/parallel.h'
      include '../RDINP/allinp.h'

      integer  iunit
      type(json_value),pointer :: m5


      ! root
      call json_value_create(m5)      ! create the value and associate the pointer
      call to_object(m5,'genfmt.json')  ! add the file name as the name of the overall structure

      call json_value_add(m5, 'mod',    5)
      call json_value_add(m5, 'vfeff',  vfeff)
      call json_value_add(m5, 'vf85e',  vf85e)
      
      call json_value_add(m5, 'mfeff',  mfeff)
      call json_value_add(m5, 'ipr5',   ipr5)
      call json_value_add(m5, 'iorder', iorder)
      call json_value_add(m5, 'critcw', dble(critcw))
      call json_value_add(m5, 'wnstar', wnstar)
     
      open(newunit=iunit, file='genfmt.json', status='REPLACE')
      call json_print(m5,iunit)
      close(iunit)

!     cleanup:
      call json_destroy(m5)


      return
      end





      subroutine json_mod6

      use json_module

      implicit double precision (a-h, o-z)
      include '../HEADERS/vers.h'
      include '../HEADERS/dim.h'
      include '../HEADERS/parallel.h'
      include '../RDINP/allinp.h'

      integer  iunit
      type(json_value),pointer :: m6


      ! root
      call json_value_create(m6)      ! create the value and associate the pointer
      call to_object(m6,'ff2x.json')  ! add the file name as the name of the overall structure

      call json_value_add(m6, 'mod',    6)
      call json_value_add(m6, 'vfeff',  vfeff)
      call json_value_add(m6, 'vf85e',  vf85e)
      
      call json_value_add(m6, 'mchi',   mchi)
      call json_value_add(m6, 'ispec',  ispec)
      call json_value_add(m6, 'idwopt', idwopt)
      call json_value_add(m6, 'ipr6',   ipr6)
      call json_value_add(m6, 'mbconv', mbconv)
      call json_value_add(m6, 'absolu', absolu)
      call json_value_add(m6, 'vrcorr', vrcorr)
      call json_value_add(m6, 'vicorr', vicorr)
      call json_value_add(m6, 's02',    s02)
      call json_value_add(m6, 'critcw', critcw)
      call json_value_add(m6, 'tk',     tk)
      call json_value_add(m6, 'thetad', thetad)
      call json_value_add(m6, 'alphat', alphat)
      call json_value_add(m6, 'thetae', thetae)
      call json_value_add(m6, 'sig2g',  sig2g)

      open(newunit=iunit, file='ff2x.json', status='REPLACE')
      call json_print(m6,iunit)
      close(iunit)

!     cleanup:
      call json_destroy(m6)


      return
      end






      subroutine json_s02

      use json_module

      implicit double precision (a-h, o-z)
      include '../HEADERS/vers.h'
      include '../HEADERS/dim.h'
      include '../HEADERS/parallel.h'
      include '../RDINP/allinp.h'

      integer  iunit
      type(json_value),pointer :: so2


      ! root
      call json_value_create(so2)      ! create the value and associate the pointer
      call to_object(so2,'s02.json')  ! add the file name as the name of the overall structure

      call json_value_add(so2, 'vfeff',  vfeff)
      call json_value_add(so2, 'vf85e',  vf85e)

      call json_value_add(so2, 'mso2conv', mso2conv)
      call json_value_add(so2, 'ipse',     ipse)
      call json_value_add(so2, 'ipsk',     ipsk)
      call json_value_add(so2, 'wsigk',    wsigk)
      call json_value_add(so2, 'cen',      cen)
      call json_value_add(so2, 'ispec',    ispec)
      call json_value_add(so2, 'ipr6',     ipr6)
      call json_value_add(so2, 'cfname',   cfname)

      
      open(newunit=iunit, file='s02.json', status='REPLACE')
      call json_print(so2,iunit)
      close(iunit)

!     cleanup:
      call json_destroy(so2)


      return
      end



      subroutine json_atoms

      use json_module

      implicit double precision (a-h, o-z)
      include '../HEADERS/vers.h'
      include '../HEADERS/dim.h'
      include '../HEADERS/parallel.h'
      include '../RDINP/allinp.h'
      double precision xx(nattx), yy(nattx), zz(nattx)

      integer  iunit
      type(json_value),pointer :: atoms


      ! root
      call json_value_create(atoms)      ! create the value and associate the pointer
      call to_object(atoms,'atoms.json')  ! add the file name as the name of the overall structure

      call json_value_add(atoms, 'vfeff',  vfeff)
      call json_value_add(atoms, 'vf85e',  vf85e)
      
      call json_value_add(atoms, 'natt', natt)
      do 10 i=1,natt
         xx(i) = ratx(1,i)
         yy(i) = ratx(2,i)
         zz(i) = ratx(3,i)
 10   continue
      call json_value_add(atoms, 'x', xx(1:natt))
      call json_value_add(atoms, 'y', yy(1:natt))
      call json_value_add(atoms, 'z', zz(1:natt))
      call json_value_add(atoms, 'iphatx', iphatx(1:natt))

      
      open(newunit=iunit, file='atoms.json', status='REPLACE')
      call json_print(atoms,iunit)
      close(iunit)

!     cleanup:
      call json_destroy(atoms)


      return
      end




      subroutine json_global(nabs)

      use json_module

      implicit double precision (a-h, o-z)
      include '../HEADERS/vers.h'
      include '../HEADERS/dim.h'
      include '../HEADERS/parallel.h'
      include '../RDINP/allinp.h'

      integer  iunit, nabs
      type(json_value),pointer :: global
      character*5 vname

      ! root
      call json_value_create(global)       ! create the value and associate the pointer
      call to_object(global,'global.json') ! add the file name as the name of the overall structure

      call json_value_add(global, 'vfeff',  vfeff)
      call json_value_add(global, 'vf85e',  vf85e)

      call json_value_add(global, 'nabs',   nabs)
      call json_value_add(global, 'iphabs', iphabs)
      call json_value_add(global, 'rclabs', rclabs)
      call json_value_add(global, 'ipol',   ipol)
      call json_value_add(global, 'ispin',  ispin)
      call json_value_add(global, 'le2',    le2)
      call json_value_add(global, 'elpty',  elpty)
      call json_value_add(global, 'angks',  angks)

      call json_value_add(global, 'evec',   evec)
      call json_value_add(global, 'xivec',  xivec)
      call json_value_add(global, 'spvec',  spvec)

      do 10 i = -1, 1
         write (vname, "(A3,I1)") "ptz", i+1
         call json_value_add(global, vname,
     1                 [dble(real(ptz(-1,i))), dble(imag(ptz(-1,i))),
     2                  dble(real(ptz( 0,i))), dble(imag(ptz( 0,i))),
     3                  dble(real(ptz( 1,i))), dble(imag(ptz( 1,i)))])
 10   continue

      open(newunit=iunit, file='global.json', status='REPLACE')
      call json_print(global,iunit)
      close(iunit)

!     cleanup:
      call json_destroy(global)


      return
      end




      subroutine json_geom(iatph, rat, iphat)

      use json_module

      implicit double precision (a-h, o-z)
      include '../HEADERS/vers.h'
      include '../HEADERS/dim.h'
      include '../HEADERS/parallel.h'
      include '../RDINP/allinp.h'
      double precision xx(nattx), yy(nattx), zz(nattx)
      double precision rat(3,natx)
      integer ib(nattx), iphat(natx), iatph(0:nphx)

      integer  iunit
      type(json_value),pointer :: geom


      ! root
      call json_value_create(geom)      ! create the value and associate the pointer
      call to_object(geom,'geom.json')  ! add the file name as the name of the overall structure
      
      call json_value_add(geom, 'vfeff',  vfeff)
      call json_value_add(geom, 'vf85e',  vf85e)

      call json_value_add(geom, 'natt',  natt)
      call json_value_add(geom, 'nph',   nph)
      call json_value_add(geom, 'iatph', iatph)
      do 10 i=1,natt
         xx(i) = rat(1,i)
         yy(i) = rat(2,i)
         zz(i) = rat(3,i)
         ib(i) = 1
 10   continue
      call json_value_add(geom, 'x',   xx(1:natt))
      call json_value_add(geom, 'y',   yy(1:natt))
      call json_value_add(geom, 'z',   zz(1:natt))
      call json_value_add(geom, 'iph', iphat(1:natt))
      call json_value_add(geom, 'ibo', ib(1:natt))

      open(newunit=iunit, file='geom.json', status='REPLACE')
      call json_print(geom,iunit)
      close(iunit)

!     cleanup:
      call json_destroy(geom)

      return
      end


      subroutine json_libpotph()

      use json_module

      implicit double precision (a-h, o-z)
      include '../HEADERS/vers.h'
      include '../HEADERS/dim.h'
      include '../HEADERS/parallel.h'
      include '../RDINP/allinp.h'
      double precision xx(nattx), yy(nattx), zz(nattx)
c      double precision rat(3,natx)

      integer  iunit
      type(json_value),pointer :: lpotph
      character*5 vname

      ! root
      call json_value_create(lpotph)        ! create the value and associate the pointer
      call to_object(lpotph,'lpotph.json')  ! add the file name as the name of the overall structure

c     content of atoms.json
      call json_value_add(lpotph, 'vfeff',  vfeff)
      call json_value_add(lpotph, 'vf85e',  vf85e)

c     TITLE
      call json_value_add(lpotph, 'ntitle',   ntitle)
      call json_value_add(lpotph, 'titles', title)

c     ATOMS
      call json_value_add(lpotph, 'natt', natt)
      do 10 i=1,natt
         xx(i) = ratx(1,i)
         yy(i) = ratx(2,i)
         zz(i) = ratx(3,i)
 10   continue
      call json_value_add(lpotph, 'x', xx(1:natt))
      call json_value_add(lpotph, 'y', yy(1:natt))
      call json_value_add(lpotph, 'z', zz(1:natt))
      call json_value_add(lpotph, 'iphatx', iphatx(1:natt))

c     POTENTIALS
      call json_value_add(lpotph, 'nph',    nph)
      call json_value_add(lpotph, 'iz',     iz)
      call json_value_add(lpotph, 'potlbl', potlbl)
      call json_value_add(lpotph, 'lmaxsc', lmaxsc)
      call json_value_add(lpotph, 'lmaxph', lmaxph)
      call json_value_add(lpotph, 'xnatph', xnatph)
      call json_value_add(lpotph, 'spinph', spinph)

c     HOLE/EDGE
      call json_value_add(lpotph, 'ihole',  ihole)

c     SCF
      call json_value_add(lpotph, 'rfms1',  dble(rfms1))
      call json_value_add(lpotph, 'lfms1',  lfms1)
      call json_value_add(lpotph, 'nscmt',  nscmt)
      call json_value_add(lpotph, 'ca1',    ca1)
      call json_value_add(lpotph, 'nmix',   nmix)
      call json_value_add(lpotph, 'ecv',    ecv)
      call json_value_add(lpotph, 'icoul',  icoul)
      
c     POLARIZATION, ELLIPTICITY
      call json_value_add(lpotph, 'ipol',   ipol)
      call json_value_add(lpotph, 'evec',   evec)
      call json_value_add(lpotph, 'elpty',  elpty)
      call json_value_add(lpotph, 'xivec',  xivec)

c     SPIN
      call json_value_add(lpotph, 'ispin',  ispin)
      call json_value_add(lpotph, 'spvec',  spvec)
      call json_value_add(lpotph, 'angks',  angks)

c     computed: ptz and gamach
      do 20 i = -1, 1
         write (vname, "(A3,I1)") "ptz", i+1
         call json_value_add(lpotph, vname,
     1                 [dble(real(ptz(-1,i))), dble(imag(ptz(-1,i))),
     2                  dble(real(ptz( 0,i))), dble(imag(ptz( 0,i))),
     3                  dble(real(ptz( 1,i))), dble(imag(ptz( 1,i)))])
 20   continue
      call json_value_add(lpotph, 'gamach', gamach)

c     EXCHANGE
      call json_value_add(lpotph, 'ixc',    ixc)
      call json_value_add(lpotph, 'vro',    vr0)
      call json_value_add(lpotph, 'vio',    vi0)
      call json_value_add(lpotph, 'ixc0',   ixc0)
      
c     AFOLP, FOLP, ION, RGRID, UNFREEZEF
      call json_value_add(lpotph, 'iafolp', iafolp)
      call json_value_add(lpotph, 'folp',   folp)
      call json_value_add(lpotph, 'xion',   xion)
      call json_value_add(lpotph, 'rgrd',   rgrd)
      call json_value_add(lpotph, 'iunf',   iunf)

c     INTERSTITIAL, JUMPRM, NOHOLE
      call json_value_add(lpotph, 'inters', inters)
      call json_value_add(lpotph, 'totvol', totvol)
      call json_value_add(lpotph, 'jumprm', jumprm)
      call json_value_add(lpotph, 'nohole', nohole)
      call json_value_add(lpotph, 'iplsmn', iPlsmn)


c     CFAVERAGE
c      call json_value_add(lpotph, 'nabs',   nabs)
c      call json_value_add(lpotph, 'iphabs', iphabs)
c      call json_value_add(lpotph, 'rclabs', rclabs)
c     MULTIPOLE
c      call json_value_add(lpotph, 'le2',    le2)

c The rest are parameters associated with:
c    XANES, FMS, PLASMON, PMBSE, TDLDA, RPHASES
c    CONTROL, PRINT 

c      call json_value_add(lpotph, 'iGrid',    iGrid)
c      call json_value_add(lpotph, 'iPlsmn',   iPlsmn)
c      call json_value_add(lpotph, 'ibasis',   ibasis)
c      call json_value_add(lpotph, 'ifxc',     ifxc)
c      call json_value_add(lpotph, 'ipmbse',   ipmbse)
c      call json_value_add(lpotph, 'ipr1',     ipr1)
c      call json_value_add(lpotph, 'ipr2',     ipr2)
c      call json_value_add(lpotph, 'ispec',    ispec)
c      call json_value_add(lpotph, 'itdlda',   itdlda)
c      call json_value_add(lpotph, 'izstd',    izstd)
c      call json_value_add(lpotph, 'l2lp',     l2lp)
c      call json_value_add(lpotph, 'lfms2',    lfms2)
c      call json_value_add(lpotph, 'lreal',    lreal)
c      call json_value_add(lpotph, 'mphase',   mphase)
c      call json_value_add(lpotph, 'mpot',     mpot)
c      call json_value_add(lpotph, 'nonlocal', nonlocal)


c      call json_value_add(lpotph, 'rfms2',  rfms2)
c      call json_value_add(lpotph, 'xkstep', xkstep)
c      call json_value_add(lpotph, 'xkmax',  xkmax)
c      call json_value_add(lpotph, 'vixan',  vixan)





      open(newunit=iunit, file='libpotph.json', status='REPLACE')
      call json_print(lpotph,iunit)
      close(iunit)

!     cleanup:
      call json_destroy(lpotph)

      return
      end
