      subroutine rexsph ( mphase, ipr2, ispec, vixan, xkstep, xkmax,
     1             gamach, rgrd,
     1             nph, lmaxph, potlbl, spinph, iatph, nat, rat, iphat,
     2             ixc, vr0, vi0, ixc0, lreal, rfms2, lfms2, l2lp,
     3             ipol, ispin, le2, angks, ptz, iPl, iGrid,
     4             izstd, ifxc, ipmbse, itdlda, nonlocal, ibasis)

      use json_module
      implicit double precision (a-h, o-z)
      logical :: found
      type(json_file) :: json   !the JSON structure read from the file:
      double precision toss
      integer,dimension(:),allocatable :: intgs
      character*80,dimension(:),allocatable :: strings
      double precision,dimension(:),allocatable :: dbpcs

      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'

cc    geom.dat
      integer  nat, iatph(0:nphx), iphat(natx), ibounc(natx)
      double precision  rat(3,natx)
cc    global.dat
c       configuration average
        integer nabs, iphabs
c       global polarization data
        integer  ipol, ispin, le2
        double precision evec(3), xivec(3), spvec(3), elpty,angks,rclabs
        complex*16 ptz(-1:1, -1:1)
cc    mod2.inp
        integer mphase, ipr2, ixc, ixc0, ispec, lreal, lfms2, l2lp, iPl, 
     &       iGrid
        double precision rgrd, gamach, xkstep, xkmax, vixan
        double precision vr0, vi0, spinph(0:nphx)
        real rfms2
        integer lmaxph(0:nphx)
        character*6  potlbl(0:nphx)
        integer izstd, ifxc, ipmbse, itdlda, nonlocal, ibasis

c     Local stuff
c      character*512 slog
c      character*80 head(nheadx)
c      dimension lhead(nheadx)

c     standard formats for string, integers and real numbers
c  10  format(a)
c  20  format (20i4)
c  30  format (6f13.5)


      call json_read_geom(nat, nph, iatph, rat, iphat, ibounc)

      call json_read_global(nabs, iphabs, rclabs, ipol, ispin, le2,
     1                      elpty, angks, evec, xivec, spvec, ptz)


      call json%load_file('xsph.json')
      if (json_failed()) then   !if there was an error reading the file
         print *, "failed to read xsph.json"
         stop
      else
         call json%get('mphase',   mphase, found)
                   if (.not. found) call bailout('mphase', 'xsph.json')
         call json%get('ipr2',   ipr2, found)
                   if (.not. found) call bailout('ipr2', 'xsph.json')
         call json%get('ixc',   ixc, found)
                   if (.not. found) call bailout('ixc', 'xsph.json')
         call json%get('ixc0',   ixc0, found)
                   if (.not. found) call bailout('ixc0', 'xsph.json')
         call json%get('ispec',   ispec, found)
                   if (.not. found) call bailout('ispec', 'xsph.json')
         call json%get('lreal',   lreal, found)
                   if (.not. found) call bailout('lreal', 'xsph.json')
         call json%get('lfms2',   lfms2, found)
                   if (.not. found) call bailout('lfms2', 'xsph.json')
         call json%get('nph',   nph, found)
                   if (.not. found) call bailout('nph', 'xsph.json')
         call json%get('l2lp',   l2lp, found)
                   if (.not. found) call bailout('l2lp', 'xsph.json')
         call json%get('iPlsmn', iPl, found)
                   if (.not. found) call bailout('iPlsmn', 'xsph.json')
         call json%get('iGrid',   iGrid, found)
                   if (.not. found) call bailout('iGrid', 'xsph.json')

         call json%get('vro',   vr0, found)
                   if (.not. found) call bailout('vr0', 'xsph.json')
         call json%get('vio',   vi0, found)
                   if (.not. found) call bailout('vi0', 'xsph.json')

         call json%get('rgrd',   rgrd, found)
                   if (.not. found) call bailout('rgrd', 'xsph.json')
         call json%get('rfms2',  toss, found)
                   if (.not. found) call bailout('rfms2', 'xsph.json')
         rfms2 = real(toss)
         call json%get('gamach',   gamach, found)
                   if (.not. found) call bailout('gamach', 'xsph.json')
         call json%get('xkstep',   xkstep, found)
                   if (.not. found) call bailout('xkstep', 'xsph.json')
         call json%get('xkmax',   xkmax, found)
                   if (.not. found) call bailout('xkmax', 'xsph.json')
         call json%get('vixan',   vixan, found)
                   if (.not. found) call bailout('vixan', 'xsph.json')

         call json%get('izstd',   izstd, found)
                   if (.not. found) call bailout('izstd', 'xsph.json')
         call json%get('ifxc',   ifxc, found)
                   if (.not. found) call bailout('ifxc', 'xsph.json')
         call json%get('ipmbse',   ipmbse, found)
                   if (.not. found) call bailout('ipmbse', 'xsph.json')
         call json%get('itdlda',   itdlda, found)
                   if (.not. found) call bailout('itdlda', 'xsph.json')
         call json%get('nonlocal',   nonlocal, found)
                if (.not. found) call bailout('nonlocal', 'xsph.json')
         call json%get('ibasis',   ibasis, found)
                   if (.not. found) call bailout('ibasis', 'xsph.json')

         call json%get('potlbl', strings, found)
                   if (.not. found) call bailout('potlbl', 'xsph.json')
         do itit = 1, nphx
c            potlbl(itit-1) = strings(itit)
            potlbl(itit-1) = strings(itit)(1:6)
         enddo
         call json%get('lmaxph', intgs, found)
                   if (.not. found) call bailout('lmaxph', 'xsph.json')
         do iph = 0, nphx
            lmaxph(iph) = intgs(iph+1)
         enddo
         call json%get('spinph', dbpcs, found)
                   if (.not. found) call bailout('spinph', 'xsph.json')
         do iph = 0, nphx
            spinph(iph) = dbpcs(iph+1)
         enddo

         call json%destroy()
      end if


c     transform to code units (bohrs and hartrees - atomic unuts)
      rfms2 = rfms2 / real(bohr)
      vr0   = vr0 / hart
      vi0   = vi0 / hart
      gamach = gamach / hart
      vixan = vixan / hart
      xkstep = xkstep * bohr
      xkmax  = xkmax  * bohr
      do i = 1,3
         do iat = 1, nat
            rat(i,iat) = rat(i,iat) / bohr
         enddo
      enddo

      return
      end
