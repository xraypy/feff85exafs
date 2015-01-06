      subroutine reff2x(mchi, ispec, ipr6, idwopt, critcw, s02, sig2g,
     1                   tk, thetad, mbconv, absolu,  !KJ added absolu 3-06
     1                   vrcorr, vicorr, alphat, thetae, nabs,
     4            elnes,ipmin,ipmax,ipstep)   !KJ added this line  1-06     

      use json_module
      implicit double precision (a-h, o-z)
      logical :: found
      type(json_file) :: json   !the JSON structure read from the file:

      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'

cc    global.dat 
c       configuration average
        integer nabs, iphabs
c       global polarization data
        integer ipol, ispin, le2
        double precision rclabs
        double precision evec(3), xivec(3), spvec(3), angks, elpty
        complex*16 ptz(-1:1, -1:1)
cc    mod6.inp
        integer  mchi, idwopt, ipr6, mbconv, absolu  !KJ added absolu 3-06
        double precision  vrcorr, vicorr, s02, tk, thetad
        double precision  alphat, thetae, sig2g

        integer elnes,ipmin,ipmax,ipstep  !KJ my variables  1-06
c     Local stuff
c      character*512 slog

c     standard formats for string, integers and real numbers
c  10  format(a)
c  20  format (20i4)
c  30  format (6f13.5)

c--json--cc    read global.inp
c--json--      open (file='global.dat', unit=3, status='unknown',iostat=ios)
c--json--c       configuration average data
c--json--        read  (3, 10) slog
c--json--        read  (3, 45) nabs, iphabs, rclabs
c--json--  45    format ( 2i8, f13.5)
c--json--      close(3)
c--json--c     read mod6.inp
c--json--      open (file='mod6.inp', unit=3, status='old',iostat=ios)
c--json--        read (3,10)  slog
c--json--        read (3,20)  mchi, ispec, idwopt, ipr6, mbconv, absolu !KJ added absolu 3-06
c--json--        read (3,10)  slog
c--json--        read (3,30)  vrcorr, vicorr, s02, critcw
c--json--        read (3,10)  slog
c--json--        read (3,30)  tk, thetad, alphat, thetae, sig2g
c--json--      close(3)

      call json_read_global(nabs, iphabs, rclabs, ipol, ispin, le2,
     1                      elpty, angks, evec, xivec, spvec, ptz)


      call json%load_file('ff2x.json')
      if (json_failed()) then   !if there was an error reading the file
         print *, "failed to read ff2x.json"
         stop
      else
         call json%get('mchi',   mchi, found)
                   if (.not. found) call bailout('mchi', 'ff2x.json')
         call json%get('ispec',   ispec, found)
                   if (.not. found) call bailout('ispec', 'ff2x.json')
         call json%get('idwopt',   idwopt, found)
                   if (.not. found) call bailout('idwopt', 'ff2x.json')
         call json%get('ipr6',   ipr6, found)
                   if (.not. found) call bailout('ipr6', 'ff2x.json')
         call json%get('mbconv',   mbconv, found)
                   if (.not. found) call bailout('mbconv', 'ff2x.json')
         call json%get('absolu',   absolu, found)
                   if (.not. found) call bailout('absolu', 'ff2x.json')
         call json%get('vrcorr',   vrcorr, found)
                   if (.not. found) call bailout('vrcorr', 'ff2x.json')
         call json%get('vicorr',   vicorr, found)
                   if (.not. found) call bailout('vicorr', 'ff2x.json')
         call json%get('s02',   s02, found)
                   if (.not. found) call bailout('s02', 'ff2x.json')
         call json%get('critcw',   critcw, found)
                   if (.not. found) call bailout('critcw', 'ff2x.json')
         call json%get('tk',   tk, found)
                   if (.not. found) call bailout('tk', 'ff2x.json')
         call json%get('thetad',   thetad, found)
                   if (.not. found) call bailout('thetad', 'ff2x.json')
         call json%get('alphat',   alphat, found)
                   if (.not. found) call bailout('alphat', 'ff2x.json')
         call json%get('thetae',   thetae, found)
                   if (.not. found) call bailout('thetae', 'ff2x.json')
         call json%get('sig2g',   sig2g, found)
                   if (.not. found) call bailout('sig2g', 'ff2x.json')
         call json%destroy()
      end if
      
      elnes  = 0
      ipstep = 1
      ipmax  = 1
      ipmin  = 1

c     transform energies to atomic units
      vrcorr = vrcorr / hart
      vicorr = vicorr / hart

      return
      end
