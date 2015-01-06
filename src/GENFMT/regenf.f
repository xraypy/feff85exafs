      subroutine regenf (mfeff, ipr5, critcw, iorder, wnstar,
     1            ipol, ispin, le2, angks, elpty, evec, xivec, ptz)

      use json_module
      implicit double precision (a-h, o-z)
      logical :: found
      type(json_file) :: json   !the JSON structure read from the file:

      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'
      
      integer elnes,ipmin,ipmax,ipstep  !KJ added these variables 1-06      

cc    global.dat 
c       configuration average
        integer nabs, iphabs
c       global polarization data
        integer  ipol, ispin, le2
        double precision evec(3), xivec(3), spvec(3), elpty,angks,rclabs
        complex*16 ptz(-1:1, -1:1)
cc    mod5.inp
        integer  mfeff, ipr5, iorder
        logical  wnstar
        double precision critcw

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
c--json--c       global polarization data
c--json--        read  (3,10)   slog
c--json--        read  (3, 50)  ipol, ispin, le2, elpty, angks
c--json--  50    format ( 3i5, 2f12.4)
c--json--        read  (3, 10) slog
c--json--        do 60 i = 1,3
c--json--          read  (3,30) evec(i), xivec(i), spvec(i)
c--json--  60    continue
c--json--        read  (3, 10) slog
c--json--        do 70 i = -1, 1
c--json--          read (3,30) a1, b1, a2, b2, a3, b3
c--json--          ptz(-1,i)= dcmplx(a1,b1)
c--json--          ptz(0,i) = dcmplx(a2,b2)
c--json--          ptz(1,i) = dcmplx(a3,b3)
c--json--  70    continue
c--json--      close(3)

      call json_read_global(nabs, iphabs, rclabs, ipol, ispin, le2,
     1                      elpty, angks, evec, xivec, spvec, ptz)


c--json--c     read mod5.inp
c--json--      open (file='mod5.inp', unit=3, status='old',iostat=ios)
c--json--        read (3,10)  slog
c--json--        read (3,180)  mfeff, ipr5, iorder, critcw, wnstar
c--json--  180   format ( 2i4, i8, f13.5, L5)
c--json--      close(3)


      call json%load_file('genfmt.json')
      if (json_failed()) then   !if there was an error reading the file
         print *, "failed to read genfmt.json"
         stop
      else
         call json%get('mfeff',   mfeff, found)
                  if (.not. found) call bailout('mfeff', 'genfmt.json')
         call json%get('ipr5',   ipr5, found)
                  if (.not. found) call bailout('ipr5', 'genfmt.json')
         call json%get('iorder',   iorder, found)
                  if (.not. found) call bailout('iorder', 'genfmt.json')
         call json%get('critcw',   critcw, found)
                  if (.not. found) call bailout('critcw', 'genfmt.json')
         call json%get('wnstar',   wnstar, found)
                  if (.not. found) call bailout('wnstar', 'genfmt.json')
         call json%destroy()
      end if

      elnes  = 0
      ipstep = 1
      ipmax  = 1
      ipmin  = 1

      return
      end
