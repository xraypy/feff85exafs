      subroutine regenf (mfeff, ipr5, critcw, iorder, wnstar,
     1            ipol, ispin, le2, angks, elpty, evec, xivec, ptz,
     2            elnes,ipmin,ipmax,ipstep)   !KJ added this line   1-06     

      implicit double precision (a-h, o-z)

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
c       global polarization data
        read  (3,10)   slog
        read  (3, 50)  ipol, ispin, le2, elpty, angks
  50    format ( 3i5, 2f12.4)
        read  (3, 10) slog
        do 60 i = 1,3
          read  (3,30) evec(i), xivec(i), spvec(i)
  60    continue
        read  (3, 10) slog
        do 70 i = -1, 1
          read (3,30) a1, b1, a2, b2, a3, b3
          ptz(-1,i)= cmplx(a1,b1)
          ptz(0,i) = cmplx(a2,b2)
          ptz(1,i) = cmplx(a3,b3)
  70    continue
      close(3)
c     read mod5.inp
      open (file='mod5.inp', unit=3, status='old',iostat=ios)
        read (3,10)  slog
        read (3,180)  mfeff, ipr5, iorder, critcw, wnstar
  180   format ( 2i4, i8, f13.5, L5)
      close(3)


c  !KJ Next section added to read ELNES variables     1-06  
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

      return
      end
