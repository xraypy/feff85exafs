      subroutine reff2x(mchi, ispec, ipr6, idwopt, critcw, s02, sig2g,
     1                   tk, thetad, mbconv, absolu,  !KJ added absolu 3-06
     1                   vrcorr, vicorr, alphat, thetae, iabs, nabs,
     4            elnes,ipmin,ipmax,ipstep)   !KJ added this line  1-06     

      implicit double precision (a-h, o-z)

      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'

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
