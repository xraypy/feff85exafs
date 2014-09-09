      subroutine reafms ( mfms, rfms2, idwopt, tk, thetad, sig2g,
     1            lmaxph, nat, iphat, rat,
     2            ipol, ispin, le2, angks, ptz,
     3            minv, rdirec, toler1, toler2,
     4            elnes,ipmin,ipmax,ipstep)   !KJ added this line 1-06     

      use json_module
      implicit double precision (a-h, o-z)
      logical :: found
      type(json_file) :: json   !the JSON structure read from the file:
      double precision toss
      integer,dimension(:),allocatable :: intgs

      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'
      
      integer elnes,ipmin,ipmax,ipstep  !KJ added these variables 1-06

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
cc    mod3.inp
        integer mfms, idwopt, minv
        integer lmaxph(0:nphx)
        real rfms2, rdirec, toler1, toler2
        double precision   tk, thetad, sig2g

c     Local stuff
c      character*512 slog
c      character*80 head(nheadx)
c      dimension lhead(nheadx)

c     standard formats for string, integers and real numbers
  10  format(a)
  20  format (20i4)
  30  format (6f13.5)

      call json_read_geom(nat, nph, iatph, rat, iphat, ibounc)

      call json_read_global(nabs, iphabs, rclabs, ipol, ispin, le2,
     1                      elpty, angks, evec, xivec, spvec, ptz)


c$$$c     Read  geom.dat file
c$$$      open (file='geom.dat', unit=3, status='old',iostat=ios)
c$$$c       read header from geom.dat
c$$$        nhead = nheadx
c$$$        call rdhead (3, nhead, head, lhead)
c$$$        nat = 0
c$$$        nph = 0
c$$$        do 40 iph = 0, nphx
c$$$  40    iatph(iph) = 0
c$$$  50    continue
c$$$           nat = nat+1
c$$$           if (nat .gt. natx)  then
c$$$              write(slog,55) ' nat, natx ', nat, natx
c$$$              call wlog(slog)
c$$$  55          format(a, 2i10)
c$$$              stop 'Bad input'
c$$$           endif
c$$$           read(3,*,end=60)  idum, (rat(j,nat),j=1,3), iphat(nat), i1b
c$$$           if (iphat(nat).gt.nph) nph = iphat(nat)
c$$$           if ( iatph(iphat(nat)).eq.0) iatph(iphat(nat)) = nat
c$$$        goto 50
c$$$  60    continue
c$$$        nat = nat-1
c$$$      close(3)
c$$$
c$$$cc    global.inp
c$$$      open (file='global.dat', unit=3, status='unknown',iostat=ios)
c$$$c       configuration average data
c$$$        read  (3, 10) slog
c$$$        read  (3, 65) nabs, iphabs, rclabs
c$$$  65    format ( 2i8, f13.5)
c$$$c       global polarization data
c$$$        read  (3,10)   slog
c$$$        read  (3, 70)  ipol, ispin, le2, elpty, angks
c$$$  70    format ( 3i5, 2f12.4)
c$$$        read  (3, 10) slog
c$$$        do 80 i = 1,3
c$$$          read  (3,30) evec(i), xivec(i), spvec(i)
c$$$  80    continue
c$$$        read  (3, 10) slog
c$$$        do 90 i = -1, 1
c$$$          read (3,30) a1, b1, a2, b2, a3, b3
c$$$          ptz(-1,i)= cmplx(a1,b1)
c$$$          ptz(0,i) = cmplx(a2,b2)
c$$$          ptz(1,i) = cmplx(a3,b3)
c$$$  90    continue
c$$$      close(3)
c$$$c     read mod3.inp
c$$$      open (file='mod3.inp', unit=3, status='old',iostat=ios)
c$$$        read (3,10)  slog
c$$$        read (3,20)  mfms, idwopt, minv
c$$$        read (3,10)  slog
c$$$        read (3,30)  rfms2, rdirec, toler1, toler2
c$$$        read (3,10)  slog
c$$$        read (3,30)  tk, thetad, sig2g
c$$$        read (3,10)  slog
c$$$        read (3,20)  (lmaxph(iph),iph=0,nph)
c$$$      close(3)

      call json%load_file('fms.json')
      if (json_failed()) then   !if there was an error reading the file
         print *, "failed to read fms.json"
         stop
      else
         call json%get('mfms',   mfms, found)
                   if (.not. found) call bailout('mfms',   'fms.json')
         call json%get('idwopt', idwopt, found)
                   if (.not. found) call bailout('idwopt', 'fms.json')
         call json%get('minv',   minv, found)
                   if (.not. found) call bailout('minv',   'fms.json')

         call json%get('rfms2',  toss, found)
                   if (.not. found) call bailout('rfms2',  'fms.json')
         rfms2  = real(toss)
         call json%get('rdirec', toss, found)
                   if (.not. found) call bailout('rdirec', 'fms.json')
         rdirec = real(toss)
         call json%get('toler1', toss, found)
                   if (.not. found) call bailout('toler1', 'fms.json')
         toler1 = real(toss)
         call json%get('toler2', toss, found)
                   if (.not. found) call bailout('toler2', 'fms.json')
         toler2 = real(toss)

         call json%get('tk',     tk, found)
                   if (.not. found) call bailout('tk',     'fms.json')
         call json%get('thetad', thetad, found)
                   if (.not. found) call bailout('thetad', 'fms.json')
         call json%get('sig2g',  sig2g, found)
                   if (.not. found) call bailout('sig2g',  'fms.json')

         call json%get('lmaxph', intgs, found)
                   if (.not. found) call bailout('lmaxph', 'fms.json')
         do 1000 iph = 0, nphx
            lmaxph(iph) = intgs(iph+1)            
 1000    continue

         call json%destroy()
      end if


      
c$$$c  !KJ Next section added to read ELNES variables 1-06     
c$$$c     read eels.inp
c$$$      elnes=0
c$$$      open(file='eels.inp',unit=3,status='old',err=900)
c$$$        read(3,*,err=900,end=900) 
c$$$        read(3,20,err=900,end=900) elnes
c$$$        read(3,*,err=900,end=900)
c$$$        read(3,*,err=900,end=900)
c$$$        read(3,*,err=900,end=900)
c$$$        read(3,20,err=900,end=900) ipmin,ipstep,ipmax
c$$$      close(3)
c$$$      goto 901
c$$$900   continue
c$$$      elnes=0
c$$$901   continue
c$$$      if(elnes.eq.0) then
c$$$        ipstep=1
c$$$        ipmax=1
c$$$        ipmin=1
c$$$      endif
c$$$               
c$$$c  !KJ end my changes

c     transform to code units (bohrs and hartrees - atomic units)
      rfms2 = rfms2 / real(bohr)
      rdirec = rdirec / real(bohr)
      do 210 iat = 1, nat
      do 210 i = 1,3
        rat(i,iat) = rat (i, iat) / bohr
 210  continue

      return
      end
