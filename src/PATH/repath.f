      subroutine repath ( ms, mpath, ipr4, pcritk, pcrith, nncrit, rmax,
     1             nlegxx, rfms2, critpw,
     2             nat, rat, iphat, ibounc,
     3             ipol, ispin, evec, xivec ,eels)  !KJ added eels 5/06

      use json_module
      implicit double precision (a-h, o-z)
      logical :: found
      type(json_file) :: json   !the JSON structure read from the file:
      double precision toss

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
cc    mod4.inp
        integer  mpath, ms, nncrit, nlegxx, ipr4
        real critpw, pcritk, pcrith,  rmax, rfms2
        
        integer eels !KJ added 5/06


c     Local stuff
      character*512 slog
      character*80 head(nheadx)
      dimension lhead(nheadx)

c     standard formats for string, integers and real numbers
  10  format(a)
  20  format (20i4)
  30  format (6f13.5)


c     Read  geom.dat file
      open (file='geom.dat', unit=3, status='old',iostat=ios)
c       read header from geom.dat
        nhead = nheadx
        call rdhead (3, nhead, head, lhead)
        nat = 0
        nph = 0
        do 40 iph = 0, nphx
  40    iatph(iph) = 0
  50    continue
         nat = nat+1
         if (nat .gt. natx)  then
           write(slog,55) ' nat, natx ', nat, natx
           call wlog(slog)
  55       format(a, 2i10)
           stop 'Bad input'
         endif
         read(3,*,end=60) idum,(rat(j,nat),j=1,3),iphat(nat),ibounc(nat)
         if (iphat(nat).gt.nph) nph = iphat(nat)
         if ( iatph(iphat(nat)).eq.0) iatph(iphat(nat)) = nat
        goto 50
  60    continue
        nat = nat-1
      close(3)

cc    global.inp
      open (file='global.dat', unit=3, status='unknown',iostat=ios)
c       configuration average data
        read  (3, 10) slog
        read  (3, 65) nabs, iphabs, rclabs
  65    format ( 2i8, f13.5)
c       global polarization data
        read  (3,10)   slog
        read  (3, 70)  ipol, ispin, le2, elpty, angks
  70    format ( 3i5, 2f12.4)
        read  (3, 10) slog
        do 80 i = 1,3
          read  (3,30) evec(i), xivec(i), spvec(i)
  80    continue
        read  (3, 10) slog
        do 90 i = -1, 1
          read (3,30) a1, b1, a2, b2, a3, b3
          ptz(-1,i)= dcmplx(a1,b1)
          ptz(0,i) = dcmplx(a2,b2)
          ptz(1,i) = dcmplx(a3,b3)
  90    continue
      close(3)
c--json--c     read mod4.inp
c--json--      open (file='mod4.inp', unit=3, status='old',iostat=ios)
c--json--        read (3,10)  slog
c--json--        read (3,20)  mpath, ms, nncrit, nlegxx, ipr4
c--json--        read (3,10)  slog
c--json--        read (3,30)  critpw, pcritk, pcrith,  rmax, rfms2
c--json--      close(3)

      call json%load_file('path.json')
      if (json_failed()) then   !if there was an error reading the file
         print *, "failed to read path.json"
         stop
      else
         call json%get('mpath',   mpath, found)
                   if (.not. found) call bailout('mpath', 'path.json')
         call json%get('ms',   ms, found)
                   if (.not. found) call bailout('ms', 'path.json')
         call json%get('nncrit',   nncrit, found)
                   if (.not. found) call bailout('nncrit', 'path.json')
         call json%get('nlegxx',   nlegxx, found)
                   if (.not. found) call bailout('nlegxx', 'path.json')
         call json%get('ipr4',   ipr4, found)
                   if (.not. found) call bailout('ipr4', 'path.json')
         call json%get('critpw',   toss, found)
                   if (.not. found) call bailout('critpw', 'path.json')
         critpw = real(toss)
         call json%get('pcritk', toss, found)
                   if (.not. found) call bailout('pcritk', 'path.json')
         pcritk = real(toss)
         call json%get('pcrith', toss, found)
                   if (.not. found) call bailout('pcrith', 'path.json')
         pcrith = real(toss)
         call json%get('rmax',   toss, found)
                   if (.not. found) call bailout('rmax', 'path.json')
         rmax = real(toss)
         call json%get('rfms2',  toss, found)
                   if (.not. found) call bailout('rfms2', 'path.json')
         rfms2 = real(toss)
      end if
      
c  !KJ Next section added to read EELS variables       5-06
c     read eels.inp
      eels=0
      open(file='eels.inp',unit=3,status='old',err=900)
        read(3,*,err=900,end=900) 
        read(3,20,err=900,end=900) eels
      close(3)
      goto 901
900   continue
      eels=0
901   continue
c  !KJ end my changes 
      

      return
      end
