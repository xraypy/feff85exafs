      subroutine rexsph ( mphase, ipr2, ispec, vixan, xkstep, xkmax,
     1             gamach, rgrd,
     1             nph, lmaxph, potlbl, spinph, iatph, nat, rat, iphat,
     2             ixc, vr0, vi0, ixc0, lreal, rfms2, lfms2, l2lp,
     3             ipol, ispin, le2, angks, ptz, iPl, iGrid,
     4             izstd, ifxc, ipmbse, itdlda, nonlocal, ibasis)

      implicit double precision (a-h, o-z)

      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'

cc    geom.dat
        integer  nat, iatph(0:nphx), iphat(natx)
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
  55          format(a, 2i10)
              stop 'Bad input'
           endif
           read(3,*,end=60)  idum, (rat(j,nat),j=1,3), iphat(nat), i1b
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
          ptz(-1,i)= cmplx(a1,b1)
          ptz(0,i) = cmplx(a2,b2)
          ptz(1,i) = cmplx(a3,b3)
  90    continue
      close(3)
c     read mod2.inp
c     Josh - added flag iPl for PLASMON card
c     Josh - added flag iGrid for user controlled grids.
      open (file='mod2.inp', unit=3, status='old',iostat=ios)
        read (3,10)  slog
        read (3,20)  mphase,ipr2,ixc,ixc0,ispec,lreal,lfms2,nph,l2lp,
     &       iPl,iGrid
        read (3,10)  slog
        read (3,30)  vr0, vi0
        read (3,10)  slog
        read (3,20)  (lmaxph(iph),iph=0,nph)
        read (3,10)  slog
        read (3,170)  (potlbl(iph),iph=0,nph)
  170   format (13a6)
        read (3,10)  slog
        read (3,30)  rgrd, rfms2, gamach, xkstep, xkmax, vixan
        read (3,30)  (spinph(iph),iph=0,nph)
        read (3,20)  izstd, ifxc, ipmbse, itdlda, nonlocal, ibasis
      close(3)

!KJ next section added for ELNES calculations 1-06
      open(3,file='eels.inp',status='old',err=100)
        read(3,*)
	read(3,20,end=100,err=100) melnes
      close(3)
      if(melnes.eq.1.and.mphase.eq.1) then
        call wlog(':INFO : rexsph reduces your polarization tensor to 
     1   the unit matrix, because eels.inp says you are doing ELNES.')
        do i=-1,1
	do j=-1,1
	ptz(i,j)=dcmplx(0,0)
	enddo
	  ptz(i,i)=dble(1)/dble(3)
	  write(*,*) (ptz(i,j),j=-1,1)
	enddo
      endif
100   continue
c !KJ end of my modifications      



c     transform to code units (bohrs and hartrees - atomic unuts)
      rfms2 = rfms2 / bohr
      vr0   = vr0 / hart
      vi0   = vi0 / hart
      gamach = gamach / hart
      vixan = vixan / hart
      xkstep = xkstep * bohr
      xkmax  = xkmax  * bohr
      do 210 i = 1,3
      do 210 iat = 1, nat
        rat(i,iat) = rat(i,iat) / bohr
 210  continue

      return
      end
