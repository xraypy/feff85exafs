      subroutine wrtall (nabs)
c     writes data stored in common blocks of allinp.h to 
c     all necessary input files for other modules.
c     version 1.0 written by Alexei Ankudinov, March 2001

c     Note: to add input variable one has to add it to the 
c        appropriate common block in allinp.h, properly initialize
c        it in subroutine iniall and modify subroutine wrtall
c        to write it to the appropriate input file.
c        (i.e. one has to make modifications in 3 places)

      implicit double precision (a-h, o-z)
      include '../HEADERS/dim.h'
      include '../HEADERS/parallel.h'
      include '../RDINP/allinp.h'

      if (.not. master) return

c     standard formats for string, integers and real numbers
  10  format(a)
  20  format (20i4)
  30  format (6f13.5)

cc    atoms.dat to be read by ffsort,
cc    that will write smaller geom.dat file
      open (file='atoms.dat', unit=3, status='unknown',iostat=ios)
        write (3, 35) natt
  35    format ('natx =  ', i7)
        write (3, 10) '    x       y        z       iph  '
        do 40  iat = 1, natt
          write(3,36) ratx(1,iat), ratx(2,iat), ratx(3,iat), iphatx(iat)
  36      format( 3f13.5, i4)
  40    continue
      close(3)

cc    global.inp
      open (file='global.dat', unit=3, status='unknown',iostat=ios)
c       configuration average data
        write (3, 10) ' nabs, iphabs - CFAVERAGE data'
        write (3, 45) nabs, iphabs, rclabs
  45    format ( 2i8, f13.5)
c       global polarization data
        write (3,10) ' ipol, ispin, le2, elpty, angks'
        write (3, 50)  ipol, ispin, le2, elpty, angks
  50    format ( 3i5, 2f12.4)
        write (3, 10) 'evec         xivec        spvec'
        do 60 i = 1,3
          write (3,30) evec(i), xivec(i), spvec(i)
  60    continue
        write (3, 10) ' polarization tensor '
        do 70 i = -1, 1
          write(3,30) dble(ptz(-1,i)), dimag(ptz(-1,i)), dble(ptz(0,i)),
     1                dimag(ptz(0,i)),  dble(ptz(1,i)), dimag(ptz(1,i))
  70    continue
      close(3)
        
cc    mod1.inp
      open (file='mod1.inp', unit=3, status='unknown',iostat=ios)
        write(3,10) 'mpot, nph, ntitle, ihole, ipr1, iafolp, ixc,ispec'
        write(3,20) mpot, nph, ntitle, ihole, ipr1, iafolp, ixc, ispec
        write(3,10) 
     1  'nmix, nohole, jumprm, inters, nscmt, icoul, lfms1, iunf'
        write(3,20)  nmix, nohole, jumprm, inters, nscmt, icoul, lfms1,
     1   iunf
        do 110 ititle = 1, ntitle
  110   write(3,10) title(ititle)
        write(3,10) 'gamach, rgrd, ca1, ecv, totvol, rfms1'
        write(3,30)  gamach, rgrd, ca1, ecv, totvol, rfms1
        write(3,10) ' iz, lmaxsc, xnatph, xion, folp'
  120   format ( 2i5, 4f13.5)
        do 130 ip = 0, nph
  130   write(3,120) iz(ip), lmaxsc(ip), xnatph(ip), xion(ip), folp(ip)
c       for OVERLAP option
        write(3,10) 'OVERLAP option: novr(iph)'
        write(3,20) ( novr(iph), iph=0,nph)
        write(3,10) ' iphovr  nnovr rovr '
  140   format ( 2i5, f13.5)
        do 150 iph = 0, nph
        do 150 iovr = 1, novr(iph)
  150   write(3,140) iphovr(iovr, iph), nnovr(iovr,iph), rovr(iovr,iph)
      close(3)

cc    mod2.inp
      open (file='mod2.inp', unit=3, status='unknown',iostat=ios)
c     Josh - added flag for PLASMON card (iPlsmn = 0, 1, or 2)
!     Josh - added flag for user difined grid (EGRID card).
        write(3,10) 'mphase,ipr2,ixc,ixc0,ispec,lreal,lfms2,nph,l2lp,
     &     iplsmn,igrid'
        write(3,20)  mphase,ipr2,ixc,ixc0,ispec,lreal,lfms2,nph,l2lp,
     &        iPlsmn, iGrid
        write(3,10) 'vr0, vi0'
        write(3,30)  vr0, vi0
        write(3,10) ' lmaxph(0:nph)'
        write(3,20)  (lmaxph(iph),iph=0,nph)
        write(3,10) ' potlbl(iph)'
        write(3,170)  (potlbl(iph),iph=0,nph)
  170   format (13a6)
        write(3,10) 'rgrd, rfms2, gamach, xkstep, xkmax, vixan'
        write(3,30)  rgrd, rfms2, gamach, xkstep, xkmax, vixan
        write(3,30)  (spinph(iph),iph=0,nph)
        write(3,20)  izstd, ifxc, ipmbse, itdlda, nonlocal, ibasis
      close(3)

cc    mod3.inp
      open (file='mod3.inp', unit=3, status='unknown',iostat=ios)
        write(3,10) 'mfms, idwopt, minv'
        write(3,20)  mfms, idwopt, minv
        write(3,10) 'rfms2, rdirec, toler1, toler2'
        write(3,30)  rfms2, rdirec, toler1, toler2
        write(3,10) 'tk, thetad, sig2g'
        write(3,30)  tk, thetad, sig2g
        write(3,10) ' lmaxph(0:nph)'
        write(3,20)  (lmaxph(iph),iph=0,nph)
      close(3)

cc    mod4.inp
      open (file='mod4.inp', unit=3, status='unknown',iostat=ios)
        write(3,10) 'mpath, ms, nncrit, nlegxx, ipr4'
        write(3,20)  mpath, ms, nncrit, nlegxx, ipr4
        write(3,10) 'critpw, pcritk, pcrith,  rmax, rfms2'
        write(3,30)  critpw, pcritk, pcrith,  rmax, rfms2
      close(3)

cc    mod5.inp
      open (file='mod5.inp', unit=3, status='unknown',iostat=ios)
        write(3,10) 'mfeff, ipr5, iorder, critcw, wnstar'
        write(3,180)  mfeff, ipr5, iorder, critcw, wnstar
  180   format ( 2i4, i8, f13.5, L5)
      close(3)

cc    mod6.inp
      open (file='mod6.inp', unit=3, status='unknown',iostat=ios)
        write(3,10) 'mchi, ispec, idwopt, ipr6, mbconv, absolu' !KJ added absolu 3-06
        write(3,20)  mchi, ispec, idwopt, ipr6, mbconv, absolu !KJ added absolu 3-06
        write(3,10) 'vrcorr, vicorr, s02, critcw'
        write(3,30)  vrcorr, vicorr, s02, critcw
        write(3,10) 'tk, thetad, alphat, thetae, sig2g'
        write(3,30)  tk, thetad, alphat, thetae, sig2g
      close(3)

c$$$cc    so2.inp - Josh Kas
c$$$      open (file='s02.inp', unit=3, status='unknown',iostat=ios)
c$$$        write(3,10) 'mso2conv, ipse, ipsk'
c$$$        write(3,20)  mso2conv, ipse, ipsk
c$$$        write(3,10) 'wsigk, cen'
c$$$        write(3,30) wsigk, cen
c$$$        write(3,10) 'ispec, ipr6'
c$$$        write(3,20)  ispec, ipr6
c$$$        write(3,10) 'cfname'
c$$$        write(3,10) cfname
c$$$      close(3)
      return
      end
