       subroutine feffdt(ntotal,iplst,nptot,ntext,text,ne,
     $      iorder, l0, rnrmav, edge, potlbl,
     $      iz,phc,ck,xk,index,
     $      nleg,deg,reff,crit,ipot,rat,achi,phchi)
c
c     removed ihole, neptsm, npot, xmu from arg list
c
c     writes feffnnnn.dat files and files.dat
c     for compatibility with the old feff
c
      implicit double precision (a-h, o-z)

      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'
      include '../HEADERS/vers.h'
      parameter (eps4 = 1.0d-4)
      parameter (eps = 1.0d-16)

      parameter (npx=15000)
      character*12 fname(npx)
      character*512 slog
      dimension iplst(npx)

c     Stuff from feff.pad, note that floating point numbers are
c     single precision
cc      character*78 string
      real rnrmav, edge
c      real xmu
cc      dimension ltext(nheadx)
      character*80 text(nheadx), lines(2*nheadx)
      character*6  potlbl(0:nphx)
      dimension iz(0:nphx)
c     central atom phase shift at l0
      complex phc(nex)
      complex ck(nex)
      real xk(nex)
      dimension index(npx)
      dimension nleg(npx)
      real deg(npx), reff(npx), crit(npx)
      dimension ipot(legtot,npx)
      real rat(3,legtot,npx)
cc      real beta(legtot,npx)
cc      real eta(legtot,npx)
cc      real ri(legtot,npx)
      real achi(nex,npx), phchi(nex,npx)
      integer istrln
c      complex*16 cchi, cfms
      dimension col1(nex), col2(nex), col3(nex), col4(nex), col5(nex)
      dimension col6(nex), col7(nex)
      real amff(nex), phff(nex)
      dimension ipoths(legtot)
      dimension thsrat(3,legtot)


      external istrln

      call wlog (' feffdt, feff.pad to feff.dat conversion: ' // vfeff
     1            // 'release ' // vf85e)

c     read feff.pad
c     Use single precision for all fp numbers in feff.pad
      do 20  itext = 1, ntext
         ltxt = istrln(text(itext))
c        text(itext) does not have carriage control
         call wlog (' ' // text(itext)(1:ltxt))
   20 continue

      write(slog,60)  nptot
   60 format (i8, ' paths to process')
      call wlog (slog)

c     make files.dat
  160 format (1x, a)
  170 format (1x, 71('-'))

c     Save filenames of feff.dat files
      open (unit=2, file='files.dat', status='unknown', iostat=ios)
      call chopen (ios, 'files.dat', 'genfmt')
c     Put phase header on top of files.dat
      do 200  itext = 1, ntext
         ltxt = istrln( text(itext))
         write(2,160)  text(itext)(1:ltxt)
  200 continue
      write(2,170)
      write(2,210)
  210 format ('    file        sig2   amp ratio    ',
     1        'deg    nlegs  r effective')
c     do each path
      call wlog ('    path     filename')

      do 500  ilist = 1, ntotal
c        find index of path
         do 410  j = 1, nptot
            if (iplst(ilist) .eq. index(j))  then
               i = j
               goto 430
            endif
  410    continue
         write(slog,420)  ilist, iplst(ilist)
  420    format (' did not find path i, iplst(i) ', 2i10)
         call wlog(slog)
  430    continue
c        Path i is the path from feff.pad that corresponds to
c        the path ilist in list.dat.  The index of the path is
c        iplst(ilist) and index(i).

c        Prepare output file feffnnnn.dat
         write(fname(i),220)  index(i)
  220    format ('feff', i4.4, '.dat')
         write(slog,230)  i, fname(i)
  230    format (i8, 5x, a)
         call wlog(slog)
c        zero is debye-waller factor column
         write(2,240) fname(i), zero, crit(i), deg(i),
     1                   nleg(i), reff(i)*bohr
  240    format(1x, a, f8.5, 2f10.3, i6, f9.4)

         ip = i
c     Write feff.dat's
         open (unit=3, file=fname(ip), status='unknown', iostat=ios)
         call chopen (ios, fname(ip), 'feffdt')
c        put header on feff.dat
c$$$         do 300  itext = 1, ntext
c$$$            ltxt = istrln(text(itext))
c$$$            write(3,160)  text(itext)(1:ltxt)
c$$$  300    continue
c$$$
c$$$         write(3,310) ip, iorder
c$$$  310    format (' Path', i5, '      icalc ', i7)
c$$$         write(3,170)
c$$$         write(3,320)  nleg(ip), deg(ip), reff(ip)*bohr, rnrmav,
c$$$     1                 edge*hart
c$$$  320    format (1x, i3, f8.3, f9.4, f10.4, f11.5,
c$$$     1           ' nleg, deg, reff, rnrmav(bohr), edge')
c$$$         write(3,330)
c$$$  330    format ('        x         y         z   pot at#')
c$$$         write(3,340)  (rat(j,nleg(ip),ip)*bohr,j=1,3),
c$$$     1                 ipot(nleg(ip),ip),
c$$$     1                 iz(ipot(nleg(ip),ip)), potlbl(ipot(nleg(ip),ip))
c$$$  340    format (1x, 3f10.4, i3, i4, 1x, a6, '   absorbing atom')
c$$$         do 360  ileg = 1, nleg(ip)-1
c$$$            write(3,350)  (rat(j,ileg,ip)*bohr,j=1,3), ipot(ileg,ip),
c$$$     1                    iz(ipot(ileg,ip)), potlbl(ipot(ileg,ip))
c$$$  350       format (1x, 3f10.4, i3, i4, 1x, a6)
c$$$  360    continue
c$$$
c$$$         write(3,370)
c$$$  370    format    ('    k   real[2*phc]   mag[feff]  phase[feff]',
c$$$     1              ' red factor   lambda     real[p]@#')

         do 900 nl = 1, nleg(ip)
            ipoths(nl) = ipot(nl, ip)
            do 910 i = 1, 3
               thsrat(i, nl) = rat(i, nl, ip)
 910        continue
 900     continue
         call fdthea(ntext, text, ip, iorder, nleg(ip), deg(ip),
     &          reff(ip), rnrmav, edge, thsrat, ipoths, iz, potlbl,
     &          nlines, lines)
         do 920 i=1, nlines
            ltxt = istrln(lines(i))
            write(3, 930) lines(i)(1:ltxt)
 920     continue
 930     format(a)

         do 1000 ie = 1, ne
            amff(ie) = achi(ie,ip)
            phff(ie) = phchi(ie,ip)
 1000    continue
         call fdtarr(ne, reff(ip), l0, amff, phff, phc, xk, ck,
     &       col1, col2, col3, col4, col5, col6, col7)
         do 1005 ie = 1, ne
            write(3,400) col1(ie), col2(ie), col3(ie), col4(ie),
     &             col5(ie), col6(ie), col7(ie)
 1005    continue

c        Make the feff.dat stuff and write it to feff.dat
c        Also write out for inspection to fort.66
c        note that dimag takes complex*16 argument, aimag takes
c        single precision complex argument.  Stuff from feff.pad
c        is single precision, cchi is complex*16
c$$$         do 450  ie = 1, ne
c$$$c           Consider chi in the standard XAFS form.  Use R = rtot/2.
c$$$            cchi = achi(ie,ip) * exp (coni*phchi(ie,ip))
c$$$            xlam = 1.0e10
c$$$            if (abs(aimag(ck(ie))) .gt. eps) xlam = 1/aimag(ck(ie))
c$$$            redfac = exp (-2 * aimag (phc(ie)))
c$$$            cdelt = 2*dble(phc(ie))
c$$$            cfms = cchi * xk(ie) * reff(ip)**2 *
c$$$     1           exp(2*reff(ip)/xlam) / redfac
c$$$            if (abs(cchi) .lt. eps)  then
c$$$               phff = 0
c$$$            else
c$$$               phff = atan2 (dimag(cchi), dble(cchi))
c$$$            endif
c$$$c           remove 2 pi jumps in phases
c$$$            if (ie .gt. 1)  then
c$$$               call pijump (phff, phffo)
c$$$               call pijump (cdelt, cdelto)
c$$$            endif
c$$$            phffo = phff
c$$$            cdelto = cdelt
c$$$c           write 1 k, momentum wrt fermi level k=sqrt(p**2-kf**2)
c$$$c                 2 central atom phase shift (real part),
c$$$c                 3 magnitude of feff,
c$$$c                 4 phase of feff,
c$$$c                 5 absorbing atom reduction factor,
c$$$c                 6 mean free path = 1/(Im (p))
c$$$c                 7 real part of local momentum p
c$$$
c$$$            write(3,400)
c$$$     1         xk(ie)/bohr,
c$$$     2         cdelt + l0*pi,
c$$$     3         abs(cfms) * bohr,
c$$$     4         phff - cdelt - l0*pi,
c$$$     5         redfac,
c$$$     6         xlam * bohr,
c$$$     7         dble(ck(ie))/bohr
c$$$
c$$$  450    continue

 400     format (1x, f6.3, 1x, 3(1pe11.4,1x),1pe10.3,1x,
     1          2(1pe11.4,1x))

c        Done with feff.dat
         close (unit=3)
  500 continue
      close (unit=2)

      return
      end
