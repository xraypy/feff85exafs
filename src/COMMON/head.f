      subroutine sthead (ntitle, title, nph, iz, rmt, rnrm,
     1     xion, ihole, ixc,
     2     vr0, vi0, gamach, xmu, xf, vint, rs,
     2     lreal,  rgrd)

c     SeT HEAD
c     This routine makes the file header, returned in head array.
c     header lines do not include a leading blank.
c     Last header line is not --------- end-of-header line

c     title lines coming into sthead include carriage control, since
c     they were read from potph.bin

      implicit double precision (a-h, o-z)

      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'
      include '../HEADERS/vers.h'

      dimension xion(0:nphx)
      dimension iz(0:nphx)
      dimension rmt(0:nphx)
      dimension rnrm(0:nphx)

      character*80 title(nheadx), store
      character*16 s1, s2

      character*10 shole(0:29)
      character*8  sout(0:7)
      data shole /'no hole',   'K  shell',  'L1 shell',  'L2 shell',
     2            'L3 shell',  'M1 shell',  'M2 shell',  'M3 shell',
     3            'M4 shell',  'M5 shell',  'N1 shell',  'N2 shell',
     4            'N3 shell',  'N4 shell',  'N5 shell',  'N6 shell',
     5            'N7 shell',  'O1 shell',  'O2 shell',  'O3 shell',
     6            'O4 shell',  'O5 shell',  'O6 shell',  'O7 shell',
     7            'P1 shell',  'P2 shell',  'P3 shell',  'P4 shell',
     8            'P5 shell',  'R1 shell'/
      data sout /'H-L exch', 'D-H exch', 'Gd state', 'DH - HL ',
     1           'DH + HL ', 'val=s+d ', 'sigmd(r)', 'sigmd=c '/


c     Fills head arrray, n = number of lines used.
c     Does not include line of dashes at the end.

      if (ntitle .ge. 1 ) then
         ii = istrln(title(1)) 
         if (ii.gt.1)  then
            write(store,100)  title(1)(1:59), vfeff//vf85e
         else
            write(store,102)  vfeff//vf85e
         endif
      else
         write(store,102)   vfeff//vf85e
      endif
  100 format( a45, t48, a30)
  102 format( t48, a30)
      title(1) = store
      nstor = 1

c     remove empty title lines
      do 120  ititle = 2, ntitle
         ii = istrln ( title (ititle) ) 
         if (ii.le.1)  goto 120
         nstor = nstor+1
         title(nstor) = title (ititle)
  120 continue
      ntitle = nstor

c     add more title lines
      if (xion(0) .ne. 0)  then
         ntitle = ntitle + 1
         write(title(ntitle),130)  iz(0), rmt(0)*bohr,
     1                    rnrm(0)*bohr, xion(0), shole(ihole)
      else
         ntitle = ntitle + 1
         write(title(ntitle),140)  iz(0), rmt(0)*bohr,
     1                    rnrm(0)*bohr, shole(ihole)
      endif
  130 format('Abs   Z=',i2, ' Rmt=',f6.3, ' Rnm=',f6.3,
     1       ' Ion=',f5.2,  1x,a10)
  140 format('Abs   Z=',i2,' Rmt=',f6.3,' Rnm=',f6.3, 1x,a10)
      if (lreal.ge.1 .or. (abs(rgrd - 0.05) .gt. 1.0e-5)) then
        ntitle = ntitle + 1
        s1 = ' '
        if (lreal.gt.1)  then
c        write(title(ntitle),144)
c 144    format ('Calculations done using only real phase shifts.')
         s1 = 'RPHASES'
        elseif (lreal.eq.1) then
c        ntitle = ntitle + 1
c        write(title(ntitle),145)
c 145    format ('Calculations done using only real self energy.')
         s1 = 'RSIGMA'
        endif
        s2 = '  '
        if (abs(rgrd - 0.05) .gt. 1.0e-5)  then
         write(s2,146)  rgrd
  146    format ('  RGRID', f7.4)
        endif
        ilen = istrln(s1)
        title(ntitle) = s1(1:ilen) // s2
      endif

      do 150  iph = 1, nph
         if (xion(iph) .ne. 0)  then
            ntitle = ntitle + 1
            write(title(ntitle),160)  iph, iz(iph),  rmt(iph)*bohr,
     1           rnrm(iph)*bohr, xion(iph)
         else
            ntitle = ntitle + 1
            write(title(ntitle),170)  iph, iz(iph),  rmt(iph)*bohr,
     1           rnrm(iph)*bohr
         endif
  150 continue
  160 format('Pot',i2,' Z=',i2,' Rmt=',f6.3,' Rnm=',f6.3,' Ion=',f5.2)
  170 format('Pot',i2,' Z=',i2,' Rmt=',f6.3,' Rnm=',f6.3)
      if (abs(vi0) .gt. 1.0e-8 .or. abs(vr0) .gt. 1.0e-8)  then
         ntitle = ntitle + 1
         write(title(ntitle),180)  gamach*hart, sout(ixc), vi0*hart,
     1                           vr0*hart
      else
         ntitle = ntitle + 1
         write(title(ntitle),190)  gamach*hart, sout(ixc)
      endif
      ntitle = ntitle + 1
  180 format('Gam_ch=',1pe9.3, 1x,a8, ' Vi=',1pe10.3, ' Vr=',1pe10.3)
  190 format('Gam_ch=',1pe9.3, 1x,a8)
  200 format('Mu=',1pe10.3, ' kf=',1pe9.3, ' Vint=',1pe10.3,
     x        ' Rs_int=',0pf6.3)
      write(title(ntitle),200)  xmu*hart, xf/bohr, vint*hart, rs

      return
      end

      subroutine wthead (io, ntitle, title)
c     Dump title lines to unit io, which must be open. 
      integer io, i, ll
      character*80 title(ntitle)

c     nice for UNIX to use with gnuplot etc.,
      do 310 i = 1, ntitle
         ll = istrln(title(i))
         write(io,300)  title(i)(1:ll)
  300    format (a)
  310 continue

      return
      end
