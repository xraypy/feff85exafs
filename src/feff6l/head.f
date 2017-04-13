      subroutine sthead (ntitle, title, ltitle, nph, iz, rmt, rnrm,
     1                  ion, ifrph, ihole, ixc,
     2                  vr0, vi0, rs0, gamach, xmu, xf, vint, rs,
     3                  nhead, lhead, head, shole, sout)

c     SeT HEAD
c     This routine makes the file header, returned in head array.
c     header lines do not include a leading blank.
c     Last header line is not --------- end-of-header line

c     title lines coming into sthead include carriage control, since
c     they were read from potph.dat

      implicit double precision (a-h, o-z)

      include 'const.h'
      include 'dim.h'

      dimension ifrph(0:nphx)
      dimension ion(0:nfrx)
      dimension iz(0:nfrx)
      dimension rmt(0:nphx)
      dimension rnrm(0:nphx)

      character*80 title(ntitle)
      parameter (nheadx = 30)
      character*80 head(nheadx)
      dimension lhead(nheadx), ltitle(ntitle)

      character*10 shole(0:9)
      character*8 sout(0:6)


      character*80 heada(nheadx)
      dimension lheada(nheadx)
      save nheada, lheada, heada
c     heada, etc., are saved for use by entry wthead

cc       common /labels/ shole, sout
      include 'vers.h'
c     character*12 vfeff, vpotph, vpaths, vgenfm, vff2ch
c     common /vers/ vfeff, vpotph, vpaths, vgenfm, vff2ch

c     FiLl head array with HEADer
c     Fills head arrray, n = number of lines used.
c     Does not include line of dashes at the end.

      nhead = 1
      if (ntitle .ge. 1  .and.  ltitle(1).gt.1)  then
         write(head(nhead),100)  title(1)(2:), vfeff, vpotph
      else
         write(head(nhead),102)  vfeff, vpotph
      endif
  100 format(a55, t56, 2a12)
  102 format(t56, 2a12)
      do 120  ititle = 2, ntitle
         if (ltitle(ititle).le.1)  goto 120
         nhead = nhead+1
         write(head(nhead),110) title(ititle)(2:)
  110    format(a79)
  120 continue
      if (ion(0) .ne. 0)  then
         nhead = nhead+1
         write(head(nhead),130)  iz(0), rmt(0)*bohr,
     1                    rnrm(0)*bohr, ion(0), shole(ihole)
      else
         nhead = nhead+1
         write(head(nhead),140)  iz(0), rmt(0)*bohr,
     1                    rnrm(0)*bohr, shole(ihole)
      endif
  130 format('Abs   Z=',i2,' Rmt=',f6.3,' Rnm=',f6.3,' Ion=',i2,1x,a10)
  140 format('Abs   Z=',i2,' Rmt=',f6.3,' Rnm=',f6.3, 1x,a10)

      do 150  iph = 1, nph
         ifr = ifrph(iph)
         if (ion(ifr) .ne. 0)  then
            nhead = nhead+1
            write(head(nhead),160)  iph, iz(ifr),  rmt(iph)*bohr,
     1           rnrm(iph)*bohr, ion(ifr)
         else
            nhead = nhead+1
            write(head(nhead),170)  iph, iz(ifr),  rmt(iph)*bohr,
     1           rnrm(iph)*bohr
         endif
  150 continue
  160 format('Pot',i2,' Z=',i2,' Rmt=',f6.3,' Rnm=',f6.3,' Ion=',i2)
  170 format('Pot',i2,' Z=',i2,' Rmt=',f6.3,' Rnm=',f6.3)
      if (abs(vi0) .gt. 1.0e-8 .or. abs(vr0) .gt. 1.0e-8)  then
         nhead = nhead+1
         write(head(nhead),180)  gamach*ryd, sout(ixc), vi0*ryd,
     1                           vr0*ryd
      else
         nhead = nhead+1
         write(head(nhead),190)  gamach*ryd, sout(ixc)
      endif
      nhead = nhead+1
  180 format('Gam_ch=',1pe9.3, 1x,a8, ' Vi=',1pe10.3, ' Vr=',1pe10.3)
  190 format('Gam_ch=',1pe9.3, 1x,a8)
  200 format('Mu=',1pe10.3, ' kf=',1pe9.3, ' Vint=',1pe10.3,
     x        ' Rs_int=',0pf6.3)
      write(head(nhead),200)  xmu*ryd, xf/bohr, vint*ryd, rs
      if (ixc .eq. 4)  then 
          nhead = nhead+1
          write(head(nhead),210)  rs0
  210     format ('Experimental DH-HL exch, rs0 = ', 1pe14.6)
      endif
      do 220  i = 1, nhead
         lhead(i) = istrln(head(i))
         heada(i) = head(i)
         lheada(i) = lhead(i)
  220 continue
      nheada = nhead

      return

      entry wthead (io)
c     Dump header to unit io, which must be open.  Add carraige control
c     to head array, which doesn't have it.

      do 310 i = 1, nheada
         ll = lheada(i)
         write(io,300)  heada(i)(1:ll)
  300    format (1x, a)
  310 continue
      end
