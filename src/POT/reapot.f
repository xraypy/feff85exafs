      subroutine reapot (mpot, rgrd, ntitle, title, ipr1, ispec,
     1           nohole, ihole, gamach, nph, iz, lmaxsc, xnatph,
     2           xion, iunf, ixc, jumprm, iafolp, folp, inters, totvol,
     3           rfms1, lfms1, nscmt, ca1, nmix, ecv, icoul,
     4           novr, iphovr, nnovr, rovr,
     5           nat, rat, iphat, iatph)
      implicit double precision (a-h, o-z)

      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'

cc    geom.dat
        integer  nat, iatph(0:nphx), iphat(natx)
        double precision  rat(3,natx)
cc    mod1.inp
        character*80 title(nheadx), head(nheadx)
        integer lhead(nheadx)
        integer mpot, nph, ntitle, ihole, ipr1, iafolp, ixc, ispec,
     1     iunf, nmix, nohole, jumprm, inters, nscmt, icoul, lfms1
        integer iz(0:nphx), lmaxsc(0:nphx)
        real rfms1
        double precision gamach, rgrd, ca1, ecv, totvol
        double precision  xnatph(0:nphx), folp(0:nphx), xion(0:nphx)
c       for OVERLAP option
        integer novr(0:nphx), iphovr(novrx,0:nphx), nnovr(novrx,0:nphx)
        double precision  rovr(novrx,0:nphx)

c     Local stuff
      character*512 slog
      character*32 s1, s2, s3

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
           if (nat .gt. natx)  then
              write(slog,55) ' nat, natx ', nat, natx
              call wlog(slog)
  55          format(a, 2i10)
              stop 'Bad input'
           endif
           nat = nat+1
           read(3,*,end=60)  idum, (rat(j,nat),j=1,3), iphat(nat), i1b
           if (iphat(nat).gt.nph) nph = iphat(nat)
           if ( iatph(iphat(nat)).eq.0) iatph(iphat(nat)) = nat
        goto 50
  60    continue
        nat = nat-1
      close(3)

c     read mod1.inp
      open (file='mod1.inp', unit=3, status='old',iostat=ios)
      call chopen (ios, 'mod1.inp', 'reapot')
        read (3,10) slog
        read (3,20) mpot, nph, ntitle, ihole, ipr1, iafolp, ixc, ispec
        read (3,10) slog
        read (3,20)  nmix, nohole, jumprm, inters, nscmt, icoul, lfms1,
     1     iunf
        do 110 ititle = 1, ntitle
  110   read (3,10) title(ititle)
        read (3,10) slog
        read (3,30)  gamach, rgrd, ca1, ecv, totvol, rfms1
        read (3,10) slog
  120   format ( 2i5, 4f13.5)
        do 130 ip = 0, nph
  130   read (3,120) iz(ip), lmaxsc(ip), xnatph(ip), xion(ip), folp(ip)
c       for OVERLAP option
        read (3,10) slog
        read (3,20) ( novr(iph), iph=0,nph)
        read (3,10) slog
  140   format ( 2i5, f13.5)
        do 150 iph = 0, nph
        do 150 iovr = 1, novr(iph)
  150   read (3,140) iphovr(iovr, iph), nnovr(iovr,iph), rovr(iovr,iph)
      close(3)

c     transform to code units (bohrs and hartrees - atomic unuts)
      rfms1 = rfms1 / bohr
      gamach = gamach / hart
      ecv   = ecv   / hart
      totvol = totvol / bohr**3
      do 210 iat = 1, nat
      do 210 i = 1,3
        rat(i,iat) = rat (i, iat) / bohr
  210 continue
      do 220 iph = 0, nph
      do 220 iovr = 1, novr(iph)
         rovr(iovr,iph) = rovr(iovr,iph) / bohr
  220 continue

c     add lines to the title
      if (mpot.eq.1) then
         ntitle = ntitle + 1
         if (nat.gt.1) then
           if (rfms1.lt.0) rfms1 = 0
           if (nscmt.gt.0) then
             write(s1, 230) nscmt, rfms1*bohr, lfms1
  230        format(' POT  SCF', i4, f8.4, i4)
           else
             write(s1, 235) 
  235        format(' POT  Non-SCF' )
           endif
         else
           write(s1, 240) 
  240      format(' POT  used OVERLAP geometry,')
         endif
         if (nohole.eq.0) then
           write(s2, 310) 
  310      format(', NO core-hole,')
         elseif (nohole.eq.2) then
           write(s2, 315) 
  315      format(', screened core-hole,')
         else
           write(s2, 320) 
  320      format(', core-hole,')
         endif
         if (iafolp.lt.0) then
           write(s3, 330) folp(0)
  330      format(' FOLP (folp(0)=', f6.3, ')' )
         else
           write(s3, 340) folp(0)
  340      format(' AFOLP (folp(0)=', f6.3, ')' )
         endif
c        concatenate 3 strings into 1
         title(ntitle) = ' '
         ilen = istrln(s1)
         istart = 1
         iend = ilen
         title(ntitle)(istart:iend) = s1(1:ilen)
         ilen = istrln(s2)
         istart = iend + 1
         iend = iend + ilen
         title(ntitle)(istart:iend) = s2(1:ilen)
         ilen = istrln(s3)
         istart = iend + 1
         iend = iend + ilen
         title(ntitle)(istart:iend) = s3(1:ilen)
      endif

      return
      end
