      subroutine reapot (mpot, rgrd, ntitle, title, ipr1, ispec,
     1           nohole, ihole, gamach, nph, iz, lmaxsc, xnatph,
     2           xion, iunf, ixc, jumprm, iafolp, folp, inters, totvol,
     3           rfms1, lfms1, nscmt, ca1, nmix, ecv, icoul,
     4           novr, iphovr, nnovr, rovr,
     5           nat, rat, iphat, iatph)

c$$$      use json_module
c$$$      implicit double precision (a-h, o-z)
c$$$      logical :: found
c$$$      character*7 vname
c$$$      type(json_file) :: json   !the JSON structure read from the file:
c$$$      double precision toss
c$$$      integer,dimension(:),allocatable :: intgs
c$$$      character*80,dimension(:),allocatable :: strings
c$$$      double precision,dimension(:),allocatable :: dbpcs

      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'

cc    geom.dat
        integer  nat, iatph(0:nphx), iphat(natx), ibounc(natx)
        double precision  rat(3,natx)
cc    mod1.inp
        character*80 title(nheadx)
c        character*80 head(nheadx)
c        integer lhead(nheadx)
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
c      character*512 slog
      character*32 s1, s2, s3

c$$$      intgs   = [INTEGER::]
c$$$      strings = [CHARACTER*80::]
c$$$      dbpcs   = [DOUBLE PRECISION::]


c     standard formats for string, integers and real numbers
c  10  format(a)
c  20  format (20i4)
c  30  format (6f13.5)


c--json--c     Read  geom.dat file
c--json--      open (file='geom.dat', unit=3, status='old',iostat=ios)
c--json--c       read header from geom.dat 
c--json--        nhead = nheadx
c--json--        call rdhead (3, nhead, head, lhead)
c--json--        nat = 0
c--json--        nph = 0
c--json--        do 40 iph = 0, nphx
c--json--  40    iatph(iph) = 0
c--json--  50    continue
c--json--           if (nat .gt. natx)  then
c--json--              write(slog,55) ' nat, natx ', nat, natx
c--json--              call wlog(slog)
c--json--  55          format(a, 2i10)
c--json--              stop 'Bad input'
c--json--           endif
c--json--           nat = nat+1
c--json--           read(3,*,end=60)  idum, (rat(j,nat),j=1,3), iphat(nat), i1b
c--json--           if (iphat(nat).gt.nph) nph = iphat(nat)
c--json--           if ( iatph(iphat(nat)).eq.0) iatph(iphat(nat)) = nat
c--json--        goto 50
c--json--  60    continue
c--json--        nat = nat-1
c--json--      close(3)

      call json_read_geom(nat, nph, iatph, rat, iphat, ibounc)

c     read mod1.inp
c--json--      open (file='mod1.inp', unit=3, status='old',iostat=ios)
c--json--      call chopen (ios, 'mod1.inp', 'reapot')
c--json--        read (3,10) slog
c--json--        read (3,20) mpot, nph, ntitle, ihole, ipr1, iafolp, ixc, ispec
c--json--        read (3,10) slog
c--json--        read (3,20)  nmix, nohole, jumprm, inters, nscmt, icoul, lfms1,
c--json--     1     iunf
c--json--        do 110 ititle = 1, ntitle
c--json--  110   read (3,10) title(ititle)
c--json--        read (3,10) slog
c--json--        read (3,30)  gamach, rgrd, ca1, ecv, totvol, rfms1
c--json--        read (3,10) slog
c--json--  120   format ( 2i5, 4f13.5)
c--json--        do 130 ip = 0, nph
c--json--  130   read (3,120) iz(ip), lmaxsc(ip), xnatph(ip), xion(ip), folp(ip)
c--json--c       for OVERLAP option
c--json--        read (3,10) slog
c--json--        read (3,20) ( novr(iph), iph=0,nph)
c--json--        read (3,10) slog
c--json--  140   format ( 2i5, f13.5)
c--json--        do 150 iph = 0, nph
c--json--        do 150 iovr = 1, novr(iph)
c--json--  150   read (3,140) iphovr(iovr, iph), nnovr(iovr,iph), rovr(iovr,iph)
c--json--      close(3)

      call json_read_pot(mpot, nph, ntitle, ihole, ipr1, iafolp,
     1       ixc, ispec, nmix, nohole, jumprm, inters, nscmt, icoul,
     2       lfms1, iunf, gamach, rgrd, ca1, ecv, totvol, rfms1,
     3       title, iz, lmaxsc, xnatph, xion, folp, novr,
     4       iphovr, nnovr, rovr )
      
c$$$      call json%load_file('pot.json')
c$$$      if (json_failed()) then   !if there was an error reading the file
c$$$         print *, "failed to read pot.json"
c$$$         stop
c$$$      else
c$$$         call json%get('mpot',   mpot, found)
c$$$                   if (.not. found) call bailout('mpot', 'pot.json')
c$$$         call json%get('nph',    nph, found)
c$$$                   if (.not. found) call bailout('nph', 'pot.json')
c$$$         call json%get('ntitle', ntitle, found)
c$$$                   if (.not. found) call bailout('ntitle', 'pot.json')
c$$$         call json%get('ihole',  ihole, found)
c$$$                   if (.not. found) call bailout('ihole', 'pot.json')
c$$$         call json%get('ipr1',   ipr1, found)
c$$$                   if (.not. found) call bailout('ipr1', 'pot.json')
c$$$         call json%get('iafolp', iafolp, found)
c$$$                   if (.not. found) call bailout('iafolp', 'pot.json')
c$$$         call json%get('ixc',    ixc, found)
c$$$                   if (.not. found) call bailout('ixc', 'pot.json')
c$$$         call json%get('ispec',  ispec, found)
c$$$                   if (.not. found) call bailout('ispec', 'pot.json')
c$$$
c$$$         call json%get('nmix',   nmix, found)
c$$$                   if (.not. found) call bailout('nmix', 'pot.json')
c$$$         call json%get('nohole', nohole, found)
c$$$                   if (.not. found) call bailout('nohole', 'pot.json')
c$$$         call json%get('jumprm', jumprm, found)
c$$$                   if (.not. found) call bailout('jumprm', 'pot.json')
c$$$         call json%get('inters', inters, found)
c$$$                   if (.not. found) call bailout('inters', 'pot.json')
c$$$         call json%get('nscmt',  nscmt, found)
c$$$                   if (.not. found) call bailout('nscmt', 'pot.json')
c$$$         call json%get('icoul',  icoul, found)
c$$$                   if (.not. found) call bailout('icoul', 'pot.json')
c$$$         call json%get('lfms1',  lfms1, found)
c$$$                   if (.not. found) call bailout('lfms1', 'pot.json')
c$$$         call json%get('iunf',   iunf, found)
c$$$                   if (.not. found) call bailout('iunf', 'pot.json')
c$$$
c$$$         call json%get('gamach', gamach, found)
c$$$                   if (.not. found) call bailout('gamach', 'pot.json')
c$$$         call json%get('rgrd',   rgrd, found)
c$$$                   if (.not. found) call bailout('rgrd', 'pot.json')
c$$$         call json%get('ca1',    ca1, found)
c$$$                   if (.not. found) call bailout('ca1', 'pot.json')
c$$$         call json%get('ecv',    ecv, found)
c$$$                   if (.not. found) call bailout('ecv', 'pot.json')
c$$$         call json%get('totvol', totvol, found)
c$$$                   if (.not. found) call bailout('totvol', 'pot.json')
c$$$         call json%get('rfms1',  toss, found)   
c$$$                   if (.not. found) call bailout('rfms1', 'pot.json')
c$$$         rfms1 = real(toss)
c$$$
c$$$         call json%get('titles', strings, found)
c$$$                   if (.not. found) call bailout('titles', 'pot.json')
c$$$         do 1000 itit = 1, nheadx
c$$$            title(itit) = strings(itit)            
c$$$ 1000    continue
c$$$         call json%get('iz', intgs, found)
c$$$                   if (.not. found) call bailout('iz', 'pot.json')
c$$$         do 1010 iph = 0, nphx
c$$$            iz(iph) = intgs(iph+1)            
c$$$ 1010    continue
c$$$         call json%get('lmaxsc', intgs, found)
c$$$                   if (.not. found) call bailout('lmaxsc', 'pot.json')
c$$$         do 1020 iph = 0, nphx
c$$$            lmaxsc(iph) = intgs(iph+1)            
c$$$ 1020    continue
c$$$         call json%get('xnatph', dbpcs, found)
c$$$                   if (.not. found) call bailout('xnatph', 'pot.json')
c$$$         do 1030 iph = 0, nphx
c$$$            xnatph(iph) = dbpcs(iph+1)            
c$$$ 1030    continue
c$$$         call json%get('xion', dbpcs, found)
c$$$                   if (.not. found) call bailout('xion', 'pot.json')
c$$$         do 1040 iph = 0, nphx
c$$$            xion(iph) = dbpcs(iph+1)            
c$$$ 1040    continue
c$$$         call json%get('folp', dbpcs, found)
c$$$                   if (.not. found) call bailout('folp', 'pot.json')
c$$$         do 1050 iph = 0, nphx
c$$$            folp(iph) = dbpcs(iph+1)            
c$$$ 1050    continue
c$$$         call json%get('novr', intgs, found)
c$$$                   if (.not. found) call bailout('novr', 'pot.json')
c$$$         do 1060 iph = 0, nphx
c$$$            novr(iph) = intgs(iph+1)            
c$$$ 1060    continue
c$$$
c$$$c        the following reconstructed all the overlap 2D arrays, see
c$$$c        RDINP/wrtjsn.f line 131 and following
c$$$         do 1200 iph = 0, nph
c$$$            write (vname, "(A6,I1)") "iphovr", iph
c$$$            call json%get(vname, intgs, found)
c$$$                      if (.not. found) call bailout(vname, 'pot.json')
c$$$            do 1220 iovr = 1, novr(iph)
c$$$               iphovr(iovr, iph) = intgs(iovr)
c$$$ 1220       continue
c$$$
c$$$            write (vname, "(A5,I1)") "nnovr", iph
c$$$            call json%get(vname, intgs, found)
c$$$                      if (.not. found) call bailout(vname, 'pot.json')
c$$$            do 1230 iovr = 1, novr(iph)
c$$$               nnovr(iovr, iph) = intgs(iovr)
c$$$ 1230       continue
c$$$
c$$$            write (vname, "(A4,I1)") "rovr", iph
c$$$            call json%get(vname, dbpcs, found)
c$$$                      if (.not. found) call bailout(vname, 'pot.json')
c$$$            do 1240 iovr = 1, novr(iph)
c$$$               rovr(iovr, iph) = dbpcs(iovr)
c$$$ 1240       continue
c$$$ 1200    continue
c$$$
c$$$         call json%destroy()
c$$$      end if

c     transform to code units (bohrs and hartrees - atomic unuts)
      rfms1 = rfms1 / real(bohr)
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
