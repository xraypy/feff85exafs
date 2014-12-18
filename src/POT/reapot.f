      subroutine reapot (mpot, rgrd, ntitle, title, ipr1, ispec,
     1           nohole, ihole, gamach, nph, iz, lmaxsc, xnatph,
     2           xion, iunf, ixc, jumprm, iafolp, folp, inters, totvol,
     3           rfms1, lfms1, nscmt, ca1, nmix, ecv, icoul,
     4           novr, iphovr, nnovr, rovr,
     5           nat, rat, iphat, iatph)

      use json_module
      implicit double precision (a-h, o-z)
      logical :: found
      character*7 vname
      type(json_file) :: json   !the JSON structure read from the file:
      double precision toss
      integer,dimension(:),allocatable :: intgs
      character*80,dimension(:),allocatable :: strings
      double precision,dimension(:),allocatable :: dbpcs

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

      intgs   = [INTEGER::]
      strings = [CHARACTER*80::]
      dbpcs   = [DOUBLE PRECISION::]


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

      call json%load_file('pot.json')
      if (json_failed()) then   !if there was an error reading the file
         print *, "failed to read pot.json"
         stop
      else
         call json%get('mpot',   mpot, found)
                   if (.not. found) call bailout('mpot', 'pot.json')
         call json%get('nph',    nph, found)
                   if (.not. found) call bailout('nph', 'pot.json')
         call json%get('ntitle', ntitle, found)
                   if (.not. found) call bailout('ntitle', 'pot.json')
         call json%get('ihole',  ihole, found)
                   if (.not. found) call bailout('ihole', 'pot.json')
         call json%get('ipr1',   ipr1, found)
                   if (.not. found) call bailout('ipr1', 'pot.json')
         call json%get('iafolp', iafolp, found)
                   if (.not. found) call bailout('iafolp', 'pot.json')
         call json%get('ixc',    ixc, found)
                   if (.not. found) call bailout('ixc', 'pot.json')
         call json%get('ispec',  ispec, found)
                   if (.not. found) call bailout('ispec', 'pot.json')

         call json%get('nmix',   nmix, found)
                   if (.not. found) call bailout('nmix', 'pot.json')
         call json%get('nohole', nohole, found)
                   if (.not. found) call bailout('nohole', 'pot.json')
         call json%get('jumprm', jumprm, found)
                   if (.not. found) call bailout('jumprm', 'pot.json')
         call json%get('inters', inters, found)
                   if (.not. found) call bailout('inters', 'pot.json')
         call json%get('nscmt',  nscmt, found)
                   if (.not. found) call bailout('nscmt', 'pot.json')
         call json%get('icoul',  icoul, found)
                   if (.not. found) call bailout('icoul', 'pot.json')
         call json%get('lfms1',  lfms1, found)
                   if (.not. found) call bailout('lfms1', 'pot.json')
         call json%get('iunf',   iunf, found)
                   if (.not. found) call bailout('iunf', 'pot.json')

         call json%get('gamach', gamach, found)
                   if (.not. found) call bailout('gamach', 'pot.json')
         call json%get('rgrd',   rgrd, found)
                   if (.not. found) call bailout('rgrd', 'pot.json')
         call json%get('ca1',    ca1, found)
                   if (.not. found) call bailout('ca1', 'pot.json')
         call json%get('ecv',    ecv, found)
                   if (.not. found) call bailout('ecv', 'pot.json')
         call json%get('totvol', totvol, found)
                   if (.not. found) call bailout('totvol', 'pot.json')
         call json%get('rfms1',  toss, found)   
                   if (.not. found) call bailout('rfms1', 'pot.json')
         rfms1 = real(toss)

         call json%get('titles', strings, found)
                   if (.not. found) call bailout('titles', 'pot.json')
         do 1000 itit = 1, nheadx
            title(itit) = strings(itit)            
 1000    continue
         call json%get('iz', intgs, found)
                   if (.not. found) call bailout('iz', 'pot.json')
         do 1010 iph = 0, nphx
            iz(iph) = intgs(iph+1)            
 1010    continue
         call json%get('lmaxsc', intgs, found)
                   if (.not. found) call bailout('lmaxsc', 'pot.json')
         do 1020 iph = 0, nphx
            lmaxsc(iph) = intgs(iph+1)            
 1020    continue
         call json%get('xnatph', dbpcs, found)
                   if (.not. found) call bailout('xnatph', 'pot.json')
         do 1030 iph = 0, nphx
            xnatph(iph) = dbpcs(iph+1)            
 1030    continue
         call json%get('xion', dbpcs, found)
                   if (.not. found) call bailout('xion', 'pot.json')
         do 1040 iph = 0, nphx
            xion(iph) = dbpcs(iph+1)            
 1040    continue
         call json%get('folp', dbpcs, found)
                   if (.not. found) call bailout('folp', 'pot.json')
         do 1050 iph = 0, nphx
            folp(iph) = dbpcs(iph+1)            
 1050    continue
         call json%get('novr', intgs, found)
                   if (.not. found) call bailout('novr', 'pot.json')
         do 1060 iph = 0, nphx
            novr(iph) = intgs(iph+1)            
 1060    continue

c        the following reconstructed all the overlap 2D arrays, see
c        RDINP/wrtjsn.f line 131 and following
         do 1200 iph = 0, nph
            write (vname, "(A6,I1)") "iphovr", iph
            call json%get(vname, intgs, found)
                      if (.not. found) call bailout(vname, 'pot.json')
            do 1220 iovr = 1, novr(iph)
               iphovr(iovr, iph) = intgs(iovr)
 1220       continue

            write (vname, "(A5,I1)") "nnovr", iph
            call json%get(vname, intgs, found)
                      if (.not. found) call bailout(vname, 'pot.json')
            do 1230 iovr = 1, novr(iph)
               nnovr(iovr, iph) = intgs(iovr)
 1230       continue

            write (vname, "(A4,I1)") "rovr", iph
            call json%get(vname, dbpcs, found)
                      if (.not. found) call bailout(vname, 'pot.json')
            do 1240 iovr = 1, novr(iph)
               rovr(iovr, iph) = dbpcs(iovr)
 1240       continue
 1200    continue

         call json%destroy()
      end if

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
