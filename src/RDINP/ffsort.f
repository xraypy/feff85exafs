      subroutine ffsort (iabs,doptz)

c KJ 1-06 : I added second input argument doptz      
      implicit double precision (a-h, o-z)

c     finds iabs-th atom of 'iphabs' type in file atoms.dat and writes
c     a smaller list of all atoms within 'rclabs' of that particular
c     absorber into 'geom.dat' file.
c      first coded by a.l.ankudinov, 1998 for CFAVERAGE card
c      modified by a.l.ankudinov, march 2001 for new i/o structure

      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'
      include '../HEADERS/vers.h'
      include '../HEADERS/parallel.h'
      include '../RDINP/allinp.h'

cc    INPUT
cc    atoms.dat
c--allinp--        integer  natt
c--allinp--        integer iphatx(nattx)
c--allinp--        double precision  ratx(3,nattx)
        logical doptz  !KJ 1-06 : call mkptz or not?    
cc    global.dat
c       configuration average
c--allinp--        integer nabs, iphabs
c       global polarization data
c--allinp--        integer  ipol, ispin, le2
c--allinp--        double precision evec(3), xivec(3), spvec(3), elpty,angks,rclabs
c--allinp--        complex*16 ptz(-1:1, -1:1)
cc    OUTPUT: geom.dat
        integer  nat
        integer iatph(0:nphx), iphat(natx), index(natx)
        double precision  rat(3,natx)

c     Local stuff
      parameter (big = 1.0e5)
      character*512 slog

      external dist

c     if (worker) go to 400

c     standard formats for string, integers and real numbers
c  10  format(a)
c  20  format (20i4)
c  30  format (6f13.5)

c--json--cc    read atoms.dat file
c--json--      open (file='atoms.dat', unit=3, status='old',iostat=ios)
c--json--        read(3, 35) slog, natt
c--json--  35    format (a8, i7)
c--json--        read  (3, 10) slog
c--json--        do 40  iat = 1, natt
c--json--          read (3,36) ratx(1,iat), ratx(2,iat), ratx(3,iat), iphatx(iat)
c--json--  36      format( 3f13.5, i4)
c--json--  40    continue
c--json--      close(3)

c  45    format ( 2i8, f13.5)
c  50    format ( 3i5, 2f12.4)
c--json--c     read global.inp
c--json--c     CFAVERAGE iphabs nabs rclabs
c--json--        open (file='global.dat', unit=3, status='old',iostat=ios)
c--json--        call chopen (ios, 'global.inp', 'ffsort')
c--json--        read (3, 10) slog
c--json--        read (3, 45) nabs, iphabs, rclabs
c--json--c       global polarization data
c--json--        read  (3,10) slog
c--json--        read  (3, 50)  ipol, ispin, le2, elpty, angks
c--json--        read  (3, 10) slog
c--json--        do 60 i = 1,3
c--json--          read  (3,30) evec(i), xivec(i), spvec(i)
c--json--  60    continue
c--json--        read  (3, 10) slog
c--json--        do 70 i = -1, 1
c--json--          read (3,30) a1, b1, a2, b2, a3, b3
c--json--          ptz(-1,i)= dcmplx(a1,b1) 
c--json--          ptz(0,i) = dcmplx(a2,b2) 
c--json--          ptz(1,i) = dcmplx(a3,b3) 
c--json--  70    continue
c--json--      close(3)

      call json_read_atoms(natt, ratx, iphatx)

      call json_read_global(nabs, iphabs, rclabs, ipol, ispin, le2,
     1                      elpty, angks, evec, xivec, spvec, ptz)


c     Find the first absorber (iphabs type) in a long list (iabs.le.0),
c     or find iabs-th atom in the list of type iphabs (iabs.gt.0)
      iatabs = 0
      icount = 0
      ifound = 0
      do 305 iat = 1, natt
         if (iphatx(iat) .eq. 0) iphatx(iat) = iphabs
         if (iphatx(iat) .eq. iphabs) icount = icount +1
         if (ifound.eq.0 .and. icount.gt.0 .and. (icount.eq.iabs .or.
     1                          (iabs.le.0 .and. icount.eq.1))) then
            iatabs = iat
            ifound =1
         endif
  305 continue

c     Make several sanity checks
      if (iatabs.eq.0 .and. natt.gt.1) then
         call wlog(' No absorbing atom (unique pot 0 or iphabs in'//
     1             ' CFAVERAGE  card) was defined.')
         call par_stop('RDINP')
      endif
      if (iphabs.eq.0 .and. icount.gt.1) then
         call wlog(' More than one absorbing atom (potential 0)')
         call wlog(' Only one absorbing atom allowed')
         call par_stop('RDINP')
      endif
      if ((icount.gt.0 .and. icount.lt.nabs) .or. nabs.le.0) then
         nabs = icount
         call wlog(' Averaging over ALL atoms of iphabs type')
      endif

c     Make absorbing atom first in the short list
      if (iatabs .ne. 0) then
         rat(1,1) = 0
         rat(2,1) = 0
         rat(3,1) = 0
         iphat(1) = 0
         index(1) = iatabs
      endif
          
c     make a smaller list of atoms from a big one
      nat = 1
      do 309 iat = 1,natt
         if (iat.ne.iatabs) then
            tmp = dist (ratx(1,iat), ratx(1,iatabs))
            if (tmp.gt.0.1 .and. tmp.le.rclabs) then
               nat = nat + 1
               if (nat.gt.natx) then
                 write (slog, 307) nat, natx
  307            format (' Number of atoms', i6, 'exceeds max allowed',
     1           ' for the pathfinder =', i6)
                 call wlog (' Use or reduce rclabs in CFAVERAGE card')
                 call wlog (' Or increase parameter natx and recompile')
                 call par_stop('RDINP')
               endif
               rat(1,nat) = ratx(1,iat)-ratx(1,iatabs)
               rat(2,nat) = ratx(2,iat)-ratx(2,iatabs)
               rat(3,nat) = ratx(3,iat)-ratx(3,iatabs)
               iphat(nat) = iphatx(iat)
               index(nat) = iat
            endif
         endif
 309  continue
c     sort atoms by distance
      do 315 iat = 1,nat-1
        r2min = rat(1,iat)**2 + rat(2,iat)**2 + rat(3,iat)**2
        imin = iat
        do 310 i = iat+1,nat
          r2 = rat(1,i)**2 + rat(2,i)**2 + rat(3,i)**2
          if (r2.lt.r2min) then
            r2min = r2
            imin = i
          endif
 310    continue
        if (imin.ne.iat) then
c         permute coordinates for atoms iat and imin
          do 311 i = 1,3
            r2 = rat(i,iat)
            rat(i,iat) = rat(i,imin)
            rat(i,imin) = r2
 311      continue
          i = iphat(iat)
          iphat(iat) = iphat(imin)
          iphat(imin) = i
          i = index(iat)
          index(iat) = index(imin)
          index(imin) = i
        endif
 315  enddo

c     rotate xyz frame for the most convinience and make
c     polarization tensor
c     make polarization tensor when z-axis is along k-vector 
      if (doptz)  !KJ I added this if-statement 1-06
     1  call mkptz( ipol, elpty, evec, xivec, ispin, spvec, nat, rat,
     2           angks, le2, ptz)
c     rewrite global.json for initial iteration to update 'ptz'
      call json_global(nabs)

c--json--c     rewrite global.inp for initial iteration to update 'ptz'
c--json--      open (file='global.dat', unit=3, status='unknown',iostat=ios)
c--json--c       configuration average data
c--json--        write (3, 10) ' nabs, iphabs - CFAVERAGE data'
c--json--        write (3, 45) nabs, iphabs, rclabs
c--json--c       global polarization data
c--json--        write (3,10) ' ipol, ispin, le2, elpty, angks'
c--json--        write (3, 50)  ipol, ispin, le2, elpty, angks
c--json--        write (3, 10) 'evec         xivec        spvec'
c--json--        do 360 i = 1,3
c--json--          write (3,30) evec(i), xivec(i), spvec(i)
c--json-- 360    continue
c--json--        write (3, 10) ' polarization tensor '
c--json--        do 370 i = -1, 1
c--json--          write(3,30) dble(ptz(-1,i)), dimag(ptz(-1,i)), dble(ptz(0,i)),
c--json--     1                dimag(ptz(0,i)),  dble(ptz(1,i)), dimag(ptz(1,i))
c--json-- 370    continue
c--json--      close(3)

c     Find model atoms for unique pots that have them
c     Use atom closest to absorber for model
      do 316  iph = 1, nphx
 316  iatph(iph) = 0
c     By construction absorbing atom is first in the list
      iatph(0) = 1
      nph = 0
      do 330  iph = 1, nphx
         rabs = big
         do 320  iat = 2, nat
            if (iph .eq. iphat(iat))  then
               tmp = dist (rat(1,iat), rat(1,1))
               if (tmp .lt. rabs)  then
c                 this is the closest so far
                  rabs = tmp
                  iatph(iph) = iat
               endif
            endif
  320    continue
         if (iatph(iph).gt.0) nph = iph
  330 continue
c     if iatph > 0, a model atom has been found.

c     Check if 2 atoms are closer together than 1.75 bohr (~.93 Ang)
      ratmin = 1.0e20
      do 480  iat = 1, nat
         do 470  jat = iat+1, nat
            rtmp = dist(rat(1,iat),rat(1,jat))
            if (rtmp .lt. ratmin)  ratmin = rtmp
            if (rtmp .lt. 1.75 * bohr)  then
               call wlog(' WARNING:  TWO ATOMS VERY CLOSE TOGETHER.' //
     1                   '  CHECK INPUT.')
               iatx = index(iat)
               jatx = index(jat)
               write(slog,'(a,2i8)') ' atoms ', iatx, jatx
               call wlog(slog)
               write(slog,'(i5,1p,3e13.5)') iatx, (ratx(i,iatx),i=1,3)
               call wlog(slog)
               write(slog,'(i5,1p,3e13.5)') jatx, (ratx(i,jatx),i=1,3)
               call wlog(slog)
               call wlog(' Run continues in case you really meant it.')
            endif
  470    continue
  480 continue

c--json--c     Write output geom.dat
c--json--      open (file='geom.dat', unit=3, status='unknown',iostat=ios)
c--json--        write (3,535) nat, nph
c--json--  535   format ('nat, nph = ', 2i5)
c--json--        write (3,516) (iatph(iph), iph=0,nph)
c--json--  516   format(16i5)
c--json--        write (3, 10) ' iat     x       y        z       iph  '
c--json--        write (3, 526)
c--json--  526   format (1x, 71('-'))
c--json--        ibounc = 1
c--json--        do 540  i = 1, nat
c--json--          write(3,536) i, rat(1,i), rat(2,i), rat(3,i), iphat(i), ibounc
c--json--  536     format(i4, 3f13.5, 2i4)
c--json--  540   continue
c--json--      close(3)
      call json_geom(iatph,rat,iphat)

c     Atoms for the pathfinder
      if (iatabs .le. 0)  then
         call wlog(' Absorbing atom coords not specified.')
         call wlog(' Cannot find multiple scattering paths.')
         call par_stop('RDINP')
      endif

c 400 call par_barrier

      return
      end
