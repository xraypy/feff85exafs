      subroutine ffsort(iabs, natt, ratx, iphatx,
     1       nabs, iphabs, rclabs, ipol, ispin, le2,
     2       elpty, angks, evec, xivec, spvec, ptz,
     3       iatph)


c***********************************************************************
c input/output:
c      ratx:   atomic positions from feff.inp on input
c              shifted and sorted by distance on output
c      iphatx: unique potentials from feff.inp on input
c              sorted by distance on output
c
c input:
c      nabs:  absorber index
c      natt:  number of atoms in cluster
c      nabs:   \
c      iphabs:  } stuff related to CFAVERAGE, ignored
c      rclabs  /
c      ipol:   flag for polarization calculation
c      ispin:  flag for spin-dep calculation
c      le2:    something related to MULTIPOLE, ignored
c      elpty:  eccentricity of elliptical light
c      angks:  angle between beam and spin 
c      evec:   polarization vector
c      xivec:  ellipticipty vecotr
c      spvec:  spin vector
c
c output
c      ptz:    polarization tensor
c      iatph:  list of example atoms for each unique potential
c***********************************************************************

c removed doptz boolean, which seems to be used only for eels, which is
c outside the scope of feff85exafs

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

      integer nat, natt, ipol, ispin, le2, iabs, nabs, iphabs
      integer iatph(0:nphx), iphat(natx), iphatx(natx), index(natx)
      double precision rclabs, elpty, angks
      double precision rat(3,natx), ratx(3,natx)
      double precision evec(3), xivec(3), spvec(3)
      complex*16 ptz(-1:1, -1:1)

c     Local stuff
      parameter (big = 1.0d5)
      character*512 slog

      external dist

c      call json_read_atoms(natt, ratx, iphatx)

c      call json_read_global(nabs, iphabs, rclabs, ipol, ispin, le2,
c     1                      elpty, angks, evec, xivec, spvec, ptz)


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
         call par_stop('libpotph')
      endif
      if (iphabs.eq.0 .and. icount.gt.1) then
         call wlog(' More than one absorbing atom (potential 0)')
         call wlog(' Only one absorbing atom allowed')
         call par_stop('libpotph')
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
            if (tmp.gt.0.1d0 .and. tmp.le.rclabs) then
               nat = nat + 1
               if (nat.gt.natx) then
                 write (slog, 307) nat, natx
  307            format (' Number of atoms', i6, 'exceeds max allowed',
     1           ' for the pathfinder =', i6)
                 call wlog (' Use or reduce rclabs in CFAVERAGE card')
                 call wlog (' Or increase parameter natx and recompile')
                 call par_stop('libpotph')
               endif
               rat(1,nat) = ratx(1,iat)-ratx(1,iatabs)
               rat(2,nat) = ratx(2,iat)-ratx(2,iatabs)
               rat(3,nat) = ratx(3,iat)-ratx(3,iatabs)
               iphat(nat) = iphatx(iat)
               index(nat) = iat
            endif
         endif
 309  continue
c     squelch a compiler warning
      i=1
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
      call mkptz(ipol, elpty, evec, xivec, ispin, spvec, nat, rat,
     1       angks, le2, ptz)

c     rewrite global.json for initial iteration to update 'ptz'
c      call json_global(nabs)


c     Find model atoms for unique pots that have them
c     Use atom closest to absorber for model
      do 316  iph = 1, nphx
 316  iatph(iph) = 0
c     By construction absorbing atom is first in the list
      iatph(0) = 1
c      nph = 0
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
c         if (iatph(iph).gt.0) nph = iph
  330 continue
c     if iatph > 0, a model atom has been found.

c     Check if 2 atoms are closer together than 1.75 bohr (~.93 Ang)
      ratmin = 1.0d20
      do 480  iat = 1, nat
         do 470  jat = iat+1, nat
            rtmp = dist(rat(1,iat),rat(1,jat))
            if (rtmp .lt. ratmin)  ratmin = rtmp
            if (rtmp .lt. 1.75d0 * bohr)  then
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

c      call json_geom(iatph,rat,iphat)


c     transfer rat back to ratx and iphat back to iphatx
      do 1000 j=1,nat
         ratx(1,j) = rat(1,j)
         ratx(2,j) = rat(2,j)
         ratx(3,j) = rat(3,j)
         iphatx(j) = iphat(j)
 1000 continue

c     Atoms for the pathfinder
      if (iatabs .le. 0)  then
         call wlog(' Absorbing atom coords not specified.')
         call wlog(' Cannot find multiple scattering paths.')
         call par_stop('libpotph')
      endif

c 400 call par_barrier

      return
      end
