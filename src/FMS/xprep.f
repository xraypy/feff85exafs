      subroutine xprep(iph0, idwopt, nat, inclus, npot,
     $            iphat, rmax, rat,
     $            izx, rnrmav, temper, thetad, sig2,
     $            minv, rdirec )

      implicit real (a-h,o-z)
      implicit integer (i-n)
c--------------------------------------------------------------------
c This subroutine prepares geometrical information for a subsequent
c call to the fms full multiple scattering package.  Information is
c passed to fms via common blocks (which live in header files during
c development).
c
c These header files are required:  dim.h xparam.h  xstruc.h
c dim.h and xparam.h must be included in the calling routine
c This routine calls wlog, so must be compiled with it
c
c This is the main file of xpreppack for use with the full multiple
c scattering package (fmspack).  The calling protocol for xpreppack
c and fmspack is;
c
c      include '../HEADERS/dim.h'
c      include 'xparam.h'
c      ...
c      call xprep(iph0, idwopt, nat, inclus, npot, iphat, rmax, rat,
c $        izx, rnrmav, temper, thetad)
c      energy loop {
c         ...
c         call fms(nsp, inclus, npot, ck, lipotx, xphase, ik, iverb, gg)
c         ... }
c
c xpreppack contains the following routines:
c   xprep:  main routine of xpreppack
c   getang: determine angles between the z axis and all pairs of atoms
c   rotxan: get all rotation matrix elements for the cluster
c   rotint: initialize arrays used in the construction of rotation
c           matrices
c   sortat: organize atoms and potentials lists for computational and
c             organizational efficiency
c   atheap: heap sort extended cluster by distance from central atom
c   xanlm:  get all legendre normalization factor
c   factst: part of legendre factor computation
c
c xpreppack currently supports use of the Debye-Waller factors for
c   estimating the effect of thermal disorder on the xanes spectrum.
c   It does this by filling a matrix with the pairwise mean square
c   displacement between atoms.  Other forms of this calculation may
c   be included in the future.  Note that it is strictly impossible
c   to correctly model disorder in the MS scattering contribution to
c   the spectrum when using the FMS technique.
c--------------------------------------------------------------------
c  input:
c     iph0:   potential index for DOS calculations (added by ala to
c             handle other-than-the-central atom) (iph0=0 for the
c             absorbing atom)
c     nat:    number of atoms in extended cluster
c     npot:   number of unique potentials in cluster
c     iphat:  (natxx) potential index for each atom in extended
c             cluster
c     rmax:   radial size of cluster
c     rat:    (3, natxx) coordiantes of each atom in extended cluster
c             as read from geometry file
c  input for correlated debye model:
c     izx:    (natxx) Z number of each atom in the cluster
c     rnrmav: average norman radius in cluster (from phase.pad)
c     temper: sample temperature
c     thetad: Debye temperature
c
c  output:
c     inclus: number of atoms in cluster (inclus <= nat)
c
c  output (all via commmon blocks in xstruc.h):
c     xphi:  matrix of angles between the z axis and pairs of atoms
c     xrat:  xyz coordinates of the atoms in the cluster, the first
c            npot+1 entries are examples of each unique potential
c     iphx:  potential indeces of each atom in the cluster, ordered
c            like xrat
c     drix:  huge matrix containing all rotation matrix elements
c            needed for computation of free elctron propagators
c     xnlm:  matrix of legendre polynomial normalization factors
c--------------------------------------------------------------------

      include '../HEADERS/dim.h'
      include '../HEADERS/const.h'
      include '../HEADERS/parallel.h'
      include 'xparam.h'
c====================================================================
c  This header file contains the structural information about the
c  cluster to be used for the full multiple scattering calculation

      common /xstruc/ xphi(nclusx,nclusx), xrat(3,nclusx),
     $            iphx(nclusx)
      save /xstruc/

c  xphi:  matrix of angles between the z axis and pairs of atoms
c  xrat:  xyz coordinates of the atoms in the cluster, the first
c         npot+1 entries are examples of each unique potential
c  iphx:  potential indeces of each atom in the cluster, ordered like
c         xrat
c********************************************************************
c**** rotation matrices for entire cluster
      complex drix
      common /rotx/ drix(-lx:lx, -lx:lx, 0:lx, 0:1, nclusx, nclusx)
      save /rotx/
c********************************************************************
c common blocks for saving rotation matrices between xanes and rotxan
       integer    jsavx, jsav, jbmagk
       parameter (jsavx = 150, roteps = 1.e-12,jbmagk=-9999)
       dimension drisav(-lx:lx,-lx:lx,0:lx,jsavx), betsav(jsavx)
       integer   ldsav(jsavx), mdsav(jsavx)
       common /rotsav/  drisav, betsav, ldsav, mdsav, jsav
       save  /rotsav/
c********************************************************************
c**** legendre polynomial normalization constants
c
      common /lnlm/ xnlm(0:lx,0:lx)
      save   /lnlm/

c********************************************************************
c**** legendre polynomial normalization constants
c
      common /xdwf/ sigsqr(nclusx,nclusx)
      save   /xdwf/
c********************************************************************
c**** save Clebsch-Gordon coefficients: <LS|J>
      dimension t3jp(0:lx, -lx:lx, 2), t3jm(0:lx, -lx:lx, 2)
      common /t3j/ t3jp, t3jm
      save   /t3j/
c********************************************************************
      integer   iphat(natxx), iphat2(natxx), izx(0:nphasx), izpair(0:2)
      dimension rat(3,natxx), rat2(3,natxx)
      double precision ra(natxx)
      character*78 line
c     sigms is written in double precision.  these are the variables
c     that it uses
      double precision dtemp, dthet, drs, dsigsq, pair(3,0:2)
      double precision sig2mx, sig2x(0:nphx,0:nphx)
c     iwarn - needed to wrtite waqrning just one time
      integer iwarn
      save iwarn
      data  iwarn /0/
      double precision cwig3j
      external cwig3j

c     Josh - logical flag for reading in dw factors
      logical readdw
c     Josh END

c  initialize geometrical arrays
      do 30 i=1,nclusx
        do 10 j=1,nclusx
          xphi(j,i) = zero
 10     continue
        do 20 j=1,3
          xrat(j,i) = zero
 20     continue
        iphx(i) = 0
 30   continue
      inclus = 0

c --- find the central atom, ipot=iph0 (iph0=0 for the absorbing atom)
      icen = 0
      do 40 i=1,nat
        iphat2(i) = iphat(i)
        if (iphat(i).eq.iph0) then
            if (icen.eq.0) then
                icen = i
            elseif (iph0.eq.0) then
                call wlog('* * * ERROR!  More than one atom '//
     $                      'in the extended cluster have ipot=0')
                call wlog('      You may only have one central atom.')
                call wlog('      Stopping in xprep.')
                call par_stop('XPREP-1')
            endif
        endif
 40   continue
c --- make sure central atom is at (0,0,0)
      do 45 i=1,nat
        rat2(1,i) = rat(1,i)-rat(1,icen)
        rat2(2,i) = rat(2,i)-rat(2,icen)
        rat2(3,i) = rat(3,i)-rat(3,icen)
 45   continue

c --- sort the atoms from extended cluster by distance from central
c     atom.
      call atheap(nat, rat2, iphat2, ra)

c --- define cluster from extended cluster by as those closer than
c     rmax to central atom
      inclus=0
      rmax2 = rmax**2
      do 50 i=1,nat
        rr = (rat2(1,i)**2 + rat2(2,i)**2 + rat2(3,i)**2)
        if (rr .gt. rmax2) then
            inclus = i-1
            goto 60
        endif
 50   continue
 60   continue
      if (inclus.eq.0) inclus=nat

c --- sanity check size of cluster
      if (inclus.gt.nclusx) then
        if (iwarn.eq.0) then
          call wlog('* * * WARNING preparing cluster for '//
     $                'FMS calculation.')
          write(line,400) inclus
 400      format('      You specified a cluster of ', i3,
     $                ' atoms for the FMS calculation.')
          call wlog(line)
          write(line,410)nclusx
          call wlog(line)
 410      format('      This exceeds the hard wired limit of ', i3,
     $                ' atoms.')
          write(line,420)nclusx
          call wlog(line)
 420      format('      The cluster size was reset to ', i3,
     $                ' and the calculation will continue.')
          iwarn = 1
        endif
        inclus = nclusx
      endif

c --- make the first few entries in xrat represent each of the
c     unique potentials, sorting around the chosen center
c     (iph0=0 for the absorbing atom)
c     call sortat(iph0, inclus, npot, iphat2, iphx, rat2, xrat)
      do 430 iat = 1, inclus
          iphx(iat) = iphat2(iat)
          xrat(1,iat) = real (rat2(1,iat))
          xrat(2,iat) = real (rat2(2,iat))
          xrat(3,iat) = real (rat2(3,iat))
 430  continue

c --- Calculate and store rotation matrix elements and phi angles
c     the k loop calculates the forward then the backward rotation
c     for an atom pair (ij). k = 0-->forward, 1-->backward
      call rotint
      lplus1 = lx+1
      mplus1 = lx+1
      do 150  i=1,inclus
        do 140 j=1,inclus
          rr = (xrat(1,i)-xrat(1,j))**2 + (xrat(2,i)-xrat(2,j))**2 +
     1         (xrat(3,i)-xrat(3,j))**2
c         if (rr.gt. rdirec**2) goto 140

          call getang(nclusx, xrat, i, j, xbeta, xphi(i,j))
          if (i.eq.j) goto 140
          do 130 k=0,1
            if (k.eq.1) xbeta = (-1) * xbeta
            call rotxan(lplus1, mplus1, xbeta, i, j, k)
 130      continue
 140    continue
 150  continue

c --- calculate spherical harmonic normalization factors
      call xanlm(lplus1,mplus1)

c --- calculate array of correlated debye waller factors
c     initialize
      do 200 iat2=1,nclusx
      do 200 iat1=1,nclusx
        sigsqr(iat1,iat2) = zero
 200  continue

c            Josh - open dwfac.dat for reading precalculated dw factors
c                   and set logical flag readdw
c             open(unit=47,file='dwfac.dat',status='unknown',iostat=ios)
c            Temporary - write output file sigsqr.dat
c             open(unit=48,file='sigsqr.dat',status='replace')
c             if(ios.le.0) readdw = .true.
c            Josh END

c     open files for sigrm and sigem
      if (idwopt.eq.1.and.master) then
         iem = 111
         open(unit=iem,file='s2_em.dat',status='unknown', iostat=ios)
         call chopen (ios, 's2_em.dat', 'sigem')
      endif
      if (idwopt.ge.1) then
         if(master) then
           irm1 =113
           open(unit=irm1,file='s2_rm2.dat',status='unknown',
     .          iostat=ios)
           call chopen (ios, 's2_rm2.dat', 'sigrm')
           irm2 = 112
           open(unit=irm2,file='s2_rm1.dat',status='unknown',
     .          iostat=ios)
           call chopen (ios, 's2_rm1.dat', 'sigrm')
         endif

c        initialize statistics for max DW factors and set rmaxem
         sig2mx=0
         do 446 iph1=0,nphx
         do 446 iph2=0,nphx
  446    sig2x(iph1,iph2) = 0
         rmaxem = 20.0
         do 447 i = 2,inclus
           rr = 0
           do 448 j=1,3
  448      rr = rr + (xrat(j,i)-xrat(j,1))**2
           rr = sqrt(rr)
           if (rr.lt.rmaxem) rmaxem = rr
  447    continue
         rmaxem = max(2.2*rmaxem, 5.0/bohr)
      endif

      npair = 0
      do 250 iat1=1,inclus-1
        do 240 iat2=iat1+1, inclus
          rr = (xrat(1,iat1)-xrat(1,iat2))**2 +
     1    (xrat(2,iat1)-xrat(2,iat2))**2 +(xrat(3,iat1)-xrat(3,iat2))**2
c         if (rr.gt. rdirec**2) goto 240

          if (idwopt.ge.0) then
             do 230 ipair=1,3
               pair(ipair,0) = dble(xrat(ipair, iat1)*bohr)
               pair(ipair,1) = dble(xrat(ipair, iat2)*bohr)
               pair(ipair,2) = dble(pair(ipair,0))
 230         continue
             izpair(0) = izx(iphx(iat1))
             izpair(1) = izx(iphx(iat2))
             izpair(2) = izpair(0)
             dtemp = dble(temper)
             dthet = dble(thetad)
             drs   = dble(rnrmav)
             ipath0=0
             if (idwopt.eq.0) then
c               use CD model
                call sigms(dtemp,dthet,drs,2,2,pair,izpair,dsigsq)
             elseif (idwopt.eq.1) then
                xr12 = (xrat(1,iat1) - xrat(1,iat2))**2
                yr12 = (xrat(2,iat1) - xrat(2,iat2))**2
                zr12 = (xrat(3,iat1) - xrat(3,iat2))**2
                rr12 = sqrt(xr12 +yr12 +zr12)
                if (rr12.le.rmaxem) then
c                  use EM method
                   npair = npair + 1
                   if (mod(npair,100).eq.0) then
                      write (line, 337) npair
  337                 format('    Doing DW factors via EM method for',
     1                ' the pair number ', i5)
                      call wlog(line)
                   endif
                   call sigem
     1             (sig2mx,sig2x,iem,dtemp,ipath0,2,pair,dsigsq)
                else
c                  use RM method
                   call sigrm
     1             (sig2mx,sig2x,irm1,irm2,dtemp,ipath0,2,pair,dsigsq)
                endif
             else
c               use RM
                call sigrm
     1          (sig2mx,sig2x,irm1,irm2,dtemp,ipath0,2,pair,dsigsq)
             endif
             sigsqr(iat1,iat2) = real(dsigsq)
c            Josh - Temporary. write to sigsqr.dat
c             write(48,*) iat1, iat2, sigsqr(iat1,iat2)
          endif
          sigsqr(iat1,iat2) = sigsqr(iat1,iat2) + sig2
          sigsqr(iat2,iat1) = sigsqr(iat1,iat2)
 240    continue
 250  continue

c     Josh - if file dwfac.dat exists, read lines
c            iat1   iat2   dwfac
c            from it and set sigsqr(iat1,iat2) = dwfac
c      if(readdw) then
c         do i = 1, 1000
c            call RdCmt(47,'#*cC')
c            read(47,*, end=255) iat1, iat2, sigsqr(iat1,iat2)
c         end do
c 255     continue
c      end if
c      close(47)
c     Josh END

c     close output for sigem sigrm
      if (master) then
        if (idwopt.eq.1) close (unit=iem)
        if (idwopt.ge.1) close (unit=irm1)
        if (idwopt.ge.1) close (unit=irm2)
      endif

c     Calculate Clebsch-Gordon coefficients <LS|J>
      do 70  l1 = 0, lx
      do 70  mm = -l1, l1
      do 70  isp1 = 1, 2
        j1 = 2 * l1
        j2 = 1
        j3p = j1 + 1
        j3m = j1 - 1
        m1 = 2*mm
        m2 = 2*isp1 - 3
c  j = l+1/2
        t3jp( l1, mm, isp1) = sqrt( j3p + 1.0e0 ) *
     1                 real( cwig3j( j1, j2, j3p, m1, m2, 2) )
        if (mod( (j2-j1-m1-m2)/2 , 2) .ne.0)
     1          t3jp( l1, mm, isp1) = - t3jp( l1, mm, isp1)

c  j = l-1/2
        t3jm( l1, mm, isp1) = sqrt( j3m + 1.0e0 ) *
     1                 real( cwig3j( j1, j2, j3m, m1, m2, 2) )
        if (mod( (j2-j1-m1-m2)/2 , 2) .ne.0)
     1          t3jm( l1, mm, isp1) = - t3jm( l1, mm, isp1)
  70  continue

      return
c  end subroutine xprep
      end
