      subroutine yprep(iph0, nat, inclus, iphat, rmax, rat)
c
c  removed izx, nopt, rdirec from arg list, unused or use commented out
c
c    yprep is the same as xprep for negative idwopt
c    simlifies calls in SCF and LDOS where DW factors should not enter

      implicit real (a-h,o-z)
      implicit integer (i-n)

      include '../HEADERS/dim.h'
      include '../HEADERS/const.h'
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
c
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

c  end of xstruc.h
c********************************************************************
      integer   iphat(natxx), iphat2(natxx)
c      integer izpair(0:2), izx(0:nphasx)
      dimension rat(3,natxx), rat2(3,natxx)
      double precision ra(natxx)
      character*78 line
c     sigms is written in double precision.  these are the variables
c     that it uses
c      double precision dtemp, dthet, drs, dsigsq, pair(3,0:2)
c      double precision sig2mx, sig2x(0:nphx,0:nphx)
c     iwarn - needed to wrtite waqrning just one time
      integer iwarn
      save iwarn
      data  iwarn /0/

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
                call par_stop('YPREP-1')
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
        if (rr.gt.rmax2) then
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
          rr = (xrat(1,i)-xrat(1,j))**2 + (xrat(2,i)-xrat(2,j))**2
     1       + (xrat(3,i)-xrat(3,j))**2
c         if (rr.gt.rdirec**2) goto 140

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

      do 200 iat2=1,nclusx
      do 200 iat1=1,nclusx
        sigsqr(iat1,iat2) = zero
 200  continue

      return
      end
