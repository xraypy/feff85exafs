      subroutine atheap(nat, rat, iphat, ra)

c--------------------------------------------------------------
c  copyright 1993 university of washington         bruce ravel
c  modified by alexei ankudinov in march 1999
c--------------------------------------------------------------
      implicit real (a-h,o-z)
      implicit integer (i-n)
c      implicit double precision (a-h,o-z)
c-------------------------------------------------------------------
c  heapsort adapted from numerical recipes.  sort atoms by distance.
c  all the pesky little do loops are for transferring rows
c  of temp into toss.
c-------------------------------------------------------------------
c  alexei ankudinov: needed to avoid unnecessary permutations when atoms
c  are at the same distance from the central atom, in order to comply 
c  feff document: the sample atom should be the nearest to absorber or
c  first in the list among equidistant
c  Add small contribution 10**-8 * number to the sorting variable ra
c  in order to achieve this.
c-------------------------------------------------------------------
c  natx:   dimension parameter from calling program
c-------------------------------------------------------------------
      dimension rat(3, nat), toss(3), iphat(nat)
      double precision ra(nat), dum 

      if (nat.lt.2) return

      l=0
      do 10 i=1,nat
         ra(i) = dble( rat(1,i)**2 + rat(2,i)**2 + rat(3,i)**2 ) +
     1           i*1.d-8
c        small addition at to prefer the old ordering
         if (l.eq.0 .and.i.gt.1) then
             if (ra(i).lt.ra(i-1)) l=1
         endif
  10  continue
c     check if array is already in order
      if (l.eq.0) return

      l  = nat/2+1
      ir = nat
 110  continue
         if (l.gt.1) then
            l = l-1
            do 120 index=1,3
               toss(index)=rat(index,l)
 120        continue
            itoss = iphat(l)
            dum = ra(l)
         else
            do 130 index=1,3
               toss(index) = rat(index,ir)
 130        continue
            itoss = iphat(ir)
            dum = ra(ir)
            do 140 index=1,3
               rat(index,ir) = rat(index,1)
 140        continue
            iphat(ir) = iphat(1)
            ra(ir) = ra(1)
            ir=ir-1
            if (ir.eq.1) then
               do 150 index=1,3
                  rat(index,1)=toss(index)
 150           continue
               iphat(1) = itoss
               ra(1) = dum
c              sort is finished
               goto 300
            endif
         endif
         i=l
         j=l+l

 160     if (j.le.ir) then
            if (j.lt.ir) then
               if ( ra(j) .lt. ra(j+1) ) then
                  j  = j + 1
               endif
            endif

            if ( dum .lt. ra(j) ) then
               do 170 index=1,3
                  rat(index,i) = rat(index,j)
 170           continue
               iphat(i) = iphat(j)
               ra(i) = ra(j)
               i=j
               j=j+j
            else
               j=ir+1
            endif
            goto 160
         endif

         do 180 index=1,3
            rat(index,i) = toss(index)
 180     continue
         iphat(i) = itoss
         ra(i) = dum

      goto 110
 300  continue

      return
c end subroutine atheap
      end
      subroutine getang(nclusx, rat, i, j, theta, phi)

c------------------------------------------------------------------
c  determine theta and phi polar angles of the vector between two
c  atom positions
c
c  inputs
c    rat:   (3,nclusx) x,y,z of all atoms in cluster
c    i, j:  indices of atoms at ends of vector Ri-Rj
c
c  outputs
c    theta: polar angle theta of vector Ri-Rj
c    phi:   polar angle phi of vector Ri-Rj
c------------------------------------------------------------------

      implicit real (a-h,o-z)
      implicit integer (i-n)

c       include 'dim.h'
c       include 'xparam.h'
      dimension rat(3,nclusx)
      parameter(tiny=1.e-7, zero=0.e0, pi=3.141592654)

      x = rat(1,i) - rat(1,j)
      y = rat(2,i) - rat(2,j)
      z = rat(3,i) - rat(3,j)
      r = sqrt(x**2 + y**2 + z**2)

c  this fails to calculate phi correctly for, as an example,
c  x=0.5e-7 and y=2e-7.  However, those numbers are below the
c  precision of the numbers stored in potph.bin.

      phi = zero
      theta  = zero
      if (i.ne.j) then
c           phi = atan2(y,x)
c        all of these conditionals will do the work for a machine that
c        cannot correctly handle a zero value for the second argument
c        of atan2
          if (abs(x).lt.tiny) then
              if (abs(y).lt.tiny) then
                  phi = zero
              elseif (y.gt.tiny) then
                  phi = pi/2
              else
                  phi = -pi/2
              endif
          else
              phi = atan2(y,x)
          endif
          if (r.gt.tiny) then
            if (z.le.-r) then
             theta = pi
            elseif ( z.lt.r) then
             theta = acos(z/r)
            endif
          endif
      endif

      return
c  end subroutine getang
      end


c====================================================================
      subroutine rotxan (lxp1, mxp1, betax, i, j, k)
      implicit real (a-h,o-z)

c     input:  lxp1, mxp1: lmax+1 & mmax+1, largest L states in matrix
c             betax is the rotation angle
c             i and j are the indeces of the atoms, thus denote
c                 which pair of atoms this is the rotation matrix for
c             k=0 for forward rotation, k=1 for backward rotation
c     output: drix(L,k,j,i) in common /rotx/
c+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+
c     adapted by BR from rot3i, version for genfmt by SIZ
c        new data structure for rotation matrices to accomodate
c        xanes calculation
c+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+
c     subroutine rot3 calculates rotation matrices for l = 0,lxp1-1

c     subroutine rot3 calculates the beta dependence of rotation
c     matrix elements using recursion of an iterated version of
c     formula (4.4.1) in edmonds.
c
c     first written:(september 17,1986) by j. mustre
c     version 2  (17 sep 86)
c     version 3  (22 feb 87) modified by j. rehr
c     version for genfmt, modified by s. zabinsky, Sept 1991
c     Initialized dri0.  Some elements may be used before being
c        initialized elsewhere -- rot3i needs to be carefully
c        checked.  S. Zabinsky, April 1993
c
c******************** warning****************************************
c     lxx must be at least lxp1 or overwriting will occur
c     nmax must be at least nm or overwriting will occur
c--------------------------------------------------------------------
c     notation dri0(l,m,n) =  drot_i(l'm'n')
c     l = l'+1, n' = n-l, m' = m-l, primes denoting subscripts
c     thus dri0(1,1,1) corresponds to the rotation matrix with
c     l' = 0, and n' and m' = 0; dri0(3,5,5) : l' = 2,n' = 2,m' = 2.
c--------------------------------------------------------------------

      include '../HEADERS/dim.h'
      include 'xparam.h'
c       include 'xstruc.h'
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
c#mn{
c common blocks for saving rotation matrices between xanes and rotxan
       integer    jsavx, jsav, jbmagk
       parameter (jsavx = 150, roteps = 1.e-12,jbmagk=-9999)
       dimension drisav(-lx:lx,-lx:lx,0:lx,jsavx), betsav(jsavx)
       integer   ldsav(jsavx), mdsav(jsavx)
       common /rotsav/  drisav, betsav, ldsav, mdsav, jsav
       save  /rotsav/
c#mn}

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
      parameter (one = 1, zero = 0)
c      needed for commented out diagnostic file
c      logical open
      parameter(lxx=24)
      parameter (pi = 3.14159 26535 89793 23846 26433e0)
      complex coni, dum 
      parameter (coni = (0,1))

c     dri0 is larger than needed for genfmt, but necessary for
c     this calculation algorithm.  Copy result into smaller
c     dri arrays (in common) at end of this routine.
      dimension  dri0 (lxx+1, 2*lxx+1, 2*lxx+1)

c#mn{
c  check whether a rotation matrix for this {beta(ileg),lxp1,mxp1} has
c  been calculated and saved.  If so, just use the saved value
       do 90 isav = 1, jsav
          if (betsav(isav).eq.jbmagk) go to 95
          if ((lxp1.eq.ldsav(isav)).and.(mxp1.eq.mdsav(isav)).and.
     $         (abs(betax-betsav(isav)).le.roteps) ) then
cc             print*, 'using drisav for ', isav, betax, lxp1, mxp1
             do 85 il = 0, lx
             do 85 m1 = -il, il
             do 85 m2 = -il, il
               drix(m2,m1,il,k,j,i)=cmplx(drisav(m2,m1,il,isav),zero)
 85          continue
             go to 770
          end if
 90    continue
 95    continue
c#mn}


c     initialize dri0
      do 150 in = 1, 2*lxx+1
        do 150 im = 1, 2*lxx+1
          do 150 il = 1, lxx+1
            dri0(il,im,in) = zero
 150  continue

      nm  = mxp1
      ndm = lxp1+nm-1
      xc  = cos(betax/2)
      xs  = sin(betax/2)
      s   = sin(betax)
      dri0(1,1,1) =  1
      dri0(2,1,1) =  xc**2
      dri0(2,1,2) =  s/sqrt(2*one)
      dri0(2,1,3) =  xs**2
      dri0(2,2,1) = -dri0(2,1,2)
      dri0(2,2,2) =  cos(betax)
      dri0(2,2,3) =  dri0(2,1,2)
      dri0(2,3,1) =  dri0(2,1,3)
      dri0(2,3,2) = -dri0(2,2,3)
      dri0(2,3,3) =  dri0(2,1,1)
      do 230  l = 3, lxp1
        ln = 2*l - 1
        lm = 2*l - 3
        if (ln .gt. ndm)  ln = ndm
        if (lm .gt. ndm)  lm = ndm
        do 220  n = 1, ln
          do 210  m = 1, lm
            t1   = (2*l-1-n) * (2*l-2-n)
            t    = (2*l-1-m) * (2*l-2-m)
            f1   = sqrt(t1/t)
            f2   = sqrt( (2*l-1-n) * (n-1) / t )
            t3   = (n-2) * (n-1)
            f3   = sqrt(t3/t)
            dlnm = f1 * xc**2 * dri0(l-1,n,m)
            if (n-1 .gt. 0) dlnm = dlnm - f2*s*dri0(l-1,n-1,m)
            if (n-2 .gt. 0) dlnm = dlnm + f3*xs**2*dri0(l-1,n-2,m)
            dri0(l,n,m) = dlnm
            if (n .gt. (2*l-3))
     $                  dri0(l,m,n) = (-1)**(n-m) * dri0(l,n,m)
 210      continue
          if (n .gt. (2*l-3)) then
              dri0(l,2*l-2,2*l-2) =  dri0(l,2,2)
              dri0(l,2*l-1,2*l-2) = -dri0(l,1,2)
              dri0(l,2*l-2,2*l-1) = -dri0(l,2,1)
              dri0(l,2*l-1,2*l-1) =  dri0(l,1,1)
          endif
 220    continue
 230  continue


c     initialize drix
      do 310 il = 0, lx
      do 310 m1 = -lx, lx
      do 310 m2 = -lx, lx
        drix(m2,m1,il,k,j,i) = cmplx(zero,zero)
        drix(m2,m1,il,k,i,i) = cmplx(zero,zero)
 310  continue

c     Copy result into drix(...,k,j,i) in /rotx/
      do 390  il = 1, lxp1
        mmx = min (il-1, mxp1-1)
        do 380  m1 = -mmx, mmx
        do 380  m2 = -mmx, mmx
          drix(m2, m1, il-1, k, j, i)=cmplx(dri0(il,m1+il,m2+il),zero)
 380    continue
 390  continue
c#mn{
c      save dri if there's room
       if (jsav.lt.jsavx) then
          jsav = jsav + 1
cc          print*, 'saving dri to ',  jsav, betax, lxp1, mxp1
          betsav(jsav) = betax
          ldsav(jsav)  = lxp1
          mdsav(jsav)  = mxp1
          do 720 il = 0, lx
          do 720 m1 = -il, il
          do 720 m2 = -il, il
            drisav(m2,m1,il,jsav) = real(drix(m2,m1,il,k,j,i))
 720      continue
       else
cc          print*, 'not saving dri to ',  betax, lxp1, mxp1
       end if
 770   continue
c#mn}

c-----test sum rule on d
c       if (idbg(1).eq.1) then
c           inquire(file='rotmat.dat', opened=open)
c           if (.not.open) then
c               iun = nxtunt(25)
c               open (iun,file='rotmat.dat',status='unknown')
c           endif
c           write(iun,*)'  '
c           write(iun,*)'atom #s : ',i,j
c           write(iun,*)  ' il, im, sum, beta'
c           write(iun,*) ' (drix(il,im,in,k,j,i),in = -il,il)'
c           do 880 il = 0,lxp1-1
c             do 870 im = -il,il
c               sum = 0
c               do 850 in = -il,il
c                 term = drix(in,im,il,k,j,i)
c                 sum = sum+term**2
c  850           continue
c               write(iun,860) il,im,sum,betax
c               write(iun,862) (drix(in,im,il,k,j,i),in = -il,il)
c  860          format(2i3,1x,f16.12,1x,f8.4)
c  862          format(5f14.6)
c  870         continue
c  880       continue
c c          close(iun)
c       endif
c-----end test------------------------

        do 920 il = 0, lx
        do 920 m1 = -il, il
          dum = coni * m1 * (xphi(i,j)-pi)
          if (k.eq.1) dum = -dum
          dum = exp( dum )
          do 910 m2 = -il, il
            if (k.eq.1) then
              drix(m2,m1,il,k,j,i) = drix(m2,m1,il,k,j,i) * dum
            else
              drix(m1,m2,il,k,j,i) = drix(m1,m2,il,k,j,i) * dum
            endif
 910       continue
 920     continue

      return
c  end subroutine rotxan
      end
c====================================================================
c#mn{
       subroutine rotint
       implicit real (a-h,o-z)
       include '../HEADERS/dim.h'
       include 'xparam.h'
c        include 'xstruc.h'
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
c#mn{
c common blocks for saving rotation matrices between xanes and rotxan
       integer    jsavx, jsav, jbmagk
       parameter (jsavx = 150, roteps = 1.e-12,jbmagk=-9999)
       dimension drisav(-lx:lx,-lx:lx,0:lx,jsavx), betsav(jsavx)
       integer   ldsav(jsavx), mdsav(jsavx)
       common /rotsav/  drisav, betsav, ldsav, mdsav, jsav
       save  /rotsav/
c#mn}

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
c initialize /rotsav/
       jsav = 0
       do 100 js = 1, jsavx
          betsav(js) = jbmagk
          ldsav(js)  = 0
          mdsav(js)  = 0
          do 90 il  = 0, lx
             do 80 m1 = -lx, lx
                do 70 m2 = -lx, lx
                   drisav(m2,m1,il,js) = 0
 70             continue
 80          continue
 90       continue
 100   continue
       return
c#mn}
       end
      subroutine sortat(iph0, nat, npot, iphat, iphx, rat, xrat)

      implicit real (a-h,o-z)
      implicit integer (i-n)
c--------------------------------------------------------------------
c  this subroutine sorts the atoms in xrat such that the first npot
c  entries are each a representative atom of a unique potential.  This
c  will mean that the upper left corner of the full MS matrix will
c  contain all of the information needed to compute the fine structure
c  and all of the electron densities.
c  NOTA BENE:  the atoms *must* have already been sorted by radial
c    distance!
c--------------------------------------------------------------------
c  input:
c    iph0:    potential index for central atom in LDOS (added by ala)
c                (iph0=0 for absorbing atom as the central atom)
c     nat:    number of atoms in cluster
c    npot:    number of unique potentials in cluster
c    iphat:   (nclusx) potential index of each atom in cluster as read
c             from geometry file
c    rat:     (3, nclusx) coordinates of each atom in cluster as read
c             from geometry file
c  output:
c    iphx:    (nclusx) potential index of each atom in cluster sorted
c             so that the first npot+1 entries are examples of each
c             ipot
c    xrat:    (3, nclusx) coordinates of each atom in cluster sorted
c             so that the first npot+1 entries are examples of each
c             ipot
c--------------------------------------------------------------------

      include '../HEADERS/dim.h'
      include 'xparam.h'
      dimension rat(3,natxx), xrat(3,nclusx)
      integer   iphat(natxx), iphx(nclusx), ipoint
      dimension ipoint(0:nphasx)
      integer iph0, ip, ilast

      do 10 i=0,nphasx
        ipoint(i) = 0
 10   continue
      do 30 ic=1,nat
        iphx(ic) = iphat(ic)
        do 20 ix=1,3
          xrat(ix,ic) = rat(ix,ic)
 20     continue
 30   continue

c     (iph0=0 for absorbing atom as the central atom)
      if (iphx(1).ne.iph0) then
          call wlog('* * * ERROR in sortat * * *')
          call wlog('            The first atom in xrat is not '//
     $                'the central atom.')
          call wlog('            Complain to Bruce immediately!')
          call par_stop('SORTAT-1')
      endif

c       if (idbg(4).eq.1) print*,'SORTAT: nat,npot: ',nat,npot
c       if (idbg(4).eq.1) print*,'SORTAT: xcen,ycen,zcen: ',
c      $            xcen,ycen,zcen

c --- find the example of each unique potential that is closest to the
c     central atom.  This will presumably be well within the cluster
c     that was used to compute the overlapped potentials
      ipoint(iph0) = 1
      do 150 ip=0,npot
        if (ip .ne. iph0) then
          do 130 iat=2,nat
            if (iphx(iat).eq.ip .and. ipoint(ip).eq.0) then
                ipoint(ip) = iat
c                print*,'>>>>> ip, ipoint(ip)', ip, ipoint(ip)
            endif
 130      continue
        endif
 150  continue

c --- now swap the first few atoms with the atoms found above
      do 200 ip=0,npot

c ----- some potentials might not be in the xanes cluster
        if (ipoint(ip).eq.0) goto 200
c ----- don't swap two potentials if examples live in the first npot
c       entries
        if (ipoint(ip).le.ip+1) goto 200

        xx  = xrat(1,1+ip)
        yy  = xrat(2,1+ip)
        zz  = xrat(3,1+ip)
        iph = iphx(1+ip)

        xrat(1,1+ip) = xrat(1,ipoint(ip))
        xrat(2,1+ip) = xrat(2,ipoint(ip))
        xrat(3,1+ip) = xrat(3,ipoint(ip))
        iphx(1+ip)  = iphx(ipoint(ip))

        xrat(1,ipoint(ip)) = xx
        xrat(2,ipoint(ip)) = yy
        xrat(3,ipoint(ip)) = zz
        iphx(ipoint(ip))  = iph

c       added by ala
c       check that substituted atom was not some ip example
c          ???BR Jan 16 1998???
        do 190 ipp = ip+1, npot
          if (ipoint(ipp).eq.ip+1) ipoint(ipp) = ipoint(ip)
  190   continue
c       set the correct pointer to ip example
        ipoint(ip) = ip+1

 200  continue

c     added by ala
c     Notice that fms will take the last atom of given type ip
c     from first npot atoms in the list as an example for ip.
c     Make more permutaions if necesary.
      ilast = -1
      nmin = min (npot+1, nat)
      do 210 ip = 0, npot
        if (ipoint(ip).ne.0) then
          do 205 iat = 1,nmin
  205     if (iphx(iat).eq.ip) ilast = iat

          if (ilast.ne.ipoint(ip)) then
            xx  = xrat(1,ilast)
            yy  = xrat(2,ilast)
            zz  = xrat(3,ilast)

            xrat(1,ilast)= xrat(1,ipoint(ip))
            xrat(2,ilast)= xrat(2,ipoint(ip))
            xrat(3,ilast)= xrat(3,ipoint(ip))

            xrat(1,ipoint(ip)) = xx
            xrat(2,ipoint(ip)) = yy
            xrat(3,ipoint(ip)) = zz
c           now ipoint(ip) = ilast, but don't need ipoint anymore
          endif
        endif
  210 continue

c       if (idbg(4).eq.1) then
c           do 220 i=1,npot+1
c             print *,i,xrat(1,i),xrat(2,i),xrat(3,i),iphx(i)
c  220      continue
c       endif
      return
c  end subroutine sortat
      end
      subroutine xanlm(lmaxp1,mmaxp1)

c------------------------------------------------------------------
c  calculate and store all of the legendre polynomial normalization
c  factors needed in the problem
c     xnlm= sqrt ((2l+1)(l-m)!/(l+m)!)
c  see, for instance, Arfken section 12.6.  Note that this lacks the
c  factor of sqrt(4*pi)
c
c  inputs:
c     lmaxp1, nmaxp1:  maximun l and m considered in the problem +1
c                      i.e. lmaxp1 = l_max+1
c
c  outputs:
c     all normalization factors passed in common /xnlm/
c------------------------------------------------------------------

      implicit real(a-h,o-z)
c       parameter(ltot=6,mtot=3)

      include '../HEADERS/dim.h'
      include 'xparam.h'
c       include 'xstruc.h'
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

      common/afctr/afac,flzero,flg(0:50)
c      common/afctr/afac,flzero,flg(0:210)
c      common/afctr/afac,flzero,flg(0:110) vax change

      call xfctst
      do 50 il=1,lmaxp1
        mmxp1 = min(mmaxp1,il)
        do 40 im=1,mmxp1
          l    = il-1
          m    = im-1
          cnlm = (2*l+1) * flg(l-m) / flg(l+m)
          cnlm = sqrt(cnlm) * afac**m
          xnlm(m,l) = cnlm

 40     continue
 50   continue
      return
c  end subroutine xlm
      end


      subroutine xfctst
c  same as feff's factst, but with a different name
      implicit real (a-h,o-z)
c     program for s3j and s6j symbols obtained from
c     steve younger of n.b.s.   modified by j.b.mann
c--------------------------------------------------------------------
c     a set to 1/64 to prevent overflow on vax
c     range on  flg set to 0:210, rather than flg(210)
c--------------------------------------------------------------------
cBR   This allows calculation of a large factorial (~100) without
cBR   overflow problems -- factor in a power of a small number then
cBR   factor it out
c--------------------------------------------------------------------
      common /afctr/ a, flzero, flg(0:50)
c      common /afctr/ a, flzero, flg(0:210)
      a=0.03125
c     a=0.015625
      flzero = 1.0
      flg(0) = 1.0
      flg(1) = a
      do 10 i=2,50
        flg(i) = flg(i-1) * i * a
 10   continue
      return
      end



c====================================================================
