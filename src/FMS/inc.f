c///////////////////////////////////////////////////////////////////////
c FEFF PROGRAMS (referred below as a System)
c Copyright (c) 1986-2002, University of Washington.
c 
c END-USER LICENSE 
c 
c A signed End-user License Agreement from the University of Washington
c Office of Technology Transfer is required to use these programs and
c subroutines.
c 
c See the URL: http://leonardo.phys.washington.edu/feff/
c 
c USE RESTRICTIONS:
c 
c 1. The End-user agrees that neither the System, nor any of its
c components shall be used as the basis of a commercial product, and
c that the System shall not be rewritten or otherwise adapted to
c circumvent the need for obtaining additional license rights.
c Components of the System subject to other license agreements are
c excluded from this restriction.
c
c 2. Modification of the System is permitted, e.g., to facilitate
c its performance by the End-user. Use of the System or any of its
c components for any purpose other than that specified in this Agreement
c requires prior approval in writing from the University of Washington.
c
c 3. The license granted hereunder and the licensed System may not be
c assigned, sublicensed, or otherwise transferred by the End-user.  
c
c 4. The End-user shall take reasonable precautions to ensure that
c neither the System nor its components are copied, or transferred out
c side of his/her current academic or government affiliated laboratory
c or disclosed to parties other than the End-user.
c 
c 5. In no event shall the End-user install or provide this System
c on any computer system on which the End-user purchases or sells
c computer-related services.
c 
c 6. Nothing in this agreement shall be construed as conferring rights
c to use in advertising, publicity, or otherwise any trademark or the
c names of the System or the UW.   In published accounts of the use or
c application of FEFF the System should be referred to  by this name,
c with an appropriate literature reference:
c 
c FEFF8: A.L. Ankudinov, B. Ravel, J.J. Rehr, and S.D. Conradson,
c        Phys. Rev. B 58, pp. 7565-7576 (1998).
c
c LIMITATION OF LIABILITY:
c
c 1.   THE UW MAKES NO WARRANTIES , EITHER EXPRESSED OR IMPLIED, AS TO
c THE CONDITION OF THE SYSTEM, ITS MERCHANTABILITY, OR ITS FITNESS FOR
c ANY PARTICULAR PURPOSE.  THE END-USER AGREES TO ACCEPT THE SYSTEM
c 'AS IS' AND IT IS UNDERSTOOD THAT THE UW IS NOT OBLIGATED TO PROVIDE
c MAINTENANCE, IMPROVEMENTS, DEBUGGING OR SUPPORT OF ANY KIND.
c
c 2. THE UW SHALL NOT BE LIABLE FOR ANY DIRECT, INDIRECT, SPECIAL,
c INCIDENTAL OR CONSEQUENTIAL DAMAGES SUFFERED BY THE END-USER OR ANY
c OTHER PARTIES FROM THE USE OF THE SYSTEM.
c
c 3.  The End-user agrees to indemnify the UW for liability resulting
c from the use of the System by End-user. The End-user and the UW each
c agree to hold the other harmless for their own negligence.
c
c TITLE:
c
c 1.  Title patent, copyright and trademark rights to the System are
c retained by the UW. The End-user shall take all reasonable precautions
c to preserve these rights.
c 
c 2.  The UW reserves the right to license or grant any other rights to
c the System to other persons or entities.
c///////////////////////////////////////////////////////////////////////
c License is applicable for routines below, until otherwise specified.
c
      subroutine fmsie( iph0, nph, lipotx, ie, em, eref, ph, iz,
     1                 rfms, lfms, nat, iphat, rath, gtr)

c     full multiple scattering code for single energy point
c     written by a.ankudinov 06.1997 using earlier written subroutines
c     coded by b.ravel
c     modified by a.ankudinov 2001 for new matrix inversion algorithms
c     Feb. 2002, a.ankudinov: fixed logic for MPI calculations
c       lfms=0  - extended system calculataions (e.g. crystal)
c       lfms=1  - small system calculations (e.g. molecule)
c       lfms=2  - same as 1 for MPI run (forces call yprep)

      implicit double precision (a-h, o-z)
      include '../HEADERS/dim.h'

c     input
      dimension iphat(natx), rath(3,natx)
      real rat(3,natx), rfms, rdirec, toler1, toler2
      real rpart,aipart
      integer nph
      dimension iz(0:nphx)
      complex*16 ph(lx+1, 0:nphx)

c     work space
      integer iph0
      complex*16 em, eref
      character*512 slog
      logical lcalc
      dimension lcalc(0:lx)
c     fms staff
      integer lipotx(0:nphx)
      complex gg(nspx*(lx+1)**2, nspx*(lx+1)**2, 0:nphx)
      complex gtr(0:lx, 0:nphx)
      complex xphase(nspx, -lx:lx, 0:nphx), ck(nspx)
      complex*16 dck
      complex conis
      parameter (conis = (0,1))
      real  temper, thetax, sig2
      save

      if (rfms .le. 0.0) goto 900

c     set default (LU) inv method
      minv = 0
      rdirec = 2*rfms
      toler1 = 0.e0
      toler2 = 0.e0

      do 30 iat=1,nat
      do 30 j=1,3
   30 rat(j,iat) = real (rath(j,iat))

c     transform to single precision
      temper =0.0e0
      thetax =0.0e0
      sig2  = 0.0e0

c      it will be nice to call yprep once for all energy points,
c      fix later, and now call it every time
      if (ie.eq.1 .or. lfms.eq.0 .or. lfms.eq.2) 
     1  call yprep(iph0, nat, inclus, nph, iphat, rfms, rat,
     2     iz, rdirec )

      if (inclus.gt.1) then

cc     call fms for a cluster around central atom
       if (ie.eq.1) then
          write (slog,35) inclus, iph0
  35      format ('        Doing FMS for a cluster of ',i3,
     1    ' atoms around iph = ',i2)
          call wlog (slog)
       endif

       dck=sqrt(2*(em-eref))
       rpart  = dble(dck)
       aipart = real(dimag(dck))
       ck(1) = cmplx(rpart, aipart)
       do 1020 ipp = 0,nph
         do 1010 ill = -lipotx(ipp), lipotx(ipp)
           rpart  = dble( ph( 1+abs(ill), ipp))
           aipart = dimag(ph( 1+abs(ill), ipp)) 
           xphase(1, ill, ipp) = cmplx(rpart, aipart)
 1010    continue
 1020  continue
       iverb=0
       if (ie.eq.1) iverb = 1
       nsp = 1
       ispin = 0
       do 1011 ill = 0,lx
 1011  lcalc(ill) = .true.
       call fms(lfms, nsp, ispin, inclus, nph, ck, lipotx, xphase, ie,
     1  iverb, minv, rdirec, toler1, toler2, lcalc, gg)

c      make ck= i, since coni is c*16
       do 1030 ip=0,nph
         if (lfms.ne.0 .or. ip.eq.iph0) then
           do 1040 il=0,lipotx(ip)
             ix = il**2
             do 1050 im=1,2*il+1
               gtr(il, ip) = gtr(il, ip) + gg(ix+im,ix+im,ip)
 1050        continue
             gtr(il,ip)= gtr(il,ip)*
     1            exp(2*conis*xphase(1,il,ip))/(2*il+1)
 1040      continue
         endif
 1030  continue
      endif

 900  continue
      return
      end
      subroutine fms(lfms, nsp, ispin, inclus, npot, ck, lipotx, xphase,
     1   ik, iverb, minv, rdirec, toler1, toler2, lcalc, gg)

      implicit real (a-h,o-z)
      implicit integer (i-n)
c--------------------------------------------------------------------
c  compute full multiple scattering within some cluster at some
c  energy 
c  This uses the LU decomposition package from LAPACK.  Driver
c  routines: cgetrf (decomposition), cgetrs (backsubstitution)
c  coded by b.ravel
c  modified by a.l.ankudinov to include spin and SO interactions
c  feb 2000
c
c  dim.h and xparam.h must be included in the calling routine
c
c  most of the information needed by this package is set into common
c  blocks the companion package xprep.  In that package, the lists of
c  atomic coordinates and potential indeces are organized so that the
c  first npot+1 entries are examples of each of the unique potentials.
c  Consequently, only the upper left hand corner of the FMS matrix
c  need be recomposed to get the set of submatrices necessary to
c  compute chi for each type of atom in the cluster.
c
c  See subroutine fmstot.f for an example of decoding the output of this
c  subroutine. The third index of gg refers to the unique potential with
c  element 0 being the absorbing atom.  
c  The first two indeces are related to the |lms> state by the
c  formula:
c       nsp=1, no spin indeces
c       lm  = ( l**2 + 1 ) + ( l + m )
c            thus {1..(lx+1)^2} ==>
c            {0,0 1,-1 1,0 1,1 2,-2 2,-1 2,0 2,1 2,2 ...
c                   lx,lx-1 lx,lx}
c       nsp=2, with spin indeces
c       lms  = 2*( l**2 + 1 ) + 2*( l + m ) + (s-1/2)
c            thus {1...2*(lx+1)^2} ==>
c            {0, 0,-1/2  0. 0,1/2
c             1,-1,-1/2  1,-1,1/2  1,0,-1/2  1,0,1/2  1,1,-1/2 1,1,1/2
c             2,-2,-1/2  2,-2,1/2  2,-1,-1/2 2,-1,1/2 ...    lx,lx,1/2}
c
c  The calling protocol for xpreppack and fmspack is;
c          include 'dim.h'
c          include 'xparam.h'
c          ...
c          call xprep(nat, inclus, npot, iphat, rmax, rat,
c     $            xnrm, izx, temper, thetad)
c          energy loop {
c             ...
c             call fms(nsp, inclus, npot, ck, lipotx, xphase,
c                      ik, iverb, gg)
c             ... }
c
c  fmspack contains the following routines:
c    fms.f:     main routine of fmspack
c    kets.f:    compute all state kets for current energy
c    xclmz.f:   compute hankle-like polynomials for current energy
c    xgllm.f:   compute z-axis propagators for current energy
c    cgetrf.f:  LU decomposition of matrix
c    cgetrs.f:  backsubstitution of LU decomposed matrix
c    lu_misc.f: various routines called by LU package
c
c---------------------------------------------------------------------
c  input
c    nsp:    1) no spin indeces 2) with spin indeces
c    inclus: number of atoms in cluster
c    npot:   number of unique potentials in cluster
c    ck:     complex momentum of current energy point
c    lipotx: (0:nphasx) max l for each unique potential
c    xphase: (0:lx, 0:nphasx) single complex array of partial wave
c            phase shifts for each unique potential
c    ik:     current energy grid index, used for run-time messages
c    iverb:  do nothing when iverb <= 0
c            1  => write a message about grid point and matrix size
c
c  passed in common from xprep package (xstruc.h)
c    xrat:   (3,nclusx) array of coordinates with first npot+1
c            elements each a unique potential
c    xphi:   (nclusx, nclusx) angles between z axis and vectors
c            connecting the atoms in the cluster
c    iphx:   (nclusx) potential index of each atom in the cluster
c    drix:   huge matrix containing all rotation matrix elements
c            needed for computation of free electron propagators
c    xnlm:   matrix of legendre polynomial normalization factors
c    xpsile: matrix containing wave functions for hybridization
c            calculation
c    sigsqr: (nclusx,nclusx) matrix of pair-wise mean square
c            displacements about interatomic distances.  Currently only
c            calculated by the correlated debye model.
c
c  output
c    gg:  (nsp*lx**2, nsp*lx**2, 0:nphasx) submatrix spanning the entire
c          angular momentum basis for each unique potential

      include '../HEADERS/dim.h'
      include '../HEADERS/parallel.h'
      include 'xparam.h'
c====================================================================
c  This common block contains the structural information about the
c  cluster to be used for the full multiple scattering calculation
c  xphi:  matrix of angles between the z axis and pairs of atoms
c  xrat:  xyz coordinates of the atoms in the cluster, the first
c         npot+1 entries are examples of each unique potential
c  iphx:  potential indeces of each atom in the cluster, ordered like
c         xrat
      common /xstruc/ xphi(nclusx,nclusx), xrat(3,nclusx),
     $            iphx(nclusx)
      save /xstruc/
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
c********************************************************************
c**** save Clebsch-Gordon coefficients: <LS|J>
      dimension t3jp(0:lx, -lx:lx, 2), t3jm(0:lx, -lx:lx, 2)
      common /t3j/ t3jp, t3jm
      save   /t3j/

c********************************************************************
      parameter (pi = 3.14159 26535 89793 23846 26433e0)
      parameter (bohr = 0.529 177 249e0)
      parameter (one = 1, zero = 0)
      complex coni
      parameter (coni = (0,1))
      complex   term, prefac, gllmz, ck(nspx)
      complex   clm(lx+2, 2*lx+3), xclm(0:lx, 0:lx, nclusx, nclusx,nspx)
      complex   xrho( nclusx, nclusx, nspx)
      integer   lipotx(0:nphasx)

c********************************************************************
c**** array of state kets at current energy
      common /stkets/ lrstat(4, istatx), istate
      save   /stkets/
      complex   xphase(nspx, -lx:lx, 0:nphasx)
      complex   tmatrx(nspx, istatx)
c     big work matrices
      complex   g0(istatx,istatx), g0t(istatx,istatx)
      logical lcalc
      dimension lcalc(0:lx)
c     return matrix containing info about each unique potential
      complex   gg(nspx*(lx+1)**2, nspx*(lx+1)**2, 0:nphasx)

      integer i0(0:nphx)
      character*3  cerr, dec
      character*13 trans
      character*75 messg

 400  format(i4)

      do 10 i=0,nphx
        if (lipotx(i).le.0)  lipotx(i) = lx
        if (lipotx(i).gt.lx) lipotx(i) = lx
        i0(i) = -1
 10   continue
c     initialize gg to zero
      do 20 i = 0, nphasx
        do 18 j = 1, nspx*(lx+1)**2
          do 16 k = 1, nspx*(lx+1)**2
            gg( k, j, i) = cmplx( zero, zero)
 16       continue
 18     continue
 20   continue

      if (lfms.eq.0) then
        ipi = iphx(1)
        ipf = iphx(1)
      else
        ipi = 0
        ipf = npot
      endif
c --- get basis kets; output array 'lrstat' passed through common
      call getkts(nsp, inclus, npot, iphx, lipotx, i0)

c --- sanity check for i0(ip)
      do 30 ip = ipi, ipf
        if (i0(ip) .lt. 0) then
          call wlog (' Cannot find all representative atoms')
          call wlog (' Increase FMS radius and rerun.')
          call par_stop(' In subroutine FMS')
        endif
  30  continue

c --- runtime message if requested
      if (iverb.gt.0 .and. minv.eq.0) then
         dec = 'LUD'
         write(messg, 4010)this_process,dec, ik, istate
 4010    format('  ',i3,'   FMS matrix (', a, ') at point ', i3,
     $               ', number of state kets =', i4)
         call wlog(messg)
      endif

c --- get all c_lm(z) values for this energy, i,j sum over all atom
c     pairs xrho and xclm are symmetric in ij
c+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
c  nota bene, in the code for setting the clmz, the indexing starts
c  at 1 rather than 0.  To my mind, that is confusing, so here I
c  reindex when I copy from clm to xclm.  See the note about this in
c  subroutine xclmz
c+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
      lplus1 = lx+1
      mplus1 = lx+1
      do 140  i=1,inclus
        do 130 j=1,i

c ------- get and store rho for this pair of atoms   bohr units on
c         r and ck
          r   = zero
          do 100 ix=1,3
            r = r + (xrat(ix,i) - xrat(ix,j))**2
 100      continue
          r   = sqrt(r)

          do 125 isp = 1,nsp
             xrho(i,j,isp) = ck(isp) * r
             xrho(j,i,isp) = xrho(i,j,isp)

c ------- store the c_lm(z) for all the rhos at this energy
c            xclm(i,j) = xclm(j,i) by symmetry
             if (i.ne.j) call xclmz(lplus1,mplus1,xrho(i,j,isp),clm)
             do 120 ll = 0,lx
               do 110 mm = 0,lx
                 if (i.eq.j) then
                     xclm(mm,ll,j,i,isp) = cmplx(zero,zero)
                 else
                     xclm(mm,ll,j,i,isp) = clm(ll+1,mm+1)
                     xclm(mm,ll,i,j,isp) = clm(ll+1,mm+1)
                 endif
 110           continue
 120         continue
 125      continue

 130    continue
 140  continue

c --- fill the G0 and T matrices for this energy
      rdir2 = rdirec**2
      do 220 ist1=1,istate
        iat1 = lrstat(1, ist1)
        l1   = lrstat(2, ist1)
        m1   = lrstat(3, ist1)
        isp1 = lrstat(4, ist1)

        do 210 ist2=1,istate
          iat2 = lrstat(1, ist2)
          l2   = lrstat(2, ist2)
          m2   = lrstat(3, ist2)
          isp2 = lrstat(4, ist2)

          rr = (xrat(1,iat1)-xrat(1,iat2))**2 +
     1    (xrat(2,iat1)-xrat(2,iat2))**2 +(xrat(3,iat1)-xrat(3,iat2))**2

c                               equation 9 in Rehr, Albers
c                               <LR| G |L'R'>

          if (iat1.eq.iat2) then
c             same atom: G=0, calculate T-matrix 
              g0(ist1,ist2)     = cmplx(zero,zero)
c             notice that T is tri-diagonal, due to conservation of
c             total momentum.(will be broken by nonspherical potential)
c             --- potential index for this atom
              iph = iphx(iat1)
            if (nsp.eq.1.and.ispin.eq.0) then
              if (ist1.eq.ist2) tmatrx(1, ist1) =
     $                    ( exp(2*coni*xphase(isp1,l1,iph)) - one )
     $                    / (2*coni) 
            else
              if (ist1.eq.ist2) then
c                set spin index for t3jm and t3jp
                 is = isp1
                 if (nsp.eq.1) then
c                  special case
                   is = 1
                   if (ispin.gt.0) is = 2
                 endif

c                diagonal matrix element
                 tmatrx(1, ist1) =
     $                    ( exp(2*coni*xphase(isp1,l1,iph)) - one )
     $                    / (2*coni) * t3jm (l1, m1, is)**2  + 
     $                    ( exp(2*coni*xphase(isp1,-l1,iph)) - one )
     $                    / (2*coni) * t3jp (l1, m1, is)**2 
              elseif (nsp.eq.2.and.l1.eq.l2.and.m1+isp1.eq.m2+isp2) then
c                same orb. mom. and total momentum projections conserved
c                calculate off-diagonal T-matrix element
c                tmatrx(2, ist1) = here only if nspx equal to 2
                 tmatrx(nsp, ist1) =
     $             ( exp(2*coni*xphase(isp1, l1,iph)) - one +
     $               exp(2*coni*xphase(isp2, l1,iph)) - one ) / (4*coni) 
     1             * t3jm (l1, m1, isp1) * t3jm (l1, m2, isp2)  + 
     $             ( exp(2*coni*xphase(isp1,-l1,iph)) - one +
     $               exp(2*coni*xphase(isp2,-l1,iph)) - one ) / (4*coni) 
     1             * t3jp (l1, m1, isp1) * t3jp (l1, m2, isp2)
              endif
            endif
          elseif (isp1.eq.isp2 .and. rr.le.rdir2) then
c           different atoms, same spin: T=0, calculate G
            g0(ist1,ist2) = cmplx(zero,zero)
            do 200 mu=-l1,l1
c             --- third arg in drix: 0==>beta, 1==>-beta
              muabs = abs(mu)
              call xgllm(muabs, ist1, ist2, lrstat,
     1                   xclm(0,0,1,1,isp1), gllmz )
              g0(ist1,ist2) = g0(ist1,ist2) +
     2             drix(mu,m1,l1,1,iat2,iat1) *  gllmz *
     3             drix(m2,mu,l2,0,iat2,iat1)
 200        continue
            prefac = exp(coni*xrho(iat1,iat2,isp1)) /
     $                  xrho(iat1,iat2,isp1)
c           use correlated debye model, sigsqr is in AA^2
            prefac = prefac * exp(-1 * sigsqr(iat1,iat2) *
     $                  ck(isp1)**2 / bohr**2)
            g0(ist1,ist2) = prefac * g0(ist1,ist2)
          else
c           different atoms, different spins:T=G=0
            g0(ist1,ist2) = cmplx(zero,zero)
          endif

c -----   end of loops over states
 210    continue
 220  continue

      if (minv.eq.0) then
         call gglu ( nsp, i0, ipi, ipf, lipotx, g0, tmatrx, g0t, gg)
      elseif (minv.eq.1) then
         dec = 'VdV'
         call ggbi ( nsp, i0, ipi, ipf, lipotx, g0, tmatrx, g0t, gg,
     1              toler1, toler2, lcalc, msord)
      elseif (minv.eq.2) then
         dec = 'LLU'
         call ggrm ( nsp, i0, ipi, ipf, lipotx, g0, tmatrx, g0t, gg,
     1              toler1, toler2, lcalc, msord)
      elseif (minv.eq.3) then
         dec = 'GMS'
         call gggm ( nsp, i0, ipi, ipf, lipotx, g0, tmatrx, g0t, gg,
     1              toler1, toler2, lcalc, msord)
      else
         dec = 'TF'
         call ggtf ( nsp, i0, ipi, ipf, lipotx, g0, tmatrx, g0t, gg,
     1              toler1, toler2, lcalc, msord)
      endif
      if (minv.ne.0) then
         write(messg, 410)this_process,dec, ik, istate, msord
 410     format('  ',i3,'. Iterative FMS (', a, ') at point ', i3,
     $               '; matrix size =', i4,'; MS order =',i5)
         call wlog(messg)
      endif

      return
      end
c--------------------------------------------------------------------
      subroutine getkts(nsp, nat, npot, iphx, lipotx, i0)

      implicit real (a-h,o-z)
      implicit integer (i-n)
c--------------------------------------------------------------------
c  construct state kets |iat,l,m> at this energy
c--------------------------------------------------------------------
c  input
c    nat:    number of atoms in cluster
c    npot:   number of unique potentials
c    iphx:   (nclusx) potential index of each atom in the cluster
c    lipotx: (nphasx) maximum angular momentum to consider for each
c            ipot
c  output
c   (istate: number of states  ---  passed in kets.h)
c    i0:     index shift for each potential representative
c   (lrstat: (4, istatx) state kets |iat,l,m> --- passed in kets.h)
c--------------------------------------------------------------------
      include '../HEADERS/dim.h'
      include 'xparam.h'
c********************************************************************
c**** array of state kets at current energy
      common /stkets/ lrstat(4, istatx), istate
      save   /stkets/
      integer   lipotx(0:nphasx), iphx(nclusx), i0(0:nphx)

      istate = 0
      do 120 iat=1,nat
        ip = iphx(iat)
c       i0(ip) - index for the ip-representative atom
c       need for simple find of states for ip-representative.
        if (i0(ip).lt.0) i0(ip) = istate
        lim = min(lx, lipotx(ip))
        do 110 l=0,lim
          do 100 m = -l, l
          do 100 isp = 1, nsp
            istate = istate + 1
            if (istate.gt.istatx) then
                call wlog('Exceeded maximum number of LR states.'//
     $                      '  Stopping')
                call par_stop('GETKTS-1')
            endif
            lrstat(1,istate) = iat
            lrstat(2,istate) = l
            lrstat(3,istate) = m
            lrstat(4,istate) = isp
 100      continue
 110    continue
 120  continue

      return
c end subroutine kets
      end
c    ----------------------------------------------------------------
      subroutine xclmz(lmaxp1,mmaxp1,rho,clm)
      implicit real(a-h,o-z)

c     calculates energy dependent factors needed in subroutine gllm
c     c(il,im) = c_l^(m)z**m/m!=c_lm             by recursion
c     c_l+1,m  = c_l-1,m-(2l+1)z(c_l,m-c_l,m-1)  l ne m
c     c_m,m    = (-z)**m (2m)!/(2**m m!)         with z=1/i rho
c
c  input:
c    lmaxp1, mmaxp1:  largest angular momentum under consideration + 1
c    rho:  distance between atoms * complex momentum at this energy
c          point
c  output:
c    clm(lx+1,lx+1):  Hankle-like polynomials from RA

      include   '../HEADERS/dim.h'
      include   'xparam.h'
      parameter (one = 1, zero = 0)
      complex coni
      parameter (coni = (0,1))
      parameter (ltotb=lx+1,mtotb=ltotb,ntotb=ltotb,mntot=mtotb+ntotb)
      complex z, cmm, clm(ltotb+1,mntot+1), rho

      cmm  = cmplx(one, zero)
      z    = (-coni)/rho

      clm(1,1) = cmplx(one,zero)
      clm(2,1) = clm(1,1) - z

      lmax = lmaxp1-1
      do 20 il=2,lmax
        clm(il+1,1) = clm(il-1,1) - (z * (2*il-1) * clm(il,1))
 20   continue
c+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
c  nota bene:  the 2l-1 factor above is correct, even though in Rehr,
c  Albers equation 4 appears with a 2l+1.  The reason has to do with
c  the indexing.  in RA the subscripts on the c's start at 0.  In this
c  piece of code, the subscripts start at 1.  If you sub l-1 for
c      l, 2l+1 --> 2l-1
c+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+


      mmxp1 = min(lmaxp1, mmaxp1)
      do 40 im=2,mmxp1
        m    = im-1
        imp1 = im+1
        cmm  = (-cmm) * (2*m-1) * z
        clm(im,im)   = cmm
        clm(imp1,im) = cmm * (2*m+1) * (1-im*z)
        do 30 il=imp1,lmax
          clm(il+1,im) = clm(il-1,im) - (2*il-1) * z *
     $                               ( clm(il,im)+clm(il,m) )
c           l = il-1
c           clm(il+1,im) = clm(l,im) - (2*l+1) * z *
c      $                               ( clm(il,im)+clm(il,m) )
 30     continue
 40   continue

      return
c  end subroutine xclmz
      end
      subroutine xgllm(mu, ist1, ist2, lrstat, xclm, gllmz)
c--------------------------------------------------------------------
c  this calculates equations 11,12 from Rehr, Albers PRB v.41,#12,
c  p.8139,  the output is the G term in equation 9 from that paper
c
c  input:
c    mu:         abs val of magnetic state in sum in eqn 11 RA, mu>=0
c    ist1, ist2: state indices of mat. elem., first index of lrstat
c    lrstat:     (4,istatx,nkmin:nex) array of LR states
c    xclm:       (0:lx,0:lx,nclusx,nclusx) array of c_lm(z) for
c                present energy value
c  output:
c    gllmz:      g_ll'^|m|(z), for present state & energy, eqn 11 RA
c--------------------------------------------------------------------
c  this requires that N_lm normalization factors and c_lm(z)
c  polynomials have already been calculated.
c--------------------------------------------------------------------
      implicit real(a-h,o-z)
      implicit integer (i-n)

      include '../HEADERS/dim.h'
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

      parameter (zero=0.e0)
      integer    lrstat(4, istatx)
      complex xclm(0:lx, 0:lx, nclusx, nclusx), sum, gllmz
      complex gam, gamtl

      iat1     = lrstat(1, ist1)
      l1       = lrstat(2, ist1)
      iat2     = lrstat(1, ist2)
      l2       = lrstat(2, ist2)
      numax    = min(l1, l2-mu)

      sum      = cmplx(zero, zero)
      do 100 nu=0,numax
        mn    = mu+nu

c       bug for xnlm with nspx=2
        gamtl = (2*l1+1) * xclm(nu, l1, iat2, iat1) / xnlm(mu, l1)
        gam   = (-1)**mu * xclm(mn, l2, iat2, iat1) * xnlm(mu, l2)

        sum   = sum + gamtl * gam
 100  continue

      gllmz = sum

      return
c  end subroutine gllm
      end

c====================================================================
      subroutine gglu( nsp, i0, ipi, ipf, lipotx, g0, tmatrx, g0t, gg)

      implicit real (a-h,o-z)
      implicit integer (i-n)
c  output
c    gg:  (nsp*lx**2, nsp*lx**2, 0:nphasx) submatrix spanning the entire
c          angular momentum basis for each unique potential

      include '../HEADERS/dim.h'
      include 'xparam.h'
      integer  i0 (0:nphx),  lipotx(0:nphx)

      parameter (one = 1, zero = 0)
      integer   ipiv(istatx)

c**** array of state kets at current energy
      common /stkets/ lrstat(4, istatx), istate

      complex   tmatrx(nspx, istatx)
c     big work matrices
      complex   g0( istatx, istatx), g0t( istatx, istatx)
      complex   g0s( istatx, nspx*(lx+1)**2)
c     return matrix containing info about each unique potential
      complex   gg(nspx*(lx+1)**2, nspx*(lx+1)**2, 0:nphasx)

      character*3  cerr
      character*13 trans

 400  format(i4)


c -------------------- LU gg
c     multiply T and G0 matrices together, construct g0t = 1 - G0*T
c     notice that the signs below for g0t ARE correct since 1 is the
c     unit matrix
c     since t is tri-diagonal, this product can be computed in n^2 time
c     also fill up some work matrices for use in eigenvalue and
c     determinant calculations and elsewhere
      do 320 icol = 1,istate
        do 310 irow = 1,istate
c         T diagonal contribution
          g0t(irow, icol) = - g0(irow, icol) * tmatrx(1, icol)
c         T off-diagonal contribution
          l1   = lrstat(2,icol)
          m1   = lrstat(3,icol)
          isp1 = lrstat(4,icol)
          m2 = m1+isp1
          if (nsp.eq.2 .and. m2.gt.-l1+1 .and. m2.lt.l1+2) then
             ist2 = icol + (-1)**isp1
             g0t(irow, icol) = g0t(irow, icol)
     1                   - g0(irow, ist2) * tmatrx(nsp, icol)
          endif
 310    continue

        g0t(icol, icol) = g0t(icol, icol) + one

 320  continue

c --- invert matrix by LU decomposition
c     call cgetrf from lapack.  this performs an LU decomposition on
c     the matrix g0t = 1 - g0*T
      call cgetrf( istate, istate, g0t, istatx, ipiv, info )
      if (info.lt.0) then
          call wlog('    *** Error in cgetrf when computing G')
          write(cerr,400)abs(info)
          call wlog('        Argument #'//cerr//
     $                ' had an illegal value.')
      elseif (info.gt.0) then
          call wlog('    *** Error in cgetrf when computing G')
          write(cerr,400)info
          call wlog('        g0t('//cerr// ','//cerr//
     $                ') is exactly 0 -- '//
     $                'this matrix cannot be decomposed.')
      endif

c     now we want g_c = (g0t)^-1 * g0.  Rather than calculating
c     the inverse of g0t from the LU decomposition, we can compute
c     g0t^-1 * g0 directly by backsubstituting the columns of G0.
c     See sect. 2.3 in Numerical Recipes or LAPACK Users' Guide
c     sect. 2.3

c     third arg in number of output columns, istate for full
c     matrix, ipart(ik) for just the parts of the matrix needed
c     to contruct fine structure + DOS functions

      do 620 ip=ipi, ipf
        ipart = nsp*(lipotx(ip)+1)**2
        do 590 is1 = 1, istate
        do 590 is2 = 1, ipart
          g0s(is1,is2) = g0(is1, is2 + i0(ip))
  590   continue

        trans = 'NotTransposed'
        call cgetrs(trans, istate, ipart, g0t, istatx,
     $                ipiv, g0s, istatx, info)
        if (info.lt.0) then
            call wlog('    *** Error in cgetrf')
            write(cerr,400) abs(info)
            call wlog('        Argument #'//cerr//
     $              ' had an invalid value.')
        endif

c **** at this point g0s contains the full MS ****

c  pack FMS matrix into an nsp*(lx+1)^2 x nsp*(lx+1)^2 matrix for each
c  ipot

        do 600 is2=1,ipart
        do 600 is1=1,ipart
          gg( is1, is2, ip) = g0s( is1+i0(ip), is2)
 600    continue
 620  continue

      return
      end
      subroutine ggbi( nsp, i0, ipi, ipf, lipotx, g0, tmatrx, g0t, gg,
     1                 toler1, toler2, lcalc, msord)

      implicit real (a-h,o-z)
      implicit integer (i-n)
c  output
c    gg:  (nsp*lx**2, nsp*lx**2, 0:nphasx) submatrix spanning the entire
c          angular momentum basis for each unique potential
c     BiCGStab algorithm: Saad, Iterative methods for ..., p. 220 (1996)

      include '../HEADERS/dim.h'
      include 'xparam.h'
      integer  i0 (0:nphx),  lipotx(0:nphx)

      parameter (one = 1, zero = 0)
      complex coni
      parameter (coni = (0,1))

c**** array of state kets at current energy
      common /stkets/ lrstat(4, istatx), istate

      complex   tmatrx(nspx, istatx)
c     big work matrices
      complex   g0( istatx, istatx), g0t( istatx, istatx)
      logical lcalc
      dimension lcalc(0:lx)

c     Lanczos method variables
      complex xvec( istatx), yvec(istatx), avec(istatx), asve(istatx)
      complex rvec(istatx), pvec(istatx), svec( istatx)
      complex aa, dd, aw, wa, ww
      complex del, delp, omega, chi, psi
c     return matrix containing info about each unique potential
      complex   gg(nspx*(lx+1)**2, nspx*(lx+1)**2, 0:nphasx)

c      notice that in gglu we invert (1-Gt), but here (1-tG).
c     multiply T and G0 matrices together, construct g0t = 1 - T*G0
c     notice that the signs below for g0t ARE correct since 1 is the
c     unit matrix
c     since t is tri-diagonal, this product can be computed in n^2 time
c     also fill up some work matrices for use in eigenvalue and
c     determinant calculations and elsewhere
c     cycle over dimensions of matrix g0t
      do 10 icol = 1,istatx
      do 10 irow = 1,istatx
 10   g0t(irow,icol) = 0

      do 30 icol = 1,istate
        do 20 irow = 1,istate
c         T diagonal contribution T(irow, irow)
          if ( abs( g0(irow, icol)) .gt. toler2 )
     1    g0t(irow,icol)=g0t(irow,icol) - tmatrx(1,irow) * g0(irow,icol) 

c         T off-diagonal contribution T(ist2, irow) in tmatr(2,irow)
c         T off-diagonal contribution T(irow, ist2) in tmatr(2,ist2)
          l1   = lrstat(2,irow)
          m1   = lrstat(3,irow)
          isp1 = lrstat(4,irow)
          m2 = m1+isp1
          if (nsp.eq.2 .and. m2.gt.-l1+1 .and. m2.lt.l1+2) then
c            spin-flip contribution
             ist2 = irow + (-1)**isp1
             if ( abs( g0(ist2, icol)) .gt. toler2)
     1       g0t(irow, icol) = g0t(irow, icol)
     2                   - tmatrx(nsp, ist2) * g0(ist2, icol) 
          endif
 20     continue

        g0t(icol, icol) = g0t(icol, icol) + one
 30   continue

      do 920 ip=ipi, ipf
        ipart = nsp*(lipotx(ip)+1)**2
        do 910 is1 = 1, ipart
          is2 = is1+i0(ip)
          l1   = lrstat(2,is2)
          if (.not.lcalc(l1)) goto 910

c         start first tier with xvec=0
          istart = -1
          msord = 0
          do 40 is = 1, istate
          avec(is) = 0
  40      xvec(is) = 0

c         RESTART here if necessary
  50      continue
          istart = istart+1

          if (istart.gt.0) call matvec( istatx,istate,g0t,xvec,avec,1)
          do 90 is = 1,istate
 90       rvec(is) = - avec(is)
c         rvec = bvec - g0t*xvec , in our case bvec(is) = delta_{is,is2}
          rvec(is2) = rvec(is2) + 1
cc          Check convergence criteria: |r_n+1| < tol
            ipass = 1
            do 92 is = 1, istate
              if ( abs(real(rvec(is))).gt.toler1) goto 93
              if ( abs(aimag(rvec(is))).gt.toler1) goto 93
 92         continue
            ipass = 0
 93         continue
            if (ipass.eq.0) goto 700

          do 95 is = 1, istate
 95       pvec(is) = rvec(is)
          call matvec( istatx,istate,g0t,pvec,avec,1)
          msord = msord + 1

c         choose yvec that del and delp close to one
          call cdot( istatx, istate, avec, avec, aa)
          call cdot( istatx, istate, rvec, avec, wa)
          aw = real(wa) - coni* aimag(wa)
          call cdot( istatx, istate, rvec, rvec, ww)
          dd = aa*ww - aw*wa
          if (abs(dd/aa/ww) .lt.1.e-8) then
            do 96 is = 1,istate
  96        yvec(is) = rvec(is) / ww
          else
            ww = ( ww - aw ) / dd
            aa = ( wa - aa) / dd
            do 97 is = 1,istate
  97        yvec(is) = rvec(is) * aa + avec(is) * ww
          endif
          call cdot( istatx, istate, yvec, rvec, del)

c         it seems ran out of precision for nit>150
          nitx = 30
          do 500 nit = 0, nitx
            call cdot( istatx, istate, yvec, avec, delp)
            omega = del / delp

            do 120 is = 1, istate
 120        svec(is) = rvec(is) - omega * avec(is)
cc          Check convergence criteria: |s_n+1| < tol
            ipass = 1
            do 122 is = 1, istate
              if ( abs(real(svec(is))).gt.toler1) goto 123
              if ( abs(aimag(svec(is))).gt.toler1) goto 123
 122        continue
            ipass = 0
 123        continue
            if (ipass.eq.0)  then
              do 124 is = 1, istate
 124          xvec(is) = xvec(is) + omega*pvec(is)
              goto 700
            endif

            call matvec( istatx,istate,g0t,svec,asve,1)
            msord = msord + 1
            call cdot( istatx, istate, asve, asve, aa)
            call cdot( istatx, istate, asve, svec, wa)
            chi = wa / aa
            do 125 is = 1, istate
 125        xvec(is) = xvec(is) + omega*pvec(is) + chi*svec(is)

            do 130 is = 1, istate
 130        rvec(is) = svec(is) - chi* asve(is)

cc          Check convergence criteria: |r_n+1| < tol
            ipass = 1
            do 370 is = 1, istate
              if ( abs(real(rvec(is))).gt.toler1) goto 380
              if ( abs(aimag(rvec(is))).gt.toler1) goto 380
 370        continue
            ipass = 0
 380        continue
            if (ipass.eq.0) goto 700

c           prepare for next iteration
            call cdot( istatx, istate, yvec, rvec, del)
            psi = del / (delp * chi)

            do 135 is = 1, istate
 135        pvec(is) = rvec(is) + psi * (pvec(is)-chi*avec(is))
            call matvec( istatx,istate,g0t,pvec,avec,1)
            msord = msord + 1

 500      continue
c         restart since ran out of iterations
          goto 50

c         exit if tolerance has been achieved
 700      continue
c         print*, ' BI iterations:', nit + istart*nitx
c         end of BI iterations

c         at this point xvec = (1-tG)**-1 * bvec  with chosen tolerance
c         pack FMS matrix into an nsp*(lx+1)^2 x nsp*(lx+1)^2 matrix 
c         for each ipot
          do 800 is2=1,ipart
            gg( is2, is1, ip) = zero
            do 790 is = 1,istate
              gg( is2, is1, ip) = gg( is2, is1, ip) +
     1        g0( is2+i0(ip), is) * xvec(is)
 790        continue
 800      continue

 910    continue
 920  continue

      return
      end
      subroutine ggrm( nsp, i0, ipi, ipf, lipotx, g0, tmatrx, g0t, gg,
     1                 toler1, toler2, lcalc, msord)

      implicit real (a-h,o-z)
      implicit integer (i-n)
c  output
c    gg:  (nsp*lx**2, nsp*lx**2, 0:nphasx) submatrix spanning the entire
c          angular momentum basis for each unique potential

      include '../HEADERS/dim.h'
      include 'xparam.h'
      integer  i0 (0:nphx),  lipotx(0:nphx)

      parameter (one = 1, zero = 0)
      complex coni
      parameter (coni = (0,1))

c**** array of state kets at current energy
      common /stkets/ lrstat(4, istatx), istate

      complex   tmatrx(nspx, istatx)
c     big work matrices
      complex   g0( istatx, istatx), g0t( istatx, istatx)
      logical lcalc
      dimension lcalc(0:lx)

      complex xvec( istatx), xket( istatx), xbra( istatx)
      complex xketp(istatx), xbrap(istatx)
      complex zvec(istatx), rvec(istatx), svec(istatx)
      complex tket(istatx), tbra(istatx)
      double precision  dum1, dum2
      complex alphac, betac, aa, bb, yy, aac, bbc, gamma
      real alpha, beta
c     return matrix containing info about each unique potential
      complex   gg(nspx*(lx+1)**2, nspx*(lx+1)**2, 0:nphasx)

c      notice that in gglu we invert (1-Gt), but here (1-tG).
c     multiply T and G0 matrices together, construct g0t = 1 - T*G0
c     notice that the signs below for g0t ARE correct since 1 is the
c     unit matrix
c     since t is tri-diagonal, this product can be computed in n^2 time
c     also fill up some work matrices for use in eigenvalue and
c     determinant calculations and elsewhere
c     cycle over dimensions of matrix g0t
      do 10 icol = 1,istatx
      do 10 irow = 1,istatx
 10   g0t(irow,icol) = 0

      do 30 icol = 1,istate
        do 20 irow = 1,istate
c         T diagonal contribution T(irow, irow)
          if ( abs( g0(irow, icol)) .gt. toler2 )
     1    g0t(irow,icol)=g0t(irow,icol) - tmatrx(1,irow) * g0(irow,icol) 

c         T off-diagonal contribution T(ist2, irow) in tmatr(2,irow)
c         T off-diagonal contribution T(irow, ist2) in tmatr(2,ist2)
          l1   = lrstat(2,irow)
          m1   = lrstat(3,irow)
          isp1 = lrstat(4,irow)
          m2 = m1+isp1
          if (nsp.eq.2 .and. m2.gt.-l1+1 .and. m2.lt.l1+2) then
c            spin-flip contribution
             ist2 = irow + (-1)**isp1
             if ( abs( g0(ist2, icol)) .gt. toler2)
     1       g0t(irow, icol) = g0t(irow, icol)
     2                   - tmatrx(nsp, ist2) * g0(ist2, icol) 
          endif
 20     continue

        g0t(icol, icol) = g0t(icol, icol) + one
 30   continue

      do 920 ip=ipi, ipf
        ipart = nsp*(lipotx(ip)+1)**2
        do 910 is1 = 1, ipart
          is2 = is1+i0(ip)
          l1   = lrstat(2,is2)
          if (.not.lcalc(l1)) goto 910

c         start first tier with xvec=0
          istart = -1
          msord=0
          do 40 is = 1, istate
          rvec(is) = 0
  40      xvec(is) = 0

c         RESTART here if necessary
  50      continue
          istart = istart+1

          if (istart.gt.0) call matvec( istatx,istate,g0t,xvec,rvec,1)
c         rvec = g0t*xvec - bvec, in our case bvec(is) = delta_{is,is2}
          rvec(is2) = rvec(is2) - 1
          do 90 is = 1,istate
 90       xket(is) = - rvec(is)
          call cdot( istatx, istate, xket, xket, bb)
          if (abs(bb).eq.0) goto 700

          xfnorm = 1.e0 / real(dble(bb))
          do 91 is = 1, istate
 91       xbra(is) = xket(is) * xfnorm
c         |t> = A |n> ; |n> - xket, |n-1> - xketp
          call matvec ( istatx, istate, g0t, xket, tket, 1)
          msord = msord + 1
          call cdot( istatx, istate, xbra, tket, aa)
          aac = real(aa) - coni*aimag(aa)
          bb = 0
          bbc= 0
          betac = aa
          yy = 1
c         initialize vectors
          do 110 is = 1,istate
            xketp(is) = 0
            xbrap(is) = 0
            zvec(is) = xket(is)
            xvec(is) = xvec(is) + zvec(is)/betac
 110      continue

          do 120 is = 1, istate
 120      svec(is) = tket(is)
          do 130 is = 1, istate
 130      rvec(is) = rvec(is) + svec(is) / betac

c         it seems ran out of precision for nit>150
          nitx = 100
          do 300 nit = 1, nitx
c           use recursion method to calculate a_n+1, b_n, |n+1>, <n+1|
            do 140 is = 1, istate
 140        tket(is) = tket(is) - aa*xket(is) - bb*xketp(is)
            call matvec ( istatx, istate, g0t, xbra, tbra, 2)
            do 150 is = 1, istate
 150        tbra(is) = tbra(is) - aac*xbra(is) - bbc*xbrap(is)
            call cdot( istatx, istate, tbra, tket, bb)
            if (abs(bb).eq.0) goto 700

            bb = sqrt (bb)
            bbc = real(bb) - coni*aimag(bb)
            do 160 is = 1, istate
              xketp(is) = xket(is)
              xbrap(is) = xbra(is)
 160        continue
            do 170 is = 1, istate
              xket(is) = tket(is) / bb
              xbra(is) = tbra(is) / bbc
 170        continue
            call matvec ( istatx, istate, g0t, xket, tket, 1)
            msord = msord + 1
            call cdot( istatx, istate, xbra, tket, aa)
            aac = real(aa) - coni*aimag(aa)
            
c           update iterative solution xvec, 
c           and residual rvec = g0t*xvec - |1>
            alphac = bb / betac
            do 210 is = 1, istate
 210        zvec(is) = xket(is) - alphac * zvec(is)
            do 220 is = 1, istate
 220        svec(is) = tket(is) - alphac * svec(is)

            betac = aa - alphac*bb
            yy = - alphac * yy
            gamma = yy / betac
            do 230 is = 1, istate
 230        xvec(is) = xvec(is) + gamma * zvec(is)
            do 240 is = 1, istate
 240        rvec(is) = rvec(is) + gamma * svec(is)

cc          Check convergence criteria: | rvec | < tol
c           call vecvec( istatx, istate, rvec, rvec, dum2)
c           if (dum2.le.tol) goto 700
cc          Check convergence criteria: | rvec | < tol
            ipass = 1
            do 250 is = 1, istate
              if ( abs(real(rvec(is))).gt.toler1) goto 260
              if ( abs(aimag(rvec(is))).gt.toler1) goto 260
 250        continue
            ipass = 0
 260        continue
            if (ipass.eq.0) goto 700

 300      continue
c         restart since ran out of iterations
          goto 50

c         exit if tolerance has been achieved
 700      continue
c         end of RM iterations

c         at this point xvec = (1-tG)**-1 * bvec  with chosen tolerance
c         pack FMS matrix into an nsp*(lx+1)^2 x nsp*(lx+1)^2 matrix 
c         for each ipot
          do 800 is2=1,ipart
            gg( is2, is1, ip) = zero
            do 790 is = 1,istate
              gg( is2, is1, ip) = gg( is2, is1, ip) +
     1        g0( is2+i0(ip), is) * xvec(is)
 790        continue
 800      continue

 910    continue
 920  continue

      return
      end

      subroutine cdot ( istatx, istate, abra, aket, cc)
c     dot product of two vectors
c     notice that we keep bra  vector as it's complex conjugate
c     thus need to conjugate abra here
      implicit real (a-h,o-z)
      implicit integer (i-n)
      complex coni
      parameter (coni = (0,1))
      complex abra, aket, cc, aa
      dimension abra(istatx), aket(istatx)

      cc = 0
      do 10 is = 1,istate
        aa = real(abra(is)) - coni*aimag(abra(is))
        cc = cc + aa * aket(is)
 10   continue
      return
      end

      subroutine vecvec ( istatx, istate, avec, bvec, cc)
c     dot product of two vectors
      implicit real (a-h,o-z)
      implicit integer (i-n)
      complex avec, bvec
      double precision cc, aa, bb
      dimension avec(istatx), bvec(istatx)

      cc = 0
      do 10 is = 1,istate
        aa = dble(real(avec(is))) * dble(real(bvec(is)))
        bb = dble(aimag(avec(is))) * dble(aimag(bvec(is)))
        cc = cc + aa + bb
 10   continue
      return
      end

      subroutine matvec (istatx, istate, amat, bvec, cvec, itrans)
c     itrans = 1  cvec = amat * bvec
c     itrans = 2  cvec = amat^+ * bvec
c     itrans = 3  cvec = amat^T * bvec
      implicit real (a-h,o-z)
      implicit integer (i-n)
      complex coni, aa
      parameter (coni = (0,1))
      complex amat, bvec, cvec
      dimension amat(istatx, istatx), bvec(istatx), cvec(istatx)

c     initialize cvec
      do 10 is = 1,istatx
 10   cvec(is) = 0

c     cycle over dimensions of amat
      do 20 icol = 1,istate
      do 20 irow = 1,istate
        if (itrans.eq.1) then
          cvec(irow) = cvec(irow) + amat(irow, icol) * bvec(icol)
        elseif(itrans.eq.2) then
          aa = real(amat(irow, icol)) - coni*aimag(amat(irow, icol))
          cvec(icol) = cvec(icol) + aa * bvec(irow)
        else
          cvec(icol) = cvec(icol) + amat(irow, icol) * bvec(irow)
        endif
 20   continue

      return
      end
      subroutine gggm( nsp, i0, ipi, ipf, lipotx, g0, tmatrx, g0t, gg,
     1                 toler1, toler2, lcalc, msord)

      implicit real (a-h,o-z)
      implicit integer (i-n)
c  output
c    gg:  (nsp*lx**2, nsp*lx**2, 0:nphasx) submatrix spanning the entire
c          angular momentum basis for each unique potential
c     Lanczos algorithm: Graves-Morris,Salam, Num.Algor.21,p.213(1999)

      include '../HEADERS/dim.h'
      include 'xparam.h'
      integer  i0 (0:nphx),  lipotx(0:nphx)

      parameter (one = 1, zero = 0)
      complex coni
      parameter (coni = (0,1))

c**** array of state kets at current energy
      common /stkets/ lrstat(4, istatx), istate

      complex   tmatrx(nspx, istatx)
c     big work matrices
      complex   g0( istatx, istatx), g0t( istatx, istatx)
      logical lcalc
      dimension lcalc(0:lx)

c     Lanczos method variables
      complex xvec( istatx), wvec(istatx), x0(istatx), x1(istatx)
      complex avec(istatx), bvec(istatx)
      complex r0(istatx), r1(istatx), t0(istatx), t1(istatx)
      complex aa, dd, aw, wa, ww, e0, e1, alpha, beta, theta, q0, q1
c     return matrix containing info about each unique potential
      complex   gg(nspx*(lx+1)**2, nspx*(lx+1)**2, 0:nphasx)

c      notice that in gglu we invert (1-Gt), but here (1-tG).
c     multiply T and G0 matrices together, construct g0t = 1 - T*G0
c     notice that the signs below for g0t ARE correct since 1 is the
c     unit matrix
c     since t is tri-diagonal, this product can be computed in n^2 time
c     also fill up some work matrices for use in eigenvalue and
c     determinant calculations and elsewhere
c     cycle over dimensions of matrix g0t
      do 10 icol = 1,istatx
      do 10 irow = 1,istatx
 10   g0t(irow,icol) = 0

      do 30 icol = 1,istate
        do 20 irow = 1,istate
c         T diagonal contribution T(irow, irow)
          if ( abs( g0(irow, icol)) .gt. toler2 )
     1    g0t(irow,icol)=g0t(irow,icol) + tmatrx(1,irow) * g0(irow,icol) 

c         T off-diagonal contribution T(ist2, irow) in tmatr(2,irow)
c         T off-diagonal contribution T(irow, ist2) in tmatr(2,ist2)
          l1   = lrstat(2,irow)
          m1   = lrstat(3,irow)
          isp1 = lrstat(4,irow)
          m2 = m1+isp1
          if (nsp.eq.2 .and. m2.gt.-l1+1 .and. m2.lt.l1+2) then
c            spin-flip contribution
             ist2 = irow + (-1)**isp1
             if ( abs( g0(ist2, icol)) .gt. toler2)
     1       g0t(irow, icol) = g0t(irow, icol)
     2                   + tmatrx(nsp, ist2) * g0(ist2, icol) 
          endif
 20     continue

c       g0t(icol, icol) = g0t(icol, icol) + one
 30   continue

      do 920 ip=ipi, ipf
        ipart = nsp*(lipotx(ip)+1)**2
        do 910 is1 = 1, ipart
          is2 = is1+i0(ip)
          l1   = lrstat(2,is2)
          if (.not.lcalc(l1)) goto 910

c         start first tier with xvec=0
          istart = -1
          msord = 0
          do 40 is = 1, istate
          bvec(is) = 0
  40      xvec(is) = 0
c         rvec = bvec - A*xvec , in our case bvec(is) = delta_{is,is2}
          bvec(is2) = 1

c         RESTART here if necessary
  50      continue
          istart = istart+1

          if (istart.gt.0) then
            do 60 is = 1, istate
  60        xvec(is) = xvec(is) + x0(is) / q0
            call matvec( istatx,istate,g0t,xvec,avec,1)
            do 70 is = 1, istate
  70        bvec(is) = avec(is) - xvec(is)
            bvec(is2) = bvec(is2) + 1
          endif
          do 80 is = 1,istate
 80       r0(is) = bvec(is)
          do 90 is = 1,istate
 90       x0(is) = 0
          do 95 is = 1, istate
 95       x1(is) = bvec(is)
          call matvec( istatx,istate,g0t,bvec,r1,1)
          msord = msord + 1

c         choose wvec that del and delp close to one
          call cdot( istatx, istate, r0, r0, ww)
          call cdot( istatx, istate, r1, r1, aa)
          call cdot( istatx, istate, r0, r1, wa)
          aw = real(wa) - coni* aimag(wa)
          dd = aa*ww - aw*wa
          if (abs(dd/aa/ww) .lt.1.e-8) then
            do 96 is = 1,istate
  96        wvec(is) = r0(is) / ww
          else
            ww = ( ww - aw ) / dd
            aa = ( wa - aa) / dd
            do 97 is = 1,istate
  97        wvec(is) = r0(is) * aa + r1(is) * ww
          endif
c         update dot products to avoid round off errors
          call cdot( istatx, istate, wvec, r0, e0)
          call cdot( istatx, istate, wvec, r1, e1)
          q0 = 1
          q1 = 1

c         it seems ran out of precision for nit>150
          nitx = 10
          do 500 nit = 1, nitx
            tol = toler1 * abs(q1) /10
cc          Check convergence criteria: |r1| < tol / 10
cc          so mostly code will not exit here
            ipass = 1
            do 98 is = 1, istate
              if ( abs(real(r1(is))).gt.tol) goto 99
              if ( abs(aimag(r1(is))).gt.tol) goto 99
  98        continue
            ipass = 0
  99        continue
            if (ipass.eq.0) then
              do 100 is = 1, istate
 100          xvec(is) = xvec(is) + x1(is) / q1
              goto 700
            endif

            alpha = e1 / e0
            do 130 is = 1, istate
 130        t0(is) = r1(is) - alpha* r0(is)
            call matvec( istatx,istate,g0t,t0,t1,1)
            msord = msord + 1

            call cdot( istatx, istate, t0, t1, wa)
            call cdot( istatx, istate, t0, t0, ww)
            call cdot( istatx, istate, t1, t1, aa)
            aw = real(wa) - coni* aimag(wa)
            theta = (wa - aa) / (ww - aw)

            do 145 is = 1, istate
 145        r0(is) = t1(is) - theta * t0(is)
            dd = 1- theta
            do 150 is = 1, istate
 150        x0(is) = t0(is) + dd * (x1(is) - alpha*x0(is))
            q0 = dd * (q1 - alpha*q0)
            tol = toler1 * abs(q0)

cc          Check convergence criteria: |r0| < tol
            ipass = 1
            do 370 is = 1, istate
              if ( abs(real(r0(is))).gt.tol) goto 380
              if ( abs(aimag(r0(is))).gt.tol) goto 380
 370        continue
            ipass = 0
 380        continue
            if (ipass.eq.0) then 
              do 390 is = 1, istate
 390          xvec(is) = xvec(is) + x0(is) / q0
              goto 700
            endif

c           prepare for next iteration
            call cdot( istatx, istate, wvec, r0, e0)
            beta = e0 / e1
            do 255 is = 1, istate
 255        t0(is) = r0(is) - beta * r1(is)
            call matvec( istatx,istate,g0t,t0,avec,1)
            msord = msord + 1
            dd = beta * theta
            do 260 is = 1, istate
 260        r1(is) = avec(is) + dd * r1(is)
            call cdot( istatx, istate, wvec, r1, e1)

            dd = beta * (1-theta)
            do 270 is = 1, istate
 270        x1(is) = x0(is) - dd * x1(is) + t0(is)
            q1 = q0 - (1-theta) * beta * q1
 500      continue
c         restart since ran out of iterations
          goto 50

c         exit if tolerance has been achieved
 700      continue
c         end of GM iterations

c         at this point xvec = (1-tG)**-1 * bvec  with chosen tolerance
c         pack FMS matrix into an nsp*(lx+1)^2 x nsp*(lx+1)^2 matrix 
c         for each ipot
          do 800 is2=1,ipart
            gg( is2, is1, ip) = zero
            do 790 is = 1,istate
              gg( is2, is1, ip) = gg( is2, is1, ip) +
     1        g0( is2+i0(ip), is) * xvec(is)
 790        continue
 800      continue

 910    continue
 920  continue

      return
      end
      subroutine ggtf( nsp, i0, ipi, ipf, lipotx, g0, tmatrx, g0t, gg,
     1                 toler1, toler2, lcalc, msord)

      implicit real (a-h,o-z)
      implicit integer (i-n)
c  output
c    gg:  (nsp*lx**2, nsp*lx**2, 0:nphasx) submatrix spanning the entire
c          angular momentum basis for each unique potential
c     TFQMR: Saad, Iterative Methods for Sparse Matrices, p.225 (1996).

      include '../HEADERS/dim.h'
      include 'xparam.h'
      integer  i0 (0:nphx),  lipotx(0:nphx)

      parameter (one = 1, zero = 0)
      complex coni
      parameter (coni = (0,1))

c**** array of state kets at current energy
      common /stkets/ lrstat(4, istatx), istate

      complex   tmatrx(nspx, istatx)
c     big work matrices
      complex   g0( istatx, istatx), g0t( istatx, istatx)
      logical lcalc
      dimension lcalc(0:lx)

      complex xvec(istatx), uvec(istatx), avec(istatx), wvec(istatx)
      complex dvec(istatx), rvec(istatx), vvec(istatx)
      complex alpha, beta, aa, rho, eta
      real tau, nu, cm, err
c     return matrix containing info about each unique potential
      complex   gg(nspx*(lx+1)**2, nspx*(lx+1)**2, 0:nphasx)

c      notice that in gglu we invert (1-Gt), but here (1-tG).
c     multiply T and G0 matrices together, construct g0t = 1 - T*G0
c     notice that the signs below for g0t ARE correct since 1 is the
c     unit matrix
c     since t is tri-diagonal, this product can be computed in n^2 time
c     also fill up some work matrices for use in eigenvalue and
c     determinant calculations and elsewhere
c     cycle over dimensions of matrix g0t
      do 10 icol = 1,istatx
      do 10 irow = 1,istatx
 10   g0t(irow,icol) = 0

      do 30 icol = 1,istate
        do 20 irow = 1,istate
c         T diagonal contribution T(irow, irow)
          if ( abs( g0(irow, icol)) .gt. toler2 )
     1    g0t(irow,icol)=g0t(irow,icol) - tmatrx(1,irow) * g0(irow,icol) 

c         T off-diagonal contribution T(ist2, irow) in tmatr(2,irow)
c         T off-diagonal contribution T(irow, ist2) in tmatr(2,ist2)
          l1   = lrstat(2,irow)
          m1   = lrstat(3,irow)
          isp1 = lrstat(4,irow)
          m2 = m1+isp1
          if (nsp.eq.2 .and. m2.gt.-l1+1 .and. m2.lt.l1+2) then
c            spin-flip contribution
             ist2 = irow + (-1)**isp1
             if ( abs( g0(ist2, icol)) .gt. toler2)
     1       g0t(irow, icol) = g0t(irow, icol)
     2                   - tmatrx(nsp, ist2) * g0(ist2, icol) 
          endif
 20     continue

        g0t(icol, icol) = g0t(icol, icol) + one
 30   continue

      do 920 ip=ipi, ipf
        ipart = nsp*(lipotx(ip)+1)**2
        do 910 is1 = 1, ipart
          is2 = is1+i0(ip)
          l1   = lrstat(2,is2)
          if (.not.lcalc(l1)) goto 910

c         start first tier with xvec=0
          istart = -1
          msord = 0
          do 40 is = 1, istate
          rvec(is) = 0
          avec(is) = 0
  40      xvec(is) = 0

c         RESTART here if necessary
  50      continue
          istart = istart+1

          if (istart.gt.0) call matvec( istatx,istate,g0t,xvec,avec,1)
          do 90 is = 1,istate
 90       uvec(is) = - avec(is)
c         uvec = bvec - g0t*xvec , in our case bvec(is) = delta_{is,is2}
          uvec(is2) = uvec(is2) + 1
          call matvec( istatx,istate,g0t,uvec,avec,1)
          msord = msord + 1
          do 95 is = 1, istate
 95       wvec(is) = uvec(is)
          do 96 is = 1, istate
 96       vvec(is) = avec(is)
          do 97 is = 1, istate
 97       dvec(is) = 0
          call cdot( istatx, istate, uvec, uvec, aa)
          tau = sqrt(real(aa))
          nu = 0
          eta = 0
c         choose rvec = uvec /aa so that dot products are about 1
          do 98 is = 1, istate
 98       rvec(is) = uvec(is) / aa
          rho = 1

c         it seems ran out of precision for nit>150
          nitx = 20
          do 300 nit = 0, nitx
            if (mod(nit,2).eq.0) then
              call cdot( istatx, istate, rvec, vvec, aa)
              alpha = rho / aa
            else
              call matvec( istatx,istate,g0t,uvec,avec,1)
              msord = msord + 1
            endif

            do 115 is = 1, istate
 115        wvec(is) = wvec(is) - alpha * avec(is)

            aa = nu**2 * eta / alpha
            do 120 is = 1, istate
 120        dvec(is) = uvec(is) + aa * dvec(is)

            call cdot( istatx, istate, wvec, wvec, aa)
            nu = sqrt(real(aa)) / tau
            cm = 1 / sqrt(1+nu**2)
            tau = tau * nu * cm
            eta = cm**2 * alpha
            do 140 is = 1, istate
 140        xvec(is) = xvec(is) + eta * dvec(is)

cc          Check convergence criteria: | rvec | < tol
            err = (1.e0 + nit) / istate
            err = tau * sqrt(err) * 10
            if ( abs(err).lt.toler1) goto 700

            if (mod(nit,2) .ne.0) then
              aa = rho
              call cdot( istatx, istate, rvec, wvec, rho)
              beta = rho / aa
              do 210 is = 1, istate
 210          uvec(is) = wvec(is) + beta * uvec(is)

              do 215 is = 1, istate
 215          vvec(is) = beta * ( avec(is) + beta * vvec(is))
              call matvec( istatx,istate,g0t,uvec,avec,1)
              msord = msord + 1
              do 220 is = 1, istate
 220          vvec(is) = avec(is) + vvec(is)
            else
              do 230 is = 1, istate
 230          uvec(is) = uvec(is) - alpha * vvec(is)
            endif
 300      continue
c         restart since ran out of iterations
          goto 50

c         exit if tolerance has been achieved
 700      continue
c         end of TFQMR iterations

c         at this point xvec = (1-tG)**-1 * bvec  with chosen tolerance
c         pack FMS matrix into an nsp*(lx+1)^2 x nsp*(lx+1)^2 matrix 
c         for each ipot
          do 800 is2=1,ipart
            gg( is2, is1, ip) = zero
            do 790 is = 1,istate
              gg( is2, is1, ip) = gg( is2, is1, ip) +
     1        g0( is2+i0(ip), is) * xvec(is)
 790        continue
 800      continue

 910    continue
 920  continue

      return
      end
      subroutine fmstot( rclust, idwopt, tk, thetad, sigma2,
     1                  lmaxph, nat, iphat, ratdbl,
     2                  ipol, ispin, le2, angks, ptz,
     3                  minv, rdirec, toler1, toler2,
     4                  elnes,ipmin,ipmax,ipstep)   !KJ added this line  1-06     
c     Calculates FMS contribution to absorption
c     uses Bruce Ravel subroutine to do FMS in self-consistency loop
c     notice that it can do FMS with polarization dependence and
c     always include l-->l-1 cahnnel.
c     written by alexei ankudinov 06.1997

      implicit double precision (a-h, o-z)
      include '../HEADERS/dim.h'
      include '../HEADERS/const.h'
      include '../HEADERS/parallel.h'
      real*8 wall_commend, wall_commst
      parameter (npadx=8)
      parameter (iblock=1)

c     input
      dimension iphat(natx),  ratdbl(3,natx)
      real rat(3,natx), rclust, rdirec, toler1, toler2
      real rpart,aipart, rnrmax, thetax, temper, sig2
      integer ne, ne1, ne3,  nph, ihole
      dimension iz(0:nphx)
      integer elnes,ipmin,ipmax,ipstep  !KJ added variables 1-06      
c     polarization data
      complex*16 ptz
      dimension ptz(-1:1,-1:1)

c     work space
      complex*16 ph(nex, -ltot:ltot, nspx, 0:nphx)
      dimension lmax(nex, 0:nphx)
c     complex energy grid emg is decomposed into em and eref to have
c     the same structure in phase.bin
      complex*16 em(nex), eref(nex, nspx)
      character*6  potlbl(0:nphx)
      character*512 slog
c     fms staff
      integer lmaxph(0:nphx)
      integer map(nex)
      complex gg(nspx*(lx+1)**2, nspx*(lx+1)**2, 0:nphx)
      complex gtr(ipmin:ipmax,nex),gtrloc(ipmin:ipmax,nex) !KJ I added ip index 1-06
      complex xphase(nspx, -lx:lx, 0:nphx), ck(nspx)
      complex*16 dck, bmat
      dimension kind(8), lind(8)
      logical ltrace, lcalc
      dimension lcalc(0:lx)
!KJ commented out  1-06    dimension bmat(-lx:lx,0:1,8, -lx:lx,0:1,8)
      dimension bmat(-lx:lx,0:1,8, -lx:lx,0:1,8,ipmin:ipmax) !KJ added last index
      complex*16 rkk(nex,8,nspx)
      complex*16 bmat0(-lx:lx,0:1,8, -lx:lx,0:1,8)  !KJ new variable 1-06      
      complex*16 dum(nex)
      integer  npot, nsp, ie, iverb, lfms
      integer ip,nip !KJ 1-06 added this variable - just local index
      integer i1,i2,i3,i4,i5,i6  !KJ for stupid f77 3-06
      integer jsize !KJ added for MPI communication      
      save

      call setkap (ihole,ikap,linit)

      do 10 ie = 1,nex
      do 10 ip = ipmin,ipmax !KJ I added this line
  10  gtr(ie,ip) = 0  !KJ added ip
  
      do 12 iph = 0,nphx
      do 12 ill = -lx, lx
      do 12 isp =  1,nspx
  12   xphase(isp, ill, iph) = 0

c     need less data than rphbin.f provides, also dimensions of ph
c     array are different.
      call rdxsph (ne, ne1, ne3, nph, ihole, rnrmav, xmu, edge,
     1           ik0, em, eref, iz, potlbl, ph, rkk, lmax, lmaxp1)
      call setkap (ihole, kinit, linit)
      npot = nph
      if (rclust.le.0.0) goto 100

      write(slog,15) ' Number of energy points = ', ne
  15  format(a, i3)
      call wlog(slog)

      do 20 iat=1,nat
      do 20 j=1,3
  20  rat(j,iat) = real (ratdbl(j,iat))

c     transform to single precision
      rnrmax = real(rnrmav)
      temper = real(tk)
      thetax = real(thetad)
      sig2 = real(sigma2)
      idwopx = idwopt

      iph0 = 0
      call xprep(iph0, idwopx, nat, inclus, npot, iphat, rclust,
     1 rat, iz, rnrmax, temper, thetax, sig2, minv, rdirec)

      if (inclus.gt.1) then
cc      call fms for a cluster around central atom
        write (slog,25) inclus, iph0
  25    format 
     1     (' Doing FMS for a cluster of ',i3,' atoms around iph = ',i2)
        call wlog (slog)
        call wlog (' Please, wait (updates every 20 points) ...')

        nsp = 1
        if (abs(ispin).eq.1 ) nsp = nspx
c       if (abs(ispin).eq.1) nsp = 2
c       if (nsp.gt.nspx) call par_stop(' FMS: increase nspx and rerun.')

        ltrace = .false.
!KJ 1-06  I added the do-loop around the call to bcoef for ELNES calcul.
            do i1=1,8
	    do i2=0,1
	    do i3=-lx,lx
	    do i4=1,8
	    do i5=0,1
	    do i6=-lx,lx
	    do ip=ipmin,ipmax
	      bmat(i6,i5,i4,i3,i2,i1,ip)=dcmplx(0,0)
	    enddo
	    enddo
	    enddo
	    enddo
	    enddo
	    enddo
	    enddo
        do ip=ipmin,ipmax,ipstep
	        open(7,file='ptz.dat',form='formatted',status='unknown')
            if (elnes.eq.1) call iniptz(ptz,ip,2)  !KJ Only change ptz for ELNES !!
            call bcoef(kinit, ipol, ptz, le2, ltrace, ispin, angks, 
     1           kind, lind, bmat0)
                do i=-1,1
		write(7,'(i3,6f10.3)') ip,(ptz(i,i1),i1=-1,1)
		enddo
            do i1=1,8
	    do i2=0,1
	    do i3=-lx,lx
	    do i4=1,8
	    do i5=0,1
	    do i6=-lx,lx
	      bmat(i6,i5,i4,i3,i2,i1,ip)=bmat0(i6,i5,i4,i3,i2,i1)
	    enddo
	    enddo
	    enddo
	    enddo
	    enddo
	    enddo
!!                 bmat(:,:,:,:,:,:,ip)=bmat0(:,:,:,:,:,:)
        enddo
c !KJ end my changes to the call to bcoef        

        istart = this_process*iblock + 1
        isize = 0
        do 90 ii=istart,ne,numprocs*iblock
           do 90 ie = ii, MIN(ii+iblock-1,ne) 
            if (worker) par_type = 3
          do 30  isp = 1, nsp
            dck=sqrt(2*(em(ie)-eref(ie,isp)))
            rpart  = real( dble(dck))
            aipart = real(dimag(dck))
            ck(isp) = cmplx(rpart, aipart)
  30      continue
          do 35 ipp = 0,nph
          do 35 isp = 1, nsp
          do 35 ill = -lmaxph(ipp), lmaxph(ipp)
              rpart  = dble( ph(ie, ill, isp, ipp))
              aipart = dimag(ph(ie, ill, isp, ipp)) 
              xphase(isp, ill, ipp) = cmplx(rpart, aipart)
  35      continue
          iverb=0
          if (mod(ie,20) .eq. 0) iverb=1

          lfms = 0
          do 68  k1 = 0, lx
  68      lcalc(k1) = .false.
          do 69  k1 = 1,8
  69      if (lind(k1).ge.0.and.lind(k1).le.lx) lcalc(lind(k1)) = .true.
          call fms(lfms, nsp, ispin, inclus, npot, ck, lmaxph, xphase,
     1         ie, iverb, minv, rdirec, toler1, toler2, lcalc, gg)

          do 70  k1 = 1,8
          do 70  is1 = 1,nsp
          do 70  k2 = 1,8
          do 70  is2 = 1,nsp
            ind = ms1 + 2 * ms2
            ix1 = nsp * ( lind(k1)**2 +  lind(k1) )
            ix2 = nsp * ( lind(k2)**2 +  lind(k2) )
            ms1 = is1 - 1
            ms2 = is2 - 1
            if (lind(k2).ge.0 .and. lind(k1).ge.0) then
             do 60  m1=-lind(k1), lind(k1)
             do 60  m2=-lind(k2), lind(k2)
             
             
c !KJ I added the do-loop around the calculation of gtr for elnes
c  calculations.  1-06
    
c !KJ original call :
c !KJ              gtr(ie)=gtr(ie) + gg(ix1+nsp*m1+is1, ix2+nsp*m2+is2,0) *
c !KJ     1        bmat(m2,ms2,k2, m1,ms1,k1) * rkk(ie,k1,is1)*rkk(ie,k2,is2)
               do ip=ipmin,ipmax,ipstep
                 gtr(ip,ie) = gtr(ip,ie) +
     1                        gg(ix1+nsp*m1+is1,ix2+nsp*m2+is2,0) *
     2                        bmat(m2,ms2,k2, m1,ms1,k1,ip)
     3                      * rkk(ie,k1,is1)*rkk(ie,k2,is2)
               enddo
c !KJ end my changes
               

 60          continue

            endif
 70       continue
          
          if (worker) then
            isize = isize + 1
	    do ip=ipmin,ipmax
            gtrloc(ip,isize) = gtr(ip,ie)  !KJ added :,  1-06
	    enddo
          endif
 90     continue
        if (worker) par_type = 2
      endif

      nip=ipmax-ipmin+1 !KJ 1-06
      jsize=isize*nip !KJ

      if (numprocs .gt. 1) then
        call seconds(wall_commst)
c Collect gtr
        if (worker) then
c-- Send pointers for gtr buffer to master
          call par_send_int(jsize,1,0,this_process) !KJ isize -> jsize 1-06
c-- Send buffer
          if (jsize .ne. 0) 
     .      call par_send_cmplx(gtrloc,jsize,0,this_process) !Josh Kas isize -> jsize
        else
c Set up map of gtr to processor 
          do is = 1,numprocs-1
            istart = is*iblock + 1 
            do ii=istart,ne,numprocs*iblock
              do ie = ii, MIN(ii+iblock-1,ne) 
                 map(ie) = is
              enddo
            enddo
          enddo
          do i = 1,numprocs-1
c-- Receive pointers for gtr buffer from i
            call par_recv_int(jsize,1,i,i) !KJ isize -> jsize  1-06
c-- Receive buffer from i
            if (jsize .ne. 0) then
              call par_recv_cmplx(gtrloc,jsize,i,i) !Josh Kas isize -> jsize
              indx = 1
              do j = 1,ne
                if (map(j) .eq. i) then
c                 Josh Kas - changed gtrloc(ipmin:ipmax,indx)
                  do ip=ipmin,ipmax
                  gtr(ip,j) = gtrloc(ip,indx) !indx:indx+nip-1) !KJ used to be gtr(j)=gtrloc(indx) 
		  enddo
                  indx = indx + 1  !KJ replaced 1 by nip   1-06 - Josh Kas - nip back to 1
                endif
              enddo
            endif
          enddo
        endif
        call seconds(wall_commend)
        wall_comm = wall_comm + wall_commend - wall_commst
      endif

c     write fms.bin


 100  continue
      if (master) then
        open (unit=1, file='fms.bin', status='unknown', iostat=ios)
c        write down title line
         write(1,105) rclust*bohr
 105     format('FMS rfms=', f7.4)
         write(1, 110) ne, ne1, ne3,  nph, npadx, nip  !KJ added nip 1-06
 110     format(6(1x,i3))  !KJ changed 5 to 6
         do ip=ipmin,ipmax      !KJ I added the do loop
            do 120 ie = 1, ne
               aa = dble ( real ( gtr(ip,ie) ) ) !KJ I added ip index
               bb = dble (aimag ( gtr(ip,ie) ) ) !KJ ditto
               dum(ie) = dcmplx (aa,bb) !KJ ditto
 120        continue
            call wrpadx(1, npadx, dum(1), ne)
         enddo                  !KJ end of my loop
         close (unit=1)
      endif
      call par_barrier

      return
      end
      subroutine reafms ( mfms, rfms2, idwopt, tk, thetad, sig2g,
     1            lmaxph, nat, iphat, rat,
     2            ipol, ispin, le2, angks, ptz,
     3            minv, rdirec, toler1, toler2,
     4            elnes,ipmin,ipmax,ipstep)   !KJ added this line 1-06     

      implicit double precision (a-h, o-z)

      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'
      
      integer elnes,ipmin,ipmax,ipstep  !KJ added these variables 1-06

cc    geom.dat
        integer  nat, iatph(0:nphx), iphat(natx)
        double precision  rat(3,natx)
cc    global.dat 
c       configuration average
        integer nabs, iphabs
c       global polarization data
        integer  ipol, ispin, le2
        double precision evec(3), xivec(3), spvec(3), elpty,angks,rclabs
        complex*16 ptz(-1:1, -1:1)
cc    mod3.inp
        integer mfms, idwopt, minv
        integer lmaxph(0:nphx)
        real rfms2, rdirec, toler1, toler2
        double precision   tk, thetad, sig2g

c     Local stuff
      character*512 slog
      character*80 head(nheadx)
      dimension lhead(nheadx)

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
           nat = nat+1
           if (nat .gt. natx)  then
              write(slog,55) ' nat, natx ', nat, natx
              call wlog(slog)
  55          format(a, 2i10)
              stop 'Bad input'
           endif
           read(3,*,end=60)  idum, (rat(j,nat),j=1,3), iphat(nat), i1b
           if (iphat(nat).gt.nph) nph = iphat(nat)
           if ( iatph(iphat(nat)).eq.0) iatph(iphat(nat)) = nat
        goto 50
  60    continue
        nat = nat-1
      close(3)

cc    global.inp
      open (file='global.dat', unit=3, status='unknown',iostat=ios)
c       configuration average data
        read  (3, 10) slog
        read  (3, 65) nabs, iphabs, rclabs
  65    format ( 2i8, f13.5)
c       global polarization data
        read  (3,10)   slog
        read  (3, 70)  ipol, ispin, le2, elpty, angks
  70    format ( 3i5, 2f12.4)
        read  (3, 10) slog
        do 80 i = 1,3
          read  (3,30) evec(i), xivec(i), spvec(i)
  80    continue
        read  (3, 10) slog
        do 90 i = -1, 1
          read (3,30) a1, b1, a2, b2, a3, b3
          ptz(-1,i)= cmplx(a1,b1)
          ptz(0,i) = cmplx(a2,b2)
          ptz(1,i) = cmplx(a3,b3)
  90    continue
      close(3)
c     read mod3.inp
      open (file='mod3.inp', unit=3, status='old',iostat=ios)
        read (3,10)  slog
        read (3,20)  mfms, idwopt, minv
        read (3,10)  slog
        read (3,30)  rfms2, rdirec, toler1, toler2
        read (3,10)  slog
        read (3,30)  tk, thetad, sig2g
        read (3,10)  slog
        read (3,20)  (lmaxph(iph),iph=0,nph)
      close(3)

      
c  !KJ Next section added to read ELNES variables 1-06     
c     read eels.inp
      elnes=0
      open(file='eels.inp',unit=3,status='old',err=900)
        read(3,*,err=900,end=900) 
        read(3,20,err=900,end=900) elnes
        read(3,*,err=900,end=900)
        read(3,*,err=900,end=900)
        read(3,*,err=900,end=900)
        read(3,20,err=900,end=900) ipmin,ipstep,ipmax
      close(3)
      goto 901
900   continue
      elnes=0
901   continue
      if(elnes.eq.0) then
        ipstep=1
        ipmax=1
        ipmin=1
      endif
               
c  !KJ end my changes

c     transform to code units (bohrs and hartrees - atomic units)
      rfms2 = rfms2 / bohr
      rdirec = rdirec / bohr
      do 210 iat = 1, nat
      do 210 i = 1,3
        rat(i,iat) = rat (i, iat) / bohr
 210  continue

      return
      end
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
c     rnrmav: average norman radius in cluster (from pahse.bin)
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
      parameter(zero=0.e0)
      parameter (bohr = 0.529 177 249e0)
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
      subroutine yprep(iph0, nat, inclus, npot, iphat, rmax, rat,
     $            izx, rdirec)
c    yprep is the same as xprep for negative idwopt
c    simlifies calls in SCF and LDOS where DW factors should not enter

      implicit real (a-h,o-z)
      implicit integer (i-n)

      include '../HEADERS/dim.h'
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
      parameter(zero=0.e0)
      parameter (bohr = 0.529 177 249e0)
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
