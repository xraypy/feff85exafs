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
c     sub-program exchange
      program ffmod3
c     subroutine ffmod3

c     full multiple scattering code (inversion of big matrix)
c     written by a.ankudinov 2000 using earlier written subroutines
c     coded by b.ravel
c     modified by a.ankudinov 2001 for new matrix inversion algorithms
c     and new I/O structure

c     INPUT:  geom.dat, global.dat and mod3.inp
c     OUTPUT:  fms.bin

      implicit double precision (a-h, o-z)

c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}
c={../HEADERS/parallel.h
      integer par_type, this_process, numprocs, my_rank
      logical master, worker, parallel_run
      real*8 wall_comm, time_comm
      common /timing/ wall_comm, time_comm
      common /parallel/ numprocs, my_rank, this_process, 
     .          master, worker, parallel_run, par_type
c= ../HEADERS/parallel.h}
      real*8 wall_start, wall_end
      dimension lmaxph(0:nphx), rat(3,natx), iphat(natx)
      complex*16 ptz(-1:1, -1:1)

      real  rfms2, rdirec, toler1, toler2
      integer elnes,ipstep,ipmin,ipmax  !KJ added these variables 1-06
      character(2000) debug ! Josh Kas added for debugging.
      
      call par_begin

c     Initialize clock
      call seconds(wall_start)
      wall_comm = 0.0

c     open the log file, unit 11.  See subroutine wlog.
      if (master) then
        open (unit=11, file='log3.dat', status='unknown', iostat=ios)
        call chopen (ios, 'log3.dat', 'feff')
      else
        par_type = 2
      endif

c     read  INPUT: files geom.dat, global.dat and mod3.inp
      call reafms(mfms, rfms2, idwopt, tk, thetad, sig2g,
     1            lmaxph, nat, iphat, rat,
     2            ipol, ispin, le2, angks, ptz,
     3            minv, rdirec, toler1, toler2,
     4            elnes,ipmin,ipmax,ipstep)   !KJ added this line    1-06     

      if (mfms.eq.1 )  then
c        first do fms inside rfms2
         call fmstot(rfms2, idwopt, tk, thetad, sig2g, 
     1            lmaxph, nat, iphat, rat,
     2            ipol, ispin, le2, angks, ptz,
     3            minv, rdirec, toler1, toler2,
     4            elnes,ipmin,ipmax,ipstep)   !KJ added this line     1-06     
         call wlog(' Done with module 3: FMS.')
      endif

c     OUTPUT: gtr.bin file for the next modules

      if (master) then
         close (unit=11)
      end if

c--   Time at end of run
      call seconds(wall_end)
      if (master .and. parallel_run) then
        write (6,*) 'total time     ', wall_end - wall_start
        write (6,*) 'communicate time', wall_comm
      endif

      call par_end

c     sub-program exchange
      stop
c     return

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

c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}
c={../HEADERS/parallel.h
      integer par_type, this_process, numprocs, my_rank
      logical master, worker, parallel_run
      real*8 wall_comm, time_comm
      common /timing/ wall_comm, time_comm
      common /parallel/ numprocs, my_rank, this_process, 
     .          master, worker, parallel_run, par_type
c= ../HEADERS/parallel.h}
c={xparam.h
c nphasx MUST be the same as nphx, the maximum number of unique
c        potentials
c natxx MUST be the same as natx, the maximum number of atoms in the
c       extendeed cluster
c nexx MUST be the same as nex, the maximum number of energy points
      parameter (nphasx=nphx)
      parameter (natxx=natx)
      parameter (nexx=nex)
      parameter (istatx=(lx+1)**2*nclusx*nspx)
      parameter (nkmin=1)
c     parameter (nkmin=-9)
c= xparam.h}
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
c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}
c={xparam.h
c nphasx MUST be the same as nphx, the maximum number of unique
c        potentials
c natxx MUST be the same as natx, the maximum number of atoms in the
c       extendeed cluster
c nexx MUST be the same as nex, the maximum number of energy points
      parameter (nphasx=nphx)
      parameter (natxx=natx)
      parameter (nexx=nex)
      parameter (istatx=(lx+1)**2*nclusx*nspx)
      parameter (nkmin=1)
c     parameter (nkmin=-9)
c= xparam.h}
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

c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}
c={xparam.h
c nphasx MUST be the same as nphx, the maximum number of unique
c        potentials
c natxx MUST be the same as natx, the maximum number of atoms in the
c       extendeed cluster
c nexx MUST be the same as nex, the maximum number of energy points
      parameter (nphasx=nphx)
      parameter (natxx=natx)
      parameter (nexx=nex)
      parameter (istatx=(lx+1)**2*nclusx*nspx)
      parameter (nkmin=1)
c     parameter (nkmin=-9)
c= xparam.h}
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

c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}
c={xparam.h
c nphasx MUST be the same as nphx, the maximum number of unique
c        potentials
c natxx MUST be the same as natx, the maximum number of atoms in the
c       extendeed cluster
c nexx MUST be the same as nex, the maximum number of energy points
      parameter (nphasx=nphx)
      parameter (natxx=natx)
      parameter (nexx=nex)
      parameter (istatx=(lx+1)**2*nclusx*nspx)
      parameter (nkmin=1)
c     parameter (nkmin=-9)
c= xparam.h}
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

c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}
c={xparam.h
c nphasx MUST be the same as nphx, the maximum number of unique
c        potentials
c natxx MUST be the same as natx, the maximum number of atoms in the
c       extendeed cluster
c nexx MUST be the same as nex, the maximum number of energy points
      parameter (nphasx=nphx)
      parameter (natxx=natx)
      parameter (nexx=nex)
      parameter (istatx=(lx+1)**2*nclusx*nspx)
      parameter (nkmin=1)
c     parameter (nkmin=-9)
c= xparam.h}
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

c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}
c={xparam.h
c nphasx MUST be the same as nphx, the maximum number of unique
c        potentials
c natxx MUST be the same as natx, the maximum number of atoms in the
c       extendeed cluster
c nexx MUST be the same as nex, the maximum number of energy points
      parameter (nphasx=nphx)
      parameter (natxx=natx)
      parameter (nexx=nex)
      parameter (istatx=(lx+1)**2*nclusx*nspx)
      parameter (nkmin=1)
c     parameter (nkmin=-9)
c= xparam.h}
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

c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}
c={xparam.h
c nphasx MUST be the same as nphx, the maximum number of unique
c        potentials
c natxx MUST be the same as natx, the maximum number of atoms in the
c       extendeed cluster
c nexx MUST be the same as nex, the maximum number of energy points
      parameter (nphasx=nphx)
      parameter (natxx=natx)
      parameter (nexx=nex)
      parameter (istatx=(lx+1)**2*nclusx*nspx)
      parameter (nkmin=1)
c     parameter (nkmin=-9)
c= xparam.h}
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

c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}
c={xparam.h
c nphasx MUST be the same as nphx, the maximum number of unique
c        potentials
c natxx MUST be the same as natx, the maximum number of atoms in the
c       extendeed cluster
c nexx MUST be the same as nex, the maximum number of energy points
      parameter (nphasx=nphx)
      parameter (natxx=natx)
      parameter (nexx=nex)
      parameter (istatx=(lx+1)**2*nclusx*nspx)
      parameter (nkmin=1)
c     parameter (nkmin=-9)
c= xparam.h}
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

c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}
c={xparam.h
c nphasx MUST be the same as nphx, the maximum number of unique
c        potentials
c natxx MUST be the same as natx, the maximum number of atoms in the
c       extendeed cluster
c nexx MUST be the same as nex, the maximum number of energy points
      parameter (nphasx=nphx)
      parameter (natxx=natx)
      parameter (nexx=nex)
      parameter (istatx=(lx+1)**2*nclusx*nspx)
      parameter (nkmin=1)
c     parameter (nkmin=-9)
c= xparam.h}
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
c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}
c={../HEADERS/const.h
      parameter (pi = 3.14159 26535 89793 23846 26433d0)
      parameter (one = 1, zero = 0)
      parameter (third = one/3)
      parameter (raddeg = 180 / pi)
      complex*16 coni
      parameter (coni = (0,1))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919 158 292 677 512 811d0)

      parameter (bohr = 0.529 177 249d0, ryd  = 13.605 698d0)
      parameter (hart = 2 * ryd)
      parameter (alpinv = 137.035 989 56d0)
c     fine structure alpha
      parameter (alphfs = 1 / alpinv)
c= ../HEADERS/const.h}
c={../HEADERS/parallel.h
      integer par_type, this_process, numprocs, my_rank
      logical master, worker, parallel_run
      real*8 wall_comm, time_comm
      common /timing/ wall_comm, time_comm
      common /parallel/ numprocs, my_rank, this_process, 
     .          master, worker, parallel_run, par_type
c= ../HEADERS/parallel.h}
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

c={../HEADERS/const.h
      parameter (pi = 3.14159 26535 89793 23846 26433d0)
      parameter (one = 1, zero = 0)
      parameter (third = one/3)
      parameter (raddeg = 180 / pi)
      complex*16 coni
      parameter (coni = (0,1))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919 158 292 677 512 811d0)

      parameter (bohr = 0.529 177 249d0, ryd  = 13.605 698d0)
      parameter (hart = 2 * ryd)
      parameter (alpinv = 137.035 989 56d0)
c     fine structure alpha
      parameter (alphfs = 1 / alpinv)
c= ../HEADERS/const.h}
c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}
      
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

c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}
c={../HEADERS/parallel.h
      integer par_type, this_process, numprocs, my_rank
      logical master, worker, parallel_run
      real*8 wall_comm, time_comm
      common /timing/ wall_comm, time_comm
      common /parallel/ numprocs, my_rank, this_process, 
     .          master, worker, parallel_run, par_type
c= ../HEADERS/parallel.h}
c={xparam.h
c nphasx MUST be the same as nphx, the maximum number of unique
c        potentials
c natxx MUST be the same as natx, the maximum number of atoms in the
c       extendeed cluster
c nexx MUST be the same as nex, the maximum number of energy points
      parameter (nphasx=nphx)
      parameter (natxx=natx)
      parameter (nexx=nex)
      parameter (istatx=(lx+1)**2*nclusx*nspx)
      parameter (nkmin=1)
c     parameter (nkmin=-9)
c= xparam.h}
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

c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}
c={xparam.h
c nphasx MUST be the same as nphx, the maximum number of unique
c        potentials
c natxx MUST be the same as natx, the maximum number of atoms in the
c       extendeed cluster
c nexx MUST be the same as nex, the maximum number of energy points
      parameter (nphasx=nphx)
      parameter (natxx=natx)
      parameter (nexx=nex)
      parameter (istatx=(lx+1)**2*nclusx*nspx)
      parameter (nkmin=1)
c     parameter (nkmin=-9)
c= xparam.h}
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
c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}
c={xparam.h
c nphasx MUST be the same as nphx, the maximum number of unique
c        potentials
c natxx MUST be the same as natx, the maximum number of atoms in the
c       extendeed cluster
c nexx MUST be the same as nex, the maximum number of energy points
      parameter (nphasx=nphx)
      parameter (natxx=natx)
      parameter (nexx=nex)
      parameter (istatx=(lx+1)**2*nclusx*nspx)
      parameter (nkmin=1)
c     parameter (nkmin=-9)
c= xparam.h}
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

c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}
c={xparam.h
c nphasx MUST be the same as nphx, the maximum number of unique
c        potentials
c natxx MUST be the same as natx, the maximum number of atoms in the
c       extendeed cluster
c nexx MUST be the same as nex, the maximum number of energy points
      parameter (nphasx=nphx)
      parameter (natxx=natx)
      parameter (nexx=nex)
      parameter (istatx=(lx+1)**2*nclusx*nspx)
      parameter (nkmin=1)
c     parameter (nkmin=-9)
c= xparam.h}
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

c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}
c={xparam.h
c nphasx MUST be the same as nphx, the maximum number of unique
c        potentials
c natxx MUST be the same as natx, the maximum number of atoms in the
c       extendeed cluster
c nexx MUST be the same as nex, the maximum number of energy points
      parameter (nphasx=nphx)
      parameter (natxx=natx)
      parameter (nexx=nex)
      parameter (istatx=(lx+1)**2*nclusx*nspx)
      parameter (nkmin=1)
c     parameter (nkmin=-9)
c= xparam.h}
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
      subroutine cgetrf( M, N, A, LDA, IPIV, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX            A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  CGETRF computes an LU factorization of a general M-by-N matrix A
!  using partial pivoting with row interchanges.
!
!  The factorization has the form
!     A = P * L * U
!  where P is a permutation matrix, L is lower triangular with unit
!  diagonal elements (lower trapezoidal if m > n), and U is upper
!  triangular (upper trapezoidal if m < n).
!
!  This is the right-looking Level 3 BLAS version of the algorithm.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) COMPLEX array, dimension (LDA,N)
!          On entry, the M-by-N matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, U(i,i) is exactly zero. The
!                factorization has been completed, but the factor U
!                is exactly singular, and division by zero will occur
!                if it is used to solve a system of equations.
!
!  ==================================================================
!
!     .. Parameters ..
      COMPLEX            ONE
      PARAMETER          ( ONE = ( 1.0E+0, 0.0E+0 ) )
!     ..
!     .. Local Scalars ..
      INTEGER            I, IINFO, J, JB, NB
!     ..
!     .. External Subroutines ..
      EXTERNAL           CGEMM, CGETF2, CLASWP, CTRSM, XERBLA
!     ..
!     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CGETRF', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 )                                          &
     &   RETURN
!
!     Determine the block size for this environment.
!
      NB = ILAENV( 1, 'CGETRF', ' ', M, N, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) THEN
!
!        Use unblocked code.
!
         CALL CGETF2( M, N, A, LDA, IPIV, INFO )
      ELSE
!
!        Use blocked code.
!
         DO 20 J = 1, MIN( M, N ), NB
            JB = MIN( MIN( M, N )-J+1, NB )
!
!           Factor diagonal and subdiagonal blocks and test for exact
!           singularity.
!
            CALL CGETF2( M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO )
!
!           Adjust INFO and the pivot indices.
!
            IF( INFO.EQ.0 .AND. IINFO.GT.0 )                            &
     &         INFO = IINFO + J - 1
            DO 10 I = J, MIN( M, J+JB-1 )
               IPIV( I ) = J - 1 + IPIV( I )
   10       CONTINUE
!
!           Apply interchanges to columns 1:J-1.
!
            CALL CLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )
!
            IF( J+JB.LE.N ) THEN
!
!              Apply interchanges to columns J+JB:N.
!
               CALL CLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1,     &
     &                      IPIV, 1 )
!
!              Compute block row of U.
!
             CALL CTRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB,   &
     &                   N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ),   &
     &                     LDA )
               IF( J+JB.LE.M ) THEN
!
!                 Update trailing submatrix.
!
                CALL CGEMM( 'No transpose', 'No transpose', M-J-JB+1,   &
     &                       N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA,     &
     &                       A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ),   &
     &                       LDA )
               END IF
            END IF
   20    CONTINUE
      END IF
      RETURN
!
!     End of CGETRF
!
      END

      SUBROUTINE CGETF2( M, N, A, LDA, IPIV, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX            A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  CGETF2 computes an LU factorization of a general m-by-n matrix A
!  using partial pivoting with row interchanges.
!
!  The factorization has the form
!     A = P * L * U
!  where P is a permutation matrix, L is lower triangular with unit
!  diagonal elements (lower trapezoidal if m > n), and U is upper
!  triangular (upper trapezoidal if m < n).
!
!  This is the right-looking Level 2 BLAS version of the algorithm.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) COMPLEX array, dimension (LDA,N)
!          On entry, the m by n matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value
!          > 0: if INFO = k, U(k,k) is exactly zero. The factorization
!               has been completed, but the factor U is exactly
!               singular, and division by zero will occur if it is
!               used to solve a system of equations.
!
!  ==================================================================
!
!     .. Parameters ..
      COMPLEX            ONE, ZERO
      PARAMETER          ( ONE = ( 1.0E+0, 0.0E+0 ),                    &
     &                   ZERO = ( 0.0E+0, 0.0E+0 ) )
!     ..
!     .. Local Scalars ..
      INTEGER            J, JP
!     ..
!     .. External Functions ..
      INTEGER            ICAMAX
      EXTERNAL           ICAMAX
!     ..
!     .. External Subroutines ..
      EXTERNAL           CGERU, CSCAL, CSWAP, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CGETF2', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 )                                          &
     &   RETURN
!
      DO 10 J = 1, MIN( M, N )
!
!        Find pivot and test for singularity.
!
         JP = J - 1 + ICAMAX( M-J+1, A( J, J ), 1 )
         IPIV( J ) = JP
         IF( A( JP, J ).NE.ZERO ) THEN
!
!           Apply the interchange to columns 1:N.
!
            IF( JP.NE.J )                                               &
     &         CALL CSWAP( N, A( J, 1 ), LDA, A( JP, 1 ), LDA )
!
!           Compute elements J+1:M of J-th column.
!
            IF( J.LT.M )                                                &
     &         CALL CSCAL( M-J, ONE / A( J, J ), A( J+1, J ), 1 )
!
         ELSE IF( INFO.EQ.0 ) THEN
!
            INFO = J
         END IF
!
         IF( J.LT.MIN( M, N ) ) THEN
!
!           Update trailing submatrix.
!
            CALL CGERU( M-J, N-J, -ONE, A( J+1, J ), 1, A( J, J+1 ),    &
     &                  LDA, A( J+1, J+1 ), LDA )
         END IF
   10 CONTINUE
      RETURN
!
!     End of CGETF2
!
      END
      SUBROUTINE cgetrs( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, LDA, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX            A( LDA, * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  CGETRS solves a system of linear equations
!     A * X = B,  A**T * X = B,  or  A**H * X = B
!  with a general N-by-N matrix A using the LU factorization computed
!  by CGETRF.
!
!  Arguments
!  =========
!
!  TRANS   (input) CHARACTER*1
!          Specifies the form of the system of equations:
!          = 'N':  A * X = B     (No transpose)
!          = 'T':  A**T * X = B  (Transpose)
!          = 'C':  A**H * X = B  (Conjugate transpose)
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  A       (input) COMPLEX array, dimension (LDA,N)
!          The factors L and U from the factorization A = P*L*U
!          as computed by CGETRF.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  IPIV    (input) INTEGER array, dimension (N)
!          The pivot indices from CGETRF; for 1<=i<=N, row i of the
!          matrix was interchanged with row IPIV(i).
!
!  B       (input/output) COMPLEX array, dimension (LDB,NRHS)
!          On entry, the right hand side matrix B.
!          On exit, the solution matrix X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  ==================================================================
!
!     .. Parameters ..
      COMPLEX            ONE
      PARAMETER          ( ONE = ( 1.0E+0, 0.0E+0 ) )
!     ..
!     .. Local Scalars ..
      LOGICAL            NOTRAN
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           CLASWP, CTRSM, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.        &
     &    LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CGETRS', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 .OR. NRHS.EQ.0 )                                       &
     &   RETURN
!
      IF( NOTRAN ) THEN
!
!        Solve A * X = B.
!
!        Apply row interchanges to the right hand sides.
!
         CALL CLASWP( NRHS, B, LDB, 1, N, IPIV, 1 )
!
!        Solve L*X = B, overwriting B with X.
!
         CALL CTRSM( 'Left', 'Lower', 'No transpose', 'Unit', N, NRHS,  &
     &               ONE, A, LDA, B, LDB )
!
!        Solve U*X = B, overwriting B with X.
!
         CALL CTRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N,    &
     &               NRHS, ONE, A, LDA, B, LDB )
      ELSE
!
!        Solve A**T * X = B  or A**H * X = B.
!
!        Solve U'*X = B, overwriting B with X.
!
         CALL CTRSM( 'Left', 'Upper', TRANS, 'Non-unit', N, NRHS, ONE,  &
     &               A, LDA, B, LDB )
!
!        Solve L'*X = B, overwriting B with X.
!
         CALL CTRSM( 'Left', 'Lower', TRANS, 'Unit', N, NRHS, ONE, A,   &
     &               LDA, B, LDB )
!
!        Apply row interchanges to the solution vectors.
!
         CALL CLASWP( NRHS, B, LDB, 1, N, IPIV, -1 )
      END IF
!
      RETURN
!
!     End of CGETRS
!
      END
      SUBROUTINE XERBLA( SRNAME, INFO )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER*6        SRNAME
      INTEGER            INFO
!     ..
!
!  Purpose
!  =======
!
!  XERBLA  is an error handler for the LAPACK routines.
!  It is called by an LAPACK routine if an input parameter has an
!  invalid value.  A message is printed and execution stops.
!
!  Installers may consider modifying the STOP statement in order to
!  call system-specific exception-handling facilities.
!
!  Arguments
!  =========
!
!  SRNAME  (input) CHARACTER*6
!          The name of the routine which called XERBLA.
!
!  INFO    (input) INTEGER
!          The position of the invalid parameter in the parameter list
!          of the calling routine.
!
! ==================================================================
!
!     .. Executable Statements ..
!
      WRITE( *, FMT = 9999 )SRNAME, INFO
!
      STOP
!
 9999 FORMAT( ' ** On entry to ', A6, ' parameter number ',I2,' had ',  &
     &      'an illegal value' )
!
!     End of XERBLA
!
      END
      subroutine  cswap (n,cx,incx,cy,incy)
!
!     interchanges two vectors.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      complex cx(*),cy(*),ctemp
      integer i,incx,incy,ix,iy,n
!
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
!
!       code for unequal increments or equal increments not equal
!         to 1
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        ctemp = cx(ix)
        cx(ix) = cy(iy)
        cy(iy) = ctemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!
!       code for both increments equal to 1
   20 do 30 i = 1,n
        ctemp = cx(i)
        cx(i) = cy(i)
        cy(i) = ctemp
   30 continue
      return
      end
      subroutine  cscal(n,ca,cx,incx)
!
!     scales a vector by a constant.
!     jack dongarra, linpack,  3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      complex ca,cx(*)
      integer i,incx,n,nincx
!
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
!
!        code for increment not equal to 1
!
      nincx = n*incx
      do 10 i = 1,nincx,incx
        cx(i) = ca*cx(i)
   10 continue
      return
!
!        code for increment equal to 1
!
   20 do 30 i = 1,n
        cx(i) = ca*cx(i)
   30 continue
      return
      end
      SUBROUTINE CGERU ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
!     .. Scalar Arguments ..
      COMPLEX            ALPHA
      INTEGER            INCX, INCY, LDA, M, N
!     .. Array Arguments ..
      COMPLEX            A( LDA, * ), X( * ), Y( * )
!     ..
!
!  Purpose
!  =======
!
!  CGERU  performs the rank 1 operation
!
!     A := alpha*x*y' + A,
!
!  where alpha is a scalar, x is an m element vector, y is an n
!  element vector and A is an m by n matrix.
!
!  Parameters
!  ==========
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of the
!           matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - COMPLEX         .
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  X      - COMPLEX          array of dimension at least
!           ( 1 + ( m - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the m
!           element vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  Y      - COMPLEX          array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCY ) ).
!           Before entry, the incremented array Y must contain the n
!           element vector y.
!           Unchanged on exit.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!  A      - COMPLEX          array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients. On exit, A is
!           overwritten by the updated matrix.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as
!           declared in the calling (sub) program. LDA must be
!           at least max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
!     .. Parameters ..
      COMPLEX            ZERO
      PARAMETER        ( ZERO = ( 0.0E+0, 0.0E+0 ) )
!     .. Local Scalars ..
      COMPLEX            TEMP
      INTEGER            I, INFO, IX, J, JY, KX
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF     ( M.LT.0 )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 7
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'CGERU ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )               &
     &   RETURN
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
      IF( INCY.GT.0 )THEN
         JY = 1
      ELSE
         JY = 1 - ( N - 1 )*INCY
      END IF
      IF( INCX.EQ.1 )THEN
         DO 20, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               DO 10, I = 1, M
                  A( I, J ) = A( I, J ) + X( I )*TEMP
   10          CONTINUE
            END IF
            JY = JY + INCY
   20    CONTINUE
      ELSE
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( M - 1 )*INCX
         END IF
         DO 40, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               IX   = KX
               DO 30, I = 1, M
                  A( I, J ) = A( I, J ) + X( IX )*TEMP
                  IX        = IX        + INCX
   30          CONTINUE
            END IF
            JY = JY + INCY
   40    CONTINUE
      END IF
!
      RETURN
!
!     End of CGERU .
!
      END
      SUBROUTINE CLASWP( N, A, LDA, K1, K2, IPIV, INCX )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      INTEGER            INCX, K1, K2, LDA, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX            A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  CLASWP performs a series of row interchanges on the matrix A.
!  One row interchange is initiated for each of rows K1 through
!  K2 of A.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.
!
!  A       (input/output) COMPLEX array, dimension (LDA,N)
!          On entry, the matrix of column dimension N to which the row
!          interchanges will be applied.
!          On exit, the permuted matrix.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.
!
!  K1      (input) INTEGER
!          The first element of IPIV for which a row interchange will
!          be done.
!
!  K2      (input) INTEGER
!          The last element of IPIV for which a row interchange will
!          be done.
!
!  IPIV    (input) INTEGER array, dimension (M*abs(INCX))
!          The vector of pivot indices.  Only the elements in
!          positions K1 through K2 of IPIV are accessed.
!          IPIV(K) = L implies rows K and L are to be interchanged.
!
!  INCX    (input) INTEGER
!          The increment between successive values of IPIV.  If IPIV
!          is negative, the pivots are applied in reverse order.
!
! ==================================================================
!
!     .. Local Scalars ..
      INTEGER            I, IP, IX
!     ..
!     .. External Subroutines ..
      EXTERNAL           CSWAP
!     ..
!     .. Executable Statements ..
!
!     Interchange row I with row IPIV(I) for each of rows K1
!     through K2.
!
      IF( INCX.EQ.0 )                                                   &
     &   RETURN
      IF( INCX.GT.0 ) THEN
         IX = K1
      ELSE
         IX = 1 + ( 1-K2 )*INCX
      END IF
      IF( INCX.EQ.1 ) THEN
         DO 10 I = K1, K2
            IP = IPIV( I )
            IF( IP.NE.I )                                               &
     &         CALL CSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
   10    CONTINUE
      ELSE IF( INCX.GT.1 ) THEN
         DO 20 I = K1, K2
            IP = IPIV( IX )
            IF( IP.NE.I )                                               &
     &         CALL CSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
            IX = IX + INCX
   20    CONTINUE
      ELSE IF( INCX.LT.0 ) THEN
         DO 30 I = K2, K1, -1
            IP = IPIV( IX )
            IF( IP.NE.I )                                               &
     &         CALL CSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
            IX = IX + INCX
   30    CONTINUE
      END IF
!
      RETURN
!
!     End of CLASWP
!
      END
      SUBROUTINE CTRSM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA,A,LDA,   &
     &                   B, LDB )
!     .. Scalar Arguments ..
      CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, LDA, LDB
      COMPLEX            ALPHA
!     .. Array Arguments ..
      COMPLEX            A( LDA, * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  CTRSM  solves one of the matrix equations
!
!     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
!
!  where alpha is a scalar, X and B are m by n matrices, A is a unit,
!  or non-unit,  upper or lower triangular matrix  and  op( A )  is
!  one of
!
!     op( A ) = A   or   op( A ) = A'   or   op( A ) = conjg( A' ).
!
!  The matrix X is overwritten on B.
!
!  Parameters
!  ==========
!
!  SIDE   - CHARACTER*1.
!           On entry, SIDE specifies whether op( A ) appears on the
!           left or right of X as follows:
!
!              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
!
!              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
!
!           Unchanged on exit.
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the matrix A is an upper
!           or lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be used
!           in the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n'   op( A ) = A.
!
!              TRANSA = 'T' or 't'   op( A ) = A'.
!
!              TRANSA = 'C' or 'c'   op( A ) = conjg( A' ).
!
!           Unchanged on exit.
!
!  DIAG   - CHARACTER*1.
!           On entry, DIAG specifies whether or not A is unit
!           triangular as follows:
!
!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of B. M must
!           be at least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of B.  N
!           must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - COMPLEX         .
!           On entry,  ALPHA specifies the scalar  alpha. When alpha
!           is zero then  A is not referenced and  B need not be set
!           before entry.
!           Unchanged on exit.
!

!  A - COMPLEX array of DIMENSION ( LDA, k ), where k is m
!           when SIDE = 'L' or 'l' and is n when SIDE = 'R' or 'r'.
!           Before entry with UPLO = 'U' or 'u', the leading k by k
!           upper triangular part of the array A must contain the
!           upper triangular matrix and the strictly lower triangular
!           part of A is not referenced.
!           Before entry  with  UPLO = 'L' or 'l',  the  leading  k
!           by k lower triangular part of the array  A must contain
!           the lower triangular matrix  and the strictly upper
!           triangular part of A is not referenced.
!           Note that when  DIAG = 'U' or 'u',  the diagonal elements
!           of A  are not referenced either,  but are assumed to be
!           unity.  Unchanged on exit.

!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as
!           declared in the calling (sub) program.  When
!           SIDE = 'L' or 'l'  then LDA  must be at least
!           max( 1, m ),  when  SIDE = 'R' or 'r'
!           then LDA must be at least max( 1, n ).
!           Unchanged on exit.
!
!  B      - COMPLEX          array of DIMENSION ( LDB, n ).
!           Before entry,  the leading  m by n part of the array
!           B must contain  the  right-hand  side  matrix  B,  and
!           on exit  is overwritten by the solution matrix  X.
!
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as
!           declared in  the  calling  (sub)  program.   LDB  must
!           be  at  least max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 3 Blas routine.
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC          CONJG, MAX
!     .. Local Scalars ..
      LOGICAL            LSIDE, NOCONJ, NOUNIT, UPPER
      INTEGER            I, INFO, J, K, NROWA
      COMPLEX            TEMP
!     .. Parameters ..
      COMPLEX            ONE
      PARAMETER        ( ONE  = ( 1.0E+0, 0.0E+0 ) )
      COMPLEX            ZERO
      PARAMETER        ( ZERO = ( 0.0E+0, 0.0E+0 ) )
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      LSIDE  = LSAME( SIDE  , 'L' )
      IF( LSIDE )THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      NOCONJ = LSAME( TRANSA, 'T' )
      NOUNIT = LSAME( DIAG  , 'N' )
      UPPER  = LSAME( UPLO  , 'U' )
!
      INFO   = 0
      IF(      ( .NOT.LSIDE                ).AND.                       &
     &         ( .NOT.LSAME( SIDE  , 'R' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.UPPER                ).AND.                       &
     &         ( .NOT.LSAME( UPLO  , 'L' ) )      )THEN
         INFO = 2
      ELSE IF( ( .NOT.LSAME( TRANSA, 'N' ) ).AND.                       &
     &         ( .NOT.LSAME( TRANSA, 'T' ) ).AND.                       &
     &         ( .NOT.LSAME( TRANSA, 'C' ) )      )THEN
         INFO = 3
      ELSE IF( ( .NOT.LSAME( DIAG  , 'U' ) ).AND.                       &
     &         ( .NOT.LSAME( DIAG  , 'N' ) )      )THEN
         INFO = 4
      ELSE IF( M  .LT.0               )THEN
         INFO = 5
      ELSE IF( N  .LT.0               )THEN
         INFO = 6
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 9
      ELSE IF( LDB.LT.MAX( 1, M     ) )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'CTRSM ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( N.EQ.0 )                                                      &
     &   RETURN
!
!     And when  alpha.eq.zero.
!
      IF( ALPHA.EQ.ZERO )THEN
         DO 20, J = 1, N
            DO 10, I = 1, M
               B( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
         RETURN
      END IF
!
!     Start the operations.
!
      IF( LSIDE )THEN
         IF( LSAME( TRANSA, 'N' ) )THEN
!
!           Form  B := alpha*inv( A )*B.
!
            IF( UPPER )THEN
               DO 60, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 30, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   30                CONTINUE
                  END IF
                  DO 50, K = M, 1, -1
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT )                                    &
     &                     B( K, J ) = B( K, J )/A( K, K )
                        DO 40, I = 1, K - 1
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   40                   CONTINUE
                     END IF
   50             CONTINUE
   60          CONTINUE
            ELSE
               DO 100, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 70, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   70                CONTINUE
                  END IF
                  DO 90 K = 1, M
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT )                                    &
     &                     B( K, J ) = B( K, J )/A( K, K )
                        DO 80, I = K + 1, M
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   80                   CONTINUE
                     END IF
   90             CONTINUE
  100          CONTINUE
            END IF
         ELSE
!
!           Form  B := alpha*inv( A' )*B
!           or    B := alpha*inv( conjg( A' ) )*B.
!
            IF( UPPER )THEN
               DO 140, J = 1, N
                  DO 130, I = 1, M
                     TEMP = ALPHA*B( I, J )
                     IF( NOCONJ )THEN
                        DO 110, K = 1, I - 1
                           TEMP = TEMP - A( K, I )*B( K, J )
  110                   CONTINUE
                        IF( NOUNIT )                                    &
     &                     TEMP = TEMP/A( I, I )
                     ELSE
                        DO 120, K = 1, I - 1
                           TEMP = TEMP - CONJG( A( K, I ) )*B( K, J )
  120                   CONTINUE
                        IF( NOUNIT )                                    &
     &                     TEMP = TEMP/CONJG( A( I, I ) )
                     END IF
                     B( I, J ) = TEMP
  130             CONTINUE
  140          CONTINUE
            ELSE
               DO 180, J = 1, N
                  DO 170, I = M, 1, -1
                     TEMP = ALPHA*B( I, J )
                     IF( NOCONJ )THEN
                        DO 150, K = I + 1, M
                           TEMP = TEMP - A( K, I )*B( K, J )
  150                   CONTINUE
                        IF( NOUNIT )                                    &
     &                     TEMP = TEMP/A( I, I )
                     ELSE
                        DO 160, K = I + 1, M
                           TEMP = TEMP - CONJG( A( K, I ) )*B( K, J )
  160                   CONTINUE
                        IF( NOUNIT )                                    &
     &                     TEMP = TEMP/CONJG( A( I, I ) )
                     END IF
                     B( I, J ) = TEMP
  170             CONTINUE
  180          CONTINUE
            END IF
         END IF
      ELSE
         IF( LSAME( TRANSA, 'N' ) )THEN
!
!           Form  B := alpha*B*inv( A ).
!
            IF( UPPER )THEN
               DO 230, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 190, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  190                CONTINUE
                  END IF
                  DO 210, K = 1, J - 1
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 200, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  200                   CONTINUE
                     END IF
  210             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 220, I = 1, M
                        B( I, J ) = TEMP*B( I, J )
  220                CONTINUE
                  END IF
  230          CONTINUE
            ELSE
               DO 280, J = N, 1, -1
                  IF( ALPHA.NE.ONE )THEN
                     DO 240, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  240                CONTINUE
                  END IF
                  DO 260, K = J + 1, N
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 250, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  250                   CONTINUE
                     END IF
  260             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 270, I = 1, M
                       B( I, J ) = TEMP*B( I, J )
  270                CONTINUE
                  END IF
  280          CONTINUE
            END IF
         ELSE
!
!           Form  B := alpha*B*inv( A' )
!           or    B := alpha*B*inv( conjg( A' ) ).
!
            IF( UPPER )THEN
               DO 330, K = N, 1, -1
                  IF( NOUNIT )THEN
                     IF( NOCONJ )THEN
                        TEMP = ONE/A( K, K )
                     ELSE
                        TEMP = ONE/CONJG( A( K, K ) )
                     END IF
                     DO 290, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  290                CONTINUE
                  END IF
                  DO 310, J = 1, K - 1
                     IF( A( J, K ).NE.ZERO )THEN
                        IF( NOCONJ )THEN
                           TEMP = A( J, K )
                        ELSE
                           TEMP = CONJG( A( J, K ) )
                        END IF
                        DO 300, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  300                   CONTINUE
                     END IF
  310             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 320, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  320                CONTINUE
                  END IF
  330          CONTINUE
            ELSE
               DO 380, K = 1, N
                  IF( NOUNIT )THEN
                     IF( NOCONJ )THEN
                        TEMP = ONE/A( K, K )
                     ELSE
                        TEMP = ONE/CONJG( A( K, K ) )
                     END IF
                     DO 340, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  340                CONTINUE
                  END IF
                  DO 360, J = K + 1, N
                     IF( A( J, K ).NE.ZERO )THEN
                        IF( NOCONJ )THEN
                           TEMP = A( J, K )
                        ELSE
                           TEMP = CONJG( A( J, K ) )
                        END IF
                        DO 350, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  350                   CONTINUE
                     END IF
  360             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 370, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  370                CONTINUE
                  END IF
  380          CONTINUE
            END IF
         END IF
      END IF
!
      RETURN
!
!     End of CTRSM .
!
      END
      SUBROUTINE CGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA,B,LDB,  &
     &                   BETA, C, LDC )
!     .. Scalar Arguments ..
      CHARACTER*1        TRANSA, TRANSB
      INTEGER            M, N, K, LDA, LDB, LDC
      COMPLEX            ALPHA, BETA
!     .. Array Arguments ..
      COMPLEX            A( LDA, * ), B( LDB, * ), C( LDC, * )
!     ..
!
!  Purpose
!  =======
!
!  CGEMM  performs one of the matrix-matrix operations
!
!     C := alpha*op( A )*op( B ) + beta*C,
!
!  where  op( X ) is one of
!
!     op( X ) = X   or   op( X ) = X'   or   op( X ) = conjg( X' ),
!
!  alpha and beta are scalars, and A, B and C are matrices,
!  with op( A ) an m by k matrix,  op( B )  a  k by n matrix and
!  C an m by n matrix.
!
!  Parameters
!  ==========
!
!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be
!           used in the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n',  op( A ) = A.
!
!              TRANSA = 'T' or 't',  op( A ) = A'.
!
!              TRANSA = 'C' or 'c',  op( A ) = conjg( A' ).
!
!           Unchanged on exit.
!
!  TRANSB - CHARACTER*1.
!           On entry, TRANSB specifies the form of op( B ) to be
!           used in the matrix multiplication as follows:
!
!              TRANSB = 'N' or 'n',  op( B ) = B.
!
!              TRANSB = 'T' or 't',  op( B ) = B'.
!
!              TRANSB = 'C' or 'c',  op( B ) = conjg( B' ).
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry,  M  specifies  the number  of rows  of the
!           matrix op( A )  and of the  matrix  C.  M  must  be at
!           least  zero. Unchanged on exit.
!
!  N      - INTEGER.
!           On entry,  N  specifies the number  of columns of the
!           matrix op( B ) and the number of columns of the matrix C.
!           N must be at least zero.
!           Unchanged on exit.
!
!  K      - INTEGER.
!           On entry,  K  specifies  the number of columns of the
!           matrix op( A ) and the number of rows of the matrix
!           op( B ). K must be at least  zero.
!           Unchanged on exit.
!
!  ALPHA  - COMPLEX         .
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A   - COMPLEX          array of DIMENSION ( LDA, ka ), where ka is
!        k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
!        Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
!        part of the array  A  must contain the matrix  A,  otherwise
!        the leading  k by m  part of the array  A  must contain  the
!        matrix A.
!        Unchanged on exit.
!
!  LDA - INTEGER.
!        On entry, LDA specifies the first dimension of A as declared
!        in the calling (sub) program. When  TRANSA = 'N' or 'n' then
!        LDA must be at least  max( 1, m ), otherwise  LDA must be at
!        least  max( 1, k ).
!        Unchanged on exit.
!
!  B   - COMPLEX          array of DIMENSION ( LDB, kb ), where kb is
!        n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
!        Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
!        part of the array  B  must contain the matrix  B,  otherwise
!        the leading  n by k  part of the array  B  must contain  the
!        matrix B.
!        Unchanged on exit.
!
!  LDB - INTEGER.
!        On entry, LDB specifies the first dimension of B as declared
!        in the calling (sub) program. When  TRANSB = 'N' or 'n' then
!        LDB must be at least  max( 1, k ), otherwise  LDB must be at
!        least  max( 1, n ).
!        Unchanged on exit.
!
!  BETA   - COMPLEX         .
!        On entry,  BETA  specifies the scalar  beta.  When  BETA  is
!        supplied as zero then C need not be set on input.
!        Unchanged on exit.
!
!  C      - COMPLEX          array of DIMENSION ( LDC, n ).
!        Before entry, the leading  m by n  part of the array  C must
!        contain the matrix  C,  except when  beta  is zero, in which
!        case C need not be set on entry.
!        On exit, the array  C  is overwritten by the  m by n  matrix
!        ( alpha*op( A )*op( B ) + beta*C ).
!
!  LDC    - INTEGER.
!        On entry, LDC specifies the first dimension of C as declared
!        in  the  calling  (sub)  program.   LDC  must  be  at  least
!        max( 1, m ).
!        Unchanged on exit.
!
!
!  Level 3 Blas routine.
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC          CONJG, MAX
!     .. Local Scalars ..
      LOGICAL            CONJA, CONJB, NOTA, NOTB
      INTEGER            I, INFO, J, L, NCOLA, NROWA, NROWB
      COMPLEX            TEMP
!     .. Parameters ..
      COMPLEX            ONE
      PARAMETER        ( ONE  = ( 1.0E+0, 0.0E+0 ) )
      COMPLEX            ZERO
      PARAMETER        ( ZERO = ( 0.0E+0, 0.0E+0 ) )
!     ..
!     .. Executable Statements ..
!
!  Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
!  conjugated or transposed, set  CONJA and CONJB  as true if  A  and
!  B  respectively are to be  transposed but  not conjugated  and set
!  NROWA, NCOLA and  NROWB  as the number of rows and  columns  of  A
!  and the number of rows of  B  respectively.
!
      NOTA  = LSAME( TRANSA, 'N' )
      NOTB  = LSAME( TRANSB, 'N' )
      CONJA = LSAME( TRANSA, 'C' )
      CONJB = LSAME( TRANSB, 'C' )
      IF( NOTA )THEN
         NROWA = M
         NCOLA = K
      ELSE
         NROWA = K
         NCOLA = M
      END IF
      IF( NOTB )THEN
         NROWB = K
      ELSE
         NROWB = N
      END IF
!
!     Test the input parameters.
!
      INFO = 0
      IF(      ( .NOT.NOTA                 ).AND.                       &
     &         ( .NOT.CONJA                ).AND.                       &
     &         ( .NOT.LSAME( TRANSA, 'T' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.NOTB                 ).AND.                       &
     &         ( .NOT.CONJB                ).AND.                       &
     &         ( .NOT.LSAME( TRANSB, 'T' ) )      )THEN
         INFO = 2
      ELSE IF( M  .LT.0               )THEN
         INFO = 3
      ELSE IF( N  .LT.0               )THEN
         INFO = 4
      ELSE IF( K  .LT.0               )THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 8
      ELSE IF( LDB.LT.MAX( 1, NROWB ) )THEN
         INFO = 10
      ELSE IF( LDC.LT.MAX( 1, M     ) )THEN
         INFO = 13
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'CGEMM ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.                                  &
     &    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE )))   &
     &   RETURN
!
!     And when  alpha.eq.zero.
!
      IF( ALPHA.EQ.ZERO )THEN
         IF( BETA.EQ.ZERO )THEN
            DO 20, J = 1, N
               DO 10, I = 1, M
                  C( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               DO 30, I = 1, M
                  C( I, J ) = BETA*C( I, J )
   30          CONTINUE
   40       CONTINUE
         END IF
         RETURN
      END IF
!
!     Start the operations.
!
      IF( NOTB )THEN
         IF( NOTA )THEN
!
!           Form  C := alpha*A*B + beta*C.
!
            DO 90, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 50, I = 1, M
                     C( I, J ) = ZERO
   50             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 60, I = 1, M
                     C( I, J ) = BETA*C( I, J )
   60             CONTINUE
               END IF
               DO 80, L = 1, K
                  IF( B( L, J ).NE.ZERO )THEN
                     TEMP = ALPHA*B( L, J )
                     DO 70, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
   70                CONTINUE
                  END IF
   80          CONTINUE
   90       CONTINUE
         ELSE IF( CONJA )THEN
!
!           Form  C := alpha*conjg( A' )*B + beta*C.
!
            DO 120, J = 1, N
               DO 110, I = 1, M
                  TEMP = ZERO
                  DO 100, L = 1, K
                     TEMP = TEMP + CONJG( A( L, I ) )*B( L, J )
  100             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  110          CONTINUE
  120       CONTINUE
         ELSE
!
!           Form  C := alpha*A'*B + beta*C
!
            DO 150, J = 1, N
               DO 140, I = 1, M
                  TEMP = ZERO
                  DO 130, L = 1, K
                     TEMP = TEMP + A( L, I )*B( L, J )
  130             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  140          CONTINUE
  150       CONTINUE
         END IF
      ELSE IF( NOTA )THEN
         IF( CONJB )THEN
!
!           Form  C := alpha*A*conjg( B' ) + beta*C.
!
            DO 200, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 160, I = 1, M
                     C( I, J ) = ZERO
  160             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 170, I = 1, M
                     C( I, J ) = BETA*C( I, J )
  170             CONTINUE
               END IF
               DO 190, L = 1, K
                  IF( B( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*CONJG( B( J, L ) )
                     DO 180, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  180                CONTINUE
                  END IF
  190          CONTINUE
  200       CONTINUE
         ELSE
!
!           Form  C := alpha*A*B'          + beta*C
!
            DO 250, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 210, I = 1, M
                     C( I, J ) = ZERO
  210             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 220, I = 1, M
                     C( I, J ) = BETA*C( I, J )
  220             CONTINUE
               END IF
               DO 240, L = 1, K
                  IF( B( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*B( J, L )
                     DO 230, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  230                CONTINUE
                  END IF
  240          CONTINUE
  250       CONTINUE
         END IF
      ELSE IF( CONJA )THEN
         IF( CONJB )THEN
!
!           Form  C := alpha*conjg( A' )*conjg( B' ) + beta*C.
!
            DO 280, J = 1, N
               DO 270, I = 1, M
                  TEMP = ZERO
                  DO 260, L = 1, K
                  TEMP = TEMP + CONJG( A( L, I ) )*CONJG( B( J, L ) )
  260             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  270          CONTINUE
  280       CONTINUE
         ELSE
!
!           Form  C := alpha*conjg( A' )*B' + beta*C
!
            DO 310, J = 1, N
               DO 300, I = 1, M
                  TEMP = ZERO
                  DO 290, L = 1, K
                     TEMP = TEMP + CONJG( A( L, I ) )*B( J, L )
  290             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  300          CONTINUE
  310       CONTINUE
         END IF
      ELSE
         IF( CONJB )THEN
!
!           Form  C := alpha*A'*conjg( B' ) + beta*C
!
            DO 340, J = 1, N
               DO 330, I = 1, M
                  TEMP = ZERO
                  DO 320, L = 1, K
                     TEMP = TEMP + A( L, I )*CONJG( B( J, L ) )
  320             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  330          CONTINUE
  340       CONTINUE
         ELSE
!
!           Form  C := alpha*A'*B' + beta*C
!
            DO 370, J = 1, N
               DO 360, I = 1, M
                  TEMP = ZERO
                  DO 350, L = 1, K
                     TEMP = TEMP + A( L, I )*B( J, L )
  350             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  360          CONTINUE
  370       CONTINUE
         END IF
      END IF
!
      RETURN
!
!     End of CGEMM .
!
      END

      INTEGER          FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3,  &
     &                 N4 )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER*( * )    NAME, OPTS
      INTEGER            ISPEC, N1, N2, N3, N4
!     ..
!
!  Purpose
!  =======
!
!  ILAENV is called from the LAPACK routines to choose problem-
!  dependent parameters for the local environment.  See ISPEC for
!  a description of the parameters.
!
!  This version provides a set of parameters which should give good,
!  but not optimal, performance on many of the currently available
!  computers.  Users are encouraged to modify this subroutine to set
!  the tuning parameters for their particular machine using the option
!  and problem size information in the arguments.
!
!  This routine will not function correctly if it is converted to all
!  lower case.  Converting it to all upper case is allowed.
!
!  Arguments
!  =========
!
!  ISPEC   (input) INTEGER
!       Specifies the parameter to be returned as the value of
!       ILAENV.
!       = 1: the optimal blocksize; if this value is 1, an unblocked
!            algorithm will give the best performance.
!       = 2: the minimum block size for which the block routine
!            should be used; if the usable block size is less than
!            this value, an unblocked routine should be used.
!       = 3: the crossover point (in a block routine, for N less
!            than this value, an unblocked routine should be used)
!       = 4: the number of shifts, used in the nonsymmetric
!            eigenvalue routines
!       = 5: the minimum column dimension for blocking to be used;
!            rectangular blocks must have dimension at least k by m,
!            where k is given by ILAENV(2,...) and m by ILAENV(5,...)
!       = 6: the crossover point for the SVD (when reducing an m by n
!            matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
!            this value, a QR factorization is used first to reduce
!            the matrix to a triangular form.)
!       = 7: the number of processors
!       = 8: the crossover point for the multishift QR and QZ methods
!            for nonsymmetric eigenvalue problems.
!
!  NAME    (input) CHARACTER*(*)
!          The name of the calling subroutine, in either upper case or
!          lower case.
!
!  OPTS    (input) CHARACTER*(*)
!          The character options to the subroutine NAME, concatenated
!          into a single character string.  For example, UPLO = 'U',
!          TRANS = 'T', and DIAG = 'N' for a triangular routine would
!          be specified as OPTS = 'UTN'.
!
!  N1      (input) INTEGER
!  N2      (input) INTEGER
!  N3      (input) INTEGER
!  N4      (input) INTEGER
!          Problem dimensions for the subroutine NAME; these may not
!          all be required.
!
! (ILAENV) (output) INTEGER
!          >= 0: the value of the parameter specified by ISPEC
!          < 0:  if ILAENV = -k, the k-th argument had an illegal
!                value.
!
!  Further Details
!  ===============
!
!  The following conventions have been used when calling ILAENV from
!  the LAPACK routines:
!  1)  OPTS is a concatenation of all of the character options to
!      subroutine NAME, in the same order that they appear in the
!      argument list for NAME, even if they are not used in
!      determining the value of the parameter specified by ISPEC.
!  2)  The problem dimensions N1, N2, N3, N4 are specified in the
!      order that they appear in the argument list for NAME.  N1 is
!      used first, N2 second, and so on, and unused problem dimensions
!      are passed a value of -1.
!  3)  The parameter value returned by ILAENV is checked for validity
!      in the calling subroutine.  For example, ILAENV is used to
!      retrieve the optimal blocksize for STRTRI as follows:
!
!      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
!      IF( NB.LE.1 ) NB = MAX( 1, N )
!
!  ==================================================================
!
!     .. Local Scalars ..
      LOGICAL            CNAME, SNAME
      CHARACTER*1        C1
      CHARACTER*2        C2, C4
      CHARACTER*3        C3
      CHARACTER*6        SUBNAM
      INTEGER            I, IC, IZ, NB, NBMIN, NX
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          CHAR, ICHAR, INT, MIN, REAL
!     ..
!     .. Executable Statements ..
!
      GO TO ( 100, 100, 100, 400, 500, 600, 700, 800 ) ISPEC
!
!     Invalid value for ISPEC
!
      ILAENV = -1
      RETURN
!
  100 CONTINUE
!
!     Convert NAME to upper case if the first character is lower case.
!
      ILAENV = 1
      SUBNAM = NAME
      IC = ICHAR( SUBNAM( 1:1 ) )
      IZ = ICHAR( 'Z' )
      IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
!
!        ASCII character set
!
         IF( IC.GE.97 .AND. IC.LE.122 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 10 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.97 .AND. IC.LE.122 )                           &
     &            SUBNAM( I:I ) = CHAR( IC-32 )
   10       CONTINUE
         END IF
!
      ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
!
!        EBCDIC character set
!
         IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.                         &
     &       ( IC.GE.145 .AND. IC.LE.153 ) .OR.                         &
     &       ( IC.GE.162 .AND. IC.LE.169 ) ) THEN
            SUBNAM( 1:1 ) = CHAR( IC+64 )
            DO 20 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.                   &
     &             ( IC.GE.145 .AND. IC.LE.153 ) .OR.                   &
     &             ( IC.GE.162 .AND. IC.LE.169 ) )                      &
     &            SUBNAM( I:I ) = CHAR( IC+64 )
   20       CONTINUE
         END IF
!
      ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
!
!        Prime machines:  ASCII+128
!
         IF( IC.GE.225 .AND. IC.LE.250 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 30 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.225 .AND. IC.LE.250 )                          &
     &            SUBNAM( I:I ) = CHAR( IC-32 )
   30       CONTINUE
         END IF
      END IF
!
      C1 = SUBNAM( 1:1 )
      SNAME = C1.EQ.'S' .OR. C1.EQ.'D'
      CNAME = C1.EQ.'C' .OR. C1.EQ.'Z'
      IF( .NOT.( CNAME .OR. SNAME ) )                                   &
     &   RETURN
      C2 = SUBNAM( 2:3 )
      C3 = SUBNAM( 4:6 )
      C4 = C3( 2:3 )
!
      GO TO ( 110, 200, 300 ) ISPEC
!
  110 CONTINUE
!
!     ISPEC = 1:  block size
!
!     In these examples, separate code is provided for setting NB for
!     real and complex.  We assume that NB will take the same value in
!     single or double precision.
!
      NB = 1
!
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.    &
     &            C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'PO' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NB = 1
         ELSE IF( SNAME .AND. C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            NB = 64
         ELSE IF( C3.EQ.'TRD' ) THEN
            NB = 1
         ELSE IF( C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.         &
     &          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.         &
     &          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.         &
     &          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.         &
     &          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.         &
     &          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.         &
     &          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.         &
     &          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.         &
     &          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( C2.EQ.'GB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'PB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'TR' ) THEN
         IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'LA' ) THEN
         IF( C3.EQ.'UUM' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'ST' ) THEN
         IF( C3.EQ.'EBZ' ) THEN
            NB = 1
         END IF
      END IF
      ILAENV = NB
      RETURN
!
  200 CONTINUE
!
!     ISPEC = 2:  minimum block size
!
      NBMIN = 2
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.         &
     &       C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 8
            ELSE
               NBMIN = 8
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.         &
     &          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.         &
     &          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.         &
     &          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.         &
     &          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.         &
     &          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.         &
     &          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.         &
     &          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.         &
     &          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      END IF
      ILAENV = NBMIN
      RETURN
!
  300 CONTINUE
!
!     ISPEC = 3:  crossover point
!
      NX = 0
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.         &
     &       C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NX = 1
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NX = 1
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.         &
     &          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.         &
     &          C4.EQ.'BR' ) THEN
               NX = 128
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.         &
     &          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.         &
     &          C4.EQ.'BR' ) THEN
               NX = 128
            END IF
         END IF
      END IF
      ILAENV = NX
      RETURN
!
  400 CONTINUE
!
!     ISPEC = 4:  number of shifts (used by xHSEQR)
!
      ILAENV = 6
      RETURN
!
  500 CONTINUE
!
!     ISPEC = 5:  minimum column dimension (not used)
!
      ILAENV = 2
      RETURN
!
  600 CONTINUE
!
!     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
!
      ILAENV = INT( REAL( MIN( N1, N2 ) )*1.6E0 )
      RETURN
!
  700 CONTINUE
!
!     ISPEC = 7:  number of processors (not used)
!
      ILAENV = 1
      RETURN
!
  800 CONTINUE
!
!     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
!
      ILAENV = 50
      RETURN
!
!     End of ILAENV
!
      END

      integer function icamax(n,cx,incx)
!
!     finds the index of element having max. absolute value.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      complex cx(*)
      real smax
      integer i,incx,ix,n
      complex zdum
      real cabs1
      cabs1(zdum) = abs(real(zdum)) + abs(aimag(zdum))
!
      icamax = 0
      if( n.lt.1 .or. incx.le.0 ) return
      icamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
!
!        code for increment not equal to 1
!
      ix = 1
      smax = cabs1(cx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(cabs1(cx(ix)).le.smax) go to 5
         icamax = i
         smax = cabs1(cx(ix))
    5    ix = ix + incx
   10 continue
      return
!
!        code for increment equal to 1
!
   20 smax = cabs1(cx(1))
      do 30 i = 2,n
         if(cabs1(cx(i)).le.smax) go to 30
         icamax = i
         smax = cabs1(cx(i))
   30 continue
      return
      end

      LOGICAL          FUNCTION LSAME( CA, CB )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     January 31, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          CA, CB
!     ..
!
!  Purpose
!  =======
!
!  LSAME returns .TRUE. if CA is the same letter as CB regardless of
!  case.
!
!  Arguments
!  =========
!
!  CA      (input) CHARACTER*1
!  CB      (input) CHARACTER*1
!          CA and CB specify the single characters to be compared.
!
! ==================================================================
!
!     .. Intrinsic Functions ..
      INTRINSIC          ICHAR
!     ..
!     .. Local Scalars ..
      INTEGER            INTA, INTB, ZCODE
!     ..
!     .. Executable Statements ..
!
!     Test if the characters are equal
!
      LSAME = CA.EQ.CB
      IF( LSAME )                                                       &
     &   RETURN
!
!     Now test for equivalence if both characters are alphabetic.
!
      ZCODE = ICHAR( 'Z' )
!
!     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
!     machines, on which ICHAR returns a value with bit 8 set.
!     ICHAR('A') on Prime machines returns 193 which is the same as
!     ICHAR('A') on an EBCDIC machine.
!
      INTA = ICHAR( CA )
      INTB = ICHAR( CB )
!
      IF( ZCODE.EQ.90 .OR. ZCODE.EQ.122 ) THEN
!
!        ASCII is assumed - ZCODE is the ASCII code of either lower or
!        upper case 'Z'.
!
         IF( INTA.GE.97 .AND. INTA.LE.122 ) INTA = INTA - 32
         IF( INTB.GE.97 .AND. INTB.LE.122 ) INTB = INTB - 32
!
      ELSE IF( ZCODE.EQ.233 .OR. ZCODE.EQ.169 ) THEN
!
!        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower
!        or upper case 'Z'.
!
         IF( INTA.GE.129 .AND. INTA.LE.137 .OR.                         &
     &       INTA.GE.145 .AND. INTA.LE.153 .OR.                         &
     &       INTA.GE.162 .AND. INTA.LE.169 ) INTA = INTA + 64
         IF( INTB.GE.129 .AND. INTB.LE.137 .OR.                         &
     &       INTB.GE.145 .AND. INTB.LE.153 .OR.                         &
     &       INTB.GE.162 .AND. INTB.LE.169 ) INTB = INTB + 64
!
      ELSE IF( ZCODE.EQ.218 .OR. ZCODE.EQ.250 ) THEN
!
!        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
!        plus 128 of either lower or upper case 'Z'.
!
         IF( INTA.GE.225 .AND. INTA.LE.250 ) INTA = INTA - 32
         IF( INTB.GE.225 .AND. INTB.LE.250 ) INTB = INTB - 32
      END IF
      LSAME = INTA.EQ.INTB
!
!     RETURN
!
!     End of LSAME
!
      END
!       SUBROUTINE CTRTRI( UPLO, DIAG, N, A, LDA, INFO )
! *
! *  -- LAPACK routine (version 2.0) --
! *     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
! *     Courant Institute, Argonne National Lab, and Rice University
! *     September 30, 1994
! *
! *     .. Scalar Arguments ..
!       CHARACTER          DIAG, UPLO
!       INTEGER            INFO, LDA, N
! *     ..
! *     .. Array Arguments ..
!       COMPLEX            A( LDA, * )
! *     ..
! *
! *  Purpose
! *  =======
! *
! *  CTRTRI computes the inverse of a complex upper or lower
! *  triangular matrix A.
! *
! *  This is the Level 3 BLAS version of the algorithm.
! *
! *  Arguments
! *  =========
! *
! *  UPLO    (input) CHARACTER*1
! *          = 'U':  A is upper triangular;
! *          = 'L':  A is lower triangular.
! *
! *  DIAG    (input) CHARACTER*1
! *          = 'N':  A is non-unit triangular;
! *          = 'U':  A is unit triangular.
! *
! *  N       (input) INTEGER
! *          The order of the matrix A.  N >= 0.
! *
! *  A       (input/output) COMPLEX array, dimension (LDA,N)
! *       On entry, the triangular matrix A.  If UPLO = 'U', the
! *       leading N-by-N upper triangular part of the array A
! *       contains the upper triangular matrix, and the strictly lower
! *       triangular part of A is not referenced.  If UPLO = 'L', the
! *       leading N-by-N lower triangular part of the array A contains
! *       the lower triangular matrix, and the strictly upper
! *       triangular part of A is not referenced.  If DIAG = 'U', the
! *       diagonal elements of A are also not referenced and are
! *       assumed to be 1.
! *       On exit, the (triangular) inverse of the original matrix, in
! *       the same storage format.
! *
! *  LDA     (input) INTEGER
! *          The leading dimension of the array A.  LDA >= max(1,N).
! *
! *  INFO    (output) INTEGER
! *       = 0: successful exit
! *       < 0: if INFO = -i, the i-th argument had an illegal value
! *       > 0: if INFO = i, A(i,i) is exactly zero.  The triangular
! *            matrix is singular and its inverse can not be computed.
! *
! *  ================================================================
! *
! *     .. Parameters ..
!       COMPLEX            ONE, ZERO
!       PARAMETER          ( ONE = ( 1.0E+0, 0.0E+0 ),
!      $                   ZERO = ( 0.0E+0, 0.0E+0 ) )
! *     ..
! *     .. Local Scalars ..
!       LOGICAL            NOUNIT, UPPER
!       INTEGER            J, JB, NB, NN
! *     ..
! *     .. External Functions ..
!       LOGICAL            LSAME
!       INTEGER            ILAENV
!       EXTERNAL           LSAME, ILAENV
! *     ..
! *     .. External Subroutines ..
!       EXTERNAL           CTRMM, CTRSM, CTRTI2, XERBLA
! *     ..
! *     .. Intrinsic Functions ..
!       INTRINSIC          MAX, MIN
! *     ..
! *     .. Executable Statements ..
! *
! *     Test the input parameters.
! *
!       INFO = 0
!       UPPER = LSAME( UPLO, 'U' )
!       NOUNIT = LSAME( DIAG, 'N' )
!       IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
!          INFO = -1
!       ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
!          INFO = -2
!       ELSE IF( N.LT.0 ) THEN
!          INFO = -3
!       ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
!          INFO = -5
!       END IF
!       IF( INFO.NE.0 ) THEN
!          CALL XERBLA( 'CTRTRI', -INFO )
!          RETURN
!       END IF
! *
! *     Quick return if possible
! *
!       IF( N.EQ.0 )
!      $   RETURN
! *
! *     Check for singularity if non-unit.
! *
!       IF( NOUNIT ) THEN
!          DO 10 INFO = 1, N
!             IF( A( INFO, INFO ).EQ.ZERO )
!      $         RETURN
!    10    CONTINUE
!          INFO = 0
!       END IF
! *
! *     Determine the block size for this environment.
! *
!       NB = ILAENV( 1, 'CTRTRI', UPLO // DIAG, N, -1, -1, -1 )
!       IF( NB.LE.1 .OR. NB.GE.N ) THEN
! *
! *        Use unblocked code
! *
!          CALL CTRTI2( UPLO, DIAG, N, A, LDA, INFO )
!       ELSE
! *
! *        Use blocked code
! *
!          IF( UPPER ) THEN
! *
! *           Compute inverse of upper triangular matrix
! *
!             DO 20 J = 1, N, NB
!                JB = MIN( NB, N-J+1 )
! *
! *              Compute rows 1:j-1 of current block column
! *
!             CALL CTRMM( 'Left', 'Upper', 'No transpose', DIAG, J-1,
!      $                  JB, ONE, A, LDA, A( 1, J ), LDA )
!             CALL CTRSM( 'Right', 'Upper', 'No transpose', DIAG, J-1,
!      $                  JB, -ONE, A( J, J ), LDA, A( 1, J ), LDA )
! *
! *              Compute inverse of current diagonal block
! *
!             CALL CTRTI2( 'Upper', DIAG, JB, A( J, J ), LDA, INFO )
!    20       CONTINUE
!          ELSE
! *
! *           Compute inverse of lower triangular matrix
! *
!             NN = ( ( N-1 ) / NB )*NB + 1
!             DO 30 J = NN, 1, -NB
!                JB = MIN( NB, N-J+1 )
!                IF( J+JB.LE.N ) THEN
! *
! *                 Compute rows j+jb:n of current block column
! *
!                CALL CTRMM( 'Left', 'Lower', 'No transpose', DIAG,
!      $                     N-J-JB+1, JB, ONE, A( J+JB, J+JB ), LDA,
!      $                     A( J+JB, J ), LDA )
!                CALL CTRSM( 'Right', 'Lower', 'No transpose', DIAG,
!      $                     N-J-JB+1, JB, -ONE, A( J, J ), LDA,
!      $                     A( J+JB, J ), LDA )
!                END IF
! *
! *              Compute inverse of current diagonal block
! *
!             CALL CTRTI2( 'Lower', DIAG, JB, A( J, J ), LDA, INFO )
!    30       CONTINUE
!          END IF
!       END IF
! *
!       RETURN
! *
! *     End of CTRTRI
! *
!       END
      SUBROUTINE CTRTI2( UPLO, DIAG, N, A, LDA, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          DIAG, UPLO
      INTEGER            INFO, LDA, N
!     ..
!     .. Array Arguments ..
      COMPLEX            A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  CTRTI2 computes the inverse of a complex upper or lower triangular
!  matrix.
!
!  This is the Level 2 BLAS version of the algorithm.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          Specifies whether the matrix A is upper or lower
!          triangular.
!          = 'U':  Upper triangular
!          = 'L':  Lower triangular
!
!  DIAG    (input) CHARACTER*1
!          Specifies whether or not the matrix A is unit triangular.
!          = 'N':  Non-unit triangular
!          = 'U':  Unit triangular
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input/output) COMPLEX array, dimension (LDA,N)
!       On entry, the triangular matrix A.  If UPLO = 'U', the
!       leading n by n upper triangular part of the array A contains
!       the upper triangular matrix, and the strictly lower
!       triangular part of A is not referenced.  If UPLO = 'L', the
!       leading n by n lower triangular part of the array A contains
!       the lower triangular matrix, and the strictly upper
!       triangular part of A is not referenced.  If DIAG = 'U', the
!       diagonal elements of A are also not referenced and are
!       assumed to be 1.
!
!       On exit, the (triangular) inverse of the original matrix, in
!       the same storage format.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value
!
!  ==================================================================
!
!     .. Parameters ..
      COMPLEX            ONE
      PARAMETER          ( ONE = ( 1.0E+0, 0.0E+0 ) )
!     ..
!     .. Local Scalars ..
      LOGICAL            NOUNIT, UPPER
      INTEGER            J
      COMPLEX            AJJ
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           CSCAL, CTRMV, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      NOUNIT = LSAME( DIAG, 'N' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CTRTI2', -INFO )
         RETURN
      END IF
!
      IF( UPPER ) THEN
!
!        Compute inverse of upper triangular matrix.
!
         DO 10 J = 1, N
            IF( NOUNIT ) THEN
               A( J, J ) = ONE / A( J, J )
               AJJ = -A( J, J )
            ELSE
               AJJ = -ONE
            END IF
!
!           Compute elements 1:j-1 of j-th column.
!
            CALL CTRMV( 'Upper', 'No transpose', DIAG, J-1, A, LDA,     &
     &                  A( 1, J ), 1 )
            CALL CSCAL( J-1, AJJ, A( 1, J ), 1 )
   10    CONTINUE
      ELSE
!
!        Compute inverse of lower triangular matrix.
!
         DO 20 J = N, 1, -1
            IF( NOUNIT ) THEN
               A( J, J ) = ONE / A( J, J )
               AJJ = -A( J, J )
            ELSE
               AJJ = -ONE
            END IF
            IF( J.LT.N ) THEN
!
!              Compute elements j+1:n of j-th column.
!
               CALL CTRMV( 'Lower', 'No transpose', DIAG, N-J,          &
     &                     A( J+1, J+1 ), LDA, A( J+1, J ), 1 )
               CALL CSCAL( N-J, AJJ, A( J+1, J ), 1 )
            END IF
   20    CONTINUE
      END IF
!
      RETURN
!
!     End of CTRTI2
!
      END
      SUBROUTINE CTRMV ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
!     .. Scalar Arguments ..
      INTEGER            INCX, LDA, N
      CHARACTER*1        DIAG, TRANS, UPLO
!     .. Array Arguments ..
      COMPLEX            A( LDA, * ), X( * )
!     ..
!
!  Purpose
!  =======
!
!  CTRMV  performs one of the matrix-vector operations
!
!     x := A*x,   or   x := A'*x,   or   x := conjg( A' )*x,
!
!  where x is an n element vector and  A is an n by n unit, or
!  non-unit, upper or lower triangular matrix.
!
!  Parameters
!  ==========
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the matrix is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANS  - CHARACTER*1.
!           On entry, TRANS specifies the operation to be performed as
!           follows:
!
!              TRANS = 'N' or 'n'   x := A*x.
!
!              TRANS = 'T' or 't'   x := A'*x.
!
!              TRANS = 'C' or 'c'   x := conjg( A' )*x.
!
!           Unchanged on exit.
!
!  DIAG   - CHARACTER*1.
!           On entry, DIAG specifies whether or not A is unit
!           triangular as follows:
!
!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.
!
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  A      - COMPLEX          array of DIMENSION ( LDA, n ).
!        Before entry with  UPLO = 'U' or 'u', the leading n by n
!        upper triangular part of the array A must contain the upper
!        triangular matrix and the strictly lower triangular part of
!        A is not referenced.
!        Before entry with UPLO = 'L' or 'l', the leading n by n
!        lower triangular part of the array A must contain the lower
!        triangular matrix and the strictly upper triangular part of
!        A is not referenced.
!        Note that when  DIAG = 'U' or 'u', the diagonal elements of
!        A are not referenced either, but are assumed to be unity.
!        Unchanged on exit.
!
!  LDA    - INTEGER.
!        On entry, LDA specifies the first dimension of A as declared
!        in the calling (sub) program. LDA must be at least
!        max( 1, n ).
!        Unchanged on exit.
!
!  X      - COMPLEX          array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the n
!           element vector x. On exit, X is overwritten with the
!           tranformed vector x.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
!     .. Parameters ..
      COMPLEX            ZERO
      PARAMETER        ( ZERO = ( 0.0E+0, 0.0E+0 ) )
!     .. Local Scalars ..
      COMPLEX            TEMP
      INTEGER            I, INFO, IX, J, JX, KX
      LOGICAL            NOCONJ, NOUNIT
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC          CONJG, MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF     ( .NOT.LSAME( UPLO , 'U' ).AND.                            &
     &         .NOT.LSAME( UPLO , 'L' )      )THEN
         INFO = 1
      ELSE IF( .NOT.LSAME( TRANS, 'N' ).AND.                            &
     &         .NOT.LSAME( TRANS, 'T' ).AND.                            &
     &         .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 2
      ELSE IF( .NOT.LSAME( DIAG , 'U' ).AND.                            &
     &         .NOT.LSAME( DIAG , 'N' )      )THEN
         INFO = 3
      ELSE IF( N.LT.0 )THEN
         INFO = 4
      ELSE IF( LDA.LT.MAX( 1, N ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'CTRMV ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( N.EQ.0 )                                                      &
     &   RETURN
!
      NOCONJ = LSAME( TRANS, 'T' )
      NOUNIT = LSAME( DIAG , 'N' )
!
!     Set up the start point in X if the increment is not unity. This
!     will be  ( N - 1 )*INCX  too small for descending loops.
!
      IF( INCX.LE.0 )THEN
         KX = 1 - ( N - 1 )*INCX
      ELSE IF( INCX.NE.1 )THEN
         KX = 1
      END IF
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
      IF( LSAME( TRANS, 'N' ) )THEN
!
!        Form  x := A*x.
!
         IF( LSAME( UPLO, 'U' ) )THEN
            IF( INCX.EQ.1 )THEN
               DO 20, J = 1, N
                  IF( X( J ).NE.ZERO )THEN
                     TEMP = X( J )
                     DO 10, I = 1, J - 1
                        X( I ) = X( I ) + TEMP*A( I, J )
   10                CONTINUE
                     IF( NOUNIT )                                       &
     &                  X( J ) = X( J )*A( J, J )
                  END IF
   20          CONTINUE
            ELSE
               JX = KX
               DO 40, J = 1, N
                  IF( X( JX ).NE.ZERO )THEN
                     TEMP = X( JX )
                     IX   = KX
                     DO 30, I = 1, J - 1
                        X( IX ) = X( IX ) + TEMP*A( I, J )
                        IX      = IX      + INCX
   30                CONTINUE
                     IF( NOUNIT )                                       &
     &                  X( JX ) = X( JX )*A( J, J )
                  END IF
                  JX = JX + INCX
   40          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 60, J = N, 1, -1
                  IF( X( J ).NE.ZERO )THEN
                     TEMP = X( J )
                     DO 50, I = N, J + 1, -1
                        X( I ) = X( I ) + TEMP*A( I, J )
   50                CONTINUE
                     IF( NOUNIT )                                       &
     &                  X( J ) = X( J )*A( J, J )
                  END IF
   60          CONTINUE
            ELSE
               KX = KX + ( N - 1 )*INCX
               JX = KX
               DO 80, J = N, 1, -1
                  IF( X( JX ).NE.ZERO )THEN
                     TEMP = X( JX )
                     IX   = KX
                     DO 70, I = N, J + 1, -1
                        X( IX ) = X( IX ) + TEMP*A( I, J )
                        IX      = IX      - INCX
   70                CONTINUE
                     IF( NOUNIT )                                       &
     &                  X( JX ) = X( JX )*A( J, J )
                  END IF
                  JX = JX - INCX
   80          CONTINUE
            END IF
         END IF
      ELSE
!
!        Form  x := A'*x  or  x := conjg( A' )*x.
!
         IF( LSAME( UPLO, 'U' ) )THEN
            IF( INCX.EQ.1 )THEN
               DO 110, J = N, 1, -1
                  TEMP = X( J )
                  IF( NOCONJ )THEN
                     IF( NOUNIT )                                       &
     &                  TEMP = TEMP*A( J, J )
                     DO 90, I = J - 1, 1, -1
                        TEMP = TEMP + A( I, J )*X( I )
   90                CONTINUE
                  ELSE
                     IF( NOUNIT )                                       &
     &                  TEMP = TEMP*CONJG( A( J, J ) )
                     DO 100, I = J - 1, 1, -1
                        TEMP = TEMP + CONJG( A( I, J ) )*X( I )
  100                CONTINUE
                  END IF
                  X( J ) = TEMP
  110          CONTINUE
            ELSE
               JX = KX + ( N - 1 )*INCX
               DO 140, J = N, 1, -1
                  TEMP = X( JX )
                  IX   = JX
                  IF( NOCONJ )THEN
                     IF( NOUNIT )                                       &
     &                  TEMP = TEMP*A( J, J )
                     DO 120, I = J - 1, 1, -1
                        IX   = IX   - INCX
                        TEMP = TEMP + A( I, J )*X( IX )
  120                CONTINUE
                  ELSE
                     IF( NOUNIT )                                       &
     &                  TEMP = TEMP*CONJG( A( J, J ) )
                     DO 130, I = J - 1, 1, -1
                        IX   = IX   - INCX
                        TEMP = TEMP + CONJG( A( I, J ) )*X( IX )
  130                CONTINUE
                  END IF
                  X( JX ) = TEMP
                  JX      = JX   - INCX
  140          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 170, J = 1, N
                  TEMP = X( J )
                  IF( NOCONJ )THEN
                     IF( NOUNIT )                                       &
     &                  TEMP = TEMP*A( J, J )
                     DO 150, I = J + 1, N
                        TEMP = TEMP + A( I, J )*X( I )
  150                CONTINUE
                  ELSE
                     IF( NOUNIT )                                       &
     &                  TEMP = TEMP*CONJG( A( J, J ) )
                     DO 160, I = J + 1, N
                        TEMP = TEMP + CONJG( A( I, J ) )*X( I )
  160                CONTINUE
                  END IF
                  X( J ) = TEMP
  170          CONTINUE
            ELSE
               JX = KX
               DO 200, J = 1, N
                  TEMP = X( JX )
                  IX   = JX
                  IF( NOCONJ )THEN
                     IF( NOUNIT )                                       &
     &                  TEMP = TEMP*A( J, J )
                     DO 180, I = J + 1, N
                        IX   = IX   + INCX
                        TEMP = TEMP + A( I, J )*X( IX )
  180                CONTINUE
                  ELSE
                     IF( NOUNIT )                                       &
     &                  TEMP = TEMP*CONJG( A( J, J ) )
                     DO 190, I = J + 1, N
                        IX   = IX   + INCX
                        TEMP = TEMP + CONJG( A( I, J ) )*X( IX )
  190                CONTINUE
                  END IF
                  X( JX ) = TEMP
                  JX      = JX   + INCX
  200          CONTINUE
            END IF
         END IF
      END IF
!
      RETURN
!
!     End of CTRMV .
!
      END
      subroutine caxpy(n,ca,cx,incx,cy,incy)
!
!     constant times a vector plus a vector.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      complex cx(*),cy(*),ca
      integer i,incx,incy,ix,iy,n
!
      if(n.le.0)return
      if (abs(real(ca)) + abs(aimag(ca)) .eq. 0.0 ) return
      if(incx.eq.1.and.incy.eq.1)go to 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        cy(iy) = cy(iy) + ca*cx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!
!        code for both increments equal to 1
!
   20 do 30 i = 1,n
        cy(i) = cy(i) + ca*cx(i)
   30 continue
      return
      end
      subroutine  ccopy(n,cx,incx,cy,incy)
!
!     copies a vector, x, to a vector, y.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      complex cx(*),cy(*)
      integer i,incx,incy,ix,iy,n
!
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        cy(iy) = cx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!
!        code for both increments equal to 1
!
   20 do 30 i = 1,n
        cy(i) = cx(i)
   30 continue
      return
      end

      SUBROUTINE CGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX,           &
     &                   BETA, Y, INCY )
!     .. Scalar Arguments ..
      COMPLEX            ALPHA, BETA
      INTEGER            INCX, INCY, LDA, M, N
      CHARACTER*1        TRANS
!     .. Array Arguments ..
      COMPLEX            A( LDA, * ), X( * ), Y( * )
!     ..
!
!  Purpose
!  =======
!
!  CGEMV  performs one of the matrix-vector operations
!
!     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,   or
!
!     y := alpha*conjg( A' )*x + beta*y,
!
!  where alpha and beta are scalars, x and y are vectors and A is an
!  m by n matrix.
!
!  Parameters
!  ==========
!
!  TRANS  - CHARACTER*1.
!           On entry, TRANS specifies the operation to be performed as
!           follows:
!
!              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
!
!              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
!
!              TRANS = 'C' or 'c'   y := alpha*conjg( A' )*x + beta*y.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - COMPLEX         .
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - COMPLEX          array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit.
!
!  X      - COMPLEX          array of DIMENSION at least
!           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
!           Before entry, the incremented array X must contain the
!           vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  BETA   - COMPLEX         .
!           On entry, BETA specifies the scalar beta. When BETA is
!           supplied as zero then Y need not be set on input.
!           Unchanged on exit.
!
!  Y      - COMPLEX          array of DIMENSION at least
!           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
!           Before entry with BETA non-zero, the incremented array Y
!           must contain the vector y. On exit, Y is overwritten by the
!           updated vector y.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
!     .. Parameters ..
      COMPLEX            ONE
      PARAMETER        ( ONE  = ( 1.0E+0, 0.0E+0 ) )
      COMPLEX            ZERO
      PARAMETER        ( ZERO = ( 0.0E+0, 0.0E+0 ) )
!     .. Local Scalars ..
      COMPLEX            TEMP
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY, LENX, LENY
      LOGICAL            NOCONJ
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC          CONJG, MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF     ( .NOT.LSAME( TRANS, 'N' ).AND.                            &
     &         .NOT.LSAME( TRANS, 'T' ).AND.                            &
     &         .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 1
      ELSE IF( M.LT.0 )THEN
         INFO = 2
      ELSE IF( N.LT.0 )THEN
         INFO = 3
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'CGEMV ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.                                  &
     &    ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )                   &
     &   RETURN
!
      NOCONJ = LSAME( TRANS, 'T' )
!
!     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
!     up the start points in  X  and  Y.
!
      IF( LSAME( TRANS, 'N' ) )THEN
         LENX = N
         LENY = M
      ELSE
         LENX = M
         LENY = N
      END IF
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( LENX - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( LENY - 1 )*INCY
      END IF
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
!     First form  y := beta*y.
!
      IF( BETA.NE.ONE )THEN
         IF( INCY.EQ.1 )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, LENY
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, LENY
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, LENY
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, LENY
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO )                                               &
     &   RETURN
      IF( LSAME( TRANS, 'N' ) )THEN
!
!        Form  y := alpha*A*x + y.
!
         JX = KX
         IF( INCY.EQ.1 )THEN
            DO 60, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  DO 50, I = 1, M
                     Y( I ) = Y( I ) + TEMP*A( I, J )
   50             CONTINUE
               END IF
               JX = JX + INCX
   60       CONTINUE
         ELSE
            DO 80, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IY   = KY
                  DO 70, I = 1, M
                     Y( IY ) = Y( IY ) + TEMP*A( I, J )
                     IY      = IY      + INCY
   70             CONTINUE
               END IF
               JX = JX + INCX
   80       CONTINUE
         END IF
      ELSE
!
!        Form  y := alpha*A'*x + y  or  y := alpha*conjg( A' )*x + y.
!
         JY = KY
         IF( INCX.EQ.1 )THEN
            DO 110, J = 1, N
               TEMP = ZERO
               IF( NOCONJ )THEN
                  DO 90, I = 1, M
                     TEMP = TEMP + A( I, J )*X( I )
   90             CONTINUE
               ELSE
                  DO 100, I = 1, M
                     TEMP = TEMP + CONJG( A( I, J ) )*X( I )
  100             CONTINUE
               END IF
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  110       CONTINUE
         ELSE
            DO 140, J = 1, N
               TEMP = ZERO
               IX   = KX
               IF( NOCONJ )THEN
                  DO 120, I = 1, M
                     TEMP = TEMP + A( I, J )*X( IX )
                     IX   = IX   + INCX
  120             CONTINUE
               ELSE
                  DO 130, I = 1, M
                     TEMP = TEMP + CONJG( A( I, J ) )*X( IX )
                     IX   = IX   + INCX
  130             CONTINUE
               END IF
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  140       CONTINUE
         END IF
      END IF
!
      RETURN
!
!     End of CGEMV .
!
      END
      double precision function dcabs1(z)
      double complex z,zz
      double precision t(2)
      equivalence (zz,t(1))
      zz = z
      dcabs1 = dabs(t(1)) + dabs(t(2))
      return
      end
      SUBROUTINE DGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB,&
     &                   BETA, C, LDC )
!     .. Scalar Arguments ..
      CHARACTER*1        TRANSA, TRANSB
      INTEGER            M, N, K, LDA, LDB, LDC
      DOUBLE PRECISION   ALPHA, BETA
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * )
!     ..
!
!  Purpose
!  =======
!
!  DGEMM  performs one of the matrix-matrix operations
!
!     C := alpha*op( A )*op( B ) + beta*C,
!
!  where  op( X ) is one of
!
!     op( X ) = X   or   op( X ) = X',
!
!  alpha and beta are scalars, and A, B and C are matrices, with op( A )
!  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
!
!  Parameters
!  ==========
!
!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n',  op( A ) = A.
!
!              TRANSA = 'T' or 't',  op( A ) = A'.
!
!              TRANSA = 'C' or 'c',  op( A ) = A'.
!
!           Unchanged on exit.
!
!  TRANSB - CHARACTER*1.
!           On entry, TRANSB specifies the form of op( B ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSB = 'N' or 'n',  op( B ) = B.
!
!              TRANSB = 'T' or 't',  op( B ) = B'.
!
!              TRANSB = 'C' or 'c',  op( B ) = B'.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry,  M  specifies  the number  of rows  of the  matrix
!           op( A )  and of the  matrix  C.  M  must  be at least  zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry,  N  specifies the number  of columns of the matrix
!           op( B ) and the number of columns of the matrix C. N must be
!           at least zero.
!           Unchanged on exit.
!
!  K      - INTEGER.
!           On entry,  K  specifies  the number of columns of the matrix
!           op( A ) and the number of rows of the matrix op( B ). K must
!           be at least  zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
!           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
!           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
!           part of the array  A  must contain the matrix  A,  otherwise
!           the leading  k by m  part of the array  A  must contain  the
!           matrix A.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
!           LDA must be at least  max( 1, m ), otherwise  LDA must be at
!           least  max( 1, k ).
!           Unchanged on exit.
!
!  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
!           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
!           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
!           part of the array  B  must contain the matrix  B,  otherwise
!           the leading  n by k  part of the array  B  must contain  the
!           matrix B.
!           Unchanged on exit.
!
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
!           LDB must be at least  max( 1, k ), otherwise  LDB must be at
!           least  max( 1, n ).
!           Unchanged on exit.
!
!  BETA   - DOUBLE PRECISION.
!           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
!           supplied as zero then C need not be set on input.
!           Unchanged on exit.
!
!  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
!           Before entry, the leading  m by n  part of the array  C must
!           contain the matrix  C,  except when  beta  is zero, in which
!           case C need not be set on entry.
!           On exit, the array  C  is overwritten by the  m by n  matrix
!           ( alpha*op( A )*op( B ) + beta*C ).
!
!  LDC    - INTEGER.
!           On entry, LDC specifies the first dimension of C as declared
!           in  the  calling  (sub)  program.   LDC  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 3 Blas routine.
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     .. Local Scalars ..
      LOGICAL            NOTA, NOTB
      INTEGER            I, INFO, J, L, NCOLA, NROWA, NROWB
      DOUBLE PRECISION   TEMP
!     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Executable Statements ..
!
!     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
!     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
!     and  columns of  A  and the  number of  rows  of  B  respectively.
!
      NOTA  = LSAME( TRANSA, 'N' )
      NOTB  = LSAME( TRANSB, 'N' )
      IF( NOTA )THEN
         NROWA = M
         NCOLA = K
      ELSE
         NROWA = K
         NCOLA = M
      END IF
      IF( NOTB )THEN
         NROWB = K
      ELSE
         NROWB = N
      END IF
!
!     Test the input parameters.
!
      INFO = 0
      IF(      ( .NOT.NOTA                 ).AND.                       &
     &         ( .NOT.LSAME( TRANSA, 'C' ) ).AND.                       &
     &         ( .NOT.LSAME( TRANSA, 'T' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.NOTB                 ).AND.                       &
     &         ( .NOT.LSAME( TRANSB, 'C' ) ).AND.                       &
     &         ( .NOT.LSAME( TRANSB, 'T' ) )      )THEN
         INFO = 2
      ELSE IF( M  .LT.0               )THEN
         INFO = 3
      ELSE IF( N  .LT.0               )THEN
         INFO = 4
      ELSE IF( K  .LT.0               )THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 8
      ELSE IF( LDB.LT.MAX( 1, NROWB ) )THEN
         INFO = 10
      ELSE IF( LDC.LT.MAX( 1, M     ) )THEN
         INFO = 13
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGEMM ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.                                  &
     &    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) ) &
     &   RETURN
!
!     And if  alpha.eq.zero.
!
      IF( ALPHA.EQ.ZERO )THEN
         IF( BETA.EQ.ZERO )THEN
            DO 20, J = 1, N
               DO 10, I = 1, M
                  C( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               DO 30, I = 1, M
                  C( I, J ) = BETA*C( I, J )
   30          CONTINUE
   40       CONTINUE
         END IF
         RETURN
      END IF
!
!     Start the operations.
!
      IF( NOTB )THEN
         IF( NOTA )THEN
!
!           Form  C := alpha*A*B + beta*C.
!
            DO 90, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 50, I = 1, M
                     C( I, J ) = ZERO
   50             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 60, I = 1, M
                     C( I, J ) = BETA*C( I, J )
   60             CONTINUE
               END IF
               DO 80, L = 1, K
                  IF( B( L, J ).NE.ZERO )THEN
                     TEMP = ALPHA*B( L, J )
                     DO 70, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
   70                CONTINUE
                  END IF
   80          CONTINUE
   90       CONTINUE
         ELSE
!
!           Form  C := alpha*A'*B + beta*C
!
            DO 120, J = 1, N
               DO 110, I = 1, M
                  TEMP = ZERO
                  DO 100, L = 1, K
                     TEMP = TEMP + A( L, I )*B( L, J )
  100             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  110          CONTINUE
  120       CONTINUE
         END IF
      ELSE
         IF( NOTA )THEN
!
!           Form  C := alpha*A*B' + beta*C
!
            DO 170, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 130, I = 1, M
                     C( I, J ) = ZERO
  130             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 140, I = 1, M
                     C( I, J ) = BETA*C( I, J )
  140             CONTINUE
               END IF
               DO 160, L = 1, K
                  IF( B( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*B( J, L )
                     DO 150, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  150                CONTINUE
                  END IF
  160          CONTINUE
  170       CONTINUE
         ELSE
!
!           Form  C := alpha*A'*B' + beta*C
!
            DO 200, J = 1, N
               DO 190, I = 1, M
                  TEMP = ZERO
                  DO 180, L = 1, K
                     TEMP = TEMP + A( L, I )*B( J, L )
  180             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  190          CONTINUE
  200       CONTINUE
         END IF
      END IF
!
      RETURN
!
!     End of DGEMM .
!
      END
      SUBROUTINE DGER  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
!     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA
      INTEGER            INCX, INCY, LDA, M, N
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
!     ..
!
!  Purpose
!  =======
!
!  DGER   performs the rank 1 operation
!
!     A := alpha*x*y' + A,
!
!  where alpha is a scalar, x is an m element vector, y is an n element
!  vector and A is an m by n matrix.
!
!  Parameters
!  ==========
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  X      - DOUBLE PRECISION array of dimension at least
!           ( 1 + ( m - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the m
!           element vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  Y      - DOUBLE PRECISION array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCY ) ).
!           Before entry, the incremented array Y must contain the n
!           element vector y.
!           Unchanged on exit.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients. On exit, A is
!           overwritten by the updated matrix.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
!     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, J, JY, KX
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF     ( M.LT.0 )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 7
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGER  ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )               &
     &   RETURN
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
      IF( INCY.GT.0 )THEN
         JY = 1
      ELSE
         JY = 1 - ( N - 1 )*INCY
      END IF
      IF( INCX.EQ.1 )THEN
         DO 20, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               DO 10, I = 1, M
                  A( I, J ) = A( I, J ) + X( I )*TEMP
   10          CONTINUE
            END IF
            JY = JY + INCY
   20    CONTINUE
      ELSE
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( M - 1 )*INCX
         END IF
         DO 40, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               IX   = KX
               DO 30, I = 1, M
                  A( I, J ) = A( I, J ) + X( IX )*TEMP
                  IX        = IX        + INCX
   30          CONTINUE
            END IF
            JY = JY + INCY
   40    CONTINUE
      END IF
!
      RETURN
!
!     End of DGER  .
!
      END
      subroutine  dscal(n,da,dx,incx)
!
!     scales a vector by a constant.
!     uses unrolled loops for increment equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      double precision da,dx(*)
      integer i,incx,m,mp1,n,nincx
!
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
!
!        code for increment not equal to 1
!
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dx(i) = da*dx(i)
   10 continue
      return
!
!        code for increment equal to 1
!
!
!        clean-up loop
!
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dx(i) = da*dx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
   50 continue
      return
      end
      subroutine  dswap (n,dx,incx,dy,incy)
!
!     interchanges two vectors.
!     uses unrolled loops for increments equal one.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
!
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
!
!       code for unequal increments or equal increments not equal
!         to 1
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dx(ix)
        dx(ix) = dy(iy)
        dy(iy) = dtemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!
!       code for both increments equal to 1
!
!
!       clean-up loop
!
   20 m = mod(n,3)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
   30 continue
      if( n .lt. 3 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,3
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
        dtemp = dx(i + 1)
        dx(i + 1) = dy(i + 1)
        dy(i + 1) = dtemp
        dtemp = dx(i + 2)
        dx(i + 2) = dy(i + 2)
        dy(i + 2) = dtemp
   50 continue
      return
      end
      SUBROUTINE DTRSM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, &
     &                   B, LDB )
!     .. Scalar Arguments ..
      CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, LDA, LDB
      DOUBLE PRECISION   ALPHA
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  DTRSM  solves one of the matrix equations
!
!     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
!
!  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
!  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
!
!     op( A ) = A   or   op( A ) = A'.
!
!  The matrix X is overwritten on B.
!
!  Parameters
!  ==========
!
!  SIDE   - CHARACTER*1.
!           On entry, SIDE specifies whether op( A ) appears on the left
!           or right of X as follows:
!
!              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
!
!              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
!
!           Unchanged on exit.
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the matrix A is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n'   op( A ) = A.
!
!              TRANSA = 'T' or 't'   op( A ) = A'.
!
!              TRANSA = 'C' or 'c'   op( A ) = A'.
!
!           Unchanged on exit.
!
!  DIAG   - CHARACTER*1.
!           On entry, DIAG specifies whether or not A is unit triangular
!           as follows:
!
!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of B. M must be at
!           least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of B.  N must be
!           at least zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
!           zero then  A is not referenced and  B need not be set before
!           entry.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m
!           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
!           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
!           upper triangular part of the array  A must contain the upper
!           triangular matrix  and the strictly lower triangular part of
!           A is not referenced.
!           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
!           lower triangular part of the array  A must contain the lower
!           triangular matrix  and the strictly upper triangular part of
!           A is not referenced.
!           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
!           A  are not referenced either,  but are assumed to be  unity.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
!           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
!           then LDA must be at least max( 1, n ).
!           Unchanged on exit.
!
!  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).
!           Before entry,  the leading  m by n part of the array  B must
!           contain  the  right-hand  side  matrix  B,  and  on exit  is
!           overwritten by the solution matrix  X.
!
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in  the  calling  (sub)  program.   LDB  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 3 Blas routine.
!
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     .. Local Scalars ..
      LOGICAL            LSIDE, NOUNIT, UPPER
      INTEGER            I, INFO, J, K, NROWA
      DOUBLE PRECISION   TEMP
!     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      LSIDE  = LSAME( SIDE  , 'L' )
      IF( LSIDE )THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      NOUNIT = LSAME( DIAG  , 'N' )
      UPPER  = LSAME( UPLO  , 'U' )
!
      INFO   = 0
      IF(      ( .NOT.LSIDE                ).AND.                       &
     &         ( .NOT.LSAME( SIDE  , 'R' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.UPPER                ).AND.                       &
     &         ( .NOT.LSAME( UPLO  , 'L' ) )      )THEN
         INFO = 2
      ELSE IF( ( .NOT.LSAME( TRANSA, 'N' ) ).AND.                       &
     &         ( .NOT.LSAME( TRANSA, 'T' ) ).AND.                       &
     &         ( .NOT.LSAME( TRANSA, 'C' ) )      )THEN
         INFO = 3
      ELSE IF( ( .NOT.LSAME( DIAG  , 'U' ) ).AND.                       &
     &         ( .NOT.LSAME( DIAG  , 'N' ) )      )THEN
         INFO = 4
      ELSE IF( M  .LT.0               )THEN
         INFO = 5
      ELSE IF( N  .LT.0               )THEN
         INFO = 6
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 9
      ELSE IF( LDB.LT.MAX( 1, M     ) )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DTRSM ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( N.EQ.0 )                                                      &
     &   RETURN
!
!     And when  alpha.eq.zero.
!
      IF( ALPHA.EQ.ZERO )THEN
         DO 20, J = 1, N
            DO 10, I = 1, M
               B( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
         RETURN
      END IF
!
!     Start the operations.
!
      IF( LSIDE )THEN
         IF( LSAME( TRANSA, 'N' ) )THEN
!
!           Form  B := alpha*inv( A )*B.
!
            IF( UPPER )THEN
               DO 60, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 30, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   30                CONTINUE
                  END IF
                  DO 50, K = M, 1, -1
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT )                                    &
     &                     B( K, J ) = B( K, J )/A( K, K )
                        DO 40, I = 1, K - 1
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   40                   CONTINUE
                     END IF
   50             CONTINUE
   60          CONTINUE
            ELSE
               DO 100, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 70, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   70                CONTINUE
                  END IF
                  DO 90 K = 1, M
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT )                                    &
     &                     B( K, J ) = B( K, J )/A( K, K )
                        DO 80, I = K + 1, M
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   80                   CONTINUE
                     END IF
   90             CONTINUE
  100          CONTINUE
            END IF
         ELSE
!
!           Form  B := alpha*inv( A' )*B.
!
            IF( UPPER )THEN
               DO 130, J = 1, N
                  DO 120, I = 1, M
                     TEMP = ALPHA*B( I, J )
                     DO 110, K = 1, I - 1
                        TEMP = TEMP - A( K, I )*B( K, J )
  110                CONTINUE
                     IF( NOUNIT )                                       &
     &                  TEMP = TEMP/A( I, I )
                     B( I, J ) = TEMP
  120             CONTINUE
  130          CONTINUE
            ELSE
               DO 160, J = 1, N
                  DO 150, I = M, 1, -1
                     TEMP = ALPHA*B( I, J )
                     DO 140, K = I + 1, M
                        TEMP = TEMP - A( K, I )*B( K, J )
  140                CONTINUE
                     IF( NOUNIT )                                       &
     &                  TEMP = TEMP/A( I, I )
                     B( I, J ) = TEMP
  150             CONTINUE
  160          CONTINUE
            END IF
         END IF
      ELSE
         IF( LSAME( TRANSA, 'N' ) )THEN
!
!           Form  B := alpha*B*inv( A ).
!
            IF( UPPER )THEN
               DO 210, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 170, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  170                CONTINUE
                  END IF
                  DO 190, K = 1, J - 1
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 180, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  180                   CONTINUE
                     END IF
  190             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 200, I = 1, M
                        B( I, J ) = TEMP*B( I, J )
  200                CONTINUE
                  END IF
  210          CONTINUE
            ELSE
               DO 260, J = N, 1, -1
                  IF( ALPHA.NE.ONE )THEN
                     DO 220, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  220                CONTINUE
                  END IF
                  DO 240, K = J + 1, N
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 230, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  230                   CONTINUE
                     END IF
  240             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 250, I = 1, M
                       B( I, J ) = TEMP*B( I, J )
  250                CONTINUE
                  END IF
  260          CONTINUE
            END IF
         ELSE
!
!           Form  B := alpha*B*inv( A' ).
!
            IF( UPPER )THEN
               DO 310, K = N, 1, -1
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( K, K )
                     DO 270, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  270                CONTINUE
                  END IF
                  DO 290, J = 1, K - 1
                     IF( A( J, K ).NE.ZERO )THEN
                        TEMP = A( J, K )
                        DO 280, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  280                   CONTINUE
                     END IF
  290             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 300, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  300                CONTINUE
                  END IF
  310          CONTINUE
            ELSE
               DO 360, K = 1, N
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( K, K )
                     DO 320, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  320                CONTINUE
                  END IF
                  DO 340, J = K + 1, N
                     IF( A( J, K ).NE.ZERO )THEN
                        TEMP = A( J, K )
                        DO 330, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  330                   CONTINUE
                     END IF
  340             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 350, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  350                CONTINUE
                  END IF
  360          CONTINUE
            END IF
         END IF
      END IF
!
      RETURN
!
!     End of DTRSM .
!
      END
      integer function idamax(n,dx,incx)
!
!     finds the index of element having max. absolute value.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      double precision dx(*),dmax
      integer i,incx,ix,n
!
      idamax = 0
      if( n.lt.1 .or. incx.le.0 ) return
      idamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
!
!        code for increment not equal to 1
!
      ix = 1
      dmax = dabs(dx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(dabs(dx(ix)).le.dmax) go to 5
         idamax = i
         dmax = dabs(dx(ix))
    5    ix = ix + incx
   10 continue
      return
!
!        code for increment equal to 1
!
   20 dmax = dabs(dx(1))
      do 30 i = 2,n
         if(dabs(dx(i)).le.dmax) go to 30
         idamax = i
         dmax = dabs(dx(i))
   30 continue
      return
      end
      integer function izamax(n,zx,incx)
!
!     finds the index of element having max. absolute value.
!     jack dongarra, 1/15/85.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      double complex zx(*)
      double precision smax
      integer i,incx,ix,n
      double precision dcabs1
!
      izamax = 0
      if( n.lt.1 .or. incx.le.0 )return
      izamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
!
!        code for increment not equal to 1
!
      ix = 1
      smax = dcabs1(zx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(dcabs1(zx(ix)).le.smax) go to 5
         izamax = i
         smax = dcabs1(zx(ix))
    5    ix = ix + incx
   10 continue
      return
!
!        code for increment equal to 1
!
   20 smax = dcabs1(zx(1))
      do 30 i = 2,n
         if(dcabs1(zx(i)).le.smax) go to 30
         izamax = i
         smax = dcabs1(zx(i))
   30 continue
      return
      end
      subroutine  zcopy(n,zx,incx,zy,incy)
!
!     copies a vector, x, to a vector, y.
!     jack dongarra, linpack, 4/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      double complex zx(*),zy(*)
      integer i,incx,incy,ix,iy,n
!
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        zy(iy) = zx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!
!        code for both increments equal to 1
!
   20 do 30 i = 1,n
        zy(i) = zx(i)
   30 continue
      return
      end
      SUBROUTINE ZGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB,&
     &                   BETA, C, LDC )
!     .. Scalar Arguments ..
      CHARACTER*1        TRANSA, TRANSB
      INTEGER            M, N, K, LDA, LDB, LDC
      COMPLEX*16         ALPHA, BETA
!     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), B( LDB, * ), C( LDC, * )
!     ..
!
!  Purpose
!  =======
!
!  ZGEMM  performs one of the matrix-matrix operations
!
!     C := alpha*op( A )*op( B ) + beta*C,
!
!  where  op( X ) is one of
!
!     op( X ) = X   or   op( X ) = X'   or   op( X ) = conjg( X' ),
!
!  alpha and beta are scalars, and A, B and C are matrices, with op( A )
!  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
!
!  Parameters
!  ==========
!
!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n',  op( A ) = A.
!
!              TRANSA = 'T' or 't',  op( A ) = A'.
!
!              TRANSA = 'C' or 'c',  op( A ) = conjg( A' ).
!
!           Unchanged on exit.
!
!  TRANSB - CHARACTER*1.
!           On entry, TRANSB specifies the form of op( B ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSB = 'N' or 'n',  op( B ) = B.
!
!              TRANSB = 'T' or 't',  op( B ) = B'.
!
!              TRANSB = 'C' or 'c',  op( B ) = conjg( B' ).
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry,  M  specifies  the number  of rows  of the  matrix
!           op( A )  and of the  matrix  C.  M  must  be at least  zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry,  N  specifies the number  of columns of the matrix
!           op( B ) and the number of columns of the matrix C. N must be
!           at least zero.
!           Unchanged on exit.
!
!  K      - INTEGER.
!           On entry,  K  specifies  the number of columns of the matrix
!           op( A ) and the number of rows of the matrix op( B ). K must
!           be at least  zero.
!           Unchanged on exit.
!
!  ALPHA  - COMPLEX*16      .
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - COMPLEX*16       array of DIMENSION ( LDA, ka ), where ka is
!           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
!           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
!           part of the array  A  must contain the matrix  A,  otherwise
!           the leading  k by m  part of the array  A  must contain  the
!           matrix A.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
!           LDA must be at least  max( 1, m ), otherwise  LDA must be at
!           least  max( 1, k ).
!           Unchanged on exit.
!
!  B      - COMPLEX*16       array of DIMENSION ( LDB, kb ), where kb is
!           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
!           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
!           part of the array  B  must contain the matrix  B,  otherwise
!           the leading  n by k  part of the array  B  must contain  the
!           matrix B.
!           Unchanged on exit.
!
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
!           LDB must be at least  max( 1, k ), otherwise  LDB must be at
!           least  max( 1, n ).
!           Unchanged on exit.
!
!  BETA   - COMPLEX*16      .
!           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
!           supplied as zero then C need not be set on input.
!           Unchanged on exit.
!
!  C      - COMPLEX*16       array of DIMENSION ( LDC, n ).
!           Before entry, the leading  m by n  part of the array  C must
!           contain the matrix  C,  except when  beta  is zero, in which
!           case C need not be set on entry.
!           On exit, the array  C  is overwritten by the  m by n  matrix
!           ( alpha*op( A )*op( B ) + beta*C ).
!
!  LDC    - INTEGER.
!           On entry, LDC specifies the first dimension of C as declared
!           in  the  calling  (sub)  program.   LDC  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 3 Blas routine.
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC          DCONJG, MAX
!     .. Local Scalars ..
      LOGICAL            CONJA, CONJB, NOTA, NOTB
      INTEGER            I, INFO, J, L, NCOLA, NROWA, NROWB
      COMPLEX*16         TEMP
!     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER        ( ONE  = ( 1.0D+0, 0.0D+0 ) )
      COMPLEX*16         ZERO
      PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) )
!     ..
!     .. Executable Statements ..
!
!     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
!     conjugated or transposed, set  CONJA and CONJB  as true if  A  and
!     B  respectively are to be  transposed but  not conjugated  and set
!     NROWA, NCOLA and  NROWB  as the number of rows and  columns  of  A
!     and the number of rows of  B  respectively.
!
      NOTA  = LSAME( TRANSA, 'N' )
      NOTB  = LSAME( TRANSB, 'N' )
      CONJA = LSAME( TRANSA, 'C' )
      CONJB = LSAME( TRANSB, 'C' )
      IF( NOTA )THEN
         NROWA = M
         NCOLA = K
      ELSE
         NROWA = K
         NCOLA = M
      END IF
      IF( NOTB )THEN
         NROWB = K
      ELSE
         NROWB = N
      END IF
!
!     Test the input parameters.
!
      INFO = 0
      IF(      ( .NOT.NOTA                 ).AND.                       &
     &         ( .NOT.CONJA                ).AND.                       &
     &         ( .NOT.LSAME( TRANSA, 'T' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.NOTB                 ).AND.                       &
     &         ( .NOT.CONJB                ).AND.                       &
     &         ( .NOT.LSAME( TRANSB, 'T' ) )      )THEN
         INFO = 2
      ELSE IF( M  .LT.0               )THEN
         INFO = 3
      ELSE IF( N  .LT.0               )THEN
         INFO = 4
      ELSE IF( K  .LT.0               )THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 8
      ELSE IF( LDB.LT.MAX( 1, NROWB ) )THEN
         INFO = 10
      ELSE IF( LDC.LT.MAX( 1, M     ) )THEN
         INFO = 13
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'ZGEMM ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.                                  &
     &    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) ) &
     &   RETURN
!
!     And when  alpha.eq.zero.
!
      IF( ALPHA.EQ.ZERO )THEN
         IF( BETA.EQ.ZERO )THEN
            DO 20, J = 1, N
               DO 10, I = 1, M
                  C( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               DO 30, I = 1, M
                  C( I, J ) = BETA*C( I, J )
   30          CONTINUE
   40       CONTINUE
         END IF
         RETURN
      END IF
!
!     Start the operations.
!
      IF( NOTB )THEN
         IF( NOTA )THEN
!
!           Form  C := alpha*A*B + beta*C.
!
            DO 90, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 50, I = 1, M
                     C( I, J ) = ZERO
   50             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 60, I = 1, M
                     C( I, J ) = BETA*C( I, J )
   60             CONTINUE
               END IF
               DO 80, L = 1, K
                  IF( B( L, J ).NE.ZERO )THEN
                     TEMP = ALPHA*B( L, J )
                     DO 70, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
   70                CONTINUE
                  END IF
   80          CONTINUE
   90       CONTINUE
         ELSE IF( CONJA )THEN
!
!           Form  C := alpha*conjg( A' )*B + beta*C.
!
            DO 120, J = 1, N
               DO 110, I = 1, M
                  TEMP = ZERO
                  DO 100, L = 1, K
                     TEMP = TEMP + DCONJG( A( L, I ) )*B( L, J )
  100             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  110          CONTINUE
  120       CONTINUE
         ELSE
!
!           Form  C := alpha*A'*B + beta*C
!
            DO 150, J = 1, N
               DO 140, I = 1, M
                  TEMP = ZERO
                  DO 130, L = 1, K
                     TEMP = TEMP + A( L, I )*B( L, J )
  130             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  140          CONTINUE
  150       CONTINUE
         END IF
      ELSE IF( NOTA )THEN
         IF( CONJB )THEN
!
!           Form  C := alpha*A*conjg( B' ) + beta*C.
!
            DO 200, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 160, I = 1, M
                     C( I, J ) = ZERO
  160             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 170, I = 1, M
                     C( I, J ) = BETA*C( I, J )
  170             CONTINUE
               END IF
               DO 190, L = 1, K
                  IF( B( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*DCONJG( B( J, L ) )
                     DO 180, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  180                CONTINUE
                  END IF
  190          CONTINUE
  200       CONTINUE
         ELSE
!
!           Form  C := alpha*A*B'          + beta*C
!
            DO 250, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 210, I = 1, M
                     C( I, J ) = ZERO
  210             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 220, I = 1, M
                     C( I, J ) = BETA*C( I, J )
  220             CONTINUE
               END IF
               DO 240, L = 1, K
                  IF( B( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*B( J, L )
                     DO 230, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  230                CONTINUE
                  END IF
  240          CONTINUE
  250       CONTINUE
         END IF
      ELSE IF( CONJA )THEN
         IF( CONJB )THEN
!
!           Form  C := alpha*conjg( A' )*conjg( B' ) + beta*C.
!
            DO 280, J = 1, N
               DO 270, I = 1, M
                  TEMP = ZERO
                  DO 260, L = 1, K
                     TEMP = TEMP +                                      &
     &                      DCONJG( A( L, I ) )*DCONJG( B( J, L ) )
  260             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  270          CONTINUE
  280       CONTINUE
         ELSE
!
!           Form  C := alpha*conjg( A' )*B' + beta*C
!
            DO 310, J = 1, N
               DO 300, I = 1, M
                  TEMP = ZERO
                  DO 290, L = 1, K
                     TEMP = TEMP + DCONJG( A( L, I ) )*B( J, L )
  290             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  300          CONTINUE
  310       CONTINUE
         END IF
      ELSE
         IF( CONJB )THEN
!
!           Form  C := alpha*A'*conjg( B' ) + beta*C
!
            DO 340, J = 1, N
               DO 330, I = 1, M
                  TEMP = ZERO
                  DO 320, L = 1, K
                     TEMP = TEMP + A( L, I )*DCONJG( B( J, L ) )
  320             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  330          CONTINUE
  340       CONTINUE
         ELSE
!
!           Form  C := alpha*A'*B' + beta*C
!
            DO 370, J = 1, N
               DO 360, I = 1, M
                  TEMP = ZERO
                  DO 350, L = 1, K
                     TEMP = TEMP + A( L, I )*B( J, L )
  350             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  360          CONTINUE
  370       CONTINUE
         END IF
      END IF
!
      RETURN
!
!     End of ZGEMM .
!
      END
      SUBROUTINE ZGERU ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
!     .. Scalar Arguments ..
      COMPLEX*16         ALPHA
      INTEGER            INCX, INCY, LDA, M, N
!     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), X( * ), Y( * )
!     ..
!
!  Purpose
!  =======
!
!  ZGERU  performs the rank 1 operation
!
!     A := alpha*x*y' + A,
!
!  where alpha is a scalar, x is an m element vector, y is an n element
!  vector and A is an m by n matrix.
!
!  Parameters
!  ==========
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - COMPLEX*16      .
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  X      - COMPLEX*16       array of dimension at least
!           ( 1 + ( m - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the m
!           element vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  Y      - COMPLEX*16       array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCY ) ).
!           Before entry, the incremented array Y must contain the n
!           element vector y.
!           Unchanged on exit.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!  A      - COMPLEX*16       array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients. On exit, A is
!           overwritten by the updated matrix.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
!     .. Parameters ..
      COMPLEX*16         ZERO
      PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) )
!     .. Local Scalars ..
      COMPLEX*16         TEMP
      INTEGER            I, INFO, IX, J, JY, KX
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF     ( M.LT.0 )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 7
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'ZGERU ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )               &
     &   RETURN
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
      IF( INCY.GT.0 )THEN
         JY = 1
      ELSE
         JY = 1 - ( N - 1 )*INCY
      END IF
      IF( INCX.EQ.1 )THEN
         DO 20, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               DO 10, I = 1, M
                  A( I, J ) = A( I, J ) + X( I )*TEMP
   10          CONTINUE
            END IF
            JY = JY + INCY
   20    CONTINUE
      ELSE
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( M - 1 )*INCX
         END IF
         DO 40, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               IX   = KX
               DO 30, I = 1, M
                  A( I, J ) = A( I, J ) + X( IX )*TEMP
                  IX        = IX        + INCX
   30          CONTINUE
            END IF
            JY = JY + INCY
   40    CONTINUE
      END IF
!
      RETURN
!
!     End of ZGERU .
!
      END
      subroutine  zscal(n,za,zx,incx)
!
!     scales a vector by a constant.
!     jack dongarra, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      double complex za,zx(*)
      integer i,incx,ix,n
!
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
!
!        code for increment not equal to 1
!
      ix = 1
      do 10 i = 1,n
        zx(ix) = za*zx(ix)
        ix = ix + incx
   10 continue
      return
!
!        code for increment equal to 1
!
   20 do 30 i = 1,n
        zx(i) = za*zx(i)
   30 continue
      return
      end
      subroutine  zswap (n,zx,incx,zy,incy)
!
!     interchanges two vectors.
!     jack dongarra, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      double complex zx(*),zy(*),ztemp
      integer i,incx,incy,ix,iy,n
!
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
!
!       code for unequal increments or equal increments not equal
!         to 1
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        ztemp = zx(ix)
        zx(ix) = zy(iy)
        zy(iy) = ztemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!
!       code for both increments equal to 1
   20 do 30 i = 1,n
        ztemp = zx(i)
        zx(i) = zy(i)
        zy(i) = ztemp
   30 continue
      return
      end
      SUBROUTINE ZTRSM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, &
     &                   B, LDB )
!     .. Scalar Arguments ..
      CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, LDA, LDB
      COMPLEX*16         ALPHA
!     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  ZTRSM  solves one of the matrix equations
!
!     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
!
!  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
!  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
!
!     op( A ) = A   or   op( A ) = A'   or   op( A ) = conjg( A' ).
!
!  The matrix X is overwritten on B.
!
!  Parameters
!  ==========
!
!  SIDE   - CHARACTER*1.
!           On entry, SIDE specifies whether op( A ) appears on the left
!           or right of X as follows:
!
!              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
!
!              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
!
!           Unchanged on exit.
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the matrix A is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n'   op( A ) = A.
!
!              TRANSA = 'T' or 't'   op( A ) = A'.
!
!              TRANSA = 'C' or 'c'   op( A ) = conjg( A' ).
!
!           Unchanged on exit.
!
!  DIAG   - CHARACTER*1.
!           On entry, DIAG specifies whether or not A is unit triangular
!           as follows:
!
!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of B. M must be at
!           least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of B.  N must be
!           at least zero.
!           Unchanged on exit.
!
!  ALPHA  - COMPLEX*16      .
!           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
!           zero then  A is not referenced and  B need not be set before
!           entry.
!           Unchanged on exit.
!
!  A      - COMPLEX*16       array of DIMENSION ( LDA, k ), where k is m
!           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
!           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
!           upper triangular part of the array  A must contain the upper
!           triangular matrix  and the strictly lower triangular part of
!           A is not referenced.
!           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
!           lower triangular part of the array  A must contain the lower
!           triangular matrix  and the strictly upper triangular part of
!           A is not referenced.
!           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
!           A  are not referenced either,  but are assumed to be  unity.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
!           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
!           then LDA must be at least max( 1, n ).
!           Unchanged on exit.
!
!  B      - COMPLEX*16       array of DIMENSION ( LDB, n ).
!           Before entry,  the leading  m by n part of the array  B must
!           contain  the  right-hand  side  matrix  B,  and  on exit  is
!           overwritten by the solution matrix  X.
!
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in  the  calling  (sub)  program.   LDB  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 3 Blas routine.
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC          DCONJG, MAX
!     .. Local Scalars ..
      LOGICAL            LSIDE, NOCONJ, NOUNIT, UPPER
      INTEGER            I, INFO, J, K, NROWA
      COMPLEX*16         TEMP
!     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER        ( ONE  = ( 1.0D+0, 0.0D+0 ) )
      COMPLEX*16         ZERO
      PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) )
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      LSIDE  = LSAME( SIDE  , 'L' )
      IF( LSIDE )THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      NOCONJ = LSAME( TRANSA, 'T' )
      NOUNIT = LSAME( DIAG  , 'N' )
      UPPER  = LSAME( UPLO  , 'U' )
!
      INFO   = 0
      IF(      ( .NOT.LSIDE                ).AND.                       &
     &         ( .NOT.LSAME( SIDE  , 'R' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.UPPER                ).AND.                       &
     &         ( .NOT.LSAME( UPLO  , 'L' ) )      )THEN
         INFO = 2
      ELSE IF( ( .NOT.LSAME( TRANSA, 'N' ) ).AND.                       &
     &         ( .NOT.LSAME( TRANSA, 'T' ) ).AND.                       &
     &         ( .NOT.LSAME( TRANSA, 'C' ) )      )THEN
         INFO = 3
      ELSE IF( ( .NOT.LSAME( DIAG  , 'U' ) ).AND.                       &
     &         ( .NOT.LSAME( DIAG  , 'N' ) )      )THEN
         INFO = 4
      ELSE IF( M  .LT.0               )THEN
         INFO = 5
      ELSE IF( N  .LT.0               )THEN
         INFO = 6
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 9
      ELSE IF( LDB.LT.MAX( 1, M     ) )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'ZTRSM ', INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( N.EQ.0 )                                                      &
     &   RETURN
!
!     And when  alpha.eq.zero.
!
      IF( ALPHA.EQ.ZERO )THEN
         DO 20, J = 1, N
            DO 10, I = 1, M
               B( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
         RETURN
      END IF
!
!     Start the operations.
!
      IF( LSIDE )THEN
         IF( LSAME( TRANSA, 'N' ) )THEN
!
!           Form  B := alpha*inv( A )*B.
!
            IF( UPPER )THEN
               DO 60, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 30, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   30                CONTINUE
                  END IF
                  DO 50, K = M, 1, -1
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT )                                    &
     &                     B( K, J ) = B( K, J )/A( K, K )
                        DO 40, I = 1, K - 1
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   40                   CONTINUE
                     END IF
   50             CONTINUE
   60          CONTINUE
            ELSE
               DO 100, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 70, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   70                CONTINUE
                  END IF
                  DO 90 K = 1, M
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT )                                    &
     &                     B( K, J ) = B( K, J )/A( K, K )
                        DO 80, I = K + 1, M
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   80                   CONTINUE
                     END IF
   90             CONTINUE
  100          CONTINUE
            END IF
         ELSE
!
!           Form  B := alpha*inv( A' )*B
!           or    B := alpha*inv( conjg( A' ) )*B.
!
            IF( UPPER )THEN
               DO 140, J = 1, N
                  DO 130, I = 1, M
                     TEMP = ALPHA*B( I, J )
                     IF( NOCONJ )THEN
                        DO 110, K = 1, I - 1
                           TEMP = TEMP - A( K, I )*B( K, J )
  110                   CONTINUE
                        IF( NOUNIT )                                    &
     &                     TEMP = TEMP/A( I, I )
                     ELSE
                        DO 120, K = 1, I - 1
                           TEMP = TEMP - DCONJG( A( K, I ) )*B( K, J )
  120                   CONTINUE
                        IF( NOUNIT )                                    &
     &                     TEMP = TEMP/DCONJG( A( I, I ) )
                     END IF
                     B( I, J ) = TEMP
  130             CONTINUE
  140          CONTINUE
            ELSE
               DO 180, J = 1, N
                  DO 170, I = M, 1, -1
                     TEMP = ALPHA*B( I, J )
                     IF( NOCONJ )THEN
                        DO 150, K = I + 1, M
                           TEMP = TEMP - A( K, I )*B( K, J )
  150                   CONTINUE
                        IF( NOUNIT )                                    &
     &                     TEMP = TEMP/A( I, I )
                     ELSE
                        DO 160, K = I + 1, M
                           TEMP = TEMP - DCONJG( A( K, I ) )*B( K, J )
  160                   CONTINUE
                        IF( NOUNIT )                                    &
     &                     TEMP = TEMP/DCONJG( A( I, I ) )
                     END IF
                     B( I, J ) = TEMP
  170             CONTINUE
  180          CONTINUE
            END IF
         END IF
      ELSE
         IF( LSAME( TRANSA, 'N' ) )THEN
!
!           Form  B := alpha*B*inv( A ).
!
            IF( UPPER )THEN
               DO 230, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 190, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  190                CONTINUE
                  END IF
                  DO 210, K = 1, J - 1
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 200, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  200                   CONTINUE
                     END IF
  210             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 220, I = 1, M
                        B( I, J ) = TEMP*B( I, J )
  220                CONTINUE
                  END IF
  230          CONTINUE
            ELSE
               DO 280, J = N, 1, -1
                  IF( ALPHA.NE.ONE )THEN
                     DO 240, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  240                CONTINUE
                  END IF
                  DO 260, K = J + 1, N
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 250, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  250                   CONTINUE
                     END IF
  260             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 270, I = 1, M
                       B( I, J ) = TEMP*B( I, J )
  270                CONTINUE
                  END IF
  280          CONTINUE
            END IF
         ELSE
!
!           Form  B := alpha*B*inv( A' )
!           or    B := alpha*B*inv( conjg( A' ) ).
!
            IF( UPPER )THEN
               DO 330, K = N, 1, -1
                  IF( NOUNIT )THEN
                     IF( NOCONJ )THEN
                        TEMP = ONE/A( K, K )
                     ELSE
                        TEMP = ONE/DCONJG( A( K, K ) )
                     END IF
                     DO 290, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  290                CONTINUE
                  END IF
                  DO 310, J = 1, K - 1
                     IF( A( J, K ).NE.ZERO )THEN
                        IF( NOCONJ )THEN
                           TEMP = A( J, K )
                        ELSE
                           TEMP = DCONJG( A( J, K ) )
                        END IF
                        DO 300, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  300                   CONTINUE
                     END IF
  310             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 320, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  320                CONTINUE
                  END IF
  330          CONTINUE
            ELSE
               DO 380, K = 1, N
                  IF( NOUNIT )THEN
                     IF( NOCONJ )THEN
                        TEMP = ONE/A( K, K )
                     ELSE
                        TEMP = ONE/DCONJG( A( K, K ) )
                     END IF
                     DO 340, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  340                CONTINUE
                  END IF
                  DO 360, J = K + 1, N
                     IF( A( J, K ).NE.ZERO )THEN
                        IF( NOCONJ )THEN
                           TEMP = A( J, K )
                        ELSE
                           TEMP = DCONJG( A( J, K ) )
                        END IF
                        DO 350, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  350                   CONTINUE
                     END IF
  360             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 370, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  370                CONTINUE
                  END IF
  380          CONTINUE
            END IF
         END IF
      END IF
!
      RETURN
!
!     End of ZTRSM .
!
      END
      SUBROUTINE CLACON( N, V, X, EST, KASE )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      INTEGER            KASE, N
      REAL               EST
!     ..
!     .. Array Arguments ..
      COMPLEX            V( N ), X( N )
!     ..
!
!  Purpose
!  =======
!
!  CLACON estimates the 1-norm of a square, complex matrix A.
!  Reverse communication is used for evaluating matrix-vector products.
!
!  Arguments
!  =========
!
!  N      (input) INTEGER
!         The order of the matrix.  N >= 1.
!
!  V      (workspace) COMPLEX array, dimension (N)
!         On the final return, V = A*W,  where  EST = norm(V)/norm(W)
!         (W is not returned).
!
!  X      (input/output) COMPLEX array, dimension (N)
!         On an intermediate return, X should be overwritten by
!               A * X,   if KASE=1,
!               A' * X,  if KASE=2,
!         where A' is the conjugate transpose of A, and CLACON must be
!         re-called with all the other parameters unchanged.
!
!  EST    (output) REAL
!         An estimate (a lower bound) for norm(A).
!
!  KASE   (input/output) INTEGER
!         On the initial call to CLACON, KASE should be 0.
!         On an intermediate return, KASE will be 1 or 2, indicating
!         whether X should be overwritten by A * X  or A' * X.
!         On the final return from CLACON, KASE will again be 0.
!
!  Further Details
!  ======= =======
!
!  Contributed by Nick Higham, University of Manchester.
!  Originally named CONEST, dated March 16, 1988.
!
!  Reference: N.J. Higham, "FORTRAN codes for estimating the one-norm of
!  a real or complex matrix, with applications to condition estimation",
!  ACM Trans. Math. Soft., vol. 14, no. 4, pp. 381-396, December 1988.
!
!  Last modified:  April, 1999
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER            ITMAX
      PARAMETER          ( ITMAX = 5 )
      REAL               ONE, TWO
      PARAMETER          ( ONE = 1.0E0, TWO = 2.0E0 )
      COMPLEX            CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0E0, 0.0E0 ),                    &
     &                   CONE = ( 1.0E0, 0.0E0 ) )
!     ..
!     .. Local Scalars ..
      INTEGER            I, ITER, J, JLAST, JUMP
      REAL               ABSXI, ALTSGN, ESTOLD, SAFMIN, TEMP
!     ..
!     .. External Functions ..
      INTEGER            ICMAX1
      REAL               SCSUM1, SLAMCH
      EXTERNAL           ICMAX1, SCSUM1, SLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL           CCOPY
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, AIMAG, CMPLX, REAL
!     ..
!     .. Save statement ..
      SAVE
!     ..
!     .. Executable Statements ..
!
      SAFMIN = SLAMCH( 'Safe minimum' )
      IF( KASE.EQ.0 ) THEN
         DO 10 I = 1, N
            X( I ) = CMPLX( ONE / REAL( N ) )
   10    CONTINUE
         KASE = 1
         JUMP = 1
         RETURN
      END IF
!
      GO TO ( 20, 40, 70, 90, 120 )JUMP
!
!     ................ ENTRY   (JUMP = 1)
!     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X.
!
   20 CONTINUE
      IF( N.EQ.1 ) THEN
         V( 1 ) = X( 1 )
         EST = ABS( V( 1 ) )
!        ... QUIT
         GO TO 130
      END IF
      EST = SCSUM1( N, X, 1 )
!
      DO 30 I = 1, N
         ABSXI = ABS( X( I ) )
         IF( ABSXI.GT.SAFMIN ) THEN
            X( I ) = CMPLX( REAL( X( I ) ) / ABSXI,                     &
     &               AIMAG( X( I ) ) / ABSXI )
         ELSE
            X( I ) = CONE
         END IF
   30 CONTINUE
      KASE = 2
      JUMP = 2
      RETURN
!
!     ................ ENTRY   (JUMP = 2)
!     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY CTRANS(A)*X.
!
   40 CONTINUE
      J = ICMAX1( N, X, 1 )
      ITER = 2
!
!     MAIN LOOP - ITERATIONS 2,3,...,ITMAX.
!
   50 CONTINUE
      DO 60 I = 1, N
         X( I ) = CZERO
   60 CONTINUE
      X( J ) = CONE
      KASE = 1
      JUMP = 3
      RETURN
!
!     ................ ENTRY   (JUMP = 3)
!     X HAS BEEN OVERWRITTEN BY A*X.
!
   70 CONTINUE
      CALL CCOPY( N, X, 1, V, 1 )
      ESTOLD = EST
      EST = SCSUM1( N, V, 1 )
!
!     TEST FOR CYCLING.
      IF( EST.LE.ESTOLD )                                               &
     &   GO TO 100
!
      DO 80 I = 1, N
         ABSXI = ABS( X( I ) )
         IF( ABSXI.GT.SAFMIN ) THEN
            X( I ) = CMPLX( REAL( X( I ) ) / ABSXI,                     &
     &               AIMAG( X( I ) ) / ABSXI )
         ELSE
            X( I ) = CONE
         END IF
   80 CONTINUE
      KASE = 2
      JUMP = 4
      RETURN
!
!     ................ ENTRY   (JUMP = 4)
!     X HAS BEEN OVERWRITTEN BY CTRANS(A)*X.
!
   90 CONTINUE
      JLAST = J
      J = ICMAX1( N, X, 1 )
      IF( ( ABS( X( JLAST ) ).NE.ABS( X( J ) ) ) .AND.                  &
     &    ( ITER.LT.ITMAX ) ) THEN
         ITER = ITER + 1
         GO TO 50
      END IF
!
!     ITERATION COMPLETE.  FINAL STAGE.
!
  100 CONTINUE
      ALTSGN = ONE
      DO 110 I = 1, N
         X( I ) = CMPLX( ALTSGN*( ONE+REAL( I-1 ) / REAL( N-1 ) ) )
         ALTSGN = -ALTSGN
  110 CONTINUE
      KASE = 1
      JUMP = 5
      RETURN
!
!     ................ ENTRY   (JUMP = 5)
!     X HAS BEEN OVERWRITTEN BY A*X.
!
  120 CONTINUE
      TEMP = TWO*( SCSUM1( N, X, 1 ) / REAL( 3*N ) )
      IF( TEMP.GT.EST ) THEN
         CALL CCOPY( N, X, 1, V, 1 )
         EST = TEMP
      END IF
!
  130 CONTINUE
      KASE = 0
      RETURN
!
!     End of CLACON
!
      END
      SUBROUTINE ZGBTF2( M, N, KL, KU, AB, LDAB, IPIV, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      INTEGER            INFO, KL, KU, LDAB, M, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         AB( LDAB, * )
!     ..
!
!  Purpose
!  =======
!
!  ZGBTF2 computes an LU factorization of a complex m-by-n band matrix
!  A using partial pivoting with row interchanges.
!
!  This is the unblocked version of the algorithm, calling Level 2 BLAS.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  KL      (input) INTEGER
!          The number of subdiagonals within the band of A.  KL >= 0.
!
!  KU      (input) INTEGER
!          The number of superdiagonals within the band of A.  KU >= 0.
!
!  AB      (input/output) COMPLEX*16 array, dimension (LDAB,N)
!          On entry, the matrix A in band storage, in rows KL+1 to
!          2*KL+KU+1; rows 1 to KL of the array need not be set.
!          The j-th column of A is stored in the j-th column of the
!          array AB as follows:
!          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
!
!          On exit, details of the factorization: U is stored as an
!          upper triangular band matrix with KL+KU superdiagonals in
!          rows 1 to KL+KU+1, and the multipliers used during the
!          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
!          See below for further details.
!
!  LDAB    (input) INTEGER
!          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!          > 0: if INFO = +i, U(i,i) is exactly zero. The factorization
!               has been completed, but the factor U is exactly
!               singular, and division by zero will occur if it is used
!               to solve a system of equations.
!
!  Further Details
!  ===============
!
!  The band storage scheme is illustrated by the following example, when
!  M = N = 6, KL = 2, KU = 1:
!
!  On entry:                       On exit:
!
!      *    *    *    +    +    +       *    *    *   u14  u25  u36
!      *    *    +    +    +    +       *    *   u13  u24  u35  u46
!      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
!     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
!     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
!     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *
!
!  Array elements marked * are not used by the routine; elements marked
!  + need not be set on entry, but are required by the routine to store
!  elements of U, because of fill-in resulting from the row
!  interchanges.
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ),                    &
     &                   ZERO = ( 0.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      INTEGER            I, J, JP, JU, KM, KV
!     ..
!     .. External Functions ..
      INTEGER            IZAMAX
      EXTERNAL           IZAMAX
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZGERU, ZSCAL, ZSWAP
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     KV is the number of superdiagonals in the factor U, allowing for
!     fill-in.
!
      KV = KU + KL
!
!     Test the input parameters.
!
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KL.LT.0 ) THEN
         INFO = -3
      ELSE IF( KU.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDAB.LT.KL+KV+1 ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGBTF2', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 )                                          &
     &   RETURN
!
!     Gaussian elimination with partial pivoting
!
!     Set fill-in elements in columns KU+2 to KV to zero.
!
      DO 20 J = KU + 2, MIN( KV, N )
         DO 10 I = KV - J + 2, KL
            AB( I, J ) = ZERO
   10    CONTINUE
   20 CONTINUE
!
!     JU is the index of the last column affected by the current stage
!     of the factorization.
!
      JU = 1
!
      DO 40 J = 1, MIN( M, N )
!
!        Set fill-in elements in column J+KV to zero.
!
         IF( J+KV.LE.N ) THEN
            DO 30 I = 1, KL
               AB( I, J+KV ) = ZERO
   30       CONTINUE
         END IF
!
!        Find pivot and test for singularity. KM is the number of
!        subdiagonal elements in the current column.
!
         KM = MIN( KL, M-J )
         JP = IZAMAX( KM+1, AB( KV+1, J ), 1 )
         IPIV( J ) = JP + J - 1
         IF( AB( KV+JP, J ).NE.ZERO ) THEN
            JU = MAX( JU, MIN( J+KU+JP-1, N ) )
!
!           Apply interchange to columns J to JU.
!
            IF( JP.NE.1 )                                               &
     &         CALL ZSWAP( JU-J+1, AB( KV+JP, J ), LDAB-1,              &
     &                     AB( KV+1, J ), LDAB-1 )
            IF( KM.GT.0 ) THEN
!
!              Compute multipliers.
!
               CALL ZSCAL( KM, ONE / AB( KV+1, J ), AB( KV+2, J ), 1 )
!
!              Update trailing submatrix within the band.
!
               IF( JU.GT.J )                                            &
     &            CALL ZGERU( KM, JU-J, -ONE, AB( KV+2, J ), 1,         &
     &                        AB( KV, J+1 ), LDAB-1, AB( KV+1, J+1 ),   &
     &                        LDAB-1 )
            END IF
         ELSE
!
!           If pivot is zero, set INFO to the index of the pivot
!           unless a zero pivot has already been found.
!
            IF( INFO.EQ.0 )                                             &
     &         INFO = J
         END IF
   40 CONTINUE
      RETURN
!
!     End of ZGBTF2
!
      END
      SUBROUTINE ZGBTRF( M, N, KL, KU, AB, LDAB, IPIV, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      INTEGER            INFO, KL, KU, LDAB, M, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         AB( LDAB, * )
!     ..
!
!  Purpose
!  =======
!
!  ZGBTRF computes an LU factorization of a complex m-by-n band matrix A
!  using partial pivoting with row interchanges.
!
!  This is the blocked version of the algorithm, calling Level 3 BLAS.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  KL      (input) INTEGER
!          The number of subdiagonals within the band of A.  KL >= 0.
!
!  KU      (input) INTEGER
!          The number of superdiagonals within the band of A.  KU >= 0.
!
!  AB      (input/output) COMPLEX*16 array, dimension (LDAB,N)
!          On entry, the matrix A in band storage, in rows KL+1 to
!          2*KL+KU+1; rows 1 to KL of the array need not be set.
!          The j-th column of A is stored in the j-th column of the
!          array AB as follows:
!          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
!
!          On exit, details of the factorization: U is stored as an
!          upper triangular band matrix with KL+KU superdiagonals in
!          rows 1 to KL+KU+1, and the multipliers used during the
!          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
!          See below for further details.
!
!  LDAB    (input) INTEGER
!          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!          > 0: if INFO = +i, U(i,i) is exactly zero. The factorization
!               has been completed, but the factor U is exactly
!               singular, and division by zero will occur if it is used
!               to solve a system of equations.
!
!  Further Details
!  ===============
!
!  The band storage scheme is illustrated by the following example, when
!  M = N = 6, KL = 2, KU = 1:
!
!  On entry:                       On exit:
!
!      *    *    *    +    +    +       *    *    *   u14  u25  u36
!      *    *    +    +    +    +       *    *   u13  u24  u35  u46
!      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
!     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
!     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
!     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *
!
!  Array elements marked * are not used by the routine; elements marked
!  + need not be set on entry, but are required by the routine to store
!  elements of U because of fill-in resulting from the row interchanges.
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ),                    &
     &                   ZERO = ( 0.0D+0, 0.0D+0 ) )
      INTEGER            NBMAX, LDWORK
      PARAMETER          ( NBMAX = 64, LDWORK = NBMAX+1 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, I2, I3, II, IP, J, J2, J3, JB, JJ, JM, JP,  &
     &                   JU, K2, KM, KV, NB, NW
      COMPLEX*16         TEMP
!     ..
!     .. Local Arrays ..
      COMPLEX*16         WORK13( LDWORK, NBMAX ),                       &
     &                   WORK31( LDWORK, NBMAX )
!     ..
!     .. External Functions ..
      INTEGER            ILAENV, IZAMAX
      EXTERNAL           ILAENV, IZAMAX
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZCOPY, ZGBTF2, ZGEMM, ZGERU, ZLASWP,   &
     &                   ZSCAL, ZSWAP, ZTRSM
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     KV is the number of superdiagonals in the factor U, allowing for
!     fill-in
!
      KV = KU + KL
!
!     Test the input parameters.
!
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KL.LT.0 ) THEN
         INFO = -3
      ELSE IF( KU.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDAB.LT.KL+KV+1 ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGBTRF', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 )                                          &
     &   RETURN
!
!     Determine the block size for this environment
!
      NB = ILAENV( 1, 'ZGBTRF', ' ', M, N, KL, KU )
!
!     The block size must not exceed the limit set by the size of the
!     local arrays WORK13 and WORK31.
!
      NB = MIN( NB, NBMAX )
!
      IF( NB.LE.1 .OR. NB.GT.KL ) THEN
!
!        Use unblocked code
!
         CALL ZGBTF2( M, N, KL, KU, AB, LDAB, IPIV, INFO )
      ELSE
!
!        Use blocked code
!
!        Zero the superdiagonal elements of the work array WORK13
!
         DO 20 J = 1, NB
            DO 10 I = 1, J - 1
               WORK13( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
!
!        Zero the subdiagonal elements of the work array WORK31
!
         DO 40 J = 1, NB
            DO 30 I = J + 1, NB
               WORK31( I, J ) = ZERO
   30       CONTINUE
   40    CONTINUE
!
!        Gaussian elimination with partial pivoting
!
!        Set fill-in elements in columns KU+2 to KV to zero
!
         DO 60 J = KU + 2, MIN( KV, N )
            DO 50 I = KV - J + 2, KL
               AB( I, J ) = ZERO
   50       CONTINUE
   60    CONTINUE
!
!        JU is the index of the last column affected by the current
!        stage of the factorization
!
         JU = 1
!
         DO 180 J = 1, MIN( M, N ), NB
            JB = MIN( NB, MIN( M, N )-J+1 )
!
!           The active part of the matrix is partitioned
!
!              A11   A12   A13
!              A21   A22   A23
!              A31   A32   A33
!
!           Here A11, A21 and A31 denote the current block of JB columns
!           which is about to be factorized. The number of rows in the
!           partitioning are JB, I2, I3 respectively, and the numbers
!           of columns are JB, J2, J3. The superdiagonal elements of A13
!           and the subdiagonal elements of A31 lie outside the band.
!
            I2 = MIN( KL-JB, M-J-JB+1 )
            I3 = MIN( JB, M-J-KL+1 )
!
!           J2 and J3 are computed after JU has been updated.
!
!           Factorize the current block of JB columns
!
            DO 80 JJ = J, J + JB - 1
!
!              Set fill-in elements in column JJ+KV to zero
!
               IF( JJ+KV.LE.N ) THEN
                  DO 70 I = 1, KL
                     AB( I, JJ+KV ) = ZERO
   70             CONTINUE
               END IF
!
!              Find pivot and test for singularity. KM is the number of
!              subdiagonal elements in the current column.
!
               KM = MIN( KL, M-JJ )
               JP = IZAMAX( KM+1, AB( KV+1, JJ ), 1 )
               IPIV( JJ ) = JP + JJ - J
               IF( AB( KV+JP, JJ ).NE.ZERO ) THEN
                  JU = MAX( JU, MIN( JJ+KU+JP-1, N ) )
                  IF( JP.NE.1 ) THEN
!
!                    Apply interchange to columns J to J+JB-1
!
                     IF( JP+JJ-1.LT.J+KL ) THEN
!
                        CALL ZSWAP( JB, AB( KV+1+JJ-J, J ), LDAB-1,     &
     &                              AB( KV+JP+JJ-J, J ), LDAB-1 )
                     ELSE
!
!                       The interchange affects columns J to JJ-1 of A31
!                       which are stored in the work array WORK31
!
                        CALL ZSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1,   &
     &                              WORK31( JP+JJ-J-KL, 1 ), LDWORK )
                        CALL ZSWAP( J+JB-JJ, AB( KV+1, JJ ), LDAB-1,    &
     &                              AB( KV+JP, JJ ), LDAB-1 )
                     END IF
                  END IF
!
!                 Compute multipliers
!
                  CALL ZSCAL( KM, ONE / AB( KV+1, JJ ), AB( KV+2, JJ ), &
     &                        1 )
!
!                 Update trailing submatrix within the band and within
!                 the current block. JM is the index of the last column
!                 which needs to be updated.
!
                  JM = MIN( JU, J+JB-1 )
                  IF( JM.GT.JJ )                                        &
     &               CALL ZGERU( KM, JM-JJ, -ONE, AB( KV+2, JJ ), 1,    &
     &                           AB( KV, JJ+1 ), LDAB-1,                &
     &                           AB( KV+1, JJ+1 ), LDAB-1 )
               ELSE
!
!                 If pivot is zero, set INFO to the index of the pivot
!                 unless a zero pivot has already been found.
!
                  IF( INFO.EQ.0 )                                       &
     &               INFO = JJ
               END IF
!
!              Copy current column of A31 into the work array WORK31
!
               NW = MIN( JJ-J+1, I3 )
               IF( NW.GT.0 )                                            &
     &            CALL ZCOPY( NW, AB( KV+KL+1-JJ+J, JJ ), 1,            &
     &                        WORK31( 1, JJ-J+1 ), 1 )
   80       CONTINUE
            IF( J+JB.LE.N ) THEN
!
!              Apply the row interchanges to the other blocks.
!
               J2 = MIN( JU-J+1, KV ) - JB
               J3 = MAX( 0, JU-J-KV+1 )
!
!              Use ZLASWP to apply the row interchanges to A12, A22, and
!              A32.
!
               CALL ZLASWP( J2, AB( KV+1-JB, J+JB ), LDAB-1, 1, JB,     &
     &                      IPIV( J ), 1 )
!
!              Adjust the pivot indices.
!
               DO 90 I = J, J + JB - 1
                  IPIV( I ) = IPIV( I ) + J - 1
   90          CONTINUE
!
!              Apply the row interchanges to A13, A23, and A33
!              columnwise.
!
               K2 = J - 1 + JB + J2
               DO 110 I = 1, J3
                  JJ = K2 + I
                  DO 100 II = J + I - 1, J + JB - 1
                     IP = IPIV( II )
                     IF( IP.NE.II ) THEN
                        TEMP = AB( KV+1+II-JJ, JJ )
                        AB( KV+1+II-JJ, JJ ) = AB( KV+1+IP-JJ, JJ )
                        AB( KV+1+IP-JJ, JJ ) = TEMP
                     END IF
  100             CONTINUE
  110          CONTINUE
!
!              Update the relevant part of the trailing submatrix
!
               IF( J2.GT.0 ) THEN
!
!                 Update A12
!
                  CALL ZTRSM( 'Left', 'Lower', 'No transpose', 'Unit',  &
     &                        JB, J2, ONE, AB( KV+1, J ), LDAB-1,       &
     &                        AB( KV+1-JB, J+JB ), LDAB-1 )
!
                  IF( I2.GT.0 ) THEN
!
!                    Update A22
!
                     CALL ZGEMM( 'No transpose', 'No transpose', I2, J2,&
     &                           JB, -ONE, AB( KV+1+JB, J ), LDAB-1,    &
     &                           AB( KV+1-JB, J+JB ), LDAB-1, ONE,      &
     &                           AB( KV+1, J+JB ), LDAB-1 )
                  END IF
!
                  IF( I3.GT.0 ) THEN
!
!                    Update A32
!
                     CALL ZGEMM( 'No transpose', 'No transpose', I3, J2,&
     &                           JB, -ONE, WORK31, LDWORK,              &
     &                           AB( KV+1-JB, J+JB ), LDAB-1, ONE,      &
     &                           AB( KV+KL+1-JB, J+JB ), LDAB-1 )
                  END IF
               END IF
!
               IF( J3.GT.0 ) THEN
!
!                 Copy the lower triangle of A13 into the work array
!                 WORK13
!
                  DO 130 JJ = 1, J3
                     DO 120 II = JJ, JB
                        WORK13( II, JJ ) = AB( II-JJ+1, JJ+J+KV-1 )
  120                CONTINUE
  130             CONTINUE
!
!                 Update A13 in the work array
!
                  CALL ZTRSM( 'Left', 'Lower', 'No transpose', 'Unit',  &
     &                        JB, J3, ONE, AB( KV+1, J ), LDAB-1,       &
     &                        WORK13, LDWORK )
!
                  IF( I2.GT.0 ) THEN
!
!                    Update A23
!
                     CALL ZGEMM( 'No transpose', 'No transpose', I2, J3,&
     &                           JB, -ONE, AB( KV+1+JB, J ), LDAB-1,    &
     &                           WORK13, LDWORK, ONE, AB( 1+JB, J+KV ), &
     &                           LDAB-1 )
                  END IF
!
                  IF( I3.GT.0 ) THEN
!
!                    Update A33
!
                     CALL ZGEMM( 'No transpose', 'No transpose', I3, J3,&
     &                           JB, -ONE, WORK31, LDWORK, WORK13,      &
     &                           LDWORK, ONE, AB( 1+KL, J+KV ), LDAB-1 )
                  END IF
!
!                 Copy the lower triangle of A13 back into place
!
                  DO 150 JJ = 1, J3
                     DO 140 II = JJ, JB
                        AB( II-JJ+1, JJ+J+KV-1 ) = WORK13( II, JJ )
  140                CONTINUE
  150             CONTINUE
               END IF
            ELSE
!
!              Adjust the pivot indices.
!
               DO 160 I = J, J + JB - 1
                  IPIV( I ) = IPIV( I ) + J - 1
  160          CONTINUE
            END IF
!
!           Partially undo the interchanges in the current block to
!           restore the upper triangular form of A31 and copy the upper
!           triangle of A31 back into place
!
            DO 170 JJ = J + JB - 1, J, -1
               JP = IPIV( JJ ) - JJ + 1
               IF( JP.NE.1 ) THEN
!
!                 Apply interchange to columns J to JJ-1
!
                  IF( JP+JJ-1.LT.J+KL ) THEN
!
!                    The interchange does not affect A31
!
                     CALL ZSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1,      &
     &                           AB( KV+JP+JJ-J, J ), LDAB-1 )
                  ELSE
!
!                    The interchange does affect A31
!
                     CALL ZSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1,      &
     &                           WORK31( JP+JJ-J-KL, 1 ), LDWORK )
                  END IF
               END IF
!
!              Copy the current column of A31 back into place
!
               NW = MIN( I3, JJ-J+1 )
               IF( NW.GT.0 )                                            &
     &            CALL ZCOPY( NW, WORK31( 1, JJ-J+1 ), 1,               &
     &                        AB( KV+KL+1-JJ+J, JJ ), 1 )
  170       CONTINUE
  180    CONTINUE
      END IF
!
      RETURN
!
!     End of ZGBTRF
!
      END
      SUBROUTINE ZGETF2( M, N, A, LDA, IPIV, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  ZGETF2 computes an LU factorization of a general m-by-n matrix A
!  using partial pivoting with row interchanges.
!
!  The factorization has the form
!     A = P * L * U
!  where P is a permutation matrix, L is lower triangular with unit
!  diagonal elements (lower trapezoidal if m > n), and U is upper
!  triangular (upper trapezoidal if m < n).
!
!  This is the right-looking Level 2 BLAS version of the algorithm.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
!          On entry, the m by n matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value
!          > 0: if INFO = k, U(k,k) is exactly zero. The factorization
!               has been completed, but the factor U is exactly
!               singular, and division by zero will occur if it is used
!               to solve a system of equations.
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ),                    &
     &                   ZERO = ( 0.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      INTEGER            J, JP
!     ..
!     .. External Functions ..
      INTEGER            IZAMAX
      EXTERNAL           IZAMAX
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZGERU, ZSCAL, ZSWAP
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGETF2', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 )                                          &
     &   RETURN
!
      DO 10 J = 1, MIN( M, N )
!
!        Find pivot and test for singularity.
!
         JP = J - 1 + IZAMAX( M-J+1, A( J, J ), 1 )
         IPIV( J ) = JP
         IF( A( JP, J ).NE.ZERO ) THEN
!
!           Apply the interchange to columns 1:N.
!
            IF( JP.NE.J )                                               &
     &         CALL ZSWAP( N, A( J, 1 ), LDA, A( JP, 1 ), LDA )
!
!           Compute elements J+1:M of J-th column.
!
            IF( J.LT.M )                                                &
     &         CALL ZSCAL( M-J, ONE / A( J, J ), A( J+1, J ), 1 )
!
         ELSE IF( INFO.EQ.0 ) THEN
!
            INFO = J
         END IF
!
         IF( J.LT.MIN( M, N ) ) THEN
!
!           Update trailing submatrix.
!
            CALL ZGERU( M-J, N-J, -ONE, A( J+1, J ), 1, A( J, J+1 ),    &
     &                  LDA, A( J+1, J+1 ), LDA )
         END IF
   10 CONTINUE
      RETURN
!
!     End of ZGETF2
!
      END
      SUBROUTINE ZGETRF( M, N, A, LDA, IPIV, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  ZGETRF computes an LU factorization of a general M-by-N matrix A
!  using partial pivoting with row interchanges.
!
!  The factorization has the form
!     A = P * L * U
!  where P is a permutation matrix, L is lower triangular with unit
!  diagonal elements (lower trapezoidal if m > n), and U is upper
!  triangular (upper trapezoidal if m < n).
!
!  This is the right-looking Level 3 BLAS version of the algorithm.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
!          On entry, the M-by-N matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
!                has been completed, but the factor U is exactly
!                singular, and division by zero will occur if it is used
!                to solve a system of equations.
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      INTEGER            I, IINFO, J, JB, NB
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZGEMM, ZGETF2, ZLASWP, ZTRSM
!     ..
!     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGETRF', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 )                                          &
     &   RETURN
!
!     Determine the block size for this environment.
!
      NB = ILAENV( 1, 'ZGETRF', ' ', M, N, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) THEN
!
!        Use unblocked code.
!
         CALL ZGETF2( M, N, A, LDA, IPIV, INFO )
      ELSE
!
!        Use blocked code.
!
         DO 20 J = 1, MIN( M, N ), NB
            JB = MIN( MIN( M, N )-J+1, NB )
!
!           Factor diagonal and subdiagonal blocks and test for exact
!           singularity.
!
            CALL ZGETF2( M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO )
!
!           Adjust INFO and the pivot indices.
!
            IF( INFO.EQ.0 .AND. IINFO.GT.0 )                            &
     &         INFO = IINFO + J - 1
            DO 10 I = J, MIN( M, J+JB-1 )
               IPIV( I ) = J - 1 + IPIV( I )
   10       CONTINUE
!
!           Apply interchanges to columns 1:J-1.
!
            CALL ZLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )
!
            IF( J+JB.LE.N ) THEN
!
!              Apply interchanges to columns J+JB:N.
!
               CALL ZLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1,     &
     &                      IPIV, 1 )
!
!              Compute block row of U.
!
               CALL ZTRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB, &
     &                     N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ), &
     &                     LDA )
               IF( J+JB.LE.M ) THEN
!
!                 Update trailing submatrix.
!
                  CALL ZGEMM( 'No transpose', 'No transpose', M-J-JB+1, &
     &                        N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA,    &
     &                        A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ),  &
     &                        LDA )
               END IF
            END IF
   20    CONTINUE
      END IF
      RETURN
!
!     End of ZGETRF
!
      END
      SUBROUTINE ZGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, LDA, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  ZGETRS solves a system of linear equations
!     A * X = B,  A**T * X = B,  or  A**H * X = B
!  with a general N-by-N matrix A using the LU factorization computed
!  by ZGETRF.
!
!  Arguments
!  =========
!
!  TRANS   (input) CHARACTER*1
!          Specifies the form of the system of equations:
!          = 'N':  A * X = B     (No transpose)
!          = 'T':  A**T * X = B  (Transpose)
!          = 'C':  A**H * X = B  (Conjugate transpose)
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  A       (input) COMPLEX*16 array, dimension (LDA,N)
!          The factors L and U from the factorization A = P*L*U
!          as computed by ZGETRF.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  IPIV    (input) INTEGER array, dimension (N)
!          The pivot indices from ZGETRF; for 1<=i<=N, row i of the
!          matrix was interchanged with row IPIV(i).
!
!  B       (input/output) COMPLEX*16 array, dimension (LDB,NRHS)
!          On entry, the right hand side matrix B.
!          On exit, the solution matrix X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
!     ..
!     .. Local Scalars ..
      LOGICAL            NOTRAN
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLASWP, ZTRSM
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.        &
     &    LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGETRS', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 .OR. NRHS.EQ.0 )                                       &
     &   RETURN
!
      IF( NOTRAN ) THEN
!
!        Solve A * X = B.
!
!        Apply row interchanges to the right hand sides.
!
         CALL ZLASWP( NRHS, B, LDB, 1, N, IPIV, 1 )
!
!        Solve L*X = B, overwriting B with X.
!
         CALL ZTRSM( 'Left', 'Lower', 'No transpose', 'Unit', N, NRHS,  &
     &               ONE, A, LDA, B, LDB )
!
!        Solve U*X = B, overwriting B with X.
!
         CALL ZTRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N,    &
     &               NRHS, ONE, A, LDA, B, LDB )
      ELSE
!
!        Solve A**T * X = B  or A**H * X = B.
!
!        Solve U'*X = B, overwriting B with X.
!
         CALL ZTRSM( 'Left', 'Upper', TRANS, 'Non-unit', N, NRHS, ONE,  &
     &               A, LDA, B, LDB )
!
!        Solve L'*X = B, overwriting B with X.
!
         CALL ZTRSM( 'Left', 'Lower', TRANS, 'Unit', N, NRHS, ONE, A,   &
     &               LDA, B, LDB )
!
!        Apply row interchanges to the solution vectors.
!
         CALL ZLASWP( NRHS, B, LDB, 1, N, IPIV, -1 )
      END IF
!
      RETURN
!
!     End of ZGETRS
!
      END
      SUBROUTINE ZLASWP( N, A, LDA, K1, K2, IPIV, INCX )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      INTEGER            INCX, K1, K2, LDA, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  ZLASWP performs a series of row interchanges on the matrix A.
!  One row interchange is initiated for each of rows K1 through K2 of A.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.
!
!  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
!          On entry, the matrix of column dimension N to which the row
!          interchanges will be applied.
!          On exit, the permuted matrix.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.
!
!  K1      (input) INTEGER
!          The first element of IPIV for which a row interchange will
!          be done.
!
!  K2      (input) INTEGER
!          The last element of IPIV for which a row interchange will
!          be done.
!
!  IPIV    (input) INTEGER array, dimension (M*abs(INCX))
!          The vector of pivot indices.  Only the elements in positions
!          K1 through K2 of IPIV are accessed.
!          IPIV(K) = L implies rows K and L are to be interchanged.
!
!  INCX    (input) INTEGER
!          The increment between successive values of IPIV.  If IPIV
!          is negative, the pivots are applied in reverse order.
!
!  Further Details
!  ===============
!
!  Modified by
!   R. C. Whaley, Computer Science Dept., Univ. of Tenn., Knoxville, USA
!
! =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I, I1, I2, INC, IP, IX, IX0, J, K, N32
      COMPLEX*16         TEMP
!     ..
!     .. Executable Statements ..
!
!     Interchange row I with row IPIV(I) for each of rows K1 through K2.
!
      IF( INCX.GT.0 ) THEN
         IX0 = K1
         I1 = K1
         I2 = K2
         INC = 1
      ELSE IF( INCX.LT.0 ) THEN
         IX0 = 1 + ( 1-K2 )*INCX
         I1 = K2
         I2 = K1
         INC = -1
      ELSE
         RETURN
      END IF
!
      N32 = ( N / 32 )*32
      IF( N32.NE.0 ) THEN
         DO 30 J = 1, N32, 32
            IX = IX0
            DO 20 I = I1, I2, INC
               IP = IPIV( IX )
               IF( IP.NE.I ) THEN
                  DO 10 K = J, J + 31
                     TEMP = A( I, K )
                     A( I, K ) = A( IP, K )
                     A( IP, K ) = TEMP
   10             CONTINUE
               END IF
               IX = IX + INCX
   20       CONTINUE
   30    CONTINUE
      END IF
      IF( N32.NE.N ) THEN
         N32 = N32 + 1
         IX = IX0
         DO 50 I = I1, I2, INC
            IP = IPIV( IX )
            IF( IP.NE.I ) THEN
               DO 40 K = N32, N
                  TEMP = A( I, K )
                  A( I, K ) = A( IP, K )
                  A( IP, K ) = TEMP
   40          CONTINUE
            END IF
            IX = IX + INCX
   50    CONTINUE
      END IF
!
      RETURN
!
!     End of ZLASWP
!
      END
      SUBROUTINE DGETF2( M, N, A, LDA, IPIV, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1992
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  DGETF2 computes an LU factorization of a general m-by-n matrix A
!  using partial pivoting with row interchanges.
!
!  The factorization has the form
!     A = P * L * U
!  where P is a permutation matrix, L is lower triangular with unit
!  diagonal elements (lower trapezoidal if m > n), and U is upper
!  triangular (upper trapezoidal if m < n).
!
!  This is the right-looking Level 2 BLAS version of the algorithm.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the m by n matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value
!          > 0: if INFO = k, U(k,k) is exactly zero. The factorization
!               has been completed, but the factor U is exactly
!               singular, and division by zero will occur if it is used
!               to solve a system of equations.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            J, JP
!     ..
!     .. External Functions ..
      INTEGER            IDAMAX
      EXTERNAL           IDAMAX
!     ..
!     .. External Subroutines ..
      EXTERNAL           DGER, DSCAL, DSWAP, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETF2', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 )                                          &
     &   RETURN
!
      DO 10 J = 1, MIN( M, N )
!
!        Find pivot and test for singularity.
!
         JP = J - 1 + IDAMAX( M-J+1, A( J, J ), 1 )
         IPIV( J ) = JP
         IF( A( JP, J ).NE.ZERO ) THEN
!
!           Apply the interchange to columns 1:N.
!
            IF( JP.NE.J )                                               &
     &         CALL DSWAP( N, A( J, 1 ), LDA, A( JP, 1 ), LDA )
!
!           Compute elements J+1:M of J-th column.
!
            IF( J.LT.M )                                                &
     &         CALL DSCAL( M-J, ONE / A( J, J ), A( J+1, J ), 1 )
!
         ELSE IF( INFO.EQ.0 ) THEN
!
            INFO = J
         END IF
!
         IF( J.LT.MIN( M, N ) ) THEN
!
!           Update trailing submatrix.
!
            CALL DGER( M-J, N-J, -ONE, A( J+1, J ), 1, A( J, J+1 ), LDA,&
     &                 A( J+1, J+1 ), LDA )
         END IF
   10 CONTINUE
      RETURN
!
!     End of DGETF2
!
      END
      SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  DGETRF computes an LU factorization of a general M-by-N matrix A
!  using partial pivoting with row interchanges.
!
!  The factorization has the form
!     A = P * L * U
!  where P is a permutation matrix, L is lower triangular with unit
!  diagonal elements (lower trapezoidal if m > n), and U is upper
!  triangular (upper trapezoidal if m < n).
!
!  This is the right-looking Level 3 BLAS version of the algorithm.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the M-by-N matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
!                has been completed, but the factor U is exactly
!                singular, and division by zero will occur if it is used
!                to solve a system of equations.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, IINFO, J, JB, NB
!     ..
!     .. External Subroutines ..
      EXTERNAL           DGEMM, DGETF2, DLASWP, DTRSM, XERBLA
!     ..
!     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETRF', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 )                                          &
     &   RETURN
!
!     Determine the block size for this environment.
!
      NB = ILAENV( 1, 'DGETRF', ' ', M, N, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) THEN
!
!        Use unblocked code.
!
         CALL DGETF2( M, N, A, LDA, IPIV, INFO )
      ELSE
!
!        Use blocked code.
!
         DO 20 J = 1, MIN( M, N ), NB
            JB = MIN( MIN( M, N )-J+1, NB )
!
!           Factor diagonal and subdiagonal blocks and test for exact
!           singularity.
!
            CALL DGETF2( M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO )
!
!           Adjust INFO and the pivot indices.
!
            IF( INFO.EQ.0 .AND. IINFO.GT.0 )                            &
     &         INFO = IINFO + J - 1
            DO 10 I = J, MIN( M, J+JB-1 )
               IPIV( I ) = J - 1 + IPIV( I )
   10       CONTINUE
!
!           Apply interchanges to columns 1:J-1.
!
            CALL DLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )
!
            IF( J+JB.LE.N ) THEN
!
!              Apply interchanges to columns J+JB:N.
!
               CALL DLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1,     &
     &                      IPIV, 1 )
!
!              Compute block row of U.
!
               CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB, &
     &                     N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ), &
     &                     LDA )
               IF( J+JB.LE.M ) THEN
!
!                 Update trailing submatrix.
!
                  CALL DGEMM( 'No transpose', 'No transpose', M-J-JB+1, &
     &                        N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA,    &
     &                        A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ),  &
     &                        LDA )
               END IF
            END IF
   20    CONTINUE
      END IF
      RETURN
!
!     End of DGETRF
!
      END
      SUBROUTINE DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, LDA, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  DGETRS solves a system of linear equations
!     A * X = B  or  A' * X = B
!  with a general N-by-N matrix A using the LU factorization computed
!  by DGETRF.
!
!  Arguments
!  =========
!
!  TRANS   (input) CHARACTER*1
!          Specifies the form of the system of equations:
!          = 'N':  A * X = B  (No transpose)
!          = 'T':  A'* X = B  (Transpose)
!          = 'C':  A'* X = B  (Conjugate transpose = Transpose)
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
!          The factors L and U from the factorization A = P*L*U
!          as computed by DGETRF.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  IPIV    (input) INTEGER array, dimension (N)
!          The pivot indices from DGETRF; for 1<=i<=N, row i of the
!          matrix was interchanged with row IPIV(i).
!
!  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
!          On entry, the right hand side matrix B.
!          On exit, the solution matrix X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            NOTRAN
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           DLASWP, DTRSM, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.        &
     &    LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETRS', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 .OR. NRHS.EQ.0 )                                       &
     &   RETURN
!
      IF( NOTRAN ) THEN
!
!        Solve A * X = B.
!
!        Apply row interchanges to the right hand sides.
!
         CALL DLASWP( NRHS, B, LDB, 1, N, IPIV, 1 )
!
!        Solve L*X = B, overwriting B with X.
!
         CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', N, NRHS,  &
     &               ONE, A, LDA, B, LDB )
!
!        Solve U*X = B, overwriting B with X.
!
         CALL DTRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N,    &
     &               NRHS, ONE, A, LDA, B, LDB )
      ELSE
!
!        Solve A' * X = B.
!
!        Solve U'*X = B, overwriting B with X.
!
         CALL DTRSM( 'Left', 'Upper', 'Transpose', 'Non-unit', N, NRHS, &
     &               ONE, A, LDA, B, LDB )
!
!        Solve L'*X = B, overwriting B with X.
!
         CALL DTRSM( 'Left', 'Lower', 'Transpose', 'Unit', N, NRHS, ONE,&
     &               A, LDA, B, LDB )
!
!        Apply row interchanges to the solution vectors.
!
         CALL DLASWP( NRHS, B, LDB, 1, N, IPIV, -1 )
      END IF
!
      RETURN
!
!     End of DGETRS
!
      END
      SUBROUTINE DLASWP( N, A, LDA, K1, K2, IPIV, INCX )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      INTEGER            INCX, K1, K2, LDA, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  DLASWP performs a series of row interchanges on the matrix A.
!  One row interchange is initiated for each of rows K1 through K2 of A.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the matrix of column dimension N to which the row
!          interchanges will be applied.
!          On exit, the permuted matrix.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.
!
!  K1      (input) INTEGER
!          The first element of IPIV for which a row interchange will
!          be done.
!
!  K2      (input) INTEGER
!          The last element of IPIV for which a row interchange will
!          be done.
!
!  IPIV    (input) INTEGER array, dimension (M*abs(INCX))
!          The vector of pivot indices.  Only the elements in positions
!          K1 through K2 of IPIV are accessed.
!          IPIV(K) = L implies rows K and L are to be interchanged.
!
!  INCX    (input) INTEGER
!          The increment between successive values of IPIV.  If IPIV
!          is negative, the pivots are applied in reverse order.
!
!  Further Details
!  ===============
!
!  Modified by
!   R. C. Whaley, Computer Science Dept., Univ. of Tenn., Knoxville, USA
!
! =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I, I1, I2, INC, IP, IX, IX0, J, K, N32
      DOUBLE PRECISION   TEMP
!     ..
!     .. Executable Statements ..
!
!     Interchange row I with row IPIV(I) for each of rows K1 through K2.
!
      IF( INCX.GT.0 ) THEN
         IX0 = K1
         I1 = K1
         I2 = K2
         INC = 1
      ELSE IF( INCX.LT.0 ) THEN
         IX0 = 1 + ( 1-K2 )*INCX
         I1 = K2
         I2 = K1
         INC = -1
      ELSE
         RETURN
      END IF
!
      N32 = ( N / 32 )*32
      IF( N32.NE.0 ) THEN
         DO 30 J = 1, N32, 32
            IX = IX0
            DO 20 I = I1, I2, INC
               IP = IPIV( IX )
               IF( IP.NE.I ) THEN
                  DO 10 K = J, J + 31
                     TEMP = A( I, K )
                     A( I, K ) = A( IP, K )
                     A( IP, K ) = TEMP
   10             CONTINUE
               END IF
               IX = IX + INCX
   20       CONTINUE
   30    CONTINUE
      END IF
      IF( N32.NE.N ) THEN
         N32 = N32 + 1
         IX = IX0
         DO 50 I = I1, I2, INC
            IP = IPIV( IX )
            IF( IP.NE.I ) THEN
               DO 40 K = N32, N
                  TEMP = A( I, K )
                  A( I, K ) = A( IP, K )
                  A( IP, K ) = TEMP
   40          CONTINUE
            END IF
            IX = IX + INCX
   50    CONTINUE
      END IF
!
      RETURN
!
!     End of DLASWP
!
      END
      REAL             FUNCTION SCSUM1( N, CX, INCX )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      INTEGER            INCX, N
!     ..
!     .. Array Arguments ..
      COMPLEX            CX( * )
!     ..
!
!  Purpose
!  =======
!
!  SCSUM1 takes the sum of the absolute values of a complex
!  vector and returns a single precision result.
!
!  Based on SCASUM from the Level 1 BLAS.
!  The change is to use the 'genuine' absolute value.
!
!  Contributed by Nick Higham for use with CLACON.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The number of elements in the vector CX.
!
!  CX      (input) COMPLEX array, dimension (N)
!          The vector whose elements will be summed.
!
!  INCX    (input) INTEGER
!          The spacing between successive values of CX.  INCX > 0.
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I, NINCX
      REAL               STEMP
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS
!     ..
!     .. Executable Statements ..
!
      SCSUM1 = 0.0E0
      STEMP = 0.0E0
      IF( N.LE.0 )                                                      &
     &   RETURN
      IF( INCX.EQ.1 )                                                   &
     &   GO TO 20
!
!     CODE FOR INCREMENT NOT EQUAL TO 1
!
      NINCX = N*INCX
      DO 10 I = 1, NINCX, INCX
!
!        NEXT LINE MODIFIED.
!
         STEMP = STEMP + ABS( CX( I ) )
   10 CONTINUE
      SCSUM1 = STEMP
      RETURN
!
!     CODE FOR INCREMENT EQUAL TO 1
!
   20 CONTINUE
      DO 30 I = 1, N
!
!        NEXT LINE MODIFIED.
!
         STEMP = STEMP + ABS( CX( I ) )
   30 CONTINUE
      SCSUM1 = STEMP
      RETURN
!
!     End of SCSUM1
!
      END
      INTEGER          FUNCTION ICMAX1( N, CX, INCX )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      INTEGER            INCX, N
!     ..
!     .. Array Arguments ..
      COMPLEX            CX( * )
!     ..
!
!  Purpose
!  =======
!
!  ICMAX1 finds the index of the element whose real part has maximum
!  absolute value.
!
!  Based on ICAMAX from Level 1 BLAS.
!  The change is to use the 'genuine' absolute value.
!
!  Contributed by Nick Higham for use with CLACON.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The number of elements in the vector CX.
!
!  CX      (input) COMPLEX array, dimension (N)
!          The vector whose elements will be summed.
!
!  INCX    (input) INTEGER
!          The spacing between successive values of CX.  INCX >= 1.
!
! =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I, IX
      REAL               SMAX
      COMPLEX            ZDUM
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, REAL
!     ..
!     .. Statement Functions ..
      REAL               CABS1
!     ..
!     .. Statement Function definitions ..
!
!     NEXT LINE IS THE ONLY MODIFICATION.
      CABS1( ZDUM ) = ABS( REAL( ZDUM ) )
!     ..
!     .. Executable Statements ..
!
      ICMAX1 = 0
      IF( N.LT.1 )                                                      &
     &   RETURN
      ICMAX1 = 1
      IF( N.EQ.1 )                                                      &
     &   RETURN
      IF( INCX.EQ.1 )                                                   &
     &   GO TO 30
!
!     CODE FOR INCREMENT NOT EQUAL TO 1
!
      IX = 1
      SMAX = CABS1( CX( 1 ) )
      IX = IX + INCX
      DO 20 I = 2, N
         IF( CABS1( CX( IX ) ).LE.SMAX )                                &
     &      GO TO 10
         ICMAX1 = I
         SMAX = CABS1( CX( IX ) )
   10    CONTINUE
         IX = IX + INCX
   20 CONTINUE
      RETURN
!
!     CODE FOR INCREMENT EQUAL TO 1
!
   30 CONTINUE
      SMAX = CABS1( CX( 1 ) )
      DO 40 I = 2, N
         IF( CABS1( CX( I ) ).LE.SMAX )                                 &
     &      GO TO 40
         ICMAX1 = I
         SMAX = CABS1( CX( I ) )
   40 CONTINUE
      RETURN
!
!     End of ICMAX1
!
      END
      INTEGER          FUNCTION IEEECK( ISPEC, ZERO, ONE )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1998
!
!     .. Scalar Arguments ..
      INTEGER            ISPEC
      REAL               ONE, ZERO
!     ..
!
!  Purpose
!  =======
!
!  IEEECK is called from the ILAENV to verify that Infinity and
!  possibly NaN arithmetic is safe (i.e. will not trap).
!
!  Arguments
!  =========
!
!  ISPEC   (input) INTEGER
!          Specifies whether to test just for inifinity arithmetic
!          or whether to test for infinity and NaN arithmetic.
!          = 0: Verify infinity arithmetic only.
!          = 1: Verify infinity and NaN arithmetic.
!
!  ZERO    (input) REAL
!          Must contain the value 0.0
!          This is passed to prevent the compiler from optimizing
!          away this code.
!
!  ONE     (input) REAL
!          Must contain the value 1.0
!          This is passed to prevent the compiler from optimizing
!          away this code.
!
!  RETURN VALUE:  INTEGER
!          = 0:  Arithmetic failed to produce the correct answers
!          = 1:  Arithmetic produced the correct answers
!
!     .. Local Scalars ..
      REAL               NAN1, NAN2, NAN3, NAN4, NAN5, NAN6, NEGINF,    &
     &                   NEGZRO, NEWZRO, POSINF
!     ..
!     .. Executable Statements ..
      IEEECK = 1
!
      POSINF = ONE / ZERO
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      NEGINF = -ONE / ZERO
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      NEGZRO = ONE / ( NEGINF+ONE )
      IF( NEGZRO.NE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      NEGINF = ONE / NEGZRO
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      NEWZRO = NEGZRO + ZERO
      IF( NEWZRO.NE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      POSINF = ONE / NEWZRO
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      NEGINF = NEGINF*POSINF
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      POSINF = POSINF*POSINF
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
!
!
!
!
!     Return if we were only asked to check infinity arithmetic
!
      IF( ISPEC.EQ.0 )                                                  &
     &   RETURN
!
      NAN1 = POSINF + NEGINF
!
      NAN2 = POSINF / NEGINF
!
      NAN3 = POSINF / POSINF
!
      NAN4 = POSINF*ZERO
!
      NAN5 = NEGINF*NEGZRO
!
      NAN6 = NAN5*0.0
!
      IF( NAN1.EQ.NAN1 ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      IF( NAN2.EQ.NAN2 ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      IF( NAN3.EQ.NAN3 ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      IF( NAN4.EQ.NAN4 ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      IF( NAN5.EQ.NAN5 ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      IF( NAN6.EQ.NAN6 ) THEN
         IEEECK = 0
         RETURN
      END IF
!
      RETURN
      END
      REAL             FUNCTION SLAMCH( CMACH )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992 
!
!     .. Scalar Arguments ..
      CHARACTER          CMACH
!     ..
!
!  Purpose
!  =======
!
!  SLAMCH determines single precision machine parameters.
!
!  Arguments
!  =========
!
!  CMACH   (input) CHARACTER*1
!          Specifies the value to be returned by SLAMCH:
!          = 'E' or 'e',   SLAMCH := eps
!          = 'S' or 's ,   SLAMCH := sfmin
!          = 'B' or 'b',   SLAMCH := base
!          = 'P' or 'p',   SLAMCH := eps*base
!          = 'N' or 'n',   SLAMCH := t
!          = 'R' or 'r',   SLAMCH := rnd
!          = 'M' or 'm',   SLAMCH := emin
!          = 'U' or 'u',   SLAMCH := rmin
!          = 'L' or 'l',   SLAMCH := emax
!          = 'O' or 'o',   SLAMCH := rmax
!
!          where
!
!          eps   = relative machine precision
!          sfmin = safe minimum, such that 1/sfmin does not overflow
!          base  = base of the machine
!          prec  = eps*base
!          t     = number of (base) digits in the mantissa
!          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
!          emin  = minimum exponent before (gradual) underflow
!          rmin  = underflow threshold - base**(emin-1)
!          emax  = largest exponent before overflow
!          rmax  = overflow threshold  - (base**emax)*(1-eps)
!
! =====================================================================
!
!     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            FIRST, LRND
      INTEGER            BETA, IMAX, IMIN, IT
      REAL               BASE, EMAX, EMIN, EPS, PREC, RMACH, RMAX, RMIN,&
     &                   RND, SFMIN, SMALL, T
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           SLAMC2
!     ..
!     .. Save statement ..
      SAVE               FIRST, EPS, SFMIN, BASE, T, RND, EMIN, RMIN,   &
     &                   EMAX, RMAX, PREC
!     ..
!     .. Data statements ..
      DATA               FIRST / .TRUE. /
!     ..
!     .. Executable Statements ..
!
      IF( FIRST ) THEN
         FIRST = .FALSE.
         CALL SLAMC2( BETA, IT, LRND, EPS, IMIN, RMIN, IMAX, RMAX )
         BASE = BETA
         T = IT
         IF( LRND ) THEN
            RND = ONE
            EPS = ( BASE**( 1-IT ) ) / 2
         ELSE
            RND = ZERO
            EPS = BASE**( 1-IT )
         END IF
         PREC = EPS*BASE
         EMIN = IMIN
         EMAX = IMAX
         SFMIN = RMIN
         SMALL = ONE / RMAX
         IF( SMALL.GE.SFMIN ) THEN
!
!           Use SMALL plus a bit, to avoid the possibility of rounding
!           causing overflow when computing  1/sfmin.
!
            SFMIN = SMALL*( ONE+EPS )
         END IF
      END IF
!
      IF( LSAME( CMACH, 'E' ) ) THEN
         RMACH = EPS
      ELSE IF( LSAME( CMACH, 'S' ) ) THEN
         RMACH = SFMIN
      ELSE IF( LSAME( CMACH, 'B' ) ) THEN
         RMACH = BASE
      ELSE IF( LSAME( CMACH, 'P' ) ) THEN
         RMACH = PREC
      ELSE IF( LSAME( CMACH, 'N' ) ) THEN
         RMACH = T
      ELSE IF( LSAME( CMACH, 'R' ) ) THEN
         RMACH = RND
      ELSE IF( LSAME( CMACH, 'M' ) ) THEN
         RMACH = EMIN
      ELSE IF( LSAME( CMACH, 'U' ) ) THEN
         RMACH = RMIN
      ELSE IF( LSAME( CMACH, 'L' ) ) THEN
         RMACH = EMAX
      ELSE IF( LSAME( CMACH, 'O' ) ) THEN
         RMACH = RMAX
      END IF
!
      SLAMCH = RMACH
      RETURN
!
!     End of SLAMCH
!
      END
!
!***********************************************************************
!
      SUBROUTINE SLAMC1( BETA, T, RND, IEEE1 )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      LOGICAL            IEEE1, RND
      INTEGER            BETA, T
!     ..
!
!  Purpose
!  =======
!
!  SLAMC1 determines the machine parameters given by BETA, T, RND, and
!  IEEE1.
!
!  Arguments
!  =========
!
!  BETA    (output) INTEGER
!          The base of the machine.
!
!  T       (output) INTEGER
!          The number of ( BETA ) digits in the mantissa.
!
!  RND     (output) LOGICAL
!          Specifies whether proper rounding  ( RND = .TRUE. )  or
!          chopping  ( RND = .FALSE. )  occurs in addition. This may not
!          be a reliable guide to the way in which the machine performs
!          its arithmetic.
!
!  IEEE1   (output) LOGICAL
!          Specifies whether rounding appears to be done in the IEEE
!          'round to nearest' style.
!
!  Further Details
!  ===============
!
!  The routine is based on the routine  ENVRON  by Malcolm and
!  incorporates suggestions by Gentleman and Marovich. See
!
!     Malcolm M. A. (1972) Algorithms to reveal properties of
!        floating-point arithmetic. Comms. of the ACM, 15, 949-951.
!
!     Gentleman W. M. and Marovich S. B. (1974) More on algorithms
!        that reveal properties of floating point arithmetic units.
!        Comms. of the ACM, 17, 276-277.
!
! =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            FIRST, LIEEE1, LRND
      INTEGER            LBETA, LT
      REAL               A, B, C, F, ONE, QTR, SAVEC, T1, T2
!     ..
!     .. External Functions ..
      REAL               SLAMC3
      EXTERNAL           SLAMC3
!     ..
!     .. Save statement ..
      SAVE               FIRST, LIEEE1, LBETA, LRND, LT
!     ..
!     .. Data statements ..
      DATA               FIRST / .TRUE. /
!     ..
!     .. Executable Statements ..
!
      IF( FIRST ) THEN
         FIRST = .FALSE.
         ONE = 1
!
!        LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BETA,
!        IEEE1, T and RND.
!
!        Throughout this routine  we use the function  SLAMC3  to ensure
!        that relevant values are  stored and not held in registers,  or
!        are not affected by optimizers.
!
!        Compute  a = 2.0**m  with the  smallest positive integer m such
!        that
!
!           fl( a + 1.0 ) = a.
!
         A = 1
         C = 1
!
!+       WHILE( C.EQ.ONE )LOOP
   10    CONTINUE
         IF( C.EQ.ONE ) THEN
            A = 2*A
            C = SLAMC3( A, ONE )
            C = SLAMC3( C, -A )
            GO TO 10
         END IF
!+       END WHILE
!
!        Now compute  b = 2.0**m  with the smallest positive integer m
!        such that
!
!           fl( a + b ) .gt. a.
!
         B = 1
         C = SLAMC3( A, B )
!
!+       WHILE( C.EQ.A )LOOP
   20    CONTINUE
         IF( C.EQ.A ) THEN
            B = 2*B
            C = SLAMC3( A, B )
            GO TO 20
         END IF
!+       END WHILE
!
!        Now compute the base.  a and c  are neighbouring floating point
!        numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and so
!        their difference is beta. Adding 0.25 to c is to ensure that it
!        is truncated to beta and not ( beta - 1 ).
!
         QTR = ONE / 4
         SAVEC = C
         C = SLAMC3( C, -A )
         LBETA = C + QTR
!
!        Now determine whether rounding or chopping occurs,  by adding a
!        bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a.
!
         B = LBETA
         F = SLAMC3( B / 2, -B / 100 )
         C = SLAMC3( F, A )
         IF( C.EQ.A ) THEN
            LRND = .TRUE.
         ELSE
            LRND = .FALSE.
         END IF
         F = SLAMC3( B / 2, B / 100 )
         C = SLAMC3( F, A )
         IF( ( LRND ) .AND. ( C.EQ.A ) )                                &
     &      LRND = .FALSE.
!
!        Try and decide whether rounding is done in the  IEEE  'round to
!        nearest' style. B/2 is half a unit in the last place of the two
!        numbers A and SAVEC. Furthermore, A is even, i.e. has last  bit
!        zero, and SAVEC is odd. Thus adding B/2 to A should not  change
!        A, but adding B/2 to SAVEC should change SAVEC.
!
         T1 = SLAMC3( B / 2, A )
         T2 = SLAMC3( B / 2, SAVEC )
         LIEEE1 = ( T1.EQ.A ) .AND. ( T2.GT.SAVEC ) .AND. LRND
!
!        Now find  the  mantissa, t.  It should  be the  integer part of
!        log to the base beta of a,  however it is safer to determine  t
!        by powering.  So we find t as the smallest positive integer for
!        which
!
!           fl( beta**t + 1.0 ) = 1.0.
!
         LT = 0
         A = 1
         C = 1
!
!+       WHILE( C.EQ.ONE )LOOP
   30    CONTINUE
         IF( C.EQ.ONE ) THEN
            LT = LT + 1
            A = A*LBETA
            C = SLAMC3( A, ONE )
            C = SLAMC3( C, -A )
            GO TO 30
         END IF
!+       END WHILE
!
      END IF
!
      BETA = LBETA
      T = LT
      RND = LRND
      IEEE1 = LIEEE1
      RETURN
!
!     End of SLAMC1
!
      END
!
!***********************************************************************
!
      SUBROUTINE SLAMC2( BETA, T, RND, EPS, EMIN, RMIN, EMAX, RMAX )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      LOGICAL            RND
      INTEGER            BETA, EMAX, EMIN, T
      REAL               EPS, RMAX, RMIN
!     ..
!
!  Purpose
!  =======
!
!  SLAMC2 determines the machine parameters specified in its argument
!  list.
!
!  Arguments
!  =========
!
!  BETA    (output) INTEGER
!          The base of the machine.
!
!  T       (output) INTEGER
!          The number of ( BETA ) digits in the mantissa.
!
!  RND     (output) LOGICAL
!          Specifies whether proper rounding  ( RND = .TRUE. )  or
!          chopping  ( RND = .FALSE. )  occurs in addition. This may not
!          be a reliable guide to the way in which the machine performs
!          its arithmetic.
!
!  EPS     (output) REAL
!          The smallest positive number such that
!
!             fl( 1.0 - EPS ) .LT. 1.0,
!
!          where fl denotes the computed value.
!
!  EMIN    (output) INTEGER
!          The minimum exponent before (gradual) underflow occurs.
!
!  RMIN    (output) REAL
!          The smallest normalized number for the machine, given by
!          BASE**( EMIN - 1 ), where  BASE  is the floating point value
!          of BETA.
!
!  EMAX    (output) INTEGER
!          The maximum exponent before overflow occurs.
!
!  RMAX    (output) REAL
!          The largest positive number for the machine, given by
!          BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating point
!          value of BETA.
!
!  Further Details
!  ===============
!
!  The computation of  EPS  is based on a routine PARANOIA by
!  W. Kahan of the University of California at Berkeley.
!
! =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            FIRST, IEEE, IWARN, LIEEE1, LRND
      INTEGER            GNMIN, GPMIN, I, LBETA, LEMAX, LEMIN, LT,      &
     &                   NGNMIN, NGPMIN
      REAL               A, B, C, HALF, LEPS, LRMAX, LRMIN, ONE, RBASE, &
     &                   SIXTH, SMALL, THIRD, TWO, ZERO
!     ..
!     .. External Functions ..
      REAL               SLAMC3
      EXTERNAL           SLAMC3
!     ..
!     .. External Subroutines ..
      EXTERNAL           SLAMC1, SLAMC4, SLAMC5
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
!     ..
!     .. Save statement ..
      SAVE               FIRST, IWARN, LBETA, LEMAX, LEMIN, LEPS, LRMAX,&
     &                   LRMIN, LT
!     ..
!     .. Data statements ..
      DATA               FIRST / .TRUE. / , IWARN / .FALSE. /
!     ..
!     .. Executable Statements ..
!
      IF( FIRST ) THEN
         FIRST = .FALSE.
         ZERO = 0
         ONE = 1
         TWO = 2
!
!        LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values of
!        BETA, T, RND, EPS, EMIN and RMIN.
!
!        Throughout this routine  we use the function  SLAMC3  to ensure
!        that relevant values are stored  and not held in registers,  or
!        are not affected by optimizers.
!
!        SLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1.
!
         CALL SLAMC1( LBETA, LT, LRND, LIEEE1 )
!
!        Start to find EPS.
!
         B = LBETA
         A = B**( -LT )
         LEPS = A
!
!        Try some tricks to see whether or not this is the correct  EPS.
!
         B = TWO / 3
         HALF = ONE / 2
         SIXTH = SLAMC3( B, -HALF )
         THIRD = SLAMC3( SIXTH, SIXTH )
         B = SLAMC3( THIRD, -HALF )
         B = SLAMC3( B, SIXTH )
         B = ABS( B )
         IF( B.LT.LEPS )                                                &
     &      B = LEPS
!
         LEPS = 1
!
!+       WHILE( ( LEPS.GT.B ).AND.( B.GT.ZERO ) )LOOP
   10    CONTINUE
         IF( ( LEPS.GT.B ) .AND. ( B.GT.ZERO ) ) THEN
            LEPS = B
            C = SLAMC3( HALF*LEPS, ( TWO**5 )*( LEPS**2 ) )
            C = SLAMC3( HALF, -C )
            B = SLAMC3( HALF, C )
            C = SLAMC3( HALF, -B )
            B = SLAMC3( HALF, C )
            GO TO 10
         END IF
!+       END WHILE
!
         IF( A.LT.LEPS )                                                &
     &      LEPS = A
!
!        Computation of EPS complete.
!
!        Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3)).
!        Keep dividing  A by BETA until (gradual) underflow occurs. This
!        is detected when we cannot recover the previous A.
!
         RBASE = ONE / LBETA
         SMALL = ONE
         DO 20 I = 1, 3
            SMALL = SLAMC3( SMALL*RBASE, ZERO )
   20    CONTINUE
         A = SLAMC3( ONE, SMALL )
         CALL SLAMC4( NGPMIN, ONE, LBETA )
         CALL SLAMC4( NGNMIN, -ONE, LBETA )
         CALL SLAMC4( GPMIN, A, LBETA )
         CALL SLAMC4( GNMIN, -A, LBETA )
         IEEE = .FALSE.
!
         IF( ( NGPMIN.EQ.NGNMIN ) .AND. ( GPMIN.EQ.GNMIN ) ) THEN
            IF( NGPMIN.EQ.GPMIN ) THEN
               LEMIN = NGPMIN
!            ( Non twos-complement machines, no gradual underflow;
!              e.g.,  VAX )
            ELSE IF( ( GPMIN-NGPMIN ).EQ.3 ) THEN
               LEMIN = NGPMIN - 1 + LT
               IEEE = .TRUE.
!            ( Non twos-complement machines, with gradual underflow;
!              e.g., IEEE standard followers )
            ELSE
               LEMIN = MIN( NGPMIN, GPMIN )
!            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
!
         ELSE IF( ( NGPMIN.EQ.GPMIN ) .AND. ( NGNMIN.EQ.GNMIN ) ) THEN
            IF( ABS( NGPMIN-NGNMIN ).EQ.1 ) THEN
               LEMIN = MAX( NGPMIN, NGNMIN )
!            ( Twos-complement machines, no gradual underflow;
!              e.g., CYBER 205 )
            ELSE
               LEMIN = MIN( NGPMIN, NGNMIN )
!            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
!
         ELSE IF( ( ABS( NGPMIN-NGNMIN ).EQ.1 ) .AND.                   &
     &            ( GPMIN.EQ.GNMIN ) ) THEN
            IF( ( GPMIN-MIN( NGPMIN, NGNMIN ) ).EQ.3 ) THEN
               LEMIN = MAX( NGPMIN, NGNMIN ) - 1 + LT
!            ( Twos-complement machines with gradual underflow;
!              no known machine )
            ELSE
               LEMIN = MIN( NGPMIN, NGNMIN )
!            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
!
         ELSE
            LEMIN = MIN( NGPMIN, NGNMIN, GPMIN, GNMIN )
!         ( A guess; no known machine )
            IWARN = .TRUE.
         END IF
!**
! Comment out this if block if EMIN is ok
         IF( IWARN ) THEN
            FIRST = .TRUE.
            WRITE( 6, FMT = 9999 )LEMIN
         END IF
!**
!
!        Assume IEEE arithmetic if we found denormalised  numbers above,
!        or if arithmetic seems to round in the  IEEE style,  determined
!        in routine SLAMC1. A true IEEE machine should have both  things
!        true; however, faulty machines may have one or the other.
!
         IEEE = IEEE .OR. LIEEE1
!
!        Compute  RMIN by successive division by  BETA. We could compute
!        RMIN as BASE**( EMIN - 1 ),  but some machines underflow during
!        this computation.
!
         LRMIN = 1
         DO 30 I = 1, 1 - LEMIN
            LRMIN = SLAMC3( LRMIN*RBASE, ZERO )
   30    CONTINUE
!
!        Finally, call SLAMC5 to compute EMAX and RMAX.
!
         CALL SLAMC5( LBETA, LT, LEMIN, IEEE, LEMAX, LRMAX )
      END IF
!
      BETA = LBETA
      T = LT
      RND = LRND
      EPS = LEPS
      EMIN = LEMIN
      RMIN = LRMIN
      EMAX = LEMAX
      RMAX = LRMAX
!
      RETURN
!
 9999 FORMAT( / / ' WARNING. The value EMIN may be incorrect:-',        &
     &      '  EMIN = ', I8, /                                          &
     &      ' If, after inspection, the value EMIN looks',              &
     &      ' acceptable please comment out ',                          &
     &      / ' the IF block as marked within the code of routine',     &
     &      ' SLAMC2,', / ' otherwise supply EMIN explicitly.', / )
!
!     End of SLAMC2
!
      END
!
!***********************************************************************
!
      REAL             FUNCTION SLAMC3( A, B )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      REAL               A, B
!     ..
!
!  Purpose
!  =======
!
!  SLAMC3  is intended to force  A  and  B  to be stored prior to doing
!  the addition of  A  and  B ,  for use in situations where optimizers
!  might hold one of these in a register.
!
!  Arguments
!  =========
!
!  A, B    (input) REAL
!          The values A and B.
!
! =====================================================================
!
!     .. Executable Statements ..
!
      SLAMC3 = A + B
!
      RETURN
!
!     End of SLAMC3
!
      END
!
!***********************************************************************
!
      SUBROUTINE SLAMC4( EMIN, START, BASE )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      INTEGER            BASE, EMIN
      REAL               START
!     ..
!
!  Purpose
!  =======
!
!  SLAMC4 is a service routine for SLAMC2.
!
!  Arguments
!  =========
!
!  EMIN    (output) EMIN
!          The minimum exponent before (gradual) underflow, computed by
!          setting A = START and dividing by BASE until the previous A
!          can not be recovered.
!
!  START   (input) REAL
!          The starting point for determining EMIN.
!
!  BASE    (input) INTEGER
!          The base of the machine.
!
! =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I
      REAL               A, B1, B2, C1, C2, D1, D2, ONE, RBASE, ZERO
!     ..
!     .. External Functions ..
      REAL               SLAMC3
      EXTERNAL           SLAMC3
!     ..
!     .. Executable Statements ..
!
      A = START
      ONE = 1
      RBASE = ONE / BASE
      ZERO = 0
      EMIN = 1
      B1 = SLAMC3( A*RBASE, ZERO )
      C1 = A
      C2 = A
      D1 = A
      D2 = A
!+    WHILE( ( C1.EQ.A ).AND.( C2.EQ.A ).AND.
!    $       ( D1.EQ.A ).AND.( D2.EQ.A )      )LOOP
   10 CONTINUE
      IF( ( C1.EQ.A ) .AND. ( C2.EQ.A ) .AND. ( D1.EQ.A ) .AND.         &
     &    ( D2.EQ.A ) ) THEN
         EMIN = EMIN - 1
         A = B1
         B1 = SLAMC3( A / BASE, ZERO )
         C1 = SLAMC3( B1*BASE, ZERO )
         D1 = ZERO
         DO 20 I = 1, BASE
            D1 = D1 + B1
   20    CONTINUE
         B2 = SLAMC3( A*RBASE, ZERO )
         C2 = SLAMC3( B2 / RBASE, ZERO )
         D2 = ZERO
         DO 30 I = 1, BASE
            D2 = D2 + B2
   30    CONTINUE
         GO TO 10
      END IF
!+    END WHILE
!
      RETURN
!
!     End of SLAMC4
!
      END
!
!***********************************************************************
!
      SUBROUTINE SLAMC5( BETA, P, EMIN, IEEE, EMAX, RMAX )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      LOGICAL            IEEE
      INTEGER            BETA, EMAX, EMIN, P
      REAL               RMAX
!     ..
!
!  Purpose
!  =======
!
!  SLAMC5 attempts to compute RMAX, the largest machine floating-point
!  number, without overflow.  It assumes that EMAX + abs(EMIN) sum
!  approximately to a power of 2.  It will fail on machines where this
!  assumption does not hold, for example, the Cyber 205 (EMIN = -28625,
!  EMAX = 28718).  It will also fail if the value supplied for EMIN is
!  too large (i.e. too close to zero), probably with overflow.
!
!  Arguments
!  =========
!
!  BETA    (input) INTEGER
!          The base of floating-point arithmetic.
!
!  P       (input) INTEGER
!          The number of base BETA digits in the mantissa of a
!          floating-point value.
!
!  EMIN    (input) INTEGER
!          The minimum exponent before (gradual) underflow.
!
!  IEEE    (input) LOGICAL
!          A logical flag specifying whether or not the arithmetic
!          system is thought to comply with the IEEE standard.
!
!  EMAX    (output) INTEGER
!          The largest exponent before overflow
!
!  RMAX    (output) REAL
!          The largest machine floating-point number.
!
! =====================================================================
!
!     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
!     ..
!     .. Local Scalars ..
      INTEGER            EXBITS, EXPSUM, I, LEXP, NBITS, TRY, UEXP
      REAL               OLDY, RECBAS, Y, Z
!     ..
!     .. External Functions ..
      REAL               SLAMC3
      EXTERNAL           SLAMC3
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MOD
!     ..
!     .. Executable Statements ..
!
!     First compute LEXP and UEXP, two powers of 2 that bound
!     abs(EMIN). We then assume that EMAX + abs(EMIN) will sum
!     approximately to the bound that is closest to abs(EMIN).
!     (EMAX is the exponent of the required number RMAX).
!
      LEXP = 1
      EXBITS = 1
   10 CONTINUE
      TRY = LEXP*2
      IF( TRY.LE.( -EMIN ) ) THEN
         LEXP = TRY
         EXBITS = EXBITS + 1
         GO TO 10
      END IF
      IF( LEXP.EQ.-EMIN ) THEN
         UEXP = LEXP
      ELSE
         UEXP = TRY
         EXBITS = EXBITS + 1
      END IF
!
!     Now -LEXP is less than or equal to EMIN, and -UEXP is greater
!     than or equal to EMIN. EXBITS is the number of bits needed to
!     store the exponent.
!
      IF( ( UEXP+EMIN ).GT.( -LEXP-EMIN ) ) THEN
         EXPSUM = 2*LEXP
      ELSE
         EXPSUM = 2*UEXP
      END IF
!
!     EXPSUM is the exponent range, approximately equal to
!     EMAX - EMIN + 1 .
!
      EMAX = EXPSUM + EMIN - 1
      NBITS = 1 + EXBITS + P
!
!     NBITS is the total number of bits needed to store a
!     floating-point number.
!
      IF( ( MOD( NBITS, 2 ).EQ.1 ) .AND. ( BETA.EQ.2 ) ) THEN
!
!        Either there are an odd number of bits used to store a
!        floating-point number, which is unlikely, or some bits are
!        not used in the representation of numbers, which is possible,
!        (e.g. Cray machines) or the mantissa has an implicit bit,
!        (e.g. IEEE machines, Dec Vax machines), which is perhaps the
!        most likely. We have to assume the last alternative.
!        If this is true, then we need to reduce EMAX by one because
!        there must be some way of representing zero in an implicit-bit
!        system. On machines like Cray, we are reducing EMAX by one
!        unnecessarily.
!
         EMAX = EMAX - 1
      END IF
!
      IF( IEEE ) THEN
!
!        Assume we are on an IEEE machine which reserves one exponent
!        for infinity and NaN.
!
         EMAX = EMAX - 1
      END IF
!
!     Now create RMAX, the largest machine number, which should
!     be equal to (1.0 - BETA**(-P)) * BETA**EMAX .
!
!     First compute 1.0 - BETA**(-P), being careful that the
!     result is less than 1.0 .
!
      RECBAS = ONE / BETA
      Z = BETA - ONE
      Y = ZERO
      DO 20 I = 1, P
         Z = Z*RECBAS
         IF( Y.LT.ONE )                                                 &
     &      OLDY = Y
         Y = SLAMC3( Y, Z )
   20 CONTINUE
      IF( Y.GE.ONE )                                                    &
     &   Y = OLDY
!
!     Now multiply by BETA**EMAX to get RMAX.
!
      DO 30 I = 1, EMAX
         Y = SLAMC3( Y*BETA, ZERO )
   30 CONTINUE
!
      RMAX = Y
      RETURN
!
!     End of SLAMC5
!
      END
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
c---------------------------------------------------------------------
c     program sigms.f
c
c     calculates debye-waller factors for each multiple
c     scattering path using Debye-Model correlations
c
c     files:  input  pathd_all.dat  multiple scattering path data
c             output fort.3  sig**2 vs path
c                    fort.2  long output
c
c     version 1  (29 july 91)
c
c     coded by j. rehr
c     path data from s. zabinsky
c
c     modified to use pdata.inp, Dec 1991, siz
c     Subroutine version, Dec 1991, siz
c
c---------------------------------------------------------------------

      subroutine sigms (tk, thetad, rs, nlegx, nleg, rat, iz, sig2)
c               tk temperature in degrees K
c               thetad debye temp in degrees K
c               rs=wigner seitz or norman radius in bohr, averaged
c                  over entire problem
c                  (4pi/3)*rs**3 = sum( (4pi/3)rnrm**3 ) / N
c                  (sum is over all atoms in the problem)
c               nlegx used in dimensions of rat and iz
c               nleg nlegs in path
c               rat positions of each atom in path
c               iz atomic number of each atom in path
c               NB Units of distance in this routine
c                  are angstroms, including sig**2.  rs is in bohr.
c               sig2 is output debye waller factor

      implicit double precision (a-h,o-z)

c={../HEADERS/const.h
      parameter (pi = 3.14159 26535 89793 23846 26433d0)
      parameter (one = 1, zero = 0)
      parameter (third = one/3)
      parameter (raddeg = 180 / pi)
      complex*16 coni
      parameter (coni = (0,1))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919 158 292 677 512 811d0)

      parameter (bohr = 0.529 177 249d0, ryd  = 13.605 698d0)
      parameter (hart = 2 * ryd)
      parameter (alpinv = 137.035 989 56d0)
c     fine structure alpha
      parameter (alphfs = 1 / alpinv)
c= ../HEADERS/const.h}

c     nlegx is max number of atoms in any one path
      dimension rat(3,0:nlegx)
      dimension iz(0:nlegx)
c#mn
       external dist

c      parameters
c               x = k_d*R   (distance parameter)
c               R distance in angstroms
c               y = hbar omegad/kT = thetad/t
c               thetad debye temp in degrees K
c               tk temperature in degrees K
c               k_d = (6*pi**2 N/V) = debye wave number
c               N/V=1/(4pi/3rs**3)
c               rs=wigner seitz or norman radius in bohr
c               ami, amj masses at sites i and j in amu
c               I = int_0^1 (y/x) dw sin(wx)coth(wy/2)

c     Note:  There are nleg atoms including the central atom
c            index 0 and index nleg both refer to central atom,
c            which makes special code unnecessary later.

      sigtot=0
      do 800 il=1,nleg
      do 800 jl=il,nleg

c        calculate r_i-r_i-1 and r_j-r_j-1

         rij = dist (rat(1,il), rat(1,jl))
         call corrfn (rij, cij, thetad, tk, iz(il), iz(jl), rs)
         sig2ij=cij

         rimjm = dist (rat(1,il-1), rat(1,jl-1))
         call corrfn (rimjm, cimjm, thetad, tk, iz(il-1), iz(jl-1), rs)
         sig2ij=sig2ij+cimjm

         rijm = dist (rat(1,il), rat(1,jl-1))
         call corrfn (rijm, cijm, thetad, tk, iz(il), iz(jl-1), rs)
         sig2ij=sig2ij-cijm

         rimj = dist (rat(1,il-1), rat(1,jl))
         call corrfn (rimj, cimj, thetad, tk, iz(il-1), iz(jl), rs)
         sig2ij=sig2ij-cimj

         riim = dist (rat(1,il), rat(1,il-1))
         rjjm = dist (rat(1,jl), rat(1,jl-1))

         ridotj=(rat(1,il)-rat(1,il-1))*(rat(1,jl)-rat(1,jl-1))+
     1          (rat(2,il)-rat(2,il-1))*(rat(2,jl)-rat(2,jl-1))+
     2          (rat(3,il)-rat(3,il-1))*(rat(3,jl)-rat(3,jl-1))
         ridotj=ridotj/(riim*rjjm)

c        double count i .ne. j  terms
         if(jl.ne.il) sig2ij=2*sig2ij
         sig2ij=sig2ij*ridotj
         sigtot=sigtot+sig2ij

  800 continue
      sig2=sigtot/4

c     sig2 is in bohr**2, just as we wanted for ff2chi
      return
      end



      subroutine corrfn(rij,cij,thetad,tk,iz1,iz2,rsavg)
c     subroutine calculates correlation function
c     c(ri,rj)=<xi xj> in the Debye approximation
c
c             =(1/N)sum_k exp(ik.(Ri-Rj))(1/sqrt(mi*mj))*
c              (hbar/2w_k)*coth(beta hbar w_k/2)
c             = (3kT/mu w_d**2)*sqrt(mu**2/mi*mj)*I
c
c      parameters
c               x = k_d*R   (distance parameter)
c               R distance in angstroms
c               y = hbar omegad/kT = thetad/t
c               thetad debye temp in degrees K
c               tk temperature in degrees K
c               k_d = (6*pi**2 N/V) = debye wave number
c               N/V=1/(4pi/3rs**3)
c               rs=wigner seitz or norman radius in bohr
c               ami, amj masses at sites i and j in amu
c               I = int_0^1 (y/x) dw sin(wx)coth(wy/2)
c
c      solution by numerical integration
c
      implicit double precision (a-h, o-z)
      common /xy/ x, yinv

c={../HEADERS/const.h
      parameter (pi = 3.14159 26535 89793 23846 26433d0)
      parameter (one = 1, zero = 0)
      parameter (third = one/3)
      parameter (raddeg = 180 / pi)
      complex*16 coni
      parameter (coni = (0,1))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919 158 292 677 512 811d0)

      parameter (bohr = 0.529 177 249d0, ryd  = 13.605 698d0)
      parameter (hart = 2 * ryd)
      parameter (alpinv = 137.035 989 56d0)
c     fine structure alpha
      parameter (alphfs = 1 / alpinv)
c= ../HEADERS/const.h}

c     con=hbar**2/kB*amu)*10**20   in ang**2 units
c     hbar = 1.054 572 666 e-34, amu = 1.660 540 e-27, 
c     kB = 1.380 6581 d-23
      parameter (con = 48.508 459 393 094)
c#mn
       external atwtd

c     external fn
c     rij=2.55
c     tk=295
c     thetad=315
c     ami=amj=63.55 at wt for Cu
c     rs=2.7

      ami=atwtd(iz1)
      amj=atwtd(iz2)
      rs=rsavg
c     thetad in degrees K, t temperature in degrees K
c     y=thetad/tk
      yinv=tk/thetad
      xkd=(9*pi/2)**(third)/(rs*bohr)
      fac=(3/2.)*con/(thetad*sqrt(ami*amj))
      rj=rij
      x=xkd*rj
c     call numerical integration
      call bingrt (grater, eps, nx)
      cij=fac*grater
      return
      end
      double precision function fn(w)
      implicit double precision (a-h,o-z)
      common/xy/x,yinv
c     fn=(sin(wx)/x)*coth(wy/2)
c     change code to allow t=0 without bombing
c     fn=2/y
      fn=2*yinv
      if(w.lt.1.e-20) return
      fac=w
      if(x.gt.0.) fac=sin(w*x)/x
      emwy=0.
      if(yinv.gt.0.0125) emwy=exp(-w/yinv)
      emwy=exp(-w/yinv)
      fn=fac*(1+emwy)/(1-emwy)
      return
      end
c-----------------------------------------------
      subroutine bingrt (b, eps, n)
c     subroutine calculates integrals between [0,1]
c      b = int_0^1 f(z) dz
c     by trapezoidal rule and binary refinement
c     (romberg integration)
c     coded by j rehr (10 Feb 92)
c     see, e.g., numerical recipes for discussion
c     and a much fancier version
c-----------------------------------------------
c     del=dz  itn=2**n tol=1.e-5
c     starting values
      implicit double precision (a-h,o-z)
      common /xy/x,yinv
      character*512 slog
c     external fn
c     error is approximately 2**(-2n) ~ 10**(-.6n)
c     so nmax=10 implies an error of 1.e-6
      parameter(nmax = 10, tol = 1.e-5)
      parameter(zero=0, one=1)
      n=0
      itn=1
      del=1.
      bn=(fn(zero)+fn(one))/2
      bo=bn
 10   continue
c     nth iteration
c     b_n+1=(b_n)/2+deln*sum_0^2**n f([2n-1]deln)
      n=n+1
      if(n.gt.nmax) go to 40
      del=del/2
      sum=0.
      do 20 i=1, itn
      zi=(2*i-1)*del
 20   sum=sum+fn(zi)
c     bnp1=b_n+1 is current value of integral
      bnp1=bn/2+del*sum
c     cancel leading error terms b=[4b-bn]/3
c     note: this is the first term in the
c     neville table - remaining errors were
c     found too small to justify the added code
      b=(4*bnp1-bn)/3
      eps=abs((b-bo)/b)
      if(eps.lt.tol) goto 60
      bn=bnp1
      bo=b
      itn=itn*2
      goto 10
 40   write(slog,50) n,itn, b,eps
      call wlog(slog)
 50   format(' not converged, n,itn,b,eps=',
     1  2i4,2e14.6)
      return
 60   continue
      return
      end
c---------------------------------------------------------------------
c     program sigem
c
c     calculate the Debye-Waller factors for each MS path
c     using the equation-of-motion methods
c
c     input files:  feff.inp and spring.inp
c
c     version 2  ( January 99)
c
c     coded by  A. Poiarkova
c
c---------------------------------------------------------------------
c  References:  
c             for the EM method: Phys. Rev. B , 59, p.948, 1999
c     also see dissertation 
c        "X-ray Absorption Fine Structure Debye-Waller Factors"
c         by Anna V. Poiarkova
c
c---------------------------------------------------------------------
c         tk temperature in degrees K
c         nleg  nlegs in path
c         rat   positions of each atom in path
c         NB Units of distance in this routine
c            are angstroms, including sig2. 
c         sig2 is output DW factor
c
c---------------------------------------------------------------------
      subroutine sigem (sig2mx, sig2x, iem, tk, ipath, nleg, rat, sig2)
      implicit double precision (a-h, o-z)

c={dwpar.h
c-*-fortran-*-
c nlegx1 MUST be the same as legtot, the maximum number of scattering
c       legs in a path
c nphx1 MUST be the same as nphx, the maximum number of atomic species

      parameter (natxdw = 200)
      parameter (nlegx1 = 9)
      parameter (nphx1=7)

c= dwpar.h}
c={../HEADERS/parallel.h
      integer par_type, this_process, numprocs, my_rank
      logical master, worker, parallel_run
      real*8 wall_comm, time_comm
      common /timing/ wall_comm, time_comm
      common /parallel/ numprocs, my_rank, this_process, 
     .          master, worker, parallel_run, par_type
c= ../HEADERS/parallel.h}

c feff parameters (from dim.h):
c     parameter (legtot=9) 
c     parameter (nphx = 7)

      parameter (nphx = nphx1)
      parameter (natx = natxdw)

c local parameters:
      parameter (amu0  = 1.660 54)
      parameter (pi = 3.14159 26535 89793 23846 26433d0)
      parameter (nwx = 700)

      double precision sig2mx, sig2x(0:nphx,0:nphx)
      dimension iphat(natx), izph(0:nphx)

c variables shared with rdspr.f:
      dimension rat1(3,natx), iz(natx)
      dimension dm(3,3,natx,natx)
      dimension rnn(3,natx,natx)
      dimension nnl(natx,natx)

c local variables:
      dimension rat(3,0:nlegx1)
      dimension nconv(0:nlegx1)
      dimension q0(3,natx)
      dimension gr(nwx), w(nwx)
      dimension nq0(0:nlegx1) 
      dimension uu(3,natx), up(3,natx), ff(3,natx)

      character*30  fname
      parameter (ntitx1 = 10)
      character*71  title(ntitx1)
      dimension ltit(ntitx1)
c     character*80  titlep(ntitx1)

      character*512 slog
      logical iem_open

      save 
      data nsigc /0/
c-------------------------------------------------------------

      inquire(unit=iem,opened=iem_open)
      if (nsigc.eq.0) then

c Read coordinates and potentials from feff.inp
      call dwrdin (rat1, iphat, izph, natom,
     1            ntitle, title, ltit)

      if (natom.gt.natx) natom=natx
      do 5 iat=1, natom
         iz(iat) = izph(iphat(iat))
         if (iphat(iat).eq.0) i0=iat
  5   continue

c Read spring.inp and build dynamical matrix
      call rdspr(rat1, iz, natom, i0, 
     1           dm, rnn, 
     1           acut, res, wmax, dosfit, zshell, w0,
     1           rintr, iprdos, nnl)

            write(slog,7)
   7        format(2x,'Calculating Debye-Waller factors via EMM...')
            call wlog(slog)
            write(slog,9)
   9        format(2x,'This might take a while.')
            call wlog(slog)

      if (ipath.ne.0.and.iem_open) then
c           Echo title cards to s2_em.dat
            do 10  i = 1, ntitle
               write(iem,12)  title(i)(1:ltit(i))
  10        continue
  12        format (1x, a)
            write(iem,17) tk, natom
  17        format(1x,'temperature =',f7.2,2x,'N_at =',i4)
            write(iem,19)
  19        format (1x, 71('-'))
            write(iem,25)
            write(slog,25)
            call wlog(slog)
  25        format(3x,'ipath',4x,'nleg',3x,'sig2',5x,
     1            'mu_ipath',2x,'check0(%)')
      endif

c Integration parameters:
      wmaxx=sqrt(zshell)
      dt=2.*pi/wmaxx/15.
c top limit in t integration:
      cutoff=2.*sqrt(2.*acut)/res/wmaxx 
      nstep=cutoff/dt
      xlam=acut/(cutoff)**2
      wl=0.0000001
c top limit in w integration:
      wm=wmax*wmaxx 
      dw=0.01 
      nw=(wm-wl)/dw + 1
      if (nw .gt. nwx) then
          nw = nwx
          dw = (wm-wl)/(nw -1)
      endif
      nfit = dosfit*nw/20.

      endif
c------------------------------------

cc    Open path input file (unit in) and read title.  Use unit 2.
c     ntitle2 = 5
c     open(unit=2,file='paths.dat',status='old', iostat=ios)
c     call chopen (ios, 'paths.dat', 'sigem')
c     call rdhead (2, ntitle2, titlep, ltit)
cc    if (ntitle2 .le. 0)  then
cc       titlep(1) = ' '
cc    endif

c 84  continue
c     read(2,*,end=1010) ipath, nleg
c     skip label (x y z ipot rleg beta eta) and read the path
c     read(2,*)
      do 78 ileg=0,nleg
c        read(2,*,end=1010) (rat(j,ileg),j=1,3)
         nconv(ileg)=0
  78  continue

      do 88 n=1,3
         aa = rat(n,nleg)
         do 87 i=0,(nleg-2)
            j=nleg-i
            rat(n,j)=rat(n,j-1)
  87     continue
         rat(n,1)=aa
  88  continue
      do 89 i=1,nleg
         nq0(i)=0.
  89  continue

c nconv converts # of an atom in the nleg list of coordinates (rat) to
c its # in the full list of all atomic coordinates (rat1)
      do 94 i=1,natom
         do 91 n=1,3
   91    q0(n,i)=0.
         do 95 jl=1,nleg
            m=0
            do 93 n=1,3
               l=nint(100.*rat(n,jl))
               l1=nint(100.*rat1(n,i))
               if (abs(l-l1).le.1) m=m+1
   93       continue
            if (m.eq.3) then
              nconv(jl)=i
            endif
   95    continue
   94 continue
c     check that all path atoms are found
      do 96 jl=1,nleg
        if (nconv(jl).eq.0) then
           print*,' did not find atom jl=', jl
           print*, rat(1,jl),rat(2,jl),rat(3,jl)
           call par_stop('SIGREM-1')
        endif
  96  continue

      atmu=0.
      iq0=0
      nconv(0)=nconv(nleg)
      do 100 il=1,nleg
         l=nconv(il)
         do 101 jq=1,iq0
  101    if(nq0(jq).eq.l) go to 102
         iq0=iq0+1
         nq0(iq0)=l
  102    continue
         nq0x=iq0
         i=nconv(il)
         im=nconv(il-1)
         ip=nconv(il+1)
c        if (il.eq.1) im=nconv(nleg)
         if (il.eq.nleg) ip=nconv(1)
         atmass=atwtd(iz(i))
      do 100 n=1,3
         atmu=atmu + 0.25*( rnn(n,i,im)+rnn(n,i,ip) )**2 /atmass
  100 continue
      atmu=1./atmu
      icount=1
  108 continue
      icount= icount+1
      if (icount.gt.10) call par_stop('SIGREM-2')

      do 115 i=1,natom
      do 115 n=1,3
  115 q0(n,i)=0.

c Build initial state vector |Q_j(0)> for the current path
      do 116 n=1,3
      do 116 il=1,nleg
         i=nconv(il)
         im=nconv(il-1)
         ip=nconv(il+1)
         if (il.eq.1) im=nconv(nleg)
         if (il.eq.nleg) ip=nconv(1)
         atmass=atwtd(iz(i))
         q0(n,i)=q0(n,i)+sqrt(atmu/atmass)*(rnn(n,im,i)-rnn(n,i,ip))/2.
  116 continue

c make sure it's normalized <Q_j(0)|Q_j(0)>=1
      q0q0=0.
      do 120 iq0=1,nq0x
      i=nq0(iq0)
      do 120 n=1,3
         q0q0=q0q0+q0(n,i)*q0(n,i)
 120  continue
      p00=nint(q0q0*1000.d0)/1000.d0
      if (abs(p00-1.d0).gt.5.d-4) then
         atmu=atmu/q0q0
         go to 108
      endif

c     to get THz units:
      wnorm=100.*w0/sqrt(amu0*10.) 
c*** moments
      a0=0.
      do 132 il=1,nq0x
      do 132 im=1,nq0x
         l=nq0(il)
         m=nq0(im)
      do 132 n1=1,3
      do 132 n2=1,3
         a0 = a0 + q0(n1,l)*dm(n1,n2,l,m)*q0(n2,m)/w0/w0
  132 continue
      a0=wnorm*sqrt(a0)

      do 125 kw=1, nwx
         gr(kw) = 0.
  125 w(kw) = (wl+(kw-1)*dw)

c  make file prdennnnn.dat
      if (master.and.ipath.ne.0.and.ipath.le.iprdos) then 
         write(fname,130)  ipath
  130    format('prden', i4.4, '.dat')
         open (unit=25, file=fname, status='unknown',iostat=ios)
         call chopen (ios, fname, 'sigem')
         do 134  i = 1, ntitle
            write(25,136)  title(i)(1:ltit(i))
  134    continue
         write(25,135) natom
  135    format('#',1x,'N_at =', i4)
  136    format ('#',1x, a)
         write(25,138)
  138    format ('#',1x, 71('-'))
      endif


c  set initial conditions
      do 150 i=1, natom
      do 150 n=1,3
         uu(n,i)=q0(n,i)
         up(n,i)=uu(n,i)
  150 continue

c Solve 3*natom equations of motion and find projected VDOS (gr)
      dt2=dt*dt
      t=dt/2.
      do 200 kstep = 1, nstep
c        damping factor:
         e1=exp(-xlam*t*t) 
         xat=0.
         do 167 i=1, natom
         do 167 n=1,3
  167    xat = xat + uu(n,i)*q0(n,i)
         xat=xat*e1
         do 170 kw=1, nw
  170    gr(kw) = gr(kw) + xat*cos(w(kw)*t)*dt
         if(kstep.eq.nstep) go to 200

         do 175 i=1,natom
         do 175 n=1,3
  175    ff(n,i)=0.

         do 180 i=1,natom
            jn=1
  185       if (nnl(i,jn).ne.0) then
               j=nnl(i,jn)
               am=w0*w0
               do 187 n1=1,3
               do 187 n2=1,3
                  ff(n1,i)=ff(n1,i)-dm(n1,n2,i,j)*uu(n2,j)/am
                  if(i.ne.j) ff(n1,j)=ff(n1,j)-dm(n1,n2,j,i)*uu(n2,i)/am
  187          continue
               jn = jn + 1
               go to 185
            endif
  180    continue

         do 199 i=1,natom
         do 199 n=1,3
            put=2.*uu(n,i)-up(n,i)+dt2*ff(n,i)
            up(n,i)=uu(n,i)
            uu(n,i)=put
  199    continue

  200 t=t+dt

      afit = 0.
      if (nfit.ne.0) then
         if (w(nfit).ne.0.) afit=gr(nfit)/(w(nfit)**4)
      endif

c fit vibr.density to A*w^4, for low w
      do 225 kw=1, nfit
         gr(kw)=afit*w(kw)**4
  225 continue

c Normalization of the pr.density of modes 
c (it's the 2/pi factor which was left out till now with,
c perhaps, a small diffrence due to the fit) 
      gr(nw)=0.
      if (gr(1).lt.0.) gr(1)=0.
      xx=(gr(1)+gr(nw))*dw/2.
      do 247 kw=2, (nw-1)
         if (gr(kw).lt.0.) gr(kw)=0.
  247 xx = xx + gr(kw)*dw
      cn1=1./xx

      if (master.and.ipath.ne.0.and.ipath.le.iprdos) then
c to get THz units:
         wnorm=100.*w0/sqrt(amu0*10.) 
         write(25,349) ipath, nleg
  349    format('#',2x,'ipath =',i3,2x,'nleg =',i2)
         write(25,350)
c 350    format(1h#,6x,5hw,THz,18x,6hrho(w))
  350    format(1h#,6x,'cm^-1',18x,6hrho(w))
         do 370 kw=1,nw
            write(25,360) w(kw)*wnorm*100./6./pi, gr(kw)*cn1/wnorm
c           write(25,360) w(kw)*wnorm, gr(kw)*cn1/wnorm
  360       format(2x,f10.3,15x,f10.7)
  370    continue
         close (unit=25)
      endif

      wt=tk/187.64/w0
      ccc=cn1
      check0 = abs((2./pi - cn1)/(2./pi))
      check0=check0*100.
      coef = ccc*0.5*0.2587926/atmu/w0
c integrate over w to get sig2
      cth=0.
      s2=0.
c     gr(1)=0.
      do 400 kw=2, (nw-1)
         cth = 1./tanh( w(kw)/(2.*wt) )
         s2 = s2 + coef*gr(kw)*cth*dw/w(kw)
  400 continue
      sig2 = s2

      if (ipath.ne.0.and.iem_open) then
         write(iem,473) ipath, nleg, sig2,atmu,check0
         write(slog,473) ipath, nleg, sig2,atmu,check0
         call wlog(slog)
  473    format(4x,i3,4x,i3,4x,f7.5,3x,f7.3,4x,f5.2)
      endif

      nsigc = nsigc + 1
c1000 go to 84
c1010 continue
c     close (unit=2)
      if (sig2.gt.1.0) then
        sig2 = 1.0d0
        call wlog (' WARNING: Found sig**2>1. Set sig2=1. ')
        write (slog,1011) nconv(1), nconv(2)
 1011   format('          Possible zero ferquency modes with atoms', i4,
     1   ' or', i4)
         call wlog(slog)
         call wlog('          Check springs.inp')
      endif
      if (check0.gt.5.0) then
        write (slog,*) ' WARNING: Failed check0 test:missing VDOS.',
     1  ' Reduce dosfit and/or increase wmax.'
        call wlog(slog)
      endif

c     update maximum DW factors
      if (sig2.gt.sig2mx) sig2mx=sig2
      if (sig2.gt.sig2x( iphat(nconv(1)),  iphat(nconv(2)) )) then
         sig2x( iphat(nconv(1)),  iphat(nconv(2)) ) = sig2
         sig2x( iphat(nconv(2)),  iphat(nconv(1)) ) = sig2
      endif

      return
      end
c----------------------------------------------------
      subroutine dwrdin (rat, iphat, izph, nat,
     1            ntitle, title, ltit)

c     Read feff.inp for sigem.f
c     (here we need only coordinates and potentials)

      implicit double precision (a-h, o-z)

c={dwpar.h
c-*-fortran-*-
c nlegx1 MUST be the same as legtot, the maximum number of scattering
c       legs in a path
c nphx1 MUST be the same as nphx, the maximum number of atomic species

      parameter (natxdw = 200)
      parameter (nlegx1 = 9)
      parameter (nphx1=7)

c= dwpar.h}

c feff parameters:
c     parameter (nphx = 7)

      parameter (nphx = nphx1)
      parameter (natx = natxdw)

      dimension iphat(natx)
      dimension rat(3,natx)
      dimension iatph(0:nphx)

      character*6  potlbl(0:nphx)

c     Local stuff
      character*150  line
      parameter (nwordx = 20)
      character*20 words(nwordx)

      parameter (ntitx = 10)
      character*71  title(ntitx)
      dimension ltit(ntitx)
      dimension izph(0:nphx)
      logical iscomm
      parameter (nssx = 16)

      parameter (big = 1.0e5)
      character*512 slog

   10 format (a)
   20 format (bn, i15)
   30 format (bn, f15.0)

cc    initialize things

      ntitle = 0

      nat = 0
      do 100  iat = 1, natx
         iphat(iat) = -1
  100 continue

      nph = 0
      do 110  iph = 0, nphx
         iatph(iph) = 0
         izph(iph) = 0
         potlbl(iph) = ' '
  110 continue

c     Open feff.inp, the input file we're going to read
      open (unit=1, file='feff.inp', status='old', iostat=ios)
      call chopen (ios, 'feff.inp', 'rdinp')

c     tokens  0 if not a token
c             1 if ATOM (ATOMS)
c             7 if TITL (TITLE)
c            10 if DEBY (DEBYE)
c            13 if PRIN (PRINT)
c            14 if POTE (POTENTIALS)
c            -1 if END  (end)
c     mode flag  0 ready to read a keyword card
c                1 reading atom positions
c                2 reading overlap instructions for unique pot
c                3 reading unique potential definitions

      mode = 0
  200 read(1,10,iostat=ios)  line
         if (ios .lt. 0)  line='END'
         call triml (line)
         if (iscomm(line))  goto 200
         nwords = nwordx
         call bwords (line, nwords, words)
         itok = itoken (words(1),'feff.inp')

c        process the card using current mode
  210    continue

         if (mode .eq. 0)  then
            if (itok .eq. 1)  then
c              ATOM
c              Following lines are atom postions, one per line
               mode = 1

            elseif (itok .eq. 7)  then
c              TITLE title...
               ntitle = ntitle + 1
               if (ntitle .le. ntitx)  then
                  title(ntitle) = line(6:)
                  call triml (title(ntitle))
               else
                  call wlog(' Too many title lines, title ignored')
                  call wlog(' ' // line(1:71))
               endif
               mode = 0

c           elseif (itok .eq. 10)  then
cc             DEBYE  temp debye-temp
cc                  temps in kelvin
cc                  These add to any sig2 from SIG2 card or files.dat
c              read(words(2),30,err=900)  tk
c              read(words(3),30,err=900)  thetad
c              idwopt=0
c              read(words(4),20,err=900)  idwopt
c              mode = 0
            elseif (itok .eq. 14)  then
c              POTENTIALS
c              Following lines are unique potential defs, 1 per line
               mode = 3

            elseif (itok .eq. -1)  then
cc             END
               goto 220
            else
               mode = 0
c *            write(slog,'(1x,a)') line(1:70)
c *            call wlog(slog)
c *            write(slog,'(1x,a)') words(1)
c *            call wlog(slog)
c *            write(slog,'(a,i8)') ' Token ', itok
c *            call wlog(slog)
c *            call wlog(' Keyword unrecognized.')
c *            call wlog(' See FEFF document -- some old features')
c *            call wlog(' are no longer available.')
c *            call par_stop('DWRDIN-1')
            endif
         elseif (mode .eq. 1)  then
            if (itok .ne. 0)  then
cc             We're done reading atoms.
cc             Change mode and process current card.
               mode = 0
               goto 210
            endif
            nat = nat+1
            if (nat .gt. natx)  then
               write(slog,'(1x,a,i5)') 'Too many atoms, max is ', natx
               call wlog(slog)
               write(slog,'(1x,a,i5,a)') 'Only', natx,
     1       ' atoms will be considered in the DW factor calculations.'
               call wlog(slog)
               nat = nat-1
               mode = 0
               goto 210
c              call par_stop('DWRDIN-2')
            endif
            if (nat.le.natx) then
               read(words(1),30,err=900)  rat(1,nat)
               read(words(2),30,err=900)  rat(2,nat)
               read(words(3),30,err=900)  rat(3,nat)
               read(words(4),20,err=900)  iphat(nat)
               if (iphat(nat).eq.0) iat0 = nat
            else
               mode = 0
               goto 210
            endif
         elseif (mode .eq. 3)  then
            if (itok .ne. 0)  then
cc             We're done reading unique potential definitions
cc             Change mode and process current card.
               mode = 0
               goto 210
            endif
            read(words(1),20,err=900)  iph
            if (iph .lt. 0  .or.  iph .gt. nphx)  then
               write(slog,'(a,i8)') 
     1             'Unique potentials must be between 0 and ',
     1             nphx
               call wlog(slog)
               write(slog,'(i8,a)') iph, ' not allowed'
               call wlog(slog)
               write(slog,'(1x,a)') line(1:71)
               call wlog(slog)
               call par_stop('DWRDIN-3')
            endif
            read(words(2),20,err=900)  izph(iph)
cc          No potential label if user didn't give us one
cc          Default set above is potlbl=' '
            if (nwords .ge. 3)  potlbl(iph) = words(3)
         else
            write(slog,'(a,i8)') 
     .        'DWRDIN-4: Mode unrecognized, mode ', mode
c           call wlog(slog)
            call par_stop(slog)
         endif
      goto 200
  220 continue

cc    We're done reading the input file, close it.
      close (unit=1)
            if (nat .gt. natx)  then
               write(slog,'(a,i8)') 'Too many atoms for DW calculations,
     1         max is ', natx
               call wlog(slog)
               write(slog,'(a,i8,a)') 'Only atoms up to #',natx,
     1         '  will be considered'
               call wlog(slog)
            endif

      do 250 iat = 1, nat
      do 250 i = 1,3
        if (iat.ne. iat0) rat(i,iat) = rat(i,iat) - rat(i,iat0)
 250  continue
      do 251 i = 1,3
 251  rat(i,iat0) = 0.d0

cc    Find out how many unique potentials we have
      nph = 0
      do 300  iph = nphx, 0, -1
         if (izph(iph) .gt. 0)  then
            nph = iph
            goto 301
         endif
  300 continue
  301 continue
cc    Must have central atom
      if (izph(0) .le. 0)  then
         call wlog(' No absorbing atom (unique pot 0) was defined.')
         call par_stop('DWRDIN-5')
      endif
cc    Find central atom (only 1 permitted)
      iatabs = -1
      do 400  iat = 1, nat
         if (iphat(iat) .eq. 0)  then
            if (iatabs .lt. 0)  then
               iatabs = iat
            else
               call wlog(' More than one absorbing atom (potential 0)')
               call wlog(' Only one absorbing atom allowed')
               call par_stop('DWRDIN-6')
            endif
         endif
  400 continue

cc    Then find model atoms for unique pots that have them
cc    Use atom closest to absorber for model
      do 330  iph = 0, nphx
         rabs = big
         do 320  iat = 1, nat
            if (iph .eq. iphat(iat))  then
               tmp = dist (rat(1,iat), rat(1,iatabs))
               if (tmp .lt. rabs)  then
cc                this is the closest so far
                  rabs = tmp
                  iatph(iph) = iat
               endif
            endif
  320    continue
  330 continue
cc    if iatph > 0, a model atom has been found.

      if (ntitle .le. 0)  then
         ntitle = 1
         title(1) = 'Null title'
      endif
      do 490  i = 1, ntitle
         ltit(i) = istrln (title(i))
  490 continue

      return

  900 continue
      call wlog(' Error reading input, bad line follows:')
      write(slog,'(1x,a)') line(1:71)
      call wlog(slog)
      call par_stop('DWRDIN-7 fatal error.')

c      return
      end

c----------------------------------------------------------
      subroutine rdspr(rat1, iz, natom, i0, 
     1           dm, rnn, 
     1           acut, res, wmax, dosfit, zshell, w0, 
     1           rintr, iprdos, nnl)

c     Read spring.inp for multiple scattering feff and
c     build dynamical matrix.

      implicit double precision (a-h, o-z)

c={dwpar.h
c-*-fortran-*-
c nlegx1 MUST be the same as legtot, the maximum number of scattering
c       legs in a path
c nphx1 MUST be the same as nphx, the maximum number of atomic species

      parameter (natxdw = 200)
      parameter (nlegx1 = 9)
      parameter (nphx1=7)

c= dwpar.h}
c={../HEADERS/parallel.h
      integer par_type, this_process, numprocs, my_rank
      logical master, worker, parallel_run
      real*8 wall_comm, time_comm
      common /timing/ wall_comm, time_comm
      common /parallel/ numprocs, my_rank, this_process, 
     .          master, worker, parallel_run, par_type
c= ../HEADERS/parallel.h}

c feff parameters:

c     parameter (nphx = nphx1)
      parameter (natx = natxdw)

c new local parameters:
      parameter (nangx = 7*natx)
      parameter (nsprx = 40)
      parameter (nshx = 100)

c variables shared with sigem.f:
      dimension rat1(3,natx), iz(natx)
      dimension dm(3,3,natx,natx)
      dimension rnn(3,natx,natx)
      dimension nnl(natx,natx)

c local variables:
      dimension rshell(natx,0:nshx)
      dimension nspr(2,nsprx), drij(natx,natx)
      dimension str(natx,natx)
      dimension ang(nangx), dang(nangx)
      dimension nang(3,nangx)
      dimension dmstr(3,3,natx,natx),dma(3,3,natx,natx)
      dimension si(3), sj(3), sk(3)

      character*150  line
      parameter (nwordx = 20)
      character*20  words(nwordx)
      character*512 slog

      logical iscomm

   10 format (a)
   20 format (bn, i15)
   30 format (bn, f15.0)

c initialize things

      do 40 i=1, natom
      do 40 j=1, natom
         str(i,j)=0.
         drij(i,j)=0.02
   40 continue
      do 47 na=1,nangx
         ang(na)=0.
         dang(na)=0.
      do 47 m=1,3
         nang(m,na)=0
   47 continue

      do 50 ispr=1, nsprx
      do 50 n=1,2
   50 nspr(n,ispr)=0

      acut=3.
      res=0.05
      dosfit=0.
      wmax=1.
      na=1
      nintr=0
      strx=10000.
      ispr=1
      iprdos = 0
      ddrij=0.02
      ddang=0.02

      open(unit=1,file='spring.inp',status='old', iostat=ios)
      call chopen (ios, 'spring.inp', 'rdspr')

c     tokens  0 if not a token
c             1 if STRE (STRETCHES)
c             2 if ANGL (ANGLES)
c             3 if VDOS
c             4 if PRDOS 
c            -1 if END  (end)
c     mode flag  0 ready to read a keyword card
c                1 reading stretches
c                2 reading angle-bends
 
      mode = 0
  200 read(1,10,iostat=ios)  line
         if (ios .lt. 0)  line='END'
         call triml (line)
         if (iscomm(line))  goto 200
         nwords = nwordx
         call bwords (line, nwords, words)
         itok = itoken (words(1),'spring.inp')

c        process the card using current mode
  210    continue

         if (mode .eq. 0)  then
            if (itok .eq. 1)  then
c              STRE
c              Following lines are stretches, one per line
c              read(words(2),20,err=900)  nintr
               mode = 1
            elseif (itok .eq. 2)  then
c              ANGL
c              Following are ...
               mode = 2
            elseif (itok .eq. 3)  then
c              VDOS
c              VDOS  resolution, a_cut, wmax, dosfit
c               0 - do not run modules, 1 - run module
               read(words(2),30,err=900)  res
               read(words(3),30,err=900)  wmax
               read(words(4),30,err=900)  dosfit
               if (nwords.gt.4) then
                   read(words(5),30,err=900)  acut
               endif
               mode = 0
            elseif (itok .eq. 4)  then
c              PRINT  iprdos
c              to print or not to print prdennnnn.dat files;
c              if the card is present, these files will be
c              printed for paths 1 through iprdos
               iprdos = 1 
               read(words(2),20,err=900)  iprdos
               mode = 0
            elseif (itok .eq. -1)  then
c              END
               goto 220
            else
               write(slog,'(1x,a)') line(1:71)
               call wlog(slog)
               write(slog,'(1x,a)') words(1)
               call wlog(slog)
               write(slog,'(a,i8)') ' Token ', itok
             call wlog(slog)
               call wlog(' Keyword unrecognized.')
               call wlog(' See FEFF document -- some old features')
               call wlog(' are no longer available.')
             call par_stop('RDSPR-1')
            endif
         elseif (mode .eq. 1)  then
            if (itok .ne. 0)  then
c              We're done reading stretches
c              Change mode and process current card.
               mode = 0
               goto 210
            endif
            read(words(1),20,err=900) ii
            i=ii+1
            call chekin (ii, natom, line)
            read(words(2),20,err=900) jj
            j=jj+1
            call chekin (jj, natom, line)
            read(words(3),30,err=900) str(i,j)
            if (str(i,j).lt.strx) then
               strx=str(i,j)
               ix=i
               jx=j
            endif
            nspr(1,ispr)=i
            nspr(2,ispr)=j
            ispr=ispr+1
            read(words(4),30,err=900) ddrij
            drij(i,j) = abs(ddrij)/100.
            drij(j,i) = abs(ddrij)/100.
         elseif (mode .eq. 2)  then
            if (itok .ne. 0)  then
c              We're done reading angle-bends
c              Change mode and process current card.
               mode = 0
               goto 210
            endif
            read(words(1),20,err=900) ii
            i=ii+1
            call chekin (ii, natom, line)
            read(words(2),20,err=900) jj
            j=jj+1
            call chekin (jj, natom, line)
            read(words(3),20,err=900) kk
            k=kk+1
            call chekin (kk, natom, line)
            read(words(4),30,err=900) ang(na)
            nang(1,na)=i
            nang(2,na)=j
            nang(3,na)=k
            read(words(5),30,err=900) ddang
            dang(na) = abs(ddang)/100.
            na=na+1
         else
            write(slog,'(a,i8)') 'Mode unrecognized, mode ', mode
            call wlog(slog)
            call par_stop('RDSPR-2')
         endif
      goto 200
  220 continue

c     We're done reading the input file, close it.
      close (unit=1)
      nax=na-1

c     write statistics on found bonds and angles into spring.dat 
      if (master) then
        open (unit=2, file='spring.dat', status='unknown',iostat=ios)
        call chopen (ios, 'spring.dat', 'spring')
        write(2,*) ' Statistics on spring constants in spring.inp.'
        write(2,*) '   STRETCHES  i  j  aa   found_number'  
      endif

c find all stretching bonds
      do 321 jspr=1, (ispr-1)
         icnt=0
         i=nspr(1,jspr)
         j=nspr(2,jspr)
         aa = str(i,j)
         ddrij=drij(i,j)
         ip=iz(i)
         jp=iz(j)
         rij = dist (rat1(1,i), rat1(1,j))
         if (aa.eq.0.) go to 321
         do 320 k=1, natom
            do 320 l=k+1, natom
               kp=iz(k)
               lp=iz(l)
               rkl = dist (rat1(1,k), rat1(1,l))
               comp = abs(rij/rkl - 1.)
               if (comp.gt.ddrij) go to 320
               if (ip.ne.kp.or.jp.ne.lp) then
                  if (ip.ne.lp.or.jp.ne.kp) go to 320
               endif
               str(k,l) = aa
               str(l,k) = aa
calex       to check the bonds, that were found
calex        print*, k,l,aa
                icnt = icnt+1
  320    continue
         str(j,i) = aa
         if (master) write (2,*) i-1, j-1, aa, icnt
  321 continue
      if (master) write(2,*) '   BENDS   i  j  k   aa   found_number'  

c find all bending angles
      naxx=nax
      do 323 na=1,nax
         icnt=1
         i=nang(1,na)
         j=nang(2,na)
         k=nang(3,na)
         ddrij=drij(i,j)
         ddrkj=drij(k,j)
         ip=iz(i)
         jp=iz(j)
         kp=iz(k)
         call coss(rat1(1,i),rat1(1,j),rat1(1,k),cosijk)
         rij = dist (rat1(1,i), rat1(1,j))
         rkj = dist (rat1(1,k), rat1(1,j))
         aa=ang(na)
c        print*, na, i,j,k, aa
         do 326 ii=1, natom
         do 326 jj=1, natom
            if (ii.eq.jj) go to 326
            rrij=dist (rat1(1,ii), rat1(1,jj))
            do 322 kk=ii+1, natom
               if (kk.eq.jj) go to 322
               rrkj=dist (rat1(1,kk), rat1(1,jj))
               comp1 = abs(rrij/rij - 1.)
               comp2 = abs(rrkj/rkj - 1.)
               if (comp1.gt.ddrij.or.comp2.gt.ddrkj) then
                  comp1 = abs(rrkj/rij - 1.)
                  comp2 = abs(rrij/rkj - 1.)
                  if (comp1.gt.ddrij.or.comp2.gt.ddrkj) go to 322
               endif
               iip=iz(ii)
               jjp=iz(jj)
               kkp=iz(kk)
            if (iip.ne.ip.or.jjp.ne.jp.or.kkp.ne.kp) then
               if (kkp.ne.ip.or.jjp.ne.jp.or.iip.ne.kp) go to 322
            endif
               call coss(rat1(1,ii),rat1(1,jj),rat1(1,kk),cssijk)
               if (dacos(cosijk).eq.0.) go to 322
                  comp = abs( dacos(cssijk)/dacos(cosijk) -1.)
               if (comp.ge.dang(na)) go to 322
               do 324 na1=1,naxx
                  ii1=nang(1,na1)
                  jj1=nang(2,na1)
                  kk1=nang(3,na1)
               if (ii.eq.ii1.and.jj.eq.jj1.and.kk.eq.kk1) go to 322
               if (kk.eq.ii1.and.jj.eq.jj1.and.ii.eq.kk1) go to 322
 324           continue
               naxx=naxx+1
               ang(naxx)=aa
               nang(1,naxx)=ii
               nang(2,naxx)=jj
               nang(3,naxx)=kk
calex          to check the bends, that were found
c              print*, naxx, ii,jj,kk,aa
               icnt = icnt + 1
               if (naxx.eq.nangx) goto 333
 322        continue
 326     continue
         if (master) write (2,*) i-1, j-1, k-1, aa, icnt
 323  continue
 333  continue

      if (master) close (unit=2)

      do 325 i=1, natom
      do 325 nshell=0, nshx
         rshell(i,nshell)=0.
 325  continue

c find shells
      rintr=0.
      nintr=1
      do 330 i=1, natom
      nshell=0
      do 335 j=1, natom
         if (j.eq.i) go to 332
         if (nshell.gt.nshx) go to 332
         rij = dist (rat1(1,i), rat1(1,j))
         ddrij=drij(i,j)
         ncount=0
         do 331 ish=0, nshell
            b = real(rshell(i,ish))
            dif=1.
            if (b.ne.0.) dif = abs(rij -b)/b
            if (dif.le.ddrij) ncount=ncount+1 
 331     continue
         if (ncount.eq.0) then
            nshell = nshell + 1
            if (str(i,j).ne.0.and.rij.gt.rintr) rintr=rij
            rshell(i,nshell) = rij
         endif
 332     do 335 n=1,3
         rnn(n,i,j)=0.
         do 335 m=1,3
            dmstr(n,m,i,j)=0.
            dma(n,m,i,j)=0.
            dm(n,m,i,j)=0.
 335  continue
c sort rshell into ascending numerical order
c and find maximum order of interacting neighbor nintr
      do 342 jsh=2,nshell 
         aa = rshell(i,jsh)
         do 341 ish=jsh-1,1,-1
            if(rshell(i,ish).le.aa) go to 340
            rshell(i,ish+1)=rshell(i,ish)
 341     continue
         ish=0
 340     rshell(i,ish+1) = aa
         if (aa.le.rintr.and.(ish+1).gt.nintr) nintr = ish+1
 342  continue
 330  continue

      zshell=0.
      i1=0
      do 350 i=1,natom
         do 352 in=1,natx
 352     nnl(i,in)=0
         do 350 j=i+1,natom
            dx=rat1(1,j)-rat1(1,i)
            dy=rat1(2,j)-rat1(2,i)
            dz=rat1(3,j)-rat1(3,i)
            dr=sqrt(dx*dx+dy*dy+dz*dz)
            rnn(1,i,j)=dx/dr
            rnn(2,i,j)=dy/dr
            rnn(3,i,j)=dz/dr
            rnn(1,j,i)=-rnn(1,i,j)
            rnn(2,j,i)=-rnn(2,i,j)
            rnn(3,j,i)=-rnn(3,i,j)
            rrij = abs( dr/rshell(1,1) -1.)
c           if (i.eq.1.and.rrij.le.drij(i,j)) zshell=zshell+1
            if (i.eq.i0.and.rrij.le.drij(i,j)) then
               zshell=zshell+1
               if (i1.eq.0.and.str(i,j).ne.0.) i1 = j
            endif
  350 continue

c Build dynm. matrix for angle bends 
      nan=0
      do 355 na=1,naxx
c        print*,na, naxx, nangx
         i=nang(1,na)
         j=nang(2,na)
         k=nang(3,na)
         if (i.eq.j.or.j.eq.k) go to 355
         if(ang(na).eq.0.) go to 355
         nan=nan+1
         rij = dist(rat1(1,i),rat1(1,j))
         rkj = dist(rat1(1,k),rat1(1,j))
         if (rij.gt.rintr.or.rkj.gt.rintr) go to 355
         call sang (i, j, k, rat1, si, sj, sk)
         do 357 n1=1,3
         do 357 n2=1,3
            dma(n1,n2,i,j)=dma(n1,n2,i,j)+ang(na)*si(n1)*sj(n2)
            dma(n1,n2,j,k)=dma(n1,n2,j,k)+ang(na)*sj(n1)*sk(n2)
            dma(n1,n2,i,k)=dma(n1,n2,i,k)+ang(na)*si(n1)*sk(n2)

            dma(n1,n2,i,i)=dma(n1,n2,i,i)+ang(na)*si(n1)*si(n2)
            dma(n1,n2,j,j)=dma(n1,n2,j,j)+ang(na)*sj(n1)*sj(n2)
            dma(n1,n2,k,k)=dma(n1,n2,k,k)+ang(na)*sk(n1)*sk(n2)

            dma(n2,n1,j,i)=dma(n1,n2,i,j)
            dma(n2,n1,k,j)=dma(n1,n2,j,k)
            dma(n2,n1,k,i)=dma(n1,n2,i,k)
  357    continue
  355 continue

c Build dynm. matrix for stretches 
      do 375 l=1,natom
      do 375 m=l,natom
         do 373 n1=1,3
            x2=str(l,m)*rnn(n1,l,m)
         do 373 n2=1,3
         dmi=0.
         if (l.eq.m) then
            do 377 i=1,natom
               if (real(str(i,m)).eq.0.) go to 377
               dmi = dmi + str(i,m)*rnn(n1,i,m)*rnn(n2,i,m)
 377        continue
         endif
         dmstr(n1,n2,l,m) = dmi - x2*rnn(n2,l,m)
         dmstr(n2,n1,m,l)=dmstr(n1,n2,l,m)
 373  continue
 375  continue

c Add two dynm. matrices D_str+D_ang
      lnx=0
      do 380 i=1,natom
      ami=sqrt(atwtd(iz(i)))
      in=0
      do 380 j=1,natom
      amj=sqrt(atwtd(iz(j)))
      sumdm=0.
      do 381 n1=1,3
      do 381 n2=1,3
         dmdm=dma(n1,n2,i,j)+dmstr(n1,n2,i,j)
         sumdm=sumdm+abs(dmdm)
         dm(n1,n2,i,j)=dmdm/ami/amj
 381  continue
      if (real(sumdm).ne.0.and.i.le.j) then
         in=in+1
         nnl(i,in)=j
      endif
      if (in.ge.lnx) lnx=in
 380  continue

      atmu = 1./(1./atwtd(iz(i0)) + 1./atwtd(iz(i1)))
      a0=0.
      do 450 i=1,2
      do 450 j=1,2
      if (i.eq.1) l=i0
      if (i.eq.2) l=i1
      if (j.eq.1) m=i0
      if (j.eq.2) m=i1
      do 450 n1=1,3
      do 450 n2=1,3
         fact = (-1)**(i+j)*atmu
         atmass = 1./atwtd(iz(l))/atwtd(iz(m))
         a0 = a0 + fact*sqrt(atmass)*rnn(n1,i0,i1)*
     1             dm(n1,n2,l,m)*rnn(n2,i0,i1)
  450 continue
c     effective freq. for the 1st shell:
      w0=sqrt(a0) 
      if (w0.eq.0.) then
         atmux = 1./(1./atwtd(iz(ix)) + 1./atwtd(iz(jx)))
         w0=sqrt(strx/atmux)
      endif

      return

  900 continue
      call wlog(' Error reading input, bad line follows:')
      write(slog,'(1x,a)') line(1:71)
      call wlog(slog)
      call par_stop('RDSPR fatal error.')

cc      return
      end

c---------------------------------------------------

      subroutine sang (i, j, k, rat1, si, sj, sk)

c*calculates coeficients s  (Sint=sum {si*ui}) connecting internal
c*coordinate delta_phi(ijk)=(valence ijk angle bend) with ui (atomic
c*displacements)

      implicit double precision (a-h, o-z)

c={dwpar.h
c-*-fortran-*-
c nlegx1 MUST be the same as legtot, the maximum number of scattering
c       legs in a path
c nphx1 MUST be the same as nphx, the maximum number of atomic species

      parameter (natxdw = 200)
      parameter (nlegx1 = 9)
      parameter (nphx1=7)

c= dwpar.h}

      parameter (natx = natxdw)

      dimension rat1(3,natx), rji(3), rjk(3)
      dimension eji(3), ejk(3), ej(3)
      dimension si(3), sk(3), sj(3)

      dji=0.
      djk=0.
      dik=0.
      do 905 m = 1, 3
         rji(m) = rat1(m,i) - rat1(m,j)
         rjk(m) = rat1(m,k) - rat1(m,j)
         dji = dji + rji(m)**2
         djk = djk + rjk(m)**2
         dik = dik + ( rat1(m,k) - rat1(m,i) )**2
         si(m) = 0.
         sj(m) = 0.
         sk(m) = 0.
  905 continue
      dji = sqrt(dji)
      djk = sqrt(djk)
      dik = sqrt(dik)

      dotj=0.
      do 910 m = 1, 3
         eji(m) = rji(m)/dji
         ejk(m) = rjk(m)/djk
         dotj = dotj + eji(m) * ejk(m)
  910 continue
c     ri = dji
c     rk = djk
c     rj = sqrt(dji*djk)
c     rj = dik
      rj=1.
      call vect (eji, ejk, ej, sinj)
      do 920 m = 1, 3
         si(m) = rj*(dotj * eji(m) - ejk(m))/dji/sinj
         sk(m) = rj*(dotj * ejk(m) - eji(m))/djk/sinj
         sj(m) = rj*((dji - djk * dotj)*eji(m) +
     1              (djk - dji * dotj)*ejk(m))/dji/djk/sinj
  920 continue
      return
      end

c-----------------------------------------------------------
      subroutine vect (v1, v2, v3, sin12)

c*calculates vector product v3 = [v1 x v2] and sin of the angle

      implicit double precision (a-h, o-z)

      dimension v1(3), v2(3), v3(3)

      v3(1) = v1(2)*v2(3) - v1(3)*v2(2)
      v3(2) = v1(3)*v2(1) - v1(1)*v2(3)
      v3(3) = v1(1)*v2(2) - v1(2)*v2(1)
      d1 = 0.
      d2 = 0.
      d3 = 0.
      do 990 m = 1, 3
         d1 = d1 + v1(m)**2
         d2 = d2 + v2(m)**2
         d3 = d3 + v3(m)**2
  990 continue
      sin12 = sqrt(d3/d1/d2)
      return
      end

c-----------------------------------------------------------
      subroutine coss (v1,v2,v3,cos12)
c* calculates cos between two vectors v1-v2 and v3-v2
      implicit double precision (a-h, o-z)
      dimension v1(3), v2(3), v3(3)

      vv1=0.
      vv2=0.
      scal=0.
      do 995 m=1,3
         vv1=vv1+(v1(m)-v2(m))**2
         vv2=vv1+(v3(m)-v2(m))**2
         scal=scal+(v1(m)-v2(m))*(v3(m)-v2(m))
  995 continue
      cos12=scal/vv1/vv2
      return
      end

c-----------------------------------------------------------
      subroutine chekin (i, natom, line)
      character*150  line
      character*512  slog

            if (i .gt. (natom-1) .or. i .lt. 0) then
               write(slog,'(a,i8)')
     1             'the atomic indexes must be between 0 and ',
     1             (natom - 1)
               call wlog(slog)
               write(slog,'(i8,a)') i, ' not allowed'
               call wlog(slog)
               write(slog,'(1x,a)') line(1:71)
               call wlog(slog)
               call par_stop('RDSPR')
            endif
        return
        end

c---------------------------------------------------------------------
c     program sigrm
c
c     calculate the Debye-Waller factors for each MS path
c     using the recursion method
c
c     input files:  feff.inp and spring.inp
c
c     version 2  ( January 99)
c
c     coded by  A. Poiarkova
c
c---------------------------------------------------------------------
c  References:
c             for the RM: J. Synchrotron Rad., 1999 (to bu published)
c     also see dissertation
c        "X-ray Absorption Fine Structure Debye-Waller Factors"
c         by Anna V. Poiarkova
c
c---------------------------------------------------------------------
c         tk temperature in degrees K
c         nleg  nlegs in path
c         rat   positions of each atom in path
c         NB Units of distance in this routine
c            are angstroms, including sig2.
c         sig2 is output DW factor
c
c---------------------------------------------------------------------
      subroutine sigrm (sig2mx, sig2x,ir1, ir2, tk,ipath,nleg,rat,sig2)
      implicit double precision (a-h, o-z)

c={dwpar.h
c-*-fortran-*-
c nlegx1 MUST be the same as legtot, the maximum number of scattering
c       legs in a path
c nphx1 MUST be the same as nphx, the maximum number of atomic species

      parameter (natxdw = 200)
      parameter (nlegx1 = 9)
      parameter (nphx1=7)

c= dwpar.h}

c feff parameters (from dim.h):
c     parameter (legtot=9) 
c     parameter (nphx = 7)

      parameter (nphx = nphx1)
      parameter (natx = natxdw)

c local parameters:
      parameter (amu0  = 1.660 54)
      double precision sig2mx, sig2x(0:nphx,0:nphx)
      dimension iphat(natx), izph(0:nphx)

c variables shared with rdspr.f:
      dimension rat1(3,natx), iz(natx)
      dimension dm(3,3,natx,natx)
      dimension rnn(3,natx,natx)
      dimension nnl(natx,natx)

c local variables:
      dimension rat(3,0:nlegx1)
      dimension nconv(0:nlegx1)
      dimension q0(3,natx)
c   list of atoms |0>=|Q>:
      dimension nq0(0:nlegx1)  
c   state |1>=D|Q>:
      dimension q1(3,natx)  
c   list of atoms in |1>:
      dimension nq1(natx)   

c     character*30  fname
      parameter (ntitx1 = 10)
      character*71  title(ntitx1)
      dimension ltit(ntitx1)
c     character*80  titlep(ntitx1)

      character*512 slog

      logical ir1_open, ir2_open

      save 
      data nsigc /0/
c-------------------------------------------------------------

      inquire(unit=ir1,opened=ir1_open)
      inquire(unit=ir2,opened=ir2_open)

      if (nsigc.eq.0) then
c        Read coordinates and potentials from feff.inp
         call dwrdin (rat1, iphat, izph, natom,
     1            ntitle, title, ltit)

         if (natom.gt.natx) natom=natx
         do 5 iat=1, natom
            iz(iat) = izph(iphat(iat))
            if (iphat(iat).eq.0) i0=iat
  5      continue

         write(slog,7)
   7     format(2x,'Calculating Debye-Waller factors via RM...')
         call wlog(slog)
         write(slog,9)
   9     format(2x,'This might take a while.')
         call wlog(slog)

c        Read spring.inp and build dynamical matrix
         call rdspr(rat1, iz, natom, i0, 
     1           dm, rnn, 
     1           acut, res, wmax, dosfit, zshell, w0,
     1           rintr, iprdos, nnl)

         if (ipath.ne.0) then
	    if(ir1_open) then
c             Echo title cards to s2_rm2.dat
              do 10  i = 1, ntitle
                write(ir1,12)  title(i)(1:ltit(i))
  10          continue
  12          format (1x, a)
              write(ir1,17) tk, natom
  17          format(1x,'temperature =',f7.2,2x,'N_at =',i4)
              write(ir1,19)
  19          format (1x, 71('-'))
              write(ir1,25)
              write(slog,25)
  25          format(1x,'ipath',2x,'nleg',4x,'sig2',3x,'mu_ipath',4x,
     1          'w_1',6x,'w_2',7x,'A1',5x,'A2')
              call wlog(slog)
            endif
            if (iprdos.ne.0.and.ir2_open) then
c              Echo title cards to s2_rm1.dat
               do 30  i = 1, ntitle
                  write(ir2,12)  title(i)(1:ltit(i))
  30           continue
               write(ir2,17) tk, natom
               write(ir2,19)
               write(ir2,35)
  35           format(1x,'ipath',2x,'nleg',4x,'sig2',3x,'mu_ipath',
     1              4x,'w_e')
            endif
         endif
      endif
      nsigc = nsigc + 1
c---- end of first time reading -------

cc    Open path input file (unit in) and read title.  Use unit 2.
c     ntitle2 = 5
c     open(unit=2,file='paths.dat',status='old', iostat=ios)
c     call chopen (ios, 'paths.dat', 'sigrm')
c     call rdhead (2, ntitle2, titlep, ltit)
cc    if (ntitle2 .le. 0)  then
cc       titlep(1) = ' '
cc    endif

c 84  continue
c     read(2,*,end=1010) ipath, nleg
c     skip label (x y z ipot rleg beta eta) and read the path
c     read(2,*)
      do 78 ileg=0,nleg
c        read(2,*,end=1010) (rat(j,ileg),j=1,3)
         nconv(ileg)=0
  78  continue

      do 88 n=1,3
         aa = rat(n,nleg)
         do 87 i=0,(nleg-2)
            j=nleg-i
            rat(n,j)=rat(n,j-1)
  87     continue
         rat(n,1)=aa
  88  continue
      do 89 i=1,nleg
         nq0(i)=0.
  89  continue

c nconv converts # of an atom in the nleg list of coordinates (rat) to
c its # in the full list of all atomic coordinates (rat1)
      do 94 i=1,natom
         do 91 n=1,3
            q1(n,i)=0.
   91    q0(n,i)=0.
            do 95 jl=1,nleg
            m=0
            do 93 n=1,3
               l=nint(100.*rat(n,jl))
               l1=nint(100.*rat1(n,i))
               if (abs(l-l1).le.1) m=m+1
   93       continue
            if (m.eq.3) then
              nconv(jl)=i
              go to 95
            endif
   95    continue
   94 continue

      atmu=0.
      inn=1
      nq1(inn)=1
      iq0=0
      nconv(0)=nconv(nleg)
      do 100 il=1,nleg
         l=nconv(il)
         do 101 jq=1,iq0
  101    if(nq0(jq).eq.l) go to 102
         iq0=iq0+1
         nq0(iq0)=l
  102    continue
         do 105 ii=1,natom
            a = dist(rat1(1,ii),rat1(1,l))
            if (a.le.rintr) then
               do 103 jn=1,inn
  103          if(nq1(jn).eq.ii) go to 105
               inn=inn+1
               nq1(inn)=ii
            endif
  105    continue
         nq1x=inn
         nq0x=iq0
         i=nconv(il)
         im=nconv(il-1)
         ip=nconv(il+1)
c        if (il.eq.1) im=nconv(nleg)
         if (il.eq.nleg) ip=nconv(1)
         atmass=atwtd(iz(i))
      do 100 n=1,3
         atmu=atmu + 0.25*( rnn(n,i,im)+rnn(n,i,ip) )**2 /atmass
  100 continue
      atmu=1./atmu
  108 continue

      do 115 i=1,natom
      do 115 n=1,3
  115 q0(n,i)=0.

c Build initial state vector |Q_j(0)> for the current path
      do 116 n=1,3
      do 116 il=1,nleg
         i=nconv(il)
         im=nconv(il-1)
         ip=nconv(il+1)
         if (il.eq.1) im=nconv(nleg)
         if (il.eq.nleg) ip=nconv(1)
         atmass=atwtd(iz(i))
         q0(n,i)=q0(n,i)+sqrt(atmu/atmass)*(rnn(n,im,i)-rnn(n,i,ip))/2.
  116 continue

c make sure it's normalized <Q_j(0)|Q_j(0)>=1
      q0q0=0.
      do 120 iq0=1,nq0x
      i=nq0(iq0)
      do 120 n=1,3
         q0q0=q0q0+q0(n,i)*q0(n,i)
 120  continue
      p00=nint(q0q0*1000.)/1000.
      if (abs(p00-1.d0).gt.5.d-4) then
         atmu=atmu/q0q0
         go to 108
      endif

c     to get THz units:
      wnorm=100.*w0/sqrt(amu0*10.) 
c*** moments
      a0=0.
      do 132 il=1,nq0x
      do 132 im=1,nq0x
         l=nq0(il)
         m=nq0(im)
      do 132 n1=1,3
      do 132 n2=1,3
         a0 = a0 + q0(n1,l)*dm(n1,n2,l,m)*q0(n2,m)/w0/w0
  132 continue
      we=wnorm*sqrt(a0)
      if (we.lt.1) then
c        recursion method is inapplicable, use statistics to set sig2
         sig2 = sig2x ( iphat(nconv(1)),  iphat(nconv(1)) )
         if (sig2.lt.1.d-6) sig2 = sig2mx
         return
      endif

      do 137 iset=1,nq1x
         i=nq1(iset)
      do 137 n1=1,3
            q1i=0.
            do 138 im=1,nq0x
               m=nq0(im)
            do 138 n2=1,3
               q1i=q1i+dm(n1,n2,i,m)*q0(n2,m)/w0/w0
  138       continue
         q1(n1,i) = q1i - a0*q0(n1,i)
  137 continue

      b0=0.
      do 139 i=1,natom
      do 139 n1=1,3
         b0=b0+q1(n1,i)*q1(n1,i)
  139 continue

      a1=0.
      do 150 iset=1,nq1x
         i=nq1(iset)
         do 150 n1=1,3
         q2=0.
         do 151 jset=1,nq1x
         j=nq1(jset)
            do 151 n2=1,3
               q2 = q2 + dm(n1,n2,i,j)*q1(n2,j)/w0/w0
  151       continue
            a1 = a1 + q1(n1,i)*q2
  150 continue

      a0=a0*wnorm**2
      a1=a1/b0
      a1=a1*wnorm**2
      b0=b0*wnorm**4

c** recursion sigma^2
      dd = (a0+a1)**2 - 4.*(a0*a1-b0)
      x1 = (a0+a1+sqrt(dd))/2.
      x2 = (a0+a1-sqrt(dd))/2.
      aa2 = (a1-x2)/(x1-x2)
c     aa2 = (a1-x2)/(x1-x2)*9./8.
      aa1 = (x1-a1)/(x1-x2)
      w1 = sqrt(x1)
      w2 = sqrt(x2)
      s1 = 3.1746/(atmu*w1*tanh(w1*7.6383/2./tk))
      s2 = 3.1746/(atmu*w2*tanh(w2*7.6383/2./tk))
      sigma2 = aa1*s1+aa2*s2
      sig2e = 3.1746/(atmu*we*tanh(we*7.6383/2./tk))

      if (ipath.ne.0) then
         write(slog,250) ipath,nleg,sigma2,atmu,w1,w2,aa1,aa2
         call wlog(slog)
         if (ir1_open) then
           write(ir1,250) ipath,nleg,sigma2,atmu,w1,w2,aa1,aa2
  250      format(1x,i3,4x,i1,3x,f9.5,2x,f7.3,2x,f7.2,2x,f7.2,
     1       4x,f5.3,2x,f5.3)
         endif
         if (iprdos.ne.0.and.ir2_open) then
            write(ir2,260) ipath,nleg,sig2e,atmu,we
  260       format(1x,i3,4x,i3,3x,f7.5,2x,f7.3,2x,f7.2)
         endif
      endif
      sig2 = sigma2

c     update maximum DW factors
      if (sig2.gt.sig2mx) sig2mx=sig2
      if (sig2.gt.sig2x( iphat(nconv(1)),  iphat(nconv(2)) )) then
         sig2x( iphat(nconv(1)),  iphat(nconv(2)) ) = sig2
         sig2x( iphat(nconv(2)),  iphat(nconv(1)) ) = sig2
      endif

      return
      end
c----------------------------------------------------
      subroutine sigte3 (iz1,iz2, sig2, alphat, thetad, reff, sig1,sig3)
c     single scattering only.

c     input: sig2
c     iz1, iz2 are iz at central atom and neighbor
c     alphat coeef of thermal expansion at high T
c     reff

c     output: sig1 sig3
      implicit double precision (a-h, o-z)
      real reff

c     con=hbar**2/kB*amu)*10**20   in ang**2 units
c     hbar = 1.054 572 666 e-34, amu = 1.660 540 e-27, 
c     kB = 1.380 6581 d-23
      parameter (con = 48.508 459 393 094)
      parameter (hbar = 1.054 572 666 e-34)
      parameter (amu = 1.660 540 e-27)
      parameter (xkb = 1.380 6581 e-23)

      ami=atwtd(iz1)*amu
      amj=atwtd(iz2)*amu

c     reduced mass
      xmu = 1 / (1/ami + 1/amj)
c     Einstein frequency
      omega = (2 * xkb * thetad) / (3 * hbar)
      xks = xmu * omega**2
      xk3 = xks**2 * reff * alphat / (3 * xkb)
      sig02 = hbar * omega / xks
      sig1 = -3 * (xk3 / xks) * sig2
      sig3 = 2 - (4.0/3.0) * (sig02 / sig2)**2
      sig3 = sig3 * (sig1 * sig2)

      return
      end
      subroutine sigm3(sig1, sig2, sig3, tk, alphat, thetae)
c     using correlated Einstein-model with a morse potential
c     Nguyen Van Hung & J.J.Rehr Phys. Rev. B 56 , 43

      implicit double precision (a-h, o-z)
      real sig02,  sig01, z
c     dimension alphat=[1/anstroems]
      parameter (bohr = 0.529 177 249d0)
      parameter(three= 3)
      parameter(four= 4 )
      parameter(fourthird= four/three)
      parameter(threequater= three/four)
              
      alphat= alphat * bohr  
      z=exp(- thetae/tk)
      sig02= (1-z)/ (1+z) * sig2
      sig01 = alphat * sig02 * threequater
      sig1 = sig01 * sig2 / sig02
      sig3 = (2- fourthird * (sig02/sig2) **2)* sig1 * sig2

      return
      end
c---------------------------------------------------------------------
c     program sigcl.f
c
c     calculates debye-waller factors for each multiple
c     scattering path using Debye-Model correlations
c
c     files:  input  pathd_all.dat  multiple scattering path data
c             output fort.3  sig**2 vs path
c                    fort.2  long output
c
c     version 1  (29 july 91)
c
c     coded by j. rehr
c     path data from s. zabinsky
c
c     modified to use pdata.inp, Dec 1991, siz
c     Subroutine version, Dec 1991, siz
c     Hacked for classical Debye Model, 7/2006, Kevin Jorissen  !KJ
c
c---------------------------------------------------------------------

      subroutine sigcl (tk, thetad, rs, nlegx, nleg, rat, iz, sig2)
c               tk temperature in degrees K
c               thetad debye temp in degrees K
c               rs=wigner seitz or norman radius in bohr, averaged
c                  over entire problem
c                  (4pi/3)*rs**3 = sum( (4pi/3)rnrm**3 ) / N
c                  (sum is over all atoms in the problem)
c               nlegx used in dimensions of rat and iz
c               nleg nlegs in path
c               rat positions of each atom in path
c               iz atomic number of each atom in path
c               NB Units of distance in this routine
c                  are angstroms, including sig**2.  rs is in bohr.
c               sig2 is output debye waller factor

      implicit double precision (a-h,o-z)

c={../HEADERS/const.h
      parameter (pi = 3.14159 26535 89793 23846 26433d0)
      parameter (one = 1, zero = 0)
      parameter (third = one/3)
      parameter (raddeg = 180 / pi)
      complex*16 coni
      parameter (coni = (0,1))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919 158 292 677 512 811d0)

      parameter (bohr = 0.529 177 249d0, ryd  = 13.605 698d0)
      parameter (hart = 2 * ryd)
      parameter (alpinv = 137.035 989 56d0)
c     fine structure alpha
      parameter (alphfs = 1 / alpinv)
c= ../HEADERS/const.h}

c     nlegx is max number of atoms in any one path
      dimension rat(3,0:nlegx)
      dimension iz(0:nlegx)
c#mn
       external dist

c      parameters
c               x = k_d*R   (distance parameter)
c               R distance in angstroms
c               y = hbar omegad/kT = thetad/t
c               thetad debye temp in degrees K
c               tk temperature in degrees K
c               k_d = (6*pi**2 N/V) = debye wave number
c               N/V=1/(4pi/3rs**3)
c               rs=wigner seitz or norman radius in bohr
c               ami, amj masses at sites i and j in amu
c               I = int_0^1 (y/x) dw sin(wx)coth(wy/2)

c     Note:  There are nleg atoms including the central atom
c            index 0 and index nleg both refer to central atom,
c            which makes special code unnecessary later.

      sigtot=0
      do 800 il=1,nleg
      do 800 jl=il,nleg

c        calculate r_i-r_i-1 and r_j-r_j-1

         rij = dist (rat(1,il), rat(1,jl))
         call corrfn2 (rij, cij, thetad, tk, iz(il), iz(jl), rs)
         sig2ij=cij

         rimjm = dist (rat(1,il-1), rat(1,jl-1))
         call corrfn2 (rimjm, cimjm, thetad, tk, iz(il-1), iz(jl-1), rs)
         sig2ij=sig2ij+cimjm

         rijm = dist (rat(1,il), rat(1,jl-1))
         call corrfn2 (rijm, cijm, thetad, tk, iz(il), iz(jl-1), rs)
         sig2ij=sig2ij-cijm

         rimj = dist (rat(1,il-1), rat(1,jl))
         call corrfn2 (rimj, cimj, thetad, tk, iz(il-1), iz(jl), rs)
         sig2ij=sig2ij-cimj

         riim = dist (rat(1,il), rat(1,il-1))
         rjjm = dist (rat(1,jl), rat(1,jl-1))

         ridotj=(rat(1,il)-rat(1,il-1))*(rat(1,jl)-rat(1,jl-1))+
     1          (rat(2,il)-rat(2,il-1))*(rat(2,jl)-rat(2,jl-1))+
     2          (rat(3,il)-rat(3,il-1))*(rat(3,jl)-rat(3,jl-1))
         ridotj=ridotj/(riim*rjjm)

c        double count i .ne. j  terms
         if(jl.ne.il) sig2ij=2*sig2ij
         sig2ij=sig2ij*ridotj
         sigtot=sigtot+sig2ij

  800 continue
      sig2=sigtot/4

c     sig2 is in bohr**2, just as we wanted for ff2chi
      return
      end



      subroutine corrfn2(rij,cij,thetad,tk,iz1,iz2,rsavg)
c     subroutine calculates correlation function
c     c(ri,rj)=<xi xj> in the Debye approximation
c
c             =(1/N)sum_k exp(ik.(Ri-Rj))(1/sqrt(mi*mj))*
c              (hbar/2w_k)*coth(beta hbar w_k/2)
c             = (3kT/mu w_d**2)*sqrt(mu**2/mi*mj)*I
c
c      parameters
c               x = k_d*R   (distance parameter)
c               R distance in angstroms
c               y = hbar omegad/kT = thetad/t
c               thetad debye temp in degrees K
c               tk temperature in degrees K
c               k_d = (6*pi**2 N/V) = debye wave number
c               N/V=1/(4pi/3rs**3)
c               rs=wigner seitz or norman radius in bohr
c               ami, amj masses at sites i and j in amu
c               I = int_0^1 (y/x) dw sin(wx)coth(wy/2)
c
c      solution by numerical integration
c
      implicit double precision (a-h, o-z)
      common /xy/ x, yinv

c={../HEADERS/const.h
      parameter (pi = 3.14159 26535 89793 23846 26433d0)
      parameter (one = 1, zero = 0)
      parameter (third = one/3)
      parameter (raddeg = 180 / pi)
      complex*16 coni
      parameter (coni = (0,1))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919 158 292 677 512 811d0)

      parameter (bohr = 0.529 177 249d0, ryd  = 13.605 698d0)
      parameter (hart = 2 * ryd)
      parameter (alpinv = 137.035 989 56d0)
c     fine structure alpha
      parameter (alphfs = 1 / alpinv)
c= ../HEADERS/const.h}

c     con=hbar**2/kB*amu)*10**20   in ang**2 units
c     hbar = 1.054 572 666 e-34, amu = 1.660 540 e-27, 
c     kB = 1.380 6581 d-23
      parameter (con = 48.508 459 393 094)
c#mn
       external atwtd

c     external fn
c     rij=2.55
c     tk=295
c     thetad=315
c     ami=amj=63.55 at wt for Cu
c     rs=2.7

      ami=atwtd(iz1)
      amj=atwtd(iz2)
      rs=rsavg
c     thetad in degrees K, t temperature in degrees K
c     y=thetad/tk
      yinv=tk/thetad
      xkd=(9*pi/2)**(third)/(rs*bohr)
      fac=(3/2.)*con/(thetad*sqrt(ami*amj))
      rj=rij
      x=xkd*rj
c     call numerical integration
      call bingrt2 (grater, eps, nx)
      cij=fac*grater
      return
      end
      double precision function fn2(w)
      implicit double precision (a-h,o-z)
      common/xy/x,yinv
c     fn=(sin(wx)/x)*coth(wy/2)
c     change code to allow t=0 without bombing
c     fn=2/y
      fn2=2*yinv
      if(w.lt.1.e-20) return
      fac=w
      if(x.gt.0.) fac=sin(w*x)/x
!      emwy=0.
!      if(yinv.gt.0.0125) emwy=exp(-w/yinv)
!      emwy=exp(-w/yinv)
!      fn2=fac*(1+emwy)/(1-emwy)
      fn2=fac*2*yinv/w  !KJ coup de grace
      return
      end
c-----------------------------------------------
      subroutine bingrt2 (b, eps, n)
c     subroutine calculates integrals between [0,1]
c      b = int_0^1 f(z) dz
c     by trapezoidal rule and binary refinement
c     (romberg integration)
c     coded by j rehr (10 Feb 92)
c     see, e.g., numerical recipes for discussion
c     and a much fancier version
c-----------------------------------------------
c     del=dz  itn=2**n tol=1.e-5
c     starting values
      implicit double precision (a-h,o-z)
      common /xy/x,yinv
      character*512 slog
c     external fn
c     error is approximately 2**(-2n) ~ 10**(-.6n)
c     so nmax=10 implies an error of 1.e-6
      parameter(nmax = 10, tol = 1.e-5)
      parameter(zero=0, one=1)
      n=0
      itn=1
      del=1.
      bn=(fn2(zero)+fn2(one))/2
      bo=bn
 10   continue
c     nth iteration
c     b_n+1=(b_n)/2+deln*sum_0^2**n f([2n-1]deln)
      n=n+1
      if(n.gt.nmax) go to 40
      del=del/2
      sum=0.
      do 20 i=1, itn
      zi=(2*i-1)*del
 20   sum=sum+fn2(zi)
c     bnp1=b_n+1 is current value of integral
      bnp1=bn/2+del*sum
c     cancel leading error terms b=[4b-bn]/3
c     note: this is the first term in the
c     neville table - remaining errors were
c     found too small to justify the added code
      b=(4*bnp1-bn)/3
      eps=abs((b-bo)/b)
      if(eps.lt.tol) goto 60
      bn=bnp1
      bo=b
      itn=itn*2
      goto 10
 40   write(slog,50) n,itn, b,eps
      call wlog(slog)
 50   format(' not converged, n,itn,b,eps=',
     1  2i4,2e14.6)
      return
 60   continue
      return
      end
c///////////////////////////////////////////////////////////////////////
c Distribution:  FEFF_MATH 1.0
c Copyright (c) [2002] University of Washington
c 
c This software was prepared in part with US Government Funding under
c DOE contract DE-FG03-97ER45623.

c Redistribution and use of this Distribution in source and binary
c formats, with or without modification is permitted, provided the 
c following conditions are met:
c 
c Redistributions must retain the above notices and the following list
c of conditions and disclaimer;
c 
c Modified formats carry the marking
c     "Based on or developed using Distribution: FEFF_MATH 1.0
c      FEFF_MATH 1.0 Copyright (c) [2002] University of Washington"
c 
c Recipient acknowledges the right of the University of Washington to
c prepare uses of this Distribution and its modifications that may be
c substantially similar or functionally equivalent to
c Recipient-prepared modifications.
c
c Recipient and anyone obtaining access to the Distribution through
c recipient's actions accept all risk associated with possession and
c use of the Distribution.
c
c THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESS OR IMPLIED
c WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
c MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
c IN NO EVENT SHALL THE UNIVERSITY OF WASHINGTON OR CONTRIBUTORS TO THE
c DISTRIBUTION BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
c EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
c PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
c REVENUE; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
c LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
c NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
c SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
c///////////////////////////////////////////////////////////////////////
c License is applicable for routines below, until otherwise specified.
c
      subroutine bcoef(kinit, ipol, ptz, le2, ltrace, ispin, angks,
     1                 kind, lind, bmat)
c     written by alexei ankudinov; march 2000
c     calculate bmat: the energy independent sum over polarization and
c     angular momenta indices
c     bmat = \sum_{p,p', all m_j} <LS|J><J|R|J1><J1|\alpha_p exp(i kz)|I>
c                    ptz(p,p') 
c            <I|\alpha_p'^* exp(-i kz) J2><J2|R'|J'><J'|L'S'>
c     where R is rotation from spin vector to x-ray k-vector
c     and R' is rotation back
c     see Eq.10 and 11 in Ankudinov,Rehr, Phys.Rev.B (accepted),
c     Theory of solid state contribution to the x-ray elastic scattering
c     aditional rotation matrices are needed when x-ray k-vector
c     is not along the spin-axis (see rotations in rdinp)

c     more precisely it is
c     bmat(l1 l1' j l ml ms; l2 l2' j' l' ml' ms') =
c        (-)**(j-j'+l2'+1) i**(l'-l) \sum_{p,p',mi,m1,mj,m2,mj'}
c        <LS|J>   r^j_{m1,mj}(angks)   3j( j l1 i ; -m1 p mi)
c        (-p)**(l1+l1'+1) ptz(p,p') (-p')**(l2+l2'+1) 
c        3j( j' l2 i ; -m2  p' mi)   r^j'_{m2,mj'}(angks)   <J'|L'S'>
c     where l1 l1' are set by the multipole moment(E1-l1=1,l1'=0;
c     E2-l1=2,l1'=1; M1-l1=1,l1'=1; etc.;
c     j and l define quantum number kappa and for each multipole moment
c     Only few final kappa's are allowed and  it is convinient
c     to denote (l1 l1' j l) by one index 'k'
c     thus  k=1-8 to include both E1 and E2 transitions;
c     ml and ms are projections of orbital and spin moments.

c     bmat  is used to calculate absorption fine structure (chi) via
c       chi = \sum_{k ms,k' ms'}  rkk(k,ms)  rkk(k',ms')
c       \sum_{ml,ml'}  bmat(k ml ms; k' ml' ms')  G_(l' ml' ms';l ml ms)
c     where sum over spins can be moved from first sum to second for
c     spin independent systems. The above expression is suitable for FMS
c     and for MS expansion on can use Eq.15 in RA paper to obtain
c     expression for the termination   matrix
c     T_{lam1 ms,lamN ms'} = \sum_{k k'} rkk(k,ms) rkk(k',ms')
c       \sum_{ml,ml'}  bmat(k ml ms; k' ml' ms') gam(l,lam1,rho1,ms)
c        gamtl(l',lamN,rhoN,ms')
c     Notice that for spin-dependent systems the scattering F matrices
c     in RA paper also should have additional spin indices. In genfmt
c     we currently neglect spin-flip processes which simplifies
c     calculations with MS expansion. (T and F are diagonal in ms,ms')
       
c     This subroutine is written for general spin-dependent asymmetric
c     system and arbitrary polarization tenzor. The symmetry of the 
c     system and polarization tenzor can be used
c     to speed up FMS or MS calculations in appropriate subroutines.
c     (see comments in subroutines mpprmp, fmstot)

c     input:
c       kinit - kappa for initial orbital
c       ipol - polarization type measurement
c       ptz  - polarization tensor (needed only for ipol=1 case)
c       le2  - sets which multipole moments to include (see mkptz)
c       ltrace- .true. for xsect.f, where need to perform trace over ml
c       angks - angle between k-vector and spin-vector 

c     output
c       lind  - orb.mom.(kappa)  needed in fmstot only (for indexing)
c       bmat  - energy independent matrix to calculate absorption 
c       in many cases bmat is diagonal due to the choice of xyz frame,
c       but for general case full 16*(2*lx+1)*16*(2*lx+1) matrix is kept

      implicit double precision (a-h,o-z)
c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}
      complex*16 coni
      parameter (coni = (0,1))

c     need only parameter lx to set max orb momentum
      complex*16 ptz, bmat, pmat, tmat
      dimension ptz(-1:1,-1:1),  bmat(-lx:lx,0:1,8, -lx:lx,0:1,8)
c       to include all possible dipole and quadrupole transitions 
c       final kp, and kpp have 8 possibilities
      logical ltrace

c     local staff
      dimension  t3j( 8, 0:1, -lx:lx+1), x3j(8, -1:1, -lx:lx+1)
c     qmat = <J2|R'|J'><J'|L'S'> - diagonal in kappa index
      dimension qmat( -lx:lx+1, -lx:lx, 0:1, 8)
c     pmat = <J1|\alpha_j exp(i kz)|I> ptz <I|\alpha_k^* exp(-i kz)|J2>
      dimension pmat( -lx:lx+1, 8, -lx:lx+1, 8)
c     tmat = pmat*qmat ; bmat = qmat^T*tmat
      dimension tmat( -lx:lx+1, 8, -lx:lx, 0:1, 8)
c     total and orbital momenta for 8 possible final kappa
      dimension jind(8), lind(8), kind(8)

      external cwig3j

      do 10 i6 = 1, 8
      do 10 i5 = 0 ,1
      do 10 i4 = -lx,lx
      do 10 i3 = 1, 8
      do 10 i2 = 0 ,1
      do 10 i1 = -lx,lx
         bmat( i1, i2, i3, i4, i5, i6) = 0
  10  continue

c     3 dipole transitions
      do 20 k=-1,1
         kap=kinit+k
         if (k.eq.0) kap=-kap
         jkap = abs(kap)
         lkap = kap
         if (kap.le.0) lkap = abs(kap) -1
c        check that orbital momentum does not exceed max allowed
         if (lkap .gt. lx) then
c          set final j and l to unphysical values
           jkap = 0
           lkap = -1 
           kap = 0
         endif
         jind(k+2) = jkap
         lind(k+2) = lkap
         kind(k+2) = kap
  20  continue

c     include 5 quadrupole or 3 mag.dipole  transitions
      do 120 k=-2,2
         jkap = abs(kinit) + k
         if (jkap.le.0) jkap = 0
         kap= jkap
         if (kinit.lt.0 .and. abs(k).ne.1) kap=-jkap
         if (kinit.gt.0 .and. abs(k).eq.1) kap=-jkap
         lkap = kap
         if(kap.le.0) lkap = - kap - 1
         if (lkap.gt.lx .or. le2.eq.0
     1                  .or. (le2.eq.1 .and. abs(k).eq.2)) then
c           set unphysical jkap and lkap to make shorter calculations
            jkap = 0
            lkap = -1
            kap = 0
         endif
         jind(k+6) = jkap
         lind(k+6) = lkap
         kind(k+6) = kap
 120  continue

      if (ipol.eq.0) then
c       polarization average case; bmat is diagonal and simple
        do 100 k = 1, 8
        do 100 ms = 0 ,1
        do 100 ml = -lind(k), lind(k)
c         i2 = (2*l1+1) , where l1 is defined by multipole moment
          i2 = 3
          if (le2.eq.2 .and. k.gt.3) i2 = 5
          bmat(ml,ms,k, ml,ms,k) = 0.5d0 / (2*lind(k)+1.d0) / i2
          if (k.le.3) bmat(ml,ms,k, ml,ms,k) = - bmat(ml,ms,k, ml,ms,k)
 100    continue
      else
c       more complicated bmat for linear(ipol=1) and circular(ipol=2)
c       polarizations
c       Put 3j factors in x3j and t3j. t3j are multiplied by
c       sqrt(2*j'+1) for  further convinience.
        do 30  mp=-lx,lx+1
        do 30  ms=0,1
        do 30  k1=1,8
  30    t3j(k1,ms,mp) = 0.0d0
        do 40  mp=-lx,lx+1
        do 40  ms=-1,1
        do 40  k1=1,8
  40      x3j(k1,ms,mp) = 0.0d0

        do 70  k1 = 1,8
        do 70  mp = -jind(k1)+1,jind(k1)
          do 50 ms=0,1
            j1 = 2 * lind(k1)
            j2 = 1
            j3 = 2 * jind(k1) - 1
            m1 = 2*(mp-ms)
            m2 = 2*ms - 1
            t3j(k1,ms,mp)=sqrt(j3+1.0d0) * cwig3j(j1,j2,j3,m1,m2,2)
            if (mod( (j2-j1-m1-m2)/2 , 2) .ne.0) 
     1          t3j(k1,ms,mp) = - t3j(k1,ms,mp)
c           t3j(m0,i)    are Clebsch-Gordon coefficients
  50      continue
          do 60 i=-1,1
            j1 = 2 * jind(k1) - 1
            j2 = 2
            if (k1.gt.3 .and. le2.eq.2) j2 = 4
            j3 = 2 * abs(kinit) - 1
            m1 = -2*mp + 1
            m2 = 2*i
            x3j(k1,i,mp)= cwig3j(j1,j2,j3,m1,m2,2)
  60      continue
  70    continue

c       calculate qmat
        do 220 i=1,8
        do 220 ms=0,1
        do 220 ml= -lind(i), lind(i)
        do 220 mj= -jind(i)+1, jind(i)
          mp = ml+ms
          jj = 2*jind(i) - 1
          mmj = 2*mj - 1
          mmp = 2*mp - 1
          value = rotwig(angks, jj, mmj, mmp, 2)
          qmat(mj,ml,ms,i) = value * t3j(i,ms,mp)
 220    continue

c       calculate pmat
        do 240 i2 = 1,8
        do 240 m2 = -jind(i2)+1, jind(i2)
        do 240 i1 = 1,8
        do 240 m1 = -jind(i1)+1, jind(i1)
          pmat(m1,i1,m2,i2) = 0
          if (abs(m2-m1).le.2) then
            do 230 j=-1,1
            do 230 i=-1,1
c             check that initial moment is the same
              if (m1-i.eq.m2-j) then
                is = 1
c               (-p) factors for M1 transitions
                if (le2.eq.1 .and. i.gt.0 .and. i1.gt.3) is = -is
                if (le2.eq.1 .and. j.gt.0 .and. i2.gt.3) is = -is
                pmat(m1,i1,m2,i2) = pmat(m1,i1,m2,i2) +
     1          is * x3j(i1,i,m1) * ptz(i,j) * x3j(i2,j,m2)
              endif
 230        continue
c           multiply by (-)^(j-j'+l2'+1) i**(l'-l) factor
c           additional (-) is from Eq.10 (-2*ck)
            is = 1
            if (mod(jind(i1)-jind(i2), 2) .ne.0) is = -is
            if (i2.le.3) is = -is
            pmat(m1,i1,m2,i2) = pmat(m1,i1,m2,i2) * is
     1           * coni**(lind(i2)-lind(i1))
          endif
 240    continue

c       calculate tmat = pmat*qmat
        do 270 i1=1,8
        do 270 ms=0,1
        do 270 ml=-lind(i1), lind(i1)
        do 270 i2=1,8
        do 270 mj=-jind(i2)+1, jind(i2)
          tmat(mj,i2, ml,ms,i1) = 0
          do 260 mp = -jind(i1)+1, jind(i1)
            tmat(mj,i2, ml,ms,i1) = tmat(mj,i2, ml,ms,i1)+
     1           pmat(mj,i2,mp,i1) * qmat(mp,ml,ms,i1)
 260      continue
 270    continue
         
c       calculate bmat = qmat^T * tmat
        do 300 i1=1,8
        do 300 ms1=0,1
        do 300 ml1=-lind(i1), lind(i1)
        do 300 i2=1,8
        do 300 ms2=0,1
        do 300 ml2=-lind(i2), lind(i2)
          bmat(ml2,ms2,i2, ml1,ms1,i1) = 0
          do 280 mj=-jind(i2)+1, jind(i2)
            bmat(ml2,ms2,i2, ml1,ms1,i1) = bmat(ml2,ms2,i2, ml1,ms1,i1)+
     1      qmat(mj,ml2,ms2,i2) * tmat(mj,i2,ml1,ms1,i1) 
 280      continue
 300    continue
c       end of ipol=1,2 cases
      endif 

      if (ltrace) then
c       need to trace bmat over ml for xsect.f
        do 390 i1 = 1, 8
        do 390 ms1 = 0,1
        do 390 i2 = 1, 8
        do 390 ms2 = 0,1
          if (lind(i1).ne.lind(i2) .or. ms1.ne.ms2) then
               bmat(0,ms2,i2, 0,ms1,i1) = 0
          else
             do 360 ml = 1, lind(i1)
               bmat(0,ms1,i2, 0,ms1,i1) =  bmat(0,ms1,i2, 0,ms1,i1) +
     1         bmat(-ml,ms1,i2, -ml,ms1,i1) + bmat(ml,ms1,i2, ml,ms1,i1)
 360         continue
          endif
 390    continue
      endif

      if (ispin .eq. 0) then
c       G(Ls,L's') is spin diagonal; trace over spin
        do 480 i1 = 1, 8
        do 480 i2 = 1, 8
        do 480 ml1 = -lind(i1), lind(i1)
        do 480 ml2 = -lind(i2), lind(i2)
           bmat(ml2,0,i2, ml1,0,i1) =   bmat(ml2,0,i2, ml1,0,i1) +
     1                                  bmat(ml2,1,i2, ml1,1,i1)
 480    continue
      elseif (ispin.eq.2 .or. (ispin.eq.1 .and. nspx.eq.1)) then
c       move spin up part into the position of spin-down
        do 490 i1 = 1, 8
        do 490 i2 = 1, 8
        do 490 ml1 = -lind(i1), lind(i1)
        do 490 ml2 = -lind(i2), lind(i2)
           bmat(ml2,0,i2, ml1,0,i1) =   bmat(ml2,1,i2, ml1,1,i1)
 490    continue

      endif

      return
      end
      subroutine besjn (x, jl, nl)

c-----------------------------------------------------------------------
c
c     purpose:  to calculate the spherical bessel functions jl and nl
c               for l = 0 to 30 (no offset)
c
c     arguments:
c       x = argument of jl and nl
c       jl = jl bessel function (abramowitz conventions)
c       nl = nl bessel function (abramowitz yl conventions)
c            Note that this array nl = abramowitz yl.
c       jl and nl must be dimensioned 
c            complex*16 jl(ltot+2), nl(ltot+2), with ltot defined in 
c            dim.h.
c
c     notes:  jl and nl should be calculated at least to 10 place
c             accuracy for the range 0<x<100 according to spot
c             checks with tables
c
c     error messages written with PRINT statement.
c
c     first coded by r. c. albers on 14 dec 82
c
c     version 3
c
c     last modified: 27 jan 83 by r. c. albers
c     dimension of jl,nl changed from 31 to 26  (10 aug 89) j. rehr
c     modified again, siz, June 1992
c
c-----------------------------------------------------------------------

      implicit double precision (a-h, o-z)
c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}

      complex*16 x
      complex*16 jl(ltot+2), nl(ltot+2)
      complex*16 cjl(ltot+2), sjl(ltot+2), cnl(ltot+2), snl(ltot+2)

      complex*16 xjl,xnl,asx,acx
      complex*16 xi,xi2,xi3,xi4,xi5,xi6,xi7,xi8,xi9,xi10,xi11

      parameter (xcut = 1.d0, xcut1 = 7.51d0, xcut2 = 5.01d0)

      if (dble(x) .le. 0)  stop 'Re(x) is .le. zero in besjn'

      lmaxp1 = ltot+2

      if (dble(x) .lt. xcut .and. abs(dimag(x)) .lt. xcut)  then
c        case Re(x) < 1, just use series expansion
         do 10 il = 1,lmaxp1
            l = il-1
            ifl = 0
            call bjnser (x,l,xjl,xnl,ifl)
            jl(il) = xjl
            nl(il) = xnl
   10    continue

      elseif (dble(x) .lt. xcut1 .and. abs(dimag(x)) .lt. xcut1)  then

c        case 1 <= Re(x) < 7.5

         call bjnser (x,lmaxp1-1,xjl,xnl,1)
         jl(lmaxp1) = xjl

         call bjnser (x,lmaxp1-2,xjl,xnl,1)
         jl(lmaxp1-1) = xjl

         if (dble(x) .lt. xcut2 .and. abs(dimag(x)) .lt. xcut2)  then
c           Re(x) < 5
            call bjnser (x,0,xjl,xnl,2)
            nl(1) = xnl
            call bjnser (x,1,xjl,xnl,2)
            nl(2) = xnl
         else
c           Re(x) >= 5
            asx = sin(x)
            acx = cos(x)
            xi = 1 / x
            xi2 = xi**2
            nl(1) = -acx*xi
            nl(2) = -acx*xi2 - asx*xi
         endif

c        Use recursion relation 10.1.19 to get nl and jl
         do 50 lp1 = 3, lmaxp1
            l = lp1 - 2
            tlxp1 = 2*l + 1
            nl(lp1) = tlxp1 * nl(lp1-1) / x  -  nl(lp1-2)
   50    continue

         do 60 lxx = 3,lmaxp1
            lp1 = lmaxp1+1-lxx
            l = lp1-1
            tlxp3 = 2*l + 3
            jl(lp1) = tlxp3 * jl(lp1+1) / x  -  jl(lp1+2)
   60    continue

      else
c        case Re(x) > 7.5
c        Use AS 10.1.8 and 10.1.9, sjl=P, qjl=Q, note that AS formulae
c        use cos (z - n*pi/2), etc., so cos and sin terms get a bit
c        scrambled (mod 4) here, since n is integer.  These are hard-
c        coded into the terms below.
         xi = 1 / x
         xi2  = xi*xi
         xi3  = xi*xi2
         xi4  = xi*xi3
         xi5  = xi*xi4
         xi6  = xi*xi5
         xi7  = xi*xi6
         xi8  = xi*xi7
         xi9  = xi*xi8
         xi10 = xi*xi9
         xi11 = xi*xi10

         sjl(1) = xi
         sjl(2) = xi2
         sjl(3) = 3*xi3 - xi
         sjl(4) = 15*xi4 - 6*xi2
         sjl(5) = 105*xi5 - 45*xi3 + xi
         sjl(6) = 945*xi6 - 420*xi4 + 15*xi2
         sjl(7) = 10395*xi7 - 4725*xi5 + 210*xi3 - xi
         sjl(8) = 135135*xi8 - 62370*xi6 + 3150*xi4 - 28*xi2
         sjl(9) = 2027025*xi9 - 945945*xi7 + 51975*xi5 
     1            - 630*xi3 + xi
         sjl(10) = 34459425*xi10 - 16216200*xi8 + 945945*xi6 
     1            - 13860*xi4 + 45*xi2
         sjl(11) = 654729075*xi11 - 310134825*xi9 + 18918900*xi7 
     1            - 315315*xi5 + 1485*xi3 - xi
         cjl(1) = 0
         cjl(2) = -xi
         cjl(3) = -3*xi2
         cjl(4) = -15*xi3 + xi
         cjl(5) = -105*xi4 + 10*xi2
         cjl(6) = -945*xi5 + 105*xi3 - xi
         cjl(7) = -10395*xi6 + 1260*xi4 - 21*xi2
         cjl(8) = -135135*xi7 + 17325*xi5 - 378*xi3 + xi
         cjl(9) = -2027025*xi8 + 270270*xi6 - 6930*xi4 + 36*xi2
         cjl(10) = -34459425*xi9 + 4729725*xi7 - 135135*xi5 
     1             + 990*xi3 - xi
         cjl(11) = -654729075*xi10 + 91891800*xi8 - 2837835*xi6 
     1             + 25740*xi4 - 55*xi2
         do 80 ie = 1,11
            snl(ie) = cjl(ie)
            cnl(ie) = -sjl(ie)
   80    continue
         do 90 lp1 = 12,lmaxp1
            l = lp1-2
            tlxp1 = float(2*l+1)
            sjl(lp1) = tlxp1*xi*sjl(lp1-1)-sjl(lp1-2)
            cjl(lp1) = tlxp1*xi*cjl(lp1-1)-cjl(lp1-2)
            snl(lp1) = tlxp1*xi*snl(lp1-1)-snl(lp1-2)
            cnl(lp1) = tlxp1*xi*cnl(lp1-1)-cnl(lp1-2)
   90    continue
         asx = sin(x)
         acx = cos(x)
         do 110 lp1 = 1,lmaxp1
            jl(lp1) = asx*sjl(lp1)+acx*cjl(lp1)
            nl(lp1) = asx*snl(lp1)+acx*cnl(lp1)
  110    continue
      endif

      return
      end
      subroutine besjh (x, lbmax, jl, hl)

c-----------------------------------------------------------------------
c
c     purpose:  to calculate the spherical bessel functions jl and hl
c               for l = 0 to lbmax (no offset)
c
c     arguments:
c       x = argument of jl and nl
c       lbmax
c       jl = jl bessel function (abramowitz conventions)
c       hl = hl^+ bessel function (messiah conventions) for Im x >=0
c       hl = hl^- bessel function (messiah conventions) for Im x < 0
c       jl and hl must be dimensioned 
c            complex*16 jl(0:lbmax), hl(0:lbmax), 
c
c     notes:  jl and hl should be calculated at least to 10 place
c             accuracy for the range 0<x<100 according to spot
c             checks with tables
c
c     error messages written with PRINT statement.
c
c     first coded by r. c. albers on 14 dec 82
c
c     version 3
c
c     last modified: 27 jan 83 by r. c. albers
c     dimension of jl,nl changed from 31 to 26  (10 aug 89) j. rehr
c     modified again, siz, June 1992
c     rewritten for jl and hl by a.l. ankudinov feb 2000
c
c-----------------------------------------------------------------------

      implicit double precision (a-h, o-z)
c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}

      complex*16 x
      complex*16 jl(0:lbmax), nl(ltot+2)
      complex*16 hl(0:lbmax)
      complex*16 cjl(ltot+2), sjl(ltot+2)

      complex*16 xjl,xnl,asx,acx, epx
      complex*16 xi,xi2,xi3,xi4,xi5,xi6,xi7,xi8,xi9,xi10,xi11

      parameter (xcut = 1.d0, xcut1 = 7.51d0, xcut2 = 5.01d0)
      complex*16 coni
      parameter (coni=(0,1))

      if (dble(x) .lt. 0)  stop 'Re(x) is .lt. zero in besjh'

      lmax = min(lbmax, ltot+1)
      lmaxp1 = lmax + 1

      if (dble(x) .lt. xcut .and. abs(dimag(x)) .lt. xcut)  then
c        case Re(x) < 1, just use series expansion
         do 10 ll = 0,lmax
            ifl = 0
            call bjnser (x,ll,xjl,xnl,ifl)
            jl(ll) = xjl
            hl(ll) = -xnl + coni*xjl
   10    continue

      elseif (dble(x) .lt. xcut1 .and. abs(dimag(x)) .lt. xcut1)  then

c        case 1 <= Re(x) < 7.5

         call bjnser (x,lmax,xjl,xnl,1)
         jl(lmax) = xjl

         call bjnser (x,lmax-1,xjl,xnl,1)
         jl(lmax-1) = xjl

         if (dble(x) .lt. xcut2 .and. abs(dimag(x)) .lt. xcut2)  then
c           Re(x) < 5
            call bjnser (x,0,xjl,xnl,2)
            nl(1) = xnl
            call bjnser (x,1,xjl,xnl,2)
            nl(2) = xnl
         else
c           Re(x) >= 5
            asx = sin(x)
            acx = cos(x)
            xi = 1 / x
            xi2 = xi**2
            nl(1) = -acx*xi
            nl(2) = -acx*xi2 - asx*xi
         endif

c        Use recursion relation 10.1.19 to get nl and jl
         do 50 lp1 = 3, lmaxp1
            l = lp1 - 2
            tlxp1 = 2*l + 1
            nl(lp1) = tlxp1 * nl(lp1-1) / x  -  nl(lp1-2)
   50    continue

         do 60 lxx = 3,lmaxp1
            lp1 = lmaxp1+1-lxx
            l = lp1-1
            tlxp3 = 2*l + 3
            jl(l) = tlxp3 * jl(l+1) / x  -  jl(l+2)
   60    continue

         do 65 il = 1, lmaxp1
            l = il - 1
            hl(l) = -nl(il) + coni*jl(l)
   65    continue

      else
c        case Re(x) > 7.5
c        Use AS 10.1.8 and 10.1.9, sjl=P, qjl=Q, note that AS formulae
c        use cos (z - n*pi/2), etc., so cos and sin terms get a bit
c        scrambled (mod 4) here, since n is integer.  These are hard-
c        coded into the terms below.
         xi = 1 / x
         xi2  = xi*xi
         xi3  = xi*xi2
         xi4  = xi*xi3
         xi5  = xi*xi4
         xi6  = xi*xi5
         xi7  = xi*xi6
         xi8  = xi*xi7
         xi9  = xi*xi8
         xi10 = xi*xi9
         xi11 = xi*xi10

         sjl(1) = xi
         sjl(2) = xi2
         sjl(3) = 3*xi3 - xi
         sjl(4) = 15*xi4 - 6*xi2
         sjl(5) = 105*xi5 - 45*xi3 + xi
         sjl(6) = 945*xi6 - 420*xi4 + 15*xi2
         sjl(7) = 10395*xi7 - 4725*xi5 + 210*xi3 - xi
         sjl(8) = 135135*xi8 - 62370*xi6 + 3150*xi4 - 28*xi2
         sjl(9) = 2027025*xi9 - 945945*xi7 + 51975*xi5 
     1            - 630*xi3 + xi
         sjl(10) = 34459425*xi10 - 16216200*xi8 + 945945*xi6 
     1            - 13860*xi4 + 45*xi2
         sjl(11) = 654729075*xi11 - 310134825*xi9 + 18918900*xi7 
     1            - 315315*xi5 + 1485*xi3 - xi
         cjl(1) = 0
         cjl(2) = -xi
         cjl(3) = -3*xi2
         cjl(4) = -15*xi3 + xi
         cjl(5) = -105*xi4 + 10*xi2
         cjl(6) = -945*xi5 + 105*xi3 - xi
         cjl(7) = -10395*xi6 + 1260*xi4 - 21*xi2
         cjl(8) = -135135*xi7 + 17325*xi5 - 378*xi3 + xi
         cjl(9) = -2027025*xi8 + 270270*xi6 - 6930*xi4 + 36*xi2
         cjl(10) = -34459425*xi9 + 4729725*xi7 - 135135*xi5 
     1             + 990*xi3 - xi
         cjl(11) = -654729075*xi10 + 91891800*xi8 - 2837835*xi6 
     1             + 25740*xi4 - 55*xi2
         do 90 lp1 = 12,lmaxp1
            l = lp1-2
            tlxp1 = float(2*l+1)
            sjl(lp1) = tlxp1*xi*sjl(lp1-1)-sjl(lp1-2)
            cjl(lp1) = tlxp1*xi*cjl(lp1-1)-cjl(lp1-2)
   90    continue
         asx = sin(x)
         acx = cos(x)
         if (dimag(x).ge. 0.d0) then
           epx = exp(coni*x)
         else 
           epx = exp(-coni*x)
         endif
         do 110 ll = 0,lmax
            lp1 = ll + 1
            jl(ll) = asx*sjl(lp1)+acx*cjl(lp1)
            if (dimag(x).ge. 0.d0) then
              hl(ll) = (sjl(lp1)+coni*cjl(lp1)) * epx
            else
              hl(ll) = (sjl(lp1)-coni*cjl(lp1)) * epx
            endif
  110    continue
      endif

      return
      end
      subroutine bjnser (x, l, jl, nl, ifl)

c-----------------------------------------------------------------------
c
c     subroutine: bjnser (x,l,jl,nl,ifl)
c
c     purpose:  to calculate the spherical bessel functions jl and nl
c
c     arguments:
c       x = argument of jl and nl
c       l = l value calculated (no offset)
c       jl = jl bessel function (abramowitz conventions)
c       nl = nl bessel function (abramowitz yl conventions)
c       ifl = 0 return both jl and nl
c             1 return jl only
c             2 return nl only
c
c     notes:  jl and nl are calculated by a series
c             expansion according to 10.1.2 and 10.1.3
c             in abramowitz and stegun (ninth printing),
c             page 437
c
c             error msgs written with PRINT statements.
c
c     first coded by r. c. albers on 26 jan 83
c
c     version 2
c
c     last modified: 27 jan 83 by r. c. albers
c
c-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

      complex*16 x,u,ux,del,pj,pn
      complex*16 jl,nl

      character*512 slog

      parameter (niter = 160, tol = 1.e-15)

      if (l .lt. 0) then
         call wlog(' l .lt. 0 in bjnser')
         stop 'bjnser 1'
      endif
      if (dble(x).lt. 0.) then
         write(slog,30) x
         call wlog(slog)
   30    format (' x = ', 1p, 2e14.6, ' is .le. 0 in bjnser')
         stop 'bjnser 2'
      endif

      lp1 = l+1
      u = x**2 / 2

c     make djl = 1 * 3 * 5 * ... * (2*l+1),
c          dnl = 1 * 3 * 5 * ... * (2*l-1)
      djl = 1
      fac = -1
      do 50 il = 1, lp1
         fac = fac + 2
         djl = fac * djl
   50 continue
      dnl = djl / (2*l+1)


      if (ifl .eq. 2)   goto 90
c     make jl
c     pj is term in { } in 10.1.2, del is last factor in the series
c     convergence test is (last factor)/(total term) <= tol
      pj = 1
      nf = 1
      nfac = 2*l + 3
      den = nfac
      sgn = -1
      ux = u
      do 60 il = 1, niter
         del = sgn*ux / den
         pj = pj + del
         trel = abs (del / pj)
         if (trel .le. tol)  goto 80
         sgn = -sgn
         ux = u*ux
         nf = nf+1
         nfac = nfac+2
         den = nf * nfac * den
   60 continue
      stop  'jl does not converge in bjnser'
   80 jl = pj * (x**l) / djl

   90 if (ifl.eq.1) return
c     make nl
c     pn is term in { } in 10.1.3, del is last factor in the series
c     convergence test is (last factor)/(total term) <= tol
      pn = 1
      nf = 1
      nfac = 1 - 2*l
      den = nfac
      sgn = -1
      ux = u
      do 100  il = 1, niter
         del = sgn * ux / den
         pn = pn + del
         trel = abs (del / pn)
         if (trel .le. tol) goto 120
         sgn = -sgn
         ux = u*ux
         nf = nf+1
         nfac = nfac+2
         den = nf * nfac * den
  100 continue
      stop  'nl does not converge in bjnser'
  120 nl = -pn * dnl / (x**lp1)

      return
      end
      subroutine conv(omega,xsec,ne1,vicorr)
c     multiply xsec by theta(omega-efermi) and
c     convolute xsec(omega) with  xloss/((omega-omega0)**2+xloss**2)/pi
c     the result is xsec0(omega0)

      implicit double precision (a-h, o-z)
c={../HEADERS/const.h
      parameter (pi = 3.14159 26535 89793 23846 26433d0)
      parameter (one = 1, zero = 0)
      parameter (third = one/3)
      parameter (raddeg = 180 / pi)
      complex*16 coni
      parameter (coni = (0,1))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919 158 292 677 512 811d0)

      parameter (bohr = 0.529 177 249d0, ryd  = 13.605 698d0)
      parameter (hart = 2 * ryd)
      parameter (alpinv = 137.035 989 56d0)
c     fine structure alpha
      parameter (alphfs = 1 / alpinv)
c= ../HEADERS/const.h}
c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}
      dimension  omega(nex)
      complex*16 xsec(nex), xsec0(nex), xsecdx

      complex*16 conv1
      external conv1

      do 100 ie = 1,ne1
         xsec0(ie) = 0.0d0
         omega0 = omega(ie)
c        Add one more point to correct for the finite grid
c        at large energies. Use linear interpolation.
         dx = max( omega(ne1) - omega(ne1-1), 50*vicorr)
         xlast = omega(ne1)+dx
         dx = dx / ( omega(ne1) - omega(ne1-1))
         xsecdx = xsec(ne1)+ (xsec(ne1)-xsec(ne1-1)) * dx

c        first interval
         do 50  i = 1, ne1-1
            xsec0(ie) = xsec0(ie) +
     1      conv1(omega(i),omega(i+1),xsec(i),xsec(i+1),omega0,vicorr)
  50     continue
c        last interval
         xsec0(ie) = xsec0(ie) +
     1   conv1(omega(ne1),xlast,xsec(ne1),xsecdx,omega0,vicorr)
         xsec0(ie) = xsec0(ie) /real(pi)
  100 continue
      do 200 ie = 1, ne1
  200 xsec(ie) = xsec0(ie)

      return
      end

      complex*16 function conv1(x1,x2,y1,y2,x0,xloss)
c     convolution of function 1/(omega-omega0-i*xloss)/pi
c     makes linear interpolation for function between x1,x2 and
c     takes advantage that the integral can be taken analytically.
      implicit double precision (a-h, o-z)
      complex*16  y1, y2, t, coni,dum, a, b
      parameter (coni = (0.0,1.0))

      d = (x2-x1) / 2.0
      a = dble(y2-y1) / 2.0
      b = dble(y2+y1) / 2.0
      t = d / ( (x1+x2)/2 - x0 - coni*xloss )
      if (abs(t) .ge. 0.1) then
         dum = 2.0*a + (b - a/t) * log((1+t)/(1-t))
      else
         dum = 2.0*b*(t+t**3 / 3.0) - 2.0/3.0 * a*t**2
      endif
      conv1 = dimag (dum)

      d = (x2-x1) / 2.0
      a = dimag(y2-y1) / 2.0
      b = dimag(y2+y1) / 2.0
      t = d / ( (x1+x2)/2 - x0 - coni*xloss )
      if (abs(t) .ge. 0.1) then
         dum = 2.0*a + (b - a/t) * log((1+t)/(1-t))
      else
         dum = 2.0*b*(t+t**3 / 3.0) - 2.0/3.0 * a*t**2
      endif
      conv1 = conv1 + coni* dimag( dum)

      return
      end
      subroutine cpl0 (x, pl0, lmaxp1)
      implicit double precision (a-h, o-z)

c-----------------------------------------------------------------------
c
c     cpl0:  Calculate associated legendre polynomials p_l0(x)
c            by recursion.
c            Adapted from aslgndr.
c
c     first written: (25 june 86) by j. j. rehr
c
c     version 1 (25 june 86) (aslgndr)
c     version 2 (March, 1992) siz
c
c-----------------------------------------------------------------------

      dimension pl0 (lmaxp1)

      lmax = lmaxp1-1

c     calculate legendre polynomials p_l0(x) up to l=lmax
      pl0(1) = 1
      pl0(2) = x
      do 10  il = 2, lmax
         l = il-1
         pl0(il+1) = ( (2*l+1)*x*pl0(il) - l*pl0(l) ) / il
   10 continue

      return
      end
      subroutine csomm (dr,dp,dq,dpas,da,m,np)
c Modified to use complex p and q.  SIZ 4/91
c integration by the method of simpson of (dp+dq)*dr**m from 
c 0 to r=dr(np)
c dpas=exponential step;
c for r in the neighborhood of zero (dp+dq)=cte*r**da
c **********************************************************************
      implicit double precision (a-h,o-z)
      dimension dr(*)
      complex*16  dp(*),dq(*),da,dc
      mm=m+1
      d1=da+mm
      da=0.0
      db=0.0
      do 70 i=1,np
      dl=dr(i)**mm
      if (i.eq.1.or.i.eq.np) go to 10
      dl=dl+dl
      if ((i-2*(i/2)).eq.0) dl=dl+dl
   10 dc=dp(i)*dl
      da=da+dc
      dc=dq(i)*dl
      da=da+dc
   70 continue
      da=dpas*da/3
      dd=exp(dpas)-1.0
      db=d1*(d1+1.0)*dd*exp((d1-1.0)*dpas)
      db=dr(1)*(dr(2)**m)/db
      dd=(dr(1)**mm)*(1.0+1.0/(dd*(d1+1.0)))/d1
      da=da+dd*(dp(1)+dq(1))-db*(dp(2)+dq(2))
      return
      end
      subroutine csomm2 (dr,dp,dpas,da,rnrm,np)
c Modified to use complex p and q.  SIZ 4/91
c Modified to use double simpson integration ALA 3/97
c integration by the method of simpson of dp*dr from 
c 0 to r=rnrm  with proper end corrections
c dpas=exponential step;
c for r in the neighborhood of zero dp=cte*r**da
c **********************************************************************
      implicit double precision (a-h,o-z)
      dimension dr(*)
      complex*16  dp(*),da,dc

      d1=dble(da)+1
      da=0.0
      db=0.0
c      np-2=inrm -point of grid just below rnrm
      a1=log(rnrm/dr(np-2)) / dpas
      a2=a1**2/8.0d0
      a3=a1**3/12.0d0
      do 70 i=1,np
         if (i.eq.1) then
            dc=dp(i) *dr(i)*9.0d0/24.0d0
         elseif (i.eq.2) then
            dc=dp(i) *dr(i)*28.0d0/24.0d0
         elseif (i.eq.3) then
            dc=dp(i)*dr(i)*23.0d0/24.0d0
         elseif (i.eq.np-3) then
            dc=dp(i)*dr(i)*(25.0d0/24.0d0-a2+a3)
         elseif (i.eq.np-2) then
            dc=dp(i)*dr(i)*(0.5d0+a1-3*a2-a3)
         elseif (i.eq.np-1) then
            dc=dp(i)*dr(i)*(-1.0d0/24.0d0+5*a2-a3)
         elseif (i.eq.np) then
            dc=dp(i)*dr(i)*(-a2+a3)
         else
c           like trapesoidal rule
            dc=dp(i)*dr(i)
         endif
         da=da+dc
   70 continue
      da=dpas*da

c     add initial point (r=0) correction
      dd=exp(dpas)-1.0
      db=d1*(d1+1.0)*dd*exp((d1-1.0)*dpas)
      db=dr(1)/db
      dd=(dr(1))*(1.0+1.0/(dd*(d1+1.0)))/d1
      da=da+dd*dp(1)-db*dp(2)
      return
      end
      double precision function cwig3j (j1,j2,j3,m1,m2,ient)
c     wigner 3j coefficient for integers  (ient=1)
c                         or semiintegers (ient=2)
c     other arguments should be multiplied by ient
 
      implicit double precision (a-h,o-z)
      parameter (idim = 58)
      character*512 slog
c     dimensions  modified for larger arguments by ala 12.12.94
      dimension al(idim+1),m(12)
      save ini, al
      data ini/1/
c     idim-1 is the largest argument of factorial to calculate

      m3=-m1-m2
      if (ini) 1,21,1
c        initialisation of the log's of the factorials
 1    ini=0
      al(1)=0.0d 00
      do 11 i=1,idim
         b=i
 11      al(i+1)=al(i)+ log(b)
 21   cwig3j=0.0d 00
      if (((ient-1)*(ient-2)).ne.0) go to 101
      ii=ient+ient
c        test triangular inequalities, parity and maximum values of m
      if (( abs(m1)+ abs(m2)).eq.0.and.mod(j1+j2+j3,ii).ne.0) go to 99
      m(1)=j1+j2-j3
      m(2)=j2+j3-j1
      m(3)=j3+j1-j2
      m(4)=j1+m1
      m(5)=j1-m1
      m(6)=j2+m2
      m(7)=j2-m2
      m(8)=j3+m3
      m(9)=j3-m3
      m(10)=j1+j2+j3+ient
      m(11)=j2-j3-m1
      m(12)=j1-j3+m2
      do 41 i=1,12
         if (i.gt.10) go to 31
         if (m(i).lt.0) go to 99
 31      if (mod(m(i),ient).ne.0) go to 101
         m(i)=m(i)/ient
         if (m(i).gt.idim) go to 101
 41   continue

c        calculate 3j coefficient
      max0= max(m(11),m(12),0)+1
      min0= min(m(1),m(5),m(6))+1
      isig=1
      if (mod(max0-1,2).ne.0) isig=-isig
      c=-al(m(10)+1)
      do 61 i=1,9
 61   c=c+al(m(i)+1)
      c=c/2.0d 00
      do 71 i=max0,min0
      j=2-i
      b=al(i)+al(j+m(1))+al(j+m(5))+al(j+m(6))+al(i-m(11))+al(i-m(12))
      cwig3j=cwig3j+isig* exp(c-b)
 71   isig=-isig
      if (mod(j1-j2-m3,ii).ne.0) cwig3j=-cwig3j
 99   return
 101     write(slog,'(a,6i5)') 'error in cwig3j ',j1,j2,j3,m1,m2,ient
         call wlog(slog)
      stop
      end
      double precision function determ(array,nord,nrows)
c
c     calculate determinate of a square matrix
c        (from bevington "data reduction and error analysis
c         for the physical sciences" pg 294)
c     array: matrix to be analyzed
c     nord: order of matrix
c     nrows:  first dimension of matrix in calling routine
c
      double precision array(nrows,nrows)
      determ = 1.
      do 150 k=1,nord
c
c
        if (array(k,k).ne.0) go to 130
        do 100 j=k,nord
          if (array(k,j).ne.0) go to 110
  100   continue
        determ = 0.
        go to 160
c
  110   do 120 i=k,nord
          saved = array(i,j)
          array(i,j) = array(i,k)
  120   array(i,k) = saved
        determ = -determ
c
  130   determ = determ*array(k,k)
        if (k.ge.nord) go to 150
        k1 = k+1
        do 140 i=k1,nord
          do 140 j=k1,nord
  140   array(i,j) = array(i,j)-array(i,k)*array(k,j)/array(k,k)
  150 continue
  160 return
c end double precision function determ
      end
      double precision function dist (r0, r1)
c     find distance between cartesian points r0 and r1
      implicit double precision (a-h, o-z)
      dimension r0(3), r1(3)
      dist = 0
      do 10  i = 1, 3
         dist = dist + (r0(i) - r1(i))**2
   10 continue
      dist = sqrt (dist)
      return
      end
      double precision function rotwig (beta, jj, m1, m2, ient)
c     uses Wigner formula (Messiah eq.C.72) to calculate rotation matrix
c     for integers  (ient=1)  or semiintegers (ient=2)
c     other arguments (except beta) should be multiplied by ient
 
      implicit double precision (a-h,o-z)
      parameter (idim = 58)
c     dimensions  modified for larger arguments by ala 12.12.94
      dimension al(idim+1),m(12)
      save ini, al
      data ini/1/
c     idim-1 is the largest argument of factorial to calculate

      if (((ient-1)*(ient-2)).ne.0) stop ' Illegal ient in rotwig.'

      if (ini.eq.1) then
c       initialisation of the log's of the factorials
        ini=0
        al(1)=0.0d 00
        do 11 i=1,idim
           b=i
 11        al(i+1)=al(i)+ log(b)
      endif
      rotwig = 0.d0

      if ( m1.ge.0 .and. abs(m1).ge.abs(m2)) then
         m1p = m1 
         m2p = m2
         betap = beta
         isign = 1
      elseif (m2.ge.0 .and. abs(m2).ge.abs(m1)) then
         m1p = m2
         m2p = m1
         betap = - beta
         isign = 1
      elseif (m1.le.0 .and. abs(m1).ge.abs(m2)) then
         m1p = - m1
         m2p = - m2
         betap = beta
         isign = (-1)**( (m1-m2)/ient ) 
      else
         m1p = - m2
         m2p = - m1
         betap = - beta
         isign = (-1)**( (m2-m1)/ient ) 
      endif

      temp = 0.d0
      zeta = cos ( betap / 2.d0 )
      eta  = sin ( betap / 2.d0 )
      do 100 it = m1p - m2p, jj - m2p, ient
        m(1) = 1 + (jj+m1p) / ient
        m(2) = 1 + (jj-m1p) / ient
        m(3) = 1 + (jj+m2p) / ient
        m(4) = 1 + (jj-m2p) / ient
        m(5) = 1 + (jj+m1p-it) / ient
        m(6) = 1 + (jj-m2p-it) / ient
        m(7) = 1 + it / ient
        m(8) = 1 + (m2p-m1p+it) / ient
        m(9)  = (2*jj+m1p-m2p-2*it) / ient 
        m(10) = (2*it-m1p+m2p) / ient 
        factor = 0.d0
        do 110 i = 1,4
  110     factor = factor + al(m(i))/2.d0 - al(m(i+4))
c       special cases to resolve 0.d0**0 problem (a.ankudinov, may 2001)
        if (m(10).eq.0 .and. m(9).eq.0) then
          temp = temp + (-1)**(it/ient)*exp(factor)
        elseif (m(10).eq.0) then
          temp = temp + (-1)**(it/ient)*zeta**m(9)*exp(factor)
        elseif (m(9).eq.0) then
          temp = temp + (-1)**(it/ient)*eta**m(10)*exp(factor)
        else
c         general expression
          temp = temp+ (-1)**(it/ient)*zeta**m(9)*eta**m(10)*exp(factor)
        endif
  100 continue

      rotwig = isign * temp
     
      return
      end
      subroutine phamp (rmt, pu, qu, ck, jl, nl, jlp, nlp, ikap,
     1                  ph, amp)
c     calculate phase shift at mt radius
c     needs to calculate atan of complex variable (coded below)
      implicit double precision (a-h, o-z)
c={../HEADERS/const.h
      parameter (pi = 3.14159 26535 89793 23846 26433d0)
      parameter (one = 1, zero = 0)
      parameter (third = one/3)
      parameter (raddeg = 180 / pi)
      complex*16 coni
      parameter (coni = (0,1))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919 158 292 677 512 811d0)

      parameter (bohr = 0.529 177 249d0, ryd  = 13.605 698d0)
      parameter (hart = 2 * ryd)
      parameter (alpinv = 137.035 989 56d0)
c     fine structure alpha
      parameter (alphfs = 1 / alpinv)
c= ../HEADERS/const.h}
c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}
      external besjn, atan2c

      complex*16 pu, qu, ck,  jl, nl, jlp, nlp, ph, amp
      complex*16 xkr, a, b, factor

c     initialize staff
      xkr = ck*rmt
      isign=1
      if (ikap.lt.0) isign = -1
      a = ck*alphfs
      factor = isign*a/(1+sqrt(1+a**2))

c     find a and b that pu = rmt*(a*jl+b*nl), qu=factor*rmt*(a*jlp+b*nlp)
      a = isign*ck*xkr* (pu*nlp - qu*nl/factor)
      b = isign*ck*xkr* (qu*jl/factor - pu*jlp)

c     pu =  amp * rmt * (jl*cos(ph) - nl*sin(ph))
c     qu =  amp * rmt * (jlp*cos(ph) - nlp*sin(ph)) * factor
c     tan(ph) = - b/a
      b = -b
      call atan2c ( a, b, amp, ph)

      return
      end
      subroutine atancc(temp, phx)
c     phx=atan(temp), for complex numbers
      implicit double precision (a-h, o-z)
      complex*16 temp, phx

      xx = dble (temp)
      yy = dimag(temp)
      if (xx .ne. 0)  then
         alph = (1 - xx**2 - yy**2)
         alph = sqrt(alph**2 + 4*xx**2) - alph
         alph = alph / (2 * xx)
         alph = atan (alph)
      else
         alph = 0
      endif
      beta = (xx**2 + (yy+1)**2) / (xx**2 + (yy-1)**2)
      beta = log(beta) / 4
      phx = dcmplx (alph, beta)

      return
      end

      subroutine atan2c(a, b, ampl, phx)
c     for complex a, b find complex ampl, phx such that:
c     a= ampl*cos(phx)  and  b= ampl*sin(phx)
c     phx=atan(b/a)
      implicit double precision (a-h, o-z)
      parameter (pi = 3.14159 26535 89793 23846 26433d0)
      complex*16 a, b, ampl, phx, temp

      aa = abs(a)
      bb = abs(b)
      if (aa+bb.eq. 0) then
         ampl=0.d0
         phx =0.d0
      elseif ( aa.gt.bb) then
         temp = b/a
         call atancc ( temp, phx)
         ampl = a / cos(phx)
      else
         temp = a/b
         call atancc ( temp, phx)
         phx = pi / 2 - phx
         ampl = b/sin(phx)
      endif

      if (dble(ampl).lt. 0.d0) then
         ampl = -ampl
         phx = phx + pi
      endif

      return
      end
      subroutine exjlnl (z, l, jl, nl)

c     purpose:  to calculate the spherical bessel functions jl and nl
c               for l = 0 to 6  using exact analytic expression
c
c     arguments:
c       z = argument of jl and nl
c       l = integer order of spherical bessel function
c       jl = jl bessel function (abramowitz conventions)
c       nl = nl bessel function (abramowitz yl conventions)
c            Note that this nl = abramowitz yl.
c
c       analytic expressions from abramowitz 10.1.11 and 10.1.12
c       recurrence relation to get analytic j4,n4  eqns 10.1.19-22 ala

      implicit double precision (a-h, o-z)

      complex*16 z, jl, nl

      complex*16 cosz, sinz

c     Exact formulae unstable for very small z, so use series
c     expansion there.  Limit of .3 chosen for 9 digit agreement.
      if (abs(z) .lt. 0.3)  then
         call bjnser (z, l, jl, nl, 0)
      else
c        use analytic formulae
         cosz = cos(z)
         sinz = sin(z)

         if (l .eq. 0)  then
            jl =  sinz / z
            nl = -cosz / z

         elseif (l .eq. 1)  then
            jl =  sinz/z**2 - cosz/z
            nl = -cosz/z**2 - sinz/z

         elseif (l .eq. 2)  then
            jl = ( 3/z**3 - 1/z)*sinz - 3*cosz/z**2
            nl = (-3/z**3 + 1/z)*cosz - 3*sinz/z**2

         elseif (l .eq. 3)  then
            jl = ( 15/z**4 - 6/z**2)*sinz + (-15/z**3 + 1/z)*cosz
            nl = (-15/z**4 + 6/z**2)*cosz + (-15/z**3 + 1/z)*sinz

         elseif (l .eq. 4)  then
            jl = ( 105/z**5 - 45/z**3 + 1/z )*sinz + 
     1                ( -105/z**4 + 10/z**2 )*cosz
            nl = (-105/z**5 + 45/z**3 - 1/z )*cosz + 
     1                ( -105/z**4 + 10/z**2 )*sinz

         elseif (l .eq. 5)  then
            jl = ( 945/z**6 - 420/z**4 + 15/z**2 )*sinz + 
     1              ( -945/z**5 + 105/z**3 - 1/z )*cosz
            nl = (-945/z**6 + 420/z**4 - 15/z**2 )*cosz + 
     1              ( -945/z**5 + 105/z**3 - 1/z )*sinz

         elseif (l .eq. 6)  then
            jl = ( 10395/z**7 - 4725/z**5 + 210/z**3 - 1/z )*sinz + 
     1              ( -10395/z**6 + 1155/z**4 - 21/z**2 )*cosz
            nl = (-10395/z**7 + 4725/z**5 - 210/z**3 + 1/z )*cosz + 
     1              ( -10395/z**6 + 1155/z**4 - 21/z**2 )*sinz

         else
            stop 'exjlnl, l out of range'
         endif
      endif

      return
      end
      subroutine polint( xa, ya, n, x, y, dy)
c     draws a polynimial P(x) of order (n-1) through n points.
c     returns y = P(x) and dy - estimate of the error
c     adapted  from numerical recipies in fortran by Press et al.

      implicit double precision (a-h,o-z)
      integer n, nmax
      parameter (nmax=4)
      dimension xa(nmax), ya(nmax), c(nmax), d (nmax)

      ns = 1
      dif = abs (x-xa(1))
      do 10 i=1,n
         dift = abs(x-xa(i))
         if (dift.lt.dif) then
            ns = i
            dif = dift
         endif
         c(i) = ya(i)
         d(i) = ya(i)
  10  continue
      y = ya(ns)
      ns = ns-1
      do 30 m=1,n-1
         do 20 i=1,n-m
            ho = xa(i)-x
            hp = xa(i+m)-x
            w = c(i+1) - d(i)
            den = ho-hp
            if (den.eq.0) pause 'failure in polint'
            den = w/den
            d(i) = hp*den
            c(i) = ho*den
  20     continue
         if (2*ns .lt. n-m) then
            dy = c(ns+1)
         else
            dy = d(ns)
            ns = ns-1
         endif
         y = y + dy
  30  continue

      return
      end
      function sdist (r0, r1)
c     find distance squared between cartesian points r0 and r1
c     single precision
      dimension r0(3), r1(3)
      sdist = 0
      do 10  i = 1, 3
         sdist = sdist + (r0(i) - r1(i))**2
   10 continue
      sdist = sqrt(sdist)
      return
      end
      subroutine somm (dr,dp,dq,dpas,da,m,np)
c
c integration by the method of simpson of (dp+dq)*dr**m from
c 0 to r=dr(np)
c dpas=exponential step;
c for r in the neighborhood of zero (dp+dq)=cte*r**da
c **********************************************************************
      implicit double precision (a-h,o-z)
      dimension dr(np), dp(np), dq(np)
      mm=m+1
      d1=da+mm
      da=0.0
      db=0.0
      do 70 i=1,np
      dl=dr(i)**mm
      if (i.eq.1.or.i.eq.np) go to 10
      dl=dl+dl
      if ((i-2*(i/2)).eq.0) dl=dl+dl
   10 dc=dp(i)*dl
      if (dc) 20,40,30
   20 db=db+dc
      go to 40
   30 da=da+dc
   40 dc=dq(i)*dl
      if (dc) 50,70,60
   50 db=db+dc
      go to 70
   60 da=da+dc
   70 continue
      da = dpas * (da + db) / 3.0
      dc=exp(dpas)-1.0
      db=d1*(d1+1.0)*dc*exp((d1-1.0)*dpas)
      db=dr(1)*(dr(2)**m)/db
      dc=(dr(1)**mm)*(1.0+1.0/(dc*(d1+1.0)))/d1
      da=da+dc*(dp(1)+dq(1))-db*(dp(2)+dq(2))
      return
      end
      subroutine somm2 (dr,dp,dpas,da,rnrm,m,np)
c Modified to use complex p and q.  SIZ 4/91
c Modified to use double simpson integration ALA 3/97
c integration by the method of simpson of dp*dr from 
c 0 to r=rnrm  with proper end corrections
c dpas=exponential step;
c for r in the neighborhood of zero dp=cte*r**da
c **********************************************************************
      implicit double precision (a-h,o-z)
      dimension dr(*)
      dimension  dp(*)

      mm = m + 1
      d1=dble(da)+mm
      da=0.0
      db=0.0
c      np-2=inrm -point of grid just below rnrm
      a1=log(rnrm/dr(np-2)) / dpas
      a2=a1**2/8.0d0
      a3=a1**3/12.0d0
      do 70 i=1,np
         if (i.eq.1) then
            dc=dp(i) *dr(i)**mm*9.0d0/24.0d0
         elseif (i.eq.2) then
            dc=dp(i) *dr(i)**mm*28.0d0/24.0d0
         elseif (i.eq.3) then
            dc=dp(i)*dr(i)**mm*23.0d0/24.0d0
         elseif (i.eq.np-3) then
            dc=dp(i)*dr(i)**mm*(25.0d0/24.0d0-a2+a3)
         elseif (i.eq.np-2) then
            dc=dp(i)*dr(i)**mm*(0.5d0+a1-3*a2-a3)
         elseif (i.eq.np-1) then
            dc=dp(i)*dr(i)**mm*(-1.0d0/24.0d0+5*a2-a3)
         elseif (i.eq.np) then
            dc=dp(i)*dr(i)**mm*(-a2+a3)
         else
c           like trapesoidal rule
            dc=dp(i)*dr(i)**mm
         endif
         da=da+dc
   70 continue
      da=dpas*da

c     add initial point (r=0) correction
      dd=exp(dpas)-1.0
      db=d1*(d1+1.0)*dd*exp((d1-1.0)*dpas)
      db=dr(1)*(dr(2)**m)/db
      dd=(dr(1)**mm)*(1.0+1.0/(dd*(d1+1.0)))/d1
      da=da+dd*dp(1)-db*dp(2)
      return
      end
      subroutine strap (x, y, n, sum)

c     Trapeziodal integration of y(x), result in sum
c     SINGLE PRECISION
c     modified by ala to handle cases for E<Efermi
c     sum only positive numbers

      dimension x(n), y(n)

      sum = y(1) * abs(x(2) - x(1))
      do 10  i = 2, n-1
         sum = sum + y(i) * abs(x(i+1) - x(i-1))
   10 continue
      sum = sum + y(n) * abs(x(n) - x(n-1))
      sum = sum/2

      return
      end
c     interpolation and extrapolation by m-th order polynomial
c     maximum m = 3. Change nmax if needed.
c     Input x and y arrays, returns y value y0 at requested x value x0.
c     Dies on error.

      subroutine terp (x, y, n, m, x0, y0)
      implicit double precision (a-h, o-z)

      dimension x(n), y(n)

c     Find out between which x points x0 lies
      i = locat (x0, n, x)
      k = min( max(i-m/2,1) , n-m )
      call polint( x(k), y(k), m+1, x0, y0, dy)

      return
      end

      function locat (x, n, xx)
      integer  u, m, n
      double precision x, xx(n)

c     Binary search for index of grid point immediately below x.
c     Array xx required to be monotonic increasing.
c     Returns
c     0            x <  xx(1)
c     1            x =  xx(1)
c     i            x =  xx(i)
c     n            x >= xx(n)

      locat = 0
      u = n+1

   10 if (u-locat .gt. 1)  then
         m = (u + locat) / 2
         if (x .lt. xx(m))  then
            u = m
         else
            locat = m
         endif
         goto 10
      endif

      return
      end


c     These routines, terp1 and locat1, are special versions to
c     be used with ff2chi, which uses some single and some double
c     precision.  They are the same as the routines in terp.f.

      subroutine terp1 (x, y, n, x0, y0)
      implicit double precision (a-h, o-z)

      real x(n), y(n)

c     Find out between which x points x0 lies
      i = locat1 (x0, n, x)
c     if i < 1, set i=1, if i > n-1, set i=n-1
      i = max (i, 1)
      i = min (i, n-1)

      if (x(i+1) - x(i) .eq. 0)  stop 'TERP-1'

      y0 = y(i) +  (x0 - x(i)) * (y(i+1) - y(i)) / (x(i+1) - x(i))

      return
      end

      function locat1 (x, n, xx)
      integer  u, m, n
      double precision x
      real xx(n)

c     Binary search for index of grid point immediately below x.
c     Array xx required to be monotonic increasing.
c     Returns
c     0            x <  xx(1)
c     1            x =  xx(1)
c     i            x =  xx(i)
c     n            x >= xx(n)

      locat1 = 0
      u = n+1

   10 if (u-locat1 .gt. 1)  then
         m = (u + locat1) / 2
         if (x .lt. xx(m))  then
            u = m
         else
            locat1 = m
         endif
         goto 10
      endif

      return
      end
c     interpolation and extrapolation by m-th order polynomial
c     maximum m = 3. Change nmax if needed.
c     Input x and y arrays, returns y value y0 at requested x value x0.
c     Dies on error.

      subroutine terpc (x, y, n, m, x0, y0)
      implicit double precision (a-h, o-z)

      complex*16 y, y0, dy
      dimension x(n), y(n)

c     Find out between which x points x0 lies
      i = locat (x0, n, x)
      k = min( max(i-m/2,1) , n-m )
      call polinc( x(k), y(k), m+1, x0, y0, dy)

      return
      end

      subroutine polinc( xa, ya, n, x, y, dy)
c     draws a polynimial P(x) of order (n-1) through n points.
c     returns y = P(x) and dy - estimate of the error
c     adapted  from numerical recipies in fortran by Press et al.

      implicit double precision (a-h,o-z)
      complex*16 ya,y,dy,c,d,w,den
      integer n, nmax
      parameter (nmax=4)
      dimension xa(nmax), ya(nmax), c(nmax), d (nmax)

      ns = 1
      dif = abs (x-xa(1))
      do 10 i=1,n
         dift = abs(x-xa(i))
         if (dift.lt.dif) then
            ns = i
            dif = dift
         endif
         c(i) = ya(i)
         d(i) = ya(i)
  10  continue
      y = ya(ns)
      ns = ns-1
      do 30 m=1,n-1
         do 20 i=1,n-m
            ho = xa(i)-x
            hp = xa(i+m)-x
            w = c(i+1) - d(i)
            den = ho-hp
            if (den.eq.0) stop 'failure in polint'
            den = w/den
            d(i) = hp*den
            c(i) = ho*den
  20     continue
         if (2*ns .lt. n-m) then
            dy = c(ns+1)
         else
            dy = d(ns)
            ns = ns-1
         endif
         y = y + dy
  30  continue

      return
      end
      subroutine trap (x, y, n, sum)
      implicit double precision (a-h, o-z)

c     Trapeziodal integration of y(x), result in sum

      dimension x(n), y(n)

      sum = y(1) * (x(2) - x(1))
      do 10  i = 2, n-1
         sum = sum + y(i) * (x(i+1) - x(i-1))
   10 continue
      sum = sum + y(n) * (x(n) - x(n-1))
      sum = sum/2

      return
      end
      SUBROUTINE CQdrtc(Coef,Sol,NSol)
c     Combutes the zeros of a quadratic polynomial
ccccccccccccccccccccccccccccccccccccccccccccccccc            
c     Input
c     Coef - array of coefficients
      COMPLEX*16 Coef(3)
ccccccccccccccccccccccccccccccccccccccccccccccccc
c     Output
c     Sol  - Array of solutions
c     NSol - # of solutions (only one if Coef(1) = 0 etc.)
c     NSol = -1 means a and b are zero
      COMPLEX*16 Sol(2)
      INTEGER NSol
ccccccccccccccccccccccccccccccccccccccccccccccccc
c     Local Variables
      COMPLEX*16 q, Sqrt
      DOUBLE PRECISION Sgn

      IF(Coef(1).eq.0.d0) THEN
         IF(Coef(2).eq.0.d0) THEN
            NSol = -1
            RETURN
         ELSE
            NSol = 1
            Sol(1) = -Coef(3)/Coef(2)
         END IF
      ELSE
         NSol = 2
         Root = Sqrt(Coef(2)**2-4.d0*Coef(1)*Coef(3))
         Sgn  = SIGN(DBLE(CONJG(Coef(2))*Root),1.d0)
         q    = -0.5d0*(Coef(2) + Sgn*Root)
         
         Sol(1) = q/Coef(1)
         Sol(2) = Coef(3)/q
      END IF

      RETURN
      END


      SUBROUTINE CCubic(Coef,Sol,NSol)
c     Combutes the zeros of a cubic polynomial
ccccccccccccccccccccccccccccccccccccccccccccccccc            
c     Input
c     Coef - array of coefficients
      COMPLEX*16 Coef(4)
ccccccccccccccccccccccccccccccccccccccccccccccccc
c     Output
c     Sol  - Array of solutions
c     NSol - # of solutions (only one if Coef(1) = 0 etc.)
c     NSol = -1 means a, b, and c are zero
      COMPLEX*16 Sol(4)
      INTEGER NSol
ccccccccccccccccccccccccccccccccccccccccccccccccc
c     Local Variables
      COMPLEX*16 P1, P2, Q, R, Coef2(3), a, b, c
      DOUBLE PRECISION Sgn, Theta
c     PARAMETERS
      COMPLEX*16 I
      PARAMETER(I = (0.d0, 1.d0))
      DOUBLE PRECISION Pi
      PARAMETER(Pi = 3.141592653589793238462643d0)

      IF(Coef(1).eq.0.d0) THEN
         Coef2(1) = Coef(2)
         Coef2(2) = Coef(3)
         Coef2(3) = Coef(4)         
         CALL CQdrtc(Coef2,Sol,NSol)
      ELSE
         a = Coef(2)/Coef(1)
         b = Coef(3)/Coef(1)
         c = Coef(4)/Coef(1)
         NSol = 3
         Q = (a**2 - 3.d0*b)/9.d0
         R = (2.d0*a**3 - 9.d0*a*b + 27.d0*c)/54.d0

         IF(((DIMAG(Q).eq.0.d0).and.(DIMAG(R).eq.0.d0)).and.
     &        (DIMAG(R**2).lt.DIMAG(Q**3))) THEN
            Theta = ACOS (DBLE(R/SQRT(Q**3)))
            Sol(1) = -2*SQRT(Q)*Cos(Theta/3.d0) - a/3.d0
            Sol(2) = -2*SQRT(Q)*Cos((Theta+2.d0*Pi)/3.d0) - a/3.d0
            Sol(3) = -2*SQRT(Q)*Cos((Theta-2.d0*Pi)/3.d0) - a/3.d0
         ELSE
            Sgn = SIGN(1.d0, DBLE(CONJG(R)*SQRT(R**2-Q**3)))
            P1 = -(R + Sgn*SQRT(R**2-Q**3))**(1.d0/3.d0)
            IF(P1.eq.0.d0) THEN
               P2 = 0.d0
            ELSE
               P2 = Q/P1
            END IF
            Sol(1) = (P1 + P2) - a/3.d0
            Sol(2) = -0.5d0*(P1 + P2) - a/3.d0 +
     &           I*SQRT(3.d0)/2.d0*(P1-P2)
            Sol(3) = -0.5d0*(P1 + P2) - a/3.d0 -
     &           I*SQRT(3.d0)/2.d0*(P1-P2)
         END IF
      END IF

      RETURN
      END
      
c///////////////////////////////////////////////////////////////////////
c Distribution:  COMMON 1.0
c Copyright (c) [2002] University of Washington
c 
c This software was prepared in part with US Government Funding under
c DOE contract DE-FG03-97ER45623.

c Redistribution and use of this Distribution in source and binary
c formats, with or without modification is permitted, provided the 
c following conditions are met:
c 
c Redistributions must retain the above notices and the following list
c of conditions and disclaimer;
c 
c Modified formats carry the marking
c     "Based on or developed using Distribution: COMMON 1.0
c      COMMON 1.0 Copyright (c) [2002] University of Washington"
c 
c Recipient acknowledges the right of the University of Washington to
c prepare uses of this Distribution and its modifications that may be
c substantially similar or functionally equivalent to
c Recipient-prepared modifications.
c
c Recipient and anyone obtaining access to the Distribution through
c recipient's actions accept all risk associated with possession and
c use of the Distribution.
c
c THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESS OR IMPLIED
c WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
c MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
c IN NO EVENT SHALL THE UNIVERSITY OF WASHINGTON OR CONTRIBUTORS TO THE
c DISTRIBUTION BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
c EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
c PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
c REVENUE; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
c LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
c NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
c SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
c///////////////////////////////////////////////////////////////////////
c License is applicable for routines below, until otherwise specified.
c
      subroutine chopen (ios, fname, mod)
c     Writes error msg and stops if error in ios flag from open
c     statement.  fname is filename, mod is module with failed open.
      character*(*) fname, mod
      character*512 slog

c     open successful
      if (ios .le. 0)  return

c     error opening file, tell user and die.
      i = istrln(fname)
      j = istrln(mod)
      write(slog,100)  fname(1:i), mod(1:j)
      call wlog(slog)

  100 format (' Error opening file, ', a, 
     2        ' in module ', a)

      call wlog(' Fatal error')
      call par_stop('CHOPEN')
      end
      subroutine fixdsp (dxorg, dxnew, dgc0, dpc0, dgcx, dpcx, jnew)

c     This fixes up the dirac spinor components (dgc and dpc) from ATOM
c     for the xsect code.

      implicit double precision (a-h, o-z)

c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}
c={../HEADERS/const.h
      parameter (pi = 3.14159 26535 89793 23846 26433d0)
      parameter (one = 1, zero = 0)
      parameter (third = one/3)
      parameter (raddeg = 180 / pi)
      complex*16 coni
      parameter (coni = (0,1))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919 158 292 677 512 811d0)

      parameter (bohr = 0.529 177 249d0, ryd  = 13.605 698d0)
      parameter (hart = 2 * ryd)
      parameter (alpinv = 137.035 989 56d0)
c     fine structure alpha
      parameter (alphfs = 1 / alpinv)
c= ../HEADERS/const.h}

      dimension dgc0(251), dpc0(251)
      dimension dgcx(nrptx), dpcx(nrptx)

      dimension xorg(nrptx), xnew(nrptx)

      parameter (xx00 = 8.8)

c     statement functions to do indexing.  delta is 'dx' for current
c     grid.  jjj is index of grid point immediately before 'r'
      xxx(j) = -xx00 + (j-1)*delta
      rrr(j) = exp (-xx00 + (j-1)*delta)
      jjj(r) = (log(r) + xx00) / delta + 1

c     Use linear interpolation in x whether necessary or not.  If
c     new grid is same as old, it shouldn't make any difference.

c     relation between x, r, and j.  xx00 = 8.8 for all grids
c     in this version, change it if more flexibility is necessary.
c     xx = -xx00 + (j-1)*delta
c     rr = exp (xx)
c     jj = (log(r) + xx00) / delta + 1; this is j immediately BELOW r

c     The dgc and dpc arrays are zero beyond a certain point, usually
c     inside the muffin tin radius.  Find this distance.
      do 100  i = 251, 1, -1
         if ( abs(dgc0(i)) .ge. 1.0d-11 .or. 
     1        abs(dpc0(i)) .ge. 1.0d-11 )  then
            imax = i
            goto 16
         endif
  100 continue
      call wlog(' Should never see this line from sub fixdsp')
   16 continue
c     jmax is the first point where both dpc and dgc are zero in
c     the original grid
      jmax = imax + 1
      if (jmax.gt.251) jmax = 251

      delta = dxorg
      do 10  j = 1, jmax
         xorg(j) = xxx(j)
   10 continue
      rmax = rrr(jmax)

c     How far out do we go in the new grid?  To the last new grid
c     point before jmax.  Everything will be zero beyond jmax.
      delta = dxnew
      jnew = jjj(rmax)
      do 20  j = 1, jnew
         xnew(j) = xxx(j)
   20 continue

c     interpolate to new grid using x, only inside of rmax
      do 30  j = 1, jnew
         call terp (xorg, dgc0,  jmax, 3, xnew(j), dgcx(j))
         call terp (xorg, dpc0,  jmax, 3, xnew(j), dpcx(j))
   30 continue

c     and zero the arrays past rmax
      do 32  j = jnew+1, nrptx
         dgcx(j) = 0
         dpcx(j) = 0
   32 continue

      return
      end
      subroutine fixdsx (iph, dxorg, dxnew, dgc, dpc, dgcn, dpcn)

c     This fixes up the dirac spinor components (dgc and dpc) from ATOM
c     for the xsect and phase codes.

      implicit double precision (a-h, o-z)

c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}
c={../HEADERS/const.h
      parameter (pi = 3.14159 26535 89793 23846 26433d0)
      parameter (one = 1, zero = 0)
      parameter (third = one/3)
      parameter (raddeg = 180 / pi)
      complex*16 coni
      parameter (coni = (0,1))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919 158 292 677 512 811d0)

      parameter (bohr = 0.529 177 249d0, ryd  = 13.605 698d0)
      parameter (hart = 2 * ryd)
      parameter (alpinv = 137.035 989 56d0)
c     fine structure alpha
      parameter (alphfs = 1 / alpinv)
c= ../HEADERS/const.h}

      dimension dgc(251,30,0:nphx+1), dpc(251,30,0:nphx+1)
      dimension dgcn(nrptx,30), dpcn(nrptx,30)

      dimension xorg(nrptx), xnew(nrptx)

      parameter (xx00 = 8.8)

c     statement functions to do indexing.  delta is 'dx' for current
c     grid.  jjj is index of grid point immediately before 'r'
      xxx(j) = -xx00 + (j-1)*delta
      rrr(j) = exp (-xx00 + (j-1)*delta)
      jjj(r) = (log(r) + xx00) / delta + 1

c     Use linear interpolation in x whether necessary or not.  If
c     new grid is same as old, it shouldn't make any difference.

c     relation between x, r, and j.  xx00 = 8.8 for all grids
c     in this version, change it if more flexibility is necessary.
c     xx = -xx00 + (j-1)*delta
c     rr = exp (xx)
c     jj = (log(r) + xx00) / delta + 1; this is j immediately BELOW r

c     The dgc and dpc arrays are zero beyond a certain point, usually
c     inside the muffin tin radius.  Find this distance.

      delta = dxorg
      do 10  j = 1, 251
         xorg(j) = xxx(j)
   10 continue

      delta = dxnew
      do 20  j = 1, nrptx
         xnew(j) = xxx(j)
   20 continue

      do 200 iorb = 1, 30
         imax = 0
         do 100  i = 251, 1, -1
            if ( abs(dgc(i,iorb,iph)) .ge. 1.0d-11 .or. 
     1           abs(dpc(i,iorb,iph)) .ge. 1.0d-11 )  then
               imax = i
               goto 16
            endif
  100    continue
   16    continue
         if (imax .eq. 0) then
            jnew = 0
            goto 35
         endif
c        jmax is the first point where both dpc and dgc are zero in
c        the original grid
         jmax = imax + 1
         if (jmax .gt. 251) jmax = 251

         delta = dxorg
         rmax = rrr(jmax)

c        How far out do we go in the new grid?  To the last new grid
c        point before jmax.  Everything will be zero beyond jmax.
         delta = dxnew
         jnew = jjj(rmax)

c        interpolate to new grid using x, only inside of rmax
         do 30  j = 1, jnew
            call terp(xorg,dgc(1,iorb,iph),jmax,3, xnew(j),dgcn(j,iorb))
            call terp(xorg,dpc(1,iorb,iph),jmax,3, xnew(j),dpcn(j,iorb))
   30    continue

c        and zero the arrays past rmax
   35    do 40  j = jnew+1, nrptx
            dgcn(j,iorb) = 0
            dpcn(j,iorb) = 0
   40    continue
  200 continue

      return
      end
      subroutine fixvar (rmt, edens, vtot, dmag,
     1                   vint, rhoint, dxorg, dxnew, jumprm,
     2                   vjump, ri, vtotph, rhoph, dmagx)

      implicit double precision (a-h, o-z)

c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}
c={../HEADERS/const.h
      parameter (pi = 3.14159 26535 89793 23846 26433d0)
      parameter (one = 1, zero = 0)
      parameter (third = one/3)
      parameter (raddeg = 180 / pi)
      complex*16 coni
      parameter (coni = (0,1))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919 158 292 677 512 811d0)

      parameter (bohr = 0.529 177 249d0, ryd  = 13.605 698d0)
      parameter (hart = 2 * ryd)
      parameter (alpinv = 137.035 989 56d0)
c     fine structure alpha
      parameter (alphfs = 1 / alpinv)
c= ../HEADERS/const.h}


      dimension edens(251), vtot (251), dmag(251)
      dimension vtotph(nrptx), rhoph(nrptx), dmagx(nrptx)
      dimension ri(nrptx)
      dimension xorg(nrptx), xnew(nrptx)

      parameter (xx00 = 8.8)

c     statement functions to do indexing.  delta is 'dx' for current
c     grid.  jjj is index of grid point immediately before 'r'
      xxx(j) = -xx00 + (j-1)*delta
      rrr(j) = exp (-xx00 + (j-1)*delta)
      jjj(r) = (log(r) + xx00) / delta + 1

c     PHASE needs
c     vtot = total potential including gs xcorr, no r**2
c     edens = rho, charge density, no factor of 4*pi, no r**2
c     From overlapping, vtot = potential only, ok as is
c                       edens = density*4*pi, so fix this here.
c     ri = r grid through imt+1

c     Only values inside the muffin tin are used, except that XCPOT
c     (in PHASE) uses values at imt+1 and requires these to be the
c     interstitial values.  So set the last part of the arrays to
c     interstitial values...

c     Use linear interpolation in x whether necessary or not.  If
c     new grid is same as old, it shouldn't make any difference.

c     relation between x, r, and j.  xx00 = 8.8 for all grids
c     in this version, change it if more flexibility is necessary.
      
c     xx = -xx00 + (j-1)*delta
c     rr = exp (xx)
c     jj = (log(r) + xx00) / delta + 1; this is j immediately BELOW r

      delta = dxorg
      jmtorg = jjj(rmt)
      jriorg = jmtorg + 1
      jrior1 = jriorg + 1
      do 10  j = 1, jrior1
         xorg(j) = xxx(j)
   10 continue

      delta = dxnew
      jmtnew = jjj(rmt)
      jrinew = jmtnew + 1
      jrine1 = jrinew + 1
      do 20  j = 1, jrine1
         xnew(j) = xxx(j)
   20 continue

c     interpolate to new grid using x, only inside of muffintin
c     jri (first interstitial point) must be set to interstitial value
      do 30  j = 1, jrinew
         call terp (xorg, vtot,  jriorg, 3, xnew(j), vtotph(j))
         call terp (xorg, edens, jrior1, 3, xnew(j), rhoph(j))
         call terp (xorg, dmag,  jrior1, 3, xnew(j), dmagx(j))
   30 continue

      if (jumprm .eq. 1) then
         xmt = log(rmt)
         call terp (xorg, vtot,  jriorg, 3, xmt, vmt)
         vjump = vint - vmt
      endif
      if (jumprm .gt. 0) then
         do 90  j = 1, jrinew
            vtotph(j) = vtotph(j) + vjump
   90    continue
      endif

      delta = dxnew
      do 180  j = 1, nrptx
         ri(j) = rrr(j)
  180 continue
      do 190  j = 1, jrinew
         rhoph(j) = rhoph(j)/(4*pi)
  190 continue
      do 200  j = jrinew+1, nrptx
         vtotph(j) = vint
         rhoph(j) = rhoint/(4*pi)
c fix later : need to calculate interstitial dmint
c      want interpolation beyond mt also
         dmagx(j) = 0.0d0
  200 continue

      return
      end
      subroutine getorb (iz, ihole, xion, iunf, norb, norbco, iorb,
     1                  iholep, nqn, nk, xnel, xnval, xmag)
c     Gets orbital data for chosen element.  Input is:
c       iz - atomic number of desired element,
c       ihole - index of core-hole orbital
c       xion  - ionicity (usually zero)
c     other arguments are output.
c       norb - total number of orbitals
c       norbco - number of core orbitals
c       iorb - index of orbital for making projections (last occupied)
c       iholep - index of core hole orbital in compacted list
c       nqn - principal quantum number for each orbital
c       nk - quantum number kappa for each orbital
c       xnel - occupation for each orbital
c       xnval - valence occupation for each orbital
c       xmag - spin magnetization for each orbital
c     Feel free to change occupation numbers for element of interest.
c     ival(i) is necessary only for partly nonlocal exchange model.
c     iocc(i) and ival(i) can be fractional
c     But you have to keep the sum of iocc(i) equal to nuclear charge.
c     Also ival(i) should be equal to iocc(i) or zero. 
c     Otherwise you have to change this subroutine or contact authors 
c     for help.

      implicit double precision (a-h, o-z)

c     Written by Steven Zabinsky, July 1989
c     modified (20 aug 1989)  table increased to at no 99
c     Recipe for final state configuration is changed. Valence
c     electron occupations are added. ala 17.1.1996

c     Table for each element has occupation of the various levels.
c     The order of the levels in each array is:

c     element  level     principal qn (nqn), kappa qn (nk)
c           1  1s        1  -1
c           2  2s        2  -1
c           3  2p1/2     2   1
c           4  2p3/2     2  -2
c           5  3s        3  -1
c           6  3p1/2     3   1
c           7  3p3/2     3  -2
c           8  3d3/2     3   2
c           9  3d5/2     3  -3
c          10  4s        4  -1
c          11  4p1/2     4   1
c          12  4p3/2     4  -2
c          13  4d3/2     4   2
c          14  4d5/2     4  -3
c          15  4f5/2     4   3
c          16  4f7/2     4  -4
c          17  5s        5  -1
c          18  5p1/2     5   1
c          19  5p3/2     5  -2
c          20  5d3/2     5   2
c          21  5d5/2     5  -3
c          22  5f5/2     5   3
c          23  5f7/2     5  -4
c          24  6s        6  -1
c          25  6p1/2     6   1
c          26  6p3/2     6  -2
c          27  6d3/2     6   2
c          28  6d5/2     6  -3
c          29  7s        7  -1

      dimension nqn(30), nk(30), xnel(30), xnval(30), xmag(30)
      dimension kappa (29)
      real iocc, ival, ispn
      dimension iocc (100, 29), ival (100, 29), ispn (100, 29)
      dimension nnum (29), iorb(-4:3)
      character*512 slog

c     kappa quantum number for each orbital
c     k = - (j + 1/2)  if l = j - 1/2
c     k = + (j + 1/2)  if l = j + 1/2
      data kappa /-1,-1, 1,-2,-1,   1,-2, 2,-3,-1,   1,-2, 2,-3, 3,
     1            -4,-1, 1,-2, 2,  -3, 3,-4,-1, 1,  -2, 2,-3,-1/

c     principal quantum number (energy eigenvalue)
      data nnum  /1,2,2,2,3,  3,3,3,3,4,  4,4,4,4,4,
     1            4,5,5,5,5,  5,5,5,6,6,  6,6,6,7/

c     occupation of each level for z = 1, 99
      data (iocc( 1,i),i=1,29)  /1,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival( 1,i),i=1,29)  /1,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn( 1,i),i=1,29)  /1,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc( 2,i),i=1,29)  /2,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival( 2,i),i=1,29)  /2,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn( 2,i),i=1,29)  /1,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc( 3,i),i=1,29)  /2,1,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival( 3,i),i=1,29)  /0,1,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn( 3,i),i=1,29)  /0,1,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc( 4,i),i=1,29)  /2,2,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival( 4,i),i=1,29)  /0,2,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn( 4,i),i=1,29)  /0,1,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc( 5,i),i=1,29)  /2,2,1,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival( 5,i),i=1,29)  /0,2,1,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn( 5,i),i=1,29)  /0,0,1,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

c     data (iocc( 6,i),i=1,29)  /2,2,2,0,0,  0,0,0,0,0,  0,0,0,0,0,
      data (iocc( 6,i),i=1,29)  /2,1,2,1,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
c     data (ival( 6,i),i=1,29)  /0,2,2,0,0,  0,0,0,0,0,  0,0,0,0,0,
      data (ival( 6,i),i=1,29)  /0,1,2,1,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn( 6,i),i=1,29)  /0,0,1,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc( 7,i),i=1,29)  /2,2,2,1,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival( 7,i),i=1,29)  /0,2,2,1,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn( 7,i),i=1,29)  /0,0,0,1,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc( 8,i),i=1,29)  /2,2,2,2,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival( 8,i),i=1,29)  /0,2,2,2,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn( 8,i),i=1,29)  /0,0,0,1,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc( 9,i),i=1,29)  /2,2,2,3,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival( 9,i),i=1,29)  /0,2,2,3,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn( 9,i),i=1,29)  /0,0,0,1,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(10,i),i=1,29)  /2,2,2,4,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(10,i),i=1,29)  /0,0,2,4,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(10,i),i=1,29)  /0,0,0,1,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(11,i),i=1,29)  /2,2,2,4,1,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(11,i),i=1,29)  /0,0,2,4,1,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(11,i),i=1,29)  /0,0,0,0,1,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(12,i),i=1,29)  /2,2,2,4,1,  1,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(12,i),i=1,29)  /0,0,0,0,1,  1,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(12,i),i=1,29)  /0,0,0,0,1,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(13,i),i=1,29)  /2,2,2,4,2,  1,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(13,i),i=1,29)  /0,0,0,0,2,  1,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(13,i),i=1,29)  /0,0,0,0,0,  1,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(14,i),i=1,29)  /2,2,2,4,2,  2,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(14,i),i=1,29)  /0,0,0,0,2,  2,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(14,i),i=1,29)  /0,0,0,0,0,  1,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(15,i),i=1,29)  /2,2,2,4,2,  2,1,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(15,i),i=1,29)  /0,0,0,0,2,  2,1,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(15,i),i=1,29)  /0,0,0,0,0,  0,1,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(16,i),i=1,29)  /2,2,2,4,2,  2,2,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(16,i),i=1,29)  /0,0,0,0,2,  2,2,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(16,i),i=1,29)  /0,0,0,0,0,  0,1,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(17,i),i=1,29)  /2,2,2,4,2,  2,3,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(17,i),i=1,29)  /0,0,0,0,2,  2,3,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(17,i),i=1,29)  /0,0,0,0,0,  0,1,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(18,i),i=1,29)  /2,2,2,4,2,  2,4,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(18,i),i=1,29)  /0,0,0,0,2,  2,4,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(18,i),i=1,29)  /0,0,0,0,0,  0,1,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(19,i),i=1,29)  /2,2,2,4,2,  2,4,0,0,1,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(19,i),i=1,29)  /0,0,0,0,2,  2,4,0,0,1,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(19,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,1,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(20,i),i=1,29)  /2,2,2,4,2,  2,4,0,0,1,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(20,i),i=1,29)  /0,0,0,0,0,  2,4,0,0,1,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(20,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,1,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(21,i),i=1,29)  /2,2,2,4,2,  2,4,1,0,1,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(21,i),i=1,29)  /0,0,0,0,0,  2,4,1,0,1,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(21,i),i=1,29)  /0,0,0,0,0,  0,0,1,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(22,i),i=1,29)  /2,2,2,4,2,  2,4,2,0,1,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(22,i),i=1,29)  /0,0,0,0,0,  2,4,2,0,1,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(22,i),i=1,29)  /0,0,0,0,0,  0,0,2,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(23,i),i=1,29)  /2,2,2,4,2,  2,4,3,0,1,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(23,i),i=1,29)  /0,0,0,0,0,  2,4,3,0,1,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(23,i),i=1,29)  /0,0,0,0,0,  0,0,3,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(24,i),i=1,29)  /2,2,2,4,2,  2,4,4,0,1,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(24,i),i=1,29)  /0,0,0,0,0,  2,4,4,0,1,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(24,i),i=1,29)  /0,0,0,0,0,  0,0,4,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(25,i),i=1,29)  /2,2,2,4,2,  2,4,4,1,1,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(25,i),i=1,29)  /0,0,0,0,0,  0,0,4,1,1,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(25,i),i=1,29)  /0,0,0,0,0,  0,0,4,1,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(26,i),i=1,29)  /2,2,2,4,2,  2,4,4,2,1,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(26,i),i=1,29)  /0,0,0,0,0,  0,0,4,2,1,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(26,i),i=1,29)  /0,0,0,0,0,  0,0,2,2,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(27,i),i=1,29)  /2,2,2,4,2,  2,4,4,3,1,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(27,i),i=1,29)  /0,0,0,0,0,  0,0,4,3,1,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(27,i),i=1,29)  /0,0,0,0,0,  0,0,0,3,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(28,i),i=1,29)  /2,2,2,4,2,  2,4,4,4,1,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(28,i),i=1,29)  /0,0,0,0,0,  0,0,4,4,1,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(28,i),i=1,29)  /0,0,0,0,0,  0,0,0,1,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(29,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,1,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(29,i),i=1,29)  /0,0,0,0,0,  0,0,4,6,1,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(29,i),i=1,29)  /0,0,0,0,0,  0,0,0,1,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(30,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,1,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(30,i),i=1,29)  /0,0,0,0,0,  0,0,4,6,1,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(30,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,1,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(31,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(31,i),i=1,29)  /0,0,0,0,0,  0,0,4,6,2,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(31,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(32,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(32,i),i=1,29)  /0,0,0,0,0,  0,0,4,6,2,  2,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(32,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(33,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,1,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(33,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,2,  2,1,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(33,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,1,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(34,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,2,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(34,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,2,  2,2,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(34,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,1,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(35,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,3,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(35,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,2,  2,3,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(35,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,1,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(36,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(36,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,2,  2,4,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(36,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,1,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(37,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,0,0,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(37,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,2,  2,4,0,0,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(37,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(38,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,0,0,0,
     1                           0,2,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(38,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  2,4,0,0,0,
     1                           0,2,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(38,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(39,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,1,0,0,
     1                           0,2,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(39,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  2,4,1,0,0,
     1                           0,2,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(39,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,1,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(40,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,2,0,0,
     1                           0,2,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(40,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  2,4,2,0,0,
     1                           0,2,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(40,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,2,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(41,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,0,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(41,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  2,4,4,0,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(41,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,3,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(42,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,1,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(42,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  2,4,4,1,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(42,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,4,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(43,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,1,0,
     1                           0,2,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(43,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,4,1,0,
     1                           0,2,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(43,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,4,1,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(44,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,3,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(44,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,4,3,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(44,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,2,2,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(45,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,4,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(45,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,4,4,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(45,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,2,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(46,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(46,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,4,6,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(46,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,1,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(47,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(47,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,4,6,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(47,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(48,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,2,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(48,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,4,6,0,
     1                           0,2,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(48,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(49,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,2,1,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(49,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,4,6,0,
     1                           0,2,1,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(49,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,1,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(50,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,2,2,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(50,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,4,6,0,
     1                           0,2,2,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(50,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,1,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(51,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,2,2,1,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(51,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,4,6,0,
     1                           0,2,2,1,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(51,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,1,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(52,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,2,2,2,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(52,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,2,2,2,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(52,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,1,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(53,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,2,2,3,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(53,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,2,2,3,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(53,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,1,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(54,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,2,2,4,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(54,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,2,2,4,0,  0,0,0,0,0,  0,0,0,0/
      data (ispn(54,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,1,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(55,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,2,2,4,0,  0,0,0,1,0,  0,0,0,0/
      data (ival(55,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,2,2,4,0,  0,0,0,1,0,  0,0,0,0/
      data (ispn(55,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,1,0,  0,0,0,0/

      data (iocc(56,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,2,2,4,0,  0,0,0,2,0,  0,0,0,0/
      data (ival(56,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,2,2,4,0,  0,0,0,2,0,  0,0,0,0/
      data (ispn(56,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,1,0,  0,0,0,0/

      data (iocc(57,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,2,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ival(57,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,2,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ispn(57,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,1,  0,0,0,0,0,  0,0,0,0/

      data (iocc(58,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,1,
     1                           0,2,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ival(58,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,1,
     1                           0,0,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ispn(58,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,1,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(59,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,2,
     1                           0,2,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ival(59,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,2,
     1                           0,0,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ispn(59,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,2,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(60,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,3,
     1                           0,2,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ival(60,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,3,
     1                           0,0,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ispn(60,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,3,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(61,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,4,
     1                           0,2,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ival(61,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,4,
     1                           0,0,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ispn(61,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,4,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(62,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,5,
     1                           0,2,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ival(62,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,5,
     1                           0,0,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ispn(62,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,5,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(63,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           0,2,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ival(63,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,6,
     1                           0,0,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ispn(63,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,6,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(64,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           1,2,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ival(64,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,6,
     1                           1,0,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ispn(64,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,6,
     1                           1,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(65,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           2,2,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ival(65,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,6,
     1                           2,0,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ispn(65,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,4,
     1                           2,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(66,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           3,2,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ival(66,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,6,
     1                           3,0,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ispn(66,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,2,
     1                           3,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(67,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           4,2,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ival(67,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,6,
     1                           4,0,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ispn(67,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           4,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(68,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           5,2,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ival(68,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,6,
     1                           5,0,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ispn(68,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           3,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(69,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           6,2,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ival(69,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,6,
     1                           6,0,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ispn(69,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           2,0,0,0,0,  0,0,0,2,0,  0,0,0,0/

      data (iocc(70,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           7,2,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ival(70,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,6,
     1                           7,0,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ispn(70,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           1,0,0,0,0,  0,0,0,0,0,  0,0,0,0/

      data (iocc(71,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ival(71,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,6,
     1                           8,0,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ispn(71,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,1,  0,0,0,0,0,  0,0,0,0/

      data (iocc(72,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,2,  0,0,0,2,0,  0,0,0,0/
      data (ival(72,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,6,
     1                           8,0,2,4,2,  0,0,0,2,0,  0,0,0,0/
      data (ispn(72,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,2,  0,0,0,0,0,  0,0,0,0/

      data (iocc(73,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,3,  0,0,0,2,0,  0,0,0,0/
      data (ival(73,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,6,
     1                           8,0,0,0,3,  0,0,0,2,0,  0,0,0,0/
      data (ispn(73,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,3,  0,0,0,0,0,  0,0,0,0/

      data (iocc(74,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,3,  1,0,0,2,0,  0,0,0,0/
c    1                           8,2,2,4,4,  0,0,0,2,0,  0,0,0,0/
      data (ival(74,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,6,
     1                           8,0,0,0,3,  1,0,0,2,0,  0,0,0,0/
c    1                           8,0,0,0,4,  0,0,0,2,0,  0,0,0,0/
      data (ispn(74,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,4,  0,0,0,0,0,  0,0,0,0/

      data (iocc(75,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  1,0,0,2,0,  0,0,0,0/
      data (ival(75,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,4,  1,0,0,2,0,  0,0,0,0/
      data (ispn(75,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,4,  1,0,0,0,0,  0,0,0,0/

      data (iocc(76,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  2,0,0,2,0,  0,0,0,0/
      data (ival(76,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,4,  2,0,0,2,0,  0,0,0,0/
      data (ispn(76,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,2,  2,0,0,0,0,  0,0,0,0/

      data (iocc(77,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  3,0,0,2,0,  0,0,0,0/
      data (ival(77,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,4,  3,0,0,2,0,  0,0,0,0/
      data (ispn(77,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  3,0,0,0,0,  0,0,0,0/

      data (iocc(78,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  5,0,0,1,0,  0,0,0,0/
      data (ival(78,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,4,  5,0,0,1,0,  0,0,0,0/
      data (ispn(78,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  2,0,0,0,0,  0,0,0,0/

      data (iocc(79,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,1,0,  0,0,0,0/
      data (ival(79,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,4,  6,0,0,1,0,  0,0,0,0/
      data (ispn(79,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  1,0,0,0,0,  0,0,0,0/

      data (iocc(80,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,0,  0,0,0,0/
      data (ival(80,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,4,  6,0,0,2,0,  0,0,0,0/
      data (ispn(80,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,1,0,  0,0,0,0/

      data (iocc(81,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,1,  0,0,0,0/
      data (ival(81,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,4,  6,0,0,2,1,  0,0,0,0/
      data (ispn(81,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,1,  0,0,0,0/

      data (iocc(82,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,2,  0,0,0,0/
      data (ival(82,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,4,  6,0,0,2,2,  0,0,0,0/
      data (ispn(82,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,1,  0,0,0,0/

      data (iocc(83,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,2,  1,0,0,0/
      data (ival(83,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,4,  6,0,0,2,2,  1,0,0,0/
      data (ispn(83,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  1,0,0,0/

      data (iocc(84,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,2,  2,0,0,0/
      data (ival(84,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,4,  6,0,0,2,2,  2,0,0,0/
      data (ispn(84,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  1,0,0,0/

      data (iocc(85,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,2,  3,0,0,0/
      data (ival(85,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,2,2,  3,0,0,0/
      data (ispn(85,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  1,0,0,0/

      data (iocc(86,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,2,  4,0,0,0/
      data (ival(86,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,2,2,  4,0,0,0/
      data (ispn(86,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  1,0,0,0/

      data (iocc(87,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,2,  4,0,0,1/
      data (ival(87,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,2,2,  4,0,0,1/
      data (ispn(87,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,1/

      data (iocc(88,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,2,  4,0,0,2/
      data (ival(88,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,2,  4,0,0,2/
      data (ispn(88,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,1/

      data (iocc(89,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,2,  4,1,0,2/
      data (ival(89,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,2,  4,1,0,2/
      data (ispn(89,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,1,0,0/

      data (iocc(90,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,2,  4,2,0,2/
      data (ival(90,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,2,  4,2,0,2/
      data (ispn(90,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,2,0,0/

      data (iocc(91,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,2,0,2,2,  4,1,0,2/
      data (ival(91,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,2,0,0,2,  4,1,0,2/
      data (ispn(91,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,2,0,0,0,  0,0,0,0/

      data (iocc(92,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,3,0,2,2,  4,1,0,2/
      data (ival(92,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,3,0,0,2,  4,1,0,2/
      data (ispn(92,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,1.5,0,0,0,  0,0,0,0/

      data (iocc(93,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,4,0,2,2,  4,1,0,2/
      data (ival(93,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,4,0,0,2,  4,1,0,2/
      data (ispn(93,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,4,0,0,0,  0,0,0,0/

      data (iocc(94,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,6,0,2,2,  4,0,0,2/
      data (ival(94,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,6,0,0,2,  4,0,0,2/
      data (ispn(94,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,5,0,0,0,  0,0,0,0/

      data (iocc(95,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,6,1,2,2,  4,0,0,2/
      data (ival(95,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,6,1,0,2,  4,0,0,2/
      data (ispn(95,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,5,1,0,0,  0,0,0,0/

      data (iocc(96,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,6,2,2,2,  4,0,0,2/
      data (ival(96,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,6,2,0,2,  4,0,0,2/
      data (ispn(96,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,5,2,0,0,  0,0,0,0/

      data (iocc(97,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,6,3,2,2,  4,0,0,2/
      data (ival(97,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,6,3,0,2,  4,0,0,2/
      data (ispn(97,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,3,3,0,0,  0,0,0,0/

      data (iocc(98,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,6,4,2,2,  4,0,0,2/
      data (ival(98,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,6,4,0,2,  4,0,0,2/
      data (ispn(98,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,1,4,0,0,  0,0,0,0/

      data (iocc(99,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,6,5,2,2,  4,0,0,2/
      data (ival(99,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,6,5,0,2,  4,0,0,2/
      data (ispn(99,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,4,0,0,  0,0,0,0/

      data (iocc(100,i),i=1,29) /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,6,6,2,2,  4,0,0,2/
      data (ival(100,i),i=1,29) /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,6,6,0,2,  4,0,0,2/
      data (ispn(100,i),i=1,29) /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,3,0,0,  0,0,0,0/

      if (iz .lt. 1  .or.  iz .gt. 99)  then
    8    format(' Atomic number ', i5, ' not available.')
         write(slog,8)  iz
         call wlog(slog)
         call par_stop('GETORB-0')
      endif

      ion = nint(xion)
      delion=xion-ion
      index = iz - ion
      ilast = 0
      iscr = 0
      iion = 0
      iholep = ihole

c     find last occupied orbital (ilast) and iion for delion.ge.0
      do 30 i=29,1,-1
         if (iion.eq.0 .and. iocc(index,i).gt.delion) iion=i
         if (ilast.eq.0 .and. iocc(index,i).gt.0) ilast=i
 30   continue

      if (ihole.gt.0) then
         if ( iocc(index,ihole) .lt. 1 ) then
           call wlog(' Cannot remove an electron from this level')
           call par_stop('GETORB-1')
         endif
      endif
      if (ihole.eq.ilast) then 
         if ( iocc(index,ihole)-delion.lt.1) then
           call wlog(' Cannot remove an electron from this level')
           call par_stop('GETORB-1')
        endif
      endif

c        the recipe for final state atomic configuration is changed
c        from iz+1 prescription, since sometimes it changed occupation
c        numbers in more than two orbitals. This could be consistent
c        only with s02=0.0. New recipe remedy this deficiency.

c     find where to put screening electron
      index1 = index + 1
      do 10  i = 1, 29
 10   if (iscr.eq.0 .and. (iocc(index1,i)-iocc(index,i)).gt.0.5) iscr=i
c     special case of hydrogen like ion
c     if (index.eq.1) iscr=2

c     find where to add or subtract charge delion (iion).
c     if (delion .ge. 0) then
c        removal of electron charge
c        iion is already found
      if (delion .lt. 0) then
c        addition of electron charge
         iion = iscr
c        except special cases
         if (ihole.ne.0 .and.
     1       iocc(index,iscr)+1-delion.gt.2*abs(kappa(iscr))) then
             iion = ilast
             if (ilast.eq.iscr .or. iocc(index,ilast)-delion.gt.
     1                          2*abs(kappa(ilast)) ) iion = ilast + 1
         endif
      endif

      norb = 0
      do 19 i=-4, 3
 19   iorb(i) = 0
      do 20  i = 1, 29
         if (iocc(index,i).gt.0 .or. (i.eq.iscr .and. ihole.gt.0)
     1       .or. (i.eq.iion .and. iocc(index,i)-delion.gt.0))  then
            if (i.ne.ihole .or. iocc(index,i).ge.1) then
               norb = norb + 1
               nqn(norb) = nnum(i)
               nk(norb)  = kappa(i)
               xnel(norb) = iocc(index,i)
               if (i.eq.ihole) then
                  xnel(norb) = xnel(norb) - 1
                  iholep = norb
               endif
               if (i.eq.iscr .and. ihole.gt.0)  xnel(norb)=xnel(norb)+1
               xnval(norb)= ival(index,i)
               if ((kappa(i).eq.-4 .or. kappa(i).eq.3) .and. iunf.eq.0)
     1           xnval(norb) = 0
               xmag(norb) = ispn(index,i)
               iorb(nk(norb)) = norb
               if (i.eq.ihole .and. xnval(norb).ge.1)
     1                         xnval(norb) = xnval(norb) - 1
               if (i.eq.iscr .and. ihole.gt.0) 
     1                         xnval(norb) = xnval(norb) + 1
               if (i.eq.iion)  xnel(norb) = xnel(norb) - delion
               if (i.eq.iion)  xnval(norb) = xnval(norb) - delion
            endif
         endif
   20 continue
      norbco = norb

c     check that all occupation numbers are within limits
      do 50 i = 1, norb
         if ( xnel(i).lt.0 .or.  xnel(i).gt.2*abs(nk(i)) .or.
     1       xnval(i).lt.0 .or. xnval(i).gt.2*abs(nk(i)) ) then
            write (slog,55) i
   55       format(' error in getorb.f. Check occupation number for ',
     1      i3, '-th orbital. May be a problem with ionicity.')
            call wlog(slog)
            call par_stop('GETORB-99')
         endif
  50  continue
c      do 60 i=1,norb
c60    xnval(i) = 0.0d0
c60    xnval(i) = xnel(i)
            
      return
      end
      double precision function getxk (e)
      implicit double precision (a-h, o-z)

c     Make xk from energy(in Hartrees) as
c          k =  sqrt(2*e)  for e > 0  (above the edge)
c          k = -sqrt(-2*e)  for e < 0  (below the edge)
      getxk = sqrt(abs(2*e))
      if (e .lt. 0)  getxk = - getxk
      return
      end
      subroutine sthead (ntitle, title, nph, iz, rmt, rnrm,
     1                  xion, ihole, ixc,
     2                  vr0, vi0, gamach, xmu, xf, vint, rs,
     2                  nohole, lreal,  rgrd)

c     SeT HEAD
c     This routine makes the file header, returned in head array.
c     header lines do not include a leading blank.
c     Last header line is not --------- end-of-header line

c     title lines coming into sthead include carriage control, since
c     they were read from potph.bin

      implicit double precision (a-h, o-z)

c={../HEADERS/const.h
      parameter (pi = 3.14159 26535 89793 23846 26433d0)
      parameter (one = 1, zero = 0)
      parameter (third = one/3)
      parameter (raddeg = 180 / pi)
      complex*16 coni
      parameter (coni = (0,1))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919 158 292 677 512 811d0)

      parameter (bohr = 0.529 177 249d0, ryd  = 13.605 698d0)
      parameter (hart = 2 * ryd)
      parameter (alpinv = 137.035 989 56d0)
c     fine structure alpha
      parameter (alphfs = 1 / alpinv)
c= ../HEADERS/const.h}
c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}
c={../HEADERS/vers.h
      character*12 vfeff
c                       123456789012  
      parameter (vfeff='Feff 8.50   ')
c= ../HEADERS/vers.h}

      dimension xion(0:nphx)
      dimension iz(0:nphx)
      dimension rmt(0:nphx)
      dimension rnrm(0:nphx)

      character*80 title(nheadx), store
      character*16 s1, s2

      character*10 shole(0:29)
      character*8  sout(0:7)
      data shole /'no hole',   'K  shell',  'L1 shell',  'L2 shell',
     2            'L3 shell',  'M1 shell',  'M2 shell',  'M3 shell',
     3            'M4 shell',  'M5 shell',  'N1 shell',  'N2 shell',
     4            'N3 shell',  'N4 shell',  'N5 shell',  'N6 shell',
     5            'N7 shell',  'O1 shell',  'O2 shell',  'O3 shell',
     6            'O4 shell',  'O5 shell',  'O6 shell',  'O7 shell',
     7            'P1 shell',  'P2 shell',  'P3 shell',  'P4 shell',
     8            'P5 shell',  'R1 shell'/
      data sout /'H-L exch', 'D-H exch', 'Gd state', 'DH - HL ',
     1           'DH + HL ', 'val=s+d ', 'sigmd(r)', 'sigmd=c '/


c     Fills head arrray, n = number of lines used.
c     Does not include line of dashes at the end.

      if (ntitle .ge. 1 ) then
         ii = istrln(title(1)) 
         if (ii.gt.1)  then
            write(store,100)  title(1)(1:), vfeff
         else
            write(store,102)  vfeff
         endif
      else
         write(store,102)   vfeff
      endif
  100 format( a55, t66, a12)
  102 format( t66, a12)
      title(1) = store
      nstor = 1

c     remove empty title lines
      do 120  ititle = 2, ntitle
         ii = istrln ( title (ititle) ) 
         if (ii.le.1)  goto 120
         nstor = nstor+1
         title(nstor) = title (ititle)
  120 continue
      ntitle = nstor

c     add more title lines
      if (xion(0) .ne. 0)  then
         ntitle = ntitle + 1
         write(title(ntitle),130)  iz(0), rmt(0)*bohr,
     1                    rnrm(0)*bohr, xion(0), shole(ihole)
      else
         ntitle = ntitle + 1
         write(title(ntitle),140)  iz(0), rmt(0)*bohr,
     1                    rnrm(0)*bohr, shole(ihole)
      endif
  130 format('Abs   Z=',i2, ' Rmt=',f6.3, ' Rnm=',f6.3,
     1       ' Ion=',f5.2,  1x,a10)
  140 format('Abs   Z=',i2,' Rmt=',f6.3,' Rnm=',f6.3, 1x,a10)
c     if (nohole.ge.0)  then
c        ntitle = ntitle + 1
c        write(title(ntitle),142)
c 142    format ('Calculations done with no core hole.')
c     endif
      if (lreal.ge.1 .or. (abs(rgrd - 0.05) .gt. 1.0e-5)) then
        ntitle = ntitle + 1
        s1 = ' '
        if (lreal.gt.1)  then
c        write(title(ntitle),144)
c 144    format ('Calculations done using only real phase shifts.')
         s1 = 'RPHASES'
        elseif (lreal.eq.1) then
c        ntitle = ntitle + 1
c        write(title(ntitle),145)
c 145    format ('Calculations done using only real self energy.')
         s1 = 'RSIGMA'
        endif
        s2 = '  '
        if (abs(rgrd - 0.05) .gt. 1.0e-5)  then
         write(s2,146)  rgrd
  146    format ('  RGRID', f7.4)
        endif
        ilen = istrln(s1)
        title(ntitle) = s1(1:ilen) // s2
      endif

      do 150  iph = 1, nph
         if (xion(iph) .ne. 0)  then
            ntitle = ntitle + 1
            write(title(ntitle),160)  iph, iz(iph),  rmt(iph)*bohr,
     1           rnrm(iph)*bohr, xion(iph)
         else
            ntitle = ntitle + 1
            write(title(ntitle),170)  iph, iz(iph),  rmt(iph)*bohr,
     1           rnrm(iph)*bohr
         endif
  150 continue
  160 format('Pot',i2,' Z=',i2,' Rmt=',f6.3,' Rnm=',f6.3,' Ion=',f5.2)
  170 format('Pot',i2,' Z=',i2,' Rmt=',f6.3,' Rnm=',f6.3)
      if (abs(vi0) .gt. 1.0e-8 .or. abs(vr0) .gt. 1.0e-8)  then
         ntitle = ntitle + 1
         write(title(ntitle),180)  gamach*hart, sout(ixc), vi0*hart,
     1                           vr0*hart
      else
         ntitle = ntitle + 1
         write(title(ntitle),190)  gamach*hart, sout(ixc)
      endif
      ntitle = ntitle + 1
  180 format('Gam_ch=',1pe9.3, 1x,a8, ' Vi=',1pe10.3, ' Vr=',1pe10.3)
  190 format('Gam_ch=',1pe9.3, 1x,a8)
  200 format('Mu=',1pe10.3, ' kf=',1pe9.3, ' Vint=',1pe10.3,
     x        ' Rs_int=',0pf6.3)
      write(title(ntitle),200)  xmu*hart, xf/bohr, vint*hart, rs

      return
      end

      subroutine wthead (io, ntitle, title)
c     Dump title lines to unit io, which must be open. 
      integer io, i, ll
      character*80 title(ntitle)

c     nice for UNIX to use with gnuplot etc.,
      do 310 i = 1, ntitle
         ll = istrln(title(i))
         write(io,300)  title(i)(1:ll)
  300    format (a)
  310 continue

      return
      end
      function itoken (word,flname)
c     chars in word assumed upper case, left justified
c     returns 0 if not a token, otherwise returns token

      character*(*) word
      character*4   w
      character*20 flname
      integer itoken

      w = word(1:4)
      call upper(w)
      
c     Tokens for feff.inp
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (flname(1:8).eq.'feff.inp') then
         if     (w .eq. 'ATOM')  then
            itoken = 1
         elseif (w .eq. 'HOLE')  then
            itoken = 2
         elseif (w .eq. 'OVER')  then
            itoken = 3
         elseif (w .eq. 'CONT')  then
            itoken = 4
         elseif (w .eq. 'EXCH')  then
            itoken = 5
         elseif (w .eq. 'ION ')  then
            itoken = 6
         elseif (w .eq. 'TITL')  then
            itoken = 7
         elseif (w .eq. 'FOLP')  then
            itoken = 8
         elseif (w .eq. 'RPAT' .or. w .eq. 'RMAX')  then
            itoken = 9
         elseif (w .eq. 'DEBY')  then
            itoken = 10
         elseif (w .eq. 'RMUL')  then
            itoken = 11
         elseif (w .eq. 'SS  ')  then
            itoken = 12
         elseif (w .eq. 'PRIN')  then
            itoken = 13
         elseif (w .eq. 'POTE')  then
            itoken = 14
         elseif (w .eq. 'NLEG')  then
            itoken = 15
         elseif (w .eq. 'CRIT')  then
            itoken = 16
         elseif (w .eq. 'NOGE')  then
            itoken = 17
         elseif (w .eq. 'IORD')  then
            itoken = 18
         elseif (w .eq. 'PCRI')  then
            itoken = 19
         elseif (w .eq. 'SIG2')  then
            itoken = 20
         elseif (w .eq. 'XANE')  then
            itoken = 21
         elseif (w .eq. 'CORR')  then
            itoken = 22
         elseif (w .eq. 'AFOL')  then
            itoken = 23
         elseif (w .eq. 'EXAF')  then
            itoken = 24
         elseif (w .eq. 'POLA')  then
            itoken = 25
         elseif (w .eq. 'ELLI')  then
            itoken = 26
         elseif (w .eq. 'RGRI')  then
            itoken = 27
         elseif (w .eq. 'RPHA')  then
            itoken = 28
         elseif (w .eq. 'NSTA')  then
            itoken = 29
         elseif (w .eq. 'NOHO')  then
            itoken = 30
         elseif (w .eq. 'SIG3')  then
            itoken = 31
         elseif (w .eq. 'JUMP')  then
            itoken = 32
         elseif (w .eq. 'MBCO')  then
            itoken = 33
         elseif (w .eq. 'SPIN')  then
            itoken = 34
         elseif (w .eq. 'EDGE')  then
            itoken = 35
         elseif (w .eq. 'SCF ')  then
            itoken = 36
         elseif (w .eq. 'FMS ')  then
            itoken = 37
         elseif (w .eq. 'LDOS')  then
            itoken = 38
         elseif (w .eq. 'INTE')  then
            itoken = 39
         elseif (w .eq. 'CFAV')  then
            itoken = 40
         elseif (w .eq. 'S02 ')  then
            itoken = 41
         elseif (w .eq. 'XES ')  then
            itoken = 42
         elseif (w .eq. 'DANE')  then
            itoken = 43
         elseif (w .eq. 'FPRI')  then
            itoken = 44
         elseif (w .eq. 'RSIG')  then
            itoken = 45
         elseif (w .eq. 'XNCD')  then
            itoken = 46
         elseif (w .eq. 'XMCD')  then
            itoken = 46
         elseif (w .eq. 'MULT')  then
            itoken = 47
         elseif (w .eq. 'UNFR')  then
            itoken = 48
         elseif (w .eq. 'TDLD')  then
            itoken = 49
         elseif (w .eq. 'PMBS')  then
            itoken = 50
         elseif (w .eq. 'PLAS')  then
            itoken = 51
         elseif (w .eq. 'SO2C')  then
            itoken = 52
         elseif (w .eq. 'SELF')  then
            itoken = 53
         elseif (w .eq. 'SFSE')  then
            itoken = 54
         elseif (w .eq. 'RCONV') then
            itoken = 55
         elseif (w .eq. 'ELNE') then !KJ new card for EELS 1-06
            itoken = 56
         elseif (w .eq. 'EXEL') then !KJ new card for EELS 1-06
            itoken = 57
         elseif (w .eq. 'MAGI') then !KJ new card for EELS 1-06
            itoken = 58
         elseif (w .eq. 'ABSO') then !KJ new card 3-06
            itoken = 59    
         elseif (w .eq. 'EGRI')  then !Josh Kas
            itoken = 60
         elseif (w .eq. 'END ')  then
            itoken = -1            
         else
            itoken = 0
         endif
      elseif (flname(1:10).eq.'spring.inp') then
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     These tokens are for spring.inp (input for eq of motion method)
         if (w .eq. 'STRE')  then
            itoken = 1
         elseif (w .eq. 'ANGL')  then
            itoken = 2
         elseif (w .eq. 'VDOS')  then
            itoken = 3
         elseif (w .eq. 'PRDOS') then
            itoken = 4
         elseif (w .eq. 'END ')  then
            itoken = -1            
         else
            itoken = 0
         endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      endif
      return
      end


c====================================================================
      integer function nxtunt(iunit)

c  this function returns the value of the next unopened unit number
c  equal to or larger than iunit.  it will return neither unit numbers
c  0, 5, or 6 nor a negative unit number
c $Id: nxtunt.f,v 1.1.1.1 2006/01/12 06:37:42 hebhop Exp $
c $Log: nxtunt.f,v $
c Revision 1.1.1.1  2006/01/12 06:37:42  hebhop
c New version of feff. feff8.5 (Extension of feff8.4)
c Includes:
c 	1) All feff8.4 capabilities.
c 	2) Screened core hole (calculation of W).
c 	3) Multiple pole self energy calculation.
c 	4) Convolution with spectral function.
c New cards and options:
c 	1) NOHOLE 2      (screened hole)
c 	2) PLASMON ipl   (multiple pole self energy)
c 	3) SO2CONV       (convolve output with spectral function)
c 	4) SELF          (print on shell self energy as calculated by Luke)
c 	5) SFSE k0        (print off shell self energy Sigma(k0,e) )
c
c Revision 1.1.1.1  2000/02/11 02:23:58  alex
c Initialize feff82
c
c Revision 1.10  1999/04/02 21:32:47  newville
c cleaned up nxtunt (matt)
c
c Revision 1.9  1999/02/11 20:08:08  alex
c x39 version: dim.h + misc. small changes
c
c Revision 1.8  1998/12/29 23:59:07  alex
c feff8x35 version
c
c Revision 1.7  1998/11/19 03:23:11  alex
c feff8x32 version
c
c Revision 1.6  1998/10/26 14:11:16  ravel
c no comments beyond column 71
c
c Revision 1.5  1998/10/18 21:47:51  alex
c feff8x30 version implements Broyden algorithm for self-consistency
c
c Revision 1.4  1998/02/24 18:31:37  ravel
c I should really be more careful.  This is the last commitment done
c      cright.
c
c Revision 1.1.1.1  1997/04/27 20:18:03  ravel
c Initial import of xanes sources, version 0.37
c
c Revision 1.1  1996/06/23 16:05:02  bruce
c Initial revision
c

       integer iunit
       logical open

       nxtunt = max(1, iunit) - 1
 10    continue
       nxtunt = nxtunt + 1
       if ((nxtunt.eq.5).or.(nxtunt.eq.6)) nxtunt = 7
       inquire (unit=nxtunt, opened=open)
       if (open) go to 10
       return
c  end integer function nxtunt
       end

c====================================================================
c     Periodic table of the elements
c     Written by Steven Zabinsky, Feb 1992.  Deo Soli Gloria

c     atwts(iz)  single precision fn, returns atomic weight
c     atwtd(iz)  double precision fn, returns atomic weight
c     atsym(iz)  character*2 fn, returns atomic symbol

      double precision function atwtd (iz)
      double precision weight
      common /atwtco/ weight(103)
      atwtd = weight(iz)
      return
      end

      real function atwts (iz)
      double precision weight
      common /atwtco/ weight(103)
      atwts = weight(iz)
      return
      end

      character*2 function atsym (iz)
      character*2 sym
      common /atsyco/ sym(103)
      atsym = sym(iz)
      return
      end

      block data prtbbd
c     PeRiodic TaBle Block Data

c     Atomic weights from inside front cover of Ashcroft and Mermin.

      double precision weight
      common /atwtco/ weight(103)

      character*2 sym
      common /atsyco/ sym(103)

      data weight /
     1   1.0079, 4.0026, 6.941,  9.0122, 10.81,   12.01,
     2   14.007, 15.999, 18.998, 20.18,  22.9898, 24.305,
     3   26.982, 28.086, 30.974, 32.064, 35.453,  39.948,
     4   39.09,  40.08,  44.956, 47.90,  50.942,  52.00,
     5   54.938, 55.85,  58.93,  58.71,  63.55,   65.38,
     6   69.72,  72.59,  74.922, 78.96,  79.91,   83.80,
     7   85.47,  87.62,  88.91,  91.22,  92.91,   95.94,
     8   98.91,  101.07, 102.90, 106.40, 107.87,  112.40,
     9   114.82, 118.69, 121.75, 127.60, 126.90,  131.30,
     x   132.91, 137.34, 138.91, 140.12, 140.91,  144.24,
     1   145,    150.35, 151.96, 157.25, 158.92,  162.50,
     2   164.93, 167.26, 168.93, 173.04, 174.97,  178.49,
     3   180.95, 183.85, 186.2,  190.20, 192.22,  195.09,
     4   196.97, 200.59, 204.37, 207.19, 208.98,  210,
     5   210,    222,    223,    226,    227,     232.04,
     6   231,    238.03, 237.05, 244,    243,     247,
     7   247,    251,    254,    257,    256,     254,
     8   257/

      data sym /  'H', 'He','Li','Be','B', 'C', 'N', 'O', 'F', 'Ne',
     1            'Na','Mg','Al','Si','P', 'S', 'Cl','Ar','K', 'Ca',
     2            'Sc','Ti','V', 'Cr','Mn','Fe','Co','Ni','Cu','Zn',
     3            'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y', 'Zr',
     4            'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',
     5            'Sb','Te','I', 'Xe','Cs','Ba','La','Ce','Pr','Nd',
     6            'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
     7            'Lu','Hf','Ta','W', 'Te','Os','Ir','Pt','Au','Hg',
     8            'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',
     9            'Pa','U', 'Np','Pu','Am','Cm','Bk','Cf','Es','Fm',
     x            'Md','No','Lw'/

      end
      subroutine pijump (ph, old)
      implicit double precision (a-h, o-z)

c     removes jumps of 2*pi in phases

c     ph = current value of phase (may be modified on output, but
c          only by multiples of 2*pi)
c     old = previous value of phase

c={../HEADERS/const.h
      parameter (pi = 3.14159 26535 89793 23846 26433d0)
      parameter (one = 1, zero = 0)
      parameter (third = one/3)
      parameter (raddeg = 180 / pi)
      complex*16 coni
      parameter (coni = (0,1))
c     kf = fa/rs with fa = (9*pi/4)**third, see Ash&Merm, pg 37
      parameter (fa = 1.919 158 292 677 512 811d0)

      parameter (bohr = 0.529 177 249d0, ryd  = 13.605 698d0)
      parameter (hart = 2 * ryd)
      parameter (alpinv = 137.035 989 56d0)
c     fine structure alpha
      parameter (alphfs = 1 / alpinv)
c= ../HEADERS/const.h}
      parameter (twopi = 2 * pi)
      dimension xph(3)

      xph(1) = ph - old
      jump =  (abs(xph(1))+ pi) / twopi
      xph(2) = xph(1) - jump*twopi
      xph(3) = xph(1) + jump*twopi


      xphmin = min (abs(xph(1)), abs(xph(2)), abs(xph(3)))
      isave = 0
      do 10  i = 1, 3
         if (abs (xphmin - abs(xph(i))) .le. 0.01)  isave = i
   10 continue
      if (isave .eq. 0)  then
         call par_stop('pijump')
      endif

      ph = old + xph(isave)

      return
      end
      subroutine rdhead (io, nhead, head, lhead)
      implicit double precision (a-h, o-z)

c     Reads title line(s) from unit io.  Returns number of lines
c     read.  If more than nheadx lines, skips over them.  End-of-header
c     marker is a line of 1 blank, 71 '-'s.
c     lhead is length of each line w/o trailing blanks.
c     header lines returned will have 1st space on line blank for
c     carriage control

      character*80 head(nhead)
      dimension lhead(nhead)
      character*80  line

      n = 0
      nheadx = nhead
      nhead = 0
   10 read(io,20)  line
   20    format(a)
         if (line(4:11) .eq. '--------')  goto 100
         n = n+1
         if (n .le. nheadx)  then
            head(n) = line
            lhead(n) = istrln(head(n))
            nhead = n
         endif
      goto 10
  100 continue
      return
      end
      subroutine rdpot ( ntitle, title, rnrmav, xmu, vint, rhoint,
     1                  emu, s02, erelax, wp, ecv,rs,xf, qtotel, 
     2                  imt, rmt, inrm, rnrm, folp, folpx, xnatph,
     3                  dgc0, dpc0, dgc, dpc, adgc, adpc,
     3                  edens, vclap, vtot, edenvl, vvalgs, dmag, xnval,
     4                  eorb, kappa, iorb, qnrm, xnmues, nohole, ihole,
     5                  inters, totvol, iafolp, xion, iunf, iz, jumprm)
c  opens pot.bin file and reads following information
c  General:
c     ntitle - number of title lines
c     title  - title itself
c     emu    - edge position (x-ray energy for final state at Fermi level)
c  Muffin-tin geometry
c     rmt    - muffin-tin radii
c     imt    - index of radial grid just below muffin-tin radii
c     rnrm   - Norman radii
c     inrm   - index of radial grid just below Norman radii
c     rnrmav - average Norman radius
c     folp   - overlap parameter for rmt
c     folpx  - maximum value for folp
c     xnatph - number of atoms of each potential type
c  Atomic orbitals info (Dirac spinors)
c     dgc0   - upper component for initial orbital
c     dpc0   - lower component for initial orbital
c     dgc    - upper components for all atomic orbitals
c     dpc    - lower components for all atomic orbitals
c     adgc   - development coefficient for upper components
c     adpc   - development coefficient for lower components
c     xnval  - number of valence electrons for each atomic orbital
c              used for core-valence separation and non-local exchange
c     eorb  - atomic enrgy of each orbital for the absorber
c  Electron density information 
c     rhoint - interstitial density
c     rs     - r_s estimate from rhoint (4/3 r_s**3 * rhoint = 1)
c     xf     - estimate of momentum at Fermi level from rhoint
c     edens  - total electron density
c     edenvl - density from valence electrons
c     dmag   - density for spin-up minus density for spin-down
c     qtotel - total charge of a cluster
c     qnrm   - charge accumulated inside Norman sphere as result of SCF
c     xnmues - occupation numbers of valence orbitals from SCF procedure
c  Potential information
c     xmu    - Fermi level position
c     ecv    - core-valence separation energy
c     vint   - muffin-tin zero energy (interstitial potential)
c     vclap  - Coulomb potential
c     vtot   - vclap + xc potential from edens
c     vvalgs - vclap + xc potential from edenvl (EXCHANGE 5 model)
c  Specific data for convolution with excitation spectrum (see mbconv)
c     s02    - many body reduction factor S_0^2 
c     erelax - estimate of relaxation energy = efrozen - emu, where
c              efrozen is edge position estimate from Koopmans theorem
c     wp     - estimate of plasmon frequency from rhoint

      implicit double precision (a-h, o-z)
c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}
      dimension imt(0:nphx), rmt(0:nphx), inrm(0:nphx),  rnrm(0:nphx)
      dimension folp(0:nphx), folpx(0:nphx), dgc0(251), dpc0(251)
      dimension dgc(251, 30, 0:nphx), dpc(251, 30, 0:nphx)
      dimension adgc(10, 30, 0:nphx), adpc(10, 30, 0:nphx)
      dimension edens(251, 0:nphx), vclap(251, 0:nphx)
      dimension vtot(251, 0:nphx), edenvl(251, 0:nphx)
      dimension vvalgs(251, 0:nphx), dmag(251, 0:nphx)
      dimension xnval(30,0:nphx), qnrm(0:nphx), xnmues(0:lx,0:nphx)
      dimension eorb(30), kappa(30)
      dimension iorb(-4:3,0:nphx), iz(0:nphx), xion(0:nphx)
      dimension xnatph(0:nphx)

      character*80 title(nheadx)

      dimension dum(13)

  10  format(a)
   20 format (bn, i15)

      open (unit=3, file='pot.bin', status='old')
      read(3,30) ntitle, nph, npadx, nohole, ihole, inters, iafolp,
     1            jumprm, iunf
  30  format(9(1x,i4))
c     nph and npadx are not passed to calling subroutine
      do 133  i  = 1, ntitle
         read(3,10) title(i)
         call triml(title(i))
  133 continue
c     Misc double precision stuff from pot.bin
      call rdpadd(3, npadx, dum(1), 13)
      rnrmav = dum(1)
      xmu    = dum(2)
      vint   = dum(3)
      rhoint = dum(4)
      emu    = dum(5)
      s02    = dum(6)
      erelax = dum(7)
      wp     = dum(8)
      ecv    = dum(9)
      rs     = dum(10)
      xf     = dum(11)
      qtotel = dum(12)
      totvol = dum(13)

c     read imt
      read (3, 40) (imt(i),i=0,nph)
  40  format(20(1x,i4))
      call rdpadd(3, npadx, rmt(0), nph+1)
c     read inrm
      read (3, 40) (inrm(i),i=0,nph)
      read (3, 40) (iz(i),i=0,nph)
      read (3, 40) (kappa(i),i=1,30)
      call rdpadd(3, npadx, rnrm(0), nph+1)
      call rdpadd(3, npadx, folp(0), nph+1)
      call rdpadd(3, npadx, folpx(0), nph+1)
      call rdpadd(3, npadx, xnatph(0), nph+1)
      call rdpadd(3, npadx, xion(0), nph+1)
      call rdpadd(3, npadx, dgc0(1), 251)
      call rdpadd(3, npadx, dpc0(1), 251)
      call rdpadd(3, npadx, dgc(1,1,0), 251*30*(nph+1) )
      call rdpadd(3, npadx, dpc(1,1,0), 251*30*(nph+1) )
      call rdpadd(3, npadx, adgc(1,1,0), 10*30*(nph+1) )
      call rdpadd(3, npadx, adpc(1,1,0), 10*30*(nph+1) )
      call rdpadd(3, npadx, edens(1,0), 251*(nph+1) )
      call rdpadd(3, npadx, vclap(1,0), 251*(nph+1) )
      call rdpadd(3, npadx, vtot(1,0), 251*(nph+1) )
      call rdpadd(3, npadx, edenvl(1,0), 251*(nph+1) )
      call rdpadd(3, npadx, vvalgs(1,0), 251*(nph+1) )
      call rdpadd(3, npadx, dmag(1,0), 251*(nph+1) )
      call rdpadd(3, npadx, xnval(1,0), 30*(nph+1) )
      call rdpadd(3, npadx, eorb(1), 30)
      do 50 iph=0,nph
 50   read (3, 60) (iorb(i,iph),i=-4,3)
 60   format(8(1x,i2))
      call rdpadd(3, npadx, qnrm(0), nph+1 )
      nn = (lx+1)*(nph+1)
      call rdpadd(3, npadx, xnmues(0,0), nn )
      close (unit=3)

      return
      end
      subroutine rdxsph ( ne, ne1, ne3, nph, ihole, rnrmav,xmu,edge,
     1               ik0, em, eref, iz, potlbl, ph, rkk, lmax, lmaxp1)
      implicit double precision (a-h, o-z)
c     reads file 'phase.bin' 
c  Energy grid information
c     em   - complex energy grid
c     eref - V_int + i*gamach/2 + self-energy correction
c     ne   - total number of points in complex energy grid
c     ne1  - number of points on main horizontal axis
c     ne2  - number of points on vertical vertical axis ne2=ne-ne1-ne3
c     ne3  - number of points on auxilary horizontal axis (need for f')
c     xmu  - Fermi energy
c     edge - x-ray frequency for final state at Fermi level
c     ik0  - grid point index at Fermi level
c  Potential type information
c     nph - number of potential types
c     iz  - charge of nuclei (atomic number)
c     potlbl - label for each potential type
c     lmax - max orb momentum for each potential type
c     ihole - index of core-hole orbital for absorber (iph=0)
c     rnrmav - average Norman radius (used in headers only)
c  Main output of xsect and phases module (except that in xsect.bin)
c     ph  - complex scattering phase shifts
c     rkk - complex multipole matrix elements

c={../HEADERS/dim.h
c      maximum number of atoms for FMS. Reduce nclusx if you need
c      smaller executable.
      parameter (nclusx=100)
c      maximum number of atoms for tdlda module.
      parameter (nclxtd=100)
c      max number of spins: 1 for spin average; 2 for spin-dep
      parameter (nspx=1)
c      max number of atoms in problem for the pathfinder
      parameter (natx =1000)
c      max number of atoms in problem for the rdinp and ffsort
      parameter (nattx =1000)
c      max orbital momentum for FMS module.
      parameter (lx=4)
c      max number of unique potentials (potph) (nphx must be ODD to
c      avoid compilation warnings about alignment in COMMON blocks)
      parameter (nphx = 11)
c      max number of ang mom (arrays 1:ltot+1)
      parameter (ltot = 24)
c      Loucks r grid used through overlap and in phase work arrays
      parameter (nrptx = 1251)
c      Number of energy points genfmt, etc.
      parameter (nex = 150)
c      Max number of distinct lambda's for genfmt
c      15 handles iord 2 and exact ss
      parameter (lamtot=15)
c      vary mmax and nmax independently
      parameter (mtot=4, ntot=2)
c      max number of path atoms, used in path finder, NOT in genfmt
      parameter (npatx = 8)
c      matches path finder, used in GENFMT
      parameter (legtot=npatx+1)
c      max number of overlap shells (OVERLAP card)
      parameter (novrx=8)
c      max number of header lines
      parameter (nheadx=30)
c      max number of poles that can be used to
c      model epsilon^-1 for HL multipole self energy
      parameter (MxPole=1000)
c= ../HEADERS/dim.h}

      character*6  potlbl
      dimension  potlbl(0:nphx)

      complex*16 ph(nex,-ltot:ltot,nspx,0:nphx), eref(nex,nspx), em(nex)
      complex*16 rkk(nex,8,nspx)
      dimension lmax0(0:nphx), lmax(nex,0:nphx)
      dimension iz(0:nphx)
c     kinit, linit, ilinit,  - initial state kappa and ang. mom.
c     lmaxp1  -largest lmax in problem + 1

c     phmin is min value to use for |phase shift|
      parameter (phmin = 1.0d-7)

c     Local staff
c     use temp to write ph, rkk, since ne < nex
      complex*16 temp(nex*(2*ltot+1))
      dimension dum(3)

      open (unit=1, file='phase.bin', status='old', iostat=ios)
      call chopen (ios, 'phase.bin', 'rdxsph')

      read(1,10) nsp, ne, ne1, ne3, nph, ihole, ik0, npadx
  10  format (8(1x,i4))

      call rdpadd(1, npadx, dum(1), 3)
      rnrmav = dum(1)
      xmu    = dum(2)
      edge   = dum(3)

      call rdpadx(1, npadx, em(1), ne)
c     call rdpadx(1, npadx, eref(1), ne)
      call rdpadx (1, npadx, temp(1), ne*nsp)
      ii = 0
      do 60 isp = 1, nsp
      do 60 ie=1, ne
        ii = ii + 1
        eref (ie, isp) = temp(ii)
  60  continue

      do 80  iph = 0, nph
         read(1, 20)  lmax0(iph), iz(iph), potlbl(iph)
  20     format(2(1x,i3), 1x, a6)

         do 75 isp = 1,nsp 
            ii = ne * (2*lmax0(iph)+1)
            call rdpadx (1, npadx, temp(1), ii )
            ii = 0
            do 70  ie = 1, ne
            do 70  ll = -lmax0(iph), lmax0(iph)
               ii = ii+ 1
               ph(ie,ll,isp,iph) = temp(ii)
   70       continue
   75    continue
   80 continue

      call rdpadx (1, npadx, temp(1), ne*8*nsp)
      ii = 0
      do 90 isp = 1,nsp 
      do 90 kdif = 1, 8
      do 90 ie=1, ne
        ii = ii + 1
        rkk (ie, kdif, isp) = temp(ii)
  90  continue

      close (unit=1)

c     make additional data for output
      lmaxp1 = 0
      do 180  iph = 0, nph
      do 180  ie = 1, ne
c        Set lmax to include only non-zero phases
         do 160  il =  lmax0(iph), 0, -1
            lmax(ie,iph) = il
            if (abs(sin(ph(ie, il, 1, iph))) .gt. phmin .or.
     3          abs(sin(ph(ie, il,nsp,iph))) .gt. phmin)  goto 161
  160    continue
  161    continue
         if (lmax(ie,iph)+1 .gt. lmaxp1)  lmaxp1 = lmax(ie,iph)+1
  180 continue

      return
      end
      subroutine setkap(ihole, kinit, linit)
      implicit double precision (a-h, o-z)

c     Set initial state ang mom and quantum number kappa
c     ihole  initial state from ihole    
c     1      K    1s      L=0 -> linit=0 
c     2      LI   2s      L=0 -> linit=0
c     3      LII  2p 1/2  L=1 -> linit=1
c     4      LIII 2p 3/2  L=1 -> linit=1
c     5+     etc.
      if (ihole.le. 2 .or. ihole.eq. 5 .or. ihole.eq.10 .or.
     1    ihole.eq.17 .or. ihole.eq.24 .or. ihole.eq.27)  then
c        hole in s state
         linit = 0
         kinit = -1
      elseif (ihole.eq. 3 .or. ihole.eq. 6 .or. ihole.eq.11 .or.
     1        ihole.eq.18 .or. ihole.eq.25 .or. ihole.eq.30)  then
c        hole in p 1/2 state
         linit = 1
         kinit = 1
      elseif (ihole.eq. 4 .or. ihole.eq. 7 .or. ihole.eq.12 .or.
     1        ihole.eq.19 .or. ihole.eq.26)  then
c        hole in p 3/2 state
         linit = 1
         kinit = -2
      elseif (ihole.eq. 8 .or. ihole.eq.13 .or.
     1        ihole.eq.20 .or. ihole.eq.27)  then
c        hole in d 3/2 state
         linit = 2
         kinit = 2
      elseif (ihole.eq. 9 .or. ihole.eq.14 .or.
     1        ihole.eq.21 .or. ihole.eq.28)  then
c        hole in d 5/2 state
         linit = 2
         kinit = -3
      elseif (ihole.eq.15 .or. ihole.eq.22)  then
c        hole in  f 5/2 state
         linit = 3
         kinit = 3
      elseif (ihole.eq.16 .or. ihole.eq.23)  then
c        hole in  f 7/2 state
         linit = 3
         kinit = -4
      else
c        some unknown hole
         call par_stop('invalid hole number in setkap')
      endif

      return
      end
C FUNCTION ISTRLN (STRING)  Returns index of last non-blank
C                           character.  Returns zero if string is
C                           null or all blank.

      FUNCTION ISTRLN (STRING)
      CHARACTER*(*)  STRING
      CHARACTER BLANK, TAB
      PARAMETER (BLANK = ' ', TAB = '	')
C     there is a tab character here  ^

C  -- If null string or blank string, return length zero.
      ISTRLN = 0
      IF (STRING (1:1) .EQ. CHAR(0))  RETURN
      IF (STRING .EQ. ' ')  RETURN

C  -- Find rightmost non-blank character.
      ILEN = LEN (STRING)
      DO 20  I = ILEN, 1, -1
         IF (STRING(I:I).NE.BLANK .AND. STRING(I:I).NE.TAB)  GOTO 30
   20 CONTINUE
   30 ISTRLN = I

      RETURN
      END
C SUBROUTINE TRIML (STRING)  Removes leading blanks.

      SUBROUTINE TRIML (STRING)
      CHARACTER*(*)  STRING
      CHARACTER*200  TMP
      CHARACTER BLANK, TAB
      PARAMETER (BLANK = ' ', TAB = '	')
C     there is a tab character here  ^

      JLEN = ISTRLN (STRING)

C  -- All blank and null strings are special cases.
      IF (JLEN .EQ. 0)  RETURN

C  -- FInd first non-blank char
      DO 10  I = 1, JLEN
         IF (STRING(I:I).NE.BLANK .AND. STRING(I:I).NE.TAB)  GOTO 20
   10 CONTINUE
   20 CONTINUE

C  -- If I is greater than JLEN, no non-blanks were found.
      IF (I .GT. JLEN)  RETURN

C  -- Remove the leading blanks.
      TMP = STRING (I:)
      STRING = TMP
      RETURN
      END
C SUBROUTINE UPPER (STRING)  Changes a-z to upper case.

      SUBROUTINE UPPER (STRING)
      CHARACTER*(*)  STRING

      JLEN = ISTRLN (STRING)

      DO 10  I = 1, JLEN
         IC = ICHAR (STRING (I:I))
         IF ((IC .LT. 97)  .OR.  (IC .GT. 122))  GOTO 10
         STRING (I:I) = CHAR (IC - 32)
   10 CONTINUE

      RETURN
      END
C SUBROUTINE LOWER (STRING)  Changes A-Z to lower case.

      SUBROUTINE LOWER (STRING)
      CHARACTER*(*)  STRING

      JLEN = ISTRLN (STRING)

      DO 10  I = 1, JLEN
         IC = ICHAR (STRING (I:I))
         IF ((IC .LT. 65) .OR.  (IC .GT. 90))  GOTO 10
         STRING (I:I) = CHAR (IC + 32)
   10 CONTINUE

      RETURN
      END
C***********************************************************************
C
      SUBROUTINE BWORDS (S, NWORDS, WORDS)
C
C     Breaks string into words.  Words are seperated by one or more
C     blanks or tabs, or a comma and zero or more blanks.
C
C     ARGS        I/O      DESCRIPTION
C     ----        ---      -----------
C     S            I       CHAR*(*)  String to be broken up
C     NWORDS      I/O      Input:  Maximum number of words to get
C                          Output: Number of words found
C     WORDS(NWORDS) O      CHAR*(*) WORDS(NWORDS)
C                          Contains words found.  WORDS(J), where J is
C                          greater then NWORDS found, are undefined on
C                          output.
C
C      Written by:  Steven Zabinsky, September 1984
C      Tab char added July 1994.
C
C**************************  Deo Soli Gloria  **************************

C  -- No floating point numbers in this routine.
      IMPLICIT INTEGER (A-Z)

      CHARACTER*(*) S, WORDS(NWORDS)

      CHARACTER BLANK, COMMA, TAB
      PARAMETER (BLANK = ' ', COMMA = ',', TAB = '	')
C     there is a tab character here               ^.

C  -- BETW    .TRUE. if between words
C     COMFND  .TRUE. if between words and a comma has already been found
      LOGICAL BETW, COMFND

C  -- Maximum number of words allowed
      WORDSX = NWORDS

C  -- SLEN is last non-blank character in string
      SLEN = ISTRLN (S)

C  -- All blank string is special case
      IF (SLEN .EQ. 0)  THEN
         NWORDS = 0
         RETURN
      ENDIF

C  -- BEGC is beginning character of a word
      BEGC = 1
      NWORDS = 0

      BETW   = .TRUE.
      COMFND = .TRUE.

      DO 10  I = 1, SLEN
         IF (S(I:I) .EQ. BLANK .OR. S(I:I) .EQ. TAB)  THEN
            IF (.NOT. BETW)  THEN
               NWORDS = NWORDS + 1
               WORDS (NWORDS) = S (BEGC : I-1)
               BETW = .TRUE.
               COMFND = .FALSE.
            ENDIF
         ELSEIF (S(I:I) .EQ. COMMA)  THEN
            IF (.NOT. BETW)  THEN
               NWORDS = NWORDS + 1
               WORDS (NWORDS) = S(BEGC : I-1)
               BETW = .TRUE.
            ELSEIF (COMFND)  THEN
               NWORDS = NWORDS + 1
               WORDS (NWORDS) = BLANK
            ENDIF
            COMFND = .TRUE.
         ELSE
            IF (BETW)  THEN
               BETW = .FALSE.
               BEGC = I
            ENDIF
         ENDIF

         IF (NWORDS .GE. WORDSX)  RETURN

   10 CONTINUE

      IF (.NOT. BETW  .AND.  NWORDS .LT. WORDSX)  THEN
         NWORDS = NWORDS + 1
         WORDS (NWORDS) = S (BEGC :SLEN)
      ENDIF

      RETURN
      END
      SUBROUTINE UNTAB (STRING)
C REPLACE TABS WITH BLANKS :    TAB IS ASCII DEPENDENT
      INTEGER        ITAB , I, ILEN, ISTRLN
      PARAMETER      (ITAB = 9)
      CHARACTER*(*)  STRING, TAB*1
      EXTERNAL ISTRLN
      TAB  = CHAR(ITAB)
      ILEN = MAX(1, ISTRLN(STRING))
 10   CONTINUE 
        I = INDEX(STRING(:ILEN), TAB ) 
        IF (I .NE. 0) THEN
            STRING(I:I) = ' '
            GO TO 10
        END IF
      RETURN
C END SUBROUTINE UNTAB
      END

      logical function iscomm (line)
c     returns true if line is a comment or blank line, false otherwise
c#mn{ rewritten to allow ";*%#" as comment characters
       character*(*) line
       iscomm = ((line.eq.' ').or.(index(';*%#',line(1:1)).ne.0))
c#mn}
      return
      end
      subroutine str2dp(str,dpval,ierr)
c  return dp number "dpval" from character string "str"
c  if str cannot be a number, ierr < 0 is returned.
      character*(*) str, fmt*15 
      double precision dpval
      integer  ierr , lenmax
      parameter ( lenmax = 40)
      logical  isnum
      external isnum
      ierr = -99
      if (isnum(str)) then
         ierr = 0
         write(fmt, 10) min(lenmax, len(str))
 10      format('(bn,f',i3,'.0)')
         read(str, fmt, err = 20, iostat=ierr) dpval
      end if    
      if (ierr.gt.0) ierr = -ierr
      return
 20   continue
      ierr = -98
      return
c end subroutine str2dp
      end
      subroutine str2re(str,val,ierr)
c  return real from character string "str"
      character*(*) str 
      double precision dpval
      real     val
      integer  ierr
      call str2dp(str,dpval,ierr)
      if (ierr.eq.0) val = dpval
      return
c end subroutine str2re
      end
      subroutine str2in(str,intg,ierr)
c  return integer from character string "str"
c  returns ierr = 1 if value was clearly non-integer
      character*(*) str 
      double precision val, tenth
      parameter (tenth = 1.d-1)
      integer  ierr, intg
      call str2dp(str,val,ierr)
      if (ierr.eq.0) then
         intg = int(val)
         if ((abs(intg - val) .gt. tenth))  ierr = 1
       end if
      return
c end subroutine str2in
      end
       logical function isnum (string)
c  tests whether a string can be a number. not foolproof!
c  to return true, string must contain:
c    - only characters in  'deDE.+-, 1234567890' (case is checked)
c    - no more than one 'd' or 'e' 
c    - no more than one '.'
c  matt newville
       character*(*)  string,  number*20
c note:  layout and case of *number* is important: do not change!
       parameter (number = 'deDE.,+- 1234567890')
       integer   iexp, idec, i, j, istrln
       external  istrln
       iexp  = 0
       idec  = 0
       isnum = .false. 
       do 100  i = 1, max(1, istrln(string))
          j = index(number,string(i:i))
          if (j.le.0)               go to 200
          if((j.ge.1).and.(j.le.4)) iexp = iexp + 1
          if (j.eq.5)               idec = idec + 1
 100   continue
c  every character in "string" is also in "number".  so, if there are 
c  not more than one exponential and decimal markers, it's a number
       if ((iexp.le.1).and.(idec.le.1)) isnum = .true.
 200   continue
       return
c  end logical function isnum
       end
      subroutine wlog (string)
      character*(*) string

c={../HEADERS/parallel.h
      integer par_type, this_process, numprocs, my_rank
      logical master, worker, parallel_run
      real*8 wall_comm, time_comm
      common /timing/ wall_comm, time_comm
      common /parallel/ numprocs, my_rank, this_process, 
     .          master, worker, parallel_run, par_type
c= ../HEADERS/parallel.h}

c     This output routine is used to replace the PRINT statement
c     for output that "goes to the terminal", or to the log file.
c     If you use a window based system, you can modify this routine
c     to handle the running output elegantly.
c     Handle carriage control in the string you pass to wlog.
c
c     The log file is also written here, hard coded here.

c     The log file is unit 11.  The log file is opened in the
c     main program, program feff.

c     make sure not to write trailing blanks

   10 format (a)

c     Suppress output in sequential loops
      if (par_type .eq. 2) return

      il = istrln (string)
      if (il .eq. 0)  then
         print10
         if (par_type .ne. 3) write(11,10)
      else
         print10, string(1:il)
         if (par_type .ne. 3) write(11,10) string(1:il)
      endif
      return
      end
      subroutine lblank (string)
      character*(*) string
c     add a leading blank, useful for carriage control
      string = ' ' // string
      return
      end
      double precision function xx (j)
      implicit double precision (a-h, o-z)
c     x grid point at index j, x = log(r), r=exp(x)
      parameter (delta = 0.050 000 000 000 000)
      parameter (c88   = 8.800 000 000 000 000)
c     xx = -8.8 + (j-1)*0.05
      xx = -c88 + (j-1)*delta
      return
      end

      double precision function rr(j)
      implicit double precision (a-h, o-z)
c     r grid point at index j
      rr = exp (xx(j))
      return
      end

      function ii(r)
      implicit double precision (a-h, o-z)
c     index of grid point immediately below postion r
      parameter (delta = 0.050 000 000 000 000)
      parameter (c88   = 8.800 000 000 000 000)
c     ii = (log(r) + 8.8) / 0.05 + 1
      ii = (log(r) + c88) / delta + 1
      return
      end
c
c PAD library:   Packed Ascii Data 
c   these routines contain code for handling packed-ascii-data  
c   (pad) arrays for writing printable character strings that 
c   represent real or complex scalars and arrays to a file.
c
c routines included in padlib are (dp==double precision):
c   wrpadd     write a dp array as pad character strings
c   wrpadx     write a dp complex array as pad character strings
c   rdpadr     read a pad character array as a real array
c   rdpadd     read a pad character array as a dp  array
c   rdpadc     read a pad character array as a complex array
c   rdpadx     read a pad character array as a dp complex array
c   pad        internal routine to convert dp number to pad string
c   unpad      internal routine to pad string to dp number
c
c routines not included, but required by padlib:
c     triml, istrln, wlog
c
c//////////////////////////////////////////////////////////////////////
c Copyright (c) 1997--2001 Matthew Newville, The University of Chicago
c Copyright (c) 1992--1996 Matthew Newville, University of Washington
c
c Permission to use and redistribute the source code or binary forms of
c this software and its documentation, with or without modification is
c hereby granted provided that the above notice of copyright, these
c terms of use, and the disclaimer of warranty below appear in the
c source code and documentation, and that none of the names of The
c University of Chicago, The University of Washington, or the authors
c appear in advertising or endorsement of works derived from this
c software without specific prior written permission from all parties.
c
c THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
c EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
c MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
c IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
c CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
c TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
c SOFTWARE OR THE USE OR OTHER DEALINGS IN THIS SOFTWARE.
c//////////////////////////////////////////////////////////////////////
c License is applicable for routines below, until otherwise specified.
c
       subroutine wrpadd(iout,npack,array,npts)
c
c write a dp array to a file in packed-ascii-data format
c
c inputs:  [ no outputs / no side effects ]
c   iout   unit to write to (assumed open)
c   npack  number of characters to use (determines precision)
c   array  real array 
c   npts   number of array elements to read
c notes:
c   real number converted to packed-ascii-data string using pad
c={padlib.h
c padlib.h -*-fortran-*-
c  header of parameters for packed-ascii-data (pad) routines
       implicit none
       character cpadr, cpadi, cpadc
       integer   maxlen, ibase, ioff, ihuge, ibas2
       double precision  ten, tenlog, huge, tiny, one, zero, base
       parameter(cpadr = '!', cpadc = '$', cpadi = '%')
       parameter(ibase = 90, ioff = 37, ihuge = 38, maxlen = 82)
       parameter(ibas2 = ibase/2, zero=0d0, one=1.d0, ten = 10.d0)
       parameter(tenlog= 2.302585092994045684d0)
       parameter(huge = ten**ihuge, tiny = one/huge)
       parameter(base = ibase*one)
c
c= padlib.h}
       integer    iout, npack, npts, mxl, js, i
       character  str*128
       double precision array(*), xr
       js  = 0
       str = ' '
       mxl = maxlen - npack + 1
       do 20 i = 1, npts
          js = js+npack
          xr = array(i)
          call pad(xr, npack, str(js-npack+1:js))
          if ((js.ge.mxl).or.(i.eq.npts)) then
             write(iout,100) cpadr, str(1:js)
             js = 0
          end if
 20    continue
       return
 100   format(a1,a)
       end
c --padlib--
       subroutine wrpadx(iout,npack,array,npts)
c write complex*16 array as pad string
c={padlib.h
c padlib.h -*-fortran-*-
c  header of parameters for packed-ascii-data (pad) routines
       implicit none
       character cpadr, cpadi, cpadc
       integer   maxlen, ibase, ioff, ihuge, ibas2
       double precision  ten, tenlog, huge, tiny, one, zero, base
       parameter(cpadr = '!', cpadc = '$', cpadi = '%')
       parameter(ibase = 90, ioff = 37, ihuge = 38, maxlen = 82)
       parameter(ibas2 = ibase/2, zero=0d0, one=1.d0, ten = 10.d0)
       parameter(tenlog= 2.302585092994045684d0)
       parameter(huge = ten**ihuge, tiny = one/huge)
       parameter(base = ibase*one)
c
c= padlib.h}
       integer    iout, npack, npts, mxl, js, i
       complex*16 array(*)
       character  str*128
       double precision xr, xi
       js = 0
       str  = ' '
       mxl  = maxlen - 2 * npack + 1
       do 20 i = 1, npts
          js = js  + 2 * npack
          xr = dble(array(i))
          xi = dimag(array(i))
          call pad(xr, npack, str(js-2*npack+1:js-npack))
          call pad(xi, npack, str(js-npack+1:js))
          if ((js.ge.mxl).or.(i.eq.npts)) then
             write(iout,100) cpadc, str(1:js)
             js = 0
          end if
 20    continue
       return
 100   format(a1,a)
       end
c --padlib--
       subroutine wrpadr(iout,npack,array,npts)
c
c write a real array to a file in packed-ascii-data format
c
c inputs:  [ no outputs / no side effects ]
c   iout   unit to write to (assumed open)
c   npack  number of characters to use (determines precision)
c   array  real array 
c   npts   number of array elements to read
c notes:
c   real number converted to packed-ascii-data string using pad
c={padlib.h
c padlib.h -*-fortran-*-
c  header of parameters for packed-ascii-data (pad) routines
       implicit none
       character cpadr, cpadi, cpadc
       integer   maxlen, ibase, ioff, ihuge, ibas2
       double precision  ten, tenlog, huge, tiny, one, zero, base
       parameter(cpadr = '!', cpadc = '$', cpadi = '%')
       parameter(ibase = 90, ioff = 37, ihuge = 38, maxlen = 82)
       parameter(ibas2 = ibase/2, zero=0d0, one=1.d0, ten = 10.d0)
       parameter(tenlog= 2.302585092994045684d0)
       parameter(huge = ten**ihuge, tiny = one/huge)
       parameter(base = ibase*one)
c
c= padlib.h}
       integer    iout, npack, npts, mxl, js, i
       character  str*128
       real    array(*)
       double precision xr
       js  = 0
       str = ' '
       mxl = maxlen - npack + 1
       do 20 i = 1, npts
          js = js+npack
          xr = dble(array(i))
          call pad(xr, npack, str(js-npack+1:js))
          if ((js.ge.mxl).or.(i.eq.npts)) then
             write(iout,100) cpadr, str(1:js)
             js = 0
          end if
 20    continue
       return
 100   format(a1,a)
       end
c --padlib--
       subroutine wrpadc(iout,npack,array,npts)
c write complex (*8) array as pad string
c={padlib.h
c padlib.h -*-fortran-*-
c  header of parameters for packed-ascii-data (pad) routines
       implicit none
       character cpadr, cpadi, cpadc
       integer   maxlen, ibase, ioff, ihuge, ibas2
       double precision  ten, tenlog, huge, tiny, one, zero, base
       parameter(cpadr = '!', cpadc = '$', cpadi = '%')
       parameter(ibase = 90, ioff = 37, ihuge = 38, maxlen = 82)
       parameter(ibas2 = ibase/2, zero=0d0, one=1.d0, ten = 10.d0)
       parameter(tenlog= 2.302585092994045684d0)
       parameter(huge = ten**ihuge, tiny = one/huge)
       parameter(base = ibase*one)
c
c= padlib.h}
       integer    iout, npack, npts, mxl, js, i
       complex    array(*)
       character  str*128
       double precision xr, xi
       js = 0
       str  = ' '
       mxl  = maxlen - 2 * npack + 1
       do 20 i = 1, npts
          js = js  + 2 * npack
          xr = dble(array(i))
          xi = aimag(array(i))
          call pad(xr, npack, str(js-2*npack+1:js-npack))
          call pad(xi, npack, str(js-npack+1:js))
          if ((js.ge.mxl).or.(i.eq.npts)) then
             write(iout,100) cpadc, str(1:js)
             js = 0
          end if
 20    continue
       return
 100   format(a1,a)
       end
c --padlib--
       subroutine rdpadd(iou,npack,array,npts)
c read dparray from packed-ascii-data file
c arguments:
c   iou    unit to read from (assumed open)                   (in)
c   npack  number of characters to use (determines precision) (in)
c   array  real array                                         (out)
c   npts   number of array elements to read / number read     (in/out)
c notes:
c   packed-ascii-data string converted to real array using  unpad
c={padlib.h
c padlib.h -*-fortran-*-
c  header of parameters for packed-ascii-data (pad) routines
       implicit none
       character cpadr, cpadi, cpadc
       integer   maxlen, ibase, ioff, ihuge, ibas2
       double precision  ten, tenlog, huge, tiny, one, zero, base
       parameter(cpadr = '!', cpadc = '$', cpadi = '%')
       parameter(ibase = 90, ioff = 37, ihuge = 38, maxlen = 82)
       parameter(ibas2 = ibase/2, zero=0d0, one=1.d0, ten = 10.d0)
       parameter(tenlog= 2.302585092994045684d0)
       parameter(huge = ten**ihuge, tiny = one/huge)
       parameter(base = ibase*one)
c
c= padlib.h}
       integer iou, npack, npts, ndline, i, istrln, ipts, iread
       double precision    array(*), unpad , tmp
       character  ctest, ccomp
       character  str*128
       external  unpad, istrln, iread
       ccomp = cpadr
       ipts = 0
 10    continue 
          i = iread(iou, str)
          if (i.lt.0) go to 50
          call triml(str)
          ctest  = str(1:1)
          str    = str(2:)
          ndline = i/npack
          if ((ctest.ne.ccomp).or.(ndline.le.0)) go to 200
          do 30 i = 1, ndline
             ipts  = ipts + 1
             tmp   = unpad(str(1-npack+i*npack:i*npack),npack)
             array(ipts) = tmp
             if (ipts.ge.npts) go to 50
 30       continue 
          go to 10
 50    continue 
       return
 200   continue
       call wlog (' -- Read_PAD error:  bad data at line:')
       i = istrln(str)
       call wlog (str(:i))
       stop ' -- fatal error in reading PAD data file -- '
       end
c --padlib--
       subroutine rdpadr(iou,npack,array,npts)
c read real array from packed-ascii-data file
c arguments:
c   iou    unit to read from (assumed open)                   (in)
c   npack  number of characters to use (determines precision) (in)
c   array  real array                                         (out)
c   npts   number of array elements to read / number read     (in/out)
c notes:
c   packed-ascii-data string converted to real array using  unpad
c={padlib.h
c padlib.h -*-fortran-*-
c  header of parameters for packed-ascii-data (pad) routines
       implicit none
       character cpadr, cpadi, cpadc
       integer   maxlen, ibase, ioff, ihuge, ibas2
       double precision  ten, tenlog, huge, tiny, one, zero, base
       parameter(cpadr = '!', cpadc = '$', cpadi = '%')
       parameter(ibase = 90, ioff = 37, ihuge = 38, maxlen = 82)
       parameter(ibas2 = ibase/2, zero=0d0, one=1.d0, ten = 10.d0)
       parameter(tenlog= 2.302585092994045684d0)
       parameter(huge = ten**ihuge, tiny = one/huge)
       parameter(base = ibase*one)
c
c= padlib.h}
       integer iou, npack, npts, ndline, i, istrln, ipts, iread
       real    array(*)
       double precision unpad , tmp
       character  ctest, ccomp
       character  str*128
       external  unpad, istrln, iread
       ccomp = cpadr
       ipts = 0
 10    continue 
          i = iread(iou, str)
          if (i.lt.0) go to 50
          call triml(str)
          ctest  = str(1:1)
          str    = str(2:)
          ndline = i/npack
          if ((ctest.ne.ccomp).or.(ndline.le.0)) go to 200
          do 30 i = 1, ndline
             ipts  = ipts + 1
             tmp   = unpad(str(1-npack+i*npack:i*npack),npack)
             array(ipts) = real(tmp)
             if (ipts.ge.npts) go to 50
 30       continue 
          go to 10
 50    continue 
       return
 200   continue
       call wlog (' -- Read_PAD error:  bad data at line:')
       i = istrln(str)
       call wlog (str(:i))
       stop ' -- fatal error in reading PAD data file -- '
       end
c --padlib--
       subroutine rdpadc(iou,npack,array,npts)
c read complex array from packed-ascii-data file
c arguments:
c   iou    unit to read from (assumed open)                  (in)
c   npack  number of characters to use (determines precision)(in)
c   array  complex array                                     (out)
c   npts   number of array elements to read / number read    (in/out)
c notes:
c   packed-ascii-data string converted to real array using  unpad
c={padlib.h
c padlib.h -*-fortran-*-
c  header of parameters for packed-ascii-data (pad) routines
       implicit none
       character cpadr, cpadi, cpadc
       integer   maxlen, ibase, ioff, ihuge, ibas2
       double precision  ten, tenlog, huge, tiny, one, zero, base
       parameter(cpadr = '!', cpadc = '$', cpadi = '%')
       parameter(ibase = 90, ioff = 37, ihuge = 38, maxlen = 82)
       parameter(ibas2 = ibase/2, zero=0d0, one=1.d0, ten = 10.d0)
       parameter(tenlog= 2.302585092994045684d0)
       parameter(huge = ten**ihuge, tiny = one/huge)
       parameter(base = ibase*one)
c
c= padlib.h}
       integer iou, npack,npts, ndline, i, istrln, ipts, np, iread
       double precision  unpad, tmpr, tmpi
       complex  array(*)
       character  ctest, ccomp
       character  str*128
       external  unpad, istrln, iread
       ccomp = cpadc
       ipts = 0
       np   = 2 * npack
 10    continue 
          i = iread(iou, str)
          if (i.lt.0) go to 50
          call triml(str)
          ctest  = str(1:1)
          str    = str(2:)
          ndline = i / np
          if ((ctest.ne.ccomp).or.(ndline.le.0)) go to 200
          do 30 i = 1, ndline
             ipts = ipts + 1
             tmpr = unpad(str(1-np+i*np:-npack+i*np),npack)
             tmpi = unpad(str(1-npack+i*np:i*np),npack)
             array(ipts) = cmplx(tmpr, tmpi)
             if (ipts.ge.npts) go to 50
 30       continue 
          go to 10
 50    continue 
       return
 200   continue
       call wlog (' -- Read_PAD error:  bad data at line:')
       i = istrln(str)
       call wlog (str(:i))
       stop ' -- fatal error in reading PAD data file -- '
       end
       subroutine rdpadx(iou,npack,array,npts)
c read complex*16 array from packed-ascii-data file
c arguments:
c   iou    unit to read from (assumed open)                  (in)
c   npack  number of characters to use (determines precision)(in)
c   array  complex array                                     (out)
c   npts   number of array elements to read / number read    (in/out)
c notes:
c   packed-ascii-data string converted to real array using  unpad
c={padlib.h
c padlib.h -*-fortran-*-
c  header of parameters for packed-ascii-data (pad) routines
       implicit none
       character cpadr, cpadi, cpadc
       integer   maxlen, ibase, ioff, ihuge, ibas2
       double precision  ten, tenlog, huge, tiny, one, zero, base
       parameter(cpadr = '!', cpadc = '$', cpadi = '%')
       parameter(ibase = 90, ioff = 37, ihuge = 38, maxlen = 82)
       parameter(ibas2 = ibase/2, zero=0d0, one=1.d0, ten = 10.d0)
       parameter(tenlog= 2.302585092994045684d0)
       parameter(huge = ten**ihuge, tiny = one/huge)
       parameter(base = ibase*one)
c
c= padlib.h}
       integer iou, npack,npts, ndline, i, istrln, ipts, np, iread
       double precision  unpad, tmpr, tmpi
       complex*16  array(*)
       character  ctest, ccomp
       character  str*128
       external  unpad, istrln, iread
       ccomp = cpadc
       ipts = 0
       np   = 2 * npack
 10    continue 
          i = iread(iou, str)
          if (i.lt.0) go to 50
          call triml(str)
          ctest  = str(1:1)
          str    = str(2:)
          ndline = i / np
          if ((ctest.ne.ccomp).or.(ndline.le.0)) go to 200
          do 30 i = 1, ndline
             ipts = ipts + 1
             tmpr = unpad(str(1-np+i*np:-npack+i*np),npack)
             tmpi = unpad(str(1-npack+i*np:i*np),npack)
             array(ipts) = cmplx(tmpr, tmpi)
             if (ipts.ge.npts) go to 50
 30       continue 
          go to 10
 50    continue 
       return
 200   continue
       call wlog (' -- Read_PAD error:  bad data at line:')
       i = istrln(str)
       call wlog (str(:i))
       stop ' -- fatal error in reading PAD data file -- '
       end

c --padlib--
       subroutine pad(xreal,npack,str)
c  convert dp number *xreal* to packed-ascii-data string *str*
c={padlib.h
c padlib.h -*-fortran-*-
c  header of parameters for packed-ascii-data (pad) routines
       implicit none
       character cpadr, cpadi, cpadc
       integer   maxlen, ibase, ioff, ihuge, ibas2
       double precision  ten, tenlog, huge, tiny, one, zero, base
       parameter(cpadr = '!', cpadc = '$', cpadi = '%')
       parameter(ibase = 90, ioff = 37, ihuge = 38, maxlen = 82)
       parameter(ibas2 = ibase/2, zero=0d0, one=1.d0, ten = 10.d0)
       parameter(tenlog= 2.302585092994045684d0)
       parameter(huge = ten**ihuge, tiny = one/huge)
       parameter(base = ibase*one)
c
c= padlib.h}
       integer  iexp, itmp, isgn, i, npack, j
       double precision xreal, xwork, xsave,onem, tenth
       parameter (onem  =  0.99999999997d0)
       parameter (tenth =  0.099999999994d0)
       character str*(*)
c
       str      = ' '
       xsave    = min(huge, max(-huge, xreal))
       isgn     = 1
       if (xsave.le.0) isgn = 0
c
       xwork    = dabs( xsave )
       iexp     = 0
       if ((xwork.lt.huge).and.(xwork.gt.tiny))  then
          iexp  =   1 + int(log(xwork) / tenlog  )
       else if (xwork.ge.huge) then
          iexp  = ihuge
          xwork = one
       else if (xwork.le.tiny)  then
          xwork = zero
       end if
c force xwork between ~0.1 and ~1
c note: this causes a loss of precision, but 
c allows backward compatibility
       xwork    = xwork / (ten ** iexp)
 20    continue
       if (xwork.ge.one) then
          xwork = xwork * 0.100000000000000d0
          iexp  = iexp + 1
       else if (xwork.le.tenth) then
          xwork = xwork * ten
          iexp  = iexp - 1
       endif
       if (xwork.ge.one) go to 20

       itmp     = int ( ibas2 * xwork ) 
       str(1:1) = char(iexp  + ioff + ibas2 )
       str(2:2) = char( 2 * itmp + isgn + ioff)
       xwork    = xwork * ibas2 - itmp
       if (npack.gt.2) then
          do 100 i = 3, npack
             itmp     = int( base * xwork + 1.d-9)
             str(i:i) = char(itmp + ioff)
             xwork    = xwork * base - itmp
 100      continue
       end if
       if (xwork.ge.0.5d0) then
          i = itmp + ioff + 1
          if (i.le.126) then
             str(npack:npack)= char(i)
          else 
             j = ichar(str(npack-1:npack-1))
             if (j.lt.126) then
                str(npack-1:npack-1) = char(j+1)
                str(npack:npack)     = char(37)
             endif 
          endif
       endif
       return
       end
c --padlib--
       double precision function unpad(str,npack)
c
c  convert packed-ascii-data string *str* to dp number *unpad*
c={padlib.h
c padlib.h -*-fortran-*-
c  header of parameters for packed-ascii-data (pad) routines
       implicit none
       character cpadr, cpadi, cpadc
       integer   maxlen, ibase, ioff, ihuge, ibas2
       double precision  ten, tenlog, huge, tiny, one, zero, base
       parameter(cpadr = '!', cpadc = '$', cpadi = '%')
       parameter(ibase = 90, ioff = 37, ihuge = 38, maxlen = 82)
       parameter(ibas2 = ibase/2, zero=0d0, one=1.d0, ten = 10.d0)
       parameter(tenlog= 2.302585092994045684d0)
       parameter(huge = ten**ihuge, tiny = one/huge)
       parameter(base = ibase*one)
c
c= padlib.h}
       double precision sum
       integer   iexp, itmp, isgn, i, npack
       character str*(*)
       unpad = zero
       if (npack.le.2) return
       iexp  =     (ichar(str(1:1)) - ioff   ) - ibas2
       isgn  = mod (ichar(str(2:2)) - ioff, 2) * 2 - 1
       itmp  =     (ichar(str(2:2)) - ioff   ) / 2
       sum   = dble(itmp/(base*base))
c       do 100 i = 3, npack
c          sum = sum + dble(ichar(str(i:i)) - ioff) / base**i
c 100   continue
       do 100 i = npack, 3, -1
          sum = sum + dble(ichar(str(i:i)) - ioff) / base**i
 100   continue
       unpad = 2 * isgn * ibase * sum * (ten ** iexp)
cc       print*, sum, iexp,unpad
       return
       end
c --padlib--
c end of pad library
c ----------
       integer function iread(lun,string)
c
c generalized internal read:
c    read a string the next line of an opened file 
c    unit, returning the real length of string
c 
c inputs:   
c   lun     opened file unit number
c outputs:
c   string  string read from file
c returns:
c   iread   useful length of string, as found from 
c                  sending string to 'sclean' to 
c                  remove non-printable characters
c                   and then istrln  
c           or
c              -1   on 'end-of-file'
c              -2   on 'error'
c
c copyright (c) 1999  Matthew Newville
       implicit none
       character*(*) string
       integer    lun, istrln
       external   istrln
       string = ' '
 10    format(a)
       read(lun, 10, end = 40, err = 50) string
       call sclean(string)
       iread = istrln(string)
       return
 40    continue 
       string = ' '
       iread = -1
       return
 50    continue 
       string = ' '
       iread = -2
       return
       end
       subroutine sclean(str) 
c
c  clean a string, especially for strings passed between 
c  different file systems, or from C functions:
c
c   1. characters in the range char(0), or char(10)...char(15) 
c      are interpreted as end-of-line characters, so that all
c      remaining characters are explicitly blanked.
c   2. all other characters below char(31) (including tab) are
c      replaced by a single blank
c
c  this is mostly useful when getting a string generated by a C 
c  function and for handling dos/unix/max line-endings.
c
c copyright (c) 1999  Matthew Newville
       character*(*) str, blank*1
       parameter (blank = ' ')
       integer i,j,is
       do 20 i = 1, len(str)
          is = ichar(str(i:i))
          if ((is.eq.0) .or. ((is.ge.10) .and. (is.le.15))) then
             do 10 j= i, len(str)
                str(j:j) = blank
 10          continue
             return
          elseif (is.le.31)  then
             str(i:i)  = blank
          end if
 20    continue 
       return
c end subroutine sclean
       end

      SUBROUTINE rdcmt(iUnt,Cmt)
      INTEGER iUnt, i1
      CHARACTER(300) line
      CHARACTER(4) Cmt
      CHARACTER TmpCmt(4), ch
      LOGICAL CmtLin

      CmtLin = .true.
      DO i1 = 1, 4
         TmpCmt(i1) = Cmt(i1:i1)
      END DO
 5    CONTINUE
      READ(iUnt,*,END=10) ch
      DO i1 = 1, 4
         IF(ch.eq.TmpCmt(i1)) goto 5
      END DO
      
      BACKSPACE(iUnt)
      
 10   CONTINUE
      
      RETURN
      END
      subroutine setgam (iz, ihole, gamach)

c     Sets gamach, core hole lifetime.  Data comes from graphs in
c     K. Rahkonen and K. Krause,
c     Atomic Data and Nuclear Data Tables, Vol 14, Number 2, 1974.
c     output gamach is in eV

      implicit double precision (a-h, o-z)

      dimension zh(8,16), gamh(8,16)

      dimension zk(8), gamkp(8)
      parameter (ryd  = 13.605 698d0)
      parameter (hart = 2*ryd)
      character*512 slog


c     Note that 0.99 replaces 1.0, 95.1 replaces 95.0 to avoid roundoff
c     trouble.
c     Gam arrays contain the gamma values.
c     We will take log10 of the gamma values so we can do linear
c     interpolation from a log plot.

      data  zh   / 0.99,  10.0, 20.0, 40.0, 50.0, 60.0, 80.0, 95.1,
     2              0.99, 18.0, 22.0, 35.0, 50.0, 52.0, 75.0,  95.1,
     3              0.99,  17.0, 28.0, 31.0, 45.0, 60.0,  80.0, 95.1,
     4              0.99,  17.0, 28.0, 31.0, 45.0, 60.0,  80.0, 95.1,
     5              0.99,  20.0, 28.0, 30.0, 36.0, 53.0,  80.0, 95.1,
     6              0.99,  20.0, 22.0, 30.0, 40.0, 68.0,  80.0, 95.1,
     7              0.99,  20.0, 22.0, 30.0, 40.0, 68.0,  80.0, 95.1,
     8              0.99,  36.0, 40.0, 48.0, 58.0, 76.0,  79.0, 95.1,
     9              0.99,  36.0, 40.0, 48.0, 58.0, 76.0,  79.0, 95.1,
     *              0.99,  30.0, 40.0, 47.0, 50.0, 63.0,  80.0, 95.1,
     1              0.99,  40.0, 42.0, 49.0, 54.0, 70.0,  87.0, 95.1,
     2              0.99,  40.0, 42.0, 49.0, 54.0, 70.0,  87.0, 95.1,
     3              0.99,  40.0, 50.0, 55.0, 60.0, 70.0,  81.0, 95.1,
     4              0.99,  40.0, 50.0, 55.0, 60.0, 70.0,  81.0, 95.1,
     5              0.99,  71.0, 73.0, 79.0, 86.0, 90.0,  95.0,100.0,
     6              0.99,  71.0, 73.0, 79.0, 86.0, 90.0,  95.0,100.0/

      data  gamh / 0.02,  0.28, 0.75,  4.8, 10.5, 21.0, 60.0, 105.0,
     2              0.07,  3.9,  3.8,  7.0,  6.0,  3.7,  8.0,  19.0,
     3              0.001, 0.12,  1.4,  0.8,  2.6,  4.1,   6.3, 10.5,
     4              0.001, 0.12, 0.55,  0.7,  2.1,  3.5,   5.4,  9.0,
     5              0.001,  1.0,  2.9,  2.2,  5.5, 10.0,  22.0, 22.0,
     6              0.001,0.001,  0.5,  2.0,  2.6, 11.0,  15.0, 16.0,
     7              0.001,0.001,  0.5,  2.0,  2.6, 11.0,  10.0, 10.0,
     8              0.0006,0.09, 0.07, 0.48,  1.0,  4.0,   2.7,  4.7,
     9              0.0006,0.09, 0.07, 0.48, 0.87,  2.2,   2.5,  4.3,
     *              0.001,0.001,  6.2,  7.0,  3.2, 12.0,  16.0, 13.0,
     1              0.001,0.001,  1.9, 16.0,  2.7, 13.0,  13.0,  8.0,
     2              0.001,0.001,  1.9, 16.0,  2.7, 13.0,  13.0,  8.0,
     3              0.001,0.001, 0.15,  0.1,  0.8,  8.0,   8.0,  5.0,
     4              0.001,0.001, 0.15,  0.1,  0.8,  8.0,   8.0,  5.0,
     5              0.001,0.001, 0.05, 0.22,  0.1, 0.16,   0.5,  0.9,
     6              0.001,0.001, 0.05, 0.22,  0.1, 0.16,   0.5,  0.9/

c     Since feff8 can be called any number of times . ALA

      if (ihole .le. 0)  then
         gamach = 0
         write(slog,'(a,1pe13.5)') ' No hole in SETGAM, gamach = ', 
     1                             gamach
         call wlog(slog)
         return
      endif
      if (ihole .gt. 16)  then
         call wlog(' This version of FEFF will set gamach = 0.1 eV ' //
     1             ' for O1 and higher hole')
         call wlog(' You can use CORRECTIONS card  to set ' //
     1   ' gamach = 0.1 + 2*vicorr ')
c        stop 'SETGAM-2'
      endif

      zz = iz
      if (ihole .le. 16)  then
         do 10  i = 1, 8
            gamkp(i) = log10 (gamh(i,ihole))
            zk(i) = zh(i,ihole)
   10    continue
         call terp (zk, gamkp, 8, 1, zz, gamach)
      else
c     include data from the tables later.
c     Now gamach=0.1eV for any O-hole for any element.
         gamach = -1.0
      endif

c     Change from log10 (gamma) to gamma
      gamach = 10.0 ** gamach


      return
      end
      subroutine iniptz(ptz,iptz,modus)
        !KJ This routine rewrites the ptz-matrix.

      implicit none
c     which polarization tensor to create      
      integer iptz
c     two ways of working
      integer modus
c     the polarization tensor
      complex*16 ptz(-1:1,-1:1)

      complex*16 zero,one,coni
      parameter (zero=(0,0),one=(1,0),coni=(0,1))
      integer i,j
      real*8 unity(3,3)



      if (iptz.lt.1.or.iptz.gt.10) then
          call wlog('Inieln sees weird iptz - returns without 
     1      changing ptz - danger of calculating nonsense !!')
      endif


      do i=1,3
      do j=1,3
	unity(i,j)=dble(0)
      enddo
      unity(i,i)=dble(1)/dble(3)
      enddo
      do i=-1,1
      do j=-1,1
        ptz(i,j)=zero
      enddo
      enddo
      

      if (modus.eq.1) then
! work in spherical coordinates

         if(iptz.eq.10) then
        do i=-1,1
	do j=-1,1
	   ptz(i,j)=unity(i+2,j+2)
        enddo
	enddo
         else
            i=(iptz-1)/3+1  ! row index
            j=iptz-3*(i-1)  ! column index
            i=i-2 !shift from 1..3 to -1..1
            j=j-2
            ptz(i,j)=one
         endif  


      elseif(modus.eq.2) then
! work in carthesian coordinates      

      if (iptz.eq.10) then ! orientation averaged spectrum
        do i=-1,1
	do j=-1,1
	   ptz(i,j)=unity(i+2,j+2)
        enddo
	enddo
	
        elseif (iptz.eq.1) then   ! x x*
          ptz(1,1)=one/dble(2)
        ptz(-1,-1)=one/dble(2)
        ptz(-1,1)=-one/dble(2)
        ptz(1,-1)=-one/dble(2)
        elseif (iptz.eq.5) then ! y y*
          ptz(1,1)=one/dble(2)
        ptz(-1,-1)=one/dble(2)
        ptz(-1,1)=one/dble(2)
        ptz(1,-1)=one/dble(2)
         elseif (iptz.eq.9) then ! z z*
        ptz(0,0)=one
        elseif (iptz.eq.2) then ! x y*
          ptz(1,1)=one*coni/dble(2)
        ptz(-1,-1)=-one*coni/dble(2)
        ptz(-1,1)=-one*coni/dble(2)
        ptz(1,-1)=one*coni/dble(2)
        elseif (iptz.eq.4) then ! x* y
          ptz(1,1)=-one*coni/dble(2)
        ptz(-1,-1)=one*coni/dble(2)
        ptz(-1,1)=-one*coni/dble(2)
        ptz(1,-1)=one*coni/dble(2)
        elseif (iptz.eq.3) then ! x z*
          ptz(-1,0)=one/dsqrt(dble(2))
        ptz(1,0)=-one/dsqrt(dble(2))
        elseif (iptz.eq.7) then ! x* z
          ptz(0,-1)=one/dsqrt(dble(2))
        ptz(0,1)=-one/dsqrt(dble(2))
        elseif (iptz.eq.6) then ! y z*
          ptz(-1,0)=-one*coni/dsqrt(dble(2))
        ptz(1,0)=-one*coni/dsqrt(dble(2))
        elseif (iptz.eq.8) then ! y* z
          ptz(0,-1)=one*coni/dsqrt(dble(2))
        ptz(0,1)=one*coni/dsqrt(dble(2))
      endif
      
      
        else
          stop 'alien modus in inieln'
        endif


      return
      end


C From HDK@psuvm.psu.edu Thu Dec  8 15:27:16 MST 1994
C 
C The following was converted from Algol recursive to Fortran iterative
C by a colleague at Penn State (a long time ago - Fortran 66, please
C excuse the GoTo's). The following code also corrects a bug in the
C Quicksort algorithm published in the ACM (see Algorithm 402, CACM,
C Sept. 1970, pp 563-567; also you younger folks who weren't born at
C that time might find interesting the history of the Quicksort
C algorithm beginning with the original published in CACM, July 1961,
C pp 321-322, Algorithm 64). Note that the following algorithm sorts
C integer data; actual data is not moved but sort is affected by sorting
C a companion index array (see leading comments). The data type being
C sorted can be changed by changing one line; see comments after
C declarations and subsequent one regarding comparisons(Fortran
C 77 takes care of character comparisons of course, so that comment
C is merely historical from the days when we had to write character
C compare subprograms, usually in assembler language for a specific
C mainframe platform at that time). But the following algorithm is
C good, still one of the best available.


      SUBROUTINE QSORTI (ORD,N,A)
C
C==============SORTS THE ARRAY A(I),I=1,2,...,N BY PUTTING THE
C   ASCENDING ORDER VECTOR IN ORD.  THAT IS ASCENDING ORDERED A
C   IS A(ORD(I)),I=1,2,...,N; DESCENDING ORDER A IS A(ORD(N-I+1)),
C   I=1,2,...,N .  THIS SORT RUNS IN TIME PROPORTIONAL TO N LOG N .
C
C
C     ACM QUICKSORT - ALGORITHM #402 - IMPLEMENTED IN FORTRAN 66 BY
C                                 WILLIAM H. VERITY, WHV@PSUVM.PSU.EDU
C                                 CENTER FOR ACADEMIC COMPUTING
C                                 THE PENNSYLVANIA STATE UNIVERSITY
C                                 UNIVERSITY PARK, PA.  16802
C
      IMPLICIT INTEGER (A-Z)
C
      DIMENSION ORD(N),POPLST(2,20)
      DOUBLE PRECISION X,XX,Z,ZZ,Y
C
C     TO SORT DIFFERENT INPUT TYPES, CHANGE THE FOLLOWING
C     SPECIFICATION STATEMENTS; FOR EXAMPLE, FOR FORTRAN CHARACTER
C     USE THE FOLLOWING:  CHARACTER *(*) A(N)
C
      DOUBLE PRECISION A(N)
C
      NDEEP=0
      U1=N
      L1=1
      DO 1  I=1,N
    1 ORD(I)=I
    2 IF (U1.LE.L1) THEN
         RETURN
      END IF
C
    3 L=L1
      U=U1
C
C PART
C
    4 P=L
      Q=U
C     FOR CHARACTER SORTS, THE FOLLOWING 3 STATEMENTS WOULD BECOME
C     X = ORD(P)
C     Z = ORD(Q)
C     IF (A(X) .LE. A(Z)) GO TO 2
C
C     WHERE "CLE" IS A LOGICAL FUNCTION WHICH RETURNS "TRUE" IF THE
C     FIRST ARGUMENT IS LESS THAN OR EQUAL TO THE SECOND, BASED ON "LEN"
C     CHARACTERS.
C
      X=A(ORD(P))
      Z=A(ORD(Q))
      IF (X.LE.Z) GO TO 5
      Y=X
      X=Z
      Z=Y
      YP=ORD(P)
      ORD(P)=ORD(Q)
      ORD(Q)=YP
    5 IF (U-L.LE.1) GO TO 15
      XX=X
      IX=P
      ZZ=Z
      IZ=Q
C
C LEFT
C
    6 P=P+1
      IF (P.GE.Q) GO TO 7
      X=A(ORD(P))
      IF (X.GE.XX) GO TO 8
      GO TO 6
    7 P=Q-1
      GO TO 13
C
C RIGHT
C
    8 Q=Q-1
      IF (Q.LE.P) GO TO 9
      Z=A(ORD(Q))
      IF (Z.LE.ZZ) GO TO 10
      GO TO 8
    9 Q=P
      P=P-1
      Z=X
      X=A(ORD(P))
C
C DIST
C
   10 IF (X.LE.Z) GO TO 11
      Y=X
      X=Z
      Z=Y
      IP=ORD(P)
      ORD(P)=ORD(Q)
      ORD(Q)=IP
   11 IF (X.LE.XX) GO TO 12
      XX=X
      IX=P
   12 IF (Z.GE.ZZ) GO TO 6
      ZZ=Z
      IZ=Q
      GO TO 6
C
C OUT
C
   13 CONTINUE
      IF (.NOT.(P.NE.IX.AND.X.NE.XX)) GO TO 14
      IP=ORD(P)
      ORD(P)=ORD(IX)
      ORD(IX)=IP
   14 CONTINUE
      IF (.NOT.(Q.NE.IZ.AND.Z.NE.ZZ)) GO TO 15
      IQ=ORD(Q)
      ORD(Q)=ORD(IZ)
      ORD(IZ)=IQ
   15 CONTINUE
      IF (U-Q.LE.P-L) GO TO 16
      L1=L
      U1=P-1
      L=Q+1
      GO TO 17
   16 U1=U
      L1=Q+1
      U=P-1
   17 CONTINUE
      IF (U1.LE.L1) GO TO 18
C
C START RECURSIVE CALL
C
      NDEEP=NDEEP+1
      POPLST(1,NDEEP)=U
      POPLST(2,NDEEP)=L
      GO TO 3
   18 IF (U.GT.L) GO TO 4
C
C POP BACK UP IN THE RECURSION LIST
C
      IF (NDEEP.EQ.0) GO TO 2
      U=POPLST(1,NDEEP)
      L=POPLST(2,NDEEP)
      NDEEP=NDEEP-1
      GO TO 18
C
C END SORT
C END QSORT
C
      END
c///////////////////////////////////////////////////////////////////////
c PAR Subroutines
c Written by J. Sims, NIST, 2001

c This software was developed at the National Institute of Standards
c and Technology by employees of the Federal Government in the course
c of their official duties. Pursuant to title 17 Section 105 of
c the United States Code this software is not subject to copyright
c protection and is in the public domain. PAR is an experimental
c system. NIST assumes no responsibility whatsoever for its use by
c other parties, and makes no guarantees, expressed or implied, about 
c its quality, reliability, or any other characteristic. We would
c appreciate acknowledgement if the software is used.

c This software can be redistributed and/or modified freely provided
c that any derivative works bear some notice that they are derived from
c it, and any modified versions bear some notice that they have been
c modified.
c///////////////////////////////////////////////////////////////////////
c License is applicable for routines below, until otherwise specified.
c
c  **************************************************
c  Parallel feff8 routines
c  Jim Sims
c  **************************************************

      subroutine par_begin
c  **************************************************
c  Initializations for parallel version(s)
c  **************************************************

c={../HEADERS/parallel.h
      integer par_type, this_process, numprocs, my_rank
      logical master, worker, parallel_run
      real*8 wall_comm, time_comm
      common /timing/ wall_comm, time_comm
      common /parallel/ numprocs, my_rank, this_process, 
     .          master, worker, parallel_run, par_type
c= ../HEADERS/parallel.h}
c={mpif.h
c={mpif.h
      include 'mpif.h'
c= not found}
c= not found}

c-- So cvd or dbx can attach to a running process
c     call sleep(30) 

      call MPI_INIT(ierrorflag)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierrorflag)
      call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierrorflag)
      this_process = my_rank

      par_type = 0
      parallel_run = .true.
c-- The following variable will be used for IO that should only be
c-- done in one process.
      master = (my_rank .eq. 0)

      worker = (.not. master)
      if (worker) par_type = 1

c     write(6,*) 'this process = ',this_process, ' worker = ',worker

      if (master) write(6,*) 'Number of processors = ',numprocs

      return
      end

      subroutine par_stop (string)
c  **************************************************
c  Abnormal termination of the parallel session
c  **************************************************
c={../HEADERS/parallel.h
      integer par_type, this_process, numprocs, my_rank
      logical master, worker, parallel_run
      real*8 wall_comm, time_comm
      common /timing/ wall_comm, time_comm
      common /parallel/ numprocs, my_rank, this_process, 
     .          master, worker, parallel_run, par_type
c= ../HEADERS/parallel.h}
c={mpif.h
c={mpif.h
      include 'mpif.h'
c= not found}
c= not found}
c     For abnormal exits 
c     If open, close unit = 11
c     Go to the barrier that workers are sitting at
c     Then everyone will call par_end and stop
      logical is_open
      character*(*) string

      inquire(unit=11,opened=is_open)
      if (is_open) then
        call wlog(string)
        close(unit=11)
      else if (string .ne. ' ') then
	print *,string
	print *,'Abnormal termination on processor ',this_process
      endif
      call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierrorflag)

      stop ' '
      end

      subroutine par_end
c  **************************************************
c  Terminate the parallel session
c  **************************************************
      call MPI_FINALIZE(ierrorflag)
      return
      end

      subroutine par_barrier
c  **************************************************
c  Calls mpi_barrier
c  **************************************************
c={mpif.h
c={mpif.h
      include 'mpif.h'
c= not found}
c= not found}
      call MPI_BARRIER(MPI_COMM_WORLD,ierrorflag)
      return
      end

      subroutine par_send_int(buf,count,dest,tag)
c  **************************************************
c  Call mpi_send for integer arrays
c  **************************************************
      integer count,dest,tag
      integer buf(*)
c={mpif.h
c={mpif.h
      include 'mpif.h'
c= not found}
c= not found}
      call MPI_SEND(buf,count,MPI_INTEGER,dest,tag,
     .              MPI_COMM_WORLD,ierrorflag)
      return
      end

      subroutine par_send_cmplx(buf,count,dest,tag)
c  **************************************************
c  Call mpi_send for complex arrays
c  **************************************************
      integer count,dest,tag
      complex buf(*)
c={mpif.h
c={mpif.h
      include 'mpif.h'
c= not found}
c= not found}
      call MPI_SEND(buf,count,MPI_COMPLEX,dest,tag,
     .              MPI_COMM_WORLD,ierrorflag)
      return
      end

      subroutine par_send_dc(buf,count,dest,tag)
c  **************************************************
c  Call mpi_send for double_complex arrays
c  **************************************************
      integer count,dest,tag
      complex*16 buf(*)
c={mpif.h
c={mpif.h
      include 'mpif.h'
c= not found}
c= not found}
      call MPI_SEND(buf,count,MPI_DOUBLE_COMPLEX,dest,tag,
     .              MPI_COMM_WORLD,ierrorflag)
      return
      end

      subroutine par_recv_int(buf,count,source,tag)
c  **************************************************
c  Call mpi_recv for integer arrays
c  **************************************************
      integer count,source,tag
      integer buf(*)
c={mpif.h
c={mpif.h
      include 'mpif.h'
c= not found}
c= not found}
      integer istat(mpi_status_size)
      call MPI_RECV(buf,count,MPI_INTEGER,source,tag,
     .              MPI_COMM_WORLD,istat,ierrorflag)
      return
      end

      subroutine par_recv_cmplx(buf,count,source,tag)
c  **************************************************
c  Call mpi_recv for complex arrays
c  **************************************************
      integer count,source,tag
      complex buf(*)
c={mpif.h
c={mpif.h
      include 'mpif.h'
c= not found}
c= not found}
      integer istat(mpi_status_size)
      call MPI_RECV(buf,count,MPI_COMPLEX,source,tag,
     .              MPI_COMM_WORLD,istat,ierrorflag)
      return
      end

      subroutine par_recv_dc(buf,count,source,tag)
c  **************************************************
c  Call mpi_recv for double complex arrays
c  **************************************************
      integer count,source,tag
      complex*16 buf(*)
c={mpif.h
c={mpif.h
      include 'mpif.h'
c= not found}
c= not found}
      integer istat(mpi_status_size)
      call MPI_RECV(buf,count,MPI_DOUBLE_COMPLEX,source,tag,
     .              MPI_COMM_WORLD,istat,ierrorflag)
      return
      end

      subroutine par_bcast_int(buf,count,source)
c  **************************************************
c  Call mpi_bcast for integer arrays
c  **************************************************
      integer count,source
      integer buf(*)
c={mpif.h
c={mpif.h
      include 'mpif.h'
c= not found}
c= not found}
      call MPI_BCAST(buf,count,MPI_INTEGER,source,
     .              MPI_COMM_WORLD,ierrorflag)
      return
      end

      subroutine par_bcast_cmplx(buf,count,source)
c  **************************************************
c  Call mpi_bcast for complex arrays
c  **************************************************
      integer count,source
      complex buf(*)
c={mpif.h
c={mpif.h
      include 'mpif.h'
c= not found}
c= not found}
      call MPI_BCAST(buf,count,MPI_COMPLEX,source,
     .              MPI_COMM_WORLD,ierrorflag)
      return
      end

      subroutine par_bcast_dc(buf,count,source)
c  **************************************************
c  Call mpi_bcast for double_complex arrays
c  **************************************************
      integer count,source
      complex*16 buf(*)
c={mpif.h
c={mpif.h
      include 'mpif.h'
c= not found}
c= not found}
      call MPI_BCAST(buf,count,MPI_DOUBLE_COMPLEX,source,
     .              MPI_COMM_WORLD,ierrorflag)
      return
      end

      subroutine MPE_DECOMP1D( n, num_procs, myid, s, e )
c  ******************************************************
c  A routine for producing a decomposition of a 1-d 
c  array when given a number of processors.  It may 
c  be used in "direct" product decomposition.  The 
c  values returned assume a "global" domain in [1:n]
c  ******************************************************
c  MPE_Decomp1d - Compute a balanced decomposition of
c  a 1-D array
c  ******************************************************
c  Input Parameters:
c  n  - Length of the array
c  num_procs - Number of processors in decomposition
c  myid  - Rank of this processor in the decomposition 
c  (0 <= rank < size)
c  ******************************************************
c  Output Parameters:
c  s,e - Array my_particles are s:e, with the original 
c  array considered as 1:n.  
c  ******************************************************

      integer n, num_procs, myid, s, e
      integer nloc, deficit
 
      nloc  = n / num_procs
      s       = myid * nloc + 1
      deficit = mod(n,num_procs)
      s       = s + min(myid,deficit)
      if (myid .lt. deficit) then
        nloc = nloc + 1
      endif
      e = s + nloc - 1
      if (e .gt. n .or. myid .eq. num_procs-1) e = n

      return
      end

      SUBROUTINE SECONDS( W)
c  ***************************************************
c  SECONDS returns the wall clock times for a process
c  in seconds.
c  ***************************************************

      real*8 W, MPI_Wtime

      W = MPI_Wtime()

      RETURN
      END
