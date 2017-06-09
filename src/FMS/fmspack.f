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
      include '../HEADERS/const.h'
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
      complex   prefac, gllmz, ck(nspx)
c      complex term
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
c      character*3  cerr
      character*3  dec
c      character*13 trans
      character*75 messg

c 400  format(i4)

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
      call getkts(nsp, inclus, iphx, lipotx, i0)

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
      subroutine getkts(nsp, nat, iphx, lipotx, i0)

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
