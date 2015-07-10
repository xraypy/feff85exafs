      subroutine pot (verbse, rgrd, nohole,
     $       inters, totvol, ecv0, nscmt, nmix, ntitle, title,
     $       nat, nph, ihole, iafolp, ixc, iphat, rat, iatph, xnatph,
     $       novr, iphovr, nnovr, rovr, folp0, xion, iunf, iz, ipr1,
     $       ispec, jumprm, lmaxsc, icoul, ca1, rfms1, lfms1,
     $
     -       rnrmav, xmu, vint, rhoint,
     1       emu, s02, erelax, wp, rs, xf, qtotel,
     2       imt, rmt, inrm, rnrm, folpx,
     3       dgc0, dpc0, dgc, dpc, adgc, adpc,
     4       edens, vclap, vtot, edenvl, vvalgs, dmag, xnval,
     5       eorb, kappa, iorb, qnrm, xnmues, nhtmp
     6       )
c                                                                   !     return stuff for use in wrpot & xsph
c gamach

c##****f* feff85exafs/pot
c##  NAME
c##     pot
c##
c##  SYNOPSIS
c##     pot(rgrd, nohole,
c##    $    inters, totvol, ecv0, nscmt, nmix, ntitle, title,
c##    $    nat, nph, ihole, iafolp, ixc, iphat, rat, iatph, xnatph,
c##    $    novr, iphovr, nnovr, rovr, folp0, xion, iunf, iz, ipr1,
c##    $    ispec, jumprm, lmaxsc, icoul, ca1, rfms1, lfms1,
c##    $    rnrmav, xmu, vint, rhoint,
c##    $    emu, s02, erelax, wp, rs, xf, qtotel,
c##    $    imt, rmt, inrm, rnrm, folpx,
c##    $    dgc0, dpc0, dgc, dpc, adgc, adpc,
c##    $    edens, vclap, vtot, edenvl, vvalgs, dmag, xnval,
c##    $    eorb, kappa, iorb, qnrm, xnmues, nhtmp)
c##
c##  FUNCTION
c##     Compute potentials for an input atomic cluster, returning data needed to compute phase shifts
c##
c##
c##  INPUTS
c##     verbse      boolean flag, true=write screen messages
c##     rgrd        radial grid spacing (RGRD)
c##     nohole      flag to compute with no hole on absorber (NOHOLE)
c##     inters      parameter of the interstitial calculation (INTERSTITIAL)
c##     totvol      parameter of the interstitial calculation (INTERSTITIAL)
c##     ecv0        core/valence separation (EXCHANGE)
c##     nscmt       max number of self-conmsistentcy iterations (SCF)
c##     nmix        number of iterations before switching to Broyden (SCF)
c##     ntitle      number of title lines (TITLE)
c##     title       title lines (TITLE)
c##     nat         number of atoms in cluster (ATOMS)
c##     nph         number of unique potentials in cluster (POTENTIALS)
c##     ihole       edge index, 1=K, 4=L3, etc (EDGE/HOLE)
c##     iafolp      flag indicating automatic overlapping (AFOLP)
c##     ixc         exchange index (EXCHANGE)
c##     iphat       unique potential indeces of atoms in cluster (ATOMS)
c##     rat         cartesian coordinates of atoms in cluster (ATOMS)
c##     iatph       example of each unique potential in cluster
c##     xnatph      stoichiometries of each potential index (POTENTIALS)
c##     novr        number of overlap shells for unique potential (obsolete in feff85exafs)
c##     iphovr      unique potential for this overlap shell (obsolete in feff85exafs)
c##     nnovr       number of atoms in overlap shell (obsolete in feff85exafs)
c##     rovr        r for overlap shell (obsolete in feff85exafs)
c##     folp0       overlap factor for rmt calculation (FOLP/AFOLP)
c##     xion        ioniziation of each potential (ION)
c##     iunf        flag to unfreeze f eectrons (UNFREEZEF)
c##     iz          Z numbers of each potential (POTENTIALS)
c##     ipr1        print flag (not used in library)
c##     ispec       spectroscopy index, 1=EXAFS (always set to 1 in feff85exafs)
c##     jumprm      flag to remove jumps at muffin tin radii (JUMPRM)
c##     lmaxsc      l max for SCF for each potential (POTENTIALS)
c##     icoul       obsolete param. for handling Coulomb potential (SCF)
c##     ca1         self-consistency convergence accelerator (SCF)
c##     rfms1       cluster radius for self-consistent calculation (SCF)
c##     lfms1       0=solid, 1=molecule (SCF)
c##     rnrmav      average Norman radius in cluster
c##     xmu         Fermi level in hartrees
c##     vint        interstitial energy
c##     rhoint      interstitial density * 4 * pi
c##
c##
c##  RESULT
c##     emu         ionization energy in adiabatic approximation
c##     s02         computed S02
c##     erelax      ionization energy - emu
c##     wp          plasmon frequency in hartrees
c##     rs          interstitial density parameter (see istprm.f and fermi.f)
c##     xf          interstital fermi momentum
c##     qtotel      total number of e in a cluster
c##     imt         r mesh index just inside rmt
c##     rmt         muffin tin radii
c##     inrm        r mesh index just inside rnorman
c##     rnrm        Norman radii
c##     folpx       overlap factor for rmt calculation
c##     dgc0        additional data needed for relativistic version
c##     dpc0        additional data needed for relativistic version
c##     dgc         additional data needed for relativistic version
c##     dpc         additional data needed for relativistic version
c##     adgc        additional data needed for relativistic version
c##     adpc        additional data needed for relativistic version
c##     edens       core density
c##     vclap       overlapped Coulomb potential
c##     vtot        overlapped Coulomb potential
c##     edenvl      additional data needed for relativistic version
c##     vvalgs      additional data needed for relativistic version
c##     dmag        data for single configuration Dirac-Fock atom code
c##     xnval       data for single configuration Dirac-Fock atom code
c##     eorb        data for single configuration Dirac-Fock atom code
c##     kappa       data for single configuration Dirac-Fock atom code
c##     iorb        data for single configuration Dirac-Fock atom code
c##     qnrm        Norman radii after ovrlp.f
c##     xnmues      
c##     nhtmp       holds nohole value
c##    
c##
c##  NOTES
c##    Call wrpot to write the pot.pad file.  Call reapot to read the potph.json (potph.inp) file.
c##    See ffmod1.f for the use of this subroutine in the conventional stand-alone program.
c##
c##  BUGS
c##    Report bugs and other issues at https://github.com/xraypy/feff85exafs
c##
c##  LICENSE
c##    See src/HEADERS/license.h for the terms of the parts of feff85exafs derived directly from Feff
c##
c##    nxjson.c and nxjson.h are Copyright (c) 2013 Yaroslav Stavnichiy <yarosla@gmail.com>.  See
c##    https://bitbucket.org/yarosla/nxjson/src
c##
c##    The C wrapper around this subroutine is released to the public domain
c##
c##  SEE ALSO
c##    The pot and xsph stand-alone programs, the libfeffphases.c wrapper
c##
c##****



c     Cluster code -- multiple shell single scattering version of FEFF
c     This program (or subroutine) calculates potentials and phase
c     shifts for unique potentials specifed by atoms and overlap cards.
c
c     Input files:  potph.inp    input data, atoms, overlaps, etc.
c     Output:       phases.bin   phase shifts for use by the rest of the
c                                program
c                   xxx.dat      various diagnostics

      implicit double precision (a-h, o-z)
      include '../HEADERS/const.h'
      include '../HEADERS/parallel.h'
      include '../HEADERS/dim.h'
      Parameter (Maxprocs = 1)

      logical verbse

c     Notes:
c        nat    number of atoms in problem
c        nph    number of unique potentials
c        ihole  hole code of absorbing atom
c        iph=0 for central atom

c     Specific atom input data
c     iphat - given specific atom, which unique pot?
      dimension iphat(natx)
c     rat(3,natx)  -  cartesian coords of specific atom
      dimension rat(3,natx)
      real rfms1

c     Unique potential input data
c     iatph(0:nphx) - given unique pot, which atom is model?
c                   (0 if none specified for this unique pot)
      dimension iatph(0:nphx)
c     xnatph(0:nphx) - given unique pot, how many atoms are there
c                      of this type? (used for interstitial calc)
      dimension xnatph(0:nphx)

c     folp(0:nphx)  - overlap factor for rmt calculation
      dimension folp(0:nphx), folp0(0:nphx), folpx(0:nphx)
c     novr(0:nphx)  - number of overlap shells for unique pot
      dimension novr(0:nphx)
c     iphovr(novrx,0:nphx)  - unique pot for this overlap shell
      dimension iphovr(novrx,0:nphx)
c     nnovr(novrx,0:nphx)  -  number of atoms in overlap shell
      dimension nnovr(novrx,0:nphx)
c     rovr(novrx,0:nphx)   -  r for overlap shell
      dimension rovr(novrx,0:nphx)

c     Free atom data
c     xion(0:nphx)  - ionicity, input
      dimension xion(0:nphx)
c     iz(0:nphx)  -   atomic number, input
      dimension iz(0:nphx)

c     ATOM output
c     Note that ATOM output is dimensioned 251, all other r grid
c     data is set to nrptx, currently 250
c     rho(251,0:nphx+1)     -   density*4*pi
      dimension rho(251,0:nphx+1)
c     vcoul(251,0:nphx+1)   -   coulomb potential
      dimension vcoul(251,0:nphx+1)
      dimension dr(251), drho(251), dvcoul(251)

c     Overlap calculation results
c     overlapped density*4*pi
      dimension edens(251,0:nphx)
c     overlapped coul pot
      dimension vclap(251,0:nphx), vclapp(251,0:nphx)
c     overlapped total potential
      dimension vtot (251,0:nphx)

c     Muffin tin calculation results
c     r mesh index just inside rmt
      dimension imt(0:nphx)
c     r mesh index just inside rnorman
      dimension inrm(0:nphx)
c     muffin tin radius
      dimension rmt(0:nphx)
c     norman radius
      dimension rnrm(0:nphx), qnrm(0:nphx), qold(0:nphx), lmaxsc(0:nphx)
      dimension xnmues(0:lx,0:nphx)
      character*80 title(nheadx)

      logical ok

      complex gtr((lx + 1) * (nphx + 1) * Maxprocs)
      complex*16 xrhoce((lx + 1) * (nphx + 1) * Maxprocs)
      complex*16 xrhole((lx + 1) * (nphx + 1) * Maxprocs)
      complex*16 yrhoce(251 * (nphx + 1) * Maxprocs )
      complex*16 yrhole(251 * (lx + 1) * (nphx + 1) * Maxprocs )
c     need irregular solution for complex potential. fix later
      dimension dgc0(251), dpc0(251)

c     additioal data needed for relativistic version
      dimension dgc(251,30,0:nphx+1), dpc(251,30,0:nphx+1)
      dimension adgc(10,30,0:nphx+1), adpc(10,30,0:nphx+1)
      dimension rhoval(251,0:nphx+1), edenvl(251,0:nphx)
      dimension vvalgs (251,0:nphx)

c     nrx = max number of r points for phase and xsect r grid
c      parameter (nrx = nrptx)
      dimension ri(nrptx)
      dimension dmag(251,0:nphx+1), xnvmu(0:lx,0:nphx+1)
      dimension xnval(30,0:nphx+1), norb(0:nphx+1), eorb(30,0:nphx+1)
      dimension kappa(30,0:nphx+1), iorb(-4:3,0:nphx+1)
c       criteria for self-consistency
      parameter (tolq = 1.D-3)
      parameter (tolmu = 3.D-3)
      logical lpass
      character*512 slog
c     Josh use nhtmp to save nohole value
      integer nhtmp

      do 4 i=1,30
         do 2 j=0,nphx+1
            kappa(i,j) = 0
            eorb(i,j) = 0
 2       continue
 4    continue

   10 format (4x, a, i5)

c     Josh - for now if nohole=2 reset to 0 so that regular nohole
c     potential is used
      nhtmp = nohole
      if (nohole.eq.2) nohole = 0
c     Josh

c     variables ecv0 and folp0 serve as input only; do not change them
c     since it will change file feff.ior content
c     ecv and folp are passed through pot.pad to next modules.
      ecv = ecv0
      do 12 i = 0, nph
   12 folp(i) = folp0(i)

      call inipot (dgc, dpc, edenvl, vvalgs, xnmues)

c     increase the length of hydrogen bonds for potential only
      call moveh (nat, iphat, iz, rat)

      nfree = 1
      do 17 i=0,nph
        if (abs(xion(i)) .gt. 1.d-3) nfree = 2
  17  continue

c     Free atom potentials and densities
c     Final state is (usually) with a core hole, initial state is 
c     w/o a corehole.
c     NB wsatom is needed in SUMAX, if changed here, change it there
c     wsatom = 15
c     do not save spinors
c     Call twice if any of xion.neq.0 ( first time with xion=0 to set
c     rnrm)

      do 99 ifree = 1, nfree

      ispinr = 0
      etfin  = 0
      do 20  iph = 0, nph
         if (verbse) then
            write(slog,10) 
     1             'free atom potential and density for atom type', iph
            call wlog(slog)
         endif
c        Include corehole if absorber (unless user says nohole)
         if (iph .eq. 0)  then
            itmp = ihole
         else
            itmp = 0
         endif
         if (nohole.ge.0 .and. iph.eq.0) then
           xionp = xion(0)
           if (nfree.eq.2 .and. ifree.eq.1) xionp = 0
           call scfdat ( ipr1, nph+1, nph, iz(0), itmp, xionp, iunf,
     1     vcoul(1,nph+1), rho(1,nph+1), dmag(1,nph+1), rhoval(1,nph+1),
     2     ispinr, dgc0, dpc0, dgc, dpc, adgc, adpc, 
     3     s02, efrozn, et, xnvmu(0,nph+1),
     4     xnval(1,nph+1), iorb(-4,nph+1), norb(nph+1),
     5     eorb(1,nph+1), kappa(1,nph+1) )
         else
           xionp = xion(iph)
           if (nfree.eq.2 .and. ifree.eq.1) xionp = 0
           call scfdat ( ipr1, iph, nph, iz(iph), itmp, xionp, iunf,
     1         vcoul(1,iph), rho(1,iph), dmag(1,iph), rhoval(1,iph),
     2         ispinr, dgc0, dpc0, dgc, dpc, adgc, adpc, 
     3         s02, efrozn, et, xnvmu(0,iph),
     4         xnval(1,iph), iorb(-4,iph), norb(iph),
     5         eorb(1,iph), kappa(1,iph) )
         endif
c        etfin is absorbing atom final state total energy, see nohole
c           case below
         if (iph .eq. 0) etfin = et
   20 continue

      if (verbse) then
         write(slog,10) 'initial state energy'
         call wlog(slog)
      endif
c     Save initial state energy and spinors for core hole orbital,
c     do not save potentials, except for nohole.
      ispinr = ihole
      itmp = 0
      if (nohole.ge.0) then
         iph = 0
         xionp = xion(iph)
         if (nfree.eq.2 .and. ifree.eq.1) xionp = 0
         call scfdat ( ipr1, iph, nph, iz(iph), itmp, xionp, iunf,
     1         vcoul(1,iph), rho(1,iph), dmag(1,iph), rhoval(1,iph),
     2         ispinr, dgc0, dpc0, dgc, dpc, adgc, adpc, 
     3         s02, efrozn, etinit, xnvmu(0,iph),
     4         xnval(1,iph), iorb(-4,iph), norb(iph),
     5         eorb(1,iph), kappa(1,iph) )
      else
         xionp = xion(0)
         if (nfree.eq.2 .and. ifree.eq.1) xionp = 0
         call scfdat ( ipr1, nph+1, nph, iz(0), itmp, xionp, iunf,
     1     vcoul(1,nph+1), rho(1,nph+1), dmag(1,nph+1), rhoval(1,nph+1),
     2     ispinr, dgc0, dpc0, dgc, dpc, adgc, adpc, 
     3     s02, efrozn, etinit, xnvmu(0,nph+1),
     4     xnval(1,nph+1), iorb(-4,nph+1), norb(nph+1),
     5     eorb(1,nph+1), kappa(1,nph+1) )
      endif

c     testing new potential for the final state. ala
      hx = 0.05d0
      x0 = -8.8d0
      if (nohole.gt.0) then
         idim = 251
         do 30 i = 1,idim
  30     dr(i) = exp(x0+hx*(i-1))
         if (nohole.eq.1) then
            do 40 i = 1,idim
  40        drho(i) = dgc0(i)**2 + dpc0(i)**2
         else
            do 50 i = 1,idim
               drho(i)=dr(i)**2 *
     1         (rho(i,0)-rhoval(i,0)-rho(i,nph+1)+rhoval(i,nph+1))
  50        continue
         endif
         call potslw ( dvcoul, drho, dr, hx,idim)
         do 60 i=1,idim
c           drho(i) = drho(i)/ dr(i)**2
c           use 1/2 of core-hole as in transition state
            drho(i) = drho(i)/2.0d0/ dr(i)**2
  60     continue
      else
         do 70 i=1,251
            drho(i) = 0
            dvcoul(i) = 0
  70     continue
      endif

c     etinit is absorbing atom initial state (no hole)
c     efrozn is ionization energy with frozen orbitals (koopman's
c      theorem)
c     etfin-etinit is ionization energy in adiabatic approximation
      erelax= -efrozn - ( etfin - etinit)
      emu = etfin - etinit

c     Overlap potentials and densitites
      do 90  iph = 0, nph
         if (verbse) then
            write(slog,10)
     1      'overlapped potential and density for unique potential', iph
            call wlog(slog)
         endif
         call ovrlp (iph, iphat, rat, iatph, novr, iphovr,
     1               nnovr, rovr, iz, nat, rho, dmag,
     2               rhoval, vcoul, edens, edenvl, vclap, qnrm)
         if (iph.eq.0) emu = emu - vclap(1,0)+vcoul(1,0)
   90 continue

      if (ifree.eq.1) then
c       Set the Norman radii 
        do 92 iph =0, nph
   92   rnrm(iph) = qnrm(iph)
      endif

   99 continue
c  end of free atom calculations (might be done twice if ION used)

cc new patch
c     itest = 1
c     if (itest.eq.1) then
cc      use orbitals with core-hole for initial orbitals
cc      orthogonaliztion problem for NRIXS calculations
c       do i = 1, 251
c         dgc0(i) = dgc(i,0)
c         dpc0(i) = dpc(i,0)
c       enddo
c     endif
cc end new patch
     
c     Find total charges for istprm
c     qtotel - total number of e in a cluster
      qtotel = 0
      do 80 iph = 0,nph
         qtotel = qtotel + (iz(iph)-xion(iph)) * xnatph(iph)
  80  continue
c     photoelectron moves out of the system
c     do not remove now since we are putting screening electron back


c     Find muffin tin radii, add gsxc to potentials, and find
c     interstitial parameters
      if (verbse) then
         write(slog,10) 'muffin tin radii and interstitial parameters'
         call wlog(slog)
      endif

      rmt(0) = -1
      xmu = 100.d0
      if (iafolp.ge.0) then
        do 101 iph=0,nph
          folpx(iph) = folp(iph)
          folp(iph) = 1
  101   continue
      endif
        
      idmag = 0
      call istprm (nph, nat, iphat, rat, iatph, xnatph,
     1            novr, iphovr, nnovr, rovr, folp, folpx, iafolp,
     1            edens, edenvl, idmag,
     2            dmag, vclap, vtot, vvalgs, imt, inrm, rmt, rnrm,
     3            ixc, rhoint,vint, rs, xf, xmu, xmunew,
     4            rnrmav, qtotel, inters, totvol)
      xmu = xmunew

c     Automatic max reasonable overlap
      if (iafolp .ge. 0)  then
         call afolp (verbse, nph, nat, iphat, rat, iatph, xnatph,
     1               novr, iphovr, nnovr, rovr, folp, folpx, iafolp,
     1               edens, edenvl,
     2               dmag, vclap, vtot, vvalgs, imt, inrm, rmt, rnrm,
     3               ixc, rhoint,vint, rs, xf, xmu, xmunew,
     4               rnrmav, qtotel, inters, totvol)
         xmu =xmunew
      endif

c     wp is plasmon frequency in hart
      wp = sqrt(12.d0*rs/fa**4) * xf**2 / 2.d0

c     Phase shift calculation
c     Atom r grid
      dx = 0.05d0
      x0 = 8.8d0

c     Find self-consistent muffin-tin potential.
      do 105 iph=0,nph
         qnrm(iph) = 0
         qold(iph) = 0
  105 continue

  100 continue
      if (nscmt.gt.0 .or. (ispec.ne.0 .and. ispec.lt.4)) call corval
     1                 (verbse, ecv, xnvmu, eorb, norb, xnval,
     1                  kappa, rgrd, nohole,
     2                  nph, edens, edenvl, vtot, vvalgs,
     3                  rmt, rnrm, ixc, rhoint, vint, jumprm,
     4                  x0, ri, dx, xion, iunf, iz,
     5                  adgc, adpc, dgc, dpc, ihole, lmaxsc)


c     find a total number of valence electrons
c     xntot - required number of valence electrons below fermi level
c     xnvmu(iph) = xnvmu(iph)-xion(iph)
c     xnvmu - number of valence electron within norman sphere
      xntot=0.0d0
      do 120 iph=0,nph
         xnvmup = 0
         do 110  i = 0,lx
  110    xnvmup = xnvmup + xnvmu(i,iph)
c x35 and earlier   xntot = xntot + xnatph(iph)*(xnvmup+xion(iph))
         xntot = xntot + xnatph(iph) * xnvmup
  120 continue

c     need to update vxcval in case if the core-valence separation was
c     made in subroutine corval. Need vxcval only for nonlocal exchange.
      if (mod(ixc,10).ge.5) then
         call  istprm (nph, nat, iphat, rat, iatph, xnatph,
     1               novr, iphovr, nnovr, rovr, folp, folpx, iafolp,
     1               edens, edenvl, idmag,
     2               dmag, vclap, vtot, vvalgs, imt, inrm, rmt, rnrm,
     3               ixc, rhoint,vint, rs, xf, xmu, xmunew,
     5               rnrmav, qtotel, inters, totvol)
         xmunew = xmu
      endif

      if (verbse) then
         write(slog,130) xmu*hart
 130     format('    : mu_old= ',f9.3)
         call wlog(slog)
      endif

c     do first nmix iterations with mixing scheme. Need for f-elements.
  140 nmix=nmix-1

c     number of processors for parallel execution
      npr = numprocs
      do 200 iscmt =1,nscmt
c        need to store coulomb potential
         do 145 ip=0,nph
         do 145 ir=1,251
  145    vclapp(ir,ip) = vclap(ir,ip)

         if (npr.le.1) then
           call scmt (verbse, iscmt, ecv, nph, nat, vclap, edens,
     1                edenvl, vtot, vvalgs, rmt, rnrm, qnrm,
     2                ixc, rhoint, vint, xmunew, jumprm,
     3                xntot, xnvmu, xnval,
     4                x0, ri, dx, xnatph, xion, iunf, iz,
     5                adgc, adpc, dgc,dpc, ihole,
     7                rat, iatph, iphat, lmaxsc, rhoval, xnmues, ok,
     8                rgrd, nohole, nscmt, icoul, ca1, rfms1, lfms1)
         else
           call scmtmp(verbse, npr, iscmt, ecv, nph, nat, vclap, edens,
     1                edenvl, vtot, vvalgs, rmt, rnrm, qnrm,
     2                ixc, rhoint, vint, xmunew, jumprm,
     3                xntot, xnvmu, xnval,
     4                x0, ri, dx, xnatph, xion, iunf, iz,
     5                adgc, adpc, dgc,dpc, ihole,
     7                rat, iatph, iphat, lmaxsc, rhoval, xnmues, ok,
     8                rgrd, nohole, nscmt, icoul, ca1, rfms1, lfms1,
     9                gtr, xrhole, xrhoce, yrhole, yrhoce )
         endif

         if (.not. ok) goto 100
c        if need to change core-valence separation then
c        start scmt loop all over again

c        write out Fermi level and charge transfers 
c        and do tests of self-consistency
         lpass = .true.
         if (iscmt.lt.nscmt .and. iscmt.le.3) lpass =.false.
         if (verbse) then
            write (slog,150)   xmunew*hart
 150        format (' mu_new= ', f9.3)
            call wlog(slog)
         endif
         if (abs (xmunew - xmu) .gt. tolmu) lpass = .false.
         xmu = xmunew
c        print out charge 
         if (verbse) call wlog(' Charge transfer:  iph  charge(iph) ')
         do 170 iph=0,nph
            if (verbse) then
               write (slog,180) iph, -qnrm(iph) + xion(iph)
               call wlog(slog)
            endif
            if (abs(qnrm(iph)-qold(iph)).gt.tolq) lpass = .false.
            qold(iph) = qnrm(iph)

c           check self-consistency of charges
            sum = -qnrm(iph)
            do 160 il=0,lx
  160       sum = sum + xnmues(il,iph) - xnvmu(il,iph)
            if (abs(sum).gt.0.05d0) lpass = .false.
  170    continue
  180    format('     ',i3, 2f9.3)

c        recalculate core density (edens) here. fix later. ala
c        call scfdat
c        for now use the old core density
         if (iscmt.eq.nscmt .or. lpass) then
c           restore  total density from previous iteration
            do 190 ip=0,nph
              do 185 ir=1,251
c                need total density for istprm
                 edens(ir,ip) = edens(ir,ip)-rhoval(ir,ip)+edenvl(ir,ip)
                 vclap(ir,ip) = vclapp(ir,ip)
  185         continue
c             remember the reported charge transfer
              qnrm(ip) = -qnrm(ip) + xion(ip)
  190       continue
c           exit self-consistency loop
            goto 210
         else
c           update valence density
            do 195 ip=0,nph
            do 195 ir=1,251
c              need total density for istprm
               edenvl(ir,ip) = rhoval(ir,ip)
  195       continue
         endif

         call  istprm (nph, nat, iphat, rat, iatph, xnatph,
     1               novr, iphovr, nnovr, rovr, folp, folpx, iafolp,
     1               edens, edenvl, idmag,
     2               dmag, vclap, vtot, vvalgs, imt, inrm, rmt, rnrm,
     3               ixc, rhoint,vint, rs, xf, xmu, xmunew,
     5               rnrmav, qtotel, inters, totvol)
         xmunew = xmu
         if (nmix.gt.0) goto 140

  200 continue
c     suspicious exit: run out of iterations (iscmt=nscmt)

c     right exit from the loop: self-consistency is achieved
  210 continue
      if (worker) go to 400

      if (nohole.gt.0) then
c        testing new final state potential
         do 220 j = 1,251
  220    edens(j,0) = edens(j,0) - drho(j)
         
c        notice that vclap is actually for the next iteration
c        in SCMT loop, thus vclap may be wrong if self-consistency
c        has not been reached
         do 230 j = 1,251
  230    vclap(j,0) = vclap(j,0) - dvcoul(j)

         call  istprm (nph, nat, iphat, rat, iatph, xnatph,
     1      novr, iphovr, nnovr, rovr, folp, folpx, iafolp,
     1      edens, edenvl, idmag,
     2      dmag, vclap, vtot, vvalgs, imt, inrm, rmt, rnrm,
     3      ixc, rhoint,vint, rs, xf, xmu, xmunew,
     5      rnrmav, qtotel, inters, totvol)
      endif

c    correct the excitation energy
c     emu = emu -vclap(1,0) + vcoul(1,0) done also above
c     emu = emu+xmu  should be done in principle but leads
c     to worse estimate of edge position. fix later. ala

      if (ipr1 .ge. 2)  then
         call wpot (nph, edens, imt, inrm,
     1              rho, vclap, vcoul, vtot, ntitle, title)
      endif


  400 call par_barrier

      return
      end
