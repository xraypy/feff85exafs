c     Josh - argument iPl has been added to arguments of xsect
      subroutine xsect (ipr2, dx, x0, ri, ne, ne1, ik0, em, edge,
     1                  ihole, emu, corr, dgc0, dpc0, jnew,
     2                  ixc, lreal, rmt, rnrm, xmu,
     2                  iPl,
     3                  vtot, vvalgs, edens, dmag, edenvl,
     4                  dgcn, dpcn, adgc, adpc, xsec, xsnorm, rkk,
     5                  iz, xion, iunf, xnval,
     5                  izstd, iorb, l2lp,
     6                  ipol, ispin, le2, angks, ptz)

c gamach, eorb, ifxc, kappa

c     right know the same self-energy is used for calculation
c     of the central atom part (xsec) and dipole m.e. for
c     scattering (rkk). You may want to run xsect separately
c     for xsec and for rkk, if you want to use different self-energy
c     for central and scattering parts.  ala. fix later

      implicit double precision (a-h, o-z)

c     INPUT
c     dx, x0, ri(nr)
c                  Loucks r-grid, ri=exp((i-1)*dx-x0)
c     ne, em(ne)   number of energy points, real energy grid
c     edge         chemical potential (energy for k=0)
c     ihole        hole code
c     emu          position of chemical potential in absorption specrum
c     dgc0(nr)     dirac upper component, ground state hole orbital
c     dpc0(nr)     dirac lower component, ground state hole orbital
c     ixc          0  Hedin-Lunqist + const real & imag part
c                  1  Dirac-Hara + const real & imag part
c                  2  ground state + const real & imag part
c                  3  Dirac-Hara + HL imag part + const real & imag part
c                  5  Dirac-Fock exchange with core electrons +
c                     ixc=0 for valence electron density
c     lreal        logical, true for real phase shifts only
c     rmt          r muffin tin
c     xmu          fermi level
c     vi0          const imag part to add to complex potential
c     gamach       core hole lifetime
c     vtot(nr)     total potential, including gsxc, final state
c     edens(nr)    density, hole orbital, final state
c     dmag(251)     density magnetization
c     edenvl      valence charge density
c     dgcn(dpcn)   large (small) dirac components for central atom
c     adgc(adpc)   their development coefficients
c
c     OUTPUT
c     xsec(ne)    atomic absorption cross section to multiply \chi
c                 (atomic background for XMCD)
c     xsnorm(ne)  atomic  absorption cross section (norm for XMCD)
c     rkk(ne, 8)  normalized reduced matrix elements for construction
c                 of termination matrix in genfmt.

      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'
      complex*16 ptz
      dimension ptz(-1:1, -1:1)

      complex*16 em(nex)
      dimension ri(nrptx), vtot(nrptx), edens(nrptx),dmag(nrptx)
      dimension dgc0(nrptx), dpc0(nrptx), vvalgs(nrptx), edenvl(nrptx)
      dimension dgcn(nrptx,30), dpcn(nrptx,30)
c      dimension eorb(30), kappa(30)
      dimension adgc(10,30), adpc(10,30), xnval(30), iorb(-4:3)
      complex*16 rkk(nex, 8), xsec(nex)
      complex*16 bmat(-lx:lx,0:1,8, -lx:lx,0:1,8)
      dimension kind(8), lind(8)
      dimension xsnorm(nex)

      dimension xp(nrptx), xq(nrptx)

c     work space for xcpot
      dimension vxcrmu(nrptx), vxcimu(nrptx), gsrel(nrptx)
      dimension vvxcrm(nrptx), vvxcim(nrptx)

c     work space for fovrg
      complex*16 p(nrptx), q(nrptx), pn(nrptx), qn(nrptx), fscf(nrptx)
      complex*16 pp(nrptx), qp(nrptx), pnp(nrptx), qnp(nrptx)
c     storage for calculation of cross term (SPIN 1 only)
      complex*16 xrcold(nrptx) , xncold(nrptx)

      complex*16  p2, ck, xkmt
c      complex*16  xkmtp, xm1, xm2, xm3, xm4, yvec(nrptx,1)
      complex*16  pu, qu, dum1, factor
      complex*16  xfnorm, xirf, xirf1
      complex*16  temp, aa, bb, cc, rkk1, rkk0, phold
      complex*16  phx(8), ph0
      complex*16  eref

      complex*16 jl,jlp1,nl,nlp1
      complex*16  v(nrptx), vval(nrptx)
      complex*16  xrc(nrptx), xnc(nrptx)
      character*512 slog
      logical ltrace
c     nesvi:  
      complex*16 xrhoce(nex), xrhopr(nex), chia(nex), cchi(nex)
      dimension omega1(nex), bf(0:2, nrptx)

      dimension pat(nrptx),qat(nrptx)
      complex*16 intr(nrptx),var(nrptx) 
c     to pass energy levels and projected DOS
c      dimension neg(30), rhoj(nex,30)
      dimension eng(nex, 30)
c     Josh - Added iPl switch for PLASMON card
c          - and WpCorr = Wi/Wp, Gamma, AmpFac
c          - to describe Im[eps^-1]
      integer iPl, ipole
      double precision WpCorr(MxPole), Gamma(MxPole), AmpFac(MxPole)
c     Josh END

c     explicitly intialize some things
      rkk0  = (0.,0.)
      rkk1  = (0.,0.)
      phold = (0.,0.)
      
      do 5 i=1,MxPole
         WpCorr(1) = -1.d30
 5    continue

      call setkap(ihole, kinit, linit)
c      PRINT*, 'dx=',dx
c     set imt and jri (use general Loucks grid)
c     rmt is between imt and jri (see function ii(r) in file xx.f)
      imt = int((log(rmt) + x0) / dx)  +  1
      jri = imt+1
      jri1 = jri+1
      if (jri1 .gt. nrptx)  call par_stop('jri .gt. nrptx in phase')

c     nesvi: define jnrm
      inrm = int((log(rnrm) + x0) / dx) + 1
      jnrm = inrm + 1

c     We'll need <i|i> later to normalize dipole matrix elements
c     <i|r|f>.  NB, dgc and dpc are r*wave_fn, so use '0' in somm to
c     get integral  psi**2 r**2 dr.
c     Square the dgc0 and dpc0 arrays before integrating.
c     <i|i> == xinorm.
c     dgc and dpc should be normalized <i|i>=1, check this here
      do 10  i = 1, nrptx
         xp(i) = dpc0(i)**2
         xq(i) = dgc0(i)**2
  10  continue
c     nb, xinorm is used for exponent on input to somm
      xinorm = 2*linit + 2
      call somm (ri, xp, xq, dx, xinorm, 0, jnrm)
      del = abs (abs(xinorm) - 1)
      if (del .gt. 1.e-2) then
         write(slog,'(a,i8,1p2e13.5)') ' ihole, xinorm ', ihole , xinorm
         call wlog(slog)
c        if using real phase shifts, don't expect great results
         if (lreal.lt.2)  then
           call wlog(' There may be convergence problems.')
           call wlog(' Xinorm should be 1. If you set the RGRID, '//
     1               'minor interpolation errors ')
           call wlog(' that will not affect final results may occur')
         endif
      endif

c     use ixc for testing
      index = ixc
c       Always use ground state self energy for xsection, quick fix
c       JJR, Jan 93
c       change for testing broadened plasmon pole 6/93
c       index = 2
c   ALA found that it is better to use index=ixc and real part of 
c   self-energy for atomic xsection. 12/96
      ltrace = .true.
      call bcoef(kinit, ipol, ptz, le2, ltrace, ispin, angks, 
     1           kind, lind, bmat)
c     set spin index to use bmat
      isp = 0
      if (ispin.eq.1) isp = nspx - 1

c     zero rkk and phx
      do 20 ie = 1,nex
      do 20 k1 = 1,8
 20   rkk(ie,k1) = 0
      do 30 k1 = 1,8
 30   phx(k1) = 0

      ifirst = 0
c     Josh - if PLASMON card is set, and using HL exc,
c          - read pole information from epsinv.dat
      IF( (iPl.gt.0).and.(ixc.eq.0) ) THEN
         open(file='exc.dat', unit=47, status='old',iostat=ios)
         call chopen (ios, 'exc.dat', 'ffmod2(xsect)')
         DO ipole = 1, MxPole
            call rdcmt(47,'#*cC')
            read(47,*,END=35) WpCorr(ipole), Gamma(ipole), AmpFac(ipole)
            Gamma(ipole)  = Gamma(ipole)/hart
            WpCorr(ipole) = (WpCorr(ipole)/hart) /
     &           SQRT(3.d0/((3 / (4*pi*edens(jri+1))) ** third)**3)
         END DO
 35      CONTINUE
         WpCorr(ipole) = -1.d30
         CLOSE(47)
      END IF
c$$$      IF(ixc.eq.0) THEN
c$$$c        Write wp as calculated from density to sigma.dat
c$$$         open(file='mpse.dat', unit=45, status='replace',iostat=ios)
c$$$         call chopen (ios, 'sigma.dat', 'ffmod2(xsect)')
c$$$         write(45,*) '# ', 'rs      wp(Hartrees)'
c$$$         write(45,*) '# ', (3 / (4*pi*edens(jri+1))) ** third,
c$$$     &        SQRT(3.d0/((3 / (4*pi*edens(jri+1))) ** third)**3)*hart
c$$$         write(45,'(a)')
c$$$     &        '# E-EFermi (eV)   Re[Sigma(E)] (eV)   Im[Sigma(E)] (eV)'
c$$$     &        // '   Re[Z]   Im[Z]   Mag[Z]   Phase[Z]   Lambda(E) (/A)'
c$$$      END IF
c     Josh END
      
      do 400 ie = 1, ne
         iph = 0
c        Josh - xcpot now has new arguments:
c             - iPl, WpCorr, Gamma, AmpFac         
         call xcpot (iph, ie, index, lreal, ifirst, jri,
     1               em(ie), xmu,
     2               vtot, vvalgs, edens, dmag, edenvl,
     3               eref, v, vval, iPl, WpCorr, AmpFac,
     4               vxcrmu, vxcimu, gsrel, vvxcrm, vvxcim,rnrm)
c         call xcpot (iph, ie, index, lreal, ifirst, jri,
c     1               em(ie), xmu,
c     2               vtot, vvalgs, edens, dmag, edenvl,
c     3               eref, v, vval, iPl, WpCorr, Gamma, AmpFac,
c     4               vxcrmu, vxcimu, gsrel, vvxcrm, vvxcim,rnrm)

c       set the method to calculate atomic cross section
c       p2 is (complex momentum)**2 referenced to energy dep xc
        p2 = em(ie) - eref
        p2f = edge - dble(eref)
        ck = sqrt (2*p2 + (p2*alphfs)**2)
        xkmt = rmt * ck

        if (mod(index,10) .lt. 5) then
           ncycle = 0
        else
c          fix later . may be ncycle can be less
           ncycle = 3
        endif
        omega = (dble(em(ie)) - edge) + emu
        omega = max (omega, 0.001d0 / hart)
c       nesvi: add omega1(ie)- need it later
        omega1(ie) = omega

c       remember the bessel functions for multipole matrix elements
        xk0 = omega * alphfs
        ilast = jnrm+6
        if (ilast.lt.jnew) ilast = jnew
        if (ilast.gt.nrptx) ilast = nrptx
        do 50 i = 1, ilast
          temp = xk0 * ri(i)
          if (abs(temp).lt.1.d0) then
c           use series expansion
            do 40 ll = 0,2
              call bjnser(temp,ll, xirf, dum1,1)
              bf(ll,i) = dble(xirf)
 40         continue
          else
c           use formula
            x = dble(temp)
            sinx = sin(x)
            cosx = cos(x)
            bf(0,i) = sinx/x
            bf(1,i) = sinx/x**2 - cosx/x
            bf(2,i) = sinx*(3/x**3-1/x) - 3*cosx/x**2
          endif
 50     continue

c       notice for spin-dep case xsnorm and xsec are spin-dep
c       and kept separately (see call to xsect in subroutine xsph)
        xsnorm(ie) = 0 
        xsec(ie) = 0
        if (dble(em(ie)).lt.-10.d0) goto 400
        if (dimag(p2).le.0.d0 .and. dble(p2).le.0.d0) goto 400

c       matrix elements for multipole (E1,E2,M1) transitions
c       The terms up to (k/c)^2 for absorption are kept.
c       L3 edge: kdif=1 (3d5/2)      kdif=2 (3d3/2), kdif=3(4s)
c       L2 edge: kdif=1 (no transition), 2 (4s),      3 (3d3/2)
        do 350 mult = 0, 2
          if (mult.eq.0) then
c           E1 transitions
            kx = 1
            ks = 2
          else
c           M1 transitions
            kx = 1
            ks = 6
c           E2 transitions
            if (mult.eq.2) kx = 2
          endif 
c         skip unnecessary calculations
          if (mult.gt.0 .and. (mult.ne.le2)) goto 350
 
c         set ilast larger than jri for better interpolation for pu
c         also need 5 points after jri for irregular solution
          ilast = jnrm + 6
          if (ilast.lt.jnew) ilast = jnew

cc        calculate screened dipole field
          ww = dble(emu+p2-edge)
          if (mult.eq.0 .and. izstd.gt.0) then
c            if (ie.eq.1) call correorb(iz, ihole, rmt, jri, dx,ri,
c     1                   p2f,edge, v, dgcn, dpcn, adgc, adpc,
c     2                   eorb, neg, eng, rhoj, kappa, norbp)
            maxsize = 1
            matsize = 0
            sfun = 1.d0
c            call phiscf (ifxc, rmt, ilast, jri, p2, p2f, emu, dx,
c    1                  ri, v, edens, dgcn, dpcn, adgc, adpc,
c     2                  iz, ihole, neg, eng, rhoj,kappa, norbp, fscf,
c     3                  yvec, maxsize, matsize, sfun)
            wse = dble(p2-eng(1,ihole))
          else
            do 159 i = 1, nrptx 
  159       fscf(i) = 1.d0
            wse = ww
          endif
      
c         ww = 1
c         ww = wse / ww
          ww = sqrt(wse/ww)

          do 300 kdif = -kx, kx
            if (omega.le.0.0) goto 300
            ind = kdif + ks
            ikap = kind (ind)
            if (ikap .eq. 0) goto 300
c           use l2lp =0 to include both transitions l-->l+/-1
c           if (l2lp.ne.0) only dipole transitions are calculated.
c            l-->l+1 transitions
            if (l2lp.eq.1 .and. ((kinit.lt.0 .and. ind.ge.3) .or.
     1          (kinit.gt.0 .and. ind.ne.3)) ) goto 300
c            l-->l-1 transitions
            if (l2lp.eq.-1 .and. ((kinit.lt.0 .and. ind.ne.3) .or.
     1          (kinit.gt.0 .and. ind.ge.3)) ) goto 300

            iold = 0
            ic3=0
c           start cycle  do ic3=0,1
c           return for ic3=1 calculations only for |ispin|=1
c           where the central atom cross terms are needed
  100       continue

            irr = -1
c           ic3p=1 to test K2Cr2O7  L3 XES 
            ic3p = ic3
            call dfovrg ( ncycle, ikap, rmt, ilast, jri, p2, dx,
     1      ri, v, vval, dgcn, dpcn, adgc, adpc,
     2               xnval, pu, qu, p, q,
     3               iz, ihole, xion, iunf, irr, ic3p)
            lfin = lind (ind)
            ilp = lfin - 1
            if (ikap .lt. 0) ilp = lfin + 1
            call exjlnl (xkmt, lfin, jl, nl)
            call exjlnl (xkmt, ilp, jlp1, nlp1)
            call phamp(rmt,pu,qu, ck, jl,nl,jlp1,nlp1, ikap, ph0,temp)

            sign = -1.0
            if (ikap.gt.0) sign = 1.0
            factor = ck*alphfs 
            factor = sign * factor/(1+sqrt(1+factor**2))
            dum1 = 1/ sqrt(1+factor**2)
            xfnorm = 1 / temp *dum1
c           normalization factor
c           xfnorm = dum1*rmt*(jl*cos(delta) - nl*sin(delta))/ Rl(rmt)
c           dum1 is relativistic correction to normalization
c           normalize regular solution
            do 130  i = 1,ilast
              p(i)=p(i)*xfnorm
              q(i)=q(i)*xfnorm
  130       continue

cc          calculate xirf including fscf - TDLDA result
            do 140 id = 1, 2
              if (id.eq.1) then
                do 121 j = 1,ilast 
                  pp(j)  = p(j)*dble(fscf(j))
                  qp(j)  = q(j)*dble(fscf(j))
  121           continue
              else
                do 122 j = 1,ilast
                  pp(j)  = p(j)*dimag(fscf(j))
                  qp(j)  = q(j)*dimag(fscf(j))
  122           continue
              endif
              ifl = 1
              if (izstd.gt.0) ifl = -1
              xirf1 = 0
              call radint(ifl, mult, bf, kinit, dgc0,dpc0, ikap, pp,qp,
     1        pn,qn,ri,dx, ilast,iold, xrc,xnc, xrcold,xncold, xirf1)
c             if (ifl.lt.0) xirf1 = xirf1 * xk0 * ww
              if (ifl.lt.0) xirf1 = xirf1 * xk0 
              if (id.eq.1) then
                xirf = xirf1
              else
                if (abs(xirf) .eq. 0.d0) then
                  xirf = xirf1
                elseif (abs(xirf1) .eq. 0.d0) then
                  xirf = xirf
                elseif (abs(xirf1) .lt. abs(xirf)) then
                  dum = abs(xirf1) / abs(xirf)
                  xirf = xirf * sqrt(1.d0 + dum**2)
                else
                  dum = abs(xirf) / abs(xirf1)
                  xirf = xirf1 * sqrt(1.d0 + dum**2)
                endif
              endif
  140       continue

c           note that for real potential  xirf is real or reduced matrix
c           element for dipole transition is pure imaginary.
            if (ic3.eq.0) then
c              can remember only E2 or M1 matrix elements
               if (mult.eq.0 .or. le2.eq.mult) then
                 rkk(ie,ind)=xirf 
                 phx(ind) = ph0
               endif
c              for f' want to include both E2 and M1 for xsnorm and xsec
c              but now only one of them is included (fix later)
               xsnorm(ie)=xsnorm(ie) +
     1         ( dble(xirf)**2 + dimag(xirf)**2 )/ (2*kx+1)
               aa =  - coni*xirf**2
               xsec(ie) = xsec(ie) -  aa * bmat(0,isp,ind, 0,isp,ind)
            elseif (iold.eq.1) then
                rkk1=xirf
                phold = ph0
            elseif (iold.eq.2) then
                rkk0=xirf
            endif

c           get irregular solution and atomic cross-section xsec
c           find irregular solution

            if(dimag(em(ie)).gt.0.d0) then
              irr = 1
c             set pu, qu - initial condition for irregular solution 
              pu = (nl*cos(ph0)+jl*sin(ph0)) *rmt * dum1
              qu=(nlp1*cos(ph0)+jlp1*sin(ph0))*factor *rmt * dum1
             
c             test on bessel functions
c             if (ikap.gt.0) print*,'test1',xkmt**2*(jl*nlp1-nl*jlp1)

              call dfovrg (ncycle, ikap, rmt, ilast, jri, p2, dx,
     1              ri, v,vval, dgcn, dpcn, adgc, adpc,
     1              xnval, pu, qu, pn, qn,
     1              iz, ihole, xion, iunf, irr, ic3p)
cc            set N- irregular solution , which is outside
cc            N=(nlp1*cos(ph0)+jlp1*sin(ph0))*factor *rmt * dum1
cc            N = i*R - H*exp(i*ph0)
              temp = exp(coni*ph0)
              do i = 1, ilast
                pn(i) = coni * p(i) - temp * pn(i)
                qn(i) = coni * q(i) - temp * qn(i)
              enddo
            else
              do 150 i = 1, ilast
                pn(i) = 0
                qn(i) = 0
  150         continue
            endif

c           combine regular and irregular solution into the
c           central atom absorption coefficient xsec (mu = dimag(xsec))
c           thus for real energy dimag(xsec)=xsnorm

c           also include TDLDA effects
            do 170 id = 1, 2
              if (id.eq.1) then
                do 131 j = 1,ilast
                  pp(j)  = p(j)*dble(fscf(j))
                  qp(j)  = q(j)*dble(fscf(j))
                  pnp(j)  = pn(j)*dble(fscf(j))
                  qnp(j)  = qn(j)*dble(fscf(j))
  131           continue
              else
                do 132 j = 1,ilast
                  pp(j)  = p(j)*dimag(fscf(j))
                  qp(j)  = q(j)*dimag(fscf(j))
                  pnp(j)  = pn(j)*dimag(fscf(j))
                  qnp(j)  = qn(j)*dimag(fscf(j))
  132           continue
              endif

c           TDLDA theory is written for the r-form of matrix elements
c           so one might want to use ifl=-1,-2 for these calculations
c           on the other hand want ifl=1,2 for DANES calculations
c           since it is more reliable at high energies and gives
c           better results for Cu test.
              ifl = 2
              if (izstd.gt.0) ifl = -2

              call radint(ifl,mult, bf, kinit, dgc0, dpc0, ikap, pp, qp,
     1            pnp, qnp, ri,dx, ilast,iold, xrc, xnc, xrcold, xncold,
     2            xirf1)
              if (ifl.lt.0) xirf1 = xirf1 * xk0**2 * ww**2
              if (id.eq.1) then
                xirf = xirf1
              else
                if (abs(xirf) .eq. 0.d0) then
                  xirf = xirf1
                elseif (abs(xirf1) .eq. 0.d0) then
                  xirf = xirf
                elseif (abs(xirf1) .lt. abs(xirf)) then
                  dum = abs(xirf1) / abs(xirf)
                  xirf = xirf * sqrt(1.d0 + dum**2)
                else
                  dum = abs(xirf) / abs(xirf1)
                  xirf = xirf1 * sqrt(1.d0 + dum**2)
                endif
              endif
  170       continue

            if (ic3.eq.0) then
               xsec(ie) = xsec(ie) - xirf * bmat(0,isp,ind, 0,isp,ind)
            endif

c           ------start of density of states part------------- 
c           added by nesvi 07/12/00
c
c           Calculate rhoc00 and rho_projected on 
c           the same grid as xsect. Need this to calculate the smooth
c           atomic ratio rho_0/mu_0 or rho_proj/mu_0.              
c           The atomic functions are normalized to 1 inside Norman radius.
c           This procedure can be called 'Renormalized atomic sphere method".
c           It gives very reasonable numbers for hole counts. In order to
c           get Mulliken counts one can extend integration limits till very
c           large R, but it's currently not recommended because of the problems
c           with the wave function's tails above Rnm.
 

            jproj =  iorb(ikap)
            if (jproj.eq.0 .and. ikap.lt.0) jproj = iorb(-ikap-1)
            kdif1 = -1
            if(kinit.gt.0) kdif1 =  1
                
            if (kdif .eq. kdif1 .and. ic3.eq.0 .and. jproj.gt.0) then
c              calculate rhoc00 (rho_0)

               temp = (2*lfin+1.0d0) / (1+factor**2) /pi *4*ck /hart
               do 500 i = 1, ilast
                 xrc(i) = pn(i)*p(i) - coni*p(i)*p(i) 
     1                   + qn(i)*q(i) - coni*q(i)*q(i)
  500          continue    
               xirf = 1
c              integration is till Norman radius, not Rmt as in xsect
               i0 = jnrm + 1
               call csomm2 (ri, xrc, dx, xirf, rnrm, i0)
               xrhoce(ie) = - xirf * temp
            
c              calculate rho_projected:              

c              pat, qat - atomic functions that we make projection on.
               do 510 i=1,nrptx
                 pat(i) = dgcn(i,jproj)
                 qat(i) = dpcn(i,jproj)
  510          continue

c     normalize pat and qat in the Norman radius sphere: <n|n>=1,
c     (renormalized atomic sphere method)
     
               do 520  i = 1, ilast
                  xp(i) = pat(i)**2 + qat(i)**2
                  xq(i) = 0
  520          continue
c     nb, xinorm is used for exponent on input to somm 
               xinorm = 2*lfin + 2
               call somm2 (ri, xp, dx, xinorm, rnrm, 0, i0)
c              call somm (ri, xp, xq, dx, xinorm, 0, jnrm)
      
               xinorm = sqrt(xinorm)
               do 530 i=1,nrptx
                  pat(i) = pat(i) / xinorm
                  qat(i) = qat(i) / xinorm
  530          continue
  
c              calculate overlap integral between f and atomic function
c              (integral Rl(r)*Psi_at(r)dr from 0 till r') 
c              intr(i) is that overlap integral. Later it
c              will be multiplied by pr(i)*Psi_at(r') and integrated 
c              till Norman radius.

               do 540 i=1,ilast
                  var(i)=pat(i)*p(i)+qat(i)*q(i)
c                 factor of 2 -integration r< r>  -->2 r r'
  540          continue

c              integration by trapezoid method
               intr(1)=var(1)*ri(1)
               do 550 i=2,ilast
                  intr(i)=intr(i-1)+ (var(i)+var(i-1))*(ri(i)-ri(i-1))
  550          continue 


c         now calculate rho_projected - xrhopr
               temp = (2*lfin+1.0d0) / (1+factor**2) /pi *4*ck /hart
c              temp = abs(ikap) / (1+factor**2) /pi *4*ck /hart
               do 560  i = 1, ilast
                 xrc(i) = pn(i)*pat(i)*intr(i)+ 
     1                    qn(i)*qat(i)*intr(i)
                 xrc(i) = xrc(i) - coni*(p(i)*pat(i)*intr(i) + 
     1                    q(i)*qat(i)*intr(i))
  560          continue

               xirf =  1
               call csomm2 (ri, xrc, dx, xirf, rnrm, i0)
               xrhopr(ie) = - xirf * temp
    
            endif
c           ----------end of density of states part---


            if (iold.gt.0) then
c             calculate cross term contribution to XMCD
c             in both cases coupling between neighbors 
c             need to remove SO interaction (ic3=1) in order
c             to avoid unphysical peak in Gd XMCD. a.l. ankudinov
              k1 = ind - 1
              if (k1.ge.1 .and.k1.le.8) then
              if (lind(k1).eq.lind(ind) .and. lind(k1).gt.0) then
                aa = exp( coni*(ph0 - phold))
                bb = 1/aa
                cc = - ( bmat(0,isp,k1, 0,isp,ind) +
     1                 bmat(0,isp,ind, 0,isp,k1) ) / 2.d0
                xsec(ie) = xsec(ie) - coni * rkk1 * rkk0 * (bb+aa) * cc
cc              combine regular and irregular solution into the
cc              central atom absorption coefficient (mu=dimag(xsec))
cc              thus for real energy dimag(xsec)=xsnorm
                call radint (3, mult, bf, kinit, dgc0, dpc0, ikap, p, q,
     1            pn, qn, ri, dx, ilast, iold, xrc, xnc, xrcold, xncold,
     2            xirf)
                xsec(ie) = xsec(ie) + xirf * cc * bb
  
                call radint (4, mult, bf, kinit, dgc0, dpc0, ikap, p, q,
     1            pn, qn, ri, dx, ilast, iold, xrc, xnc, xrcold, xncold,
     2             xirf)
                xsec(ie) = xsec(ie) + xirf * cc * aa
              endif
              endif
            endif
cc          end of |ispin=1| cross term calculations

c           prepare for ic3=1 cross term calculations if needed
            if (ic3.eq.0 .and. abs(ispin).eq.1) then
              iold = 0
              if (ind.lt.8 .and. lind(ind).gt.0) then
                k1 = ind + 1
                if (lind(k1).eq.lind(ind)) iold = 1
              endif
              if (ind.gt.1 .and. lind(ind).gt.0) then
                k1 = ind - 1
                if (lind(k1).eq.lind(ind)) iold = 2
              endif
c             need to remove SO interaction to calculate cross term
c             big effect for Gd XMCD calculations
              if (iold.gt.0) then
c               repeat calculation for current kdif with SO turned off
                ic3 = 1
                goto 100
              endif
            endif

  300     continue
  350   continue

        if (omega.gt.0.0) then
c         prefac = (8 * pi / 3)  * alphfs * omega  -- nonrelativistic
c         relativistic is (for alpha form)
          prefac = 4 * pi * alpinv / omega * bohr**2
          xsnorm(ie) =  xsnorm(ie) * prefac * 2*abs(ck) 
          xnorm= sqrt( xsnorm(ie) )
          xsec(ie) = xsec(ie) * prefac* 2*ck

c         put complex sqrt(prefactor) into reduced matrix elements rkk
          ck = sqrt ( prefac * (2*ck))
c         guarantee that we have the right root
          if (dimag(ck) .lt. 0) ck = -ck
c         add central atom phase shift here. 
          do 360 kdif = 1 , 8
 360      rkk(ie,kdif)= rkk(ie,kdif) * ck/xnorm * exp(coni*phx(kdif))
        endif
 400  continue
c     end of energy cycle

c     Josh - Close sigma.dat
      close(45)
c     Josh END

      if (ipr2.ge.3) then
c       calculate mu_0/rho_0 for XMCD normalization.
        do 410 ie=1,ne
           chia(ie) = 0
  410   continue
        vrcorr = 0
        vicorr = 0
        call xscorr(1, em, ne1, ne, ik0, xrhoce,xsnorm,chia,
     1     vrcorr, vicorr, cchi)
        do 420 ie = 1, ne1
            xrhoce(ie)  = coni* dimag(xrhoce(ie)+cchi(ie))
  420   continue
        call xscorr(1, em, ne1, ne, ik0, xrhopr,xsnorm,chia,
     1     vrcorr, vicorr, cchi)
        do 425 ie = 1, ne1
            xrhopr(ie)  = coni* dimag(xrhopr(ie)+cchi(ie))
  425   continue    
        call xscorr(1, em, ne1, ne, ik0, xsec,xsnorm,chia,
     1     vrcorr, vicorr, cchi)
        do 430 ie = 1, ne1
            cchi(ie)  = coni* dimag(xsec(ie)+cchi(ie))
  430   continue

        open(unit=3,file='ratio.dat',status='unknown', iostat=ios)
        open(unit=4,file='ratiop.dat',status='unknown', iostat=ios)
c       normalize to xsec at 50 ev above edge
        edg50 = emu +50.0 / hart
        call terp (omega1, xsnorm, ne1, 1, edg50, xsedge)
        write(3,440) xsedge, emu * hart 
  440   format ('# Normalization factor:', e12.4,
     1     ' Angstrom**2. Fermi level at ', f7.1, ' eV.')
        write(3,450)
  450   format ('#   Energy      rho_0        mu_0       rho_0/mu_0 ')
     
        write(4,440) xsedge, emu * hart 
        write(4,455)
  455   format ('#   Energy      rho_proj      mu_0      rho_proj/mu_0',
     1   '    mu_deloc ')

        do 470 ie=1,ne1 
           if (dimag(cchi(ie)).eq.0.d0 .and. ie.lt.ik0) then
              cchi(ie)=cchi(ik0)
              xrhoce(ie)=xrhoce(ik0)
              xrhopr(ie)=xrhopr(ik0)
           endif
           ratio = dimag(xrhoce(ie)) / dimag(cchi(ie)) * xsedge
           ratiop = dimag(xrhopr(ie)) / dimag(cchi(ie)) * xsedge

           write(3,460)  dble(em(ie))*hart, dimag(xrhoce(ie)),
     1          dimag(cchi(ie))/xsedge, ratio*corr
c          corr is the ratio N_av/N_j, responsible for difference in
c          counts due to variation of wave function due to spin-orbit
  460      format(f12.6, 2x, e12.6,2x,e12.6,2x,e12.6,1x,e12.6)      
           write(4,465)  dble(em(ie))*hart, dimag(xrhopr(ie)),
     1          dimag(cchi(ie))/xsedge, ratiop,
     2          dimag(xrhoce(ie)-xrhopr(ie))/ratio 
c     also write contribution to mu_0 from delocalized states defined as
c     (rho-rho_proj)/ratio 
  465      format(f12.6, 2x, e12.6,2x,e12.6,2x,e12.6,1x,e12.6,2x,e12.6)    
      
  470   continue    
        close(unit=3)
        close(unit=4)
      endif 

      return
      end
