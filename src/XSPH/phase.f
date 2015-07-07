c     Josh - argument iPl has been added to arguments of xsect
      subroutine phase (iph, dx, x0, ri, ne, ne1, ne3, em,
     1                  ixc, nsp, lreal, rmt,rnrm, xmu,
     2                  iPl,
     2                  vtot, vvalgs, edens, dmag, edenvl,
     3                  dgcn, dpcn, adgc, adpc, eref, ph, lmax,
     2                  iz, ihole, xion, iunf, xnval, ispin)
c     2                  vi0, iPl, gamach,

      implicit double precision (a-h, o-z)

c     INPUT
c     iph          unique pot index (used for messages only)
c     dx, x0, ri(nr)
c                  Loucks r-grid, ri=exp((i-1)*dx-x0)
c     ne, em(ne)   number of energy points, real energy grid
c     ixc        0  Hedin-Lunqist + const real & imag part
c                  1  Dirac-Hara + const real & imag part
c                  2  ground state + const real & imag part
c                  3  Dirac-Hara + HL imag part + const real & imag part
c                  4, 5, 6, see rdinp or xcpot
c     lreal        1 for real self energy and 2 for real phase shifts 
c     rmt          r muffin tin
c     xmu          fermi level
c     vi0          const imag part to add to complex potential
c     gamach       core hole lifetime
c     vtot(nr)     total potential, including gsxc
c     vvalgs(nr)   overlap Coulomb+gsxc potential for valence electrons
c     edens(nr)    density
c     dmag(nr)     density magnetization
c     edenvl(nr)  valence charge density
c     dgcn(dpcn)   large (small) dirac components for 'iph' atom
c     adgc(adpc)   their development coefficients
c
c     OUTPUT
c     eref(ne)     complex energy reference including energy dep xc
c     ph(nex,ltot+1) complex scattering phase shifts
c     lmax         max l (lmax = kmax*rmt)

      include '../HEADERS/dim.h'
      include '../HEADERS/const.h'

      complex*16 em(nex)
      dimension  ri(nrptx), vtot(nrptx), edens(nrptx)
      dimension  dmag(nrptx), vvalgs(nrptx), edenvl(nrptx)
      dimension  adgc(10,30), adpc(10,30), xnval(30)
      dimension  dgcn(nrptx,30), dpcn(nrptx,30)
      complex*16  eref(nex)
      complex*16  ph(nex,-ltot:ltot)
      integer ispin

c     work space for xcpot
      dimension   vxcrmu(nrptx), vxcimu(nrptx), gsrel(nrptx)
      dimension   vvxcrm(nrptx), vvxcim(nrptx)
c     p and q were needed in xsect to calc. matrix elements.
      complex*16 p(nrptx), q(nrptx)

      complex*16  p2, ck, xkmt, temp, pu, qu
      complex*16 jl(ltot+2), nl(ltot+2)
c      complex*16 xkmtp, nlp(ltot+2), jlp(ltot+2)

      complex*16 v(nrptx), vval(nrptx)
      character*512 slog
c     Josh - Added iPl switch for PLASMON card
c          - and WpCorr = Wi/Wp, Gamma, AmpFac
c          - to describe Im[eps^-1]
      integer iPl, ipole
      double precision WpCorr(MxPole), Gamma(MxPole), AmpFac(MxPole),
     &     rnrm
c     Josh END

c{#mn: g77 (and other compilers) have an intrinsic function besjn, 
c      so besjn should be declared  external 
         external besjn
c#mn}

      do 5 i=1,MxPole
         WpCorr(1) = -1.d30
 5    continue


      clz = 0.d0

c     zero phase shifts (some may not be set below)
      xkmax = 0
      ne12 = ne - ne3
      do 100  ie = 1, ne
         do 90  il = -ltot, ltot
            ph(ie,il) = 0
   90    continue
         if (ie.le.ne12 .and. xkmax.lt.dble(em(ie))) xkmax= dble(em(ie))
  100 continue
      xkmax = sqrt(xkmax * 2)

c     Use kmax to find accurate l-points
c     limit l, lmax = prefac* kmax * rmt
c     prefac is set not to have warning message for Cu metal for kmax=20
      prefac = 0.7d0
      lmax = int(prefac * rmt * xkmax)
      lmax = max(lmax, 5)
      if (lmax.gt.ltot) then
        ik = nint( ltot / rmt / bohr / prefac )
        write (slog, 110) ik
  110   format('      Phase shift calculation is accurate to k=', i2)
        call wlog(slog)
        write (slog, 120)
  120   format('      See FEFF document to increase the range.')
        call wlog(slog)
      endif
      lmax = min (lmax, ltot)
c     set imt and jri (use general Loucks grid)
c     rmt is between imt and jri (see function ii(r) in file xx.f)
      imt = int((log(rmt) + x0) / dx)  +  1
      jri = imt+1
      jri1 = jri+1
      if (jri1 .gt. nrptx)  call par_stop('jri .gt. nrptx in phase')

      ifirst = 0
      index = ixc
c     Josh - if PLASMON card is set, and using HL exc,
c          - read pole information from epsinv.dat
      IF( (iPl.gt.0).and.(ixc.eq.0) ) THEN
         open(file='exc.dat', unit=47, status='old',iostat=ios)
         call chopen(ios,'exc.dat','ffmod2(phase)')
         DO ipole = 1, MxPole
            call rdcmt(47,'#*cC')
            read(47,*,END=125) WpCorr(ipole), Gamma(ipole),
     &           AmpFac(ipole)
            Gamma(ipole)  = Gamma(ipole)/hart
            WpCorr(ipole) = (WpCorr(ipole)/hart) /
     &           SQRT(3.d0/((3 / (4*pi*edens(jri+1))) ** third)**3)
         END DO
 125     CONTINUE
         WpCorr(ipole) = -1.d30
         CLOSE(47)
      END IF
c$$$      IF(ixc.eq.0) THEN
c$$$c     Write wp as calculated from density to sigma.dat
c$$$         open(file='mpse.dat', unit=45, status='replace', iostat=ios)
c$$$         call chopen(ios, 'sigma.dat', 'ffmod2(phase)')
c$$$         write(45,*) '# ', 'rs      wp(eV)'
c$$$         write(45,*) '# ', (3 / (4*pi*edens(jri+1))) ** third, 
c$$$     &        SQRT(3.d0/((3 / (4*pi*edens(jri+1))) ** third)**3)*hart
c$$$         write(45,*) '# mu (eV)'
c$$$         write(45,*) '# ', xmu
c$$$         write(45,'(a)') 
c$$$     &         '# E-EFermi (eV)   Re[Sigma(E)] (eV)   Im[Sigma(E)] (eV)'
c$$$     &       // '   Re[Z]   Im[Z]   Mag[Z]   Phase[Z]   Lambda(E) (/A)'
c$$$      END IF
c     Josh END
      
c     calculate phase shifts
      do 220 ie = 1, ne12

c        Josh - xcpot now has new arguments:
c             - iPl, WpCorr, Gamma, AmpFac         
         call xcpot (iph, ie, index, lreal, ifirst, jri,
     1               em(ie), xmu,
     2               vtot, vvalgs, edens, dmag, edenvl,
     3               eref(ie), v, vval, iPl, WpCorr, AmpFac,
     4               vxcrmu, vxcimu, gsrel, vvxcrm, vvxcim, rnrm)
c         call xcpot (iph, ie, index, lreal, ifirst, jri,
c     1               em(ie), xmu,
c     2               vtot, vvalgs, edens, dmag, edenvl,
c     3               eref(ie), v, vval, iPl, WpCorr, Gamma, AmpFac,
c     4               vxcrmu, vxcimu, gsrel, vvxcrm, vvxcim, rnrm)

         if (dble(em(ie)).lt.-10.d0 .or. dble(em(ie)) .gt.3.d2) goto 220
c        p2 is (complex momentum)**2 referenced to energy dep xc
c        notice that constant Im part (gamach/2+vi0) is cancelled,
c        since it is also present in v and vval.
         p2 = em(ie) - eref(ie) 
         if (lreal.gt.1 .and. ie.le.ne1) p2 = dble(p2)
         ck =  sqrt (2*p2+ (p2*alphfs)**2)
         xkmt = rmt * ck
         if (dble(p2).le.0.d0 .and. dimag(p2) .le.0.d0) goto 220

c{#mn  see note above about declaring besjn as external
c#mn}
         call besjn (xkmt, jl, nl)

         if (mod(ixc,10) .lt. 5) then
             ncycle = 0
         else
             ncycle = 3
         endif

         do 210  ll = -lmax, lmax
            il = abs(ll) +  1
c           nonlocal exchange is unstable for high il.
c           need to do integrals instead of diff. eq. fix later
c           use local xc for high il
            if (il*dx.gt.0.50) then
               ncycle=0
            endif

c  v should be V_N+V_COUL+V_XCtotal-V_mt, vval= V_N+V_COUL+V_XCVAL-V_mt
            ikap = ll - 1
            if ( ll.gt.0 ) ikap=ll
            ilp = il + 1
            if (ikap.gt.0) ilp = il - 1
            ic3 = 0

            if(nsp.eq.1 .and. ispin.eq.0) then
c              remove spin-orbit interaction
c              otherwise, get wrong results e.g. for Pt metal
               if (ll.ne.0) ic3 = 1
               ikap = -1 - abs(ll)
               ilp = il + 1
            endif

c_lz  add term (C L_z) (p.32 of Ankoudinov's thesis) 
c     currently just add constant potential only within mt radius
c     keep intersitial level the same
c OPC for U for jj coupling
c           if (ll.eq.3 .and. iph.eq.1) then
c              clz = -0.5d0 / hart
c              if (ispin.lt.0) clz = -clz
c              do 180 i = 1, jri
c                 v(i) = v(i) + clz
c                 vval(i) = vval(i) + clz
c180           continue
c           endif
c OPC for U for LS coupling
            if (abs(ll).eq.3 .and. iph.eq.1 .and. ispin.eq.1) then
               clz = -0.0d0 / hart
               if (ikap.lt.0) clz = -clz
               do 180 i = 1, jri
                  v(i) = v(i) + clz
                  vval(i) = vval(i) + clz
 180           continue
            endif

c           never use irr=0, only positive or negative
            irr = -1
            call dfovrg (ncycle, ikap, rmt, jri, jri, p2, dx,
     1               ri, v,vval, dgcn, dpcn, adgc, adpc,
     1               xnval, pu, qu, p, q,
     1               iz, ihole, xion, iunf, irr, ic3)

c        restore potential for clz=0
c OPC for U for jj coupling
c           if (ll.eq.3 .and. iph.eq.1) then
c OPC for U for LS coupling
            if (abs(ll).eq.3 .and. iph.eq.1 .and. ispin.eq.1) then
               do 190 i = 1, jri
                  v(i) = v(i) - clz
                  vval(i) = vval(i) - clz
 190           continue
            endif
            call phamp (rmt, pu, qu, ck, jl(il), nl(il),
     1                  jl(ilp), nl(ilp), ikap, ph(ie,ll), temp)

c           cut phaseshift calculation if they become too small
            if (abs(ph(ie,ll)) .lt. 1.0e-6 .and. ll.ge.4)  goto 220
c           new cut function introduced by Rivas
            if(abs(exp((0,2)*ph(ie,ll))-1.).lt.1.0e-5) ph(ie,ll)=0
            if (abs(ph(ie,ll)) .lt. 1.0e-5 .and. ll.ge.4)  goto 220

  210    continue
  220 continue
c     Josh - Close sigma.dat
      close(45)
c     Josh END

      do 230 ie = ne12+1, ne
  230 eref(ie) = eref(ne1)

      return
      end
