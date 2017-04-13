      subroutine xcpot (iph, ie, index, lreal, ifirst, jri,
     1                  em, xmu,
     2                 vtot, vvalgs, densty, dmag, denval,
     3                  eref, v, vval, iPl, WpCorr, AmpFac,
     4                  vxcrmu, vxcimu, gsrel, vvxcrm, vvxcim,rnrm)
c      subroutine xcpot (iph, ie, index, lreal, ifirst, jri,
c     1                  em, xmu,
c     2                 vtot, vvalgs, densty, dmag, denval,
c     3                  eref, v, vval, iPl, WpCorr, Gamma, AmpFac,
c     4                  vxcrmu, vxcimu, gsrel, vvxcrm, vvxcim,rnrm)

      implicit double precision (a-h, o-z)
c     calculate self-energy correction
c     first coded j. mustre de leon
c     last modified a.ankudinov 1996 for non-local self-energies
c     Ankudinov, Rehr, J. Physique IV, vol. 7, C2-121, (1997).

c     INPUT
c     iph, ie used only for debug and labels.
c     index       0  Hedin-Lunqvist + const real & imag part
c                 1  Dirac-Hara + const real & imag part
c                 2  ground state + const real & imag part
c                 3  Dirac-Hara + HL imag part + const real & imag part
c                 4  See rdinp for comment
c     lreal       not equal zero for real self energy
c     ifirst      first entry flag, set to zero before first call for
c                 each unique potential, see vxcrmu and vxcimu below
c     jri         index of first interstitial point in current
c                 Loucks r grid
c     em          current energy grid point
c     xmu         fermi level
c     vi0         const imag part to subtract from potential
c     gamach      core hole lifetime
c     vtot(nr)    total potential (coulomb and gs exchange corr)
c     vvalgs(nr)  total coulomb + gs xc potential from valence electrons
c     densty(nr)  electron density
c     dmag(nr)    density magnetization
c     denval(nr)  valence electron density
c     iPl         Control for many pole self energy (Josh)
c
c     OUTPUT
c     eref        complex energy reference for current energy
c     v(nr)       complex potential including energy dep xc
c     vval(nr)    as above,but xc from valence electrons only
c     em          current energy
c
c     WORKSPACE
c     vxcrmu and vxcimu are calculated only on first entry for a
c     particular unique potential, re-used on subsequent entries.
c     vxcrmu(nr)  real part of xc at fermi level
c     vxcimu(nr)  imag part of xc at fermi level
c     gsrel(nr) ratio of gs xc potentials with and without magnetization
c     vvxcrm(nr)  real part of xc at fermi level from valence electrons
c     vvxcim(nr)  imag part of xc at fermi level from valence electrons


      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'

      dimension   vtot(nrptx), vvalgs(nrptx), densty(nrptx)
      dimension   dmag(nrptx), denval(nrptx)
      complex*16  em, eref, v(nrptx), vval(nrptx)
      dimension   vxcrmu(nrptx), vxcimu(nrptx)
      dimension   vvxcrm(nrptx), vvxcim(nrptx), gsrel(nrptx)

c     Josh added variables:
c     ZRnrm      - renormalization constant
c     ZTmp       - Temp var for ZRnrm
c     csig       - control for using many pole self energy
c     NRPts      - Number of points to use inside atom.
c                  Other points are linearly interpolated.
c     WpCorr     - Array of frequencies for many pole self energy.
c     rsTmp      - Temp var for rs
c     WpTmp      - Temp var for Wp
c     AmpFac     - g_i (pole strengths)
c     Gamma      - pole broadening
c     RsInt      - Rs in the intersitial
c     DRs        - RsCore - RsInt
c     delrHL     - Re[delta sigma] for many pole self energy
c     deliHL     - Im[delta sigma] for many pole self energy
c     Rs1(NRPts) - Array of Rs points for interpolation
      complex*16  delta, deltav, ZRnrm, ZTemp
      character*512 slog
      logical csig, rdmpse
      integer NRPts, lastPl
      parameter (tol=0.0004d0)
      parameter (NRPts=10)
      double precision WpCorr(MxPole), AmpFac(MxPole)
c      double precision  Gamma(MxPole), rsTmp, WpTmp
      double precision RsInt, DRs, delrHL(NRPts), deliHL(NRPts),
     &     Rs1(NRPts), dx, x0, volume, totvol, rnrm, ri, riInt, omp,
     &     ompmax
      complex*16 delavg
c    Josh END

c     First calculate vxc to correct the local momentum dispersion
c     relation, delta = vxc(e,k) - vxc(mu,k), and
c               p^2 = k^2 -mu + kf^2 - delta.
c     In jr theory, v(e,r) = vcoul(r) + vxc(e,r) =
c                          = vcoul(r) + vxcgs(r) + delta(e,r).

c    at jri potential is smooth continuation of potential to r(jri)
c    at this point potential jumps to interstitial value at jri+1
c     Atom r grid
      dx = 0.05d0
      x0 = 8.8d0
      totvol = 4.d0/3.d0*pi*rnrm**3
      delavg = 0.d0
c      totvol = 0.d0
      ZRnrm = 0.d0
      csig=.false.
      ompm1 = 0.d0
      jri1 = jri + 1
c      rsTmp = RsCorr
      xfval=0
      delvr = 0
      delvi = 0
      lastPl = 0
      nmax=1
      nul=0
      ibp = index / 10
      ixc = mod(index,10)
      ixcTmp=ixc
      DO i = 1, MxPole
         IF(WpCorr(i).le.0.d0) then
            lastPl = i-1
            GOTO 5
         END IF
      END DO
   5  CONTINUE
      if((ixc.eq.0).and.(iPl.gt.0)) then
         csig=.true.
      end if
      if (ixc .eq. 2 .or. dble(em).le.xmu)  then
         do 10  i = 1, jri1
            v(i) = vtot(i)
            vval(i) = vvalgs(i)
   10    continue
c        Ground state exchange, no self energy calculation
         goto 888
      endif

c     Josh - Added CSigma to calculate HL Sigma with broadening and such.
c     Calculate Rs at the core and interstitial densities.
      if (densty(jri1).le.0) then
         RsInt =10
      else
         RsInt = (3 / (4*pi*densty(jri1))) ** third
      endif
      if (densty(1).le.0) then
         rscore =101.d0
      else
         rscore = (3 / (4*pi*densty(1))) ** third
      endif
      DRs = (RsInt-rscore)/(NRPts-1)
      omp = SQRT(3.d0/RsInt**3)*hart
      ompmax=omp*WpCorr(lastPl)
c     Now calculate delta sigma as a function of Rs and energy
      if (csig) then
         do i= NRPts, 1, -1
            rdmpse = .false.
            delrHL(i) = 0.d0
            deliHL(i) = 0.d0
            Rs1(i)=rscore+DBLE(i-1)*DRs

c           If iPl > 1, use renormalization, else not
c           If iPl = 2, use Sigma(r) = Sigma[Wp(r)*Wp/Wp(RsInt)]
c           If iPl = 3, use Sigma(r) = 0 outside of intersitial region
c                       actually linearly interpolates to zero at the
c                       first RPt.
c           If iPl > 3, use Sigma(r) = Sigma(RsInt) (Sigma as a bulk property)
            if(iPl.gt.1) then
               if((iPl.eq.2).or.(i.eq.NRPts)) then
c                  call CSigZ(em,xmu,Rs1(i),delrHL(i),deliHL(i),ZTemp,
c     &                 WpCorr,Gamma,AmpFac)
                  call CSigZ(em,xmu,Rs1(i),delrHL(i),deliHL(i),ZTemp,
     &                 WpCorr,AmpFac)
               elseif(iPl.eq.3) then
                  delrHL(i) = 0.d0
                  deliHL(i) = 0.d0
               else
                  delrHL(i) = delrHL(NRPts)
               end if
               if(i.eq.NRPts) ZRnrm = ZTemp
            else
c               call CSigma(em,xmu,Rs1(i),delrHL(i),deliHL(i),WpCorr,
c     &              Gamma,AmpFac)
               call CSigma(em,xmu,Rs1(i),delrHL(i),deliHL(i),WpCorr,
     &              AmpFac)
            end if
c            Josh Kas - Write self energy to mpse.bin for fast processing later
c            write(23,'(6f30.10)') DBLE(em), Rs1(i), delrHL(i),
c     &           deliHL(i), Ztemp
c 17         continue
c           debugging output of deltaSigma(em, rs)
c           write(44,'(6f30.10)') DBLE(em), Rs1(i), delrHL(i), deliHL(i),
c     &           dble(ZRnrm), dimag(ZRnrm)
         end do
c        write(44,*)
      end if
c     END Josh

c     Add the self energy correction
      do 20  i =  jri1,1,-1
         ri = exp((i-1)*dx - x0)
         if(i.eq.jri1) then
            riInt = ri
         end if
         niter = 0
         if (densty(i).le.0) then
            rs =10
         else
            rs = (3 / (4*pi*densty(i))) ** third
         endif
c         write(22,*) 1.d0*exp(dble(i)*0.01), densty(i)
c        Josh - If csigma is turned on, interpolate onto rs.
c        Then skip to 15 (skip other calculations and self
c        consistency)
         if(csig) then
            omp = SQRT(3.d0/rs**3)*hart
            if(iPl.ge.4) then
               delr = delrHL(NRPts)
               deli = deliHL(NRPts)
            else
               call terp (Rs1, delrHL, NRPts, 1, rs, delr)
               call terp (Rs1, deliHL, NRPts, 1, rs, deli)
            end if
            if((iPl.ne.5).or.(omp.lt.ompmax)) then
               goto 15
            end if
         end if
c        END Josh

c        xf = 1.9191.../rs
         xf = fa / rs
         rsm = rs / (1+dmag(i))**third
         xfm = fa / rsm

         if (ixc.eq.5) then
            if ( denval(i) .gt. 0.00001) then
               rsval = (3 / (4*pi*denval(i))) ** third
               if (rsval.gt.10.0) rsval=10.0
            else
               rsval = 10.0
            endif
            xfval = fa / rsval
         elseif (ixc.ge.6) then
            if (densty(i) .le. denval(i) ) then
               rscore = 101.0
            else
               rscore = (3 / (4*pi*(densty(i)-denval(i)))) ** third
            endif
         endif

         if (ifirst .eq. 0)  then
c           vxc_mu indep of energy, calc only once
c           Calculate vxc at fermi level e = mu, j.m. 1/12/89
            xk = xf * 1.00001
            gsrel(i) = 1.0d0
            if (ixc .lt. 5) then
              call sigma(ixc, ibp,rs,rscore,xk,vxcrmu(i),vxcimu(i))
              if (index .eq. 0) then
c  do not need 4 following lines for gs difference in potential
c                xmag = 1.0d0+ dmag(i)
c                call vbh(rs,xmag,v1)
c                call vbh(rs, 1.0d0,v0)
c                if (v0 .ne. 0) gsrel(i) = v1/v0
              endif
            else
              call sigma(nul,ibp, rs, rscore,xk,vxcrmu(i),vxcimu(i))
            endif
            if (ixc.eq.5 ) then
               xkpp = xfval * 1.00001
               call sigma
     1         (ixc, ibp, rsval, rscore, xkpp, vvxcrm(i),vvxcim(i))
               if (ixc.eq.5 .and. i.eq.jri1) then
                  vvxcrm(jri1) =  vxcrmu(jri1)
                  vvxcim(jri1) =  vxcimu(jri1)
               endif
            elseif (ixc .ge. 6) then
               call sigma
     1         (ixc, ibp, rs, rscore, xk, vvxcrm(i), vvxcim(i))
               if (ixc.eq.6 .and. i.eq.jri1) then
                  vvxcrm(jri1) =  vxcrmu(jri1)
                  vvxcim(jri1) =  vxcimu(jri1)
               endif
            else
               vvxcrm(i) = 0.0d0
               vvxcim(i) = 0.0d0
            endif
         endif

c        xk2 is the local momentum squared, p^2 = k^2 - 2*mu + kf^2,
c        k^2 represents energy measured from vacuum.
c        See formula 2.15 in Lee and Beni's paper with the last 2
c        terms neglected.  (complete reference?)
         xk2 = 2 * (dble(em) - xmu) + xf**2
         xk = sqrt(xk2)
         xkm2 = 2 * (dble(em) - xmu) + xfm**2
c        quick fix
         if (xkm2.lt.0) xkm2=xk2
         xkm = sqrt(xkm2)

c        find \delta_1
         if (ixc .lt. 5) then
            call sigma (ixc, ibp, rs, rscore, xk, vxcr, vxci)
         else
            call sigma (nul, ibp, rs, rscore, xk, vxcr, vxci)
         endif
         del1r = gsrel(i) * (vxcr - vxcrmu(i))

c        Correct local momentum according to the formula
c        p^2 = k^2 - 2*mu + kf^2 - 2*delta.  Note that imag part
c        of delta is ignored, since xk2 is a real quantity.

c        find xk(em) by iterative solution of dyson equation
  50     continue
         xk2 = 2*(dble(em) - xmu - del1r) + xf**2
         if (xk2 .lt. 0)  then
            write(slog,'(1pe13.5, 3i8, a)')
     1         xk2, i, ie, iph, ' xk2, i, ie, iph'
            call wlog(slog)
            call wlog(' em, xf**2, xmu, delta')
            write(slog,'(1p, 5e13.5)') dble(em), xf**2, xmu, del1r
            call wlog(slog)
            call par_stop('XCPOT-2')
         endif
         xk = sqrt (xk2)

c        calculate \delta_2 and \delta_v,2 with the corrected
c        local momentum
         call sigma (ixc, ibp, rs, rscore, xk, vxcr, vxci)
c        delta corrected calculated with new local momentum
         delr = gsrel(i) * (vxcr - vxcrmu(i))
         deli = vxci-vxcimu(i)

         if (ixc.ge.5 .and. i.eq.jri1 .and. xk.gt.xf) then
            if (ixc.eq.5 .or. ixc.eq.6) then
               delvr = delr
               delvi = deli
            endif
         endif

         if (niter.lt.nmax) then
            del1r=delr
            niter=niter+1
            go to 50
         endif

         if (ixc .ge. 5 .and. i.lt.jri1 .and. xk.gt.xf) then
            if (ixc.eq.5) then
               xkpp=sqrt(xk**2-xf**2+xfval**2)
               call sigma (ixc, ibp, rsval,rscore,xkpp,vxcvr,vxcvi)
            else
               call sigma (ixc, ibp, rs, rscore, xk, vxcvr, vxcvi)
            endif
            delvr = vxcvr-vvxcrm(i)
            delvi = vxcvi-vvxcim(i)
         endif

c        Josh - Skip SC loop if CSigma is called. CSigma calculates self consistently.
 15      continue

         delta = dcmplx(delr,deli)

c	 Josh - write out delta sigma at interstitial level to sigma.dat.
c         if(i.eq.jri1) then
c            write(45,'(X,20e14.6)') (DBLE(em) - xmu)*hart, delr*hart,
c     &                        deli*hart, DBLE(ZRnrm), DIMAG(ZRnrm),
c     &                        SQRT(DBLE(ZRnrm)**2+DIMAG(ZRnrm)**2),
c     &                        ATAN2(DIMAG(ZRnrm),DBLE(ZRnrm)),
c     &                        SQRT(DBLE(em-xmu)/2.d0)/ABS(deli)*bohr
c         end if
c	 Josh END

         if (ixc .eq. 5) delta = dcmplx(delr,delvi)
         v(i) = vtot(i) + delta
         if (ixc .ge. 5) then
            deltav = dcmplx(delvr,delvi)
            vval(i) = vvalgs(i) + deltav
         endif
         if(i.eq.jri1) then
            volume = 0.d0
         elseif(i.eq.jri) then
            volume = 4.d0*pi/3.d0*(rnrm**3 - exp(3.d0*((i-1)*dx - x0)))
         else
            volume = 4.d0*pi/3.d0*exp(3.d0*((i-1)*dx-x0))*(exp(3.d0*dx)
     &           - 1.d0)
         end if
         if(volume.lt.0.d0) volume = 0.d0
         omp = SQRT(3.d0/rs**3)
c         write(39,'(I5,20f30.10)') i, dble(em-xmu), volume/totvol,
c     &        exp((i-1)*dx - x0)*bohr,
c     &        rnrm*bohr, densty(i), denval(i), dble(volume*delta),
c     &        dimag(volume*delta), omp*hart, omp-ompm1
c         ompm1 = omp
c         delavg = delavg + volume*delta
 20   continue
c      write(39,*)

      ifirst = 1
      delavg = delavg/totvol
c      write(38,'(X,20e14.6)') (DBLE(em) - xmu)*hart, dble(delavg)*hart,
c     &     dimag(delavg)*hart,
c     &     SQRT(DBLE(em-xmu)/2.d0)/ABS(dimag(delavg))*bohr,
c     &     totvol
c     Reference the potential with respect to mt potential, ie,
c     first interstitial point.  v(jri1) = 0

c     Note that the reference does not contain the core hole lifetime
c     since the total atomic potential should have it. However in the
c     perturbation  deltav = v - vmt it cancels out.
c     ( deltav = vat - igamma - (vatmt-igamma) ).

 888  eref = v(jri1)
      do 910 i = 1, jri1
  910 v(i) = v(i) - eref
      if (ixc.ge.5) then
         do 920 i = 1, jri1
  920    vval(i) = vval(i) - eref
      else
         do 930 i = 1, jri1
  930    vval(i) = v(i)
      endif

c     Real self energy, zero imag part
      if (lreal.gt.0)  then
         do 950  i = 1, jri1
            v(i) = dble(v(i))
            if (ixc.gt.4)  vval(i) = dble(vval(i))
  950    continue
         eref = dble(eref)
      endif

      return
      end

      subroutine sigma (ixc, ibp, rs, rscore, xk, vr, vi)
      implicit double precision (a-h, o-z)

      if ((ixc.eq.0 .or. ixc.ge.5) .and. ibp .eq. 0) then
         call rhl (rs, xk, vr, vi)
      elseif ((ixc.eq.0.or. ixc.ge.5) .and. ibp .eq. 1) then
         call rhlbp (rs, xk, vr, vi)
      elseif (ixc .eq. 1) then
         vi = 0
         call edp(rs,xk,vr)
      elseif (ixc .eq. 3) then
         call edp(rs,xk,vr)
         call imhl (rs,xk,vi,icusp)
      endif

      if (ixc .ge. 6) then
         call edp(rscore,xk,vrp)
         vr = vr - vrp
      endif

      return
      end
