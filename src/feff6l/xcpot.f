      subroutine xcpot (iph, ie, nr, index, ifirst, jri,
     1                  em, xmu, vi0, rs0, gamach,
     2                  vr, densty,
     3                  eref, v,
     4                  vxcrmu, vxcimu)

      implicit double precision (a-h, o-z)

c     INPUT
c     iph, ie used only for debug and labels.
c     nr          number of points in current Loucks r-grid
c     index       0  Hedin-Lunqvist + const real & imag part
c                 1  Dirac-Hara + const real & imag part
c                 2  ground state + const real & imag part
c                 3  Dirac-Hara + HL imag part + const real & imag part
c                 4  See rdinp for comment
c     ifirst      first entry flag, set to zero before first call for
c                 each unique potential, see vxcrmu and vxcimu below
c     jri         index of first interstitial point in current
c                 Loucks r grid
c     em          current energy grid point
c     xmu         fermi level
c     vi0         const imag part to subtract from potential
c     rs0         user input density cutoff, index=4 only
c     gamach      core hole lifetime
c     vr(nr)      total potential (coulomb and gs exchange corr)
c     densty(nr)  electron density
c
c     OUTPUT
c     eref        complex energy reference for current energy
c     v(nr)       complex potential including energy dep xc
c
c     WORKSPACE
c     vxcrmu and vxcimu are calculated only on first entry for a
c     particular unique potential, re-used on subsequent entries.
c     vxcrmu(nr)  real part of xc at fermi level
c     vxcimu(nr)  imag part of xc at fermi level
c
c     This subroutine uses atomic (hartree) units for energy,
c     phase uses rydbergs.  All inputs to and outputs from xcpot are
c     in rydbergs.  (Factor of 2 to convert from one to the other.)


      include 'const.h'
      include 'dim.h'

      dimension   vr(nr), densty(nr)
      complex*16  eref, v(nr)
      dimension   vxcrmu(nr), vxcimu(nr)
       character*128 messag
      complex*16  delta

c     First calculate vxc to correct the local momentum dispersion
c     relation, delta = vxc(e,k) - vxc(mu,k), and
c               p^2 = k^2 -mu + kf^2 - delta.
c     In jr theory, v(e,r) = vcoul(r) + vxc(e,r) =
c                          = vcoul(r) + vxcgs(r) + delta(e,r).

      if (index .eq. 2)  then
c        Ground state exchange, no self energy calculation
         do 10  i = 1, jri
            v(i) = vr(i)
   10    continue
      else
c        Add the self energy correction
         do 20  i = 1, jri
            rs = (3 / (4*pi*densty(i))) ** third
c           xf = 1.9191.../rs
            xf = fa / rs

c           xk2 is the local momentum squared, p^2 = k^2 - mu + kf^2,
c           k^2 represents energy measured from vacuum.
c           See formula 2.15 in Lee and Beni's paper with the last 2
c           terms neglected.  (complete reference?)
            xk2 = em + xf**2 - xmu

            if (xk2 .lt. 0)  then
               call echo(' error at XCPOT: xk2<0')
               call echo('  i, jri, rs, densty(i):')
               write(messag,'(2x,2i8,2g15.6)') i, jri, rs, densty(i)
               call echo(messag)
               call echo('  xf, fa, em, xmu, xk2:')
               write(messag,'(2x,5g15.6)')  xf, fa, em, xmu, xk2
               call echo(messag)
               call fstop(' at XCPOT-1')
            endif
            xk = sqrt(xk2)
            if (index .eq. 0)  call rhl(rs,xk,vxcr,vxci)
            if (index .eq. 1)  call edp(rs,xk,vi0,vxcr,vxci)
            if (index .eq. 3)  then
               call edp(rs,xk,vi0,vxcr,vxci)
               call imhl(rs,xk,vxci,icusp)
            elseif (index .eq. 4)  then
               rstmp = (1/rs**3 - 1/rs0**3) ** (-third)
               call edp(rstmp,xk,vi0,vxcr1,vxci1)
               call rhl(rs0,xk,vxcr2,vxci2)
               vxcr = vxcr1 + vxcr2
               vxci = vxci1 + vxci2
            endif

            if (ifirst .eq. 0)  then
c              vxc_mu indep of energy, calc only once
c              Calculate vxc at fermi level e = mu, j.m. 1/12/89
               xk = xf * 1.00001
               if (index .eq. 0) call rhl(rs,xk,vxcrmu(i),vxcimu(i))
               if (index .eq. 1) call edp(rs,xk,vi0,vxcrmu(i),vxcimu(i))
               if (index .eq. 3) then
                  call edp(rs,xk,vi0,vxcrmu(i),vxcimu(i))
                  call imhl (rs,xk,vxcimu(i),icusp)
               elseif (index .eq. 4)  then
                  rstmp = (1/rs**3 - 1/rs0**3) ** (-third)
                  call edp(rstmp,xk,vi0,vxcr1,vxci1)
                  call rhl(rs0,xk,vxcr2,vxci2)
                  vxcrmu(i) = vxcr1 + vxcr2
                  vxcimu(i) = vxci1 + vxci2
               endif
            endif

            delta = dcmplx (vxcr-vxcrmu(i), vxci-vxcimu(i))

c           Correct local momentum according to the formula
c           p^2 = k^2 - mu + kf^2 - delta.  Note that imag part
c           of delta is ignored, since xk2 is a real quantity.
            xk2 = em + xf**2 - xmu - delta
            if (xk2 .lt. 0)  then
               call echo(' error at XCPOT: xk2<0')
               call echo('   i, ie, iph, xk2, em, xf**2,'//
     $              ' xmu, delta') 
               write(messag,'(2x,3i7,5g15.6)')
     $              i, ie, iph, xk2, em, xf**2, xmu, delta
               call echo(messag)
               call fstop(' at XCPOT-2')
            endif
            xk = sqrt (xk2)

c           recalculate vxc(e,k) and vxc(mu,k) with the corrected
c           local momentum
            if (index .eq. 0)  call rhl(rs,xk,vxcr, vxci)
            if (index .eq. 1)  call edp(rs,xk,vi0,vxcr,vxci)
            if (index .eq. 3)  then
               call edp(rs,xk,vi0,vxcr,vxci)
               call imhl (rs,xk,vxci,icusp)
            elseif (index .eq. 4)  then
               rstmp = (1/rs**3 - 1/rs0**3) ** (-third)
               call edp(rstmp,xk,vi0,vxcr1,vxci1)
               call rhl(rs0,xk,vxcr2,vxci2)
               vxcr = vxcr1 + vxcr2
               vxci = vxci1 + vxci2
            endif

c           delta corrected calculated with new local momentum
            delta = dcmplx (vxcr-vxcrmu(i), vxci-vxcimu(i))

c           Note multiplication by 2 in the exchange correlation part to
c           to convert it to rydberg units.
   19       continue
            v(i) = vr(i) + 2*delta

   20    continue
      endif

c     Reference the potential with respect to mt potential, ie,
c     first interstitial point.  v(jri) = 0

c     Note that the reference does not contain the core hole lifetime
c     since the total atomic potential should have it. However in the
c     perturbation  deltav = v - vmt it cancels out.
c     ( deltav = vat - igamma - (vatmt-igamma) ).

      eref = v(jri)
      do 11  i = 1, jri
         v(i) = v(i) - eref
   11 continue

c     igamma added to the reference so that k^2 = E - Eref, where
c     Eref = Vat(mt) - igamma / 2
      eref = eref - coni * gamach / 2

c     Add const imag part
      eref = eref - coni * vi0

      ifirst = 1
      return
      end
