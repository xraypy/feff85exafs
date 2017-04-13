      subroutine phase (iph, nr, dx, x0, ri, ne, em, edge,
     1                  index, rmt, xmu, vi0, rs0, gamach,
     2                  vtot, edens,
     3                  eref, ph, lmax)

      implicit double precision (a-h, o-z)

c     INPUT
c     iph          unique pot index (used for messages only)
c     nr, dx, x0, ri(nr)
c                  Loucks r-grid, ri=exp((i-1)*dx-x0)
c     ne, em(ne)   number of energy points, real energy grid
c     edge         energy for k=0 (note, edge=xmu-vr0)
c     index        0  Hedin-Lunqist + const real & imag part
c                  1  Dirac-Hara + const real & imag part
c                  2  ground state + const real & imag part
c                  3  Dirac-Hara + HL imag part + const real & imag part
c                  4, 5, 6, see rdinp or xcpot
c     rmt          r muffin tin
c     xmu          fermi level
c     vi0          const imag part to add to complex potential
c     rs0          user input density cutoff, used only with ixc=4
c     gamach       core hole lifetime
c     vtot(nr)     total potential, including gsxc
c     edens(nr)    density
c
c     OUTPUT
c     eref(ne)     complex energy reference including energy dep xc
c     ph(nex,ltot+1) complex scattering phase shifts
c     lmax         max l (lmax = kmax*rmt)

      include 'dim.h'
      parameter (bohr = 0.529 177 249, ryd  = 13.605 698)

      dimension   ri(nr), em(nex), vtot(nr), edens(nr)
      complex*16  eref(nex)
      complex*16  ph(nex,ltot+1)

c     work space for xcpot
      dimension   vxcrmu(nrptx), vxcimu(nrptx)
c     work space for fovrg
      complex*16 p(nrptx), q(nrptx), ps(nrptx), qs(nrptx), vm(nrptx)

      complex*16  p2, xkmt, temp, dny, pu, qu
      complex*16 jl(ltot+2), nl(ltot+2)
      complex*16 v(nrptx)
      external besjn

c     zero phase shifts (some may not be set below)
      do 100  ie = 1, ne
         do 90  il = 1, ltot+1
            ph(ie,il) = 0
   90    continue
  100 continue

c     limit l, lmax = kmax * rmt
c     lmax = rmt * sqrt(em(ne)-edge)
c     Use kmax = 20 so we get enough l-points even if kmax is small
      lmax = rmt * (20 * bohr)
      lmax = min (lmax, ltot)

c     set imt and jri (use general Loucks grid)
c     rmt is between imt and jri (see function ii(r) in file xx.f)
      imt = (log(rmt) + x0) / dx  +  1
      jri = imt+1
      if (jri .gt. nr)  call fstop(' at PHASE: jri > nr')
c     xmt is floating point version of imt, so that
c     rmt = (exp (x-1)*dx - x0).  xmt used in fovrg
      xmt = (log(rmt) + x0) / dx  +  1

      ifirst = 0
c     calculate phase shifts
      do 220 ie = 1, ne

         ihard = 0
         call xcpot (iph, ie, nr, index, ifirst, jri,
     1               em(ie), xmu, vi0, rs0, gamach,
     2               vtot, edens,
     3               eref(ie), v,
     4               vxcrmu, vxcimu)

c        fovrg needs v in form pot*r**2
         do 120  i = 1, jri
            v(i) = v(i) * ri(i)**2
  120    continue

c        p2 is (complex momentum)**2 referenced to energy dep xc
         p2 = em(ie) - eref(ie)
         xkmt = rmt * sqrt (p2)
         call besjn (xkmt, jl, nl)

         do 210  il = 1, lmax+1
            l = il - 1

            call fovrg (il, ihard, rmt, xmt, jri, p2, 
     1                  nr, dx, ri, v, dny,
     1                  pu, qu, p, q, ps, qs, vm)


            temp = (jl(il)*(dny-l) + xkmt*jl(il+1))  /
     1             (nl(il)*(dny-l) + xkmt*nl(il+1))
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
            beta = (xx**2 + (yy+1)**2) /
     1             (xx**2 + (yy-1)**2)
            beta = log(beta) / 4

            ph(ie,il) = dcmplx (alph, beta)

c           cut phaseshift calculation if they become too small
            if (abs(ph(ie,il)) .lt. 1.0e-6)  goto 220

  210    continue

  220 continue


c     Warn user if fovrg failed ihard test.
      if (ihard .gt. 0)  then
         call echo('  Hard tests failed in fovrg.')
         call echo(' Muffin-tin radius may be too large;'//
     1               ' coordination number too small.')
      endif

      return
      end
