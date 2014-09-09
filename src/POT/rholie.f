      subroutine rholie ( ri05, nr05, dx, x0, ri, em,
     2                  ixc, rmt, rnrm,
     3                  vtot, vvalgs, xnval, dgcn, dpcn, eref,
     4                  adgc, adpc, xrhole, xrhoce, yrhole, yrhoce, ph,
     i                  iz, xion, iunf, ihole, lmaxsc)

      implicit double precision (a-h, o-z)

c     INPUT
c     dx, x0, ri(nr)
c                  Loucks r-grid, ri=exp((i-1)*dx-x0)
c     ne, em(ne)   number of energy points,  complex energy grid
c     ixc          0  Hedin-Lunqist + const real & imag part
c                  1  Dirac-Hara + const real & imag part
c                  2  ground state + const real & imag part
c                  3  Dirac-Hara + HL imag part + const real & imag part
c                  5  Dirac-Fock exchange with core electrons +
c                     ixc=0 for valence electron density
c     rmt          r muffin tin
c     rnrm         r norman
c     vtot(nr)     total potential, including gsxc, final state
c     dgcn(dpcn)   large (small) dirac components for central atom
c     adgc(adpc)   their development coefficients
c
c     OUTPUT
c     xrhole(0:lx)  integral over r of density function
c     xrhoce(0:lx)  the same integral for embedded atom only
c     yrhole(251,0:lx)   density function
c     yrhoce(251)        density function for embedded atom


      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'

c     max number allowed in xsect r-grid
      parameter (nrx = nrptx)

c     output
      complex*16  xrhole(0:lx)
      complex*16  xrhoce(0:lx)
      complex*16  yrhole(251,0:lx), yrhoce(251)
      complex*16  ph(lx+1)

      dimension ri(nrptx), ri05(251)
      dimension  vtot(nrptx), vvalgs(nrptx)
      complex*16 vtotc(nrptx), vvalc(nrptx)
      dimension xnval(30), dgcn(nrptx,30), dpcn(nrptx,30)
      dimension adgc(10,30), adpc(10,30)

c     energy grid in complex e-plane
      complex*16 em, eref

c     work space for dfovrg: regular and irregular solutions
      complex*16 pr(nrx), qr(nrx), pn(nrx), qn(nrx)

      complex*16  p2, xkmt, ck
c      complex*16  xck
      complex*16  pu, qu
      complex*16  xfnorm, xirf
      complex*16  temp,  phx, tempc

      complex*16 jl,jlp1,nl,nlp1
      complex*16  xpc(nrx)

c     initialize
      lmax=lmaxsc
      if (lmax.gt.lx) lmax = lx
      if (iz.le.4) lmax=2
      if (iz.le.2) lmax=1
      do 20 i = 1, nrptx
         vtotc(i)=vtot(i)
         vvalc(i)= vvalgs(i)
  20  continue
c     set imt and jri (use general Loucks grid)
c     rmt is between imt and jri (see function ii(r) in file xx.f)
      imt  = int((log(rmt) + x0) / dx)  +  1
      jri  = imt+1
      if (jri .gt. nrptx)  call par_stop('jri .gt. nrptx in phase')
      inrm = int((log(rnrm) + x0) / dx)  +  1
      jnrm = inrm+1

c     set limits for tabulations
      nr05= int((log(rnrm) + x0) / 0.05d0) + 5
      if (nr05.gt.251) nr05 = 251
c     ilast is the last integration point
c     it is larger than jnrm for better interpolations
      ilast = nint( (nr05-1) *0.05d0 / dx ) + 1
      if (ilast.gt.nrptx) ilast=nrptx

      do 10 lll = 0, lx
      do 10 j = 1, 251
         yrhole(j,lll) = 0
  10  continue
      do 30 j = 1, 251
  30  yrhoce(j) = 0

c     p2 is 0.5*(complex momentum)**2 referenced to energy dep xc
c     need hartree units for dfovrg
      p2 = em - eref
      if (mod(ixc,10) .lt. 5) then
        ncycle = 0
      else
        ncycle = 3
      endif
      ck = sqrt(2*p2 + (p2*alphfs)**2)
      xkmt = rmt * ck

      do 200 lll=0,lx
        if (lll.gt.lmax) then
           ph(lll+1) = 0
           xrhoce(lll) = 0
           xrhole(lll) = 0
           do 110 i = 1,251
  110      yrhole(i,lll) = 0
           goto 200
        endif

c       may want to use ihole=0 for new screening. 
c       don't want ro use it now
c       ihole = 0
        ikap = -1-lll
        irr = -1
        ic3 = 1
        if (lll.eq.0) ic3 = 0
        call dfovrg ( ncycle, ikap, rmt, ilast, jri, p2, dx,
     $                ri, vtotc, vvalc, dgcn, dpcn, adgc, adpc,
     $                xnval, pu, qu, pn, qn,
     $                iz, ihole, xion, iunf, irr, ic3)
        call exjlnl (xkmt, lll, jl, nl)
        call exjlnl (xkmt, lll+1, jlp1, nlp1)
        call phamp (rmt, pu, qu, ck,  jl, nl, jlp1, nlp1, ikap,
     1                  phx, temp)
        ph(lll+1)=phx

c     Normalize final state  at rmt to
c     rmt*(jl*cos(delta) - nl*sin(delta))
        xfnorm = 1 / temp
c     normalize regular solution
        do 133  i = 1,ilast
          pr(i)=pn(i)*xfnorm
          qr(i)=qn(i)*xfnorm
  133   continue

c      find irregular solution
        irr = 1
        pu = ck*alphfs
        pu = - pu/(1+sqrt(1+pu**2))
c       set pu, qu - initial condition for irregular solution at ilast
c       qu=(nlp1*cos(phx)+jlp1*sin(phx))*pu *rmt
c       pu = (nl*cos(phx)+jl*sin(phx)) *rmt
        qu=(nlp1*cos(phx)+jlp1*sin(phx))*pu *rmt 
        pu = (nl*cos(phx)+jl*sin(phx)) *rmt 

        call dfovrg (ncycle, ikap, rmt, ilast, jri, p2, dx,
     1              ri, vtotc,vvalc, dgcn, dpcn, adgc, adpc,
     1              xnval, pu, qu, pn, qn,
     1              iz, ihole, xion, iunf, irr, ic3)
cc      set N- irregular solution , which is outside
cc      N=(nlp1*cos(ph0)+jlp1*sin(ph0))*factor *rmt * dum1
cc      N = i*R - H*exp(i*ph0)
        temp = exp(coni*phx)
c       calculate wronskian
        qu = 2 * alpinv * temp * ( pn(jri)*qr(jri) - pr(jri)*qn(jri) )
        qu = 1 /qu / ck
c       qu should be close to 1
        do i = 1, ilast
          pn(i) = coni * pr(i) - temp * pn(i)*qu
          qn(i) = coni * qr(i) - temp * qn(i)*qu
        enddo

c     ATOM,  dgc0 is large component, ground state hole orbital
c     .      dpc0 is small component, ground state hole orbital
c     FOVRG, p    is large component, final state photo electron
c     .      q    is small component, final state photo electron

            
c    combine all constant factors to temp
c    add relativistic correction to normalization and factor 2*lll+1
        pu = ck*alphfs
        pu = - pu/(1+sqrt(1+pu**2))
        temp = (2*lll+1.0d0)/(1+pu**2) /pi *ck * 2
c    also scale by appropriate step in complex energy
        do 190  i = 1, ilast
          xpc(i) = pr(i) * pr(i) + qr(i) * qr(i) 
 190    continue
          
        do 191 ir=1,nr05
           call terpc(ri, xpc, ilast, 3, ri05(ir), tempc)
           tempc = tempc * temp
           yrhole(ir,lll)= tempc
 191    continue

        xirf = lll*2 + 2
c       i0 should be less or equal to  ilast
        i0=jnrm+1
        call csomm2 (ri, xpc, dx, xirf, rnrm, i0)
c       print out xirf for Bruce
        xrhole(lll) = xirf*temp

c     only central atom contribution needs irregular solution
        do 195  i = 1, ilast
          xpc(i) = pn(i)*pr(i)-coni*pr(i)*pr(i)
     1           + qn(i)*qr(i)-coni*qr(i)*qr(i)
c         yrhoce(i)=yrhoce(i) - temp*xpc(i)
 195    continue
        do 196 ir=1,nr05
           call terpc(ri, xpc, ilast, 3, ri05(ir), tempc)
           yrhoce(ir)=yrhoce(ir) - temp*tempc
 196    continue

        xirf =  1
        call csomm2 (ri, xpc, dx, xirf, rnrm, i0)
        xrhoce(lll) =  - xirf* temp
 200  continue 

      return
      end
