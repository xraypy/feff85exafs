      subroutine rholsz ( dx, x0, ri, em,
     2                  ixc, rmt, rnrm,
     3                  vtot, vvalgs, xnval, dgcn, dpcn, eref,
     4                  adgc, adpc, xrhole, xrhoce, ph,
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


      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'

c     max number allowed in xsect r-grid
      parameter (nrx = nrptx)

c     output
      complex*16  xrhole(-4:3,-4:3)
      complex*16  xrhoce(-4:3, -4:3)
      complex*16  ph(lx+1)

      dimension ri(nrptx)
c      dimension ri05(251)
      dimension  vtot(nrptx), vvalgs(nrptx)
      complex*16 vtotc(nrptx), vvalc(nrptx)
      dimension xnval(30), dgcn(nrptx,30), dpcn(nrptx,30)
      dimension adgc(10,30), adpc(10,30)

c     energy grid in complex e-plane
      complex*16 em, eref

c     work space for dfovrg: regular and irregular solutions
      complex*16 pr(nrx,2,2), qr(nrx,2,2), pn(nrx,2,2), qn(nrx,2,2)

      complex*16  p2, xkmt, ck
c      complex*16 xck
      complex*16  pu, qu
      complex*16  xfnorm, xirf, xmp, xpm
      complex*16  temp,  phx, phm(2,2), factor

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

      do 10 i = -4,3
      do 10 j = -4,3
         xrhole(i,j) = 0
         xrhoce(i,j) = 0
  10  continue
      do 15 i=1,lx+1
  15  ph(i) = 0

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

      do 200 lll=0,lmax
        do 199 jd = 0,1
          ikap = (lll+jd)* (-1)**jd
          if (ikap.eq.0) goto 199

          ilp = lll + 1
          if (ikap.gt.0) ilp = lll - 1
          im = 1+ jd

          do 150 j = 1, 2
            ic3 = j-1
            if (lll.eq.0 .and. ic3.eq.1) goto 150

            irr = -1
            call dfovrg ( ncycle, ikap, rmt, ilast, jri, p2, dx,
     $                ri, vtotc, vvalc, dgcn, dpcn, adgc, adpc,
     $                xnval, pu, qu, pn(1,im,j), qn(1,im,j),
     $                iz, ihole, xion, iunf, irr, ic3)
            
            call exjlnl (xkmt, lll, jl, nl)
            call exjlnl (xkmt, ilp, jlp1, nlp1)
            call phamp (rmt, pu, qu, ck,  jl, nl, jlp1, nlp1, ikap,
     1                  phx, temp)
            if (lll.eq.0)  ph(1)=phx
            phm(im,j) = phx

c           Normalize final state  at rmt to
c           rmt*(jl*cos(delta) - nl*sin(delta))
            xfnorm = 1 / temp
c           normalize regular solution
            do 133  i = 1,ilast
              pr(i,im,j)=pn(i,im,j)*xfnorm
              qr(i,im,j)=qn(i,im,j)*xfnorm
  133       continue

c          find irregular solution
            irr = 1
            pu = ck*alphfs
            factor = pu/(1+sqrt(1+pu**2))
            if (ikap.lt.0) factor = -factor
c           set pu, qu - initial condition for irregular solution at ilast
c           qu=(nlp1*cos(phx)+jlp1*sin(phx))*pu *rmt
c           pu = (nl*cos(phx)+jl*sin(phx)) *rmt
            qu=(nlp1*cos(phx)+jlp1*sin(phx))* factor *rmt 
            pu = (nl*cos(phx)+jl*sin(phx)) *rmt 

            call dfovrg (ncycle, ikap, rmt, ilast, jri, p2, dx,
     1              ri, vtotc,vvalc, dgcn, dpcn, adgc, adpc,
     1              xnval, pu, qu, pn(1,im,j), qn(1,im,j),
     1              iz, ihole, xion, iunf, irr, ic3)
cc            set N- irregular solution , which is outside
cc            N=(nlp1*cos(ph0)+jlp1*sin(ph0))*factor *rmt * dum1
cc            N = i*R - H*exp(i*ph0)
              temp = exp(coni*phx)
              do i = 1, ilast
                pn(i,im,j) = coni * pr(i,im,j) - temp * pn(i,im,j)
                qn(i,im,j) = coni * qr(i,im,j) - temp * qn(i,im,j)
              enddo

 150      continue

c         combine all constant factors to temp
c         add relativistic correction to normaliz. and factor 2*lll+1
          temp = 2*ck / (1+factor**2) / pi
  
c         ic3 = 0, j= ic3+1
          j = 1
c         calculate diagonal radial integrals R(k1,k1) - xrhoce and xrhole
            do 190  i = 1, ilast
              xpc(i) = pr(i,im,j) **2 + qr(i,im,j) **2
 190        continue
            xirf = lll*2 + 2
c           i0 should be less or equal to  ilast
            i0=jnrm+1
            call csomm2 (ri, xpc, dx, xirf, rnrm, i0)
            xrhole(ikap,ikap) =xirf*temp*exp(coni*(phm(im,j)+phm(im,j)))

c         only central atom contribution needs irregular solution
            do 195  i = 1, ilast
              xpc(i) = pn(i,im,j)*pr(i,im,j)+ qn(i,im,j) *qr(i,im,j)
              xpc(i) = xpc(i) - coni*(pr(i,im,j)**2 + qr(i,im,j)**2)
 195        continue
            xirf =  1
            call csomm2 (ri, xpc, dx, xirf, rnrm, i0)
            xrhoce(ikap,ikap) = - xirf * temp

c         calculate cross terms
          if (ikap.lt.-1) then
            k1 = ikap + 2*lll + 1
            do 290  i = 1, ilast
              xpc(i) = pr(i,1,j) * pr(i,2,j) + qr(i,1,j) * qr(i,2,j) 
 290        continue
            xirf = lll*2 + 2
c           i0 should be less or equal to  ilast
            i0=jnrm+1
            call csomm2 (ri, xpc, dx, xirf, rnrm, i0)
            xrhole (ikap, k1) = xirf*temp* exp(coni*(phm(1,j)+phm(2,j)))
            xrhole (k1, ikap) = xrhole (ikap, k1)

c           ic3 = 1, j= ic3+1
            j = 2
            xpm =  exp(coni*(phm(1,j)-phm(2,j))) / 2
            xmp =  exp(coni*(phm(2,j)-phm(1,j))) / 2
            do 295  i = 1, ilast
              xpc(i) = (pn(i,1,j)*pr(i,2,j)+ qn(i,1,j) *qr(i,2,j)) * xmp
     1               + (pn(i,2,j)*pr(i,1,j)+ qn(i,2,j) *qr(i,1,j)) * xpm
              xpc(i) = xpc(i) - coni*(xpm+xmp) *
     1                 (pr(i,1,j)*pr(i,2,j) + qr(i,1,j)*qr(i,2,j))
 295        continue
            xirf =  1
            call csomm2 (ri, xpc, dx, xirf, rnrm, i0)
            xrhoce(ikap,k1) = - xirf * temp
            xrhoce(k1,ikap) =  xrhoce(ikap,k1)
          endif
 199    continue 
 200  continue 

c     calculate phase shift in old way (ic3=1) test new one
c     which is commented out above later
      do 300 lll = 1,lmax
          im = 1
          ikap = -lll-1
          irr = -1
          ic3 = 1
          call dfovrg ( ncycle, ikap, rmt, ilast, jri, p2, dx,
     $                ri, vtotc, vvalc, dgcn, dpcn, adgc, adpc,
     $                xnval, pu, qu, pr(1,im,1), qr(1,im,1),
     $                iz, ihole, xion, iunf, irr, ic3)
            
          call exjlnl (xkmt, lll, jl, nl)
          call exjlnl (xkmt, lll+1, jlp1, nlp1)
          call phamp (rmt, pu, qu, ck,  jl, nl, jlp1, nlp1, ikap,
     1                  phx, temp)
          ph(1+lll)=phx
 300  continue

      return
      end
