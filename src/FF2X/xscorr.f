      subroutine xscorr(ispec, emxs ,ne1, ne, ik0, xsec, xsnorm, chia,
     1                  vrcorr, vicorr, cchi)
c     convolute xmu(E)=xsec+xsnorm*chia with lorentzian using 
c     calculations in the complex energy plane

c     Input: ispec - type of spectroscopy
c       emxs - complex energy grid
c       ne1 - number of points on horizonatal axis 
c       ne - total number of points (ne-ne1) points on vertical axis
c       ik0 - Fermi level index on horizontal axis
c       xsec, xsnorm, chia - give function f in complex energy plain
c           xmu(ie) = xsec + xsnorm*chia
c       vrcorr = correction for the shift of the Fermi level
c       vicorr = 0 (disabled)
c     Output: cchi(w) - result of convolution for w = dble(emxs)
c       cchi(w) = \int_C dE xmu(E)*xloss/pi/((E-w)**2+xloss**2) = 
c       xmu(w+i*xloss)* [1/2+atan(w-efermi/xloss)/pi] +
c       \int_C dE ff(E)*xloss/pi/((E-w)**2+xloss**2)
c       where ff(E)=xmu(E)-xmu(w+i*xloss) for w<efermi we use 
c       xmu(efermi+i*xloss) instead of xmu(w+i*xloss);
c       contour C starts at efermi, goes vertically to efermi+i*xloss 
c       and then goes horizontally to infinity + i*xloss

      implicit double precision (a-h, o-z)
      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'

      dimension  xsnorm(nex), omega(nex)
      complex*16 emxs(nex), xsec(nex), chia(nex), cchi(nex) 
      complex*16 xmu(nex), aa, bb, f1, f2, ff(nex), xmu0
c      complex*16 c1
      parameter (eps4 = 1.0d-4)
      complex*16 ec(nex), fc(nex), e1,e2, z1,z2, corr
      complex*16 lorenz
      external lorenz, astep

      ne2 = ne-ne1
      efermi = dble(emxs(ne)) 
      xloss = dimag(emxs(1))
      vicorr = 0

c     xmu - analytic function in complex energy plain
      do  ie = 1,ne
        xmu (ie) = xsec(ie) + xsnorm(ie)*chia(ie)
      enddo
c     real frequencies
      do ie = 1, ne1
        omega(ie) = dble(emxs(ie))
      enddo

      if (abs(vrcorr).gt.eps4) then
c       account for the fermi level shift
        bb = xmu(ik0)
        efermi = efermi - vrcorr
        call terpc(omega, xmu ,ne1, 1, efermi, bb)

c       shift the vertical axis
        do ie = 1, ne2
          emxs(ne1+ie) = emxs(ne1+ie) - vrcorr
        enddo

c       rescale values on vertical axis
        bb = bb/xmu(ik0)
        do ie = ne1+1, ne
          xmu(ie) = xmu (ie) * bb 
        enddo
      else
        bb = 1
      endif

c     construct the integration countur C
      nc = 0
c     start with points on vertical axis below xloss
      do ie = 1,ne2
        if (dimag(emxs(ne1+ie)).lt.xloss) then
          nc = nc+1
          ec(nc) = emxs(ne1+ie)
          fc(nc) = xmu(ne1+ie)
        endif
      enddo
c     add corner at efermi + xloss*i
      nc = nc+1
      ic0 = nc
      if (abs(vrcorr).gt.eps4) then
        ec(nc) = efermi + coni*xloss
        fc(nc) = bb * xmu(ik0)
      else
        ec(nc) = emxs(ik0)
        fc(nc) = xmu(ik0)
      endif
c     add points on horizontal axis above efermi
      if (ispec.ne.2) then
        do ie = 1,ne1
          if (dble(emxs(ie))-efermi.gt.eps4) then
            nc = nc+1
            ec(nc) = emxs(ie)
            fc(nc) = xmu(ie)
          endif
        enddo
      else
c       ispec=2 - emission calculations- need points below E_fermi
        do ie = ne1,1,-1
          if (efermi-dble(emxs(ie)).gt.eps4) then
            nc = nc+1
            ec(nc) = emxs(ie)
            fc(nc) = xmu(ie)
          endif
        enddo
      endif
c     endo of countour construction
              
c     cycle over frequency points 
      do ie = 1, ne1
        if (omega(ie).ge.efermi) then
          xmu0 = xmu(ie)
          if (ispec.eq.2) xmu0 = xmu(ik0)*bb
        else
          xmu0 = xmu(ik0)*bb
          if (ispec.eq.2) xmu0 = xmu(ie)
        endif
        e1 = omega(ie) + coni*xloss
        e2 = omega(ie) - coni*xloss
        do ic = 1, nc 
          ff(ic) = fc(ic) - xmu0
        enddo
        dele = omega(ie) - efermi
        cchi(ie) = xmu0 * astep( xloss, dele)
        if (ispec.eq.2) cchi(ie) = xmu0 - cchi(ie)
        corr = 0

        if (abs(dele).lt.eps4) dele = 0.0d0
        w1 = dimag(ec(1))
        w2 = dimag(ec(2))
        w3 = dimag(ec(3))
        ip =0

c       add half matsubara pole contribution
c       equivalent to integral from efermi to efermi+i*w1
        corr = corr + lorenz(xloss,w1,dele)*ff(1) *coni*w1
        if (nc0.gt.3) then
c       add sommerfeld correction (correction for derivative)
c         corr = corr + coni * w1**2 / 6   / (w3-w2) *
c    2   (lorenz(xloss,w3,dele)*ff(3)-lorenz(xloss,w2,dele)*ff(2))
        endif


c       cycle over contour points 
        do ic = 1,nc-1
c         perform integration over contour from efermi+i*2*w1 to 
c         efermi+i*xloss; linear interpolation of ff between  z1 and z2
          z1 = ec(ic)
          z2 = ec(ic+1)
c         if (ic.eq.1) z1 = efermi+coni*2*w1
          f1 = ff(ic)
          f2 = ff(ic+1)
c         if (ic.eq.1) f1 = (f1*(z2-z1) + f2*(z1-ec(ic))) / (z2-ec(ic))
c         add correction from pole above real axis
          aa = 0
          if (abs(z1-e1).gt.eps4 .and. abs(z2-e1).gt.eps4) then
            aa = log((z2-e1)/(z1-e1)) *(f1*(z2-e1)+f2*(e1-z1))
c           z1 or z2 equal to e1; in this case corr is exactly zero
          endif
c         second pole 
          aa = aa - log((z2-e2)/(z1-e2)) *(f1*(z2-e2)+f2*(e2-z1))
          corr = corr + aa/ (z2-z1) /2/pi/coni
        enddo
c       end of cycle over contour points
        if (ispec.eq.2) corr = -corr
c       if (ispec.eq.2) corr = 0

        cchi(ie) = cchi(ie) +  corr
c       return the result of convolution minus bare value
        cchi(ie) = cchi(ie) - xmu(ie)
      enddo
c     end of cycle over frequency points

c     restore the input energy mesh
      if (abs(vrcorr).gt.eps4) then
        do  ie = ne1+1, ne
          emxs(ie) = emxs(ie) + vrcorr
        enddo
      endif

      return
      end

      complex*16 function lorenz (xloss, w, dele)
      implicit double precision (a-h, o-z)
      include '../HEADERS/const.h'
c     ifp is dummy now. correspond to ifp=0 in old code
c     can remove it and change calls to lorenz in other routines

      lorenz = xloss /pi / (xloss**2+(coni*w-dele)**2)

      return
      end

      double precision function astep ( xloss, dele)
      implicit double precision (a-h, o-z)
      include '../HEADERS/const.h'

      astep = 0.5d0 + atan(dele/xloss) /pi
      if (astep.lt.0.d0) astep = 0.d0
      if (astep.gt.1.d0) astep = 1.d0

      return
      end
