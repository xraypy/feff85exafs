      subroutine ff2g (gtr, iph, ie, ilast, xrhoce, xrhole, xrhocp,
     1             ee, ep, yrhole, yrhoce,  yrhocp, rhoval, xnmues,
     2             xnatph, xntot, iflr, iflrp, fl, fr, iunf)
      implicit double precision (a-h, o-z)
c     the main output is l-dos in xrhoce, and valence density 
c     of states at distance r
c     in yrhoce, which at the input are only embedded atom quantities

      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'

      complex*16 xrhoce(0:lx,0:nphx), xrhocp(0:lx,0:nphx)
      complex*16 yrhoce(251), yrhocp(251)
      complex*16 xrhole(0:lx)
      complex*16 yrhole(251,0:lx)
      complex gtr(0:lx)
      complex*16 ee, ep, del, der, fl, fr
      complex*16 cchi(0:lx)
      dimension xnmues(0:lx), rhoval(251)

c     chi from fms is contained in gtr
      do j = 0, lx
         cchi(j) =  dble( real( gtr(j) )) + coni* dble(aimag( gtr(j) ))
      enddo

      do il = 0,lx
         xrhoce(il, iph)=xrhoce(il, iph)+
     $        cchi(il)*xrhole(il) 
         if (ie.eq.1) xrhocp(il,iph) = xrhoce(il,iph)
      enddo

      del = ee-ep
      der = del
c     if iflr=1 add/subtract integral from point to real axis
c     factor 2 below comes from spin degeneracy
      if (iflr.eq.1) der = der - coni * 2 * dimag(ee)
      if (iflrp.eq.1) del = del + coni * 2 * dimag(ep)
      do il = 0, lx
         if (il.le.2 .or. iunf.ne.0) then
            fl = fl + 2 * xrhocp(il,iph) * xnatph
            fr = fr + 2 * xrhoce(il,iph) * xnatph
            xnmues(il) = xnmues(il) + 
     $           dimag( xrhoce(il,iph) * der + xrhocp(il,iph) * del)
            xntot = xntot + xnmues(il) * xnatph
         endif
      enddo

cc    calculate r-dependent l-dos for later use
      do il = 0,lx
         do ir = 1,ilast
            if (il.le.2 .or. iunf.ne.0) then
               yrhoce(ir) = yrhoce(ir) + cchi(il)*yrhole(ir,il)
               if (ie.eq.1) yrhocp(ir) = yrhoce(ir)
            endif
         enddo
      enddo

      do ir = 1, ilast
         rhoval(ir) = rhoval(ir) + dimag(yrhoce(ir)*der+yrhocp(ir)*del)
      enddo

      return
      end
