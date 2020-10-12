      subroutine snlm (lmaxp1, mmaxp1, xnlm)
      implicit double precision(a-h,o-z)

c     Set nlm, legendre normalization factors, xnlm in common /nlm/
c     Calculates legendre norm factors
c     xnlm= sqrt ((2l+1)(l-m)!/(l+m)!)

      include '../HEADERS/dim.h'
c     include 'nlm.h'
      dimension xnlm(ltot+1,mtot+1)


c     flg(i) = i! * afac**i, set in factst
      dimension flg(0:210)

      call factst (afac, flg)

c     initialize xnlm explicitly
      do il = 1, ltot+1
         do im = 1, mtot+1
            xnlm(il,im) = 0
         enddo
      enddo

      do il = 1, lmaxp1
         mmxp1 = min (mmaxp1, il)
         do im = 1, mmxp1
            l = il-1
            m = im-1
            cnlm = (2*l+1) * flg(l-m) / flg(l+m)
            cnlm = sqrt(cnlm) * afac**m
            xnlm(il,im) = cnlm
         enddo
      enddo

      return
      end
      subroutine factst (afac, flg)
      implicit double precision (a-h,o-z)

c     FACTorial SeT, flg(i) = i! * afac**i
      dimension flg(0:210)

c     afac = 1/64 works with double precision on a VAX
      afac = 1./64.

      flg(0) = 1
      flg(1) = afac

      do i = 2, 210
         flg(i) = flg(i-1) * i * afac
      enddo
      return
      end
