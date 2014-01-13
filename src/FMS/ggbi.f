      subroutine ggbi( nsp, i0, ipi, ipf, lipotx, g0, tmatrx, g0t, gg,
     1                 toler1, toler2, lcalc, msord)

      implicit real (a-h,o-z)
      implicit integer (i-n)
c  output
c    gg:  (nsp*lx**2, nsp*lx**2, 0:nphasx) submatrix spanning the entire
c          angular momentum basis for each unique potential
c     BiCGStab algorithm: Saad, Iterative methods for ..., p. 220 (1996)

      include '../HEADERS/dim.h'
      include 'xparam.h'
      integer  i0 (0:nphx),  lipotx(0:nphx)

      parameter (one = 1, zero = 0)
      complex coni
      parameter (coni = (0,1))

c**** array of state kets at current energy
      common /stkets/ lrstat(4, istatx), istate

      complex   tmatrx(nspx, istatx)
c     big work matrices
      complex   g0( istatx, istatx), g0t( istatx, istatx)
      logical lcalc
      dimension lcalc(0:lx)

c     Lanczos method variables
      complex xvec( istatx), yvec(istatx), avec(istatx), asve(istatx)
      complex rvec(istatx), pvec(istatx), svec( istatx)
      complex aa, dd, aw, wa, ww
      complex del, delp, omega, chi, psi
c     return matrix containing info about each unique potential
      complex   gg(nspx*(lx+1)**2, nspx*(lx+1)**2, 0:nphasx)

c      notice that in gglu we invert (1-Gt), but here (1-tG).
c     multiply T and G0 matrices together, construct g0t = 1 - T*G0
c     notice that the signs below for g0t ARE correct since 1 is the
c     unit matrix
c     since t is tri-diagonal, this product can be computed in n^2 time
c     also fill up some work matrices for use in eigenvalue and
c     determinant calculations and elsewhere
c     cycle over dimensions of matrix g0t
      do 10 icol = 1,istatx
      do 10 irow = 1,istatx
 10   g0t(irow,icol) = 0

      do 30 icol = 1,istate
        do 20 irow = 1,istate
c         T diagonal contribution T(irow, irow)
          if ( abs( g0(irow, icol)) .gt. toler2 )
     1    g0t(irow,icol)=g0t(irow,icol) - tmatrx(1,irow) * g0(irow,icol) 

c         T off-diagonal contribution T(ist2, irow) in tmatr(2,irow)
c         T off-diagonal contribution T(irow, ist2) in tmatr(2,ist2)
          l1   = lrstat(2,irow)
          m1   = lrstat(3,irow)
          isp1 = lrstat(4,irow)
          m2 = m1+isp1
          if (nsp.eq.2 .and. m2.gt.-l1+1 .and. m2.lt.l1+2) then
c            spin-flip contribution
             ist2 = irow + (-1)**isp1
             if ( abs( g0(ist2, icol)) .gt. toler2)
     1       g0t(irow, icol) = g0t(irow, icol)
     2                   - tmatrx(nsp, ist2) * g0(ist2, icol) 
          endif
 20     continue

        g0t(icol, icol) = g0t(icol, icol) + one
 30   continue

      do 920 ip=ipi, ipf
        ipart = nsp*(lipotx(ip)+1)**2
        do 910 is1 = 1, ipart
          is2 = is1+i0(ip)
          l1   = lrstat(2,is2)
          if (.not.lcalc(l1)) goto 910

c         start first tier with xvec=0
          istart = -1
          msord = 0
          do 40 is = 1, istate
          avec(is) = 0
  40      xvec(is) = 0

c         RESTART here if necessary
  50      continue
          istart = istart+1

          if (istart.gt.0) call matvec( istatx,istate,g0t,xvec,avec,1)
          do 90 is = 1,istate
 90       rvec(is) = - avec(is)
c         rvec = bvec - g0t*xvec , in our case bvec(is) = delta_{is,is2}
          rvec(is2) = rvec(is2) + 1
cc          Check convergence criteria: |r_n+1| < tol
            ipass = 1
            do 92 is = 1, istate
              if ( abs(real(rvec(is))).gt.toler1) goto 93
              if ( abs(aimag(rvec(is))).gt.toler1) goto 93
 92         continue
            ipass = 0
 93         continue
            if (ipass.eq.0) goto 700

          do 95 is = 1, istate
 95       pvec(is) = rvec(is)
          call matvec( istatx,istate,g0t,pvec,avec,1)
          msord = msord + 1

c         choose yvec that del and delp close to one
          call cdot( istatx, istate, avec, avec, aa)
          call cdot( istatx, istate, rvec, avec, wa)
          aw = real(wa) - coni* aimag(wa)
          call cdot( istatx, istate, rvec, rvec, ww)
          dd = aa*ww - aw*wa
          if (abs(dd/aa/ww) .lt.1.e-8) then
            do 96 is = 1,istate
  96        yvec(is) = rvec(is) / ww
          else
            ww = ( ww - aw ) / dd
            aa = ( wa - aa) / dd
            do 97 is = 1,istate
  97        yvec(is) = rvec(is) * aa + avec(is) * ww
          endif
          call cdot( istatx, istate, yvec, rvec, del)

c         it seems ran out of precision for nit>150
          nitx = 30
          do 500 nit = 0, nitx
            call cdot( istatx, istate, yvec, avec, delp)
            omega = del / delp

            do 120 is = 1, istate
 120        svec(is) = rvec(is) - omega * avec(is)
cc          Check convergence criteria: |s_n+1| < tol
            ipass = 1
            do 122 is = 1, istate
              if ( abs(real(svec(is))).gt.toler1) goto 123
              if ( abs(aimag(svec(is))).gt.toler1) goto 123
 122        continue
            ipass = 0
 123        continue
            if (ipass.eq.0)  then
              do 124 is = 1, istate
 124          xvec(is) = xvec(is) + omega*pvec(is)
              goto 700
            endif

            call matvec( istatx,istate,g0t,svec,asve,1)
            msord = msord + 1
            call cdot( istatx, istate, asve, asve, aa)
            call cdot( istatx, istate, asve, svec, wa)
            chi = wa / aa
            do 125 is = 1, istate
 125        xvec(is) = xvec(is) + omega*pvec(is) + chi*svec(is)

            do 130 is = 1, istate
 130        rvec(is) = svec(is) - chi* asve(is)

cc          Check convergence criteria: |r_n+1| < tol
            ipass = 1
            do 370 is = 1, istate
              if ( abs(real(rvec(is))).gt.toler1) goto 380
              if ( abs(aimag(rvec(is))).gt.toler1) goto 380
 370        continue
            ipass = 0
 380        continue
            if (ipass.eq.0) goto 700

c           prepare for next iteration
            call cdot( istatx, istate, yvec, rvec, del)
            psi = del / (delp * chi)

            do 135 is = 1, istate
 135        pvec(is) = rvec(is) + psi * (pvec(is)-chi*avec(is))
            call matvec( istatx,istate,g0t,pvec,avec,1)
            msord = msord + 1

 500      continue
c         restart since ran out of iterations
          goto 50

c         exit if tolerance has been achieved
 700      continue
c         print*, ' BI iterations:', nit + istart*nitx
c         end of BI iterations

c         at this point xvec = (1-tG)**-1 * bvec  with chosen tolerance
c         pack FMS matrix into an nsp*(lx+1)^2 x nsp*(lx+1)^2 matrix 
c         for each ipot
          do 800 is2=1,ipart
            gg( is2, is1, ip) = zero
            do 790 is = 1,istate
              gg( is2, is1, ip) = gg( is2, is1, ip) +
     1        g0( is2+i0(ip), is) * xvec(is)
 790        continue
 800      continue

 910    continue
 920  continue

      return
      end
