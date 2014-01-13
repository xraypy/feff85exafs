      subroutine gggm( nsp, i0, ipi, ipf, lipotx, g0, tmatrx, g0t, gg,
     1                 toler1, toler2, lcalc, msord)

      implicit real (a-h,o-z)
      implicit integer (i-n)
c  output
c    gg:  (nsp*lx**2, nsp*lx**2, 0:nphasx) submatrix spanning the entire
c          angular momentum basis for each unique potential
c     Lanczos algorithm: Graves-Morris,Salam, Num.Algor.21,p.213(1999)

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
      complex xvec( istatx), wvec(istatx), x0(istatx), x1(istatx)
      complex avec(istatx), bvec(istatx)
      complex r0(istatx), r1(istatx), t0(istatx), t1(istatx)
      complex aa, dd, aw, wa, ww, e0, e1, alpha, beta, theta, q0, q1
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
     1    g0t(irow,icol)=g0t(irow,icol) + tmatrx(1,irow) * g0(irow,icol) 

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
     2                   + tmatrx(nsp, ist2) * g0(ist2, icol) 
          endif
 20     continue

c       g0t(icol, icol) = g0t(icol, icol) + one
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
          bvec(is) = 0
  40      xvec(is) = 0
c         rvec = bvec - A*xvec , in our case bvec(is) = delta_{is,is2}
          bvec(is2) = 1

c         RESTART here if necessary
  50      continue
          istart = istart+1

          if (istart.gt.0) then
            do 60 is = 1, istate
  60        xvec(is) = xvec(is) + x0(is) / q0
            call matvec( istatx,istate,g0t,xvec,avec,1)
            do 70 is = 1, istate
  70        bvec(is) = avec(is) - xvec(is)
            bvec(is2) = bvec(is2) + 1
          endif
          do 80 is = 1,istate
 80       r0(is) = bvec(is)
          do 90 is = 1,istate
 90       x0(is) = 0
          do 95 is = 1, istate
 95       x1(is) = bvec(is)
          call matvec( istatx,istate,g0t,bvec,r1,1)
          msord = msord + 1

c         choose wvec that del and delp close to one
          call cdot( istatx, istate, r0, r0, ww)
          call cdot( istatx, istate, r1, r1, aa)
          call cdot( istatx, istate, r0, r1, wa)
          aw = real(wa) - coni* aimag(wa)
          dd = aa*ww - aw*wa
          if (abs(dd/aa/ww) .lt.1.e-8) then
            do 96 is = 1,istate
  96        wvec(is) = r0(is) / ww
          else
            ww = ( ww - aw ) / dd
            aa = ( wa - aa) / dd
            do 97 is = 1,istate
  97        wvec(is) = r0(is) * aa + r1(is) * ww
          endif
c         update dot products to avoid round off errors
          call cdot( istatx, istate, wvec, r0, e0)
          call cdot( istatx, istate, wvec, r1, e1)
          q0 = 1
          q1 = 1

c         it seems ran out of precision for nit>150
          nitx = 10
          do 500 nit = 1, nitx
            tol = toler1 * abs(q1) /10
cc          Check convergence criteria: |r1| < tol / 10
cc          so mostly code will not exit here
            ipass = 1
            do 98 is = 1, istate
              if ( abs(real(r1(is))).gt.tol) goto 99
              if ( abs(aimag(r1(is))).gt.tol) goto 99
  98        continue
            ipass = 0
  99        continue
            if (ipass.eq.0) then
              do 100 is = 1, istate
 100          xvec(is) = xvec(is) + x1(is) / q1
              goto 700
            endif

            alpha = e1 / e0
            do 130 is = 1, istate
 130        t0(is) = r1(is) - alpha* r0(is)
            call matvec( istatx,istate,g0t,t0,t1,1)
            msord = msord + 1

            call cdot( istatx, istate, t0, t1, wa)
            call cdot( istatx, istate, t0, t0, ww)
            call cdot( istatx, istate, t1, t1, aa)
            aw = real(wa) - coni* aimag(wa)
            theta = (wa - aa) / (ww - aw)

            do 145 is = 1, istate
 145        r0(is) = t1(is) - theta * t0(is)
            dd = 1- theta
            do 150 is = 1, istate
 150        x0(is) = t0(is) + dd * (x1(is) - alpha*x0(is))
            q0 = dd * (q1 - alpha*q0)
            tol = toler1 * abs(q0)

cc          Check convergence criteria: |r0| < tol
            ipass = 1
            do 370 is = 1, istate
              if ( abs(real(r0(is))).gt.tol) goto 380
              if ( abs(aimag(r0(is))).gt.tol) goto 380
 370        continue
            ipass = 0
 380        continue
            if (ipass.eq.0) then 
              do 390 is = 1, istate
 390          xvec(is) = xvec(is) + x0(is) / q0
              goto 700
            endif

c           prepare for next iteration
            call cdot( istatx, istate, wvec, r0, e0)
            beta = e0 / e1
            do 255 is = 1, istate
 255        t0(is) = r0(is) - beta * r1(is)
            call matvec( istatx,istate,g0t,t0,avec,1)
            msord = msord + 1
            dd = beta * theta
            do 260 is = 1, istate
 260        r1(is) = avec(is) + dd * r1(is)
            call cdot( istatx, istate, wvec, r1, e1)

            dd = beta * (1-theta)
            do 270 is = 1, istate
 270        x1(is) = x0(is) - dd * x1(is) + t0(is)
            q1 = q0 - (1-theta) * beta * q1
 500      continue
c         restart since ran out of iterations
          goto 50

c         exit if tolerance has been achieved
 700      continue
c         end of GM iterations

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
