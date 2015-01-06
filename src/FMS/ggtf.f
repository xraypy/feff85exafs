      subroutine ggtf( nsp, i0, ipi, ipf, lipotx, g0, tmatrx, g0t, gg,
     1                 toler1, toler2, lcalc, msord)

      implicit real (a-h,o-z)
      implicit integer (i-n)
c  output
c    gg:  (nsp*lx**2, nsp*lx**2, 0:nphasx) submatrix spanning the entire
c          angular momentum basis for each unique potential
c     TFQMR: Saad, Iterative Methods for Sparse Matrices, p.225 (1996).

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

      complex xvec(istatx), uvec(istatx), avec(istatx), wvec(istatx)
      complex dvec(istatx), rvec(istatx), vvec(istatx)
      complex alpha, beta, aa, rho, eta
      real tau, nu, cm, err
c     return matrix containing info about each unique potential
      complex   gg(nspx*(lx+1)**2, nspx*(lx+1)**2, 0:nphasx)

      alpha = (1., 0.)

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
          rvec(is) = 0
          avec(is) = 0
  40      xvec(is) = 0

c         RESTART here if necessary
  50      continue
          istart = istart+1

          if (istart.gt.0) call matvec( istatx,istate,g0t,xvec,avec,1)
          do 90 is = 1,istate
 90       uvec(is) = - avec(is)
c         uvec = bvec - g0t*xvec , in our case bvec(is) = delta_{is,is2}
          uvec(is2) = uvec(is2) + 1
          call matvec( istatx,istate,g0t,uvec,avec,1)
          msord = msord + 1
          do 95 is = 1, istate
 95       wvec(is) = uvec(is)
          do 96 is = 1, istate
 96       vvec(is) = avec(is)
          do 97 is = 1, istate
 97       dvec(is) = 0
          call cdot( istatx, istate, uvec, uvec, aa)
          tau = sqrt(real(aa))
          nu = 0
          eta = 0
c         choose rvec = uvec /aa so that dot products are about 1
          do 98 is = 1, istate
 98       rvec(is) = uvec(is) / aa
          rho = 1

c         it seems ran out of precision for nit>150
          nitx = 20
          do 300 nit = 0, nitx
            if (mod(nit,2).eq.0) then
              call cdot( istatx, istate, rvec, vvec, aa)
              alpha = rho / aa
            else
              call matvec( istatx,istate,g0t,uvec,avec,1)
              msord = msord + 1
            endif

            do 115 is = 1, istate
 115        wvec(is) = wvec(is) - alpha * avec(is)

            aa = nu**2 * eta / alpha
            do 120 is = 1, istate
 120        dvec(is) = uvec(is) + aa * dvec(is)

            call cdot( istatx, istate, wvec, wvec, aa)
            nu = sqrt(real(aa)) / tau
            cm = 1 / sqrt(1+nu**2)
            tau = tau * nu * cm
            eta = cm**2 * alpha
            do 140 is = 1, istate
 140        xvec(is) = xvec(is) + eta * dvec(is)

cc          Check convergence criteria: | rvec | < tol
            err = (1.e0 + nit) / istate
            err = tau * sqrt(err) * 10
            if ( abs(err).lt.toler1) goto 700

            if (mod(nit,2) .ne.0) then
              aa = rho
              call cdot( istatx, istate, rvec, wvec, rho)
              beta = rho / aa
              do 210 is = 1, istate
 210          uvec(is) = wvec(is) + beta * uvec(is)

              do 215 is = 1, istate
 215          vvec(is) = beta * ( avec(is) + beta * vvec(is))
              call matvec( istatx,istate,g0t,uvec,avec,1)
              msord = msord + 1
              do 220 is = 1, istate
 220          vvec(is) = avec(is) + vvec(is)
            else
              do 230 is = 1, istate
 230          uvec(is) = uvec(is) - alpha * vvec(is)
            endif
 300      continue
c         restart since ran out of iterations
          goto 50

c         exit if tolerance has been achieved
 700      continue
c         end of TFQMR iterations

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
