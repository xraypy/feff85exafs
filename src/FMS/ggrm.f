      subroutine ggrm( nsp, i0, ipi, ipf, lipotx, g0, tmatrx, g0t, gg,
     1                 toler1, toler2, lcalc, msord)

      implicit real (a-h,o-z)
      implicit integer (i-n)
c  output
c    gg:  (nsp*lx**2, nsp*lx**2, 0:nphasx) submatrix spanning the entire
c          angular momentum basis for each unique potential

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

      complex xvec( istatx), xket( istatx), xbra( istatx)
      complex xketp(istatx), xbrap(istatx)
      complex zvec(istatx), rvec(istatx), svec(istatx)
      complex tket(istatx), tbra(istatx)
c      double precision  dum1, dum2
      complex alphac, betac, aa, bb, yy, aac, bbc, gamma
c      real alpha, beta
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
          msord=0
          do 40 is = 1, istate
          rvec(is) = 0
  40      xvec(is) = 0

c         RESTART here if necessary
  50      continue
          istart = istart+1

          if (istart.gt.0) call matvec( istatx,istate,g0t,xvec,rvec,1)
c         rvec = g0t*xvec - bvec, in our case bvec(is) = delta_{is,is2}
          rvec(is2) = rvec(is2) - 1
          do 90 is = 1,istate
 90       xket(is) = - rvec(is)
          call cdot( istatx, istate, xket, xket, bb)
          if (abs(bb).eq.0) goto 700

          xfnorm = 1.e0 / real(dble(bb))
          do 91 is = 1, istate
 91       xbra(is) = xket(is) * xfnorm
c         |t> = A |n> ; |n> - xket, |n-1> - xketp
          call matvec ( istatx, istate, g0t, xket, tket, 1)
          msord = msord + 1
          call cdot( istatx, istate, xbra, tket, aa)
          aac = real(aa) - coni*aimag(aa)
          bb = 0
          bbc= 0
          betac = aa
          yy = 1
c         initialize vectors
          do 110 is = 1,istate
            xketp(is) = 0
            xbrap(is) = 0
            zvec(is) = xket(is)
            xvec(is) = xvec(is) + zvec(is)/betac
 110      continue

          do 120 is = 1, istate
 120      svec(is) = tket(is)
          do 130 is = 1, istate
 130      rvec(is) = rvec(is) + svec(is) / betac

c         it seems ran out of precision for nit>150
          nitx = 100
          do 300 nit = 1, nitx
c           use recursion method to calculate a_n+1, b_n, |n+1>, <n+1|
            do 140 is = 1, istate
 140        tket(is) = tket(is) - aa*xket(is) - bb*xketp(is)
            call matvec ( istatx, istate, g0t, xbra, tbra, 2)
            do 150 is = 1, istate
 150        tbra(is) = tbra(is) - aac*xbra(is) - bbc*xbrap(is)
            call cdot( istatx, istate, tbra, tket, bb)
            if (abs(bb).eq.0) goto 700

            bb = sqrt (bb)
            bbc = real(bb) - coni*aimag(bb)
            do 160 is = 1, istate
              xketp(is) = xket(is)
              xbrap(is) = xbra(is)
 160        continue
            do 170 is = 1, istate
              xket(is) = tket(is) / bb
              xbra(is) = tbra(is) / bbc
 170        continue
            call matvec ( istatx, istate, g0t, xket, tket, 1)
            msord = msord + 1
            call cdot( istatx, istate, xbra, tket, aa)
            aac = real(aa) - coni*aimag(aa)
            
c           update iterative solution xvec, 
c           and residual rvec = g0t*xvec - |1>
            alphac = bb / betac
            do 210 is = 1, istate
 210        zvec(is) = xket(is) - alphac * zvec(is)
            do 220 is = 1, istate
 220        svec(is) = tket(is) - alphac * svec(is)

            betac = aa - alphac*bb
            yy = - alphac * yy
            gamma = yy / betac
            do 230 is = 1, istate
 230        xvec(is) = xvec(is) + gamma * zvec(is)
            do 240 is = 1, istate
 240        rvec(is) = rvec(is) + gamma * svec(is)

cc          Check convergence criteria: | rvec | < tol
c           call vecvec( istatx, istate, rvec, rvec, dum2)
c           if (dum2.le.tol) goto 700
cc          Check convergence criteria: | rvec | < tol
            ipass = 1
            do 250 is = 1, istate
              if ( abs(real(rvec(is))).gt.toler1) goto 260
              if ( abs(aimag(rvec(is))).gt.toler1) goto 260
 250        continue
            ipass = 0
 260        continue
            if (ipass.eq.0) goto 700

 300      continue
c         restart since ran out of iterations
          goto 50

c         exit if tolerance has been achieved
 700      continue
c         end of RM iterations

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

      subroutine cdot ( istatx, istate, abra, aket, cc)
c     dot product of two vectors
c     notice that we keep bra  vector as it's complex conjugate
c     thus need to conjugate abra here
      implicit real (a-h,o-z)
      implicit integer (i-n)
      complex coni
      parameter (coni = (0,1))
      complex abra, aket, cc, aa
      dimension abra(istatx), aket(istatx)

      cc = 0
      do 10 is = 1,istate
        aa = real(abra(is)) - coni*aimag(abra(is))
        cc = cc + aa * aket(is)
 10   continue
      return
      end

      subroutine vecvec ( istatx, istate, avec, bvec, cc)
c     dot product of two vectors
      implicit real (a-h,o-z)
      implicit integer (i-n)
      complex avec, bvec
      double precision cc, aa, bb
      dimension avec(istatx), bvec(istatx)

      cc = 0
      do 10 is = 1,istate
        aa = dble(real(avec(is))) * dble(real(bvec(is)))
        bb = dble(aimag(avec(is))) * dble(aimag(bvec(is)))
        cc = cc + aa + bb
 10   continue
      return
      end

      subroutine matvec (istatx, istate, amat, bvec, cvec, itrans)
c     itrans = 1  cvec = amat * bvec
c     itrans = 2  cvec = amat^+ * bvec
c     itrans = 3  cvec = amat^T * bvec
      implicit real (a-h,o-z)
      implicit integer (i-n)
      complex coni, aa
      parameter (coni = (0,1))
      complex amat, bvec, cvec
      dimension amat(istatx, istatx), bvec(istatx), cvec(istatx)

c     initialize cvec
      do 10 is = 1,istatx
 10   cvec(is) = 0

c     cycle over dimensions of amat
      do 20 icol = 1,istate
      do 20 irow = 1,istate
        if (itrans.eq.1) then
          cvec(irow) = cvec(irow) + amat(irow, icol) * bvec(icol)
        elseif(itrans.eq.2) then
          aa = real(amat(irow, icol)) - coni*aimag(amat(irow, icol))
          cvec(icol) = cvec(icol) + aa * bvec(irow)
        else
          cvec(icol) = cvec(icol) + amat(irow, icol) * bvec(irow)
        endif
 20   continue

      return
      end
