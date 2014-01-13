      subroutine gglu( nsp, i0, ipi, ipf, lipotx, g0, tmatrx, g0t, gg)

      implicit real (a-h,o-z)
      implicit integer (i-n)
c  output
c    gg:  (nsp*lx**2, nsp*lx**2, 0:nphasx) submatrix spanning the entire
c          angular momentum basis for each unique potential

      include '../HEADERS/dim.h'
      include 'xparam.h'
      integer  i0 (0:nphx),  lipotx(0:nphx)

      parameter (one = 1, zero = 0)
      integer   ipiv(istatx)

c**** array of state kets at current energy
      common /stkets/ lrstat(4, istatx), istate

      complex   tmatrx(nspx, istatx)
c     big work matrices
      complex   g0( istatx, istatx), g0t( istatx, istatx)
      complex   g0s( istatx, nspx*(lx+1)**2)
c     return matrix containing info about each unique potential
      complex   gg(nspx*(lx+1)**2, nspx*(lx+1)**2, 0:nphasx)

      character*3  cerr
      character*13 trans

 400  format(i4)


c -------------------- LU gg
c     multiply T and G0 matrices together, construct g0t = 1 - G0*T
c     notice that the signs below for g0t ARE correct since 1 is the
c     unit matrix
c     since t is tri-diagonal, this product can be computed in n^2 time
c     also fill up some work matrices for use in eigenvalue and
c     determinant calculations and elsewhere
      do 320 icol = 1,istate
        do 310 irow = 1,istate
c         T diagonal contribution
          g0t(irow, icol) = - g0(irow, icol) * tmatrx(1, icol)
c         T off-diagonal contribution
          l1   = lrstat(2,icol)
          m1   = lrstat(3,icol)
          isp1 = lrstat(4,icol)
          m2 = m1+isp1
          if (nsp.eq.2 .and. m2.gt.-l1+1 .and. m2.lt.l1+2) then
             ist2 = icol + (-1)**isp1
             g0t(irow, icol) = g0t(irow, icol)
     1                   - g0(irow, ist2) * tmatrx(nsp, icol)
          endif
 310    continue

        g0t(icol, icol) = g0t(icol, icol) + one

 320  continue

c --- invert matrix by LU decomposition
c     call cgetrf from lapack.  this performs an LU decomposition on
c     the matrix g0t = 1 - g0*T
      call cgetrf( istate, istate, g0t, istatx, ipiv, info )
      if (info.lt.0) then
          call wlog('    *** Error in cgetrf when computing G')
          write(cerr,400)abs(info)
          call wlog('        Argument #'//cerr//
     $                ' had an illegal value.')
      elseif (info.gt.0) then
          call wlog('    *** Error in cgetrf when computing G')
          write(cerr,400)info
          call wlog('        g0t('//cerr// ','//cerr//
     $                ') is exactly 0 -- '//
     $                'this matrix cannot be decomposed.')
      endif

c     now we want g_c = (g0t)^-1 * g0.  Rather than calculating
c     the inverse of g0t from the LU decomposition, we can compute
c     g0t^-1 * g0 directly by backsubstituting the columns of G0.
c     See sect. 2.3 in Numerical Recipes or LAPACK Users' Guide
c     sect. 2.3

c     third arg in number of output columns, istate for full
c     matrix, ipart(ik) for just the parts of the matrix needed
c     to contruct fine structure + DOS functions

      do 620 ip=ipi, ipf
        ipart = nsp*(lipotx(ip)+1)**2
        do 590 is1 = 1, istate
        do 590 is2 = 1, ipart
          g0s(is1,is2) = g0(is1, is2 + i0(ip))
  590   continue

        trans = 'NotTransposed'
        call cgetrs(trans, istate, ipart, g0t, istatx,
     $                ipiv, g0s, istatx, info)
        if (info.lt.0) then
            call wlog('    *** Error in cgetrf')
            write(cerr,400) abs(info)
            call wlog('        Argument #'//cerr//
     $              ' had an invalid value.')
        endif

c **** at this point g0s contains the full MS ****

c  pack FMS matrix into an nsp*(lx+1)^2 x nsp*(lx+1)^2 matrix for each
c  ipot

        do 600 is2=1,ipart
        do 600 is1=1,ipart
          gg( is1, is2, ip) = g0s( is1+i0(ip), is2)
 600    continue
 620  continue

      return
      end
