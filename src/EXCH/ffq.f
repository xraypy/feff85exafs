      double precision function ffq (q, ef, xk, wp, alph)
      implicit double precision (a-h,o-z)

c     input:  q, wp, alph, ef, xk
c             q is dimensionless, normalized to fermi momentum
c             xk is momentum in invBohrs
c     output: ffq only

      wq = sqrt (wp**2 + alph*q**2 + q**4)
      ffq = (wp+wq)/(q**2) + alph/(2*wp)
      ffq = ((ef*wp) / (4*xk))  * log(ffq)

      return
      end
