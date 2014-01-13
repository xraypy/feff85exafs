      subroutine cofcon (a,b,p,q)
c     acceleration of the convergence in the iterative process
c     b is the part of final iteration n is a function of the error (p)
c     (p) at iteration n and the error (q) at the iteration n-1.
c     if the product p*q is positive  b is increased by 0.1
c                        zero b is unchanged
c                        negative b is decreased by 0.1
c     b is between 0.1 and 0.9
c                a = 1. - b
c     ** at the end makes q=p
c
      implicit double precision (a-h,o-z)

      if (p*q)  11,31,21
 11   if (b .ge. 0.2) b = b - 0.1
      go to 31

 21   if (b .le. 0.8) b = b + 0.1

 31   a = 1.0 - b
      q=p
      return
      end
