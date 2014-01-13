      double precision function aprdev (a,b,l)
c     the result of this function is the coefficient for the term of 
c     power (l-1) for the product of two polynomes, whose coefficients
c     are in rows a and b 
 
      implicit double precision (a-h,o-z)
      dimension a(10),b(10)
 
      aprdev=0.0d 00
      do 11 m=1,l
 11      aprdev=aprdev+a(m)*b(l+1-m)
      return
      end
