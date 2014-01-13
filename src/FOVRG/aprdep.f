      double precision function aprdep (a,b,l)
c     need to be in library for ATOM and PHASE; renamed aprdev
c     the result of this function is the coefficient for the term of 
c     power (l-1) for the product of two polynomes, whose coefficients
c     are in rows a and b 
 
      implicit double precision (a-h,o-z)
      dimension a(10),b(10)
 
      aprdep=0.0d 00
      do 11 m=1,l
 11      aprdep=aprdep+a(m)*b(l+1-m)
      return
      end
