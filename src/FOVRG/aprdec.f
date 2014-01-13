      complex*16  function aprdec(ala,bla,lla)
c     the result of this function is the coefficient for the term of
c     power (l-1) for the product of two polynomes, whose coefficients
c     are in rows a and b
 
      implicit double precision (a-h, o-z)
      complex*16 ala (10)
      integer lla
      dimension bla(10)
 
      aprdec = (0.0d0, 0.0d0)
      do 11 m = 1, lla
 11      aprdec = aprdec + ala(m) * bla(lla+1-m)
      return
      end
