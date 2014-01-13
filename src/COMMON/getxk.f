      double precision function getxk (e)
      implicit double precision (a-h, o-z)

c     Make xk from energy(in Hartrees) as
c          k =  sqrt(2*e)  for e > 0  (above the edge)
c          k = -sqrt(-2*e)  for e < 0  (below the edge)
      getxk = sqrt(abs(2*e))
      if (e .lt. 0)  getxk = - getxk
      return
      end
