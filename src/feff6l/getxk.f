      double precision function getxk (e)
      implicit double precision (a-h, o-z)

c     Make xk from energy as
c          k =  sqrt( e)  for e > 0  (above the edge)
c          k = -sqrt(-e)  for e < 0  (below the edge)

      getxk = sqrt(abs(e))
      if (e .lt. 0)  getxk = - getxk
      return
      end
