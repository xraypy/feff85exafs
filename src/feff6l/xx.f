      double precision function xx (j)
      implicit double precision (a-h, o-z)
c     x grid point at index j, x = log(r), r=exp(x)
      parameter (delta = 0.050 000 000 000 000)
      parameter (c88   = 8.800 000 000 000 000)
c     xx = -8.8 + (j-1)*0.05
      xx = -c88 + (j-1)*delta
      return
      end

      double precision function rr(j)
      implicit double precision (a-h, o-z)
c     r grid point at index j
      rr = exp (xx(j))
      return
      end

      function ii(r)
      implicit double precision (a-h, o-z)
c     index of grid point immediately below postion r
      parameter (delta = 0.050 000 000 000 000)
      parameter (c88   = 8.800 000 000 000 000)
c     ii = (log(r) + 8.8) / 0.05 + 1
      ii = (log(r) + c88) / delta + 1
      return
      end
