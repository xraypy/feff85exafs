      double precision function exchee (d,dr)
      implicit double precision (a-h,o-z)
      save
c jm if density= 0,make exchange energy equal to zero
      if (d .eq. 0.0) then
      exchee=0.0
      else
      x=(3.0*dr**2/d)**.333333333333/30.0
      rx=1.0/x
      exchee=.02520*(x**3*log(1.0+rx)+x*.50-x**2-1.0/3.0-0.2020129
     1 2*rx)
      endif
      return
      end
