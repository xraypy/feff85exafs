      double precision function dalp (d1,d2,d3,d4)
      implicit double precision (a-h,o-z)
      save
c
c procedure of pratt to accelerate the convergence
c d1=initial (n-1);   d2=final (n-1);   d3=initial (n);   d4=final (n);
c **********************************************************************
      if ((d1+d4).eq.(d2+d3)) go to 10
      d=(d4-d2)/((d1+d4)-(d2+d3))
      if (d.lt.0.0) go to 20
      if (d.lt.0.5) go to 30
   10 d=0.5
      go to 30
   20 d=0.0
   30 dalp=d
      return
      end
