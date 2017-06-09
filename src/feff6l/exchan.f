      double precision function exchan (d,dr,dexv)
      implicit double precision (a-h,o-z)
      save
c  dexv=0.0, hedin-barth corr. and exch. potential
c  dexv.ne. 0.0, dexv*slater exchange potential
c  d=4pi*rho*r^2 , radial density for r=dr
c  this function calculates exch=-r*Vexch
c  105.27578=32*(pi^2)/3
c  comments added by j. mustre 8/27/87
      if (dexv.eq.0.0) go to 10
      exchan=3.0*dexv*((dr*d/105.27578)**(1.0/3.0))
      return
   10 continue
      rrs=(d/(3.0*dr**2))**.33333333333
      exchan=+0.5*(1.22177412*rrs+.0504*log(30.0*rrs+1.0))*dr
      return
      end
