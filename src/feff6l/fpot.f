      double precision function fpot (r,z,wa)
      implicit double precision (a-h,o-z)
      save
c
c thomas fermi potential at the point r; z=atomic number
c wa=number of electrons-z-1
c **********************************************************************
      wc=sqrt((r*(z+wa)**(1.0/3.0))/0.88530)
      wd=wc*(0.601120*wc+1.810610)+1.0
      we=wc*(wc*(wc*(wc*(0.04793*wc+0.21465)+0.77112)+1.39515)+1
     1 .81061)+1.0
      wc=(z+wa)*(wd/we)**2-wa
      fpot=-wc/r
      return
      end
