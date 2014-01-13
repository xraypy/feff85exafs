      double precision function dentfa (dr,dz,ch)
c     analitical approximation of potential is created for electrons in
c     thomas-fermi model for atom or free ion. dr distance from nucleus
c     with charge dz  
c     ch=ionicity = number of electrons-dz-1
      implicit double precision (a-h,o-z)
 
      dentfa=0.0d 00
      if ((dz+ch).lt.1.0d-04) return
      w=dr*(dz+ch)**(1./3.)
      w=sqrt(w/0.8853)
      t=w*(0.60112*w+1.81061)+1.
      w=w*(w*(w*(w*(0.04793*w+0.21465)+0.77112)+1.39515)+1.81061)+1
      dentfa=(dz+ch)*(1.0d 00-(t/w)**2)/dr
      return
      end
