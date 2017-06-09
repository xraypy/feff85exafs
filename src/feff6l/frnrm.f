      subroutine frnrm (rho, iz, rnrm)
      implicit double precision (a-h, o-z)
      include 'dim.h'
      dimension rho(nrptx)

c     finds norman radius

c     Need overlapped densities.  We'll get them in the form
c     4*pi*density = rho.  Also need z of atom

c     Then integrate out to the point where the integral of
c     4*pi*density*r**2 is equal to iz
      sum = 0.0
      do 10  i = 1, nrptx-1
         fr = rho(i+1) * rr(i+1)**3
         fl = rho(i)   * rr(i)**3
         sumsav = sum
         sum = sum + 0.025*(fr+fl)
         if (sum .ge. iz)  then
            inrm = i+1
            goto 20
         endif
   10 continue
       call echo( ' FRNRM Could not integrate enough charge')
       call echo( '       to reach required z.')
       call fstop(' at FRNRM-1')
   20 continue
c     inrm is too big, subtract one from irnm and interpolate
c     to get correct value
      inrm = inrm - 1
      deltaq = iz - sumsav
      fr = rho(inrm+1) * rr(inrm+1)**3
      fl = rho(inrm)   * rr(inrm)**3
c     dipas is delta i * 0.05
      dipas = 2*deltaq / (fl + fr)
      rnrm = rr(inrm)*(1 + dipas)

      return
      end
