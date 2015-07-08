      subroutine frnrm (rho, iz, rnrm)
      implicit double precision (a-h, o-z)
      include '../HEADERS/dim.h'
      dimension rho(nrptx)
      dimension xpc(251), ri(251)
c#mn
      external rr
      x=0

c     finds norman radius

c     Need overlapped densities.  We'll get them in the form
c     4*pi*density = rho.  Also need z of atom

c     Then integrate out to the point where the integral of
c     4*pi*density*r**2 is equal to iz
      sum= (9*rho(1)*rr(1)**3+28*rho(2)*rr(2)**3+23*rho(3)*rr(3)**3)/480
c     add initial point (r=0) correction (see subroutine somm2)
      dpas = 0.05d0
      d1 = 3.0d0
      dd=exp(dpas)-1.0d0
      db=d1*(d1+1.0d0)*dd*exp((d1-1.0d0)*dpas)
      db=rr(1)/db
      dd=rr(1)*(1.0d0+1.0d0/(dd*(d1+1.0d0)))/d1
      sum = sum + dd*rho(1)*rr(1)**2 - db*rho(2)*rr(2)**2

      fl = rho(4) *rr(4)**3
      fr = rho(5) *rr(5)**3
      frr = rho(6) *rr(6)**3
      sum = sum + (25*fl + 12 *fr -frr)/480
      do 10  i = 7, nrptx
         fll = fl
         fl = fr
         fr = frr
         frr = rho(i) * rr(i)**3
         sumsav = sum
         sum = sum + (13*(fr+fl) -fll -frr)/480
         if (sum .ge. iz)  then
            inrm = i-2
            x= (iz-sumsav)/(sum-sumsav)
            goto 20
         endif
   10 continue
      call wlog(' FRNRM Could not integrate enough charge to reach' //
     1          ' required z.')
      call par_stop('FRNRM-1')
   20 continue
      rnrm = rr(inrm)*(1 + x*0.05d0)
     
c     add next order correction ALA 3/97
        dx05 = 0.05d0
        x0 = 8.8d0
        jnrm =  int((log(rnrm) + x0) / dx05)  +  2
        i0=jnrm+1
        xirf = 2
        do 710 ir = 1, jnrm+2
           ri(ir) = rr(ir)
           xpc(ir) = rho(ir)*ri(ir)**2
  710   continue

        call somm2 (ri, xpc, dx05, xirf, rnrm,0,i0)
c       dq is how many new electrons are within norman sphere
        dn1 = xirf-iz
        x2 = x - dn1/((1-x)*xpc(inrm) + x*xpc(inrm+1))
        if (abs(x2-x).gt.0.0001d0) then
          xirf = 2
          rnrm = rr(inrm)*(1 + x2*0.05d0)
          call somm2 (ri, xpc, dx05, xirf, rnrm,0,i0)
          dn2 = xirf-iz
c         Newton-Raphson methof to find zeroes
          x = x2 - dn2 * (x2-x)/(dn2-dn1)
        endif
        rnrm = rr(inrm)*(1 + x*0.05d0)

      return
      end
