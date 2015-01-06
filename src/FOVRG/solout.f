c      subroutine solout(en, fl, agi, api, kap, rmt,
c     1                  jri, max0, ic3, vm, iwkb)
      subroutine solout(en, fl, agi, api, kap,
     1                  jri, max0, ic3, vm, iwkb)
c                  resolution of the dirac equation
c                   p' - kap*p/r = - ( en/cl-v )*g - eg/r
c                   g' + kap*g/r = ( 2*cl+en/cl-v )*p + ep/r
c at the origin v approximately is -z/(r*cl) due to the point nucleus
c en one-electron energy in atomic units and negative
c fl power of the first term in development at the origin
c agi (api) initial values of the first development coefficient
c at the origin of the large(small)component
c kap quantum number kappa
c max0 the last point of tabulation of the wave function
 
      implicit double precision (a-h,o-z)
      include '../HEADERS/dim.h'
      include '../HEADERS/const.h'
      parameter (npi=6, test=1.0d+5)
      parameter (ccl=2*alpinv, csq=ccl**2 )
      complex*16 en,agi,api
      complex*16 gg,ag,gp,ap,dv,av,eg,ceg,ep,cep, vm(nrptx)
      common/comdic/cl,dz,gg(nrptx),ag(10),gp(nrptx),ap(10),dv(nrptx),
     1   av(10),eg(nrptx),ceg(10),ep(nrptx),cep(10)

      complex*16 ec,eph,f,g

c gg,gp -output, dv,eg,ep - input
c
c cl speed of light (approximately 137.037 in atomic units)
c dz nuclear charge
c gg (gp) large (small) component
c dv direct potential (v)     eg and ep exchange potentials
c ag,ap,av,ceg and cep are respectively the
c development coefficients for gg,gp,dv,eg and ep
c
      common/tabtec/hx,dr(nrptx),test1,test2,ndor,np,nes,method,idim
c hx exponential step
c dr radial mesh
c test1,test2,nes,method are dummy.
c  ndor number of terms for the developments at the origin
c np maximum number of the tabulation points
c idim dimension of the block dr


c{#mn: g77 chokes on taking real of a double-precision complex
c      if (real(av(1)).lt.0.0d 00.and.kap.gt.0) api=-agi*(kap+fl)/av(1)
c      if (real(av(1)).lt.0.0d 00.and.kap.lt.0) api=-agi*av(1)/(kap-fl)
      if (dble(av(1)).lt.0.0d 00.and.kap.gt.0) api=-agi*(kap+fl)/av(1)
      if (dble(av(1)).lt.0.0d 00.and.kap.lt.0) api=-agi*av(1)/(kap-fl)
c#mn}
      ec=en/cl
      ag(1)=agi
      ap(1)=api
      do 115 i=2,ndor
         ag(i)=ceg(i-1)
 115     ap(i)=cep(i-1)
c     integration of the inhomogenious system
c     no need in normalization, since we can use 
c     normalization agi=ag(1)=const
 
c            solution of the inhomogenios dirac equation
c gg gp initially exch. terms, at the time of return are wave functions
c ag and ap development coefficients of  gg and gp
c en one-electron energy
c fl power of the first development term at the origin
c agi (api) initial values of the first development coefficients
c at the origin of a large (small) component
 
c     initial values for the outward integration
      if (ic3.eq.0) then
c       Desclaux power expansion
         do 35 j=2,ndor
            k=j-1
            a=fl+kap+k
            b=fl-kap+k
            eph=a*b+av(1)*av(1)
            f=(ec+ccl)*ap(k)+ap(j)
            g=ec*ag(k)+ag(j)
            do 31 i=1,k
               f=f-av(i+1)*ap(j-i)
 31            g=g-av(i+1)*ag(j-i)
 
            ag(j)=(b*f+av(1)*g)/eph
 35         ap(j)=(av(1)*f-a*g)/eph

         do  41 i = 1,1
            gg(i)=0.0d 00
            gp(i)=0.0d 00
         do 41 j=1,ndor
            a=fl+j-1
            b=dr(i)**a
            gg(i)=gg(i)+b*ag(j)
 41         gp(i)=gp(i)+b*ap(j)
      else
c        see fovrg.f in feff6, be aware of different units
         twoz = -dble(av(1)) * 2.0*cl
         rat1 = twoz/ccl
         rat2 = rat1**2
         rat3 = csq/twoz
         il = -kap
         if (kap.gt.0) il = kap+1
         l0 = il-1
         ag(1) = agi
         if (twoz.le.0.0) then
            ap(1) = -ec/(2.0*il+1.0)*dr(1)*ag(1)
            ag(2) = 0.0
            ap(2) = 0.0
            ag(3) = 0.0
            ap(3) = 0.0
         else
            ap(1) = (fl-il)*rat3*ag(1)
            ag(2) = (3.0*fl-rat2)/(2.0*fl+1.0) * ag(1)
            ap(2)= rat3*( (fl -l0)*ag(2) - ag(1) ) -ap(1)
            ag(3)=( (fl+3.0*il)*ag(2) - 3.0*l0*ag(1) + 
     1      (fl+il+3.0)/rat3*ap(2) ) /(fl+1.0)/4.0
            ap(3)=( rat3*(2.0*l0*(fl+2.0-il)-l0-rat2)*ag(2)
     1      - 3.0*l0*rat3*(fl+2.0-il)*ag(1) + (fl+3.0-2.0*il-rat2)
     2      *ap(2) ) /(fl+1.0)/4.0
            ap(1) = ap(1)/ccl
            ag(2)= ag(2)*rat3
            ap(2)= ap(2)*rat3/ccl
            ag(3)= ag(3)*rat3**2
            ap(3)= ap(3)*rat3**2/ccl
         endif
         gg(1)=dr(1)**fl * (ag(1)+dr(1)*(ag(2)+dr(1)*ag(3)))
         gp(1)=dr(1)**fl * (ap(1)+dr(1)*(ap(2)+dr(1)*ap(3)))
      endif

      i0=1
      iflat = min ( jri, iwkb)
      call intout (en, i0, kap, iflat, ic3, vm)

      do 100 i = iflat, max0-1
         if (i.eq.iwkb) then
            eph = cl* ( 3*dv(iwkb+1) - dv(iwkb+2)) /2
            if (iwkb.eq.jri-1) eph=  cl* (dv(i) + dv(i+1)) /2
         else
            eph = cl* (dv(i) + dv(i+1)) /2
         endif
         if (ic3.gt.0 .and. i.lt.jri) then
           rav = (dr(i)+dr(i+1)) / 2
           ec = rav**3 * ( ccl+ (en - eph) / cl )**2
           eph = eph + ic3 * cl / ec * (vm(i) + vm(i+1)) / 2
         endif
         call flatv( dr(i), dr(i+1), gg(i), gp(i), en, eph, kap,
     1               gg(i+1), gp(i+1))
  100 continue

      return
      end
