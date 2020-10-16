      subroutine nucdev (av,dr,dv,dz,hx,nuc,np,ndor,dr1)
c        * construction of nuclear potential *
c av coefficients of the development at the origin of nuclear potential
c dr  tabulation points
c dv  nuclear potential
c dz  nuclear charge
c hx  exponential step
c nuc index of the nuclear radius
c np  number of tabulation points
c ndor number of the coefficients for development at the origin
c the declared below arguments are saved, dr1 is the first

      implicit double precision (a-h,o-z)
      dimension av(10),dr(251),dv(251),at(251)

c    specify atomic mass and thickness of nuclear shell
c a atomic mass (negative or null for the point charge)
c epai parameter of the fermi density distribution
c (negative or null for uniform distribution), which is
c       cte / (1. + exp((r-rn)/epai) )
c with nuclear radius rn= 2.2677e-05 * (a**(1/3))

c calculate radial mesh
      a = 0.0
      epai = 0.0

      if (a.le.1.0d-01) then
         nuc=1
      else
         a=dz*(a**(1./3.))*2.2677d-05
         b=a/ exp(hx*(nuc-1))
         if (b.le.dr1) then
            dr1=b
         else
            b=log(a/dr1)/hx
            nuc=3+2*int(b/2.0)
            if (nuc.ge.np) call par_stop('dr1 too small')
c           index of atomic radius larger than dimension of dr
            dr1=a*exp(-(nuc-1)*hx)
         endif
      endif

      dr(1)=dr1/dz
      do l=2,np
         dr(l)=dr(1)* exp(hx*(l-1))
      enddo
      if (ndor.lt.5) then
c       * it should be at least 5 development coefficients
         call par_stop
     .     ('stopped in programm nucdev, ndor should be > 4.')
c        stop
      endif
c  calculate nuclear potential on calculated radial mesh
      do i=1,ndor
         av(i)=0.0d00
      enddo
      if (epai.le.0.0) then
         do i=1,np
            dv(i)=-dz/dr(i)
         enddo
         if (nuc.le.1) then
            av(1)=-dz
         else
            av(2)=-3.0d 00*dz/(dr(nuc)+dr(nuc))
            av(4)=-av(2)/(3.0d 00*dr(nuc)*dr(nuc))
            l=nuc-1
            do i=1,l
               dv(i)=av(2)+av(4)*dr(i)*dr(i)
            enddo
         endif
      else
         b= exp(-dr(nuc)/epai)
         b=1.0d 00/(1.0d 00+b)
         av(4)=b
         av(5)=epai*b*(b-1.0d 00)
         if (ndor.gt.5) then
            at(1)=1.0d 00
            at(2)=1.0d 00
            nf=1
            do i=6,ndor
               n=i-4
               nf=n*nf
               dv(1)=n*at(1)
               n1=n+1
               dv(n1)=1.0d 00
               do j=2,n
                  dv(j)=(n-j+2)*at(j-1)+(n-j+1)*at(j)
               enddo
               do j=1,n1
                  m=n+1-j
                  l=1
                  if (mod(j,2).eq.0) l=-l
                  av(i)=av(i)+l*dv(j)*(b**m)
                  at(j)=dv(j)
               enddo
               av(i)=b*av(i)*(epai**n)/nf
            enddo
         endif
         do i=1,np
            b=1.0d 00+ exp((dr(i)-dr(nuc))/epai)
            if ((b*av(4)).gt.1.0d+15) go to 51
            dv(i)=dr(i)*dr(i)*dr(i)/b
            l=i
         enddo
 51      continue
         if (l.ge.(np-1)) l=np-2
         k=l+1
         do i=k,np
            dv(i)=0.0d 00
         enddo
         at(1)=0.0d 00
         at(2)=0.0d 00
         k=2
         do i=4,ndor
            k=k+1
            do j=1,2
               at(j)=at(j)+av(i)*(dr(j)**k)/k
            enddo
            av(i)=av(i)/(k*(k-1))
            av(2)=av(2)+av(i)*(dr(1)**k)
         enddo
         a=hx/2.4d+01
         b=a*1.3d+01
         k=l+1
         do i=3,k
            at(i)=at(i-1)+b*(dv(i-1)+dv(i))-a*(dv(i-2)+dv(i+1))
         enddo
         dv(l)=at(l)
         do i=k,np
            dv(i)=dv(l)
         enddo
         e= exp(hx)
         c=1.0d 00/(e*e)
         i=l-1
 83      continue
         dv(i)=dv(i+1)/e+b*(at(i+1)/e+at(i))-a*(at(i+2)*c+at(i-1)*e)
         i=i-1
         if (i .gt. 1) go to 83
         dv(1)=dv(3)*c+hx*(at(1)+4.0d 00*at(2)/e+at(3)*c)/3.0d 00
         av(2)=(av(2)+dv(1))/dr(1)
         a=-dz/dv(l)
         do i=4,ndor
            av(i)=-a*av(i)
         enddo
         av(2)=a*av(2)
         do i=1,np
            dv(i)=a*dv(i)/dr(i)
         enddo
      endif

      return
      end
