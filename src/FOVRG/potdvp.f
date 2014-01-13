      subroutine potdvp
c     this programm uses aprdep,multrk,yzkrdf
c     to calculate potential development coefficients
      implicit double precision (a-h,o-z)
      include '../HEADERS/dim.h'
      common/dff/ cg(nrptx,30), cp(nrptx,30), bg(10,30), bp(10,30),
     1              fl(30), fix(30), ibgp
      complex*16 dg,ag,dp,ap,dv,av,eg,ceg,ep,cep
      common/comdic/cl,dz,dg(nrptx),ag(10),dp(nrptx),ap(10),dv(nrptx),
     2         av(10),eg(nrptx),ceg(10),ep(nrptx),cep(10)
c     dg,dp to get data from yzkrdf, dv,eg,ep -output for soldir
      common/itescf/testy,rap(2),teste,nz,norb,norbsc
      common/ratom1/xnel(30),en(30),scc(30),scw(30),sce(30),
     1nq(30),kap(30),nmax(30)
      common/snoyac/dvn(nrptx),anoy(10),nuc
      common/tabtec/hx,dr(nrptx),test1,test2,ndor,np,nes,method,idim
      dimension bgj(10),bpj(10)
c#mn
       external aprdep

      do 9 i=1,10
 9       av(i)=anoy(i)
 
c     calculate density development coefficients
      do 31 i=1,ndor
 31   ag(i)=0.0d 00
      do 51 j=1,norb-1
         do 33 i = 1,10
            bgj(i) = bg(i,j)
 33         bpj(i) = bp(i,j)
         n=2* abs(kap(j))
         l=ndor+2-n
         if (l.le.0) go to 51
         do 41 i=1,l
            m=n-2+i
 41         ag(m)=ag(m)+xnel(j)*(aprdep(bgj,bgj,i)+
     1            aprdep(bpj,bpj,i))*fix(j)**2
 51   continue

c     transform density coefficients into ones for potential
      ap(1)=0.0d 00 
      do 15 i=1,ndor
         ag(i)=ag(i)/(i+2)/(i+1)
         ap(1)=ap(1)+ag(i)*dr(1)**(i+1)
 15   continue

      do 61 i=1,ndor
         l=i+3
         if (l.gt.ndor) go to 61
         av(l)=av(l)-ag(i)
 61   continue
c     av(2)=avoy(2) + ap(1)+(vxcvzl(1)-dvn(1)) in order 
c     to have sum av(i)*dr(1)**(i-2)=vxcval(1)
      av(2)=av(2)+ap(1)
 
c addition of nuclear potential and division of potentials and
c       their development limits by speed of light
      do 527 i=1,10
 527     av(i)=av(i)/cl
      return
      end
