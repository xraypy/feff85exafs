      subroutine potrdf (ia)
c        this programm uses akeato(bkeato),aprdev,multrk,yzkrdf
      implicit double precision (a-h,o-z)
      common cg(251,30), cp(251,30), bg(10,30), bp(10,30),
     1        fl(30), fix(30), ibgp
      common/comdir/cl,dz,dg(251),ag(10),dp(251),ap(10),dv(251),av(10),
     2              eg(251),ceg(10),ep(251),cep(10)
c     dg,dp to get data from yzkrdf, dv,eg,ep -output for soldir
      dimension at(251),bt(251)
      common/itescf/testy,rap(2),teste,nz,norb,norbsc
      common/ratom1/xnel(30),en(30),scc(30),scw(30),sce(30),
     1nq(30),kap(30),nmax(30)
      common/scrhf1/eps(435),nre(30),ipl
      common/snoyau/dvn(251),anoy(10),nuc
      common/tabtes/hx,dr(251),test1,test2,ndor,np,nes,method,idim
      dimension bgj(10),bpj(10)
c#mn
       external akeato, bkeato, aprdev

      do i=1,ndor
         cep(i)=0.0d 00
         ceg(i)=0.0d 00
         av(i)=anoy(i)
      enddo
      do i=1,idim
         at(i)=0.0d 00
         bt(i)=0.0d 00
         ep(i)=0.0d 00
         eg(i)=0.0d 00
         dv(i)=0.0d 00
      enddo

c     coulomb terms
      jia=2* abs(kap(ia))-1
      k=0
 21   continue
      do i=1,idim
         dg(i)=0.0d 00
      enddo
      do i=1,ndor
         ag(i)=0.0d 00
      enddo
      max0=0
      do j=1,norb
         do i = 1,10
            bgj(i) = bg(i,j)
            bpj(i) = bp(i,j)
         enddo
         m=2* abs(kap(j))-1
         if (k.le.m) then
            a=akeato(ia,j,k)/xnel(ia)
            if (a.ne.0.0d 00) then
               m=nmax(j)
               do i=1,m
                  dg(i)=dg(i)+a*(cg(i,j)*cg(i,j)+cp(i,j)*cp(i,j))
               enddo
               n=2* abs(kap(j))-k
               l=ndor+2-n
               if (l.gt.0) then
c     quick fix of development coefficients
                  a = a * fix(j)**2
                  do i=1,l
                     m=n-2+i
                     ag(m)=ag(m)+a*(aprdev(bgj,bgj,i)+
     1                    aprdev(bpj,bpj,i))
                  enddo
               endif
            endif
         endif

         max0= max(max0,nmax(j))
      enddo
      call yzkrdf (0,max0,k)
      do i=1,ndor
         l=k+i+3
         if (l.le.ndor) then
            av(l)=av(l)-ag(i)
         endif
      enddo
      do i=1,idim
         dv(i)=dv(i)+dg(i)
      enddo
      k=k+2
      if (k.le.ndor) av(k)=av(k)+ap(1)
      if (k.lt.jia) go to 21

c     exchange terms
      if (method.ne.0) then
         do j=1,norb
            if (j.ne.ia) then
               max0=nmax(j)
               jj=2* abs(kap(j))-1
               kma=(jj+jia)/2
               k= abs(jj-kma)
               if ((kap(j)*kap(ia)).lt.0) k=k+1

 111           continue
               a=bkeato(j,ia,k)/xnel(ia)
               if (a.ne.0.0d 00) then
                  call yzkrdf (j,ia,k)
                  do i=1,max0
                     eg(i)=eg(i)+a*dg(i)*cg(i,j)
                     ep(i)=ep(i)+a*dg(i)*cp(i,j)
                  enddo
                  n=k+1+ abs(kap(j))- abs(kap(ia))
                  if (n.le.ndor) then
                     do i=n,ndor
                        ceg(i)=ceg(i)+bg(i+1-n,j)*a*ap(1)*fix(j)/fix(ia)
                        cep(i)=cep(i)+bp(i+1-n,j)*a*ap(1)*fix(j)/fix(ia)
                     enddo
                  endif
                  i=2* abs(kap(j))+1
                  if (i.le.ndor) then
                     do ix = 1,10
                        bgj(ix) = bg(ix,j)
                        bpj(ix) = bp(ix,j)
                     enddo
                     do n=i,ndor
                        ceg(n)=ceg(n)-a*aprdev(ag,bgj,n+1-i)*fix(j)**2
                        cep(n)=cep(n)-a*aprdev(ag,bpj,n+1-i)*fix(j)**2
                     enddo
                  endif
               endif
               k=k+2
               if (k.le.kma) go to 111
            endif
         enddo
      endif

      if (ipl.ne.0) then
         do j=1,norbsc
            if (kap(j).ne.kap(ia).or.j.eq.ia) go to 481
            if (nre(j).lt.0.and.nre(ia).lt.0) go to 481
            m= max(j,ia)
            i= min(j,ia)+((m-1)*(m-2))/2
            a=eps(i)*xnel(j)
            max0=nmax(j)
            do i=1,max0
               at(i)=at(i)+a*cg(i,j)
               bt(i)=bt(i)+a*cp(i,j)
            enddo
            do i=1,ndor
               ceg(i)=ceg(i)+bg(i,j)*a
               cep(i)=cep(i)+bp(i,j)*a
            enddo
 481        continue
         enddo
      endif

c addition of nuclear potential and division of potentials and
c       their development limits by speed of light
      do i=1,ndor
         av(i)=av(i)/cl
         cep(i)=cep(i)/cl
         ceg(i)=ceg(i)/cl
      enddo
      do i=1,idim
         dv(i)=(dv(i)/dr(i)+dvn(i))/cl
         ep(i)=(ep(i)+bt(i)*dr(i))/cl
         eg(i)=(eg(i)+at(i)*dr(i))/cl
      enddo
      return
      end
