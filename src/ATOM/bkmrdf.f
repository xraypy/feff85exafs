      subroutine bkmrdf (i,j,k)
c     angular coefficients for the breit term. i and j are the numbers
c     of orbitals and  k is the value of k in uk(1,2)
c        this programm uses cwig3j
c     coefficients for magnetic interaction  are in cmag
c     and those for retarded term are in cret
c     the order correspond to -1 0 and +1

      implicit double precision (a-h,o-z)
      common/ratom1/xnel(30),en(30),scc(30),scw(30),sce(30),
     1nq(30),kap(30),nmax(30)
      common/tabre/cmag(3),cret(3)
c#mn
       external cwig3j

      do l=1,3
        cmag(l)=0.0d 00
        cret(l)=0.0d 00
      enddo
      ji=2* abs(kap(i))-1
      jj=2* abs(kap(j))-1
      kam=kap(j)-kap(i)
      l=k-1
      do m=1,3
         if (l.ge.0) then
            a=cwig3j(ji,jj,l+l,-1,1,2)**2
            if (a.ne.0.0d00) then
               c=l+l+1
               if (m.eq.2) then
                  d=k*(k+1)
                  cm=(kap(i)+kap(j))**2
                  cz=cm
                  cp=cm
               else
                  if (m.lt.2) then
                     cm=(kam+k)**2
                     cz=kam*kam-k*k
                     cp=(k-kam)**2
                     n=k
                  else if (m.gt.2) then
                     cm=(kam-l)**2
                     cz=kam*kam-l*l
                     cp=(kam+l)**2
                     n=l
                     c=-c
                  endif
                  l1=l+1
                  am=(kam-l)*(kam+l1)/c
                  az=(kam*kam+l*l1)/c
                  ap=(l+kam)*(kam-l1)/c
                  d=n*(k+k+1)
                  c= abs(c)*d
                  if (c.ne.0.0d 00) c=n/c
                  cret(1)=cret(1)+a*(am-c*cm)
                  cret(2)=cret(2)+(a+a)*(az-c*cz)
                  cret(3)=cret(3)+a*(ap-c*cp)
               endif
            endif
            if (d.ne.0.0d 00) then
               a=a/d
               cmag(1)=cmag(1)+cm*a
               cmag(2)=cmag(2)+cz*(a+a)
               cmag(3)=cmag(3)+cp*a
            endif
         endif
      l=l+1
      enddo
      return
      end
