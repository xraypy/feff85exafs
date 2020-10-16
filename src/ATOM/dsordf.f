      double precision function dsordf (i,j,n,jnd,a)
c              * calculation of diff. integrals*
c        integration by simpson method of the   hg*(r**n)
c        hg(l)=cg(l,i)*cg(l,j)+cp(l,i)*cp(l,j)  if jnd=1
c        hg=expression above multiplied by  dg  if jnd=-1
c        hg(l)=cg(l,i)*cp(l,j)                  if jnd=2
c        hg=expression above multiplied by  dg  if jnd=-2
c        hg(l)=dg(l)*cg(l,i)+dp(l)*cp(l,j)      if jnd=3
c        hg(l)=dg(l)*dg(l)+dp(l)*dp(l)          if jnd=4
c        hg is constructed by calling program   if jnd>=5
c                  cg(l,i)  large component of the orbital i
c                  cp(l,j)  small component of the orbital j
c        a is such that dg,dp or hg following the case
c        behave at the origin as cte*r**a
c        the integration is made as far as dr(j) for jnd>3
c
c        the development limits at the origin (used for calculation
c        of integral form 0 to dr(1) ) of functions dg,dp and hg are
c        supposed to be in blocks ag,ap and chg respectively
c        this program uses  aprdev
c
      implicit double precision (a-h,o-z)
      common cg(251,30), cp(251,30), bg(10,30), bp(10,30),
     1         fl(30), fix(30), ibgp
      common/comdir/cl,dz,dg(251),ag(10),dp(251),ap(10),bidcom(783)
      dimension hg(251),chg(10)
      common/ratom1/xnel(30),en(30),scc(30),scw(30),sce(30),
     1nq(30),kap(30),nmax(30)
      common/tabtes/hx,dr(251),test1,test2,ndor,np,nes,method,idim
      dimension bgi(10),bgj(10),bpi(10),bpj(10)

c        construction of the array hg
      b=0
      if (jnd.gt.3) then
         max0=j
         b=a
         go to 101
      endif
      max0= min(nmax(i),nmax(j))
      do l= 1,ibgp
        bgi(l) = bg(l,i)
        bgj(l) = bg(l,j)
        bpi(l) = bp(l,i)
        bpj(l) = bp(l,j)
      enddo

      if (abs(jnd).le.2) then
         if (abs(jnd).lt.2) then
            do l=1,max0
               hg(l)=cg(l,i)*cg(l,j)+cp(l,i)*cp(l,j)
            enddo
            do l=1,ndor
               chg(l)=aprdev(bgi,bgj,l)+aprdev(bpi,bpj,l)
            enddo
         else
            do l=1,max0
               hg(l)=cg(l,i)*cp(l,j)
            enddo
            do l=1,ndor
               chg(l)=aprdev(bgi,bpj,l)
            enddo
         endif
         b=fl(i)+fl(j)
         if (jnd.le.0) then
            do l=1,max0
               hg(l)=hg(l)*dg(l)
            enddo
            do l=1,ndor
               ap(l)=chg(l)
            enddo
            b=b+a
            do l=1,ndor
               chg(l)=aprdev(ap,ag,l)
            enddo
         endif
         go to 301
      endif
 101  continue

      if (jnd.eq.4) then
         do l=1,max0
            hg(l)=dg(l)*dg(l)+dp(l)*dp(l)
         enddo
         b=b+b
         do l=1,ndor
            chg(l)=aprdev(ag,ag,l)+aprdev(ap,ap,l)
         enddo
      else if (jnd.lt.4) then
         do l=1,max0
            hg(l)=dg(l)*cg(l,i)+dp(l)*cp(l,j)
         enddo
         b=a+fl(i)
         do  l=1,ndor
            chg(l)=aprdev(bgi,ag,l)+aprdev(bpj,ap,l)
         enddo
      endif
c     integration of the hg
 301  continue


      dsordf=0.0d 00
      io=n+1
      do l=1,max0
         hg(l)=hg(l)*(dr(l)**io)
      enddo
      do l=2,max0,2
         dsordf=dsordf+hg(l)+hg(l)+hg(l+1)
      enddo
      dsordf=hx*(dsordf+dsordf+hg(1)-hg(max0))/3.0d 00
c        integral from 0 to dr(1)
      b=b+n
      do l=1,ndor
         b=b+1.0d 00
         dsordf=dsordf+chg(l)*(dr(1)**b)/b
      enddo
      return
      end
