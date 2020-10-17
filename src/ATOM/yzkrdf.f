      subroutine yzkrdf (i,j,k)
c       * calculate  function yk *
c yk = r * integral of f(s)*uk(r,s)
c uk(r,s) = rinf**k/rsup**(k+1)   rinf=min(r,s)   rsup=max(r,s)
c f(s)=cg(s,i)*cg(s,j)+cp(s,i)*cp(s,j)      if nem=0
c f(s)=cg(s,i)*cp(s,j)                      if nem is non zero
c f(s) is constructed by the calling programm  if i < or =0
c in the last case a function f (lies in the block dg) is supposedly
c tabulated untill point dr(j), and its' devlopment coefficients
c at the origin are in ag and the power in r of the first term is k+2

c the output functions yk and zk are in the blocks dp and dg.
c at the origin  yk = cte * r**(k+1) - developement limit,
c cte lies in ap(1) and development coefficients in ag.
c        this programm uses aprdev and yzkteg
 
      implicit double precision (a-h,o-z)
      common cg(251,30), cp(251,30), bg(10,30), bp(10,30),
     1         fl(30), fix(30), ibgp
      common/comdir/cl,dz,dg(251),ag(10),dp(251),ap(10),bidcom(783)
      dimension chg(10)
      common/ratom1/xnel(30),en(30),scc(30),scw(30),sce(30),
     1nq(30),kap(30),nmax(30)
      common/tabtes/hx,dr(251),test1,test2,ndor,np,nes,method,idim
      common/inelma/nem
      dimension bgi(10),bgj(10),bpi(10),bpj(10)
c#mn
       external aprdev
 
       if (i.le.0) then
          ap(1)=k+2
          id=j
       else
c     construction of the function f
          do l= 1,ibgp
             bgi(l) = bg(l,i)
             bgj(l) = bg(l,j)
             bpi(l) = bp(l,i)
             bpj(l) = bp(l,j)
          enddo
          id= min(nmax(i),nmax(j))
          ap(1)=fl(i)+fl(j)
          if (nem.eq.0) then
             do l=1,id
                dg(l)=cg(l,i)*cg(l,j)+cp(l,i)*cp(l,j)
             enddo
             do l=1,ndor
                ag(l)=aprdev(bgi,bgj,l)+aprdev(bpi,bpj,l)
             enddo
          else
             do l=1,id
                dg(l)=cg(l,i)*cp(l,j)
             enddo
             do l=1,ndor
                ag(l)=aprdev(bgi,bpj,l)
             enddo
          endif
       endif
      call yzkteg (dg,ag,dp,chg,dr,ap(1),hx,k,ndor,id,idim)
      return
      end
