c      subroutine yzkrdc (i,k,flps,ps,qs,aps,aqs,p2, norb)
      subroutine yzkrdc (i,k,flps,ps,qs,aps,aqs)
c       * calculate  function yk *
c yk = r * integral of f(s)*uk(r,s)
c uk(r,s) = rinf**k/rsup**(k+1)   rinf=min(r,s)   rsup=max(r,s)
c j=norb for photoelectron
c f(s)=cg(s,i)*cg(s,j)+cp(s,i)*cp(s,j)
c f(s) is constructed by the calling programm  if i < or =0
c in the last case a function f (lies in the block dg) is supposedly
c tabulated untill point dr(j), and its' devlopment coefficients
c at the origin are in ag and the power in r of the first term is k+2

c the output functions yk and zk are in the blocks dp and dg.
c at the origin  yk = cte * r**(k+1) - developement limit,
c cte lies in ap(1) and development coefficients in ag.
c        this programm uses aprdec and yzktec
 
      implicit double precision (a-h,o-z)
      include '../HEADERS/dim.h'
      complex*16 aprdec,dyzk
c     complex*16 p2
c     complex*16 a1,a2,b1,b2,coni
c     complex*16 xck, temp, ck, phx
c      parameter (coni=(0.d0,1.d0))
      complex*16 ps(nrptx),qs(nrptx),aps(10),aqs(10)
      common/dff/cg(nrptx,30), cp(nrptx,30), bg(10,30), bp(10,30),
     1             fl(30), fix(30), ibgp
      complex*16 dg,ag,dp,ap,bidcom, chg(10)
      common/comdic/cl,dz,dg(nrptx),ag(10),dp(nrptx),ap(10),
     1   bidcom(3*nrptx+30)
      common/ratom1/xnel(30),en(30),scc(30),scw(30),sce(30),
     1   nq(30),kap(30),nmax(30)
      common/tabtec/hx,dr(nrptx),test1,test2,ndor,np,nes,method,idim
      dimension bgi(10),bpi(10)
c#mn
       external aprdec
 
c     construction of the function f
      do l= 1,ibgp
        bgi(l) = bg(l,i)
        bpi(l) = bp(l,i)
      enddo
      id=min(nmax(i),np)
      ap(1)=fl(i)+flps
      do l=1,id
         dg(l)=cg(l,i)*ps(l)+cp(l,i)*qs(l)
      enddo
      do l = id+1,idim
         dg(l) = 0.0d0
      enddo
      do l=1,ndor
         ag(l) = aprdec(aps,bgi,l) + aprdec(aqs,bpi,l)
      enddo

      dyzk = 0

      call yzktec (dg,ag,dp,chg,dr,ap(1),hx,k,ndor,id,idim, dyzk)
      return
      end

