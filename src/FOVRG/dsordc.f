      complex*16 function dsordc(j,a,dg,dp,ag,ap)
c              * calculation of overlap integrals*
c        integration by simpson method of the   hg*(r**0)
c        hg(l)=dg(l)*cg(l,j)+dp(l)*cp(l,j)
c                cg,cp(l,j)  orbital j
c        a is such that dg,dp or hg following the case
c        behave at the origin as cte*r**a
c        the development limits at the origin (used for calculation
c        of integral form 0 to dr(1) ) of functions dg,dp and hg are
c        supposed to be in blocks ag,ap and chg respectively
c        this program uses   aprdec
c
      implicit double precision (a-h,o-z)
      include '../HEADERS/dim.h'
      complex*16 aprdec
      common/dff/ cg(nrptx,30), cp(nrptx,30), bg(10,30), bp(10,30),
     1              fl(30), fix(30), ibgp
      complex*16 dg(nrptx),ag(10),dp(nrptx),ap(10)
      complex*16 hg(nrptx),chg(10)
c     common/ratom1/xnel(30),en(30),scc(30),scw(30),sce(30),
c    1   nq(30),kap(30),nmax(30)
      common/tabtec/hx,dr(nrptx),test1,test2,ndor,np,nes,method,idim
      dimension bgj(10),bpj(10)

c        construction of the array hg
      do  15 l= 1,ibgp
        bgj(l) = bg(l,j)
 15     bpj(l) = bp(l,j)

      do 221 l=1,idim
 221  hg(l)=dg(l)*cg(l,j)+dp(l)*cp(l,j)
      b=a+fl(j)
      do 241 l=1,ndor
 241     chg(l) = aprdec(ag,bgj,l) + aprdec(ap,bpj,l)
 
c        integration of the hg
      dsordc = (0.0d0, 0.0d0)
      do 305 l=1,idim
 305     hg(l)=hg(l)*dr(l)
      do 311 l=2,idim,2
 311     dsordc=dsordc+hg(l)+hg(l)+hg(l+1)
      dsordc=hx*(dsordc+dsordc+hg(1)-hg(idim))/3.0d0
c        integral from 0 to dr(1)
      do 331 l=1,ndor
         b=b+1.0d 00
 331     dsordc=dsordc+chg(l)*(dr(1)**b)/b
      return
      end
