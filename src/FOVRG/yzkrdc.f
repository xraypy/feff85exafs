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
      do  5 l= 1,ibgp
        bgi(l) = bg(l,i)
  5     bpi(l) = bp(l,i)
      id=min(nmax(i),np)
      ap(1)=fl(i)+flps
      do 11 l=1,id
 11   dg(l)=cg(l,i)*ps(l)+cp(l,i)*qs(l)
      do 12 l = id+1,idim
 12    dg(l) = 0.0d0
      do 21 l=1,ndor
 21   ag(l) = aprdec(aps,bgi,l) + aprdec(aqs,bpi,l)

      dyzk = 0
c     if (id .ge. nmax(norb)) then
c        id = nmax(norb)-1
c        ck0 = log(cg(id,i)/cg(id+1,i))  / (dr(id+1)-dr(id))
c        ck = sqrt(2*p2)
c        xck = ck/cl
c        xck = -xck/(1+sqrt(1+xck**2))
c        temp = -ps(id+1) / qs(id+1) *xck
c        xx = dble (temp)
c        yy = dimag(temp)
c        if (xx .ne. 0)  then
c            alph = (1 - xx**2 - yy**2)
c            alph = sqrt(alph**2 + 4*xx**2) - alph
c            alph = alph / (2 * xx)
c            alph = atan (alph)
c        else
c            alph = 0
c        endif
c        beta = (xx**2 + (yy+1)**2) / (xx**2 + (yy-1)**2)
c        beta = log(beta) / 4

c        phx = dcmplx (alph, beta)
c        a1 =   ps(id+1) / sin(phx)
c        a2 = - qs(id+1) / cos(phx)
c        xck=ck*dr(id+1)
c        phx = phx -xck
c        a1 = a1*cg(id+1,i)/2/coni
c        a2 = a2*cp(id+1,i)/2
c        b1=exp(coni*phx) * (a1 - a2)
c        b2=exp(-coni*phx) * (-a1 - a2)
c        xck = (ck0 - coni*ck)*dr(id+1)
c        n = k +1
c        dyzk = dyzk + b1*exp(-xck)/xck
c        dyzk = dyzk + b1*expint(n,xck)
c        xck = (ck0 + coni*ck)*dr(id+1)
c        dyzk = dyzk + b2*exp(-xck)/xck
c        dyzk = dyzk + b2*expint(n,xck)
c        dyzk = dyzk*dr(id+1)
c     endif

      call yzktec (dg,ag,dp,chg,dr,ap(1),hx,k,ndor,id,idim, dyzk)
      return
      end

c     complex*16 function expint(n,x)
c     implicit double precision (a-h,o-z)
c     integer n, maxit
c     complex*16 x, b, c, d, h, del, fact, zero
c     parameter (zero=(0.d0,0.d0))
c     parameter (maxit=100, eps=1.d-7, fpmin=1.d-30, euler=.5772156649)

c     nm1 = n - 1
c     if (n.lt.0 .or. (dble(x).lt.0.d0 .and. dimag(x).eq.0.d0) .or.
c    1     (x.eq.zero .and. (n.eq.0.or.n.eq.1))) then
c        call par_stop('Bad arguments in expint')
c     elseif (n.eq.0) then
c        expint = exp(-x) / x
c     elseif (x.eq.0) then
c        expint = 1.d0 /nm1
c     elseif (dble(x).gt.1) then
c        b = x + n
c        c = 1/fpmin
c        d = 1/b
c        h = d
c        do 10 i=1,maxit
c           a = -i*(nm1+i)
c           b = b + 2
c           d = 1 / (a*d+b)
c           c = b + a/c
c           del = c*d
c           h = h*del
c           if (abs(del-1) .lt. eps) then
c              expint = h * exp(-x)
c              return
c           endif
c 10     continue
c        call par_stop(' continued fraction failed in expint')
c     else
c        if (nm1.ne.0) then
c           expint = 1/nm1
c        else
c           expint = -log(x) - euler
c        endif
c        fact = 1
c        do 30 i=1,maxit
c           fact = - fact *x / i
c           if (i.ne.nm1) then
c              del = - fact / (i-nm1)
c           else
c              psi = - euler
c              do 20 ii=1,nm1
c                 psi = psi + 1.d0 / ii
c 20           continue
c              del = fact*(-log(x)+psi)
c           endif
c           expint = expint + del
c           if (abs(del).lt.abs(expint)*eps) return
c 30     continue
c        call par_stop('series failed in expint')
c     endif
c     return
c     end
