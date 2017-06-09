c---------------------------------------------------------------------
c     program sigms.f
c
c     calculates debye-waller factors for each multiple
c     scattering path using Debye-Model correlations
c
c     files:  input  pathd_all.dat  multiple scattering path data
c             output fort.3  sig**2 vs path
c                    fort.2  long output
c
c     version 1  (29 july 91)
c
c     coded by j. rehr
c     path data from s. zabinsky
c
c     modified to use pdata.inp, Dec 1991, siz
c     Subroutine version, Dec 1991, siz
c
c---------------------------------------------------------------------

      subroutine sigms (tk, thetad, rs, nlegx, nleg, rat, iz, sig2)
c               tk temperature in degrees K
c               thetad debye temp in degrees K
c               rs=wigner seitz or norman radius in bohr, averaged
c                  over entire problem
c                  (4pi/3)*rs**3 = sum( (4pi/3)rnrm**3 ) / N
c                  (sum is over all atoms in the problem)
c               nlegx used in dimensions of rat and iz
c               nleg nlegs in path
c               rat positions of each atom in path (in bohr)
c               iz atomic number of each atom in path
c               NB Units of distance in this routine
c                  are angstroms, including sig**2
c               sig2 is output, debye waller factor in bohr**-2

      implicit double precision (a-h,o-z)

      include 'const.h'

c     nlegx is max number of atoms in any one path
      dimension rat(3,0:nlegx)
      dimension iz(0:nlegx)

c      parameters
c               x = k_d*R   (distance parameter)
c               R distance in angstroms
c               y = hbar omegad/kT = thetad/t
c               thetad debye temp in degrees K
c               tk temperature in degrees K
c               k_d = (6*pi**2 N/V) = debye wave number
c               N/V=1/(4pi/3rs**3)
c               rs=wigner seitz or norman radius in bohr
c               ami, amj masses at sites i and j in amu
c               I = int_0^1 (y/x) dw sin(wx)coth(wy/2)

c     Note:  There are nleg atoms including the central atom
c            index 0 and index nleg both refer to central atom,
c            which makes special code unnecessary later.
      sum = 0
      ntot = 0

      sigtot=0
      do 800 il=1,nleg
      do 800 jl=il,nleg

c        calculate r_i-r_i-1 and r_j-r_j-1

         rij = dist (rat(1,il), rat(1,jl))
         call corrfn (rij, cij, thetad, tk, iz(il), iz(jl), rs)
         sig2ij=cij

         rimjm = dist (rat(1,il-1), rat(1,jl-1))
         call corrfn (rimjm, cimjm, thetad, tk, iz(il-1), iz(jl-1), rs)
         sig2ij=sig2ij+cimjm

         rijm = dist (rat(1,il), rat(1,jl-1))
         call corrfn (rijm, cijm, thetad, tk, iz(il), iz(jl-1), rs)
         sig2ij=sig2ij-cijm

         rimj = dist (rat(1,il-1), rat(1,jl))
         call corrfn (rimj, cimj, thetad, tk, iz(il-1), iz(jl), rs)
         sig2ij=sig2ij-cimj

         riim = dist (rat(1,il), rat(1,il-1))
         rjjm = dist (rat(1,jl), rat(1,jl-1))

         ridotj=(rat(1,il)-rat(1,il-1))*(rat(1,jl)-rat(1,jl-1))+
     1          (rat(2,il)-rat(2,il-1))*(rat(2,jl)-rat(2,jl-1))+
     2          (rat(3,il)-rat(3,il-1))*(rat(3,jl)-rat(3,jl-1))
         ridotj=ridotj/(riim*rjjm)

c        double count i .ne. j  terms
         if(jl.ne.il) sig2ij=2*sig2ij
         sig2ij=sig2ij*ridotj
         sigtot=sigtot+sig2ij

  800 continue
      sig2=sigtot/4

c     sig2 is in bohr**2, just as we wanted for ff2chi
      return
      end



      subroutine corrfn(rij,cij,thetad,tk,iz1,iz2,rsavg)
c     subroutine calculates correlation function
c     c(ri,rj)=<xi xj> in the Debye approximation
c
c             =(1/N)sum_k exp(ik.(Ri-Rj))(1/sqrt(mi*mj))*
c              (hbar/2w_k)*coth(beta hbar w_k/2)
c             = (3kT/mu w_d**2)*sqrt(mu**2/mi*mj)*I
c
c      parameters
c               x = k_d*R   (distance parameter)
c               R distance in angstroms
c               y = hbar omegad/kT = thetad/t
c               thetad debye temp in degrees K
c               tk temperature in degrees K
c               k_d = (6*pi**2 N/V) = debye wave number
c               N/V=1/(4pi/3rs**3)
c               rs=wigner seitz or norman radius in bohr
c               ami, amj masses at sites i and j in amu
c               I = int_0^1 (y/x) dw sin(wx)coth(wy/2)
c
c      solution by numerical integration
c
      implicit double precision (a-h, o-z)
      common /xy/ x, yinv

      include 'const.h'
       external atwtd, dist
c     con=hbar**2/kB*amu)*10**20   in ang**2 units
c     hbar = 1.054 572 666 e-34, amu = 1.660 540 e-27, 
c     kB = 1.380 6581 d-23
      parameter (con = 48.508 459 393 094)

c     external fn
c     rij=2.55
c     tk=295
c     thetad=315
c     ami=amj=63.55 at wt for Cu
c     rs=2.7

      ami=atwtd(iz1)
      amj=atwtd(iz2)
      rs=rsavg
c     thetad in degrees K, t temperature in degrees K
c     y=thetad/tk
      yinv=tk/thetad
      xkd=(9*pi/2)**(third)/(rs*bohr)
      fac=(3/2.)*con/(thetad*sqrt(ami*amj))
      rj=rij
      x=xkd*rj
c     call numerical integration
      call bingrt (grater, eps, nx)
      cij=fac*grater
      return
      end
      double precision function fn(w)
      implicit double precision (a-h,o-z)
      common/xy/x,yinv
c     fn=(sin(wx)/x)*coth(wy/2)
c     change code to allow t=0 without bombing
c     fn=2/y
      fn=2*yinv
      if(w.lt.1.e-20) return
      fac=w
      if(x.gt.0.) fac=sin(w*x)/x
      emwy=0.
      if(yinv.gt.0.0125) emwy=exp(-w/yinv)
      emwy=exp(-w/yinv)
      fn=fac*(1+emwy)/(1-emwy)
      return
      end
c-----------------------------------------------
      subroutine bingrt (b, eps, n)
c     subroutine calculates integrals between [0,1]
c      b = int_0^1 f(z) dz
c     by trapezoidal rule and binary refinement
c     (romberg integration)
c     coded by j rehr (10 Feb 92)
c     see, e.g., numerical recipes for discussion
c     and a much fancier version
c-----------------------------------------------
c     del=dz  itn=2**n tol=1.e-5
c     starting values
      implicit double precision (a-h,o-z)
       character messag*128
      common /xy/x,yinv
c     external fn
c     error is approximately 2**(-2n) ~ 10**(-.6n)
c     so nmax=10 implies an error of 1.e-6
      parameter(nmax = 10, tol = 1.e-5)
      parameter(zero=0, one=1)
      n=0
      itn=1
      del=1.
      bn=(fn(zero)+fn(one))/2
      bo=bn
 10   continue
c     nth iteration
c     b_n+1=(b_n)/2+deln*sum_0^2**n f([2n-1]deln)
      n=n+1
      if(n.gt.nmax) go to 40
      del=del/2
      sum=0.
      do 20 i=1, itn
      zi=(2*i-1)*del
 20   sum=sum+fn(zi)
c     bnp1=b_n+1 is current value of integral
      bnp1=bn/2+del*sum
c     cancel leading error terms b=[4b-bn]/3
c     note: this is the first term in the
c     neville table - remaining errors were
c     found too small to justify the added code
      b=(4*bnp1-bn)/3
      eps=abs((b-bo)/b)
      if(eps.lt.tol) goto 60
      bn=bnp1
      bo=b
      itn=itn*2
      goto 10
 40   continue 
       write(messag, 50) n,itn, b,eps
       call echo(messag)
 50   format(' not converged, n,itn,b,eps=',
     1  2i4,2e14.6)
      return
 60   continue
c     print70, n, itn, b, eps
c70   format(' n,itn,b,eps=' 2i4,2e16.8)
      return
      end
