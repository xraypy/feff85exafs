c      subroutine vlda(ia, xnval,srho, srhovl,vtrho, ilast, idfock)
      subroutine vlda(xnval,srho, srhovl,vtrho, ilast, idfock)
c    this program calculates xc potential, using core-vlaence separation
c    discussed in ankuodinov's thesis.  
c    written by alexei ankoudinov. 11.07.96
      implicit double precision (a-h,o-z)
      include '../HEADERS/const.h'
      dimension xnval(30), srho (251), srhovl(251), vtrho(251)
      common cg(251,30), cp(251,30), bg(10,30), bp(10,30),
     1        fl(30), fix(30), ibgp
      common/comdir/cl,dz,dg(251),ag(10),dp(251),ap(10),dv(251),av(10),
     2              eg(251),ceg(10),ep(251),cep(10)
      common/itescf/testy,rap(2),teste,nz,norb,norbsc
      common/ratom1/xnel(30),en(30),scc(30),scw(30),sce(30),
     1nq(30),kap(30),nmax(30)
      common/tabtes/hx,dr(251),test1,test2,ndor,np,nes,method,idim
 
      do 10 i = 1, 251
        srho(i)   = zero
        srhovl(i) = zero
 10   continue 
c  find total and valence densities. Remove self-interaction if SIC
      do 50 j = 1, norb
         a = xnel(j)
         b = xnval(j)
c     use to test SIC
c       if (j .eq. ia) a=a-1.0d0
c       if (j .eq. ia) b=b-1.0d0
         do 50 i = 1, nmax(j)
            srho(i)   = srho(i)   + a * (cg(i,j)**2+cp(i,j)**2)
            srhovl(i) = srhovl(i) + b * (cg(i,j)**2+cp(i,j)**2)
 50   continue 

c  constract lda potential. Put your favorite model into vbh.f.
c  exch=5,6 correspond to 2 ways of core-valence separation of V_xc.
      rhoc = zero
      do 90 i = 1, 251
        rho = srho(i) / (dr(i)**2)
        if (idfock.eq.5) then
c          for exch=5 valence density*4*pi
           rhoc = srhovl(i) / (dr(i)**2)
        elseif (idfock.eq.6) then
c          for exch=6 core density*4*pi
           rhoc = (srho(i)-srhovl(i)) / (dr(i)**2)
        elseif (idfock.eq.1) then
           rhoc = zero
        elseif (idfock.eq.2) then
           rhoc = srho(i) / (dr(i)**2)
        else
            call par_stop(' undefined idfock in subroutine vlda')
        endif

        if (rho .gt. zero ) then
           rs = (rho/3)**(-third)
           rsc = 101.d0
           if (rhoc .gt.zero) rsc = (rhoc/3)**(-third)
           xm = one
c          vbh and edp in Hartrees
           if (idfock.eq.5 .or. idfock.eq.2) then
c             for exch=5, 2
              call vbh(rsc, xm, vxcvl)
           elseif (idfock.eq.6) then
c             for exch=6
              call vbh(rs, xm, vvbh)
                 xf = fa/rs
              call edp(rsc,xf,vdh)
              vxcvl = vvbh - vdh
           elseif (idfock.eq.1) then
c          for pure Dirac-Fock
              vxcvl = zero
           endif

c   contribution to the total energy from V_xc:=\int d^3 r V_xc * rho/2
           if (ilast.gt.0) vtrho (i) = vtrho(i) +
     1         vxcvl * srho(i)
c    1         vxcvl * xnel(ia)*(cg(i,ia)**2+cp(i,ia)**2)
c           use to test SIC

c  add to the total potential and correct it's development coefficients
           if (i.eq.1) av(2) = av(2) +vxcvl/cl
           dv(i) = dv(i) +vxcvl/cl
        endif
 90   continue

      return
      end
