      subroutine potsl (dv,d,dp,dr,dpas,dexv,z,np,ion,icut,dvn)
c
c coulomb potential uses a 4-point integration method
c dv=potential;  d=density;  dp=bloc de travail; dr=radial mesh;
c dpas=exponential step; dexv=multiplicative coefficient for the exchang
c z=atomic number;  np=number of points; ion=z-number of electrons
c if icut is zero one corrects the potential by -(ion+1)/r
c **********************************************************************
      implicit double precision (a-h,o-z)
      save
      dimension dv(251), d(251), dp(251), dr(251), dvn(251)
      das=dpas/24.0
      do 10 i=1,np
   10 dv(i)=d(i)*dr(i)
      dlo=exp(dpas)
      dlo2=dlo*dlo
      dp(2)=dr(1)*(d(2)-d(1)*dlo2)/(12.0*(dlo-1.0))
      dp(1)=dv(1)/3.0-dp(2)/dlo2
      dp(2)=dv(2)/3.0-dp(2)*dlo2
      j=np-1
      do 20 i=3,j
   20 dp(i)=dp(i-1)+das*(13.0*(dv(i)+dv(i-1))-(dv(i-2)+dv(i+1)))
      dp(np)=dp(j)
      dv(j)=dp(j)
      dv(np)=dp(j)
      do 30 i=3,j
      k=np+1-i
   30 dv(k)=dv(k+1)/dlo+das*(13.0*(dp(k+1)/dlo+dp(k))-(dp(k+2)/dlo2+dp
     1 (k-1)*dlo))
      dv(1)=dv(3)/dlo2+dpas*(dp(1)+4.0*dp(2)/dlo+dp(3)/dlo2)/3.0
      dlo=-(ion+1)
      do 40 i=1,np
      dvn(i)=dv(i)/dr(i)
      dv(i)=dv(i)-(z+exchan(d(i),dr(i),dexv))
      if (icut.ne.0) go to 40
      if (dv(i).gt.dlo) dv(i)=dlo
   40 dv(i)=dv(i)/dr(i)
      return
      end
