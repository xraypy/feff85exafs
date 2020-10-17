      subroutine potslw (dv,d,dr,dpas,np)
c
c coulomb potential uses a 4-point integration method
c dv=potential;  d=density;  dp=bloc de travail; dr=radial mesh
c dpas=exponential step;
c np=number of points
c **********************************************************************

      implicit double precision (a-h,o-z)
      dimension dv(251), d(251), dp(251), dr(251)

      das=dpas/24.0
      do i=1,np
         dv(i)=d(i)*dr(i)
      enddo
      dlo=exp(dpas)
      dlo2=dlo*dlo
      dp(2)=dr(1)*(d(2)-d(1)*dlo2)/(12.0*(dlo-1.0))
      dp(1)=dv(1)/3.0-dp(2)/dlo2
      dp(2)=dv(2)/3.0-dp(2)*dlo2
      j=np-1
      do i=3,j
         dp(i)=dp(i-1)+das*(13.0*(dv(i)+dv(i-1))-(dv(i-2)+dv(i+1)))
      enddo
      dp(np)=dp(j)
      dv(j)=dp(j)
      dv(np)=dp(j)
      do i=3,j
         k=np+1-i
         dv(k)=dv(k+1)/dlo+das*
     $        (13.0*(dp(k+1)/dlo+dp(k)) - (dp(k+2)/dlo2 + dp(k-1)*dlo))
      enddo
      dv(1)=dv(3)/dlo2+dpas*(dp(1)+4.0*dp(2)/dlo+dp(3)/dlo2)/3.0
      do i=1,np
         dv(i)=dv(i)/dr(i)
      enddo
      return
      end
