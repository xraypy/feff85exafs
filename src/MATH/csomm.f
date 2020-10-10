      subroutine csomm (dr,dp,dq,dpas,da,m,np)
c Modified to use complex p and q.  SIZ 4/91
c integration by the method of simpson of (dp+dq)*dr**m from
c 0 to r=dr(np)
c dpas=exponential step;
c for r in the neighborhood of zero (dp+dq)=cte*r**da
c **********************************************************************
      implicit double precision (a-h,o-z)
      dimension dr(*)
      double precision d1, db, dd
      complex*16 dp(*),dq(*),da,dc
      mm=m+1
      d1=dble(da+mm)
      da=0.d0
      db=0.d0
      do i=1,np
         dl=dr(i)**mm
         if (i.eq.1.or.i.eq.np) go to 10
         dl=dl+dl
         if ((i-2*(i/2)).eq.0) dl=dl+dl
 10      dc=dp(i)*dl
         da=da+dc
         dc=dq(i)*dl
         da=da+dc

      enddo
      da=dpas*da/3
      dd=exp(dpas)-1.d0
      db=d1*(d1+1.d0)*dd*exp((d1-1.d0)*dpas)
      db=dr(1)*(dr(2)**m)/db
      dd=(dr(1)**mm)*(1.d0+1.d0/(dd*(d1+1.d0)))/d1
      da=da+dd*(dp(1)+dq(1))-db*(dp(2)+dq(2))
      return
      end
