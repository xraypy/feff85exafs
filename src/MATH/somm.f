      subroutine somm (dr,dp,dq,dpas,da,m,np)
c
c integration by the method of simpson of (dp+dq)*dr**m from
c 0 to r=dr(np)
c dpas=exponential step;
c for r in the neighborhood of zero (dp+dq)=cte*r**da
c **********************************************************************
      implicit double precision (a-h,o-z)
      dimension dr(np), dp(np), dq(np)
      mm=m+1
      d1=da+mm
      da=0.d0
      db=0.d0
      do i=1,np
         dl=dr(i)**mm
         if (i.eq.1.or.i.eq.np) go to 10
         dl=dl+dl
         if ((i-2*(i/2)).eq.0) dl=dl+dl
 10      continue
         dc=dp(i)*dl
         if (dc.lt.0) then
            db=db+dc
         else if (dc.gt.0) then
            da=da+dc
         endif
         dc=dq(i)*dl
         if (dc.lt.0) then
            db=db+dc
         else if (dc.gt.0) then
            da=da+dc
         endif
      enddo
      da = dpas * (da + db) / 3.d0
      dc=exp(dpas)-1.d0
      db=d1*(d1+1.d0)*dc*exp((d1-1.d0)*dpas)
      db=dr(1)*(dr(2)**m)/db
      dc=(dr(1)**mm)*(1.d0+1.d0/(dc*(d1+1.d0)))/d1
      da=da+dc*(dp(1)+dq(1))-db*(dp(2)+dq(2))
      return
      end
