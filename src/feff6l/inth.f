      subroutine inth (dp,dq,dv,dr)
c
c integration by the 5-point method of adams for the large
c component dp and the small component dq at the point dr;
c dv being the potential at this point
c **********************************************************************
      implicit double precision (a-h,o-z)
      save
      common /ps1/ dep(5), deq(5), db, dvc, dsal, dk, dm
c
c dep,deq the derivatives of dp and dq; db=energy/dvc;
c dvc=speed of light in atomic units; dsal=2.*dvc; dk=kappa quantum numb
c dm=exponential step/720.
c dkoef1=405./502., dkoef2=27./502.
c **********************************************************************
      data dkoef1 /.9462151394422310/, dkoef2 /.5378486055776890d-1/
      dpr=dp+dm*((251.0*dep(1)+2616.0*dep(3)+1901.0*dep(5))-(1274.0
     1 *dep(2)+2774.0*dep(4)))
      dqr=dq+dm*((251.0*deq(1)+2616.0*deq(3)+1901.0*deq(5))-(1274.0
     1 *deq(2)+2774.0*deq(4)))
      do 10 i=2,5
      dep(i-1)=dep(i)
   10 deq(i-1)=deq(i)
      dsum=(db-dv/dvc)*dr
      dep(5)=-dk*dpr+(dsal*dr+dsum)*dqr
      deq(5)=dk*dqr-dsum*dpr
      dp=dp+dm*((106.0*dep(2)+646.0*dep(4)+251.0*dep(5))-(19.0*dep(1
     1 )+264.0*dep(3)))
      dq=dq+dm*((106.0*deq(2)+646.0*deq(4)+251.0*deq(5))-(19.0*deq(1
     1 )+264.0*deq(3)))
      dp=dkoef1*dp+dkoef2*dpr
      dq=dkoef1*dq+dkoef2*dqr
      dep(5)=-dk*dp+(dsal*dr+dsum)*dq
      deq(5)=dk*dq-dsum*dp
      return
      end
