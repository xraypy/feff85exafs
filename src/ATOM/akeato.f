      double precision function akeato(i,j,k)
c     angular coefficient by the direct coulomb integral fk for orbitals
c     i and j
      implicit double precision (a-h,o-z)
      common/mulabk/afgk
      dimension afgk(30,30,0:3)
 
c     afgk angular coefficients by integrales fk and gk
c        coefficient of integral fk(i;j) is in  afgk(min,max)
c        and that of integral gk(i;j) is in  afgk(max,min)
c        max=max(i,j) min=min(i,j)
 
      if (i .le. j) then 
         akeato=afgk(i,j,k/2)
      else
         akeato=afgk(j,i,k/2)
      endif
      return

      entry bkeato (i,j,k)
c     angular coefficient at the exchange coulomb integral gk
 
      bkeato=0.0d0
      if (i .lt. j) then
         bkeato=afgk(j,i,k/2)
      elseif (i.gt.j) then
         bkeato=afgk(i,j,k/2)
      endif
      return
      end
