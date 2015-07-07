      subroutine moveh (nat, iphat, iz, rath)
c    Increase length of bonds with hydrogen atoms
c    Move Hydrogens for potentials. Otherwise MT geometry is screwed up.
      implicit double precision (a-h, o-z)
      include '../HEADERS/dim.h'

c     input is everything; output is modified  atomic coordinates (rath) 
c     nat is number of atoms in cluster
      dimension iphat(natx),  iz(0:nphx)
      dimension rath(3,natx)

      do 970 iat = 1, nat
        if (iz(iphat(iat)) .eq. 1) then
c         find the nearest atom A, units for rat are bohr.
          rah = 100
          ia = 0
          do i = 1,nat
              rattmp = dist(rath(1,iat), rath(1,i) )
              if (rattmp.lt. rah .and. i.ne. iat) then
                 ia = i
                 rah = rattmp
              endif
          enddo
          if (iz(iphat(ia)).eq.1) goto 970

c         set max distance as function of rah ( set by calculations
c         for H2O and GeH4)
          ratmax = rah + 4.d0/rah**2 

c        find shortest AB bond (neither A or B are H)
          rab = 10
          ib = 0
          do i = 1,nat
              rattmp = dist(rath(1,ia), rath(1,i))
              if (i.ne.ia .and. iz(iphat(i)).ne.1 .and.
     1            rab.gt.rattmp) then
                 rab = rattmp
                 ib = i
              endif
          enddo
          if (rab.lt.ratmax) ratmax = 0.95d0*rab + 0.05d0*rah
          if (rah .gt. ratmax) goto 970

c         increase rah to ratmax and check that A is still closest to H
          ratmin = rah
  960     do i = 1,3
           rath(i,iat)=rath(i,ia)+ratmax/ratmin*(rath(i,iat)-rath(i,ia))
          enddo
          rbh = 10
          ib = 0
          do i = 1,nat
              rattmp = dist(rath(1,iat), rath(1,i))
              if (i.ne.iat .and. rbh.gt.rattmp) then
                 rbh = rattmp
                 ib = i
              endif
          enddo

          if (ia.ne.ib) then
             rab = dist(rath(1,ia),rath(1,ib))
             rattmp = ratmax*rab**2/(ratmax**2+rab**2-rbh**2)
             ratmin = ratmax
             ratmax = 0.95d0*rattmp +0.05d0*rah
             goto 960
          endif
        endif
  970 continue

      return
      end
