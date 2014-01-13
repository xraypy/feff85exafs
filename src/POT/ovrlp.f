      subroutine ovrlp (iph, iphat, rat, iatph, novr, iphovr,
     1                nnovr, rovr, iz, nat, rho, dmag,
     2                rhoval, vcoul, edens, edenvl, vclap, rnrm)

c     Overlaps coulomb potentials and electron densities for current
c     unique potential
      implicit double precision (a-h, o-z)

      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'

      dimension iphat(natx)
      dimension rat(3,natx)
      dimension iatph(0:nphx)
      dimension novr(0:nphx)
      dimension iphovr(novrx,0:nphx)
      dimension nnovr(novrx,0:nphx)
      dimension rovr(novrx,0:nphx)
      dimension iz(0:nphx)
      dimension rho(251,0:nphx+1), dmag(251,0:nphx+1)
      dimension vcoul(251,0:nphx+1), rhoval(251,0:nphx+1)
      dimension edens(251,0:nphx), edenvl(251,0:nphx)
      dimension vclap(251,0:nphx)
      dimension rnrm(0:nphx)
c#mn
       external dist

c     start with free atom values for current atom
      do 100  i = 1, 251
         vclap(i,iph) = vcoul(i,iph)
         edens(i,iph) = rho  (i,iph)
         
cc       investigate effect of central atom spin only
c        if (iph.ge.1) dmag(i,iph) = 0.0

         edenvl(i,iph) = rhoval  (i,iph)
  100 continue

      if (novr(iph) .gt. 0)  then
         do 104  iovr = 1, novr(iph)
            rnn  = rovr(iovr,iph)
            ann  = nnovr(iovr,iph)
            infr = iphovr(iovr,iph)
            call sumax (rnn, ann, vcoul(1,infr), vclap(1,iph))
            call sumax (rnn, ann, rho  (1,infr), edens(1,iph))
            call sumax (rnn, ann, rho  (1,infr), edenvl(1,iph))
  104    continue
      else
c        Do overlapping from geometry with model atom iat
         iat = iatph(iph)

c        overlap with all atoms within r overlap max (rlapx)
c        12 au = 6.35 ang  This number pulled out of a hat...
         rlapx = 12
c        inat is Index of Neighboring ATom
         do 110  inat = 1, nat
c           don't overlap atom with itself
            if (inat .eq. iat)  goto 110

c           if neighbor is too far away, don't overlap it
            rnn = dist (rat(1,inat), rat(1,iat))
            if (rnn .gt. rlapx)  goto 110

            infr = iphat(inat)
            call sumax (rnn, one, vcoul(1,infr), vclap(1,iph))
            call sumax (rnn, one, rho  (1,infr), edens(1,iph))
            call sumax (rnn, one, rho  (1,infr), edenvl(1,iph))
cala        call sumax (rnn, one, rhoval(1,infr), edenvl(1,iph))
  110       continue
      endif

c     set norman radius
      call frnrm (edens(1,iph), iz(iph), rnrm(iph))

c     remember ratio dmag/edens , not dmag itself
      do 200 i = 1,251
        if (edens(i,iph) .gt. 0.d0) then
          dmag(i,iph) = dmag(i,iph) / edens(i,iph)
        else
          dmag(i,iph) = 0.d0
        endif
 200  continue

      return
      end
