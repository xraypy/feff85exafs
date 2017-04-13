      subroutine ovrlp (iph, iphat, rat, iatph, ifrph, novr,
     1                  iphovr, nnovr, rovr, iz, nat, rho, vcoul,
     2                  edens, vclap, rnrm)

c     Overlaps coulomb potentials and electron densities for current
c     unique potential
      implicit double precision (a-h, o-z)

      include 'const.h'
      include 'dim.h'

      dimension iphat(natx)
      dimension rat(3,natx)
      dimension iatph(0:nphx)
      dimension ifrph(0:nphx)
      dimension novr(0:nphx)
      dimension iphovr(novrx,0:nphx)
      dimension nnovr(novrx,0:nphx)
      dimension rovr(novrx,0:nphx)
      dimension iz(0:nfrx)
      dimension rho(251,0:nfrx)
      dimension vcoul(251,0:nfrx)
      dimension edens(nrptx,0:nphx)
      dimension vclap(nrptx,0:nphx)
      dimension rnrm(0:nphx)

c     find out which free atom we're dealing with
      ifr = ifrph(iph)

c     start with free atom values for current atom
      do 100  i = 1, 250
         vclap(i,iph) = vcoul(i,ifr)
         edens(i,iph) = rho  (i,ifr)
  100 continue

      if (novr(iph) .gt. 0)  then
         do 104  iovr = 1, novr(iph)
            rnn  = rovr(iovr,iph)
            ann  = nnovr(iovr,iph)
            infr = ifrph(iphovr(iovr,iph))
            call sumax (250, rnn, ann, vcoul(1,infr), vclap(1,iph))
            call sumax (250, rnn, ann, rho  (1,infr), edens(1,iph))
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

            infr = ifrph(iphat(inat))
            call sumax (250, rnn, one, vcoul(1,infr), vclap(1,iph))
            call sumax (250, rnn, one, rho  (1,infr), edens(1,iph))
  110       continue
      endif

c     set norman radius
      call frnrm (edens(1,iph), iz(ifr), rnrm(iph))

      return
      end
