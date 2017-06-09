      subroutine istprm (nph, nat, iphat, rat, iatph, xnatph,
     1                   novr, iphovr, nnovr, rovr, folp, edens,
     2                   vclap, vtot, imt, inrm, rmt, rnrm, 
     2                   rhoint,
     3                   vint, rs, xf, xmu, rnrmav, intclc)

c     Finds interstitial parameters, rmt, vint, etc.
      implicit double precision (a-h, o-z)

      include 'const.h'
      include 'dim.h'

      dimension iphat(natx)
      dimension rat(3,natx)
      dimension iatph(0:nphx)
      dimension xnatph(0:nphx)
      dimension novr(0:nphx)
      dimension iphovr(novrx,0:nphx)
      dimension nnovr(novrx,0:nphx)
      dimension rovr(novrx,0:nphx)
      dimension folp(0:nphx)
      dimension edens(nrptx,0:nphx)
      dimension vclap(nrptx,0:nphx)
      dimension vtot (nrptx,0:nphx)
      dimension imt(0:nphx)
      dimension inrm(0:nphx)
      dimension rmt(0:nphx)
      dimension rnrm(0:nphx)

       character*128 messag
c     intclc = 0, average evenly over all atoms
c              1, weight be lorentzian, 1 / (1 + 3*x**2), x = r/rnn,
c                 r   = distance to central atom,
c                 rnn = distance of near neighbor to central atom

c Find muffin tin radii.  We'll find rmt based on norman prescription,
c ie, rmt(i) = R * folp * rnrm(i) / (rnrm(i) + rnrm(j)),
c a simple average
c based on atoms i and j.  We average the rmt's from each pair of
c atoms, weighting by the volume of the lense shape formed by the
c overlap of the norman spheres.
c NB, if folp=1, muffin tins touch without overlap, folp>1 gives
c overlapping muffin tins.
c
c rnn is distance between sphere centers
c rnrm is the radius of the norman sphere
c xl_i is the distance to the plane containing the circle of the
c    intersection
c h_i  = rnrm_i - xl_i is the height of the ith atom's part of
c    the lense
c vol_i = (pi/3)*(h_i**2 * (3*rnrm_i - h_i))
c
c xl_i = (rnrm_i**2 - rnrm_j**2 + rnn**2) / (2*rnn)

      do 140  iph = 0, nph
         voltot = 0
         rmtavg = 0
         if (novr(iph) .gt. 0)  then
c           Overlap explicitly defined by overlap card
            do 124  iovr = 1, novr(iph)
               rnn  = rovr(iovr,iph)
               inph = iphovr(iovr,iph)
c              Don't avg if norman spheres don't overlap
               if (rnrm(iph)+rnrm(inph) .le. rnn)  goto 124
               voltmp = calcvl (rnrm(iph), rnrm(inph), rnn)
               voltmp = voltmp + calcvl (rnrm(inph), rnrm(iph), rnn)
               rmttmp = rnn * folp(iph) * rnrm(iph) /
     1                  (rnrm(iph) + rnrm(inph))
               ntmp = nnovr(iovr,iph)
               rmtavg = rmtavg + rmttmp*voltmp*ntmp
               voltot = voltot + voltmp*ntmp
  124       continue
         else
            iat = iatph(iph)
            do 130  inat = 1, nat
               if (inat .eq. iat)  goto 130
               rnn = dist (rat(1,inat), rat(1,iat))
               inph = iphat(inat)
c              Don't avg if norman spheres don't overlap
               if (rnrm(iph)+rnrm(inph) .lt. rnn)  goto 130
               voltmp = calcvl (rnrm(iph), rnrm(inph), rnn)
               voltmp = voltmp + calcvl (rnrm(inph), rnrm(iph), rnn)
               rmttmp = rnn * folp(iph) * rnrm(iph) /
     1                  (rnrm(iph) + rnrm(inph))
               rmtavg = rmtavg + rmttmp*voltmp
               voltot = voltot + voltmp
  130       continue
         endif
         if (rmtavg .le. 0)  then
            write(messag, 132) iat, iph
            call echo(messag)
  132       format (' WARNING: NO ATOMS CLOSE ENOUGH TO OVERLAP ATOM',
     1              i5, ',  UNIQUE POT', i5, '!!', /,
     2              '          Rmt set to Rnorman.  May be error in ',
     3              'input file.')
            rmt(iph) = rnrm(iph)
         else
            rmt(iph) = rmtavg / voltot
         endif
  140 continue

c     Need potential with ground state xc, put it into vtot
      do 160  iph = 0, nph
         call sidx (edens(1,iph), 250, rmt(iph), rnrm(iph),
     1              imax, imt(iph), inrm(iph))
         do 150  i = 1, imax
            rs = (edens(i,iph)/3)**(-third)
c           vhedbr from Von Barth Hedin paper, 1971
            vhedbr = -1.22177412/rs - 0.0504*log(30/rs + 1)
            vtot(i,iph) = vclap(i,iph) + vhedbr
  150    continue
  160 continue

c     What to do about interstitial values?
c     Calculate'em for all atoms, print'em out for all unique pots along
c     with derivative quantities, like fermi energy, etc.
c     Interstitial values will be average over all atoms in problem.

c     rnrmav is averge norman radius,
c     (4pi/3)rnrmav**3 = (sum((4pi/3)rnrm(i)**3)/n, sum over all atoms
c     in problem
      rnrmav = 0
      xn = 0
      rs = 0
      vint   = 0
      rhoint = 0
c     volint is total interstitial volume
      volint = 0

      do 170  iph = 0, nph
c        Use all atoms
         call istval (vtot(1,iph), edens(1,iph), rmt(iph), imt(iph),
     2                rnrm(iph), inrm(iph), vintx, rhintx, ierr)
c        if no contribution to interstitial region, skip this unique pot
         if (ierr .ne. 0)  goto 170
         call fermi (rhintx, vintx, xmu, rs, xf)
c        (factor 4pi/3 cancel in numerator and denom, so leave out)
         volx = (rnrm(iph)**3 - rmt(iph)**3)
         if (volx .le. 0)  goto 170
         volint = volint + volx * xnatph(iph)
         vint   = vint   + vintx * volx * xnatph(iph)
         rhoint = rhoint + rhintx* volx * xnatph(iph)
  170 continue
c     If no contribution to interstitial from any atom, die.
      if (volint .le. 0)  then
         call echo(' No interstitial density.  Check input file.')
         call fstop(' at ISTPRM')
      endif
      vint   = vint   / volint
      rhoint = rhoint / volint
      call fermi (rhoint, vint, xmu, rs, xf)
      do 180  iph = 0, nph
         rnrmav = rnrmav + xnatph(iph) * rnrm(iph)**3
         xn = xn + xnatph(iph)
  180 continue
      rnrmav = (rnrmav/xn) ** third


      return
      end

      double precision function calcvl (r1, r2, r)
      implicit double precision (a-h, o-z)
      include 'const.h'
      xl = (r1**2 - r2**2 + r**2) / (2*r)
      h = r1 - xl
      calcvl = (pi/3) * h**2 * (3*r1 - h)
      return
      end
