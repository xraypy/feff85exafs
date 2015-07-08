      subroutine istprm ( nph, nat, iphat, rat, iatph, xnatph,
     1                novr, iphovr, nnovr, rovr, folp, folpx, iafolp,
     1                edens, edenvl, idmag,
     2                dmag, vclap, vtot, vvalgs, imt, inrm, rmt, rnrm, 
     2                ixc, rhoint, vint, rs, xf, xmu, xmunew,
     3                rnrmav, qtotel, inters, totvol)

c     Finds interstitial parameters, rmt, vint, etc.
      implicit double precision (a-h, o-z)

      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'

      dimension iphat(natx)
      dimension rat(3,natx)
      dimension iatph(0:nphx)
      dimension xnatph(0:nphx)
      dimension novr(0:nphx)
      dimension iphovr(novrx,0:nphx)
      dimension nnovr(novrx,0:nphx)
      dimension rovr(novrx,0:nphx)
      dimension folp(0:nphx), folpx(0:nphx)
      dimension edens(251,0:nphx), edenvl(251,0:nphx)
      dimension dmag(251,0:nphx+1)
      dimension vclap(251,0:nphx)
      dimension vtot (251,0:nphx), vvalgs (251,0:nphx)
      dimension imt(0:nphx)
      dimension inrm(0:nphx)
      dimension rmt(0:nphx)
      dimension rnrm(0:nphx)
      parameter (big = 5000)
      character*512 slog
      logical lnear
      dimension lnear(0:nphx), inn(0:nphx), rnnmin(0:nphx)
c#mn
      external dist

c     work space for linear algebra
      dimension ri(251)
      parameter (novp=40)
      complex cmovp(novp*(nphx+1)+1,novp*(nphx+1)+1)
      integer ipiv(novp*(nphx+1)+1)
      save lnear

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

c     find rmt from rnrm only on first call of istprm (rmt(0)=-1)
      if (rmt(0).le.0.0d0) then
      do 10 iph=0,nph
  10  lnear(iph)=.false.
      do 140  iph = 0, nph
         voltot = 0
         rmtavg = 0
         inrm(iph) = ii(rnrm(iph))
         if (novr(iph) .gt. 0)  then
c           Overlap explicitly defined by overlap card
            rnear = big
            inters = mod(inters,6)
c           use Norman prescription only in this case

            do 124  iovr = 1, novr(iph)
               rnn  = rovr(iovr,iph)
               inph = iphovr(iovr,iph)
               if (rnn .le. rnear) then
                  rnear = rnn
                  rnnmin(iph) = rnn
                  inn(iph) = inph
               endif
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
            rnear = big
            rmt(iph) = big
            do 130  inat = 1, nat
               if (inat .eq. iat)  goto 130
               rnn = dist (rat(1,inat), rat(1,iat))
               inph = iphat(inat)
               if (rnn .le. rnear) then
                  rnear = rnn
                  rnnmin(iph) = rnn
                  inn(iph) = inph
               endif
c              Don't avg if norman spheres don't overlap
               if (rnrm(iph)+rnrm(inph) .lt. rnn)  goto 130

               if (inters.lt.6) then
c                Norman prescription
                 voltmp = calcvl (rnrm(iph), rnrm(inph), rnn)
                 voltmp = voltmp + calcvl (rnrm(inph), rnrm(iph), rnn)
                 rmttmp = rnn * folp(iph) * rnrm(iph) /
     1                  (rnrm(iph) + rnrm(inph))
                 rmtavg = rmtavg + rmttmp*voltmp
                 voltot = voltot + voltmp
               else
c                Matching point prescription
                 do 125 i=inrm(iph),1,-1
                   j=ii(rnn-rnrm(iph))
                   if (vclap(i,iph).le.vclap(j,inph)) then
                     d1 = (vclap(i+1,iph)-vclap(i,iph))/(rr(i+1)-rr(i))
                     d2 =(vclap(j,inph)-vclap(j-1,inph))/(rr(j)-rr(j-1))
                     rmtavg = rr(i) + 
     1               (vclap(j,inph)+d2*(rnn-rr(i)-rr(j))-vclap(i,iph))
     2               /(d1+d2)
                     goto 127
c                    exit from the loop
                   endif
  125            continue
  127            continue
                 if (rmtavg.lt.rmt(iph)) rmt(iph) = rmtavg
               endif
  130       continue
         endif

c        special situation if rnrm is too close or larger than
c        the nearest neighbor distance
         if (rnrm(iph).ge.rnear) lnear(iph) = .true.

         if (rmtavg .le. 0)  then
            write(slog,132) iat, iph
            call wlog(slog)
  132       format (' WARNING: NO ATOMS CLOSE ENOUGH TO OVERLAP ATOM',
     1              i5, ',  UNIQUE POT', i5, '!!  ', 
     2              'Rmt set to Rnorman.  May be error in ',
     3              'input file.')
            rmt(iph) = rnrm(iph)
         elseif(inters.lt.6) then
c           Norman prescription
            rmt(iph) = rmtavg / voltot
            if (rmt(iph) .ge. rnear)  then
c              print*,iph, rmt(iph), rnear
               call wlog(' Rmt >= distance to nearest neighbor.  ' //
     1            'Not physically, meaningful.')
               call wlog(' FEFF may crash.  Look for error in ATOM '//
     1            'list or OVERLAP cards.')
            endif
            if (rnrm(iph) .ge. rnear) then
              imax = ii(rnear) - 1
c             begin until loop
 133            if (vclap(imax,iph).lt.vclap(imax+1,iph)) goto 134
                imax = imax-1
                goto 133
c             end of until loop
 134          continue
              rmt(iph) = exp(xx(imax)) - 0.0001d0
            endif
         endif

  140 continue

c     set maximum value for folp(iph) if AFOLP is in use
c     LMTO lore says no more than 15% overlap
c     do 144 iph = 0, nph
c 144 folpx(iph) = 1.15
c     already done in pot.f

      do 145 iph = 0, nph
         if (iafolp. gt. 0 ) then
            temp = 0.2d0 + 0.8d0 * rnrm(iph) / rmt(iph)
         else
            temp = 0.3d0 + 0.7d0 * rnrm(iph) / rmt(iph)
         endif
         if (temp.lt.folpx(iph)) folpx(iph) = temp
         temp = rnnmin(iph)/rmt(iph)/1.06d0
         if (temp.lt.folpx(iph)) folpx(iph) = temp
         temp = exp( -(novp-3)*0.05d0)
c      make sure that with given folpx(iph) the construction
c      of the overlapping matrix in movrlp will not fail
         if (lnear(iph)) then
c           lnear=.true. only when hydrogens are present in the system.
c           want to scale both rmt for iph and inn, so that overlapping
c           matrix calculations will not fail
            temp = rnnmin(iph) / (rmt(iph)*1.05d0 + temp*rmt(inn(iph)))
            if (temp.lt.folpx(iph)) folpx(iph) = temp
            if (temp.lt.folpx(inn(iph))) folpx(inn(iph)) = temp
         else
            temp = (rnnmin(iph) - rnrm(iph))/ (temp*rmt(inn(iph)))
            if (temp.lt.folpx(inn(iph))) folpx(inn(iph)) = temp
         endif
  145 continue

      endif
c     end of finding rmt from rnrm on first call of istprm.

c     Need potential with ground state xc, put it into vtot
      do 160  iph = 0, nph
         call sidx (edens(1,iph), 250, rmt(iph), rnrm(iph),
     1              imax, imt(iph), inrm(iph))
         do 150  i = 1, imax
            if (edens(i,iph).le.0) then
             if(mod(i,10).eq.0) then
               write(slog, 149) 'negative dens ', i,iph
  149          format (a, 2i3)
               call wlog(slog)
             endif
             rs = 100
             xmag=1.0d0
            else
              rs = (edens(i,iph)/3)**(-third)
c     spin dependent xc potential for ground state from Von Barth, Hedin
c     J.Phys.C:Solid State Phys., 5, 1629 (1972).
c     xmag/2 -fraction of spin up or down, depending on sign in renorm.f
c     put xmag = 1.0 to calculate cmd with external potential difference
              xmag = 1.0d0 + idmag*dmag(i,iph)
            endif
c           wrong for ferromagnets, need to overlap dmag(i)

c           vvbh from Von Barth Hedin paper, 1971
            call vbh(rs,xmag,vvbh)
            vtot(i,iph) = vclap(i,iph) + vvbh

            if (mod(ixc,10).eq.5) then
              rsval = 10.0d0
              if (edenvl(i,iph) .gt. 0.00001d0) 
     1           rsval = (edenvl(i,iph)/3)**(-third)
              if (rsval.gt.10.0d0) rsval = 10.0d0
              xmagvl = 1.0d0 + idmag * dmag(i,iph) 
     1                      * edens(i,iph) / edenvl(i,iph)
              call vbh(rsval,xmagvl,vvbhvl)
              vvalgs(i,iph) = vclap(i,iph) + vvbhvl
            elseif (mod(ixc,10) .ge. 6) then
              if (edens(i,iph).le.edenvl(i,iph)) then
                 rscore =101.0d0
              else
                 rscore = ((edens(i,iph)-edenvl(i,iph)) / 3)**(-third)
              endif
              rsmag = (edens(i,iph)*(1+idmag*dmag(i,iph)) / 3)**(-third)
              xfmag = fa/rsmag
              call edp(rscore,xfmag,vrdh)
              vvalgs(i,iph) = vclap(i,iph) + vvbh - vrdh
            else
              vvalgs(i,iph) = 0.d0
            endif
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
c     volint is total interstitial volume
      volint = 0
      do 180  iph = 0, nph
         rnrmav = rnrmav + xnatph(iph) * rnrm(iph)**3
         volint=volint-xnatph(iph) * rmt(iph)**3
         xn = xn + xnatph(iph)
  180 continue
      if (totvol.le.0.0d0) then
         volint=4*pi/3 *(volint+rnrmav)
      else
         volint=4*pi/3 *volint + totvol
      endif
c     volume of lenses from overlapping mt spheres is added in movrlp.
      rnrmav = (rnrmav/xn) ** third

      rs = 0
      vint   = 0
      rhoint = 0
      rsval = 0

      call movrlp(nph, nat, iphat, rat, iatph, xnatph,
     1            novr, iphovr, nnovr, rovr,
     2            imt, rmt, rnrm, ri, lnear,
     3            cmovp,ipiv, volint,inters)

c     If no contribution to interstitial from any atom, die.
      if (volint .le. 0)  then
         call wlog(' No interstitial density.  Check input file.')
         call par_stop('ISTPRM')
      endif

c     find interstitial density

      call ovp2mt(nph, edens, 0, qtotel, ri, xnatph, lnear,
     1            inrm, imt, rnrm, rmt, cmovp,ipiv, rhoint,inters)
      rhoint = 4*pi * rhoint / volint

      if (ixc.ge.5) then
c        find valence potential inside mt sphere (vintvl -dummy)
         call ovp2mt(nph, vvalgs, 1, qtotel, ri, xnatph, lnear,
     1           inrm, imt, rnrm, rmt, cmovp, ipiv, vintvl,inters)
      endif

c     find potential inside mt sphere and vint
      call ovp2mt(nph, vtot, 1, qtotel, ri, xnatph, lnear,
     1            inrm, imt, rnrm, rmt, cmovp, ipiv, vint,inters)

      if (vint.ge.xmu) then
        write(slog,'(a)')
     1  ' WARNING:interstitial level found above Fermi level'
        call wlog(slog)
        write(slog,'(a)')
     1  '  Results may be unreliable. See manual for details'
        call wlog(slog)
        vint = xmu - 0.05d0
        call ovp2mt(nph, vtot, 2, qtotel, ri, xnatph, lnear,
     1            inrm, imt, rnrm, rmt, cmovp, ipiv, vint,inters)
      endif
      call fermi (rhoint, vint, xmunew, rs, xf)

      return
      end

      double precision function calcvl (r1, r2, r)
      implicit double precision (a-h, o-z)
      include '../HEADERS/const.h'
      xl = (r1**2 - r2**2 + r**2) / (2*r)
      h = r1 - xl
      calcvl = (pi/3) * h**2 * (3*r1 - h)
      return
      end
