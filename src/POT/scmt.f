      subroutine scmt (verbse, iscmt, ecv, nph, nat, vclap,
     1                edens, edenvl, vtot, vvalgs, rmt, rnrm,qnrm,
     2                ixc, rhoint, vint, xmu, jumprm,
     3                xnferm, xnvmu, xnval,
     4                x0, ri, dx, xnatph, xion, iunf, iz,
     5                adgc, adpc, dgc,dpc, ihole,
     7                rat,iatph,iphat, lmaxsc, rhoval, xnmues, ok,
     8                rgrd, nohole, nscmt, icoul, ca1, rfms1, lfms1)

c     Finds new Fermi level (xmu), electron counts (qnrm)
c     and new valence densities (rhoval).

      implicit double precision (a-h, o-z)
      include '../HEADERS/dim.h'
      include '../HEADERS/const.h'

      logical verbse
c     input
      dimension dmagx(nrptx), dmag0(251)
      dimension vclap(251,0:nphx)
      dimension vtot (251,0:nphx), vvalgs (251,0:nphx)
      dimension xnvmu(0:lx,0:nphx+1), rmt(0:nphx),rnrm(0:nphx)
      dimension xnval (30,0:nphx)
      dimension qnrm(0:nphx), dq(0:nphx)
      dimension ri(nrptx), ri05(251), nr05(0:nphx)
      dimension xnatph(0:nphx), iz(0:nphx), xion(0:nphx)
      dimension rat(3,natx),iatph(0:nphx),iphat(natx), lmaxsc(0:nphx)
      real  rfms1
c     input and output
      dimension edens(251,0:nphx), edenvl(251,0:nphx)
      dimension rhoval(251,0:nphx+1)

c     work space
      dimension xnmues(0:lx, 0:nphx)
      complex gtr(0:lx, 0:nphx)
      dimension dum(nrptx), vtotph(nrptx),vvalph(nrptx)
      dimension dgc(251,30,0:nphx+1), dpc(251,30,0:nphx+1)
      dimension adgc(10,30,0:nphx+1), adpc(10,30,0:nphx+1)
      dimension dgcn(nrptx,30), dpcn(nrptx,30)
      complex*16 xrhoce(0:lx,0:nphx), xrhocp(0:lx,0:nphx)
      complex*16 xrhole(0:lx,0:nphx)
      complex*16 yrhoce(251,0:nphx), yrhocp(251,0:nphx)
      complex*16 yrhole(251,0:lx,0:nphx)
      complex*16 ph(lx+1, 0:nphx)
      integer iph
c     complex energy grid emg is decomposed into em and eref
      parameter (negx = 80)
      complex*16 emg(negx), em, eref, ee, ep, fl, fr, fxa
c     nflrx should be odd and defines the max of Im energy for
c     the countour
      parameter (nflrx = 17)
      dimension step(nflrx)
c     stuff from feff.f for rdinp, pathfinder and genfmt
      logical upok, ok
c      logical wnstar
c     Following passed to pathfinder, which is single precision.
      character*512 slog
      integer ient
      data ient /0/

c     save staff from rdinp, so no need to call it again
      save   ri05, ient

      upok = .false.
      idir = 1
      ient = ient + 1
      if (ient.eq.1) then
         do 15 i= 1,251
  15     ri05(i) = exp (-8.8d0+0.05d0*(i-1))
      endif

      if (verbse) then
         write (slog,10) iscmt, nscmt
 10      format('              SCF ITERATION NUMBER',i3,'  OUT OF',i3)
         call wlog(slog)

         call wlog (' Calculating energy and space dependent l-DOS.')
         call wlog (' It takes time ...')
      endif

c     initialize new valence density
      do 16 iph=0,nphx
      do 16 ir=1,251
  16  rhoval(ir,iph) = 0

c     polarization average in scmt and ldos

      call grids (ecv, xmu, negx, neg, emg, step, nflrx)

c     ie - is number of energy points calculated
      ie = 0
      xndifp = 0
      xndif  = 0
      ee = emg(1)
      ep = dble(ee)
      do 22 iph=0,nphx
      do 22 il=0,lx
        xrhoce(il,iph) = 0
        xnmues(il,iph) = 0
  22  continue
      do 23 iph=0,nphx
      do 23  ir = 1,251
  23  yrhoce(ir,iph) = 0
      iflr = nflrx
      iflrp = nflrx

c     Start the cycle over energy points (ie)
  25  continue
      ie = ie + 1

      do 29 iph = 0,nph
        do 860 il = 0,lx
  860     xrhocp(il,iph) = xrhoce(il,iph)
        do 870 i = 1,251
  870     yrhocp(i,iph) = yrhoce(i,iph)
  29  continue

      if (ie.eq.1 .or. mod(ie,20).eq.0) then
         if (verbse) then
            write(slog,30) ie, dble(ee)*hart
 30         format('     point # ', i3, '  energy = ', f7.3)
            call wlog(slog)
         endif
      endif

      do 100  iph = 0, nph

         do 35 i=1, 251
  35     dmag0(i) = 0.d0
cc       use spin-unpolarized case to get SCF. set dmagx to zero
cc       may want to replace dmag0 with dmag(1,iph) for spin-dependent
cc       extension of SCF procedure.
         call fixvar (rmt(iph),edens(1,iph),vtot(1,iph),dmag0,
     1                vint, rhoint, dx, rgrd, jumprm,
     2                vjump, ri, vtotph, dum, dmagx)
         if (mod(ixc,10) .ge.5) then
            if (jumprm .gt. 0) jumprm = 2
            call fixvar (rmt(iph), edenvl(1,iph), vvalgs(1,iph),
     1                dmag0, vint, rhoint, dx, rgrd , jumprm,
     2                vjump, ri, vvalph, dum, dmagx)
            if (jumprm .gt. 0) jumprm = 1
         endif

         call fixdsx (iph, dx, rgrd , dgc, dpc, dgcn, dpcn)
        jri = int((log(rmt(iph)) + x0) / rgrd) + 2
        jri1 = jri+1
        eref = vtotph(jri1)
        do 40 i = 1, jri1
  40    vtotph(i) = vtotph(i) - dble(eref)
        if (ixc.ge.5) then
           do 50 i = 1, jri1
  50       vvalph(i) = vvalph(i) - dble(eref)
        else
           do 60 i = 1, jri1
  60       vvalph(i) = vtotph(i)
        endif

         itmp = 0
         if (iph.eq.0 .and. nohole.lt.0) itmp = ihole
         call rholie( ri05, nr05(iph), rgrd, x0, ri, ee,
     2           ixc, rmt(iph), rnrm(iph),
     3           vtotph, vvalph, xnval(1,iph), dgcn, dpcn, eref,
     4           adgc(1,1,iph), adpc(1,1,iph), xrhole(0,iph),
     5          xrhoce(0,iph),yrhole(1,0,iph),yrhoce(1,iph),ph(1,iph),
     6           iz(iph), xion(iph), iunf, itmp,lmaxsc(iph))
  100 continue

c     Write out phases for fmsie
c     transform neg,emg to em,ne,eref first
       em= dble(ee)
       eref=dble(eref)-coni*dimag(ee)

cc    call fms for a cluster around central atom
       do 195 iph0 = 0,nph
       do 195 il = 0, lx
  195  gtr(il,iph0) = 0
      if (rfms1 .gt. 0) then
        if (lfms1 .ne. 0) then
          iph0 = 0
          call fmsie(verbse, iph0, nph, lmaxsc, ie,  em, eref, ph,
     1                rfms1, lfms1, nat, iphat, rat, gtr)
        else
          do 190 iph0 = 0, nph
  190     call fmsie(verbse, iph0, nph, lmaxsc, ie, em, eref, ph,
     1                rfms1, lfms1, nat, iphat, rat, gtr)
        endif
      endif

      xntot = 0
      fl = 0
      fr = 0
      do 300 iph = 0,nph
c       calculate density and integrated number of electrons in each
c       channel for each type of atoms density, etc., find xntot.
        call ff2g (gtr(0,iph), iph, ie, nr05(iph), xrhoce,
     1     xrhole(0,iph), xrhocp, ee, ep,
     2     yrhole(1,0,iph), yrhoce(1,iph), yrhocp(1,iph), rhoval(1,iph),
     3     xnmues(0,iph), xnatph(iph), xntot, iflr, iflrp, fl, fr, iunf)
  300 continue

      if (ie.ne.1) xndifp = xndif
      xndif = xntot - xnferm

c     decide on next energy point; there are nflrx floors, defined
c     by the magnitude of Im part. Each floor has it's height and
c     horizontal step to search for Fermi level associated with it.
c     The driver below will decide whether to go left or right on
c     the current floor, go one floor up or down.

      if ((ie.lt.neg .and. ient.gt.1) .or.
     1    (ient.eq.1.and.ie.lt.nflrx)) then
         ep = ee
         ee = emg(ie+1)
         if (ie.eq.neg-1) then
c          reset iflr variables
           iflrp = 2
           iflr  = 1
         endif
         goto 25
      elseif (ient.eq.1 .and. ie.eq.nflrx) then
         upok = .false.
         idir = 1
         if (xntot.gt. xnferm) idir = -1
         ep = ee
         ee = ee + idir * step(iflr)
         goto 25
      elseif (ient.gt.1 .and. ie.eq.neg) then
         upok = .true.
         iflrp = 1
         iflr  = 1
         idir = -1
         if (xntot.lt. xnferm) idir = 1
         ep = ee
         ee = ee + idir * step(iflr)
         goto 25
      else
c       check if the fermi level is found
        if (iflrp.eq.1 .and. iflr.eq.1 .and.
     1                xndifp*xndif .le. 0.d0) then
c          Fermi level is found ; do not goto 25
           if (xndif.eq.0) then
              xmunew = dble(ee)
              a=0
           else
              a = xndif/(xndif-xndifp)
              do 220 i = 1,4
                fxa = a*fl + (1-a)*fr
                bb = dimag(dcmplx((ep-ee)*(fr+fxa)/2 +
     $               coni*imag(ee)*(fr-fl)))
                xndif1 = xndif + a * bb
                a = a - xndif1 / bb
  220         continue
              xmunew = dble((1-a)*ee+a*ep)
           endif

c          add end cap corrections to the configuration and density
c          factor 2 for spin degeneracy
           do 250 iph = 0,nph
              do 230 il = 0,lx
               if (il.le.2 .or. iunf.ne.0) then
                fl = xrhocp(il,iph) * 2
                fr = xrhoce(il,iph) * 2
                fxa = a*fl + (1-a)*fr
                bb = dimag(dcmplx((ep-ee)*(fr+fxa)/2 +
     1               coni*dimag(dcmplx(ee))*(fr-fl)))
                xnmues(il,iph) = xnmues(il,iph) + a * bb
               endif
  230         continue
              do 240 ir = 1,nr05(iph)
                fl = yrhocp(ir,iph) * 2
                fr = yrhoce(ir,iph) * 2
                fxa = a*fl + (1-a)*fr
                bb = dimag(dcmplx((ep-ee)*(fr+fxa)/2 +
     1               coni*dimag(dcmplx(ee))*(fr-fl)))
                rhoval(ir,iph) = rhoval(ir,iph) + a * bb
  240         continue
  250      continue
        else
c          continue search ; goto 25 eventually
           if (iflr.eq.iflrp) then
c            previous step was gorizontal
             if (xndifp*xndif.le.0) then
c               need to step down
                upok =.false.
                iflrp = iflr
                iflr = iflr - 1
                ep = ee
                ee = dble(ee) + coni*4*step(iflr)
             elseif (abs(xndif).gt.10.d0*abs(xndif-xndifp)
     1          .and. upok) then
c               need to go up one floor since too far from fermi level
                iflrp = iflr
                if (iflr.lt.nflrx) then
                  iflr = iflr+1
                  ep = ee
                  ee = dble(ee) +  coni*4*step(iflr)
                else
                  ep = ee
                  ee = ee + idir* step(iflr)
                endif
             else
c               keep the same floor and direction
                ep = ee
                ee = ee + idir* step(iflr)
             endif
           else
c            previous step was up or down (vertical)
c            check the direction of search
             idir = -1
             if (xndif.lt.0) idir = 1
             iflrp = iflr
             ep = ee
             ee = ee + idir* step(iflr)
           endif
           goto 25
        endif
      endif
c     END of the loop over energy in comlex plane.
c     new fermi level and densities are calculated.

c     report configuration; repeat iteration if found bad counts.
      ok = .true.
      if (verbse) then
         call wlog('  Electronic configuration')
         call wlog('   iph    il      N_el')
      endif
 310  format (2i6, f9.3)
      do 320 ip= 0,nph
      do 320 il = 0,lx
         if (verbse) then
            write (slog,310) ip,il,xnmues(il,ip)
            call wlog(slog)
         endif
c        check that occupation numbers are consistent with those
c        set in getorb.f
         diff = abs(xnmues(il,ip) - xnvmu(il,ip))
         if (diff.gt.13.1d0 .or. (il.eq.2 .and. diff.gt. 9.1d0) .or.
     1   (il.eq.1 .and. diff.gt.5.1d0) .or.
     2   (il.eq.0 .and. diff.gt.1.95d0)) then
            call wlog (' Found bad counts.')
            write (slog,311) xnvmu(il,ip)
  311       format('  Occupation number in getorb is ', f9.3)
            call wlog(slog)
            call wlog ('  Will repeat this iteration ')
c            if (ient.gt.1) ok = .false.
         endif
 320  continue

c     if (.not. ok) then will restart SCF loop
      if (ok) then
         xmu = xmunew
c        find rhoval via Broyden algorithm
         call broydn( iscmt, ca1, nph, xnvmu,
     1         nr05 , xnatph, rnrm, qnrm, edenvl, rhoval, dq)

c        calculate new vclap - overlap coulomb potential
         call coulom (icoul, nph, nr05 , rhoval, edenvl, edens,
     2     nat, rat, iatph, iphat, rnrm, dq, iz, vclap)

c       update array edens
        do 350 ip=0,nph
           do 330 ir=1,nr05 (ip)
             edens(ir,ip)=edens(ir,ip)-edenvl(ir,ip)+rhoval(ir,ip)
  330      continue
           do 340 ir=nr05 (ip)+1,251
             edens(ir,ip)=0.0d0
             edenvl(ir,ip)=0.0d0
  340      continue
  350   continue
      endif

      return
      end
