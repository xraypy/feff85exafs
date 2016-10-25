      subroutine scmtmp (verbse, npr, iscmt, ecv, nph, nat, vclap,
     2                edens, edenvl, vtot, vvalgs, rmt, rnrm,qnrm,
     2                ixc, rhoint, vint, xmu, jumprm,
     3                xnferm, xnvmu, xnval,
     4                x0, ri, dx, xnatph, xion, iunf, iz,
     5                adgc, adpc, dgc,dpc, ihole,
     7                rat,iatph,iphat, lmaxsc, rhoval, xnmues, ok,
     8                rgrd, nohole, nscmt, icoul, ca1, rfms1, lfms1,
     9                gtr, xrhole, xrhoce, yrhole, yrhoce )

c     Finds new Fermi level (xmu), electron counts (qnrm)
c     and new valence densities (rhoval).

      implicit double precision (a-h, o-z)
      include '../HEADERS/dim.h'
      include '../HEADERS/const.h'
      include '../HEADERS/parallel.h'
      real*8 wall_commend, wall_commst
      Parameter (Maxprocs = 1)

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
      dimension dum(nrptx), vtotph(nrptx),vvalph(nrptx)
      dimension dgc(251,30,0:nphx+1), dpc(251,30,0:nphx+1)
      dimension adgc(10,30,0:nphx+1), adpc(10,30,0:nphx+1)
      dimension dgcn(nrptx,30), dpcn(nrptx,30)
      complex*16 ph(lx+1, 0:nphx)
      complex*16 xrhocp(0:lx,0:nphx), yrhocp(251,0:nphx)

c     special dimension for MPI
c     Maxprocs = max number of processors for parallel execution
c--   This is set in parallel.h
c     npr - actual number of processors is passed to this subroutine
      complex gtr(0:lx, 0:nph, Maxprocs)
      complex*16 xrhoce(0:lx,0:nph,Maxprocs)
      complex*16 xrhole(0:lx,0:nph,Maxprocs)
      complex*16 yrhoce(251,0:nph,Maxprocs)
      complex*16 yrhole(251,0:lx,0:nph,Maxprocs)

      integer iph
c     complex energy grid emg is decomposed into em and eref
      parameter (negx = 80)
      complex*16 emg(negx), em, eref, ee, ep, fl, fr, fxa
c     nflrx should be odd and defines the max of Im energy for
c     the countour
      parameter (nflrx = 17)
      dimension step(nflrx)
c     stuff from feff.f for rdinp, pathfinder and genfmt
      logical ok
c      logical wnstar
c     Following passed to pathfinder, which is single precision.
      character*512 slog
      integer ient
      data ient /0/

c     save stuff from rdinp, so no need to call it again
      save   ri05, ient

      xndif = 0.0d0
      xndifp = 0.0d0
      ient = ient + 1
      if (ient.eq.1) then
         xmu = -0.25d0
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
      ietot0 = 1
      ee = emg(1)
      ep = dble(ee)
      do 21 ipr=1,Maxprocs
      do 21 iph=0,nph
      do 21 il=0,lx
  21    xrhoce(il,iph, ipr) = 0
      do 22 iph=0,nphx
      do 22 il=0,lx
  22    xnmues(il,iph) = 0
      do 23 ipr=1,Maxprocs
      do 23 iph=0,nph
      do 23  ir = 1,251
  23  yrhoce(ir,iph,ipr) = 0
      iflr = nflrx
      iflrp = nflrx

      nproc = min(npr, Maxprocs)
      n1 = 1
      n2 = min(neg, nproc)

c     Start the cycle over energy points (ie)
  25  continue

c     slow loop for MPI execution
      ie = this_process + n1
      ietot = ietot0 + this_process
      if (ie .gt. n2) go to 200
        ipr = 1 + ie - n1
        if (worker) par_type = 3

c       print *,'process n1 n2 ietot',this_process,n1,n2,ietot

        if (ietot.eq.1 .or. mod(ietot,20).eq.0) then
           if (verbse) then
              write(slog,30) ietot, dble(emg(ie))*hart
 30           format('     point # ', i3, '  energy = ', f7.3)
              call wlog(slog)
           endif
        endif

        do 100  iph = 0, nph
          do 35 i=1, 251
  35      dmag0(i) = 0.d0
cc        use spin-unpolarized case to get SCF. set dmagx to zero
cc        may want to replace dmag0 with dmag(1,iph) for spin-dependent
cc        extension of SCF procedure.
          call fixvar (rmt(iph),edens(1,iph),vtot(1,iph),dmag0,
     1                  vint, rhoint, dx, rgrd, jumprm,
     2                  vjump, ri, vtotph, dum, dmagx)
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
  40      vtotph(i) = vtotph(i) - dble(eref)
          if (ixc.ge.5) then
            do 50 i = 1, jri1
  50        vvalph(i) = vvalph(i) - dble(eref)
          else
            do 60 i = 1, jri1
  60        vvalph(i) = vtotph(i)
          endif

           itmp = 0
           if (iph.eq.0 .and. nohole.lt.0) itmp = ihole
           call rholie( ri05, nr05(iph), rgrd, x0, ri, emg(ie), ixc,
     2           rmt(iph), rnrm(iph), vtotph, vvalph, xnval(1,iph),
     3           dgcn, dpcn, eref, adgc(1,1,iph), adpc(1,1,iph),
     4           xrhole(0,iph,ipr), xrhoce(0,iph,ipr),
     5           yrhole(1,0,iph,ipr), yrhoce(1,iph,ipr),
     6           ph(1,iph), iz(iph), xion(iph), iunf, itmp,lmaxsc(iph))
  100   continue

c       Write out phases for fmsie
c       transform neg,emg to em,ne,eref first
        em= dble(emg(ie))
        eref=dble(eref)-coni*dimag(emg(ie))

cc      call fms for a cluster around central atom
        do 115 iph0 = 0,nph
        do 115 il = 0, lx
  115   gtr(il,iph0,ipr) = 0
        if (rfms1 .gt. 0) then
          if (lfms1 .ne. 0) then
            iph0 = 0
c           set logic to call yprep on every processor
            lfms = lfms1
            if (ietot0.eq.1) lfms = 2
            call fmsie(verbse, iph0, nph, lmaxsc, ietot, em, eref, ph,
     1           rfms1, lfms, nat, iphat, rat, gtr(0,0,ipr))
          else
            do 190 iph0 = 0, nph
  190       call fmsie(verbse, iph0, nph, lmaxsc, ietot, em, eref, ph,
     1           rfms1, lfms1, nat, iphat, rat, gtr(0,0,ipr))
          endif
        endif
  200 continue
c     end of slow loop for MPI execution

      ietot0 = ietot0 + n2 - n1 + 1
      if (worker) par_type = 2

      ixl = (lx + 1) * (nph + 1)
      ixly = ixl * 251
      ixlc = (nph + 1) * 251
      if (nproc .gt. 1) then
        call seconds(wall_commst)
        if (worker .and. (ie .le. n2)) then
c-- Send pointers for gtr buffer to master
          call par_send_int(ixl,1,0,this_process)
          call par_send_int(ixly,1,0,this_process)
          call par_send_int(ixlc,1,0,this_process)
c-- Send buffer
          if (ixl .ne. 0) then
            call par_send_cmplx(gtr(0,0,ipr),ixl,0,this_process)
            call par_send_dc(xrhoce(0,0,ipr),ixl, 0, this_process)
            call par_send_dc(xrhole(0,0,ipr),ixl, 0, this_process)
          endif
          if (ixly .ne. 0)
     .      call par_send_dc(yrhole(1,0,0,ipr),ixly, 0, this_process)
          if (ixlc .ne. 0)
     .      call par_send_dc(yrhoce(1,0,ipr),ixlc, 0, this_process)
        else if (master) then
          do i = 1,n2-n1
c-- Receive pointers for gtr buffer from i
            call par_recv_int(ixl,1,i,i)
            call par_recv_int(ixly,1,i,i)
            call par_recv_int(ixlc,1,i,i)
c-- Receive buffer from i
            if (ixl .ne. 0) then
              call par_recv_cmplx(gtr(0,0,i+1),ixl,i,i)
              call par_recv_dc(xrhoce(0,0,i+1),ixl,i,i)
              call par_recv_dc(xrhole(0,0,i+1),ixl,i,i)
            endif
            if (ixly .ne. 0)
     .        call par_recv_dc(yrhole(1,0,0,i+1),ixly,i,i)
            if (ixlc .ne. 0)
     .        call par_recv_dc(yrhoce(1,0,i+1),ixlc,i,i)
          enddo
        endif
c-- Broadcast gtr
c-- Needed here since we aren't done yet
        ilen = ixl * (n2 - n1 + 1)
        ileny = ilen * 251
        ilenc = (nph + 1) * (n2 - n1 + 1) * 251
        call par_bcast_cmplx(gtr(0,0,1),ilen,0)
        call par_bcast_dc(xrhoce(0,0,1),ilen,0)
        call par_bcast_dc(xrhole(0,0,1),ilen,0)
        call par_bcast_dc(yrhole(1,0,0,1),ileny,0)
        call par_bcast_dc(yrhoce(1,0,1),ilenc,0)
        call seconds(wall_commend)
        wall_comm = wall_comm + wall_commend - wall_commst
      endif

c     fast loop (does not need parallel execution)
c     uses results of the above loop to find Fermi level
c     and to decide on next set of energy points
      do 300 ie = n1, n2
        ipr = 1+ ie -n1
        ee = emg(ie)

        if (ie.eq.1 .and. iflrp.ne.1) then
c         the absolutely first point on energy grid
          do 206 iph = 0,nph
          do 206 il = 0,lx
  206     xrhocp(il,iph) = xrhoce(il,iph, ipr)
          do 207 iph = 0,nph
          do 207 i = 1,251
  207     yrhocp(i,iph) = yrhoce(i,iph, ipr)
        endif

        xntot = 0
        if (ie.eq.neg .and. iflrp.gt.1) iflr = 1
        fl = 0
        fr = 0
        do 210 iph = 0,nph
c         calculate density and integrated number of electrons in each
c         channel for each type of atoms density, etc., find xntot.
          call ff2g (gtr(0,iph,ipr), iph,ie, nr05(iph), xrhoce(0,0,ipr), 
     1      xrhole(0,iph,ipr), xrhocp, ee, ep, yrhole(1,0,iph,ipr),
     2      yrhoce(1,iph,ipr),yrhocp(1,iph),rhoval(1,iph),
     3      xnmues(0,iph), xnatph(iph), xntot, iflr, iflrp, fl, fr,iunf)
  210   continue

c       check whether Fermi level is found between points n1 and n2
c       and decide on next set of energy points;
        if (ie.ne.1 .or. iflrp.eq.1) xndifp = xndif
        xndif = xntot - xnferm
c       if (master) print*,'xndif = ', xndif, 'xntot = ',xntot

c       check if the fermi level is found
        if ( iflr.eq.1) then
          if (xndifp*xndif .le. 0.d0) then
c         Fermi level is found ; exit from energy loop
             if (xndif.eq.0) then
               xmunew = dble(emg(ie))
               a=0
             else
               a = xndif/(xndif-xndifp)
               do 220 i = 1,4
                 fxa = a*fl + (1-a)*fr
                 bb = dimag(dcmplx((ep-ee)*(fr+fxa)/2 +
     1                coni*dimag(ee)*(fr-fl)))
                 xndif1 = xndif + a * bb
                 a = a - xndif1 / bb
  220          continue
               xmunew = dble((1-a)*ee+a*ep)
             endif

c            add end cap corrections to the configuration and density
c            factor 2 for spin degeneracy
             do 250 iph = 0,nph
               do 230 il = 0,lx
                if (il.le.2 .or. iunf.ne.0) then
                 fl = xrhocp(il,iph) * 2
                 fr = xrhoce(il,iph,ipr) * 2
                 fxa = a*fl + (1-a)*fr
                 bb = dimag(dcmplx((ep-ee)*(fr+fxa)/2 +
     $                coni*dimag(ee)*(fr-fl)))
                 xnmues(il,iph) = xnmues(il,iph) + a * bb
                endif
  230          continue
               do 240 ir = 1,nr05(iph)
                 fl = yrhocp(ir,iph) * 2
                 fr = yrhoce(ir,iph,ipr) * 2
                 fxa = a*fl + (1-a)*fr
                 bb = dimag(dcmplx((ep-ee)*(fr+fxa)/2 +
     $                coni*dimag(ee)*(fr-fl)))
                 rhoval(ir,iph) = rhoval(ir,iph) + a * bb
  240          continue
  250        continue

c            exit from the energy loop
             goto 305
          endif
        endif
        ep = emg(ie)
        do 256 iph = 0,nph
        do 256 il = 0,lx
  256   xrhocp(il,iph) = xrhoce(il,iph, ipr)
        do 257 iph = 0,nph
        do 257 i = 1,251
  257   yrhocp(i,iph) = yrhoce(i,iph, ipr)

 300  continue

      if (n2.lt.neg .and. iflrp.gt.1) then
        n1 = n2+1
        n2 = min(neg, n2+nproc)
      else
c       set direction of search
        iflr = 1
        iflrp = 1
        idir = -1
        if (xndif.lt.0) idir = 1
        n1 = 1
        n2 = min(nproc, negx)
        do 303 ie = n1, n2
 303    emg(ie) = ep+ idir*step(iflr) * ie
      endif
      goto 25

c     END of the loop over energy in complex plane.
c     new fermi level and densities are calculated.
 305  continue

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
            if (ient.gt.1) ok = .false.
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
