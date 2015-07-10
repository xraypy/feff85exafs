      subroutine corval (verbse,
     1       ecv, xnvmu, eorb, norb, xnval, kappa, rgrd,
     2       nohole, nph, edens, edenvl, vtot, vvalgs,
     3       rmt, rnrm, ixc, rhoint, vint, jumprm,
     4       x0, ri, dx, xion, iunf, iz,
     5       adgc, adpc, dgc, dpc, ihole, lmaxsc)

c     Finds the core-valence separation for the cluster of atoms.
c     written by ala 10 1998

c     Input: necessary atomic data and the muffin-tin potential data
c     Output:
c          xnvmu - number of valence atoms for each channel
c          ecv   - core-valence separation energy
c     Algorithm:
c       definite valence electron - above -20 eV;
c       definite core electrons   - below -70 ev;
c       first find suspicious points in LDOS (central atom only)
c       between -20 and -70, which are written in eldos array
c       After sorting, the lowest valence state is found and
c       all core states above this energy are reassigned to valence.
c       The "ecv" should be between the lowest valence energy and
c       the highest core level. Also it should be far enough 
c       (see variable tol) from both of the above levels and V_int.
c       If fails to find "ecv" for a given core-valence separation,
c       then the highest core level is reassigned to valence and
c       attempt to find "ecv' is repeated.

      implicit double precision (a-h, o-z)
      include '../HEADERS/dim.h'
      include '../HEADERS/const.h'

c     input
      dimension dmagx(nrptx), dmag0(251)
      dimension vtot (251,0:nphx), vvalgs (251,0:nphx)
      dimension xnvmu(0:lx,0:nphx+1), rmt(0:nphx),rnrm(0:nphx)
      dimension ri(nrptx), ri05(251)
      dimension iz(0:nphx), xion(0:nphx),xnval(30,0:nphx)
      dimension norb(0:nphx), kappa(30,0:nphx), iiorb(0:lx,0:nphx)
      dimension eorb(30,0:nphx), eldos(0:lx,0:nphx)
      dimension lmaxsc(0:nphx), ival(0:lx, 0:nphx), ifound(0:lx)
c     input and output
      dimension edens(251,0:nphx), edenvl(251,0:nphx)

c     work space
      dimension dum(nrptx), vtotph(nrptx),vvalph(nrptx)
      dimension dgc(251,30,0:nphx+1), dpc(251,30,0:nphx+1)
      dimension adgc(10,30,0:nphx+1), adpc(10,30,0:nphx+1)
      dimension dgcn(nrptx,30), dpcn(nrptx,30)
      complex*16 xrhoce(0:lx)
      complex*16 xrhole(0:lx)
      complex*16 yrhoce(251)
      complex*16 yrhole(251,0:lx)
      complex*16 ph(lx+1)
c     dimension 32 = (0:lx)*(0:nphx)
      dimension en(32)
      integer ll(32), ip(32), icv(32)
      complex*16 emg, eref, eimag
      dimension xp(0:lx), xpeak(0:lx)
c     stuff from feff.f for rdinp, pathfinder and genfmt
c     Following passed to pathfinder, which is single precision.
      character*512 slog
      logical ok, verbse

      do 5 i=1,32
         en(i)  = 0.0d0
         ll(i)  = 0
         ip(i)  = 0
         icv(i) = 0
 5    continue

      if (verbse) then
         write (slog,10) 
 10      format('              Core-valence separation ')
         call wlog(slog)
      endif

c     initialize staff
      do 15 i= 1,251
        dmag0(i) = 0.d0
  15  ri05(i) = exp (-8.8d0+0.05d0*(i-1))
      do 20 iph = 0, nphx
      do 20 il = 0, lx
         eldos(il, iph) = 0
         iiorb(il, iph) = 0
         ival(il, iph) = 0
  20  continue

      tol = 5.0d0/hart
      if (vint - ecv.lt.tol) ecv = vint - tol
      elow = -70.0d0/hart
      ehigh = -20.0d0/hart
      eimag = coni*1.5d0/hart
c     make energy step about 0.5 eV
      ne = 1 + nint((ehigh-elow)*2*hart)
      de = (ehigh-elow)/(ne-1)

c     find out problematic energies for core-valence separation
      do 100 iph = 0, nph
      do 100 iorb = 1, norb(iph)
        if (eorb(iorb,iph).lt.ehigh-tol.and.eorb(iorb,iph).gt.elow) then
          lll = -kappa(iorb,iph) - 1
          if (lll.lt.0) lll = kappa(iorb,iph)
c        skip in special case for Hf,Lu,Ta; treat f-electrons as valence
c        or as core according to UNFREEZEF
          if((iz(iph).ge.71.and.iz(iph).le.73) .and. lll.eq.3) goto 100
          if(iunf.eq.0 .and. lll.eq.3) goto 100

          eldos(lll,iph) = eorb(iorb,iph)
          ival(lll,iph) = 1
          if (xnval(iorb,iph).lt. 0.1d0) ival(lll,iph)=-1
          iiorb(lll,iph) = iorb
        endif
  100 continue

      do 500  iph = 0, nph
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
  40     vtotph(i) = vtotph(i) - dble(eref)
         if (ixc.ge.5) then
           do 50 i = 1, jri1
  50       vvalph(i) = vvalph(i) - dble(eref)
         else
           do 60 i = 1, jri1
  60       vvalph(i) = vtotph(i)
         endif
         itmp = 0
         if (iph.eq.0 .and. nohole.lt.0) itmp = ihole

         xx = dimag(eimag)
         nfound = 0
         do 80 il = 0,lx
           xpeak(il) = (2*il+1.d0)/(6*xx*pi)
           xp(il) = 0
           ifound(il) = 1
           if (ival(il,iph).ne.0) ifound(il) = 0
           nfound = nfound + ifound(il)
  80     continue
         if (nfound .eq. lx+1) goto 500

c        start the search for suspicious maxima in LDOS for iph
         ie = 0
  200    ie=ie + 1
            emg = elow + de*(ie-1) + eimag
            call rholie( ri05, nr05, rgrd, x0, ri, emg,
     2           ixc, rmt(iph), rnrm(iph),
     3           vtotph, vvalph, xnval(1,iph), dgcn, dpcn, eref,
     4           adgc(1,1,iph), adpc(1,1,iph), xrhole,
     5           xrhoce, yrhole, yrhoce, ph,
     6           iz(iph), xion(iph), iunf, itmp,lmaxsc(iph))

c           find the suspicious peaks on ldos and correct the energy
            nfound = 0
            do 400 il = 0, lx
               if (ival(il,iph).ne.0 .and. ifound(il).eq.0) then
c                suspicious ldos; find the first peak in ldos that
c                contains more than 1 electron is not found yet
                 xx = dimag(xrhoce(il))
                 if ((ie.eq.ne .or. xx.lt.xp(il)) .and.
     1                xp(il).gt.xpeak(il)) then
                   ifound(il) = 1
                   eldos(il,iph) = elow + de*(ie-2)
c      print*,iph,' approx count is ',xp(il)*pi*dimag(eimag),' in l=',il
                 else
                   xp(il) = xx
                 endif
               endif
               nfound = nfound + ifound(il)
  400       continue
         if (nfound.lt.lx+1 .and. ie.lt.ne) goto 200

         if (nfound.lt.lx+1) then 
            call wlog ('WARNING: fatal error in subroutine corval. Try')
            call wlog ('  to reduce ca1 in SCF card. If does not help,')
            call wlog ('SEND bug report to AUTHORS')
            call par_stop('CORVAL-1')
         endif
  500 continue

c     arrange suspicious levels in order
      ne = 0
      do 600 iph = 0,nph
      do 600  il = 0, lx
         if (eldos(il,iph) .lt. 0) then
            ne = ne + 1
c           find in which position to put the new energy
            inew = ne
            do 580 ie = 1,ne-1
               if (en(ie).gt.eldos(il,iph) .and. inew.eq.ne) inew = ie
  580       continue
            do 590 ie = ne-1,inew, -1
               en(ie+1) = en(ie)
               icv(ie+1) = icv(ie)
               ll(ie+1) = ll(ie)
               ip(ie+1) = ip(ie)
  590       continue
            en(inew) = eldos(il,iph)
            icv(inew) = ival(il,iph)
            ll(inew) = il
            ip(inew) = iph
         endif
  600 continue

c     goto exit if there is no suspicious points
      if (ne.eq.0) goto 999
 
c     find the highest core and lowest valence energies
      ic = 0
      iv = ne + 1
      do 700 ie = 1,ne
         if (icv(ie).eq.-1) then
            ic = ie
         else
            if (ie.lt.iv) iv = ie
         endif
  700 continue

c     change assignment from core to valence, if core state above lowest
c     valence
      do 720 ie=iv+1,ic
        if (icv(ie).lt.0) then
           iph = ip(ie)
           icv(ie) = 1
           ival(ll(ie),iph) = 1
c          update occupation number
           xnvmu(ll(ie), iph) = xnvmu(ll(ie), iph) + 4*ll(ie)+2
c          update valence density
           iorb = iiorb(ll(ie),iph)

           do 710 ir = 1,251
             edenvl(ir,iph) =  edenvl(ir,iph) + 2*(ll(ie)+1)*
     1       (dgc(ir,iorb,iph)**2 + dpc(ir,iorb,iph)**2)/ri05(ir)**2
             if (ll(ie).ne.0) then
               edenvl(ir,iph) =  edenvl(ir,iph) + 2*ll(ie)*
     1         (dgc(ir,iorb-1,iph)**2+dpc(ir,iorb-1,iph)**2)/ri05(ir)**2
             endif
  710      continue
        endif
  720 continue
      ic = iv - 1

c     check if suggested ecv is between core and valence
      ok = .false.
      if (ic.gt. 0) then
        if (iv.le.ne) then
          if (ecv-en(ic).gt.tol .and. en(iv)-ecv.gt.tol) ok = .true.
        else
          if (ecv-en(ic).gt.tol) ok = .true.
        endif
      else
        if (iv.le.ne) then
          if (en(iv)-ecv.gt.tol) ok = .true.
        endif
      endif
      if (ok) goto 999

  800 ecv = vint - tol
      if (iv.le.ne) ecv = min(ecv,en(iv)-tol)
      if (ic.eq.0) goto 899
      if (ecv-en(ic).gt.tol) goto 899

c     need to reassign the last core state to valence
      ic = ic - 1
      iv = iv - 1
      icv(iv) = 1
      ival(ll(iv),ip(iv)) = 1
      xnvmu(ll(iv),ip(iv)) =  xnvmu(ll(iv),ip(iv)) + 4*ll(iv)+2
c     update valence density
      iph = ip(iv)
      iorb = iiorb(ll(iv),iph)
      do 810 ir = 1,251
        edenvl(ir,iph) =  edenvl(ir,iph)+ 2*(ll(iv)+1)*
     1  (dgc(ir,iorb,iph)**2 + dpc(ir,iorb,iph)**2)/ri05(ir)**2
        if (ll(iv).ne.0) then
          edenvl(ir,iph) =  edenvl(ir,iph)+ 2*ll(iv)*
     1    (dgc(ir,iorb-1,iph)**2+dpc(ir,iorb-1,iph)**2)/ri05(ir)**2
        endif
  810 continue
      go to 800

899   continue
c     update the core valence separation in array xnval
c     need to do that for second call of 'corval' and for ixc=5,6
      do 900  ie = iv, ne
         iph = ip(ie)
         lll = ll(ie)
         iorb = iiorb(lll,iph)
         if (xnval(iorb,iph).lt.0.1d0) then
            xnval(iorb,iph) = 2*lll+2
            if (lll.gt.0) xnval(iorb-1,iph) = 2*lll
         endif
  900 continue

999   continue
      return
      end
