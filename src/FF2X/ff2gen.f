      subroutine rdxbin (s02p, erelax, wp, edgep, s02, gamach, ne1, ik0,
     1                   emxs, omega, xkxs, xsnorm, xsec, nxsec, mbconv,
     2                   title, ntitle)

      implicit double precision (a-h, o-z)
      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'

c     header from xsect.bin
      dimension ltitle(nheadx)
      character*80 title(nheadx)
      complex*16 emxs(nex), xsec(nex)
      dimension omega(nex), xkxs(nex), xsnorm(nex)

      double precision er(nex), ei(nex), xsn(nex)
      double precision col4(nex), col5(nex)

c$$$      open (unit=8, file='xsect.bin', status='old', iostat=ios)
c$$$c     read xsect.bin
c$$$      ntitle = nheadx
c$$$      call rdhead (8, ntitle, title, ltitle)
c$$$c     read method for xsec calculation
c$$$      read(8,*)  s02p, erelax, wp, edgep, emu
c$$$      if (mbconv .gt.0 .or. s02.le.0.1) s02=s02p
c$$$c     read gamach (in eV) for use in atan at absorption edge
c$$$c     and convert to code units
c$$$      read(8,*)  gamach, ne1, ik0
c$$$      gamach = gamach / hart
c$$$c     skip label and read after it
c$$$      read(8,*)
c$$$      i = 1
c$$$  300    read(8,*,end=310)  ereal, eimag, xsnorm(i), dum1, dum2
c$$$         xsec(i) = dum1 + coni*dum2
c$$$c        xsect.bin is in eV and invA, convert to code units here
c$$$         emxs(i) = (ereal + coni*eimag) / hart
c$$$         xkxs(i) = getxk (dble(emxs(i)) - edgep)
c$$$         omega(i) = dble(emxs(i)) - edgep + emu
c$$$         nxsec = i
c$$$         i = i + 1
c$$$         if (i.le.nex) goto 300
c$$$  310 continue
c$$$      close(unit=8)

      call read_xsect(ntitle, title, s02p, erelax, wp, edgep, emu,
     1                gamach, nxsec, ne1, ik0,
     2                er, ei, xsn, col4, col5)

      if (mbconv .gt.0 .or. s02.le.0.1) s02=s02p
      gamach = gamach / hart
      do 1000 i=1,nxsec
         xsec(i) = col4(i) + coni*col5(i)
         emxs(i) = (er(i) + coni*ei(i)) / hart
         xkxs(i) = getxk (dble(emxs(i)) - edgep)
         omega(i) = dble(emxs(i)) - edgep + emu
         xsnorm(i) = xsn(i)
 1000 continue
      do 1010 i=1,ntitle
         ltitle(i) = istrln(title(i))
 1010 continue

      return
      end

      subroutine wrhead (iunit, nhead, head, dwcorr, s02,
     1  tk, thetad, sig2g, alphat, vrcorr, vicorr, critcw)

      implicit double precision (a-h, o-z)

      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'
      include '../HEADERS/vers.h'
      parameter (eps4 = 1.0d-4)
      character*80  head(nheadx)
      logical dwcorr
      character*2 coment
      parameter (coment='# ')

c     write miscellanious staff into headers
c     add feff verdion to the first line
      ll = istrln(head(1))
      if (ll .lt. 0)  then
        head(1)= 'Untitled'
        ll = istrln(head(1))
      endif
      write(iunit,310) coment, head(1)(1:), vfeff//vf85e
  310 format (a2, a45, t48, a30)

c     the rest of the title
      do 330  ihead = 2, nhead
         ll = istrln(head(ihead))
         if (ll .gt. 0)  then
            write(iunit,320) coment, head(ihead)(1:ll)
         endif
  320    format (a2, a)
  330 continue
      if (dwcorr)  then
         write(iunit,340)  coment, s02, tk, thetad, sig2g
  340    format (a2,' S02=', f5.3, '  Temp=', f7.2,'  Debye_temp=',f7.2,
     1        '  Global_sig2=', f8.5)
      else
         write(iunit,341)  coment, s02, sig2g
  341    format (a2, ' S02=', f5.3,
     1   '                                        Global_sig2=', f8.5)
      endif
      if (alphat .gt. zero)  then
         write(iunit,321)  coment, alphat
  321    format (a2, ' 1st and 3rd cumulants, alphat = ', 1pe20.4)
      endif

      if (abs(vrcorr).ge.eps4 .or. abs(vicorr).ge.eps4)  then
         write(iunit,342) coment, vrcorr*hart, vicorr*hart
      endif
  342 format (a2, ' Energy zero shift, vr, vi ', 1p, 2e14.5)

      if (critcw .gt. 0)  write(iunit,350) coment, critcw
  350 format (a2, ' Curved wave amplitude ratio filter ', f7.3, '%')
      write(iunit,360) coment
  360 format (a2, '    file         sig2 tot  cw amp ratio   deg',
     1        '  nlegs   reff  inp sig2')
c     stop writing misc. staff to files

      return
      end


      subroutine dwadd (ntotal,nptot,idwopt,ip,index,crit,critcw,sig2g,
     1  sig2u, dwcorr, rnrmav, nleg, deg, reff, iz, ipot, rat,tk,thetad,
     2  alphat, thetae, mbconv, s02, ne1, ck, achi, phchi, nkx, xk, xk0,
     3  xkp, cchi, iabs, nabs, ispec, ipr4, nhead,
     4  head, vrcorr, vicorr, nused)
      implicit double precision (a-h, o-z)

      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'
      parameter (eps4 = 1.0d-4)
      character*80  head(nheadx)
      parameter (npx=15000)
c     indices of paths to do, read from list.dat
      dimension ip(npx)
      real sig2u(npx)

      parameter (nfinex = 601)
      complex*16 cchi(nfinex), ccpath(nfinex), ccc, ckp
c     to keep Im part of cchi 11.18.97 ala
      complex*16 dw, dw1, dw3
      dimension xkp(nfinex), xk0(nfinex)

      logical dwcorr
      dimension rattmp(3,0:legtot)
      dimension iztmp(0:legtot)
      character*512 slog
      character*12 fname
      real rnrmav
      dimension iz(0:nphx)
c     central atom phase shift at l0
      complex ck(nex)
      real xk(nex)
      dimension index(npx)
      dimension nleg(npx)
      real deg(npx), reff(npx), crit(npx)
      dimension ipot(legtot,npx)
      real rat(3,legtot,npx)
      real achi(nex,npx), phchi(nex,npx)
      dimension sig2x(0:nphx, 0:nphx)
      character*2 coment
      parameter (coment='# ')

c     Keep stats on paths used to make chi
      nused = 0
      xkref = dble(ck(1)**2) - xk(1)*abs(xk(1))

c     open the files for sigrm and sigem
      if (idwopt.eq.1) then
         iem = 111
         open(unit=iem,file='s2_em.dat',status='unknown', iostat=ios)
         call chopen (ios, 's2_em.dat', 'sigem')
      elseif (idwopt.eq.2) then
         irm1 =111
         open(unit=irm1,file='s2_rm2.dat',status='unknown', iostat=ios)
         call chopen (ios, 's2_rm2.dat', 'sigrm')
         irm2 = 112
         open(unit=irm2,file='s2_rm1.dat',status='unknown', iostat=ios)
         call chopen (ios, 's2_rm1.dat', 'sigrm')
      endif
      if (alphat .gt. 0) then
        icum = 113
        open(unit=icum, file='cum.dat', status='unknown', iostat=ios)
        call chopen (ios, 'cum.dat', 'sig3')
        Write(icum, 363)
  363  format('# first and third icumulant for single scattering paths')
        write(icum,364) thetae, alphat
        write(icum,365)
  364   format ('# Einstein-Temp. =',f9.2 ,'   ', 'alpha=',f9.5)
  365   format ('#       file   sig1    sig2    sig3 ')
      endif

      if (idwopt.ge.1) then
c        initialize statistics for max DW for sigrm
         sig2mx=0
         do 400 iph1=0,nphx
         do 400 iph2=0,nphx
  400    sig2x(iph1, iph2) = 0
      endif


c     cycle over all paths in the list
      do 560  ilist = 1, ntotal
c        find index of path
         do 410  j = 1, nptot
            if (ip(ilist) .eq. index(j))  then
               ipath = j
               goto 430
            endif
  410    continue
         write(slog,420)  ilist, ip(ilist)
  420    format (' did not find path i, ip(i) ', 2i10)
         call wlog(slog)
  430    continue
c        Path ipath is the path from feff.pad that corresponds to
c        the path ilist in list.dat.  The index of the path is
c        ip(ilist) and index(ipath).

c        Use this path if it passes critcw filter
         if (crit(ipath) .lt. critcw)  goto 550

c        do debye-waller factors, get sig2d from correlated debye
c        model if required
c        A note about units:  sig2g, sig2u() and sig2d are all in
c        Angs**2.  Convert to code units after we've summed them.
         sig2 = sig2g + sig2u(ilist)
         if (dwcorr .and. idwopt.ge.0)  then
c           note that stuff from feff.pad is single precision and
c           mostly in multidim. arrays.  sigms is double precision
c           and its arrays are dimensioned for a single path, so
c           use tmp variables to call it.  tk, thetad and sig2d are
c           all dp, and therefore OK.  Also note that sigms takes
c           inputs in angstroms, except for rs which is in bohr.
            rs = rnrmav
            do 460  ileg = 1, nleg(ipath)
               iztmp(ileg) = iz(ipot(ileg,ipath))
               do 450  j = 1, 3
                  rattmp(j,ileg) = rat(j,ileg,ipath) * bohr
  450          continue
  460       continue
            iztmp(0) = iztmp(nleg(ipath))
            do 470  j = 1,3
               rattmp(j,0) = rattmp(j,nleg(ipath))
  470       continue
            if (idwopt.eq.0) then
c             use CD model
              call sigms (tk, thetad, rs, legtot, nleg(ipath),
     1                  rattmp, iztmp, sig2d)
            elseif (idwopt.eq.1) then
c             use EM method
              call sigem
     1        (sig2mx,sig2x,iem,tk,ipath,nleg(ipath),rattmp,sig2d)
            elseif (idwopt.eq.3) then  !KJ 7/06 added this section
c             use CL model
              call sigcl (tk, thetad, rs, legtot, nleg(ipath),
     1                  rattmp, iztmp, sig2d)
            else
c             use RM
              call sigrm
     1        (sig2mx,sig2x,irm1,irm2,tk,ipath,nleg(ipath),rattmp,sig2d)
            endif
            sig2 = sig2 + sig2d
         endif
         sig2 = sig2 / (bohr**2)

c        Do first and third cumulants
         sig1 = 0
         sig3 = 0
         if (alphat .gt. zero  .and. nleg(ipath) .eq. 2)  then
           if (thetae.le.0.d0) then
c            call sig3  to get sig1 and sig3 for single scattering paths
c           use reff(ipath) for r, note that reff is single precision
             iz1 = iztmp(nleg(ipath))
             iz2 = iztmp(1)
             call sigte3(iz1, iz2, sig2, alphat, thetad, reff(ipath),
     1                   sig1, sig3)
           else
c            this gets sig1 and sig3 for single scattering paths
c            using Morse potential
             call sigm3(sig1, sig2, sig3, tk, alphat, thetae)
           endif
           write(icum,475) index(ipath),  sig1 * bohr,
     1                 sig2*(bohr**2), sig3*(bohr**3)
  475      format( i10,f9.5,f9.5,' ',f9.7)
         endif

c        put the debye-waller factor and other cumulants into
c        achi and phchi
         if (mbconv .gt. 0) s02 = 1.0
         do 480  i = 1, ne1
            dw = exp(-2 * sig2 * ck(i)**2)
            dw1 = exp (2 * coni * ck(i) * sig1)
            dw3 = exp ((-4 * coni * ck(i)**3 * sig3) / 3)
            dw = dw * dw1 * dw3
            phdw = 0.0
            if (abs(dw).gt.0) phdw = atan2 (dimag(dw), dble(dw))
            achi(i,ipath) = achi(i,ipath) *
     1           real(abs(dw) * s02 * deg(ipath))
            phchi(i,ipath) = phchi(i,ipath) + real(phdw)
  480    continue
c        make sure no 2pi jumps in phase
         do 490  i = 2, ne1
c           phchi is single precision, so use tmp variables
            curr = phchi (i, ipath)
            old = phchi (i-1, ipath)
            call pijump (curr, old)
            phchi (i, ipath) = real(curr)
  490    continue

         do 500  ik = 1, nkx
            call terp1 (xk, achi(1,ipath),  ne1, xk0(ik), achi0)
            call terp1 (xk, phchi(1,ipath), ne1, xk0(ik), phchi0)
            ccpath(ik) =
     1        achi0 * exp (coni * (2 * xk0(ik) * reff(ipath) + phchi0))
c           note that this already includes s02, deg, sig2, etc.
c           sum total complex chi
            cchi(ik) = cchi(ik) + ccpath(ik)
  500    continue
         nused = nused + 1

         if (iabs.eq.nabs) then
c           Put path into chi.dat, xmu.dat as required
            if (abs(sig2u(ilist)) .gt. 0.000001)  then
              write(3,515)  coment, index(ipath), sig2*(bohr**2),
     1          crit(ipath), deg(ipath), nleg(ipath), reff(ipath)*bohr,
     2          sig2u(ilist)
              write(8,515)  coment, index(ipath), sig2*(bohr**2),
     1          crit(ipath), deg(ipath), nleg(ipath), reff(ipath)*bohr,
     2            sig2u(ilist)
            else
              write(3,515) coment, index(ipath), sig2*(bohr**2),
     1          crit(ipath), deg(ipath), nleg(ipath), reff(ipath)*bohr
              write(8,515) coment, index(ipath), sig2*(bohr**2),
     1          crit(ipath), deg(ipath), nleg(ipath), reff(ipath)*bohr
            endif
  515       format(a2, 1x, i10, 5x, f9.5, 2f10.2, i6, f9.4, f9.5)
         endif

c        write out a chinnnn.dat for this path, if necessary.
         if (ipr4 .eq. 2 .and. iabs.eq.nabs .and. ispec.eq.0)  then
c           make filename chipnnnn.dat
            write(fname,520)  index(ipath)
  520       format('chip', i4.4, '.dat')
            open (unit=9, file=fname, status='unknown',iostat=ios)
            call chopen (ios, fname, 'ff2chi')
            do 530  ihead = 1, nhead
               lhead = istrln(head(ihead))
               if (lhead .gt. 0)  then
                  write(9,320) head(ihead)(1:lhead)
  320             format (a)
               endif
  530       continue
            if (dwcorr)  then
               write(9,340)  s02, tk, thetad, sig2g
  340          format (' S02', f7.3, '  Temp', f8.2,'  Debye temp',f8.2,
     1        '  Global sig2', f9.5)
            else
               write(9,341)  s02, sig2g
  341          format (' S02', f7.3,
     1      '                                        Global sig2', f9.5)
            endif
            if (alphat .gt. zero)  then
               write(9,321)  alphat
  321          format (' 1st and 3rd cumulants, alphat = ', 1pe20.4)
            endif

            if (abs(vrcorr).ge.eps4 .or. abs(vicorr).ge.eps4)  then
               write(9,342)  vrcorr, vicorr
  342          format (' Energy zero shift, vr, vi ', 1p, 2e14.5)
            endif
            write(9,*) 'Debye-waller factor ', sig2, sig3

            write(9,610)
  610       format (1x, 71('-'))
            write(9,535)
  535       format ('       k         chi           mag          ',
     1              'phase        phase-2kr  @#')
            do 540  i = 1, nkx
               ckp = sqrt (xkp(i)*abs(xkp(i)) + xkref)
c              it would be better to use interpolation for ckp
c              fix later if complaints about chipnnn.dat files, ala
               xlam0 =  - dimag(ckp)
               ccc = ccpath(i) * exp(2 * reff(ipath) * xlam0)
               phase = 0
               if (abs(ccc) .gt. 0)  phase=atan2(dimag(ccc), dble(ccc))
               if (i .gt. 1)  call pijump (phase, phase0)
               phase0 = phase
               write(9,630)  xkp(i)/bohr, dimag(ccc), abs(ccc), phase,
     1                       phase-2*xk0(i)*reff(ipath)
  630          format (1x, f10.4, 3x, 4(1pe13.6,1x))
  540       continue
            close (unit=9)
         endif

  550    continue
  560 continue

c     close files opened for sigem and sigrem
      if (idwopt.eq.1) then
        close (unit=iem)
      elseif (idwopt.eq.2) then
        close (unit=irm1)
        close (unit=irm2)
      endif
      if (alphat .gt. 0) then
        close (unit=icum)
      endif

      return
      end
