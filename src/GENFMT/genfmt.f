      subroutine genfmt (ipr5, critcw, iorder, wnstar,
     1       ipol, ispin, le2, angks, elpty, evec, xivec, ptz)
      implicit double precision (a-h, o-z)

c     altered by matt newville (jan 1999):
c     format of feff.pad changed to packed-ascii, and all writes changed.
c     altered by alex ankudinov(feb 2000); disabled use of paths in LDOS.

      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'
      include '../HEADERS/vers.h'

c+----------------------------------------------------------------------
c     removing local common blocks, replacing them with explicit passing
c     of the various data srtuctures
c+----------------------------------------------------------------------
c     include 'clmz.h'
      complex*16 clmi(ltot+1,mtot+ntot+1,legtot)
c     include 'fmatrx.h'
      complex*16 fmati(lamtot,lamtot,legtot)
c     include 'lambda.h'
c     .   mlam(lamtot),     !mu for each lambda
c     .   nlam(lamtot),     !nu for each lambda
c     .   lamx,             !max lambda in problem
c     .   laml0x,           !max lambda for vectors involving absorbing atom
c     .   mmaxp1, nmax      !max mu in problem + 1, max nu in problem
      integer mlam(lamtot), nlam(lamtot), lamx, laml0x, mmaxp1, nmax
c     include 'nlm.h'
      dimension xnlm(ltot+1,mtot+1)
c     include 'rotmat.h'
      dimension dri(ltot+1,2*mtot+1,2*mtot+1,legtot+1)
c     include 'pdata.h'
      character*80 text(5)
      character*6  potlbl(0:nphx)
      complex*16 ph(nex,-ltot:ltot,0:nphx), eref(nex), em(nex)
      double precision rat(3,0:legtot+1)
      double precision ri(legtot), beta(legtot+1), eta(0:legtot+1)
      double precision deg, rnrmav, xmu, edge
      integer lmax(nex,0:nphx), ipot(0:legtot), iz(0:nphx), ltext(5)
      integer nsc, nleg, npot, ne, ik0, ipath, ihole
      integer kinit, linit, ilinit, lmaxp1, ntext
c     common /pdata/ ph(nex,-ltot:ltot,0:nphx), !complex phase shifts ipot=0
c     .  eref(nex),                             !complex energy reference
c     .  rat(3,0:legtot+1),                     !position of each atom, code units(bohr)
c     .  em(nex),                               !energy mesh
c     .  ri(legtot), beta(legtot+1), eta(0:legtot+1), !r, beta, eta for each leg
c     .  deg, rnrmav, xmu, edge,                !(output only)
c     .  lmax(nex,0:nphx),                      !max l with non-zero phase for each energy
c     .  ipot(0:legtot),                        !potential for each atom in path
c     .  iz(0:nphx),                            !atomic number (output only)
c     .  ltext (5),                             !length of each string
c     .  nsc, nleg,                             !nscatters, nlegs (nleg = nsc+1)
c     .  npot, ne,                              !number of potentials, energy points
c     .  ik0,                                   !index of energy grid corresponding to k=0 (edge)
c     .  ipath, ihole,                          !index of current path  and hole (output only)
c     .  kinit, linit, ilinit,                  ! initial state kappa and ang. mom.
c     .  lmaxp1,                                !largest lmax in problem + 1
c     .  ntext                                  !number of text  lines




      dimension   evec(3), xivec(3)
      complex*16  ptz
      dimension   ptz(-1:1, -1:1), lind(8)

      complex*16  rho(legtot), pmati(lamtot,lamtot,2)
      complex*16  pllp, ptrac, srho, prho, cfac
      complex*16  cchi(nex)
      complex*16  rkk(nex,8)
      complex*16  rkk2(nex,8,nspx)
      complex*16  eref2(nex,nspx), ph4, bmati
      dimension   ph4(nex,-ltot:ltot, nspx, 0:nphx)
      dimension   bmati(-mtot:mtot, 8, -mtot:mtot, 8)
      dimension   xk(nex), ckmag(nex)
      complex*16  ck(nex)
      dimension   eps1(3), eps2(3), vec1(3), vec2(3)

      character*128 string
      character*256 phpad
      character*512 slog
      logical done, wnstar

c      padlib staff
      double precision phff(nex), amff(nex),  xkr(nex)
      integer  mpadx
      parameter (mpadx = 8)
      character*75 wfmt, atsym*2
      external atsym, cwig3j, istrln

c     Input flags:
c     iorder, order of approx in f-matrix expansion (see setlam)
c             (normal use, 2.  Do ss exactly regardless of iorder)

c     used for divide-by-zero and trig tests
      parameter (eps = 1.0d-16)
      external getxk, xstar

c+---------------------------------------------------------------------
c begin intialization

      phpad = 'phase.pad'
      call genfmt_prep(phpad, ispin,
c     arguments for rdxsph
     &       ne, ne1, ne3, npot, ihole, rnrmav,
     &       xmu, edge, ik0, ixc, rs, vint,
     &       em, eref2, iz, potlbl, ph4, rkk2, lmax, lmaxp1,
c     arguments for setkap
     &       kinit, linit, ilinit,
c     argument for snlm (also a return)
     &       xnlm,
c     things set in genfmt_prep
     &       eref, ph, xk, ck, ckmag, xkr,
     &       nsp, ll, npath, ntotal, nused, xportx)

c end of intialization
c+---------------------------------------------------------------------

c     Open path input file (unit in) and read text .  Use unit 1.
      ntext  = 5
      open (unit=1, file='paths.dat', status='old', iostat=ios)
      call chopen (ios, 'paths.dat', 'genfmt')
      call rdhead (1, ntext , text, ltext)
      if (ntext  .le. 0)  then
         text (1) = ' '
      endif
c
c     Save indices of paths for use by ff2chi
      open (unit=2, file='list.dat', status='unknown', iostat=ios)
      call chopen (ios, 'list.dat', 'genfmt')
c     Put phase header on top of list.dat
      call wthead (2, ntext , text )
      write(2, 125)
 125  format (1x, 71('-'))
      write(2, 135)
 135  format('  pathindex     sig2   amp ratio    ',
     1       'deg    nlegs  r effective')

c     Open nstar.dat if necessary
      if (wnstar)  then
         open(unit=4,file='nstar.dat', status='unknown', iostat=ios)
         call chopen (ios, 'nstar.dat', 'genfmt')
         write(4,'(1x,a,f8.4)' ) ' polarization', evec
         write(4,'(1x,a)' ) ' npath     n*'
      endif

c     Set crit0 for keeping feff.dat's
      if (ipr5 .le. 0)  crit0 = 2*critcw/3
c     Make a header for the running messages.
      write(slog, 155) critcw
 155  format ('    Curved wave chi amplitude ratio', f7.2, '%')
      call wlog(slog)
      if (ipr5 .le. 0)  then
         write(slog,165) crit0
         call wlog(slog)
      endif
 165  format('    Discard feff.dat for paths with cw ratio <',
     1       f7.2, '%')
      write(slog,195)
 195  format ('    path  cw ratio     deg    nleg  reff')
      call wlog(slog)

c     open feff.pad for storing path info
c     for now, use double precision.  After it's working, try
c     single precision.
c     Use single precision for all fp numbers in feff.pad
c     !KJ replaced by next lines      open (unit=3, file='feff.pad', status='unknown', iostat=ios)
c     !KJ idem      call chopen (ios, 'feff.pad', 'genfmt')
      open (unit=3, file='feff.pad', status='unknown', iostat=ios) !KJ f1  1-06
      call chopen (ios, 'feff.pad', 'genfmt') !KJ introduced f1
c     put label line in feff.pad so other programs know it really
c     is a feff.pad file
      string = '#_feff.pad v03: ' // vfeff // vf85e
      jstr   = istrln(string)
      write(3, '(a)')  string(1:jstr)

c     save stuff that is the same for all paths
c     header, ck, central atom phase shifts
      write(3, '(a2,6(1x,i4))') '#_', npot, ne, mpadx

c     Misc stuff from phase.pad and genfmt call
 345  format(a2,3(1x,i7), 3(1x,g14.7))
      write(3, 345) '#&', ihole, iorder, ilinit, rnrmav, xmu, edge
 395  format('(',i3,'(1x,a6),',i3,'(1x,i3))')
      write(wfmt, 395) npot+1, npot+1
      write(string,wfmt) (potlbl(i),i=0,npot) , (iz(i),i=0,npot)
      jstr = istrln(string)
      write(3, '(a2,a)') '#@',string(:jstr)

c     Central atom phase shifts
      call wrpadx(3,mpadx, ph(1,ll, 0),ne)
      call wrpadx(3,mpadx, ck,ne)
      call wrpadd(3,mpadx, xkr,ne)

c     While not done, read path, find feff.
 1000 continue

c     Read current path
      call rdpath (1, done, ipol, potlbl, rat, ri, beta, eta,
     &       deg, ipot, nsc, nleg, npot, ipath)

      icalc = iorder
      if (.not.done)  then
         npath = npath + 1
         ntotal = ntotal + 1
         if (wnstar)  then
c           should be ipol=1
            do 1150 ic =1,3
               vec1(ic) = rat(ic,1) - rat(ic,0)
               vec2(ic) = rat(ic,nleg-1) - rat(ic,0)
               eps1(ic) = evec(ic)
 1150       continue
            if (elpty.ne.0.0) then
               eps2(1) = xivec(2)*evec(3)-xivec(3)*evec(2)
               eps2(2) = xivec(3)*evec(1)-xivec(1)*evec(3)
               eps2(3) = xivec(1)*evec(2)-xivec(2)*evec(1)
            endif
            ndeg = nint (deg)
            xxstar = xstar(eps1, eps2, vec1, vec2, ndeg, elpty,
     &             ilinit)
            write(4,'(1x,i6,f10.3)')  npath, xxstar
         endif

c        Need reff
         reff = 0
         do 1200  i = 1, nleg
            reff = reff + ri(i)
 1200    continue
         reff = reff/2

c        Set lambda for low k
         call setlam(icalc, 1, beta, nsc, nleg, ilinit,
     &          mlam, nlam, lamx, laml0x, mmaxp1, nmax)

c        Calculate and store rotation matrix elements
         do 1300  isc = 1, nleg
            call rot3i (lmaxp1, mmaxp1, isc, beta, dri)
 1300    continue
         if (ipol.gt.0)  then
c           one more rotation in polarization case
c           NEED MORE rot3j FOR CENTRAL ATOM ( l \pm 1 )
            call rot3i (ilinit+1, ilinit+1, nleg+1, beta, dri)
         endif

c        Start cycle over spin
         do ie = 1, ne
            cchi(ie) = 0
         enddo

         do 6000 is = 1, nsp
            if (nsp.eq.1) then
               call mmtr(bmati, ipol, ispin, le2, angks, ptz, lind,
     &                dri, eta, nsc, nleg, kinit, ilinit)
            else
               call mmtr(bmati, ipol, is, le2, angks, ptz, lind,
     &                dri, eta, nsc, nleg, kinit, ilinit)
            endif
            do 110 ie = 1, ne
               eref(ie) = eref2(ie,is)
 110        continue
            do 120 iph = 0, npot
               do 122 ie = 1, ne
                  do 124 il = -lmax(ie, iph), lmax(ie, iph)
                     ph(ie,il, iph) = ph4(ie, il, is, iph)
 124              continue
 122           continue
 120        continue
            do 130 ie = 1, ne
               do 132 kdif = 1, 8
                  rkk(ie,kdif) = rkk2(ie,kdif,is)
 132           continue
 130        continue
            do 140 ie = 1, ne
               ck(ie) = sqrt (2* (em(ie) - eref(ie)))
 140        continue

c           Big energy loop
            do 5000  ie = 1, ne
c              complex rho
               do 2010  ileg = 1, nleg
                  rho(ileg) = ck(ie) * ri(ileg)
 2010          continue
c              if ck is zero, xafs is undefined.  Make it zero and jump
c              to end of calc part of loop.
               if (abs(ck(ie)) .le. eps)  then
                  cchi(ie) = cchi(ie) + 0
                  write(slog,2055)  ie, ck(ie)
 2055             format (' genfmt: ck=0.  ie, ck(ie)',i5,1p,2e14.5)
                  call wlog(slog)
                  goto 4990
               endif
c              Calculate and store spherical wave factors c_l^(m)z^m/m!
c              in a matrix clmi(il,im,ileg), ileg=1...nleg.
c              Result is that common /clmz/ is updated for use by fmtrxi.
c
c              zero clmi arrays
               do 2100  ileg = 1, legtot
                  do 2102 im = 1, mtot+ntot+1
                     do 2104  il = 1, ltot+1
                        clmi(il,im,ileg) = 0
 2104                continue
 2102             continue
 2100          continue
               mnmxp1 = mmaxp1 + nmax
               do 2150  ileg = 1, nleg
                  isc0 = ileg-1
                  if (isc0.eq.0) isc0=nleg
                  isc1 = ileg
                  lxp1 = max (lmax(ie,ipot(isc0))+1,
     &                   lmax(ie,ipot(isc1))+1)
                  mnp1 = min (lxp1, mnmxp1)
                  call sclmz (rho, lxp1, mnp1, ileg, clmi)
 2150          continue

c              Calculate and store scattering matrices fmati.
c              First matrix
               call fmtrxi(lamx, laml0x, ie, 2, 1,
     &                clmi, mlam, nlam, xnlm, dri,
     &                ph, eta, lmax, ipot, fmati)
c              Last matrix if needed
               if (nleg .gt. 2)  then
                  call fmtrxi(laml0x, lamx, ie, nleg, nleg-1,
     &                   clmi, mlam, nlam, xnlm, dri,
     &                   ph, eta, lmax, ipot, fmati)
               endif
c              Intermediate scattering matrices
               do 2200  ilegp = 2, nsc-1
                  ileg = ilegp + 1
                  call fmtrxi(lamx, lamx, ie, ileg, ilegp,
     &                   clmi, mlam, nlam, xnlm, dri,
     &                   ph, eta, lmax, ipot, fmati)
 2200          continue

c              Big matrix multiplication loops.
c              Calculates trace of matrix product
c              M(1,N) * f(N,N-1) * ... * f(3,2) * f(2,1), as in reference.
c              We will calculate the trace over lambda_N, working from
c              right to left.
c              Use only 2 pmati arrays, alternating indp (index p)
c              1 and 2.

c              to start f(2,1) -> pmat(1)
               indp = 1
               do 2250 lmp = 1, laml0x
                  do 2252 lm = 1, lamx
                     pmati(lm,lmp,indp)= fmati(lm,lmp,1)
 2252             continue
 2250          continue

c              f(N,N-1) * ... * f(3,2) * [f(2,1)]
c              Term in [] is pmat(1)
               do 2900 isc = 2, nleg-1
c                 indp is current p matrix, indp0 is previous p matrix
                  indp = 2 - mod(isc,2)
                  indp0 = 1 + mod(indp,2)
                  do 2850  lmp = 1, laml0x
                     do 2852  lm = 1, lamx
                        pllp=0
                        do 2800 lmi = 1, lamx
                           pllp = pllp +
     1                            fmati(lm,lmi,isc)*pmati(lmi,lmp,indp0)
 2800                   continue
                        pmati(lm,lmp,indp) = pllp
 2852                continue
 2850             continue
 2900          continue

c              srho=sum pr(i), prho = prod pr(i)
               srho=0
               prho=1
               do 3200  ileg = 1, nleg
                  srho = srho + rho(ileg)
                  prho = prho * rho(ileg)
 3200          continue

c              Termination matrix, fmati(...,nleg)
c              Polarization enters only this matrix
c              this will fill fmati(...,nleg) (NO LONGER in common /fmtrxi/)
               call mmtrxi(rkk, laml0x, bmati, ie, 1, nleg,lind,
     &                clmi, mlam, nlam, xnlm, eta, fmati)

c              Final trace over matrix
               ptrac=0
               do 4400  lm = 1, laml0x
                  do 4402  lmp = 1, laml0x
                     ptrac = ptrac + fmati(lm,lmp,nleg) *
     &                      pmati(lmp,lm,indp)
 4402             continue
 4400          continue

c              Calculate xafs
c              Complex chi (without 2kr term)
c              ipot(nleg) is central atom
c              cdel1 = exp(2*coni*ph(ie,ilinit+1,0))
c              central atom phase shift are included in normalized
c              reduced matrix elements rkk(....)
               cfac = exp(coni*(srho-2*xk(ie)*reff)) / prho

c              now factor 1/(2*l0+1) is inside termination matrix
c              cchi(ie) = ptrac * cfac/(2*l0+1)
               if (nsp.eq.2 .and. is.eq.1) cfac = -cfac
               cchi(ie) = cchi(ie) + ptrac * cfac
c              write(7,5) xk(ie), -12*dimag(cchi(ie)*exp(coni*2*xk(ie)*reff))
c              5              format (3f13.5)

c              When ck(ie)=0, xafs is set to zero.  Calc above undefined.
c              Jump to here from ck(ie)=0 test above.
 4990          continue

 5000       continue
c           end of energy loop
 6000    continue
c        end of loop over spins

c        compute the importance factor of this path
         call import(ne1, nsp, ik0, reff, deg, ckmag, em, eref2,
     &          cchi, xportx, crit)

c        Write path data to feff.pad if we need it.
         if (ipr5 .ge. 1  .or.  crit .ge. crit0)  then
c           write path info
 7225       format('(i6,1x,i3,1x,f7.3,1x,f11.7,1x,e15.4,',i3,
     &             '(1x,i2))')  !KJ  1-06
c           7225         format('(i6,1x,i3,1x,f7.3,1x,f11.7,1x,f9.4,',i3,'(1x,i2))') !KJ original code
            write(wfmt, 7225) nleg
            write(string,wfmt) ipath, nleg, deg, reff*bohr,
     $             crit, (ipot(i),i=1, nleg)
            jstr = istrln(string)
            write(3,'(a2,a)') '##',string(:jstr)
            call wrpadd(3,mpadx, rat(1,1),3*nleg)
            call wrpadd(3,mpadx, beta,nleg)
            call wrpadd(3,mpadx, eta,nleg)
            call wrpadd(3,mpadx, ri,nleg)
            phffo = 0
            do 7700  ie = 1, ne
               phff(ie) = 0
               if (abs(cchi(ie)) .ge. eps) then
                  phff(ie) = atan2 (dimag(cchi(ie)),
     $                   dble(cchi(ie)))
               end if

c              remove 2 pi jumps in phase
               if (ie.gt.1) call pijump (phff(ie), phffo)
               phffo    = phff(ie)
               amff(ie) = dble(abs(cchi(ie)))
 7700       continue
            call wrpadd(3,mpadx, amff,ne)
            call wrpadd(3,mpadx, phff,ne)

c           Put feff.dat and stuff into list.dat
c           zero is debye-waller factor column
            write(2,8215) ipath, zero, crit, deg, nleg, reff*bohr
 8215       format(1x, i8, f12.5, 2f10.3, i6, f9.4) !KJ original code
c           8215         format(1x, i8, f12.5, e15.4, f10.3, i6, f9.4) !KJ  1-06

c           Tell user about the path we just did
            write(slog, 8225) ipath, crit, deg, nleg, reff*bohr
 8225       format (3x, i4, 2f10.3, i6, f9.4) !KJ original code
c           8225         format (3x, i4, e15.4, f10.3, i6, f9.4) !KJ    1-06
            call wlog(slog)
            nused = nused+1
         else
c           path unimportant, tell user
            write(slog, 8235) ipath, crit, deg, nleg, reff*bohr
 8235       format (3x, i4, 2f10.3, i6, f9.4, ' neglected')
            call wlog(slog)
         endif
c        goto next path
         goto 1000
c        done with loop over paths
      end if
c     close paths.dat, list.dat, feff.pad, nstar.dat
      close (unit=1)
      close (unit=2)
      close (unit=3)
      if (wnstar) close (unit=4)
      write(slog,'(1x,i4,a,i4,a)') nused,' paths kept, ',
     $       ntotal,' examined.'
      call wlog(slog)

      return
      end
