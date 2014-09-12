      program onepath

      implicit double precision (a-h, o-z)

      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'

c+---------------------------------------------------------------------
c  declarations for regenf
      double precision evec(3), xivec(3)
      complex*16 ptz(-1:1, -1:1)
      integer  mfeff, ipr5, iorder
      logical  wnstar
      double precision critcw, angks, elpty


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
c      character*80 text(5)
      character*6  potlbl(0:nphx)
      complex*16 ph(nex,-ltot:ltot,0:nphx), eref(nex), em(nex)
      complex*16 caps(nex)
      double precision rat(3,0:legtot+1)
      double precision ri(legtot), beta(legtot+1), eta(0:legtot+1)
      double precision deg, rnrmav, xmu, edge
      integer lmax(nex,0:nphx), ipot(0:legtot), iz(0:nphx)
c      integer ltext(5), ntext
      integer nsc, nleg, npot, ne, ik0, ihole
      integer kinit, linit, ilinit, lmaxp1
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


c+----------------------------------------------------------------------
      complex*16  lind(8)

      complex*16  rho(legtot), pmati(lamtot,lamtot,2)
      complex*16  pllp, ptrac, srho, prho, cfac
      complex*16  cchi(nex)
      complex*16  rkk(nex,8)
      complex*16  rkk2(nex,8,nspx)
      complex*16  eref2(nex,nspx), ph4, bmati
      dimension   ph4(nex,-ltot:ltot, nspx, 0:nphx)
      dimension   bmati(-mtot:mtot, 8, -mtot:mtot, 8)
      dimension   xk(nex), ckmag(nex)
      complex*16  ck(nex), ckp
      dimension   ffmag(nex)
      dimension   eps1(3), eps2(3), vec1(3), vec2(3)

      character*512 slog
      integer ntit
      character*80 titles(nheadx)
      
c      padlib staff
      double precision phff(nex), amff(nex),  xkr(nex)
      integer  mpadx
      parameter (mpadx = 8)
      character*2 atsym
      external atsym, cwig3j, istrln

c     Input flags:
c     iorder, order of approx in f-matrix expansion (see setlam)
c             (normal use, 2.  Do ss exactly regardless of iorder)

c+----------------------------------------------------------------------
c  from feffdt.f
      character*12 fname
      complex*16 ccchi, cfms


c     used for divide-by-zero and trig tests
      parameter (eps = 1.0e-16)
      external getxk, xstar




c                 genfmt.json and global.json
      call regenf(mfeff, ipr5, critcw, iorder, wnstar,
     1            ipol, ispin, le2, angks, elpty, evec, xivec, ptz)



      call genfmt_prep(ispin,
c     arguments for rdxsph
     &       ne, ne1, ne3, npot, ihole, rnrmav,
     &       xmu, edge, ik0,
     &       em, eref2, iz, potlbl, ph4, rkk2, lmax, lmaxp1,
c     arguments for setkap
     &       kinit, linit, ilinit,
c     argument for snlm (an output)
     &       xnlm,
c     things set in genfmt_prep
     &       eref, ph, xk, ck, ckmag, xkr,
     &       nsp, ll, npath, ntotal, nused, xportx)

c     central atoms phase shifts
      do 10 ie=1,ne
         caps(ie) = ph(ie, ll, 0)
 10   continue

c+----------------------------------------------------------------------
c  this section is cut-n-pasted from genfmt
c  this is the loop over paragraphs in the paths.dat file
c  the call to rdpath is replaced by the reading of the onepath.json file (for now)
c+----------------------------------------------------------------------

      call json_read_onepath(ipol, index, nleg, nsc, deg, rat, ipot,
     &       ri, beta, eta)
c     this return ri, beta, eta, rat (like ri, but with 0th and (n++1)th atom
c                 ipath, deg, nleg

      call read_titles(ntit, titles)

      icalc = iorder
      npath = npath + 1
      ntotal = ntotal + 1
      if (wnstar)  then
c        should be ipol=1
         do 1150 ic =1,3
            vec1(ic) = rat(ic,1) - rat(ic,0)
            vec2(ic) = rat(ic,nleg-1) - rat(ic,0)
            eps1(ic) = evec(ic)
 1150    continue
         if (elpty.ne.0.0) then
            eps2(1) = xivec(2)*evec(3)-xivec(3)*evec(2)
            eps2(2) = xivec(3)*evec(1)-xivec(1)*evec(3)
            eps2(3) = xivec(1)*evec(2)-xivec(2)*evec(1)
         endif
         ndeg = nint (deg)
         xxstar = xstar(eps1, eps2, vec1, vec2, ndeg, elpty,
     &          ilinit)
c        write(4,'(1x,i6,f10.3)')  npath, xxstar
      endif

c     Need reff in code units
      reff = 0
      do 1200  i = 1, nleg
         reff = reff + ri(i)
 1200 continue 
      reff = reff/2

c     Set lambda for low k
      call setlam(icalc, 1, beta, nsc, nleg, ilinit,
     &       mlam, nlam, lamx, laml0x, mmaxp1, nmax)

c     Calculate and store rotation matrix elements
      do 1300  isc = 1, nleg
         call rot3i (lmaxp1, mmaxp1, isc, beta, dri)
 1300 continue
      if (ipol.gt.0)  then
c        one more rotation in polarization case
c        NEED MORE rot3j FOR CENTRAL ATOM ( l \pm 1 )
         call rot3i (ilinit+1, ilinit+1, nleg+1, beta, dri)
      endif 

c     Start cycle over spin
      do ie = 1, ne
         cchi(ie) = 0
      enddo

      do 6000 is = 1, nsp
         if (nsp.eq.1) then
            call mmtr(bmati, ipol, ispin, le2, angks, ptz, lind,
     &             dri, eta, nsc, nleg, kinit, ilinit)
         else
            call mmtr(bmati, ipol, is, le2, angks, ptz, lind,
     &             dri, eta, nsc, nleg, kinit, ilinit)
         endif
         do 110 ie = 1, ne
            eref(ie) = eref2(ie,is)
 110     continue
         do 120 iph = 0, npot
            do 122 ie = 1, ne
               do 124 il = -lmax(ie, iph), lmax(ie, iph)
                  ph(ie,il, iph) = ph4(ie, il, is, iph)
 124           continue
 122        continue
 120     continue
         do 130 ie = 1, ne
            do 132 kdif = 1, 8
               rkk(ie,kdif) = rkk2(ie,kdif,is)
 132        continue
 130     continue
         do 140 ie = 1, ne
            ck(ie) = sqrt (2* (em(ie) - eref(ie)))
 140     continue

c        Big energy loop
         do 5000  ie = 1, ne
c           complex rho
            do 2010  ileg = 1, nleg
               rho(ileg) = ck(ie) * ri(ileg)
 2010       continue
c           if ck is zero, xafs is undefined.  Make it zero and jump
c           to end of calc part of loop.
            if (abs(ck(ie)) .le. eps)  then
               cchi(ie) = cchi(ie) + 0
               write(slog,2055)  ie, ck(ie)
 2055          format (' genfmt: ck=0.  ie, ck(ie)',i5,1p,2e14.5)
               call wlog(slog)
               goto 4990
            endif
c           Calculate and store spherical wave factors c_l^(m)z^m/m!
c           in a matrix clmi(il,im,ileg), ileg=1...nleg.
c           Result is that common /clmz/ is updated for use by fmtrxi.
c           
c           zero clmi arrays
            do 2100  ileg = 1, legtot
               do 2102 im = 1, mtot+ntot+1
                  do 2104  il = 1, ltot+1
                     clmi(il,im,ileg) = 0
 2104             continue
 2102          continue
 2100       continue
            mnmxp1 = mmaxp1 + nmax
            do 2150  ileg = 1, nleg
               isc0 = ileg-1
               if (isc0.eq.0) isc0=nleg
               isc1 = ileg
               lxp1 = max (lmax(ie,ipot(isc0))+1,
     &                lmax(ie,ipot(isc1))+1)
               mnp1 = min (lxp1, mnmxp1)
               call sclmz (rho, lxp1, mnp1, ileg, clmi)
 2150       continue

c           Calculate and store scattering matrices fmati.
c           First matrix
            call fmtrxi(lamx, laml0x, ie, 2, 1,
     &             clmi, mlam, nlam, xnlm, dri,
     &             ph, eta, lmax, ipot, fmati)
c           Last matrix if needed
            if (nleg .gt. 2)  then
               call fmtrxi(laml0x, lamx, ie, nleg, nleg-1,
     &                clmi, mlam, nlam, xnlm, dri,
     &                ph, eta, lmax, ipot, fmati)
            endif
c           Intermediate scattering matrices
            do 2200  ilegp = 2, nsc-1
               ileg = ilegp + 1
               call fmtrxi(lamx, lamx, ie, ileg, ilegp,
     &                clmi, mlam, nlam, xnlm, dri,
     &                ph, eta, lmax, ipot, fmati)
 2200       continue

c           Big matrix multiplication loops.
c           Calculates trace of matrix product
c           M(1,N) * f(N,N-1) * ... * f(3,2) * f(2,1), as in reference.
c           We will calculate the trace over lambda_N, working from
c           right to left.
c           Use only 2 pmati arrays, alternating indp (index p)
c           1 and 2.

c           to start f(2,1) -> pmat(1)
            indp = 1
            do 2250 lmp = 1, laml0x
               do 2252 lm = 1, lamx
                  pmati(lm,lmp,indp)= fmati(lm,lmp,1)
 2252          continue
 2250       continue

c           f(N,N-1) * ... * f(3,2) * [f(2,1)]
c           Term in [] is pmat(1)
            do 2900 isc = 2, nleg-1
c              indp is current p matrix, indp0 is previous p matrix
               indp = 2 - mod(isc,2)
               indp0 = 1 + mod(indp,2)
               do 2850  lmp = 1, laml0x
                  do 2852  lm = 1, lamx
                     pllp=0
                     do 2800 lmi = 1, lamx
                        pllp = pllp +
     1                         fmati(lm,lmi,isc)*pmati(lmi,lmp,indp0)
 2800                continue
                     pmati(lm,lmp,indp) = pllp
 2852             continue
 2850          continue
 2900       continue

c           srho=sum pr(i), prho = prod pr(i)
            srho=0
            prho=1
            do 3200  ileg = 1, nleg
               srho = srho + rho(ileg)
               prho = prho * rho(ileg)
 3200       continue

c           Termination matrix, fmati(...,nleg)
c           Polarization enters only this matrix
c           this will fill fmati(...,nleg) (NO LONGER in common /fmtrxi/)
            call mmtrxi(rkk, laml0x, bmati, ie, 1, nleg,lind,
     &             clmi, mlam, nlam, xnlm, eta, fmati)

c           Final trace over matrix
            ptrac=0
            do 4400  lm = 1, laml0x
               do 4402  lmp = 1, laml0x
                  ptrac = ptrac + fmati(lm,lmp,nleg) *
     &                   pmati(lmp,lm,indp)
 4402          continue
 4400       continue

c           Calculate xafs
c           Complex chi (without 2kr term)
c           ipot(nleg) is central atom
c           cdel1 = exp(2*coni*ph(ie,ilinit+1,0))
c           central atom phase shift are included in normalized
c           reduced matrix elements rkk(....)
            cfac = exp(coni*(srho-2*xk(ie)*reff)) / prho

c           now factor 1/(2*l0+1) is inside termination matrix
c           cchi(ie) = ptrac * cfac/(2*l0+1)
            if (nsp.eq.2 .and. is.eq.1) cfac = -cfac
            cchi(ie) = cchi(ie) + ptrac * cfac
c           write(7,5) xk(ie), -12*dimag(cchi(ie)*exp(coni*2*xk(ie)*reff))
c           5              format (3f13.5)

c           When ck(ie)=0, xafs is set to zero.  Calc above undefined.
c           Jump to here from ck(ie)=0 test above.
 4990       continue

 5000    continue
c        end of energy loop
 6000 continue
c     end of loop over spins


c     Make importance factor, deg*(integral (|chi|*d|p|))
c     make ffmag (|chi|)
c     xport   importance factor
      do 6810  ie = 1, ne1
         if (nsp.eq.2) then
            eref(ie) = (eref2(ie,1) + eref2(ie,nsp)) /2
c           !KJ eref(ie) = (eref2(ie,1) + eref2(ie,2)) /2
            ck(ie) = sqrt (2* (em(ie) - eref(ie)))
         endif
         ckp = ck(ie)
         xlam0 = dimag(ck(ie)) - dimag(ckp)
         ffmag(ie) = abs( cchi(ie) * exp(2*reff*xlam0) )
 6810 continue

c     integrate from edge (ik0) to ne
      nemax = ne1 - ik0 + 1
      call trap (ckmag(ik0), ffmag(ik0), nemax, xport)
      xport = abs(deg*xport)
      if (xportx.le.0)  xportx = xport
      crit = 100 * xport / xportx
c     use line below to disable importance factor (e.g. for dichroism)
c     crit = crit0+1

      phffo = 0
      do 7700  ie = 1, ne
         phff(ie) = 0
         if (abs(cchi(ie)) .ge. eps) then
            phff(ie) = atan2 (dimag(cchi(ie)), dble(cchi(ie)))
         end if
         
c        remove 2 pi jumps in phase
         if (ie.gt.1) call pijump (phff(ie), phffo)
         phffo    = phff(ie)
         amff(ie) = dble(abs(cchi(ie)))
 7700 continue


c+----------------------------------------------------------------------
c  the following get stored in feff.bin:
c        ipath, nleg, deg, reff (*bohr), crit, ipot(1, nleg)
c        rat
c        beta
c        eta
c        ri
c        amff
c        phff
c+----------------------------------------------------------------------
c  instead, we'll skip straight to feffdt where the stuff from feff.bin
c  has been read and is written out to the form of a feffNNNN.dat file
c+----------------------------------------------------------------------

c     find index of path
      ip = index

c     Path i is the path from feff.bin that corresponds to
c     the path ilist in list.dat.  The index of the path is
c     iplst(ilist) and index(i).

c     Prepare output file feffnnnn.dat
      write(fname,220)  ip
 220  format ('f3ff', i4.4, '.dat')
      write(slog,230)  ip, fname
 230  format (i8, 5x, a)
      call wlog(slog)

c     Write feff.dat's
      open (unit=3, file=fname, status='unknown', iostat=ios)
      call chopen (ios, fname, 'onepath')

c     put header on feff.dat
      do 300  itext = 1, ntit
         ltxt = istrln(titles(itext))
         write(3,160)  titles(itext)(1:ltxt)
 300  continue
 160  format (1x, a)

      write(3,310) ip, iorder
 310  format (' Path', i5, '      icalc ', i7)
      write(3,170)
 170  format (1x, 71('-'))
      write(3,320)  nleg, deg, reff*bohr, rnrmav, 
     1       edge*hart
 320  format (1x, i3, f8.3, f9.4, f10.4, f11.5, 
     1       ' nleg, deg, reff, rnrmav(bohr), edge')
      write(3,330)
 330  format ('        x         y         z   pot at#')
      write(3,340)  (rat(j,nleg)*bohr,j=1,3), 
     1       ipot(nleg),
     1       iz(ipot(nleg)), potlbl(ipot(nleg))
 340  format (1x, 3f10.4, i3, i4, 1x, a6, '   absorbing atom')
      do 360  ileg = 1, nleg-1
         write(3,350)  (rat(j,ileg)*bohr,j=1,3), ipot(ileg),
     1          iz(ipot(ileg)), potlbl(ipot(ileg))
 350     format (1x, 3f10.4, i3, i4, 1x, a6)
 360  continue

      write(3,370)
 370  format    ('    k   real[2*phc]   mag[feff]  phase[feff]',
     1       ' red factor   lambda     real[p]@#')

c     Make the feff.dat stuff and write it to feff.dat
c     Also write out for inspection to fort.66
c     note that dimag takes complex*16 argument, aimag takes
c     single precision complex argument.  Stuff from feff.bin
c     is single precision, cchi is complex*16
      do 450  ie = 1, ne1
c        Consider chi in the standard XAFS form.  Use R = rtot/2.
         ccchi = amff(ie) * exp (coni*phff(ie))
         xlam = 1.0e10
         if (abs(aimag(ck(ie))) .gt. eps) xlam = 1/aimag(ck(ie))
         redfac = exp (-2 * aimag (caps(ie)))
         cdelt = 2*dble(caps(ie))
         cfms = ccchi * xk(ie) * reff**2 *
     1          exp(2*reff/xlam) / redfac
         if (abs(ccchi) .lt. eps)  then
            phfff = 0
         else
            phfff = atan2 (dimag(ccchi), dble(ccchi))
         endif
c        remove 2 pi jumps in phases
         if (ie .gt. 1)  then
            call pijump (phfff, phfffo)
            call pijump (cdelt, cdelto)
         endif
         phfffo = phfff
         cdelto = cdelt
c        write 1 k, momentum wrt fermi level k=sqrt(p**2-kf**2)
c        2 central atom phase shift (real part),
c        3 magnitude of feff,
c        4 phase of feff,
c        5 absorbing atom reduction factor,
c        6 mean free path = 1/(Im (p))
c        7 real part of local momentum p

         write(3,400)
     1          xk(ie)/bohr,
     2          cdelt + ilinit*pi,
     3          abs(cfms) * bohr,
     4          phfff - cdelt - ilinit*pi,
     5          redfac,
     6          xlam * bohr,
     7          dble(ck(ie))/bohr
 400     format (1x, f6.3, 1x, 3(1pe11.4,1x),1pe10.3,1x,
     1          2(1pe11.4,1x))

 450  continue

c     Done with feff.dat
      close (unit=3)

      end

