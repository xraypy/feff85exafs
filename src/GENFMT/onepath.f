      subroutine onepath(index, nleg, deg, iorder,
     &       ipot, rat,
     &       ipol, evec, elpty, xivec,
     &       innnn, ijson, ri, beta, eta,
     &       ne1,col1,col2,col3,col4,col5,col6,col7)

      implicit double precision (a-h, o-z)

c+---------------------------------------------------------------------
c  compute a single path, generating the F matrix then returning the 
c  information contained in a feffNNNN.dat file
c
c  INPUT:
c    index:    path index                            integer
c    nleg:     number of legs in path                integer
c    deg:      path degeneracy                       double
c    iorder:   order of approximation in genfmt      integer
c    ipot:     array of unique potentials            integer(legtot)
c    rat:      cartesian coordinates of scatterers   double(3,0:legtot+1)
c    ipol:     flag to do polarization               integer
c    evec:     polarization vector                   double(3)
c    elpty:    ellipticity                           double
c    xivec:    direction of travel                   double(3)
c    innnn:    flag to write feffNNNN.dat file       integer
c    ijson:    flag to write feffNNNN.json file      integer
c
c    also requires a phase.bin file from an earlier run of xsph
c
c  OUTPUT
c    ri:       leg lengths                           double(legtot)
c    beta:     beta angles                           double(legtot)
c    eta:      eta angles                            double(legtot)
c    ne:       number of k-grid points               integer
c    col1:     k-grid                                double(nex)
c    col2:     central atom phase shifts             double(nex)
c    col3:     magnitude of F_eff                    double(nex)
c    col4:     phase of F_eff                        double(nex)
c    col5:     reduction factor                      double(nex)
c    col6:     mean free path                        double(nex)
c    col7:     real partof complex momentum          double(nex)
c+---------------------------------------------------------------------

      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'

c+---------------------------------------------------------------------
c     parameters related to the call to regenf
c
c     Input flags:
c     iorder, order of approx in f-matrix expansion (see setlam)
c             (normal use, 2.  Do ss exactly regardless of iorder)
c+---------------------------------------------------------------------
      double precision evec(3), xivec(3), spvec(3)
      complex*16 ptz(-1:1, -1:1)
      integer  iorder
c     integer  mfeff, ipr5
      logical  wnstar
      double precision angks, elpty
c     double precision critcw
      logical nnnn, json
      integer innnn, ijson

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
      complex caps(nex)
      double precision rat(3,0:legtot+1), rathea(3,legtot)
      double precision ri(legtot), beta(legtot+1), eta(0:legtot+1)
      double precision deg, rnrmav, xmu, edge
      integer lmax(nex,0:nphx), ipot(0:legtot), ipthea(legtot),
     &       iz(0:nphx)
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
c     parameters used in the code taken from GENFMT/genfmt.f
c+----------------------------------------------------------------------
      complex*16  rho(legtot), pmati(lamtot,lamtot,2)
      complex*16  pllp, ptrac, srho, prho, cfac
      complex*16  cchi(nex), rkk(nex,8), rkk2(nex,8,nspx)
      complex*16  eref2(nex,nspx), ph4(nex,-ltot:ltot, nspx, 0:nphx)
      complex*16  bmati(-mtot:mtot, 8, -mtot:mtot, 8)
      complex*16  ck(nex), lind(8)
      dimension   xk(nex), ckmag(nex)
c     ckp and ffmag are used to compute importance factor
c     complex*16  ckp
c     dimension   ffmag(nex)
      dimension   eps1(3), eps2(3), vec1(3), vec2(3)

      character*512 slog
      integer ntit
      character*80 titles(nheadx), lines(2*nheadx)
      
c+----------------------------------------------------------------------
c     parameters related to using padlib
c+----------------------------------------------------------------------
      real phff(nex), amff(nex)
      double precision xkr(nex)
      integer  mpadx
      parameter (mpadx = 8)

c+----------------------------------------------------------------------
c     parameters related to using fdtarr.f and fdthea.f
c+----------------------------------------------------------------------
      character*12 fname
      character*13 fjson
      dimension col1(nex), col2(nex), col3(nex), col4(nex), col5(nex)
      dimension col6(nex), col7(nex)
      real sxk(nex)
      complex sck(nex)


c     used for divide-by-zero and trig tests
      parameter (eps = 1.0e-16)
      external xstar

      dimension atarr(3,natx)
      do 5 i=1,natx
         atarr(1, i) = 0
         atarr(2, i) = 0
         atarr(3, i) = 0
 5    continue
c     atarr is a dummy array used to call mkptz
c     CAUTION: atom coordinates may have been changed by Feff for some
c     funny polarization or ellipticity.  need a test case of funny
c     pol/ell

      wnstar = .false.
c+----------------------------------------------------------------------
c     read genfmt.json and global.json
c keep at input: iorder, ipol, evec, elpty, xivec
c+----------------------------------------------------------------------
c      call regenf(mfeff, ipr5, critcw, iorder, wnstar,
c     &       ipol, ispin, le2, angks, elpty, evec, xivec, ptz)



c+----------------------------------------------------------------------
c     initialize everything needed for the genfmt calculation
c+----------------------------------------------------------------------
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

c+----------------------------------------------------------------------
c     pull out the central atom phase shifts
c+----------------------------------------------------------------------
      do 10 ie=1,ne
         caps(ie) = cmplx(ph(ie, ll, 0))
 10   continue

c+----------------------------------------------------------------------
c     read the input JSON file for this program: onepath.json
c     return ri, beta, eta, rat (like ri, but with 0th and (n++1)th atom)
c+----------------------------------------------------------------------
      le2      = 0
      ispin    = 0
      spvec(1) = 0
      spvec(2) = 0
      spvec(3) = 0
c     call json_read_onepath(index, iorder, ipol,
c    &       nleg, deg, rat, ipot, elpty, evec, xivec, nnnn, json)
      call json_read_geometry(nleg, rat, ipot)
      call pathgeom(nleg, nsc, ipol, rat, ipot, ri, beta, eta)
      call mkptz(ipol, elpty, evec, xivec, ispin, spvec, natx, atarr,
     &       angks, le2, ptz)

      nnnn = .false.
      if (innnn .gt. 0) nnnn=.true.
      json = .false.
      if (ijson .gt. 0) json=.true.

c+----------------------------------------------------------------------
c     fetch the standard output header lines from xsect.json
c+----------------------------------------------------------------------
      call read_titles(ntit, titles)

c+----------------------------------------------------------------------
c  this section is cut-n-pasted from genfmt
c  this is the loop over paragraphs in the paths.dat file
c  the call to rdpath is replaced by the reading of the onepath.json file (for now)
c+----------------------------------------------------------------------
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
      call setlam(iorder, 1, beta, nsc, nleg, ilinit,
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
         do 510 ie = 1, ne
            eref(ie) = eref2(ie,is)
 510     continue
         do 520 iph = 0, npot
            do 522 ie = 1, ne
               do 524 il = -lmax(ie, iph), lmax(ie, iph)
                  ph(ie,il, iph) = ph4(ie, il, is, iph)
 524           continue
 522        continue
 520     continue
         do 530 ie = 1, ne
            do 532 kdif = 1, 8
               rkk(ie,kdif) = rkk2(ie,kdif,is)
 532        continue
 530     continue
         do 540 ie = 1, ne
            ck(ie) = sqrt (2* (em(ie) - eref(ie)))
 540     continue

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


c+----------------------------------------------------------------------
c     compute the importance factor of this path
c+----------------------------------------------------------------------
c     call import(ne1, nsp, ik0, deg, ckmag, em, eref2,
c     &       cchi, xportx, crit)


c+----------------------------------------------------------------------
c     compute mag and phase arrays for F_eff, set single precision
c     arrays for xk and ck
c+----------------------------------------------------------------------
      phffo = 0
      do 15 ie = 1, ne
         phff(ie) = 0
         if (abs(cchi(ie)) .ge. eps) then
            phff(ie) = real(atan2 (dimag(cchi(ie)), dble(cchi(ie))))
         end if
         
c        remove 2 pi jumps in phase
         if (ie.gt.1) call pijump (dble(phff(ie)), phffo)
         phffo    = dble(phff(ie))
         amff(ie) = real(abs(cchi(ie)))
         sxk(ie)  = real(xk(ie))
         sck(ie)  = cmplx(ck(ie))
 15   continue


c+----------------------------------------------------------------------
c  the following get stored in feff.bin for each path:
c        ipath, nleg, deg, reff (*bohr), crit, ipot(1, nleg)
c        rat beta eta ri amff phff
c
c  instead, we'll skip straight to the chore performed in feffdt where
c  the stuff from feff.bin has been read and is written out to the form
c  of a feffNNNN.dat file
c+----------------------------------------------------------------------

c+----------------------------------------------------------------------
c        compute the columns of feffNNNN.dat
c+----------------------------------------------------------------------
      call fdtarr(ne1, real(reff), ilinit, amff, phff, caps, sxk,sck,
     &       col1, col2, col3, col4, col5, col6, col7)


      if (nnnn) then
c        Prepare output file feffnnnn.dat
         write(fname,20)  index
 20      format ('f3ff', i4.4, '.dat')
         write(slog,30)  index, fname
 30      format (i8, 5x, a)
         call wlog(slog)

c        Write feff.dat's
         open (unit=3, file=fname, status='unknown', iostat=ios)
         call chopen (ios, fname, 'onepath')


c+----------------------------------------------------------------------
c        write out the feffNNNN.dat header
c+----------------------------------------------------------------------
         do 36 il=1,legtot
            ipthea(il) = ipot(il)
            do 33 ix=1,3
               rathea(ix,il) = rat(ix,il)
 33         continue
 36      continue
         call fdthea(ntit, titles, index, iorder, nleg, real(deg),
     &          real(reff), real(rnrmav), real(edge), rathea, ipthea,
     &          iz, potlbl, nlines, lines)
         do 40 i=1, nlines
            write(3, 50)lines(i)
 40      continue
 50      format(a)

c+----------------------------------------------------------------------
c        write out the feffNNNN.dat columns
c+----------------------------------------------------------------------
         do 60 ie = 1, ne1
            write(3,70) col1(ie), col2(ie), col3(ie), col4(ie),
     &             col5(ie), col6(ie), col7(ie)

 60      continue
 70      format (1x, f6.3, 1x, 3(1pe11.4,1x),1pe10.3,1x,
     1          2(1pe11.4,1x))

c        Done with feff.dat
         close (unit=3)
      end if
c     end of conditional for writing feffNNNN.dat


c+----------------------------------------------------------------------
c     write out a JSON file with the same information as feffNNNN.dat
c+----------------------------------------------------------------------
      if (json) then
         write(fjson,80)  index
 80      format ('feff', i4.4, '.json')
         write(slog,90)  index, fjson
 90      format (i8, 5x, a)
         call wlog(slog)

         call json_nnnn(fjson, ntit, titles, rat, ipot, ri, beta, eta,
     &          index, iorder, nleg, deg, reff, rnrmav, edge,
     &          ne1, col1, col2, col3, col4, col5, col6, col7)

      end if
c     end of conditional for writing feffNNNN.json

      end

