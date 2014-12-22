      subroutine prcrit (neout, nncrit, ik0out, cksp, fbeta, ckspc, 
     1                   fbetac, potlb0, xlam, xlamc)
      implicit double precision (a-h, o-z)

c     Prepare fbeta arrays, etc., for pathfinder criteria
c
c     Note that path finder is single precision, so be sure that
c     things are correct precision in calls and declarations!
c     See declarations below for details.
c     
c     Inputs:  Reads phase.pad
c     Output:  neout   'ne', number of energy grid points
c              ik0out  index of energy grid with k=0
c              cksp    |p| at each energy grid point in single precision
c              fbeta   |f(beta)| for each angle, npot, energy point, sp
c              ckspc   |p| at each necrit point in single precision
c              fbetac  |f(beta)| for each angle, npot, nncrit point, sp
c              potlb0  unique potential labels
c              xlam    mean free path for each energy point in Ang, sp
c              xlamc   mean free path for each nncrit point in Ang, sp

      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'
      character*6  potlbl
      dimension  potlbl(0:nphx)

c     staff originally kept in common blocks of pdata.h
      complex*16 ph, eref, em
      double precision rnrmav, xmu, edge
      dimension ph( nex, ltot+1, 0:nphx), eref(nex), em(nex),
     1 lmax(nex,0:nphx), iz(0:nphx)

c     Output variables SINGLE PRECISION for use with path finder.
c     BE CAREFUL!!
      parameter (necrit=9, nbeta=40)
      real fbetac(-nbeta:nbeta,0:nphx,necrit), ckspc(necrit)
      real fbeta(-nbeta:nbeta,0:nphx,nex), cksp(nex)
      real xlamc(necrit)
      real xlam(nex)
      character*6  potlb0(0:nphx)

c     Local variables
      complex*16 cfbeta, tl, cktmp
      dimension dcosb(-nbeta:nbeta)
      dimension pl(ltot+1)
      dimension iecrit(necrit)
      parameter (eps = 1.0e-16)
      complex*16 rkk(nex,8,nspx), eref2(nex,nspx)
      complex*16 ph4(nex, -ltot:ltot, nspx, 0:nphx)

      character*256 phpad

c     Need stuff from phase.pad
c     Read phase calculation input, data returned via commons
      phpad = 'phase.pad'
      call rdxsph (phpad, ne, ne1, ne3, npot, ihole,
     1     rnrmav, xmu, edge, ik0, ixc, rs, vint,
     2     em, eref2, iz, potlbl, ph4, rkk, lmax, lmaxp1  )
 
      do 10 ie = 1, ne
  10  eref(ie) = eref2(ie,1)
      do 20 iph = 0, npot
      do 20 ie = 1, ne
      do 20 il = 0, lmax(ie, iph)
  20  ph(ie,il+1, iph) = ph4(ie, -il, 1, iph)

      neout = ne1
      ik0out = ik0
      do 40  i = 0, nphx
         potlb0(i) = potlbl(i)
   40 continue

c     |p| at each energy point (path finder uses invA, convert here)
c     Also make mfp (xlam) in Ang
      do 100  ie = 1, ne
         cktmp = sqrt (2*(em(ie) - eref(ie)))
         cksp(ie) = real(dble (cktmp) / bohr)
c        xlam code lifted from genfmt
         xlam(ie) = 1.0e10
         if (abs(dimag(cktmp)) .gt. eps) xlam(ie) = real(1/dimag(cktmp))
         xlam(ie) = xlam(ie) * real(bohr)
  100 continue

c     Make the cos(beta)'s
c     Grid is from -40 to 40, 81 points from -1 to 1, spaced .025
      do 200  ibeta = -nbeta, nbeta
         dcosb(ibeta) = 0.025 * ibeta
  200 continue
c     watch out for round-off error
      dcosb(-nbeta) = -1
      dcosb(nbeta)  =  1

c     make fbeta (f(beta) for all energy points
      do 280  ibeta = -nbeta, nbeta
         call cpl0 (dcosb(ibeta), pl, lmaxp1)
         do 260  iii = 0, npot
            do 250  ie = 1, ne
               cfbeta = 0
               do 245  il = 1, lmax(ie,iii)+1
                  tl = (exp (2*coni*ph(ie,il,iii)) - 1) / (2*coni)
                  cfbeta = cfbeta + tl*pl(il)*(2*il-1)
  245          continue
               fbeta(ibeta,iii,ie) = real( abs(cfbeta) )
  250       continue
  260    continue
  280 continue

c     Make similar arrays for only the icrit points

c     Use 9 points at k=0,1,2,3,4,6,8,10,12 invA
c     See phmesh for energy gid definition.  These seem to work fine, 
c     and results aren't too sensitive to choices of k.  As few as 4
c     points work well (used 0,3,6,9), but time penalty for 9 points
c     is small and increased safety seems to be worth it.
      iecrit(1) = ik0
      iecrit(2) = ik0 + 5
      iecrit(3) = ik0 + 10
      iecrit(4) = ik0 + 15
      iecrit(5) = ik0 + 20
      iecrit(6) = ik0 + 30
      iecrit(7) = ik0 + 34
      iecrit(8) = ik0 + 38
      iecrit(9) = ik0 + 40

c     make sure that we have enough energy grid points to use all
c     9 iecrits
      nncrit = 0
      do 290  ie = 1, necrit
         if (iecrit(ie) .gt. ne)  goto 295
         nncrit = ie
  290 continue
  295 continue
      if (nncrit .eq. 0) call par_stop('bad nncrit in prcrit')
            

      do 320  icrit = 1, nncrit
         ie = iecrit(icrit)
         ckspc(icrit) = cksp(ie)
         xlamc(icrit) = xlam(ie)
         do 310  ibeta = -nbeta, nbeta
            do 300  iii = 0, npot
               fbetac(ibeta,iii,icrit) = fbeta(ibeta,iii,ie)
  300       continue
  310    continue
  320 continue

      return
      end
