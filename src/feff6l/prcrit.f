      subroutine prcrit (neout, nncrit, ik0out, cksp, fbeta, ckspc, 
     1                   fbetac, potlb0)
      implicit double precision (a-h, o-z)

c     Prepare fbeta arrays, etc., for pathfinder criteria
c
c     Note that path finder is single precision, so be sure that
c     things are correct precision in calls and declarations!
c     See declarations below for details.
c     
c     Inputs:  Reads phase.bin
c     Output:  neout   'ne', number of energy grid points
c              ik0out  index of energy grid with k=0
c              cksp    |p| at each energy grid point in single precision
c              fbeta   |f(beta)| for each angle, npot, energy point, sp
c              ckspc   |p| at each necrit point in single precision
c              fbetac  |f(beta)| for each angle, npot, nncrit point, sp
c              potlb0  unique potential labels

      include 'const.h'
      include 'dim.h'
      include 'pdata.h'

c     Output variables SINGLE PRECISION for use with path finder.
c     BE CAREFUL!!
      parameter (necrit=9, nbeta=40)
      real fbetac(-nbeta:nbeta,0:npotx,necrit), ckspc(necrit)
      real fbeta(-nbeta:nbeta,0:npotx,nex), cksp(nex)
      character*6  potlb0(0:npotx)
      character*128 messag

c     Local variables
      complex*16 cfbeta, tl
      dimension dcosb(-nbeta:nbeta)
      dimension pl(ltot+1)
      dimension iecrit(necrit)

cc      print*, 'prcrit 01 ne = ' , ne
c     Need stuff from phase.bin
c     Read phase calculation input, data returned via commons
      open (unit=1, file='phase.bin', status='old',
     1      access='sequential', form='unformatted', iostat=ios)
      call chopen (ios, 'phase.bin', 'prcrit')
      call rphbin (1)
      close (unit=1)
c     Pass out ne, ik0, potlbl (from rphbin via /pdata/)
      neout = ne
cc      print*, 'prcrit 02 ne = ' , ne
      ik0out = ik0
      do 40  i = 0, npotx
         potlb0(i) = potlbl(i)
   40 continue

c     |p| at each energy point (path finder uses invA, convert here)
      do 100  ie = 1, ne
         cksp(ie) = abs (sqrt (em(ie) - eref(ie))) / bohr
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
               fbeta(ibeta,iii,ie) = abs(cfbeta)
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
       if (nncrit .eq. 0) call fstop(' at PRCRIT: bad nncrit')
       write(messag,'(1x,a,i7)') ' nncrit in prcrit ', nncrit
       call echo(messag)

      do 320  icrit = 1, nncrit
         ie = iecrit(icrit)
         ckspc(icrit) = cksp(ie)
         do 310  ibeta = -nbeta, nbeta
            do 300  iii = 0, npot
               fbetac(ibeta,iii,icrit) = fbeta(ibeta,iii,ie)
  300       continue
  310    continue
  320 continue

      return
      end
