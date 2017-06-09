      program feff
c
c  EXAFS only version of FEFF6 
c  Modified by Matt Newville from original code by John Rehr
c  see LICENSE for copying details
c
       implicit none
       include 'vers.h'
       include 'dim.h'
       include 'const.h'

       integer ntitx
       parameter (ntitx = 16)
       character*128  title(ntitx), tmpstr*32, fname*64
       integer       ltit(ntitx)

       integer ibeta, ipotn, ik0, ipot, ios, ne, ie
       integer nlegx, nncrit, isporb, ipotnn, i, nlegxx
       integer ipr2, ipr3, ipr4,  iorder, icsig
       integer mphase, mpath, mfeff, mchi, ms, ntitle
       double precision s02, tk, thetad, sig2g, critcw
       double precision  angle, cosb, vicorr, vrcorr
c     Following passed to pathfinder, which is single precision.
c     Be careful to always declare these!
       integer necrit, nbeta
      parameter (necrit=9, nbeta=40)
      real fbetac(-nbeta:nbeta,0:npotx,necrit), ckspc(necrit)
      real fbeta(-nbeta:nbeta,0:npotx,nex), cksp(nex)
      real rmax, critpw, pcritk, pcrith
      character*6  potlbl(0:npotx)
      character*128 inpfil, lfile

      character*10 shole(0:9)
      character*8 sout(0:6)

      integer istat, il, iox, istrln
      external istrln

   10 format (1x, a)
ccc      vfeff = 'Feff 6L.02'


       call sca_init
       call echo_init
       call open_echofile('feff.run')
       call fstop_init('feff.err')
       call echo(vfeff)
       

       shole(0) = 'no hole'
       shole(1) = 'K shell'
       shole(2) = 'LI shell'
       shole(3) = 'LII shell'
       shole(4) = 'LIII shell'
       shole(5) = 'MI shell'
       shole(6) = 'MII shell'
       shole(7) = 'MIII shell'
       shole(8) = 'MIV shell'
       shole(9) = 'MV shell'

       sout(0) = 'H-L exch'
       sout(1) = 'D-H exch'
       sout(2) = 'Gd state'
       sout(3) = 'DH - HL '
       sout(4) = 'DH + HL '
       sout(5) = 'HLnoimag'
       sout(6) = 'Gs HL   '

       call get_inpfile('feff.inp',inpfil,istat)

       call rdinp (inpfil,
     $      mphase, mpath, mfeff, mchi, ms,
     1      ntitle, title, ltit,
     2      critcw, 
     1      ipr2, ipr3, ipr4,
     1      s02, tk, thetad, sig2g,
     1      nlegxx,
     1      rmax, critpw, pcritk, pcrith, nncrit,
     2      icsig, iorder, vrcorr, vicorr, isporb)

      do 20  i = 1, ntitle
         call echo(title(i)(1:ltit(i)))
   20 continue

      if (mphase .eq. 1)  then
         call echo( 'Calculating potentials and phases...')
         call potph (isporb, shole, sout)
         open (unit=1, file='potph.dat', status='old', iostat=ios)
         call chopen (ios, 'potph.dat', 'feff')
         close (unit=1, status='delete')
      endif

      if (ms.eq.1  .and.  mpath.eq.1)  then


         call echo('Preparing plane wave scattering amplitudes...')
cc         print*, 'feff-> prcrit  ne, nncrit = ' , ne, nncrit
         call prcrit (ne, nncrit, ik0, cksp, fbeta, ckspc, 
     1                fbetac, potlbl)
cc         print*, 'after prcrit  ne, nncrit = ' , ne, nncrit

c        Dump out fbetac for central atom and first pot
         if (ipr2 .ge. 3 .and. ipr2.ne.5)  then
            do 260  ipot = 0, 1
               do 250  ie = 1, nncrit
                  write(fname,200)  ie, ipot
  200             format ('fbeta', i1, 'p', i1, '.dat')
                  open (unit=1, file=fname)
                  write(1,210)  ipot, ie, ckspc(ie)
  210             format ('# ipot, ie, ckspc(ie) ', 2i5, 1pe20.6, /
     1                    '#  angle(degrees), fbeta/|p|,  fbeta')
                  do 240  ibeta = -nbeta, nbeta
                     cosb = .025 * ibeta
                     if (cosb .gt.  1)  cosb =  1
                     if (cosb .lt. -1)  cosb = -1
                     angle = acos (cosb)
                     write(1,230)  angle*raddeg, 
     1                  fbetac(ibeta,ipot,ie)/ckspc(ie),
     2                  fbetac(ibeta,ipot,ie)
  230                format (f10.4, 1p, 2e15.6)
  240             continue
                  close (unit=1)
  250          continue
  260       continue
         endif

         call echo('Searching for paths...')
         call paths (ckspc, fbetac, pcritk, pcrith, nncrit,
     1               rmax, nlegxx, ipotnn)

         call echo('Eliminating path degeneracies...')
         call pathsd (ckspc, fbetac, ne, ik0, cksp, fbeta,
     1                critpw, ipotnn, ipr2, 
     1                pcritk, pcrith, nncrit, potlbl)

         if (ipr2 .lt. 2)  then
            open (unit=1, file='geom.dat', status='old')
            call chopen (ios, 'geom.dat', 'feff')
            close (unit=1, status='delete')
         endif
      endif

      if (mfeff .eq. 1)  then
         call echo('Calculating EXAFS parameters...')
         call genfmt (ipr3, critcw, sig2g, iorder)
      endif

c      if (mchi .eq. 1)  then
c         call echo('Calculating chi...')
c         call ff2chi (ipr4, critcw, s02, tk, thetad, icsig,
c     1        vrcorr, vicorr)
c      endif

       call echo('Feff done.  Have a nice day.')
       end
