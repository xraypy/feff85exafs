      subroutine setlam (icalc, ie)
      implicit double precision (a-h, o-z)

c     Set lambda array based on icalc and ie
c     icalc  what to do
c      0     i0, ss exact
c      1     i1, ss exact
c      2     i2, ss exact
c     10     cute algorithm
c     <0     do exactly as told, decode as:
c               icalc = -(nmax + 100*mmax + 10 000*(iord+1))
c               Note that iord=0 <=> nmax=mmax=0, so use
c                  icalc = -10 000 for this case.
c               iord = 2*nmax + mmax, so if you want iord to control,
c               set nmax and mmax large enough-- if you want nmax and
c               mmax to control, set iord = 2*nmax + mmax...

c     inputs: ie used for cute algorithm
c             nsc used from /pdata/ to recognize ss paths
c     output: variables in /lambda/ set

      include 'const.h'
      include 'dim.h'
      include 'lambda.h'
      include 'pdata.h'
      dimension mlam0(lamtot), nlam0(lamtot)
      character messag*128
c     one degree in radians
      parameter (onedeg = .01745329252)

c     Set iord, nmax and mmax based on icalc
      if (icalc .lt. 0)  then
c        decode it and do what user wants
         icode = -icalc
         nmax = mod(icode,100)
         mmax = mod(icode,10000)/100
         iord = icode/10000 -1
      elseif (nsc .eq. 1)  then
         mmax = il0-1
         nmax = il0-1
         iord = 2*nmax + mmax
      elseif (icalc .lt. 10)  then
         iord = icalc
         mmax = iord
         nmax = iord/2
      elseif (icalc .eq. 10)  then
c        do cute algorithm
c        set mmax = L0 if straight line path, otherwise set mmax = 3
         mmax = il0-1
         do 10  ileg = 1, nleg
            mag1 = abs(beta(ileg))
            mag2 = abs(mag1 - pi)
c           if beta is not 0 or pi, path is non-linear
            if (mag1.gt.onedeg .and. mag2.gt.onedeg) mmax = 3
   10    continue
c        Set nmax based on ie and l0.
c        k <= 12 invA (ie=41)  nmax = L0
c        k >= 13 invA (ie=42)  nmax =  9
         nmax = il0-1
         if (ie .ge. 42)  nmax = 9
         iord = 2*nmax + mmax
      else
         write(messag,'(1x,a,i7)') '  icalc out of range ', icalc
         call echo(messag)
         call fstop('at SETLAM')
      endif

c-----construct index lambda (lam), (mu, nu) = mlam(lam), nlam(lam)
c     lamtot, ntot, mtot are maximum lambda, mu and nu to consider
c     Use ...0 for making indices, then sort into arrays with no
c     trailing 0 so laml0x is minimimized. (note: this is a crude
c     n**2 sort -- can 'improve' to nlog_2(n) if necessary)
      lam = 0
      do 20 in = 1, nmax+1
         n = in - 1
         do 20  im = 1, mmax+1
            m = im-1
            jord = 2*n+m
            if (jord .gt. iord)  goto 20
            if (lam .ge. lamtot)  then
               call echo('Lambda array filled, some order lost')
               goto 21
            endif
            lam = lam+1
            mlam0(lam) = -m
            nlam0(lam) = n
            if (m .eq. 0)  goto 20
            if (lam .ge. lamtot)  then
               call echo('Lambda array filled, some order lost')
               goto 21
            endif
            lam = lam+1
            mlam0(lam) = m
            nlam0(lam) = n
   20 continue
   21 continue
      lamx=lam
c     lamx must be less than lamtot
       if (lamx .gt. lamtot)
     $     call fstop(' at SETLAM lamx > lamtot')

c     laml0x is biggest lam for non-zero fmatrix, also set mmax and nmax
c     Sort mlam0 and nlam0 to use min possible laml0x
      lam = 0
      do 30  lam0 = 1, lamx
         if ((nlam0(lam0).le.l0) .and. (iabs(mlam0(lam0)).le.l0)) then
            lam = lam+1
            nlam(lam) = nlam0(lam0)
            mlam(lam) = mlam0(lam0)
            nlam0(lam0) = -1
         endif
   30 continue
      laml0x = lam
      do 40  lam0 = 1, lamx
         if (nlam0(lam0) .ge. 0)  then
            lam = lam+1
            nlam(lam) = nlam0(lam0)
            mlam(lam) = mlam0(lam0)
         endif
   40 continue

      mmaxp1 = 0
      nmax = 0
      do 50  lam = 1, lamx
         if (mlam(lam)+1 .gt. mmaxp1)  mmaxp1 = mlam(lam)+1
         if (nlam(lam) .gt. nmax)  nmax = nlam(lam)
   50 continue

      if (nmax.gt.ntot .or. mmaxp1.gt.mtot+1)  then
         call echo(' mmaxp1, nmax, mtot, ntot icalc')

         write(messag,'(3x,5i8)') 
     1        mmaxp1, nmax, mtot, ntot, icalc
         call echo(messag)

         call fstop(' at SETLAM')
      endif

      return
      end
