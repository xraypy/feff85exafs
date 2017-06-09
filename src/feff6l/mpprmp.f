      subroutine mpprmp (npat, ipat, xp, yp, zp)

c     make path parameters,  xp, yp,zp for each atom for a given
c     path.

c     Input is list of atoms (npat, ipat(npat)), output are
c     x,y,z coord. of path in standard frame of reference
c     (see comments in timrep.f or here below)

      include 'dim.h'
      include 'pola.h'
      double precision  ro2, norm, zvec, xvec, yvec, ri, xp1, yp1, zp1
      dimension ipat(npatx+1), zvec(3), xvec(3), yvec(3)

      common /atoms/ rat(3,0:natx), ipot(0:natx), i1b(0:natx)

      dimension xp(npatx), yp(npatx), zp(npatx)
      dimension xp1(npatx), yp1(npatx), zp1(npatx)
      dimension ri(3,npatx)

      parameter (eps4 = 1.0E-4)

c        get the atoms in this path
c        we actually have them already via the ipat array

c     initialize staff
      do 10 j = 1, npatx
         xp(j) = 0
         yp(j) = 0
         zp(j) = 0
         xp1(j) = 0
         yp1(j) = 0
         zp1(j) = 0
   10 continue
      nleg = npat + 1
      do 20  j = 1, npat
      do 20  i = 1, 3
         ri(i,j) = rat(i,ipat(j)) - rat(i,0)
   20 continue
      do 30  j = nleg, npatx
      do 30  i = 1, 3
         ri(i,j) = 0
   30 continue
      do 40 i =1, 3
         xvec(i) = 0.0
         yvec(i) = 0.0
         zvec(i) = 0.0
   40 continue

      if (.not. pola) then
c        z-axis along first leg
         norm = ri(1,1)*ri(1,1)+ri(2,1)*ri(2,1)+ri(3,1)*ri(3,1)
         norm = sqrt(norm)
         do 140 i = 1, 3
           zvec(i) = ri(i,1)/norm
  140    continue
      else
c        z-axis in direction of polarization
         do 120 i = 1, 3
           zvec(i) = evec(i)
  120    continue
      endif

      do 160 j = 1,npat
      do 160 i = 1, 3
        zp1(j) = zp1(j) + zvec(i)*ri(i,j)
  160 continue

      num = 1
      if (.not. pola) then
c        first nonzero z-coord. is already positive
         goto 240
      endif
  200 continue
      if (abs(zp1(num)) .gt. eps4) then
         if (zp1(num) .lt. 0.0) then
c           inverse all z-coordinates and zvec, if 
c           first nonzero z-coordinate is negative 
            do 210 j = 1, 3
               zvec(j) = - zvec(j)
  210       continue
            do 220 j = 1, npat
               zp1(j) = - zp1(j)
  220       continue
         endif
         goto 240
      endif
      num = num +1
      if (num .lt. nleg) then
         goto 200
      endif
c     here first nonzero z-coordinate is positive
  240 continue

      num = 1
  300 continue
      ro2 = 0.0
      do 310 i =1, 3
         ro2 = ro2 + ri(i,num)*ri(i,num)
  310 continue
c     looking for first atom which is not on z-axis
      ro2 = ro2 - zp1(num)*zp1(num)
      ro2 = sqrt(abs(ro2))
      if (ro2 .ge. eps4) then
c     if atom not on the z-axis then
         if (elpty .eq. 0.0) then
c           if not elliptical polarization then
c           choose x-axis so that x-coord. positive and y=0.
            do 320 i = 1, 3
               xvec(i) = ri(i,num) - zvec(i)*zp1(num)
  320       continue
            do 330 i = 1, 3
               xvec(i) = xvec(i)/ro2
  330       continue
         else
c           if elliptical polarization then
c           choose x-axis along incident beam
            do 350 i =1, 3
               xvec(i) = ivec(i)
  350       continue
         endif
         yvec(1) = zvec(2)*xvec(3) - zvec(3)*xvec(2)
         yvec(2) = zvec(3)*xvec(1) - zvec(1)*xvec(3)
         yvec(3) = zvec(1)*xvec(2) - zvec(2)*xvec(1)
         goto 390
      endif
      num = num + 1
      if (num .lt. nleg) then
         goto 300
      endif
  390 continue

c     calculate x,y coord for each atom in chosen frame of reference
      do 400 j = 1, npat
      do 400 i =1,3
         xp1(j) = xp1(j) + xvec(i)*ri(i,j)
         yp1(j) = yp1(j) + yvec(i)*ri(i,j)
  400 continue

      if ( elpty .ne. 0.0) then
c        if no polarization or linear polarization then first nonzero
c        x-coordinate is already positive, no need to check it.
         num = 1
  500    continue
         if (abs(xp1(num)) .ge. eps4) then
            if (xp1(num) .lt. 0.0) then
               do 510 j = 1, npat
                  xp1(j) = - xp1(j)
  510          continue
            endif
            goto 520
         endif
         num = num + 1
         if (num .lt. nleg) then
            goto 500
         endif
  520    continue
      endif

      num = 1
  570 continue
c     inverse all y-coordinates if first nonzero y-coord is negative
      if (abs(yp1(num)) .ge. eps4) then
         if (yp1(num) .lt. 0.0) then
            do 580 j = 1, npat
               yp1(j) = - yp1(j)
  580       continue
         endif
         goto 590
      endif
      num = num + 1
      if (num .lt. nleg) then
         goto 570
      endif
  590 continue

      do 595 j = 1, npat
        xp(j) = xp1(j)
        yp(j) = yp1(j)
        zp(j) = zp1(j)
  595 continue
c     now xp,yp,zp represent the path in standard order
      return
      end
