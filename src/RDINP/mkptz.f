      subroutine mkptz (ipol, elpty, evec, xivec, ispin, spvec,nat,rat,
     1                  angks, le2, ptz)
c     choose new right handed frame of reference with z along spvec,
c      y along (xivec cross spvec); simpler choice if one of them is 0.
c     get all vectors in new frame and
c     makes polarization tensor ptz when z is rotated along k-vector

c     input:
c     ipol = 0  random k-vector orientation in 3d; ptz(i,j)=\delta_{i,j}
c     ipol = 1 for polarizion vector eps and it's  complex conjugate epc
c        ptz(j,i) = 0.5 [(eps(-i))^* eps(-j) + (epc(-i))^* epc(-j)]
c        notice that complex conjugation and taking i-th component
c        are non commuting operations. (eps(-i))^* = (-)^i (epc(i))
c     ipol = 2 ptz(i,j)= i*\delta_{i,j}
c     elpty - ellipticity (optional for ipol=1)
c     xivec - direction of x-ray propagation
c     ispin - type of spin calculations
c        0 - spin independent
c        -1,1 - spin dependent potential
c        2 - calculations with spin-up potential
c       -2 - calculations with spin-down potential
c     spvec - direction of spin vector (along z at the output)
c     nat - number of atoms
c     rat - xyz cordinates of atoms (changed due to the rotations)

c     output:
c     angks - angle between k-vector and spin-vector (0-pi)
c     le2   - 0-only E1, 1-E1+M1, 2-E1+E2, 3-E1+E2+M1 transitions
c     ptz   - polarization tensor

      implicit double precision (a-h, o-z)

c     all input and output through common area /pol/
      include '../HEADERS/const.h'
      dimension evec(3), xivec(3), spvec(3), rat(3,nat)
      complex*16 ptz
      dimension ptz(-1:1, -1:1)

c     addittonal local stuff to create polarization tensor ptz(i,j)
      dimension e2(3)
      complex*16  e(3),eps,epc
      dimension eps(-1:1),epc(-1:1)
      character*512 slog

c     make z axis along propagation (XIVEC).
c     le2=0 - only E1 transitions; le2=1 - E1+M1; le2=2 - E1+E2 
      rr = xivec(1)**2 + xivec(2)**2 + xivec(3)**2
      if (rr.eq.0) then
        angks = 0
c       special case when xivec is not specified
        if (ipol.eq.1) then
c         need to know xivec for E2 and M1 transitions
c         leave only E1 contribution
          if (le2.ne.0) call wlog(
     1    '  Can do only E1 transitions. Specify k-vector for M1 or E2')
          le2 = 0
        else
c         for polarization average of circular dichroizm
          if (ispin.ne.0) then
c           spin-dependent case
            do 10 i = 1,3
  10        xivec(i) = spvec(i)
            rr = xivec(1)**2 + xivec(2)**2 + xivec(3)**2
          endif
        endif
      endif
            
              
      if (rr.gt.0) then
         rsp = sqrt(rr)
         rr = xivec(1)**2 + xivec(2)**2
         if ( rr.ne.0 .or. xivec(3).lt.0) then
           if (rr.eq. 0) then
             cst = - 1
             snt = 0
             csf = 1
             snf = 0
           else
c            rotation is defined by angles theta and fi
             rr = sqrt(rr)
             cst = xivec(3) / rsp
             snt = rr / rsp
             csf = xivec(1) / rr
             snf = xivec(2) / rr
           endif
c          rotate all vectors
           do 20 i = 1, nat
 20        call rotate (rat(1,i), cst, snt, csf, snf)
           call rotate (evec, cst, snt, csf, snf)
           call rotate (xivec, cst, snt, csf, snf)
           call rotate (spvec, cst, snt, csf, snf)
         endif
      endif


c     initialize ptz
      do 30 i=-1,1
      do 30 j=-1,1
 30   ptz(j,i) = 0

c     make ptz in the frame when z is along xivec, except ipol=0
      if (ipol .eq. 0) then
         do 40 i=-1,1
 40      ptz(i,i) = 1.d0 /3.d0
      elseif (ipol .eq. 2) then
         ptz( 1, 1) =  1.d0
         ptz(-1,-1) = -1.d0
      elseif (ipol .eq. 1) then
c       Normalize polarization vector
        x = sqrt (evec(1)**2 + evec(2)**2 + evec(3)**2)
        if (x .le. 0.000001) then
         call wlog(' STOP  Polarization vector of almost zero length')
         call wlog(' Correct POLARIZATION card')
         call par_stop('MKPTZ-1')
        endif
        do 50  i = 1, 3
         evec(i) = evec(i) / x
  50    continue
        x = sqrt (xivec(1)**2 + xivec(2)**2 + xivec(3)**2)
        if (x .gt. 0) then
c         run elliptical polarization code
          do 60  i = 1, 3
            xivec(i) = xivec(i) / x
  60      continue
          x = evec(1)*xivec(1)+evec(2)*xivec(2)+evec(3)*xivec(3)
          if (abs(x) .gt. 0.9) then
            call wlog(' polarization')
            write(slog,292)  (evec(i), i=1,3)
            call wlog(slog)
            call wlog(' incidence')
            write(slog,292) (xivec(i), i=1,3)
            call wlog(slog)
            call wlog(' dot product')
            write(slog,292)  x
            call wlog(slog)
  292       format (5x, 1p, 2e13.5)
            call wlog(' STOP polarization almost parallel' //
     1                ' to the incidence')
            call wlog(' Correct ELLIPTICITY and POLARIZATION cards')
            call par_stop('MKPTZ-2')
          endif
          if (x .ne. 0.0) then
c           if xivec not normal to evec then make in normal, keeping the
c           plane based on two vectors
            call wlog(' Changing polarization vector!')
            call wlog(' Incidence is not normal to polarization.')
            call wlog(' Check your input for errors. Run continues.')
            do 70  i = 1,3
              evec(i) = evec(i) - x*xivec(i)
  70        continue
            x = sqrt (evec(1)**2 + evec(2)**2 + evec(3)**2)
            do 80   i = 1, 3
               evec(i) = evec(i) / x
  80        continue
          endif
        else
c         elpty cannot be used with xivec=0
          elpty = 0.0
        endif 
     
        e2(1) = xivec(2)*evec(3)-xivec(3)*evec(2)
        e2(2) = xivec(3)*evec(1)-xivec(1)*evec(3)
        e2(3) = xivec(1)*evec(2)-xivec(2)*evec(1)
        do 90   i = 1,3
          e(i) = (evec(i)+elpty*e2(i)*coni)
  90    continue 
        eps(-1) =  (e(1)-coni*e(2))/sqrt(2.0)
        eps(0)  =   e(3)
        eps(1)  = -(e(1)+coni*e(2))/sqrt(2.0)
        do 100  i = 1,3
          e(i) = (evec(i)-elpty*e2(i)*coni)
  100   continue 
        epc(-1) =  (e(1)-coni*e(2))/sqrt(2.0)
        epc(0)  =   e(3)
        epc(1)  = -(e(1)+coni*e(2))/sqrt(2.0)
        do 110 i = -1,1
        do 110 j = -1,1
c         ptz(j,i) = (-1.0)**i * epc(i)*eps(-j) / (1+elpty**2)
c         above - true polarization tensor for given ellipticity, 
c         below - average over left and right in order to have
c         path reversal simmetry
          ptz(j,i) = ((-1.0)**i)*(epc(i)*eps(-j)+eps(i)*epc(-j))
     1               /(1+elpty**2)/2.0
  110   continue
      endif
c     end of making polarization tensor

      angks = 0


c     second rotate so that z parrallel to spin
c     note that new y-axis is normal to spin AND incidence vector
c     which simplifies further expression for rotation matrix
      rr = spvec(1)**2 + spvec(2)**2 + spvec(3)**2
      if (rr.gt.0) then
         rsp = sqrt(rr)
         rr = spvec(1)**2 + spvec(2)**2
         if ( rr.ne.0 .or. spvec(3).lt.0) then
           if (rr.eq. 0) then
             cst = - 1
             snt = 0
             csf = 1
             snf = 0
             angks = pi
           else
c            rotation is defined by angles theta and fi
             rr = sqrt(rr)
             cst = spvec(3) / rsp
             snt = rr / rsp
             csf = spvec(1) / rr
             snf = spvec(2) / rr
             angks = acos( cst)
           endif
c          rotate all vectors
           do 120 i = 1, nat
 120       call rotate (rat(1,i), cst, snt, csf, snf)
           call rotate (evec, cst, snt, csf, snf)
           call rotate (xivec, cst, snt, csf, snf)
         endif
      endif

      return
      end

      subroutine rotate (vec, cst, snt, csf, snf)
      implicit double precision (a-h, o-z)
c     rotates vector to a new coordinate system
c     Euler angles: alpha=phi, beta=theta, gamma=0
      dimension vec(3), temp (3)

      temp(1) = vec(1)*cst*csf + vec(2)*cst*snf - vec(3)*snt
      temp(2) = -vec(1)*snf + vec(2)*csf
      temp(3) = vec(1)*csf*snt + vec(2)*snt*snf + vec(3)*cst
      do 10 i = 1,3
  10  vec(i) = temp(i)

      return
      end
