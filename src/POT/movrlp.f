      subroutine movrlp ( nph, nat, iphat, rat, iatph, xnatph,
     1                novr, iphovr, nnovr, rovr,
     2                imt, rmt, rnrm, ri, lnear,
     3                cmovp, ipiv, volint, inters)

c     Constructs overlap matrix based on geometry of overlapped
c     muffin-tin spheres. Uses LU decomposition for inversion of matrix
c     
      implicit double precision (a-h, o-z)

      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'

      dimension iphat(natx)
      dimension rat(3,natx)
      dimension iatph(0:nphx)
      dimension xnatph(0:nphx)
      dimension novr(0:nphx)
      dimension iphovr(novrx,0:nphx)
      dimension nnovr(novrx,0:nphx)
      dimension rovr(novrx,0:nphx)
      dimension imt(0:nphx)
      dimension rmt(0:nphx)
      dimension rnrm(0:nphx)
      logical lnear
      dimension lnear(0:nphx)
c     local
      character*512 slog
c     work space for linear algebra
      dimension ri(251)
      parameter (novp=40)
      complex cmovp(novp*(nphx+1)+1,novp*(nphx+1)+1)
      real bmat(nphx+1,novp*(nphx+1))
      integer ipiv(novp*(nphx+1)+1)
c#mn
       external dist, ii

      ntmp = 0
      iat0 = -999

c     get ipot and irav from inters
      ipot = mod(inters,2)
      irav = (inters-ipot) / 2
      do 20 i=1,251
  20  ri(i)=exp(-8.8d0+(i-1)*0.05d0)
      exphx=exp(0.025d0)

c     initiallly cmovp is a unit matrix up to ncp
      ncp = novp*(nph+1)+1
      do 30 i2=1,ncp
      do 30 i1=1,ncp
        cmovp(i1,i2) = 0.
        if ( i1.eq.i2 ) cmovp(i1,i2) = 1.
        if (i2.eq.ncp) cmovp(i1,i2) = 0.01
  30  continue
      do 40 i2=1,ncp-1
      do 40 i1=1,nph+1
        bmat (i1,i2) = 0.e0
  40  continue
      xn = 0.d0

      do 200 ip1=0,nph
        if (novr(ip1) .gt. 0 ) then
           nlast = novr(ip1)
        else
           iat0 = iatph(ip1)
           ntmp = 1
           nlast = nat
        endif
        if (irav.eq.1) then
          rav = (rmt(ip1) + rnrm(ip1)) / 2
        elseif (irav.eq.0) then
          rav =  rnrm(ip1)
        else
          rav=ri(imt(ip1)+1)
        endif
        if (lnear(ip1)) rav=ri(imt(ip1)+1)

        do 190 iat = 1,nlast
          if (novr(ip1) .gt. 0 ) then
             ntmp = nnovr(iat,ip1)
             ip2 = iphovr(iat,ip1)
             rnn = rovr(iat,ip1)
          else
            if (iat.eq.iat0) goto 190
            ip2 = iphat(iat)
            rnn = dist (rat(1,iat0), rat(1,iat))
          endif

c         correct for double counting volume and area
          if (rnn .lt. rmt(ip1)+rmt(ip2)) then
c            correct interstitial volume
             volint = volint + xnatph(ip1) * ntmp *
     1       (calcvl( rmt(ip1), rmt(ip2), rnn) +
     2       calcvl(rmt(ip1), rmt(ip2), rnn)) / 2.d0
          endif

c         using expression for vtot(jri) ,(jri=i1)
c         first fill  matrix bmat
          ix1 = ip1+1

          if (rav+rmt(ip2) .le. rnn) goto 100
          imin2 = ii( rnn-rav )
          if (imt(ip2)-imin2 .ge. novp-1) then
             write(slog,132) ip1
  132        format(' FOLP for POTENTIAL type ',i3,' is too big.')
             call wlog (slog)
             write(slog,'(a)') ' Reduce overlap using FOLP and rerun'
             call wlog (slog)
             call par_stop('MOVRLP-1')
          endif
          imin2=imt(ip2)-novp+1

          do 80 i2 = imin2,imt(ip2)
             r1=ri(i2)/exphx
             r2=ri(i2)*exphx
             if (i2.eq.imt(ip2)) r2=rmt(ip2)
             if (i2.eq.imt(ip2))   r1=(r1+2*ri(imt(ip2))-rmt(ip2))/2.d0
             if (i2.eq.imt(ip2)-1) r2=(r2+2*ri(imt(ip2))-rmt(ip2))/2.d0
             if (r2+rav .lt. rnn) goto 80
             if (r1+rav .lt. rnn) then
c               use linear interpolation between cases xr=0, xr=1
                xr = (rnn-rav-r1)/ (r2-r1)
                r1 = rnn-rav   
                temp =  (r2**2 - r1**2) / (4*rnn*rav) * ntmp
                ind2=i2+1
                if (i2.eq.imt(ip2))  ind2=i2-1
                xr = xr * (r2-ri(i2)) / (ri(ind2)-ri(i2))
                ix2 = ip2*novp + i2 - imin2 + 1
                bmat (ix1,ix2) = bmat (ix1,ix2) + real(temp*(1-xr))
                ix2 = ip2*novp + ind2 - imin2 + 1
                bmat (ix1,ix2)=bmat (ix1,ix2) + real(temp*xr)
             else
                temp = (r2**2 - r1**2) / (4*rnn*rav   ) * ntmp
                ix2 = ip2*novp + i2 - imin2 + 1
                bmat (ix1,ix2) = bmat (ix1,ix2) + real( temp)
             endif
  80      continue

c         using expression for vtot(i) ,(i<jri)
c         construct matrix  cmovp
 100      if (rmt(ip1)+rmt(ip2) .le. rnn) goto 190

          imin1=ii(rnn-rmt(ip2))
          imin2=ii(rnn-rmt(ip1))
          if (imt(ip1)-imin1.ge.novp-1 .or. imt(ip2)-imin2.ge.novp-1) 
     1               call par_stop('tell authors to INCREASE NOVP')
          imin1=imt(ip1)-novp+1
          imin2=imt(ip2)-novp+1

          do 180 i1 = imin1,imt(ip1)
            ri1=ri(i1)/exphx
            ri2=ri(i1)*exphx
            if (i1.eq.imt(ip1)) ri2=rmt(ip1)
            if (i1.eq.imt(ip1)) ri1=(ri1+2*ri(imt(ip1))-rmt(ip1))/2.d0
            if (i1.eq.imt(ip1)-1)
     1                         ri2=(ri2+2*ri(imt(ip1))-rmt(ip1))/2.d0
            ix1 = i1-imin1+1  + ip1*novp
            do 170 i2 = imin2,imt(ip2)
              r1=ri(i2)/exphx
              r2=ri(i2)*exphx
              if (i2.eq.imt(ip2)) r2=rmt(ip2)
              if (i2.eq.imt(ip2))   r1=(r1+2*ri(imt(ip2))-rmt(ip2))/2.d0
              if (i2.eq.imt(ip2)-1) r2=(r2+2*ri(imt(ip2))-rmt(ip2))/2.d0
              if (r2+ri2.lt.rnn) goto 170

c             calculate volume of intersection
              temp = calcvl(ri2,r2,rnn) + calcvl(r2,ri2,rnn)
              if (ri1+r2.gt.rnn)
     1          temp = temp - calcvl(ri1,r2,rnn) - calcvl(r2,ri1,rnn)
              if (ri2+r1.gt.rnn)
     1          temp = temp - calcvl(ri2,r1,rnn) - calcvl(r1,ri2,rnn)
              if (ri1+r1.gt.rnn)
     1          temp = temp + calcvl(ri1,r1,rnn) + calcvl(r1,ri1,rnn)
c             volume of intersection (temp) should be devided by volume
c             volume between spheres ri1 and ri2
              temp=temp / ( 4.d0/3.d0*pi * (ri2**3-ri1**3) ) * ntmp

              if (r1+ri2.lt.rnn) then
c               use linear interpolation between cases xr=0, xr=1
                xr = (rnn-ri(i1)-r1)/ (r2-r1)

                ind2=i2+1
                if (i2.eq.imt(ip2))  ind2=i2-1
                xr = xr * (r2-ri(i2)) / (ri(ind2)-ri(i2))
                ix2 = i2-imin2+1 + ip2*novp
                cmovp(ix1,ix2)=cmovp(ix1,ix2) 
     1                              +cmplx (real(temp*(1-xr)))
                ix2 = ind2-imin2+1 + ip2*novp
                cmovp(ix1,ix2)=cmovp(ix1,ix2) 
     1                               +cmplx (real(temp*xr))
                r1=rnn-ri2
              else
                ix1 = i1-imin1+1 + ip1*novp
                ix2 = i2-imin2+1 + ip2*novp
                cmovp(ix1,ix2)=cmovp(ix1,ix2)  +cmplx (real(temp))
              endif
 170        continue
 180      continue

 190     continue
         xn = xn + xnatph(ip1)
  200 continue

c     using matrix bmat fill in the last row of matrix cmvovp
c     this is additional equation to find Vint.
c     switch to local equation from average over all atoms
      if (ipot .eq. 0) then
         do 260 iph=0, nph
c          xn may differ from nat, if atom list have more natx atoms
c          see rdinp.f
           aa = xnatph(iph)/xn
           do 250 ix1 = 1, ncp-1
  250      cmovp(ncp,ix1) = cmovp(ncp, ix1) +
     &             cmplx(real(aa)*bmat(iph+1,ix1))
  260    continue
      else  
         iph=0
         do 270 ix1 = 1, ncp-1
  270    cmovp(ncp,ix1) = cmovp(ncp, ix1) + bmat(iph+1,ix1)
      endif

c --- invert matrices by LU decomposition
c     call cgetrf from lapack.  this performs an LU decomposition on
c     the matrix 
      istatx=novp*(nphx+1) + 1
      call cgetrf( ncp, ncp, cmovp, istatx, ipiv, info )
      if (info.ne.0) then
          call wlog('    *** Error in cgetrf when computing cmovp')
      endif

c     have to check that the last was not permuted, otherwise
c     the density calculation will be wrong
c     this is also why we put 0.01 in last column and not 1.0
      if (ipiv(ncp).ne.ncp) 
     .  call par_stop('illegal permutation in ipiv ')

      return
      end
