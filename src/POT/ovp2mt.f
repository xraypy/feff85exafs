      subroutine ovp2mt( nph, vtot, lrewr, qtot,ri,xnatph,lnear,
     1             inrm, imt, rnrm, rmt, cmovp, ipiv, vint, inters)
c  INPUT: nph - number of diferent potentials
c   vtot(i,iph) - potential OR density at point i for potential iph
c   lrewr       - if lrewr .gt. 0 potential will be overwritten
c                   density is never overwritten (lcoul.lt.0)
c                  lrewr=0 density calculation
c                  lrewr=1 potential calculation, vint estimated
c                  lrewr=2 potential calculation, vint is fixed
c   lcoul       -  .gt.0  (potential only) calculate charge for each iph
c                  .eq.0  (potential only) flat interstitial potential 
c                  .lt.0  (density only) calc charge inside MT spheres
c   qtot       -  for density only, total electron charge of cluster
c   ri         -  loucks radial grid
c   xnatph     -  number of atoms of type iph in the cluster
c   cmovp      -  LU decomposed overlapped matrix from movrlp.f
c   ipiv       -  pivoting indices for matrix cmovp
c  OUTPUT
c    vtot    if lrewr.gt.0  decomposed overlapped potential
c            if lrewr.le.0  old prescription for potential inside MT
c              spheres or don't want to overwrite densities
c    vint    mt zero level for potentials; charge outside mt spheres for
c            densities

      implicit double precision (a-h, o-z)

      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'

      dimension vtot(251,0:nphx), xnatph(0:nphx)
      dimension inrm(0:nphx), imt(0:nphx), rmt(0:nphx), rnrm(0:nphx)
      dimension vtotav(0:nphx)
c     work space for linear algebra
      parameter (novp=40)
      complex cmovp(novp*(nphx+1)+1,novp*(nphx+1)+1)
      complex cvovp(novp*(nphx+1)+1)
      integer ipiv(novp*(nphx+1)+1)
      dimension  ri(251)
      character*13 trans
      dimension  crho(251)
      logical lnear
      dimension lnear(0:nphx)
cpot      character*30 fname

c      get ipot and irav from inters
      ipot = mod(inters,2)
      irav = (inters-ipot) / 2
c     prepare cvovp and bvec from vtot
      ncp=0
      do 25 ip1=0,nph
      do 25 i=1,novp
        ncp = ncp + 1
        ix1 = imt(ip1)-novp + i
        cvovp(ncp)= real( vtot(ix1,ip1) )
       if (lrewr.eq.2) cvovp(ncp) = cvovp(ncp) - real(vint)
  25  continue
      do 27 ip1=0,nph
         if (irav .eq. 1) then
           rav = (rmt(ip1) + rnrm(ip1)) / 2
         elseif(irav.eq.0) then
           rav =  rnrm(ip1)
         else
           rav = ri(imt(ip1)+1)
         endif
         if (lnear(ip1)) rav = ri(imt(ip1)+1)
         call terp(ri,vtot(1,ip1),inrm(ip1)+2,3,rav,vtotav(ip1))
  27  continue
      istatx=novp*(nphx+1)+1
      trans = 'NotTransposed'
      nrhs = 1

c     find parameters for interstitial potential
      if (lrewr.gt.0) then
c        dealing with potentials
         if (lrewr.eq.1) then
c           additional equation to find vint
            ncp = ncp + 1
            cvovp(ncp) = 0
            bsum = 0
c           switch from average equation for vint to the local one
            nphlst = 0
            if (ipot .eq. 0) nphlst = nph
            do 430 iph=0,nphlst
               cvovp(ncp) = cvovp(ncp) + real(vtotav(iph)*xnatph(iph))
               bsum = bsum + xnatph(iph)
  430       continue
            cvovp(ncp) = cvovp(ncp) / real(bsum)
         endif

         call cgetrs(trans, ncp, nrhs, cmovp, istatx,
     $               ipiv, cvovp, istatx, info)
         if (info.lt.0) then
             call par_stop('    *** Error in cgetrf')
c            stop
         endif

         if (lrewr.eq.1) vint = dble(real(cvovp(ncp))) /100.0d0

c        rewrite vtot
         do 550 iph=0,nph
 
cpot  to write out ovp tot pot and it's mt approxim, comment out cpot
cpot         write(fname,172)  iph
cpot  172    format('potp', i2.2, '.dat')
cpot         open (unit=1, file=fname, status='unknown', iostat=ios)
cpot         call chopen (ios, fname, 'wpot')

            do 500 i=1,novp
              index1=imt(iph)-novp + i
              index2=i+novp*iph

cpot            write(1,176) i, ri(index1), 
cpot     1             vtot(index1,iph),  dble(real(cvovp(index2)))+vint
cpot  176       format (1x, i4, 1p, 3e12.4)

              vtot(index1,iph) = dble(real(cvovp(index2)))+vint
  500       continue

cpot         close (unit=1)

c           use second order extrapolation
            j=imt(iph)+1
            call terp (ri,vtot(1,iph),imt(iph),2,ri(j),vtot(j,iph))
            do 505 j=imt(iph)+2, 251
  505       vtot(j,iph) = vint
  550    continue
      else
c        dealing with  density calculations. vint  is the total
c        charge inside mt spheres.
c        Divided by interstitial volume in istprm

         call cgetrs(trans, ncp, nrhs, cmovp, istatx,
     $            ipiv, cvovp, istatx, info)
         if (info.lt.0) then
             call par_stop('    *** Error in cgetrf')
c            stop
         endif

         vint = 0
         do 450 iph=0,nph
            do 440 i=1,imt(iph)+2
               if (i.lt.imt(iph)-novp+1) then
                 crho(i) =  vtot(i,iph)*ri(i)**2
               elseif (i.le. imt(iph)) then
                 ix1 = novp*iph +i-imt(iph)+novp
                 crho(i) = dble(real(cvovp(ix1))) * ri(i)**2
c                crho(i) =  vtot(i,iph)*ri(i)**2
               else
                 call terp(ri,crho,imt(iph),2,ri(i), crho(i) )
               endif
  440       continue
            np = imt(iph) + 2
            cdum = 0
            dpas = 0.05d0
            call somm2 (ri,crho,dpas,cdum,rmt(iph),0,np)
            vint = vint + xnatph(iph) * cdum
  450    continue
         vint=qtot-vint
      endif

      return
      end
