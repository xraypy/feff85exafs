      subroutine feffloss(nph, iz, xnatph, rnrm, npoles, eps0,
     1       write_opcons, write_loss, write_exc, verbose,
     2       wpcorr, gamma, ampfac, delta)

      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'
      include '../HEADERS/parallel.h'

      integer epsmax
      parameter(epsmax = 700)
      
      integer iz(0:nphx), npoles
      character*2 Components(0:nphx)
      double precision NumDens(0:nphx), xnatph(0:nphx)
      double precision rnrm(0:nphx), eps0, sumrl, xNElec, gamma, csumrl
      double precision VTot
      logical write_opcons, write_loss, write_exc, verbose

      double precision wpcorr(MxPole), delta(MxPole), ampfac(MxPole)

      double precision thiseps(epsmax,3,0:nphx)
      double precision energy(epsmax), loss(epsmax)

      character*12 EpsFile
      character*2 atsym
      external atsym

      if (verbose) print*, 'nph, rnrm', nph, rnrm(0:nph)

c     Find the number density of each component.
      VTot = 0.d0
      do 10 iph = 0, nph
         VTot = VTot + xnatph(iph)*4.d0/3.d0*pi*(rnrm(iph)*bohr)**3
         NumDens(iph) = -1.d0
 10   continue

      do 20 iph = 0, nph
         if(NumDens(iph).lt.0.d0) NumDens(iph) = xnatph(iph)/VTot
c        test with numDens = 1. Should give same loss as input file for a single file.
c         if(.false.) then
c            if(iph.ne.0) then
c               NumDens(iph) = 1.d0
c            else
c               NumDens(iph) = 0.d0
c            end if
c         end if
         Components(iph) = atsym(iz(iph))
 20   continue

c     Get opcons{Element}.dat from database.
      do 34 iph = 0, nph
         do 32 i1 = 1,epsmax
            do 30 i2 = 1,3
               thiseps(i1,i2,iph) = 0.d0
 30         continue
 32      continue
         EpsFile = 'opcons' // atsym(iz(iph)) // '.dat'
         if (verbose) print*, Components(iph), NumDens(iph), EpsFile
         call epsdb(iz(iph), thiseps(1,1,iph))
         if (write_opcons)
     &          call write_eps(iz(iph), thiseps(1,1,iph), EpsFile)
 34   continue

      NComps = nph + 1
      call AddEps(NComps, NumDens, thiseps, NDataTot,
     &       energy, loss)
  
c     Open output file and write data.
      if (write_loss) then
         open(file='loss.dat', unit=60, status='replace',iostat=ios)
         call chopen (ios, 'loss.dat', 'libfeffloss')
         write(60,'(A)') '# E(eV)    Loss'
         do 40 i1 = 1, NDataTot
            write(60,*), energy(i1), loss(i1)
 40      continue
      end if
  
c     ! sumrl, xNElec, gamma, and csumrl -- a bit cryptic
      sumrl = 1.d0
      xNElec = 1.d0
      gamma = 0.01
      csumrl= xNElec/sumrl
      do 50 i1 = 1, NDataTot
         loss(i1) = loss(i1)*csumrl
 50   continue

c     getomi finds poles and weights
      call getomi(energy, loss, NDataTot, NPoles,
     $       wpcorr, ampfac, delta, eps0)

      if (write_exc) then
         open(unit=13,file='exc.dat',status='replace')
         write(13,'(A33,I4,A6)') '# Loss function represented with ',
     &          NPoles, ' poles'
         write(13,'(A23,f8.4)') '# Dielectric constant: ', eps0
         do 60 i1 = 1, NPoles
            write(13,'(4f20.10)'), wpcorr(i1), wpcorr(i1)*gamma,
     &              ampfac(i1), Delta(i1)
 60      continue
         close(13)
      end if

      return
      end
