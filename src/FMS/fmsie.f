      subroutine fmsie(verbse, iph0, nph, lipotx, ie, em, eref, ph,
     1                 rfms, lfms, nat, iphat, rath, gtr)

c     full multiple scattering code for single energy point
c     written by a.ankudinov 06.1997 using earlier written subroutines
c     coded by b.ravel
c     modified by a.ankudinov 2001 for new matrix inversion algorithms
c     Feb. 2002, a.ankudinov: fixed logic for MPI calculations
c       lfms=0  - extended system calculataions (e.g. crystal)
c       lfms=1  - small system calculations (e.g. molecule)
c       lfms=2  - same as 1 for MPI run (forces call yprep)

      implicit double precision (a-h, o-z)
      include '../HEADERS/dim.h'

c     input
      dimension iphat(natx), rath(3,natx)
      real rat(3,natx), rfms, rdirec, toler1, toler2
      real rpart,aipart
      integer nph
c      dimension iz(0:nphx)
      complex*16 ph(lx+1, 0:nphx)

c     work space
      integer iph0
      complex*16 em, eref
      character*512 slog
      logical lcalc, verbse
      dimension lcalc(0:lx)
c     fms staff
      integer lipotx(0:nphx)
      complex gg(nspx*(lx+1)**2, nspx*(lx+1)**2, 0:nphx)
      complex gtr(0:lx, 0:nphx)
      complex xphase(nspx, -lx:lx, 0:nphx), ck(nspx)
      complex*16 dck
      complex conis
      parameter (conis = (0.0, 1.0))
      real  temper, thetax, sig2
      save

      if (rfms .le. 0.0) goto 900

c     set default (LU) inv method
      minv = 0
      rdirec = 2*rfms
      toler1 = 0.e0
      toler2 = 0.e0

      do 30 iat=1,nat
      do 30 j=1,3
   30 rat(j,iat) = real (rath(j,iat))

c     transform to single precision
      temper =0.0e0
      thetax =0.0e0
      sig2  = 0.0e0

c      it will be nice to call yprep once for all energy points,
c      fix later, and now call it every time
      if (ie.eq.1 .or. lfms.eq.0 .or. lfms.eq.2)
     1  call yprep(iph0, nat, inclus, iphat, rfms, rat)

      if (inclus.gt.1) then

cc     call fms for a cluster around central atom
       if (ie.eq.1 .and. verbse) then
          write (slog,35) inclus, iph0
  35      format ('        Doing FMS for a cluster of ',i3,
     1    ' atoms around iph = ',i2)
          call wlog (slog)
       endif

       dck=sqrt(2*(em-eref))
       rpart  = real(dble(dck))
       aipart = real(dimag(dck))
       ck(1) = cmplx(rpart, aipart)
       do 1020 ipp = 0,nph
         do 1010 ill = -lipotx(ipp), lipotx(ipp)
           rpart  = real(dble (ph( 1+abs(ill), ipp)))
           aipart = real(dimag(ph( 1+abs(ill), ipp)))
           xphase(1, ill, ipp) = cmplx(rpart, aipart)
 1010    continue
 1020  continue
       iverb=0
       if (ie.eq.1) iverb = 1
       if (.not. verbse) iverb = 0
       nsp = 1
       ispin = 0
       do 1011 ill = 0,lx
 1011  lcalc(ill) = .true.
       call fms(lfms, nsp, ispin, inclus, nph, ck, lipotx, xphase, ie,
     1  iverb, minv, rdirec, toler1, toler2, lcalc, gg)

c      make ck= i, since coni is c*16
       do 1030 ip=0,nph
         if (lfms.ne.0 .or. ip.eq.iph0) then
           do 1040 il=0,lipotx(ip)
             ix = il**2
             do 1050 im=1,2*il+1
               gtr(il, ip) = gtr(il, ip) + gg(ix+im,ix+im,ip)
 1050        continue
             gtr(il,ip)= gtr(il,ip)*
     1            exp(2*conis*xphase(1,il,ip))/(2*il+1)
 1040      continue
         endif
 1030  continue
      endif

 900  continue
      return
      end
