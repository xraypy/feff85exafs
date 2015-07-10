      subroutine fmssz(verbse, iph0, ie, em, eref, ph, nph,
     1           rfms, lfms, nat, iphat, rath, amat, lipotx, gctr, gtr)
c     uses Bruce Ravel subroutine to do FMS in self-consistency loop
c     written by alexei ankudinov 06.1997

      implicit double precision (a-h, o-z)
      include '../HEADERS/dim.h'

c     input
      dimension iphat(natx), rath(3,natx)
      real rat(3,natx), rfms
      real rpart,aipart
      integer nph
c      dimension iz(0:nphx)

c     work space
      complex*16 ph(lx+1, 0:nphx)
c      integer iph
      complex*16 em, eref
      character*512 slog
c     fms staff
      integer lipotx(0:nphx)
      complex gg(nspx*(lx+1)**2, nspx*(lx+1)**2, 0:nphx)
      complex gtr(2,2, 3,0:lx, 0:nphx)
      real amat(-lx:lx,2,2, 3,0:lx), gctr(2,2, 3,0:lx, 0:nphx)
      complex xphase(nspx, -lx:lx, 0:nphx), ck(nspx)
      complex*16 dck
      real  rdirec, toler1, toler2
      logical lcalc, verbse
      dimension lcalc(0:lx)
      save

      if (rfms .gt. 0) then
        do 25 iat=1,nat
        do 25 j=1,3
   25   rat(j,iat) = real (rath(j,iat))

c       transform to single precision
        minv = 0
        rdirec = 2*rfms
        toler1 = 0.e0
        toler2 = 0.e0


c       it will be nice to call yprep once for all energy points,
c       fix later, and now call it every time
        if (ie.eq.1 .or. lfms.eq.0) 
     1    call yprep(iph0, nat, inclus, iphat, rfms, rat)

        if (inclus.gt.1) then
cc        call fms for a cluster around central atom
          if (ie.eq.1 .and. verbse) then
             write (slog,35) inclus, iph0
  35         format ('        Doing FMS for a cluster of ',i3,
     1       ' atoms around iph = ',i2)
             call wlog (slog)
          endif

          dck=sqrt(2*(em-eref))
          rpart  = real(dble(dck))
          aipart = real(dimag(dck))
          ck(1) = cmplx(rpart, aipart)
          do 50 ipp = 0,nph
            do 40 ill = -lipotx(ipp), lipotx(ipp)
              rpart  = real(dble( ph(abs(ill)+1,ipp)))
              aipart = real(dimag(ph(abs(ill)+1,ipp)))
              xphase(1, ill, ipp) = cmplx(rpart, aipart)
  40        continue
  50      continue
          iverb=0
          if (ie.eq.1) iverb = 1
c         neglect spin-flip processes (fix later for ispin=1)
          nsp = 1
          ispin = 0
          do 55 ill = 0, lx
  55      lcalc(ill) = .true.
          call fms(lfms, nsp, ispin, inclus, nph, ck, lipotx, xphase,ie,
     1     iverb, minv, rdirec, toler1, toler2, lcalc,gg)
        endif
      endif

      do 200 ip=0,nph

        if (lfms.ne.0 .or. ip.eq.iph0) then
          do 190 lpp =0,lipotx(ip)
             ix1 = lpp**2 
             do 170 im=1,2*lpp+1
c              now cycle over gtr dimensions
               do 100 iop = 1,3
               do 100 i2 = 1,2
               do 100 i1 = 1,2
                 if (rfms.gt.0 .and. inclus.gt.0) gtr(i1,i2,iop,lpp,ip)= 
     1             gtr(i1,i2,iop,lpp,ip) + amat(im-lpp-1,i1,i2,iop,lpp)
     2             * gg(ix1+im,ix1+im,ip)
                 gctr(i1, i2, iop,lpp,ip)= gctr(i1, i2, iop,lpp,ip)
     1             + amat(im-lpp-1,i1,i2,iop,lpp)
 100           continue
 170         continue
 190      continue
        endif
 200  continue

      return
      end
