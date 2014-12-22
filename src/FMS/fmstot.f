      subroutine fmstot( rclust, idwopt, tk, thetad, sigma2,
     1                  lmaxph, nat, iphat, ratdbl,
     2                  ipol, ispin, le2, angks, ptz,
     3                  minv, rdirec, toler1, toler2,
     4                  elnes,ipmin,ipmax,ipstep)   !KJ added this line  1-06     
c     Calculates FMS contribution to absorption
c     uses Bruce Ravel subroutine to do FMS in self-consistency loop
c     notice that it can do FMS with polarization dependence and
c     always include l-->l-1 cahnnel.
c     written by alexei ankudinov 06.1997

      implicit double precision (a-h, o-z)
      include '../HEADERS/dim.h'
      include '../HEADERS/const.h'
      include '../HEADERS/parallel.h'
      real*8 wall_commend, wall_commst
      parameter (npadx=8)
      parameter (iblock=1)

c     input
      dimension iphat(natx),  ratdbl(3,natx)
      real rat(3,natx), rclust, rdirec, toler1, toler2
      real rpart,aipart, rnrmax, thetax, temper, sig2
      integer ne, ne1, ne3,  nph, ihole
      dimension iz(0:nphx)
      integer elnes,ipmin,ipmax,ipstep  !KJ added variables 1-06      
c     polarization data
      complex*16 ptz
      dimension ptz(-1:1,-1:1)

c     work space
      complex*16 ph(nex, -ltot:ltot, nspx, 0:nphx)
      dimension lmax(nex, 0:nphx)
c     complex energy grid emg is decomposed into em and eref to have
c     the same structure in phase.pad
      complex*16 em(nex), eref(nex, nspx)
      character*6  potlbl(0:nphx)
      character*(*) phpad
      character*512 slog
c     fms staff
      integer lmaxph(0:nphx)
      integer map(nex)
      complex gg(nspx*(lx+1)**2, nspx*(lx+1)**2, 0:nphx)
      complex gtr(ipmin:ipmax,nex),gtrloc(ipmin:ipmax,nex) !KJ I added ip index 1-06
      complex xphase(nspx, -lx:lx, 0:nphx), ck(nspx)
      complex*16 dck, bmat
      dimension kind(8), lind(8)
      logical ltrace, lcalc
      dimension lcalc(0:lx)
!KJ commented out  1-06    dimension bmat(-lx:lx,0:1,8, -lx:lx,0:1,8)
      dimension bmat(-lx:lx,0:1,8, -lx:lx,0:1,8,ipmin:ipmax) !KJ added last index
      complex*16 rkk(nex,8,nspx)
      complex*16 bmat0(-lx:lx,0:1,8, -lx:lx,0:1,8)  !KJ new variable 1-06      
      complex*16 dum(nex)
      integer  npot, nsp, ie, iverb, lfms
      integer ip,nip !KJ 1-06 added this variable - just local index
      integer i1,i2,i3,i4,i5,i6  !KJ for stupid f77 3-06
      integer jsize !KJ added for MPI communication      
      save

      call setkap (ihole,ikap,linit)

      do 10 ie = 1,nex
      do 10 ip = ipmin,ipmax !KJ I added this line
  10  gtr(ie,ip) = 0  !KJ added ip
  
      do 14 iph = 0,nphx
         do 13 ill = -lx, lx
            do 12 isp =  1,nspx
               xphase(isp, ill, iph) = 0
 12         continue
 13      continue
 14   continue

c     need less data than rphpad.f provides, also dimensions of ph
c     array are different.
      phpad = 'phase.pad'
      call rdxsph (phpad, ne, ne1, ne3, nph, ihole, rnrmav, xmu, edge,
     1     ik0, ixc, rs, vint,
     2     em, eref, iz, potlbl, ph, rkk, lmax, lmaxp1)
      call setkap (ihole, kinit, linit)
      npot = nph
      if (rclust.le.0.0) goto 100

      write(slog,15) ' Number of energy points = ', ne
  15  format(a, i3)
      call wlog(slog)

      do 20 iat=1,nat
      do 20 j=1,3
  20  rat(j,iat) = real (ratdbl(j,iat))

c     transform to single precision
      rnrmax = real(rnrmav)
      temper = real(tk)
      thetax = real(thetad)
      sig2 = real(sigma2)
      idwopx = idwopt

      iph0 = 0
      call xprep(iph0, idwopx, nat, inclus, npot, iphat, rclust,
     1 rat, iz, rnrmax, temper, thetax, sig2, minv, rdirec)

      if (inclus.gt.1) then
cc      call fms for a cluster around central atom
        write (slog,25) inclus, iph0
  25    format 
     1     (' Doing FMS for a cluster of ',i3,' atoms around iph = ',i2)
        call wlog (slog)
        call wlog (' Please, wait (updates every 20 points) ...')

        nsp = 1
        if (abs(ispin).eq.1 ) nsp = nspx
c       if (abs(ispin).eq.1) nsp = 2
c       if (nsp.gt.nspx) call par_stop(' FMS: increase nspx and rerun.')

        ltrace = .false.
!KJ 1-06  I added the do-loop around the call to bcoef for ELNES calcul.
            do i1=1,8
            do i2=0,1
            do i3=-lx,lx
            do i4=1,8
            do i5=0,1
            do i6=-lx,lx
            do ip=ipmin,ipmax
              bmat(i6,i5,i4,i3,i2,i1,ip)=dcmplx(0,0)
            enddo
            enddo
            enddo
            enddo
            enddo
            enddo
            enddo
        do ip=ipmin,ipmax,ipstep
                open(7,file='ptz.dat',form='formatted',status='unknown')
            if (elnes.eq.1) call iniptz(ptz,ip,2)  !KJ Only change ptz for ELNES !!
            call bcoef(kinit, ipol, ptz, le2, ltrace, ispin, angks, 
     1           kind, lind, bmat0)
                do i=-1,1
                write(7,'(i3,6f10.3)') ip,(ptz(i,i1),i1=-1,1)
                enddo
            do i1=1,8
            do i2=0,1
            do i3=-lx,lx
            do i4=1,8
            do i5=0,1
            do i6=-lx,lx
              bmat(i6,i5,i4,i3,i2,i1,ip)=bmat0(i6,i5,i4,i3,i2,i1)
            enddo
            enddo
            enddo
            enddo
            enddo
            enddo
!!                 bmat(:,:,:,:,:,:,ip)=bmat0(:,:,:,:,:,:)
        enddo
c !KJ end my changes to the call to bcoef        

        istart = this_process*iblock + 1
        isize = 0
        do 90 ii=istart,ne,numprocs*iblock
           do 90 ie = ii, MIN(ii+iblock-1,ne) 
            if (worker) par_type = 3
          do 30  isp = 1, nsp
            dck=sqrt(2*(em(ie)-eref(ie,isp)))
            rpart  = real( dble(dck))
            aipart = real(dimag(dck))
            ck(isp) = cmplx(rpart, aipart)
  30      continue
          do 35 ipp = 0,nph
          do 35 isp = 1, nsp
          do 35 ill = -lmaxph(ipp), lmaxph(ipp)
              rpart  = dble( ph(ie, ill, isp, ipp))
              aipart = dimag(ph(ie, ill, isp, ipp)) 
              xphase(isp, ill, ipp) = cmplx(rpart, aipart)
  35      continue
          iverb=0
          if (mod(ie,20) .eq. 0) iverb=1

          lfms = 0
          do 68  k1 = 0, lx
  68      lcalc(k1) = .false.
          do 69  k1 = 1,8
  69      if (lind(k1).ge.0.and.lind(k1).le.lx) lcalc(lind(k1)) = .true.
          call fms(lfms, nsp, ispin, inclus, npot, ck, lmaxph, xphase,
     1         ie, iverb, minv, rdirec, toler1, toler2, lcalc, gg)

          do 70  k1 = 1,8
          do 70  is1 = 1,nsp
          do 70  k2 = 1,8
          do 70  is2 = 1,nsp
            ind = ms1 + 2 * ms2
            ix1 = nsp * ( lind(k1)**2 +  lind(k1) )
            ix2 = nsp * ( lind(k2)**2 +  lind(k2) )
            ms1 = is1 - 1
            ms2 = is2 - 1
            if (lind(k2).ge.0 .and. lind(k1).ge.0) then
             do 60  m1=-lind(k1), lind(k1)
             do 60  m2=-lind(k2), lind(k2)
             
             
c !KJ I added the do-loop around the calculation of gtr for elnes
c  calculations.  1-06
    
c !KJ original call :
c !KJ              gtr(ie)=gtr(ie) + gg(ix1+nsp*m1+is1, ix2+nsp*m2+is2,0) *
c !KJ     1        bmat(m2,ms2,k2, m1,ms1,k1) * rkk(ie,k1,is1)*rkk(ie,k2,is2)
               do ip=ipmin,ipmax,ipstep
                 gtr(ip,ie) = gtr(ip,ie) +
     1                        gg(ix1+nsp*m1+is1,ix2+nsp*m2+is2,0) *
     2                        bmat(m2,ms2,k2, m1,ms1,k1,ip)
     3                      * rkk(ie,k1,is1)*rkk(ie,k2,is2)
               enddo
c !KJ end my changes
               

 60          continue

            endif
 70       continue
          
          if (worker) then
            isize = isize + 1
            do ip=ipmin,ipmax
            gtrloc(ip,isize) = gtr(ip,ie)  !KJ added :,  1-06
            enddo
          endif
 90     continue
        if (worker) par_type = 2
      endif

      nip=ipmax-ipmin+1 !KJ 1-06
      jsize=isize*nip !KJ

      if (numprocs .gt. 1) then
        call seconds(wall_commst)
c Collect gtr
        if (worker) then
c-- Send pointers for gtr buffer to master
          call par_send_int(jsize,1,0,this_process) !KJ isize -> jsize 1-06
c-- Send buffer
          if (jsize .ne. 0) 
     .      call par_send_cmplx(gtrloc,jsize,0,this_process) !Josh Kas isize -> jsize
        else
c Set up map of gtr to processor 
          do is = 1,numprocs-1
            istart = is*iblock + 1 
            do ii=istart,ne,numprocs*iblock
              do ie = ii, MIN(ii+iblock-1,ne) 
                 map(ie) = is
              enddo
            enddo
          enddo
          do i = 1,numprocs-1
c-- Receive pointers for gtr buffer from i
            call par_recv_int(jsize,1,i,i) !KJ isize -> jsize  1-06
c-- Receive buffer from i
            if (jsize .ne. 0) then
              call par_recv_cmplx(gtrloc,jsize,i,i) !Josh Kas isize -> jsize
              indx = 1
              do j = 1,ne
                if (map(j) .eq. i) then
c                 Josh Kas - changed gtrloc(ipmin:ipmax,indx)
                  do ip=ipmin,ipmax
                  gtr(ip,j) = gtrloc(ip,indx) !indx:indx+nip-1) !KJ used to be gtr(j)=gtrloc(indx) 
                  enddo
                  indx = indx + 1  !KJ replaced 1 by nip   1-06 - Josh Kas - nip back to 1
                endif
              enddo
            endif
          enddo
        endif
        call seconds(wall_commend)
        wall_comm = wall_comm + wall_commend - wall_commst
      endif

c     write fms.bin


 100  continue
      if (master) then
        open (unit=1, file='fms.bin', status='unknown', iostat=ios)
c        write down title line
         write(1,105) rclust*bohr
 105     format('FMS rfms=', f7.4)
         write(1, 110) ne, ne1, ne3,  nph, npadx, nip  !KJ added nip 1-06
 110     format(6(1x,i3))  !KJ changed 5 to 6
         do ip=ipmin,ipmax      !KJ I added the do loop
            do 120 ie = 1, ne
               aa = dble ( real ( gtr(ip,ie) ) ) !KJ I added ip index
               bb = dble (aimag ( gtr(ip,ie) ) ) !KJ ditto
               dum(ie) = dcmplx (aa,bb) !KJ ditto
 120        continue
            call wrpadx(1, npadx, dum(1), ne)
         enddo                  !KJ end of my loop
         close (unit=1)
      endif
      call par_barrier

      return
      end
