      subroutine fpf0 ( iz, iholep, srho, dr, hx,
     1     dgc0, dpc0, dgc, dpc, 
     2     eatom, xnel, norb, eorb, kappa)
c      everything is input. output is written in fpf0.dat
c      to be read by ff2afs.f to get scattering amplitude

      use json_module
      implicit double precision (a-h,o-z)
      integer  iunit
      type(json_value),pointer :: fp
      double precision,dimension(81) :: qval, xirfval

      include '../HEADERS/dim.h'
      include '../HEADERS/const.h'
      include '../HEADERS/parallel.h'
c     save central atom dirac components, see comments below.
      dimension dgc0(251), dpc0(251)
      dimension dgc(251, 30, 0:nphx), dpc(251, 30, 0:nphx)
      dimension xnel(30), eorb(30), kappa(30)
      dimension srho(251), dr(251), xpc(251), xqc(251)
      logical open_16
 
c     output arrays
      dimension enosc(13), oscstr(13), index(13)
 
      if (master) then
        open (unit=16, file='fpf0.dat', status='unknown', iostat=ios)
        fpcorr =  -(iz/82.5)**2.37
        write (16,*)  ' atom Z = ', iz
        write (16,10)  eatom *alphfs**2 *5/3, fpcorr,
     1        ' total energy part of fprime - 5/3*E_tot/mc**2'
  10    format (2(1pe19.5), a)
        open_16 = .true.
        
        call json_value_create(fp)
        call to_object(fp,'fpf0.json')
        call json_value_add(fp, 'iz', iz)
        call json_value_add(fp, 'comment',
     1       'total energy part of fprime - 5/3*E_tot/mc**2')
        call json_value_add(fp, 'etot', eatom *alphfs**2 *5/3)
        call json_value_add(fp, 'fpcorr', fpcorr)
      else
        open_16 = .false.
      endif

c     get oscillator strengths
      do 20 i=1,13
        oscstr(i)=0.d0
        enosc(i)=0.d0
  20  continue
      enosc(1)= eorb(iholep)
      index(1)= iholep
      kinit = kappa(iholep)
      oscstr(1) = 2*abs(kinit)
c     always will use first spot to represent initial state
      nosc=1
      np = 251
      xmult1 = 0.0d0
      xmult2 = 0.0d0
      do 30 iorb =1, norb
        if (xnel(iorb) .gt.0.d0) then
c         it is core orbital, check if it satisfies dipole selection
          jkap = kappa(iorb)
          if (jkap+kinit.eq.0 .or. abs(jkap-kinit).eq.1) then
             nosc = nosc+1
c            calculate reduced dipole matrix element
             kdif= jkap-kinit
             if (abs(kdif).gt.1) kdif=0
c            xirf = <i |p| f> relativistic version of dipole m.e.
c            from Grant,Advan.Phys.,v.19,747(1970) eq. 6.30, using
c            Messiah's "Q.M." appendices to reduce 9j,3j symbols
c            to simple coefficients xmult1,2. ala 12.12.95
             twoj = 2.0d0*abs(kinit) - 1.0d0
             if (kdif.eq.-1 .and. kinit.gt.0) then
                xmult1 = 0.0d0
                xmult2 = sqrt(2.0d0 * (twoj+1)*(twoj-1)/twoj )
             elseif (kdif.eq.-1 .and. kinit.lt.0) then
                xmult1 = 0.0d0
                xmult2 = - sqrt(2.0d0 * (twoj+1)*(twoj+3)/(twoj+2) )
             elseif (kdif.eq. 0 .and. kinit.gt.0) then
                xmult1 = - sqrt( (twoj+1)*twoj/(twoj+2) )
                xmult2 = - sqrt( (twoj+1)*(twoj+2)/twoj )
             elseif (kdif.eq. 0 .and. kinit.lt.0) then
                xmult1 = sqrt( (twoj+1)*(twoj+2)/twoj )
                xmult2 = sqrt( (twoj+1)*twoj/(twoj+2) )
             elseif (kdif.eq. 1 .and. kinit.gt.0) then
                xmult1 = sqrt(2.0d0 * (twoj+1)*(twoj+3)/(twoj+2) )
                xmult2 = 0.0d0
             elseif (kdif.eq. 1 .and. kinit.lt.0) then
                xmult1 = - sqrt(2.0d0 * (twoj+1)*(twoj-1)/twoj )
                xmult2 = 0.0d0
             endif
             xk0 = abs(eorb(iorb)-eorb(iholep)) * alphfs
             do 190  i = 1, np
                xj0 = sin(xk0*dr(i))/(xk0*dr(i))
                xpc(i) = (xmult1*dgc0(i)*dpc(i,iorb,0)+
     1            xmult2*dpc0(i)*dgc(i,iorb,0)) * xj0
                xqc(i) = 0.0d0
  190        continue
c            xirf=lfin+linit+2
             xirf=2
             call somm (dr, xpc, xqc, hx, xirf, 0, np)
             oscstr(nosc) = xirf**2/3.0d0 
             enosc(nosc) = eorb(iorb)
             index(nosc) = iorb
          endif
        endif
  30  continue

c     write down information about oscillators
      if(open_16) then
        write(16, *) nosc
        do 210 i=1,nosc
          write(16,220) oscstr(i), enosc(i), index(i)
  220     format ( f9.5, f12.3, i4)
  210   continue
        call json_value_add(fp, 'nosc',   nosc)
        call json_value_add(fp, 'oscstr', oscstr(1:nosc))
        call json_value_add(fp, 'enosc',  enosc(1:nosc))
        call json_value_add(fp, 'index',  index(1:nosc))
      endif

c     calculate and write out f0(Q) on grid delq=0.5 Angstorm**(-1)
      dq=0.5*bohr 
      do 300 iq = 1,81
         xk0 = dq*(iq-1)
c        srho is 4*pi*density 
         do 560  i = 1, np
            xj0 = 1.d0
            if(iq.gt.1) xj0 = sin(xk0*dr(i))/(xk0*dr(i))
            xpc(i) = srho(i) * (dr(i)**2) *xj0
            xqc(i) = 0.d0
  560    continue
         xirf = 2.d0
         call somm (dr, xpc, xqc, hx, xirf, 0, np)
         if (open_16) then
            write (16, 570) 0.5*(iq-1), xirf
            qval(iq) = 0.5*(iq-1)
            xirfval(iq) = xirf
         endif
  570    format ( f5.1, 1x, f9.4)
  300 continue
      if (open_16) then
        call json_value_add(fp, 'qval',    qval)
        call json_value_add(fp, 'xirfval', xirfval)         
        open(newunit=iunit, file='fpf0.json', status='REPLACE')
        call json_print(fp,iunit)
        close(iunit)
        call json_destroy(fp)
      endif

      if (open_16) close(unit=16)

      return
      end
