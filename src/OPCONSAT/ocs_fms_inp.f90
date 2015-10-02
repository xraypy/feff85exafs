!=======================================================================
!     FMS
!=======================================================================

module fms_inp
  use ldos_inp
  use global_inp,only: ldecmx
  implicit none
  character(*),parameter,private :: filename='fms.inp'
  integer mfms, idwopt, ipr3 !ipr3 is currently dummy - not in fms.inp
  real rprec
  !KJ rprec seems to be bogus input, i.e. not used anywhere in entire FEFF90.  Set to 0 and kept here for compatibility.
  double precision   tk, thetad, sig2g

contains

  subroutine fms_write
    implicit none
    integer iph
    open (file=filename, unit=3, status='unknown')
    write(3,10) 'mfms, idwopt, minv'
    write(3,20)  mfms, idwopt, minv
    write(3,10) 'rfms2, rdirec, toler1, toler2'
    write(3,30)  rfms2, rdirec, toler1, toler2
    write(3,10) 'tk, thetad, sig2g'
    write(3,30)  tk, thetad, sig2g
    write(3,10) ' lmaxph(0:nph)'
    write(3,20)  (lmaxph(iph),iph=0,nph)
    !KJ 7-09 Next 2 lines for feff8q
    write(3,'(a24)') ' the number of decomposition channels'
    write(3,'(i5)') ldecmx
    close(3)
    ! standard formats for string, integers and real numbers
10  format(a)
20  format (20i4)
30  format (9f13.5)
  end subroutine fms_write

  subroutine fms_read
    integer iph
    open (file=filename, unit=3, status='unknown')
    read(3,*) ; read(3,*)  mfms, idwopt, minv
    read(3,*) ; read(3,*)  rfms2, rdirec, toler1, toler2
    read(3,*) ; read(3,*)  tk, thetad, sig2g
    read(3,*) ; read(3,*)  (lmaxph(iph),iph=0,nph)
    !KJ 7-09 Next line for feff8q
    read(3,*) ; read(3,*) ldecmx
    close(3)
  end subroutine fms_read

  subroutine fms_init
    mfms = 1
    idwopt = -1
    sig2g = 0.d0
    thetad = 0.d0
    tk = 0.d0
    ipr3 = 0
    rprec = 0.e0
  end subroutine fms_init

end module fms_inp
