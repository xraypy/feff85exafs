!=======================================================================
!     FF2X
!=======================================================================

module ff2x_inp
  use global_inp
  use xsph_inp
  use fms_inp
  use genfmt_inp
  implicit none
  character(*),parameter,private :: filename='ff2x.inp'
  integer  mchi, ipr6, mbconv, absolu !KJ added absolu 3-06
  double precision  vrcorr, vicorr, s02, alphat, thetae


contains

  subroutine ff2x_write
    integer i
    open (file=filename, unit=3, status='unknown')
    write(3,10) 'mchi, ispec, idwopt, ipr6, mbconv, absolu, iGammaCH' !KJ added absolu 3-06
    write(3,20)  mchi, ispec, idwopt, ipr6, mbconv, absolu, iGammaCH !KJ added absolu 3-06
    write(3,10) 'vrcorr, vicorr, s02, critcw'
    write(3,30)  vrcorr, vicorr, s02, critcw
    write(3,10) 'tk, thetad, alphat, thetae, sig2g'
    write(3,30)  tk, thetad, alphat, thetae, sig2g
    !KJ 7-09 next 4 lines for feff8q
    write(3,10) 'momentum transfer'
    write(3, '(3f13.5)') (xivec(i),i=1,3)
    write(3,'(a24)') ' the number of decomposition channels'
    write(3,'(i5)') ldecmx
    close(3)
    ! standard formats for string, integers and real numbers
10  format(a)
20  format (20i4)
30  format (9f13.5)
  end subroutine ff2x_write

  subroutine ff2x_read
    integer i
    open (file=filename, unit=3, status='old')
    read(3,*) ; read(3,*)  mchi, ispec, idwopt, ipr6, mbconv, absolu, iGammaCH !KJ added absolu 3-06
    read(3,*) ; read(3,*)  vrcorr, vicorr, s02, critcw
    read(3,*) ; read(3,*)  tk, thetad, alphat, thetae, sig2g
    read(3,*) ; read(3, *) (xivec(i),i=1,3)
    read(3,*) ; read(3,*) ldecmx
    close(3)
  end subroutine ff2x_read

  subroutine ff2x_init
    absolu=0  !KJ 3-06 for ABSOLUTE card
    mchi = 1
    ipr6 = 0
    mbconv = 0
    vicorr = 0.d0
    vrcorr = 0.d0
    s02 = 1.d0
    alphat = 0.d0
    thetae = 0.d0
  end subroutine ff2x_init

end module ff2x_inp
