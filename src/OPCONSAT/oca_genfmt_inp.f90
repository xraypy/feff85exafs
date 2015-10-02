!=======================================================================
!     GENFMT
!=======================================================================

module genfmt_inp
  use global_inp
  implicit none
  character(*),parameter,private :: filename='genfmt.inp'
  integer  mfeff, ipr5, iorder
  logical  wnstar
  double precision critcw

contains

  subroutine genfmt_write
    open (file=filename, unit=3, status='unknown')
    write(3,10) 'mfeff, ipr5, iorder, critcw, wnstar'
    write(3,180)  mfeff, ipr5, iorder, critcw, wnstar
180 format ( 2i4, i8, f13.5, L5)
    !KJ 7-09 Next 2 lines for feff8q
    write(3,'(a24)') ' the number of decomposition channels'
    write(3,'(i5)') ldecmx
    close(3)
    ! standard formats for string, integers and real numbers
10  format(a)
    ! 20  format (20i4)
    ! 30  format (9f13.5)
  end subroutine genfmt_write

  subroutine genfmt_read
    open (file=filename, unit=3, status='old')
    read(3,*) ; read(3,*)  mfeff, ipr5, iorder, critcw, wnstar
    !KJ 7-09 Next line for feff8q
    read(3,*) ; read(3,*) ldecmx
    close(3)
  end subroutine genfmt_read

  subroutine genfmt_init
    mfeff = 1
    ipr5 = 0
    iorder = 2
    wnstar = .false.
    critcw = 4*1.d0
  end subroutine genfmt_init

end module genfmt_inp
