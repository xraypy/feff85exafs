!=======================================================================
!     SFCONV
!=======================================================================

module sfconv_inp
  use global_inp, only : ispec
  use ff2x_inp, only : ipr6
  implicit none
  character(*),parameter,private :: filename='sfconv.inp'
  integer  msfconv, ipse, ipsk
  double precision wsigk, cen
  character(12) cfname

contains

  subroutine sfconv_write
    !c    sfconv.inp - Josh Kas
    open (file=filename, unit=3, status='unknown')
    write(3,10) 'msfconv, ipse, ipsk'
    write(3,20)  msfconv, ipse, ipsk
    write(3,10) 'wsigk, cen'
    write(3,30) wsigk, cen
    write(3,10) 'ispec, ipr6'
    write(3,20)  ispec, ipr6
    write(3,10) 'cfname'
    write(3,10) cfname
    close(3)
    ! standard formats for string, integers and real numbers
10  format(a)
20  format (20i4)
30  format (9f13.5)
  end subroutine sfconv_write

  subroutine sfconv_read
    open (file=filename, unit=3, status='old')
    read (3,*) ; read (3,*)  msfconv, ipse, ipsk
    read (3,*) ; read (3,*)  wsigk, cen
    read (3,*) ; read (3,*)  ispec, ipr6
    read (3,*) ; read (3,*)  cfname
    close(3)
  end subroutine sfconv_read

  subroutine sfconv_init
    msfconv = 0 ! Josh Kas
    ipse = 0
    ipsk = 0
    wsigk = 0.d0 ! Josh Kas
    cen = 0.d0 ! Josh Kas
    cfname = 'NULL'
  end subroutine sfconv_init

end module sfconv_inp
