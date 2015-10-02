!=======================================================================
!     PATHS
!=======================================================================

module paths_inp
  use ldos_inp
  implicit none
  character(*),parameter,private :: filename='paths.inp'
  integer  mpath, ms, nncrit, nlegxx, ipr4, ica  !KJ added ica 6-06
  !KJ nncrit seems to be bogus input, i.e. not set in rdinp at all ; fully internal to PATH.  Set to 0 and kept here for compatibility.
  real critpw, pcritk, pcrith,  rmax

contains

  subroutine paths_write
    implicit none
    open (file=filename, unit=3, status='unknown')
    write(3,10) 'mpath, ms, nncrit, nlegxx, ipr4'
    write(3,20)  mpath, ms, nncrit, nlegxx, ipr4
    write(3,10) 'critpw, pcritk, pcrith,  rmax, rfms2'
    write(3,30)  critpw, pcritk, pcrith,  rmax, rfms2
    write(3,10) 'ica' !KJ 6-06
    write(3,20)  ica  !KJ 6-06
    close(3)
    ! standard formats for string, integers and real numbers
10  format(a)
20  format (20i4)
30  format (9f13.5)
  end subroutine paths_write

  subroutine paths_read
    open (file=filename, unit=3, status='old')
    read(3,*) ; read(3,*)  mpath, ms, nncrit, nlegxx, ipr4
    read(3,*) ; read(3,*)  critpw, pcritk, pcrith,  rmax, rfms2
    read(3,*) ; read(3,*)  ica  !KJ 6-06
    close(3)
  end subroutine paths_read

  subroutine paths_init
    mpath = 1
    ms = 0
    ipr4 = 0
    ica=-1 !KJ 6-06
    critpw = 2.5*1.e0
    pcritk = 0.e0
    pcrith = 0.e0
    rmax = -1 * 1.e0
    nlegxx = 10
    nncrit = 0
  end subroutine paths_init

end module paths_inp
