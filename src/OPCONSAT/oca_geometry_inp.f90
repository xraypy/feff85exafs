!=======================================================================
!     GEOMETRY
!=======================================================================

module geometry_inp
  use dimsmod, only: nattx
  implicit none
  !c    atoms.dat
  integer  natt
  integer iphatx(nattx)
  double precision  ratx(3,nattx)

contains

  subroutine geometry_write_atoms
    integer iat
    double precision distance
    !c    atoms.dat to be read by ffsort, which will write smaller geom.dat file
    open (file='atoms.dat', unit=3, status='unknown')
    write (3, 35) natt
35  format ('natx =  ', i7)
    write (3, 10) '    x       y        z       iph  '
    do iat = 1, natt
       distance=dsqrt((ratx(1,iat)-ratx(1,1))**2+(ratx(2,iat)-ratx(2,1))**2+(ratx(3,iat)-ratx(3,1))**2) ! core hole should be in position 1 by now
       write(3,36) ratx(1,iat), ratx(2,iat), ratx(3,iat), iphatx(iat), distance
36     format( 3f13.5, i4, f13.5)
    enddo
    close(3)
10  format(a)
    ! 20  format (20i4)
    ! 30  format (9f13.5)
  end subroutine geometry_write_atoms

  subroutine geometry_init
    natt = 0
    iphatx(:) = -1
    ratx(:,:) = 0.d0 !KJ added 7-09
  end subroutine geometry_init

end module geometry_inp
