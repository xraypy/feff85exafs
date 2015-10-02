!=======================================================================
!     COMPTON
!=======================================================================


module compton_inp
  implicit none
  character(*),parameter,private :: filename='compton.inp'

  ! spatial and momentum grid parameters
  integer :: ns, nphi, nz, nzp, npq
  real :: smax, phimax, zmax, zpmax, pqmax

  ! flags
  logical :: do_compton, do_rhozzp
  logical :: force_jzzp
  integer run_compton_module

  ! apodization function type
  integer :: window
  real :: window_cutoff

  real :: temperature
  logical ::  set_chemical_potential
  real :: chemical_potential

  integer, parameter :: WINDOW_STEP = 0, WINDOW_HANNING = 1
contains
  subroutine compton_write
    if (do_compton .or. do_rhozzp) then
       run_compton_module=1
    else
       run_compton_module=0
    endif
    open (file=filename, unit=3, status='unknown')
    write(3,10) 'run compton module?'
    write(3,*)  run_compton_module 
    write(3,10) 'pqmax, npq'
    write(3,*) pqmax, npq
    write(3,10) 'ns, nphi, nz, nzp'
    write(3,20) ns, nphi, nz, nzp
    write(3,10) 'smax, phimax, zmax, zpmax'
    write(3,30) smax, phimax, zmax, zpmax
    write(3,10) 'jpq? rhozzp? force_recalc_jzzp?'
    write(3,*) do_compton, do_rhozzp, force_jzzp
    write(3,10) 'window_type (0=Step, 1=Hann), window_cutoff'
    write(3,*) window, window_cutoff
    write(3,10) 'temperature (in eV)'
    write(3,30) temperature
    write(3,10) 'set_chemical_potential? chemical_potential(eV)'
    write(3,*) set_chemical_potential, chemical_potential
    close(3)
    ! standard formats for string, integers, real numbers
10  format(a)
20  format (20i4)
30  format (9f13.5)
  end subroutine compton_write

  subroutine compton_read
    open (file=filename, unit=3, status='old',err=100)
    read(3,*,end=100,err=100) ; read(3,*,end=100,err=100) run_compton_module
    read(3,*,end=100,err=100) ; read(3,*,end=100,err=100) pqmax, npq
    read(3,*,end=100,err=100) ; read(3,*,end=100,err=100) ns, nphi, nz, nzp
    read(3,*,end=100,err=100) ; read(3,*,end=100,err=100) smax, phimax, zmax, zpmax
    read(3,*,end=100,err=100) ; read(3,*,end=100,err=100) do_compton, do_rhozzp, force_jzzp
    read(3,*,end=100,err=100) ; read(3,*,end=100,err=100) window, window_cutoff
    read(3,*,end=100,err=100) ; read(3,*,end=100,err=100) temperature
    read(3,*,end=100,err=100) ; read(3,*,end=100,err=100) set_chemical_potential, chemical_potential
    close(3)
    return
100 run_compton_module=0  ! no compton.inp -> don't do compton
    return
  end subroutine compton_read

  subroutine compton_init
    real, parameter :: pi  = 3.1415926535897932384626433832795
    ns   = 32
    nphi = 32
    nz   = 32
    nzp  = 144

    smax   = 0
    phimax = 2*pi
    zmax   = 0
    zpmax  = 10.0

    npq   = 1000
    pqmax = 5.0

    do_compton     = .false.
    do_rhozzp  = .false.
    force_jzzp = .false.
    run_compton_module=0

    window = WINDOW_HANNING
    window_cutoff = 0

    temperature = 0.0
    set_chemical_potential = .false.
    chemical_potential = 0
  end subroutine compton_init
end module compton_inp
