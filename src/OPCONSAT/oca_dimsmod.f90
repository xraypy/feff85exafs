!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: m_dimsmod.f90,v $:
! $Revision: 1.30 $
! $Author: jorissen $
! $Date: 2012/06/29 01:05:24 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module DimsMod
  ! This module contains dimensions for data arrays

  ! The file in which dimensions current to the calculation are saved :
  character*20, parameter :: dimFName = '.dimensions.dat'
  ! Set the following according to max. available memory on your system :
  ! The hardcoded limit on cluster size that can NEVER be exceeded :
  integer, parameter :: nclusxhardlimit = 3000
  ! The hardcoded upper limit on l-values that can NEVER be exceeded :
  integer, parameter :: lxhardlimit = 20
  ! The hardcoded upper limit on the number of potentials that can NEVER be exceeded :
  integer, parameter :: nphxhardlimit = 31
  private dimFName, lxhardlimit, nclusxhardlimit !! Meaning no other code can change these ...

  integer,parameter :: nclxtd = 100     ! Maximum number of atoms for tdlda module.
  integer,parameter :: nspx   = 1       ! Max number of spins: 1 for spin average; 2 for spin-dep
  integer,parameter :: natx   = 1000    ! Max number of atoms in problem for the pathfinder and ffsort
  integer,parameter :: nattx  = 1000   ! Max number of atoms in problem for the rdinp
  integer,parameter :: nphx   = 11      ! Max number of unique potentials (potph)
  integer,parameter :: ltot   = 24      ! Max number of ang mom (arrays 1:ltot+1)
  integer,parameter :: nrptx  = 1251    ! Loucks r grid used through overlap and in phase work arrays
  integer,parameter :: nex    = 150     ! Number of energy points genfmt, etc.
  integer,parameter :: lamtot = 15      ! Max number of distinct lambda's for genfmt 15 handles iord 2 and exact ss
  integer,parameter :: mtot   = 4       ! Vary mmax and nmax independently
  integer,parameter :: ntot   = 2 
  integer,parameter :: npatx  = 8       ! Max number of path atoms, used  in path finder, NOT in genfmt
  integer,parameter :: legtot = npatx+1 ! Matches path finder, used in GENFMT
  integer,parameter :: novrx  = 8       ! Max number of overlap shells (OVERLAP card)
  integer,parameter :: nheadx = 30      ! Max number of header lines
  !integer,parameter :: nheadx = 20+nphxhardlimit      ! Max number of header lines !KJ 7-09 added term to accomodate large systems in xsect.bin header
  integer,parameter :: MxPole = 1000    ! Max number of poles that can be used to model epsilon^-1 for HL multipole self energy
  integer,parameter :: nwordx = max(100,2+2*nphxhardlimit)     ! An infuriatingly stupid parameter that shows up in a few places. KJ added 7-09.  used to be 20 - must be at least 2*(1+nphx) for feff.bin header.
  integer,parameter :: novp = 40 ! For istprm, movrlp, ovp2mt - an atom list cutoff that should be high enough to include one atom of each potential type.  Added 2-2011 !KJ

  ! NON PARAMETER STATEMENTS
  integer :: nclusx    ! Maximum number of atoms for FMS.
  integer :: lx = 4       ! Max orbital momentum for FMS module.

  ! OLD XPARAM.H MODULE
  integer,parameter :: natxx = natx
  integer,parameter :: nexx = nex
  integer,parameter :: nkmin = 1
  integer,parameter :: nphasx = nphx
  integer :: istatx

contains

  subroutine write_dimensions(nclusxuserlimit,lxuserlimit)
    implicit none
    ! Write dimension data to a file
    integer,intent(in) :: nclusxuserlimit,lxuserlimit
    integer :: ios  ! IO Status

    !   3/ Apply hardcoded dimension limits
    if(nclusxuserlimit.gt.0) then
       nclusx=min(nclusx,nclusxuserlimit)
    else
       nclusx=min(nclusx,nclusxhardlimit)
    endif
    if(lxuserlimit.ge.0) then
       lx=min(lx,lxuserlimit)
    else
       lx=min(lx,lxhardlimit)
    endif
    open(10,FILE=trim(dimFName),STATUS='unknown',FORM='formatted',IOSTAT=ios)
    call chopen(ios,trim(dimFName),'dimsmod')
    if (ios.ne.0) stop "Error writing dimensions.dat.  Quiting."

    write(10,*) nclusx,lx
    close(10)
  end subroutine write_dimensions


  subroutine init_dimensions
    implicit none
    ! Read dimensions from file
    integer :: ios  ! IO Status
    open(10,FILE=trim(dimFName),STATUS='old',FORM='formatted',IOSTAT=ios)
    call chopen(ios,trim(dimFName),'dimsmod')
    read(10,*) nclusx,lx
    close(10)
    ! OLD XPARAM DIMENSIONS
    istatx=(lx+1)**2*nclusx*nspx
  end subroutine init_dimensions


end module DimsMod

