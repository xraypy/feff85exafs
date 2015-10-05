program opconsAt
  use atoms_inp
  use potential_inp
  use geometry_inp
  use constants
  use AtomicPotIO
  use opcons_inp
  use par
  use errorfile
  implicit none
  integer iph, NComps ! , nnn, iph2, iX, iY, iZ2, iAt, iAt2
  integer, allocatable :: izComp(:), nAtComp(:)
  character(2), allocatable :: Components(:)
  character(12), allocatable :: EpsFiles(:), thisfile
  character(2),external :: GetElement
  real(8),allocatable :: rnrm(:)
  real(8) VTot ! , rNN, r, rCnt(3), point(3), x2, y2, z2, dX
  ! REAL(8),PARAMETER :: RTol  = (3.d0)**2
  ! INTEGER,PARAMETER :: NGrid = 200

  integer i1, i2
  integer, parameter :: epsmax = 700
  real(8) thiseps(epsmax,3)
  
  !KJ 1-2012:
  call par_begin
  if (worker) go to 400
  call OpenErrorfileAtLaunch('opconsat')

  !CALL opcons_read
  !IF(.NOT.run_opcons) STOP
  !if (.not. run_opcons) goto 400
  call opcons_init

  call atoms_read
  call potential_read

  allocate(izComp(0:nph), Components(0:nph), nAtComp(0:nph), EpsFiles(0:nph), rnrm(0:nph), thisfile)

  call ReadAtomicPots(nph, rnrm)
  print*, 'nph, rnrm', nph, rnrm

  ! Find the number density of each component.
  VTot = 0.d0
  do iph = 0, nph
     VTot = VTot + xnatph(iph)*4.d0/3.d0*pi*(rnrm(iph)*bohr)**3
  end do

  do iph = 0, nph
     if(NumDens(iph).lt.0.d0) NumDens(iph) = xnatph(iph)/VTot
     ! test with numDens = 1. Should give same loss as input file for a single file.
     if(.false.) then
        if(iph.ne.0) then
           NumDens(iph) = 1.d0
        else
           NumDens(iph) = 0.d0
        end if
     end if
     Components(iph) = GetElement(iz(iph))
  end do

  ! Get opcons{Element}.dat from database.
  !PRINT*, iz, nph
  do iph = 0, nph
     do i1 = 1,epsmax
        do i2 = 1,3
           thiseps(i1,i2) = 0.d0
        end do
     end do
     EpsFiles(iph) = 'opcons' // trim(adjustl(Components(iph))) // '.dat'
     print*, Components(iph), NumDens(iph), EpsFiles(iph)
     call epsdb(iz(iph), thiseps)
     call write_eps(iz(iph), thiseps, EpsFiles(iph))
  end do

  NComps = nph + 1
  call AddEps(EpsFiles,NComps,NumDens,print_eps)

  deallocate(izComp, Components, nAtComp, EpsFiles, rnrm, thisfile)
  
  !KJ 1-2012:
400 call par_barrier
  call par_end
  call WipeErrorfileAtFinish
  stop


end program opconsAt
