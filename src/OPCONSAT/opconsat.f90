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
  character(2), allocatable :: Components(:)
  ! character(12), allocatable :: EpsFiles(:)
  character(12), allocatable :: thisfile
  character(2),external :: GetElement
  real(8),allocatable :: rnrm(:)
  real(8) VTot ! , rNN, r, rCnt(3), point(3), x2, y2, z2, dX
  ! REAL(8),PARAMETER :: RTol  = (3.d0)**2
  ! INTEGER,PARAMETER :: NGrid = 200

  character(12) EpsFile
  integer i1, i2, NDataTot, iUnit, NPoles
  integer, parameter :: epsmax = 700
  real(8) thiseps(epsmax,3,0:nphx)
  real(8) energy(epsmax), loss(epsmax)
  real(8) g(10000), gamma, omi(10000), Delta(10000), eps0, sumrl, xNElec, csumrl
  character ch

  logical excinp
  logical write_loss, write_opcons, write_exc, verbose
  write_loss   = .true.
  write_opcons = .true.
  write_exc    = .true.
  verbose      = .true.
  NPoles       = 100
  eps0         = -1.d0
  
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

  allocate(Components(0:nph), rnrm(0:nph), thisfile)
  ! allocate(EpsFiles(0:nph))

  
  call ReadAtomicPots(nph, rnrm)
  if (verbose) print*, 'nph, rnrm', nph, rnrm

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
  do iph = 0, nph
     do i1 = 1,epsmax
        do i2 = 1,3
           thiseps(i1,i2,iph) = 0.d0
        end do
     end do
     EpsFile = 'opcons' // trim(adjustl(Components(iph))) // '.dat'
     if (verbose) print*, Components(iph), NumDens(iph), EpsFile
     call epsdb(iz(iph), thiseps(:,:,iph))
     if (write_opcons) call write_eps(iz(iph), thiseps(:,:,iph), EpsFile)
  end do

  NComps = nph + 1
  call AddEps(NComps, NumDens, thiseps, NDataTot, energy, loss)
  
  ! Open output file and write data.
  if (write_loss) then
     call OpenFl('loss.dat')
     call GetIOFileInfo('loss.dat', UnitNumber = iUnit)
     write(iUnit,'(A)') '# E(eV)    Loss'
     do i1 = 1, NDataTot
        write(iUnit,*), energy(i1), loss(i1)
     end do
     call CloseFl('loss.dat')
  end if
  
  deallocate(Components, rnrm, thisfile)


  !! sumrl, xNElec, gamma, and csumrl -- a bit cryptic
  sumrl = 1.d0
  xNElec = 1.d0
  inquire( file='exc.inp', exist=excinp ) 
  if (excinp) then
     open(unit=61, file='exc.inp', status='old')
     read(61,*) NPoles
     read(61,*) ch
     if (ch.eq.'y'.or.ch.eq.'Y') then
        eps0 = -1.d0
     else
        read(61,*) ch
        if(ch.eq.'n'.or.ch.eq.'N') then
           eps0 = -2.d0
        else
           read(61,*) eps0
        end if
     end if
     close(61)
  else
     print*, '# Enter number of poles:'
     read*, NPoles
     print*, 'Is this a metal? (y/n)'
     read*, ch
     if(ch.eq.'y'.or.ch.eq.'Y') then
        eps0 = -1.d0
     else
        print*, 'Would you like to set the dielectric constant? (y/n)'
        read*, ch
        if(ch.eq.'n'.or.ch.eq.'N') then
           eps0 = -2.d0
        else
       
           ! This input can be used to correct the dielectric constant,
           ! which is related to the inverse moment.
           ! Use eps0 = -2 to ignore this correction.
           print*, '# Enter dielectric constant: '
           read*, eps0
        end if
     end if
     if (verbose) print*, eps0
  end if
  
  gamma = 0.01
  csumrl= xNElec/sumrl
  do i1 = 1, NDataTot
     loss(i1) = loss(i1)*csumrl
  end do

  ! getomi finds poles and weights
  call getomi(energy, loss, NDataTot, NPoles, omi, g, Delta, eps0)

  if (write_exc) then
     open(unit=13,file='exc.dat',status='replace')
     write(13,'(A33,I4,A6)') '# Loss function represented with ', NPoles, ' poles'
     write(13,'(A23,f8.4)') '# Dielectric constant: ', eps0
     do i1 = 1, NPoles
        write(13,'(4f20.10)'), omi(i1), omi(i1)*gamma, g(i1), Delta(i1)
     end do
     close(13)
  end if
  
  !KJ 1-2012:
400 call par_barrier
  call par_end
  call WipeErrorfileAtFinish
  stop


end program opconsAt
