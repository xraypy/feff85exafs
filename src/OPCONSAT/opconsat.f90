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
  character(12), allocatable :: EpsFiles(:)
  character(2),external :: GetElement
  real(8),allocatable :: rnrm(:)
  real(8) VTot ! , rNN, r, rCnt(3), point(3), x2, y2, z2, dX
  ! REAL(8),PARAMETER :: RTol  = (3.d0)**2
  ! INTEGER,PARAMETER :: NGrid = 200

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

  allocate(izComp(0:nph), Components(0:nph), nAtComp(0:nph), EpsFiles(0:nph), rnrm(0:nph))

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
  call epsdb(iz,nph+1)

  ! Get the file names.
  do iph = 0, nph
     EpsFiles(iph) = 'opcons' // trim(adjustl(Components(iph))) // '.dat'
     print*, Components(iph), NumDens(iph), EpsFiles(iph)
  end do
  NComps = nph + 1
  call AddEps(EpsFiles,NComps,NumDens,print_eps)

  !KJ 1-2012:
400 call par_barrier
  call par_end
  call WipeErrorfileAtFinish
  stop


end program opconsAt


      subroutine bwords_nc (s, nwords, words)
!
!     Breaks string into words.  Words are seperated by one or more
!     blanks or tabs, but not a comma.
!
!     ARGS        I/O      DESCRIPTION
!     ----        ---      -----------
!     S            I       CHAR*(*)  String to be broken up
!     NWORDS      I/O      Input:  Maximum number of words to get
!                          Output: Number of words found
!     WORDS(NWORDS) O      CHAR*(*) WORDS(NWORDS)
!                          Contains words found.  WORDS(J), where J is
!                          greater than NWORDS found, are undefined on
!                          output.
!
!      Written by:  Steven Zabinsky, September 1984
!      Tab char added July 1994.
!
!**************************  Deo Soli Gloria  **************************

!  -- No floating point numbers in this routine.
      implicit integer (a-z)

      character*(*) s, words(nwords)

      character blank, tab
      parameter (blank = ' ', tab = '	')
!     there is a tab character here               ^.

!  -- BETW    .TRUE. if between words
!     COMFND  .TRUE. if between words and a comma has already been found
      logical betw, comfnd

!  -- Maximum number of words allowed
      wordsx = nwords

!  -- SLEN is last non-blank character in string
      slen = istrln (s)

!  -- All blank string is special case
      if (slen .eq. 0)  then
         nwords = 0
         return
      endif

!  -- BEGC is beginning character of a word
      begc = 1
      nwords = 0

      betw   = .true.
      comfnd = .true.

      do 10  i = 1, slen
         if (s(i:i) .eq. blank .or. s(i:i) .eq. tab)  then
            if (.not. betw)  then
               nwords = nwords + 1
               words (nwords) = s (begc : i-1)
               betw = .true.
               comfnd = .false.
            endif
         else
            if (betw)  then
               betw = .false.
               begc = i
            endif
         endif

         if (nwords .ge. wordsx)  return

   10 continue

      if (.not. betw  .and.  nwords .lt. wordsx)  then
         nwords = nwords + 1
         words (nwords) = s (begc :slen)
      endif
 
      return
      end
