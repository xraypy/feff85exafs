!=======================================================================
!     ATOMS
!=======================================================================

module atoms_inp
  ! The geom.dat file
  use dimsmod,only: natx,nphx,nheadx
  implicit none
  integer nat, nph, iatph(0:nphx), iphat(natx), ibounc(natx)
  ! ibounc is currently set to 1 for all atoms in ffsort.  Path uses it.  Probably discontinued variable but ah well. !KJ
  double precision  rat(3,natx)
  character(*),parameter,private :: filename='geom.dat'
  !		iphat(natx)  -  given specific atom, which unique pot?
  !		rat(3,natx)  -  cartesian coords of specific atom
  !		iatph(0:nphx)  - given unique pot, which atom is model?
  !                      (0 if none specified for this unique pot)

contains

  subroutine atoms_read
    !		Read  geom.dat file
    implicit none
    character*512 slog
    character*80 head(nheadx)
    integer lhead(nheadx),j1,nhead ! ,j
    real*8 rdum1(3)
    integer idum1,idum2

    call json_read_geom(nat, nph, iatph, rat, iphat, ibounc)
    
!     open (file=filename, unit=3, status='old')
!     !			read header
!     nhead = nheadx
!     call rdhead (3, nhead, head, lhead)
!     nat = 0
!     nph = 0
!     iatph(:)=0
! 50  continue
!     !KJ I switched up statements below so that code doesn't falsely abort when nat=natx.
!     !KJ			nat = nat+1
!     if (nat .gt. natx)  then
!        write(slog,'(a, 2i10)') ' nat, natx ', nat, natx
!        call wlog(slog)
!        stop 'Bad input'
!     endif
!     read(3,*,end=60) j1,rdum1(1:3),idum1,idum2
!     nat = nat+1
!     rat(1:3,nat)=rdum1(1:3)
!     iphat(nat)=idum1
!     ibounc(nat)=idum2
!     !KJ			read(3,*,end=60)  j1, (rat(j,nat),j=1,3), iphat(nat), ibounc(nat) !KJ j2  !KJ put ibounc back in for program PATH
!     if (iphat(nat).gt.nph) nph = iphat(nat)
!     if ( iatph(iphat(nat)).eq.0) iatph(iphat(nat)) = nat
!     goto 50
! 60  continue
!     !KJ			nat = nat-1
!     close(3)
  end subroutine atoms_read

end module atoms_inp
