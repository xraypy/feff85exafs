!=======================================================================
!     RECIPROCAL
!=======================================================================

module reciprocal_inp
  !     k-space variables :
  use controls  !KJ 8/06
  use struct, nphstr => nph
  use kklist,only: nkp,usesym,nkx,nky,nkz,ktype
  use strfacs,only: streta,strrmax,strgmax,init_strfacs
  implicit none
  integer icorehole
  real*8 streimag ! additional broadening for calculation KKR structure factors ; not recommended
  character(*),parameter,private :: filename='reciprocal.inp'

contains

  subroutine reciprocal_write
    !KJ next file added 8/06
    integer i
    open (file=filename, unit=3, status='unknown')
    !       in which space are we?
    write(3,10) 'spacy'
    write(3,20) spacy
    if(spacy.eq.0) then
       write(3,10) 'lattice vectors  (in A, in Carthesian coordinates)'
       write(3,30) a1
       write(3,30) a2
       write(3,30) a3
       write(3,10) 'Volume scaling factor (A^3); eimag; core hole'
       write(3,30) dble(-1),dble(0),dble(1)
       write(3,10) 'lattice type  (P,I,F,R,B,CXY,CYZ,CXZ)'
       write(3,10) latticename
       write(3,10) '#atoms in unit cell ; position absorber ; corehole?'
       write(3,20) nats,absorber,icorehole
       write(3,10) '# k-points total/x/y/z ; ktype; use symmetry?'
       write(3,*) nkp,nkx,nky,nkz,ktype,usesym  ! format line 20 limits integer to 4 positions - not enough for nkp!
       write(3,10) 'ppos'
       do i=1,nats
          write(3,30) ppos(:,i)
       enddo
       write(3,10) 'ppot'
       !KJ bugfix 5/2012: It's important not to use formatting when there are more atoms than fit on one line!!			   
       write(3,*) ppot
       write(3,10) 'streta,strgmax,strrmax'
       write(3,30) streta,strgmax,strrmax
    endif
    close(3)
    ! standard formats for string, integers and real numbers
10  format(a)
20  format (20i4)
30  format (9f13.5)
  end subroutine reciprocal_write

  subroutine reciprocal_read(celvin)
    use struct, nphstr => nph
    integer i
    real*8,intent(out) :: celvin
    open (3,file=filename,status='unknown',err=167)
    read(3,*,end=167,err=167)
    read(3,*,end=167,err=167) spacy
    if(spacy.eq.0) then
       read(3,*) ; read(3,*) a1(:)
       read(3,*) a2(:)
       read(3,*) a3(:)
       read(3,*) ; read(3,*) celvin,streimag,cholestrength
       read(3,*) ; read(3,*) latticename
       lattice=latticename(1:1)
       read(3,*) ; read(3,*) nats,absorber,icorehole
       read(3,*) ; read(3,*) nkp,nkx,nky,nkz,ktype,usesym
       read(3,*)
       !Careful: the next statement used to be "if size(ppot).eq.0".  However, on ifort size(ppot)=0 but on gfortran it =1!!
       !Hence the new instruction.
       !I wish if(allocated(ppot)) would work here; I don't understand why it doesn't.
       if(size(ppot).lt.nats) call init_struct(nats) !KJ 7-09 bugfix call this only once ; I can't seem to use "allocated(ppos)" here?
       do i=1,nats
          read(3,*) ppos(:,i)
       enddo
       read(3,*) ; read(3,*) ppot
       read(3,*) ; read(3,*) streta,strgmax,strrmax
       if(icorehole.eq.1) then
          corehole=.true.
       else
          corehole=.false.
       endif
    endif
    return
167 spacy=1
    return
  end subroutine reciprocal_read

  subroutine reciprocal_init
    call init_controls
    call init_strfacs
    icorehole = 1  ! use core hole
    streimag = dble(0) ! no extra broadening for KKR struc factors
    cholestrength = dble(1) ! don't mess with core hole
  end subroutine reciprocal_init

end module reciprocal_inp
