!=======================================================================
!     OPCONS
!=======================================================================
MODULE opcons_inp
  USE dimsmod
  LOGICAL run_opcons, print_eps
  REAL(8) NumDens(0:nphx)
  character(*),parameter,private :: filename='opcons.inp'

CONTAINS
  SUBROUTINE opcons_init
    run_opcons = .FALSE.
    print_eps  = .FALSE.
    NumDens(:) = -1.d0
  END SUBROUTINE opcons_init

  SUBROUTINE opcons_write
    ! INTEGER iph

    OPEN(FILE=filename,UNIT=8,STATUS='REPLACE')

    WRITE(8,'(A)') 'run_opcons'
    WRITE(8,*) run_opcons
    WRITE(8,'(A)') 'print_eps'
    WRITE(8,*) print_eps
    WRITE(8,'(A)') 'NumDens(0:nphx)'
    WRITE(8,*) NumDens(0:nphx)

    CLOSE(8)
  END SUBROUTINE opcons_write

  SUBROUTINE opcons_read
    ! INTEGER iph

    OPEN(FILE=filename,UNIT=8,STATUS='OLD')

    READ(8,*)
    READ(8,*) run_opcons
    READ(8,*)
    READ(8,*) print_eps
    READ(8,*)
    READ(8,*) NumDens(0:nphx)
  END SUBROUTINE opcons_read

END MODULE opcons_inp
