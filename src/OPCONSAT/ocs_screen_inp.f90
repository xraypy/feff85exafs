!=======================================================================
!     SCREEN
!=======================================================================

module screen_inp
  use atoms_inp,only: nph
  implicit none

  TYPE ScreenInputVars
     integer ner, nei, maxl, irrh, iend, lfxc, nrptx0
     double precision emin, emax, eimax, ermin, rfms
  END TYPE ScreenInputVars

  character(*),parameter,private :: filename='screen.inp'
  TYPE(ScreenInputVars) ScreenI

contains

  subroutine screen_write
    open(unit=3,file=filename,status='unknown')
    write(3,*) 'ner',ScreenI%ner
    write(3,*) 'nei',ScreenI%nei
    write(3,*) 'maxl',ScreenI%maxl
    write(3,*) 'irrh',ScreenI%irrh
    write(3,*) 'iend',ScreenI%iend
    write(3,*) 'lfxc',ScreenI%lfxc
    write(3,*) 'emin',ScreenI%emin
    write(3,*) 'emax',ScreenI%emax
    write(3,*) 'eimax',ScreenI%eimax
    write(3,*) 'ermin',ScreenI%ermin
    write(3,*) 'rfms',ScreenI%rfms
    write(3,*) 'nrptx0',ScreenI%nrptx0
    close(3)
    return
  end subroutine screen_write

  subroutine screen_inp_parse(str,vars)
    implicit none
    character*3,intent(in) :: str
    real*8,intent(in) ::  vars
    if (str .eq. 'ner') then
       ScreenI%ner   = int(vars)
    elseif (str .eq. 'nei') then
       ScreenI%nei   = int(vars)
    elseif (str .eq. 'max') then
       ScreenI%maxl  = int(vars)
    elseif (str .eq. 'irr') then
       ScreenI%irrh  = int(vars)
    elseif (str .eq. 'ien') then
       ScreenI%iend  = int(vars)
    elseif (str .eq. 'lfx') then
       ScreenI%lfxc  = int(vars)
    elseif (str .eq. 'emi') then
       ScreenI%emin  = vars
    elseif (str .eq. 'ema') then
       ScreenI%emax  = vars
    elseif (str .eq. 'eim') then
       ScreenI%eimax = vars
    elseif (str .eq. 'erm') then
       ScreenI%ermin = vars
    elseif (str .eq. 'rfm') then
       ScreenI%rfms  = vars
    elseif (str .eq. 'nrp')then
       ScreenI%nrptx0  = int(vars)
    else
       call wlog("Unrecognized keyword submitted to screen.inp in SCREEN_INP_PARSE ; aborting.")
       stop
    endif
    return
  end subroutine screen_inp_parse

  subroutine screen_inp_parse_and_write(str,vars)
    !KJ No longer used (1-2012).  Used in a previous version of feff.
    implicit none
    character*3,intent(in) :: str
    real*8,intent(in) ::  vars
    open(unit=3,file=filename,status='unknown',access='append')
    if (str .eq. 'ner') then
       ScreenI%ner   = int(vars)
       write(3,*) 'ner',ScreenI%ner
    elseif (str .eq. 'nei') then
       ScreenI%nei   = int(vars)
       write(3,*) 'nei',ScreenI%nei
    elseif (str .eq. 'max') then
       ScreenI%maxl  = int(vars)
       write(3,*) 'maxl',ScreenI%maxl
    elseif (str .eq. 'irr') then
       ScreenI%irrh  = int(vars)
       write(3,*) 'irrh',ScreenI%irrh
    elseif (str .eq. 'ien') then
       ScreenI%iend  = int(vars)
       write(3,*) 'iend',ScreenI%iend
    elseif (str .eq. 'lfx') then
       ScreenI%lfxc  = int(vars)
       write(3,*) 'lfxc',ScreenI%lfxc
    elseif (str .eq. 'emi') then
       ScreenI%emin  = vars
       write(3,*) 'emin',ScreenI%emin
    elseif (str .eq. 'ema') then
       ScreenI%emax  = vars
       write(3,*) 'emax',ScreenI%emax
    elseif (str .eq. 'eim') then
       ScreenI%eimax = vars
       write(3,*) 'eimax',ScreenI%eimax
    elseif (str .eq. 'erm') then
       ScreenI%ermin = vars
       write(3,*) 'ermin',ScreenI%ermin
    elseif (str .eq. 'rfm') then
       ScreenI%rfms  = vars
       write(3,*) 'rfms',ScreenI%rfms
    elseif (str .eq. 'nrp')then
       ScreenI%nrptx0  = int(vars)
       write(3,*) 'nrptx0',ScreenI%nrptx0
    else
       call wlog("Unrecognized keyword submitted to screen.inp ; aborting.")
       stop
    endif
    close(3)
    return
  end subroutine screen_inp_parse_and_write


  subroutine screen_read
    ! Reads screen.inp.  This routine is set up a little different from its brothers in the other input modules.
    ! This is to keep it compatible with situations where there either is no screen.inp file (in which case defaults are used for all variables),
    ! and with situations where screen.inp contains only the variables for which non-default values are specified.
    ! This is because I've only added mandatory screen.inp files being written by rdinp now 1-2012.  KJ
    integer i
    character*8 strs
    character*3 str
    double precision vars
    call screen_init  !KJ set defaults in case screen.inp doesn't exist!
    open (file=filename, unit=3, status='old', err=60)
    !KJ			read (3,*)  !KJ 11-2011 removing header line from screen.inp because it is incompatible with screen_inp_and_parse above.
    do i = 1, 12
       read(3,*,end=60)  strs, vars
       str = strs(1:3)
       if (str .eq. 'ner') ScreenI%ner   = nint(vars)
       if (str .eq. 'nei') ScreenI%nei   = nint(vars)
       if (str .eq. 'max') ScreenI%maxl  = nint(vars)
       if (str .eq. 'irr') ScreenI%irrh  = nint(vars)
       if (str .eq. 'ien') ScreenI%iend  = nint(vars)
       if (str .eq. 'lfx') ScreenI%lfxc  = nint(vars)
       if (str .eq. 'emi') ScreenI%emin  = vars
       if (str .eq. 'ema') ScreenI%emax  = vars
       if (str .eq. 'eim') ScreenI%eimax = vars
       if (str .eq. 'erm') ScreenI%ermin = vars
       if (str .eq. 'rfm') ScreenI%rfms  = vars
       if (str .eq. 'nrp') ScreenI%nrptx0  = nint(vars)
    end do
60  continue
    close(3)
    return
  end subroutine screen_read

  subroutine screen_init
    ScreenI%ner   = 40
    ScreenI%nei   = 20
    ScreenI%maxl  = 4
    ScreenI%irrh  = 1
    ScreenI%iend  = 0
    ScreenI%emin  = -40.0d0 !KJ This and next 3 values are in eV ; converted to Ha at a later point in the code (screen/rdgeom.f90)
    ScreenI%emax  = 0.0d0
    ScreenI%eimax = 2.0d0
    ScreenI%ermin = 0.001d0
    ScreenI%lfxc  = 0
    ScreenI%rfms  = 4.0d0
    ScreenI%nrptx0 = 251
  end subroutine screen_init

end module screen_inp
