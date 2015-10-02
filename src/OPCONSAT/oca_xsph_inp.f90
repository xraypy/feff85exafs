!=======================================================================
!     XSPH
!=======================================================================

module xsph_inp
  use dimsmod
  use global_inp
  use potential_inp
  use ldos_inp
  implicit none
  character(*),parameter,private :: filename='xsph.inp'
  integer mphase, ipr2, ixc0, lreal, iPlsmn
  integer iGammaCH, iGrid, NPoles
  character*6  potlbl(0:nphx)
  !		potlbl(0:nphx)    -   label for user convienence
  double precision xkstep, xkmax, vixan, vr0, vi0, Eps0, EGap
  integer izstd, ifxc, ipmbse, itdlda, nonlocal, ibasis
  !		!KJ for the energy grid card EGRID :
  integer iegrid,egrid3a
  real*8 egrid3b,egrid3c
  character*100 egridfile
  logical lopt

contains

  subroutine xsph_write
    integer iph
    open (file=filename, unit=3, status='unknown')
    !     Josh - added flag for PLASMON card (iPlsmn = 0, 1, or 2)
    write(3,10) 'mphase,ipr2,ixc,ixc0,ispec,lreal,lfms2,nph,l2lp,iPlsmn,NPoles,iGammaCH,iGrid'
    write(3,20)  mphase,ipr2,ixc,ixc0,ispec,lreal,lfms2,nph,l2lp,   &
         &        iPlsmn, NPoles, iGammaCH, iGrid
    write(3,10) 'vr0, vi0'
    write(3,30)  vr0, vi0
    write(3,10) ' lmaxph(0:nph)'
    write(3,20)  (lmaxph(iph),iph=0,nph)
    write(3,10) ' potlbl(iph)'
    write(3,170)  (potlbl(iph),iph=0,nph)
170 format (13a6)
    write(3,10) 'rgrd, rfms2, gamach, xkstep, xkmax, vixan, Eps0, EGap'
    write(3,30)  rgrd, rfms2, gamach, xkstep, xkmax, vixan, Eps0, EGap
    write(3,30)  (spinph(iph),iph=0,nph)
    write(3,20)  izstd, ifxc, ipmbse, itdlda, nonlocal, ibasis
    ! Commented out by Fer
    ! The following lines are commented out because they are not being read
    ! in rexsph (commented out by JK). This screws up anything that comes after
    ! them in mod2.inp (for example, the ChSh parameters that I'm including.
    !!KJ next lines contain EGRID variables ; added 01-07
    !        write(3,10) 'iegrid,egrid3a,egrid3b,egrid3c'
    !          write(3,'(2i4,2f13.5)') iegrid,egrid3a,egrid3b,egrid3c !format statement is a mix of 20 and 30
    !          write(3,10) 'egridfile'
    !          write(3,10) egridfile
    !!KJ
    ! Added by Fer
    ! Correction of the excitation energy for chemical shifts
    write(3,10) 'ChSh_Type:'
    write(3,20) ChSh_Type
    !KJ 7-09 Next 2 lines for feff8q
    write(3,'(a)') ' the number of decomposition channels ; only used for nrixs'
    write(3,'(i5)') ldecmx
    write(3,'(a)') 'lopt'
    write(3,*) lopt
    close(3)
    ! standard formats for string, integers and real numbers
10  format(a)
20  format (20i4)
30  format (20f13.5)
  end subroutine xsph_write

  subroutine xsph_read
    integer iph
    open (file=filename, unit=3, status='old')
    read(3,*) ; read(3,*)  mphase,ipr2,ixc,ixc0,ispec,lreal,lfms2,nph,l2lp,iPlsmn, NPoles, iGammaCH, iGrid
    read(3,*) ; read(3,*)  vr0, vi0
    read(3,*) ; read(3,*)  (lmaxph(iph),iph=0,nph)
    read(3,*) ; read(3,'(13a6)')  (potlbl(iph),iph=0,nph)
    read(3,*) ; read(3,*)  rgrd, rfms2, gamach, xkstep, xkmax, vixan, Eps0, EGap
    read(3,*)  (spinph(iph),iph=0,nph)
    read(3,*)  izstd, ifxc, ipmbse, itdlda, nonlocal, ibasis
    !!KJ next lines contain EGRID variables ; added 01-07
    !          read(3,*) ; read(3,'(2i4,2f13.5)') iegrid,egrid3a,egrid3b,egrid3c !format statement is a mix of 20 and 30
    !          read(3,*) ; read(3,10) egridfile
    read(3,*) ; read(3,*) ChSh_Type
    !KJ 7-09 Next 2 lines for feff8q
    read(3,*) ; read(3,*) ldecmx
    read(3,*) ; read(3,*) lopt
    close(3)
  end subroutine xsph_read

  subroutine xsph_init
    lopt = .false.
    izstd = 0
    ifxc = 0
    ipmbse = 0
    itdlda = 0
    nonlocal = 0
    ibasis = 0
    potlbl(0:nphx) = ' '
    mphase = 1
    ipr2 = 0
    ixc0 = -1
    lreal = 0
    iPlsmn = 0 ! Josh Kas
    NPoles = 100 ! JJK 3/9/2010
    EGap = 0.d0 ! JJK 4/2010
    iGammaCH = 0
    iGrid = 0
    vr0 = 0.d0
    vi0 = 0.d0
    xkmax = 20*1.d0
    xkstep = 0.07*1.d0
    vixan = 0.d0
    iegrid=0 !KJ for EGRID card 1-07
    egridfile=' '
    egrid3a=0
    egrid3b=dble(0)
    egrid3c=dble(0)
  end subroutine xsph_init


end module xsph_inp
