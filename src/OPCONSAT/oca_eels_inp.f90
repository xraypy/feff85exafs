!=======================================================================
!     EELS
!=======================================================================

module eels_inp
  !		Beam direction in crystal frame of feff.inp (a.u.)
  use global_inp,only: xivec
  implicit none
  character(*),parameter,private :: filename='eels.inp'
  !		Beam energy in eV :
  real*8 ebeam
  !		Convergence semiangle in rad :
  real*8 aconv
  !		Collection semiangle in rad :
  real*8 acoll
  !		Integration mesh for q-vectors (radial/angular mesh size)
  integer nqr,nqf
  !		Detector position ; angles in rad w.r.t. x and y directions
  real*8 thetax,thetay
  !       what kind of q-mesh : uniform (U), logarithmic (L), or one dimensional logarithmic (1)
  !       not currently in eels.inp/feff.inp
  character*1      qmodus
  !		Parameter for logarithmic mesh - not currently in eels.inp/feff.inp
  real*8 th0
  !		Make magic angle plot if magic=1
  integer        magic
  !		Evaluate magic angle at this energy point
  real*8        emagic
  !		Orientation sensitive?
  integer        aver
  !		Do we have cross-terms?
  integer        cross
  !		Do we do anything at all?
  integer        eels
  !		How many spectra to combine
  integer ipmin,ipmax,ipstep ,nip
  !		Where do we take input from :
  integer iinput           
  !		Which column? - to be replaced by more advanced switch
  integer spcol
  !       Relativistic calculation or not?  Converted into logical inside eels-module.
  integer relat

contains

  subroutine eels_write
    open (file=filename, unit=3, status='unknown')
    write(3,10) 'calculate ELNES?'
    write(3,20) eels
    write(3,10) 'average? relativistic? cross-terms? Which input?'
    write(3,20) aver, relat, cross, iinput, spcol
    write(3,10) 'polarizations to be used ; min step max'
    write(3,20) ipmin,ipstep,ipmax
    write(3,10) 'beam energy in eV'
    write(3,30) ebeam
    write(3,10) 'beam direction in arbitrary units'
    write(3,30) xivec
    write(3,10) 'collection and convergence semiangle in rad'
    write(3,30) acoll,aconv
    write(3,10) 'qmesh - radial and angular grid size'
    write(3,20) nqr,nqf
    write(3,10) 'detector positions - two angles in rad'
    write(3,30) thetax,thetay
    write(3,10) 'calculate magic angle if magic=1'
    write(3,20) magic
    write(3,10) 'energy for magic angle - eV above threshold'
    write(3,30) emagic
    close(3)
    ! standard formats for string, integers and real numbers
10  format(a)
20  format (20i4)
30  format (9f13.5)
  end subroutine eels_write

  subroutine eels_read
    open (file=filename, unit=3, status='old',err=100)
    read(3,*) ; read(3,*,end=100,err=100) eels
    read(3,*) ; read(3,*,err=209) aver, relat, cross, iinput,spcol ; goto 210
209 iinput=1;spcol=4;relat=1;cross=1;aver=0  !restore defaults - this construction for older files.
210 read(3,*) ; read(3,*) ipmin,ipstep,ipmax  !KJ this un-Kevin-like construction for older files ...
    nip=1+((ipmax-ipmin)/ipstep)
    read(3,*) ; read(3,*) ebeam
    read(3,*) ; read(3,*) xivec
    read(3,*) ; read(3,*) acoll,aconv
    read(3,*) ; read(3,*) nqr,nqf
    read(3,*) ; read(3,*) thetax,thetay
    read(3,*) ; read(3,*) magic
    read(3,*) ; read(3,*) emagic
    close(3)
    return
100 eels = 0  ; ipmin=1 ; ipmax=1 ; ipstep=1 ! no eels.inp -> don't do eels
    return
  end subroutine eels_read


  subroutine eels_init !default values for everything (except xivec)
    ebeam  = 0.d0
    aconv  = 0.d0
    acoll  = 0.d0
    nqr    = 0
    nqf    = 0
    magic  = 0
    emagic = 0.d0
    eels   = 0
    relat  = 1
    cross  = 1
    aver   = 0
    thetax = 0.d0
    thetay = 0.d0
    ipmin  = 1
    ipmax  = 1
    nip    = 1
    ipstep = 1
    iinput = 1  ! xmu.dat - files from ff2x
    spcol  = 4   ! xmu.dat - use spectrum mu(omega)
    qmodus = 'U'  !  U for uniform grid 
    th0    = 0.d0
  end subroutine eels_init

end module eels_inp
