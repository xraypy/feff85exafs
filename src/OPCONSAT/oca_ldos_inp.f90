!=======================================================================
!     LDOS
!=======================================================================

module ldos_inp
  use atoms_inp,only: nph
  use potential_inp,only: ixc, rgrd
  use global_inp,only: ispin
  use dimsmod,only : nphx
  implicit none
  character(*),parameter,private :: filename='ldos.inp'
  integer mldos, lfms2, minv, lmaxph(0:nphx)
  double precision emin, emax, eimag
  integer neldos
  real rdirec, toler1, toler2, rfms2
  logical save_g0, save_compton_info ! BAM 2/2012

contains

  subroutine ldos_write
    integer iph
    open (file=filename, unit=3, status='unknown')
    write(3,10) 'mldos, lfms2, ixc, ispin, minv, neldos'
    write(3,20)  mldos, lfms2, ixc, ispin, minv, neldos
    write(3,10) 'rfms2, emin, emax, eimag, rgrd'
    write(3,30)  rfms2, emin, emax, eimag, rgrd
    write(3,10) 'rdirec, toler1, toler2'
    write(3,30)  rdirec, toler1, toler2
    write(3,10) ' lmaxph(0:nph)'
    write(3,20)  (lmaxph(iph),iph=0,nph)
    write(3,10) 'save_g0? save_compton_info?'
    write(3,*)  save_g0, save_compton_info
    close(3)
    ! standard formats for string, integers and real numbers
10  format(a)
20  format (20i4)
30  format (9f13.5)
  end subroutine ldos_write

  subroutine ldos_read
    integer iph
    open (file=filename, unit=3, status='old')
    read(3,*) ; read(3,*)  mldos, lfms2, ixc, ispin, minv, neldos
    read(3,*) ; read(3,*)  rfms2, emin, emax, eimag, rgrd
    read(3,*) ; read(3,*)  rdirec, toler1, toler2
    read(3,*) ; read(3,*)  (lmaxph(iph),iph=0,nph)
    read(3,*) ; read(3,*)  save_g0, save_compton_info
    close(3)
  end subroutine ldos_read

  subroutine ldos_init
    mldos = 0
    lfms2 = 0
    minv = 0
    emax = 0.d0
    emin = 1000*1.d0
    eimag = -1*1.d0
    neldos = 101
    rfms2 = -1 * 1.e0
    rdirec = -1 * 1.e0
    toler1 = 1.e-3
    toler2 = 1.e-3
    lmaxph(:) = 0
    save_g0 = .false.
    save_compton_info = .false.
  end subroutine ldos_init

end module ldos_inp
