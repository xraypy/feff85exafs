!=======================================================================
!     POTENTIAL
!=======================================================================

module potential_inp
  use dimsmod, only: nheadx, nphx, novrx
  use global_inp, only: ispec
  use atoms_inp, only : nph
  implicit none
  character(*),parameter,private :: filename='mod1.inp'

  character*80 title(nheadx)
  integer mpot, ntitle, ihole, ipr1, iafolp, iunf,             &
       nmix, nohole, jumprm, inters, nscmt, icoul, lfms1, ixc
  integer iz(0:nphx)
  !		iz(0:nphx)    - atomic number, input
  integer lmaxsc(0:nphx)
  real rfms1
  double precision gamach, rgrd, ca1, ecv, totvol
  double precision  xnatph(0:nphx), folp(0:nphx), spinph(0:nphx)
  !		xnatph(0:nphx) - given unique pot, how many atoms are there
  !                      of this type? (used for interstitial calc)
  !		folp(0:nphx) -  overlap factor for rmt calculation
  double precision  xion(0:nphx)
  !		xion(0:nphx)  - ionicity, input
  logical ExternalPot
  !     for OVERLAP option
  logical StartFromFile
  ! read potential from pot.bin file and start from there
  integer novr(0:nphx), iphovr(novrx,0:nphx), nnovr(novrx,0:nphx)
  double precision  rovr(novrx,0:nphx)
  !		novr(0:nphx) -  number of overlap shells for unique pot
  !		iphovr(novrx,0:nphx) -  unique pot for this overlap shell
  !		nnovr(novrx,0:nphx) -   number of atoms in overlap shell
  !		rovr(novrx,0:nphx)  -   r for overlap shell
  ! Added by Fer
  ! Used to correct the excitation energy for chemical shifts
  integer  ChSh_Type
  integer configtype !KJ 12-2010 : which method for choosing atomic configuration?
  double precision corval_emin  !KJ 12-2012 defines energy window for search for core-valence separation energy.

  !       criteria for self-consistency
  real*8,parameter :: tolmu = 1.D-3  ! Fermi level (Ha)
  real*8,parameter :: tolq = 1.D-3   ! net charge on atom iph (e)
  real*8,parameter :: tolqp = 2.D-4  ! partial charge (e.g. l=1) on atom iph (e)
  real*8,parameter :: tolsum = 0.05  ! total valence charge in Norman sphere compared to formal valence charge


contains

  subroutine potential_write
    integer ititle,ip,iph,iovr
    open (file=filename, unit=3, status='unknown')
    write(3,10) 'mpot, nph, ntitle, ihole, ipr1, iafolp, ixc,ispec'
    write(3,20) mpot, nph, ntitle, ihole, ipr1, iafolp, ixc, ispec
    write(3,10) 'nmix, nohole, jumprm, inters, nscmt, icoul, lfms1, iunf'
    write(3,20)  nmix, nohole, jumprm, inters, nscmt, icoul, lfms1, iunf
    do ititle = 1, ntitle
       write(3,10) title(ititle)
    enddo
    write(3,10) 'gamach, rgrd, ca1, ecv, totvol, rfms1, corval_emin'
    write(3,30)  gamach, rgrd, ca1, ecv, totvol, rfms1, corval_emin
    write(3,10) ' iz, lmaxsc, xnatph, xion, folp'
120 format ( 2i5, 4f13.5)
    do ip = 0, nph
       write(3,120) iz(ip), lmaxsc(ip), xnatph(ip), xion(ip), folp(ip)
    enddo
    write(3,10) 'ExternalPot switch, StartFromFile switch'
    write(3,*) ExternalPot,StartFromFile
    !       for OVERLAP option
    write(3,10) 'OVERLAP option: novr(iph)'
    write(3,20) ( novr(iph), iph=0,nph)
    write(3,10) ' iphovr  nnovr rovr '
140 format ( 2i5, f13.5)
    do iph = 0, nph
       do iovr = 1, novr(iph)
          write(3,140) iphovr(iovr, iph), nnovr(iovr,iph), rovr(iovr,iph)
       enddo
    enddo
    ! Added by Fer
    ! Correction of the excitation energy for chemical shifts
    write(3,10) 'ChSh_Type:'
    write(3,20) ChSh_Type
    write(3,10) 'ConfigType:'
    write(3,20) configtype
    close(3)
    ! standard formats for string, integers and real numbers
10  format(a)
20  format (20i4)
30  format (9f13.5)
  end subroutine potential_write

  subroutine potential_read
    ! integer ititle,ip,iph,iovr

    call json_read_pot(mpot, nph, ntitle, ihole, ipr1, iafolp, &
         ixc, ispec, nmix, nohole, jumprm, inters, nscmt, icoul, &
         lfms1, iunf, gamach, rgrd, ca1, ecv, totvol, rfms1, &
         title, iz, lmaxsc, xnatph, xion, folp, novr, &
         iphovr, nnovr, rovr )
    
    ! open (file=filename, unit=3, status='old')
    ! read(3,*) ; read(3,*) mpot, nph, ntitle, ihole, ipr1, iafolp, ixc, ispec
    ! read(3,*) ; read(3,*)  nmix, nohole, jumprm, inters, nscmt, icoul, lfms1, iunf
    ! do ititle = 1, ntitle
    !    read(3,*) title(ititle)
    ! enddo
    ! read(3,*) ; read(3,*)  gamach, rgrd, ca1, ecv, totvol, rfms1 !, corval_emin
    ! read(3,*)
    ! do ip = 0, nph
    !    read(3,*) iz(ip), lmaxsc(ip), xnatph(ip), xion(ip), folp(ip)
    ! enddo
    ! !read(3,*) ; read(3,*) ExternalPot, StartFromFile
    ! read(3,*) ; read(3,*) (novr(iph), iph=0,nph)
    ! read(3,*)
    ! do iph = 0, nph
    !    do iovr = 1, novr(iph)
    !       read(3,*) iphovr(iovr, iph), nnovr(iovr,iph), rovr(iovr,iph)
    !    enddo
    ! enddo
    ! !read(3,*) ; read(3,*) ChSh_Type
    ! !read(3,*,end=55) ; read(3,*,end=55) configtype
    ! ! 55  continue
    ! close(3)
  end subroutine potential_read

  subroutine potential_init
    title(:) = ' '
    mpot = 1
    nph = 0
    ntitle = 0
    ihole = 1
    ipr1 = 0
    iafolp = 0
    iunf = 0
    nmix = 1
    nohole = -1
    jumprm = 0
    inters = 0
    nscmt = 0
    icoul = 0
    ixc = 0
    lfms1 = 0
    iz(:) = -1
    lmaxsc(:) = 0
    rfms1 = -1 * 1.e0
    ca1 = 0.d0
    ecv = -40*1.d0
    rgrd = 0.05 * 1.d0
    totvol = 0.d0
    gamach = 0.d0 !initialized later by setgam
    xnatph(:) = 0.d0
    spinph(:) = 0.d0
    xion(:) = 0.d0
    folp(:) = 1.d0
    ExternalPot = .false.
    StartFromFile = .false. !KJ added 12-10
    novr(:) = 0
    iphovr(:,:)=0 !KJ added 7-09
    nnovr(:,:)=0 !KJ
    rovr(:,:) = 0.d0 !KJ
    ChSh_Type = 0 !Fer : standard feff
    configtype=1 !KJ 12-2010 standard feff9
    corval_emin=-70.d0 ! eV
  end subroutine potential_init

end module potential_inp
