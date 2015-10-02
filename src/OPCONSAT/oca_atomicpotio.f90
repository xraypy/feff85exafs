!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: m_atomicpotio.f90,v $:
! $Revision: 1.16 $
! $Author: jorissen $
! $Date: 2010/12/14 00:22:37 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE AtomicPotIO
  USE ErrorMod
  USE IOMod
  USE dimsmod
  IMPLICIT NONE

  CHARACTER(3),PRIVATE :: FileType
  !PARAMETER(FileType = 'PAD')
  PARAMETER(FileType = 'TXT')


CONTAINS
  SUBROUTINE WriteAtomicPots(nph, iz, ihole, rho, dmag, rhoval, vcoul, dgc0, dpc0, dgc, dpc, &
       & adgc, adpc, erelax, emu, xnvmu, xnval, norb, eorb, drho, dvcoul, iphat,    &
       & rat, iatph, novr, iphovr, nnovr, rovr, nat, edens, edenvl, vclap, rnrm,    &
       & kappa, iorb, s02)

    ! Scalar data
    INTEGER,INTENT(IN) :: nph, nat, ihole
    DOUBLE PRECISION, INTENT(IN) :: emu, erelax, s02

    ! 1D data
    INTEGER,INTENT(IN) :: iz(:), norb(:), iphat(:), iatph(:), novr(:)
    DOUBLE PRECISION,INTENT(IN) :: dgc0(:), dpc0(:), drho(:), dvcoul(:), rnrm(:)

    ! 2D data
    INTEGER,INTENT(IN) :: iphovr(:,:), nnovr(:,:), kappa(:,:), iorb(:,:)
    DOUBLE PRECISION,INTENT(IN) :: rho(:,:), dmag(:,:), rhoval(:,:), vcoul(:,:), &
         & xnvmu(:,:), xnval(:,:), eorb(:,:), rat(:,:), rovr(:,:), edens(:,:),   &
         & edenvl(:,:), vclap(:,:)

    ! 3D data
    DOUBLE PRECISION,INTENT(IN) :: dgc(:,:,:), dpc(:,:,:), adgc(:,:,:), adpc(:,:,:)




    CHARACTER*80 Headers(20)
    CHARACTER*20 ColumnLabels(20)

    INTEGER iph

    ! Initialize stuff.
    Headers(:) = ' '
    ColumnLabels(:) = ' '

    ! Define some headers.
    Headers(1) = 'This file contains information about the free atom potentials.'
    Headers(2) = 'nat    - number of atoms.'
    Headers(3) = 'nph    - number of unique potentials'      
    Headers(3) = 'ihole  - hole index'      
    Headers(4) = 'erelax - relaxation energy for each occupied orbital.'
    Headers(5) = 'emu - edge energy for each occupied orbital.'
    Headers(6) = 's02 - many body amplitude reduction.'
    ColumnLabels(1) = 'nph'
    ColumnLabels(2) = 'nat'
    ColumnLabels(3) = 'ihole'
    ColumnLabels(4) = 'erelax'
    ColumnLabels(5) = 'emu'
    ColumnLabels(6) = 's02'
    ! Write headers along with the first line and columnlabels.
    ! Write scalar data.
    CALL WriteData('apot.bin', Int1 = nph, Int2 = nat, Int3 = ihole, Double4 = erelax, Double5 = emu, &
         & Double6 = s02, Headers = Headers, ColumnLabels = ColumnLabels, FileType = FileType, &
         ForceNewSection = .TRUE. ) !KJ 7-09 added ForceNewSection bc otherwise doesn't work on my Windows pc
    Headers(:) = ' '
    ColumnLabels(:) = ' '

    ! Write iz, iatph, novr, rnrm.
    Headers(1) = 'iz(0:nphx)    - atomic number for each unique potential'
    Headers(2) = 'iatph(0:nphx) - given unique pot, which atom is model?'
    Headers(3) = 'novr(0:nphx)  - number of overlap shells for each unique pot'
    Headers(4) = 'rnrm(0:nphx)  - norman radius for each unique potential.'
    ColumnLabels(1) = 'iz'
    ColumnLabels(2) = 'iatph'
    ColumnLabels(3) = 'novr'
    ColumnLabels(4) = 'rnrm'
    CALL WriteArrayData('apot.bin', Int1 = iz, Int2 = iatph, Int3 = novr, Double4 = rnrm, &
         & Headers = Headers, ColumnLabels = ColumnLabels, FileType = FileType,           &
         & ForceNewSection = .TRUE.)
    Headers(:) = ' '
    ColumnLabels(:) = ' '

    ! Write norb
    Headers(1) = 'norb(0:nphx+1) - number of occupied orbitals for each unique potential'
    CALL WriteArrayData('apot.bin', Int1 = norb, Headers = Headers, &
         & FileType = FileType, ForceNewSection = .TRUE.)

    ! Write iphat
    Headers(1) = 'iphat(natx) - given specific atom, which unique pot?'
    CALL WriteArrayData('apot.bin', Int1 = iphat, Headers = Headers, &
         & FileType = FileType, ForceNewSection = .TRUE.)

    ! Write dgc0, dpc0, drho, and dvcoul
    Headers(1) = 'dgc0   - upper component of core hole orbital'
    Headers(2) = 'dpc0   - lower component of core hole orbital'
    Headers(3) = 'drho   - core hole density.'
    Headers(4) = 'dvcoul - core hole coulomb potential.'
    ColumnLabels(1) = 'dgc0'
    ColumnLabels(2) = 'dpc0'
    ColumnLabels(3) = 'drho'
    ColumnLabels(4) = 'dvcoul'
    CALL WriteArrayData('apot.bin', Double1 = dgc0, Double2 = dpc0, Double3 = drho, &
         & Double4 = dvcoul, Headers = Headers, ColumnLabels = ColumnLabels,          &
         & FileType = FileType, ForceNewSection = .TRUE.)
    Headers(:) = ' '
    ColumnLabels(:) = ' '

    ! Write iphovr.
    Headers(1) = 'iphovr(novrx,0:nphx) - unique pot for each overlap shell'
    CALL Write2D('apot.bin', iphovr, Headers = Headers, FileType = FileType)

    ! Write nnovr.
    Headers(1) = 'nnovr(novrx,0:nphx) - number of atoms in overlap shell'
    CALL Write2D('apot.bin', nnovr, Headers = Headers, FileType = FileType)

    ! Write rho.
    Headers(1) = 'rho(r,0:nphx+1) - atomic density for each unique potential'
    Headers(2) = '                  nph+1 holds final state potential for absorber'
    CALL Write2D('apot.bin', rho, Headers = Headers, FileType = FileType)
    Headers(:) = ' '

    ! Write dmag.
    Headers(1) = 'dmag(r,nph+1) - ?'
    CALL Write2D('apot.bin', dmag, Headers = Headers, FileType = FileType)

    ! Write rhoval.
    Headers(1) = 'rhoval(r,nph+1) - atomic valence density for each unique potential.'
    CALL Write2D('apot.bin', rhoval, Headers = Headers, FileType = FileType)

    ! Write vcoul.
    Headers(1) = 'vcoul(r,nph) - coulomb potential for each unique potential.'
    CALL Write2D('apot.bin', vcoul, Headers = Headers, FileType = FileType)

    ! Write xnvmu.
    Headers(1) = 'xnvmu(0:lx,0:nphx+1) - number of val. el, w/in norman sphere for each channel'
    CALL Write2D('apot.bin', xnvmu, Headers = Headers, FileType = FileType)

    ! Write xnval.
    Headers(1) = 'xnval(30,nph) - occupation of each orbital.'
    CALL Write2D('apot.bin', xnval, Headers = Headers, FileType = FileType)

    ! Write eorb.
    Headers(1) = 'eorb(norb,iph)'
    CALL Write2D('apot.bin', eorb, Headers = Headers, FileType = FileType)

    ! Write rat.
    Headers(1) = 'rat(3,nat) - cartesian coordinates for each atom.'
    CALL Write2D('apot.bin', rat, Headers = Headers, FileType = FileType)

    ! Write rovr
    Headers(1) = 'rovr(novrx,0:nphx) - r for overlap shell'
    CALL Write2D('apot.bin', rovr, Headers = Headers, FileType = FileType)

    ! Write edens.
    Headers(1) = 'edens(r,0:nphx) - overlapped density for each unique potential.'
    CALL Write2D('apot.bin', edens, Headers = Headers, FileType = FileType)

    ! Write edenvl
    Headers(1) = 'edenvl(r,0:nphx) - overlapped density of valence electrons?'
    CALL Write2D('apot.bin', edenvl, Headers = Headers, FileType = FileType)


    ! Write vclap.
    Headers(1) = 'vclap(r,0:nphx) - overlapped coulomb potential.'
    CALL Write2D('apot.bin', vclap, Headers = Headers, FileType = FileType)

    ! Write kappa.
    Headers(1) = 'kappa(norb,0:nph) - quntum number kappa for each orbital and potential.'
    CALL Write2D('apot.bin', kappa, Headers = Headers, FileType = FileType)

    ! Write iorb.
    Headers(1) = 'iorb(-4:3,0:nphx+1) - last occupied orbital of a particular kappa or 0 if none.'
    CALL Write2D('apot.bin', iorb, Headers = Headers, FileType = FileType)

    ! Write dgc.
    Headers(1) = 'dgc(r,30,nph) - upper component of each obital for each unique potential.'
    CALL WriteData('apot.bin',Headers = Headers)
    DO iph = 1, SIZE(dgc,3)
       IF(iph.eq.1) THEN
          Headers(1) = 'dgc(r,30,nph) - upper component of each obital for each unique potential.'
          Headers(2) = 'dgc(r,norb,' // ACHAR(iph+48) // ')'
       ELSEIF(iph.lt.10) THEN
          Headers(1) = 'dgc(r,norb,' // ACHAR(iph+48) // ')'
       ELSEIF(iph.lt.100) THEN
          Headers(1) = 'dgc(r,norb,' // ACHAR((iph/10)+48) // ACHAR(mod(iph,10)+48) // ')'
       ELSE
          stop 'ERROR iph too large in WriteAtomicPots'
       END IF
       CALL Write2D('apot.bin', dgc(:,:,iph), Headers = Headers, FileType = FileType)
    END DO

    ! Write dpc.
    Headers(1) = 'dpc(r,30,nph) - lower component of each obital and unique potential.'
    CALL WriteData('apot.bin',Headers = Headers)
    DO iph = 1, SIZE(dpc,3)
       IF(iph.eq.1) THEN
          Headers(1) = 'dpc(r,30,nph) - lower component of each obital and unique potential.'
          Headers(2) = 'dpc(r,norb,' // ACHAR(iph+48) // ')'
       ELSEIF(iph.lt.10) THEN
          Headers(1) = 'dpc(r,norb,' // ACHAR(iph+48) // ')'
       ELSEIF(iph.lt.100) THEN
          Headers(1) = 'dpc(r,norb,' // ACHAR((iph/10)+48) // ACHAR(mod(iph,10)+48) // ')'
       ELSE
          stop 'ERROR iph too large in WriteAtomicPots'
       END IF
       CALL Write2D('apot.bin', dpc(:,:,iph), Headers = Headers, FileType = FileType)
    END DO

    ! Write adgc.
    DO iph = 1, SIZE(adgc,3)
       IF(iph.eq.1) THEN
          Headers(1) = 'adgc(r,30,nph) - upper devel. coeficients for each obital and unique potential.'
          WRITE(Headers(2),'(A,I10,A)') 'adgc(r,norb,', iph, ')'
       ELSE
          WRITE(Headers(1),'(A,I10,A)') 'adgc(r,norb,', iph, ')'
       END IF
       CALL Write2D('apot.bin', adgc(:,:,iph), Headers = Headers, FileType = FileType)
    END DO

    ! Write adpc.
    Headers(1) = 'adpc(r,30,nph) - lower devel. coeficients for each obital and unique potential.'
    CALL WriteData('apot.bin',Headers = Headers)
    DO iph = 1, SIZE(adpc,3)
       IF(iph.eq.1) THEN
          Headers(1) = 'adpc(r,30,nph) - lower devel. coeficients for each obital and unique potential.'
          WRITE(Headers(2),'(A,I10,A)') 'adpc(r,norb,', iph, ')'
       ELSE
          WRITE(Headers(1),'(A,I10,A)') 'adpc(r,norb,', iph, ')'
       END IF
       CALL Write2D('apot.bin', adpc(:,:,iph), Headers = Headers, FileType = FileType)
    END DO

    ! Close the file so that we can read it later.
    CALL CloseFl('apot.bin')
  END SUBROUTINE WriteAtomicPots


  SUBROUTINE ReadAtomicPots(nph, rnrm)
    INTEGER nph
    INTEGER imt(0:nphx), rmt(0:nphx), inrm(0:nphx)
    REAL(8) folp(0:nphx), folpx(0:nphx), dgc0(251), dpc0(251)
    REAL(8) dgc(251, 30, 0:nphx), dpc(251, 30, 0:nphx)
    REAL(8) adgc(10, 30, 0:nphx), adpc(10, 30, 0:nphx)
    REAL(8) edens(251, 0:nphx), vclap(251, 0:nphx)
    REAL(8) vtot(251, 0:nphx), edenvl(251, 0:nphx)
    REAL(8) vvalgs(251, 0:nphx), dmag(251, 0:nphx)
    REAL(8) xnval(30,0:nphx), qnrm(0:nphx), xnmues(0:lx,0:nphx)
    REAL(8) eorb(30), kappa(30), rnrm(0:nphx)
    INTEGER iorb(-4:3,0:nphx), iz(0:nphx), xion(0:nphx), iafolp, ihole, &
         &        inters, iunf, jumprm, nohole, ntitle
    REAL(8) xnatph(0:nphx), ecv, emu, erelax, qtotel, rhoint, rnrmav,   &
         &        rs, s02, totvol, vint, wp, xf, xmu

    character*80 title(nheadx)

    CALL rdpot ( ntitle, title, rnrmav, xmu, vint, rhoint, &
         emu, s02, erelax, wp, ecv,rs,xf, qtotel,          &
         imt, rmt, inrm, rnrm, folp, folpx, xnatph,        &
         dgc0, dpc0, dgc, dpc, adgc, adpc,                 &
         edens, vclap, vtot, edenvl, vvalgs, dmag, xnval,  &
         eorb, kappa, iorb, qnrm, xnmues, nohole, ihole,   &
         inters, totvol, iafolp, xion, iunf, iz, jumprm)
  END SUBROUTINE ReadAtomicPots

  ! Hack for Fer's mt potentials
  SUBROUTINE WriteExternalPot(vtot, vint, edens, rhoint, rat, xmu, imt, rmt)
    DOUBLE PRECISION vtot(:,:), vint, edens(:,:), rhoint, rat(:,:), xmu, rmt(:), RadialGrid(251), xNElOut
    DOUBLE PRECISION,ALLOCATABLE :: ratTranspose(:,:)
    INTEGER nat, NRPts, i1, i2, imt(:), ncoord
    INTEGER,ALLOCATABLE :: iat(:)
    CHARACTER(30) FileName

    ncoord = 3
    FileName = 'extpot.aip'
    nat = 2
    ! Read number of atoms
    CALL WriteData(FileName, Int1 = nat, ForceNewSection = .TRUE.)

    ! Read atomic numbers
    ALLOCATE(iat(nat))
    iat(1) = 6
    iat(2) = 8
    CALL WriteArrayData(FileName, Int1 = iat, ForceNewSection = .TRUE.)

    ! Read coordinates
    ALLOCATE(ratTranspose(nat,3))
    DO i1 = 1, nat
       DO i2 = 1, 3
          ratTranspose(i1,i2) = rat(i2, i1)
       END DO
    END DO
    CALL Write2D(FileName, ratTranspose)

    ! Read number of radial points
    NRPts = 210
    CALL WriteData(FileName, Int1 = NRPts, ForceNewSection = .TRUE.)

    ! Read Interstitial potential
    CALL WriteData(FileName, Double1 = vint, ForceNewSection = .TRUE.)

    ! Read Fermi energy
    CALL WriteData(FileName, Double1 = xmu, ForceNewSection = .TRUE.)

    ! Read number of electrons outside muffin tins.
    xNElOut = 1
    CALL WriteData(FileName, Double1 = xNElOut, ForceNewSection = .TRUE.)

    ! Read muffin tin radii and jri.
    CALL WriteArrayData(FileName, Int1 = imt, Double2 = rmt, ForceNewSection = .TRUE.)

    ! Read radial grid, vtot, and edens
    DO i1 = 1, nat
       DO i2 = 1, NRPts
          CALL WriteData(FileName, Double1 = RadialGrid(i2), Double2 = vtot(i2,i1), Double3 = edens(i2,i1))
       END DO
    END DO

    CALL CloseFl(FileName)
  END SUBROUTINE WriteExternalPot


  ! Hack for Fer's mt potentials
  SUBROUTINE ReadExternalPot(vtot, vint, edens, rhoint, rat, xmu, imt, rmt)
    USE Mtdp
    USE IOFiles

    DOUBLE PRECISION vtot(:,:), vint, edens(:,:), rhoint, rat(:,:), xmu, rmt(:) ! , RadialGrid(251), vTmp , xNElOut
    DOUBLE PRECISION :: rmt2(SIZE(rmt)), vtot2(SIZE(vtot,1),SIZE(vtot,2)), edens2(SIZE(edens,1),SIZE(edens,2)) 
    INTEGER nat, NRPts, i1, i2, imt(:), ncoord, iSort(200), iunit ! , iError
    INTEGER :: imt2(SIZE(imt)) ! , iz(SIZE(imt))
    CHARACTER(30) PotFile, SortFile, mtdpFile
    TYPE(Mtdp_Data_Type) :: Mtdp_Data
    ncoord = 3
    PotFile  = 'extpot.aip'
    SortFile = 'sort.aip'
    mtdpFile = 'GeCl4.04.dft.mtdp'

    ! Read mtdp file.
    CALL OpenFl(mtdpFile,FileStatus='OLD',FileAction='READ')
    CALL GetIOFileInfo(mtdpFile, UnitNumber = iunit)
    CALL Read_Mtdp(iunit,Mtdp_Data)

    ! Number of points in the radial grid
    NRPts = Mtdp_Data%nR

    ! Number of atoms
    nat = Mtdp_Data%nAt

    ! For now atomic numbers are not needed.

    ! Atomic coordinates - ignore for now
    !      rat(:,:) = Mtdp_Data%At_XYZ

    ! Muffin-Tin Radii
    rmt2(:nat) = Mtdp_Data%At_R(:)

    ! Index of the radii on the radial grid
    imt2(:nat) = Mtdp_Data%At_iR(:)

    ! Electron density inside each Muffin-Tin
    edens2(:,:nat) = Mtdp_Data%At_Den(:,:)

    ! Potential (H+XC) inside each Muffin-Tin
    vtot2(:,:nat) = Mtdp_Data%At_Pot(:,:)

    ! Empty spheres

    ! muffin tin radii
    rmt2(nat+1:Mtdp_Data%nESph+nat) = Mtdp_Data%ESph_R(:)

    ! muffin tin index
    imt2(nat+1:Mtdp_Data%nESph+nat) = Mtdp_Data%ESph_iR(:)

    ! Electron density inside each empty sphere
    edens2(:,nat+1:Mtdp_Data%nESph+nat) = Mtdp_Data%ESph_Den(:,:)

    ! Potential (H+XC) inside each empty sphere
    vtot2(:,nat+1:Mtdp_Data%nESph+nat) = Mtdp_Data%ESph_Pot(:,:)

    PRINT*, 'Old vint = ', vint
    ! Interstitial potential
    vint = Mtdp_Data%V_Int
    PRINT*, 'Enter vint: '
    READ*, vint    

    ! HOMO energy
    xmu = (Mtdp_Data%V_HOMO + Mtdp_Data%V_LUMO)/2.d0
    !      xmu = Mtdp_Data%V_HOMO
    nat = nat + Mtdp_Data%nESph
    ! Read isort out of sort.aip. This tells which potentials go with which atoms in feff.
    iSort(:) = -1
    ! PRINT*, 'Reading iSort'
    CALL ReadArrayData(SortFile, Int1 = iSort)
    ! iSort(i1) is the unique potential given the i1th potential defined in the extpot.aip file.
    ! iSort can be zero, but arrays here are defined from 1, so we need to shift.
    DO i1 = 1, SIZE(iSort)
       iSort(i1) = iSort(i1) + 1
    END DO

    ! Fill vtot, edens, rmt, and imt with the proper values based on iSort
    DO i1 = 1, SIZE(iSort)
       IF(iSort(i1).gt.0) THEN
          ! Fill rmt and imt
          IF(iSort(i1).gt.nat) CALL Error('ERROR: Number of potentials defined '// &
               & 'in sort.aip is greater than number defined in extpot.aip.')
          rmt(iSort(i1)) = rmt2(i1)
          imt(iSort(i1)) = imt2(i1)

          ! Fill vtot and edens
          DO i2 = 1, NRPts
             vtot(i2,iSort(i1)) = vtot2(i2,i1)
             edens(i2,iSort(i1)) = edens2(i2,i1)
          END DO
       END IF
    END DO

    ! If vint is too close to zero, reset it.
    !      IF(vint.gt.-0.1d0) vint = -0.1d0
    DO i1 = 1, nat
       DO i2 = NRPts + 1, 251
          vtot(i2,i1) = vint
          edens(i2,i1) = rhoint
       END DO
    END DO

    !      DO i1 = 1, 251
    !         vtot(i2,nat+1) = vtot(NRPts,1)
    !         edens(i2,nat+1) = edens(NRPts,1)
    !      END DO

    CALL CloseFl(SortFile)
    CALL CloseFl(mtdpFile)
    ! PRINT*, 'Done reading external potential.'
  END SUBROUTINE ReadExternalPot

  !     SUBROUTINE SplitMtDP(PotFiles,mtdpFile)
  !       USE Mtdp
  !       USE IOFiles
  !       ! This subroutine reads an mtdp file and splits it into nph pot files.
  !       ! Input is the array of potential file names.
  !       CHARACTER PotFiles(:)

  !       INTEGER iunit, iat
  !       CHARACTER AtomLabel

  !       CALL OpenFl(mtdpFile)
  !       CALL GetIOFileInfo(UnitNumber = iunit)
  !       CALL Read_Mtdp(iunit,Mtdp_Data)

  !       DO iat = 1, Mtdp_Data%nAt
  !          ! Number of points in the radial grid
  !          Mtdp_Data%nR

  !          CALL WriteData



  SUBROUTINE ReadExternalPotWien2k   ! (vtot, vint, edens, rhoint, rat, xmu, imt, rmt)

  END SUBROUTINE ReadExternalPotWien2k



END MODULE AtomicPotIO
