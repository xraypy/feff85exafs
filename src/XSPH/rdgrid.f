      SUBROUTINE RdGrid(em,ne,nGrid,iGridType,GridMin,GridMax,GridStep,
     &     nGridMax,nex)
!     Read data from grid.inp.
!     file should have lines with the following format:
!     
!     Grid_Type    GridMin    GridMax    GridStep
!     
!     Grid_Type can be any of the following (case insensitive):
!     egrid    - (regular in energy)
!     kgrid    - (regular in k)
!     expgrid  - (exponential in energy)
!     usergrid - (read a grid from the file)
!
!     Energy, and k are given relative to the edge.
!     For egrid and expgrid, GridMin, GridMax, and GridStep are
!     given in (eV).
!     For kgrid, units are inverse angstroms.
!     usergrid is a special case and is followed by one energy
!     point per line, i.e.
!            usergrid
!            -1.01
!            -0.55
!            10.01
!              .!              .
!              .
!
!     More than one grid may be specified, and grids can
!     overlap, for example:
!
!     egrid -10 10 0.1
!     kgrid  0  15 0.5
!
!     will make overlapping grids going from E = -10 eV to
!     k = 15 Angstrom**(-1). Up to 10 different grids can be
!     defined.
!     If the 'last' keyword is used in the GridMin field, i.e.
!        expgrid  last  100
!     the specified grid will start where the last grid ended.
!     This is usefull when defining non-overlapping k/e grids.
!     Comments lines have #,!,c, or * at the beginning.
      include '../HEADERS/const.h'
!     Input:
!     nGridMax - max number of grids that can be defined.
!     nex      - max number of energy points
      INTEGER nGridMax, nex

!     Output:
!     nGrid     - number of grids defined in file
!     iGridType - Type of grid. (0 = user, 1 = energy, 2 = k, 3 = exponential)
!     GridMin   - Minimum value of grid.
!     GridMax   - Maximum value of grid.
!     GridStep  - Step size.
!     ne        - number of energy points
!     em(nex)        - energy grid
      INTEGER nGrid, iGridType(nGridMax), ne
      DOUBLE PRECISION GridMin(nGridMax), GridMax(nGridMax),
     &     GridStep(nGridMax)
      COMPLEX*16 em(nex)

!     Local Variables:
!     ios      - i/o error flag
!     iUGrid   - unit number for grid.inp
!     RealE    - real part of energy
!     ImagE    - imaginary part of energy
!     nWords   - number of words in the line
!     Words(4) - array of words
!     line     - string to hold line
!     ieMin    - index of minimum of user defined grid
!     ieMax    - index of max of user defined frid
      INTEGER ios, iUGrid, nWords
c      INTEGER ieMin, ieMax
      DOUBLE PRECISION RealE, ImagE
      CHARACTER(20) Words(10)
      CHARACTER(100) line

!     Loop Variables:
      INTEGER i1, i2

!     Externals
      LOGICAL isnum
      EXTERNAL isnum

      iUGrid = 22
      OPEN(unit=iUGrid,file='grid.inp',status='old',iostat=ios)
      CALL CHOPEN(ios, 'grid.inp', 'xsph')
      
      DO nGrid = 1, nGridMax
!        Read comment lines
         CALL rdcmt(iUGrid,'#!*C')
!        Read data line into string variable "line" and change to
!        lowercase.
         READ(iUGrid,'(A)',END=5) line
c         CALL lower(line)
!        bwords breaks line into words which are then passed
!        back in Words array
         nWords = 4
         CALL untab(line)
         CALL bwords(line,nWords,Words)

!        Set iGridType
         IF(Words(1).eq.'usergrid') THEN
            iGridType(nGrid) = 0
         ELSEIF(Words(1).eq.'egrid') THEN
            iGridType(nGrid) = 1
c            IF(nwords.ne.4) CALL GridError('Error in grid.inp', line)
         ELSEIF(Words(1).eq.'kgrid') THEN
            iGridType(nGrid) = 2
         ELSEIF(Words(1).eq.'expgrid') THEN
            iGridType(nGrid) = 3
         END IF
         
         IF(iGridType(nGrid).ne.0) THEN
            IF(Words(2).eq.'last') THEN
               ! Set the grid minimum to the max of the last grid.
               IF(nGrid.gt.1) THEN
                  CALL SetGridMin(GridMin,GridMax,GridStep,iGridType,
     &                 nGrid)
               ELSE
                  GridMin(1) = 0.d0
               END IF
            ELSE
               READ(Words(2),*) GridMin(nGrid)
            END IF
            READ(Words(3),*) GridMax(nGrid)
            READ(Words(4),*) GridStep(nGrid)
         END IF

         IF(iGridType(nGrid).eq.0) THEN
!        User defined points: read from file.
            DO i2 = 1, nex
               ! Read comments
               CALL rdcmt(iUGrid,'#!*C')
               ! Read line
               READ(iUGrid,'(A)',END=5) line
               nwords = 2
               ! break line into words
               CALL untab(line)
               CALL bwords(line,nWords,Words)
               ! if first word is number, Real(em) = num
               IF(isnum(Words(1))) THEN
                  READ(Words(1),*) RealE
                  ! if second word exists and is a num, Im(em) = num
                  IF((nWords.ge.2).and.isnum(Words(2)))
     &                 READ(Words(2),*) ImagE                  
                  em(i2) = (RealE + coni*ImagE)
                  ne = ne + 1
               ! If first word is not a number, exit loop and read line again.   
               ELSE
                  ! Set GridMax and GridMin for reference
                  GridMin(nGrid) = DBLE(em(ne - i2 + 1))
                  GridMax(nGrid) = DBLE(em(ne))

                  BACKSPACE(iUGrid)
                  EXIT
               END IF
            END DO
         END IF
      END DO
 5    CONTINUE
      nGrid = nGrid - 1

      DO i1 = 1, nGrid
         IF(iGridType(i1).eq.2) THEN
!     k-Grid. Set units to bohr**(-1)
            GridMin(i1) = GridMin(i1)*bohr
            GridMax(i1) = GridMax(i1)*bohr
            GridStep(i1) = GridStep(i1)*bohr
         ELSE
!     e-grid. Set units to hartrees
            GridMin(i1) = GridMin(i1)/hart
            GridMax(i1) = GridMax(i1)/hart
            GridStep(i1) = GridStep(i1)/hart
         END IF
      END DO
      DO i1 = 1, ne
         em(i1) = em(i1)/hart
      END DO

      CLOSE(iUGrid)
      RETURN
      END

      SUBROUTINE SetGridMin(GridMin, GridMax, GridStep, iGridType,
     &     nGrid)
!     This sets the minimum of the current grid to the maximum of the last grid + GridStep
      include '../HEADERS/const.h'
!     Input:
!     GridMin   - array that holds grid minima
!     GridMax   - array that holds grid maxima
!     GridStep  - array of steps
!     iGridType - array of grid types
!     nGrid     - current grid
      INTEGER nGrid
      INTEGER iGridType(nGrid)
      DOUBLE PRECISION GridMin(nGrid), GridMax(nGrid), GridStep(nGrid)

!     Output: GridMin(nGrid) (minimum of current grid.

      IF((iGridType(nGrid).ne.2).and.(iGridType(nGrid-1).ne.2).or.
     &    (iGridType(nGrid).eq.iGridType(nGrid-1))) THEN
!     If neither grid is a k grid, or if both are k-grid, just set the minimum to the previous
!     maximum.
         GridMin(nGrid) = GridMax(nGrid-1) + GridStep(nGrid)
      ELSEIF(iGridType(nGrid).eq.2) THEN
!     If current grid is k, kmin = sqrt(2*emax)
         GridMin(nGrid) =
     &        SQRT(2*GridMax(nGrid-1)/hart)/bohr + GridStep(nGrid)
      ELSE
!     If current grid is e, emin = k**2/2
         GridMin(nGrid) =
     &        (GridMax(nGrid-1)*bohr)**2/2*hart + GridStep(nGrid)
      END IF

      RETURN
      END

      SUBROUTINE GridError(message, line)
      CHARACTER(300) message, line
      
      CALL wlog(message)
      CALL wlog(line)
      STOP

      RETURN
      END
