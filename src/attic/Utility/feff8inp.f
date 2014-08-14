      program feff8i
c
c     reads feff.inp file and writes feff8.inp file in new format
      implicit double precision (a-h, o-z)

      character*150  line, linep
      parameter (nwordx = 12)
      character*15 words(nwordx)

      logical iscomm

   10 format (a)
   20 format (bn, i15)
   30 format (bn, f15.0)

      iscmt =0
      ifms  =0

c     Open feff.inp, the input file we're going to read
      open (unit=1, file='feff.inp', status='old', iostat=ios)
      call chopen (ios, 'feff.inp', 'rdinp')

c     Open feff8.inp, the output file we're going to write
      open (unit=2, file='feff8.inp', status='unknown')

c     tokens  0 if not a token
c             4 if CONT (CONTROL)
c             9 if RMAX
c            13 if PRIN (PRINT)
c            36 if SCMT do self-consistency loop
c            37 if FMS  use FMS for cluster of the size rfms
c            38 if LDOS print out l-dos for specified energy range
c            40 if NABS to do configuration average

  200 read(1,10,iostat=ios)  line
         if (ios .lt. 0)  goto 900

         linep = line
         call triml (line)
         if (iscomm(line)) then
            jlen = istrln (linep)
            write (2, 10) linep (1:jlen)
            goto 200
         endif

         nwords = nwordx
         call bwords (line, nwords, words)
         itok = itoken (words(1))

         if (itok .eq. 4)  then
c           CONTROL  mphase, mpath, mfeff, mchi
c            0 - do not run modules, 1 - run module
            if (nwords.eq.5) then
c              feff7 input file
               read(words(2),20,err=900)  mpot
               mphase = mpot
               mfms = mpot
               read(words(3),20,err=900)  mpath
               read(words(4),20,err=900)  mfeff
               read(words(5),20,err=900)  mchi
            else
c              feff8 input file
               read(words(2),20,err=900)  mpot
               read(words(3),20,err=900)  mphase
               read(words(4),20,err=900)  mfms
               read(words(5),20,err=900)  mpath
               read(words(6),20,err=900)  mfeff
               read(words(7),20,err=900)  mchi
            endif
            write (2,100) ' CONTROL', 
     1        mpot, mphase, mfms, mpath, mfeff, mchi
  100       format (a, 6(2x,i1))
         elseif (itok .eq. 9)  then
c           RMAX  rmax (max r for ss and pathfinder)
            read(words(2),30,err=900)  rmax
            write (2, 110) ' RPATH ', rmax
  110       format (a, f8.3)
         elseif (itok .eq. 13)  then
c           PRINT  ipr1  ipr2  ipr3  ipr4 ipr5 ipr6
            if (nwords.eq.5) then
c              feff7 input file
               read(words(2),20,err=900)  ipr1
               ipr2 = ipr1
               ipr3 = ipr1
               read(words(3),20,err=900)  ipr4
               read(words(4),20,err=900)  ipr5
               read(words(5),20,err=900)  ipr6
            else
c              feff8 input file
               read(words(2),20,err=900)  ipr1
               read(words(3),20,err=900)  ipr2
               read(words(4),20,err=900)  ipr3
               read(words(5),20,err=900)  ipr4
               read(words(6),20,err=900)  ipr5
               read(words(7),20,err=900)  ipr6
            endif
            write (2, 100) ' PRINT',ipr1,ipr2,ipr3,ipr4,ipr5,ipr6
         elseif (itok .eq. 36)  then
c           SCMT   rfms [ lfms nscmt  icoul ca1 ]
c           number of cycles, mode of calculating coulomb potential,
c           convergence accelerator
            nscmt = 30
            ca1 = 0.2d0
            icoul = 0
            ecv = 0
            if (nwords.gt.1) read(words(2),20,err=900)  nscmt
            if (nwords.gt.2) read(words(3),30,err=900)  ca1
            if (nwords.gt.3) read(words(4),30,err=900)  ecv
            if (ifms.eq.1) write (2,140) ' SCF ',rfms1, lfms1,
     1         nscmt, ca1, icoul, ecv
            iscmt = 1
         elseif (itok .eq. 37)  then
c           FMS   rfms2  (lfms2)
c           radius of the cluster to do FMS
            read(words(2),30,err=900)  rfms1
            read(words(3),30,err=900)  rfms2
            lfms1=0
            lfms2=0
            if (nwords.gt.3) then
               read(words(4),20,err=900)  lfms1
               lfms1 = 1 - lfms1
            endif
            if (nwords.gt.4) then
               read(words(5),20,err=900)  lfms2
               lfms2 = 1 - lfms2
            endif
            write (2, 130) ' FMS ', rfms2, lfms2
  130       format ( a, f7.3,1x, i2)
            if (iscmt.eq.1) write (2,140) ' SCF ',rfms1, lfms1,
     1         nscmt, ca1, icoul, ecv
  140       format (a, f7.3,1x,i1,1x,i3,1x,f7.3,1x,i2,1x,f6.1)
            ifms = 1
         elseif (itok .eq. 40) then
c           NABS  iphabs nabs rclabs
            read(words(2),20,err=900)  iphabs
            read(words(3),20,err=900)  nabs
            read(words(4),30,err=900)  rclabs
            write (2, 120) ' CAVERAGE ', iphabs, nabs, rclabs
  120       format ( a, i2, i6, f7.3)
         else
            jlen = istrln(linep)
            write(2,10) linep(1:jlen)
         endif
      goto 200

  900 continue
c     We're done reading the input file, close it.
      close (unit=1)
      close (unit=2)

      stop
      end


      function itoken (word)
c     chars in word assumed upper case, left justified
c     returns 0 if not a token, otherwise returns token

      character*(*) word
      character*4   w

      w = word(1:4)
      call upper(w)
      if (w .eq. 'CONT')  then
         itoken = 4
      elseif (w .eq. 'RMAX')  then
         itoken = 9
      elseif (w .eq. 'PRIN')  then
         itoken = 13
      elseif (w .eq. 'SCMT')  then
         itoken = 36
      elseif (w .eq. 'FMS ')  then
         itoken = 37
      elseif (w .eq. 'NABS')  then
         itoken = 40
      else
         itoken = 0
      endif
      return
      end

      logical function iscomm (line)
c     returns true if line is a comment or blank line, false otherwise
      character*(*) line
      iscomm = .false.
      if (istrln(line).le.0  .or.  line(1:1).eq.'*')  iscomm = .true.
      return
      end

      subroutine warnex (string)
      implicit double precision (a-h, o-z)
c     This prints a warning message if the user is using an
c     expert option.
      character*(*) string

      call wlog(string)
      call wlog(' Expert option, please read documentation ' //
     1          'carefully and check your results.')
      return
      end

      subroutine wlog (string)
      character*(*) string

c     This output routine is ued to replace the PRINT statement
c     for output that "goes to the terminal", or to the log file.
c     If you use a window based system, you can modify this routine
c     to handle the running output elegantly.
c     Handle carriage control in the string you pass to wlog.
c
c     The log file is also written here, hard coded here.

c     The log file is unit 11.  The log file is opened in the
c     main program, program feff.

c     make sure not to write trailing blanks


   10 format (a)

      il = istrln (string)
      if (il .eq. 0)  then
         print10
         write(11,10)
      else
         print10, string(1:il)
         write(11,10) string(1:il)
      endif
      return
      end
      subroutine lblank (string)
      character*(*) string
c     add a leading blank, useful for carriage control
      string = ' ' // string
      return
      end
      subroutine chopen (ios, fname, mod)
c     Writes error msg and stops if error in ios flag from open
c     statement.  fname is filename, mod is module with failed open.
      character*(*) fname, mod
      character*512 slog

c     open successful
      if (ios .le. 0)  return

c     error opening file, tell user and die.
      i = istrln(fname)
      j = istrln(mod)
      write(slog,100)  fname(1:i), mod(1:j)
      call wlog(slog)

  100 format (' Error opening file, ', a, 
     2        ' in module ', a)

      call wlog(' Fatal error')
      stop 'CHOPEN'
      end
C FUNCTION ISTRLN (STRING)  Returns index of last non-blank
C                           character.  Returns zero if string is
C                           null or all blank.

      FUNCTION ISTRLN (STRING)
      CHARACTER*(*)  STRING
      CHARACTER BLANK, TAB
      PARAMETER (BLANK = ' ', TAB = '	')
C     there is a tab character here  ^

C  -- If null string or blank string, return length zero.
      ISTRLN = 0
      IF (STRING (1:1) .EQ. CHAR(0))  RETURN
      IF (STRING .EQ. ' ')  RETURN

C  -- Find rightmost non-blank character.
      ILEN = LEN (STRING)
      DO 20  I = ILEN, 1, -1
         IF (STRING(I:I).NE.BLANK .AND. STRING(I:I).NE.TAB)  GOTO 30
   20 CONTINUE
   30 ISTRLN = I

      RETURN
      END
C SUBROUTINE TRIML (STRING)  Removes leading blanks.

      SUBROUTINE TRIML (STRING)
      CHARACTER*(*)  STRING
      CHARACTER*200  TMP
      CHARACTER BLANK, TAB
      PARAMETER (BLANK = ' ', TAB = '	')
C     there is a tab character here  ^

      JLEN = ISTRLN (STRING)

C  -- All blank and null strings are special cases.
      IF (JLEN .EQ. 0)  RETURN

C  -- FInd first non-blank char
      DO 10  I = 1, JLEN
         IF (STRING(I:I).NE.BLANK .AND. STRING(I:I).NE.TAB)  GOTO 20
   10 CONTINUE
   20 CONTINUE

C  -- If I is greater than JLEN, no non-blanks were found.
      IF (I .GT. JLEN)  RETURN

C  -- Remove the leading blanks.
      TMP = STRING (I:)
      STRING = TMP
      RETURN
      END
C SUBROUTINE UPPER (STRING)  Changes a-z to upper case.

      SUBROUTINE UPPER (STRING)
      CHARACTER*(*)  STRING

      JLEN = ISTRLN (STRING)

      DO 10  I = 1, JLEN
         IC = ICHAR (STRING (I:I))
         IF ((IC .LT. 97)  .OR.  (IC .GT. 122))  GOTO 10
         STRING (I:I) = CHAR (IC - 32)
   10 CONTINUE

      RETURN
      END
C SUBROUTINE LOWER (STRING)  Changes A-Z to lower case.

      SUBROUTINE LOWER (STRING)
      CHARACTER*(*)  STRING

      JLEN = ISTRLN (STRING)

      DO 10  I = 1, JLEN
         IC = ICHAR (STRING (I:I))
         IF ((IC .LT. 65) .OR.  (IC .GT. 90))  GOTO 10
         STRING (I:I) = CHAR (IC + 32)
   10 CONTINUE

      RETURN
      END
C***********************************************************************
C
      SUBROUTINE BWORDS (S, NWORDS, WORDS)
C
C     Breaks string into words.  Words are seperated by one or more
C     blanks or tabs, or a comma and zero or more blanks.
C
C     ARGS        I/O      DESCRIPTION
C     ----        ---      -----------
C     S            I       CHAR*(*)  String to be broken up
C     NWORDS      I/O      Input:  Maximum number of words to get
C                          Output: Number of words found
C     WORDS(NWORDS) O      CHAR*(*) WORDS(NWORDS)
C                          Contains words found.  WORDS(J), where J is
C                          greater then NWORDS found, are undefined on
C                          output.
C
C      Written by:  Steven Zabinsky, September 1984
C      Tab char added July 1994.
C
C**************************  Deo Soli Gloria  **************************

C  -- No floating point numbers in this routine.
      IMPLICIT INTEGER (A-Z)

      CHARACTER*(*) S, WORDS(NWORDS)

      CHARACTER BLANK, COMMA, TAB
      PARAMETER (BLANK = ' ', COMMA = ',', TAB = '	')
C     there is a tab character here               ^.

C  -- BETW    .TRUE. if between words
C     COMFND  .TRUE. if between words and a comma has already been found
      LOGICAL BETW, COMFND

C  -- Maximum number of words allowed
      WORDSX = NWORDS

C  -- SLEN is last non-blank character in string
      SLEN = ISTRLN (S)

C  -- All blank string is special case
      IF (SLEN .EQ. 0)  THEN
         NWORDS = 0
         RETURN
      ENDIF

C  -- BEGC is beginning character of a word
      BEGC = 1
      NWORDS = 0

      BETW   = .TRUE.
      COMFND = .TRUE.

      DO 10  I = 1, SLEN
         IF (S(I:I) .EQ. BLANK .OR. S(I:I) .EQ. TAB)  THEN
            IF (.NOT. BETW)  THEN
               NWORDS = NWORDS + 1
               WORDS (NWORDS) = S (BEGC : I-1)
               BETW = .TRUE.
               COMFND = .FALSE.
            ENDIF
         ELSEIF (S(I:I) .EQ. COMMA)  THEN
            IF (.NOT. BETW)  THEN
               NWORDS = NWORDS + 1
               WORDS (NWORDS) = S(BEGC : I-1)
               BETW = .TRUE.
            ELSEIF (COMFND)  THEN
               NWORDS = NWORDS + 1
               WORDS (NWORDS) = BLANK
            ENDIF
            COMFND = .TRUE.
         ELSE
            IF (BETW)  THEN
               BETW = .FALSE.
               BEGC = I
            ENDIF
         ENDIF

         IF (NWORDS .GE. WORDSX)  RETURN

   10 CONTINUE

      IF (.NOT. BETW  .AND.  NWORDS .LT. WORDSX)  THEN
         NWORDS = NWORDS + 1
         WORDS (NWORDS) = S (BEGC :SLEN)
      ENDIF

      RETURN
      END
