      subroutine fdtxdi(ntext, text, ip, iorder, nleg, deg, reff,
     &       rnrmav, edge, rat, ipot, iz, potlbl, nlines, lines)

c+---------------------------------------------------------------------
c     write a feffNNNN.dat equivalent in the XDI format
c     all header information is treated as XDI metadata in the
c     "feff8l" namespace, save for the absorber element and edge,
c     extracted from the header line for potential 0 and presented
c     as Element.edge and Element.symbol
c+---------------------------------------------------------------------
c     "Based on or developed using Distribution: FEFF8.5L
c      Copyright (c) [2013] University of Washington"
c
C  See ../HEADERS/license.h for full llicense information
c+---------------------------------------------------------------------

      implicit double precision (a-h, o-z)
      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'
      include '../HEADERS/vers.h'
      parameter (npx=15000)

      double precision rat(3,legtot)
      dimension ipot(legtot)
      dimension iz(0:nphx)
      character*6  potlbl(0:nphx)
      character*2  atsym
      character*10 v
      character*80 text(nheadx), lines(2*nheadx), buff, thisln
      real deg, reff, rnrmav, edge

      parameter (nwordx = 20)
      character*20 words(nwordx)
      
      external istrln, triml, atsym

      nlines = 1
 10   format('# XDI/1.0 feff8l/', a10)
      v = vf85e
      call triml(v)
      write(thisln,10) v
      lines(nlines) = thisln
      nlines = nlines+1

c     absorber line, eg: Abs   Z=29 Rmt= 0.000 Rnm= 0.000 K  shell
      buff = text(2)
      nwords = nwordx
      call bword2(buff, nwords, words)
 20   format('# Element.edge:    ', a2)
      write(thisln,20) words(8)
      lines(nlines) = thisln
      nlines = nlines+1

 22   format('# Element.symbol:  ', a2)
      read(words(3), *) i
      write(thisln,22) atsym(i)
      lines(nlines) = thisln
      nlines = nlines+1

      lines(nlines) = '# Column.1:        k inverse Angstrom'
      nlines = nlines+1
      lines(nlines) = '# Column.2:        real_phc'
      nlines = nlines+1
      lines(nlines) = '# Column.3:        magnitude_feff'
      nlines = nlines+1
      lines(nlines) = '# Column.4:        phase_feff'
      nlines = nlines+1
      lines(nlines) = '# Column.5:        reduction_factor'
      nlines = nlines+1
      lines(nlines) = '# Column.6:        lambda'
      nlines = nlines+1
      lines(nlines) = '# Column.7:        real_p'
      nlines = nlines+1
      
 40   format('# feff8l.path_index: ', i4)
      write(thisln,40) ip
      lines(nlines) = thisln
      nlines = nlines+1

 200  format('# feff8l.nleg:        ', i1)
      write(thisln,200) nleg
      lines(nlines) = thisln
      nlines = nlines+1

 210  format('# feff8l.degen:       ', f6.3)
      write(thisln,210) deg
      lines(nlines) = thisln
      nlines = nlines+1
      
 220  format('# feff8l.reff:        ', f7.4, ' Angstrom')
      write(thisln,220) reff*bohr
      lines(nlines) = thisln
      nlines = nlines+1
      
 50   format('# feff8l.iorder:     ', i2)
      write(thisln,50) iorder
      lines(nlines) = thisln
      nlines = nlines+1

c     still on Abs line
 60   format('# feff8l.abs_radii:   ', a5, '/', a5)
      write(thisln,60) words(5), words(7)
      lines(nlines) = thisln
      nlines = nlines+1
      
      do 199 i=3,ntext
         buff = text(i)
         nwords = nwordx
         call bword2(buff, nwords, words)
         if (words(1) .eq. 'Gam_ch') then

 100        format('# feff8l.gamma_ch:    ', a10, ' eV')
            write(thisln, 100) words(2)
            lines(nlines) = thisln
            nlines = nlines+1
            
 110        format('# feff8l.exchange:    ', a10)
            write(thisln, 110) words(3)
            lines(nlines) = thisln
            nlines = nlines+1
            
         elseif (words(1) .eq. 'Mu') then

 120        format('# feff8l.mu:          ', a10, ' eV')
            write(thisln, 120) words(2)
            lines(nlines) = thisln
            nlines = nlines+1
            
 130        format('# feff8l.kf:          ', a10,
     1             ' inverse Angstrom')
            write(thisln, 130) words(4)
            lines(nlines) = thisln
            nlines = nlines+1
            
 140        format('# feff8l.vint:        ', a10, ' eV')
            write(thisln, 140) words(6)
            lines(nlines) = thisln
            nlines = nlines+1
            
 150        format('# feff8l.rs_int:      ', a10, ' Angstrom')
            write(thisln, 150) words(8)
            lines(nlines) = thisln
            nlines = nlines+1

c        this is a POT line
         elseif (words(1) .eq. 'Pot')  then

 160        format('# feff8l.pot', a1, ':        ', a2)
            read(words(4), *) j
            write(thisln, 160) words(2), atsym(j)
            lines(nlines) = thisln
            nlines = nlines+1

 170        format('# feff8l.pot', a1, '_radii:  ', a5, '/', a5)
            write(thisln,170) words(2), words(6), words(8)
            lines(nlines) = thisln
            nlines = nlines+1
            
         end if
 199  continue


 230  format('# feff8l.rnrmav:      ', f7.4, ' Angstrom')
      write(thisln,230) rnrmav
      lines(nlines) = thisln
      nlines = nlines+1
      
 240  format('# feff8l.edge:        ', f9.5, ' eV')
      write(thisln,240) edge*hart
      lines(nlines) = thisln
      nlines = nlines+1


 300  format('# feff8l.atom0:       ', 3f10.4, i3, i4, 1x, a6)
      write(thisln,300)  (rat(j,nleg)*bohr,j=1,3), ipot(nleg),
     1       iz(ipot(nleg)), potlbl(ipot(nleg))
      lines(nlines) = thisln
      nlines = nlines+1
      
      do 399  ileg = 1, nleg-1
 310     format('# feff8l.atom', i1,
     1          ':       ', 3f10.4, i3, i4, 1x, a6)
         write(thisln,310)  ileg, (rat(j,ileg)*bohr,j=1,3), ipot(ileg),
     1          iz(ipot(ileg)), potlbl(ipot(ileg))
         lines(nlines) = thisln
         nlines = nlines+1
 399  continue

      lines(nlines) = '# -----'
      nlines = nlines+1
      write(thisln,370)
 370  format    ('#   k   real[2*phc]   mag[feff]  phase[feff]',
     1       ' red factor   lambda     real[p]')
      lines(nlines) = thisln

      
      return
      end
      



      SUBROUTINE BWORD2 (S, NWORDS, WORDS)
C
C     Breaks string into words.  Words are seperated by one or more
C     blanks or tabs, or a comma or an equal sign and zero or more blanks.
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
C      Equal char added by BR September 2015
C
C**************************  Deo Soli Gloria  **************************

C  -- No floating point numbers in this routine.
      IMPLICIT INTEGER (A-Z)

      CHARACTER*(*) S, WORDS(NWORDS)

      CHARACTER BLANK, COMMA, EQUAL, TAB
      PARAMETER (BLANK = ' ', COMMA = ',', EQUAL = '=', TAB = '	')
C     there is a tab character here                            ^.

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
         ELSEIF (S(I:I) .EQ. COMMA .OR. S(I:I) .EQ. EQUAL)  THEN
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
      
