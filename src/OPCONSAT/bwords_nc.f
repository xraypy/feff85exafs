      subroutine bwords_nc (s, nwords, words)
c
c     Breaks string into words.  Words are seperated by one or more
c     blanks or tabs, but not a comma.  This is just like the normal
c     bwords, excpet for the comma
c
c     ARGS        I/O      DESCRIPTION
c     ----        ---      -----------
c     S            I       CHAR*(*)  String to be broken up
c     NWORDS      I/O      Input:  Maximum number of words to get
c                          Output: Number of words found
c     WORDS(NWORDS) O      CHAR*(*) WORDS(NWORDS)
c                          Contains words found.  WORDS(J), where J is
c                          greater than NWORDS found, are undefined on
c                          output.
c
c      Written by:  Steven Zabinsky, September 1984
c      Tab char added July 1994.
c
c**************************  Deo Soli Gloria  **************************

c  -- No floating point numbers in this routine.
      implicit integer (a-z)

      character*(*) s, words(nwords)

      character blank, tab
      parameter (blank = ' ', tab = '	')
c     there is a tab character here               ^.

c  -- BETW    .TRUE. if between words
c     COMFND  .TRUE. if between words and a comma has already been found
      logical betw, comfnd

c  -- Maximum number of words allowed
      wordsx = nwords

c  -- SLEN is last non-blank character in string
      slen = istrln (s)

c  -- All blank string is special case
      if (slen .eq. 0)  then
         nwords = 0
         return
      endif

c  -- BEGC is beginning character of a word
      begc = 1
      nwords = 0

      betw   = .true.
      comfnd = .true.

      do 10  i = 1, slen
         if (s(i:i) .eq. blank .or. s(i:i) .eq. tab)  then
            if (.not. betw)  then
               nwords = nwords + 1
               words (nwords) = s (begc : i-1)
               betw = .true.
               comfnd = .false.
            endif
         else
            if (betw)  then
               betw = .false.
               begc = i
            endif
         endif

         if (nwords .ge. wordsx)  return

   10 continue

      if (.not. betw  .and.  nwords .lt. wordsx)  then
         nwords = nwords + 1
         words (nwords) = s (begc :slen)
      endif

      return
      end
