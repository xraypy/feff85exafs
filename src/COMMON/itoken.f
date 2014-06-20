      function itoken (word,flname)
c     chars in word assumed upper case, left justified
c     returns 0 if not a token, otherwise returns token

      character*(*) word
      character*4   w
      character*20 flname
      integer itoken

      itoken = 0
      w = word(1:4)
      call upper(w)
      
c     Tokens for feff.inp
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (flname(1:8).eq.'feff.inp') then
         if     (w .eq. 'ATOM')  then
            itoken = 1
         elseif (w .eq. 'HOLE')  then
            itoken = 2
         elseif (w .eq. 'OVER')  then
            itoken = 3
         elseif (w .eq. 'CONT')  then
            itoken = 4
         elseif (w .eq. 'EXCH')  then
            itoken = 5
         elseif (w .eq. 'ION ')  then
            itoken = 6
         elseif (w .eq. 'TITL')  then
            itoken = 7
         elseif (w .eq. 'FOLP')  then
            itoken = 8
         elseif (w .eq. 'RPAT' .or. w .eq. 'RMAX')  then
            itoken = 9
         elseif (w .eq. 'DEBY')  then
            itoken = 10
         elseif (w .eq. 'RMUL')  then
            itoken = 11
         elseif (w .eq. 'SS  ')  then
            itoken = 12
         elseif (w .eq. 'PRIN')  then
            itoken = 13
         elseif (w .eq. 'POTE')  then
            itoken = 14
         elseif (w .eq. 'NLEG')  then
            itoken = 15
         elseif (w .eq. 'CRIT')  then
            itoken = 16
         elseif (w .eq. 'NOGE')  then
            itoken = 17
         elseif (w .eq. 'IORD')  then
            itoken = 18
         elseif (w .eq. 'PCRI')  then
            itoken = 19
         elseif (w .eq. 'SIG2')  then
            itoken = 20
         elseif (w .eq. 'XANE')  then
            itoken = 21
         elseif (w .eq. 'CORR')  then
            itoken = 22
         elseif (w .eq. 'AFOL')  then
            itoken = 23
         elseif (w .eq. 'EXAF')  then
            itoken = 24
         elseif (w .eq. 'POLA')  then
            itoken = 25
         elseif (w .eq. 'ELLI')  then
            itoken = 26
         elseif (w .eq. 'RGRI')  then
            itoken = 27
         elseif (w .eq. 'RPHA')  then
            itoken = 28
         elseif (w .eq. 'NSTA')  then
            itoken = 29
         elseif (w .eq. 'NOHO')  then
            itoken = 30
         elseif (w .eq. 'SIG3')  then
            itoken = 31
         elseif (w .eq. 'JUMP')  then
            itoken = 32
         elseif (w .eq. 'MBCO')  then
            itoken = 33
         elseif (w .eq. 'SPIN')  then
            itoken = 34
         elseif (w .eq. 'EDGE')  then
            itoken = 35
         elseif (w .eq. 'SCF ')  then
            itoken = 36
         elseif (w .eq. 'FMS ')  then
            itoken = 37
         elseif (w .eq. 'LDOS')  then
            itoken = 38
         elseif (w .eq. 'INTE')  then
            itoken = 39
         elseif (w .eq. 'CFAV')  then
            itoken = 40
         elseif (w .eq. 'S02 ')  then
            itoken = 41
         elseif (w .eq. 'XES ')  then
            itoken = 42
         elseif (w .eq. 'DANE')  then
            itoken = 43
         elseif (w .eq. 'FPRI')  then
            itoken = 44
         elseif (w .eq. 'RSIG')  then
            itoken = 45
         elseif (w .eq. 'XNCD')  then
            itoken = 46
         elseif (w .eq. 'XMCD')  then
            itoken = 46
         elseif (w .eq. 'MULT')  then
            itoken = 47
         elseif (w .eq. 'UNFR')  then
            itoken = 48
         elseif (w .eq. 'TDLD')  then
            itoken = 49
         elseif (w .eq. 'PMBS')  then
            itoken = 50
         elseif (w .eq. 'PLAS')  then
            itoken = 51
         elseif (w .eq. 'SO2C')  then
            itoken = 52
         elseif (w .eq. 'SELF')  then
            itoken = 53
         elseif (w .eq. 'SFSE')  then
            itoken = 54
         elseif (w .eq. 'RCONV') then
            itoken = 55
         elseif (w .eq. 'ELNE') then !KJ new card for EELS 1-06
            itoken = 56
         elseif (w .eq. 'EXEL') then !KJ new card for EELS 1-06
            itoken = 57
         elseif (w .eq. 'MAGI') then !KJ new card for EELS 1-06
            itoken = 58
         elseif (w .eq. 'ABSO') then !KJ new card 3-06
            itoken = 59    
         elseif (w .eq. 'EGRI')  then !Josh Kas
            itoken = 60
         elseif (w .eq. 'END ')  then
            itoken = -1            
         else
            itoken = 0
         endif
      elseif (flname(1:10).eq.'spring.inp') then
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     These tokens are for spring.inp (input for eq of motion method)
         if (w .eq. 'STRE')  then
            itoken = 1
         elseif (w .eq. 'ANGL')  then
            itoken = 2
         elseif (w .eq. 'VDOS')  then
            itoken = 3
         elseif (w .eq. 'PRDOS') then
            itoken = 4
         elseif (w .eq. 'END ')  then
            itoken = -1            
         else
            itoken = 0
         endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      endif
      return
      end


