      subroutine chopen (ios, fname, mod)
c     Writes error msg and stops if error in ios flag from open
c     statement.  fname is filename, mod is module with failed open.
      character*(*) fname, mod, s*512
       integer il, istrln
c     open successful
      if (ios .le. 0)  return

c     error opening file, tell user and die.
       il = istrln(fname)
       write(s,100) fname(1:il), mod
 100   format ('cannot open file "',
     $      a, '" in module ', a)
       call fstop(s)
       end
