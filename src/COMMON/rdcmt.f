      subroutine rdcmt(iunt, cmt)
      integer iunt, i1
c      character(300) line
      character(4) cmt
      character tmpcmt(4), ch
      logical cmtlin

      cmtlin = .true.
      do i1 = 1, 4
         tmpcmt(i1) = cmt(i1:i1)
      end do
 5    continue
      read(iunt,*,end=10) ch
      do i1 = 1, 4
         if(ch.eq.tmpcmt(i1)) goto 5
      end do
      
 10   continue
      backspace(iunt)
      
      return
      end
