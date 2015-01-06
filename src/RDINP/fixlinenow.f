      subroutine fixlinenow(w,nw)
      implicit none 
      integer nwordx     
      parameter (nwordx = 20)
      character*20 w(nwordx),k
      integer nw,nwold,i
      logical iscomm
      external iscomm
      
      nwold=nw
      do i=1,nw
        k=w(i)
        call untab(k)
        call triml(k)
        w(i)=k
!       write(*,'(a1,a20,a1,i1)') '-',w(i),'-',iscomm(w(i))
        if (iscomm(w(i))) then
          nw=i-1
          exit
        endif
      enddo
      
      if (nw.ne.nwold) then
        do i=nw+1,nwold
        
          w(i)='                    '
        enddo
      endif
      
      
      return
      end
      
