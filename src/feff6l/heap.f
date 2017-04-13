c     These heap routines maintain a heap (array h) and an index
c     array (array ih) used to keep other data associated with the heap
c     elements.

      subroutine hup (h, ih, n)
c     heap is in order except for last element, which is new and must
c     be bubbled through to its proper location
c     new element is at i, j = index of parent
      integer  n,i,j
      integer  ih(n)
      dimension h(n)

      i = n

   10 j = i/2
c     if no parent, we're at the top of the heap, and done
      if (j .eq. 0)  return
      if (h(i) .lt. h(j))  then
         call swap (h(i), h(j))
         call iswap (ih(i), ih(j))
         i = j
         goto 10
      endif
      return
      end

      subroutine hdown (h, ih, n)
c     h is in order, except that 1st element has been replaced.
c     Bubble it down to its proper location.  New element is i,
c     children are j and k.

      integer  n,i,j,k
      integer  ih(n)
      dimension h(n)

      i = 1

   10 continue
      j = 2*i
      k = j + 1

c     if j > n, new element is at bottom, we're done
      if (j .gt. n)  return
c     handle case where new element has only one child
      if (k .gt. n)  k = j

      if (h(j) .gt. h(k))  j = k
c     j is now index of smallest of children

      if (h(i) .gt. h(j))  then
         call swap (h(i), h(j))
         call iswap (ih(i), ih(j))
         i = j
         goto 10
      endif

      return
      end

      subroutine swap (a, b)
      t = a
      a = b
      b = t
      return
      end

      subroutine iswap (i, j)
      integer  i,j,k
      k = i
      i = j
      j = k
      return
      end
