      double precision function determ(array,nord,nrows)
c
c     calculate determinate of a square matrix
c        (from bevington "data reduction and error analysis
c         for the physical sciences" pg 294)
c     array: matrix to be analyzed
c     nord: order of matrix
c     nrows:  first dimension of matrix in calling routine
c
      double precision array(nrows,nrows)
      double precision saved
      determ = 1.d0
      do k=1,nord
        if (array(k,k).ne.0) go to 130
        do j=k,nord
           if (array(k,j).ne.0) go to 110
        enddo
        determ = 0.d0
        go to 160
c
 110    continue
        do i=k,nord
           saved = array(i,j)
           array(i,j) = array(i,k)
           array(i,k) = saved
        enddo
        determ = -determ
c
 130    continue
        determ = determ*array(k,k)
        if (k.lt.nord) then
           k1 = k+1
           do i=k1,nord
              do j=k1,nord
                 array(i,j)=array(i,j)-array(i,k)*array(k,j)/array(k,k)
              enddo
           enddo
        endif
      enddo
 160  return
c     end double precision function determ
      end
