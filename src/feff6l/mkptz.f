      subroutine mkptz
c     makes polarization temsor ptz if necessary
      implicit double precision (a-h, o-z)

c     all input and output through common area /pol/
      include 'pola.h'
      include 'const.h'

c     addittonal local stuff to create polarization tensor ptz(i,j)
      real e2(3)
      complex*16  e(3),eps,epc
      dimension eps(-1:1),epc(-1:1)
       character*128 messag

c     Begin to make polarization tensor
c     Normalize polarization vector
      x = sqrt (evec(1)**2 + evec(2)**2 + evec(3)**2)
      if (x .eq. 0.0) then
         call fstop(
     $        ' at MKPTZ: Polarization vector has zero length')
      endif
      do 290  i = 1, 3
         evec(i) = evec(i) / x
  290 continue
      if (elpty .eq. 0.0) then
c        run linear polarization code
         do 291 i = 1, 3
            ivec(i) = 0.0
  291    continue
      endif
      x = sqrt (ivec(1)**2 + ivec(2)**2 + ivec(3)**2)
      if (x .gt. 0) then
c        run elliptical polarization code
         do 293  i = 1, 3
            ivec(i) = ivec(i) / x
  293    continue
         x = evec(1)*ivec(1)+evec(2)*ivec(2)+evec(3)*ivec(3)
 2935    format(1x,a,3g14.5)
         if (abs(x) .gt. 0.9) then
            write(messag,2935) ' polarization  ',
     $           evec(1), evec(2), evec(3)
            call echo(messag)
            write(messag,2935) ' incidence  ',
     $           ivec(1), ivec(2), ivec(3)
            call fstop(' at MKPTZ: Polarization '//
     $           'almost parallel to the incidence')
         endif
         if (x .ne. 0.0) then
c          if ivec not normal to evec then make in normal, keeping the
c          plane based on two vectors
           do 294 i = 1,3
              ivec(i) = ivec(i) - x*evec(i)
  294      continue
           x = sqrt (ivec(1)**2 + ivec(2)**2 + ivec(3)**2)
           do 295  i = 1, 3
              ivec(i) = ivec(i) / x
  295      continue
         endif
      else
         elpty = 0.0
      endif 
     
      e2(1) = ivec(2)*evec(3)-ivec(3)*evec(2)
      e2(2) = ivec(3)*evec(1)-ivec(1)*evec(3)
      e2(3) = ivec(1)*evec(2)-ivec(2)*evec(1)
      do 296  i = 1,3
        e(i) = (evec(i)+elpty*e2(i)*coni)
  296 continue 
      eps(-1) =  (e(1)-coni*e(2))/sqrt(2.0)
      eps(0)  =   e(3)
      eps(1)  = -(e(1)+coni*e(2))/sqrt(2.0)
      do 297  i = 1,3
        e(i) = (evec(i)-elpty*e2(i)*coni)
  297 continue 
      epc(-1) =  (e(1)-coni*e(2))/sqrt(2.0)
      epc(0)  =   e(3)
      epc(1)  = -(e(1)+coni*e(2))/sqrt(2.0)
      do 298 i = -1,1
      do 298 j = -1,1
c        ptz(i,j) = ((-1.0)**i)*epc(-i)*eps(j)/(1+elpty**2)
c       above - true polarization tensor for given ellipticity, 
c       below - average over left and right in order to have
c       path reversal simmetry
        ptz(i,j) = ((-1.0)**i)*(epc(-i)*eps(j)+eps(-i)*epc(j))
     1               /(1+elpty**2)/2.0
  298 continue
c     end of making polarization tensor

      return
      end
