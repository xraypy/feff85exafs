c     make e mesh for phase
c     input:  iprint, ispec, edge, vi0, gamach, xkmax, xkstep
c     output: ne, ne1, em(ne), ik0 [grid point with k=0]
c             ne -  total number of points in array em
c             ne1 - number of points on horizontal grid 

      subroutine phmesh (iprint, ispec, edge, emu, vi0, gamach,
     1                  xkmax, xkstep, vixan, ne, ne1, em, ik0, ne3)
      implicit double precision (a-h, o-z)
      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'
      complex*16 em(nex), tempc

c     see arrays.h
c#mn
      external getxk

c     nemax - max number of points on horizontal axis
      xloss = gamach/2 + vi0
      if (xloss.lt.0) xloss = 0
      xvert = max(xloss, 0.02/hart)
      xloss = xvert
      aa = 0.5d0
      ne3 = 0
      xim = xloss*aa
      if (vixan.gt.0.0001) xim = vixan
      ik0 = 0

      if (ispec.le.3)  then
c        make energy mesh for XANES with FMS
c        around fermi level step is regular in energy (xloss/2)
c        and regular in k at high energies

c        10 points below Fermi level
         nemax = 10
c        dk = 0.14*bohr
         dk = 2*xkstep
         n1 = int (xim/2/dk**2)         
         n2 = int ( sqrt(n1*2*xim) / dk )
         if ( (dk*(n2+1))**2 .gt. (n1+1)*2*xim ) n1 = n1+1
         n1 = min (n1,nemax)
         do 10 i = 1, n1
  10     em(nemax+1-i) = -xim*i + edge + coni*xloss
         nb = nemax-n1
         do 20 i = 1, nb
  20     em(nb + 1 -i) = -(dk*(n2+i))**2/2 + edge + coni*xloss
         nmin = nemax
         ik0 = nemax+1
      endif

      if (ispec .gt. 0 .and. ispec.le.3)  then
c        make energy mesh for XANES with FMS
c        around fermi level step is regular in energy (xloss/2)
c        and regular in k at high energies
c        90 points above Fermi level
         nemax = 100 - nemax 
         n1 = int (xim/2/xkstep**2)
         n2 = int ( sqrt(n1*2*xim) / xkstep )
         n1 = n1 + 1
         if ( (xkstep*(n2+1))**2 .gt. n1*2*xim ) n1 = n1+1
         n1 = min (n1,nemax)
         if (ispec.ne.2) then
            nb = int (xkmax**2 /xim/2) + 1
         else
            nb = int (abs(edge - xkmax/bohr/hart) /xim) + 1
         endif
         if (nb .le. n1) n1 = nb
         do 30 i = 1, n1
  30     em(nmin+i) = xim*(i-1)
         if (ispec.ne.2) then
            nb = int( xkmax / xkstep)  - n2
         else
            nb = int( sqrt(abs(2*(edge-xkmax/bohr/hart))) / xkstep) - n2
         endif
         nb = min(nb, nemax-n1)
         nb = max(nb,0)
         do 40 i = 1, nb
  40     em(nmin+n1+i) = (xkstep*(n2+i))**2 /2
         ne1 = nmin+n1+nb
         do 50 i = ik0, ne1
  50     em(i) = em(i) + edge + coni*xloss

      elseif (ispec.eq.4) then
c        grid for atomic f' calculation regular in energy
         nemax = 100
         emin = xkmax / bohr /hart
         emax = xkstep / bohr / hart
         ne = 1
         emin = emin - emu + edge
         emax  = emax - emu + edge
         em(1) = emin
         if (emin .lt. emax) then
            if (vixan.le.0.d0) vixan = (emax-emin) / (nemax-1)

  85        ne = ne + 1
            em(ne) = em(ne-1) + vixan
            if ( ne.lt.nemax .and. dble(em(ne)).lt.emax) goto 85
         endif

         ne1 = ne
         nemax = nex-ne
         if (nemax.gt.100) nemax=100
         de = 3.d0 /hart
         elimit = min (2.0d5/hart, 20*emu)
         elimit = max (elimit, 1.0d3/hart)
         elimit = elimit - emu

         ne2 = 0
         ne3 = nemax
         ne = ne1+ne2+ne3
         em(ne1+1) = edge
         do 88 i = 1,ne3-1
            dep = 0
            if (dble(em(ne1+i)).gt.0) 
     1             dep=dble(em(ne1+i) * 
     2             (exp( log( elimit/em(ne1+i) ) / (ne3-i) ) -1))
            if (dep.lt.de) dep = de
            em(ne1+i+1) = em(ne1+i) + dep
  88     continue
      else
c        energy mesh for EXAFS or XANES without FMS
c        20 pts (0 le k le 1.9, delk=0.1 ang(-1) )
c        20 pts (2 le k le 5.8, delk=0.2 ang(-1) )
c         9 pts (6 le k le 10., delk=0.5 ang(-1) )
c        10 pts (11 le k le 20.0, delk=1.0 ang(-1) )
         ne = 0
         if (ispec.lt.0) ne = 10
         nemax = 100
         delk = bohr/10
         do 111 i=1,20
            tempk=(i-1)*delk
            ne = ne+1
            em(ne)=tempk**2/2 +edge + coni*xloss
            if (i.eq.1)  ik0 = ne
  111    continue
         delk = bohr/5
         n2 = 20
         do 112 i=1,n2
            tempk=2*bohr + (i-1)*delk
            ne = ne+1
            em(ne)=tempk**2/2 +edge + coni*xloss
  112    continue
         delk = bohr/2
         do 113 i=1,9
            tempk=6*bohr + (i-1)*delk
            ne = ne+1
            em(ne)=tempk**2/2 +edge + coni*xloss
  113    continue
         delk=bohr
         do 114 i=1,10
            tempk=11*bohr + (i-1)*delk
            ne = ne+1
            em(ne)=tempk**2/2 +edge + coni*xloss
  114    continue

c        while loop
  115    if (tempk. lt.xkmax) then
            tempk = tempk + delk
            ne = ne+1
            em(ne)=tempk**2/2 +edge + coni*xloss
            goto 115
         endif

         ne = min (ne, nemax)
         ne1 = ne
      endif

      if (ispec.le.3)  then
c        make the vertical grid in energy plane
c        first point is at 0.005 ev, second at 0.01 ev and
c        exponential grid with step 0.4 after that up to 50 eV
         tempk = 0.005d0/hart
         em(ne1+1) = edge + coni*tempk
         tempk = tempk*2
         em(ne1+2) = edge + coni*tempk
c        chose delk that point edge+coni*xloss is in the middle of step
c        delk = 0.6 is ok for Cu K edge, but needs more testing
         delk = 0.4d0
         n1 = nint ( log(xloss/tempk)/delk - 0.5d0)
         if (n1.le.0) n1 = 1
         bb = exp(delk)
         aa = 2*xloss /(1+bb)
         aa = aa/bb**n1
         if (aa.le. tempk) aa = aa*bb
cJosh         if (aa.le.tempk .or. aa.ge. xloss) 
cJosh     .     call par_stop(' Bad mesh in phmesh')
c        delk = log (xloss/tempk) /(n1+0.5)
c        n1 = nint( log(1000/hart/tempk) / delk )
c        n1 = nint( log(50/hart/aa) / delk )
         ee = min(50.d0/hart,20.d0*xloss)
c         print*, ee, aa, xloss, n1
         n1 = nint( log(ee/aa) / delk )
         do 60 i = 0, n1
  60     em(ne1+3+i) = edge +coni*aa*exp(delk*i)
         ne = ne1 + n1 + 3

c        for DANES need additional points
         if (abs(ispec).eq.3) then
            ne3 = min(nex,150) - ne
            em(ne+1) = dble(2*em(ne1)-em(ne1-1))
            dk = log(7.d4/dble(em(ne+1))) / (ne3-1)
            dk = exp(dk)
            do 80 i = 1, ne3-1
  80        em(ne+i+1)= em(ne+i)*dk
            do 90 i = 1, ne3
  90        em(ne+i)= em(ne+i)+coni*1.d-8
            ne = ne + ne3
            
         endif
      endif

c     need to reverse order for horizontal grid for XES
      if (ispec.eq.2) then
         do 150 ie = 1, ne1
  150    em(ie) = 2*(edge + coni*xloss) - em(ie)
         np = ne1 / 2
         do 160 ie=1,np
            ip = ne1+1-ie
            tempc = em(ie)
            em(ie) = em(ip)
            em(ip) = tempc
  160    continue
         ik0 = ne1+1-ik0
      endif

      if (iprint .ge. 3)  then
         open (unit=44, file='emesh.dat', status='unknown')
         write(44,*) 'edge, bohr, edge*hart ', edge, bohr, edge*hart
         write(44,*) 'ispec, ik0 ', ispec, ik0
         write(44,*) 'ie, em(ie)*hart, xk(ie)'
         do 230  ie = 1, ne
           write(44,220) ie, dble(em(ie))*hart,
     1                   getxk(dble(em(ie))-edge)/bohr
  220      format (i5, 3f20.5)
  230    continue
         close (unit=44)
      endif

      return
      end
