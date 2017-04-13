c SUBROUTINE SUMAX (NPTS, RN, ANN, AA2, AASUM)
c This is a version of the subroutine sumax found on page 110 of
c Louck's book.  It performs eq 3.22, using simpson's rule and
c taking advantage of the logarithmic grid so that sum f(r)*dr becomes
c sum over f(r)*r*(0.05).  Linear interpolation is used at the end
c caps.  This version does not sum over 14 shells of identical
c atoms, instead it averages the contribution of one or more atoms
c of type 2 at the location of atom 1.  Louck's description (except
c for his integration algorithm) is very clear.
c
c input:  npts      number of points to consider
c         rn        distance from atom 1 to atom 2 in au
c         ann       number of type 2 atoms to add to atom 1, can
c                   be fractional
c         aa2(i)    potential or density at atom 2
c output: aasum(i)  spherically summed contribution added into this
c                   array so that sumax can be called repeatedly
c                   and the overlapped values summed into aasum
c
c Note that this routine requires that all position data be on a
c grid  rr(j) = exp (-8.8d0 + (j-1)*0.05d0), which is the grid
c used by Louck, and also used by ATOM if nuclear options not used.
c
c Coded by Steven Zabinsky, December 1989
c Modified for FEFF cluster code, August 1990, siz
c Bug fixed, May 1991, SIZ
c Another bug fixed, Mar 1992, SIZ
c
c T.L.Louck, "Augmented Plane Wave Method", W.A.Benjamin, Inc., 1967

      subroutine sumax (npts, rn, ann, aa2, aasum)
      implicit double precision (a-h, o-z)
      parameter (nptx=250)
      dimension aa2(nptx), aasum(nptx)
      dimension stor(nptx)
       external xx
c     jjchi     index beyond which aa2 is zero
c     jtop      index just below distance to neighbor
c               aasum is calculated only up to index jtop

c     Wigner-Seitz radius is set to 15 in ATOM.
      rws = 15
      jjchi = ii(rws)
      jtop  = ii(rn)

      topx = xx(jjchi)

      do 120  i = 1, jtop
         x = xx(i)
         xint = 0.0
         et = exp(x)
         blx = log(rn-et)
         if (blx .ge. topx)  goto 119
         jbl = 2.0+20.0*(blx+8.8)
         if (jbl .lt. 1)  jbl=1
         if (jbl .ge. 2)  then
c           use linear interp to make end cap near center of neighbor
            xjbl = jbl
            xbl = 0.05 * (xjbl-1.0) - 8.8
            g = xbl-blx
            xint = xint+0.5*g*(aa2(jbl)*(2.0-20.0*g)*exp(2.0*xbl)
     1             +20.0*g*aa2(jbl-1)*exp(2.0*(xbl-0.05)))
         endif
         tlx = log(rn+et)
         if (tlx .ge. topx)  then
            jtl = jjchi
            go to 90
         endif
         jtl = 1.0 + 20.0*(tlx+8.8)
         if (jtl .lt. jbl)  then
c           handle peculiar special case at center of atom 1
            fzn = aa2(jtl)*exp(2.0*(xbl-0.05))
            fz3 = aa2(jbl)*exp(2.0*xbl)
            fz2 = fzn+20.0*(fz3-fzn)*(tlx-xbl+0.05)
            fz1 = fzn+20.0*(fz3-fzn)*(blx-xbl+0.05)
            xint = 0.5*(fz1+fz2)*(tlx-blx)
            go to 119
         endif
         xjtl = jtl
         xtl = 0.05*(xjtl-1.0)-8.8
         c = tlx-xtl
         xint = xint+0.5*c*(aa2(jtl)*(2.0-20.0*c)
     1         *exp(2.0*xtl)+aa2(jtl+1)*20.0*c
     2         *exp(2.0*(xtl+0.05)))

   90    if (jtl .gt. jbl)  then
  100       xint = xint+0.5*(aa2(jbl)*exp(2.0*xbl)+aa2(jbl+1)
     1             *exp(2.0*(xbl+0.05)))*0.05
            jbl = jbl+1
            if (jbl .lt. jtl) then
               xbl = xbl+0.05
               go to 100
            endif
         endif
  119    stor(i) = 0.5*xint*ann/(rn*et)
  120 continue

      do 190  i = 1, jtop
         aasum(i) = aasum(i) + stor(i)
  190 continue

      return
      end
