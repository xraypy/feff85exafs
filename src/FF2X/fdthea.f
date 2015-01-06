      subroutine fdthea(ntext, text, ip, iorder, nleg, deg, reff,
     &       rnrmav, edge, rat, ipot, iz, potlbl, nlines, lines)

c+---------------------------------------------------------------------
c     "Based on or developed using Distribution: FEFF8.5L
c      Copyright (c) [2013] University of Washington"
c
C  See ../HEADERS/license.h for full llicense information
c+---------------------------------------------------------------------
      implicit double precision (a-h, o-z)
      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'
      parameter (npx=15000)

      double precision rat(3,legtot)
      dimension ipot(legtot)
      dimension iz(0:nphx)
      character*6  potlbl(0:nphx)
      character*80 text(nheadx), lines(2*nheadx), thisln
      real deg, reff, rnrmav, edge
      
      external istrln

      nlines = 1
      do 300  itext = 1, ntext
         ltxt = istrln(text(itext))
         write(thisln,160)  text(itext)(1:ltxt)
         lines(nlines) = thisln
         nlines = nlines+1
 300  continue
 160  format (1x, a)

      write(thisln,310) ip, iorder
 310  format (' Path', i5, '      icalc ', i7)
      lines(nlines) = thisln
      nlines = nlines+1

      write(thisln,170)
  170 format (1x, 71('-'))
      lines(nlines) = thisln
      nlines = nlines+1

      write(thisln,320)  nleg, deg, reff*bohr, rnrmav, 
     1       edge*hart
 320  format (1x, i3, f8.3, f9.4, f10.4, f11.5, 
     1       ' nleg, deg, reff, rnrmav(bohr), edge')
      lines(nlines) = thisln
      nlines = nlines+1

      write(thisln,330)
 330  format ('        x         y         z   pot at#')
      lines(nlines) = thisln
      nlines = nlines+1

      write(thisln,340)  (rat(j,nleg)*bohr,j=1,3), ipot(nleg),
     1       iz(ipot(nleg)), potlbl(ipot(nleg))
 340  format (1x, 3f10.4, i3, i4, 1x, a6, '   absorbing atom')
      lines(nlines) = thisln
      nlines = nlines+1

      do 360  ileg = 1, nleg-1
         write(thisln,350)  (rat(j,ileg)*bohr,j=1,3), ipot(ileg),
     1          iz(ipot(ileg)), potlbl(ipot(ileg))
 350     format (1x, 3f10.4, i3, i4, 1x, a6)
         lines(nlines) = thisln
         nlines = nlines+1
 360  continue
      
      write(thisln,370)
 370  format    ('    k   real[2*phc]   mag[feff]  phase[feff]',
     1       ' red factor   lambda     real[p]@#')
      lines(nlines) = thisln

      return
      end
