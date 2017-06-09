      subroutine rpotph (io, nhead0, head0, lhead0,
     1             nat, nph, nfr, ihole, gamach, iafolp, intclc,
     1             ixc, vr0, vi0, rs0, iphat, rat, iatph, ifrph, 
     1             xnatph, novr,
     2             iphovr, nnovr, rovr, folp, ion, iz, iprint, 
     2             ixanes, nemax, xkmin, xkmax, potlbl)
      implicit double precision (a-h, o-z)

c     Notes:
c        nat   number of atoms in problem
c        nph   number of unique potentials
c        nfr   number of unique free atoms
c        ihole hole code of absorbing atom
c        iph=0 for central atom
c        ifr=0 for central atom
c        xkmin, xkmax  min and max energy mesh points to consider

      include 'dim.h'

      character*(*) head0(nhead0)
      dimension lhead0(nhead0)

c     End of line comments removed -- see include file arrays.h for
c     comments.
c     Specific atom input data
      dimension iphat(natx)
      dimension rat(3,natx)

c     Unique potential input data
      dimension iatph(0:nphx)
      dimension ifrph(0:nphx)
      dimension xnatph(0:nphx)
      character*6  potlbl(0:nphx)

      dimension folp(0:nphx)
      dimension novr(0:nphx)
      dimension iphovr(novrx,0:nphx)
      dimension nnovr(novrx,0:nphx)
      dimension rovr(novrx,0:nphx)

c     Free atom data
      dimension ion(0:nfrx)
      dimension iz(0:nfrx)

c     read and save header from old file, has carriage control char
      head0(1) = ' '
      call rdhead (io, nhead0, head0, lhead0)
      read(io,*) ihole, gamach, iprint, iafolp, intclc
      read(io,*) ixc, vr0, vi0, rs0
      read(io,*) ixanes, nemax, xkmin, xkmax
      read(io,*) nfr
      do 710  ifr = 0, nfr
         read(io,*)  index, iz(ifr), ion(ifr)
  710 continue
      read(io,*) nat
      do 720  iat = 1, nat
         read(io,*) index, iphat(iat), (rat(j,iat),j=1,3)
  720 continue
      read(io,*) nph
      do 740  iph = 0, nph
         read(io,*) index, iatph(iph), ifrph(iph), xnatph(iph),
     1                folp(iph), novr(iph)
         read(io,*) potlbl(iph)
         do 730  iovr = 1, novr(iph)
            read(io,*) iphovr(iovr,iph), nnovr(iovr,iph),
     1                   rovr(iovr,iph)
  730    continue
  740 continue

      return
      end
