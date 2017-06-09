      subroutine wphase (nph, em, eref, lmax, ne, ph)

c     Writes phase data to file PHASExx.DAT for each shell

      implicit double precision (a-h, o-z)

      include 'dim.h'

      complex*16 eref(nex)
      complex*16 ph(nex,ltot+1,0:nphx)
      dimension em(nex)
      dimension lmax(0:nphx)
      character*30  fname

c     Dump phase data, eref and complex phase for each shell
      do 260  iph = 0, nph
c        prepare file for shell's phase data
         write(fname,242)  iph
  242    format('phase', i2.2, '.dat')
         open (unit=1, file=fname, status='unknown', iostat=ios)
         call chopen (ios, fname, 'wphase')
         call wthead  (1)
c        write out unique pot and lmax
         write(1,244)   iph, lmax(iph), ne
  244    format (1x, 3i4, '   unique pot,  lmax, ne ')
c        for each energy
c        ie, em, eref, p=sqrt(em-eref)
c        ph array to ltot+1, 5 values per line
         do 250  ie = 1, ne
            xp = sqrt(em(ie) - eref(ie))
            write(1,246)  ie, em(ie), eref(ie), sqrt(em(ie)-eref(ie))
  246       format ('   ie        energy      re(eref)',
     1              '      im(eref)',
     2              '         re(p)         im(p)', /,
     3              1x, i4, 1p, 5e14.6)
            write(1,248)  (ph(ie,ll,iph), ll=1,lmax(iph)+1)
  248       format (1x, 1p, 4e14.6)
  250    continue
         close(unit=1)
  260 continue

      return
      end
