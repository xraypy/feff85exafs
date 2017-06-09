      subroutine wpot (nph, edens, ifrph, imt, inrm,
     1                 rho, vclap, vcoul, vtot)

c     Writes potentials to file name POTxx.DAT for each unique pot.

      implicit double precision (a-h, o-z)

      include 'const.h'
      include 'dim.h'

      dimension ifrph(0:nphx)
      dimension rho(251,0:nfrx)
      dimension vcoul(251,0:nfrx)
      dimension edens(nrptx,0:nphx)
      dimension vclap(nrptx,0:nphx)
      dimension vtot (nrptx,0:nphx)
      dimension imt(0:nphx)
      dimension inrm(0:nphx)

      character*30 fname

c     note units --
c     potentials in rydbergs, so that v * 13.6 -> eV
c     density in #/(bohr)**3, so rho * e / (.529)**3 -> e/(Ang)**3

      do 180  iph = 0, nph
         ifr = ifrph(iph)
c        prepare file for unique potential data
         write(fname,172)  iph
  172    format('pot', i2.2, '.dat')
         open (unit=1, file=fname, status='unknown', iostat=ios)
         call chopen (ios, fname, 'wpot')
         call wthead(1)
         write(1,173)  iph, imt(iph), inrm(iph)
  173    format (1x, 3i4, '  Unique potential, I_mt, I_norman.',
     1          '    Following data in atomic units.')
         write(1,*) ' ifr ', ifr
         write(1,174)
  174    format ('   i      r         vcoul        rho',
     1           '     ovrlp vcoul  ovrlp vtot  ovrlp rho')
         do 178  i = 1, nrptx
            write(1,176) i, rr(i), vcoul(i,ifr), rho(i,ifr)/(4*pi),
     1                vclap(i,iph), vtot(i,iph), edens(i,iph)/(4*pi)
  176       format (1x, i3, 1p, 6e12.4)
  178    continue
         close(unit=1)
  180 continue

      return
      end
