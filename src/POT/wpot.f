      subroutine wpot (nph, edens, imt, inrm,
     1                 rho, vclap, vcoul, vtot, ntitle, title)

c     Writes potentials to file name POTxx.DAT for each unique pot.

      implicit double precision (a-h, o-z)

      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'

      dimension rho(251,0:nphx+1)
      dimension vcoul(251,0:nphx+1)
      dimension edens(251,0:nphx)
      dimension vclap(251,0:nphx)
      dimension vtot (251,0:nphx)
      dimension imt(0:nphx)
      dimension inrm(0:nphx)
      character*80 title(ntitle)

      character*30 fname
c#mn
       external rr

c     note units --
c     potentials in hartrees, so that v * 27.2 -> eV
c     density in #/(bohr)**3, so rho * e / (.529)**3 -> e/(Ang)**3

      do 180  iph = 0, nph
c        prepare file for unique potential data
         write(fname,172)  iph
  172    format('pot', i2.2, '.dat')
         open (unit=1, file=fname, status='unknown', iostat=ios)
         call chopen (ios, fname, 'wpot')
         call wthead(1, ntitle, title)
         write(1,173)  iph, imt(iph), inrm(iph)
  173    format (1x, 3i4, '  Unique potential, I_mt, I_norman.',
     1          '    Following data in atomic units.')
         write(1,*) ' iph ', iph
         write(1,174)
  174    format ('   i      r         vcoul        rho',
     1           '     ovrlp vcoul  ovrlp vtot  ovrlp rho')
c        need some limit here, 1250 points is silly.  Use
c        r <= 38, which gives 249 points with usual rgrid
         do 178  i = 1, 251
            if (rr(i) .gt. 38)  goto 179
            write(1,176) i, rr(i), vcoul(i,iph), rho(i,iph)/(4*pi),
     1                vclap(i,iph), vtot(i,iph), edens(i,iph)/(4*pi)
  176       format (1x, i4, 1p, 6e12.4)
  178    continue
  179    continue
         close(unit=1)
  180 continue

      return
      end
