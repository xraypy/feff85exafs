      subroutine istval (vtot, rholap, rmt, imt, rws, iws, vint, rhoint,
     1                   ierr)

c     This subroutine calculates interstitial values of v and rho
c     for an overlapped atom.  Inputs are everything except vint and
c     rhoint, which are returned.  vtot includes ground state xc.
c     rhoint is form density*4*pi, same as rholap
c
c     ierr = 0, normal exit
c          =-1, rmt=rws, no calculation possible

      implicit double precision (a-h, o-z)

      include 'dim.h'
      parameter (delta = 0.050 000 000 000 000)

      dimension vtot (nrptx)
      dimension rholap (nrptx)

c     Integrations are done in x (r = exp(x), see Louck's grid)
c     Trapezoidal rule, end caps use linear interpolation.
c     imt is grid point immediately below rmt, etc.
c     We will integrate over spherical shell and divide by volume of
c     shell, so leave out factor 4pi, vol = r**3/3, not 4pi*r**3/3,
c     similarly leave out 4pi in integration.

c     If rmt and rws are the same, cannot contribute to interstitial
c     stuff, set error flag
      vol = (rws**3 - rmt**3) / 3
      if (vol .le. 0)  then
         ierr = -1
         return
      endif
      ierr = 0

c     Calculation of vint including exchange correlation
c     Trapezoidal rule from imt+1 to iws
      vint = 0
      do 100  i = imt, iws-1
         fr = rr(i+1)**3 * vtot(i+1)
         fl = rr(i)**3   * vtot(i)
         vint = vint + (fr+fl)*delta/2
  100 continue
c     End cap at rws (rr(iws) to rws)
      xws = log (rws)
      xiws = xx(iws)
      g = xws - xiws
      fr = rr(iws+1)**3 * vtot(iws+1)
      fl = rr(iws)**3   * vtot(iws)
      vint = vint + (g/2) * ( (2-(g/delta))*fl + (g/delta)*fr)
c     End cap at rmt (rmt to rr(imt+1))
      xmt = log (rmt)
      ximt = xx(imt)
      g = xmt - ximt
      fr = rr(imt+1)**3 * vtot(imt+1)
      fl = rr(imt)**3   * vtot(imt)
      vint = vint - (g/2) * ( (2-(g/delta))*fl + (g/delta)*fr)
      vint = vint / vol

c     Calculation of rhoint
c     Trapezoidal rule from imt+1 to iws
      rhoint = 0
      do 200  i = imt, iws-1
         fr = rr(i+1)**3 * rholap(i+1)
         fl = rr(i)**3   * rholap(i)
         rhoint = rhoint + (fr+fl)*delta/2
  200 continue
c     End cap at rws (rr(iws) to rws)
      xws = log (rws)
      xiws = xx(iws)
      g = xws - xiws
      fr = rr(iws+1)**3 * rholap(iws+1)
      fl = rr(iws)**3   * rholap(iws)
      rhoint = rhoint + (g/2) * ( (2-(g/delta))*fl + (g/delta)*fr)
c     End cap at rmt (rmt to rr(imt+1))
      xmt = log (rmt)
      ximt = xx(imt)
      g = xmt - ximt
      fr = rr(imt+1)**3 * rholap(imt+1)
      fl = rr(imt)**3   * rholap(imt)
      rhoint = rhoint - (g/2) * ( (2-(g/delta))*fl + (g/delta)*fr)
      rhoint = rhoint / vol

      return
      end
