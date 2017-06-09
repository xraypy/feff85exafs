c***********************************************************************
c
c     this subroutine calculates the ' energy dependent
c     exchange-correlation potential' (or 'dirac- hara potential')
c     ref.: paper by s.h.chou, j.j.rehr, e.a.stern, e.r.davidson (1986)
c
c     inputs:    rs in a.u.
c                xk momentum in a.u.
c                vi0 constant imaginary part in rydbergs
c     outputs:   vr --- dirac potential (Hartrees)
c                vi --- constant imag part of the potential (Hartrees)
c     written by j. mustre 8/31/87
c**********************************************************************

      subroutine edp (rs, xk, vi0, vr, vi)
      implicit double precision (a-h, o-z)
      include 'const.h'

      xf = fa / rs

c     p = sqrt (k^2 + kf^2) is the local momentum, and x = p / kf
c     Reference formula 23 in Role of Inelastic effects in EXAFS
c     by Rehr and Chou. EXAFS1 conference editted by Bianconi.
c     x is local momentum in units of fermi momentum

      x = xk / xf
      x = x + 1.0e-5
c     set to fermi level if below fermi level
      if (x .lt. 1.00001) x = 1.00001
      c = abs( (1+x) / (1-x) )
      c = log(c)
      vr = - (xf/pi) * (1 + c * (1-x**2) / (2*x))

c     Note vi=vi0/2 to have both real and imaginary part in hartrees
c     to be consistent with  other subroutines.
      vi = vi0 / 2
      return
      end
