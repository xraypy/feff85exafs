      subroutine afolp (verbse, nph, nat, iphat, rat, iatph, xnatph,
     1                novr, iphovr, nnovr, rovr, folp, folpx, iafolp,
     1                edens, edenvl,
     2                dmag, vclap, vtot, vvalgs, imt, inrm, rmt, rnrm, 
     2                ixc, rhoint, vint, rs, xf, xmu, xmunew,
     3                rnrmav, qtotel, inters, totvol)

c     find folp(iph) automatically and recalculates
c     interstitial parameters, rmt, vint, etc.
c     written by ala 11.97
      implicit double precision (a-h, o-z)

      include '../HEADERS/const.h'
      include '../HEADERS/dim.h'

      logical verbse

      dimension iphat(natx)
      dimension rat(3,natx)
      dimension iatph(0:nphx)
      dimension xnatph(0:nphx)
      dimension novr(0:nphx)
      dimension iphovr(novrx,0:nphx)
      dimension nnovr(novrx,0:nphx)
      dimension rovr(novrx,0:nphx)
      dimension folp(0:nphx), folpx(0:nphx)
      dimension edens(251,0:nphx), edenvl(251,0:nphx)
      dimension dmag(251,0:nphx+1)
      dimension vclap(251,0:nphx)
      dimension vtot (251,0:nphx), vvalgs (251,0:nphx)
      dimension imt(0:nphx)
      dimension inrm(0:nphx)
      dimension rmt(0:nphx), rmtx(0:nphx)
      dimension rnrm(0:nphx)
      character*512 slog

      do 5 iph=0,nph
         rmtx(iph) = rmt(iph) / folp(iph)
   5  continue

      if (verbse)
     1call wlog('    : ipot, Norman radius, Muffin tin radius, Overlap')
      if (iafolp.ge.0) then
         do 400  iph = 0, nph
c          old algorithm for automatic overlap
c          folp(iph) = 1 + 0.7*(rnrm(iph)/rmt(iph) - 1)
           folp(iph) = folpx(iph)
           rmt(iph) = folp(iph) * rmtx(iph)

           if (verbse) then
 398          format(i10, 1p, 3e13.5)
              write(slog,398) iph, rnrm(iph)*bohr, rmt(iph)*bohr,
     1               folp(iph)
              call wlog(slog)
           endif
  400    continue

         idmag = 0
         call istprm (nph, nat, iphat, rat, iatph, xnatph,
     1               novr, iphovr, nnovr, rovr, folp, folpx, iafolp,
     1               edens, edenvl, idmag,
     2               dmag, vclap, vtot, vvalgs, imt, inrm, rmt, rnrm,
     3               ixc, rhoint,vint, rs, xf, xmu, xmunew,
     4               rnrmav, qtotel, inters, totvol)

      endif

      return
      end
