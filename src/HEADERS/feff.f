      program feff
      integer iabs, nabs
      logical ceels !KJ added ceels 5-06

      call rdinp(nabs,ceels)
      call ffmod1
      call ffmod2
      do 900 iabs = 1, nabs
        if (nabs.gt.1) call ffsort (iabs,ceels) !KJ 5-6 added second argument
        call ffmod4
        call ffmod5
        call ffmod6(iabs)
        call ffmod9
 900  continue
      stop
      end
