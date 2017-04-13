      subroutine rphbin (in)
      implicit double precision (a-h, o-z)

c     Reads input from unit in.  Returns (via /pdata/)
c       energy mesh (ne, em and eref),
c       ph (npot, lmax, lmaxp1, ph),
c       final state (l0, il0)
c
c     phmin is min value to use for |phase shift|

      include 'dim.h'
      include 'pdata.h'

      parameter (phmin = 1.0e-8)

c     These header lines do not include carriage control
      read(in) ntext
      do 62  i = 1, ntext
         read(in) text(i)
         read(in) ltext(i)
   62 continue
      read(in) ne, npot, ihole, rnrmav, xmu, edge, ik0
      read(in) (em(ie),ie=1,ne)
      read(in) (eref(ie),ie=1,ne)
      lmaxp1 = 0
      do 80  iph = 0, npot
         read(in) lmax0, iz(iph)
         read(in) potlbl(iph)
         do 70  ie = 1, ne
            read(in)  (ph(ie,ll,iph), ll=1,lmax0+1)
            lmax(ie,iph) = 0
c           Set lmax to include only non-zero phases
            do 60  il = 1, lmax0+1
               if (abs(ph(ie,il,iph)) .lt. phmin)  goto 61
               lmax(ie,iph) = il-1
   60       continue
   61       continue
            if (lmax(ie,iph)+1 .gt. lmaxp1)  lmaxp1 = lmax(ie,iph)+1
   70    continue
   80 continue

c-----l0 is angular momentum of final state
c     Selection rule says that final state has angmom = l_init+1
c     ihole  initial state from ihole         final state
c     1      K    1s      L=0 -> linit=0   L0=1 -> lfinal=1
c     2      LI   2s      L=0 -> linit=0   L0=1 -> lfinal=1
c     3      LII  2p 1/2  L=1 -> linit=1   L0=2 -> lfinal=2
c     4      LIII 2p 3/2  L=1 -> linit=1   L0=2 -> lfinal=2
c     5+     M -- think about this later...
      if (ihole .le. 2)  then
c        hole in s state (1s or 2s)
         linit = 0
         lfinal = 1
      elseif (ihole .le. 4)  then
c        hole in p state (2p 1/2  or  2p 3/2)
         linit = 1
         lfinal = 2
      else
c        some m hole, n=3, could go to d state
         call echo(' Feff6L only does K and L shells')
         call fstop(' at RPHBIN')
      endif
      l0 = lfinal
      il0 = l0 + 1

      return
      end
