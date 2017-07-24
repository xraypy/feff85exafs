c     Note that leg nleg is the leg ending at the central atom, so that
c     ipot(nleg) is central atom potential, rat(nleg) position of
c     central atom.
c     Central atom has ipot=0
c     For later convience, rat(,0) and ipot(0) refer to the central
c     atom, and are the same as rat(,nleg), ipot(nleg).

c     text and title arrays include carriage control
      character*80 text, title
      character*6  potlbl

c common /str/:
c text header from potph
c title from paths.dat
c potential labels for output

      common /str/ text(40), title(5), potlbl(0:npotx)

      complex*16 ph, eref


c     common /pdata/:
c ph(nex,ltot+1,0:npotx),	!complex phase shifts,
c                          !central atom ipot=0
c rat(3,0:legtot+1),		!position of each atom, code units(bohr)
c eref(nex),		!complex energy reference
c em(nex),		!energy mesh
c ri(legtot), beta(legtot+1), eta(0:legtot+1), !r, beta, eta for each leg
c deg, rnrmav, xmu, edge,	!(output only)
c lmax(nex,0:npotx),	!max l with non-zero phase for each energy
c ipot(0:legtot),	!potential for each atom in path
c iz(0:npotx),	!atomic number (output only)
c ltext(40), ltitle(5),	!length of each string
c nsc, nleg,	!nscatters, nlegs (nleg = nsc+1)
c npot, ne,	!number of potentials, energy points
c ik0,		!index of energy grid corresponding to k=0 (edge)
c ipath, 	!index of current path (output only)
c ihole,	!(output only)
c l0, il0,	!lfinal and lfinal+1 (used for indices)
c lmaxp1,	!largest lmax in problem + 1
c ntext, ntitle	!number of text and title lines

      common /pdata/ ph(nex,ltot+1,0:npotx), rat(3,0:legtot+1),
     $     eref(nex), em(nex), ri(legtot), beta(legtot+1),
     $     eta(0:legtot+1), deg, rnrmav, xmu, edge, lmax(nex,0:npotx),
     $     ipot(0:legtot), iz(0:npotx), ltext(40), ltitle(5), nsc, nleg,
     $     npot, ne, ik0, ipath, ihole, l0, il0, lmaxp1, ntext, ntitle
