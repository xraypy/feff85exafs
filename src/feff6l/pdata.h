c     Note that leg nleg is the leg ending at the central atom, so that
c     ipot(nleg) is central atom potential, rat(nleg) position of 
c     central atom.
c     Central atom has ipot=0
c     For later convience, rat(,0) and ipot(0) refer to the central
c     atom, and are the same as rat(,nleg), ipot(nleg).

c     text and title arrays include carriage control
      character*80 text, title
      character*6  potlbl
      common /str/ text(40),	!text header from potph
     1             title(5),	!title from paths.dat
     1             potlbl(0:npotx)	! potential labels for output

      complex*16 ph, eref
      common /pdata/
     1 ph(nex,ltot+1,0:npotx),	!complex phase shifts,
     1					!central atom ipot=0
     1 rat(3,0:legtot+1),		!position of each atom, code units(bohr)
     1 eref(nex),		!complex energy reference
     1 em(nex),		!energy mesh
     1 ri(legtot), beta(legtot+1), eta(0:legtot+1), !r, beta, eta for each leg
     1 deg, rnrmav, xmu, edge,	!(output only)
     1 lmax(nex,0:npotx),	!max l with non-zero phase for each energy
     1 ipot(0:legtot),	!potential for each atom in path
     1 iz(0:npotx),	!atomic number (output only)
     1 ltext(40), ltitle(5),	!length of each string
     1 nsc, nleg,	!nscatters, nlegs (nleg = nsc+1)
     1 npot, ne,	!number of potentials, energy points
     1 ik0,		!index of energy grid corresponding to k=0 (edge)
     1 ipath, 	!index of current path (output only)
     1 ihole,	!(output only)
     1 l0, il0,	!lfinal and lfinal+1 (used for indices)
     1 lmaxp1,	!largest lmax in problem + 1
     1 ntext, ntitle	!number of text and title lines
