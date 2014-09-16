c -*- fortran -*-
c     Note that leg nleg is the leg ending at the central atom, so that
c     ipot(nleg) is central atom potential, rat(nleg) position of 
c     central atom.
c     Central atom has ipot=0
c     For later convience, rat(,0) and ipot(0) refer to the central
c     atom, and are the same as rat(,nleg), ipot(nleg).

c     text arrays 
      character*80 text 
      character*6  potlbl
c     text from paths.dat, potential labels
      common /str/ text (5), potlbl(0:nphx)

      complex*16 ph, eref, em
      double precision rnrmav, xmu, edge
c     common /pdata/ ph(nex,-ltot:ltot,0:nphx), !complex phase shifts ipot=0
c    1 eref(nex),		!complex energy reference
c    1 rat(3,0:legtot+1),	!position of each atom, code units(bohr)
c    1 em(nex),		!energy mesh
c    1 ri(legtot), beta(legtot+1), eta(0:legtot+1), !r, beta, eta for each leg
c    1 deg, rnrmav, xmu, edge,	!(output only)
c    1 lmax(nex,0:nphx),	!max l with non-zero phase for each energy
c    1 ipot(0:legtot),	!potential for each atom in path
c    1 iz(0:nphx),	!atomic number (output only)
c    1 ltext (5),	!length of each string
c    1 nsc, nleg,	!nscatters, nlegs (nleg = nsc+1)
c    1 npot, ne,	!number of potentials, energy points
c    1 ik0,		!index of energy grid corresponding to k=0 (edge)
c    1 ipath, ihole, 	!index of current path  and hole (output only)
c    1 kinit, linit, ilinit,  ! initial state kappa and ang. mom.
c    1 lmaxp1,	!largest lmax in problem + 1
c    1 ntext 	!number of text  lines

      common /pdata/ ph(nex,-ltot:ltot,0:nphx), eref(nex), em(nex),
     1 rat(3,0:legtot+1),  
     1 ri(legtot), beta(legtot+1), eta(0:legtot+1),
     1 deg, rnrmav, xmu, edge, lmax(nex,0:nphx), ipot(0:legtot),
     1 iz(0:nphx), ltext (5),
     1 nsc, nleg, npot, ne, ik0, ipath, ihole, 
     1 kinit, linit, ilinit,
     1 lmaxp1, ntext 
