c -*- fortran -*-
c     common /lambda/  
c    4   mlam(lamtot), 	!mu for each lambda
c    5   nlam(lamtot),	!nu for each lambda
c    1   lamx, 		!max lambda in problem
c    2   laml0x, 	!max lambda for vectors involving absorbing atom
c    3   mmaxp1, nmax 	!max mu in problem + 1, max nu in problem
      common /lambda/ mlam(lamtot), nlam(lamtot), lamx, laml0x,
     1                mmaxp1, nmax
