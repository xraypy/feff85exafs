c  common /lambda/  :
c     mlam(lamtot), 	!mu for each lambda
c     nlam(lamtot),	!nu for each lambda
c     lamx, 		!max lambda in problem
c     laml0x, 	!max lambda for vectors involving absorbing atom
c     mmaxp1, nmax 	!max mu in problem + 1, max nu in problem

      common /lambda/ mlam(lamtot), nlam(lamtot),
     $ lamx, laml0x, mmaxp1, nmax
