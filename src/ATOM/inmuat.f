      subroutine inmuat (ihole, xionin, iunf, xnval, iholep, xmag, iorb)
      implicit double precision (a-h,o-z)
      dimension xnval(30), xmag(30), iorb(-4:3)
      common/itescf/testy,rap(2),teste,nz,norb,norbsc
c the meaning of common variables is described below
      common/ratom1/xnel(30),en(30),scc(30),scw(30),sce(30),
     1nq(30),kap(30),nmax(30)
c en one-electron energies
c scc factors for acceleration of convergence
c scw precisions of wave functions
c sce precisions of one-electron energies
c nmax number of tabulation points for orbitals
      common/scrhf1/eps(435),nre(30),ipl
c eps non diagonal lagrange parameters
c nre distingue: - the shell is closed (nre <0)
c                  the shell is open (nre>0)
c                - the orbitals in the integral rk if abs(nre) > or =2
c ipl define the existence of lagrange parameters (ipl>0)
      common/snoyau/dvn(251),anoy(10),nuc
c dvn nuclear potential
c anoy development coefficients at the origin of nuclear potential
c this development is supposed to be written anoy(i)*r**(i-1)
c nuc index of nuclear radius (nuc=1 for point charge)
      common/tabtes/hx,dr(251),test1,test2,ndor,np,nes,method,idim
      data nucm/11/,nesn/50/,ideps/435/

      ndor=10
      testy=1.0d-05
c testy precision for the wave functions
      teste=5.0d-06
c teste precision for the one-electron energies
      rap(1)=100.
      rap(2)=10.
c rap tests of precision for soldir
      do 10 i = 1, 30
        en(i) = 0.d0
        xmag(i) = 0
  10  xnval(i) = 0

      call getorb (nz, ihole, xionin, iunf, norb, norbsc, iorb,
     1             iholep, nq, kap, xnel, xnval, xmag)
      xk=0
      do 411 i=1,norb
 411  xk=xk+xnel(i)
      if ( abs(nz-xionin-xk) .gt. 0.001) then
         call par_stop('check number of electrons in getorb.f')
c        stop
      endif
      norbsc=norb
c nz atomic number     noi ionicity (nz-number of electrons)
c norb number of orbitals
c xnel(i) number of electrons on orbital i.
c first norbsc orbitals will be determined selfconsistently,
c the rest of orbitals are orthogonolized if iorth is non null,
c and their energies are those on cards if iene is non null
c or otherwise are the values obtained from solving dirac equation
      nes=nesn
c nes number of attempts in program soldir
      nuc=nucm
c nuc number of points inside nucleus (11 by default)
      do 171 i=1,ideps
 171  eps(i)=0.0d 00

      idim = 251
      if (mod(idim,2) .eq. 0) idim=idim-1

      ipl=0
c if ipl non null, it permits a repartition of tabulation points
c and certain precision tests.
      do 401 i=1,norb
         nre(i)=-1
         llq= abs(kap(i))
         l=llq+llq
         if (kap(i).lt.0) llq=llq-1
         if (llq.lt.0.or.llq.ge.nq(i).or.llq.gt.3) then
            call par_stop('kappa out of range, check getorb.f')
c           stop
         endif
         nmax(i)=idim
         scc(i)=0.3
         if (xnel(i) .lt. l)  nre(i)=1
         if (xnel(i) .lt. 0.5)  scc(i)=1.0
         do 385 j=1,i-1
            if (kap(j).ne.kap(i)) go to 385
            if (nre(j).gt.0.or.nre(i).gt.0) ipl=ipl+1
 385     continue
 401  continue
      return
      end
