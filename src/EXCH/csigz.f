      SUBROUTINE CSigZ(Energy, Mu, Rs, ReSig, ImSig, ZTot, WpScl,
     &     AmpFac)
c      SUBROUTINE CSigZ(Energy, Mu, Rs, ReSig, ImSig, ZTot, WpScl, Gamma,
c     &     AmpFac)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Written by Josh Kas
c     This subroutine calculates the self energy Sigma(k(Energy),Energy)
c     based on an electron gas model of epsilon^-1.
c      
c     Solve: k0**2 = 2*Energy - 2*Mu -
c                    2*(Sigma(k0,Energy)-Sigma(kFermi,EFermi))
c            
c     Steps:
c
c            1. k0**2  = 2*(Energy-Mu) + SigmaF (SigmaF is self energy at fermi level).
c            2. Sigma0 = Sigma(k0,Energy)
c            3. Find derivative w.r.t. E dSgdE
c            4. k1**2  = 
c                  k0**2 - 2*(Sigma0-SigmaF)/(1-dSgdE)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      include '../HEADERS/const.h'
c     Parameters:
c     MxPole - Maximum number of poles
      INTEGER MxPole
      PARAMETER(MxPole=1000)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Input:
c     Energy - Energy at which to evaluate Sigma
c     Mu     - Fermi energy as calculated by FEFF.
c     Rs     - R sub s (sphere of radius Rs holds charge e)
c     WpScl  - Scale Wp in interstitial by WpScl
c     Gamma  - Use broadening Gamma when calculating Sigma
c     AmpFac - Use amplitude AmpFac for plasmon pole.
      DOUBLE PRECISION Rs, WpScl(MxPole), AmpFac(MxPole), Mu
c      DOUBLE PRECISION Gamma(MxPole),
      COMPLEX*16 Energy
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Output:
c     ReSig  - Re[Sigma(k,e)]
c     ImSig  - Im[Sigma(k,e)]
c     ZTot   - Renormalization factor Z = 1/(1-dSgdE)
c     Note: Atomic units are used.
      DOUBLE PRECISION ReSig, ImSig
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Local Variables:
c     kFermi - Fermi momentum calculated from Rs via electron gas
c              approximation.
c     EFermi - Fermi energy = kFermi^2/2
c     Wp     - Electron gas plasmon frequency.
c     Gam    - broadening for broadened pole.
c     ckF    - complex variable to store kFermi
c     ck0    - complex momentum
c     SigmaF - Self energy at the fermi energy and fermi momentum
c              Does not include the Hartree Fock part.
c     Sigma0 - Single pole self energy evaluated at
c              Wp = Wp(ipole)/Wp(R_Interstitial)*Wp(Rs)
c     dSgdE  - Derivative of sigma w.r.t. Energy
c     ZTot   - Renormalization factor Z = 1/(1-dSgdE)
c     RelEn  - Energy relative to the fermi energy from FEFF      
c              Relen = Energy - Mu + EFermi
c     SigTot - Total many pole deltaSigma = Sigma(E,k(E))-Sigma(EF,kF)
c     DelHF  - Sigma_HartreeFock(k) - Sigma_HartreeFock(kF)
      DOUBLE PRECISION kFermi, EFermi, Wp, Gam
      COMPLEX*16 ckF, ck0, SigmaF, Sigma0, dSgdE, ZTot, RelEn, RelEnP,
     &     SigTot, DelHF, SigmaP
c     Loop variables
      INTEGER i1, i2

c     Parameters:
      DOUBLE PRECISION DPZero
      PARAMETER(DPZero = 0.d0)
      INTEGER MxIter
      PARAMETER(MxIter = 1)

c     Externals:
c     Sigma1 - calculates the energy dependent part of self energy
c              for a single pole.
c     dSigma - calculates derivative of self energy w.r.t energy.
c     HFExt  - calculates Hartree Fock exchange
      COMPLEX*16 Sigma1, dSigma, HFExc
      EXTERNAL Sigma1, dSigma, HFExc
      
c     Initialization
      ZTot = 0.d0
      kFermi = fa/Rs
      EFermi = kFermi*kFermi/2.d0
      SigTot=0.d0
      dSgdE = 0.d0
      SigmaF = 0.d0
      Gam = 0.d0

c     Loop1: Start self consistency loop.
c     This does not seem to work, so MxIter = 1
      DO i2 = 1, MxIter

c        Loop2: Loop over poles to get SigmaF
         DO i1 = 1, MxPole
            IF(WpScl(i1).lt.0) GOTO 5
            
c           Wp is in Hartrees            
            Wp = SQRT(3.d0/Rs**3)*WpScl(i1)
c            Gam = Gamma(i1)
            
c           find Sigma_Fermi (SigmaF)
            ckF = kFermi*1.00001d0
            RelEn = EFermi
            SigmaF = SigmaF + Sigma1(ckF,RelEn,Wp,Gam,AmpFac(i1),
     &           kFermi, EFermi)
         END DO
c        End Loop2
 5       CONTINUE
         
         dsgdE = 0.d0
c        Loop3: Loop over poles         
         DO i1 = 1, MxPole
            IF(WpScl(i1).lt.0) GOTO 10
c           Wp is in Hartrees
            Wp = SQRT(3.d0/Rs**3)*WpScl(i1)
c            Gam = Gamma(i1)
c           Start with ck0=Sqrt[Re(Energy)-Mu+EFermi]
            RelEn = DBLE(Energy) - Mu + EFermi
            ck0 = SQRT(2.d0*DBLE(RelEn))
            
c           Find Sigma0 = Sigma(ck0,E); ck0=SQRT(2*(Energy-Mu))
            Sigma0 = Sigma1(ck0,RelEn,Wp,Gam,AmpFac(i1),kFermi,
     &           EFermi)
                  
            RelEnP = RelEn*0.001d0
            SigmaP = Sigma1(ck0,RelEnP,Wp,Gam,AmpFac(i1),kFermi,
     &           EFermi)

c            write(71,*) DBLE(RelEn-EFermi), DBLE(dSgdE), DIMAG(dSgDE)
            dSgdE = dSgdE + (SigmaP - Sigma0)/(RelEnP-RelEn)
c            dSgdE = dSgdE + 
c     &          dSigma(ck0,RelEn,Wp,Gam,AmpFac(i1),kFermi,EFermi)
c           Uncomment the following line to print derivative            
c            write(72,*) DBLE(RelEn-EFermi), DBLE(dSgdE2), DIMAG(dSgDE2)         

c           SigTot is sum of poles
            SigTot = SigTot + Sigma0
            
c        End Loop3: loop over poles.            
         END DO         
 10      CONTINUE
         
c     End Loop1: self-consistency loop      
      END DO


c     Add Hartree Fock part of delta Sigma.
      DelHF = HFExc(ck0,EFermi,kFermi) - HFExc(ckF,EFermi,kFermi)
      SigTot = SigTot - SigmaF + DelHF

c     Form ZTot and return Re and Im parts of Sigma.
      ZTot = 1.d0/(1.d0-dSgdE)
c      ZTot = 1.d0
      SigTot = ZTot*(SigTot)
      
      ReSig = DBLE(SigTot)
      ImSig = DIMAG(SigTot)
      
      RETURN
      END

