      SUBROUTINE CSigma(Energy, Mu, Rs, ReSig, ImSig, WpScl,!
     &     AmpFac)
c      SUBROUTINE CSigma(Energy, Mu, Rs, ReSig, ImSig, WpScl, Gamma,
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
c            3. k1**2  = 
c                  k0**2 - 2*(Sigma0-SigmaF)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      include '../HEADERS/const.h'
c     Parameters:
c     MxPole - Maximum number of poles      
      INTEGER MxPole
      PARAMETER(MxPole=1000)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Input:
c     Energy - Energy at which to evaluate Sigma
c     Mu     - Fermi energy as calculated by FEFF
c     Rs     - R sub s (sphere of radius Rs holds charge e)
c     WpScl  - Scale Wp in interstitial by WpScl
c     Gamma  - Use broadening Gamma when calculating Sigma
c     AmpFac - Use amplitude AmpFac for plasmon pole.
c     Note: Atomic units are used.
      DOUBLE PRECISION Rs, WpScl(MxPole), AmpFac(MxPole), Mu
c      DOUBLE PRECISION Gamma(MxPole)
      COMPLEX*16 Energy
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Output:
c     ReSig  - Re[Sigma(Energy,k(Energy))]
c     ImSig  - Im[Sigma(Energy,k(Energy))]
      DOUBLE PRECISION ReSig, ImSig
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
      COMPLEX*16 ckF, ck0, SigmaF, Sigma0,
     &     RelEn, SigTot, DelHF
      INTEGER i1, i2

c     Parameters:
      DOUBLE PRECISION DPZero, h
      PARAMETER(DPZero = 0.d0, h = 1.d-3)
      INTEGER MxIter

      PARAMETER(MxIter = 1)

c     Externals:
      COMPLEX*16 Sigma1, dSigma, HFExc
      EXTERNAL Sigma1, dSigma, HFExc
      
c     Initialization
      kFermi = fa/Rs
      EFermi = kFermi*kFermi/2.d0
      SigTot=0.d0
      SigmaF = 0.d0 
      Gam = 0.d0
      ckF = 0.d0

c     Loop1: Start self consistency loop.
      DO i2 = 1, MxIter
c        Loop2: Loop over poles to find SigmaF
         DO i1 = 1, MxPole
            IF(WpScl(i1).lt.-1000.d0) GOTO 5
            
c           Wp is in Hartrees
            Wp = SQRT(3.d0/rs**3)*WpScl(i1)
            
c           find Sigma_Fermi (SigmaF)
            ckF = kFermi*1.00001d0
            RelEn = EFermi
            SigmaF = SigmaF + Sigma1(ckF,RelEn,Wp,Gam,AmpFac(i1),
     &           kFermi,EFermi)
         END DO
 5       CONTINUE
         
c        Loop3: Loop over poles
         DO i1 = 1, MxPole
            IF(WpScl(i1).lt.-1000.d0) GOTO 10
c           Wp is in Hartrees
            Wp = SQRT(3.d0/rs**3)*WpScl(i1)

c           Start with ck0=Sqrt[Re(Energy)-Mu+EFermi]
            RelEn = DBLE(Energy) - Mu + EFermi
            ck0 = SQRT(2.d0*DBLE(RelEn))
            
c           Find Sigma0
            Sigma0 = Sigma1(ck0,RelEn,Wp,Gam,AmpFac(i1),kFermi,
     &           EFermi)
            
            SigTot = SigTot + Sigma0
            
c        End loop over poles.            
         END DO
 10      CONTINUE
         
c     End self-consistency loop
      END DO

c     Form delta sigma and retur.n
      SigTot = SigTot - SigmaF
      DelHF = HFExc(ck0,EFermi,kFermi) - HFExc(ckF,EFermi,kFermi)
      SigTot = SigTot + DelHF
      
      ReSig = DBLE(SigTot)
      ImSig = DIMAG(SigTot)
      
      RETURN
      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      FUNCTION Sigma1(ck,Energy,Wi,Gamma,Amp,kFermi,EFermi)
c     Written by Josh Kas
c     Function Sigma calculates the energy dependent part
c     of Sigma(ck,Energy) from H.L. electron gas model.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      include '../HEADERS/const.h'
c     Input:
c     ck     - complex momentum
c     Energy - Energy
c     Wi     - Plasmon pole energy
c     Gamma  - Broadening of plasmon pole.
c     Amp    - Amplitude of plasmon pole.
c     kFermi - Fermi momentum
c     EFermi - Fermi energy
c              This is used when calculating dSigma/dE
      DOUBLE PRECISION Wi, Gamma, Amp, kFermi, EFermi
      COMPLEX*16 ck, Energy
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Output:
c     Sigma(ck,Energy)
      COMPLEX*16 Sigma1
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Local Variables:
c     NSing  - Number of singularities of integrand (used in cgratr)
c     NCalls - Number of fcn calls used in cgratr
c     MaxNR  - Max number of regions used in cgratr
c     DPPar  - Array of double precision parameters passed to cgratr
c              to be used in the functions r1, r2, and r3
c     CPar   - Array of complex parameters passed to cgratr
c              to be used in the functions r1, r2, and r3
c     Limit1 - Lower limit of integration
c     Limit2 - Upper limit of integration
c     HLInt1 - Integral of r2 (first integral in eq. 13 of H.L.)
c     HLInt2 - Integral of r1 (second integral in eq. 13 of H.L.)
c     HLInt3 - Integral of r1 or r3 (3rd or 4th integral)
c     XSing  - Array of singularities of integrand (used by cgratr)      
      INTEGER NSing, NCalls, MaxNR
      DOUBLE PRECISION DPPar(10), Beta
      COMPLEX*16 CPar(10), Limit1, Limit2, HLInt1, HLInt2, HLInt3,
     &     XSing(20)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Loop variables:
      INTEGER i1, i2
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Parameters:
c     ZeroPl - lower bound Limit1
c     Inf    - Upper bound of Limit2
c     AbsErr - absolute error used by cgratr
c     RelErr - Relative error used by cgratr
c     Error  - used for error codes by cgratr      
      DOUBLE PRECISION ZeroPl, Inf, AbsErr, RelErr, Error
      PARAMETER(ZeroPl = 1.d-5, Inf = 1.d2)
      PARAMETER(AbsErr = 1.d-5, RelErr = 1.d-4)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Externals:
c     Externals:
c     cgratr      - integration routine
c     dr1,dr2,dr3 - functions to integrate
c     HFExc       - Calculates Hartree Fock exchange      
      COMPLEX*16 cgratr, r1, r2, r3, HFExc
      EXTERNAL cgratr, r1, r2, r3, HFExc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c     Initialization:
      NSing  = 0
      NCalls = 0
      MaxNR  = 0
c     DPPar is array of dp parameters to evaluate functions in cgratr.
c     Everything is in dimensionless units.
c     1. xwg
      DPPar(1) = (Wi/EFermi)
c     2. xgam
      DPPar(2) = gamma/EFermi
c     3. xe
      DPPar(3) = DBLE(Energy)/EFermi
c     4. xeg (gap energy)
      DPPar(4) = 0.d0
c     CPar is array of complex parameters to evaluate functions in cgratr.
c     ck in dimensionless units.
      CPar(1) = ck/kFermi
c     2. ce (complex energy)
      CPar(2) = Energy/EFermi
      
c     Josh - This is a possible fix for functions that overlap zero by a large
c     amount so that Wp does not equal Wi. 
*      Beta = 1.d0*Gamma*SQRT(2.d0/(Wi**2+Gamma**2))
c      Beta  = 0.9d0*Gamma*Gamma/(Wi**2+Gamma**2)
c      Wp =  2.d0*Gamma*LOG(Gamma/Beta) + 2.d0*ATAN2(Wi,Gamma) -
c     &     Gamma*LOG(Gamma**2 + Wi**2)
c      Wp = SQRT(Wi*Wp)
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      DO i1 = 1, 1
         IF(i1.eq.2) THEN
            DPPar(1) = 0.d0
            DPPar(2) = 0.9d0*Gamma/EFermi
         END IF         
c     Calculate integrals in eq. 13 of H.L.
c     1)
c     Integral { ln[((ck+q)**2-E+Wq)/((ck-q)**2-E+Wq)] }
c     from ck + kFermi to Inf.
         Limit1 = ck/kFermi+1.d0
         Limit2 = Inf
      
c     Find singularities in r2
         iFcn = 2
         CALL FndSng(Limit1, Limit2, NSing, XSing, DPPar, CPar, iFcn)
      
c     Calculate integral
         HLInt1 = cgratr(r2,DPPar,CPar,Limit1,Limit2,AbsErr,RelErr,
     &        NSing, XSing,Error,NCalls,MaxNR)
         DO i2 = 1, NSing
            XSing(i2) = (1.d0,0)*XSing(i2)
         END DO

c     2)
c     Integral { ln[(kFermi**2-E-Wq)/(kFermi**2-E+Wq)*
c                     ((ck+q)**2-E+Wq)/((ck-q)**2-E-Wq)] }
c     From ck - kFermi to ck + kFermi
         Limit1 = MAX(ABS(DBLE(ck)/kFermi-1.d0), ZeroPl)
         Limit2 = ck/kFermi+1.d0

c     Find singularities in r1
         IFcn = 1
         CALL FndSng(Limit1, Limit2, NSing, XSing, DPPar, CPar, iFcn)

c     Calculate integral
         HLInt2 = cgratr(r1,DPPar,CPar,Limit1,Limit2,AbsErr,RelErr,
     &        NSing, XSing,Error,NCalls,MaxNR)
         DO i2 = 1, NSing
            XSing(i2) = (1.d0,0)*XSing(i2)
         END DO      
            
c     3)
c     Theta(kFermi-Re(ck)) *
c     Integral { ln[((ck+q)**2-E-Wq)/((ck-q)**2-E-Wq)] } +
c     Theta(Re(ck)-kFermi) *
c     Integral { ln[((ck+q)**2-E+Wq)/((ck-q)**2-E+Wq)] }
c     (Integrals from 0 to kFermi - k and 0 to k - kFermi)
         Limit1 = ZeroPl
         Limit2 = ABS(DBLE(ck)/kfermi-1.d0)

c     If ck = kFermi, HLInt3 = 0
         IF((ABS(DBLE(ck)-kFermi).lt.ZeroPl).or.
     &        (DBLE(Limit2).le.DBLE(Limit1))) THEN
            HLInt3 = 0.d0
c     If ck < kFermi, HLint3 = Integral { ln[((ck+q)**2-E-Wq)/((ck-q)**2-E-Wq)] }
         ELSEIF(DBLE(ck).lt.kFermi) THEN
            Limit2 = 1.d0 - DBLE(ck)/kFermi
            
c     Find singularities in r3
            iFcn = 3         
            CALL FndSng(Limit1, Limit2, NSing, XSing, DPPar, CPar, iFcn)
c     Calculate integral
            HLInt3 = cgratr(r3,DPPar,CPar,Limit1,Limit2,AbsErr,RelErr,
     &           NSing,XSing,Error,NCalls,MaxNR)
            DO i2 = 1, NSing
               XSing(i2) = (1.d0,0)*XSing(i2)
            END DO
         
c     Else ck > kFermi, HLint3 = Integral { ln[((ck+q)**2-E+Wq)/((ck-q)**2-E+Wq)] }
         ELSE
            Limit2 = DBLE(ck)/kFermi - 1.d0
            
c     Find singularities in r2
            iFcn = 2
            CALL FndSng(Limit1, Limit2, NSing, XSing, DPPar, CPar, iFcn)
            
c     Calculate integral
            HLInt3 = cgratr(r2,DPPar,CPar,Limit1,Limit2,AbsErr,RelErr,
     &           NSing,XSing,Error,NCalls,MaxNR)
            DO i2 = 1, NSing
               XSing(i2) = (1.d0,0)*XSing(i2)
            END DO
         END IF
                           
         IF(i1.eq.1) THEN
            Sigma1 = - Amp*Wi*(Wi-coni*Gamma)/(2.d0*pi*EFermi*ck)*
     &           (HLInt1 + HLInt2 + HLInt3)
         ELSE
            Sigma1 = Sigma1 - Amp*Beta*
     &           Wi*(coni*Gamma)/(pi*EFermi*ck)*
     &           (HLInt1 + HLInt2 + HLInt3)
         END IF
      END DO
      
      RETURN
      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      FUNCTION dSigma(ck,Energy,Wi,Gamma,Amp,kFermi,EFermi)
c     Written by Josh Kas
c     Function dSigma calculates dSigma(ck,Energy)/dE from H.L.
c     electron gas model.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      include '../HEADERS/const.h'
c     Input:
c     ck     - complex momentum
c     Energy - Energy
c     Wi     - Plasmon pole energy
c     Gamma  - Broadening of plasmon pole.
c     Amp    - Amplitude of plasmon pole.
c     kFermi - Fermi momentum
c     EFermi - Fermi energy
c              This is used when calculating dSigma/dE
      DOUBLE PRECISION Wi, Gamma, Amp, kFermi, EFermi
      COMPLEX*16 ck, Energy
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Output:
c     Sigma(ck,Energy)
      COMPLEX*16 dSigma
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Local Variables:
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Local Variables:
c     Wp     - Plasmon frequency, for broadened poles, Wp is not the
c              same as Wi.
c     Beta   - g_i for the negative weight pole at zero. 
c              Adding the negative weight pole at zero corrects the
c              diverging sum rule for epsilon^-1. This is irrelevant
c              for unbroadened poles.
c     NSing  - Number of singularities of integrand (used in cgratr)
c     NCalls - Number of fcn calls used in cgratr
c     MaxNR  - Max number of regions used in cgratr
c     DPPar  - Array of double precision parameters passed to cgratr
c              to be used in the functions r1, r2, and r3
c     CPar   - Array of complex parameters passed to cgratr
c              to be used in the functions r1, r2, and r3
c     Limit1 - Lower limit of integration
c     Limit2 - Upper limit of integration
c     HLInt1 - Integral of dr2 (derivative of first integral in eq. 13 of H.L.)
c     HLInt2 - Integral of dr1 (derivative second integral in eq. 13 of H.L.)
c     HLInt3 - Integral of dr1 or dr3 (derivative of 3rd or 4th integral)
c     XSing  - Array of singularities of integrand (used by cgratr)      
      INTEGER NSing, NCalls, MaxNR
      DOUBLE PRECISION DPPar(10), Beta
c      DOUBLE PRECISION Wp
      COMPLEX*16 CPar(10), Limit1, Limit2, HLInt1, HLInt2, HLInt3,
     &     XSing(20)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c     Loop Variables:
      INTEGER i1, i2
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c     Parameters:
c     ZeroPl - lower bound Limit1
c     Inf    - Upper bound of Limit2
c     AbsErr - absolute error used by cgratr
c     RelErr - Relative error used by cgratr
c     Error  - used for error codes by cgratr
      DOUBLE PRECISION ZeroPl, Inf, AbsErr, RelErr, Error
      PARAMETER(ZeroPl = 1.d-5, Inf = 1.d2)
      PARAMETER(AbsErr = 1.d-5, RelErr = 1.d-4)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Externals:
c     cgratr      - integration routine
c     dr1,dr2,dr3 - functions to integrate
c     HFExc       - Calculates Hartree Fock exchange
      COMPLEX*16 cgratr, dr1, dr2, dr3, HFExc
      EXTERNAL cgratr, dr1, dr2, dr3, HFExc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc            
c     Initialization:
      NSing  = 0
      NCalls = 0
      MaxNR  = 0
c     DPPar is array of dp parameters to evaluate functions in cgratr.
c     Everything is in dimensionless units.
c     1. xwg
      DPPar(1) = Wi/EFermi
c     2. xgam
      DPPar(2) = Gamma/EFermi
c     3. xe
      DPPar(3) = dble(Energy/EFermi)
c     4. xeg
      DPPar(4) = 0.d0
c     CPar is array of complex parameters to evaluate functions in cgratr.
c     ck in dimensionless units.
      CPar(1) = ck/kFermi
c     2. ce (complex energy)
      CPar(2) = Energy/EFermi + coni*DPPar(2)
      
c     Josh - This is a possible fix for functions that overlap zero by a large
c     amount so that Wp does not equal Wi. 
c     Wp= pi*Wi/2 + Wi*ArcTan[Wi/Gamma] - Gamma*Log[Beta] + Gamma*Log[Gamma] - 
c               1/2*Gamma*Log[Wi**2 + Gamma**2]
c      Beta = 1.d0*Gamma*SQRT(2.d0/(Wi**2+Gamma**2))
c      Beta  = 0.9d0*Gamma*Gamma/(Wi**2+Gamma**2)
c      Wp =  2*Gamma*LOG(Gamma/Beta) + 2*ATAN2(Wi,Gamma) -
c     &     Gamma*LOG(Gamma**2 + Wi**2)
c      Wp = SQRT(Wi*Wp)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      DO i1 = 1, 1
         IF(i1.eq.2) THEN
            DPPar(1) = 0.d0
            DPPar(2) = 0.9d0*Gamma/EFermi
         END IF         
         
c     Calculate derivatives of integrals in eq. 13 of H.L.      
c     1)
c     Integral d/dE{ ln[((ck+q)**2-E+Wq)/((ck-q)**2-E+Wq)] }
c     from ck + kFermi to Inf.      
         Limit1 = ck/kFermi+1.d0
         Limit2 = Inf
c     Find singularities in dr2
         iFcn = 2
         CALL FndSng(Limit1, Limit2, NSing, XSing, DPPar, CPar, iFcn)
         
c     Calculate integral
         HLInt1 = cgratr(dr2,DPPar,CPar,Limit1,Limit2,AbsErr,RelErr,
     &        NSing, XSing,Error,NCalls,MaxNR)
         DO i2 = 1, NSing
            XSing(i2) = (1.d0,0)*XSing(i2)
         END DO
         
c     2)
c     Integral d/dE{ ln[(kFermi**2-E-Wq)/(kFermi**2-E+Wq)*
c     ((ck+q)**2-E+Wq)/((ck-q)**2-E-Wq)] }
c     From ck - kFermi to ck + kFermi
         Limit1 = MAX(ABS(DBLE(ck)/kFermi-1.d0), ZeroPl)
         Limit2 = ck/kFermi+1.d0
         
c     Find singularities in dr1
         iFcn = 1
         CALL FndSng(Limit1, Limit2, NSing, XSing, DPPar, CPar, iFcn)
         HLInt2 = cgratr(dr1,DPPar,CPar,Limit1,Limit2,AbsErr,RelErr,
     &        NSing, XSing,Error,NCalls,MaxNR)
         DO i2 = 1, NSing
            XSing(i2) = (1.d0,0)*XSing(i2)
         END DO
         
c     3)
c     Theta(kFermi-Re(ck)) *
c     Integral { ln[((ck+q)**2-E-Wq)/((ck-q)**2-E-Wq)] +
c     Theta(Re(ck)-kFermi) *
c     Integral { ln[((ck+q)**2-E+Wq)/((ck-q)**2-E+Wq)]
c     (Integrals from 0 to kFermi - k and 0 to k - kFermi)
         Limit1 = ZeroPl
         Limit2 = ABS(DBLE(ck)/kfermi-1.d0)
         
c     If ck = kFermi, HLInt3 = 0      
         IF((ABS(DBLE(ck)-kFermi).lt.ZeroPl).or.
     &        (DBLE(Limit2).le.DBLE(Limit1))) THEN         
            HLInt3 = 0.d0
c     If ck < kFermi, HLint3 = Integral d/dE{ ln[((ck+q)**2-E-Wq)/((ck-q)**2-E-Wq)] }          
         ELSEIF(DBLE(ck).lt.kFermi) THEN
            Limit2 = 1.d0 - DBLE(ck)/kFermi
            
c     Find singularities in r3        
            iFcn = 3
            CALL FndSng(Limit1, Limit2, NSing, XSing, DPPar, CPar, iFcn)
            
c     Calculate integral         
            HLInt3 = cgratr(dr3,DPPar,CPar,Limit1,Limit2,AbsErr,RelErr,
     &           NSing,XSing,Error,NCalls,MaxNR)
            DO i2 = 1, NSing
               XSing(i2) = (1.d0,0)*XSing(i2)
            END DO
c     Else ck > kFermi, HLint3 = Integral d/dE{ ln[((ck+q)**2-E+Wq)/((ck-q)**2-E+Wq)] }         
         ELSE
            Limit2 = DBLE(ck)/kFermi - 1.d0
            
c     Find singularities in r2         
            iFcn = 2
            CALL FndSng(Limit1, Limit2, NSing, XSing, DPPar, CPar, iFcn)
            
c     Calculate integral         
            HLInt3 = cgratr(dr2,DPPar,CPar,Limit1,Limit2,AbsErr,RelErr,
     &           NSing,XSing,Error,NCalls,MaxNR)
            DO i2 = 1, NSing
               XSing(i2) = (1.d0,0)*XSing(i2)
            END DO
         END IF

         IF(i1.eq.1) THEN
            dSigma = - Amp*Wi*(Wi-coni*Gamma)/(2.d0*pi*EFermi*ck)*
     &           (HLInt1 + HLInt2 + HLInt3)
         ELSE
            dSigma = dSigma - Amp*Beta*
     &           Wi*(coni*Gamma)/(pi*EFermi*ck)*
     &           (HLInt1 + HLInt2 + HLInt3)
         END IF
      END DO
      
      RETURN
      END      

      FUNCTION HFExc(ckIn, EFermi, kFermi)
c     returns dirac-hara hartree-fock exchange
c     ck - complex momentum in units of kFermi
      include '../HEADERS/const.h'
      COMPLEX*16 ckIn, ck, HFExc, c
      DOUBLE PRECISION EFermi, kFermi
      ck = ckIn/kFermi
      c=-2.d0*EFermi/(pi*kFermi)
      IF(ABS(ck-1.d0).le.0.00001d0) THEN
         HFExc = c
      ELSE
         HFExc = c*(1.d0+(1.d0/ck-ck)* log( (1.d0+ck)/(ck-1.d0) )/2.d0)
      END IF
      RETURN
      END
      
c****************************************************************************
c     the following function routines are used for evaluating integals and
c     their derivatives.
c****************************************************************************
      complex*16 function r1(q,dppar,cpar)
      implicit double precision (a-h,o-z)
      include '../HEADERS/const.h'
      
c     Input:
      double precision dppar(4)
c     dppar contains:
c     xe, xeg, xwg, xgam
      complex*16 cpar(2)

c     Local Variables:
      complex*16 fq,fqq,fiq,a1,a2,a3,a4,t1,t2,q,ck,xe
      external fq

      ck=CPar(1)
c     3. xe
      xe = CPar(2)
c     4. xeg
      xeg = DPPar(4)
c     print*, 'call fq(q),ck', q, ck
      fqq=fq(q,dppar)
c     print*,'fqq=', fqq
      fiq=1./(q*fqq)
      a1=1.d0-xeg-xe-fqq - coni*1.d-10
      a2=(ck+q)**2-xe+fqq - coni*1.d-10
      a3=(ck-q)**2-xe-fqq - coni*1.d-10
      a4=1.d0+xeg-xe+fqq - coni*1.d-10
      t1=(a1*a2)
      t2=(a3*a4)
c     print*,'a1,a2,a3,a4,t1,t2',a1,a2,a3,a4,t1,t2
c      t1=t1/t2
      r1=fiq*(log(a1)+log(a2)-log(a3)-log(a4))
c     Test with r=E
c      r1=xe
c     print*,'r1 return to cgratr', r1
      return
      end
c****************************************************************************
      complex*16 function dr1(q,dppar,cpar)
      implicit double precision (a-h,o-z)
      include '../HEADERS/const.h'
      
c     Input:
      double precision dppar(4)
c     dppar contains:
c     xe, xeg, xwg, xgam
      complex*16 cpar(2)

c     Local Variables:
      complex*16 fq,fqq,fiq,a1,a2,a3,a4,q,ck,xe
      external fq

      ck=CPar(1)
c     3. xe
      xe = CPar(2)
c     4. xeg
      xeg = DPPar(4)
c     print*, 'call fq(q),ck', q, ck
      fqq=fq(q,dppar)
c     print*,'fqq=', fqq
      fiq=1./(q*fqq)
      a1=1.d0-xeg-xe-fqq - coni*1.d-10
      a2=(ck+q)**2-xe+fqq - coni*1.d-10
      a3=(ck-q)**2-xe-fqq - coni*1.d-10
      a4=1.d0+xeg-xe+fqq - coni*1.d-10

      dr1 = -fiq*(1.d0/a1+1.d0/a2-1.d0/a3-1.d0/a4)
c     Test with r=E
c      dr1=1.d0
c      write(51,*) dble(q), dble(dr1)
c     print*,'r1 return to cgratr', r1
      return
      end
c**********************************************************************
      complex*16 function r2(q,dppar,cpar)
      implicit double precision (a-h,o-z)
      include '../HEADERS/const.h'
      
c     Input:
      double precision dppar(4)
c     dppar contains:
c     xe, xeg, xwg, xgam
      complex*16 cpar(2)

c     Local Variables:
      complex*16 fqq,fq,fiq,a1,a2,q,ck,xe
      ck = CPar(1)
c     3. xe
      xe = CPar(2)
c     4. xeg
c      xeg = DPPar(4)
      
      fqq=fq(q,dppar)
      fiq=1.d0/(q*fqq)
      a1=((ck+q)**2-xe+fqq) - coni*1.d-10
      a2=((ck-q)**2-xe+fqq) - coni*1.d-10
      r2=fiq*(log(a1)-log(a2))
c     Test with r=E
c      r2=xe      
      return
      end
c**********************************************************************
      complex*16 function dr2(q,dppar,cpar)
      implicit double precision (a-h,o-z)
      include '../HEADERS/const.h'
      
c     Input:
      double precision dppar(4)
c     dppar contains:
c     xe, xeg, xwg, xgam
      complex*16 cpar(2)

c     Local Variables:
      complex*16 fqq,fq,fiq,a1,a2,q,ck,xe
      ck = CPar(1)
c     3. xe
      xe = CPar(2)
c     4. xeg
c      xeg = DPPar(4)
      
      fqq=fq(q,dppar)
      fiq=1.d0/(q*fqq)
      a1=((ck+q)**2-xe+fqq) - coni*1.d-10
      a2=((ck-q)**2-xe+fqq) - coni*1.d-10
      dr2=-fiq*(1.d0/a1-1.d0/a2)
c     Test with r=E
c      dr2=1.d0      
c      write(52,*) dble(q), dble(dr2)
      return
      end  
c**********************************************************************
      complex*16 function r3(q,dppar,cpar)
      implicit double precision (a-h,o-z)
      include '../HEADERS/const.h'

c     Input:
      double precision dppar(4)
c     dppar contains:
c     xe, xeg, xwg, xgam
      complex*16 cpar(2)

c     Local Variables:
      complex*16 fqq,fq,fiq,a1,a2,q,ck,xe
      ck = CPar(1)
c     3. xe
      xe = CPar(2)
c     4. xeg
c      xeg = DPPar(4)
c     valid only for k<kf, q<kf-k ?
      fqq=fq(q,dppar)
      fiq=1.d0/(q*fqq)
      a1=( (ck+q)**2-xe-fqq) - coni*1.d-10
      a2=( (ck-q)**2-xe-fqq) - coni*1.d-10
      r3=fiq*(log(a1) - log(a2))
c     Test with r=E
c      r3=xe      
      return
      end
c**********************************************************************
      complex*16 function dr3(q,dppar,cpar)
      implicit double precision (a-h,o-z)
      include '../HEADERS/const.h'

c     Input:
      double precision dppar(4)
c     dppar contains:
c     xe, xeg, xwg, xgam
      complex*16 cpar(2)

c     Local Variables:
      complex*16 fqq,fq,fiq,a1,a2,q,ck,xe
      ck = CPar(1)
c     3. xe
      xe = CPar(2)
c     4. xeg
c      xeg = DPPar(4)
c     valid only for k<kf, q<kf-k ?
      fqq=fq(q,dppar)
      fiq=1.d0/(q*fqq)
      a1=( (ck+q)**2-xe-fqq) - coni*1.d-10
      a2=( (ck-q)**2-xe-fqq) - coni*1.d-10
      dr3=-fiq*(1.d0/a1-1.d0/a2)
c     Test with r=E
c      dr3=1.d0
c      write(53,*) dble(q), dble(dr2)
      return
      end     
c**********************************************************************
      complex*16 function fq(q,dppar)
      implicit double precision (a-h,o-z)
      include '../HEADERS/const.h'
      complex*16 q
      double precision dppar(4)

c     1. xwg
      xwg = DPPar(1)
c     2. xgam
      xgam = DPPar(2)
      
c     Here I am going to change the dispersion relation to 
c     wq = wp + 1/2 * q**2 
c     This makes calculation of broadened poles easier. 
c      fq = xwg-coni*xgam + q**2
      
c     fq(q)=w1(q)=sqrt(w1**2+((omega(q)-omega_p)/omega_f)**2)
c     omega(q)**2=omega_p**2+omega_g**2(q)
c     fq(q)=xwg+a2*q**2+a4*q**4    xwg=(w1/ef)**2
c     electron gas parameters xwg=wp**2 a2=4/3, a4=1
c     uncomment the following 4 lines to use the old dispersion relation.
      a2=4.d0/3.d0
      a4=1.d0
      fq=(xwg-coni*xgam)**2 + a2*q*q + a4*q**4
      fq=sqrt(fq)
      return
      end

               
c*********************************************************************
c   This is Steve White's rewrite of Mike Teter's integration routine.  
c   Modified by J. Rehr for complex integration.
c   The following is a listing of the arguments in the initial function 
c   statement:
c     fn    -- routine requires external function statement in MAIN
c     xmin  -- lower limit
c     xmax  -- upper limit
c     abr   -- absolute tolerable error
c     rlr   -- relative tolerable error
c     nsing -- number of singularities or regions requiring 
c     special attention
c     xsing -- array of locations of singularities or endpoints
c     of special regions
c     error -- output for routine error messages
c     numcal-- the number of times fn was called
c     maxns -- the maximum number of regions being considered simultaneously
c     function cgratr(fn,xmin,xmax,abr,rlr,nsing,xsing,error,numcal,maxns)
c     fn declared double precision
c     double precision function cgratr(fn,xmin,xmax,abr,rlr,
c     fn declared complex*16
      
      complex*16 function cgratr(fn,dppar,cpar,xmin,xmax,abr,rlr,
     &     nsing,xsing,error,numcal,maxns)
      implicit double precision (a-h,o-z)
      parameter (mx=1500)
      integer nsing
      complex*16 fn,value,valu,fval(3,mx),xmax,xmin,del,del1
      complex*16 xleft(mx), xsing(20), cpar(10)
      double precision dppar(10)
      external fn
c     dimension xleft(mx),fval(3,mx),dx(3),wt(3)
      dimension wt9(9),dx(3),wt(3)
c     dimension xsing(20)
      logical atsing
      save dx,wt,wt9
      data dx/0.1127016653792583  ,0.5  ,0.8872983346207417  /
      data wt/0.277777777777777778  ,0.4444444444444444444  ,
     1     0.2777777777777777778  /
      data wt9/0.0616938806304841571  ,0.108384229110206161  ,
     1     0.0398463603260281088  ,0.175209035316976464  ,
     2     0.229732989232610220  ,0.175209035316976464  ,
     3     0.0398463603260281088  ,0.108384229110206161  ,
     4     0.0616938806304841571  /
c     nstack is the number of different intervals into which the 
c     integration region is currently divided. The number of regions can
c     grow if more accuracy is needed by dividing the right-most region
c     into three regions. The number of regions shrinks when the integral
c     over the right-most region is accurate enough, in which case that
c     integral is added to the total (stored in cgratr) and the region
c     is removed from consideration (and a new region is the right-most).
      nstack=nsing+1
      maxns = nstack
      error=0.  
      cgratr=0.  
c     The array xleft stores the boundary points of the regions.
c     The singular points just govern the initial placement of the regions.
      xleft(1)=xmin
      xleft(nsing+2)=xmax
      if(nsing.gt.0) then
         do 9 j=1,nsing
            xleft(j+1)=xsing(j)
 9       continue
      endif
c     For each region, calculate the function and store at three selected points.
      do 1 jj=1,nstack
         del=xleft(jj+1)-xleft(jj)
c     print*, 'fn call j= ,'
         do 1 j=1,3
c     print*, 'fn call in cgratr j= ',j
            fval(j,jj)=fn(xleft(jj)+del*dx(j),dppar,cpar)
 1    continue
c     print*, 'output of fn call, fval(j,jj)',fval(j,jj)
      numcal = nstack * 3
 6    continue
      if(nstack+3.ge.mx) then
         write(*,*) ' TOO MANY REGIONS'
         stop 0006
      endif
c     Divide the rightmost region into three subregions.  
      del=xleft(nstack+1)-xleft(nstack)
      xleft(nstack+3)=xleft(nstack+1)
      xleft(nstack+1)=xleft(nstack)+del*dx(1)*2.
      xleft(nstack+2)=xleft(nstack+3)-del*dx(1)*2.
c     The three data points already found for the region become the 
c     middle data points (number 2 in first index of fval) for each region.
      fval(2,nstack+2)=fval(3,nstack)
      fval(2,nstack+1)=fval(2,nstack)
      fval(2,nstack)=fval(1,nstack)
c     Now do the integral over the right-most region in two different ways-
c     a three point integral (valu) over each of the three subregions 
c     and a more accurate nine-point integral (value) over whole region.
c     valu is used only for the error estimate.
      icount=0
      value=0.  
      valu=0.  
      do 3 j=nstack,nstack+2
         del1=xleft(j+1)-xleft(j)
c     print*, 'fn call 2'
         fval(1,j)=fn(xleft(j)+dx(1)*del1,dppar,cpar)
         fval(3,j)=fn(xleft(j)+dx(3)*del1,dppar,cpar)
c     print*, 'fn call 2'
         numcal = numcal + 2
         do 5 k=1,3
            icount=icount+1
            value=value+wt9(icount)*fval(k,j)*del
            valu=valu+fval(k,j)*wt(k)*del1
 5       continue
 3    continue
      dif=abs(value-valu)
c     If the following condition is true, add in this integral to the total,
c     and reduce the number of regions under consideration.
      frac = dble(del / (xmax - xmin))
      atsing = .false.
      if(frac .le. 1.0e-8) atsing = .true.
      if(dif .le. abr*frac .or. dif.le.rlr*abs(value) .or. 
     1     (atsing .and. 
     2     (frac .le. 1.0e-15 .or. dif .le. abr*0.1  ))) then
c     The following commented out line is Teeter's old error criterion.
c     if(dif.le.abr.or.dif.le.rlr*abs(value))then
         cgratr=cgratr+value
         error=error+abs(dif)
         nstack=nstack-1
c     If no more regions, we are done.
         if(nstack.le.0) return
      else
c     If the integration is insufficiently accurate, make each of the 
c     three subregions of the right-most region into regions.
c     On next pass the right-most of these is the new current region.
         nstack=nstack+2
         maxns = max(maxns,nstack)
      endif
      go to 6
      end
