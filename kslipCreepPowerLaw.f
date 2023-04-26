C   Christos Skamniotis & Nicolo Grilli
C   University of Oxford
C   December 2021 
C
C   Power law plasticity coupled with empirical creep law for tertiary creep:

      subroutine kslipCreepPowerLaw(xNorm,xDir,tau,signtau,tauc,
     +  dtime,nSys,iphase,currentTemperature,Lp,Lcr,tmat,gammaDot,gammaDotPlastic,
     +  gammaDotCreep,cubicslip,creep,softening,usvars,nsvars)
      
      INCLUDE 'ABA_PARAM.INC'
	  
	  ! number of slip system
      integer, intent(in):: nSys
        ! activation flag for cubic slip (additional 6 systems activated when loading is not along 001)
      INTEGER,intent(in) :: cubicslip
        ! activation flag for tertiary creep
      INTEGER,intent(in) :: creep
        ! activation flag for tertiary creep
      INTEGER,intent(in) :: softening
      
       ! number of Abaqus state variables
      INTEGER,intent(in) :: nsvars
      
      ! Abaqus state variables
      REAL*8,intent(in) :: usvars(nsvars)

	  ! phase
      integer, intent(in):: iphase

      ! slip directions and normals	  
      real*8, intent(in) :: xNorm(nSys,3),xDir(nSys,3)
	  
	  ! resolved shear stress and critical resolved shear stress
	  ! and sign of the resolved shear stress
	  ! tauc is positive by definition
      real*8, intent(in) :: tau(nSys), tauc(nSys), signtau(nSys)
	  
	  ! time step
      real*8, intent(in) :: dtime	  
	  
	  ! Temperature in C 
	  real*8, intent(in) :: currentTemperature
       
	  ! plastic velocity gradient
       real*8, intent(out) :: Lp(3,3)
       
        ! creep velocity gradient
        real*8, intent(out) :: Lcr(3,3)
	  
	  ! and its derivative with respect to the stress
	  real*8, intent(out) :: tmat(6,6)
	  
	  ! plastic strain rate on each slip system
	  real*8, intent(out) :: gammaDotPlastic(nSys)
        
        ! creep strain rate on each slip system
	  real*8, intent(out) :: gammaDotCreep(nSys)
        
        ! total inelastic strain rate on each slip system
	  real*8, intent(out) :: gammaDot(nSys)
        
        ! Gas constant (J*mol/K)
	  real*8, parameter :: R = 8.314462
      
        
******************************************
** The following parameters must be set **

*** RATE DEPENDENT PLASTICITY (thermally activated glide)

      ! reference strain rate (1/s)
	  real*8, parameter :: ref_gammaDot = 7.071136e-5
      ! rate sensitivity multiplier (1/Kelvin)
	  real*8, parameter :: slopeM = -0.036273
      ! rate sensitivity constant (-/-)
        real*8, parameter :: constantM = 65.33827694

*** TERTIARY CREEP (dislocation climb & damage)
!!!!  Initial creep rate constants        
      ! Activation energy for creep (J/mol)
      real*8, parameter :: Qo = 440.0e3
      ! reference rate (1/s)
      real*8, parameter :: ao = 2.0e7 
      ! stress multiplier (1/MPa)
      real*8, parameter :: bo = 5.0e-2
!!!!  Climb/damage constants  
      ! Activation energy for damage (J/mol)
      real*8, parameter :: QD = 170.0e3
      ! reference rate (1/s)
      real*8, parameter :: aD = 3.0e-1
      ! stress multiplier (1/MPa)
      real*8, parameter :: bD = 3.5e-2
!!!  decreasing strain softening at low Temp
      real*8, parameter :: TrefaD = 1203
      real*8, parameter :: multiaD = 15.0
!!!  increase stress sensitivity at high Temp
      real*8, parameter :: Trefbo = 1173.15
      real*8, parameter :: multibo = 5.0
!!!  creep strenghtening due to rafting at high Temp
      real*8, parameter :: Trefao = 1213.15
      real*8, parameter :: multiao = 22.0      
      
      ! prefactor parameters used to avoid divergence at low tau values
      ! multiplier
      real*8, parameter :: mult = 0.7
      ! power
      real*8, parameter :: pow = 20.0
	  
**       End of parameters to set       **
******************************************

        ! slip system index
        integer :: i
	  
	  ! Schmid tensor and its transpose
	  real*8 :: SNij(3,3), NSij(3,3)
	  
	  ! Schmid tensor and its transpose in Voigt notation
	  real*8 :: sni(6), nsi(6)
	  
	  ! higher order Schmid tensor in Voigt notation
	  real*8 :: SNNS(6,6)
	  
	  ! temporary slip normal and slip direction
	  real*8 :: tempNorm(3), tempDir(3)
	  
	  ! temporary variable to calculate the Jacobian
	  real*8 :: result1
	 
	  ! Jacobian
        real*8 :: result4(6,6)
        
        ! RSS/CRSS ratio
	  real*8 :: tau_ratio

	  ! rate sensitivity
	  real*8 :: mpower
        
        ! prefactor
	  real*8 :: prefactor
        
        ! creep constants varying with temperature
	  real*8 :: aDVar, aoVar, boVar
        
        real*8 :: Temperature
        
        Temperature = CurrentTemperature + 273.15
C
C  *** CALCULATE LP AND THE DERIVATIVE OF PLASTIC STRAIN INCREMENT WITH 
C   RESPECT TO THE STRESS DEFINED AS tmat***
C
C
        tmat = 0.0
	  Lp = 0.0
        Lcr = 0.0
	  result4 = 0.0
        mpower = slopeM * Temperature + constantM
        
c     compute creep constants for current temp
      aDVar=aD*exp(multiaD*(Temperature-TrefaD)/TrefaD)
      aoVar=ao*exp(multiao*(Trefao-Temperature)/Trefao)
      boVar=bo*exp(multibo*(Temperature-Trefbo)/Trefbo)
        
c      WRITE(*,*) "aDvar", aDVar
c      WRITE(*,*) "aovar", aoVar
c      WRITE(*,*) "bovar", boVar
      
      
	  ! contribution to Lp of all slip systems
      do i=1,nSys
	  
        tau_ratio=tau(i)/tauc(i)

       if  (tau_ratio > 0.0) then
            
           ! strain rate due to thermally activated glide (rate dependent plasticity)
           gammaDotPlastic(i) = signtau(i)*ref_gammaDot*tau_ratio**mpower
            
           ! calculate derivative d ( gammaDot(i) ) / d ( tau(i) )
           result1 = ref_gammaDot*mpower*(1/tauc(i))*(tau_ratio**(mpower-1))
          
        !  add tertiary creep rate
        if (creep == 1) then  
           
           prefactor=(mult*tau(i)/(1+mult*tau(i)))**pow
          
           gammaDotCreep(i) =  signtau(i)*prefactor*
     +                                    (  aoVar*exp(boVar*tau(i)-Qo/(R*Temperature))+
     +                    abs(usvars(89+i))*softening*aDVar*exp(bD*tau(i)-QD/(R*Temperature))  )
          
          result1 = result1 + prefactor*(aoVar*boVar*exp(boVar*tau(i)-Qo/(R*Temperature)) +
     +                  abs(usvars(89+i))*aDVar*softening*bD*exp(bD*tau(i)-QD/(R*Temperature)))+
     +              pow*(((mult*tau(i))/(1+mult*tau(i)))**(pow-1))*mult/((1+mult*tau(i))**2)* 
     +                     (  aoVar*exp(boVar*tau(i)-Qo/(R*Temperature))+
     +                    abs(usvars(89+i))*aDVar*softening*exp(bD*tau(i)-QD/(R*Temperature))  )
          
         else
            gammaDotCreep(i) = 0.0
          
         end if
        
        ! calculate SNNS
          tempNorm = xNorm(i,:)
          tempDir = xDir(i,:)
          SNij = spread(tempDir,2,3)*spread(tempNorm,1,3)
          NSij = spread(tempNorm,2,3)*spread(tempDir,1,3)
          call KGMATVEC6(SNij,sni)         
          call KGMATVEC6(NSij,nsi) 
          SNNS = spread(sni,2,6)*spread(nsi,1,6)
  
        ! contribution to Jacobian
          result4 = result4 + dtime*result1*SNNS     
		  
          gammaDot(i) = gammaDotPlastic(i) + gammaDotCreep(i) 
          
        ! plastic (inelastic) velocity gradient contribution
          Lp = Lp + gammaDot(i)*SNij
          
        ! creep velocity gradient contribution
          Lcr = Lcr + gammaDotCreep(i)*SNij
		
        else   
         
          gammaDotPlastic(i) = 0.0
          gammaDotCreep(i) = 0.0
          gammaDot(i) = 0.0
		  
        end if

      end do
      
      tmat = 0.5*(result4+transpose(result4))
      
      return
      end 