C   Christos Skamniotis
C   University of Oxford
C   March 2021 
c
      subroutine FatigueLife(tau,signtau,time,dtime,TEMP,DTEMP,
     +          nSys,gammaDot,gammaDotPlastic,gammaDotCreep,
     +          cubicslip,creep,usvars,nsvars)
      
      INCLUDE 'ABA_PARAM.INC'
	  
	  ! number of slip system
      integer, intent(in):: nSys
        ! activation flag for cubic slip (additional 6 systems activated when loading is not along 001)
      INTEGER,intent(in) :: cubicslip
        ! activation flag for tertiary creep
      INTEGER,intent(in) :: creep
        ! number of Abaqus state variables
      INTEGER,intent(in) :: nsvars
        ! Abaqus state variables
      REAL*8,intent(inout) :: usvars(nsvars)
        ! resolved shear stress
      real*8, intent(in) :: tau(nSys)
        ! sign of resolved shear stress
      real*8, intent(in) :: signtau(nSys)
        ! time: time(1) = step time; time(2) = total time
      REAL*8,intent(in) :: time(2)
	  ! time step
      real*8, intent(in) :: dtime
        ! temperature step
      real*8, intent(in) :: DTEMP	
	  ! Temperature in Kelvin 
      real*8, intent(in) :: TEMP
	  ! plastic gamma rate on each slip system
      real*8, intent(in) :: gammaDotPlastic(nSys)
        ! creep gamma rate on each slip system
      real*8, intent(in) :: gammaDotCreep(nSys)
        ! total inelastic gamma rate on each slip system
      real*8, intent(in) :: gammaDot(nSys)
        ! Gas constant (J*mol/K)
      real*8, parameter :: R = 8.314462
        ! macro stress and strain used to calculate princ components
      real*8 :: stresses(6), strains(6), Creepstrains(6)
        ! princ components
      real*8 :: princ(3), princDir(3,3)

      
c  PARAMETERS
******************************************
c  LOADING CYCLE
      real*8, parameter ::  period=5000
      
      
      
      
c CREEP FAILURE STRAIN
      real*8, parameter :: FailStrain = 0.1

c FATIGUE LIFE MODEL        
c constants obtained by converting de-Nf data to gamma-Nf data from:
c A. Scholz, Y.Wang, S. Linn, C. Berger, R. Znajda. Modeling of mechanical properties of alloy CMSX-4
      real*8, parameter :: mult1oct = 132399.61
      real*8, parameter :: mult2oct = -3873.78
      real*8, parameter :: mult1cube = 46300.47
      real*8, parameter :: mult2cube = -1323.28
c Fatigue-creep Interaction Coefficient (calibrated)       
      real*8, parameter :: Interaction = 0
        
******************************************
      ! slip system index
      integer :: i, iCrit, id
      ! number of cycle simulated
      integer :: k
      ! other
      real*8 :: gammaCyclic, creepInc, maxVal, maxCreepstrain
        
      
******************************************
c calculation for loading cycle
      k=INT(time(2)/period)  ! k=0 for first cycle
      
      
      
******************************************
C   CALCULATIONS OF PRINC COMPONENTS, DIRECTION of Smax & max princ creep strain    

c   determine princ stress
      stresses=usvars(48:53)
      CALL SPRIND(stresses, princ, princDir, 1, 3, 3)  ! princ values & directions
      usvars(146:148)=princ
      maxVal=1e-20   ! initial small value
      do i=1,3
        if (abs(princ(i)) > maxVal) then
           maxVal = abs(princ(i))
           id=i
        endif
      end do    
      usvars(149)=princ(id)   ! current critical princ stress                  
      usvars(150:152)=maxVal*princDir(id,:)   ! components of critical princ stress
c      usvars(150:152)=abs(princDir(id,:))     ! magnitudes of direction cosines
          
c  determine total princ strains      
      strains=usvars(11:16)
      CALL SPRINC(strains, princ, 2, 3, 3)   ! only princ values
      usvars(153:155)=princ
      maxVal=1e-20   ! initial small value
      do i=1,3
        if (abs(princ(i)) > maxVal) then
           maxVal = abs(princ(i))
           id=i
        endif
      end do    
      usvars(156)=princ(id)   ! current critical princ strain
      
c  determine princ creep strains      
      Creepstrains=usvars(57:62)
      CALL SPRINC(Creepstrains, princ, 2, 3, 3)   ! only princ values
      maxVal=1e-20   ! initial small value
      do i=1,3
        if (abs(princ(i)) > maxVal) then
           maxVal = abs(princ(i))
           id=i
        endif
      end do    
      maxCreepstrain=princ(id)   ! current critical princ creep strain
      
*****************************************      
      
      
      
      
*************************************     
CCCCCCC  caclulations for each slip system
        maxVal=1e-20   ! initial small value
        creepInc=0.0
      
        do i=1,nSys
            
          ! signed gamma in 18 slip systems (90-107)          
          usvars(89+i) = usvars(89+i) + gammaDot(i)*dtime
          
          ! accumulated gamma in 18 slip systems (108-125)
          usvars(107+i) = usvars(107+i) + abs(gammaDot(i))*dtime
          
          ! find critical system
          if (usvars(107+i) > maxVal) then
             
             maxVal= usvars(107+i)
             iCrit=i
             
          end if   
      
          ! additive accumulated creep gamma increment (accounting all slip systems)
          creepInc=creepInc + abs(gammaDotCreep(i))*dtime
          
        end do
          
        
        
c   identify beggining of cycle when k changes value
         if (k > usvars(126)) then  
             
           usvars(137) = usvars(127)  !  accumulated plastic gamma per cycle for critical system 
           usvars(138) = usvars(128)  !  accumulated creep gamma per cycle for critical system
           usvars(139) = usvars(129)  !  accumulated creep gamma per cycle accounting all systems
           usvars(140) = usvars(130)  !  accumulated max princ macro creep strain per cycle
c      CREEP DAMAGE per cycle
           usvars(141) = usvars(139)*0.4082/FailStrain  !  creep damage per cycle based on all slip systems
           usvars(142) = usvars(140)/FailStrain  !  creep damage per cycle based on max macro creep strain
           
           ! Assign zero values to accumulated plastic and creep strains & damage per cycle for calculation in the next cycle
           usvars(127:130) = 0
      
         end if
         
         
C updating variables (NOT OUTPUT)         
c  cycle number identifier      
         usvars(126)=k
! update accumulated plastic gamma per cycle (for critical slip system)
          usvars(127) = usvars(127) + abs(gammaDotPlastic(iCrit))*dtime
         
! update accumulated creep gamma per cycle (for critical slip system)
          usvars(128) = usvars(128) + abs(gammaDotCreep(iCrit))*dtime
          
! additive accumulated creep gamma (accounting all slip systems)
          usvars(129) = usvars(129) +creepInc
          
! additive accumulated max creep macro strain
          usvars(130) = usvars(130) +abs(maxCreepstrain  -usvars(131))    ! note 131 is the old value
          
! update old creep macro strain          
         usvars(131)=maxCreepstrain
          
         
          
C OUTPUT variables         
          ! critical system
          usvars(132) = iCrit
          ! signed tau value for critical system
          usvars(133) = signtau(iCrit)*tau(iCrit)
          ! gamma for critical system
          usvars(134) = usvars(89+iCrit)
          if (usvars(134) == 1) then
             usvars(134) = 0.0
          endif   
          ! accumulated gamma for critical system -- USED FOR CHOOSING FATIGUE LOCATION
          usvars(135) = maxVal
          ! additive accumulated gamma (accounting all slip systems) -- USED FOR CHOOSING CREEP LOCATION
          usvars(136) = sum(usvars(108:125))       

         
          
c  FATIGUE DAMAGE per cycle
         gammaCyclic=usvars(137)+usvars(138)
         if (iCrit <= 12) then  
           usvars(143) = 1/ ( mult1oct * exp(gammaCyclic*mult2oct))   ! octahedral slip
         else
           usvars(143) = 1/ (mult1cube * exp(gammaCyclic*mult2cube))    ! cubic slip
         endif
         
       
C  TOTAL DAMAGE per cycle
! based on gamma creep
      usvars(144) = usvars(143) + usvars(141) + Interaction*usvars(143)*usvars(141)/(usvars(143)+usvars(141))
! based on macro creep strain
      usvars(145) = usvars(143) + usvars(142) + Interaction*usvars(143)*usvars(142)/(usvars(143)+usvars(142))
      
      
      
      return
      end
     
      
     
     
      