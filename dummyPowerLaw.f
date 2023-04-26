C   Christos Skamniotis & Nicolo Grilli
C   University of Oxford
C   December 2021 
C
C   Power law plasticity coupled with empirical creep law for tertiary creep:

      subroutine dummyPowerLaw(xNorm,xDir,tau,signtau,tauc,
     +  dtime,nSys,Lp,
     +  tmat,gammaDot)
      
      implicit none
	  
	  ! number of slip system
      integer, intent(in):: nSys

      ! slip directions and normals	  
      real*8, intent(in) :: xNorm(nSys,3),xDir(nSys,3)
	  
	  ! resolved shear stress and critical resolved shear stress
	  ! and sign of the resolved shear stress
	  ! tauc is positive by definition
       real*8, intent(in) :: tau(nSys), tauc(nSys), signtau(nSys)
	  
	  ! time step
       real*8, intent(in) :: dtime	  
	  
	  ! plastic velocity gradient
       real*8, intent(out) :: Lp(3,3)
	  
	  ! and its derivative with respect to the stress
	  real*8, intent(out) :: tmat(6,6)
	  
	  ! plastic strain rate on each slip system
	  real*8, intent(out) :: gammaDot(nSys)
        
******************************************
        ! reference strain rate (1/s)
	  real*8, parameter :: ref_gammaDot = 1e-2
	  
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

C
C  *** CALCULATE LP AND THE DERIVATIVE OF PLASTIC STRAIN INCREMENT WITH 
C   RESPECT TO THE STRESS DEFINED AS tmat***
C
        tmat = 0.0
	  Lp = 0.0
	  result4 = 0.0
        mpower = 20.0
        
	  ! contribution to Lp of all slip systems
      do i=1,nSys
	  
       tau_ratio=tau(i)/tauc(i)

       if  (tau_ratio > 0.0) then
            
            ! strain rate due to thermally activated glide (rate dependent plasticity)
            gammaDot(i) = signtau(i)*ref_gammaDot*tau_ratio**mpower
            
            ! calculate derivative d ( gammaDot(i) ) / d ( tau(i) )
            result1 = ref_gammaDot*mpower*(1/tauc(i))*(tau_ratio**(mpower-1))
        
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
          
        ! plastic velocity gradient contribution
          Lp = Lp + gammaDot(i)*SNij
        !  WRITE(*,*) "Lp", Lp
		
        else   
         
          gammaDot(i) = 0.0
		  
        end if

      end do
      
      !WRITE(*,*) "Lp", Lp
      
      tmat = 0.5*(result4+transpose(result4))
             
      return
      end 