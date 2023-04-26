********************************************
** KHARDENING updates the state variables **
**  that determine mechanical hardening   **
********************************************
c
      SUBROUTINE kHardening(pdot,p,plasStrainrate,dtime,slip,gammaDot,
     + nSys,irradiate,tauSolute,gammast,rhossd,iphase,rhofor,
     + Temperature,rhosub,twinon,twinvolfrac,nTwin,
     + nTwinStart,nTwinEnd,gammaTwinDot)
c
      INCLUDE 'ABA_PARAM.INC'
c
      ! number of slip systems
      INTEGER,intent(in) :: nSys
c
      ! plastic strain rate
      REAL*8,intent(in) :: plasStrainRate(3,3)
c
      ! time increment
      REAL*8,intent(in) :: dtime
c
      ! plastic shear rate on slip systems
      ! and absolute value
      REAL*8,intent(in) :: gammadot(nSys)
c
      ! activate irradiation effect
      INTEGER,intent(in) :: irradiate
c
      ! prefactor for SSD evolution equation
      REAL*8,intent(in) :: gammast 
c
      ! crystal type
      INTEGER,intent(in) :: iphase
c
      ! Current temperature
      REAL*8,intent(in) :: Temperature
c
      ! twin systems activation flag
      INTEGER,intent(in) :: twinon
c
      ! total number of twin systems
      INTEGER,intent(in) :: nTwin
c
      ! the active twins are the ones in the
      ! interval [nTwinStart,nTwinEnd] in the
      ! twin system file
      INTEGER,intent(in) :: nTwinStart,nTwinEnd
c     
      ! twin rate
      REAL*8,intent(in) :: gammaTwinDot(nTwin)
c
      ! Von Mises invariant plastic strain rate
      REAL*8,intent(inout) :: p
c
      ! crystallographic slip
      ! needed for model with irradiation
      REAL*8,intent(inout) :: slip
c
      ! increase in tauc due to solute force
      REAL*8,intent(inout) :: tauSolute
c
      ! total SSD density
      REAL*8,intent(inout) :: rhossd
c
      ! forest dislocation density
      REAL*8,intent(inout) :: rhofor(nSys)
c
      ! substructure dislocation density
      REAL*8,intent(inout) :: rhosub 
c
      ! twin volume fraction
      REAL*8,intent(inout) :: twinvolfrac(nTwin)
c
      ! Von Mises invariant plastic strain rate
      REAL*8,intent(out) :: pdot
c
      ! absolute value of the plastic shear rate on slip systems
      REAL*8,dimension(nSys) :: tgammadot
c
      INTEGER :: I
c      
      ! calculate Von Mises invariant plastic strain rate
      pdot=sqrt(2./3.*sum(plasStrainrate*plasStrainrate)) 
c      
      p = p + pdot*dtime
c
      ! update crystallographic slip
      slip = slip + sum(abs(gammaDot))*dtime 
c       
      ! update solute force
      if (irradiate == 1) then           
        tauSolute = 750.0*exp(-slip/0.025)
      end if         
c
    !=========================================================================
    ! SSD Evolution
c     
      rhossd = rhossd + (gammast*pdot*dtime)
c      
    !=========================================================================
c
    !=========================================================================
    ! Dislocation density evolution (explicit Euler time integration)
    ! alpha-Uranium model
c
      if (iphase == 5) then
c
        ! slip independent of twin fraction
        DO I=1,nSys
          tgammaDot(I) = abs(gammaDot(I))
        END DO
c
        ! forest dislocations evolution
        ! using constants from calibration of tensile bar 3
        ! using twin-slip interaction model
        rhofor(1) = rhofor(1) + 43.2*max(sqrt(rhofor(1))-
     +(0.17100+2.6093e-03*Temperature)*rhofor(1),0.0)*tgammaDot(1)*dtime
        rhofor(2) = rhofor(2) + 6320.0*max(sqrt(rhofor(2))-
     +(0.25650+5.8708e-04*Temperature)*rhofor(2),0.0)*tgammaDot(2)*dtime
        rhofor(3) = rhofor(3) + 0.24*max(sqrt(rhofor(3))-
     +(0.11718+1.9289e-04*Temperature)*rhofor(3),0.0)*tgammaDot(3)*dtime
        rhofor(4) = rhofor(4) + 0.24*max(sqrt(rhofor(4))-
     +(0.11718+1.9289e-04*Temperature)*rhofor(4),0.0)*tgammaDot(4)*dtime
        rhofor(5) = rhofor(5) + 800.0*max(sqrt(rhofor(5))-
     +(0.12+1.5e-05*Temperature)*rhofor(5),0.0)*tgammaDot(5)*dtime
        rhofor(6) = rhofor(6) + 800.0*max(sqrt(rhofor(6))-
     +(0.12+1.5e-05*Temperature)*rhofor(6),0.0)*tgammaDot(6)*dtime
        rhofor(7) = rhofor(7) + 800.0*max(sqrt(rhofor(7))-
     +(0.12+1.5e-05*Temperature)*rhofor(7),0.0)*tgammaDot(7)*dtime
        rhofor(8) = rhofor(8) + 800.0*max(sqrt(rhofor(8))-
     +(0.12+1.5e-05*Temperature)*rhofor(8),0.0)*tgammaDot(8)*dtime
c
        ! substructure dislocations evolution
        rhosub = rhosub + 0.216*(17.545+0.26771*Temperature)*
     +            rhofor(1)*sqrt(rhosub)*tgammaDot(1)*dtime
c
        ! twin volume fraction evolution
        ! dot(f)^beta = dot(gamma)^beta / gamma^twin
        ! see Kalidindi JMPS 1998 (equations 32 and 33a)
        ! if twinvolfractot = 1.0 => all gammaTwinDot are zero => no evolution
        ! McCabe 2010: (130) twin induces a total strain 0.299
        if (twinon == 1) then ! twin active
          twinvolfrac(nTwinStart) = twinvolfrac(nTwinStart) + 
     +              gammaTwinDot(nTwinStart)*dtime/0.299
          twinvolfrac(nTwinEnd) = twinvolfrac(nTwinEnd) + 
     +              gammaTwinDot(nTwinEnd)*dtime/0.299
        end if
c
      end if
c
      RETURN
c
      END
