********************************************************
** KCRSS calculates the CRSS of slip and twin systems **
********************************************************
c
      SUBROUTINE kCRSS(iphase,tauc,nSys,G12,burgerv,gndtot,irradiate,
     + tauSolute,gndcut,rhofor,rhosub,Temperature,homogtwin,
     + nTwinStart,nTwinEnd,twinvolfrac,tauctwin,nTwin,TwinIntegral,
     + twinvolfractotal,twinon)
c
      INCLUDE 'ABA_PARAM.INC'
c
      ! crystal type
      INTEGER,intent(in) :: iphase               
c
      ! number of slip systems
      INTEGER,intent(in) :: nSys
c
      ! total number of twin systems
      INTEGER,intent(in) :: nTwin
c
      ! shear modulus for Taylor dislocation law
      REAL*8,intent(in) :: G12
c
      ! Burgers vectors
      REAL*8,intent(in) :: burgerv(nSys)
c
      ! scalar total GND density
      REAL*8,intent(in) :: gndtot
c
      ! activate irradiation effect
      INTEGER,intent(in) :: irradiate
c
      ! increase in tauc due to solute force
      REAL*8,intent(in) :: tauSolute
c
      ! GND density (immobile)
      REAL*8,intent(in) :: gndcut(nSys)
c
      ! forest dislocation density
      REAL*8,intent(in) :: rhofor(nSys)
c
      ! substructure dislocation density
      REAL*8,intent(in) :: rhosub 
c
      ! Current temperature
      REAL*8,intent(in) :: Temperature
c
      ! homogenize twin model
      INTEGER,intent(in) :: homogtwin
c
      ! the active twins are the ones in the
      ! interval [nTwinStart,nTwinEnd] in the
      ! twin system file
      INTEGER,intent(in) :: nTwinStart,nTwinEnd
c	  
	  ! twin systems activation flag
      INTEGER,intent(in) :: twinon
c
      ! twin volume fraction
      REAL*8,intent(in) :: twinvolfrac(nTwin)
c
      ! average of the twin volume fraction
      ! over the neighbourhood
      ! two twin systems
      REAL*8,intent(in) :: TwinIntegral(nTwin)
c
      ! total twin volume fraction 
      REAL*8,intent(in) :: twinvolfractotal
c
      ! critical resolved shear stress of slip systems
      REAL*8,intent(inout) :: tauc(nSys)
c
      ! critical resolved shear stress of twin systems
      REAL*8,intent(inout) :: tauctwin(nTwin)
c
      INTEGER :: i
c
      ! check crystal type
      if (iphase == 1) then
c
          ! Taylor dislocation law
          tauc = tauc + 0.0065*G12*(burgerv(1))*sqrt(gndtot)
c         
          if (irradiate == 1) then
              tauc = tauc + tauSolute
          end if
c
      else if (iphase == 2) then
c
          ! Taylor dislocation law
c          tauc = tauc + 0.32*G12*(burgerv(1))*sqrt(gndcut)
           tauc = tauc   ! simple case of perfect plasticity (Christos)
c
      else if (iphase == 5) then 
c
          ! alpha-Uranium model with forest and substructure dislocations
          ! R.J. McCabe, L. Capolungo, P.E. Marshall, C.M. Cady, C.N. Tomé
          ! Deformation of wrought uranium: Experiments and modeling
          ! Acta Materialia 58 (2010) 5447–5459
          tauc(1) = tauc(1) + 19.066 * sqrt(rhofor(1)) + 1.8218 * 
     +         sqrt(rhosub) * log(1.0 / (burgerv(1) * sqrt(rhosub)))
          tauc(2) = tauc(2) + 18.832 * sqrt(rhofor(2)) + 1.7995 * 
     +           sqrt(rhosub) * log(1.0 / (burgerv(2) * sqrt(rhosub)))
          tauc(3) = tauc(3) + 54.052 * sqrt(rhofor(3)) + 5.1650 * 
     +           sqrt(rhosub) * log(1.0 / (burgerv(3) * sqrt(rhosub)))
          tauc(4) = tauc(4) + 54.052 * sqrt(rhofor(4)) + 5.1650 * 
     +           sqrt(rhosub) * log(1.0 / (burgerv(4) * sqrt(rhosub)))
          tauc(5) = tauc(5) + 123.357 * sqrt(rhofor(5)) + 11.7875 * 
     +           sqrt(rhosub) * log(1.0 / (burgerv(5) * sqrt(rhosub)))
          tauc(6) = tauc(6) + 123.357 * sqrt(rhofor(6)) + 11.7875 * 
     +           sqrt(rhosub) * log(1.0 / (burgerv(6) * sqrt(rhosub)))
          tauc(7) = tauc(7) + 123.357 * sqrt(rhofor(7)) + 11.7875 * 
     +           sqrt(rhosub) * log(1.0 / (burgerv(7) * sqrt(rhosub)))
          tauc(8) = tauc(8) + 123.357 * sqrt(rhofor(8)) + 11.7875 * 
     +           sqrt(rhosub) * log(1.0 / (burgerv(8) * sqrt(rhosub)))
c
          ! Zecevic 2016 temperature dependence
          tauc(1) = tauc(1) * exp(-(Temperature-293.0)/140.0)
          tauc(2) = tauc(2) * exp(-(Temperature-293.0)/140.0)
          tauc(3) = tauc(3) * exp(-(Temperature-293.0)/140.0)
          tauc(4) = tauc(4) * exp(-(Temperature-293.0)/140.0)
          tauc(5) = tauc(5) * exp(-(Temperature-293.0)/140.0)
          tauc(6) = tauc(6) * exp(-(Temperature-293.0)/140.0)
          tauc(7) = tauc(7) * exp(-(Temperature-293.0)/140.0)
          tauc(8) = tauc(8) * exp(-(Temperature-293.0)/140.0)
c
          ! Daniel, Lesage, 1971 minimum value as minima
          tauc(1) = max(tauc(1),4.0)
          tauc(2) = max(tauc(2),4.0)
          tauc(3) = max(tauc(3),4.0)
          tauc(4) = max(tauc(4),4.0)
          tauc(5) = max(tauc(5),4.0)
          tauc(6) = max(tauc(6),4.0)
          tauc(7) = max(tauc(7),4.0)
          tauc(8) = max(tauc(8),4.0)
c
          if (twinon == 1) then ! twin active
c
            if (homogtwin == 1) then ! homogenized twin model
c
	      ! add cross hardening of one twin system on the other
              DO i=nTwinStart,nTwinEnd
                tauctwin(i) = tauctwin(i) + 96.79*twinvolfractotal
              END DO
c
            else ! discrete twin model
c
              ! add hardening in the nucleation stage
              ! 50% is the critical twin volume fraction at which the
              ! softest value is reached
              DO i=nTwinStart,nTwinEnd
                if (twinvolfrac(i) < 0.5) then
                  tauctwin(i) = tauctwin(i) + 37.5*(0.5-twinvolfrac(i))
                  tauctwin(i) = tauctwin(i) + 2000.0*TwinIntegral(i)
                else ! twinvolfrac(i) > 0.5
                  tauctwin(i) = tauctwin(i) + 37.5*(twinvolfrac(i)-0.5)
                  tauctwin(i) = tauctwin(i) + 2000.0*TwinIntegral(i)
                end if
              END DO
c
              ! local interaction between twin systems
              ! only when threshold is passed
              if (twinvolfrac(nTwinEnd) > 0.5) then
                tauctwin(nTwinStart) = tauctwin(nTwinStart) + 
     +       200.0*twinvolfrac(nTwinEnd)
              end if
              if (twinvolfrac(nTwinStart) > 0.5) then
                tauctwin(nTwinEnd) = tauctwin(nTwinEnd) + 
     +       200.0*twinvolfrac(nTwinStart)
              end if
c
            end if ! choice of twin model (homogeneous/discrete)
c
          end if ! twin active
c
      end if ! check crystal type
c
      RETURN
c
      END
