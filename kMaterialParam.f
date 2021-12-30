********************************************
** KMATERIAL sets the material constants  **
**      for elasticity and plasticity     **
********************************************
c
      SUBROUTINE kMaterialParam(iphase,caratio,compliance,G12,thermat,
     + gammast,burgerv,nSys,tauc,screwplanes,Temperature,
     + tauctwin,nTwin,twinon,nTwinStart,nTwinEnd,TwinIntegral,
     + singlecrystal, cubicslip)
c
      ! crystal type
      INTEGER,intent(in) :: iphase
c      
CHRISTOS START
      INTEGER,intent(in) :: singlecrystal
      INTEGER,intent(in) :: cubicslip
CHRISTOS END            
c
      ! number of slip systems
      INTEGER,intent(in) :: nSys
c
      ! current temperature
      REAL*8,intent(in) :: Temperature
c
      ! total number of twin systems
      INTEGER,intent(in) :: nTwin
c
      ! twin systems activation flag
      INTEGER,intent(in) :: twinon
c
      ! the active twins are the ones in the
      ! interval [nTwinStart,nTwinEnd] in the
      ! twin system file
      INTEGER,intent(in) :: nTwinStart,nTwinEnd
c
      ! c/a ratio for hcp crystals
      REAL*8,intent(out) :: caratio
c
      ! elastic compliance matrix in the crystal reference frame 
      REAL*8,intent(out) :: compliance(6,6)
c
      ! shear modulus for Taylor dislocation law
      REAL*8,intent(out) :: G12
c
      ! thermal eigenstrain to model thermal expansion
      REAL*8,intent(out) :: thermat(3,3)
c
      ! prefactor for SSD evolution equation
      REAL*8,intent(out) :: gammast 
c
      ! Burgers vectors
      REAL*8,intent(out) :: burgerv(nSys)
c
      ! critical resolved shear stress of slip systems
      REAL*8,intent(out) :: tauc(nSys)
c
      ! number of screw planes 
      integer,intent(out) :: screwplanes
c
      ! critical resolved shear stress of twin systems
      REAL*8,intent(out) :: tauctwin(nTwin)
c
******************************************
** The following parameters must be set **
c
      ! material name
      ! more materials can be added
      ! by modifying this subroutine
      ! materials available are the following:
      ! 'zirconium', 'alphauranium', 'tungsten', 'copper', 'carbide',
      ! 'olivine', 'CMSX4' (Nickel alloy)
      character(len=*), parameter :: matname = 'CMSX4' 
c
**       End of parameters to set       **
******************************************
c
      ! Burgers vector scalars
      REAL*8 :: burger1, burger2
c
      ! Elastic constants scalars
      REAL*8 :: E1, E2, E3, G13, v12, v13
      REAL*8 :: G23, v23
CHRISTOS START      
      REAL*8 :: C11, C12
CHRISTOS END      
c
      ! hcp: basal critical resolved shear stress
      REAL*8 :: XTAUC1
c
      ! hcp: prismatic critical resolved shear stress
      REAL*8 :: XTAUC2
c
      ! hcp: prismatic critical resolved shear stress
      REAL*8 :: XTAUC4
c
      ! thermal expansion coefficients
      REAL*8 :: alpha1, alpha2, alpha3
c
      INTEGER :: i, j
c
      ! select crystal type       
      SELECT CASE(iphase)
c
      case(0) !hcp
c
        SELECT CASE(matname)
c
        case('zirconium')
c
        ! Burgers vectors (um)
        burger1 = 2.28E-4   
        burger2 = 4.242E-4 ! sqrt(a^2 + c^2)                      
        caratio = 1.57
c
        ! elastic constants (MPa units)
        E1 = 289.38E3
        E3 = 335.17E3
        G12 = 132.80E3
        G13 = 162.50E3
        v12 = 0.09 
        v13 = 0.04
c
        ! CRSS (MPa units)
        XTAUC1 = 15.2 ! basal
        XTAUC2 = 67.7 ! prismatic
        XTAUC4 = 2000.0 ! pyramidal
c
        ! thermal expansion coefficients
        alpha1 = 9.5D-6
        alpha2 = alpha1
        alpha3 = 0.5895*alpha1
c
        ! prefactor for SSD evolution equation
       	gammast = 0.0
c
        case default
        WRITE(*,*)"Not sure what material"
        END SELECT
c
      ! elastic constants based on crystal symmetry
      E2 = E1
      G23 = G13
      v23 = v13
c
      ! assign Burgers vector scalars
      burgerv(1:6) = burger1
      burgerv(7:12) = burger2
c
      ! assign CRSS
      tauc(1:3) = XTAUC1
      tauc(4:6) = XTAUC2
      tauc(7:12) = XTAUC4
c
      case(1) !bcc
c
        SELECT CASE(matname)
c
        case('tungsten')
c
        ! CRSS (MPa units)
        XTAUC1 = 360.0
c
        ! Burgers vectors (um)
        burger1 = 2.74E-4
c
        ! elastic constants (MPa units)
        E1 = 421E3
        G12 = 164.4E3
        v12 = 0.28
c
        ! thermal expansion coefficients
        alpha1 = 9.5e-6
        alpha2 = alpha1
        alpha3 = 0.5895*alpha1
c
        ! prefactor for SSD evolution equation
       	gammast = 0.0
c
        case default
        WRITE(*,*)"Not sure what material"
        END SELECT
c
      ! elastic constants based on crystal symmetry
      E2 = E1
      E3 = E1
      v13 = v12
      v23 = v12
      G13 = G12
      G23 = G12
c    
      ! assign Burgers vector scalars
      burgerv = burger1
c
      ! assign CRSS: same for all slip systems
      tauc = XTAUC1
c     
      case(2) !fcc
c
        SELECT CASE(matname)
c
        case('copper')
c
        ! CRSS (MPa units)
        XTAUC1 = 20.0
c
        ! Burgers vectors (um)
        burger1 = 2.55E-4
c
        ! elastic constants (MPa units)
        E1 = 66.69E3
        v12 = 0.4189
        G12 = 75.4E3
c
        ! thermal expansion coefficient
        alpha1 = 13.0e-6
c
CHRISTOS START         
        case('CMSX4')
c
        ! CRSS (MPa units)
        if (Temperature .LE. 850.0) then   !  Celcius units
           tauc=-0.000000001051*Temperature**4 
     +                 + 0.000001644382*Temperature**3
     +                -0.000738679333*Temperature**2 + 0.128385617901*Temperature 
     +                +446.547978926622
           if (cubicslip == 1) then
             tauc(13:18)=-0.000000001077*Temperature**4 
     +                + 0.000001567820*Temperature**3
     +                -0.000686532147*Temperature**2 - 0.074981918833*Temperature 
     +                +571.706771689334
           end if
        else
           tauc=-1.1707*Temperature + 1478.9
           if (cubicslip == 1) then
             tauc(13:18)=-0.9097*Temperature + 1183
           end if
        end if
c
c
        ! temperature dependent stiffness constants (MPa units)
        if (Temperature .LE. 800.0) then   !  Celcius units
           C11=-40.841*Temperature +251300
           C12=-14.269*Temperature +160965
        else
       C11=0.111364*Temperature**2-295.136*Temperature+382827
       C12=-0.000375*Temperature**3+1.3375*Temperature**2
     +              -1537.5*Temperature+716000
        end if  
        G12=-0.00002066*Temperature**3+0.021718*Temperature**2
     +       -38.3179*Temperature+129864
c        C11=199E3 (temperature idnependent case - for verification)
c        C12=141E3
c        G12=93E3
        E1 = (C11-C12)*(C11+2*C12)/(C11+C12)
        v12 = E1*C12/((C11-C12)*(C11+2*C12))
c
        ! temperature dependent thermal expansion coefficient
        alpha1 = 9.119E-9*Temperature +1.0975E-5
c       alpha1 = 20.1E-6  (temperature idnependent case - for verification)
c
        ! prefactor for SSD evolution equation
        gammast = 0.0
CHRISTOS END        
c
        case default
        WRITE(*,*)"Not sure what material"
        END SELECT
c
      screwplanes = 2.0
c      
      ! elastic constants based on crystal symmetry
      E2 = E1
      E3 = E1
      v13 = v12
      v23 = v12
      G13 = G12
      G23 = G12
      ! isotropic thermal expansion                  
      alpha2 = alpha1
      alpha3 = alpha1
c      
      ! assign Burgers vector scalars
      burgerv = 2.55E-4
c
c
      case(3) ! carbide
c
        SELECT CASE(matname)
c
        case('carbide')
c
        ! CRSS (MPa units)
        XTAUC1 = 2300.0
c
        ! Burgers vectors (um)
        burger1 = 3.5072e-4
c
        ! elastic constants (MPa units)
        E1 = 207.0E+4
        v12 = 0.28
c
        ! thermal expansion coefficients
        alpha1 = 4.5e-6
        alpha2 = alpha1
        alpha3 = alpha1
c
        case default
        WRITE(*,*)"Not sure what material"
        END SELECT
c	  
      ! elastic constants based on crystal symmetry      
      E3 = E1
      v13 = v12
      G12 = E1/(2.0*(1.0+v12)) 
      G13 = G12
      G23 = G12
c    
      ! assign Burgers vector scalars
      burgerv = burger1
c
      ! assign CRSS: same for all slip systems
      tauc = XTAUC1
c      
      case(4) ! olivine
c
        SELECT CASE(matname)
c
        case('olivine')
c
        ! CRSS (MPa units)
        XTAUC1 = 2000.0
c
        ! Burgers vectors (um)
        burger1 = 2.74E-4
c
        ! elastic constants (MPa units)
        E1 =  421E3
        v12 = 0.28
        G12 = 164.4E3
c
        ! thermal expansion coefficients
        alpha1 = 4.5e-6
        alpha2 = alpha1
        alpha3 = alpha1
c
        case default
        WRITE(*,*)"Not sure what material"
        END SELECT
c
      ! elastic constants based on crystal symmetry       
      E2 = E1
      E3 = E1
      v13 = v12
      v23 = v12
      G13 = G12
      G23 = G12
c
      ! assign Burgers vector scalars
      burgerv = burger1
c
      ! assign CRSS
      tauc(1:7) = XTAUC1*(/2.0,1.0,1.0,1.5,4.0,4.0,1.0/)
c
      case(5) !alpha-Uranium
c
        SELECT CASE(matname)
c
        case('alphauranium')
c
        ! Burgers vectors (micrometre)
        burgerv(1:2) = 2.85e-4
        burgerv(3:4) = 6.51e-4
        burgerv(5:8) = 11.85e-4
c
        ! Constant factor tau_0^alpha for CRSS (MPa)
        ! Calhoun 2013 values
        ! with temperature dependence as in Zecevic 2016
        ! added after hardening is included
        tauc(1) = 24.5
        tauc(2) = 85.5
        tauc(3) = 166.5
        tauc(4) = 166.5
        tauc(5) = 235.0
        tauc(6) = 235.0
        tauc(7) = 235.0
        tauc(8) = 235.0
c
        ! Daniel 1971, Figure 9
        ! softest value (MPa)
        if (twinon == 1) then
          tauctwin(nTwinStart) = 6.25
          tauctwin(nTwinEnd) = 6.25
        end if
c
        ! elastic moduli (MPa) and Poissons ratios
        ! see PRS Literature Review by Philip Earp
        ! Short Crack Propagation in Uranium, an Anisotropic Polycrystalline Metal
        ! temperature dependence according to Daniel 1971
        E1 = 203665.987780 * (1.0 - 0.000935*(Temperature-293.0))
        E2 = 148588.410104 * (1.0 - 0.000935*(Temperature-293.0))
        E3 = 208768.267223 * (1.0 - 0.000935*(Temperature-293.0))
        v12 = 0.242363
        v13 = -0.016293
        v23 = 0.387816
        G12 = 74349.442379 * (1.0 - 0.000935*(Temperature-293.0))
        G13 = 73421.439060 * (1.0 - 0.000935*(Temperature-293.0))
        G23 = 124378.109453 * (1.0 - 0.000935*(Temperature-293.0))
c
        ! define thermal expansion coefficients as a function of temperature
        ! Lloyd, Barrett, 1966
        ! Thermal expansion of alpha Uranium
        ! Journal of Nuclear Materials 18 (1966) 55-59
        alpha1 = 24.22e-6 - 9.83e-9 * Temperature + 46.02e-12 * 
     +            Temperature * Temperature
        alpha2 = 3.07e-6 + 3.47e-9 * Temperature - 38.45e-12 * 
     +            Temperature * Temperature
        alpha3 = 8.72e-6 + 37.04e-9 * Temperature + 9.08e-12 * 
     +            Temperature * Temperature
c
        case default
        WRITE(*,*)"Not sure what material"
        END SELECT
c
      case default
      WRITE(*,*)"Not sure what crystal type"
      END SELECT
c
C     *** SET UP ELASTIC STIFFNESS MATRIX IN LATTICE SYSTEM ***  
C     notation as in: https://en.wikipedia.org/wiki/Hooke%27s_law 
C     However, the Voigt notation in Abaqus for strain and stress vectors
C     has the following order:
C     stress = (sigma11,sigma22,sigma33,tau12 tau13 tau23)
C     strain = (epsilon11,epsilon22,epsilon33,gamma12,gamma13,gamma23)
c
      compliance(1,1:3) = (/1./E1,-v12/E1,-v13/E1/)
      compliance(2,2:3) =         (/1./E2,-v23/E2/)
      compliance(3,3:3) =                 (/1./E3/)
      compliance(4,4:4) =                       (/1./G12/)
      compliance(5,5:5) =                       (/1./G13/)
      compliance(6,6:6) =                       (/1./G23/)
c
      ! symmetrize compliance matrix
      DO i=2,6
         DO j=1,i-1
            compliance(i,j)=compliance(j,i)
         END DO
      END DO
c
      ! define thermal eigenstrain in the lattice system
      thermat(1,1) = alpha1 
      thermat(2,2) = alpha2 
      thermat(3,3) = alpha3
c
      RETURN
c
      END

