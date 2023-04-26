***********************************************************
** Crystal Plasticity UMAT 
**  University of Oxford   
** Developed by Nicolo Grilli 2020,
** rewrite of UMAT by Ed Tarleton which was based on UEL by Fionn Dunne 2007
**********************************************************
c
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
c
      INCLUDE 'ABA_PARAM.INC'
c
#include <SMAAspUserSubroutines.hdr>
c
      CHARACTER*80 CMNAME
c
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
c
******************************************
** The following parameters must be set **
c      
      ! activate debug mode with Visual Studio
      ! 0 = off ; 1 = on
      integer, parameter :: debug = 0 
c
      ! activate GND calculations
      ! 0 = off ; 1 = on
      integer, parameter :: gndon = 0

      ! activate twin systems
      ! 0 = off ; 1 = on
      integer, parameter :: twinon = 0
c      
      ! total number of twins (active/inactive)
      integer, parameter :: TotalNTwin = 2
c
      ! the active twins are the ones in the
      ! interval [nTwinStart,nTwinEnd] in the
      ! twin system file
      integer, parameter :: nTwinStart = 1
      integer, parameter :: nTwinEnd = 2
c
      ! activate irradiation effect
      ! 0 = off ; 1 = on 
      integer, parameter :: irradiate = 0
      
c Activate cubic slip systems for single crystal FCC
      integer, parameter :: cubicslip = 1
      
c Activate dynamic incrementation for cyclic plasticity-creep analysis
      integer, parameter :: dynIncrementation = 1     
      
CCCCCC incrementation parameters for cyclic plasticity-creep analysis

      ! max allowable stress increment
      real*8, parameter :: maxStressInc = 50.0
      ! max time increment when max stress increment is exceeded
      real*8, parameter :: maxDtStress = 10.0
      real*8 :: ratioDtStress
      
      ! muliplier and exponent related to plastic strain rate
      
      real*8, parameter :: multiDt = 0.011
      real*8, parameter :: powerDt = -0.65
      real*8 :: ratioDtpStrain
CCCCCCC
**       End of parameters to set       **
******************************************
c
      ! dimension of the space
      PARAMETER (M=3,N=3)     

      ! Gauss points coordinates of the parent C3D8 element
      DIMENSION gauss(8,3), gausscoords(3,8)
c
      ! identity matrix
      DIMENSION xI(3,3)
c
      ! coordinates of the nodes of the parent element
      DIMENSION xnat(8,3)
c
      ! curl of the plastic deformation gradient
      DIMENSION curlfp(8,9)
c
      ! plastic deformation gradient
      DIMENSION Fp(8,9)
c
      ! stress components
      real*8 :: stressvec(6)
c
      ! stress increment
      real*8 :: dstressinc(6)
c
      ! strain increment
      real*8 :: dtotstran(6)
c
      ! total strain
      real*8 :: totstran(6)
c
      ! time derivative of the deformation gradient
      real*8 :: Fdot(3,3)
c
      ! inverse of the deformation gradient
      real*8 :: invF(3,3)
c
      ! deformation gradient at previous increment
      real*8 :: F(3,3)
c
      ! velocity gradient
      real*8 :: L(3,3)
c
      ! Von Mises invariant plastic strain rate
      real*8 :: pdot
c
      ! Von Mises stress
      real*8 :: vms
c
      integer :: i, j, K, debugWait
c
      ! number of slip systems
      integer :: nSys
c
      ! total number of twins
      integer :: nTwin
c
      ! number of screw systems
      integer :: ns
c
      ! number of active slip systems considered
      integer :: L0=12 ! HCP
      integer :: L1=12 ! BCC
      integer :: L2 ! FCC: variable because cubic systems can be included
      integer :: L4=7  ! Olivine
      integer :: LalphaUranium=8 ! alpha-uranium
c      
      ! crystal type
      ! first material constant in the input file
      ! 0 = HCP
      ! 1 = BCC
      ! 2 = FCC
      ! 3 = Carbide
      ! 4 = Olivine
      ! 5 = Orthorombic
      integer :: iphase
c
      ! integral of the twin volume fraction in a 1 um cylinder around
      ! the twin plane normal is calculated to determine the CRSS
      ! of the twin system in the discrete twin model
      real*8 :: TwinIntegral(TotalNTwin)
c
      ! position of the Gauss point in the C3D8 parent element
      PARAMETER ( xgauss = 0.577350269189626 ) 
c
      include 'mycommon.f'
c
      ! arrays to store element index and IP index
      ! of the neighbouring IPs
      integer, target :: Twin1UpDownEl(ArrayNUpDown) 
      integer, target :: Twin1UpDownIP(ArrayNUpDown) 
      !first twin sys storage
      integer, target :: Twin2UpDownEl(ArrayNUpDown) 
      integer, target :: Twin2UpDownIP(ArrayNUpDown) 
      ! second twin sys storage
c     
      pointer(ptrTwin1UpDownEl,Twin1UpDownEl)
      pointer(ptrTwin1UpDownIP,Twin1UpDownIP)
c
      pointer(ptrTwin2UpDownEl,Twin2UpDownEl)
      pointer(ptrTwin2UpDownIP,Twin2UpDownIP)
c      
      ! wait here to attach Visual Studio debugger
      do while (debug == 1)
          debugWait = 1
      end do
c
      ! if sinh( ) in the slip law has probably blown up 
      ! then try again with smaller dt
      if(any(DFGRD1 /= DFGRD1)) then
          call MutexLock( 1 )      ! lock Mutex #1
              pnewdt = 0.5
              write(*,*) "*** WARNING DFGRD1  = NaN: noel, npt, time: ",
     1        noel, npt, time
          call MutexUnlock( 1 )      ! lock Mutex #1
          return
      end if   
c
C increase slip systems if cubic slip is activated     
      if (cubicslip == 0) then
        L2=12
      else  
        L2=18
      end if
 
      ! read crystal type from input file
      ! first material constant
      iphase = int(props(1))      
c
      SELECT CASE(iphase)
      CASE(0) !hcp
      nSys = L0
      ns = 9 ! number of screw systems
c            
      case(1) !bcc
      nSys = L1
      ns = 4
c     
      case(2) !fcc
      nSys = L2
      ns = 6
c      
      case(3) !carbide
      nSys = L2
      ns = 0
c      
      case(4) !olivine
      nSys = L4
      ns = 2
c
      case(5) !alpha-Uranium
      nSys = 8
      ns = 0
c
      case default      
      WRITE(*,*)"Not sure what crystal type. Material constants."
      END SELECT
c	  
	  ! assign number of twins to initialize
	  ! arrays with the right size in kmat
	  nTwin = TotalNTwin
c
C     WRITE PROPS INTO SVARS TO INITIALIZE (ONCE ONLY),
c
      if (kinc <= 1 .and. kstep==1) then
c
         STATEV = 0.
c      
         ! read rotation matrix from the material constants
         ! 2 to 10 in the input file
         ! and assign to state variables 1 to 9
         ! order of the components in the input file must be
         ! R11, R12, R13, R21, R22, R23, R31, R32, R33
         do i=1,3
            do j=1,3
             STATEV(j+(i-1)*3) = props(j+1+((i-1)*3))
            end do
         end do
c         call lapinverse(T,M,info,Tinv)
c         do i=1,3
c           do j=1,3
c             STATEV(j+(i-1)*3) = Tinv(i,j)
c           end do
c         end do
c
         ! Initialize plastic deformation gradient
         ! to identity
	   STATEV(81) = 1.0
	   STATEV(85) = 1.0
         STATEV(89) = 1.0
c
         ! initialize cumulative plastic slip to zero
         STATEV(35) = 0.0
c      
         ! initialize hardening from defects
         if (irradiate == 1) then
           STATEV(36) = 750.0 
         else 
           STATEV(36) = 0.0 
         end if
c
         ! initialize sessile SSD density
         STATEV(54) = 0.01
c
      ! initialize dislocation density and twin phase field
      call kRhoTwinInit(nTwin,M,NSTATV,STATEV,twinon,
     +     nTwinStart,nTwinEnd,iphase,COORDS,nSys) 
c
         ! initialize plastic deformation gradient
         call MutexLock( 2 )      ! lock Mutex #2
             kFp(noel,npt,1:9) = 0.0
         call MutexUnlock( 2 )   ! unlock Mutex #2
c
         ! initialize gauss points coordinates 
         call MutexLock( 3 )      ! lock Mutex #3
             kgausscoords(noel,npt,1:3) = coords(1:3)
         call MutexUnlock( 3 )   ! unlock Mutex #3
c      
         ! initialize curl of plastic deformation gradient
         call MutexLock( 4 )      ! lock Mutex #4    
             kcurlFp(noel,npt,1:9) = 0.0
         call MutexUnlock( 4 )   ! unlock Mutex #4
c
         ! initialize number of neighbouring IPs
         call MutexLock( 7 )      ! lock Mutex #7
             kNoNeighbours(noel,npt,1:2) = 0
         call MutexUnlock( 7 )   ! unlock Mutex #7
c      
         ! initialize grain index
         ! it is the constant 11 in the
         ! material constants of the input file
         call MutexLock( 8 )      ! lock Mutex #8
             kGrainIndex(noel,npt) = props(11)
             NUpDownExceeded = 0
         call MutexUnlock( 8 )   ! unlock Mutex #8
c
         ! initialize maximum stress of the cohesive law
         ! and effective opening of fracture surface
         call MutexLock( 13 )      ! lock Mutex #13
             kSigma0(noel,npt) = 1300.0
             kDeltaEff(noel,npt) = 0.0
         call MutexUnlock( 13 )   ! unlock Mutex #13
c
         ! define local thread variables to store the element and IP
         ! indices that are within a length l0
         ! from the current noel and npt
         if (twinon == 1) then
           ptrTwin1UpDownEl = SMAIntArrayCreate(1,ArrayNUpDown,0)
           ptrTwin1UpDownIP = SMAIntArrayCreate(2,ArrayNUpDown,0)
           ptrTwin2UpDownEl = SMAIntArrayCreate(3,ArrayNUpDown,0)
           ptrTwin2UpDownIP = SMAIntArrayCreate(4,ArrayNUpDown,0)
         end if
c
      end if ! END OF WRITE PROPS INTO SVARS TO INITIALIZE (ONCE ONLY)
c
      ! access the arrays TwinUpDownEl, TwinUpDownIP
      ! and assign them to the pointers
      if (twinon == 1) then
        ptrTwin1UpDownEl = SMAIntArrayAccess(1)
        ptrTwin1UpDownIP = SMAIntArrayAccess(2)
        ptrTwin2UpDownEl = SMAIntArrayAccess(3)
        ptrTwin2UpDownIP = SMAIntArrayAccess(4)
      end if
c
! calculation of Twin1UpDownEl, Twin1UpDownIP
! calculation of Twin2UpDownEl, Twin2UpDownIP
! at the second increment when kgausscoords
! have been assigned for all elements and IPs
      if (kinc == 2 .and. kstep == 1 .and. twinon == 1) then
c
        call MutexLock( 9 )      ! lock Mutex #9
c
      call kfindneighbourhood(noel,npt,COORDS,STATEV,NSTATV,
     + Twin1UpDownEl,Twin1UpDownIP,Twin2UpDownEl,Twin2UpDownIP,
     + nTwinStart,nTwinEnd,nTwin)
c
        call MutexUnlock( 9 )   ! unlock Mutex #9
c
      end if 
      ! END OF CALCULATION OF TwinUpDownEl & TwinUpDownIP (ONCE ONLY)
c
      ! define identity matrix
      xI=0.            
      DO I=1,3
         xI(I,I)=1.
      END DO
c
C   DETERMINE DEFORMATION AND VELOCITY GRADIENTS
      stressvec = 0.0
      F=0.  
      invF = 0. 
      Fdot = 0.
      L=0.
      F = DFGRD0
      Fdot = (DFGRD1-DFGRD0)/DTIME
      CALL lapinverse(F,3,info,invF) 
      L = matmul(Fdot,invF)
c
      ! assign curl of plastic deformation gradient
      ! to state variables for output
c      DO i=1,9
c          statev(37+i) = kcurlfp(noel,npt,i)
c      END DO  
      
c
! calculate integral of the twin volume fraction
! for the two twin systems
! in the up/down region
      TwinIntegral(1:nTwin) = 0.0
c
      if (twinon == 1) then ! active twin
c
        if (kinc >= 2 .or. kstep > 1) then 
! TwinUpDownEl, TwinUpDownIP already assigned
c
      call ktwinvolfracintegral(noel,npt,TwinIntegral,
     + Twin1UpDownEl,Twin1UpDownIP,Twin2UpDownEl,Twin2UpDownIP,
     + nTwinStart,nTwinEnd,nTwin)
c
        end if ! check increment
c
      end if ! active twin
      
      
CCCccc   Dynamic incrementation for cyclic analysis of CMSX-4
      if (dynIncrementation == 1) then
        ! suggested incr ratio related to plastic strain (continuous)
        if (kinc >= 1) then
           ratioDtpStrain = (multiDt*STATEV(32)**powerDt)/dtime
        end if 
      
        ! upper bound of incr ratio used when stress increment is exceeded (non-continuous)
        if (STATEV(33) > maxStressInc) then 
           ratioDtStress=maxDtStress/dtime
        else
           ratioDtStress=10000.0   ! arbitrary high value
        end if
      
        ! use the minimum of the 2 criteria
        pnewdt=min(ratioDtpStrain,ratioDtStress)
      end if   
ccccccccccccccccccccccccccccccccccccccc      
c      
C   CALL KMAT FOR MATERIAL BEHAVIOUR - Global stiffness matrix C
C      and stress 
c
      call kmat(dtime,NSTATV,STATEV,xI,NOEL,NPT,time,F,
     +    L,iphase,irradiate,DDSDDE,stressvec,dstressinc,totstran,
     +    dtotstran,TEMP,DTEMP,vms,pdot,pnewdt,gndon,nSys,nTwin,ns,
     +    coords,TwinIntegral,nTwinStart,nTwinEnd,twinon,cubicslip)
c
      ! store twin phase field
      ! in common block variable
      call MutexLock( 14 )      ! lock Mutex #14
          kTwinVolFrac(noel,npt,1) = STATEV(106+nTwinStart)
          kTwinVolFrac(noel,npt,2) = STATEV(106+nTwinEnd)
      call MutexUnlock( 14 )   ! unlock Mutex #14
c
      ! store UEL variable in STATEV
      ! for output

c        STATEV(124) = kDeltaEff(noel,npt)
c        STATEV(125) = kSigma0(noel,npt)

C CHRISTOS - print deformation gradient
c      DO i=1,3
c       DO j=1,3 
c        statev(37+j+(i-1)*3) = F(i,j)
c       END DO
c      END DO 

c
C    RECOVER stress for calculating residual force
      DO K=1,6
           stress(K)=STATEV(47+K)
      END DO
c
      ! GND calculations
      if (gndon == 1) then
c
          call MutexLock( 5 )      ! lock Mutex #5            
              kFp(noel,npt,1:9)= statev(81:89)
          call MutexUnlock( 5 )    ! unlock Mutex #5 
c       
          IF (npt == 8 ) THEN ! update curl Fp
c        
C=======================================================================   
C   SPECIFY GAUSS POINT LOCAL COORDS, AND WEIGHTING CONSTANTS
c      
          xnat(1,:) = (/-1.,1.,1./)
          xnat(2,:) = (/-1.,-1.,1./)
          xnat(3,:) = (/-1.,1.,-1./)
          xnat(4,:) = (/-1.,-1.,-1./)
          xnat(5,:) = (/1.,1.,1./)
          xnat(6,:) = (/1.,-1.,1./)
          xnat(7,:) = (/1.,1.,-1./)
          xnat(8,:) = (/1.,-1.,-1./)      
c     
          do i = 1,8
            gauss(i,:)=(/xnat(i,1),xnat(i,2),xnat(i,3)/)*xgauss
          end do       
c         
          DO i=1,3         
            gausscoords(i,1:8) = kgausscoords(noel,1:8,i)
          END DO             
c
          Fp = kFp(noel,1:8,1:9)
C=======================================================================   
C   A FULL INTEGRATION GRADIENT SCHEME
c
          CALL kcurlET(curlFp,Fp,xnat,gauss,gausscoords) 
 ! faster and cleaned up version but same result as kcurl
c      
          call MutexLock( 6 )      ! lock Mutex #6
            kcurlFp(noel,1:8,1:9) = curlFp(1:8,1:9)     
          call MutexUnlock( 6 )      ! unlock Mutex #6  
c    
          END IF ! 8 IP check
      end if
c          
      RETURN
c
      END
c
      include 'kmat.f'
      include 'uexternaldb.f'
      include 'kdirns.f'
      include 'ktwinrot.f'
      include 'kslip0.f'
      include 'kslip5ET.f'
      include 'kslip6ET.f'
      include 'kslipPowerLaw.f'
      include 'kgndl2ET.f'
      include 'kcurlET.f'
      include 'kshapes.f'
      include 'utils.f'
      include 'UEL.for'
      include 'ktwinneighbourhood.f'
      include 'kRhoTwinInit.f'
      include 'kMaterialParam.f'
      include 'kCRSS.f'
      include 'kHardening.f'
CHRISTOS START
      include 'MPC.f'
      include 'DISPfull.f'
      include 'DLOADfull.f'
      include 'NickelSuperalloy.f'
      include 'kslipCreepPowerLaw.f'
      include 'dummyPowerLaw'
CHRISTOS END  

            
            
