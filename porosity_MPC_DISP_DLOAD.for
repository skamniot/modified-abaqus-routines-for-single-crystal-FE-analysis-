      SUBROUTINE MPC(UE,A,JDOF,MDOF,N,JTYPE,X,U,UINIT,MAXDOF,LMPC,
     1               KSTEP,KINC,TIME,NT,NF,TEMP,FIELD,LTRAN,TRAN)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION UE(MDOF), A(MDOF,MDOF,N), JDOF(MDOF,N), X(6,N),
     1          U(MAXDOF,N), UINIT(MAXDOF,N), TIME(2), TEMP(NT,N),
     2          FIELD(NF,NT,N)
c      
      if (JTYPE.eq.1) then  
         JDOF(1,1)=3
         JDOF(1,2)=3     
         A(1,1,1) = 1
         A(1,1,2) =-1
      else if (JTYPE.eq.2) then
         JDOF(1,1)=3
         JDOF(1,2)=3     
         A(1,1,1) = 1
         A(1,1,2) = -1
      else if (JTYPE.eq.3) then
         JDOF(1,1)=1
         JDOF(1,2)=1     
         A(1,1,1) = 1
         A(1,1,2) =-1        
      else              
        JDOF(1,1)=1
        JDOF(2,1)=2
        JDOF(1,2)=1
        JDOF(2,2)=2        
        JDOF(1,3)=1
        A(1,1,1) = 1
        A(1,1,2) =-1
        A(1,1,3) =-1
        A(2,2,1) = 1
        A(2,2,2) =-1            
      end if
C  Give the value of dependent degree of freedom by solving the constraint equation about u1
c     print *, 'UE_default', UE(1)
      RETURN
      END
c
c


      SUBROUTINE  DISP(U,KSTEP,KINC,TIME,NODE,NOEL,JDOF,COORDS)
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION U(3),TIME(2),COORDS(3)
c  
      real*8, parameter ::  t_ramp=1    ! cases 1,2,3,4   (ramp duration)
c      
c        
c     CHOOSE PARAMETERS FOR DOUBLE PLATE
      real*8, parameter ::  Tmax=1050, DT=150, DT_h=0.5, T0=20
      real*8, parameter ::  tc=1
C
c     FIXED PARAMETERS FOR DOUBLE PLATE
      real*8, parameter ::  DT_c=0, H=1, th=1
c
      real*8 :: period
      real*8 :: Tmedh, Tmedc, Tmin
      real*8 :: Tmax_t, Tmedh_t, Tmedc_t, Tmin_t
c      
c     for non-linear geometry the original coords must be used     
      REAL(8), DIMENSION (100000, 2)::COORDS0
      COMMON /myblock/ COORDS0
      IF (KINC .EQ. 1) THEN
           COORDS0(NODE,1)= COORDS(2)
      END IF

      DT_ped=1-DT_h-DT_c
      Tmedh=Tmax-DT*DT_h
      Tmedc=Tmedh-DT*DT_ped
      Tmin=Tmedc-DT*DT_c
      if (TIME(2) .LE. t_ramp) then
          Tmax_t=((Tmax-T0)/t_ramp)*(TIME(2)) +T0
          Tmedh_t=((Tmedh-T0)/t_ramp)*(TIME(2)) +T0
          Tmedc_t=((Tmedc-T0)/t_ramp)*(TIME(2)) +T0
          Tmin_t=((Tmin-T0)/t_ramp)*(TIME(2)) +T0
      else
          Tmax_t=Tmax
          Tmedh_t=Tmedh
          Tmedc_t=Tmedc
          Tmin_t=Tmin
      end if    
c
c   PROVIDE U(1) field
      if (COORDS0(NODE,1) .GE. tc+H) then
         U(1)=Tmedh_t + 
     1     (Tmax_t-Tmedh_t)*(COORDS0(NODE,1)-tc-H)/th
      elseif (COORDS0(NODE,1) .GE. tc) then
         U(1)=Tmedc_t + 
     1     (Tmedh_t-Tmedc_t)*(COORDS0(NODE,1)-tc)/H
      else    
         U(1)=Tmin_t + 
     1     (Tmedc_t-Tmin_t)*COORDS0(NODE,1)/tc
      end if  
c          
      RETURN
      END
c
c


      SUBROUTINE DLOAD(F,KSTEP,KINC,TIME,NOEL,NPT,LAYER,KSPT,
     1 COORDS,JLTYP,SNAME)
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION TIME(2), COORDS (3)
      CHARACTER*80 SNAME
c      
      REAL, PARAMETER::  p_nom=0
      REAL, PARAMETER::  tc=1
c      
c     FIXED
      REAL, PARAMETER:: th=1, W=2.8, D_film=0.1, D_imp=0.5
      REAL, PARAMETER :: pi = 3.1415927,  gamma=60
      REAL S_film, S_imp, S_nom, gamma_rad, p_net
c
c     calculate p_net
      gamma_rad=gamma*pi/180
      S_nom=(tc+th)*W
      S_film=D_film*th/((cos(gamma_rad))**2)
      S_imp=D_imp*tc
      p_net=p_nom*S_nom/(S_nom-2*S_film-S_imp)
c      
      F=-p_net
c          
      RETURN
      END