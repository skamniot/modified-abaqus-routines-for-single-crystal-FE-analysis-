      SUBROUTINE MPC(UE,A,JDOF,MDOF,N,JTYPE,X,U,UINIT,MAXDOF,LMPC,
     1               KSTEP,KINC,TIME,NT,NF,TEMP,FIELD,LTRAN,TRAN)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION UE(MDOF), A(MDOF,MDOF,N), JDOF(MDOF,N), X(6,N),
     1          U(MAXDOF,N), UINIT(MAXDOF,N), TIME(2), TEMP(NT,N),
     2          FIELD(NF,NT,N)
C     added variables 
c      REAL, PARAMETER::  xmax=30, zmax =20
C     REAL, PARAMETER:: PRECISION = 1E-04, pi = 3.1415927     
c     
c     type 4,5        
      if (N.eq.4) then  
        JDOF(1,1)=1
        JDOF(2,1)=2
        JDOF(3,1)=3
        JDOF(1,2)=1
        JDOF(2,2)=2
        JDOF(3,2)=3        
        JDOF(1,3)=1
        JDOF(1,4)=3
        A(1,1,1) = 1
        A(1,1,2) =-1
        A(1,1,3) =-1       
        A(2,2,1) = 1
        A(2,2,2) =-1        
        A(3,3,1) = 1
        A(3,1,4) =-1 
c    type 3            
      elseif (N.eq.3) then               
        JDOF(1,1)=1
        JDOF(2,1)=2
        JDOF(3,1)=3
        JDOF(1,2)=1
        JDOF(2,2)=2
        JDOF(3,2)=3        
        JDOF(1,3)=1
        A(1,1,1) = 1
        A(1,1,2) =-1
        A(1,1,3) =-1      
        A(2,2,1) = 1
        A(2,2,2) =-1        
        A(3,3,1) = 1
        A(3,3,2) =-1     
      else if (JTYPE.eq.3) then
        JDOF(1,1)=3
        JDOF(1,2)=3     
        A(1,1,1) = 1
        A(1,1,2) =-1
      else if (JTYPE.eq.2) then
        JDOF(1,1)=2
        JDOF(1,2)=2   
        JDOF(2,1)=3
        JDOF(2,2)=3     
        A(1,1,1) = 1
        A(1,1,2) =-1
        A(2,2,1) = 1
        A(2,2,2) =-1
      else
        JDOF(1,1)=3
        JDOF(1,2)=3     
        A(1,1,1) = 1
        A(1,1,2) =-1  
      end if
C  Give the value of dependent degree of freedom by solving the constraint equation about u1
c      print *, 'UE_default', UE(1)
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
c     VARIABLE PARAMETERS
      REAL, PARAMETER::  DT=0, tc=1
c      
c     FIXED PARAMETERS
      REAL, PARAMETER::   Tmax=1100, Rp=0.5, H=1.5, th=1, Rc=0
      REAL DTh, DTp, DTc, Rh
C
      Rh=1-Rc-Rp
      DTh=DT*Rh
      DTp=DT*Rp
      DTc=DT*Rc
c      COORDS(2)=COORDS(2)+tc+H/2
      if (COORDS(2)-0.5 .GE. tc+H) then
          U(1)=Tmax-DTh + DTh*(COORDS(2)-tc-H)/th
      elseif (COORDS(2)-0.5 .GE. tc) then
          U(1)=Tmax-DTh-DTp + DTp*(COORDS(2)-tc)/H
      else    
          U(1)=Tmax-DTh-DTp-DTc +DTc*(COORDS(2))/tc
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
      REAL, PARAMETER::  p_nom=225
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