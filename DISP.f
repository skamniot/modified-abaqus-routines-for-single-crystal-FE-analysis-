      SUBROUTINE  DISP(U,KSTEP,KINC,TIME,NODE,NOEL,JDOF,COORDS)
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION U(3),TIME(2),COORDS(3)
c     
c     VARIABLE PARAMETERS
      REAL, PARAMETER::   Tmax=1000, DT=300, DT_ped=0.5
      REAL, PARAMETER::   tc=1
c     
c     FIXED PARAMETERS
      REAL, PARAMETER:: DT_c=0, T0=20, H=1, th=1
c    
      REAL Tmedh, Tmedc, Tmin
c      
c     for non-linear geometry      
c      REAL(8), DIMENSION (30000, 1)::COORDS0
c      COMMON /myblock/ COORDS0
c      IF (KINC .EQ. 1) THEN
c         COORDS0(NODE,1)= COORDS(2)
c      END IF
c
      Tmedh=Tmax-DT*(1-DT_ped-DT_c)
      Tmedc=Tmedh-DT*DT_ped
      Tmin=Tmedc-DT*DT_c  
c      
      if (COORDS(2) .GE. tc+H) then
         U(1)=Tmedh + 
     1     (Tmax-Tmedh)*(COORDS(2)-tc-H)/th
      elseif (COORDS(2) .GE. tc) then
         U(1)=Tmedc + 
     1     (Tmedh-Tmedc)*(COORDS(2)-tc)/H
      else    
         U(1)=Tmin + 
     1     (Tmedc-Tmin)*COORDS(2)/tc
      end if    
c          
      RETURN
      END