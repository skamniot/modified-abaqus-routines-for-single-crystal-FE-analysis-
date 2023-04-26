      SUBROUTINE  DISP(U,KSTEP,KINC,TIME,NODE,NOEL,JDOF,COORDS)
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION U(3),TIME(2),COORDS(3)
c
c     VARIABLE PARAMETERS
      REAL, PARAMETER::  DT=100, tc=1
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