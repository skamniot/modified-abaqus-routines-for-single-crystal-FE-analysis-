c     The routine applies a temperature field that cycles with time 
c     on a double wall geometry
      SUBROUTINE  DISP(U,KSTEP,KINC,TIME,NODE,NOEL,JDOF,COORDS)
c
      INCLUDE 'ABA_PARAM.INC'
c
      DIMENSION U(3),TIME(2),COORDS(3)
c     
c     VARIABLE PARAMETERS
      REAL, PARAMETER::  Tmax=1100, DT=250, DT_h=0.5
      REAL, PARAMETER::  lag_up=1, lag_down=1
      REAL, PARAMETER::  tc=1
c
c     FIXED PARAMETERS
      REAL, PARAMETER::  DT_c=0, T0=20, H=1, th=1
      REAL, PARAMETER::  T_cycle=2, t_up=0.02, t_down=0.02, t_rest=0.5
c
c     
      REAL heatup, dwell, cooldown, rest, period, t_steady
      REAL Tmedh, Tmedc, Tmin, DT_ped, COORD_R
      REAL Tmax_t, Tmedh_t, Tmedc_t, Tmin_t, k  
c      
c     for non-linear geometry the original coords must be used     
      REAL(8), DIMENSION (40000, 1)::COORDS0
      COMMON /myblock/ COORDS0
      IF (KINC .EQ. 1) THEN
         COORDS0(NODE,1)= COORDS(2)
      END IF
c
c
      period=T_cycle+t_rest 
      t_steady=T_cycle-t_up-t_down
      heatup= t_up*lag_up
      cooldown=t_down*lag_down
      dwell=t_steady+(1-lag_up)*t_up
      rest=t_rest+(1-lag_down)*t_down
      k=INT(TIME(2)/period)
c
      DT_ped=1-DT_h-DT_c
      Tmedh=Tmax-DT*DT_h
      Tmedc=Tmedh-DT*DT_ped
      Tmin=Tmedc-DT*DT_c
      if (TIME(2) .LE. k*period + heatup) then
          Tmax_t=((Tmax-T0)/heatup)*(TIME(2)-k*period) +T0
          Tmedh_t=((Tmedh-T0)/heatup)*(TIME(2)-k*period) +T0
          Tmedc_t=((Tmedc-T0)/heatup)*(TIME(2)-k*period) +T0
          Tmin_t=((Tmin-T0)/heatup)*(TIME(2)-k*period) +T0
      else if (TIME(2) .LE. k*period + heatup+dwell) then
          Tmax_t=Tmax
          Tmedh_t=Tmedh
          Tmedc_t=Tmedc
          Tmin_t=Tmin
      else if (TIME(2) .LE. k*period + heatup+dwell+cooldown) then
          Tmax_t=Tmax-((Tmax-T0)/cooldown)*
     1                (TIME(2)-k*period-heatup-dwell)
          Tmedh_t=Tmedh-((Tmedh-T0)/cooldown)*
     1               (TIME(2)-k*period-heatup-dwell)
          Tmedc_t=Tmedc-((Tmedc-T0)/cooldown)*
     1               (TIME(2)-k*period-heatup-dwell)
          Tmin_t=Tmin-((Tmin-T0)/cooldown)*
     1               (TIME(2)-k*period-heatup-dwell)
      else    
          Tmax_t=T0
          Tmedh_t=T0
          Tmedc_t=T0
          Tmin_t=T0
      end if
c
c   PROVIDE U(1) TO ABAQUS
c
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
