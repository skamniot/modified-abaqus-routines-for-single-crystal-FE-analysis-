c     The routine applies a temperature field that cycles with time 
c     on a double plate geometry
      SUBROUTINE  DISP(U,KSTEP,KINC,TIME,NODE,NOEL,JDOF,COORDS)
c
      INCLUDE 'ABA_PARAM.INC'
c
      DIMENSION U(3),TIME(2),COORDS(3)
c
*************************************************************************
c     Determine the type of displacement controlled loading:  
      ! 1 = ramp displacement
      ! 10 = cyclic displacement
      ! 2 = ramp uniform temperature
      ! 20 = cyclic uniform temperature
      ! 3 = ramp temperature field (here applied in double plate geometry)
      ! 30 = cyclic temperature field (here applied in double plate geometry)
      integer, parameter :: loadtype = 1
*************************************************************************
c    
      real*8, parameter ::  t_ramp=0.05    ! cases 1,2,3
c      
      real*8, parameter ::  T_cycle=0.25, t_up=0.1, t_down=0.1, t_rest=0.05  ! cases 10, 20, 30
c      
      real*8, parameter ::  Umax = 0.1   ! max displacement - cases 1,10
c      
      real*8, parameter ::  Temp_max = 400   ! max temperature - cases 2,20
c        
c     CHOOSE PARAMETERS FOR DOUBLE PLATE (for cases 3, 30)
      real*8, parameter ::  Tmax=1100, DT=100, DT_h=0.5
      real*8, parameter ::  lag_up=1, lag_down=1   ! out-of-phase factors for combined force-temp loading
      real*8, parameter ::  tc=1
c      
c     FIXED PARAMETERS FOR DOUBLE PLATE (for cases 3, 30)
      real*8, parameter ::  DT_c=0, T0=20, H=1.5, th=1
c
c     
      real*8 :: heatup, dwell, cooldown, rest, period, t_steady
      real*8 :: Tmedh, Tmedc, Tmin, DT_ped, COORD_R
      real*8 :: Tmax_t, Tmedh_t, Tmedc_t, Tmin_t
      integer :: k
c      
c     for non-linear geometry the original coords must be used     
      REAL(8), DIMENSION (40000, 1)::COORDS0
      COMMON /myblock/ COORDS0
      IF (KINC .EQ. 1) THEN
         COORDS0(NODE,1)= COORDS(2)
      END IF
c
c     calculation of loading time variables
      period=T_cycle+t_rest 
      t_steady=T_cycle-t_up-t_down
      heatup= t_up*lag_up
      cooldown=t_down*lag_down
      dwell=t_steady+(1-lag_up)*t_up
      rest=t_rest+(1-lag_down)*t_down
      k=INT(TIME(2)/period)
c
c
c       
      SELECT CASE(loadtype)
c          
      CASE(1) !ramp displacement
c
      if (TIME(2) .LE. t_ramp) then
          U(1)=(Umax/t_ramp)*TIME(2)
      else 
          U(1)=Umax
      endif
c
c
c      
      case(10) !cyclic displacement
c          
      if (TIME(2) .LE. k*period + t_up) then
          U(1)=(Umax/t_up)*(TIME(2)-k*period)
      else if (TIME(2) .LE. k*period + t_up+t_steady) then
          U(1)=Umax
      else if (TIME(2) .LE. k*period + t_up+t_steady+t_down) then
          U(1)=Umax-(Umax/t_down)*
     1                (TIME(2)-k*period-t_up-t_steady)
      else    
          U(1)=0.0
      end if
c      
c
c      
      case(2) !ramp temperature
c          
      if (TIME(2) .LE. t_ramp) then
          U(1)=(Temp_max/t_ramp)*TIME(2)
      else 
          U(1)=Temp_max
      endif
c
c
c      
      case(20) !cyclic temperature
c          
      if (TIME(2) .LE. k*period + heatup) then
          U(1)=(Temp_max/heatup)*(TIME(2)-k*period)
      else if (TIME(2) .LE. k*period + heatup+dwell) then
          U(1)=Temp_max
      else if (TIME(2) .LE. k*period + heatup+dwell+cooldown) then
          U(1)=Temp_max-(Temp_max/cooldown)*
     1                (TIME(2)-k*period-heatup-dwell)
      else    
          U(1)=0.0
      end if
c
c
c     
      case(3) ! ramp temperature field
c
      DT_ped=1-DT_h-DT_c
      Tmedh=Tmax-DT*DT_h
      Tmedc=Tmedh-DT*DT_ped
      Tmin=Tmedc-DT*DT_c
      if (TIME(2) .LE. t_ramp) then
          Tmax_t=((Tmax-T0)/t_ramp)*(TIME(2)-k*period) +T0
          Tmedh_t=((Tmedh-T0)/t_ramp)*(TIME(2)-k*period) +T0
          Tmedc_t=((Tmedc-T0)/t_ramp)*(TIME(2)-k*period) +T0
          Tmin_t=((Tmin-T0)/t_ramp)*(TIME(2)-k*period) +T0
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
c
c      
      case(30) ! cyclic temperature field
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
c
      case default      
      WRITE(*,*)"Not sure what type of displacement controlled loading"
      END SELECT
c
c
      RETURN
      END
