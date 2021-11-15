c     The routine applies uniform traction load that cycles with time
c     The traction is applied on wall surfaces that are reduced by holes
c     Therefore for a given nominal traction p_nom the code computes the 
c     net traction that needs to be applied on the net wall sections to 
c     transmit a certain force - useful in turbine blade applications 
c     where the load represents centrifugal stresses
      SUBROUTINE DLOAD(F,KSTEP,KINC,TIME,NOEL,NPT,LAYER,KSPT,
     1 COORDS,JLTYP,SNAME)
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION TIME(2), COORDS (3)
      CHARACTER*80 SNAME
c      
*************************************************************************
c     Determine loading type: 
!     1 = ramp traction   
!     10 = cyclic traction
!     2 = ramp traction on double plate model  
!     20 = cyclic traction on double plate model
      integer, parameter :: loadtype = 1
*************************************************************************
c     Maximum traction 
      REAL, PARAMETER::   p_nom=0.0
c      
      REAL, PARAMETER::   tc=1    ! for double plate model 
c      
      REAL, PARAMETER::  t_ramp=1    ! for loadtype = 1,2
c      
!     for lodtype = 10, 20
      REAL, PARAMETER::  T_cycle=2, t_up=0.02, t_down=0.02, t_rest=0.5 
c      
      REAL period, k, t_steady, p_t, p_net
c      
c     Geometric dimensions for double plate model     
      REAL, PARAMETER:: th=1, W=5, D_film=0.6, D_imp=0.6, film_angle=60
      REAL S_film, S_imp, S_nom, film_angle_rad
c
c
c     loading time variables
        k=INT(TIME(2)/T_cycle)
        t_steady=T_cycle-t_up-t_down
        period=T_cycle+t_rest
        k=INT(TIME(2)/period)
c
c     calculate p_net for double plate model
      film_angle_rad=film_angle*pi/180
      S_nom=(tc+th)*W
      S_film=D_film*th/((cos(film_angle_rad))**2)
      S_imp=D_imp*tc
      p_net=p_nom*S_nom/(S_nom-2*S_film-S_imp)
c
c
      SELECT CASE(loadtype)
c          
      CASE(1) !ramp traction
c
      if (TIME(2) .LE. t_ramp) then
          p_t=(p_nom/t_ramp)*TIME(2)
      else
          p_t=p_nom
      end if
      F=-p_t
c          
      CASE(10) !cyclic traction
c
      if (TIME(2) .LE. k*period + t_up) then
          p_t=(p_nom/t_up)*(TIME(2)-k*period)
      else if (TIME(2) .LE. k*period + t_up+t_steady) then
          p_t=p_nom
      else if (TIME(2) .LE. k*period + t_up+t_steady+t_down) then
          p_t=p_nom-(p_nom/t_down)*
     1        (TIME(2)-k*period-t_up-t_steady)
      else    
          p_t=0.0
      end if
      F=-p_t
c          
      CASE(2) !ramp traction on double plate  
c
      if (TIME(2) .LE. t_ramp) then
          p_t=(p_net/t_ramp)*TIME(2)
      else
          p_t=p_net
      end if
      F=-p_t
c    
      CASE(20) !ramp traction on double plate  
c 
      if (TIME(2) .LE. k*period + t_up) then
          p_t=(p_net/t_up)*(TIME(2)-k*period)
        else if (TIME(2) .LE. k*period + t_up+t_steady) then
          p_t=p_net
        else if (TIME(2) .LE. k*period + t_up+t_steady+t_down) then
          p_t=p_net-(p_net/t_down)*
     1        (TIME(2)-k*period-t_up-t_steady)
        else    
          p_t=0.0
        end if
      F=-p_t
c        
      case default      
      WRITE(*,*)"Not sure what type of displacement controlled loading"
      END SELECT
c      
c          
      RETURN
      END
