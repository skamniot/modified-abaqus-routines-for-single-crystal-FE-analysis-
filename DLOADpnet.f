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
!     2 = ramp traction for double plate model  
!     20 = cyclic traction for double plate model
      integer, parameter :: loadtype = 20
*************************************************************************
c     net traction at minimum double wall section 
      REAL, PARAMETER::   p_net=0
      
c     porosity for double plate model     
      REAL, PARAMETER:: fi=0.1
c      
      REAL, PARAMETER::  t_ramp=1    ! for loadtype = 1,2
c      
!     for lodtype = 10, 20
      real*8, parameter ::  T_cycle=10000.0, t_up=500, t_down=500, t_rest=5000.0
c      
      REAL period, k, t_steady, p_t, p_applied, factor
c
c     loading time variables
        k=INT(TIME(2)/T_cycle)
        t_steady=T_cycle-t_up-t_down
        period=T_cycle+t_rest
        k=INT(TIME(2)/period)
c
c     calculate p_net for double plate model
      if (fi .EQ. 0.02) then
        factor=0.957
      elseif (fi .EQ. 0.05) then
        factor=0.944
      elseif (fi .EQ. 0.1) then
        factor=0.934
      elseif (fi .EQ. 0.2) then
        factor=0.923
      elseif (fi .EQ. 0.3) then
        factor=0.917
      end if
      
      p_applied=p_net*factor
      
c
      SELECT CASE(loadtype)
c          
      CASE(1) !ramp traction
c
      if (TIME(2) .LE. t_ramp) then
          p_t=(p_net/t_ramp)*TIME(2)
      else
          p_t=p_net
      end if
      F=-p_t
c          
      CASE(10) !cyclic traction
c
      if (TIME(2) .LE. k*period + t_up) then
          p_t=(p_net/t_up)*(TIME(2)-k*period)
      else if (TIME(2) .LE. k*period + t_up+t_steady) then
          p_t=p_net
      else if (TIME(2) .LE. k*period + t_up+t_steady+t_down) then
          p_t=p_net-(p_nom/t_down)*
     1        (TIME(2)-k*period-t_up-t_steady)
      else    
          p_t=0.0
      end if
      F=-p_t
c          
      CASE(2) !ramp traction on double plate  
c  
      if (TIME(2) .LE. t_ramp) then
          p_t=(p_applied/t_ramp)*TIME(2)
      else
          p_t=p_applied
      end if
      F=-p_t
c    
      CASE(20) !ramp traction on double plate  
c 
      if (TIME(2) .LE. k*period + t_up) then
          p_t=(p_applied/t_up)*(TIME(2)-k*period)
        else if (TIME(2) .LE. k*period + t_up+t_steady) then
          p_t=p_applied
        else if (TIME(2) .LE. k*period + t_up+t_steady+t_down) then
          p_t=p_applied-(p_net/t_down)*
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
