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
c     Maximum nominal traction at steady state of cycle
      REAL, PARAMETER::   p_nom=100
c      
      REAL, PARAMETER::   tc=1     
c      
      REAL, PARAMETER:: T_cycle=2, t_up=0.02, t_down=0.02, t_rest=0.5 
      REAL period, k, t_steady, p_t
c     Geometric dimensions      
      REAL, PARAMETER:: th=1, W=5, D_film=0.6, D_imp=0.6, film_angle=60
      REAL S_film, S_imp, S_nom, film_angle_rad, p_net
c
c     calculate p_net at steady state (maximum)
      film_angle_rad=film_angle*pi/180
      S_nom=(tc+th)*W
      S_film=D_film*th/((cos(film_angle_rad))**2)
      S_imp=D_imp*tc
      p_net=p_nom*S_nom/(S_nom-2*S_film-S_imp)
c         
      k=INT(TIME(2)/T_cycle)
      t_steady=T_cycle-t_up-t_down
        period=T_cycle+t_rest
        k=INT(TIME(2)/period)
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
c        
c   Provide F to ABAQUS
        F=-p_t
c          
      RETURN
      END
