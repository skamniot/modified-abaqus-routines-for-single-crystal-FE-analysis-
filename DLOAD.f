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