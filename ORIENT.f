      SUBROUTINE ORIENT(T,NOEL,NPT,LAYER,KSPT,COORDS,BASIS,
     1 ORNAME,NNODES,CNODES,JNNUM)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 ORNAME
C
      DIMENSION T(3,3),COORDS(3),BASIS(3,3),CNODES(3,NNODES)
      DIMENSION JNNUM(NNODES)
c
! Transformation matrix
c      REAL*8,intent(out) :: T(3,3)
c      
      REAL, PARAMETER:: pi = 3.141592653
c     Define euler rotation angles between global system and crystal system IN DEGREES     
      REAL, PARAMETER::  fii=45, theta=45, psii=0
      REAL  th, fi, psi
c
C     inversion of angles psi & fi to convert from our intrinsic rotation
c     and the extrinsic rotation considered by abaqus
      psi=fii*pi/180
      th=theta*pi/180
      fi=psii*pi/180
      T(1,1)=cos(fi)*cos(psi)-sin(fi)*cos(th)*sin(psi)
      T(1,2)=sin(fi)*cos(psi)+cos(fi)*cos(th)*sin(psi)
      T(1,3)=sin(th)*sin(psi)
      T(2,1)=-sin(psi)*cos(fi)-sin(fi)*cos(th)*cos(psi)
      T(2,2)=-sin(fi)*sin(psi)+cos(fi)*cos(th)*cos(psi)
      T(2,3)=sin(th)*cos(psi)
      T(3,1)=sin(fi)*sin(th)
      T(3,2)=-sin(th)*cos(fi)
      T(3,3)=cos(th)
c
      RETURN
      END