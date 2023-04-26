      SUBROUTINE MPC(UE,A,JDOF,MDOF,N,JTYPE,X,U,UINIT,MAXDOF,LMPC,
     1               KSTEP,KINC,TIME,NT,NF,TEMP,FIELD,LTRAN,TRAN)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION UE(MDOF), A(MDOF,MDOF,N), JDOF(MDOF,N), X(6,N),
     1          U(MAXDOF,N), UINIT(MAXDOF,N), TIME(2), TEMP(NT,N),
     2          FIELD(NF,NT,N)
c      
      if (JTYPE.eq.1) then  
         JDOF(1,1)=3
         JDOF(1,2)=3     
         A(1,1,1) = 1
         A(1,1,2) =-1
      else if (JTYPE.eq.2) then
         JDOF(1,1)=3
         JDOF(1,2)=3     
         A(1,1,1) = 1
         A(1,1,2) = -1
      else if (JTYPE.eq.3) then
         JDOF(1,1)=1
         JDOF(1,2)=1     
         A(1,1,1) = 1
         A(1,1,2) =-1        
      else              
        JDOF(1,1)=1
        JDOF(2,1)=2
        JDOF(1,2)=1
        JDOF(2,2)=2        
        JDOF(1,3)=1
        A(1,1,1) = 1
        A(1,1,2) =-1
        A(1,1,3) =-1
        A(2,2,1) = 1
        A(2,2,2) =-1            
      end if
C  Give the value of dependent degree of freedom by solving the constraint equation about u1
c     print *, 'UE_default', UE(1)
      RETURN
      END