!Least Squares Density Minimisation
      subroutine kgndl2ET(curlFp,xNin,xDin,tau,burgerv,iphase,ne,ns,
     + screwplanes,jelem,kint,time,gndall,gndcut,gndmob)
c              
!      implicit real*8(a-h,o-z)
      implicit none
c
      integer :: i,m,n,ialloc,info
      real*8  :: det, costheta, sumtau
      real*8,parameter :: zero = 1.0e-6
c      
      integer,intent(in):: ne,iphase,jelem,kint,ns,screwplanes
      real*8,intent(in) :: time(2)
      real*8,intent(in) :: curlFp(3,3),xNin(ne,3),xDin(ne,3),tau(ne),
     + burgerv(ne)   
c      
      real*8,dimension(3,3) :: xrot,btdyad,bsdyad
      real*8,dimension(3)   :: tempn,temps,tempt
      real*8,dimension(9)   :: btvec,bsvec,gv
      real*8,dimension(ne+ns,3) :: xLine
      real*8, dimension(ne,ne+ns) :: tdotn
c      
      !real*8,dimension(:,:),allocatable :: A,A1,A2,Ainv
      !real*8,dimension(:),allocatable :: xvec,absgnds
      real*8,dimension(9,ne+ns) :: A
      real*8,dimension(ne+ns,9) :: Ainv            
      real*8,dimension(9,9) :: A1, A2      
      real*8,dimension(ne+ns):: xvec
      !real*8,dimension(:),allocatable :: xvec
      real*8, dimension(ne) :: gnde
      real*8, dimension(ns) :: gnds
c      
      integer,dimension(ns) :: screw
      integer, dimension(ne) :: stype
      integer, dimension(ne, screwplanes-1) :: sgroup
c     
      character (len=*), parameter :: fmt2 = "(24(' ',(I2,1X)))",
     + fmt3="(3(' ',(ES11.3,1X)))",
     + fmt9 = "(9(' ',(ES11.3,1X)))",fmt33 = "(33(' ',(ES11.3,1X)))"
c          
      real*8,intent(out) :: gndall(ne+ns), gndcut(ne), gndmob    
c      
      integer :: j
c      
      !EXTERNAL DGELSD
c    
      !select case(iphase) 
      !case (0); ns = 9; ne=12; !HCP       
      !case (1); ns = 4; ne=24 !BCC     
      !case (2); ns = 6; ne=12 !FCC    
      !case (4); ns = 2; ne=7 !olivine 
      !end select 
c     
      !allocate(absgnds(ns), STAT=ialloc)
c      
      select case(iphase) !build up the screw slip systems
          case(0) !HCP
          !<a> slip
          screw(1) = 1 
          screw(2) = 2
          screw(3) = 3
          !<c+a>
          screw(4) = 7
          screw(5) = 8
          screw(6) = 9
          screw(7) = 10
          screw(8) = 11
          screw(9) = 12          
c          
          case(1) ! BCC
          !a/2<111>
          screw(1) = 1           
          screw(2) = 2
          screw(3) = 3
          screw(4) = 6
c          
          case(2) ! FCC
          screw(1) = 1           
          screw(2) = 2
          screw(3) = 3
          screw(4) = 4
          screw(5) = 6
          screw(6) = 8               
c          
          case(4) ! olivine
          screw(1) = 1
          screw(2) = 3       
c          
      end select
c              
c      
      gnde = 0.; gnds = 0.; gndall = 0.0; gndcut = 0.0; gndmob =0.0
c      
c            
      !Compute the geometric dislocation tensor
      gv = reshape(curlFp,(/9/))!This happens by column.
c                 
      !===BEGIN LONG IF
      if(maxval(abs(gv)) <= zero) then 
      !Trying to handle the situation when G = 0.
         gnde = 0.; gnds = 0.
      else         
c            
      m = 9; n = ne+ns      
      !allocate(A(m,n),Ainv(n,m),xvec(n),STAT=ialloc)
      A = 0.; Ainv = 0.; xvec=0. 
c      
      !m < n  !right inverse
       !allocate(A1(m,m),A2(m,m),STAT=ialloc)
       A1=0.; A2=0.
c      
      !Construct the matrix of dyadics
      !Edge
c     
      do i = 1,ne	      
          tempn = xNin(i,:)
          temps = xDin(i,:)   
          CALL KVECPROD(temps,tempn,tempt)
c                    
          btdyad = spread(tempt,2,3)*spread(temps,1,3)*burgerv(i) 
c         
          A(:,i) = reshape(btdyad,(/9/))              
          xLine(i,:) = tempt     
      end do
c      
      do i = 1,ns          
          j = screw(i) 
          temps = xDin(j,:)  
          tempt = temps
c        
          btdyad = spread(tempt,2,3)*spread(temps,1,3)*burgerv(i)     
c         
          A(:,ne+i) = reshape(btdyad,(/9/))    
          xLine(ne+i,:) = tempt
c
      end do
c      
c           
!      if ((jelem == 7910 .or. jelem == 8010) .and. kint == 6) then
!         write(6,*) "A, ne",ubound(A,1),ubound(A,2),ne
!         write(6,fmt33) (A(k,:),k=1,ubound(A,1))  !Print rows
!         write(6,*)  
!      end if 
c
!!      call ksvd(A,m,n,Ainv)
!      
c     
      ! m < n right inverse
       A1 = matmul(A,transpose(A)) ![9x9] = [9xn][nx9]
       call lapinverse(A1,m,info,A2)
       if(info /= 0) write(*,*) "inverse failure: A1 in kgndl2"
       Ainv = matmul(transpose(A),A2) ![nx9]=[nx9][9x9]
c     
      !three-ifs solution
      xvec = matmul(Ainv,gv) ![nx1] = [nx9][9x1]
c      
      do i=1,ne+ns !ubound(xvec,1)
         if (xvec(i) /= xvec(i)) then       
             xvec(i) = 0.0  
           write(*,*) "error in kgndl2 jelem,kint,time",jelem,kint,time
         end if
      end do
c  
c      
c         
c      
      !Density to output
      gndall = sqrt(xvec*xvec) !sqrt(gnde*gnde +gnds*gnds) 
c      
      !deallocate(A,A1,A2,Ainv,xvec)  
c      
      where(gndall > 1.0E7) gndall = 1.0E7 
      !catching infinities. 1.0e7 for microns, 1.0e19 for metres.
c       
c            
      do i =1,ne           
          tempn = xNin(i,:)
          do j=1,n
               tempt = xLine(j,:)          
               CALL KDOTPROD(tempt,tempn,costheta)
               tdotn(i,j) = costheta
          end do
      end do
c      
      gndcut= matmul(tdotn,gndall)
c      
      gnde = gndall(1:ne)
      gnds = gndall(ne:ne+ns)
c           
c      
      !do i =1,ne
      !    j = stype(i) 
      !    
      !    sumtau = tau(i)
      !    do n=1,screwplanes-1 ! other possible slip planes which screw could glide on
      !        sumtau = sumtau + tau(sgroup(i,n))              
      !    end do
      !    
      !    gndmob(i) = gnde(i) + gnds(j)*tau(i)/sumtau
      !end do
c      
      !Ben's experimental data storage groups
      !absgnde = sqrt(gnde*gnde); absgnds = sqrt(gnds*gnds)
      !if(iphase == 0) then
      !!edge
      !abasedge = sum(absgnde(1:3))
      !aprismedge = sum(absgnde(4:6))
      !apyramedge = sum(absgnde(7:12))
      !capyramedge = sum(absgnde(13:18))
      !!screw
      !ascrew = sum(absgnds(1:3))
      !capyramscrew = sum(absgnds(4:9))
      !end if
      !
      !===END LONG IF
c            
      end if
c      
c      
c      
      return
      end
     