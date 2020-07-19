!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE tree_module

      use setup_module

      IMPLICIT NONE

! r8 is 8-byte (double precision) real

      INTEGER,PARAMETER :: r8=SELECTED_REAL_KIND(12)

! global variables for taylor expansions

      INTEGER :: torderlim!,torder
      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:) :: cf,cf1,cf2
      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:,:,:) :: b1
      REAL(KIND=r8) :: xyzminmax(6)
      REAL(KIND=r8) :: strs1,strs2,strs3,strs4,strs5,strs6

! global variables for Tricubic

      INTEGER :: trcord,trcsize 
      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:,:) :: BinvT,BinvTM ! Transpose of the inverse of B
      REAL(KIND=r8),PARAMETER :: bfac = 0.037037037037037037037037_r8
! global variables to track tree levels 
 
      INTEGER :: minlevel,maxlevel!,lattlim,maxparnode
      
! global variables used when computing potential/force

      INTEGER :: orderoffset
      REAL(KIND=r8),DIMENSION(3) :: tarpos
      REAL(KIND=r8) :: thetasq!,mactheta

! global variables for postition and charge storage

      INTEGER,ALLOCATABLE,DIMENSION(:)  :: orderarr
      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:) :: txx,tyy,tzz,tchge

! node pointer and node type declarations

      TYPE tnode_pointer
           TYPE(tnode), POINTER :: p_to_tnode
      END TYPE tnode_pointer
      TYPE tnode
           INTEGER          :: numpar,ibeg,iend
           REAL(KIND=r8)    :: x_min,y_min,z_min
           REAL(KIND=r8)    :: x_max,y_max,z_max
           REAL(KIND=r8)    :: x_mid,y_mid,z_mid
           REAL(KIND=r8)    :: radius,sqradius,aspect
           REAL(KIND=r8)    :: dx,dy,dz,rdx,rdy,rdz
           INTEGER          :: level,num_children,exist_ms
           REAL(KIND=r8)    :: ccp(1:16,1:3)
           REAL(KIND=r8),DIMENSION(:,:,:),POINTER :: ms
           REAL(KIND=r8),DIMENSION(:),POINTER :: trcms
           REAL(KIND=r8),DIMENSION(:),POINTER :: trcmsx
           REAL(KIND=r8),DIMENSION(:),POINTER :: trcmsy
           REAL(KIND=r8),DIMENSION(:),POINTER :: trcmsz
           TYPE(tnode_pointer), DIMENSION(8) :: child
      END TYPE tnode

      save txx,tyy,tzz,tchge
      save strs1,strs2,strs3,strs4,strs5,strs6

      CONTAINS
!!!!!!!!!!!!!!!
      SUBROUTINE alloc_tree_arrays(idnode)
                      
      IMPLICIT NONE
!
! SETUP allocates and initializes arrays needed for the Taylor expansion.
! Also, global variables are set and the Cartesian coordinates of
! the smallest box containing the particles is determined. The particle
! postions and charges are copied so that they can be restored upon exit.
!
      INTEGER,INTENT(IN) :: idnode 

! local variables

      INTEGER :: err,i,j
      REAL(KIND=r8) :: t1

! global integers and reals:  TORDER, TORDERLIM and THETASQ

      orderoffset=1
      torderlim=torder+orderoffset
      thetasq=mactheta*mactheta

      trcord = 3
      trcsize = (trcord+1)**3

! allocate Tricubic variables

      ALLOCATE(BinvT(1:trcsize,1:trcsize),
     x         BinvTM(1:trcsize,1:trcsize), STAT=err)
      IF(err .NE. 0) THEN
        WRITE(6,*) 'Error allocating Tricubic B-inverse matrix'
        STOP
      END IF

! Read in B-inverse matrix
      IF(idnode .EQ. 0) THEN 
         OPEN(unit=80,file='Bnew2Inv_Matrix.txt')
         READ(80,*)BinvT ! Read in the tranpose of B-inverse 
!                     !((BinvT(i,j), j=1,trcsize),i=1,trcsize)
         BinvT = bfac*BinvT
         CLOSE(80)
      END IF

      IF(idnode .EQ. 0) THEN 
         OPEN(unit=81,file='BinvMatrix.txt')
         READ(81,*)BinvTM ! Read in the tranpose of B-inverse
         CLOSE(81)
      END IF

!      DO i=1,trcsize
!         WRITE(50,*)BinvT(i,:)
!      END DO

! allocate global Taylor expansion variables

      ALLOCATE(cf(0:torder),cf1(torder+1),cf2(torder+1),
     x        b1(-2:torderlim,-2:torderlim,-2:torderlim), STAT=err)
      IF (err .NE. 0) THEN
         WRITE(6,*) 'Error allocating Taylor variables! '
         STOP
      END IF

! initialize arrays for Taylor sums and coeffs

      DO i=0,torder
         cf(i)=-(REAL(i,KIND=r8)+1.0_r8)
      END DO
      DO i=1,torderlim
         t1=1.0_r8/REAL(i,KIND=r8)
         cf1(i)=2.0_r8-t1!-2.0_r8
         cf2(i)=t1-1.0_r8
      END DO

! Allocate arrays and copy variables. Also create and initialize orderarr.

      ALLOCATE(txx(mxatms),tyy(mxatms),tzz(mxatms), 
     x          tchge(mxatms),orderarr(mxatms),STAT=err)
      IF (err .NE. 0) THEN
         WRITE(6,*) 'Error allocating copy variables! '
         STOP
      END IF    

      RETURN
      END SUBROUTINE alloc_tree_arrays 
!!!!!!!!!!!!!
      RECURSIVE SUBROUTINE CREATE_TREE(p,ibeg,iend,x,y,z,q,shrink,
     x                                 xyzmm,level,arrdim)
      IMPLICIT NONE
!
! CREATE_TREE recursively create the tree structure. Node P is
! input, which contains particles indexed from IBEG to IEND. After
! the node parameters are set subdivision occurs if IEND-IBEG+1 > MAXPARNODE.
! Real array XYZMM contains the min and max values of the coordinates
! of the particle in P, thus defining the box.   
!
      TYPE(tnode),POINTER :: p
      INTEGER,INTENT(IN) :: ibeg,iend,shrink,level,arrdim
      REAL(KIND=r8),DIMENSION(arrdim),INTENT(INOUT) :: x,y,z,q
      REAL(KIND=r8),DIMENSION(6),INTENT(IN) :: xyzmm

! local variables

      REAL(KIND=r8) :: x_mid,y_mid,z_mid,xl,yl,zl,lmax,t1,t2,t3
      INTEGER, DIMENSION(8,2) :: ind
      REAL(KIND=r8), DIMENSION(6,8) :: xyzmms
      INTEGER :: i,j,limin,limax,err,loclev,numposchild
      REAL(KIND=r8), DIMENSION(6) ::  lxyzmm
     
! allocate pointer 

      ALLOCATE(p,STAT=err)
      IF (err .NE. 0) THEN
         WRITE(6,*) 'Error allocating pointer! '
         STOP
      END IF

! set node fields: number of particles, exist_ms
! and xyz bounds 

      p%numpar=iend-ibeg+1
      p%exist_ms=0

      IF (shrink .EQ. 1) THEN   
         p%x_min=MINVAL(x(ibeg:iend))
         p%x_max=MAXVAL(x(ibeg:iend))
         p%y_min=MINVAL(y(ibeg:iend))
         p%y_max=MAXVAL(y(ibeg:iend))
         p%z_min=MINVAL(z(ibeg:iend))
         p%z_max=MAXVAL(z(ibeg:iend))
      ELSE
         p%x_min=xyzmm(1)
         p%x_max=xyzmm(2)
         p%y_min=xyzmm(3)
         p%y_max=xyzmm(4)
         p%z_min=xyzmm(5)
         p%z_max=xyzmm(6)        
      END IF

! compute aspect ratio

      xl=p%x_max-p%x_min
      yl=p%y_max-p%y_min
      zl=p%z_max-p%z_min

      lmax=MAX(xl,yl,zl)
      t1=lmax
      t2=MIN(xl,yl,zl)
      IF (t2 .NE. 0.0_r8) THEN
         p%aspect=t1/t2
      ELSE
         p%aspect=0.0_r8
      END IF

! midpoint coordinates , RADIUS and SQRADIUS 

      p%x_mid=(p%x_max+p%x_min)/2.0_r8
      p%y_mid=(p%y_max+p%y_min)/2.0_r8
      p%z_mid=(p%z_max+p%z_min)/2.0_r8
      t1=p%x_max-p%x_mid
      t2=p%y_max-p%y_mid
      t3=p%z_max-p%z_mid
      p%sqradius=t1*t1+t2*t2+t3*t3
      p%radius=SQRT(p%sqradius)

! side length, inverse of side length
      p%dx = xl 
      p%dy = yl 
      p%dz = zl 

      p%rdx = 1.0_r8/xl
      p%rdy = 1.0_r8/yl
      p%rdz = 1.0_r8/zl

! sixteen cluster evaluation points
      p%ccp(1,1)=p%x_min
      p%ccp(1,2)=p%y_min
      p%ccp(1,3)=p%z_min
      p%ccp(2,1)=p%x_max
      p%ccp(2,2)=p%y_min
      p%ccp(2,3)=p%z_min
      p%ccp(3,1)=p%x_min
      p%ccp(3,2)=p%y_max
      p%ccp(3,3)=p%z_min
      p%ccp(4,1)=p%x_max
      p%ccp(4,2)=p%y_max
      p%ccp(4,3)=p%z_min

      p%ccp(5,1)=p%x_min
      p%ccp(5,2)=p%y_min
      p%ccp(5,3)=p%z_max
      p%ccp(6,1)=p%x_max
      p%ccp(6,2)=p%y_min
      p%ccp(6,3)=p%z_max
      p%ccp(7,1)=p%x_min
      p%ccp(7,2)=p%y_max
      p%ccp(7,3)=p%z_max
      p%ccp(8,1)=p%x_max
      p%ccp(8,2)=p%y_max
      p%ccp(8,3)=p%z_max

      p%ccp(9,1) =p%x_min+0.25_r8*xl
      p%ccp(9,2) =p%y_min+0.25_r8*yl
      p%ccp(9,3) =p%z_min+0.25_r8*zl
      p%ccp(10,1)=p%x_max-0.25_r8*xl
      p%ccp(10,2)=p%ccp(9,2)
      p%ccp(10,3)=p%ccp(9,3)
      p%ccp(11,1)=p%ccp(9,1) 
      p%ccp(11,2)=p%y_max-0.25_r8*yl
      p%ccp(11,3)=p%ccp(9,3)
      p%ccp(12,1)=p%ccp(10,1)
      p%ccp(12,2)=p%ccp(11,2)
      p%ccp(12,3)=p%ccp(9,3)

      p%ccp(13,1)=p%ccp(9,1)
      p%ccp(13,2)=p%ccp(9,2)
      p%ccp(13,3)=p%z_max-0.25_r8*zl
      p%ccp(14,1)=p%ccp(10,1)
      p%ccp(14,2)=p%ccp(9,2)
      p%ccp(14,3)=p%ccp(13,3)
      p%ccp(15,1)=p%ccp(9,1)
      p%ccp(15,2)=p%ccp(11,2)
      p%ccp(15,3)=p%ccp(13,3)
      p%ccp(16,1)=p%ccp(10,1)
      p%ccp(16,2)=p%ccp(11,2)
      p%ccp(16,3)=p%ccp(13,3)


! set particle limits, tree level of node, and nullify children pointers

      p%ibeg=ibeg
      p%iend=iend
      p%level=level
      IF (maxlevel .LT. level) THEN
         maxlevel=level
      END IF
      p%num_children=0
      DO i=1,8
         NULLIFY(p%child(i)%p_to_tnode)
      END DO

      IF (p%exist_ms .EQ. 0) THEN
          ALLOCATE(p%ms(0:torder,0:torder,0:torder),STAT=err)
          IF (err .NE. 0) THEN
             WRITE(6,*) 'Error allocating node moments! '
             STOP
          END IF
          CALL COMP_MS(p,x,y,z,q,arrdim)

          ALLOCATE(p%trcms(1:trcsize),p%trcmsx(1:trcsize), 
     x           p%trcmsy(1:trcsize),p%trcmsz(1:trcsize),STAT=err)
          IF (err .NE. 0) THEN
             WRITE(6,*) 'Error allocating node moments! '
             STOP
          END IF
          CALL TRICUBIC_MOMENTS(p,x,y,z,q,arrdim)

          p%exist_ms=1
      END IF   

      IF (p%numpar .GT. maxparnode) THEN
!
! set IND array to 0 and then call PARTITION routine.  IND array holds indices
! of the eight new subregions. Also, setup XYZMMS array in case SHRINK=1
!
         xyzmms(1,1)=p%x_min
         xyzmms(2,1)=p%x_max
         xyzmms(3,1)=p%y_min
         xyzmms(4,1)=p%y_max
         xyzmms(5,1)=p%z_min
         xyzmms(6,1)=p%z_max
         ind(1,1)=ibeg
         ind(1,2)=iend
         x_mid=p%x_mid
         y_mid=p%y_mid
         z_mid=p%z_mid

         CALL PARTITION_8(x,y,z,q,xyzmms,xl,yl,zl,lmax,numposchild,
     x                    x_mid,y_mid,z_mid,ind,arrdim)
!
! create children if indicated and store info in parent
!
         loclev=level+1
         DO i=1,numposchild
            IF (ind(i,1) .LE. ind(i,2)) THEN
               p%num_children=p%num_children+1
               lxyzmm=xyzmms(:,i)
               CALL CREATE_TREE(p%child(p%num_children)%p_to_tnode,
     x                          ind(i,1),ind(i,2),x,y,z,q,shrink,
     x                          lxyzmm,loclev,arrdim)
            END IF
         END DO
      ELSE
         IF (level .LT. minlevel) THEN
            minlevel=level
         END IF
      END IF   

      END SUBROUTINE CREATE_TREE      
!!!!!!!!!!!!!!!
      SUBROUTINE PARTITION_8(x,y,z,q,xyzmms,xl,yl,zl,lmax,numposchild,
     x                       x_mid,y_mid,z_mid,ind,arrdim)
      IMPLICIT NONE

! PARTITION_8 determines the particle indices of the eight sub boxes
! containing the particles after the box defined by particles I_BEG
! to I_END is divided by its midpoints in each coordinate direction.
! The determination of the indices is accomplished by the subroutine
! PARTITION. A box is divided in a coordinate direction as long as the
! resulting aspect ratio is not too large. This avoids the creation of
! "narrow" boxes in which Talyor expansions may become inefficient.
! On exit the INTEGER array IND (dimension 8 x 2) contains
! the indice limits of each new box (node) and NUMPOSCHILD the number 
! of possible children.  If IND(J,1) > IND(J,2) for a given J this indicates
! that box J is empty.
!
      INTEGER, INTENT(IN) :: arrdim
      REAL(KIND=r8),DIMENSION(arrdim),INTENT(INOUT) :: x,y,z,q
      INTEGER, DIMENSION(8,2),INTENT(INOUT) :: ind
      REAL(KIND=r8),DIMENSION(6,8),INTENT(INOUT) :: xyzmms
      REAL(KIND=r8), INTENT(IN) :: x_mid,y_mid,z_mid,xl,yl,zl,lmax
      INTEGER,INTENT(INOUT) :: numposchild

! local variables

      INTEGER :: temp_ind,i
      REAL(KIND=r8) :: critlen

      numposchild=1
      critlen=lmax/sqrt(2.0_r8)

      IF (xl .GE. critlen) THEN
         CALL PARTITION(x,y,z,q,orderarr,ind(1,1),ind(1,2),
     x                  x_mid,temp_ind,arrdim)

         ind(2,1)=temp_ind+1
         ind(2,2)=ind(1,2)
         ind(1,2)=temp_ind
         xyzmms(:,2)=xyzmms(:,1)
         xyzmms(2,1)=x_mid
         xyzmms(1,2)=x_mid
         numposchild=2*numposchild
      END IF 
 
      IF (yl .GE. critlen) THEN
         DO i=1,numposchild
            CALL PARTITION(y,x,z,q,orderarr,ind(i,1),ind(i,2),
     x                   y_mid,temp_ind,arrdim)
            ind(numposchild+i,1)=temp_ind+1
            ind(numposchild+i,2)=ind(i,2)
            ind(i,2)=temp_ind
            xyzmms(:,numposchild+i)=xyzmms(:,i)
            xyzmms(4,i)=y_mid
            xyzmms(3,numposchild+i)=y_mid
         END DO
         numposchild=2*numposchild
      END IF

      IF (zl .GE. critlen) THEN
         DO i=1,numposchild
            CALL PARTITION(z,x,y,q,orderarr,ind(i,1),ind(i,2),
     x                     z_mid,temp_ind,arrdim)
            ind(numposchild+i,1)=temp_ind+1
            ind(numposchild+i,2)=ind(i,2)
            ind(i,2)=temp_ind
            xyzmms(:,numposchild+i)=xyzmms(:,i)
            xyzmms(6,i)=z_mid
            xyzmms(5,numposchild+i)=z_mid
         END DO
         numposchild=2*numposchild
      END IF

      RETURN 
      END SUBROUTINE PARTITION_8
!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE TREE_COMPFP_PBC(p,x,y,z,q,tpeng,tvir,boxl,epsq,
     x                            arrdim,idnode,mxnode)

      use config_module
      use utility_module

      IMPLICIT NONE

! TREE_COMPF is the driver routine which calls COMPF_TREE for each
! particle, setting the global variable TARPOS before the call. 
! The current target particle's x coordinate and charge are changed
! so that it does not interact with itself. P is the root node of the tree. 
!

      INTEGER,INTENT(IN) :: arrdim,idnode,mxnode
      TYPE(tnode),POINTER :: p  
      REAL(KIND=r8),DIMENSION(arrdim),INTENT(INOUT) :: x,q
      REAL(KIND=r8),DIMENSION(arrdim),INTENT(IN) :: y,z
      REAL(KIND=r8),INTENT(INOUT) :: tpeng,tvir
      REAL(KIND=r8),INTENT(IN)    :: boxl(3),epsq

! local variables

      INTEGER :: i,j,ii,nx,ny,nz,atm0,atm1
      REAL(KIND=r8),DIMENSION(3) :: f,flocal,fff
      REAL(KIND=r8) :: penglocal,peng,vir,virlocal,tempx,tempq
      REAL(KIND=r8) :: boxlx,boxly,boxlz,nxboxl,nyboxl,nzboxl
      REAL(KIND=r8) :: fx,fy,fz

      boxlx=cell(1); boxly=cell(5); boxlz=cell(9)

      atm0=(idnode*arrdim)/mxnode+1
      atm1=((idnode+1)*arrdim)/mxnode

      tpeng=0.0_r8; tvir=0.0_r8
      fff=0.0_r8
      DO i=atm0,atm1
         f=0.0_r8
         peng=0.0_r8
         vir=0.0_r8
         strs1=0.0_r8; strs2=0.0_r8; strs3=0.0_r8
         strs4=0.0_r8; strs5=0.0_r8; strs6=0.0_r8
 
         tarpos(1)=x(i)
         tarpos(2)=y(i)
         tarpos(3)=z(i)

         DO nz=-lattlim,-1
            nzboxl=REAL(nz,KIND=r8)*boxlz

            DO ny=-lattlim,lattlim
               nyboxl=REAL(ny,KIND=r8)*boxly

               DO nx=-lattlim,lattlim
                  nxboxl=REAL(nx,KIND=r8)*boxlx                  

                 CALL COMPFP_TREE_PBC(p,penglocal,virlocal,flocal,
     x                  x,y,z,q,arrdim,nxboxl,nyboxl,nzboxl,idnode)

                  f=f+flocal
                  
                  peng=peng+penglocal     
                  vir=vir+virlocal
                  
               END DO
            END DO
         END DO

         DO nz=1,lattlim
            nzboxl=REAL(nz,KIND=r8)*boxlz

            DO ny=-lattlim,lattlim
               nyboxl=REAL(ny,KIND=r8)*boxly

               DO nx=-lattlim,lattlim
                  nxboxl=REAL(nx,KIND=r8)*boxlx

                  CALL COMPFP_TREE_PBC(p,penglocal,virlocal,flocal,
     x                  x,y,z,q,arrdim,nxboxl,nyboxl,nzboxl,idnode)

                  f=f+flocal

                  peng=peng+penglocal
                  vir=vir+virlocal

               END DO
            END DO
         END DO

         nz=0
         nzboxl=REAL(nz,KIND=r8)*boxlz

         DO ny=-lattlim,-1
            nyboxl=REAL(ny,KIND=r8)*boxly

            DO nx=-lattlim,lattlim
               nxboxl=REAL(nx,KIND=r8)*boxlx

               CALL COMPFP_TREE_PBC(p,penglocal,virlocal,flocal,
     x               x,y,z,q,arrdim,nxboxl,nyboxl,nzboxl,idnode)

               f=f+flocal

               peng=peng+penglocal
               vir=vir+virlocal

            END DO
         END DO

         DO ny=1,lattlim
            nyboxl=REAL(ny,KIND=r8)*boxly

            DO nx=-lattlim,lattlim
               nxboxl=REAL(nx,KIND=r8)*boxlx

               CALL COMPFP_TREE_PBC(p,penglocal,virlocal,flocal,
     x               x,y,z,q,arrdim,nxboxl,nyboxl,nzboxl,idnode)
     
               f=f+flocal

               peng=peng+penglocal
               vir=vir+virlocal

            END DO
         END DO

         nz=0; ny=0
         nzboxl=REAL(nz,KIND=r8)*boxlz
         nyboxl=REAL(ny,KIND=r8)*boxly

         DO nx=-lattlim,-1
            nxboxl=REAL(nx,KIND=r8)*boxlx


            CALL COMPFP_TREE_PBC(p,penglocal,virlocal,flocal,
     x            x,y,z,q,arrdim,nxboxl,nyboxl,nzboxl,idnode)
     
            f=f+flocal

            peng=peng+penglocal
            vir=vir+virlocal

         END DO

         DO nx=1,lattlim
            nxboxl=REAL(nx,KIND=r8)*boxlx

            CALL COMPFP_TREE_PBC(p,penglocal,virlocal,flocal,
     x            x,y,z,q,arrdim,nxboxl,nyboxl,nzboxl,idnode)
     
            f=f+flocal

            peng=peng+penglocal
            vir=vir+virlocal

         END DO

         nz=0; ny=0; nx=0
         nzboxl=REAL(nz,KIND=r8)*boxlz
         nyboxl=REAL(ny,KIND=r8)*boxly
         nxboxl=REAL(nx,KIND=r8)*boxlx

         tempx=x(i)
         tempq=q(i)

         q(i)=0.0_r8
         x(i)=tempx+1000.0_r8

         CALL COMPFP_TREE_PBC(p,penglocal,virlocal,flocal,
     x         x,y,z,q,arrdim,nxboxl,nyboxl,nzboxl,idnode)
     
         f=f+flocal

         peng=peng+0.5_r8*penglocal
         vir=vir+0.5_r8*virlocal

         x(i)=tempx
         q(i)=tempq

         ii=orderarr(i)

         tempq=r4pie0*tempq/epsq
         fx=-tempq*f(1)
         fy=-tempq*f(2)
         fz=-tempq*f(3)

         fxx(ii)=fx !-tempq*f(1)
         fyy(ii)=fy !-tempq*f(2)
         fzz(ii)=fz !-tempq*f(3)
         tpeng=tpeng+tempq*peng
         tvir=tvir-tempq*vir

         fff(1)=fff(1)+fx
         fff(2)=fff(2)+fy
         fff(3)=fff(3)+fz

!     complete stress tensor
        
         stress(1)=stress(1)-tempq*strs1
         stress(2)=stress(2)-tempq*strs2
         stress(3)=stress(3)-tempq*strs3
         stress(4)=stress(4)-tempq*strs2
         stress(5)=stress(5)-tempq*strs4
         stress(6)=stress(6)-tempq*strs5
         stress(7)=stress(7)-tempq*strs3
         stress(8)=stress(8)-tempq*strs5
         stress(9)=stress(9)-tempq*strs6

      END DO

c     remove COM drift arising from treecode approximations

      if(mxnode.gt.1)call gdsum(fff,3,buffer)

     
      fff(1)=fff(1)/REAL(arrdim,KIND=r8)
      fff(2)=fff(2)/REAL(arrdim,KIND=r8)
      fff(3)=fff(3)/REAL(arrdim,KIND=r8)

      do i=atm0,atm1

        fxx(i)=fxx(i)-fff(1)
        fyy(i)=fyy(i)-fff(2)
        fzz(i)=fzz(i)-fff(3)

      enddo


      RETURN
      END SUBROUTINE TREE_COMPFP_PBC
!!!!!!!!!!!!!!
      RECURSIVE SUBROUTINE COMPFP_TREE_PBC(p,peng,vir,f,x,y,z,q,arrdim,
     x                                      nxboxl,nyboxl,nzboxl,idnode)
     
      IMPLICIT NONE

      INTEGER,INTENT(IN) :: arrdim,idnode
      TYPE(tnode),POINTER :: p      
      REAL(KIND=r8),INTENT(INOUT) :: peng,vir
      REAL(KIND=r8),DIMENSION(3),INTENT(INOUT) :: f
      REAL(KIND=r8),DIMENSION(arrdim),INTENT(IN) :: x,y,z,q
      REAL(KIND=r8),INTENT(IN) :: nxboxl,nyboxl,nzboxl

! local variables

      REAL(KIND=r8) :: tx,ty,tz,distsq,t1,t2,penglocal,virlocal
      REAL(KIND=r8) :: cfj,cfk,mtmp,flocal(3),fflocal(3)
      REAL(KIND=r8), DIMENSION(1:trcsize) ::bvec,bdvecx,bdvecy,bdvecz 
      INTEGER :: i,j,k,j1,k1,err

! determine DISTSQ for MAC test

      tx=tarpos(1)-p%x_mid-nxboxl
      ty=tarpos(2)-p%y_mid-nyboxl
      tz=tarpos(3)-p%z_mid-nzboxl
      distsq=tx*tx+ty*ty+tz*tz

! intialize potential energy and force 

      vir=0.0_r8
      peng=0.0_r8
      f=0.0_r8

! If MAC is accepted and there is more than 1 particle in the 
! box use the expansion for the approximation.

      IF ((p%sqradius .LT. distsq*thetasq) .AND.
     x   (p%sqradius .NE. 0.0_r8)) THEN
         penglocal=0.0_r8; virlocal=0.0_r8   
         flocal=0.0_r8
         
c      CALL COMPFP_DIRECT_PBC(penglocal,flocal,p%ibeg,p%iend,
c     x                      x,y,z,q,arrdim,nxboxl,nyboxl,nzboxl)
c      WRITE(500,*)tarpos
c      WRITE(500,*)penglocal,flocal
c      WRITE(500,*)
c         penglocal=0.0_r8; virlocal=0.0_r8
c         flocal=0.0_r8
c         fflocal=0.0_r8

         IF (taylortree) THEN  ! Taylor treecode
            CALL COMP_TCOEFF_PBC(tx,ty,tz,distsq)
            DO k=0,torder
               cfk=cf(k); k1=k+1
               DO j=0,torder-k
                  cfj=cf(j); j1=j+1
                  DO i=0,torder-k-j
                     mtmp     =p%ms(i,j,k)
                     flocal(1)=flocal(1)+cf(i)*b1(i+1,j,k)*mtmp
                     flocal(2)=flocal(2)+cfj  *b1(i,j1,k) *mtmp
                     flocal(3)=flocal(3)+cfk  *b1(i,j,k1) *mtmp
                     peng     =peng+b1(i,j,k)*mtmp
                  END DO
               END DO
            END DO
c            WRITE(500,*)tarpos
c            WRITE(500,*)peng,flocal
c            WRITE(500,*)
c            STOP
         ELSE IF (tricubictree) THEN   ! Tricubic treecode
            bvec=0.0_r8; bdvecx=0.0_r8; bdvecy=0.0_r8; bdvecz=0.0_r8
            IF(marsden)THEN
               CALL kerEVAL_M(p,bvec,bdvecx,bdvecy,bdvecz,
     x                      nxboxl,nyboxl,nzboxl)
               peng     =peng+DOT_PRODUCT(bvec,p%trcms)
               flocal(1)=flocal(1)-DOT_PRODUCT(bvec,p%trcmsx)
               flocal(2)=flocal(2)-DOT_PRODUCT(bvec,p%trcmsy)
               flocal(3)=flocal(3)-DOT_PRODUCT(bvec,p%trcmsz)

               fflocal(1)=fflocal(1)-DOT_PRODUCT(bdvecx,p%trcms)
               fflocal(2)=fflocal(2)-DOT_PRODUCT(bdvecy,p%trcms)
               fflocal(3)=fflocal(3)-DOT_PRODUCT(bdvecz,p%trcms)

c               WRITE(500,*)tarpos
c               WRITE(500,*)peng,flocal
c               WRITE(500,*)peng,fflocal
c               WRITE(500,*)
c               STOP
             ELSE
                CALL kerEVAL_B2(p,bvec,bdvecx,bdvecy,bdvecz,
     x                      nxboxl,nyboxl,nzboxl)
               peng     =peng+DOT_PRODUCT(bvec,p%trcms)
               flocal(1)=flocal(1)-DOT_PRODUCT(bvec,p%trcmsx)
               flocal(2)=flocal(2)-DOT_PRODUCT(bvec,p%trcmsy)
               flocal(3)=flocal(3)-DOT_PRODUCT(bvec,p%trcmsz)

               fflocal(1)=fflocal(1)-DOT_PRODUCT(bdvecx,p%trcms)
               fflocal(2)=fflocal(2)-DOT_PRODUCT(bdvecy,p%trcms)
               fflocal(3)=fflocal(3)-DOT_PRODUCT(bdvecz,p%trcms)

c               WRITE(500,*)tarpos
c               WRITE(500,*)peng,flocal
c               WRITE(500,*)peng,fflocal
c               WRITE(500,*)
c               STOP
             END IF
         ELSE
           IF(idnode .eq. 0)WRITE(nrite,*)"TREECODE UNSPECIFIED"
           STOP
         END IF

         f=f+flocal
         vir=vir+flocal(1)*tx+flocal(2)*ty+flocal(3)*tz

!     calculate stress tensor
                
         strs1=strs1+tx*flocal(1)
         strs2=strs2+tx*flocal(2)
         strs3=strs3+tx*flocal(3)
         strs4=strs4+ty*flocal(2)
         strs5=strs5+ty*flocal(3)
         strs6=strs6+tz*flocal(3)

      ELSE

! If MAC fails check to see if the are children. If not, perform direct 
! calculation.  If there are children, call routine recursively for each.
!
      IF (p%num_children .EQ. 0) THEN
         CALL COMPFP_DIRECT_PBC(penglocal,flocal,p%ibeg,p%iend,
     x                      x,y,z,q,arrdim,nxboxl,nyboxl,nzboxl)
         vir=penglocal
         peng=penglocal
         f=flocal
      ELSE
         DO i=1,p%num_children
          CALL COMPFP_TREE_PBC(p%child(i)%p_to_tnode,penglocal,virlocal,
     x                flocal,x,y,z,q,arrdim,nxboxl,nyboxl,nzboxl,idnode)
            vir=vir+virlocal
            peng=peng+penglocal
            f=f+flocal
         END DO  
      END IF 
      END IF

      RETURN
      END SUBROUTINE COMPFP_TREE_PBC

!!!!!!!!!!!!!!
      SUBROUTINE kerEVAL_B2(p,bvec,bdvecx,bdvecy,bdvecz,
     x                   nxboxl,nyboxl,nzboxl)

      IMPLICIT NONE

! kerEVAL_B2 evaluates the function and its derivatives at the sixteen points, i.e.
! The eight corners and the eight internal points

      TYPE(tnode),POINTER :: p
      REAL(KIND=r8),DIMENSION(1:trcsize),INTENT(INOUT) :: bvec,bdvecx
      REAL(KIND=r8),DIMENSION(1:trcsize),INTENT(INOUT) :: bdvecy,bdvecz
      REAL(KIND=r8),INTENT(IN) :: nxboxl,nyboxl,nzboxl

! local variables

      INTEGER :: i,i16,i32,i48
      REAL(KIND=r8) :: tx,ty,tz,rrr,rr2,rr3,rr5

      DO i=1,16

        i16=i+16; i32=i+32; i48=i+48

        tx=tarpos(1)-p%ccp(i,1)-nxboxl
        ty=tarpos(2)-p%ccp(i,2)-nyboxl
        tz=tarpos(3)-p%ccp(i,3)-nzboxl

        rr2 = 1.0_r8/(tx*tx + ty*ty + tz*tz)
        rrr = SQRT(rr2)
        rr3 = rrr*rr2
        rr5 = rr2*rr3

        bvec(i)=rrr
        bvec(i16)=tx*rr3*p%dx
        bvec(i32)=ty*rr3*p%dy
        bvec(i48)=tz*rr3*p%dz

        bdvecx(i16)=(3.0_r8*tx*tx*rr5-rr3)*p%dx*p%dx
        bdvecx(i32)=3.0_r8*tx*ty*rr5*p%dx*p%dy
        bdvecx(i48)=3.0_r8*tx*tz*rr5*p%dx*p%dz

        bdvecy(i32)=(3.0_r8*ty*ty*rr5-rr3)*p%dy*p%dy
        bdvecy(i48)=3.0_r8*ty*tz*rr5*p%dy*p%dz

        bdvecz(i48)=(3.0_r8*tz*tz*rr5-rr3)*p%dz*p%dz

      END DO

      bdvecx(1:16) =bvec(17:32)
      
      bdvecy(1:16) =bvec(33:48)
      bdvecy(17:32)=bdvecx(33:48)

      bdvecz(1:16) =bvec(49:64)
      bdvecz(17:32)=bdvecx(49:64)
      bdvecz(33:48)=bdvecy(49:64)
 
      RETURN

      END SUBROUTINE kerEVAL_B2
!!!!!!!!!!!!!!
      SUBROUTINE kerEVAL_M(p,bvec,bdvecx,bdvecy,bdvecz,
     x                   nxboxl,nyboxl,nzboxl)

      IMPLICIT NONE

! kerEVAL_M evaluates the function and its derivatives at the eight corners

      TYPE(tnode),POINTER :: p
      REAL(KIND=r8),DIMENSION(1:trcsize),INTENT(INOUT) :: bvec,bdvecx
      REAL(KIND=r8),DIMENSION(1:trcsize),INTENT(INOUT) :: bdvecy,bdvecz
      REAL(KIND=r8),INTENT(IN) :: nxboxl,nyboxl,nzboxl

! local variables

      INTEGER :: i,i8,i16,i24,i32,i40,i48,i56
      REAL(KIND=r8) :: tx,ty,tz,rrr,rr2,rr3,rr5,rr7
      
      DO i=1,8
 
        i8=i+8; i16=i+16; i24=i+24; i32=i+32; 
        i40=i+40; i48=i+48; i56=i+56

        tx=tarpos(1)-p%ccp(i,1)-nxboxl
        ty=tarpos(2)-p%ccp(i,2)-nyboxl
        tz=tarpos(3)-p%ccp(i,3)-nzboxl

        rr2 = 1.0_r8/(tx*tx + ty*ty + tz*tz)
        rrr = SQRT(rr2)
        rr3 = rrr*rr2
        rr5 = rr2*rr3

        bvec(i)=rrr
        bvec(i8) =tx*rr3*p%dx
        bvec(i16)=ty*rr3*p%dy
        bvec(i24)=tz*rr3*p%dz
        bvec(i32)=3.0_r8*ty*bvec(i8)*rr2*p%dy
        bvec(i40)=3.0_r8*tz*bvec(i8)*rr2*p%dz
        bvec(i48)=3.0_r8*tz*bvec(i16)*rr2**p%dz
        bvec(i56)=5.0_r8*tz*bvec(i32)*rr2*p%dz

c        bdvecx(i16)=(tx*tx*rr5-rr3)*p%dx*p%dx
c        bdvecx(i32)=tx*ty*rr5*p%dx*p%dy
c        bdvecx(i48)=tx*tz*rr5*p%dx*p%dz

c        bdvecy(i32)=(ty*ty*rr5-rr3)*p%dy*p%dy
c        bdvecy(i48)=ty*tz*rr5*p%dy*p%dz

c        bdvecz(i48)=(tz*tz*rr5-rr3)*p%dz*p%dz

      END DO

c      bdvecx(1:16) =bvec(17:32)

c      bdvecy(1:16) =bvec(33:48)
c      bdvecy(17:32)=bdvecx(33:48)

c      bdvecz(1:16) =bvec(49:64)
c      bdvecz(17:32)=bdvecx(49:64)
c      bdvecz(33:48)=bdvecy(49:64)

      RETURN

      END SUBROUTINE kerEVAL_M


!!!!!!!!!!!!!!
      SUBROUTINE COMPFP_DIRECT_PBC(peng,f,ibeg,iend,x,y,z,q,arrdim,
     x                             nxboxl,nyboxl,nzboxl)
     
      IMPLICIT NONE

! COMPF_DIRECT directly computes the force on the current target
! particle determined by the global variable TARPOS.

      INTEGER,INTENT(IN) :: ibeg,iend,arrdim
      REAL(KIND=r8),DIMENSION(3),INTENT(OUT) :: f
      REAL(KIND=r8),INTENT(OUT) :: peng
      REAL(KIND=r8),DIMENSION(arrdim),INTENT(IN) :: x,y,z,q
      REAL(KIND=r8),INTENT(IN) :: nxboxl,nyboxl,nzboxl

! local variables

      INTEGER :: i
      REAL(KIND=r8) :: dist,t1,tx,ty,tz,qi,fx,fy,fz

      peng=0.0_r8
      f=0.0_r8
      DO i=ibeg,iend
         tx=x(i)+nxboxl-tarpos(1)
         ty=y(i)+nyboxl-tarpos(2)
         tz=z(i)+nzboxl-tarpos(3)
         qi=q(i)
         t1=1.0_r8/(tx*tx+ty*ty+tz*tz)

         dist=SQRT(t1)
         peng=peng+qi*dist

         t1=qi*dist*t1
         fx=tx*t1 
         fy=ty*t1
         fz=tz*t1

         f(1)=f(1)+fx
         f(2)=f(2)+fy
         f(3)=f(3)+fz

!     calculate stress tensor

         strs1=strs1+tx*fx
         strs2=strs2+tx*fy
         strs3=strs3+tx*fz
         strs4=strs4+ty*fy
         strs5=strs5+ty*fz
         strs6=strs6+tz*fz

      END DO  

      RETURN
      END SUBROUTINE COMPFP_DIRECT_PBC

!!!!!!!!
      SUBROUTINE COMP_TCOEFF_PBC(tx,ty,tz,distsq)
      IMPLICIT NONE

! COMP_TCOEFF computes the Taylor coefficients of the potential
! using a recurrence formula.  The center of the expansion is the
! midpoint of the node P.  TARPOS and TORDERLIM are globally defined.  
!
      REAL(KIND=r8),INTENT(IN) :: tx,ty,tz,distsq
      
! local variables

      REAL(KIND=r8) :: fac,sqfac,t1
      INTEGER :: k1,k2,k3,n,k21,k31,k22,k32

! setup variables

      fac=1.0_r8/(distsq)
      sqfac=SQRT(fac)

      b1=0.0_r8

! k1=0, k2=0, k3=0

      b1(0,0,0)=sqfac

!      k2=0, k3=0

      do k1=1,torderlim
         b1(k1,0,0)=fac*(cf1(k1)*tx*b1(k1-1,0,0)+ cf2(k1)*b1(k1-2,0,0))
      end do

!      k3=0
      do k2=1,torderlim
         k21=k2-1; k22=k2-2
         do k1=0,torderlim-k2
            n=k1+k2
            b1(k1,k2,0)=fac*(cf1(n)*(tx*b1(k1-1,k2,0)+ty*b1(k1,k21,0))+
     x                        cf2(n)*(b1(k1-2,k2,0)+b1(k1,k22,0)) )
     
         end do
      end do

      do k3=1,torderlim
         k31=k3-1; k32=k3-2
         do k2=0,torderlim-k3
            k21=k2-1; k22=k2-2
            do k1=0,torderlim-k3-k2
            n=k1+k2+k3
            b1(k1,k2,k3)=fac*(cf1(n)*(tx*b1(k1-1,k2,k3)+ 
     x                          ty*b1(k1,k21,k3)+tz*b1(k1,k2,k31))+
     x          cf2(n)*(b1(k1-2,k2,k3)+b1(k1,k22,k3)+b1(k1,k2,k32)))
            end do
         end do
       end do

      RETURN
      END SUBROUTINE COMP_TCOEFF_PBC    
!!!!!!!!!!!!!!!
      SUBROUTINE COMP_MS(p,x,y,z,q,arrdim)
      IMPLICIT NONE

! COMP_MS computes the moments for node P needed in the Taylor approximation
!
      INTEGER,INTENT(IN) :: arrdim
      TYPE(tnode),POINTER :: p 
      REAL(KIND=r8),DIMENSION(arrdim),INTENT(IN) :: x,y,z,q

! local variables

      INTEGER :: i,k1,k2,k3
      REAL(KIND=r8) :: dx,dy,dz,tx,ty,tz,qi
     
      p%ms=0.0_r8
      DO i=p%ibeg,p%iend
         dx=x(i)-p%x_mid
         dy=y(i)-p%y_mid
         dz=z(i)-p%z_mid
         qi=q(i)
         tz=1.0_r8
         DO k3=0,torder
            ty=tz
            DO k2=0,torder-k3
               tx=ty 
               DO k1=0,torder-k3-k2
                  p%ms(k1,k2,k3)=p%ms(k1,k2,k3)+qi*tx
                  tx=tx*dx
              END DO
              ty=ty*dy
            END DO
            tz=tz*dz
         END DO
      END DO
         
      RETURN
      END SUBROUTINE COMP_MS
!!!!!!!!!!!!!!!!!!!!!!!!!!
   
      SUBROUTINE TRICUBIC_MOMENTS(p,x,y,z,q,arrdim)
      IMPLICIT NONE

! TRICUBIC_MOMENTS computes the moments for node P needed in the Tricubic approximation
!
      INTEGER,INTENT(IN) :: arrdim
      TYPE(tnode),POINTER :: p
      REAL(KIND=r8),DIMENSION(arrdim),INTENT(IN) :: x,y,z,q

! localc variables

      INTEGER       :: id,i,j,k,i1,i4j,kk
      REAL(KIND=r8) :: tx,ty,tz,rtx,rty,rtz,qi,sij,s
      REAL(KIND=r8) :: mm(1:trcsize)
      LOGICAL       :: txzero,tyzero,tzzero
      mm=0.0_r8; 
      p%trcms=0.0_r8; p%trcmsx=0.0_r8; p%trcmsy=0.0_r8; p%trcmsz=0.0_r8;
      DO  id=p%ibeg,p%iend

          tx=(x(id)-p%x_min)*p%rdx
          ty=(y(id)-p%y_min)*p%rdy
          tz=(z(id)-p%z_min)*p%rdz

          txzero=(abs(tx) .LT. 1.0e-8)
          tyzero=(abs(ty) .LT. 1.0e-8)
          tzzero=(abs(tz) .LT. 1.0e-8)
          rtx=Merge(0.0_r8,1.0_r8/tx,txzero)
          rty=Merge(0.0_r8,1.0_r8/ty,tyzero)
          rtz=Merge(0.0_r8,1.0_r8/tz,tzzero)

          qi=q(id)
          DO i=0,trcord
   
             i1  = i+1
             sij = qi
             DO j=0,trcord
                
                i4j = i1 + 4 * j
                s   = sij
                DO k=0,trcord
   
                   kk = i4j + 16 * k
                   
                   p%trcms(kk) =p%trcms(kk) + s
                   p%trcmsx(kk)=p%trcmsx(kk)+REAL(i,kind=r8)*s*rtx*p%rdx
                   p%trcmsy(kk)=p%trcmsy(kk)+REAL(j,kind=r8)*s*rty*p%rdy
                   p%trcmsz(kk)=p%trcmsz(kk)+REAL(k,kind=r8)*s*rtz*p%rdz

                   s = s*tz
                END DO

                sij = sij*ty
             END DO

             qi = qi*tx

          END DO

       END DO     

! Multiplication by B-inverse transpose

       IF(marsden)THEN
          p%trcms = MATMUL(BinvTM,p%trcms)
          p%trcmsx= MATMUL(BinvTM,p%trcmsx)
          p%trcmsy= MATMUL(BinvTM,p%trcmsy)
          p%trcmsz= MATMUL(BinvTM,p%trcmsz)
       ELSE
          p%trcms = MATMUL(BinvT,p%trcms)
          p%trcmsx= MATMUL(BinvT,p%trcmsx)
          p%trcmsy= MATMUL(BinvT,p%trcmsy)
          p%trcmsz= MATMUL(BinvT,p%trcmsz)
       END IF
       RETURN
       END SUBROUTINE TRICUBIC_MOMENTS
!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE CLEANUP(p)
      IMPLICIT NONE

! CLEANUP deallocates allocated global variables and then
! calls recursive routine REMOVE_NODE to delete the tree.
!
      TYPE(tnode),POINTER :: p      

! local variables
  
      INTEGER :: err

      CALL REMOVE_NODE(p)
      DEALLOCATE(p, STAT=err)
      IF (err .NE. 0) THEN
         WRITE(6,*) 'Error deallocating root node! '
         STOP
      END IF 
      NULLIFY(p)         

      RETURN
      END SUBROUTINE CLEANUP
!!!!!!!!!!!
      RECURSIVE SUBROUTINE REMOVE_NODE(p)
      IMPLICIT NONE

! REMOVE_NODE recursively removes each node from the
! tree and deallocates its memory for MS array if it
! exits.
!
      TYPE(tnode),POINTER :: p 

! local variables

      INTEGER :: i,err

      IF (p%exist_ms .EQ. 1) THEN
         DEALLOCATE(p%ms,p%trcms,p%trcmsx,p%trcmsy,p%trcmsz,STAT=err)
         IF (err .NE. 0) THEN
            WRITE(6,*) 'Error deallocating node MS! '
            STOP
         END IF               
      END IF

      IF (p%num_children .GT. 0) THEN
          DO i=1,p%num_children
            CALL REMOVE_NODE(p%child(i)%p_to_tnode)
            DEALLOCATE(p%child(i)%p_to_tnode,STAT=err)
            IF (err .NE. 0) THEN
               WRITE(6,*) 'Error deallocating node child! '
               STOP
            END IF                           
          END DO
      END IF 

      RETURN                
      END SUBROUTINE REMOVE_NODE      

      END MODULE tree_module
!!!!!!!!!!!!!!!
      SUBROUTINE TREECODE(imcon,tpeng,tvir,epsq,idnode,mxnode,nstep)

      use config_module
      use setup_module
      use tree_module
      use utility_module

      IMPLICIT NONE

      INTEGER,INTENT(IN) :: imcon,idnode,mxnode,nstep
      REAL(KIND=r8),INTENT(IN) :: epsq
      REAL(KIND=r8),INTENT(OUT) :: tpeng,tvir

! local variables

      TYPE(tnode),POINTER :: troot
      INTEGER :: i,level,err,shrink

      REAL(KIND=r8)      :: boxl(3) 


      shrink=1

      txx=xxx
      tyy=yyy
      tzz=zzz
      tchge=chge

      call images(imcon,0,1,mxatms,cell,txx,tyy,tzz)

! find bounds of Cartesian box enclosing the particles

      xyzminmax(1)=MINVAL(txx(1:mxatms))
      xyzminmax(2)=MAXVAL(txx(1:mxatms))
      xyzminmax(3)=MINVAL(tyy(1:mxatms))
      xyzminmax(4)=MAXVAL(tyy(1:mxatms))
      xyzminmax(5)=MINVAL(tzz(1:mxatms))
      xyzminmax(6)=MAXVAL(tzz(1:mxatms))

! nullify pointer to root of tree (TROOT) and create tree

      NULLIFY(troot)  

! set global variables to track tree levels during construction

      level=0
      minlevel=50000
      maxlevel=0

      boxl(1)=cell(1); boxl(2)=cell(5); boxl(3)=cell(9)

      DO i=1,mxatms
         orderarr(i)=i
      END DO


      CALL CREATE_TREE(troot,1,mxatms,txx,tyy,tzz,tchge,shrink,
     x                xyzminmax,level,mxatms)


      CALL TREE_COMPFP_PBC(troot,txx,tyy,tzz,tchge,tpeng,tvir,boxl,
     x                     epsq,mxatms,idnode,mxnode)


! Call CLEANUP to deallocate global variables and tree structure.

      CALL CLEANUP(troot)

      END SUBROUTINE TREECODE
!!!!!!!!!!!!!
      SUBROUTINE PARTITION(a,b,c,q,indarr,ibeg,iend,val,midind,arrdim)
      IMPLICIT NONE
!
! PARTITION determines the index MIDIND, after partitioning
! in place the  arrays A,B,C and Q,  such that 
! A(IBEG:MIDIND) <= VAL and  A(MIDIND+1:IEND) > VAL. 
! If on entry IBEG > IEND or  A(IBEG:IEND) > VAL then MIDIND
! is returned as IBEG-1. 
! 
      INTEGER,PARAMETER :: r8=SELECTED_REAL_KIND(12)
      INTEGER, INTENT(IN) :: arrdim,ibeg,iend
      REAL(KIND=r8),DIMENSION(arrdim),INTENT(INOUT) :: a,b,c,q
      INTEGER,DIMENSION(arrdim),INTENT(INOUT) :: indarr   
      INTEGER, INTENT(INOUT) :: midind   
      REAL(KIND=r8) val

! local variables

      REAL(KIND=r8) ta,tb,tc,tq
      INTEGER lower,upper,tind

      IF (ibeg .LT. iend) THEN

! temporarily store IBEG entries and set A(IBEG)=VAL for 
! the partitoning algorithm.  

         ta=a(ibeg)
         tb=b(ibeg)
         tc=c(ibeg)
         tq=q(ibeg)
         tind=indarr(ibeg)
         a(ibeg)=val 
         upper=ibeg
         lower=iend

         DO WHILE (upper .NE. lower)
            DO WHILE ((upper .LT. lower) .AND. (val .LT. a(lower)))
                  lower=lower-1
            END DO
            IF (upper .NE. lower) THEN
               a(upper)=a(lower)
               b(upper)=b(lower)
               c(upper)=c(lower)
               q(upper)=q(lower)
               indarr(upper)=indarr(lower)
            END IF
            DO WHILE ((upper .LT. lower) .AND. (val .GE. a(upper)))
                  upper=upper+1
            END DO
            IF (upper .NE. lower) THEN
               a(lower)=a(upper)
               b(lower)=b(upper)
               c(lower)=c(upper)
               q(lower)=q(upper)
               indarr(lower)=indarr(upper)
            END IF
         END DO
         midind=upper

! replace TA in position UPPER and change MIDIND if TA > VAL 

         IF (ta .GT. val) THEN
            midind=upper-1
         END IF
         a(upper)=ta
         b(upper)=tb
         c(upper)=tc
         q(upper)=tq
         indarr(upper)=tind

      ELSEIF (ibeg .EQ. iend) THEN
         IF (a(ibeg) .LE. val) THEN
            midind=ibeg
         ELSE
            midind=ibeg-1
         END IF
      ELSE
         midind=ibeg-1
      END IF

      RETURN
      END SUBROUTINE PARTITION
!!!!!!!!!!!!!!!!!!
      subroutine tree_excl(iatm,ik,engcpe,vircpe,epsq)

!***********************************************************************
!     
!     dl_poly subroutine for calculating exclusion corrections to
!     the treecode method
!     
!     parallel replicated data version
!     
!     
!***********************************************************************
      
      use setup_module
      use exclude_module

      implicit none

      integer iatm,ik,m,jatm
      real(8) engcpe,vircpe,epsq,fx,fy,fz,strs1,strs2,strs3
      real(8) strs5,strs6,strs9,chgea,chgprd,rrr,rsq,coul,fcoul

      
!     initialise stress tensor accumulators

      strs1 = 0.d0
      strs2 = 0.d0
      strs3 = 0.d0
      strs5 = 0.d0
      strs6 = 0.d0
      strs9 = 0.d0

!     initialise potential energy and virial
      
      engcpe=0.d0
      vircpe=0.d0     

!     start of primary loop for forces evaluation
      
      chgea=chge(iatm)/epsq*r4pie0
      
      if(abs(chgea).gt.1.d-10)then
        
        do m=1,nexatm(ik)
!     atomic index and charge product
          
          jatm=lexatm(ik,m)
          chgprd=chgea*chge(jatm)
          
          if(abs(chgprd).gt.1.d-10)then

!     calculate interatomic distance
            
            rsq=xdf(m)**2+ydf(m)**2+zdf(m)**2
            rrr = sqrt(rsq)

!     calculate potential energy and virial
            
            coul = chgprd/rrr
            engcpe = engcpe - coul
            
!     calculate forces
            
            fcoul = coul/rsq
            fx = fcoul*xdf(m)
            fy = fcoul*ydf(m)
            fz = fcoul*zdf(m)
            
            fxx(iatm) = fxx(iatm) - fx
            fyy(iatm) = fyy(iatm) - fy
            fzz(iatm) = fzz(iatm) - fz

            fxx(jatm) = fxx(jatm) + fx
            fyy(jatm) = fyy(jatm) + fy
            fzz(jatm) = fzz(jatm) + fz

!     calculate stress tensor
            
            strs1 = strs1 - xdf(m)*fx
            strs2 = strs2 - xdf(m)*fy
            strs3 = strs3 - xdf(m)*fz
            strs5 = strs5 - ydf(m)*fy
            strs6 = strs6 - ydf(m)*fz
            strs9 = strs9 - zdf(m)*fz
            
          endif
          
        enddo

!     virial
        
        vircpe=vircpe-engcpe

!     complete stress tensor
        
        stress(1) = stress(1) + strs1
        stress(2) = stress(2) + strs2
        stress(3) = stress(3) + strs3
        stress(4) = stress(4) + strs2
        stress(5) = stress(5) + strs5
        stress(6) = stress(6) + strs6
        stress(7) = stress(7) + strs3
        stress(8) = stress(8) + strs6
        stress(9) = stress(9) + strs9

      endif
      
      return
      end subroutine tree_excl



