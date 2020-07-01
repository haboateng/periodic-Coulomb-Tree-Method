     program kgfac_par
!=====================================================
! This program computes the Kirkwood G-factor  
! from a dl_poly run of water
! 
! This is the parallel version
! September 12, 2019
!=====================================================

     implicit none

     include "comms.inc"

     integer, parameter :: r8=selected_real_kind(12)

     real(kind=r8),allocatable,dimension(:) ::dpsum,dpsumAll
     real(kind=r8),allocatable,dimension(:) :: atmw,chge
     real(kind=r8),allocatable,dimension(:) :: xxx,yyy,zzz
     real(kind=r8) :: xi,yi,zi,xj,yj,zj,tstep,tmp,bgrid
     real(kind=r8) :: imux,imuy,imuz,jmux,jmuy,jmuz,boxlen,rgrid 
     real(kind=r8), parameter :: mulen = 0.59467146_r8
     integer :: nmols,natms,nsteps,nsteps1,nfiles,err,filenum
     integer :: atmnum,i,j,k,m,jj,i1,i2,i3,j1,j2,j3,ngrid

     character(len=20) :: filename,atmname
     character(len=80) :: header

! Parameters for parallel processing
     integer :: idnode, mxnode, ierr,fnum1,fnum2,iatm1,iatm2

!Begin MPI Communications

     call initcomms()
     call gsync()

!Find the id of processor

     call machine(idnode,mxnode)

     if (idnode==0) then
        read(*,*)nmols,nsteps,tstep,boxlen,ngrid,nfiles
     end if

     call MPI_BCAST(nmols,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
     call gsync()
     call MPI_BCAST(nsteps,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
     call gsync()
     call MPI_BCAST(tstep,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
     call gsync()
     call MPI_BCAST(boxlen,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
     call gsync()
     call MPI_BCAST(ngrid,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
     call gsync()
     call MPI_BCAST(nfiles,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
     call gsync()

     bgrid = boxlen/REAL(ngrid,KIND=r8)
     rgrid = REAL(ngrid,KIND=r8)/boxlen
     natms=nmols*3
     nsteps1=nsteps-1
     imols1 = (idnode*imols)/mxnode + 1
     imols2 = ((idnode+1)*imols)/mxnode

     allocate(dpsum(0:ngrid),dpsumAll(0:ngrid),dpcumsum(0:ngrid),atmw(natms),&
              chge(natms),xxx(natms),yyy(natms),zzz(natms),stat=err) 

     if(err.ne.0)then
       write(6,*)'error allocating compute arrays'
       stop
     end if

      dpsum = 0.0_r8; dpsumAll = 0.0_r8; dpcumsum = 0.0_r8
      dpsum(0)=1.0_r8

! read in data and sum up

     imols1 = (idnode*nmols)/mxnode + 1
     imols2 = ((idnode+1)*nmols)/mxnode
     do filenum=1,nfiles
        if(idnode.eq.0)then
           write(filename,'("HISTORY",I0)')filenum
           write(6,*)filenum, filename

           open(unit=82,form='unformatted',file=filename,status='old',action='read')
           read(82)
           read(82)
           read(82)
           read(82)atmw
           read(82)chge
        end if

        do k=0,nsteps1
           if(idnode.eq.0)then
              read(82)
              read(82)
              read(82)xxx
              read(82)yyy
              read(82)zzz
           end if


           call MPI_BCAST(xxx,natms,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
           call gsync()
           call MPI_BCAST(yyy,natms,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
           call gsync()
           call MPI_BCAST(zzz,natms,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
           call gsync()

! read in positions and velocities and find the center of mass of
! each molecule as well as the velocity of the center of mass

           do i=imols1,imols2

              i1=3*(i-1)+1; i2=i1+1; i3=i1+2;

! coordinate of oxygen atom of molecule i
              xi=xxx(i1); yi=yyy(i1); zi=zzz(i1);

! unit dipole vector of molecule i
              imux=(xi-0.5_r8*(xxx(i2)+xxx(i3)))/mulen;
              imuy=(yi-0.5_r8*(yyy(i2)+yyy(i3)))/mulen;
              imuz=(zi-0.5_r8*(zzz(i2)+zzz(i3)))/mulen;

              
              do j = 1,i-1

                 j1=3*(j-1)+1; j2=j1+1; j3=j1+2

! coordinate of oxygen atom of molecule j
                 xj=xxx(j1); yj=yyy(j1); zj=zzz(j1)

! unit dipole vector of molecule j
                 dx=xi-xj; dy=yi-yj; dx=zi-zj

                 rrr = sqrt(dx*dx + dy*dy + dz*dz)
                 kk  = CEILING(rrr*rgrid)

                 jmux=(xj-0.5_r8*(xxx(j2)+xxx(j3)))/mulen;
                 jmuy=(yj-0.5_r8*(yyy(j2)+yyy(j3)))/mulen;
                 jmuz=(zj-0.5_r8*(zzz(j2)+zzz(j3)))/mulen;

                 temp=imux*jmux+imuy+jmuy+imuz*jmuz
                 dpsum(kk)=dpsum(kk)+temp
              end do

              do j = i+1,nmols

                 j1=3*(j-1)+1; j2=j1+1; j3=j1+2

! coordinate of oxygen atom of molecule j
                 xj=xxx(j1); yj=yyy(j1); zj=zzz(j1)

! unit dipole vector of molecule j
                 dx=xi-xj; dy=yi-yj; dx=zi-zj

                 rrr = sqrt(dx*dx + dy*dy + dz*dz)
                 kk  = CEILING(rrr*rgrid)

                 jmux=(xj-0.5_r8*(xxx(j2)+xxx(j3)))/mulen;
                 jmuy=(yj-0.5_r8*(yyy(j2)+yyy(j3)))/mulen;
                 jmuz=(zj-0.5_r8*(zzz(j2)+zzz(j3)))/mulen;

                 temp=imux*jmux+imuy+jmuy+imuz*jmuz
                 dpsum(kk)=dpsum(kk)+temp
              end do

           end do
          
        end do

        if(idnode.eq.0) close(82)

     end do

! Sum up dpsum and send to master

     call MPI_REDUCE(dpsum,dpsumAll,ngrid,MPI_DOUBLE_PRECISION,&
           MPI_SUM,0,MPI_COMM_WORLD,ierr)

! find cumulative sum, average over number of molecules, steps and files
! and write out Kirkwood G-factor data 

     if(idnode.eq.0)then
        temp = real(nmols*nsteps*nfiles)

        open(unit=83,file='gfac.txt',status='unknown',action='write',&
             position='rewind')

        do i=0,ngrid
           write(83,13)bgrid*real(i,kind=r8),sum(dpsumAll(0:i))/temp 
        end do

        close(83)
     end if

 13   FORMAT(3(2X,E24.16))

     deallocate(dpsum,dpsumAll)
     deallocate(xxx,yyy,zzz)

     call exitcomms()

     end program kgfac_par

