      module pair_module

c***********************************************************************
c     
c     dl_poly module for defining atom pair data
c     copyright - daresbury laboratory
c     author    - w. smith    mar 2004
c     
c***********************************************************************

      use setup_module
      use error_module

      implicit none

      integer, allocatable :: ilist(:),jlist(:)
      real(8), allocatable :: xdf(:),ydf(:),zdf(:)
      real(8), allocatable :: rsqdf(:)

      save ilist,jlist,xdf,ydf,zdf,rsqdf

      contains

      subroutine alloc_pair_arrays(idnode)

      implicit none

      integer, parameter :: nnn=6

      integer i,fail,idnode
      dimension fail(nnn)

      allocate (ilist(mxxdf),stat=fail(1))
      allocate (jlist(mxxdf),stat=fail(2))
      allocate (xdf(mxxdf),stat=fail(3))
      allocate (ydf(mxxdf),stat=fail(4))
      allocate (zdf(mxxdf),stat=fail(5))
      allocate (rsqdf(mxxdf),stat=fail(6))

      do i=1,nnn
        if(fail(i).gt.0)call error(idnode,1940)
      enddo

      end subroutine alloc_pair_arrays

      end module pair_module
