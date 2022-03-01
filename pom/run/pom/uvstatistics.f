!_______________________________________________________________________
      subroutine uvstatistics(iflg)

      implicit none
      include 'pom.h'

      integer iflg

      integer i,j,k

      integer init 
      data init /0/
      save init

      real riprint

      if(init.eq.0) then

      do k=1,kb
        do j=1,jm
          do i=1,im
            um1(i,j,k)=0.
            vm1(i,j,k)=0.
            um2(i,j,k)=0.
            vm2(i,j,k)=0.
            umx(i,j,k)=0.
            vmx(i,j,k)=0.
          end do
        end do
      end do
      init=1

      end if

      if(iflg.eq.0) then

      riprint=1.e0/real(iprint)
      do k=1,kbm1
        do j=1,jm
          do i=1,im
            um1(i,j,k)=um1(i,j,k)+u(i,j,k)*riprint
            vm1(i,j,k)=vm1(i,j,k)+v(i,j,k)*riprint
            um2(i,j,k)=um2(i,j,k)+u(i,j,k)*u(i,j,k)*riprint
            vm2(i,j,k)=vm2(i,j,k)+v(i,j,k)*v(i,j,k)*riprint
          end do
        end do
      end do

      do k=1,kbm1
        do j=1,jm
          do i=1,im
            if(abs(u(i,j,k)).gt.abs(umx(i,j,k)))then
              umx(i,j,k)=u(i,j,k)
            end if
            if(abs(v(i,j,k)).gt.abs(vmx(i,j,k)))then
              vmx(i,j,k)=v(i,j,k)
            end if
          end do
        end do
      end do

      else

      do k=1,kb
        do j=1,jm
          do i=1,im
            um1(i,j,k)=0.
            vm1(i,j,k)=0.
            um2(i,j,k)=0.
            vm2(i,j,k)=0.
            umx(i,j,k)=0.
            vmx(i,j,k)=0.
          end do
        end do
      end do

      end if

      return
      end
