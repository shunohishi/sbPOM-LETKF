! bounds_forcing_frs.f

!_______________________________________________________________________
      subroutine bcond_frs(idx)
! apply open boundary conditions by nesting with
! the flow relaxation scheme

      implicit none
      include 'pom.h'
      integer idx
      integer i,j,k
      real cfrs

      if(idx.eq.4) then

! temperature and salinity boundary conditions (using uf and vf,
! respectively)

        do i=2,4
          do k=1,kbm1
            do j=1,jm
!              ! east
              if(n_east.eq.-1 .and. fsm(im,j).eq.1) then
                cfrs=1.0-tanh(0.5*(i-1))
                uf(im-i+1,j,k)=fsm(im-i+1,j)*
     $            (cfrs*tref(im-i+1,j,k)+(1.0-cfrs)*uf(im-i+1,j,k))
                vf(im-i+1,j,k)=fsm(im-i+1,j)*
     $            (cfrs*sref(im-i+1,j,k)+(1.0-cfrs)*vf(im-i+1,j,k))
              end if
!              ! west
              if(n_west.eq.-1 .and. fsm(1,j).eq.1) then
                cfrs=1.0-tanh(0.5*(i-1))
                uf(i,j,k)=fsm(i,j)*
     $            (cfrs*tref(i,j,k)+(1.0-cfrs)*uf(i,j,k))
                vf(i,j,k)=fsm(i,j)*
     $            (cfrs*sref(i,j,k)+(1.0-cfrs)*vf(i,j,k))
              end if
            end do
          end do
        end do

        do j=2,4
          do k=1,kbm1
            do i=1,im
              ! north
              if(n_north.eq.-1 .and. fsm(i,jm).eq.1) then
                cfrs=1.0-tanh(0.5*(j-1))
                uf(i,jm-j+1,k)=fsm(i,jm-j+1)*
     $           (cfrs*tref(i,jm-j+1,k)+(1.0-cfrs)*uf(i,jm-j+1,k))
                vf(i,jm-j+1,k)=fsm(i,jm-j+1)*
     $           (cfrs*sref(i,jm-j+1,k)+(1.0-cfrs)*vf(i,jm-j+1,k))
              end if
              ! south
              if(n_south.eq.-1 .and. fsm(i,1).eq.1) then
                cfrs=1.0-tanh(0.5*(j-1))
                uf(i,j,k)=fsm(i,j)*
     $           (cfrs*tref(i,j,k)+(1.0-cfrs)*uf(i,j,k))
                vf(i,j,k)=fsm(i,j)*
     $           (cfrs*sref(i,j,k)+(1.0-cfrs)*vf(i,j,k))
              end if
            end do
          end do
        end do

!        write(6,*) 'in frs uf tref ', uf(2,jm/2,1), tref(2,jm/2,1)
!        write(6,*) 'in frs vf sref ', vf(2,jm/2,1), sref(2,jm/2,1)
      endif

      return

      end

