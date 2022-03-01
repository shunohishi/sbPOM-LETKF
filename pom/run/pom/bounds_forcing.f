! bounds_forcing.f

! spcify variable boundary conditions, atmospheric forcing, restoring
!_______________________________________________________________________
      subroutine bcond_orig(idx)
! apply open boundary conditions
! closed boundary conditions are automatically enabled through
! specification of the masks, dum, dvm and fsm, in which case the open
! boundary conditions, included below, will be overwritten
      implicit none
      include 'pom.h'
      integer idx
      integer i,j,k
      real ga,u1,wm
      real hmax

      if(idx.eq.1) then

! External (2-D) elevation boundary conditions
        do j=1,jm
          if(n_west.eq.-1) elf(1,j)=elf(2,j)
          if(n_east.eq.-1) elf(im,j)=elf(imm1,j)
        end do

        do i=1,im
          if(n_south.eq.-1) elf(i,1)=elf(i,2)
          if(n_north.eq.-1) elf(i,jm)=elf(i,jmm1)
        end do

        do j=1,jm
          do i=1,im
            elf(i,j)=elf(i,j)*fsm(i,j)
          end do
        end do

        return

      else if(idx.eq.2) then

! external (2-D) velocity boundary conditions
        do j=2,jmm1
          ! west
          if(n_west.eq.-1) then
            uaf(2,j)=uabw(j)-rfw*sqrt(grav/h(2,j))*(el(2,j)-elw(j))
            uaf(2,j)=ramp*uaf(2,j)
            uaf(1,j)=uaf(2,j)
            vaf(1,j)=0.e0
          end if

          ! east
          if(n_east.eq.-1) then
            uaf(im,j)=uabe(j)
     $                     +rfe*sqrt(grav/h(imm1,j))*(el(imm1,j)-ele(j))
            uaf(im,j)=ramp*uaf(im,j)
            vaf(im,j)=0.e0
          end if
        end do

        do i=2,imm1
          ! south
          if(n_south.eq.-1) then
            vaf(i,2)=vabs(i)-rfs*sqrt(grav/h(i,2))*(el(i,2)-els(i))
            vaf(i,2)=ramp*vaf(i,2)
            vaf(i,1)=vaf(i,2)
            uaf(i,1)=0.e0
          end if

          ! north
          if(n_north.eq.-1) then
            vaf(i,jm)=vabn(i)
     $                     +rfn*sqrt(grav/h(i,jmm1))*(el(i,jmm1)-eln(i))
            vaf(i,jm)=ramp*vaf(i,jm)
            uaf(i,jm)=0.e0
          end if
        end do

        do j=1,jm
          do i=1,im
            uaf(i,j)=uaf(i,j)*dum(i,j)
            vaf(i,j)=vaf(i,j)*dvm(i,j)
          end do
        end do

        return

      else if(idx.eq.3) then

! internal (3-D) velocity boundary conditions
! radiation conditions
! smoothing is used in the direction tangential to the boundaries
        hmax=4500.e0

        do k=1,kbm1
          do j=2,jmm1
            ! east
            if(n_east.eq.-1) then
              ga=sqrt(h(im,j)/hmax)
              uf(im,j,k)=ga*(.25e0*u(imm1,j-1,k)+.5e0*u(imm1,j,k)
     $                       +.25e0*u(imm1,j+1,k))
     $                    +(1.e0-ga)*(.25e0*u(im,j-1,k)+.5e0*u(im,j,k)
     $                      +.25e0*u(im,j+1,k))
              vf(im,j,k)=0.e0
            end if
            ! west
            if(n_west.eq.-1) then
              ga=sqrt(h(1,j)/hmax)
              uf(2,j,k)=ga*(.25e0*u(3,j-1,k)+.5e0*u(3,j,k)
     $                      +.25e0*u(3,j+1,k))
     $                   +(1.e0-ga)*(.25e0*u(2,j-1,k)+.5e0*u(2,j,k)
     $                     +.25e0*u(2,j+1,k))
              uf(1,j,k)=uf(2,j,k)
              vf(1,j,k)=0.e0
            end if
          end do
        end do

        do k=1,kbm1
          do i=2,imm1
            ! south
            if(n_south.eq.-1) then
              ga=sqrt(h(i,1)/hmax)
              vf(i,2,k)=ga*(.25e0*v(i-1,3,k)+.5e0*v(i,3,k)
     $                      +.25e0*v(i+1,3,k))
     $                   +(1.e0-ga)*(.25e0*v(i-1,2,k)+.5e0*v(i,2,k)
     $                     +.25e0*v(i+1,2,k))
              vf(i,1,k)=vf(i,2,k)
              uf(i,1,k)=0.e0
            end if
            ! north
            if(n_north.eq.-1) then
              ga=sqrt(h(i,jm)/hmax)
              vf(i,jm,k)=ga*(.25e0*v(i-1,jmm1,k)+.5e0*v(i,jmm1,k)
     $                       +.25e0*v(i+1,jmm1,k))
     $                    +(1.e0-ga)*(.25e0*v(i-1,jm,k)+.5e0*v(i,jm,k)
     $                      +.25e0*v(i+1,jm,k))
              uf(i,jm,k)=0.e0
            end if
          end do
        end do

        do k=1,kbm1
          do j=1,jm
            do i=1,im
              uf(i,j,k)=uf(i,j,k)*dum(i,j)
              vf(i,j,k)=vf(i,j,k)*dvm(i,j)
            end do
          end do
        end do

        return

      else if(idx.eq.4) then

! temperature and salinity boundary conditions (using uf and vf,
! respectively)
        do k=1,kbm1
          do j=1,jm
            ! east
            if(n_east.eq.-1) then
              u1=2.e0*u(im,j,k)*dti/(dx(im,j)+dx(imm1,j))
              if(u1.le.0.e0) then
                uf(im,j,k)=t(im,j,k)-u1*(tbe(j,k)-t(im,j,k))
                vf(im,j,k)=s(im,j,k)-u1*(sbe(j,k)-s(im,j,k))
              else
                uf(im,j,k)=t(im,j,k)-u1*(t(im,j,k)-t(imm1,j,k))
                vf(im,j,k)=s(im,j,k)-u1*(s(im,j,k)-s(imm1,j,k))
                if(k.ne.1.and.k.ne.kbm1) then
                  wm=.5e0*(w(imm1,j,k)+w(imm1,j,k+1))*dti
     $                /((zz(k-1)-zz(k+1))*dt(imm1,j))
                  uf(im,j,k)=uf(im,j,k)-wm*(t(imm1,j,k-1)-t(imm1,j,k+1))
                  vf(im,j,k)=vf(im,j,k)-wm*(s(imm1,j,k-1)-s(imm1,j,k+1))
                endif
              end if
            end if

            ! west
            if(n_west.eq.-1) then
              u1=2.e0*u(2,j,k)*dti/(dx(1,j)+dx(2,j))
              if(u1.ge.0.e0) then
                uf(1,j,k)=t(1,j,k)-u1*(t(1,j,k)-tbw(j,k))
                vf(1,j,k)=s(1,j,k)-u1*(s(1,j,k)-sbw(j,k))
              else
                uf(1,j,k)=t(1,j,k)-u1*(t(2,j,k)-t(1,j,k))
                vf(1,j,k)=s(1,j,k)-u1*(s(2,j,k)-s(1,j,k))
                if(k.ne.1.and.k.ne.kbm1) then
                  wm=.5e0*(w(2,j,k)+w(2,j,k+1))*dti
     $                /((zz(k-1)-zz(k+1))*dt(2,j))
                  uf(1,j,k)=uf(1,j,k)-wm*(t(2,j,k-1)-t(2,j,k+1))
                  vf(1,j,k)=vf(1,j,k)-wm*(s(2,j,k-1)-s(2,j,k+1))
                end if
              end if
            end if
          end do

          do i=1,im
            ! south
            if(n_south.eq.-1) then
              u1=2.e0*v(i,2,k)*dti/(dy(i,1)+dy(i,2))
              if(u1.ge.0.e0) then
                uf(i,1,k)=t(i,1,k)-u1*(t(i,1,k)-tbs(i,k))
                vf(i,1,k)=s(i,1,k)-u1*(s(i,1,k)-sbs(i,k))
              else
                uf(i,1,k)=t(i,1,k)-u1*(t(i,2,k)-t(i,1,k))
                vf(i,1,k)=s(i,1,k)-u1*(s(i,2,k)-s(i,1,k))
                if(k.ne.1.and.k.ne.kbm1) then
                  wm=.5e0*(w(i,2,k)+w(i,2,k+1))*dti
     $                /((zz(k-1)-zz(k+1))*dt(i,2))
                  uf(i,1,k)=uf(i,1,k)-wm*(t(i,2,k-1)-t(i,2,k+1))
                  vf(i,1,k)=vf(i,1,k)-wm*(s(i,2,k-1)-s(i,2,k+1))
                end if
              end if
            end if

            ! north
            if(n_north.eq.-1) then
              u1=2.e0*v(i,jm,k)*dti/(dy(i,jm)+dy(i,jmm1))
              if(u1.le.0.e0) then
                uf(i,jm,k)=t(i,jm,k)-u1*(tbn(i,k)-t(i,jm,k))
                vf(i,jm,k)=s(i,jm,k)-u1*(sbn(i,k)-s(i,jm,k))
              else
                uf(i,jm,k)=t(i,jm,k)-u1*(t(i,jm,k)-t(i,jmm1,k))
                vf(i,jm,k)=s(i,jm,k)-u1*(s(i,jm,k)-s(i,jmm1,k))
                if(k.ne.1.and.k.ne.kbm1) then
                  wm=.5e0*(w(i,jmm1,k)+w(i,jmm1,k+1))*dti
     $                /((zz(k-1)-zz(k+1))*dt(i,jmm1))
                  uf(i,jm,k)=uf(i,jm,k)-wm*(t(i,jmm1,k-1)-t(i,jmm1,k+1))
                  vf(i,jm,k)=vf(i,jm,k)-wm*(s(i,jmm1,k-1)-s(i,jmm1,k+1))
                end if
              end if
            end if
          end do
        end do

        do k=1,kbm1
          do j=1,jm
            do i=1,im
              uf(i,j,k)=uf(i,j,k)*fsm(i,j)
              vf(i,j,k)=vf(i,j,k)*fsm(i,j)
            end do
          end do
        end do

        return

      else if(idx.eq.5) then

! vertical velocity boundary conditions
        do k=1,kbm1
          do j=1,jm
            do i=1,im
              w(i,j,k)=w(i,j,k)*fsm(i,j)
            end do
          end do
        end do

        return

      else if(idx.eq.6) then

! q2 and q2l boundary conditions

        do k=1,kb
          do j=1,jm
            ! west
            if(n_west.eq.-1) then
              u1=2.e0*u(2,j,k)*dti/(dx(1,j)+dx(2,j))
              if(u1.ge.0.e0) then
                uf(1,j,k)=q2(1,j,k)-u1*(q2(1,j,k)-small)
                vf(1,j,k)=q2l(1,j,k)-u1*(q2l(1,j,k)-small)
              else
                uf(1,j,k)=q2(1,j,k)-u1*(q2(2,j,k)-q2(1,j,k))
                vf(1,j,k)=q2l(1,j,k)-u1*(q2l(2,j,k)-q2l(1,j,k))
              end if
            end if

            ! east
            if(n_east.eq.-1) then
              u1=2.e0*u(im,j,k)*dti/(dx(im,j)+dx(imm1,j))
              if(u1.le.0.e0) then
                uf(im,j,k)=q2(im,j,k)-u1*(small-q2(im,j,k))
                vf(im,j,k)=q2l(im,j,k)-u1*(small-q2l(im,j,k))
              else
                uf(im,j,k)=q2(im,j,k)-u1*(q2(im,j,k)-q2(imm1,j,k))
                vf(im,j,k)=q2l(im,j,k)-u1*(q2l(im,j,k)-q2l(imm1,j,k))
              end if
            end if
          end do

          do i=1,im
            ! south
            if(n_south.eq.-1) then
              u1=2.e0*v(i,2,k)*dti/(dy(i,1)+dy(i,2))
              if(u1.ge.0.e0) then
                uf(i,1,k)=q2(i,1,k)-u1*(q2(i,1,k)-small)
                vf(i,1,k)=q2l(i,1,k)-u1*(q2l(i,1,k)-small)
              else
                uf(i,1,k)=q2(i,1,k)-u1*(q2(i,2,k)-q2(i,1,k))
                vf(i,1,k)=q2l(i,1,k)-u1*(q2l(i,2,k)-q2l(i,1,k))
              end if
            end if

            ! north
            if(n_north.eq.-1) then
              u1=2.e0*v(i,jm,k)*dti/(dy(i,jm)+dy(i,jmm1))
              if(u1.le.0.e0) then
                uf(i,jm,k)=q2(i,jm,k)-u1*(small-q2(i,jm,k))
                vf(i,jm,k)=q2l(i,jm,k)-u1*(small-q2l(i,jm,k))
              else
                uf(i,jm,k)=q2(i,jm,k)-u1*(q2(i,jm,k)-q2(i,jmm1,k))
                vf(i,jm,k)=q2l(i,jm,k)-u1*(q2l(i,jm,k)-q2l(i,jmm1,k))
              end if
            end if
          end do
        end do

        do k=1,kb
          do j=1,jm
            do i=1,im
              uf(i,j,k)=uf(i,j,k)*fsm(i,j)+1.e-10
              vf(i,j,k)=vf(i,j,k)*fsm(i,j)+1.e-10
            end do
          end do
        end do

        return

      endif

      end

!_______________________________________________________________________
      subroutine bcond(idx)
! apply open boundary conditions by nesting
! closed boundary conditions are automatically enabled through
! specification of the masks, dum, dvm and fsm, in which case the open
! boundary conditions, included below, will be overwritten
      implicit none
      include 'pom.h'
      integer idx
      integer i,j,k
      real ga,u1,wm
      real hmax

      if(idx.eq.1) then

! External (2-D) elevation boundary conditions
        do j=1,jm
          if(n_west.eq.-1) elf(1,j)=elw(j)
          if(n_east.eq.-1) elf(im,j)=ele(j)
        end do

        do i=1,im
          if(n_south.eq.-1) elf(i,1)=els(i)
          if(n_north.eq.-1) elf(i,jm)=eln(i)
        end do

        do j=1,jm
          do i=1,im
            elf(i,j)=elf(i,j)*fsm(i,j)
          end do
        end do

        return

      else if(idx.eq.2) then

! external (2-D) velocity boundary conditions
        do j=2,jmm1
          ! west
          if(n_west.eq.-1) then
            uaf(2,j)=uabw(j)-rfw*sqrt(grav/h(2,j))*(el(2,j)-elw(j))
            uaf(2,j)=ramp*uaf(2,j)
            uaf(1,j)=uaf(2,j)
            vaf(1,j)=vabw(j)
          end if

          ! east
          if(n_east.eq.-1) then
            uaf(im,j)=uabe(j)
     $                     +rfe*sqrt(grav/h(imm1,j))*(el(imm1,j)-ele(j))
            uaf(im,j)=ramp*uaf(im,j)
            vaf(im,j)=vabe(j)
          end if
        end do

        do i=2,imm1
          ! south
          if(n_south.eq.-1) then
            vaf(i,2)=vabs(i)-rfs*sqrt(grav/h(i,2))*(el(i,2)-els(i))
            vaf(i,2)=ramp*vaf(i,2)
            vaf(i,1)=vaf(i,2)
            uaf(i,1)=uabs(i)
          end if

          ! north
          if(n_north.eq.-1) then
            vaf(i,jm)=vabn(i)
     $                     +rfn*sqrt(grav/h(i,jmm1))*(el(i,jmm1)-eln(i))
            vaf(i,jm)=ramp*vaf(i,jm)
            uaf(i,jm)=uabn(i)
          end if
        end do

        do j=1,jm
          do i=1,im
            uaf(i,j)=uaf(i,j)*dum(i,j)
            vaf(i,j)=vaf(i,j)*dvm(i,j)
          end do
        end do

        return

      else if(idx.eq.3) then

! internal (3-D) velocity boundary conditions
! radiation conditions
! smoothing is used in the direction tangential to the boundaries
        hmax=4500.e0

        do k=1,kbm1
          do j=2,jmm1
            ! east
            if(n_east.eq.-1) then
              uf(im,j,k)=ube(j,k)
              vf(im,j,k)=vbe(j,k)
            end if
            ! west
            if(n_west.eq.-1) then
              uf(2,j,k)=ubw(j,k)
              uf(1,j,k)=uf(2,j,k)
              vf(1,j,k)=vbw(j,k)
            end if
          end do
        end do

        do k=1,kbm1
          do i=2,imm1
            ! south
            if(n_south.eq.-1) then
              vf(i,2,k)=vbs(i,k)
              vf(i,1,k)=vf(i,2,k)
              uf(i,1,k)=ubs(i,k)
            end if
            ! north
            if(n_north.eq.-1) then
              vf(i,jm,k)=vbn(i,k)
              uf(i,jm,k)=ubn(i,k)
            end if
          end do
        end do

        do k=1,kbm1
          do j=1,jm
            do i=1,im
              uf(i,j,k)=uf(i,j,k)*dum(i,j)
              vf(i,j,k)=vf(i,j,k)*dvm(i,j)
            end do
          end do
        end do

        return

      else if(idx.eq.4) then

! temperature and salinity boundary conditions (using uf and vf,
! respectively)
        do k=1,kbm1
          do j=1,jm
            ! east
            if(n_east.eq.-1) then
              uf(im,j,k)=tbe(j,k) 
              vf(im,j,k)=sbe(j,k) 
            end if

            ! west
            if(n_west.eq.-1) then
              uf(1,j,k)=tbw(j,k)
              vf(1,j,k)=sbw(j,k)
            end if
          end do

          do i=1,im
            ! south
            if(n_south.eq.-1) then
              uf(i,1,k)=tbs(i,k)
              vf(i,1,k)=sbs(i,k)
            end if

            ! north
            if(n_north.eq.-1) then
              uf(i,jm,k)=tbn(i,k)
              vf(i,jm,k)=sbn(i,k)
            end if
          end do
        end do

        do k=1,kbm1
          do j=1,jm
            do i=1,im
              uf(i,j,k)=uf(i,j,k)*fsm(i,j)
              vf(i,j,k)=vf(i,j,k)*fsm(i,j)
            end do
          end do
        end do

        return

      else if(idx.eq.5) then

! vertical velocity boundary conditions
        do k=1,kbm1
          do j=1,jm
            do i=1,im
              w(i,j,k)=w(i,j,k)*fsm(i,j)
            end do
          end do
        end do

        return

      else if(idx.eq.6) then

! q2 and q2l boundary conditions

        do k=1,kb
          do j=1,jm
            ! west
            if(n_west.eq.-1) then
              u1=2.e0*u(2,j,k)*dti/(dx(1,j)+dx(2,j))
              if(u1.ge.0.e0) then
                uf(1,j,k)=q2(1,j,k)-u1*(q2(1,j,k)-small)
                vf(1,j,k)=q2l(1,j,k)-u1*(q2l(1,j,k)-small)
              else
                uf(1,j,k)=q2(1,j,k)-u1*(q2(2,j,k)-q2(1,j,k))
                vf(1,j,k)=q2l(1,j,k)-u1*(q2l(2,j,k)-q2l(1,j,k))
              end if
            end if

            ! east
            if(n_east.eq.-1) then
              u1=2.e0*u(im,j,k)*dti/(dx(im,j)+dx(imm1,j))
              if(u1.le.0.e0) then
                uf(im,j,k)=q2(im,j,k)-u1*(small-q2(im,j,k))
                vf(im,j,k)=q2l(im,j,k)-u1*(small-q2l(im,j,k))
              else
                uf(im,j,k)=q2(im,j,k)-u1*(q2(im,j,k)-q2(imm1,j,k))
                vf(im,j,k)=q2l(im,j,k)-u1*(q2l(im,j,k)-q2l(imm1,j,k))
              end if
            end if
          end do

          do i=1,im
            ! south
            if(n_south.eq.-1) then
              u1=2.e0*v(i,2,k)*dti/(dy(i,1)+dy(i,2))
              if(u1.ge.0.e0) then
                uf(i,1,k)=q2(i,1,k)-u1*(q2(i,1,k)-small)
                vf(i,1,k)=q2l(i,1,k)-u1*(q2l(i,1,k)-small)
              else
                uf(i,1,k)=q2(i,1,k)-u1*(q2(i,2,k)-q2(i,1,k))
                vf(i,1,k)=q2l(i,1,k)-u1*(q2l(i,2,k)-q2l(i,1,k))
              end if
            end if

            ! north
            if(n_north.eq.-1) then
              u1=2.e0*v(i,jm,k)*dti/(dy(i,jm)+dy(i,jmm1))
              if(u1.le.0.e0) then
                uf(i,jm,k)=q2(i,jm,k)-u1*(small-q2(i,jm,k))
                vf(i,jm,k)=q2l(i,jm,k)-u1*(small-q2l(i,jm,k))
              else
                uf(i,jm,k)=q2(i,jm,k)-u1*(q2(i,jm,k)-q2(i,jmm1,k))
                vf(i,jm,k)=q2l(i,jm,k)-u1*(q2l(i,jm,k)-q2l(i,jmm1,k))
              end if
            end if
          end do
        end do

        do k=1,kb
          do j=1,jm
            do i=1,im
              uf(i,j,k)=uf(i,j,k)*fsm(i,j)+1.e-10
              vf(i,j,k)=vf(i,j,k)*fsm(i,j)+1.e-10
            end do
          end do
        end do

        return

      endif

      end

!_______________________________________________________________________
      subroutine bcondorl(idx)
! this is an optional subroutine replacing  bcond and using Orlanski's
! scheme (J. Comp. Phys. 21, 251-269, 1976), specialized for the
! seamount problem
      implicit none
      include 'pom.h'
      integer idx
      real cl,denom
      real ar,eps
      integer i,j,k

      if(idx.eq.1) then
      
! external (2-D) elevation boundary conditions
        do  j=1,jm
          if(n_west.eq.-1) elf(1,j)=elf(2,j)
          if(n_east.eq.-1) elf(im,j)=elf(imm1,j)
        end do

        do j=1,jm
          do i=1,im
            elf(i,j)=elf(i,j)*fsm(i,j)
          end do
        end do

        return

      else if(idx.eq.2) then

! external (2-D) velocity  boundary conditions
        do j=2,jmm1
          ! east
          if(n_east.eq.-1) then
            denom=(uaf(im-1,j)+uab(im-1,j)-2.e0*ua(im-2,j))
            if(denom.eq.0.0e0)denom=0.01e0
            cl=(uab(im-1,j)-uaf(im-1,j))/denom
            if(cl.gt.1.e0) cl=1.e0
            if(cl.lt.0.e0) cl=0.e0
            uaf(im,j)=(uab(im,j)*(1.e0-cl)+2.e0*cl*ua(im-1,j))
     $                  /(1.e0+cl)
            vaf(im,j)=0.e0
          end if

          ! west
          if(n_west.eq.-1) then
            denom=(uaf(3,j)+uab(3,j)-2.e0*ua(4,j))
            if(denom.eq.0.0e0)denom=0.01e0
            cl=(uab(3,j)-uaf(3,j))/denom
            if(cl.gt.1.e0) cl=1.e0
            if(cl.lt.0.e0) cl=0.e0
            uaf(2,j)=(uab(2,j)*(1.e0-cl)+2.e0*cl*ua(3,j))
     $                 /(1.e0+cl)
            uaf(1,j)=uaf(2,j)
            vaf(1,j)=0.e0
          end if
        end do

        do j=1,jm
          do i=1,im
            uaf(i,j)=uaf(i,j)*dum(i,j)
            vaf(i,j)=vaf(i,j)*dvm(i,j)
          end do
        end do

        return

      else if(idx.eq.3) then

! internal (3-D) velocity boundary conditions

        do k=1,kbm1
          do j=2,jmm1
            ! east
            if(n_east.eq.-1) then
              denom=(uf(im-1,j,k)+ub(im-1,j,k)-2.e0*u(im-2,j,k))
              if(denom.eq.0.e0)denom=0.01e0
              cl=(ub(im-1,j,k)-uf(im-1,j,k))/denom
              if(cl.gt.1.e0) cl=1.e0
              if(cl.lt.0.e0) cl=0.e0
              uf(im,j,k)=(ub(im,j,k)*(1.e0-cl)+2.e0*cl*u(im-1,j,k))
     $                    /(1.e0+cl)
              vf(im,j,k)=0.e0
            end if

            ! west
            if(n_west.eq.-1) then
              denom=(uf(3,j,k)+ub(3,j,k)-2.e0*u(4,j,k))
              if(denom.eq.0.e0)denom=0.01e0
              cl=(ub(3,j,k)-uf(3,j,k))/denom
              if(cl.gt.1.e0) cl=1.e0
              if(cl.lt.0.e0) cl=0.e0
              uf(2,j,k)=(ub(2,j,k)*(1.e0-cl)+2.e0*cl*u(3,j,k))
     $                   /(1.e0+cl)
              uf(1,j,k)=uf(2,j,k)
              vf(1,j,k)=0.e0
            end if
          end do
        end do

        do k=1,kbm1
          do j=1,jm
            do i=1,im
              uf(i,j,k)=uf(i,j,k)*dum(i,j)
              vf(i,j,k)=vf(i,j,k)*dvm(i,j)
            end do
          end do
        end do

        return

      else if(idx.eq.4) then

! temperature and salinity boundary conditions (using uf and vf,
! respectively)
        do k=1,kbm1
          do j=1,jm
            ! east
            if(n_east.eq.-1) then
              ube(j,k)=ub(im,j,k)
              denom=(uf(im-1,j,k)+tb(im-1,j,k)-2.e0*t(im-2,j,k))
              if(denom.eq.0.e0) denom=0.01e0
              cl=(tb(im-1,j,k)-uf(im-1,j,k))/denom
              if(cl.gt.1.e0) cl=1.e0
              if(cl.lt.0.e0) cl=0.e0
              uf(im,j,k)=(tb(im,j,k)*(1.e0-cl)+2.e0*cl*t(im-1,j,k))
     $                    /(1.e0+cl)
              if(cl.eq.0.e0.and.ube(j,k).le.0.e0) uf(im,j,k)=tbe(j,k)

              denom=(vf(im-1,j,k)+sb(im-1,j,k)-2.e0*s(im-2,j,k))
              if(denom.eq.0.e0) denom=0.01e0
              cl=(sb(im-1,j,k)-vf(im-1,j,k))/denom
              if(cl.gt.1.e0) cl=1.e0
              if(cl.lt.0.e0) cl=0.e0
              vf(im,j,k)=(sb(im,j,k)*(1.e0-cl)+2.e0*cl*s(im-1,j,k))
     $                    /(1.e0+cl)
              if(cl.eq.0.e0.and.ube(j,k).le.0.e0) vf(im,j,k)=sbe(j,k)
            end if

            ! west
            if(n_west.eq.-1) then
              ubw(j,k)=ub(2,j,k)
              denom=(uf(2,j,k)+tb(2,j,k)-2.e0*t(3,j,k))
              if(denom.eq.0.e0) denom=0.01e0
              cl=(tb(2,j,k)-uf(2,j,k))/denom
              if(cl.gt.1.e0) cl=1.e0
              if(cl.lt.0.e0) cl=0.e0
              uf(1,j,k)=(tb(1,j,k)*(1.e0-cl)+2.e0*cl*t(2,j,k))/(1.e0+cl)
              if(cl.eq.0.e0.and.ubw(j,k).ge.0.e0) uf(1,j,k)=tbw(j,k)

              denom=(vf(2,j,k)+sb(2,j,k)-2.e0*s(3,j,k))
              if(denom.eq.0.e0) denom=0.01e0
              cl=(sb(2,j,k)-vf(2,j,k))/denom
              if(cl.gt.1.e0) cl=1.e0
              if(cl.lt.0.e0) cl=0.e0
              vf(1,j,k)=(sb(1,j,k)*(1.e0-cl)+2.e0*cl*s(2,j,k))/(1.e0+cl)
              if(cl.eq.0.e0.and.ubw(j,k).ge.0.e0) vf(1,j,k)=sbw(j,k)
            end if
          end do
        end do

        do k=1,kbm1
          do j=1,jm
            do i=1,im
              uf(i,j,k)=uf(i,j,k)*fsm(i,j)
              vf(i,j,k)=vf(i,j,k)*fsm(i,j)
            end do
          end do
        end do

        return

      else if(idx.eq.5) then

! vertical velocity boundary conditions
        do k=1,kbm1
          do j=1,jm
            do i=1,im
              w(i,j,k)=w(i,j,k)*fsm(i,j)
            end do
          end do
        end do

        return

      else if(idx.eq.6) then

! q2 and q2l boundary conditions
        do k=1,kb
          do j=1,jm
            if(n_east.eq.-1) then
              uf(im,j,k)=1.e-10
              vf(im,j,k)=1.e-10
            end if
            if(n_west.eq.-1) then
              uf(1,j,k)=1.e-10
              vf(1,j,k)=1.e-10
            end if
          end do

          do j=1,jm
            do i=1,im
              uf(i,j,k)=uf(i,j,k)*fsm(i,j)
              vf(i,j,k)=vf(i,j,k)*fsm(i,j)
            end do
          end do
        end do

        return

      endif

      end

