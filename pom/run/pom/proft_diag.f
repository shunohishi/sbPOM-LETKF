!_______________________________________________________________________
      subroutine proft_diag(f,wfsurf,fsurf,nbc,dh)
! solves for vertical diffusion of temperature and salinity using method
! described by Richmeyer and Morton (1967)
! note: wfsurf and swrad are negative values when water column is
! warming or salt is being added
      implicit none
      include 'pom.h'
      real f(im,jm,kb),wfsurf(im,jm)
      real fsurf(im,jm),dh(im,jm)

! diag
      include 'pom_diag.h'
      
      integer nbc
      real rad(im,jm,kb),r(5),ad1(5),ad2(5)
      integer i,j,k,ki

      integer iglb,jglb,iloc,jloc,nproc

! irradiance parameters after Paulson and Simpson (1977)
!       ntp               1      2       3       4       5
!   Jerlov type           i      ia      ib      ii     iii
      data r   /       .58e0,  .62e0,  .67e0,  .77e0,  .78e0 /
      data ad1 /       .35e0,  .60e0,  1.0e0,  1.5e0,  1.4e0 /
      data ad2 /       23.e0,  20.e0,  17.e0,  14.e0,  7.9e0 /

! surface boundary condition:
!       nbc   prescribed    prescribed   short wave
!             temperature      flux      penetration
!             or salinity               (temperature
!                                           only)
!        1        no           yes           no
!        2        no           yes           yes
!        3        yes          no            no
!        4        yes          no            yes
! note that only 1 and 3 are allowed for salinity.

! diag
      do k=1,kb
        do j=1,jm
          do i=1,im
            vdifterm(i,j,k)=f(i,j,k)
          end do
        end do
      end do

! the following section solves the equation
!     dti2*(kh*f')'-f=-fb
      do j=1,jm
        do i=1,im
          dh(i,j)=h(i,j)+etf(i,j)
        end do
      end do

      do k=2,kbm1
        do j=1,jm
          do i=1,im
            a(i,j,k-1)=-dti2*(kh(i,j,k)+umol)
     $                  /(dz(k-1)*dzz(k-1)*dh(i,j)*dh(i,j))
            c(i,j,k)=-dti2*(kh(i,j,k)+umol)
     $                  /(dz(k)*dzz(k-1)*dh(i,j)*dh(i,j))
          end do
        end do
      end do

! calculate penetrative radiation. At the bottom any unattenuated
! radiation is deposited in the bottom layer
      do k=1,kb
        do j=1,jm
          do i=1,im
            rad(i,j,k)=0.e0
          end do
        end do
      end do

      if(nbc.eq.2.or.nbc.eq.4) then
        do k=1,kbm1
          do j=1,jm
            do i=1,im
              rad(i,j,k)=swrad(i,j)
     $                    *(r(ntp)*exp(z(k)*dh(i,j)/ad1(ntp))
     $                      +(1.e0-r(ntp))*exp(z(k)*dh(i,j)/ad2(ntp)))
            end do
          end do
        end do
      endif

      if(nbc.eq.1) then

        do j=1,jm
          do i=1,im
            ee(i,j,1)=a(i,j,1)/(a(i,j,1)-1.e0)
            gg(i,j,1)=-dti2*wfsurf(i,j)/(-dz(1)*dh(i,j))-f(i,j,1)
            gg(i,j,1)=gg(i,j,1)/(a(i,j,1)-1.e0)
          end do
        end do

      else if(nbc.eq.2) then

        do j=1,jm
          do i=1,im
            ee(i,j,1)=a(i,j,1)/(a(i,j,1)-1.e0)
            gg(i,j,1)=dti2*(wfsurf(i,j)+rad(i,j,1)-rad(i,j,2))
     $                 /(dz(1)*dh(i,j))
     $                   -f(i,j,1)
            gg(i,j,1)=gg(i,j,1)/(a(i,j,1)-1.e0)
          end do
        end do

      else if(nbc.eq.3.or.nbc.eq.4) then

        do j=1,jm
          do i=1,im
            ee(i,j,1)=0.e0
            gg(i,j,1)=fsurf(i,j)
          end do
        end do

      endif

! debug
      iglb=(im_global-2)/2
      jglb=(jm_global-2)/2
      nproc=int((iglb-2)/(im-2))+int((jglb-2)/(jm-2))*nproc_x

! debug
!      if(my_task.eq.nproc .and. mod(iint-1,1).eq.0 .and.
!     $   nbc.eq.2) then
!        iloc=iglb-mod(nproc,nproc_x)*(im-2)
!        jloc=jglb-(nproc/nproc_x)*(jm-2)
!        write(99,*) 'in proft iint ',iint
!        write(99,*) 'iloc jloc nproc ',iloc,jloc,nproc
!        write(99,*) 'iglb jglb ',iglb,jglb
!        write(99,*) 'gg ee ',gg(iloc,jloc,1),ee(iloc,jloc,1)
!        write(99,*) 'rad swrad ',rad(iloc,jloc,1),swrad(iloc,jloc)
!        write(99,*) 'f wfsurf ',f(iloc,jloc,1),wfsurf(iloc,jloc)
!        write(99,*) 'wusurf wvsurf '
!     $   ,wusurf(iloc,jloc),wvsurf(iloc,jloc)
!        write(99,*) 'wtsurf wssurf ',wtsurf(iloc,jloc),wssurf(iloc,jloc)
!        call flush(99)
!      end if

      do k=2,kbm2
        do j=1,jm
          do i=1,im
            gg(i,j,k)=1.e0/(a(i,j,k)+c(i,j,k)*(1.e0-ee(i,j,k-1))-1.e0)
            ee(i,j,k)=a(i,j,k)*gg(i,j,k)
            gg(i,j,k)=(c(i,j,k)*gg(i,j,k-1)-f(i,j,k)
     $                 +dti2*(rad(i,j,k)-rad(i,j,k+1))
     $                   /(dh(i,j)*dz(k)))
     $                 *gg(i,j,k)
          end do
        end do
      end do

! bottom adiabatic boundary condition
      do j=1,jm
        do i=1,im
          f(i,j,kbm1)=(c(i,j,kbm1)*gg(i,j,kbm2)-f(i,j,kbm1)
     $                 +dti2*(rad(i,j,kbm1)-rad(i,j,kb))
     $                   /(dh(i,j)*dz(kbm1)))
     $                 /(c(i,j,kbm1)*(1.e0-ee(i,j,kbm2))-1.e0)
        end do
      end do

      do k=2,kbm1
        ki=kb-k
        do j=1,jm
          do i=1,im
          f(i,j,ki)=(ee(i,j,ki)*f(i,j,ki+1)+gg(i,j,ki))
          end do
        end do
      end do

! diag

      do j=2,jmm1
        do i=2,imm1
          sfcterm(i,j)=-wfsurf(i,j)/(dz(1)*dh(i,j))
        end do
      end do

      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            radterm(i,j,k)=-(rad(i,j,k)-rad(i,j,k+1))
     $                      /(dh(i,j)*dz(k))
            gvdterm(i,j,k)=(f(i,j,k)-vdifterm(i,j,k))/dti2
            vdifterm(i,j,k)=(f(i,j,k)-vdifterm(i,j,k))/dti2
     $                     -radterm(i,j,k)
            if(k.eq.1) then
               vdifterm(i,j,k)=vdifterm(i,j,k)-sfcterm(i,j)
            end if
          end do
        end do
      end do

      return
      end
