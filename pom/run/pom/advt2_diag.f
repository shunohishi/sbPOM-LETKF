!_______________________________________________________________________
      subroutine advt2_diag(fb,f,fclim,ff,xflux,yflux)
! integrate conservative scalar equations
! this is a first-order upstream scheme, which reduces implicit
! diffusion using the Smolarkiewicz iterative upstream scheme with an
! antidiffusive velocity
! it is based on the subroutines of Gianmaria Sannino (Inter-university
! Computing Consortium, Rome, Italy) and Vincenzo Artale (Italian
! National Agency for New Technology and Environment, Rome, Italy)
      implicit none
      include 'pom.h'
      real fb(im,jm,kb),f(im,jm,kb),fclim(im,jm,kb),ff(im,jm,kb)
      real xflux(im,jm,kb),yflux(im,jm,kb)
      real fbmem(im,jm,kb),eta(im,jm)
      real xmassflux(im,jm,kb),ymassflux(im,jm,kb),zwflux(im,jm,kb)

! diag
      include 'pom_diag.h'

      integer i,j,k,itera

! diag
      do k=1,kb
        do j=1,jm
          do i=1,im
            xadvterm(i,j,k)=0.e0
            yadvterm(i,j,k)=0.e0
            wadvterm(i,j,k)=0.e0
          end do
        end do
      end do

! calculate horizontal mass fluxes
      do k=1,kb
        do j=1,jm
          do i=1,im
            xmassflux(i,j,k)=0.e0
            ymassflux(i,j,k)=0.e0
          end do
        end do
      end do

      do k=1,kbm1
        do j=2,jmm1
          do i=2,im
            xmassflux(i,j,k)=0.25e0*(dy(i-1,j)+dy(i,j))
     $                             *(dt(i-1,j)+dt(i,j))*u(i,j,k)
          end do
        end do

        do j=2,jm
          do i=2,imm1
            ymassflux(i,j,k)=0.25e0*(dx(i,j-1)+dx(i,j))
     $                             *(dt(i,j-1)+dt(i,j))*v(i,j,k)
          end do
        end do
      end do

      do j=1,jm
        do i=1,im
          fb(i,j,kb)=fb(i,j,kbm1)
          eta(i,j)=etb(i,j)
        end do
      end do

      do k=1,kb
        do j=1,jm
          do i=1,im
            zwflux(i,j,k)=w(i,j,k)
            fbmem(i,j,k)=fb(i,j,k)
          end do
        end do
      end do

!! debug
!      do k=1,kb-1
!      do j=1,jm
!      do i=1,im
!        write(97,*) 'fbmem ',my_task,i,j,k,fbmem(i,j,k)
!      end do
!      end do
!      end do
!      do k=1,kb-1
!      do j=1,jm
!      do i=1,im
!        write(98,*) 'xmassflux ',my_task,i,j,k,xmassflux(i,j,k)
!      end do
!      end do
!      end do
!      do k=1,kb-1
!      do j=1,jm
!      do i=1,im
!        write(99,*) 'ymassflux ',my_task,i,j,k,ymassflux(i,j,k)
!      end do
!      end do
!      end do
!      call mpi_finalize
!      stop

! start Smolarkiewicz scheme
      do itera=1,nitera

! upwind advection scheme
        do k=1,kbm1
          do j=2,jm
            do i=2,im
              xflux(i,j,k)=0.5e0
     $                      *((xmassflux(i,j,k)+abs(xmassflux(i,j,k)))
     $                        *fbmem(i-1,j,k)+
     $                        (xmassflux(i,j,k)-abs(xmassflux(i,j,k)))
     $                        *fbmem(i,j,k))

              yflux(i,j,k)=0.5e0
     $                      *((ymassflux(i,j,k)+abs(ymassflux(i,j,k)))
     $                        *fbmem(i,j-1,k)+
     $                        (ymassflux(i,j,k)-abs(ymassflux(i,j,k)))
     $                        *fbmem(i,j,k))
            end do
          end do
        end do

        do j=2,jmm1
          do i=2,imm1
            zflux(i,j,1)=0.e0
            if(itera.eq.1) zflux(i,j,1)=w(i,j,1)*f(i,j,1)*art(i,j)
            zflux(i,j,kb)=0.e0
          end do
        end do

        do k=2,kbm1
          do j=2,jmm1
            do i=2,imm1
              zflux(i,j,k)=0.5e0
     $                      *((zwflux(i,j,k)+abs(zwflux(i,j,k)))
     $                       *fbmem(i,j,k)+
     $                        (zwflux(i,j,k)-abs(zwflux(i,j,k)))
     $                       *fbmem(i,j,k-1))
              zflux(i,j,k)=zflux(i,j,k)*art(i,j)
            end do
          end do
        end do

! add net advective fluxes and step forward in time
        do k=1,kbm1
          do j=2,jmm1
            do i=2,imm1
              ff(i,j,k)=xflux(i+1,j,k)-xflux(i,j,k)
     $                 +yflux(i,j+1,k)-yflux(i,j,k)
     $                 +(zflux(i,j,k)-zflux(i,j,k+1))/dz(k)
              ff(i,j,k)=(fbmem(i,j,k)*(h(i,j)+eta(i,j))*art(i,j)
     $                   -dti2*ff(i,j,k))/((h(i,j)+etf(i,j))*art(i,j))
            end do
          end do
        end do
        ! next line added on 22-Jul-2009 by Raffaele Bernardello
        call exchange3d_mpi(ff(:,:,1:kbm1),im,jm,kbm1)

! diag
        do k=1,kbm1
          do j=2,jmm1
            do i=2,imm1
              xadvterm(i,j,k)=xadvterm(i,j,k)
     $                 -(xflux(i+1,j,k)-xflux(i,j,k))
     $                  /((h(i,j)+etf(i,j))*art(i,j))
              yadvterm(i,j,k)=yadvterm(i,j,k)
     $                 -(yflux(i,j+1,k)-yflux(i,j,k))
     $                  /((h(i,j)+etf(i,j))*art(i,j))
              wadvterm(i,j,k)=wadvterm(i,j,k)
     $                 -(zflux(i,j,k)-zflux(i,j,k+1))/dz(k)
     $                  /((h(i,j)+etf(i,j))*art(i,j))
            end do
          end do
        end do

! calculate antidiffusion velocity
        call smol_adif(xmassflux,ymassflux,zwflux,ff)

        do j=1,jm
          do i=1,im
            eta(i,j)=etf(i,j)
            do k=1,kb
              fbmem(i,j,k)=ff(i,j,k)
            end do
          end do
        end do

! end of Smolarkiewicz scheme
      end do

! add horizontal diffusive fluxes
      do k=1,kb
        do j=1,jm
          do i=1,im
            fb(i,j,k)=fb(i,j,k)-fclim(i,j,k)
          end do
        end do
      end do

      do k=1,kbm1
        do j=2,jm
          do i=2,im
            xmassflux(i,j,k)=0.5e0*(aam(i,j,k)+aam(i-1,j,k))
            ymassflux(i,j,k)=0.5e0*(aam(i,j,k)+aam(i,j-1,k))
          end do
        end do
      end do

      do k=1,kbm1
        do j=2,jm
          do i=2,im
           xflux(i,j,k)=-xmassflux(i,j,k)*(h(i,j)+h(i-1,j))*tprni
     $                   *(fb(i,j,k)-fb(i-1,j,k))*dum(i,j)
     $                   *(dy(i,j)+dy(i-1,j))*0.5e0/(dx(i,j)+dx(i-1,j))
           yflux(i,j,k)=-ymassflux(i,j,k)*(h(i,j)+h(i,j-1))*tprni
     $                   *(fb(i,j,k)-fb(i,j-1,k))*dvm(i,j)
     $                   *(dx(i,j)+dx(i,j-1))*0.5e0/(dy(i,j)+dy(i,j-1))
          end do
        end do
      end do

      do k=1,kb
        do j=1,jm
          do i=1,im
            fb(i,j,k)=fb(i,j,k)+fclim(i,j,k)
          end do
        end do
      end do

! add net horizontal fluxes and step forward in time
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            ff(i,j,k)=ff(i,j,k)-dti2*(xflux(i+1,j,k)-xflux(i,j,k)
     $                               +yflux(i,j+1,k)-yflux(i,j,k))
     $                           /((h(i,j)+etf(i,j))*art(i,j))
          end do
        end do
      end do

! diag
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            xdifterm(i,j,k)=-(xflux(i+1,j,k)-xflux(i,j,k))
     $                       /((h(i,j)+etf(i,j))*art(i,j))
            ydifterm(i,j,k)=-(yflux(i,j+1,k)-yflux(i,j,k))
     $                       /((h(i,j)+etf(i,j))*art(i,j))
          end do
        end do
      end do

      return
      end
