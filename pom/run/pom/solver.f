! solver.f

! main subroutines for POM
!_______________________________________________________________________
! kii2b
!      subroutine advave(curv2d)
      subroutine advave
! calculate horizontal advection and diffusion
      implicit none
      include 'pom.h'
      real curv2d(im,jm)
      integer i,j

! u-advection and diffusion

! advective fluxes
      do j=1,jm
        do i=1,im
          advua(i,j)=0.e0
        end do
      end do

      do j=2,jm
        do i=2,imm1
          fluxua(i,j)=.125e0*((d(i+1,j)+d(i,j))*ua(i+1,j)
     $                       +(d(i,j)+d(i-1,j))*ua(i,j))
     $                      *(ua(i+1,j)+ua(i,j))
        end do
      end do

      do j=2,jm
        do i=2,im
          fluxva(i,j)=.125e0*((d(i,j)+d(i,j-1))*va(i,j)
     $                       +(d(i-1,j)+d(i-1,j-1))*va(i-1,j))
     $                      *(ua(i,j)+ua(i,j-1))
        end do
      end do

! add viscous fluxes
      do j=2,jm
        do i=2,imm1
          fluxua(i,j)=fluxua(i,j)
     $                 -d(i,j)*2.e0*aam2d(i,j)*(uab(i+1,j)-uab(i,j))
     $                   /dx(i,j)
        end do
      end do

      do j=2,jm
        do i=2,im
          tps(i,j)=.25e0*(d(i,j)+d(i-1,j)+d(i,j-1)+d(i-1,j-1))
     $              *(aam2d(i,j)+aam2d(i,j-1)
     $                +aam2d(i-1,j)+aam2d(i-1,j-1))
     $              *((uab(i,j)-uab(i,j-1))
     $                 /(dy(i,j)+dy(i-1,j)+dy(i,j-1)+dy(i-1,j-1))
     $               +(vab(i,j)-vab(i-1,j))
     $                 /(dx(i,j)+dx(i-1,j)+dx(i,j-1)+dx(i-1,j-1)))
          fluxua(i,j)=fluxua(i,j)*dy(i,j)
          fluxva(i,j)=(fluxva(i,j)-tps(i,j))*.25e0
     $                 *(dx(i,j)+dx(i-1,j)+dx(i,j-1)+dx(i-1,j-1))
        end do
      end do
      call exchange2d_mpi(fluxua,im,jm)

      do j=2,jmm1
        do i=2,imm1
          advua(i,j)=fluxua(i,j)-fluxua(i-1,j)
     $                +fluxva(i,j+1)-fluxva(i,j)
        end do
      end do

      call exchange2d_mpi(advua,im,jm) !2018.08

! v-advection and diffusion
      do j=1,jm
        do i=1,im
          advva(i,j)=0.e0
        end do
      end do

! advective fluxes
      do j=2,jm
        do i=2,im
          fluxua(i,j)=.125e0*((d(i,j)+d(i-1,j))*ua(i,j)
     $                       +(d(i,j-1)+d(i-1,j-1))*ua(i,j-1))
     $                      *(va(i-1,j)+va(i,j))
        end do
      end do

      do j=2,jmm1
        do i=2,im
          fluxva(i,j)=.125e0*((d(i,j+1)+d(i,j))*va(i,j+1)
     $                       +(d(i,j)+d(i,j-1))*va(i,j))
     $                      *(va(i,j+1)+va(i,j))
        end do
      end do

! add viscous fluxes
      do j=2,jmm1
        do i=2,im
          fluxva(i,j)=fluxva(i,j)
     $                 -d(i,j)*2.e0*aam2d(i,j)*(vab(i,j+1)-vab(i,j))
     $                   /dy(i,j)
        end do
      end do

      do j=2,jm
        do i=2,im
          fluxva(i,j)=fluxva(i,j)*dx(i,j)
          fluxua(i,j)=(fluxua(i,j)-tps(i,j))*.25e0
     $                 *(dy(i,j)+dy(i-1,j)+dy(i,j-1)+dy(i-1,j-1))
        end do
      end do
      call exchange2d_mpi(fluxva,im,jm)

      do j=2,jmm1
        do i=2,imm1
          advva(i,j)=fluxua(i+1,j)-fluxua(i,j)
     $                +fluxva(i,j)-fluxva(i,j-1)
        end do
      end do

      call exchange2d_mpi(advva,im,jm) !2018.08

      if(mode.eq.2) then

        do j=2,jmm1
          do i=2,imm1
            wubot(i,j)=-0.5e0*(cbc(i,j)+cbc(i-1,j))
     $                  *sqrt(uab(i,j)**2
     $                        +(.25e0*(vab(i,j)+vab(i,j+1)
     $                                 +vab(i-1,j)+vab(i-1,j+1)))**2)
     $                  *uab(i,j)
          end do
        end do

        call exchange2d_mpi(wubot,im,jm) !2018.07.28

        do j=2,jmm1
          do i=2,imm1
            wvbot(i,j)=-0.5e0*(cbc(i,j)+cbc(i,j-1))
     $                  *sqrt(vab(i,j)**2
     $                        +(.25e0*(uab(i,j)+uab(i+1,j)
     $                                +uab(i,j-1)+uab(i+1,j-1)))**2)
     $                  *vab(i,j)
          end do
        end do

        call exchange2d_mpi(wvbot,im,jm) !2018.07.28

        do j=2,jmm1
          do i=2,imm1
            curv2d(i,j)=.25e0
     $                   *((va(i,j+1)+va(i,j))*(dy(i+1,j)-dy(i-1,j))
     $                    -(ua(i+1,j)+ua(i,j))*(dx(i,j+1)-dx(i,j-1)))
     $                   /(dx(i,j)*dy(i,j))
          end do
        end do

        call exchange2d_mpi(curv2d,im,jm)

        do j=2,jmm1
          if(n_west.eq.-1) then
          do i=3,imm1
            advua(i,j)=advua(i,j)-aru(i,j)*.25e0
     $                  *(curv2d(i,j)*d(i,j)
     $                    *(va(i,j+1)+va(i,j))
     $                    +curv2d(i-1,j)*d(i-1,j)
     $                    *(va(i-1,j+1)+va(i-1,j)))
          end do
          else
          do i=2,imm1
            advua(i,j)=advua(i,j)-aru(i,j)*.25e0
     $                  *(curv2d(i,j)*d(i,j)
     $                    *(va(i,j+1)+va(i,j))
     $                    +curv2d(i-1,j)*d(i-1,j)
     $                    *(va(i-1,j+1)+va(i-1,j)))
          end do
          end if
        end do

        call exchange2d_mpi(advua,im,jm) !2018.08

        do i=2,imm1
          if(n_south.eq.-1) then
          do j=3,jmm1
            advva(i,j)=advva(i,j)+arv(i,j)*.25e0
     $                  *(curv2d(i,j)*d(i,j)
     $                    *(ua(i+1,j)+ua(i,j))
     $                    +curv2d(i,j-1)*d(i,j-1)
     $                    *(ua(i+1,j-1)+ua(i,j-1)))
          end do
          else
          do j=2,jmm1
            advva(i,j)=advva(i,j)+arv(i,j)*.25e0
     $                  *(curv2d(i,j)*d(i,j)
     $                    *(ua(i+1,j)+ua(i,j))
     $                    +curv2d(i,j-1)*d(i,j-1)
     $                    *(ua(i+1,j-1)+ua(i,j-1)))
          end do
          end if
        end do

        call exchange2d_mpi(advva,im,jm) !2018.08

      endif

      return
      end

!_______________________________________________________________________
! kii2b
!      subroutine advct(xflux,yflux,curv)
      subroutine advct
! calculate the horizontal portions of momentum advection well in
! advance of their use in advu and advv so that their vertical integrals
! (created in the main program) may be used in the external (2-D) mode
! calculation
      implicit none
      include 'pom.h'
      real xflux(im,jm,kb),yflux(im,jm,kb)
      real curv(im,jm,kb)
      real dtaam
      integer i,j,k

      do k=1,kb
        do j=1,jm
          do i=1,im
            curv(i,j,k)=0.e0
            advx(i,j,k)=0.e0
            xflux(i,j,k)=0.e0
            yflux(i,j,k)=0.e0
          end do
        end do
      end do

      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            curv(i,j,k)=.25e0*((v(i,j+1,k)+v(i,j,k))
     $                         *(dy(i+1,j)-dy(i-1,j))
     $                         -(u(i+1,j,k)+u(i,j,k))
     $                         *(dx(i,j+1)-dx(i,j-1)))
     $                       /(dx(i,j)*dy(i,j))
          end do
        end do
      end do

      call exchange3d_mpi(curv(:,:,1:kbm1),im,jm,kbm1)

! calculate x-component of velocity advection

! calculate horizontal advective fluxes
      do k=1,kbm1
        do j=1,jm
          do i=2,imm1
            xflux(i,j,k)=.125e0*((dt(i+1,j)+dt(i,j))*u(i+1,j,k)
     $                           +(dt(i,j)+dt(i-1,j))*u(i,j,k))
     $                         *(u(i+1,j,k)+u(i,j,k))
          end do
        end do
      end do

      do k=1,kbm1
        do j=2,jm
          do i=2,im
            yflux(i,j,k)=.125e0*((dt(i,j)+dt(i,j-1))*v(i,j,k)
     $                           +(dt(i-1,j)+dt(i-1,j-1))*v(i-1,j,k))
     $                         *(u(i,j,k)+u(i,j-1,k))
          end do
        end do
      end do

! add horizontal diffusive fluxes
      do k=1,kbm1
        do j=2,jm
          do i=2,imm1
            xflux(i,j,k)=xflux(i,j,k)
     $                    -dt(i,j)*aam(i,j,k)*2.e0
     $                    *(ub(i+1,j,k)-ub(i,j,k))/dx(i,j)
            dtaam=.25e0*(dt(i,j)+dt(i-1,j)+dt(i,j-1)+dt(i-1,j-1))
     $             *(aam(i,j,k)+aam(i-1,j,k)
     $               +aam(i,j-1,k)+aam(i-1,j-1,k))
            yflux(i,j,k)=yflux(i,j,k)
     $                    -dtaam*((ub(i,j,k)-ub(i,j-1,k))
     $                            /(dy(i,j)+dy(i-1,j)
     $                              +dy(i,j-1)+dy(i-1,j-1))
     $                            +(vb(i,j,k)-vb(i-1,j,k))
     $                            /(dx(i,j)+dx(i-1,j)
     $                              +dx(i,j-1)+dx(i-1,j-1)))

            xflux(i,j,k)=dy(i,j)*xflux(i,j,k)
            yflux(i,j,k)=.25e0*(dx(i,j)+dx(i-1,j)
     $                          +dx(i,j-1)+dx(i-1,j-1))*yflux(i,j,k)
          end do
        end do
      end do
      call exchange3d_mpi(xflux(:,:,1:kbm1),im,jm,kbm1)

! do horizontal advection
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            advx(i,j,k)=xflux(i,j,k)-xflux(i-1,j,k)
     $                   +yflux(i,j+1,k)-yflux(i,j,k)
          end do
        end do
      end do

      do k=1,kbm1
        do j=2,jmm1
          if(n_west.eq.-1) then
          do i=3,imm1
            advx(i,j,k)=advx(i,j,k)
     $                   -aru(i,j)*.25e0
     $                     *(curv(i,j,k)*dt(i,j)
     $                        *(v(i,j+1,k)+v(i,j,k))
     $                       +curv(i-1,j,k)*dt(i-1,j)
     $                        *(v(i-1,j+1,k)+v(i-1,j,k)))
          end do
          else
          do i=2,imm1
            advx(i,j,k)=advx(i,j,k)
     $                   -aru(i,j)*.25e0
     $                     *(curv(i,j,k)*dt(i,j)
     $                        *(v(i,j+1,k)+v(i,j,k))
     $                       +curv(i-1,j,k)*dt(i-1,j)
     $                        *(v(i-1,j+1,k)+v(i-1,j,k)))
          end do
          end if
        end do
      end do

! calculate y-component of velocity advection

      do k=1,kb
        do j=1,jm
          do i=1,im
            advy(i,j,k)=0.e0
            xflux(i,j,k)=0.e0
            yflux(i,j,k)=0.e0
          end do
        end do
      end do

! calculate horizontal advective fluxes
      do k=1,kbm1
        do j=2,jm
          do i=2,im
            xflux(i,j,k)=.125e0*((dt(i,j)+dt(i-1,j))*u(i,j,k)
     $                           +(dt(i,j-1)+dt(i-1,j-1))*u(i,j-1,k))
     $                         *(v(i,j,k)+v(i-1,j,k))
          end do
        end do
      end do

      do k=1,kbm1
        do j=2,jmm1
          do i=1,im
            yflux(i,j,k)=.125e0*((dt(i,j+1)+dt(i,j))*v(i,j+1,k)
     $                           +(dt(i,j)+dt(i,j-1))*v(i,j,k))
     $                         *(v(i,j+1,k)+v(i,j,k))
          end do
        end do
      end do

! add horizontal diffusive fluxes
      do k=1,kbm1
        do j=2,jmm1
          do i=2,im
            dtaam=.25e0*(dt(i,j)+dt(i-1,j)+dt(i,j-1)+dt(i-1,j-1))
     $             *(aam(i,j,k)+aam(i-1,j,k)
     $               +aam(i,j-1,k)+aam(i-1,j-1,k))
            xflux(i,j,k)=xflux(i,j,k)
     $                    -dtaam*((ub(i,j,k)-ub(i,j-1,k))
     $                            /(dy(i,j)+dy(i-1,j)
     $                              +dy(i,j-1)+dy(i-1,j-1))
     $                            +(vb(i,j,k)-vb(i-1,j,k))
     $                            /(dx(i,j)+dx(i-1,j)
     $                              +dx(i,j-1)+dx(i-1,j-1)))
            yflux(i,j,k)=yflux(i,j,k)
     $                    -dt(i,j)*aam(i,j,k)*2.e0
     $                    *(vb(i,j+1,k)-vb(i,j,k))/dy(i,j)

            xflux(i,j,k)=.25e0*(dy(i,j)+dy(i-1,j)
     $                          +dy(i,j-1)+dy(i-1,j-1))*xflux(i,j,k)
            yflux(i,j,k)=dx(i,j)*yflux(i,j,k)
          end do
        end do
      end do
      call exchange3d_mpi(yflux(:,:,1:kbm1),im,jm,kbm1)

! do horizontal advection
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            advy(i,j,k)=xflux(i+1,j,k)-xflux(i,j,k)
     $                   +yflux(i,j,k)-yflux(i,j-1,k)
          end do
        end do
      end do

      do k=1,kbm1
        do i=2,imm1
          if(n_south.eq.-1) then
          do j=3,jmm1
            advy(i,j,k)=advy(i,j,k)
     $                   +arv(i,j)*.25e0
     $                     *(curv(i,j,k)*dt(i,j)
     $                        *(u(i+1,j,k)+u(i,j,k))
     $                       +curv(i,j-1,k)*dt(i,j-1)
     $                        *(u(i+1,j-1,k)+u(i,j-1,k)))
          end do
          else
          do j=2,jmm1
            advy(i,j,k)=advy(i,j,k)
     $                   +arv(i,j)*.25e0
     $                     *(curv(i,j,k)*dt(i,j)
     $                        *(u(i+1,j,k)+u(i,j,k))
     $                       +curv(i,j-1,k)*dt(i,j-1)
     $                        *(u(i+1,j-1,k)+u(i,j-1,k)))
          end do
          end if
        end do
      end do

      return
      end

!_______________________________________________________________________
! kii2b
!      subroutine advq(qb,q,qf,xflux,yflux)
       subroutine advq(qb,q,qf)
! calculates horizontal advection and diffusion, and vertical advection
! for turbulent quantities
      implicit none
      include 'pom.h'
      real qb(im,jm,kb),q(im,jm,kb),qf(im,jm,kb)
      real xflux(im,jm,kb),yflux(im,jm,kb)
      integer i,j,k

! do horizontal advection
      do k=2,kbm1
        do j=2,jm
          do i=2,im
            xflux(i,j,k)=.125e0*(q(i,j,k)+q(i-1,j,k))
     $                    *(dt(i,j)+dt(i-1,j))*(u(i,j,k)+u(i,j,k-1))
            yflux(i,j,k)=.125e0*(q(i,j,k)+q(i,j-1,k))
     $                    *(dt(i,j)+dt(i,j-1))*(v(i,j,k)+v(i,j,k-1))
          end do
        end do
      end do

! do horizontal diffusion
      do k=2,kbm1
        do j=2,jm
          do i=2,im
            xflux(i,j,k)=xflux(i,j,k)
     $                    -.25e0*(aam(i,j,k)+aam(i-1,j,k)
     $                            +aam(i,j,k-1)+aam(i-1,j,k-1))
     $                          *(h(i,j)+h(i-1,j))
     $                          *(qb(i,j,k)-qb(i-1,j,k))*dum(i,j)
     $                          /(dx(i,j)+dx(i-1,j))
            yflux(i,j,k)=yflux(i,j,k)
     $                    -.25e0*(aam(i,j,k)+aam(i,j-1,k)
     $                            +aam(i,j,k-1)+aam(i,j-1,k-1))
     $                          *(h(i,j)+h(i,j-1))
     $                          *(qb(i,j,k)-qb(i,j-1,k))*dvm(i,j)
     $                          /(dy(i,j)+dy(i,j-1))
            xflux(i,j,k)=.5e0*(dy(i,j)+dy(i-1,j))*xflux(i,j,k)
            yflux(i,j,k)=.5e0*(dx(i,j)+dx(i,j-1))*yflux(i,j,k)
          end do
        end do
      end do

! kii2b
      call exchange3d_mpi(xflux(:,:,1:kbm1),im,jm,kbm1)
      call exchange3d_mpi(yflux(:,:,1:kbm1),im,jm,kbm1)

! do vertical advection, add flux terms, then step forward in time
      do k=2,kbm1
        do j=2,jmm1
          do i=2,imm1
            qf(i,j,k)=(w(i,j,k-1)*q(i,j,k-1)-w(i,j,k+1)*q(i,j,k+1))
     $                 *art(i,j)/(dz(k)+dz(k-1))
     $                 +xflux(i+1,j,k)-xflux(i,j,k)
     $                 +yflux(i,j+1,k)-yflux(i,j,k)
            qf(i,j,k)=((h(i,j)+etb(i,j))*art(i,j)
     $                 *qb(i,j,k)-dti2*qf(i,j,k))
     $                /((h(i,j)+etf(i,j))*art(i,j))
          end do
        end do
      end do

      return
      end

!_______________________________________________________________________
      subroutine advt1(fb,f,fclim,ff)
! integrate conservative scalar equations
! this is centred scheme, as originally provide in POM (previously
! called advt)
!     2019.04 S.Ohishi real --> double precision

      implicit none
      include 'pom.h'

      integer i,j,k

      real fb(im,jm,kb),f(im,jm,kb),fclim(im,jm,kb),ff(im,jm,kb)

      double precision dfb(im,jm,kb),df(im,jm,kb)
      double precision dfclim(im,jm,kb),dff(im,jm,kb)

      double precision xflux(im,jm,kb),yflux(im,jm,kb),zflux(im,jm,kb)
      double precision xyzflux(im,jm,kb)

      double precision dfmx(im,jm,kb),dtmx(im,jm),dxm(im,jm) !f(i-1)--- fmx(i)/u(i)--- f(i) 
      double precision dfmy(im,jm,kb),dtmy(im,jm),dym(im,jm) !f(j-1)--- fmy(j)/v(j)--- f(j) 
      double precision dfmz(im,jm,kb) !f(k-1)---fmz(k)/w(k)---f(k)
      double precision dtmt(im,jm) !h+efb --- dtmt --- h+etf

!     debug
      double precision tmp(im,jm,kb)

!     Initialization
      do k=1,kb
         do j=1,jm
            do i=1,im
               dfb(i,j,k)=dble(fb(i,j,k))
               df(i,j,k)=dble(f(i,j,k))
               dfclim(i,j,k)=dble(fclim(i,j,k))
               ff(i,j,k)=0.e0
               dff(i,j,k)=0.d0
            end do
         end do
      end do

      do j=1,jm
         do i=1,im
            df(i,j,kb)=df(i,j,kbm1)
            dfb(i,j,kb)=dfb(i,j,kbm1)
         end do
      end do

      if(budget == 1)then
         do k=1,kb
            do j=1,jm
               do i=1,im
                  xadvterm(i,j,k)=0.e0
                  yadvterm(i,j,k)=0.e0
                  zadvterm(i,j,k)=0.e0
                  advterm(i,j,k)=0.e0
                  xdifterm(i,j,k)=0.e0
                  ydifterm(i,j,k)=0.e0
               end do
            end do
         end do
      end if

!     Intermidiate value
      do k=1,kbm1
         do j=2,jm
            do i=2,im
               dfmx(i,j,k)=0.5d0*(df(i-1,j,k)+df(i,j,k))
               dfmy(i,j,k)=0.5d0*(df(i,j-1,k)+df(i,j,k))
            end do
         end do
      end do

      do k=2,kbm1
         do j=2,jmm1
            do i=2,imm1
               dfmz(i,j,k)=0.5d0*(df(i,j,k-1)+df(i,j,k))
            end do
         end do
      end do
      
      do j=2,jmm1
         do i=2,imm1
            dfmz(i,j,1)=df(i,j,1)
            dfmz(i,j,kb)=0.d0
         end do
      end do

      do j=2,jm
         do i=2,im
            dtmx(i,j)=0.5d0*(dble(dt(i-1,j))+dble(dt(i,j)))
            dtmy(i,j)=0.5d0*(dble(dt(i,j-1))+dble(dt(i,j)))
            dtmt(i,j)=dble(h(i,j))
     $           +0.5d0*(dble(etb(i,j))+dble(etf(i,j)))
            dxm(i,j)=0.5d0*(dble(dx(i-1,j))+dble(dx(i,j)))
            dym(i,j)=0.5d0*(dble(dy(i,j-1))+dble(dy(i,j)))
         end do
      end do

!---  do horizontal fluxes
      do k=1,kbm1
         do j=2,jm
            do i=2,im               
               xflux(i,j,k)=
     $              dtmx(i,j)*dfmx(i,j,k)*dble(u(i,j,k))*dym(i,j) !UDf*dy
               yflux(i,j,k)=
     $              dtmy(i,j)*dfmy(i,j,k)*dble(v(i,j,k))*dxm(i,j) !VDf*dx
            end do
         end do
      end do

!     do vertical advection
      do k=1,kb
         do j=2,jmm1
            do i=2,imm1
               zflux(i,j,k)=dfmz(i,j,k)*dble(w(i,j,k))*dble(art(i,j)) !Wf*dxdy
            end do
         end do
      end do

! add advection and then step forward in time
      do k=1,kbm1
         do j=2,jmm1
            do i=2,imm1
               xyzflux(i,j,k)=
     $              (xflux(i+1,j,k)-xflux(i,j,k))
     $              +(yflux(i,j+1,k)-yflux(i,j,k))
     $              +(zflux(i,j,k)-zflux(i,j,k+1))/dble(dz(k))
               dff(i,j,k)=
     $              (dfb(i,j,k)*(dble(h(i,j))+dble(etb(i,j))) !Dn-1*Tn-1*dxdy
     $              *dble(art(i,j))
     $              -1.d0*dble(dti2)*xyzflux(i,j,k)) !adv
     $              /((dble(h(i,j))+dble(etf(i,j)))*dble(art(i,j))) !/(Dn+1*dxdy)
            end do
         end do
      end do

!     exchanged3d_mpi
      call exchange3d_mpi_dble(dff(:,:,1:kbm1),im,jm,kbm1)
      
      if(budget == 1)then
         
         do k=1,kbm1
            do j=2,jmm1
               do i=2,imm1
                  
                  xadvterm(i,j,k)=
     $                 -1.d0*(xflux(i+1,j,k)-xflux(i,j,k)) !-d(UDT)*dy
     $                 +0.5d0*(dfmx(i+1,j,k)+dfmx(i,j,k))
     $                 *(dtmx(i+1,j)*dble(u(i+1,j,k))*dym(i+1,j)
     $                 -dtmx(i,j)*dble(u(i,j,k))*dym(i,j)) !Td(UD)*dy
                  xadvterm(i,j,k)=dble(fsm(i,j))*xadvterm(i,j,k)
     $                 /((dble(h(i,j))+dble(etf(i,j)))*dble(art(i,j))) !-UdT/dx = (-d(UDT)*dy+Td(UD)*dy)/(D*dxdy)
                  
                  yadvterm(i,j,k)=
     $                 -1.d0*(yflux(i,j+1,k)-yflux(i,j,k)) !-d(VDT)*dx
     $                 +0.5d0*(dfmy(i,j+1,k)+dfmy(i,j,k))
     $                 *(dtmy(i,j+1)*dble(v(i,j+1,k))*dxm(i,j+1)
     $                 -dtmy(i,j)*dble(v(i,j,k))*dxm(i,j)) !Td(VD)*dx
                  yadvterm(i,j,k)=dble(fsm(i,j))*yadvterm(i,j,k)
     $                 /((dble(h(i,j))+dble(etf(i,j)))*dble(art(i,j))) !-VdT/dy = (-d(VDT)*dx+Td(VD)*dx)/(D*dxdy)
                  
                  zadvterm(i,j,k)=
     $                 -1.d0*(zflux(i,j,k)-zflux(i,j,k+1)) !-d(wT)*dxdy
     $                 +0.5d0*(dfmz(i,j,k+1)+dfmz(i,j,k))
     $                 *(w(i,j,k)-w(i,j,k+1))*dble(art(i,j)) !Tdw*dxdy
                  zadvterm(i,j,k)=dble(fsm(i,j))*zadvterm(i,j,k)
     $                 /((dble(h(i,j))+dble(etf(i,j)))
     $                 *dble(art(i,j)*dz(k))) !-wdT/dz=(-d(wT)*dxdy+Tdw*dxdy)/(D*dxdy*dz)

                  advterm(i,j,k)=
     $                 dble(fsm(i,j))*(dff(i,j,k)-dfb(i,j,k))/dble(dti2)

               end do
            end do
         end do
      end if
      
!--- add horizontal diffusive fluxes
      do k=1,kb
         do j=1,jm
            do i=1,im
               dfb(i,j,k)=dfb(i,j,k)-dfclim(i,j,k)
            end do
         end do
      end do

      do k=1,kbm1
         do j=2,jm
            do i=2,im
               xflux(i,j,k)=
     $              -0.5d0*(aam(i,j,k)+aam(i-1,j,k))
     $              *(h(i,j)+h(i-1,j))*tprni
     $              *(dfb(i,j,k)-dfb(i-1,j,k))*dum(i,j)
     $              *(dy(i,j)+dy(i-1,j))*0.5d0/(dx(i,j)+dx(i-1,j))
               yflux(i,j,k)=
     $              -0.5d0*(aam(i,j,k)+aam(i,j-1,k))
     $              *(h(i,j)+h(i,j-1))*tprni
     $              *(dfb(i,j,k)-dfb(i,j-1,k))*dvm(i,j)
     $              *(dx(i,j)+dx(i,j-1))*0.5d0/(dy(i,j)+dy(i,j-1))
            end do
         end do
      end do

      do k=1,kbm1
         do j=2,jmm1
            do i=2,imm1
               dff(i,j,k)=dff(i,j,k)
     $              -1.d0*dble(dti2)
     $              *(xflux(i+1,j,k)-xflux(i,j,k)
     $              +yflux(i,j+1,k)-yflux(i,j,k))
     $              /((dble(h(i,j))+dble(etf(i,j)))*dble(art(i,j)))
            end do
         end do
      end do

      do k=1,kb
         do j=1,jm
            do i=1,im
               dfb(i,j,k)=dfb(i,j,k)+dfclim(i,j,k)
            end do
         end do
      end do
      
!     conserve horizontal diffusion
      if(budget == 1)then
         do k=1,kbm1
            do j=2,jmm1
               do i=2,imm1
                  xdifterm(i,j,k)=-(xflux(i+1,j,k)-xflux(i,j,k))
     $                 /((h(i,j)+etf(i,j))*art(i,j))
                  ydifterm(i,j,k)=-(yflux(i,j+1,k)-yflux(i,j,k))
     $                 /((h(i,j)+etf(i,j))*art(i,j))
               end do
            end do
         end do
      end if
      
!     exchanged3d_mpi
      call exchange3d_mpi_dble(dff(:,:,1:kbm1),im,jm,kbm1)
      
      if(budget == 1)then
         call exchange3d_mpi(xadvterm(:,:,1:kbm1),im,jm,kbm1)
         call exchange3d_mpi(yadvterm(:,:,1:kbm1),im,jm,kbm1)
         call exchange3d_mpi(zadvterm(:,:,1:kbm1),im,jm,kbm1)
         call exchange3d_mpi(advterm(:,:,1:kbm1),im,jm,kbm1)
         call exchange3d_mpi(xdifterm(:,:,1:kbm1),im,jm,kbm1)
         call exchange3d_mpi(ydifterm(:,:,1:kbm1),im,jm,kbm1)
      end if
      
!     dff --> ff
      do k=1,kb
         do j=1,jm
            do i=1,im
               ff(i,j,k)=real(dff(i,j,k))
            end do
         end do
      end do

      return
      end

!_______________________________________________________________________
      subroutine advt2(fb,f,fclim,ff)
!     integrate conservative scalar equations
!     this is a first-order upstream scheme, which reduces implicit
!     diffusion using the Smolarkiewicz iterative upstream scheme with an
!     antidiffusive velocity
!     it is based on the subroutines of Gianmaria Sannino (Inter-university
!     Computing Consortium, Rome, Italy) and Vincenzo Artale (Italian
!     National Agency for New Technology and Environment, Rome, Italy)
!     S.Ohishi modified 2019.04

      implicit none
      include 'pom.h'
      
      integer i,j,k,itera
      
      real fb(im,jm,kb),f(im,jm,kb),fclim(im,jm,kb),ff(im,jm,kb)
      double precision xflux(im,jm,kb),yflux(im,jm,kb)
      double precision zflux(im,jm,kb),xyzflux(im,jm,kb)
      double precision dfbmem(im,jm,kb),eta(im,jm)
      double precision xmassflux(im,jm,kb),ymassflux(im,jm,kb)
      double precision zwflux(im,jm,kb)
      
      double precision dfb(im,jm,kb),df(im,jm,kb)
      double precision dfclim(im,jm,kb),dff(im,jm,kb)
     
      double precision dxm(im,jm),dym(im,jm)
      double precision dtmx(im,jm),dtmy(im,jm)
 
      double precision tmp(im,jm,kb)
      
!       Initialization
      do k=1,kb
         do j=1,jm
            do i=1,im
               tmp(i,j,k)=0.d0
               dfb(i,j,k)=dble(fb(i,j,k))
               df(i,j,k)=dble(f(i,j,k))
               dfclim(i,j,k)=dble(fclim(i,j,k))
            end do
         end do
      end do
            
!       Preparation
      do j=1,jm
         do i=1,im
            dfb(i,j,kb)=dfb(i,j,kbm1)
         end do
      end do

      do j=2,jm
         do i=2,im
            dxm(i,j)=0.5d0*(dble(dx(i,j-1))+dble(dx(i,j)))
            dym(i,j)=0.5d0*(dble(dy(i-1,j))+dble(dy(i,j)))
            dtmx(i,j)=0.5d0*(dble(dt(i-1,j))+dble(dt(i,j)))
            dtmy(i,j)=0.5d0*(dble(dt(i,j-1))+dble(dt(i,j)))
         end do
      end do            
      
!     calculate horizontal mass fluxes
      do k=1,kb
         do j=1,jm
            do i=1,im
               xmassflux(i,j,k)=0.d0
               ymassflux(i,j,k)=0.d0
            end do
         end do
      end do
      
      do k=1,kbm1
         do j=2,jmm1
            do i=2,im
!               xmassflux(i,j,k)=
!     $              0.25d0*(dble(dy(i-1,j))+dble(dy(i,j)))
!     $              *(dble(dt(i-1,j))+dble(dt(i,j)))*dble(u(i,j,k))
               xmassflux(i,j,k)=dym(i,j)*dtmx(i,j)*dble(u(i,j,k))
            end do
         end do
         
         do j=2,jm
            do i=2,imm1
!               ymassflux(i,j,k)=
!     $              0.25d0*(dble(dx(i,j-1))+dble(dx(i,j)))
!     $              *(dble(dt(i,j-1))+dble(dt(i,j)))*dble(v(i,j,k))
               ymassflux(i,j,k)=dxm(i,j)*dtmy(i,j)*dble(v(i,j,k))
            end do
         end do
      end do
      
      do j=1,jm
         do i=1,im
            dfb(i,j,kb)=dfb(i,j,kbm1)
            eta(i,j)=dble(etb(i,j))
         end do
      end do
      
      do k=1,kb
         do j=1,jm
            do i=1,im
               zwflux(i,j,k)=dble(w(i,j,k))
               dfbmem(i,j,k)=dble(fb(i,j,k))
            end do
         end do
      end do
      
!     Conservation
      if(budget == 1)then
         do k=1,kb
            do j=1,jm
               do i=1,im
                  ff(i,j,k)=0.e0
                  dff(i,j,k)=0.d0
                  xadvterm(i,j,k)=0.e0
                  yadvterm(i,j,k)=0.e0
                  zadvterm(i,j,k)=0.e0
                  xdifterm(i,j,k)=0.e0
                  ydifterm(i,j,k)=0.e0
               end do
            end do
         end do
      end if
      
!     start Smolarkiewicz scheme
      do itera=1,nitera
         
!     upwind advection scheme
         do k=1,kbm1
            do j=2,jm
               do i=2,im
                  xflux(i,j,k)=0.5d0
     $                 *((xmassflux(i,j,k)+dabs(xmassflux(i,j,k)))
     $                 *dfbmem(i-1,j,k)+
     $                 (xmassflux(i,j,k)-dabs(xmassflux(i,j,k)))
     $                 *dfbmem(i,j,k))
                  
                  yflux(i,j,k)=0.5d0
     $                 *((ymassflux(i,j,k)+dabs(ymassflux(i,j,k)))
     $                 *dfbmem(i,j-1,k)+
     $                 (ymassflux(i,j,k)-dabs(ymassflux(i,j,k)))
     $                 *dfbmem(i,j,k))
               end do
            end do
         end do
         
         do j=2,jmm1
            do i=2,imm1
               zflux(i,j,1)=0.d0
               if(itera.eq.1)then
!                  zflux(i,j,1)=dble(w(i,j,1))*df(i,j,1)*dble(art(i,j))
                  zflux(i,j,1)=dble(w(i,j,1))*dfb(i,j,1)*dble(art(i,j))
               end if
               zflux(i,j,kb)=0.d0
            end do
         end do
         
         do k=2,kbm1
            do j=2,jmm1
               do i=2,imm1
                  zflux(i,j,k)=0.5d0
     $                 *((zwflux(i,j,k)+dabs(zwflux(i,j,k)))
     $                 *dfbmem(i,j,k)+
     $                 (zwflux(i,j,k)-dabs(zwflux(i,j,k)))
     $                 *dfbmem(i,j,k-1))
                  zflux(i,j,k)=zflux(i,j,k)*dble(art(i,j))
               end do
            end do
         end do
         
!     add net advective fluxes and step forward in time
         do k=1,kbm1
            do j=2,jmm1
               do i=2,imm1
                  
!     Original             
!     ff(i,j,k)=xflux(i+1,j,k)-xflux(i,j,k)
!     $                 +yflux(i,j+1,k)-yflux(i,j,k)
!     $                 +(zflux(i,j,k)-zflux(i,j,k+1))/dz(k)
!     ff(i,j,k)=(fbmem(i,j,k)*(h(i,j)+eta(i,j))*art(i,j)
!     $                   -dti2*ff(i,j,k))/((h(i,j)+etf(i,j))*art(i,j))
                  
                  xyzflux(i,j,k)=
     $                 xflux(i+1,j,k)-xflux(i,j,k)
     $                 +yflux(i,j+1,k)-yflux(i,j,k)
     $                 +(zflux(i,j,k)-zflux(i,j,k+1))/dz(k)
                  dff(i,j,k)=
     $ (dfbmem(i,j,k)*(dble(h(i,j))+dble(eta(i,j)))*dble(art(i,j))
     $ -dble(dti2)*dble(xyzflux(i,j,k)))
     $              /((dble(h(i,j))+dble(etf(i,j)))*dble(art(i,j)))
                                    
               end do
            end do
         end do

! next line added on 22-Jul-2009 by Raffaele Bernardello
         call exchange3d_mpi_dble(dff(:,:,1:kbm1),im,jm,kbm1)

!     Conservation
         if(budget == 1 .and. itera == 1)then

            do k=1,kbm1
               do j=2,jmm1
                  do i=2,imm1
                                          
                     xadvterm(i,j,k)=
     $                    -1.d0*(xflux(i+1,j,k)-xflux(i,j,k)) !-d(UDT)*dy
     $                    +dfb(i,j,k)
     $                    *(dble(u(i+1,j,k))*dtmx(i+1,j)*dym(i+1,j)
     $                    -dble(u(i,j,k))*dtmx(i,j)*dym(i,j)) !Td(UD)*dy
                     xadvterm(i,j,k)=
     $                    dble(fsm(i,j))*xadvterm(i,j,k)/dble(art(i,j)) 
!     -UdT/dx*D
                     xadvterm(i,j,k)=
     $                    xadvterm(i,j,k)/(dble(h(i,j))+dble(etf(i,j))) 
!     -UdT/dx
                     
                     yadvterm(i,j,k)=
     $                    -1.d0*(yflux(i,j+1,k)-yflux(i,j,k)) !-d(VDT)*dx
     $                    +dfb(i,j,k)
     $                    *(dble(v(i,j+1,k))*dtmy(i,j+1)*dxm(i,j+1) 
     $                    -dble(v(i,j,k))*dtmy(i,j)*dxm(i,j)) !Td(VD)*dx
                     yadvterm(i,j,k)=
     $                    dble(fsm(i,j))*yadvterm(i,j,k)/dble(art(i,j)) 
!     -VdT/dy*D
                     yadvterm(i,j,k)=
     $                    yadvterm(i,j,k)/(dble(h(i,j))+dble(etf(i,j))) 
!     -VdT/dy
                     
                     zadvterm(i,j,k)=
     $                    -1.d0*(zflux(i,j,k)-zflux(i,j,k+1)) !-d(wT)*dxdy
     $                    +dfb(i,j,k)
     $                    *(dble(w(i,j,k))-dble(w(i,j,k+1)))
     $                    *dble(art(i,j)) !Tdw*dxdy
                     zadvterm(i,j,k)=
     $                    dble(fsm(i,j))*zadvterm(i,j,k)
     $                    /(dble(art(i,j))*dble(dz(k))) !-wdT/dz
                     zadvterm(i,j,k)=
     $                    zadvterm(i,j,k)/(dble(h(i,j))+dble(etf(i,j))) 
!     -w/D dT/dz
                     
                     advterm(i,j,k)=
     $                    dble(fsm(i,j))*(dff(i,j,k)-dfb(i,j,k))
     $                    /dble(dti2)
                     
                  end do
               end do
            end do
         else if(budget == 1)then
            do k=1,kbm1
               do j=2,jmm1
                  do i=2,imm1
                     
                     xadvterm(i,j,k)=xadvterm(i,j,k)
     $                    -1.d0*(xflux(i+1,j,k)-xflux(i,j,k))
     $                    /((dble(h(i,j))+dble(etf(i,j)))
     $                    *dble(art(i,j)))
                     
                     yadvterm(i,j,k)=yadvterm(i,j,k)
     $                    -1.d0*(yflux(i,j+1,k)-yflux(i,j,k))         
     $                    /((dble(h(i,j))+dble(etf(i,j)))
     $                    *dble(art(i,j)))
                     
                     zadvterm(i,j,k)=zadvterm(i,j,k)
     $                    -1.d0*(zflux(i,j,k)-zflux(i,j,k+1))
     $                    /((dble(h(i,j))+dble(etf(i,j)))*dble(dz(k))
     $                    *dble(art(i,j)))
                     
                     advterm(i,j,k)=
     $                    dble(fsm(i,j))*(dff(i,j,k)-dfb(i,j,k))
     $                    /dble(dti2)
                     
                  end do
               end do
            end do            
         end if
         
         if(budget == 1)then
            call exchange3d_mpi(xadvterm(:,:,1:kbm1),im,jm,kbm1)
            call exchange3d_mpi(yadvterm(:,:,1:kbm1),im,jm,kbm1)
            call exchange3d_mpi(zadvterm(:,:,1:kbm1),im,jm,kbm1)
            call exchange3d_mpi(advterm(:,:,1:kbm1),im,jm,kbm1)
         end if

!     calculate antidiffusion velocity
         call smol_adif_dble(xmassflux,ymassflux,zwflux,dff)
         
         do j=1,jm
            do i=1,im
               eta(i,j)=etf(i,j)
               do k=1,kb
                  dfbmem(i,j,k)=dff(i,j,k)
               end do
            end do
         end do
         
!     end of Smolarkiewicz scheme
      end do
      
      do k=1,kb
         do j=1,jm
            do i=1,im
               ff(i,j,k)=real(dff(i,j,k))
               tmp(i,j,k)=dff(i,j,k)
            end do
         end do
      end do
            
! add horizontal diffusive fluxes
      do k=1,kb
         do j=1,jm
            do i=1,im
               dfb(i,j,k)=dfb(i,j,k)-dfclim(i,j,k)
            end do
         end do
      end do
      
      do k=1,kbm1
         do j=2,jm
            do i=2,im
               xflux(i,j,k)=
     $              -0.5d0*(dble(aam(i,j,k))+dble(aam(i-1,j,k)))
     $              *(dble(h(i,j))+dble(h(i-1,j)))*dble(tprni)
     $              *(dfb(i,j,k)-dfb(i-1,j,k))*dble(dum(i,j))
     $              *(dble(dy(i,j))+dble(dy(i-1,j)))*0.5d0
     $              /(dble(dx(i,j))+dble(dx(i-1,j)))
               yflux(i,j,k)=
     $              -0.5d0*(dble(aam(i,j,k))+dble(aam(i,j-1,k)))
     $              *(dble(h(i,j))+dble(h(i,j-1)))*dble(tprni)
     $              *(dfb(i,j,k)-dfb(i,j-1,k))*dble(dvm(i,j))
     $              *(dble(dx(i,j))+dble(dx(i,j-1)))*0.5d0
     $              /(dble(dy(i,j))+dble(dy(i,j-1)))
            end do
         end do
      end do
      
      do k=1,kb
         do j=1,jm
            do i=1,im
               dfb(i,j,k)=dfb(i,j,k)+dfclim(i,j,k)
            end do
         end do
      end do
      
!     add net horizontal fluxes and step forward in time
      do k=1,kbm1
         do j=2,jmm1
            do i=2,imm1
               dff(i,j,k)=dff(i,j,k)
     $              -dble(dti2)*(xflux(i+1,j,k)-xflux(i,j,k)
     $              +yflux(i,j+1,k)-yflux(i,j,k))
     $              /((dble(h(i,j))+dble(etf(i,j)))*dble(art(i,j)))
            end do
         end do
      end do

!     Conservation
      if(budget == 1)then
         do k=1,kbm1
            do j=2,jmm1
               do i=2,imm1
                  xdifterm(i,j,k)=-(xflux(i+1,j,k)-xflux(i,j,k))
     $                 /((h(i,j)+etf(i,j))*art(i,j))
                  ydifterm(i,j,k)=-(yflux(i,j+1,k)-yflux(i,j,k))
     $                 /((h(i,j)+etf(i,j))*art(i,j))
               end do
            end do
         end do
         call exchange3d_mpi(xdifterm(:,:,1:kbm1),im,jm,kbm1)
         call exchange3d_mpi(ydifterm(:,:,1:kbm1),im,jm,kbm1)
      end if

      do k=1,kb
         do j=1,jm
            do i=1,im
               ff(i,j,k)=real(dff(i,j,k))
            end do
         end do
      end do

      return
      end

!_______________________________________________________________________
      subroutine advt1_original(fb,f,fclim,ff)
! integrate conservative scalar equations
! this is centred scheme, as originally provide in POM (previously
! called advt)
      implicit none
      include 'pom.h'
      real fb(im,jm,kb),f(im,jm,kb),fclim(im,jm,kb),ff(im,jm,kb)
      real xflux(im,jm,kb),yflux(im,jm,kb),zflux(im,jm,kb)
      integer i,j,k

      do j=1,jm
        do i=1,im
           f(i,j,kb)=f(i,j,kbm1)
           fb(i,j,kb)=fb(i,j,kbm1)
        end do
      end do

! do advective fluxes
      do k=1,kbm1
        do j=2,jm
          do i=2,im
            xflux(i,j,k)=.25e0*((dt(i,j)+dt(i-1,j))
     $                          *(f(i,j,k)+f(i-1,j,k))*u(i,j,k))
            yflux(i,j,k)=.25e0*((dt(i,j)+dt(i,j-1))
     $                          *(f(i,j,k)+f(i,j-1,k))*v(i,j,k))
          end do
        end do
      end do

! add diffusive fluxes
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
            xflux(i,j,k)=xflux(i,j,k)
     $                    -.5e0*(aam(i,j,k)+aam(i-1,j,k))
     $                         *(h(i,j)+h(i-1,j))*tprni
     $                         *(fb(i,j,k)-fb(i-1,j,k))*dum(i,j)
     $                         /(dx(i,j)+dx(i-1,j))
            yflux(i,j,k)=yflux(i,j,k)
     $                    -.5e0*(aam(i,j,k)+aam(i,j-1,k))
     $                         *(h(i,j)+h(i,j-1))*tprni
     $                         *(fb(i,j,k)-fb(i,j-1,k))*dvm(i,j)
     $                         /(dy(i,j)+dy(i,j-1))
            xflux(i,j,k)=.5e0*(dy(i,j)+dy(i-1,j))*xflux(i,j,k)
            yflux(i,j,k)=.5e0*(dx(i,j)+dx(i,j-1))*yflux(i,j,k)
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

! do vertical advection
      do j=2,jmm1
        do i=2,imm1
          zflux(i,j,1)=f(i,j,1)*w(i,j,1)*art(i,j)
          zflux(i,j,kb)=0.e0
        end do
      end do

      do k=2,kbm1
        do j=2,jmm1
          do i=2,imm1
            zflux(i,j,k)=.5e0*(f(i,j,k-1)+f(i,j,k))*w(i,j,k)*art(i,j)
          end do
        end do
      end do

! add net horizontal fluxes and then step forward in time
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            ff(i,j,k)=xflux(i+1,j,k)-xflux(i,j,k)
     $                 +yflux(i,j+1,k)-yflux(i,j,k)
     $                 +(zflux(i,j,k)-zflux(i,j,k+1))/dz(k)

            ff(i,j,k)=(fb(i,j,k)*(h(i,j)+etb(i,j))*art(i,j)
     $                 -dti2*ff(i,j,k))
     $                 /((h(i,j)+etf(i,j))*art(i,j))
          end do
        end do
      end do

      return
      end

!_______________________________________________________________________
      subroutine advt2_original(fb,f,fclim,ff)
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
      real xflux(im,jm,kb),yflux(im,jm,kb),zflux(im,jm,kb)
      real fbmem(im,jm,kb),eta(im,jm)
      real xmassflux(im,jm,kb),ymassflux(im,jm,kb),zwflux(im,jm,kb)
      integer i,j,k,itera

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

      return
      end

!_______________________________________________________________________
      subroutine advu
! do horizontal and vertical advection of u-momentum, and includes
! coriolis, surface slope and baroclinic terms
      implicit none
      include 'pom.h'
      integer i,j,k

! do vertical advection
      do k=1,kb
        do j=1,jm
          do i=1,im
            uf(i,j,k)=0.e0
          end do
        end do
      end do

      do k=2,kbm1
        do j=1,jm
          do i=2,im
            uf(i,j,k)=.25e0*(w(i,j,k)+w(i-1,j,k))
     $                     *(u(i,j,k)+u(i,j,k-1))
          end do
        end do
      end do

! combine horizontal and vertical advection with coriolis, surface
! slope and baroclinic terms
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            uf(i,j,k)=advx(i,j,k)
     $                 +(uf(i,j,k)-uf(i,j,k+1))*aru(i,j)/dz(k)
     $                 -aru(i,j)*.25e0
     $                   *(cor(i,j)*dt(i,j)
     $                      *(v(i,j+1,k)+v(i,j,k))
     $                     +cor(i-1,j)*dt(i-1,j)
     $                       *(v(i-1,j+1,k)+v(i-1,j,k)))
     $                 +grav*.125e0*(dt(i,j)+dt(i-1,j))
     $                   *(egf(i,j)-egf(i-1,j)+egb(i,j)-egb(i-1,j)
     $                     +(e_atmos(i,j)-e_atmos(i-1,j))*2.e0)
     $                   *(dy(i,j)+dy(i-1,j))
     $                 +drhox(i,j,k)
          end do
        end do
      end do

!  step forward in time
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            uf(i,j,k)=((h(i,j)+etb(i,j)+h(i-1,j)+etb(i-1,j))
     $                 *aru(i,j)*ub(i,j,k)
     $                 -2.e0*dti2*uf(i,j,k))
     $                /((h(i,j)+etf(i,j)+h(i-1,j)+etf(i-1,j))
     $                  *aru(i,j))
          end do
        end do
      end do

      return
      end

!_______________________________________________________________________
      subroutine advv
! do horizontal and vertical advection of v-momentum, and includes
! coriolis, surface slope and baroclinic terms
      implicit none
      include 'pom.h'
      integer i,j,k

! do vertical advection
      do k=1,kb
        do j=1,jm
          do i=1,im
            vf(i,j,k)=0.e0
          end do
        end do
      end do

      do k=2,kbm1
        do j=2,jm
          do i=1,im
            vf(i,j,k)=.25e0*(w(i,j,k)+w(i,j-1,k))
     $                     *(v(i,j,k)+v(i,j,k-1))
          end do
        end do
      end do

! combine horizontal and vertical advection with coriolis, surface
! slope and baroclinic terms
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            vf(i,j,k)=advy(i,j,k)
     $                 +(vf(i,j,k)-vf(i,j,k+1))*arv(i,j)/dz(k)
     $                 +arv(i,j)*.25e0
     $                   *(cor(i,j)*dt(i,j)
     $                      *(u(i+1,j,k)+u(i,j,k))
     $                     +cor(i,j-1)*dt(i,j-1)
     $                       *(u(i+1,j-1,k)+u(i,j-1,k)))
     $                 +grav*.125e0*(dt(i,j)+dt(i,j-1))
     $                   *(egf(i,j)-egf(i,j-1)+egb(i,j)-egb(i,j-1)
     $                     +(e_atmos(i,j)-e_atmos(i,j-1))*2.e0)
     $                   *(dx(i,j)+dx(i,j-1))
     $                 +drhoy(i,j,k)
          end do
        end do
      end do

! step forward in time
      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            vf(i,j,k)=((h(i,j)+etb(i,j)+h(i,j-1)+etb(i,j-1))
     $                 *arv(i,j)*vb(i,j,k)
     $                 -2.e0*dti2*vf(i,j,k))
     $                /((h(i,j)+etf(i,j)+h(i,j-1)+etf(i,j-1))
     $                  *arv(i,j))
          end do
        end do
      end do

      return
      end

!_______________________________________________________________________
      subroutine baropg
! calculate  baroclinic pressure gradient
      implicit none
      include 'pom.h'
      integer i,j,k

!      do k=1,kb
!        do j=1,jm
!          do i=1,im
!            rho(i,j,k)=rho(i,j,k)-rmean(i,j,k)
!          end do
!        end do
!      end do

      call exchange3d_mpi(rho(:,:,1:kbm1),im,jm,kbm1) ! 2018.07.28

! calculate x-component of baroclinic pressure gradient
      do j=2,jmm1
        do i=2,imm1
          drhox(i,j,1)=.5e0*grav*(-zz(1))*(dt(i,j)+dt(i-1,j))
     $                  *(rho(i,j,1)-rho(i-1,j,1))
        end do
      end do

      do k=2,kbm1
        do j=2,jmm1
          do i=2,imm1
            drhox(i,j,k)=drhox(i,j,k-1)
     $                    +grav*.25e0*(zz(k-1)-zz(k))
     $                      *(dt(i,j)+dt(i-1,j))
     $                      *(rho(i,j,k)-rho(i-1,j,k)
     $                        +rho(i,j,k-1)-rho(i-1,j,k-1))
     $                    +grav*.25e0*(zz(k-1)+zz(k))
     $                      *(dt(i,j)-dt(i-1,j))
     $                      *(rho(i,j,k)+rho(i-1,j,k)
     $                        -rho(i,j,k-1)-rho(i-1,j,k-1))
          end do
        end do
      end do

      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            drhox(i,j,k)=.25e0*(dt(i,j)+dt(i-1,j))
     $                        *drhox(i,j,k)*dum(i,j)
     $                        *(dy(i,j)+dy(i-1,j))
          end do
        end do
      end do

! calculate y-component of baroclinic pressure gradient
      do j=2,jmm1
        do i=2,imm1
          drhoy(i,j,1)=.5e0*grav*(-zz(1))*(dt(i,j)+dt(i,j-1))
     $                  *(rho(i,j,1)-rho(i,j-1,1))
        end do
      end do

      do k=2,kbm1
        do j=2,jmm1
          do i=2,imm1
            drhoy(i,j,k)=drhoy(i,j,k-1)
     $                    +grav*.25e0*(zz(k-1)-zz(k))
     $                      *(dt(i,j)+dt(i,j-1))
     $                      *(rho(i,j,k)-rho(i,j-1,k)
     $                        +rho(i,j,k-1)-rho(i,j-1,k-1))
     $                    +grav*.25e0*(zz(k-1)+zz(k))
     $                      *(dt(i,j)-dt(i,j-1))
     $                      *(rho(i,j,k)+rho(i,j-1,k)
     $                        -rho(i,j,k-1)-rho(i,j-1,k-1))
          end do
        end do
      end do

!      i = 2
!      j = 2
!      k = 2
!      write(6,*) '1 in baropg i j k drhoy ', i, j, k, drhoy(i,j,k)
!      write(6,*) '1 in baropg i j k rho ', rho(i,j,k), rho(i,j-1,k)
!      write(6,*) 'dt dvm dx ramp ',dt(i,j),dvm(i,j),dx(i,j),ramp

      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            drhoy(i,j,k)=.25e0*(dt(i,j)+dt(i,j-1))
     $                        *drhoy(i,j,k)*dvm(i,j)
     $                        *(dx(i,j)+dx(i,j-1))
          end do
        end do
      end do

      do k=1,kb
        do j=2,jmm1
          do i=2,imm1
            drhox(i,j,k)=ramp*drhox(i,j,k)
            drhoy(i,j,k)=ramp*drhoy(i,j,k)
          end do
        end do
      end do

!      do k=1,kb
!        do j=1,jm
!          do i=1,im
!            rho(i,j,k)=rho(i,j,k)+rmean(i,j,k)
!          end do
!        end do
!      end do

!      i = 2
!      j = 2
!      k = 2
!      write(6,*) '2 in baropg i j k drhoy ', i, j, k, drhoy(i,j,k)

      return
      end

!_______________________________________________________________________
      subroutine bcondrad(idx)
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

! eExternal (2-D) elevation boundary conditions
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
      subroutine baropg_mcc
! calculate  baroclinic pressure gradient
! 4th order correction terms, following McCalpin
      implicit none
      include 'pom.h'
      integer i,j,k
      real d4(im,jm),ddx(im,jm),drho(im,jm,kb),rhou(im,jm,kb)
      real rho4th(0:im,0:jm,kb),d4th(0:im,0:jm)

!      do k=1,kb
!        do j=1,jm
!          do i=1,im
!            rho(i,j,k)=rho(i,j,k)-rmean(i,j,k)
!          end do
!        end do
!      end do

      call exchange3d_mpi(rho(:,:,1:kbm1),im,jm,kbm1) ! 2018.07.28

! convert a 2nd order matrices to special 4th order
! special 4th order case
      call order2d_mpi(d,d4th,im,jm)
      call order3d_mpi(rho,rho4th,im,jm,kb)

! compute terms correct to 4th order
      do i=1,im
        do j=1,jm
          ddx(i,j)=0.
          d4(i,j)=0.
        end do
      end do
      do k=1,kb
        do j=1,jm
          do i=1,im
            rhou(i,j,k)=0.
            drho(i,j,k)=0.
          end do
        end do
      end do

! compute DRHO, RHOU, DDX and D4
      do j=1,jm
        do i=2,im
          do k=1,kbm1
            drho(i,j,k)=(rho(i,j,k)-rho(i-1,j,k))*dum(i,j)
            rhou(i,j,k)=0.5*(rho(i,j,k)+rho(i-1,j,k))*dum(i,j)
          end do
          ddx(i,j)=(d(i,j)-d(i-1,j))*dum(i,j)
          d4(i,j)=.5*(d(i,j)+d(i-1,j))*dum(i,j)
        end do
      end do

      if(n_west.eq.-1) then
        do j=1,jm
          do i=3,imm1
            do k=1,kbm1
              drho(i,j,k)=drho(i,j,k) - (1./24.)*
     $                    (dum(i+1,j)*(rho(i+1,j,k)-rho(i,j,k))-
     $                    2*(rho(i,j,k)-rho(i-1,j,k))+
     $                    dum(i-1,j)*(rho(i-1,j,k)-rho(i-2,j,k)))
              rhou(i,j,k)=rhou(i,j,k) + (1./16.)*
     $                    (dum(i+1,j)*(rho(i,j,k)-rho(i+1,j,k))+
     $                    dum(i-1,j)*(rho(i-1,j,k)-rho(i-2,j,k)))
            end do
            ddx(i,j)=ddx(i,j)-(1./24.)*
     $               (dum(i+1,j)*(d(i+1,j)-d(i,j))-
     $               2*(d(i,j)-d(i-1,j))+
     $               dum(i-1,j)*(d(i-1,j)-d(i-2,j)))
            d4(i,j)=d4(i,j)+(1./16.)*
     $              (dum(i+1,j)*(d(i,j)-d(i+1,j))+
     $              dum(i-1,j)*(d(i-1,j)-d(i-2,j)))
          end do
        end do
      else
        do j=1,jm
          do i=2,imm1
            do k=1,kbm1
              drho(i,j,k)=drho(i,j,k) - (1./24.)*
     $                   (dum(i+1,j)*(rho(i+1,j,k)-rho(i,j,k))-
     $                    2*(rho(i,j,k)-rho(i-1,j,k))+
     $                    dum(i-1,j)*(rho(i-1,j,k)-rho4th(i-2,j,k)))
              rhou(i,j,k)=rhou(i,j,k) + (1./16.)*
     $                    (dum(i+1,j)*(rho(i,j,k)-rho(i+1,j,k))+
     $                    dum(i-1,j)*(rho(i-1,j,k)-rho4th(i-2,j,k)))
            end do
            ddx(i,j)=ddx(i,j)-(1./24.)*
     $               (dum(i+1,j)*(d(i+1,j)-d(i,j))-
     $               2*(d(i,j)-d(i-1,j))+
     $               dum(i-1,j)*(d(i-1,j)-d4th(i-2,j)))
            d4(i,j)=d4(i,j)+(1./16.)*
     $              (dum(i+1,j)*(d(i,j)-d(i+1,j))+
     $              dum(i-1,j)*(d(i-1,j)-d4th(i-2,j)))
          end do
        end do
      end if

! calculate x-component of baroclinic pressure gradient
      do j=2,jmm1
        do i=2,imm1
          drhox(i,j,1)=grav*(-zz(1))*d4(i,j)*drho(i,j,1)
        end do
      end do

      do k=2,kbm1
        do j=2,jmm1
          do i=2,imm1
            drhox(i,j,k)=drhox(i,j,k-1)
     $                   +grav*0.5e0*dzz(k-1)*d4(i,j)
     $                   *(drho(i,j,k-1)+drho(i,j,k))
     $                   +grav*0.5e0*(zz(k-1)+zz(k))*ddx(i,j)
     $                   *(rhou(i,j,k)-rhou(i,j,k-1))
          end do
        end do
      end do

      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            drhox(i,j,k)=.25e0*(dt(i,j)+dt(i-1,j))
     $                        *drhox(i,j,k)*dum(i,j)
     $                        *(dy(i,j)+dy(i-1,j))
          end do
        end do
      end do

! compute terms correct to 4th order
      do i=1,im
        do j=1,jm
          ddx(i,j)=0.
          d4(i,j)=0.
        end do
      end do
      do k=1,kb
        do j=1,jm
          do i=1,im
            rhou(i,j,k)=0.
            drho(i,j,k)=0.
          end do
        end do
      end do

! compute DRHO, RHOU, DDX and D4
      do j=2,jm
        do i=1,im
          do k=1,kbm1
            drho(i,j,k)=(rho(i,j,k)-rho(i,j-1,k))*dvm(i,j)
            rhou(i,j,k)=.5*(rho(i,j,k)+rho(i,j-1,k))*dvm(i,j)
          end do
          ddx(i,j)=(d(i,j)-d(i,j-1))*dvm(i,j)
          d4(i,j)=.5*(d(i,j)+d(i,j-1))*dvm(i,j)
        end do
      end do

      if(n_south.eq.-1) then
        do j=3,jmm1
          do i=1,im
            do k=1,kbm1
              drho(i,j,k)=drho(i,j,k)-(1./24.)*
     $                    (dvm(i,j+1)*(rho(i,j+1,k)-rho(i,j,k))-
     $                    2*(rho(i,j,k)-rho(i,j-1,k))+
     $                    dvm(i,j-1)*(rho(i,j-1,k)-rho(i,j-2,k)))
              rhou(i,j,k)=rhou(i,j,k)+(1./16.)*
     $                    (dvm(i,j+1)*(rho(i,j,k)-rho(i,j+1,k))+
     $                    dvm(i,j-1)*(rho(i,j-1,k)-rho(i,j-2,k)))
            end do
            ddx(i,j)=ddx(i,j)-(1./24)*
     $               (dvm(i,j+1)*(d(i,j+1)-d(i,j))-
     $               2*(d(i,j)-d(i,j-1))+
     $               dvm(i,j-1)*(d(i,j-1)-d(i,j-2)))
            d4(i,j)=d4(i,j)+(1./16.)*
     $              (dvm(i,j+1)*(d(i,j)-d(i,j+1))+
     $              dvm(i,j-1)*(d(i,j-1)-d(i,j-2)))
          end do
        end do
      else
        do j=2,jmm1
          do i=1,im
            do k=1,kbm1
              drho(i,j,k)=drho(i,j,k)-(1./24.)*
     $                    (dvm(i,j+1)*(rho(i,j+1,k)-rho(i,j,k))-
     $                    2*(rho(i,j,k)-rho(i,j-1,k))+
     $                    dvm(i,j-1)*(rho(i,j-1,k)-rho4th(i,j-2,k)))
              rhou(i,j,k)=rhou(i,j,k)+(1./16.)*
     $                    (dvm(i,j+1)*(rho(i,j,k)-rho(i,j+1,k))+
     $                    dvm(i,j-1)*(rho(i,j-1,k)-rho4th(i,j-2,k)))
            end do
            ddx(i,j)=ddx(i,j)-(1./24)*
     $               (dvm(i,j+1)*(d(i,j+1)-d(i,j))-
     $               2*(d(i,j)-d(i,j-1))+
     $               dvm(i,j-1)*(d(i,j-1)-d4th(i,j-2)))
            d4(i,j)=d4(i,j)+(1./16.)*
     $              (dvm(i,j+1)*(d(i,j)-d(i,j+1))+
     $              dvm(i,j-1)*(d(i,j-1)-d4th(i,j-2)))
          end do
        end do
      end if

! calculate y-component of baroclinic pressure gradient
      do j=2,jmm1
        do i=2,imm1
          drhoy(i,j,1)=grav*(-zz(1))*d4(i,j)*drho(i,j,1)
        end do
      end do

      do k=2,kbm1
        do j=2,jmm1
          do i=2,imm1
            drhoy(i,j,k)=drhoy(i,j,k-1)
     $                   +grav*0.5e0*dzz(k-1)*d4(i,j)
     $                   *(drho(i,j,k-1)+drho(i,j,k))
     $                   +grav*0.5e0*(zz(k-1)+zz(k))*ddx(i,j)
     $                   *(rhou(i,j,k)-rhou(i,j,k-1))
          end do
        end do
      end do

      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            drhoy(i,j,k)=.25e0*(dt(i,j)+dt(i,j-1))
     $                        *drhoy(i,j,k)*dvm(i,j)
     $                        *(dx(i,j)+dx(i,j-1))
          end do
        end do
      end do

      do k=1,kb
        do j=2,jmm1
          do i=2,imm1
            drhox(i,j,k)=ramp*drhox(i,j,k)
            drhoy(i,j,k)=ramp*drhoy(i,j,k)
          end do
        end do
      end do

!      do k=1,kb
!        do j=1,jm
!          do i=1,im
!            rho(i,j,k)=rho(i,j,k)+rmean(i,j,k)
!          end do
!        end do
!      end do

      return
      end
!_______________________________________________________________________
      subroutine baropg_thiem

! calculate baroclinic pressure gradient (drhox,drhoy) 
! following Thiem and Berntsen (2006; Ocean Modelling, 12, 140-156)
!
!     2018.05 Created by S.Ohishi
!

      implicit none
      include 'pom.h'
      integer i,j,k

      real,parameter :: w1=3.e0,w2=1.e0 !Weightning coefficient(3-4:1)
!     real,parameter :: w1=1.e0,w2=0.e0 !Weightning coefficient(3-4:1)
      
      real ww2x(im,jm),ww2y(im,jm)
!     real del_rho

!     Subtract rmean for truncation error
!      do k=1,kb
!         do j=1,jm
!            do i=1,im
!               rho(i,j,k)=rho(i,j,k)-rmean(i,j,k)
!            end do
!         end do
!      end do

      call exchange3d_mpi(rho(:,:,1:kbm1),im,jm,kbm1) ! 2018.07.28

      do j=2,jmm1
         do i=2,imm1

            ww2x(i,j)=w2*fsm(i  ,j  )*fsm(i-1,j-1)
     &            *fsm(i-1,j  )*fsm(i  ,j-1)
     &            *fsm(i  ,j+1)
     &            *fsm(i-1,j+1)

            ww2y(i,j)=w2*fsm(i  ,j  )*fsm(i-1,j-1)
     &           *fsm(i-1,j  )*fsm(i  ,j-1)
     &           *fsm(i+1,j  )
     &           *fsm(i+1,j-1)

         end do
      end do

!     k=1
      do j=2,jmm1
         do i=2,imm1

            drhox(i,j,1)=
     &           grav*(-zz(1))
     &           *(w1*0.5e0*(dt(i,j)+dt(i-1,j)) !(i,j)-(i-1,j)
     &           *(rho(i,j,1)-rho(i-1,j,1))
     &           +ww2x(i,j)*0.125e0
     &           *((dt(i,j)+dt(i-1,j-1)) !(i,j)-(i-1,j-1)
     &           *(rho(i,j,1)-rho(i-1,j-1,1))
     &           -(dt(i-1,j)+dt(i,j-1)) !(i-1,j)-(i,j-1)
     &           *(rho(i-1,j,1)-rho(i,j-1,1))
     &           +(dt(i,j+1)+dt(i-1,j)) !(i,j+1)-(i-1,j)
     &           *(rho(i,j+1,1)-rho(i-1,j,1))
     &           -(dt(i-1,j+1)+dt(i,j)) !(i-1,j+1)-(i,j)
     &           *(rho(i-1,j+1,1)-rho(i,j,1))))
     &           /(w1+ww2x(i,j))

            drhoy(i,j,1)=
     &           grav*(-zz(1))
     &           *(w1*0.5e0*(dt(i,j)+dt(i,j-1)) !(i,j)-(i,j-1)
     &           *(rho(i,j,1)-rho(i,j-1,1))
     &           +ww2y(i,j)*0.125e0
     &           *((dt(i,j)+dt(i-1,j-1)) !(i,j)-(i-1,j-1)
     &           *(rho(i,j,1)-rho(i-1,j-1,1))
     &           +(dt(i-1,j)+dt(i,j-1)) !(i-1,j)-(i,j-1)
     &           *(rho(i-1,j,1)-rho(i,j-1,1))
     &           +(dt(i+1,j)+dt(i,j-1)) !(i+1,j)-(i,j-1)
     &           *(rho(i+1,j,1)-rho(i,j-1,1))
     &           +(dt(i,j)+dt(i+1,j-1)) !(i,j)-(i+1,j-1)
     &           *(rho(i,j,1)-rho(i+1,j-1,1))))
     &           /(w1+ww2y(i,j))

         end do
      end do

!     k=2-kbm1 integration
      do k=2,kbm1
         do j=2,jmm1
            do i=2,imm1
               
               drhox(i,j,k)=drhox(i,j,k-1)
     &              +0.25e0*grav*
     &              (w1*((zz(k-1)-zz(k))*(dt(i,j)+dt(i-1,j)) !(i,j)-(i-1,j)
     &              *(rho(i,j,k-1)-rho(i-1,j,k-1)
     &              +rho(i,j,k  )-rho(i-1,j,k  ))
     &              -(zz(k-1)+zz(k))*(dt(i,j)-dt(i-1,j))
     &              *(rho(i,j,k-1)+rho(i-1,j,k-1)
     &              -rho(i,j,k  )-rho(i-1,j,k  )))
     &              +ww2x(i,j)*0.25e0
     &              *(((zz(k-1)-zz(k))*(dt(i,j)+dt(i-1,j-1)) !(i,j)-(i-1,j-1)
     &              *(rho(i,j,k-1)-rho(i-1,j-1,k-1)
     &              +rho(i,j,k  )-rho(i-1,j-1,k  ))
     &              -(zz(k-1)+zz(k))*(dt(i,j)-dt(i-1,j-1))
     &              *(rho(i,j,k-1)+rho(i-1,j-1,k-1)
     &              -rho(i,j,k  )-rho(i-1,j-1,k  )))
     &              -((zz(k-1)-zz(k))*(dt(i-1,j)+dt(i,j-1)) !(i-1,j)-(i,j-1)
     &              *(rho(i-1,j,k-1)-rho(i,j-1,k-1)
     &              +rho(i-1,j,k  )-rho(i,j-1,k  ))
     &              -(zz(k-1)+zz(k))*(dt(i-1,j)-dt(i,j-1))
     &              *(rho(i-1,j,k-1)+rho(i,j-1,k-1)
     &              -rho(i-1,j,k  )-rho(i,j-1,k  )))
     &              +((zz(k-1)-zz(k))*(dt(i,j+1)+dt(i-1,j)) !(i,j+1)-(i-1,j)
     &              *(rho(i,j+1,k-1)-rho(i-1,j,k-1)
     &              +rho(i,j+1,k  )-rho(i-1,j,k  ))
     &              -(zz(k-1)+zz(k))*(dt(i,j+1)-dt(i-1,j))
     &              *(rho(i,j+1,k-1)+rho(i-1,j,k-1)
     &              -rho(i,j+1,k  )-rho(i-1,j,k  )))
     &              -((zz(k-1)-zz(k))*(dt(i-1,j+1)+dt(i,j)) !(i-1,j+1)-(i,j)
     &              *(rho(i-1,j+1,k-1)-rho(i,j,k-1)
     &              +rho(i-1,j+1,k  )-rho(i,j,k  ))
     &              -(zz(k-1)+zz(k))*(dt(i-1,j+1)-dt(i,j))
     &              *(rho(i-1,j+1,k-1)+rho(i,j,k-1)
     &              -rho(i-1,j+1,k  )-rho(i,j,k  )))))
     &              /(w1+ww2x(i,j))

               drhoy(i,j,k)=drhoy(i,j,k-1)
     &              +0.25e0*grav*
     &              (w1*((zz(k-1)-zz(k))*(dt(i,j)+dt(i,j-1)) !(i,j)-(i,j-1)
     &              *(rho(i,j,k-1)-rho(i,j-1,k-1)
     &              +rho(i,j,k  )-rho(i,j-1,k  ))
     &              -(zz(k-1)+zz(k))*(dt(i,j)-dt(i,j-1))
     &              *(rho(i,j,k-1)+rho(i,j-1,k-1)
     &              -rho(i,j,k  )-rho(i,j-1,k  )))
     &              +ww2y(i,j)*0.25e0
     &              *(((zz(k-1)-zz(k))*(dt(i,j)+dt(i-1,j-1)) !(i,j)-(i-1,j-1)
     &              *(rho(i,j,k-1)-rho(i-1,j-1,k-1)
     &              +rho(i,j,k  )-rho(i-1,j-1,k  ))
     &              -(zz(k-1)+zz(k))*(dt(i,j)-dt(i-1,j-1))
     &              *(rho(i,j,k-1)+rho(i-1,j-1,k-1)
     &              -rho(i,j,k  )-rho(i-1,j-1,k  )))
     &              +((zz(k-1)-zz(k))*(dt(i-1,j)+dt(i,j-1)) !(i-1,j)-(i,j-1)
     &              *(rho(i-1,j,k-1)-rho(i,j-1,k-1)
     &              +rho(i-1,j,k  )-rho(i,j-1,k  ))
     &              -(zz(k-1)+zz(k))*(dt(i-1,j)-dt(i,j-1))
     &              *(rho(i-1,j,k-1)+rho(i,j-1,k-1)
     &              -rho(i-1,j,k  )-rho(i,j-1,k  )))
     &              +((zz(k-1)-zz(k))*(dt(i+1,j)+dt(i,j-1)) !(i+1,j)-(i,j-1)
     &              *(rho(i+1,j,k-1)-rho(i,j-1,k-1)
     &              +rho(i+1,j,k  )-rho(i,j-1,k  ))
     &              -(zz(k-1)+zz(k))*(dt(i+1,j)-dt(i,j-1))
     &              *(rho(i+1,j,k-1)+rho(i,j-1,k-1)
     &              -rho(i+1,j,k  )-rho(i,j-1,k  )))
     &              +((zz(k-1)-zz(k))*(dt(i,j)+dt(i+1,j-1)) !(i,j)-(i+1,j-1)
     &              *(rho(i,j,k-1)-rho(i+1,j-1,k-1)
     &              +rho(i,j,k  )-rho(i+1,j-1,k  ))
     &              -(zz(k-1)+zz(k))*(dt(i,j)-dt(i+1,j-1))
     &              *(rho(i,j,k-1)+rho(i+1,j-1,k-1)
     &              -rho(i,j,k  )-rho(i+1,j-1,k  )))))
     &              /(w1+ww2y(i,j))

            end do
         end do
      end do

      do k=1,kbm1
         do j=2,jmm1
            do i=2,imm1
               
               drhox(i,j,k)=dum(i,j)*0.25e0*(dt(i,j)+dt(i-1,j))
     &              *drhox(i,j,k)*(dy(i,j)+dy(i-1,j))*ramp

               drhoy(i,j,k)=dvm(i,j)*0.25e0*(dt(i,j)+dt(i,j-1))
     &              *drhoy(i,j,k)*(dx(i,j)+dx(i,j-1))*ramp

            end do
         end do
      end do

!     Add rmean
!      do k=1,kb
!         do j=1,jm
!            do i=1,im
!               rho(i,j,k)=rho(i,j,k)+rmean(i,j,k)
!            end do
!         end do
!      end do

!     debug
!      write(6,*) "drhox in baropg_thiem:",
!     &     dum(2,2),east_u(2,2),north_u(2,2),drhox(2,2,2)
!      write(6,*) "drhoy in baropg_thiem:",
!     &     dvm(2,2),east_v(2,2),north_v(2,2),drhoy(2,2,2)

      return
      end

!---

      function del_rho(i1,j1,i2,j2,k)

!     Not used in baropg_thiem
      implicit none
      include 'pom.h'

      integer i1,j1,i2,j2,k
      real del_rho

      if(k==1)then
         del_rho=(dt(i1,j1)+dt(i2,j2))
     &        *(rho(i1,j1,1)-rho(i2,j2,1))
      else
         del_rho=
     &        (zz(k-1)-zz(k))*(dt(i1,j1)+dt(i2,j2))
     &        *(rho(i1,j1,k-1)-rho(i2,j2,k-1)
     &         +rho(i1,j1,k  )-rho(i2,j2,k  ))
     &        -(zz(k-1)+zz(k))*(dt(i1,j1)-dt(i2,j2))
     &        *(rho(i1,j1,k-1)+rho(i2,j2,k-1)
     &         -rho(i1,j1,k  )-rho(i2,j2,k  ))

      end if

      return
      end function

!_______________________________________________________________________
      subroutine dens(si,ti,rhoo)
! calculate (density-1000.)/rhoref.
! see: Mellor, G.L., 1991, J. Atmos. Oceanic Tech., 609-611
! note: if pressure is not used in dens, buoyancy term (boygr) in
! subroutine profq must be changed (see note in subroutine profq)
      implicit none
      include 'pom.h'
      real si(im,jm,kb),ti(im,jm,kb),rhoo(im,jm,kb)
      integer i,j,k
      real cr,p,rhor,sr,tr,tr2,tr3,tr4

      do k=1,kbm1
        do j=1,jm
          do i=1,im

            tr=ti(i,j,k)+tbias
            sr=si(i,j,k)+sbias
            tr2=tr*tr
            tr3=tr2*tr
            tr4=tr3*tr

! approximate pressure in units of bars
            p=grav*rhoref*(-zz(k)* h(i,j))*1.e-5

            rhor=-0.157406e0+6.793952e-2*tr
     $           -9.095290e-3*tr2+1.001685e-4*tr3
     $           -1.120083e-6*tr4+6.536332e-9*tr4*tr

            rhor=rhor+(0.824493e0-4.0899e-3*tr
     $               +7.6438e-5*tr2-8.2467e-7*tr3
     $               +5.3875e-9*tr4)*sr
     $               +(-5.72466e-3+1.0227e-4*tr
     $               -1.6546e-6*tr2)*abs(sr)**1.5
     $               +4.8314e-4*sr*sr

            cr=1449.1e0+.0821e0*p+4.55e0*tr-.045e0*tr2
     $                 +1.34e0*(sr-35.e0)
            rhor=rhor+1.e5*p/(cr*cr)*(1.e0-2.e0*p/(cr*cr))

            rhoo(i,j,k)=rhor/rhoref*fsm(i,j)

          end do
        end do
      end do

      return
      end

! kii2b
!_______________________________________________________________________
      subroutine profq
! solve for q2 (twice the turbulent kinetic energy), q2l (q2 x turbulent
! length scale), km (vertical kinematic viscosity) and kh (vertical
! kinematic diffusivity), using a simplified version of the level 2 1/2
! model of Mellor and Yamada (1982)
! in this version, the Craig-Banner sub-model whereby breaking wave tke
! is injected into the surface is included. However, we use an
! analytical solution to the near surface tke equation to solve for q2
! at the surface giving the same result as C-B diffusion. The new scheme
! is simpler and more robust than the latter scheme
      implicit none
      include 'pom.h'
      real a(im,jm,kb),c(im,jm,kb)
      real ee(im,jm,kb),gg(im,jm,kb)
      real sm(im,jm,kb),sh(im,jm,kb)
      real cc(im,jm,kb)
      real gh(im,jm,kb),boygr(im,jm,kb),dh(im,jm),stf(im,jm,kb)
      real prod(im,jm,kb)
      real a1,a2,b1,b2,c1
      real coef1,coef2,coef3,coef4,coef5
      real const1,e1,e2,ghc
      real p,sef,sp,tp
      real l0(im,jm)
      real cbcnst,surfl,shiw
      real utau2(im,jm)
      real df0,df1,df2
      integer i,j,k,ki

      data a1,b1,a2,b2,c1/0.92e0,16.6e0,0.74e0,10.1e0,0.08e0/
      data e1/1.8e0/,e2/1.33e0/
      data sef/1.e0/
      data cbcnst/100./surfl/2.e5/shiw/0.0/

!      write(6,*) 'in profq '
!      call flush(6)

      do j=1,jm
        do i=1,im
          dh(i,j)=h(i,j)+etf(i,j)
        end do
      end do

      do k=1,kb
        do j=1,jm
          do i=1,im
            a(i,j,k)=0.e0
            c(i,j,k)=0.e0
            ee(i,j,k)=0.e0
            gg(i,j,k)=0.e0
          end do
        end do
      end do

      do k=2,kbm1
        do j=1,jm
          do i=1,im
            a(i,j,k)=-dti2*(kq(i,j,k+1)+kq(i,j,k)+2.e0*umol)*.5e0
     $                /(dzz(k-1)*dz(k)*dh(i,j)*dh(i,j))
            c(i,j,k)=-dti2*(kq(i,j,k-1)+kq(i,j,k)+2.e0*umol)*.5e0
     $                /(dzz(k-1)*dz(k-1)*dh(i,j)*dh(i,j))
          end do
        end do
      end do

!      write(6,*) 'in profq 0'
!      call flush(6)

! the following section solves the equation:
!     dti2*(kq*q2')' - q2*(2.*dti2*dtef+1.) = -q2b

! surface and bottom boundary conditions
      const1=(16.6e0**(2.e0/3.e0))*sef

! initialize fields that are not calculated on all boundaries
! but are later used there
      do i=1,im
        do j=1,jm
          l0(i,j)=0.
         end do
      end do
      do i=1,im
        do j=1,jm
          do k=1,kb
            boygr(i,j,k)=0.
            prod(i,j,k)=0.
          end do
        end do
      end do

      do j=1,jmm1
        do i=1,imm1
          utau2(i,j)=sqrt((.5e0*(wusurf(i,j)+wusurf(i+1,j)))**2
     $                   +(.5e0*(wvsurf(i,j)+wvsurf(i,j+1)))**2)
          uf(i,j,kb)=sqrt((.5e0*(wubot(i,j)+wubot(i+1,j)))**2
     $                   +(.5e0*(wvbot(i,j)+wvbot(i,j+1)))**2)*const1
        end do
      end do
      call exchange2d_mpi(utau2,im,jm)
      call exchange2d_mpi(uf(:,:,kb),im,jm)

!      write(6,*) 'in profq a'
!      call flush(6)

      do j=1,jm
        do i=1,im
! wave breaking energy- a variant of Craig & Banner (1994)
! see Mellor and Blumberg, 2003.
          ee(i,j,1)=0.e0
          gg(i,j,1)=(15.8*cbcnst)**(2./3.)*utau2(i,j)
! surface length scale following Stacey (1999).
          l0(i,j)=surfl*utau2(i,j)/grav
        end do
      end do

! calculate speed of sound squared
      do k=1,kb
        do j=1,jm
          do i=1,im
            cc(i,j,k)=0.
          end do
        end do
      end do
      do k=1,kbm1
        do j=1,jm
          do i=1,im
            tp=t(i,j,k)+tbias
            sp=s(i,j,k)+sbias
! calculate pressure in units of decibars
            p=grav*rhoref*(-zz(k)*h(i,j))*1.e-4
            cc(i,j,k)=1449.1e0+.00821e0*p+4.55e0*tp-.045e0*tp**2
     $                 +1.34e0*(sp-35.0e0)
            cc(i,j,k)=cc(i,j,k)
     $                 /sqrt((1.e0-.01642e0*p/cc(i,j,k))
     $                   *(1.e0-0.40e0*p/cc(i,j,k)**2))
          end do
        end do
      end do

! calculate buoyancy gradient
      do k=2,kbm1
        do j=1,jm
          do i=1,im
            q2b(i,j,k)=abs(q2b(i,j,k))
            q2lb(i,j,k)=abs(q2lb(i,j,k))
            boygr(i,j,k)=grav*(rho(i,j,k-1)-rho(i,j,k))
     $                    /(dzz(k-1)*h(i,j))
! *** note: comment out next line if dens does not include pressure
     $      +(grav**2)*2.e0/(cc(i,j,k-1)**2+cc(i,j,k)**2)
          end do
        end do
      end do

      do k=2,kbm1
        do j=1,jm
          do i=1,im
            l(i,j,k)=abs(q2lb(i,j,k)/q2b(i,j,k))
            if(z(k).gt.-0.5) l(i,j,k)=max(l(i,j,k),kappa*l0(i,j))
            gh(i,j,k)=(l(i,j,k)**2)*boygr(i,j,k)/q2b(i,j,k)
            gh(i,j,k)=min(gh(i,j,k),.028e0)
          end do
        end do
      end do

      do j=1,jm
        do i=1,im
          l(i,j,1)=kappa*l0(i,j)
          l(i,j,kb)=0.e0
          gh(i,j,1)=0.e0
          gh(i,j,kb)=0.e0
        end do
      end do

! calculate production of turbulent kinetic energy:
      do k=2,kbm1
        do j=2,jmm1
          do i=2,imm1
            prod(i,j,k)=km(i,j,k)*.25e0*sef
     $                   *((u(i,j,k)-u(i,j,k-1)
     $                      +u(i+1,j,k)-u(i+1,j,k-1))**2
     $                     +(v(i,j,k)-v(i,j,k-1)
     $                      +v(i,j+1,k)-v(i,j+1,k-1))**2)
     $                   /(dzz(k-1)*dh(i,j))**2
! add shear due to internal wave field
     $             -shiw*km(i,j,k)*boygr(i,j,k)
            prod(i,j,k)=prod(i,j,k)+kh(i,j,k)*boygr(i,j,k)
          end do
        end do
      end do
      call exchange3d_mpi(prod(:,:,2:kbm1),im,jm,kbm2)

!      write(6,*) 'in profq b'
!      call flush(6)

! note: Richardson # dep. dissipation correction (Mellor, 2001; Ezer,
! 2000), depends on ghc the critical number (empirical -6 to -2) to
! increase mixing
      ghc=-6.0e0
      do k=1,kb
        do j=1,jm
          do i=1,im
            stf(i,j,k)=1.e0
! It is unclear yet if diss. corr. is needed when surf. waves are included.
!           if(gh(i,j,k).lt.0.e0)
!    $        stf(i,j,k)=1.0e0-0.9e0*(gh(i,j,k)/ghc)**1.5e0
!           if(gh(i,j,k).lt.ghc) stf(i,j,k)=0.1e0
            dtef(i,j,k)=sqrt(abs(q2b(i,j,k)))*stf(i,j,k)
     $                   /(b1*l(i,j,k)+small)
          end do
        end do
      end do

      do k=2,kbm1
        do j=1,jm
          do i=1,im
            gg(i,j,k)=1.e0/(a(i,j,k)+c(i,j,k)*(1.e0-ee(i,j,k-1))
     $                      -(2.e0*dti2*dtef(i,j,k)+1.e0))
            ee(i,j,k)=a(i,j,k)*gg(i,j,k)
            gg(i,j,k)=(-2.e0*dti2*prod(i,j,k)+c(i,j,k)*gg(i,j,k-1)
     $                 -uf(i,j,k))*gg(i,j,k)
          end do
        end do
      end do

      do k=1,kbm1
        ki=kb-k
        do j=1,jm
          do i=1,im
            uf(i,j,ki)=ee(i,j,ki)*uf(i,j,ki+1)+gg(i,j,ki)
          end do
        end do
      end do

! the following section solves the equation:
!     dti2(kq*q2l')' - q2l*(dti2*dtef+1.) = -q2lb
      do j=1,jm
        do i=1,im
          vf(i,j,1)=0.
          vf(i,j,kb)=0.
          ee(i,j,2)=0.e0
          gg(i,j,2)=-kappa*z(2)*dh(i,j)*q2(i,j,2)
          vf(i,j,kb-1)=kappa*(1+z(kbm1))*dh(i,j)*q2(i,j,kbm1)
        end do
      end do
      do k=2,kbm1
        do j=1,jm
          do i=1,im
            dtef(i,j,k)=dtef(i,j,k)
     $                   *(1.e0+e2*((1.e0/abs(z(k)-z(1))
     $                               +1.e0/abs(z(k)-z(kb)))
     $                                *l(i,j,k)/(dh(i,j)*kappa))**2)
          end do
        end do
      end do
      do k=3,kbm1
        do j=1,jm
          do i=1,im
            gg(i,j,k)=1.e0/(a(i,j,k)+c(i,j,k)*(1.e0-ee(i,j,k-1))
     $                      -(dti2*dtef(i,j,k)+1.e0))
            ee(i,j,k)=a(i,j,k)*gg(i,j,k)
            gg(i,j,k)=(dti2*(-prod(i,j,k)*l(i,j,k)*e1)
     $                 +c(i,j,k)*gg(i,j,k-1)-vf(i,j,k))*gg(i,j,k)
          end do
        end do
      end do

      do k=1,kb-2
        ki=kb-k
        do j=1,jm
          do i=1,im
            vf(i,j,ki)=ee(i,j,ki)*vf(i,j,ki+1)+gg(i,j,ki)
          end do
        end do
      end do
! the following is to counter the problem of the ratio of two small
! numbers (l = q2l/q2) or one number becoming negative. Two options are
! included below. In this application, the second option, l was less
! noisy when uf or vf is small
      do k=2,kbm1
        do j=1,jm
          do i=1,im
!           if(uf(i,j,k).le.small.or.vf(i,j,k).le.small) then
!             uf(i,j,k)=small
!             vf(i,j,k)=0.1*dt(i,j)*small
!           end if
          uf(i,j,k)=abs(uf(i,j,k))
          vf(i,j,k)=abs(vf(i,j,k))
          end do
        end do
      end do

! the following section solves for km and kh
      coef4=18.e0*a1*a1+9.e0*a1*a2
      coef5=9.e0*a1*a2

! note that sm and sh limit to infinity when gh approaches 0.0288
      do k=1,kb
        do j=1,jm
          do i=1,im
            coef1=a2*(1.e0-6.e0*a1/b1*stf(i,j,k))
            coef2=3.e0*a2*b2/stf(i,j,k)+18.e0*a1*a2
            coef3=a1*(1.e0-3.e0*c1-6.e0*a1/b1*stf(i,j,k))
            sh(i,j,k)=coef1/(1.e0-coef2*gh(i,j,k))
            sm(i,j,k)=coef3+sh(i,j,k)*coef4*gh(i,j,k)
            sm(i,j,k)=sm(i,j,k)/(1.e0-coef5*gh(i,j,k))
          end do
        end do
      end do

! there are 2 options for kq which, unlike km and kh, was not derived by
! Mellor and Yamada but was purely empirical based on neutral boundary
! layer data. The choice is whether or not it should be subject to the
! stability factor, sh. Generally, there is not a great difference in
! output
      do k=1,kb
        do j=1,jm
          do i=1,im
            prod(i,j,k)=l(i,j,k)*sqrt(abs(q2(i,j,k)))
            kq(i,j,k)=(prod(i,j,k)*.41e0*sh(i,j,k)+kq(i,j,k))*.5e0
!            kq(i,j,k)=(prod(i,j,k)*.20+kq(i,j,k))*.5e0
            km(i,j,k)=(prod(i,j,k)*sm(i,j,k)+km(i,j,k))*.5e0
            kh(i,j,k)=(prod(i,j,k)*sh(i,j,k)+kh(i,j,k))*.5e0
          end do
        end do
      end do

!      write(6,*) 'in profq c'
!      call flush(6)

! cosmetics: make boundr. values as interior (even if not used, printout
! may show strange values)
      do k=1,kb
        do i=1,im
          if(n_north.eq.-1) then
            km(i,jm,k)=km(i,jmm1,k)
            kh(i,jm,k)=kh(i,jmm1,k)
            kq(i,jm,k)=kq(i,jmm1,k)
          end if
          if(n_south.eq.-1) then
            km(i,1,k)=km(i,2,k)
            kh(i,1,k)=kh(i,2,k)
            kq(i,1,k)=kq(i,2,k)
          end if
        end do
        do j=1,jm
          if(n_east.eq.-1) then
            km(im,j,k)=km(imm1,j,k)
            kh(im,j,k)=kh(imm1,j,k)
            kq(im,j,k)=kq(imm1,j,k)
          end if
          if(n_west.eq.-1) then
            km(1,j,k)=km(2,j,k)
            kh(1,j,k)=kh(2,j,k)
            kq(1,j,k)=kq(2,j,k)
          end if
        end do
      end do

      do k=1,kb
        do i=1,im
          do j=1,jm
            km(i,j,k)=km(i,j,k)*fsm(i,j)
            kh(i,j,k)=kh(i,j,k)*fsm(i,j)
            kq(i,j,k)=kq(i,j,k)*fsm(i,j)
          end do
        end do
      end do

      call exchange3d_mpi(km,im,jm,kb) ! 2018.07
      call exchange3d_mpi(kh,im,jm,kb) ! 2018.07
      call exchange3d_mpi(kq,im,jm,kb) ! 2018.08
      call exchange3d_mpi(l,im,jm,kb) !2018.08

      return
      end

!_______________________________________________________________________
      subroutine proft(f,wfsurf,fsurf,nbc)
! solves for vertical diffusion of temperature and salinity using method
! described by Richmeyer and Morton (1967)
! note: wfsurf and swrad are negative values when water column is
! warming or salt is being added
      implicit none
      include 'pom.h'
      real f(im,jm,kb),wfsurf(im,jm)
      real fsurf(im,jm),dh(im,jm)
      integer nbc
      real*8 a(im,jm,kb),c(im,jm,kb)
      real*8 ee(im,jm,kb),gg(im,jm,kb)
      real rad(im,jm,kb),r(5),ad1(5),ad2(5)
      integer i,j,k,ki

      real fb(im,jm,kb) !for conservation (budget)

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

!     Conservation
      if(budget == 1)then
         do k=1,kb
            do j=1,jm
               do i=1,im
                  fb(i,j,k)=f(i,j,k)
                  radterm(i,j,k)=0.e0
               end do
            end do
         end do
      end if

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
!        do k=1,kbm1
! rad(kb) passes into the bottom without
! heating lower layer or is reflected from bottom without been adsorbed.
! 2013.10.05 Miyazawa
        do k=1,kb
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

!     Conservation
      if(budget == 1)then
         
         do j=1,jm
            do i=1,im
               sfcterm(i,j)=-wfsurf(i,j)/(dz(1)*dh(i,j))
            end do
         end do

         do k=1,kbm1
            do j=1,jm
               do i=1,im
                  radterm(i,j,k)=-(rad(i,j,k)-rad(i,j,k+1))
     $                 /(dh(i,j)*dz(k))
                  zdifterm(i,j,k)=(f(i,j,k)-fb(i,j,k))/dti2
     $                 -radterm(i,j,k)
               end do
            end do
         end do

         do j=1,jm
            do i=1,im
               zdifterm(i,j,1)=zdifterm(i,j,1)-sfcterm(i,j)
            end do
         end do

      end if

      return
      end

!_______________________________________________________________________
      subroutine profu
! solves for vertical diffusion of x-momentum using method described by
! Richmeyer and Morton (1967)
! note: wusurf has the opposite sign to the wind speed
      implicit none
      include 'pom.h'
      real*8 a(im,jm,kb),c(im,jm,kb)
      real*8 ee(im,jm,kb),gg(im,jm,kb)
      real dh(im,jm)
      integer i,j,k,ki

! the following section solves the equation
!   dti2*(km*u')'-u=-ub
      do j=1,jm
        do i=1,im
          dh(i,j)=1.e0
        end do
      end do

      do j=2,jm
        do i=2,im
          dh(i,j)=(h(i,j)+etf(i,j)+h(i-1,j)+etf(i-1,j))*.5e0
        end do
      end do

      do k=1,kb
        do j=2,jm
          do i=2,im
            c(i,j,k)=(km(i,j,k)+km(i-1,j,k))*.5e0
          end do
        end do
      end do

      do k=2,kbm1
        do j=1,jm
          do i=1,im
            a(i,j,k-1)=-dti2*(c(i,j,k)+umol)
     $                  /(dz(k-1)*dzz(k-1)*dh(i,j)*dh(i,j))
            c(i,j,k)=-dti2*(c(i,j,k)+umol)
     $                /(dz(k)*dzz(k-1)*dh(i,j)*dh(i,j))
          end do
        end do
      end do

      do j=1,jm
        do i=1,im
          ee(i,j,1)=a(i,j,1)/(a(i,j,1)-1.e0)
          gg(i,j,1)=(-dti2*wusurf(i,j)/(-dz(1)*dh(i,j))
     $               -uf(i,j,1))
     $               /(a(i,j,1)-1.e0)
        end do
      end do

      do k=2,kbm2
        do j=1,jm
          do i=1,im
            gg(i,j,k)=1.e0/(a(i,j,k)+c(i,j,k)*(1.e0-ee(i,j,k-1))-1.e0)
            ee(i,j,k)=a(i,j,k)*gg(i,j,k)
            gg(i,j,k)=(c(i,j,k)*gg(i,j,k-1)-uf(i,j,k))*gg(i,j,k)
          end do
        end do
      end do

      do j=2,jmm1
        do i=2,imm1
          tps(i,j)=0.5e0*(cbc(i,j)+cbc(i-1,j))
     $              *sqrt(ub(i,j,kbm1)**2
     $                +(.25e0*(vb(i,j,kbm1)+vb(i,j+1,kbm1)
     $                         +vb(i-1,j,kbm1)+vb(i-1,j+1,kbm1)))**2)
          uf(i,j,kbm1)=(c(i,j,kbm1)*gg(i,j,kbm2)-uf(i,j,kbm1))
     $                  /(tps(i,j)*dti2/(-dz(kbm1)*dh(i,j))-1.e0
     $                    -(ee(i,j,kbm2)-1.e0)*c(i,j,kbm1))
          uf(i,j,kbm1)=uf(i,j,kbm1)*dum(i,j)
        end do
      end do

      do k=2,kbm1
        ki=kb-k
        do j=2,jmm1
          do i=2,imm1
            uf(i,j,ki)=(ee(i,j,ki)*uf(i,j,ki+1)+gg(i,j,ki))*dum(i,j)
          end do
        end do
      end do

      do j=2,jmm1
        do i=2,imm1
          wubot(i,j)=-tps(i,j)*uf(i,j,kbm1)
        end do
      end do
      call exchange2d_mpi(wubot,im,jm)

      return
      end

!_______________________________________________________________________
      subroutine profv
! solves for vertical diffusion of x-momentum using method described by
! Richmeyer and Morton (1967)
! note: wvsurf has the opposite sign to the wind speed
      implicit none
      include 'pom.h'
      real*8 a(im,jm,kb),c(im,jm,kb)
      real*8 ee(im,jm,kb),gg(im,jm,kb)
      real dh(im,jm)
      integer i,j,k,ki

! the following section solves the equation
!     dti2*(km*u')'-u=-ub

      do j=1,jm
        do i=1,im
          dh(i,j)=1.e0
        end do
      end do

      do j=2,jm
        do i=2,im
          dh(i,j)=.5e0*(h(i,j)+etf(i,j)+h(i,j-1)+etf(i,j-1))
        end do
      end do

      do k=1,kb
        do j=2,jm
          do i=2,im
            c(i,j,k)=(km(i,j,k)+km(i,j-1,k))*.5e0
          end do
        end do
      end do

      do k=2,kbm1
        do j=1,jm
          do i=1,im
            a(i,j,k-1)=-dti2*(c(i,j,k)+umol)
     $                  /(dz(k-1)*dzz(k-1)*dh(i,j)*dh(i,j))
            c(i,j,k)=-dti2*(c(i,j,k)+umol)
     $                /(dz(k)*dzz(k-1)*dh(i,j)*dh(i,j))
          end do
        end do
      end do

      do j=1,jm
        do i=1,im
          ee(i,j,1)=a(i,j,1)/(a(i,j,1)-1.e0)
          gg(i,j,1)=(-dti2*wvsurf(i,j)/(-dz(1)*dh(i,j))-vf(i,j,1))
     $               /(a(i,j,1)-1.e0)
        end do
      end do

      do k=2,kbm2
        do j=1,jm
          do i=1,im
            gg(i,j,k)=1.e0/(a(i,j,k)+c(i,j,k)*(1.e0-ee(i,j,k-1))-1.e0)
            ee(i,j,k)=a(i,j,k)*gg(i,j,k)
            gg(i,j,k)=(c(i,j,k)*gg(i,j,k-1)-vf(i,j,k))*gg(i,j,k)
          end do
        end do
      end do

      do j=2,jmm1
        do i=2,imm1
          tps(i,j)=0.5e0*(cbc(i,j)+cbc(i,j-1))
     $              *sqrt((.25e0*(ub(i,j,kbm1)+ub(i+1,j,kbm1)
     $                            +ub(i,j-1,kbm1)+ub(i+1,j-1,kbm1)))**2
     $                    +vb(i,j,kbm1)**2)
          vf(i,j,kbm1)=(c(i,j,kbm1)*gg(i,j,kbm2)-vf(i,j,kbm1))
     $                  /(tps(i,j)*dti2/(-dz(kbm1)*dh(i,j))-1.e0
     $                    -(ee(i,j,kbm2)-1.e0)*c(i,j,kbm1))
          vf(i,j,kbm1)=vf(i,j,kbm1)*dvm(i,j)
        end do
      end do

      do k=2,kbm1
        ki=kb-k
        do j=2,jmm1
          do i=2,imm1
            vf(i,j,ki)=(ee(i,j,ki)*vf(i,j,ki+1)+gg(i,j,ki))*dvm(i,j)
          end do
        end do
      end do

      do j=2,jmm1
        do i=2,imm1
          wvbot(i,j)=-tps(i,j)*vf(i,j,kbm1)
        end do
      end do
      call exchange2d_mpi(wvbot,im,jm)

      return
      end

!_______________________________________________________________________
      subroutine smol_adif(xmassflux,ymassflux,zwflux,ff)
! calculate the antidiffusive velocity used to reduce the numerical
! diffusion associated with the upstream differencing scheme
! this is based on a subroutine of Gianmaria Sannino (Inter-university
! Computing Consortium, Rome, Italy) and Vincenzo Artale (Italian
! National Agency for New Technology and Environment, Rome, Italy)
      implicit none
      include 'pom.h'
      real ff(im,jm,kb)
      real xmassflux(im,jm,kb),ymassflux(im,jm,kb),zwflux(im,jm,kb)
      real mol,abs_1,abs_2
      real value_min,epsilon
      real udx,u2dt,vdy,v2dt,wdz,w2dt
      integer i,j,k
      parameter (value_min=1.e-9,epsilon=1.0e-14)

! apply temperature and salinity mask
      do k=1,kb
        do i=1,im
          do j=1,jm
            ff(i,j,k)=ff(i,j,k)*fsm(i,j)
          end do
        end do
      end do

! recalculate mass fluxes with antidiffusion velocity
      do k=1,kbm1
        do j=2,jmm1
          do i=2,im
            if(ff(i,j,k).lt.value_min.or.
     $         ff(i-1,j,k).lt.value_min) then
              xmassflux(i,j,k)=0.e0
            else
              udx=abs(xmassflux(i,j,k))
              u2dt=dti2*xmassflux(i,j,k)*xmassflux(i,j,k)*2.e0
     $              /(aru(i,j)*(dt(i-1,j)+dt(i,j)))
              mol=(ff(i,j,k)-ff(i-1,j,k))
     $             /(ff(i-1,j,k)+ff(i,j,k)+epsilon)
              xmassflux(i,j,k)=(udx-u2dt)*mol*sw
              abs_1=abs(udx)
              abs_2=abs(u2dt)
              if(abs_1.lt.abs_2) xmassflux(i,j,k)=0.e0
            end if
          end do
        end do
      end do

      do k=1,kbm1
        do j=2,jm
          do i=2,imm1
            if(ff(i,j,k).lt.value_min.or.
     $         ff(i,j-1,k).lt.value_min) then
              ymassflux(i,j,k)=0.e0
            else
             vdy=abs(ymassflux(i,j,k))
             v2dt=dti2*ymassflux(i,j,k)*ymassflux(i,j,k)*2.e0
     $             /(arv(i,j)*(dt(i,j-1)+dt(i,j)))
             mol=(ff(i,j,k)-ff(i,j-1,k))
     $            /(ff(i,j-1,k)+ff(i,j,k)+epsilon)
             ymassflux(i,j,k)=(vdy-v2dt)*mol*sw
             abs_1=abs(vdy)
             abs_2=abs(v2dt)
             if(abs_1.lt.abs_2) ymassflux(i,j,k)=0.e0
            end if
          end do
        end do
      end do

      do k=2,kbm1
        do j=2,jmm1
          do i=2,imm1
            if(ff(i,j,k).lt.value_min.or.
     $         ff(i,j,k-1).lt.value_min) then
              zwflux(i,j,k)=0.e0
            else
              wdz=abs(zwflux(i,j,k))
              w2dt=dti2*zwflux(i,j,k)*zwflux(i,j,k)/(dzz(k-1)*dt(i,j))
              mol=(ff(i,j,k-1)-ff(i,j,k))
     $             /(ff(i,j,k)+ff(i,j,k-1)+epsilon)
              zwflux(i,j,k)=(wdz-w2dt)*mol*sw
              abs_1=abs(wdz)
              abs_2=abs(w2dt)
              if(abs_1.lt.abs_2)zwflux(i,j,k)=0.e0
            end if
          end do
        end do
      end do

      return
      end

!_______________________________________________________________________
      subroutine smol_adif_dble(xmassflux,ymassflux,zwflux,ff)
! calculate the antidiffusive velocity used to reduce the numerical
! diffusion associated with the upstream differencing scheme
! this is based on a subroutine of Gianmaria Sannino (Inter-university
! Computing Consortium, Rome, Italy) and Vincenzo Artale (Italian
! National Agency for New Technology and Environment, Rome, Italy)
!     S. Ohishi 2019.03 real --> double precision
      implicit none
      include 'pom.h'
      double precision ff(im,jm,kb)
      double precision xmassflux(im,jm,kb),ymassflux(im,jm,kb)
      double precision zwflux(im,jm,kb)
      double precision mol,abs_1,abs_2
      double precision value_min,epsilon
      double precision udx,u2dt,vdy,v2dt,wdz,w2dt
      integer i,j,k
      parameter (value_min=1.d-9,epsilon=1.0d-14)

! apply temperature and salinity mask
      do k=1,kb
         do i=1,im
            do j=1,jm
               ff(i,j,k)=ff(i,j,k)*dble(fsm(i,j))
            end do
         end do
      end do
      
!     recalculate mass fluxes with antidiffusion velocity
      do k=1,kbm1
         do j=2,jmm1
            do i=2,im
               if(ff(i,j,k).lt.value_min.or.
     $              ff(i-1,j,k).lt.value_min) then
                  xmassflux(i,j,k)=0.d0
               else
                  udx=dabs(xmassflux(i,j,k))
                  u2dt=
     $                 dble(dti2)*xmassflux(i,j,k)*xmassflux(i,j,k)*2.d0
     $                 /(dble(aru(i,j))*(dble(dt(i-1,j))+dble(dt(i,j))))
                  mol=(ff(i,j,k)-ff(i-1,j,k))
     $                 /(ff(i-1,j,k)+ff(i,j,k)+epsilon)
                  xmassflux(i,j,k)=dble((udx-u2dt)*mol*sw)
                  abs_1=dabs(udx)
                  abs_2=dabs(u2dt)
                  if(abs_1.lt.abs_2) xmassflux(i,j,k)=0.d0
               end if
            end do
         end do
      end do
      
      do k=1,kbm1
         do j=2,jm
            do i=2,imm1
               if(ff(i,j,k).lt.value_min.or.
     $              ff(i,j-1,k).lt.value_min) then
                  ymassflux(i,j,k)=0.d0
               else
                  vdy=dabs(ymassflux(i,j,k))
                  v2dt=
     $                 dble(dti2)*ymassflux(i,j,k)*ymassflux(i,j,k)*2.d0
     $                 /(dble(arv(i,j))*(dble(dt(i,j-1))+dble(dt(i,j))))
                  mol=(ff(i,j,k)-ff(i,j-1,k))
     $                 /(ff(i,j-1,k)+ff(i,j,k)+epsilon)
                  ymassflux(i,j,k)=dble((vdy-v2dt)*mol*sw)
                  abs_1=dabs(vdy)
                  abs_2=dabs(v2dt)
                  if(abs_1.lt.abs_2) ymassflux(i,j,k)=0.d0
               end if
            end do
         end do
      end do
      
      do k=2,kbm1
         do j=2,jmm1
            do i=2,imm1
               if(ff(i,j,k).lt.value_min.or.
     $              ff(i,j,k-1).lt.value_min) then
                  zwflux(i,j,k)=0.d0
               else
                  wdz=dabs(zwflux(i,j,k))
                  w2dt=
     $                 dble(dti2)*zwflux(i,j,k)*zwflux(i,j,k)
     $                 /(dble(dzz(k-1))*dble(dt(i,j)))
                  mol=(ff(i,j,k-1)-ff(i,j,k))
     $                 /(ff(i,j,k)+ff(i,j,k-1)+epsilon)
                  zwflux(i,j,k)=dble((wdz-w2dt)*mol*sw)
                  abs_1=dabs(wdz)
                  abs_2=dabs(w2dt)
                  if(abs_1.lt.abs_2)zwflux(i,j,k)=0.d0
               end if
            end do
         end do
      end do
      
      return
      end

!_______________________________________________________________________
! kii2b
!      subroutine vertvl(xflux,yflux)
      subroutine vertvl
! calculates vertical velocity
      implicit none
      include 'pom.h'
      real xflux(im,jm,kb),yflux(im,jm,kb)
      integer i,j,k

! reestablish boundary conditions
      do k=1,kbm1
        do j=2,jm
          do i=2,im
            xflux(i,j,k)=.25e0*(dy(i,j)+dy(i-1,j))
     $                    *(dt(i,j)+dt(i-1,j))*u(i,j,k)
          end do
        end do
      end do

      do k=1,kbm1
        do j=2,jm
          do i=2,im
            yflux(i,j,k)=.25e0*(dx(i,j)+dx(i,j-1))
     $                    *(dt(i,j)+dt(i,j-1))*v(i,j,k)
          end do
        end do
      end do

! note: if one wishes to include freshwater flux, the surface velocity
! should be set to vflux(i,j). See also change made to 2-D volume
! conservation equation which calculates elf
        do j=2,jmm1
          do i=2,imm1
            w(i,j,1)=0.5*(vfluxb(i,j)+vfluxf(i,j))
          end do
        end do

      do k=1,kbm1
        do j=2,jmm1
          do i=2,imm1
            w(i,j,k+1)=w(i,j,k)
     $                +dz(k)*((xflux(i+1,j,k)-xflux(i,j,k)
     $                        +yflux(i,j+1,k)-yflux(i,j,k))
     $                        /(dx(i,j)*dy(i,j))
     $                        +(etf(i,j)-etb(i,j))/dti2)
          end do
        end do
      end do

      return
      end

! kii2b
!_______________________________________________________________________
      subroutine realvertvl
! calculates real vertical velocity (wr)
      implicit none
      include 'pom.h'
      integer i,j,k
      real dxr,dxl,dyt,dyb
      do k=1,kb
        do j=1,jm
          do i=1,im
            wr(i,j,k)=0.
          end do
        end do
      end do

      do k=1,kbm1
        do j=1,jm
          do i=1,im
            tps(i,j)=zz(k)*dt(i,j) + et(i,j)
          end do
        end do
        do j=2,jmm1
          do i=2,imm1
            dxr=2.0/(dx(i+1,j)+dx(i,j))
            dxl=2.0/(dx(i,j)+dx(i-1,j))
            dyt=2.0/(dy(i,j+1)+dy(i,j))
            dyb=2.0/(dy(i,j)+dy(i,j-1))
            wr(i,j,k)=0.5*(w(i,j,k)+w(i,j,k+1))+0.5*
     $                (u(i+1,j,k)*(tps(i+1,j)-tps(i,j))*dxr+
     $                 u(i,j,k)*(tps(i,j)-tps(i-1,j))*dxl+
     $                 v(i,j+1,k)*(tps(i,j+1)-tps(i,j))*dyt+
     $                 v(i,j,k)*(tps(i,j)-tps(i,j-1))*dyb)
     $                +(1.0+zz(k))*(etf(i,j)-etb(i,j))/dti2
          end do
        end do
      end do

      call exchange3d_mpi(wr(:,:,1:kbm1),im,jm,kbm1)

      do k=1,kb
        do i=1,im
          if(n_south.eq.-1) wr(i,1,k)=wr(i,2,k)
          if(n_north.eq.-1) wr(i,jm,k)=wr(i,jmm1,k)
        end do
      end do
      do k=1,kb
        do j=1,jm
          if(n_west.eq.-1) wr(1,j,k)=wr(2,j,k)
          if(n_east.eq.-1) wr(im,j,k)=wr(imm1,j,k)
        end do
      end do

      do k=1,kbm1
        do j=1,jm
          do i=1,im
            wr(i,j,k)=fsm(i,j)*wr(i,j,k)
          end do
        end do
      end do

      return
      end

