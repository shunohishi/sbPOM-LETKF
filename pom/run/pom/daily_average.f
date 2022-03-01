!     daily_average.f 

!     average variables within 1 day
!_________________________________________________________________

      subroutine daily_average

      implicit none
      include 'pom.h'

      integer i,j
      real duvm(im,jm)

      do j=1,jm
         do i=1,im
            duvm(i,j)=dum(i,j)*dvm(i,j)
         end do
      end do

!SSH
      call calculate_daily_average(iint,iprint,im_local,jm_local,1,
     &     el,el_dave,fsm)
!Latent heat flux
      call calculate_daily_average(iint,iprint,im_local,jm_local,1,
     &     lhf,lhf_dave,fsm)
!Sensible heat flux
      call calculate_daily_average(iint,iprint,im_local,jm_local,1,
     &     shf,shf_dave,fsm)
!Longwave radiation
      call calculate_daily_average(iint,iprint,im_local,jm_local,1,
     &     lwr,lwr_dave,fsm)
!Shortwave radiation
      call calculate_daily_average(iint,iprint,im_local,jm_local,1,
     &     swr,swr_dave,fsm)
!zonal wind
      call calculate_daily_average(iint,iprint,im_local,jm_local,1,
     &     windu,windu_dave,fsm)
!meridional wind
      call calculate_daily_average(iint,iprint,im_local,jm_local,1,
     &     windv,windv_dave,fsm)
!wind speed
      call calculate_daily_average(iint,iprint,im_local,jm_local,1,
     &     winds,winds_dave,fsm)
!zonal wind stress
      call calculate_daily_average(iint,iprint,im_local,jm_local,1,
     &     tauu,tauu_dave,dum)
!meridional wind stress
      call calculate_daily_average(iint,iprint,im_local,jm_local,1,
     &     tauv,tauv_dave,dvm)
!magnitude of wind stress
      call calculate_daily_average(iint,iprint,im_local,jm_local,1,
     &     taus,taus_dave,duvm)
!air humidity
      call calculate_daily_average(iint,iprint,im_local,jm_local,1,
     &     qa,qa_dave,fsm)
!surface saturated humidity
      call calculate_daily_average(iint,iprint,im_local,jm_local,1,
     &     qs,qs_dave,fsm)
!air temperature
      call calculate_daily_average(iint,iprint,im_local,jm_local,1,
     &     ta,ta_dave,fsm)
! freshwater flux      
      if(lssf == 1)then
!     Evaporation
         call calculate_daily_average(iint,iprint,im_local,jm_local,1,
     &        evap,evap_dave,fsm)
!     Evaporation flux
         call calculate_daily_average(iint,iprint,im_local,jm_local,1,
     &        eflux,eflux_dave,fsm)
!     Precipitation
         call calculate_daily_average(iint,iprint,im_local,jm_local,1,
     &        prep,prep_dave,fsm)
!     Precipitation flux
         call calculate_daily_average(iint,iprint,im_local,jm_local,1,
     &        pflux,pflux_dave,fsm)
!     River discharge
         call calculate_daily_average(iint,iprint,im_local,jm_local,1,
     &        river,river_dave,fsm)
!     River discharge flux
         call calculate_daily_average(iint,iprint,im_local,jm_local,1,
     &        rflux,rflux_dave,fsm)         
      end if

!     S
      call calculate_daily_average(iint,iprint,im_local,jm_local,kb,
     &     s,s_dave,fsm)
!     T
      call calculate_daily_average(iint,iprint,im_local,jm_local,kb,
     &     t,t_dave,fsm)
!     U
      call calculate_daily_average(iint,iprint,im_local,jm_local,kb,
     &     u,u_dave,dum)
!     V
      call calculate_daily_average(iint,iprint,im_local,jm_local,kb,
     &     v,v_dave,dvm)
!     W
      call calculate_daily_average(iint,iprint,im_local,jm_local,kb,
     &     wr,wr_dave,fsm)

!     When convervation each term of budget equations
      if(budget == 1)then
!     T tendency
         call calculate_daily_average(iint,iprint,im_local,jm_local,kb,
     &        dtdt,dtdt_dave,fsm)
!     T Advection
         call calculate_daily_average(iint,iprint,im_local,jm_local,kb,
     &        txadv,txadv_dave,fsm)
         call calculate_daily_average(iint,iprint,im_local,jm_local,kb,
     &        tyadv,tyadv_dave,fsm)
         call calculate_daily_average(iint,iprint,im_local,jm_local,kb,
     &        tzadv,tzadv_dave,fsm)
         call calculate_daily_average(iint,iprint,im_local,jm_local,kb,
     &        tadv,tadv_dave,fsm)
!     T Diffusion
         call calculate_daily_average(iint,iprint,im_local,jm_local,kb,
     &        txdif,txdif_dave,fsm)
         call calculate_daily_average(iint,iprint,im_local,jm_local,kb,
     &        tydif,tydif_dave,fsm)
         call calculate_daily_average(iint,iprint,im_local,jm_local,kb,
     &        tzdif,tzdif_dave,fsm)
!     Suface heat flux
         call calculate_daily_average(iint,iprint,im_local,jm_local,1,
     &        tsfc,tsfc_dave,fsm)
!     Shortwave penetration
         call calculate_daily_average(iint,iprint,im_local,jm_local,kb,
     &        qz,qz_dave,fsm)

!     S Tendency
         call calculate_daily_average(iint,iprint,im_local,jm_local,kb,
     &        dsdt,dsdt_dave,fsm)
!     S Advection
         call calculate_daily_average(iint,iprint,im_local,jm_local,kb,
     &        sxadv,sxadv_dave,fsm)
         call calculate_daily_average(iint,iprint,im_local,jm_local,kb,
     &        syadv,syadv_dave,fsm)
         call calculate_daily_average(iint,iprint,im_local,jm_local,kb,
     &        szadv,szadv_dave,fsm)
         call calculate_daily_average(iint,iprint,im_local,jm_local,kb,
     &        sadv,sadv_dave,fsm)
!     S Diffusion
         call calculate_daily_average(iint,iprint,im_local,jm_local,kb,
     &        sxdif,sxdif_dave,fsm)
         call calculate_daily_average(iint,iprint,im_local,jm_local,kb,
     &        sydif,sydif_dave,fsm)
         call calculate_daily_average(iint,iprint,im_local,jm_local,kb,
     &        szdif,szdif_dave,fsm)
!     Suface freshwater flux
         call calculate_daily_average(iint,iprint,im_local,jm_local,1,
     &        ssfc,ssfc_dave,fsm)

         if(lroff)then
!     Temperature round off
            call calculate_daily_average(iint,iprint,
     &           im_local,jm_local,kb,
     &           troff,troff_dave,fsm)
!     Salinity round off
            call calculate_daily_average(iint,iprint,
     &           im_local,jm_local,kb,
     &           sroff,sroff_dave,fsm)
         end if

!     Temprature nuding
         call calculate_daily_average(iint,iprint,
     &        im_local,jm_local,kb,
     &        tnudge,tnudge_dave,fsm)

!     Salinity nuding
         call calculate_daily_average(iint,iprint,
     &        im_local,jm_local,kb,
     &        snudge,snudge_dave,fsm)

!     Kinematic Visosity (aam) & Vertical Diffusivity (kh)
         call calculate_daily_average(iint,iprint,im_local,jm_local,kb,
     &        aam,aam_dave,fsm)
         call calculate_daily_average(iint,iprint,im_local,jm_local,kb,
     &        kh,kh_dave,fsm)
      end if

      end subroutine

!__________________________________________________________________

      subroutine calculate_daily_average(iint,nint,im,jm,km,dat,ave,fsm)

      implicit none
      integer iint,nint         !time step
      integer i,j,k,im,jm,km    !the number of grid
      real dat(im,jm,km)        !original data
      real ave(im,jm,km)        !sum -> average
      real fsm(im,jm)           !mask

      if(mod(iint,nint) == 1)then
         ave(:,:,:)=0.
      end if

      do k=1,km
         do j=1,jm
            do i=1,im
               ave(i,j,k)=ave(i,j,k)+fsm(i,j)*dat(i,j,k)
            end do
         end do
      end do

      if(mod(iint,nint) == 0)then
         do k=1,km
            do j=1,jm
               do i=1,im
                  ave(i,j,k)=fsm(i,j)*ave(i,j,k)/real(nint)
               end do
            end do
         end do
      end if

      end subroutine
