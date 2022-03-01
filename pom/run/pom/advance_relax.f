! advance.f

! advance POM

!_______________________________________________________________________
      subroutine advance
! advance POM 1 step in time
      implicit none
      include 'pom.h'

! debug
      integer iglb,jglb,iloc,jloc,nproc
      integer i,j,k
      integer iuvmax,juvmax,kuvmax
      integer iuvmin,juvmin,kuvmin
      integer itmax,jtmax,ktmax
      integer itmin,jtmin,ktmin
      integer ismax,jsmax,ksmax
      integer ismin,jsmin,ksmin
      integer ielmax,jelmax,kelmax
      integer ielmin,jelmin,kelmin
      real uvmax,uvmin,tmax,tmin,smax,smin,elmax,elmin

      if(my_task .eq. master_task .and. mod(iint,100) .eq. 0)
     $     write(6,*) "start advance routine at ",iint 

! get time
      call get_time

! set time dependent surface boundary conditions
      call surface_forcing

! set time dependent lateral boundary conditions
      if(llbc .eq. 1)then
         call lateral_forcing_mclim
      elseif(llbc .eq. 2)then
         call lateral_forcing
      else
         write(*,"(/a)") "POM terminated with error: llbc"
        call finalize_mpi
        stop
      endif

! set time dependent reference temperature and salinity
      if(llbc .eq. 1)then
         call set_tsdata_mclim
      elseif(llbc .eq. 2)then
         call set_tsdata
      else
         write(*,"(/a)") "POM terminated with error: llbc"
         call finalize_mpi
         stop
      endif

! set monthly climatological data
      call set_tsclim

! set lateral viscosity
      call lateral_viscosity

! form vertical averages of 3-D fields for use in external (2-D) mode
      call mode_interaction

! external (2-D) mode calculation
      do iext=1,isplit
        call mode_external
      end do

!     Check NaN
      if(lnan)then
         call check_nan(1,"uaf",uaf)
         call check_nan(1,"vaf",vaf)
         call check_nan(1,"elf",elf)
         call check_nan(1,"etf",etf)
      end if

! internal (3-D) mode calculation
      call mode_internal

!     Restoring variable at sea bottom during assimilation
!     S.Ohishi 2018.12
      if(assim == 1 .or. assim == 2)then
!         call restore_bottom_ts
      end if

! print section
      call print_section

! uv point sampling
!      call uvstore

! uv statistics
!      call uvstatistics(0)

! daily average (S.Ohishi 2018.07)
      if(ldave == 1)then
         call daily_average
      end if

      if(netcdf_file.ne.'nonetcdf' .and. mod(iint,iprint).eq.0)then
         call write_output_netcdf
      end if
      
! uv statistics clear 
!      if(mod(iint,iprint).eq.0) call uvstatistics(1)

! write restart
      if(mod(iint,irestart).eq.0)then
!     if(my_task .eq. master_task) write(6,*) iint,"write restart"
         call write_restart_netcdf
      end if

      if(mod(iint,iens) .eq. 0)then
         call write_ens_netcdf
      end if

! check CFL condition
!      if(my_task .eq. master_task)
!     $     write(6,*) iint,"check CFL condition"
      call check_velocity

      if(my_task .eq. master_task .and. mod(iint,100) .eq. 0)then 
         write(6,*) "finish advance routine at ",iint 
      end if

      return
      end

!_______________________________________________________________________
      subroutine get_time
! return the model time
      implicit none
      include 'pom.h'

      time=dti*float(iint)/86400.e0+time0
      if(iint.ge.iswtch) iprint=nint(prtd2*24.e0*3600.e0/dti)
!      iprint=iint
      if(lramp) then
        ramp=time/period
        if(ramp.gt.1.e0) ramp=1.e0
      else
        ramp=1.e0
      end if
      return
      end

!_______________________________________________________________________
      subroutine surface_forcing
! set time dependent surface boundary conditions
      implicit none
      include 'pom.h'
      integer i,j
      real tatm,satm

      call surface_airseaflux
      call surface_freshwaterflux

      do j=2,jmm1
        do i=2,imm1

! wind stress
! value is negative for westerly or southerly winds. The wind stress
! should be tapered along the boundary to suppress numerically induced
! oscilations near the boundary (Jamart and Ozer, JGR, 91, 10621-10631)
!          wusurf(i,j)=0.e0
!          wvsurf(i,j)=0.e0

          e_atmos(i,j)=0.e0
          vfluxf(i,j)=0.e0

! set w(i,j,1)=vflux(i,j).ne.0 if one wishes non-zero flow across
! the sea surface. See calculation of elf(i,j) below and subroutines
! vertvl, advt1 (or advt2). If w(1,j,1)=0, and, additionally, there
! is no net flow across lateral boundaries, the basin volume will be
! constant; if also vflux(i,j).ne.0, then, for example, the average
! salinity will change and, unrealistically, so will total salt
          w(i,j,1)=vfluxf(i,j)

! set wtsurf to the sensible heat, the latent heat (which involves
! only the evaporative component of vflux) and the long wave
! radiation
!          wtsurf(i,j)=0.e0
!
! set swrad to the short wave radiation
!          swrad(i,j)=0.e0

! to account for change in temperature of flow crossing the sea
! surface (generally quite small compared to latent heat effect)
!          tatm=t(i,j,1)+tbias    ! an approximation
!          wtsurf(i,j)=wtsurf(i,j)+vfluxf(i,j)*(tatm-t(i,j,1)-tbias)
!          wtsurf(i,j)=50./4.1876e06*(tclim(i,j,1)-t(i,j,1))
!
! set the salinity of water vapor/precipitation which enters/leaves
! the atmosphere (or e.g., an ice cover)
!          satm=0.e0
!          wssurf(i,j)=            vfluxf(i,j)*(satm-s(i,j,1)-sbias)

        end do
      end do

      return
      end

!_______________________________________________________________________
      subroutine lateral_forcing_mclim
!     set monthly climatological lateral boudary forcing
!     Created by 2018.08 S.Ohishi

      implicit none
      include 'pom.h'

      integer i,j,k
      integer julday_time,year,month,day,isecflg
      real timepre,sec,ratio

      timepre=dti*float(iint-1)/86400.e0+time0
      julday_time=int(julday_start+timepre+0.001)
      call caldat(julday_time,year,month,day)

      if(iint .eq. 1)then
         if(day .lt. 15)then
            month=month-1
            if(month .lt. 1) month=12
         end if
         call read_lbc_mclim_netcdf(month)
      else
         sec=(timepre-int(timepre))*86400.
         isecflg=int(abs(sec))
         if(day .eq. 15 .and. isecflg .eq. 0)then
            call read_lbc_mclim_netcdf(month)
         end if
      end if

      call ratioclim(ratio,julday_time,timepre)
      if(mod(iint-1,100) .eq. 0 .and. my_task .eq. master_task)then
         write(6,*) 'lbcclim ratio:',ratio,"at ",iint
      end if

      do j=2,jmm1
        ele(j)=(1.-ratio)*ele0(j)+ratio*ele1(j)
        uabe(j)=(1.-ratio)*uabe0(j)+ratio*uabe1(j)
        vabe(j)=(1.-ratio)*vabe0(j)+ratio*vabe1(j)
        elw(j)=(1.-ratio)*elw0(j)+ratio*elw1(j)
        uabw(j)=(1.-ratio)*uabw0(j)+ratio*uabw1(j)
        vabw(j)=(1.-ratio)*vabw0(j)+ratio*vabw1(j)
      end do       

      do k=1,kbm1
        do j=2,jmm1
          tbe(j,k)=(1.-ratio)*tbe0(j,k)+ratio*tbe1(j,k)
          sbe(j,k)=(1.-ratio)*sbe0(j,k)+ratio*sbe1(j,k)
          ube(j,k)=(1.-ratio)*ube0(j,k)+ratio*ube1(j,k)
          vbe(j,k)=(1.-ratio)*vbe0(j,k)+ratio*vbe1(j,k)
          tbw(j,k)=(1.-ratio)*tbw0(j,k)+ratio*tbw1(j,k)
          sbw(j,k)=(1.-ratio)*sbw0(j,k)+ratio*sbw1(j,k)
          ubw(j,k)=(1.-ratio)*ubw0(j,k)+ratio*ubw1(j,k)
          vbw(j,k)=(1.-ratio)*vbw0(j,k)+ratio*vbw1(j,k)
        end do       
      end do

      do i=2,imm1
        eln(i)=(1.-ratio)*eln0(i)+ratio*eln1(i)
        uabn(i)=(1.-ratio)*uabn0(i)+ratio*uabn1(i)
        vabn(i)=(1.-ratio)*vabn0(i)+ratio*vabn1(i)
        els(i)=(1.-ratio)*els0(i)+ratio*els1(i)
        uabs(i)=(1.-ratio)*uabs0(i)+ratio*uabs1(i)
        vabs(i)=(1.-ratio)*vabs0(i)+ratio*vabs1(i)
      end do       

      do k=1,kbm1
        do i=2,imm1
          tbn(i,k)=(1.-ratio)*tbn0(i,k)+ratio*tbn1(i,k)
          sbn(i,k)=(1.-ratio)*sbn0(i,k)+ratio*sbn1(i,k)
          ubn(i,k)=(1.-ratio)*ubn0(i,k)+ratio*ubn1(i,k)
          vbn(i,k)=(1.-ratio)*vbn0(i,k)+ratio*vbn1(i,k)
          tbs(i,k)=(1.-ratio)*tbs0(i,k)+ratio*tbs1(i,k)
          sbs(i,k)=(1.-ratio)*sbs0(i,k)+ratio*sbs1(i,k)
          ubs(i,k)=(1.-ratio)*ubs0(i,k)+ratio*ubs1(i,k)
          vbs(i,k)=(1.-ratio)*vbs0(i,k)+ratio*vbs1(i,k)
        end do
      end do
      
      return
      end
      
!_______________________________________________________________________
      subroutine lateral_forcing
! set time dependent surface boundary conditions
      implicit none
      include 'pom.h'
      real time0_data,ratio
      integer timing,i,j,k

      if(iint.eq.1) call read_lbc_netcdf
   
      time0_data=lbctime_julday-julday_start
      call readtiming(timing,ratio,lbctime_dayint,time0_data)
      if(timing.eq.1) then
        call read_lbc_netcdf
        if(my_task.eq.master_task) then
          write(6,*) 'lbc timing ratio time ',timing,ratio
     $               ,dti*float(iint-1)/86400.e0+time0
        end if
      end if

      do j=2,jmm1
        ele(j)=(1.-ratio)*ele0(j)+ratio*ele1(j)
        uabe(j)=(1.-ratio)*uabe0(j)+ratio*uabe1(j)
        vabe(j)=(1.-ratio)*vabe0(j)+ratio*vabe1(j)
        elw(j)=(1.-ratio)*elw0(j)+ratio*elw1(j)
        uabw(j)=(1.-ratio)*uabw0(j)+ratio*uabw1(j)
        vabw(j)=(1.-ratio)*vabw0(j)+ratio*vabw1(j)
      end do       

      do k=1,kbm1
        do j=2,jmm1
          tbe(j,k)=(1.-ratio)*tbe0(j,k)+ratio*tbe1(j,k)
          sbe(j,k)=(1.-ratio)*sbe0(j,k)+ratio*sbe1(j,k)
          ube(j,k)=(1.-ratio)*ube0(j,k)+ratio*ube1(j,k)
          vbe(j,k)=(1.-ratio)*vbe0(j,k)+ratio*vbe1(j,k)
          tbw(j,k)=(1.-ratio)*tbw0(j,k)+ratio*tbw1(j,k)
          sbw(j,k)=(1.-ratio)*sbw0(j,k)+ratio*sbw1(j,k)
          ubw(j,k)=(1.-ratio)*ubw0(j,k)+ratio*ubw1(j,k)
          vbw(j,k)=(1.-ratio)*vbw0(j,k)+ratio*vbw1(j,k)
        end do       
      end do

      do i=2,imm1
        eln(i)=(1.-ratio)*eln0(i)+ratio*eln1(i)
        uabn(i)=(1.-ratio)*uabn0(i)+ratio*uabn1(i)
        vabn(i)=(1.-ratio)*vabn0(i)+ratio*vabn1(i)
        els(i)=(1.-ratio)*els0(i)+ratio*els1(i)
        uabs(i)=(1.-ratio)*uabs0(i)+ratio*uabs1(i)
        vabs(i)=(1.-ratio)*vabs0(i)+ratio*vabs1(i)
      end do       

      do k=1,kbm1
        do i=2,imm1
          tbn(i,k)=(1.-ratio)*tbn0(i,k)+ratio*tbn1(i,k)
          sbn(i,k)=(1.-ratio)*sbn0(i,k)+ratio*sbn1(i,k)
          ubn(i,k)=(1.-ratio)*ubn0(i,k)+ratio*ubn1(i,k)
          vbn(i,k)=(1.-ratio)*vbn0(i,k)+ratio*vbn1(i,k)
          tbs(i,k)=(1.-ratio)*tbs0(i,k)+ratio*tbs1(i,k)
          sbs(i,k)=(1.-ratio)*sbs0(i,k)+ratio*sbs1(i,k)
          ubs(i,k)=(1.-ratio)*ubs0(i,k)+ratio*ubs1(i,k)
          vbs(i,k)=(1.-ratio)*vbs0(i,k)+ratio*vbs1(i,k)
        end do
      end do

      return
      end

!_______________________________________________________________________
      
      subroutine set_tsdata_mclim
!     monthly climatology temperature and salinity
!     2018.08.20 S.Ohishi

      implicit none
      include 'pom.h'

      integer i,j,k
      integer julday_time,year,month,day,isecflg
      real timepre,sec,ratio

      timepre=dti*float(iint-1)/86400.e0+time0
      julday_time=int(julday_start+timepre+0.001)
      call caldat(julday_time,year,month,day)

      if(iint .eq. 1)then
         if(day .lt. 15)then
            month=month-1
            if(month .lt. 1) month=12
         end if
         call read_tsdata_mclim_netcdf(month)
      else
         sec=(timepre-int(timepre))*86400.
         isecflg=int(abs(sec))
         if(day .eq. 15 .and. isecflg .eq. 0)then
            call read_tsdata_mclim_netcdf(month)
         end if
      end if

      call ratioclim(ratio,julday_time,timepre)
      if(mod(iint-1,100) .eq. 0 .and. my_task .eq. master_task)then
         write(6,*) 'tsdataclim ratio:',ratio,"at ",iint
      end if
      
      do k=1,kb
         do j=1,jm
            do i=1,im
            tref(i,j,k)=(1.-ratio)*tref0(i,j,k)
     $                   +ratio*tref1(i,j,k)
            sref(i,j,k)=(1.-ratio)*sref0(i,j,k)
     $                   +ratio*sref1(i,j,k)
            end do
         end do
      end do

      return
      end

!_______________________________________________________________________
      subroutine set_tsdata
! set time dependent reference temperature and salinity
      implicit none
      include 'pom.h'
      real time0_data,ratio
      integer timing,i,j,k

      if(iint.eq.1) call read_tsdata_netcdf
   
      time0_data=tsdatatime_julday-julday_start
      call readtiming(timing,ratio,tsdatatime_dayint,time0_data)
      if(timing.eq.1) then
        call read_tsdata_netcdf
        if(my_task.eq.master_task) then
          write(6,*) 'tsdata timing ratio time ',timing,ratio
     $               ,dti*float(iint-1)/86400.e0+time0
        end if
      end if

      do k=1,kb
        do j=1,jm
          do i=1,im
            tref(i,j,k)=(1.-ratio)*tref0(i,j,k)
     $                   +ratio*tref1(i,j,k)
            sref(i,j,k)=(1.-ratio)*sref0(i,j,k)
     $                   +ratio*sref1(i,j,k)
          end do
        end do       
      end do

!      if(my_task.eq.master_task) then
!        i = im/2
!        j = jm/2
!        k = 1
!        write(6,*) 'tref0 tref1 tref '
!     $   ,tref0(i,j,k),tref1(i,j,k),tref(i,j,k)
!        write(6,*) 'sref0 sref1 sref '
!     $   ,sref0(i,j,k),sref1(i,j,k),sref(i,j,k)
!      end if

      return
      end

!_______________________________________________________________________
      subroutine set_tsclim
! set monthly climatological temperature and salinity data
      implicit none
      include 'pom.h'
      real timepre,sec,ratio
      integer timing,i,j,k
      integer julday_time,year,month,day,isecflg

      timepre=dti*float(iint-1)/86400.e0+time0
      julday_time=int(julday_start+timepre+0.001)
      call caldat(julday_time,year,month,day)
     
      if(iint.eq.1) then
        if(day.lt.15) then
          month=month-1
          if(month.lt.1) month=12
        end if
        call read_tsclim_monthly_netcdf(month)
      else
        sec=(timepre-int(timepre))*86400.
        isecflg=int(abs(sec))
        if(day.eq.15.and.isecflg.eq.0) then
          call read_tsclim_monthly_netcdf(month)
        end if
      end if   

      call ratioclim(ratio,julday_time,timepre)
      if(mod(iint-1,100).eq.0.and.my_task.eq.master_task) then
        write(6,*) 'tsclim ratio:', ratio ,"at ",iint
      end if

      do k=1,kb
        do j=1,jm
          do i=1,im
            tclimm(i,j,k)=(1.-ratio)*tclim0(i,j,k)
     $                   +ratio*tclim1(i,j,k)
            sclimm(i,j,k)=(1.-ratio)*sclim0(i,j,k)
     $                   +ratio*sclim1(i,j,k)
          end do
        end do       
      end do

      return
      end

!_______________________________________________________________________
      subroutine lateral_viscosity
! set the lateral viscosity
      implicit none
      include 'pom.h'
      integer i,j,k
! if mode=2 then initial values of aam2d are used. If one wishes
! to use Smagorinsky lateral viscosity and diffusion for an
! external (2-D) mode calculation, then appropiate code can be
! adapted from that below and installed just before the end of the
! "if(mode.eq.2)" loop in subroutine advave

! calculate Smagorinsky lateral viscosity:
! ( hor visc = horcon*dx*dy*sqrt((du/dx)**2+(dv/dy)**2
!                                +.5*(du/dy+dv/dx)**2) )
      if(mode.ne.2) then

! kii2b
!        call advct(a,c,ee)
        call advct

        if (npg .eq. 1) then
          call baropg
        else if (npg .eq. 2) then
          call baropg_mcc
        else if (npg .eq. 3) then
          call baropg_thiem
        else
          error_status=1
          write(6,'(/''Error: invalid value for npg'')')
        end if

        do k=1,kbm1
          do j=2,jmm1
            do i=2,imm1
              aam(i,j,k)=horcon*dx(i,j)*dy(i,j)
     $                    *sqrt( ((u(i+1,j,k)-u(i,j,k))/dx(i,j))**2
     $                          +((v(i,j+1,k)-v(i,j,k))/dy(i,j))**2
     $                    +.5e0*(.25e0*(u(i,j+1,k)+u(i+1,j+1,k)
     $                                 -u(i,j-1,k)-u(i+1,j-1,k))
     $                    /dy(i,j)
     $                    +.25e0*(v(i+1,j,k)+v(i+1,j+1,k)
     $                           -v(i-1,j,k)-v(i-1,j+1,k))
     $                    /dx(i,j)) **2)
     $                 +aamadd
            end do
          end do
        end do

! create sponge zones
!        do k=1,kbm1
!          do j=2,jmm1
!            do i=2,imm1
!              aam(i,j,k)=aam(i,j,k)+1000*exp(-(j_global(j)-2)*1.5)
!     $                    +1000*exp((j_global(j)-jm_global+1)*1.5)
!            end do
!          end do
!        end do
        call exchange3d_mpi(aam(:,:,1:kbm1),im,jm,kbm1)
      end if

      return
      end

!_______________________________________________________________________
      subroutine mode_interaction
! form vertical averages of 3-D fields for use in external (2-D) mode
      implicit none
      include 'pom.h'
      integer i,j,k

      if(mode.ne.2) then

        do j=1,jm
          do i=1,im
            adx2d(i,j)=0.e0
            ady2d(i,j)=0.e0
            drx2d(i,j)=0.e0
            dry2d(i,j)=0.e0
            aam2d(i,j)=0.e0
          end do
        end do

        do k=1,kbm1
          do j=1,jm
            do i=1,im
              adx2d(i,j)=adx2d(i,j)+advx(i,j,k)*dz(k)
              ady2d(i,j)=ady2d(i,j)+advy(i,j,k)*dz(k)
              drx2d(i,j)=drx2d(i,j)+drhox(i,j,k)*dz(k)
              dry2d(i,j)=dry2d(i,j)+drhoy(i,j,k)*dz(k)
              aam2d(i,j)=aam2d(i,j)+aam(i,j,k)*dz(k)
            end do
          end do
        end do

        call exchange2d_mpi(drx2d,im,jm) !2018.08
        call exchange2d_mpi(dry2d,im,jm) !2018.08
        call exchange2d_mpi(aam2d,im,jm) !2018.08

! kii2b
!        call advave(tps)
        call advave

        do j=1,jm
          do i=1,im
            adx2d(i,j)=adx2d(i,j)-advua(i,j)
            ady2d(i,j)=ady2d(i,j)-advva(i,j)
          end do
        end do

        call exchange2d_mpi(adx2d,im,jm) !2018.08
        call exchange2d_mpi(ady2d,im,jm) !2018.08

      end if

      do j=1,jm
        do i=1,im
          egf(i,j)=el(i,j)*ispi
        end do
      end do

      do j=1,jm
        do i=2,im
          utf(i,j)=ua(i,j)*(d(i,j)+d(i-1,j))*isp2i
        end do
      end do
      do j=2,jm
        do i=1,im
          vtf(i,j)=va(i,j)*(d(i,j)+d(i,j-1))*isp2i
        end do
      end do

      call exchange2d_mpi(utf,im,jm) !2018.08
      call exchange2d_mpi(vtf,im,jm) !2018.08

      return
      end

!_______________________________________________________________________
      subroutine mode_external
! calculate the external (2-D) mode
      implicit none
      include 'pom.h'
      integer i,j

      do j=2,jm
        do i=2,im
          fluxua(i,j)=.25e0*(d(i,j)+d(i-1,j))
     $                 *(dy(i,j)+dy(i-1,j))*ua(i,j)
          fluxva(i,j)=.25e0*(d(i,j)+d(i,j-1))
     $                 *(dx(i,j)+dx(i,j-1))*va(i,j)
        end do
      end do

! NOTE addition of surface freshwater flux, w(i,j,1)=vflux, compared
! with pom98.f. See also modifications to subroutine vertvl
      do j=2,jmm1
        do i=2,imm1
          elf(i,j)=elb(i,j)
     $              +dte2*(-(fluxua(i+1,j)-fluxua(i,j)
     $                      +fluxva(i,j+1)-fluxva(i,j))/art(i,j)
     $                      -vfluxf(i,j))
        end do
      end do

      if(ltide) then
        call tide_bcond(1)
      else
        call bcond(1)
      end if

!     S.Ohishi (2018.12)
!      if(assim == 1)then
!         call cal_el_iau
!      end if

      call exchange2d_mpi(elf,im,jm)

! kii2b
!      if(mod(iext,ispadv).eq.0) call advave(tps)
      if(mod(iext,ispadv).eq.0) call advave

      do j=2,jmm1
        do i=2,im
          uaf(i,j)=adx2d(i,j)+advua(i,j)
     $              -aru(i,j)*.25e0
     $                *(cor(i,j)*d(i,j)*(va(i,j+1)+va(i,j))
     $                 +cor(i-1,j)*d(i-1,j)*(va(i-1,j+1)+va(i-1,j)))
     $              +.25e0*grav*(dy(i,j)+dy(i-1,j))
     $                *(d(i,j)+d(i-1,j))
     $                *((1.e0-2.e0*alpha)
     $                   *(el(i,j)-el(i-1,j))
     $                  +alpha*(elb(i,j)-elb(i-1,j)
     $                         +elf(i,j)-elf(i-1,j))
     $                  +e_atmos(i,j)-e_atmos(i-1,j))
     $              +drx2d(i,j)+aru(i,j)*(wusurf(i,j)-wubot(i,j))
        end do
      end do

      do j=2,jmm1
        do i=2,im
          uaf(i,j)=((h(i,j)+elb(i,j)+h(i-1,j)+elb(i-1,j))
     $                *aru(i,j)*uab(i,j)
     $              -4.e0*dte*uaf(i,j))
     $             /((h(i,j)+elf(i,j)+h(i-1,j)+elf(i-1,j))
     $                 *aru(i,j))
        end do
      end do

      do j=2,jm
        do i=2,imm1
          vaf(i,j)=ady2d(i,j)+advva(i,j)
     $              +arv(i,j)*.25e0
     $                *(cor(i,j)*d(i,j)*(ua(i+1,j)+ua(i,j))
     $                 +cor(i,j-1)*d(i,j-1)*(ua(i+1,j-1)+ua(i,j-1)))
     $              +.25e0*grav*(dx(i,j)+dx(i,j-1))
     $                *(d(i,j)+d(i,j-1))
     $                *((1.e0-2.e0*alpha)*(el(i,j)-el(i,j-1))
     $                  +alpha*(elb(i,j)-elb(i,j-1)
     $                         +elf(i,j)-elf(i,j-1))
     $                  +e_atmos(i,j)-e_atmos(i,j-1))
     $              +dry2d(i,j)+arv(i,j)*(wvsurf(i,j)-wvbot(i,j))
        end do
      end do

      do j=2,jm
        do i=2,imm1
          vaf(i,j)=((h(i,j)+elb(i,j)+h(i,j-1)+elb(i,j-1))
     $                *vab(i,j)*arv(i,j)
     $              -4.e0*dte*vaf(i,j))
     $             /((h(i,j)+elf(i,j)+h(i,j-1)+elf(i,j-1))
     $                 *arv(i,j))
        end do
      end do

      if(ltide) then
        call tide_bcond(2)
      else 
        call bcond(2)
      end if

      call exchange2d_mpi(uaf,im,jm)
      call exchange2d_mpi(vaf,im,jm)

      if(iext.eq.(isplit-2))then
        do j=1,jm
          do i=1,im
            etf(i,j)=.25e0*smoth*elf(i,j)
          end do
        end do

      else if(iext.eq.(isplit-1)) then

        do j=1,jm
          do i=1,im
            etf(i,j)=etf(i,j)+.5e0*(1.-.5e0*smoth)*elf(i,j)
          end do
        end do

      else if(iext.eq.isplit) then

        do j=1,jm
          do i=1,im
            etf(i,j)=(etf(i,j)+.5e0*elf(i,j))*fsm(i,j)
          end do
        end do

      end if

! apply filter to remove time split
      do j=1,jm
        do i=1,im
          ua(i,j)=ua(i,j)+.5e0*smoth*(uab(i,j)-2.e0*ua(i,j)+uaf(i,j))
          va(i,j)=va(i,j)+.5e0*smoth*(vab(i,j)-2.e0*va(i,j)+vaf(i,j))
          el(i,j)=el(i,j)+.5e0*smoth*(elb(i,j)-2.e0*el(i,j)+elf(i,j))
          elb(i,j)=el(i,j)
          el(i,j)=elf(i,j)
          d(i,j)=h(i,j)+el(i,j)
          uab(i,j)=ua(i,j)
          ua(i,j)=uaf(i,j)
          vab(i,j)=va(i,j)
          va(i,j)=vaf(i,j)
        end do
      end do

      if(iext.ne.isplit) then
        do j=1,jm
          do i=1,im
            egf(i,j)=egf(i,j)+el(i,j)*ispi
          end do
        end do
        do j=1,jm
          do i=2,im
            utf(i,j)=utf(i,j)+ua(i,j)*(d(i,j)+d(i-1,j))*isp2i
          end do
        end do
        do j=2,jm
          do i=1,im
            vtf(i,j)=vtf(i,j)+va(i,j)*(d(i,j)+d(i,j-1))*isp2i
          end do
        end do
       end if

      call exchange2d_mpi(utf,im,jm) !2018.08
      call exchange2d_mpi(vtf,im,jm) !2018.08

      return
      end

!_______________________________________________________________________
      subroutine mode_internal
!     calculate the internal (3-D) mode
      implicit none
      include 'pom.h'
      integer i,j,k
      real dxr,dxl,dyt,dyb

      if((iint.ne.1.or.time0.ne.0.e0).and.mode.ne.2) then

!     adjust u(z) and v(z) such that depth average of (u,v) = (ua,va)
         do j=1,jm
            do i=1,im
               tps(i,j)=0.e0
            end do
         end do

         do k=1,kbm1
            do j=1,jm
               do i=1,im
                  tps(i,j)=tps(i,j)+u(i,j,k)*dz(k)
               end do
            end do
         end do

         do k=1,kbm1
            do j=1,jm
               do i=2,im
                  u(i,j,k)=(u(i,j,k)-tps(i,j))+
     $                 (utb(i,j)+utf(i,j))/(dt(i,j)+dt(i-1,j))
               end do
            end do
         end do
         
         do j=1,jm
            do i=1,im
               tps(i,j)=0.e0
            end do
         end do

         do k=1,kbm1
            do j=1,jm
               do i=1,im
                  tps(i,j)=tps(i,j)+v(i,j,k)*dz(k)
               end do
            end do
         end do

         do k=1,kbm1
            do j=2,jm
               do i=1,im
                  v(i,j,k)=(v(i,j,k)-tps(i,j))+
     $                 (vtb(i,j)+vtf(i,j))/(dt(i,j)+dt(i,j-1))
               end do
            end do
         end do

         call exchange3d_mpi(u(:,:,1:kbm1),im,jm,kbm1) ! 2018.07.27
         call exchange3d_mpi(v(:,:,1:kbm1),im,jm,kbm1) ! 2018.07.27

!     calculate w from u, v, dt (h+et), etf and etb
         call vertvl

         if(ltide) then
            call tide_bcond(5)
         else
            call bcond(5)
         end if

         call exchange3d_mpi(w,im,jm,kb)


!     set uf and vf to zero
         do k=1,kb
            do j=1,jm
               do i=1,im
                  uf(i,j,k)=0.e0
                  vf(i,j,k)=0.e0
               end do
            end do
         end do

!     calculate q2f and q2lf using uf, vf, a and c as temporary variables
         call advq(q2b,q2,uf)
         if(lmynnf) then
            call profq_mynnf
         else
            call advq(q2lb,q2l,vf)
            call profq
         end if

         if(ltide) then
            call tide_bcond(6)
         else
            call bcond(6)
         end if

         call exchange3d_mpi(uf(:,:,2:kbm1),im,jm,kbm2)
         if(.not.lmynnf) call exchange3d_mpi(vf(:,:,2:kbm1),im,jm,kbm2)

         if(lmynnf) then

            do k=1,kb
               do j=1,jm
                  do i=1,im
                     q2(i,j,k)=q2(i,j,k)
     $                    +.5e0*smoth*(uf(i,j,k)+q2b(i,j,k)
     $                    -2.e0*q2(i,j,k))
                     q2b(i,j,k)=q2(i,j,k)
                     q2(i,j,k)=uf(i,j,k)
                  end do
               end do
            end do

         else

            do k=1,kb
               do j=1,jm
                  do i=1,im
                     q2(i,j,k)=q2(i,j,k)
     $                    +.5e0*smoth*(uf(i,j,k)+q2b(i,j,k)
     $                    -2.e0*q2(i,j,k))
                     q2l(i,j,k)=q2l(i,j,k)
     $                    +.5e0*smoth*(vf(i,j,k)+q2lb(i,j,k)
     $                    -2.e0*q2l(i,j,k))
                     q2b(i,j,k)=q2(i,j,k)
                     q2(i,j,k)=uf(i,j,k)
                     q2lb(i,j,k)=q2l(i,j,k)
                     q2l(i,j,k)=vf(i,j,k)
                  end do
               end do
            end do

         end if

!     calculate tf and sf using uf, vf, a and c as temporary variables
         if(mode.ne.4)then

!     Temperature advection
            if(nadv.eq.1 .and. budget .eq. 1)then
               call advt1(tb,t,tclim,uf)
            else if(nadv .eq. 1)then
               call advt1_original(tb,t,tclim,uf)
            else if(nadv .eq. 2 .and. budget .eq. 1)then
               call advt2(tb,t,tclim,uf)
            else if(nadv .eq. 2)then
               call advt2_original(tb,t,tclim,uf)
            else
               error_status=1
               write(6,'(/''Error: invalid value for nadv'')')
            end if

            if(budget .eq. 1)then
               do k=1,kb
                  do j=1,jm
                     do i=1,im
                        txadv(i,j,k)=xadvterm(i,j,k)*86400.e0 !*86400.:[deg C/sec] -> [deg C/day]
                        tyadv(i,j,k)=yadvterm(i,j,k)*86400.e0
                        tzadv(i,j,k)=zadvterm(i,j,k)*86400.e0
                        tadv(i,j,k)=advterm(i,j,k)*86400.e0
                        txdif(i,j,k)=xdifterm(i,j,k)*86400.e0
                        tydif(i,j,k)=ydifterm(i,j,k)*86400.e0
                     end do
                  end do
               end do
            end if

!     Salinity advection
            if(nadv .eq. 1 .and. budget .eq. 1)then
               call advt1(sb,s,sclim,vf)
            else if(nadv .eq. 1)then
               call advt1_original(sb,s,sclim,vf)
            else if(nadv .eq. 2 .and. budget .eq. 1)then
               call advt2(sb,s,sclim,vf)
            else if(nadv .eq. 2)then
               call advt2_original(sb,s,sclim,vf)
            else
               error_status=1
               write(6,'(/''Error: invalid value for nadv'')')
            end if

            if(budget .eq. 1)then
               do k=1,kb
                  do j=1,jm
                     do i=1,im
                        sxadv(i,j,k)=xadvterm(i,j,k)*86400.e0
                        syadv(i,j,k)=yadvterm(i,j,k)*86400.e0
                        szadv(i,j,k)=zadvterm(i,j,k)*86400.e0
                        sadv(i,j,k)=advterm(i,j,k)*86400.e0
                        sxdif(i,j,k)=xdifterm(i,j,k)*86400.e0
                        sydif(i,j,k)=ydifterm(i,j,k)*86400.e0
                     end do
                  end do
               end do
            end if
            
!     Temperature diffusion
            call proft(uf,wtsurf,tsurf,nbct)
            if(budget == 1)then
               do k=1,kb
                  do j=1,jm
                     do i=1,im
                        tzdif(i,j,k)=zdifterm(i,j,k)*86400.e0
                        qz(i,j,k)=radterm(i,j,k)*86400.e0
                     end do
                  end do
               end do
               do j=1,jm
                  do i=1,im
                     tsfc(i,j)=sfcterm(i,j)*86400.e0
                  end do
               end do
            end if

!     Salinity diffusion
            call proft(vf,wssurf,ssurf,nbcs)
            if(budget == 1)then
               do k=1,kb
                  do j=1,jm
                     do i=1,im
                        szdif(i,j,k)=zdifterm(i,j,k)*86400.e0
                     end do
                  end do
               end do
               do j=1,jm
                  do i=1,im
                     ssfc(i,j)=sfcterm(i,j)*86400.e0
                  end do
               end do              
            end if
            
            if(ltide) then
               call tide_bcond(4)
            else
               call bcond(4)
            end if
            
            if(lfrs) call bcond_frs(4)
            
!     S.Ohishi(2018.12)
            if(assim == 1)then
               call ts_iau
            end if

!     S.Ohishi(2020.04)
            call ts_nudging

!     S.Ohishi (2020.02)
            if(lroff)then
               call ts_roff
            end if


            call exchange3d_mpi(uf(:,:,1:kbm1),im,jm,kbm1)
            call exchange3d_mpi(vf(:,:,1:kbm1),im,jm,kbm1)

!     S.Ohishi(2019.04)
            if(budget == 1)then

!     k=1,kbm1
               do k=1,kbm1
                  do j=1,jm
                     do i=1,im
                        dtdt(i,j,k)=(uf(i,j,k)-tb(i,j,k))/dti2*86400.e0
                        dsdt(i,j,k)=(vf(i,j,k)-sb(i,j,k))/dti2*86400.e0
                     end do
                  end do
               end do
               
!     k=kb
               do j=1,jm
                  do i=1,im
                     dtdt(i,j,kb)=0.e0
                     dsdt(i,j,kb)=0.e0
                  end do
               end do
               
            end if
            
!     Asselin filter (Asselin 1972)
            do k=1,kb
               do j=1,jm
                  do i=1,im
                     t(i,j,k)=t(i,j,k)
     $                    +.5e0*smoth*(uf(i,j,k)+tb(i,j,k)
     $                    -2.e0*t(i,j,k))
                     s(i,j,k)=s(i,j,k)
     $                    +.5e0*smoth*(vf(i,j,k)+sb(i,j,k)
     $                    -2.e0*s(i,j,k))
                     tb(i,j,k)=t(i,j,k)
                     t(i,j,k)=uf(i,j,k)
                     sb(i,j,k)=s(i,j,k)
                     s(i,j,k)=vf(i,j,k)
                  end do
               end do
            end do

!     Bottom value S.Ohishi 2019.04
            do j=1,jm
               do i=1,im
                  t(i,j,kb)=0.e0
                  tb(i,j,kb)=0.e0
                  s(i,j,kb)=0.e0
                  sb(i,j,kb)=0.e0
               end do
            end do                          

            call dens(s,t,rho)
            
         end if

!     Check NaN
         if(lnan)then
            call check_nan(kb,"tf",uf)
            call check_nan(kb,"sf",vf)
         end if
         
!     write(6,*) 'after t and s'
!     call flush(6)

!     calculate uf and vf
         call advu
         call advv

         call profu
         call profv

         if(ltide) then
            call tide_bcond(3)
         else
            call bcond(3)
         end if

!     S.Ohishi(2019.12)
         if(assim == 1)then
            call uv_iau
         end if

         call exchange3d_mpi(uf(:,:,1:kbm1),im,jm,kbm1)
         call exchange3d_mpi(vf(:,:,1:kbm1),im,jm,kbm1)

!     write(6,*) 'after caluculate uf and vf'
!     call flush(6)

         do j=1,jm
            do i=1,im
               tps(i,j)=0.e0
            end do
         end do

         do k=1,kbm1
            do j=1,jm
               do i=1,im
                  tps(i,j)=tps(i,j)
     $                 +(uf(i,j,k)+ub(i,j,k)-2.e0*u(i,j,k))*dz(k)
               end do
            end do
         end do

         do k=1,kbm1
            do j=1,jm
               do i=1,im
                  u(i,j,k)=u(i,j,k)
     $                 +.5e0*smoth*(uf(i,j,k)+ub(i,j,k)
     $                 -2.e0*u(i,j,k)-tps(i,j))
               end do
            end do
         end do

         do j=1,jm
            do i=1,im
               tps(i,j)=0.e0
            end do
         end do

         do k=1,kbm1
            do j=1,jm
               do i=1,im
                  tps(i,j)=tps(i,j)
     $                 +(vf(i,j,k)+vb(i,j,k)-2.e0*v(i,j,k))*dz(k)
               end do
            end do
         end do

         do k=1,kbm1
            do j=1,jm
               do i=1,im
                  v(i,j,k)=v(i,j,k)
     $                 +.5e0*smoth*(vf(i,j,k)+vb(i,j,k)
     $                 -2.e0*v(i,j,k)-tps(i,j))
               end do
            end do
         end do

         do k=1,kb
            do j=1,jm
               do i=1,im
                  ub(i,j,k)=u(i,j,k)
                  u(i,j,k)=uf(i,j,k)
                  vb(i,j,k)=v(i,j,k)
                  v(i,j,k)=vf(i,j,k)
               end do
            end do
         end do

      end if

      do j=1,jm
         do i=1,im
            egb(i,j)=egf(i,j)
            etb(i,j)=et(i,j)
            et(i,j)=etf(i,j)
            dt(i,j)=h(i,j)+et(i,j)
            utb(i,j)=utf(i,j)
            vtb(i,j)=vtf(i,j)
            vfluxb(i,j)=vfluxf(i,j)
         end do
      end do

      call realvertvl

!     Check NaN
      if(lnan)then
         call check_nan(kb,"uf",uf)
         call check_nan(kb,"vf",vf)
      end if

      return
      end

!_______________________________________________________________________
      subroutine print_section
! print output
      implicit none
      include 'pom.h'
      real atot,darea,dvol,eaver,saver,taver,vtot,tsalt
      integer i,j,k

      if(mod(iint,iprint).eq.0) then

! print time
        if(my_task.eq.master_task) write(6,'(/
     $    ''**********************************************************''
     $    /''time ='',f12.4,'', iint ='',i8,'', iext ='',i8,
     $    '', iprint ='',i8)') time,iint,iext,iprint

! check for errors
        call sum0d_mpi(error_status,master_task)
        call bcast0d_mpi(error_status,master_task)
        if(error_status.ne.0) then
          if(my_task.eq.master_task) write(*,'(/a)')
     $                                       'POM terminated with error'
          call finalize_mpi
          stop
        end if

! local averages
        vtot=0.e0
        atot=0.e0
        taver=0.e0
        saver=0.e0
        eaver=0.e0
        do k=1,kbm1
          do j=1,jm
            do i=1,im
              darea=dx(i,j)*dy(i,j)*fsm(i,j)
              dvol=darea*dt(i,j)*dz(k)
              vtot=vtot+dvol
              taver=taver+tb(i,j,k)*dvol
              saver=saver+sb(i,j,k)*dvol
            end do
          end do
        end do
!
        do j=1,jm
         do i=1,im
            darea=dx(i,j)*dy(i,j)*fsm(i,j)
            atot=atot+darea
            eaver=eaver+et(i,j)*darea
          end do
        end do

!        write(6,*) "vtot:",vtot
!        write(6,*) "taver:",taver
!        write(6,*) "eaver:",eaver

        if(vtot .ne. vtot)then
           write(*,"(/a)") "POM terminated with error: vtot"
           call finalize_mpi
           stop
        else if(vtot==0.)then
           taver=0.
           saver=0.
        else
           taver=taver/vtot
           saver=saver/vtot
        end if

        if(atot==0.)then
           eaver=0.
        else
           eaver=eaver/atot
        end if

        tsalt=(saver+sbias)*vtot

! print averages
! global averages requiere to transfer high amounts of data between
! processor - therefore, only local average for master_task is printed
        if(my_task.eq.master_task) write(6,'(/''vtot = '',e16.7,
     $    ''   atot = '',e16.7,''  eaver ='',e16.7/''taver ='',e16.7,
     $    ''   saver ='',e16.7,''  tsalt ='',e16.7)')
     $    vtot,atot,eaver,taver,saver,tsalt

      end if

      return
      end

!_______________________________________________________________________
      subroutine ts_nudging
!     Temperature & Salinity nudges toward tsclim monthly climatology
!     Created by S.Ohishi 2020.04

      implicit none
      include 'pom.h'

      integer i,j,k
      
      real time_scale

!     Temperature
      do k=1,kb-1

         if(k == 1)then
            time_scale=ts_nudge
         else
            time_scale=ti_nudge
         end if

         if(time_scale > 0.)then
            time_scale=dti2/(time_scale*86400.)
         else
            tnudge(:,:,k)=0.
            cycle
         end if

         do j=1,jm
            do i=1,im
               tnudge(i,j,k)=(tclimm(i,j,k)-uf(i,j,k))*time_scale
               uf(i,j,k)=uf(i,j,k)+tnudge(i,j,k)
            end do
         end do
      end do

!     Salinity
      do k=1,kb-1

         if(k == 1)then
            time_scale=ss_nudge
         else
            time_scale=si_nudge
         end if

         if(time_scale > 0.)then
            time_scale=dti2/(time_scale*86400.)
         else
            snudge(:,:,k)=0.
            cycle
         end if

         do j=1,jm
            do i=1,im
               snudge(i,j,k)=(sclimm(i,j,k)-vf(i,j,k))*time_scale
               vf(i,j,k)=vf(i,j,k)+snudge(i,j,k)
            end do
         end do
      end do

      return
      end

!_______________________________________________________________________
      subroutine ts_roff
!     created by S.Ohishi 2020.02 

      implicit none
      include 'pom.h'
      
      integer i,j,k
      
      real,parameter :: tmin=-1.8e0,tmax=40.e0
      real,parameter :: smin=30.e0,smax=40.e0

      do k=1,kb-1
         do j=1,jm
            do i=1,im

               if(uf(i,j,k) == 0.e0) cycle

               if(uf(i,j,k) < tmin)then
                  troff(i,j,k)=uf(i,j,k)-tmin
                  uf(i,j,k)=tmin
               else if(tmax < uf(i,j,k))then
                  troff(i,j,k)=tmax-uf(i,j,k)
                  uf(i,j,k)=tmax
               end if        

            end do
         end do
      end do

      do k=1,kb-1
         do j=1,jm
            do i=1,im

               if(vf(i,j,k) == 0.e0) cycle

               if(vf(i,j,k) < smin)then
                  sroff(i,j,k)=vf(i,j,k)-smin
                  vf(i,j,k)=smin
               else if(smax < vf(i,j,k))then
                  sroff(i,j,k)=smax-vf(i,j,k)
                  vf(i,j,k)=smax     
               endif

            end do
         end do
      end do

      return
      end subroutine
!_______________________________________________________________________
      subroutine check_velocity
! check if velocity condition is violated
      implicit none
      include 'pom.h'
      real vamax,atot,darea,dvol,eaver,saver,taver,vtot,tsalt
      integer i,j,k
      integer imax,jmax

      vamax=0.e0

      do j=1,jm
        do i=1,im
          if(abs(vaf(i,j)).ge.vamax) then
            vamax=abs(vaf(i,j))
            imax=i
            jmax=j
          end if
        end do
      end do

      if(vamax.gt.vmaxl) then
        if(my_task.eq.master_task.and.error_status.eq.0) write(6,'(/
     $    ''Error: velocity condition violated''/''time ='',f9.4,
     $    '', iint ='',i8,'', iext ='',i8,'', iprint ='',i8,/
     $    ''vamax ='',e12.3,''   imax,jmax ='',2i5)')
     $    time,iint,iext,iprint,vamax,imax,jmax
        error_status=1
      end if

      return
      end

!_____________________________________________________________________________

      subroutine ts_iau

!     IAU method 
!     Relaxation of TS toward analysis TS

      implicit none
      include 'pom.h'

      integer i,j,k

      do k=1,kb
         do j=1,jm
            do i=1,im
               uf(i,j,k)=uf(i,j,k)+2.*t_iau(i,j,k)/real(iend)
               vf(i,j,k)=vf(i,j,k)+2.*s_iau(i,j,k)/real(iend)
!               uf(i,j,k)=uf(i,j,k)+t_iau(i,j,k)/real(iend)
!               vf(i,j,k)=vf(i,j,k)+s_iau(i,j,k)/real(iend)
            end do
         end do
      end do

      return
      end
      
!________________________________________________________________________________

      subroutine uv_iau

!     IAU method 
!     Relaxation of TS toward analysis UV

      implicit none
      include 'pom.h'

      integer i,j,k

      do k=1,kb
         do j=1,jm
            do i=1,im
               uf(i,j,k)=uf(i,j,k)+2.*u_iau(i,j,k)/real(iend)
               vf(i,j,k)=vf(i,j,k)+2.*v_iau(i,j,k)/real(iend)
!               uf(i,j,k)=uf(i,j,k)+u_iau(i,j,k)/real(iend)
!               vf(i,j,k)=vf(i,j,k)+v_iau(i,j,k)/real(iend)
            end do
         end do
      end do

      return
      end
      
!____________________________________________________________________________________

      subroutine check_nan(ndep,var,dat)

      implicit none
      include 'pom.h'

      integer i,j,idep
      integer ndep

      real dat(im,jm,ndep)
      
      character(10) var,cint

      do idep=1,ndep
         do j=2,jmm1
            do i=2,imm1
               if(dat(i,j,idep) .ne. dat(i,j,idep))then
                  write(cint,'(i10.10)') iint
                  write(*,"(a,3i6)") "Position: ",i,j,idep
                  write(*,"(a)") "Time step: "//cint
                  write(*,'(i6,e10.2)') 
     $                 0,dat(i,j,idep)
                  write(*,'(i6,3e10.2)') 
     $                 1,uaf(i,j),adx2d(i,j),advua(i,j)
                  write(*,'(i6,3e10.2)') 
     $                 2,aru(i,j)
                  write(*,'(i6,4e10.2)') 
     $                 3,cor(i,j),d(i,j),va(i,j+1),va(i,j)
                  write(*,'(i6,4e10.2)') 
     $                 4,cor(i-1,j),d(i-1,j),va(i-1,j+1),va(i-1,j)
                  write(*,'(i6,2e10.2)') 
     $                 5,dy(i,j),dy(i-1,j)
                  write(*,'(i6,2e10.2)') 
     $                 6,d(i,j),d(i-1,j)
                  write(*,'(i6,2e10.2)') 
     $                 7,el(i,j),el(i-1,j)
                  write(*,'(i6,2e10.2)') 
     $                 8,elb(i,j),elb(i-1,j)
                  write(*,'(i6,2e10.2)') 
     $                 9,elf(i,j),elf(i-1,j)
                  write(*,'(i6,2e10.2)') 
     $                 10,e_atmos(i,j),e_atmos(i-1,j)
                  write(*,'(i6,4e10.2)') 
     $                 11,drx2d(i,j),aru(i,j),wusurf(i,j),wubot(i,j)
                  write(*,'(i6,2e10.2)') 
     $                 12,windu(i,j),windv(i,j)
                  write(*,'(i6,2e10.2)') 
     $                 13,wusurf(i,j),wvsurf(i,j)
                  write(*,"(a)") 
     $                 "POM terminated with error: "//trim(var)
                  call finalize_mpi
                  stop                  
               end if
            end do
         end do
      end do

      return
      end
