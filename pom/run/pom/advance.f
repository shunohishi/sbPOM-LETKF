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

! get time
      write(*,*) "get time"
      write(6,*) "get time"
      write(98,*) "get time"
      write(99,*) "get time"
      call get_time

! set time dependent surface boundary conditions
      call surface_forcing

! set time dependent lateral boundary conditions
      call lateral_forcing

! set time dependent reference temperature and salinity
      call set_tsdata

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
!      write(6,*) 'after mode_external'
!      call flush(6)

! debug
!      iglb=(im_global-2)/2
!      jglb=(jm_global-2)/2
!      nproc=int((iglb-2)/(im-2))+int((jglb-2)/(jm-2))*nproc_x

! debug
!      if(my_task.eq.nproc .and. mod(iint-1,1).eq.0) then 
!        iloc=iglb-mod(nproc,nproc_x)*(im-2)
!        jloc=jglb-(nproc/nproc_x)*(jm-2)
!        write(99,*) 'in advance iint ',iint
!        write(99,*) 'iloc jloc nproc ',iloc,jloc,nproc
!        write(99,*) 'iglb jglb ',iglb,jglb
!        write(99,*) 'u ',u(iloc,jloc,1:kb-1)
!        write(99,*) 'v ',v(iloc,jloc,1:kb-1)
!        write(99,*) 't ',t(iloc,jloc,1:kb-1)
!        write(99,*) 's ',s(iloc,jloc,1:kb-1)
!        write(99,*) 'el elb ',el(iloc,jloc),elb(iloc,jloc)
!        write(99,*) 'egf ',egf(iloc,jloc)
!        write(99,*) 'uab vab ',uab(iloc,jloc),vab(iloc,jloc)
!        write(99,*) 'ua va ',ua(iloc,jloc),va(iloc,jloc)
!        write(99,*) 'utf/d vtf/d '
!     $             ,utf(iloc,jloc)/dt(iloc,jloc)
!     $             ,vtf(iloc,jloc)/dt(iloc,jloc)
!        call flush(99)
!      end if

!      uvmax=0.
!      uvmin=99.
!      tmax=0.
!      tmin=99.
!      smax=0.
!      smin=99.
!      elmax=0.
!      elmin=99.

!      do i=1,im
!      do j=1,jm
!      do k=1,kb-1
!        if(u(i,j,k).gt.uvmax .and. fsm(i,j).eq.1.) then
!          uvmax=u(i,j,k)
!          iuvmax=i
!          juvmax=j
!          kuvmax=k
!        end if
!        if(u(i,j,k).lt.uvmin .and. fsm(i,j).eq.1.) then
!          uvmin=u(i,j,k)
!          iuvmin=i
!          juvmin=j
!          kuvmin=k
!        end if
!        if(v(i,j,k).gt.uvmax .and. fsm(i,j).eq.1.) then
!          uvmax=v(i,j,k)
!          iuvmax=i
!          juvmax=j
!          kuvmax=k
!        end if
!        if(v(i,j,k).lt.uvmin .and. fsm(i,j).eq.1.) then
!          uvmin=v(i,j,k)
!          iuvmin=i
!          juvmin=j
!          kuvmin=k
!        end if
!        if(t(i,j,k).gt.tmax .and. fsm(i,j).eq.1.) then
!          tmax=t(i,j,k)
!          itmax=i
!          jtmax=j
!          ktmax=k
!        end if
!        if(t(i,j,k).lt.tmin .and. fsm(i,j).eq.1.) then
!          tmin=t(i,j,k)
!          itmin=i
!          jtmin=j
!          ktmin=k
!        end if
!        if(s(i,j,k).gt.smax .and. fsm(i,j).eq.1.) then
!          smax=s(i,j,k)
!          ismax=i
!          jsmax=j
!          ksmax=k
!        end if
!        if(s(i,j,k).lt.smin .and. fsm(i,j).eq.1.) then
!          smin=s(i,j,k)
!          ismin=i
!          jsmin=j
!          ksmin=k
!        end if
!      end do
!      end do
!      end do
!      do i=1,im
!      do j=1,jm
!        if(el(i,j).gt.elmax .and. fsm(i,j).eq.1.) then
!          elmax=el(i,j)
!          ielmax=i
!          jelmax=j
!        end if
!        if(el(i,j).lt.elmin .and. fsm(i,j).eq.1.) then
!          elmin=el(i,j)
!          ielmin=i
!          jelmin=j
!        end if
!      end do
!      end do
!      write(99,*) 'in advance iint my_task ',iint,my_task
!      write(99,*) 'uvmax ',my_task,uvmax,iuvmax,juvmax,kuvmax
!      write(99,*) 'uvmin ',my_task,uvmin,iuvmin,juvmin,kuvmin
!      write(99,*) 'tmax ',my_task,tmax,itmax,jtmax,ktmax
!      write(99,*) 'tmin ',my_task,tmin,itmin,jtmin,ktmin
!      write(99,*) 'smax ',my_task,smax,ismax,jsmax,ksmax
!      write(99,*) 'smin ',my_task,smin,ismin,jsmin,ksmin
!      write(99,*) 'elmax ',my_task,elmax,ielmax,jelmax
!      write(99,*) 'elmin ',my_task,elmin,ielmin,jelmin

!      call flush(99)

! internal (3-D) mode calculation
! diag

!      write(6,*) 'before mode_internal'
!      call flush(6)

      call mode_internal
!      call mode_internal_diag

!      write(6,*) 'after mode_internal'
!      call flush(6)

! print section
      call print_section

! write output
      if(netcdf_file.ne.'nonetcdf' .and. mod(iint,iprint).eq.0)
     $                                         call write_output_netcdf

! write restart
      if(mod(iint,irestart).eq.0) call write_restart_netcdf

! check CFL condition
      call check_velocity

      return
      end

!_______________________________________________________________________
      subroutine get_time
! return the model time
      implicit none
      include 'pom.h'

      time=dti*float(iint)/86400.e0+time0
      if(iint.ge.iswtch) iprint=nint(prtd2*24.e0*3600.e0/dti)
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
        write(6,*) 'tsclim ratio', ratio 
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

        if (npg.eq.1) then
          call baropg
        else if (npg.eq.2) then
          call baropg_mcc
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

! kii2b
!        call advave(tps)
        call advave

        do j=1,jm
          do i=1,im
            adx2d(i,j)=adx2d(i,j)-advua(i,j)
            ady2d(i,j)=ady2d(i,j)-advva(i,j)
          end do
        end do

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

!      call bcond(1)
      call tide_bcond(1)

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

!      call bcond(2)
      call tide_bcond(2)

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

      return
      end

!_______________________________________________________________________
      subroutine mode_internal
! calculate the internal (3-D) mode
      implicit none
      include 'pom.h'
      integer i,j,k
      real dxr,dxl,dyt,dyb

      if((iint.ne.1.or.time0.ne.0.e0).and.mode.ne.2) then

! adjust u(z) and v(z) such that depth average of (u,v) = (ua,va)
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

! calculate w from u, v, dt (h+et), etf and etb
! kii2b
!        call vertvl(a,c)
        call vertvl

!        call bcond(5)
        call tide_bcond(5)

        call exchange3d_mpi(w,im,jm,kb)

!        write(6,*) 'after w'
!        call flush(6)

! set uf and vf to zero
        do k=1,kb
          do j=1,jm
            do i=1,im
              uf(i,j,k)=0.e0
              vf(i,j,k)=0.e0
            end do
          end do
        end do

! calculate q2f and q2lf using uf, vf, a and c as temporary variables
! kii2b
!        call advq(q2b,q2,uf,a,c)
!        call advq(q2lb,q2l,vf,a,c)
!        call profq(a,c,tps,dtef)
        call advq(q2b,q2,uf)
        call advq(q2lb,q2l,vf)

!        write(6,*) 'after advq'
!        call flush(6)

        call profq

!        write(6,*) 'after profq'
!        call flush(6)

!        call bcond(6)
        call tide_bcond(6)

!        write(6,*) 'after bcond(6)'
!        call flush(6)

        call exchange3d_mpi(uf(:,:,2:kbm1),im,jm,kbm2)
        call exchange3d_mpi(vf(:,:,2:kbm1),im,jm,kbm2)

!        write(6,*) 'after q2f and q2lf'
!        call flush(6)

        do k=1,kb
          do j=1,jm
            do i=1,im
              q2(i,j,k)=q2(i,j,k)
     $                   +.5e0*smoth*(uf(i,j,k)+q2b(i,j,k)
     $                                -2.e0*q2(i,j,k))
              q2l(i,j,k)=q2l(i,j,k)
     $                   +.5e0*smoth*(vf(i,j,k)+q2lb(i,j,k)
     $                                -2.e0*q2l(i,j,k))
              q2b(i,j,k)=q2(i,j,k)
              q2(i,j,k)=uf(i,j,k)
              q2lb(i,j,k)=q2l(i,j,k)
              q2l(i,j,k)=vf(i,j,k)
            end do
          end do
        end do

! calculate tf and sf using uf, vf, a and c as temporary variables
        if(mode.ne.4) then
          if(nadv.eq.1) then
! kii2b
!            call advt1(tb,t,tclim,uf,a,c)
!            call advt1(sb,s,sclim,vf,a,c)
            call advt1(tb,t,tclim,uf)
            call advt1(sb,s,sclim,vf)
          else if(nadv.eq.2) then
! kii2b
!            call advt2(tb,t,tclim,uf,a,c)
!            call advt2(sb,s,sclim,vf,a,c)
            call advt2(tb,t,tclim,uf)
            call advt2(sb,s,sclim,vf)
          else
            error_status=1
            write(6,'(/''Error: invalid value for nadv'')')
          end if

! kii2b
!          call proft(uf,wtsurf,tsurf,nbct,tps)
!          call proft(vf,wssurf,ssurf,nbcs,tps)
          call proft(uf,wtsurf,tsurf,nbct)
          call proft(vf,wssurf,ssurf,nbcs)

!          call bcond(4)
          call tide_bcond(4)
          if(lfrs) call bcond_frs(4)

          call exchange3d_mpi(uf(:,:,1:kbm1),im,jm,kbm1)
          call exchange3d_mpi(vf(:,:,1:kbm1),im,jm,kbm1)

          do k=1,kb
            do j=1,jm
              do i=1,im
                t(i,j,k)=t(i,j,k)
     $                    +.5e0*smoth*(uf(i,j,k)+tb(i,j,k)
     $                                 -2.e0*t(i,j,k))
                s(i,j,k)=s(i,j,k)
     $                    +.5e0*smoth*(vf(i,j,k)+sb(i,j,k)
     $                                 -2.e0*s(i,j,k))
                tb(i,j,k)=t(i,j,k)
                t(i,j,k)=uf(i,j,k)
                sb(i,j,k)=s(i,j,k)
                s(i,j,k)=vf(i,j,k)
              end do
            end do
          end do

          call dens(s,t,rho)

        end if

!        write(6,*) 'after t and s'
!        call flush(6)

! calculate uf and vf
        call advu
        call advv
        call profu
        call profv

!        call bcond(3)
        call tide_bcond(3)

        call exchange3d_mpi(uf(:,:,1:kbm1),im,jm,kbm1)
        call exchange3d_mpi(vf(:,:,1:kbm1),im,jm,kbm1)

!        write(6,*) 'after caluculate uf and vf'
!        call flush(6)

        do j=1,jm
          do i=1,im
            tps(i,j)=0.e0
          end do
        end do

        do k=1,kbm1
          do j=1,jm
            do i=1,im
              tps(i,j)=tps(i,j)
     $                  +(uf(i,j,k)+ub(i,j,k)-2.e0*u(i,j,k))*dz(k)
            end do
          end do
        end do

        do k=1,kbm1
          do j=1,jm
            do i=1,im
              u(i,j,k)=u(i,j,k)
     $                  +.5e0*smoth*(uf(i,j,k)+ub(i,j,k)
     $                               -2.e0*u(i,j,k)-tps(i,j))
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
     $                  +(vf(i,j,k)+vb(i,j,k)-2.e0*v(i,j,k))*dz(k)
            end do
          end do
        end do

        do k=1,kbm1
          do j=1,jm
            do i=1,im
              v(i,j,k)=v(i,j,k)
     $                  +.5e0*smoth*(vf(i,j,k)+vb(i,j,k)
     $                               -2.e0*v(i,j,k)-tps(i,j))
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

! kii2b
! calculate real w as wr
!      do k=1,kb
!        do j=1,jm
!          do i=1,im
!            wr(i,j,k)=0.
!          end do
!        end do
!      end do
!
!      do k=1,kbm1
!        do j=1,jm
!          do i=1,im
!            tps(i,j)=zz(k)*dt(i,j) + et(i,j)
!          end do
!        end do
!        do j=2,jmm1
!          do i=2,imm1
!            dxr=2.0/(dx(i+1,j)+dx(i,j))
!            dxl=2.0/(dx(i,j)+dx(i-1,j))
!            dyt=2.0/(dy(i,j+1)+dy(i,j))
!            dyb=2.0/(dy(i,j)+dy(i,j-1))
!            wr(i,j,k)=0.5*(w(i,j,k)+w(i,j,k+1))+0.5*
!     $                (u(i+1,j,k)*(tps(i+1,j)-tps(i,j))*dxr+
!     $                 u(i,j,k)*(tps(i,j)-tps(i-1,j))*dxl+
!     $                 v(i,j+1,k)*(tps(i,j+1)-tps(i,j))*dyt+
!     $                 v(i,j,k)*(tps(i,j)-tps(i,j-1))*dyb)
!     $                +(1.0+zz(k))*(etf(i,j)-etb(i,j))/dti2
!          end do
!        end do
!      end do
!
!      call exchange3d_mpi(wr(:,:,1:kbm1),im,jm,kbm1)
!
!      do k=1,kb
!        do i=1,im
!          if(n_south.eq.-1) wr(i,1,k)=wr(i,2,k)
!          if(n_north.eq.-1) wr(i,jm,k)=wr(i,jmm1,k)
!        end do
!      end do
!      do k=1,kb
!        do j=1,jm
!          if(n_west.eq.-1) wr(1,j,k)=wr(2,j,k)
!          if(n_east.eq.-1) wr(im,j,k)=wr(imm1,j,k)
!        end do
!      end do
!
!      do k=1,kbm1
!        do j=1,jm
!          do i=1,im
!            wr(i,j,k)=fsm(i,j)*wr(i,j,k)
!          end do
!        end do
!      end do
! kii2b
      call realvertvl

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
     $    /''time ='',f9.4,'', iint ='',i8,'', iext ='',i8,
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

        do j=1,jm
          do i=1,im
            darea=dx(i,j)*dy(i,j)*fsm(i,j)
            atot=atot+darea
            eaver=eaver+et(i,j)*darea
          end do
        end do

        taver=taver/vtot
        saver=saver/vtot
        eaver=eaver/atot
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

