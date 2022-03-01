! initialize.f

! initialize POM: define constant, read initial values, grid and initial
! conditions

!_______________________________________________________________________
      subroutine initialize
! initialize POM
      implicit none
      include 'pom.h'
      integer i,j,k
      integer year,month,day,julday
      character*4 cyear
      character*2 cmonth,cday

! initialize the MPI execution environment and create communicator for
!internal POM communications
      call initialize_mpi

! distribute the model domain across processors
      call distribute_mpi

! read input values and define constants
      call read_input

! initialize arrays for safety (may be overwritten)
      call initialize_arrays

! read in grid data
      call read_grid

! read initial and lateral boundary conditions
      call initial_conditions

! update initial conditions and set the remaining conditions
      call update_initial

! calculate the bottom friction coefficient
      call bottom_friction

! read restart data from a previous run
      if(nread_rst.ne.0) call read_restart_netcdf
      if(nread_rst.eq.2) call restart_tsuv

! read iau data from letkf analysis (S.Ohishi 2018.12)
      if(assim .eq. 1) call read_iau_netcdf
      
! write grid and initial conditions
!      if(netcdf_file.ne.'nonetcdf') call write_output_netcdf
! debug
!      call write_restart_netcdf
!      call finalize_mpi
!      stop

! set start julian day
      cyear=time_start(1:4)
      read(cyear,*) year
      cmonth=time_start(6:7)
      read(cmonth,*) month
      cday=time_start(9:10)
      read(cday,*) day
      julday_start=julday(year,month,day)

! check for errors
      call sum0d_mpi(error_status,master_task)
      call bcast0d_mpi(error_status,master_task)
!      if(error_status.ne.0) then
!        if(my_task.eq.master_task) write(*,'(/a)')
!     $                                       'POM terminated with error'
!        call finalize_mpi
!        stop
!      end if
      if(my_task.eq.master_task) write(*,'(/a)') 'End of initialization'

      return
      end

!_______________________________________________________________________
      subroutine read_input
! read input values and defines constants
      implicit none
      include 'pom.h'
      namelist/pom_nml/ title,netcdf_file,mode,assim,
     $     nadv,nitera,sw,npg,dte,
     $     isplit,time_start,nread_rst,read_rst_file,read_iau_file,
     $     write_rst,write_rst_file,
     $     write_ens,write_ens_file,
     $     budget,days,prtd1,prtd2,
     $     swtch,
     $     ts_nudge,ti_nudge,ss_nudge,si_nudge

! read input namelist
      open(idevinout,file='pom.nml',status='old')
      read(idevinout,nml=pom_nml)
      close(idevinout)
  
! Input of filenames and constants

! Logical for inertial ramp (.true. if inertial ramp to be applied
! to wind stress and baroclinic forcing, otherwise .false.)
      lramp=.false.

! Reference density (recommended values: 1025 for seawater,
! 1000 for freswater; S.I. units):
      rhoref=1025.e0

! Temperature bias (deg. C)
      tbias=0.e0

! Salinity bias
      sbias=0.e0

! gravity constant (S.I. units)
      grav=9.806e0

! von Karman's constant
      kappa=0.4e0

! Bottom roughness (metres)
      z0b=.01e0

! Minimum bottom friction coeff.
      cbcmin=.0025e0

! Maximum bottom friction coeff.
      cbcmax=1.e0

! Smagorinsky diffusivity coeff.
      horcon=0.1e0
      aamadd = 0.

! Inverse horizontal turbulent Prandtl number (ah/am; dimensionless):
! NOTE that tprni=0.e0 yields zero horizontal diffusivity!
      tprni=.2e0

! Background viscosity used in subroutines profq, proft, profu and
! profv (S.I. units):
      umol=2.e-5

! Maximum magnitude of vaf (used in check that essentially tests
! for CFL violation):
      vmaxl=10.e0

! Maximum allowable value of:
!   <difference of depths>/<sum of depths>
! for two adjacent cells (dimensionless). This is used in subroutine
! slpmax. If >= 1, then slpmax is not applied:
      slmax=2.e0

! Water type, used in subroutine proft.
!    ntp    Jerlov water type
!     1            i
!     2            ia
!     3            ib
!     4            ii
!     5            iii
      ntp=2

! Surface temperature boundary condition, used in subroutine proft:
!    nbct   prescribed    prescribed   short wave
!           temperature      flux      penetration
!     1        no           yes           no
!     2        no           yes           yes
!     3        yes          no            no
!     4        yes          no            yes
      nbct=2

! Surface salinity boundary condition, used in subroutine proft:
!    nbcs   prescribed    prescribed
!            salinity      flux
!     1        no           yes
!     3        yes          no
! NOTE that only 1 and 3 are allowed for salinity.
      nbcs=1

!     2018.06 added by S.Ohishi
!     2020.04 modified by S.Ohishi
! Surface salinity forcing:
!     lssf  prescribed
!             E-P-R
!       1      yes
!       2      no
      lssf=1

! 2018.08 added by S.Ohishi
! lateral boudary forcing:
!     llbc   monthly climatology     daily  
!      1           yes                no
!      2           no                 yes
      llbc=1

!2018.08 added by S.Ohishi
! output:
!    ldave   daily average      instanteneous
!      1         yes                 no
!      2         no                  yes
      ldave=1

! for wind stress formula by Mellor and Blumberg (2004)
!      lmbws=.true.

!2018.09 S.Ohishi
! turbulent heat flux & wind stress:
! #2-5: use https://github.com/brodeau/aerobulk, Brodeau et al. (2017) 
!  lthf_ws
!   1   Kondo(1975) & Mellor and Blumberg (2004)
!   2   coare 3.0
!   3   coare 3.5
!   4   ncar
!   5   ecmwf
      lthf_ws=3
      
! Step interval during which external (2-D) mode advective terms are
! not updated (dimensionless):
      ispadv=5

! Constant in temporal filter used to prevent solution splitting
! (dimensionless):
      smoth=0.10e0

! Weight used for surface slope term in external (2-D) dynamic
! equation (a value of alpha = 0.e0 is perfectly acceptable, but the
! value, alpha=.225e0 permits a longer time step):
      alpha=0.225e0

! Initial value of aam:
      aam_init=500.e0

! Flow Relaxation Scheme
      lfrs=.false.

! for surface heat flux & wind stress relative to ocean current
      lrtvf=.true.

! for MYNNF
      lmynnf=.true.

! seto2 experiments
! for tide
      ltide=.false.

! Surface wave breaking (Mellor and Blumberg, 2004)
      lwbrk=.false.
! for internal wave breaking parameterization
      liwbrk=.false.

! S.Ohishi 2020.04
!for temperature/salinity round off
      lroff=.true.

! S.Ohishi 2020.10
! Check NaN
      lnan=.true.

! End of input of constants

! calculate some constants
      small=1.e-9           ! Small value
      pi=atan(1.e0)*4.e0    ! PI

      dti=dte*float(isplit)
      dte2=dte*2
      dti2=dti*2

      iend=max0(nint(days*24.e0*3600.e0/dti),2)
      iprint=nint(prtd1*24.e0*3600.e0/dti)
      iswtch=nint(swtch*24.e0*3600.e0/dti)
      irestart=nint(write_rst*24.e0*3600.e0/dti)
      iens=nint(write_ens*24.e0*3600.e0/dti)

      ispi=1.e0/float(isplit)
      isp2i=1.e0/(2.e0*float(isplit))

! initialise time
      time0=0.e0
      time=0.e0

! print initial summary
      if(my_task.eq.master_task) then
        write(6,'(/'' title      = '',a40)') title
        write(6,'(/'' mode       = '',i10)') mode
        write(6,'(/'' assim       = '',i10)') assim
        write(6,'('' nadv       = '',i10)') nadv
        write(6,'('' nitera     = '',i10)') nitera
        write(6,'('' sw         = '',f10.4)') sw
        write(6,'('' npg         = '',i10)') npg
        write(6,'('' nread_rst  = '',i10)') nread_rst
        write(6,'('' write_rst  = '',f10.4)') write_rst
        write(6,'('' irestart   = '',i10)') irestart
        write(6,'('' write_ens  = '',f10.4)') write_ens
        write(6,'('' iens       = '',i10)') iens
        write(6,'('' dte        = '',f10.2)') dte
        write(6,'('' dti        = '',f10.1)') dti
        write(6,'('' isplit     = '',i10)') isplit
        write(6,'('' time_start = '',a26)') time_start
        write(6,'('' days       = '',f10.4)') days
        write(6,'('' iend       = '',i10)') iend
        write(6,'('' prtd1      = '',f10.4)') prtd1
        write(6,'('' iprint     = '',i10)') iprint
        write(6,'('' prtd2      = '',f10.4)') prtd2
        write(6,'('' swtch      = '',f10.2)') swtch
        write(6,'('' iswtch     = '',i10)') iswtch
        write(6,'('' lramp      = '',l10)') lramp
        write(6,'('' lfrs       = '',l10)') lfrs
        write(6,'('' lrtvf      = '',l10)') lrtvf
        write(6,'('' lmynnf     = '',l10)') lmynnf
        write(6,'('' rhoref     = '',f10.3)') rhoref
        write(6,'('' tbias      = '',f10.3)') tbias
        write(6,'('' sbias      = '',f10.3)') sbias
        write(6,'('' grav       = '',f10.4)') grav
        write(6,'('' kappa      = '',f10.4)') kappa
        write(6,'('' z0b        = '',f10.6)') z0b
        write(6,'('' cbcmin     = '',f10.6)') cbcmin
        write(6,'('' cbcmax     = '',f10.6)') cbcmax
        write(6,'('' horcon     = '',f10.3)') horcon
        write(6,'('' aamadd     = '',f10.3)') aamadd
        write(6,'('' tprni      = '',f10.4)') tprni
        write(6,'('' umol       = '',f10.4)') umol
        write(6,'('' vmaxl      = '',f10.4)') vmaxl
        write(6,'('' slmax      = '',f10.4)') slmax
        write(6,'('' ntp        = '',i10)') ntp
        write(6,'('' nbct       = '',i10)') nbct
        write(6,'('' nbcs       = '',i10)') nbcs
        write(6,'('' budget     = '',i10)') budget       
        write(6,'('' ldave      = '',i10)') ldave       
        write(6,'('' llbc       = '',i10)') llbc       
        write(6,'('' lssf       = '',i10)') lssf        
        write(6,'('' lthf_ws    = '',i10)') lthf_ws
        write(6,'('' ispadv     = '',i10)') ispadv
        write(6,'('' smoth      = '',f10.4)') smoth
        write(6,'('' alpha      = '',f10.4)') alpha
        write(6,'('' ts_nudge   = '',f10.4)') ts_nudge
        write(6,'('' ti_nudge   = '',f10.4)') ti_nudge
        write(6,'('' ss_nudge   = '',f10.4)') ss_nudge
        write(6,'('' si_nudge   = '',f10.4)') si_nudge
! seto2 experiments
        write(6,'('' ltide      = '',l10)') ltide
        write(6,'('' lwbrk      = '',l10)') lwbrk
        write(6,'('' liwbrk     = '',l10)') liwbrk
        write(6,'('' lroff      = '',l10)') lroff
        write(6,'('' lnan       = '',l10)') lnan

      end if

      return
      end

!_______________________________________________________________________
      subroutine initialize_arrays
! initialize arrays for safety
      implicit none
      include 'pom.h'
      integer i,j,k

! boundary arrays
      do i=1,im
        vabn(i)=0.e0
        vabs(i)=0.e0
        eln(i)=0.e0
        els(i)=0.e0
        do k=1,kb
          vbn(i,k)=0.e0
          vbs(i,k)=0.e0
          tbn(i,k)=0.e0
          tbs(i,k)=0.e0
          sbn(i,k)=0.e0
          sbs(i,k)=0.e0
        end do
      end do

      do j=1,jm
        uabe(j)=0.e0
        uabw(j)=0.e0
        ele(j)=0.e0
        elw(j)=0.e0
        do k=1,kb
          ube(j,k)=0.e0
          ubw(j,k)=0.e0
          tbe(j,k)=0.e0
          tbw(j,k)=0.e0
          sbe(j,k)=0.e0
          sbw(j,k)=0.e0
        end do
      end do

! 2-D and 3-D arrays
      do j=1,jm
        do i=1,im
          uab(i,j)=0.e0
          vab(i,j)=0.e0
          elb(i,j)=0.e0
          etb(i,j)=0.e0
          e_atmos(i,j)=0.e0
          vfluxb(i,j)=0.e0
          vfluxf(i,j)=0.e0
          wusurf(i,j)=0.e0
          wvsurf(i,j)=0.e0
          wtsurf(i,j)=0.e0
          wssurf(i,j)=0.e0
          swrad(i,j)=0.e0
          drx2d(i,j)=0.e0
          dry2d(i,j)=0.e0
        end do
      end do

      do k=1,kbm1
        do j=1,jm
          do i=1,im
            ub(i,j,k)=0.e0
            vb(i,j,k)=0.e0
          end do
        end do
      end do

      return
      end

!_______________________________________________________________________
      subroutine read_grid
! set up vertical and horizontal grid, topography, areas and masks
      implicit none
      include 'pom.h'
      integer i,j,k
      real deg2rad

! degrees to radians
      deg2rad=pi/180.

! read grid
      call read_grid_netcdf

! derived vertical grid variables
      do k=1,kb-1
        dz(k)=z(k)-z(k+1)
        dzz(k)=zz(k)-zz(k+1)
      end do
      dz(kb)=0.
      dzz(kb)=0.

! print vertical grid information
      if(my_task.eq.master_task) then
        write(6,'(/2x,a,7x,a,9x,a,9x,a,9x,a)') 'k','z','zz','dz','dzz'
        do k=1,kb
          write(6,'(1x,i5,4f10.3)') k,z(k),zz(k),dz(k),dzz(k)
        end do
      end if

! set up Coriolis parameter
        do j=1,jm
          do i=1,im
            cor(i,j)=2.*7.29e-5*sin(north_e(i,j)*deg2rad)
          end do
        end do

! inertial period for temporal filter
      period=(2.e0*pi)/abs(cor(im/2,jm/2))/86400.e0
!      period=60.

! calculate areas of "t" and "s" cells
      do j=1,jm
        do i=1,im
          art(i,j)=dx(i,j)*dy(i,j)
        end do
      end do

! calculate areas of "u" and "v" cells
      do j=2,jm
        do i=2,im
          aru(i,j)=.25e0*(dx(i,j)+dx(i-1,j))*(dy(i,j)+dy(i-1,j))
          arv(i,j)=.25e0*(dx(i,j)+dx(i,j-1))*(dy(i,j)+dy(i,j-1))
        end do
      end do
      call exchange2d_mpi(aru,im,jm)
      call exchange2d_mpi(arv,im,jm)

      if (n_west.eq.-1) then
        do j=1,jm
          aru(1,j)=aru(2,j)
          arv(1,j)=arv(2,j)
        end do
      end if

      if (n_south.eq.-1) then
        do i=1,im
          aru(i,1)=aru(i,2)
          arv(i,1)=arv(i,2)
        end do
      end if

      do i=1,im
        do j=1,jm
          d(i,j)=h(i,j)+el(i,j)
          dt(i,j)=h(i,j)+et(i,j)
          end do
      end do

      return
      end

!_______________________________________________________________________
      subroutine initial_conditions
! set up initial and lateral boundary conditions
      implicit none
      include 'pom.h'
      integer i,j,k
      real elejmid,elwjmid
      real tb2(im,jm,kb),sb2(im,jm,kb)

      call read_initial_tsuv_netcdf(kb,tb,sb,ub,vb)
      do k=1,kbm1
        do j=1,jm
          do i=1,im
            tb(i,j,k)=tb(i,j,k)-tbias
            sb(i,j,k)=sb(i,j,k)-sbias
          end do
        end do
      end do

! density
      call dens(sb,tb,rho)

      do k=1,kbm1
        do j=1,jm
          do i=1,im
            t(i,j,k)=tb(i,j,k)
            s(i,j,k)=sb(i,j,k)
            tclim(i,j,k)=tb(i,j,k)
            sclim(i,j,k)=sb(i,j,k)
          end do
        end do
      end do

! read climatological/mean temperature and salinity 
      call read_tsclim_netcdf
      do k=1,kbm1
        do j=1,jm
          do i=1,im
            tclim(i,j,k)=tclim(i,j,k)-tbias
            sclim(i,j,k)=sclim(i,j,k)-sbias
!            tmean(i,j,k)=tmean(i,j,k)-tbias
!            smean(i,j,k)=smean(i,j,k)-sbias
          end do
        end do
      end do

! velocity
      do j=2,jmm1
        uabw(j)=0.
        uabe(j)=0.
        do k=1,kbm1
          uabw(j)=uabw(j)+ub(2,j,k)*dz(k)
          uabe(j)=uabe(j)+ub(im,j,k)*dz(k)
        end do
      end do

      do k=1,kbm1
        do j=1,jm
          ube(j,k)=ub(im,j,k)
          ubw(j,k)=ub(2,j,k)
        end do
      end do

!      do j=1,jm
!        write(6,*) 'j uabw uabe ',j,uabw(j),uabe(j)
!      end do
!      do j=1,jm
!        write(6,*) 'j ubw ube ',j,ubw(j,1),ube(j,1)
!      end do

! lateral viscosity: add a sponge layer at the east (downstream) end
      if(n_east.eq.-1) then
        do j=1,jm
          do i=1,im
           aam2d(i,j)=100.e0*(1.+10.e0*exp(-(im-i)/10.e0))
          end do
        end do
      end if

! generated horizontally averaged density field (in this application,
! the initial condition for density is a function of z (the vertical
! cartesian coordinate) -- when this is not so, make sure that rmean
! has been area averaged before transfer to sigma coordinates)
!      call dens(smean,tmean,rmean)
!      do k=1,kbm1
!        do j=1,jm
!          do i=1,im
!            rmean(i,j,k)=rho(i,j,k)
!          end do
!        end do
!      end do

! set lateral boundary conditions, for use in subroutine bcond (in the
! seamount problem, the east and west boundaries are open, while the
! south and north boundaries are closed through the specification of the
! masks fsm, dum and dvm)
      rfe=1.0
      rfw=1.0
      rfn=1.0
      rfs=1.0

      do j=2,jmm1
!        uabw(j)=uab(2,j)
!        uabe(j)=uab(imm1,j)
! set geostrophically conditioned elevations at the boundaries
!        ele(j)=ele(j-1)-cor(imm1,j)*uab(imm1,j)/grav*dy(imm1,j-1)
!        elw(j)=elw(j-1)-cor(2,j)*uab(2,j)/grav*dy(2,j-1)
        ele(j)=ele(j-1)-cor(imm1,j)*uabe(j)/grav*dy(imm1,j-1)
        elw(j)=elw(j-1)-cor(2,j)*uabw(j)/grav*dy(2,j-1)
      end do

! adjust boundary elevations so that they are zero in the middle of the
! channel
      elejmid=ele(jmm1/2)
      elwjmid=elw(jmm1/2)
      do j=2,jmm1
        ele(j)=(ele(j)-elejmid)*fsm(im,j)
        elw(j)=(elw(j)-elwjmid)*fsm(2,j)
      end do

! set thermodynamic boundary conditions (for the seamount problem, and
! other possible applications, lateral thermodynamic boundary conditions
! are set equal to the initial conditions and are held constant
! thereafter - users may create variable boundary conditions)
      do k=1,kbm1
        do j=1,jm
          tbe(j,k)=tb(im,j,k)
          tbw(j,k)=tb(1,j,k)
          sbe(j,k)=sb(im,j,k)
          sbw(j,k)=sb(1,j,k)
        end do
        do i=1,im
          tbn(i,k)=tb(i,jm,k)
          tbs(i,k)=tb(i,1,k)
          sbn(i,k)=sb(i,jm,k)
          sbs(i,k)=sb(i,1,k)
        end do
      end do

      return
      end

!_______________________________________________________________________
      subroutine update_initial
! update the initial conditions and set the remaining initial conditions
      implicit none
      include 'pom.h'
      integer i,j,k

      do i=1,im
        do j=1,jm
          ua(i,j)=uab(i,j)
          va(i,j)=vab(i,j)
          el(i,j)=elb(i,j)
          et(i,j)=etb(i,j)
          etf(i,j)=et(i,j)
          d(i,j)=h(i,j)+el(i,j)
          dt(i,j)=h(i,j)+et(i,j)
          w(i,j,1)=vfluxf(i,j)
        end do
      end do

      do k=1,kb
        do j=1,jm
          do i=1,im
            l(i,j,k)=0.1*dt(i,j)
            q2b(i,j,k)=small
            q2lb(i,j,k)=l(i,j,k)*q2b(i,j,k)
            kh(i,j,k)=l(i,j,k)*sqrt(q2b(i,j,k))
            km(i,j,k)=kh(i,j,k)
            kq(i,j,k)=kh(i,j,k)
            aam(i,j,k)=aam_init
          end do
        end do
      end do

      do k=1,kbm1
         do j=1,jm      
            do i=1,im
               q2(i,j,k)=q2b(i,j,k)
               q2l(i,j,k)=q2lb(i,j,k)
               t(i,j,k)=tb(i,j,k)
               s(i,j,k)=sb(i,j,k)
               u(i,j,k)=ub(i,j,k)
               v(i,j,k)=vb(i,j,k)
            end do
         end do
      end do
      
!     S.Ohishi 2018.08.20
      if (npg.eq.1) then
         call baropg
      else if (npg.eq.2) then
         call baropg_mcc
      else if (npg .eq. 3) then
         call baropg_thiem
      else
         error_status=1
         write(6,'(/''Error: invalid value for npg'')')
      end if
      
!      call baropg

      do k=1,kbm1
        do j=1,jm
          do i=1,im
            drx2d(i,j)=drx2d(i,j)+drhox(i,j,k)*dz(k)
            dry2d(i,j)=dry2d(i,j)+drhoy(i,j,k)*dz(k)
          end do
        end do
      end do

      return
      end

!_______________________________________________________________________
      subroutine restart_tsuv
! restart using tsuv data 
      implicit none
      include 'pom.h'
      integer i,j,k

!      do i=1,im
!        do j=1,jm
!          uab(i,j)=ua(i,j)
!          vab(i,j)=va(i,j)
!          elb(i,j)=el(i,j)
!          etb(i,j)=el(i,j)
!          et(i,j)=el(i,j)
!          etf(i,j)=el(i,j)
!          d(i,j)=h(i,j)+el(i,j)
!          dt(i,j)=h(i,j)+el(i,j)
!        end do
!      end do

      do k=1,kbm1
         do j=1,jm
            do i=1,im
               tb(i,j,k)=t(i,j,k)
               sb(i,j,k)=s(i,j,k)
               ub(i,j,k)=u(i,j,k)
               vb(i,j,k)=v(i,j,k)
            end do
         end do
      end do

!      call baropg

!      do k=1,kbm1
!        do j=1,jm
!          do i=1,im
!            drx2d(i,j)=drx2d(i,j)+drhox(i,j,k)*dz(k)
!            dry2d(i,j)=dry2d(i,j)+drhoy(i,j,k)*dz(k)
!          end do
!        end do
!      end do

      return
      end

!_______________________________________________________________________
      subroutine bottom_friction
! calculate the bottom friction coefficient
      implicit none
      include 'pom.h'
      integer i,j
! calculate bottom friction
      do j=1,jm
         do i=1,im
            cbc(i,j)=(kappa/log((1.+zz(kbm1))*h(i,j)/z0b))**2
            cbc(i,j)=max(cbcmin,cbc(i,j))
!     if the following is invoked, then it is probable that the wrong
!     choice of z0b or vertical spacing has been made:
            cbc(i,j)=min(cbcmax,cbc(i,j))
         end do
      end do
      
      return
      end
      
!_______________________________________________________________________
      subroutine ztosig(zs,tb,zz,h,t,im,jm,ks,kb,
     $     n_west,n_east,n_south,n_north)
! interpolate vertically
      implicit none
      integer im,jm,ks,kb
      real zs(ks),tb(im,jm,ks),zz(kb),h(im,jm),t(im,jm,kb),tin(ks),
     $                                                  tout(kb),zzh(kb)
      integer n_west,n_east,n_south,n_north
      real tmax
      integer i,j,k

      do j=2,jm-1
         do i=2,im-1
            if (h(i,j).gt.1.0) then
!     special interp on z-lev for cases of no data because h smoothing
               do k=1,ks
                  tin(k)=tb(i,j,k)
                  if (zs(k).le.h(i,j) .and. tin(k).lt.0.01) then
                     tmax=amax1(tb(i-1,j,k),tb(i+1,j,k),
     $                    tb(i,j-1,k),tb(i,j+1,k))
                     tin(k)=tmax
                  endif
                  if (tin(k).lt.0.01 .and. k.ne.1) tin(k)=tin(k-1)
               end do
               
               do k=1,kb
                  zzh(k)=-zz(k)*h(i,j)
               end do
               
!     vertical spline interp
               call splinc(zs,tin,ks,2.e30,2.e30,zzh,tout,kb)
               
               do k=1,kb
                  t(i,j,k)=tout(k)
               end do
               
            end if
         end do
      end do
      call exchange3d_mpi(t,im,jm,kb)
      
!     boundaries
      do k=1,kb
         do j=1,jm
            if(n_west.eq.-1) t(1,j,k)=t(2,j,k)
            if(n_east.eq.-1) t(im,j,k)=t(im-1,j,k)
         end do
         do i=1,im
            if(n_south.eq.-1) t(i,1,k)=t(i,2,k)
            if(n_north.eq.-1) t(i,jm,k)=t(i,jm-1,k)
         end do
      end do
      
      return
      end
      
!_______________________________________________________________________
      subroutine splinc(x,y,n,yp1,ypn,xnew,ynew,m)
! interpolate using splines
      parameter (nmax=300)
      dimension x(n),y(n),y2(nmax),u(nmax),xnew(m),ynew(m)

      if (yp1.gt..99e30) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif

      do i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     $      /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      end do

      if (ypn.gt..99e30) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif

      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
      end do

      do i=1,m
        call splint(x,y,y2,n,xnew(i),ynew(i))
      end do

      return
      end

!_______________________________________________________________________
      subroutine splint(xa,ya,y2a,n,x,y)
      dimension xa(n),ya(n),y2a(n)

      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
        goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.)  then
        error_staus=1
        write(6,'(/a)') 'Error: bad xa input in splint'
      end if
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+
     $      ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
      return
      end

