! pom.h

! contain parameters for the model domain and for the decomposition
! local domain size, and common POM variables

!_______________________________________________________________________
! Grid size parameters
      integer
     $  im_global      ,! number of global grid points in x
     $  jm_global      ,! number of global grid points in y
     $  kb             ,! number of grid points in z
     $  im_local       ,! number of local grid points in x
     $  jm_local       ,! number of local grid points in y
     $  n_proc          ! number of processors

! im, jm in ssize.h equal im_local, jm_local, respectively

! Correct values for im_local and jm_local are found using
!   n_proc=(im_global-2)/(im_local-2)*(jm_global-2)/(jm_local-2)
! Values higher than necessary will not cause the code to fail, but
! will allocate more memory than is necessary. Value that are too low
! will cause the code to exit

      parameter(
     $  im_global=722 ,
     $  jm_global=386 ,
     $  kb=50         ,
     $  im_local=362  ,
     $  jm_local=194  ,
     $  n_proc=4     )

      integer idevinout     ! device number for output/input files
      parameter(
     $  idevinout=73   )

!_______________________________________________________________________
! Efective grid size
      integer
     $  im             ,! number of grid points used in local x domains
     $  imm1           ,! im-1
     $  imm2           ,! im-2
     $  jm             ,! number of grid points used in local y domains
     $  jmm1           ,! jm-1
     $  jmm2           ,! jm-2
     $  kbm1           ,! kb-1
     $  kbm2            ! kb-2

! Note that im and jm may be different between local domains
! im and jm are equal to or lower than im_local and jm_local, depending
! on the use of correct values for im_local and jm_local

      common/blksiz/
     $  im             ,
     $  imm1           ,
     $  imm2           ,
     $  jm             ,
     $  jmm1           ,
     $  jmm2           ,
     $  kbm1           ,
     $  kbm2

!_______________________________________________________________________
! Parallel variables
      integer
     $  my_task        ,! actual parallel processor ID
     $  master_task    ,! master processor ID
     $  pom_comm       ,! POM model MPI group communicator
     $  i_global       ,! global i index for each point in local domain
     $  j_global       ,! global j index for each point in local domain
     $  n_west         ,! western parallel processor ID
     $  n_east         ,! eastern parallel processor ID
     $  n_south        ,! southern parallel processor ID
     $  n_north        ,! northern parallel processor ID
     $  nproc_x        ,! number of processors in x
     $  nproc_y         ! number of processors in y

      common/blkpar/
     $  my_task        ,
     $  master_task    ,
     $  pom_comm       ,
     $  i_global(im_local),
     $  j_global(jm_local),
     $  n_west         ,
     $  n_east         ,
     $  n_south        ,
     $  n_north        ,
     $  nproc_x        ,
     $  nproc_y        

!_______________________________________________________________________
! Scalars
      integer
     $  iint           ,
     $  iprint         ,! interval in iint at which variables are printed
     $  mode           ,! calculation mode
     $  assim          ,! assimilation (0: No assimilation, 1:assimilation)
     $  ntp            ,! water type
     $  iend           ,! total internal mode time steps
     $  iext           ,
     $  ispadv         ,! step interval for updating external advective terms
     $  isplit         ,! dti/dte
     $  iswtch         ,! time step interval to switch from prtd1 to prtd2
     $  nadv           ,! advection scheme
     $  nbct           ,! surface temperature boundary condition
     $  nbcs           ,! surface salinity boundary condition
     $  budget         ,! conservation each term in budget equations 
     $  ldave          ,! daily average output
     $  llbc           ,! Lateral boudary condition
     $  lssf           ,! fresh water flux
     $  lthf_ws        ,! turbulent heat flux & wind stress algorithm
     $  nitera         ,! number of iterations for Smolarkiewicz scheme
     $  npg            ,! pressure gradient scheme
     $  nread_rst      ,! index to start from restoart file
     $  irestart       ,! restart flag
     $  iens           ,! ensemble flag
     $  error_status   ,! error flag
     $  julday_start   ,! julian day of start time
     $  lbctime_len    ,! total number of lateral boundary data
     $  lbctime_julday ,! start time of lateral boundary data
     $  atmtime_len    ,! total number of atmospheric condition data
     $  atmtime_julday ,! start time of atmospheric condition data
     $  ffluxtime_len    ,! total number of fresh water flux condition data
     $  ffluxtime_julday ,! start time of fresh water flux condition data
     $  tsdatatime_len ,! total number of t/s reference data
     $  tsdatatime_julday ! start time of t/s reference data

      real
     $  alpha          ,! weight for surface slope term in external eq
     $  dte            ,! external (2-D) time step (s)
     $  dti            ,! internal (3-D) time step (s)
     $  dti2           ,! 2*dti
     $  grav           ,! gravity constant (S.I. units)
     $  kappa          ,! von Karman's constant'
     $  pi             ,! pi
     $  ramp           ,! inertial ramp
     $  rfe            ,! flag for eastern open boundary (see bcond)
     $  rfn            ,! flag for northern open boundary (see bcond)
     $  rfs            ,! flag for southern open boundary (see bcond)
     $  rfw            ,! flag for western open boundary (see bcond)
     $  rhoref         ,! reference density
     $  sbias          ,! salinity bias
     $  slmax          ,
     $  small          ,! small value
     $  tbias          ,! temperature bias
     $  time           ,! model time (days)
     $  tprni          ,! inverse horizontal turbulent Prandtl number
     $  umol           ,! background viscosity
     $  vmaxl          ,! max vaf used to test for model blow-up
     $  write_rst      ,
     $  write_ens      ,
     $  aam_init       ,! initial value of aam
     $  cbcmax         ,! maximum bottom friction coefficient
     $  cbcmin         ,! minimum bottom friction coefficient
     $  days           ,! run duration in days
     $  dte2           ,! 2*dte
     $  horcon         ,! smagorinsky diffusivity coefficient
     $  aamadd         ,! Additional viscosity constant
     $  ispi           ,! dte/dti
     $  isp2i          ,! dte/(2*dti)
     $  period         ,! inertial period
     $  prtd1          ,! initial print interval (days)
     $  prtd2          ,! final print interval (days)
     $  smoth          ,! constant to prevent solution splitting
     $  sw             ,! smoothing parameter for Smolarkiewicz scheme
     $  swtch          ,! time to switch from prtd1 to prtd2
     $  time0          ,! initial time (days)
     $  z0b            ,! bottom roughness
     $  atmtime_dayint ,! time interval of atmospheric data
     $  ffluxtime_dayint   ,! time interval of fresh water flux data
     $  lbctime_dayint ,! time interval of lateral boundary data
     $  tsdatatime_dayint  ,! time interval of t/s reference data
     $  ts_nudge       ,! SST nuding timescale [day]
     $  ti_nudge       ,! T nuding timescale [day]
     $  ss_nudge       ,! SSS nuding timescale [day]
     $  si_nudge        ! S nuding timescale [day]

      common/blkcon/
     $  alpha          ,
     $  dte            ,
     $  dti            ,
     $  dti2           ,
     $  grav           ,
     $  kappa          ,
     $  pi             ,
     $  ramp           ,
     $  rfe            ,
     $  rfn            ,
     $  rfs            ,
     $  rfw            ,
     $  rhoref         ,
     $  sbias          ,
     $  slmax          ,
     $  small          ,
     $  tbias          ,
     $  time           ,
     $  tprni          ,
     $  umol           ,
     $  vmaxl          ,
     $  write_rst      ,
     $  write_ens      ,
     $  iint           ,
     $  iprint         ,
     $  mode           ,
     $  assim          ,
     $  ntp            ,
     $  aam_init       ,
     $  cbcmax         ,
     $  cbcmin         ,
     $  days           ,
     $  dte2           ,
     $  horcon         ,
     $  aamadd         , 
     $  ispi           ,
     $  isp2i          ,
     $  period         ,
     $  prtd1          ,
     $  prtd2          ,
     $  smoth          ,
     $  sw             ,
     $  swtch          ,
     $  time0          ,
     $  z0b            ,
     $  iend           ,
     $  iext           ,
     $  ispadv         ,
     $  isplit         ,
     $  iswtch         ,
     $  nadv           ,
     $  nbct           ,
     $  nbcs           ,
     $  budget         ,
     $  ldave          ,
     $  llbc           ,
     $  lssf           ,
     $  lthf_ws        ,
     $  nitera         ,
     $  npg            ,
     $  nread_rst      ,
     $  irestart       ,
     $  iens           ,
     $  error_status   ,
     $  julday_start   ,
     $  lbctime_len    ,
     $  lbctime_julday ,
     $  lbctime_dayint ,
     $  atmtime_len    ,
     $  atmtime_julday ,
     $  atmtime_dayint ,
     $  ffluxtime_len    ,
     $  ffluxtime_julday ,
     $  ffluxtime_dayint ,
     $  tsdatatime_len    ,
     $  tsdatatime_julday ,
     $  tsdatatime_dayint ,
     $  ts_nudge       ,
     $  ti_nudge       ,
     $  ss_nudge       ,
     $  si_nudge

!_______________________________________________________________________
! 1-D arrays
      real
     $  dz             ,! z(k)-z(k+1)
     $  dzz            ,! zz(k)-zz(k+1)
     $  z              ,! sigma coordinate from z=0 (surface) to z=-1 (bottom)
     $  zz              ! sigma coordinate, intermediate between z

      common/blk1d/ 
     $  dz(kb)         ,
     $  dzz(kb)        ,
     $  z(kb)          ,
     $  zz(kb)

!_______________________________________________________________________
! 2-D arrays
      real
     $  aam2d          ,! vertical average of aam
     $  advua          ,! sum of the 2nd, 3rd and 4th terms in eq (18)
     $  advva          ,! sum of the 2nd, 3rd and 4th terms in eq (19)
     $  adx2d          ,! vertical integral of advx
     $  ady2d          ,! vertical integral of advy
     $  art            ,! cell area centered on T grid points
     $  aru            ,! cell area centered on U grid points
     $  arv            ,! cell area centered on V grid points
     $  cbc            ,! bottom friction coefficient
     $  cor            ,! coriolis parameter
     $  d              ,! h+el
     $  drx2d          ,! vertical integral of drhox
     $  dry2d          ,! vertical integral of drhoy
     $  dt             ,! h+et
     $  dum            ,! mask for u velocity
     $  dvm            ,! mask for v velocity
     $  dx             ,! grid spacing in x
     $  dy             ,! grid spacing in y
     $  east_c         ,! horizontal coordinate of cell corner points in x
     $  east_e         ,! horizontal coordinate of elevation points in x
     $  east_u         ,! horizontal coordinate of U points in x
     $  east_v         ,! horizontal coordinate of V points in x
     $  e_atmos        ,! atmospheric pressure
     $  egb            ,! surface elevation use for pressure gradient at time n-1
     $  egf            ,! surface elevation use for pressure gradient at time n+1
     $  el             ,! surface elevation used in the external mode at time n
     $  elb            ,! surface elevation used in the external mode at time n-1
     $  elf            ,! surface elevation used in the external mode at time n+1
     $  el_dave        ,! surface elevation used in the external mode averaged 1 day
     $  el_iau         ,! surface elevation difference for iau analysis (anal-model)
     $  et             ,! surface elevation used in the internal mode at time n
     $  etb            ,! surface elevation used in the internal mode at time n-1
     $  etf            ,! surface elevation used in the internal mode at time n+1
     $  fluxua         ,
     $  fluxva         ,
     $  fsm            ,! mask for scalar variables
     $  h              ,! bottom depth
     $  north_c        ,! horizontal coordinate of cell corner points in y
     $  north_e        ,! horizontal coordinate of elevation points in y
     $  north_u        ,! horizontal coordinate of U points in y
     $  north_v        ,! horizontal coordinate of V points in y
     $  psi            ,
     $  rot            ,! rotation angle
     $  ssurf          ,
     $  swrad          ,! short wave radiation incident on the ocean surface
     $  lhf            ,!latent heat flux [W/m^2]
     $  lhf_dave       ,!latent heat flux [W/m^2] averaged within 1day
     $  shf            ,!sensible heat flux [W/m^2]
     $  shf_dave       ,!sensible heat flux [W/m^2] averaged within 1day
     $  swr            ,!shortwave radiation [W/m^2]
     $  swr_dave       ,!shortwave radiation [W/m^2] averaged within 1day
     $  lwr            ,!longwave radiation [W/m^2]
     $  lwr_dave       ,!longwave radiation [W/m^2] averaged within 1day
     $  vfluxb         ,! volume flux through water column surface at time n-1
     $  vfluxf         ,! volume flux through water column surface at time n+1
     $  tps            ,
     $  tsurf          ,
     $  ua             ,! vertical mean of u at time n
     $  uab            ,! vertical mean of u at time n-1
     $  uaf            ,! vertical mean of u at time n+1
     $  ua_iau         ,! difference of vertical mean of u for iau analysis
     $  utb            ,! ua time averaged over the interval dti at time n-1
     $  utf            ,! ua time averaged over the interval dti at time n+1
     $  va             ,! vertical mean of v at time n
     $  vab            ,! vertical mean of v at time n-1
     $  vaf            ,! vertical mean of v at time n+1
     $  va_iau         ,! difference of vertical mean of v for iau analysis
     $  vtb            ,! va time averaged over the interval dti at time n-1
     $  vtf            ,! va time averaged over the interval dti at time n+1
     $  wssurf         ,! <ws(0)> salinity flux at the surface
     $  wtsurf         ,! <wt(0)> temperature flux at the surface
     $  wubot          ,! x-momentum flux at the bottom
     $  wusurf         ,! <wu(0)> momentum flux at the surface
     $  wvbot          ,! y-momentum flux at the bottom
     $  wvsurf         ,! <wv(0)> momentum flux at the surface
     $  sfcterm        ,! surface forcing term
     $  tsfcf          ,! temperature surface forcing term
     $  tsfc           ,!
     $  tsfcb          ,!
     $  tsfc_dave      ,! temperature surface forcing averaged within 1day
     $  ssfcf          ,! salinity surface forcing term
     $  ssfc           ,!
     $  ssfcb          ,!
     $  ssfc_dave       ! salinity surface forcing averaged within 1day

      common/blk2d/
     $  aam2d(im_local,jm_local)   ,
     $  advua(im_local,jm_local)   ,
     $  advva(im_local,jm_local)   ,
     $  adx2d(im_local,jm_local)   ,
     $  ady2d(im_local,jm_local)   ,
     $  art(im_local,jm_local)     ,
     $  aru(im_local,jm_local)     ,
     $  arv(im_local,jm_local)     ,
     $  cbc(im_local,jm_local)     ,
     $  cor(im_local,jm_local)     ,
     $  d(im_local,jm_local)       ,
     $  drx2d(im_local,jm_local)   ,
     $  dry2d(im_local,jm_local)   ,
     $  dt(im_local,jm_local)      ,
     $  dum(im_local,jm_local)     ,
     $  dvm(im_local,jm_local)     ,
     $  dx(im_local,jm_local)      ,
     $  dy(im_local,jm_local)      ,
     $  east_c(im_local,jm_local)  ,
     $  east_e(im_local,jm_local)  ,
     $  east_u(im_local,jm_local)  ,
     $  east_v(im_local,jm_local)  ,
     $  e_atmos(im_local,jm_local) ,
     $  egb(im_local,jm_local)     ,
     $  egf(im_local,jm_local)     ,
     $  el(im_local,jm_local)      ,
     $  elb(im_local,jm_local)     ,
     $  elf(im_local,jm_local)     ,
     $  el_dave(im_local,jm_local) ,
     $  el_iau(im_local,jm_local)  ,
     $  et(im_local,jm_local)      ,
     $  etb(im_local,jm_local)     ,
     $  etf(im_local,jm_local)     ,
     $  fluxua(im_local,jm_local)  ,
     $  fluxva(im_local,jm_local)  ,
     $  fsm(im_local,jm_local)     ,
     $  h(im_local,jm_local)       ,
     $  north_c(im_local,jm_local) ,
     $  north_e(im_local,jm_local) ,
     $  north_u(im_local,jm_local) ,
     $  north_v(im_local,jm_local) ,
     $  psi(im_local,jm_local)     ,
     $  rot(im_local,jm_local)     ,
     $  ssurf(im_local,jm_local)   ,
     $  swrad(im_local,jm_local)   ,
     $  lhf(im_local,jm_local)     ,
     $  lhf_dave(im_local,jm_local),
     $  shf(im_local,jm_local)     ,
     $  shf_dave(im_local,jm_local),
     $  lwr(im_local,jm_local)     ,
     $  lwr_dave(im_local,jm_local),
     $  swr(im_local,jm_local)     ,
     $  swr_dave(im_local,jm_local),
     $  vfluxb(im_local,jm_local)  ,
     $  vfluxf(im_local,jm_local)  ,
     $  tps(im_local,jm_local)     ,
     $  tsurf(im_local,jm_local)   ,
     $  ua(im_local,jm_local)      ,
     $  uab(im_local,jm_local)     ,
     $  uaf(im_local,jm_local)     ,
     $  ua_iau(im_local,jm_local)  ,
     $  utb(im_local,jm_local)     ,
     $  utf(im_local,jm_local)     ,
     $  va(im_local,jm_local)      ,
     $  vab(im_local,jm_local)     ,
     $  vaf(im_local,jm_local)     ,
     $  va_iau(im_local,jm_local)  ,
     $  vtb(im_local,jm_local)     ,
     $  vtf(im_local,jm_local)     ,
     $  wssurf(im_local,jm_local)  ,
     $  wtsurf(im_local,jm_local)  ,
     $  wubot(im_local,jm_local)   ,
     $  wusurf(im_local,jm_local)  ,
     $  wvbot(im_local,jm_local)   ,
     $  wvsurf(im_local,jm_local)  ,
     $  sfcterm(im_local,jm_local) ,! surface forcing term
     $  tsfcf(im_local,jm_local)   ,! temperature surface forcing term
     $  tsfc(im_local,jm_local)    ,!
     $  tsfcb(im_local,jm_local)   ,!
     $  tsfc_dave(im_local,jm_local) ,! temperature surface forcing averaged within 1day
     $  ssfcf(im_local,jm_local)   ,! salinity surface forcing term
     $  ssfc(im_local,jm_local)    ,!
     $  ssfcb(im_local,jm_local)   ,!
     $  ssfc_dave(im_local,jm_local) ! salinity surface forcing averaged within 1day

!_______________________________________________________________________
! 3-D arrays
      real
     $  aam            ,! horizontal kinematic viscosity
     $  aam_dave       ,! horizontal kinematic viscosity averaged within 1 day
     $  advx           ,! x-horizontal advection and diffusion terms
     $  advy           ,! y-horizontal advection and diffusion terms
     $  drhox          ,! x-component of the internal baroclinic pressure
     $  drhoy          ,! y-component of the internal baroclinic pressure
     $  dtef           ,
     $  kh             ,! vertical diffusivity
     $  kh_dave        ,! vertical diffusivity averaged within 1 day
     $  km             ,! vertical kinematic viscosity
     $  kq             ,
     $  l              ,! turbulence length scale
     $  q2b            ,! twice the turbulent kinetic energy at time n-1
     $  q2             ,! twice the turbulent kinetic energy at time n
     $  q2lb           ,! q2 x l at time n-1
     $  q2l            ,! q2 x l at time n
     $  rho            ,! density
     $  sb             ,! salinity at time n-1
     $  sclim          ,! horizontally averaged salinity
     $  s              ,! salinity at time n
     $  s_dave         ,! salinity averaged within 1 day
     $  s_iau          ,! salinity difference for iau analysis
     $  dsdt           ,! salinity tendency in 1 day
     $  dsdt_dave      ,! salinity tendency in 1 day averaged within 1day
     $  tb             ,! temperature at time n-1
     $  tclim          ,! horizontally averaged temperature
     $  t              ,! temperature at time n
     $  t_dave         ,! temperature averaged within 1 day
     $  t_iau          ,! temperature difference for iau analysis
     $  dtdt           ,! temperature tendency in 1 day
     $  dtdt_dave      ,!                               averaged within 1 day
     $  ub             ,! horizontal velocity in x at time n-1
     $  uf             ,! horizontal velocity in x at time n+1
     $  u              ,! horizontal velocity in x at time n
     $  u_dave         ,! horizontal velocity in x averaged within 1 day
     $  u_iau          ,! horizontal velocity difference in x for iau analysis 
     $  vb             ,! horizontal velocity in y at time n-1
     $  vf             ,! horizontal velocity in y at time n+1
     $  v              ,! horizontal velocity in y at time n
     $  v_dave         ,! horizontal velocity in y averaged within 1 day
     $  v_iau          ,! horizontal velocity difference in y for iau analysis 
     $  w              ,! sigma coordinate vertical velocity
     $  wr             ,! real (z coordinate) vertical velocity
     $  wr_dave        ,! real (z coordinate) vertical velocity averaged within 1 day
! conservation each term
     $  xadvterm       ,! advection term - x direction
     $  yadvterm       ,!                - y direction
     $  zadvterm       ,!                - z direction
     $  advterm        ,! advection term - x+y+z direction
     $  xdifterm       ,! diffusion term - x direction
     $  ydifterm       ,!                - y direction
     $  zdifterm       ,!                - z direction
     $  radterm        ,! radiation term
     $  txadv          ,! temperature advection term -x direction
     $  tyadv          ,!                            -y direction
     $  tzadv          ,!                            -z direction 
     $  tadv           ,!                            -total 
     $  txdif          ,! temperature diffusion term -x direction
     $  tydif          ,! 
     $  tzdif          ,!
     $  qz             ,! shortwave penetration term
     $  troff          ,! temperature round off
     $  tnudge         ,! temperature nuding
     $  txadv_dave     ,! daily average of temperature advection term
     $  tyadv_dave     ,! 
     $  tzadv_dave     ,! 
     $  tadv_dave      ,! 
     $  txdif_dave     ,!                              diffusion term
     $  tydif_dave     ,! 
     $  tzdif_dave     ,!
     $  qz_dave        ,! daily average of shortwave penetration term
     $  troff_dave     ,! daily temperature round off
     $  tnudge_dave    ,! daily temperature nuding
     $  sxadv          ,! salinity advection term -x direction
     $  syadv          ,!                         -y direction
     $  szadv          ,!                         -z direction
     $  sadv           ,!                         -total
     $  sxdif          ,! salinity diffusion term - xdirection
     $  sydif          ,!                         - ydirection
     $  szdif          ,!                         - zdirection
     $  sroff          ,! salinity round off
     $  snudge         ,! salinity nudging
     $  sxadv_dave     ,! daily average of salinity advection term
     $  syadv_dave     ,! 
     $  szadv_dave     ,! 
     $  sadv_dave      ,! 
     $  sxdif_dave     ,!                           diffusion term
     $  sydif_dave     ,! 
     $  szdif_dave     ,!
     $  sroff_dave     ,! daily salinify round off
     $  snudge_dave     !daily salinity nuding

      common/blk3d/
     $  aam(im_local,jm_local,kb)  ,
     $  aam_dave(im_local,jm_local,kb)  ,
     $  advx(im_local,jm_local,kb) ,
     $  advy(im_local,jm_local,kb) ,
     $  drhox(im_local,jm_local,kb),
     $  drhoy(im_local,jm_local,kb),
     $  dtef(im_local,jm_local,kb) ,
     $  kh(im_local,jm_local,kb)   ,
     $  kh_dave(im_local,jm_local,kb)   ,
     $  km(im_local,jm_local,kb)   ,
     $  kq(im_local,jm_local,kb)   ,
     $  l(im_local,jm_local,kb)    ,
     $  q2b(im_local,jm_local,kb)  ,
     $  q2(im_local,jm_local,kb)   ,
     $  q2lb(im_local,jm_local,kb) ,
     $  q2l(im_local,jm_local,kb)  ,
     $  rho(im_local,jm_local,kb)  ,
     $  sb(im_local,jm_local,kb)   ,
     $  sclim(im_local,jm_local,kb),
     $  s(im_local,jm_local,kb)    ,
     $  s_dave(im_local,jm_local,kb),
     $  s_iau(im_local,jm_local,kb),
     $  dsdt(im_local,jm_local,kb) ,
     $  dsdt_dave(im_local,jm_local,kb),
     $  tb(im_local,jm_local,kb)   ,
     $  tclim(im_local,jm_local,kb),
     $  t(im_local,jm_local,kb)    ,
     $  t_dave(im_local,jm_local,kb),
     $  t_iau(im_local,jm_local,kb),
     $  dtdt(im_local,jm_local,kb) ,
     $  dtdt_dave(im_local,jm_local,kb) ,
     $  ub(im_local,jm_local,kb)   ,
     $  uf(im_local,jm_local,kb)   ,
     $  u(im_local,jm_local,kb)    ,
     $  u_dave(im_local,jm_local,kb),
     $  u_iau(im_local,jm_local,kb),
     $  vb(im_local,jm_local,kb)   ,
     $  vf(im_local,jm_local,kb)   ,
     $  v(im_local,jm_local,kb)    ,
     $  v_dave(im_local,jm_local,kb),
     $  v_iau(im_local,jm_local,kb),
     $  w(im_local,jm_local,kb)    ,
     $  wr(im_local,jm_local,kb)   ,
     $  wr_dave(im_local,jm_local,kb),
     $  xadvterm(im_local,jm_local,kb)    ,
     $  yadvterm(im_local,jm_local,kb)    ,
     $  zadvterm(im_local,jm_local,kb)    ,
     $  advterm(im_local,jm_local,kb)     ,     
     $  xdifterm(im_local,jm_local,kb)    ,
     $  ydifterm(im_local,jm_local,kb)    ,
     $  zdifterm(im_local,jm_local,kb)    ,
     $  radterm(im_local,jm_local,kb)     ,
     $  txadv(im_local,jm_local,kb)       ,
     $  tyadv(im_local,jm_local,kb)       ,
     $  tzadv(im_local,jm_local,kb)       ,
     $  tadv(im_local,jm_local,kb)        ,
     $  txdif(im_local,jm_local,kb)       ,
     $  tydif(im_local,jm_local,kb)       ,
     $  tzdif(im_local,jm_local,kb)       ,
     $  qz(im_local,jm_local,kb)          ,
     $  troff(im_local,jm_local,kb)       ,
     $  tnudge(im_local,jm_local,kb)      ,
     $  txadv_dave(im_local,jm_local,kb)  ,
     $  tyadv_dave(im_local,jm_local,kb)  ,
     $  tzadv_dave(im_local,jm_local,kb)  ,
     $  tadv_dave(im_local,jm_local,kb)   ,
     $  txdif_dave(im_local,jm_local,kb)  ,
     $  tydif_dave(im_local,jm_local,kb)  ,
     $  tzdif_dave(im_local,jm_local,kb)  ,
     $  qz_dave(im_local,jm_local,kb)     ,
     $  troff_dave(im_local,jm_local,kb)  ,
     $  tnudge_dave(im_local,jm_local,kb) ,
     $  sxadv(im_local,jm_local,kb)       ,
     $  syadv(im_local,jm_local,kb)       ,
     $  szadv(im_local,jm_local,kb)       ,
     $  sadv(im_local,jm_local,kb)        ,
     $  sxdif(im_local,jm_local,kb)       ,
     $  sydif(im_local,jm_local,kb)       ,
     $  szdif(im_local,jm_local,kb)       ,
     $  sroff(im_local,jm_local,kb)       ,
     $  snudge(im_local,jm_local,kb)        ,
     $  sxadv_dave(im_local,jm_local,kb)    ,
     $  syadv_dave(im_local,jm_local,kb)    ,
     $  szadv_dave(im_local,jm_local,kb)    ,
     $  sadv_dave(im_local,jm_local,kb)     ,
     $  sxdif_dave(im_local,jm_local,kb)    ,
     $  sydif_dave(im_local,jm_local,kb)    ,
     $  szdif_dave(im_local,jm_local,kb)    ,
     $  sroff_dave(im_local,jm_local,kb)    ,
     $  snudge_dave(im_local,jm_local,kb)

!_______________________________________________________________________
! 1 and 2-D boundary value arrays
      real
     $  ele,ele0,ele1    ,! elevation at the eastern open boundary
     $  eln,eln0,eln1    ,! elevation at the northern open boundary
     $  els,els0,els1    ,! elevation at the southern open boundary
     $  elw,elw0,elw1    ,! elevation at the western open boundary
     $  sbe,sbe0,sbe1    ,! salinity at the eastern open boundary
     $  sbn,sbn0,sbn1    ,! salinity at the northern open boundary
     $  sbs,sbs0,sbs1    ,! salinity at the southern open boundary
     $  sbw,sbw0,sbw1    ,! salinity at the western open boundary
     $  tbe,tbe0,tbe1    ,! temperature at the eastern open boundary
     $  tbn,tbn0,tbn1    ,! temperature at the northern open boundary
     $  tbs,tbs0,tbs1    ,! temperature at the southern open boundary
     $  tbw,tbw0,tbw1    ,! temperature at the western open boundary
     $  uabe,uabe0,uabe1 ,! vertical mean of u at the eastern open boundary
     $  uabn,uabn0,uabn1 ,! vertical mean of u at the northern open boundary
     $  uabs,uabs0,uabs1 ,! vertical mean of u at the southern open boundary
     $  uabw,uabw0,uabw1 ,! vertical mean of u at the western open boundary
     $  ube,ube0,ube1    ,! u at the eastern open boundary
     $  ubn,ubn0,ubn1    ,! u at the northern open boundary
     $  ubs,ubs0,ubs1    ,! u at the southern open boundary
     $  ubw,ubw0,ubw1    ,! u at the western open boundary
     $  vabe,vabe0,vabe1 ,! vertical mean of v at the eastern open boundary
     $  vabn,vabn0,vabn1 ,! vertical mean of v at the northern open boundary
     $  vabs,vabs0,vabs1 ,! vertical mean of v at the southern open boundary
     $  vabw,vabw0,vabw1 ,! vertical mean of v at the western open boundary
     $  vbe,vbe0,vbe1    ,! v at the eastern open boundary
     $  vbn,vbn0,vbn1    ,! v at the northern open boundary
     $  vbs,vbs0,vbs1    ,! v at the southern open boundary
     $  vbw,vbw0,vbw1     ! v at the western open boundary

      common/bdry/
     $  ele(jm_local),ele0(jm_local),ele1(jm_local)         ,
     $  eln(im_local),eln0(im_local),eln1(im_local)         ,        
     $  els(im_local),els0(im_local),els1(im_local)         ,
     $  elw(jm_local),elw0(jm_local),elw1(jm_local)         ,
     $  sbe(jm_local,kb),sbe0(jm_local,kb),sbe1(jm_local,kb),
     $  sbn(im_local,kb),sbn0(im_local,kb),sbn1(im_local,kb),
     $  sbs(im_local,kb),sbs0(im_local,kb),sbs1(im_local,kb),
     $  sbw(jm_local,kb),sbw0(jm_local,kb),sbw1(jm_local,kb),
     $  tbe(jm_local,kb),tbe0(jm_local,kb),tbe1(jm_local,kb),
     $  tbn(im_local,kb),tbn0(im_local,kb),tbn1(im_local,kb),
     $  tbs(im_local,kb),tbs0(im_local,kb),tbs1(im_local,kb),
     $  tbw(jm_local,kb),tbw0(jm_local,kb),tbw1(jm_local,kb),
     $  uabe(jm_local),uabe0(jm_local),uabe1(jm_local)      ,
     $  uabn(im_local),uabn0(im_local),uabn1(im_local)      ,
     $  uabs(im_local),uabs0(im_local),uabs1(im_local)      ,
     $  uabw(jm_local),uabw0(jm_local),uabw1(jm_local)      ,
     $  ube(jm_local,kb),ube0(jm_local,kb),ube1(jm_local,kb),
     $  ubn(im_local,kb),ubn0(im_local,kb),ubn1(im_local,kb),
     $  ubs(im_local,kb),ubs0(im_local,kb),ubs1(im_local,kb),
     $  ubw(jm_local,kb),ubw0(jm_local,kb),ubw1(jm_local,kb),
     $  vabe(jm_local),vabe0(jm_local),vabe1(jm_local)      ,
     $  vabn(im_local),vabn0(im_local),vabn1(im_local)      ,
     $  vabs(im_local),vabs0(im_local),vabs1(im_local)      ,
     $  vabw(jm_local),vabw0(jm_local),vabw1(jm_local)      ,
     $  vbe(jm_local,kb),vbe0(jm_local,kb),vbe1(jm_local,kb),
     $  vbn(im_local,kb),vbn0(im_local,kb),vbn1(im_local,kb),
     $  vbs(im_local,kb),vbs0(im_local,kb),vbs1(im_local,kb),
     $  vbw(jm_local,kb),vbw0(jm_local,kb),vbw1(jm_local,kb)

!_______________________________________________________________________
! Character variables
      character*26
     $  time_start      ! date and time of start of initial run of model

      character*40
     $  source         ,
     $  title

      character*120
     $  netcdf_file    ,
     $  read_rst_file  ,
     $  read_iau_file  ,
     $  write_rst_file ,
     $  write_ens_file

      common/blkchar/
     $  time_start     ,
     $  source         ,
     $  title          ,
     $  netcdf_file    ,
     $  read_rst_file  ,
     $  read_iau_file  ,
     $  write_rst_file ,
     $  write_ens_file

!_______________________________________________________________________
! Logical variables
      logical lramp
      logical lfrs
      logical lmsm
      logical lrtvf
      logical lmynnf
      logical lmbws
      logical liwbrk
      logical lroff
      logical lnan

! seto2 experiments
      logical ltide
      logical lwbrk

      common/blklog/ lramp,lfrs,
     $               lrtvf,lmsm,lmynnf,lmbws,liwbrk,
     $               lroff,lnan,
! seto2 experiments
     $               ltide,lwbrk
!_______________________________________________________________________
! Additional 2-D arrays
      real
     $  windu,windu0,windu1,windu_dave  ,! east-west component of wind
     $  windv,windv0,windv1,windv_dave  ,! north-south component of wind
     $  winds,winds0,winds1,winds_dave  ,! wind speed
     $  tauu,tauu_dave                  ,! east-west component of wind stress
     $  tauv,tauv_dave                  ,! north-south component of wind stress
     $  taus,taus_dave                  ,! magnitude of wind stress
     $  airt,airt0,airt1                ,! air temperature [K]
     $  ta,ta_dave                      ,! air temperature [degree C]
     $  airh,airh0,airh1                ,! air humidity [g/g]
     $  qa,qa_dave, ! air humidity [g/kg] 
     $  qs,qs_dave, ! surface saturated specific humidity [g/kg]
     $  swrad0,swrad1        ,! short wave radiation
     $  cloud,cloud0,cloud1  ,! cloud fraction
     $  slp,slp0,slp1        ,! sea level pressure [Pa]
     $  evap,evap0,evap1,evap_dave     ,! evaporation (E) [mm/day]
     $  prep,prep0,prep1,prep_dave     ,! precipitation (P) 
     $  river,river0,river1,river_dave ,! river discharge (R)
     $  fflux,fflux0,fflux1            ,! fresh water flux (E-P-R)
     $  eflux,eflux_dave               ,! evaporation flux [m/s]
     $  pflux,pflux_dave               ,! precipitation flux
     $  rflux,rflux_dave                ! river discharge flux

      common/addblk2d/
     $  windu(im_local,jm_local)   ,
     $  windu0(im_local,jm_local)  ,
     $  windu1(im_local,jm_local)  ,
     $  windu_dave(im_local,jm_local),
     $  windv(im_local,jm_local)   ,
     $  windv0(im_local,jm_local)  ,
     $  windv1(im_local,jm_local)  ,
     $  windv_dave(im_local,jm_local),
     $  winds(im_local,jm_local)   ,
     $  winds0(im_local,jm_local)  ,
     $  winds1(im_local,jm_local)  ,
     $  winds_dave(im_local,jm_local),
     $  tauu(im_local,jm_local)     ,
     $  tauu_dave(im_local,jm_local),
     $  tauv(im_local,jm_local)     ,
     $  tauv_dave(im_local,jm_local),
     $  taus(im_local,jm_local)     ,
     $  taus_dave(im_local,jm_local),
     $  airt(im_local,jm_local)    ,
     $  airt0(im_local,jm_local)   ,
     $  airt1(im_local,jm_local)   ,
     $  ta(im_local,jm_local)      ,
     $  ta_dave(im_local,jm_local) ,
     $  airh(im_local,jm_local)    ,
     $  airh0(im_local,jm_local)   ,
     $  airh1(im_local,jm_local)   ,
     $  qa(im_local,jm_local)      ,
     $  qa_dave(im_local,jm_local) ,
     $  qs(im_local,jm_local)      ,
     $  qs_dave(im_local,jm_local) ,
     $  swrad0(im_local,jm_local)  ,
     $  swrad1(im_local,jm_local)  ,
     $  cloud(im_local,jm_local)   ,
     $  cloud0(im_local,jm_local)  ,
     $  cloud1(im_local,jm_local)  ,
     $  slp(im_local,jm_local)     ,
     $  slp0(im_local,jm_local)    ,
     $  slp1(im_local,jm_local)    ,
     $  evap(im_local,jm_local)    ,
     $  evap0(im_local,jm_local)   ,
     $  evap1(im_local,jm_local)   ,
     $  evap_dave(im_local,jm_local),
     $  prep(im_local,jm_local)     ,
     $  prep0(im_local,jm_local)    ,
     $  prep1(im_local,jm_local)    ,
     $  prep_dave(im_local,jm_local),
     $  river(im_local,jm_local)     ,
     $  river0(im_local,jm_local)    ,
     $  river1(im_local,jm_local)    ,
     $  river_dave(im_local,jm_local),
     $  fflux(im_local,jm_local)     ,
     $  fflux0(im_local,jm_local)    ,
     $  fflux1(im_local,jm_local)    , 
     $  eflux(im_local,jm_local)     ,
     $  eflux_dave(im_local,jm_local),
     $  pflux(im_local,jm_local)     ,
     $  pflux_dave(im_local,jm_local),
     $  rflux(im_local,jm_local)     ,
     $  rflux_dave(im_local,jm_local)
!_______________________________________________________________________
! Additional 3-D arrays
      real
     $  tclimm, tclim0, tclim1      ,! monthly climatological temperature
     $  sclimm, sclim0, sclim1      ,! monthly climatological salinity
     $  tref, tref0, tref1          ,! reference temperature
     $  sref, sref0, sref1           ! reference salinity

      common/addblk3d/
     $  tclimm(im_local,jm_local,kb)  ,
     $  tclim0(im_local,jm_local,kb)  ,
     $  tclim1(im_local,jm_local,kb)  ,
     $  sclimm(im_local,jm_local,kb)  ,
     $  sclim0(im_local,jm_local,kb)  ,
     $  sclim1(im_local,jm_local,kb)  ,
     $  tref(im_local,jm_local,kb)   ,
     $  tref0(im_local,jm_local,kb)  ,
     $  tref1(im_local,jm_local,kb)  ,
     $  sref(im_local,jm_local,kb)   ,
     $  sref0(im_local,jm_local,kb)  ,
     $  sref1(im_local,jm_local,kb)
