MODULE common_setting

  IMPLICIT NONE
  PUBLIC

  !---common.f90
  !-----------------------------------------------------------------------
  ! Variable size definitions
  !-----------------------------------------------------------------------
  INTEGER,PARAMETER :: r_size=kind(0.0d0)
  INTEGER,PARAMETER :: r_dble=kind(0.0d0)
  INTEGER,PARAMETER :: r_sngl=kind(0.0e0)

  !-----------------------------------------------------------------------
  ! Constants
  !-----------------------------------------------------------------------
  REAL(r_size),PARAMETER :: pi=3.1415926535d0
  REAL(r_size),PARAMETER :: gg=9.81d0
  REAL(r_size),PARAMETER :: rd=287.0d0
  REAL(r_size),PARAMETER :: cp=1005.7d0
  REAL(r_size),PARAMETER :: re=6371.3d3
  REAL(r_size),PARAMETER :: r_omega=7.292d-5
  REAL(r_size),PARAMETER :: t0c=273.15d0
  REAL(r_size),PARAMETER :: undef=-9.99d33

  !---common_letkf.f90
  !=======================================================================
  !  LEKF Model Independent Parameters
  !=======================================================================
  INTEGER,PARAMETER :: nbv=100    !Ensemble size

  !---common_pom.f90
  !-----------------------------------------------------------------------
  ! General parameters
  !-----------------------------------------------------------------------
  INTEGER,PARAMETER :: nlon=722 !longitude
  INTEGER,PARAMETER :: nlat=386 !latitude
  INTEGER,PARAMETER :: nlev=50  !depth
  INTEGER,PARAMETER :: nv3d=4 ! u,v,t,s
  INTEGER,PARAMETER :: nv2d=3 ! zeta,ubar,vbar
  INTEGER,PARAMETER :: iv3d_u=1
  INTEGER,PARAMETER :: iv3d_v=2
  INTEGER,PARAMETER :: iv3d_t=3
  INTEGER,PARAMETER :: iv3d_s=4
  INTEGER,PARAMETER :: iv2d_z=1
  INTEGER,PARAMETER :: iv2d_ubar=2
  INTEGER,PARAMETER :: iv2d_vbar=3
  INTEGER,PARAMETER :: nij0=nlon*nlat
  INTEGER,PARAMETER :: nlevall=nlev*nv3d+nv2d
  INTEGER,PARAMETER :: ngpv=nij0*nlevall
  REAL(r_size),PARAMETER :: hnosea=1.

  !---common_obs_pom.f90
  INTEGER,PARAMETER :: id_u_obs=2819
  INTEGER,PARAMETER :: id_v_obs=2820
  INTEGER,PARAMETER :: id_t_obs=3073
  INTEGER,PARAMETER :: id_z_obs=2567
  INTEGER,PARAMETER :: id_s_obs=3332

  !---letkf_obs.f90
  INTEGER,PARAMETER :: nslots=1 ! number of time slots for 4D-LETKF
  INTEGER,PARAMETER :: nbslot=1 ! basetime slot

  !Horizontal localization
  LOGICAL,PARAMETER :: DL=.false. !Dual localization
  !Single Localization
  REAL(r_size),PARAMETER :: sigma_obs=300.0d3 ! horizontal [m]
  !Dual Localization
  REAL(r_size),PARAMETER :: sigma_obs_large=1000.d3 !large-scale [m]
  REAL(r_size),PARAMETER :: sigma_obs_small=300.d3 !small-scale [m]
  REAL(r_size),PARAMETER :: g_sigma=300.d3 !E-folding scale in Gauss filter

  !Vertical
  REAL(r_size),PARAMETER :: sigma_obsv=100.0d0 ! vertical [m] 
  !*minus for no v localization

  !Gross error Check using innovation (y-Hx) S.Ohishi 2018.10
  LOGICAL,PARAMETER :: LGE_IO=.false.
  REAL(r_size),PARAMETER :: sstqc_min=-5.d0,sstqc_max=5.d0 !SST
  REAL(r_size),PARAMETER :: tqc_min=-5.d0,tqc_max=5.d0 !Temperature
  REAL(r_size),PARAMETER :: sssqc_min=-1.d0,sssqc_max=1.d0 !SSS
  REAL(r_size),PARAMETER :: sqc_min=-2.d0,sqc_max=2.d0 !Salinity
  REAL(r_size),PARAMETER :: hqc_min=-1.d0,hqc_max=1.d0 !SSH

  !Obs. Error
  LOGICAL,PARAMETER :: OBS_ERR=.false. !Use the following observation error:
  REAL(r_size),PARAMETER :: sst_err=1.5d0   !SST[degree C]
  REAL(r_size),PARAMETER :: t_err=1.5d0     !T[degree C]
  REAL(r_size),PARAMETER :: sss_err=0.3d0  !SSS[psu]
  REAL(r_size),PARAMETER :: s_err=0.3d0    !S[psu] 
  REAL(r_size),PARAMETER :: ssh_err=0.2d0  !SSH[m]
  REAL(r_size),PARAMETER :: u_err=0.2d0    !U[m/s]
  REAL(r_size),PARAMETER :: v_err=0.2d0    !V[m/s]
  LOGICAL,PARAMETER :: OBS_SN=.false. !Use S/N ration for Pf and R
  REAL(r_size),PARAMETER :: sst_sn=1.2d0 !SST
  REAL(r_size),PARAMETER :: t_sn=1.4d0   !T
  REAL(r_size),PARAMETER :: sss_sn=1.2d0 !SSS
  REAL(r_size),PARAMETER :: s_sn=1.4d0   !S 
  REAL(r_size),PARAMETER :: ssh_sn=1.2d0 !SSH
  REAL(r_size),PARAMETER :: u_sn=1.4d0   !U
  REAL(r_size),PARAMETER :: v_sn=1.4d0   !V
  
  REAL(r_size),PARAMETER :: sigma_obst=3.0d0 ! slots

  !Adaptive Background & Observation Error Inflation (AOEI) S.Ohishi
  LOGICAL,PARAMETER :: AOEI=.true.
  REAL(r_size),PARAMETER :: ALPHA_AOEI=1.d0
  !0.d0 <= ALPHA_AOEI <= 1.d0:
  ! (1-alpha) * sigma_b**2 + alpha * sigma_o**2 = B + R = d^o_b**2 

  !Round off for analysis S. Ohishi 2020.03
  LOGICAL,PARAMETER :: ROFF=.true.
  !Lower & Upper limit for Analysis T & S
  REAL(r_size),PARAMETER :: tmin=-1.8d0,tmax=40.d0 
  REAL(r_size),PARAMETER :: smin=30.d0,smax=40.d0

  !---letkf_tools.f90
  ! Multiplicative inflation: cov_infl_mul
  ! > 0: globally constant covariance inflation
  ! < 0: 3D inflation values input from a GPV file "infl_mul.grd"
  REAL(r_size),PARAMETER :: cov_infl_mul = 1.00d0

  !Additive inflation: sp_infl_add
  REAL(r_size),PARAMETER :: sp_infl_add = 0.d0 !additive inflation

  !RTPP/RTPS
  REAL(r_size),PARAMETER :: ALPHA_RTPP   = 0.9d0
  REAL(r_size),PARAMETER :: ALPHA_RTPS   = 0.0d0 

END MODULE common_setting
