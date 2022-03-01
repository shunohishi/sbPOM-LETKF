c
c Tidal module for OCMs
c Sergey Varlamov, 2001-2011
c
c See example and test program at the end of this file.
c
c All values are in cgs units: cm, cm^2/s
c
c Modifications 20110624:
c
c   1. Added support for OTIS ATLAS solutions that incorporate local solutuions
c      into the regional (PO) and global solution (tpxo7)
c   2. Earlier (20100615) was added support for OTIS 2009 solutions.
c
c Modification 20081210:
c   1. Changed "nx_full" to "nx_full-1" in input tidal grid definition (bug)
c   2. Added support for fes2004 data
c
c Modification 20070312:
c   1. Add support for the NAO regional tidal UV data
C
c Modification 20050601:
c   1. Support for 4 long-term tidal harmonics added: 'mm ','mf ','ssa','sa '
c   2. If number of requested harmonics is zero, program checks the presence
c      of data files for all 20 harmonics and use all available;
c      absent harmonics names are reported
c   3. Test program was not modified yet for new version!
c
c
c Directives for activating test mode
c
!#ifdef VFD
 !MS$  DEFINE _TDTEST_=1
!#endif

!#define _TDTEST_  1
!#define CONVERT2BIGENDIAN
C--------------------------------------------------------------------------

      MODULE TIDE_SERV

        IMPLICIT NONE

        TYPE TOPO_G
C
C         REQUEST IS SPECIFIED AT THE INITIALIZATION STEP
C         AND IS MODIFIED TO SELECTED COVERING AREA.
C         VALUES MUST BE ORDERED SO THAT MIN < MAX
C
          REAL    XMIN, XMAX
          REAL    YMIN, YMAX
C
C         NEXT PARAMETERS ARE READ FROM THE TIDE FILE HEADERS.
C         SAVED ARE INTERPOLATION COEFFICIENTS FOR THE REGION
C         SPECIFIED ABOUT THAT COULD BE SHIFTED
C         TO THE CLOSEST EXTERNAL GRID BOUNDARIES.
C         GRID NOW IS REGULAR WITH FIXED STEPS ONLY.
C
          REAL    X0_FULL, X1_FULL
          REAL    Y0_FULL, Y1_FULL
          REAL    DX, DY
          INTEGER NX_FULL, NY_FULL
          INTEGER NX, NY
          INTEGER GLOBX
        END TYPE

        INTEGER,PARAMETER::       TIDE_NUM_MODELS = 9,
     *                            TIDE_NAO99_GLOBAL = 1,
     *                            TIDE_NAO99_JAPAN_SSH = 2,
     *                            TIDE_NAO99_JAPAN_SSH_UVH = 3,
     *                            TIDE_LEGOS_2002 = 4,
     *                            TIDE_FES2004 = 5,
     *                            TIDE_OTIS_GLOBAL = 6,   ! other OTIS models could come after it
     *                            TIDE_OTIS_GLOBAL_ATLAS = 7,
     *                            TIDE_OTIS_PO = 8,
     *                            TIDE_OTIS_PO_ATLAS = 9
        TYPE TIDE_D
          CHARACTER(255)          PATH
c
c     1 for global NAO SSH,
c     2 for regional NAO around Japan (SSH only),
c     3 for regional NAO around Japan (SSH and UV),
c     4 for global LEGOS 0.25 degree SSH
c     5 for global LEGOS-FES2004 0.125 degree SSH
c
          INTEGER                 TIDE_MODEL
          INTEGER::               NUM_TH=0        !num of harmonics to estimate (input+by admittance)
          INTEGER::               ENABLE_ADM=1
          INTEGER::               NUM_TH_INPUT=0  !num of harmonics to enter
          INTEGER::               NUM_REGIONS=0
          TYPE(TOPO_G), POINTER:: GRID(:) => null()
          CHARACTER(4), POINTER:: HARMONIC_NAME(:) => null()
        END TYPE
C
C       FIXED DATA FILES extentions:
C                                 '.nao  ' FOR NAO-GLOBAL,
C                                 '_j.nao' FOR REGIONAL SSH,
C                                 'vfield.' file name FOR REGIONAL UV
C                                 '.legos2002' FOR LEGOS
C
        CHARACTER(*),PARAMETER::  DEPTH_FILE_NAOJ = 'bathy_j_5min.data'
        CHARACTER(17),PARAMETER:: FILE_EXT(TIDE_NUM_MODELS) =
     *   (/'.nao             ','_j.nao           ','vfield.          ',
     *     '.legos2002       ','_fes2004.asc     ',
     *     '_tpxo7.2         ','_tpxo7_atlas     ',
     *     '.PO_2009.out     ','.PO_2010atlas.out'/)
C
C     POINTERS TO REGIONS AND GRIDS ARE INTRODUCED
C     AS THE NUMBER OF REGIONS IS NOT FIXED - IT DEPENDS
C     FROM NUMBER OF OPEN BOUNDARY SEGMENTS AND THEIR LOCATION,
C     WHERE TIDES HAVE TO BE DEFINED. IF SPATIAL TIDAL FORCING
C     FOR ALL MODEL DOMAIN IS REQUESTED - SINGLE REGION OF MODEL
C     DOMAIN SIZE COULD BE THE BEST CHOISE.
C     SAME IS FOR HARMONICS NAMES.
C
C************  TABLE OF FIXED COEFFICIENTS FOR TIDAL FORMULAES **********C
C                                                                        C
C       NAO - regional and global models harmonic names
C       and LEGOS harmonics that could be downloaded
C       and simulated 'directly'
C                                                                        C
        INTEGER,PARAMETER::  NHARM_ALL = 22

        CHARACTER(4),DIMENSION(NHARM_ALL),PARAMETER:: TW_NAMES=(/
     *        'q1  ','o1  ','m1  ','p1  ','k1  ',
     *        'j1  ','oo1 ','2n2 ','mu2 ','n2  ',
     *        'nu2 ','m2  ','l2  ','t2  ','s2  ',
     *        'k2  ','mm  ','mf  ','ssa ','sa  ',
     *        'mtm ','msqm'/)
!       next list is used as equivalent to the above list
!       for file names generation only if fes2004 model is used
        CHARACTER(4),DIMENSION(NHARM_ALL),PARAMETER:: TW_NAMES_FES2004
     *     =(/'Q1  ','O1  ','M1  ','P1  ','K1  ',
     *        'J1  ','OO1 ','2N2 ','Mu2 ','N2  ',
     *        'Nu2 ','M2  ','L2  ','T2  ','S2  ',
     *        'K2  ','Mm  ','Mf  ','SSa ','Sa  ',
     *        'Mtm ','MSqm'/)
!
! with legos2004 need to extend number of harmonics:
! infer additional constituents by admittance
! Values from legos2004,
! duirnal from Richard Ray perth2 program
! semiduirnal (from N2-M2-K2) by spline admittances [see GRL 18[5]:845-848,1991]
!
        INTEGER,PARAMETER:: ADM_TD_I(3,NHARM_ALL) = RESHAPE(
     *    (/1,0,0, 2,0,0, 3,0,0, 4,0,0, 5,0,0,
     *      2,5,0, 2,5,0, !j1,oo1 => f(o1,k1)
     *      8,0,0,
     *      10,12,16,     !mu2 => f(n2,m2,k2)
     *      10,0,0,
     *      10,12,16,     !nu2 => f(n2,m2,k2)
     *      12,0,0,
     *      10,12,16,     !l2 => f(n2,m2,k2)
     *      10,12,16,     !t2 => f(n2,m2,k2)
     *      15,0,0, 16,0,0, 17,0,0, 18,0,0, 19,0,0, 20,0,0,
     *      21,0,0, 22,0,0/), (/3,NHARM_ALL/))

        DOUBLE PRECISION,PARAMETER:: ADM_TD_C(3,NHARM_ALL) = RESHAPE(
     *    (/1.,0.,0., 1.,0.,0., 1.,0.,0., 1.,0.,0., 1.,0.,0.,
     *      -0.0389,0.0836,0., -0.0431,0.0613,0., !j1,oo1
     *      1.,0.,0.,
     *      0.351535557706,-0.046278307672,0.069439968323,  !mu2 => f(n2,m2,k2)
     *      1.,0.,0.,
     *      0.156878802427,0.006755704028,-0.006104695053,  !nu2 => f(n2,m2,k2)
     *      1.,0.,0.,
     *      -0.051653455134,0.027869916824,0.077137765667,  !l2 => f(n2,m2,k2)
     *      -0.020101177502,0.008331518844,0.180480173707,  !t2 => f(n2,m2,k2)
     *      1.,0.,0., 1.,0.,0., 1.,0.,0., 1.,0.,0., 1.,0.,0., 1.,0.,0.,
     *      1.,0.,0., 1.,0.,0./), (/3,NHARM_ALL/))

!        INTEGER,PARAMETER:: NHARM_LEGOS2004 = 33
!        INTEGER,DIMENSION(NHARM_LEGOS2004):: ACTIVE_LEGOS2004 =
!     *  (//)
!       Up to 15 could be entered from LEGOS2004 model,
!       plus 18 short-term harmonics could be estimated by addmittance
!       Q1 O1 K1 2N2 N2 M2 K2 S2 P1 M4 S1 Mf Mm Mtm MSqm
!
!        INTEGER,DIMENSION(NHARM_ALL),PARAMETER:: CNUM_IN_LEGOS2004 =
!     *(/0,1,25,8,2,31,32,3,16,4,15,5,17,18,7,6,12,11,-1,-1,13,14/)
!
!       Are some of long-period waves directly entered?
!       if zeroes - estimate equilibrium tide, -1 -only this harmonic:
!                                  Sa, Ssa, Mm, Mf, Mtm, Msqm
!        integer:: LP_DEFINED(6) = (/-1,0,0,0,0,0/)
!         INTEGER,DIMENSION(NHARM_ALL),PARAMETER:: NUM_IN_LP_DEFINED =
!     *  (/-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,3,4,2,1,5,6/)

C
C*************************************************************************
c
C        COEFFICIENTS FOR NODAL FACTORS F AND U,
C        AI - FOR COS(N), BI - FOR SIN(N), Japan Tidal Tables, 1992
C        200811 - corrected by Edition 1.0 of IH[ydrographic]O
C        Harmonic Constituents with Nodal Corrections
C
        DOUBLE PRECISION,DIMENSION(8),PARAMETER::
C
C        TYPE::  'o1 '   'k1 '   'j1 '   'oo1'   'm2 '   'k2 '
c                'mm '   'mf '
C
     *   A1 = (/ 1.0176, 1.0060, 1.1029, 1.1027, 1.0004, 1.0246,
     *           1.0,    1.084/),
     *   A2 = (/ 0.1871, 0.1150, 0.1676, 0.6504,-0.0373, 0.2863,
     *          -0.13,   0.415/),
     *   A3 = (/-0.0147,-0.0088,-0.0170, 0.0317, 0.0002, 0.0083,
     *           0.0013, 0.039/),
     *   A4 = (/ 0.0014, 0.0006, 0.0016,-0.0014, 0.0,   -0.0015,
     *           0.0,    0.0/),
     *   B1 = (/10.80,  -8.86, -12.94, -36.68,  -2.14, -17.74,
     *           0.0,  -23.7/),
     *   B2 = (/-1.34,   0.68,   1.34,   4.02,   0.0,    0.68,
     *           0.0,    2.7/),
     *   B3 = (/ 0.19,  -0.07,  -0.19,  -0.57,   0.0,   -0.04,
     *           0.0,   -0.4/)
C
C        NODAL FACTOR COEFFICIENT'S NUMBER TO USE FOR EACH HARMONIC
C        NEGATIVE VALUES CORRESPOND TO SPECIAL CASES:
C            -1 :: M1 TYPE
C            -2 :: L2 TYPE
C            -3 :: F=1, U=0
C            -5 :: f=F(5), u=-U(5) [M2]
C
        INTEGER,DIMENSION(NHARM_ALL),PARAMETER:: I_NOD =
!  q1 o1 m1 p1 k1 j1 oo1 2n2 mu2 n2 nu2 m2 l2 t2 s2 k2 mm mf ssa sa mtm msqm
     *  (/1,1,-1,-3,2,3,4,5,5,5,5,5,-2,-3,-3,6,7,8,-3,-3,7,-5/)
C
C       COEFFICIENTS FOR V0
C
        DOUBLE PRECISION,DIMENSION(NHARM_ALL),PARAMETER::
     *TIDE_TYPE=
     * (/ 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2,2,2, 0,0,0,0, 0, 0/),
     *C2=
     * (/-3,-2,-1, 0, 0, 1, 2,-4,-4,-3,-3,-2,-1, 0,0,0, 1,2,0,0, 3, 4/),
     *C3=
     * (/ 1, 1, 1,-1, 1, 1, 1, 2, 4, 2, 4, 2, 2,-1,0,2, 0,0,2,1, 0,-2/),
     *C4=
     * (/ 1, 0, 0, 0, 0,-1, 0, 2, 0, 1,-1, 0,-1, 0,0,0,-1,0,0,0,-1, 0/),
     *C5=
     * (/270,270,90,270,90,90,90,0,0,0,0,0,180,283,0,0, 0,0,0,0, 0, 0/)
C
C       EQUILIBRIUM TIDE HEIGHT
C       (FROM A.SCHILLER (2004), OCEAN MODELLING (6) PP.31-49,
C       AFTER SCHWIDERSKI, 1980)
C
c        REAL, DIMENSION(NHARM_ALL),PARAMETER:: H_EQ_CM =
c             Q1      O1       M1       P1       K1       J1
c     *    ( / 1.9256, 10.0574, 0.0,     4.6843,  14.1565, 0.0,
c             OO1     2N2      MU2      N2       NU2      M2
c     *        0.0,    0.0,     0.0,     4.6398,  0.0,     24.2334,
c             L2      T2       S2       K2
c     *        0.0,    0.0,     11.2841, 3.0704/)
c
c       From Kantha & Clayson, 2000, based on Dietrich et al. (1980)
c                                    and Desai (1996):
        REAL, DIMENSION(NHARM_ALL),PARAMETER:: H_EQ_CM =
c          Q1      O1       M1       P1       K1       J1
     *  (/1.9387, 10.1266, 0.7965, 4.7129,  14.2408, 0.7965,
c          OO1     2N2      MU2      N2       NU2      M2
     *    0.4361, 0.6184,  0.7463,  4.6735,  0.8877,  24.4102,
c          L2      T2       S2       K2       Mm        Mf
     *    0.6899, 0.6636,  11.3572, 3.0875,  2.2191,  4.2017,
C          Ssa       Sa      Mtm       MSqm
     *    1.9542,   0.3104, 0.80425, 0.12931/)
C
C       EQUILIBRIUM TIDE HEIGHT FROM TAMURA POTENTIAL AMPLITUDE
C       REFERRING TO EPOCH 2000. Translatin value used: 26.70374427
C       TAMURA, Y. 1987: A HARMONIC DEVELOPMENT OF THE TIDE
C       GENERATING POTENTIAL. BULLETIN D'INFORMATIONS
C       MAREES TERRESTRES NO. 99, BRUXELLES 1987.
C
        REAL,PARAMETER:: TP2CM = 26.70374427
        REAL, DIMENSION(NHARM_ALL),PARAMETER:: TP_T87 =
C          Q1        O1        M1        P1        K1        J1
     *  (/0.072136, 0.376763, 0.029631, 0.175307, 0.529876, 0.029630,
C          OO1       2N2       MU2       N2        NU2       M2
     *    0.016212, 0.023009, 0.027768, 0.173881, 0.033027, 0.908184,
C          L2        T2        S2        K2        Mm        Mf
     *    0.025670, 0.024701, 0.422535, 0.114860, 0.082569, 0.156303,
C          Ssa       Sa        Mtm       MSqm
     *    0.072732, 0.011549, 0.029926, 0.004779/)
C
C      Eq of motion must replace: eta -> eta_geocentral - H_EQ,
c      that with known loading tide (inv to eta) could be:
c
c      (eta-eta_load) - alpha*H_EQ = beta*eta - alpha*H_EQ
c
c      beta is used if is unknown loading tide in solid earth;
c      alpha counts for results of potential increase due to the
c      elastic solid earth tide and self-attraction together with
c      decrease caused by resulting displacement from the earth center.
c
c
c       From Kantha & Clayson, 2000, based on Dietrich et al. (1980)
c                                    and Desai (1996):
        REAL, DIMENSION(NHARM_ALL),PARAMETER:: TP_ALPHA =
C             Q1      O1       M1       P1       K1       J1
     *   (/   0.6946, 0.695,   0.6962,  0.7059,  0.7364,  0.6911,
C             OO1     2N2      MU2      N2       NU2      M2
     *        0.6925, 0.693,   0.693,   0.693,   0.693,   0.693,
C             L2      T2       S2       K2       Mm       Mf
     *        0.693,  0.693,   0.693,   0.693,   0.693,   0.693,
C             Ssa     Sa        Mtm       MSqm
     *        0.693,  0.693,   0.693,   0.693/)

        REAL, DIMENSION(NHARM_ALL),PARAMETER:: TP_BETA =
C             Q1      O1       M1       P1       K1       J1
     *   (/   0.940, 0.940,    0.940,  0.940,  0.940,  0.940,
C             OO1     2N2      MU2      N2       NU2      M2
     *        0.940, 0.953,   0.953,   0.953,   0.953,   0.953,
C             L2      T2       S2       K2       Mm       Mf
     *        0.953,  0.953,   0.953,   0.953,   0.953,   0.953,
C             Ssa     Sa        Mtm       MSqm
     *        0.953,  0.953,   0.953,  0.953/)

      INTEGER,PRIVATE:: I

      DOUBLE PRECISION,PARAMETER::  SIGMA(NHARM_ALL) = (/(
     *        TIDE_TYPE(I)*15D0  + C2(I)*0.54901652D0 +
     *        C3(I)*0.04106864D0 + C4(I)*0.00464181D0,I=1,NHARM_ALL)/)
C
C**********************************************************************
C
      INTEGER, SAVE:: NUM_IN_USER(NHARM_ALL) = (/(0,i=1,NHARM_ALL)/)  !DEFINED FOR EACH IN NHARM_ALL
      INTEGER, SAVE:: NUM_IN_TMOD(NHARM_ALL) = (/(0,i=1,NHARM_ALL)/)  !DEFINED FOR EACH IN NUM_TH
      INTEGER, SAVE:: NUM_IN_INPUT(NHARM_ALL) = (/(0,i=1,NHARM_ALL)/) !DEFINED FOR EACH IN NHARM_ALL
      INTEGER, SAVE:: EQV_TIDE(NHARM_ALL) = (/(1,i=1,NHARM_ALL)/)     !DEFINED FOR EACH IN NHARM_ALL

      REAL, ALLOCATABLE, SAVE::   TIDE_CT(:,:,:,:,:,:),
     *                            TIDE_X(:,:,:), TIDE_Y(:,:,:)
      INTEGER,SAVE::              TIDE_ISTART = 0
      INTEGER,SAVE::              TIDE_NXM = 0, TIDE_NYM = 0

      PRIVATE      !! everything local except ...
      PUBLIC :: TOPO_G, TIDE_D, NHARM_ALL, TW_NAMES,
     *          ALLOC_TIDE_GRID, INIT_SL_TIDAL_HARMONICS,
     *          GET_SL_TIDAL_HARMONICS, FREE_TIDE_HARMONICS,
     *          INIT_TIDAL_EQ_PARAM, GET_TIDAL_EQ_PARAM,FREE_TIDE,
     *          ALLOC_TIDE_HARMONICS,TIDE_NUM_MODELS,TIDE_NAO99_GLOBAL,
     *          TIDE_NAO99_JAPAN_SSH,TIDE_NAO99_JAPAN_SSH_UVH,
     *          TIDE_LEGOS_2002,H_EQ_CM,TP_BETA,TP_ALPHA,TIDE_TYPE,
     *          NUM_IN_TMOD, NUM_IN_USER, NUM_IN_INPUT,TIDE_FES2004,
     *          EQV_TIDE,INIT_EQV_TIDE,
     *          TIDE_OTIS_GLOBAL, TIDE_OTIS_GLOBAL_ATLAS,
     *          TIDE_OTIS_PO, TIDE_OTIS_PO_ATLAS


      CONTAINS

C**********************************************************************C
C
      SUBROUTINE INIT_EQV_TIDE( TD, IRET, LOG )
C
c       Call it first if want all user-specified harmonics to be used.
c       If tides are requested on input stage list could be modified.
C       If list of harmonics is empty - use all available to estimate
C       equivalent tide, else only specified.
C
      IMPLICIT NONE
      TYPE(TIDE_D),     INTENT(INOUT):: TD
      INTEGER,          INTENT(OUT  ):: IRET
      INTEGER,OPTIONAL, INTENT(IN   ):: LOG

      INTEGER NF_LOG, NUM_TH, I, J, M

      IRET = 0
      IF(PRESENT(LOG))THEN
        NF_LOG = LOG
      ELSE
        NF_LOG = 0
      ENDIF

      NUM_TH = TD%NUM_TH
      EQV_TIDE = 1
      IF(NUM_TH <= 0)RETURN !save default units
      EQV_TIDE = 0
      DO M = 1, NUM_TH
        J=0
        DO I = 1, NHARM_ALL
          IF(TRIM(TD%HARMONIC_NAME(M)) == TRIM(TW_NAMES(I))) J = I
        ENDDO
        IF(J == 0)THEN
          IF(NF_LOG > 0 )WRITE(NF_LOG,*)
     *      ' *** ERROR, HARMONIC NAME NOT FROM TMOD, IGNORE: ',
     *      TRIM(TD%HARMONIC_NAME(M))
          CYCLE
        ENDIF
        EQV_TIDE(J) = 1
      ENDDO
      IF(MAXVAL(EQV_TIDE)==0)THEN
        IRET = 1
        IF(NF_LOG > 0 )WRITE(NF_LOG,*)
     *' * ERROR EQV_TIDE: NO REQUESTED HARMONIC(S) IN TMOD'
      ENDIF

      END SUBROUTINE INIT_EQV_TIDE
C
C**********************************************************************C

C**********************************************************************C
C
      SUBROUTINE INIT_TIDAL_EQ_PARAM( TD, IRET, LOG, IO )
C
C**********************************************************************C
C
C    1) DEFINE THE CORRESPONDENCE OF REQUESTED TIDAL HARMONICS WITH
C       HARMONICS AVAILABLE FROM THE TIDAL DATA FILES
C    2) IF NO HARMONICS SPECIFIED - USE ALL AVAILABLE
C    3) CHECK THE PRESENCE OF ORIGINAL TIDAL DATA FILES
C
C     IRET = 0 - NORMAL EXECUTION, ANOTHER VALUE - ERROR                 C
C
C    SERGEY M. VARLAMOV, 1991-2009
C    MODIFICATION - JUNE 2010
C                                                                      C
C**********************************************************************C
      IMPLICIT NONE
C
C    ARGUMENTS FOR ALL ENTRIES::
C
      TYPE(TIDE_D),     INTENT(INOUT):: TD
      INTEGER,          INTENT(OUT  ):: IRET
      INTEGER,OPTIONAL, INTENT(IN   ):: LOG, IO
C
C     LOCAL AUTOMATIC VARIABLES::
C
      INTEGER   I, J, K, M, IERR, NF_LOG, IMODEL, I_ADM, NUM_TH,
     *          N_SYM_EXT, N_SYM_PATH, NHI, NF_IO
      LOGICAL   EXST
      CHARACTER FILE_NAME*255
      CHARACTER(4),DIMENSION(NHARM_ALL):: NAMES_IN_OTIS_FILE
      integer no,mo,nco
      real theta_lim(2),phi_lim(2)

      IRET = 0

      IF(PRESENT(LOG))THEN
        NF_LOG = LOG
      ELSE
        NF_LOG = 0
      ENDIF

      IF(PRESENT(IO))THEN
        NF_IO = IO
      ELSE
        NF_IO = 81
      ENDIF
C
C     IF TIDE HARMONICS ARE NOT DEFINED - USE ALL AVAILABLE::
C
      NUM_TH = TD%NUM_TH

      IF(NUM_TH <= 0)THEN
        NUM_TH = NHARM_ALL
        CALL ALLOC_TIDE_HARMONICS(TD, NUM_TH, IERR)

        IF(IERR .NE. 0)THEN
          IRET = 61
          IF(NF_LOG > 0 )WRITE(NF_LOG,*)
     *      ' *** PROBLEMS ALLOCATING MEMORY BY '//
     *      'ALLOC_TIDE_HARMONICS'
          RETURN
        ENDIF

        DO M = 1, NUM_TH
          TD%HARMONIC_NAME(M) = TW_NAMES(M)
        ENDDO
      ENDIF
      IF(NUM_TH < NHARM_ALL)THEN
        IF(TD%ENABLE_ADM > 0)THEN
          IF(NF_LOG > 0 )WRITE(NF_LOG,*)
     *      ' *** REQUESTED LIMITED NUMBER OF HARMONICS, '//
     *      'DISABLE ESTIMATION BY ADMITTANCE: ',NUM_TH,NHARM_ALL
          TD%ENABLE_ADM = 0
        ENDIF
      ENDIF
c
c     even if model has u&v data (type==2), check ssh files only
c
      IMODEL = TD%TIDE_MODEL
      IF(TD%TIDE_MODEL == TIDE_NAO99_JAPAN_SSH_UVH)IMODEL=IMODEL-1

      N_SYM_EXT = LEN_TRIM(FILE_EXT(IMODEL))
C
C     ON SOME SYSTEMS EMPTY STRINGS ARE INITIALIZED BY "0"s,
C     AND TRIM RETURNS FULL STRING.
C
      N_SYM_PATH = LEN_TRIM(TD%PATH)
      IF(N_SYM_PATH > 0 .AND. ICHAR(TD%PATH(1:1)) == 0) N_SYM_PATH=0
      IF(N_SYM_PATH > 255 - N_SYM_EXT - 1)THEN
        IRET = 3
        IF(NF_LOG > 0 )WRITE(NF_LOG,*)
     *    ' *** TOO LONG PATH TO TIDAL DATA:', TD%TIDE_MODEL,
     *    TRIM(TD%PATH)
        RETURN
      ENDIF
      IF(N_SYM_PATH == 0)TD%PATH = "." !LOCAL DIRECTORY

      NUM_IN_TMOD = 0
      NUM_IN_USER = 0
      NUM_IN_INPUT = 0
      NHI = 0

      if(TD%TIDE_MODEL >= TIDE_OTIS_GLOBAL)then
!
!       single file keeps all harmonics
!
        IF(IMODEL==TIDE_OTIS_PO .or. IMODEL==TIDE_OTIS_PO_ATLAS)THEN
          FILE_NAME = TRIM(TD%PATH) // '/' //
     *    'hf' // TRIM(FILE_EXT(IMODEL))
        ELSE
          FILE_NAME = TRIM(TD%PATH) // '/' //
     *    'h' // TRIM(FILE_EXT(IMODEL))
        ENDIF
C
C       CHECK THAT FILE EXISTS
C
        INQUIRE(FILE = FILE_NAME, EXIST = EXST)
        IF(.NOT. EXST)THEN
          IRET = 69
          IF(NF_LOG > 0)WRITE(NF_LOG,*)
     *    ' *** DO NOT FOUND HARMONICS FILE: ',
     *    TRIM(FILE_NAME)," FOR OTIS MODEL ",TD%TIDE_MODEL
          RETURN
        endif

#ifdef CONVERT2BIGENDIAN
        OPEN(UNIT=NF_IO,FILE=FILE_NAME,STATUS='OLD',FORM='UNFORMATTED',
     *    ACTION='READ', convert='big_endian', IOSTAT=IERR)
#else
        OPEN(UNIT=NF_IO,FILE=FILE_NAME,STATUS='OLD',FORM='UNFORMATTED',
     *    ACTION='READ', IOSTAT=IERR)
#endif

        IF(IERR.NE.0) THEN
          IRET = 67
          IF(NF_LOG > 0 )WRITE(NF_LOG,*)
     *      ' *** CAN NOT OPEN TIDAL DATA FILE:', TRIM(FILE_NAME)
          RETURN
        ENDIF
        read(NF_IO) no,mo,nco
        if(nco <1 .or. nco > NHARM_ALL)then
          IRET = 67
          IF(NF_LOG > 0 )WRITE(NF_LOG,*)
     *      ' *** Number of harmonics in otis file ', TRIM(FILE_NAME),
     *      'exceeds limit: ',nco, NHARM_ALL
          RETURN
        endif
        rewind(NF_IO)
        read(NF_IO) no,mo,nco,theta_lim,phi_lim,
     *              NAMES_IN_OTIS_FILE(1:nco)
        close(NF_IO)
        NAMES_IN_OTIS_FILE(nco+1:NHARM_ALL)="    "
      endif

      DO M = 1, NUM_TH
C
C       CREATE THE TIDAL DATA FILE NAME FOR HARMONIC REQUESTED
C
        J=0
        DO I = 1, NHARM_ALL
          IF(TRIM(TD%HARMONIC_NAME(M)) == TRIM(TW_NAMES(I))) J = I
        ENDDO
        IF(J == 0)THEN
          IF(NF_LOG > 0 )WRITE(NF_LOG,*)
     *      ' *** ERROR, HARMONIC NAME NOT FROM TMOD, IGNORE: ',
     *      TRIM(TD%HARMONIC_NAME(M))
          CYCLE
        ENDIF
        if(TD%TIDE_MODEL >= TIDE_OTIS_GLOBAL)then
          J = 0
#ifdef _TDTEST_
          write(*,*)"Look for <",TRIM(TD%HARMONIC_NAME(M)),">"
#endif
          DO I = 1, nco
#ifdef _TDTEST_
          write(*,*)"Check ",I,J," <",TRIM(NAMES_IN_OTIS_FILE(I)),">"
#endif
            IF(TRIM(TD%HARMONIC_NAME(M)) ==
     *         TRIM(NAMES_IN_OTIS_FILE(I))) J = I
          ENDDO

          IF(J == 0)THEN
            IF(NF_LOG > 0 )WRITE(NF_LOG,*)
     *    ' *** WARNING, HARMONIC NAME DO NOT FOUND IN OTIS, IGNORE: ',
     *        TRIM(TD%HARMONIC_NAME(M))
            CYCLE
#ifdef _TDTEST_
          else
            write(*,*)"<",TRIM(TD%HARMONIC_NAME(M)), "> in OTIS ",J
#endif
          ENDIF
          EXST = .TRUE.
        else

          IF(IMODEL == TIDE_FES2004)THEN
            FILE_NAME = TRIM(TD%PATH) // '/' //
     *      TRIM(TW_NAMES_FES2004(J)) // TRIM(FILE_EXT(IMODEL))
          ELSE
            FILE_NAME = TRIM(TD%PATH) // '/' //
     *      TRIM(TD%HARMONIC_NAME(M)) // TRIM(FILE_EXT(IMODEL))
          ENDIF
C
C         CHECK THAT FILE EXISTS
C
          INQUIRE(FILE = FILE_NAME, EXIST = EXST)
        endif
        IF(.NOT. EXST)THEN
          NUM_IN_INPUT(M) = 0 !IT WAS REQUESTED BUT TRESE IS NO THIS FILE
          IF(NF_LOG > 0 )WRITE(NF_LOG,*)
     *      ' *** DO NOT FOUND TIDAL DATA, IGNORE: ',
     *      TRIM(TD%HARMONIC_NAME(M)),", ",TRIM(FILE_NAME)
        ELSE
          NHI = NHI + 1
          NUM_IN_INPUT(M) = NHI
        ENDIF
      ENDDO

      DO M = 1, NUM_TH
        J=0
        DO I = 1, NHARM_ALL
          IF(TRIM(TD%HARMONIC_NAME(M)) == TRIM(TW_NAMES(I))) J = I
        ENDDO
C
C       CHECK FOR DUBLICATE NAMES
C
        K = 0
        DO I = 1,M-1
          IF(IABS(NUM_IN_TMOD(I)) == J)K=1
        ENDDO

        IF(J > 0 .AND. K == 0)THEN
          NUM_IN_USER(J) = M
          IF(NUM_IN_INPUT(M) == 0 )THEN
            IF(TD%ENABLE_ADM == 0)THEN
              NUM_IN_TMOD(M) = 0
              NUM_IN_USER(J) = 0
            ELSE
              NUM_IN_TMOD(M) = -J
            ENDIF
          ELSE
            NUM_IN_TMOD(M) = J
          ENDIF
        ELSE
          IRET = 65
          IF(NF_LOG > 0)WRITE(NF_LOG,*)
     *    ' *** NON-EXISTING OR DUBLICATE TMOD HARMONIC IS REQUESTED: ',
     *    TRIM(TD%HARMONIC_NAME(M)),M,J,I
          RETURN
        ENDIF
      ENDDO
C
C     TEST POSSIBILITY TO USE ADMITTANCE
C
      IF( TD%ENABLE_ADM /= 0 )THEN
        DO M = 1, NUM_TH
          J = NUM_IN_TMOD(M)
          IF(J < 0)THEN
C
C           COULD BE ESTIMATED BY ADMITTANCE ?
C
            J = -J
            NUM_IN_TMOD(M) = 0
            I = ADM_TD_I(1,J)
            IF(I /= J)THEN !POINTS ON OTHER THEN THIS HARMONIC
              K = NUM_IN_USER(I)
              IF (K > 0)THEN
                IF(NUM_IN_INPUT(K) > 0)THEN !WAS ENTERED, CONTINUE
                  I = ADM_TD_I(2,J) !IT IS > 0, DO NOT CHECK
                  K = NUM_IN_USER(I)
                  IF (K > 0)THEN
                    IF(NUM_IN_INPUT(K) > 0)THEN
                      I = ADM_TD_I(3,J)
                      IF(I > 0)THEN      !NEED 3RD BASE HARMONIC
                        K = NUM_IN_USER(I)
                        IF (K > 0)THEN
                          IF(NUM_IN_INPUT(K) > 0)THEN
                           NUM_IN_TMOD(M) = -J
                          ENDIF
                        ENDIF
                      ELSE
                        NUM_IN_TMOD(M) = -J
                      ENDIF
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF
            ENDIF
          ENDIF
        ENDDO
      ENDIF
C
C     COMPACT WORKING ARRAYS
C
      J = NUM_TH
      I_ADM = 0
      M = 1
      DO WHILE(M <= J)
        IF(NUM_IN_TMOD(M) == 0)THEN
          J = J - 1
          IF(J <= 0)EXIT
          DO I = M,J
            TD%HARMONIC_NAME(I)=TD%HARMONIC_NAME(I+1)
            NUM_IN_TMOD(I) = NUM_IN_TMOD(I+1)
            NUM_IN_INPUT(I) = NUM_IN_INPUT(I+1)
          ENDDO
!         clean tale
          TD%HARMONIC_NAME(J+1) = " "
          NUM_IN_TMOD(J+1) = 0
          NUM_IN_INPUT(J+1) = 0

          DO K = 1,NHARM_ALL
            IF(NUM_IN_USER(K) == M)THEN
              NUM_IN_USER(K) = 0
            ELSEIF(NUM_IN_USER(K) > M)THEN
              NUM_IN_USER(K) = NUM_IN_USER(K)-1
            ENDIF
          ENDDO
          CYCLE
        ELSE IF(NUM_IN_TMOD(M) < 0)THEN
          I_ADM = I_ADM + 1
        ENDIF
        M = M + 1
      ENDDO
      NUM_TH = J
      IF(NUM_TH < 1)THEN
        IF(NF_LOG > 0 )WRITE(NF_LOG,*)
     *  ' *** DO NOT FOUND ANY INPUT TIDAL DATA FILES!'
        IRET = -2
        TD%NUM_TH = 0
        TD%NUM_TH_INPUT = 0
        RETURN
      ELSE
        TD%NUM_TH = NUM_TH
        TD%NUM_TH_INPUT = NUM_TH - I_ADM
      ENDIF

      IF(NF_LOG > 0 )WRITE(NF_LOG,*)
     *' *** FOUND ',TD%NUM_TH," DATA FILES WITH TIDAL HARMONICS"

      RETURN

      END SUBROUTINE  INIT_TIDAL_EQ_PARAM
C
C**********************************************************************C
C
      SUBROUTINE GET_TIDAL_EQ_PARAM(TD, GLB_UTC_TIME, PH, F,
     *            PHA, FA)
C
C**********************************************************************C
C     DEFINE TIME DEPENDANT TIDE PARAMETERS
C
C    1) ESTIMATE TIDAL EQUATION  NODAL PARAMETERS. IT IS DONE
C       AFTER TIME CHANGES FOR MORE THEN MAX_UPDATE_PHASE_TIME_INTERVAL SEC
C    2) CALCULATE THE QUICK PHASE OF ALL REQUESTED TIDAL HARMONICS

C     NODAL_FACTOR == F
C     PHASE_WITH_NODAL CORRECTION  == PH
C
C     H_TIDE(I) = F(I)*TH_AMP(I)*COS((PH(I)-TH_PHS(I))*RAD_IN_D)
C
C     WHERE I IS HARMONIC NUMBER, TH_AMP(I) AND TH_PHS(I) ARE
C     CONSTANT TIDAL HARMONIC AMPLITUDE AND PHASE IN THE FIXED POINT
C     UNDER CONSIDERATION
C
C*********************************************************************
        USE COM_SERV, ONLY:TDATE,TSEC,TDT,HOUR,RAD_IN_D,D_IN_RAD,DAY
        IMPLICIT NONE
        REAL, INTENT(INOUT)::             F(*),
     *                                    PH(*),
     *                                    FA(NHARM_ALL),
     *                                    PHA(NHARM_ALL)
        TYPE(TIDE_D),     INTENT(IN)::    TD
        DOUBLE PRECISION, INTENT(IN   ):: GLB_UTC_TIME
C
C       FIXED PARAMETERS::
C
        REAL,            PARAMETER::  EPS = 1E-6
        DOUBLE PRECISION,PARAMETER::  MAX_UPDATE_PHASE_TIME_INTERVAL =
     *                                DAY/48D0
C
C     LOCAL SAVED VARIABLES::
C
        DOUBLE PRECISION, SAVE::      V0(NHARM_ALL),
     *                                U(NHARM_ALL)
        DOUBLE PRECISION,SAVE::       LAST_UPDATE_PHASE_TIME =
     *                                    -999D0*DAY,
     *                                REF_TIME = -999D99
C
C     LOCAL AUTOMATIC VARIABLES::
C
        INTEGER           I,J,M,NUM_TH
        DOUBLE PRECISION  T, D, Y, LL, S, H, P, N, N2, N3, P2,
     *                    PP, FCOSU_L2, FSINU_L2, F_L2, V, UU,
     *                    U_L2, F_M1, U_M1, FCOSU_M1, FSINU_M1,
     *                    F_MM, U_MM, COSN, COS2N, COS3N, COSP2mN,
     *                    SINN, SIN2N, SIN3N
        TYPE (TDATE)::    DT
        DOUBLE PRECISION  UPDATE_PHASE_PASSED_TIME,
     *                    START_YEAR_TIME

        NUM_TH=TD%NUM_TH

        UPDATE_PHASE_PASSED_TIME =
     *    GLB_UTC_TIME-LAST_UPDATE_PHASE_TIME

        IF(UPDATE_PHASE_PASSED_TIME >
     *     MAX_UPDATE_PHASE_TIME_INTERVAL)THEN

          LAST_UPDATE_PHASE_TIME = GLB_UTC_TIME
C
C         CURRENT UTC CALENDAR DATE::
C
          DT = TDT(GLB_UTC_TIME)
C
C         START OF DAY GLOBAL TIME IN SEC::
C
          DT%IHOUR=0
          DT%IMIN=0
          DT%ISEC=0
          REF_TIME = TSEC(DT)
C----------------------------------------------------------
C         NODAL PHASE CORRECTION
C----------------------------------------------------------
          Y=DT%IYEAR
C
C         START OF YEAR
C
          DT%IMONTH=1
          DT%IDAY=1
          START_YEAR_TIME = TSEC(DT)

          D = INT((REF_TIME - START_YEAR_TIME)/DAY)
C
C         NUMBER OF LEAP YEARS SINCE 2000 (NEGATIVE, IF BEFORE 2000):
C
          LL= INT((Y+3)/4)-500

          S = 211.728D0+129.38471D0*(Y-2000D0)+13.176396D0*(D+LL)
          S = MOD(S,360D0)
          H = 279.974D0-0.23871D0*(Y-2000D0)+0.985647D0*(D+LL)
          H = MOD(H,360D0)
          P = 83.298D0+40.66229D0*(Y-2000D0)+0.111404D0*(D+LL)
          P = MOD(P,360D0)
          N = 125.071D0-19.32812D0*(Y-2000D0)-0.052954D0*(D+LL)

          N  = MOD(N,360D0)
          N2 = MOD(N*2D0,360D0)*RAD_IN_D
          N3 = MOD(N*3D0,360D0)*RAD_IN_D
          N  = N*RAD_IN_D

          PP = MOD(P,360D0)*RAD_IN_D
          P2 = MOD(P*2D0,360D0)*RAD_IN_D

          COSN = COS(N)
          COS2N = COS(N2)
          COS3N = COS(N3)
          COSP2mN = COS(P2-N)
          SINN = SIN(N)
          SIN2N = SIN(N2)
          SIN3N = SIN(N3)

          FCOSU_L2 = 1. - 0.2505D0*COS(P2) - 0.1102D0*COSP2mN -
     *                    0.0156D0*COS(P2-N2) - 0.0370D0*COSN
          FSINU_L2 = 0. - 0.2505D0*SIN(P2) - 0.1102D0*SIN(P2-N) -
     *                    0.0156D0*SIN(P2-N2) - 0.0370D0*SINN
          F_L2 = SQRT(FCOSU_L2*FCOSU_L2 + FSINU_L2*FSINU_L2)

          IF(ABS(FCOSU_L2) < EPS)THEN
            IF(ABS(FSINU_L2) < EPS)THEN
              F_L2=0D0
              U_L2=0D0
            ELSEIF(FSINU_L2 <= -EPS)THEN
              U_L2=-90D0
            ELSE
              U_L2=90D0
            ENDIF
          ELSE
            U_L2 = ATAN(FSINU_L2/FCOSU_L2)*D_IN_RAD
            IF(FCOSU_L2 <= -EPS)U_L2=U_L2+180D0
          ENDIF

          FCOSU_M1 = 2D0*COS(PP) + 0.4D0*COS(PP-N)
          FSINU_M1 =    SIN(PP) + 0.2D0*SIN(PP-N)
          F_M1 = SQRT(FCOSU_M1*FCOSU_M1 + FSINU_M1*FSINU_M1)
!         definition from Harm. Constants Specification, 2006
          F_MM = 1D0-0.1311*COSN+0.0538*COS(P2)+0.0205*COSP2mN
          U_MM = 0D0

          IF(ABS(FCOSU_M1) < EPS)THEN
            IF(ABS(FSINU_M1) < EPS)THEN
              F_M1=0D0
              U_M1=0D0
            ELSEIF(FSINU_M1 <= -EPS)THEN
              U_M1=-90D0
            ELSE
              U_M1=90D0
            ENDIF
          ELSE
            U_M1 = ATAN(FSINU_M1/FCOSU_M1)*D_IN_RAD
            IF(FCOSU_M1 <= -EPS)U_M1=U_M1+180D0
          ENDIF

          DO M = 1,NHARM_ALL
            V = C2(M)*S + C3(M)*H + C4(M)*P+C5(M)
            V0(M) = MOD(V, 360D0)
            I = I_NOD(M)
            IF(I==-5)I=5
            IF(I > 0)THEN
              FA(M)=A1(I) + A2(I)*COSN + A3(I)*COS2N +
     *             A4(I)*COS3N
              UU = B1(I)*SINN+B2(I)*SIN2N+B3(I)*SIN3N
              U(M)=MOD(UU,360D0)
            ELSEIF(I == -3)THEN
              FA(M)=1D0
              U(M)=0D0
            ELSEIF(I == -1)THEN
C
C             M1 TYPE
C
              FA(M)=F_M1
              U(M)=U_M1
            ELSEIF(I == -2)THEN
C
C             L2 TYPE
C
              FA(M)=F_L2
              U(M)=U_L2
C            ELSE !COMMENT FOR PARALLEL
CC
CC              ERROR
CC
C               IRET = 71
C               IF(NF_LOG > 0)WRITE(NF_LOG,*)
C     *           ' *** INCORRECT TYPE FOR NODAL FACTORS',
C     *           TD%HARMONIC_NAME(M)
C               RETURN
            ENDIF
            IF(I_NOD(M)==-5)U(M)=-U(M)

          ENDDO !HARMONICS CYCLE
        ENDIF !END SLOW PHASE PARAMETERS DEFINITION
C----------------------------------------------------------
C       TOTAL PHASE CALCULATION, IS DONE IN ANY CASE
C       TO UPDATE THE FORMAL ARGUMENT IF IT WAS MODIFIED BY
C       USER OUTSIDE OF SUBROUTINE:
C----------------------------------------------------------
C
C       TIME IN HOURS PASSED AFTER THE REFERENCE DATE/TIME
C
        T = (GLB_UTC_TIME - REF_TIME)/HOUR

        DO M=1,NHARM_ALL
          V = MOD(SIGMA(M)*T, 360D0)
          PHA(M) = MOD(V0(M)+U(M)+V, 360D0)
        ENDDO

!        WRITE(*,'(a,(10f8.2))')'F: u=',u
!        WRITE(*,'(a,(10f8.2))')'F: V0=',V0
!        WRITE(*,'(a,(10f8.2))')'F: F=',FA
!        WRITE(*,'(a,(10f8.2))')'F: F=',F_MM

!        write(*,*)"Tide, define F, NUM_TH=",NUM_TH

        IF(NUM_TH >= 1)THEN
          DO M = 1,NUM_TH
            J = IABS(NUM_IN_TMOD(M))
            F(M)=FA(J)
            PH(M) = PHA(J)
          ENDDO
        ENDIF
        RETURN

      END SUBROUTINE GET_TIDAL_EQ_PARAM
C
C**********************************************************************C

C**********************************************************************C
C
      SUBROUTINE FREE_TIDE_HARMONICS(TD, LOG, FREE_ALL)
C
C**********************************************************************C
        IMPLICIT NONE
        TYPE(TIDE_D), INTENT(INOUT)::   TD
        INTEGER,OPTIONAL,INTENT(IN)::   LOG, FREE_ALL
        INTEGER IT,I
        IT = 0
        IF(ALLOCATED(TIDE_CT))THEN
          DEALLOCATE(TIDE_CT, STAT = I)
          IF(I .NE. 0)IT=IT+1
        ENDIF
        IF(ALLOCATED(TIDE_X))THEN
          DEALLOCATE(TIDE_X, STAT = I)
          IF(I .NE. 0)IT=IT+2
        ENDIF
        IF(ALLOCATED(TIDE_Y))THEN
          DEALLOCATE(TIDE_Y, STAT = I)
          IF(I .NE. 0)IT=IT+4
        ENDIF
        IF(PRESENT(FREE_ALL))THEN
          IF(FREE_ALL > 0)THEN
            IF(ASSOCIATED(TD%GRID))THEN
              DEALLOCATE(TD%GRID, STAT = I)
              NULLIFY(TD%GRID)
              IF(I .NE. 0) IT=IT+8
            ENDIF
            IF(ASSOCIATED(TD%HARMONIC_NAME))THEN
              DEALLOCATE(TD%HARMONIC_NAME, STAT = I)
              NULLIFY(TD%HARMONIC_NAME)
              IF(I .NE. 0)IT=IT+16
            ENDIF
          ENDIF
        ENDIF
        IF(IT .NE. 0)THEN
          IF(PRESENT(LOG))THEN
            IF(LOG > 0) WRITE (UNIT=LOG,FMT=*)
     *      ' *** DEALLOCATION ERROR IN FREE_TIDE_HARMONICS:',IT
            ENDIF
        ENDIF

      END SUBROUTINE FREE_TIDE_HARMONICS
C
C**********************************************************************C

C**********************************************************************C
C
      SUBROUTINE FREE_TIDE(TD, LOG)
C
C**********************************************************************C
        IMPLICIT                      NONE
        TYPE(TIDE_D),     INTENT(INOUT):: TD
        INTEGER,OPTIONAL, INTENT(IN   ):: LOG
        IF(PRESENT(LOG))THEN
          CALL FREE_TIDE_HARMONICS(TD, LOG = LOG, FREE_ALL = 1)
        ELSE
          CALL FREE_TIDE_HARMONICS(TD, FREE_ALL = 1)
        ENDIF
      END SUBROUTINE FREE_TIDE
C
C**********************************************************************

C**********************************************************************C
C
      SUBROUTINE INIT_SL_TIDAL_HARMONICS
     * ( TD, IRET, IO, LOG, DEBUG, IO2, VOLUME_FLUX)
C                                                                      C
C**********************************************************************C
C
C    1) OPEN TIDAL DATA FILES
C    2) DEFINE THE TIDAL MODEL GRID AND
C    3) POSITION OF REQUESTED SUB-GRIDS
C    4) IF NO SUB-GRIDS SPECIFIED - TAKE ALL TIDAL MODEL GRID AS SINGLE
C       SUB-GRID
C    5) ALLOCATE THE SUB-GRIDS SPECIFICATION ARRAYS X&Y, WORKING ARRAYS,
C       AND INTERPOLATION COEFFICIENT ARRAYS
C    6) ENTER THE TIDAL HARMONICS AMPLITUDE AND PHASE FROM
C       TIDAL SEA LEVEL DATA FILE
C    7) ESTIMATE THE INTERPOLATION COEFFICIENTS FOR THEM
C       (SIN/COS TRANSFORMATION IS USED FOR PHASE)
C
C     HARMONIC NAME MUST BE SAME AS USED BY MATSUMOTO, IN LOWER CASE
C     TRAILING BLANKS ARE IMPORTANT
C
C     IRET = 0 - NORMAL EXECUTION, ANOTHER VALUE - ERROR
C
C    MODIFICATION - JULY, 2001
C    MODIFICATION - JANUARY 2005
C    MODIFICATION - JANUARY 2010
C    MODIFICATION - JUNE 2010
C
C    SERGEY M. VARLAMOV
C
C**********************************************************************C
C
      USE COM_SERV, ONLY: RAD_IN_D, EXTRAPOLATE_BY_MASK, GET_BSPLINES_S
      IMPLICIT NONE
C
C     ARGUMENTS
C
      TYPE(TIDE_D), INTENT(INOUT)::   TD
      INTEGER,      INTENT(OUT  )::   IRET
      INTEGER,OPTIONAL,INTENT(IN)::   IO, IO2, LOG, DEBUG,
     *                                VOLUME_FLUX
C
C     PARAMETERS
C
      REAL,PARAMETER::            EPS = 1E-3
C
C     LOCAL
C
      REAL, ALLOCATABLE:: MASK(:,:), UV(:,:,:), MASK_F(:,:,:)
      DOUBLE PRECISION, ALLOCATABLE:: AMP(:,:,:),  PHS(:,:,:),
     * MDLDPTH(:,:)

      INTEGER  NUM_TH, NUM_REGIONS
      INTEGER  NX, NY, NX_FULL, NY_FULL, K_MIN, K_MAX,
     *         NH, NR, NP, I, J, K, INEG, IERR,
     *         NXYM, NXC, NF_IO, NF_LOG, NF_IO2, NHI,
     *         IT, GLOBX, II, JJ, IDEBUG, IVOLUME_FLUX, NUM_TH_INPUT,
     *         NPARAM,IM,JM,IM1,JM1,IMODEL
      REAL DX, DY, ADX, ADY, DDX, DDY
      REAL XX, YY, XMINF, XMAXF, YMINF, D,
     *     OZ, OZ_MAX, CR, CCR, A, P, X0, Y0

      DOUBLE PRECISION:: X0_FULL, Y0_FULL, Y1_FULL, X1_FULL, UNDEF,
     *                   DXX, DYY, RMDLDPTN
      CHARACTER::        FILE_NAME*255
      CHARACTER::        FILE_TOPO*255
C
C     INITIALIZE ONCE ONLY, SINGLE INITIALIZATION COPY ALLOWED::
C
      IRET = 0
      IF(PRESENT(LOG))THEN
        NF_LOG = LOG
      ELSE
        NF_LOG = 0
      ENDIF

      IF(TD%TIDE_MODEL < 1 .OR. TD%TIDE_MODEL > TIDE_NUM_MODELS)THEN
        IRET = 1
        IF(NF_LOG > 0 )WRITE(NF_LOG,*)
     *    ' *** UNSUPPORTED TIDAL MODEL NUMBER:',TD%TIDE_MODEL,
     *    TRIM(TD%PATH)
        RETURN
      ENDIF

      IF(PRESENT(VOLUME_FLUX))THEN
        IVOLUME_FLUX = VOLUME_FLUX
      ELSE
        IVOLUME_FLUX = 0
      ENDIF
      IF(PRESENT(DEBUG))THEN
        IDEBUG = DEBUG
      ELSE
        IDEBUG = 0
      ENDIF
      IF(PRESENT(IO))THEN
        NF_IO = IO
      ELSE
        NF_IO = 81
      ENDIF
      IF(PRESENT(IO2))THEN
        NF_IO2 = IO2
      ELSE
        NF_IO2 = 83
      ENDIF

      NPARAM=1
      IF(TD%TIDE_MODEL == TIDE_NAO99_JAPAN_SSH_UVH .or.
     *   TD%TIDE_MODEL >= TIDE_OTIS_GLOBAL) NPARAM = 3

      NUM_TH      = TD%NUM_TH
      NUM_REGIONS = TD%NUM_REGIONS

!      write(*,*)"Tide initialization: num_th,NUM_REGIONS=",
!     *           NUM_TH,NUM_REGIONS
      NUM_TH_INPUT = TD%NUM_TH_INPUT
      NHI = 0
      DO NH = 1, NUM_TH
        IF(NUM_IN_TMOD(NH) < 0)THEN
          CYCLE
        ENDIF
        NHI = NHI + 1
C
C       CREATE THE TIDAL DATA FILE NAME FOR HARMONIC REQUESTED
C
c       if model has u&v data (type==2), first get ssh grid and data:
c
        IMODEL=TD%TIDE_MODEL
        IF(TD%TIDE_MODEL == TIDE_NAO99_JAPAN_SSH_UVH)IMODEL=IMODEL-1

        if(IMODEL >= TIDE_OTIS_GLOBAL)then
!
!       single file keeps all harmonics
!
          IF(IMODEL==TIDE_OTIS_PO .or. IMODEL==TIDE_OTIS_PO_ATLAS)THEN
            FILE_NAME = TRIM(TD%PATH) // '/' //
     *      'hf' // TRIM(FILE_EXT(IMODEL))
          ELSE
            FILE_NAME = TRIM(TD%PATH) // '/' //
     *      'h' // TRIM(FILE_EXT(IMODEL))
          ENDIF
#ifdef CONVERT2BIGENDIAN
          OPEN(UNIT=NF_IO,FILE=FILE_NAME,STATUS='OLD',
     *      FORM='UNFORMATTED', ACTION='READ',
     *      CONVERT='big_endian', IOSTAT=IERR)
#else
          OPEN(UNIT=NF_IO,FILE=FILE_NAME,STATUS='OLD',
     *      FORM='UNFORMATTED', ACTION='READ', IOSTAT=IERR)
#endif
        ELSE
          IF(iMODEL == TIDE_FES2004)THEN
            FILE_NAME = TRIM(TD%PATH) // '/' //
     *      TRIM(TW_NAMES_FES2004(NUM_IN_TMOD(NH))) //
     *      TRIM(FILE_EXT(IMODEL))
          ELSE
            FILE_NAME = TRIM(TD%PATH) // '/' //
     *      TRIM(TD%HARMONIC_NAME(NH)) // TRIM(FILE_EXT(IMODEL))
          ENDIF
          OPEN(UNIT=NF_IO,FILE=FILE_NAME,STATUS='OLD',
     *      FORM='FORMATTED', ACTION='READ', IOSTAT=IERR)
        ENDIF
        IF(IERR.NE.0) THEN
          IRET = 5
          IF(NF_LOG > 0 )WRITE(NF_LOG,*)
     *      ' *** CAN NOT OPEN TIDAL DATA FILE:', TRIM(FILE_NAME)
          RETURN
        ENDIF

        IF(NHI == 1)THEN
C
C         DEFINE THE TIDAL MODEL GRID PARAMETERS
C         FROM FIRST HARMONIC USED::
C
          IF(TD%TIDE_MODEL <= TIDE_NAO99_JAPAN_SSH_UVH)THEN
             CALL DEF_NAO_GRID(NF_IO, DXX, DYY, X0_FULL, X1_FULL,
     *         Y0_FULL, Y1_FULL, GLOBX, NX_FULL, NY_FULL, IERR)
          ELSE IF(TD%TIDE_MODEL == TIDE_LEGOS_2002)THEN
            CALL DEF_LEGOS_GRID(DXX, DYY, X0_FULL, X1_FULL,
     *        Y0_FULL, Y1_FULL, GLOBX, NX_FULL, NY_FULL)
          ELSE IF(TD%TIDE_MODEL == TIDE_FES2004)THEN
            CALL DEF_FES2004_GRID(NF_IO, DXX, DYY, X0_FULL, X1_FULL,
     *        Y0_FULL, Y1_FULL, GLOBX, NX_FULL, NY_FULL, IERR)
          ELSE IF(IMODEL >= TIDE_OTIS_GLOBAL .and.
     *            IMODEL <= TIDE_NUM_MODELS)then
            CALL DEF_OTIS_GRID(NF_IO, DXX, DYY, X0_FULL, X1_FULL,
     *        Y0_FULL, Y1_FULL, GLOBX, NX_FULL, NY_FULL, IERR)
          ELSE
            IERR = -1
            IF(NF_LOG > 0 )WRITE(NF_LOG,*)
     *      ' *** UNKNOWN TIDAL MODEL TYPE ',TD%TIDE_MODEL
          ENDIF
          IF( IERR .NE. 0)THEN
            IRET = 7
            IF(NF_LOG > 0 )WRITE(NF_LOG,*)
     *      ' *** ERROR READING TIDE FILE HEADER FROM ',
     *      TRIM(FILE_NAME)
            RETURN
          ENDIF
C
C         CHECK THE GRID CONSISTENCY, ::
C
          IF( DABS(Y1_FULL-Y0_FULL - DYY*(NY_FULL-1)) > EPS .OR.
     *        DABS(X1_FULL-X0_FULL - DXX*(NX_FULL-1)) > EPS )THEN
            IRET = 9
            IF(NF_LOG > 0 )then
              WRITE(NF_LOG,*)
     *        ' *** ERROR IN TIDE GRID DEFINITION FROM '//
     *        'FILE HEADER', IERR, TRIM(FILE_NAME)
              write(nf_log,*)Y1_FULL,Y0_FULL,DYY,NY_FULL
              write(nf_log,*)x1_FULL,x0_FULL,Dxx,Nx_FULL,globx
            endif
            RETURN
          ENDIF

#ifdef _TDTEST_
          WRITE(*,*)"GRID: ", TRIM(FILE_NAME)
          write(*,*)Y1_FULL,Y0_FULL,DYY,NY_FULL
          write(*,*)x1_FULL,x0_FULL,Dxx,Nx_FULL,globx
#endif
          IF(NUM_REGIONS <= 0) THEN
C
C           USER DID NOT SPECIFIED, USE (SINGLE) ALL TIDAL MODEL REGION::
C
            NUM_REGIONS = 1
            CALL ALLOC_TIDE_GRID(TD, NUM_REGIONS, IERR)
            IF(IERR .NE. 0)THEN
              IRET = 11
              IF(NF_LOG > 0 )WRITE(NF_LOG,*)
     *        ' *** PROBLEMS ALLOCATING MEMORY BY '//
     *        'ALLOC_TIDE_GRID FOR MAX REGION'
              RETURN
            ENDIF
C
C           SPECIFY REQUEST AS FOR ALL ORIGINAL DATA
C
            TD%GRID(1)%XMIN = MIN(X0_FULL,X1_FULL)
            TD%GRID(1)%XMAX = MAX(X0_FULL,X1_FULL)
            TD%GRID(1)%YMIN = MIN(Y0_FULL,Y1_FULL)
            TD%GRID(1)%YMAX = MAX(Y0_FULL,Y1_FULL)

          ENDIF

          DX = DXX
          DY = DYY
          ADX = ABS(DX)
          ADY = ABS(DY)
          DDX = EPS*ADX
          DDY = EPS*ADY

          DO NR = 1, NUM_REGIONS
C
C           COPY TO TIDE INFO STRUCTURE::
C
            TD%GRID(NR)%X0_FULL = X0_FULL
            TD%GRID(NR)%X1_FULL = X1_FULL
            TD%GRID(NR)%Y0_FULL = Y0_FULL
            TD%GRID(NR)%Y1_FULL = Y1_FULL
            TD%GRID(NR)%DX = DX
            TD%GRID(NR)%DY = DY
            TD%GRID(NR)%NX_FULL = NX_FULL
            TD%GRID(NR)%NY_FULL = NY_FULL
            TD%GRID(NR)%GLOBX = GLOBX

          ENDDO

!      write(*,*)"Tide Grid Before Adj: ",1,TD%GRID(1)
C
C         ANALYZE THE REGION REQUESTS RELATIVE TO DATA REGION
C         all are done for the ssh grid only,
c         ignore here shifted u,v grids (in this version)
C
          DO NR = 1, NUM_REGIONS
            IF(GLOBX > 0)THEN
C
C             IT IS GLOBAL X MODEL, ANY X REQUEST IS 'OK',
C             JUST ALIGN REQUEST WITH GRID - SHIFT 'OUTSIDE'
C             RESULT COULD BE MORE THEN 360 DEGREES ZONE - OK!
C
              INEG = 0
              D = TD%GRID(NR)%XMIN - X0_FULL
              IF(D < 0.)INEG = 1
              K_MIN = INT(D/ADX)-1-INEG
              TD%GRID(NR)%XMIN = X0_FULL + ADX*K_MIN
              INEG = 0
              D = TD%GRID(NR)%XMAX - X0_FULL
              IF(D < 0.)INEG = 1
              K_MAX = INT(D/ADX)+2-INEG
              TD%GRID(NR)%XMAX = X0_FULL + ADX*K_MAX

            ELSE
C
C             LIMITED COVERAGE DATA. REQUEST IS ASSUMED TO BE IN
C             'NORMAL ORDER', INPUT DATA COULD BE IN ANY, RE-ORDER
C
              IF(X0_FULL > X1_FULL)THEN
                XMINF = X1_FULL
                XMAXF = X0_FULL
              ELSE
                XMINF = X0_FULL
                XMAXF = X1_FULL
              ENDIF
C
C             FIND CYCLIC (360 DEGREE) SHIFT OF INPUT DATA
C             THAT WILL RESULT IN LARGEST OVERLAPPING
C
              OZ_MAX = -1.
              CCR = 0.
              DO I = 1,11
                CR = 360.*(I-6)
                OZ = MIN(XMAXF+CR,TD%GRID(NR)%XMAX) -
     *               MAX(XMINF+CR,TD%GRID(NR)%XMIN)
                IF(OZ > OZ_MAX)THEN
                  OZ_MAX = OZ
                  CCR = CR
                ENDIF
              ENDDO
C              IF(OZ_MAX <= 2.*ADX)THEN ! REQUESTED AREA SMALL, DISABLE LATER
C                IRET = -3
C                RETURN
C              ENDIF

              XMINF = XMINF + CCR
              INEG = 0
              D = TD%GRID(NR)%XMIN - XMINF
              IF(D < 0.)INEG = 1
              K_MIN = INT(D/ADX)-1-INEG
              IF(K_MIN < 0)K_MIN = 0
              IF(K_MIN > NX_FULL-1)K_MIN = NX_FULL-1
              TD%GRID(NR)%XMIN = XMINF + ADX*K_MIN

!              write(*,*)"X:OZ_MAX,CCR,INEG,K_MIN=",OZ_MAX,CCR,INEG,K_MIN

              INEG = 0
              D = TD%GRID(NR)%XMAX - XMINF
              IF(D < 0.)INEG = 1
              K_MAX = INT(D/ADX)+2-INEG
              IF(K_MAX < 0)K_MAX = 0
              IF(K_MAX > NX_FULL-1)K_MAX = NX_FULL-1
              TD%GRID(NR)%XMAX = XMINF + ADX*K_MAX

!              write(*,*)"X,INEG,K_MAX=",INEG,K_MAX

            ENDIF
C
C     ************ Y DOMAIN **************
C
            IF(Y0_FULL > Y1_FULL)THEN
              YMINF = Y1_FULL
            ELSE
              YMINF = Y0_FULL
            ENDIF

            INEG = 0
            D = TD%GRID(NR)%YMIN - YMINF
            IF(D < 0.)INEG = 1
            K_MIN = INT(D/ADY)-1-INEG
            IF(K_MIN < 0)K_MIN = 0
            IF(K_MIN > NY_FULL-1)K_MIN = NY_FULL-1
            TD%GRID(NR)%YMIN = YMINF + ADY*K_MIN

!            write(*,*)"Y:INEG,K_MIN=",INEG,K_MIN

            INEG = 0
            D = TD%GRID(NR)%YMAX - YMINF
            IF(D < 0.)INEG = 1
            K_MAX = INT(D/ADY)+2-INEG
            IF(K_MAX < 0)K_MAX = 0
            IF(K_MAX > NY_FULL-1)K_MAX = NY_FULL-1
            TD%GRID(NR)%YMAX = YMINF + ADY*K_MAX

!            write(*,*)"Y:INEG,K_MAX=",INEG,K_MAX

            TD%GRID(NR)%NX =
     *        NINT((TD%GRID(NR)%XMAX-TD%GRID(NR)%XMIN)/ADX)+1
            TD%GRID(NR)%NY =
     *        NINT((TD%GRID(NR)%YMAX-TD%GRID(NR)%YMIN)/ADY)+1

            IF(TD%GRID(NR)%NX < 4 .OR. TD%GRID(NR)%NY < 4 )THEN
C
C             SMALL REGION, SPLINE INTERPOLATION METHOD
C             CAN NOT BE APPLIED::
C
              IF(NF_LOG > 0 )WRITE(NF_LOG,*)
     *        ' *** TIDAL DATA DO NOT COVER WELL SUBREGION ', NR,
     *        TD%GRID(NR)
              IRET = 15
              RETURN
            ENDIF
            TIDE_NXM = MAX(TD%GRID(NR)%NX, TIDE_NXM )
            TIDE_NYM = MAX(TD%GRID(NR)%NY, TIDE_NYM )
          ENDDO ! END OF NR REGIONS ANALYSIS

!      write(*,*)"Tide Grid AFTER Adj: ",1,TD%GRID(1)

          NXYM    = MAX(TIDE_NXM, TIDE_NYM)

!      write(*,*)"Tide param allocation:",TIDE_NXM,TIDE_NYM,NPARAM,NUM_TH,
!     * NUM_REGIONS,NX_FULL,NY_FULL,NXYM,2*TIDE_NXM*TIDE_NYM+2*NXYM+1
C
C         ALLOCATE THE WORKING ARRAYS
C
          ALLOCATE(
     *      TIDE_CT(2, TIDE_NXM, 2*TIDE_NYM, NPARAM*2,
     *              NUM_TH_INPUT, NUM_REGIONS),
     *      TIDE_X(TIDE_NXM, NUM_REGIONS,NPARAM),
     *      TIDE_Y(TIDE_NYM, NUM_REGIONS,NPARAM),
     *      MASK(TIDE_NXM,TIDE_NYM),
     *      UV(TIDE_NXM,TIDE_NYM,2),
     *      AMP(NX_FULL,NY_FULL,NPARAM),
     *      PHS(NX_FULL,NY_FULL,NPARAM),
     *      MDLDPTH(NX_FULL,NY_FULL),
     *      MASK_F(NX_FULL,NY_FULL,NPARAM),
     *      STAT = IERR)
          IF( IERR .NE. 0)THEN
            IRET = 17
            IF(NF_LOG > 0 )WRITE(NF_LOG,*)
     *      ' *** ERROR ALLOCATING MEMORY FOR TIDAL INFO ARRAYS'
            RETURN
          ENDIF

          TIDE_CT=0.
          TIDE_X=0.
          TIDE_Y=0.

          DO NR = 1, NUM_REGIONS
            X0 = TD%GRID(NR)%XMIN
            Y0 = TD%GRID(NR)%YMIN
            NX = TD%GRID(NR)%NX
            NY = TD%GRID(NR)%NY

            DO I = 1, NX
              TIDE_X(I, NR, 1:NPARAM) = X0 + (I-1)*ADX
            ENDDO
            DO J = 1, NY
              TIDE_Y(J, NR, 1:NPARAM) = Y0 + (J-1)*ADY
            ENDDO
            IF(NPARAM > 1)THEN
              TIDE_X(1:NX,NR,2) = TIDE_X(1:NX,NR,2)-0.5*ADX  !U
              TIDE_Y(1:NY,NR,3) = TIDE_Y(1:NY,NR,3)-0.5*ADY  !V

!       write(*,*)'grid 1:',TIDE_X(1,NR,1),TIDE_X(NX,NR,1),TIDE_Y(1,NR,1),TIDE_Y(NY,NR,1)
!       write(*,*)'grid 2:',TIDE_X(1,NR,2),TIDE_X(NX,NR,2),TIDE_Y(1,NR,2),TIDE_Y(NY,NR,2)
!       write(*,*)'grid 3:',TIDE_X(1,NR,3),TIDE_X(NX,NR,3),TIDE_Y(1,NR,3),TIDE_Y(NY,NR,3)

            ENDIF
          ENDDO
C
C         CYCLIC FACTOR FOR GLOBAL MODEL
C
          NXC = NINT(360./ADX)

          IF(TD%TIDE_MODEL == TIDE_NAO99_JAPAN_SSH_UVH .AND.
     *       IVOLUME_FLUX > 0)THEN
C
C           ENTER NAOJ DEPTH IN CM FOR ESTIMATION TIDAL VOLUME FLUXES
C
            FILE_TOPO = TRIM(TD%PATH) // '/' //
     *      TRIM(DEPTH_FILE_NAOJ)

            OPEN(UNIT=NF_IO2,FILE=FILE_TOPO,STATUS='OLD',
     *      FORM='FORMATTED',ACTION='READ',IOSTAT=IERR)

            IF(IERR.NE.0) THEN
              IRET = 6
              IF(NF_LOG > 0 )WRITE(NF_LOG,*)
     *        ' *** CAN NOT OPEN NAOJ DEPTH FILE:', TRIM(FILE_TOPO)
              RETURN
            ENDIF

            MDLDPTH = 0D0

            CALL RD_NAO_MDLDPTH(NF_IO2, MDLDPTH, NX_FULL, NY_FULL,
     *        X0_FULL, Y0_FULL, DXX, DYY, IERR)

            CLOSE(UNIT=NF_IO2)

            IF( IERR .NE. 0)THEN
              IRET = 20
              IF(NF_LOG > 0 )WRITE(NF_LOG,*)
     *        ' *** ERROR READING TIDAL INFO: ',IERR, TRIM(FILE_TOPO)
              RETURN
            ENDIF
          ENDIF ! DEPTH INPUT FOR NAOJ
        ENDIF ! END OF NH == 1, GRID DEFINITION AND ARRAY ALLOCATIONS

        MASK_F = 0.
        AMP = 0D0
        PHS = 0D0
C
C       DATA INPUT FOR HARMONIC NH:
C       AMP IN CM, PHASE IN DEGREES (FOR GREENWICH MERIDIAN)
C
        IF(TD%TIDE_MODEL <= TIDE_NAO99_JAPAN_SSH_UVH)THEN
          CALL RD_NAO_MAP(NF_IO, AMP, PHS, NX_FULL, NY_FULL,
     *        UNDEF, IERR)
        ELSE IF(TD%TIDE_MODEL == TIDE_LEGOS_2002)THEN
          CALL RD_LEGOS_MAP(NF_IO, AMP, PHS, NX_FULL, NY_FULL,
     *        UNDEF, IERR)
        ELSE IF(TD%TIDE_MODEL == TIDE_FES2004)THEN
          CALL RD_FES2004_MAP(NF_IO, AMP, PHS, NX_FULL, NY_FULL,
     *        UNDEF, IERR)
        ELSE IF(TD%TIDE_MODEL >= TIDE_OTIS_GLOBAL .and.
     *          TD%TIDE_MODEL <= TIDE_NUM_MODELS)then
          CALL RD_H_OTIS(NF_IO, AMP, PHS, NX_FULL, NY_FULL,
     *        TD%HARMONIC_NAME(NH), UNDEF, IERR)
        ENDIF
!        write(*,*)"input done, close ",NF_IO,trim(FILE_NAME)
        CLOSE(NF_IO)

        IF( IERR .NE. 0)THEN
          IRET = 20
          IF(NF_LOG > 0 )WRITE(NF_LOG,*)
     *    ' *** ERROR READING TIDAL INFO: ',IERR, TRIM(FILE_NAME)
          RETURN
        ENDIF

        IF(IDEBUG > 0 .AND. NF_LOG > 0)THEN
          WRITE(NF_LOG,*)' === SSH TIDAL HARMONIC ENTERED: ',
     *    TRIM(TD%HARMONIC_NAME(NH)),", ",TRIM(FILE_NAME)
        ENDIF

        DO J = 1, NY_FULL
          DO I = 1, NX_FULL
            IF (AMP(I,J,1) < UNDEF)THEN
              IF(AMP(I,J,1) < 0.)THEN
                IF(NF_LOG > 0)WRITE(NF_LOG,*)
     *          ' *** WARNING: TIDE SSH AMPLITUDE IS NEGATIVE: ',
     *          AMP(I,J,1),I,J,NH,TRIM(FILE_NAME)
              ELSE
                MASK_F(I,J,1) = 1.
              ENDIF
            endif
          ENDDO  ! X COLUMNS CYCLE END
        ENDDO  ! Y LINES CYCLE END

        IF( IERR > 0)THEN
          IRET = 21
          IF(NF_LOG > 0 )WRITE(NF_LOG,*)
     *    ' *** NEGATIVE AMPLITUDE OF REQUESTED AREA: ',IERR
          RETURN
        ENDIF

        IF(TD%TIDE_MODEL == TIDE_NAO99_JAPAN_SSH_UVH)THEN
C
C         ENTER ADDITIONAL U&V INFORMATION
C
          FILE_NAME = TRIM(TD%PATH) // '/' //
     *    TRIM(FILE_EXT(TD%TIDE_MODEL)) //
     *    TRIM(TD%HARMONIC_NAME(NH))
C
C         OPEN THE FILE FOR READING DATA
C
          OPEN(UNIT=NF_IO,FILE=FILE_NAME,STATUS='OLD',FORM='FORMATTED',
     *      ACTION='READ', IOSTAT=IERR)

          IF(IERR.NE.0) THEN
            IRET = 6
            IF(NF_LOG > 0 )WRITE(NF_LOG,*)
     *      ' *** CAN NOT OPEN TIDAL U&V DATA FILE:', TRIM(FILE_NAME)
            RETURN
          ENDIF

          CALL RD_NAO_UV(NF_IO, AMP(1,1,2), PHS(1,1,2), MASK_F(1,1,2),
     *      NX_FULL, NY_FULL, X0_FULL, Y0_FULL, DXX, DYY, IERR)

          IF( IERR .NE. 0)THEN
            IRET = 7
            IF(NF_LOG > 0 )WRITE(NF_LOG,*)
     *    ' *** ERROR READING TIDAL INFO: ',IERR, TRIM(FILE_NAME)
            RETURN
          ENDIF

          CLOSE(UNIT=NF_IO)

          IF(IDEBUG > 0 .AND. NF_LOG > 0)THEN
            WRITE(NF_LOG,*)' === U&V TIDAL HARMONICS ENTERED: ',
     *      TRIM(TD%HARMONIC_NAME(NH)),", ",
     *      TRIM(FILE_NAME)
          ENDIF
          IF( IVOLUME_FLUX > 0 )THEN
C
C           ESTIMATE VOLUME FLUXES INSTEAD OF VELOCITY
c           USING NAO DEPTH DATA AND ASSUMING C-GRID
C
            IF(DX > 0)THEN
              IM=1
            ELSE
              IM=-1
            ENDIF
            IF(DY > 0)THEN
              JM=1
            ELSE
              JM=-1
            ENDIF

            DO J=1,NY_FULL
              JM1=J-JM
              IF(JM1 < 1)JM1=1
              IF(JM1 > NY_FULL)JM1=NY_FULL
              DO I=1,NX_FULL
                IM1=I-IM
                IF(IM1 < 1)IM1=1
                IF(IM1 > NX_FULL)IM1=NX_FULL
                IF(MASK_F(I,J,2) > 0.5)THEN
                  RMDLDPTN =0.5*(MDLDPTH(I,J)+MDLDPTH(IM1,J))
                  IF(RMDLDPTN < EPS)THEN
                    IF(NF_LOG > 0 )WRITE(NF_LOG,*)
     * ' = WARNING: TIDAL MODEL DEPTH ZERO AROUND U POINT',I,J
                  ENDIF
                  AMP(I,J,2)=AMP(I,J,2)*RMDLDPTN
                  IF(MASK_F(I,J,1) < EPS)THEN
                    IF(NF_LOG > 0 )WRITE(NF_LOG,*)
     * ' = WARNING: SSH UNDEFINED AROUND U POINT',I,J
                  ENDIF
                ENDIF
                IF(MASK_F(I,J,3) > 0.5)THEN
                  RMDLDPTN =0.5*(MDLDPTH(I,J)+MDLDPTH(I,JM1))
                  IF(RMDLDPTN < EPS)THEN
                    IF(NF_LOG > 0 )WRITE(NF_LOG,*)
     * ' = WARNING: TIDAL MODEL DEPTH ZERO AROUND V POINT',I,J
                  ENDIF
                  AMP(I,J,3)=AMP(I,J,3)*RMDLDPTN
                  IF(MASK_F(I,J,1) < EPS)THEN
                    IF(NF_LOG > 0 )WRITE(NF_LOG,*)
     * ' = WARNING: SSH UNDEFINED AROUND V POINT',I,J
                  ENDIF
                ENDIF
              ENDDO
            ENDDO
          ENDIF
        ENDIF

        IF(TD%TIDE_MODEL >= TIDE_OTIS_GLOBAL)then

          IF( IVOLUME_FLUX < 1 )THEN
            IRET = 7
            IF(NF_LOG > 0 )WRITE(NF_LOG,*)
     *      ' *** DO not realized yet input of u,v'//
     *      ' from OTIS (only UH, VH)'
            RETURN
          ENDIF
C
C         ENTER ADDITIONAL U&V INFORMATION,
!         single file keeps all harmonics
!
          IF(IMODEL==TIDE_OTIS_PO .or. IMODEL==TIDE_OTIS_PO_ATLAS)THEN
            FILE_NAME = TRIM(TD%PATH) // '/' //
     *        'uv' // TRIM(FILE_EXT(IMODEL))
          ELSE
            FILE_NAME = TRIM(TD%PATH) // '/' //
     *        'u' // TRIM(FILE_EXT(IMODEL))
          ENDIF
#ifdef CONVERT2BIGENDIAN
          OPEN(UNIT=NF_IO,FILE=FILE_NAME,STATUS='OLD',
     *      FORM='UNFORMATTED', ACTION='READ',
     *      CONVERT='big_endian', IOSTAT=IERR)
#else
          OPEN(UNIT=NF_IO,FILE=FILE_NAME,STATUS='OLD',
     *      FORM='UNFORMATTED', ACTION='READ', IOSTAT=IERR)
#endif
          IF(IERR.NE.0) THEN
            IRET = 6
            IF(NF_LOG > 0 )WRITE(NF_LOG,*)
     *      ' *** CAN NOT OPEN TIDAL U&V DATA FILE:', TRIM(FILE_NAME)
            RETURN
          ENDIF

          CALL RD_UV_OTIS(NF_IO, AMP(1,1,2), PHS(1,1,2), MASK_F(1,1,2),
     *      NX_FULL, NY_FULL, TD%HARMONIC_NAME(NH), IERR)

          IF( IERR .NE. 0)THEN
            IRET = 7
            IF(NF_LOG > 0 )WRITE(NF_LOG,*)
     *    ' *** ERROR READING TIDAL INFO: ',IERR, TRIM(FILE_NAME)
            RETURN
          ENDIF

          CLOSE(UNIT=NF_IO)

          IF(IDEBUG > 0 .AND. NF_LOG > 0)THEN
            WRITE(NF_LOG,*)' === U&V TIDAL HARMONICS ENTERED: ',
     *      TRIM(TD%HARMONIC_NAME(NH)),", ",TRIM(FILE_NAME)
          ENDIF
        ENDIF

!        IF(NPARAM > 1)THEN
!       write(*,*)"****** Outside of cycle ********"
!       nr=1
!       NX = TD%GRID(NR)%NX
!       NY = TD%GRID(NR)%NY
!       write(*,*)'grid 1:',TIDE_X(1,NR,1),TIDE_X(NX,NR,1),TIDE_Y(1,NR,1),TIDE_Y(NY,NR,1)
!       write(*,*)'grid 2:',TIDE_X(1,NR,2),TIDE_X(NX,NR,2),TIDE_Y(1,NR,2),TIDE_Y(NY,NR,2)
!       write(*,*)'grid 3:',TIDE_X(1,NR,3),TIDE_X(NX,NR,3),TIDE_Y(1,NR,3),TIDE_Y(NY,NR,3)
!        ENDIF

        DO NR = 1, NUM_REGIONS
          NX = TD%GRID(NR)%NX
          NY = TD%GRID(NR)%NY
          X0 = TD%GRID(NR)%XMIN
          Y0 = TD%GRID(NR)%YMIN

          DO NP = 1, NPARAM
C
C           GIVE INITIAL LAND VALUES TO MASK AND ZERO TO AMP/PHASE::
C
            MASK = 0
            UV   = 0
            IERR = 0
            DO J = 1, NY_FULL
              YY = Y0_FULL + (J-1)*DY !LATITUDE
              JJ = NINT((YY - Y0)/ADY) + 1
              IF( JJ > 0 .AND. JJ <= NY )THEN
                DO I = 1, NX_FULL
                  IF (MASK_F(I,J,NP) > 0.5)THEN
                    XX = X0_FULL + (I-1)*DX
                    A = AMP(I,J,NP)
                    P = PHS(I,J,NP)*RAD_IN_D
                    II = NINT(MOD((XX - X0),360.)/ADX) + 1
                    IT = II
                    DO WHILE(IT <= NX)
                      IF(IT > 0)THEN
c                        IF(A < 0. .and. np==1)IERR = IERR + 1
                        IF(A < 0.)IERR = IERR + 1
                        UV(IT,JJ,1) = A*SIN(P)
                        UV(IT,JJ,2) = A*COS(P)
c                        UV(IT,JJ,3) = A
                        MASK(IT,JJ) = 1.
                      ENDIF
                      IT = IT + NXC
                    ENDDO
                    IT = II - NXC  !IT ELIMINATES DUBLICATE ASSIGNMENT
                    DO WHILE(IT > 0)
                      IF(IT <= NX)THEN
c                        IF(A < 0. .and. np==1)IERR = IERR + 1
                        IF(A < 0. )IERR = IERR + 1
                        UV(IT,JJ,1) = A*SIN(P)
                        UV(IT,JJ,2) = A*COS(P)
c                        UV(IT,JJ,3) = A
                        MASK(IT,JJ) = 1.
                      ENDIF
                      IT = IT - NXC
                    ENDDO
                  ENDIF
                ENDDO  ! X COLUMNS CYCLE END
              ENDIF
            ENDDO  ! Y LINES CYCLE END

            IF( IERR > 0)THEN
              IRET = 21
              IF(NF_LOG > 0 )WRITE(NF_LOG,*)
     *        ' *** NEGATIVE AMPLITUDE OF REQUESTED AREA: ',IERR
              RETURN
            ENDIF
C
C           EXTRAPOLATE TIDE AMPLITUDE*sin:cos(PHASE) CLOSE TO THE LAND
C           IT NEEDS TO FILL SOME ISLAND POINTS IN ORIGINAL TIDAL MODEL,
C           WHICH ARE SEA POINTS IN USED O.C. MODEL
C
            DO K=1,2
              CALL EXTRAPOLATE_BY_MASK
     *        ( UV(1,1,K), MASK, NX, NY, TIDE_NXM, IERR )
              IF(IERR. NE. 0) THEN
                IRET = 23
                IF(NF_LOG > 0 )WRITE(NF_LOG,*)
     *          ' *** ERROR IN TIDE_HARMONICS -> MASK CORRECTION: ',
     *          IERR,K
                RETURN
              ENDIF
C
C             CALCULATE BICUBIC SPLINE COEFFICIENTS
C
!      write(*,*)'Do splines for K,NP,NR,NH=',K,NP,NR,NH,NX,NY,TIDE_NXM,
!     *(NP-1)*3+K
!      write(*,*)'grid:',K,TIDE_X(1,NR,NP),TIDE_X(NX,NR,NP),TIDE_Y(1,NR,NP),TIDE_Y(NY,NR,NP)

              CALL GET_BSPLINES_S
     *        ( UV(1,1,K), TIDE_X(1,NR,NP), NX, TIDE_Y(1,NR,NP), NY,
     *          TIDE_CT(1,1,1,(NP-1)*2+K,NHI,NR), TIDE_NXM, IERR)
              IF(IERR .NE. 0)THEN
                IRET = 20 +K*2
                IF(NF_LOG > 0 )WRITE (UNIT=NF_LOG,FMT=*)
     *          ' *** ERROR CALCULATING SPLINE COEFFICIENTS'//
     *          ' FOR TIDAL DATA, K,NP,NR,NH,NHI,ERR: ',
     *          K,NP,NR,NH,NHI,IERR
                RETURN
              ENDIF
!      write(*,*)'grid after:',K,TIDE_X(1,NR,NP),TIDE_X(NX,NR,NP),
!     *TIDE_Y(1,NR,NP),TIDE_Y(NY,NR,NP)
            ENDDO ! K: SIN-COS
          ENDDO ! NP: NPARAM
        ENDDO ! NR: NUM_REGIONS
      ENDDO !NH: NUM_HARMONICS

      TIDE_ISTART = 1

      IT = 0
      IF(ALLOCATED(MASK))THEN
        DEALLOCATE(MASK, STAT = I)
        IF(I .NE. 0)IT=IT+1
      ENDIF
      IF(ALLOCATED(UV))THEN
        DEALLOCATE(UV, STAT = I)
        IF(I .NE. 0)IT=IT+2
      ENDIF
      IF(ALLOCATED(AMP))THEN
        DEALLOCATE(AMP, STAT = I)
        IF(I .NE. 0)IT=IT+4
      ENDIF
      IF(ALLOCATED(PHS))THEN
        DEALLOCATE(PHS, STAT = I)
        IF(I .NE. 0)IT=IT+8
      ENDIF
      IF(ALLOCATED(MDLDPTH))THEN
        DEALLOCATE(MDLDPTH, STAT = I)
        IF(I .NE. 0)IT=IT+16
      ENDIF
      IF(ALLOCATED(MASK_F))THEN
        DEALLOCATE(MASK_F, STAT = I)
        IF(I .NE. 0)IT=IT+32
      ENDIF
      IF(IT .NE. 0)THEN
        IF(NF_LOG > 0) WRITE (UNIT=NF_LOG,FMT=*)
     *  ' *** DEALLOCATION ERROR IN INIT_SL_TIDAL_HARMONICS:',IT
      ENDIF
      RETURN
      END SUBROUTINE INIT_SL_TIDAL_HARMONICS
C
C**********************************************************************

C**********************************************************************C
C
      SUBROUTINE GET_SL_TIDAL_HARMONICS( NHRM, TD, PT, HAMP, HPHASE,
     *                              HAMP_SIN_PHASE, HAMP_COS_PHASE,
     *                              IRET, LOG, CUT_EXTREMUMS, PARAM)
C
C************************************************************************
c      PARAM = 1 is default for the ssh, 2 - u, 3 - v (if model has them)
C
C       ESTIMATE THE REQUESTED TIDAL HARMONIC AMPLITUDE AND
C       PHASE (DEGREES) IN POINT SPECIFIED
C
C************************************************************************
      USE COM_SERV, ONLY: TPOINT, D_IN_RAD, BCSEVL_EX
      IMPLICIT NONE
C
C     ARGUMENTS FOR ALL ENTRIES::
C
      TYPE(TPOINT), INTENT(IN   )::   PT
      TYPE(TIDE_D), INTENT(INOUT)::   TD
      REAL,         INTENT(OUT  )::   HAMP, HPHASE
      REAL,         INTENT(OUT  )::   HAMP_SIN_PHASE, HAMP_COS_PHASE
      INTEGER,      INTENT(IN   )::   NHRM
      INTEGER,      INTENT(OUT  )::   IRET
      INTEGER,OPTIONAL,INTENT(IN)::   LOG,CUT_EXTREMUMS, PARAM

      REAL,PARAMETER::            EPS = 1E-3

      INTEGER  NX, NY, NR, I, K, IERR, nr_best, NF_LOG, NHI, IH,
     *         ICUT_EXTREMUMS, IPARAM, JM
      REAL     ADX, ADY, DDX, DDY
      REAL     XX, YY, XM, A, X0, Y0, X1, Y1,
     *         dxy_best, xr_best, dx2, dy1, dy2, dxy

      IF(PRESENT(LOG))THEN
        NF_LOG = LOG
      ELSE
        NF_LOG = 0
      ENDIF

      IF(PRESENT(CUT_EXTREMUMS))THEN
        ICUT_EXTREMUMS = CUT_EXTREMUMS
      ELSE
        ICUT_EXTREMUMS = 0
      ENDIF

      IF(PRESENT(PARAM))THEN
        IPARAM = PARAM
      ELSE
        IPARAM = 1
      ENDIF

      IF(TIDE_ISTART <= 0)THEN
        IF(NF_LOG > 0)WRITE(NF_LOG,*)
     *  ' *** CALL GET_SL_TIDAL_HARMONICS BEFORE INITIALIZATION'
        IRET = 41
        RETURN
      ENDIF

      IF(IPARAM < 1 .or. IPARAM > 3 .or.
     *  (IPARAM > 1 .and.
     *   (TD%TIDE_MODEL /= TIDE_NAO99_JAPAN_SSH_UVH .and.
     *   (TD%TIDE_MODEL < TIDE_OTIS_GLOBAL))))THEN
        IF(NF_LOG > 0)WRITE(NF_LOG,*)
     *  ' *** GET_SL_TIDAL_HARMONICS, PARAM(1) NOT IN MODEL(2): ',
     *  IPARAM, TD%TIDE_MODEL
        IRET = 42
        RETURN
      ENDIF
C
C    CHECK THE INDEX NH OF REQUESTED HARMONIC
C
      IF(NHRM < 1 .OR. NHRM > TD%NUM_TH)THEN
        IF(NF_LOG > 0)WRITE(NF_LOG,*)
     *  ' *** INDEX OF REQUESTED TIDAL HARMONIC EXCEEDS LIMITS: ',
     *  NHRM,', MAX = ',TD%NUM_TH
        IRET = 43
        RETURN
      ENDIF

      XX = PT%X
      YY = PT%Y
C
C     WHAT REGION POINT BELONGS TO? Search "best" (05.2007)
c     For multiple parameters grids could be shifted, min-max are not valid
c
      nr_best = 0
      adx = abs(TD%GRID(1)%DX)
      ady = abs(TD%GRID(1)%DY)
      DDX = EPS*ADX
      DDY = EPS*ADY

      dxy_best = -(adx+ady)
      xr_best = dxy_best
      DO NR = 1, TD%NUM_REGIONS

        NX = TD%GRID(NR)%NX
        NY = TD%GRID(NR)%NY

        X0 = TIDE_X(1,NR,IPARAM)
        Y0 = TIDE_Y(1,NR,IPARAM)
        X1 = TIDE_X(NX,NR,IPARAM)
        Y1 = TIDE_Y(NY,NR,IPARAM)

        XM = MOD(XX - X0, 360.)
        DO WHILE(XM < -DDX)
          XM = XM + 360.
        ENDDO

        IF( XM >= -DDX          .AND.
     *      XM <= X1 + DDX - X0 .AND.
     *      YY >= Y0 - DDY      .AND.
     *      YY <= Y1 + DDY           )THEN
C
C         IT IS IN REGION NR, CORRECT POSITION AND EXIT CYCLE::
C         IF POINT WAS "VERY CLOSE" TO BOUNDARY FROM OUTSIDE -
C         SHIFT IT TO BOUNDARY::
C
          dx2=x1-x0-xm
          dy1=yy-y0
          dy2=y1-yy
          dxy=min(xm,dx2,dy1,dy2)
          if(dxy > dxy_best)then
            dxy_best=dxy
            nr_best=nr
            xr_best=xm
          endif
        ENDIF
      ENDDO

      IF( NR_best == 0 )THEN
        IF( NF_LOG > 0 )then
          WRITE (UNIT=NF_LOG,FMT=*)
     *    ' *** POINT TO DEFINE TIDE OUT OF SELECTED DOMAINS: ',
     *    XX,YY,TD%NUM_REGIONS
          WRITE (UNIT=NF_LOG,FMT=*)
     *      TD%GRID(1)%X0_FULL,TD%GRID(1)%X1_FULL,
     *      TD%GRID(1)%Y0_FULL,TD%GRID(1)%Y1_FULL,
     *      TD%GRID(1)%dx,TD%GRID(1)%dy,
     *      TD%GRID(1)%nx_FULL,TD%GRID(1)%ny_FULL
        endif
        IRET = 45
        RETURN
      ENDIF
C
C     IT IS IN REGION NR, CORRECT POSITION.
C     IF POINT WAS "VERY CLOSE" TO BOUNDARY FROM OUTSIDE -
C      SHIFT IT TO BOUNDARY::
C
      NR = nr_best
      NX = TD%GRID(NR)%NX
      NY = TD%GRID(NR)%NY

      X0 = TIDE_X(1,NR,IPARAM)
      Y0 = TIDE_Y(1,NR,IPARAM)
      X1 = TIDE_X(NX,NR,IPARAM)
      Y1 = TIDE_Y(NY,NR,IPARAM)

      XX = Xr_best + X0
      IF(XX < X0) XX = X0
      IF(XX > X1) XX = X1
      IF(YY < Y0) YY = Y0
      IF(YY > Y1) YY = Y1

      JM = IABS(NUM_IN_TMOD(NHRM))

      HAMP_SIN_PHASE = 0.
      HAMP_COS_PHASE = 0.

      DO I = 1,3
        IF(ADM_TD_I(I,JM) > 0)THEN
          IH = NUM_IN_USER(ADM_TD_I(I,JM))
          NHI = NUM_IN_INPUT(IH)
          DO K = 1,2
        CALL BCSEVL_EX(TIDE_X(1,NR,IPARAM),NX,TIDE_Y(1,NR,IPARAM),NY,
     *    TIDE_CT(1,1,1,(IPARAM-1)*2+K,NHI,NR),TIDE_NXM,XX,YY,A,
     *    IERR,ICUT_EXTREMUMS)

        IF(IERR .NE. 0)THEN
          IRET = 50+K*2
          IF(NF_LOG > 0 )THEN
             WRITE (UNIT=NF_LOG,FMT=*)
     *    ' *** ERROR IN SPLINE INTERPOLATION, TIDE: ', IERR,K
             WRITE (UNIT=NF_LOG,FMT=*)nx,ny,nr,iparam
             WRITE (UNIT=NF_LOG,FMT=*)TIDE_NXM,TIDE_NYM,TD%NUM_REGIONS
             WRITE (UNIT=NF_LOG,FMT=*)xx
             WRITE (UNIT=NF_LOG,FMT=*)yy
             WRITE (UNIT=NF_LOG,FMT=*)TIDE_X(1:nx,nr,iparam)
             WRITE (UNIT=NF_LOG,FMT=*)TIDE_Y(1:ny,nr,iparam)
          endif
          RETURN
        ENDIF
        IF(K==1)THEN
          HAMP_SIN_PHASE = HAMP_SIN_PHASE + A*ADM_TD_C(I,JM)
        ELSE
          HAMP_COS_PHASE = HAMP_COS_PHASE + A*ADM_TD_C(I,JM)
        ENDIF
          ENDDO
        ENDIF
      ENDDO

      HPHASE = ATAN2(HAMP_SIN_PHASE,HAMP_COS_PHASE)*D_IN_RAD
      HAMP = SQRT(HAMP_SIN_PHASE*HAMP_SIN_PHASE +
     *            HAMP_COS_PHASE*HAMP_COS_PHASE)

      RETURN
      END SUBROUTINE GET_SL_TIDAL_HARMONICS
C
C**********************************************************************

C**********************************************************************
C
      SUBROUTINE ALLOC_TIDE_GRID(TD, NR, MERR)
C
C**********************************************************************
        IMPLICIT                      NONE
        TYPE(TIDE_D), INTENT(INOUT):: TD
        INTEGER,      INTENT(IN   ):: NR
        INTEGER,      INTENT(OUT  ):: MERR
        IF(TD%NUM_REGIONS > 0)THEN
          MERR = -1
        ELSE IF(ASSOCIATED(TD%GRID)) THEN
          MERR = -2
        ELSE
          ALLOCATE(TD%GRID(NR), STAT=MERR)
          IF(MERR == 0) TD%NUM_REGIONS = NR
        ENDIF
      END SUBROUTINE ALLOC_TIDE_GRID
c
C**********************************************************************

C**********************************************************************
C
      SUBROUTINE ALLOC_TIDE_HARMONICS(TD, NH, MERR)
C
C**********************************************************************
        IMPLICIT                      NONE
        TYPE(TIDE_D), INTENT(INOUT):: TD
        INTEGER,      INTENT(IN   ):: NH
        INTEGER,      INTENT(OUT  ):: MERR

        IF(TD%NUM_TH > 0)THEN
          MERR = -1
        ELSE IF(ASSOCIATED(TD%HARMONIC_NAME))THEN
          MERR = -2
        ELSE
          ALLOCATE(TD%HARMONIC_NAME(NH), STAT = MERR)
          IF(MERR == 0) TD%NUM_TH = NH
        ENDIF
      END SUBROUTINE ALLOC_TIDE_HARMONICS
C
C**********************************************************************C

C----NAO INPUT --------------------------------------------------------
C
C MODIFICATION OF THE ORIGINAL PROGRAM NAO2XYAP OF
C KOJI MATSUMOTO, JAPAN NAO, 1999.06.08
C
C APPLICABLE FOR BOTH LOCAL AROUND JAPAN (1/12 DEGREE) AND
C GLOBAL (1/2 DEGREE) COMPONENT FILES
C DEFINE THE MODEL GRID FROM FILE HEADER,
C READ THE SEA LEVEL TIDAL HARMONICS AMPLITUDE AND PHAZE
C
C    VERSION 3.00
C
C    LAST MODIFICATION - JULY, 2001
C
C    SERGEY M. VARLAMOV
C
C$MAIN
C--------------------------------------------------------------------------
C
C      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C      PARAMETER (MMAX   = 721, NMAX   = 541)
C      PARAMETER (IUNT05 = 5)
C
C      DIMENSION AMP(MMAX,NMAX), PHS(MMAX,NMAX)
C
C      CALL RDMAP(IUNT05, AMP   , PHS   , DX    , DY    ,
C     +           MMAX  , NMAX  , XMIN  , YMAX  , MEND  ,
C     +           NEND  , UNDEF                          )
C
C      DO N = 1,NEND
C         Y = YMAX - DBLE(N-1)*DY
C         DO M = 1,MEND
C            X = XMIN + DBLE(M-1)*DX
C
C            IF (AMP(M,N).LT.UNDEF) THEN
CCC               WRITE(6,'(4F9.4)')X,Y,AMP(M,N),PHS(M,N)
C            ENDIF
C
C         ENDDO
C      ENDDO
C
C 19   STOP
C      END
C
C$DEF_NAOTIDE_GRID
C------------------------------------------------------------
      SUBROUTINE DEF_NAO_GRID( IU, DX, DY , X0, X1, Y0, Y1,
     *                     GLOBX, MEND, NEND, IRET )
C----------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(OUT):: X0, X1, Y0, Y1, DX, DY
      INTEGER,          INTENT(IN ):: IU
      INTEGER,          INTENT(OUT):: GLOBX, MEND, NEND, IRET
C------------------------------------------------------------
C    DX AND DY ARE ABS VALUES IN NAO DATA.
C    DATA ORDER IN FILES AND IN RETURN ARRAYS IS FROM NORTH (YMAX)
C    TO SOUTH (YMIN).
C    I AM COUNTING FOR IT RETURNING X0, X1, Y0, Y1, DX AND DY
C    JUST FOR THIS DATA SET.
C    GLOBX IS ZERO FOR REGIONAL GRIDS AND 1 OR 2 FOR GLOBAL,
C    WITH GLOBX = 2 IF FIRS AND LAST LONG IN FILE ARE SAME (REPEARTED)
C
C    S.VARLAMOV, JANUARY 2005
C
C------------------------------------------------------------
C
C    LOCALS::
C
      CHARACTER        FMT*6, WAVE*3, NAME*50, DATE*50
      INTEGER          IDEFF
      DOUBLE PRECISION AUNIT, PUNIT, D_X
C
C     INVERSE ORDER OF Y BOUNDS IN CALL RDHEAD !!!
C
      CALL RDHEAD( IU, NAME, WAVE, DATE, X0, X1, Y1, Y0, DX, DY,
     *             MEND, NEND, IDEFF, FMT, AUNIT, PUNIT, IRET)

      IF(IRET .GT. 0)RETURN

      DY = -DY
C
C     ASSUME NON-ZERO INCREMENT
C
      D_X = X1 - X0
      IF(DX > 0D0)THEN
        IF(D_X < 0.01*DX ) X1 = X1 + 360D0
      ELSE
        IF(D_X > 0.01*DX ) X1 = X1 - 360D0
      ENDIF

      D_X = DABS(X1 - X0)
      GLOBX = 1
      IF(D_X > 360D0 - 0.01*DABS(DX))THEN
        GLOBX = 2
      ELSE IF(D_X + DABS(DX)*1.01 < 360D0) THEN
        GLOBX = 0
      ENDIF
C
C     RETURN BACK TO THE START OF FILE::
C
      REWIND(IU)

      END SUBROUTINE DEF_NAO_GRID
C
C$RDCMP
C----------------------------------------------------------------
      SUBROUTINE RDCMP(IU    , AMP   , PHS   , MEND  , NEND  ,
     +                 MMAX  , NMAX  , FMT   , AUNIT , PUNIT , IRET )
C----------------------------------------------------------------
C
C     VSM: IT READS AND SAVE DATA IN THE SAME ORDER
C     AS IT IS RECORDED IN FILES
C
      IMPLICIT                        NONE
      INTEGER, PARAMETER::            KC = 10
C
      CHARACTER,        INTENT(IN ):: FMT*6
      INTEGER,          INTENT(IN ):: IU, MEND, NEND, MMAX, NMAX
      INTEGER,          INTENT(OUT):: IRET
      DOUBLE PRECISION, INTENT(IN ):: AUNIT, PUNIT
      DOUBLE PRECISION, INTENT(OUT):: AMP(MMAX,NMAX), PHS(MMAX,NMAX)
C
C     LOCALS::
C
      INTEGER, ALLOCATABLE::  IAMP(:), IPHS(:)
      INTEGER::               KEND, KREM, N, K, M1, M2, M, I, IT

      ALLOCATE( IAMP(MEND), IPHS(MEND), STAT = N)
      IF( N .NE. 0)THEN
        IRET = 1
        GOTO 100
      ENDIF

      AMP = 0D0
      PHS = 0D0
C
C -----< READING LOOP >-----
C
      KEND = MEND/KC
      KREM = MOD(MEND,KC)
      DO N = 1,NEND
         DO K = 1,KEND
            M1 = (K-1)*KC + 1
            M2 = M1 + KC - 1
            READ(IU, FMT, IOSTAT = IRET) (IAMP(M),M = M1, M2)
            IF(IRET .GT. 0)GOTO 100
         ENDDO
         IF (KREM.NE.0) THEN
            M1 = KEND*KC + 1
            M2 = KEND*KC + KREM
            READ(IU, FMT, IOSTAT = IRET) (IAMP(M),M=M1,M2)
            IF(IRET .GT. 0)GOTO 100
         ENDIF
         DO K = 1,KEND
            M1 = (K-1)*KC + 1
            M2 = M1 + KC - 1
            READ(IU, FMT, IOSTAT = IRET) (IPHS(M),M=M1,M2)
            IF(IRET .GT. 0)GOTO 100
         ENDDO
         IF (KREM.NE.0) THEN
            M1 = KEND*KC + 1
            M2 = KEND*KC + KREM
            READ(IU, FMT, IOSTAT = IRET) (IPHS(M),M=M1,M2)
            IF(IRET .GT. 0)GOTO 100
         ENDIF
         DO M = 1,MEND
            AMP(M,N) = DBLE(IAMP(M))*AUNIT  ! IN CENTIMETERS
            PHS(M,N) = DBLE(IPHS(M))*PUNIT  ! IN DEGREES
         ENDDO
      ENDDO

      IRET = 0

100   CONTINUE
      IT = 0
      I = 0
      IF(ALLOCATED( IAMP))DEALLOCATE( IAMP, STAT = I)
      IF(I .NE. 0)THEN
        IT=IT+1
        I = 0
      ENDIF
      IF(ALLOCATED( IPHS))DEALLOCATE( IPHS, STAT = I)
      IF(I .NE. 0)THEN
        IT=IT+2
        I = 0
      ENDIF
      IF(IT .NE. 0)THEN
C        WRITE (*,FMT=*)
C     *  ' *** DEALLOCATION ERROR IN RDCMP:',IT
        IRET = IT
      ENDIF

      END SUBROUTINE RDCMP
C
C$RDHEAD
C----------------------------------------------------------------
      SUBROUTINE RDHEAD(IU    , NAME  , WAVE  , DATE  , XMIN  ,
     +                  XMAX  , YMIN  , YMAX  , DX    , DY    ,
     +                  MEND  , NEND  , IDEFF , FMT   , AUNIT ,
     +                  PUNIT , IRET                           )
C----------------------------------------------------------------
      IMPLICIT                        NONE
      CHARACTER,        INTENT(OUT):: FMT*6, WAVE*3, NAME*50, DATE*50
      DOUBLE PRECISION, INTENT(OUT):: XMIN, XMAX, YMIN, YMAX, DX, DY,
     *                                AUNIT, PUNIT
      INTEGER,          INTENT(IN ):: IU
      INTEGER,          INTENT(OUT):: MEND, NEND, IDEFF, IRET
C
C    LOCALS::
C
      INTEGER::          IDX, IDY
C
      READ(IU,'(13X,A50)', IOSTAT = IRET) NAME
      READ(IU,'(13X,A3 )', IOSTAT = IRET) WAVE
      READ(IU,'(13X,F5.3,19X,F4.2)', IOSTAT = IRET) AUNIT,PUNIT
      READ(IU,'(13X,A50)', IOSTAT = IRET) DATE
      READ(IU,'(7X,F7.2,3(9X,F7.2))', IOSTAT = IRET) XMIN,XMAX,YMIN,YMAX
      READ(IU,'(12X,I2,14X,I2,2(9X,I7))',IOSTAT = IRET)IDX,IDY,MEND,NEND
      READ(IU,'(16X,I6,11X,A6)', IOSTAT = IRET) IDEFF, FMT
      IF(IRET .GT. 0)RETURN
      IF (IDX.EQ.50) THEN
         DX = 0.5D0
         DY = 0.5D0
      ELSE
         DX = 1.D0/DBLE(IDX)
         DY = 1.D0/DBLE(IDY)
      ENDIF

      IRET = 0
      END SUBROUTINE RDHEAD

C$RDMAP
C------------------------------------------------------------
      SUBROUTINE RD_NAO_MAP(IU, AMP, PHS, MMAX, NMAX, UNDEF, IRET)
C------------------------------------------------------------
C
C     REMOVED RETURN VALUES THAT HAVE TO BE DEFINED MORE ACCURATELY
C     WITH THE CALL TO DO DEF_NAOTIDE_GRID
C
C------------------------------------------------------------
      IMPLICIT                        NONE
      INTEGER,          INTENT(IN ):: IU, MMAX, NMAX
      INTEGER,          INTENT(OUT):: IRET
      DOUBLE PRECISION, INTENT(OUT):: AMP(MMAX, NMAX), PHS(MMAX, NMAX),
     *                                UNDEF
C
C    LOCALS::
C
      CHARACTER        FMT*6, WAVE*3, NAME*50,DATE*50
      INTEGER          IDEFF, MEND, NEND
      DOUBLE PRECISION XMAX, YMIN, AUNIT, PUNIT, XMIN, YMAX, DX, DY

      CALL RDHEAD(IU    , NAME  , WAVE  , DATE  , XMIN  ,
     +            XMAX  , YMIN  , YMAX  , DX    , DY    ,
     +            MEND  , NEND  , IDEFF , FMT   , AUNIT ,
     +            PUNIT , IRET                           )

      IF(IRET .NE. 0)RETURN
C
C     CHECK THE CURRENT FILE GRID SIZE WITH ALLOCATED::
C
      IF(NEND > NMAX .OR. MEND > MMAX)THEN
        IRET = -1
        RETURN
      ENDIF

CCCCCC
C      WRITE(6,6002) NAME
C 6002 FORMAT('   MODEL NAME   = ',A50)
C      WRITE(6,6003) DATE
C 6003 FORMAT('   CREATED DATE = ',A50)
C      WRITE(6,101) XMIN, XMAX, YMIN, YMAX
C      WRITE(6,102) DX, DY, MEND, NEND
C      WRITE(6,103) AUNIT, PUNIT
C 101  FORMAT(4X,'XMIN = ',F6.2,', XMAX = ',F6.2,', YMIN = ',F6.2,
C     +          ', YMAX = ',F6.2)
C 102  FORMAT(4X,'DX   = ',F6.3,', DY =   ',F6.3,', MEND = ',I4,
C     +          '  , NEND = ',I4)
C 103  FORMAT(4X,'AMPLITUDE UNIT = ',F6.3,', PHASE UNIT =  ',F6.3)
CCCCCC

      UNDEF = DBLE(IDEFF)*AUNIT

      CALL RDCMP(IU    , AMP   , PHS   , MEND  , NEND  ,
     +           MMAX  , NMAX  , FMT   , AUNIT , PUNIT , IRET )

      END SUBROUTINE RD_NAO_MAP
C
C------------------------------------------------------------
      SUBROUTINE RD_NAO_UV(IU, AMP, PHS, MASK, MMAX, NMAX,
     *  X0, Y0, DX, DY, IRET)
C------------------------------------------------------------
      IMPLICIT                        NONE
      INTEGER,          INTENT(IN ):: IU, MMAX, NMAX
      INTEGER,          INTENT(OUT):: IRET
      DOUBLE PRECISION, INTENT(OUT):: AMP(MMAX, NMAX, 2),
     *                                PHS(MMAX, NMAX, 2)
c
C     NEXT ARE IN DEGREES, C-GRID'S SSH-GRID, NOT U EITHER V GRIDS
c     as it is written in the NAO u&v output
c
      DOUBLE PRECISION, INTENT(IN)::  X0, Y0, DX, DY
      REAL, INTENT(OUT)::             MASK(MMAX, NMAX, 2)
C
C    LOCALS::
C
      DOUBLE PRECISION XX,YY,UA,UP,VA,VP
      INTEGER          IX, JY
      DOUBLE PRECISION,parameter:: EPS=1D-6

      MASK = 0.
      AMP  = 0D0
      PHS  = 0D0
      IRET = 0

c      write(*,*)"RD_UV x0 y0 dx dy ",X0, Y0, DX, DY
C
C     READ FROM NAO MATSUMOTO TIDAL CURRENTS ORIGINAL FILES:
C     AMP IN CM/SEC, PHASE IN DEGREES (GREENWICH)
C
      DO WHILE(.TRUE.)
        READ(UNIT=IU, FMT=731, END=730, ERR=1)XX,YY,UA,UP,VA,VP
        IX=NINT((XX-X0)/DX)+1
        JY=NINT((YY-Y0)/DY)+1
        IF(IX.GE.1.AND.IX.LE.MMAX.AND.
     *     JY.GE.1.AND.JY.LE.NMAX)THEN
          AMP(IX,JY,1)=UA
          AMP(IX,JY,2)=VA
          PHS(IX,JY,1)=UP
          PHS(IX,JY,2)=VP
          IF( ABS(UA) > EPS .or. ABS(UP) > EPS ) MASK(IX,JY,1)=1.
          IF( ABS(VA) > EPS .or. ABS(VP) > EPS ) MASK(IX,JY,2)=1.
        ELSE
          IRET = 2
          RETURN
        ENDIF
      ENDDO
731   FORMAT(2F8.3,4D11.4)
730   CONTINUE
      RETURN
1     IRET = 1
      END SUBROUTINE RD_NAO_UV
C
C------------------------------------------------------------
      SUBROUTINE RD_NAO_MDLDPTH(IU, MDLDPTH, MMAX, NMAX,
     *  X0, Y0, DX, DY, IRET)
C------------------------------------------------------------
C     Y0 IN NAO IS NORTH, DY IS NEGATIVE, DATA IN ARRAYS GO
C     FROM NORT TO SOUTH... DO SAME FOR THE DEPTH...
C     READ DEPTH, REPLACE SMALL BY USED MIN AND CONVERT TO CM
C------------------------------------------------------------
      IMPLICIT                        NONE
      INTEGER,          INTENT(IN ):: IU, MMAX, NMAX
      INTEGER,          INTENT(OUT):: IRET
      DOUBLE PRECISION, INTENT(OUT):: MDLDPTH(MMAX, NMAX)
c
C     NEXT ARE IN DEGREES, C-GRID'S SSH-GRID
c
      DOUBLE PRECISION, INTENT(IN)::  X0, Y0, DX, DY
C
C    LOCALS::
C
      DOUBLE PRECISION LON0,LAT0,DXYR,MDLDPTH_MIN,YY
      INTEGER          I,J,JY,NX,NY
      DOUBLE PRECISION,parameter:: EPS=1D-6
      CHARACTER(256)   TMP1,TMP2

      MDLDPTH = 0D0
      IRET = 0

c      write(*,*)"RD_DEPTH x0 y0 dx dy ",X0, Y0, DX, DY
C
C     READ FROM NAO MATSUMOTO MODEL ORIGINAL DEPTH FILE:
C     MDLDPTH IN M, POSITIVE FOR SEA DEPTH, MIN USED IS 20M
C     (REPLACE SMALL POSITIVE); CHANGE TO CM
C
      READ(UNIT=IU,FMT=*,IOSTAT=IRET)LON0,LAT0,DXYR,NX,NY,MDLDPTH_MIN

      IF(IRET /= 0)THEN
        IRET = 1
        RETURN
      ELSE IF( DABS(LON0-X0) > EPS .OR.
     *    NX > MMAX .OR. NY > NMAX .OR.
     *    DABS(1D0/DXYR-DX) > EPS .OR.
     *    DABS(1D0/DXYR-DABS(DY)) > EPS )THEN
        IRET = 2
        RETURN
      ENDIF

      DO J=1,NY
        READ(UNIT=IU, FMT=*, IOSTAT=IRET)TMP1,TMP2,YY

        JY=NINT((YY-Y0)/DY)+1

c      write(*,*)"RD_DEPTH J Jy ",j, jY

        IF(JY .GE. 1 .AND. JY .LE. NY)THEN
          READ(UNIT=IU, FMT='(20f5.0)', IOSTAT=IRET)
     *     (MDLDPTH(I,JY),I=1,NX)

          IF(IRET /= 0)THEN
            IRET = 3
            RETURN
          ENDIF

          DO I=1,NX
            IF(MDLDPTH(I,JY) <= 0D0)THEN
              MDLDPTH(I,JY)=0D0
            ELSE
              IF(MDLDPTH(I,JY) < MDLDPTH_MIN)THEN
                MDLDPTH(I,JY)=MDLDPTH_MIN*100D0
              ELSE
                MDLDPTH(I,JY)=MDLDPTH(I,JY)*100D0
              ENDIF
            ENDIF
          ENDDO
        ELSE
          IRET = 4
          RETURN
        ENDIF
      ENDDO

      END SUBROUTINE RD_NAO_MDLDPTH
C
C----------------------< END OF NAO INPUT >----------------------

C----LEGOS 2002 INPUT -------------------------------------------
C
C    LOW LEVEL SUBROUTINES TO READ LEGOS-2002 0.25 DEGREE TIDAL DATA
C    METEO-FRANCE DPREVI-MAR LEGOS DATA IN ASCII FORMAT ARE WRITTEN
C    FROM SOUTH POLE TO NORTH POLE (BOTH INCLUDED AND BOTH UBDEFINED),
C    ONE LINE FOR EACH LONGITUDE, TOTAL 721 PAIRS OF VALUES (AMP-PHS).
C    LONGITUDE IS FROM ZERO TO 360, FIRST AND LAST LINES SAME
C    (0 AND 360 DUBLICATED), TOTAL 1441 LINES OF DATA
C
C    LAST MODIFICATION - MAY, 2005
C
C    SERGEY M. VARLAMOV
C
C------------------------------------------------------------
      SUBROUTINE DEF_LEGOS_GRID( DX, DY, X0, X1, Y0, Y1,
     *                     GLOBX, MEND, NEND)
C----------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(OUT):: X0, X1, Y0, Y1, DX, DY
      INTEGER,          INTENT(OUT):: GLOBX, MEND, NEND
C------------------------------------------------------------
C    GLOBX IS ZERO FOR REGIONAL GRIDS AND 1 OR 2 FOR GLOBAL,
C    WITH GLOBX = 2 IF FIRS AND LAST LONG IN FILE ARE SAME (REPEATED)
C
C    S.VARLAMOV, MAY, 2005
C
C------------------------------------------------------------
C
C    LOCALS::
C
      DOUBLE PRECISION D_X

      DX = 0.25
      DY = 0.25
      X0 = 0.0
      X1 = 360.0
      Y0 = -90.0
      Y1 = 90.0
      MEND = 1441
      NEND = 721
C
C     ASSUME NON-ZERO INCREMENT
C
      D_X = X1 - X0
      IF(DX > 0D0)THEN
        IF(D_X < 0.01*DX ) X1 = X1 + 360D0
      ELSE
        IF(D_X > 0.01*DX ) X1 = X1 - 360D0
      ENDIF

      D_X = DABS(X1 - X0)
      GLOBX = 1
      IF(D_X > 360D0 - 0.01*DABS(DX))THEN
        GLOBX = 2
      ELSE IF(D_X + DABS(DX)*1.01 < 360D0) THEN
        GLOBX = 0
      ENDIF

      END SUBROUTINE DEF_LEGOS_GRID

C----------------------------------------------------------------
      SUBROUTINE RDCMP_LEGOS(IU    , AMP   , PHS   , MEND  , NEND  ,
     +                       MMAX  , NMAX  , IRET )
C----------------------------------------------------------------
      IMPLICIT                        NONE
      INTEGER,          INTENT(IN ):: IU, MEND, NEND, MMAX, NMAX
      INTEGER,          INTENT(OUT):: IRET
      DOUBLE PRECISION, INTENT(OUT):: AMP(MMAX,NMAX), PHS(MMAX,NMAX)
C
C     LOCALS::
C
      REAL, ALLOCATABLE::  RAMP_PHS(:)
      INTEGER::            N, M

      ALLOCATE( RAMP_PHS(NEND*2), STAT = N)
      IF( N .NE. 0)THEN
        IRET = 1
        GOTO 100
      ENDIF
C
C -----< READING LOOP >-----
C
      DO M = 1,MEND
        READ(IU, *, IOSTAT = IRET) RAMP_PHS
        IF(IRET > 0)GOTO 100
        DO N = 1,NEND
          AMP(M,N) = DBLE(RAMP_PHS(N*2-1))  ! IN CENTIMETERS
          PHS(M,N) = DBLE(RAMP_PHS(N*2  ))  ! IN DEGREES
        ENDDO
      ENDDO

      IRET = 0

100   CONTINUE
      IF(ALLOCATED(RAMP_PHS))DEALLOCATE(RAMP_PHS, STAT = N)
      IF(N .NE. 0)THEN
C        WRITE (*,FMT=*)
C     *  ' *** DEALLOCATION ERROR IN RDCMP_LEGOS:',N
        IRET = 2
      ENDIF

      END SUBROUTINE RDCMP_LEGOS

C------------------------------------------------------------
      SUBROUTINE RD_LEGOS_MAP(IU, AMP, PHS, MMAX, NMAX, UNDEF, IRET)
C------------------------------------------------------------
C------------------------------------------------------------
      IMPLICIT                        NONE
      INTEGER,          INTENT(IN ):: IU, MMAX, NMAX
      INTEGER,          INTENT(OUT):: IRET
      DOUBLE PRECISION, INTENT(OUT):: AMP(MMAX, NMAX), PHS(MMAX, NMAX),
     *                                UNDEF
C
C    LOCALS::
C
      INTEGER          GLOBX, MEND, NEND
      DOUBLE PRECISION X0, Y0, X1, Y1, DX, DY

      UNDEF = 999.
      CALL DEF_LEGOS_GRID( DX, DY, X0, X1, Y0, Y1,
     *                     GLOBX, MEND, NEND)

      CALL RDCMP_LEGOS(IU    , AMP   , PHS   , MEND  , NEND  ,
     +                 MMAX  , NMAX  , IRET )

      END SUBROUTINE RD_LEGOS_MAP
C
C ----------------------< END OF LEGOS 2002 INPUT >----------------------

C--- LEGOS FES2004 INPUT-------------------------------------------------
C
C DEFINE THE FES2004 MODEL GRID FROM FILE HEADER,
C READ THE SEA LEVEL TIDAL HARMONICS AMPLITUDE AND PHAZE
C
C    DECEMBER, 2008
C
C    SERGEY M. VARLAMOV
C
C$MAIN
C--------------------------------------------------------------------------

C------------------------------------------------------------
      SUBROUTINE DEF_FES2004_GRID( IU, DX, DY , X0, X1, Y0, Y1,
     *                     GLOBX, MEND, NEND, IRET )
C----------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(OUT):: X0, X1, Y0, Y1, DX, DY
      INTEGER,          INTENT(IN ):: IU
      INTEGER,          INTENT(OUT):: GLOBX, MEND, NEND, IRET
C------------------------------------------------------------
C    DX AND DY ARE ABS VALUES IN NAO DATA.
C    DATA ORDER IN FILES AND IN RETURN ARRAYS IS FROM NORTH (YMAX)
C    TO SOUTH (YMIN).
C    I AM COUNTING FOR IT RETURNING X0, X1, Y0, Y1, DX AND DY
C    JUST FOR THIS DATA SET.
C    GLOBX IS ZERO FOR REGIONAL GRIDS AND 1 OR 2 FOR GLOBAL,
C    WITH GLOBX = 2 IF FIRS AND LAST LONG IN FILE ARE SAME (REPEARTED)
C
C    S.VARLAMOV, JANUARY 2005
C
C------------------------------------------------------------
C
C    LOCALS::
C
      DOUBLE PRECISION AUNDEF, PUNDEF, D_X
C
C     INVERSE ORDER OF Y BOUNDS IN CALL RDHEAD !!!
C
!      write(*,*)"Define fes2004 grid", iu

      CALL RDHEAD_FES2004( IU, X0, X1, Y0, Y1, DX, DY,
     *             MEND, NEND, AUNDEF, PUNDEF, IRET)

      IF(IRET .GT. 0)RETURN
C
C     ASSUME NON-ZERO INCREMENT
C
      D_X = X1 - X0
      IF(DX > 0D0)THEN
        IF(D_X < 0.01*DX ) X1 = X1 + 360D0
      ELSE
        IF(D_X > 0.01*DX ) X1 = X1 - 360D0
      ENDIF

      D_X = DABS(X1 - X0)
      GLOBX = 1
      IF(D_X > 360D0 - 0.01*DABS(DX))THEN
        GLOBX = 2
      ELSE IF(D_X + DABS(DX)*1.01 < 360D0) THEN
        GLOBX = 0
      ENDIF
C
C     RETURN BACK TO THE START OF FILE::
C
      REWIND(IU)

      END SUBROUTINE DEF_FES2004_GRID
C
C$RDCMP
C----------------------------------------------------------------
      SUBROUTINE RDCMP_FES2004(IU    , AMP   , PHS   , MEND  , NEND  ,
     +                 MMAX  , NMAX  , IRET )
C----------------------------------------------------------------
C
C     VSM: IT READS AND SAVE DATA IN THE SAME ORDER
C     AS IT IS RECORDED IN FILES
C
      IMPLICIT                        NONE
      INTEGER, PARAMETER::            KC = 30
C
      INTEGER,          INTENT(IN ):: IU, MEND, NEND, MMAX, NMAX
      INTEGER,          INTENT(OUT):: IRET
      DOUBLE PRECISION, INTENT(OUT):: AMP(MMAX,NMAX), PHS(MMAX,NMAX)
C
C     LOCALS::
C
      REAL, ALLOCATABLE::     RAMP(:), RPHS(:)
      INTEGER::               KEND, KREM, N, K, M1, M2, M, I, IT

      ALLOCATE( RAMP(MEND), RPHS(MEND), STAT = N)
      IF( N .NE. 0)THEN
        IRET = 1
        GOTO 100
      ENDIF

      AMP = 0D0
      PHS = 0D0
C
C -----< READING LOOP >-----
C
      REWIND(IU)
      DO N = 1,5
        READ(IU,*)
      ENDDO

      KEND = MEND/KC
      KREM = MOD(MEND,KC)
      DO N = 1,NEND
         DO K = 1,KEND
            M1 = (K-1)*KC + 1
            M2 = M1 + KC - 1
!      write(*,*)"input 1: ",n,k,kc,m1,m2,MEND
            READ(IU, *) (RAMP(M),M = M1,M2)
            READ(IU, *) (RPHS(M),M = M1,M2)
         ENDDO
         IF (KREM.NE.0) THEN
            M1 = KEND*KC + 1
            M2 = KEND*KC + KREM
!      write(*,*)"input 2: ",n,k,kc,m1,m2,MEND
            READ(IU, *) (RAMP(M),M=M1,M2)
            READ(IU, *) (RPHS(M),M=M1,M2)
         ENDIF
         DO M = 1,MEND
            AMP(M,N) = DBLE(RAMP(M))  ! IN CENTIMETERS
            PHS(M,N) = DBLE(RPHS(M))  ! IN DEGREES
         ENDDO
      ENDDO
!      write(*,*)"input done"
      IRET = 0

100   CONTINUE
      IT = 0
      I = 0
      IF(ALLOCATED( RAMP))DEALLOCATE( RAMP, STAT = I)
      IF(I .NE. 0)THEN
        IT=IT+1
        I = 0
      ENDIF
      IF(ALLOCATED( RPHS))DEALLOCATE( RPHS, STAT = I)
      IF(I .NE. 0)THEN
        IT=IT+2
        I = 0
      ENDIF
      IF(IT .NE. 0)THEN
C        WRITE (*,FMT=*)
C     *  ' *** DEALLOCATION ERROR IN RDCMP:',IT
        IRET = IT
      ENDIF

      END SUBROUTINE RDCMP_FES2004
C
C$RDHEAD
C----------------------------------------------------------------
      SUBROUTINE RDHEAD_FES2004(IU    , XMIN  , XMAX  ,
     * YMIN  , YMAX  ,  DX    , DY    , MEND  , NEND  ,
     * AUNDEF , PUNDEF, IRET)
C----------------------------------------------------------------
      IMPLICIT                        NONE
      DOUBLE PRECISION, INTENT(OUT):: XMIN, XMAX, YMIN, YMAX, DX, DY,
     *                                AUNDEF, PUNDEF
      INTEGER,          INTENT(IN ):: IU
      INTEGER,          INTENT(OUT):: MEND, NEND, IRET
!      real t1
C
      IRET = 0
!      write(*,*)"read header: ",iu
      REWIND(IU)
      READ(IU,*) XMIN,XMAX
!      write(*,*)"read header: ",XMIN, XMAX
      READ(IU,*) YMIN,YMAX
      READ(IU,*) DX,DY
      READ(IU,*) MEND,NEND
      READ(IU,*) AUNDEF, PUNDEF

!      write(*,*)"read header: ",XMIN, XMAX, YMIN, YMAX, DX, DY,
!     *AUNDEF, PUNDEF,MEND, NEND
c
c     start to read next record; else gfortran cause SEG error on data import
c
!      READ(IU,*) t1
!      write(*,*)"read values: ",t1
!      BACKSPACE(IU)

      END SUBROUTINE RDHEAD_FES2004

C$RDMAP
C------------------------------------------------------------
      SUBROUTINE RD_FES2004_MAP(IU, AMP, PHS, MMAX, NMAX, UNDEF, IRET)
C------------------------------------------------------------
C
C     REMOVED RETURN VALUES THAT HAVE TO BE DEFINED MORE ACCURATELY
C     WITH THE CALL TO DO DEF_NAOTIDE_GRID
C
C------------------------------------------------------------
      IMPLICIT                        NONE
      INTEGER,          INTENT(IN ):: IU, MMAX, NMAX
      INTEGER,          INTENT(OUT):: IRET
      DOUBLE PRECISION, INTENT(OUT):: AMP(MMAX, NMAX), PHS(MMAX, NMAX),
     *                                UNDEF
C
C    LOCALS::
C
      INTEGER          MEND, NEND
      DOUBLE PRECISION XMAX, YMIN, AUNDEF, PUNDEF, XMIN, YMAX, DX, DY

      CALL RDHEAD_FES2004(IU    , XMIN  ,
     +            XMAX  , YMIN  , YMAX  , DX    , DY    ,
     +            MEND  , NEND  , AUNDEF ,
     +            PUNDEF , IRET                           )

      IF(IRET .NE. 0)RETURN
C
C     CHECK THE CURRENT FILE GRID SIZE WITH ALLOCATED::
C
      IF(NEND > NMAX .OR. MEND > MMAX)THEN
        IRET = -1
        RETURN
      ENDIF

      UNDEF = AUNDEF

      CALL RDCMP_FES2004(IU    , AMP   , PHS   , MEND  , NEND  ,
     +           MMAX  , NMAX  , IRET )

      END SUBROUTINE RD_FES2004_MAP

C----------------------------------------------------------------
      SUBROUTINE DEF_OTIS_GRID( IU, DX, DY , X0, X1, Y0, Y1,
     *                     GLOBX, MEND, NEND, IRET )
C----------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(OUT):: X0, X1, Y0, Y1, DX, DY
      INTEGER,          INTENT(IN ):: IU
      INTEGER,          INTENT(OUT):: GLOBX, MEND, NEND, IRET
C------------------------------------------------------------
! From the OTIS readme file and code:
!
!The lats and lons (in degrees) for individual grid nodes are:
!
!For h-nodes   lat(i,j) =  theta_lim(1) + (j-.5)*dy
!              lon(i,j) =  phi_lim(1) + (i-.5)*dx
!
!For U-nodes:  lat(i,j) =  theta_lim(1) + (j-.5)*dy
!              lon(i,j) =  phi_lim(1) + (i- 1)*dx
!
!For V-nodes:  lat(i,j) =  theta_lim(1) + (j- 1)*dy
!              lon(i,j) =  phi_lim(1) + (i-.5)*dx
!
!Here dx,dy is resolution in longitude and latitude correspondingly.
!        dtheta = (th_lim(2)-th_lim(1))/m
!        dphi = (ph_lim(2)-ph_lim(1))/n
!------------------------------------------------------------
C    GLOBX IS ZERO FOR REGIONAL GRIDS AND 1 OR 2 FOR GLOBAL,
C    WITH GLOBX = 2 IF FIRST AND LAST LONG IN FILE ARE SAME (REPEARTED)
C
C    S.VARLAMOV, JUNE 2010
C------------------------------------------------------------
C
C    LOCALS::
C
      DOUBLE PRECISION D_X
      integer nco
      real theta_lim(2),phi_lim(2)

      IRET = 0

      read(IU) MEND, NEND, nco, theta_lim, phi_lim
!     start OTIS SPECIFIC, see above
      DX = (phi_lim(2)-phi_lim(1))/MEND
      DY = (theta_lim(2)-theta_lim(1))/NEND

      X0 = phi_lim(1)+0.5*DX
      X1 = phi_lim(2)-0.5*DX
      Y0 = theta_lim(1)+0.5*DY
      Y1 = theta_lim(2)-0.5*DY
!     end OTIS SPECIFIC
C
C     ASSUME NON-ZERO INCREMENT
C
      D_X = X1 - X0
      IF(DX > 0D0)THEN
        IF(D_X < 0.01*DX ) X1 = X1 + 360D0
      ELSE
        IF(D_X > 0.01*DX ) X1 = X1 - 360D0
      ENDIF

      D_X = DABS(X1 - X0)
      GLOBX = 1
      IF(D_X > 360D0 - 0.01*DABS(DX))THEN
        GLOBX = 2
      ELSE IF(D_X + DABS(DX)*1.01 < 360D0) THEN
        GLOBX = 0
      ENDIF
C
C     RETURN BACK TO THE START OF FILE::
C
      REWIND(IU)

      END SUBROUTINE DEF_OTIS_GRID

C----------------------------------------------------------------
      SUBROUTINE RD_H_OTIS(IU    , AMP   , PHS   , MEND  , NEND  ,
     +                     HARMONIC_NAME, UNDEF, IRET )
C----------------------------------------------------------------
C
C     VSM: IT READS AND SAVE DATA IN THE SAME ORDER
C     AS IT IS RECORDED IN FILES
C
      USE COM_SERV, ONLY: D_IN_RAD
      IMPLICIT                        NONE
      INTEGER,          INTENT(IN ):: IU, MEND, NEND
      INTEGER,          INTENT(OUT):: IRET
      DOUBLE PRECISION, INTENT(OUT):: UNDEF,
     *                                AMP(MEND,NEND),PHS(MEND,NEND)
      CHARACTER*4,      INTENT(IN ):: HARMONIC_NAME
      CHARACTER*4       HARMONIC_NAMES(NHARM_ALL)
C
C     LOCALS::
C
      INTEGER::               N, IC, K, M, I, IT
      complex, allocatable::  ho(:,:)
      integer no,mo,nco
      real theta_lim(2),phi_lim(2)

      AMP = 0D0
      PHS = 0D0
      UNDEF = 990d0
      IRET = 0

      ALLOCATE( ho(MEND,NEND), STAT = N)

      IF( N .NE. 0)THEN
        IRET = 1
        GOTO 105
      ENDIF
C
C -----< READ >-----
C
      REWIND(IU)
!     read header
      READ(IU) no,mo,nco,theta_lim,phi_lim,
     *              HARMONIC_NAMES(1:nco)
      if(no /= MEND .or. mo /= NEND .or. nco > NHARM_ALL)then
        IRET = 2
        goto 105
      endif
#ifdef _TDTEST_
      write(*,*) "read ",no,mo,nco,theta_lim,phi_lim,
     *              HARMONIC_NAMES(1:nco)
      write(*,*) "requested ",HARMONIC_NAME
#endif

      do ic=1,nco
#ifdef _TDTEST_
      write(*,*) "read ",ic
#endif
       read(IU)ho
       if(HARMONIC_NAMES(IC)==HARMONIC_NAME)exit
      enddo
      if(ic > nco)then
        IRET = 3
        goto 105
      endif

      DO N = 1,NEND
         DO M = 1,MEND
            AMP(M,N) = ABS(ho(M,N))*100D0  ! IN CENTIMETERS
            if(AMP(M,N) > 0d0)then
              PHS(M,N) = 
     * D_IN_RAD*ATAN2(-AIMAG(ho(M,N)),REAL(ho(M,N)))  ! IN DEGREES
            else
              AMP(M,N) = 999D0
              PHS(M,N) = 999D0
            endif
         ENDDO
      ENDDO
!      write(*,*)"input done"

105   CONTINUE
      IT = 0
      I = 0
      IF(ALLOCATED( ho))DEALLOCATE( ho, STAT = I)
      IF(I .NE. 0)THEN
        IT=IT+1
        I = 0
      ENDIF
      IF(IT .NE. 0)THEN
C        WRITE (*,FMT=*)
C     *  ' *** DEALLOCATION ERROR IN RDCMP:',IT
        IRET = IT
      ENDIF

      END SUBROUTINE RD_H_OTIS

C----------------------------------------------------------------
      SUBROUTINE RD_UV_OTIS(IU, AMP, PHS, MASK, MEND, NEND,
     +                     HARMONIC_NAME, IRET )
C----------------------------------------------------------------
C
C     VSM: IT READS AND SAVE DATA IN THE SAME ORDER
C     AS IT IS RECORDED IN FILES
C
      USE COM_SERV, ONLY: D_IN_RAD
      IMPLICIT                        NONE
      INTEGER,          INTENT(IN ):: IU, MEND, NEND
      INTEGER,          INTENT(OUT):: IRET
      DOUBLE PRECISION, INTENT(OUT):: AMP(MEND,NEND,2),PHS(MEND,NEND,2)
      CHARACTER*4,      INTENT(IN ):: HARMONIC_NAME
      CHARACTER*4       HARMONIC_NAMES(NHARM_ALL)
      REAL, INTENT(OUT)::             MASK(MEND,NEND,2)
C
C     LOCALS::
C
      INTEGER::               N, IC, M, I, IT
      complex, allocatable::  uvo(:,:,:)
      integer no,mo,nco
      real theta_lim(2),phi_lim(2)

      IRET = 0
      AMP = 0D0
      PHS = 0D0
      MASK = 0.

      ALLOCATE( uvo(2,MEND,NEND), STAT = N)

      IF( N .NE. 0)THEN
        IRET = 1
        GOTO 106
      ENDIF
C
C -----< READ >-----
C
      REWIND(IU)
!     read header
      READ(IU) no,mo,nco,theta_lim,phi_lim,
     *              HARMONIC_NAMES(1:nco)
      if(no /= MEND .or. mo /= NEND .or. nco > NHARM_ALL)then
        IRET = 2
        goto 106
      endif
      do ic=1,nco
       read(IU)uvo
       if(HARMONIC_NAMES(IC)==HARMONIC_NAME)exit
      enddo
      if(ic > nco)then
        IRET = 3
        goto 106
      endif

      DO N = 1,NEND
         DO M = 1,MEND
            AMP(M,N,1) = ABS(uvo(1,M,N))
            if(AMP(M,N,1) > 0d0)then
              AMP(M,N,1) = AMP(M,N,1)*10000D0  ! IN CENTIMETERS**2/s
              PHS(M,N,1) = D_IN_RAD*
     *          ATAN2(-AIMAG(uvo(1,M,N)),REAL(uvo(1,M,N)))  ! IN DEGREES
              MASK(M,N,1) = 1.
            endif
            AMP(M,N,2) = ABS(uvo(2,M,N))
            if(AMP(M,N,2) > 0d0)then
              AMP(M,N,2) = AMP(M,N,2)*10000D0  ! IN CENTIMETERS**2/s
              PHS(M,N,2) = D_IN_RAD*
     *          ATAN2(-AIMAG(uvo(2,M,N)),REAL(uvo(2,M,N)))  ! IN DEGREES
              MASK(M,N,2) = 1.
            endif
         ENDDO
      ENDDO
!      write(*,*)"input done"
      IRET = 0

106   CONTINUE
      IT = 0
      I = 0
      IF(ALLOCATED( uvo))DEALLOCATE( uvo, STAT = I)
      IF(I .NE. 0)THEN
        IT=IT+1
        I = 0
      ENDIF
      IF(IT .NE. 0)THEN
C        WRITE (*,FMT=*)
C     *  ' *** DEALLOCATION ERROR IN RDCMP:',IT
        IRET = IT
      ENDIF

      END SUBROUTINE RD_UV_OTIS
C
C
C----------------------< END OF FES2004 INPUT >----------------------
      END MODULE

#ifdef _TDTEST_

!#ifdef VFD
!MS$  IF DEFINED (_TDTEST_)
!#endif

C***********************************************************
C     TEST PROGRAM - ONE POINT, LIMITED ERROR CHECKING ETC.
#define USE_TOTAL_GRID

      PROGRAM TEST_TIDE
C***********************************************************
C     (1) INITIALIZATION:
C-----------------------------------------------------------
        USE COM_SERV, ONLY: TPOINT, TDATE, TSEC, TDT, RAD_IN_D
        USE TIDE_SERV
!        USE TOPO_SERV, ONLY: TOPO_D,INIT_TOPO,DEPTH,FREE_TOPOGR
        IMPLICIT NONE
        TYPE(TIDE_D)::  SEA_TIDE
C-----------------------------------------------------------
C     (1.0) estimate fixed geographical parameters for your model
C-----------------------------------------------------------
C         HERE: NUM_TH = SEA_TIDE%NUM_TH, it could be modified
c                                         if some files do not exist
C               NPT  - NUMBER OF POINTS TO FIND TIDAL harmonics,
C               NLON - NUMBER OF POINTS IN EAST-WEST DIRECTION
C                      to estimate tidal forcing,
C               NLAT - NUMBER OF POINTS IN SOUTH - NORTH DIRECTION
C                      to estimate tidal forcing,
        INTEGER,PARAMETER:: NPT = 1,
     *                      NLON = NPT,
     *                      NLAT = NPT,
     *                      NF_IO = 33,
     *                      NF_LOG = 44,
     *                      TIDE_DEBUG_PRINT = 1
c        INTEGER,PARAMETER::  NUM_TH = 1
c        CHARACTER(3)::       NAMES(NUM_TH) =(/"m2 "/)
        INTEGER::   NUM_TH = 0
#ifdef USE_TOTAL_GRID
        INTEGER,PARAMETER:: NUM_REGIONS = 0
#else
        INTEGER,PARAMETER:: NUM_REGIONS = 3  !ok - checked 0,1,3
        REAL:: XMIN_VALUES(NUM_REGIONS) = (/-15., 130., -15. /),
     *         XMAX_VALUES(NUM_REGIONS) = (/2.,   150.,   5. /),
     *         YMIN_VALUES(NUM_REGIONS) = (/42.,  30.,   52. /),
     *         YMAX_VALUES(NUM_REGIONS) = (/52.,  40.,   62. /)   !OK - checked 1,3
#endif
!
!       Saint peter port, France, compare with http://www.shom.fr, "Predict..."
!
!!        TYPE(TPOINT)::   PT = TPOINT(-2.5, 49.45, 0.0)
!        TYPE(TPOINT)::   PT = TPOINT(-7.688, 59.195, 0.0)
        TYPE(TPOINT)::   PT = TPOINT(126.0, 32.0, 0.0)

        REAL,allocatable:: TH_AMP(:,:), TH_PHS(:,:),
     *                     TH_AMP_SIN(:,:), TH_AMP_COS(:,:),
     *                     H_EQ_LAT(:,:),
     *                     FA(:),PHA(:),
     *                     F(:),PH(:),STDPH(:),CTDPH(:)
        INTEGER            ITIDE_TYPE

        INTEGER M,IERR,NR,M_EQ_MAX,JT,iparam,num_param
        REAL    H_TMP,H_EQV,TD_BETA,H_EQVL,FF,PP,TT,TMP,TM

        TYPE (TDATE)::    DT
        DOUBLE PRECISION  GLB_UTC_TIME,UTC_STOP_TIME,TIME_STEP

C-----------------------------------------------------------

        OPEN(NF_LOG, FILE = "TEST_TIDE.LOG")

        DT = TDATE(1983,1,1,0,0,0)
        GLB_UTC_TIME = TSEC(DT)

        DT = TDATE(1983,1,2,0,0,0)
        UTC_STOP_TIME = TSEC(DT)

        TIME_STEP = 3600D0

C      ! OR INVERSE
C
C      DT = TDT(GLB_UTC_TIME)
C

C-----------------------------------------------------------
C     (1.1) define data set for input of tidal harmonics
C-----------------------------------------------------------
        SEA_TIDE%PATH = "/home/vsm/OTIS/OTPS/DATA"
!        SEA_TIDE%PATH = "/home/vsm/FES2004/data"
!        SEA_TIDE%PATH = "/home/vsm/RIAMOM/NAO"
!        SEA_TIDE%PATH = "/home/G10201/vsm/RIAMOM/NAO"
c
c     TIDE_NAO99_GLOBAL = 1        for global NAO99 SSH,
c     TIDE_NAO99_JAPAN_SSH = 2     for regional NAO99 around Japan (SSH only),
c     TIDE_NAO99_JAPAN_SSH_UVH = 3 for regional NAO99 around Japan (SSH and UV),
c     TIDE_LEGOS_2002 = 4          for global LEGOS2002 0.25 degree SSH
c     TIDE_FES2004 = 5             for global FES2004 0.125 degree SSH
!     TIDE_OTIS_GLOBAL             for global OTIS 0.25 SSH and UH,VH
c
!        SEA_TIDE%TIDE_MODEL = TIDE_FES2004
!        SEA_TIDE%TIDE_MODEL = TIDE_OTIS_GLOBAL
        SEA_TIDE%TIDE_MODEL = TIDE_OTIS_PO
!        SEA_TIDE%TIDE_MODEL = TIDE_NAO99_JAPAN_SSH_UVH
c       initialize other members
        SEA_TIDE%NUM_REGIONS = -1
        NULLIFY(SEA_TIDE%GRID)
        NULLIFY(SEA_TIDE%HARMONIC_NAME)
C-----------------------------------------------------------
C     (1.2) define what tidal harmonics to use
C-----------------------------------------------------------

C        CALL ALLOC_TIDE_HARMONICS(SEA_TIDE, NUM_TH, IERR)         !MOD
C        WRITE(NF_LOG,*)" @@@ ALLOC_TIDE_HARMONICS RETURN ", IERR
C        SEA_TIDE%HARMONIC_NAME(1:NUM_TH) = NAMES(1:NUM_TH)
C
C        alternative for all available harmonics (NHARM_ALL)
C
         SEA_TIDE%NUM_TH = NUM_TH
         SEA_TIDE%ENABLE_ADM=1 !estimate if possible some tide components by addmittance

C-----------------------------------------------------------
C     (1.3) define sub-regions for interpolation of tidal constituents;
C           if only astronomical tidal parameters are necessary -
C           use 0 and skip initialization and usage of related subs
C-----------------------------------------------------------
#ifndef USE_TOTAL_GRID
        CALL ALLOC_TIDE_GRID(SEA_TIDE, NUM_REGIONS, IERR)         !MOD
        WRITE(NF_LOG,*)" @@@ ALLOC_TIDE_GRID RETURN ", IERR
     *  ," for ",NUM_REGIONS
        WRITE(*,*)" @@@ ALLOC_TIDE_GRID RETURN ", IERR
     *  ," for ",NUM_REGIONS
        if(ierr .ne. 0)stop

          SEA_TIDE%GRID(1:NUM_REGIONS)%XMIN = XMIN_VALUES(1:NUM_REGIONS)
          SEA_TIDE%GRID(1:NUM_REGIONS)%XMAX = XMAX_VALUES(1:NUM_REGIONS)
          SEA_TIDE%GRID(1:NUM_REGIONS)%YMIN = YMIN_VALUES(1:NUM_REGIONS)
          SEA_TIDE%GRID(1:NUM_REGIONS)%YMAX = YMAX_VALUES(1:NUM_REGIONS)
!#else
C    ! ALL AVAILABLE GRID WILL BE USED, IT IS NOT THE BEST SOLUTION
C    ! FOR GLOBAL DATA...
#endif
C-----------------------------------------------------------
C     (1.4) initialize fixed astronomic parameters
C-----------------------------------------------------------
        write(*,*)"INIT_TIDAL_EQ_PARAM:"
        CALL INIT_TIDAL_EQ_PARAM( SEA_TIDE, IERR, LOG=NF_LOG, IO=NF_IO)         !MOD

        WRITE(NF_LOG,*)" @@@ INIT_TIDAL_EQ_PARAM RETURN ", IERR
        WRITE(*,*)" @@@ INIT_TIDAL_EQ_PARAM RETURN ", IERR,
     *  SEA_TIDE%NUM_TH, SEA_TIDE%NUM_TH_INPUT,
     *  SEA_TIDE%HARMONIC_NAME(1:SEA_TIDE%NUM_TH)
        WRITE(*,*)" @@@ INIT_TIDAL_EQ_PARAM NUM_IN_USER: ",
     *  NUM_IN_USER
        WRITE(*,*)" @@@ INIT_TIDAL_EQ_PARAM NUM_IN_TMOD: ",
     *  NUM_IN_TMOD
        WRITE(*,*)" @@@ INIT_TIDAL_EQ_PARAM NUM_IN_TMOD: ",
     *  NUM_IN_INPUT

        if(ierr .ne. 0)stop
C-----------------------------------------------------------
C     (1.5) read tidal harmonics and
C           estimate interpolation polinoms                      !MOD
C-----------------------------------------------------------
        write(*,*)"INIT_SL_TIDAL_HARMONICS:"
        CALL INIT_SL_TIDAL_HARMONICS                             !MOD
     *     ( SEA_TIDE, IERR, IO=NF_IO, LOG=NF_LOG,
     *       DEBUG = TIDE_DEBUG_PRINT,VOLUME_FLUX=1)

        WRITE(NF_LOG,*)" @@@ INIT_SL_TIDAL_HARMONICS RETURN ", IERR
        WRITE(*,*)" @@@ INIT_SL_TIDAL_HARMONICS RETURN ", IERR
        if(ierr .ne. 0)stop
        NUM_TH = SEA_TIDE%NUM_TH
c
c       allocate necessary working arrays
c
        write(*,*)"allocate:"
        allocate(
     *  TH_AMP(NPT,NUM_TH),
     *  TH_PHS(NPT,NUM_TH),
     *  TH_AMP_SIN(NPT,NUM_TH),
     *  TH_AMP_COS(NPT,NUM_TH),
     *  H_EQ_LAT(NHARM_ALL,NLAT),
     *  FA(NHARM_ALL),
     *  PHA(NHARM_ALL),
     *  F(NUM_TH),
     *  PH(NUM_TH),
     *  STDPH(NUM_TH),
     *  CTDPH(NUM_TH), stat=ierr)
        WRITE(NF_LOG,*)" @@@ allocate RETURN ", IERR
        WRITE(*,*)" @@@ allocate RETURN ", IERR
        if(ierr .ne. 0)stop


        NR = SEA_TIDE%NUM_REGIONS
        WRITE(NF_LOG,*)" @@@ TIDE_GRID NR ", NR
        WRITE(NF_LOG,*)" @@@ TIDE_GRIDS   ",
     *   ( SEA_TIDE%GRID(M)%XMIN,
     *     SEA_TIDE%GRID(M)%XMAX,
     *     SEA_TIDE%GRID(M)%YMIN,
     *     SEA_TIDE%GRID(M)%YMAX, M=1,NR)

C-----------------------------------------------------------
C
C     (1.6.1) define fixed amplitude and phase of tidal harmonics
C             in given points, if CUT_EXTR > 0 - no amplitude values
C             exceeding original gridded values will be allowed
C             in interpolated amplitude, silent.
C-----------------------------------------------------------
        write(*,*)"GET_SL_TIDAL_HARMONICS: "
        if(SEA_TIDE%TIDE_MODEL >= TIDE_OTIS_GLOBAL .or.
!     *       SEA_TIDE%TIDE_MODEL == TIDE_OTIS_PO .or.
     *       SEA_TIDE%TIDE_MODEL == TIDE_NAO99_JAPAN_SSH_UVH)then
          num_param = 3
        else
          num_param = 1
        endif

        DO M = 1, SEA_TIDE%NUM_TH
          do iparam = 1,num_param
          write(*,*)"GET_SL_TIDAL_HARMONICS for param: ",iparam,
     *    M," of ", SEA_TIDE%NUM_TH
          CALL GET_SL_TIDAL_HARMONICS( M, SEA_TIDE, PT,
     *         TH_AMP(1,M), TH_PHS(1,M),
     *         TH_AMP_SIN(1,M), TH_AMP_COS(1,M),
     *         IERR, LOG=NF_LOG, CUT_EXTREMUMS=0,PARAM=iparam)

        WRITE(NF_LOG,*)" @@@ GET_SL_TIDAL_HARMONICS RETURN ", IERR,
     *  ", ",TRIM(TW_NAMES(IABS(NUM_IN_TMOD(M)))),", param ",iparam
        WRITE(NF_LOG,*)" @@@ GET_SL_TIDAL_HARMONICS RESULT ", PT,
     *    M, TH_AMP(1,M), TH_PHS(1,M)
        WRITE(*,*)" @@@ GET_SL_TIDAL_HARMONICS RETURN ", IERR,
     *  ", ",TRIM(TW_NAMES(IABS(NUM_IN_TMOD(M)))),", param ",iparam
        WRITE(*,*)" @@@ GET_SL_TIDAL_HARMONICS RESULT ", PT,
     *    M, TH_AMP(1,M), TH_PHS(1,M)
          if(ierr .ne. 0)stop
          enddo
        ENDDO
        write(*,*)"FREE_TIDE_HARMONICS: "
        CALL FREE_TIDE_HARMONICS(SEA_TIDE, LOG=NF_LOG)                !MOD
        WRITE(NF_LOG,*)" @@@ FREE_TIDE_HARMONICS"
        WRITE(*,*)" @@@ FREE_TIDE_HARMONICS"
C-----------------------------------------------------------
C     (1.6.2) get fixed equivalent tide height in local arrays
C             and estimate fixed geographical shapes for tidal forcing
C-----------------------------------------------------------

        write(*,*)"EQV TIDE 1: "
        IF(SEA_TIDE%NUM_TH == 1)THEN
          M_EQ_MAX = 1          !JUST ONE HARMONIC
        ELSE
          M_EQ_MAX = NHARM_ALL  !ALL DEFINED
        ENDIF
        DO M = 1,M_EQ_MAX
          JT = M
          IF(M_EQ_MAX == 1) JT = IABS(NUM_IN_TMOD(1))   !JUST ONE HARMONIC

          ITIDE_TYPE = NINT(TIDE_TYPE(JT))
          IF(ITIDE_TYPE == 1)THEN
              !DUIRNAL
            H_EQ_LAT(M,1) = H_EQ_CM(JT)*TP_ALPHA(JT)*
     *                      SIN(PT%Y*2.*RAD_IN_D)
          ELSE IF(ITIDE_TYPE == 2)THEN
              !SEMI-DUIRNAL
            H_EQ_LAT(M,1) = H_EQ_CM(JT)*TP_ALPHA(JT)*
     *                      COS(PT%Y*RAD_IN_D)**2
          ELSE
              !LONG-PERIOD
            H_EQ_LAT(M,1) = H_EQ_CM(JT)*TP_ALPHA(JT)*
     *                      (1.5*COS(PT%Y*RAD_IN_D)**2 - 1.)
          ENDIF
        ENDDO

        WRITE(*,*)" @@@ EQUIV_TIDE_LAT ",
     *    H_EQ_LAT(1:M_EQ_MAX,1)

C***********************************************************
C
C      (2) ESTIMATION OF TIDAL PARAMETERS
C          FOR GIVEN UTC TIME IN GIVEN POINT
C          GLB_UTC_TIME must BE SPECIFIED AS DOUBLE PRECISION
C          VALUE OF SECONDS PASSED FROM 00:00 01.01.1870
C          AS DEFINED BY FUNCTION TSEC FROM  CURRENT UTC CALENDAR DATE:
C***********************************************************
        DT = TDATE(1950,1,1,0,0,0)

        TM = 0.
        DO WHILE(GLB_UTC_TIME < UTC_STOP_TIME+1D0)

          CALL GET_TIDAL_EQ_PARAM(SEA_TIDE,GLB_UTC_TIME, PH, F,
     *            PHA, FA)
!          if(ierr .ne. 0)stop
          DO M = 1,SEA_TIDE%NUM_TH
            CTDPH(M)=COS(PH(M)*RAD_IN_D)
            STDPH(M)=SIN(PH(M)*RAD_IN_D)
          ENDDO
C-----------------------------------------------------------
C      (2.1) tidal sea level anomaly in point *
C-----------------------------------------------------------
          H_TMP = 0.
          DO M = 1,SEA_TIDE%NUM_TH
            H_TMP = H_TMP +
c     *        F(M)*TH_AMP(1,M)*COS((PH(M)-TH_PHS(1,M))*RAD_IN_D)
     *        F(M)*(TH_AMP_SIN(1,M)*STDPH(M)+TH_AMP_COS(1,M)*CTDPH(M))
          ENDDO

          IF(TM < 1.)THEN
            write(*,*)"GET_TIDAL_EQ_PARAM: ",tm
            WRITE(NF_LOG,*)" @@@ GET_TIDAL_EQ_PARAM RETURN ", IERR
            WRITE(NF_LOG,*)" @@@ GET_TIDAL_EQ_PARAM RESULT ",
     *        TDT(GLB_UTC_TIME),
     *        PH(1:SEA_TIDE%NUM_TH), F(1:SEA_TIDE%NUM_TH)
            WRITE(NF_LOG,*)" @@@ TIDAL [SEA LEVEL-UVH] ANOMALY IS ",
     *        H_TMP
            WRITE(*,*)" @@@ GET_TIDAL_EQ_PARAM RETURN ", IERR
            WRITE(*,*)" @@@ GET_TIDAL_EQ_PARAM RESULT ", DT,
     *        PH(1:SEA_TIDE%NUM_TH), F(1:SEA_TIDE%NUM_TH)
            WRITE(*,*)" @@@ TIDAL [SEA LEVEL-UVH] ANOMALY IS ", H_TMP

            WRITE(NF_LOG,*)" TIME, LEVEL"
          ENDIF

C-----------------------------------------------------------
C      (2.2) equlibrium tide sea level anomaly in point *
C-----------------------------------------------------------
          H_EQV = 0.
          H_EQVL = 0.
          DO M=1,M_EQ_MAX
            JT = M
            IF(M_EQ_MAX == 1)THEN !SINGLE HARMONIC
              JT = IABS(NUM_IN_TMOD(M))
              TD_BETA = TP_BETA(JT)
            ELSE
              TD_BETA = 0.95
            ENDIF
            FF = FA(JT)
            PP = PHA(JT)
            TT = TIDE_TYPE(JT)
            TMP = FF*H_EQ_LAT(M,1)*
     *              COS((PP+TT*PT%X)*RAD_IN_D)
            H_EQV = H_EQV + TMP
! test single long-period wave::
!            IF(NINT(TT) == 0 .and. m == 20)then
            IF(NINT(TT) == 0)then
              H_EQVL = H_EQVL + TMP
!              write(*,*)'F:',M,FF,PP,TT,H_EQ_LAT(M,1),COS(PP*RAD_IN_D)
            endif
          ENDDO

          WRITE(NF_LOG,*)TM/3600.,", ",H_TMP,", ",H_EQV,", ",
     *    H_EQVL
          TM = TM + TIME_STEP
          GLB_UTC_TIME = GLB_UTC_TIME + TIME_STEP
        ENDDO

C***********************************************************
c
c    (3) clean-up after all work done
c
C***********************************************************

          CALL FREE_TIDE(SEA_TIDE, LOG=NF_LOG)
          WRITE(NF_LOG,*)" @@@ FREE_TIDE "
          WRITE(*,*)" @@@ FREE_TIDE "
          deallocate(
     *  TH_AMP,
     *  TH_PHS,
     *  TH_AMP_SIN,
     *  TH_AMP_COS,
     *  H_EQ_LAT,
     *  FA,
     *  PHA,
     *  F,
     *  PH,
     *  STDPH,
     *  CTDPH, stat=ierr)

          CLOSE(NF_LOG)
          END PROGRAM
!#ifdef VFD
!MS$  ENDIF
!#endif

#endif
