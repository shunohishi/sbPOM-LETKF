
c
c     NEC SX-8 fortran compiler do not support integer*1 type,
c     but character*1 work as c unsigned int.
c
#ifdef use_int1
#if use_int1 == 0
#undef use_int1
#endif
#endif

C--------------------------------------------------------------------------
      MODULE COM_SERV
C--------------------------------------------------------------------------
C
C      SOME COMMON PHYSICAL PARAMETERS, VARIABLES DEFINITIONS,
C      INTERFACES AND SERVIVE (TIME) PROCEDURES
C      INCLUDE DATA STRUCTURE DEFINITIONS AND SOME CONSTANTS
C      FOR GRIB DATA ACCESS AND OIL SPILL PROBLEMS
C
C      VERSION 4.00
C
C      MODIFICATION - NOVEMBER, 2002
C
C      MOVED "DAYS" IN PARAMETERS AND BUG FOR DATE RESTORATION FOR LARGE
C      YEAR NUMBERS (AFTER YEAR 4796!) REMOVED
C
C      MODIFIED DEC 28 2004:
C
c      COUNT MORE CORRECTLY LEEP YEARS, SKIP 100s EXEPT 400s.
C      - I MET PROBLEM WITH GRADS ON FEB 29, 2100:
C      NO SUCH DATE IN THE GREGORIAN CALENDAR.
C      IT MAKES THE USE OF PREVIOUSLY SAVED DOUBLE TIME VALUES
C      INCOMPATIBLE WITH NEW VERSION. TO AVOID CONFLICTS,
C      MODULE NAME WAS CHANGED
C
C      20061125 - moved here StrLowCase & StrUpCase subroutines
C      20070907 - added DATE2FNAME subroutine
C
C      SERGEY M. VARLAMOV
C
C--------------------------------------------------------------------------
          IMPLICIT NONE
          SAVE

          TYPE TDATE
              INTEGER IYEAR
              INTEGER IMONTH
              INTEGER IDAY
              INTEGER IHOUR
              INTEGER IMIN
              INTEGER ISEC
          END TYPE

          TYPE TPOINT
              REAL X
              REAL Y
              REAL Z
          END TYPE

          TYPE TSECTION
              TYPE(TPOINT)  INI_POINT
              TYPE(TPOINT)  FIN_POINT
              CHARACTER(4)  CODE
C              integer       reserved_for_alignment
          END TYPE

          TYPE TSPOINT
              TYPE(TPOINT)  POINT_POS
              CHARACTER(4)  CODE
C              integer       reserved_for_alignment
          END TYPE
C
C         SOME GLOBAL PARAMETERS, IN METRIC SYSTEM
C
          REAL, PARAMETER::
     *                      PI=3.1415926535897932384626433832795,
     *                      D_IN_RAD = 180./PI,
     *                      RAD_IN_D = PI/180.,
     *                      HOUR     = 3600.,
     *                      DAY      = HOUR * 24.
C
C          EARTH RADIUS & G, M & M/S^2
C
          REAL, PARAMETER:: R_EARTH = 6.371E6,
     *                      GRAVITY = 9.81
C
C         RETURN PARAMETERS DEFINITIONS
C         FOR VECTOR INTERPOLATION SUBROUTINES -
C         LINEAR COMPONENTS INTERPOLATION,
C         VECTOR ROTATION INTERPOLATION OR
C         SMALL VALUE
C
          INTEGER, PARAMETER::    INTERP_ZERO = 1,
     *                            INTERP_COMP = 2,
     *                            INTERP_VECT = 3
C
C         PROVIDE HERE CONSTANTS FOR INTERPOLATION
C         METHOD TO BE USED BY USER
C         FOR SPLINE, LINEAR AND NEAREST POINT METHODS
C         THEY ARE EQUEL TO NUMBER OF POINTS USED
C         FOR TIME INTERPOLATION.
C
          INTEGER, PARAMETER::    INTERP_SPLINE  = 4,
     *                            INTERP_INVDIST = 3,
     *                            INTERP_LINEAR  = 2,
     *                            INTERP_NEAREST = 1

          CHARACTER( * ), PRIVATE, PARAMETER ::
     *      LOWER_CASE = 'abcdefghijklmnopqrstuvwxyz'
          CHARACTER( * ), PRIVATE, PARAMETER ::
     *      UPPER_CASE = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

          INTERFACE

            INTEGER(4) FUNCTION DATE2NAME
     *      (YEAR,MONTH,DAY,HOUR,MIN,YDAY,STRLEN,PATH,PATH_FORMAT)
              INTEGER(4),INTENT(IN)::  YEAR,MONTH,DAY,HOUR,MIN,
     *                                 YDAY,STRLEN
#ifdef use_int1
              INTEGER(1),INTENT(OUT):: PATH(*)
              INTEGER(1),INTENT(IN)::  PATH_FORMAT(*)
#else
              CHARACTER,INTENT(OUT):: PATH(*)
              CHARACTER,INTENT(IN)::  PATH_FORMAT(*)
#endif
            END FUNCTION

            INTEGER(4) FUNCTION   NAME2DATE
     *      (YEAR,MONTH,DAY,HOUR,MIN,PATH,PATH_FORMAT)
              INTEGER(4),INTENT(OUT):: YEAR,MONTH,DAY,HOUR,MIN
#ifdef use_int1
              INTEGER(1),INTENT(IN)::  PATH(*),PATH_FORMAT(*)
#else
              CHARACTER,INTENT(IN)::  PATH(*),PATH_FORMAT(*)
#endif
            END FUNCTION

            INTEGER(4) FUNCTION MODIFY_PATH_MASK_M(PATH,MAXLEN)
              INTEGER(4),INTENT(IN):: MAXLEN
#ifdef use_int1
              INTEGER(1),INTENT(INOUT)::  PATH(*)
#else
              CHARACTER,INTENT(INOUT)::  PATH(*)
#endif
            END FUNCTION

              SUBROUTINE IBCCCU (F,X,NX,Y,NY,C,IC,WK,IER)
                  IMPLICIT NONE
                  INTEGER, INTENT(IN )::      NX,NY,IC
                  INTEGER, INTENT(OUT)::      IER
                  REAL,    INTENT(IN )::      F(IC,*), X(*), Y(*)
                  REAL,    INTENT(OUT)::      C(2,IC,*)
                  REAL,    INTENT(OUT)::      WK(*)
              END SUBROUTINE

              SUBROUTINE GET_BSPLINES (F,X,NX,Y,NY,C,IC,IER)
                IMPLICIT NONE
                INTEGER, INTENT(IN )::    NX,NY,IC
                INTEGER, INTENT(OUT)::    IER
                REAL,    INTENT(IN )::    F(IC,*),X(*),Y(*)
                REAL,    INTENT(OUT)::    C(2,IC,*)
              END SUBROUTINE

              SUBROUTINE GET_BSPLINES_S (F,X,NX,Y,NY,C,IC,IER)
                IMPLICIT NONE
                INTEGER, INTENT(IN )::    NX,NY,IC
                INTEGER, INTENT(OUT)::    IER
                REAL,    INTENT(IN )::    F(IC,*),X(*),Y(*)
                REAL,    INTENT(OUT)::    C(2,IC,*)
              END SUBROUTINE

              SUBROUTINE DBCEVL (X,NX,Y,NY,C,IC,XL,YL,PDS,IER)
                  IMPLICIT NONE
                  INTEGER, INTENT(IN )::      NX,NY,IC
                  INTEGER, INTENT(OUT)::      IER
                  REAL,    INTENT(IN )::      X(*),Y(*)
                  REAL,    INTENT(IN )::      C(2,IC,*),XL,YL
                  REAL,    INTENT(OUT)::      PDS(*)
              END SUBROUTINE

              SUBROUTINE BCSEVL_S (X,NX,Y,NY,C,IC,XL,YL,ZL,IER)
                  IMPLICIT NONE
                  INTEGER, INTENT(IN ):: NX,NY,IC
                  INTEGER, INTENT(OUT):: IER
                  REAL,    INTENT(IN ):: X(*),Y(*),C(2,IC,*),XL,YL
                  REAL,    INTENT(OUT):: ZL
              END SUBROUTINE

              SUBROUTINE BCSEVL_V (X,NX,Y,NY,C,IC,XL,NXL,YL,NYL,
     *                             PDS,NXD,IER)
                  IMPLICIT NONE
                  INTEGER, INTENT(IN )::    NX,NY,IC,NXL,NYL,NXD
                  INTEGER, INTENT(OUT)::    IER
                  REAL,    INTENT(IN )::    X(*),Y(*),C(2,IC,*),
     *                                      XL(*),YL(*)
                  REAL,    INTENT(OUT)::    PDS(NXD,*)
              END SUBROUTINE

              SUBROUTINE BCSEVL_VE (X,NX,Y,NY,C,IC,XL,NXL,YL,NYL,
     *                             PDS,NXD,IER)
                  IMPLICIT NONE
                  INTEGER, INTENT(IN )::    NX,NY,IC,NXL,NYL,NXD
                  INTEGER, INTENT(OUT)::    IER
                  REAL,    INTENT(IN )::    X(*),Y(*),C(2,IC,*),
     *                                      XL(*),YL(*)
                  REAL,    INTENT(OUT)::    PDS(NXD,*)
              END SUBROUTINE

              SUBROUTINE BCSEVL_VM(X,NX,Y,NY,C,IC,XL,NXL,YL,NYL,
     *                    PDS,NXD,MASK,UNDEF,IER)
                IMPLICIT NONE
                INTEGER, INTENT(IN )::  NX,NY,IC,NXL,NYL,NXD
                INTEGER, INTENT(OUT)::  IER
                REAL,    INTENT(IN )::  X(*),Y(*),C(2,IC,*),
     *                                  XL(*),YL(*),
     *                                  MASK(NXD,*), UNDEF
                REAL,    INTENT(OUT)::  PDS(NXD,*)
              END SUBROUTINE

              SUBROUTINE BCSEVL_EX(X,NX,Y,NY,C,IC,XL,YL,ZL,IER,CUT_EXTR)
                  IMPLICIT NONE
                  INTEGER, INTENT(IN ):: NX,NY,IC,CUT_EXTR
                  INTEGER, INTENT(OUT):: IER
                  REAL,    INTENT(IN ):: X(*),Y(*),C(2,IC,*),XL,YL
                  REAL,    INTENT(OUT):: ZL
              END SUBROUTINE

              SUBROUTINE BCSEVL_PVEXP(X,NX,Y,NY,C,IC,PT,NP,IER,
     *                                GLOBX,UNDEF,CUT_EXTR,MASK,VALUE)
                IMPLICIT NONE
                TYPE TPOINT
                  REAL X
                  REAL Y
                  REAL Z
                END TYPE
                INTEGER, INTENT(IN)::           NX,NY,IC,NP
                REAL, INTENT(IN)::              X(*),Y(*),C(2,IC,*)
                TYPE(TPOINT), INTENT(INOUT)::   PT(*)
                INTEGER, INTENT(OUT)::          IER
                INTEGER, OPTIONAL, INTENT(IN):: GLOBX, CUT_EXTR
                LOGICAL, OPTIONAL, INTENT(IN):: MASK(*)
                REAL, OPTIONAL, INTENT(IN)::    UNDEF
                REAL, OPTIONAL, INTENT(OUT)::   VALUE(*)
              END SUBROUTINE BCSEVL_PVEXP


              SUBROUTINE VECTOR_LININTERP(U, V, R, WX, WY,
     *                   ZERO_WIND, MAX_ANGLE_ROTATE, METHOD)
C
C                  MAX_ANGLE_ROTATE - IN RADIANS,
C                  ZERO_WIND - IN USED VECTOR UNITS;
C
                  IMPLICIT NONE
                  REAL,    INTENT( IN)::      U(*), V(*), R, ZERO_WIND,
     *                                        MAX_ANGLE_ROTATE
                  REAL,    INTENT(OUT)::      WX, WY
                  INTEGER, INTENT(OUT)::      METHOD
              END SUBROUTINE

              SUBROUTINE VECTOR_SPLINEINTERP(U, V, T, R, WX, WY,
     *            ZERO_WIND, MAX_ANGLE_ROTATE, METHOD)
                  IMPLICIT NONE
                  REAL,    INTENT( IN)::    U(*), V(*), R, T(*),
     *                                      ZERO_WIND, MAX_ANGLE_ROTATE
                  REAL,    INTENT(OUT)::    WX, WY
                  INTEGER, INTENT(OUT)::    METHOD
C
C             METHOD IS ONE OF PREDIFINED CONSTANTS:
C                  INTERP_COMP (COMPONENTS),
C                  INTERP_VECT (INTERPOLATE AMPLITUDE AND PHASE) OR
C                  INTERP_ZERO (SMALL COMPONENTS WERE INTERPOLATED)
C
              END SUBROUTINE

              SUBROUTINE INVDIST_2DINTERPOL(X,NX,Y,NY,XL,YL,Z,NXDIM,ZL,
     *                                      IER)
C
C      SIMPLE INTERPOLATION - INVERSE DISTANCE WEIGHTED
C      FROM 4 CLOSEST GRID POINT VALUES.
C
                  IMPLICIT NONE
                  INTEGER, INTENT(IN )::      NX,NY,NXDIM
                  INTEGER, INTENT(OUT)::      IER
                  REAL,    INTENT(IN )::      X(*),Y(*),Z(NXDIM,*),XL,YL
                  REAL,    INTENT(OUT)::      ZL

              END SUBROUTINE

              SUBROUTINE NEARPT_2DINTERPOL (X,NX,Y,NY,XL,YL,Z,NXDIM,ZL,
     *                                      IER)
C
C                  SIMPLE INTERPOLATION FROM CLOSEST GRID POINT.
C
                  IMPLICIT NONE
                  INTEGER, INTENT(IN )::      NX,NY,NXDIM
                  INTEGER, INTENT(OUT)::      IER
                  REAL,    INTENT(IN )::      X(*),Y(*),Z(NXDIM,*),XL,YL
                  REAL,    INTENT(OUT)::      ZL

              END SUBROUTINE

              SUBROUTINE SPLINE_INTERPOLATE_1D(X,Y,XX,YY,IERR)
                  IMPLICIT NONE
                  INTEGER, INTENT(OUT)::      IERR
                  REAL, INTENT(IN)::          X(*),Y(*),XX
                  REAL, INTENT(OUT)::         YY
              END SUBROUTINE

              SUBROUTINE      GRADIENTS_FOR_1D_SPLINE_GLM
     *            (X, Y, GRAD, NP, INDT, IGRAD, IERR)
                  IMPLICIT NONE
                  INTEGER, INTENT(IN)::   NP, INDT(*), IGRAD
                  INTEGER, INTENT(OUT)::  IERR
                  REAL, INTENT(IN)::      X(*),Y(NP,*)
                  REAL, INTENT(OUT)::     GRAD(NP,*)
              END SUBROUTINE

              SUBROUTINE GRADIENTS_FOR_1D_SPLINE_GLM_M
     *        (X, Y, GRAD, NX, NY, NXDIM, NYDIM, INDT, IGRAD, IERR)
                  IMPLICIT        NONE
                  INTEGER, INTENT(IN)::   NX, NY, NXDIM, NYDIM,
     *                                    INDT(*), IGRAD
                  INTEGER, INTENT(OUT)::  IERR
                  REAL, INTENT(IN)::      X(*),Y(NXDIM, NYDIM,*)
                  REAL, INTENT(OUT)::     GRAD(NXDIM, NYDIM,*)
              END SUBROUTINE

              SUBROUTINE      GRAD_1D_SPLINE_GLM_CYCL
     *            (X, Y, GRAD, NP, INDT, IGRAD, NCYCLE, IERR)
                  IMPLICIT NONE
                  INTEGER, INTENT(IN)::   NP, INDT(*), IGRAD, NCYCLE
                  INTEGER, INTENT(OUT)::  IERR
                  REAL, INTENT(IN)::      X(*),Y(NP,*)
                  REAL, INTENT(OUT)::     GRAD(NP,*)
              END SUBROUTINE

              SUBROUTINE SPLINE_INTERPOLATE_1D_GLM
     *            (X, Y, GRAD, NP, INDT, IGRAD, XX, YY, IERR)
                  IMPLICIT NONE
                  INTEGER, INTENT(IN)::   NP, INDT(*), IGRAD
                  INTEGER, INTENT(OUT)::  IERR
                  REAL, INTENT(IN)::      X(*),Y(NP,*),GRAD(NP,*),XX
                  REAL, INTENT(OUT)::     YY(*)
              END SUBROUTINE

              SUBROUTINE SPLINE_INTERPOLATE_1D_GLM_M
     *        (X,Y,GRAD,NX,NY,NXDIM,NYDIM,INDT,IGRAD,XX,YY,IERR)
                  IMPLICIT NONE
                  INTEGER, INTENT(IN)::   NX,NY,NXDIM,NYDIM,
     *                                    INDT(*),IGRAD
                  INTEGER, INTENT(OUT)::  IERR
                  REAL, INTENT(IN)::      X(*),Y(NXDIM, NYDIM,*),
     *                                    GRAD(NXDIM, NYDIM, *),XX
                  REAL, INTENT(OUT)::     YY(NXDIM, *)
              END SUBROUTINE

              SUBROUTINE SPLINE_INTERPOLATE
     *            (XIN,YIN,NIN,XOUT,YOUT,NOUT,UNDEF,IRET)
                  IMPLICIT NONE
                  INTEGER, INTENT(IN)::   NIN,NOUT
                  INTEGER, INTENT(OUT)::  IRET
                  REAL, INTENT(IN)::      XIN(*),YIN(*),
     *                                    XOUT(*),UNDEF
                  REAL, INTENT(OUT)::     YOUT(*)
              END SUBROUTINE

              SUBROUTINE LINEAR_INTERPOLATE
     *            (XIN,YIN,NIN,XOUT,YOUT,NOUT,UNDEF)
                  IMPLICIT NONE
                  INTEGER, INTENT(IN)::   NIN,NOUT
                  REAL, INTENT(IN)::      XIN(*),YIN(*),
     *                                    XOUT(*),UNDEF
                  REAL, INTENT(OUT)::     YOUT(*)
              END SUBROUTINE

              SUBROUTINE EXTRAPOLATE_BY_MASK(F,MASK,NX,NY,NX_D,IRET)
                  IMPLICIT NONE
                  INTEGER, INTENT(IN )::      NX, NY, NX_D
                  INTEGER, INTENT(OUT)::      IRET
                  REAL,  INTENT(INOUT)::      F(NX_D,*)
                  REAL,    INTENT(IN )::      MASK(NX_D,*)
              END SUBROUTINE

              SUBROUTINE EXTRAPOLATE_BY_MASK_V5(F,MASK,NX,NY,NX_F,IRET)
                  IMPLICIT NONE
                  INTEGER, INTENT(IN )::      NX, NY, NX_F
                  INTEGER, INTENT(OUT)::      IRET
                  REAL,  INTENT(INOUT)::      F(NX_F,*)
                  REAL,    INTENT(IN )::      MASK(NX,*)
              END SUBROUTINE

              SUBROUTINE MASK_REPLACE_MISSING(IRET, DAT, NX, NY, NXDIM,
     *            UNDEF_MIN, UNDEF_MAX, MASK)
                  IMPLICIT NONE
                  INTEGER, INTENT(IN )::      NX, NY, NXDIM
                  INTEGER, INTENT(OUT)::      IRET
                  REAL,  INTENT(INOUT)::      DAT(NXDIM,*)
                  REAL,    INTENT(IN ),OPTIONAL::
     *                                        MASK(NX,*),
     *                                        UNDEF_MIN, UNDEF_MAX
              END SUBROUTINE

              SUBROUTINE MASK_REPLACE_MISSING_M
     *        (IRET, R_MEAN, DAT, NX, NY, NXDIM,
     *        UNDEF_MIN, UNDEF_MAX, MASK)
                IMPLICIT NONE
                INTEGER, INTENT(IN   )::  NX, NY, NXDIM
                INTEGER, INTENT(OUT  )::  IRET
                REAL,    INTENT(INOUT)::  DAT(NXDIM,*)
                REAL,    INTENT(OUT)::    R_MEAN
                REAL,    INTENT(IN   ),OPTIONAL::
     *                                    MASK(NX,*),
     *                                    UNDEF_MIN, UNDEF_MAX

              END SUBROUTINE

              SUBROUTINE RANLUX(RVEC, LENV)
                  IMPLICIT NONE
                  INTEGER, INTENT( IN)::      LENV
                  REAL,    INTENT(OUT)::      RVEC(*)
              END SUBROUTINE

              SUBROUTINE RANLUX0(RVEC, LENV)
                  IMPLICIT NONE
                  INTEGER, INTENT( IN)::      LENV
                  REAL,    INTENT(OUT)::      RVEC(*)
              END SUBROUTINE

              SUBROUTINE LUBKSB(A,N,NP,INDX,B)
                  IMPLICIT NONE
                  INTEGER, INTENT(IN   )::    N,NP, INDX(*)
                  REAL,    INTENT(IN   )::    A(NP,*)
                  REAL,    INTENT(INOUT)::    B(*)
              END SUBROUTINE

              SUBROUTINE LUDCMP(A, N, NP, INDX, D, ERR)
                  IMPLICIT NONE
                  INTEGER, INTENT(IN   )::    N,NP
                  INTEGER, INTENT(OUT  )::    ERR
                  INTEGER, INTENT(INOUT)::    INDX(*)
                  REAL,    INTENT(OUT  )::    D
                  REAL,    INTENT(INOUT)::    A(NP,*)
              END SUBROUTINE

          END INTERFACE

        CONTAINS
C***************************************************************************C
C                                                                           C
          DOUBLE PRECISION FUNCTION TSEC(TM)
C                                                                           C
C             CALCULATE SECONDS FROM 00:00:00 01.01.1880                    C
C                                                                           C
C             CHECK THAT YEAR IS AFTER 1879 AND MONTH BETWEEN 1 AND 12.     C
C             DAY NUMBER, HOURS, MINUTES AND SECONDS ARE NOT CHECKED -      C
C             THEY ARE ADDED TO THE TIME AS IS. IT COULD BE USED TO         C
C             ESTIMATE DATE FROM TIME SHIFT ON FIXED NUMBER OF DAYS, HOURS, C
C             MINUTES AND SECONDS. JUST CHANGE DATE, CALL TSEC AND          C
C             TRANSFORM IT INTO CORRECT DATE BY CALLING TDT.                C
C                                                                           C
C***************************************************************************C

              IMPLICIT NONE
              TYPE(TDATE), INTENT(IN)::      TM

              DOUBLE PRECISION,PARAMETER::
     *            MN = 60.,     !SEC
     *            H = MN*60.,   !3600 SEC
     *            D = H*24.,    !86400 SEC
     *            Y = D*365.,   !31536000 SEC
     *            DAYS(12) = (/31.,28.,31.,30.,31.,30.,
     *                         31.,31.,30.,31.,30.,31./)
              INTEGER IYPASSED,IYP4,IYP100,IYP400,I,IM

              TSEC = 0D0
              IYPASSED = TM%IYEAR-1880
              IM = TM%IMONTH-1
C
C             CHECK THAT YEAR IS AFTER 1879 AND MONTH BETWEEN 1 AND 12.
C
              IF(IYPASSED < 0 .OR. IM < 0 .OR. IM >= 12)RETURN

              IYP4   = INT((TM%IYEAR-1877)/4)
              IYP100 = INT((TM%IYEAR-1801)/100)
              IYP400 = INT((TM%IYEAR-1601)/400)

              TSEC = DBLE(IYPASSED)*Y + DBLE(IYP4-IYP100+IYP400)*D

              IF( IM > 0)THEN
                  DO I = 1, IM
                      TSEC = TSEC + DAYS(I)*D
                  ENDDO
C
C                 IM == 2 => MARCH; FEBRUARY PASSED, WAS IT A LEAP YEAR?
C
                  IF(IM >= 2 .AND. MOD(TM%IYEAR,4) == 0)THEN
C
C                   SOME YEARS ARE SKIPPED FROM LEAP:
C
                    IF(MOD(TM%IYEAR,100) /= 0 .OR.
     *                 MOD(TM%IYEAR,400) == 0      )THEN
                      TSEC = TSEC+D
                    ENDIF
                  ENDIF
              ENDIF

              TSEC = DBLE(TM%IDAY-1)*D + TSEC
              TSEC = DBLE(TM%IHOUR)*H  + TSEC
              TSEC = DBLE(TM%IMIN)*MN + TSEC
              TSEC = DBLE(TM%ISEC)     + TSEC
          END FUNCTION TSEC
c
C***************************************************************************C
C
          FUNCTION TDT(TM) RESULT (TD)
C
C             CALCULATE THE DATE FROM SECONDS PASSED AFTER
C             00:00:00 01.01.1880
C             IF TM IS NEGATIVE - RETURN DATE(1880,1,1,0,0,0)
C
C             ON SOME SYSTEMS, LIKE LINUX INTEL COMPILER ON PC,
C             ROUNDING PROBLEMS GENERATED OUTPUT WITH 60 SEC:
C             TO SOLVE THIS PROBLEM, KEEP ACCURACY UP TO 1 SEC,
C             ADD 0.01 SEC FOR ALIGNMENT (JANUARY 2004).
C
C***************************************************************************C
              IMPLICIT NONE
              TYPE(TDATE):: TD
              DOUBLE PRECISION, INTENT(IN):: TM

              DOUBLE PRECISION PREV, PASSED, TMW, ADD_FEBRUARY
              INTEGER IYEAR, IMONTH, IDAY, IHOUR, IMIN, I_PW,
     *                IYP4, IYP100, IYP400, ICYEAR

              DOUBLE PRECISION,PARAMETER::
     *            MN = 60.,     !SEC
     *            H = MN*60.,   !3600 SEC
     *            D = H*24.,    !86400 SEC
     *            Y = D*365.,   !31536000 SEC
     *            DAYS(12) = (/31.,28.,31.,30.,31.,30.,
     *                         31.,31.,30.,31.,30.,31./)

              TD   = TDATE(1880,1,1,0,0,0)

              I_PW = PRECISION(TMW)
              TMW  = 9.*10.**I_PW  ! ESTIMATION FOR LARGEST INTEGER PRESENTED IN TM!

              IF( TM >= 0D0 .AND. TM <= TMW)THEN

                TMW = TM + 1D-2
C
C               APROXIMATELY ESTIMATE NUMBER OF YEARS PASSED,
C               (UNDERESTIMATE IN ANY CASE AS 1 IS SUBTRACTED)::
C
                IYEAR = INT(TMW/(Y+D*0.25))-1
                IF(IYEAR > 0)THEN
                  ICYEAR = 1880 + IYEAR
                  IYP4   = INT((ICYEAR-1877)/4)
                  IYP100 = INT((ICYEAR-1801)/100)
                  IYP400 = INT((ICYEAR-1601)/400)
                  PASSED = DBLE(IYEAR)*Y + DBLE(IYP4-IYP100+IYP400)*D
                ELSE
                  IYEAR  = 0
                  PASSED = 0D0
                ENDIF
                PREV = PASSED
C
C               FIND TRUE NUMBER OF CALENDAR YEARS PASSED
C               THE CYCLE WILL BE PASSED AT LEAST ONE TIME
C
                DO WHILE( TMW >= PASSED)
                  PREV  = PASSED
                  PASSED = PASSED + Y
                  ICYEAR = 1880 + IYEAR
                  IF( MOD(ICYEAR,4) == 0 .AND.
     *               (MOD(ICYEAR,100) /= 0 .OR.
     *                MOD(ICYEAR,400) == 0 ))THEN
                    PASSED = PASSED + D
                  ENDIF
                  IYEAR = IYEAR + 1 !IT WILL BE NEXT YEAR
                ENDDO
                TMW = TMW - PREV
                IYEAR = IYEAR - 1

                ICYEAR = 1880 + IYEAR
                ADD_FEBRUARY = 0D0
                IF( MOD(ICYEAR,4) == 0 .AND.
     *             (MOD(ICYEAR,100) /= 0 .OR.
     *              MOD(ICYEAR,400) == 0 ))ADD_FEBRUARY = D
C
C               FOR MONTHS GO DIRECTLY -
C               DO NOT ESTIMATE APPROXIMATE START VALUE::
C
                PREV   = 0D0
                PASSED = 0D0
                IMONTH = 0
                DO I_PW=1,12
                  PASSED = PASSED + DAYS(I_PW)*D
                  IF(I_PW == 2) PASSED = PASSED + ADD_FEBRUARY
c                  IF(PASSED > TMW)THEN
                  IF(IMONTH < 1 .AND. PASSED > TMW)THEN
                      TMW = TMW - PREV
                      IMONTH = I_PW
C                      EXIT
                  ENDIF
                  PREV = PASSED
                ENDDO
C
C               IT IS IMPOSSIBLE IN NORMAL CASE, BUT...::
C
                IF(IMONTH <= 12) THEN
                    TD%IYEAR  = IYEAR + 1880
                    TD%IMONTH = IMONTH

                    IDAY      = INT(TMW/D)
                    TD%IDAY   = IDAY+1
                    TMW       = TMW-D*IDAY

                    IHOUR     = INT(TMW/H)
                    TD%IHOUR  = IHOUR
                    TMW       = TMW-H*IHOUR

                    IMIN      = INT(TMW/MN)
                    TD%IMIN   = IMIN
                    TMW       = TMW-MN*IMIN

                    TD%ISEC   = TMW
                  ENDIF
              ENDIF
          END FUNCTION TDT

          SUBROUTINE CLEAN_STR_C(STRING)
            CHARACTER(*),INTENT(INOUT):: STRING
            STRING = CHAR(0)
          END SUBROUTINE

          REAL FUNCTION CYCL_COR(LON_BASE,LON_TEST)
            IMPLICIT NONE
            REAL, INTENT(IN):: LON_BASE,LON_TEST
            REAL D1, D2, D, CR, CCR2
            INTEGER I
            REAL, PARAMETER:: EPS = 1E-7
C
C           IF ABS DIFFERRENCE IS 180+I*360 - IT IS
C           NOT FULLY DEFINED CASE.
C           IN THIS CASE TAKE 'I' SUCH THAT DIFFERENCE
C           IS SMALEST (180) AND HAVE SIGN AS FOR INPUT
C           DATA.
C
            D1 = HUGE(D1)
            D2 = D1
            CYCL_COR = 0.
            CCR2 = 0.
            DO I = 1,11
              CR = 360.*(I-6)
              D = ABS(LON_BASE-(LON_TEST+CR))
              IF(D <= D1+EPS)THEN !COULD BE +180, -180
                D2 = D1
                D1 = D
                CCR2 = CYCL_COR
                CYCL_COR = CR !ADD THIS CYCLIC CORRECTION VALUE
              ENDIF
            ENDDO
            IF(ABS(D1-D2) <= EPS)THEN  !COULD BE +180, -180
              D = (LON_BASE-LON_TEST)
              D2 = (LON_BASE-(LON_TEST+CCR2))
              IF(D*D2 > 0.) CYCL_COR = CCR2
            ENDIF

          END FUNCTION
C----------------------------------------------------------------------
C       GREATEST COMMON DIVIDER, EUCLID ALGORITHM
C       RETURN ZERO IF ONE OF NUMBER IS ZERO,
C       GREATEST POSITIVE INTEGER DIVIDER IN OTHER CASE
C----------------------------------------------------------------------
        INTEGER FUNCTION GR_COM_DIV(A,B)
          IMPLICIT NONE
          INTEGER, INTENT(IN):: A,B
          INTEGER CC,AA,BB
          AA = ABS(A)
          BB = ABS(B)
          CC = MIN(AA,BB) !SMALLER OF TWO
          IF(CC <= 1)THEN
            GR_COM_DIV = CC
            RETURN
          ENDIF
          IF(AA > BB) BB = AA   !GREATER VALUE
          DO WHILE(CC .NE. 0)
            AA = BB
            BB = CC
            CC = MOD(AA,BB)
          ENDDO
          GR_COM_DIV = BB
        END FUNCTION

********************************************************************************************
        subroutine StrUpCase ( String )
********************************************************************************************
        CHARACTER( * ), INTENT( INOUT ) :: String
        INTEGER i, n

!       Loop over string elements
        DO i = 1, LEN( String )
!         Find location of letter in lower case constant string
          n = INDEX( LOWER_CASE, String( i:i ) )
!         If current substring is a lower case letter, make it upper case
          IF ( n /= 0 ) String( i:i ) = UPPER_CASE( n:n )
        END DO
      END subroutine StrUpCase

********************************************************************************************
        subroutine StrLowCase ( String )
********************************************************************************************
        CHARACTER( * ), INTENT( INOUT ) :: String
        INTEGER :: i, n

        DO i = 1, LEN( String )
          n = INDEX( UPPER_CASE, String( i:i ) )
          IF ( n /= 0 ) String( i:i ) = LOWER_CASE( n:n )
        END DO
      END subroutine StrLowCase

********************************************************************************************
      END MODULE
C--------------------------------------------------------------------------

