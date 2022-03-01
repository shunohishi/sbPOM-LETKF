C-----------------------------------------------------------------------
C
C    BICUBIC SPLINE EVALUATORS
C
C    RESTRUCTURIZATION OF CODE FOR FORTRAN 90
C    INDEX SEARCHING
C    BASED ON IMSL ROUTINE DBCEVL
C
C-----------------------------------------------------------------------
C
C   ARGUMENTS    X      - VECTOR OF LENGTH NX. (INPUT) X MUST BE
C                           ORDERED SO THAT X(I) .LT. X(I+1) FOR
C                           I=1,...,NX-1.
C                NX     - NUMBER OF ELEMENTS IN X. (INPUT) NX MUST BE
C                           .GE. 2.
C                Y      - VECTOR OF LENGTH NY. (INPUT) Y MUST BE
C                           ORDERED SO THAT Y(J) .LT. Y(J+1) FOR
C                           J=1,...,NY-1.
C                NY     - NUMBER OF ELEMENTS IN Y. (INPUT) NY MUST BE
C                           .GE. 2.
C                         NOTE - THE COORDINATE PAIRS (X(I),Y(J)), FOR
C                           I=1,...,NX AND J=1,...,NY, GIVE THE POINTS
C                           WHERE THE FUNCTION VALUES ARE DEFINED.
C                C      - ARRAY OF SPLINE COEFFICIENTS. (INPUT)
C                           C IS OF DIMENSION 2 BY NX BY 2 BY NY.
C                           THE SPLINE COEFFICIENTS CAN BE COMPUTED BY
C                           IMSL SUBROUTINE IBCCCU.
C                           (NOTE - C IS TREATED INTERNALLY AS A
C                            2 BY NX BY 2*NY ARRAY BECAUSE CERTAIN
C                            ENVIRONMENTS DO NOT PERMIT QUADRUPLY-
C                            DIMENSIONED ARRAYS.  IN THESE
C                            ENVIRONMENTS THE CALLING PROGRAM MAY
C                            DIMENSION C IN THE SAME MANNER.)
C                IC     - SECOND DIMENSION OF ARRAY C EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT
C                           (INPUT).  IC MUST BE .GE. NX.
C                XL,YL  - (XL,YL) IS THE POINT AT WHICH THE SPLINE
C                           IS TO BE EVALUATED. (INPUT)
C                ZL      - VALUE OF THE BICUBIC SPLINE, S(X,Y),
C                           EVALUATED AT X=XL AND Y=YL. (OUTPUT)
C                XL,YL  - (XL,YL) IS THE POINT AT WHICH THE MIXED
C                           PARTIAL DERIVATIVES OF THE SPLINE ARE TO BE
C                           EVALUATED. (INPUT)
C                PDS    - VECTOR OF LENGTH 6 CONTAINING THE PARTIAL
C                           DERIVATIVES OF THE BICUBIC SPLINE, S(X,Y),
C                           EVALUATED AT X=XL AND Y=YL. (OUTPUT)
C                             PDS(1) = S(XL,YL)
C                             PDS(2) = DS/DX
C                             PDS(3) = DS/DY
C                             PDS(4) = D(DS/DX)/DY
C                             PDS(5) = D(DS/DX)/DX
C                             PDS(6) = D(DS/DY)/DY.
C                IER    - ERROR PARAMETER. (OUTPUT)
C

C-----------------------------------------------------------------------
C
      SUBROUTINE DBCEVL (X,NX,Y,NY,C,IC,XL,YL,PDS,IER)
      IMPLICIT   NONE
C
C                            SPECIFICATIONS FOR ARGUMENTS
C
      INTEGER, INTENT(IN )::    NX,NY,IC
      INTEGER, INTENT(OUT)::    IER
      REAL,    INTENT(IN )::    X(NX),Y(NY),C(2,IC,2*NY),XL,YL
      REAL,    INTENT(OUT)::    PDS(6)
C
C                            SPECIFICATIONS FOR LOCAL VARIABLES
C
      INTEGER                 I,K,KM1,KP1,KP2,LXPL,LX,LY,L,LXP1,IX,IY,J
      REAL                    HX,HY,SVX(2),SV(2),U,V,R,
     *                        SUX(2),SUY(2),SU(2),SXY(2)
      REAL,PARAMETER::        EPS = 1.0E-7
C
C     STATEMENT-FUNCTION::
C
C      REAL SPLN0,S0,SH,SP0,SPH,H,D
C     *  ,SPLN1,SPLN2
C      SPLN0(S0,SH,SP0,SPH,H,D) = S0+D*(H*SP0+D*(3.*(SH-S0)-
C     *            (SPH+2.*SP0)*H+D*(2.*(S0-SH)+(SPH+SP0)*H)))
c      SPLN1(S0,SH,SP0,SPH,H,D) = SP0+D*(6.*(SH-S0)/H-2.*
c     *            (SPH+2.*SP0)+3.*D*(2.*(S0-SH)/H+(SPH+SP0)))
C      SPLN2(S0,SH,SP0,SPH,H,D) = 6.*(SH-S0)/(H*H)-2.*
C     *            (SPH+2.*SP0)/H+D*(2.*(S0-SH)/(H*H)+(SPH+SP0)/H)*6.
C
C     PRELIMINARY INDEXES ESTIMATION WITH ASSUMPTION
C    OF CONSTANT STEPS IS INTRODUCED::
C
      R = X(NX)-X(1)
      IF(ABS(R) < EPS)THEN
        IER = 1
        RETURN
      ENDIF

      R = (XL-X(1))/R
      IF(R < 0.0 .OR. R > 1.0)THEN
        IER = 2
        RETURN
      ENDIF

      IX = INT(R*(NX-1))+1
      IF( IX == NX ) IX = NX-1

      IF( XL > X(IX)) THEN
        DO I=IX+1,NX-1
          IF (XL <= X(I)) EXIT
        ENDDO
        LX = I-1
      ELSE ! (XL <= X(IX))
        DO I=IX-1,2,-1
          IF (XL > X(I)) EXIT
        ENDDO
        LX = I
        IF( LX == 0 ) LX = 1
      ENDIF

      R = Y(NY)-Y(1)
      IF(ABS(R) < EPS)THEN
        IER = 3
        RETURN
      ENDIF

      R = (YL-Y(1))/R
      IF(R < 0.0 .OR. R > 1.0)THEN
        IER = 4
        RETURN
      ENDIF

      IY = INT(R*(NY-1)) + 1
      IF( IY == NY ) IY = NY - 1

      IF( YL > Y(IY)) THEN
        DO I=IY+1,NY-1
          IF (YL < Y(I)) EXIT
        ENDDO
        LY = I-1
      ELSE ! (YL .LE. Y(IY))
        DO I=IY-1,2,-1
          IF (YL > Y(I)) EXIT
        ENDDO
        LY = I
        IF( LY == 0 ) LY = 1
      ENDIF

      LXP1 = LX+1
      HX = X(LXP1)-X(LX)
      HY = Y(LY+1)-Y(LY)
      U = (XL-X(LX))/HX
      V = (YL-Y(LY))/HY
      K = 2*LY
      KP1 = K+1
      KP2 = K+2
      KM1 = K-1

      DO L=1,2
          LXPL = LX-1+L
          I = 2*(LY-1+L)
          J = I-1
          SUX(L) = SPLN1(C(1,LX,J),C(1,LXP1,J),C(2,LX,J),
     *        C(2,LXP1,J),HX,U)
          SXY(L) = SPLN1(C(1,LX,I),C(1,LXP1,I),C(2,LX,I),
     *        C(2,LXP1,I),HX,U)
          SU(L) = SPLN0(C(1,LX,J),C(1,LXP1,J),C(2,LX,J),
     *        C(2,LXP1,J),HX,U)
          SUY(L) = SPLN0(C(1,LX,I),C(1,LXP1,I),C(2,LX,I),
     *        C(2,LXP1,I),HX,U)
          SV(L) = SPLN0(C(1,LXPL,KM1),C(1,LXPL,KP1),C(1,LXPL,K),
     *        C(1,LXPL,KP2),HY,V)
          SVX(L) = SPLN0(C(2,LXPL,KM1),C(2,LXPL,KP1),C(2,LXPL,K),
     *        C(2,LXPL,KP2),HY,V)
      ENDDO

      PDS(1) = SPLN0(SV(1),SV(2),SVX(1),SVX(2),HX,U)
      PDS(2) = SPLN1(SV(1),SV(2),SVX(1),SVX(2),HX,U)
      PDS(3) = SPLN1(SU(1),SU(2),SUY(1),SUY(2),HY,V)
      PDS(4) = SPLN1(SUX(1),SUX(2),SXY(1),SXY(2),HY,V)
      PDS(5) = SPLN2(SV(1),SV(2),SVX(1),SVX(2),HX,U)
      PDS(6) = SPLN2(SU(1),SU(2),SUY(1),SUY(2),HY,V)

      IER = 0

      CONTAINS

        REAL FUNCTION SPLN0(S0,SH,SP0,SPH,H,D)
          REAL,INTENT(IN):: S0,SH,SP0,SPH,H,D
          SPLN0 = S0+D*(H*SP0+D*(3.*(SH-S0)-
     *            (SPH+2.*SP0)*H+D*(2.*(S0-SH)+(SPH+SP0)*H)))
        END FUNCTION SPLN0

        REAL FUNCTION SPLN1(S0,SH,SP0,SPH,H,D)
           REAL,INTENT(IN):: S0,SH,SP0,SPH,H,D
           SPLN1 = SP0+D*(6.*(SH-S0)/H-2.*
     *            (SPH+2.*SP0)+3.*D*(2.*(S0-SH)/H+(SPH+SP0)))
        END FUNCTION SPLN1

        REAL FUNCTION SPLN2(S0,SH,SP0,SPH,H,D)
          REAL,INTENT(IN):: S0,SH,SP0,SPH,H,D
          SPLN2 = 6.*(SH-S0)/(H*H)-2.*
     *            (SPH+2.*SP0)/H+D*(2.*(S0-SH)/(H*H)+(SPH+SP0)/H)*6.
        END FUNCTION SPLN2

      END SUBROUTINE DBCEVL
C
C-----------------------------------------------------------------------
C
C     SIMPLIFIED VERSION, SPLINE VALUES ONLY
C
C-----------------------------------------------------------------------
C
      SUBROUTINE BCSEVL_S (X,NX,Y,NY,C,IC,XL,YL,ZL,IER)
      IMPLICIT                NONE
C
C                             SPECIFICATIONS FOR ARGUMENTS
C
      INTEGER, INTENT(IN )::  NX,NY,IC
      INTEGER, INTENT(OUT)::  IER
      REAL,    INTENT(IN )::  X(NX),Y(NY),C(2,IC,2*NY),XL,YL
      REAL,    INTENT(OUT)::  ZL
C
C                             SPECIFICATIONS FOR LOCAL VARIABLES
C
      INTEGER                 K,KM1,KP1,KP2,LX,LY,LXP1
      REAL                    HX,HY,U,V

      IER = 1
      IF(XL < X(1) .OR. XL > X(NX) .OR.
     *   YL < Y(1) .OR. YL > Y(NY))RETURN

      IER = 0
      DO K=1,NX
        IF(XL < X(K)) EXIT
      ENDDO
      LX = K-1
      IF(LX == NX )LX = NX-1

      DO K=1,NY
        IF(YL < Y(K)) EXIT
      ENDDO
      LY = K-1
      IF(LY == NY )LY = NY-1

      LXP1 = LX+1
      HX = X(LXP1)-X(LX)
      HY = Y(LY+1)-Y(LY)
      U = (XL-X(LX))/HX
      V = (YL-Y(LY))/HY
      K = 2*LY
      KP1 = K+1
      KP2 = K+2
      KM1 = K-1

      ZL = SPLN0(
     *      SPLN0(C(1,LX,KM1),C(1,LX,KP1),C(1,LX,K),
     *            C(1,LX,KP2),HY,V),
     *      SPLN0(C(1,LXP1,KM1),C(1,LXP1,KP1),C(1,LXP1,K),
     *            C(1,LXP1,KP2),HY,V),
     *      SPLN0(C(2,LX,KM1),C(2,LX,KP1),C(2,LX,K),
     *            C(2,LX,KP2),HY,V),
     *      SPLN0(C(2,LXP1,KM1),C(2,LXP1,KP1),C(2,LXP1,K),
     *            C(2,LXP1,KP2),HY,V),
     *      HX,U)


      CONTAINS

        REAL FUNCTION SPLN0(S0,SH,SP0,SPH,H,D)
          REAL,INTENT(IN):: S0,SH,SP0,SPH,H,D
          SPLN0 = S0+D*(H*SP0+D*(3.*(SH-S0)-
     *           (SPH+2.*SP0)*H+D*(2.*(S0-SH)+(SPH+SP0)*H)))
        END FUNCTION SPLN0

      END SUBROUTINE BCSEVL_S
C
C-----------------------------------------------------------------------
C
C     Could optionally cut extremums
C
C-----------------------------------------------------------------------
C

      SUBROUTINE BCSEVL_EX(X,NX,Y,NY,C,IC,XL,YL,ZL,IER,CUT_EXTR)
      IMPLICIT                NONE
C
C                            SPECIFICATIONS FOR ARGUMENTS
C
      INTEGER, INTENT(IN )::    NX,NY,IC,CUT_EXTR
      REAL,    INTENT(IN )::    X(NX),Y(NY),C(2,IC,2*NY),XL,YL
      INTEGER, INTENT(OUT)::    IER
      REAL,    INTENT(OUT)::    ZL
C
C                            SPECIFICATIONS FOR LOCAL VARIABLES
C
      INTEGER                 K,KM1,KP1,KP2,LX,LY,LXP1
      REAL                    U,V,HX,HY,CMX,CMN,C1,C2,C3,C4,ZZ

      IER = 1
      IF(XL < X(1) .OR. XL > X(NX) .OR.
     *   YL < Y(1) .OR. YL > Y(NY))RETURN

      IER = 0

      DO K=1,NX
        IF(XL < X(K)) EXIT
      ENDDO
      LX = K-1
      IF(LX == NX )LX = NX-1

      DO K=1,NY
        IF(YL < Y(K)) EXIT
      ENDDO
      LY = K-1
      IF(LY == NY )LY = NY-1

      LXP1 = LX+1
      HX = X(LXP1)-X(LX)
      HY = Y(LY+1)-Y(LY)
      U = (XL-X(LX))/HX
      V = (YL-Y(LY))/HY
      K = 2*LY
      KP1 = K+1
      KP2 = K+2
      KM1 = K-1
      C1=C(1,LX,KM1)
      C2=C(1,LXP1,KM1)
      C3=C(1,LX,KP1)
      C4=C(1,LXP1,KP1)

      ZZ = SPLN0(
     *      SPLN0(C1,C3,C(1,LX,K),C(1,LX,KP2),HY,V),
     *      SPLN0(C2,C4,C(1,LXP1,K),C(1,LXP1,KP2),HY,V),
     *      SPLN0(C(2,LX,KM1),C(2,LX,KP1),C(2,LX,K),
     *            C(2,LX,KP2),HY,V),
     *      SPLN0(C(2,LXP1,KM1),C(2,LXP1,KP1),C(2,LXP1,K),
     *            C(2,LXP1,KP2),HY,V),
     *      HX,U)
      IF(CUT_EXTR .NE. 0)THEN
        CMN=MIN(C1,C2,C3,C4)
        CMX=MAX(C1,C2,C3,C4)
        IF(ZZ < CMN)THEN
          ZZ = CMN
        ELSEIF(ZZ > CMX)THEN
          ZZ = CMX
        ENDIF
      ENDIF
      ZL = ZZ

      CONTAINS

        REAL FUNCTION SPLN0(S0,SH,SP0,SPH,H,D)
          REAL,INTENT(IN):: S0,SH,SP0,SPH,H,D
          SPLN0 = S0+D*(H*SP0+D*(3.*(SH-S0)-
     *            (SPH+2.*SP0)*H+D*(2.*(S0-SH)+(SPH+SP0)*H)))
        END FUNCTION SPLN0

      END SUBROUTINE BCSEVL_EX
C
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C
C    VERSION OPTIMIZED FOR VECTOR PROCESSORS (FUJITSU VPP).
C    FOR SCALAR PROCESSING PLEASE USE MORE CLEAR
C    DBCEVL(INCLUDE PARTIAL DERIVATIVES ESTIMATIONS)
C    OR BCSEVL_S
C
C    S. VARLAMOV
c    Last revised: DECEMBER 2003
C
C-----------------------------------------------------------------------
C
C                XL      - VECTOR OF LENGTH NXL. (INPUT) XL MUST BE
C                           ALL INSIDE OF X(1)<=XL(I)<=X(NX) FOR
C                           I=1,...,NXL-1.
C                NXL     - NUMBER OF ELEMENTS IN XL. (INPUT)
C                YL      - VECTOR OF LENGTH NYL. (INPUT) YL MUST BE
C                           ALL INSIDE OF Y(1)<=YL(J)<=Y(NY) FOR
C                           J=1,...,NYL-1.
C                NYL     - NUMBER OF ELEMENTS IN YL. (INPUT)
C                         NOTE - THE COORDINATE PAIRS (XL(I),YL(J)), FOR
C                           I=1,...,NXL AND J=1,...,NYL, GIVE THE POINTS
C                           WHERE THE THE BICUBIC SPLINE IS EVALUATED.
C                ZL     - ARRAY OF BICUBIC SPLINE VALUES, S(X,Y),
C                           EVALUATED AT X=XL(I) AND Y=YL(I). (OUTPUT)
C                           DIMENSION ZL(NXD,NYL), ZL(I,J) = S(XL(I),YL(J))
C                NXD    - FIRST DIMENSION OF ARRAY ZL EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT
C                           (INPUT).  NXD MUST BE .GE. NXL.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE BCSEVL_V(X,NX,Y,NY,C,IC,XL,NXL,YL,NYL,ZL,NXD,IER)
      IMPLICIT                NONE
C                             ARGUMENTS
      INTEGER, INTENT(IN )::  NX,NY,IC,NXL,NYL,NXD
      INTEGER, INTENT(OUT)::  IER
      REAL,    INTENT(IN )::  X(NX),Y(NY),C(2,IC,2*NY),XL(NXL),YL(NYL)
      REAL,    INTENT(OUT)::  ZL(NXD,NYL)
C                             LOCAL
      INTEGER                 I,J,K,L,KM1,KP1,KP2,LX,LY,LXP1
      REAL                    HX,HY,U,V
      REAL,PARAMETER::        EPS = 1.0E-7
cC
cC     STATEMENT-FUNCTION::
cC
c      REAL SPLN0,S0,SH,SP0,SPH,H,D
c
c      SPLN0(S0,SH,SP0,SPH,H,D) = S0+D*(H*SP0+D*(3.*(SH-S0)-
c     *            (SPH+2.*SP0)*H+D*(2.*(S0-SH)+(SPH+SP0)*H)))
c
      IER = 0
      IF( X(NX)-X(1) < EPS)                    IER = IER + 1
      IF( Y(NY)-Y(1) < EPS)                    IER = IER + 2
      IF( XL(1) < X(1) .OR. XL(1) > X(NX))     IER = IER + 4
      IF( XL(NXL) < X(1) .OR. XL(NXL) > X(NX)) IER = IER + 8
      IF(NXL > 3)THEN
        I = NXL/2
        IF(XL(I) < X(1) .OR. XL(I) > X(NX))    IER = IER + 16
      ENDIF
      IF( YL(1) < Y(1) .OR. YL(1) > Y(NY))     IER = IER + 32
      IF( YL(NYL) < Y(1) .OR. YL(NYL) > Y(NY)) IER = IER + 64
      IF(NYL > 3)THEN
        I = NYL/2
        IF(YL(I) < Y(1) .OR. YL(I) > Y(NY))    IER = IER + 128
      ENDIF
      IF( IER > 0) RETURN

      DO J = 1,NYL
        LY = 1
        DO L = 1,NY-1
          IF(YL(J)>=Y(L) .AND. YL(J)<=Y(L+1))THEN
            LY = L
          ENDIF
        ENDDO

        HY = Y(LY+1)-Y(LY)
        V = (YL(J)-Y(LY))/HY
        K = 2*LY
        KP1 = K+1
        KP2 = K+2
        KM1 = K-1

        DO I = 1,NXL

          LX = 1
          DO L = 1,NX-1
            IF(XL(I)>=X(L) .AND. XL(I)<=X(L+1))THEN
              LX = L
            ENDIF
          ENDDO

          LXP1 = LX + 1
          HX = X(LXP1) - X(LX)
          U = (XL(I) - X(LX))/HX

          ZL(I,J) = SPLN0(
     *      SPLN0(C(1,LX,KM1),C(1,LX,KP1),C(1,LX,K),
     *            C(1,LX,KP2),HY,V),
     *      SPLN0(C(1,LXP1,KM1),C(1,LXP1,KP1),C(1,LXP1,K),
     *            C(1,LXP1,KP2),HY,V),
     *      SPLN0(C(2,LX,KM1),C(2,LX,KP1),C(2,LX,K),
     *            C(2,LX,KP2),HY,V),
     *      SPLN0(C(2,LXP1,KM1),C(2,LXP1,KP1),C(2,LXP1,K),
     *            C(2,LXP1,KP2),HY,V),
     *      HX,U)
        ENDDO
      ENDDO

      CONTAINS

        REAL FUNCTION SPLN0(S0,SH,SP0,SPH,H,D)
          REAL,INTENT(IN):: S0,SH,SP0,SPH,H,D
          SPLN0 = S0+D*(H*SP0+D*(3.*(SH-S0)-
     *            (SPH+2.*SP0)*H+D*(2.*(S0-SH)+(SPH+SP0)*H)))
        END FUNCTION SPLN0


      END SUBROUTINE BCSEVL_V
C
C-----------------------------------------------------------------------
C
C     VERSION WHEN MASK IS USED
C
C
      SUBROUTINE BCSEVL_VM(X,NX,Y,NY,C,IC,XL,NXL,YL,NYL,ZL,NXD,MASK,
     *                     UNDEF,IER)
C
      IMPLICIT                NONE
C                             ARGUMENTS
      INTEGER, INTENT(IN )::  NX,NY,IC,NXL,NYL,NXD
      INTEGER, INTENT(OUT)::  IER
      REAL,    INTENT(IN )::  X(NX),Y(NY),C(2,IC,2*NY),XL(NXL),YL(NYL),
     *                        MASK(NXD,NYL), UNDEF
      REAL,    INTENT(OUT)::  ZL(NXD,NYL)
C                             LOCAL
      INTEGER                 I,J,K,L,KM1,KP1,KP2,LX,LY,LXP1
      REAL                    HX,HY,U,V
      REAL,PARAMETER::        EPS = 1.0E-7
cC
cC     STATEMENT-FUNCTION::
cC
c      REAL SPLN0,S0,SH,SP0,SPH,H,D
c
c      SPLN0(S0,SH,SP0,SPH,H,D) = S0+D*(H*SP0+D*(3.*(SH-S0)-
c     *            (SPH+2.*SP0)*H+D*(2.*(S0-SH)+(SPH+SP0)*H)))
c
      IER = 0
      IF( X(NX)-X(1) < EPS)                    IER = IER + 1
      IF( Y(NY)-Y(1) < EPS)                    IER = IER + 2
      IF( XL(1) < X(1) .OR. XL(1) > X(NX))     IER = IER + 4
      IF( XL(NXL) < X(1) .OR. XL(NXL) > X(NX)) IER = IER + 8
      IF(NXL > 3)THEN
        I = NXL/2
        IF(XL(I) < X(1) .OR. XL(I) > X(NX))    IER = IER + 16
      ENDIF
      IF( YL(1) < Y(1) .OR. YL(1) > Y(NY))     IER = IER + 32
      IF( YL(NYL) < Y(1) .OR. YL(NYL) > Y(NY)) IER = IER + 64
      IF(NYL > 3)THEN
        I = NYL/2
        IF(YL(I) < Y(1) .OR. YL(I) > Y(NY))    IER = IER + 128
      ENDIF
      IF( IER > 0) RETURN

      DO J=1,NYL
        LY = 1
        DO L = 1,NY-1
          IF(YL(J)>=Y(L) .AND. YL(J)<=Y(L+1))THEN
            LY=L
          ENDIF
        ENDDO

        HY = Y(LY+1)-Y(LY)
        V = (YL(J)-Y(LY))/HY
        K = 2*LY
        KP1 = K+1
        KP2 = K+2
        KM1 = K-1

        DO I=1,NXL
          IF(MASK(I,J)>0.)THEN
            LX = 1
            DO L = 1,NX-1
              IF(XL(I)>=X(L) .AND. XL(I)<=X(L+1))THEN
                LX=L
              ENDIF
            ENDDO

            LXP1 = LX+1
            HX = X(LXP1)-X(LX)
            U = (XL(I)-X(LX))/HX

            ZL(I,J) = SPLN0(
     *      SPLN0(C(1,LX,KM1),C(1,LX,KP1),C(1,LX,K),
     *            C(1,LX,KP2),HY,V),
     *      SPLN0(C(1,LXP1,KM1),C(1,LXP1,KP1),C(1,LXP1,K),
     *            C(1,LXP1,KP2),HY,V),
     *      SPLN0(C(2,LX,KM1),C(2,LX,KP1),C(2,LX,K),
     *            C(2,LX,KP2),HY,V),
     *      SPLN0(C(2,LXP1,KM1),C(2,LXP1,KP1),C(2,LXP1,K),
     *            C(2,LXP1,KP2),HY,V),
     *      HX,U)
          ELSE
            ZL(I,J) = UNDEF
          ENDIF
        ENDDO
      ENDDO

      CONTAINS

        REAL FUNCTION SPLN0(S0,SH,SP0,SPH,H,D)
          REAL,INTENT(IN):: S0,SH,SP0,SPH,H,D
          SPLN0 = S0+D*(H*SP0+D*(3.*(SH-S0)-
     *            (SPH+2.*SP0)*H+D*(2.*(S0-SH)+(SPH+SP0)*H)))
        END FUNCTION SPLN0

      END SUBROUTINE BCSEVL_VM
C
C-----------------------------------------------------------------------
C
C     Vectorized with [unconditional] limitation of extremums
C
C-----------------------------------------------------------------------
C
      SUBROUTINE BCSEVL_VE(X,NX,Y,NY,C,IC,XL,NXL,YL,NYL,ZL,NXD,IER)
      IMPLICIT                NONE
C                             ARGUMENTS
      INTEGER, INTENT(IN )::  NX,NY,IC,NXL,NYL,NXD
      INTEGER, INTENT(OUT)::  IER
      REAL,    INTENT(IN )::  X(NX),Y(NY),C(2,IC,2*NY),XL(NXL),YL(NYL)
      REAL,    INTENT(OUT)::  ZL(NXD,NYL)
C                             LOCAL
      INTEGER                 I,J,K,L,KM1,KP1,KP2,LX,LY,LXP1
      REAL                    HX,HY,U,V,C1,C2,C3,C4,ZZ,CMN,CMX
      REAL,PARAMETER::        EPS = 1.0E-7
cC
cC     STATEMENT-FUNCTION::
cC
c      REAL SPLN0,S0,SH,SP0,SPH,H,D
c
c      SPLN0(S0,SH,SP0,SPH,H,D) = S0+D*(H*SP0+D*(3.*(SH-S0)-
c     *            (SPH+2.*SP0)*H+D*(2.*(S0-SH)+(SPH+SP0)*H)))
c
      IER = 0
      IF( X(NX)-X(1) < EPS)                    IER = IER + 1
      IF( Y(NY)-Y(1) < EPS)                    IER = IER + 2
      IF( XL(1) < X(1) .OR. XL(1) > X(NX))     IER = IER + 4
      IF( XL(NXL) < X(1) .OR. XL(NXL) > X(NX)) IER = IER + 8
      IF(NXL > 3)THEN
        I = NXL/2
        IF(XL(I) < X(1) .OR. XL(I) > X(NX))    IER = IER + 16
      ENDIF
      IF( YL(1) < Y(1) .OR. YL(1) > Y(NY))     IER = IER + 32
      IF( YL(NYL) < Y(1) .OR. YL(NYL) > Y(NY)) IER = IER + 64
      IF(NYL > 3)THEN
        I = NYL/2
        IF(YL(I) < Y(1) .OR. YL(I) > Y(NY))    IER = IER + 128
      ENDIF
      IF( IER > 0) RETURN


      DO J = 1,NYL
        LY = 1
        DO L = 1,NY-1
          IF(YL(J)>=Y(L) .AND. YL(J)<=Y(L+1))THEN
            LY = L
          ENDIF
        ENDDO

        HY = Y(LY+1)-Y(LY)
        V = (YL(J)-Y(LY))/HY
        K = 2*LY
        KP1 = K+1
        KP2 = K+2
        KM1 = K-1

        DO I = 1,NXL

          LX = 1
          DO L = 1,NX-1
            IF(XL(I)>=X(L) .AND. XL(I)<=X(L+1))THEN
              LX = L
            ENDIF
          ENDDO

          LXP1 = LX + 1
          HX = X(LXP1) - X(LX)
          U = (XL(I) - X(LX))/HX

          C1=C(1,LX,KM1)
          C2=C(1,LXP1,KM1)
          C3=C(1,LX,KP1)
          C4=C(1,LXP1,KP1)

          ZZ = SPLN0(
     *      SPLN0(C1,C3,C(1,LX,K),
     *            C(1,LX,KP2),HY,V),
     *      SPLN0(C2,C4,C(1,LXP1,K),
     *            C(1,LXP1,KP2),HY,V),
     *      SPLN0(C(2,LX,KM1),C(2,LX,KP1),C(2,LX,K),
     *            C(2,LX,KP2),HY,V),
     *      SPLN0(C(2,LXP1,KM1),C(2,LXP1,KP1),C(2,LXP1,K),
     *            C(2,LXP1,KP2),HY,V),
     *      HX,U)

          CMN=MIN(C1,C2,C3,C4)
          CMX=MAX(C1,C2,C3,C4)
          IF(ZZ < CMN)THEN
            ZZ = CMN
          ELSEIF(ZZ > CMX)THEN
            ZZ = CMX
          ENDIF
          ZL(I,J) = ZZ

        ENDDO
      ENDDO

      CONTAINS

        REAL FUNCTION SPLN0(S0,SH,SP0,SPH,H,D)
          REAL,INTENT(IN):: S0,SH,SP0,SPH,H,D
          SPLN0 = S0+D*(H*SP0+D*(3.*(SH-S0)-
     *            (SPH+2.*SP0)*H+D*(2.*(S0-SH)+(SPH+SP0)*H)))
        END FUNCTION SPLN0


      END SUBROUTINE BCSEVL_VE
C
C-----------------------------------------------------------------------
C
C     VERSION WHEN MASK IS USED
C
C
      SUBROUTINE BCSEVL_VME(X,NX,Y,NY,C,IC,XL,NXL,YL,NYL,ZL,NXD,MASK,
     *                    UNDEF,IER)
      IMPLICIT                NONE
C                             ARGUMENTS
      INTEGER, INTENT(IN )::  NX,NY,IC,NXL,NYL,NXD
      INTEGER, INTENT(OUT)::  IER
      REAL,    INTENT(IN )::  X(NX),Y(NY),C(2,IC,2*NY),XL(NXL),YL(NYL),
     *                        MASK(NXD,NYL), UNDEF
      REAL,    INTENT(OUT)::  ZL(NXD,NYL)
C                             LOCAL
      INTEGER                 I,J,K,L,KM1,KP1,KP2,LX,LY,LXP1
      REAL                    HX,HY,U,V,C1,C2,C3,C4,ZZ,CMN,CMX
      REAL,PARAMETER::        EPS = 1.0E-7
cC
cC     STATEMENT-FUNCTION::
cC
c      REAL SPLN0,S0,SH,SP0,SPH,H,D
c
c      SPLN0(S0,SH,SP0,SPH,H,D) = S0+D*(H*SP0+D*(3.*(SH-S0)-
c     *            (SPH+2.*SP0)*H+D*(2.*(S0-SH)+(SPH+SP0)*H)))
c
      IER = 0
      IF( X(NX)-X(1) < EPS)                    IER = IER + 1
      IF( Y(NY)-Y(1) < EPS)                    IER = IER + 2
      IF( XL(1) < X(1) .OR. XL(1) > X(NX))     IER = IER + 4
      IF( XL(NXL) < X(1) .OR. XL(NXL) > X(NX)) IER = IER + 8
      IF(NXL > 3)THEN
        I = NXL/2
        IF(XL(I) < X(1) .OR. XL(I) > X(NX))    IER = IER + 16
      ENDIF
      IF( YL(1) < Y(1) .OR. YL(1) > Y(NY))     IER = IER + 32
      IF( YL(NYL) < Y(1) .OR. YL(NYL) > Y(NY)) IER = IER + 64
      IF(NYL > 3)THEN
        I = NYL/2
        IF(YL(I) < Y(1) .OR. YL(I) > Y(NY))    IER = IER + 128
      ENDIF
      IF( IER > 0) RETURN

      DO J=1,NYL
        LY = 1
        DO L = 1,NY-1
          IF(YL(J)>=Y(L) .AND. YL(J)<=Y(L+1))THEN
            LY=L
          ENDIF
        ENDDO

        HY = Y(LY+1)-Y(LY)
        V = (YL(J)-Y(LY))/HY
        K = 2*LY
        KP1 = K+1
        KP2 = K+2
        KM1 = K-1

        DO I=1,NXL
          IF(MASK(I,J)>0.)THEN
            LX = 1
            DO L = 1,NX-1
              IF(XL(I)>=X(L) .AND. XL(I)<=X(L+1))THEN
                LX=L
              ENDIF
            ENDDO

            LXP1 = LX + 1
            HX = X(LXP1) - X(LX)
            U = (XL(I) - X(LX))/HX

            C1=C(1,LX,KM1)
            C2=C(1,LXP1,KM1)
            C3=C(1,LX,KP1)
            C4=C(1,LXP1,KP1)


            ZZ = SPLN0(
     *      SPLN0(C1,C3,C(1,LX,K),
     *            C(1,LX,KP2),HY,V),
     *      SPLN0(C2,C4,C(1,LXP1,K),
     *            C(1,LXP1,KP2),HY,V),
     *      SPLN0(C(2,LX,KM1),C(2,LX,KP1),C(2,LX,K),
     *            C(2,LX,KP2),HY,V),
     *      SPLN0(C(2,LXP1,KM1),C(2,LXP1,KP1),C(2,LXP1,K),
     *            C(2,LXP1,KP2),HY,V),
     *      HX,U)

            CMN=MIN(C1,C2,C3,C4)
            CMX=MAX(C1,C2,C3,C4)
            IF(ZZ < CMN)THEN
              ZZ = CMN
            ELSEIF(ZZ > CMX)THEN
              ZZ = CMX
            ENDIF
            ZL(I,J) = ZZ
          ELSE
            ZL(I,J) = UNDEF
          ENDIF
        ENDDO
      ENDDO

      CONTAINS

        REAL FUNCTION SPLN0(S0,SH,SP0,SPH,H,D)
          REAL,INTENT(IN):: S0,SH,SP0,SPH,H,D
          SPLN0 = S0+D*(H*SP0+D*(3.*(SH-S0)-
     *            (SPH+2.*SP0)*H+D*(2.*(S0-SH)+(SPH+SP0)*H)))
        END FUNCTION SPLN0

      END SUBROUTINE BCSEVL_VME

C-----------------------------------------------------------------------
C
C    BICUBIC SPLINE EVALUATOR WITH CUTTING OFF EXTREMUMS
C
C    VERSION OPTIMIZED FOR VECTOR PROCESSORS (FUJITSU VPP).
C-----------------------------------------------------------------------
C      CALL BCSEVL_PVEXP
C     *  (X,NX,Y,NY,C,IC,PT,NP,IER,GLOBX,UNDEF,CUT_EXTR,MASK)
C-----------------------------------------------------------------------
C   ARGUMENTS    PT      - TPOINT TYPE VECTOR OF LENGTH NP. (INPUT)
C                NP      - NUMBER OF ELEMENTS IN PT. (INPUT)
C                Z       - ARRAY OF BICUBIC SPLINE VALUES, S(PT),
C                            EVALUATED AT PT (PT%Z IS IGNORED) (OUTPUT)
C                            DIMENSION Z(PT), Z(PT) = UNDEF
C                            IF POINT IS OUT OF X,Y DOMAIN
C (OPTIONAL)     GLOBX   - NON ZERO IF X IS A LONGITUDE ON GLOBE, IN THIS
C                           CASE ALL X VALUES MUST BE IN DEGREES
C                           (PERIOD 360 ASSUMED)
C (OPTIONAL)     UNDEF   - VALUE ASSIGNED TO POINTS OUTSIDE OF
C                           INTERPOLATION DOMAIN
C (OPTIONAL)     CUT_EXTR - IF NON ZERO, EXTREMUMS EXCEEDING 4 CLOSEST
C                           POINT VALUES ARE CUTOFF
C (OPTIONAL)     MASK     - LOGICAL, DIMENSION(NP). IF PROVIDED, ESTIMATIONS
C                           ARE GOING ONLY WHERE MASK VALUES ARE TRUE
C
C                IER     - ERROR PARAMETER. (OUTPUT)
C
C-----------------------------------------------------------------------
      SUBROUTINE BCSEVL_PVEXP
     * (X,NX,Y,NY,C,IC,PT,NP,IER,GLOBX,UNDEF,CUT_EXTR,MASK,VALUE)
C-----------------------------------------------------------------------
c      USE COM_SERV, ONLY: TPOINT
      IMPLICIT NONE

      TYPE TPOINT
          REAL X
          REAL Y
          REAL Z
      END TYPE
C-----------------------------------------------------------------------
C     ARGUMENTS
C-----------------------------------------------------------------------
      INTEGER,     INTENT(IN):: NX, NY, IC, NP
      REAL,        INTENT(IN):: X(NX), Y(NY), C(2,IC,2*NY)
      TYPE(TPOINT),INTENT(INOUT):: PT(NP)
      INTEGER, INTENT(OUT):: IER
      INTEGER, OPTIONAL, INTENT(IN):: GLOBX, CUT_EXTR
      LOGICAL, OPTIONAL, INTENT(IN):: MASK(NP)
      REAL,    OPTIONAL, INTENT(IN):: UNDEF
      REAL,    OPTIONAL, INTENT(OUT):: VALUE(NP)
C-----------------------------------------------------------------------
C     LOCAL
C-----------------------------------------------------------------------
      INTEGER          K,KM1,KP1,KP2,LX,LY,LXP1,IP,GLOBXL,CUT_EXTRL
      REAL             HX,HY,U,V,ZZ,C1,C2,C3,C4,CMN,CMX,UNDEFL,UNDEFCRIT
      REAL,PARAMETER:: EPS = 1.0E-7
      LOGICAL          UNDEFINE
      REAL             XX(NP),YY(NP),LVALUE(NP)
      INTEGER          INDI(NP),INDJ(NP)

      IER = 0
      IF(X(NX)-X(1) < EPS)IER = IER + 1
      IF(Y(NY)-Y(1) < EPS)IER = IER + 2
      IF(IER > 0) RETURN
C-----------------------------------------------------------------------
C     DEFAULT VALUES
C-----------------------------------------------------------------------
      GLOBXL = 0
      CUT_EXTRL = 0

      IF(PRESENT(GLOBX))GLOBXL=GLOBX
      IF(PRESENT(CUT_EXTR))CUT_EXTRL=CUT_EXTR
      IF(PRESENT(UNDEF))THEN
        UNDEFL=UNDEF
        UNDEFINE = .TRUE.
        UNDEFCRIT = 0.  !NOT USED
      ELSE
        UNDEFL = HUGE(UNDEFL)
        UNDEFCRIT = UNDEFL-1.
        UNDEFINE = .FALSE.
      ENDIF

      IF(GLOBXL == 0)THEN
        XX = PT%X
      ELSE
        XX = MOD(PT%X-X(1),360.)
        WHERE(XX < 0.) XX = XX + 360.
        XX = XX + X(1)
      ENDIF
      WHERE(ABS(XX-X(1))<EPS)XX=X(1)
      WHERE(ABS(XX-X(NX))<EPS)XX=X(NX)
      YY = PT%Y
      WHERE(ABS(YY-Y(1))<EPS)YY=Y(1)
      WHERE(ABS(YY-Y(NY))<EPS)YY=Y(NY)
C-----------------------------------------------------------------------
C     IF NP >> NX, NY NEXT WAY LOOKS MOST EFFICIENTLY VECTORIZED ON VPP:
C-----------------------------------------------------------------------
      INDI = 0
      INDJ = 0
      DO LX = 1,NX-1
        WHERE(INDI == 0 .AND. XX >= X(LX) .AND. XX <= X(LX+1))
     *        INDI = LX
      ENDDO
      DO LY = 1,NY-1
        WHERE(INDJ == 0 .AND. YY >= Y(LY) .AND. YY <= Y(LY+1))
     *        INDJ = LY
      ENDDO
      IF(PRESENT(MASK))THEN
        WHERE(.NOT. MASK)INDI = 0
      ENDIF

      DO IP = 1,NP
        LY = INDJ(IP)
        LX = INDI(IP)
        IF(LX*LY == 0)THEN
          LVALUE(IP) = UNDEFL
        ELSE
          HY = Y(LY+1)-Y(LY)
          V = (YY(IP)-Y(LY))/HY
          K = 2*LY
          KP1 = K+1
          KP2 = K+2
          KM1 = K-1
          LXP1 = LX + 1
          HX = X(LXP1) - X(LX)
          U = (XX(IP) - X(LX))/HX
          C1=C(1,LX,KM1)
          C2=C(1,LXP1,KM1)
          C3=C(1,LX,KP1)
          C4=C(1,LXP1,KP1)

          ZZ = SPLN0(
     *      SPLN0(C1,C3,C(1,LX,K),C(1,LX,KP2),HY,V),
     *      SPLN0(C2,C4,C(1,LXP1,K),C(1,LXP1,KP2),HY,V),
     *      SPLN0(C(2,LX,KM1),C(2,LX,KP1),C(2,LX,K),
     *            C(2,LX,KP2),HY,V),
     *      SPLN0(C(2,LXP1,KM1),C(2,LXP1,KP1),C(2,LXP1,K),
     *            C(2,LXP1,KP2),HY,V),
     *      HX,U)
          IF(CUT_EXTRL .NE. 0)THEN
            CMN=MIN(C1,C2,C3,C4)
            CMX=MAX(C1,C2,C3,C4)
            IF(ZZ < CMN)THEN
              ZZ = CMN
            ELSEIF(ZZ > CMX)THEN
              ZZ = CMX
            ENDIF
          ENDIF
          LVALUE(IP) = ZZ
        ENDIF
      ENDDO

      IF(PRESENT(VALUE))THEN
        IF(UNDEFINE)THEN
          VALUE = LVALUE
        ELSE
C-----------------------------------------------------------------------
C         KEEP MASKED OR OUT OF DOMAIN VALUES NOT MODIFIED
C-----------------------------------------------------------------------
          WHERE(LVALUE < UNDEFCRIT)VALUE = LVALUE
        ENDIF
      ELSE
        IF(UNDEFINE)THEN
          PT(:)%Z = LVALUE(:)
        ELSE
C-----------------------------------------------------------------------
C         KEEP MASKED OR OUT OF DOMAIN VALUES NOT MODIFIED
C-----------------------------------------------------------------------
          WHERE(LVALUE(:) < UNDEFCRIT)PT(:)%Z = LVALUE(:)
        ENDIF
      ENDIF

      CONTAINS

        REAL FUNCTION SPLN0(S0,SH,SP0,SPH,H,D)
          REAL,INTENT(IN):: S0,SH,SP0,SPH,H,D
          SPLN0 = S0+D*(H*SP0+D*(3.*(SH-S0)-
     *            (SPH+2.*SP0)*H+D*(2.*(S0-SH)+(SPH+SP0)*H)))
        END FUNCTION SPLN0

      END SUBROUTINE BCSEVL_PVEXP
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C
C   BICUBIC SPLINE TWO-DIMENSIONAL COEFFICIENT CALCULATOR
C
C    RESTRUCTURIZATION OF CODE FOR FORTRAN 90
C    BASED ON IMSL ROUTINE IBCCCU
C
C-----------------------------------------------------------------------
C
C   USAGE               - CALL GET_BSPLINES (F,X,NX,Y,NY,C,IC,IER)
C
C   ARGUMENTS    F      - NX BY NY MATRIX CONTAINING THE FUNCTION
C                           VALUES. (INPUT) F(I,J) IS THE FUNCTION VALUE
C                           AT THE POINT (X(I),Y(J)) FOR I=1,...,NX AND
C                           J=1,...,NY.
C                X      - VECTOR OF LENGTH NX. (INPUT) X MUST BE
C                           ORDERED SO THAT X(I) .LT. X(I+1) FOR
C                           I=1,...,NX-1.
C                NX     - NUMBER OF ELEMENTS IN X. (INPUT) NX MUST BE
C                           .GE. 4.
C                Y      - VECTOR OF LENGTH NY. (INPUT) Y MUST BE
C                           ORDERED SO THAT Y(J) .LT. Y(J+1) FOR
C                           J=1,...,NY-1.
C                NY     - NUMBER OF ELEMENTS IN Y. (INPUT) NY MUST BE
C                           .GE. 4.
C                         NOTE - THE COORDINATE PAIRS (X(I),Y(J)), FOR
C                           I=1,...,NX AND J=1,...,NY, GIVE THE POINTS
C                           WHERE THE FUNCTION VALUES F(I,J) ARE
C                           DEFINED.
C                C      - ARRAY OF SPLINE COEFFICIENTS. (OUTPUT)
C                           C IS OF DIMENSION 2 BY NX BY 2 BY NY.
C                           AT THE POINT (X(I),Y(J))
C                             C(1,I,1,J) = S
C                             C(2,I,1,J) = DS/DX
C                             C(1,I,2,J) = DS/DY
C                             C(2,I,2,J) = D(DS/DX)/DY
C                           WHERE S(X,Y) IS THE SPLINE APPROXIMATION.
C                           (NOTE - C IS TREATED INTERNALLY AS A
C                             2 BY NX BY 2*NY ARRAY BECAUSE CERTAIN
C                             ENVIRONMENTS DO NOT PERMIT QUADRUPLY-
C                             DIMENSIONED ARRAYS.  IN THESE
C                             ENVIRONMENTS THE CALLING PROGRAM MAY
C                             DIMENSION C IN THE SAME MANNER.)
C                IC     - ROW DIMENSION OF MATRIX F AND SECOND
C                           DIMENSION OF ARRAY C EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT.
C                           (INPUT). IC MUST BE .GE. NX.
C                IER    - ERROR PARAMETER. (OUTPUT)
C-----------------------------------------------------------------------
C
      SUBROUTINE GET_BSPLINES (F,X,NX,Y,NY,C,IC,IER)
C
C                             SPECIFICATIONS FOR ARGUMENTS
C
      IMPLICIT                NONE
      INTEGER, INTENT(IN )::    NX,NY,IC
      INTEGER, INTENT(OUT)::    IER
      REAL,    INTENT(IN )::    F(IC,NY),X(NX),Y(NY)
      REAL,    INTENT(OUT)::    C(2,IC,2*NY)
!
!     automatic working array
!
      REAL     WK(2*NX*NY)

      INTERFACE
        SUBROUTINE GET_SPLINES (TAU,GTAU,N,M,VS,IC1,IC2,IER)
          INTEGER, INTENT(IN)::   N,M,IC1,IC2
          INTEGER, INTENT(OUT)::  IER
          REAL, INTENT(IN)::      TAU(N),GTAU(IC1,M)
          REAL, INTENT(OUT)::     VS(IC2,2,N)
        END SUBROUTINE
      END INTERFACE

      IER = 0
      IF(IC < NX)IER = IER + 1
      IF(NX < 4) IER = IER + 2
      IF(NY < 4) IER = IER + 4
      IF(IER > 0)RETURN

      CALL GET_SPLINES(X,F,NX,NY,WK,IC,NY,IER)
      IF(IER == 0)THEN
        CALL GET_SPLINES(Y,WK,NY,2*NX,C,NY,2*IC,IER)
      ENDIF

      END SUBROUTINE GET_BSPLINES
C-----------------------------------------------------------------------
      SUBROUTINE GET_BSPLINES_S (F,X,NX,Y,NY,C,IC,IER)
C
C                             SPECIFICATIONS FOR ARGUMENTS
C
      IMPLICIT                NONE
      INTEGER, INTENT(IN )::    NX,NY,IC
      INTEGER, INTENT(OUT)::    IER
      REAL,    INTENT(IN )::    F(IC,NY),X(NX),Y(NY)
      REAL,    INTENT(OUT)::    C(2,IC,2*NY)
!
!     automatic working array
!
      REAL     WK(2*NX*NY)

      INTERFACE
        SUBROUTINE GET_SPLINES_S (TAU,GTAU,N,M,VS,IC1,IC2,IER)
          INTEGER, INTENT(IN)::   N,M,IC1,IC2
          INTEGER, INTENT(OUT)::  IER
          REAL, INTENT(IN)::      TAU(N),GTAU(IC1,M)
          REAL, INTENT(OUT)::     VS(IC2,2,N)
        END SUBROUTINE
      END INTERFACE

      IER = 0
      IF(IC < NX)IER = IER + 1
      IF(NX < 5) IER = IER + 2
      IF(NY < 5) IER = IER + 4
      IF(IER > 0)RETURN

      CALL GET_SPLINES_S(X,F,NX,NY,WK,IC,NY,IER)
      IF(IER == 0)THEN
        CALL GET_SPLINES_S(Y,WK,NY,2*NX,C,NY,2*IC,IER)
      ENDIF

      END SUBROUTINE GET_BSPLINES_S

C-----------------------------------------------------------------------
!#define BSPLINE_ERRORS
C-----------------------------------------------------------------------
      SUBROUTINE GET_SPLINES_S(TAU,GTAU,N,M,VS,IC1,IC2,IER)
C
C     FORALL IS NOT VECTORIZED ON FUJITSU VPP
C
      IMPLICIT NONE
      INTEGER, INTENT(IN)::   N,M,IC1,IC2
      INTEGER, INTENT(OUT)::  IER
      REAL, INTENT(IN)::      TAU(N),GTAU(IC1,M)
      REAL, INTENT(OUT)::     VS(IC2,2,N)
!
!     automatic working array
!
      REAL  W(N,2)

      INTEGER I,J,K,NM1,NM2,NM3,JM1,JP1,II
      REAL    C1,C2,AA,BB,CC,DD,GG,DTAU1,DTAU2,G,H,RATIO,U

      IER = 20
      IF(N < 5)RETURN
      IER = 21

      NM1 = N-1
      NM2 = N-2
      NM3 = N-3

!      W(NM1,1) = 0.
!      W(N,1) = 0.

      W(2,1) = TAU(3)-TAU(1)
      W(NM2,1) = TAU(N)-TAU(NM2)
      DTAU1 = TAU(2)-TAU(1)
      RATIO = DTAU1/W(2,1)
      W(1,2) = (RATIO-1.)*(RATIO-1.)
      W(1,1) = RATIO*(RATIO-1.)
      C1 = RATIO*(2.*RATIO-3.)

      DO J=3,NM3
        W(J,1) = TAU(J+1)-TAU(J)
      ENDDO
#ifdef BSPLINE_ERRORS
      IF(ANY(W(2:NM2,1)<= 0.0))RETURN
      IER = 22
#endif
      DO J=2,NM3
        W(J,2) = 2.*(W(J,1)+W(J+1,1))-W(J+1,1)*W(J-1,1)/W(J-1,2)
      ENDDO
#ifdef BSPLINE_ERRORS
      IF(ANY(W(1:NM3,2) == 0.0))RETURN
      IER = 23
#endif
      DTAU2 = TAU(N-1)-TAU(NM2)
      RATIO = DTAU2/W(NM2,1)
      G = -(RATIO-1.)*(RATIO-1.)/W(NM3,2)
      C2 = RATIO*(2.*RATIO-3.)

      W(NM2,2) = RATIO*(RATIO-1.)
      W(NM2,2) = G*W(NM3,1)+W(NM2,2)

#ifdef BSPLINE_ERRORS
      IF(W(NM2,2) == 0.0))RETURN
      IER = 24
#endif

!$OMP  PARALLEL DEFAULT(NONE) IF (N*M > 2000)
!$OMP& SHARED(VS,TAU,GTAU,W,M,N,NM1,NM2,NM3,DTAU1,DTAU2,C1,C2,G)
!$OMP&  PRIVATE(I,J,K,II,AA,BB,CC,DD,GG,H,U,JM1,JP1)
!$OMP DO
      DO K=1,M
        VS(K,1,1) = GTAU(1,K)
        VS(K,1,NM2) = GTAU(N,K)
      ENDDO
!$OMP END DO
!$OMP DO
!      DO I=1,NM2
!        II = I + 1
!        IF(I == 1)THEN
!          II = 1
!        ELSEIF(I == NM2)THEN
!          II = N
!        ENDIF
!        DO K=1,M
!          VS(K,1,I) = GTAU(II,K)
!        ENDDO
!      ENDDO
      DO I=2,NM3
        DO K=1,M
          VS(K,1,I) = GTAU(I+1,K)
        ENDDO
      ENDDO
!$OMP END DO
!$OMP DO
      DO I=2,NM2
        DO K=1,M
          VS(K,2,I) = (VS(K,1,I)-VS(K,1,I-1))/W(I,1)
        ENDDO
      ENDDO
!$OMP END DO
!$OMP DO
      DO K=1,M
        VS(K,2,1) = (GTAU(2,K)-GTAU(1,K))/DTAU1+VS(K,2,2)*C1
      ENDDO
!$OMP END DO

      DO I=2,NM3
        CC = 3.*W(I,1)
        DD = 3.*W(I+1,1)
        GG = -W(I+1,1)/W(I-1,2)
!$OMP DO
        DO K=1,M
          VS(K,2,I) = GG*VS(K,2,I-1)+CC*VS(K,2,I+1)+DD*VS(K,2,I)
        ENDDO
!$OMP END DO NOWAIT
      ENDDO
!$OMP BARRIER

!$OMP DO
      DO K=1,M
        CC = VS(K,2,NM2)*C2+
     *                (GTAU(NM1,K)-VS(K,1,NM3))/DTAU2
        VS(K,2,NM2) = (G*VS(K,2,NM3)+CC)/W(NM2,2)
      ENDDO
!$OMP END DO

      DO J = NM3,1,-1
!$OMP DO
        DO K=1,M
          VS(K,2,J) = (VS(K,2,J)-W(J,1)*VS(K,2,J+1))/W(J,2)
        ENDDO
!$OMP END DO NOWAIT
      ENDDO
!$OMP BARRIER

!$OMP DO
      DO J=1,2
        DO K=1,M
          VS(K,J,N) = VS(K,J,NM2)
        ENDDO
      ENDDO
!$OMP END DO

      DO I=NM1,2,-1
!$OMP DO
        DO J=1,2
          DO K=1,M
            VS(K,J,I) = VS(K,J,I-1)
          ENDDO
        ENDDO
!$OMP END DO
      ENDDO
!$OMP BARRIER

      DO J=2,NM1,NM3  !j=2,j=NM1 only
        JM1 = J-1     !1,NM2 only
        JP1 = J+1     !3,N only
C-----------------------------------------------------------------------
C       NEXT CORRECTIONS ARE REQUIRED FOR SMALL AREA N <= 4 ONLY
C-----------------------------------------------------------------------
!        IF (JM1 == 2) JM1 = 1
!        IF (JP1 == NM1) JP1 = N
C-----------------------------------------------------------------------
        H = TAU(JP1)-TAU(JM1)
        U = TAU(J)-TAU(JM1)
!$OMP DO
        DO K=1,M
          AA = VS(K,1,JM1)
          BB = VS(K,2,JM1)
          CC = (3.*(VS(K,1,JP1)-VS(K,1,JM1))/H-(VS(K,2,JP1)+
     *          2.*VS(K,2,JM1)))/H
          DD = (2.*(VS(K,1,JM1)-VS(K,1,JP1))/H+(VS(K,2,JP1)+
     *          VS(K,2,JM1)))/(H*H)
          VS(K,1,J) = AA+U*(BB+U*(CC+DD*U))
          VS(K,2,J) = BB+U*(2.*CC+3.*DD*U)
        ENDDO
!$OMP END DO
      ENDDO
!$OMP END PARALLEL

      IER = 0

      END SUBROUTINE GET_SPLINES_S
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C
      SUBROUTINE GET_SPLINES (TAU,GTAU,N,M,VS,IC1,IC2,IER)
C
C     FORALL IS NOT VECTORIZED ON FUJITSU VPP
C
      IMPLICIT NONE
      INTEGER, INTENT(IN)::   N,M,IC1,IC2
      INTEGER, INTENT(OUT)::  IER
      REAL, INTENT(IN)::      TAU(N),GTAU(IC1,M)
      REAL, INTENT(OUT)::     VS(IC2,2,N)
!
!     automatic working array
!
      REAL  W(N,2)

      INTEGER I,JM1,JP1,J,K,NM1,NM2,NM3
      REAL    AA,BB,C1,C2,CC,DD,DTAU,G,H,RATIO,U,XILIM

      NM1 = N-1
      NM2 = N-2
      NM3 = N-3
      IER = 11

      W(2,1) = TAU(3)-TAU(1)
      XILIM = TAU(1)
      IF(NM3 >= 2)THEN
        XILIM = TAU(NM2)
        DO I=3,NM2
          W(I,1) = TAU(I+1)-TAU(I)
        ENDDO
      ENDIF
      W(NM2,1) = TAU(N)-XILIM

      IF(ANY(W(2:NM2,1)<= 0.0))RETURN
      IER = 12
!$OMP PARALLEL DO SHARED(VS,GTAU,M,N,NM2) PRIVATE(K)
      do k=1,m
        VS(k,1,1) = GTAU(1,k)
        VS(k,1,NM2) = GTAU(N,k)
      enddo
!$OMP END PARALLEL DO
      IF (NM3 >= 2)THEN
!$OMP  PARALLEL DO PRIVATE(I,K)
!$OMP& SHARED(NM3,M,VS,GTAU)
        DO I=2,NM3
          DO K=1,M
            VS(K,1,I) = GTAU(I+1,K)
          ENDDO
        ENDDO
!$OMP END PARALLEL DO
      ENDIF

!$OMP PARALLEL DO PRIVATE(I,K)
!$OMP&    SHARED(NM2,M,VS,W)
      DO I = 2,NM2
        DO K=1,M
          VS(K,2,I) = (VS(K,1,I)-VS(K,1,I-1))/W(I,1)
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

      DTAU = TAU(2)-TAU(1)
      RATIO = DTAU/W(2,1)
      W(1,2) = (RATIO-1.)*(RATIO-1.)
      W(1,1) = RATIO*(RATIO-1.)
      C1 = RATIO*(2.*RATIO-3.)

!$OMP  PARALLEL DO PRIVATE(K)
!$OMP& SHARED(M,VS,GTAU,DTAU,C1)
      DO K=1,M
        VS(K,2,1) = (GTAU(2,K)-GTAU(1,K))/DTAU+VS(K,2,2)*C1
      ENDDO
!$OMP END PARALLEL DO

      IF (NM3 >= 2) THEN
        DO I=2,NM3
          C1 = 3.*W(I,1)
          C2 = 3.*W(I+1,1)
          G = -W(I+1,1)/W(I-1,2)
          W(I,2) = G*W(I-1,1)+2.*(W(I,1)+W(I+1,1))
!$OMP  PARALLEL DO PRIVATE(K)
!$OMP& SHARED(M,I,VS,G,C1,C2)
          DO K=1,M
            VS(K,2,I) = G*VS(K,2,I-1)+C1*VS(K,2,I+1)+C2*VS(K,2,I)
          ENDDO
!$OMP END PARALLEL DO
        ENDDO
      ENDIF

      DTAU = TAU(N-1)-XILIM

      IF(ANY(W(1:NM3,2) == 0.0))RETURN
      IER = 13

      RATIO = DTAU/W(NM2,1)
      G = -(RATIO-1.)*(RATIO-1.)/W(NM3,2)
      W(NM2,2) = RATIO*(RATIO-1.)

      C1 = RATIO*(2.*RATIO-3.)
!$OMP  PARALLEL DO PRIVATE(K)
!$OMP& SHARED(M,VS,GTAU,C1,NM1,NM2,NM3,DTAU)
      DO K=1,M
        VS(K,2,NM2) = VS(K,2,NM2)*C1+
     *                (GTAU(NM1,K)-VS(K,1,NM3))/DTAU
      ENDDO
!$OMP END PARALLEL DO
      W(NM2,2) = G*W(NM3,1)+W(NM2,2)

      IF(W(NM2,2) == 0.0)RETURN
      IER = 14

!$OMP  PARALLEL DO PRIVATE(K)
!$OMP& SHARED(M,VS,G,W,NM2,NM3)
      DO K=1,M
        VS(K,2,NM2) = (G*VS(K,2,NM3)+VS(K,2,NM2))/W(NM2,2)
      ENDDO
!$OMP END PARALLEL DO
      DO J = NM3,1,-1
!$OMP  PARALLEL DO PRIVATE(K)
!$OMP& SHARED(M,J,VS,W)
        DO K=1,M
          VS(K,2,J) = (VS(K,2,J)-W(J,1)*VS(K,2,J+1))/W(J,2)
        ENDDO
!$OMP END PARALLEL DO
      ENDDO

!$OMP  PARALLEL DO PRIVATE(J,K)
!$OMP& SHARED(M,VS,N,NM2)
      DO J=1,2
        DO K=1,M
          VS(K,J,N) = VS(K,J,NM2)
        ENDDO
      ENDDO
!$OMP END PARALLEL DO
      DO I=NM1,2,-1
!$OMP  PARALLEL DO PRIVATE(J,K)
!$OMP& SHARED(M,VS,I)
        DO J=1,2
          DO K=1,M
            VS(K,J,I) = VS(K,J,I-1)
          ENDDO
        ENDDO
!$OMP END PARALLEL DO
      ENDDO

      DO J=2,NM1,NM3  !j=2,j=NM1 only
        JM1 = J-1
        JP1 = J+1
C-----------------------------------------------------------------------
C       NEXT CORRECTION IS VALID FOR SMALL AREA N <= 4 ONLY
C-----------------------------------------------------------------------
        IF (JM1 == 2) JM1 = 1
        IF (JP1 == NM1) JP1 = N
        H = TAU(JP1)-TAU(JM1)
        U = TAU(J)-TAU(JM1)
!$OMP  PARALLEL DO PRIVATE(K,AA,BB,CC,DD)
!$OMP& SHARED(VS,M,J,JM1,JP1,H,U)
        DO K=1,M
          AA = VS(K,1,JM1)
          BB = VS(K,2,JM1)
          CC = (3.*(VS(K,1,JP1)-VS(K,1,JM1))/H-(VS(K,2,JP1)+
     *          2.*VS(K,2,JM1)))/H
          DD = (2.*(VS(K,1,JM1)-VS(K,1,JP1))/H+(VS(K,2,JP1)+
     *          VS(K,2,JM1)))/(H*H)
          VS(K,1,J) = AA+U*(BB+U*(CC+DD*U))
          VS(K,2,J) = BB+U*(2.*CC+3.*DD*U)
        ENDDO
!$OMP END PARALLEL DO
      ENDDO

      IER = 0

      END SUBROUTINE GET_SPLINES
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C
C   BICUBIC SPLINE TWO-DIMENSIONAL COEFFICIENT CALCULATOR
C
C    RESTRUCTURIZATION OF CODE FOR FORTRAN 90
C    BASED ON IMSL ROUTINE IBCCCU
C
C-----------------------------------------------------------------------
C
C   USAGE               - CALL IBCCCU (F,X,NX,Y,NY,C,IC,WK,IER)
C
C   ARGUMENTS    F      - NX BY NY MATRIX CONTAINING THE FUNCTION
C                           VALUES. (INPUT) F(I,J) IS THE FUNCTION VALUE
C                           AT THE POINT (X(I),Y(J)) FOR I=1,...,NX AND
C                           J=1,...,NY.
C                X      - VECTOR OF LENGTH NX. (INPUT) X MUST BE
C                           ORDERED SO THAT X(I) .LT. X(I+1) FOR
C                           I=1,...,NX-1.
C                NX     - NUMBER OF ELEMENTS IN X. (INPUT) NX MUST BE
C                           .GE. 4.
C                Y      - VECTOR OF LENGTH NY. (INPUT) Y MUST BE
C                           ORDERED SO THAT Y(J) .LT. Y(J+1) FOR
C                           J=1,...,NY-1.
C                NY     - NUMBER OF ELEMENTS IN Y. (INPUT) NY MUST BE
C                           .GE. 4.
C                         NOTE - THE COORDINATE PAIRS (X(I),Y(J)), FOR
C                           I=1,...,NX AND J=1,...,NY, GIVE THE POINTS
C                           WHERE THE FUNCTION VALUES F(I,J) ARE
C                           DEFINED.
C                C      - ARRAY OF SPLINE COEFFICIENTS. (OUTPUT)
C                           C IS OF DIMENSION 2 BY NX BY 2 BY NY.
C                           AT THE POINT (X(I),Y(J))
C                             C(1,I,1,J) = S
C                             C(2,I,1,J) = DS/DX
C                             C(1,I,2,J) = DS/DY
C                             C(2,I,2,J) = D(DS/DX)/DY
C                           WHERE S(X,Y) IS THE SPLINE APPROXIMATION.
C                           (NOTE - C IS TREATED INTERNALLY AS A
C                             2 BY NX BY 2*NY ARRAY BECAUSE CERTAIN
C                             ENVIRONMENTS DO NOT PERMIT QUADRUPLY-
C                             DIMENSIONED ARRAYS.  IN THESE
C                             ENVIRONMENTS THE CALLING PROGRAM MAY
C                             DIMENSION C IN THE SAME MANNER.)
C                IC     - ROW DIMENSION OF MATRIX F AND SECOND
C                           DIMENSION OF ARRAY C EXACTLY AS
C                           SPECIFIED IN THE DIMENSION STATEMENT.
C                           (INPUT). IC MUST BE .GE. NX.
C                WK     - WORK VECTOR OF LENGTH
C                           2*NX*NY+2*MAX(NX,NY)
C                IER    - ERROR PARAMETER. (OUTPUT)
C-----------------------------------------------------------------------
C
      SUBROUTINE IBCCCU (F,X,NX,Y,NY,C,IC,WK,IER)
C
C                             SPECIFICATIONS FOR ARGUMENTS
C
      IMPLICIT                NONE
      INTEGER, INTENT(IN )::    NX,NY,IC
      INTEGER, INTENT(OUT)::    IER
      REAL,    INTENT(IN )::    F(IC,NY),X(NX),Y(NY)
      REAL,    INTENT(OUT)::    C(2,IC,2*NY),WK(2*NX*NY+2*MAX(NX,NY))
C
C                             SPECIFICATIONS FOR LOCAL VARIABLES
C
      INTEGER::                IWK

      INTERFACE
          SUBROUTINE IBCDCU (TAU,GTAU,N,M,W,VS,IC1,IC2,IER)
            INTEGER, INTENT(IN)::   N,M,IC1,IC2
            INTEGER, INTENT(OUT)::  IER
            REAL, INTENT(IN)::      TAU(N),GTAU(IC1,M)
            REAL, INTENT(OUT)::     W(N,2),VS(IC2,2,N)
        END SUBROUTINE
      END INTERFACE
C
C                             FIRST EXECUTABLE STATEMENT
C
      IF (IC .LT. NX)THEN
        IER = 1
        RETURN
      ENDIF

      IF (NX .LT. 4) THEN
        IER = 2
        RETURN
      ENDIF

      IF (NY .LT. 4) THEN
        IER = 3
        RETURN
      ENDIF

      IWK = 2*NY*NX
      CALL IBCDCU(X,F,NX,NY,WK(IWK+1),WK,IC,NY,IER)

      IF (IER .GT. 0) RETURN
      CALL IBCDCU(Y,WK,NY,2*NX,WK(IWK+1),C,NY,2*IC,IER)

      END SUBROUTINE IBCCCU
C
C-----------------------------------------------------------------------
C
      SUBROUTINE IBCDCU (TAU,GTAU,N,M,W,VS,IC1,IC2,IER)
C
C                             SPECIFICATIONS FOR ARGUMENTS
C
      IMPLICIT                NONE
      INTEGER, INTENT(IN)::   N,M,IC1,IC2
      INTEGER, INTENT(OUT)::  IER
      REAL, INTENT(IN)::      TAU(N),GTAU(IC1,M)
      REAL, INTENT(OUT)::     W(N,2),VS(IC2,2,N)
C
C                             SPECIFICATIONS FOR LOCAL VARIABLES
C
      INTEGER                 I,JJ,JM1,JP1,J,K,LIM,LP1,NM1
      REAL                    AA,BB,C1,C2,CC,DD,DTAU,G,H,RATIO,U,XILIM
C
C                             FIRST EXECUTABLE STATEMENT
C
      LIM = N-3
      NM1 = N-1
      LP1 = LIM+1

      IER = 11

      W(2,1) = TAU(3)-TAU(1)
      IF (W(2,1) .LE. 0.0) RETURN

      XILIM = TAU(1)
      IF (LIM.GE.2)THEN
          XILIM = TAU(N-2)
CC!$PAR DOSERIAL
          DO I=2,LIM
              J = I+1
              W(J,1) = TAU(I+2)-TAU(J)
              IF (W(J,1).LE.0.0) RETURN
          ENDDO
      ENDIF
      W(LP1,1) = TAU(N)-XILIM
      IF (W(LP1,1).LE.0.0) RETURN

      IF (LIM.GE.2)THEN  !PARALLELIZED CASE
C-----PARALLELIZATION---------------------------------------------------
C
C     OpenMP STYLE
C
CC!$OMP PARALLEL DEFAULT(NONE)
CC!$OMP DO PRIVATE(K)
CC!$OMP&    SHARED(M, VS, GTAU)
C
C     SUN STYLE
C
CC!$PAR DOALL PRIVATE(K)
CC!$PAR&    SHARED(M, VS, GTAU)
        DO K=1,M
          VS(K,1,1) = GTAU(1,K)
        ENDDO
CC!$OMP END DO NOWAIT

CC!$OMP DO PRIVATE(K)
CC!$OMP& SHARED(M,N,VS,GTAU)
C
C     SUN STYLE
C
CC!$PAR DOALL PRIVATE(K)
CC!$PAR&    SHARED(M,N,VS,GTAU)
        DO K=1,M
          VS(K,1,LP1) = GTAU(N,K)
        ENDDO
CC!$OMP END DO NOWAIT

CC!$OMP DO PRIVATE(I,J,K)
CC!$OMP& SHARED(LIM,M,VS,GTAU)
C
C     SUN STYLE
C
CC!$PAR DOALL PRIVATE(I,J,K)
CC!$PAR&    SHARED(LIM,M,VS,GTAU)
CC!$PAR&    READONLY(LIM)
        DO I=2,LIM
          J = I+1
          DO K=1,M
            VS(K,1,I) = GTAU(J,K)
          ENDDO
        ENDDO
CC!$OMP END DO NOWAIT
CC!$OMP END PARALLEL
      ELSE
CC!$OMP PARALLEL DEFAULT(NONE)
CC!$OMP DO PRIVATE(K)
CC!$OMP&    SHARED(M, VS, GTAU)
C
C     SUN STYLE
C
CC!$PAR DOALL PRIVATE(K)
CC!$PAR&    SHARED(M, VS, GTAU)
C
        DO K=1,M
          VS(K,1,1) = GTAU(1,K)
        ENDDO
CC!$OMP END DO NOWAIT

CC!$OMP DO PRIVATE(K)
CC!$OMP& SHARED(M,N,VS,GTAU)
C
C     SUN STYLE
C
CC!$PAR DOALL PRIVATE(K)
CC!$PAR&    SHARED(M,N,VS,GTAU)
C
        DO K=1,M
          VS(K,1,LP1) = GTAU(N,K)
        ENDDO
CC!$OMP END DO NOWAIT
CC!$OMP END PARALLEL
      ENDIF
C
C     CAN PARALLELIZE EXTERNAL LOOP
C
CC!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(I,K)
CC!$OMP&    SHARED(LP1,M,VS,W)
C
C     SUN STYLE
C
CC!$PAR DOALL PRIVATE(I,K)
CC!$PAR&    SHARED(LP1,M,VS,W)
CC!$PAR&    READONLY(LP1)
C
      DO I = 2,LP1
        DO K=1,M
          VS(K,2,I) = (VS(K,1,I)-VS(K,1,I-1))/W(I,1)
        ENDDO
      ENDDO

      DTAU = TAU(2)-TAU(1)
      RATIO = DTAU/W(2,1)
      W(1,2) = (RATIO-1.)*(RATIO-1.)
      W(1,1) = RATIO*(RATIO-1.)
      C1 = RATIO*(2.*RATIO-3.)

CC!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(K)
CC!$OMP&    SHARED(M,VS,GTAU,DTAU,C1)
C
C     SUN STYLE
C
CC!$PAR DOALL PRIVATE(K)
CC!$PAR&    SHARED(M,VS,GTAU,DTAU,C1)
CC!$PAR&    READONLY(DTAU,C1)
C
      DO K=1,M
        VS(K,2,1) = (GTAU(2,K)-GTAU(1,K))/DTAU+VS(K,2,2)*C1
      ENDDO

      IF (LIM.GE.2) THEN
C
C     CAN NOT PARALLELIZE EXTERNAL LOOP
C     - RECURENT SEQUANTIAL UPDATE OF VS IS GOING
C     INTERNAL ONLY - SMALL GAIN (??)
C
C     SUN STYLE
C
CC!$PAR DOSERIAL
        DO I=2,LIM
          J = I+1
          JJ = I-1
          G = -W(J,1)/W(JJ,2)
          C1 = 3.*W(I,1)
          C2 = 3.*W(J,1)
          W(I,2) = G*W(JJ,1)+2.*(W(I,1)+W(J,1))
C
C     SUN STYLE
C
CC!$PAR DOALL PRIVATE(K)
CC!$PAR&    SHARED(M,I,J,JJ,VS,G,C1,C2)
CC!$PAR&    READONLY(I,J,JJ,G,C1,C2)
C
          DO K=1,M
            VS(K,2,I) = G*VS(K,2,JJ)+C1*VS(K,2,J)+C2*VS(K,2,I)
          ENDDO
        ENDDO
      ENDIF

      DTAU = TAU(N-1)-XILIM
      RATIO = DTAU/W(LP1,1)
      G = -(RATIO-1.)*(RATIO-1.)/W(LIM,2)
      W(LP1,2) = RATIO*(RATIO-1.)
      C1 = RATIO*(2.*RATIO-3.)
C
C     SUN STYLE
C
CC!$PAR DOALL PRIVATE(K)
CC!$PAR&    SHARED(VS,M,N,LP1,LIM,DTAU,C1)
CC!$PAR&    READONLY(LP1,LIM,DTAU,C1)
C
      DO K=1,M
        VS(K,2,LP1) = (GTAU(N-1,K)-VS(K,1,LIM))/DTAU+VS(K,2,LP1)*C1
      ENDDO

      W(LP1,2) = G*W(LIM,1)+W(LP1,2)
C
C     SUN STYLE
C
CC!$PAR DOALL PRIVATE(K)
CC!$PAR&    SHARED(VS,W,M,LP1,LIM,G)
CC!$PAR&    READONLY(LP1,LIM,G)
C
      DO K=1,M
          VS(K,2,LP1) = (G*VS(K,2,LIM)+VS(K,2,LP1))/W(LP1,2)
      ENDDO
C
C     CAN NOT PARALLELIZE EXTERNAL LOOP
C     - RECURENT SEQUANTIAL UPDATE OF VS IS GOING
C     INTERNAL ONLY - SMALL GAIN (??)
C
C     SUN STYLE
C
CC!$PAR DOSERIAL
      DO J = LIM,1,-1
        JJ = J+1
C
C     SUN STYLE
C
CC!$PAR DOALL PRIVATE(K)
CC!$PAR&    SHARED(VS,W,M,J,JJ)
CC!$PAR&    READONLY(J,JJ)
C
        DO K=1,M
          VS(K,2,J) = (VS(K,2,J)-W(J,1)*VS(K,2,JJ))/W(J,2)
        ENDDO
      ENDDO
C
C     CAN NOT PARALLELIZE EXTERNAL LOOP
C     - RECURENT SEQUANTIAL UPDATE OF VS IS GOING
C     INTERNAL ONLY - SMALL GAIN (??)
C
C     SUN STYLE
C
CC!$PAR DOSERIAL
      DO JJ=1,N
        J = N+1-JJ
        JM1 = J-1
        IF (J.EQ.N) JM1 = J-2
        IF (J.EQ.1) JM1 = J
C
C     SUN STYLE
C
CC!$PAR DOALL PRIVATE(K)
CC!$PAR&    SHARED(VS,M,J,JM1)
CC!$PAR&    READONLY(J,JM1)
C
        DO K=1,M
          VS(K,1,J) = VS(K,1,JM1)
          VS(K,2,J) = VS(K,2,JM1)
        ENDDO
      ENDDO
C
C     NEXT GIVES ONLY TWO J VALUES:2 AND N-1=NM1
C
CC!$PAR DOSERIAL
      DO J=2,NM1,LIM
        JM1 = J-1
        JP1 = J+1
C       NEXT CORRECTION IS VALID FOR SMALL AREA N <= 4 ONLY
        IF (JM1.EQ.2) JM1 = 1
        IF (JP1.EQ.NM1) JP1 = N

        H = TAU(JP1)-TAU(JM1)
        U = TAU(J)-TAU(JM1)
C
C     SUN STYLE
C
CC!$PAR DOALL PRIVATE(K,AA,BB,CC,DD)
CC!$PAR&    SHARED(VS,M,J,JM1,JP1,H,U)
CC!$PAR&    READONLY(J,JM1,JP1,H,U)
C
        DO K=1,M
          AA = VS(K,1,JM1)
          BB = VS(K,2,JM1)
          CC = (3.*(VS(K,1,JP1)-VS(K,1,JM1))/H-(VS(K,2,JP1)+
     *          2.*VS(K,2,JM1)))/H
          DD = (2.*(VS(K,1,JM1)-VS(K,1,JP1))/H+(VS(K,2,JP1)+
     *          VS(K,2,JM1)))/(H*H)
          VS(K,1,J) = AA+U*(BB+U*(CC+DD*U))
          VS(K,2,J) = BB+U*(2.*CC+3.*DD*U)
        ENDDO
      ENDDO

      IER = 0

      END SUBROUTINE IBCDCU
C
C******************************************************************
C
!#define TEST_SPLINES

#ifdef TEST_SPLINES
      PROGRAM MAIN
      implicit none
      INTEGER, PARAMETER::    NX = 7, NY = 7, IC = NX
      INTEGER  IER,i,j
      REAL     F(IC,NY),X(NX),Y(NY)
      REAL     C(2,IC,2*NY)

      c=0.
      do j=1,ny
      y(j)=j
      do i=1,nx
      x(i)=i
      f(i,j)=(i**2 + J**2)*0.1
      enddo
      enddo
      write(*,*)"Old version"
      call GET_BSPLINES (F,X,NX,Y,NY,C,IC,IER)
      write(*,*)"iret=",ier
      do j=1,ny
      write(*,"(i3,10(f8.3,f8.3))")j,(C(1,I,j),C(2,I,j),i=1,nx)
      write(*,"(i3,10(f8.3,f8.3))")ny+j,(C(1,I,ny+j),C(2,I,ny+j),i=1,nx)
      enddo
      c=0.
      write(*,*)"New version"
      call GET_BSPLINES_S (F,X,NX,Y,NY,C,IC,IER)
      write(*,*)"iret=",ier
      do j=1,ny
      write(*,"(i3,10(f8.3,f8.3))")j,(C(1,I,j),C(2,I,j),i=1,nx)
      write(*,"(i3,10(f8.3,f8.3))")ny+j,(C(1,I,ny+j),C(2,I,ny+j),i=1,nx)
      enddo

      END PROGRAM
#endif
