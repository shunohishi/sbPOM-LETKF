
!#define TEST_MASK
C******************************************************************
c      SUBROUTINE EXTRAPOLATE_BY_MASK   ( F, MASK, NX, NY, NX_D, IRET )
c      SUBROUTINE EXTRAPOLATE_BY_MASK_V5( F, MASK, NX, NY, NX_F, IRET )
c
c      SUBROUTINE MASK_REPLACE_MISSING
c     *  (IRET, DAT, NX, NY, NXDIM, UNDEF_MIN, UNDEF_MAX, MASK)
c      SUBROUTINE MASK_REPLACE_MISSING_M
c     *  (IRET, R_MEAN, DAT, NX, NY, NXDIM, UNDEF_MIN, UNDEF_MAX, MASK)
c
c     Sergey Varlamov, 2000-2008
c
C******************************************************************

      SUBROUTINE EXTRAPOLATE_BY_MASK( F, MASK, NX, NY, NX_D, IRET )

C******************************************************************
C
C     ACCORDING TO MASK VALUE CHANGE THE LAND VALUES TO INVERSE
C     OF DISTANCE WEIGHTED VALUES OF NEAREST SEA POINTS.
C     IF NO SEA POINT AT R-MAX GRID STEPS  DISTANCE, THE MEAN VALUE
C     IS USED. IF ALL POINTS ARE LAND POINTS - RETURN ERROR.
C
C     S. VARLAMOV
C
C     LAST REVISION:  MAY 10, 2001
C
C*******************************************************************

      IMPLICIT                  NONE
      INTEGER, INTENT(IN )::    NX, NY, NX_D
      INTEGER, INTENT(OUT)::    IRET
      REAL,  INTENT(INOUT)::    F(NX_D,NY)
      REAL,    INTENT(IN )::    MASK(NX_D,NY)

      REAL,PARAMETER::          R_MAX = 6., RINV = 0.1/R_MAX
      INTEGER                   I, J, I1, J1, I2, J2 ,JJ, II, R
      REAL                      RS, RR, R_MEAN, W
      DOUBLE PRECISION          FSUM
      IRET = 0
C
C     ESTIMATE THE FIELD MEAN VALUE FOR 'FAR' EXTRAPOLATION
C
      RR = 0.
!      RS = 0.
      FSUM = 0D0
!$OMP  PARALLEL DO PRIVATE(i,j) REDUCTION(+:RR,FSUM)
!$OMP& SHARED(NY,NX,MASK,F)
      DO J=1,NY
          DO I=1,NX
              IF(MASK(I,J) > 0.5)THEN !DO NOT LAND POINT
                RR = RR + 1.
                FSUM = FSUM + F(I,J)
!                W=(F(I,J)-RS)/RR
!                RS=RS+W
              ENDIF
          ENDDO
      ENDDO
!$OMP END PARALLEL DO
      IF(RR < 1.)THEN
        IRET = 1
C       WRITE(*,*)' ENTIRE FIELD UNDEFINED!'
        RETURN
      ENDIF
      R_MEAN = FSUM/RR

      DO J=1,NY
          DO I=1,NX
              IF(MASK(I,J) < 0.5)THEN !LAND POINT
C
C                INTERPOLATE AS WEIGHTED VALUE FROM NEAREST SEA POINTS
C
                  R=0
                  DO WHILE (.TRUE.)
                      R=R+1
                      IF(R > R_MAX)THEN
C                         'FAR' INTERPOLATION, NOT IMPORTANT
                          F(I,J)=R_MEAN
                          EXIT
                      ENDIF

                      I1=I-R
                      IF(I1 < 1)I1=1
                      I2=I+R
                      IF(I2 > NX)I2=NX
                      J1=J-R
                      IF(J1 < 1)J1=1
                      J2=J+R
                      IF(J2 > NY)J2=NY
                      RS=0. !SUM OF WEIGHTS
                      W=0.  !SUM OF VALUES
                      DO JJ = J1,J2
                           DO II = I1,I2
                               IF(MASK(II,JJ) < 0.5)CYCLE !LAND POINT
                               RR=(JJ-J)*(JJ-J)+(II-I)*(II-I)
                               RR=SQRT(RR)
                               !CHECK THE DISTANCE AS R+0.5
                               IF(RR > R+0.5)CYCLE
                               RR=1./RR
                               RS = RS + RR
                               W = W + RR*F(II,JJ)
                           ENDDO
                      ENDDO
                      IF(RS >= RINV)THEN
                        F(I,J) = W/RS
                        EXIT
                      ENDIF
                  ENDDO
              ENDIF
          ENDDO
      ENDDO

      END SUBROUTINE EXTRAPOLATE_BY_MASK

      SUBROUTINE EXTRAPOLATE_BY_MASK_V5( F, MASK, NX, NY, NX_F, IRET )

C******************************************************************
C
C     ACCORDING TO MASK VALUE CHANGE THE LAND VALUES TO INVERSE
C     OF DISTANCE WEIGHTED VALUES OF NEAREST SEA POINTS.
C     IF NO SEA POINT AT R-MAX GRID STEPS  DISTANCE, THE MEAN VALUE
C     IS USED. IF ALL POINTS ARE LAND POINTS - RETURN ERROR.
C
C     S. VARLAMOV
C
C     LAST REVISION:  MAY 10, 2001
C
C*******************************************************************

      IMPLICIT                  NONE
      INTEGER, INTENT(IN   )::  NX, NY, NX_F
      INTEGER, INTENT(OUT  )::  IRET
      REAL,    INTENT(INOUT)::  F(NX_F,NY)
      REAL,    INTENT(IN   )::  MASK(NX,NY)
      INTEGER,PARAMETER::       IJ_MAX = 4
      REAL,PARAMETER::          W0 = 1E-7, R_MAX = IJ_MAX+0.5
      INTEGER                   I, J, I1, J1, I2, J2 ,JJ, II
      REAL                      R_MEAN, Q
      DOUBLE PRECISION          RS, W, QR


      IRET = 0
C
C     ESTIMATE THE FIELD MEAN VALUE FOR 'FAR' EXTRAPOLATION
C
      W = 0D0
      RS = 0D0
!$OMP  PARALLEL DO PRIVATE(i,j) REDUCTION(+:W,RS)
!$OMP& SHARED(NY,NX,MASK,F)
      DO J=1,NY
        DO I=1,NX
          IF(MASK(I,J) > 0.5)THEN !SEA POINT
            W = W + 1D0
            RS = RS + F(I,J)
          ENDIF
        ENDDO
      ENDDO
!$OMP END PARALLEL DO
      IF(W < 1D0)THEN
        IRET = 1
C       WRITE(*,*)' ENTIRE FIELD UNDEFINED!'
        RETURN
      ENDIF

      R_MEAN = RS/W

!$OMP  PARALLEL DO PRIVATE(i,j,j1,j2,i1,i2,ii,jj,W,RS,Q,QR)
!$OMP& SHARED(NY,NX,MASK,F,R_MEAN)
      DO J=1,NY
C        J1 = MAX(J-IJ_MAX,1 )
C        J2 = MIN(J+IJ_MAX,NY)
        J1 = J-IJ_MAX
        J2 = J+IJ_MAX
        DO I=1,NX
          IF(MASK(I,J) < 0.5)THEN !LAND POINT
C
C           INTERPOLATE AS WEIGHTED VALUE FROM NEAREST SEA POINTS
C
C            I1 = MAX(I-IJ_MAX,1 )
C            I2 = MIN(I+IJ_MAX,NX)
            I1 = I - IJ_MAX
            I2 = I + IJ_MAX
            W  = W0         !SUM OF WEIGHTS
            RS = R_MEAN*W0  !SUM OF VALUES
            DO JJ = J1,J2
              DO II = I1,I2
                IF(II>0 .AND. JJ>0 .AND. II<=NX .AND. JJ<=NY)THEN
                  IF(MASK(II,JJ) > 0.5)THEN !SEA POINT
                    Q = SQRT(FLOAT((JJ-J)**2+(II-I)**2))
                    !CHECK THE DISTANCE
                    IF(Q < R_MAX)THEN
                      QR=1./Q
                      RS = RS + QR*F(II,JJ)
                      W = W + QR
                    ENDIF
                  ENDIF
                ENDIF
              ENDDO
            ENDDO
            F(I,J) = RS/W
          ENDIF
        ENDDO
      ENDDO
!$OMP END PARALLEL DO
      END SUBROUTINE EXTRAPOLATE_BY_MASK_V5


      SUBROUTINE MASK_REPLACE_MISSING
     *  (IRET, DAT, NX, NY, NXDIM, UNDEF_MIN, UNDEF_MAX, MASK)


C******************************************************************
C
C     ACCORDING TO MASK VALUE AND LIMITS FOR FUNCTION VALUES
C     CHANGE BAD AND/OR MISSING VALUES TO THE INVERSE
C     OF DISTANCE-WEIGHTED VALUES OF NEAREST SEA POINTS.
C     IF NO SEA POINT AT R-MAX GRID STEPS  DISTANCE, THE MEAN VALUE
C     IS USED.
C     IF ALL POINTS ARE LAND POINTS - RETURN ERROR=2.
C
C     S. VARLAMOV
C
C     LAST REVISION:  JULY 27, 2004
C
C*******************************************************************

      IMPLICIT                  NONE
      INTEGER, INTENT(IN   )::  NX, NY, NXDIM
      INTEGER, INTENT(OUT  )::  IRET
      REAL,    INTENT(INOUT)::  DAT(NXDIM,NY)
      REAL,    INTENT(IN   ),OPTIONAL::
     *                          MASK(NX,NY), UNDEF_MIN, UNDEF_MAX
      INTEGER,PARAMETER::       IJ_MAX = 4
      REAL,PARAMETER::          W0 = 1E-7, R_MAX = IJ_MAX+0.5
      INTEGER                   I, J, I1, J1, I2, J2 ,JJ, II
      REAL                      R_MEAN, Q
      DOUBLE PRECISION          RS, W, QR
      LOGICAL                   TMP_MASK(NX,NY)

      IRET = 0
      TMP_MASK = .TRUE.   !ALL SEA
C
C     ESTIMATE THE FIELD MEAN VALUE FOR 'FAR' EXTRAPOLATION
C
      IF(PRESENT(UNDEF_MIN))THEN
!$OMP PARALLEL WORKSHARE SHARED(DAT,NX,NY,TMP_MASK,UNDEF_MIN)
        WHERE(DAT(1:NX,1:NY) <= UNDEF_MIN)
          TMP_MASK(:,:) = .FALSE.   !LAND OR MISSING
        ENDWHERE
!$OMP END PARALLEL WORKSHARE
      ENDIF
      IF(PRESENT(UNDEF_MAX))THEN
!$OMP PARALLEL WORKSHARE SHARED(DAT,NX,NY,TMP_MASK,UNDEF_MAX)
        WHERE(DAT(1:NX,1:NY) >= UNDEF_MAX)
          TMP_MASK(:,:) = .FALSE.   !LAND OR MISSING
        ENDWHERE
!$OMP END PARALLEL WORKSHARE
      ENDIF
      IF(PRESENT(MASK))THEN
!$OMP PARALLEL WORKSHARE SHARED(MASK,TMP_MASK)
        WHERE(MASK < 0.5)      !LAND POINT
          TMP_MASK = .FALSE.   !LAND OR MISSING
        ENDWHERE
!$OMP END PARALLEL WORKSHARE
      ENDIF

      W = 0D0
      RS = 0D0
!$OMP  PARALLEL DO PRIVATE(i,j) REDUCTION(+:W,RS)
!$OMP& SHARED(NY,NX,TMP_MASK,DAT)
      DO J=1,NY
        DO I=1,NX
          IF(TMP_MASK(I,J))THEN !SEA POINT
            W = W + 1D0
            RS = RS + DAT(I,J)
          ENDIF
        ENDDO
      ENDDO
!$OMP END PARALLEL DO
      IF(W < 1D0)THEN
        IRET = 2
C       WRITE(*,*)' ENTIRE FIELD UNDEFINED!'
        RETURN
      ENDIF

      R_MEAN = RS/W

!$OMP  PARALLEL DO PRIVATE(i,j,j1,j2,i1,i2,ii,jj,W,RS,Q,QR)
!$OMP& SHARED(NY,NX,TMP_MASK,DAT)
      DO J=1,NY
        J1 = J-IJ_MAX
        J2 = J+IJ_MAX
        DO I=1,NX
          IF(.NOT. TMP_MASK(I,J))THEN !LAND OR MISSING DATA POINT
C
C           INTERPOLATE AS WEIGHTED VALUE FROM NEAREST SEA POINTS
C
            I1 = I - IJ_MAX
            I2 = I + IJ_MAX
            W  = W0         !SUM OF WEIGHTS
            RS = R_MEAN*W0  !SUM OF VALUES
            DO JJ = J1,J2
              DO II = I1,I2
                IF(II>0 .AND. JJ>0 .AND. II<=NX .AND. JJ<=NY)THEN
                  IF(TMP_MASK(II,JJ))THEN !SEA POINT
                    Q = SQRT(FLOAT((JJ-J)**2+(II-I)**2))
                    !CHECK THE DISTANCE
                    IF(Q < R_MAX)THEN
                      QR=1./Q
                      RS = RS + QR*DAT(II,JJ)
                      W = W + QR
                    ENDIF
                  ENDIF
                ENDIF
              ENDDO
            ENDDO
            DAT(I,J) = RS/W
          ENDIF   !LAND OR MISSING DATA POINT
        ENDDO
      ENDDO
!$OMP END PARALLEL DO
      END SUBROUTINE MASK_REPLACE_MISSING


      SUBROUTINE MASK_REPLACE_MISSING_M
     *  (IRET, R_MEAN, DAT, NX, NY, NXDIM, UNDEF_MIN, UNDEF_MAX, MASK)

C******************************************************************
C
C     ACCORDING TO MASK VALUE AND LIMITS FOR FUNCTION VALUES
C     CHANGE BAD AND/OR MISSING VALUES TO THE INVERSE
C     OF DISTANCE-WEIGHTED VALUES OF NEAREST SEA POINTS.
C     IF NO SEA POINT AT R-MAX GRID STEPS  DISTANCE, THE MEAN VALUE
C     IS USED. IF ALL POINTS ARE LAND POINTS - RETURN ERROR.
C
C     S. VARLAMOV
C
C     LAST REVISION:  JULY 27, 2004
C
C*******************************************************************
#ifdef TEST_MASK
#ifdef _OPENMP
      USE OMP_LIB
#endif
#endif
      IMPLICIT                  NONE
      INTEGER, INTENT(IN   )::  NX, NY, NXDIM
      INTEGER, INTENT(OUT  )::  IRET
!      REAL,    INTENT(INOUT)::  DAT(NXDIM,NY)
      REAL,    INTENT(INOUT)::  DAT(NXDIM,*)
      REAL,    INTENT(OUT)::    R_MEAN
      REAL,    INTENT(IN   ),OPTIONAL::
     *                          MASK(NX,*), UNDEF_MIN, UNDEF_MAX
      INTEGER,PARAMETER::       IJ_MAX = 4
      REAL,PARAMETER::          W0 = 1E-7, R_MAX = IJ_MAX+0.5
      INTEGER                   I, J, I1, J1, I2, J2 ,JJ, II
      REAL                      Q
      DOUBLE PRECISION          RS, W, RQ, WQ, QR
      LOGICAL                   TMP_MASK(NX,NY)

      IRET = 0
      R_MEAN = HUGE(Q)
      W = 0D0
      RS = 0D0

!$OMP  PARALLEL DEFAULT(NONE) IF (NX*NY >= 400)
!$OMP& SHARED(UNDEF_MIN,UNDEF_MAX,TMP_MASK,MASK)
!$OMP& SHARED(DAT,NX,NY,W,RS,R_MEAN)
!$OMP& PRIVATE(i,j,j1,j2,i1,i2,ii,jj,WQ,RQ,Q,QR)
#ifdef TEST_MASK
#ifdef _OPENMP
      I = OMP_GET_NUM_THREADS()
      write(*,*)"Run on ",I,"threads"
#endif
#endif
!$OMP WORKSHARE
      TMP_MASK = .TRUE.   !ALL SEA
!$OMP END WORKSHARE
      IF(PRESENT(UNDEF_MIN))THEN
!$OMP WORKSHARE
        WHERE(DAT(1:NX,1:NY) <= UNDEF_MIN)
          TMP_MASK(1:NX,1:NY) = .FALSE.   !LAND OR MISSING
        ENDWHERE
!$OMP END WORKSHARE
      ENDIF
      IF(PRESENT(UNDEF_MAX))THEN
!$OMP WORKSHARE
        WHERE(DAT(1:NX,1:NY) >= UNDEF_MAX)
          TMP_MASK(1:NX,1:NY) = .FALSE.   !LAND OR MISSING
        ENDWHERE
!$OMP END WORKSHARE
      ENDIF
      IF(PRESENT(MASK))THEN
!$OMP WORKSHARE
        WHERE(MASK(1:NX,1:NY) < 0.5)      !LAND POINT
          TMP_MASK(1:NX,1:NY) = .FALSE.   !LAND OR MISSING
        ENDWHERE
!$OMP END WORKSHARE
      ENDIF
C
C     ESTIMATE THE FIELD MEAN VALUE FOR 'FAR' EXTRAPOLATION
C
!$OMP DO
!$OMP& REDUCTION(+:W,RS)
      DO J=1,NY
        DO I=1,NX
          IF(TMP_MASK(I,J))THEN !SEA POINT
            W = W + 1D0
            RS = RS + DAT(I,J)
          ENDIF
        ENDDO
      ENDDO
!$OMP END DO

      IF(W < 1D0)THEN
        R_MEAN = 0.
!        IRET = 2
!        RETURN
      ELSE
        R_MEAN = RS/W
      ENDIF
!$OMP DO
      DO J=1,NY
        J1 = J-IJ_MAX
        J2 = J+IJ_MAX
        DO I=1,NX
          IF(.NOT. TMP_MASK(I,J))THEN !LAND OR MISSING DATA POINT
C
C           INTERPOLATE AS WEIGHTED VALUE FROM THE NEAREST SEA POINTS
C
            I1 = I - IJ_MAX
            I2 = I + IJ_MAX
            WQ  = W0         !SUM OF WEIGHTS
            RQ = R_MEAN*W0   !SUM OF VALUES
            DO JJ = J1,J2
              DO II = I1,I2
                IF(II>0 .AND. JJ>0 .AND. II<=NX .AND. JJ<=NY)THEN !SEA POINT
                  IF( TMP_MASK(II,JJ) )THEN !SEA POINT
                    Q = SQRT(FLOAT((JJ-J)**2+(II-I)**2))
                    !CHECK THE DISTANCE
                    IF(Q < R_MAX)THEN
                      QR=1./Q
                      RQ = RQ + QR*DAT(II,JJ)
                      WQ = WQ + QR
                    ENDIF
                  ENDIF
                ENDIF
              ENDDO
            ENDDO
            DAT(I,J) = RQ/WQ
          ENDIF   !LAND OR MISSING DATA POINT
        ENDDO
      ENDDO
!$OMP END DO
!$OMP END PARALLEL

      END SUBROUTINE MASK_REPLACE_MISSING_M

C
C******************************************************************
C
#ifdef TEST_MASK
      PROGRAM MAIN
      implicit none
      INTEGER, PARAMETER::    NX = 21, NY = 21, NXDIM = NX
      INTEGER  IRET,i,j
      REAL F(NX,NY)
      REAL R_MEAN
      REAL MASK(NX,NY), UNDEF_MIN, UNDEF_MAX
      INTERFACE
              SUBROUTINE MASK_REPLACE_MISSING_M
     *        (IRET, R_MEAN, DAT, NX, NY, NXDIM,
     *        UNDEF_MIN, UNDEF_MAX, MASK)
                IMPLICIT                  NONE
                INTEGER, INTENT(IN   )::  NX, NY, NXDIM
                INTEGER, INTENT(OUT  )::  IRET
                REAL,    INTENT(INOUT)::  DAT(NXDIM,*)
                REAL,    INTENT(OUT)::    R_MEAN
                REAL,    INTENT(IN   ),OPTIONAL::
     *                          MASK(NX,NY), UNDEF_MIN, UNDEF_MAX
              END SUBROUTINE
      END INTERFACE

      do j=1,ny
      do i=1,nx
      f(i,j)=(i**2 + J**2)*0.01
      enddo
      enddo
      write(*,*)"Old version"
      CALL MASK_REPLACE_MISSING_M
     *  (IRET, R_MEAN, F, NX, NY, NXDIM, UNDEF_MIN=2.)
!, UNDEF_MAX, MASK)

      write(*,*)"iret=",IRET,R_MEAN

      END PROGRAM
#endif
