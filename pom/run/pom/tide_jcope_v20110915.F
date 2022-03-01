
      module jpom_tide
c
c       Now in this module is placed all information
c       related to tidal forcing for the JPOM.
c       Sergey Varlamov, 2009.01.08
c
c       Version January 29, 2009
c
c       modified by Yasumasa Miyazawa, 2012.02
c
        USE TIDE_SERV
        IMPLICIT NONE
c
c       solution specific to JPOM structure, define model size:
c
        include 'ssize.h'
!        private:: im, jm, km, rlons0, rlats0, rlon_step, rlat_step
        private:: im, jm, km


        INTEGER,PRIVATE,SAVE:: NF_LOG=6
        REAL,SAVE::
     *         ZN_TIDE(IM,JM),          !tidal sea surface anomaly from harmonics
     *         UVHC_TIDE(IM,JM,2),      !tidal volume fluxes from harmonics
     *         ZN_EQTR(IM,JM),          !astronomic equivalent tidal sea level anomaly
     *         ZN_MEAN(IM,JM),          !model time-filtered sea level
     *         UV_MEAN(IM,JM,2),           !model time-filtered 2D velocity field
     *         TTS_MEAN(IM,JM,KM,2),    !model time-filtered 3D temperature and salinity fields
     *         TUV_MEAN(IM,JM,KM,2)     !model time-filtered 3D velocity fields
        REAL,SAVE::
     *         ZN_TIDE_i(IM,JM,2),          !tidal sea surface anomaly from harmonics
     *         UVHC_TIDE_i(IM,JM,2,2),      !tidal volume fluxes from harmonics
     *         ZN_EQTR_i(IM,JM,2)           !astronomic equivalent tidal sea level anomaly

        REAL,ALLOCATABLE,PRIVATE,SAVE::
     *         TH_AS(:,:,:,:), TH_AC(:,:,:,:), !(IM,JM,NUM_TH,NUM_TPARAM)
     *         FS_LON(:,:),FC_LON(:,:),        !(0:2,IM)
     *         FLAT(:,:)                       !(NHARM_ALL,JM)

        TYPE(TIDE_D),SAVE::  SEA_TIDE = TIDE_D(
     *    "/home/G10201/vsm/TIDES", !PATH         (in)
     *    TIDE_NAO99_JAPAN_SSH_UVH, !TIDE_MODEL   (in)
     *    0,                        !NUM_TH       (inout)
     *    1,                        !ENABLE_ADM   (in)
     *    0,                        !NUM_TH_INPUT (out)
     *    0,                        !NUM_REGIONS  (inout)
     *    null(),                   !GRID(:)      (inout)
     *    null()                    !HARMONIC_NAME(:) (inout)
     *  )
C
C       SOME EXTRA PARAMETERS THAT USUALLY ARE DEFINED IN THE CALLING CODE
C
        REAL,SAVE:: R_DAYS_AVERAGE = 50./24.,  !20090702 changed from 48 to 50 h, vsm
     *              R_DAYS_DAMP = 2E0,
     *              CM_TO_USER_UNITS = 0.01,
     *              R_TD_AVERAGE = 0.
        DOUBLE PRECISION,PRIVATE,SAVE::
     *              UTC_TIME_TIDES_FROM=0D0, UTC_AVERAGE_FROM=0D0,
     *              TSEC_AVERAGE=0D0, TSEC_DAMP=1D0
        INTEGER,SAVE::
     *    TIDE_DEBUG_PRINT=1
        LOGICAL,SAVE::
     *    L_TIDE_BD=.false., !DOES IT NEEDS ESTIMATE TIDE IN MODEL (BOUNDARY) POINTS
     *    L_TIDE_FC=.false., !DOES IT NEEDS ESTIMATE EQUIVALENT ASTRONOMIC TIDE
     *    L_TIDE_UVH=.false. !DOES IT NEEDS ESTIMATE TIDAL VOLUME FLUXES (if possible)
        
! miyazawa 12.02.06
        REAL,PARAMETER:: START_TIDE_DELAY_DAYS=0.0
        DOUBLE PRECISION,SAVE:: UTC_TM_JDS, UTC_TM_JDS_START_TIDE,
     *   timejdso,timeincrs_mn,timeincrs_hf,time0jdso_hf,time0jdso_mn,
     *   d_tmp1,d_tmp2,timejdrtr,tide_io_time
        DOUBLE PRECISION,SAVE:: dd_tmp1=-1,dd_tmp2=-1

        PRIVATE !all except PUBLIC
        PUBLIC  ZN_TIDE,ZN_EQTR,ZN_MEAN,L_TIDE_BD,L_TIDE_FC,L_TIDE_UVH,
     *          JPOM_INIT_TIDE,GET_TIDE,FREE_TIDES,UV_MEAN,
     *          WRITE_TIDE_MEAN,READ_TIDE_MEAN,SEA_TIDE,
     *          UVHC_TIDE,GET_TIDE2,null_init_tide,
     *          init_td_mean,update_td_mean,R_DAYS_DAMP,
     *          R_DAYS_AVERAGE,CM_TO_USER_UNITS,TTS_MEAN,TUV_MEAN,
     *          TIDE_NAO99_GLOBAL,TIDE_NAO99_JAPAN_SSH,
     *          TIDE_NAO99_JAPAN_SSH_UVH,TIDE_LEGOS_2002,TIDE_FES2004,
     *          TIDE_OTIS_GLOBAL,TIDE_OTIS_PO,
     *          TIDE_OTIS_GLOBAL_ATLAS,TIDE_OTIS_PO_ATLAS,
     *          START_TIDE_DELAY_DAYS,UTC_TM_JDS,UTC_TM_JDS_START_TIDE,
     *          timejdso,timeincrs_mn,timeincrs_hf,
     *          time0jdso_hf,time0jdso_mn,
     *          d_tmp1,d_tmp2,timejdrtr,tide_io_time,
     *          dd_tmp1,dd_tmp2

      contains
C
C**********************************************************************C
C
      subroutine null_init_tide(DT)
        implicit none
        REAL,INTENT(IN),OPTIONAL:: DT
        ZN_TIDE=0
        UVHC_TIDE=0
        ZN_EQTR=0
        ZN_MEAN=0
        UV_MEAN=0
        TTS_MEAN=0
        TUV_MEAN=0
        ZN_TIDE_i=0
        UVHC_TIDE_i=0
        ZN_EQTR_i=0
        IF(PRESENT(DT))THEN
          TSEC_AVERAGE = 86400D0*R_DAYS_AVERAGE
          R_TD_AVERAGE  = DT/TSEC_AVERAGE
        endif
      end subroutine
C
      subroutine init_td_mean(UTC_TM_JDS,DT,NFL,ZN,UA,VA,T,S,U,V)
        use time_serv
        implicit none
        REAL,INTENT(IN),OPTIONAL:: ZN(IM,JM),UA(IM,JM),VA(IM,JM)
        REAL,INTENT(IN),OPTIONAL:: T(IM,JM,KM),S(IM,JM,KM)
        REAL,INTENT(IN),OPTIONAL:: U(IM,JM,KM),V(IM,JM,KM)
        REAL,INTENT(IN):: DT
        INTEGER,INTENT(IN):: NFL
        DOUBLE PRECISION, INTENT(IN ):: UTC_TM_JDS
        INTEGER IRET

        CALL CNV_TIME(UTC_TM_JDS,UTC_AVERAGE_FROM,NFL,IRET)
        IF(IRET /= 0)THEN
          WRITE(NFL,*)' *** TIDE AND POM CALENDARS DIVERGE! '
        ENDIF

        TSEC_AVERAGE = 86400D0*R_DAYS_AVERAGE
        R_TD_AVERAGE  = DT/TSEC_AVERAGE
        WRITE(NFL,*)'init_td_mean: time filtering scale is ',
     *  R_DAYS_AVERAGE,
     *  ' days, each step weight is ',R_TD_AVERAGE
        IF(PRESENT(ZN) .AND. PRESENT(UA) .AND. PRESENT(VA))THEN
          ZN_MEAN = ZN
          UV_MEAN(:,:,1) = UA(:,:)
          UV_MEAN(:,:,2) = VA(:,:)
          WRITE(NFL,*)'*** Initialized tide mean field'
        ELSE
          ZN_MEAN = 0.0
          UV_MEAN = 0.0
          WRITE(NFL,*)'*** Nullified tide mean field'
        ENDIF
        IF(PRESENT(T) .AND. PRESENT(S))THEN
          TTS_MEAN(:,:,:,1) = T(:,:,:)
          TTS_MEAN(:,:,:,2) = S(:,:,:)
          WRITE(NFL,*)'*** Initialized 3D mean TS fields'
        ELSE
          TTS_MEAN = 0.0
          WRITE(NFL,*)'*** Nullified 3D mean TS fields'
        ENDIF
        IF(PRESENT(U) .AND. PRESENT(V))THEN
          TUV_MEAN(:,:,:,1) = U(:,:,:)
          TUV_MEAN(:,:,:,2) = V(:,:,:)
          WRITE(NFL,*)'*** Initialized 3D mean UV fields'
        ELSE
          TUV_MEAN = 0.0
          WRITE(NFL,*)'*** Nullified 3D mean UV fields'
        ENDIF
      end subroutine init_td_mean

      subroutine update_td_mean(UTC_TM_JDS,NFL,ZN,UA,VA,T,S,U,V)
        use time_serv
        implicit none
        REAL,INTENT(IN):: ZN(IM,JM),UA(IM,JM),VA(IM,JM)
        REAL,INTENT(IN),OPTIONAL:: T(IM,JM,KM),S(IM,JM,KM)
        REAL,INTENT(IN),OPTIONAL:: U(IM,JM,KM),V(IM,JM,KM)
        Integer,intent(in):: NFL
        DOUBLE PRECISION, INTENT(IN ):: UTC_TM_JDS

        DOUBLE PRECISION  UTC_TIME, DRAMP
        integer iret
        real rb

        CALL CNV_TIME(UTC_TM_JDS,UTC_TIME,NFL,IRET)
        IF(IRET /= 0)THEN
          WRITE(NFL,*)'update_td_mean: TIDE AND POM CALENDARS DIVERGE!'
          RETURN
        ENDIF
        IF(TSEC_AVERAGE <= 0D0 .or. R_TD_AVERAGE <= 0.)then
          WRITE(NFL,*)'update_td_mean: ',
     *    'do not initialized filter parameters!'
          RETURN
        ENDIF
        DRAMP = UTC_TIME-UTC_AVERAGE_FROM
        RB=0
        IF(DRAMP < -2D0)THEN
          WRITE(NFL,*)'update_td_mean: ',
     *    'TIDE ESTIMATION TIME IS BEFORE THE INITIALIZATION TIME'
          RETURN
        ELSE IF(DRAMP < TSEC_AVERAGE)then
          RB = R_TD_AVERAGE*TSEC_AVERAGE/DRAMP
          if(RB > 1.0) RB = 1.
        else
          RB = R_TD_AVERAGE
        endif
        ZN_MEAN = (ZN - ZN_MEAN)*RB + ZN_MEAN
        UV_MEAN(:,:,1) = (UA(:,:)- UV_MEAN(:,:,1))*RB + UV_MEAN(:,:,1)
        UV_MEAN(:,:,2) = (VA(:,:)- UV_MEAN(:,:,2))*RB + UV_MEAN(:,:,2)
        IF(PRESENT(T) .AND. PRESENT(S))THEN
          TTS_MEAN(:,:,:,1) =
     *      (T(:,:,:)-TTS_MEAN(:,:,:,1))*RB + TTS_MEAN(:,:,:,1)
          TTS_MEAN(:,:,:,2) =
     *      (S(:,:,:)-TTS_MEAN(:,:,:,2))*RB + TTS_MEAN(:,:,:,2)
        ENDIF
        IF(PRESENT(U) .AND. PRESENT(V))THEN
          TUV_MEAN(:,:,:,1) =
     *      (U(:,:,:)-TUV_MEAN(:,:,:,1))*RB + TUV_MEAN(:,:,:,1)
          TUV_MEAN(:,:,:,2) =
     *      (V(:,:,:)-TUV_MEAN(:,:,:,2))*RB + TUV_MEAN(:,:,:,2)
        ENDIF
      end subroutine update_td_mean
C
C**********************************************************************C
C
      SUBROUTINE JPOM_INIT_TIDE(UTC_TM_JDS,IM_,JM_,KM_,XTGR,YTGR,
     *                          XUGR,YUGR,XVGR,YVGR,
     *                          NF_LG_,NF_IO,IRET,MT_FILE)
C**********************************************************************C
C
C     INIT - ESTIMATE THE TIME-INDEPENDENT AMPLITUDE AND PHASE OF USED TIDAL
C            HARMONICS IN ALL OPEN (TIDAL FORCED) BOUNDARY POINTS
C            AND FOR ALL MODEL DOMAIN;
C
C     IMPORTANT!
C
c     USE OCEAN TIDES WITH THE LOADING TIDES ESTIMATED
c     (NOT GEO-CENTRIC [*_gc.nao] AND
C     NOT RADIAL LOADING [*_rload.nao])
C
C     IRET = 0 - NORMAL EXECUTION, ANOTHER VALUE - ERROR
C
C     SERGEY VARLAMOV
C
C     JULY 2001
C
C**********************************************************************C

      USE COM_SERV, ONLY: TPOINT, RAD_IN_D, TDATE, TSEC, TDT
      use time_serv

!      USE       TIDE_SERV
!#ifdef USE_PORT_LIB
!      USE PORT_LIB_NAME, ONLY: FLUSH
!#endif
      IMPLICIT  NONE

      INTEGER,          INTENT(OUT):: IRET
      INTEGER,          INTENT(IN)::  IM_,JM_,NF_LG_,NF_IO,KM_
      DOUBLE PRECISION, INTENT(IN ):: XTGR(IM), YTGR(JM),
     *                                XUGR(IM), YUGR(JM),
     *                                XVGR(IM), YVGR(JM)
      DOUBLE PRECISION, INTENT(IN ):: UTC_TM_JDS
      CHARACTER*(*),INTENT(IN )::     MT_FILE
C
C ******** LOCAL SAVED VARIABLES::
C
C     ARRAYS FOR TIME-VARIABLE TIDAL PHASE AND AMPLITUDE FACTOR,
C     ARRAYS TO KEEP CONSTANT TIDAL HARMONICS IN ALL MODEL POINTS:
C
      INTEGER,SAVE::              NUM_TH = 0,
     *                            NUM_TPARAM = 1
!     *                            ,ISTART = 0
C
C ******** LOCAL AUTOMATIC VALUES::
C
      INTEGER                   I,J,M,IERR,
     *                          IT,IPAR,JT,ITIDE_TYPE
!,I1,I2,J1,J2,L,
      REAL                      HAMP, HPHASE, HAMPSIN, HAMPCOS
      TYPE (TPOINT)::           PTH
      REAL, PARAMETER::         EPS = 1.0E-5
      LOGICAL::                 INIT_MEAN, LT_BD, LT_FC, LT_UH,
     *                          MTF_EXIST=.FALSE.
      INTEGER                   IDEVRW, IMT, JMT, TIDE_MODEL
      CHARACTER(120)::          HRMNC_NAMES = " "
!
!     initialize from input and assume calling filter each EXT mode step
!
      NF_LOG = NF_LG_

      WRITE(NF_LOG,*)'*** L_TIDE_BD = ',L_TIDE_BD
      WRITE(NF_LOG,*)'*** L_TIDE_FC = ',L_TIDE_FC
      WRITE(NF_LOG,*)'*** L_TIDE_UVH = ',L_TIDE_UVH
      IF(R_TD_AVERAGE > 0.)THEN
        WRITE(NF_LOG,*)'*** Time filtering scale is ',
     *  R_DAYS_AVERAGE,
     *  ' days, each step weight is ',R_TD_AVERAGE
      ENDIF
C-----------------------------------------------------------------
C     Do not used temporaly, but present in the code
!      REAL      XUGR(1), YUGR(1)
C-----------------------------------------------------------------
      IRET = 0
      ZN_TIDE = 0
      UVHC_TIDE = 0
      ZN_EQTR = 0

      IF((.NOT. L_TIDE_BD) .AND. (.NOT. L_TIDE_FC))RETURN

      IF(IM_ /= IM .OR. JM_ /= JM  .OR. KM_ /= KM)THEN
        WRITE(NF_LOG,*)' *** ERROR (JPOM_INIT_TIDE): '//
     *  'IM AND JM DIMENSION VALUES DIFFER IN MODULE AND OUTSIDE: ',
     *  IM_,IM,JM_,JM,KM_,KM
        IRET = 1
        RETURN
      ENDIF

!      IF(ISTART .NE. 0)THEN
!        IRET = -1
!        WRITE(NF_LOG,*)' *** WARNING: MULTIPLE INITIALIZATION '//
!     *    'OF TIDE PARAMETERS, DO NOT RECOMMENDED!'//
!     *    ' N_INI = ',ISTART + 1
C
C       IF NEED MULTIPLE INITIALIZATION (NEXT DELETE PREVIOUS!)
C       MANUALLY FREE MEMORY BY:
C       CALL FREE_TIDES
C       IF IT WILL NOT BE DONE - NEXT CALL TO INIT_TIDE_PARAM
C       WILL RETURN AN ERROR...
C
C       DO NOT STOP HERE, CONTINUE EXECUTION
C
!      ENDIF
!      ISTART = ISTART + 1

      CALL CNV_TIME(UTC_TM_JDS,UTC_TIME_TIDES_FROM,NF_LOG,IRET)
      IF(IRET /= 0)THEN
        IRET = 5
        WRITE(NF_LOG,*)' *** TIDE AND POM CALENDARS DIVERGE! '
      ENDIF
      TSEC_DAMP = R_DAYS_DAMP*86400D0
      WRITE(NF_LOG,*)'*** Fixed damping for ',R_DAYS_DAMP,
     *  ' days or ',TSEC_DAMP,' secs'
!
!     DEFAULT MODEL is NAO-JAPAN
!
      if(SEA_TIDE%TIDE_MODEL<=0)then
        SEA_TIDE%TIDE_MODEL = TIDE_NAO99_JAPAN_SSH
      endif
C-----------------------------------------------------------
C     (1.2) define what tidal harmonics to use
C-----------------------------------------------------------

!      CALL ALLOC_TIDE_HARMONICS(SEA_TIDE, NUM_TH, IERR)
!      AND SPECIFY NECESSARY BY NAMES EXACTLY AS IN TIDE_SERV MODULE
!      or
!      SEA_TIDE%NUM_TH = 0  ! ALL AVAILABLE WILL BE ALLOCATED
!
      INQUIRE(EXIST=MTF_EXIST, FILE=MT_FILE)

      IF(MTF_EXIST)THEN
        CALL OPEN_FILE(IDEVRW,MT_FILE,'UNFORMATTED READ SEQUENTIAL')
        READ(IDEVRW,IOSTAT=M)LT_BD,LT_FC,LT_UH,
     *    TIDE_MODEL, NUM_TH, NUM_TPARAM, IMT, JMT, I

        IF( M /= 0 .or.
     *     (.not. LT_BD .and. L_TIDE_BD) .or.           !no ssh info
     *     (.not. LT_FC .and. L_TIDE_FC) .or.           !no h_eqv info
     *     (.not. LT_UH .and. L_TIDE_UVH) .or.           !no uvh info
     *     (TIDE_MODEL /= SEA_TIDE%TIDE_MODEL) .OR.
     *     (IMT /= IM) .OR. (JMT /= JM) .OR. (I /= NHARM_ALL) .OR.
     *     ((SEA_TIDE%NUM_TH > 0) .AND. (NUM_TH > SEA_TIDE%NUM_TH)))THEN

          WRITE(NF_LOG,*) 'BINARY TIDAL INPUT FILE EXISTS',
     *    ' BUT DO NOT CONFORM TO MODEL CONFIGURATION.'
          WRITE(NF_LOG,*) ' WILL BE RE-CREATED WITH THE SAME NAME: ',
     *    TRIM(MT_FILE)
          WRITE(NF_LOG,*)'MODEL: ',IM,JM,NHARM_ALL,
     *      SEA_TIDE%NUM_TH,SEA_TIDE%TIDE_MODEL
          WRITE(NF_LOG,*)'FILE:  ',IMT,JMT,I,NUM_TH,
     *      TIDE_MODEL,NUM_TPARAM

          MTF_EXIST = .FALSE.
          CALL RELEASE_UNIT(IDEVRW)
        ELSE
          SEA_TIDE%NUM_TH = NUM_TH
          READ(IDEVRW)NUM_IN_TMOD,EQV_TIDE
!          next is not necessary also could be done...
!          CALL ALLOC_TIDE_HARMONICS(SEA_TIDE, NUM_TH, IERR)
!          IF(IERR /= 0)THEN
!            WRITE(NF_LOG,*)' *** ERROR (MODULE JPOM_TIDE): '//
!     *      'ALLOCATE TIDE HARMONIC NAMES'
!            IRET = 1
!            RETURN
!          ENDIF
!          SEA_TIDE%HARMONIC_NAME(1:NUM_TH) = HARMONIC_NAME(1:NUM_TH)
          WRITE(NF_LOG,*)'READ FROM ',trim(MT_FILE),':'
          WRITE(NF_LOG,*) LT_BD,LT_FC,LT_UH,
     *                    TIDE_MODEL,NUM_TH,NUM_TPARAM,
     *                    IMT,JMT,I
          WRITE(NF_LOG,*)'NUM_IN_MODULE: ',
     *     (NUM_IN_TMOD(i),i=1,NHARM_ALL)
          WRITE(NF_LOG,*)'EQV TIDE: ', (EQV_TIDE(i),i=1,NHARM_ALL)
          if(NUM_TH > 0 .and. NUM_TH <= NHARM_ALL)then
            WRITE(NF_LOG,*)'GOT SEA LEVEL/SEA CURRENT HARMONICS: ',
     *      ((TW_NAMES(IABS(NUM_IN_TMOD(i)))//" "),i=1,NUM_TH)
          endif
          j = 0
          HRMNC_NAMES = "NO"
          do i = 1,NHARM_ALL
            if( EQV_TIDE(i) > 0 )then
               write(HRMNC_NAMES(j*5+1:),*) TW_NAMES(i)
               j = j + 1
            endif
          enddo
          if(j > 0) WRITE(NF_LOG,*)'BODY FORCE FIXED TO: ',
     *              trim(HRMNC_NAMES)
        ENDIF
      ENDIF

      IF(L_TIDE_FC)THEN  !NEED EQUIVALENT TIDE HEIGHT
        IF(.not. MTF_EXIST)THEN
!         if user defined limited number of harmonics by names
!         take this only for the astronomic body force
          CALL INIT_EQV_TIDE(SEA_TIDE, IERR, LOG=NF_LOG)
          IF(IERR /= 0)THEN
            IRET = 5
            WRITE(NF_LOG,*)
     *      ' *** ERROR INITIALIZING EQV TIDE INFO',IERR
            RETURN
          ENDIF
        ENDIF
        ALLOCATE(
     *    FS_LON(0:2,IM),
     *    FC_LON(0:2,IM),
     *    FLAT(NHARM_ALL,JM),
     *    STAT = IERR)

        IF(IERR.NE.0)THEN
          IRET = 6
          WRITE(NF_LOG,*)
     *    ' *** ERROR(1) ALLOCATING MEMORY IN INIT_TIDE '
          RETURN
        ENDIF

        FS_LON = 0.
        FC_LON = 0.
        FLAT = 0.

        DO M = 1,NHARM_ALL
          IF(EQV_TIDE(M) == 0)CYCLE
          ITIDE_TYPE = NINT(TIDE_TYPE(M))
          IF(ITIDE_TYPE == 1)THEN
            !DUIRNAL
            FLAT(M,1:JM) = H_EQ_CM(M)*TP_ALPHA(M)*
     *                      SIN(YTGR(1:JM)*2.*RAD_IN_D)
          ELSE IF(ITIDE_TYPE == 2)THEN
            !SEMI-DUIRNAL
            FLAT(M,1:JM) = H_EQ_CM(M)*TP_ALPHA(M)*
     *                      COS(YTGR(1:JM)*RAD_IN_D)**2
          ELSE
            !LONG-PERIOD
            FLAT(M,1:JM) = H_EQ_CM(M)*TP_ALPHA(M)*
     *                      (1.5*COS(YTGR(1:JM)*RAD_IN_D)**2 - 1.)
          ENDIF
        ENDDO
        DO I = 1,IM
C         *** LONG-PERIOD ***
          FS_LON(0,I) = 0.
          FC_LON(0,I) = 1.
C         *** DUIRNAL ***
          FS_LON(1,I) = SIN(XTGR(I)*RAD_IN_D)
          FC_LON(1,I) = COS(XTGR(I)*RAD_IN_D)
C         *** SEMI-DUIRNAL ***
          FS_LON(2,I) = SIN(XTGR(I)*2.*RAD_IN_D)
          FC_LON(2,I) = COS(XTGR(I)*2.*RAD_IN_D)
        ENDDO
      ENDIF

      IF(.NOT. L_TIDE_BD)THEN
        IF(MTF_EXIST)THEN
          CALL RELEASE_UNIT(IDEVRW)
        ELSE
          CALL OPEN_FILE
     *      (IDEVRW,MT_FILE,'UNFORMATTED WRITE SEQUENTIAL')
          WRITE(IDEVRW)L_TIDE_BD,L_TIDE_FC,L_TIDE_UVH,
     *      SEA_TIDE%TIDE_MODEL,SEA_TIDE%NUM_TH,
     *      NUM_TPARAM,IM,JM,NHARM_ALL
!          HARMONIC_NAME(1:NUM_TH) = SEA_TIDE%HARMONIC_NAME(1:NUM_TH)
          WRITE(IDEVRW)NUM_IN_TMOD,EQV_TIDE
          CALL RELEASE_UNIT(IDEVRW)

          WRITE(NF_LOG,*)'WROTE TO ',trim(MT_FILE),':'
          WRITE(NF_LOG,*)L_TIDE_BD,L_TIDE_FC,L_TIDE_UVH,
     *                   SEA_TIDE%TIDE_MODEL,SEA_TIDE%NUM_TH,
     *                   NUM_TPARAM,IM,JM,NHARM_ALL
          WRITE(NF_LOG,*)'NUM_IN_MODULE: ', (NUM_IN_TMOD(i),i=1,NUM_TH)
          WRITE(NF_LOG,*)'EQV TIDE: ', (EQV_TIDE(i),i=1,NHARM_ALL)
!          WRITE(NF_LOG,*)'WROTE SEA LEVEL HARMONICS: ',
!     *       ((TW_NAMES(IABS(NUM_IN_TMOD(i)))//" "),i=1,NUM_TH)
          j = 0
          DO I = 1,NHARM_ALL
            IF( EQV_TIDE(i) > 0 )THEN
               WRITE(HRMNC_NAMES(j*5+1:),*) TW_NAMES(i)
               j = j + 1
            ENDIF
          ENDDO
!          WRITE(NF_LOG,*)"NO SEA LEVEL DATA WRITTEN!"
          IF(j > 0) WRITE(NF_LOG,*)'BODY FORCE FIXED TO: ',
     *              TRIM(HRMNC_NAMES)
        ENDIF
        RETURN
      ENDIF

      IF(MTF_EXIST)THEN
        ALLOCATE(
     *      TH_AS(IM,JM,NUM_TH,NUM_TPARAM),
     *      TH_AC(IM,JM,NUM_TH,NUM_TPARAM),
     *      STAT = IERR)
        IF(IERR.NE.0)THEN
            IRET = 6
            WRITE(NF_LOG,*)
     *      ' *** ERROR ALLOCATING MEMORY IN INIT_TIDE '
            RETURN
        ENDIF
        READ(IDEVRW)TH_AS,TH_AC,I
        IF(I /= 123456)THEN
            IRET = 6
            WRITE(NF_LOG,*)
     *      ' *** ERROR READING MODEL TIDAL HARMONICS FROM ',
     *      TRIM(MT_FILE)
            RETURN
        ENDIF
        CALL RELEASE_UNIT(IDEVRW)
        RETURN
      ENDIF
!
!     First time initialization
!
!     NPT = MAX(IM+1,JM+1)
      NUM_TPARAM = 1
      IF((SEA_TIDE%TIDE_MODEL==TIDE_NAO99_JAPAN_SSH_UVH .or.
     *   SEA_TIDE%TIDE_MODEL>=TIDE_OTIS_GLOBAL) .and. L_TIDE_UVH)then
        NUM_TPARAM = 3
      else
        L_TIDE_UVH=.false.
      endif
!
!     avoid using of all tide model grid for undefined region,
!     use OCM domain instead:
!
      IF(SEA_TIDE%NUM_REGIONS <= 0)THEN
        NULLIFY(SEA_TIDE%GRID)
        CALL ALLOC_TIDE_GRID(SEA_TIDE, 1, IERR)
        if(ierr /= 0)then
          WRITE(NF_LOG,*)' *** error (module jpom_tide): '//
     *    'allocate tide definition grid'
          IRET = 1
          return
        endif

        SEA_TIDE%GRID(1)%XMIN = XTGR(1) - 3*(XTGR(2)-XTGR(1))
        SEA_TIDE%GRID(1)%XMAX = XTGR(IM) + 3*(XTGR(IM)-XTGR(IM-1))
        SEA_TIDE%GRID(1)%YMIN = YTGR(1) - 3*(YTGR(2)-YTGR(1))
        SEA_TIDE%GRID(1)%YMAX = YTGR(JM) + 3*(YTGR(JM)-YTGR(JM-1))
      ENDIF
C-----------------------------------------------------------
C     (1.4) initialize fixed astronomic parameters
C-----------------------------------------------------------
      CALL INIT_TIDAL_EQ_PARAM(SEA_TIDE, IERR, LOG=NF_LOG)
C
C     NEGATIVE VALUE IS WARNING: MULTIPLE INITIALIZAION WAS DONE,
C     OR NO FILES FOUND PASS IT::
C
      IF(IERR > 0 .OR. (IERR < -1 .AND. L_TIDE_BD))THEN
        IRET = 6
        WRITE(NF_LOG,*)
     *  ' *** ERROR INITIALIZING TIDE INFO FROM INIT_TIDE ',
     *  IERR, L_TIDE_BD, L_TIDE_FC
        RETURN
      ENDIF
C-----------------------------------------------------------
C     (1.5) read tidal harmonics and
C           estimate interpolation polinoms                      !MOD
C-----------------------------------------------------------
      CALL INIT_SL_TIDAL_HARMONICS
     *( SEA_TIDE, IERR, IO=NF_IO, LOG=NF_LOG,
     *  DEBUG = TIDE_DEBUG_PRINT, VOLUME_FLUX = 1)
      IF(IERR .NE. 0)THEN
        IRET = 7
        WRITE(NF_LOG,*)
     *  ' *** ERROR IN INIT_SL_TIDAL_HARMONICS: ', IERR
        RETURN
      ENDIF

      NUM_TH = SEA_TIDE%NUM_TH

c      write(NF_LOG,*)"Tide dims: ",NUM_TPARAM,NUM_TH,NPT,NHARM_ALL
c      write(*,*)"tide dims: ",NUM_TPARAM,NUM_TH,NPT,NHARM_ALL
c      call flush(nf_log)

      ALLOCATE(
     *  TH_AS(IM,JM,NUM_TH,NUM_TPARAM),
     *  TH_AC(IM,JM,NUM_TH,NUM_TPARAM),
     *  STAT = IERR)

      IF(IERR.NE.0)THEN
        IRET = 6
        WRITE(NF_LOG,*)
     *  ' *** ERROR ALLOCATING MEMORY IN INIT_TIDE '
        RETURN
      ENDIF

      TH_AS = 0.
      TH_AC = 0.
C-----------------------------------------------------------
C     DEFINE TIDAL HARMONICS IN ALL MODEL POINTS
C
        DO M=1,NUM_TH
C-----------------------------------------------------------
!          DO NBD = 1, NUM_OBSEG
!            IF(OBSEG(NBD)%NPP > 0 .AND. OBSEG(NBD)%TIDE)THEN
              DO IPAR=1,NUM_TPARAM
                DO J = 1,JM
                  IF(IPAR==1)THEN       !SEA LEVEL
                    PTH%Y  = YTGR(J)
                  ELSE IF(IPAR==2)THEN  !U velocity component
                    PTH%Y  = YUGR(J)
                  else                  !V velocity component
                    PTH%Y  = YVGR(J)
                  ENDIF
                  DO I = 1,IM
                    IF(IPAR==1)THEN      !SEA LEVEL
                      PTH%X = XTGR(I)
                    ELSE IF(IPAR==2)THEN !U velocity component
                      PTH%X = XUGR(I)
                    ELSE                 !V velocity component
                      PTH%X = XVGR(I)
                    ENDIF

                    CALL GET_SL_TIDAL_HARMONICS
     *              ( M, SEA_TIDE, PTH, HAMP, HPHASE,
     *                HAMPSIN, HAMPCOS,
     *                IERR, LOG=NF_LOG, CUT_EXTREMUMS=0,PARAM=IPAR)

                    IF(IERR .NE. 0)THEN
                      IRET = 3
                      WRITE(NF_LOG,*)
     *                ' *** ERROR IN GET_SL_TIDAL_HARMONICS FOR ',
     *                SEA_TIDE%HARMONIC_NAME(M), IERR
                      RETURN
                    ENDIF
C
C                   SPECIAL INDEX IS POSSIBLE BECAUSE ONE CYCLE
C                   IS ALWAYS "ZERO" LENGTH::
C
!                    L=J-J1+I-I1+1
!                    TH_AMP(L,M,NBD,IPAR) = HAMP
!                    TH_PHS(L,M,NBD,IPAR) = HPHASE
                    TH_AS(I,J,M,IPAR) = HAMPSIN
                    TH_AC(I,J,M,IPAR) = HAMPCOS
                  ENDDO
                ENDDO
              ENDDO
!            ENDIF
!          ENDDO
C-----------------------------------------------------------
        ENDDO !CYCLE M ENDS (HARMONIC)
C
C      DEBUG PRINT
C      
      IF(TIDE_DEBUG_PRINT > 0)THEN
        WRITE(NF_LOG,*)' NON-ZERO TIDAL HARMONICS DEFINED:'
        WRITE(NF_LOG,*)'   I   J    AMP  PHASE'
        DO M=1,NUM_TH
          WRITE(NF_LOG,*)'                      ',
     *    SEA_TIDE%HARMONIC_NAME(M)
!          DO J = 1, JM
!            DO I = 1, IM
!C              IF(TH_AMP(I,J,M,1) < EPS)CYCLE
!C              WRITE(NF_LOG,'(2I4,3(F12.2,F8.2))') I,J,
!C     *        (TH_AMP(I,J,M,IPAR),TH_PHS(I,J,M,IPAR),
!C     *        IPAR=1,NUM_TPARAM)
!              IF(TH_AS(I,J,M,1)*TH_AS(I,J,M,1)+
!     *           TH_AC(I,J,M,1)*TH_AC(I,J,M,1) < EPS)CYCLE
!              WRITE(NF_LOG,'(2I4,3(F12.2,F8.2))') I,J,
!     *        (TH_AS(I,J,M,IPAR),TH_AC(I,J,M,IPAR),
!     *        IPAR=1,NUM_TPARAM)
!            ENDDO
!          ENDDO
        ENDDO
      ENDIF
C
C     DO NOT NEED MORE ARRAYS OF INTERPOLATION POLINOMS:
C
!      IF(TIDE_DEBUG_PRINT == 0)
!     *  CALL FREE_TIDE_HARMONICS(SEA_TIDE,LOG=NF_LOG)
      CALL FREE_TIDE_HARMONICS(SEA_TIDE,LOG=NF_LOG)
!
!     save binary model tide version
!
      CALL OPEN_FILE
     *  (IDEVRW,MT_FILE,'UNFORMATTED WRITE SEQUENTIAL')
      WRITE(IDEVRW)L_TIDE_BD,L_TIDE_FC,L_TIDE_UVH,
     *  SEA_TIDE%TIDE_MODEL,SEA_TIDE%NUM_TH,
     *  NUM_TPARAM,IM,JM,NHARM_ALL
!      HARMONIC_NAME(1:NUM_TH) = SEA_TIDE%HARMONIC_NAME(1:NUM_TH)
      WRITE(IDEVRW)NUM_IN_TMOD,EQV_TIDE
      WRITE(IDEVRW)TH_AS,TH_AC,123456
      CALL RELEASE_UNIT(IDEVRW)

      WRITE(NF_LOG,*)'WROTE TO ',trim(MT_FILE),':'
      WRITE(NF_LOG,*)L_TIDE_BD,L_TIDE_FC,
     *               SEA_TIDE%TIDE_MODEL,SEA_TIDE%NUM_TH,
     *               NUM_TPARAM,IM,JM,NHARM_ALL
      WRITE(NF_LOG,*)'NUM_IN_MODULE: ', (NUM_IN_TMOD(i),i=1,NUM_TH)
      WRITE(NF_LOG,*)'EQV TIDE: ', (EQV_TIDE(i),i=1,NHARM_ALL)
      WRITE(NF_LOG,*)'WROTE SEA LEVEL HARMONICS: ',
     *  ((TW_NAMES(IABS(NUM_IN_TMOD(i)))//" "),i=1,NUM_TH)
      j = 0
      DO I = 1,NHARM_ALL
        IF( EQV_TIDE(I) > 0 )THEN
           WRITE(HRMNC_NAMES(j*5+1:),*) TW_NAMES(i)
           j = j + 1
        ENDIF
      ENDDO
      IF(j > 0) WRITE(NF_LOG,*)'BODY FORCE FIXED TO: ',
     *          TRIM(HRMNC_NAMES)

      END SUBROUTINE JPOM_INIT_TIDE

C**********************************************************************C
      SUBROUTINE GET_TIDE(UTC_TM_JDS)
C**********************************************************************C
C     USE GET_TIDE PARAM TO ESTIMATE TIME-DEPENDANT FACTORS
C     OF TIDAL EQUATIONS FOR ALL HARMONICS AT THE UTC_TIME
C     AND CALCULATE RESULTING TOTAL TIDAL SEA LEVEL IN
C     ALL OPEN (TIDAL FORCED) BOUNDARY POINTS AND FOR MODEL DOMAIN

      USE       COM_SERV, ONLY: RAD_IN_D
      use time_serv
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN ):: UTC_TM_JDS

      INTEGER  IRET,NUM_TH,IPAR,I,J,M,JT,ITIDE_TYPE,NUM_TPARAM
      REAL,SAVE::
     *         F(NHARM_ALL)=   (/(0,i=1,NHARM_ALL)/),
     *         TDPH(NHARM_ALL)=(/(0,i=1,NHARM_ALL)/),
     *         FA(NHARM_ALL)=  (/(0,i=1,NHARM_ALL)/),
     *         PHA(NHARM_ALL)= (/(0,i=1,NHARM_ALL)/)
      REAL     CTDPH(NHARM_ALL), STDPH(NHARM_ALL),RF(NHARM_ALL)
      DOUBLE PRECISION::
     *         GLB_UTC_TIME, DRAMP
      REAL     TRAMP,H_TMP,FF,PP,COS_PH,SIN_PH,TD_BETA
!     *         ,RB
C
C     INITIALIZE ARRAYS BY ZEROES::
C
C     skip all errors analysis
C
      CALL CNV_TIME(UTC_TM_JDS,GLB_UTC_TIME,NF_LOG,IRET)
      DRAMP = GLB_UTC_TIME-UTC_TIME_TIDES_FROM

!      write(NF_LOG,*)'Test get_tide, UTC_TM_JDS=',
!     *  UTC_TM_JDS,DRAMP,TSEC_DAMP,DRAMP/TSEC_DAMP
!
!     return if called before tides are declared to start in initialization
!
      IF(DRAMP <= 0D0)RETURN
      TRAMP = 1E0
      IF(DRAMP < TSEC_DAMP) TRAMP = DRAMP/TSEC_DAMP

      CALL GET_TIDAL_EQ_PARAM(SEA_TIDE, GLB_UTC_TIME,
     *                        TDPH, F, PHA, FA)

      IF(L_TIDE_BD)THEN
        NUM_TH = SEA_TIDE%NUM_TH
        RF = F*TRAMP*CM_TO_USER_UNITS

        DO M = 1,NUM_TH
          CTDPH(M) = COS(TDPH(M)*RAD_IN_D)
          STDPH(M) = SIN(TDPH(M)*RAD_IN_D)
        ENDDO

        if(L_TIDE_UVH)THEN
          NUM_TPARAM = 3
        else
          NUM_TPARAM = 1
        endif

        DO IPAR=1,NUM_TPARAM
          IF(IPAR == 2)RF = RF*CM_TO_USER_UNITS  !UVH needs convert both UV and H!
          DO J = 1,JM
            DO I = 1,IM
              H_TMP = 0.
              DO M = 1,NUM_TH
                H_TMP = H_TMP +
     *          RF(M)*(TH_AS(I,J,M,IPAR)*STDPH(M)+
     *                 TH_AC(I,J,M,IPAR)*CTDPH(M))
              ENDDO
              IF(IPAR==1)THEN  !SEA LEVEL
                ZN_TIDE(I,J)  = H_TMP
              ELSE
                UVHC_TIDE(I,J,IPAR-1)  = H_TMP
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      IF(L_TIDE_FC)THEN
C-----------------------------------------------------------
C      equlibrium tide sea level anomaly in point *
C-----------------------------------------------------------
        RF = FA*TRAMP*CM_TO_USER_UNITS

        ZN_EQTR = 0.
        JT=0
        DO M = 1,NHARM_ALL
          IF(EQV_TIDE(M) == 1)THEN
            TD_BETA = TP_BETA(M)
            JT = JT + 1
            FF = RF(M)
            PP = PHA(M)
            ITIDE_TYPE = NINT(TIDE_TYPE(M))

            COS_PH = COS(PP*RAD_IN_D)
            SIN_PH = SIN(PP*RAD_IN_D)

            DO J=1,JM
              DO I=1,IM
                ZN_EQTR(I,J) = ZN_EQTR(I,J) +
     *          FF*FLAT(M,J)*(COS_PH*FC_LON(ITIDE_TYPE,I)-
     *                        SIN_PH*FS_LON(ITIDE_TYPE,I))
              ENDDO
            ENDDO
          ENDIF
        ENDDO
c
C       Using next ZN_EQTR please replace in external U-V equations
C       ZN to (ZN - ZN_EQTR) to simulate tidal forcing
C
!vsm: temporary disable feedback as unclear with the short waves...
!        if(present(zn))then
!          IF(JT > 1) TD_BETA = 0.95
!          RB = (1.-TD_BETA)*TRAMP
!          ZN_EQTR = ZN_EQTR + ZN*RB
!        endif
      ENDIF

!      write(NF_LOG,*)'Test get_tide, Z_EQV,Z_TIDE=',
!     *  maxval(ZN_EQTR),minval(ZN_EQTR),
!     *  maxval(ZN_TIDE),minval(ZN_TIDE)

      END SUBROUTINE GET_TIDE
C
C**********************************************************************C
C
      SUBROUTINE FREE_TIDES()
        IMPLICIT NONE
        integer IT,I
        IT = 0
        I = 0
        IF(ALLOCATED(TH_AS))DEALLOCATE(TH_AS,STAT=I)
        IF(I .NE. 0)THEN
          IT=IT+1
          I = 0
        ENDIF
        IF(ALLOCATED(TH_AC))DEALLOCATE(TH_AC,STAT=I)
        IF(I .NE. 0)THEN
          IT=IT+2
          I = 0
        ENDIF
        IF(ALLOCATED(FS_LON))DEALLOCATE(FS_LON,STAT=I)
        IF(I .NE. 0)THEN
          IT=IT+4
          I = 0
        ENDIF
        IF(ALLOCATED(FC_LON))DEALLOCATE(FC_LON,STAT=I)
        IF(I .NE. 0)THEN
          IT=IT+4
          I = 0
        ENDIF
        IF(ALLOCATED(FLAT))DEALLOCATE(FLAT,STAT=I)
        IF(I .NE. 0)THEN
          IT=IT+128
          I = 0
        ENDIF
        IF(IT .NE. 0)THEN
          IF(NF_LOG.GT.0) WRITE (UNIT=NF_LOG,FMT=*)
     *    ' *** DEALLOCATION ERROR IN FREE_TIDE:',IT
        ENDIF

        CALL FREE_TIDE(SEA_TIDE, LOG=NF_LOG)

      END SUBROUTINE FREE_TIDES
C
C**********************************************************************C
C
      subroutine write_tide_mean
     &( tide_mean_file,TIME_JDS_START,TIME_JDS_END,idev )
c
      implicit none
      integer,intent(in)::          idev
      double precision,intent(in):: TIME_JDS_START,TIME_JDS_END
      character*(*),intent(in)::    tide_mean_file
      integer,parameter:: chk=123456
c
      integer idevrw
      logical it_exists
!
!     Backup previous restart files if these exist
!     to avoid problems if process would be killed
!     when writing
!
      inquire(FILE=tide_mean_file, EXIST=it_exists)
      if(it_exists)then
        call system("mv "//trim(tide_mean_file)//" "//
     *              trim(tide_mean_file)//".bak")
      endif

      call open_file
     *  (idevrw,tide_mean_file,'unformatted write sequential')
      write(idevrw) TIME_JDS_START,TIME_JDS_END,im,jm
      write(idevrw) zn_mean,UV_mean,TTS_MEAN,TUV_MEAN,chk
      call release_unit(idevrw)
      write(idev,*) 'Wrote tide mean fields for timejd = ',
     *              TIME_JDS_END/86400D0
      return
      end subroutine write_tide_mean
C
C**********************************************************************C
C
      subroutine read_tide_mean
     &( tide_mean_file,TIME_JDS_START,TIME_JDS_END,idev )
c
      use time_serv
      implicit none
      integer,intent(in)::           idev
      double precision,intent(out):: TIME_JDS_START,TIME_JDS_END
      character*(*),intent(in)::     tide_mean_file
c
      integer idevrw
      integer i,j,chk

      call open_file
     & (idevrw,tide_mean_file,'unformatted read sequential check')

      read(idevrw,err=902) TIME_JDS_START,TIME_JDS_END,i,j
      if(i /= im .or. j /= jm )then
        write(idev,*)
     &  'read_tide_mean: wrong dimensions IM,JM in input file: ',
     *  i,j,im,jm,TIME_JDS_START/86400D0,TIME_JDS_END/86400D0
        goto 902
      endif

      read(idevrw,err=902)  zn_mean,UV_mean,TTS_MEAN,TUV_MEAN,chk
      if(chk /= 123456)then
        write(idev,*) 'read_tide_mean: wrong end-marker: ',chk
        goto 902
      endif
      call release_unit(idevrw)
!
!     transfer averaging start time to the module global variable
!
      CALL CNV_TIME(TIME_JDS_START,UTC_AVERAGE_FROM,idev,I)
c
      write(idev,*)
     * 'read_tide_mean: entered time-filtered fields for timejd=',
     * TIME_JDS_START/86400D0,TIME_JDS_END/86400D0
c
      return

 902  continue
      write(0,*)'$$$ tide mean fields read error ',trim(tide_mean_file)
      stop

      end subroutine read_tide_mean

C**********************************************************************C
      SUBROUTINE GET_TIDE2(UTC_TM_JDS1,UTC_TM_JDS2)
C**********************************************************************C
C     USE GET_TIDE PARAM TO ESTIMATE TIME-DEPENDANT FACTORS
C     OF TIDAL EQUATIONS FOR ALL HARMONICS AT THE UTC_TIME
C     AND CALCULATE RESULTING TOTAL TIDAL SEA LEVEL IN
C     ALL OPEN (TIDAL FORCED) BOUNDARY POINTS AND FOR MODEL DOMAIN

      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN):: UTC_TM_JDS1
      DOUBLE PRECISION, INTENT(IN), OPTIONAL:: UTC_TM_JDS2
      DOUBLE PRECISION,PARAMETER:: t_time_conf=2D0
      DOUBLE PRECISION,SAVE:: t_time_i(2) = (/-1d99,-1d99/)
      INTEGER,SAVE:: t_ind_i(2) = (/1,2/)
      REAL RRTide

        if(.not. present(UTC_TM_JDS2))then
!
!         Linear interpolate or,
!         if interpolation impossible - define for single fixed time
!
          if( UTC_TM_JDS1 < t_time_i(1)-t_time_conf .or.
     *        UTC_TM_JDS1 > t_time_i(2)+t_time_conf .or.
     *        t_time_i(1) < 0d0 .or. t_time_i(2) < 0d0 .or.
     *        t_time_i(2) - t_time_i(1) < t_time_conf)then
            call GET_TIDE(UTC_TM_JDS1)
            write(NF_LOG,*)"GET_TIDE2: GET instead of interpolation:",
     &      t_time_i,UTC_TM_JDS1
          else
            RRTide=(UTC_TM_JDS1-t_time_i(1))/(t_time_i(2)-t_time_i(1))
            ZN_TIDE(:,:)= ZN_TIDE_i(:,:,t_ind_i(1)) +
     *      RRTide*(ZN_TIDE_i(:,:,t_ind_i(2))-ZN_TIDE_i(:,:,t_ind_i(1)))
            ZN_EQTR(:,:)= ZN_EQTR_i(:,:,t_ind_i(1)) +
     *      RRTide*(ZN_EQTR_i(:,:,t_ind_i(2))-ZN_EQTR_i(:,:,t_ind_i(1)))
            if(L_TIDE_UVH)THEN
              UVHC_TIDE(:,:,:)= UVHC_TIDE_i(:,:,:,t_ind_i(1)) +
     *        RRTide*(UVHC_TIDE_i(:,:,:,t_ind_i(2))-
     *          UVHC_TIDE_i(:,:,:,t_ind_i(1)))
            endif
          endif
          return
        endif
!
!       Enter 1 or 2 cases - output field timing unclear - interpolate!
!
        call GET_TIDE(UTC_TM_JDS2)
        t_ind_i(1) = t_ind_i(2)
        t_ind_i(2) = 3-t_ind_i(2)
        t_time_i(1)=t_time_i(2)
        t_time_i(2) = UTC_TM_JDS2
        ZN_TIDE_i(:,:,t_ind_i(2))=ZN_TIDE(:,:)
        ZN_EQTR_i(:,:,t_ind_i(2))=ZN_EQTR(:,:)
        if(L_TIDE_UVH)THEN
          UVHC_TIDE_i(:,:,:,t_ind_i(2))=UVHC_TIDE(:,:,:)
        endif

        if(DABS(UTC_TM_JDS1-t_time_i(1)) >= t_time_conf)then
          call GET_TIDE(UTC_TM_JDS1)
          t_time_i(1) = UTC_TM_JDS1
          ZN_TIDE_i(:,:,t_ind_i(1))=ZN_TIDE(:,:)
          ZN_EQTR_i(:,:,t_ind_i(1))=ZN_EQTR(:,:)
          if(L_TIDE_UVH)THEN
            UVHC_TIDE_i(:,:,:,t_ind_i(1))=UVHC_TIDE(:,:,:)
          endif
          write(NF_LOG,*)"GET_TIDE2: GET instead of exchange:",
     &    t_time_i,UTC_TM_JDS1
        endif

      end subroutine GET_TIDE2

C**********************************************************************C
      end module jpom_tide
C**********************************************************************C

