
      module time_serv
        DOUBLE PRECISION,PRIVATE,SAVE::
     *              CAL_DIF=-208026662400D0
        LOGICAL,PRIVATE,SAVE:: FIRST_CALL = .true.
        contains
C**********************************************************************C
      DOUBLE PRECISION FUNCTION UTC_TM_FROM_JDS(UTC_TM_JDS)
        IMPLICIT  NONE
        DOUBLE PRECISION, INTENT(IN ):: UTC_TM_JDS
        UTC_TM_FROM_JDS = UTC_TM_JDS + CAL_DIF
      END FUNCTION

      DOUBLE PRECISION FUNCTION UTC_JDS_FROM_TM(UTC_TIME)
        IMPLICIT  NONE
        DOUBLE PRECISION, INTENT(IN ):: UTC_TIME
        UTC_JDS_FROM_TM = UTC_TIME - CAL_DIF
      END FUNCTION
C**********************************************************************C
C
      SUBROUTINE CAL_TIME(UTC_TM_JDS,Y,M,D,H,N,S)
C
C     CHANGE POM CALENDAR TIME TO TIDE INTERNAL CALENDAR TIME
C
C**********************************************************************C
      USE       COM_SERV, ONLY: TDATE, TSEC, TDT
      IMPLICIT  NONE
      DOUBLE PRECISION, INTENT(IN ):: UTC_TM_JDS
      INTEGER, INTENT(OUT)::          Y,M,D,H,N,S

      TYPE(TDATE)::                   UTC_DATE
      DOUBLE PRECISION::              UTC_TIME

      UTC_TIME = UTC_TM_JDS + CAL_DIF
      UTC_DATE = TDT(UTC_TIME)
      Y = UTC_DATE%IYEAR
      M = UTC_DATE%IMONTH
      D = UTC_DATE%IDAY
      H = UTC_DATE%IHOUR
      N = UTC_DATE%IMIN
      S = UTC_DATE%ISEC

      END SUBROUTINE CAL_TIME
C
C**********************************************************************C
C
      SUBROUTINE CNV_TIME(UTC_TM_JDS,UTC_TIME,NFL,IRET)
C
C     CHANGE POM CALENDAR TIME TO TIDE INTERNAL CALENDAR TIME
C
C**********************************************************************C
      USE       COM_SERV, ONLY: TDATE, TSEC, TDT
      use julian_calender, only: caldat
      IMPLICIT  NONE
      DOUBLE PRECISION, INTENT(IN ):: UTC_TM_JDS
      DOUBLE PRECISION, INTENT(OUT):: UTC_TIME
      INTEGER, INTENT(OUT)::          IRET
      INTEGER, INTENT(IN )::          NFL

      TYPE(TDATE)::                   UTC_DATE
      INTEGER I

      IRET = 0
!     DIFFERENCE BETWEEN THE UTC_TIME AND UTC_TM_JDS
!     MUST BE CONSTANT, ELSE DIFFERENT CALENDARS ARE USED
      IF(FIRST_CALL)THEN
        I = INT(UTC_TM_JDS/86400D0)
        UTC_DATE = TDATE(0,0,0,0,0,0)
        CALL CALDAT(I,UTC_DATE%IYEAR,UTC_DATE%IMONTH,UTC_DATE%IDAY)
        UTC_TIME = TSEC(UTC_DATE)+(UTC_TM_JDS-I*86400D0)
        IF(DABS(CAL_DIF - UTC_TIME + UTC_TM_JDS) > 2D0)THEN
          IRET = 1
          WRITE(NFL,*)' *** TIDE/METEO AND POM CALENDARS DIVERGE! ',
     *    CAL_DIF, UTC_TIME, UTC_TM_JDS, UTC_TIME-UTC_TM_JDS
        ENDIF
        CAL_DIF = UTC_TIME - UTC_TM_JDS
        FIRST_CALL = .false.
!        write(NFL,*)"TIDES INITIALIZED FROM ",tdt(UTC_TIME)
      ELSE
        UTC_TIME = UTC_TM_JDS + CAL_DIF
      ENDIF
      END SUBROUTINE CNV_TIME
C
C**********************************************************************C

      end module

