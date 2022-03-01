      program main

      implicit none

      integer iend, iint, julday_start, julday_time, julday
      integer iyear, imonth, iday, ihour, iprint
      real time0, time, dti
      character*4 cyear
      character*2 cmonth, cday, chour

      real*8 timeday

      iend = 100000
      dti = 0.5 * 20.
      iprint = int(3600./dti)
      time0 = 300.
      julday_start = julday(2012,06,01)

      do iint = 1, iend
          time = dti * float(iint) / 86400.e0 + time0 
          julday_time=julday_start+int(time)
!          julday_time=int(julday_start+time+0.001)
          call caldat(julday_time,iyear,imonth,iday)
!          ihour=nint((time-int(time))*24.)
          timeday=dble(dti)*dble(iint)/86400.d0+dble(time0)
          ihour=int((timeday+0.00001 - int(timeday))*24.)
          if(ihour.eq.24) ihour = 0
          write(cyear,'(i4.4)') iyear
          write(cmonth,'(i2.2)') imonth
          write(cday,'(i2.2)') iday
          write(chour,'(i2.2)') ihour
          if(mod(iint,iprint).eq.0)
     &     write(6,*) 
     &     iint, time, cyear//'/'//cmonth//'/'//cday//' '//chour
      end do

      stop
      end

!_______________________________________________________________________
      INTEGER FUNCTION JULDAY(IYYY,MONTH,DD)
C
C ********** SUBROUTINE DESCRIPTION:
C
C FINDS THE JULIAN DAY FROM A DATE.
C
C ********** ORIGINAL AUTHOR AND DATE:
C
C PRESS,FLANNERY,TEUKOLSKY,VETTERLING 1986.
C NUMERICAL RECIPES
C
C ********** REVISION HISTORY:
C
C
C ********** ARGUMENT DEFINITIONS:
C
      implicit none
      INTEGER IYYY,MONTH,DD
C
C NAME   IN/OUT DESCRIPTION
C
C IYYY     I    YEAR
C MONTH    I    MONTH (1 TO 12)
C DD       I    DAY OF MONTH
C JULDAY   O    JULIAN DAY
C
C ********** COMMON BLOCKS:
C
C NONE
C
C ********** LOCAL PARAMETER DEFINITIONS:
C
      INTEGER IGREG
      PARAMETER (IGREG = 15 + 31*(10 + 12*1582))
C
C ********** LOCAL VARIABLE DEFINITIONS:
C
      INTEGER JY,JM,JA
C
C NAME   DESCRIPTION
C
C
C ********** OTHER ROUTINES AND FUNCTIONS CALLED:
C
C INT    - INTRINSIC TRUNCATE
C
C---+67--1----+----2----+----3----+----4----+----5----+----6----+----7--
C
      IF (IYYY .LT. 0) IYYY = IYYY + 1
      IF (MONTH .GT. 2) THEN
        JY = IYYY
        JM = MONTH + 1
      ELSE
        JY = IYYY - 1
        JM = MONTH + 13
      ENDIF
      JULDAY = INT(365.25*JY) + INT(30.6001*JM) + DD + 1720995
      IF (DD + 31*(MONTH + 12*IYYY) .GE. IGREG) THEN
        JA = INT(0.01*JY)
        JULDAY = JULDAY + 2 - JA + INT(0.25*JA)
      ENDIF
      RETURN
      END

!_______________________________________________________________________
      SUBROUTINE CALDAT(JULIAN,IYYY,MONTH,DD)
C
C ********** SUBROUTINE DESCRIPTION:
C
C GIVEN THE JULIAN DAY, RETURNS THE YEAR, MONTH AND DAY OF MONTH.
C
C ********** ORIGINAL AUTHOR AND DATE:
C
C PRESS,FLANNERY,TEUKOLSKY,VETTERLING 1986.
C NUMERICAL RECIPES
C
C ********** REVISION HISTORY:
C
C
C ********** ARGUMENT DEFINITIONS:
C
      implicit none
      INTEGER JULIAN,IYYY,MONTH,DD
C
C NAME   IN/OUT DESCRIPTION
C
C JULIAN   I    THE JULIAN DAY
C IYYY     O    THE YEAR
C MONTH    O    THE MONTH (1 TO 12)
C DD       O    THE DAY OF THE MONTH
C
C ********** COMMON BLOCKS:
C
C NONE
C
C ********** LOCAL PARAMETER DEFINITIONS:
C
      INTEGER IGREG
      PARAMETER (IGREG=2299161)
C
C ********** LOCAL VARIABLE DEFINITIONS:
C
      INTEGER JALPHA,JA,JB,JC,JD,JE
C
C NAME   DESCRIPTION
C
C
C ********** OTHER ROUTINES AND FUNCTIONS CALLED:
C
C
C---+67--1----+----2----+----3----+----4----+----5----+----6----+----7--
C
C
      IF (JULIAN .GE. IGREG) THEN
        JALPHA = INT(((JULIAN - 1867216) - 0.25)/36524.25)
        JA = JULIAN + 1 + JALPHA - INT(0.25*JALPHA)
      ELSE
        JA = JULIAN
      ENDIF
      JB = JA + 1524
      JC = INT(6680. + ((JB - 2439870) - 122.1)/365.25)
      JD = 365*JC + INT(0.25*JC)
      JE = INT((JB - JD)/30.6001)
      DD = JB - JD - INT(30.6001*JE)
      MONTH = JE - 1
      IF (MONTH .GT. 12) MONTH = MONTH - 12
      IYYY = JC - 4715
      IF (MONTH .GT. 2) IYYY = IYYY - 1
      IF (IYYY .LE. 0) IYYY = IYYY - 1
      RETURN
      END
