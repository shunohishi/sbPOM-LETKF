
!#define TEST_CALDAT_D
!
      module julian_calender
        implicit none
        TYPE iJUL
              INTEGER:: YEAR=0
              INTEGER:: MONTH=1
              INTEGER:: DAY=1
              INTEGER:: HOUR=0
              INTEGER:: MIN=0
              INTEGER:: SEC=0
        END TYPE
        TYPE sJUL
              CHARACTER:: YEAR*4="0000"
              CHARACTER:: MONTH*2="01"
              CHARACTER:: DAY*2="01"
              CHARACTER:: HOUR*2="00"
              CHARACTER:: MIN*2="00"
              CHARACTER:: SEC*2="00"
        END TYPE
!vsm: add 0.5 secs in days for time adjustment to avoid rounding problem
        REAL,PARAMETER:: adj_isec = 0.5     !vsm: add 0.5 secs to adjust time
        DOUBLE PRECISION,PARAMETER:: adj_dsec = adj_isec/86400d0

      contains

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
      INTEGER,intent(in):: IYYY,MONTH,DD
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
      INTEGER IGREG,iYYYi
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
      iYYYi = iYYY
      IF (IYYY .LT. 0) IYYYi = IYYY + 1
      IF (MONTH .GT. 2) THEN
        JY = IYYYi
        JM = MONTH + 1
      ELSE
        JY = IYYYi - 1
        JM = MONTH + 13
      ENDIF
      JULDAY = INT(365.25*JY) + INT(30.6001*JM) + DD + 1720995
      IF (DD + 31*(MONTH + 12*IYYYi) .GE. IGREG) THEN
        JA = INT(0.01*JY)
        JULDAY = JULDAY + 2 - JA + INT(0.25*JA)
      ENDIF
      RETURN
      END FUNCTION

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
      INTEGER,intent(in):: JULIAN
      INTEGER,intent(out):: IYYY,MONTH,DD
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
      END SUBROUTINE

! vsm: Using real*4 for julday do not permit to recover even hours from it!
!
      SUBROUTINE CALDAT_D(JULIAN_D,IYYY,MONTH,DD,HH,MM,SEC)
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
      INTEGER,INTENT(OUT):: IYYY,MONTH,DD,HH,MM,SEC
      DOUBLE PRECISION, INTENT(IN):: JULIAN_D
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
      INTEGER, PARAMETER:: IGREG=2299161
C
C ********** LOCAL VARIABLE DEFINITIONS:
C
      INTEGER JALPHA,JA,JB,JC,JD,JE
      INTEGER JULIAN
      DOUBLE PRECISION JULIAN_RA, JULIAN_H
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
      JULIAN_RA = JULIAN_D + adj_dsec
      JULIAN = INT(JULIAN_RA)
      JULIAN_H = JULIAN_RA - JULIAN
      HH = int(JULIAN_H*24d0)
      MM = int((JULIAN_H-HH/24d0)*1440d0)
      SEC = int((JULIAN_H-HH/24d0-MM/1440d0)*86400d0)
!      write(*,*)"W: ",JULIAN_D,JULIAN_RA,JULIAN,JULIAN_H,HH,MM,SEC

      IF (JULIAN .GE. IGREG) THEN
        JALPHA = INT(((JULIAN - 1867216) - 0.25D0)/36524.25D0)
        JA = JULIAN + 1 + JALPHA - INT(0.25D0*JALPHA)
      ELSE
        JA = JULIAN
      ENDIF
      JB = JA + 1524
      JC = INT(6680D0 + ((JB - 2439870) - 122.1D0)/365.25D0)
      JD = 365*JC + INT(0.25D0*JC)
      JE = INT((JB - JD)/30.6001D0)
      DD = JB - JD - INT(30.6001D0*JE)
      MONTH = JE - 1
      IF (MONTH .GT. 12) MONTH = MONTH - 12
      IYYY = JC - 4715
      IF (MONTH .GT. 2) IYYY = IYYY - 1
      IF (IYYY .LE. 0) IYYY = IYYY - 1
      RETURN
      END SUBROUTINE

      SUBROUTINE CALDAT_IS(JULIAN_D,Ti,Ts)
C
        implicit none
        DOUBLE PRECISION, INTENT(IN):: JULIAN_D
        TYPE (iJUL),INTENT(OUT):: Ti
        TYPE (sJUL),OPTIONAL,INTENT(OUT):: Ts

        call CALDAT_D(JULIAN_D,Ti%YEAR,Ti%MONTH,Ti%DAY,
     *                         Ti%HOUR,Ti%MIN,Ti%SEC)
        if(present(Ts))then
           write(Ts%YEAR,'(i4.4)') Ti%YEAR
           write(Ts%MONTH,'(i2.2)') Ti%MONTH
           write(Ts%DAY,'(i2.2)') Ti%DAY
           write(Ts%HOUR,'(i2.2)') Ti%HOUR
           write(Ts%MIN,'(i2.2)') Ti%MIN
           write(Ts%SEC,'(i2.2)') Ti%SEC
        endif
      END SUBROUTINE
      END MODULE

#ifdef TEST_CALDAT_D
      program main
      use julian_calender
!        INTEGER:: JULDAY,IYYY=2007,MONTH=7,DD=30,HH=23,MM=59,SEC=59,JULDAYi
!        INTEGER:: JULDAY,IYYY=3007,MONTH=1,DD=1,HH=0,MM=0,SEC=0,JULDAYi
!        INTEGER:: JULDAY,IYYY=3009,MONTH=7,DD=30,HH=4,MM=25,SEC=59,JULDAYi
        INTEGER:: JULDAY,IYYY=-2009,MONTH=7,DD=30,HH=4,MM=25,SEC=59,
     *    JULDAYi
        DOUBLE PRECISION JULIAN_R
        type(iJUL)Ti
        type(sJUL)Ts
        character*(*),parameter:: date_frm =
     *'(a,i4.4,".",i2.2,".",i2.2," ",i2.2,":",i2.2,":",i2.2,'//
     *'1x,i8,1x,i8)'

        JULDAYi = JULDAY(IYYY,MONTH,DD)
        JULIAN_R = JULDAYi + (HH + (MM+SEC/60D0)/60D0)/24D0

        write(*,*)'get julday_i: ',
     *  JULDAYi,IYYY,MONTH,DD,HH,MM,SEC,JULIAN_R
        IYYY=0
        MONTH=0
        DD=0
        HH=03
        MM=0
        SEC=0
        call CALDAT_D(JULIAN_R,IYYY,MONTH,DD,HH,MM,SEC)

        write(*,*)'Inverted: ',JULDAYi,IYYY,MONTH,DD,HH,MM,SEC,JULIAN_R

        call CALDAT_IS(JULIAN_R,Ti)
        write(*,*)'Get: ',Ti
        call CALDAT_IS(JULIAN_R,Ti,Ts)
        write(*,*)'Get2: ',Ti
        write(*,*)'Get3: ',Ts

        write(*,date_frm)'Date: ',Ti,345,12345678

      end program

#endif
