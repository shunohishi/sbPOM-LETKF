!_______________________________________________________________________
      subroutine readtiming(timing,ratio,dayint_data,time0_data)
! determin read timing of forcing data
      implicit none
      include 'pom.h'
! arguments
      integer timing
      real ratio,dayint_data,time0_data
! local
      real dev,dif,timepre

      timepre=dti*float(iint-1)/86400.e0+time0
      dev=(timepre-time0_data)*86400.
     $    /(dayint_data*86400.)
      dif=dev-nint(dev)

      timing=0
      if(abs(dif)+1.e-6.lt.dti/86400./dayint_data) then
        timing=1
      end if

      ratio=dif
      if(abs(dif)+1.e-6.lt.dti/86400./dayint_data) then
        ratio=int(dif)
      else
        if(ratio.lt.0.) ratio=1.+ratio
      end if

      return
      end

!_______________________________________________________________________
      subroutine ratioclim(ratio,julday_time,time)
!     julday,time: [day]
      implicit none
! arguments
      integer julday_time
      real ratio,time
! local
      integer iyy0,iyy1,imm0,imm1,idd0,idd1
      integer julday0,julday1,julday
!
      call caldat(julday_time,iyy0,imm0,idd0)
      if(idd0.lt.15) then
        imm0=imm0-1
        if(imm0.lt.1) then
          imm0=12
          iyy0=iyy0-1
        end if
      end if
      imm1=imm0+1
      iyy1=iyy0
      if(imm1.gt.12) then
        imm1=1
        iyy1=iyy0+1
      end if
      julday0=julday(iyy0,imm0,15)
      julday1=julday(iyy1,imm1,15)
!
      ratio=(julday_time+time-int(time)-julday0)
     $      /(julday1-julday0)
!
      return
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
