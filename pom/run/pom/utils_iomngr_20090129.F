#include "sysdep.h"
#include "stdfio.h"
c-----------------------------------------------------------------------
c     utility package
c
c     require:
c       nothing
c
c     provide:
#if defined NOT_HAVE_LEN_TRIM
c       function len_trim
#endif
c       function find_substring
#if defined NOT_HAVE_ADJUSTL
c       subroutine char_adjustl
#endif
c       subroutine get_date_and_time
c       subroutine char_lowercase
c       subroutine char_uppercase
c       subroutine message_dump
c       subroutine string2int4
c       subroutine string2real4
#if ! defined __GNUC__
c       subroutine string2int8
c       subroutine string2real8
#endif
c       subroutine get_digit_of_integer
c       subroutine get_digit_of_real
c       subroutine get_format_int
c       subroutine get_format_real
c       function nearest_index
c       function heaviside
c
c     Note: meaning of continuation character is as follows
c       I: input
c       O: output
c       M: input/output (modified)
c       W: work area
c       &: no meaning (general continuation of statement)
c
c     Creator: Takashi Kagimoto (kagimoto@jamstec.go.jp)
c     Create:  2004/April/22
c     Modify:  2004/May/10
c        (add function, nearest_index)
c     Modify:  2004/May/14
c        (add function, heaviside)
c-----------------------------------------------------------------------

#if defined NOT_HAVE_LEN_TRIM
c-----------------------------------------------------------------------
c     name:        len_trim
c     description: return length of a given string after removing
c                 trailing blanks
c
c     N.B.: len_trim() is an intrinsic function in Fortran90 and also
c         an extension of FORTRAN77 in some compilers
c-----------------------------------------------------------------------
      function len_trim(
     I     string)
      implicit none
c---- arguments
      character*(*) string
c---- return value
      integer*4     len_trim
c---- local
      integer*4     i

      len_trim = 0
      do i = len(string), 1, -1
        if (string(i:i) .ne. ' ') then
          len_trim = i
          goto 999
        endif
      enddo
 999  continue

      return
      end
#endif
#if defined NOT_HAVE_ADJUSTL
      subroutine char_adjustl(
     M       str)
      character*(*) str
      integer*4     n, nst

      do n = 1, len(str)
        if (str(n:n) .ne. ' ') goto 10
      enddo
 10   continue
      nst = n
      do n = nst, len(str)
        str(n-nst+1:n-nst+1) = str(n:n)
      enddo

      return
      end
#endif
c-----------------------------------------------------------------------
c     name:        get_date_and_time
c     description: get the current date and time in the following form
c                    hr:mi:se mo/dy/yyyy
c                  where
c                    hr: hour   (2 digits)
c                    mi: minute (2 digits)
c                    se: second (2 digits)
c                    mo: month  (2 digits)
c                    dy: day    (2 digits)
c                    yyyy: year (4 digits)
c
c     N.B.: This routine highly depends upon the system as well as
c          compiler that you used, thus you MUST customize this for
c          the specific platform. Currently SGI MIPS Pro F77 and F90,
c          Sun Workshop and Intel Fortran Compiler are supported.
c     N.B.: If you use Intel Fortran Compiler ver.7.x, you need to
c          link this with portability library using -Vaxlib option.
c-----------------------------------------------------------------------
      subroutine get_date_and_time(
     O     cur_time)
      implicit none
c---- arguments
      character*(*) cur_time
c---- local
      integer*4     year, month, day, hour, minute, second
#if defined __INTEL_COMPILER
      integer*4     milisec
#endif
      character*128 fmt
!#if defined __sgi
!      character*10  cdate, ctime
!#elif defined __GNUC__
!      character*10  cdate
!      character*8   ctime
!#endif
!#if (defined __sgi && defined _LANGUAGE_FORTRAN77) || defined __GNUC__
#if (defined __sgi && defined _LANGUAGE_FORTRAN77)
      integer*4     n
      character*3   cmon(12)
      data          cmon   /'jan', 'feb', 'mar', 'apr', 'may', 'jun',
     &                      'jul', 'aug', 'sep', 'oct', 'nov', 'dec'/
#endif
#if defined __sun
      integer*4     iarray(3)
#endif
!vsm: default
      character*10  cdate, ctime

      fmt = '(i2.2,":",i2.2,":",i2.2,1x,i2.2,"/",i2.2,"/",i4.4)'

      if (len(cur_time) .lt. 20) then
         call message_dump('get_date_and_time', 'E',
     &        'Length of the string of the argument must be >= 20')
      endif

!#if defined __sgi && defined _LANGUAGE_FORTRAN90
!      call date_and_time(cdate, ctime)
!      read(cdate, '(i4,i2,i2)') year, month, day
!      read(ctime(1:6), '(i2,i2,i2)') hour, minute, second
!vsm #elif (defined __sgi && defined _LANGUAGE_FORTRAN77) || defined __GNUC__
#if (defined __sgi && defined _LANGUAGE_FORTRAN77)
      call date(cdate)
      call time(ctime)
      read(ctime, '(i2,1x,i2,1x,i2)') hour, minute, second
      read(cdate, '(i2,5x,i2)') day, year
      year = year + 2000
      call char_lowercase(cdate)
      do n = 1, 12
        if (cdate(4:6) .eq. cmon(n)) then
          month = n
          goto 99
        endif
      enddo
 99   continue
#elif defined __sun
      call idate(iarray)
      day   = iarray(1)
      month = iarray(2)
      year  = iarray(3)
      call itime(iarray)
      hour   = iarray(1)
      minute = iarray(2)
      second = iarray(3)
#elif defined __INTEL_COMPILER
      call getdat(year, month, day)
      call gettim(hour, minute, second, milisec)
#else
      call date_and_time(cdate, ctime)
      read(cdate, '(i4,i2,i2)') year, month, day
      read(ctime(1:6), '(i2,i2,i2)') hour, minute, second
#endif

      write(cur_time, fmt) hour, minute, second, month, day, year

      return
      end

c-----------------------------------------------------------------------
c     name:        char_lowercase
c     description: change a given string into its lower case
c-----------------------------------------------------------------------
      subroutine char_lowercase(
     M     string)
      implicit none
c---- arguments
      character*(*) string
c---- local
      integer*4   ica, icz, i, idx, ila, iskip
      character*1 char

      integer*4   ichar
      intrinsic   ichar
#if defined NOT_HAVE_LEN_TRIM
      integer*4   len_trim
      external    len_trim
#endif

      ica = ichar('A')
      ila = ichar('a')
      iskip = ila - ica
      icz = ichar('Z')
      do i = 1, len_trim(string)
        idx = ichar(string(i:i))
        if (ica .le. idx .and. idx .le. icz)
     &    string(i:i) = char(idx+iskip)
      enddo

      return
      end

c-----------------------------------------------------------------------
c     name:        char_uppercase
c     description: change a given string into its upper case
c-----------------------------------------------------------------------
      subroutine char_uppercase(
     M     string)
      implicit none
c---- arguments
      character*(*) string
c---- local
      integer*4   ila, ilz, i, idx, ica, iskip
      character*1 char

      integer*4   ichar
      intrinsic   ichar
#if defined NOT_HAVE_LEN_TRIM
      integer*4   len_trim
      external    len_trim
#endif

      ila = ichar('a')
      ica = ichar('A')
      iskip = ica - ila
      ilz = ichar('z')
      do i = 1, len_trim(string)
        idx = ichar(string(i:i))
        if (ila .le. idx .and. idx .le. ilz)
     &    string(i:i) = char(idx+iskip)
      enddo

      return
      end

c-----------------------------------------------------------------------
c     name:        message_dump
c     description: print out the message as standard error
c                 The first argument is either 'E' or 'W'.
c                  E: fatal error -> program will be stopped
c                  W: warning     -> program is continued
c-----------------------------------------------------------------------
      subroutine message_dump(
     I     caller, type, message)
      implicit none
c---- arguments
      character*1   type
      character*(*) caller, message
c---- local
      character*1   c

#if defined NOT_HAVE_TRIM
      character*512 string
# if defined NOT_HAVE_LEN_TRIM
      integer*4   len_trim
      external    len_trim
# endif
#endif

      c = type
      call char_lowercase(c)

      if (c .eq. 'e') then
#if defined NOT_HAVE_TRIM
        string = '*** ERROR ('
     &        // caller(1:len_trim(caller))
     &        // ') '
     &        // message(1:len_trim(message))
        write(STDERR,'(a)') string(1:len_trim(string))
#else
        write(STDERR,'(a)')
     &    '*** ERROR ('//trim(caller)//') '//trim(message)
#endif
        stop
      else if (c .eq. 'w') then
        write(STDERR,'(a)')
#if defined NOT_HAVE_TRIM
        string = '*** WARN ('
     &        // caller(1:len_trim(caller))
     &        // ') '
     &        // message(1:len_trim(message))
        write(STDERR,'(a)') string(1:len_trim(string))
#else
        write(STDERR,'(a)')
     &    '*** WARN ('//trim(caller)//') '//trim(message)
#endif
      endif

      return
      end

c-----------------------------------------------------------------------
c     name:        find_substring
c     description: if substring is a part of string, return TRUE
c                 otherwise return FALSE
c-----------------------------------------------------------------------
      function find_substring(
     I     string, substring)
      implicit none
c---- arguments
      character*(*) string, substring
c---- return value
      logical       find_substring
c---- local
#if defined NOT_HAVE_TRIM && defined NOT_HAVE_LEN_TRIM
      integer*4   len_trim
      external    len_trim
#endif

      find_substring = .false.
#if defined NOT_HAVE_TRIM
      if (index(string(1:len_trim(string)),
     &          substring(1:len_trim(substring))) .ne. 0) then
#else
      if (index(trim(string), trim(substring)) .ne. 0) then
#endif
        find_substring = .true.
      endif

      return
      end

c-----------------------------------------------------------------------
c     name:        string2int4
c     description: return single precision integer value which is
c                 read from the string
c-----------------------------------------------------------------------
      subroutine string2int4(
     I     string,
     O     ival)
      implicit none
c---- arguments
      character*(*) string
      integer*4     ival
c---- local
      integer*4     int4max
      parameter    (int4max  = z'7fffffff')
      character*128 ctmp
      integer*4     length_str, ndigit
      character*16  fmt
#if defined NOT_HAVE_LEN_TRIM
      integer*4     len_trim
      external      len_trim
#endif

#if defined NOT_HAVE_ADJUSTL
      ctmp = string
      call char_adjustl(ctmp)
#else
      ctmp = adjustl(string)
#endif
      length_str = len_trim(ctmp)
      if (length_str .eq. 0) then
        ival = int4max
        return
      endif

      call get_digit_of_integer(ctmp, length_str, ndigit)

      if (ndigit .eq. 0) then
        call message_dump('string2int4', 'E', 'Data scan is failed.')
      endif

      call get_format_int(ndigit, fmt)

      read(ctmp(1:ndigit), fmt) ival

      return
      end

c-----------------------------------------------------------------------
c     name:        string2real4
c     description: return single precision float value which is
c                 read from the string
c-----------------------------------------------------------------------
      subroutine string2real4(
     I     string,
     O     rval)
      implicit none
c---- arguments
      character*(*) string
      real*4        rval
c---- local
      real*4        real4max
      parameter    (real4max = z'7f7fffff')
      character*128 ctmp
      integer*4     length_str, ndigit_all, ndigit_dec, ndigit_exp
      integer*4     iexp
      character*16  fmt
#if defined NOT_HAVE_LEN_TRIM
      integer*4     len_trim
      external      len_trim
#endif

#if defined NOT_HAVE_ADJUSTL
      ctmp = string
      call char_adjustl(ctmp)
#else
      ctmp = adjustl(string)
#endif
      length_str = len_trim(ctmp)
      if (length_str .eq. 0) then
        rval = real4max
        return
      endif

      call get_digit_of_real(ctmp, length_str, ndigit_all, ndigit_dec,
     &     ndigit_exp)

      if (ndigit_all .eq. 0) then
        call message_dump('string2real4', 'E', 'Data scan is failed.')
      endif

      call get_format_real(ndigit_all, ndigit_dec, fmt)

      if (ndigit_exp .eq. 0) then
        read(ctmp(1:ndigit_all), fmt) rval
      else
        read(ctmp(1:ndigit_all-ndigit_exp-1), fmt) rval
        write(fmt,'("(i",i1,")")') ndigit_exp
        read(ctmp(ndigit_all-ndigit_exp+1:ndigit_all), fmt) iexp
        rval = rval * real(10.0)**(iexp)
      endif

      return
      end

#if ! defined __GNUC__
c-----------------------------------------------------------------------
c     name:        string2int8
c     description: return double precision integer value which is
c                 read from the string
c-----------------------------------------------------------------------
      subroutine string2int8(
     I     string,
     O     ival)
      implicit none
c---- arguments
      character*(*) string
      integer*8     ival
c---- local
      integer*8     int8max
      parameter    (int8max  = z'7fffffffffffffff')
      character*128 ctmp
      integer*4     length_str, ndigit
      character*16  fmt
# if defined NOT_HAVE_LEN_TRIM
      integer*4     len_trim
      external      len_trim
# endif

# if defined NOT_HAVE_ADJUSTL
      ctmp = string
      call char_adjustl(ctmp)
# else
      ctmp = adjustl(string)
# endif
      length_str = len_trim(ctmp)
      if (length_str .eq. 0) then
        ival = int8max
        return
      endif

      call get_digit_of_integer(ctmp, length_str, ndigit)

      if (ndigit .eq. 0) then
        call message_dump('string2int8', 'E', 'Data scan is failed.')
      endif

      call get_format_int(ndigit, fmt)

      read(ctmp(1:ndigit), fmt) ival

      return
      end

c-----------------------------------------------------------------------
c     name:        string2real8
c     description: return double precision float value which is
c                 read from the string
c-----------------------------------------------------------------------
      subroutine string2real8(
     I     string,
     O     rval)
      implicit none
c---- arguments
      character*(*) string
      real*8        rval
c---- local
      real*8        real8max
      parameter    (real8max = z'7fefffffffffffff')
      character*128 ctmp
      integer*4     length_str, ndigit_all, ndigit_dec, ndigit_exp
      integer*4     iexp
      character*16  fmt
# if defined NOT_HAVE_LEN_TRIM
      integer*4     len_trim
      external      len_trim
# endif

# if defined NOT_HAVE_ADJUSTL
      ctmp = string
      call char_adjustl(ctmp)
# else
      ctmp = adjustl(string)
# endif
      length_str = len_trim(ctmp)
      if (length_str .eq. 0) then
        rval = real8max
        return
      endif

      call get_digit_of_real(ctmp, length_str, ndigit_all, ndigit_dec,
     &     ndigit_exp)

      if (ndigit_all .eq. 0) then
        call message_dump('string2real8', 'E', 'Data scan is failed.')
      endif

      call get_format_real(ndigit_all, ndigit_dec, fmt)

      if (ndigit_exp .eq. 0) then
        read(ctmp(1:ndigit_all), fmt) rval
      else
        read(ctmp(1:ndigit_all-ndigit_exp-1), fmt) rval
        write(fmt,'("(i",i1,")")') ndigit_exp
        read(ctmp(ndigit_all-ndigit_exp+1:ndigit_all), fmt) iexp
        rval = rval * dble(10.0)**(iexp)
      endif

      return
      end
#endif

c-----------------------------------------------------------------------
c     name:        get_digit_of_integer
c     description: return the number of digits of an integer value in
c                 the string
c-----------------------------------------------------------------------
      subroutine get_digit_of_integer(
     I     string, length,
     O     ndigit)
      implicit none
c---- arguments
      character*(*) string
      integer*4     length
      integer*4     ndigit
c---- local
      integer*4     i, ic0, ic9
      logical       lsign

      lsign = .false.
      ndigit = 0
      ic0 = ichar('0')
      ic9 = ichar('9')

      do i = 1, length
        if (string(i:i) .eq. ' ') goto 99
        if (string(i:i) .eq. '-' .or. string(i:i) .eq. '+') then
          if (.not. lsign) then
            ndigit = ndigit + 1
            lsign  = .true.
            goto 10
          else
            call message_dump('get_digit_of_integer', 'E',
     &           'Encounter +/- sign more than twice.')
          endif
        endif
        if (ic0 .le. ichar(string(i:i)) .and.
     &      ichar(string(i:i)) .le. ic9) then
          ndigit = ndigit + 1
        endif
 10     continue
      enddo
 99   continue

      return
      end

c-----------------------------------------------------------------------
c     name:        get_digit_of_real
c     description: return the number of digits of an floating point
c                 value in the string
c-----------------------------------------------------------------------
      subroutine get_digit_of_real(
     I     string, length,
     O     ndigit_all, ndigit_dec, ndigit_exp)
      implicit none
c---- arguments
      character*(*) string
      integer*4     length
      integer*4     ndigit_all, ndigit_dec, ndigit_exp
c---- local
      integer*4     i, ic0, ic9
      logical       lsign, lexp

      lsign = .false.
      lexp  = .false.
      ndigit_all = 0
      ndigit_dec = 0
      ndigit_exp = 0
      ic0 = ichar('0')
      ic9 = ichar('9')

      do i = 1, length
        if (string(i:i) .eq. ' ') goto 99
        ndigit_all = ndigit_all + 1
        if (string(i:i) .eq. '.') then
          ndigit_dec = 0
          goto 10
        else if (string(i:i) .eq. 'e' .or. string(i:i) .eq. 'E' .or.
     &           string(i:i) .eq. 'd' .or. string(i:i) .eq. 'D') then
          lexp  = .true.
          lsign = .false.
          goto 10
        else if (string(i:i) .eq. '-' .or. string(i:i) .eq. '+') then
          if (.not. lsign) then
            if (.not. lexp) then
              ndigit_dec = ndigit_dec + 1
            else
              ndigit_exp = ndigit_exp + 1
            endif
            lsign = .true.
            goto 10
          else
            call message_dump('get_digit_of_real', 'E',
     &           'Encounter +/- sign more than twice.')
          endif
        else if (ic0 .le. ichar(string(i:i)) .and.
     &           ichar(string(i:i)) .le. ic9) then
          if (.not. lexp) then
            ndigit_dec = ndigit_dec + 1
          else
            ndigit_exp = ndigit_exp + 1
          endif
        else
          ndigit_all = 0
          goto 99
        endif
 10     continue
      enddo
 99   continue

      return
      end

c-----------------------------------------------------------------------
c     name:        get_format_int
c     description: return appropriate format to read an integer value
c                 from the string
c
c     N.B.: maximum integer value has at most 10 digits in single
c          precision and 19 digits in double precision. But the routine
c          does not take care of this in constructing the format. If
c          the number of digits exceeds 100, this will put the error
c          message into the standard error and stop the process.
c-----------------------------------------------------------------------
      subroutine get_format_int(
     I     ndigit,
     O     fmt)
      implicit none
c---- arguments
      integer*4     ndigit
      character*(*) fmt

      if (ndigit .lt. 10) then
        write(fmt,'("(i",i1,")")') ndigit
      else if (ndigit .lt. 100) then
        write(fmt,'("(i",i2,")")') ndigit
      else
        call message_dump('get_format_int', 'E',
     &       'The digit must be less than 100.')
      endif

      return
      end

c-----------------------------------------------------------------------
c     name:        get_format_real
c     description: return appropriate format to read a real value from
c                 the string
c
c     N.B.: Like as get_format_int(), the routine does not take care of
c          maximum digit of real value in constructing the format. If
c          the number of digits exceeds 100, this will put the error
c          message into the standard error and stop the process.
c-----------------------------------------------------------------------
      subroutine get_format_real(
     I     ndigit_all, ndigit_dec,
     O     fmt)
      implicit none
c---- arguments
      integer*4     ndigit_all, ndigit_dec
      character*(*) fmt

      if (ndigit_all .lt. 10 .and. ndigit_dec .lt. 10) then
        write(fmt,'("(f",i1,".",i1,")")') ndigit_all, ndigit_dec
      else if (ndigit_all .lt. 100 .and. ndigit_dec .lt. 10) then
        write(fmt,'("(f",i2,".",i1,")")') ndigit_all, ndigit_dec
      else if (ndigit_all .lt. 100 .and. ndigit_dec .lt. 100) then
        write(fmt,'("(f",i2,".",i2,")")') ndigit_all, ndigit_dec
      else
        call message_dump('get_format_real', 'E',
     &       'The digit must be less than 100.')
      endif

      return
      end

c-----------------------------------------------------------------------
c     name:        nearest_index
c     description: return the nearest index of array to a given value
c
c     N.B.: sequence of number given with an array must be ascending or
c         descending order. The order is specified with the fourth
c         argument of the function, 'ascend' or 'descend'. If the given
c         value is between two neighbouring numbers in the array, the
c         return value is selected according to the fifth argument,
c         'small' or 'large'. The former means that the selected index
c         is smaller than the other, and vice verse.
c     N.B.: If a given value is out-of-bounds, 0 is returned.
c-----------------------------------------------------------------------
      function nearest_index(
     I     val, seq, isize, order, side)
      implicit none
c---- arguments
      integer*4     isize
      real*4        val, seq(isize)
      character*(*) order, side
c---- return value
      integer*4     nearest_index
c---- local
      integer*4     i
      logical       lascend, lsmall
      character*16  string
c
c     check arguments (pass 1)
c
      string = order
      call char_lowercase(string)
      if (string .eq. 'ascend') then
        lascend = .true.
      else if (string .eq. 'descend') then
        lascend = .false.
      else
        call message_dump('nearest_index', 'E',
     &       'order must be ascend or descend.')
      endif
c
c     check arguments (pass 2)
c
      string = side
      call char_lowercase(string)
      if (string .eq. 'small') then
        lsmall = .true.
      else if (string .eq. 'large') then
        lsmall = .false.
      else
        call message_dump('nearest_index', 'E',
     &       'side must be small or large.')
      endif
c
c     look for the nearest index
c
      if (lascend) then
c------ <ascending order> ------
        if (seq(1)  .ge. seq(isize)) then
          call message_dump('nearest_index', 'E',
     &         'array is not arranged in the ascending order.')
        endif
c------ (out-of-bounds) ------
        if (val .lt. seq(1) .or. seq(isize) .lt. val) then
          nearest_index = 0
          return
        endif
c------ (in-bounds) ------
        if (lsmall) then
          do i = 1, isize-1
            if (seq(i) .le. val .and. val .lt. seq(i+1)) goto 10
          enddo
 10       continue
        else
          do i = 2, isize
            if (seq(i-1) .lt. val .and. val .le. seq(i)) goto 20
          enddo
 20       continue
        endif
      else
c------ <descending order> ------
        if (seq(1) .le. seq(isize)) then
          call message_dump('nearest_index', 'E',
     &         'array is not arranged in the descending order.')
        endif
c------ (out-of-bounds) ------
        if (val .gt. seq(1) .or. seq(isize) .gt. val) then
          nearest_index = 0
          return
        endif
c------ (in-bounds) ------
        if (lsmall) then
          do i = 1, isize-1
            if (seq(i) .ge. val .and. val .gt. seq(i+1)) goto 30
          enddo
 30       continue
        else
          do i = 2, isize
            if (seq(i-1) .gt. val .and. val .ge. seq(i)) goto 40
          enddo
 40       continue
        endif
      endif

      nearest_index = i

      return
      end

c-----------------------------------------------------------------------
c
c     Heaviside step function
c
c        heaviside = 1         if a >= 0
c        heaviside = 0         if a <  0
c
      function heaviside(
     I          a)
      implicit none
c---- arguments
      real*4      a
c---- return value
      integer*4   heaviside

      heaviside = int(0.5*(1.0+sign(1.0,a)))

      return
      end
