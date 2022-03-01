#include "sysdep.h"
#include "stdfio.h"
c-----------------------------------------------------------------------
c     I/O manager is a module-like group of routines
c
c     N.B.: the use of entry statement is old-fashioned, so that the
c          program structure should be change into that with module
c          written in Fortran90/95/2000 in the future, because the
c          entry statement is obsolete in such the standards.
c
c     require:
c        char_lowercase
c        find_substring
c        string2int4
c        message_dump
c
c     provide:
c        subroutine iomanager
c        entry init_iomngr
c        entry get_unitnumber
c        entry release_unit
c        entry query_ioinit
c        entry put_filename
c        subroutine open_file
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
c     Modify:  2004/June/09
c        Sun Workshop FORTRAN77 and SGI MIPS Pro FORTRAN77 does not
c       support 'position' item in open statement, and thus 'append'
c       mode should be specified with 'access' item.
c-----------------------------------------------------------------------
      subroutine iomanager
      implicit none
c---- arguments
      integer*4     fd
      logical       qinit
      character*(*) filename
c---- local
      integer*4     unit_start, unit_end, n
      parameter    (unit_start = 7, unit_end = 99)
c---- static
      logical       lunit_used(unit_start:unit_end)
      character*128 unit_name(unit_start:unit_end)
      logical       linit
      data          linit/ .false. /
      save          lunit_used, unit_name, linit

      return
c-----------------------------------------------------------------------
c     initialize I/O manager
c-----------------------------------------------------------------------
      entry init_iomngr

      do n = unit_start, unit_end
         lunit_used(n) = .false.
         unit_name(n)  = ''
      enddo
      linit = .true.

      return
c-----------------------------------------------------------------------
c     get unused unit number and set the corresponding flag bit to be
c     TRUE (1) such that the unit is shown to be currently used.
c-----------------------------------------------------------------------
      entry get_unitnumber(
     O     fd)

      do n = unit_start, unit_end
         if (.not. lunit_used(n)) goto 10
      enddo

      call message_dump('get_unitnumber', 'E',
     &     'No more unit is available.')

 10   continue
      fd = n
      lunit_used(fd) = .true.

      return
c-----------------------------------------------------------------------
c     save file name into the corresponding array
c-----------------------------------------------------------------------
      entry put_filename(
     I     filename, fd)

      unit_name(fd) = filename

      return
c-----------------------------------------------------------------------
c     close the specified unit and set the flag bit to be FALSE (0) such
c     that the unit is shown to be available.
c-----------------------------------------------------------------------
      entry release_unit(
     I     fd)

      if (.not. linit) then
        call message_dump('release_unit', 'E',
     &       'init_iomngr() is not called.')
      endif
      close(unit=fd)
      lunit_used(fd) = .false.
      unit_name(fd)  = ''

      return
c-----------------------------------------------------------------------
      entry query_ioinit(
     O     qinit)
      qinit = linit

      return
      end

c-----------------------------------------------------------------------
c     name:        open_file
c     description: open a file with specified options such as FORM,
c                 ACCESS, ACTION, RECL and so on.
c-----------------------------------------------------------------------
      subroutine open_file(
     O     fd,
     I     filename, options)
      implicit none
c---- argument
      integer*4,intent(out)::    fd
      character*(*),intent(in):: filename, options
c---- local
      character*256 string, message
      character*16  form, access, position, action, reclen_unit
      integer*4     reclen, wordlen, n1, n2
      logical       lexist, lcheck, linit
!vsm-s:
      integer io_stat
#if ! (defined __GNUC__ || (defined __sgi && defined _LANGUAGE_FORTRAN77) \
                        || (defined __sun && defined __SUNPRO_F77))
      real*4::      dummy=1.0
#endif
      logical       find_substring
      external      find_substring
#if defined NOT_HAVE_LEN_TRIM
      integer*4     len_trim
      external      len_trim
#endif
#if defined DEBUG
      character*128 fmt_dir, fmt_seq
      character*1   quote

      quote = char(ichar("'"))
      fmt_dir = '("open(unit=",i2,", file='//quote//'",a,"'//quote//
     &          ', form='//quote//'",a,"'//quote//
     &          ', access='//quote//'",a,"'//quote//
     &          ', recl=",XXX,", action='//
     &          quote//'",a,"'//quote//')")'
      fmt_seq = '("open(unit=",i2,", file='//quote//'",a,"'//quote//
     &          ', form='//quote//'",a,"'//quote//
     &          ', access='//quote//'",a,"'//quote//
     &          ', action='//quote//'",a,"'//quote//')")'
#endif

      call query_ioinit(linit)
      if (.not. linit) then
         call message_dump('open_file', 'E',
     &        'init_iomngr() is not called.')
      endif
c
c     default value
c
      form     = 'formatted'
      access   = 'sequential'
      position = 'rewind'
      action   = 'readwrite'
      lcheck   = .false.

#if defined __GNUC__ || (defined __sun && defined __SUNPRO_F77)
      wordlen = 4
#elif defined __sgi && defined _LANGUAGE_FORTRAN77
      wordlen = 1
#else
      inquire(iolength=wordlen) dummy
#endif
      if (wordlen .eq. 1) then
         reclen_unit = 'words'
      else
         reclen_unit = 'bytes'
      endif

      string = options
      call char_lowercase(string)
c
c     check options for open statement
c
      if (find_substring(string, 'append')) then
        position = 'append'
#if (defined __sun && defined __SUNPRO_F77) || (defined __sgi && defined _LANGUAGE_FORTRAN77)
        access='append'
#endif
      endif

      if (find_substring(string, 'unformatted')) then
c
c     binary file
c
        form = 'unformatted'

        if (find_substring(string, 'direct')) then
c
c       direct access format
c
          access = 'direct'

          n1 = index(string, "=") + 1
          do n2 = n1, len_trim(string)
            if (string(n2:n2) .ne. ' ') goto 10
          enddo
 10       continue
          n1 = n2
          do n2 = n1+1, len_trim(string)
            if (string(n2:n2) .eq. ' ') goto 20
          enddo
 20       continue
          n2 = n2-1
          call string2int4(string(n1:n2), reclen)

          if (find_substring(string, 'bytes')) then
            if (reclen_unit .eq. 'words') then
              reclen = reclen / wordlen
            endif
          else if (find_substring(string, 'words')) then
            if (reclen_unit .eq. 'bytes') then
              reclen = reclen * wordlen
            endif
          endif
        endif
      endif
c
c     strict specification of the action
c
      if (find_substring(string, 'readwrite')) then
        action = 'readwrite'
      else if (find_substring(string, 'read')) then
        action = 'read'
      else if (find_substring(string, 'write')) then
        action = 'write'
      endif
c
c     lcheck=.true.: check file existence
c
      if (find_substring(string, 'check')) then
        lcheck = .true.
      endif
c
c   get unit number
c
      call get_unitnumber(fd)
      call put_filename(filename, fd)

      if (lcheck) then
        inquire(exist=lexist, file=filename)
        if (.not. lexist) then
#if defined NOT_HAVE_TRIM
          message = 'File, '//filename(1:len_trim(filename))//
     &         ', is not found.'
#else
          message = 'File, '//trim(filename)//', is not found.'
#endif
          call message_dump('open_file', 'E', message)
        endif
      endif
      if (access .eq. 'direct') then
#if defined DEBUG
        n1 = index(fmt_dir, 'XXX')
        n2 = log10(real(reclen))+1
        if (n2 .lt. 10) then
          write(fmt_dir(n1:n1+2), '(a2,i1)') ' I',n2
        else if (n2 .lt. 100) then
          write(fmt_dir(n1:n1+2), '(a1,i2)') 'I',n2
        else
          call message_dump('open_file', 'E',
     &         'The digit of reclen must be less than 100.')
        endif
# if defined NOT_HAVE_TRIM
        write(STDOUT,fmt_dir) fd,
     &    filename(1:len_trim(filename)),
     &    form(1:len_trim(form)),
     &    access(1:len_trim(access)),
     &    reclen,
     &    action(1:len_trim(action))
# else
        write(STDOUT,fmt_dir) fd,
     &    trim(filename),
     &    trim(form),
     &    trim(access),
     &    reclen,
     &    trim(action)
# endif
#endif
#if defined __GNUC__ || (defined __sgi && defined _LANGUAGE_FORTRAN77)
        open(unit=fd, file=filename, form=form, access=access,
     &       recl=reclen, iostat=io_stat)
#else
        open(unit=fd, file=filename, form=form, access=access,
     &       recl=reclen, action=action, iostat=io_stat)
#endif
      else
#if defined DEBUG
# if defined NOT_HAVE_TRIM
        write(STDOUT,fmt_seq) fd,
     &    filename(1:len_trim(filename)),
     &    form(1:len_trim(form)),
     &    access(1:len_trim(access)),
     &    action(1:len_trim(action))
# else
        write(STDOUT,fmt_seq) fd,
     &    trim(filename),
     &    trim(form),
     &    trim(access),
     &    trim(action)
# endif
#endif
#if defined __GNUC__ || (defined __sgi && defined _LANGUAGE_FORTRAN77) \
                     || (defined __sun && defined __SUNPRO_F77)
        open(unit=fd, file=filename, form=form, access=access,
     *       iostat=io_stat)
#else
        open(unit=fd, file=filename, form=form, access=access,
     &       action=action, position=position, iostat=io_stat)
#endif
      endif
      if( io_stat .ne. 0) write(*,*)" Failed to open ",trim(filename)

      return
      end
