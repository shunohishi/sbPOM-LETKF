! pom.f

! main program

      program pom
      implicit none
      include 'pom.h'
      integer ndd,nhh,nmm,nss
      double precision realtime,realtime0

      realtime0 = realtime()

! initialize model
      call initialize

! initialize tide
      if(ltide) call tide_initialize

! main loop
      do iint=1,iend
        if(ltide) call tide_advance ! update time dependent tidal parameters
        if(mod(iint,20).eq.0) call diag_w(6)
        call get_time
        call advance ! advance model
        if(error_status.ne.0) then
          if(my_task.eq.master_task) write(*,'(/a)')
     $                                       'POM terminated with error'
          call finalize_mpi
          stop
        end if
      end do

! computer time usage
      if (my_task.eq.master_task) then
        realtime0 = realtime() - realtime0
        ndd = realtime0/86400.
        nhh = (realtime0-ndd*86400.)/3600.
        nmm = (realtime0-ndd*86400.-nhh*3600.)/60.
        nss = realtime0-ndd*86400.-nhh*3600.-nmm*60.
        write (*,'(/a,i3,a,i2,a,i2,a,i2,a)') 'Elapsed time usage = ',
     &        ndd,' days ',nhh,' hours ',nmm,' minutes ',nss,' seconds'
        write (*,'(/a)') 'POM terminated successfully '
      end if

! finalize mpi
      call finalize_mpi

      if (my_task.eq.master_task) then
         write(6,*) "***normal end***"
      end if

      stop
      end

!_______________________________________________________________________
      real*8 function realtime()
      integer count, count_rate, count_max
      call system_clock(i,j,k)
      realtime=i/dble(j)
      return
      end
