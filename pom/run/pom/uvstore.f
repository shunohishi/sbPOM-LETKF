!_______________________________________________________________________
      subroutine uvstore

      implicit none
      include 'pom.h'

      integer iglb,jglb,iloc,jloc,nproc
      integer m,k,isec,iyear,imonth,iday,ihour,julday_time
      integer mstep
      data mstep/0/
      save mstep
      integer init
      data init/0/
      save init

      real ptdepth
      double precision timeday,timelocal

      character*4 cyear
      character*2 cmonth,cday,chour

      if(init.eq.0) then
        if(mstepmax.lt.iprint) then
          if(my_task.eq.master_task) then
            write(*,'(/a)') 'in uvstore: mstepmax < iprint'
            write(*,'(/a)') 'POM terminated with error'
          end if
          call finalize_mpi
          stop
        end if
        init=1
      end if

      mstep = mstep + 1

      iglb=341 ! 139.144 E izu12
      jglb=461 !  34.2546N
      nproc=int((iglb-2)/(im-2))+int((jglb-2)/(jm-2))*nproc_x

      if(my_task.eq.nproc) then 
        iloc=iglb-mod(nproc,nproc_x)*(im-2)
        jloc=jglb-(nproc/nproc_x)*(jm-2)
        do k=1,kbm1
          upt(k,mstep) = u(iloc,jloc,k)
          vpt(k,mstep) = v(iloc,jloc,k)
        end do   
        if(mstep.eq.iprint) then
!          julday_time=int(julday_start+time+0.001)
          julday_time=julday_start+int(time)
          call caldat(julday_time,iyear,imonth,iday)
!          timeday=dble(dti)*dble(iint)/86400.d0+dble(time0)
!          ihour=int((timeday+0.00001 - int(timeday))*24.)
          ihour=nint((time-int(time))*24.)
          if(ihour.eq.24) ihour = 0
          write(cyear,'(i4.4)') iyear
          write(cmonth,'(i2.2)') imonth
          write(cday,'(i2.2)') iday
          write(chour,'(i2.2)') ihour
          open(idevinout,file='UVPT_'//cyear//cmonth//cday//chour,
     $            form='formatted')
          do m=1,iprint
            isec=int(dti*m)
            do k=1,kbm1
              ptdepth=abs((h(iloc,jloc)+etf(iloc,jloc))*zz(k))
              write(idevinout,'((i5.5),x,i3.3,x,f7.2,2(x,f8.5))') 
     $          isec,k,ptdepth,upt(k,m),vpt(k,m)
            end do
          end do
          close(idevinout)
        end if
      end if

      if(mstep.eq.iprint) then
        mstep = 0
      end if

      return
      end
