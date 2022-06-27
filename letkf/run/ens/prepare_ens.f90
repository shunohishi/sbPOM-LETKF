program main

  implicit none

  integer iyr,imon,iday
  integer nmem
  character(100) dir

  open(1,file="ens_info.txt",status="old")
  read(1,'(a)') dir
  read(1,*) iyr,imon,iday
  read(1,*) nmem
  close(1)

  call io_time(dir,iyr,imon,iday,nmem)

end program main

!-------------------------------------------------------------
! Read & Write Time in 001
!-------------------------------------------------------------

subroutine io_time(dir,iyr,imon,iday,nmem)

  implicit none
  include "netcdf.inc"

  integer iyr,imon,iday
  integer imem,nmem

  integer status,access
  integer ncid,varid

  real time

  character(4) year
  character(3) mem
  character(2) month,day
  character(100) dir,filename
  
  write(year,'(i4.4)') iyr
  write(month,'(i2.2)') imon
  write(day,'(i2.2)') iday

  do imem=1,nmem

     write(mem,'(i3.3)') imem
     filename=trim(dir)//"/"//mem//"/restart."//year//month//day//".nc"
     status=access(trim(filename)," ")

     if(status /= 0)then
        write(*,*) "***Error: Not found "//trim(filename)
        stop
     end if

     if(imem == 1)then
        
        status=nf_open(trim(filename),nf_nowrite,ncid)
        
        status=nf_inq_varid(ncid,"time",varid)
        status=nf_get_var_real(ncid,varid,time)
        
        status=nf_close(ncid)

     else

        status=nf_open(trim(filename),nf_write,ncid)
        
        status=nf_inq_varid(ncid,"time",varid)
        status=nf_put_var_real(ncid,varid,time)
        
        status=nf_close(ncid)
                
     end if

  end do

end subroutine io_time
