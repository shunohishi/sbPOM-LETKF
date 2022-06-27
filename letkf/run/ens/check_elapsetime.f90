program main

  integer(8) itime

  !Get residual time [s]
  call chkelapse(itime)

  ![s] --> [min]
  itime=nint(dble(itime)/60.d0)
  
  open(1,file="check_elapsetime.dat",status="replace")
  write(1,'(i6)') itime
  close(1)

end program main

