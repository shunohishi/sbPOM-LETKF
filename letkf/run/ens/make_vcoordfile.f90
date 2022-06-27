program main

  implicit none

  integer status,system

  integer imem,nmem
  integer inode,enode !each node for each job
  integer eproc !each processor for each job
  integer ipn,pn !Processor number
  integer inum

  character(100) dirname,filename
  character(3) mem
  
  !Read ensemble run information
  open(1,file="vcoordfile_info.txt",status="old")
  read(1,'(a)') dirname
  read(1,*) nmem
  read(1,*) enode
  read(1,*) eproc
  close(1)
  status=system("rm -f vcoordfile_info.txt")

  if(mod(eproc,enode) /= 0)then
     write(*,*) "***Error: mod(eproc,enode) /= 0"
     stop
  endif
  
  pn=eproc/enode
  inum=0

  do imem=1,nmem

     write(mem,'(i3.3)') imem
     filename=trim(dirname)//"/"//mem//"/vcoordfile"

     open(1,file=trim(filename),status="replace")
     do inode=1,enode
        do ipn=1,pn
           write(1,'(a1,i0,a1)') "(",inum,")"
        end do
        inum=inum+1
     end do
     close(1)

  end do

end program main
