!-------------------------------------------------------------------------
! Make Mean & Spread |
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
! Created by S.Ohishi 2019.10
! Modified by S.Ohishi 2020.04 for nudging term
!-------------------------------------------------------------------------
program main

  implicit none

  integer iyr,imon,iday,nday
  integer syr,smon,sday
  integer im,jm,km
  integer imem,nmem
  integer iswitch_budget
  integer ivar
  integer status,system

  real,allocatable :: time(:)
  real,allocatable :: z(:),zz(:)
  real,allocatable :: east_u(:,:),east_v(:,:),east_e(:,:)
  real,allocatable :: north_u(:,:),north_v(:,:),north_e(:,:)
  real,allocatable :: h(:,:),fsm(:,:),dum(:,:),dvm(:,:)
  real,allocatable :: dat2d(:,:),mean2d(:,:),sprd2d(:,:),pass2d(:,:),miss2d(:,:)
  real,allocatable :: dat3d(:,:,:),mean3d(:,:,:),sprd3d(:,:,:),pass3d(:,:,:),miss3d(:,:,:)

  character(100) dir,region
  character(10) var
  character(4) year
  character(2) month

  !Read ensemble run information
  open(1,file="ens_info.txt",status="old")
  read(1,'(a)') dir
  read(1,'(a)') region
  read(1,*) syr,smon,sday
  read(1,*) iyr,imon,nday
  read(1,*) nmem
  read(1,*) iswitch_budget
  close(1)
  status=system("rm -f ens_info.txt")

  write(year,'(i4.4)') iyr
  write(month,'(i2.2)') imon

  !Read grid infomation (im,jm,km)
  ! write(*,*) "Read grid information"
  call read_ngrid(dir,region,im,jm,km)
  allocate(time(nday))
  allocate(z(km),zz(km))
  allocate(east_u(im,jm),east_v(im,jm),east_e(im,jm))
  allocate(north_u(im,jm),north_v(im,jm),north_e(im,jm))
  allocate(h(im,jm),fsm(im,jm),dum(im,jm),dvm(im,jm))
  allocate(dat2d(im,jm),mean2d(im,jm),sprd2d(im,jm),pass2d(im,jm),miss2d(im,jm))
  allocate(dat3d(im,jm,km),mean3d(im,jm,km),sprd3d(im,jm,km),pass3d(im,jm,km),miss3d(im,jm,km))

  !Read model information(time,...,dvm)
  ! write(*,*) "Read information"
  call read_info(dir,region,nday,im,jm,km, &
       & time,z,zz,east_u,east_v,east_e,north_u,north_v,north_e, &
       & h,fsm,dum,dvm)
  
  !Create netcdf file
  ! write(*,*) "Create netcdf"
  call create_nc("mean",dir,region,syr,smon,sday,im,jm,km,iswitch_budget)
  call create_nc("sprd",dir,region,syr,smon,sday,im,jm,km,iswitch_budget)
  
  !Write model information(z,...,dvm)
  ! write(*,*) "Write information"
  call write_info("mean",dir,region,nday,im,jm,km, &
       & time,z,zz,east_u,east_v,east_e,north_u,north_v,north_e, &
       & h,fsm,dum,dvm)
  call write_info("sprd",dir,region,nday,im,jm,km, &
       & time,z,zz,east_u,east_v,east_e,north_u,north_v,north_e, &
       & h,fsm,dum,dvm)

  do ivar=1,48

     !Read variable name
     call read_var(ivar,var)
     write(*,*) "Calculate Mean/Sprd of "//trim(var)//" in mesp_ens.out"

     do iday=1,nday
                
        do imem=1,nmem
           ! write(*,*) "mean",imem,"/",nmem
           if(1 <= ivar .and. ivar <= 22)then
              ! write(*,*) "Read ensemble"
              call read_ens(var,dir,region,imem,iday,im,jm,1,dat2d)
              ! write(*,*) "Calculate Mean"
              call ensemble_mean(imem,nmem,im,jm,1,dat2d,mean2d,pass2d,miss2d)
           elseif(23 <= ivar .and. ivar <= 48)then
              ! write(*,*) "Read ensemble"
              call read_ens(var,dir,region,imem,iday,im,jm,km,dat3d)
              ! write(*,*) "Calculate Mean"
              call ensemble_mean(imem,nmem,im,jm,km,dat3d,mean3d,pass3d,miss3d)
           end if
        end do

        do imem=1,nmem
           ! write(*,*) "sprd",imem,"/",nmem
           if(1 <= ivar .and. ivar <= 22)then
              ! write(*,*) "Read ensemble"
              call read_ens(var,dir,region,imem,iday,im,jm,1,dat2d)
              ! write(*,*) "Calculate Sprd"
              call ensemble_sprd(imem,nmem,im,jm,1,dat2d,mean2d,sprd2d,pass2d,miss2d)
           elseif(23 <= ivar .and. ivar <= 48)then
              ! write(*,*) "Read ensemble"
              call read_ens(var,dir,region,imem,iday,im,jm,km,dat3d)
              ! write(*,*) "Calculate Sprd"
              call ensemble_sprd(imem,nmem,im,jm,km,dat3d,mean3d,sprd3d,pass3d,miss3d)
           end if
        end do

        if(1 <= ivar .and. ivar <= 22)then
           ! write(*,*) "Write Mean"
           call write_ens("mean",var,dir,region,iday,im,jm,1,mean2d)
           ! write(*,*) "Write Spread"
           call write_ens("sprd",var,dir,region,iday,im,jm,1,sprd2d)
        elseif(23 <= ivar .and. ivar <= 48)then
           ! write(*,*) "Write Mean"
           call write_ens("mean",var,dir,region,iday,im,jm,km,mean3d)
           ! write(*,*) "Write Spread"
           call write_ens("sprd",var,dir,region,iday,im,jm,km,sprd3d)
        end if
        
     end do
  end do

  deallocate(time)
  deallocate(z,zz)
  deallocate(east_u,east_v,east_e)
  deallocate(north_u,north_v,north_e)
  deallocate(h,fsm,dum,dvm)
  deallocate(dat2d,mean2d,sprd2d,pass2d,miss2d)
  deallocate(dat3d,mean3d,sprd3d,pass3d,miss3d)

  write(*,*) "=== End Mean & Sprd ==="

end program

!----------------------------------------------------------------------
! Read variable name |
!----------------------------------------------------------------------

subroutine read_var(ivar,var)

  implicit none

  integer ivar
  character(10) var

  !2D
  if(ivar == 1) var="el"
  if(ivar == 2) var="lhf"
  if(ivar == 3) var="shf"
  if(ivar == 4) var="lwr"
  if(ivar == 5) var="swr"
  if(ivar == 6) var="windu"
  if(ivar == 7) var="windv"
  if(ivar == 8) var="winds"
  if(ivar == 9) var="tauu"
  if(ivar == 10) var="tauv"
  if(ivar == 11) var="taus"
  if(ivar == 12) var="qa"
  if(ivar == 13) var="qs"
  if(ivar == 14) var="ta"
  if(ivar == 15) var="evap"
  if(ivar == 16) var="prep"
  if(ivar == 17) var="river"
  if(ivar == 18) var="eflux"
  if(ivar == 19) var="pflux"
  if(ivar == 20) var="rflux"
  if(ivar == 21) var="tsfc"
  if(ivar == 22) var="ssfc"

  !3D
  if(ivar == 23) var="u"
  if(ivar == 24) var="v"
  if(ivar == 25) var="w"
  if(ivar == 26) var="t"
  if(ivar == 27) var="s"
  if(ivar == 28) var="dtdt"
  if(ivar == 29) var="txadv"
  if(ivar == 30) var="tyadv"
  if(ivar == 31) var="tzadv"
  if(ivar == 32) var="tadv"
  if(ivar == 33) var="txdif"
  if(ivar == 34) var="tydif"
  if(ivar == 35) var="tzdif"
  if(ivar == 36) var="qz"
  if(ivar == 37) var="tnudge"
  if(ivar == 38) var="dsdt"
  if(ivar == 39) var="sxadv"
  if(ivar == 40) var="syadv"
  if(ivar == 41) var="szadv"
  if(ivar == 42) var="sadv"
  if(ivar == 43) var="sxdif"
  if(ivar == 44) var="sydif"
  if(ivar == 45) var="szdif"
  if(ivar == 46) var="snudge"
  if(ivar == 47) var="aam"
  if(ivar == 48) var="kh"

end subroutine read_var

!-----------------------------------------------------------------------
! Read Ensemble |
!-----------------------------------------------------------------------

subroutine read_ngrid(dir,region,im,jm,km)

  implicit none
  include "netcdf.inc"

  integer im,jm,km
  integer status,access
  integer ncid,varid,dimid

  real time

  character(100) dir,region
  character(200) filename

  filename=trim(dir)//"/001/out/"//trim(region)//".nc"
  
  status=access(trim(filename)," ")
  if(status /= 0)then
     write(*,*) "***Error: Not Found "//trim(filename)
     stop
  end if

  status=nf_open(trim(filename),nf_nowrite,ncid)

  status=nf_inq_dimid(ncid,"x",dimid)
  status=nf_inq_dimlen(ncid,dimid,im)

  status=nf_inq_dimid(ncid,"y",dimid)
  status=nf_inq_dimlen(ncid,dimid,jm)

  status=nf_inq_dimid(ncid,"z",dimid)
  status=nf_inq_dimlen(ncid,dimid,km)

  status=nf_close(ncid)

end subroutine read_ngrid

!-----------------------------------------

subroutine read_info(dir,region,nday,im,jm,km, &
     & time,z,zz,east_u,east_v,east_e,north_u,north_v,north_e, &
     & h,fsm,dum,dvm)

  implicit none
  include "netcdf.inc"

  integer nday
  integer im,jm,km
  integer status,access
  integer ncid,varid
  
  real time(nday)
  real z(km),zz(km)
  real east_u(im,jm),east_v(im,jm),east_e(im,jm)
  real north_u(im,jm),north_v(im,jm),north_e(im,jm)
  real h(im,jm),fsm(im,jm),dum(im,jm),dvm(im,jm)

  character(100) dir,region
  character(200) filename

  filename=trim(dir)//"/001/out/"//trim(region)//".nc"

  status=access(trim(filename)," ")
  if(status /= 0)then
     write(*,*) "***Error: Not Found "//trim(filename)
     stop
  end if

  status=nf_open(trim(filename),nf_nowrite,ncid)

  status=nf_inq_varid(ncid,"time",varid)
  status=nf_get_var_real(ncid,varid,time)

  status=nf_inq_varid(ncid,"z",varid)
  status=nf_get_var_real(ncid,varid,z)

  status=nf_inq_varid(ncid,"zz",varid)
  status=nf_get_var_real(ncid,varid,zz)

  status=nf_inq_varid(ncid,"east_u",varid)
  status=nf_get_var_real(ncid,varid,east_u)

  status=nf_inq_varid(ncid,"east_v",varid)
  status=nf_get_var_real(ncid,varid,east_v)

  status=nf_inq_varid(ncid,"east_e",varid)
  status=nf_get_var_real(ncid,varid,east_e)

  status=nf_inq_varid(ncid,"north_u",varid)
  status=nf_get_var_real(ncid,varid,north_u)

  status=nf_inq_varid(ncid,"north_v",varid)
  status=nf_get_var_real(ncid,varid,north_v)

  status=nf_inq_varid(ncid,"north_e",varid)
  status=nf_get_var_real(ncid,varid,north_e)

  status=nf_inq_varid(ncid,"h",varid)
  status=nf_get_var_real(ncid,varid,h)

  status=nf_inq_varid(ncid,"fsm",varid)
  status=nf_get_var_real(ncid,varid,fsm)

  status=nf_inq_varid(ncid,"dum",varid)
  status=nf_get_var_real(ncid,varid,dum)

  status=nf_inq_varid(ncid,"dvm",varid)
  status=nf_get_var_real(ncid,varid,dvm)

  status=nf_close(ncid)

end subroutine read_info

!------------------------------------------------------------

subroutine read_ens(var,dir,region,imem,iday,im,jm,km,dat)

  implicit none
  include "netcdf.inc"

  integer imem
  integer iday
  integer im,jm,km
  integer status,access
  integer ncid,varid

  real dat(im,jm,km)

  character(3) mem
  character(10) var
  character(100) dir,region
  character(200) filename

  write(mem,'(i3.3)') imem  

  filename=trim(dir)//"/"//mem//"/out/"//trim(region)//".nc"
  
  status=access(trim(filename)," ")
  if(status /= 0)then
     write(*,*) "***Error: Not Found "//trim(filename)
     stop
  end if
  
  status=nf_open(trim(filename),nf_nowrite,ncid)

  status=nf_inq_varid(ncid,trim(var),varid)
  if(status == nf_noerr)then
     if(km == 1)then
        status=nf_get_vara_real(ncid,varid,(/1,1,iday/),(/im,jm,1/),dat)
     else
        status=nf_get_vara_real(ncid,varid,(/1,1,1,iday/),(/im,jm,km,1/),dat)
     end if
  else
     dat(:,:,:)=0.
  end if

  status=nf_close(ncid)

end subroutine read_ens

!----------------------------------------------------------------------
!  Calculate Ensemble Mean |
!----------------------------------------------------------------------

subroutine ensemble_mean(imem,nmem,im,jm,km,dat,mean,pass,miss)

  implicit none

  integer imem,nmem
  integer i,j,k,im,jm,km
  
  real dat(im,jm,km)
  real mean(im,jm,km),pass(im,jm,km),miss(im,jm,km)

  if(imem == 1)then
     mean(:,:,:)=0.
     pass(:,:,:)=0.
     miss(:,:,:)=0.
  end if

  do k=1,km
     do j=1,jm
        do i=1,im
           if(dat(i,j,k) == 0.)then
              miss(i,j,k)=miss(i,j,k)+1.
           else
              mean(i,j,k)=mean(i,j,k)+dat(i,j,k)
              pass(i,j,k)=pass(i,j,k)+1.
           end if
        end do
     end do
  end do

  if(imem == nmem)then
     do k=1,km
        do j=1,jm
           do i=1,im
              if(pass(i,j,k) == 0.)then
                 mean(i,j,k)=0.
              else
                 mean(i,j,k)=mean(i,j,k)/pass(i,j,k)
              end if
           end do
        end do
     end do
  end if

end subroutine ensemble_mean

!-----------------------------

subroutine ensemble_sprd(imem,nmem,im,jm,km,dat,mean,sprd,pass,miss)

  implicit none

  integer imem,nmem
  integer i,j,k,im,jm,km
  
  real dat(im,jm,km)
  real mean(im,jm,km),sprd(im,jm,km),pass(im,jm,km),miss(im,jm,km)

  if(imem == 1)then
     sprd(:,:,:)=0.
     pass(:,:,:)=0.
     miss(:,:,:)=0.
  end if

  do k=1,km
     do j=1,jm
        do i=1,im
           if(dat(i,j,k) == 0.)then
              miss(i,j,k)=miss(i,j,k)+1.
           else
              sprd(i,j,k)=sprd(i,j,k)+(dat(i,j,k)-mean(i,j,k))**2.
              pass(i,j,k)=pass(i,j,k)+1.
           end if
        end do
     end do
  end do

  if(imem == nmem)then
     do k=1,km
        do j=1,jm
           do i=1,im
              if(pass(i,j,k) == 0. .or. pass(i,j,k) == 1.)then
                 sprd(i,j,k)=0.
              else
                 sprd(i,j,k)=sqrt(sprd(i,j,k)/(pass(i,j,k)-1.))
              end if
           end do
        end do
     end do
  end if

end subroutine ensemble_sprd

!----------------------------------------------------------------------
! Create Netcdf |
!----------------------------------------------------------------------

subroutine def_var_netcdf(ncid,name,nvdims,vdims,varid,long_name,units,coords,lcoords)

  implicit none
  include "netcdf.inc"

  integer ncid,nvdims,varid
  integer vdims(nvdims)
  integer status,length
  
  character(*) name,long_name,units,coords

  logical lcoords

  status=nf_def_var(ncid,trim(name),nf_float,nvdims,vdims,varid)

  length=len(trim(long_name))
  status=nf_put_att_text(ncid,varid,'long_name',length,trim(long_name))

  length=len(trim(units))
  status=nf_put_att_text(ncid,varid,'units',length,trim(units))

  if(lcoords)then
     length=len(trim(coords))
     status=nf_put_att_text(ncid,varid,'coordinates',length,trim(coords))
  end if

end subroutine def_var_netcdf

!--------------------------------

subroutine create_nc(mesp,dir,region,iyr,imon,iday,im,jm,km,iswitch_budget)

  implicit none
  include "netcdf.inc"

  integer iyr,imon,iday
  integer im,jm,km
  integer status,system,access
  integer ncid,dimid,varid
  integer vdims(4)
  integer x_dimid,y_dimid,z_dimid,t_dimid
  integer length
  integer iswitch_budget !0:No, 1: Budget term
  integer iswitch_ff !0: No, 1: freshwater flux
  integer iswitch_tnudge !0: No, 1: T nudging
  integer iswitch_snudge !0: No, 1: S nudging

  character(200) filename
  character(100) dir,region
  character(4) year,mesp
  character(2) month,day

  write(year,'(i4.4)') iyr
  write(month,'(i2.2)') imon
  write(day,'(i2.2)') iday

  !Check freshwater flux & Nudging
  iswitch_ff=0
  iswitch_tnudge=0
  iswitch_snudge=0

  filename=trim(dir)//"/001/out/"//trim(region)//".nc"
  status=nf_open(trim(filename),nf_nowrite,ncid)
  status=nf_inq_varid(ncid,"evap",varid)
  if(status == nf_noerr)iswitch_ff=1
  status=nf_inq_varid(ncid,"tnudge",varid)
  if(status == nf_noerr)iswitch_tnudge=1
  status=nf_inq_varid(ncid,"snudge",varid)
  if(status == nf_noerr)iswitch_snudge=1
  status=nf_close(ncid)

  !Write filename
  filename=trim(dir)//"/"//trim(mesp)//"/out/"//trim(region)//".nc"

  !NF_CREATE
  status=access(trim(filename)," ")
  if(status == 0)then
     return
  endif
  
  !  status=system("rm -f "//trim(filename))
  status=nf_create(trim(filename),ior(nf_noclobber,nf_64bit_offset),ncid)

  !GLOBAL ATTRIBUTE
  length=len(trim(region))
  status=nf_put_att_text(ncid,nf_global,'title',length,trim(region))
  length=len('output file')
  status=nf_put_att_text(ncid,nf_global,'description',length,'output file')

  !TIME
  status=nf_def_dim(ncid,'time',nf_unlimited,t_dimid)
  vdims(1)=t_dimid
  call def_var_netcdf(ncid,'time',1,vdims,varid, &
       & 'time','days since '//year//'-'//month//'-'//day//' 00:00:00 +00:00', &
       & ' ',.false.)

  !Z
  status=nf_def_dim(ncid,'z',km,z_dimid)
  vdims(1)=z_dimid
  call def_var_netcdf(ncid,'z',1,vdims,varid, &
       & 'sigma of cell face','sigma_level', &
       & ' ',.false.)
  length=len(trim('ocean_sigma_coordinate'))
  status=nf_put_att_text(ncid,varid,'standard_name', &
       & length,'ocean_sigma_coordinate')
  length=len(trim('sigma: z eta: el depth: h'))
  status=nf_put_att_text(ncid,varid,'formula_terms', &
       & length,'sigma: z eta: el depth: h')

  !ZZ
  call def_var_netcdf(ncid,'zz',1,vdims,varid, &
       & 'sigma of cell centre','sigma_level', &
       & ' ',.false.)
  length=len(trim('ocean_sigma_coordinate'))
  status=nf_put_att_text(ncid,varid,'standard_name', &
       & length,'ocean_sigma_coordinate')
  length=len(trim('sigma: zz eta: el depth: h'))
  status=nf_put_att_text(ncid,varid,'formula_terms', &
       & length,'sigma: zz eta: el depth: h')

  !Y
  status=nf_def_dim(ncid,'y',jm,y_dimid)

  !X
  status=nf_def_dim(ncid,'x',im,x_dimid)

  !Position
  vdims(1)=x_dimid
  vdims(2)=y_dimid

  call def_var_netcdf(ncid,'east_u',2,vdims,varid, &
       & 'easting of u-points','metre', &
       & 'east_u north_u',.true.)
  call def_var_netcdf(ncid,'east_v',2,vdims,varid, &
       & 'easting of v-points','metre', &
       & 'east_v north_v',.true.)
  call def_var_netcdf(ncid,'east_e',2,vdims,varid, &
       & 'easting of elevation points','metre', &
       & 'east_e north_e',.true.)
  call def_var_netcdf(ncid,'north_u',2,vdims,varid, &
       & 'northing of u-points','metre', &
       & 'east_u north_u',.true.)
  call def_var_netcdf(ncid,'north_v',2,vdims,varid, &
       & 'northing of v-points','metre', &
       & 'east_v north_v',.true.)
  call def_var_netcdf(ncid,'north_e',2,vdims,varid, &
       & 'northing of elevation points','metre', &
       & 'east_e north_e',.true.)
  call def_var_netcdf(ncid,'h',2,vdims,varid, &
       & 'undisturbed water depth','metre', &
       & 'east_e north_e',.true.)
  call def_var_netcdf(ncid,'fsm',2,vdims,varid, &
       & 'free surface mask','dimensionless', &
       & 'east_e north_e',.true.)
  call def_var_netcdf(ncid,'dum',2,vdims,varid, &
       & 'u-velocity mask','dimensionless', &
       & 'east_u north_u',.true.)
  call def_var_netcdf(ncid,'dvm',2,vdims,varid, &
       & 'v-velocity mask','dimensionless', &
       & 'east_v north_v',.true.)
  
  !2D variable
  vdims(1)=x_dimid
  vdims(2)=y_dimid
  vdims(3)=t_dimid
  
  call def_var_netcdf(ncid,'el',3,vdims,varid, &
       & 'surface elevation','metre', &
       & 'east_e north_e',.true.)
  call def_var_netcdf(ncid,'lhf',3,vdims,varid, &
       & 'latent heat flux','W/m^2', &
       & 'east_e north_e',.true.)
  call def_var_netcdf(ncid,'shf',3,vdims,varid, &
       & 'sensible heat flux','W/m^2', &
       & 'east_e north_e',.true.)
  call def_var_netcdf(ncid,'lwr',3,vdims,varid, &
       & 'longwave radiation','W/m^2', &
       & 'east_e north_e',.true.)
  call def_var_netcdf(ncid,'swr',3,vdims,varid, &
       & 'shortwave radiation','W/m^2', &
       & 'east_e north_e',.true.)
  call def_var_netcdf(ncid,'windu',3,vdims,varid, &
       & 'zonal wind','m/s', &
       & 'east_e north_e',.true.)
  call def_var_netcdf(ncid,'windv',3,vdims,varid, &
       & 'meridional wind','m/s', &
       & 'east_e north_e',.true.)
  call def_var_netcdf(ncid,'winds',3,vdims,varid, &
       & 'wind speed','m/s', &
       & 'east_e north_e',.true.)
  call def_var_netcdf(ncid,'tauu',3,vdims,varid, &
       & 'zonal wind stress','N/m^2', &
       & 'east_u north_u',.true.)
  call def_var_netcdf(ncid,'tauv',3,vdims,varid, &
       & 'meridional wind stress','N/m^2', &
       & 'east_v north_v',.true.)
  call def_var_netcdf(ncid,'taus',3,vdims,varid, &
       & 'magnitude of wind stress','N/m^2', &
       & 'east_v north_v',.true.)
  call def_var_netcdf(ncid,'qa',3,vdims,varid, &
       & 'air specific humidity','g/kg', &
       & 'east_e north_e',.true.)
  call def_var_netcdf(ncid,'qs',3,vdims,varid, &
       & 'surface saturated specific humidity','g/kg', &
       & 'east_e north_e',.true.)
  call def_var_netcdf(ncid,'ta',3,vdims,varid, &
       & 'air temperature','degree C', &
       & 'east_e north_e',.true.)
  
  if(iswitch_ff == 1)then
     call def_var_netcdf(ncid,'evap',3,vdims,varid, &
                   'evaporation','mm/day', &
                   'east_e north_e',.true.)
     call def_var_netcdf(ncid,'prep',3,vdims,varid, &
                   'precipitation','mm/day', &
                   'east_e north_e',.true.)
     call def_var_netcdf(ncid,'river',3,vdims,varid, &
                   'river discharge','mm/day', &
                   'east_e north_e',.true.)
     call def_var_netcdf(ncid,'eflux',3,vdims,varid, &
                   'evaporation flux','m/s', &
                   'east_e north_e',.true.)
     call def_var_netcdf(ncid,'pflux',3,vdims,varid, &
                   'precipitation flux','m/s', &
                   'east_e north_e',.true.)
     call def_var_netcdf(ncid,'rflux',3,vdims,varid, &
                   'river discharge flux','m/s', &
                   'east_e north_e',.true.)
  end if
  
  !3D variable
  vdims(1)=x_dimid
  vdims(2)=y_dimid
  vdims(3)=z_dimid
  vdims(4)=t_dimid
  
  call def_var_netcdf(ncid,'u',4,vdims,varid, &
       & 'x-velocity','metre/sec', &
       & 'east_u north_u zz',.true.)
  call def_var_netcdf(ncid,'v',4,vdims,varid, &
       & 'y-velocity','metre/sec', &
       & 'east_v north_v zz',.true.)
  call def_var_netcdf(ncid,'w',4,vdims,varid, &
       & 'z-velocity','metre/sec', &
       & 'east_e north_e z',.true.)
  call def_var_netcdf(ncid,'t',4,vdims,varid, &
       & 'potential temperature','degree C', &
       & 'east_e north_e zz',.true.)
  call def_var_netcdf(ncid,'s',4,vdims,varid, &
       & 'salinity','psu', &
       & 'east_e north_e zz',.true.)

  if(iswitch_budget == 1)then
     call def_var_netcdf(ncid,'dtdt',4,vdims,varid, &
          & 'temperature tendency','degree C/day', &
          & 'east_e north_e zz',.true.)
     call def_var_netcdf(ncid,'txadv',4,vdims,varid, &
          & 'x-temp advection','degree C/day', &
          & 'east_e north_e zz',.true.)
     call def_var_netcdf(ncid,'tyadv',4,vdims,varid, &
          & 'y-temp advection','degree C/day', &
          & 'east_e north_e zz',.true.)
     call def_var_netcdf(ncid,'tzadv',4,vdims,varid, &
          & 'z-temp advection','degree C/day', &
          & 'east_e north_e zz',.true.)
     call def_var_netcdf(ncid,'tadv',4,vdims,varid, &
          & 'temp advection','degree C/day', &
          & 'east_e north_e zz',.true.)
     call def_var_netcdf(ncid,'txdif',4,vdims,varid, &
          & 'x-temp diffusion','degree C/day', &
          & 'east_e north_e zz',.true.)
     call def_var_netcdf(ncid,'tydif',4,vdims,varid, &
          & 'y-temp diffusion','degree C/day', &
          & 'east_e north_e zz',.true.)
     call def_var_netcdf(ncid,'tzdif',4,vdims,varid, &
          & 'z-temp diffusion','degree C/day', &
          & 'east_e north_e zz',.true.)
     call def_var_netcdf(ncid,'qz',4,vdims,varid, &
          & 'shortwave penetration','degree C/day', &
          & 'east_e north_e zz',.true.)
     call def_var_netcdf(ncid,'dsdt',4,vdims,varid, &
          & 'salinity tendency','psu/day', &
          & 'east_e north_e zz',.true.)
     call def_var_netcdf(ncid,'sxadv',4,vdims,varid, &
          & 'x-sal advection','psu/day', &
          & 'east_e north_e zz',.true.)
     call def_var_netcdf(ncid,'syadv',4,vdims,varid, &
          & 'y-sal advection','psu/day', &
          & 'east_e north_e zz',.true.)
     call def_var_netcdf(ncid,'szadv',4,vdims,varid, &
          & 'z-sal advection','psu/day', &
          & 'east_e north_e zz',.true.)
     call def_var_netcdf(ncid,'sadv',4,vdims,varid, &
          & 'sal advection','psu/day', &
          & 'east_e north_e zz',.true.)
     call def_var_netcdf(ncid,'sxdif',4,vdims,varid, &
          & 'x-sal diffusion','psu/day', &
          & 'east_e north_e zz',.true.)
     call def_var_netcdf(ncid,'sydif',4,vdims,varid, &
          & 'y-sal diffusion','psu/day', &
          & 'east_e north_e zz',.true.)
     call def_var_netcdf(ncid,'szdif',4,vdims,varid, &
          & 'z-sal diffusion','psu/day', &
          & 'east_e north_e zz',.true.)
     call def_var_netcdf(ncid,'aam',4,vdims,varid, &
          & 'Am (Ah=0.2*Am)','meter^2/sec', &
          & 'east_e north_e zz',.true.)
     call def_var_netcdf(ncid,'kh',4,vdims,varid, &
          & 'vertical diffusivity','meter^2/sec', &
          & 'east_e north_e zz',.true.)           

     if(iswitch_tnudge == 1)then
        call def_var_netcdf(ncid,'tnudge',4,vdims,varid, &
             &           'temperature nudging','degree C/day', &
             &           'east_e north_e zz',.true.)
     end if

     if(iswitch_snudge == 1)then
        call def_var_netcdf(ncid,'snudge',4,vdims,varid, &
             &           'salinity nudging','psu/day', &
             &           'east_e north_e zz',.true.)
     end if

  end if

  !END Definition
  status=nf_enddef(ncid)

  status=nf_close(ncid)

end subroutine create_nc

!----------------------------------------------------------------------
! Write Netcdf |
!----------------------------------------------------------------------

subroutine write_info(mesp,dir,region,nday,im,jm,km, &
     & time,z,zz,east_u,east_v,east_e,north_u,north_v,north_e, &
     & h,fsm,dum,dvm)

  implicit none
  include "netcdf.inc"

  integer iday,nday
  integer im,jm,km
  integer status,access
  integer ncid,varid
  
  real time(nday)
  real z(km),zz(km)
  real east_u(im,jm),east_v(im,jm),east_e(im,jm)
  real north_u(im,jm),north_v(im,jm),north_e(im,jm)
  real h(im,jm),fsm(im,jm),dum(im,jm),dvm(im,jm)

  character(4) mesp
  character(100) dir,region
  character(200) filename

  filename=trim(dir)//"/"//trim(mesp)//"/out/"//trim(region)//".nc"

  status=access(trim(filename)," ")
  if(status /= 0)then
     write(*,*) "***Error: Not Found "//trim(filename)
     stop
  end if

  status=nf_open(trim(filename),nf_write,ncid)

  status=nf_inq_varid(ncid,"time",varid)
  status=nf_put_vara_real(ncid,varid,(/1/),(/nday/),time)

  status=nf_inq_varid(ncid,"z",varid)
  status=nf_put_var_real(ncid,varid,z)

  status=nf_inq_varid(ncid,"zz",varid)
  status=nf_put_var_real(ncid,varid,zz)

  status=nf_inq_varid(ncid,"east_u",varid)
  status=nf_put_var_real(ncid,varid,east_u)

  status=nf_inq_varid(ncid,"east_v",varid)
  status=nf_put_var_real(ncid,varid,east_v)

  status=nf_inq_varid(ncid,"east_e",varid)
  status=nf_put_var_real(ncid,varid,east_e)

  status=nf_inq_varid(ncid,"north_u",varid)
  status=nf_put_var_real(ncid,varid,north_u)

  status=nf_inq_varid(ncid,"north_v",varid)
  status=nf_put_var_real(ncid,varid,north_v)

  status=nf_inq_varid(ncid,"north_e",varid)
  status=nf_put_var_real(ncid,varid,north_e)

  status=nf_inq_varid(ncid,"h",varid)
  status=nf_put_var_real(ncid,varid,h)

  status=nf_inq_varid(ncid,"fsm",varid)
  status=nf_put_var_real(ncid,varid,fsm)

  status=nf_inq_varid(ncid,"dum",varid)
  status=nf_put_var_real(ncid,varid,dum)

  status=nf_inq_varid(ncid,"dvm",varid)
  status=nf_put_var_real(ncid,varid,dvm)

  status=nf_close(ncid)

end subroutine write_info

!-------------------------------------------

subroutine write_ens(mesp,var,dir,region,iday,im,jm,km,dat)

  implicit none
  include "netcdf.inc"

  integer iday
  integer im,jm,km
  integer status,access
  integer ncid,varid

  real dat(im,jm,km)

  character(4) mesp
  character(10) var
  character(100) dir,region
  character(200) filename

  filename=trim(dir)//"/"//trim(mesp)//"/out/"//trim(region)//".nc"
  
  status=access(trim(filename)," ")
  if(status /= 0)then
     write(*,*) "***Error: Not Found "//trim(filename)
     stop
  end if
  
  status=nf_open(trim(filename),nf_write,ncid)

  status=nf_inq_varid(ncid,trim(var),varid)
  if(status == nf_noerr)then
     if(km == 1)then
        status=nf_put_vara_real(ncid,varid,(/1,1,iday/),(/im,jm,1/),dat)
     else
        status=nf_put_vara_real(ncid,varid,(/1,1,1,iday/),(/im,jm,km,1/),dat)
     end if
  end if

  status=nf_close(ncid)

end subroutine write_ens
