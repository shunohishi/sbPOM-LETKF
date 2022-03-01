!_______________________________________________________________________
      subroutine write_diag_netcdf

! write diagnosys data

      implicit none
# include "mpif.h"
# include "netcdf.inc"
      include 'pom.h'
      character*120 netcdf_out_file,str_tmp
!      integer nprint
      integer nday,nhour
      integer time_dimid,x_dimid,y_dimid,z_dimid
      integer z_varid,zz_varid,dx_varid,dy_varid,east_c_varid,
     $        east_e_varid,east_u_varid,east_v_varid,north_c_varid,
     $        north_e_varid,north_u_varid,north_v_varid,rot_varid,
     $        h_varid,fsm_varid,dum_varid,dvm_varid
      integer time_varid,sfcf_varid,xadv_varid,yadv_varid,
     $        wadv_varid,xdif_varid,ydif_varid,vdif_varid,
     $        rad_varid,gvd_varid,temp_varid,
     $        windu_varid,windv_varid,u_varid,v_varid
      integer ncid,status
      integer vdims(4)
      integer length

! -- local
      real glb2d(im_global,jm_global),glb3d(im_global,jm_global,kb)

! debug
!      character*1 cha_my_task
      integer i,j

! diag
      include 'pom_diag.h'
      character*100 netcdf_diag_file
      netcdf_diag_file='hbgt'

! create netcdf file

      if(my_task.eq.master_task) then

!        nprint=(iint+time0*86400/dti)/iprint
!        write(netcdf_out_file,'(''out/'',a,''.'',i4.4,''.nc'')') 
!     $                                  trim(netcdf_diag_file),nprint
        nday=(iint*dti+time0*86400)/86400.
        nhour=((iint*dti+time0*86400)-nday*86400.)/3600.
        write(netcdf_out_file,
     $      '(''out/'',a,''.'',i4.4,''.'',i2.2,''.nc'')') 
     $       trim(netcdf_diag_file),nday,nhour
        write(*,'(/''writing file '',a)') trim(netcdf_out_file)
        status=nf_create(netcdf_out_file,nf_clobber,ncid)
        call handle_error_netcdf('nf_create: '//netcdf_out_file,
     $                          status,nf_noerr)

! define global attributes
        length=len(trim(title))
        status=nf_put_att_text(ncid,nf_global,'title',length,
     $                          trim(title))
        call handle_error_netcdf('nf_put_att_text',status,nf_noerr)
        str_tmp='output file'
        length=len(trim(str_tmp))
        status=nf_put_att_text(ncid,nf_global,'description',length,
     $                         trim(str_tmp))
        call handle_error_netcdf('nf_put_att_text',status,nf_noerr)

! define dimensions
        length=1
        status=nf_def_dim(ncid,'time',length,time_dimid)
        call handle_error_netcdf('nf_def_dim: time',status,nf_noerr)
        length=kb
        status=nf_def_dim(ncid,'z',length,z_dimid)
        call handle_error_netcdf('nf_def_dim: z',status,nf_noerr)
        length=jm_global
        status=nf_def_dim(ncid,'y',length,y_dimid)
        call handle_error_netcdf('nf_def_dim: y',status,nf_noerr)
        length=im_global
        status=nf_def_dim(ncid,'x',length,x_dimid)
        call handle_error_netcdf('nf_def_dim: x',status,nf_noerr)

! define variables and their attributes
        vdims(1)=time_dimid
        str_tmp='days since '//time_start
        call def_var_netcdf(ncid,'time',1,vdims,time_varid,
     $                     'time',str_tmp,
     $                     ' ',.false.)

        vdims(1)=z_dimid
        call def_var_netcdf(ncid,'z',1,vdims,z_varid,
     $                     'sigma of cell face','sigma_level',
     $                     ' ',.false.)
        length=22
        status=nf_put_att_text(ncid,z_varid,'standard_name',
     $                          length,'ocean_sigma_coordinate')
        length=26
        status=nf_put_att_text(ncid,z_varid,'formula_terms',
     $                          length,'sigma: z eta: elb depth: h')
        call handle_error_netcdf('nf_put_att_text: z',status,nf_noerr)

        call def_var_netcdf(ncid,'zz',1,vdims,zz_varid,
     $                     'sigma of cell centre','sigma_level',
     $                     ' ',.false.)
        length=22
        status=nf_put_att_text(ncid,zz_varid,'standard_name',
     $                         length,'ocean_sigma_coordinate')
        length=27
        status=nf_put_att_text(ncid,zz_varid,'formula_terms',
     $                          length,'sigma: zz eta: elb depth: h')
        call handle_error_netcdf('nf_put_att_text: zz',status,nf_noerr)

        vdims(1)=x_dimid
        vdims(2)=y_dimid
        call def_var_netcdf(ncid,'dx',2,vdims,dx_varid,
     $                     'grid increment in x','metre',
     $                     'east_e north_e',.true.)
        call def_var_netcdf(ncid,'dy',2,vdims,dy_varid,
     $                     'grid increment in y','metre',
     $                     'east_e north_e',.true.)
        call def_var_netcdf(ncid,'east_u',2,vdims,east_u_varid,
     $                     'easting of u-points','metre',
     $                     'east_u north_u',.true.)
        call def_var_netcdf(ncid,'east_v',2,vdims,east_v_varid,
     $                     'easting of v-points','metre',
     $                     'east_v north_v',.true.)
        call def_var_netcdf(ncid,'east_e',2,vdims,east_e_varid,
     $                     'easting of elevation points','metre',
     $                     'east_e north_e',.true.)
        call def_var_netcdf(ncid,'east_c',2,vdims,east_c_varid,
     $                     'easting of cell corners','metre',
     $                     'east_c north_c',.true.)
        call def_var_netcdf(ncid,'north_u',2,vdims,north_u_varid,
     $                     'northing of u-points','metre',
     $                     'east_u north_u',.true.)
        call def_var_netcdf(ncid,'north_v',2,vdims,north_v_varid,
     $                     'northing of v-points','metre',
     $                     'east_v north_v',.true.)
        call def_var_netcdf(ncid,'north_e',2,vdims,north_e_varid,
     $                     'northing of elevation points','metre',
     $                     'east_e north_e',.true.)
        call def_var_netcdf(ncid,'north_c',2,vdims,north_c_varid,
     $                     'northing of cell corners','metre',
     $                     'east_c north_c',.true.)
        call def_var_netcdf(ncid,'rot',2,vdims,rot_varid,
     $                    'Rotation angle of x-axis wrt. east','degree',
     $                     'east_e north_e',.true.)
        call def_var_netcdf(ncid,'h',2,vdims,h_varid,
     $                     'undisturbed water depth','metre',
     $                     'east_e north_e',.true.)
        call def_var_netcdf(ncid,'fsm',2,vdims,fsm_varid,
     $                     'free surface mask','dimensionless',
     $                     'east_e north_e',.true.)
        call def_var_netcdf(ncid,'dum',2,vdims,dum_varid,
     $                     'u-velocity mask','dimensionless',
     $                     'east_u north_u',.true.)
        call def_var_netcdf(ncid,'dvm',2,vdims,dvm_varid,
     $                     'v-velocity mask','dimensionless',
     $                     'east_v north_v',.true.)
     
        vdims(1)=x_dimid
        vdims(2)=y_dimid
        vdims(3)=time_dimid
        call def_var_netcdf(ncid,'sfcterm',3,vdims,sfcf_varid,
     $                     'surface heat flux','K/sec',
     $                     'east_e north_e',.true.)
        call def_var_netcdf(ncid,'windu',3,vdims,windu_varid,
     $                     'surface eastward wind','m/sec',
     $                     'east_u north_u',.true.)
        call def_var_netcdf(ncid,'windv',3,vdims,windv_varid,
     $                     'surface northward wind','m/sec',
     $                     'east_v north_v',.true.)

        vdims(1)=x_dimid
        vdims(2)=y_dimid
        vdims(3)=z_dimid
        vdims(4)=time_dimid
        call def_var_netcdf(ncid,'xadvterm',4,vdims,xadv_varid,
     $                     'x-advection','K/sec',
     $                     'east_e north_e zz',.true.)
        call def_var_netcdf(ncid,'yadvterm',4,vdims,yadv_varid,
     $                     'y-advection','K/sec',
     $                     'east_e north_e zz',.true.)
        call def_var_netcdf(ncid,'wadvterm',4,vdims,wadv_varid,
     $                     'w-advection','K/sec',
     $                     'east_e north_e zz',.true.)
        call def_var_netcdf(ncid,'xdifterm',4,vdims,xdif_varid,
     $                     'x-diffusion','K/sec',
     $                     'east_e north_e zz',.true.)
        call def_var_netcdf(ncid,'ydifterm',4,vdims,ydif_varid,
     $                     'y-diffusion','K/sec',
     $                     'east_e north_e zz',.true.)
        call def_var_netcdf(ncid,'vdifterm',4,vdims,vdif_varid,
     $                     'v-diffusion','K/sec',
     $                     'east_e north_e zz',.true.)
        call def_var_netcdf(ncid,'radterm',4,vdims,rad_varid,
     $                     'radiation','K/sec',
     $                     'east_e north_e zz',.true.)
        call def_var_netcdf(ncid,'gvdterm',4,vdims,gvd_varid,
     $                     'v-diffusion+heat','K/sec',
     $                     'east_e north_e zz',.true.)
        call def_var_netcdf(ncid,'t',4,vdims,temp_varid,
     $                     't','K',
     $                     'east_e north_e zz',.true.)
        call def_var_netcdf(ncid,'u',4,vdims,u_varid,
     $                     'u','m/sec',
     $                     'east_u north_u zz',.true.)
        call def_var_netcdf(ncid,'v',4,vdims,v_varid,
     $                     'v','m/sec',
     $                     'east_v north_v zz',.true.)

! end definitions
        status=nf_enddef(ncid)
        call handle_error_netcdf('nf_enddef: diag',status,nf_noerr)

      end if

! write data

      if(my_task.eq.master_task) then
        status=nf_put_var_real(ncid,time_varid,time)
        call handle_error_netcdf('nf_put_var_real:time',status,
     $                           nf_noerr)
        status=nf_put_var_real(ncid,z_varid,z)
        call handle_error_netcdf('nf_put_var_real: z',status,nf_noerr)
        status=nf_put_var_real(ncid,zz_varid,zz)
        call handle_error_netcdf('nf_put_var_real: zz',status,nf_noerr)
      end if

      call merge2d(glb2d,dx)
      if(my_task.eq.master_task) then
        status=nf_put_var_real(ncid,dx_varid,glb2d)
        call handle_error_netcdf('nf_put_var_real: dx',status,nf_noerr)
      end if

      call merge2d(glb2d,dy)
      if(my_task.eq.master_task) then
        status=nf_put_var_real(ncid,dy_varid,glb2d)
        call handle_error_netcdf('nf_put_var_real: dy',status,nf_noerr)
      end if

      call merge2d(glb2d,east_u)
      if(my_task.eq.master_task) then
        status=nf_put_var_real(ncid,east_u_varid,glb2d)
        call handle_error_netcdf('nf_put_var_real: east_u',status,
     $                             nf_noerr)
      end if

      call merge2d(glb2d,east_v)
      if(my_task.eq.master_task) then
        status=nf_put_var_real(ncid,east_v_varid,glb2d)
        call handle_error_netcdf('nf_put_var_real: east_v',status,
     $                             nf_noerr)
      end if

      call merge2d(glb2d,east_e)
      if(my_task.eq.master_task) then
        status=nf_put_var_real(ncid,east_e_varid,glb2d)
        call handle_error_netcdf('nf_put_var_real: east_e',status,
     $                             nf_noerr)
      end if

      call merge2d(glb2d,east_c)
      if(my_task.eq.master_task) then
        status=nf_put_var_real(ncid,east_c_varid,glb2d)
        call handle_error_netcdf('nf_put_var_real: east_c',status,
     $                              nf_noerr)
      end if

      call merge2d(glb2d,north_u)
      if(my_task.eq.master_task) then
        status=nf_put_var_real(ncid,north_u_varid,glb2d)
        call handle_error_netcdf('nf_put_var_real: north_u',status,
     $                             nf_noerr)
      end if

      call merge2d(glb2d,north_v)
      if(my_task.eq.master_task) then
        status=nf_put_var_real(ncid,north_v_varid,glb2d)
        call handle_error_netcdf('nf_put_var_real: north_v',status,
     $                             nf_noerr)
      end if

      call merge2d(glb2d,north_e)
      if(my_task.eq.master_task) then
        status=nf_put_var_real(ncid,north_e_varid,glb2d)
        call handle_error_netcdf('nf_put_var_real: north_e',status,
     $                           nf_noerr)
      end if

      call merge2d(glb2d,north_c)
      if(my_task.eq.master_task) then
        status=nf_put_var_real(ncid,north_c_varid,glb2d)
        call handle_error_netcdf('nf_put_var_real: north_c',status,
     $                           nf_noerr)
      end if

      call merge2d(glb2d,rot)
      if(my_task.eq.master_task) then
        status=nf_put_var_real(ncid,rot_varid,glb2d)
        call handle_error_netcdf('nf_put_var_real: rot',status,nf_noerr)
      end if

      call merge2d(glb2d,h)
      if(my_task.eq.master_task) then
        status=nf_put_var_real(ncid,h_varid,glb2d)
        call handle_error_netcdf('nf_put_var_real: h',status,nf_noerr)
      end if

      call merge2d(glb2d,fsm)
      if(my_task.eq.master_task) then
        status=nf_put_var_real(ncid,fsm_varid,glb2d)
        call handle_error_netcdf('nf_put_var_real: fsm',status,nf_noerr)
      end if

      call merge2d(glb2d,dum)
      if(my_task.eq.master_task) then
        status=nf_put_var_real(ncid,dum_varid,glb2d)
        call handle_error_netcdf('nf_put_var_real: dum',status,nf_noerr)
      end if

      call merge2d(glb2d,dvm)
      if(my_task.eq.master_task) then
        status=nf_put_var_real(ncid,dvm_varid,glb2d)
        call handle_error_netcdf('nf_put_var_real: dvm',status,nf_noerr)
      end if

      call merge2d(glb2d,sfcmean)
      if(my_task.eq.master_task) then
        status=nf_put_var_real(ncid,sfcf_varid,glb2d)
        call handle_error_netcdf('nf_put_var_real: sfcmean',status,
     $                           nf_noerr)
      end if

      call merge2d(glb2d,windu)
      if(my_task.eq.master_task) then
        status=nf_put_var_real(ncid,windu_varid,glb2d)
        call handle_error_netcdf('nf_put_var_real: windu',status,
     $                           nf_noerr)
      end if

      call merge2d(glb2d,windv)
      if(my_task.eq.master_task) then
        status=nf_put_var_real(ncid,windv_varid,glb2d)
        call handle_error_netcdf('nf_put_var_real: windv',status,
     $                           nf_noerr)
      end if

      call merge3d(glb3d,xadvmean)
      if(my_task.eq.master_task) then
        status=nf_put_var_real(ncid,xadv_varid,glb3d)
        call handle_error_netcdf
     $   ('nf_put_var_real: xadv',status,nf_noerr)
      end if

      call merge3d(glb3d,yadvmean)
      if(my_task.eq.master_task) then
        status=nf_put_var_real(ncid,yadv_varid,glb3d)
        call handle_error_netcdf
     $    ('nf_put_var_real: yadv',status,nf_noerr)
      end if

      call merge3d(glb3d,wadvmean)
      if(my_task.eq.master_task) then
        status=nf_put_var_real(ncid,wadv_varid,glb3d)
        call handle_error_netcdf
     $    ('nf_put_var_real: wadv',status,nf_noerr)
      end if

      call merge3d(glb3d,xdifmean)
      if(my_task.eq.master_task) then
        status=nf_put_var_real(ncid,xdif_varid,glb3d)
        call handle_error_netcdf
     $    ('nf_put_var_real: xdif',status,nf_noerr)
      end if

      call merge3d(glb3d,ydifmean)
      if(my_task.eq.master_task) then
        status=nf_put_var_real(ncid,ydif_varid,glb3d)
        call handle_error_netcdf
     $    ('nf_put_var_real: ydif',status,nf_noerr)
      end if

      call merge3d(glb3d,vdifmean)
      if(my_task.eq.master_task) then
        status=nf_put_var_real(ncid,vdif_varid,glb3d)
        call handle_error_netcdf
     $    ('nf_put_var_real: vdif',status,nf_noerr)
      end if

      call merge3d(glb3d,radmean)
      if(my_task.eq.master_task) then
        status=nf_put_var_real(ncid,rad_varid,glb3d)
        call handle_error_netcdf
     $    ('nf_put_var_real: rad',status,nf_noerr)
      end if

      call merge3d(glb3d,gvdmean)
      if(my_task.eq.master_task) then
        status=nf_put_var_real(ncid,gvd_varid,glb3d)
        call handle_error_netcdf
     $    ('nf_put_var_real: gvd',status,nf_noerr)
      end if

      call merge3d(glb3d,uf)
      if(my_task.eq.master_task) then
        status=nf_put_var_real(ncid,temp_varid,glb3d)
        call handle_error_netcdf
     $    ('nf_put_var_real: t',status,nf_noerr)
      end if

      call merge3d(glb3d,u)
      if(my_task.eq.master_task) then
        status=nf_put_var_real(ncid,u_varid,glb3d)
        call handle_error_netcdf
     $    ('nf_put_var_real: u',status,nf_noerr)
      end if

      call merge3d(glb3d,v)
      if(my_task.eq.master_task) then
        status=nf_put_var_real(ncid,v_varid,glb3d)
        call handle_error_netcdf
     $    ('nf_put_var_real: v',status,nf_noerr)
      end if

! close file
      if(my_task.eq.master_task) then
        status=nf_close(ncid)
        call handle_error_netcdf('nf_close: diag',status,nf_noerr)
      end if

! debug 
!      sfcmean=my_task
!      call merge2d(glb2d,sfcmean)
!      if(my_task.eq.master_task) then
!      open(88,file='output',form='formatted')
!        do j=1,jm_global
!         write(88,*) j,(glb2d(i,j),i=1,im_global)
!        end do
!        close(88)
!      end if
!      call finalize_mpi
!      stop

      return
      end
