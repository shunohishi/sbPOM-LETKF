! io_netcdf.F

! NetCDF I/O + MPI

!_______________________________________________________________________
      subroutine def_var_netcdf(ncid,name,nvdims,vdims,varid,
     $                          long_name,units,coords,lcoords)
# include "mpif.h"
# include "netcdf.inc"
      integer vdims(nvdims)
      integer ncid,nvdims,varid
      logical lcoords
      character*(*) name,long_name,units,coords
      integer status
      integer length
! define variable
      status=nf_def_var(ncid,name,nf_float,nvdims,vdims,varid)
      call handle_error_netcdf('nf_put_att_text',status,nf_noerr)
! define attributes
      length=len(trim(long_name))
      status=nf_put_att_text(ncid,varid,'long_name',
     $                          length,trim(long_name))
      call handle_error_netcdf('nf_put_att_text',status,nf_noerr)

      length=len(trim(units))
      status=nf_put_att_text(ncid,varid,'units',
     $                          length,trim(units))
      call handle_error_netcdf('nf_put_att_text',status,nf_noerr)
! add coordinates attribute, if necessary
      if(lcoords) then
        length=len(trim(coords))
        status=nf_put_att_text(ncid,varid,'coordinates',
     $                            length,trim(coords))
        call handle_error_netcdf('nf_put_att_text',status,nf_noerr)
      endif
      return
      end

!_______________________________________________________________________
      subroutine handle_error_netcdf(routine,status,nf_noerr)
      implicit none
      include 'pom.h'
      integer status,nf_noerr
      character*(*) routine
      if(status.ne.nf_noerr) then
        error_status=1
        if(my_task.eq.master_task) write(*,'(/a,a)')
     $                                  'Error: NetCDF routine ',routine
      end if
      return
      end

!_______________________________________________________________________
      subroutine write_output_netcdf
!     write output data
      implicit none
# include "mpif.h"
# include "netcdf.inc"
      include 'pom.h'
      character*120 netcdf_out_file,str_tmp
      integer nprint
      integer time_dimid,x_dimid,y_dimid,z_dimid
      integer z_varid,zz_varid,dx_varid,dy_varid,east_c_varid,
     $     east_e_varid,east_u_varid,east_v_varid,north_c_varid,
     $     north_e_varid,north_u_varid,north_v_varid,rot_varid,
     $     h_varid,fsm_varid,dum_varid,dvm_varid
      integer time_varid,uab_varid,vab_varid,el_varid,
     $     lhf_varid,shf_varid,lwr_varid,swr_varid,
     $     windu_varid,windv_varid,winds_varid,
     $     tauu_varid,tauv_varid,taus_varid,
     $     qa_varid,qs_varid,ta_varid,
     $     evap_varid,prep_varid,river_varid,
     $     eflux_varid,pflux_varid,rflux_varid,
     $     u_varid,v_varid,w_varid,t_varid,s_varid,rho_varid
!     Convervation
      integer dtdt_varid,dsdt_varid,
     $     txadv_varid,tyadv_varid,tzadv_varid,tadv_varid,
     $     txdif_varid,tydif_varid,tzdif_varid,
     $     tsfc_varid,qz_varid,
     $     sxadv_varid,syadv_varid,szadv_varid,sadv_varid,
     $     sxdif_varid,sydif_varid,szdif_varid,
     $     ssfc_varid,
     $     aam_varid,kh_varid
!     IAU
      integer tiau_varid,siau_varid
!     Round Off
      integer troff_varid,sroff_varid
!     Nudging
      integer tnudge_varid,snudge_varid

      integer ncid,status
      integer vdims(4)
      integer length

      real glb2d(im_global,jm_global),glb3d(im_global,jm_global,kb)

      integer i,j,k

!     create netcdf file
      if(my_task.eq.master_task)then

         nprint=nint(iint/real(iprint))
!     nprint=1
!     nprint=nint((iint+time0*86400/dti)/iprint)
         write(netcdf_out_file,'(''out/'',a,''.nc'')') trim(netcdf_file)
         write(*,'(/''writing file '',a)') trim(netcdf_out_file)

      end if
         
      if(my_task .eq. master_task .and. nprint .eq. 1)then
            
         status=nf_create(netcdf_out_file,
     $        ior(nf_noclobber,nf_64bit_offset),ncid)
         call handle_error_netcdf('nf_create: '//netcdf_out_file,
     $        status,nf_noerr)
         
!     define global attributes
         length=len(trim(title))
         status=nf_put_att_text(ncid,nf_global,'title',length,
     $        trim(title))
         call handle_error_netcdf('nf_put_att_text',status,nf_noerr)
         str_tmp='output file'
         length=len(trim(str_tmp))
         status=nf_put_att_text(ncid,nf_global,'description',length,
     $        trim(str_tmp))
         call handle_error_netcdf('nf_put_att_text',status,nf_noerr)
         
!     define dimensions
         status=nf_def_dim(ncid,'time',nf_unlimited,time_dimid)
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
         
!     define variables and their attributes
         vdims(1)=time_dimid
         str_tmp='days since '//time_start
         call def_var_netcdf(ncid,'time',1,vdims,time_varid,
     $        'time',str_tmp,
     $        ' ',.false.)
         
         vdims(1)=z_dimid
         call def_var_netcdf(ncid,'z',1,vdims,z_varid,
     $        'sigma of cell face','sigma_level',
     $        ' ',.false.)
         length=22
         status=nf_put_att_text(ncid,z_varid,'standard_name',
     $        length,'ocean_sigma_coordinate')
         call handle_error_netcdf
     $        ('nf_put_att_text: z',status,nf_noerr)
         length=26
         status=nf_put_att_text(ncid,z_varid,'formula_terms',
     $        length,'sigma: z eta: el depth: h')
         call handle_error_netcdf
     $        ('nf_put_att_text: z',status,nf_noerr)
         
         call def_var_netcdf(ncid,'zz',1,vdims,zz_varid,
     $        'sigma of cell centre','sigma_level',
     $        ' ',.false.)
         length=22
         status=nf_put_att_text(ncid,zz_varid,'standard_name',
     $        length,'ocean_sigma_coordinate')
         call handle_error_netcdf
     $        ('nf_put_att_text: zz',status,nf_noerr)
         length=27
         status=nf_put_att_text(ncid,zz_varid,'formula_terms',
     $        length,'sigma: zz eta: el depth: h')
         call handle_error_netcdf
     $        ('nf_put_att_text: zz',status,nf_noerr)
         
         vdims(1)=x_dimid
         vdims(2)=y_dimid
         call def_var_netcdf(ncid,'east_u',2,vdims,east_u_varid,
     $        'easting of u-points','metre',
     $        'east_u north_u',.true.)
         call def_var_netcdf(ncid,'east_v',2,vdims,east_v_varid,
     $        'easting of v-points','metre',
     $        'east_v north_v',.true.)
         call def_var_netcdf(ncid,'east_e',2,vdims,east_e_varid,
     $        'easting of elevation points','metre',
     $        'east_e north_e',.true.)
         call def_var_netcdf(ncid,'north_u',2,vdims,north_u_varid,
     $        'northing of u-points','metre',
     $        'east_u north_u',.true.)
         call def_var_netcdf(ncid,'north_v',2,vdims,north_v_varid,
     $        'northing of v-points','metre',
     $        'east_v north_v',.true.)
         call def_var_netcdf(ncid,'north_e',2,vdims,north_e_varid,
     $        'northing of elevation points','metre',
     $        'east_e north_e',.true.)
         call def_var_netcdf(ncid,'h',2,vdims,h_varid,
     $        'undisturbed water depth','metre',
     $        'east_e north_e',.true.)
         call def_var_netcdf(ncid,'fsm',2,vdims,fsm_varid,
     $        'free surface mask','dimensionless',
     $        'east_e north_e',.true.)
         call def_var_netcdf(ncid,'dum',2,vdims,dum_varid,
     $        'u-velocity mask','dimensionless',
     $        'east_u north_u',.true.)
         call def_var_netcdf(ncid,'dvm',2,vdims,dvm_varid,
     $        'v-velocity mask','dimensionless',
     $        'east_v north_v',.true.)
         
         vdims(1)=x_dimid
         vdims(2)=y_dimid
         vdims(3)=time_dimid
         call def_var_netcdf(ncid,'el',3,vdims,el_varid,
     $        'surface elevation','metre',
     $        'east_e north_e',.true.)
!     Heat fluxes (2018.07 S.Ohishi)
         call def_var_netcdf(ncid,'lhf',3,vdims,lhf_varid,
     $        'latent heat flux','W/m^2',
     $        'east_e north_e',.true.)
         call def_var_netcdf(ncid,'shf',3,vdims,shf_varid,
     $        'sensible heat flux','W/m^2',
     $        'east_e north_e',.true.)
         call def_var_netcdf(ncid,'lwr',3,vdims,lwr_varid,
     $        'longwave radiation','W/m^2',
     $        'east_e north_e',.true.)
         call def_var_netcdf(ncid,'swr',3,vdims,swr_varid,
     $        'shortwave radiation','W/m^2',
     $        'east_e north_e',.true.)
         call def_var_netcdf(ncid,'windu',3,vdims,windu_varid,
     $        'zonal wind','m/s',
     $        'east_e north_e',.true.)
         call def_var_netcdf(ncid,'windv',3,vdims,windv_varid,
     $        'meridional wind','m/s',
     $        'east_e north_e',.true.)
         call def_var_netcdf(ncid,'winds',3,vdims,winds_varid,
     $        'wind speed','m/s',
     $        'east_e north_e',.true.)
         call def_var_netcdf(ncid,'tauu',3,vdims,tauu_varid,
     $        'zonal wind stress','N/m^2',
     $        'east_u north_u',.true.)
         call def_var_netcdf(ncid,'tauv',3,vdims,tauv_varid,
     $        'meridional wind stress','N/m^2',
     $        'east_v north_v',.true.)
         call def_var_netcdf(ncid,'taus',3,vdims,taus_varid,
     $        'magnitude of wind stress','N/m^2',
     $        'east_v north_v',.true.)
         call def_var_netcdf(ncid,'qa',3,vdims,qa_varid,
     $        'air specific humidity','g/kg',
     $        'east_e north_e',.true.)
         call def_var_netcdf(ncid,'qs',3,vdims,qs_varid,
     $        'surface saturated specific humidity','g/kg',
     $        'east_e north_e',.true.)
         call def_var_netcdf(ncid,'ta',3,vdims,ta_varid,
     $        'air temperature','degree C',
     $        'east_e north_e',.true.)
         
         if(lssf .eq. 1)then
            call def_var_netcdf(ncid,'evap',3,vdims,evap_varid,
     $           'evaporation','mm/day',
     $           'east_e north_e',.true.)
            call def_var_netcdf(ncid,'prep',3,vdims,prep_varid,
     $           'precipitation','mm/day',
     $           'east_e north_e',.true.)
            call def_var_netcdf(ncid,'river',3,vdims,river_varid,
     $           'river discharge','mm/day',
     $           'east_e north_e',.true.)
            call def_var_netcdf(ncid,'eflux',3,vdims,eflux_varid,
     $           'evaporation flux','m/s',
     $           'east_e north_e',.true.)
            call def_var_netcdf(ncid,'pflux',3,vdims,pflux_varid,
     $           'precipitation flux','m/s',
     $           'east_e north_e',.true.)
            call def_var_netcdf(ncid,'rflux',3,vdims,rflux_varid,
     $           'river discharge flux','m/s',
     $           'east_e north_e',.true.)
         end if
         
         if(budget .eq. 1)then
            call def_var_netcdf(ncid,'tsfc',3,vdims,tsfc_varid,
     $           'surface heatflux forcing','degree C/day',
     $           'east_e north_e',.true.)

            call def_var_netcdf(ncid,'ssfc',3,vdims,ssfc_varid,
     $           'surface freshwater forcing','psu/day',
     $           'east_e north_e',.true.)
         end if
         
         vdims(1)=x_dimid
         vdims(2)=y_dimid
         vdims(3)=z_dimid
         vdims(4)=time_dimid
         call def_var_netcdf(ncid,'u',4,vdims,u_varid,
     $        'x-velocity','metre/sec',
     $        'east_u north_u zz',.true.)
         call def_var_netcdf(ncid,'v',4,vdims,v_varid,
     $        'y-velocity','metre/sec',
     $        'east_v north_v zz',.true.)
         call def_var_netcdf(ncid,'w',4,vdims,w_varid,
     $        'z-velocity','metre/sec',
     $        'east_e north_e z',.true.)
         call def_var_netcdf(ncid,'t',4,vdims,t_varid,
     $        'potential temperature','degree C',
     $        'east_e north_e zz',.true.)
         call def_var_netcdf(ncid,'s',4,vdims,s_varid,
     $        'salinity','psu',
     $        'east_e north_e zz',.true.)
         if(budget .eq. 1)then
            call def_var_netcdf(ncid,'dtdt',4,vdims,dtdt_varid,
     $           'temperature tendency','degree C/day',
     $           'east_e north_e zz',.true.)
            call def_var_netcdf(ncid,'txadv',4,vdims,txadv_varid,
     $           'x-temp advection','degree C/day',
     $           'east_e north_e zz',.true.)
            call def_var_netcdf(ncid,'tyadv',4,vdims,tyadv_varid,
     $           'y-temp advection','degree C/day',
     $           'east_e north_e zz',.true.)
            call def_var_netcdf(ncid,'tzadv',4,vdims,tzadv_varid,
     $           'z-temp advection','degree C/day',
     $           'east_e north_e zz',.true.)
            call def_var_netcdf(ncid,'tadv',4,vdims,tadv_varid,
     $           'temp advection','degree C/day',
     $           'east_e north_e zz',.true.)
            call def_var_netcdf(ncid,'txdif',4,vdims,txdif_varid,
     $           'x-temp diffusion','degree C/day',
     $           'east_e north_e zz',.true.)
            call def_var_netcdf(ncid,'tydif',4,vdims,tydif_varid,
     $           'y-temp diffusion','degree C/day',
     $           'east_e north_e zz',.true.)
            call def_var_netcdf(ncid,'tzdif',4,vdims,tzdif_varid,
     $           'z-temp diffusion','degree C/day',
     $           'east_e north_e zz',.true.)
            call def_var_netcdf(ncid,'qz',4,vdims,qz_varid,
     $           'shortwave penetration','degree C/day',
     $           'east_e north_e zz',.true.)
            call def_var_netcdf(ncid,'dsdt',4,vdims,dsdt_varid,
     $           'salinity tendency','psu/day',
     $           'east_e north_e zz',.true.)
            call def_var_netcdf(ncid,'sxadv',4,vdims,sxadv_varid,
     $           'x-sal advection','psu/day',
     $           'east_e north_e zz',.true.)
            call def_var_netcdf(ncid,'syadv',4,vdims,syadv_varid,
     $           'y-sal advection','psu/day',
     $           'east_e north_e zz',.true.)
            call def_var_netcdf(ncid,'szadv',4,vdims,szadv_varid,
     $           'z-sal advection','psu/day',
     $           'east_e north_e zz',.true.)
            call def_var_netcdf(ncid,'sadv',4,vdims,sadv_varid,
     $           'sal advection','psu/day',
     $           'east_e north_e zz',.true.)
            call def_var_netcdf(ncid,'sxdif',4,vdims,sxdif_varid,
     $           'x-sal diffusion','psu/day',
     $           'east_e north_e zz',.true.)
            call def_var_netcdf(ncid,'sydif',4,vdims,sydif_varid,
     $           'y-sal diffusion','psu/day',
     $           'east_e north_e zz',.true.)
            call def_var_netcdf(ncid,'szdif',4,vdims,szdif_varid,
     $           'z-sal diffusion','psu/day',
     $           'east_e north_e zz',.true.)
            call def_var_netcdf(ncid,'aam',4,vdims,aam_varid,
     $           'Am (Ah=0.2*Am)','meter^2/sec',
     $           'east_e north_e zz',.true.)
            call def_var_netcdf(ncid,'kh',4,vdims,kh_varid,
     $           'vertical diffusivity','meter^2/sec',
     $           'east_e north_e zz',.true.)           
         end if

         if(budget == 1 .and. assim == 1.)then
            call def_var_netcdf(ncid,'tiau',4,vdims,tiau_varid,
     $           'temperature for IAU method','deg C',
     $           'east_e north_e zz',.true.)           
            call def_var_netcdf(ncid,'siau',4,vdims,siau_varid,
     $           'salinity for IAU method','psu',
     $           'east_e north_e zz',.true.)
         end if

         if(budget == 1 .and. lroff)then
            call def_var_netcdf(ncid,'troff',4,vdims,troff_varid,
     $           'temperature round off','degree C/day',
     $           'east_e north_e zz',.true.)
            call def_var_netcdf(ncid,'sroff',4,vdims,sroff_varid,
     $           'salinity round off','psu/day',
     $           'east_e north_e zz',.true.)
         end if

         if(budget == 1)then
            call def_var_netcdf(ncid,'tnudge',4,vdims,tnudge_varid,
     $           'temperature nudging','degree C/day',
     $           'east_e north_e zz',.true.)
         end if

         if(budget == 1)then
            call def_var_netcdf(ncid,'snudge',4,vdims,snudge_varid,
     $           'salinity nudging','psu/day',
     $           'east_e north_e zz',.true.)
         end if
         
!     end definitions
         status=nf_enddef(ncid)
         call handle_error_netcdf
     $        ('nf_enddef: output',status,nf_noerr)

      end if

!     write data ------------------------------------------

      if(my_task .eq. master_task .and. nprint .ne. 1)then

         status=nf_open(netcdf_out_file,nf_write,ncid)
         call handle_error_netcdf('nf_open',status,nf_noerr)

      end if

      if(my_task .eq. master_task)then
         status=nf_inq_varid(ncid,"time",time_varid)
         call handle_error_netcdf('nf_inq_varid: t'
     $        ,status,nf_noerr)

         status=nf_put_vara_real(ncid,time_varid,
     $        (/nprint/),(/1/),time-1.)
         call handle_error_netcdf('nf_put_var_real:time',status,
     $        nf_noerr)
      end if
      
      if(my_task.eq.master_task .and. iint/iprint .eq. 1) then
         status=nf_inq_varid(ncid,"z",z_varid)
         call handle_error_netcdf('nf_inq_varid: z'
     $        ,status,nf_noerr)
         status=nf_put_var_real(ncid,z_varid,z)
         call handle_error_netcdf('nf_put_var_real: z',status,nf_noerr)
         status=nf_inq_varid(ncid,"zz",zz_varid)
      call handle_error_netcdf('nf_inq_varid: zz'
     $                          ,status,nf_noerr)
         status=nf_put_var_real(ncid,zz_varid,zz)
         call handle_error_netcdf('nf_put_var_real: zz',status,nf_noerr)
      end if
      
      call merge2d(glb2d,east_u)
      if(my_task.eq.master_task .and. iint/iprint .eq. 1) then
         status=nf_inq_varid(ncid,"east_u",east_u_varid)
         call handle_error_netcdf('nf_inq_varid: east_u'
     $        ,status,nf_noerr)
         status=nf_put_var_real(ncid,east_u_varid,glb2d)
         call handle_error_netcdf('nf_put_var_real: east_u',status,
     $        nf_noerr)
      end if
      
      call merge2d(glb2d,east_v)
      if(my_task.eq.master_task .and. iint/iprint .eq. 1) then
         status=nf_inq_varid(ncid,"east_v",east_v_varid)
         call handle_error_netcdf('nf_inq_varid: east_v'
     $        ,status,nf_noerr)
         status=nf_put_var_real(ncid,east_v_varid,glb2d)
         call handle_error_netcdf('nf_put_var_real: east_v',status,
     $        nf_noerr)
      end if
      
      call merge2d(glb2d,east_e)
      if(my_task.eq.master_task .and. iint/iprint .eq. 1) then
         status=nf_inq_varid(ncid,"east_e",east_e_varid)
         call handle_error_netcdf('nf_inq_varid: east_e'
     $        ,status,nf_noerr)
         status=nf_put_var_real(ncid,east_e_varid,glb2d)
         call handle_error_netcdf('nf_put_var_real: east_e',status,
     $        nf_noerr)
      end if
            
      call merge2d(glb2d,north_u)
      if(my_task.eq.master_task .and. iint/iprint .eq. 1) then
         status=nf_inq_varid(ncid,"north_u",north_u_varid)
         call handle_error_netcdf('nf_inq_varid: north_u'
     $        ,status,nf_noerr)
         status=nf_put_var_real(ncid,north_u_varid,glb2d)
         call handle_error_netcdf('nf_put_var_real: north_u',status,
     $        nf_noerr)
      end if
      
      call merge2d(glb2d,north_v)
      if(my_task.eq.master_task .and. iint/iprint .eq. 1) then
         status=nf_inq_varid(ncid,"north_v",north_v_varid)
         call handle_error_netcdf('nf_inq_varid: north_v'
     $        ,status,nf_noerr)
         status=nf_put_var_real(ncid,north_v_varid,glb2d)
         call handle_error_netcdf('nf_put_var_real: north_v',status,
     $        nf_noerr)
      end if
      
      call merge2d(glb2d,north_e)
      if(my_task.eq.master_task .and. iint/iprint .eq. 1) then
         status=nf_inq_varid(ncid,"north_e",north_e_varid)
         call handle_error_netcdf('nf_inq_varid: north_e'
     $        ,status,nf_noerr)
         status=nf_put_var_real(ncid,north_e_varid,glb2d)
         call handle_error_netcdf('nf_put_var_real: north_e',status,
     $        nf_noerr)
      end if
            
      call merge2d(glb2d,h)
      if(my_task.eq.master_task .and. iint/iprint .eq. 1) then
         status=nf_inq_varid(ncid,"h",h_varid)
      call handle_error_netcdf('nf_inq_varid: h'
     $                          ,status,nf_noerr)
         status=nf_put_var_real(ncid,h_varid,glb2d)
         call handle_error_netcdf('nf_put_var_real: h',status,nf_noerr)
      end if
      
      call merge2d(glb2d,fsm)
      if(my_task.eq.master_task .and. iint/iprint .eq. 1) then
         status=nf_inq_varid(ncid,"fsm",fsm_varid)
         call handle_error_netcdf('nf_inq_varid: fsm'
     $        ,status,nf_noerr)
         status=nf_put_var_real(ncid,fsm_varid,glb2d)
         call handle_error_netcdf
     $        ('nf_put_var_real: fsm',status,nf_noerr)
      end if
      
      call merge2d(glb2d,dum)
      if(my_task.eq.master_task .and. iint/iprint .eq. 1) then
         status=nf_inq_varid(ncid,"dum",dum_varid)
         call handle_error_netcdf('nf_inq_varid: dum'
     $        ,status,nf_noerr)
         status=nf_put_var_real(ncid,dum_varid,glb2d)
         call handle_error_netcdf
     $        ('nf_put_var_real: dum',status,nf_noerr)
      end if
      
      call merge2d(glb2d,dvm)
      if(my_task.eq.master_task .and. iint/iprint .eq. 1) then
         status=nf_inq_varid(ncid,"dvm",dvm_varid)
         call handle_error_netcdf('nf_inq_varid: dvm'
     $        ,status,nf_noerr)
         status=nf_put_var_real(ncid,dvm_varid,glb2d)
         call handle_error_netcdf
     $        ('nf_put_var_real: dvm',status,nf_noerr)
      end if
            
      if(ldave == 1)then
         call merge2d(glb2d,el_dave)
      else if(ldave == 2)then
         call merge2d(glb2d,el)
      end if
      if(my_task.eq.master_task) then
         status=nf_inq_varid(ncid,"el",el_varid)
         call handle_error_netcdf('nf_inq_varid: el'
     $        ,status,nf_noerr)
         status=nf_put_vara_real(ncid,el_varid,
     $        (/1,1,nprint/),(/im_global,jm_global,1/),glb2d)
         call handle_error_netcdf('nf_put_var_real: el',status,
     $        nf_noerr)
      end if

      if(ldave == 1)then
         call merge2d(glb2d,lhf_dave)
      else if(ldave == 2)then
         call merge2d(glb2d,lhf)
      end if
      if(my_task.eq.master_task) then
         status=nf_inq_varid(ncid,"lhf",lhf_varid)
         call handle_error_netcdf('nf_inq_varid: lhf'
     $        ,status,nf_noerr)
         status=nf_put_vara_real(ncid,lhf_varid,
     $        (/1,1,nprint/),(/im_global,jm_global,1/),glb2d)
         call handle_error_netcdf('nf_put_var_real: lhf',status,
     $        nf_noerr)
      end if

      if(ldave == 1)then
         call merge2d(glb2d,shf_dave)
      else if(ldave == 2)then
         call merge2d(glb2d,shf)
      end if
      if(my_task.eq.master_task) then
         status=nf_inq_varid(ncid,"shf",shf_varid)
         call handle_error_netcdf('nf_inq_varid: shf'
     $        ,status,nf_noerr)
         status=nf_put_vara_real(ncid,shf_varid,
     $        (/1,1,nprint/),(/im_global,jm_global,1/),glb2d)

         call handle_error_netcdf('nf_put_var_real: shf',status,
     $        nf_noerr)
      end if

      if(ldave == 1)then
         call merge2d(glb2d,lwr_dave)
      else if(ldave == 2)then
         call merge2d(glb2d,lwr)
      end if
      if(my_task.eq.master_task) then
         status=nf_inq_varid(ncid,"lwr",lwr_varid)
         call handle_error_netcdf('nf_inq_varid: lwr'
     $        ,status,nf_noerr)
         status=nf_put_vara_real(ncid,lwr_varid,
     $        (/1,1,nprint/),(/im_global,jm_global,1/),glb2d)
         call handle_error_netcdf('nf_put_var_real: lwr',status,
     $        nf_noerr)
      end if

      if(ldave == 1)then
         call merge2d(glb2d,swr_dave)
      else if(ldave == 2)then
         call merge2d(glb2d,swr)
      end if
      if(my_task.eq.master_task) then
         status=nf_inq_varid(ncid,"swr",swr_varid)
         call handle_error_netcdf('nf_inq_varid: swr'
     $        ,status,nf_noerr)
         status=nf_put_vara_real(ncid,swr_varid,
     $        (/1,1,nprint/),(/im_global,jm_global,1/),glb2d)
         call handle_error_netcdf('nf_put_var_real: swr',status,
     $        nf_noerr)
      end if

      if(ldave == 1)then
         call merge2d(glb2d,windu_dave)
      else if(ldave == 2)then
         call merge2d(glb2d,windu)
      end if
      if(my_task.eq.master_task) then
         status=nf_inq_varid(ncid,"windu",windu_varid)
         call handle_error_netcdf('nf_inq_varid: windu'
     $        ,status,nf_noerr)
         status=nf_put_vara_real(ncid,windu_varid,
     $        (/1,1,nprint/),(/im_global,jm_global,1/),glb2d)
         call handle_error_netcdf('nf_put_var_real: windu',status,
     $        nf_noerr)
      end if

      if(ldave == 1)then
         call merge2d(glb2d,windv_dave)
      else if(ldave == 2)then
         call merge2d(glb2d,windv)
      end if
      if(my_task.eq.master_task) then
         status=nf_inq_varid(ncid,"windv",windv_varid)
         call handle_error_netcdf('nf_inq_varid: windv'
     $        ,status,nf_noerr)
         status=nf_put_vara_real(ncid,windv_varid,
     $        (/1,1,nprint/),(/im_global,jm_global,1/),glb2d)
         call handle_error_netcdf('nf_put_var_real: windv',status,
     $        nf_noerr)
      end if

      if(ldave == 1)then
         call merge2d(glb2d,winds_dave)
      else if(ldave == 2)then
         call merge2d(glb2d,winds)
      end if
      if(my_task.eq.master_task) then
         status=nf_inq_varid(ncid,"winds",winds_varid)
         call handle_error_netcdf('nf_inq_varid: winds'
     $        ,status,nf_noerr)
         status=nf_put_vara_real(ncid,winds_varid,
     $        (/1,1,nprint/),(/im_global,jm_global,1/),glb2d)
         call handle_error_netcdf('nf_put_var_real: winds',status,
     $        nf_noerr)
      end if

      if(ldave == 1)then
         call merge2d(glb2d,tauu_dave)
      else if(ldave == 2)then
         call merge2d(glb2d,tauu)
      end if
      if(my_task.eq.master_task) then
         status=nf_inq_varid(ncid,"tauu",tauu_varid)
         call handle_error_netcdf('nf_inq_varid: tauu'
     $        ,status,nf_noerr)
         status=nf_put_vara_real(ncid,tauu_varid,
     $        (/1,1,nprint/),(/im_global,jm_global,1/),glb2d)
         call handle_error_netcdf('nf_put_var_real: tauu',status,
     $        nf_noerr)
      end if

      if(ldave == 1)then
         call merge2d(glb2d,tauv_dave)
      else if(ldave == 2)then
         call merge2d(glb2d,tauv)
      end if
      if(my_task.eq.master_task) then
         status=nf_inq_varid(ncid,"tauv",tauv_varid)
         call handle_error_netcdf('nf_inq_varid: tauv'
     $        ,status,nf_noerr)
         status=nf_put_vara_real(ncid,tauv_varid,
     $        (/1,1,nprint/),(/im_global,jm_global,1/),glb2d)
         call handle_error_netcdf('nf_put_var_real: tauv',status,
     $        nf_noerr)
      end if

      if(ldave == 1)then
         call merge2d(glb2d,taus_dave)
      else if(ldave == 2)then
         call merge2d(glb2d,taus)
      end if
      if(my_task.eq.master_task) then
         status=nf_inq_varid(ncid,"taus",taus_varid)
         call handle_error_netcdf('nf_inq_varid: taus'
     $        ,status,nf_noerr)
         status=nf_put_vara_real(ncid,taus_varid,
     $        (/1,1,nprint/),(/im_global,jm_global,1/),glb2d)
         call handle_error_netcdf('nf_put_var_real: taus',status,
     $        nf_noerr)
      end if

      if(ldave == 1)then
         call merge2d(glb2d,qa_dave)
      else if(ldave == 2)then
         call merge2d(glb2d,qa)
      end if
      if(my_task.eq.master_task) then
         status=nf_inq_varid(ncid,"qa",qa_varid)
         call handle_error_netcdf('nf_inq_varid: qa'
     $        ,status,nf_noerr)
         status=nf_put_vara_real(ncid,qa_varid,
     $        (/1,1,nprint/),(/im_global,jm_global,1/),glb2d)
         call handle_error_netcdf('nf_put_var_real: qa',status,
     $        nf_noerr)
      end if

      if(ldave == 1)then
         call merge2d(glb2d,qs_dave)
      else if(ldave == 2)then
         call merge2d(glb2d,qs)
      end if
      if(my_task.eq.master_task) then
         status=nf_inq_varid(ncid,"qs",qs_varid)
         call handle_error_netcdf('nf_inq_varid: qs'
     $        ,status,nf_noerr)
         status=nf_put_vara_real(ncid,qs_varid,
     $        (/1,1,nprint/),(/im_global,jm_global,1/),glb2d)
         call handle_error_netcdf('nf_put_var_real: qs',status,
     $        nf_noerr)
      end if

      if(ldave == 1)then
         call merge2d(glb2d,ta_dave)
      else if(ldave == 2)then
         call merge2d(glb2d,ta)
      end if
      if(my_task.eq.master_task) then
         status=nf_inq_varid(ncid,"ta",ta_varid)
         call handle_error_netcdf('nf_inq_varid: ta'
     $        ,status,nf_noerr)
         status=nf_put_vara_real(ncid,ta_varid,
     $        (/1,1,nprint/),(/im_global,jm_global,1/),glb2d)
         call handle_error_netcdf('nf_put_var_real: ta',status,
     $        nf_noerr)
      end if

      if(lssf .eq. 1)then

         if(ldave == 1)then
            call merge2d(glb2d,evap_dave)
         else if(ldave == 2)then
            call merge2d(glb2d,evap)
         end if
         if(my_task.eq.master_task) then
            status=nf_inq_varid(ncid,"evap",evap_varid)
            call handle_error_netcdf('nf_inq_varid: evap'
     $           ,status,nf_noerr)
            status=nf_put_vara_real(ncid,evap_varid,
     $           (/1,1,nprint/),(/im_global,jm_global,1/),glb2d)
            call handle_error_netcdf('nf_put_var_real: evap',status,
     $           nf_noerr)
         end if

         if(ldave == 1)then
            call merge2d(glb2d,prep_dave)
         else if(ldave == 2)then
            call merge2d(glb2d,prep)
         end if
         if(my_task.eq.master_task) then
            status=nf_inq_varid(ncid,"prep",prep_varid)
            call handle_error_netcdf('nf_inq_varid: prep'
     $           ,status,nf_noerr)
            status=nf_put_vara_real(ncid,prep_varid,
     $           (/1,1,nprint/),(/im_global,jm_global,1/),glb2d)
            call handle_error_netcdf('nf_put_var_real: prep',status,
     $           nf_noerr)
         end if

         if(ldave == 1)then
            call merge2d(glb2d,river_dave)
         else if(ldave == 2)then
            call merge2d(glb2d,river)
         end if
         if(my_task.eq.master_task) then
            status=nf_inq_varid(ncid,"river",river_varid)
            call handle_error_netcdf('nf_inq_varid: river'
     $           ,status,nf_noerr)
            status=nf_put_vara_real(ncid,river_varid,
     $           (/1,1,nprint/),(/im_global,jm_global,1/),glb2d)
            call handle_error_netcdf('nf_put_var_real: river',status,
     $           nf_noerr)
         end if

         if(ldave == 1)then
            call merge2d(glb2d,eflux_dave)
         else if(ldave == 2)then
            call merge2d(glb2d,eflux)
         end if
         if(my_task.eq.master_task) then
            status=nf_inq_varid(ncid,"eflux",eflux_varid)
            call handle_error_netcdf('nf_inq_varid: eflux'
     $           ,status,nf_noerr)
            status=nf_put_vara_real(ncid,eflux_varid,
     $           (/1,1,nprint/),(/im_global,jm_global,1/),glb2d)
            call handle_error_netcdf('nf_put_var_real: eflux',status,
     $           nf_noerr)
         end if

         if(ldave == 1)then
            call merge2d(glb2d,pflux_dave)
         else if(ldave == 2)then
            call merge2d(glb2d,pflux)
         end if
         if(my_task.eq.master_task) then
            status=nf_inq_varid(ncid,"pflux",pflux_varid)
            call handle_error_netcdf('nf_inq_varid: pflux'
     $           ,status,nf_noerr)
            status=nf_put_vara_real(ncid,pflux_varid,
     $           (/1,1,nprint/),(/im_global,jm_global,1/),glb2d)
            call handle_error_netcdf('nf_put_var_real: pflux',status,
     $           nf_noerr)
         end if

         if(ldave == 1)then
            call merge2d(glb2d,rflux_dave)
         else if(ldave == 2)then
            call merge2d(glb2d,rflux)
         end if
         if(my_task.eq.master_task) then
            status=nf_inq_varid(ncid,"rflux",rflux_varid)
            call handle_error_netcdf('nf_inq_varid: rflux'
     $           ,status,nf_noerr)
            status=nf_put_vara_real(ncid,rflux_varid,
     $           (/1,1,nprint/),(/im_global,jm_global,1/),glb2d)
            call handle_error_netcdf('nf_put_var_real: rflux',status,
     $           nf_noerr)
         end if

      end if

      if(budget .eq. 1)then

         if(ldave == 1)then
            call merge2d(glb2d,tsfc_dave)
         else if(ldave == 2)then
            call merge2d(glb2d,tsfc)
         end if
         if(my_task.eq.master_task)then
            status=nf_inq_varid(ncid,"tsfc",tsfc_varid)
            call handle_error_netcdf('nf_inq_varid: tsfc'
     $           ,status,nf_noerr)
            status=nf_put_vara_real(ncid,tsfc_varid,
     $           (/1,1,nprint/),(/im_global,jm_global,1/),glb2d)
            call handle_error_netcdf('nf_put_var_real: tsfc',status,
     $           nf_noerr)
         end if

         if(ldave == 1)then
            call merge2d(glb2d,ssfc_dave)
         else if(ldave == 2)then
            call merge2d(glb2d,ssfc)
         end if
         if(my_task.eq.master_task)then
            status=nf_inq_varid(ncid,"ssfc",ssfc_varid)
            call handle_error_netcdf('nf_inq_varid: ssfc'
     $           ,status,nf_noerr)
            status=nf_put_vara_real(ncid,ssfc_varid,
     $           (/1,1,nprint/),(/im_global,jm_global,1/),glb2d)
            call handle_error_netcdf('nf_put_var_real: ssfc',status,
     $           nf_noerr)
         end if

      end if

      if(ldave == 1)then
         call merge3d(glb3d,u_dave)
      else if(ldave == 2)then
         call merge3d(glb3d,u)
      end if
      if(my_task.eq.master_task)then
         status=nf_inq_varid(ncid,"u",u_varid)
         call handle_error_netcdf('nf_inq_varid: u'
     $        ,status,nf_noerr)
         status=nf_put_vara_real(ncid,u_varid,
     $        (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
         call handle_error_netcdf('nf_put_var_real: u',status,nf_noerr)
      end if

      if(ldave == 1)then
         call merge3d(glb3d,v_dave)
      else if(ldave == 2)then
         call merge3d(glb3d,v)
      end if
      if(my_task.eq.master_task)then
         status=nf_inq_varid(ncid,"v",v_varid)
         call handle_error_netcdf('nf_inq_varid: v'
     $        ,status,nf_noerr)
         status=nf_put_vara_real(ncid,v_varid,
     $        (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
         call handle_error_netcdf('nf_put_var_real: v',status,nf_noerr)
      end if

      if(ldave == 1)then
         call merge3d(glb3d,wr_dave)
      else if(ldave == 2)then
         call merge3d(glb3d,wr)
      end if
      if(my_task.eq.master_task)then
         status=nf_inq_varid(ncid,"w",w_varid)
         call handle_error_netcdf('nf_inq_varid: w'
     $        ,status,nf_noerr)
         status=nf_put_vara_real(ncid,w_varid,
     $        (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
         call handle_error_netcdf('nf_put_var_real: w',status,nf_noerr)
      end if

      if(ldave == 1)then
         call merge3d(glb3d,t_dave)
      else if(ldave == 2)then
         call merge3d(glb3d,t)
      end if
      if(my_task.eq.master_task)then
         status=nf_inq_varid(ncid,"t",t_varid)
         call handle_error_netcdf('nf_inq_varid: t'
     $        ,status,nf_noerr)
         status=nf_put_vara_real(ncid,t_varid,
     $        (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
         call handle_error_netcdf('nf_put_var_real: t',status,nf_noerr)
      end if

      if(ldave == 1)then
         call merge3d(glb3d,s_dave)
      else if(ldave == 2)then
         call merge3d(glb3d,s)
      end if
      if(my_task.eq.master_task)then
         status=nf_inq_varid(ncid,"s",s_varid)
         call handle_error_netcdf('nf_inq_varid: s'
     $        ,status,nf_noerr)
         status=nf_put_vara_real(ncid,s_varid,
     $        (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
         call handle_error_netcdf('nf_put_var_real: s',status,nf_noerr)
      end if

!     Convervation
      if(budget == 1)then

!     T tendency
         if(ldave == 1)then
            call merge3d(glb3d,dtdt_dave)
         else if(ldave == 2)then
            call merge3d(glb3d,dtdt)
         end if
         if(my_task.eq.master_task)then
            status=nf_inq_varid(ncid,"dtdt",dtdt_varid)
            call handle_error_netcdf('nf_inq_varid: dtdt'
     $           ,status,nf_noerr)
            status=nf_put_vara_real(ncid,dtdt_varid,
     $           (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
            call handle_error_netcdf('nf_put_var_real: dtdt',
     $           status,nf_noerr)
         end if

!     T Advection
         if(ldave == 1)then
            call merge3d(glb3d,txadv_dave)
         else if(ldave == 2)then
            call merge3d(glb3d,txadv)
         end if
         if(my_task.eq.master_task)then
            status=nf_inq_varid(ncid,"txadv",txadv_varid)
            call handle_error_netcdf('nf_inq_varid: txadv'
     $           ,status,nf_noerr)
            status=nf_put_vara_real(ncid,txadv_varid,
     $           (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
            call handle_error_netcdf('nf_put_var_real: txadv',
     $           status,nf_noerr)
         end if

         if(ldave == 1)then
            call merge3d(glb3d,tyadv_dave)
         else if(ldave == 2)then
            call merge3d(glb3d,tyadv)
         end if
         if(my_task.eq.master_task)then
            status=nf_inq_varid(ncid,"tyadv",tyadv_varid)
            call handle_error_netcdf('nf_inq_varid: tyadv'
     $           ,status,nf_noerr)
            status=nf_put_vara_real(ncid,tyadv_varid,
     $           (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
            call handle_error_netcdf('nf_put_var_real: tyadv',
     $           status,nf_noerr)
         end if

         if(ldave == 1)then
            call merge3d(glb3d,tzadv_dave)
         else if(ldave == 2)then
            call merge3d(glb3d,tzadv)
         end if
         if(my_task.eq.master_task)then
            status=nf_inq_varid(ncid,"tzadv",tzadv_varid)
            call handle_error_netcdf('nf_inq_varid: tzadv'
     $           ,status,nf_noerr)
            status=nf_put_vara_real(ncid,tzadv_varid,
     $           (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
            call handle_error_netcdf('nf_put_var_real: tzadv',
     $           status,nf_noerr)
         end if

         if(ldave == 1)then
            call merge3d(glb3d,tadv_dave)
         else if(ldave == 2)then
            call merge3d(glb3d,tadv)
         end if
         if(my_task.eq.master_task)then
            status=nf_inq_varid(ncid,"tadv",tadv_varid)
            call handle_error_netcdf('nf_inq_varid: tadv'
     $           ,status,nf_noerr)
            status=nf_put_vara_real(ncid,tadv_varid,
     $           (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
            call handle_error_netcdf('nf_put_var_real: tadv',
     $           status,nf_noerr)
         end if

!     T Diffusion
         if(ldave == 1)then
            call merge3d(glb3d,txdif_dave)
         else if(ldave == 2)then
            call merge3d(glb3d,txdif)
         end if
         if(my_task.eq.master_task)then
            status=nf_inq_varid(ncid,"txdif",txdif_varid)
            call handle_error_netcdf('nf_inq_varid: txdif'
     $           ,status,nf_noerr)
            status=nf_put_vara_real(ncid,txdif_varid,
     $           (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
            call handle_error_netcdf('nf_put_var_real: txdif',
     $           status,nf_noerr)
         end if

         if(ldave == 1)then
            call merge3d(glb3d,tydif_dave)
         else if(ldave == 2)then
            call merge3d(glb3d,tydif)
         end if
         if(my_task.eq.master_task)then
            status=nf_inq_varid(ncid,"tydif",tydif_varid)
            call handle_error_netcdf('nf_inq_varid: tydif'
     $           ,status,nf_noerr)
            status=nf_put_vara_real(ncid,tydif_varid,
     $           (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
            call handle_error_netcdf('nf_put_var_real: tydif',
     $           status,nf_noerr)
         end if

         if(ldave == 1)then
            call merge3d(glb3d,tzdif_dave)
         else if(ldave == 2)then
            call merge3d(glb3d,tzdif)
         end if
         if(my_task.eq.master_task)then
            status=nf_inq_varid(ncid,"tzdif",tzdif_varid)
            call handle_error_netcdf('nf_inq_varid: tzdif'
     $           ,status,nf_noerr)
            status=nf_put_vara_real(ncid,tzdif_varid,
     $           (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
            call handle_error_netcdf('nf_put_var_real: tzdif',
     $           status,nf_noerr)
         end if

         if(ldave == 1)then
            call merge3d(glb3d,qz_dave)
         else if(ldave == 2)then
            call merge3d(glb3d,qz)
         end if
         if(my_task.eq.master_task)then
            status=nf_inq_varid(ncid,"qz",qz_varid)
            call handle_error_netcdf('nf_inq_varid: qz'
     $           ,status,nf_noerr)
            status=nf_put_vara_real(ncid,qz_varid,
     $           (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
            call handle_error_netcdf('nf_put_var_real: qz',
     $           status,nf_noerr)
         end if

!     S tendency
         if(ldave == 1)then
            call merge3d(glb3d,dsdt_dave)
         else if(ldave == 2)then
            call merge3d(glb3d,dsdt)
         end if
         if(my_task.eq.master_task)then
            status=nf_inq_varid(ncid,"dsdt",dsdt_varid)
            call handle_error_netcdf('nf_inq_varid: dsdt'
     $           ,status,nf_noerr)
            status=nf_put_vara_real(ncid,dsdt_varid,
     $           (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
            call handle_error_netcdf('nf_put_var_real: dsdt',
     $           status,nf_noerr)
         end if

!     S Advection
         if(ldave == 1)then
            call merge3d(glb3d,sxadv_dave)
         else if(ldave == 2)then
            call merge3d(glb3d,sxadv)
         end if
         if(my_task.eq.master_task)then
            status=nf_inq_varid(ncid,"sxadv",sxadv_varid)
            call handle_error_netcdf('nf_inq_varid: sxadv'
     $        ,status,nf_noerr)
            status=nf_put_vara_real(ncid,sxadv_varid,
     $           (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
            call handle_error_netcdf('nf_put_var_real: sxadv',
     $           status,nf_noerr)
         end if

         if(ldave == 1)then
            call merge3d(glb3d,syadv_dave)
         else if(ldave == 2)then
            call merge3d(glb3d,syadv)
         end if
         if(my_task.eq.master_task)then
            status=nf_inq_varid(ncid,"syadv",syadv_varid)
            call handle_error_netcdf('nf_inq_varid: syadv'
     $           ,status,nf_noerr)
            status=nf_put_vara_real(ncid,syadv_varid,
     $           (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
            call handle_error_netcdf('nf_put_var_real: syadv',
     $           status,nf_noerr)
         end if

         if(ldave == 1)then
            call merge3d(glb3d,szadv_dave)
         else if(ldave == 2)then
            call merge3d(glb3d,szadv)
         end if
         if(my_task.eq.master_task)then
            status=nf_inq_varid(ncid,"szadv",szadv_varid)
            call handle_error_netcdf('nf_inq_varid: szadv'
     $           ,status,nf_noerr)
            status=nf_put_vara_real(ncid,szadv_varid,
     $           (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
            call handle_error_netcdf('nf_put_var_real: szadv',
     $           status,nf_noerr)
         end if

         if(ldave == 1)then
            call merge3d(glb3d,sadv_dave)
         else if(ldave == 2)then
            call merge3d(glb3d,sadv)
         end if
         if(my_task.eq.master_task)then
            status=nf_inq_varid(ncid,"sadv",sadv_varid)
            call handle_error_netcdf('nf_inq_varid: sadv'
     $           ,status,nf_noerr)
            status=nf_put_vara_real(ncid,sadv_varid,
     $           (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
            call handle_error_netcdf('nf_put_var_real: sadv',
     $           status,nf_noerr)
         end if

!     S Diffusion
         if(ldave == 1)then
            call merge3d(glb3d,sxdif_dave)
         else if(ldave == 2)then
            call merge3d(glb3d,sxdif)
         end if
         if(my_task.eq.master_task)then
            status=nf_inq_varid(ncid,"sxdif",sxdif_varid)
            call handle_error_netcdf('nf_inq_varid: sxdif'
     $           ,status,nf_noerr)
            status=nf_put_vara_real(ncid,sxdif_varid,
     $           (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
            call handle_error_netcdf('nf_put_var_real: sxdif',
     $           status,nf_noerr)
         end if

         if(ldave == 1)then
            call merge3d(glb3d,sydif_dave)
         else if(ldave == 2)then
            call merge3d(glb3d,sydif)
         end if
         if(my_task.eq.master_task)then
            status=nf_inq_varid(ncid,"sydif",sydif_varid)
            call handle_error_netcdf('nf_inq_varid: sydif'
     $           ,status,nf_noerr)
            status=nf_put_vara_real(ncid,sydif_varid,
     $           (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
            call handle_error_netcdf('nf_put_var_real: sydif',
     $           status,nf_noerr)
         end if

         if(ldave == 1)then
            call merge3d(glb3d,szdif_dave)
         else if(ldave == 2)then
            call merge3d(glb3d,szdif)
         end if
         if(my_task.eq.master_task)then
            status=nf_inq_varid(ncid,"szdif",szdif_varid)
            call handle_error_netcdf('nf_inq_varid: szdif'
     $           ,status,nf_noerr)
            status=nf_put_vara_real(ncid,szdif_varid,
     $           (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
            call handle_error_netcdf('nf_put_var_real: szdif',
     $           status,nf_noerr)
         end if

!     AAM
         if(ldave == 1)then
            call merge3d(glb3d,aam_dave)
         else if(ldave == 2)then
            call merge3d(glb3d,aam)
         end if
         if(my_task.eq.master_task)then
            status=nf_inq_varid(ncid,"aam",aam_varid)
            call handle_error_netcdf('nf_inq_varid: aam'
     $           ,status,nf_noerr)
            status=nf_put_vara_real(ncid,aam_varid,
     $           (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
            call handle_error_netcdf('nf_put_var_real: aam',
     $           status,nf_noerr)
         end if

!     KH
         if(ldave == 1)then
            call merge3d(glb3d,kh_dave)
         else if(ldave == 2)then
            call merge3d(glb3d,kh)
         end if
         if(my_task.eq.master_task)then
            status=nf_inq_varid(ncid,"kh",kh_varid)
            call handle_error_netcdf('nf_inq_varid: kh'
     $           ,status,nf_noerr)
            status=nf_put_vara_real(ncid,kh_varid,
     $           (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
            call handle_error_netcdf('nf_put_var_real: kh',
     $           status,nf_noerr)
         end if

      end if

!     IAU
      if(budget == 1 .and. assim == 1)then

!         do k=1,kb
!            do j=1,jm
!               do i=1,im
!                  t_iau(i,j,k)=0.5e0*t_iau(i,j,k)
!               end do
!            end do
!         end do

         call merge3d(glb3d,t_iau)

         if(my_task.eq.master_task)then
            status=nf_inq_varid(ncid,"tiau",tiau_varid)
            call handle_error_netcdf('nf_inq_varid: tiau'
     $           ,status,nf_noerr)
            status=nf_put_vara_real(ncid,tiau_varid,
     $           (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
            call handle_error_netcdf('nf_put_var_real: tiau',
     $           status,nf_noerr)
         end if

!         do k=1,kb
!            do j=1,jm
!               do i=1,im
!                  s_iau(i,j,k)=0.5e0*s_iau(i,j,k)
!               end do
!            end do
!         end do

         call merge3d(glb3d,s_iau)

         if(my_task.eq.master_task)then
            status=nf_inq_varid(ncid,"siau",siau_varid)
            call handle_error_netcdf('nf_inq_varid: siau'
     $           ,status,nf_noerr)
            status=nf_put_vara_real(ncid,siau_varid,
     $           (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
            call handle_error_netcdf('nf_put_var_real: siau',
     $           status,nf_noerr)
         end if

      end if

!     Rounf off
      if(budget == 1 .and. lroff)then

         if(ldave == 1)then
            call merge3d(glb3d,troff_dave)
         else if(ldave == 2)then
            call merge3d(glb3d,troff)
         end if
         if(my_task.eq.master_task)then
            status=nf_inq_varid(ncid,"troff",troff_varid)
            call handle_error_netcdf('nf_inq_varid: troff'
     $           ,status,nf_noerr)
            status=nf_put_vara_real(ncid,troff_varid,
     $           (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
            call handle_error_netcdf('nf_put_var_real: troff',
     $           status,nf_noerr)
         end if

         if(ldave == 1)then
            call merge3d(glb3d,sroff_dave)
         else if(ldave == 2)then
            call merge3d(glb3d,sroff)
         end if
         if(my_task.eq.master_task)then
            status=nf_inq_varid(ncid,"sroff",sroff_varid)
            call handle_error_netcdf('nf_inq_varid: sroff'
     $           ,status,nf_noerr)
            status=nf_put_vara_real(ncid,sroff_varid,
     $           (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
            call handle_error_netcdf('nf_put_var_real: sroff',
     $           status,nf_noerr)
         end if

      end if

!     Nuding
      if(budget == 1)then

         if(ldave == 1)then
            call merge3d(glb3d,tnudge_dave)
         else if(ldave == 2)then
            call merge3d(glb3d,tnudge)
         end if

         if(my_task.eq.master_task)then
            status=nf_inq_varid(ncid,"tnudge",tnudge_varid)
            call handle_error_netcdf('nf_inq_varid: tnudge'
     $           ,status,nf_noerr)
            status=nf_put_vara_real(ncid,tnudge_varid,
     $           (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
            call handle_error_netcdf('nf_put_var_real: tnudge',
     $           status,nf_noerr)
         end if

      end if

      if(budget == 1)then

         if(ldave == 1)then
            call merge3d(glb3d,snudge_dave)
         else if(ldave == 2)then
            call merge3d(glb3d,snudge)
         end if

         if(my_task.eq.master_task)then
            status=nf_inq_varid(ncid,"snudge",snudge_varid)
            call handle_error_netcdf('nf_inq_varid: snudge'
     $           ,status,nf_noerr)
            status=nf_put_vara_real(ncid,snudge_varid,
     $           (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
            call handle_error_netcdf('nf_put_var_real: snudge',
     $           status,nf_noerr)
         end if

      end if

!     close file
      if(my_task.eq.master_task)then
         status=nf_close(ncid)
         call handle_error_netcdf('nf_close: output',status,nf_noerr)
      end if

      return
      end

!_______________________________________________________________________
      subroutine write_restart_netcdf
! write data for a seamless restart
      implicit none
# include "mpif.h"
# include "netcdf.inc"
      include 'pom.h'
      character*120 netcdf_out_file,str_tmp
      integer nprint
      integer ncid,status
      integer time_dimid,x_dimid,y_dimid,z_dimid
      integer time_varid,wubot_varid,wvbot_varid,aam2d_varid,ua_varid,
     $        uab_varid,va_varid,vab_varid,el_varid,elb_varid,et_varid,
     $        etb_varid,egb_varid,utb_varid,vtb_varid,u_varid,ub_varid,
     $        w_varid,v_varid,vb_varid,t_varid,tb_varid,s_varid,
     $        sb_varid,rho_varid,adx2d_varid,ady2d_varid,advua_varid,
     $        advva_varid,km_varid,kh_varid,kq_varid,l_varid,q2_varid,
     $        q2b_varid,aam_varid,q2l_varid,q2lb_varid
      integer vdims(4)
      integer length
      real glb2d(im_global,jm_global),glb3d(im_global,jm_global,kb)

      if(my_task.eq.master_task) then
         
!     create netcdf restart file
!     nprint=(iint+time0*86400/dti)/irestart
!     write(netcdf_out_file,'(''out/'',a,''.'',i5.5,''.nc'')') 
!     $                                       trim(write_rst_file),nprint
         
         nprint=nint(iint/real(irestart))
         netcdf_out_file="out/"//trim(write_rst_file)//".nc"
         write(*,'(/''writing file '',a)') trim(netcdf_out_file)
         
      end if
      
!     debug
      if(my_task .eq. master_task .and. nprint .eq. 1)then
         write(6,*) "Start to set definition"
      endif

      if(my_task.eq.master_task .and. nprint .eq. 1)then      

         status=nf_create(netcdf_out_file, 
     $        ior(nf_noclobber,nf_64bit_offset), ncid)
         call handle_error_netcdf('nf_create',status,nf_noerr)
         
!     define global attributes
         length=len(trim(title))
         status=nf_put_att_text(ncid,nf_global,'title',length,
     $        trim(title))
         call handle_error_netcdf('nf_put_att_text',status,nf_noerr)
         
         str_tmp='restart file'
         length=len(trim(str_tmp))
         status=nf_put_att_text(ncid,nf_global,'description',length,
     $        trim(str_tmp))
         call handle_error_netcdf('nf_put_att_text',status,nf_noerr)
         
!     define dimensions
         status=nf_def_dim(ncid,'time',nf_unlimited,time_dimid)
         call handle_error_netcdf('nf_def_dim',status,nf_noerr)
         
         length=kb
         status=nf_def_dim(ncid,'z',length,z_dimid)
         call handle_error_netcdf('nf_def_dim',status,nf_noerr)
         
         length=jm_global
         status=nf_def_dim(ncid,'y',length,y_dimid)
         call handle_error_netcdf('nf_def_dim',status,nf_noerr)
         
         length=im_global
         status=nf_def_dim(ncid,'x',length,x_dimid)
         call handle_error_netcdf('nf_def_dim',status,nf_noerr)
         
!     define variables and their attributes:
         vdims(1)=time_dimid
         str_tmp='days since '//time_start
         call def_var_netcdf(ncid,'time',1,vdims,time_varid,
     $        'time',str_tmp,' ',.false.)
         
         vdims(1)=x_dimid
         vdims(2)=y_dimid
         vdims(3)=time_dimid
         call def_var_netcdf(ncid,'wubot',3,vdims,wubot_varid,
     $        'x-momentum flux at the bottom',
     $        'metre^2/sec^2',
     $        'east_u north_u',.true.)
         call def_var_netcdf(ncid,'wvbot',3,vdims,wvbot_varid,
     $        'y-momentum flux at the bottom',
     $        'metre^2/sec^2',
     $        'east_v north_v',.true.)
         call def_var_netcdf(ncid,'aam2d',3,vdims,aam2d_varid,
     $        'vertical average of aam',
     $        'metre^2/sec',
     $        'east_e north_e',.true.)
         call def_var_netcdf(ncid,'ua',3,vdims,ua_varid,
     $        'vertical mean of u',
     $        'metre/sec',
     $        'east_u north_u',.true.)
         call def_var_netcdf(ncid,'uab',3,vdims,uab_varid,
     $        'vertical mean of u at time -dt',
     $        'metre/sec',
     $        'east_u north_u',.true.)
         call def_var_netcdf(ncid,'va',3,vdims,va_varid,
     $        'vertical mean of v',
     $        'metre/sec',
     $        'east_v north_v',.true.)
         call def_var_netcdf(ncid,'vab',3,vdims,vab_varid,
     $        'vertical mean of v at time -dt',
     $        'metre/sec',
     $        'east_v north_v',.true.)
         call def_var_netcdf(ncid,'el',3,vdims,el_varid,
     $        'surface elevation in external mode',
     $        'metre',
     $        'east_e north_e',.true.)
         call def_var_netcdf(ncid,'elb',3,vdims,elb_varid,
     $        'surface elevation in external mode at -dt',
     $        'metre',
     $        'east_e north_e',.true.)
         call def_var_netcdf(ncid,'et',3,vdims,et_varid,
     $        'surface elevation in internal mode',
     $        'metre',
     $        'east_e north_e',.true.)
         call def_var_netcdf(ncid,'etb',3,vdims,etb_varid,
     $        'surface elevation in internal mode at -dt',
     $        'metre',
     $        'east_e north_e',.true.)
         call def_var_netcdf(ncid,'egb',3,vdims,egb_varid,
     $        'surface elevation for pres. grad. at -dt',
     $        'metre',
     $        'east_e north_e',.true.)
         call def_var_netcdf(ncid,'utb',3,vdims,utb_varid,
     $        'ua time averaged over dti',
     $        'metre/sec',
     $        'east_u north_u',.true.)
         call def_var_netcdf(ncid,'vtb',3,vdims,vtb_varid,
     $        'va time averaged over dti',
     $        'metre/sec',
     $        'east_v north_v',.true.)
         call def_var_netcdf(ncid,'adx2d',3,vdims,adx2d_varid,
     $        'vertical integral of advx',
     $        '-',
     $        'east_u north_u',.true.)
         call def_var_netcdf(ncid,'ady2d',3,vdims,ady2d_varid,
     $        'vertical integral of advy',
     $        '-',
     $        'east_v north_v',.true.)
         call def_var_netcdf(ncid,'advua',3,vdims,advua_varid,
     $        'sum of 2nd, 3rd and 4th terms in eq (18)',
     $        '-',
     $        'east_u north_u',.true.)
         call def_var_netcdf(ncid,'advva',3,vdims,advva_varid,
     $        'sum of 2nd, 3rd and 4th terms in eq (19)',
     $        '-',
     $        'east_v north_v',.true.)
         
         vdims(1)=x_dimid
         vdims(2)=y_dimid
         vdims(3)=z_dimid
         vdims(4)=time_dimid
         call def_var_netcdf(ncid,'u',4,vdims,u_varid,
     $        'x-velocity',
     $        'metre/sec',
     $        'east_u north_u zz',.true.)
         call def_var_netcdf(ncid,'ub',4,vdims,ub_varid,
     $        'x-velocity at time -dt',
     $        'metre/sec',
     $        'east_u north_u zz',.true.)
         call def_var_netcdf(ncid,'v',4,vdims,v_varid,
     $        'y-velocity',
     $        'metre/sec',
     $        'east_v north_v zz',.true.)
         call def_var_netcdf(ncid,'vb',4,vdims,vb_varid,
     $        'y-velocity at time -dt',
     $        'metre/sec',
     $        'east_v north_v zz',.true.)
         call def_var_netcdf(ncid,'w',4,vdims,w_varid,
     $        'sigma-velocity',
     $        'metre/sec',
     $        'east_e north_e zz',.true.)
         call def_var_netcdf(ncid,'t',4,vdims,t_varid,
     $        'potential temperature',
     $        'degree C',
     $        'east_e north_e zz',.true.)
         call def_var_netcdf(ncid,'tb',4,vdims,tb_varid,
     $        'potential temperature at time -dt',
     $        'degree C',
     $        'east_e north_e zz',.true.)
         call def_var_netcdf(ncid,'s',4,vdims,s_varid,
     $        'salinity x rho / rhoref',
     $        'PSS',
     $        'east_e north_e zz',.true.)
         call def_var_netcdf(ncid,'sb',4,vdims,sb_varid,
     $        'salinity x rho / rhoref at time -dt',
     $        'PSS',
     $        'east_e north_e zz',.true.)
         call def_var_netcdf(ncid,'rho',4,vdims,rho_varid,
     $        '(density-1000)/rhoref',
     $        'dimensionless',
     $        'east_e north_e zz',.true.)
         call def_var_netcdf(ncid,'km',4,vdims,km_varid,
     $        'vertical kinematic viscosity',
     $        'metre^2/sec',
     $        'east_e north_e zz',.true.)
         call def_var_netcdf(ncid,'kh',4,vdims,kh_varid,
     $        'vertical diffusivity',
     $        'metre^2/sec',
     $        'east_e north_e zz',.true.)
         call def_var_netcdf(ncid,'kq',4,vdims,kq_varid,
     $        'kq',
     $        'metre^2/sec',
     $        'east_e north_e zz',.true.)
         call def_var_netcdf(ncid,'l',4,vdims,l_varid,
     $        'turbulence length scale',
     $        '-',
     $        'east_e north_e zz',.true.)
         call def_var_netcdf(ncid,'q2',4,vdims,q2_varid,
     $        'twice the turbulent kinetic energy',
     $        'metre^2/sec^2',
     $        'east_e north_e zz',.true.)
         call def_var_netcdf(ncid,'q2b',4,vdims,q2b_varid,
     $        'twice the turbulent kinetic energy at -dt',
     $        'metre^2/sec^2',
     $        'east_e north_e zz',.true.)
         call def_var_netcdf(ncid,'aam',4,vdims,aam_varid,
     $        'horizontal kinematic viscosity',
     $        'metre^2/sec',
     $        'east_e north_e zz',.true.)
         call def_var_netcdf(ncid,'q2l',4,vdims,q2l_varid,
     $        'q2 x l',
     $        'metre^3/sec^2',
     $        'east_e north_e zz',.true.)
         call def_var_netcdf(ncid,'q2lb',4,vdims,q2lb_varid,
     $        'q2 x l at time -dt',
     $        'metre^3/sec^2',
     $        'east_e north_e zz',.true.)
                  
!     end definitions
         status=nf_enddef(ncid)
         call handle_error_netcdf('nf_enddef',status,nf_noerr)
         
      end if

!     debug
      if(my_task .eq. master_task .and. nprint .eq. 1)then
         write(6,*) "Finish to set definition"
      endif

      if(my_task.eq.master_task .and. nprint .ne. 1) then      
         
         status=nf_open(netcdf_out_file,nf_write,ncid)
         call handle_error_netcdf('nf_open',status,nf_noerr)
         
      end if

!     debug
      if(my_task .eq. master_task) write(6,*) "Start to set write"

! write data

!     debug
      if(my_task .eq. master_task) write(6,*) "Write 2d"

      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,time_varid,(/nprint/),(/1/),time)
        call handle_error_netcdf('nf_put_vara_real',status,nf_noerr)
      end if

      call merge2d(glb2d,wubot)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,wubot_varid,
     $        (/1,1,nprint/),(/im_global,jm_global,1/),glb2d)
        call handle_error_netcdf('nf_put_vara_real',status,nf_noerr)
      end if

      call merge2d(glb2d,wvbot)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,wvbot_varid,
     $        (/1,1,nprint/),(/im_global,jm_global,1/),glb2d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

      call merge2d(glb2d,aam2d)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,aam2d_varid,
     $        (/1,1,nprint/),(/im_global,jm_global,1/),glb2d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

      call merge2d(glb2d,ua)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,ua_varid,
     $        (/1,1,nprint/),(/im_global,jm_global,1/),glb2d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

      call merge2d(glb2d,uab)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,uab_varid,
     $        (/1,1,nprint/),(/im_global,jm_global,1/),glb2d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

      call merge2d(glb2d,va)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,va_varid,
     $        (/1,1,nprint/),(/im_global,jm_global,1/),glb2d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

      call merge2d(glb2d,vab)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,vab_varid,
     $        (/1,1,nprint/),(/im_global,jm_global,1/),glb2d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

      call merge2d(glb2d,el)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,el_varid,
     $        (/1,1,nprint/),(/im_global,jm_global,1/),glb2d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

      call merge2d(glb2d,elb)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,elb_varid,
     $        (/1,1,nprint/),(/im_global,jm_global,1/),glb2d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

      call merge2d(glb2d,et)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,et_varid,
     $        (/1,1,nprint/),(/im_global,jm_global,1/),glb2d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

      call merge2d(glb2d,etb)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,etb_varid,
     $        (/1,1,nprint/),(/im_global,jm_global,1/),glb2d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

      call merge2d(glb2d,egb)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,egb_varid,
     $        (/1,1,nprint/),(/im_global,jm_global,1/),glb2d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

      call merge2d(glb2d,utb)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,utb_varid,
     $        (/1,1,nprint/),(/im_global,jm_global,1/),glb2d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

      call merge2d(glb2d,vtb)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,vtb_varid,
     $        (/1,1,nprint/),(/im_global,jm_global,1/),glb2d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

      call merge2d(glb2d,adx2d)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,adx2d_varid,
     $        (/1,1,nprint/),(/im_global,jm_global,1/),glb2d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

      call merge2d(glb2d,ady2d)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,ady2d_varid,
     $        (/1,1,nprint/),(/im_global,jm_global,1/),glb2d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

      call merge2d(glb2d,advua)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,advua_varid,
     $        (/1,1,nprint/),(/im_global,jm_global,1/),glb2d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

      call merge2d(glb2d,advva)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,advva_varid,
     $        (/1,1,nprint/),(/im_global,jm_global,1/),glb2d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

!     debug
      if(my_task .eq. master_task) write(6,*) "Write 3d u"

      call merge3d(glb3d,u)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,u_varid,
     $        (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

      call merge3d(glb3d,ub)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,ub_varid,
     $        (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

!     debug
      if(my_task .eq. master_task) write(6,*) "Write 3d v"

      call merge3d(glb3d,v)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,v_varid,
     $        (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

      call merge3d(glb3d,vb)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,vb_varid,
     $        (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

!     debug
      if(my_task .eq. master_task) write(6,*) "Write 3d w"

      call merge3d(glb3d,w)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,w_varid,
     $        (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

!     debug
      if(my_task .eq. master_task) write(6,*) "Write 3d t"

      call merge3d(glb3d,t)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,t_varid,
     $        (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

      call merge3d(glb3d,tb)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,tb_varid,
     $        (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

!     debug
      if(my_task .eq. master_task) write(6,*) "Write 3d s"

      call merge3d(glb3d,s)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,s_varid,
     $        (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

      call merge3d(glb3d,sb)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,sb_varid,
     $        (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

!     debug
      if(my_task .eq. master_task) write(6,*) "Write 3d rho"

      call merge3d(glb3d,rho)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,rho_varid,
     $        (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

!     debug
      if(my_task .eq. master_task) write(6,*) "Write 3d km"

      call merge3d(glb3d,km)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,km_varid,
     $        (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

!     debug
      if(my_task .eq. master_task) write(6,*) "Write 3d kh"

      call merge3d(glb3d,kh)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,kh_varid,
     $        (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

!     debug
      if(my_task .eq. master_task) write(6,*) "Write 3d kq"

      call merge3d(glb3d,kq)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,kq_varid,
     $        (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

!     debug
      if(my_task .eq. master_task) write(6,*) "Write 3d l"

      call merge3d(glb3d,l)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,l_varid,
     $        (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

!     debug
      if(my_task .eq. master_task) write(6,*) "Write 3d q"

      call merge3d(glb3d,q2)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,q2_varid,
     $        (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

      call merge3d(glb3d,q2b)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,q2b_varid,
     $        (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

      call merge3d(glb3d,aam)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,aam_varid,
     $        (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

      call merge3d(glb3d,q2l)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,q2l_varid,
     $        (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

      call merge3d(glb3d,q2lb)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,q2lb_varid,
     $        (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

! close file
      if(my_task.eq.master_task) then
        status=nf_close(ncid)
        call handle_error_netcdf('nf_close',status,nf_noerr)
      end if

!     debug
      if(my_task .eq. master_task) write(6,*) "Finish to set write"

      return
      end
!_______________________________________________________________________
      subroutine write_ens_netcdf
! write data for a seamless restart
! Created by S.Ohishi
! almost same as write_restart_netcdf, but for filename
      implicit none
# include "mpif.h"
# include "netcdf.inc"
      include 'pom.h'
      character*120 netcdf_out_file,str_tmp
      integer nprint
      integer ncid,status
      integer time_dimid,x_dimid,y_dimid,z_dimid
      integer time_varid,wubot_varid,wvbot_varid,aam2d_varid,ua_varid,
     $        uab_varid,va_varid,vab_varid,el_varid,elb_varid,et_varid,
     $        etb_varid,egb_varid,utb_varid,vtb_varid,u_varid,ub_varid,
     $        w_varid,v_varid,vb_varid,t_varid,tb_varid,s_varid,
     $        sb_varid,rho_varid,adx2d_varid,ady2d_varid,advua_varid,
     $        advva_varid,km_varid,kh_varid,kq_varid,l_varid,q2_varid,
     $        q2b_varid,aam_varid,q2l_varid,q2lb_varid
      integer vdims(4)
      integer length
      real glb2d(im_global,jm_global),glb3d(im_global,jm_global,kb)

      if(my_task.eq.master_task)then
         
!     create netcdf restart file
!     nprint=(iint+time0*86400/dti)/irestart
!     write(netcdf_out_file,'(''out/'',a,''.'',i5.5,''.nc'')') 
!     $                                       trim(write_rst_file),nprint
         
         nprint=1
         write(netcdf_out_file,'(''out/'',a,''.'',i2.2,''.nc'')') 
     $        trim(write_ens_file),nint(iint/real(iprint))
         write(*,'(/''writing file '',a)') trim(netcdf_out_file)
         
      end if
      
      if(my_task.eq.master_task .and. nprint .eq. 1)then      

         status=nf_create(netcdf_out_file, 
     $        ior(nf_noclobber,nf_64bit_offset), ncid)
         call handle_error_netcdf('nf_create',status,nf_noerr)
         
!     define global attributes
         length=len(trim(title))
         status=nf_put_att_text(ncid,nf_global,'title',length,
     $        trim(title))
         call handle_error_netcdf('nf_put_att_text',status,nf_noerr)
         
         str_tmp='restart file'
         length=len(trim(str_tmp))
         status=nf_put_att_text(ncid,nf_global,'description',length,
     $        trim(str_tmp))
         call handle_error_netcdf('nf_put_att_text',status,nf_noerr)
         
!     define dimensions
         status=nf_def_dim(ncid,'time',nf_unlimited,time_dimid)
         call handle_error_netcdf('nf_def_dim',status,nf_noerr)
         
         length=kb
         status=nf_def_dim(ncid,'z',length,z_dimid)
         call handle_error_netcdf('nf_def_dim',status,nf_noerr)
         
         length=jm_global
         status=nf_def_dim(ncid,'y',length,y_dimid)
         call handle_error_netcdf('nf_def_dim',status,nf_noerr)
         
         length=im_global
         status=nf_def_dim(ncid,'x',length,x_dimid)
         call handle_error_netcdf('nf_def_dim',status,nf_noerr)
         
!     define variables and their attributes:
         vdims(1)=time_dimid
         str_tmp='days since '//time_start
         call def_var_netcdf(ncid,'time',1,vdims,time_varid,
     $        'time',str_tmp,' ',.false.)
         
         vdims(1)=x_dimid
         vdims(2)=y_dimid
         vdims(3)=time_dimid
         call def_var_netcdf(ncid,'wubot',3,vdims,wubot_varid,
     $        'x-momentum flux at the bottom',
     $        'metre^2/sec^2',
     $        'east_u north_u',.true.)
         call def_var_netcdf(ncid,'wvbot',3,vdims,wvbot_varid,
     $        'y-momentum flux at the bottom',
     $        'metre^2/sec^2',
     $        'east_v north_v',.true.)
         call def_var_netcdf(ncid,'aam2d',3,vdims,aam2d_varid,
     $        'vertical average of aam',
     $        'metre^2/sec',
     $        'east_e north_e',.true.)
         call def_var_netcdf(ncid,'ua',3,vdims,ua_varid,
     $        'vertical mean of u',
     $        'metre/sec',
     $        'east_u north_u',.true.)
         call def_var_netcdf(ncid,'uab',3,vdims,uab_varid,
     $        'vertical mean of u at time -dt',
     $        'metre/sec',
     $        'east_u north_u',.true.)
         call def_var_netcdf(ncid,'va',3,vdims,va_varid,
     $        'vertical mean of v',
     $        'metre/sec',
     $        'east_v north_v',.true.)
         call def_var_netcdf(ncid,'vab',3,vdims,vab_varid,
     $        'vertical mean of v at time -dt',
     $        'metre/sec',
     $        'east_v north_v',.true.)
         call def_var_netcdf(ncid,'el',3,vdims,el_varid,
     $        'surface elevation in external mode',
     $        'metre',
     $        'east_e north_e',.true.)
         call def_var_netcdf(ncid,'elb',3,vdims,elb_varid,
     $        'surface elevation in external mode at -dt',
     $        'metre',
     $        'east_e north_e',.true.)
         call def_var_netcdf(ncid,'et',3,vdims,et_varid,
     $        'surface elevation in internal mode',
     $        'metre',
     $        'east_e north_e',.true.)
         call def_var_netcdf(ncid,'etb',3,vdims,etb_varid,
     $        'surface elevation in internal mode at -dt',
     $        'metre',
     $        'east_e north_e',.true.)
         call def_var_netcdf(ncid,'egb',3,vdims,egb_varid,
     $        'surface elevation for pres. grad. at -dt',
     $        'metre',
     $        'east_e north_e',.true.)
         call def_var_netcdf(ncid,'utb',3,vdims,utb_varid,
     $        'ua time averaged over dti',
     $        'metre/sec',
     $        'east_u north_u',.true.)
         call def_var_netcdf(ncid,'vtb',3,vdims,vtb_varid,
     $        'va time averaged over dti',
     $        'metre/sec',
     $        'east_v north_v',.true.)
         call def_var_netcdf(ncid,'adx2d',3,vdims,adx2d_varid,
     $        'vertical integral of advx',
     $        '-',
     $        'east_u north_u',.true.)
         call def_var_netcdf(ncid,'ady2d',3,vdims,ady2d_varid,
     $        'vertical integral of advy',
     $        '-',
     $        'east_v north_v',.true.)
         call def_var_netcdf(ncid,'advua',3,vdims,advua_varid,
     $        'sum of 2nd, 3rd and 4th terms in eq (18)',
     $        '-',
     $        'east_u north_u',.true.)
         call def_var_netcdf(ncid,'advva',3,vdims,advva_varid,
     $        'sum of 2nd, 3rd and 4th terms in eq (19)',
     $        '-',
     $        'east_v north_v',.true.)
         
         vdims(1)=x_dimid
         vdims(2)=y_dimid
         vdims(3)=z_dimid
         vdims(4)=time_dimid
         call def_var_netcdf(ncid,'u',4,vdims,u_varid,
     $        'x-velocity',
     $        'metre/sec',
     $        'east_u north_u zz',.true.)
         call def_var_netcdf(ncid,'ub',4,vdims,ub_varid,
     $        'x-velocity at time -dt',
     $        'metre/sec',
     $        'east_u north_u zz',.true.)
         call def_var_netcdf(ncid,'v',4,vdims,v_varid,
     $        'y-velocity',
     $        'metre/sec',
     $        'east_v north_v zz',.true.)
         call def_var_netcdf(ncid,'vb',4,vdims,vb_varid,
     $        'y-velocity at time -dt',
     $        'metre/sec',
     $        'east_v north_v zz',.true.)
         call def_var_netcdf(ncid,'w',4,vdims,w_varid,
     $        'sigma-velocity',
     $        'metre/sec',
     $        'east_e north_e zz',.true.)
         call def_var_netcdf(ncid,'t',4,vdims,t_varid,
     $        'potential temperature',
     $        'degree C',
     $        'east_e north_e zz',.true.)
         call def_var_netcdf(ncid,'tb',4,vdims,tb_varid,
     $        'potential temperature at time -dt',
     $        'degree C',
     $        'east_e north_e zz',.true.)
         call def_var_netcdf(ncid,'s',4,vdims,s_varid,
     $        'salinity x rho / rhoref',
     $        'PSS',
     $        'east_e north_e zz',.true.)
         call def_var_netcdf(ncid,'sb',4,vdims,sb_varid,
     $        'salinity x rho / rhoref at time -dt',
     $        'PSS',
     $        'east_e north_e zz',.true.)
         call def_var_netcdf(ncid,'rho',4,vdims,rho_varid,
     $        '(density-1000)/rhoref',
     $        'dimensionless',
     $        'east_e north_e zz',.true.)
         call def_var_netcdf(ncid,'km',4,vdims,km_varid,
     $        'vertical kinematic viscosity',
     $        'metre^2/sec',
     $        'east_e north_e zz',.true.)
         call def_var_netcdf(ncid,'kh',4,vdims,kh_varid,
     $        'vertical diffusivity',
     $        'metre^2/sec',
     $        'east_e north_e zz',.true.)
         call def_var_netcdf(ncid,'kq',4,vdims,kq_varid,
     $        'kq',
     $        'metre^2/sec',
     $        'east_e north_e zz',.true.)
         call def_var_netcdf(ncid,'l',4,vdims,l_varid,
     $        'turbulence length scale',
     $        '-',
     $        'east_e north_e zz',.true.)
         call def_var_netcdf(ncid,'q2',4,vdims,q2_varid,
     $        'twice the turbulent kinetic energy',
     $        'metre^2/sec^2',
     $        'east_e north_e zz',.true.)
         call def_var_netcdf(ncid,'q2b',4,vdims,q2b_varid,
     $        'twice the turbulent kinetic energy at -dt',
     $        'metre^2/sec^2',
     $        'east_e north_e zz',.true.)
         call def_var_netcdf(ncid,'aam',4,vdims,aam_varid,
     $        'horizontal kinematic viscosity',
     $        'metre^2/sec',
     $        'east_e north_e zz',.true.)
         call def_var_netcdf(ncid,'q2l',4,vdims,q2l_varid,
     $        'q2 x l',
     $        'metre^3/sec^2',
     $        'east_e north_e zz',.true.)
         call def_var_netcdf(ncid,'q2lb',4,vdims,q2lb_varid,
     $        'q2 x l at time -dt',
     $        'metre^3/sec^2',
     $        'east_e north_e zz',.true.)
                  
!     end definitions
         status=nf_enddef(ncid)
         call handle_error_netcdf('nf_enddef',status,nf_noerr)
         
      end if

      if(my_task.eq.master_task .and. nprint .ne. 1) then      
         
         status=nf_open(netcdf_out_file,nf_write,ncid)
         call handle_error_netcdf('nf_open',status,nf_noerr)
         
      end if

! write data

      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,time_varid,(/nprint/),(/1/),time)
        call handle_error_netcdf('nf_put_vara_real',status,nf_noerr)
      end if

      call merge2d(glb2d,wubot)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,wubot_varid,
     $        (/1,1,nprint/),(/im_global,jm_global,1/),glb2d)
        call handle_error_netcdf('nf_put_vara_real',status,nf_noerr)
      end if

      call merge2d(glb2d,wvbot)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,wvbot_varid,
     $        (/1,1,nprint/),(/im_global,jm_global,1/),glb2d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

      call merge2d(glb2d,aam2d)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,aam2d_varid,
     $        (/1,1,nprint/),(/im_global,jm_global,1/),glb2d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

      call merge2d(glb2d,ua)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,ua_varid,
     $        (/1,1,nprint/),(/im_global,jm_global,1/),glb2d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

      call merge2d(glb2d,uab)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,uab_varid,
     $        (/1,1,nprint/),(/im_global,jm_global,1/),glb2d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

      call merge2d(glb2d,va)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,va_varid,
     $        (/1,1,nprint/),(/im_global,jm_global,1/),glb2d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

      call merge2d(glb2d,vab)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,vab_varid,
     $        (/1,1,nprint/),(/im_global,jm_global,1/),glb2d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

      call merge2d(glb2d,el)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,el_varid,
     $        (/1,1,nprint/),(/im_global,jm_global,1/),glb2d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

      call merge2d(glb2d,elb)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,elb_varid,
     $        (/1,1,nprint/),(/im_global,jm_global,1/),glb2d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

      call merge2d(glb2d,et)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,et_varid,
     $        (/1,1,nprint/),(/im_global,jm_global,1/),glb2d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

      call merge2d(glb2d,etb)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,etb_varid,
     $        (/1,1,nprint/),(/im_global,jm_global,1/),glb2d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

      call merge2d(glb2d,egb)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,egb_varid,
     $        (/1,1,nprint/),(/im_global,jm_global,1/),glb2d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

      call merge2d(glb2d,utb)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,utb_varid,
     $        (/1,1,nprint/),(/im_global,jm_global,1/),glb2d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

      call merge2d(glb2d,vtb)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,vtb_varid,
     $        (/1,1,nprint/),(/im_global,jm_global,1/),glb2d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

      call merge2d(glb2d,adx2d)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,adx2d_varid,
     $        (/1,1,nprint/),(/im_global,jm_global,1/),glb2d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

      call merge2d(glb2d,ady2d)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,ady2d_varid,
     $        (/1,1,nprint/),(/im_global,jm_global,1/),glb2d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

      call merge2d(glb2d,advua)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,advua_varid,
     $        (/1,1,nprint/),(/im_global,jm_global,1/),glb2d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

      call merge2d(glb2d,advva)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,advva_varid,
     $        (/1,1,nprint/),(/im_global,jm_global,1/),glb2d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

      call merge3d(glb3d,u)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,u_varid,
     $        (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

      call merge3d(glb3d,ub)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,ub_varid,
     $        (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

      call merge3d(glb3d,v)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,v_varid,
     $        (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

      call merge3d(glb3d,vb)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,vb_varid,
     $        (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

      call merge3d(glb3d,w)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,w_varid,
     $        (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

      call merge3d(glb3d,t)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,t_varid,
     $        (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

      call merge3d(glb3d,tb)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,tb_varid,
     $        (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

      call merge3d(glb3d,s)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,s_varid,
     $        (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

      call merge3d(glb3d,sb)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,sb_varid,
     $        (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

      call merge3d(glb3d,rho)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,rho_varid,
     $        (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

      call merge3d(glb3d,km)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,km_varid,
     $        (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

      call merge3d(glb3d,kh)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,kh_varid,
     $        (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

      call merge3d(glb3d,kq)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,kq_varid,
     $        (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

      call merge3d(glb3d,l)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,l_varid,
     $        (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

      call merge3d(glb3d,q2)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,q2_varid,
     $        (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

      call merge3d(glb3d,q2b)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,q2b_varid,
     $        (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

      call merge3d(glb3d,aam)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,aam_varid,
     $        (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

      call merge3d(glb3d,q2l)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,q2l_varid,
     $        (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

      call merge3d(glb3d,q2lb)
      if(my_task.eq.master_task) then
        status=nf_put_vara_real(ncid,q2lb_varid,
     $        (/1,1,1,nprint/),(/im_global,jm_global,kb,1/),glb3d)
        call handle_error_netcdf('nf_put_var_real',status,nf_noerr)
      end if

! close file
      if(my_task.eq.master_task) then
        status=nf_close(ncid)
        call handle_error_netcdf('nf_close',status,nf_noerr)
      end if

      return
      end

!_______________________________________________________________________
      subroutine read_grid_netcdf
! read grid data
      implicit none
# include "mpif.h"
# include "netcdf.inc"
      include 'pom.h'
      character*120 netcdf_grid_file
      integer z_varid,zz_varid,dx_varid,dy_varid,east_c_varid,
     $        east_e_varid,east_u_varid,east_v_varid,north_c_varid,
     $        north_e_varid,north_u_varid,north_v_varid,rot_varid,
     $        h_varid,fsm_varid,dum_varid,dvm_varid
      integer ncid,status
      integer vdims(2)
      integer length
      integer start(2),edge(2)
      real glb2d(im_global,jm_global),glb3d(im_global,jm_global,kb)
! debug
!      integer i,j
!      character*1 cha_my_task

! open netcdf file
      write(netcdf_grid_file,'(''in/'',a,''.grid.nc'')') 
     $                                                 trim(netcdf_file)
      write(*,'(/''reading file '',a)')
     $  trim(netcdf_grid_file)
      status=nf_open(netcdf_grid_file,nf_nowrite,ncid)
      call handle_error_netcdf('nf_open: '//netcdf_grid_file,
     $                          status,nf_noerr)
     
! get variables
      status=nf_inq_varid(ncid,'z',z_varid)
      call handle_error_netcdf('nf_inq_varid: z'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'zz',zz_varid)
      call handle_error_netcdf('nf_inq_varid: zz'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'dx',dx_varid)
      call handle_error_netcdf('nf_inq_varid: dx'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'dy',dy_varid)
      call handle_error_netcdf('nf_inq_varid: dy'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'east_u',east_u_varid)
      call handle_error_netcdf('nf_inq_varid: east_u'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'east_v',east_v_varid)
      call handle_error_netcdf('nf_inq_varid: east_v'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'east_e',east_e_varid)
      call handle_error_netcdf('nf_inq_varid: east_e'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'east_c',east_c_varid)
      call handle_error_netcdf('nf_inq_varid: east_c'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'north_u',north_u_varid)
      call handle_error_netcdf('nf_inq_varid: north_u'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'north_v',north_v_varid)
      call handle_error_netcdf('nf_inq_varid: north_v'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'north_e',north_e_varid)
      call handle_error_netcdf('nf_inq_varid: north_e'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'north_c',north_c_varid)
      call handle_error_netcdf('nf_inq_varid: north_c'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'rot',rot_varid)
      call handle_error_netcdf('nf_inq_varid: rot'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'h',h_varid)
      call handle_error_netcdf('nf_inq_varid: h'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'fsm',fsm_varid)
      call handle_error_netcdf('nf_inq_varid: fsm'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'dum',dum_varid)
      call handle_error_netcdf('nf_inq_varid: dum'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'dvm',dvm_varid)
      call handle_error_netcdf('nf_inq_varid: dvm'
     $                          ,status,nf_noerr)

! get data
      status=nf_get_var_real(ncid,z_varid,z)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      status=nf_get_var_real(ncid,zz_varid,zz)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)

      status=nf_get_var_real(ncid,dx_varid,glb2d)
      call handle_error_netcdf('nf_get_var_real ',
     $                          status,nf_noerr)
      call scatter2d(dx,glb2d)

!      write(cha_my_task,'(i1.1)') my_task
!      write(*,*) cha_my_task//'/output'
!      open(88,file=cha_my_task//'/output',form='formatted')
!      do j=1,jm
!        write(88,*) (dx(i,j),i=1,im)
!      end do
!      close(88)
!      call mpi_finalize
!      stop
      
      status=nf_get_var_real(ncid,dy_varid,glb2d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter2d(dy,glb2d)

      status=nf_get_var_real(ncid,east_u_varid,glb2d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter2d(east_u,glb2d)

      status=nf_get_var_real(ncid,east_v_varid,glb2d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter2d(east_v,glb2d)

      status=nf_get_var_real(ncid,east_e_varid,glb2d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter2d(east_e,glb2d)

      status=nf_get_var_real(ncid,east_c_varid,glb2d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter2d(east_c,glb2d)

      status=nf_get_var_real(ncid,north_u_varid,glb2d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter2d(north_u,glb2d)

      status=nf_get_var_real(ncid,north_v_varid,glb2d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter2d(north_v,glb2d)

      status=nf_get_var_real(ncid,north_e_varid,glb2d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter2d(north_e,glb2d)

      status=nf_get_var_real(ncid,north_c_varid,glb2d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter2d(north_c,glb2d)

      status=nf_get_var_real(ncid,rot_varid,glb2d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter2d(rot,glb2d)

      status=nf_get_var_real(ncid,h_varid,glb2d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter2d(h,glb2d)

      status=nf_get_var_real(ncid,fsm_varid,glb2d)
      call handle_error_netcdf('nf_get_var_real',
     $                         status,nf_noerr)
      call scatter2d(fsm,glb2d)

      status=nf_get_var_real(ncid,dum_varid,glb2d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter2d(dum,glb2d)

      status=nf_get_var_real(ncid,dvm_varid,glb2d)
      call handle_error_netcdf('nf_get_var_real',
     $                         status,nf_noerr)
      call scatter2d(dvm,glb2d)

! close file:
      status=nf_close(ncid)
      call handle_error_netcdf('nf_close: grid',status,nf_noerr)

      return
      end

!_______________________________________________________________________
      subroutine read_restart_netcdf
! read data for a seamless restart
      implicit none
# include "mpif.h"
# include "netcdf.inc"
      include 'pom.h'
      character*120 netcdf_in_file
      integer ncid,status
      integer time_varid,wubot_varid,wvbot_varid,aam2d_varid,ua_varid,
     $        uab_varid,va_varid,vab_varid,el_varid,elb_varid,et_varid,
     $        etb_varid,egb_varid,utb_varid,vtb_varid,u_varid,ub_varid,
     $        w_varid,v_varid,vb_varid,t_varid,tb_varid,s_varid,
     $        sb_varid,rho_varid,adx2d_varid,ady2d_varid,advua_varid,
     $        advva_varid,km_varid,kh_varid,kq_varid,l_varid,q2_varid,
     $        q2b_varid,aam_varid,q2l_varid,q2lb_varid
      integer vdims(3)
      integer length
      integer i,j,k
      real glb2d(im_global,jm_global),glb3d(im_global,jm_global,kb)

! open netcdf restart file

      write(netcdf_in_file,'(''in/'',a)') trim(read_rst_file)
      if(my_task.eq.master_task) write(*,'(/''reading file '',a)')
     $ trim(netcdf_in_file)
      status=nf_open(netcdf_in_file,nf_nowrite,ncid)
      call handle_error_netcdf('nf_open: '//trim(read_rst_file),
     $                          status,nf_noerr)

! get variables
      status=nf_inq_varid(ncid,'time',time_varid)
      call handle_error_netcdf('nf_inq_varid: time'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'wubot',wubot_varid)
      call handle_error_netcdf('nf_inq_varid: wubot'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'wvbot',wvbot_varid)
      call handle_error_netcdf('nf_inq_varid: wvbot'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'aam2d',aam2d_varid)
      call handle_error_netcdf('nf_inq_varid: aam2d'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'ua',ua_varid)
      call handle_error_netcdf('nf_inq_varid: ua'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'uab',uab_varid)
      call handle_error_netcdf('nf_inq_varid: uab'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'va',va_varid)
      call handle_error_netcdf('nf_inq_varid: va'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'vab',vab_varid)
      call handle_error_netcdf('nf_inq_varid: vab'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'el',el_varid)
      call handle_error_netcdf('nf_inq_varid: el'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'elb',elb_varid)
      call handle_error_netcdf('nf_inq_varid: elb'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'et',et_varid)
      call handle_error_netcdf('nf_inq_varid: et'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'etb',etb_varid)
      call handle_error_netcdf('nf_inq_varid: etb'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'egb',egb_varid)
      call handle_error_netcdf('nf_inq_varid: egb'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'utb',utb_varid)
      call handle_error_netcdf('nf_inq_varid: utb'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'vtb',vtb_varid)
      call handle_error_netcdf('nf_inq_varid: vtb'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'u',u_varid)
      call handle_error_netcdf('nf_inq_varid: u'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'ub',ub_varid)
      call handle_error_netcdf('nf_inq_varid: ub'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'v',v_varid)
      call handle_error_netcdf('nf_inq_varid: v'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'vb',vb_varid)
      call handle_error_netcdf('nf_inq_varid: vb'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'w',w_varid)
      call handle_error_netcdf('nf_inq_varid: w'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'t',t_varid)
      call handle_error_netcdf('nf_inq_varid: t'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'tb',tb_varid)
      call handle_error_netcdf('nf_inq_varid: tb'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'s',s_varid)
      call handle_error_netcdf('nf_inq_varid: s'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'sb',sb_varid)
      call handle_error_netcdf('nf_inq_varid: sb'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'rho',rho_varid)
      call handle_error_netcdf('nf_inq_varid: rho'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'adx2d',adx2d_varid)
      call handle_error_netcdf('nf_inq_varid: adx2d'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'ady2d',ady2d_varid)
      call handle_error_netcdf('nf_inq_varid: ady2d'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'advua',advua_varid)
      call handle_error_netcdf('nf_inq_varid: advua'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'advva',advva_varid)
      call handle_error_netcdf('nf_inq_varid: advva'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'km',km_varid)
      call handle_error_netcdf('nf_inq_varid: km'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'kh',kh_varid)
      call handle_error_netcdf('nf_inq_varid: kh'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'kq',kq_varid)
      call handle_error_netcdf('nf_inq_varid: kq'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'l',l_varid)
      call handle_error_netcdf('nf_inq_varid: l'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'q2',q2_varid)
      call handle_error_netcdf('nf_inq_varid: q2'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'q2b',q2b_varid)
      call handle_error_netcdf('nf_inq_varid: q2b'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'aam',aam_varid)
      call handle_error_netcdf('nf_inq_varid: aam'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'q2l',q2l_varid)
      call handle_error_netcdf('nf_inq_varid: q2l'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'q2lb',q2lb_varid)
      call handle_error_netcdf('nf_inq_varid: q2lb'
     $                          ,status,nf_noerr)

! get data
      status=nf_get_var_real(ncid,time_varid,time0)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)

      status=nf_get_var_real(ncid,wubot_varid,glb2d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter2d(wubot,glb2d)

      status=nf_get_var_real(ncid,wvbot_varid,glb2d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter2d(wvbot,glb2d)

      status=nf_get_var_real(ncid,aam2d_varid,glb2d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter2d(aam2d,glb2d)

      status=nf_get_var_real(ncid,ua_varid,glb2d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter2d(ua,glb2d)

      status=nf_get_var_real(ncid,uab_varid,glb2d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter2d(uab,glb2d)

      status=nf_get_var_real(ncid,va_varid,glb2d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter2d(va,glb2d)

      status=nf_get_var_real(ncid,vab_varid,glb2d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter2d(vab,glb2d)

      status=nf_get_var_real(ncid,el_varid,glb2d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter2d(el,glb2d)

      status=nf_get_var_real(ncid,elb_varid,glb2d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter2d(elb,glb2d)

      status=nf_get_var_real(ncid,et_varid,glb2d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter2d(et,glb2d)

      status=nf_get_var_real(ncid,etb_varid,glb2d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter2d(etb,glb2d)

      status=nf_get_var_real(ncid,egb_varid,glb2d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter2d(egb,glb2d)

      status=nf_get_var_real(ncid,utb_varid,glb2d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter2d(utb,glb2d)

      status=nf_get_var_real(ncid,vtb_varid,glb2d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter2d(vtb,glb2d)

      status=nf_get_var_real(ncid,adx2d_varid,glb2d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter2d(adx2d,glb2d)

      status=nf_get_var_real(ncid,ady2d_varid,glb2d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter2d(ady2d,glb2d)

      status=nf_get_var_real(ncid,advua_varid,glb2d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter2d(advua,glb2d)

      status=nf_get_var_real(ncid,advva_varid,glb2d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter2d(advva,glb2d)

      status=nf_get_var_real(ncid,u_varid,glb3d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter3d(u,glb3d)

      status=nf_get_var_real(ncid,ub_varid,glb3d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter3d(ub,glb3d)

      status=nf_get_var_real(ncid,v_varid,glb3d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter3d(v,glb3d)

      status=nf_get_var_real(ncid,vb_varid,glb3d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter3d(vb,glb3d)

      status=nf_get_var_real(ncid,w_varid,glb3d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter3d(w,glb3d)

      status=nf_get_var_real(ncid,t_varid,glb3d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter3d(t,glb3d)

      status=nf_get_var_real(ncid,tb_varid,glb3d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter3d(tb,glb3d)

      status=nf_get_var_real(ncid,s_varid,glb3d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter3d(s,glb3d)

      status=nf_get_var_real(ncid,sb_varid,glb3d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter3d(sb,glb3d)

      status=nf_get_var_real(ncid,rho_varid,glb3d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter3d(rho,glb3d)

      status=nf_get_var_real(ncid,km_varid,glb3d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter3d(km,glb3d)

      status=nf_get_var_real(ncid,kh_varid,glb3d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter3d(kh,glb3d)

      status=nf_get_var_real(ncid,kq_varid,glb3d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter3d(kq,glb3d)

      status=nf_get_var_real(ncid,l_varid,glb3d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter3d(l,glb3d)

      status=nf_get_var_real(ncid,q2_varid,glb3d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter3d(q2,glb3d)

      status=nf_get_var_real(ncid,q2b_varid,glb3d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter3d(q2b,glb3d)

      status=nf_get_var_real(ncid,aam_varid,glb3d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter3d(aam,glb3d)

      status=nf_get_var_real(ncid,q2l_varid,glb3d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter3d(q2l,glb3d)

      status=nf_get_var_real(ncid,q2lb_varid,glb3d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter3d(q2lb,glb3d)

! close file
      status=nf_close(ncid)
      call handle_error_netcdf('nf_close',status,nf_noerr)

!     for IAU
      if(assim == 3)then
         do j=1,jm
            do i=1,im
               ua(i,j)=uab(i,j)
               va(i,j)=vab(i,j)
               el(i,j)=elb(i,j)
            end do
         end do
         do k=1,kb
            do j=1,jm
               do i=1,im
                  u(i,j,k)=ub(i,j,k)
                  v(i,j,k)=vb(i,j,k)
                  t(i,j,k)=tb(i,j,k)
                  s(i,j,k)=sb(i,j,k)
               end do
            end do
         end do
      end if
      
! initialize elevation
      do j=1,jm
        do i=1,im
          d(i,j)=h(i,j)+el(i,j)
          dt(i,j)=h(i,j)+et(i,j)
        end do
      end do

! initialize time
      time=time0

      return
      end

!_______________________________________________________________________
      subroutine read_iau_netcdf
! read data for a seamless restart
      implicit none
# include "mpif.h"
# include "netcdf.inc"
      include 'pom.h'
      character*120 netcdf_in_file
      integer ncid,status
      integer ua_varid,uab_varid,va_varid,vab_varid,el_varid,elb_varid,
     $     u_varid,ub_varid,v_varid,vb_varid,
     $     t_varid,tb_varid,s_varid,sb_varid
      integer vdims(3)
      integer length
      integer i,j,k
      real glb2d(im_global,jm_global),glb3d(im_global,jm_global,kb)
      real ua_iau1(im,jm),ua_iau2(im,jm)
      real va_iau1(im,jm),va_iau2(im,jm)
      real el_iau1(im,jm),el_iau2(im,jm)
      real t_iau1(im,jm,kb),t_iau2(im,jm,kb)
      real s_iau1(im,jm,kb),s_iau2(im,jm,kb)
      real u_iau1(im,jm,kb),u_iau2(im,jm,kb)
      real v_iau1(im,jm,kb),v_iau2(im,jm,kb)

! open netcdf restart file

      write(netcdf_in_file,'(''in/'',a)') trim(read_iau_file)
      if(my_task.eq.master_task) write(*,'(/''reading file '',a)')
     $     trim(netcdf_in_file)
      status=nf_open(netcdf_in_file,nf_nowrite,ncid)
      call handle_error_netcdf('nf_open: '//trim(read_iau_file),
     $                          status,nf_noerr)

! get variables
      status=nf_inq_varid(ncid,'ua',ua_varid)
      call handle_error_netcdf('nf_inq_varid: ua'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'uab',uab_varid)
      call handle_error_netcdf('nf_inq_varid: uab'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'va',va_varid)
      call handle_error_netcdf('nf_inq_varid: va'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'vab',vab_varid)
      call handle_error_netcdf('nf_inq_varid: vab'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'el',el_varid)
      call handle_error_netcdf('nf_inq_varid: el'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'elb',elb_varid)
      call handle_error_netcdf('nf_inq_varid: elb'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'u',u_varid)
      call handle_error_netcdf('nf_inq_varid: u'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'ub',ub_varid)
      call handle_error_netcdf('nf_inq_varid: ub'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'v',v_varid)
      call handle_error_netcdf('nf_inq_varid: v'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'vb',vb_varid)
      call handle_error_netcdf('nf_inq_varid: vb'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'t',t_varid)
      call handle_error_netcdf('nf_inq_varid: t'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'tb',tb_varid)
      call handle_error_netcdf('nf_inq_varid: tb'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'s',s_varid)
      call handle_error_netcdf('nf_inq_varid: s'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'sb',sb_varid)
      call handle_error_netcdf('nf_inq_varid: sb'
     $                          ,status,nf_noerr)

! get data
      status=nf_get_var_real(ncid,ua_varid,glb2d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter2d(ua_iau1,glb2d)

      status=nf_get_var_real(ncid,uab_varid,glb2d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter2d(ua_iau2,glb2d)

      status=nf_get_var_real(ncid,va_varid,glb2d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter2d(va_iau1,glb2d)

      status=nf_get_var_real(ncid,vab_varid,glb2d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter2d(va_iau2,glb2d)

      status=nf_get_var_real(ncid,el_varid,glb2d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter2d(el_iau1,glb2d)

      status=nf_get_var_real(ncid,elb_varid,glb2d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter2d(el_iau2,glb2d)

      status=nf_get_var_real(ncid,u_varid,glb3d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter3d(u_iau1,glb3d)

      status=nf_get_var_real(ncid,ub_varid,glb3d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter3d(u_iau2,glb3d)

      status=nf_get_var_real(ncid,v_varid,glb3d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter3d(v_iau1,glb3d)

      status=nf_get_var_real(ncid,vb_varid,glb3d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter3d(v_iau2,glb3d)

      status=nf_get_var_real(ncid,t_varid,glb3d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter3d(t_iau1,glb3d)

      status=nf_get_var_real(ncid,tb_varid,glb3d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter3d(t_iau2,glb3d)

      status=nf_get_var_real(ncid,s_varid,glb3d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter3d(s_iau1,glb3d)

      status=nf_get_var_real(ncid,sb_varid,glb3d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter3d(s_iau2,glb3d)

! close file
      status=nf_close(ncid)
      call handle_error_netcdf('nf_close',status,nf_noerr)

!     Calculate variable_iau (anal-model)
      do j=1,jm
         do i=1,im
            ua_iau(i,j)=ua_iau2(i,j)-ua_iau1(i,j)
            va_iau(i,j)=va_iau2(i,j)-va_iau1(i,j)
            el_iau(i,j)=el_iau2(i,j)-el_iau1(i,j)
         end do
      end do

      do k=1,kb
         do j=1,jm
            do i=1,im
               t_iau(i,j,k)=t_iau2(i,j,k)-t_iau1(i,j,k)
               s_iau(i,j,k)=s_iau2(i,j,k)-s_iau1(i,j,k)
               u_iau(i,j,k)=u_iau2(i,j,k)-u_iau1(i,j,k)
               v_iau(i,j,k)=v_iau2(i,j,k)-v_iau1(i,j,k)
            end do
         end do
      end do

      return
      end

!_______________________________________________________________________
      subroutine read_initial_ts_netcdf(k,temp,salt)
! read initial temperature and salinity
      implicit none
# include "mpif.h"
# include "netcdf.inc"
      include 'pom.h'
      integer k
      real temp(im,jm,k),salt(im,jm,k)
      real glb3d(im_global,jm_global,kb)
      character*120 netcdf_ic_file
      integer ncid,status
      integer t_varid,s_varid
      integer vdims(3)
   
! open netcdf ic file
      write(netcdf_ic_file,'(''in/'',a,''.ic.nc'')') trim(netcdf_file)
      if(my_task.eq.master_task) write(*,'(/''reading file '',a)')
     $  trim(netcdf_ic_file)
      status=nf_open(netcdf_ic_file,nf_nowrite,ncid)
      call handle_error_netcdf('nf_open: '//netcdf_ic_file,
     $                         status,nf_noerr)

! get variables
      status=nf_inq_varid(ncid,'t',t_varid)
      call handle_error_netcdf('nf_inq_varid: t'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'s',s_varid)
      call handle_error_netcdf('nf_inq_varid: s'
     $                          ,status,nf_noerr)

! get data

      status=nf_get_var_real(ncid,t_varid,glb3d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter3d(t,glb3d)

      status=nf_get_var_real(ncid,s_varid,glb3d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter3d(s,glb3d)

! close file
      status=nf_close(ncid)
      call handle_error_netcdf('nf_close',status,nf_noerr)

      return
      end

!_______________________________________________________________________
      subroutine read_initial_tsuv_netcdf(k,temp,salt,xvel,yvel)
! read initial temperature, salinity, x-velocity, and y-velocity
      implicit none
# include "mpif.h"
# include "netcdf.inc"
      include 'pom.h'
      integer k
      real temp(im,jm,k),salt(im,jm,k)
      real xvel(im,jm,k),yvel(im,jm,k)
      real glb3d(im_global,jm_global,kb)
      character*120 netcdf_ic_file
      integer ncid,status
      integer t_varid,s_varid
      integer u_varid,v_varid
      integer vdims(3)
! debug
!      character*1 cha_my_task
!      integer i,j

! open netcdf ic file
      write(netcdf_ic_file,'(''in/'',a,''.ic.nc'')') trim(netcdf_file)
      if(my_task.eq.master_task) write(*,'(/''reading file '',a)')
     $  trim(netcdf_ic_file)
      status=nf_open(netcdf_ic_file,nf_nowrite,ncid)
      call handle_error_netcdf('nf_open: '//netcdf_ic_file,
     $                          status,nf_noerr)

! get variables
      status=nf_inq_varid(ncid,'t',t_varid)
      call handle_error_netcdf('nf_inq_varid: t'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'s',s_varid)
      call handle_error_netcdf('nf_inq_varid: s'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'u',u_varid)
      call handle_error_netcdf('nf_inq_varid: u'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'v',v_varid)
      call handle_error_netcdf('nf_inq_varid: v'
     $                          ,status,nf_noerr)

! get data
      status=nf_get_var_real(ncid,t_varid,glb3d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter3d(temp,glb3d)

      status=nf_get_var_real(ncid,s_varid,glb3d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter3d(salt,glb3d)

      status=nf_get_var_real(ncid,u_varid,glb3d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter3d(xvel,glb3d)

      status=nf_get_var_real(ncid,v_varid,glb3d)
      call handle_error_netcdf('nf_get_var_real',
     $                          status,nf_noerr)
      call scatter3d(yvel,glb3d)

! close file
      status=nf_close(ncid)
      call handle_error_netcdf('nf_close',status,nf_noerr)

      return
      end

!_______________________________________________________________________
      subroutine read_lbc_mclim_netcdf(month)
!     read monthly climatology of lateral boundary condition
!     Created by 2018.08 S.Ohishi
      
      implicit none
# include "mpif.h"
# include "netcdf.inc"
      include 'pom.h'

      integer month
      character*120 netcdf_lbc_file
      integer ncid,status,varid
      integer start(3),count(3)
      
      real glb2dew(jm_global),glb2dns(im_global)
      real glb3dew(jm_global,kb),glb3dns(im_global,kb)
      
!     open netcdf file
      write(netcdf_lbc_file,'(''in/'',a,''.lbc.nc'')') trim(netcdf_file)
      if(my_task .eq. master_task)then 
         write(*,'(/''reading file '',a)') trim(netcdf_lbc_file)
      endif
      status=nf_open(netcdf_lbc_file,nf_nowrite,ncid)
      call handle_error_netcdf('nf_open: '//netcdf_lbc_file,
     $     status,nf_noerr)

! check month

      if(my_task .eq. master_task) write(6,*) 'lbc month ',month
      if(month .le. 0 .or. month .gt. 12) then
         if(my_task .eq. master_task) write(*,'(/a)')
     $        'wrong month in read_lbc_mclim_netcdf'
         call finalize_mpi
         stop
      end if

!     specify 2d array shape at wester/eastern boudary

      start(1)=1
      count(1)=jm_global
      start(2)=month
      count(2)=1

!     eastern boudary
      status=nf_inq_varid(ncid,'ele',varid)
      call handle_error_netcdf('nf_inq_varid: ele',status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2dew)
      call handle_error_netcdf('nf_get_vara_real: ele',status,nf_noerr)
      call scatter_bnd_ew(ele0,glb2dew,1)
      
      status=nf_inq_varid(ncid,'uabe',varid)
      call handle_error_netcdf('nf_inq_varid: uabe',status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2dew)
      call handle_error_netcdf('nf_get_vara_real: uabe',status,nf_noerr)
      call scatter_bnd_ew(uabe0,glb2dew,1)
      
      status=nf_inq_varid(ncid,'vabe',varid)
      call handle_error_netcdf('nf_inq_varid: vabe',status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2dew)
      call handle_error_netcdf('nf_get_vara_real: vabe',status,nf_noerr)
      call scatter_bnd_ew(vabe0,glb2dew,1)
      
!     western boudary
      status=nf_inq_varid(ncid,'elw',varid)
      call handle_error_netcdf('nf_inq_varid: elw',status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2dew)
      call handle_error_netcdf('nf_get_vara_real: elw',status,nf_noerr)
      call scatter_bnd_ew(elw0,glb2dew,1)
      
      status=nf_inq_varid(ncid,'uabw',varid)
      call handle_error_netcdf('nf_inq_varid: uabw',status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2dew)
      call handle_error_netcdf('nf_get_vara_real: uabw',status,nf_noerr)
      call scatter_bnd_ew(uabw0,glb2dew,1)
      
      status=nf_inq_varid(ncid,'vabw',varid)
      call handle_error_netcdf('nf_inq_varid: vabw',status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2dew)
      call handle_error_netcdf('nf_get_vara_real: vabw',status,nf_noerr)
      call scatter_bnd_ew(vabw0,glb2dew,1)
      
!     specify 3d array shape at wester/eastern boudary
      
      start(1)=1
      count(1)=jm_global
      start(2)=1
      count(2)=kb
      start(3)=month
      count(3)=1
      
!     eastern boudary
      status=nf_inq_varid(ncid,'tbe',varid)
      call handle_error_netcdf('nf_inq_varid: tbe',status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dew)
      call handle_error_netcdf('nf_get_vara_real: tbe',status,nf_noerr)
      call scatter_bnd_ew(tbe0,glb3dew,kb)
      
      status=nf_inq_varid(ncid,'sbe',varid)
      call handle_error_netcdf('nf_inq_varid: sbe',status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dew)
      call handle_error_netcdf('nf_get_vara_real: sbe',status,nf_noerr)
      call scatter_bnd_ew(sbe0,glb3dew,kb)
      
      status=nf_inq_varid(ncid,'ube',varid)
      call handle_error_netcdf('nf_inq_varid: ube',status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dew)
      call handle_error_netcdf('nf_get_vara_real: ube',status,nf_noerr)
      call scatter_bnd_ew(ube0,glb3dew,kb)
      
      status=nf_inq_varid(ncid,'vbe',varid)
      call handle_error_netcdf('nf_inq_varid: vbe',status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dew)
      call handle_error_netcdf('nf_get_vara_real: vbe',status,nf_noerr)
      call scatter_bnd_ew(vbe0,glb3dew,kb)
      
!     western boudary
      status=nf_inq_varid(ncid,'tbw',varid)
      call handle_error_netcdf('nf_inq_varid: tbw',status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dew)
      call handle_error_netcdf('nf_get_vara_real: tbw',status,nf_noerr)
      call scatter_bnd_ew(tbw0,glb3dew,kb)
      
      status=nf_inq_varid(ncid,'sbw',varid)
      call handle_error_netcdf('nf_inq_varid: sbw',status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dew)
      call handle_error_netcdf('nf_get_vara_real: sbw',status,nf_noerr)
      call scatter_bnd_ew(sbw0,glb3dew,kb)
      
      status=nf_inq_varid(ncid,'ubw',varid)
      call handle_error_netcdf('nf_inq_varid: ubw',status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dew)
      call handle_error_netcdf('nf_get_vara_real: ubw',status,nf_noerr)
      call scatter_bnd_ew(ubw0,glb3dew,kb)
      
      status=nf_inq_varid(ncid,'vbw',varid)
      call handle_error_netcdf('nf_inq_varid: vbw',status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dew)
      call handle_error_netcdf('nf_get_vara_real: vbw',status,nf_noerr)
      call scatter_bnd_ew(vbw0,glb3dew,kb)
      
!     specify 2d array shape at northern/southern boudary
      
      start(1)=1
      count(1)=im_global
      start(2)=month
      count(2)=1
      
!     northern boudary
      status=nf_inq_varid(ncid,'eln',varid)
      call handle_error_netcdf('nf_inq_varid: eln',status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2dns)
      call handle_error_netcdf('nf_get_var_real: eln',status,nf_noerr)
      call scatter_bnd_ns(eln0,glb2dns,1)
      
      status=nf_inq_varid(ncid,'uabn',varid)
      call handle_error_netcdf('nf_inq_varid: uabn',status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2dns)
      call handle_error_netcdf('nf_get_var_real: uabn',status,nf_noerr)
      call scatter_bnd_ns(uabn0,glb2dns,1)
      
      status=nf_inq_varid(ncid,'vabn',varid)
      call handle_error_netcdf('nf_inq_varid: vabn',status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2dns)
      call handle_error_netcdf('nf_get_var_real: vabn',status,nf_noerr)
      call scatter_bnd_ns(vabn0,glb2dns,1)
      
!     southern boudary
      status=nf_inq_varid(ncid,'els',varid)
      call handle_error_netcdf('nf_inq_varid: els',status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2dns)
      call handle_error_netcdf('nf_get_var_real: els',status,nf_noerr)
      call scatter_bnd_ns(els0,glb2dns,1)
      
      status=nf_inq_varid(ncid,'uabs',varid)
      call handle_error_netcdf('nf_inq_varid: uabs',status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2dns)
      call handle_error_netcdf('nf_get_var_real: uabs',status,nf_noerr)
      call scatter_bnd_ns(uabs0,glb2dns,1)
      
      status=nf_inq_varid(ncid,'vabs',varid)
      call handle_error_netcdf('nf_inq_varid: vabs',status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2dns)
      call handle_error_netcdf('nf_get_var_real: vabs',status,nf_noerr)
      call scatter_bnd_ns(vabs0,glb2dns,1)
      
!     specify 3d array shape at northern/southern boudary
      
      start(1)=1
      count(1)=im_global
      start(2)=1
      count(2)=kb
      start(3)=month
      count(3)=1
      
!     northern boudary
      status=nf_inq_varid(ncid,'tbn',varid)
      call handle_error_netcdf('nf_inq_varid: tbn',status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dns)
      call handle_error_netcdf('nf_get_vara_real: tbn',status,nf_noerr)
      call scatter_bnd_ns(tbn0, glb3dns,kb)
      
      status=nf_inq_varid(ncid,'sbn',varid)
      call handle_error_netcdf('nf_inq_varid: sbn',status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dns)
      call handle_error_netcdf('nf_get_vara_real: sbn',status,nf_noerr)
      call scatter_bnd_ns(sbn0, glb3dns,kb)
      
      status=nf_inq_varid(ncid,'ubn',varid)
      call handle_error_netcdf('nf_inq_varid: ubn',status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dns)
      call handle_error_netcdf('nf_get_vara_real: ubn',status,nf_noerr)
      call scatter_bnd_ns(ubn0, glb3dns,kb)
      
      status=nf_inq_varid(ncid,'vbn',varid)
      call handle_error_netcdf('nf_inq_varid: vbn',status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dns)
      call handle_error_netcdf('nf_get_vara_real: vbn',status,nf_noerr)
      call scatter_bnd_ns(vbn0, glb3dns,kb)
      
!     southern boudary
      status=nf_inq_varid(ncid,'tbs',varid)
      call handle_error_netcdf('nf_inq_varid: tbs',status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dns)
      call handle_error_netcdf('nf_get_vara_real: tbs',status,nf_noerr)
      call scatter_bnd_ns(tbs0, glb3dns,kb)
      
      status=nf_inq_varid(ncid,'sbs',varid)
      call handle_error_netcdf('nf_inq_varid: sbs',status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dns)
      call handle_error_netcdf('nf_get_vara_real: sbs',status,nf_noerr)
      call scatter_bnd_ns(sbs0, glb3dns,kb)
      
      status=nf_inq_varid(ncid,'ubs',varid)
      call handle_error_netcdf('nf_inq_varid: ubs',status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dns)
      call handle_error_netcdf('nf_get_vara_real: ubs',status,nf_noerr)
      call scatter_bnd_ns(ubs0, glb3dns,kb)
      
      status=nf_inq_varid(ncid,'vbs',varid)
      call handle_error_netcdf('nf_inq_varid: vbs',status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dns)
      call handle_error_netcdf('nf_get_vara_real: vbs',status,nf_noerr)
      call scatter_bnd_ns(vbs0, glb3dns,kb)
         
!     next month

      month=month+1
      if(month .gt. 12) month=1
      if(my_task .eq. master_task) write(6,*) 'lbcmclim month ', month

!     specify 2d array shape at wester/eastern boudary

      start(1)=1
      count(1)=jm_global
      start(2)=month
      count(2)=1

!     eastern boudary
      status=nf_inq_varid(ncid,'ele',varid)
      call handle_error_netcdf('nf_inq_varid: ele',status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2dew)
      call handle_error_netcdf('nf_get_vara_real: ele',status,nf_noerr)
      call scatter_bnd_ew(ele1,glb2dew,1)
      
      status=nf_inq_varid(ncid,'uabe',varid)
      call handle_error_netcdf('nf_inq_varid: uabe',status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2dew)
      call handle_error_netcdf('nf_get_vara_real: uabe',status,nf_noerr)
      call scatter_bnd_ew(uabe1,glb2dew,1)
      
      status=nf_inq_varid(ncid,'vabe',varid)
      call handle_error_netcdf('nf_inq_varid: vabe',status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2dew)
      call handle_error_netcdf('nf_get_vara_real: vabe',status,nf_noerr)
      call scatter_bnd_ew(vabe1,glb2dew,1)
      
!     western boudary
      status=nf_inq_varid(ncid,'elw',varid)
      call handle_error_netcdf('nf_inq_varid: elw',status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2dew)
      call handle_error_netcdf('nf_get_vara_real: elw',status,nf_noerr)
      call scatter_bnd_ew(elw1,glb2dew,1)
      
      status=nf_inq_varid(ncid,'uabw',varid)
      call handle_error_netcdf('nf_inq_varid: uabw',status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2dew)
      call handle_error_netcdf('nf_get_vara_real: uabw',status,nf_noerr)
      call scatter_bnd_ew(uabw1,glb2dew,1)
      
      status=nf_inq_varid(ncid,'vabw',varid)
      call handle_error_netcdf('nf_inq_varid: vabw',status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2dew)
      call handle_error_netcdf('nf_get_vara_real: vabw',status,nf_noerr)
      call scatter_bnd_ew(vabw1,glb2dew,1)
      
!     specify 3d array shape at wester/eastern boudary
      
      start(1)=1
      count(1)=jm_global
      start(2)=1
      count(2)=kb
      start(3)=month
      count(3)=1
      
!     eastern boudary
      status=nf_inq_varid(ncid,'tbe',varid)
      call handle_error_netcdf('nf_inq_varid: tbe',status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dew)
      call handle_error_netcdf('nf_get_vara_real: tbe',status,nf_noerr)
      call scatter_bnd_ew(tbe1,glb3dew,kb)
      
      status=nf_inq_varid(ncid,'sbe',varid)
      call handle_error_netcdf('nf_inq_varid: sbe',status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dew)
      call handle_error_netcdf('nf_get_vara_real: sbe',status,nf_noerr)
      call scatter_bnd_ew(sbe1,glb3dew,kb)
      
      status=nf_inq_varid(ncid,'ube',varid)
      call handle_error_netcdf('nf_inq_varid: ube',status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dew)
      call handle_error_netcdf('nf_get_vara_real: ube',status,nf_noerr)
      call scatter_bnd_ew(ube1,glb3dew,kb)
      
      status=nf_inq_varid(ncid,'vbe',varid)
      call handle_error_netcdf('nf_inq_varid: vbe',status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dew)
      call handle_error_netcdf('nf_get_vara_real: vbe',status,nf_noerr)
      call scatter_bnd_ew(vbe1,glb3dew,kb)
      
!     western boudary
      status=nf_inq_varid(ncid,'tbw',varid)
      call handle_error_netcdf('nf_inq_varid: tbw',status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dew)
      call handle_error_netcdf('nf_get_vara_real: tbw',status,nf_noerr)
      call scatter_bnd_ew(tbw1,glb3dew,kb)
      
      status=nf_inq_varid(ncid,'sbw',varid)
      call handle_error_netcdf('nf_inq_varid: sbw',status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dew)
      call handle_error_netcdf('nf_get_vara_real: sbw',status,nf_noerr)
      call scatter_bnd_ew(sbw1,glb3dew,kb)
      
      status=nf_inq_varid(ncid,'ubw',varid)
      call handle_error_netcdf('nf_inq_varid: ubw',status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dew)
      call handle_error_netcdf('nf_get_vara_real: ubw',status,nf_noerr)
      call scatter_bnd_ew(ubw1,glb3dew,kb)
      
      status=nf_inq_varid(ncid,'vbw',varid)
      call handle_error_netcdf('nf_inq_varid: vbw',status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dew)
      call handle_error_netcdf('nf_get_vara_real: vbw',status,nf_noerr)
      call scatter_bnd_ew(vbw1,glb3dew,kb)
      
!     specify 2d array shape at northern/southern boudary
      
      start(1)=1
      count(1)=im_global
      start(2)=month
      count(2)=1
      
!     northern boudary
      status=nf_inq_varid(ncid,'eln',varid)
      call handle_error_netcdf('nf_inq_varid: eln',status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2dns)
      call handle_error_netcdf('nf_get_var_real: eln',status,nf_noerr)
      call scatter_bnd_ns(eln1,glb2dns,1)
      
      status=nf_inq_varid(ncid,'uabn',varid)
      call handle_error_netcdf('nf_inq_varid: uabn',status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2dns)
      call handle_error_netcdf('nf_get_var_real: uabn',status,nf_noerr)
      call scatter_bnd_ns(uabn1,glb2dns,1)
      
      status=nf_inq_varid(ncid,'vabn',varid)
      call handle_error_netcdf('nf_inq_varid: vabn',status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2dns)
      call handle_error_netcdf('nf_get_var_real: vabn',status,nf_noerr)
      call scatter_bnd_ns(vabn1,glb2dns,1)
      
!     southern boudary
      status=nf_inq_varid(ncid,'els',varid)
      call handle_error_netcdf('nf_inq_varid: els',status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2dns)
      call handle_error_netcdf('nf_get_var_real: els',status,nf_noerr)
      call scatter_bnd_ns(els1,glb2dns,1)
      
      status=nf_inq_varid(ncid,'uabs',varid)
      call handle_error_netcdf('nf_inq_varid: uabs',status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2dns)
      call handle_error_netcdf('nf_get_var_real: uabs',status,nf_noerr)
      call scatter_bnd_ns(uabs1,glb2dns,1)
      
      status=nf_inq_varid(ncid,'vabs',varid)
      call handle_error_netcdf('nf_inq_varid: vabs',status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2dns)
      call handle_error_netcdf('nf_get_var_real: vabs',status,nf_noerr)
      call scatter_bnd_ns(vabs1,glb2dns,1)
      
!     specify 3d array shape at northern/southern boudary
      
      start(1)=1
      count(1)=im_global
      start(2)=1
      count(2)=kb
      start(3)=month
      count(3)=1
      
!     northern boudary
      status=nf_inq_varid(ncid,'tbn',varid)
      call handle_error_netcdf('nf_inq_varid: tbn',status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dns)
      call handle_error_netcdf('nf_get_vara_real: tbn',status,nf_noerr)
      call scatter_bnd_ns(tbn1, glb3dns,kb)
      
      status=nf_inq_varid(ncid,'sbn',varid)
      call handle_error_netcdf('nf_inq_varid: sbn',status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dns)
      call handle_error_netcdf('nf_get_vara_real: sbn',status,nf_noerr)
      call scatter_bnd_ns(sbn1, glb3dns,kb)
      
      status=nf_inq_varid(ncid,'ubn',varid)
      call handle_error_netcdf('nf_inq_varid: ubn',status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dns)
      call handle_error_netcdf('nf_get_vara_real: ubn',status,nf_noerr)
      call scatter_bnd_ns(ubn1, glb3dns,kb)
      
      status=nf_inq_varid(ncid,'vbn',varid)
      call handle_error_netcdf('nf_inq_varid: vbn',status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dns)
      call handle_error_netcdf('nf_get_vara_real: vbn',status,nf_noerr)
      call scatter_bnd_ns(vbn1, glb3dns,kb)
      
!     southern boudary
      status=nf_inq_varid(ncid,'tbs',varid)
      call handle_error_netcdf('nf_inq_varid: tbs',status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dns)
      call handle_error_netcdf('nf_get_vara_real: tbs',status,nf_noerr)
      call scatter_bnd_ns(tbs1, glb3dns,kb)
      
      status=nf_inq_varid(ncid,'sbs',varid)
      call handle_error_netcdf('nf_inq_varid: sbs',status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dns)
      call handle_error_netcdf('nf_get_vara_real: sbs',status,nf_noerr)
      call scatter_bnd_ns(sbs1, glb3dns,kb)
      
      status=nf_inq_varid(ncid,'ubs',varid)
      call handle_error_netcdf('nf_inq_varid: ubs',status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dns)
      call handle_error_netcdf('nf_get_vara_real: ubs',status,nf_noerr)
      call scatter_bnd_ns(ubs1, glb3dns,kb)
      
      status=nf_inq_varid(ncid,'vbs',varid)
      call handle_error_netcdf('nf_inq_varid: vbs',status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dns)
      call handle_error_netcdf('nf_get_vara_real: vbs',status,nf_noerr)
      call scatter_bnd_ns(vbs1, glb3dns,kb)

!     close file
      status=nf_close(ncid)
      call handle_error_netcdf('nf_close',status,nf_noerr)

      return
      end

!_______________________________________________________________________
      subroutine read_lbc_netcdf

! read lateral boundary condition
      implicit none
# include "mpif.h"
# include "netcdf.inc"
      include 'pom.h'
      character*120 netcdf_lbc_file
      integer ncid,status,varid
      integer lbctime_dimid
      character*30 lbctime_name
      character*4 cyear
      character*2 cmonth,cday
      integer n, start(3),count(3)
      real lbctime,lbctime0,time0_lbc
      real glb2dew(jm_global),glb2dns(im_global)
      real glb3dew(jm_global,kb),glb3dns(im_global,kb)
      integer year,month,day,julday
      integer init
      save init
      data init/0/
! debug
      character*1 cha_my_task
      integer i,j
   
! open netcdf lbc file
      write(netcdf_lbc_file,'(''in/'',a,''.lbc.nc'')') trim(netcdf_file)
      if(my_task.eq.master_task) write(*,'(/''reading file '',a)')
     $  trim(netcdf_lbc_file)
      status=nf_open(netcdf_lbc_file,nf_nowrite,ncid)
      call handle_error_netcdf('nf_open: '//netcdf_lbc_file,
     $                         status,nf_noerr)

! get time
      if(init.eq.0) then
 
        status=nf_inq_dimid(ncid,'t',lbctime_dimid)
        call handle_error_netcdf('nf_inq_dimid: t'
     $                          ,status,nf_noerr)
        status=nf_inq_dimlen(ncid,lbctime_dimid,lbctime_len)

        status=nf_inq_varid(ncid,'lbctime',varid)
        call handle_error_netcdf('nf_inq_varid: lbctime'
     $                          ,status,nf_noerr)

        status=nf_get_att_text(ncid,varid,'units',lbctime_name)
        call handle_error_netcdf('nf_inq_att_text: lbctime'
     $                          ,status,nf_noerr)

        do n=1,2
          start(1)=n
          count(1)=1
          status=nf_get_vara_real(ncid,varid
     $                           ,start,count,lbctime)
          call handle_error_netcdf('nf_get_vara_real: lbctime'
     $                          ,status,nf_noerr)
          if(n.eq.2) lbctime_dayint=lbctime-lbctime0
          lbctime0=lbctime
        end do
      
        cyear=lbctime_name(12:15)
        read(cyear,*) year
        cmonth=lbctime_name(17:18)
        read(cmonth,*) month
        cday=lbctime_name(20:21)
        read(cday,*) day
        lbctime_julday = int(julday(year,month,day)+0.001)
 
        call caldat(int(lbctime_julday+(lbctime_len-1)*lbctime_dayint)
     $           ,year,month,day)

        if(my_task.eq.master_task) then
          write(6,*) 'lbc tlength ',lbctime_len
          write(6,*) 'lbc dayint',lbctime_dayint
          write(6,*) 'lbc start_time ', lbctime_name
          write(6,*) 'lbc end_time ',year,month,day 
        end if

        init=1

! close file
        status=nf_close(ncid)
        call handle_error_netcdf('nf_close',status,nf_noerr)

        return

      end if

      time0_lbc=lbctime_julday-julday_start
      if(my_task.eq.master_task) 
     $ write(6,*) 'lbctime_julday ',lbctime_julday
      if(my_task.eq.master_task) 
     $ write(6,*) 'julday_start ',julday_start

! get number on step-0

      n=(time-dti/86400.-time0_lbc)/lbctime_dayint+1
      if(my_task.eq.master_task) write(6,*) 'time ',time-dti/86400.
      if(my_task.eq.master_task) write(6,*) 'time0_lbc ',time0_lbc
      if(my_task.eq.master_task) 
     $ write(6,*) 'lbctime_dayint ',lbctime_dayint
      if(my_task.eq.master_task) write(6,*) 'lbc step ',n

      if(n.le.0) then
        if(my_task.eq.master_task) write(*,'(/a)')
     $                                     'wrong time of boundary data'
        call finalize_mpi
        stop
      end if

! get variables on step-0
      
      start(1)=1
      count(1)=jm_global
      start(2)=n
      count(2)=1

      status=nf_inq_varid(ncid,'ele',varid)
      call handle_error_netcdf('nf_inq_varid: ele'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2dew)
      call handle_error_netcdf('nf_get_var_real: ele',
     $                          status,nf_noerr)
      call scatter_bnd_ew(ele0,glb2dew,1)

      status=nf_inq_varid(ncid,'uabe',varid)
      call handle_error_netcdf('nf_inq_varid: uabe'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2dew)
      call handle_error_netcdf('nf_get_var_real: uabe',
     $                          status,nf_noerr)
      call scatter_bnd_ew(uabe0,glb2dew,1)

      status=nf_inq_varid(ncid,'vabe',varid)
      call handle_error_netcdf('nf_inq_varid: vabe'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2dew)
      call handle_error_netcdf('nf_get_var_real: vabe',
     $                          status,nf_noerr)
      call scatter_bnd_ew(vabe0,glb2dew,1)

      status=nf_inq_varid(ncid,'elw',varid)
      call handle_error_netcdf('nf_inq_varid: elw'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2dew)
      call handle_error_netcdf('nf_get_var_real: elw',
     $                          status,nf_noerr)
      call scatter_bnd_ew(elw0,glb2dew,1)

      status=nf_inq_varid(ncid,'uabw',varid)
      call handle_error_netcdf('nf_inq_varid: uabw'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2dew)
      call handle_error_netcdf('nf_get_var_real: uabw',
     $                          status,nf_noerr)
      call scatter_bnd_ew(uabw0,glb2dew,1)

      status=nf_inq_varid(ncid,'vabw',varid)
      call handle_error_netcdf('nf_inq_varid: vabw'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2dew)
      call handle_error_netcdf('nf_get_var_real: vabw',
     $                          status,nf_noerr)
      call scatter_bnd_ew(vabw0,glb2dew,1)

      start(1)=1
      count(1)=jm_global
      start(2)=1
      count(2)=kb
      start(3)=n
      count(3)=1

      status=nf_inq_varid(ncid,'tbe',varid)
      call handle_error_netcdf('nf_inq_varid: tbe'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dew)
      call handle_error_netcdf('nf_get_var_real: tbe',
     $                          status,nf_noerr)
      call scatter_bnd_ew(tbe0,glb3dew,kb)

      status=nf_inq_varid(ncid,'sbe',varid)
      call handle_error_netcdf('nf_inq_varid: sbe'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dew)
      call handle_error_netcdf('nf_get_var_real: sbe',
     $                          status,nf_noerr)
      call scatter_bnd_ew(sbe0,glb3dew,kb)

      status=nf_inq_varid(ncid,'ube',varid)
      call handle_error_netcdf('nf_inq_varid: ube'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dew)
      call handle_error_netcdf('nf_get_var_real: ube',
     $                          status,nf_noerr)
      call scatter_bnd_ew(ube0,glb3dew,kb)

      status=nf_inq_varid(ncid,'vbe',varid)
      call handle_error_netcdf('nf_inq_varid: vbe'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dew)
      call handle_error_netcdf('nf_get_var_real: vbe',
     $                          status,nf_noerr)
      call scatter_bnd_ew(vbe0,glb3dew,kb)

      status=nf_inq_varid(ncid,'tbw',varid)
      call handle_error_netcdf('nf_inq_varid: tbw'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dew)
      call handle_error_netcdf('nf_get_var_real: tbw',
     $                          status,nf_noerr)
      call scatter_bnd_ew(tbw0,glb3dew,kb)

      status=nf_inq_varid(ncid,'sbw',varid)
      call handle_error_netcdf('nf_inq_varid: sbw'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dew)
      call handle_error_netcdf('nf_get_var_real: sbw',
     $                          status,nf_noerr)
      call scatter_bnd_ew(sbw0,glb3dew,kb)

      status=nf_inq_varid(ncid,'ubw',varid)
      call handle_error_netcdf('nf_inq_varid: ubw'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dew)
      call handle_error_netcdf('nf_get_var_real: ubw',
     $                          status,nf_noerr)
      call scatter_bnd_ew(ubw0,glb3dew,kb)

      status=nf_inq_varid(ncid,'vbw',varid)
      call handle_error_netcdf('nf_inq_varid: vbw'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dew)
      call handle_error_netcdf('nf_get_var_real: vbw',
     $                          status,nf_noerr)
      call scatter_bnd_ew(vbw0,glb3dew,kb)

      start(1)=1
      count(1)=im_global
      start(2)=n
      count(2)=1

      status=nf_inq_varid(ncid,'eln',varid)
      call handle_error_netcdf('nf_inq_varid: eln'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2dns)
      call handle_error_netcdf('nf_get_var_real: eln',
     $                          status,nf_noerr)
      call scatter_bnd_ns(eln0,glb2dns,1)

      status=nf_inq_varid(ncid,'uabn',varid)
      call handle_error_netcdf('nf_inq_varid: uabn'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2dns)
      call handle_error_netcdf('nf_get_var_real: uabn',
     $                          status,nf_noerr)
      call scatter_bnd_ns(uabn0,glb2dns,1)

      status=nf_inq_varid(ncid,'vabn',varid)
      call handle_error_netcdf('nf_inq_varid: vabn'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2dns)
      call handle_error_netcdf('nf_get_var_real: vabn',
     $                          status,nf_noerr)
      call scatter_bnd_ns(vabn0,glb2dns,1)

      status=nf_inq_varid(ncid,'els',varid)
      call handle_error_netcdf('nf_inq_varid: els'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2dns)
      call handle_error_netcdf('nf_get_var_real: els',
     $                          status,nf_noerr)
      call scatter_bnd_ns(els0,glb2dns,1)

      status=nf_inq_varid(ncid,'uabs',varid)
      call handle_error_netcdf('nf_inq_varid: uabs'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2dns)
      call handle_error_netcdf('nf_get_var_real: uabs',
     $                          status,nf_noerr)
      call scatter_bnd_ns(uabs0,glb2dns,1)

      status=nf_inq_varid(ncid,'vabs',varid)
      call handle_error_netcdf('nf_inq_varid: vabs'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2dns)
      call handle_error_netcdf('nf_get_var_real: vabs',
     $                          status,nf_noerr)
      call scatter_bnd_ns(vabs0,glb2dns,1)

      start(1)=1
      count(1)=im_global
      start(2)=1
      count(2)=kb
      start(3)=n
      count(3)=1

      status=nf_inq_varid(ncid,'tbn',varid)
      call handle_error_netcdf('nf_inq_varid: tbn'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dns)
      call handle_error_netcdf('nf_get_var_real: tbn',
     $                          status,nf_noerr)
      call scatter_bnd_ns(tbn0,glb3dns,kb)

      status=nf_inq_varid(ncid,'sbn',varid)
      call handle_error_netcdf('nf_inq_varid: sbn'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dns)
      call handle_error_netcdf('nf_get_var_real: sbn',
     $                          status,nf_noerr)
      call scatter_bnd_ns(sbn0,glb3dns,kb)

      status=nf_inq_varid(ncid,'ubn',varid)
      call handle_error_netcdf('nf_inq_varid: ubn'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dns)
      call handle_error_netcdf('nf_get_var_real: ubn',
     $                          status,nf_noerr)
      call scatter_bnd_ns(ubn0,glb3dns,kb)

      status=nf_inq_varid(ncid,'vbn',varid)
      call handle_error_netcdf('nf_inq_varid: vbn'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dns)
      call handle_error_netcdf('nf_get_var_real: vbn',
     $                          status,nf_noerr)
      call scatter_bnd_ns(vbn0,glb3dns,kb)

      status=nf_inq_varid(ncid,'tbs',varid)
      call handle_error_netcdf('nf_inq_varid: tbs'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dns)
      call handle_error_netcdf('nf_get_var_real: tbs',
     $                          status,nf_noerr)
      call scatter_bnd_ns(tbs0,glb3dns,kb)

      status=nf_inq_varid(ncid,'sbs',varid)
      call handle_error_netcdf('nf_inq_varid: sbs'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dns)
      call handle_error_netcdf('nf_get_var_real: sbs',
     $                          status,nf_noerr)
      call scatter_bnd_ns(sbs0,glb3dns,kb)

      status=nf_inq_varid(ncid,'ubs',varid)
      call handle_error_netcdf('nf_inq_varid: ubs'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dns)
      call handle_error_netcdf('nf_get_var_real: ubs',
     $                          status,nf_noerr)
      call scatter_bnd_ns(ubs0,glb3dns,kb)

      status=nf_inq_varid(ncid,'vbs',varid)
      call handle_error_netcdf('nf_inq_varid: vbs'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dns)
      call handle_error_netcdf('nf_get_var_real: vbs',
     $                          status,nf_noerr)
      call scatter_bnd_ns(vbs0,glb3dns,kb)

! get number on step-1

      n=(time-dti/86400.+lbctime_dayint-time0_lbc)/lbctime_dayint+1
      if(my_task.eq.master_task) write(6,*) 'lbc step ',n

      if(n.le.0) then
        if(my_task.eq.master_task) write(*,'(/a)')
     $                                     'wrong time of boundary data'
        call finalize_mpi
        stop
      end if

! get variables on step-1
      
      start(1)=1
      count(1)=jm_global
      start(2)=n
      count(2)=1

      status=nf_inq_varid(ncid,'ele',varid)
      call handle_error_netcdf('nf_inq_varid: ele'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2dew)
      call handle_error_netcdf('nf_get_var_real: ele',
     $                          status,nf_noerr)
      call scatter_bnd_ew(ele1,glb2dew,1)

      status=nf_inq_varid(ncid,'uabe',varid)
      call handle_error_netcdf('nf_inq_varid: uabe'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2dew)
      call handle_error_netcdf('nf_get_var_real: uabe',
     $                          status,nf_noerr)
      call scatter_bnd_ew(uabe1,glb2dew,1)

      status=nf_inq_varid(ncid,'vabe',varid)
      call handle_error_netcdf('nf_inq_varid: vabe'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2dew)
      call handle_error_netcdf('nf_get_var_real: vabe',
     $                          status,nf_noerr)
      call scatter_bnd_ew(vabe1,glb2dew,1)

      status=nf_inq_varid(ncid,'elw',varid)
      call handle_error_netcdf('nf_inq_varid: elw'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2dew)
      call handle_error_netcdf('nf_get_var_real: elw',
     $                          status,nf_noerr)
      call scatter_bnd_ew(elw1,glb2dew,1)

      status=nf_inq_varid(ncid,'uabw',varid)
      call handle_error_netcdf('nf_inq_varid: uabw'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2dew)
      call handle_error_netcdf('nf_get_var_real: uabw',
     $                          status,nf_noerr)
      call scatter_bnd_ew(uabw1,glb2dew,1)

      status=nf_inq_varid(ncid,'vabw',varid)
      call handle_error_netcdf('nf_inq_varid: vabw'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2dew)
      call handle_error_netcdf('nf_get_var_real: vabw',
     $                          status,nf_noerr)
      call scatter_bnd_ew(vabw1,glb2dew,1)

      start(1)=1
      count(1)=jm_global
      start(2)=1
      count(2)=kb
      start(3)=n
      count(3)=1

      status=nf_inq_varid(ncid,'tbe',varid)
      call handle_error_netcdf('nf_inq_varid: tbe'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dew)
      call handle_error_netcdf('nf_get_var_real: tbe',
     $                          status,nf_noerr)
      call scatter_bnd_ew(tbe1,glb3dew,kb)

      status=nf_inq_varid(ncid,'sbe',varid)
      call handle_error_netcdf('nf_inq_varid: sbe'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dew)
      call handle_error_netcdf('nf_get_var_real: sbe',
     $                          status,nf_noerr)
      call scatter_bnd_ew(sbe1,glb3dew,kb)

      status=nf_inq_varid(ncid,'ube',varid)
      call handle_error_netcdf('nf_inq_varid: ube'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dew)
      call handle_error_netcdf('nf_get_var_real: ube',
     $                          status,nf_noerr)
      call scatter_bnd_ew(ube1,glb3dew,kb)

      status=nf_inq_varid(ncid,'vbe',varid)
      call handle_error_netcdf('nf_inq_varid: vbe'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dew)
      call handle_error_netcdf('nf_get_var_real: vbe',
     $                          status,nf_noerr)
      call scatter_bnd_ew(vbe1,glb3dew,kb)

      status=nf_inq_varid(ncid,'tbw',varid)
      call handle_error_netcdf('nf_inq_varid: tbw'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dew)
      call handle_error_netcdf('nf_get_var_real: tbw',
     $                          status,nf_noerr)
      call scatter_bnd_ew(tbw1,glb3dew,kb)

      status=nf_inq_varid(ncid,'sbw',varid)
      call handle_error_netcdf('nf_inq_varid: sbw'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dew)
      call handle_error_netcdf('nf_get_var_real: sbw',
     $                          status,nf_noerr)
      call scatter_bnd_ew(sbw1,glb3dew,kb)

      status=nf_inq_varid(ncid,'ubw',varid)
      call handle_error_netcdf('nf_inq_varid: ubw'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dew)
      call handle_error_netcdf('nf_get_var_real: ubw',
     $                          status,nf_noerr)
      call scatter_bnd_ew(ubw1,glb3dew,kb)

      status=nf_inq_varid(ncid,'vbw',varid)
      call handle_error_netcdf('nf_inq_varid: vbw'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dew)
      call handle_error_netcdf('nf_get_var_real: vbw',
     $                          status,nf_noerr)
      call scatter_bnd_ew(vbw1,glb3dew,kb)

      start(1)=1
      count(1)=im_global
      start(2)=n
      count(2)=1

      status=nf_inq_varid(ncid,'eln',varid)
      call handle_error_netcdf('nf_inq_varid: eln'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2dns)
      call handle_error_netcdf('nf_get_var_real: eln',
     $                          status,nf_noerr)
      call scatter_bnd_ns(eln1,glb2dns,1)

      status=nf_inq_varid(ncid,'uabn',varid)
      call handle_error_netcdf('nf_inq_varid: uabn'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2dns)
      call handle_error_netcdf('nf_get_var_real: uabn',
     $                          status,nf_noerr)
      call scatter_bnd_ns(uabn1,glb2dns,1)

      status=nf_inq_varid(ncid,'vabn',varid)
      call handle_error_netcdf('nf_inq_varid: vabn'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2dns)
      call handle_error_netcdf('nf_get_var_real: vabn',
     $                          status,nf_noerr)
      call scatter_bnd_ns(vabn1,glb2dns,1)

      status=nf_inq_varid(ncid,'els',varid)
      call handle_error_netcdf('nf_inq_varid: els'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2dns)
      call handle_error_netcdf('nf_get_var_real: els',
     $                          status,nf_noerr)
      call scatter_bnd_ns(els1,glb2dns,1)

      status=nf_inq_varid(ncid,'uabs',varid)
      call handle_error_netcdf('nf_inq_varid: uabs'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2dns)
      call handle_error_netcdf('nf_get_var_real: uabs',
     $                          status,nf_noerr)
      call scatter_bnd_ns(uabs1,glb2dns,1)

      status=nf_inq_varid(ncid,'vabs',varid)
      call handle_error_netcdf('nf_inq_varid: vabs'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2dns)
      call handle_error_netcdf('nf_get_var_real: vabs',
     $                          status,nf_noerr)
      call scatter_bnd_ns(vabs1,glb2dns,1)

      start(1)=1
      count(1)=im_global
      start(2)=1
      count(2)=kb
      start(3)=n
      count(3)=1

      status=nf_inq_varid(ncid,'tbn',varid)
      call handle_error_netcdf('nf_inq_varid: tbn'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dns)
      call handle_error_netcdf('nf_get_var_real: tbn',
     $                          status,nf_noerr)
      call scatter_bnd_ns(tbn1,glb3dns,kb)

      status=nf_inq_varid(ncid,'sbn',varid)
      call handle_error_netcdf('nf_inq_varid: sbn'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dns)
      call handle_error_netcdf('nf_get_var_real: sbn',
     $                          status,nf_noerr)
      call scatter_bnd_ns(sbn1,glb3dns,kb)

      status=nf_inq_varid(ncid,'ubn',varid)
      call handle_error_netcdf('nf_inq_varid: ubn'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dns)
      call handle_error_netcdf('nf_get_var_real: ubn',
     $                          status,nf_noerr)
      call scatter_bnd_ns(ubn1,glb3dns,kb)

      status=nf_inq_varid(ncid,'vbn',varid)
      call handle_error_netcdf('nf_inq_varid: vbn'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dns)
      call handle_error_netcdf('nf_get_var_real: vbn',
     $                          status,nf_noerr)
      call scatter_bnd_ns(vbn1,glb3dns,kb)

      status=nf_inq_varid(ncid,'tbs',varid)
      call handle_error_netcdf('nf_inq_varid: tbs'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dns)
      call handle_error_netcdf('nf_get_var_real: tbs',
     $                          status,nf_noerr)
      call scatter_bnd_ns(tbs1,glb3dns,kb)

      status=nf_inq_varid(ncid,'sbs',varid)
      call handle_error_netcdf('nf_inq_varid: sbs'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dns)
      call handle_error_netcdf('nf_get_var_real: sbs',
     $                          status,nf_noerr)
      call scatter_bnd_ns(sbs1,glb3dns,kb)

      status=nf_inq_varid(ncid,'ubs',varid)
      call handle_error_netcdf('nf_inq_varid: ubs'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dns)
      call handle_error_netcdf('nf_get_var_real: ubs',
     $                          status,nf_noerr)
      call scatter_bnd_ns(ubs1,glb3dns,kb)

      status=nf_inq_varid(ncid,'vbs',varid)
      call handle_error_netcdf('nf_inq_varid: vbs'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3dns)
      call handle_error_netcdf('nf_get_var_real: vbs',
     $                          status,nf_noerr)
      call scatter_bnd_ns(vbs1,glb3dns,kb)

! close file
      status=nf_close(ncid)
      call handle_error_netcdf('nf_close',status,nf_noerr)

      return
      end

!_______________________________________________________________________
      subroutine read_atm_netcdf

! read atmospheric condition 
      implicit none
# include "mpif.h"
# include "netcdf.inc"
      include 'pom.h'
      character*120 netcdf_atm_file
      integer ncid,status,varid
      integer atmtime_dimid
      character*30 atmtime_name
      character*4 cyear
      character*2 cmonth,cday
      integer n, start(3),count(3)
      real atmtime,atmtime0,time0_atm
      real glb2d(im_global,jm_global)
      integer year,month,day,julday
      integer init
      save init
      data init/0/
! debug
      character*1 cha_my_task
      integer i,j
   
! open netcdf atm file
      write(netcdf_atm_file,'(''in/'',a,''.atm.nc'')') trim(netcdf_file)
      if(my_task.eq.master_task) write(*,'(/''reading file '',a)')
     $  trim(netcdf_atm_file)
      status=nf_open(netcdf_atm_file,nf_nowrite,ncid)
      call handle_error_netcdf('nf_open: '//netcdf_atm_file,
     $                         status,nf_noerr)

! get time
      if(init.eq.0) then
 
        status=nf_inq_dimid(ncid,'t',atmtime_dimid)
        call handle_error_netcdf('nf_inq_dimid: t'
     $                          ,status,nf_noerr)
        status=nf_inq_dimlen(ncid,atmtime_dimid,atmtime_len)

        status=nf_inq_varid(ncid,'atmtime',varid)
        call handle_error_netcdf('nf_inq_varid: atmtime'
     $                          ,status,nf_noerr)

        status=nf_get_att_text(ncid,varid,'units',atmtime_name)
        call handle_error_netcdf('nf_inq_att_text: atmtime'
     $                          ,status,nf_noerr)

        do n=1,2
          start(1)=n
          count(1)=1
          status=nf_get_vara_real(ncid,varid
     $                           ,start,count,atmtime)
          call handle_error_netcdf('nf_get_vara_real: atmtime'
     $                          ,status,nf_noerr)
          if(n.eq.2) atmtime_dayint=atmtime-atmtime0
          atmtime0=atmtime
        end do
      
        cyear=atmtime_name(12:15)
        read(cyear,*) year
        cmonth=atmtime_name(17:18)
        read(cmonth,*) month
        cday=atmtime_name(20:21)
        read(cday,*) day
        atmtime_julday = int(julday(year,month,day)+0.001)
 
        call caldat(int(atmtime_julday+(atmtime_len-1)*atmtime_dayint)
     $           ,year,month,day)

        if(my_task.eq.master_task) then
          write(6,*) 'atm tlength ',atmtime_len
          write(6,*) 'atm dayint',atmtime_dayint
          write(6,*) 'atm start_time ', atmtime_name
          write(6,*) 'atm end_time ',year,month,day 
        end if

        init=1

! close file
        status=nf_close(ncid)
        call handle_error_netcdf('nf_close',status,nf_noerr)

        return

      end if

      time0_atm=atmtime_julday-julday_start

! get number on step-0

      n=(time-dti/86400.-time0_atm)/atmtime_dayint+1

      if(my_task.eq.master_task) write(6,*) 'atm step ',n

      if(n.le.0) then
        if(my_task.eq.master_task) write(*,'(/a)')
     $                              'wrong time of atmospheric data'
        call finalize_mpi
        stop
      end if

! get variables on step-0
      
      start(1)=1
      count(1)=im_global
      start(2)=1
      count(2)=jm_global
      start(3)=n
      count(3)=1

      status=nf_inq_varid(ncid,'windu',varid)
      call handle_error_netcdf('nf_inq_varid: windu'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2d)
      call handle_error_netcdf('nf_get_var_real: windu',
     $                          status,nf_noerr)
      call scatter2d(windu0,glb2d)

      status=nf_inq_varid(ncid,'windv',varid)
      call handle_error_netcdf('nf_inq_varid: windv'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2d)
      call handle_error_netcdf('nf_get_var_real: windv',
     $                          status,nf_noerr)
      call scatter2d(windv0,glb2d)

      status=nf_inq_varid(ncid,'airt',varid)
      call handle_error_netcdf('nf_inq_varid: airt'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2d)
      call handle_error_netcdf('nf_get_var_real: airt',
     $                          status,nf_noerr)
      call scatter2d(airt0,glb2d)

      status=nf_inq_varid(ncid,'airh',varid)
      call handle_error_netcdf('nf_inq_varid: airh'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2d)
      call handle_error_netcdf('nf_get_var_real: airh',
     $                          status,nf_noerr)
      call scatter2d(airh0,glb2d)

      status=nf_inq_varid(ncid,'swrad',varid)
      call handle_error_netcdf('nf_inq_varid: swrad'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2d)
      call handle_error_netcdf('nf_get_var_real: swrad',
     $                          status,nf_noerr)
      call scatter2d(swrad0,glb2d)

      status=nf_inq_varid(ncid,'cloud',varid)
      call handle_error_netcdf('nf_inq_varid: cloud'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2d)
      call handle_error_netcdf('nf_get_var_real: cloud',
     $                          status,nf_noerr)
      call scatter2d(cloud0,glb2d)

      status=nf_inq_varid(ncid,'slp',varid)
      call handle_error_netcdf('nf_inq_varid: slp'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2d)
      call handle_error_netcdf('nf_get_var_real: slp',
     $                          status,nf_noerr)
      call scatter2d(slp0,glb2d)

! get number on step-1

      n=(time-dti/86400.+atmtime_dayint-time0_atm)/atmtime_dayint+1
      if(my_task.eq.master_task) write(6,*) 'atm step ',n

      if(n.le.0) then
        if(my_task.eq.master_task) write(*,'(/a)')
     $                              'wrong time of atmospheric data'
        call finalize_mpi
        stop
      end if

! get variables on step-1

      start(1)=1
      count(1)=im_global
      start(2)=1
      count(2)=jm_global
      start(3)=n
      count(3)=1

      status=nf_inq_varid(ncid,'windu',varid)
      call handle_error_netcdf('nf_inq_varid: windu'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2d)
      call handle_error_netcdf('nf_get_var_real: windu',
     $                          status,nf_noerr)
      call scatter2d(windu1,glb2d)

      status=nf_inq_varid(ncid,'windv',varid)
      call handle_error_netcdf('nf_inq_varid: windv'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2d)
      call handle_error_netcdf('nf_get_var_real: windv',
     $                          status,nf_noerr)
      call scatter2d(windv1,glb2d)

      status=nf_inq_varid(ncid,'airt',varid)
      call handle_error_netcdf('nf_inq_varid: airt'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2d)
      call handle_error_netcdf('nf_get_var_real: airt',
     $                          status,nf_noerr)
      call scatter2d(airt1,glb2d)

      status=nf_inq_varid(ncid,'airh',varid)
      call handle_error_netcdf('nf_inq_varid: airh'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2d)
      call handle_error_netcdf('nf_get_var_real: airh',
     $                          status,nf_noerr)
      call scatter2d(airh1,glb2d)

      status=nf_inq_varid(ncid,'swrad',varid)
      call handle_error_netcdf('nf_inq_varid: swrad'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2d)
      call handle_error_netcdf('nf_get_var_real: swrad',
     $                          status,nf_noerr)
      call scatter2d(swrad1,glb2d)

      status=nf_inq_varid(ncid,'cloud',varid)
      call handle_error_netcdf('nf_inq_varid: cloud'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2d)
      call handle_error_netcdf('nf_get_var_real: cloud',
     $                          status,nf_noerr)
      call scatter2d(cloud1,glb2d)

      status=nf_inq_varid(ncid,'slp',varid)
      call handle_error_netcdf('nf_inq_varid: slp'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2d)
      call handle_error_netcdf('nf_get_var_real: slp',
     $                          status,nf_noerr)
      call scatter2d(slp1,glb2d)

! close file
      status=nf_close(ncid)
      call handle_error_netcdf('nf_close',status,nf_noerr)

      return
      end

!_______________________________________________________________________
      subroutine read_fflux_netcdf

! read fresh water flux condition
!     2018.06.06 S.Ohishi
!     2018.08.20 S.Ohishi
      implicit none
# include "mpif.h"
# include "netcdf.inc"
      include 'pom.h'
      character*120 netcdf_fflux_file
      integer ncid,status,varid
      integer ffluxtime_dimid
      character*30 ffluxtime_name
      character*4 cyear
      character*2 cmonth,cday
      integer n, start(3),count(3)
      real ffluxtime,ffluxtime0,time0_fflux
      real glb2d(im_global,jm_global)
      integer year,month,day,julday
      integer init
      save init
      data init/0/
! debug
      character*1 cha_my_task
      integer i,j
   
! open netcdf fresh water flux file
      write(netcdf_fflux_file,'(''in/'',a,''.fflux.nc'')') 
     $     trim(netcdf_file)
      if(my_task.eq.master_task) write(*,'(/''reading file '',a)')
     $  trim(netcdf_fflux_file)
      status=nf_open(netcdf_fflux_file,nf_nowrite,ncid)
      call handle_error_netcdf('nf_open: '//netcdf_fflux_file,
     $                         status,nf_noerr)

! get time
      if(init.eq.0) then
 
        status=nf_inq_dimid(ncid,'t',ffluxtime_dimid)
        call handle_error_netcdf('nf_inq_dimid: t'
     $                          ,status,nf_noerr)
        status=nf_inq_dimlen(ncid,ffluxtime_dimid,ffluxtime_len)

        status=nf_inq_varid(ncid,'ffluxtime',varid)
        call handle_error_netcdf('nf_inq_varid: ffluxtime'
     $                          ,status,nf_noerr)

        status=nf_get_att_text(ncid,varid,'units',ffluxtime_name)
        call handle_error_netcdf('nf_inq_att_text: ffluxtime'
     $                          ,status,nf_noerr)

        do n=1,2
          start(1)=n
          count(1)=1
          status=nf_get_vara_real(ncid,varid
     $                           ,start,count,ffluxtime)
          call handle_error_netcdf('nf_get_vara_real: ffluxtime'
     $                          ,status,nf_noerr)
          if(n.eq.2) ffluxtime_dayint=ffluxtime-ffluxtime0
          ffluxtime0=ffluxtime
        end do
      
        cyear=ffluxtime_name(12:15)
        read(cyear,*) year
        cmonth=ffluxtime_name(17:18)
        read(cmonth,*) month
        cday=ffluxtime_name(20:21)
        read(cday,*) day
        ffluxtime_julday = int(julday(year,month,day)+0.001)
 
        call caldat(
     $       int(ffluxtime_julday+(ffluxtime_len-1)*ffluxtime_dayint)
     $       ,year,month,day)

        if(my_task.eq.master_task) then
          write(6,*) 'fflux tlength ',ffluxtime_len
          write(6,*) 'fflux dayint',ffluxtime_dayint
          write(6,*) 'fflux start_time ', ffluxtime_name
          write(6,*) 'fflux end_time ',year,month,day 
        end if

        init=1

! close file
        status=nf_close(ncid)
        call handle_error_netcdf('nf_close',status,nf_noerr)

        return

      end if

      time0_fflux=ffluxtime_julday-julday_start

! get number on step-0

      n=(time-dti/86400.-time0_fflux)/ffluxtime_dayint+1

      if(my_task.eq.master_task) write(6,*) 'fflux step ',n

      if(n.le.0) then
        if(my_task.eq.master_task) write(*,'(/a)')
     $                              'wrong time of freshwater flux data'
        call finalize_mpi
        stop
      end if

! get variables on step-0
      
      start(1)=1
      count(1)=im_global
      start(2)=1
      count(2)=jm_global
      start(3)=n
      count(3)=1

!     Evap
!      status=nf_inq_varid(ncid,'evap',varid)
!      call handle_error_netcdf('nf_inq_varid: evap'
!     $                          ,status,nf_noerr)
!      status=nf_get_vara_real(ncid,varid,start,count,glb2d)
!      call handle_error_netcdf('nf_get_var_real: evap',
!     $                          status,nf_noerr)
!      call scatter2d(evap0,glb2d)

!     Prep
      status=nf_inq_varid(ncid,'prep',varid)
      call handle_error_netcdf('nf_inq_varid: prep'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2d)
      call handle_error_netcdf('nf_get_var_real: prep',
     $                          status,nf_noerr)
      call scatter2d(prep0,glb2d)

!     River discharge
      status=nf_inq_varid(ncid,'river',varid)
      call handle_error_netcdf('nf_inq_varid: river'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2d)
      call handle_error_netcdf('nf_get_var_real: river',
     $                          status,nf_noerr)
      call scatter2d(river0,glb2d)

!     Fresh water flux (E-P-R)
!      status=nf_inq_varid(ncid,'fflux',varid)
!      call handle_error_netcdf('nf_inq_varid: fflux'
!     $                          ,status,nf_noerr)
!      status=nf_get_vara_real(ncid,varid,start,count,glb2d)
!      call handle_error_netcdf('nf_get_var_real: fflux',
!     $                          status,nf_noerr)
!      call scatter2d(fflux0,glb2d)
      
! get number on step-1

      n=(time-dti/86400.+ffluxtime_dayint-time0_fflux)
     $     /ffluxtime_dayint+1
      if(my_task.eq.master_task) write(6,*) 'fflux step ',n

      if(n.le.0) then
        if(my_task.eq.master_task) write(*,'(/a)')
     $        'wrong time of fresh water flux data'
        call finalize_mpi
        stop
      end if

! get variables on step-1

      start(1)=1
      count(1)=im_global
      start(2)=1
      count(2)=jm_global
      start(3)=n
      count(3)=1

!     Evap
!      status=nf_inq_varid(ncid,'evap',varid)
!      call handle_error_netcdf('nf_inq_varid: evap'
!     $                          ,status,nf_noerr)
!      status=nf_get_vara_real(ncid,varid,start,count,glb2d)
!      call handle_error_netcdf('nf_get_var_real: evap',
!     $                          status,nf_noerr)
!      call scatter2d(evap1,glb2d)

!     Prep
      status=nf_inq_varid(ncid,'prep',varid)
      call handle_error_netcdf('nf_inq_varid: prep'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2d)
      call handle_error_netcdf('nf_get_var_real: prep',
     $                          status,nf_noerr)
      call scatter2d(prep1,glb2d)

!     River discharge
      status=nf_inq_varid(ncid,'river',varid)
      call handle_error_netcdf('nf_inq_varid: river'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb2d)
      call handle_error_netcdf('nf_get_var_real: river',
     $                          status,nf_noerr)
      call scatter2d(river1,glb2d)

!     Fresh water flux (E-P-R)
!      status=nf_inq_varid(ncid,'fflux',varid)
!      call handle_error_netcdf('nf_inq_varid: fflux'
!     $                          ,status,nf_noerr)
!      status=nf_get_vara_real(ncid,varid,start,count,glb2d)
!      call handle_error_netcdf('nf_get_var_real: fflux',
!     $                          status,nf_noerr)
!      call scatter2d(fflux1,glb2d)

! close file
      status=nf_close(ncid)
      call handle_error_netcdf('nf_close',status,nf_noerr)

      return
      end

!_______________________________________________________________________
      subroutine read_tsclim_netcdf(month)

! read temperature and salinity climatological/mean data
      implicit none
# include "mpif.h"
# include "netcdf.inc"
      include 'pom.h'
      integer month
      character*120 netcdf_tsclim_file
      integer ncid,status,varid
      integer start(3),count(3)
      real glb3d(im_global,jm_global,kb)

! open netcdf file
      write(netcdf_tsclim_file,'(''in/'',a,''.tsclim.nc'')') 
     $                         trim(netcdf_file)
      if(my_task.eq.master_task) write(*,'(/''reading file '',a)')
     $  trim(netcdf_tsclim_file)
      status=nf_open(netcdf_tsclim_file,nf_nowrite,ncid)
      call handle_error_netcdf('nf_open: '//netcdf_tsclim_file,
     $                         status,nf_noerr)

! specify array shape
      start(1)=1
      count(1)=im_global
      start(2)=1
      count(2)=jm_global
      start(3)=1
      count(3)=kb

! read data
      status=nf_inq_varid(ncid,'tclim',varid)
      call handle_error_netcdf('nf_inq_varid: tclim'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3d)
      call handle_error_netcdf('nf_get_var_real: tclim',
     $                          status,nf_noerr)
      call scatter3d(tclim,glb3d)      

      status=nf_inq_varid(ncid,'sclim',varid)
      call handle_error_netcdf('nf_inq_varid: sclim'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3d)
      call handle_error_netcdf('nf_get_var_real: sclim',
     $                          status,nf_noerr)
      call scatter3d(sclim,glb3d)      

!      status=nf_inq_varid(ncid,'tmean',varid)
!      call handle_error_netcdf('nf_inq_varid: tmean'
!     $                          ,status,nf_noerr)
!      status=nf_get_vara_real(ncid,varid,start,count,glb3d)
!      call handle_error_netcdf('nf_get_var_real: tmean',
!     $                          status,nf_noerr)
!      call scatter3d(tmean,glb3d)      

!      status=nf_inq_varid(ncid,'smean',varid)
!      call handle_error_netcdf('nf_inq_varid: smean'
!     $                          ,status,nf_noerr)
!      status=nf_get_vara_real(ncid,varid,start,count,glb3d)
!      call handle_error_netcdf('nf_get_var_real: smean',
!     $                          status,nf_noerr)
!      call scatter3d(smean,glb3d)      

! close file
      status=nf_close(ncid)
      call handle_error_netcdf('nf_close',status,nf_noerr)

      return
      end

!_______________________________________________________________________
      subroutine read_tsclim_monthly_netcdf(month)

! read temperature and salinity climatological data
      implicit none
# include "mpif.h"
# include "netcdf.inc"
      include 'pom.h'
      integer month
      character*120 netcdf_tsclim_file
      integer ncid,status,varid
      integer n,start(4),count(4)
      real glb3d(im_global,jm_global,kb)

! open netcdf atm file
      write(netcdf_tsclim_file,'(''in/'',a,''.tsclim.nc'')') 
     $                         trim(netcdf_file)
      if(my_task.eq.master_task) write(*,'(/''reading file '',a)')
     $  trim(netcdf_tsclim_file)
      status=nf_open(netcdf_tsclim_file,nf_nowrite,ncid)
      call handle_error_netcdf('nf_open: '//netcdf_tsclim_file,
     $                         status,nf_noerr)

! check month

      if(my_task.eq.master_task) write(6,*) 'tsclim month ',month
      if(month.le.0.or.month.gt.12) then
        if(my_task.eq.master_task) write(*,'(/a)')
     $                            'wrong month in read_tsclim_netcdf'
        call finalize_mpi
        stop
      end if

! get variables on step-0

      start(1)=1
      count(1)=im_global
      start(2)=1
      count(2)=jm_global
      start(3)=1
      count(3)=kb
      start(4)=month
      count(4)=1

      status=nf_inq_varid(ncid,'tclim',varid)
      call handle_error_netcdf('nf_inq_varid: tclim'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3d)
      call handle_error_netcdf('nf_get_var_real: tclim',
     $                          status,nf_noerr)
      call scatter3d(tclim0,glb3d)      

      status=nf_inq_varid(ncid,'sclim',varid)
      call handle_error_netcdf('nf_inq_varid: sclim'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3d)
      call handle_error_netcdf('nf_get_var_real: sclim',
     $                          status,nf_noerr)
      call scatter3d(sclim0,glb3d)      

! get month on step-1

      month=month+1
      if(month.gt.12) month=1
      if(my_task.eq.master_task) write(6,*) 'tsclim month ',month

! get variables on step-1

      start(1)=1
      count(1)=im_global
      start(2)=1
      count(2)=jm_global
      start(3)=1
      count(3)=kb
      start(4)=month
      count(4)=1

      status=nf_inq_varid(ncid,'tclim',varid)
      call handle_error_netcdf('nf_inq_varid: tclim'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3d)
      call handle_error_netcdf('nf_get_var_real: tclim',
     $                          status,nf_noerr)
      call scatter3d(tclim1,glb3d)      

      status=nf_inq_varid(ncid,'sclim',varid)
      call handle_error_netcdf('nf_inq_varid: sclim'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3d)
      call handle_error_netcdf('nf_get_var_real: sclim',
     $                          status,nf_noerr)
      call scatter3d(sclim1,glb3d)      

! close file
      status=nf_close(ncid)
      call handle_error_netcdf('nf_close',status,nf_noerr)

      return
      end

!_______________________________________________________________________
      subroutine read_restart_tsuv_netcdf
! read data for a seamless restart
      implicit none
      return
      end

!_______________________________________________________________________

      subroutine read_tsdata_mclim_netcdf(month)

!     read temperature salinity data
      implicit none
# include "mpif.h"
# include "netcdf.inc"
      include 'pom.h'

      integer i,j
      integer ncid,status,varid
      integer month
      integer start(4),count(4)

      real glb3d(im_global,jm_global,kb)

      character*120 netcdf_tsdata_file

!     open ntecdf file
      write(netcdf_tsdata_file,'(''in/'',a,''.tsdata.nc'')') 
     $     trim(netcdf_file)
      if(my_task .eq. master_task)then
         write(*,'(''reading file '',a)') trim(netcdf_tsdata_file)
      end if
      status=nf_open(netcdf_tsdata_file,nf_nowrite,ncid)
      call handle_error_netcdf('nf_open: '//netcdf_tsdata_file,
     $     status,nf_noerr)

!     present month
      
      start(1)=1
      count(1)=im_global
      start(2)=1
      count(2)=jm_global
      start(3)=1
      count(3)=kb
      start(4)=month
      count(4)=1

      status=nf_inq_varid(ncid,'temp',varid)
      call handle_error_netcdf('nf_inq_varid: temp',status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3d)
      call handle_error_netcdf('nf_get_vara_real: temp',status,nf_noerr)
      call scatter3d(tref0,glb3d)

      status=nf_inq_varid(ncid,'sal',varid)
      call handle_error_netcdf('nf_inq_varid: sal',status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3d)
      call handle_error_netcdf('nf_get_vara_real: sal',status,nf_noerr)
      call scatter3d(sref0,glb3d)

!     next month

      month=month+1
      if(month .gt. 12) month=1
      if(my_task .eq. master_task) write(6,*) 'tsdata month ',month

      start(1)=1
      count(1)=im_global
      start(2)=1
      count(2)=jm_global
      start(3)=1
      count(3)=kb
      start(4)=month
      count(4)=1

      status=nf_inq_varid(ncid,'temp',varid)
      call handle_error_netcdf('nf_inq_varid: temp',status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3d)
      call handle_error_netcdf('nf_get_vara_real: temp',status,nf_noerr)
      call scatter3d(tref1,glb3d)

      status=nf_inq_varid(ncid,'sal',varid)
      call handle_error_netcdf('nf_inq_varid: sal',status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3d)
      call handle_error_netcdf('nf_get_vara_real: sal',status,nf_noerr)
      call scatter3d(sref1,glb3d)

!     close file
      status=nf_close(ncid)
      call handle_error_netcdf('nf_close',status,nf_noerr)
      
      return
      end

!_______________________________________________________________________
      subroutine read_tsdata_netcdf

! read temperature salinity data
      implicit none
# include "mpif.h"
# include "netcdf.inc"
      include 'pom.h'
      character*120 netcdf_tsdata_file
      integer ncid,status,varid
      integer tsdatatime_dimid
      character*30 tsdatatime_name
      character*4 cyear
      character*2 cmonth,cday
      integer n, start(4),count(4)
      real tsdatatime,tsdatatime0,time0_tsdata
      real glb3d(im_global,jm_global,kb)
      integer year,month,day,julday
      integer init
      save init
      data init/0/
! debug
      character*1 cha_my_task
      integer i,j
   
! open netcdf tsdata file
      write(netcdf_tsdata_file,'(''in/'',a,''.tsdata.nc'')') 
     $ trim(netcdf_file)
      if(my_task.eq.master_task) write(*,'(/''reading file '',a)')
     $  trim(netcdf_tsdata_file)
      status=nf_open(netcdf_tsdata_file,nf_nowrite,ncid)
      call handle_error_netcdf('nf_open: '//netcdf_tsdata_file,
     $                         status,nf_noerr)

! get time
      if(init.eq.0) then
 
        status=nf_inq_dimid(ncid,'t',tsdatatime_dimid)
        call handle_error_netcdf('nf_inq_dimid: t'
     $                          ,status,nf_noerr)
        status=nf_inq_dimlen(ncid,tsdatatime_dimid,tsdatatime_len)

        status=nf_inq_varid(ncid,'tsdatatime',varid)
        call handle_error_netcdf('nf_inq_varid: tsdatatime'
     $                          ,status,nf_noerr)

        status=nf_get_att_text(ncid,varid,'units',tsdatatime_name)
        call handle_error_netcdf('nf_inq_att_text: tsdatatime'
     $                          ,status,nf_noerr)

        do n=1,2
          start(1)=n
          count(1)=1
          status=nf_get_vara_real(ncid,varid
     $                           ,start,count,tsdatatime)
          call handle_error_netcdf('nf_get_vara_real: tsdatatime'
     $                          ,status,nf_noerr)
          if(n.eq.2) tsdatatime_dayint=tsdatatime-tsdatatime0
          tsdatatime0=tsdatatime
        end do
      
        cyear=tsdatatime_name(12:15)
        read(cyear,*) year
        cmonth=tsdatatime_name(17:18)
        read(cmonth,*) month
        cday=tsdatatime_name(20:21)
        read(cday,*) day
        tsdatatime_julday = int(julday(year,month,day)+0.001)
 
        call caldat(
     $     int(tsdatatime_julday+(tsdatatime_len-1)*tsdatatime_dayint)
     $    ,year,month,day)

        if(my_task.eq.master_task) then
          write(6,*) 'tsdata tlength ',tsdatatime_len
          write(6,*) 'tsdata dayint',tsdatatime_dayint
          write(6,*) 'tsdata start_time ', tsdatatime_name
          write(6,*) 'tsdata end_time ',year,month,day 
        end if

        init=1

! close file
        status=nf_close(ncid)
        call handle_error_netcdf('nf_close',status,nf_noerr)

        return

      end if

      time0_tsdata=tsdatatime_julday-julday_start
      if(my_task.eq.master_task) 
     $ write(6,*) 'tsdatatime_julday ',tsdatatime_julday
      if(my_task.eq.master_task) 
     $ write(6,*) 'julday_start ',julday_start

! get number on step-0

      n=(time-dti/86400.-time0_tsdata)/tsdatatime_dayint+1
      if(my_task.eq.master_task) write(6,*) 'time ',time-dti/86400.
      if(my_task.eq.master_task) write(6,*) 'time0_tsdata ',time0_tsdata
      if(my_task.eq.master_task) 
     $ write(6,*) 'tsdatatime_dayint ',tsdatatime_dayint
      if(my_task.eq.master_task) write(6,*) 'tsdata step ',n

      if(n.le.0) then
        if(my_task.eq.master_task) write(*,'(/a)')
     $                                     'wrong time of ts data'
        call finalize_mpi
        stop
      end if

! get variables on step-0
      
      start(1)=1
      count(1)=im_global
      start(2)=1
      count(2)=jm_global
      start(3)=1
      count(3)=kb
      start(4)=n
      count(4)=1

!      status=nf_inq_varid(ncid,'t',varid)
      status=nf_inq_varid(ncid,'temp',varid)
      call handle_error_netcdf('nf_inq_varid: temp'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3d)
      call handle_error_netcdf('nf_get_var_real: temp',
     $                          status,nf_noerr)
      call scatter3d(tref0,glb3d)

!      status=nf_inq_varid(ncid,'s',varid)
      status=nf_inq_varid(ncid,'sal',varid)
      call handle_error_netcdf('nf_inq_varid: sal'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3d)
      call handle_error_netcdf('nf_get_var_real: sal',
     $                          status,nf_noerr)
      call scatter3d(sref0,glb3d)

! get number on step-1

      n=(time-dti/86400.+tsdatatime_dayint-time0_tsdata)
     $   /tsdatatime_dayint+1
      if(my_task.eq.master_task) write(6,*) 'tsdata step ',n

      if(n.le.0) then
        if(my_task.eq.master_task) write(*,'(/a)')
     $                                     'wrong time of ts data'
        call finalize_mpi
        stop
      end if

! get variables on step-1
      
      start(1)=1
      count(1)=im_global
      start(2)=1
      count(2)=jm_global
      start(3)=1
      count(3)=kb
      start(4)=n
      count(4)=1

!      status=nf_inq_varid(ncid,'t',varid)
      status=nf_inq_varid(ncid,'temp',varid)
      call handle_error_netcdf('nf_inq_varid: temp'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3d)
      call handle_error_netcdf('nf_get_var_real: temp',
     $                          status,nf_noerr)
      call scatter3d(tref1,glb3d)

!      status=nf_inq_varid(ncid,'s',varid)
      status=nf_inq_varid(ncid,'sal',varid)
      call handle_error_netcdf('nf_inq_varid: sal'
     $                          ,status,nf_noerr)
      status=nf_get_vara_real(ncid,varid,start,count,glb3d)
      call handle_error_netcdf('nf_get_var_real: sal',
     $                          status,nf_noerr)
      call scatter3d(sref1,glb3d)

! close file
      status=nf_close(ncid)
      call handle_error_netcdf('nf_close',status,nf_noerr)

      return
      end















