MODULE common_obs_pom
  !=======================================================================
  !
  ! [PURPOSE:] Observational procedures
  !
  ! [HISTORY:]
  !   01/23/2009 Takemasa MIYOSHI  created
  !   02/03/2009 Takemasa MIYOSHI  modified for ROMS
  !   02/02/2010 Yasumasa MIYAZAWA  modified for POM
  !
  !=======================================================================
  !$USE OMP_LIB
  USE common_setting
  USE common
  USE common_pom

  IMPLICIT NONE
  PUBLIC

CONTAINS
  !-----------------------------------------------------------------------
  ! Transformation from model variables to an observation
  !-----------------------------------------------------------------------
  SUBROUTINE Trans_XtoY(elm,ri,rj,rlev,v3d,v2d,yobs)
    IMPLICIT NONE
    REAL(r_size),INTENT(IN) :: elm
    REAL(r_size),INTENT(IN) :: ri,rj,rlev
    REAL(r_size),INTENT(IN) :: v3d(nlon,nlat,nlev,nv3d)
    REAL(r_size),INTENT(IN) :: v2d(nlon,nlat,nv2d)
    REAL(r_size),INTENT(OUT) :: yobs
    REAL(r_size) :: wk1(1),wk2(1),depth(nlev)
    INTEGER :: k

    !Add Ohishi 2018.10.11
    IF(phi0(NINT(ri),NINT(rj)) == hnosea)THEN
       yobs = undef
       write(6,*) "nosea at ",NINT(ri),NINT(rj)
       return
    END IF

    SELECT CASE (NINT(elm))
    CASE(id_u_obs) ! U
       yobs = v3d(NINT(ri),NINT(rj),nlev,iv3d_u) ! only surface
    CASE(id_v_obs) ! V
       yobs = v3d(NINT(ri),NINT(rj),nlev,iv3d_v) ! only surface
    CASE(id_t_obs) ! T
       wk1(1) = rlev
       ! added by Miyazawa 2011.03.03
       ! modified by Miyazawa 2012.07.02
       ! modified by Ohishi 2018.10.11
       CALL calc_depth(v2d(NINT(ri),NINT(rj),iv2d_z),phi0(NINT(ri),NINT(rj)),depth)
       !       IF(nint(rlev) == 0) THEN
       IF(depth(nlev) <= rlev .and. rlev <= 0.) THEN
          yobs = v3d(NINT(ri),NINT(rj),nlev,iv3d_t) ! SST
       ELSE
          !      write(6,*) 'obschk rlev depths ', rlev, depth(1), depth(nlev)
          IF(depth(1) <= rlev .and. rlev < depth(nlev)) THEN ! depth(1) < ... < depth(nlev) < 0
             CALL com_interp_spline(nlev,depth,v3d(NINT(ri),NINT(rj),:,iv3d_t),1,wk1,wk2)
             yobs = wk2(1)
          ELSE
             yobs = undef
             write(6,*) "rlev, depth(1), depth(nlev) at undef:",rlev,depth(1),depth(nlev)
          END IF
       END IF
    CASE(id_s_obs) ! S
       wk1(1) = rlev
       ! added by Miyazawa 2011.03.03
       ! modified by Miyazawa 2012.07.02
       ! modified by Ohishi 2018.10.11
       CALL calc_depth(v2d(NINT(ri),NINT(rj),iv2d_z),phi0(NINT(ri),NINT(rj)),depth)
!       IF(nint(rlev) == 0) THEN
       IF(depth(nlev) <= rlev .and. rlev <= 0.) THEN
          yobs = v3d(NINT(ri),NINT(rj),nlev,iv3d_s) ! SSS
       ELSE
          !      write(6,*) 'obschk rlev depths ', rlev, depth(1), depth(nlev)
          IF(depth(1) <= rlev .and. rlev < depth(nlev)) THEN ! depth(1) < ... < depth(nlev) < 0
             CALL com_interp_spline(nlev,depth,v3d(NINT(ri),NINT(rj),:,iv3d_s),1,wk1,wk2)
             yobs = wk2(1)
          ELSE
             yobs = undef
             write(6,*) "rlev, depth(1), depth(nlev) at undef:",rlev,depth(1),depth(nlev)
          END IF
       END IF

    CASE(id_z_obs) ! Z
       yobs = v2d(NINT(ri),NINT(rj),iv2d_z)
    END SELECT

    RETURN
  END SUBROUTINE Trans_XtoY
  !-----------------------------------------------------------------------
  ! Interpolation
  !-----------------------------------------------------------------------
  SUBROUTINE itpl_2d(var,ri,rj,var5)
    IMPLICIT NONE
    REAL(r_size),INTENT(IN) :: var(nlon,nlat)
    REAL(r_size),INTENT(IN) :: ri
    REAL(r_size),INTENT(IN) :: rj
    REAL(r_size),INTENT(OUT) :: var5
    REAL(r_size) :: ai,aj
    INTEGER :: i,j

    i = CEILING(ri)
    ai = ri - REAL(i-1,r_size)
    j = CEILING(rj)
    aj = rj - REAL(j-1,r_size)

    IF(i <= nlon) THEN
       var5 = var(i-1,j-1) * (1-ai) * (1-aj) &
            & + var(i  ,j-1) *    ai  * (1-aj) &
            & + var(i-1,j  ) * (1-ai) *    aj  &
            & + var(i  ,j  ) *    ai  *    aj
    ELSE
       var5 = var(i-1,j-1) * (1-ai) * (1-aj) &
            & + var(1  ,j-1) *    ai  * (1-aj) &
            & + var(i-1,j  ) * (1-ai) *    aj  &
            & + var(1  ,j  ) *    ai  *    aj
    END IF

    RETURN
  END SUBROUTINE itpl_2d

  SUBROUTINE itpl_3d(var,ri,rj,rk,var5)
    IMPLICIT NONE
    REAL(r_size),INTENT(IN) :: var(nlon,nlat,nlev)
    REAL(r_size),INTENT(IN) :: ri
    REAL(r_size),INTENT(IN) :: rj
    REAL(r_size),INTENT(IN) :: rk
    REAL(r_size),INTENT(OUT) :: var5
    REAL(r_size) :: ai,aj,ak
    INTEGER :: i,j,k

    i = CEILING(ri)
    ai = ri - REAL(i-1,r_size)
    j = CEILING(rj)
    aj = rj - REAL(j-1,r_size)
    k = CEILING(rk)
    ak = rk - REAL(k-1,r_size)

    IF(i <= nlon) THEN
       var5 = var(i-1,j-1,k-1) * (1-ai) * (1-aj) * (1-ak) &
            & + var(i  ,j-1,k-1) *    ai  * (1-aj) * (1-ak) &
            & + var(i-1,j  ,k-1) * (1-ai) *    aj  * (1-ak) &
            & + var(i  ,j  ,k-1) *    ai  *    aj  * (1-ak) &
            & + var(i-1,j-1,k  ) * (1-ai) * (1-aj) *    ak  &
            & + var(i  ,j-1,k  ) *    ai  * (1-aj) *    ak  &
            & + var(i-1,j  ,k  ) * (1-ai) *    aj  *    ak  &
            & + var(i  ,j  ,k  ) *    ai  *    aj  *    ak
    ELSE
       var5 = var(i-1,j-1,k-1) * (1-ai) * (1-aj) * (1-ak) &
            & + var(1  ,j-1,k-1) *    ai  * (1-aj) * (1-ak) &
            & + var(i-1,j  ,k-1) * (1-ai) *    aj  * (1-ak) &
            & + var(1  ,j  ,k-1) *    ai  *    aj  * (1-ak) &
            & + var(i-1,j-1,k  ) * (1-ai) * (1-aj) *    ak  &
            & + var(1  ,j-1,k  ) *    ai  * (1-aj) *    ak  &
            & + var(i-1,j  ,k  ) * (1-ai) *    aj  *    ak  &
            & + var(1  ,j  ,k  ) *    ai  *    aj  *    ak
    END IF

    RETURN
  END SUBROUTINE itpl_3d
  !-----------------------------------------------------------------------
  ! Monitor departure
  !-----------------------------------------------------------------------
  SUBROUTINE monit_dep(nn,elm,dep,qc)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: nn
    REAL(r_size),INTENT(IN) :: elm(nn)
    REAL(r_size),INTENT(IN) :: dep(nn)
    INTEGER,INTENT(IN) :: qc(nn)
    REAL(r_size) :: rmse_u,rmse_v,rmse_t,rmse_s,rmse_z
    REAL(r_size) :: bias_u,bias_v,bias_t,bias_s,bias_z
    INTEGER :: n,iu,iv,it,is,iz

    rmse_u = 0.0d0
    rmse_v = 0.0d0
    rmse_t = 0.0d0
    rmse_s = 0.0d0
    rmse_z = 0.0d0
    bias_u = 0.0d0
    bias_v = 0.0d0
    bias_t = 0.0d0
    bias_s = 0.0d0
    bias_z = 0.0d0
    iu = 0
    iv = 0
    it = 0
    is = 0
    iz = 0
    DO n=1,nn
       IF(qc(n) /= 1) CYCLE
       SELECT CASE(NINT(elm(n)))
       CASE(id_u_obs)
          rmse_u = rmse_u + dep(n)**2
          bias_u = bias_u + dep(n)
          iu = iu + 1
       CASE(id_v_obs)
          rmse_v = rmse_v + dep(n)**2
          bias_v = bias_v + dep(n)
          iv = iv + 1
       CASE(id_t_obs)
          rmse_t = rmse_t + dep(n)**2
          bias_t = bias_t + dep(n)
          it = it + 1
       CASE(id_s_obs)
          rmse_s = rmse_s + dep(n)**2
          bias_s = bias_s + dep(n)
          is = is + 1
       CASE(id_z_obs)
          rmse_z = rmse_z + dep(n)**2
          bias_z = bias_z + dep(n)
          iz = iz + 1
       END SELECT
    END DO
    IF(iu == 0) THEN
       rmse_u = undef
       bias_u = undef
    ELSE
       rmse_u = SQRT(rmse_u / REAL(iu,r_size))
       bias_u = bias_u / REAL(iu,r_size)
    END IF
    IF(iv == 0) THEN
       rmse_v = undef
       bias_v = undef
    ELSE
       rmse_v = SQRT(rmse_v / REAL(iv,r_size))
       bias_v = bias_v / REAL(iv,r_size)
    END IF
    IF(it == 0) THEN
       rmse_t = undef
       bias_t = undef
    ELSE
       rmse_t = SQRT(rmse_t / REAL(it,r_size))
       bias_t = bias_t / REAL(it,r_size)
    END IF
    IF(is == 0) THEN
       rmse_s = undef
       bias_s = undef
    ELSE
       rmse_s = SQRT(rmse_s / REAL(is,r_size))
       bias_s = bias_s / REAL(is,r_size)
    END IF
    IF(iz == 0) THEN
       rmse_z = undef
       bias_z = undef
    ELSE
       rmse_z = SQRT(rmse_z / REAL(iz,r_size))
       bias_z = bias_z / REAL(iz,r_size)
    END IF

    WRITE(6,'(A)') '== OBSERVATIONAL DEPARTURE ================================='
    WRITE(6,'(5A12)') 'U','V','T','SALT','ZETA'
    WRITE(6,'(5ES12.3)') bias_u,bias_v,bias_t,bias_s,bias_z
    WRITE(6,'(5ES12.3)') rmse_u,rmse_v,rmse_t,rmse_s,rmse_z
    WRITE(6,'(A)') '== NUMBER OF OBSERVATIONS TO BE ASSIMILATED ================'
    WRITE(6,'(5A12)') 'U','V','T','SALT','ZETA'
    WRITE(6,'(5I12)') iu,iv,it,is,iz
    WRITE(6,'(A)') '============================================================'

    RETURN
  END SUBROUTINE monit_dep
  !-----------------------------------------------------------------------
  ! Basic modules for observation input
  !-----------------------------------------------------------------------
  ! 
  ! Modifyied by S.Ohishi 2018.09.21
  !
  !------------------------------------

  SUBROUTINE get_nobs(cfile,nn)

    IMPLICIT NONE
    INCLUDE "netcdf.inc"

    CHARACTER(*),INTENT(IN) :: cfile
    INTEGER,INTENT(OUT) :: nn
    !  REAL(r_sngl) :: wk(6)
    REAL(r_sngl),ALLOCATABLE :: ele(:)
    !  INTEGER :: ios
    INTEGER :: iu,iv,it,is,iz
    !  INTEGER :: iunit
    INTEGER :: status,ncid,dimid,varid
    INTEGER :: n
    !  LOGICAL :: ex

    nn = 0
    iu = 0
    iv = 0
    it = 0
    is = 0
    iz = 0

    status=nf_open(cfile,nf_nowrite,ncid)
    if(status /= nf_noerr)then
       WRITE(6,'(2A)') cfile,' does not exist -- skipped'
    end if

    status=nf_inq_dimid(ncid,"nobs",dimid)
    if(status /= nf_noerr)then
       write(6,*) "***Error: nf_inq_dimid (nobs)"
       stop
    end if
    status=nf_inq_dimlen(ncid,dimid,nn)
    if(status /= nf_noerr)then
       write(6,*) "***Error: nf_inq_dimlen (nobs)"
       stop
    endif

    allocate(ele(nn))

    status=nf_inq_varid(ncid,"ele",varid)
    if(status /= nf_noerr)then
       write(6,*) "***Error: nf_inq_varid (ele)"
       stop
    end if
    status=nf_get_var_real(ncid,varid,ele)
    if(status /= nf_noerr)then
       write(6,*) "***Error: nf_get_var_real (ele)"
       stop
    end if
    status=nf_close(ncid)
    if(status /= nf_noerr)then
       write(6,*) "***Error: nf_close"  
       stop
    end if

    DO n=1,nn

       SELECT CASE(NINT(ele(n)))
       CASE(id_u_obs)
          iu = iu + 1
       CASE(id_v_obs)
          iv = iv + 1
       CASE(id_t_obs)
          it = it + 1
       CASE(id_s_obs)
          is = is + 1
       CASE(id_z_obs)
          iz = iz + 1
       END SELECT

    END DO

    deallocate(ele)

    WRITE(6,'(I10,A)') nn,' OBSERVATIONS INPUT'
    WRITE(6,'(A12,I10)') '          U:',iu
    WRITE(6,'(A12,I10)') '          V:',iv
    WRITE(6,'(A12,I10)') '          T:',it
    WRITE(6,'(A12,I10)') '       SALT:',is
    WRITE(6,'(A12,I10)') '       ZETA:',iz

    RETURN
  END SUBROUTINE get_nobs

  !-------------------------------------------------------------

  SUBROUTINE read_obs(cfile,nn,elem,rlon,rlat,rlev,odat,oerr)

    IMPLICIT NONE
    INCLUDE "netcdf.inc"

    CHARACTER(*),INTENT(IN) :: cfile
    INTEGER,INTENT(IN) :: nn
    REAL(r_size),INTENT(OUT) :: elem(nn) ! element number
    REAL(r_size),INTENT(OUT) :: rlon(nn) ! longitude
    REAL(r_size),INTENT(OUT) :: rlat(nn) ! latitude
    REAL(r_size),INTENT(OUT) :: rlev(nn) ! depth [meters]
    REAL(r_size),INTENT(OUT) :: odat(nn)
    REAL(r_size),INTENT(OUT) :: oerr(nn)
    INTEGER :: n
    INTEGER :: status,ncid,varid
    REAL(r_sngl) :: rtmp(nn)

    !open
    status=nf_open(cfile,nf_nowrite,ncid)
    if(status /= nf_noerr)then
       write(6,*) "***Error: nf_open"
       stop
    end if

    !element
    status=nf_inq_varid(ncid,"ele",varid)
    if(status /= nf_noerr)then
       write(6,*) "***Error: nf_inq_varid (ele)"
       stop
    end if
    status=nf_get_var_real(ncid,varid,rtmp)
    if(status /= nf_noerr)then
       write(6,*) "***Error: nf_get_var_real (ele)"
       stop
    end if
    do n=1,nn
       elem(n)=real(rtmp(n),r_size)
    end do

    !longitude
    status=nf_inq_varid(ncid,"lon",varid)
    if(status /= nf_noerr)then
       write(6,*) "***Error: nf_inq_varid (lon)"
       stop
    end if
    status=nf_get_var_real(ncid,varid,rtmp)
    if(status /= nf_noerr)then
       write(6,*) "***Error: nf_get_var_real (lon)"
       stop
    end if
    do n=1,nn
       rlon(n)=real(rtmp(n),r_size)
    end do

    !latitude
    status=nf_inq_varid(ncid,"lat",varid)
    if(status /= nf_noerr)then
       write(6,*) "***Error: nf_inq_varid (lat)"
       stop
    end if
    status=nf_get_var_real(ncid,varid,rtmp)
    if(status /= nf_noerr)then
       write(6,*) "***Error: nf_get_var_real (lat)"  
       stop
    end if
    do n=1,nn
       rlat(n)=real(rtmp(n),r_size)
    end do

    !level
    status=nf_inq_varid(ncid,"lev",varid)
    if(status /= nf_noerr)then
       write(6,*) "***Error: nf_inq_varid (lev)"
       stop
    end if
    status=nf_get_var_real(ncid,varid,rtmp)
    if(status /= nf_noerr)then
       write(6,*) "***Error: nf_get_var_real (lev)"
       stop
    end if
    do n=1,nn
       rlev(n)=real(rtmp(n),r_size)
    end do

    !data
    status=nf_inq_varid(ncid,"dat",varid)
    if(status /= nf_noerr)then
       write(6,*) "***Error: nf_inq_varid (dat)"
       stop
    end if
    status=nf_get_var_real(ncid,varid,rtmp)
    if(status /= nf_noerr)then
       write(6,*) "***Error: nf_get_var_real (dat)"
       stop
    end if
    do n=1,nn
       odat(n)=real(rtmp(n),r_size)
    end do

    !Error
    status=nf_inq_varid(ncid,"err",varid)
    if(status /= nf_noerr)then
       write(6,*) "***Error: nf_inq_varid (err)"
       stop
    end if
    status=nf_get_var_real(ncid,varid,rtmp)
    if(status /= nf_noerr)then
       write(6,*) "***Error: nf_get_var_real (err)"
       stop
    end if
    do n=1,nn
       oerr(n)=real(rtmp(n),r_size)
    end do

    status=nf_close(ncid)
    if(status /= nf_noerr)then
       write(6,*) "***Error: nf_close"
       stop
    end if

    IF(OBS_ERR)THEN
       DO n=1,nn
          
          SELECT CASE(NINT(elem(n)))
          CASE(id_u_obs)
             oerr(n)=u_err
          CASE(id_v_obs)
             oerr(n)=u_err             
          CASE(id_t_obs)
             IF(rlev(n) == 0.d0)THEN
                oerr(n)=sst_err
             ELSE
                oerr(n)=t_err
             ENDIF
          CASE(id_s_obs)
             IF(rlev(n) == 0.d0)THEN
                oerr(n)=sss_err
             ELSE
                oerr(n)=s_err
             ENDIF
          CASE(id_z_obs)
             oerr(n)=ssh_err
          END SELECT

       END DO

       WRITE(6,'(A14,F12.5)') 'OBS ERROR SST:',sst_err
       WRITE(6,'(A14,F12.5)') 'OBS ERROR T  :',t_err
       WRITE(6,'(A14,F12.5)') 'OBS ERROR SSS:',sss_err
       WRITE(6,'(A14,F12.5)') 'OBS ERROR S  :',s_err
       WRITE(6,'(A14,F12.5)') 'OBS ERROR U  :',u_err
       WRITE(6,'(A14,F12.5)') 'OBS ERROR V  :',v_err

    END IF

    RETURN
  END SUBROUTINE read_obs

END MODULE common_obs_pom
