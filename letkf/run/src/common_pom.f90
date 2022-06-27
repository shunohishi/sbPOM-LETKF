MODULE common_pom
!=======================================================================
!
! [PURPOSE:] Common Information for ROMS
!
! [HISTORY:]
!   10/15/2004 Takemasa Miyoshi  created
!   01/23/2009 Takemasa Miyoshi  modified
!   02/02/2009 Takemasa Miyoshi  modified for ROMS
!   02/01/2010 Yasumasa Miyazawa  modified for POM
!
!=======================================================================
!$USE OMP_LIB
  USE common_setting
  USE common
  IMPLICIT NONE
  PUBLIC

  REAL(r_size),SAVE :: lon(nlon,nlat)
  REAL(r_size),SAVE :: lat(nlon,nlat)
  REAL(r_size),SAVE :: dlon(nlon,nlat) !S.Ohishi 2019.08
  REAL(r_size),SAVE :: dlat(nlon,nlat)
  REAL(r_size),SAVE :: fsm(nlon,nlat)
  REAL(r_size),SAVE :: phi0(nlon,nlat)
  REAL(r_size),SAVE :: sigma_level0(nlev)
  REAL(r_size),SAVE :: sigma_level(nlev)
  CHARACTER(4),SAVE :: element(nv3d+nv2d)

CONTAINS
!-----------------------------------------------------------------------
! Set the parameters
!-----------------------------------------------------------------------
SUBROUTINE set_common_pom
  IMPLICIT NONE
  INCLUDE 'netcdf.inc'
  INTEGER :: ncid,istat,varid,k

  WRITE(6,'(A)') 'Hello from set_common_pom'
  !
  ! Elements
  !
  element(iv3d_u) = 'U   '
  element(iv3d_v) = 'V   '
  element(iv3d_t) = 'T   '
  element(iv3d_s) = 'SALT'
  element(nv3d+iv2d_z) = 'ZETA'
  element(nv3d+iv2d_ubar) = 'UBAR'
  element(nv3d+iv2d_vbar) = 'VBAR'
  !
  ! Lon, Lat, f, orography
  !
  WRITE(6,'(A)') '  >> accessing to file: grd.nc'
  istat = NF_OPEN('grd.nc',NF_NOWRITE,ncid)
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF OPEN ERROR'
    STOP
  END IF
  istat = NF_INQ_VARID(ncid,'east_e',varid)
  istat = NF_GET_VAR_DOUBLE(ncid,varid,lon)
  istat = NF_INQ_VARID(ncid,'north_e',varid)
  istat = NF_GET_VAR_DOUBLE(ncid,varid,lat)
  istat = NF_INQ_VARID(ncid,'h',varid)
  istat = NF_GET_VAR_DOUBLE(ncid,varid,phi0)
  istat = NF_INQ_VARID(ncid,'fsm',varid) !S.Ohishi 2020.05
  istat = NF_GET_VAR_DOUBLE(ncid,varid,fsm)
  istat = NF_INQ_VARID(ncid,'dx',varid) !S.Ohishi 2019.08
  istat = NF_GET_VAR_DOUBLE(ncid,varid,dlon) 
  istat = NF_INQ_VARID(ncid,'dy',varid) !S.Ohishi 2019.08
  istat = NF_GET_VAR_DOUBLE(ncid,varid,dlat)
  istat = NF_INQ_VARID(ncid,'zz',varid)
  istat = NF_GET_VAR_DOUBLE(ncid,varid,sigma_level0)
  DO k=1,nlev
    sigma_level(nlev-k+1)=sigma_level0(k)
  END DO
  istat = NF_CLOSE(ncid)

  RETURN
END SUBROUTINE set_common_pom
!-----------------------------------------------------------------------
! File I/O (netCDF)
!-----------------------------------------------------------------------
!-- Read a grid file ---------------------------------------------------
! 2018.10 S.Ohishi modified
SUBROUTINE read_grd(filename,v3d,v2d)
  IMPLICIT NONE
  INCLUDE 'netcdf.inc'
  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_size),INTENT(OUT) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_size),INTENT(OUT) :: v2d(nlon,nlat,nv2d)
  REAL(r_sngl) :: tmp2d(nlon,nlat),tmp3d(nlon,nlat,nlev)
  !  REAL(r_sngl) :: buf4(nlon,nlat,nlev)
  INTEGER :: i,j,k
  INTEGER :: ncid,istat,varid

  istat = NF_OPEN(filename,NF_NOWRITE,ncid)
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF OPEN ERROR'
    STOP
  END IF
  !!! z
  istat = NF_INQ_VARID(ncid,'elb',varid)
  istat = NF_GET_VARA_REAL(ncid,varid,(/1,1,1/),(/nlon,nlat,1/),tmp2d)
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF READ ERROR (zeta)'
    STOP
  END IF
  DO j=1,nlat
    DO i=1,nlon
      v2d(i,j,iv2d_z) = REAL(tmp2d(i,j),r_size)
    END DO
  END DO
  !!! ubar
  istat = NF_INQ_VARID(ncid,'uab',varid)
  istat = NF_GET_VARA_REAL(ncid,varid,(/1,1,1/),(/nlon,nlat,1/),tmp2d)
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF READ ERROR (ubar)'
    STOP
  END IF
  DO j=1,nlat
    DO i=1,nlon
      v2d(i,j,iv2d_ubar) = REAL(tmp2d(i,j),r_size)
    END DO
  END DO
  !!! vbar
  istat = NF_INQ_VARID(ncid,'vab',varid)
  istat = NF_GET_VARA_REAL(ncid,varid,(/1,1,1/),(/nlon,nlat,1/),tmp2d)
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF READ ERROR (vbar)'
    STOP
  END IF
  DO j=1,nlat
    DO i=1,nlon
      v2d(i,j,iv2d_vbar) = REAL(tmp2d(i,j),r_size)
    END DO
  END DO
  !!! u
  istat = NF_INQ_VARID(ncid,'ub',varid)
  istat = NF_GET_VARA_REAL(ncid,varid,(/1,1,1,1/),(/nlon,nlat,nlev,1/),tmp3d)
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF READ ERROR (u)'
    STOP
  END IF
  DO k=1,nlev
    DO j=1,nlat
      DO i=1,nlon
!        v3d(i,j,k,iv3d_u) = REAL(buf4(i,j,k),r_size)
! POM -> ROMS type vertical order 2010.08.31
        v3d(i,j,nlev-k+1,iv3d_u) = REAL(tmp3d(i,j,k),r_size)
      END DO
    END DO
  END DO
  !!! v
  istat = NF_INQ_VARID(ncid,'vb',varid)
  istat = NF_GET_VARA_REAL(ncid,varid,(/1,1,1,1/),(/nlon,nlat,nlev,1/),tmp3d)
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF READ ERROR (v)'
    STOP
  END IF
  DO k=1,nlev
    DO j=1,nlat
      DO i=1,nlon
!        v3d(i,j,k,iv3d_v) = REAL(buf4(i,j,k),r_size)
! POM -> ROMS type vertical order 2010.08.31
        v3d(i,j,nlev-k+1,iv3d_v) = REAL(tmp3d(i,j,k),r_size)
      END DO
    END DO
  END DO
  !!! t
  istat = NF_INQ_VARID(ncid,'tb',varid)
  istat = NF_GET_VARA_REAL(ncid,varid,(/1,1,1,1/),(/nlon,nlat,nlev,1/),tmp3d)
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF READ ERROR (temp)'
    STOP
  END IF
  DO k=1,nlev
    DO j=1,nlat
      DO i=1,nlon
!        v3d(i,j,k,iv3d_t) = REAL(buf4(i,j,k),r_size)
! POM -> ROMS type vertical order 2010.08.31
        v3d(i,j,nlev-k+1,iv3d_t) = REAL(tmp3d(i,j,k),r_size)
      END DO
    END DO
  END DO
  !!! s
  istat = NF_INQ_VARID(ncid,'sb',varid)
  istat = NF_GET_VARA_REAL(ncid,varid,(/1,1,1,1/),(/nlon,nlat,nlev,1/),tmp3d)
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF READ ERROR (salt)'
    STOP
  END IF
  DO k=1,nlev
    DO j=1,nlat
      DO i=1,nlon
!        v3d(i,j,k,iv3d_s) = REAL(buf4(i,j,k),r_size)
! POM -> ROMS type vertical order 2010.08.31
        v3d(i,j,nlev-k+1,iv3d_s) = REAL(tmp3d(i,j,k),r_size)
      END DO
    END DO
  END DO

  istat = NF_CLOSE(ncid)

  RETURN
END SUBROUTINE read_grd

SUBROUTINE read_grd4(filename,v3d,v2d,v3dl,v2dl)
  IMPLICIT NONE
  INCLUDE 'netcdf.inc'
  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_sngl),INTENT(OUT) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(OUT) :: v2d(nlon,nlat,nv2d)
  REAL(r_sngl),INTENT(OUT),OPTIONAL :: v3dl(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(OUT),OPTIONAL :: v2dl(nlon,nlat,nv2d)
  REAL(r_sngl) :: tmp2d(nlon,nlat),tmp3d(nlon,nlat,nlev)
  INTEGER :: i,j,k
  INTEGER :: ncid,istat,varid
  INTEGER :: iv2d,iv3d

  istat = NF_OPEN(filename,NF_NOWRITE,ncid)
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF OPEN ERROR'
    STOP
  END IF
  !!! z
  istat = NF_INQ_VARID(ncid,'elb',varid)
  istat = NF_GET_VARA_REAL(ncid,varid,(/1,1,1/),(/nlon,nlat,1/),tmp2d)
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF READ ERROR (zeta)'
    STOP
  END IF
  DO j=1,nlat
     DO i=1,nlon
        v2d(i,j,iv2d_z)=tmp2d(i,j)
     END DO
  END DO
  !!! ubar
  istat = NF_INQ_VARID(ncid,'uab',varid)
  istat = NF_GET_VARA_REAL(ncid,varid,(/1,1,1/),(/nlon,nlat,1/),tmp2d)
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF READ ERROR (ubar)'
    STOP
  END IF
  DO j=1,nlat
     DO i=1,nlon
        v2d(i,j,iv2d_ubar)=tmp2d(i,j)
     END DO
  END DO
  !!! vbar
  istat = NF_INQ_VARID(ncid,'vab',varid)
  istat = NF_GET_VARA_REAL(ncid,varid,(/1,1,1/),(/nlon,nlat,1/),tmp2d)
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF READ ERROR (vbar)'
    STOP
  END IF
  DO j=1,nlat
     DO i=1,nlon
        v2d(i,j,iv2d_vbar)=tmp2d(i,j)
     END DO
  END DO
  !!! u
  istat = NF_INQ_VARID(ncid,'ub',varid)
!  istat = NF_GET_VAR_REAL(ncid,varid,v3d(:,:,:,iv3d_u))
  istat = NF_GET_VARA_REAL(ncid,varid,(/1,1,1,1/),(/nlon,nlat,nlev,1/),tmp3d)
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF READ ERROR (u)'
    STOP
  END IF
  DO k=1,nlev
    DO j=1,nlat
      DO i=1,nlon
! POM -> ROMS type vertical order 2010.08.31
        v3d(i,j,nlev-k+1,iv3d_u) = tmp3d(i,j,k)
      END DO
    END DO
  END DO
  !!! v
  istat = NF_INQ_VARID(ncid,'v',varid)
!  istat = NF_GET_VAR_REAL(ncid,varid,v3d(:,:,:,iv3d_v))
  istat = NF_GET_VARA_REAL(ncid,varid,(/1,1,1,1/),(/nlon,nlat,nlev,1/),tmp3d)
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF READ ERROR (v)'
    STOP
  END IF
  DO k=1,nlev
    DO j=1,nlat
      DO i=1,nlon
! POM -> ROMS type vertical order 2010.08.31
        v3d(i,j,nlev-k+1,iv3d_v) = tmp3d(i,j,k)
      END DO
    END DO
  END DO
  !!! t
  istat = NF_INQ_VARID(ncid,'tb',varid)
!  istat = NF_GET_VAR_REAL(ncid,varid,v3d(:,:,:,iv3d_t))
  istat = NF_GET_VARA_REAL(ncid,varid,(/1,1,1,1/),(/nlon,nlat,nlev,1/),tmp3d)
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF READ ERROR (temp)'
    STOP
  END IF
  DO k=1,nlev
    DO j=1,nlat
      DO i=1,nlon
! POM -> ROMS type vertical order 2010.08.31
        v3d(i,j,nlev-k+1,iv3d_t) = tmp3d(i,j,k)
      END DO
    END DO
  END DO
  !!! s
  istat = NF_INQ_VARID(ncid,'sb',varid)
!  istat = NF_GET_VAR_REAL(ncid,varid,v3d(:,:,:,iv3d_s))
  istat = NF_GET_VARA_REAL(ncid,varid,(/1,1,1,1/),(/nlon,nlat,nlev,1/),tmp3d)
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF READ ERROR (salt)'
    STOP
  END IF
  DO k=1,nlev
    DO j=1,nlat
      DO i=1,nlon
! POM -> ROMS type vertical order 2010.08.31
        v3d(i,j,nlev-k+1,iv3d_s) = tmp3d(i,j,k)
      END DO
    END DO
  END DO

  istat = NF_CLOSE(ncid)

  !GAUSS FILTER
  IF(DL)THEN
     v2dl(:,:,:)=v2d(:,:,:)
     DO iv2d=1,nv2d
        call gauss_filter(v2dl(:,:,iv2d))
     END DO
     
     v3dl(:,:,:,:)=v3d(:,:,:,:)
     DO iv3d=1,nv3d
        DO k=1,nlat
           call gauss_filter(v3dl(:,:,k,iv3d))
        END DO
     END DO
  END IF

  RETURN
END SUBROUTINE read_grd4
!-- Write a grid file -------------------------------------------------
SUBROUTINE write_grd(filename,v3d,v2d)
  IMPLICIT NONE
  INCLUDE 'netcdf.inc'
  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_size),INTENT(IN) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nlon,nlat,nv2d)
  REAL(r_sngl) :: tmp2d(nlon,nlat),tmp3d(nlon,nlat,nlev)
  INTEGER :: i,j,k
  INTEGER :: ncid,istat,varid

  istat = NF_OPEN(filename,NF_WRITE,ncid)
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF OPEN ERROR'
    STOP
  END IF
  !!! z
  istat = NF_INQ_VARID(ncid,'elb',varid)
  DO j=1,nlat
    DO i=1,nlon
      tmp2d(i,j) = REAL(v2d(i,j,iv2d_z),r_sngl)
    END DO
  END DO
  istat = NF_PUT_VARA_REAL(ncid,varid,(/1,1,1/),(/nlon,nlat,1/),tmp2d)
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF WRITE ERROR (zeta)'
    STOP
  END IF
  !!! ubar
  istat = NF_INQ_VARID(ncid,'uab',varid)
  DO j=1,nlat
    DO i=1,nlon
      tmp2d(i,j) = REAL(v2d(i,j,iv2d_ubar),r_sngl)
    END DO
  END DO
  istat = NF_PUT_VARA_REAL(ncid,varid,(/1,1,1/),(/nlon,nlat,1/),tmp2d)
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF WRITE ERROR (ubar)'
    STOP
  END IF
  !!! vbar
  istat = NF_INQ_VARID(ncid,'vab',varid)
  DO j=1,nlat
    DO i=1,nlon
      tmp2d(i,j) = REAL(v2d(i,j,iv2d_vbar),r_sngl)
    END DO
  END DO
  istat = NF_PUT_VARA_REAL(ncid,varid,(/1,1,1/),(/nlon,nlat,1/),tmp2d)
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF WRITE ERROR (vbar)'
    STOP
  END IF
  !!! u
  istat = NF_INQ_VARID(ncid,'ub',varid)
  DO k=1,nlev
    DO j=1,nlat
      DO i=1,nlon
!        buf4(i,j,k) = REAL(v3d(i,j,k,iv3d_u),r_sngl)
! ROMS -> POM type vertical order 2010.08.31
        tmp3d(i,j,k) = REAL(v3d(i,j,nlev-k+1,iv3d_u),r_sngl)
      END DO
    END DO
  END DO
  istat = NF_PUT_VARA_REAL(ncid,varid,(/1,1,1,1/),(/nlon,nlat,nlev,1/),tmp3d)
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF WRITE ERROR (u)'
    STOP
  END IF
  !!! v
  istat = NF_INQ_VARID(ncid,'vb',varid)
  DO k=1,nlev
    DO j=1,nlat
      DO i=1,nlon
!        buf4(i,j,k) = REAL(v3d(i,j,k,iv3d_v),r_sngl)
! ROMS -> POM type vertical order 2010.08.31
        tmp3d(i,j,k) = REAL(v3d(i,j,nlev-k+1,iv3d_v),r_sngl)
      END DO
    END DO
  END DO
  istat = NF_PUT_VARA_REAL(ncid,varid,(/1,1,1,1/),(/nlon,nlat,nlev,1/),tmp3d)
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF WRITE ERROR (v)'
    STOP
  END IF
  !!! t
  istat = NF_INQ_VARID(ncid,'tb',varid)
  DO k=1,nlev
    DO j=1,nlat
      DO i=1,nlon
!        buf4(i,j,k) = REAL(v3d(i,j,k,iv3d_t),r_sngl)
! ROMS -> POM type vertical order 2010.08.31
        tmp3d(i,j,k) = REAL(v3d(i,j,nlev-k+1,iv3d_t),r_sngl)
      END DO
    END DO
  END DO
  istat = NF_PUT_VARA_REAL(ncid,varid,(/1,1,1,1/),(/nlon,nlat,nlev,1/),tmp3d)
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF WRITE ERROR (temp)'
    STOP
  END IF
  !!! s
  istat = NF_INQ_VARID(ncid,'sb',varid)
  DO k=1,nlev
    DO j=1,nlat
      DO i=1,nlon
!        buf4(i,j,k) = REAL(v3d(i,j,k,iv3d_s),r_sngl)
! ROMS -> POM type vertical order 2010.08.31
        tmp3d(i,j,k) = REAL(v3d(i,j,nlev-k+1,iv3d_s),r_sngl)
      END DO
    END DO
  END DO
  istat = NF_PUT_VARA_REAL(ncid,varid,(/1,1,1,1/),(/nlon,nlat,nlev,1/),tmp3d)
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF WRITE ERROR (salt)'
    STOP
  END IF

  istat = NF_CLOSE(ncid)

  RETURN
END SUBROUTINE write_grd

SUBROUTINE write_grd4(filename,v3d,v2d)
  IMPLICIT NONE
  INCLUDE 'netcdf.inc'
  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_sngl),INTENT(IN) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(IN) :: v2d(nlon,nlat,nv2d)
  REAL(r_sngl) :: tmp2d(nlon,nlat),tmp3d(nlon,nlat,nlev)
  !  REAL(r_sngl) :: buf4(nlon,nlat,nlev)
  INTEGER :: i,j,k
  INTEGER :: ncid,istat,varid

  istat = NF_OPEN(filename,NF_WRITE,ncid)
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF OPEN ERROR'
    STOP
  END IF
  !!! z
  DO j=1,nlat
     DO i=1,nlon
        tmp2d(i,j)=v2d(i,j,iv2d_z)
     END DO
  END DO
  istat = NF_INQ_VARID(ncid,'elb',varid)
  istat = NF_PUT_VARA_REAL(ncid,varid,(/1,1,1/),(/nlon,nlat,1/),tmp2d)
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF WRITE ERROR (zeta)'
    STOP
  END IF
  !!! ubar
  DO j=1,nlat
     DO i=1,nlon
        tmp2d(i,j)=v2d(i,j,iv2d_ubar)
     END DO
  END DO
  istat = NF_INQ_VARID(ncid,'uab',varid)
  istat = NF_PUT_VARA_REAL(ncid,varid,(/1,1,1/),(/nlon,nlat,1/),tmp2d)
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF WRITE ERROR (ubar)'
    STOP
  END IF
  !!! vbar
  DO j=1,nlat
     DO i=1,nlon
        tmp2d(i,j)=v2d(i,j,iv2d_vbar)
     END DO
  END DO
  istat = NF_INQ_VARID(ncid,'vab',varid)
  istat = NF_PUT_VARA_REAL(ncid,varid,(/1,1,1/),(/nlon,nlat,1/),tmp2d)
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF WRITE ERROR (vbar)'
    STOP
  END IF
  !!! u
  DO k=1,nlev
    DO j=1,nlat
      DO i=1,nlon
! ROMS -> POM type vertical order 2010.08.31
        tmp3d(i,j,k) = v3d(i,j,nlev-k+1,iv3d_u)
      END DO
    END DO
  END DO
  istat = NF_INQ_VARID(ncid,'ub',varid)
  istat = NF_PUT_VARA_REAL(ncid,varid,(/1,1,1,1/),(/nlon,nlat,nlev,1/),tmp3d)
!  istat = NF_PUT_VAR_REAL(ncid,varid,v3d(:,:,:,iv3d_u))
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF WRITE ERROR (u)'
    STOP
  END IF
  !!! v
  DO k=1,nlev
    DO j=1,nlat
      DO i=1,nlon
! ROMS -> POM type vertical order 2010.08.31
        tmp3d(i,j,k) = v3d(i,j,nlev-k+1,iv3d_v)
      END DO
    END DO
  END DO
  istat = NF_INQ_VARID(ncid,'vb',varid)
  istat = NF_PUT_VARA_REAL(ncid,varid,(/1,1,1,1/),(/nlon,nlat,nlev,1/),tmp3d)
!  istat = NF_PUT_VAR_REAL(ncid,varid,v3d(:,:,:,iv3d_v))
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF WRITE ERROR (v)'
    STOP
  END IF
  !!! t
  DO k=1,nlev
    DO j=1,nlat
      DO i=1,nlon
! ROMS -> POM type vertical order 2010.08.31
        tmp3d(i,j,k) = v3d(i,j,nlev-k+1,iv3d_t)
      END DO
    END DO
  END DO
  istat = NF_INQ_VARID(ncid,'tb',varid)
  istat = NF_PUT_VARA_REAL(ncid,varid,(/1,1,1,1/),(/nlon,nlat,nlev,1/),tmp3d)
!  istat = NF_PUT_VAR_REAL(ncid,varid,v3d(:,:,:,iv3d_t))
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF WRITE ERROR (temp)'
    STOP
  END IF
  !!! s
  DO k=1,nlev
    DO j=1,nlat
      DO i=1,nlon
! ROMS -> POM type vertical order 2010.08.31
        tmp3d(i,j,k) = v3d(i,j,nlev-k+1,iv3d_s)
      END DO
    END DO
  END DO
  istat = NF_INQ_VARID(ncid,'sb',varid)
  istat = NF_PUT_VARA_REAL(ncid,varid,(/1,1,1,1/),(/nlon,nlat,nlev,1/),tmp3d)
!  istat = NF_PUT_VAR_REAL(ncid,varid,v3d(:,:,:,iv3d_s))
  IF(istat /= NF_NOERR) THEN
    WRITE(6,'(A)') 'netCDF WRITE ERROR (salt)'
    STOP
  END IF

  istat = NF_CLOSE(ncid)

  RETURN
END SUBROUTINE write_grd4
!-----------------------------------------------------------------------
! Depth
!-----------------------------------------------------------------------
SUBROUTINE calc_depth(zeta,bottom,depth)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: zeta
  REAL(r_size),INTENT(IN) :: bottom
  REAL(r_size),INTENT(OUT) :: depth(nlev)
  INTEGER :: k

  DO k=1,nlev
    depth(k) = sigma_level(k)*(bottom+zeta)
  END DO

  RETURN
END SUBROUTINE calc_depth
!-----------------------------------------------------------------------
! Monitor
!-----------------------------------------------------------------------
SUBROUTINE monit_grd(v3d,v2d)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: v3d(nlon,nlat,nlev,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nlon,nlat,nv2d)
  INTEGER :: k,n

  DO k=1,nlev
    WRITE(6,'(I2,A)') k,'th level'
    DO n=1,nv3d
      WRITE(6,'(A,2ES10.2)') element(n),MAXVAL(v3d(:,:,k,n)),MINVAL(v3d(:,:,k,n))
    END DO
  END DO

  DO n=1,nv2d
    WRITE(6,'(A,2ES10.2)') element(nv3d+n),MAXVAL(v2d(:,:,n)),MINVAL(v2d(:,:,n))
  END DO

  RETURN
END SUBROUTINE monit_grd
!-----------------------------------------------------------------------
! Ensemble manipulations
!-----------------------------------------------------------------------
SUBROUTINE ensmean_grd(member,nij,v3d,v2d,v3dm,v2dm)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: member
  INTEGER,INTENT(IN) :: nij
  REAL(r_size),INTENT(IN) :: v3d(nij,nlev,member,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nij,member,nv2d)
  REAL(r_size),INTENT(OUT) :: v3dm(nij,nlev,nv3d)
  REAL(r_size),INTENT(OUT) :: v2dm(nij,nv2d)
  INTEGER :: i,k,m,n

  DO n=1,nv3d
!$OMP PARALLEL DO PRIVATE(i,k,m)
    DO k=1,nlev
      DO i=1,nij
        v3dm(i,k,n) = v3d(i,k,1,n)
        DO m=2,member
          v3dm(i,k,n) = v3dm(i,k,n) + v3d(i,k,m,n)
        END DO
        v3dm(i,k,n) = v3dm(i,k,n) / REAL(member,r_size)
      END DO
    END DO
!$OMP END PARALLEL DO
  END DO

  DO n=1,nv2d
    DO i=1,nij
      v2dm(i,n) = v2d(i,1,n)
      DO m=2,member
        v2dm(i,n) = v2dm(i,n) + v2d(i,m,n)
      END DO
      v2dm(i,n) = v2dm(i,n) / REAL(member,r_size)
    END DO
  END DO

  RETURN
END SUBROUTINE ensmean_grd

!--------------------------------------------------------------------------
! Gauss filter |
!--------------------------------------------------------------------------
!
! sigma: exp(-r^2/sigma^2) in window function
! Apply Gaussian filter within exp(-r^2/sigma^2) > 0.1
!
!--------------------------------------------------------------------------

subroutine gauss_filter(dat)

  implicit none

  real,parameter :: rmiss=0.e0
  
  integer ilon1,ilat1,ilon2,ilat2
  integer idlon(nlat),idlat
  integer pass(nlon,nlat),itmp

  real(r_size) r
  real(r_sngl) dat(nlon,nlat)
  real(r_sngl) sum1(nlon,nlat),sum2(nlon,nlat)
  real(r_sngl) tmp(nlon,nlat)  

  !---Estimate calculating idlon,idlat
  
  idlat=nint(g_sigma*dlog(10.d0)*180.d0/(pi*re*abs(lat(1,2)-lat(1,1))))

  !$OMP PARALLEL DO PRIVATE(ilat1,ilat2,itmp)
  do ilat1=1,nlat

     idlon(ilat1)=0

     do ilat2=ilat1-idlat,ilat1+idlat

        if(ilat2 < 1 .or. nlat < ilat2)cycle
           
        itmp=nint(g_sigma*dlog(10.d0)*180.d0/(pi*re*abs(lon(2,1)-lon(1,1))*cos(lat(1,ilat2)*pi/180.d0)))
        if(idlon(ilat1) < itmp)then
           idlon(ilat1)=itmp
        end if

     end do
  end do
  !$OMP END PARALLEL DO

  write(*,*) "idlat:",idlat
  write(*,*) "idlon:",maxval(idlon)
  
  !---Apply Gaussian filter
  !$OMP PARALLEL DO PRIVATE(ilon1,ilon2,ilat1,ilat2,r)
  do ilat1=1,nlat
     
     if(mod(ilat1,100) == 1)then
        write(*,*) ilat1,"/",nlat,"in Gaussian Filter"
     end if
     
     do ilon1=1,nlon
        
        sum1(ilon1,ilat1)=0.d0
        sum2(ilon1,ilat1)=0.d0
        pass(ilon1,ilat1)=0
        
        if(fsm(ilon1,ilat1) == 0.d0)then
           tmp(ilon1,ilat1)=rmiss
           cycle
        end if
        
        do ilat2=ilat1-idlat,ilat1+idlat
           
           if(ilat2 < 1 .or. nlat < ilat2)cycle
           
           do ilon2=ilon1-idlon(ilat2),ilon1+idlon(ilat2)
              
              if(ilon2 < 1 .or. nlon < ilon2)cycle

              !Calculate distance between (lon(i1),lat(j1)) and (lon(i2),lat(j2))
              call com_distll_1(dble(lon(ilon1,1)),dble(lat(1,ilat1)),dble(lon(ilon2,1)),dble(lat(1,ilat2)),r)
              
              if(r <= g_sigma*sqrt(dlog(10.d0)) .and. fsm(ilon2,ilat2) == 1.d0)then 
                 sum1(ilon1,ilat1)= &
                      & sum1(ilon1,ilat1)+dat(ilon2,ilat2)*exp(-1.d0*r*r/(g_sigma*g_sigma))
                 sum2(ilon1,ilat1)= &
                      & sum2(ilon1,ilat1)+exp(-1.d0*r*r/(g_sigma*g_sigma))
                 pass(ilon1,ilat1)=pass(ilon1,ilat1)+1
              end if
              
           end do !ilon2
        end do !ilat2
        
        if(pass(ilon1,ilat1) == 0)then
           tmp(ilon1,ilat1)=rmiss
        else
           tmp(ilon1,ilat1)=sum1(ilon1,ilat1)/sum2(ilon1,ilat1)
        end if
          
     end do !ilon1
  end do !ilat1
  !$OMP END PARALLEL DO  

  !dat <-- tmp
  dat(:,:)=tmp(:,:)
  
end subroutine gauss_filter

END MODULE common_pom
