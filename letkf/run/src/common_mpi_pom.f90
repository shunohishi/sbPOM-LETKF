MODULE common_mpi_pom
!=======================================================================
!
! [PURPOSE:] MPI procedures
!
! [ATTENTION:]
!   DO NOT COMPILE WITH BOTH INLINE EXPANSION AND OMP OPTIONS TOGETHER
!    (Use ONE if you want, but DON'T USE BOTH AT THE SAME TIME)
!
! [HISTORY:]
!   01/23/2009 Takemasa Miyoshi  created
!   02/03/2009 Takemasa Miyoshi  modified for ROMS
!   02/02/2010 Yasumasa Miyazawa  modified for POM
!
!=======================================================================
!$USE OMP_LIB
  USE common_setting
  USE common
  USE common_mpi
  USE common_pom
  IMPLICIT NONE
  PUBLIC

  INTEGER,SAVE :: nij1
  INTEGER,SAVE :: nij1max
  INTEGER,ALLOCATABLE,SAVE :: nij1node(:)
  REAL(r_size),ALLOCATABLE,SAVE :: phi1(:) 
  REAL(r_size),ALLOCATABLE,SAVE :: lon1(:),lat1(:)
  REAL(r_size),ALLOCATABLE,SAVE :: dlon1(:),dlat1(:) !S.Ohishi 2019.08
  REAL(r_size),ALLOCATABLE,SAVE :: ri1(:),rj1(:)

CONTAINS
SUBROUTINE set_common_mpi_pom
  REAL(r_sngl) :: v3dg(nlon,nlat,nlev,nv3d)
  REAL(r_sngl) :: v2dg(nlon,nlat,nv2d)
  REAL(r_size),ALLOCATABLE :: v3d(:,:,:)
  REAL(r_size),ALLOCATABLE :: v2d(:,:)
  INTEGER :: i,j,n

  WRITE(6,'(A)') 'Hello from set_common_mpi_pom'
  i = MOD(nlon*nlat,nprocs)
  nij1max = (nlon*nlat - i)/nprocs + 1
  IF(myrank < i) THEN
    nij1 = nij1max
  ELSE
    nij1 = nij1max - 1
  END IF
  WRITE(6,'(A,I3.3,A,I6)') 'MYRANK ',myrank,' number of grid points: nij1= ',nij1
  ALLOCATE(nij1node(nprocs))
  DO n=1,nprocs
    IF(n-1 < i) THEN
      nij1node(n) = nij1max
    ELSE
      nij1node(n) = nij1max - 1
    END IF
  END DO

  ALLOCATE(phi1(nij1))
  ALLOCATE(lon1(nij1))
  ALLOCATE(lat1(nij1))
  ALLOCATE(dlon1(nij1)) !S.Ohishi 2019.08
  ALLOCATE(dlat1(nij1))
  ALLOCATE(ri1(nij1))
  ALLOCATE(rj1(nij1))

  ALLOCATE(v3d(nij1,nlev,nv3d))
  ALLOCATE(v2d(nij1,nv2d))
  v3dg(:,:,1,1) = REAL(lon,r_sngl)
  v3dg(:,:,1,2) = REAL(lat,r_sngl)
  DO j=1,nlat
    DO i=1,nlon
      v3dg(i,j,2,1) = REAL(i,r_sngl)
      v3dg(i,j,2,2) = REAL(j,r_sngl)
    END DO
  END DO
  v3dg(:,:,3,1) = REAL(dlon,r_sngl) !S.Ohishi 2019.08
  v3dg(:,:,3,2) = REAL(dlat,r_sngl)
  v2dg(:,:,1) = REAL(phi0,r_sngl)

  CALL scatter_grd_mpi(0,v3dg,v2dg,v3d,v2d)

  lon1(:) = v3d(:,1,1)
  lat1(:) = v3d(:,1,2)
  ri1(:) = v3d(:,2,1)
  rj1(:) = v3d(:,2,2)
  dlon1(:) = v3d(:,3,1) !S.Ohishi 2019.08
  dlat1(:) = v3d(:,3,2)
  phi1 = v2d(:,1)

  RETURN
END SUBROUTINE set_common_mpi_pom
!-----------------------------------------------------------------------
! Scatter gridded data to processes (nrank -> all)
!-----------------------------------------------------------------------
SUBROUTINE scatter_grd_mpi(nrank,v3dg,v2dg,v3d,v2d)
  INTEGER,INTENT(IN) :: nrank
  REAL(r_sngl),INTENT(IN) :: v3dg(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(IN) :: v2dg(nlon,nlat,nv2d)
  REAL(r_size),INTENT(OUT) :: v3d(nij1,nlev,nv3d)
  REAL(r_size),INTENT(OUT) :: v2d(nij1,nv2d)
  REAL(r_sngl) :: bufs(nij1max,nlevall,nprocs)
  REAL(r_sngl) :: bufr(nij1max,nlevall)
  INTEGER :: j,k,n,ierr,n0

  IF(myrank == nrank) THEN
    j=0
    DO n=1,nv3d
      DO k=1,nlev
        j = j+1
        CALL grd_to_buf(v3dg(:,:,k,n),bufs(:,j,:))
      END DO
    END DO

    DO n=1,nv2d
      j = j+1
      CALL grd_to_buf(v2dg(:,:,n),bufs(:,j,:))
    END DO
  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  n = nij1max * nlevall
  n0= n
  CALL MPI_SCATTER(bufs,n ,MPI_REAL,&
                 & bufr,n0,MPI_REAL,nrank,MPI_COMM_WORLD,ierr)

  j=0
  DO n=1,nv3d
    DO k=1,nlev
      j = j+1
      v3d(:,k,n) = REAL(bufr(1:nij1,j),r_size)
    END DO
  END DO

  DO n=1,nv2d
    j = j+1
    v2d(:,n) = REAL(bufr(1:nij1,j),r_size)
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

  RETURN
END SUBROUTINE scatter_grd_mpi
!-----------------------------------------------------------------------
! Gather gridded data (all -> nrank)
!-----------------------------------------------------------------------
SUBROUTINE gather_grd_mpi(nrank,v3d,v2d,v3dg,v2dg)
  INTEGER,INTENT(IN) :: nrank
  REAL(r_size),INTENT(IN) :: v3d(nij1,nlev,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nij1,nv2d)
  REAL(r_sngl),INTENT(OUT) :: v3dg(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(OUT) :: v2dg(nlon,nlat,nv2d)
  REAL(r_sngl) :: bufs(nij1max,nlevall)
  REAL(r_sngl) :: bufr(nij1max,nlevall,nprocs)
  INTEGER :: j,k,n,ierr,n0

  j=0
  DO n=1,nv3d
    DO k=1,nlev
      j = j+1
      bufs(1:nij1,j) = REAL(v3d(:,k,n),r_sngl)
    END DO
  END DO

  DO n=1,nv2d
    j = j+1
    bufs(1:nij1,j) = REAL(v2d(:,n),r_sngl)
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  n = nij1max * nlevall
  n0= n
  CALL MPI_GATHER(bufs,n ,MPI_REAL,&
                & bufr,n0,MPI_REAL,nrank,MPI_COMM_WORLD,ierr)

  IF(myrank == nrank) THEN
    j=0
    DO n=1,nv3d
      DO k=1,nlev
        j = j+1
        CALL buf_to_grd(bufr(:,j,:),v3dg(:,:,k,n))
      END DO
    END DO

    DO n=1,nv2d
      j = j+1
      CALL buf_to_grd(bufr(:,j,:),v2dg(:,:,n))
    END DO
  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

  RETURN
END SUBROUTINE gather_grd_mpi
!-----------------------------------------------------------------------
! Read ensemble data and distribute to processes
!-----------------------------------------------------------------------
SUBROUTINE read_ens_mpi(file,member,v3d,v2d)
  CHARACTER(4),INTENT(IN) :: file
  INTEGER,INTENT(IN) :: member
  REAL(r_size),INTENT(OUT) :: v3d(nij1,nlev,member,nv3d)
  REAL(r_size),INTENT(OUT) :: v2d(nij1,member,nv2d)
  REAL(r_sngl) :: v3dg(nlon,nlat,nlev,nv3d)
  REAL(r_sngl) :: v2dg(nlon,nlat,nv2d)
  INTEGER :: l,n,ll,im
  CHARACTER(10) :: filename='file000.nc'

  ll = CEILING(REAL(member)/REAL(nprocs))
  DO l=1,ll
    im = myrank+1 + (l-1)*nprocs
    IF(im <= member) THEN
      WRITE(filename(1:7),'(A4,I3.3)') file,im
      WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is reading a file ',filename
      CALL read_grd4(filename,v3dg,v2dg)
    END IF

    DO n=0,nprocs-1
      im = n+1 + (l-1)*nprocs
      IF(im <= member) THEN
        CALL scatter_grd_mpi(n,v3dg,v2dg,v3d(:,:,im,:),v2d(:,im,:))
      END IF
    END DO
  END DO

  RETURN
END SUBROUTINE read_ens_mpi
!-----------------------------------------------------------------------
! Write ensemble data after collecting data from processes
!-----------------------------------------------------------------------
SUBROUTINE write_ens_mpi(file,member,v3d,v2d)
  CHARACTER(4),INTENT(IN) :: file
  INTEGER,INTENT(IN) :: member
  REAL(r_size),INTENT(IN) :: v3d(nij1,nlev,member,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nij1,member,nv2d)
  REAL(r_sngl) :: v3dg(nlon,nlat,nlev,nv3d)
  REAL(r_sngl) :: v2dg(nlon,nlat,nv2d)
  INTEGER :: l,n,ll,im
  CHARACTER(10) :: filename='file000.nc'

  ll = CEILING(REAL(member)/REAL(nprocs))
  DO l=1,ll
    DO n=0,nprocs-1
      im = n+1 + (l-1)*nprocs
      IF(im <= member) THEN
        CALL gather_grd_mpi(n,v3d(:,:,im,:),v2d(:,im,:),v3dg,v2dg)
      END IF
    END DO

    im = myrank+1 + (l-1)*nprocs
    IF(im <= member) THEN
      WRITE(filename(1:7),'(A4,I3.3)') file,im
      WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing a file ',filename
      CALL write_grd4(filename,v3dg,v2dg)
    END IF
  END DO

  RETURN
END SUBROUTINE write_ens_mpi
!-----------------------------------------------------------------------
! gridded data -> buffer
!-----------------------------------------------------------------------
SUBROUTINE grd_to_buf(grd,buf)
  REAL(r_sngl),INTENT(IN) :: grd(nlon,nlat)
  REAL(r_sngl),INTENT(OUT) :: buf(nij1max,nprocs)
  INTEGER :: i,j,m,ilon,ilat

  DO m=1,nprocs
    DO i=1,nij1node(m)
      j = m-1 + nprocs * (i-1)
      ilon = MOD(j,nlon) + 1
      ilat = (j-ilon+1) / nlon + 1
      buf(i,m) = grd(ilon,ilat)
    END DO
  END DO

  RETURN
END SUBROUTINE grd_to_buf
!-----------------------------------------------------------------------
! buffer -> gridded data
!-----------------------------------------------------------------------
SUBROUTINE buf_to_grd(buf,grd)
  REAL(r_sngl),INTENT(IN) :: buf(nij1max,nprocs)
  REAL(r_sngl),INTENT(OUT) :: grd(nlon,nlat)
  INTEGER :: i,j,m,ilon,ilat

  DO m=1,nprocs
    DO i=1,nij1node(m)
      j = m-1 + nprocs * (i-1)
      ilon = MOD(j,nlon) + 1
      ilat = (j-ilon+1) / nlon + 1
      grd(ilon,ilat) = buf(i,m)
    END DO
  END DO

  RETURN
END SUBROUTINE buf_to_grd
!-----------------------------------------------------------------------
! STORING DATA (ensemble mean and spread)
!-----------------------------------------------------------------------
SUBROUTINE write_ensmspr_mpi(file,member,v3d,v2d)
  CHARACTER(4),INTENT(IN) :: file
  INTEGER,INTENT(IN) :: member
  REAL(r_size),INTENT(IN) :: v3d(nij1,nlev,member,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nij1,member,nv2d)
  REAL(r_size) :: v3dm(nij1,nlev,nv3d)
  REAL(r_size) :: v2dm(nij1,nv2d)
  REAL(r_size) :: v3ds(nij1,nlev,nv3d)
  REAL(r_size) :: v2ds(nij1,nv2d)
  REAL(r_sngl) :: v3dg(nlon,nlat,nlev,nv3d)
  REAL(r_sngl) :: v2dg(nlon,nlat,nv2d)
  INTEGER :: i,k,m,n
  CHARACTER(10) :: filename='file000.nc'

  CALL ensmean_grd(member,nij1,v3d,v2d,v3dm,v2dm)

  CALL gather_grd_mpi(0,v3dm,v2dm,v3dg,v2dg)
  IF(myrank == 0) THEN
    WRITE(filename(1:7),'(A4,A3)') file,'_me'
    WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing a file ',filename
    CALL write_grd4(filename,v3dg,v2dg)
  END IF

  DO n=1,nv3d
!$OMP PARALLEL DO PRIVATE(i,k,m)
    DO k=1,nlev
      DO i=1,nij1
        v3ds(i,k,n) = (v3d(i,k,1,n)-v3dm(i,k,n))**2
        DO m=2,member
          v3ds(i,k,n) = v3ds(i,k,n) + (v3d(i,k,m,n)-v3dm(i,k,n))**2
        END DO
        v3ds(i,k,n) = SQRT(v3ds(i,k,n) / REAL(member-1,r_size))
      END DO
    END DO
!$OMP END PARALLEL DO
  END DO

  DO n=1,nv2d
!$OMP PARALLEL DO PRIVATE(i,k,m)
    DO i=1,nij1
      v2ds(i,n) = (v2d(i,1,n)-v2dm(i,n))**2
      DO m=2,member
        v2ds(i,n) = v2ds(i,n) + (v2d(i,m,n)-v2dm(i,n))**2
      END DO
      v2ds(i,n) = SQRT(v2ds(i,n) / REAL(member-1,r_size))
    END DO
!$OMP END PARALLEL DO
  END DO

  CALL gather_grd_mpi(0,v3ds,v2ds,v3dg,v2dg)
  IF(myrank == 0) THEN
    WRITE(filename(1:7),'(A4,A3)') file,'_sp'
    WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing a file ',filename
    CALL write_grd4(filename,v3dg,v2dg)
  END IF

  RETURN
END SUBROUTINE write_ensmspr_mpi

END MODULE common_mpi_pom
