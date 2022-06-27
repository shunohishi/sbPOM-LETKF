MODULE letkf_tools
  !=======================================================================
  !
  ! [PURPOSE:] Module for LETKF with POM
  !
  ! [HISTORY:]
  !   01/26/2009 Takemasa Miyoshi  created
  !   02/03/2009 Takemasa Miyoshi  modified for ROMS
  !   01/26/2011 Yasumasa Miyazawa  modified for POM (check 'pom' or 'POM')
  !
  !=======================================================================
  !$USE OMP_LIB
  USE common_setting
  USE common
  USE common_mpi
  USE common_pom
  USE common_mpi_pom
  USE common_letkf
  USE letkf_obs

  IMPLICIT NONE

  PRIVATE
  PUBLIC ::  das_letkf

  INTEGER,SAVE :: nobstotal

CONTAINS
  !-----------------------------------------------------------------------
  ! Data Assimilation
  !-----------------------------------------------------------------------
  SUBROUTINE das_letkf(gues3d,gues2d,anal3d,anal2d)
    IMPLICIT NONE
    CHARACTER(11) :: inflfile='infl_mul.nc'
    REAL(r_size),INTENT(INOUT) :: gues3d(nij1,nlev,nbv,nv3d) ! background ensemble
    REAL(r_size),INTENT(INOUT) :: gues2d(nij1,nbv,nv2d)      !  output: destroyed
    REAL(r_size),INTENT(OUT) :: anal3d(nij1,nlev,nbv,nv3d) ! analysis ensemble
    REAL(r_size),INTENT(OUT) :: anal2d(nij1,nbv,nv2d)
    REAL(r_size),ALLOCATABLE :: mean3d(:,:,:)
    REAL(r_size),ALLOCATABLE :: mean2d(:,:)
    REAL(r_size),ALLOCATABLE :: hdxf(:,:)
    REAL(r_size),ALLOCATABLE :: rdiag(:)
    REAL(r_size),ALLOCATABLE :: rloc(:)
    REAL(r_size),ALLOCATABLE :: dep(:)
    REAL(r_size),ALLOCATABLE :: work3d(:,:,:)
    REAL(r_size),ALLOCATABLE :: work2d(:,:)
    REAL(r_sngl),ALLOCATABLE :: work3dg(:,:,:,:)
    REAL(r_sngl),ALLOCATABLE :: work2dg(:,:,:)
    REAL(r_size),ALLOCATABLE :: depth(:,:)
    REAL(r_size) :: parm
    REAL(r_size) :: trans(nbv,nbv,nv3d+nv2d) !W+w : S.Ohishi 2019.08
    REAL(r_size) :: transw(nbv,nbv) !W: S.Ohishi 2019.08
    REAL(r_size) :: transm(nbv) !w: S.Ohishi 2019.08
    REAL(r_size) :: beta(nv3d) !for RTPS S.Ohishi 2019.08
    REAL(r_size) :: pa(nbv,nbv) !Pa: S.Ohishi 2019.08
    INTEGER :: ij,ilev,n,m,i,j,k,nobsl,ierr
    LOGICAL :: ex

    WRITE(6,'(A)') 'Hello from das_letkf'
    WRITE(6,'(A,F12.5)') 'cov_infl_mul :',cov_infl_mul !SO
    WRITE(6,'(A,F12.5)') 'sp_infl_add  :',sp_infl_add  !SO
    WRITE(6,'(A,F12.5)') 'ALPHA_RTPP   :',ALPHA_RTPP   !SO
    WRITE(6,'(A,F12.5)') 'ALPHA_RTPS   :',ALPHA_RTPS   !SO
    nobstotal = nobs
    WRITE(6,'(A,I8)') 'Target observation numbers : NOBS=',nobs
    !
    ! In case of no obs
    !
    IF(nobstotal == 0) THEN
       WRITE(6,'(A)') 'No observation assimilated'
       anal3d = gues3d
       anal2d = gues2d
       RETURN
    END IF
    !
    ! Check inflation parameter
    ! S. Ohishi 2019.08
    !
    IF((ALPHA_RTPP /= 0.d0 .or. ALPHA_RTPS /= 0.d0) .and. cov_infl_mul /= 1.d0) THEN !apply RTPP/RTPS
       WRITE(6,'(A)') '***Error: Check inflation parameter'
       STOP
    ELSE IF(ALPHA_RTPP == 0.d0 .and. ALPHA_RTPS == 0.d0 .and. cov_infl_mul == 0.d0) THEN !Not apply RTPP/RTPS
       WRITE(6,'(A)') '***Error: Check inflation parameter'
       STOP
    END IF
    !
    ! FCST PERTURBATIONS
    !
    ALLOCATE(mean3d(nij1,nlev,nv3d))
    ALLOCATE(mean2d(nij1,nv2d))
    CALL ensmean_grd(nbv,nij1,gues3d,gues2d,mean3d,mean2d)
    DO n=1,nv3d
       DO m=1,nbv
          DO k=1,nlev
             DO i=1,nij1
                gues3d(i,k,m,n) = gues3d(i,k,m,n) - mean3d(i,k,n)
             END DO
          END DO
       END DO
    END DO
    DO n=1,nv2d
       DO m=1,nbv
          DO i=1,nij1
             gues2d(i,m,n) = gues2d(i,m,n) - mean2d(i,n)
          END DO
       END DO
    END DO
    !
    ! depth
    !
    ALLOCATE(depth(nij1,nlev))
    DO i=1,nij1
       CALL calc_depth(mean2d(i,iv2d_z),phi1(i),depth(i,:))
    END DO
    !
    ! multiplicative inflation
    !
    IF(cov_infl_mul > 0.0d0) THEN ! fixed multiplicative inflation parameter
       ALLOCATE( work3d(nij1,nlev,nv3d) )
       ALLOCATE( work2d(nij1,nv2d) )
       work3d = cov_infl_mul
       work2d = cov_infl_mul
    END IF
    IF(cov_infl_mul <= 0.0d0) THEN ! 3D parameter values are read-in
       ALLOCATE( work3dg(nlon,nlat,nlev,nv3d) )
       ALLOCATE( work2dg(nlon,nlat,nv2d) )
       ALLOCATE( work3d(nij1,nlev,nv3d) )
       ALLOCATE( work2d(nij1,nv2d) )
       INQUIRE(FILE=inflfile,EXIST=ex)
       IF(ex) THEN
          IF(myrank == 0) THEN
             WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is reading.. ',inflfile
             CALL read_grd4(inflfile,work3dg,work2dg)
             !work3dg =1.0
             !work2dg =1.0
             !CALL write_grd4(inflfile,work3dg,work2dg)
             !WRITE(6,'(2A)') '!!ERROR: stop',inflfile
             !STOP
          END IF
          CALL scatter_grd_mpi(0,work3dg,work2dg,work3d,work2d)
       ELSE
          WRITE(6,'(2A)') '!!ERROR: no such file exist: ',inflfile
          STOP
          work3d = -1.0d0 * cov_infl_mul
       END IF
    END IF
    !
    ! MAIN ASSIMILATION LOOP
    !
    ALLOCATE( hdxf(1:nobstotal,1:nbv),rdiag(1:nobstotal),rloc(1:nobstotal),dep(1:nobstotal) )
    DO ilev=1,nlev
       !    WRITE(6,'(A,I3)') 'ilev = ',ilev
       DO ij=1,nij1

          CALL obs_local(ij,ilev,depth(ij,ilev),hdxf,rdiag,rloc,dep,nobsl)
          parm = work3d(ij,ilev,iv3d_t)

          IF(ALPHA_RTPP /= 0.d0)THEN
             CALL letkf_core(nobstotal,nobsl,hdxf,rdiag,rloc,dep,parm,transw,transm)!SO
             CALL weight_RTPP(nbv,transw,transm,trans(:,:,1))                           !SO
             if(nv3d+nv2d > 1)then                                                      !SO
                do n=2,nv3d+nv2d                                                        !SO
                   trans(:,:,n)=trans(:,:,1)                                            !SO
                end do                                                                  !SO
             end if                                                                     !SO
          ELSE IF(ALPHA_RTPS /= 0.d0)THEN                                               !SO
             CALL letkf_core(nobstotal,nobsl,hdxf,rdiag,rloc,dep,parm,transw,transm,pa) !SO
             CALL weight_RTPS(nbv,transw,transm,pa,gues2d(ij,:,:),gues3d(ij,ilev,:,:),trans)!SO 
          ELSE                                                                          !SO
             CALL letkf_core(nobstotal,nobsl,hdxf,rdiag,rloc,dep,parm,trans(:,:,1)) !SO
             if(nv3d+nv2d > 1)then                                                      !SO
                do n=2,nv3d+nv2d                                                        !SO
                   trans(:,:,n)=trans(:,:,1)                                            !SO
                end do                                                                  !SO
             end if                                                                     !SO
          ENDIF

          DO n=1,nv3d
             DO m=1,nbv
                anal3d(ij,ilev,m,n) = mean3d(ij,ilev,n)
                DO k=1,nbv
                   anal3d(ij,ilev,m,n) = anal3d(ij,ilev,m,n) &
                        & + gues3d(ij,ilev,k,n) * trans(k,m,n) !SO
                END DO
             END DO
          END DO
          IF(ilev == nlev) THEN
             DO n=1,nv2d
                ! pom
                !          IF(n==iv2d_ubar .OR. n==iv2d_vbar .OR. n==iv2d_hbl) THEN
                IF(n==iv2d_ubar .OR. n==iv2d_vbar) THEN
                   !
                   DO m=1,nbv
                      anal2d(ij,m,n)  = mean2d(ij,n) + gues2d(ij,m,n)
                   END DO
                ELSE
                   DO m=1,nbv
                      anal2d(ij,m,n)  = mean2d(ij,n)
                      DO k=1,nbv
                         anal2d(ij,m,n) = anal2d(ij,m,n) + gues2d(ij,k,n) * trans(k,m,n) !SO
                      END DO
                   END DO
                END IF
             END DO
          END IF
          IF(cov_infl_mul < 0.0d0) THEN
             work3d(ij,ilev,:) = parm
             IF(ilev == nlev) work2d(ij,:) = parm
          END IF
       END DO
    END DO
    DEALLOCATE(hdxf,rdiag,rloc,dep)
    IF(cov_infl_mul < 0.0d0) THEN
       CALL gather_grd_mpi(0,work3d,work2d,work3dg,work2dg)
       IF(myrank == 0) THEN
          WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing.. ',inflfile
          CALL write_grd4(inflfile,work3dg,work2dg)
       END IF
       DEALLOCATE(work3dg,work2dg,work3d,work2d)
    END IF
    !
    ! Additive inflation
    !
    IF(sp_infl_add > 0.0d0) THEN
       CALL read_ens_mpi('addi',nbv,gues3d,gues2d)
       ALLOCATE( work3d(nij1,nlev,nv3d) )
       ALLOCATE( work2d(nij1,nv2d) )
       CALL ensmean_grd(nbv,nij1,gues3d,gues2d,work3d,work2d)
       DO n=1,nv3d
          DO m=1,nbv
             DO k=1,nlev
                DO i=1,nij1
                   gues3d(i,k,m,n) = gues3d(i,k,m,n) - work3d(i,k,n)
                END DO
             END DO
          END DO
       END DO
       DO n=1,nv2d
          DO m=1,nbv
             DO i=1,nij1
                gues2d(i,m,n) = gues2d(i,m,n) - work2d(i,n)
             END DO
          END DO
       END DO

       DEALLOCATE(work3d,work2d)
       WRITE(6,'(A)') '===== Additive covariance inflation ====='
       WRITE(6,'(A,F10.4)') '  parameter:',sp_infl_add
       WRITE(6,'(A)') '========================================='
       !    parm = 0.7d0
       !    DO ilev=1,nlev
       !      parm_infl_damp(ilev) = 1.0d0 + parm &
       !        & + parm * REAL(1-ilev,r_size)/REAL(nlev_dampinfl,r_size)
       !      parm_infl_damp(ilev) = MAX(parm_infl_damp(ilev),1.0d0)
       !    END DO
       DO n=1,nv3d
          DO m=1,nbv
             DO ilev=1,nlev
                DO ij=1,nij1
                   anal3d(ij,ilev,m,n) = anal3d(ij,ilev,m,n) &
                        & + gues3d(ij,ilev,m,n) * sp_infl_add
                END DO
             END DO
          END DO
       END DO
       DO n=1,nv2d
          DO m=1,nbv
             DO ij=1,nij1
                anal2d(ij,m,n) = anal2d(ij,m,n) + gues2d(ij,m,n) * sp_infl_add
             END DO
          END DO
       END DO
    END IF

    !
    ! Round off using max,min
    ! Added by S.Ohishi 2020.04
    IF(ROFF)then
       DO m=1,nbv
          DO ilev=1,nlev
             DO ij=1,nij1
                IF(anal3d(ij,ilev,m,iv3d_t) == 0.d0 .or. anal3d(ij,ilev,m,iv3d_s) == 0.d0) CYCLE
                anal3d(ij,ilev,m,iv3d_t) = MAX(anal3d(ij,ilev,m,iv3d_t),tmin)
                anal3d(ij,ilev,m,iv3d_t) = MIN(anal3d(ij,ilev,m,iv3d_t),tmax)
                anal3d(ij,ilev,m,iv3d_s) = MAX(anal3d(ij,ilev,m,iv3d_s),smin)
                anal3d(ij,ilev,m,iv3d_s) = MIN(anal3d(ij,ilev,m,iv3d_s),smax)
             END DO
          END DO
       END DO
    END IF
    !
    ! Force non-negative s
    !
    !DO m=1,nbv
    !   DO ilev=1,nlev
    !      DO ij=1,nij1
    !         anal3d(ij,ilev,m,iv3d_s) = MAX(anal3d(ij,ilev,m,iv3d_s),0.0d0)
    !      END DO
    !   END DO
    !END DO

    DEALLOCATE(mean3d,mean2d)
    RETURN
  END SUBROUTINE das_letkf
  !-----------------------------------------------------------------------
  ! Project global observations to local
  !     (hdxf_g,dep_g,rdiag_g) -> (hdxf,dep,rdiag)
  !-----------------------------------------------------------------------
  SUBROUTINE obs_local(ij,ilev,depth,hdxf,rdiag,rloc,dep,nobsl)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: ij,ilev
    REAL(r_size),INTENT(IN) :: depth
    REAL(r_size),INTENT(OUT) :: hdxf(nobstotal,nbv)
    REAL(r_size),INTENT(OUT) :: rdiag(nobstotal)
    REAL(r_size),INTENT(OUT) :: rloc(nobstotal)
    REAL(r_size),INTENT(OUT) :: dep(nobstotal)
    INTEGER,INTENT(OUT) :: nobsl
    REAL(r_size) :: minlon,maxlon,minlat,maxlat,dist,dlev
    INTEGER,ALLOCATABLE:: nobs_use(:)
    INTEGER :: imin,imax,jmin,jmax,im,ichan
    INTEGER :: n,nn
    !
    ! INITIALIZE
    !
    IF( nobs > 0 ) THEN
       ALLOCATE(nobs_use(nobs))
    END IF
    !
    ! data search
    !
    !  S. Ohishi 2019.08
    imin = MAX(NINT(ri1(ij) - dist_zero/dlon1(ij)),1)
    imax = MIN(NINT(ri1(ij) + dist_zero/dlon1(ij)),nlon)
    jmin = MAX(NINT(rj1(ij) - dist_zero/dlat1(ij)),1)
    jmax = MIN(NINT(rj1(ij) + dist_zero/dlat1(ij)),nlat)

    !  imin = MAX(NINT(ri1(ij) - dist_zero),1)
    !  imax = MIN(NINT(ri1(ij) + dist_zero),nlon)
    !  jmin = MAX(NINT(rj1(ij) - dist_zero),1)
    !  jmax = MIN(NINT(rj1(ij) + dist_zero),nlat)

    nn = 1
    IF( nobs > 0 ) CALL obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
    nn = nn-1
    IF(nn < 1) THEN
       nobsl = 0
       RETURN
    END IF
    !
    ! CONVENTIONAL
    !
    nobsl = 0
    IF(nn > 0) THEN
       DO n=1,nn
          IF(NINT(obselm(nobs_use(n))) == id_z_obs .AND. ilev < nlev) THEN
             dlev = ABS(depth)
          ELSE IF(NINT(obselm(nobs_use(n))) /= id_z_obs) THEN
             dlev = ABS(obslev(nobs_use(n)) - depth)
          ELSE
             dlev = 0.0d0
          END IF
          IF(dist_zerov >= 0. .and. dlev > dist_zerov) CYCLE !SO for no vertical localization
!          IF(dlev > dist_zerov) CYCLE

          !S.Ohishi 2019.08
          call com_distll_1(obslon(nobs_use(n)),obslat(nobs_use(n)),lon1(ij),lat1(ij),dist)
          ! dist = SQRT((obsi(nobs_use(n))-ri1(ij))**2 &
          ! & + (obsj(nobs_use(n))-rj1(ij))**2)
          IF( dist > dist_zero ) CYCLE

          nobsl = nobsl + 1
          hdxf(nobsl,:) = obshdxf(nobs_use(n),:)
          dep(nobsl)    = obsdep(nobs_use(n))
          !
          ! Observational localization
          !
          rdiag(nobsl) = obserr(nobs_use(n))**2
          
          !S. Ohishi 2020.04
          !Add No vertical localization
          IF(sigma_obsv <= 0.)THEN
             rloc(nobsl) = EXP(-0.5d0 * (dist/sigma_obs)**2) !Non vertical loc.
          ELSE
             rloc(nobsl) = EXP(-0.5d0 * ((dist/sigma_obs)**2 + (dlev/sigma_obsv)**2)) !original
          ENDIF

       END DO
    END IF
    !
    ! DEBUG
    ! IF( ILEV == 1 .AND. ILON == 1 ) &
    ! & WRITE(6,*) 'ILEV,ILON,ILAT,NN,TVNN,NOBSL=',ilev,ij,nn,tvnn,nobsl
    !
    IF( nobsl > nobstotal ) THEN
       WRITE(6,'(A,I5,A,I5)') 'FATAL ERROR, NOBSL=',nobsl,' > NOBSTOTAL=',nobstotal
       WRITE(6,*) 'IJ,NN=',ij,nn
       STOP 99
    END IF
    !
    IF( nobs > 0 ) THEN
       DEALLOCATE(nobs_use)
    END IF
    !
    RETURN
  END SUBROUTINE obs_local

  SUBROUTINE obs_local_sub(imin,imax,jmin,jmax,nn,nobs_use)
    INTEGER,INTENT(IN) :: imin,imax,jmin,jmax
    INTEGER,INTENT(INOUT) :: nn, nobs_use(nobs)
    INTEGER :: j,n,ib,ie,ip

    DO j=jmin,jmax
       IF(imin > 1) THEN
          ib = nobsgrd(imin-1,j)+1
       ELSE
          IF(j > 1) THEN
             ib = nobsgrd(nlon,j-1)+1
          ELSE
             ib = 1
          END IF
       END IF
       ie = nobsgrd(imax,j)
       n = ie - ib + 1
       IF(n == 0) CYCLE
       DO ip=ib,ie
          IF(nn > nobs) THEN
             WRITE(6,*) 'FATALERROR, NN > NOBS', NN, NOBS
          END IF
          nobs_use(nn) = ip
          nn = nn + 1
       END DO
    END DO

    RETURN
  END SUBROUTINE obs_local_sub

  !------------------------------------------------------------------------------
  ! RTPP |
  !-------
  !
  ! W_inf = alpha*I + (1 - alpha) * W_tmp
  ! W = W_inf + w
  !
  !---------------------------------------
  ! S.Ohishi 2019.08
  !-----------------------------------------------------------------------------

  SUBROUTINE weight_RTPP(nbv,transw,transm,trans)

    IMPLICIT NONE
    INTEGER,INTENT(IN) :: nbv
    REAL(r_size),INTENT(IN) :: transw(nbv,nbv)
    REAL(r_size),INTENT(IN) :: transm(nbv)
    REAL(r_size),INTENT(OUT) :: trans(nbv,nbv)
    INTEGER :: i,j
    
    !W_tmp = (1-alpha) * W_tmp
    do j=1,nbv
       do i=1,nbv
          trans(i,j)=(1.d0-ALPHA_RTPP)*transw(i,j)
       end do
    end do

    !W_inf = alpha*I + W_tmp
    do i=1,nbv
       trans(i,i)=trans(i,i)+ALPHA_RTPP*1.d0
    end do

    !W = W_inf + w
    do j=1,nbv
       do i=1,nbv
          trans(i,j)=trans(i,j)+transm(i)
       end do
    end do

  END SUBROUTINE weight_RTPP

!-------------------------------------------------------------------
! RTPS |
!-------
!
! W_inf = [ alpha * sqrt(X_b * X_b^T)/sqrt( (nbv-1)*X_b * Pa * X_b^T) + 1-alpha ] W_tmp
! W = W_inf + w
!
!---------------------------
! S.Ohishi 2019.08
!-------------------------------------------------------------------

  SUBROUTINE weight_RTPS(nbv,transw,transm,pa,gues2d,gues3d,trans)

    IMPLICIT NONE
    INTEGER,INTENT(IN) :: nbv
    REAL(r_size),INTENT(IN) :: transw(nbv,nbv)
    REAL(r_size),INTENT(IN) :: transm(nbv)
    REAL(r_size),INTENT(IN) :: pa(nbv,nbv)
    REAL(r_size),INTENT(IN) :: gues3d(nbv,nv3d)
    REAL(r_size),INTENT(IN) :: gues2d(nbv,nv2d)
    REAL(r_size),INTENT(OUT) :: trans(nbv,nbv,nv3d+nv2d)

    INTEGER :: i,j,n

    REAL(r_size) :: spr_a,spr_b

    do n=1,nv3d+nv2d
       
       if(1 <= n .and. n <= nv3d) call spread(nbv,gues3d(:,n),pa,spr_a,spr_b)
       if(nv3d+1 <= n .and. n <= nv3d+nv2d) call spread(nbv,gues2d(:,n-nv3d),pa,spr_a,spr_b)
          
       if(spr_a > 0.d0 .and. spr_b > 0.d0)then
       do j=1,nbv
          do i=1,nbv
             !W_inf
             trans(i,j,n)=(ALPHA_RTPS*spr_b/spr_a + 1.d0 - ALPHA_RTPS)*transw(i,j)
             !W = W_inf + w
             trans(i,j,n)=trans(i,j,n)+transm(i)
          end do
       end do
       else
       do j=1,nbv
          do i=1,nbv
             trans(i,j,n)=transw(i,j)+transm(i)
          end do
       end do
       end if

    end do

  END SUBROUTINE weight_RTPS  

  !-----

  SUBROUTINE spread(nbv,gues,pa,spr_a,spr_b)

    IMPLICIT NONE
    INTEGER,INTENT(IN) :: nbv
    REAL(r_size),INTENT(IN) :: gues(nbv)
    REAL(r_size),INTENT(IN) :: pa(nbv,nbv)
    REAL(r_size),INTENT(OUT) :: spr_a,spr_b

    INTEGER :: i,j

    spr_b=0.d0
    do i=1,nbv
       spr_b=spr_b+gues(i)*gues(i)
    end do
    spr_b=sqrt(spr_b)

    spr_a=0.d0
    do j=1,nbv
       do i=1,nbv
          spr_a=spr_a+gues(i)*pa(i,j)*gues(j)
       end do
    end do
    spr_a=sqrt(real(nbv-1.d0,r_size)*spr_a)

  END SUBROUTINE spread

END MODULE letkf_tools
