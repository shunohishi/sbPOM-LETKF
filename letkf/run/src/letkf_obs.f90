MODULE letkf_obs
  !=======================================================================
  !
  ! [PURPOSE:] Observational procedures
  !
  ! [HISTORY:]
  !   01/23/2009 Takemasa MIYOSHI  created
  !   01/26/2011 Yasumasa MIYAZAWA  modified for POM (check 'pom' or 'POM')
  !
  !=======================================================================
  !$USE OMP_LIB
  USE common_setting
  USE common
  USE common_mpi
  USE common_pom
  USE common_obs_pom
  USE common_mpi_pom
  USE common_letkf

  IMPLICIT NONE
  PUBLIC

  INTEGER,SAVE :: nobs

  REAL(r_size),SAVE :: dist_zero
  REAL(r_size),SAVE :: dist_zerov
  REAL(r_size),ALLOCATABLE,SAVE :: obselm(:)
  REAL(r_size),ALLOCATABLE,SAVE :: obslon(:)
  REAL(r_size),ALLOCATABLE,SAVE :: obslat(:)
  REAL(r_size),ALLOCATABLE,SAVE :: obslev(:)
  REAL(r_size),ALLOCATABLE,SAVE :: obsdat(:)
  REAL(r_size),ALLOCATABLE,SAVE :: obserr(:)
  REAL(r_size),ALLOCATABLE,SAVE :: obsi(:)
  REAL(r_size),ALLOCATABLE,SAVE :: obsj(:)
  REAL(r_size),ALLOCATABLE,SAVE :: obsdep(:)
  REAL(r_size),ALLOCATABLE,SAVE :: obshdxf(:,:)
  INTEGER,SAVE :: nobsgrd(nlon,nlat)

CONTAINS
  !-----------------------------------------------------------------------
  ! Initialize
  !-----------------------------------------------------------------------
  SUBROUTINE set_letkf_obs
    IMPLICIT NONE
    REAL(r_size) :: v3d(nlon,nlat,nlev,nv3d)
    REAL(r_size) :: v2d(nlon,nlat,nv2d)
    ! 2016.02.07
    REAL(r_size),PARAMETER :: gross_error=10.0d0
    REAL(r_size) :: dz,tg,qg
    REAL(r_size) :: obserrt !SO
    REAL(r_size) :: dob2_sprd2 !dob(obsdep)**2 - sprd**2: S0
    REAL(r_size),ALLOCATABLE :: wk2d(:,:)
    INTEGER,ALLOCATABLE :: iwk2d(:,:)
    REAL(r_size),ALLOCATABLE :: tmpelm(:)
    REAL(r_size),ALLOCATABLE :: tmplon(:)
    REAL(r_size),ALLOCATABLE :: tmplat(:)
    REAL(r_size),ALLOCATABLE :: tmplev(:)
    REAL(r_size),ALLOCATABLE :: tmpdat(:)
    REAL(r_size),ALLOCATABLE :: tmperr(:)
    REAL(r_size),ALLOCATABLE :: tmpi(:)
    REAL(r_size),ALLOCATABLE :: tmpj(:)
    REAL(r_size),ALLOCATABLE :: tmpdep(:)
    REAL(r_size),ALLOCATABLE :: tmphdxf(:,:)
    INTEGER,ALLOCATABLE :: tmpqc0(:,:)
    INTEGER,ALLOCATABLE :: tmpqc(:)
    REAL(r_size),ALLOCATABLE :: tmp2elm(:)
    REAL(r_size),ALLOCATABLE :: tmp2lon(:)
    REAL(r_size),ALLOCATABLE :: tmp2lat(:)
    REAL(r_size),ALLOCATABLE :: tmp2lev(:)
    REAL(r_size),ALLOCATABLE :: tmp2dat(:)
    REAL(r_size),ALLOCATABLE :: tmp2err(:)
    REAL(r_size),ALLOCATABLE :: tmp2i(:)
    REAL(r_size),ALLOCATABLE :: tmp2j(:)
    REAL(r_size),ALLOCATABLE :: tmp2dep(:)
    REAL(r_size),ALLOCATABLE :: tmp2hdxf(:,:)
    INTEGER :: nobslots(nslots)
    INTEGER :: n,i,j,ierr,islot,nn,l,im
    INTEGER :: nj(0:nlat-1)
    INTEGER :: njs(1:nlat-1)
    INTEGER :: issh=0,nssh=0,isst=0,nsst=0,isss=0,nsss=0 !SO
    INTEGER :: it=0,nt=0,is=0,ns=0 !SO
    !The number of obs removed by large deviation S.Ohishi 2018.10
    CHARACTER(9) :: obsfile='obsTT.dat'
    CHARACTER(10) :: guesfile='gsTTNNN.nc'

    ! 2012.07.02 Miyazawa
    INTEGER :: numdelqc1 = 0, numdelqc2 = 0, numdelqc3 = 0

    ! 2012.09.10 Miyazawa
    REAL(r_size),ALLOCATABLE :: sprd_hdxf(:) !SO

    WRITE(6,'(A)') 'Hello from set_letkf_obs'

    !S. Ohishi 2020.04
    !Add No vertical localication
    dist_zero = sigma_obs * SQRT(10.0d0/3.0d0) * 2.0d0
    IF(sigma_obsv <= 0)THEN
       dist_zerov = sigma_obsv
    ELSE
       dist_zerov = sigma_obsv * SQRT(10.0d0/3.0d0) * 2.0d0
    ENDIF

    DO islot=1,nslots
       WRITE(obsfile(4:5),'(I2.2)') islot
       CALL get_nobs(obsfile,nobslots(islot))
    END DO
    nobs = SUM(nobslots)
    WRITE(6,'(I10,A)') nobs,' TOTAL OBSERVATIONS INPUT'

    !Add S.Ohishi for no observation
    IF(nobs == 0)THEN
       CALL monit_mean('gues')
       CALL monit_mean('anal')
       WRITE(6,'(A)') 'Not Applying LETKF because of No Observation'

       CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
       CALL finalize_mpi    
       STOP
    END IF

    !
    ! INITIALIZE GLOBAL VARIABLES
    !
    ALLOCATE( tmpelm(nobs) )
    ALLOCATE( tmplon(nobs) )
    ALLOCATE( tmplat(nobs) )
    ALLOCATE( tmplev(nobs) )
    ALLOCATE( tmpdat(nobs) )
    ALLOCATE( tmperr(nobs) )
    ALLOCATE( tmpi(nobs) )
    ALLOCATE( tmpj(nobs) )
    !  ALLOCATE( tmpk(nobs) )
    ALLOCATE( tmpdep(nobs) )
    ALLOCATE( tmphdxf(nobs,nbv) )
    ALLOCATE( tmpqc0(nobs,nbv) )
    ALLOCATE( tmpqc(nobs) )
    tmpqc0 = 0
    tmphdxf = 0.0d0
    !
    ! LOOP of timeslots
    !
    nn=0
    timeslots: DO islot=1,nslots
       IF(nobslots(islot) == 0) CYCLE
       WRITE(obsfile(4:5),'(I2.2)') islot
       CALL read_obs(obsfile,nobslots(islot),&
            & tmpelm(nn+1:nn+nobslots(islot)),tmplon(nn+1:nn+nobslots(islot)),&
            & tmplat(nn+1:nn+nobslots(islot)),tmplev(nn+1:nn+nobslots(islot)),&
            & tmpdat(nn+1:nn+nobslots(islot)),tmperr(nn+1:nn+nobslots(islot)) )
       !2(Accurate)-->1(Fast) Ohishi 2017.12.04
       CALL com_pos2ij(1,nlon,nlat,lon,lat,nobslots(islot), &
            & tmplon(nn+1:nn+nobslots(islot)),tmplat(nn+1:nn+nobslots(islot)),&
            & tmpi  (nn+1:nn+nobslots(islot)),tmpj  (nn+1:nn+nobslots(islot)))
       l=0

       DO
          im = myrank+1 + nprocs * l
          IF(im > nbv) EXIT
          WRITE(guesfile(3:7),'(I2.2,I3.3)') islot,im
          WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is reading a file ',guesfile
          CALL read_grd(guesfile,v3d,v2d)
          !$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(n)
          DO n=1,nobslots(islot)
             !
             ! observational operator
             !

             CALL Trans_XtoY(tmpelm(nn+n),&
                  & tmpi(nn+n),tmpj(nn+n),tmplev(nn+n),v3d,v2d,tmphdxf(nn+n,im))

             ! added by Miyazawa 2011.03.03
             IF(tmphdxf(nn+n,im) > undef) THEN
                tmpqc0(nn+n,im) = 1
             ELSE
                numdelqc3 = numdelqc3 + 1
                tmpqc0(nn+n,im) = 0
             END IF

          END DO
          !$OMP END PARALLEL DO
          l = l+1
       END DO
       nn = nn + nobslots(islot)
    END DO timeslots

    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    ALLOCATE(wk2d(nobs,nbv))
    wk2d = tmphdxf
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    CALL MPI_ALLREDUCE(wk2d,tmphdxf,nobs*nbv,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    DEALLOCATE(wk2d)
    ALLOCATE(iwk2d(nobs,nbv))
    iwk2d = tmpqc0
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    CALL MPI_ALLREDUCE(iwk2d,tmpqc0,nobs*nbv,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
    DEALLOCATE(iwk2d)

    !$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(n,i)
    DO n=1,nobs

       tmpqc(n) = MINVAL(tmpqc0(n,:))

       IF(tmpqc(n) /= 1) THEN
          numdelqc1 = numdelqc1 + 1
          CYCLE
       END IF

       tmpdep(n) = tmphdxf(n,1)
       DO i=2,nbv
          tmpdep(n) = tmpdep(n) + tmphdxf(n,i)
       END DO
       tmpdep(n) = tmpdep(n) / REAL(nbv,r_size) !Hx

       DO i=1,nbv
          tmphdxf(n,i) = tmphdxf(n,i) - tmpdep(n) ! Hdx
       END DO

       tmpdep(n) = tmpdat(n) - tmpdep(n) ! y-Hx

       ! SSHA QC S. Ohishi 2020.04
       IF(NINT(tmpelm(n)) == id_z_obs) THEN

          nssh=nssh+1
          IF(AOEI)THEN
          ELSE IF(tmpdep(n) < hqc_min .or. hqc_max < tmpdep(n))THEN 
             issh=issh+1
             IF(LGE_IO)THEN
                write(6,*) "---SSH---"
                write(6,*) "position(x,y,z):",n,tmplon(n),tmplat(n),tmplev(n)
                write(6,*) "data(y,Hx):",tmpdat(n),tmpdat(n)-tmpdep(n)
             END IF
             numdelqc2 = numdelqc2 + 1
             tmpqc(n) = 0
          END IF

       ! SST/Temperature QC S.Ohishi 2020.04
       ELSE IF(NINT(tmpelm(n)) == id_t_obs) THEN

          IF(tmplev(n) == 0.d0)THEN
             nsst=nsst+1
          ELSE
             nt=nt+1
          ENDIF
          
          IF(AOEI)THEN
          ELSE IF(tmplev(n) == 0.d0)THEN
             IF(tmpdep(n) < sstqc_min .or. sstqc_max < tmpdep(n))THEN
                isst=isst+1
                IF(LGE_IO)THEN
                   write(6,*) "---SST---"
                   write(6,*) "position(x,y,z):",n,tmplon(n),tmplat(n),tmplev(n)
                   write(6,*) "data(y,Hx):",tmpdat(n),tmpdat(n)-tmpdep(n)
                END IF
                numdelqc2 = numdelqc2 + 1
                tmpqc(n) = 0
             ENDIF
          ELSE IF(tmplev(n) /= 0.d0)THEN
             IF(tmpdep(n) < tqc_min .or. tqc_max < tmpdep(n))THEN
                it=it+1
                IF(LGE_IO)THEN
                   write(6,*) "---T---"
                   write(6,*) "position(x,y,z):",n,tmplon(n),tmplat(n),tmplev(n)
                   write(6,*) "data(y,Hx):",tmpdat(n),tmpdat(n)-tmpdep(n)
                END IF
                numdelqc2 = numdelqc2 + 1
                tmpqc(n) = 0
             ENDIF
          END IF

       ! SSS/Salinity QC S.Ohishi 2020.04
       ELSE IF(NINT(tmpelm(n)) == id_s_obs) THEN

          IF(tmplev(n) == 0.d0)THEN
             nsss=nsss+1
          ELSE
             ns=ns+1
          ENDIF
   
          IF(AOEI)THEN
          ELSE IF(tmplev(n) == 0.d0)THEN
             IF(tmpdep(n) < sssqc_min .or. sssqc_max < tmpdep(n))THEN
                isss=isss+1
                IF(LGE_IO)THEN
                   write(6,*) "---SSS---"
                   write(6,*) "position(x,y,z):",n,tmplon(n),tmplat(n),tmplev(n)
                   write(6,*) "data(y,Hx):",tmpdat(n),tmpdat(n)-tmpdep(n)
                END IF
                numdelqc2 = numdelqc2 + 1
                tmpqc(n) = 0
             END IF
          ELSE IF(tmplev(n) /= 0.d0)THEN
             IF(tmpdep(n) < sqc_min .or. sqc_max < tmpdep(n))THEN
                is=is+1
                IF(LGE_IO)THEN
                   write(6,*) "---S---"
                   write(6,*) "position(x,y,z):",n,tmplon(n),tmplat(n),tmplev(n)
                   write(6,*) "data(y,Hx):",tmpdat(n),tmpdat(n)-tmpdep(n)
                END IF
                numdelqc2 = numdelqc2 + 1
                tmpqc(n) = 0
             END IF
          END IF

          !Other variable
       ELSE
          IF(ABS(tmpdep(n)) > gross_error*tmperr(n)) THEN !gross error
             numdelqc2 = numdelqc2 + 1
             tmpqc(n) = 0
          END IF
       END IF

    END DO
    !$OMP END PARALLEL DO
    DEALLOCATE(tmpqc0)

    WRITE(6,'(I10,A)') SUM(tmpqc),' OBSERVATIONS TO BE ASSIMILATED'
    WRITE(6,'(I10,A)') numdelqc1,' deleted by original tmpqc values'
    WRITE(6,'(I10,A)') numdelqc2,' deleted by large deviations'
    WRITE(6,'(A,I10,I10)') "SSH:",issh,nssh
    WRITE(6,'(A,I10,I10)') "SST:",isst,nsst
    WRITE(6,'(A,I10,I10)') "SSS:",isss,nsss
    WRITE(6,'(A,I10,I10)') "T:",it,nt
    WRITE(6,'(A,I10,I10)') "S:",is,ns
    WRITE(6,'(I10,A)') numdelqc3,' deleted by undefined point'

    CALL monit_dep(nobs,tmpelm,tmpdep,tmpqc)
    !
    ! temporal observation localization
    !
    nn = 0
    DO islot=1,nslots
       tmperr(nn+1:nn+nobslots(islot)) = tmperr(nn+1:nn+nobslots(islot)) &
            & * exp(0.25d0 * (REAL(islot-nbslot,r_size) / sigma_obst)**2)
       nn = nn + nobslots(islot)
    END DO
    !
    ! SELECT OBS IN THE NODE
    !
    nn = 0
    DO n=1,nobs

       IF(tmpqc(n) /= 1) CYCLE

       nn = nn+1
       tmpelm(nn) = tmpelm(n)
       tmplon(nn) = tmplon(n)
       tmplat(nn) = tmplat(n)
       tmplev(nn) = tmplev(n)
       tmpdat(nn) = tmpdat(n)
       tmperr(nn) = tmperr(n)
       tmpi(nn) = tmpi(n)
       tmpj(nn) = tmpj(n)
       tmpdep(nn) = tmpdep(n)
       tmphdxf(nn,:) = tmphdxf(n,:)
       tmpqc(nn) = tmpqc(n)
    END DO
    nobs = nn
    WRITE(6,'(I10,A,I3.3)') nobs,' OBSERVATIONS TO BE ASSIMILATED IN MYRANK ',myrank
    !
    ! SORT
    !
    ALLOCATE( tmp2elm(nobs) )
    ALLOCATE( tmp2lon(nobs) )
    ALLOCATE( tmp2lat(nobs) )
    ALLOCATE( tmp2lev(nobs) )
    ALLOCATE( tmp2dat(nobs) )
    ALLOCATE( tmp2err(nobs) )
    ALLOCATE( tmp2i(nobs) )
    ALLOCATE( tmp2j(nobs) )
    ALLOCATE( tmp2dep(nobs) )
    ALLOCATE( tmp2hdxf(nobs,nbv) )
    ALLOCATE( obselm(nobs) )
    ALLOCATE( obslon(nobs) )
    ALLOCATE( obslat(nobs) )
    ALLOCATE( obslev(nobs) )
    ALLOCATE( obsdat(nobs) )
    ALLOCATE( obserr(nobs) )
    ALLOCATE( obsi(nobs) )
    ALLOCATE( obsj(nobs) ) 
    ALLOCATE( obsdep(nobs) )
    ALLOCATE( obshdxf(nobs,nbv) )
    ALLOCATE( sprd_hdxf(nobs) ) !SO

    nobsgrd = 0
    nj = 0
    !$OMP PARALLEL PRIVATE(i,j,n,nn)
    !$OMP DO SCHEDULE(DYNAMIC)
    DO j=1,nlat-1
       DO n=1,nobs
          IF(tmpj(n) < j .OR. j+1 <= tmpj(n)) CYCLE
          nj(j) = nj(j) + 1
       END DO
    END DO
    !$OMP END DO
    !$OMP DO SCHEDULE(DYNAMIC)
    DO j=1,nlat-1
       njs(j) = SUM(nj(0:j-1))
    END DO
    !$OMP END DO
    !$OMP DO SCHEDULE(DYNAMIC)
    DO j=1,nlat-1
       nn = 0
       DO n=1,nobs
          IF(tmpj(n) < j .OR. j+1 <= tmpj(n)) CYCLE
          nn = nn + 1
          tmp2elm(njs(j)+nn) = tmpelm(n)
          tmp2lon(njs(j)+nn) = tmplon(n)
          tmp2lat(njs(j)+nn) = tmplat(n)
          tmp2lev(njs(j)+nn) = tmplev(n)
          tmp2dat(njs(j)+nn) = tmpdat(n)
          tmp2err(njs(j)+nn) = tmperr(n)
          tmp2i(njs(j)+nn) = tmpi(n)
          tmp2j(njs(j)+nn) = tmpj(n)
          tmp2dep(njs(j)+nn) = tmpdep(n)
          tmp2hdxf(njs(j)+nn,:) = tmphdxf(n,:)
       END DO
    END DO
    !$OMP END DO
    !$OMP DO SCHEDULE(DYNAMIC)
    DO j=1,nlat-1
       IF(nj(j) == 0) THEN
          nobsgrd(:,j) = njs(j)
          CYCLE
       END IF
       nn = 0
       DO i=1,nlon
          DO n=njs(j)+1,njs(j)+nj(j)
             IF(tmp2i(n) < i .OR. i+1 <= tmp2i(n)) CYCLE
             nn = nn + 1
             obselm(njs(j)+nn) = tmp2elm(n)
             obslon(njs(j)+nn) = tmp2lon(n)
             obslat(njs(j)+nn) = tmp2lat(n)
             obslev(njs(j)+nn) = tmp2lev(n)
             obsdat(njs(j)+nn) = tmp2dat(n)
             obserr(njs(j)+nn) = tmp2err(n)
             obsi(njs(j)+nn) = tmp2i(n)
             obsj(njs(j)+nn) = tmp2j(n)
             obsdep(njs(j)+nn) = tmp2dep(n)
             obshdxf(njs(j)+nn,:) = tmp2hdxf(n,:)
          END DO
          nobsgrd(i,j) = njs(j) + nn
       END DO
       IF(nn /= nj(j)) THEN
          !$OMP CRITICAL
          WRITE(6,'(A,2I)') 'OBS DATA SORT ERROR: ',nn,nj(j)
          WRITE(6,'(F6.2,A,F6.2)') j,'< J <',j+1
          WRITE(6,'(F6.2,A,F6.2)') MINVAL(tmp2j(njs(j)+1:njs(j)+nj(j))),'< OBSJ <',MAXVAL(tmp2j(njs(j)+1:njs(j)+nj(j)))
          !$OMP END CRITICAL
       END IF
    END DO
    !$OMP END DO
    !$OMP END PARALLEL

    IF(AOEI .and. OBS_SN)THEN
       WRITE(6,'(A)') 'Not Applying LETKF because of both AOER and OBS_SN'
       CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
       CALL finalize_mpi    
       STOP
    END IF

    !Ensemble spread
    IF(AOEI .or. OBS_SN) THEN
       DO n=1,nobs
          sprd_hdxf(n) = obshdxf(n,1)**2
          DO i=2,nbv
             sprd_hdxf(n) = sprd_hdxf(n) + obshdxf(n,i)**2
          END DO
          sprd_hdxf(n) = SQRT(sprd_hdxf(n) / REAL(nbv-1,r_size))
       END DO
    END IF

    !AOEI (Minamide and Zhang 2017; Monthly Weather Review) S.Ohishi 2019.08
    !Add ALPHA_AOEI S.Ohishi 2019.10
    IF(AOEI)THEN
       DO n=1,nobs

          obserrt = obserr(n) !specified obs. err
          dob2_sprd2 = obsdep(n)*obsdep(n)-sprd_hdxf(n)*sprd_hdxf(n) !dob**2 - sprd**2

          obserr(n) = sqrt( max(obserrt*obserrt,ALPHA_AOEI*dob2_sprd2) )

       END DO

    END IF

    !Obs. Error using S/N ratio (Carton et al. 2018; J. Clim.) S.Ohishi 2020.06
    IF(OBS_SN)THEN
       
       DO n=1,nobs
          SELECT CASE(NINT(obselm(n)))
             CASE(id_u_obs)
                obserr(n)=sprd_hdxf(n)*u_sn
             CASE(id_v_obs)
                obserr(n)=sprd_hdxf(n)*v_sn
             CASE(id_t_obs)
                IF(obslev(n) == 0.d0)THEN
                   obserr(n)=sprd_hdxf(n)*sst_sn
                ELSE
                   obserr(n)=sprd_hdxf(n)*t_sn
                END IF
             CASE(id_s_obs)
                IF(obslev(n) == 0.d0)THEN
                   obserr(n)=sprd_hdxf(n)*sss_sn
                ELSE
                   obserr(n)=sprd_hdxf(n)*s_sn
                END IF
             CASE(id_z_obs)
                obserr(n)=sprd_hdxf(n)*ssh_sn
             END SELECT
       END DO

    END IF

    DEALLOCATE( tmp2elm )
    DEALLOCATE( tmp2lon )
    DEALLOCATE( tmp2lat )
    DEALLOCATE( tmp2lev )
    DEALLOCATE( tmp2dat )
    DEALLOCATE( tmp2err )
    DEALLOCATE( tmp2i )
    DEALLOCATE( tmp2j )
    DEALLOCATE( tmp2dep )
    DEALLOCATE( tmp2hdxf )
    DEALLOCATE( tmpelm )
    DEALLOCATE( tmplon )
    DEALLOCATE( tmplat )
    DEALLOCATE( tmplev )
    DEALLOCATE( tmpdat )
    DEALLOCATE( tmperr )
    DEALLOCATE( tmpi )
    DEALLOCATE( tmpj )
    DEALLOCATE( tmpdep )
    DEALLOCATE( tmphdxf )
    DEALLOCATE( tmpqc )
    DEALLOCATE( sprd_hdxf ) !SO

    RETURN
  END SUBROUTINE set_letkf_obs
  !-----------------------------------------------------------------------
  ! Monitor departure from gues/anal mean |
  !----------------------------------------
  !
  ! Add Ts,Ti,Ss,Si monitor by S.Ohishi 2020.05 
  !
  !-----------------------------------------------------------------------
  SUBROUTINE monit_mean(file)
    CHARACTER(4),INTENT(IN) :: file
    REAL(r_size) :: v3d(nlon,nlat,nlev,nv3d)
    REAL(r_size) :: v2d(nlon,nlat,nv2d)
    REAL(r_size) :: elem
    REAL(r_size) :: bias_u,bias_v,bias_t,bias_s,bias_z
    REAL(r_size) :: bias_ts,bias_ti,bias_ss,bias_si !S.Ohishi 2020.05
    REAL(r_size) :: rmse_u,rmse_v,rmse_t,rmse_s,rmse_z
    REAL(r_size) :: rmse_ts,rmse_ti,rmse_ss,rmse_si !S.Ohishi 2020.05
    REAL(r_size) :: hdxf,dep,ri,rj,rk
    INTEGER :: n,iu,iv,it,is,iz
    INTEGER :: its,iti,iss,isi !S.Ohishi 2020.05
    CHARACTER(10) :: filename='filexxx.nc'

    rmse_u  = 0.0d0
    rmse_v  = 0.0d0
    rmse_t  = 0.0d0
    rmse_s = 0.0d0
    rmse_z = 0.0d0
    rmse_ts = 0.0d0 !S.Ohishi 2020.05
    rmse_ti = 0.0d0
    rmse_ss = 0.0d0
    rmse_si = 0.0d0
    bias_u = 0.0d0
    bias_v = 0.0d0
    bias_t = 0.0d0
    bias_s = 0.0d0
    bias_z = 0.0d0
    bias_ts = 0.0d0 !S.Ohishi 2020.05
    bias_ti = 0.0d0
    bias_ss = 0.0d0
    bias_si = 0.0d0
    iu = 0
    iv = 0
    it = 0
    is = 0
    iz = 0
    its = 0 !S.Ohishi 2020.05
    iti = 0
    iss = 0
    isi = 0

    WRITE(filename(1:7),'(A4,A3)') file,'_me'
    CALL read_grd(filename,v3d,v2d)

    DO n=1,nobs
       CALL Trans_XtoY(obselm(n),obsi(n),obsj(n),obslev(n),v3d,v2d,hdxf)
       dep = obsdat(n) - hdxf
       SELECT CASE(NINT(obselm(n)))
       CASE(id_u_obs)
          rmse_u = rmse_u + dep**2
          bias_u = bias_u + dep
          iu = iu + 1
       CASE(id_v_obs)
          rmse_v = rmse_v + dep**2
          bias_v = bias_v + dep
          iv = iv + 1
       CASE(id_t_obs)
          rmse_t = rmse_t + dep**2
          bias_t = bias_t + dep
          it = it + 1
          IF(obslev(n) == 0.d0)THEN
             rmse_ts = rmse_ts + dep**2
             bias_ts = bias_ts + dep
             its = its + 1             
          ELSE
             rmse_ti = rmse_ti + dep**2
             bias_ti = bias_ti + dep
             iti = iti + 1
          ENDIF
       CASE(id_s_obs)
          rmse_s = rmse_s + dep**2
          bias_s = bias_s + dep
          is = is + 1
          IF(obslev(n) == 0.d0)THEN
             rmse_ss = rmse_ss + dep**2
             bias_ss = bias_ss + dep
             iss = iss + 1
          ELSE
             rmse_si = rmse_si + dep**2
             bias_si = bias_si + dep
             isi = isi + 1
          END IF
       CASE(id_z_obs)
          rmse_z = rmse_z + dep**2
          bias_z = bias_z + dep
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
    IF(its == 0) THEN
       rmse_ts = undef
       bias_ts = undef
    ELSE
       rmse_ts = SQRT(rmse_ts / REAL(its,r_size))
       bias_ts = bias_ts / REAL(its,r_size)
    END IF
    IF(iti == 0) THEN
       rmse_ti = undef
       bias_ti = undef
    ELSE
       rmse_ti = SQRT(rmse_ti / REAL(iti,r_size))
       bias_ti = bias_ti / REAL(iti,r_size)
    END IF
    IF(is == 0) THEN
       rmse_s = undef
       bias_s = undef
    ELSE
       rmse_s = SQRT(rmse_s / REAL(is,r_size))
       bias_s = bias_s / REAL(is,r_size)
    END IF
    IF(iss == 0) THEN
       rmse_ss = undef
       bias_ss = undef
    ELSE
       rmse_ss = SQRT(rmse_ss / REAL(iss,r_size))
       bias_ss = bias_ss / REAL(iss,r_size)
    END IF
    IF(isi == 0) THEN
       rmse_si = undef
       bias_si = undef
    ELSE
       rmse_si = SQRT(rmse_si / REAL(isi,r_size))
       bias_si = bias_si / REAL(isi,r_size)
    END IF
    IF(iz == 0) THEN
       rmse_z = undef
       bias_z = undef
    ELSE
       rmse_z = SQRT(rmse_z / REAL(iz,r_size))
       bias_z = bias_z / REAL(iz,r_size)
    END IF

    WRITE(6,'(3A)') '== PARTIAL OBSERVATIONAL DEPARTURE (',file,') =================='
    WRITE(6,'(9A12)') 'U','V','T','Ts','Tint','S','Ss','Sint','ZETA'
    WRITE(6,'(9ES12.3)') bias_u,bias_v,bias_t,bias_ts,bias_ti,bias_s,bias_ss,bias_si,bias_z
    WRITE(6,'(9ES12.3)') rmse_u,rmse_v,rmse_t,rmse_ts,rmse_ti,rmse_s,rmse_ss,rmse_si,rmse_z
    WRITE(6,'(A)') '== NUMBER OF OBSERVATIONS =================================='
    WRITE(6,'(9A12)') 'U','V','T','Ts','Tint','S','Ss','Sint','ZETA'
    WRITE(6,'(9I12)') iu,iv,it,its,iti,is,iss,isi,iz
    WRITE(6,'(A)') '============================================================'

    RETURN
  END SUBROUTINE monit_mean

END MODULE letkf_obs
