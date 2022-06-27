MODULE common_letkf
!=======================================================================
!
! [PURPOSE:] Local Ensemble Transform Kalman Filtering (LETKF)
!            Model Independent Core Module
!
! [REFERENCES:]
!  [1] Ott et al., 2004: A local ensemble Kalman filter for atmospheric
!    data assimilation. Tellus, 56A, 415-428.
!  [2] Hunt et al., 2007: Efficient Data Assimilation for Spatiotemporal
!    Chaos: A Local Ensemble Transform Kalman Filter. Physica D, 230,
!    112-126.
!
! [HISTORY:]
!  01/21/2009 Takemasa Miyoshi  Created at U. of Maryland, College Park
!
!=======================================================================
!$USE OMP_LIB
  USE common_setting
  USE common
  USE common_mtx

  IMPLICIT NONE
  PUBLIC

CONTAINS
!=======================================================================
!  Main Subroutine of LETKF Core
!   INPUT
!     nobs             : array size, but only first nobsl elements are used
!     nobsl            : total number of observation assimilated at the point
!     hdxb(nobs,nbv)   : obs operator times fcst ens perturbations
!     rdiag(nobs)      : observation error variance
!     rloc(nobs)       : localization weighting function
!     dep(nobs)        : observation departure (yo-Hxb)
!     parm_infl        : covariance inflation parameter
!   OUTPUT
!     trans(nbv,nbv) : transformation matrix (If transm exists, large W component)
!   Optional,OUTPUT !S.Ohishi 2019.08
!     tansm(nbv)   : transformation matrix (small w component)
!     pao(nbv,nbv) : covariance matrix Pa
!-----------------------------------------------------------------------
! 2019.08 S.Ohishi modified refering Dr. Kotsuki code
!=======================================================================
SUBROUTINE letkf_core(nobs,nobsl,hdxb,rdiag,rloc,dep,parm_infl,trans,transm,pao) !SO
  IMPLICIT NONE
!  INTEGER,INTENT(IN) :: nbv
  INTEGER,INTENT(IN) :: nobs
  INTEGER,INTENT(IN) :: nobsl
  REAL(r_size),INTENT(IN) :: hdxb(1:nobs,1:nbv)
  REAL(r_size),INTENT(IN) :: rdiag(1:nobs)
  REAL(r_size),INTENT(IN) :: rloc(1:nobs)
  REAL(r_size),INTENT(IN) :: dep(1:nobs)
  REAL(r_size),INTENT(INOUT) :: parm_infl
  REAL(r_size),INTENT(OUT) :: trans(nbv,nbv)
  REAL(r_size) :: hdxb_rinv(nobsl,nbv) !S.Ohishi 2019.08
  REAL(r_size) :: eivec(nbv,nbv)
  REAL(r_size) :: eival(nbv)
  REAL(r_size) :: pa(nbv,nbv)
  REAL(r_size) :: work1(nbv,nbv)
  REAL(r_size) :: work2(nbv,nobsl) !S.Ohishi 2019.08
  REAL(r_size) :: work3(nbv)
  REAL(r_size) :: rho
  REAL(r_size) :: parm(4),sigma_o,gain
  REAL(r_size),PARAMETER :: sigma_b = 0.04d0 !error stdev of parm_infl
  INTEGER :: i,j,k
  !For RTPP & RTPS !S.Ohishi 2019.08
  REAL(r_size),INTENT(OUT),OPTIONAL :: transm(nbv)
  REAL(r_size),INTENT(OUT),OPTIONAL :: pao(nbv,nbv)

  IF(nobsl == 0) THEN

    trans = 0.0d0
    DO i=1,nbv
      trans(i,i) = SQRT(parm_infl)
    END DO

!   S.Ohishi 2019.08
    IF(PRESENT(transm))THEN
        transm=0.0d0
    END IF

!   S.Ohishi 2019.08
    IF(PRESENT(pao))THEN
        pao=0.0d0
        DO i=1,nbv
         pao(i,i)=parm_infl/REAL(nbv-1,r_size)
        END DO
    END IF

    RETURN
  ELSE
!-----------------------------------------------------------------------
!  hdxb Rinv
!-----------------------------------------------------------------------
  DO j=1,nbv
    DO i=1,nobsl
      hdxb_rinv(i,j) = hdxb(i,j) / rdiag(i) * rloc(i)
    END DO
  END DO
!-----------------------------------------------------------------------
!  hdxb^T Rinv hdxb
!-----------------------------------------------------------------------
  CALL dgemm('t','n',nbv,nbv,nobsl, &
        & 1.0d0,hdxb_rinv,nobsl,hdxb(1:nobsl,:), &
        & nobsl,0.0d0,work1,nbv) !S.Ohishi 2019.08

!  DO j=1,nbv
!    DO i=1,nbv
!      work1(i,j) = hdxb_rinv(1,i) * hdxb(1,j)
!      DO k=2,nobsl
!        work1(i,j) = work1(i,j) + hdxb_rinv(k,i) * hdxb(k,j)
!      END DO
!    END DO
!  END DO

!-----------------------------------------------------------------------
!  hdxb^T Rinv hdxb + (m-1) I / rho (covariance inflation)
!-----------------------------------------------------------------------
  rho = 1.0d0 / parm_infl
  DO i=1,nbv
    work1(i,i) = work1(i,i) + REAL(nbv-1,r_size) * rho
  END DO
!-----------------------------------------------------------------------
!  eigenvalues and eigenvectors of [ hdxb^T Rinv hdxb + (m-1) I ]
!-----------------------------------------------------------------------
  CALL mtx_eigen(1,nbv,work1,eival,eivec,i)
!-----------------------------------------------------------------------
!  Pa = [ hdxb^T Rinv hdxb + (m-1) I ]inv
!-----------------------------------------------------------------------
  DO j=1,nbv
    DO i=1,nbv
      work1(i,j) = eivec(i,j) / eival(j)
    END DO
  END DO

  CALL dgemm('n','t',nbv,nbv,nbv, &
        & 1.0d0,work1,nbv,eivec,&
        & nbv,0.0d0,pa,nbv) !S.Ohishi 2019.08

!  DO j=1,nbv
!    DO i=1,nbv
!      pa(i,j) = work1(i,1) * eivec(j,1)
!      DO k=2,nbv
!        pa(i,j) = pa(i,j) + work1(i,k) * eivec(j,k)
!      END DO
!    END DO
!  END DO

!-----------------------------------------------------------------------
!  Pa hdxb_rinv^T
!-----------------------------------------------------------------------

  CALL dgemm('n','t',nbv,nobsl,nbv, &
        & 1.0d0,pa,nbv,hdxb_rinv,&
        & nobsl,0.0d0,work2,nbv) !S.Ohishi 2019.08

!  DO j=1,nobsl
!    DO i=1,nbv
!      work2(i,j) = pa(i,1) * hdxb_rinv(j,1)
!      DO k=2,nbv
!        work2(i,j) = work2(i,j) + pa(i,k) * hdxb_rinv(j,k)
!      END DO
!    END DO
!  END DO

!-----------------------------------------------------------------------
!  Pa hdxb_rinv^T dep
!-----------------------------------------------------------------------
  DO i=1,nbv
    work3(i) = work2(i,1) * dep(1)
    DO j=2,nobsl
      work3(i) = work3(i) + work2(i,j) * dep(j)
    END DO
  END DO
!-----------------------------------------------------------------------
!  T = sqrt[(m-1)Pa]
!-----------------------------------------------------------------------

  DO j=1,nbv
    rho = SQRT( REAL(nbv-1,r_size) / eival(j) )
    DO i=1,nbv
      work1(i,j) = eivec(i,j) * rho
    END DO
  END DO

  CALL dgemm('n','t',nbv,nbv,nbv, &
        & 1.0d0,work1,nbv,eivec,&
        & nbv,0.0d0,trans,nbv)

!  DO j=1,nbv
!    DO i=1,nbv
!      trans(i,j) = work1(i,1) * eivec(j,1)
!      DO k=2,nbv
!        trans(i,j) = trans(i,j) + work1(i,k) * eivec(j,k)
!      END DO
!    END DO
!  END DO

!-----------------------------------------------------------------------
!  T + Pa hdxb_rinv^T dep
!-----------------------------------------------------------------------

!S.Ohishi 2019.08
!If transm is present, return both trans and transm without adding them

  IF(PRESENT(transm))THEN
    transm=work3
  ELSE
    DO j=1,nbv
        DO i=1,nbv
        trans(i,j) = trans(i,j) + work3(i)
        END DO
    END DO
  END IF

!S.Ohishi 2019.08
  IF(PRESENT(pao)) pao = pa

!-----------------------------------------------------------------------
!  Inflation estimation
!-----------------------------------------------------------------------
  parm = 0.0d0
  DO i=1,nobsl
    parm(1) = parm(1) + dep(i)*dep(i)/rdiag(i) * rloc(i)
  END DO
  DO j=1,nbv
    DO i=1,nobsl
      parm(2) = parm(2) + hdxb_rinv(i,j) * hdxb(i,j)
    END DO
  END DO
  parm(2) = parm(2) / REAL(nbv-1,r_size)
! 2015/11/01 (masked grid: no sea)
  if(abs(parm(2)).lt.1e-10) parm(2) = 1.0
!
  parm(3) = SUM(rloc(1:nobsl))
  parm(4) = (parm(1)-parm(3))/parm(2) - parm_infl
!  sigma_o = 1.0d0/REAL(nobsl,r_size)/MAXVAL(rloc(1:nobsl))
  sigma_o = 2.0d0/parm(3)*((parm_infl*parm(2)+parm(3))/parm(2))**2
  gain = sigma_b**2 / (sigma_o + sigma_b**2)
  parm_infl = parm_infl + gain * parm(4)

! debug
!  write(6,*) 'nbv nobsl ', nbv, nobsl
!  write(6,*) 'parm(1) ', parm(1)
!  write(6,*) 'parm(2) ', parm(2)
!  write(6,*) 'parm(3) ', parm(3)
!  write(6,*) 'parm(4) ', parm(4)
!  write(6,*) 'parm_infl ', parm_infl
!  write(6,*) 'sigma_o gain ', sigma_o, gain
!  call flush(6)

  RETURN
  END IF
END SUBROUTINE letkf_core

END MODULE common_letkf
