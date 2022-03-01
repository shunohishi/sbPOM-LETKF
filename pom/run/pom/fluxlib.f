!_______________________________________________________________________
      subroutine surface_airseaflux
! set time dependent surface airsea fluxes
!     upward(sea->air) heat transport: positive wtsurf and swrad
!     upward(sea->air) salt transport: positive wssurf
!
      implicit none
      include 'pom.h'
      integer timing,i,j
      real time0_data,ratio
      real qh(im,jm),qe(im,jm),lwrad(im,jm),sstm(im,jm)
      real dqds, adj1
      parameter(dqds=10./(30.*86400.),adj1=1.0)

      integer iglb,jglb,iloc,jloc,nproc

      if(iint.eq.1) call read_atm_netcdf            

      time0_data=atmtime_julday-julday_start
      call readtiming(timing,ratio,atmtime_dayint,time0_data)
      if(timing.eq.1) then
        call read_atm_netcdf
        if(my_task.eq.master_task) then
          write(6,*) 'atm timing ratio time ',timing,ratio
     $               ,dti*float(iint-1)/86400.e0+time0
        end if
      end if

      do j=1,jm
        do i=1,im
          windu(i,j)=(1.-ratio)*windu0(i,j)+ratio*windu1(i,j)
          windv(i,j)=(1.-ratio)*windv0(i,j)+ratio*windv1(i,j)
          airt(i,j)=(1.-ratio)*airt0(i,j)+ratio*airt1(i,j)
          airh(i,j)=(1.-ratio)*airh0(i,j)+ratio*airh1(i,j)
          swrad(i,j)=(1.-ratio)*swrad0(i,j)+ratio*swrad1(i,j)
          cloud(i,j)=(1.-ratio)*cloud0(i,j)+ratio*cloud1(i,j)
          sst(i,j)=(1.-ratio)*sst0(i,j)+ratio*sst1(i,j)
        end do
      end do       

! surface heat flux

      do j=1,jm
        do i=1,im
 
! [Celcius] -> {Kelvin]
          sstm(i,j)=t(i,j,1)+tbias+273.15
          airt(i,j)=airt(i,j)+273.15
        end do
      end do

      call wind_stress(wusurf,wvsurf,windu,windv,im,jm)

      call sensibleh(qh,sstm,airt,windu,windv,im,jm,my_task)
      call latenth(qe,sstm,airt,airh,windu,windv,im,jm,my_task)
      call longrad(lwrad,sstm,airt,airh,cloud,im,jm,my_task)
 
! unit change [W/m**2] -> [K m/sec]

      do j=2,jmm1
        do i=2,imm1
          wtsurf(i,j)=
     $     (adj1*(qh(i,j)+qe(i,j))+lwrad(i,j))/(4.1876e6)*fsm(i,j)
          swrad(i,j)=-swrad(i,j)/(4.1876e6)*fsm(i,j)
        end do
      end do

! for output
!      call exchange2d_mpi(wtsurf,im,jm)
!      call exchange2d_mpi(swrad,im,jm)

! debug
!      iglb=(im_global-2)/2
!      jglb=(jm_global-2)/2
!      nproc=int((iglb-2)/(im-2))+int((jglb-2)/(jm-2))*nproc_x
!      if(my_task.eq.nproc .and. mod(iint-1,200).eq.0) then
!        iloc=iglb-mod(nproc,nproc_x)*(im-2)
!        jloc=jglb-(nproc/nproc_x)*(jm-2)
!        write(99,*) 'in surface_airseaflux iint ',iint
!        write(99,*) 'iloc jloc nproc ',iloc,jloc,nproc
!        write(99,*) 'iglb jglb ',iglb,jglb
!        write(99,*) 'qh qe ',qh(iloc,jloc),qe(iloc,jloc)
!        write(99,*) 'swrad lwrad '
!     $    ,swrad(iloc,jloc)*4.1876e6,lwrad(iloc,jloc)
!        write(99,*) 'wtsurf [W/m**2] ',wtsurf(iloc,jloc)*4.1876e6
!        write(99,*) 'wtsurf ',wtsurf(iloc,jloc)
!        call flush(99) 
!      end if
 
! surface salt flux [psu m/sec]

      do j=2,jmm1
        do i=2,imm1
          wssurf(i,j)=dqds*(s(i,j,1)-sclimm(i,j,1))*fsm(i,j)      
        end do
      end do

      return
      end

!_______________________________________________________________________
      subroutine wind_stress(wusurf,wvsurf,windu,windv,im,jm)
! set wind stress  (Large and Pond, 1981)
      implicit none
! -- arguments
      integer im,jm
      real wusurf(im,jm),wvsurf(im,jm)
      real windu(im,jm),windv(im,jm)
! -- local      
      real cofvmag(3)
      data cofvmag/4.0,10.0,26.0/
      integer i,j
      real wvel,cd,rhoa
!
      rhoa=1.2
      do j=2,jm-1
        do i=2,im-1
          wvel=sqrt(windu(i,j)*windu(i,j)+windv(i,j)*windv(i,j))
          if(wvel.lt.cofvmag(2)) then
            cd=1.14
          else if( wvel.lt.cofvmag(3) ) then
            cd=0.49+0.065*wvel
          else
            cd=0.49+0.065*cofvmag(3)
          end if
          cd=cd*1.e-3
          wusurf(i,j)=rhoa*cd*wvel*windu(i,j)
          wvsurf(i,j)=rhoa*cd*wvel*windv(i,j)
! unit and direction changes
          wusurf(i,j)=-(wusurf(i,j))*1.e-4*0
          wvsurf(i,j)=-(wvsurf(i,j))*1.e-4*0
        end do
      end do
!
      return
      end

!_______________________________________________________________________
      subroutine sensibleh(qh,ts,ta,ua,va,im,jm,my_task)
      implicit none
!
!     upward(sea->air): positive
!     ts: sea surface temperature [K]
!     ta: air temperature above 2m [k]
!
! -- arugumanets
      integer im,jm,my_task
      real qh(im,jm)
      real ts(im,jm),ta(im,jm),ua(im,jm),va(im,jm)
      real ua2m(im,jm),va2m(im,jm),w(im,jm)
! -- local
      integer i,j
      real z0,z10,z2,rkalm
      parameter(z10=10.0,z2=1.5,rkalm=0.4)
      real rhoa,cp,Chc
      parameter(rhoa=1.2,cp=1.005e03,Chc=1.1e-3)
!      real dqh(im,jm)
      real wmag(im,jm)
      real Ch10m(im,jm),Ch2m(im,jm)
      real Cd10m(im,jm),Cd2m(im,jm)

! debug
      integer iglb,jglb,iloc,jloc,nproc
      integer im_global,jm_global,nproc_x
      im_global=254
      jm_global=182
      nproc_x=(im_global-2)/(im-2)
      iglb=(im_global-2)/2
      jglb=(jm_global-2)/2
      nproc=int((iglb-2)/(im-2))+int((jglb-2)/(jm-2))*nproc_x
      iloc=iglb-mod(nproc,nproc_x)*(im-2)
      jloc=jglb-(nproc/nproc_x)*(jm-2)

      call bulkcof(cd10m,ua,va,ts,ta,wmag,im,jm,1)
      call bulkcof(ch10m,ua,va,ts,ta,wmag,im,jm,2)

!      transformation from Ch(10m) to Ch(2m)

      do j=2,jm-1
         do i=2,im-1

!      transformation from Ch(10m) to Ch(2m)

             if( cd10m(i,j).gt.0. ) then
               Cd2m(i,j) = rkalm**2
     &           *( rkalm*Cd10m(i,j)**(-0.5)-log(z10/z2) )**(-2)
             else
               Cd2m(i,j) = 0.0
             end if
             if( ch10m(i,j).gt.0. ) then
               Ch2m(i,j) = rkalm*sqrt(Cd2m(i,j))
     &            /( rkalm*sqrt( Cd10m(i,j) )/Ch10m(i,j) +log(z2/z10) )
             else
               Ch2m(i,j) = 0.0
             end if

!      transformation from U(10m) to U(2m)

             if( cd10m(i,j).gt.0. ) then
               z0 = exp( log(z10)-rkalm*( Cd10m(i,j) )**(-0.5) )
             else
               z0 = 0.0
             end if
             if( z0.gt.1.e-7 ) then
               ua2m(i,j) = ua(i,j)*log(z2/z0)/log(z10/z0)
               va2m(i,j) = va(i,j)*log(z2/z0)/log(z10/z0)
             else
               ua2m(i,j) = 0.0
               va2m(i,j) = 0.0
             end if

         enddo
      enddo


      do j=2,jm-1
         do i=2,im-1
!  Kondo(1975)
            qh(i,j) = cp*rhoa*Ch2m(i,j)
     &      *sqrt( ua2m(i,j)**2+va2m(i,j)**2 )*( ts(i,j)-ta(i,j) ) 

! debug
!            if(my_task.eq.nproc .and. i.eq.iloc .and. j.eq. jloc) then
!              write(99,*) 'in sensibleh '
!              write(99,*) 'iloc jloc nproc ',iloc,jloc,nproc
!              write(99,*) 'iglb jglb ',iglb,jglb
!              write(99,*) 'qh ',qh(iloc,jloc)
!              write(99,*) 'Ch2m ',Ch2m(iloc,jloc)
!              write(99,*) 'ua2m va2m ',ua2m(iloc,jloc),va2m(iloc,jloc)
!              write(99,*) 'ts ta ',ts(iloc,jloc),ta(iloc,jloc)
!              call flush(99)
!            end if

!            dqh(i,j) = cp*rhoa*Ch2m(i,j)
!     &      *sqrt( ua2m(i,j)**2+va2m(i,j)**2 )
             if(ua2m(i,j).eq.0.0.and.va2m(i,j).eq.0.) qh(i,j)=0.0
!             if(ua2m(i,j).eq.0.0.and.va2m(i,j).eq.0.) dqh(i,j)=0.0
         enddo
      enddo
c
      return
      end

!_______________________________________________________________________
      subroutine latenth(qe,ts,ta,qa,ua,va,im,jm,my_task)
      implicit none
!
!     upward(sea->air): positive
!     ts: sea surface temperature [K]
!     ta: air temperature above 2m [k]
!     qa: specific humidity [g/g]
!     ua,va: wind velocity above 10m
!
! -- argumants
      integer im,jm,my_task
      real Qe(im,jm)
      real ts(im,jm),ta(im,jm),qa(im,jm),ua(im,jm),va(im,jm)
! -- local
      integer i,j
      real z10,z2,rkalm
      parameter(z10=10.0,z2=1.5,rkalm=0.4)
      real rhoa,rL,Pa,Ce
      parameter(rhoa=1.2,rL=2.501e6,Pa=1013.,Ce=1.1e-3)
      real ua2m(im,jm),va2m(im,jm),wmag(im,jm)
      real Ce10m(im,jm),Ce2m(im,jm)
      real Cd10m(im,jm),Cd2m(im,jm)
!      real dQe(im,jm)
      real devide,qs,z0,esat,desat_dt,esat_const
    
! debug
      integer iglb,jglb,iloc,jloc,nproc
      integer im_global,jm_global,nproc_x
      im_global=254
      jm_global=182
      nproc_x=(im_global-2)/(im-2)
      iglb=(im_global-2)/2
      jglb=(jm_global-2)/2
      nproc=int((iglb-2)/(im-2))+int((jglb-2)/(jm-2))*nproc_x
      iloc=iglb-mod(nproc,nproc_x)*(im-2)
      jloc=jglb-(nproc/nproc_x)*(jm-2)

!     bulk coefficients

      call bulkcof(cd10m,ua,va,ts,ta,wmag,im,jm,1)
      call bulkcof(ce10m,ua,va,ts,ta,wmag,im,jm,3)

      do j=2,jm-1
         do i=2,im-1

!      transformation from Ce(10m) to Ce(2m)

             if( cd10m(i,j).gt.0. ) then
               Cd2m(i,j) = rkalm**2
     &           *( rkalm*Cd10m(i,j)**(-0.5)-log(z10/z2) )**(-2)
             else
               Cd2m(i,j) = 0.0
             end if
             if( ce10m(i,j).gt.0. ) then
               devide = rkalm*sqrt( Cd10m(i,j) )/Ce10m(i,j) +log(z2/z10) 
             else
               devide = 0.
             end if
             if( devide.gt.0. ) then
               Ce2m(i,j) = rkalm*sqrt(Cd2m(i,j))/devide
             else
               Ce2m(i,j) = 0.0
             end if

!      transformation from U(10m) to U(2m)

             if( cd10m(i,j).gt.0. ) then
               z0 = exp( log(z10)-rkalm*( Cd10m(i,j) )**(-0.5) )
               if( z0.gt.1.e-7 ) then
                 devide = log(z10/z0)
               else
                 z0 = 0.0
                 devide = 0.
               end if             
             else
               z0 = 0.0
               devide = 0.
             end if
             if( devide.gt.1.e-7 .and. z0.gt.1.e-7 ) then
               ua2m(i,j) = ua(i,j)*log(z2/z0)/devide
               va2m(i,j) = va(i,j)*log(z2/z0)/devide
             else
               Ce2m(i,j) = 0.0
               ua2m(i,j) = 0.0
               va2m(i,j) = 0.0
             end if

         end do
      end do

      do j=2,jm-1
         do i=2,im-1

!               qs = 0.622*esat(ts(i,j))/Pa
!     &             /(1.-0.378*esat(ts(i,j))/Pa)
               esat_const = 
     &          6.1078 * 10.0 ** (7.5 * (ts(i,j) - 273.15) /
     &                            (237.3 + ts(i,j) - 273.15))
               qs = 0.622*esat_const/Pa
     &             /(1.-0.378*esat_const/Pa)
               Qe(i,j) = rL*rhoa*Ce2m(i,j)
     &              *sqrt( ua2m(i,j)**2+va2m(i,j)**2 )
     &              *( qs-qa(i,j) )

! debug
!            if(my_task.eq.nproc .and. i.eq.iloc .and. j.eq. jloc) then
!              write(99,*) 'in latenth '
!              write(99,*) 'iloc jloc nproc ',iloc,jloc,nproc
!              write(99,*) 'iglb jglb ',iglb,jglb
!              write(99,*) 'qe ',qe(iloc,jloc)
!              write(99,*) 'Ce2m ',Ce2m(iloc,jloc)
!              write(99,*) 'ua2m va2m ',ua2m(iloc,jloc),va2m(iloc,jloc)
!              write(99,*) 'qs qa ts ',qs,qa(iloc,jloc),ts(iloc,jloc)
!              call flush(99)
!            end if

!               dQe(i,j) = rL*rhoa*Ce2m(i,j)
!     &              *sqrt( ua2m(i,j)**2+va2m(i,j)**2 )
!     &              * desat_dt( ts(i,j) )*0.622/Pa
         enddo
      enddo

      return
      end

!_______________________________________________________________________
!      real function esat(tk)
!
!      parameter(a=7.5, b=237.3)
!      real tk,tc,a,b
!
!      tc = tk - 273.15
!      esat = 6.1078*10.0**( a*tc/(b+tc) )
!
!      return
!      end
!
!!_______________________________________________________________________
!      real function desat_dt(tk)
!
!      parameter(a=7.5, b=237.3)
!      real tk,tc,a,b
!
!      tc = tk - 273.15
!      desat_dt = 6.1078*( 2500.- 2.4*tc ) / ( 0.4615*tk**2 )
!     &     * 10.0**( a*tc/(b+tc) )
!
!      return
!      end

!_______________________________________________________________________
      subroutine bulkcof(cf,u,v,ts,ta,w,id,jd,ind)
!
***********************************************************************
*          Bulk coefficient applied to turbulent fluxes               *
*                        Kondo(1975)                                  *
*
*          IND=1 ----> CDD   for momentum flux
*          IND=2 ----> CHD   for sensible heat flux
*          IND=3 ----> CED   for latent heat flux
*
*          coefficient        in the nutral case
*          CDD(ID,JD,MD)      CD(ID,JD,MD)
*          CHD(ID,JD,MD)      CH(ID,JD,MD)
*          CED(ID,JD,MD)      CE(ID,JD,MD)
*
*     Input data      Coefficient which depends on the wind speed
*                     in the nutral case  
*
*                     CD(ID,JD,MD)
*     W(ID,JD,MD)---->CH(ID,JD,MD)
*                     CE(ID,JD,MD) 
*
*               the stability conditions  stability parameter 
*
* TS(ID,JD,MD)    (TS-TA)<0      stable   S<0   S=S0*(|S0|/(|S0|+0.01))
* TA(ID,JD,MD)--->(TS-TA)>0    unstable   S>0   S0=(TS-TA)/W**2
*                 (TS-TA)=0      nutral   S=0   If W=0, S=0
*
************************************************************************
C
      DIMENSION W(ID,JD),TS(ID,JD),TA(ID,JD)
      DIMENSION CF(ID,JD), U(ID,JD),V(ID,JD)
      DIMENSION A(3,5),B(3,5),C(3,5),P(3,5),RR(3)
C
C---------------Coefficients of (CD, CH, CE) formula 
C
      DATA ( (A(I,J),I=1,3), J=1,5 ) 
     &     /0., 0., 0., 0.771, 0.927, 0.969, 0.867, 1.15, 1.18,
     &      1.2, 1.17, 1.196, 0., 1.652, 1.68/
      DATA ( (B(I,J),I=1,3), J=1,5 )
     &  /1.08, 1.185, 1.23, 0.0858, 0.0546, 0.0521, 0.0667, 0.01, 0.01, 
     &      0.025, 0.0075, 0.008, 0.073, -0.017, -0.016/
      DATA ( (C(I,J),I=1,3), J=1,5 )
     &     /0., 0., 0., 0., 0., 0., 0., 0., 0., 
     &      -0.00045, -0.0004, 0., 0., 0., 0./
      DATA ( (P(I,J),I=1,3), J=1,5 )
     &     /-0.15, -0.157, -0.16, 1., 1., 1., 1., 1., 1., 
     &      1., 1., 1., 1., 1., 1./
      DATA (RR(I),I=1,3) /0.47, 0.63, 0.63/
C
C--------------------Write  coefficients
C
c      WRITE(6,*) 'Coefficients in the nutral case calculations' 
c      DO 1000 I=1, 3
c      WRITE(6,100) (A(I,J),J=1,5)
c      WRITE(6,100) (B(I,J),J=1,5)
c      WRITE(6,100) (C(I,J),J=1,5)
c      WRITE(6,100) (P(I,J),J=1,5)
c 1000 CONTINUE
c  100 FORMAT(1H , 5(F8.5,1x)/)
C

      DO 2000 I=1, ID
      DO 2000 J=1, JD

       W(I,J) = sqrt( U(I,J)**2 + V(I,J)**2 )

c       IF(DMSK(I,J). NE. 0.) THEN
C
C-----------------CF0 for the nutral case
C
        IF(W(I,J).LE.2.2) THEN
         KKK=1
        ELSE
         IF(W(I,J).LE.5.) THEN
          KKK=2
         ELSE
          IF(W(I,J).LE.8.) THEN
           KKK=3
          ELSE
           IF(W(I,J).LE.25.) THEN
            KKK=4
           ELSE
            KKK=5
           END IF
          END IF
         END IF
        END IF
C 
C----------------The calculation is possible for W lager than 0.3 m/s, 
C                based on Kondo(1975) formula.
C                MSK cannot be applied for all valiables.
C
      IF(W(I,J). LT. 0.3) THEN
c       CF(I,J)=-9999.
       CF(I,J)=1.1e-3
       GO TO 2000
      END IF  
C
      CF0=1.E-3*( A(IND,KKK)+B(IND,KKK)*(W(I,J)**P(IND,KKK))+
     &            C(IND,KKK)*( W(I,J)-8. )**2 )
C
C------------------S for stability parameter
C
      IF(W(I,J). NE. 0.) THEN
         S0=( TS(I,J)-TA(I,J) )/(W(I,J)**2)
         S=S0*( ABS(S0)/(ABS(S0)+0.01) )
      ELSE
         S=0.
      END IF
C
C-----------------CF(I,J,M) for bulk coefficient
C
C-----------------for the stable case
C
       IF(S. LT. 0.) THEN
        IF(S. GE. -3.3) THEN
         CF(I,J)=CF0*( 0.1+0.03*S+0.9*EXP(4.8*S) )
        ELSE
         CF(I,J)=0.
        END IF
C
C-----------------for the unstable case
C
       ELSE
        IF(S. GT. 0.) THEN
         CF(I,J)=CF0*( 1.0+RR(IND)*S**0.5 )
C
C-----------------for the nutral case
C
        ELSE
         CF(I,J)=CF0
        END IF
       END IF
C
C----------------for the lack value
C
c       ELSE
c        CF(I,J)=-9999.
c       END IF
C
 2000 CONTINUE

C  
      RETURN
      END

!_______________________________________________________________________
      subroutine longrad(ql,ts,ta,qa,tcdc,im,jm,my_task)
      implicit none
!
!     upward(sea->air): positive
!     ts: sea surface temperature [K]
!     ta: air temperature above 2m [k]
!     qa: specific humidity [g/g]
!     tcdc: total clound cover [0:1]
!     downward: positive
!
! -- arguments
      integer im,jm,my_task
      real ql(im,jm),ts(im,jm),ta(im,jm),qa(im,jm),tcdc(im,jm)
! -- local
      integer i,j
!      real dql(im,jm)
      real ea
      real eps,sigma,bb,Pa

! debug
      integer iglb,jglb,iloc,jloc,nproc
      integer im_global,jm_global,nproc_x
      im_global=254
      jm_global=182
      nproc_x=(im_global-2)/(im-2)
      iglb=(im_global-2)/2
      jglb=(jm_global-2)/2
      nproc=int((iglb-2)/(im-2))+int((jglb-2)/(jm-2))*nproc_x
      iloc=iglb-mod(nproc,nproc_x)*(im-2)
      jloc=jglb-(nproc/nproc_x)*(jm-2)

      eps = 0.97        ! emissivity of the ocean     
      sigma = 5.670e-8  ! Stefan-Boltzmann constant [W m^-2 K^-4]
      bb = 0.8          ! linear correction factor
      Pa = 1013.        ! atmospheric pressure on sea surface

      do j=2,jm-1
         do i=2,im-1
            ea = qa(i,j)*Pa/(0.622+0.378*qa(i,j))
            ql(i,j) = eps*sigma*ts(i,j)**4
     &           *( 0.39 - 0.05*( ea )**0.5 )
     &           *( 1.0 - bb*tcdc(i,j) )
     &           + 4.0*eps*sigma*ts(i,j)**3*(ts(i,j)-ta(i,j))

! debug
!            if(my_task.eq.nproc .and. i.eq.iloc .and. j.eq. jloc) then
!              write(99,*) 'in longrad '
!              write(99,*) 'iloc jloc nproc ',iloc,jloc,nproc
!              write(99,*) 'iglb jglb ',iglb,jglb
!              write(99,*) 'ql ',ql(iloc,jloc)
!              write(99,*) 'ts ta ',ts(iloc,jloc),ta(iloc,jloc)
!              write(99,*) 'ea qa tcdc ',ea,qa(iloc,jloc),tcdc(iloc,jloc)
!              call flush(99)
!            end if

!            dql(i,j) = 4.0*eps*sigma*ts(i,j)**3
!     &           *( 0.39 - 0.05*( ea )**0.5 )
!     &           *( 1.0 - bb*tcdc(i,j) )
!     &           + 4.0*eps*sigma*ts(i,j)*3
!     &           + 3.0*4.0*eps*sigma*ts(i,j)**2*(ts(i,j)-ta(i,j))
         end do
      end do

      return
      end

