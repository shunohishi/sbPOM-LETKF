!
! Module for definition of constants and estimation of parameters
! for Mellor-Yamada-Nakanishi-Nino-Furuichi versions of level 2.5 turbulence
! clousure models. Coded by vsm, November 2012
!

! Define GEN_SIGMA variable for the generalized sigma coordinates model,
! else regular sigma coordinates model is considered.
! Different is the vertical model layers layout description used.

!#define GEN_SIGMA

      module MYNNF_lev25_2012

       implicit none

       private
!
!     Fixed parameters for MYNNF-2.5 scheme 
!
        real,parameter::
     *  g1 = 0.235,
     *  b1 = 24.0,
     *  b2 = 15.0,
     *  c2 = 0.75,
     *  c3 = 0.352,
     *  c4 = 0.0,
     *  c5 = 0.2,
     *  pr = 0.74,
     *  vk = 0.4,                              !Karman const
     *  a1 = b1*( 1.0-3.0*g1 )/6.0,
!     *  c1 = g1 -1.0/( 3.0*a1*b1**(1.0/3.0) ),  !on SX only integer power is allowed
     *  c1 = 0.13706763,
     *  a2 = a1*( g1-c1 )/( g1*pr ),
     *  g2 = b2/b1*( 1.0-c3 ) +2.0*a1/b1*( 3.0-2.0*c2 )

      real,parameter:: alp1 = 0.23    !Multiplier for L as planetary boundary layer z-scale
      real,parameter:: alp2 = 0.53    !0.53 in Furuichi et al, 2012 for LES; 1 in NN
      real,parameter:: alp3 = 1./3.7  !1/3.7 in NN; 1 in Furuichi et al, 2012, but no 
      real,parameter:: almost_zero = 1e-12 !small value for nonsingular numerics
!
!     Constants for surfase boundary layer Monin-Obukhov length scale "Lmo" estimation
!
      real,parameter:: grav = 9.808
      real,parameter:: alpha_sw = 1.7e-4  !1/K, sea water thermal expansion coefficient
      real,parameter:: beta_sw = 7.5e-4   !1/PSU, sea water thermal expansion coefficient
!
!     Traditional values for Mellor-Yamada level 2.5 model parameters:
!
      real,parameter::
     *  a1_my= 0.92,
     *  a2_my= 0.74,
     *  b1_my=16.60,
     *  b2_my=10.10,
     *  c1_my= 0.08

      public mynn_get_l,mynn_get_q2l2,mynn_get_ShSm,my_get_ShSm

      contains

      subroutine mynn_get_l(l,q2,boygr,wtsurf,wssurf,
     *                      wusurf_t,wvsurf_t,z0,
#ifdef GEN_SIGMA
     *                      z,dzz,im,jm,kb,ntp)
#else
     *                      z,dzz,dh,im,jm,kb,ntp)
#endif
!
!     Turbulence length scale estimation for MYNN model
!
      implicit none
      integer,intent(in):: im,jm,kb,ntp
      real,intent(out):: l(im,jm,kb)        !turbulence length scale, m 
      real,intent(in)::  q2(im,jm,kb)       !square of turbulence velocity scale
      real,intent(in)::  boygr(im,jm,kb)    !in POM, N**2 ~= -boygr/1.025
      real,intent(in)::  wtsurf(im,jm)      !Total_Heat_Flux/(pho*cp), K*(m/s), positive if ocean losses heat
      real,intent(in)::  wssurf(im,jm)      !Total_Salt_Flux/Rho_fresh_water, [PSU*(m/s)], -(s+sbias)*WQ*RoFWR where
!                                            WQ [kg/m**2/s] is positive mass flux for evaporation (salination),
!                                            wssurf is positive for precipitation case (desalination)
      real,intent(in):: wusurf_t(im,jm),wvsurf_t(im,jm)   !surface_stress/Rho, defined in T-points (not in U and V points as wusrf and wvsurf of POM).
      real,intent(in):: z0(im,jm)           !surface roughness length; bottom roughness is fixed as 0.01 m
#ifdef GEN_SIGMA
      real,intent(in):: z(im,jm,kb),dzz(im,jm,kb)
#else
      real,intent(in):: z(kb),dzz(kb),dh(im,jm)     !sea depth dh is required only for SIGMA layers model
#endif
!
!     Working arrays used for PBL scale,
!     inverse Monin-Obukhov scale and buoyancy flux at the surface
!
      real lt(im,jm),LmoR(im,jm),Bf(im,jm)
!
      integer i,j,k
      real qdz,zk,lb,ls,lh,lr,u_star,zn_MO,N,q,qc,rr,hf
!
!     Irradiance parameters after Paulson and Simpson, JPO, 1977, 952-956.
!     Same as in the PROFT subroutine. Theoretically for the surface buoyancy flux
!     estimation have to use only part of SWR adsorbed by ocean surface mixed layer.
!     Use ocean depth as mixed layer depth. Part of radiation reaching the bottom
!     is considered to be removed from system. Could be important to prevent over-heating
!     (proft) and over-stabilization near the surface in the shallow coastal zones.
!
!     Irradiance parameters after Paulson and Simpson (1977)
!                 NTP         =     1      2      3      4     5
!             JERLOV TYPE     =     I      IA     IB     II    III
      real,parameter:: R(5)   = (/  0.58,  0.62,  0.67,  0.77, 0.78/)
      real,parameter:: AD1(5) = (/  0.35,  0.60,  1.00,  1.50, 1.40/)
      real,parameter:: AD2(5) = (/ 23.00, 20.00, 17.00, 14.00, 7.90/)
!
!     Use lt and bf as working arrays for vertical integrals
      do j=1,jm
        do i=1,im
          lt(i,j)  = 0.
          bf(i,j) = 0.
        end do
      end do
c
      do k=2,kb-1
        do j=1,jm
          do i=1,im
#ifdef GEN_SIGMA
            qdz=sqrt(abs(q2(i,j,k)))*dzz(i,j,k-1)
            lt(i,j)=lt(i,j)+qdz*abs(z(i,j,k))
#else
            qdz=sqrt(abs(q2(i,j,k)))*dzz(k-1)*dh(i,j)
            zk=abs(z(k)*dh(i,j))
            lt(i,j)=lt(i,j)+qdz*zk
#endif
            bf(i,j)=bf(i,j)+qdz
          enddo
        enddo
      enddo
c
      do j=1,jm
        do i=1,im
c ** Length scale depending on the PBL depth
          lt(i,j)=alp1*
     *            max(lt(i,j)/max(bf(i,j),almost_zero),almost_zero)
          l(i,j,1)=vk*abs(z0(i,j))       !surface value, L -> ZERO: define it in main cycle
          l(i,j,kb)=vk*0.01              !bottom value, L -> ZERO (almost_zero)
        end do
      end do
!
!     Estimate surface buoyancy flux, momentum flux parameter u_star and
!     inverse of Monin-Obukhov length scale LmoR.
!
      do j=1,jm
        do i=1,im
!         negative sea depth
#ifdef GEN_SIGMA
          zk=z(i,j,kb)
#else
          zk=z(kb)*dh(i,j)
#endif
!         Remove radiation that reaches bottom. It could decrease
!         surface heating stabilization impact in shallow waters
!
          rr = exp(zk/ad1(ntp))*r(ntp) +
     &         exp(zk/ad2(ntp))*(1.-r(ntp))
          hf = wtsurf(i,j)
          Bf(i,j) = grav*(alpha_sw*hf-beta_sw*wssurf(i,j))  !Bf positive for the case of convection
!          u_star = sqrt(0.5*sqrt(
!     *        (wusurf(i,j)+wusurf(i+1,j))**2+
!     *        (wvsurf(i,j)+wvsurf(i,j+1))**2))
          u_star = sqrt(sqrt(wusurf_t(i,j)**2+wvsurf_t(i,j)**2))
          LmoR(i,j)= -vk*Bf(i,j)/max(u_star**3,almost_zero)
        end do
      end do
!
!     Turbulence length scale estimation
!
      do k=1,kb-1
        do j=1,jm
          do i=1,im
!
!           Length scale in the surface ML, possibly - convective ML
!
#ifdef GEN_SIGMA
            zk = vk*(abs(z(i,j,k))+z0(i,j))
#else
            zk = vk*(abs(z(k)*dh(i,j))+z0(i,j))
#endif
            zn_MO = zk/vk*LmoR(i,j)

#ifdef GEN_SIGMA
            if(abs(z(i,j,k)) < lt(i,j))then  !v20130313, consider surface conditions only in PBL
#else
            if(abs(z(k)*dh(i,j)) < lt(i,j))then
#endif

              if(zn_MO >= 1.0)then  !stabilizing surface fluxes
                ls = alp3*zk
              elseif(zn_MO >= 0.0)then  !weekly stable to neutral stratification
                ls = zk/(1.+2.7*zn_MO)
              else                      !convective instability could develope
                ls = zk*(1.-100.*zn_MO)**0.2
              endif

            else
              ls = zk
            endif
            ls = max(ls,almost_zero)
! test case: no MO impact
!            ls = max(zk,almost_zero)
!
!           Bottom mixed layer, distance from bottom. Assume zn_MO_bottom = 0
!
#ifdef GEN_SIGMA
            lh = vk*(abs(z(i,j,kb)-z(i,j,k))+0.01)
#else
            lh = vk*(abs((z(kb)-z(k))*dh(i,j))+0.01)
#endif
!            lh = max(lh,almost_zero)
!
!           Length scale limited by the buoyancy effect
!
            if(k>1 .and. boygr(i,j,k) < 0.0) then !statically stable stratification
              N = sqrt(-boygr(i,j,k)/1.025)
              q = sqrt(abs(q2(i,j,k)))
              lb=alp2*q/N
! test case: comment out MO impact terms
#ifdef GEN_SIGMA
              if(abs(z(i,j,k)) < lt(i,j))then
#else
              if(abs(dh(i,j)*z(k)) < lt(i,j))then
#endif
!v20130313, consider surface conditions only in PBL
                if( zn_MO < 0.0 )then !count for convection impact
!                 Convective velocity scale
                  qc = (Bf(i,j)*Lt(i,j))**(1./3.)
                  lb = lb*sqrt(1.+40.*qc/(Lt(i,j)*N))
                endif
              endif

              lb=max(lb,almost_zero)
            else !statically unstable vertical stratification
              lb=1./almost_zero
            endif
!
!           Length scale controlled by the smallest length scale
!           among the three length scales: lt, lb, and (ls-lh)
!
            lr = 1./lt(i,j) + 1./lb + 1./min(ls,lh)
            l(i,j,k) = max(1./lr,almost_zero)
          end do
        end do
!        write(*,*)"MYNN a1,a2,c1,g2=",a1,a2,c1,g2
!        write(*,*)'k,z,L,ls,lh,lb,lt,zn_MO,1/l_MO=',
!     *             k,zk/vk,l(1,1,k),ls,lh,lb,lt(1,1),zn_MO,LmoR(1,1)
      end do
      end subroutine mynn_get_l
!
!******************************************************************************
! Example of shear square (shear2) estimation in calling code
!      do k=2,kb-1
!        do j=1,jmm1
!          do i=1,imm1
!            shear2=0.25*((u(i,  j,k)-u(i,  j,k-1)
!     $                   +u(i+1,j,k)-u(i+1,j,k-1))**2
!     $                  +(v(i,  j,k)-v(i,j,  k-1)
!     $                   +v(i,j+1,k)-v(i,j+1,k-1))**2)/
!#ifdef GEN_SIGMA
!     $                   dzz(i,j,k-1)**2
!#else
!     $                   (dzz(k-1)*dh(i,j))**2
!#endif
!          enddo
!        enddo
!      enddo
!******************************************************************************
!
      subroutine mynn_get_q2l2(l,q2l2,boygr,shear2,im,jm,kb)
!
!     Level 2 turbulent q square (q2l2) estimation for MYNN model
!     Values are defined on the walls of the grid boxes.
!     Boundaries (surface,bottom,lateral) must be treated in the calling code
!
      implicit none
      integer,intent(in):: im,jm,kb
      real,intent(out)::   q2l2(im,jm,kb)
      real,intent(in):: l(im,jm,kb), shear2(im,jm,kb)
      real,intent(in):: boygr(im,jm,kb)        !in POM N**2 ~= -boygr/1.025
!
!     some constants
!
!     Critical Richardson number rfc:
      real, parameter:: rfc = g1/( g1+g2 )
      real, parameter:: f1 = b1*( g1-c1 ) +3.0*a2*( 1.0-c2 )*( 1.0-c5 )
     *                       +2.0*a1*( 3.0-2.0*c2 )
      real, parameter:: f2 = b1*( g1+g2 ) -3.0*a1*( 1.0-c2 )
      real, parameter:: rf1 = b1*( g1-c1 )/f1
      real, parameter:: rf2 = b1* g1/f2
      real, parameter:: smc = a1 /a2* f1/f2
      real, parameter:: shc = 3.0*a2*( g1+g2 )
c
      real, parameter:: ri1=0.5/smc
      real, parameter:: ri2=rf1*smc
      real, parameter:: ri3=4.0*rf2*smc -2.0*ri2
      real, parameter:: ri4=ri2**2
c
      integer i,j,k
      real ri, rf, sh2, sm2
c
      do k=2,kb-1
        do j=1,jm
          do i=1,im
!
!Gradient Richardson number
!
            Ri = -boygr(i,j,k)/(1.025*max(shear2(i,j,k),almost_zero))
!
!Flux Richardson number rf, here
! rfc=0.2984, rf1=0.374, rf2=0.313
! ri1=0.806 ri2=0.232 ri3=0.313 ri4=0.0537
! D < 0, sqrt is OK
!        
            Rf = min( ri1*( Ri+ri2-sqrt(Ri**2-ri3*Ri+ri4) ), rfc )
            sh2 = shc*( rfc-Rf )/( 1.0-Rf )
            sm2 = smc*( rf1-Rf )/( rf2-Rf ) * sh2
            q2l2(i,j,k) = b1*sm2*(1.0-Rf)*
     *                    max(shear2(i,j,k),almost_zero)*l(i,j,k)**2
!            write(*,*)"Ri,shc,smc,rf,rfc,rf1,rf2,sh2,sm2,l(i,j,k),"//
!     *                "ri1,ri2,ri3,ri4="
!            write(*,*) Ri,shc,smc,rf,rfc,rf1,rf2,sh2,sm2,l(i,j,k),
!     *                 ri1,ri2,ri3,ri4
          enddo
        enddo
      enddo
      q2l2(:,:,1) = 0.
      q2l2(:,:,kb) = 0.

      end subroutine mynn_get_q2l2
!
!******************************************************************************
!
      subroutine mynn_get_ShSm(sh,sm,l,q2,q2l2,boygr,shear2,im,jm,kb)
!
!     Definition of level 2.5 Nakanishi-Nino stability functions.
!
      implicit none
      integer,intent(in):: im,jm,kb
      real,intent(out)::   sh(im,jm,kb),sm(im,jm,kb)

      real,dimension(im,jm,kb),intent(in):: l, q2, q2l2
!      real,dimension(im,jm,kb),intent(in):: shear2
!start tmp  diagnostic test
      real,dimension(im,jm,kb),intent(inout):: shear2
!end tmp  diagnostic test

      real,intent(in):: boygr(im,jm,kb)        !in POM N**2 ~= -boygr/1.025
      integer i,j,k
      real ac,gh,gm,rr,ac2,a2c2,f1,f2,f3,f4,f5,d25
!
!     vsm: limiting values for non-singular numerics
!
!      real,parameter:: gh_max = 4.55e-2
!      real,parameter:: gm_max = 1e0     !1.0e-1

!      ac = 1e0
!      ac2 = ac*ac
      a2c2 = a2*(1.-c2)

      do k=1,kb
        do j=1,jm
          do i=1,im
            if(q2l2(i,j,k) > 0.0 .and.
     *        abs(q2(i,j,k)) < q2l2(i,j,k))then
              ac = sqrt(abs(q2(i,j,k))/q2l2(i,j,k))
            else
              ac = 1.0
            endif
            rr = (l(i,j,k)**2)/max(q2(i,j,k),almost_zero)
            gh = boygr(i,j,k)/1.025*rr
            gm = shear2(i,j,k)*rr
!vsm, start: introduce limitations on gh and gm to avoid singularities
!     instead of ac that is often very small in new developing unstable layers...
!            gh = min(gh,gh_max)  !d25 could have two roots, ~0.047 and 0.554 for gm==0
!            gm = min(gm,gm_max)
!end
            ac2 = ac*ac
            f1 = 1.-3.*ac2*a2*b2*(1.-c3)*gh
            f2 = 1.-9.*ac2*a1*a2c2*gh
            f3 = f1+9.*ac2*a2*a2c2*(1.-c5)*gh
            f4 = f1-12.*ac2*a1*a2c2*gh
            f5 = 6.*ac2*a1*a1*gm
            d25 = f2*f4+f5*f3
!start tmp diagnostic tests
            if(abs(d25) <= 0.0)then
              shear2(i,j,k) = -1.+d25
            endif
!end tmp  diagnostic test
            d25 = max(d25,almost_zero)
            sm(i,j,k)=ac*a1*(f3-3.*c1*f4)/d25
            sh(i,j,k)=ac*a2*(f2+3.*c1*f5)/d25
!      if(i==485 .and. j==357 .and. k==2)then
!      write(*,*)"MYNN: l,q2,q2l2,bg,sh2,sh,sm=",
!     *l(i,j,k),q2(i,j,k),q2l2(i,j,k),boygr(i,j,k),shear2(i,j,k),
!     *sh(i,j,k),sm(i,j,k)
!      endif
!      write(*,*)"MYNN gh=",gh,gm
          enddo
        enddo
      enddo
      end subroutine mynn_get_ShSm
!
!******************************************************************************
!
      subroutine my_get_ShSm(sh,sm,l,q2,boygr,im,jm,kb)
!
!  Definition of level 2.5 Mellor-Yamada stability functions.
!
C  NOTE: Richardson # dep. dissipation correction (Mellor, 2001; Ezer, 2000) 
C  disabled here; stf=1.0 as initialized above.
C  It is unclear yet if diss. corr. is needed when surf. waves are included.
!
      implicit none
      integer,intent(in):: im,jm,kb
      real,intent(out)::   sh(im,jm,kb),sm(im,jm,kb)
      real,dimension(im,jm,kb),intent(in):: l,q2,boygr  !in POM: N**2 ~= -boygr/1.025

      integer i,j,k
      real gh,coef1,coef2,coef3,coef4,coef5
!
      real, parameter:: stf=1. !stf(i,j,k), see above for dissipation correction

      coef1=a2_my*(1.-6.*a1_my/b1_my*stf)
      coef2=3.*a2_my*(b2_my/stf+6.*a1_my)
      coef3=a1_my*(1.-3.*c1_my-6.0*a1_my/b1_my*stf)

      coef4=9.*a1_my*(2.*a1_my+a2_my)
      coef5=9.*a1_my*a2_my

      do k=1,kb
        do j=1,jm
          do i=1,im
            if(k==1 .or. k==kb)then
              gh = 0.
            else
              gh=(l(i,j,k)**2)*boygr(i,j,k)/(1.025*q2(i,j,k))
!              write(*,*)"MY gh=",gh,1./coef2,1./coef5
C
C     sm and sh limit to infinity when gh approaches 0.0288, :
C
              gh=min(gh,0.028)
            endif
!
!     vsm: for large L it is possible that sm becomes negative:
!          it was found for week stable stratification and L~100m
!
            sh(i,j,k)=coef1/(1.-coef2*gh)
            sm(i,j,k)=(coef3+sh(i,j,k)*coef4*gh)/(1.-coef5*gh)
          end do
        end do
      end do

      end subroutine my_get_ShSm

      end module

