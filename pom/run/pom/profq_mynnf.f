!_______________________________________________________________________
      subroutine profq_mynnf
!     solve for q2 (twice the turbulent kinetic energy), q2l (q2 x turbulent
!     length scale), km (vertical kinematic viscosity) and kh (vertical
!     kinematic diffusivity), using a simplified version of the level 2 1/2
!     model of Mellor and Yamada (1982)
!     in this version, the Craig-Banner sub-model whereby breaking wave tke
!     is injected into the surface is included. However, we use an
!     analytical solution to the near surface tke equation to solve for q2
!     at the surface giving the same result as C-B diffusion. The new scheme
!     is simpler and more robust than the latter scheme
!     
!     including Mellor-Yamada-Nakanishi-Niino-Furuichi version
!     see lmynnf
!     
      use MYNNF_lev25_2012

      implicit none
      include 'pom.h'
      real a(im,jm,kb),c(im,jm,kb)
      real ee(im,jm,kb),gg(im,jm,kb)
      real sm(im,jm,kb),sh(im,jm,kb)
      real cc(im,jm,kb)
      real gh(im,jm,kb),boygr(im,jm,kb),dh(im,jm),stf(im,jm,kb)
      real prod(im,jm,kb)
      real a1,a2,b1,b2,c1
      real coef1,coef2,coef3,coef4,coef5
      real const1,e1,e2
      real p,sef,sp,tp
      real l0(im,jm)
      real cbcnst,surfl,shiw
      real utau2(im,jm)
      real df0,df1,df2
      integer i,j,k,ki

      data a1,b1,a2,b2,c1/0.92e0,16.6e0,0.74e0,10.1e0,0.08e0/
      data e1/1.8e0/,e2/1.33e0/
      data sef/1.e0/

!     surface wave breaking
!     data cbcnst/100./surfl/2.e5/
      data cbcnst/50./surfl/1.e5/

      if(liwbrk) then
         shiw=0.7
      else
         shiw=0.0
      end if

      do j=1,jm
         do i=1,im
            dh(i,j)=h(i,j)+etf(i,j)
         end do
      end do

      do k=1,kb
         do j=1,jm
            do i=1,im
               a(i,j,k)=0.e0
               c(i,j,k)=0.e0
               ee(i,j,k)=0.e0
               gg(i,j,k)=0.e0
            end do
         end do
      end do

      do k=2,kbm1
         do j=2,jmm1
            do i=2,imm1
               a(i,j,k)=-dti2*(kq(i,j,k+1)+kq(i,j,k)+2.e0*umol)*.5e0
     $              /(dzz(k-1)*dz(k)*dh(i,j)*dh(i,j))
               c(i,j,k)=-dti2*(kq(i,j,k-1)+kq(i,j,k)+2.e0*umol)*.5e0
     $              /(dzz(k-1)*dz(k-1)*dh(i,j)*dh(i,j))
            end do
         end do
      end do

!     write(6,*) 'in profq 0'
!     call flush(6)

!     the following section solves the equation:
!     dti2*(kq*q2')' - q2*(2.*dti2*dtef+1.) = -q2b

!     surface and bottom boundary conditions
      const1=(16.6e0**(2.e0/3.e0))*sef
      if(lmynnf) const1=(24.0e0**(2.e0/3.e0))*sef

!     initialize fields that are not calculated on all boundaries
!     but are later used there
      do i=1,im
         do j=1,jm
            l0(i,j)=0.
         end do
      end do
      do i=1,im
         do j=1,jm
            do k=1,kb
               boygr(i,j,k)=0.
               prod(i,j,k)=0.
            end do
         end do
      end do

      do j=2,jmm1
         do i=2,imm1
            utau2(i,j)=sqrt((.5e0*(wusurf(i,j)+wusurf(i+1,j)))**2
     $           +(.5e0*(wvsurf(i,j)+wvsurf(i,j+1)))**2)
            uf(i,j,kb)=sqrt((.5e0*(wubot(i,j)+wubot(i+1,j)))**2
     $           +(.5e0*(wvbot(i,j)+wvbot(i,j+1)))**2)*const1
         end do
      end do

!     call exchange2d_mpi(utau2,im,jm)
!     call exchange2d_mpi(uf(:,:,kb),im,jm)

      if(lwbrk) then

         do j=2,jmm1
            do i=2,imm1
!     wave breaking energy- a variant of Craig & Banner (1994)
!     see Mellor and Blumberg, 2004.
               ee(i,j,1)=0.e0
               gg(i,j,1)=(15.8*cbcnst)**(2./3.)*utau2(i,j)
!     surface length scale following Stacey (1999).
               l0(i,j)=surfl*utau2(i,j)/grav
            end do
         end do

      end if

!     calculate speed of sound squared
      do k=1,kb
         do j=1,jm
            do i=1,im
               cc(i,j,k)=0.
            end do
         end do
      end do
      do k=1,kbm1
         do j=2,jmm1
            do i=2,imm1
               tp=t(i,j,k)+tbias
               sp=s(i,j,k)+sbias
!     calculate pressure in units of decibars
               p=grav*rhoref*(-zz(k)*h(i,j))*1.e-4
               cc(i,j,k)=1449.1e0+.00821e0*p+4.55e0*tp-.045e0*tp**2
     $              +1.34e0*(sp-35.0e0)
               cc(i,j,k)=cc(i,j,k)
     $              /sqrt((1.e0-.01642e0*p/cc(i,j,k))
     $              *(1.e0-0.40e0*p/cc(i,j,k)**2))
            end do
         end do
      end do

!     calculate buoyancy gradient
      do k=2,kbm1
         do j=2,jmm1
            do i=2,imm1
               q2b(i,j,k)=abs(q2b(i,j,k))
               q2lb(i,j,k)=abs(q2lb(i,j,k))
               boygr(i,j,k)=grav*(rho(i,j,k-1)-rho(i,j,k))
     $              /(dzz(k-1)*h(i,j))
!     *** note: comment out next line if dens does not include pressure
     $              +(grav**2)*2.e0/(cc(i,j,k-1)**2+cc(i,j,k)**2)
            end do
         end do
      end do

      if(.not.lmynnf) then

         do k=2,kbm1
            do j=2,jmm1
               do i=2,imm1
                  l(i,j,k)=abs(q2lb(i,j,k)/q2b(i,j,k))
                  if(z(k).gt.-0.5) l(i,j,k)=max(l(i,j,k),kappa*l0(i,j))
                  gh(i,j,k)=(l(i,j,k)**2)*boygr(i,j,k)/q2b(i,j,k)
                  gh(i,j,k)=min(gh(i,j,k),.028e0)
               end do
            end do
         end do

         do j=2,jmm1
            do i=2,imm1
               l(i,j,1)=kappa*l0(i,j)
               l(i,j,kb)=0.e0
               gh(i,j,1)=0.e0
               gh(i,j,kb)=0.e0
            end do
         end do

      end if

!     calculate production of turbulent kinetic energy:
      do k=2,kbm1
         do j=2,jmm1
            do i=2,imm1
               prod(i,j,k)=km(i,j,k)*.25e0*sef
     $              *((u(i,j,k)-u(i,j,k-1)
     $              +u(i+1,j,k)-u(i+1,j,k-1))**2
     $              +(v(i,j,k)-v(i,j,k-1)
     $              +v(i,j+1,k)-v(i,j+1,k-1))**2)
     $              /(dzz(k-1)*dh(i,j))**2
!     add shear due to internal wave field
     $              -shiw*km(i,j,k)*boygr(i,j,k)
               prod(i,j,k)=prod(i,j,k)+kh(i,j,k)*boygr(i,j,k)
            end do
         end do
      end do
!     call exchange3d_mpi(prod(:,:,2:kbm1),im,jm,kbm2)

!     write(6,*) 'in profq b'
!     call flush(6)

      if(lmynnf) then

         do k=1,kb
            do j=2,jmm1
               do i=2,imm1
                  stf(i,j,k)=1.e0
                  dtef(i,j,k)=sqrt(abs(q2b(i,j,k)))*stf(i,j,k)
     $                 /(24.0e0*l(i,j,k)+small)
               end do
            end do
         end do

      else

         do k=1,kb
            do j=2,jmm1
               do i=2,imm1
                  stf(i,j,k)=1.e0
                  dtef(i,j,k)=sqrt(abs(q2b(i,j,k)))*stf(i,j,k)
     $                 /(b1*l(i,j,k)+small)
               end do
            end do
         end do

      end if

      do k=2,kbm1
         do j=2,jmm1
            do i=2,imm1
               gg(i,j,k)=1.e0/(a(i,j,k)+c(i,j,k)*(1.e0-ee(i,j,k-1))
     $              -(2.e0*dti2*dtef(i,j,k)+1.e0))
               ee(i,j,k)=a(i,j,k)*gg(i,j,k)
               gg(i,j,k)=(-2.e0*dti2*prod(i,j,k)+c(i,j,k)*gg(i,j,k-1)
     $              -uf(i,j,k))*gg(i,j,k)
            end do
         end do
      end do

      do k=1,kbm1
         ki=kb-k
         do j=2,jmm1
            do i=2,imm1
               uf(i,j,ki)=ee(i,j,ki)*uf(i,j,ki+1)+gg(i,j,ki)
            end do
         end do
      end do

      if(.not.lmynnf) then

!     the following section solves the equation:
!     dti2(kq*q2l')' - q2l*(dti2*dtef+1.) = -q2lb

         do j=2,jmm1
            do i=2,imm1
               vf(i,j,1)=0.
               vf(i,j,kb)=0.
               ee(i,j,2)=0.e0
               gg(i,j,2)=-kappa*z(2)*dh(i,j)*q2(i,j,2)
               vf(i,j,kb-1)=kappa*(1+z(kbm1))*dh(i,j)*q2(i,j,kbm1)
            end do
         end do
         do k=2,kbm1
            do j=2,jmm1
               do i=2,imm1
                  dtef(i,j,k)=dtef(i,j,k)
     $                 *(1.e0+e2*((1.e0/abs(z(k)-z(1))
     $                 +1.e0/abs(z(k)-z(kb)))
     $                 *l(i,j,k)/(dh(i,j)*kappa))**2)
               end do
            end do
         end do
         do k=3,kbm1
            do j=2,jmm1
               do i=2,imm1
                  gg(i,j,k)=1.e0/(a(i,j,k)+c(i,j,k)*(1.e0-ee(i,j,k-1))
     $                 -(dti2*dtef(i,j,k)+1.e0))
                  ee(i,j,k)=a(i,j,k)*gg(i,j,k)
                  gg(i,j,k)=(dti2*(-prod(i,j,k)*l(i,j,k)*e1)
     $                 +c(i,j,k)*gg(i,j,k-1)-vf(i,j,k))*gg(i,j,k)
               end do
            end do
         end do

         do k=1,kb-2
            ki=kb-k
            do j=2,jmm1
               do i=2,imm1
                  vf(i,j,ki)=ee(i,j,ki)*vf(i,j,ki+1)+gg(i,j,ki)
               end do
            end do
         end do

      end if

!     the following is to counter the problem of the ratio of two small
!     numbers (l = q2l/q2) or one number becoming negative. Two options are
!     included below. In this application, the second option, l was less
!     noisy when uf or vf is small
      do k=2,kbm1
         do j=2,jmm1
            do i=2,imm1
               uf(i,j,k)=abs(uf(i,j,k))
               vf(i,j,k)=abs(vf(i,j,k))
            end do
         end do
      end do

      if(lmynnf) then

!     gh = shear2
         gh(:,:,:) = 0.
         do k=2,kb
            do j=2,jmm1
               do i=2,imm1
                  if(dzz(k-1)*dh(i,j).gt.0.) then
                     gh(i,j,k)=0.25*                   
     &                    ((u(i,  j,k)-u(i,  j,k-1)           
     &                    +u(i+1,j,k)-u(i+1,j,k-1))**2       
     &                    +(v(i,  j,k)-v(i,j,  k-1)           
     &                    +v(i,j+1,k)-v(i,j+1,k-1))**2)/     
     &                    (dzz(k-1)*dh(i,j))**2
                  endif
               enddo
            enddo
         enddo      

         call mynn_get_l(l,UF,boygr,wtsurf,wssurf,     
     $        wusurf,wvsurf,l0,z,dzz,dh,im,jm,kb,ntp)
!     use vf as storage for q2l2 - level 2 diagnostic of q2
         call mynn_get_q2l2(l,vf,boygr,gh,im,jm,kb)
         call mynn_get_ShSm(sh,sm,l,uf,vf,boygr,gh,im,jm,kb)

      else

!     the following section solves for km and kh
         coef4=18.e0*a1*a1+9.e0*a1*a2
         coef5=9.e0*a1*a2

!     note that sm and sh limit to infinity when gh approaches 0.0288
         do k=1,kb
            do j=2,jmm1
               do i=2,imm1
                  coef1=a2*(1.e0-6.e0*a1/b1*stf(i,j,k))
                  coef2=3.e0*a2*b2/stf(i,j,k)+18.e0*a1*a2
                  coef3=a1*(1.e0-3.e0*c1-6.e0*a1/b1*stf(i,j,k))
                  sh(i,j,k)=coef1/(1.e0-coef2*gh(i,j,k))
                  sm(i,j,k)=coef3+sh(i,j,k)*coef4*gh(i,j,k)
                  sm(i,j,k)=sm(i,j,k)/(1.e0-coef5*gh(i,j,k))
               end do
            end do
         end do

      end if

      if(lmynnf) then

         do k=1,kb
            do j=2,jmm1
               do i=2,imm1
                  prod(i,j,k)=l(i,j,k)*sqrt(abs(uf(i,j,k)))
                  kq(i,j,k)=(prod(i,j,k)*3.0e0*sm(i,j,k)+kq(i,j,k))*.5e0
                  km(i,j,k)=(prod(i,j,k)*sm(i,j,k)+km(i,j,k))*.5e0
                  kh(i,j,k)=(prod(i,j,k)*sh(i,j,k)+kh(i,j,k))*.5e0
               end do
            end do
         end do

      else

         do k=1,kb
            do j=2,jmm1
               do i=2,imm1
                  prod(i,j,k)=l(i,j,k)*sqrt(abs(q2(i,j,k)))
                  kq(i,j,k)=(prod(i,j,k)*.41e0*sh(i,j,k)+kq(i,j,k))*.5e0
                  km(i,j,k)=(prod(i,j,k)*sm(i,j,k)+km(i,j,k))*.5e0
                  kh(i,j,k)=(prod(i,j,k)*sh(i,j,k)+kh(i,j,k))*.5e0
               end do
            end do
         end do

      end if

!     cosmetics: make boundr. values as interior (even if not used, printout
!     may show strange values)
      do k=1,kb
         do i=1,im
            if(n_north.eq.-1) then
               km(i,jm,k)=km(i,jmm1,k)
               kh(i,jm,k)=kh(i,jmm1,k)
               kq(i,jm,k)=kq(i,jmm1,k)
            end if
            if(n_south.eq.-1) then
               km(i,1,k)=km(i,2,k)
               kh(i,1,k)=kh(i,2,k)
               kq(i,1,k)=kq(i,2,k)
            end if
         end do
         do j=1,jm
            if(n_east.eq.-1) then
               km(im,j,k)=km(imm1,j,k)
               kh(im,j,k)=kh(imm1,j,k)
               kq(im,j,k)=kq(imm1,j,k)
            end if
            if(n_west.eq.-1) then
               km(1,j,k)=km(2,j,k)
               kh(1,j,k)=kh(2,j,k)
               kq(1,j,k)=kq(2,j,k)
            end if
         end do
      end do

      do k=1,kb
         do i=1,im
            do j=1,jm
               km(i,j,k)=km(i,j,k)*fsm(i,j)
               kh(i,j,k)=kh(i,j,k)*fsm(i,j)
               kq(i,j,k)=kq(i,j,k)*fsm(i,j)
            end do
         end do
      end do

      call exchange3d_mpi(km,im,jm,kb) ! 2018.07
      call exchange3d_mpi(kh,im,jm,kb) ! 2018.07
      call exchange3d_mpi(kq,im,jm,kb) ! 2018.08
      call exchange3d_mpi(l,im,jm,kb) !2018.08

      return
      end
