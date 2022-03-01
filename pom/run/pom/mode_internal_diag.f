!_______________________________________________________________________
      subroutine mode_internal_diag
! calculate the internal (3-D) mode
      implicit none
      include 'pom.h'
      integer i,j,k
      real dxr,dxl,dyt,dyb

! diag
      integer ptmx,pt
!      parameter(ptmx=5)
!      real ptlon(ptmx),ptlat(ptmx)
!      data ptlon/134.8,134.9,135.1,135.3,135.5/
!      data ptlat/ 33.8, 34.0, 33.8, 33.6, 33.4/
      parameter(ptmx=1)
      real ptlon(ptmx),ptlat(ptmx)
      data ptlon/135.5/
      data ptlat/ 33.4/

      integer imdiag
!      parameter(imdiag = 60) ! 1-hour
      parameter(imdiag = 240) ! 4-hour
!      parameter(imdiag = 1440) !  1-day
      integer iglb,jglb,iloc,jloc,nproc,kloc
      real diag_sum
      include 'pom_diag.h'

      if((iint.ne.1.or.time0.ne.0.e0).and.mode.ne.2) then

! adjust u(z) and v(z) such that depth average of (u,v) = (ua,va)
        do j=1,jm
          do i=1,im
            tps(i,j)=0.e0
          end do
        end do

        do k=1,kbm1
          do j=1,jm
            do i=1,im
              tps(i,j)=tps(i,j)+u(i,j,k)*dz(k)
            end do
          end do
        end do

        do k=1,kbm1
          do j=1,jm
            do i=2,im
              u(i,j,k)=(u(i,j,k)-tps(i,j))+
     $                 (utb(i,j)+utf(i,j))/(dt(i,j)+dt(i-1,j))
            end do
          end do
        end do

        do j=1,jm
          do i=1,im
            tps(i,j)=0.e0
          end do
        end do

        do k=1,kbm1
          do j=1,jm
            do i=1,im
              tps(i,j)=tps(i,j)+v(i,j,k)*dz(k)
            end do
          end do
        end do

        do k=1,kbm1
          do j=2,jm
            do i=1,im
              v(i,j,k)=(v(i,j,k)-tps(i,j))+
     $                 (vtb(i,j)+vtf(i,j))/(dt(i,j)+dt(i,j-1))
            end do
          end do
        end do

! calculate w from u, v, dt (h+et), etf and etb
        call vertvl(a,c)

        call bcond(5)

        call exchange3d_mpi(w,im,jm,kb)

! set uf and vf to zero
        do k=1,kb
          do j=1,jm
            do i=1,im
              uf(i,j,k)=0.e0
              vf(i,j,k)=0.e0
            end do
          end do
        end do

! calculate q2f and q2lf using uf, vf, a and c as temporary variables
        call advq(q2b,q2,uf,a,c)
        call advq(q2lb,q2l,vf,a,c)
        call profq(a,c,tps,dtef)

        call bcond(6)

        call exchange3d_mpi(uf(:,:,2:kbm1),im,jm,kbm2)
        call exchange3d_mpi(vf(:,:,2:kbm1),im,jm,kbm2)

        do k=1,kb
          do j=1,jm
            do i=1,im
              q2(i,j,k)=q2(i,j,k)
     $                   +.5e0*smoth*(uf(i,j,k)+q2b(i,j,k)
     $                                -2.e0*q2(i,j,k))
              q2l(i,j,k)=q2l(i,j,k)
     $                   +.5e0*smoth*(vf(i,j,k)+q2lb(i,j,k)
     $                                -2.e0*q2l(i,j,k))
              q2b(i,j,k)=q2(i,j,k)
              q2(i,j,k)=uf(i,j,k)
              q2lb(i,j,k)=q2l(i,j,k)
              q2l(i,j,k)=vf(i,j,k)
            end do
          end do
        end do

! calculate tf and sf using uf, vf, a and c as temporary variables
        if(mode.ne.4) then
          if(nadv.eq.1) then
            call advt1(tb,t,tclim,uf,a,c)
            call advt1(sb,s,sclim,vf,a,c)
          else if(nadv.eq.2) then
! diag
!            call advt2(tb,t,tclim,uf,a,c)
            call advt2_diag(tb,t,tclim,uf,a,c)
            call advt2(sb,s,sclim,vf,a,c)
          else
            error_status=1
            write(6,'(/''Error: invalid value for nadv'')')
          end if

! diag
!          call proft(uf,wtsurf,tsurf,nbct,tps)
          call proft_diag(uf,wtsurf,tsurf,nbct,tps)
          call proft(vf,wssurf,ssurf,nbcs,tps)

          call bcond(4)

          call exchange3d_mpi(uf(:,:,1:kbm1),im,jm,kbm1)
          call exchange3d_mpi(vf(:,:,1:kbm1),im,jm,kbm1)

! diag
          if(iint.eq.1) then

          do k=1,kb
            do j=1,jm
              do i=1,im
              xadvmean(i,j,k)=0.
              yadvmean(i,j,k)=0.
              wadvmean(i,j,k)=0.
              xdifmean(i,j,k)=0.
              ydifmean(i,j,k)=0.
              vdifmean(i,j,k)=0.
              gvdmean(i,j,k)=0.
              radmean(i,j,k)=0.
              end do
            end do
          end do

          do j=1,jm
            do i=1,im
              sfcmean(i,j)=0.
            end do
          end do

          end if

          do k=1,kb
            do j=1,jm
              do i=1,im
              xadvmean(i,j,k)=xadvmean(i,j,k)
     $              +xadvterm(i,j,k)/imdiag 
              yadvmean(i,j,k)=yadvmean(i,j,k)
     $              +yadvterm(i,j,k)/imdiag 
              wadvmean(i,j,k)=wadvmean(i,j,k)
     $              +wadvterm(i,j,k)/imdiag 
              xdifmean(i,j,k)=xdifmean(i,j,k)
     $              +xdifterm(i,j,k)/imdiag 
              ydifmean(i,j,k)=ydifmean(i,j,k)
     $              +ydifterm(i,j,k)/imdiag 
              vdifmean(i,j,k)=vdifmean(i,j,k)
     $              +vdifterm(i,j,k)/imdiag 
              gvdmean(i,j,k)=gvdmean(i,j,k)
     $              +gvdterm(i,j,k)/imdiag 
              radmean(i,j,k)=radmean(i,j,k)
     $              +radterm(i,j,k)/imdiag 
              end do
            end do
          end do

          do j=1,jm
            do i=1,im
              sfcmean(i,j)=sfcmean(i,j)+sfcterm(i,j)/imdiag
            end do
          end do

          if(mod(iint,imdiag).eq.0) then 

          call exchange2d_mpi(sfcmean,im,jm)       
          call exchange2d_mpi(windu,im,jm)       
          call exchange2d_mpi(windv,im,jm)       
          call exchange3d_mpi(xadvmean(:,:,1:kbm1),im,jm,kbm1)
          call exchange3d_mpi(yadvmean(:,:,1:kbm1),im,jm,kbm1)
          call exchange3d_mpi(wadvmean(:,:,1:kbm1),im,jm,kbm1)
          call exchange3d_mpi(xdifmean(:,:,1:kbm1),im,jm,kbm1)
          call exchange3d_mpi(ydifmean(:,:,1:kbm1),im,jm,kbm1)
          call exchange3d_mpi(radmean(:,:,1:kbm1),im,jm,kbm1)  
          call exchange3d_mpi(vdifmean(:,:,1:kbm1),im,jm,kbm1) 
          call exchange3d_mpi(gvdmean(:,:,1:kbm1),im,jm,kbm1) 
          call exchange3d_mpi(uf(:,:,1:kbm1),im,jm,kbm1) 

          call write_diag_netcdf

          do pt=1,ptmx

            iglb=nint((ptlon(pt)-133)*36+2) 
            jglb=nint((ptlat(pt)- 30)*36+2)
            kloc=1
            nproc=int((iglb-2)/(im-2))+int((jglb-2)/(jm-2))*nproc_x

            if(my_task.eq.nproc) then

            iloc=iglb-mod(nproc,nproc_x)*(im-2)
            jloc=jglb-(nproc/nproc_x)*(jm-2)
            write(99,*) 'diag_temp iint ',iint
            write(99,*) 'lon lat ',
     $         east_e(iloc,jloc),north_e(iloc,jloc)
            write(99,*) 'iloc jloc nproc ',iloc,jloc,nproc
            write(99,*) 'iglb jglb kloc ',iglb,jglb,kloc
            write(99,*) 'iint pt xadv ',
     $        iint, pt, xadvmean(iloc,jloc,kloc)
            write(99,*) 'iint pt yadv ', 
     $        iint, pt, yadvmean(iloc,jloc,kloc)
            write(99,*) 'iint pt wadv ',
     $        iint, pt, wadvmean(iloc,jloc,kloc)
            write(99,*) 'iint pt xdif ',
     $        iint, pt, xdifmean(iloc,jloc,kloc)
            write(99,*) 'iint pt ydif ',
     $        iint, pt, ydifmean(iloc,jloc,kloc)
            write(99,*) 'iint pt vdif ',
     $        iint, pt, vdifmean(iloc,jloc,kloc)
            write(99,*) 'iint pt rad  ',
     $        iint, pt, radmean(iloc,jloc,kloc)
            write(99,*) 'iint pt gvd  ',
     $        iint, pt, gvdmean(iloc,jloc,kloc)
            diag_sum = ( xadvmean(iloc,jloc,kloc)
     $                  +yadvmean(iloc,jloc,kloc)
     $                  +wadvmean(iloc,jloc,kloc)
     $                  +xdifmean(iloc,jloc,kloc)
     $                  +ydifmean(iloc,jloc,kloc)
     $                  +gvdmean(iloc,jloc,kloc)
     $                 )*dti2*imdiag
            write(99,*) 'iint pt tend*dti2  ',
     $        iint, pt, diag_sum, diag_sum/dti2/imdiag
            write(99,*) 'iint pt sfcf ',
     $        iint, pt, sfcmean(iloc,jloc)
            if(kloc.eq.1)
     $        write(99,*) 'iint pt vdif+rad+sfcf ',
     $        iint, pt, vdifterm(iloc,jloc,kloc)
     $                 +radterm(iloc,jloc,kloc)
     $                 +sfcterm(iloc,jloc)
            diag_sum = ( xadvterm(iloc,jloc,kloc)
     $                  +yadvterm(iloc,jloc,kloc)
     $                  +wadvterm(iloc,jloc,kloc)
     $                  +xdifterm(iloc,jloc,kloc)
     $                  +ydifterm(iloc,jloc,kloc)
     $                  +gvdterm(iloc,jloc,kloc)
     $                 )*dti2
            write(99,*) 'iint pt sum*dti2 ', iint, pt, diag_sum 
            write(99,*) 'i p t t+ tf ', 
     $                  iint,pt,
     $                  tb(iloc,jloc,kloc),
     $                  tb(iloc,jloc,kloc)+diag_sum, 
     $                  uf(iloc,jloc,kloc)
            endif

          end do

          do k=1,kb
            do j=1,jm
              do i=1,im
              xadvmean(i,j,k)=0.
              yadvmean(i,j,k)=0.
              wadvmean(i,j,k)=0.
              xdifmean(i,j,k)=0.
              ydifmean(i,j,k)=0.
              vdifmean(i,j,k)=0.
              gvdmean(i,j,k)=0.
              radmean(i,j,k)=0.
              end do
            end do
          end do

          do j=1,jm
            do i=1,im
              sfcmean(i,j)=0.
            end do
          end do

          end if

          do k=1,kb
            do j=1,jm
              do i=1,im
                t(i,j,k)=t(i,j,k)
     $                    +.5e0*smoth*(uf(i,j,k)+tb(i,j,k)
     $                                 -2.e0*t(i,j,k))
                s(i,j,k)=s(i,j,k)
     $                    +.5e0*smoth*(vf(i,j,k)+sb(i,j,k)
     $                                 -2.e0*s(i,j,k))
                tb(i,j,k)=t(i,j,k)
                t(i,j,k)=uf(i,j,k)
                sb(i,j,k)=s(i,j,k)
                s(i,j,k)=vf(i,j,k)
              end do
            end do
          end do

          call dens(s,t,rho)

        end if

! calculate uf and vf
        call advu
        call advv
        call profu
        call profv

        call bcond(3)

        call exchange3d_mpi(uf(:,:,1:kbm1),im,jm,kbm1)
        call exchange3d_mpi(vf(:,:,1:kbm1),im,jm,kbm1)

        do j=1,jm
          do i=1,im
            tps(i,j)=0.e0
          end do
        end do

        do k=1,kbm1
          do j=1,jm
            do i=1,im
              tps(i,j)=tps(i,j)
     $                  +(uf(i,j,k)+ub(i,j,k)-2.e0*u(i,j,k))*dz(k)
            end do
          end do
        end do

        do k=1,kbm1
          do j=1,jm
            do i=1,im
              u(i,j,k)=u(i,j,k)
     $                  +.5e0*smoth*(uf(i,j,k)+ub(i,j,k)
     $                               -2.e0*u(i,j,k)-tps(i,j))
            end do
          end do
        end do

        do j=1,jm
          do i=1,im
            tps(i,j)=0.e0
          end do
        end do

        do k=1,kbm1
          do j=1,jm
            do i=1,im
              tps(i,j)=tps(i,j)
     $                  +(vf(i,j,k)+vb(i,j,k)-2.e0*v(i,j,k))*dz(k)
            end do
          end do
        end do

        do k=1,kbm1
          do j=1,jm
            do i=1,im
              v(i,j,k)=v(i,j,k)
     $                  +.5e0*smoth*(vf(i,j,k)+vb(i,j,k)
     $                               -2.e0*v(i,j,k)-tps(i,j))
            end do
          end do
        end do

        do k=1,kb
          do j=1,jm
            do i=1,im
              ub(i,j,k)=u(i,j,k)
              u(i,j,k)=uf(i,j,k)
              vb(i,j,k)=v(i,j,k)
              v(i,j,k)=vf(i,j,k)
            end do
          end do
        end do

      end if

      do j=1,jm
        do i=1,im
          egb(i,j)=egf(i,j)
          etb(i,j)=et(i,j)
          et(i,j)=etf(i,j)
          dt(i,j)=h(i,j)+et(i,j)
          utb(i,j)=utf(i,j)
          vtb(i,j)=vtf(i,j)
          vfluxb(i,j)=vfluxf(i,j)
        end do
      end do

! calculate real w as wr
      do k=1,kb
        do j=1,jm
          do i=1,im
            wr(i,j,k)=0.
          end do
        end do
      end do

      do k=1,kbm1
        do j=1,jm
          do i=1,im
            tps(i,j)=zz(k)*dt(i,j) + et(i,j)
          end do
        end do
        do j=2,jmm1
          do i=2,imm1
            dxr=2.0/(dx(i+1,j)+dx(i,j))
            dxl=2.0/(dx(i,j)+dx(i-1,j))
            dyt=2.0/(dy(i,j+1)+dy(i,j))
            dyb=2.0/(dy(i,j)+dy(i,j-1))
            wr(i,j,k)=0.5*(w(i,j,k)+w(i,j,k+1))+0.5*
     $                (u(i+1,j,k)*(tps(i+1,j)-tps(i,j))*dxr+
     $                 u(i,j,k)*(tps(i,j)-tps(i-1,j))*dxl+
     $                 v(i,j+1,k)*(tps(i,j+1)-tps(i,j))*dyt+
     $                 v(i,j,k)*(tps(i,j)-tps(i,j-1))*dyb)
     $                +(1.0+zz(k))*(etf(i,j)-etb(i,j))/dti2
          end do
        end do
      end do

      call exchange3d_mpi(wr(:,:,1:kbm1),im,jm,kbm1)

      do k=1,kb
        do i=1,im
          if(n_south.eq.-1) wr(i,1,k)=wr(i,2,k)
          if(n_north.eq.-1) wr(i,jm,k)=wr(i,jmm1,k)
        end do
      end do
      do k=1,kb
        do j=1,jm
          if(n_west.eq.-1) wr(1,j,k)=wr(2,j,k)
          if(n_east.eq.-1) wr(im,j,k)=wr(imm1,j,k)
        end do
      end do

      do k=1,kbm1
        do j=1,jm
          do i=1,im
            wr(i,j,k)=fsm(i,j)*wr(i,j,k)
          end do
        end do
      end do

      return
      end

