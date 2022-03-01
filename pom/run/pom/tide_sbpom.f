! tide_sbpom.f

!_______________________________________________________________________
      subroutine tide_initialize
! initialize tide module

      use jpom_tide

      implicit none
      include 'pom.h'
      integer i,j,k
      character*3 cha_pe

      double precision:: xlon(im), ylat(jm)
      double precision:: xlonu(im), ylatu(jm)
      double precision:: xlonv(im), ylatv(jm)
      
      INTEGER IDEV_TIDE_DIAG,IDEV_TIDE_INPT,IRET

      logical:: L_TIDE_BD_INI, L_TIDE_FC_INI, L_TIDE_UVH_INI
     *         ,L_TIDE_INITIALIZED = .false.

      CHARACTER(80):: TIDE_DB_PATH=""
      INTEGER::       TIDE_MODEL=-1

      call init_iomngr

      IDEV_TIDE_DIAG=6
      IDEV_TIDE_INPT=6

      L_TIDE_BD=.true.
      L_TIDE_FC=.false.
      L_TIDE_UVH=.true.
      SEA_TIDE%PATH = "./TIDE"
      SEA_TIDE%TIDE_MODEL = TIDE_NAO99_JAPAN_SSH_UVH

      call null_init_tide(dti)

      L_TIDE_BD_INI = L_TIDE_BD
      L_TIDE_FC_INI = L_TIDE_FC
      L_TIDE_UVH_INI = L_TIDE_UVH

!      if(nread_rst > 0)then
!         guaranty input of tide mean restart data and
!         no tide amplitude dumping from the start of simulation
!          UTC_TM_JDS_START_TIDE = 86400D0*
!     &      (julday_start-max(R_DAYS_DAMP,abs(START_TIDE_DELAY_DAYS)))
!      else
          UTC_TM_JDS_START_TIDE =
     &      (julday_start+START_TIDE_DELAY_DAYS)*86400D0
!      endif

      xlon=east_c(:,1)
      ylat=north_c(1,:)
      xlonu=east_u(:,1)
      ylatu=north_u(1,:)
      xlonv=east_v(:,1)
      ylatv=north_v(1,:)
      write(cha_pe,'(i3.3)') my_task

      CALL JPOM_INIT_TIDE(UTC_TM_JDS_START_TIDE,IM,JM,KB,
     *      XLON,YLAT,XLONU,YLATU,XLONV,YLATV,
     *      idev_tide_diag,idev_tide_inpt,
     *      IRET,'tideinit.bin.'//cha_pe)
      if(iret /= 0)then
          write(6,*)
     *    '$$$ Error initializing tides: ',iret,my_task
          stop
      endif

      L_TIDE_INITIALIZED = .true.

      return
      end

!_______________________________________________________________________
      subroutine tide_advance
! estimate tide parameters at the start and end time of internal step;
! then linearly interpolated each external time step.

      use jpom_tide

      implicit none
      include 'pom.h'
      integer i,j,k

      d_tmp1 = DBLE(julday_start+time0)*86400D0 + DBLE(DTI*(IINT-1.))
      d_tmp2 = d_tmp1 + DBLE(DTI)
      call GET_TIDE2(d_tmp1,d_tmp2)
      tide_io_time = d_tmp1

!      do j = 464,462,-1
!      write(6,*) 
!     $  (fsm(i,j),i=56,58)
!      end do
!      stop
!      i = 2
!      j = jm/2  
!      i = 57
!      j = 463
!      write(6,*) 'time h u v ', 
!     $  time, zn_tide(i,j), 
!     $  uvhc_tide(i,j,1)/h(i,j), 
!     $  uvhc_tide(i,j,2)/h(i,j)
!      write(6,*) 'time el u v ', 
!     $  time, el(i,j), 
!     $  u(i+1,j-1,1),
!     $  v(i,j,1)

      return
      end

!_______________________________________________________________________
      subroutine tide_advance_external
! for larger time steps avoid steepnes of the sea level changes,
! do tide interpolation into the fast (external) mode cycle

      use jpom_tide

      implicit none
      include 'pom.h'
      integer i,j,k

      tide_io_time = tide_io_time + DBLE(DTE)
      call GET_TIDE2(tide_io_time)

      return
      end

!_______________________________________________________________________
      subroutine tide_bcond(idx)
! apply open boundary conditions by nesting with tide forcing

! tide
      use jpom_tide

      implicit none
      include 'pom.h'
      integer idx
      integer i,j,k
      real ga,u1,wm
      real hmax

! tide
!     L_TIDE_BD=.true.
!     L_TIDE_UVH=.true.
      real rf

      if(idx.eq.1) then

! External (2-D) elevation boundary conditions
        do j=1,jm
          if(n_west.eq.-1) elf(1,j)=elf(2,j)
          if(n_east.eq.-1) elf(im,j)=elf(im-1,j)
        end do

        do i=1,im
          if(n_south.eq.-1) elf(i,1)=elf(i,2)
          if(n_north.eq.-1) elf(i,jm)=elf(i,jm-1)
        end do

        do j=1,jm
          do i=1,im
            elf(i,j)=elf(i,j)*fsm(i,j)
          end do
        end do

        return

      else if(idx.eq.2) then

! external (2-D) velocity boundary conditions
        do j=2,jmm1
          ! west
          if(n_west.eq.-1) then
! tide
            uaf(2,j)=ramp*uabw(j)
     $               +2.e0*uvhc_tide(2,j,1)
     $               /(h(1,j)+h(2,j)+el(1,j)+el(2,j))
            rf=(ramp*elw(j)+zn_tide(2,j)-el(2,j))*sqrt(grav/h(2,j))
            if((j.eq.2 .and. dvm(2,2).eq.1.) .or.
     $         (j.eq.jmm1 .and. dvm(2,jm).eq.1.)) then
              rf=.5e0*rf
            end if
            uaf(1,j)=uaf(2,j)+rf
            vaf(1,j)=ramp*vabw(j)
          end if

          ! east
          if(n_east.eq.-1) then
! tide
            uaf(im,j)=ramp*uabe(j)
     $                +2.e0*uvhc_tide(im,j,1)
     $                /(h(imm1,j)+h(im,j)+el(imm1,j)+el(im,j))
            rf=-(ramp*ele(j)+zn_tide(imm1,j)-el(imm1,j))
     $          *sqrt(grav/h(imm1,j))
            if((j.eq.2 .and. dvm(imm1,2).eq.1.) .or.
     $         (j.eq.jmm1 .and. dvm(imm1,jm).eq.1.)) then
              rf=.5e0*rf
            end if
            uaf(im,j)=uaf(im,j)+rf
            vaf(im,j)=ramp*vabe(j)
          end if
        end do

        do i=2,imm1
          ! south
          if(n_south.eq.-1) then
! tide
            vaf(i,2)=ramp*vabs(i)
     $               +2.e0*uvhc_tide(i,2,2)
     $               /(h(i,1)+h(i,2)+el(i,1)+el(i,2))
            rf=(ramp*els(i)+zn_tide(i,2)-el(i,2))
     $         *sqrt(grav/h(i,2))
            if((i.eq.2 .and. dum(2,2).eq.1.) .or.
     $         (i.eq.imm1 .and. dum(im,2).eq.1.)) then
              rf=.5e0*rf
            end if
            vaf(i,2)=vaf(i,2)+rf
            vaf(i,1)=vaf(i,2)
            uaf(i,1)=ramp*uabs(i)
          end if

          ! north
          if(n_north.eq.-1) then
! tide
            vaf(i,jm)=ramp*vabn(i)
     $               +2.e0*uvhc_tide(i,jm,2)
     $               /(h(i,jmm1)+h(i,jm)+el(i,jmm1)+el(i,jm))
            rf=-(ramp*eln(i)+zn_tide(i,jmm1)-el(i,jmm1))
     $          *sqrt(grav/h(i,jmm1))
            if((i.eq.2 .and. dum(2,jmm1).eq.1.) .or.
     $         (i.eq.imm1 .and. dum(im,jmm1).eq.1.)) then
              rf=.5e0*rf
            end if
            vaf(i,jm)=vaf(i,jm)+rf
            uaf(i,jm)=ramp*uabn(i)
          end if
        end do

        do j=1,jm
          do i=1,im
            uaf(i,j)=uaf(i,j)*dum(i,j)
            vaf(i,j)=vaf(i,j)*dvm(i,j)
          end do
        end do

        return

      else if(idx.eq.3) then

! internal (3-D) velocity boundary conditions
! radiation conditions
! smoothing is used in the direction tangential to the boundaries
        hmax=4500.e0

        do k=1,kbm1
          do j=2,jmm1
            ! east
            if(n_east.eq.-1) then
              uf(im,j,k)=(utb(im,j)+utf(im,j))
     $                   /(dt(im,j)+dt(imm1,j))
     $                   +ramp*(ube(j,k)-uabe(j))
              vf(im,j,k)=ramp*vbe(j,k)
            end if
            ! west
            if(n_west.eq.-1) then
              uf(2,j,k)=(utb(2,j)+utf(2,j))
     $                  /(dt(2,j)+dt(1,j))
     $                  +ramp*(ubw(j,k)-uabw(j))
              uf(1,j,k)=uf(2,j,k)
              vf(1,j,k)=ramp*vbw(j,k)
            end if
          end do
        end do

        do k=1,kbm1
          do i=2,imm1
            ! south
            if(n_south.eq.-1) then
              vf(i,2,k)=(vtb(i,2)+vtf(i,2))
     $                  /(dt(i,2)+dt(i,1))
     $                  +ramp*(vbs(i,k)-vabs(i))
              vf(i,1,k)=vf(i,2,k)
              uf(i,1,k)=ramp*ubs(i,k)
            end if
            ! north
            if(n_north.eq.-1) then
              vf(i,jm,k)=(vtb(i,jm)+vtf(i,jm))
     $                  /(dt(i,jm)+dt(i,jmm1))
     $                   +ramp*(vbn(i,k)-vabn(i))
              uf(i,jm,k)=ramp*ubn(i,k)
            end if
          end do
        end do

        do k=1,kbm1
          do j=1,jm
            do i=1,im
              uf(i,j,k)=uf(i,j,k)*dum(i,j)
              vf(i,j,k)=vf(i,j,k)*dvm(i,j)
            end do
          end do
        end do

        return

      else if(idx.eq.4) then

! temperature and salinity boundary conditions (using uf and vf,
! respectively)
        do k=1,kbm1
          do j=1,jm
            ! east
            if(n_east.eq.-1) then
              uf(im,j,k)=tbe(j,k) 
              vf(im,j,k)=sbe(j,k) 
            end if

            ! west
            if(n_west.eq.-1) then
              uf(1,j,k)=tbw(j,k)
              vf(1,j,k)=sbw(j,k)
            end if
          end do

          do i=1,im
            ! south
            if(n_south.eq.-1) then
              uf(i,1,k)=tbs(i,k)
              vf(i,1,k)=sbs(i,k)
            end if

            ! north
            if(n_north.eq.-1) then
              uf(i,jm,k)=tbn(i,k)
              vf(i,jm,k)=sbn(i,k)
            end if
          end do
        end do

        do k=1,kbm1
          do j=1,jm
            do i=1,im
              uf(i,j,k)=uf(i,j,k)*fsm(i,j)
              vf(i,j,k)=vf(i,j,k)*fsm(i,j)
            end do
          end do
        end do

        return

      else if(idx.eq.5) then

! vertical velocity boundary conditions
        do k=1,kbm1
          do j=1,jm
            do i=1,im
              w(i,j,k)=w(i,j,k)*fsm(i,j)
            end do
          end do
        end do

        return

      else if(idx.eq.6) then

! q2 and q2l boundary conditions

        do k=1,kb
          do j=1,jm
            ! west
            if(n_west.eq.-1) then
              u1=2.e0*u(2,j,k)*dti/(dx(1,j)+dx(2,j))
              if(u1.ge.0.e0) then
                uf(1,j,k)=q2(1,j,k)-u1*(q2(1,j,k)-small)
                vf(1,j,k)=q2l(1,j,k)-u1*(q2l(1,j,k)-small)
              else
                uf(1,j,k)=q2(1,j,k)-u1*(q2(2,j,k)-q2(1,j,k))
                vf(1,j,k)=q2l(1,j,k)-u1*(q2l(2,j,k)-q2l(1,j,k))
              end if
            end if

            ! east
            if(n_east.eq.-1) then
              u1=2.e0*u(im,j,k)*dti/(dx(im,j)+dx(imm1,j))
              if(u1.le.0.e0) then
                uf(im,j,k)=q2(im,j,k)-u1*(small-q2(im,j,k))
                vf(im,j,k)=q2l(im,j,k)-u1*(small-q2l(im,j,k))
              else
                uf(im,j,k)=q2(im,j,k)-u1*(q2(im,j,k)-q2(imm1,j,k))
                vf(im,j,k)=q2l(im,j,k)-u1*(q2l(im,j,k)-q2l(imm1,j,k))
              end if
            end if
          end do

          do i=1,im
            ! south
            if(n_south.eq.-1) then
              u1=2.e0*v(i,2,k)*dti/(dy(i,1)+dy(i,2))
              if(u1.ge.0.e0) then
                uf(i,1,k)=q2(i,1,k)-u1*(q2(i,1,k)-small)
                vf(i,1,k)=q2l(i,1,k)-u1*(q2l(i,1,k)-small)
              else
                uf(i,1,k)=q2(i,1,k)-u1*(q2(i,2,k)-q2(i,1,k))
                vf(i,1,k)=q2l(i,1,k)-u1*(q2l(i,2,k)-q2l(i,1,k))
              end if
            end if

            ! north
            if(n_north.eq.-1) then
              u1=2.e0*v(i,jm,k)*dti/(dy(i,jm)+dy(i,jmm1))
              if(u1.le.0.e0) then
                uf(i,jm,k)=q2(i,jm,k)-u1*(small-q2(i,jm,k))
                vf(i,jm,k)=q2l(i,jm,k)-u1*(small-q2l(i,jm,k))
              else
                uf(i,jm,k)=q2(i,jm,k)-u1*(q2(i,jm,k)-q2(i,jmm1,k))
                vf(i,jm,k)=q2l(i,jm,k)-u1*(q2l(i,jm,k)-q2l(i,jmm1,k))
              end if
            end if
          end do
        end do

        do k=1,kb
          do j=1,jm
            do i=1,im
              uf(i,j,k)=uf(i,j,k)*fsm(i,j)+1.e-10
              vf(i,j,k)=vf(i,j,k)*fsm(i,j)+1.e-10
            end do
          end do
        end do

        return

      endif

      end

!_______________________________________________________________________
      subroutine diag_w(idev)
! diagnose overshuting of vertical layers by vertical motion
! assuming layers thickness is dzz

      implicit none
      integer idev

      include 'pom.h'
      real dzw(im,jm,kb)
      real eps
      parameter(eps=1E-3)

      integer i,j,k,ii,jj

      integer iwmax,jwmax,iwmin,jwmin
      real wmax,wmin

      integer idifmax,jdifmax,kdifmax
      real wdifmax

      dzw(:,:,kb-1) = 0.0
      dzw(:,:,kb) = 0.0
      do k=1,kb-2
        dzw(:,:,k)=abs(w(:,:,k+1))*DTI*2
      enddo

      ii=0
      wdifmax=0.
      do k=1,kb-2
        do j=2,jmm1
          do i=2,imm1
            if(dzw(i,j,k).gt.dzz(k)*dt(i,j)) then
              ii=ii+1
              if(abs(dzw(i,j,k)-dzz(k)*dt(i,j)).gt.wdifmax) then
                wdifmax=abs(dzw(i,j,k)-dzz(k)*dt(i,j))
                idifmax=i
                jdifmax=j
                kdifmax=k
              end if
            end if
          end do
        end do
      end do
         
      jj=0
      wmax=0.
      wmin=0.
      do j=2,jmm1
        do i=2,imm1
          if(abs(w(i,j,kb)).gt.eps) then
            jj=jj+1
            if(w(i,j,kb).gt.wmax) then
              wmax=w(i,j,kb)
              iwmax=i
              jwmax=j
            end if
            if(w(i,j,kb).lt.wmin) then
              wmin=w(i,j,kb)
              iwmin=i
              jwmin=j
            end if
          end if
        end do
      end do

      if(ii.eq.0.and.jj.eq.0) then
        return
      end if

      if(jj.gt.0) then
        write(idev,*)
     $    "$$$ in diag_w, w_bottom >",eps," in N points: ",jj
        write(idev,*)
     $    'wmax i j ',wmax,iwmax,jwmax
        write(idev,*)
     $    'wmin i j ',wmin,iwmin,jwmin
      end if

      if(ii.gt.0) then
        write(idev,*)
     $    "$$$ in diag_w, dzz < w*dt*2 in N points: ",ii
        write(idev,*)
     $    'i j k',idifmax,jdifmax,kdifmax
        write(idev,*)'w: ',w(idifmax,jdifmax,kdifmax)
        write(idev,*)'dz: ',dzz(kdifmax)*dt(idifmax,jdifmax)
      endif

      end

      
