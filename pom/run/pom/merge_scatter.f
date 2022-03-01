c
c  merge2d() merges a 2-D array into a 2-D array in rank0 process.
c
      subroutine merge2d(dest,src)

      implicit none
      include 'mpif.h'
      include 'pom.h'
      real    dest(im_global,jm_global),src(im,jm)
c     --
      integer nblk,ijloc,iloc,jloc,iglb,jglb,m,n,ierr
      real work(im*jm*n_proc)

      nblk=im*jm
      call MPI_Gather(src,nblk,MPI_REAL,
     $                work,nblk,MPI_REAL,
     $                master_task,pom_comm,ierr)

      if(my_task.eq.master_task) then

         do m=1,nblk*n_proc
            
            n=int(m/nblk)

            ijloc=mod(m,nblk)
            if(ijloc.eq.0) then
               n=n-1
               ijloc=nblk
            end if
            
            jloc=int(ijloc/im)+1
            iloc=mod(ijloc,im)
            if(iloc.eq.0) then
               jloc=jloc-1
               iloc=im
            end if
            
            iglb=iloc+mod(n,nproc_x)*(im-2)
            jglb=jloc+(n/nproc_x)*(jm-2)
            
!     Shun Ohishi 2018.07
!            if(iloc == 1 .and. 0 < mod(n,nproc_x)) cycle
!            if(iloc == im .and. mod(n,nproc_x) < nproc_x-1) cycle
!            if(jloc == 1 .and. n/nproc_x > 0) cycle
!            if(jloc == jm .and. n/nproc_x < nproc_y-1)cycle

            dest(iglb,jglb)=work(m)
            
         end do
         
      endif
      
      call MPI_Barrier(pom_comm,ierr)

      return
      end

c
c  merge3d() merges a 3D array (src) in each process into
c  a array (dest) in rank0 process.
c
      subroutine merge3d(dest,src)

      implicit none
      include 'mpif.h'
      include 'pom.h'
      real    dest(im_global,jm_global,kb),src(im,jm,kb)
c     --
      integer nblk,kloc,ijkloc,ijloc,iloc,jloc,iglb,jglb,m,n,ierr
      real work(im*jm*kb*n_proc)

      nblk=im*jm*kb
      call MPI_Gather(src,nblk,MPI_REAL,
     $     work,nblk,MPI_REAL,
     $     master_task,pom_comm,ierr)
      
      if(my_task.eq.master_task) then
         
         do m=1,nblk*n_proc

            n=int(m/nblk)

            ijkloc=mod(m,nblk)
            if(ijkloc.eq.0) then
               n=n-1
               ijkloc=nblk
            end if

            kloc=int(ijkloc/(im*jm))+1
            ijloc=mod(ijkloc,im*jm)

            if(ijloc.eq.0) then
               kloc=kloc-1
               ijloc=im*jm
            end if

            jloc=int(ijloc/im)+1

            iloc=mod(ijloc,im)
            if(iloc.eq.0) then
               jloc=jloc-1
               iloc=im
            end if

            iglb=iloc+mod(n,nproc_x)*(im-2)
            jglb=jloc+(n/nproc_x)*(jm-2)

!     Shun Ohishi 2018.07            
!            if(iloc == 1 .and. 0 < mod(n,nproc_x)) cycle
!            if(iloc == im .and. mod(n,nproc_x) < nproc_x-1) cycle
!            if(jloc == 1 .and. n/nproc_x > 0) cycle
!            if(jloc == jm .and. n/nproc_x < nproc_y-1)cycle
            
            dest(iglb,jglb,kloc)=work(m)
!     if(kloc.eq.1) write(6,*) m,n,iglb,jglb
         end do
         
      endif
      
      call MPI_Barrier(pom_comm,ierr)

      return
      end

c
c  scatter2d() scatter a global 2-D array into a local 2-D array
c
      subroutine scatter2d(dest,src)

      implicit none
      include 'mpif.h'
      include 'pom.h'
      real    dest(im,jm),src(im_global,jm_global)
c     --
      integer iloc,jloc,i,j
!      character*1 cha_my_task
!      write(cha_my_task,'(i1.1)') my_task
!      open(88,file=cha_my_task//'/scatter2d',form='formatted')

      do jloc=1,jm
        do iloc=1,im
          dest(iloc,jloc)=src(i_global(iloc),j_global(jloc))
!          write(88,*) 'g->l ',i_global(iloc),j_global(jloc),iloc,jloc
        end do
      end do

!      close(88)

      return
      end

c
c  scatter3d() scatter a global 3-D array into a local 3-D array
c
      subroutine scatter3d(dest,src)

      implicit none
      include 'mpif.h'
      include 'pom.h'
      real    dest(im,jm,kb),src(im_global,jm_global,kb)
c     --
      integer iloc,jloc,i,j,k

      do k=1,kb
        do jloc=1,jm
          do iloc=1,im
            dest(iloc,jloc,k)=src(i_global(iloc),j_global(jloc),k)
          end do
        end do
      end do

      return
      end

c
c  scatter_bnd_ew() 
c    scatter a global 2-D/3-D array into 
c    a eastern/western boundary of local 2-D array
c
      subroutine scatter_bnd_ew(dest,src,kmax)
      implicit none
      include 'mpif.h'
      include 'pom.h'
      integer kmax
      real    dest(jm,kmax),src(jm_global,kmax)
c     --
      integer jloc,j,k

      if(n_west.eq.-1.or.n_east.eq.-1) then
        do k=1,kmax
          do jloc=1,jm
            dest(jloc,k)=src(j_global(jloc),k)
          end do
        end do
      end if

      return
      end

c
c  scatter_bnd_ns() 
c    scatter a global 2-D/3-D array into 
c    a northern/southern boundary of local 2-D array
c
      subroutine scatter_bnd_ns(dest,src,kmax)
      implicit none
      include 'mpif.h'
      include 'pom.h'
      integer kmax
      real    dest(im,kmax),src(im_global,kmax)
c     --
      integer iloc,i,k

      if(n_north.eq.-1.or.n_south.eq.-1) then
        do k=1,kmax
          do iloc=1,im
            dest(iloc,k)=src(i_global(iloc),k)
          end do
        end do
      end if

      return
      end

