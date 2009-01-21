! $Id$
#include "piernik.def"

module resistivity

! Written by R.K. Pawlaszek, September 2006 - December 2006
! Modified by M. Hanasz, February 2007
! Modified by R.K. Pawlaszek, June 2007

!      use arrays
!      use grid !, only: dxmn
      use constants, only: pi, small, big
!      use mag_boundaries
      use mpi_setup
      use fluidindex, only : ibx,iby,ibz,icx,icy,icz
      use initionized, only : idni,imxi,imyi,imzi
#ifndef ISO
      use initionized, only : ieni
#endif /* ISO */
      implicit none
      real    :: cfl_resist, eta_0, eta_1, j_crit, deint_max
      integer :: eta_scale

      real :: eta_max, dt_resist, dt_eint
      real :: eta_max_proc, eta_max_all
      integer, dimension(3) :: loc_eta_max

contains

      subroutine init_resistivity
         use mpi_setup
         implicit none
         character(LEN=100) :: par_file, tmp_log_file

         namelist /RESISTIVITY/ cfl_resist, eta_0, eta_1, eta_scale, j_crit, deint_max
         
         par_file = trim(cwd)//'/problem.par'
         tmp_log_file = trim(cwd)//'/tmp.log'

         cfl_resist  =  0.4
         eta_0       =  0.0
         eta_1       =  0.0
         eta_scale   =  4
         j_crit      =  1.0e6
         deint_max   = 0.01
       
         if(proc .eq. 0) then
            open(1,file=par_file)
               read(unit=1,nml=RESISTIVITY)
            close(1)
            open(3, file=tmp_log_file, position='append')
               write(unit=3,nml=RESISTIVITY)
            close(3)

            ibuff(1) = eta_scale

            rbuff(1) = cfl_resist
            rbuff(2) = eta_0
            rbuff(3) = eta_1
            rbuff(4) = j_crit
            rbuff(5) = deint_max

            call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
            call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

         else

            call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
            call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

            eta_scale          = ibuff(1)

            cfl_resist         = rbuff(1)
            eta_0              = rbuff(2)
            eta_1              = rbuff(3)
            j_crit             = rbuff(4)
            deint_max          = rbuff(5)
         endif
         if(eta_scale < 0) then
            write(*,*) 'eta_scale must be greater or equal 0'
            call mpistop
         endif


      end subroutine init_resistivity

      subroutine compute_resist(eta,ici)
       use arrays, only : b,u
       use grid,   only : dl,xdim,ydim,zdim,nx,ny,nz,is,ie,js,je,ks,ke,nzd
       use func,   only : mshift, pshift

      implicit none
      integer,intent(in)                    :: ici
      real, dimension(nx,ny,nz), intent(inout) :: eta
      real, dimension(:,:,:), allocatable   ::  wb,etahelp
      allocate(wb(nx,ny,nz),etahelp(nx,ny,nz))

!--- square current computing in cell corner step by step

!--- current_z
!       wb(:,:,:) = (b(iby,:,:,:)-mshift(b(iby,:,:,:),xdim))/dl(xdim) &
!                   -(b(ibx,:,:,:)-mshift(b(ibx,:,:,:),ydim))/dl(ydim)
       wb(:,:,:) = (b(iby,:,:,:)-mshift(b(iby,:,:,:),xdim))/dl(xdim)
       wb = wb -   (b(ibx,:,:,:)-mshift(b(ibx,:,:,:),ydim))/dl(ydim)

       eta(:,:,:) = 0.25*( wb(:,:,:) + mshift(wb(:,:,:),zdim) )**2

     if(nzd /=1) then
!--- current_x
!       wb(:,:,:) = (b(ibz,:,:,:)-mshift(b(ibz,:,:,:),ydim))/dl(ydim) &
!                   -(b(iby,:,:,:)-mshift(b(iby,:,:,:),zdim))/dl(zdim)
       wb(:,:,:) = (b(ibz,:,:,:)-mshift(b(ibz,:,:,:),ydim))/dl(ydim)
       wb = wb -   (b(iby,:,:,:)-mshift(b(iby,:,:,:),zdim))/dl(zdim)

       eta(:,:,:) = eta(:,:,:) + 0.25*( wb(:,:,:)+mshift(wb(:,:,:),xdim) )**2
!--- current_y
!       wb(:,:,:) = (b(ibx,:,:,:)-mshift(b(ibx,:,:,:),zdim))/dl(zdim) &
!                   -(b(ibz,:,:,:)-mshift(b(ibz,:,:,:),xdim))/dl(xdim)
       wb(:,:,:) = (b(ibx,:,:,:)-mshift(b(ibx,:,:,:),zdim))/dl(zdim)
       wb = wb -   (b(ibz,:,:,:)-mshift(b(ibz,:,:,:),xdim))/dl(xdim)

       eta(:,:,:) = eta(:,:,:) +0.25*( wb(:,:,:) + mshift(wb(:,:,:),ydim))**2
     endif

!--- wb = current**2
        wb(:,:,:) = eta(:,:,:)

        eta(:,:,:) = eta_0 + eta_1*sqrt(max(0.0,eta(:,:,:)- j_crit**2 ))
!>
!! \todo Following lines are splitted into separate lines because of intel and gnu dbgs
!! shoud that be so? Is there any other solution instead splitting?
!<
        if(nzd /= 1) then
          etahelp(:,:,:)  =   mshift(eta(:,:,:),xdim)
          etahelp = etahelp + pshift(eta(:,:,:),xdim)
          etahelp = etahelp + mshift(eta(:,:,:),ydim)
          etahelp = etahelp + pshift(eta(:,:,:),ydim)
          etahelp = etahelp + mshift(eta(:,:,:),zdim)
          etahelp = etahelp + pshift(eta(:,:,:),zdim)
          etahelp = (etahelp+dble(eta_scale)*eta(:,:,:))/(6.+dble(eta_scale))
          where(eta > eta_0)
             eta = etahelp
          endwhere
        else
          where(eta > eta_0)
            eta(:,:,:) = (mshift(eta(:,:,:),xdim) + pshift(eta(:,:,:),xdim) &
                         +mshift(eta(:,:,:),ydim) + pshift(eta(:,:,:),ydim) &
                         +dble(eta_scale)*eta(:,:,:))/(4.+dble(eta_scale))
          endwhere
        endif

        eta_max_proc      = maxval(eta(is:ie,js:je,ks:ke))
        loc_eta_max       = maxloc(eta(is:ie,js:je,ks:ke))

        eta_max = eta_max_proc

        call MPI_REDUCE(eta_max_proc, eta_max_all, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, ierr)
        call MPI_BCAST (eta_max_all, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)

        eta_max = eta_max_all

#ifndef ISO
        dt_eint = deint_max * abs(minval(               &
                  ( u(ieni,is:ie,js:je,ks:ke)           &
                  - 0.5*( u(imxi,is:ie,js:je,ks:ke)**2  &
                        + u(imyi,is:ie,js:je,ks:ke)**2  &
                        + u(imzi,is:ie,js:je,ks:ke)**2 )&
                        /u(idni,is:ie,js:je,ks:ke)      &
                  - 0.5*( b(ibx,is:ie,js:je,ks:ke)**2   &
                        + b(iby,is:ie,js:je,ks:ke)**2   &
                        + b(ibz,is:ie,js:je,ks:ke)**2)) &
                        /( eta(is:ie,js:je,ks:ke)       &
                        *wb(is:ie,js:je,ks:ke)+small) ))
#endif /* ISO */

      deallocate(wb,etahelp)

      if (ici .eq. icz) then
         eta(:,:,:)=0.5*(eta(:,:,:)+pshift(eta(:,:,:),zdim))
      elseif(ici .eq. icy) then
         eta(:,:,:)=0.5*(eta(:,:,:)+pshift(eta(:,:,:),ydim))
      elseif(ici .eq. icx) then
         eta(:,:,:)=0.5*(eta(:,:,:)+pshift(eta(:,:,:),xdim))
      endif

      return

      end subroutine compute_resist

!-----------------------------------------------------------------------

      subroutine timestep_resist
        use grid, only : dxmn
        use constants, only : big

        implicit none
        real dx2,dt_resist_min

        if(eta_max .ne. 0.) then
          dx2 = dxmn**2
          dt_resist = cfl_resist*dx2/(2.*eta_max)
#ifndef ISO
          dt_resist = min(dt_resist,dt_eint)
#endif /* ISO */
        else
          dt_resist = big
        endif

        call MPI_REDUCE(dt_resist, dt_resist_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, comm, ierr)
        call MPI_BCAST (dt_resist_min, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)

        dt_resist = dt_resist_min

      end subroutine timestep_resist

!-----------------------------------------------------------------------------

  subroutine tvdd(ibi,ici,n)
    use arrays, only : b,wcu
    use start,  only : dt
    use func,   only : mshift, pshift
    use grid,   only : nx,ny,nz,dxmn,dl

    implicit none
    real                     :: di
    real,dimension(:,:,:), allocatable :: w,wm,wp,dw,eta,b1
    integer                  :: ibi,ici,n

    allocate(w(nx,ny,nz), wm(nx,ny,nz), wp(nx,ny,nz), dw(nx,ny,nz), eta(nx,ny,nz))
    allocate(b1(nx,ny,nz))

    di = dl(n)
    eta = 0.0

    call compute_resist(eta,ici)

! HALF STEP
    w(:,:,:) = (b(ibi,:,:,:)-mshift(b(ibi,:,:,:),n))/di
    w  = eta*w
    b1 = b(ibi,:,:,:)+0.5*(pshift(w,n)-w)*dt/di

! FULL STEP
    w(:,:,:) = (b1(:,:,:)-mshift(b1(:,:,:),n))/di
    w  = eta*w
    wp = 0.5*(pshift(w,n)-w)
    wm = 0.5*(w-mshift(w,n))
    dw = 0.
    where(wm*wp > 0.) dw=2.*wm*wp/(wm+wp)
    wcu = (w+dw)*dt
    deallocate(w,wm,wp,dw,eta,b1)
  end subroutine tvdd

!-------------------------------------------------------------------------------

  subroutine diffuseby_x
    use func, only : mshift, pshift
    use grid, only : xdim,nzd
    use fluidindex, only : iby,icz
    use arrays, only : wcu
    use mag_boundaries, only : bnd_emf
    implicit none

    call tvdd(iby,icz,xdim)
    call bnd_emf(wcu,'emfz','xdim')
    call bnd_emf(wcu,'emfz','ydim')
    if(nzd /= 1) call bnd_emf(wcu,'emfz','zdim')

  end subroutine diffuseby_x

  subroutine diffusebz_x
    use func, only : mshift, pshift
    use grid, only : xdim,nzd
    use fluidindex, only : ibz,icy
    use arrays, only : wcu
    use mag_boundaries, only : bnd_emf
    implicit none

    call tvdd(ibz,icy,xdim)
    call bnd_emf(wcu,'emfy','xdim')
    call bnd_emf(wcu,'emfy','ydim')
    if(nzd /= 1) call bnd_emf(wcu,'emfy','zdim')

  end subroutine diffusebz_x

  subroutine diffusebz_y
    use func, only : mshift, pshift
    use grid, only : ydim,nzd
    use fluidindex, only : ibz,icx
    use arrays, only : wcu
    use mag_boundaries, only : bnd_emf
    implicit none

    call tvdd(ibz,icx,ydim)
    call bnd_emf(wcu,'emfx','ydim')
    if(nzd /= 1) call bnd_emf(wcu,'emfx','zdim')
    call bnd_emf(wcu,'emfx','xdim')
  end subroutine diffusebz_y

  subroutine diffusebx_y
    use func, only : mshift, pshift
    use grid, only : ydim,nzd
    use fluidindex, only : ibx,icz
    use arrays, only : wcu
    use mag_boundaries, only : bnd_emf
    implicit none

    call tvdd(ibx,icz,ydim)
    call bnd_emf(wcu, 'emfz', 'ydim')
    if(nzd /= 1) call bnd_emf(wcu, 'emfz', 'zdim')
    call bnd_emf(wcu, 'emfz', 'xdim')
  end subroutine diffusebx_y

  subroutine diffusebx_z
    use func, only : mshift, pshift
    use grid, only : zdim,nzd
    use fluidindex, only : ibx,icy
    use arrays, only : wcu
    use mag_boundaries, only : bnd_emf
    implicit none

    call tvdd(ibx,icy,zdim)
    if(nzd /= 1) call bnd_emf(wcu, 'emfy', 'zdim')
    call bnd_emf(wcu, 'emfy', 'xdim')
    call bnd_emf(wcu, 'emfy', 'ydim')
  end subroutine diffusebx_z

  subroutine diffuseby_z
    use func, only : mshift, pshift
    use grid, only : zdim,nzd
    use fluidindex, only : iby,icx
    use arrays, only : wcu
    use mag_boundaries, only : bnd_emf
    implicit none

    call tvdd(iby,icx,zdim)
    if(nzd /= 1) call bnd_emf(wcu, 'emfx', 'zdim')
    call bnd_emf(wcu, 'emfx', 'xdim')
    call bnd_emf(wcu, 'emfx', 'ydim')
  end subroutine diffuseby_z

end module resistivity
