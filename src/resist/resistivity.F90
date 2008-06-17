#include "mhd.def"

module resistivity

! Written by R.K. Pawlaszek, September 2006 - December 2006
! Modified by M. Hanasz, February 2007
! Modified by R.K. Pawlaszek, June 2007

      use arrays
      use start
      use grid !, only: dxmn
      use constants, only: pi, small, big
      use mag_boundaries
      use mpi_setup
      
      real eta_max, dt_resist, dt_eint
      real eta_max_proc, eta_max_all
      integer loc_eta_max(3)

contains

      subroutine compute_resist(eta,ici)
      
      implicit none
      integer,intent(in)                    :: ici
      real, dimension(nx,ny,nz), intent(inout) :: eta
      real, dimension(:,:,:), allocatable   ::  wb
      allocate(wb(nx,ny,nz))
      
!--- square current computing in cell corner step by step

!--- current_z
        wb(:,:,:) = (b(iby,:,:,:)-cshift(b(iby,:,:,:),shift=-1,dim=xdim))/dl(xdim) &
                   -(b(ibx,:,:,:)-cshift(b(ibx,:,:,:),shift=-1,dim=ydim))/dl(ydim)

       eta(:,:,:) = (wb(:,:,:) &
             +cshift(wb(:,:,:),shift=-1,dim=zdim))**2/4.0

     if(dimensions .eq. '3d') then
!--- current_x
        wb(:,:,:) = (b(ibz,:,:,:)-cshift(b(ibz,:,:,:),shift=-1,dim=ydim))/dl(ydim) &
                   -(b(iby,:,:,:)-cshift(b(iby,:,:,:),shift=-1,dim=zdim))/dl(zdim)

       eta(:,:,:) = eta(:,:,:) + (wb(:,:,:) &
                          +cshift(wb(:,:,:),shift=-1,dim=xdim))**2/4.0
!--- current_y
        wb(:,:,:) = (b(ibx,:,:,:)-cshift(b(ibx,:,:,:),shift=-1,dim=zdim))/dl(zdim) &
                   -(b(ibz,:,:,:)-cshift(b(ibz,:,:,:),shift=-1,dim=xdim))/dl(xdim)

       eta(:,:,:) = eta(:,:,:) + (wb(:,:,:) &
                          +cshift(wb(:,:,:),shift=-1,dim=ydim))**2/4.0
     endif

!--- wb = current**2
        wb(:,:,:) = eta(:,:,:)

        eta(:,:,:) = eta_0 + eta_1*sqrt(max(0.0,eta(:,:,:)-j_crit**2))

        if(dimensions .eq. '3d') then
          where(eta > eta_0)
            eta(:,:,:) = (cshift(eta(:,:,:),shift=-1,dim=xdim) + cshift(eta(:,:,:),shift=1,dim=xdim) &
                         +cshift(eta(:,:,:),shift=-1,dim=ydim) + cshift(eta(:,:,:),shift=1,dim=ydim) &
                         +cshift(eta(:,:,:),shift=-1,dim=zdim) + cshift(eta(:,:,:),shift=1,dim=zdim) &
                         +dble(eta_scale)*eta(:,:,:))/(6.+dble(eta_scale))
          endwhere
        else
          where(eta > eta_0)
            eta(:,:,:) = (cshift(eta(:,:,:),shift=-1,dim=xdim) + cshift(eta(:,:,:),shift=1,dim=xdim) &
                         +cshift(eta(:,:,:),shift=-1,dim=ydim) + cshift(eta(:,:,:),shift=1,dim=ydim) &
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
        dt_eint = deint_max * abs(minval( &
                  ( u(iena,is:ie,js:je,ks:ke) &
                  - 0.5*(u(imxa,is:ie,js:je,ks:ke)**2  &
                        +u(imya,is:ie,js:je,ks:ke)**2  &
                        +u(imza,is:ie,js:je,ks:ke)**2) &
                        /u(idna,is:ie,js:je,ks:ke)     &
                - 0.5*(b(ibx,is:ie,js:je,ks:ke)**2  &
                      +b(iby,is:ie,js:je,ks:ke)**2  &
                      +b(ibz,is:ie,js:je,ks:ke)**2))&
                        /(eta(is:ie,js:je,ks:ke)   &
                        *wb(is:ie,js:je,ks:ke)+small))) 
#endif /* ISO */
        
      deallocate(wb)

      if (ici .eq. icz) then
         eta(:,:,:)=0.5*(eta(:,:,:)+cshift(eta(:,:,:),shift=1,dim=zdim))
      elseif(ici .eq. icy) then
         eta(:,:,:)=0.5*(eta(:,:,:)+cshift(eta(:,:,:),shift=1,dim=ydim))
      elseif(ici .eq. icx) then
         eta(:,:,:)=0.5*(eta(:,:,:)+cshift(eta(:,:,:),shift=1,dim=xdim))
      endif

      return

      end subroutine compute_resist

!-----------------------------------------------------------------------

      subroutine timestep_resist

        use constants, only : big

        implicit none
        real dx2,dt_resist_all

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

    implicit none
#ifdef SPLIT
#ifdef ORIG
    real,dimension(nx,ny,nz) :: b1
#endif /* ORIG */
#endif /* SPLIT */

    real                     :: di
    real,dimension(:,:,:), allocatable :: w,wm,wp,dw,eta
    integer                  :: ibi,ici,n

    allocate(w(nx,ny,nz), wm(nx,ny,nz), wp(nx,ny,nz), dw(nx,ny,nz), eta(nx,ny,nz))

    di = dl(n)
    eta = 0.0

    call compute_resist(eta,ici)
#ifdef SPLIT
#ifdef ORIG

! HALF STEP
    w(:,:,:) = (b(ibi,:,:,:)-cshift(b(ibi,:,:,:),shift=-1,dim=n))/di
    w  = eta*w
    b1 = b(ibi,:,:,:)+(cshift(w,shift=1,dim=n)-w)*dt/di/2.

! FULL STEP
    w(:,:,:) = (b1(:,:,:)-cshift(b1(:,:,:),shift=-1,dim=n))/di
    w  = eta*w
    wp = (cshift(w,shift=1,dim=n)-w)/2.
    wm = (w-cshift(w,shift=-1,dim=n))/2.
    dw = 0.
    where(wm*wp > 0.) dw=2.*wm*wp/(wm+wp)
    wcu = (w+dw)*dt
#endif /* ORIG */
#ifdef SSP

    w(:,:,:) = (b(ibi,:,:,:)-cshift(b(ibi,:,:,:),shift=-1,dim=n))/di
    w  = eta*w
    wp = (cshift(w,shift=1,dim=n)-w)/2.
    wm = (w-cshift(w,shift=-1,dim=n))/2.
    dw = 0.
    where(wm*wp > 0.) dw=2.*wm*wp/(wm+wp)
    wcu = (w+dw)*dt
#endif /* SSP */
#else /* SPLIT */

    w(:,:,:) = (b(ibi,:,:,:)-cshift(b(ibi,:,:,:),shift=-1,dim=n))/di
    w  = eta*w
#ifdef ORIG

    if(istep .eq. 1) then
      wa = w*dt
    elseif(istep .eq. 2) then
      wp = (cshift(w,shift=1,dim=n)-w)/2.
      wm = (w-cshift(w,shift=-1,dim=n))/2.
      dw = 0.
      where(wm*wp > 0.) dw=2.*wm*wp/(wm+wp)
      wa = (w+dw)*dt
    else
      write(*,*) 'That high integration order not implemented'
      stop
    endif
#endif /* ORIG */
#ifdef SSP
    wp = (cshift(w,shift=1,dim=n)-w)/2.
    wm = (w-cshift(w,shift=-1,dim=n))/2.
    dw = 0.
    where(wm*wp > 0.) dw=2.*wm*wp/(wm+wp)
    wa = (w+dw)*dt
#endif /* SSP */
#endif /* SPLIT */
    deallocate(w,wm,wp,dw,eta)
  end subroutine tvdd

!-------------------------------------------------------------------------------

  subroutine diffuseby_x
  
    call tvdd(iby,icz,xdim)
#ifdef SPLIT

    call bnd_emf(wcu,'emfz','xdim')
    call bnd_emf(wcu,'emfz','ydim')
    if(dimensions .eq. '3d') call bnd_emf(wcu,'emfz','zdim')
#else /* SPLIT */

    call bnd_emf(wa,'emfz','xdim')
    call bnd_emf(wa,'emfz','ydim')
    if(dimensions .eq. '3d') call bnd_emf(wa,'emfz','zdim')

    Lb(iby,:,:,:) = Lb(iby,:,:,:) - wa/dl(xdim)
    wa = cshift(wa,shift=1,dim=xdim)
    Lb(iby,:,:,:) = Lb(iby,:,:,:) + wa/dl(xdim)
    wa = cshift(wa,shift=-1,dim=xdim)
    Lb(ibx,:,:,:) = Lb(ibx,:,:,:) + wa/dl(ydim)
    wa = cshift(wa,shift=1,dim=ydim)
    Lb(ibx,:,:,:) = Lb(ibx,:,:,:) - wa/dl(ydim)
#endif /* SPLIT */

  end subroutine diffuseby_x

  subroutine diffusebz_x

    call tvdd(ibz,icy,xdim)
#ifdef SPLIT

    call bnd_emf(wcu,'emfy','xdim')
    call bnd_emf(wcu,'emfy','ydim')
    if(dimensions .eq. '3d') call bnd_emf(wcu,'emfy','zdim')
#else /* SPLIT */

    call bnd_emf(wa, 'emfy', 'xdim')
    call bnd_emf(wa, 'emfy', 'ydim')
    if(dimensions .eq. '3d') call bnd_emf(wa, 'emfy', 'zdim')

    Lb(ibz,:,:,:) = Lb(ibz,:,:,:) - wa/dl(xdim)
    wa = cshift(wa,shift=1,dim=xdim)
    Lb(ibz,:,:,:) = Lb(ibz,:,:,:) + wa/dl(xdim)
    wa = cshift(wa,shift=-1,dim=xdim)
    Lb(ibx,:,:,:) = Lb(ibx,:,:,:) + wa/dl(zdim)
    wa = cshift(wa,shift=1,dim=zdim)
    Lb(ibx,:,:,:) = Lb(ibx,:,:,:) - wa/dl(zdim)
#endif /* SPLIT */

  end subroutine diffusebz_x

  subroutine diffusebz_y

    call tvdd(ibz,icx,ydim)
#ifdef SPLIT

    call bnd_emf(wcu,'emfx','ydim')
    if(dimensions .eq. '3d') call bnd_emf(wcu,'emfx','zdim')
    call bnd_emf(wcu,'emfx','xdim')
#else /* SPLIT */

    call bnd_emf(wa, 'emfx', 'ydim')
    if(dimensions .eq. '3d') call bnd_emf(wa, 'emfx', 'zdim')
    call bnd_emf(wa, 'emfx', 'xdim')
  
    Lb(ibz,:,:,:) = Lb(ibz,:,:,:) - wa/dl(ydim)
    wa = cshift(wa,shift=1,dim=ydim)
    Lb(ibz,:,:,:) = Lb(ibz,:,:,:) + wa/dl(ydim)
    wa = cshift(wa,shift=-1,dim=ydim)
    Lb(iby,:,:,:) = Lb(iby,:,:,:) + wa/dl(zdim)
    wa = cshift(wa,shift=1,dim=zdim)
    Lb(iby,:,:,:) = Lb(iby,:,:,:) - wa/dl(zdim)
#endif /* SPLIT */
  end subroutine diffusebz_y

  subroutine diffusebx_y

    call tvdd(ibx,icz,ydim)
#ifdef SPLIT

    call bnd_emf(wcu, 'emfz', 'ydim')
    if(dimensions .eq. '3d') call bnd_emf(wcu, 'emfz', 'zdim')
    call bnd_emf(wcu, 'emfz', 'xdim')
#else /* SPLIT */

    call bnd_emf(wa, 'emfz', 'ydim')
    if(dimensions .eq. '3d') call bnd_emf(wa, 'emfz', 'zdim')
    call bnd_emf(wa, 'emfz', 'xdim')

    Lb(ibx,:,:,:) = Lb(ibx,:,:,:) - wa/dl(ydim)
    wa = cshift(wa,shift=1,dim=ydim)
    Lb(ibx,:,:,:) = Lb(ibx,:,:,:) + wa/dl(ydim)
    wa = cshift(wa,shift=-1,dim=ydim)
    Lb(iby,:,:,:) = Lb(iby,:,:,:) + wa/dl(xdim)
    wa = cshift(wa,shift=1,dim=xdim)
    Lb(iby,:,:,:) = Lb(iby,:,:,:) - wa/dl(xdim)
#endif /* SPLIT */
  end subroutine diffusebx_y

  subroutine diffusebx_z

    call tvdd(ibx,icy,zdim)
#ifdef SPLIT

    call bnd_emf(wcu, 'emfy', 'zdim')
    call bnd_emf(wcu, 'emfy', 'xdim')
    call bnd_emf(wcu, 'emfy', 'ydim')
#else /* SPLIT */

    call bnd_emf(wa, 'emfy', 'zdim')
    call bnd_emf(wa, 'emfy', 'xdim')
    call bnd_emf(wa, 'emfy', 'ydim')

    Lb(ibx,:,:,:) = Lb(ibx,:,:,:) - wa/dl(zdim)
    wa = cshift(wa,shift=1,dim=zdim)
    Lb(ibx,:,:,:) = Lb(ibx,:,:,:) + wa/dl(zdim)
    wa = cshift(wa,shift=-1,dim=zdim)
    Lb(ibz,:,:,:) = Lb(ibz,:,:,:) + wa/dl(xdim)
    wa = cshift(wa,shift=1,dim=xdim)
    Lb(ibz,:,:,:) = Lb(ibz,:,:,:) - wa/dl(xdim)
#endif /* SPLIT */
  end subroutine diffusebx_z

  subroutine diffuseby_z
    
    call tvdd(iby,icx,zdim)
#ifdef SPLIT

    call bnd_emf(wcu, 'emfx', 'zdim')
    call bnd_emf(wcu, 'emfx', 'xdim')
    call bnd_emf(wcu, 'emfx', 'ydim')
#else /* SPLIT */

    call bnd_emf(wa, 'emfx', 'zdim')
    call bnd_emf(wa, 'emfx', 'xdim')
    call bnd_emf(wa, 'emfx', 'ydim')

    Lb(iby,:,:,:) = Lb(iby,:,:,:) - wa/dl(zdim)
    wa = cshift(wa,shift=1,dim=zdim)
    Lb(iby,:,:,:) = Lb(iby,:,:,:) + wa/dl(zdim)
    wa = cshift(wa,shift=-1,dim=zdim)
    Lb(ibz,:,:,:) = Lb(ibz,:,:,:) + wa/dl(ydim)
    wa = cshift(wa,shift=1,dim=ydim)
    Lb(ibz,:,:,:) = Lb(ibz,:,:,:) - wa/dl(ydim)
#endif /* SPLIT */
  end subroutine diffuseby_z

end module resistivity
