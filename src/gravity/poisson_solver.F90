#include "mhd.def"
module poisson_solver

! Written by G. Kowal
! Adopted for this code by R. Maliszewski & M. Hanasz, May 2006

!  use start
!  use arrays
!  use grid
  use constants
  implicit none

contains

!!=====================================================================
!!
!! SUBROUTINE POISSON: solves Poisson equation
!!
  subroutine poisson
    use start, only  : bnd_xl, bnd_xr, bnd_yl, bnd_yr, nb, nxd, nyd, dz
    use arrays, only : idna 
    implicit none
  
    if( bnd_xl .eq. 'per' .and. bnd_xr .eq. 'per' .and. &
        bnd_yl .eq. 'per' .and. bnd_yr .eq. 'per'      ) then

         call poisson_xyp(u(idna,nb+1:nb+nxd,nb+1:nb+nyd,:), &
                           gp(nb+1:nb+nxd,nb+1:nb+nyd,:),dz)
  
! Boundary conditions

         gp(1:nb,:,:)              = gp(nxd+1:nxd+nb,:,:)
         gp(nxd+nb+1:nxd+2*nb,:,:) = gp(nb+1:2*nb,:,:)
         gp(:,1:nb,:)              = gp(:,nyd+1:nyd+nb,:)
         gp(:,nyd+nb+1:nyd+2*nb,:) = gp(:,nb+1:2*nb,:)
    else
      write(*,*) 'POISSON SOLVER: not implemented for boundary conditions'
      write(*,*) '                xdim: ',bnd_xl,',  ', bnd_xr
      write(*,*) '                ydim: ',bnd_yl,',  ', bnd_yr
      stop
    endif

  end subroutine poisson

!!=====================================================================
!!
!! SUBROUTINE POISSON_XYP: solves Poisson equation for periodic
!! bnd conditions in X and Y, and arbitrary bnd in Z
!!
  subroutine poisson_xyp(den, pot, dz)

    implicit none

    real, dimension(:,:,:), intent(in)  :: den
    real, dimension(:,:,:), intent(out) :: pot
    real, optional        , intent(in)  :: dz

    complex(kind=8), dimension(:,:,:), allocatable :: fft
    complex(kind=8), dimension(:,:)  , allocatable :: ctmp
    complex(kind=8), dimension(:)    , allocatable :: ztmp
    complex(kind=8), dimension(:)    , allocatable :: vl, vd, vu, vs

    real(kind=8)   , dimension(:,:)  , allocatable :: rtmp
    real(kind=8)   , dimension(:)    , allocatable :: kx, ky

    integer        , dimension(:)    , allocatable :: ipiv

    real(kind=8)    :: lambda, factor

    integer(kind=8) :: planf, planb

    integer         :: nx, ny, nz, np, p, q, k, info

    real(kind=8), parameter :: dpi = 6.2831853071795865d0
    integer     , parameter :: FFTW_ESTIMATE = 64
!
!----------------------------------------------------------------------
!

! determine input array dimension
!
  nx = size(den, 1)
  ny = size(den, 2)
  nz = size(den, 3)

! compute dimensions for complex arrays
!
    np = nx / 2 + 1

! allocate arrays
!
    allocate(fft(np,ny,nz))
    allocate(ctmp(np,ny))
    allocate(rtmp(nx,ny))
    allocate(ztmp(nz))
    allocate(kx(np))
    allocate(ky(ny))
    allocate(vd(nz))
    allocate(vl(nz-1))
    allocate(vu(nz-1))
    allocate(vs(nz-2))
    allocate(ipiv(nz))

! create plan for the forward FFT
!
    call dfftw_plan_dft_r2c_2d(planf, nx, ny, rtmp, ctmp, FFTW_ESTIMATE)

! perform forward FFT for each k
!
    factor = 1.0

    if (present(dz) .eq. .true.) factor = dz * dz

    do k = 1, nz
      rtmp(:,:)   = factor * den(:,:,k)

      call dfftw_execute(planf)

      fft(:,:,k) = ctmp(:,:) / real(nx * ny)
    enddo

! destroy plan for the forward FFT
!
    call dfftw_destroy_plan(planf)

! compute wave vectors
!
    do p = 1, np
      kx(p) = dpi * (p - 1) / nx
    enddo

    ky(1) = 0.0
    do q = 2, ny / 2 + 1
      ky(q) = dpi * (q - 1) / ny
      ky(ny+2-q) = ky(q)
    enddo

! compute eigenvalues for each p and q and solve linear system
!
    do q = 1, ny
      do p = 1, np

        lambda = - kx(p) * kx(p) - ky(q) * ky(q) - 2.0

        vd(:)  = lambda
        vl(:)  =    1.0
        vu(:)  =    1.0

! factorize the matrix
!
        call zgttrf(nz, vl, vd, vu, vs, ipiv, info)

! solve the system
!
        ztmp(:) = fft(p,q,:)

        call zgttrs('N', nz, 1, vl, vd, vu, vs, ipiv, ztmp, nz, info)

        fft(p,q,:) = ztmp(:)

      enddo
    enddo

! create plan for the backward FFT
!
    call dfftw_plan_dft_c2r_2d(planb, nx, ny, ctmp, rtmp, FFTW_ESTIMATE)

! perform inverse FFT for each k
!
    do k = 1, nz
      ctmp(:,:)   = fft(:,:,k)
      call dfftw_execute(planb)
      pot(:,:,k) = fpiG*rtmp(:,:)
    enddo

! destroy plan for the backward FFT
!
    call dfftw_destroy_plan(planb)

! deallocate rhofft
!
    if (allocated(ipiv))     deallocate(ipiv)
    if (allocated(vl))     deallocate(vl)
    if (allocated(vd))     deallocate(vd)
    if (allocated(vu))     deallocate(vu)
    if (allocated(vs))     deallocate(vs)
    if (allocated(kx))     deallocate(kx)
    if (allocated(ky))     deallocate(ky)
    if (allocated(ztmp))   deallocate(ztmp)
    if (allocated(rtmp))   deallocate(rtmp)
    if (allocated(ctmp))   deallocate(ctmp)
    if (allocated(fft))    deallocate(fft)

  end subroutine poisson_xyp

end module poisson_solver
