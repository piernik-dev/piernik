! $Id$
!
! PIERNIK Code Copyright (C) 2006 Michal Hanasz
!
!    This file is part of PIERNIK code.
!
!    PIERNIK is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    PIERNIK is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with PIERNIK.  If not, see <http://www.gnu.org/licenses/>.
!
!    Initial implementation of PIERNIK code was based on TVD split MHD code by
!    Ue-Li Pen
!        see: Pen, Arras & Wong (2003) for algorithm and
!             http://www.cita.utoronto.ca/~pen/MHD
!             for original source code "mhd.f90"
!
!    For full list of developers see $PIERNIK_HOME/license/pdt.txt
!
#include "piernik.h"
module poissonsolver
! pulled by POISSON_FFT
#ifdef POISSON_FFT

  implicit none

  private
  public :: poisson_solve

contains

!!=====================================================================
!!
!! SUBROUTINE POISSON: solves Poisson equation
!!
  subroutine poisson_solve(dens)

    use arrays,         only: u, sgp
    use dataio_pub,     only: die
    use grid,           only: x, nx, ny, nz, dz, dx, nb, nxb, nyb, nzb
    use mpisetup,       only: bnd_xl, bnd_xr, bnd_yl, bnd_yr, bnd_zl, bnd_zr

    implicit none
    real, dimension(:,:,:), intent(in)  :: dens
#ifdef SHEAR
    real, dimension(:,:,:), allocatable :: ala
#endif /* SHEAR */
    integer       :: i

    if (any([nx, ny, nz] <= 1)) call die("[poissonsolver:poisson_solve] Only 3D setups are supported") !BEWARE 2D and 1D probably need small fixes

    if ( bnd_xl .eq. 'per' .and. bnd_xr .eq. 'per' .and. &
        bnd_yl .eq. 'per' .and. bnd_yr .eq. 'per' .and. &
        bnd_zl .ne. 'per' .and. bnd_zr .ne. 'per'        ) then ! Periodic in X and Y, nonperiodic in Z

         call poisson_xyp(dens(nb+1:nb+nxb,nb+1:nb+nyb,:), &
                           sgp(nb+1:nb+nxb,nb+1:nb+nyb,:),dz)

        call die("[poissonsolver:poisson_solve] poisson_xyp called")

    elseif ( bnd_xl .eq. 'per' .and. bnd_xr .eq. 'per' .and. &
            bnd_yl .eq. 'per' .and. bnd_yr .eq. 'per' .and. &
            bnd_zl .eq. 'per' .and. bnd_zr .eq. 'per'        ) then ! Fully 3D periodic
        call poisson_xyzp(dens(nb+1:nb+nxb,nb+1:nb+nyb,nb+1:nb+nzb), sgp(nb+1:nb+nxb,nb+1:nb+nyb,nb+1:nb+nzb))

#ifdef SHEAR
    elseif ( bnd_xl .eq. 'she' .and. bnd_xr .eq. 'she' .and. &
             bnd_yl .eq. 'per' .and. bnd_yr .eq. 'per' ) then ! 2D shearing box

         if (dimensions=='3d') then
           if (.not.allocated(ala)) allocate(ala(nx-2*nb,ny-2*nb,nz-2*nb))
           ala = dens(nb+1:nb+nxb,nb+1:nb+nyb,nb+1:nb+nzb)
           call poisson_xyzp(ala(:,:,:), sgp(nb+1:nb+nxb,nb+1:nb+nyb,nb+1:nb+nzb))

           sgp(:,:,1:nb)              = sgp(:,:,nzb+1:nzb+nb)
           sgp(:,:,nzb+nb+1:nzb+2*nb) = sgp(:,:,nb+1:2*nb)

         else
            call poisson_xy2d(dens(nb+1:nb+nxb,nb+1:nb+nyb,1), &
                            sgp(nb+1:nb+nxb,      nb+1:nb+nyb,1), &
                            sgp(1:nb,             nb+1:nb+nyb,1), &
                            sgp(nxb+nb+1:nxb+2*nb,nb+1:nb+nyb,1), &
                            dx)

         endif
         sgp(:,1:nb,:)              = sgp(:,nyb+1:nyb+nb,:)
         sgp(:,nyb+nb+1:nyb+2*nb,:) = sgp(:,nb+1:2*nb,:)

         if (allocated(ala)) deallocate(ala)
#endif /* SHEAR */
    else
       call die("[poissonsolver:poisson_solve not implemented for current boundary conditions")
    endif

  end subroutine poisson_solve

!!=====================================================================
!!
!! SUBROUTINE POISSON_XY2D: solves Poisson equation for periodic
!! bnd conditions in X and Y in 2D
!! ASSUMPTION: 1) dx = dy
!!
#ifdef SHEAR
  subroutine poisson_xy2d(den, pot, lpot, rpot, dx)
    use arrays,    only: x
    use constants, only: newtong, dpi
    use grid,      only: xmin, xmax, ymin, ymax, nxb, nyb, nb
    use shear,     only: dely
    implicit none

    real, dimension(:,:), intent(in)              :: den
    real, dimension(:,:), intent(out)             :: pot
    real, dimension(:,:), intent(out)             :: lpot, rpot
    real, intent(in)                              :: dx

    complex*16, dimension(size(den,1)/2 + 1)      :: ctmp
    complex*16, dimension(size(den,1))            :: ctm2
    complex*16, dimension(size(den,1))            :: ctm3
    complex*16, dimension(size(den,1))            :: ctm4

    complex*16, dimension(size(den,1),size(den,1)/2 +1 ) :: fft,fft2

    real(kind=8)   , dimension(size(den,1))       :: rtmp
    real(kind=8)   , dimension(size(den,1)/2 + 1) :: ky
    real(kind=8)   , dimension(size(den,1))       :: kx,xx

    real(kind=8)    :: lambda, factor, St

    integer(kind=8) :: pf1, pf2, ppb1, pb2

    complex(kind=8), parameter :: j = (0, 1)

    integer*8      :: n,np, p, q, i

    integer*8     , parameter :: FFTW_ESTIMATE = 64
    integer*8     , parameter :: FFTW_FORWARD  = 1
    integer*8     , parameter :: FFTW_BACKWARD = -1
!
!----------------------------------------------------------------------
!

    St = dely / (xmax - xmin)
    St = -St * nyb /(ymax - ymin)
! determine input array dimension
!
    n = size(den, 1)
    xx(:) = x(nb+1:nxb+nb)
    if (n /= size(den,2)) stop 'nx /= ny in poisson_xy'

! compute dimensions for complex arrays
!
    np = n / 2 + 1

! compute wave vectors
!
    ky(1) = 0.0
    do p = 2, np
      ky(p) = dpi * (p-1) / n !      kx(p) = dpi * (p - 1) / nx
    enddo

    kx(1) = 0.0
    do q = 2, np
      kx(q) = dpi * (q-1) / n !      ky(q) = dpi * (q - 1) / ny
      kx(n+2-q) = kx(q)
    enddo

! create plan for the forward 1D REAL => COMPLEX  FFT
!
    call dfftw_plan_dft_r2c_1d(pf1,  size(rtmp), rtmp, ctmp, FFTW_ESTIMATE)

! perform forward FFT for each k
!
    do i = 1, n
      rtmp(:)   = den(i,:)
      call dfftw_execute(pf1)
      fft(i,:) = ctmp(:) * exp(cmplx(0.0, St*ky(:)*xx(i)))
    enddo
    call dfftw_destroy_plan(pf1)

    call dfftw_plan_dft_1d(pf2, size(ctm2), ctm2, ctm3, FFTW_FORWARD, FFTW_ESTIMATE)

!! Shift and perform final FFT
!
    do i = 1, np
      ctm2(:) = fft(:,i)
      call dfftw_execute(pf2)
      fft2(:,i) = ctm3(:)
    enddo

    call dfftw_destroy_plan(pf2)

!! compute eigenvalues for each p and q and solve
!!
    fft2(:,:)  = dpi*newtong*dx*dx * fft2(:,:)
    do q = 1, np
      do p = 1, n
        factor = cos(kx(p)) + cos(ky(q)) - 2.0
        if (factor == 0) then
           fft2(p,q) = cmplx(0.0,0.0)
        else
           fft2(p,q) = fft2(p,q) / factor
        endif
      enddo
    enddo

!! create plan for the backward FFT
!!
    call dfftw_plan_dft_1d(ppb1, n, ctm3, ctm4, -1, FFTW_ESTIMATE)

    do i = 1,np
      ctm3(:) = fft2(:,i)
      call dfftw_execute(ppb1)
!      fft(:,i) = (1./n) * ctm4(:)
      fft(:,i) = ctm4(:)
    enddo
! destroy plan for the 1st inverse FFT
!
    call dfftw_destroy_plan(ppb1)

    call dfftw_plan_dft_c2r_1d(pb2, n, ctmp, rtmp, FFTW_ESTIMATE)

    do i = 1,n
      ctmp(:)   = fft(i,:)  * exp( cmplx(0.0 , -St*ky(:)*xx(i)) )
      call dfftw_execute(pb2)
      pot(i,:)  = rtmp(:) / n / n
    enddo
    do i = 1,nb
      ctmp(:)   = fft(n-nb+i,:) * exp( cmplx(0.0, -St*ky(:)*x(i)) )
      call dfftw_execute(pb2)
      lpot(i,:) = rtmp(:) / n / n

      ctmp(:)   = fft(i,:) * exp( cmplx(0.0, -St*ky(:)*x(nxb+nb+i)) )
      call dfftw_execute(pb2)
      rpot(i,:) = rtmp(:) / n / n
    enddo
! destroy plan for the 1st inverse FFT
!
    call dfftw_destroy_plan(pb2)

    return
  end subroutine poisson_xy2d
#endif /* SHEAR */

!!=====================================================================
!!
!! SUBROUTINE POISSON_XY: solves Poisson equation for periodic
!! bnd conditions in X and Y in 2D
!! ASSUMPTION: 1) dx = dy
!!
  subroutine poisson_xy(den, pot, dx)

    use constants, only: newtong, dpi
    implicit none

    real, dimension(:,:), intent(in)  :: den
    real, dimension(:,:), intent(out) :: pot
    real, intent(in)                  :: dx

    complex(kind=8), dimension(:,:)  , allocatable :: ctmp

    real(kind=8)   , dimension(:,:)  , allocatable :: rtmp
    real(kind=8)   , dimension(:)    , allocatable :: kx, ky

    real(kind=8)    :: factor

    integer(kind=8) :: planf, planb

    integer         :: nx, ny, np, p, q

    integer     , parameter :: FFTW_ESTIMATE = 64
!
!----------------------------------------------------------------------
!

! determine input array dimension
!
  nx = size(den, 1)
  ny = size(den, 2)
  if (nx /= ny) stop 'nx /= ny in poisson_xy'

! compute dimensions for complex arrays
!
    np = nx / 2 + 1

! allocate arrays
!
    allocate(ctmp(np,ny))
    allocate(rtmp(nx,ny))
    allocate(kx(np))
    allocate(ky(ny))

! create plan for the forward FFT
!
    call dfftw_plan_dft_r2c_2d(planf, nx, ny, rtmp, ctmp, FFTW_ESTIMATE)

! perform forward FFT for each k
!

!!!    factor = dx * dx
   factor = 1.0

    rtmp(:,:)   = factor * den(:,:)
    call dfftw_execute(planf)
!    fft(:,:) = ctmp(:,:) !/ real(nx * ny)

! destroy plan for the forward FFT
!
    call dfftw_destroy_plan(planf)

! compute wave vectors
!
    do p = 1, np
      kx(p) = dpi * (p - 1) / nx
      kx(p) = dpi * (p) / nx

    enddo

    ky(1) = 0.0
    ky(1) = dpi / ny
    do q = 2, ny / 2 + 1
      ky(q) = dpi * (q - 1) / ny
      ky(q) = dpi * (q) / ny
      ky(ny+2-q) = ky(q)
    enddo

! compute eigenvalues for each p and q and solve linear system
!

    ctmp(:,:)  = dpi*newtong*dx*dx * ctmp(:,:)

    do q = 1, ny
      do p = 1, np
        factor = cos(kx(p)) + cos(ky(q)) - 2.0 !+ 1.e-10
        ctmp(p,q) = ctmp(p,q) / factor
      enddo
    enddo

! create plan for the backward FFT
!
    call dfftw_plan_dft_c2r_2d(planb, nx, ny, ctmp, rtmp, FFTW_ESTIMATE)

! perform inverse FFT for each k
!
    call dfftw_execute(planb)
    pot(:,:) = rtmp(:,:)/(nx*ny)

! destroy plan for the backward FFT
!
    call dfftw_destroy_plan(planb)

! deallocate rhofft
!
    if (allocated(kx))     deallocate(kx)
    if (allocated(ky))     deallocate(ky)
    if (allocated(rtmp))   deallocate(rtmp)
    if (allocated(ctmp))   deallocate(ctmp)

  end subroutine poisson_xy

!!=====================================================================
!!
!! SUBROUTINE POISSON_XYP: solves Poisson equation for periodic
!! bnd conditions in X and Y, and arbitrary bnd in Z
!!
  subroutine poisson_xyp(den, pot, dz)

    use constants, only: fpiG, dpi
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

    !BEWARE: This routine will work incorrectly if dx /= dz or dy /= dz (poisson_xyzp was corrected in r2124)
    if (present(dz) .eqv. .true.) factor = dz * dz

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
    if (allocated(ipiv))   deallocate(ipiv)
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

!!=====================================================================
!!
!! SUBROUTINE POISSON_XYZP: solves Poisson equation for periodic
!! bnd conditions in X, Y and Z
!!
  subroutine poisson_xyzp(den, pot)

    use constants, only: fpiG, dpi
    use grid,      only: dx, dy, dz

    implicit none

    real, dimension(:,:,:), intent(in)  :: den
    real, dimension(:,:,:), intent(out) :: pot

    complex, dimension(:,:,:), allocatable :: fft
    complex, dimension(:,:,:), allocatable :: ctmp
    real   , dimension(:,:,:)  , allocatable :: rtmp
    real   , dimension(:)    , allocatable    :: kx, ky, kz
    real    :: norm
    integer(kind=selected_int_kind(16)) :: planf, plani
    integer         :: nx, ny, nz, np, i, j, k
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
    allocate(ctmp(np,ny,nz))
    allocate(rtmp(nx,ny,nz))
    allocate(kx(nx))
    allocate(ky(ny))
    allocate(kz(nz))

    norm = 1.0 / real( nx * ny * nz )

! BEWARE: the plans can probably be reused and it might be more efficient to create them with FFTW_MEASURE
! create plan for the forward FFT
    call dfftw_plan_dft_r2c_3d(planf, nx, ny, nz, den, ctmp, FFTW_ESTIMATE)
! perform forward FFT 3D
!    call dfftw_execute(planf) ! Some fortran compilers make segfaulting code with simple dfftw_execute(plan)
    call dfftw_execute_dft_r2c(planf, den, ctmp)
! destroy plan for the forward FFT
    call dfftw_destroy_plan(planf)

! compute eigenvalues for each p, q and r and solve linear system
! ToDo: this can be done only once if we do not change arrays

    kx(:) = (cos(dpi/nx*(/(j,j=0,np-1)/))-1.)/dx**2
    ky(:) = (cos(dpi/ny*(/(j,j=0,ny-1)/))-1.)/dy**2
    kz(:) = (cos(dpi/nz*(/(j,j=0,nz-1)/))-1.)/dz**2

! Correction for 4-th order (5-point) Laplacian - seems not to be important, at least in Jeans test.
! For integral approximation of the 4-th order Laplacian replace (7.-cos(x))/6. by (13.-cos(x))/12.
!    kx(:) = kx(:) * (7.-cos(dpi/nx*(/(j,j=0,np-1)/)))/6.
!    ky(:) = ky(:) * (7.-cos(dpi/ny*(/(j,j=0,ny-1)/)))/6.
!    kz(:) = kz(:) * (7.-cos(dpi/nz*(/(j,j=0,nz-1)/)))/6.
    forall (i=1:np,j=1:ny,k=1:nz, (kx(i)+ky(j)+kz(k) - 3.0) /= 0.0)
       fft(i,j,k) = 0.5 * ctmp(i,j,k) / (kx(i)+ky(j)+kz(k) - 3.0)
    endforall

! create plan for the inverse FFT3D
    call dfftw_plan_dft_c2r_3d(plani, nx, ny, nz, fft, rtmp, FFTW_ESTIMATE)
! perform inverse FFT 3D
    call dfftw_execute(plani)
! destroy plan for the inverse FFT3D
    call dfftw_destroy_plan(plani)

    pot(:,:,:) = fpiG * rtmp(:,:,:) * norm

! deallocate temporary arrays
!
    if (allocated(kx))     deallocate(kx)
    if (allocated(ky))     deallocate(ky)
    if (allocated(kz))     deallocate(kz)
    if (allocated(rtmp))   deallocate(rtmp)
    if (allocated(ctmp))   deallocate(ctmp)
    if (allocated(fft))    deallocate(fft)

!------------------------------------------
    end subroutine poisson_xyzp

#else /* !POISSON_FFT */
#warning This should not happen. Probably the poissonsolver.F90 file is included in object directory by mistake.
#endif /* !POISSON_FFT */

end module poissonsolver
