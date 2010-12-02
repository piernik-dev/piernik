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
#include "macros.h"
!>
!! \brief (KK)
!!
!! In this module following namelist of parameters is specified:
!! \copydetails shear::init_shear
!<
module shear
! pulled by SHEAR
   implicit none
   private
   public  :: init_shear, yshift, qshear, omega, eps, delj, dely
#ifdef FFTW
   public  :: unshear_fft
#endif FFTW
   real    :: ts, dely, eps, omega, qshear, dts, ddly
   integer :: delj

contains

!>
!! \brief Routine to set parameter values from namelist SHEARING
!!
!! \n \n
!! @b SHEARING
!! \n \n
!! <table border="+1">
!! <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!! <tr><td>omega </td><td>0.0  </td><td>real value</td><td>\copydoc shear::omega </td></tr>
!! <tr><td>qshear</td><td>0.0  </td><td>real value</td><td>\copydoc shear::qshear</td></tr>
!! </table>
!! \n \n
!<
  subroutine init_shear

    use dataio_pub,     only: par_file, ierrh, namelist_errh, compare_namelist  ! QA_WARN required for diff_nml
    use dataio_pub,     only: printinfo
    use mpisetup,       only: ierr, MPI_DOUBLE_PRECISION, proc, rbuff, buffer_dim, comm

    implicit none

    namelist /SHEARING/ omega, qshear

#ifdef VERBOSE
    call printinfo("[shear:init_shear]: commencing...")
#endif /* VERBOSE */

    omega   = 0.0
    qshear  = 0.0

    if (proc == 0) then
       diff_nml(SHEARING)

       rbuff(1) = omega
       rbuff(2) = qshear

    endif

    call MPI_Bcast(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

    if (proc /= 0) then
       omega   = rbuff(1)
       qshear  = rbuff(2)
    endif

#ifdef VERBOSE
    call printinfo("[shear:init_shear]: finished. \o/")
#endif /* VERBOSE */
  end subroutine init_shear

  subroutine yshift(ts,dts)
    use arrays, only: u
    use grid,   only: dy, Lx, nyb, x, nb, ny
    implicit none
    real, intent(in) :: ts, dts
#ifdef FFTW
    integer :: i
#endif /* FFTW */

    ddly  = dts * qshear*omega*Lx
    dely  = ts  * qshear*omega*Lx
    delj  = mod(int(dely/dy),nyb)
    eps   = mod(dely,dy)/dy
#ifdef FFTW
    do i=LBOUND(u,1),UBOUND(u,1)
       u(i,:,nb+1:ny-nb,:) = unshear_fft( u(i,:,nb+1:ny-nb,:), x(:),ddly)
    enddo
    u(:,:,1:nb,:)       = u(:,:,ny-2*nb+1:ny-nb,:)
    u(:,:,ny-nb+1:ny,:) = u(:,:,nb+1:2*nb,:)
#endif /* FFTW */
  end subroutine yshift


#ifdef FFTW
  function unshear_fft(qty,x,ddy,inv)
    use constants, only: dpi
    use grid,      only: dy, nb, Lx
    implicit none
    include "fftw3.f"
    real, intent(in) :: ddy
    logical, optional               :: inv
    real, dimension(:,:,:)          :: qty
    real, dimension(:), intent(in)  :: x
    real, dimension(size(qty,1),size(qty,2),size(qty,3)) :: unshear_fft
    integer :: nx,ny,nz,p,np,q
    real    :: St

    integer(kind=8) :: planf,planb

    complex*16      , dimension(:)   , allocatable :: ctmp
    real(kind=8)    , dimension(:)   , allocatable :: rtmp
    real(kind=8)    , dimension(:)     , allocatable :: ky

    ! constants from fftw3.f
!    integer, parameter :: FFTW_ESTIMATE=64

    St = - ddy / dy / Lx
    if (.not.present(inv)) St = -St

    nx = size(qty,1)
    ny = size(qty,2)
    nz = size(qty,3)

    np = ny / 2 + 1

    if (.not.allocated(ctmp)) allocate(ctmp(np))
    if (.not.allocated(rtmp)) allocate(rtmp(ny))
    if (.not.allocated(ky)  ) allocate(  ky(np))

    ky(1) = 0.0
    do p = 2, np
      ky(p) = dpi * (p-1) / ny
    enddo

    call dfftw_plan_dft_r2c_1d(planf, ny, rtmp, ctmp, FFTW_ESTIMATE)
    call dfftw_plan_dft_c2r_1d(planb, ny, ctmp, rtmp, FFTW_ESTIMATE)

    do q = 1, nx
       do p = 1, nz
         rtmp(:)  = qty(q,:,p)
         call dfftw_execute(planf)
         ctmp(:)  = ctmp(:)*exp(cmplx( 0.0, St*ky(:)*x(q) ) )
         call dfftw_execute(planb)
         unshear_fft(q,:,p)  = rtmp(:) / (ny)
       enddo
    enddo

    call dfftw_destroy_plan(planf)
    call dfftw_destroy_plan(planb)

    if (allocated(rtmp))   deallocate(rtmp)
    if (allocated(ctmp))   deallocate(ctmp)
    if (allocated(ky)  )   deallocate(ky  )
    return

  end function unshear_fft
#endif /* FFTW */

  function unshear(qty,x,inv)
    use grid,     only: nb, xmax, xmin, nyb, dy

    logical, optional               :: inv
    real, dimension(:,:,:)          :: qty
    real, dimension(:), intent(in)  :: x
    real, dimension(size(qty,1),size(qty,2),size(qty,3)) :: unshear
    real, dimension(:,:), allocatable:: temp
    integer :: i,sg,my,nx,ny,nz,ndl
    real    :: fx,dl,ddl

    nx = size(qty,1)
    ny = size(qty,2)
    nz = size(qty,3)

    my = 3*nyb+2*nb

    fx = dely / (xmax - xmin)
    sg = -1

    if (.not.allocated(temp)) allocate(temp(my,nz))

    unshear = 0.0

    if (present(inv)) then
       fx = - fx
    endif
    do i = 1, nx
      dl  = fx * x(i)
      ndl = mod(int(dl/dy),nyb)
      ddl = mod(dl,dy)/dy

      temp(         1:  nyb+nb,:)   = qty(i,   1:nyb+nb ,:)
      temp(  nyb+nb+1:2*nyb+nb,:)   = qty(i,nb+1:nyb+nb,:)
      temp(2*nyb+nb+1:3*nyb+2*nb,:) = qty(i,nb+1:ny    ,:)

      temp = cshift(temp,dim=1,shift=ndl)

!      temp(1:nb,:) = temp(nyb+1:nyb+nb,:)          ! not needed
!      temp(nb+nyb+1:nyb+2*nb,:) = temp(nb+1:2*nb,:)

      temp(:,:) = (1.0+ddl)*(1.0-ddl) * temp(:,:) &
            - 0.5*(ddl)*(1.0-ddl) * cshift(temp(:,:),shift= sg,dim=1) &
            + 0.5*(ddl)*(1.0+ddl) * cshift(temp(:,:),shift=-sg,dim=1)

      unshear(i,nb+1:nb+nyb,:) = temp(nb+nyb+1:nb+2*nyb,:)

      unshear(i,1:nb,:)          = unshear(i,nyb+1:nyb+nb,:)
      unshear(i,nyb+nb+1:ny,:)   = unshear(i,nb+1 :2*nb,:)

!      unshear(i,:,:) = max(unshear(i,:,:), smalld)
    enddo
    if (allocated(temp)) deallocate(temp)
    return
  end function unshear
end module shear
