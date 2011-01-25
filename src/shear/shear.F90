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
   public :: csvk, delj, dely, eps, eta_gas, global_gradP, init_shear, omega, qshear, yshift
#ifdef FFTW
   public  :: unshear_fft
#endif FFTW
   real    :: ts, dely, eps, omega, qshear, dts, ddly, eta_gas, csvk
   integer :: delj

   real, dimension(:), allocatable :: global_gradP

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

    use dataio_pub,     only: par_file, ierrh, namelist_errh, compare_namelist, cmdl_nml  ! QA_WARN required for diff_nml
    use dataio_pub,     only: printinfo, die, code_progress, PIERNIK_INIT_BASE
    use mpisetup,       only: ierr, master, slave, rbuff, buffer_dim, comm
    use mpi,            only: MPI_DOUBLE_PRECISION
    use fluidindex,     only: flind

    implicit none
    integer       :: i

    namelist /SHEARING/ omega, qshear, eta_gas, csvk

    if (code_progress < PIERNIK_INIT_BASE) call die("[shear:init_shear] fluids not initialized.")

#ifdef VERBOSE
    call printinfo("[shear:init_shear]: commencing...")
#endif /* VERBOSE */

    omega   = 0.0
    qshear  = 0.0

    if (master) then
       diff_nml(SHEARING)

       rbuff(1) = omega
       rbuff(2) = qshear
       rbuff(3) = eta_gas
       rbuff(4) = csvk

    endif

    call MPI_Bcast(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

    if (slave) then
       omega   = rbuff(1)
       qshear  = rbuff(2)
       eta_gas = rbuff(3)
       csvk    = rbuff(4)
    endif

#ifdef VERBOSE
    call printinfo("[shear:init_shear]: finished. \o/")
#endif /* VERBOSE */

    allocate(global_gradP(flind%fluids))
    do i = 1, flind%fluids
       global_gradP(i) = 2.0*omega * eta_gas * flind%all_fluids(i)%cs / csvk
    enddo

  end subroutine init_shear

  subroutine yshift(ts,dts)

    use arrays,   only: u
    use grid,     only: cg
    use mpisetup, only: dom

    implicit none

    real, intent(in) :: ts, dts
#ifdef FFTW
    integer :: i
#endif /* FFTW */

    ddly  = dts * qshear*omega*dom%Lx
    dely  = ts  * qshear*omega*dom%Lx
    delj  = mod(int(dely/cg%dy), cg%nyb)
    eps   = mod(dely, cg%dy)/cg%dy
#ifdef FFTW
    do i=LBOUND(u,1),UBOUND(u,1)
       u(i,:, cg%js:cg%je,:) = unshear_fft( u(i,:, cg%js:cg%je,:), cg%x(:),ddly)
    enddo
    u(:,:,1:cg%nb,:)       = u(:,:, cg%ny-2*cg%js:cg%je,:)
    u(:,:, cg%je+1:cg%ny,:) = u(:,:, cg%js:cg%jsb,:)
#endif /* FFTW */
  end subroutine yshift

#ifdef FFTW
  function unshear_fft(qty,x,ddy,inv)

    use constants, only: dpi
    use grid,      only: cg
    use mpisetup,  only: dom

    implicit none

!    include "fftw3.f" ! this may give tons of warnings on unused parameters

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
    integer, parameter :: FFTW_ESTIMATE=64

    St = - ddy / cg%dy / dom%Lx
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

    use grid,     only: cg
    use mpisetup, only: dom

    implicit none

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

    my = 3*cg%nyb+2*cg%nb

    fx = dely / dom%Lx
    sg = -1

    if (.not.allocated(temp)) allocate(temp(my,nz))

    unshear = 0.0

    if (present(inv)) then
       fx = - fx
    endif
    do i = 1, cg%nx
      dl  = fx * x(i)
      ndl = mod(int(dl/cg%dy), cg%nyb)
      ddl = mod(dl, cg%dy)/cg%dy

      temp(         1:  cg%je,:)   = qty(i,   1:cg%je ,:)
      temp(  cg%je+1:2*cg%nyb+cg%nb,:)   = qty(i, cg%js:cg%je,:)
      temp(2*cg%nyb+cg%nb+1:3*cg%nyb+2*cg%nb,:) = qty(i, cg%js:ny    ,:)

      temp = cshift(temp,dim=1,shift=ndl)

!      temp(1:nb,:) = temp(nyb+1:nyb+nb,:)          ! not needed
!      temp(nb+nyb+1:nyb+2*nb,:) = temp(nb+1:2*nb,:)

      temp(:,:) = (1.0+ddl)*(1.0-ddl) * temp(:,:) &
            - 0.5*(ddl)*(1.0-ddl) * cshift(temp(:,:),shift= sg,dim=1) &
            + 0.5*(ddl)*(1.0+ddl) * cshift(temp(:,:),shift=-sg,dim=1)

      unshear(i, cg%js:cg%je,:) = temp(cg%je+1:cg%nb+2*cg%nyb,:)

      unshear(i,1:cg%nb,:)          = unshear(i, cg%jeb:cg%je,:)
      unshear(i, cg%je+1:ny,:)   = unshear(i, cg%js:cg%jsb,:)

!      unshear(i,:,:) = max(unshear(i,:,:), smalld)
    enddo
    if (allocated(temp)) deallocate(temp)
    return
  end function unshear
end module shear
