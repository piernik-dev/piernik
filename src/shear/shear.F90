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
!! \brief Module of shearing box routines
!!
!! In this module following namelist of parameters is specified:
!! \copydetails shear::init_shear
!<
module shear
! pulled by SHEAR
   implicit none

   private
   public :: csvk, delj, dely, eps, eta_gas, global_gradP, init_shear, omega, qshear, shear_3sweeps, shear_acc, yshift
#ifdef FFTW
   public  :: unshear_fft
#endif /* FFTW */

   include "fftw3.f"

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
!! <tr><td>omega  </td><td>0.0  </td><td>real value</td><td>\copydoc shear::omega  </td></tr>
!! <tr><td>qshear </td><td>0.0  </td><td>real value</td><td>\copydoc shear::qshear </td></tr>
!! <tr><td>eta_gas</td><td>0.0  </td><td>real value</td><td>\copydoc shear::eta_gas</td></tr>
!! <tr><td>csvk   </td><td>1.0  </td><td>real value</td><td>\copydoc shear::csvk   </td></tr>
!! </table>
!! The list is active while \b "SHEAR" is defined.
!! \n \n
!<
   subroutine init_shear

      use constants,      only: PIERNIK_INIT_GRID
      use dataio_pub,     only: nh  ! QA_WARN required for diff_nml
      use dataio_pub,     only: printinfo, die, code_progress
      use fluidindex,     only: flind
      use mpisetup,       only: master, slave, rbuff, piernik_MPI_Bcast

      implicit none
      integer       :: i

      namelist /SHEARING/ omega, qshear, eta_gas, csvk

      if (code_progress < PIERNIK_INIT_GRID) call die("[shear:init_shear] fluids not initialized.")

#ifdef VERBOSE
      call printinfo("[shear:init_shear]: commencing...")
#endif /* VERBOSE */

      omega   = 0.0
      qshear  = 0.0
      eta_gas = 0.0
      csvk    = 1.0

      if (master) then
         diff_nml(SHEARING)

         rbuff(1) = omega
         rbuff(2) = qshear
         rbuff(3) = eta_gas
         rbuff(4) = csvk

      endif

      call piernik_MPI_Bcast(rbuff)

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
         global_gradP(i) = 2.0*omega * eta_gas * flind%all_fluids(i)%fl%cs / csvk
      enddo

   end subroutine init_shear
!--------------------------------------------------------------------------------------------------

   subroutine shear_3sweeps

      use cg_leaves,           only: leaves
      use constants,           only: xdim, ydim
      use dataio_pub,          only: die
      use domain,              only: dom, is_multicg
      use fluidboundaries,     only: bnd_u
      use global,              only: t, dt
      use grid_cont,           only: grid_container

      implicit none
      type(grid_container), pointer :: cg

      cg => leaves%first%cg
      if (is_multicg) call die("[shear:shear_3sweeps] multiple grid pieces per processor not implemented yet") !nontrivial SHEAR

      if (dom%has_dir(ydim)) call yshift(t, dt)
      if (dom%has_dir(xdim)) call bnd_u(xdim, cg)
      if (dom%has_dir(ydim)) call bnd_u(ydim, cg)

   end subroutine shear_3sweeps

!--------------------------------------------------------------------------------------------------
   function shear_acc(sweep,u) result(rotacc)

      use constants,  only: xdim, ydim
      use fluidindex, only: flind, iarr_all_dn, iarr_all_my

      implicit none

      real, dimension(:,:), intent(in)        :: u
      integer(kind=4),      intent(in)        :: sweep

      real, dimension(2)                      :: df                 !< \deprecated  additional acceleration term used in streaming problem
      real, dimension(flind%fluids,size(u,2)) :: vy0
      real, dimension(flind%fluids,size(u,2)) :: rotacc
      integer :: ind

      df = 0.0
#ifdef FLUID_INTERACTIONS
      df = global_gradP       !! \deprecated BEWARE: only for backward compatibility with old streaming problem
#endif /* !FLUID_INTERACTIONS */
      where (u(iarr_all_dn,:) > 0.0)
         vy0(:,:)  = u(iarr_all_my,:)/u(iarr_all_dn,:)
      elsewhere
         vy0(:,:)  = 0.0
      endwhere
      do ind = 1, flind%fluids
!        if (sweep == xdim) then
!           rotacc(ind,:) =  2.0*omega*(vy0(ind,:) + qshear*omega*cg%x(:))
!        else if (sweep == ydim)  then
!           rotacc(ind,:) = - 2.0*omega*vy0(ind,:)          ! with global shear
!        else
!           rotacc(ind,:) = 0.0
!        endif
         if (sweep == xdim) then
            rotacc(ind,:) =  2.0*omega*vy0(ind,:) + df(ind)  ! global_gradient
         else if (sweep == ydim)  then
            rotacc(ind,:) = (qshear - 2.0)*omega*vy0(ind,:)  ! with respect to global shear (2.5D)
         else
            rotacc(ind,:) = 0.0
         endif
      enddo

   end function shear_acc

   subroutine yshift(ts,dts)

      use cg_leaves,   only: leaves
      use constants,   only: xdim
#ifdef FFTW
      use constants,   only: ydim, LO, HI
#endif /* FFTW */
      use dataio_pub,  only: die
      use domain,      only: dom, is_multicg
      use grid_cont,   only: grid_container

      implicit none

      real, intent(in) :: ts, dts

      type(grid_container), pointer :: cg
#ifdef FFTW
      integer :: i
#endif /* FFTW */

      cg => leaves%first%cg
      if (is_multicg) call die("[shear:yshift] multiple grid pieces per procesor not implemented yet") !nontrivial

      ddly  = dts * qshear*omega*dom%L_(xdim)
      dely  = ts  * qshear*omega*dom%L_(xdim)
      delj  = mod(int(dely*cg%idy), int(cg%nyb))
      eps   = mod(dely, cg%dy)*cg%idy
#ifdef FFTW
      do i=lbound(cg%u,1),ubound(cg%u,1)
         cg%u(i,:, cg%js:cg%je,:) = unshear_fft( cg%u(i,:, cg%js:cg%je,:), cg%x(:),ddly)
      enddo
      cg%u(:,:,cg%lhn(ydim, LO):cg%lh1(ydim, LO),:) = cg%u(:,:, cg%jeb:cg%je,:)
      cg%u(:,:, cg%lh1(ydim, HI):cg%lhn(ydim, HI),:) = cg%u(:,:, cg%js:cg%jsb,:)
#endif /* FFTW */

   end subroutine yshift
!--------------------------------------------------------------------------------------------------
#ifdef FFTW
   function unshear_fft(qty,x,ddy,inv)

      use cg_leaves,   only: leaves
      use constants,   only: dpi, xdim
      use dataio_pub,  only: die
      use domain,      only: dom, is_multicg
      use grid_cont,   only: grid_container

      implicit none

      !integer, parameter :: FFTW_ESTIMATE=64

      real,                   intent(in) :: ddy
      real, dimension(:),     intent(in) :: x
      real, dimension(:,:,:), intent(in) :: qty
      logical,                optional   :: inv
      real, dimension(size(qty,1),size(qty,2),size(qty,3)) :: unshear_fft
      integer :: nx,ny,nz,p,np,q
      real    :: St

      integer(kind=8)                            :: planf, planb

      complex(kind=8), dimension(:), allocatable :: ctmp
      real(kind=8),    dimension(:), allocatable :: rtmp
      real(kind=8),    dimension(:), allocatable :: ky
      type(grid_container), pointer              :: cg

      cg => leaves%first%cg
      if (is_multicg) call die("[shear:unshear_fft] multiple grid pieces per procesor not implemented yet") !nontrivial

      St = - ddy * cg%idy / dom%L_(xdim)
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
!--------------------------------------------------------------------------------------------------
   function unshear(qty,x,inv)

      use cg_leaves,   only: leaves
      use constants,   only: xdim, half
      use dataio_pub,  only: die
      use domain,      only: dom, is_multicg
      use grid_cont,   only: grid_container

      implicit none

      real, dimension(:,:,:), intent(in) :: qty
      real, dimension(:),     intent(in) :: x
      logical,                optional   :: inv
      real, dimension(size(qty,1),size(qty,2),size(qty,3)) :: unshear
      real, dimension(:,:), allocatable  :: temp
      integer :: i, sg, my, nx, ny, nz, ndl
      real    :: fx, dl, ddl
      type(grid_container), pointer :: cg

      nx = size(qty,1)
      ny = size(qty,2)
      nz = size(qty,3)

      cg => leaves%first%cg
      if (is_multicg) call die("[shear:unshear] multiple grid pieces per procesor not implemented yet") !nontrivial

      my = 3*cg%nyb+2*dom%nb

      fx = dely / dom%L_(xdim)
      sg = -1

      if (.not.allocated(temp)) allocate(temp(my,nz))

      unshear = 0.0

      if (present(inv)) fx = - fx

      do i = 1, cg%n_(xdim)
         dl  = fx * x(i)
         ndl = mod(int(dl*cg%idy), int(cg%nyb))
         ddl = mod(dl, cg%dy)*cg%idy

         temp(         1:  cg%je,:)   = qty(i,   1:cg%je ,:)
         temp(  cg%je+1:2*cg%nyb+dom%nb,:)   = qty(i, cg%js:cg%je,:)
         temp(2*cg%nyb+dom%nb+1:3*cg%nyb+2*dom%nb,:) = qty(i, cg%js:ny    ,:)

         temp = cshift(temp,dim=1,shift=ndl)

         !      temp(1:nb,:) = temp(nyb+1:nyb+nb,:)          ! not needed
         !      temp(nb+nyb+1:nyb+2*nb,:) = temp(nb+1:2*nb,:)

         temp(:,:) = (1.0+ddl)*(1.0-ddl) * temp(:,:) &
              - half*(ddl)*(1.0-ddl) * cshift(temp(:,:),shift= sg,dim=1) &
              + half*(ddl)*(1.0+ddl) * cshift(temp(:,:),shift=-sg,dim=1)

         unshear(i, cg%js:cg%je,:) = temp(cg%je+1:dom%nb+2*cg%nyb,:)

         unshear(i,1:dom%nb,:)          = unshear(i, cg%jeb:cg%je,:)
         unshear(i, cg%je+1:ny,:)   = unshear(i, cg%js:cg%jsb,:)

         !      unshear(i,:,:) = max(unshear(i,:,:), smalld)
      enddo

      if (allocated(temp)) deallocate(temp)

   end function unshear
!--------------------------------------------------------------------------------------------------
#ifdef SHEAR_BND
   subroutine bnd_shear_u(dir, cg)

      use dataio_pub, only: die
      use grid_cont,  only: grid_container

      implicit none

      integer(kind=4),               intent(in)    :: dir
      type(grid_container), pointer, intent(inout) :: cg

      call die("[shear:bnd_shear_u] Unimplemented")

      if (.false.) cg%is_old = cg%is_old .or. (dir == -1) ! suppress compiler warnings

   end subroutine bnd_shear_u
#endif /* SHEAR_BND */
end module shear
