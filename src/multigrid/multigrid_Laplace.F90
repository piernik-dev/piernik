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

!> \brief Multigrid Laplace operator selector routines

module multigrid_Laplace
! pulled by MULTIGRID && GRAV

   implicit none

   private
   public :: residual, approximate_solution_relax, vT_A_v, ord_laplacian, ord_laplacian_outer

   integer(kind=4) :: ord_laplacian          !< Laplace operator order; allowed values are 2, -4 (default) and 4 (not fully implemented)
   integer(kind=4) :: ord_laplacian_outer    !< Laplace operator order for isolated boundaries (useful as long as -4 is not fully implemented)

contains

!> \brief Select appropriate order of laplacian, depending on which solve operation we're on

   function ordL()

      use multigridvars, only: grav_bnd, bnd_givenval

      implicit none

      integer(kind=4) :: ordL

      ordL = ord_laplacian
      if (grav_bnd == bnd_givenval) ordL = ord_laplacian_outer

   end function ordL

!>
!! \brief Select requested routine for calculation of the residuum for the Poisson equation.
!!
!! \details Different orders of the Laplace operator result in different quality of the approximation of potential
!! * Second order (O_I2) is the simplest operator (3-point in 1D, 5-point in 2D and 7-point in 3D). Its error is proportional to cg%dx**2
!! * Fourth order (O_I4) has width of 5 cells in each direction (5, 9 and 13-points in 1D, 2D and 3D, respectively). Its error is proportional to cg%dx**4.
!!   The value of its coefficients depends on interpretation of variables (e.g. point values, integral over cell)
!! * Fourth order Mehrstellen (-O_I4) is compact (takes 3, 9, and 27 solution values in 1D, 2D and 3D),
!!   but requires values of adjacent source cells as well (3, 5 and 7 points in 1D, 2D and 3D).
!!   Its error is proportional to cg%dx**4 and tends to be smaller than the error of simple fourth order operator by a factor of 3..4.
!!   \todo check how the coefficients depend on interpretation of cell values (center point values vs integral over cell)
!!
!! Note that each kind of Laplace operator requires its own approximate solver (relaxation scheme, Green function for FFT) for optimal convergence.
!! \warning Relaxation is not implemented for the fourth order operator.
!<

   subroutine residual(cg_llst, src, soln, def)

      use cg_list_bnd,         only: cg_list_bnd_T
      use constants,           only: O_I2, O_I4
      use dataio_pub,          only: die
      use multigrid_Laplace2,  only: residual2
      use multigrid_Laplace4,  only: residual4
      use multigrid_Laplace4M, only: residual_Mehrstellen

      implicit none

      class(cg_list_bnd_T), intent(inout) :: cg_llst !< pointer to a list of grids for which we approximate the solution
      integer(kind=4),      intent(in) :: src     !< index of source in cg%q(:)
      integer(kind=4),      intent(in) :: soln    !< index of solution in cg%q(:)
      integer(kind=4),      intent(in) :: def     !< index of defect in cg%q(:)

      if (any(def == [ src, soln ])) call die("[multigrid_Laplace:residual] Cannot put the result into one of the input fields.") ! Use %q_copy method in such case
      select case (ordL())
         case (O_I2)
            call residual2(cg_llst, src, soln, def)
         case (O_I4)
            call residual4(cg_llst, src, soln, def)
         case (-O_I4)
            call residual_Mehrstellen(cg_llst, src, soln, def)
         case default
            call die("[multigrid_Laplace:residual] The order of Laplacian must be equal to 2, 4 or -4")
      end select

   end subroutine residual

!>
!! \brief Relaxation selector routine.
!!
!! \details The relaxation routines are the most costly routine in a serial run. Try to find optimal values for nsmoo.
!! The relaxation routines also depends a lot on communication, which may limit scalability of the multigrid.
!<

   subroutine approximate_solution_relax(curl, src, soln, nsmoo)

      use cg_level_connected,  only: cg_level_connected_T
      use constants,           only: O_I2, O_I4
      use dataio_pub,          only: die
      use multigrid_Laplace2,  only: approximate_solution_rbgs2
      use multigrid_Laplace4,  only: approximate_solution_relax4
      use multigrid_Laplace4M, only: approximate_solution_relax4M

      implicit none

      type(cg_level_connected_T), pointer, intent(inout) :: curl  !< pointer to a level for which we approximate the solution
      integer(kind=4),                     intent(in)    :: src   !< index of source in cg%q(:)
      integer(kind=4),                     intent(in)    :: soln  !< index of solution in cg%q(:)
      integer(kind=4),                     intent(in)    :: nsmoo !< Number of smoothing operations to perform

      select case (ordL())
         case (O_I2)
            call approximate_solution_rbgs2  (curl, src, soln, nsmoo)
         case (O_I4)
            call approximate_solution_relax4 (curl, src, soln, nsmoo)
         case (-O_I4)
            call approximate_solution_relax4M(curl, src, soln, nsmoo)
         case default
            call die("[multigrid_Laplace:approximate_solution_relax] The order of Laplacian must be equal to 2, 4 or -4")
      end select

   end subroutine approximate_solution_relax

!> \brief Selector for v*Laplacian(v) routine

   real function vT_A_v(var)

      use cg_leaves,          only: leaves
      use cg_list_global,     only: all_cg
      use constants,          only: O_I2, GEO_XYZ, ndims, INVALID, dsetnamelen
      use dataio_pub,         only: warn
      use domain,             only: dom
      use mpisetup,           only: master
      use multigrid_Laplace2, only: vT_A_v_2
      use named_array_list,   only: qna

      implicit none

      integer(kind=4), intent(in) :: var   !< Variable on which we want to calculate the operation

      logical, save :: firstcall = .true.
      character(len=dsetnamelen), parameter :: cg_L_n = "cg_L"
      integer(kind=4), save :: cg_L = INVALID

      if (dom%geometry_type == GEO_XYZ .and. ordL() == O_I2 .and. dom%eff_dim == ndims) then
         vT_A_v = vT_A_v_2(var)
      else
         if (firstcall) then
            if (master) call warn("[multigrid_Laplace:vT_A_v] No direct support for v*Laplacian(v) operation. Using workaround (a bit more expensive).")
            call all_cg%reg_var(cg_L_n)
            cg_L = qna%ind(cg_L_n)
         endif
         firstcall = .false.
         call leaves%set_q_value(qna%wai, 0.)
         call residual(leaves, qna%wai, var, cg_L)
         vT_A_v = -leaves%scalar_product(var, cg_L)
      endif

   end function vT_A_v

end module multigrid_Laplace
