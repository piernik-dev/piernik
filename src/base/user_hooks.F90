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
!>
!! \brief Definitions of subroutine templates for user hooks
!<
module user_hooks

   implicit none

   private
   public :: problem_customize_solution, problem_grace_passed, finalize_problem, cleanup_problem, problem_refine_derefine, &
        &    custom_emf_bnd, at_user_dims_settings, at_user_area_settings, problem_post_restart, user_vars_arr_in_restart

   interface

      subroutine no_args
         implicit none
      end subroutine no_args

      subroutine tab_args(tab)
         implicit none
         real, dimension(:,:,:), intent(inout) :: tab
      end subroutine tab_args

      subroutine indx_args(ll,lr,ch,lo)
         implicit none
         integer(kind=4), dimension(:), intent(out) :: ll, lr
         integer,         dimension(:), intent(out) :: ch
         integer(kind=8), dimension(:), intent(out) :: lo
      end subroutine indx_args

      subroutine user_area(area)
         implicit none
         integer, dimension(:), intent(out) :: area
      end subroutine user_area

   end interface

   procedure(no_args),   pointer :: problem_customize_solution => NULL()
   procedure(no_args),   pointer :: problem_grace_passed       => NULL()
   procedure(no_args),   pointer :: user_vars_arr_in_restart   => NULL()
   procedure(no_args),   pointer :: problem_post_restart       => NULL()
   procedure(no_args),   pointer :: finalize_problem           => NULL()
   procedure(no_args),   pointer :: cleanup_problem            => NULL()
   procedure(no_args),   pointer :: problem_refine_derefine    => NULL()
   procedure(tab_args),  pointer :: custom_emf_bnd             => NULL()
   procedure(indx_args), pointer :: at_user_dims_settings      => NULL()
   procedure(user_area), pointer :: at_user_area_settings      => NULL()

end module user_hooks
