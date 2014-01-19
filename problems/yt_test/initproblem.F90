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

module initproblem

   use constants, only: dsetnamelen

   implicit none

   private
   public :: read_problem_par, problem_initial_conditions, problem_pointers

   ! namelist parameters
   real :: span  !< If > 0 then do 2 refinements; use span value to scale their distance
   real :: width !< Relative size of the refined regions
   namelist /PROBLEM_CONTROL/ span, width

   ! other private data
   character(len=dsetnamelen), parameter :: const_n = "const", lin_n = "linear", quad_n = "parabolic"

contains

!> Set up some user hooks

   subroutine problem_pointers

      use dataio_user, only: user_reg_var_restart
      use user_hooks,  only: problem_refine_derefine

      implicit none

      user_reg_var_restart    => register_user_var
      problem_refine_derefine => add_a_patch

   end subroutine problem_pointers

!> \brief Read runtime parameters

   subroutine read_problem_par

      use dataio_pub, only: nh      ! QA_WARN required for diff_nml
      use mpisetup,   only: rbuff, master, slave, piernik_MPI_Bcast

      implicit none

      ! namelist default parameter values
      span  = 0.
      width = 1.

      if (master) then

         if (.not.nh%initialized) call nh%init()
         open(newunit=nh%lun, file=nh%tmp1, status="unknown")
         write(nh%lun,nml=PROBLEM_CONTROL)
         close(nh%lun)
         open(newunit=nh%lun, file=nh%par_file)
         nh%errstr=""
         read(unit=nh%lun, nml=PROBLEM_CONTROL, iostat=nh%ierrh, iomsg=nh%errstr)
         close(nh%lun)
         call nh%namelist_errh(nh%ierrh, "PROBLEM_CONTROL")
         read(nh%cmdl_nml,nml=PROBLEM_CONTROL, iostat=nh%ierrh)
         call nh%namelist_errh(nh%ierrh, "PROBLEM_CONTROL", .true.)
         open(newunit=nh%lun, file=nh%tmp2, status="unknown")
         write(nh%lun,nml=PROBLEM_CONTROL)
         close(nh%lun)
         call nh%compare_namelist()

         rbuff(1) = span
         rbuff(2) = width

      endif

      call piernik_MPI_Bcast(rbuff)

      if (slave) then

         span = rbuff(1)
         width = rbuff(2)

      endif

      call register_user_var

   end subroutine read_problem_par

!> \brief Calculate the Mandelbrot iterations for all leaves

   subroutine problem_initial_conditions

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: small
      use dataio_pub,       only: warn
      use grid_cont,        only: grid_container
      use fluidindex,       only: iarr_all_dn
      use mpisetup,         only: master
      use named_array_list, only: qna!, wna

      implicit none

      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg
      integer :: i, j, k
      real, dimension(:,:,:), pointer :: f0, f1, f2

      ! Create the initial density arrays
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         cg%b(:, :, :, :) = 0.
         cg%u(:, :, :, :) = 0.

         f0 => cg%q(qna%ind(const_n))%arr
         f1 => cg%q(qna%ind(lin_n  ))%arr
         f2 => cg%q(qna%ind(quad_n ))%arr
         if (.not. associated(f0) .or. .not. associated(f1) .or. .not. associated(f2)) then
            if (master) call warn("[initproblem:analytic_solution] Cannot store the data")
            return
         endif

         cg%u(:, :, :, :) = small
         do k = cg%ks, cg%ke
            do j = cg%js, cg%je
               do i = cg%is, cg%ie
                  f0(i, j, k) = 7+1e-5*(i*j*k)/(cg%ie*cg%je*cg%ke)
                  f1(i, j, k) = cg%x(i)    + cg%y(j)    + cg%z(k)
                  f2(i, j, k) = cg%x(i)**2 + cg%y(j)**2 + cg%z(k)**2
                  cg%u(iarr_all_dn(1), i, j, k) = f0(i, j, k)
               enddo
            enddo
         enddo

         cgl => cgl%nxt
      enddo

! Why the statement below gives "Invalid writes" under valgrind?
!      call leaves%qw_copy(qna%ind(const_n), wna%fi, iarr_all_dn(1)) ! prevent spurious FP exceptions

   end subroutine problem_initial_conditions

!> \brief Add fields for iteration count and real and imaginary coordinate of the point at the end of iterations

   subroutine register_user_var

      use cg_list_global, only: all_cg
      use constants,      only: AT_NO_B

      implicit none

      call all_cg%reg_var(const_n, restart_mode = AT_NO_B)
      call all_cg%reg_var(lin_n,   restart_mode = AT_NO_B)
      call all_cg%reg_var(quad_n,  restart_mode = AT_NO_B)

   end subroutine register_user_var

!> \brief add a small patch on the top refinement level in a quite abusive way

   subroutine add_a_patch

      use cg_level_connected, only: cg_level_connected_T
      use cg_level_finest,    only: finest
      use cg_list,            only: cg_list_element
      use constants,          only: refinement_factor, xdim, ydim, zdim, LO
      use dataio_pub,         only: printinfo
      use refinement,         only: level_max

      implicit none

      type(cg_level_connected_T), pointer :: rlev
      type(cg_list_element), pointer :: cgl

      rlev => finest%level

      if (finest%level%level_id < level_max) then
         call printinfo("[initproblem:add_a_patch] adding a patch ...")
         call finest%add_finer
         cgl => rlev%first
         do while (associated(cgl))
            if (span <= 0.) then
               call finest%level%add_patch(int(cgl%cg%n_b(:)*width, kind=8), cgl%cg%my_se(:, LO)*refinement_factor + int(cgl%cg%n_b(:)*(1.-width/2.), kind=4))
            else
               call finest%level%add_patch(int(cgl%cg%n_b(:)*width, kind=8), &
                    int([ cgl%cg%my_se(xdim:ydim, LO)*refinement_factor + cgl%cg%n_b(xdim:ydim)*(1.-width/2.+span), rlev%n_d(zdim) - cgl%cg%n_b(zdim)*width/2.], kind=8))
               call finest%level%add_patch(int(cgl%cg%n_b(:)*width, kind=8), &
                    int([ cgl%cg%my_se(xdim:ydim, LO)*refinement_factor + cgl%cg%n_b(xdim:ydim)*(1.-width/2.-span), rlev%n_d(zdim) - cgl%cg%n_b(zdim)*width/2.], kind=8))
            endif
            cgl => cgl%nxt
         enddo
      endif

   end subroutine add_a_patch

end module initproblem
