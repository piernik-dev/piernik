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

!> Brief calculate the lovely shape of the Mandelbrot set. Refine on the set to follow the interesting details.

module initproblem

   use constants, only: dsetnamelen

   implicit none

   private
   public :: read_problem_par, problem_initial_conditions, problem_pointers

   ! namelist parameters
   integer(kind=4) :: order   !< Order of mandelbrot set
   integer(kind=4) :: maxiter !< Maximum number of iterations
   logical :: smooth_map      !< Try continuous colouring
   real :: ref_thr            !< threshold for refining a grid
   real :: deref_thr          !< threshold for derefining a grid

   namelist /PROBLEM_CONTROL/  order, maxiter, smooth_map, ref_thr, deref_thr

   ! other private data
   character(len=dsetnamelen), parameter :: mand_n = "mand", re_n = "real", imag_n = "imag"

contains

!> Set up some user hooks

   subroutine problem_pointers

      use dataio_user, only: user_reg_var_restart, user_vars_hdf5
      use user_hooks,  only: problem_refine_derefine

      implicit none

      user_reg_var_restart    => register_user_var
      user_vars_hdf5          => mand_vars
      problem_refine_derefine => mark_the_set

   end subroutine problem_pointers

!> \brief Read runtime parameters

   subroutine read_problem_par

      use dataio_pub, only: nh      ! QA_WARN required for diff_nml
      use dataio_pub, only: warn, die
      use domain,     only: dom
      use mpisetup,   only: lbuff, ibuff, rbuff, master, slave, piernik_MPI_Bcast

      implicit none

      ! namelist default parameter values
      order = 2
      maxiter = 100
      smooth_map = .true.
      ref_thr = 5.
      deref_thr = 2.

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

         ibuff(1) = order
         ibuff(2) = maxiter

         lbuff(1) = smooth_map

         rbuff(1) = ref_thr
         rbuff(2) = deref_thr

      endif

      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(lbuff)
      call piernik_MPI_Bcast(rbuff)

      if (slave) then

         order      = ibuff(1)
         maxiter    = ibuff(2)

         smooth_map = lbuff(1)

         ref_thr    = rbuff(1)
         deref_thr  = rbuff(2)

      endif

      if (any(dom%has_dir(:) .neqv. [ .true., .true., .false. ])) &
           call die("[initproblem:read_problem_par] Mandelbrot is supposed to by run only with XY plane and without Z-direction present")

      if (order /= 2) then
         call warn("[initproblem:read_problem_par] Only order == 2 is supported at the moment")
         order = 2
      endif

      if (ref_thr <= deref_thr) call die("[initproblem:read_problem_par] ref_thr <= deref_thr")

      call register_user_var

   end subroutine read_problem_par

!> \brief Calculate the Mandelbrot iterations for all leaves

   subroutine problem_initial_conditions

      use cg_list,          only: cg_list_element
      use cg_leaves,        only: leaves
      use dataio_pub,       only: warn
      use grid_cont,        only: grid_container
      use fluidindex,       only: iarr_all_dn
      use mpisetup,         only: master
      use named_array_list, only: qna, wna

      implicit none

      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg
      integer :: i, j, k, nit
      real, dimension(:,:,:), pointer :: mand, r__l, imag
      real, parameter :: bailout2 = 10.
      real :: zx, zy, zt, cx, cy, rnit

      ! Create the initial density arrays
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         cg%b(:, :, :, :) = 0.
         cg%u(:, :, :, :) = 0.

         mand => cg%q(qna%ind(mand_n))%arr
         r__l => cg%q(qna%ind(re_n  ))%arr
         imag => cg%q(qna%ind(imag_n))%arr
         if (.not. associated(mand) .or. .not. associated(r__l) .or. .not. associated(imag)) then
            if (master) call warn("[initproblem:analytic_solution] Cannot store the set")
            return
         endif

         do k = cg%ks, cg%ke
            do j = cg%js, cg%je
               cy = cg%y(j)
               do i = cg%is, cg%ie
                  cx = cg%x(i)

                  zx = cx
                  zy = cy
                  nit = 0
                  do while (zx*zx + zy*zy < bailout2 .and. nit < maxiter)
                     zt = zx*zx - zy*zy + cx
                     zy = 2*zx*zy + cy
                     zx = zt
                     nit = nit + 1
                  enddo

                  rnit = nit
                  if (smooth_map .and. zx*zx + zy*zy > bailout2) rnit = rnit + 1 - log(log(sqrt(zx*zx + zy*zy)))/log(2.)

                  mand(i, j, k) = log(rnit)
                  r__l(i, j, k) = zx
                  imag(i, j, k) = zy

               enddo
            enddo
         enddo

         cgl => cgl%nxt
      enddo

      call leaves%qw_copy(qna%ind(mand_n), wna%fi, iarr_all_dn(1)) ! prevent spurious FP exceptions

   end subroutine problem_initial_conditions

!> \brief Add fields for iteration count and real and imaginary coordinate of the point at the end of iterations

   subroutine register_user_var

      use cg_list_global, only: all_cg
      use constants,      only: AT_NO_B

      implicit none

      call all_cg%reg_var(mand_n, restart_mode = AT_NO_B)
      call all_cg%reg_var(re_n,   restart_mode = AT_NO_B)
      call all_cg%reg_var(imag_n, restart_mode = AT_NO_B)

   end subroutine register_user_var

!> \brief Make some derived variables available for output

   subroutine mand_vars(var, tab, ierrh, cg)

      use grid_cont,   only: grid_container
      use named_array_list, only: qna

      implicit none

      character(len=*), intent(in)                    :: var
      real(kind=4), dimension(:,:,:), intent(inout)   :: tab
      integer, intent(inout)                          :: ierrh
      type(grid_container), pointer, intent(in)       :: cg

      ierrh = 0
      select case (trim(var))
         case ("distance", "dist") ! Supply the alternative name to comply with the old 4-letter limit
            tab(:,:,:) = real(log(sqrt(cg%q(qna%ind(re_n))%span(cg%ijkse)**2 + cg%q(qna%ind(imag_n))%span(cg%ijkse)**2)), kind=4)
         case ("angle", "ang")
            tab(:,:,:) = real(atan2(cg%q(qna%ind(imag_n))%span(cg%ijkse), cg%q(qna%ind(re_n))%span(cg%ijkse)), kind=4)
         case default
            ierrh = -1
      end select

   end subroutine mand_vars

!> \brief Mark interesting blocks for refinement

   subroutine mark_the_set

      use cg_list,          only: cg_list_element
      use cg_leaves,        only: leaves
      use fluidindex,       only: iarr_all_dn
      use named_array_list, only: wna
      use refinement,       only: ref_flag

      implicit none

      type(cg_list_element), pointer :: cgl
      real :: nitd

      !call leaves%leaf_arr3d_boundaries(qna%ind(mand_n))

      cgl => leaves%first
      do while (associated(cgl))
         if (any(cgl%cg%leafmap)) then
            ! Cannot use mand_n as long as it stays uninitialized during second call to problem_refine_derefine in update_refinement
            ! wna%fi is vital and thus automagilaclly prolonged
            ! Possible fixes:
            ! * make mand_n vital
            ! * do not call problem_refine_derefine twice in update_refinement
!!$            nitd = exp(maxval(cgl%cg%q(qna%ind(mand_n))%span(cgl%cg%ijkse), mask=cgl%cg%leafmap)) - &
!!$                 & exp(minval(cgl%cg%q(qna%ind(mand_n))%span(cgl%cg%ijkse), mask=cgl%cg%leafmap))
             nitd = exp(maxval(cgl%cg%w(wna%fi)%arr(iarr_all_dn(1), cgl%cg%is:cgl%cg%ie, cgl%cg%js:cgl%cg%je, cgl%cg%ks:cgl%cg%ke), mask=cgl%cg%leafmap)) - &
                  & exp(minval(cgl%cg%w(wna%fi)%arr(iarr_all_dn(1), cgl%cg%is:cgl%cg%ie, cgl%cg%js:cgl%cg%je, cgl%cg%ks:cgl%cg%ke), mask=cgl%cg%leafmap))
            cgl%cg%refine_flags = ref_flag( nitd >= ref_thr, nitd < deref_thr )
         endif
         cgl => cgl%nxt
      enddo

   end subroutine mark_the_set

end module initproblem
