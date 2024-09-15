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
!! \brief calculate the lovely shape of the Mandelbrot set.  Refine on the
!! set to follow the interesting details.
!!
!! The Mandelbrot problem is intended for stress testing of the AMR
!! subsystem. There are more efficient fractal generators around but
!! one may consider some fun ideas:
!! * Use some multiprecision or implement own fixed point, optimized for
!!   these calculations.
!! * Use speedup tricks like these in fast deep zoom programs.
!! * Detect interior of minibrotsfor further speedups (already sort of works
!!   as it dos not get refined too much.
!<

module initproblem

   use constants, only: dsetnamelen

   implicit none

   private
   public :: read_problem_par, problem_initial_conditions, problem_pointers

   ! namelist parameters
   integer(kind=4) :: order   !< Order of mandelbrot set
   integer(kind=4) :: maxiter !< Maximum number of iterations
   logical :: smooth_map      !< Try continuous colouring
   logical :: log_polar       !< Use polar mapping around x_polar + i * y_polar
   real :: x_polar            !< x-coordinate for polar mode
   real :: y_polar            !< y-coordinate for polar mode
   real :: c_polar            !< correct colouring with x-coordinate multiplied by this factor
   real :: ref_thr            !< threshold for refining a grid

   namelist /PROBLEM_CONTROL/  order, maxiter, smooth_map, log_polar, x_polar, y_polar, c_polar, ref_thr

   ! other private data
   character(len=dsetnamelen), parameter :: mand_n = "mand", re_n = "real", imag_n = "imag"

contains

!> Set up some user hooks

   subroutine problem_pointers

      use dataio_user, only: user_reg_var_restart
#ifdef HDF5
      use dataio_user, only: user_vars_hdf5
#endif /* HDF5 */
      use user_hooks,  only: problem_refine_derefine

      implicit none

      user_reg_var_restart    => register_user_var
#ifdef HDF5
      user_vars_hdf5          => mand_vars
#endif /* HDF5 */
      problem_refine_derefine => mark_the_set

   end subroutine problem_pointers

!> \brief Read runtime parameters

   subroutine read_problem_par

      use bcast,      only: piernik_MPI_Bcast
      use constants,  only: ydim, LO, HI, dpi
      use dataio_pub, only: warn, die, nh
      use domain,     only: dom
      use mpisetup,   only: lbuff, ibuff, rbuff, master, slave

      implicit none

      ! namelist default parameter values
      order = 2
      maxiter = 100
      smooth_map = .true.
      log_polar = .false.
      x_polar = 0.
      y_polar = 0.
      c_polar = 0.
      ref_thr = 1.

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
         lbuff(2) = log_polar

         rbuff(1) = ref_thr
         rbuff(3) = x_polar
         rbuff(4) = y_polar
         rbuff(5) = c_polar

      endif

      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(lbuff)
      call piernik_MPI_Bcast(rbuff)

      if (slave) then

         order      = ibuff(1)
         maxiter    = ibuff(2)

         smooth_map = lbuff(1)
         log_polar  = lbuff(2)

         ref_thr    = rbuff(1)
         x_polar    = rbuff(3)
         y_polar    = rbuff(4)
         c_polar    = rbuff(5)

      endif

      if (any(dom%has_dir(:) .neqv. [ .true., .true., .false. ])) &
           call die("[initproblem:read_problem_par] Mandelbrot is supposed to by run only with XY plane and without Z-direction present")

      if (order /= 2) then
         if (master) call warn("[initproblem:read_problem_par] Only order == 2 is supported at the moment")
         order = 2
      endif

      if (log_polar .and. master) then
         if (dom%edge(ydim, HI) - dom%edge(ydim, LO) < 0.999*dpi) call warn("[initproblem:read_problem_par] not covering full angle")
         if (dom%edge(ydim, HI) - dom%edge(ydim, LO) > 1.001*dpi) call warn("[initproblem:read_problem_par] covering more than full angle")
      endif

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
      real, parameter :: bailout2 = 10., min_log_mand = 0.1
      real :: zx, zy, zt, cx, cy, rnit, r, f

      ! Create the initial density arrays
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         if (.not. cg%is_old) then

            call cg%set_constant_b_field([0., 0., 0.])
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
                  do i = cg%is, cg%ie
                     if (log_polar) then
                        r = 10**cg%x(i)
                        f = cg%y(j)
                        cx = x_polar + r*cos(f)
                        cy = y_polar + r*sin(f)
                     else
                        cx = cg%x(i)
                        cy = cg%y(j)
                     endif

                     zx = cx
                     zy = cy
                     nit = 1
                     do while (zx*zx + zy*zy < bailout2 .and. nit < maxiter)
                        zt = zx*zx - zy*zy + cx
                        zy = 2*zx*zy + cy
                        zx = zt
                        nit = nit + 1
                     enddo

                     rnit = nit
                     if (smooth_map .and. zx*zx + zy*zy > bailout2) rnit = rnit + 1 - log(log(sqrt(zx*zx + zy*zy)))/log(2.)

                     if (nit >= maxiter) then
                        mand(i, j, k) = min_log_mand ! increase contrast between interior and exterior
                     else
                        mand(i, j, k) = max(min_log_mand, log(max(rnit, min_log_mand)) + c_polar * cg%x(i))
                     endif
                     r__l(i, j, k) = zx
                     imag(i, j, k) = zy

                  enddo
               enddo
            enddo

         endif
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

      character(len=*),              intent(in)    :: var
      real, dimension(:,:,:),        intent(inout) :: tab
      integer,                       intent(inout) :: ierrh
      type(grid_container), pointer, intent(in)    :: cg

      ierrh = 0
      select case (trim(var))
         case ("distance", "dist") ! Supply the alternative name to comply with the old 4-letter limit
            tab(:,:,:) = log(sqrt(cg%q(qna%ind(re_n))%span(cg%ijkse)**2 + cg%q(qna%ind(imag_n))%span(cg%ijkse)**2))
         case ("angle", "ang")
            tab(:,:,:) = atan2(cg%q(qna%ind(imag_n))%span(cg%ijkse), cg%q(qna%ind(re_n))%span(cg%ijkse))
         case default
            ierrh = -1
      end select

   end subroutine mand_vars

!> \brief Mark interesting blocks for refinement

   subroutine mark_the_set

      use cg_list,    only: cg_list_element
      use cg_leaves,  only: leaves
      use domain,     only: dom
      use fluidindex, only: iarr_all_dn

      implicit none

      type(cg_list_element), pointer :: cgl
      real :: nitd, diffmax
      integer(kind=4) :: i, j, k

      cgl => leaves%first
      do while (associated(cgl))
         associate (cg => cgl%cg)
            ! Cannot use mand_n as long as it stays uninitialized during second call to problem_refine_derefine in update_refinement
            ! wna%fi is vital and thus automagilaclly prolonged
            ! Possible fixes:
            ! * make mand_n vital
            ! * do not call problem_refine_derefine twice in update_refinement
!!$            nitd = exp(maxval(cg%q(qna%ind(mand_n))%span(cg%ijkse), mask=cg%leafmap)) - &
!!$                 & exp(minval(cg%q(qna%ind(mand_n))%span(cg%ijkse), mask=cg%leafmap))

            diffmax = -huge(1.)
            ! Look one cell beyond boundary to prevent unnecessary derefinements
            do i = cg%is-dom%D_x, cg%ie+dom%D_x
               do j = cg%js-dom%D_y, cg%je+dom%D_y
                  do k = cg%ks-dom%D_z, cg%ke+dom%D_z
                     nitd = maxval(abs(exp(cg%u(iarr_all_dn(1), i, j, k)) - [ &
                          &            exp(cg%u(iarr_all_dn(1), i+dom%D_x, j, k)), &
                          &            exp(cg%u(iarr_all_dn(1), i-dom%D_x, j, k)), &
                          &            exp(cg%u(iarr_all_dn(1), i, j+dom%D_y, k)), &
                          &            exp(cg%u(iarr_all_dn(1), i, j-dom%D_y, k)), &
                          &            exp(cg%u(iarr_all_dn(1), i, j, k+dom%D_z)), &
                          &            exp(cg%u(iarr_all_dn(1), i, j, k-dom%D_z)) ] ) )
                     if (nitd >= ref_thr) call cg%flag%set(i, j, k)
                     diffmax = max(diffmax, nitd)
                  enddo
               enddo
            enddo

         end associate
         cgl => cgl%nxt
      enddo

   end subroutine mark_the_set

end module initproblem
