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

!!$ ============================================================================
!>
!! \brief Multigrid solver initialization and cleanup
!!
!! \details This module contains initialization subroutine that reads the MULTIGRID_SOLVER namelist, performs all required allocation and initializes variables and arrays.
!! The initialization subroutine also calls secondary initialization routines for self-gravity and cosmic ray diffussion when applicable.
!! This module contains a cleanup subroutine that does all deallocations. While this is not strictly required for correct execution of the multigrid solver, it helps a lot
!! to debug memory leaks.
!<
module multigrid
! pulled by MULTIGRID

   implicit none

   private
   public :: init_multigrid, cleanup_multigrid

contains

!!$ ============================================================================
!!
!! Initializations and cleanup
!!
!>
!! \brief Routine to set parameters values from namelist MULTIGRID_SOLVER
!!
!! \n \n
!! @b MULTIGRID_SOLVER
!! \n \n
!! <table border="+1">
!! <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!! <tr><td>level_max            </td><td>1      </td><td>integer value </td><td>\copydoc init_multigrid::level_max           </td></tr>
!! <tr><td>ord_prolong          </td><td>0      </td><td>integer value </td><td>\copydoc multigridvars::ord_prolong          </td></tr>
!! <tr><td>ord_prolong_face_norm</td><td>0      </td><td>integer value </td><td>\copydoc multigridvars::ord_prolong_face_norm</td></tr>
!! <tr><td>ord_prolong_face_par </td><td>0      </td><td>integer value </td><td>\copydoc multigridvars::ord_prolong_face_par </td></tr>
!! <tr><td>stdout               </td><td>.false.</td><td>logical       </td><td>\copydoc multigridvars::stdout               </td></tr>
!! <tr><td>verbose_vcycle       </td><td>.false.</td><td>logical       </td><td>\copydoc multigridvars::verbose_vcycle       </td></tr>
!! <tr><td>do_ascii_dump        </td><td>.false.</td><td>logical       </td><td>\copydoc multigridhelpers::do_ascii_dump     </td></tr>
!! <tr><td>dirty_debug          </td><td>.false.</td><td>logical       </td><td>\copydoc multigridhelpers::dirty_debug       </td></tr>
!! </table>
!! \n \n
!<
   subroutine init_multigrid

      use cg_list_lev,         only: cg_list_level, cg_list_patch
      use constants,           only: PIERNIK_INIT_GRID, LO, HI, LONG, I_TWO, I_ONE, half, O_INJ, O_LIN, O_I2
      use dataio_pub,          only: msg, par_file, namelist_errh, compare_namelist, cmdl_nml, lun, ierrh  ! QA_WARN required for diff_nml
      use dataio_pub,          only: printinfo, warn, die, code_progress
      use decomposition,       only: divide_domain!, deallocate_pse
      use domain,              only: dom, is_uneven
      use gc_list,             only: cg_list_element, all_cg
      use grid,                only: base_lev
      use grid_cont,           only: grid_container
      use mpi,                 only: MPI_INTEGER, MPI_LOGICAL, MPI_IN_PLACE, MPI_LOR, MPI_COMM_NULL
      use mpisetup,            only: comm, ierr, proc, master, slave, nproc, FIRST, buffer_dim, ibuff, lbuff
      use multigridhelpers,    only: dirtyH, do_ascii_dump, dirty_debug, set_dirty
      use multigridvars,       only: roof, base, single_base, source_n, solution_n, defect_n, correction_n, source, solution, defect, correction, &
           &                         ord_prolong, ord_prolong_face_norm, ord_prolong_face_par, stdout, verbose_vcycle, tot_ts, is_mg_uneven
      use types,               only: cdd
#ifdef GRAV
      use multigrid_gravity,   only: init_multigrid_grav, init_multigrid_grav_post
#endif /* GRAV */
#ifdef COSM_RAYS
      use multigrid_diffusion, only: init_multigrid_diff, init_multigrid_diff_post
#endif /* COSM_RAYS */

      implicit none

      integer, parameter    :: level_incredible = 50  !< Increase this value only if your base domain contains much more than 10^15 cells in any active direction ;-)
      integer               :: level_max              !< Maximum allowed levels of base grid coarsening

      integer               :: j, g
      logical, save         :: frun = .true.          !< First run flag
      type(cg_list_element), pointer :: cgl
      type(cg_list_level),   pointer :: curl, tmpl    !< current level (a pointer sliding along the linked list) and temporary level
      type(cg_list_patch) :: patch                    !< wrapper for current level (to be passed to divide domain)
      type(grid_container),  pointer :: cg            !< current grid container

      namelist /MULTIGRID_SOLVER/ level_max, ord_prolong, ord_prolong_face_norm, ord_prolong_face_par, stdout, verbose_vcycle, do_ascii_dump, dirty_debug

      if (code_progress < PIERNIK_INIT_GRID) call die("[multigrid:init_multigrid] grid, geometry, constants or arrays not initialized") ! This check is too weak (geometry), arrays are required only for multigrid_gravity

      if (.not.frun) call die("[multigrid:init_multigrid] Called more than once.")
      frun = .false.

      ! Default values for namelist variables
      level_max             = level_incredible
      ord_prolong           = O_INJ
      ord_prolong_face_norm = O_I2
      ord_prolong_face_par  = O_INJ
      ! May all the logical parameters be .false. by default
      stdout                = .false.
      verbose_vcycle        = .false.
      do_ascii_dump         = .false.
      dirty_debug           = .false.

      if (master) then

         diff_nml(MULTIGRID_SOLVER)

         ibuff(1) = level_max
         ibuff(2) = ord_prolong
         ibuff(3) = ord_prolong_face_norm
         ibuff(4) = ord_prolong_face_par

         lbuff(1) = stdout
         lbuff(2) = verbose_vcycle
         lbuff(3) = do_ascii_dump
         lbuff(4) = dirty_debug

      endif

      call MPI_Bcast(ibuff, buffer_dim, MPI_INTEGER, FIRST, comm, ierr)
      call MPI_Bcast(lbuff, buffer_dim, MPI_LOGICAL, FIRST, comm, ierr)

      if (slave) then

         level_max             = ibuff(1)
         ord_prolong           = int(ibuff(2), kind=4)
         ord_prolong_face_norm = int(ibuff(3), kind=4)
         ord_prolong_face_par  = int(ibuff(4), kind=4)

         stdout           = lbuff(1)
         verbose_vcycle   = lbuff(2)
         do_ascii_dump    = lbuff(3)
         dirty_debug      = lbuff(4)

      endif

      stdout = stdout .and. master

      single_base = (nproc == 1)

      if (dom%eff_dim < 1 .or. dom%eff_dim > 3) call die("[multigrid:init_multigrid] Unsupported number of dimensions.")

!! \todo Make an array of subroutine pointers
#ifdef GRAV
      call init_multigrid_grav
#endif /* GRAV */
#ifdef COSM_RAYS
      call init_multigrid_diff
#endif /* COSM_RAYS */

      !! Sanity checks
      if (ord_prolong_face_norm < O_INJ) then
         if (master) call warn("[multigrid:init_multigrid] ord_prolong_face_norm < 0 is not defined, defaulting to 0")
         ord_prolong_face_norm = O_INJ
      endif

      if (ord_prolong == -1) ord_prolong = O_LIN

      if (level_max <= 0) then
         if (master) call warn("[multigrid:init_multigrid] level_max < 1: solving on a single grid may be extremely slow")
         level_max = 0
      endif

      do j = 0, level_incredible
         if (any((mod(base_lev%n_d(:), 2_LONG**(j+1)) /= 0 .or. base_lev%n_d(:)/2**(j+1) < dom%nb) .and. dom%has_dir(:))) exit
      enddo
      if (level_max > j) then
         if (master) then
            if (level_max /= level_incredible) call warn("[multigrid:init_multigrid] level_max is too big,")
            write(msg,'(a,i3)')"[multigrid:init_multigrid] Automatically set level_max = ",j
            call printinfo(msg)
         endif
         level_max = j
      endif

      roof => base_lev

      curl => roof
      do while (associated(curl))

         base => curl

         if (curl%lev <= -level_max) exit

         ! create coarser level:
         allocate(tmpl)
         call tmpl%init         ! create an empty cg list
         tmpl%lev = curl%lev -1 ! set id (number)
         curl%coarser => tmpl   ! set up connectivity
         tmpl%finer => curl
         tmpl%coarser => null()

         ! set up decomposition of coarse levels
         where (dom%has_dir(:))
            tmpl%n_d(:) = tmpl%finer%n_d / I_TWO
         elsewhere
            tmpl%n_d(:) = I_ONE
         endwhere

         if (master .and. any(tmpl%n_d(:)*2 /= tmpl%finer%n_d(:) .and. dom%has_dir(:))) then
            write(msg, '(a,3f10.1,a,i3)')"[multigrid:init_multigrid] Fractional number of domain cells: ", half*tmpl%finer%n_d(:), " at level ",tmpl%lev
            call die(msg)
         endif

         if (tmpl%lev == -level_max .and. single_base) then
            call base_on_single(tmpl)
         else

            patch%list_level => tmpl
            patch%n_d = int(tmpl%n_d, kind=4)

            if (.not. divide_domain(patch)) then
               write(msg,'(a,i4)')"[multigrid:init_multigrid] Coarse domain decomposition failed at level ",tmpl%lev
               if (master) call warn(msg)
               if (single_base) then
                  call base_on_single(tmpl)
                  level_max = -tmpl%lev
               else
                  deallocate(tmpl)
                  nullify(curl%coarser)
                  exit
                  !> \warning possible memory leak here
               endif
            endif

         endif

         call tmpl%print_segments

         do g = lbound(tmpl%pse(proc)%sel(:,:,:), dim=1), ubound(tmpl%pse(proc)%sel(:,:,:), dim=1)
            call tmpl%init_new_cg(g)
            call all_cg%add(tmpl%last%cg)
         enddo

         curl => curl%coarser
      enddo

      is_mg_uneven = is_uneven
      curl => roof
      do while (associated(curl))

         cgl => curl%first
         do while (associated(cgl))
            cg => cgl%cg

            if (any(cg%n_b(:) < dom%nb .and. dom%has_dir(:))) then
               write(msg, '(a,i1,a,3i4,2(a,i2))')"[multigrid:init_multigrid] Number of guardcells exceeds number of interior cells: ", &
                    dom%nb, " > ", cg%n_b(:), " at level ", curl%lev, ". You may try to set level_max <=", -curl%lev
               call die(msg)
            endif

            if (any(cg%n_b(:) * 2**(roof%lev - curl%lev) /= roof%first%cg%n_b(:) .and. dom%has_dir(:)) .and. (.not. associated(curl, base) .or. .not. single_base)) is_mg_uneven = .true.

            ! data storage
            if ( allocated(cg%mg%bnd_x) .or. allocated(cg%mg%bnd_y) .or. allocated(cg%mg%bnd_z)) call die("[multigrid:init_multigrid] multigrid boundary arrays already allocated")
            allocate(cg%mg%bnd_x(cg%js:cg%je, cg%ks:cg%ke, LO:HI))
            allocate(cg%mg%bnd_y(cg%is:cg%ie, cg%ks:cg%ke, LO:HI))
            allocate(cg%mg%bnd_z(cg%is:cg%ie, cg%js:cg%je, LO:HI))

            ! array initialization
            if (dirty_debug) then
               cg%mg%bnd_x(:, :, :) = dirtyH
               cg%mg%bnd_y(:, :, :) = dirtyH
               cg%mg%bnd_z(:, :, :) = dirtyH
            endif

            if (.not. associated(cg%wa)) cg%wa => cg%q(all_cg%wai)%arr ! required for CR diffusion

            cgl => cgl%nxt
         enddo

         call curl%vertical_prep

         curl => curl%coarser ! descend until null() is encountered
      enddo
      call MPI_Allreduce(MPI_IN_PLACE, is_mg_uneven, I_ONE, MPI_LOGICAL, MPI_LOR, comm, ierr)

      if ((is_mg_uneven .or. is_uneven .or. cdd%comm3d == MPI_COMM_NULL) .and. ord_prolong /= O_INJ .and. master) &
           call warn("[multigrid:init_multigrid] prolongation order /= injection may not work correctly on uneven or noncartesian domains yet.")

      call all_cg%reg_var(source_n,     ord_prolong = ord_prolong, multigrid = .true.)
      call all_cg%reg_var(solution_n,   ord_prolong = ord_prolong, multigrid = .true.)
      call all_cg%reg_var(defect_n,     ord_prolong = ord_prolong, multigrid = .true.)
      call all_cg%reg_var(correction_n, ord_prolong = ord_prolong, multigrid = .true.)

      source     = all_cg%ind(source_n)
      solution   = all_cg%ind(solution_n)
      defect     = all_cg%ind(defect_n)
      correction = all_cg%ind(correction_n)

      call set_dirty(source)
      call set_dirty(solution)
      call set_dirty(defect)
      call set_dirty(correction)

      tot_ts = 0.

#ifdef GRAV
      call init_multigrid_grav_post
#endif /* !GRAV */
#ifdef COSM_RAYS
      call init_multigrid_diff_post
#endif /* COSM_RAYS */

      ! summary
      if (master) then
         write(msg, '(a,i2,a,3i4,a)')"[multigrid:init_multigrid] Initialized ", roof%lev - base%lev, " coarse levels, coarse level resolution [ ", base%n_d(:)," ]"
         call printinfo(msg)
      endif

   end subroutine init_multigrid

!> \brief Deallocate, destroy, demolish ...

   subroutine cleanup_multigrid

      use constants,           only: I_ONE
      use dataio_pub,          only: msg, printinfo
      use cg_list_lev,         only: cg_list_level
      use grid,                only: base_lev
      use mpi,                 only: MPI_DOUBLE_PRECISION
      use mpisetup,            only: master, nproc, FIRST, LAST, comm, ierr
      use multigridvars,       only: base, tot_ts
#ifdef GRAV
      use multigrid_gravity,   only: cleanup_multigrid_grav
#endif /* GRAV */
#ifdef COSM_RAYS
      use multigrid_diffusion, only: cleanup_multigrid_diff
#endif /* COSM_RAYS */

      implicit none

      integer :: g
      real, allocatable, dimension(:) :: all_ts
      type(cg_list_level),   pointer :: curl

#ifdef GRAV
      call cleanup_multigrid_grav
#endif /* GRAV */
#ifdef COSM_RAYS
      call cleanup_multigrid_diff
#endif /* COSM_RAYS */

      !! \todo move this loop to grid::cleanup_grid
      curl => base
      do while (associated(curl) .and. .not. associated(curl, base_lev))
         call curl%delete
         curl => curl%finer
         deallocate(curl%coarser)
         if (allocated(curl%pse)) then
            do g = FIRST, LAST
               deallocate(curl%pse(g)%sel)
            enddo
            deallocate(curl%pse)
         endif
      enddo

      if (allocated(all_ts)) deallocate(all_ts)
      allocate(all_ts(FIRST:LAST))

      call MPI_Gather(tot_ts, I_ONE, MPI_DOUBLE_PRECISION, all_ts, I_ONE, MPI_DOUBLE_PRECISION, FIRST, comm, ierr)

      if (master) then
         write(msg, '(a,3(g11.4,a))')"[multigrid] Spent ", sum(all_ts)/nproc, " seconds in multigrid_solve_* (min= ",minval(all_ts)," max= ",maxval(all_ts),")."
         call printinfo(msg, .false.)
      endif

      if (allocated(all_ts)) deallocate(all_ts)

   end subroutine cleanup_multigrid

!> \brief Put the whole base level on the master CPU (\todo try a random or the last one)

   subroutine base_on_single(tmpl)

      use cg_list_lev,   only: cg_list_level
      use constants,     only: LO, HI, I_ONE
      use decomposition, only: allocate_pse
      use mpisetup,      only: nproc, FIRST

      implicit none

      type(cg_list_level), pointer, intent(inout) :: tmpl

      integer, dimension(nproc) :: n_cg

      n_cg(:) = 0
      n_cg(1) = I_ONE
      call allocate_pse(tmpl, n_cg)
      tmpl%pse(FIRST)%sel(I_ONE,  :, LO) = 0
      tmpl%pse(FIRST)%sel(I_ONE,  :, HI) = tmpl%n_d(:)-1

   end subroutine base_on_single

end module multigrid
