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

#if defined(__INTEL_COMPILER)
   !! \deprecated remove this clause as soon as Intel Compiler gets required
   !! features and/or bug fixes
   use cg_list  ! QA_WARN required for ifort-12.1.x
#endif /* __INTEL_COMPILER */

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
!! <tr><td>do_ascii_dump        </td><td>.false.</td><td>logical       </td><td>\copydoc global::do_ascii_dump               </td></tr>
!! <tr><td>dirty_debug          </td><td>.false.</td><td>logical       </td><td>\copydoc global::dirty_debug                 </td></tr>
!! </table>
!! \n \n
!<
   subroutine init_multigrid

      use cg_list,             only: cg_list_element
      use cg_list_global,      only: all_cg
      use cg_list_level,       only: cg_list_level_T, base_lev, finest, coarsest
      use constants,           only: PIERNIK_INIT_GRID, LO, HI, I_ONE, O_INJ, O_LIN, O_I2, refinement_factor, dirtyH, base_level_offset
      use dataio_pub,          only: msg, par_file, namelist_errh, compare_namelist, cmdl_nml, lun, ierrh  ! QA_WARN required for diff_nml
      use dataio_pub,          only: printinfo, warn, die, code_progress
      use cart_comm,           only: cdd
      use domain,              only: dom, is_uneven, minsize
      use global,              only: dirty_debug, do_ascii_dump
      use grid_cont,           only: grid_container
      use mpi,                 only: MPI_LOGICAL, MPI_IN_PLACE, MPI_LOR, MPI_COMM_NULL
      use mpisetup,            only: comm, mpi_err, master, slave, nproc, ibuff, lbuff, piernik_MPI_Bcast
      use multigridvars,       only: single_base, source_n, solution_n, defect_n, correction_n, source, solution, defect, correction, &
           &                         ord_prolong, ord_prolong_face_norm, ord_prolong_face_par, stdout, verbose_vcycle, tot_ts, is_mg_uneven
      use named_array,         only: qna
#ifdef GRAV
      use multigrid_gravity,   only: init_multigrid_grav, init_multigrid_grav_post
#endif /* GRAV */
#ifdef COSM_RAYS
      use multigrid_diffusion, only: init_multigrid_diff, init_multigrid_diff_post
#endif /* COSM_RAYS */

      implicit none

      integer, parameter    :: level_incredible = 50  !< Increase this value only if your base domain contains much more than 10^15 cells in any active direction ;-)
      integer(kind=4)       :: level_max              !< Maximum allowed levels of base grid coarsening

      integer(kind=4)       :: j
      logical, save         :: frun = .true.          !< First run flag
      type(cg_list_element), pointer :: cgl
      type(cg_list_level_T), pointer :: curl          !< current level (a pointer sliding along the linked list) and temporary level
      type(grid_container),  pointer :: cg            !< current grid container

      namelist /MULTIGRID_SOLVER/ level_max, ord_prolong, ord_prolong_face_norm, ord_prolong_face_par, stdout, verbose_vcycle, do_ascii_dump, dirty_debug

      if (code_progress < PIERNIK_INIT_GRID) call die("[multigrid:init_multigrid] grid, geometry, constants or arrays not initialized")
      ! This check is too weak (geometry), arrays are required only for multigrid_gravity

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

      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(lbuff)

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
         if (any((mod(base_lev%n_d(:), int(refinement_factor, kind=8)**(j+1)) /= 0 .or. base_lev%n_d(:)/refinement_factor**(j+1) < minsize(:)) .and. dom%has_dir(:))) exit
      enddo
      if (level_max > j) then
         if (master) then
            if (level_max /= level_incredible) call warn("[multigrid:init_multigrid] level_max is too big,")
            write(msg,'(a,i3)')"[multigrid:init_multigrid] Automatically set level_max = ",j
            call printinfo(msg)
         endif
         level_max = j
      endif

      curl => finest
      do while (associated(curl))

         if (curl%level_id <= -level_max) exit

         ! create coarser level:
         call curl%add_level(coarse = .true.)
         if (coarsest%level_id == -level_max .and. single_base) then
            call coarsest%add_patch(int(coarsest%n_d(:), kind=4), base_level_offset, n_pieces=I_ONE)
         else
            call coarsest%add_patch(int(coarsest%n_d(:), kind=4), base_level_offset)
         endif
         call coarsest%init_all_new_cg

         curl => curl%coarser
      enddo

      is_mg_uneven = is_uneven
      curl => finest
      do while (associated(curl))

         cgl => curl%first
         do while (associated(cgl))
            cg => cgl%cg

            if (any(cg%n_b(:) < dom%nb .and. dom%has_dir(:))) then
               write(msg, '(a,i1,a,3i4,2(a,i2))')"[multigrid:init_multigrid] Number of guardcells exceeds number of interior cells: ", &
                    dom%nb, " > ", cg%n_b(:), " at level ", curl%level_id, ". You may try to set level_max <=", -curl%level_id
               call die(msg)
            endif

            if (any(cg%n_b(:) * refinement_factor**(finest%level_id - curl%level_id) /= finest%first%cg%n_b(:) .and. dom%has_dir(:)) .and. &
                 (.not. associated(curl, coarsest) .or. .not. single_base)) is_mg_uneven = .true.

            ! data storage
            if ( allocated(cg%mg%bnd_x) .or. allocated(cg%mg%bnd_y) .or. allocated(cg%mg%bnd_z)) &
                 call die("[multigrid:init_multigrid] multigrid boundary arrays already allocated")
            allocate(cg%mg%bnd_x(cg%js:cg%je, cg%ks:cg%ke, LO:HI))
            allocate(cg%mg%bnd_y(cg%is:cg%ie, cg%ks:cg%ke, LO:HI))
            allocate(cg%mg%bnd_z(cg%is:cg%ie, cg%js:cg%je, LO:HI))

            ! array initialization
            if (dirty_debug) then
               cg%mg%bnd_x(:, :, :) = dirtyH
               cg%mg%bnd_y(:, :, :) = dirtyH
               cg%mg%bnd_z(:, :, :) = dirtyH
            endif

            if (.not. associated(cg%wa)) cg%wa => cg%q(qna%wai)%arr ! required for CR diffusion

            cgl => cgl%nxt
         enddo

         curl => curl%coarser ! descend until null() is encountered
      enddo
      call MPI_Allreduce(MPI_IN_PLACE, is_mg_uneven, I_ONE, MPI_LOGICAL, MPI_LOR, comm, mpi_err)

      if ((is_mg_uneven .or. is_uneven .or. cdd%comm3d == MPI_COMM_NULL) .and. ord_prolong /= O_INJ .and. master) &
           call warn("[multigrid:init_multigrid] prolongation order /= injection may not work correctly on uneven or noncartesian domains yet.")

      call all_cg%reg_var(source_n,     ord_prolong = ord_prolong, multigrid = .true.)
      call all_cg%reg_var(solution_n,   ord_prolong = ord_prolong, multigrid = .true.)
      call all_cg%reg_var(defect_n,     ord_prolong = ord_prolong, multigrid = .true.)
      call all_cg%reg_var(correction_n, ord_prolong = ord_prolong, multigrid = .true.)

      source     = qna%ind(source_n)
      solution   = qna%ind(solution_n)
      defect     = qna%ind(defect_n)
      correction = qna%ind(correction_n)

      call all_cg%set_dirty(source)
      call all_cg%set_dirty(solution)
      call all_cg%set_dirty(defect)
      call all_cg%set_dirty(correction)

      tot_ts = 0.

#ifdef GRAV
      call init_multigrid_grav_post
#endif /* !GRAV */
#ifdef COSM_RAYS
      call init_multigrid_diff_post
#endif /* COSM_RAYS */

      ! summary
      if (master) then
         write(msg, '(a,i2,a,3i4,a)')"[multigrid:init_multigrid] Initialized ", finest%level_id - coarsest%level_id, &
              &                      " coarse levels, coarse level resolution [ ", coarsest%n_d(:)," ]"
         call printinfo(msg)
      endif

   end subroutine init_multigrid

!> \brief Deallocate, destroy, demolish ...

   subroutine cleanup_multigrid

      use constants,           only: I_ONE
      use dataio_pub,          only: msg, printinfo
      use mpi,                 only: MPI_DOUBLE_PRECISION
      use mpisetup,            only: master, nproc, FIRST, LAST, comm, mpi_err
      use multigridvars,       only: tot_ts
#ifdef GRAV
      use multigrid_gravity,   only: cleanup_multigrid_grav
#endif /* GRAV */
#ifdef COSM_RAYS
      use multigrid_diffusion, only: cleanup_multigrid_diff
#endif /* COSM_RAYS */

      implicit none

      real, allocatable, dimension(:) :: all_ts

#ifdef GRAV
      call cleanup_multigrid_grav
#endif /* GRAV */
#ifdef COSM_RAYS
      call cleanup_multigrid_diff
#endif /* COSM_RAYS */

      if (allocated(all_ts)) deallocate(all_ts)
      allocate(all_ts(FIRST:LAST))

      call MPI_Gather(tot_ts, I_ONE, MPI_DOUBLE_PRECISION, all_ts, I_ONE, MPI_DOUBLE_PRECISION, FIRST, comm, mpi_err)

      if (master) then
         write(msg, '(a,3(g11.4,a))')"[multigrid] Spent ", sum(all_ts)/nproc, " seconds in multigrid_solve_* (min= ",minval(all_ts)," max= ",maxval(all_ts),")."
         call printinfo(msg, .false.)
      endif

      if (allocated(all_ts)) deallocate(all_ts)

   end subroutine cleanup_multigrid

end module multigrid
