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

!!$ ============================================================================
!>
!! \brief Multigrid solver initialization and cleanup
!!
!! \details This module contains initialization subroutine that reads the MULTIGRID_SOLVER namelist, performs all required allocation and initializes variables and arrays.
!! The initialization subroutine also calls secondary initialization routines for self-gravity and cosmic ray diffusion when applicable.
!! This module contains a cleanup subroutine that does all deallocations. While this is not strictly required for correct execution of the multigrid solver, it helps a lot
!! to debug memory leaks.
!<
module multigrid
! pulled by MULTIGRID

   implicit none

   private
   public :: multigrid_par, init_multigrid, cleanup_multigrid, init_multigrid_ext

   integer(kind=4), parameter :: level_incredible = 50  !< Should be enough for up to 10^15 cells in any active direction ;-)
   integer(kind=4)            :: level_depth            !< Maximum allowed levels of base grid coarsening

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
!! <tr><td>level_depth          </td><td>1      </td><td>integer value </td><td>\copydoc multigrid::init_multigrid::level_depth</td></tr>
!! <tr><td>ord_prolong          </td><td>0      </td><td>integer value </td><td>\copydoc multigridvars::ord_prolong          </td></tr>
!! <tr><td>verbose_vcycle       </td><td>.false.</td><td>logical       </td><td>\copydoc multigridvars::verbose_vcycle       </td></tr>
!! <tr><td>do_ascii_dump        </td><td>.false.</td><td>logical       </td><td>\copydoc global::do_ascii_dump               </td></tr>
!! <tr><td>dirty_debug          </td><td>.false.</td><td>logical       </td><td>\copydoc global::dirty_debug                 </td></tr>
!! </table>
!! The list is active while \b "MULTIGRID" is defined.
!! \n \n
!<
   subroutine multigrid_par

      use bcast,               only: piernik_MPI_Bcast
      use cg_list_global,      only: all_cg
      use constants,           only: PIERNIK_INIT_DOMAIN, O_LIN, O_D3, I_ZERO
      use dataio_pub,          only: warn, die, code_progress, nh
      use domain,              only: dom
      use global,              only: dirty_debug, do_ascii_dump, show_n_dirtys !< \warning: alien variables go to local namelist
      use mpisetup,            only: master, slave, nproc, ibuff, lbuff
      use multigridvars,       only: single_base, ord_prolong, verbose_vcycle, tot_ts, &
           &                         source_n, solution_n, defect_n, correction_n, source, solution, defect, correction
      use named_array_list,    only: qna
#ifdef SELF_GRAV
      use multigrid_gravity,   only: multigrid_grav_par
#endif /* SELF_GRAV */
#ifdef COSM_RAYS
      use multigrid_diffusion, only: multigrid_diff_par
#endif /* COSM_RAYS */

      implicit none

      logical, save         :: frun = .true.          !< First run flag

      namelist /MULTIGRID_SOLVER/ level_depth, ord_prolong, verbose_vcycle, do_ascii_dump, dirty_debug, show_n_dirtys

      if (code_progress < PIERNIK_INIT_DOMAIN) call die("[multigrid:multigrid_par] grid, geometry, constants or arrays not initialized")
      ! This check is too weak (geometry), arrays are required only for multigrid_gravity

      if (.not.frun) call die("[multigrid:multigrid_par] Called more than once.")
      frun = .false.

      ! Default values for namelist variables
      level_depth           = merge(level_incredible, I_ZERO, dom%eff_dim /= 0)
      ord_prolong           = O_D3
      show_n_dirtys         = 16
      ! May all the logical parameters be .false. by default
      verbose_vcycle        = .false.
      do_ascii_dump         = .false.
      dirty_debug           = .false.

      if (master) then

         if (.not.nh%initialized) call nh%init()
         open(newunit=nh%lun, file=nh%tmp1, status="unknown")
         write(nh%lun,nml=MULTIGRID_SOLVER)
         close(nh%lun)
         open(newunit=nh%lun, file=nh%par_file)
         nh%errstr=""
         read(unit=nh%lun, nml=MULTIGRID_SOLVER, iostat=nh%ierrh, iomsg=nh%errstr)
         close(nh%lun)
         call nh%namelist_errh(nh%ierrh, "MULTIGRID_SOLVER")
         read(nh%cmdl_nml,nml=MULTIGRID_SOLVER, iostat=nh%ierrh)
         call nh%namelist_errh(nh%ierrh, "MULTIGRID_SOLVER", .true.)
         open(newunit=nh%lun, file=nh%tmp2, status="unknown")
         write(nh%lun,nml=MULTIGRID_SOLVER)
         close(nh%lun)
         call nh%compare_namelist()

         ibuff(1) = level_depth
         ibuff(2) = ord_prolong
         ibuff(3) = show_n_dirtys

         lbuff(1) = verbose_vcycle
         lbuff(2) = do_ascii_dump
         lbuff(3) = dirty_debug

      endif

      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(lbuff)

      if (slave) then

         level_depth    = ibuff(1)
         ord_prolong    = int(ibuff(2), kind=4)
         show_n_dirtys  = ibuff(3)

         verbose_vcycle = lbuff(1)
         do_ascii_dump  = lbuff(2)
         dirty_debug    = lbuff(3)

      endif

      single_base = (nproc == 1)

#ifndef DEBUG
      if (dirty_debug .and. master) call warn("[multigrid:multigrid_par] dirty_debug is supposed to be set only in debugging runs. Remember to disable it in production runs")
#endif /* !DEBUG */
      if (dom%eff_dim < 0 .or. dom%eff_dim > 3) call die("[multigrid:multigrid_par] Unsupported number of dimensions.")

!! \todo Make an array of subroutine pointers
#ifdef SELF_GRAV
      call multigrid_grav_par
#endif /* SELF_GRAV */
#ifdef COSM_RAYS
      call multigrid_diff_par
#endif /* COSM_RAYS */

      !! Sanity checks
      if (ord_prolong == -1) ord_prolong = O_LIN

      tot_ts = 0.

      call all_cg%reg_var(source_n,     ord_prolong = ord_prolong, multigrid = .true.)
      call all_cg%reg_var(solution_n,   ord_prolong = ord_prolong, multigrid = .true.)
      call all_cg%reg_var(defect_n,     ord_prolong = ord_prolong, multigrid = .true.)
      call all_cg%reg_var(correction_n, ord_prolong = ord_prolong, multigrid = .true.)

      source     = qna%ind(source_n)
      solution   = qna%ind(solution_n)
      defect     = qna%ind(defect_n)
      correction = qna%ind(correction_n)

   end subroutine multigrid_par

   subroutine init_multigrid

      use cg_list,            only: cg_list_element
      use cg_level_base,      only: base
      use cg_level_coarsest,  only: coarsest
      use cg_level_connected, only: cg_level_connected_t
      use cg_level_finest,    only: finest
      use constants,          only: PIERNIK_INIT_GRID, I_ONE, refinement_factor, V_VERBOSE
      use dataio_pub,         only: printinfo, warn, die, code_progress, msg
      use domain,             only: dom, minsize
      use grid_cont,          only: grid_container
      use mpisetup,           only: master
      use multigridvars,      only: single_base
#ifdef SELF_GRAV
      use multigrid_gravity,  only: init_multigrid_grav
#endif /* SELF_GRAV */

      implicit none

      integer(kind=4)       :: j

      type(cg_list_element), pointer :: cgl
      type(cg_level_connected_t), pointer :: curl  !< current level (a pointer sliding along the linked list) and temporary level
      type(grid_container),  pointer :: cg         !< current grid container

      if (code_progress < PIERNIK_INIT_GRID) call die("[multigrid:init_multigrid] grid, geometry, constants or arrays not initialized")
      ! This check is too weak (geometry), arrays are required only for multigrid_gravity

      if (level_depth <= 0) then
         if (master) call warn("[multigrid:init_multigrid] level_depth < 1: solving on a single grid may be extremely slow")
         level_depth = 0
      endif

      do j = 0, level_incredible
         if (any((mod(base%level%l%n_d(:), int(refinement_factor, kind=8)**(j+1)) /= 0 .or. base%level%l%n_d(:)/refinement_factor**(j+1) < minsize(:)) .and. dom%has_dir(:))) exit
         if (any((mod(base%level%l%off(:), int(refinement_factor, kind=8)**(j+1)) /= 0 .and. dom%has_dir(:)))) exit
      enddo
      if (level_depth > j) then
         if (master) then
            if (level_depth /= level_incredible) call warn("[multigrid:init_multigrid] level_depth is too big,")
            write(msg,'(a,i3)')"[multigrid:init_multigrid] Automatically set level_depth = ",j
            call printinfo(msg, V_VERBOSE)
         endif
         level_depth = j
      endif

      do while (coarsest%level%l%id > -level_depth)
         call coarsest%add_coarser
      enddo

      curl => finest%level
      do while (associated(curl))

         cgl => curl%first
         do while (associated(cgl))
            cg => cgl%cg

            if (any(cg%n_b(:) < dom%nb .and. dom%has_dir(:))) then
               write(msg, '(a,i1,a,3i4,2(a,i2))')"[multigrid:init_multigrid] Number of guardcells exceeds number of interior cells: ", &
                    dom%nb, " > ", cg%n_b(:), " at level ", curl%l%id, ". You may try to set level_depth <=", -curl%l%id
               call die(msg)
            endif

            cgl => cgl%nxt
         enddo

         curl => curl%coarser ! descend until null() is encountered
      enddo

      curl => base%level%coarser
      do while (associated(curl))
         if (master) then
            if (curl%l%id == -level_depth .and. single_base) then
               call curl%add_patch(n_pieces=I_ONE)
            else
               !> \todo When there is more AMR::bsize-pieces than processes, consider forcing cartesian or noncartesian decomposition
               call curl%add_patch
            endif
         endif
         call curl%init_all_new_cg

         curl => curl%coarser
      enddo

#ifdef SELF_GRAV
      call init_multigrid_grav
#endif /* SELF_GRAV */

      ! summary
      if (master) then
         write(msg, '(a,i2,a,3i4,a)')"[multigrid:init_multigrid] Initialized ", finest%level%l%id - coarsest%level%l%id, &
              &                      " coarse levels, coarsest level resolution [ ", coarsest%level%l%n_d(:)," ]"
         call printinfo(msg, V_VERBOSE)
      endif

   end subroutine init_multigrid

!> set up pointers for cg%mg initialization

   subroutine init_multigrid_ext

      use constants,          only: dsetnamelen
      use grid_container_ext, only: cg_ext, cg_extptrs
#ifdef SELF_GRAV
      use multigrid_gravity,  only: init_multigrid_grav_ext
#endif /* SELF_GRAV */

      implicit none

      procedure(cg_ext), pointer :: mg_cg_init_p, mg_cg_cleanup_p
      character(len=dsetnamelen), parameter :: mg_ext_name = "multigrid"

      mg_cg_init_p => mg_cg_init
      mg_cg_cleanup_p => mg_cg_cleanup
      call cg_extptrs%extend(mg_cg_init_p, mg_cg_cleanup_p, mg_ext_name)

#ifdef SELF_GRAV
      call init_multigrid_grav_ext(mg_ext_name)
#endif /* SELF_GRAV */
!!$#ifdef COSM_RAYS
!!$      call init_multigrid_cr_ext(mg_ext_name)
!!$#endif /* COSM_RAYS */

   end subroutine init_multigrid_ext

!> \brief Allocate some multigrid-specific arrays

   subroutine mg_cg_init(cg)

      use constants, only: LO, HI
      use grid_cont, only: grid_container

      implicit none

      type(grid_container), pointer,  intent(inout) :: cg

      allocate(cg%mg)

      allocate(cg%mg%bnd_x(cg%js:cg%je, cg%ks:cg%ke, LO:HI))
      allocate(cg%mg%bnd_y(cg%is:cg%ie, cg%ks:cg%ke, LO:HI))
      allocate(cg%mg%bnd_z(cg%is:cg%ie, cg%js:cg%je, LO:HI))

   end subroutine mg_cg_init

!> \brief Deallocate what was allocated in mg_cg_init

   subroutine mg_cg_cleanup(cg)

      use grid_cont, only: grid_container

      implicit none

      type(grid_container), pointer,  intent(inout) :: cg

      if (allocated(cg%mg%bnd_x)) deallocate(cg%mg%bnd_x)
      if (allocated(cg%mg%bnd_y)) deallocate(cg%mg%bnd_y)
      if (allocated(cg%mg%bnd_z)) deallocate(cg%mg%bnd_z)
      deallocate(cg%mg)

   end subroutine mg_cg_cleanup

!> \brief Deallocate, destroy, demolish ...

   subroutine cleanup_multigrid

      use constants,           only: I_ONE, V_LOG
      use dataio_pub,          only: msg, printinfo
      use MPIF,                only: MPI_DOUBLE_PRECISION, MPI_COMM_WORLD
      use MPIFUN,              only: MPI_Gather
      use mpisetup,            only: master, nproc, FIRST, LAST, err_mpi
      use multigridvars,       only: tot_ts
#ifdef SELF_GRAV
      use multigrid_gravity,   only: cleanup_multigrid_grav
#endif /* SELF_GRAV */
#ifdef COSM_RAYS
      use multigrid_diffusion, only: cleanup_multigrid_diff
#endif /* COSM_RAYS */

      implicit none

      real, allocatable, dimension(:) :: all_ts

#ifdef SELF_GRAV
      call cleanup_multigrid_grav
#endif /* SELF_GRAV */
#ifdef COSM_RAYS
      call cleanup_multigrid_diff
#endif /* COSM_RAYS */

      if (allocated(all_ts)) deallocate(all_ts)
      allocate(all_ts(FIRST:LAST))

      call MPI_Gather(tot_ts, I_ONE, MPI_DOUBLE_PRECISION, all_ts, I_ONE, MPI_DOUBLE_PRECISION, FIRST, MPI_COMM_WORLD, err_mpi)

      if (master) then
         write(msg, '(a,3(g11.4,a))')"[multigrid] Spent ", sum(all_ts)/nproc, " seconds in multigrid_solve_* (min= ",minval(all_ts)," max= ",maxval(all_ts),")."
         call printinfo(msg, V_LOG)
      endif

      if (allocated(all_ts)) deallocate(all_ts)

   end subroutine cleanup_multigrid

end module multigrid
