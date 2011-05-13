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
!! <tr><td>level_max            </td><td>1      </td><td>integer value </td><td>\copydoc multigrid::level_max                </td></tr>
!! <tr><td>ord_prolong          </td><td>0      </td><td>integer value </td><td>\copydoc multigridvars::ord_prolong          </td></tr>
!! <tr><td>ord_prolong_face_norm</td><td>0      </td><td>integer value </td><td>\copydoc multigridvars::ord_prolong_face_norm</td></tr>
!! <tr><td>ord_prolong_face_par </td><td>0      </td><td>integer value </td><td>\copydoc multigridvars::ord_prolong_face_par </td></tr>
!! <tr><td>stdout               </td><td>.false.</td><td>logical       </td><td>\copydoc multigridvars::stdout               </td></tr>
!! <tr><td>verbose_vcycle       </td><td>.false.</td><td>logical       </td><td>\copydoc multigridvars::verbose_vcycle       </td></tr>
!! <tr><td>do_ascii_dump        </td><td>.false.</td><td>logical       </td><td>\copydoc multigridhelpers::do_ascii_dump     </td></tr>
!! <tr><td>dirty_debug          </td><td>.false.</td><td>logical       </td><td>\copydoc multigridhelpers::dirty_debug       </td></tr>
!! <tr><td>multidim_code_3D     </td><td>.false.</td><td>logical       </td><td>\copydoc multigridhelpers::multidim_code_3D  </td></tr>
!! </table>
!! \n \n
!<
   subroutine init_multigrid

      use grid,                only: cg, D_x, D_y, D_z
      use multigridvars,       only: lvl, plvl, roof, base, mg_nb, ngridvars, correction, single_base, &
           &                         is_external, ord_prolong, ord_prolong_face_norm, ord_prolong_face_par, stdout, verbose_vcycle, tot_ts, is_mg_uneven
      use mpi,                 only: MPI_INTEGER, MPI_LOGICAL, MPI_DOUBLE_PRECISION, MPI_IN_PLACE, MPI_LOR, MPI_MIN, MPI_MAX, MPI_COMM_NULL
      use mpisetup,            only: comm, comm3d, ierr, proc, master, slave, nproc, has_dir, buffer_dim, ibuff, lbuff, dom, eff_dim, geometry_type, is_uneven
      use multigridhelpers,    only: mg_write_log, dirtyH, do_ascii_dump, dirty_debug, multidim_code_3D
      use multigridmpifuncs,   only: mpi_multigrid_prep
      use dataio_pub,          only: warn, die, code_progress
      use constants,           only: PIERNIK_INIT_ARRAYS, xdim, ydim, zdim, GEO_RPZ, LO, HI
      use dataio_pub,          only: msg, par_file, namelist_errh, compare_namelist, cmdl_nml  ! QA_WARN required for diff_nml
#ifdef GRAV
      use multigrid_gravity,   only: init_multigrid_grav, init_multigrid_grav_post
#endif /* GRAV */
#ifdef COSM_RAYS
      use multigrid_diffusion, only: init_multigrid_diff, init_multigrid_diff_post
#endif /* COSM_RAYS */

      implicit none

      integer, parameter    :: level_min = 1          !< Base (coarsest) level number
      integer               :: level_max              !< Levels of multigrid refinement

      integer               :: ierrh, j
      logical, save         :: frun = .true.          !< First run flag
      real                  :: mb_alloc, min_m, max_m !< Allocation counter
      integer, dimension(3) :: aerr                   !> \deprecated BEWARE: hardcoded magic integer. Update when you change number of simultaneous error checks
      type(plvl), pointer   :: curl                   !> current level (a pointer sliding along the linked list)

      namelist /MULTIGRID_SOLVER/ level_max, ord_prolong, ord_prolong_face_norm, ord_prolong_face_par, stdout, verbose_vcycle, do_ascii_dump, dirty_debug, multidim_code_3D

      if (code_progress < PIERNIK_INIT_ARRAYS) call die("[multigrid:init_multigrid] grid, geometry, constants or arrays not initialized") ! This check is too weak (geometry), arrays are required only for multigrid_gravity

      if (.not.frun) call die("[multigrid:init_multigrid] Called more than once.")
      frun = .false.

      ! Default values for namelist variables
      level_max             = 1
      ord_prolong           = 0
      ord_prolong_face_norm = 1
      ord_prolong_face_par  = 0
      ! May all the logical parameters be .false. by default
      stdout                = .false.
      verbose_vcycle        = .false.
      do_ascii_dump         = .false.
      dirty_debug           = .false.
      multidim_code_3D      = .false.

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
         lbuff(5) = multidim_code_3D

      endif

      call MPI_Bcast(ibuff, buffer_dim, MPI_INTEGER, 0, comm, ierr)
      call MPI_Bcast(lbuff, buffer_dim, MPI_LOGICAL, 0, comm, ierr)

      if (slave) then

         level_max             = ibuff(1)
         ord_prolong           = ibuff(2)
         ord_prolong_face_norm = ibuff(3)
         ord_prolong_face_par  = ibuff(4)

         stdout           = lbuff(1)
         verbose_vcycle   = lbuff(2)
         do_ascii_dump    = lbuff(3)
         dirty_debug      = lbuff(4)
         multidim_code_3D = lbuff(5)

      endif

      ! mark external faces
      is_external(:, :) = .false.
      do j=xdim, zdim
         if (has_dir(j) .and. .not. dom%periodic(j)) then
            is_external(j, LO) = (cg%off(j) == 0)
            is_external(j, HI) = (cg%off(j)+cg%n_b(j) == dom%n_d(j))
         endif
      enddo

      is_mg_uneven = is_uneven
      mb_alloc = 0
      single_base = (nproc == 1)

      ngridvars = correction  !< 4 variables are required for basic use of the multigrid solver
      if (eff_dim < 1 .or. eff_dim > 3) call die("[multigrid:init_multigrid] Unsupported number of dimensions.")

      if (geometry_type == GEO_RPZ) multidim_code_3D = .true. ! temporarily

!> \todo Make an array of subroutine pointers
#ifdef GRAV
      call init_multigrid_grav
#endif /* GRAV */
#ifdef COSM_RAYS
      call init_multigrid_diff
#endif /* COSM_RAYS */

      !! Sanity checks
      if (abs(ord_prolong) > 2*mg_nb) call die("[multigrid:init_multigrid] not enough guardcells for given prolongation operator order")
      if (ord_prolong_face_norm < 0) then
         if (master) call warn("[multigrid:init_multigrid] ord_prolong_face_norm < 0 is not defined, defaulting to 0")
         ord_prolong_face_norm = 0
      endif

      if (allocated(lvl)) call die("[multigrid:init_multigrid] lvl already allocated")
      allocate(lvl(level_min:level_max), stat=aerr(1))
      if (aerr(1) /= 0) call die("[multigrid:init_multigrid] Allocation error: lvl")
      mb_alloc = mb_alloc + size(lvl)

      ! handy shortcuts
      base => lvl(level_min)
      roof => lvl(level_max)

      ! set up connections between levels
      base%level = level_min
      curl => base
      do while (associated(curl))
         if (.not. associated(curl, roof)) then
            curl%finer => lvl(curl%level + 1)
            curl%finer%level = curl%level + 1
         else
            curl%finer => null()
         endif

         if (.not. associated(curl, base)) then
            curl%coarser => lvl(curl%level - 1)
         else
            curl%coarser => null()
         endif

         curl => curl%finer
      enddo

      curl => roof
      do while (associated(curl))

         if (associated(curl, roof)) then
            curl%dom = dom       ! inherit roof from global domain
            curl%dom%nb = mg_nb  ! the multigrid solver relies typically on 2 guardcells (roof%dom%n[xyz]t is wrong here, but it is not used in multigrid)
         else
            ! set up decomposition of coarse levels
            curl%dom = curl%finer%dom
            where (has_dir(:)) curl%dom%n_d(:) = curl%dom%n_d(:) / 2
            if (any(curl%dom%n_d(:)*2 /= curl%finer%dom%n_d(:) .and. has_dir(:))) then
               write(msg, '(a,3f10.1)')"[multigrid:init_multigrid] Fractional number of domain cells: ", 0.5*curl%finer%dom%n_d(:)
               call die(msg) ! handling this would require coarse grids bigger than base grid
            endif
            if (.not. allocated(curl%dom%se)) then
               if (master) call warn("[multigrid:init_multigrid] curl%dom%se not allocated implicitly")! this should detect changes in compiler behavior
               allocate(curl%dom%se(0:nproc-1, xdim:zdim, LO:HI))
               !else curl%dom%se was implicitly allocated during assignment
            endif
            if (associated(curl, base) .and. single_base) then ! Use curl%dom%se(nproc, :, :) with care, because this is an antiparallel thing
               curl%dom%se(:,  :, LO) = 0
               curl%dom%se(0,  :, HI) = curl%dom%n_d(:)-1 ! put the whole base level on the master CPU (\todo try a random or the last one)
               curl%dom%se(1:, :, HI) = -1
            else
               curl%dom%se(:, :, LO) = (curl%finer%dom%se(:, :, LO) + 1) / 2
               curl%dom%se(:, :, HI) =  curl%finer%dom%se(:, :, HI)      / 2
            endif
            call curl%dom%set_derived ! fix what was inherited from finer level or dom and should be recalculated
         endif

         call curl%init(curl%dom)

         if (any(curl%n_b(:) < curl%nb .and. has_dir(:) .and. .not. curl%empty)) then
            write(msg, '(a,i1,a,3i4,2(a,i2))')"[multigrid:init_multigrid] Number of guardcells exceeds number of interior cells: ", &
                 curl%nb, " > ", curl%n_b(:), " at level ", curl%level, ". You may try to set level_max <=", roof%level-curl%level
            call die(msg)
         endif

         if (any(curl%n_b(:) * 2**(roof%level - curl%level) /= roof%n_b(:) .and. has_dir(:)) .and. .not. curl%empty .and. (.not. associated(curl, base) .or. .not. single_base)) is_mg_uneven = .true.

         curl%vol  = 1.

         !> \todo check if these are correctly defined for multipole solver
         curl%dxy = 1.
         curl%dxz = 1.
         curl%dyz = 1.

         if (has_dir(xdim)) then
            curl%idx2  = 1. / curl%dx**2                      ! auxiliary invariants
            curl%vol   = curl%vol * (cg%xmaxb-cg%xminb)
            curl%dxy   = curl%dxy * curl%dx
            curl%dxz   = curl%dxz * curl%dx
         else
            curl%idx2  = 0.
         endif

         if (has_dir(ydim)) then
            curl%idy2  = 1. / curl%dy**2
            curl%vol   = curl%vol * (cg%ymaxb-cg%yminb)
            curl%dxy   = curl%dxy * curl%dy
            curl%dyz   = curl%dyz * curl%dy
         else
            curl%idy2  = 0.
         endif

         if (has_dir(zdim)) then
            curl%idz2  = 1. / curl%dz**2
            curl%vol   = curl%vol * (cg%zmaxb-cg%zminb)
            curl%dxz   = curl%dxz * curl%dz
            curl%dyz   = curl%dyz * curl%dz
         else
            curl%idz2  = 0.
         endif

         curl%dvol2 = curl%dvol**2

         ! data storage
         !> \deprecated BEWARE prolong_x and %prolong_xy are used only with RBGS relaxation when ord_prolong /= 0
         if ( allocated(curl%prolong_x) .or. allocated(curl%prolong_xy) .or. allocated(curl%mgvar) ) call die("[multigrid:init_multigrid] multigrid arrays already allocated")
         allocate( curl%mgvar     (curl%nx, curl%ny,                  curl%nz,                  ngridvars), stat=aerr(1) )
         allocate( curl%prolong_x (curl%nx, curl%nyb/2+2*curl%nb, curl%nzb/2+2*curl%nb),            stat=aerr(2) )
         allocate( curl%prolong_xy(curl%nx, curl%ny,                  curl%nzb/2+2*curl%nb),            stat=aerr(3) )
         if (any(aerr(1:3) /= 0)) call die("[multigrid:init_multigrid] Allocation error: curl%*")
         if ( .not. allocated(curl%prolong_x) .or. .not. allocated(curl%prolong_xy) .or. .not. allocated(curl%mgvar) .or. &
              ((.not. allocated(curl%x) .or. .not. allocated(curl%y) .or. .not. allocated(curl%z)) .and. .not. curl%empty) ) &
              call die("[multigrid:init_multigrid] some multigrid arrays not allocated")
         mb_alloc  = mb_alloc + size(curl%prolong_x) + size(curl%prolong_xy) + size(curl%mgvar)

         if ( allocated(curl%bnd_x) .or. allocated(curl%bnd_y) .or. allocated(curl%bnd_z)) call die("[multigrid:init_multigrid] multigrid boundary arrays already allocated")
         allocate( curl%bnd_x(curl%js:curl%je, curl%ks:curl%ke, LO:HI), stat=aerr(1) )
         allocate( curl%bnd_y(curl%is:curl%ie, curl%ks:curl%ke, LO:HI), stat=aerr(2) )
         allocate( curl%bnd_z(curl%is:curl%ie, curl%js:curl%je, LO:HI), stat=aerr(3) )
         if (any(aerr(1:3) /= 0)) call die("[multigrid:init_multigrid] Allocation error: curl%bnd_?")
         mb_alloc  = mb_alloc + size(curl%bnd_x) + size(curl%bnd_y) + size(curl%bnd_z)

         ! array initialization
         if (dirty_debug) then
            curl%mgvar     (:, :, :, :) = dirtyH
            curl%prolong_x (:, :, :)    = dirtyH
            curl%prolong_xy(:, :, :)    = dirtyH
            curl%bnd_x     (:, :, :)    = dirtyH
            curl%bnd_y     (:, :, :)    = dirtyH
            curl%bnd_z     (:, :, :)    = dirtyH
         else
            curl%mgvar     (:, :, :, :) = 0.0 ! should not be necessary if dirty_debug shows nothing suspicious
         endif

         curl => curl%coarser ! descend until null() is encountered
      enddo

      if (any(dom%se(:,:,:) /= roof%dom%se(:,:,:))) call die("[multigrid:init_multigrid] dom%se or roof%dom%se corrupted") ! this should detect changes in compiler behavior

      call MPI_Allreduce(MPI_IN_PLACE, is_mg_uneven, 1, MPI_LOGICAL, MPI_LOR, comm, ierr)
      if ((is_mg_uneven .or. is_uneven .or. comm3d == MPI_COMM_NULL) .and. ord_prolong /= 0) then
         ord_prolong = 0
         if (master) call warn("[multigrid:init_multigrid] prolongation order /= injection not implemented on uneven or noncartesian domains yet.")
      endif

      call mpi_multigrid_prep

      tot_ts = 0.

#ifdef GRAV
      call init_multigrid_grav_post(mb_alloc)
#endif /* !GRAV */
#ifdef COSM_RAYS
      call init_multigrid_diff_post(mb_alloc)
#endif /* COSM_RAYS */

      ! summary
      mb_alloc = mb_alloc / 2.**17 ! sizeof(double)/2.**20
      call MPI_Allreduce(mb_alloc, min_m, 1, MPI_DOUBLE_PRECISION, MPI_MIN, comm, ierr)
      call MPI_Allreduce(mb_alloc, max_m, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm, ierr)
      if (master) then
         write(msg, '(a,i2,a,3i4,a,2(f6.1,a))')"[multigrid:init_multigrid] Initialized ", roof%level, " levels, coarse level resolution [ ", &
            base%dom%n_d(:)," ], allocated", min_m, " ..", max_m, "MiB"
         call mg_write_log(msg)
      endif

   end subroutine init_multigrid

!!$ ============================================================================
!>
!! \brief Deallocate, destroy, demolish ...
!>

   subroutine cleanup_multigrid

      use multigridvars,      only: lvl, plvl, base, tot_ts, tgt_list
      use mpisetup,           only: master, nproc, comm, ierr, has_dir
      use constants,          only: xdim, zdim, LO, HI, BND, BLK, INVALID
      use mpi,                only: MPI_DOUBLE_PRECISION
      use multigridhelpers,   only: mg_write_log
      use dataio_pub,         only: msg
#ifdef GRAV
      use multigrid_gravity,  only: cleanup_multigrid_grav
#endif /* GRAV */
#ifdef COSM_RAYS
      use multigrid_diffusion, only: cleanup_multigrid_diff
#endif /* COSM_RAYS */

      implicit none

      integer :: ib, d, g
      real, allocatable, dimension(:) :: all_ts
      integer, parameter :: nseg = 2
      type(tgt_list), dimension(nseg) :: io_tgt
      type(plvl), pointer   :: curl

#ifdef GRAV
      call cleanup_multigrid_grav
#endif /* GRAV */
#ifdef COSM_RAYS
      call cleanup_multigrid_diff
#endif /* COSM_RAYS */

      if (allocated(lvl)) then
         curl => base
         do while (associated(curl))
            if (allocated(curl%prolong_xy)) deallocate(curl%prolong_xy)
            if (allocated(curl%prolong_x))  deallocate(curl%prolong_x)
            if (allocated(curl%mgvar))      deallocate(curl%mgvar)
            if (allocated(curl%bnd_x))      deallocate(curl%bnd_x)
            if (allocated(curl%bnd_y))      deallocate(curl%bnd_y)
            if (allocated(curl%bnd_z))      deallocate(curl%bnd_z)
            if (allocated(curl%dom%se))     deallocate(curl%dom%se)
            io_tgt(1:nseg) = [ curl%f_tgt, curl%c_tgt ]
            do ib = 1, nseg
               if (allocated(io_tgt(ib)%seg)) then
                  do g = 1, ubound(io_tgt(ib)%seg, dim=1)
                     if (allocated(io_tgt(ib)%seg(g)%buf)) deallocate(io_tgt(ib)%seg(g)%buf)
                  enddo
                  deallocate(io_tgt(ib)%seg)
               endif
            enddo
            do ib = 1, curl%nb
               do d = xdim, zdim
                  if (has_dir(d)) then
                     if (curl%mmbc(d, LO, BND, ib) /= INVALID) call MPI_Type_free(curl%mmbc(d, LO, BND, ib), ierr)
                     if (curl%mmbc(d, LO, BLK, ib) /= INVALID) call MPI_Type_free(curl%mmbc(d, LO, BLK, ib), ierr)
                     if (curl%mmbc(d, HI, BLK, ib) /= INVALID) call MPI_Type_free(curl%mmbc(d, HI, BLK, ib), ierr)
                     if (curl%mmbc(d, HI, BND, ib) /= INVALID) call MPI_Type_free(curl%mmbc(d, HI, BND, ib), ierr)
                  endif
               enddo
            enddo
            call curl%cleanup
            curl => curl%finer
         enddo
         deallocate(lvl)
      endif

      if (allocated(all_ts)) deallocate(all_ts)
      allocate(all_ts(0:nproc-1))

      call MPI_Gather(tot_ts, 1, MPI_DOUBLE_PRECISION, all_ts, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)

      if (master) then
         write(msg, '(a,3(g11.4,a))')"[multigrid] Spent ", sum(all_ts)/nproc, " seconds in multigrid_solve_* (min= ",minval(all_ts)," max= ",maxval(all_ts),")."
         call mg_write_log(msg, .false.)
      endif

      if (allocated(all_ts)) deallocate(all_ts)

   end subroutine cleanup_multigrid

end module multigrid
