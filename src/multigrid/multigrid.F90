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
!! <tr><td>level_max       </td><td>1      </td><td>integer value  </td><td>\copydoc multigridvars::level_max          </td></tr>
!! <tr><td>ord_prolong     </td><td>0      </td><td>integer value  </td><td>\copydoc multigridvars::ord_prolong        </td></tr>
!! <tr><td>ord_prolong_face</td><td>0      </td><td>integer value  </td><td>\copydoc multigridvars::ord_prolong_face   </td></tr>
!! <tr><td>stdout          </td><td>.false.</td><td>logical        </td><td>\copydoc multigridvars::stdout             </td></tr>
!! <tr><td>verbose_vcycle  </td><td>.false.</td><td>logical        </td><td>\copydoc multigridvars::verbose_vcycle     </td></tr>
!! <tr><td>do_ascii_dump   </td><td>.false.</td><td>logical        </td><td>\copydoc multigridhelpers::do_ascii_dump   </td></tr>
!! <tr><td>dirty_debug     </td><td>.false.</td><td>logical        </td><td>\copydoc multigridhelpers::dirty_debug     </td></tr>
!! <tr><td>multidim_code_3D</td><td>.false.</td><td>logical        </td><td>\copydoc multigridhelpers::multidim_code_3D</td></tr>
!! </table>
!! \n \n
!<
   subroutine init_multigrid

      use grid,                only: cg, D_x, D_y, D_z
      use multigridvars,       only: lvl, dom_lvl, level_max, level_min, level_gb, roof, base, gb, mg_nb, ngridvars, correction, &
           &                         is_external, ord_prolong, ord_prolong_face, stdout, verbose_vcycle, tot_ts, is_mg_uneven
      use mpi,                 only: MPI_INTEGER, MPI_LOGICAL, MPI_IN_PLACE, MPI_LOR
      use mpisetup,            only: comm, ierr, proc, master, slave, nproc, has_dir, buffer_dim, ibuff, lbuff, dom, eff_dim, geometry_type, is_uneven
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

      integer                          :: ierrh, div, idx, j
      logical, save                    :: frun = .true.          !< First run flag
      real                             :: mb_alloc               !< Allocation counter
      integer, dimension(3)            :: aerr                   !> \deprecated BEWARE: hardcoded magic integer. Update when you change number of simultaneous error checks

      namelist /MULTIGRID_SOLVER/ level_max, ord_prolong, ord_prolong_face, stdout, verbose_vcycle, do_ascii_dump, dirty_debug, multidim_code_3D

      if (code_progress < PIERNIK_INIT_ARRAYS) call die("[multigrid:init_multigrid] grid, geometry, constants or arrays not initialized") ! This check is too weak (geometry), arrays are required only for multigrid_gravity

      if (.not.frun) call die("[multigrid:init_multigrid] Called more than once.")
      frun = .false.

      ! Default values for namelist variables
      level_max         = 1
      ord_prolong       = 0
      ord_prolong_face  = 0
      ! May all the logical parameters be .false. by default
      stdout                 = .false.
      verbose_vcycle         = .false.
      do_ascii_dump          = .false.
      dirty_debug            = .false.
      multidim_code_3D       = .false.

      if (master) then

         diff_nml(MULTIGRID_SOLVER)

         ibuff(1) = level_max
         ibuff(2) = ord_prolong
         ibuff(3) = ord_prolong_face

         lbuff(1) = stdout
         lbuff(2) = verbose_vcycle
         lbuff(3) = do_ascii_dump
         lbuff(4) = dirty_debug
         lbuff(5) = multidim_code_3D

      endif

      call MPI_Bcast(ibuff, buffer_dim, MPI_INTEGER, 0, comm, ierr)
      call MPI_Bcast(lbuff, buffer_dim, MPI_LOGICAL, 0, comm, ierr)

      if (slave) then

         level_max        = ibuff(1)
         ord_prolong      = ibuff(2)
         ord_prolong_face = ibuff(3)

         stdout           = lbuff(1)
         verbose_vcycle   = lbuff(2)
         do_ascii_dump    = lbuff(3)
         dirty_debug      = lbuff(4)
         multidim_code_3D = lbuff(5)

      endif

      is_mg_uneven = is_uneven

      ngridvars = correction  !< 4 variables are required for basic use of the multigrid solver
      if (eff_dim < 1 .or. eff_dim > 3) call die("[multigrid:init_multigrid] Unsupported number of dimensions.")

      if (geometry_type == GEO_RPZ) multidim_code_3D = .true. ! temporarily

!> \todo Make array of subroutine pointers
#ifdef GRAV
      call init_multigrid_grav
#endif /* GRAV */
#ifdef COSM_RAYS
      call init_multigrid_diff
#endif /* COSM_RAYS */

      !! Sanity checks
      if (abs(ord_prolong) > 2*mg_nb) call die("[multigrid:init_multigrid] not enough guardcells for given prolongation operator order")
      if (allocated(lvl)) call die("[multigrid:init_multigrid] lvl already allocated")
      allocate(lvl(level_gb:level_max), stat=aerr(1))                                 ! level_gb = level_min-1 contains some global base level data
      if (aerr(1) /= 0) call die("[multigrid:init_multigrid] Allocation error: lvl")
      mb_alloc = size(lvl)

      if (allocated(dom_lvl)) call die("[multigrid:init_multigrid] dom_lvl already allocated")
      allocate(dom_lvl(level_min:level_max), stat=aerr(1))                             ! level_gb = level_min-1 contains some global base level data
      if (aerr(1) /= 0) call die("[multigrid:init_multigrid] Allocation error: dom_lvl")
      mb_alloc = mb_alloc + size(dom_lvl) !> \deprecated it should be multiplied by the sizeof(dom_lvl)/sizeof(real)

      dom_lvl(level_max) = dom
      dom_lvl(level_max)%nb = mg_nb
      do idx = level_max-1, level_min, -1
         dom_lvl(idx) = dom_lvl(idx + 1)
         where (has_dir(:)) dom_lvl(idx)%n_d(:) = dom_lvl(idx)%n_d(:) / 2
         if (any(dom_lvl(idx)%n_d(:)*2 /= dom_lvl(idx + 1)%n_d(:) .and. has_dir(:))) then
            write(msg, '(a,3f10.1)')"[multigrid:init_multigrid] Fractional number of domain cells: ", 0.5*dom_lvl(idx + 1)%n_d(:)
            call die(msg)
         endif
         if (.not. allocated(dom_lvl(idx)%se)) then
            if (master) call warn("[multigrid:init_multigrid] dom_lvl%se not allocated implicitly")! this should detect changes in compiler behavior
            allocate(dom_lvl(idx)%se(0:nproc-1, xdim:zdim, LO:HI))
            !else dom_lvl%se was implicitly allocated during assignment
         endif
         dom_lvl(idx)%se(:, :, LO) = (dom_lvl(idx+1)%se(:, :, LO) + 1) / 2
         dom_lvl(idx)%se(:, :, HI) =  dom_lvl(idx+1)%se(:, :, HI)      / 2
         call dom_lvl(idx)%set_derived
      enddo
      if (any(dom%se(:,:,:) /= dom_lvl(level_max)%se(:,:,:))) call die("[multigrid:init_multigrid] dom%se or dom_lvl(level_max)%se corrupted") ! this should detect changes in compiler behavior

      !! Initialization of all regular levels (all but global base)
      !! Following loop gives us:
      !!    * shape (lvl(level_max  )) = (cg%nxb  , cg%nyb  , cg%nzb  ) + (2*cg%nb, 2*cg%nb, 2*cg%nb)
      !!    * shape (lvl(level_max-1)) = (cg%nxb/2, cg%nyb/2, cg%nzb/2) + (2*cg%nb, 2*cg%nb, 2*cg%nb)
      !!    * shape (lvl(level_max-2)) = (cg%nxb/4, cg%nyb/4, cg%nzb/4) + (2*cg%nb, 2*cg%nb, 2*cg%nb)
      !!    * ...
      !!    * shape (lvl(1)) = (cg%nxb/2**(level_max-1), cg%nyb/2**(level_max-1), cg%nzb/2**(level_max-1)) + (2*cg%nb, 2*cg%nb, 2*ncg%b)
      do idx = level_max, level_min, -1

         call lvl(idx)%init(dom_lvl(idx))

         lvl(idx)%level = idx                                      ! level number

         div = 2**(level_max -idx)                                 ! derefinement factor with respect to the top level

         if (any(lvl(idx)%n_b(:) < lvl(idx)%nb .and. has_dir(:))) then
            write(msg, '(a,i1,a,3i4,2(a,i2))')"[multigrid:init_multigrid] Number of guardcells exceeds number of interior cells: ", &
                 lvl(idx)%nb, " > ", lvl(idx)%n_b(:), " at level ", idx, ". You may try to set level_max <=", level_max-idx
            call die(msg)
         endif

         if (any(lvl(idx)%n_b(:) * div /= lvl(level_max)%n_b(:) .and. has_dir(:))) then
            write(msg, '(a,3f6.1,2(a,i2))')"[multigrid:init_multigrid] Fractional number of cells in: ", &
                 lvl(level_max)%n_b(:)/real(div), " at level ", idx, ". You may try to set level_max <=", level_max-idx
            call die(msg)
            is_mg_uneven = .true.
         endif

         lvl(idx)%vol  = 1.

         !> \todo check if these are correctly defined for multipole solver
         lvl(idx)%dxy = 1.
         lvl(idx)%dxz = 1.
         lvl(idx)%dyz = 1.

         if (has_dir(xdim)) then
            lvl(idx)%idx2  = 1. / lvl(idx)%dx**2                      ! auxiliary invariants
            lvl(idx)%vol   = lvl(idx)%vol * (cg%xmaxb-cg%xminb)
            lvl(idx)%dxy   = lvl(idx)%dxy * lvl(idx)%dx
            lvl(idx)%dxz   = lvl(idx)%dxz * lvl(idx)%dx
         else
            lvl(idx)%idx2  = 0.
         endif

         if (has_dir(ydim)) then
            lvl(idx)%idy2  = 1. / lvl(idx)%dy**2
            lvl(idx)%vol   = lvl(idx)%vol * (cg%ymaxb-cg%yminb)
            lvl(idx)%dxy   = lvl(idx)%dxy * lvl(idx)%dy
            lvl(idx)%dyz   = lvl(idx)%dyz * lvl(idx)%dy
         else
            lvl(idx)%idy2  = 0.
         endif

         if (has_dir(zdim)) then
            lvl(idx)%idz2  = 1. / lvl(idx)%dz**2
            lvl(idx)%vol   = lvl(idx)%vol * (cg%zmaxb-cg%zminb)
            lvl(idx)%dxz   = lvl(idx)%dxz * lvl(idx)%dz
            lvl(idx)%dyz   = lvl(idx)%dyz * lvl(idx)%dz
         else
            lvl(idx)%idz2  = 0.
         endif

         lvl(idx)%dvol2 = lvl(idx)%dvol**2

         ! data storage
         !> \deprecated BEWARE prolong_x and %prolong_xy are used only with RBGS relaxation when ord_prolong /= 0
         if ( allocated(lvl(idx)%prolong_x) .or. allocated(lvl(idx)%prolong_xy) .or. allocated(lvl(idx)%mgvar) ) call die("[multigrid:init_multigrid] multigrid arrays already allocated")
         allocate( lvl(idx)%mgvar     (lvl(idx)%nx, lvl(idx)%ny,                  lvl(idx)%nz,                  ngridvars), stat=aerr(1) )
         allocate( lvl(idx)%prolong_x (lvl(idx)%nx, lvl(idx)%nyb/2+2*lvl(idx)%nb, lvl(idx)%nzb/2+2*lvl(idx)%nb),            stat=aerr(2) )
         allocate( lvl(idx)%prolong_xy(lvl(idx)%nx, lvl(idx)%ny,                  lvl(idx)%nzb/2+2*lvl(idx)%nb),            stat=aerr(3) )
         if (any(aerr(1:3) /= 0)) call die("[multigrid:init_multigrid] Allocation error: lvl(idx)%*")
         if ( .not. allocated(lvl(idx)%prolong_x) .or. .not. allocated(lvl(idx)%prolong_xy) .or. .not. allocated(lvl(idx)%mgvar) .or. &
              .not. allocated(lvl(idx)%x) .or. .not. allocated(lvl(idx)%y) .or. .not. allocated(lvl(idx)%z) ) &
              call die("[multigrid:init_multigrid] some multigrid arrays not allocated")
         mb_alloc  = mb_alloc + size(lvl(idx)%prolong_x) + size(lvl(idx)%prolong_xy) + size(lvl(idx)%mgvar) + size(lvl(idx)%x)  + size(lvl(idx)%y) + size(lvl(idx)%z)

         if ( allocated(lvl(idx)%bnd_x) .or. allocated(lvl(idx)%bnd_y) .or. allocated(lvl(idx)%bnd_z)) call die("[multigrid:init_multigrid] multigrid boundary arrays already allocated")
         allocate( lvl(idx)%bnd_x(lvl(idx)%js:lvl(idx)%je, lvl(idx)%ks:lvl(idx)%ke, LO:HI), stat=aerr(1) )
         allocate( lvl(idx)%bnd_y(lvl(idx)%is:lvl(idx)%ie, lvl(idx)%ks:lvl(idx)%ke, LO:HI), stat=aerr(2) )
         allocate( lvl(idx)%bnd_z(lvl(idx)%is:lvl(idx)%ie, lvl(idx)%js:lvl(idx)%je, LO:HI), stat=aerr(3) )
         if (any(aerr(1:3) /= 0)) call die("[multigrid:init_multigrid] Allocation error: lvl(idx)%bnd_?")
         mb_alloc  = mb_alloc + size(lvl(idx)%bnd_x) + size(lvl(idx)%bnd_y) + size(lvl(idx)%bnd_z)

         ! array initialization
         if (dirty_debug) then
            lvl(idx)%mgvar     (:, :, :, :) = dirtyH
            lvl(idx)%prolong_x (:, :, :)    = dirtyH
            lvl(idx)%prolong_xy(:, :, :)    = dirtyH
            lvl(idx)%bnd_x     (:, :, :)    = dirtyH
            lvl(idx)%bnd_y     (:, :, :)    = dirtyH
            lvl(idx)%bnd_z     (:, :, :)    = dirtyH
         else
            lvl(idx)%mgvar     (:, :, :, :) = 0.0 ! should not be necessary if dirty_debug shows nothing suspicious
         endif

         ! set up connections between levels
         lvl(idx)%i_rst%proc = proc
         lvl(idx)%i_rst%se(:,:) = reshape( [ lvl(idx)%is, lvl(idx)%js, lvl(idx)%ks, lvl(idx)%ie, lvl(idx)%je, lvl(idx)%ke ], shape(lvl(idx)%i_rst%se(:,:)) )
         lvl(idx)%i_rst%nextgrid => null()
         lvl(idx)%o_rst = lvl(idx)%i_rst

         if (idx<level_max) lvl(idx)%i_rst%nextgrid => lvl(idx+1)
         if (idx>level_min) lvl(idx)%o_rst%nextgrid => lvl(idx-1)
      enddo

      call MPI_Allreduce(MPI_IN_PLACE, is_mg_uneven, 1, MPI_LOGICAL, MPI_LOR, comm, ierr)
      if ((is_mg_uneven .or. is_uneven) .and. ord_prolong /= 0) then
         ord_prolong = 0
         if (master) call warn("[multigrid:init_multigrid] prolongation order /= injection not implemented on uneven domains.")
      endif

      call mpi_multigrid_prep

      ! handy shortcuts
      base => lvl(level_min)
      roof => lvl(level_max)
      gb   => lvl(level_gb) !> \todo remove this feature. New restriction and prolongation will work with lvl(level_min) nonempty only on master

      tot_ts = 0.

      ! mark external faces
      is_external(:, :) = .false.
      do j=xdim, zdim
         if (has_dir(j) .and. .not. dom%periodic(j)) then
            is_external(j, LO) = (cg%off(j) == 0)
            is_external(j, HI) = (cg%off(j)+cg%n_b(j) == dom%n_d(j))
         endif
      enddo

#ifdef GRAV
      call init_multigrid_grav_post(mb_alloc)
#endif /* !GRAV */
#ifdef COSM_RAYS
      call init_multigrid_diff_post(mb_alloc)
#endif /* COSM_RAYS */

      ! summary
      if (master) then
         write(msg, '(a,i2,a,3(i4,a),f6.1,a)')"[multigrid:init_multigrid] Initialized ", level_max, " levels, coarsest resolution [ ", &
            lvl(level_min)%nxb, ",", lvl(level_min)%nyb, ",", lvl(level_min)%nzb, " ] per processor, allocated", mb_alloc*8./1048576., "MiB" ! sizeof(double)/2.**20
         call mg_write_log(msg)
      endif

   end subroutine init_multigrid

!!$ ============================================================================
!>
!! \brief Deallocate, destroy, demolish ...
!>

   subroutine cleanup_multigrid

      use multigridvars,      only: lvl, dom_lvl, level_gb, level_min, level_max, tot_ts, gb_cartmap
      use mpisetup,           only: master, nproc, comm3d, ierr, has_dir
      use constants,          only: xdim, zdim, LO, HI, BND, BLK
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

      integer :: i, ib, d
      real, allocatable, dimension(:) :: all_ts

#ifdef GRAV
      call cleanup_multigrid_grav
#endif /* GRAV */
#ifdef COSM_RAYS
      call cleanup_multigrid_diff
#endif /* COSM_RAYS */

      if (allocated(dom_lvl)) then
         do i=level_min, level_max
            if (allocated(dom_lvl(i)%se))     deallocate(dom_lvl(i)%se)
         enddo
         deallocate(dom_lvl)
      endif
      if (allocated(lvl)) then
         do i=level_gb, level_max
            if (allocated(lvl(i)%prolong_xy)) deallocate(lvl(i)%prolong_xy)
            if (allocated(lvl(i)%prolong_x))  deallocate(lvl(i)%prolong_x)
            if (allocated(lvl(i)%mgvar))      deallocate(lvl(i)%mgvar)
            if (allocated(lvl(i)%bnd_x))      deallocate(lvl(i)%bnd_x)
            if (allocated(lvl(i)%bnd_y))      deallocate(lvl(i)%bnd_y)
            if (allocated(lvl(i)%bnd_z))      deallocate(lvl(i)%bnd_z)

            if (i >= level_min) then
               do ib = 1, lvl(i)%nb
                  do d = xdim, zdim
                     if (has_dir(d)) then
                        call MPI_Type_free(lvl(i)%mmbc(d, LO, BND, ib), ierr)
                        call MPI_Type_free(lvl(i)%mmbc(d, LO, BLK, ib), ierr)
                        call MPI_Type_free(lvl(i)%mmbc(d, HI, BLK, ib), ierr)
                        call MPI_Type_free(lvl(i)%mmbc(d, HI, BND, ib), ierr)
                     endif
                  enddo
               enddo
               call lvl(i)%cleanup
            endif

         enddo
         deallocate(lvl)
      endif

      if (allocated(gb_cartmap))  deallocate(gb_cartmap)

      if (allocated(all_ts)) deallocate(all_ts)
      allocate(all_ts(0:nproc-1))

      call MPI_Gather(tot_ts, 1, MPI_DOUBLE_PRECISION, all_ts, 1, MPI_DOUBLE_PRECISION, 0, comm3d, ierr)

      if (master) then
         write(msg, '(a,3(g11.4,a))')"[multigrid] Spent ", sum(all_ts)/nproc, " seconds in multigrid_solve_* (min= ",minval(all_ts)," max= ",maxval(all_ts),")."
         call mg_write_log(msg, .false.)
      endif

      if (allocated(all_ts)) deallocate(all_ts)

   end subroutine cleanup_multigrid

end module multigrid
