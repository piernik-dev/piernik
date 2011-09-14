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
!! \brief Implicit diffusion multigrid solver
!!
!! \details This module contains routines and variables specific for multigrid diffusion solver.
!!
!! Note that for constant-coefficient problems one can use Poisson solver routines, including FFT solver.
!<

module multigrid_diffusion
! pulled by MULTIGRID && COSM_RAYS

   use multigridvars, only: correction, vcycle_stats
   use constants,     only: cbuff_len, ndims

   implicit none

   private
   public :: init_multigrid_diff, init_multigrid_diff_post, cleanup_multigrid_diff, multigrid_solve_diff
   public :: diff_tstep_fac, diff_explicit, diff_dt_crs_orig

   ! namelist parameters
   real               :: norm_tol                                     !< stop V-cycle iterations when the ratio of norms ||residual||/||source|| is below this value
   real               :: vcycle_abort                                 !< abort the V-cycle when lhs norm raises by this factor
   real               :: overrelax                                    !< overrealaxation factor (if < 1. then works as underrelaxation), use with care
   integer            :: max_cycles                                   !< Maximum allowed number of V-cycles
   integer            :: nsmool                                       !< smoothing cycles per call
   integer            :: nsmoob                                       !< smoothing cycles on base level;  \todo implement a convergence check
   real               :: diff_theta                                   !< 0. is explicit, 1. is fully implicit 0.5 is Crank-Nicholson
   real,    protected :: diff_tstep_fac                               !< How much we stretch timestep. Note that for diff_theta == 0. this should not be > 1.
   logical, protected :: diff_explicit                                !< If .true. then do not use multigrid for diffusion
   logical            :: allow_explicit                               !< When timestep is limited somewhere else, allow explicit calculation (should be a bit faster)
   real               :: diff_dt_crs_orig                             !< timestep calculated at timestepcosmicrays.F90, before enlarging by diff_tstep_fac
   character(len=cbuff_len) :: diff_bnd_str                           !< Type of diffusion boundary conditions. Can be "isolated", "reflecting" or "zero" (there are some aliases as well)

   ! mgvar entries for the B field
   enum, bind(C)
      enumerator :: diff_bx = correction+1                            !< index of B_x in the cg%b%arr(:,:,:,:) array
      enumerator :: diff_by                                           !< index of B_y in the cg%b%arr(:,:,:,:) array
      enumerator :: diff_bz                                           !< index of B_z in the cg%b%arr(:,:,:,:) array
   end enum
   integer, parameter, dimension(ndims) :: idiffb = [diff_bx, diff_by, diff_bz]

   ! miscellaneous
   logical, allocatable, dimension(:) :: norm_was_zero                !< Flag for suppressing repeated warnings on nonexistent CR components
   type(vcycle_stats) :: vstat                                        !< V-cycle statistics
   integer(kind=4) :: diff_extbnd                                     !< external boundary type for relaxation and computation of residuum

contains

!!$ ============================================================================
!!
!! Initialization
!!
!>
!! \brief Routine to set parameters values from namelist MULTIGRID_DIFFUSION
!!
!! \n \n
!! @b MULTIGRID_DIFFUSION
!! \n \n
!! <table border="+1">
!! <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!! <tr><td>norm_tol      </td><td>1.e-2  </td><td>real value     </td><td>\copydoc multigrid_diffusion::norm_tol      </td></tr>
!! <tr><td>vcycle_abort  </td><td>2.0    </td><td>real value     </td><td>\copydoc multigrid_diffusion::vcycle_abort  </td></tr>
!! <tr><td>max_cycles    </td><td>20     </td><td>integer value  </td><td>\copydoc multigrid_diffusion::max_cycles    </td></tr>
!! <tr><td>nsmool        </td><td>4      </td><td>integer value  </td><td>\copydoc multigrid_diffusion::nsmool        </td></tr>
!! <tr><td>nsmoob        </td><td>1      </td><td>integer value  </td><td>\copydoc multigrid_diffusion::nsmoob        </td></tr>
!! <tr><td>overrelax     </td><td>1.     </td><td>real value     </td><td>\copydoc multigrid_diffusion::overrelax     </td></tr>
!! <tr><td>diff_theta    </td><td>1.     </td><td>real value     </td><td>\copydoc multigrid_diffusion::diff_theta    </td></tr>
!! <tr><td>diff_tstep_fac</td><td>1.     </td><td>real value     </td><td>\copydoc multigrid_diffusion::diff_tstep_fac</td></tr>
!! <tr><td>diff_explicit </td><td>.false.</td><td>logical        </td><td>\copydoc multigrid_diffusion::diff_explicit </td></tr>
!! <tr><td>allow_explicit</td><td>.true. </td><td>logical        </td><td>\copydoc multigrid_diffusion::allow_explicit</td></tr>
!! <tr><td>diff_bnd_str  </td><td>"zero" </td><td>string of chars</td><td>\copydoc multigrid_diffusion::diff_bnd_str  </td></tr>
!! </table>
!! \n \n
!<
   subroutine init_multigrid_diff

      use constants,     only: GEO_XYZ, half, zero, one
      use dataio_pub,    only: par_file, ierrh, namelist_errh, compare_namelist, cmdl_nml, lun, getlun      ! QA_WARN required for diff_nml
      use dataio_pub,    only: die, warn, msg
      use domain,        only: geometry_type
      use mpi,           only: MPI_DOUBLE_PRECISION, MPI_INTEGER, MPI_LOGICAL, MPI_CHARACTER
      use mpisetup,      only: comm, ierr, master, slave, nproc, ibuff, rbuff, lbuff, cbuff, buffer_dim, FIRST
      use multigridvars, only: ngridvars, extbnd_zero, extbnd_extrapolate, extbnd_mirror, extbnd_antimirror, single_base

      implicit none

      logical, save                    :: frun = .true.          !< First run flag

      namelist /MULTIGRID_DIFFUSION/ norm_tol, vcycle_abort, max_cycles, nsmool, nsmoob, overrelax, &
           &                         diff_theta, diff_tstep_fac, diff_explicit, allow_explicit, diff_bnd_str

      if (.not.frun) call die("[multigrid_diffusion:init_multigrid_diff] Called more than once.")
      frun = .false.
      if (geometry_type /= GEO_XYZ) call die("[multigrid_gravity:init_multigrid_gravdiffusion:init_multigrid_diff] non-cartesian geometry not implemented yet.")

      ! Default values for namelist variables
      norm_tol       = 1.e-2
      vcycle_abort   = 2.    ! unused as yet
      diff_theta     = 1.
      diff_tstep_fac = 1.
      overrelax      = 1.
      max_cycles     = 20
      nsmool         = 4
      nsmoob         = 1
      diff_explicit  = .false.
      allow_explicit = .true.
      diff_bnd_str   = "zero"

      if (master) then

         diff_nml(MULTIGRID_DIFFUSION)

         rbuff(1) = norm_tol
         rbuff(2) = vcycle_abort
         rbuff(3) = diff_theta
         rbuff(4) = diff_tstep_fac
         rbuff(5) = overrelax

         ibuff(1) = max_cycles
         ibuff(2) = nsmool
         ibuff(3) = nsmoob

         lbuff(1) = diff_explicit
         lbuff(2) = allow_explicit

         cbuff(1) = diff_bnd_str

      endif

      call MPI_Bcast(cbuff, cbuff_len*buffer_dim, MPI_CHARACTER,        FIRST, comm, ierr)
      call MPI_Bcast(ibuff,           buffer_dim, MPI_INTEGER,          FIRST, comm, ierr)
      call MPI_Bcast(rbuff,           buffer_dim, MPI_DOUBLE_PRECISION, FIRST, comm, ierr)
      call MPI_Bcast(lbuff,           buffer_dim, MPI_LOGICAL,          FIRST, comm, ierr)

      if (slave) then

         norm_tol       = rbuff(1)
         vcycle_abort   = rbuff(2)
         diff_theta     = rbuff(3)
         diff_tstep_fac = rbuff(4)
         overrelax      = rbuff(5)

         max_cycles     = ibuff(1)
         nsmool         = ibuff(2)
         nsmoob         = ibuff(3)

         diff_explicit  = lbuff(1)
         allow_explicit = lbuff(2)

         diff_bnd_str   = cbuff(1)(1:len(diff_bnd_str))

      endif

      ngridvars = max(ngridvars, diff_bz)

      ! boundaries
      diff_extbnd = extbnd_zero
      select case (diff_bnd_str)
         case ("isolated", "iso", "free")
            diff_extbnd = extbnd_extrapolate
         case ("reflecting", "refl", "styrofoam")
            diff_extbnd = extbnd_mirror
         case ("zero", "cold", "antireflecting")
            diff_extbnd = extbnd_antimirror
         case default
            write(msg,'(3a)')"[multigrid_diffusion:init_multigrid_diff] Non-recognized boundary description '",diff_bnd_str,"'"
            call die(msg)
      end select

      !diffusion
      if (.not. diff_explicit) then
         if (diff_theta < zero .or. diff_theta > one) call die("[multigrid_diffusion:init_multigrid] diff_theta must lie in the range [0. .. 1.]")
         if (diff_theta < half .and. diff_tstep_fac > one .and. master) call warn("[multigrid_diffusion:init_multigrid] diff_tstep_fac > 1. for diff_theta < 0.5 might be unstable")
         ! calculate exact limit formula
         ! for diff_theta=0. stable diff_tstep_fac is 0.5 in 2D (guess: 0.333 in 3D)
         ! for diff_theta<0.5 stable diff_tstep_fac rises by 1./(1.-2.*diff_theta)
      endif

      if (overrelax /= 1 .and. master) then
         write(msg, '(a,f8.5)')"[multigrid_diffusion:init_multigrid_diff] Overrelaxation factor = ", overrelax
         call warn(msg)
      endif

      if (single_base .and. nproc > 1) then
         call warn("[multigrid_diffusion:init_multigrid_diff] single_base disabled just in case")
         single_base = (nproc == 1)
      endif

   end subroutine init_multigrid_diff

!!$ ============================================================================
!>
!! \brief Initialization - continued after allocation of everything interesting
!<

   subroutine init_multigrid_diff_post(mb_alloc)

      use dataio_pub,         only: die
      use fluidindex,         only: flind
      use multigridhelpers,   only: vcycle_stats_init
      use multigridvars,      only: vcycle_stats

      implicit none

      real, intent(inout)              :: mb_alloc               !< Allocation counter

      if (allocated(norm_was_zero)) call die("[multigrid_diffusion:init_multigrid] norm_was_zero already allocated")
      allocate(norm_was_zero(flind%crs%all))
      mb_alloc = mb_alloc + size(norm_was_zero)/2
      norm_was_zero(:) = .false.

      call vcycle_stats_init(vstat, max_cycles)
      mb_alloc = mb_alloc + 2*max_cycles

   end subroutine init_multigrid_diff_post

!!$ ============================================================================
!>
!! \brief Cleanup
!<

   subroutine cleanup_multigrid_diff

      implicit none

      if (allocated(vstat%factor)) deallocate(vstat%factor)
      if (allocated(vstat%time)) deallocate(vstat%time)
      if (allocated(norm_was_zero)) deallocate(norm_was_zero)

   end subroutine cleanup_multigrid_diff

!!$ ============================================================================
!>
!! \brief Multigrid diffusion driver. This is the only multigrid routine intended to be called from the fluidupdate module.
!! This routine is also responsible for communicating the solution to the rest of world
!<

   subroutine multigrid_solve_diff

      use constants,     only: xdim, ydim, zdim
      use crdiffusion,   only: cr_diff
      use dataio_pub,    only: halfstep, warn, printinfo, msg
      use fluidindex,    only: flind, ibx, iby, ibz
      use global,        only: dt
      use mpisetup,      only: master
      use multigridvars, only: ts, tot_ts, stdout
      use timer,         only: set_timer

      implicit none

      logical, save :: frun = .true.
      integer       :: cr_id         ! maybe we should make this variable global in the module and do not pass it as an argument?

      ts =  set_timer("multigrid_diffusion", .true.)

      if (diff_explicit .or. (allow_explicit .and. dt/diff_dt_crs_orig<1)) then

         if (frun) then
            if (master .and. diff_explicit) call warn("[multigrid_diffusion:multigrid_solve_diff] Multigrid was initialized but is not used")
            frun = .false.
         endif
         if (halfstep) then
            call cr_diff(zdim,ibz)
            call cr_diff(ydim,iby)
            call cr_diff(xdim,ibx)
         else
            call cr_diff(xdim,ibx)
            call cr_diff(ydim,iby)
            call cr_diff(zdim,ibz)
         endif

      else

         if (dt < 0.99999 * diff_dt_crs_orig * diff_tstep_fac .and. .not. halfstep .and. master) then
            write(msg,'(a,f8.3,a)')"[multigrid_diffusion:multigrid_solve_diff] Timestep limited somewhere else: dt = ",dt/diff_dt_crs_orig, " of explicit dt_crs."
            call printinfo(msg, stdout)
         endif

         ! set diffBC
         call init_b

         do cr_id = 1, flind%crs%all
            call init_source(cr_id)
            if (vstat%norm_rhs /= 0) then
               if (norm_was_zero(cr_id) .and. master) then
                  write(msg,'(a,i2,a)')"[multigrid_diffusion:multigrid_solve_diff] CR-fluid #",cr_id," is now available in measurable quantities."
                  call printinfo(msg)
               endif
               norm_was_zero(cr_id) = .false.

               call init_solution(cr_id)
               ! do substepping
               call vcycle_hg(cr_id)
               ! enddo
            else
               if (.not. norm_was_zero(cr_id) .and. master) then
                  write(msg,'(a,i2,a)')"[multigrid_diffusion:multigrid_solve_diff] Source norm of CR-fluid #",cr_id," == 0., skipping."
                  call warn(msg)
               endif
               norm_was_zero(cr_id) = .true.
            endif
         enddo

      endif
      ts = set_timer("multigrid_diffusion")
      tot_ts = tot_ts + ts

   end subroutine multigrid_solve_diff

!!$ ============================================================================
!>
!! \brief Make a local copy of source
!<

   subroutine init_source(cr_id)

      use multigridvars,      only: roof, source, defect, correction
      use initcosmicrays,     only: iarr_crs
      use grid,               only: all_cg
      use gc_list,            only: cg_list_element
      use grid_cont,          only: grid_container
      use multigridbasefuncs, only: norm_sq
      use dataio_pub,         only: die
      use multigridhelpers,   only: set_dirty, check_dirty

      implicit none

      integer, intent(in) :: cr_id !< CR component index
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg

      if (all_cg%cnt > 1) call die("[multigrid_diffusion:init_source] multiple grid pieces per procesor not implemented yet") !nontrivial plvl

      call set_dirty(source)
      call set_dirty(correction)
      call set_dirty(defect)
      ! Trick residual subroutine to initialize with: u + (1-theta) dt grad (c grad u)
      if (diff_theta /= 0.) then
         cgl => all_cg%first
         do while (associated(cgl))
            cg => cgl%cg
            roof%mgvar(roof%is:roof%ie, roof%js:roof%je, roof%ks:roof%ke, correction) = (1. -1./diff_theta) * cg%u%arr(iarr_crs(cr_id), cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke)
            roof%mgvar(roof%is:roof%ie, roof%js:roof%je, roof%ks:roof%ke, defect)     =     -1./diff_theta  * cg%u%arr(iarr_crs(cr_id), cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke)
            call residual(roof, defect, correction, source, cr_id)
            cgl => cgl%nxt
         enddo
      else
         call die("[multigrid_diffusion:init_source] diff_theta = 0 not supported.")
      endif
      call check_dirty(roof, source, "init source")

      call norm_sq(source, vstat%norm_rhs)

   end subroutine init_source

!!$ ============================================================================
!>
!! \brief Initialize solution with current CR density (no solution recycling as yet)
!<

   subroutine init_solution(cr_id)

      use dataio_pub,       only: die
      use grid,             only: all_cg
      use gc_list,          only: cg_list_element
      use grid_cont,        only: grid_container
      use initcosmicrays,   only: iarr_crs
      use multigridhelpers, only: set_dirty, check_dirty
      use multigridvars,    only: roof, solution

      implicit none

      integer, intent(in) :: cr_id !< CR component index
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg

      if (all_cg%cnt > 1) call die("[multigrid_diffusion:init_solution] multiple grid pieces per procesor not implemented yet") !nontrivial plvl

      call set_dirty(solution)
      cgl => all_cg%first
      do while (associated(cgl))
         cg => cgl%cg
         roof%mgvar(roof%is:roof%ie, roof%js:roof%je, roof%ks:roof%ke, solution) = cg%u%arr(iarr_crs(cr_id),  cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke)
         cgl => cgl%nxt
      enddo
      call check_dirty(roof, solution, "init solution")

   end subroutine init_solution

!!$ ============================================================================
!>
!! \brief Initialize magnetic field components
!!
!! \todo test what happens if we use magnetic field interpolated to the cell centers (some optimizations may then become available)
!<

   subroutine init_b

      use constants,          only: I_ONE
      use dataio_pub,         only: die
      use domain,             only: D_x, D_y, D_z
      use fluidindex,         only: ibx, iby, ibz
      use grid,               only: all_cg
      use gc_list,            only: cg_list_element
      use grid_cont,          only: grid_container
      use multigridbasefuncs, only: restrict_all
      use multigridhelpers,   only: set_dirty, check_dirty, dirty_label
      use multigridmpifuncs,  only: mpi_multigrid_bnd
      use multigridvars,      only: base, roof, extbnd_mirror, plvl

      implicit none

      integer(kind=4) :: ib
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg
      type(plvl), pointer :: curl

      if (all_cg%cnt > 1) call die("[multigrid_diffusion:init_b] multiple grid pieces per procesor not implemented yet") !nontrivial plvl

      if (diff_bx+iby-ibx /= diff_by .or. diff_bx+ibz-ibx /= diff_bz) call die("[multigrid_diffusion:init_b] Something is wrong with diff_by or diff_bz indices.")

      do ib = ibx, ibz
         call set_dirty(diff_bx+ib-ibx)
         cgl => all_cg%first
         do while (associated(cgl))
            cg => cgl%cg
            roof%mgvar(       roof%is-D_x:roof%ie+D_x, roof%js-D_y:roof%je+D_y, roof%ks-D_z:roof%ke+D_z, diff_bx+ib-ibx) = &
                 cg%b%arr(ib,   cg%is-D_x:  cg%ie+D_x,   cg%js-D_y:  cg%je+D_y,   cg%ks-D_z:  cg%ke+D_z)
            cgl => cgl%nxt
         enddo
         call restrict_all(diff_bx+ib-ibx)             ! Implement correct restriction (and probably also separate inter-process communication) routines

         curl => base
         do while (associated(curl%finer)) ! from base to one level below roof
            call mpi_multigrid_bnd(curl, diff_bx+ib-ibx, I_ONE, extbnd_mirror, .true.) !> \todo use global boundary type for B
            !>
            !! |deprecated BEWARE b is set on a staggered grid; corners should be properly set here (now they are not)
            !! the problem is that the cg%b%arr(:,:,:,:) elements are face-centered so restriction and external boundaries should take this into account
            !<
            write(dirty_label, '(a,i1)')"init b",ib
            call check_dirty(curl, diff_bx+ib-ibx, dirty_label)
            curl => curl%finer
         enddo
      enddo

   end subroutine init_b

!!$ ============================================================================
!>
!! \brief Huang-Greengard V-cycle
!<

   subroutine vcycle_hg(cr_id)

      use dataio_pub,         only: msg, warn
      use grid,               only: all_cg
      use gc_list,            only: cg_list_element
      use grid_cont,          only: grid_container
      use initcosmicrays,     only: iarr_crs
      use mpisetup,           only: master
      use multigridbasefuncs, only: norm_sq, restrict_all, prolong_level
      use multigridhelpers,   only: set_dirty, check_dirty, do_ascii_dump, numbered_ascii_dump, brief_v_log, dirty_label
      use multigridvars,      only: source, defect, solution, correction, base, roof, plvl, ts, tot_ts
      use timer,              only: set_timer

      implicit none

      integer, intent(in) :: cr_id !< CR component index

      real, parameter    :: barely_greater_than_1 = 1.05
      integer, parameter :: convergence_history = 2
      integer            :: v
      real               :: norm_lhs, norm_rhs, norm_old
      logical            :: dump_every_step
      type(plvl), pointer :: curl
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg

      write(vstat%cprefix,'("C",i1,"-")') cr_id !> \deprecated BEWARE: this is another place with 0 <= cr_id <= 9 limit
      write(dirty_label, '("md_",i1,"_dump")')  cr_id

      inquire(file = "_dump_every_step_", EXIST=dump_every_step) ! use for debug only
      do_ascii_dump = do_ascii_dump .or. dump_every_step

      call norm_sq(solution, norm_rhs)
      norm_old = norm_rhs

      do v = 0, max_cycles

         call set_dirty(defect)

         call residual(roof, source, solution, defect, cr_id)
         call norm_sq(defect, norm_lhs)
         ts = set_timer("multigrid_diffusion")
         tot_ts = tot_ts + ts

         vstat%count = v
         if (norm_lhs /= 0) then
            vstat%factor(vstat%count) = norm_old/norm_lhs
         else
            vstat%factor(vstat%count) = huge(1.0)
         endif
         vstat%time(vstat%count) = ts

         norm_old = norm_lhs

         if (dump_every_step) call numbered_ascii_dump(dirty_label, v)

         if (norm_lhs/norm_rhs <= norm_tol) exit

         if (v>convergence_history) then
            if (product(vstat%factor(v-convergence_history:v)) < barely_greater_than_1 .and. master) then
               write(msg, '(a,i3,a,g15.5)')"[multigrid_diffusion:vcycle_hg] Too slow convergence: cycle = ",v,", norm_lhs/norm_rhs = ", norm_lhs/norm_rhs
               call warn(msg)
               exit
            endif
         endif

         call restrict_all(defect)

         call set_dirty(correction)
         base%mgvar(:, :, :, correction) = 0.

         curl => base
         do while (associated(curl))
            call approximate_solution(curl, defect, correction, cr_id)
            call prolong_level(curl, correction)
            curl => curl%finer
         enddo

         call check_dirty(roof, correction, "c_residual")
         call check_dirty(roof, defect, "d_residual")
         roof%mgvar     (roof%is:roof%ie, roof%js:roof%je, roof%ks:roof%ke, solution) = &
              roof%mgvar(roof%is:roof%ie, roof%js:roof%je, roof%ks:roof%ke, solution) - &
              roof%mgvar(roof%is:roof%ie, roof%js:roof%je, roof%ks:roof%ke, correction)

      enddo

      if (dump_every_step) call numbered_ascii_dump(dirty_label)

      call check_dirty(roof, solution, "v_soln")

      if (v > max_cycles) then
         if (master .and. norm_lhs/norm_rhs > norm_tol) then
            write(msg, '(a,i3,a,g15.5)')"[multigrid_diffusion:vcycle_hg] Not enough V-cycles to achieve convergence: cycle = ",v,", norm_lhs/norm_rhs = ", norm_lhs/norm_rhs
            call warn(msg)
         endif
         v = max_cycles
      endif

      vstat%norm_final = norm_lhs/norm_rhs
      call brief_v_log(vstat)

      call norm_sq(solution, norm_rhs)
      call norm_sq(defect, norm_lhs)
!     Do we need to take care of boundaries here?
!      call mpi_multigrid_bnd(roof, solution, I_ONE, diff_extbnd)
!      cg%u%arr(iarr_crs(cr_id), is-D_x:cg%ie+D_x, cg%js-D_y:cg%je+D_y, cg%ks-D_z:cg%ke+D_z) = roof%mgvar(roof%is-D_x:roof%ie+D_x, roof%js-D_y:roof%je+D_y, roof%ks-D_z:roof%ke+D_z, solution)
      cgl => all_cg%first
      do while (associated(cgl))
         cg => cgl%cg
         cg%u%arr(iarr_crs(cr_id), cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) = roof%mgvar(roof%is:roof%ie, roof%js:roof%je, roof%ks:roof%ke, solution)
         cgl => cgl%nxt
      enddo

   end subroutine vcycle_hg

!!$ ============================================================================
!>
!! \brief Compute diffusive flux in the crdim-direction
!!
!! OPT: this routine can consume as much as half of the CPU time of the whole simulation
!! OPT: moving the loops over i, j and k inside it can speed them up by more than 20% (or much more after merging directional variants)
!! OPT: bcomp and magb are invariants of the solution and can be precomputed in init_b (requires additional 4*ndim temporary arrays)
!! OPT: *curl%idl(idir) is a measurably faster than /curl%dl(idir) but will require update of the "gold" template due to roundoff differences
!<

   subroutine diff_flux(crdim, im, soln, curl, cr_id, Keff)

      use constants,      only: xdim, ydim, zdim, ndims
      use dataio_pub,     only: die
      use domain,         only: has_dir
      use grid,           only: all_cg
      use grid_cont,      only: grid_container
      use initcosmicrays, only: K_crs_perp, K_crs_paral
      use multigridvars,  only: plvl

      implicit none

      integer(kind=4),           intent(in)  :: crdim
      integer, dimension(ndims), intent(in)  :: im           !< [first cell index, second cell index, third cell index]
      integer(kind=4),           intent(in)  :: soln         !< multigrid variable to differentiate
      type(plvl), pointer,       intent(in)  :: curl         !< level on which differentiate
      integer,                   intent(in)  :: cr_id        !< CR component index
      real, optional,            intent(out) :: Keff         !< effective diffusion coefficient for relaxation

      real                                   :: magb, fcrdif, kbm
      real                                   :: b_par, b_perp, d_par, db_perp
      integer, dimension(ndims)              :: ilm, imp, imm, ilmp, ilmm
      integer(kind=4)                        :: idir
      logical, dimension(ndims)              :: present_not_crdim
      type(grid_container), pointer          :: cg

      ilm(:) = im(:) ; ilm(crdim) = ilm(crdim) - 1
      present_not_crdim(:) = has_dir(:) .and. ( [ xdim,ydim,zdim ] /= crdim )

      cg => all_cg%first%cg
      if (all_cg%cnt > 1) call die("[multigrid_diffusion:diff_flux] multiple grid pieces per procesor not implemented yet") !nontrivial plvl

      ! Assumes has_dir(crdim)
      !> \warning *curl%idl(crdim) makes a difference
      d_par = (curl%mgvar(im(xdim), im(ydim), im(zdim), soln) - curl%mgvar(ilm(xdim), ilm(ydim), ilm(zdim), soln)) / curl%dl(crdim)
      fcrdif = K_crs_perp(cr_id) * d_par
      if (present(Keff)) Keff = K_crs_perp(cr_id)

      if (K_crs_paral(cr_id) /= 0.) then

         b_perp = 0.
         db_perp = 0.
         b_par = curl%mgvar(im(xdim), im(ydim), im(zdim), idiffb(crdim))
         magb = b_par**2

         do idir = xdim, zdim
            if (present_not_crdim(idir)) then
               imp(:) = im(:) ; imp(idir) = imp(idir) + 1 ; ilmp(:) = imp(:) ; ilmp(crdim) = ilmp(crdim) - 1
               imm(:) = im(:) ; imm(idir) = imm(idir) - 1 ; ilmm(:) = imm(:) ; ilmm(crdim) = ilmm(crdim) - 1
               b_perp = sum(curl%mgvar(ilm(xdim):imp(xdim), ilm(ydim):imp(ydim), ilm(zdim):imp(zdim), idiffb(idir)))*0.25
               magb = magb + b_perp**2
               !> \warning *curl%idl(crdim) makes a difference
               db_perp = db_perp + b_perp*((curl%mgvar(ilmp(xdim), ilmp(ydim), ilmp(zdim), soln) + curl%mgvar(imp(xdim), imp(ydim), imp(zdim), soln)) - &
                  &                        (curl%mgvar(ilmm(xdim), ilmm(ydim), ilmm(zdim), soln) + curl%mgvar(imm(xdim), imm(ydim), imm(zdim), soln))) * 0.25 / curl%dl(idir)
            endif
         enddo

         if (magb /= 0.) then
            kbm = K_crs_paral(cr_id) * b_par / magb
            fcrdif = fcrdif + kbm * (b_par*d_par + db_perp)
            if (present(Keff)) Keff = Keff + kbm * b_par
         endif

      endif

      cg%wa(im(xdim), im(ydim), im(zdim)) = fcrdif ! * diff_theta * dt / curl%dl(crdim) !> \warning *curl%idl(crdim) makes a difference

   end subroutine diff_flux

!!$ ============================================================================
!>
!! \brief 2nd order: grad (c grad)
!!
!! defect = solution - source - grad (c grad (solution))
!<

   subroutine residual(curl, src, soln, def, cr_id)

      use constants,         only: xdim, ydim, zdim, I_ONE, ndims, LO, HI
      use dataio_pub,        only: die
      use domain,            only: has_dir
      use global,            only: dt
      use grid,              only: all_cg
      use grid_cont,         only: grid_container
      use multigridhelpers,  only: check_dirty
      use multigridmpifuncs, only: mpi_multigrid_bnd
      use multigridvars,     only: plvl

      implicit none

      type(plvl), pointer, intent(in) :: curl !< level for which approximate the solution
      integer(kind=4),     intent(in) :: src   !< index of source in lvl()%mgvar
      integer(kind=4),     intent(in) :: soln  !< index of solution in lvl()%mgvar
      integer(kind=4),     intent(in) :: def   !< index of defect in lvl()%mgvar
      integer,             intent(in) :: cr_id !< CR component index

      integer                         :: i, j, k
      integer(kind=4)                 :: idir
      integer, dimension(ndims)       :: iml, imh
      type(grid_container), pointer   :: cg

      cg => all_cg%first%cg
      if (all_cg%cnt > 1) call die("[multigrid_diffusion:residual] multiple grid pieces per procesor not implemented yet") !nontrivial plvl

      call mpi_multigrid_bnd(curl, soln, I_ONE, diff_extbnd, .true.) ! corners are required for fluxes

      do k = curl%ks, curl%ke
         curl%mgvar     (curl%is:curl%ie, curl%js:curl%je, k, def)  =   &
              curl%mgvar(curl%is:curl%ie, curl%js:curl%je, k, soln) -   &
              curl%mgvar(curl%is:curl%ie, curl%js:curl%je, k, src)
      enddo

      do idir = xdim, zdim
         if (has_dir(idir)) then
            imh = curl%ijkse(:,HI) ; imh(idir) = imh(idir) + 1
            iml = curl%ijkse(:,LO) ; iml(idir) = iml(idir) + 1
            do k = curl%ks, imh(zdim)
               do j = curl%js, imh(ydim)
                  do i = curl%is, imh(xdim)
                     call diff_flux(idir, [i, j, k], soln, curl, cr_id)
                  enddo
               enddo
            enddo
            curl%mgvar     (curl%is  :curl%ie,   curl%js:curl%je, curl%ks:curl%ke, def)    = &
                 curl%mgvar(curl%is  :curl%ie,   curl%js:curl%je, curl%ks:curl%ke, def)    - &
              (       cg%wa(iml(xdim):imh(xdim), iml(ydim):imh(ydim), iml(zdim):imh(zdim)) - &
                      cg%wa(curl%is  :curl%ie,   curl%js:curl%je, curl%ks:curl%ke) ) * diff_theta * dt / curl%dl(idir)
         endif
      enddo

      call check_dirty(curl, def, "res def")

   end subroutine residual

!!$ ============================================================================
!>
!! \brief Relaxation.
!!
!! \details This is the most costly routine in a serial run. Try to find optimal values for nsmool.
!! It seems that for this particular scheme nsmoob can be set even to 1 when the base level is coarsened enough.
!! This routine also depends a lot on communication so it  may limit scalability of the multigrid.
!<

   subroutine approximate_solution(curl, src, soln, cr_id)

      use constants,         only: xdim, ydim, zdim, one, half, I_ONE, ndims
      use dataio_pub,        only: die
      use domain,            only: has_dir
      use global,            only: dt
      use grid,              only: all_cg
      use grid_cont,         only: grid_container
      use multigridmpifuncs, only: mpi_multigrid_bnd
      use multigridvars,     only: plvl, base, extbnd_donothing

      implicit none

      type(plvl), pointer, intent(in) :: curl  !< level for which approximate the solution
      integer(kind=4),     intent(in) :: src   !< index of source in lvl()%mgvar
      integer(kind=4),     intent(in) :: soln  !< index of solution in lvl()%mgvar
      integer,             intent(in) :: cr_id !< CR component index

      integer, parameter              :: RED_BLACK = 2 !< the checkerboard requires two sweeps

      integer                         :: n, i, j, k, i1, j1, k1, id, jd, kd, nsmoo
      integer(kind=4)                 :: idir
      integer, dimension(ndims)       :: im, ih
      real,    dimension(ndims)       :: idl2
      real                            :: Keff1, Keff2, dLdu, temp
      type(grid_container), pointer   :: cg

      idl2 = [curl%idx2, curl%idy2, curl%idz2]
      cg => all_cg%first%cg
      if (all_cg%cnt > 1) call die("[multigrid_diffusion:approximate_solution] multiple grid pieces per procesor not implemented yet") !nontrivial plvl

      if (associated(curl, base)) then
         nsmoo = nsmoob
      else
         nsmoo = nsmool
      endif

      i1 = curl%is; id = 1 ! mv to multigridvars, init_multigrid
      j1 = curl%js; jd = 1
      k1 = curl%ks; kd = 1
      if (has_dir(xdim)) then
         id = RED_BLACK
      else if (has_dir(ydim)) then
         jd = RED_BLACK
      else if (has_dir(zdim)) then
         kd = RED_BLACK
      endif

      do n = 1, RED_BLACK*nsmoo
         if (mod(n,2) == 1) then
            call mpi_multigrid_bnd(curl, soln, I_ONE, diff_extbnd, .true.) ! corners are required for fluxes
         else
            call mpi_multigrid_bnd(curl, soln, I_ONE, extbnd_donothing, .true.)
         endif

         if (kd == RED_BLACK) k1 = curl%ks + mod(n, RED_BLACK)
         do k = k1, curl%ke, kd
            if (jd == RED_BLACK) j1 = curl%js + mod(n+k, RED_BLACK)
            do j = j1, curl%je, jd
               if (id == RED_BLACK) i1 = curl%is + mod(n+j+k, RED_BLACK)
               do i = i1, curl%ie, id

                  temp = curl%mgvar(i, j, k, soln) - curl%mgvar(i, j, k, src)
                  dLdu = 0.

                  im = [i, j, k]

                  do idir = xdim, zdim
                     if (has_dir(idir)) then

                        ih = im ; ih(idir) = ih(idir) + 1
                        call diff_flux(idir, im, soln, curl, cr_id, Keff1)
                        call diff_flux(idir, ih, soln, curl, cr_id, Keff2)

                        temp = temp - (cg%wa(ih(xdim), ih(ydim), ih(zdim)) - cg%wa(i, j, k)) * diff_theta * dt / curl%dl(idir)
                        dLdu = dLdu - 2 * (Keff1 + Keff2) * idl2(idir)

                     endif
                  enddo

                  ! ToDo add an option to automagically fine-tune overrelax
                  curl%mgvar(i, j, k, soln) = curl%mgvar(i, j, k, soln) - overrelax * temp/(one - half * diff_theta * dt * dLdu)

               enddo
            enddo
         enddo
      enddo

   end subroutine approximate_solution

end module multigrid_diffusion
