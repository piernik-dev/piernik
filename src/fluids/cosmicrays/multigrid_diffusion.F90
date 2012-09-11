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

   use constants,        only: cbuff_len, ndims, dsetnamelen
   use multigrid_vstats, only: vcycle_stats

   implicit none

   private
   public :: init_multigrid_diff, init_multigrid_diff_post, cleanup_multigrid_diff, multigrid_solve_diff
   public :: diff_tstep_fac, diff_explicit, diff_dt_crs_orig

   ! namelist parameters
   real               :: norm_tol                                     !< stop V-cycle iterations when the ratio of norms ||residual||/||source|| is below this value
   real               :: vcycle_abort                                 !< abort the V-cycle when lhs norm raises by this factor
   real               :: overrelax                                    !< overrealaxation factor (if < 1. then works as underrelaxation), use with care
   integer(kind=4)    :: max_cycles                                   !< Maximum allowed number of V-cycles
   integer(kind=4)    :: nsmool                                       !< smoothing cycles per call
   integer(kind=4)    :: nsmoob                                       !< smoothing cycles on coarsest level;  \todo implement a convergence check
   real               :: diff_theta                                   !< 0. is explicit, 1. is fully implicit 0.5 is Crank-Nicholson
   real,    protected :: diff_tstep_fac                               !< How much we stretch timestep. Note that for diff_theta == 0. this should not be > 1.
   logical, protected :: diff_explicit                                !< If .true. then do not use multigrid for diffusion
   logical            :: allow_explicit                               !< When timestep is limited somewhere else, allow explicit calculation (should be a bit faster)
   real               :: diff_dt_crs_orig                             !< timestep calculated at timestepcosmicrays.F90, before enlarging by diff_tstep_fac
   character(len=cbuff_len) :: diff_bnd_str                           !< Type of diffusion boundary conditions. Can be "isolated", "reflecting" or "zero" (there are some aliases as well)

   ! mgvar entries for the B field
   character(len=dsetnamelen), parameter :: diff_bx_n = "diff_bx"  !< index of B_x in the cg%b(:,:,:,:) array
   character(len=dsetnamelen), parameter :: diff_by_n = "diff_by"  !< index of B_y in the cg%b(:,:,:,:) array
   character(len=dsetnamelen), parameter :: diff_bz_n = "diff_bz"  !< index of B_z in the cg%b(:,:,:,:) array
   integer, dimension(ndims) :: idiffb

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
!! The list is active while \b "COSM_RAYS" and \b "MULTIGRID" are defined.
!! \n \n
!<
   subroutine init_multigrid_diff

      use cg_list_global, only: all_cg
      use constants,      only: BND_ZERO, BND_XTRAP, BND_REF, BND_NEGREF, xdim, ydim, zdim, GEO_XYZ, half, zero, one
      use dataio_pub,     only: par_file, ierrh, namelist_errh, compare_namelist, cmdl_nml, lun      ! QA_WARN required for diff_nml
      use dataio_pub,     only: die, warn, msg
      use domain,         only: dom
      use mpisetup,       only: master, slave, nproc, ibuff, rbuff, lbuff, cbuff, piernik_MPI_Bcast
      use multigridvars,  only: single_base
      use named_array_list, only: qna

      implicit none

      logical, save :: frun = .true.          !< First run flag

      namelist /MULTIGRID_DIFFUSION/ norm_tol, vcycle_abort, max_cycles, nsmool, nsmoob, overrelax, &
           &                         diff_theta, diff_tstep_fac, diff_explicit, allow_explicit, diff_bnd_str

      if (.not.frun) call die("[multigrid_diffusion:init_multigrid_diff] Called more than once.")
      frun = .false.
      if (dom%geometry_type /= GEO_XYZ) call die("[multigrid_gravity:init_multigrid_gravdiffusion:init_multigrid_diff] non-cartesian geometry not implemented yet.")

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

      call piernik_MPI_Bcast(cbuff, cbuff_len)
      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(rbuff)
      call piernik_MPI_Bcast(lbuff)

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

      ! boundaries
      diff_extbnd = BND_ZERO
      select case (diff_bnd_str)
         case ("isolated", "iso", "free")
            diff_extbnd = BND_XTRAP
         case ("reflecting", "refl", "styrofoam")
            diff_extbnd = BND_REF
         case ("zero", "cold", "antireflecting")
            diff_extbnd = BND_NEGREF
         case default
            write(msg,'(3a)')"[multigrid_diffusion:init_multigrid_diff] Non-recognized boundary description '",diff_bnd_str,"'"
            call die(msg)
      end select

      !diffusion
      if (.not. diff_explicit) then
         if (diff_theta < zero .or. diff_theta > one) call die("[multigrid_diffusion:init_multigrid] diff_theta must lie in the range [0. .. 1.]")
         if (diff_theta < half .and. diff_tstep_fac > one .and. master) call warn("[multigrid_diffusion:init_multigrid] diff_tstep_fac > 1. for diff_theta < 0.5 might be unstable")
         ! calculate exact limit formula
         ! for diff_theta = 0. stable diff_tstep_fac is 0.5 in 2D (guess: 0.333 in 3D)
         ! for diff_theta < 0.5 stable diff_tstep_fac rises by 1./(1.-2.*diff_theta)
      endif

      if (overrelax /= 1 .and. master) then
         write(msg, '(a,f8.5)')"[multigrid_diffusion:init_multigrid_diff] Overrelaxation factor = ", overrelax
         call warn(msg)
      endif

      if (single_base .and. nproc > 1) then
         call warn("[multigrid_diffusion:init_multigrid_diff] single_base disabled just in case")
         single_base = (nproc == 1)
      endif

      call all_cg%reg_var(diff_bx_n, multigrid = .true.)
      call all_cg%reg_var(diff_by_n, multigrid = .true.)
      call all_cg%reg_var(diff_bz_n, multigrid = .true.)
      idiffb(xdim) = qna%ind(diff_bx_n)
      idiffb(ydim) = qna%ind(diff_by_n)
      idiffb(zdim) = qna%ind(diff_bz_n)

   end subroutine init_multigrid_diff

!!$ ============================================================================
!>
!! \brief Initialization - continued after allocation of everything interesting
!<

   subroutine init_multigrid_diff_post

      use dataio_pub,         only: die
      use fluidindex,         only: flind

      implicit none

      if (allocated(norm_was_zero)) call die("[multigrid_diffusion:init_multigrid] norm_was_zero already allocated")
      allocate(norm_was_zero(flind%crs%all))
      norm_was_zero(:) = .false.

      call vstat%init(max_cycles)

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
      use fluidindex,    only: flind
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
            call cr_diff(zdim)
            call cr_diff(ydim)
            call cr_diff(xdim)
         else
            call cr_diff(xdim)
            call cr_diff(ydim)
            call cr_diff(zdim)
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

      use cg_list_global, only: all_cg
#if defined(__INTEL_COMPILER)
      use cg_level_connected,  only: cg_level_connected_T  ! QA_WARN workaround for stupid INTEL compiler
#endif /* __INTEL_COMPILER */
      use cg_level_connected,  only: finest
      use dataio_pub,     only: die
      use cg_list,        only: dirty_label
      use cg_list_dataop, only: ind_val
      use cg_leaves,      only: leaves
      use initcosmicrays, only: iarr_crs
      use multigridvars,  only: source, defect, correction
      use named_array_list, only: qna, wna

      implicit none

      integer, intent(in) :: cr_id !< CR component index

      call all_cg%set_dirty(source)
      call all_cg%set_dirty(correction)
      call all_cg%set_dirty(defect)
      ! Trick residual subroutine to initialize with: u + (1-theta) dt grad (c grad u)
      if (diff_theta /= 0.) then
         call leaves%wq_copy(wna%fi, int(iarr_crs(cr_id)), qna%wai)
         call leaves%q_lin_comb( [ ind_val(qna%wai, (1. -1./diff_theta)) ], correction)
         call leaves%q_lin_comb( [ ind_val(qna%wai,     -1./diff_theta ) ], defect)
         call residual(finest, defect, correction, source, cr_id)
      else
         call die("[multigrid_diffusion:init_source] diff_theta = 0 not supported.")
      endif
      write(dirty_label, '(a,i2)')"init source#", cr_id
      call finest%check_dirty(source, dirty_label)

      vstat%norm_rhs = leaves%norm_sq(source)

   end subroutine init_source

!!$ ============================================================================
!>
!! \brief Initialize solution with current CR density (no solution recycling as yet)
!<

   subroutine init_solution(cr_id)

      use cg_list_global, only: all_cg
#if defined(__INTEL_COMPILER)
      use cg_level_connected,  only: cg_level_connected_T  ! QA_WARN workaround for stupid INTEL compiler
#endif /* __INTEL_COMPILER */
      use cg_level_connected,  only: finest
      use cg_leaves,      only: leaves
      use initcosmicrays, only: iarr_crs
      use multigridvars,  only: solution
      use named_array_list, only: wna

      implicit none

      integer, intent(in) :: cr_id !< CR component index

      call all_cg%set_dirty(solution)
      call leaves%wq_copy(wna%fi, int(iarr_crs(cr_id)), solution)
      call finest%check_dirty(solution, "init solution")

   end subroutine init_solution

!!$ ============================================================================
!>
!! \brief Initialize magnetic field components
!!
!! \todo test what happens if we use magnetic field interpolated to the cell centers (some optimizations may then become available)
!<

   subroutine init_b

      use cg_list_global, only: all_cg
      use cg_level_connected,  only: cg_level_connected_T, coarsest, finest
      use constants,      only: I_ONE, xdim, zdim, HI, LO, BND_REF
      use domain,         only: dom
      use cg_leaves,      only: leaves
      use cg_list,        only: cg_list_element, dirty_label
      use grid_cont,      only: grid_container
      use named_array,    only: p3, p4
      use named_array_list, only: wna

      implicit none

      integer(kind=4) :: ib
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg
      type(cg_level_connected_T),   pointer :: curl

      do ib = xdim, zdim
         call all_cg%set_dirty(idiffb(ib))
#if 1
         cgl => leaves%first
         do while (associated(cgl))
            cg => cgl%cg
            p3 => cg%q(idiffb(ib))%span(cg%ijkse(:,LO)-dom%D_(:),cg%ijkse(:,HI)+dom%D_(:))
            p4 => cg%w(wna%bi )%span(cg%ijkse(:,LO)-dom%D_(:),cg%ijkse(:,HI)+dom%D_(:))
            p3 = p4(ib,:,:,:)
            cgl => cgl%nxt
         enddo
#else
         ! This works well but copies all guardcells, which is not necessary
         call leaves%wq_copy(wna%bi, ib, idiffb(ib))
#endif
         call finest%restrict_to_floor_q_1var(idiffb(ib))             ! Implement correct restriction (and probably also separate inter-process communication) routines

         curl => coarsest
         do while (associated(curl%finer)) ! from coarsest to one level below finest
            call curl%arr3d_boundaries(idiffb(ib), nb = I_ONE, bnd_type = BND_REF, corners = .true.) !> \todo use global boundary type for B
            !>
            !! |deprecated BEWARE b is set on a staggered grid; corners should be properly set here (now they are not)
            !! the problem is that the cg%b(:,:,:,:) elements are face-centered so restriction and external boundaries should take this into account
            !<
            curl => curl%finer
         enddo
         write(dirty_label, '(a,i1)')"init b",ib
         call curl%check_dirty(idiffb(ib), dirty_label)
      enddo

   end subroutine init_b

!!$ ============================================================================
!>
!! \brief Huang-Greengard V-cycle
!<

   subroutine vcycle_hg(cr_id)

      use cg_list_global, only: all_cg
      use cg_level_connected,  only: cg_level_connected_T, coarsest, finest
      use dataio_pub,     only: msg, warn
      use cg_list,        only: dirty_label
      use cg_list_dataop, only: ind_val
      use global,         only: do_ascii_dump
      use cg_leaves,      only: leaves
      use initcosmicrays, only: iarr_crs
      use mpisetup,       only: master
      use multigridvars,  only: source, defect, solution, correction, ts, tot_ts
      use named_array_list, only: wna
      use timer,          only: set_timer

      implicit none

      integer, intent(in) :: cr_id !< CR component index

      real, parameter    :: barely_greater_than_1 = 1.05
      integer, parameter :: convergence_history = 2
      integer            :: v
      real               :: norm_lhs, norm_rhs, norm_old
      logical            :: dump_every_step
      type(cg_level_connected_T), pointer :: curl

      write(vstat%cprefix,'("C",i1,"-")') cr_id !> \deprecated BEWARE: this is another place with 0 <= cr_id <= 9 limit
      write(dirty_label, '("md_",i1,"_dump")')  cr_id

      inquire(file = "_dump_every_step_", EXIST=dump_every_step) ! use for debug only
      do_ascii_dump = do_ascii_dump .or. dump_every_step

      norm_lhs = 0.
      norm_rhs = leaves%norm_sq(solution)
      norm_old = norm_rhs

      do v = 0, max_cycles

         call all_cg%set_dirty(defect)

         call residual(finest, source, solution, defect, cr_id)
         norm_lhs = leaves%norm_sq(defect)
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

         if (dump_every_step) call all_cg%numbered_ascii_dump([ source, solution, defect, correction ], dirty_label, v)

         if (norm_lhs/norm_rhs <= norm_tol) exit

         if (v>convergence_history) then
            if (product(vstat%factor(v-convergence_history:v)) < barely_greater_than_1) then
               if (master) then
                  write(msg, '(a,i3,a,g15.5)')"[multigrid_diffusion:vcycle_hg] Too slow convergence: cycle = ",v,", norm_lhs/norm_rhs = ", norm_lhs/norm_rhs
                  call warn(msg)
               endif
               exit
            endif
         endif

         call finest%restrict_to_floor_q_1var(defect)

         !call all_cg%set_dirty(correction)
         call coarsest%set_q_value(correction, 0.)

         curl => coarsest
         do while (associated(curl))
            call approximate_solution(curl, defect, correction, cr_id)
            if (.not. associated(curl, finest)) call curl%prolong_q_1var(correction)
            curl => curl%finer
         enddo

         call finest%check_dirty(correction, "c_residual")
         call finest%check_dirty(defect, "d_residual")
         call leaves%q_lin_comb( [ ind_val(solution, 1.), ind_val(correction, -1.) ], solution) ! solution := solution - correction

      enddo

      if (dump_every_step) call all_cg%numbered_ascii_dump([ source, solution, defect, correction ], dirty_label)

      call finest%check_dirty(solution, "v_soln")

      if (v > max_cycles) then
         if (master .and. norm_lhs/norm_rhs > norm_tol) then
            write(msg, '(a,i3,a,g15.5)')"[multigrid_diffusion:vcycle_hg] Not enough V-cycles to achieve convergence: cycle = ",v,", norm_lhs/norm_rhs = ", norm_lhs/norm_rhs
            call warn(msg)
         endif
         v = max_cycles
      endif

      vstat%norm_final = norm_lhs/norm_rhs
      call vstat%brief_v_log

      norm_rhs = leaves%norm_sq(solution)
      norm_lhs = leaves%norm_sq(defect)
!     Do we need to take care of boundaries here?
!      call finest%arr3d_boundaries(solution, nb = I_ONE, bnd_type = diff_extbnd)
!      cg%u%span(iarr_crs(cr_id),cg%ijkse(:,LO)-dom%D_,cg%ijkse(:,HI)+dom%D_) = cg%q(solution)%span(cg%ijkse(:,LO)-dom%D_,cg%ijkse(:,HI)+dom%D_)

      call leaves%qw_copy(solution, wna%fi, int(iarr_crs(cr_id)))

   end subroutine vcycle_hg

!!$ ============================================================================
!>
!! \brief Compute diffusive flux in the crdim-direction
!!
!! OPT: this routine can consume as much as half of the CPU time of the whole simulation
!! OPT: moving the loops over i, j and k inside it can speed them up by more than 20% (or much more after merging directional variants)
!! OPT: b_par/perp and magb are invariants of the solution and can be precomputed in init_b (requires additional 4*ndim temporary arrays)
!<

   subroutine diff_flux(crdim, im, soln, cg, cr_id, Keff)

      use constants,      only: xdim, ydim, zdim, ndims, oneq
      use domain,         only: dom
      use grid_cont,      only: grid_container
      use initcosmicrays, only: K_crs_perp, K_crs_paral

      implicit none

      integer(kind=4),               intent(in)  :: crdim        !< direction in which we calculate flux
      integer, dimension(:),         intent(in)  :: im           !< [first cell index, second cell index, third cell index]
      integer,                       intent(in)  :: soln         !< multigrid variable to differentiate
      type(grid_container), pointer, intent(in)  :: cg           !< level on which differentiate
      integer,                       intent(in)  :: cr_id        !< CR component index
      real, optional,                intent(out) :: Keff         !< effective diffusion coefficient for relaxation

      real                                   :: magb, fcrdif, kbm
      real                                   :: b_par, b_perp, d_par, db
      integer, dimension(ndims)              :: ilm, imp, imm, ilmp, ilmm
      integer(kind=4)                        :: idir
      logical, dimension(ndims)              :: present_not_crdim

      ilm(:) = im(:) ; ilm(crdim) = ilm(crdim) - 1
      present_not_crdim(:) = dom%has_dir(:) .and. ( [ xdim,ydim,zdim ] /= crdim )

      ! Assumes dom%has_dir(crdim)
      !> \warning *cg%idl(crdim) makes a difference
      d_par = (cg%q(soln)%point(im) - cg%q(soln)%point(ilm)) * cg%idl(crdim)
      fcrdif = K_crs_perp(cr_id) * d_par
      if (present(Keff)) Keff = K_crs_perp(cr_id)

      if (K_crs_paral(cr_id) /= 0.) then

         b_perp = 0.
         b_par = cg%q(idiffb(crdim))%point(im)
         db = d_par * b_par
         magb = b_par**2

         do idir = xdim, zdim
            if (present_not_crdim(idir)) then
               imp(:) = im(:) ; imp(idir) = imp(idir) + 1 ; ilmp(:) = imp(:) ; ilmp(crdim) = ilmp(crdim) - 1
               imm(:) = im(:) ; imm(idir) = imm(idir) - 1 ; ilmm(:) = imm(:) ; ilmm(crdim) = ilmm(crdim) - 1
               b_perp = sum(cg%q(idiffb(idir))%span(int(ilm, kind=4), int(imp, kind=4)))*oneq
               magb = magb + b_perp**2
               !> \warning *cg%idl(crdim) makes a difference
               db = db + b_perp*((cg%q(soln)%point(ilmp) + cg%q(soln)%point(imp)) - (cg%q(soln)%point(ilmm) + cg%q(soln)%point(imm))) * oneq * cg%idl(idir)
            endif
         enddo

         if (magb /= 0.) then
            kbm = K_crs_paral(cr_id) * b_par / magb
            fcrdif = fcrdif + kbm * db
            if (present(Keff)) Keff = Keff + kbm * b_par
         endif

      endif

      cg%wa(im(xdim), im(ydim), im(zdim)) = fcrdif ! * diff_theta * dt / cg%dl(crdim) !> \warning *cg%idl(crdim) makes a difference

   end subroutine diff_flux

!!$ ============================================================================
!>
!! \brief 2nd order: grad (c grad)
!!
!! defect = solution - source - grad (c grad (solution))
!<

   subroutine residual(curl, src, soln, def, cr_id)

      use cg_level_connected,  only: cg_level_connected_T
      use constants,      only: xdim, ydim, zdim, I_ONE, ndims, LO, HI
      use domain,         only: dom
      use cg_list,        only: cg_list_element
      use cg_list_dataop, only: ind_val
      use global,         only: dt
      use grid_cont,      only: grid_container
      use named_array,    only: p3
      use named_array_list, only: qna

      implicit none

      type(cg_level_connected_T), pointer, intent(in) :: curl !< level for which approximate the solution
      integer,                      intent(in) :: src   !< index of source in cg%q(:)
      integer,                      intent(in) :: soln  !< index of solution in cg%q(:)
      integer,                      intent(in) :: def   !< index of defect in cg%q(:)
      integer,                      intent(in) :: cr_id !< CR component index

      integer                        :: i, j, k
      integer(kind=4)                :: idir
      integer, dimension(ndims)      :: iml, imh
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer  :: cg

      call curl%arr3d_boundaries(soln, nb = I_ONE, bnd_type = diff_extbnd, corners = .true.) ! corners are required for fluxes

      call curl%q_lin_comb([ ind_val(soln, 1.), ind_val(src, -1.) ], def)

      cgl => curl%first
      do while (associated(cgl))
         cg => cgl%cg
         do idir = xdim, zdim
            if (dom%has_dir(idir)) then
               imh = cg%ijkse(:,HI) ; imh(idir) = imh(idir) + 1
               iml = cg%ijkse(:,LO) ; iml(idir) = iml(idir) + 1
               do k = cg%ks, imh(zdim)
                  do j = cg%js, imh(ydim)
                     do i = cg%is, imh(xdim)
                        call diff_flux(idir, [i, j, k], soln, cg, cr_id)
                     enddo
                  enddo
               enddo

               p3 => cg%q(def)%span(cg%ijkse)
               p3 = p3 - (cg%q(qna%wai)%span(int(iml, kind=4), int(imh, kind=4)) - cg%q(qna%wai)%span(cg%ijkse) ) * diff_theta * dt * cg%idl(idir)
            endif
         enddo
         cgl => cgl%nxt
      enddo

!      call curl%check_dirty(def, "res def")

   end subroutine residual

!!$ ============================================================================
!>
!! \brief Relaxation.
!!
!! \details This is the most costly routine in a serial run. Try to find optimal values for nsmool.
!! It seems that for this particular scheme nsmoob can be set even to 1 when the coarsest level is coarsened enough.
!! This routine also depends a lot on communication so it  may limit scalability of the multigrid.
!<

   subroutine approximate_solution(curl, src, soln, cr_id)

      use cg_level_connected,  only: cg_level_connected_T, coarsest
      use constants,      only: xdim, ydim, zdim, one, half, I_ONE, ndims, BND_NONE
      use domain,         only: dom
      use cg_list,        only: cg_list_element
      use global,         only: dt
      use grid_cont,      only: grid_container
      use named_array_list, only: qna

      implicit none

      type(cg_level_connected_T), pointer, intent(in) :: curl  !< level for which approximate the solution
      integer,                      intent(in) :: src   !< index of source in cg%q(:)
      integer,                      intent(in) :: soln  !< index of solution in cg%q(:)
      integer,                      intent(in) :: cr_id !< CR component index

      integer, parameter              :: RED_BLACK = 2 !< the checkerboard requires two sweeps

      integer                         :: n, i, j, k, i1, j1, k1, id, jd, kd, nsmoo
      integer(kind=4)                 :: idir
      integer, dimension(ndims)       :: im, ih
      real                            :: Keff1, Keff2, dLdu, temp
      type(grid_container), pointer   :: cg
      type(cg_list_element), pointer  :: cgl

      if (associated(curl, coarsest)) then
         nsmoo = nsmoob
      else
         nsmoo = nsmool
      endif

      do n = 1, RED_BLACK*nsmoo
         if (mod(n,2) == 1) then
            call curl%arr3d_boundaries(soln, nb = I_ONE, bnd_type = diff_extbnd, corners = .true.) ! corners are required for fluxes
         else
            call curl%arr3d_boundaries(soln, nb = I_ONE, bnd_type = BND_NONE, corners = .true.)
         endif

         cgl => curl%first
         do while (associated(cgl))
            cg => cgl%cg

            i1 = cg%is; id = 1 ! mv to multigridvars, init_multigrid
            j1 = cg%js; jd = 1
            k1 = cg%ks; kd = 1
            if (dom%has_dir(xdim)) then
               id = RED_BLACK
            else if (dom%has_dir(ydim)) then
               jd = RED_BLACK
            else if (dom%has_dir(zdim)) then
               kd = RED_BLACK
            endif

            if (kd == RED_BLACK) k1 = cg%ks + int(mod(n+cg%off(zdim), int(RED_BLACK, kind=8)), kind=4)
            do k = k1, cg%ke, kd
               if (jd == RED_BLACK) j1 = cg%js + int(mod(n+k+sum(cg%off(ydim:zdim)), int(RED_BLACK, kind=8)), kind=4)
               do j = j1, cg%je, jd
                  if (id == RED_BLACK) i1 = cg%is + int(mod(n+j+k+sum(cg%off(xdim:zdim)), int(RED_BLACK, kind=8)), kind=4)
                  do i = i1, cg%ie, id

                     temp = cg%q(soln)%arr(i, j, k) - cg%q(src)%arr(i, j, k)
                     dLdu = 0.

                     im = [i, j, k]

                     do idir = xdim, zdim
                        if (dom%has_dir(idir)) then

                           ih = im ; ih(idir) = ih(idir) + 1
                           call diff_flux(idir, im, soln, cg, cr_id, Keff1)
                           call diff_flux(idir, ih, soln, cg, cr_id, Keff2)

                           temp = temp - (cg%q(qna%wai)%point(ih) - cg%wa(i, j, k)) * diff_theta * dt * cg%idl(idir)
                           dLdu = dLdu - 2 * (Keff1 + Keff2) * cg%idl2(idir)

                        endif
                     enddo

                     ! ToDo add an option to automagically fine-tune overrelax
                     cg%q(soln)%arr(i, j, k) = cg%q(soln)%arr(i, j, k) - overrelax * temp/(one - half * diff_theta * dt * dLdu)

                  enddo
               enddo
            enddo
            cgl => cgl%nxt
         enddo
      enddo

   end subroutine approximate_solution

end module multigrid_diffusion
