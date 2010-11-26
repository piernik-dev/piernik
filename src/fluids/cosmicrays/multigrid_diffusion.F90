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
!<
!! This module contains routines and variables specific for multigrid diffusion solver.
!!
!! Note that for constant-coefficient problems one can use Poisson solver routines, including FFT solver.
!!

module multigrid_diffusion
! pulled by MULTIGRID && COSM_RAYS
#if defined(COSM_RAYS) && defined(MULTIGRID)

   use multigridvars, only: correction, vcycle_stats
   use mpisetup,      only: cbuff_len

   implicit none

   private
   public :: init_multigrid_diff, init_multigrid_diff_post, cleanup_multigrid_diff, multigrid_solve_diff
   public :: diff_tstep_fac, diff_explicit, diff_dt_crs_orig

   ! namelist parameters
   real               :: norm_tol                                     !< stop V-cycle iterations when the ratio of norms ||residual||/||source|| is below this value
   real               :: vcycle_abort                                 !< abort the V-cycle when lhs norm raises by this factor
   integer            :: max_cycles                                   !< Maximum allowed number of V-cycles
   integer            :: nsmool                                       !< smoothing cycles per call
   integer            :: nsmoob                                       !< smoothing cycles on base level when gb_no_fft = .true. (a convergence check would be much better)
   real               :: diff_theta                                   !< 0. is explicit, 1. is fully implicit 0.5 is Crank-Nicholson
   real, protected    :: diff_tstep_fac                               !< How much we stretch timestep. Note that for diff_theta == 0. this should not be > 1.
   logical, protected :: diff_explicit                                !< If .true. then do not use multigrid for diffusion
   real               :: diff_dt_crs_orig                             !< timestep calculated at timestepcosmicrays.F90, before enlarging by diff_tstep_fac
   character(len=cbuff_len) :: diff_bnd_str                           !< Type of diffusion boundary conditions.

   ! mgvar entries for the B field
   integer, parameter :: diff_bx = correction+1, diff_by = diff_bx + 1, diff_bz = diff_by + 1 !< indices pointing to multigrid copies of the b(:,:,:,:) array

   ! miscellaneous
   logical, allocatable, dimension(:) :: norm_was_zero                !< Flag for suppressing repeated warnings on nonexistent CR components
   type(vcycle_stats) :: vstat                                        !< V-cycle statistics
   integer :: diff_extbnd                                             !< external boundary type for relaxation and computation of residuum

contains

!!$ ============================================================================
!!
!! Initialization
!!

   subroutine init_multigrid_diff

      use multigridvars,      only: ngridvars, extbnd_zero, extbnd_extrapolate, extbnd_mirror, extbnd_antimirror
      use mpisetup,           only: buffer_dim, comm, ierr, proc, ibuff, rbuff, lbuff, cbuff, MPI_DOUBLE_PRECISION, MPI_INTEGER, MPI_LOGICAL, MPI_CHARACTER
      use dataio_pub,         only: par_file, ierrh, namelist_errh, compare_namelist      ! QA_WARN required for diff_nml
      use dataio_pub,         only: die, warn, msg

      implicit none

      logical, save                    :: frun = .true.          !< First run flag

      namelist /MULTIGRID_DIFFUSION/ norm_tol, vcycle_abort, max_cycles, nsmool, nsmoob, &
           &                         diff_theta, diff_tstep_fac, diff_explicit, diff_bnd_str

      if (.not.frun) call die("[multigrid_diffusion:init_multigrid_diff] Called more than once.")
      frun = .false.

      ! Default values for namelist variables
      norm_tol       = 1.e-2
      vcycle_abort   = 2.    ! unused as yet
      diff_theta     = 1.
      diff_tstep_fac = 1.
      max_cycles     = 20
      nsmool         = 4
      nsmoob         = 1
      diff_explicit  = .false.
      diff_bnd_str   = "zero"

      if (proc == 0) then

         diff_nml(MULTIGRID_DIFFUSION)

         rbuff(1) = norm_tol
         rbuff(2) = vcycle_abort
         rbuff(3) = diff_theta
         rbuff(4) = diff_tstep_fac

         ibuff(1) = max_cycles
         ibuff(2) = nsmool
         ibuff(3) = nsmoob

         lbuff(1) = diff_explicit

         cbuff(1) = diff_bnd_str

      endif

      call MPI_Bcast(cbuff, cbuff_len*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_Bcast(ibuff,           buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_Bcast(rbuff,           buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)
      call MPI_Bcast(lbuff,           buffer_dim, MPI_LOGICAL,          0, comm, ierr)

      if (proc /= 0) then

         norm_tol       = rbuff(1)
         vcycle_abort   = rbuff(2)
         diff_theta     = rbuff(3)
         diff_tstep_fac = rbuff(4)

         max_cycles     = ibuff(1)
         nsmool         = ibuff(2)
         nsmoob         = ibuff(3)

         diff_explicit  = lbuff(1)

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
         if (diff_theta < 0. .or. diff_theta > 1.) call die("[multigrid_diffusion:init_multigrid] diff_theta must lie in the range [0. .. 1.]")
         if (diff_theta < 0.5 .and. diff_tstep_fac>1. .and. proc == 0) call warn("[multigrid_diffusion:init_multigrid] diff_tstep_fac > 1. for diff_theta < 0.5 might be unstable")
         ! calculate exact limit formula
         ! for diff_theta=0. stable diff_tstep_fac is 0.5 in 2D (guess: 0.333 in 3D)
         ! for diff_theta<0.5 stable diff_tstep_fac rises by 1./(1.-2.*diff_theta)
      endif

   end subroutine init_multigrid_diff

!!$ ============================================================================
!!
!! Initialization - continued after allocation of everything interesting
!!

   subroutine init_multigrid_diff_post(cgrid, mb_alloc)

      use types,              only: grid_container
      use dataio_pub,         only: die
      use fluidindex,         only: nvar
      use multigridhelpers,   only: vcycle_stats_init
      use multigridvars,      only: vcycle_stats

      implicit none

      type(grid_container), intent(in) :: cgrid                  !< copy of grid variables
      real, intent(inout)              :: mb_alloc               !< Allocation counter

      integer                          :: ierr

      if (.false.) ierr = cgrid%nxb ! suppress compiler warnings on unused arguments

      if (allocated(norm_was_zero)) call die("[multigrid_diffusion:init_multigrid] norm_was_zero already allocated")
      allocate(norm_was_zero(nvar%crs%all), stat=ierr)
      if (ierr /= 0) call die("[multigrid_diffusion:init_multigrid] Allocation error: norm_was_zero")
      mb_alloc = mb_alloc + size(norm_was_zero)/2
      norm_was_zero(:) = .false.

      call vcycle_stats_init(vstat, max_cycles)
      mb_alloc = mb_alloc + 2*max_cycles

   end subroutine init_multigrid_diff_post

!!$ ============================================================================
!!
!! Cleanup
!!

   subroutine cleanup_multigrid_diff

      implicit none

      if (allocated(vstat%factor)) deallocate(vstat%factor)
      if (allocated(vstat%time)) deallocate(vstat%time)
      if (allocated(norm_was_zero)) deallocate(norm_was_zero)

   end subroutine cleanup_multigrid_diff

!!$ ============================================================================
!!
!! Multigrid diffusion driver. This is the only multigrid routine intended to be called from the fluidupdate module.
!! This routine is also responsible for communicating the solution to the rest of world
!!

   subroutine multigrid_solve_diff

      use dataio_pub,         only: halfstep, warn, printinfo, msg
      use crdiffusion,        only: cr_diff_x, cr_diff_y, cr_diff_z
      use timer,              only: timer_
      use multigridvars,      only: ts, tot_ts, stdout
      use fluidindex,         only: nvar
      use mpisetup,           only: dt, proc

      implicit none

      logical, save :: frun = .true.
      integer       :: cr_id         ! maybe we should make this variable global in the module and do not pass it as an argument?

      ts =  timer_("multigrid_diffusion", .true.)
      if (diff_explicit) then
         if (frun) then
            if (proc == 0) call warn("[multigrid_diffusion:multigrid_solve_diff] Multigrid was initialized but is not used")
            frun = .false.
         endif
         if (halfstep) then
            call cr_diff_z
            call cr_diff_y
            call cr_diff_x
         else
            call cr_diff_x
            call cr_diff_y
            call cr_diff_z
         endif
      else

         if (dt < 0.99999 * diff_dt_crs_orig * diff_tstep_fac .and. .not. halfstep .and. proc == 0) then
            write(msg,'(a,f8.3,a)')"[multigrid_diffusion:multigrid_solve_diff] Timestep limited somewhere else: dt = ",dt/diff_dt_crs_orig, " of explicit dt_crs."
            call printinfo(msg, stdout)
         endif

         ! set diffBC
         call init_b

         do cr_id = 1, nvar%crs%all
            call init_source(cr_id)
            if (vstat%norm_rhs_orig /= 0) then
               if (norm_was_zero(cr_id) .and. proc == 0) then
                  write(msg,'(a,i2,a)')"[multigrid_diffusion:multigrid_solve_diff] CR-fluid #",cr_id," is now available in measurable quantities."
                  call printinfo(msg)
               endif
               norm_was_zero(cr_id) = .false.

               call init_solution(cr_id)
               ! do substepping
               call vcycle_hg(cr_id)
               ! enddo
            else
               if (.not. norm_was_zero(cr_id) .and. proc == 0) then
                  write(msg,'(a,i2,a)')"[multigrid_diffusion:multigrid_solve_diff] Source norm of CR-fluid #",cr_id," == 0., skipping."
                  call warn(msg)
               endif
               norm_was_zero(cr_id) = .true.
            endif
         enddo

      endif
      ts = timer_("multigrid_diffusion")
      tot_ts = tot_ts + ts

   end subroutine multigrid_solve_diff

!!$ ============================================================================
!!
!! Make a local copy of source
!!

   subroutine init_source(cr_id)

      use multigridvars,      only: roof, source, defect, correction
      use initcosmicrays,     only: iarr_crs
      use grid,               only: is, ie, js, je, ks, ke
      use arrays,             only: u
      use multigridbasefuncs, only: norm_sq
      use dataio_pub,         only: die
      use multigridhelpers,   only: set_dirty, check_dirty

      implicit none

      integer, intent(in) :: cr_id !< CR component index

      call set_dirty(source)
      call set_dirty(correction)
      call set_dirty(defect)
      ! Trick residual subroutine to initialize with: u + (1-theta) dt grad (c grad u)
      if (diff_theta /= 0.) then
         roof%mgvar(roof%is:roof%ie, roof%js:roof%je, roof%ks:roof%ke, correction) = (1. -1./diff_theta) * u(iarr_crs(cr_id), is:ie, js:je, ks:ke)
         roof%mgvar(roof%is:roof%ie, roof%js:roof%je, roof%ks:roof%ke, defect)     =     -1./diff_theta  * u(iarr_crs(cr_id), is:ie, js:je, ks:ke)
         call residual(roof%level, defect, correction, source, cr_id)
      else
         call die("[multigrid_diffusion:init_source] diff_theta = 0 not supported.")
      endif
      call check_dirty(roof%level, source, "init source")

      call norm_sq(source, vstat%norm_rhs_orig)

   end subroutine init_source

!!$ ============================================================================
!!
!! Initialize solution with current CR density (no solution recycling as yet)
!!

   subroutine init_solution(cr_id)

      use multigridvars,  only: roof, solution
      use grid,           only: is, ie, js, je, ks, ke
      use arrays,         only: u
      use initcosmicrays, only: iarr_crs
      use multigridhelpers,   only: set_dirty, check_dirty

      implicit none

      integer, intent(in) :: cr_id !< CR component index

      call set_dirty(solution)
      roof%mgvar(roof%is:roof%ie, roof%js:roof%je, roof%ks:roof%ke, solution) = u(iarr_crs(cr_id),  is:ie, js:je, ks:ke)
      call check_dirty(roof%level, solution, "init solution")

   end subroutine init_solution

!!$ ============================================================================
!!
!! Initialize magnetic field components
!!
!! ToDo: test what happens if we use magnetic field interpolated to the cell centers (some optimizations may then become available)
!!

   subroutine init_b

      use arrays,             only: b
      use grid,               only: is, ie, js, je, ks, ke
      use multigridvars,      only: roof, level_min, level_max, extbnd_mirror, D_x, D_y, D_z
      use multigridbasefuncs, only: restrict_all
      use multigridmpifuncs,  only: mpi_multigrid_bnd
      use fluidindex,         only: ibx, iby, ibz
      use dataio_pub,         only: die
      use multigridhelpers,   only: set_dirty, check_dirty, dirty_label

      implicit none

      integer :: ib, il

      if (diff_bx+iby-ibx /= diff_by .or. diff_bx+ibz-ibx /= diff_bz) call die("[multigrid_diffusion:init_b] Something is wrong with diff_by or diff_bz indices.")

      do ib = ibx, ibz
         call set_dirty(diff_bx+ib-ibx)
         roof%mgvar(roof%is-D_x:roof%ie+D_x, roof%js-D_y:roof%je+D_y, roof%ks-D_z:roof%ke+D_z, diff_bx+ib-ibx) = b(ib, is-D_x:ie+D_x, js-D_y:je+D_y, ks-D_z:ke+D_z)
         call restrict_all(diff_bx+ib-ibx)             ! Implement correct restriction (and probably also separate inter-process communication) routines
         do il = level_min, level_max-1
            call mpi_multigrid_bnd(il, diff_bx+ib-ibx, 1, extbnd_mirror) ! ToDo: use global boundary type for B
            !BEWARE b is set on a staggered grid; corners should be properly set here (now they are not)
            ! the problem is that the b(:,:,:,:) elements are face-centered so restriction and external boundaries should take this into account
            write(dirty_label, '(a,i1)')"init b",ib
            call check_dirty(il, diff_bx+ib-ibx, dirty_label)
         enddo
      enddo

   end subroutine init_b

!!$ ============================================================================
!!
!! Huang-Greengard V-cycle
!!

   subroutine vcycle_hg(cr_id)

      use multigridvars,      only: source, defect, solution, correction, base, roof, level_min, level_max, ts, tot_ts, stdout!, D_x, D_y, D_z
      use multigridbasefuncs, only: norm_sq, restrict_all, prolong_level
      use multigridhelpers,   only: set_dirty, check_dirty, do_ascii_dump, numbered_ascii_dump, mg_write_log, brief_v_log
!      use multigridmpifuncs,  only: mpi_multigrid_bnd
      use initcosmicrays,     only: iarr_crs
      use arrays,             only: u
      use grid,               only: is, ie, js, je, ks, ke
      use dataio_pub,         only: msg, warn
      use mpisetup,           only: proc
      use timer,              only: timer_

      implicit none

      integer, intent(in) :: cr_id !< CR component index

      real, parameter    :: barely_greater_than_1 = 1.05
      integer, parameter :: convergence_history = 2
      integer            :: v, l
      real               :: norm_lhs, norm_rhs, norm_old
      logical            :: dump_every_step

      write(vstat%cprefix,'("D",i1)') cr_id !BEWARE: this is another place with 0 <= cr_id <= 9 limit

      inquire(file = "_dump_every_step_", EXIST=dump_every_step) ! use for debug only
      do_ascii_dump = do_ascii_dump .or. dump_every_step

      call norm_sq(solution, norm_rhs)
      norm_old = norm_rhs

      do v = 0, max_cycles

         call set_dirty(defect)

         call residual(level_max, source, solution, defect, cr_id)
         call norm_sq(defect, norm_lhs)
         ts = timer_("multigrid_diffusion")
         tot_ts = tot_ts + ts

         vstat%count = v
         if (norm_lhs /= 0) then
            vstat%factor(vstat%count) = norm_old/norm_lhs
         else
            vstat%factor(vstat%count) = huge(1.0)
         endif
         vstat%time(vstat%count) = ts

         norm_old = norm_lhs

         if (dump_every_step) call numbered_ascii_dump("md_dump", v)

         if (norm_lhs/norm_rhs <= norm_tol) exit

         if (v>convergence_history) then
            if (product(vstat%factor(v-convergence_history:v)) < barely_greater_than_1 .and. proc == 0) then
               write(msg, '(a,i3,a,g15.5)')"[multigrid_diffusion:vcycle_hg] Too slow convergence: cycle = ",v,", norm_lhs/norm_rhs = ", norm_lhs/norm_rhs
               call warn(msg)
               exit
            endif
         endif

         call restrict_all(defect)

         call set_dirty(correction)
         base%mgvar(:, :, :, correction) = 0.

         do l = level_min, level_max
            call approximate_solution(l, defect, correction, cr_id)
            call prolong_level(l, correction)
         enddo

         call check_dirty(level_max, correction, "c_residual")
         call check_dirty(level_max, defect, "d_residual")
         roof%mgvar     (roof%is:roof%ie, roof%js:roof%je, roof%ks:roof%ke, solution) = &
              roof%mgvar(roof%is:roof%ie, roof%js:roof%je, roof%ks:roof%ke, solution) - &
              roof%mgvar(roof%is:roof%ie, roof%js:roof%je, roof%ks:roof%ke, correction)

      enddo

      if (dump_every_step) call numbered_ascii_dump("md_dump")

      call check_dirty(level_max, solution, "v_soln")

      if (v > max_cycles) then
         if (proc == 0 .and. norm_lhs/norm_rhs > norm_tol) then
            write(msg, '(a,i3,a,g15.5)')"[multigrid_diffusion:vcycle_hg] Not enough V-cycles to achieve convergence: cycle = ",v,", norm_lhs/norm_rhs = ", norm_lhs/norm_rhs
            call warn(msg)
         endif
         v = max_cycles
      endif

      vstat%norm_final = norm_lhs/norm_rhs
      call brief_v_log(vstat)

      call norm_sq(solution, norm_rhs)
      call norm_sq(defect, norm_lhs)
      if (proc == 0 .and. stdout) then
         write(msg,'(a,3g15.5)')"[multigrid_diffusion:vcycle_hg] norms: src, soln, defect: ",vstat%norm_rhs_orig, norm_rhs, norm_lhs
         call mg_write_log(msg)
      endif
!     Do we need to take care of boundaries here?
!      call mpi_multigrid_bnd(roof%level, solution, 1, diff_extbnd)
!      u(iarr_crs(cr_id), is-D_x:ie+D_x, js-D_y:je+D_y, ks-D_z:ke+D_z) = roof%mgvar(roof%is-D_x:roof%ie+D_x, roof%js-D_y:roof%je+D_y, roof%ks-D_z:roof%ke+D_z, solution)
      u(iarr_crs(cr_id), is:ie, js:je, ks:ke) = roof%mgvar(roof%is:roof%ie, roof%js:roof%je, roof%ks:roof%ke, solution)

   end subroutine vcycle_hg

!!$ ============================================================================
!!
!! Compute diffusive flux in the x-direction
!!

   ! BEWARE: almost replicated code (see crdiffusion.F90)
   subroutine diff_flux_x(i, j, k, soln, lev, cr_id, Keff)

      use multigridvars,     only: lvl, YDIR, ZDIR, has_dir
      use initcosmicrays,    only: K_crs_perp, K_crs_paral
      use mpisetup,          only: dt
      use arrays,            only: wa

      implicit none

      integer, intent(in) :: i, j, k      !< cell indices
      integer, intent(in) :: soln         !< multigrid variable to differentiate
      integer, intent(in) :: lev          !< level on which differentiate
      integer, intent(in) :: cr_id        !< CR component index
      real, optional, intent(out) :: Keff !< effective diffusion coefficient for relaxation

      real                :: fcrdif, decr1, decr2, decr3
      real                :: b1b, b2b, b3b, magb

      ! Assumes has_dir(XDIR)
      decr1 = (lvl(lev)%mgvar(i, j, k, soln) - lvl(lev)%mgvar(i-1, j, k, soln)) / lvl(lev)%dx
      fcrdif = K_crs_perp(cr_id) * decr1
      if (present(Keff)) Keff = K_crs_perp(cr_id)

      if (K_crs_paral(cr_id) /= 0.) then

         b1b = lvl(lev)%mgvar(i, j, k, diff_bx)

         if (has_dir(YDIR)) then
            b2b = sum(lvl(lev)%mgvar(i-1:i, j:j+1, k, diff_by))*0.25
            decr2 = ((lvl(lev)%mgvar(i-1, j+1, k, soln) + lvl(lev)%mgvar(i, j+1, k, soln))- &
                 &   (lvl(lev)%mgvar(i-1, j-1, k, soln) + lvl(lev)%mgvar(i, j-1, k, soln))) * 0.25 / lvl(lev)%dy
         else
            b2b = 0.
            decr2 = 0.
         endif

         if (has_dir(ZDIR)) then
            b3b = sum(lvl(lev)%mgvar(i-1:i, j, k:k+1, diff_bz))*0.25
            decr3 = ((lvl(lev)%mgvar(i-1, j, k+1, soln) + lvl(lev)%mgvar(i, j, k+1, soln)) - &
                 &   (lvl(lev)%mgvar(i-1, j, k-1, soln) + lvl(lev)%mgvar(i, j, k-1, soln))) * 0.25 / lvl(lev)%dz
         else
            b3b = 0.
            decr3 = 0.
         endif

         magb = b1b**2 + b2b**2 + b3b**2
         if (magb /= 0.) then
            fcrdif = fcrdif + K_crs_paral(cr_id) * b1b * (b1b*decr1 + b2b*decr2 + b3b*decr3) / magb
            if (present(Keff)) Keff = Keff +  K_crs_paral(cr_id) * b1b**2/magb
         endif

      endif

      wa(i, j, k) = fcrdif * diff_theta * dt / lvl(lev)%dx

   end subroutine diff_flux_x

!!$ ============================================================================
!!
!! Compute diffusive flux in the y-direction
!!

   subroutine diff_flux_y(i, j, k, soln, lev, cr_id, Keff)

      use multigridvars,     only: lvl, XDIR, ZDIR, has_dir
      use initcosmicrays,    only: K_crs_perp, K_crs_paral
      use mpisetup,          only: dt
      use arrays,            only: wa

      implicit none

      integer, intent(in) :: i, j, k      !< cell indices
      integer, intent(in) :: soln         !< multigrid variable to differentiate
      integer, intent(in) :: lev          !< level on which differentiate
      integer, intent(in) :: cr_id        !< CR component index
      real, optional, intent(out) :: Keff !< effective diffusion coefficient for relaxation

      real                :: fcrdif, decr1, decr2, decr3
      real                :: b1b, b2b, b3b, magb

      ! Assumes has_dir(YDIR)
      decr2 = (lvl(lev)%mgvar(i, j, k, soln) - lvl(lev)%mgvar(i, j-1, k, soln)) / lvl(lev)%dy
      fcrdif = K_crs_perp(cr_id) * decr2
      if (present(Keff)) Keff = K_crs_perp(cr_id)

      if (K_crs_paral(cr_id) /= 0.) then

         if (has_dir(XDIR)) then
            b1b = sum(lvl(lev)%mgvar(i:i+1, j-1:j, k, diff_bx))*0.25
            decr1 = ((lvl(lev)%mgvar(i+1, j-1, k, soln) + lvl(lev)%mgvar(i+1, j,  k, soln)) - &
                 &   (lvl(lev)%mgvar(i-1, j-1, k, soln) + lvl(lev)%mgvar(i-1, j,  k, soln))) * 0.25 / lvl(lev)%dx
         else
            b1b = 0.
            decr1 = 0.
         endif

         b2b = lvl(lev)%mgvar(i, j, k, diff_by)

         if (has_dir(ZDIR)) then
            b3b = sum(lvl(lev)%mgvar(i, j-1:j, k:k+1, diff_bz))*0.25
            decr3 = ((lvl(lev)%mgvar(i, j-1, k+1, soln) + lvl(lev)%mgvar(i, j,  k+1, soln)) - &
                 &   (lvl(lev)%mgvar(i, j-1, k-1, soln) + lvl(lev)%mgvar(i, j,  k-1, soln))) * 0.25 / lvl(lev)%dz
         else
            b3b = 0.
            decr3 = 0.
         endif

         magb = b1b**2 + b2b**2 + b3b**2
         if (magb /= 0.) then
            fcrdif = fcrdif + K_crs_paral(cr_id) * b2b * (b1b*decr1 + b2b*decr2 + b3b*decr3) / magb
            if (present(Keff)) Keff = Keff +  K_crs_paral(cr_id) * b2b**2/magb
         endif

      endif

      wa(i, j, k) = fcrdif * diff_theta * dt / lvl(lev)%dy

   end subroutine diff_flux_y

!!$ ============================================================================
!!
!! Compute diffusive flux in the z-direction
!!

   subroutine diff_flux_z(i, j, k, soln, lev, cr_id, Keff)

      use multigridvars,     only: lvl, XDIR, YDIR, has_dir
      use initcosmicrays,    only: K_crs_perp, K_crs_paral
      use mpisetup,          only: dt
      use arrays,            only: wa

      implicit none

      integer, intent(in) :: i, j, k      !< cell indices
      integer, intent(in) :: soln         !< multigrid variable to differentiate
      integer, intent(in) :: lev          !< level on which differentiate
      integer, intent(in) :: cr_id        !< CR component index
      real, optional, intent(out) :: Keff !< effective diffusion coefficient for relaxation

      real                :: fcrdif, decr1, decr2, decr3
      real                :: b1b, b2b, b3b, magb

      ! Assumes has_dir(ZDIR)
      decr3 = (lvl(lev)%mgvar(i, j, k, soln) - lvl(lev)%mgvar(i, j, k-1, soln)) / lvl(lev)%dz
      fcrdif = K_crs_perp(cr_id) * decr3
      if (present(Keff)) Keff = K_crs_perp(cr_id)

      if (K_crs_paral(cr_id) /= 0.) then

         if (has_dir(XDIR)) then
            b1b = sum(lvl(lev)%mgvar(i:i+1, j, k-1:k, diff_bx))*0.25
            decr1 = ((lvl(lev)%mgvar(i+1, j, k-1, soln) + lvl(lev)%mgvar(i+1, j,  k, soln)) - &
                 &   (lvl(lev)%mgvar(i-1, j, k-1, soln) + lvl(lev)%mgvar(i-1, j,  k, soln))) * 0.25 / lvl(lev)%dx
         else
            b1b = 0.
            decr1 = 0.
         endif

         if (has_dir(YDIR)) then
            b2b = sum(lvl(lev)%mgvar(i, j:j+1, k-1:k, diff_by))*0.25
            decr2 = ((lvl(lev)%mgvar(i, j+1, k-1, soln) + lvl(lev)%mgvar(i,  j+1, k, soln)) - &
                 &   (lvl(lev)%mgvar(i, j-1, k-1, soln) + lvl(lev)%mgvar(i,  j-1, k, soln))) * 0.25 / lvl(lev)%dy
         else
            b2b = 0.
            decr2 = 0.
         endif

         b3b =  lvl(lev)%mgvar(i, j, k, diff_bz)

         magb = b1b**2 + b2b**2 + b3b**2
         if (magb /= 0.) then
            fcrdif = fcrdif + K_crs_paral(cr_id) * b3b * (b1b*decr1 + b2b*decr2 + b3b*decr3) / magb
            if (present(Keff)) Keff = Keff +  K_crs_paral(cr_id) * b3b**2/magb
         endif

      endif

      wa(i, j, k) = fcrdif * diff_theta * dt / lvl(lev)%dz

   end subroutine diff_flux_z

!!$ ============================================================================
!!
!! 2nd order: grad (c grad)
!!
!! defect = solution - source - grad (c grad (solution))
!!

   subroutine residual(lev, src, soln, def, cr_id)

      use multigridvars,     only: lvl, XDIR, YDIR, ZDIR, has_dir
      use multigridmpifuncs, only: mpi_multigrid_bnd
      use initcosmicrays,    only: K_crs_perp, K_crs_paral
      use mpisetup,          only: dt
      use arrays,            only: wa
      use multigridhelpers,  only: check_dirty

      implicit none

      integer, intent(in) :: lev   !< level for which approximate the solution
      integer, intent(in) :: src   !< index of source in lvl()%mgvar
      integer, intent(in) :: soln  !< index of solution in lvl()%mgvar
      integer, intent(in) :: def   !< index of defect in lvl()%mgvar
      integer, intent(in) :: cr_id !< CR component index

      integer             :: i, j, k

      call mpi_multigrid_bnd(lev, soln, 1, diff_extbnd) ! corners are required for fluxes

      do k = lvl(lev)%ks, lvl(lev)%ke
         lvl(         lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js  :lvl(lev)%je,   k,   def)                =   &
              &   lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js  :lvl(lev)%je,   k,   soln)               -   &
              &   lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js  :lvl(lev)%je,   k,   src)
      enddo

      if (has_dir(XDIR)) then
         do k = lvl(lev)%ks, lvl(lev)%ke
            do j = lvl(lev)%js, lvl(lev)%je
               do i = lvl(lev)%is, lvl(lev)%ie+1
                  call diff_flux_x(i, j, k, soln, lev, cr_id)
               enddo
            enddo
            lvl     (lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js  :lvl(lev)%je,   k,   def) = &
                 lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js  :lvl(lev)%je,   k,   def) - &
                 ( wa(lvl(lev)%is+1:lvl(lev)%ie+1, lvl(lev)%js  :lvl(lev)%je,   k)        - &
                 & wa(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js  :lvl(lev)%je,   k) )
         enddo
      endif

      if (has_dir(YDIR)) then
         do k = lvl(lev)%ks, lvl(lev)%ke
            do j = lvl(lev)%js, lvl(lev)%je+1
               do i = lvl(lev)%is, lvl(lev)%ie
                  call diff_flux_y(i, j, k, soln, lev, cr_id)
               enddo
            enddo
            lvl     (lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js  :lvl(lev)%je,   k,   def) = &
                 lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js  :lvl(lev)%je,   k,   def) - &
                 ( wa(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js+1:lvl(lev)%je+1, k)        - &
                 & wa(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js  :lvl(lev)%je,   k) )
         enddo
      endif

      if (has_dir(ZDIR)) then
         do k = lvl(lev)%ks, lvl(lev)%ke+1
            do j = lvl(lev)%js, lvl(lev)%je
               do i = lvl(lev)%is, lvl(lev)%ie
                  call diff_flux_z(i, j, k, soln, lev, cr_id)
               enddo
            enddo
         enddo
         lvl     (lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js  :lvl(lev)%je,   lvl(lev)%ks  :lvl(lev)%ke,   def) = &
              lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js  :lvl(lev)%je,   lvl(lev)%ks  :lvl(lev)%ke,   def) - &
              ( wa(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js  :lvl(lev)%je,   lvl(lev)%ks+1:lvl(lev)%ke+1)      - &
              & wa(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js  :lvl(lev)%je,   lvl(lev)%ks  :lvl(lev)%ke  ) )
      endif

      call check_dirty(lev, def, "res def")

   end subroutine residual

!!$ ============================================================================
!!
!! Relaxation.
!!
!! This is the most costly routine in a serial run. Try to find optimal values for nsmool.
!! It seems that for this particular scheme nsmoob can be set even to 1 when the base level is coarsened enough.
!! This routine also depends a lot on communication so it  may limit scalability of the multigrid.
!!

   subroutine approximate_solution(lev, src, soln, cr_id)

      use multigridvars,      only: level_min, has_dir, XDIR, YDIR, ZDIR, lvl, extbnd_donothing
      use multigridmpifuncs,  only: mpi_multigrid_bnd
      use mpisetup,           only: dt
      use arrays,             only: wa

      implicit none

      integer, intent(in) :: lev   !< level for which approximate the solution
      integer, intent(in) :: src   !< index of source in lvl()%mgvar
      integer, intent(in) :: soln  !< index of solution in lvl()%mgvar
      integer, intent(in) :: cr_id !< CR component index

      integer, parameter :: RED_BLACK = 2 !< the checkerboard requires two sweeps

      integer :: n, i, j, k, i1, j1, k1, id, jd, kd
      integer :: nsmoo
      real    :: Keff1, Keff2, dLdu, temp

      if (lev == level_min) then
         nsmoo = nsmoob
      else
         nsmoo = nsmool
      endif

      i1 = lvl(lev)%is; id = 1 ! mv to multigridvars, init_multigrid
      j1 = lvl(lev)%js; jd = 1
      k1 = lvl(lev)%ks; kd = 1
      if (has_dir(XDIR)) then
         id = RED_BLACK
      else if (has_dir(YDIR)) then
         jd = RED_BLACK
      else if (has_dir(ZDIR)) then
         kd = RED_BLACK
      endif

      do n = 1, RED_BLACK*nsmoo
         if (mod(n,2) == 1) then
            call mpi_multigrid_bnd(lev, soln, 1, diff_extbnd) ! no corners are required here
         else
            call mpi_multigrid_bnd(lev, soln, 1, extbnd_donothing)
         endif

         if (kd == RED_BLACK) k1 = lvl(lev)%ks + mod(n, RED_BLACK)
         do k = k1, lvl(lev)%ke, kd
            if (jd == RED_BLACK) j1 = lvl(lev)%js + mod(n+k, RED_BLACK)
            do j = j1, lvl(lev)%je, jd
               if (id == RED_BLACK) i1 = lvl(lev)%is + mod(n+j+k, RED_BLACK)
               do i = i1, lvl(lev)%ie, id

                  temp = lvl(lev)%mgvar(i, j, k, soln) - lvl(lev)%mgvar(i, j, k, src)
                  dLdu = 0.

                  if (has_dir(XDIR)) then

                     call diff_flux_x(i,   j, k, soln, lev, cr_id, Keff1)
                     call diff_flux_x(i+1, j, k, soln, lev, cr_id, Keff2)

                     temp = temp - (wa(i+1, j, k) - wa(i, j, k))
                     dLdu = dLdu - 2 * (Keff1 + Keff2) * lvl(lev)%idx2

                  endif

                  if (has_dir(YDIR)) then

                     call diff_flux_y(i, j,   k, soln, lev, cr_id, Keff1)
                     call diff_flux_y(i, j+1, k, soln, lev, cr_id, Keff2)

                     temp = temp - (wa(i, j+1, k) - wa(i, j, k))
                     dLdu = dLdu - 2 * (Keff1 + Keff2) * lvl(lev)%idy2

                  endif

                  if (has_dir(ZDIR)) then

                     call diff_flux_z(i, j, k,   soln, lev, cr_id, Keff1)
                     call diff_flux_z(i, j, k+1, soln, lev, cr_id, Keff2)

                     temp = temp - (wa(i, j, k+1) - wa(i, j, k))
                     dLdu = dLdu - 2 * (Keff1 + Keff2) * lvl(lev)%idz2

                  endif

                  lvl(lev)%mgvar(i, j, k, soln) = lvl(lev)%mgvar(i, j, k, soln) - temp/(1.e0 - 0.5e0 * diff_theta * dt * dLdu)

               enddo
            enddo
         enddo
      enddo

   end subroutine approximate_solution

#else /* !(COSM_RAYS && MULTIGRID) */
#warning This should not happen. Probably the multigrid_diffusion.F90 file is included in object directory by mistake.
#endif /* !(COSM_RAYS && MULTIGRID) */

end module multigrid_diffusion
