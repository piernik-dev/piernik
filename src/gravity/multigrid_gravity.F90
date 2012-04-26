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
!! \brief Multigrid Poisson solver
!!
!! \details This module contains routines and variables specific for multigrid self-gravity solver.
!!
!! Some code pieces (low-level FFT routines) here are not really gravity-related,
!! but these are not needed for implicit CR-diffusion solver either.
!! These parts of code can be moved to other multigrid files when any other multigrid solver uses them.
!<

module multigrid_gravity
! pulled by MULTIGRID && GRAV

   use constants,     only: cbuff_len
   use multigridvars, only: vcycle_stats

   implicit none

   private
   public :: init_multigrid_grav, init_multigrid_grav_post, cleanup_multigrid_grav, multigrid_solve_grav
   public :: grav_bnd

   include "fftw3.f"
   ! constants from fftw3.f
   !   integer, parameter :: FFTW_MEASURE=0, FFTW_PATIENT=32, FFTW_ESTIMATE=64
   !   integer, parameter :: FFTW_RODFT01=8, FFTW_RODFT10=9

   ! multigrid constants
   enum, bind(C)
      enumerator :: fft_rcr = 1                                       !< type of FFT transform: full
      enumerator :: fft_dst                                           !< type of FFT transform: discrete sine
      enumerator :: fft_none=-1                                       !< type of FFT transform: none
   end enum

   ! namelist parameters
   real               :: norm_tol                                     !< stop V-cycle iterations when the ratio of norms ||residual||/||source|| is below this value
   real               :: overrelax                                    !< overrealaxation factor (if < 1. then works as underrelaxation), provided for tests
   real               :: overrelax_x                                  !< x-direction overrealaxation factor for fine tuning convergence ratio when cell spacing is not equal in all 3 directions. Use with care, patience and lots of hope.
   real               :: overrelax_y                                  !< y-direction overrealaxation factor for fine tuning convergence ratio when cell spacing is not equal in all 3 directions. Use with care, patience and lots of hope.
   real               :: overrelax_z                                  !< z-direction overrealaxation factor for fine tuning convergence ratio when cell spacing is not equal in all 3 directions. Use with care, patience and lots of hope.
   real               :: Jacobi_damp                                  !< omega factor for damped Jacobi relaxation. Jacobi_damp == 1 gives undamped method. Try 0.5 in 1D.
   real               :: vcycle_abort                                 !< abort the V-cycle when lhs norm raises by this factor
   real               :: L4_strength                                  !< strength of the 4th order terms in the Laplace operator; 0.: 2nd, 1.: 4th direct, 0.5: 4th integral
   integer            :: max_cycles                                   !< Maximum allowed number of V-cycles
   integer            :: nsmool                                       !< smoothing cycles per call
   integer            :: nsmoob                                       !< smoothing cycles on base level when cannot use FFT. (a convergence check would be much better)
   integer            :: nsmoof                                       !< FFT iterations per call
   integer            :: ord_laplacian                                !< Laplace operator order; allowed values are 2 (default) and 4 (experimental, not fully implemented)
   integer            :: ord_time_extrap                              !< Order of temporal extrapolation for solution recycling; -1 means 0-guess, 2 does parabolic interpolation
   logical            :: trust_fft_solution                           !< Bypass the V-cycle, when doing FFT on whole domain, make sure first that FFT is properly set up.
   logical            :: base_no_fft                                  !< Deny solving the base level with FFT. Can be very slow.
   logical            :: prefer_rbgs_relaxation                       !< Prefer relaxation over FFT local solver. Typically faster.
   !> \todo allow to perform one or more V-cycles with FFT method, the switch to the RBGS (may save one V-cycle in some cases)
   logical            :: fft_full_relax                               !< Perform full or boundary relaxation after local FFT solve
   logical            :: fft_patient                                  !< Spend more time in init_multigrid to find faster fft plan
   logical            :: multidim_code_3D                             !< prefer code written for any 1D and 2D configuration even in 3D for benchmarking and debugging
   character(len=cbuff_len) :: grav_bnd_str                           !< Type of gravitational boundary conditions.
   logical            :: require_FFT                                  !< .true. if we use FFT solver anywhere (and need face prolongation)

   ! boundaries
   integer            :: grav_bnd                                     !< boundary type for computational domain
   !integer           :: grav_extbnd_mode                             !< external boundary mode

   integer            :: fftw_flags = FFTW_MEASURE                    !< or FFTW_PATIENT on request

   ! solution recycling
   integer, parameter :: nold_max=3                                   !< maximum implemented extrapolation order
   integer :: nold                                                    !< number of old solutions kept for solution recycling
   type :: old_soln                                                   !< container for an old solution with its timestamp
      integer :: i_hist                                               !< index to the old solution
      real :: time                                                    !< time of the old solution
   end type old_soln
   type :: soln_history                                               !< container for a set of several old potential solutions
      type(old_soln), dimension(nold_max) :: old
      integer :: last                                                 !< index of the last stored potential
      logical :: valid                                                !< .true. when old(last) was properly initialized
    contains
      procedure :: init_history
   end type soln_history
   type(soln_history), target :: inner, outer                         !< storage for recycling the inner and outer potentials

   ! miscellaneous
   type(vcycle_stats) :: vstat                                        !< V-cycle statistics

contains

!!$ ============================================================================
!!
!! Initialization
!!
!>
!! \brief Routine to set parameters values from namelist MULTIGRID_GRAVITY
!!
!! \n \n
!! @b MULTIGRID_GRAVITY
!! \n \n
!! <table border="+1">
!! <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!! <tr><td>norm_tol              </td><td>1.e-6  </td><td>real value     </td><td>\copydoc multigrid_gravity::norm_tol              </td></tr>
!! <tr><td>vcycle_abort          </td><td>2.0    </td><td>real value     </td><td>\copydoc multigrid_gravity::vcycle_abort          </td></tr>
!! <tr><td>max_cycles            </td><td>20     </td><td>integer value  </td><td>\copydoc multigrid_gravity::max_cycles            </td></tr>
!! <tr><td>nsmool                </td><td>4      </td><td>integer value  </td><td>\copydoc multigrid_gravity::nsmool                </td></tr>
!! <tr><td>nsmoob                </td><td>100    </td><td>integer value  </td><td>\copydoc multigrid_gravity::nsmoob                </td></tr>
!! <tr><td>overrelax             </td><td>1.     </td><td>real value     </td><td>\copydoc multigrid_gravity::overrelax             </td></tr>
!! <tr><td>overrelax_x           </td><td>1.     </td><td>real value     </td><td>\copydoc multigrid_gravity::overrelax_x           </td></tr>
!! <tr><td>overrelax_y           </td><td>1.     </td><td>real value     </td><td>\copydoc multigrid_gravity::overrelax_y           </td></tr>
!! <tr><td>overrelax_z           </td><td>1.     </td><td>real value     </td><td>\copydoc multigrid_gravity::overrelax_z           </td></tr>
!! <tr><td>Jacobi_damp           </td><td>1.     </td><td>real value     </td><td>\copydoc multigrid_gravity::Jacobi_damp           </td></tr>
!! <tr><td>L4_strength           </td><td>1.0    </td><td>real value     </td><td>\copydoc multigrid_gravity::L4_strength           </td></tr>
!! <tr><td>nsmoof                </td><td>1      </td><td>integer value  </td><td>\copydoc multigrid_gravity::nsmoof                </td></tr>
!! <tr><td>ord_laplacian         </td><td>2      </td><td>integer value  </td><td>\copydoc multigrid_gravity::ord_laplacian         </td></tr>
!! <tr><td>ord_time_extrap       </td><td>1      </td><td>integer value  </td><td>\copydoc multigrid_gravity::ord_time_extrap       </td></tr>
!! <tr><td>prefer_rbgs_relaxation</td><td>.true. </td><td>logical        </td><td>\copydoc multigrid_gravity::prefer_rbgs_relaxation</td></tr>
!! <tr><td>base_no_fft           </td><td>.false.</td><td>logical        </td><td>\copydoc multigrid_gravity::base_no_fft           </td></tr>
!! <tr><td>fft_full_relax        </td><td>.false.</td><td>logical        </td><td>\copydoc multigrid_gravity::fft_full_relax        </td></tr>
!! <tr><td>fft_patient           </td><td>.false.</td><td>logical        </td><td>\copydoc multigrid_gravity::fft_patient           </td></tr>
!! <tr><td>trust_fft_solution    </td><td>.false.</td><td>logical        </td><td>\copydoc multigrid_gravity::trust_fft_solution    </td></tr>
!! <tr><td>coarsen_multipole     </td><td>1      </td><td>integer value  </td><td>\copydoc multipole::coarsen_multipole             </td></tr>
!! <tr><td>lmax                  </td><td>16     </td><td>integer value  </td><td>\copydoc multipole::lmax                          </td></tr>
!! <tr><td>mmax                  </td><td>-1     </td><td>integer value  </td><td>\copydoc multipole::mmax                          </td></tr>
!! <tr><td>ord_prolong_mpole     </td><td>-2     </td><td>integer value  </td><td>\copydoc multipole::ord_prolong_mpole             </td></tr>
!! <tr><td>use_point_monopole    </td><td>.false.</td><td>logical        </td><td>\copydoc multipole::use_point_monopole            </td></tr>
!! <tr><td>interp_pt2mom         </td><td>.false.</td><td>logical        </td><td>\copydoc multipole::interp_pt2mom                 </td></tr>
!! <tr><td>interp_mom2pot        </td><td>.false.</td><td>logical        </td><td>\copydoc multipole::interp_mom2pot                </td></tr>
!! <tr><td>multidim_code_3D      </td><td>.false.</td><td>logical        </td><td>\copydoc multigrid_gravity::multidim_code_3D      </td></tr>
!! <tr><td>grav_bnd_str          </td><td>"periodic"/"dirichlet"</td><td>string of chars</td><td>\copydoc multigrid_gravity::grav_bnd_str          </td></tr>
!! </table>
!! \n \n
!<
   subroutine init_multigrid_grav

      use constants,     only: GEO_XYZ, GEO_RPZ, BND_PER
      use dataio_pub,    only: par_file, ierrh, namelist_errh, compare_namelist, cmdl_nml, lun  ! QA_WARN required for diff_nml
      use dataio_pub,    only: msg, die, warn
      use domain,        only: dom, is_uneven, is_multicg
      use mpi,           only: MPI_CHARACTER, MPI_DOUBLE_PRECISION, MPI_INTEGER, MPI_LOGICAL, MPI_COMM_NULL
      use mpisetup,      only: buffer_dim, comm, ierr, master, slave, ibuff, cbuff, rbuff, lbuff, FIRST
      use multigridvars, only: bnd_periodic, bnd_dirichlet, bnd_isolated, bnd_invalid, single_base
      use multipole,     only: use_point_monopole, lmax, mmax, ord_prolong_mpole, coarsen_multipole, interp_pt2mom, interp_mom2pot
      use types,         only: cdd

      implicit none

      integer       :: periodic_bnd_cnt   !< counter of periodic boundaries in existing directions
      logical, save :: frun = .true.      !< First run flag

      namelist /MULTIGRID_GRAVITY/ norm_tol, vcycle_abort, max_cycles, nsmool, nsmoob, &
           &                       overrelax, overrelax_x, overrelax_y, overrelax_z, Jacobi_damp, L4_strength, nsmoof, ord_laplacian, ord_time_extrap, &
           &                       prefer_rbgs_relaxation, base_no_fft, fft_full_relax, fft_patient, trust_fft_solution, &
           &                       coarsen_multipole, lmax, mmax, ord_prolong_mpole, use_point_monopole, interp_pt2mom, interp_mom2pot, multidim_code_3D, &
           &                       grav_bnd_str

      if (.not.frun) call die("[multigrid_gravity:init_multigrid_grav] Called more than once.")
      frun = .false.

      ! Default values for namelist variables
      norm_tol               = 1.e-6
      overrelax              = 1.
      overrelax_x            = 1.
      overrelax_y            = 1.
      overrelax_z            = 1.
      Jacobi_damp            = 1.
      vcycle_abort           = 2.
      L4_strength            = 1.0

      coarsen_multipole      = 1
      if (is_uneven) coarsen_multipole = 0
      lmax                   = 16
      mmax                   = -1 ! will be automatically set to lmax unless explicitly limited in problem.par
      max_cycles             = 20
      nsmool                 = 4
      nsmoob                 = 100
      nsmoof                 = 1
      ord_laplacian          = 2
      ord_prolong_mpole      = -2
      ord_time_extrap        = 1

      use_point_monopole     = .false.
      trust_fft_solution     = .false.
      base_no_fft            = is_multicg
      prefer_rbgs_relaxation = .true. !> \warning it seems that FFT local solver is broken somewhere. ToDo: restore .false. as a default
      fft_full_relax         = .false.
      fft_patient            = .false.
      interp_pt2mom          = .false.
      interp_mom2pot         = .false.
      multidim_code_3D       = .false.

      periodic_bnd_cnt = count(dom%periodic(:) .and. dom%has_dir(:))

      if (periodic_bnd_cnt == dom%eff_dim) then
         grav_bnd_str = "periodic"
      else
         grav_bnd_str = "dirichlet"
      endif

      if (master) then

         diff_nml(MULTIGRID_GRAVITY)

         ! FIXME when ready
         if (dom%geometry_type == GEO_RPZ) then
            call warn("[multigrid_gravity:init_multigrid_grav] cylindrical geometry support is under development.")
            ! switch off FFT-related bits
            base_no_fft = .true.
            prefer_rbgs_relaxation = .true.
            ord_laplacian = 2
            L4_strength = 0.
            ! ord_prolong_mpole = 0
         else if (dom%geometry_type /= GEO_XYZ) then
            call die("[multigrid_gravity:init_multigrid_grav] non-cartesian geometry not implemented yet.")
         endif

         if (is_multicg .and. .not. base_no_fft) then
            call warn("[multigrid_gravity:init_multigrid_grav] base_no_fft set to .true. for multicg configuration")
            base_no_fft = .true.
         endif

         rbuff(1) = norm_tol
         rbuff(2) = overrelax
         rbuff(3) = overrelax_x
         rbuff(4) = overrelax_y
         rbuff(5) = overrelax_z
         rbuff(6) = Jacobi_damp
         rbuff(7) = vcycle_abort
         rbuff(8) = L4_strength

         ibuff( 1) = coarsen_multipole
         ibuff( 2) = lmax
         ibuff( 3) = mmax
         ibuff( 4) = max_cycles
         ibuff( 5) = nsmool
         ibuff( 6) = nsmoob
         ibuff( 7) = nsmoof
         ibuff( 8) = ord_laplacian
         ibuff( 9) = ord_prolong_mpole
         ibuff(10) = ord_time_extrap

         lbuff(1) = use_point_monopole
         lbuff(2) = trust_fft_solution
         lbuff(3) = base_no_fft
         lbuff(4) = prefer_rbgs_relaxation
         lbuff(5) = fft_full_relax
         lbuff(6) = fft_patient
         lbuff(7) = interp_pt2mom
         lbuff(8) = interp_mom2pot
         lbuff(9) = multidim_code_3D

         cbuff(1) = grav_bnd_str

      endif

      call MPI_Bcast(cbuff, cbuff_len*buffer_dim, MPI_CHARACTER,        FIRST, comm, ierr)
      call MPI_Bcast(ibuff,           buffer_dim, MPI_INTEGER,          FIRST, comm, ierr)
      call MPI_Bcast(rbuff,           buffer_dim, MPI_DOUBLE_PRECISION, FIRST, comm, ierr)
      call MPI_Bcast(lbuff,           buffer_dim, MPI_LOGICAL,          FIRST, comm, ierr)

      if (slave) then

         norm_tol       = rbuff(1)
         overrelax      = rbuff(2)
         overrelax_x    = rbuff(3)
         overrelax_y    = rbuff(4)
         overrelax_z    = rbuff(5)
         Jacobi_damp    = rbuff(6)
         vcycle_abort   = rbuff(7)
         L4_strength    = rbuff(8)

         coarsen_multipole = ibuff( 1)
         lmax              = ibuff( 2)
         mmax              = ibuff( 3)
         max_cycles        = ibuff( 4)
         nsmool            = ibuff( 5)
         nsmoob            = ibuff( 6)
         nsmoof            = ibuff( 7)
         ord_laplacian     = ibuff( 8)
         ord_prolong_mpole = ibuff( 9)
         ord_time_extrap   = ibuff(10)

         use_point_monopole      = lbuff(1)
         trust_fft_solution      = lbuff(2)
         base_no_fft             = lbuff(3)
         prefer_rbgs_relaxation  = lbuff(4)
         fft_full_relax          = lbuff(5)
         fft_patient             = lbuff(6)
         interp_pt2mom           = lbuff(7)
         interp_mom2pot          = lbuff(8)
         multidim_code_3D        = lbuff(9)

         grav_bnd_str   = cbuff(1)(1:len(grav_bnd_str))

      endif

      if (dom%geometry_type == GEO_RPZ) multidim_code_3D = .true. ! temporarily

      ! boundaries
      grav_bnd = bnd_invalid
      select case (grav_bnd_str)
         case ("isolated", "iso")
            grav_bnd = bnd_isolated
         case ("periodic", "per")
            if (any(dom%bnd(:,:) /= BND_PER)) &
                 call die("[multigrid_gravity:init_multigrid_grav] cannot enforce periodic boundaries for gravity on a not fully periodic domain")
            grav_bnd = bnd_periodic
         case ("dirichlet", "dir")
            grav_bnd = bnd_dirichlet
         case default
            call die("[multigrid_gravity:init_multigrid_grav] Non-recognized boundary description.")
      end select

      if (periodic_bnd_cnt == dom%eff_dim) then ! fully periodic domain
         if (grav_bnd /= bnd_periodic .and. master) call warn("[multigrid_gravity:init_multigrid_grav] Ignoring non-periodic boundary conditions for gravity on a fully periodic domain.")
         grav_bnd = bnd_periodic
      else if (periodic_bnd_cnt > 0 .and. periodic_bnd_cnt < dom%eff_dim) then
         if (master) call warn("[multigrid_gravity:init_multigrid_grav] Mixing periodic and non-periodic boundary conditions for gravity is experimental.")
         prefer_rbgs_relaxation = .true.
         base_no_fft = .true.
      endif
!!$      select case (grav_bnd)
!!$         case (bnd_periodic)
!!$            grav_extbnd_mode = extbnd_donothing
!!$         case (bnd_isolated, bnd_dirichlet, bnd_givenval)
!!$            grav_extbnd_mode = extbnd_antimirror
!!$         case default
!!$            call die("[multigrid_gravity:init_multigrid_grav] Unsupported grav_bnd.")
!!$            !grav_extbnd_mode = extbnd_donothing
!!$      end select

      if (.not. (grav_bnd == bnd_periodic .or. grav_bnd == bnd_dirichlet .or. grav_bnd == bnd_isolated) .and. .not. base_no_fft) then
         base_no_fft = .true.
         if (master) call warn("[multigrid_gravity:init_multigrid_grav] Use of FFT not allowed by current boundary type/combination.")
      endif

      ! something is a bit messed up here
      if (cdd%comm3d /= MPI_COMM_NULL .and. .not. base_no_fft) then
         base_no_fft = .true.
         if (master) call warn("[multigrid_gravity:init_multigrid_grav] cdd%comm3d disables use of FFT at base level")
      endif

      single_base = .not. base_no_fft

      if (.not. prefer_rbgs_relaxation .and. any([ overrelax, overrelax_x, overrelax_y, overrelax_z ] /= 1.)) then
         if (master) call warn("[multigrid_gravity:init_multigrid_grav] Overrelaxation is disabled for FFT local solver.")
         overrelax = 1.
         overrelax_x   = 1.
         overrelax_y   = 1.
         overrelax_z   = 1.
      endif

      if (master) then
         if ((Jacobi_damp <= 0. .or. Jacobi_damp>1.)) then
            write(msg, '(a,g12.5,a)')"[multigrid_gravity:init_multigrid_grav] Jacobi_damp = ",Jacobi_damp," is outside (0, 1] interval."
            call warn(msg)
         endif
         if (overrelax /= 1. .or. overrelax_x /= 1. .or. overrelax_y /= 1. .or. overrelax_z /= 1.) then
            write(msg, '(a,f8.5,a,3f8.5,a)')"[multigrid_gravity:init_multigrid_grav] Overrelaxation factors: global = ", overrelax, ", directional = [", overrelax_x, overrelax_y, overrelax_z, "]"
            call warn(msg)
         endif
      endif

      if (fft_patient) fftw_flags = FFTW_PATIENT

   end subroutine init_multigrid_grav

!!$ ============================================================================
!>
!! \brief Initialization - continued after allocation of everything interesting
!<

   subroutine init_multigrid_grav_post

      use constants,     only: pi, dpi, GEO_XYZ, one, zero, half, sgp_n
      use dataio_pub,    only: die, warn
      use domain,        only: dom
      use gc_list,       only: cg_list_element, all_cg
      use cg_list_lev,   only: cg_list_level
      use grid,          only: leaves
      use grid_cont,     only: grid_container
      use mpi,           only: MPI_COMM_NULL
      use mpisetup,      only: master, nproc, proc
      use multigridvars, only: roof, base, bnd_periodic, bnd_dirichlet, bnd_isolated, is_mg_uneven, need_general_pf, single_base
      use multipole,     only: init_multipole, coarsen_multipole
      use types,         only: cdd

      implicit none

      real, allocatable, dimension(:)  :: kx, ky, kz             !< FFT kernel directional components for convolution
      integer :: i, j
      type(cg_list_level), pointer :: curl
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg

      need_general_pf = cdd%comm3d == MPI_COMM_NULL .or. single_base .or. is_mg_uneven

      if (need_general_pf .and. coarsen_multipole /= 0) then
         coarsen_multipole = 0
         if (master) call warn("[multigrid_gravity:init_multigrid_grav_post] multipole coarsening on uneven domains or with cdd%comm3d == MPI_COMM_NULL is not implemented yet.")
      endif

      ! solution recycling
      ord_time_extrap = min(nold_max-1, max(-1, ord_time_extrap))
      nold = ord_time_extrap + 1
      if (nold > 0) then
         call inner%init_history(nold, "i")
         if (grav_bnd == bnd_isolated) call outer%init_history(nold, "o")
      endif

      call leaves%set_q_value(all_cg%ind(sgp_n), 0.) !Initialize all the guardcells, even those which does not impact the solution

      curl => base
      do while (associated(curl))
         ! this should work correctly also when dom%eff_dim < 3
         cgl => curl%first
         do while (associated(cgl))
            cg => cgl%cg

            cg%mg%r  = overrelax   / 2.
            cg%mg%rx = cg%dvol2 * cg%idx2
            cg%mg%ry = cg%dvol2 * cg%idy2
            cg%mg%rz = cg%dvol2 * cg%idz2
            cg%mg%r  = cg%mg%r  / (cg%mg%rx + cg%mg%ry + cg%mg%rz)
            cg%mg%rx = overrelax_x * cg%mg%rx * cg%mg%r
            cg%mg%ry = overrelax_y * cg%mg%ry * cg%mg%r
            cg%mg%rz = overrelax_z * cg%mg%rz * cg%mg%r
            cg%mg%r  = cg%mg%r  * cg%dvol2
            !>
            !! \deprecated BEWARE: some of the above invariants may be not optimally defined - the convergence ratio drops when dx /= dy or dy /= dz or dx /= dz
            !! and overrelaxation factors are required to get any convergence (often poor)
            !<
            if (prefer_rbgs_relaxation) then
               cg%mg%fft_type = fft_none
            else if (grav_bnd == bnd_periodic .and. nproc == 1) then
               cg%mg%fft_type = fft_rcr
            else if (grav_bnd == bnd_periodic .or. grav_bnd == bnd_dirichlet .or. grav_bnd == bnd_isolated) then
               cg%mg%fft_type = fft_dst
            else
               cg%mg%fft_type = fft_none
            endif

            cgl => cgl%nxt
         enddo
         curl => curl%finer
      enddo

      ! data related to local and global base-level FFT solver
      cgl => base%first
      do while (associated(cgl))
         if (base_no_fft) then
            cgl%cg%mg%fft_type = fft_none
         else
            select case (grav_bnd)
               case (bnd_periodic)
                  cgl%cg%mg%fft_type = fft_rcr
               case (bnd_dirichlet, bnd_isolated)
                  cgl%cg%mg%fft_type = fft_dst
               case default
                  cgl%cg%mg%fft_type = fft_none
                  if (master) call warn("[multigrid_gravity:init_multigrid_grav_post] base_no_fft set but no suitable boundary conditions found. Reverting to RBGS relaxation.")
            end select
         endif
         cgl => cgl%nxt
      enddo

      !special initialization of global base-level FFT-related data
      if (dom%geometry_type /= GEO_XYZ) then
         curl => base
         do while (associated(curl))
            cgl => curl%first
            do while (associated(cgl))
               if (cgl%cg%mg%fft_type /= fft_none) call die("[multigrid_gravity:init_multigrid_grav_post] FFT is not allowed in non-cartesian coordinates.")
               cgl => cgl%nxt
            enddo
            curl => curl%finer
         enddo
      endif

      require_FFT = .false.

      ! FFT solver storage and data
      curl => base
      do while (associated(curl))

         cgl => curl%first
         do while (associated(cgl))
            cg => cgl%cg
            cg%mg%planf = 0
            cg%mg%plani = 0

            if (cg%mg%fft_type /= fft_none) then

               require_FFT = .true.

               if (dom%geometry_type /= GEO_XYZ) call die("[multigrid_gravity:init_multigrid_grav_post] FFT is not allowed in non-cartesian coordinates.")

               select case (cg%mg%fft_type)
                  case (fft_rcr)
                     cg%mg%nxc = cg%nxb / 2 + 1
                  case (fft_dst)
                     cg%mg%nxc = cg%nxb
                  case default
                     call die("[multigrid_gravity:init_multigrid_grav_post] Unknown FFT type.")
               end select

               if (allocated(cg%mg%Green3D) .or. allocated(cg%mg%src)) call die("[multigrid_gravity:init_multigrid_grav_post] Green3D or src arrays already allocated")
               allocate(cg%mg%Green3D(cg%mg%nxc, cg%nyb, cg%nzb))
               allocate(cg%mg%src    (cg%nxb,    cg%nyb, cg%nzb))

               if (allocated(kx)) deallocate(kx)
               if (allocated(ky)) deallocate(ky)
               if (allocated(kz)) deallocate(kz)
               allocate(kx(cg%mg%nxc), ky(cg%nyb), kz(cg%nzb))

               select case (cg%mg%fft_type)

                  ! cg%mg%fft_norm is set such that the following sequence gives identity:
                  ! call dfftw_execute(cg%mg%planf); cg%mg%fftr(:, :, :) = cg%mg%fftr(:, :, :) * cg%mg%fft_norm ; call dfftw_execute(cg%mg%plani)

                  case (fft_rcr)
                     if (allocated(cg%mg%fft)) call die("[multigrid_gravity:init_multigrid_grav_post] fft or Green3D array already allocated")
                     allocate(cg%mg%fft(cg%mg%nxc, cg%nyb, cg%nzb))

                     cg%mg%fft_norm = one / real( product(cg%n_b(:), mask=dom%has_dir(:)) ) ! No 4 pi G factor here because the source was already multiplied by it

                     ! FFT local solver initialization for 2nd order (3-point) Laplacian
                     ! sin(k*x-d) - 2.*sin(k*x) + sin(k*x+d) = 2 * (cos(d)-1) * sin(k*x) = -4 * sin(d/2)**2 * sin(k*x)
                     ! For 4th order: a*sin(k*x) + b*(sin(k*x-d) + sin(k*x+d)) + c*(sin(k*x-2*d) + sin(k*x+2*d)), a+2*b+2*c == 0 it would be:
                     ! 4*(a+b+(a+2*b)*cos(d)) * sin(d/2)**2 * sin(k*x)
                     ! For 6th order: a*sin(k*x) + b*(sin(k*x-d) + sin(k*x+d)) + c*(sin(k*x-2*d) + sin(k*x+2*d)) + e*(sin(k*x-3*d) + sin(k*x+3*d)), a+2*b+2*c+2*e == 0 it would be:
                     ! 2*(3*a+4*b+2*c+4*(a+2*b+c)*cos(d)+2*(a+2*(b+c))*cos(2*d)) * sin(d/2)**2 * sin(k*x)
                     ! asymptotically: -d**2/2 for d<pi

                     kx(:) = cg%idx2 * (cos(dpi/cg%nxb*[( j, j=0, cg%mg%nxc-1 )]) - one)
                     ky(:) = cg%idy2 * (cos(dpi/cg%nyb*[( j, j=0, cg%nyb-1 )]) - one)
                     kz(:) = cg%idz2 * (cos(dpi/cg%nzb*[( j, j=0, cg%nzb-1 )]) - one)
                     call dfftw_plan_dft_r2c_3d(cg%mg%planf, cg%nxb, cg%nyb, cg%nzb, cg%mg%src, cg%mg%fft, fftw_flags)
                     call dfftw_plan_dft_c2r_3d(cg%mg%plani, cg%nxb, cg%nyb, cg%nzb, cg%mg%fft, cg%mg%src, fftw_flags)

                  case (fft_dst)

                     if (allocated(cg%mg%fftr)) call die("[multigrid_gravity:init_multigrid_grav_post] fftr array already allocated")
                     allocate(cg%mg%fftr(cg%mg%nxc, cg%nyb, cg%nzb))

                     cg%mg%fft_norm = one / (8. * real( product(cg%n_b(:), mask=dom%has_dir(:)) ))
                     kx(:) = cg%idx2 * (cos(pi/cg%nxb*[( j, j=1, cg%mg%nxc )]) - one)
                     ky(:) = cg%idy2 * (cos(pi/cg%nyb*[( j, j=1, cg%nyb )]) - one)
                     kz(:) = cg%idz2 * (cos(pi/cg%nzb*[( j, j=1, cg%nzb )]) - one)
                     call dfftw_plan_r2r_3d(cg%mg%planf, cg%nxb, cg%nyb, cg%nzb, cg%mg%src,  cg%mg%fftr, FFTW_RODFT10, FFTW_RODFT10, FFTW_RODFT10, fftw_flags)
                     call dfftw_plan_r2r_3d(cg%mg%plani, cg%nxb, cg%nyb, cg%nzb, cg%mg%fftr, cg%mg%src,  FFTW_RODFT01, FFTW_RODFT01, FFTW_RODFT01, fftw_flags)

                  case default
                     call die("[multigrid_gravity:init_multigrid_grav_post] Unknown FFT type.")
               end select

               ! compute Green's function for 7-point 3D discrete laplacian
               do i = 1, cg%mg%nxc
                  do j = 1, cg%nyb
                     where ( (kx(i) + ky(j) + kz(:)) /= 0 )
                        cg%mg%Green3D(i,j,:) = half * cg%mg%fft_norm / (kx(i) + ky(j) + kz(:))
                     elsewhere
                        cg%mg%Green3D(i,j,:) = zero
                     endwhere
                  enddo
               enddo

            endif

            cgl => cgl%nxt
         enddo
         curl => curl%finer
      enddo

      if (require_FFT) call mpi_multigrid_prep_grav !supplement to mpi_multigrid_prep

      cgl => roof%first
      do while (associated(cgl))
         if (cgl%cg%mg%fft_type == fft_none .and. trust_fft_solution) then
            if (master) call warn("[multigrid_gravity:init_multigrid_grav_post] cannot trust FFT solution on the roof.")
            trust_fft_solution = .false.
         endif
         cgl => cgl%nxt
      enddo

      if (allocated(kx)) deallocate(kx)
      if (allocated(ky)) deallocate(ky)
      if (allocated(kz)) deallocate(kz)

      if (grav_bnd == bnd_isolated) call init_multipole

      call vstat%init(max_cycles)

   end subroutine init_multigrid_grav_post

!> \brief Initialize structure for keeping historical potential fields

   subroutine init_history(this, nold, prefix)

      use constants,        only: singlechar, dsetnamelen, AT_IGNORE
      use gc_list,          only: all_cg
      use multigridhelpers, only: set_dirty

      implicit none

      class(soln_history),       intent(inout) :: this   !< COMMENT ME
      integer,                   intent(in)    :: nold   !< how many historic points to store
      character(len=singlechar), intent(in)    :: prefix !< 'i' for inner, 'o' for outer potential

      integer :: i
      character(len=dsetnamelen) :: hname

      do i = 1, nold
         write(hname,'(2a,i2.2)')prefix,"_h_",i
         call all_cg%reg_var(hname, .true., AT_IGNORE, multigrid = .true.)
         this%old(i) = old_soln(all_cg%ind(hname), -huge(1.0))
         call set_dirty(this%old(i)%i_hist)
      enddo
      this%valid = .false.
      this%last  = 1

   end subroutine init_history

!!$ ============================================================================
!>
!! \brief Set up communication for face prolongation
!!
!! \todo implement also prolongation of coarsened multipoles
!!
!! \todo move this to cg_list_level
!<

   subroutine mpi_multigrid_prep_grav

      use constants,     only: xdim, ydim, zdim, ndims, LO, HI, LONG, zero, one, half
      use dataio_pub,    only: warn, die
      use domain,        only: dom
      use cg_list_lev,   only: cg_list_level
      use gc_list,       only: cg_list_element
      use grid_cont,     only: pr_segment, grid_container, is_overlap
      use mpisetup,      only: proc, nproc, FIRST, LAST, procmask, inflate_req, req
      use multigridvars, only: base, ord_prolong_face_norm, need_general_pf
#ifdef DEBUG
      use constants,     only: two
      use piernikdebug,  only: aux_R
#endif /* DEBUG */

      implicit none

      integer :: d, g, j, lh, hl, l, nl
      integer(kind=8), dimension(xdim:zdim) :: ijks, per
      logical, dimension(xdim:zdim) :: dmask
      integer(kind=8), dimension(xdim:zdim, LO:HI) :: coarsened, b_layer
      type(pr_segment), pointer :: seg
      type(cg_list_level), pointer   :: curl                   !> current level (a pointer sliding along the linked list)
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg            !< current grid container
      logical :: is_internal_fine
      integer, parameter :: max_opfn = 2
      real, dimension(0:max_opfn) :: opfn_c_ff, opfn_c_cf

      if (.not. need_general_pf) return

      call inflate_req(size([LO, HI]) * 2 * nproc * ndims)

      if (ord_prolong_face_norm > max_opfn) ord_prolong_face_norm = max_opfn
      select case (ord_prolong_face_norm)
         case (0)
            opfn_c_ff(:) = [ half, zero, zero ]
            opfn_c_cf(:) = [ one,  zero, zero ]
         case (1)
#ifdef DEBUG
            ! the maximum convergence is at aux_R(1) = 0.25 +/- 0.05
            opfn_c_ff(:) = half * [ one + aux_R(1), -aux_R(1), zero ]
            ! the maximum convergence is at aux_R(2) = -0.05 +/- 0.05 and 0
            opfn_c_cf(:) = [ one + two*aux_R(2), -aux_R(2), zero ]
#else /* !DEBUG */
            opfn_c_ff(:) = [  5., -one, zero ] / 8.  ! adjusted experimentally
            opfn_c_cf(:) = [ 18.,  one, zero ] / 20. ! adjusted experimentally
#endif /* !DEBUG */
         case (2)
#ifdef DEBUG
            ! the maximum convergence is at aux_R(1) = 0.35 +/- 0.1 and aux_R(3) = 0.20 +- 0.05
            opfn_c_ff(:) = half * [ one+aux_R(1), -two*aux_R(1)+aux_R(3), aux_R(1)-aux_R(3) ]
            ! the maximum convergence is at aux_R(2) = -0.01 +/- 0.02 and aux_R(4) = -0.125 +- 0.025
            opfn_c_cf(:) = [ one + two*aux_R(2), -two*aux_R(2)+aux_R(4), aux_R(2)-aux_R(4) ]
#else /* !DEBUG */
            opfn_c_ff(:) = [ 27., -10., 3.  ] / 40. ! adjusted experimentally
            opfn_c_cf(:) = [  8., -one, one ] / 8.
#endif /* !DEBUG */
         case default
            call die("mg:mmpg opfn_c_[cf]f(:)")
      end select

      curl => base
      do while (associated(curl))

         if (ubound(curl%pse(proc)%sel(:,:,:), dim=1) > 1) call die("[multigrid_gravity:mpi_multigrid_prep_grav] Multiple blocks per process not implemented yet")

         per(:) = 0
         where (dom%periodic(:)) per(:) = curl%n_d(:)

         do d = xdim, zdim
            dmask(:) = .false.
            dmask(d) = .true.
            do lh = LO, HI
               if (dom%has_dir(d)) then
                  hl = LO+HI-lh

                  ! find coarse target for receiving data to be prolonged
                  if (associated(curl%coarser)) then
                     cgl => curl%first
                     do while (associated(cgl))
                        cg => cgl%cg
                        if (.not. cg%ext_bnd(d, lh)) then
                           ijks(:) = cg%ijkse(:, LO) - cg%off(:)  ! add this to convert absolute cell coordinates to local indices. (+nb - off(:))
                           procmask(:) = 0
                           ! two layers of cells are required for even locations (+ two layers per each interpolation order)
                           ! one layer of cells is required for odd locations (the local domain face is exactly at the centers of coarse cells)
                           coarsened(:, :) = cg%my_se(:, :)/2
                           coarsened(d, hl) = coarsened(d, lh)
                           select case (lh)
                              case (LO)
                                 if (mod(cg%off(d),    2_LONG) == 0) coarsened(d, :) = coarsened(d, :) + [ -1-ord_prolong_face_norm,   ord_prolong_face_norm ]
                              case (HI)
                                 if (mod(cg%h_cor1(d), 2_LONG) == 0) coarsened(d, :) = coarsened(d, :) + [   -ord_prolong_face_norm, 1+ord_prolong_face_norm ]
                           end select

                           do j = FIRST, LAST
                              if (ubound(curl%coarser%pse(j)%sel, dim=1) > 0) then
                                 if (is_overlap(coarsened(:,:), curl%coarser%pse(j)%sel(1, :, :), per)) procmask(j) = 1
                              endif
                           enddo
                           allocate(cg%mg%pfc_tgt(d, lh)%seg(count(procmask(:) /= 0)))

                           g = 0
                           do j = FIRST, LAST
                              if (procmask(j) /= 0) then
                                 g = g + 1
                                 if (.not. allocated(cg%mg%pfc_tgt(d, lh)%seg) .or. g>ubound(cg%mg%pfc_tgt(d, lh)%seg, dim=1)) call die("mg:mmpg pfc_tgt g>")
                                 seg => cg%mg%pfc_tgt(d, lh)%seg(g)
                                 if (allocated(seg%buf)) then
                                    call warn("mg:mmpg o seg%buf a a")
                                    deallocate(seg%buf)
                                 endif
                                 seg%proc = j
                                 ! find cross-section of own face segment with refined coarse segment
                                 b_layer(:, :) = cg%my_se(:, :)
                                 b_layer(d, hl) = b_layer(d, lh)
                                 !b_layer(d, lh) = b_layer(d, lh) + 2*lh-LO-HI ! extend to two layers of buffer

                                 where (.not. dmask(:)) ! find extents perpendicular to d
                                    seg%se(:, LO) = max(b_layer(:, LO), curl%coarser%pse(j)%sel(1, :, LO)*2  )
                                    seg%se(:, HI) = min(b_layer(:, HI), curl%coarser%pse(j)%sel(1, :, HI)*2+1)
                                 endwhere
                                 seg%se(d, :) = b_layer(d, :)
                                 allocate(seg%buf(seg%se(xdim, HI)/2-seg%se(xdim, LO)/2 + 1, &
                                      &           seg%se(ydim, HI)/2-seg%se(ydim, LO)/2 + 1, &
                                      &           seg%se(zdim, HI)/2-seg%se(zdim, LO)/2 + 1))

                                 seg%se(:, LO) = seg%se(:, LO) - cg%off(:) !+ ijks(:)
                                 seg%se(:, HI) = seg%se(:, HI) - cg%off(:) !+ ijks(:)
                              endif
                           enddo
                        endif

                        cgl => cgl%nxt
                     enddo
                  endif

                  ! find fine target(s) for sending the data to be prolonged
                  if (associated(curl%finer)) then
                     cgl => curl%first
                     do while (associated(cgl))
                        cg => cgl%cg
                        ijks(:) = cg%ijkse(:, LO) - cg%off(:)  ! add this to convert absolute cell coordinates to local indices. (+nb - off(:))
                        procmask(:) = 0
                        do j = FIRST, LAST
                           is_internal_fine = dom%periodic(d)
                           coarsened(:, :) = curl%finer%pse(j)%sel(1, :, :)/2
                           coarsened(d, hl) = coarsened(d, lh)
                           select case (lh)
                              case (LO)
                                 if (mod(curl%finer%pse(j)%sel(1, d, LO),     2_LONG) == 0) coarsened(d, :) = coarsened(d, :) + [ -1-ord_prolong_face_norm,   ord_prolong_face_norm ]
                                 is_internal_fine = is_internal_fine .or. (curl%finer%pse(j)%sel(1, d, lh) /= 0)
                              case (HI)
                                 if (mod(curl%finer%pse(j)%sel(1, d, HI) + 1, 2_LONG) == 0) coarsened(d, :) = coarsened(d, :) + [   -ord_prolong_face_norm, 1+ord_prolong_face_norm ]
                                 is_internal_fine = is_internal_fine .or. (curl%finer%pse(j)%sel(1, d, lh) + 1 < curl%finer%n_d(d))
                           end select
                           if (is_internal_fine) then
                              if (is_overlap(coarsened(:, :), cg%my_se(:, :), per)) procmask(j) = 1
                           endif
                        enddo
                        allocate(cg%mg%pff_tgt(d, lh)%seg(count(procmask(:) /= 0)))

                        g = 0
                        do j = FIRST, LAST
                           if (procmask(j) /= 0) then
                              g = g + 1
                              if (.not. allocated(cg%mg%pff_tgt(d, lh)%seg) .or. g>ubound(cg%mg%pff_tgt(d, lh)%seg, dim=1)) call die("mg:mmpg pff_tgt g>")
                              seg => cg%mg%pff_tgt(d, lh)%seg(g)
                              if (allocated(seg%buf)) then
                                 call warn("mg:mmpg o seg%buf a a")
                                 deallocate(seg%buf)
                              endif
                              seg%proc = j

                              ! find cross-section of own segment with coarsened fine face segment
                              coarsened(:, :) = curl%finer%pse(j)%sel(1, :, :)
                              coarsened(d, hl) = coarsened(d, lh)
                              coarsened(:, :) = coarsened(:, :)/2
                              where (.not. dmask(:))
                                 seg%se(:, LO) = max(cg%my_se(:, LO), coarsened(:, LO))
                                 seg%se(:, HI) = min(cg%my_se(:, HI), coarsened(:, HI))
                              endwhere
                              seg%se(d, :) = coarsened(d, :)
                              allocate(seg%buf(seg%se(xdim, HI)-seg%se(xdim, LO) + 1, &
                                   &           seg%se(ydim, HI)-seg%se(ydim, LO) + 1, &
                                   &           seg%se(zdim, HI)-seg%se(zdim, LO) + 1))

                              coarsened(:, :) = curl%finer%pse(j)%sel(1, :, :)
                              coarsened(d, hl) = coarsened(d, lh)
                              coarsened(d, lh) = coarsened(d, lh) + 2*lh-LO-HI ! extend to two layers of buffer
                              coarsened(:, :) = coarsened(:, :)/2
                              coarsened(d, :) = coarsened(d, :) + [ -ord_prolong_face_norm, ord_prolong_face_norm ]

                              seg%se(:, LO) = max(cg%my_se(:, LO), coarsened(:, LO)) + ijks(:)
                              seg%se(:, HI) = min(cg%my_se(:, HI), coarsened(:, HI)) + ijks(:)

                              coarsened(d, :) = coarsened(d, :) - [ -ord_prolong_face_norm, ord_prolong_face_norm ] ! revert broadening
                              allocate(seg%f_lay(seg%se(d, HI) - seg%se(d, LO) + 1))
                              do l = 1, size(seg%f_lay(:))
                                 seg%f_lay(l)%layer = l + int(seg%se(d, LO), kind=4) - 1
                                 nl = int(minval(abs(seg%f_lay(l)%layer - ijks(d) - coarsened(d, :))), kind=4)
                                 if (mod(curl%finer%pse(j)%sel(1, d, lh) + lh - LO, 2_LONG) == 0) then ! fine face at coarse face
                                    seg%f_lay(l)%coeff = opfn_c_ff(nl)
                                 else                                                              ! fine face at coarse center
                                    seg%f_lay(l)%coeff = opfn_c_cf(nl)
                                 endif
                              enddo

                           endif
                        enddo

                        cgl => cgl%nxt
                     enddo
                  endif

               endif
            enddo
         enddo

         curl => curl%finer
      enddo

   end subroutine mpi_multigrid_prep_grav

!!$ ============================================================================
!>
!! \brief Cleanup
!<

   subroutine cleanup_multigrid_grav

      use constants,     only: LO, HI, ndims
      use gc_list,       only: cg_list_element
      use cg_list_lev,   only: cg_list_level
      use grid_cont,     only: tgt_list
      use multipole,     only: cleanup_multipole
      use multigridvars, only: base

      implicit none

      integer :: g, ib
      integer, parameter :: nseg = 2*(HI-LO+1)*ndims
      type(tgt_list), dimension(nseg) :: io_tgt
      type(cg_list_level), pointer :: curl
      type(cg_list_element), pointer :: cgl

      call cleanup_multipole

      if (allocated(vstat%factor)) deallocate(vstat%factor)
      if (allocated(vstat%time)) deallocate(vstat%time)

      curl => base
      do while (associated(curl))
         cgl => curl%first
         do while (associated(cgl))
            if (allocated(cgl%cg%mg%fft))     deallocate(cgl%cg%mg%fft)
            if (allocated(cgl%cg%mg%fftr))    deallocate(cgl%cg%mg%fftr)
            if (allocated(cgl%cg%mg%src))     deallocate(cgl%cg%mg%src)
            if (allocated(cgl%cg%mg%Green3D)) deallocate(cgl%cg%mg%Green3D)

            if (cgl%cg%mg%planf /= 0) call dfftw_destroy_plan(cgl%cg%mg%planf)
            if (cgl%cg%mg%plani /= 0) call dfftw_destroy_plan(cgl%cg%mg%plani)

            io_tgt(1:nseg) = [ cgl%cg%mg%pfc_tgt, cgl%cg%mg%pff_tgt ]
            do ib = 1, nseg
               if (allocated(io_tgt(ib)%seg)) then
                  do g = lbound(io_tgt(ib)%seg, dim=1), ubound(io_tgt(ib)%seg, dim=1)
                     if (allocated(io_tgt(ib)%seg(g)%buf)) deallocate(io_tgt(ib)%seg(g)%buf)
                     if (allocated(io_tgt(ib)%seg(g)%f_lay)) deallocate(io_tgt(ib)%seg(g)%f_lay)
                  enddo
                  deallocate(io_tgt(ib)%seg)
               endif
            enddo
            cgl => cgl%nxt
         enddo
         curl => curl%finer
      enddo

      call dfftw_cleanup

   end subroutine cleanup_multigrid_grav

!!$ ============================================================================
!>
!! \brief This routine tries to construct first guess of potential based on previously obtained solution, if any.
!! If for some reason it is not desired to do the extrapolation (i.e. timestep varies too much) it would be good to have more control on this behavior.
!!
!! \details Quadratic extrapolation in time often gives better guess for smooth potentials, but is more risky for sharp-peaked potential fields (like moving self-bound clumps)
!! Rational extrapolation (1 + a t)(b + c t) can give 2 times better and 10 times worse guess depending on timestep.
!!
!! Set history%valid to .false. to force start from scratch.
!<

   subroutine init_solution(history)

#if defined(__INTEL_COMPILER)
      use cg_list_lev,      only: cg_list_level  ! QA_WARN workaround for stupid INTEL compiler
#endif /* __INTEL_COMPILER */
      use dataio_pub,       only: msg, die, printinfo
      use gc_list,          only: ind_val, all_cg
      use global,           only: t
      use grid,             only: leaves
      use mpisetup,         only: master
      use multigridvars,    only: roof, stdout, solution
      use multigridhelpers, only: set_dirty, check_dirty

      implicit none

      type(soln_history), intent(inout) :: history !< inner or outer potential history for recycling

      integer :: p0, p1, p2, ordt
      real, dimension(3) :: dt_fac
      enum, bind(C)
         enumerator :: RESET = -1, COPY = 0, LINEAR, QUADRATIC !! order of temporal interpolation
      end enum

      call set_dirty(solution)

      p0 = history%last
      if (nold > 0) then
         p1 = 1 + mod(p0 + nold - 2, nold) ! index of previous save
         p2 = 1 + mod(p1 + nold - 2, nold)
      else
         p1 = p0
         p2 = p0
      endif

      ! selfgrav_clump/initproblem.F90 requires monotonic time sequence t > history%old(p0)%time > history%old(p1)%time > history%old(p2)%time
      ordt = ord_time_extrap
      if (history%valid) then
         if ( history%old(p2)%time < history%old(p1)%time .and. &        ! quadratic interpolation
              history%old(p1)%time < history%old(p0)%time .and. &
              history%old(p0)%time < t) then
            ordt = min(int(QUADRATIC), ord_time_extrap)
         else if (history%old(p1)%time < history%old(p0)%time .and. &
              &   history%old(p0)%time < t) then      ! linear extrapolation
            ordt = min(int(LINEAR), ord_time_extrap)
         else                                                            ! simple recycling
            ordt = min(int(COPY), ord_time_extrap)
         endif
      else                                                               ! coldstart
         ordt = min(int(RESET), ord_time_extrap)
      endif

      select case (ordt)
         case (:RESET)
            if (master .and. ord_time_extrap > ordt) then
               write(msg, '(3a)')"[multigrid_gravity:init_solution] Clearing ",trim(vstat%cprefix),"solution."
               call printinfo(msg, stdout)
            endif
            call all_cg%set_q_value(solution, 0.)
            history%old(:)%time = -huge(1.0)
         case (COPY)
            call leaves%q_copy(history%old(p0)%i_hist, solution)
            if (master .and. ord_time_extrap > ordt) then
               write(msg, '(3a)')"[multigrid_gravity:init_solution] No extrapolation of ",trim(vstat%cprefix),"solution."
               call printinfo(msg, stdout)
            endif
         case (LINEAR)
            dt_fac(1) = (t - history%old(p0)%time) / (history%old(p0)%time - history%old(p1)%time)
            call leaves%q_lin_comb( [ ind_val(history%old(p0)%i_hist, (1.+dt_fac(1))), &
                 &                    ind_val(history%old(p1)%i_hist,    -dt_fac(1) ) ], solution )
            if (master .and. ord_time_extrap > ordt) then
               write(msg, '(3a)')"[multigrid_gravity:init_solution] Linear extrapolation of ",trim(vstat%cprefix),"solution."
               call printinfo(msg, stdout)
            endif
         case (QUADRATIC)
            dt_fac(1) = (t - history%old(p0)%time) / (history%old(p1)%time - history%old(p2)%time)
            dt_fac(2) = (t - history%old(p1)%time) / (history%old(p2)%time - history%old(p0)%time)
            dt_fac(3) = (t - history%old(p2)%time) / (history%old(p0)%time - history%old(p1)%time)
            call leaves%q_lin_comb([ ind_val(history%old(p0)%i_hist, -dt_fac(2)*dt_fac(3)), &
                 &                   ind_val(history%old(p1)%i_hist, -dt_fac(1)*dt_fac(3)), &
                 &                   ind_val(history%old(p2)%i_hist, -dt_fac(1)*dt_fac(2)) ], solution )
         case default
            call die("[multigrid_gravity:init_solution] Extrapolation order not implemented")
      end select

      call check_dirty(roof, solution, "init_soln")

   end subroutine init_solution

!!$ ============================================================================
!>
!! \brief Make a local copy of source (density) and multiply by 4 pi G
!!
!! \details Typically i_all_dens is a copy of fluidindex::iarr_all_sg.
!! Passing this as an argument allows for independent computation of the potential for several density fields if necessary.
!<

   subroutine init_source(i_all_dens)

#if defined(__INTEL_COMPILER)
      use cg_list_lev,        only: cg_list_level    ! QA_WARN workaround for stupid INTEL compiler
#endif /* __INTEL_COMPILER */
      use constants,          only: GEO_RPZ, LO, HI, xdim, ydim, zdim
      use dataio_pub,         only: die
      use domain,             only: dom
      use gc_list,            only: cg_list_element
      use grid,               only: leaves
      use grid_cont,          only: grid_container
      use multigridhelpers,   only: set_dirty, check_dirty
      use multigridvars,      only: roof, source, bnd_periodic, bnd_dirichlet, bnd_givenval
      use units,              only: fpiG
#ifdef JEANS_PROBLEM
      use problem_pub,        only: jeans_d0, jeans_mode ! hack for tests
#endif /* JEANS_PROBLEM */

      implicit none

      integer(kind=4), dimension(:), optional, intent(in) :: i_all_dens !< indices to selfgravitating fluids

      real :: fac
      integer :: i
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg

      call set_dirty(source)

      if (present(i_all_dens)) then
         cgl => leaves%first
         do while (associated(cgl))
            cg => cgl%cg
            cg%q(source)%arr(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) = fpiG * sum(cg%u(i_all_dens, cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke), dim=1)
            cgl => cgl%nxt
         enddo
      else
         if (grav_bnd /= bnd_givenval) call die("[multigrid_gravity:init_source] empty space allowed only for given value boundaries.")
         call leaves%set_q_value(source, 0.)
      endif

      cg => leaves%first%cg
      select case (grav_bnd)
         case (bnd_periodic) ! probably also bnd_neumann
            call leaves%subtract_average(source)
         case (bnd_dirichlet)
#ifdef JEANS_PROBLEM
            if (jeans_mode == 1) call leaves%q_add_val(source, - fpiG * jeans_d0) ! remove density bias
#endif /* JEANS_PROBLEM */
         case (bnd_givenval) ! convert potential into a layer of imaginary mass (subtract second derivative normal to computational domain boundary)

            cgl => leaves%first
            do while (associated(cgl))
               cg => cgl%cg
               if (cg%ext_bnd(xdim, LO)) then
                  fac = 2. * cg%idx2 / fpiG
                  if (dom%geometry_type == GEO_RPZ .and. cg%x(cg%is) /= 0.) fac = fac - 1./(cg%dx * cg%x(cg%is) * fpiG) !> BEWARE is it cg%x(ie), cg%x(ie+1) or something in the middle?
                  cg%q(source)%arr       (cg%is, cg%js:cg%je, cg%ks:cg%ke) = &
                       & cg%q(source)%arr(cg%is, cg%js:cg%je, cg%ks:cg%ke) - &
                       & cg%mg%bnd_x(            cg%js:cg%je, cg%ks:cg%ke, LO) * fac
               endif
               if (cg%ext_bnd(xdim, HI)) then
                  fac = 2. * cg%idx2 / fpiG
                  if (dom%geometry_type == GEO_RPZ .and. cg%x(cg%ie) /= 0.) fac = fac - 1./(cg%dx * cg%x(cg%ie) * fpiG) !> BEWARE is it cg%x(ie), cg%x(ie+1) or something in the middle?
                  cg%q(source)%arr       (cg%ie, cg%js:cg%je, cg%ks:cg%ke) = &
                       & cg%q(source)%arr(cg%ie, cg%js:cg%je, cg%ks:cg%ke) - &
                       & cg%mg%bnd_x(            cg%js:cg%je, cg%ks:cg%ke, HI) * fac
               endif
               if (cg%ext_bnd(ydim, LO)) then
                  if (dom%geometry_type == GEO_RPZ) then
                     do i = cg%is, cg%ie
                        if (cg%x(i) /= 0.) then
                           cg%q(source)%arr       (i, cg%js, cg%ks:cg%ke) = &
                                & cg%q(source)%arr(i, cg%js, cg%ks:cg%ke) - &
                                & cg%mg%bnd_y(i,             cg%ks:cg%ke, LO) * 2. * cg%idy2 / fpiG / cg%x(i)**2
                        endif
                     enddo
                  else
                     cg%q(source)%arr       (cg%is:cg%ie, cg%js, cg%ks:cg%ke) = &
                          & cg%q(source)%arr(cg%is:cg%ie, cg%js, cg%ks:cg%ke) - &
                          & cg%mg%bnd_y     (cg%is:cg%ie,        cg%ks:cg%ke, LO) * 2. * cg%idy2 / fpiG
                  endif
               endif
               if (cg%ext_bnd(ydim, HI)) then
                  if (dom%geometry_type == GEO_RPZ) then
                     do i = cg%is, cg%ie
                        if (cg%x(i) /= 0.) then
                           cg%q(source)%arr       (i, cg%je, cg%ks:cg%ke) = &
                                & cg%q(source)%arr(i, cg%je, cg%ks:cg%ke) - &
                                & cg%mg%bnd_y     (i,        cg%ks:cg%ke, HI) * 2. * cg%idy2 / fpiG / cg%x(i)**2
                        endif
                     enddo
                  else
                     cg%q(source)%arr       (cg%is:cg%ie, cg%je, cg%ks:cg%ke) = &
                          & cg%q(source)%arr(cg%is:cg%ie, cg%je, cg%ks:cg%ke) - &
                          & cg%mg%bnd_y     (cg%is:cg%ie,        cg%ks:cg%ke, HI) * 2. * cg%idy2 / fpiG
                  endif
               endif
               if (cg%ext_bnd(zdim, LO)) cg%q(source)%arr(cg%is:cg%ie, cg%js:cg%je, cg%ks) = &
                    &                    cg%q(source)%arr(cg%is:cg%ie, cg%js:cg%je, cg%ks) - &
                    &                    cg%mg%bnd_z     (cg%is:cg%ie, cg%js:cg%je,          LO) * 2. * cg%idz2 / fpiG
               if (cg%ext_bnd(zdim, HI)) cg%q(source)%arr(cg%is:cg%ie, cg%js:cg%je, cg%ke) = &
                    &                    cg%q(source)%arr(cg%is:cg%ie, cg%js:cg%je, cg%ke) - &
                    &                    cg%mg%bnd_z     (cg%is:cg%ie, cg%js:cg%je,          HI) * 2. * cg%idz2 / fpiG
               !> \todo compactify the above mess
               cgl => cgl%nxt
            enddo

         case default
            call die("[multigrid_gravity:init_source] Unknown boundary type")
      end select

      call check_dirty(roof, source, "init_src")

   end subroutine init_source

!!$ ============================================================================
!>
!! \brief This routine manages old copies of potential for recycling.
!<

   subroutine store_solution(history)

#if defined(__INTEL_COMPILER)
      use cg_list_lev,       only: cg_list_level   ! QA_WARN workaround for stupid INTEL compiler
#endif /* __INTEL_COMPILER */
      use domain,            only: dom
      use global,            only: t
      use grid,              only: leaves
      use multigridvars,     only: roof, bnd_isolated, bnd_givenval, solution, extbnd_extrapolate, extbnd_mirror
      use multigridmpifuncs, only: mpi_multigrid_bnd

      implicit none

      type(soln_history), intent(inout) :: history !< inner or outer potential history to store recent solution

      if (nold <= 0) return

      if (history%valid) then
         history%last = 1 + mod(history%last, nold)
      else
         history%old(:)%time = t ! prevents extrapolation too early
      endif

      call leaves%q_copy(solution, history%old(history%last)%i_hist)
      history%old(history%last)%time = t
      history%valid = .true.

      ! Update guardcells of the solution before leaving. This can be done in higher-level routines that collect all the gravity contributions, but would be less safe.
      ! Extrapolate isolated boundaries, remember that grav_bnd is messed up by multigrid_solve_*
      if (grav_bnd == bnd_isolated .or. grav_bnd == bnd_givenval) then
         call mpi_multigrid_bnd(roof, solution, dom%nb, extbnd_extrapolate)
      else
         call mpi_multigrid_bnd(roof, solution, dom%nb, extbnd_mirror)
      endif

   end subroutine store_solution

!!$ ============================================================================
!>
!! \brief Multigrid gravity driver. This is the only multigrid routine intended to be called from the gravity module.
!! This routine is also responsible for communicating the solution to the rest of world via sgp array.
!<

   subroutine multigrid_solve_grav(i_all_dens)

      use constants,     only: sgp_n
      use gc_list,       only: all_cg
      use grid,          only: leaves
      use multigridvars, only: solution, bnd_isolated, bnd_dirichlet, bnd_givenval, tot_ts, ts
      use multipole,     only: multipole_solver
      use timer,         only: set_timer

      implicit none

      integer(kind=4), dimension(:), intent(in) :: i_all_dens !< indices to selfgravitating fluids

      logical :: isolated

      ts =  set_timer("multigrid", .true.)

      isolated = (grav_bnd == bnd_isolated) !> \deprecated BEWARE: not elegant; probably there should be two global grav_bnd variables

      if (isolated) then
         grav_bnd = bnd_dirichlet
         vstat%cprefix = "Gi-"
      else
#ifdef COSM_RAYS
         vstat%cprefix = "G-"
#else /* !COSM_RAYS */
         vstat%cprefix = ""
#endif /* !COSM_RAYS */
      endif

      call init_source(i_all_dens)

      call vcycle_hg(inner)

      call leaves%q_copy(solution, all_cg%ind(sgp_n))

      if (isolated) then
         grav_bnd = bnd_givenval

         vstat%cprefix = "Go-"
         call multipole_solver
         call init_source

         call vcycle_hg(outer)

         call leaves%q_add(solution, all_cg%ind(sgp_n)) ! add solution to sgp

         grav_bnd = bnd_isolated ! restore

      endif

      ts = set_timer("multigrid")
      tot_ts = tot_ts + ts

   end subroutine multigrid_solve_grav

!!$ ============================================================================
!>
!! \brief The solver. Here we choose an adaptation of the Huang-Greengard V-cycle.
!! For more difficult problems, like variable coefficient diffusion equation a more sophisticated V-cycle may be more effective.
!<

   subroutine vcycle_hg(history)

      use cg_list_lev,        only: cg_list_level
      use constants,          only: cbuff_len
      use dataio_pub,         only: msg, die, warn, printinfo
      use grid,               only: leaves
      use mpisetup,           only: master, nproc
      use multigridhelpers,   only: set_dirty, check_dirty, do_ascii_dump, numbered_ascii_dump
      use multigridvars,      only: roof, base, source, solution, correction, defect, verbose_vcycle, bnd_periodic, stdout, tot_ts, ts
      use timer,              only: set_timer

      implicit none

      type(soln_history), intent(inout) :: history !< inner or outer potential history used for initializing first guess

      real,    parameter :: suspicious_factor = 1.05 !> \deprecated If the norm decreases too slowly then dump diagnostic output (BEWARE: this option is for tests only)
      integer            :: v
      real               :: norm_rhs, norm_lhs, norm_old, norm_lowest
      logical            :: dump_every_step, dump_result
      logical, save      :: norm_was_zero = .false.
      integer, parameter       :: fmtlen = 32
      character(len=fmtlen)    :: fmt
      character(len=cbuff_len) :: dname
      type(cg_list_level), pointer :: curl

      inquire(file = "_dump_every_step_", EXIST=dump_every_step) ! use for debug only
      inquire(file = "_dump_result_", EXIST=dump_result)

      do_ascii_dump = do_ascii_dump .or. dump_every_step .or. dump_result

      ! On single CPU use FFT if possible because it is faster. Can be disabled by prefer_rbgs_relaxation = .true.
      if (nproc == 1 .and. roof%first%cg%mg%fft_type /= fft_none) then
         call set_dirty(solution)
         call fft_solve_roof
         if (trust_fft_solution) then
            write(msg, '(3a)')"[multigrid_gravity:vcycle_hg] FFT solution trusted, skipping ", trim(vstat%cprefix), "cycle."
            call printinfo(msg, stdout)
            return
         endif
      else
         call init_solution(history)
      endif

      norm_lhs = 0.
      norm_rhs = leaves%norm_sq(source)
      norm_old = norm_rhs
      norm_lowest = norm_rhs

      if (norm_rhs == 0.) then ! empty domain => potential == 0.
         if (master .and. .not. norm_was_zero) call warn("[multigrid_gravity:vcycle_hg] No gravitational potential for an empty space.")
         norm_was_zero = .true.
         call store_solution(history)
         return
      else
         if (master .and. norm_was_zero) call warn("[multigrid_gravity:vcycle_hg] Spontaneous mass creation detected!")
         norm_was_zero = .false.
      endif

!      if (.not. history%valid .and. prefer_rbgs_relaxation) call approximate_solution(roof, source, solution) ! not necessary when init_solution called FFT
! difficult statement: for approximate_solution_fft it requires to pass a flag to use guardcells instead of prolonging faces.
! how much does it improve? (make a benchmark at some point)

      ! iterations
      do v = 0, max_cycles

         call set_dirty(defect)
         call residual(roof, source, solution, defect) !leaves
         if (grav_bnd == bnd_periodic) call leaves%subtract_average(defect)
         call check_dirty(roof, defect, "residual") !leaves

         norm_lhs = leaves%norm_sq(defect)
         ts = set_timer("multigrid")
         tot_ts = tot_ts + ts
         if (master .and. verbose_vcycle) then
            if (norm_old/norm_lhs < 1.e5) then
               fmt='(3a,i3,a,f12.9,a,f8.2,a,f7.3)'
            else
               fmt='(3a,i3,a,f12.9,a,es8.2,a,f7.3)'
            endif
            write(msg, fmt)"[multigrid_gravity] ", trim(vstat%cprefix), "Cycle:", v, " norm/rhs= ", norm_lhs/norm_rhs, " reduction factor= ", norm_old/norm_lhs, "   dt_wall= ", ts
            call printinfo(msg, stdout)
         endif

         vstat%count = v
         if (norm_lhs /= 0) then
            vstat%factor(vstat%count) = norm_old/norm_lhs
         else
            vstat%factor(vstat%count) = huge(1.0)
         endif
         vstat%time(vstat%count) = ts

         if (norm_old/norm_lhs <= suspicious_factor .or. dump_every_step .or. (norm_lhs/norm_rhs <= norm_tol .and. dump_result)) then
            write(dname,'(2a)')trim(vstat%cprefix),"mdump"
            if (dump_result .and. norm_lhs/norm_rhs <= norm_tol) then
               call numbered_ascii_dump(dname)
            else
               call numbered_ascii_dump(dname, v)
            endif
         endif

         if (norm_lhs/norm_rhs <= norm_tol) exit

         if (v<1) then ! forgive poor convergence in some first V-cycles
            norm_lowest = norm_lhs
         else
            if (norm_lhs < norm_lowest) then
               norm_lowest = norm_lhs
            else
               if (norm_lhs/norm_lowest > vcycle_abort) then
                  vstat%norm_final = norm_lhs/norm_rhs
                  if (.not. verbose_vcycle) call vstat%brief_v_log
                  call die("[multigrid_gravity:vcycle_hg] Serious nonconvergence detected.")
                  !In such case one may increase nsmool, decrease refinement depth or use FFT
               endif
            endif
         endif

         norm_old = norm_lhs

         ! the Huang-Greengard V-cycle
         call roof%restrict_to_floor_q_1var(defect)

         ! call set_dirty(correction)
         call base%set_q_value(correction, 0.)

         curl => base
         do while (associated(curl))
            call approximate_solution(curl, defect, correction)
            call check_dirty(curl, correction, "Vup relax+")
            curl => curl%finer
         enddo
         call leaves%q_add(correction, solution)

      enddo

      if (v > max_cycles) then
         if (master .and. norm_lhs/norm_rhs > norm_tol) call warn("[multigrid_gravity:vcycle_hg] Not enough V-cycles to achieve convergence.")
         v = max_cycles
      endif

      vstat%norm_final = norm_lhs/norm_rhs
      if (.not. verbose_vcycle) call vstat%brief_v_log

      call check_dirty(roof, solution, "final_solution")

      call store_solution(history)

   end subroutine vcycle_hg

!!$ ============================================================================
!>
!! \brief Calculate the residuum for the Poisson equation.
!<

   subroutine residual(curl, src, soln, def)

      use dataio_pub,  only: die
      use cg_list_lev, only: cg_list_level

      implicit none

      type(cg_list_level), pointer, intent(in) :: curl !< pointer to a level for which we approximate the solution
      integer,                      intent(in) :: src  !< index of source in cg%q(:)
      integer,                      intent(in) :: soln !< index of solution in cg%q(:)
      integer,                      intent(in) :: def  !< index of defect in cg%q(:)

      select case (ord_laplacian)
      case (2)
         call residual2(curl, src, soln, def)
      case (4)
         call residual4(curl, src, soln, def)
      case default
         call die("[multigrid_gravity:residual] The parameter 'ord_laplacian' must be 2 or 4")
      end select

   end subroutine residual

!!$ ============================================================================
!>
!! \brief 2nd order Laplacian
!<

   subroutine residual2(curl, src, soln, def)

      use constants,         only: xdim, ydim, zdim, ndims, GEO_XYZ, GEO_RPZ, zero, half, I_ONE
      use dataio_pub,        only: die
      use domain,            only: dom
      use gc_list,           only: cg_list_element
      use cg_list_lev,       only: cg_list_level
      use grid_cont,         only: grid_container
      use multigridvars,     only: extbnd_antimirror
      use multigridmpifuncs, only: mpi_multigrid_bnd

      implicit none

      type(cg_list_level), pointer, intent(in) :: curl !< pointer to a level for which we approximate the solution
      integer,                      intent(in) :: src  !< index of source in cg%q(:)
      integer,                      intent(in) :: soln !< index of solution in cg%q(:)
      integer,                      intent(in) :: def  !< index of defect in cg%q(:)

      real    :: L0, Lx, Ly, Lz, Lx1
      integer :: i, j, k
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg

      call mpi_multigrid_bnd(curl, soln, I_ONE, extbnd_antimirror) ! no corners required

      ! Possible optimization candidate: reduce cache misses (secondary importance, cache-aware implementation required)
      ! Explicit loop over k gives here better performance than array operation due to less cache misses (at least on 32^3 and 64^3 arrays)
      cgl => curl%first
      do while (associated(cgl))
         cg => cgl%cg

         ! Coefficients for a simplest 3-point Laplacian operator: [ 1, -2, 1 ]
         ! for 2D and 1D setups appropriate elements of [ Lx, Ly, Lz ] should be == 0.
         Lx = cg%idx2
         Ly = cg%idy2
         Lz = cg%idz2
         L0 = -2. * (Lx + Ly + Lz)

         select case (dom%geometry_type)
            case (GEO_XYZ)
               if (dom%eff_dim == ndims .and. .not. multidim_code_3D) then
                  do k = cg%ks, cg%ke
                     cg%q(def)%arr        (cg%is  :cg%ie,   cg%js  :cg%je,   k)         = &
                          & cg%q(src)%arr (cg%is  :cg%ie,   cg%js  :cg%je,   k)         - &
                          ( cg%q(soln)%arr(cg%is-1:cg%ie-1, cg%js  :cg%je,   k)         + &
                          & cg%q(soln)%arr(cg%is+1:cg%ie+1, cg%js  :cg%je,   k))   * Lx - &
                          ( cg%q(soln)%arr(cg%is  :cg%ie,   cg%js-1:cg%je-1, k)         + &
                          & cg%q(soln)%arr(cg%is  :cg%ie,   cg%js+1:cg%je+1, k))   * Ly - &
                          ( cg%q(soln)%arr(cg%is  :cg%ie,   cg%js  :cg%je,   k-1)       + &
                          & cg%q(soln)%arr(cg%is  :cg%ie,   cg%js  :cg%je,   k+1)) * Lz - &
                          & cg%q(soln)%arr(cg%is  :cg%ie,   cg%js  :cg%je,   k)    * L0
                  enddo
               else
                  ! In 3D this implementation can give a bit more cache misses, few times more writes and significantly more instructions executed than monolithic 3D above
                  do k = cg%ks, cg%ke
                     cg%q(def)%arr        (cg%is  :cg%ie,   cg%js  :cg%je,   k)    = &
                          & cg%q(src)%arr (cg%is  :cg%ie,   cg%js  :cg%je,   k)    - &
                          & cg%q(soln)%arr(cg%is  :cg%ie,   cg%js  :cg%je,   k)    * L0
                     if (dom%has_dir(xdim)) &
                          & cg%q(def)%arr (cg%is  :cg%ie,   cg%js  :cg%je,   k)    = &
                          & cg%q(def)%arr (cg%is  :cg%ie,   cg%js  :cg%je,   k)    - &
                          ( cg%q(soln)%arr(cg%is-1:cg%ie-1, cg%js  :cg%je,   k)    + &
                          & cg%q(soln)%arr(cg%is+1:cg%ie+1, cg%js  :cg%je,   k))   * Lx
                     if (dom%has_dir(ydim)) &
                          & cg%q(def)%arr (cg%is  :cg%ie,   cg%js  :cg%je,   k)    = &
                          & cg%q(def)%arr (cg%is  :cg%ie,   cg%js  :cg%je,   k)    - &
                          ( cg%q(soln)%arr(cg%is  :cg%ie,   cg%js-1:cg%je-1, k)    + &
                          & cg%q(soln)%arr(cg%is  :cg%ie,   cg%js+1:cg%je+1, k))   * Ly
                     if (dom%has_dir(zdim)) &
                          & cg%q(def)%arr (cg%is  :cg%ie,   cg%js  :cg%je,   k)    = &
                          & cg%q(def)%arr (cg%is  :cg%ie,   cg%js  :cg%je,   k)    - &
                          ( cg%q(soln)%arr(cg%is  :cg%ie,   cg%js  :cg%je,   k-1)  + &
                          & cg%q(soln)%arr(cg%is  :cg%ie,   cg%js  :cg%je,   k+1)) * Lz
                  enddo
               endif
            case (GEO_RPZ)
               Lx = cg%idx2
               Lz = cg%idz2
               do k = cg%ks, cg%ke
                  do j = cg%js, cg%je
                     do i = cg%is, cg%ie
                        if (cg%x(i) /= zero) then !> \todo convert Ly, Lx1 and L0 into precomputed arrays
                           Ly = cg%idy2 / cg%x(i)**2 ! cylindrical factor
                           Lx1 = half / (cg%dx * cg%x(i))
                        else
                           Ly = zero
                           Lx1 = zero
                        endif
                        L0 = -2. * (Lx + Ly + Lz)
                        cg%q(def)%arr(i, j, k) = cg%q(src)%arr (i,   j,   k)   - cg%q(soln)%arr(i,   j,   k)    * L0
                        if (dom%has_dir(xdim))   cg%q(def)%arr (i,   j,   k)   = cg%q(def) %arr(i,   j,   k)         - &
                             &                  (cg%q(soln)%arr(i+1, j,   k)   + cg%q(soln)%arr(i-1, j,   k))   * Lx - &
                             &                  (cg%q(soln)%arr(i+1, j,   k)   - cg%q(soln)%arr(i-1, j,   k))   * Lx1    ! cylindrical term
                        if (dom%has_dir(ydim))   cg%q(def)%arr (i,   j,   k)   = cg%q(def) %arr(i,   j,   k)         - &
                             &                  (cg%q(soln)%arr(i,   j+1, k)   + cg%q(soln)%arr(i,   j-1, k))   * Ly
                        if (dom%has_dir(zdim))   cg%q(def)%arr (i,   j,   k)   = cg%q(def) %arr(i,   j,   k)         - &
                             &                  (cg%q(soln)%arr(i,   j,   k+1) + cg%q(soln)%arr(i,   j,   k-1)) * Lz
                     enddo
                  enddo
               enddo
            case default
               call die("[multigrid_gravity:residual2] Unsupported geometry.")
         end select
         cgl => cgl%nxt
      enddo

   end subroutine residual2

!!$ ============================================================================
!>
!! \brief 4th order Laplacian
!!
!! \details Significantly slows down convergence, does not seem to improve quality of solution in simple tests.
!!
!! L4 = [0, 1, -2, 1, 0] + L4_strength * 1./12. * [ -1, 4, -6, 4, -1 ]
!! For integrated face fluxes in the 4th order Laplacian estimate set L4_strength = 0.5
!! For simple 5-point L4 set L4_strength = 1.0
!!
!! There also exists more compact Mehrstellen scheme.
!<

   subroutine residual4(curl, src, soln, def)

      use constants,         only: I_TWO, ndims, idm2, xdim, ydim, zdim
      use dataio_pub,        only: die, warn
      use domain,            only: dom
      use gc_list,           only: cg_list_element
      use cg_list_lev,       only: cg_list_level
      use grid_cont,         only: grid_container
      use mpisetup,          only: master
      use multigridvars,     only: bnd_givenval, extbnd_antimirror
      use multigridmpifuncs, only: mpi_multigrid_bnd
      use named_array,       only: p3

      implicit none

      type(cg_list_level), pointer, intent(in) :: curl !< pointer to a level for which we approximate the solution
      integer,                      intent(in) :: src  !< index of source in cg%q(:)
      integer,                      intent(in) :: soln !< index of solution in cg%q(:)
      integer,                      intent(in) :: def  !< index of defect in cg%q(:)

      real, parameter     :: L4_scaling = 1./12. ! with L4_strength = 1. this gives an L4 approximation for finite differences approach
      integer, parameter  :: L2w = 2             ! #layers of boundary cells for L2 operator

      real                :: c21, c41, c42 !, c20, c40
      real                :: L0, Lx1, Lx2, Ly1, Ly2, Lz1, Lz2, Lx, Ly, Lz

      logical, save       :: firstcall = .true.
      integer             :: i, j, k
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg

      if (dom%eff_dim<ndims) call die("[multigrid_gravity:residual4] Only 3D is implemented")

      if (firstcall) then
         if (master) call warn("[multigrid_gravity:residual4] residual order 4 is experimental.")
         firstcall = .false.
      endif

      call mpi_multigrid_bnd(curl, soln, I_TWO, extbnd_antimirror) ! no corners required

      c21 = 1.
      c42 = - L4_scaling * L4_strength
      c41 = c21 + 4. * L4_scaling * L4_strength
      !c20 = -2.
      !c40 = c20 - 6. * L4_strength

      cgl => curl%first
      do while (associated(cgl))
         cg => cgl%cg

         Lx1 = c41 * cg%idx2
         Ly1 = c41 * cg%idy2
         Lz1 = c41 * cg%idz2
         Lx2 = c42 * cg%idx2
         Ly2 = c42 * cg%idy2
         Lz2 = c42 * cg%idz2
         ! L0  = c40 * (cg%idx2 + cg%idy2 + cg%idz2 )
         L0 = -2. * (Lx1 + Lx2 + Ly1 + Ly2 + Lz1 + Lz2)

         !> \deprecated BEWARE: cylindrical factors go here
         p3 => cg%q(def)%span(cg%ijkse)
         p3 = cg%q(src)%span(cg%ijkse) - &
              cg%q(soln)%span(cg%ijkse-2*idm2(xdim,:,:)) * Lx2 - cg%q(soln)%span(cg%ijkse+2*idm2(xdim,:,:)) * Lx2 - &
              cg%q(soln)%span(cg%ijkse-  idm2(xdim,:,:)) * Lx1 - cg%q(soln)%span(cg%ijkse+  idm2(xdim,:,:)) * Lx1 - &
              cg%q(soln)%span(cg%ijkse-2*idm2(ydim,:,:)) * Ly2 - cg%q(soln)%span(cg%ijkse+2*idm2(ydim,:,:)) * Ly2 - &
              cg%q(soln)%span(cg%ijkse-  idm2(ydim,:,:)) * Ly1 - cg%q(soln)%span(cg%ijkse+  idm2(ydim,:,:)) * Ly1 - &
              cg%q(soln)%span(cg%ijkse-2*idm2(zdim,:,:)) * Lz2 - cg%q(soln)%span(cg%ijkse+2*idm2(zdim,:,:)) * Lz2 - &
              cg%q(soln)%span(cg%ijkse-  idm2(zdim,:,:)) * Lz1 - cg%q(soln)%span(cg%ijkse+  idm2(zdim,:,:)) * Lz1 - cg%q(soln)%span(cg%ijkse)   * L0

         ! WARNING: not optimized
         if (grav_bnd == bnd_givenval) then ! probably also in some other cases
            ! Use L2 Laplacian in two layers of cells next to the boundary because L4 seems to be incompatible with present image mass construction
            Lx = c21 * cg%idx2
            Ly = c21 * cg%idy2
            Lz = c21 * cg%idz2
            L0 = -2. * (Lx + Ly + Lz)

            do k = cg%ks, cg%ke
               do j = cg%js, cg%je
                  do i = cg%is, cg%ie
                     if ( i<cg%is+L2w .or. i>cg%ie-L2w .or. j<cg%js+L2w .or. j>cg%je-L2w .or. k<cg%ks+L2w .or. k>cg%ke-L2w) then
                        cg%q(def)%arr        (i,   j,   k)   = cg%q(src)%arr (i,   j,   k)         - &
                             ( cg%q(soln)%arr(i-1, j,   k)   + cg%q(soln)%arr(i+1, j,   k))   * Lx - &
                             ( cg%q(soln)%arr(i,   j-1, k)   + cg%q(soln)%arr(i,   j+1, k))   * Ly - &
                             ( cg%q(soln)%arr(i,   j,   k-1) + cg%q(soln)%arr(i,   j,   k+1)) * Lz - &
                             & cg%q(soln)%arr(i,   j,   k)                                    * L0
                     endif
                  enddo
               enddo
            enddo
         endif
         cgl => cgl%nxt
      enddo

   end subroutine residual4

!!$ ============================================================================
!>
!! \brief This routine has to find an approximate solution for given source field and implemented differential operator
!<

   subroutine approximate_solution(curl, src, soln)

      use cg_list_lev,        only: cg_list_level
      use multigridhelpers,   only: check_dirty
      use multigridbasefuncs, only: prolong_level
      use multigridvars,      only: base, roof, correction

      implicit none

      type(cg_list_level), pointer, intent(in) :: curl !< pointer to a level for which we approximate the solution
      integer,                      intent(in) :: src  !< index of source in cg%q(:)
      integer,                      intent(in) :: soln !< index of solution in cg%q(:)

      call check_dirty(curl, src, "approx_soln src-")

      if (associated(curl%first)) then
         if (curl%first%cg%mg%fft_type /= fft_none) then
            call approximate_solution_fft(curl, src, soln)
         else
            call check_dirty(curl, soln, "approx_soln soln-")
            call approximate_solution_rbgs(curl, src, soln)
         endif
      endif

      if (prefer_rbgs_relaxation .and. soln == correction .and. .not. associated(curl, roof)) call prolong_level(curl, correction)
      !> \deprecated BEWARE other implementations of the multigrid algorithm may be incompatible with prolongation called from here

      call check_dirty(curl, soln, "approx_soln soln+")

   end subroutine approximate_solution

!!$ ============================================================================
!>
!! \brief Red-Black Gauss-Seidel relaxation.
!!
!! \details This is the most costly routine in a serial run. Try to find optimal values for nsmool and nsmoob.
!! This routine also depends a lot on communication so it  may limit scalability of the multigrid.
!! \todo Implement convergence check on base level (not very important since we have a FFT solver for base level)
!<

   subroutine approximate_solution_rbgs(curl, src, soln)

      use constants,         only: xdim, ydim, zdim, ndims, GEO_XYZ, GEO_RPZ, I_ONE
      use dataio_pub,        only: die
      use domain,            only: dom
      use gc_list,           only: cg_list_element
      use cg_list_lev,       only: cg_list_level
      use grid_cont,         only: grid_container
      use multigridvars,     only: base, extbnd_antimirror
      use multigridhelpers,  only: dirty_debug, check_dirty, dirty_label
      use multigridmpifuncs, only: mpi_multigrid_bnd

      implicit none

      type(cg_list_level), pointer, intent(in) :: curl !< pointer to a level for which we approximate the solution
      integer,                      intent(in) :: src  !< index of source in cg%q(:)
      integer,                      intent(in) :: soln !< index of solution in cg%q(:)

      integer, parameter :: RED_BLACK = 2 !< the checkerboard requires two sweeps

      integer :: n, i, j, k, i1, j1, k1, id, jd, kd
      integer :: nsmoo
      real    :: crx, crx1, cry, crz, cr
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg

      if (associated(curl, base)) then
         nsmoo = nsmoob
      else
         nsmoo = nsmool
      endif

      if (dom%geometry_type == GEO_RPZ .and. .not. multidim_code_3D) call die("[multigrid_gravity:approximate_solution_rbgs] multidim_code_3D = .false. not implemented")

      do n = 1, RED_BLACK*nsmoo
         call mpi_multigrid_bnd(curl, soln, I_ONE, extbnd_antimirror) ! no corners are required here

         if (dirty_debug) then
            write(dirty_label, '(a,i5)')"relax soln- smoo=", n
            call check_dirty(curl, soln, dirty_label)
         endif
         cgl => curl%first
         do while (associated(cgl))
            cg => cgl%cg

            ! Possible optimization: this is the most costly part of the RBGS relaxation (instruction count, read and write data, L1 and L2 read cache miss)
            ! do n = 1, nsmoo
            !    call mpi_multigrid_bnd(curl, soln, I_ONE, extbnd_antimirror)
            !    relax single layer of red cells at all faces
            !    call mpi_multigrid_bnd(curl, soln, I_ONE, extbnd_antimirror)
            !    relax interior cells (except for single layer of cells at all faces), first red, then 1-cell behind black one.
            !    relax single layer of black cells at all faces
            ! enddo

            ! with explicit outer loops it is easier to describe a 3-D checkerboard :-)

            if (dom%eff_dim==ndims .and. .not. multidim_code_3D) then
               do k = cg%ks, cg%ke
                  do j = cg%js, cg%je
                     i1 = cg%is + mod(n+j+k, RED_BLACK)
                     if (dom%geometry_type == GEO_RPZ) then
!!$                  cg%q(soln)%arr(i1  :cg%ie  :2, j,   k) = &
!!$                       cg%mg%rx * (cg%q(soln)%arr(i1-1:cg%ie-1:2, j,   k  ) + cg%q(soln)%arr(i1+1:cg%ie+1:2, j,   k))   + &
!!$                       cg%mg%ry * (cg%q(soln)%arr(i1  :cg%ie  :2, j-1, k  ) + cg%q(soln)%arr(i1:  cg%ie:  2, j+1, k))   + &
!!$                       cg%mg%rz * (cg%q(soln)%arr(i1  :cg%ie  :2, j,   k-1) + cg%q(soln)%arr(i1:  cg%ie:  2, j,   k+1)) - &
!!$                       cg%mg%r  *  cg%q(src)%arr( i1  :cg%ie  :2, j,   k  )  + &
!!$                       cg%mg%rx * (cg%q(soln)%arr(i1+1:cg%ie+1:2, j,   k  ) - cg%q(soln)%arr(i1-1:cg%ie-1:2, j,   k)) * fac(i1:cg%ie:2)
                        call die("[multigrid_gravity:approximate_solution_rbgs] This variant of relaxation loop is not implemented for cylindrical coordinates.")
                     else
                        cg%q(soln)%arr(i1  :cg%ie  :2, j,   k) = &
                             cg%mg%rx * (cg%q(soln)%arr(i1-1:cg%ie-1:2, j,   k)   + cg%q(soln)%arr(i1+1:cg%ie+1:2, j,   k))   + &
                             cg%mg%ry * (cg%q(soln)%arr(i1  :cg%ie  :2, j-1, k)   + cg%q(soln)%arr(i1:  cg%ie:  2, j+1, k))   + &
                             cg%mg%rz * (cg%q(soln)%arr(i1  :cg%ie  :2, j,   k-1) + cg%q(soln)%arr(i1:  cg%ie:  2, j,   k+1)) - &
                             cg%mg%r  *  cg%q(src)%arr (i1  :cg%ie  :2, j,   k)
                     endif
                  enddo
               enddo
            else
               ! In 3D this variant significantly increases instruction count and also some data read
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

               if (kd == RED_BLACK) k1 = cg%ks + mod(n, RED_BLACK)
               select case (dom%geometry_type)
                  case (GEO_XYZ)
                     do k = k1, cg%ke, kd
                        if (jd == RED_BLACK) j1 = cg%js + mod(n+k, RED_BLACK)
                        do j = j1, cg%je, jd
                           if (id == RED_BLACK) i1 = cg%is + mod(n+j+k, RED_BLACK)
                           cg%q(soln)%arr                           (i1:  cg%ie  :id, j,   k)   = &
                                & (1. - Jacobi_damp)* cg%q(soln)%arr(i1  :cg%ie  :id, j,   k)   - &
                                &       Jacobi_damp * cg%q(src)%arr (i1  :cg%ie  :id, j,   k)   * cg%mg%r
                           if (dom%has_dir(xdim))     cg%q(soln)%arr(i1  :cg%ie  :id, j,   k)   = cg%q(soln)%arr(i1:  cg%ie:  id, j,   k)    + &
                                &       Jacobi_damp *(cg%q(soln)%arr(i1-1:cg%ie-1:id, j,   k)   + cg%q(soln)%arr(i1+1:cg%ie+1:id, j,   k))   * cg%mg%rx
                           if (dom%has_dir(ydim))     cg%q(soln)%arr(i1  :cg%ie  :id, j,   k)   = cg%q(soln)%arr(i1:  cg%ie:  id, j,   k)    + &
                                &       Jacobi_damp *(cg%q(soln)%arr(i1  :cg%ie  :id, j-1, k)   + cg%q(soln)%arr(i1:  cg%ie:  id, j+1, k))   * cg%mg%ry
                           if (dom%has_dir(zdim))     cg%q(soln)%arr(i1  :cg%ie  :id, j,   k)   = cg%q(soln)%arr(i1:  cg%ie:  id, j,   k)    + &
                                &       Jacobi_damp *(cg%q(soln)%arr(i1  :cg%ie  :id, j,   k-1) + cg%q(soln)%arr(i1:  cg%ie:  id, j,   k+1)) * cg%mg%rz
                        enddo
                     enddo
                  case (GEO_RPZ)
                     do k = k1, cg%ke, kd
                        if (jd == RED_BLACK) j1 = cg%js + mod(n+k, RED_BLACK)
                        do j = j1, cg%je, jd
                           if (id == RED_BLACK) i1 = cg%is + mod(n+j+k, RED_BLACK)
                           do i = i1, cg%ie, id
                              cr  = overrelax / 2.
                              crx = cg%dvol2 * cg%idx2 * cg%x(i)**2
                              cry = cg%dvol2 * cg%idy2
                              crz = cg%dvol2 * cg%idz2 * cg%x(i)**2
                              cr  = cr / (crx + cry + crz)
                              crx = overrelax_x * crx * cr
                              cry = overrelax_y * cry * cr
                              crz = overrelax_z * crz * cr
                              cr  = cr * cg%dvol2 * cg%x(i)**2

                              crx1 = 2. * cg%x(i) * cg%idx
                              if (crx1 /= 0.) crx1 = 1./crx1
                              cg%q(soln)%arr                           (i,   j,   k)   = &
                                   & (1. - Jacobi_damp)* cg%q(soln)%arr(i,   j,   k)   - &
                                   &       Jacobi_damp * cg%q(src)%arr (i,   j,   k)   * cr
                              if (dom%has_dir(xdim))     cg%q(soln)%arr(i,   j,   k)   = cg%q(soln)%arr(i,   j,   k)    + &
                                   &       Jacobi_damp *(cg%q(soln)%arr(i-1, j,   k)   + cg%q(soln)%arr(i+1, j,   k))   * crx + &
                                   &       Jacobi_damp *(cg%q(soln)%arr(i+1, j,   k)   - cg%q(soln)%arr(i-1, j,   k))   * crx * crx1
                              if (dom%has_dir(ydim))     cg%q(soln)%arr(i,   j,   k)   = cg%q(soln)%arr(i,   j,   k)    + &
                                   &       Jacobi_damp *(cg%q(soln)%arr(i,   j-1, k)   + cg%q(soln)%arr(i,   j+1, k))   * cry
                              if (dom%has_dir(zdim))     cg%q(soln)%arr(i,   j,   k)   = cg%q(soln)%arr(i,   j,   k)    + &
                                   &       Jacobi_damp *(cg%q(soln)%arr(i,   j,   k-1) + cg%q(soln)%arr(i,   j,   k+1)) * crz
                           enddo
                        enddo
                     enddo
                  case default
                     call die("[multigrid_gravity:approximate_solution_rbgs] Unsupported geometry.")
               end select
            endif
            cgl => cgl%nxt
         enddo

         if (dirty_debug) then
            write(dirty_label, '(a,i5)')"relax soln+ smoo=", n
            call check_dirty(curl, soln, dirty_label)
         endif

      enddo

   end subroutine approximate_solution_rbgs

!!$ ============================================================================
!>
!! \brief FFT given-boundary Poisson solver applied to local domain. Should require less communication than RBGS implementation.
!!
!! \todo test a configuration with wider area being subjected to FFT (sizes would no longer be 2**n) to avoid the need of relaxation
!<

   subroutine approximate_solution_fft(curl, src, soln)

      use constants,         only: LO, HI, ndims, xdim, ydim, zdim, GEO_XYZ, half, I_ONE, idm2
      use dataio_pub,        only: die, warn
      use domain,            only: dom
      use gc_list,           only: cg_list_element
      use cg_list_lev,       only: cg_list_level
      use grid_cont,         only: grid_container
      use multigridhelpers,  only: dirty_debug, check_dirty, dirtyL
      use multigridmpifuncs, only: mpi_multigrid_bnd
      use multigridvars,     only: base, extbnd_antimirror, single_base
      use named_array,       only: p3

      implicit none

      type(cg_list_level), pointer, intent(in) :: curl !< pointer to a level for which we approximate the solution
      integer,                      intent(in) :: src  !< index of source in cg%q(:)
      integer,                      intent(in) :: soln !< index of solution in cg%q(:)

      integer :: nf, n, nsmoo
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg
      integer, dimension(ndims,LO:HI) :: pdn, D_2
      integer                         :: dir

      if (curl%first%cg%mg%fft_type == fft_none) call die("[multigrid_gravity:approximate_solution_fft] unknown FFT type")

      if (dom%geometry_type /= GEO_XYZ) call die("[multigrid_gravity:approximate_solution_fft] FFT is not allowed in non-cartesian coordinates.")

      do nf = 1, nsmoof
         cgl => curl%first
         do while (associated(cgl))
            cg => cgl%cg ; p3 => cg%q(src)%span(cg%ijkse)
            cg%mg%src(:, :, :) = p3
            cgl => cgl%nxt
         enddo

         if (curl%first%cg%mg%fft_type == fft_dst) then !correct boundaries on non-periodic local domain
            if (nf == 1 .and. .not. associated(curl, base)) then
               call make_face_boundaries(curl, soln)
            else
               call mpi_multigrid_bnd(curl, soln, I_ONE, extbnd_antimirror)
               cgl => curl%first
               do while (associated(cgl))
                  cg => cgl%cg
                  if (dom%has_dir(xdim)) then
                     cg%mg%bnd_x(:, :, LO) = half* sum (cg%q(soln)%arr(cg%is-1:cg%is, cg%js:cg%je, cg%ks:cg%ke), 1)
                     cg%mg%bnd_x(:, :, HI) = half* sum (cg%q(soln)%arr(cg%ie:cg%ie+1, cg%js:cg%je, cg%ks:cg%ke), 1)
                  endif
                  if (dom%has_dir(ydim)) then
                     cg%mg%bnd_y(:, :, LO) = half* sum (cg%q(soln)%arr(cg%is:cg%ie, cg%js-1:cg%js, cg%ks:cg%ke), 2)
                     cg%mg%bnd_y(:, :, HI) = half* sum (cg%q(soln)%arr(cg%is:cg%ie, cg%je:cg%je+1, cg%ks:cg%ke), 2)
                  endif
                  if (dom%has_dir(zdim)) then
                     cg%mg%bnd_z(:, :, LO) = half* sum (cg%q(soln)%arr(cg%is:cg%ie, cg%js:cg%je, cg%ks-1:cg%ks), 3)
                     cg%mg%bnd_z(:, :, HI) = half* sum (cg%q(soln)%arr(cg%is:cg%ie, cg%js:cg%je, cg%ke:cg%ke+1), 3)
                  endif
                  cgl => cgl%nxt
               enddo
            endif

            cgl => curl%first
            do while (associated(cgl))
               cg => cgl%cg
               if (dirty_debug) then
                  if (dom%has_dir(xdim) .and. any(abs(curl%first%cg%mg%bnd_x(:, :, :)) > dirtyL)) call warn("approximate_solution_fft dirty bnd_x")
                  if (dom%has_dir(ydim) .and. any(abs(curl%first%cg%mg%bnd_y(:, :, :)) > dirtyL)) call warn("approximate_solution_fft dirty bnd_y")
                  if (dom%has_dir(zdim) .and. any(abs(curl%first%cg%mg%bnd_z(:, :, :)) > dirtyL)) call warn("approximate_solution_fft dirty bnd_z")
               endif

               if (dom%has_dir(xdim)) then
                  cg%mg%src(1,      :, :) = cg%mg%src(1,      :, :) - cg%mg%bnd_x(:, :, LO) * 2. * cg%idx2
                  cg%mg%src(cg%nxb, :, :) = cg%mg%src(cg%nxb, :, :) - cg%mg%bnd_x(:, :, HI) * 2. * cg%idx2
               endif
               if (dom%has_dir(ydim)) then
                  cg%mg%src(:, 1,      :) = cg%mg%src(:, 1,      :) - cg%mg%bnd_y(:, :, LO) * 2. * cg%idy2
                  cg%mg%src(:, cg%nyb, :) = cg%mg%src(:, cg%nyb, :) - cg%mg%bnd_y(:, :, HI) * 2. * cg%idy2
               endif
               if (dom%has_dir(zdim)) then
                  cg%mg%src(:, :, 1     ) = cg%mg%src(:, :, 1     ) - cg%mg%bnd_z(:, :, LO) * 2. * cg%idz2
                  cg%mg%src(:, :, cg%nzb) = cg%mg%src(:, :, cg%nzb) - cg%mg%bnd_z(:, :, HI) * 2. * cg%idz2
               endif
               cgl => cgl%nxt
            enddo
         endif

         call fft_convolve(curl)

         cgl => curl%first
         do while (associated(cgl))
            cg => cgl%cg
            p3 => cg%q(soln)%span(cg%ijkse) ; p3 = cg%mg%src(:, :, :)
            cgl => cgl%nxt
         enddo

         call check_dirty(curl, soln, "approx_soln fft+")

         !> \deprecated BEWARE use dom%has_dir() here in a way that does not degrade performance

         if (associated(curl, base) .and. single_base) then !> do not relax if it is a whole domain (BEWARE: simplified check)
            nsmoo = 0
         else
            nsmoo = nsmool
         endif

         !relax the boundaries
         do n = 1, nsmoo
            call mpi_multigrid_bnd(curl, soln, I_ONE, extbnd_antimirror)
            ! Possible optimization: This is a quite costly part of the local FFT solver
            cgl => curl%first
            do while (associated(cgl))
               cg => cgl%cg
               if (fft_full_relax) then
                  p3 => cg%q(soln)%span(cg%ijkse)
                  if (dom%eff_dim == ndims .and. .not. multidim_code_3D) then
                     p3 = cg%mg%rx * (cg%q(soln)%span(cg%ijkse-idm2(xdim,:,:)) + cg%q(soln)%span(cg%ijkse+idm2(xdim,:,:))) + &
                          cg%mg%ry * (cg%q(soln)%span(cg%ijkse-idm2(ydim,:,:)) + cg%q(soln)%span(cg%ijkse+idm2(ydim,:,:))) + &
                          cg%mg%rz * (cg%q(soln)%span(cg%ijkse-idm2(zdim,:,:)) + cg%q(soln)%span(cg%ijkse+idm2(zdim,:,:))) - &
                          cg%mg%r  *  cg%q(src)%span (cg%ijkse)
                  else

                     call die("[multigrid_gravity:approximate_solution_fft] fft_full_relax is allowed only for 3D at the moment")

                     ! An additional array (cg%mg%src would be good enough) is required here to assemble partial results or use red-black passes
                     p3 = - cg%mg%r * cg%q(src)%span(cg%ijkse)
                     if (dom%has_dir(xdim)) p3 = p3 + cg%mg%rx * (cg%q(soln)%span(cg%ijkse-idm2(xdim,:,:)) + cg%q(soln)%span(cg%ijkse+idm2(xdim,:,:)))
                     if (dom%has_dir(ydim)) p3 = p3 + cg%mg%ry * (cg%q(soln)%span(cg%ijkse-idm2(ydim,:,:)) + cg%q(soln)%span(cg%ijkse+idm2(ydim,:,:)))
                     if (dom%has_dir(zdim)) p3 = p3 + cg%mg%rz * (cg%q(soln)%span(cg%ijkse-idm2(zdim,:,:)) + cg%q(soln)%span(cg%ijkse+idm2(zdim,:,:)))
                  endif
               else
                  ! relax only two layers of cells (1 is  significantly worse, 3 does not improve much)
                  ! edges are relaxed twice, corners are relaxed three times which seems to be good

                  D_2 = spread(dom%D_,2,2)
                  do dir = xdim, zdim
                     if (dom%has_dir(dir)) then
                        pdn = cg%ijkse ; pdn(dir,HI) = pdn(dir,HI) + dom%D_(dir) ! -X/-Y/-Z
                        p3 => cg%q(soln)%span(pdn)
                        p3 = cg%mg%rx * (cg%q(soln)%span(pdn-idm2(xdim,:,:)*D_2) + cg%q(soln)%span(pdn+idm2(xdim,:,:)*D_2)) + &
                             cg%mg%ry * (cg%q(soln)%span(pdn-idm2(ydim,:,:)*D_2) + cg%q(soln)%span(pdn+idm2(ydim,:,:)*D_2)) + &
                             cg%mg%rz * (cg%q(soln)%span(pdn-idm2(zdim,:,:)*D_2) + cg%q(soln)%span(pdn+idm2(zdim,:,:)*D_2)) - cg%mg%r  *  cg%q(src)%span(pdn)

                        pdn = cg%ijkse ; pdn(dir,HI) = pdn(dir,LO) - dom%D_(dir) ! +X/+Y/+Z
                        p3 => cg%q(soln)%span(pdn)
                        p3 = cg%mg%rx * (cg%q(soln)%span(pdn-idm2(xdim,:,:)*D_2) + cg%q(soln)%span(pdn+idm2(xdim,:,:)*D_2)) + &
                             cg%mg%ry * (cg%q(soln)%span(pdn-idm2(ydim,:,:)*D_2) + cg%q(soln)%span(pdn+idm2(ydim,:,:)*D_2)) + &
                             cg%mg%rz * (cg%q(soln)%span(pdn-idm2(zdim,:,:)*D_2) + cg%q(soln)%span(pdn+idm2(zdim,:,:)*D_2)) - cg%mg%r  *  cg%q(src)%span(pdn)
                     endif
                  enddo

               endif

               cgl => cgl%nxt
            enddo
         enddo

         call check_dirty(curl, soln, "approx_soln relax+")

      enddo

   end subroutine approximate_solution_fft

!!$ ============================================================================
!>
!! \brief This routine prepares boundary values for local-FFT solver
!<

   subroutine make_face_boundaries(curl, soln)

      use dataio_pub,         only: warn
      use cg_list_lev,        only: cg_list_level
      use mpisetup,           only: nproc
      use multigridbasefuncs, only: zero_boundaries, prolong_faces
      use multigridvars,      only: bnd_givenval, bnd_periodic, base, single_base

      implicit none

      type(cg_list_level), pointer, intent(in) :: curl !< pointer to a level for which we approximate the solution
      integer,                      intent(in) :: soln !< index of solution in cg%q(:)

      if (grav_bnd == bnd_periodic .and. (nproc == 1 .or. (associated(curl, base) .and. single_base) ) ) then
         call zero_boundaries(curl)
      else
         if (.not. associated(curl, base)) then
            call prolong_faces(curl, soln)
         else
            if (grav_bnd /= bnd_givenval) call zero_boundaries(curl)
            call warn("m:mfb WTF?")
         endif
      endif

   end subroutine make_face_boundaries

!>
!! \brief Solve finest level if allowed (single cg and single thread)
!<

   subroutine fft_solve_roof

      use dataio_pub,    only: die
      use domain,        only: is_multicg
      use grid_cont,     only: grid_container
      use multigridvars, only: roof, source, solution
      use named_array,   only: p3

      implicit none

      type(grid_container), pointer :: cg

      if (is_multicg) call die("[multigrid_gravity:fft_solve_roof] multicg not possible")

      cg => roof%first%cg

      if (cg%mg%fft_type == fft_none) return

      p3 => cg%q(source)%span(cg%ijkse)
      cg%mg%src(:, :, :) = p3
      call fft_convolve(roof)
      p3 = cg%mg%src(:, :, :)

   end subroutine fft_solve_roof

!!$ ============================================================================
!>
!! \brief Do the FFT convolution
!<

   subroutine fft_convolve(curl)

      use dataio_pub,  only: die
      use gc_list,     only: cg_list_element
      use cg_list_lev, only: cg_list_level
      use grid_cont,   only: grid_container

      implicit none

      type(cg_list_level), pointer, intent(in) :: curl !< pointer to a level at which make the convolution

      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg

      cgl => curl%first
      do while (associated(cgl))
         cg => cgl%cg

         ! do the convolution in Fourier space; cg%mg%src(:,:,:) -> cg%mg%fft{r}(:,:,:)
         call dfftw_execute(cg%mg%planf)

         select case (cg%mg%fft_type)
            case (fft_rcr)
               cg%mg%fft  = cg%mg%fft  * cg%mg%Green3D
            case (fft_dst)
               cg%mg%fftr = cg%mg%fftr * cg%mg%Green3D
            case default
               call die("[multigrid_gravity:fft_convolve] Unknown FFT type.")
         end select

         call dfftw_execute(cg%mg%plani) ! cg%mg%fft{r}(:,:,:) -> cg%mg%src(:,:,:)
         cgl => cgl%nxt
      enddo

   end subroutine fft_convolve

end module multigrid_gravity
