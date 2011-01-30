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
#if defined(MULTIGRID) && defined(GRAV)

   use mpisetup,      only: cbuff_len
   use multigridvars, only: vcycle_stats

   implicit none

   private
   public :: init_multigrid_grav, init_multigrid_grav_post, cleanup_multigrid_grav, multigrid_solve_grav
   public :: grav_bnd

   ! multigrid constants
   integer, parameter :: fft_rcr=1                                    !< type of FFT transform: full
   integer, parameter :: fft_dst=fft_rcr+1                            !< type of FFT transform: discrete sine
   integer, parameter :: fft_none=-1                                  !< type of FFT transform: none

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
   integer            :: nsmoob                                       !< smoothing cycles on base level when gb_no_fft = .true. (a convergence check would be much better)
   integer            :: nsmoof                                       !< FFT iterations per call
   integer            :: ord_laplacian                                !< Laplace operator order; allowed values are 2 (default) and 4 (experimental, not fully implemented)
   integer            :: ord_time_extrap                              !< Order of temporal extrapolation for solution recycling; -1 means 0-guess, 2 does parabolic interpolation
   logical            :: trust_fft_solution                           !< Bypass the V-cycle, when doing FFT on whole domain, make sure first that FFT is properly set up.
   logical            :: gb_no_fft                                    !< Deny solving the base level with FFT. Can be very slow.
   logical            :: prefer_rbgs_relaxation                       !< Prefer relaxation over FFT local solver. Typically faster.
   !> \todo allow to perform one or more V-cycles with FFT method, the switch to the RBGS (may save one V-cycle in some cases)
   logical            :: fft_full_relax                               !< Perform full or boundary relaxation after local FFT solve
   logical            :: gb_solve_gather                              !< Prefer MPI_Gather over Send/Recv when solving global base level (looks a bit faster on small domains)
   logical            :: fft_patient                                  !< Spend more time in init_multigrid to find faster fft plan
   character(len=cbuff_len) :: grav_bnd_str                           !< Type of gravitational boundary conditions.

   ! boundaries
   integer                     :: grav_bnd                            !< boundary type for computational domain
   !integer                     :: grav_extbnd_mode                    !< external boundary mode

   ! global base-level FFT solver
   real,       dimension(:,:,:,:), allocatable :: gb_src_temp         !< Storage for collected base level if using gb_solve_gather

   ! constants from fftw3.f
   integer, parameter :: FFTW_MEASURE=0, FFTW_PATIENT=32, FFTW_ESTIMATE=64
   integer, parameter :: FFTW_RODFT01=8, FFTW_RODFT10=9

   integer            :: fftw_flags = FFTW_MEASURE                    !< or FFTW_PATIENT on request

   ! solution recycling
   integer, parameter :: nold_max=3                                   !< maximum implemented extrapolation order
   integer :: nold                                                    !< number of old solutions kept for solution recycling
   type :: old_soln                                                   !< container for an old solution with its timestamp
      real, dimension(:,:,:), allocatable :: soln
      real :: time
   end type old_soln
   type :: soln_history                                               !< container for a set of several old potential solutions
      type(old_soln), dimension(nold_max) :: old
      integer :: last                                                 !< index of the last stored potential
      logical :: valid                                                !< .true. when old(last) was properly initialized
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
!! <tr><td>prefer_rbgs_relaxation</td><td>.false.</td><td>logical        </td><td>\copydoc multigrid_gravity::prefer_rbgs_relaxation</td></tr>
!! <tr><td>gb_no_fft             </td><td>.false.</td><td>logical        </td><td>\copydoc multigrid_gravity::gb_no_fft             </td></tr>
!! <tr><td>gb_solve_gather       </td><td>.false.</td><td>logical        </td><td>\copydoc multigrid_gravity::gb_solve_gather       </td></tr>
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
!! <tr><td>grav_bnd_str          </td><td>"periodic"/"dirichlet"</td><td>string of chars</td><td>\copydoc multigrid_gravity::grav_bnd_str          </td></tr>
!! </table>
!! \n \n
!<
   subroutine init_multigrid_grav

      use multigridvars,      only: bnd_periodic, bnd_dirichlet, bnd_isolated, bnd_invalid, correction, mg_nb, ngridvars, periodic_bnd_cnt, non_periodic_bnd_cnt
      use multipole,          only: use_point_monopole, lmax, mmax, ord_prolong_mpole, coarsen_multipole, interp_pt2mom, interp_mom2pot
      use mpisetup,           only: buffer_dim, comm, ierr, master, slave, ibuff, cbuff, rbuff, lbuff, dom, has_dir, geometry
      use mpi,                only: MPI_CHARACTER, MPI_DOUBLE_PRECISION, MPI_INTEGER, MPI_LOGICAL
      use dataio_pub,         only: par_file, ierrh, namelist_errh, compare_namelist, cmdl_nml  ! QA_WARN required for diff_nml
      use dataio_pub,         only: msg, die, warn

      implicit none

      logical, save                    :: frun = .true.          !< First run flag

      namelist /MULTIGRID_GRAVITY/ norm_tol, vcycle_abort, max_cycles, nsmool, nsmoob, &
           &                       overrelax, overrelax_x, overrelax_y, overrelax_z, Jacobi_damp, L4_strength, nsmoof, ord_laplacian, ord_time_extrap, &
           &                       prefer_rbgs_relaxation, gb_no_fft, gb_solve_gather, fft_full_relax, fft_patient, trust_fft_solution, &
           &                       coarsen_multipole, lmax, mmax, ord_prolong_mpole, use_point_monopole, interp_pt2mom, interp_mom2pot, &
           &                       grav_bnd_str

      if (.not.frun) call die("[multigrid_gravity:init_multigrid_grav] Called more than once.")
      frun = .false.

      if (geometry /= "cartesian") call die("[multigrid_gravity:init_multigrid_grav] non-cartesian geometry not implemented yet.")

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
      gb_no_fft              = .false.
      prefer_rbgs_relaxation = .false.
      fft_full_relax         = .false.
      gb_solve_gather        = .false.
      fft_patient            = .false.
      interp_pt2mom          = .false.
      interp_mom2pot         = .false.

      if (periodic_bnd_cnt == 2*count(has_dir(:))) then
         grav_bnd_str = "periodic"
      else
         grav_bnd_str = "dirichlet"
      endif

      if (master) then

         diff_nml(MULTIGRID_GRAVITY)

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
         lbuff(3) = gb_no_fft
         lbuff(4) = prefer_rbgs_relaxation
         lbuff(5) = fft_full_relax
         lbuff(6) = gb_solve_gather
         lbuff(7) = fft_patient
         lbuff(8) = interp_pt2mom
         lbuff(9) = interp_mom2pot

         cbuff(1) = grav_bnd_str

      endif

      call MPI_Bcast(cbuff, cbuff_len*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_Bcast(ibuff,           buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_Bcast(rbuff,           buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)
      call MPI_Bcast(lbuff,           buffer_dim, MPI_LOGICAL,          0, comm, ierr)

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
         gb_no_fft               = lbuff(3)
         prefer_rbgs_relaxation  = lbuff(4)
         fft_full_relax          = lbuff(5)
         gb_solve_gather         = lbuff(6)
         fft_patient             = lbuff(7)
         interp_pt2mom           = lbuff(8)
         interp_mom2pot          = lbuff(9)

         grav_bnd_str   = cbuff(1)(1:len(grav_bnd_str))

      endif

      ngridvars = max(ngridvars, correction)

      ! boundaries
      grav_bnd = bnd_invalid
      select case (grav_bnd_str)
         case ("isolated", "iso")
            grav_bnd = bnd_isolated
         case ("periodic", "per")
            if (dom%bnd_xl_dom /= 'per' .or. dom%bnd_xr_dom /= 'per' .or. dom%bnd_yl_dom /= 'per' .or. dom%bnd_yr_dom /= 'per' .or. dom%bnd_zl_dom /= 'per' .or. dom%bnd_zr_dom /= 'per') &
                 call die("[multigrid_gravity:init_multigrid_grav] cannot enforce periodic boundaries for gravity on a not fully periodic domain")
            grav_bnd = bnd_periodic
         case ("dirichlet", "dir")
            grav_bnd = bnd_dirichlet
         case default
            call die("[multigrid_gravity:init_multigrid_grav] Non-recognized boundary description.")
      end select

      if (periodic_bnd_cnt == 2*count(has_dir(:))) then ! fully periodic domain
         if (grav_bnd /= bnd_periodic .and. master) call warn("[multigrid_gravity:init_multigrid_grav] Ignoring non-periodic boundary conditions for gravity on a fully periodic domain.")
         grav_bnd = bnd_periodic
      else if (periodic_bnd_cnt > 0 .and. non_periodic_bnd_cnt > 0) then
         if (master) call warn("[multigrid_gravity:init_multigrid_grav] Mixing periodic and non-periodic boundary conditions for gravity is experimental.")
         prefer_rbgs_relaxation = .true.
         gb_no_fft = .true.
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

      if (.not. (grav_bnd == bnd_periodic .or. grav_bnd == bnd_dirichlet .or. grav_bnd == bnd_isolated) .and. .not. gb_no_fft) then
         gb_no_fft = .true.
         if (master) call warn("[multigrid_gravity:init_multigrid_grav] Use of FFT not allowed by current boundary type/combination.")
      endif

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

      !! Sanity checks
      if (abs(ord_laplacian) > 2*mg_nb) call die("[multigrid_gravity:init_multigrid_grav] not enough guardcells for given laplacian operator order")

   end subroutine init_multigrid_grav

!!$ ============================================================================
!>
!! \brief Initialization - continued after allocation of everything interesting
!<

   subroutine init_multigrid_grav_post(mb_alloc)

      use arrays,             only: sgp
      use multigridvars,      only: lvl, roof, base, gb, level_gb, level_max, level_min, bnd_periodic, bnd_dirichlet, bnd_isolated, vcycle_stats
      use mpisetup,           only: master, nproc, psize, xdim, ydim, zdim
      use multigridhelpers,   only: vcycle_stats_init, dirty_debug, dirtyH
      use constants,          only: pi, dpi
      use dataio_pub,         only: die, warn
      use multipole,          only: init_multipole

      implicit none

      real, intent(inout)              :: mb_alloc               !< Allocation counter

      type(soln_history), pointer      :: os
      real, allocatable, dimension(:)  :: kx, ky, kz             !< FFT kernel directional components for convolution
      integer, dimension(6)            :: aerr                   !> \deprecated BEWARE: hardcoded magic integer. Update when you change number of simultaneous error checks
      integer :: i, j, idx

      do idx = level_max, level_min, -1
         ! this should work correctly also when eff_dim < 3
         lvl(idx)%r  = overrelax   / 2.
         lvl(idx)%rx = lvl(idx)%dvol2 * lvl(idx)%idx2
         lvl(idx)%ry = lvl(idx)%dvol2 * lvl(idx)%idy2
         lvl(idx)%rz = lvl(idx)%dvol2 * lvl(idx)%idz2
         lvl(idx)%r  = lvl(idx)%r  / (lvl(idx)%rx + lvl(idx)%ry + lvl(idx)%rz)
         lvl(idx)%rx = overrelax_x * lvl(idx)%rx * lvl(idx)%r
         lvl(idx)%ry = overrelax_y * lvl(idx)%ry * lvl(idx)%r
         lvl(idx)%rz = overrelax_z * lvl(idx)%rz * lvl(idx)%r
         lvl(idx)%r  = lvl(idx)%r  * lvl(idx)%dvol2

         !>
         !! \deprecated BEWARE: some of the above invariants may be not optimally defined - the convergence ratio drops when dx /= dy or dy /= dz or dx /= dz
         !! and overrelaxation factors are required to get any convergence (often poor)
         !<

         if (prefer_rbgs_relaxation) then
            lvl(idx)%fft_type = fft_none
         else if (grav_bnd == bnd_periodic .and. nproc == 1) then
            lvl(idx)%fft_type = fft_rcr
         else if (grav_bnd == bnd_periodic .or. grav_bnd == bnd_dirichlet .or. grav_bnd == bnd_isolated) then
            lvl(idx)%fft_type = fft_dst
         else
            lvl(idx)%fft_type = fft_none
         endif
      enddo

      ! solution recycling
      ord_time_extrap = min(nold_max-1, max(-1, ord_time_extrap))
      nold = ord_time_extrap + 1
      if (nold > 0) then
         do j = 1, 2 ! inner and outer arrays
            if (associated(os)) nullify(os)
            select case (j)
               case (1)
                  os => inner
               case (2)
                  if (grav_bnd == bnd_isolated) os => outer
            end select
            if (associated(os)) then
               do i = 1, nold
                  if ( allocated(os%old(i)%soln) ) call die("[multigrid_gravity:init_multigrid_grav_post] os%old(:)%soln arrays already allocated")
                  allocate( os%old(i)%soln( roof%nx, roof%ny, roof%nz), stat=aerr(i) )
                  mb_alloc = mb_alloc + size(os%old(i)%soln)
                  os%old(i)%time= -HUGE(1.0)
                  if (dirty_debug) os%old(i)%soln(:, :, :) = dirtyH
               enddo
               if (any(aerr(1:nold) /= 0)) call die("[multigrid_gravity:init_multigrid_grav_post] Allocation error: os%old(:)%soln")
               os%valid = .false.
               os%last  = 1
            endif
         enddo
      endif

      sgp(:,:,:) = 0. !Initialize all the guardcells, even those which does not impact the solution

      ! data related to local and global base-level FFT solver
      if (gb_no_fft) then
         gb%fft_type = fft_none
      else
         select case (grav_bnd)
            case (bnd_periodic)
               gb%fft_type = fft_rcr
            case (bnd_dirichlet, bnd_isolated)
               gb%fft_type = fft_dst
            case default
               gb%fft_type = fft_none
               if (master) call warn("[multigrid_gravity:init_multigrid_grav_post] gb_no_fft set but no suitable boundary conditions found. Reverting to RBGS relaxation.")
         end select
      endif

      !special initialization of global base-level FFT-related data
      if (gb%fft_type /= fft_none) then

         !> \deprecated BEWARE only small subset of gb% members is ever initialized

         gb%nxb = base%nxb * psize(xdim)
         gb%nyb = base%nyb * psize(ydim)
         gb%nzb = base%nzb * psize(zdim)

         gb%level = level_gb

         gb%idx2  = base%idx2
         gb%idy2  = base%idy2
         gb%idz2  = base%idz2

         if (gb_solve_gather) then
            if (allocated(gb_src_temp)) call die("[multigrid_gravity:init_multigrid_grav_post] gb_src_temp already allocated")
            allocate(gb_src_temp(base%nxb, base%nyb, base%nzb, 0:nproc-1), stat=aerr(1))
            if (aerr(1) /= 0) call die("[multigrid_gravity:init_multigrid_grav_post] Allocation error: gb_src_temp")
            mb_alloc = mb_alloc + size(gb_src_temp)
         endif

      endif

      ! FFT solver storage and data
      do idx = level_max, level_gb, -1

         lvl(idx)%planf = 0
         lvl(idx)%plani = 0

         if (lvl(idx)%fft_type /= fft_none) then

            select case (lvl(idx)%fft_type)
               case (fft_rcr)
                  lvl(idx)%nxc = lvl(idx)%nxb / 2 + 1
               case (fft_dst)
                  lvl(idx)%nxc = lvl(idx)%nxb
               case default
                  call die("[multigrid_gravity:init_multigrid_grav_post] Unknown FFT type.")
            end select

            if (allocated(lvl(idx)%Green3D) .or. allocated(lvl(idx)%src)) call die("[multigrid_gravity:init_multigrid_grav_post] Green3D or src arrays already allocated")
            allocate(lvl(idx)%Green3D(lvl(idx)%nxc, lvl(idx)%nyb, lvl(idx)%nzb), stat=aerr(1))
            allocate(lvl(idx)%src    (lvl(idx)%nxb, lvl(idx)%nyb, lvl(idx)%nzb), stat=aerr(2))
            if (any(aerr(1:2) /= 0)) call die("[multigrid_gravity:init_multigrid_grav_post] Allocation error: lvl(idx)%Green3D or lvl(idx)%src.")
            mb_alloc = mb_alloc + size(lvl(idx)%Green3D) + size(lvl(idx)%src)

            if (allocated(kx)) deallocate(kx)
            if (allocated(ky)) deallocate(ky)
            if (allocated(kz)) deallocate(kz)
            allocate(kx(lvl(idx)%nxc), stat=aerr(1))
            allocate(ky(lvl(idx)%nyb), stat=aerr(2))
            allocate(kz(lvl(idx)%nzb), stat=aerr(3))
            if (any(aerr(1:3) /= 0)) call die("[multigrid_gravity:init_multigrid_grav_post] Allocation error: k[xyz]")

            select case (lvl(idx)%fft_type)

              ! lvl(idx)%fft_norm is set such that the following sequence gives identity:
              ! call dfftw_execute(lvl(idx)%planf); lvl(idx)%fftr(:, :, :) = lvl(idx)%fftr(:, :, :) * lvl(idx)%fft_norm ; call dfftw_execute(lvl(idx)%plani)

               case (fft_rcr)
                  if (allocated(lvl(idx)%fft)) call die("[multigrid_gravity:init_multigrid_grav_post] fft or Green3D array already allocated")
                  allocate(lvl(idx)%fft(lvl(idx)%nxc, lvl(idx)%nyb, lvl(idx)%nzb), stat=aerr(1))
                  if (aerr(1) /= 0) call die("[multigrid_gravity:init_multigrid_grav_post] Allocation error: fft.")
                  mb_alloc  = mb_alloc + 2*size(lvl(idx)%fft)

                  lvl(idx)%fft_norm = 1. / real( lvl(idx)%nxb * lvl(idx)%nyb * lvl(idx)%nzb ) ! No 4 pi G factor here because the source was already multiplied by it

                  ! FFT local solver initialization for 2nd order (3-point) Laplacian
                  ! sin(k*x-d) - 2.*sin(k*x) + sin(k*x+d) = 2 * (cos(d)-1) * sin(k*x) = -4 * sin(d/2)**2 * sin(k*x)
                  ! For 4th order: a*sin(k*x) + b*(sin(k*x-d) + sin(k*x+d)) + c*(sin(k*x-2*d) + sin(k*x+2*d)), a+2*b+2*c == 0 it would be:
                  ! 4*(a+b+(a+2*b)*cos(d)) * sin(d/2)**2 * sin(k*x)
                  ! For 6th order: a*sin(k*x) + b*(sin(k*x-d) + sin(k*x+d)) + c*(sin(k*x-2*d) + sin(k*x+2*d)) + e*(sin(k*x-3*d) + sin(k*x+3*d)), a+2*b+2*c+2*e == 0 it would be:
                  ! 2*(3*a+4*b+2*c+4*(a+2*b+c)*cos(d)+2*(a+2*(b+c))*cos(2*d)) * sin(d/2)**2 * sin(k*x)
                  ! asymptotically: -d**2/2 for d<pi

                  kx(:) = lvl(idx)%idx2 * (cos(dpi/lvl(idx)%nxb*(/(j, j=0, lvl(idx)%nxc-1)/)) - 1.)
                  ky(:) = lvl(idx)%idy2 * (cos(dpi/lvl(idx)%nyb*(/(j, j=0, lvl(idx)%nyb-1)/)) - 1.)
                  kz(:) = lvl(idx)%idz2 * (cos(dpi/lvl(idx)%nzb*(/(j, j=0, lvl(idx)%nzb-1)/)) - 1.)
                  call dfftw_plan_dft_r2c_3d(lvl(idx)%planf, lvl(idx)%nxb, lvl(idx)%nyb, lvl(idx)%nzb, lvl(idx)%src, lvl(idx)%fft, fftw_flags)
                  call dfftw_plan_dft_c2r_3d(lvl(idx)%plani, lvl(idx)%nxb, lvl(idx)%nyb, lvl(idx)%nzb, lvl(idx)%fft, lvl(idx)%src, fftw_flags)

               case (fft_dst)

                  if (allocated(lvl(idx)%fftr)) call die("[multigrid_gravity:init_multigrid_grav_post] fftr array already allocated")
                  allocate(lvl(idx)%fftr(lvl(idx)%nxc, lvl(idx)%nyb, lvl(idx)%nzb), stat=aerr(1))
                  if (aerr(1) /= 0) call die("[multigrid_gravity:init_multigrid_grav_post] Allocation error: fftr.")
                  mb_alloc  = mb_alloc + size(lvl(idx)%fftr)

                  lvl(idx)%fft_norm = 1. / (8. * real( lvl(idx)%nxb * lvl(idx)%nyb * lvl(idx)%nzb ))
                  kx(:) = lvl(idx)%idx2 * (cos(pi/lvl(idx)%nxb*(/(j, j=1, lvl(idx)%nxc)/)) - 1.)
                  ky(:) = lvl(idx)%idy2 * (cos(pi/lvl(idx)%nyb*(/(j, j=1, lvl(idx)%nyb)/)) - 1.)
                  kz(:) = lvl(idx)%idz2 * (cos(pi/lvl(idx)%nzb*(/(j, j=1, lvl(idx)%nzb)/)) - 1.)
                  call dfftw_plan_r2r_3d(lvl(idx)%planf, lvl(idx)%nxb, lvl(idx)%nyb, lvl(idx)%nzb, lvl(idx)%src,  lvl(idx)%fftr, &
                       &                 FFTW_RODFT10, FFTW_RODFT10, FFTW_RODFT10, fftw_flags)
                  call dfftw_plan_r2r_3d(lvl(idx)%plani, lvl(idx)%nxb, lvl(idx)%nyb, lvl(idx)%nzb, lvl(idx)%fftr, lvl(idx)%src,  &
                       &                 FFTW_RODFT01, FFTW_RODFT01, FFTW_RODFT01, fftw_flags)

               case default
                  call die("[multigrid_gravity:init_multigrid_grav_post] Unknown FFT type.")
            end select

            ! compute Green's function for 7-point 3D discrete laplacian
            do i = 1, lvl(idx)%nxc
               do j = 1, lvl(idx)%nyb
                  where ( (kx(i) + ky(j) + kz(:)) /= 0 )
                     lvl(idx)%Green3D(i,j,:) = 0.5 * lvl(idx)%fft_norm / (kx(i) + ky(j) + kz(:))
                  elsewhere
                     lvl(idx)%Green3D(i,j,:) = 0.0
                  endwhere
               enddo
            enddo

         endif
      enddo

      if (roof%fft_type == fft_none .and. trust_fft_solution) then
         if (master) call warn("[multigrid_gravity:init_multigrid_grav_post] cannot trust FFT solution on the roof.")
         trust_fft_solution = .false.
      endif

      if (allocated(kx)) deallocate(kx)
      if (allocated(ky)) deallocate(ky)
      if (allocated(kz)) deallocate(kz)

      if (grav_bnd == bnd_isolated) call init_multipole(mb_alloc)

      call vcycle_stats_init(vstat, max_cycles)
      mb_alloc = mb_alloc + 2*max_cycles

   end subroutine init_multigrid_grav_post

!!$ ============================================================================
!>
!! \brief Cleanup
!<

   subroutine cleanup_multigrid_grav

      use multipole,          only: cleanup_multipole
      use multigridvars, only: lvl, level_gb, level_max

      implicit none

      integer :: i

      call cleanup_multipole

      do i = 1, nold
         if (allocated(inner%old(i)%soln)) deallocate(inner%old(i)%soln)
         if (allocated(outer%old(i)%soln)) deallocate(outer%old(i)%soln)
      enddo

      if (allocated(gb_src_temp)) deallocate(gb_src_temp)
      if (allocated(vstat%factor)) deallocate(vstat%factor)
      if (allocated(vstat%time)) deallocate(vstat%time)

      if (allocated(lvl)) then
         do i = level_gb, level_max
            if (allocated(lvl(i)%fft))        deallocate(lvl(i)%fft)
            if (allocated(lvl(i)%fftr))       deallocate(lvl(i)%fftr)
            if (allocated(lvl(i)%src))        deallocate(lvl(i)%src)
            if (allocated(lvl(i)%Green3D))    deallocate(lvl(i)%Green3D)

            if (lvl(i)%planf /= 0) call dfftw_destroy_plan(lvl(i)%planf)
            if (lvl(i)%plani /= 0) call dfftw_destroy_plan(lvl(i)%plani)
         enddo
      endif

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

      use mpisetup,         only: master, t
      use multigridhelpers, only: set_dirty, check_dirty, mg_write_log
      use dataio_pub,       only: msg, die
      use multigridvars,    only: lvl, roof, level_min, level_max, stdout, solution

      implicit none

      type(soln_history), intent(inout) :: history !< inner or outer potential history for recycling

      integer :: l, p0, p1, p2, ordt
      real, dimension(3)  :: dt_fac

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
            ordt = min(2, ord_time_extrap)
         else if (history%old(p1)%time < history%old(p0)%time .and. &
              &   history%old(p0)%time < t) then      ! linear extrapolation
            ordt = min(1, ord_time_extrap)
         else                                                            ! simple recycling
            ordt = min(0, ord_time_extrap)
         endif
      else                                                               ! coldstart
         ordt = min(-1, ord_time_extrap)
      endif

      select case (ordt)
         case (:-1)
            if (master .and. ord_time_extrap > -1) then
               write(msg, '(3a)')"[multigrid_gravity:init_solution] Clearing ",trim(vstat%cprefix),"solution."
               call mg_write_log(msg, stdout)
            endif
            do l = level_min, level_max
               lvl(l)%mgvar(:, :, :, solution) = 0.
            enddo
            history%old(:)%time = -HUGE(1.0)
         case (0)
            roof%mgvar(:, :, :, solution) = history%old(p0)%soln(:, :, :)
            if (master .and. ord_time_extrap > 0) then
               write(msg, '(3a)')"[multigrid_gravity:init_solution] No extrapolation of ",trim(vstat%cprefix),"solution."
               call mg_write_log(msg, stdout)
            endif
         case (1)
            dt_fac(1) = (t - history%old(p0)%time) / (history%old(p0)%time - history%old(p1)%time)
            roof%mgvar(:, :, :, solution) = (1. + dt_fac(1)) * history%old(p0)%soln(:, :, :) - dt_fac(1) *  history%old(p1)%soln(:, :, :)
            if (master .and. ord_time_extrap > 1) then
               write(msg, '(3a)')"[multigrid_gravity:init_solution] Linear extrapolation of ",trim(vstat%cprefix),"solution."
               call mg_write_log(msg, stdout)
            endif
         case (2)
            dt_fac(1) = (t - history%old(p0)%time) / (history%old(p1)%time - history%old(p2)%time)
            dt_fac(2) = (t - history%old(p1)%time) / (history%old(p2)%time - history%old(p0)%time)
            dt_fac(3) = (t - history%old(p2)%time) / (history%old(p0)%time - history%old(p1)%time)
            roof%mgvar(:, :, :, solution) = - dt_fac(2) * dt_fac(3) * history%old(p0)%soln(:, :, :) &
                 &                          - dt_fac(1) * dt_fac(3) * history%old(p1)%soln(:, :, :) &
                 &                          - dt_fac(1) * dt_fac(2) * history%old(p2)%soln(:, :, :)
         case default
            call die("[multigrid_gravity:init_solution] Extrapolation order not implemented")
      end select

      call check_dirty(level_max, solution, "init_soln")

   end subroutine init_solution

!!$ ============================================================================
!>
!! \brief Make a local copy of source (density) and multiply by 4 pi G
!<

   subroutine init_source(dens)

#ifdef JEANS_PROBLEM
      use problem_pub,        only: jeans_d0, jeans_mode ! hack for tests
#endif /* JEANS_PROBLEM */
      use constants,          only: fpiG
      use grid,               only: cg
      use dataio_pub,         only: die
      use multigridhelpers,   only: set_dirty, check_dirty
      use multigridbasefuncs, only: substract_average
      use multigridvars,      only: roof, source, level_max, is_external, bnd_periodic, bnd_dirichlet, bnd_givenval, &
           &                        XLO, XHI, YLO, YHI, ZLO, ZHI, LOW, HIGH

      implicit none

      real, optional, dimension(:,:,:), intent(in)  :: dens !< input source field or nothing for empty space

      call set_dirty(source)

      if (present(dens)) then
         roof%mgvar(roof%is:roof%ie, roof%js:roof%je, roof%ks:roof%ke, source) = fpiG * dens(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke)
      else
         if (grav_bnd /= bnd_givenval) call die("[multigrid_gravity:init_source] empty space allowed only for given value boundaries.")
         roof%mgvar(roof%is:roof%ie, roof%js:roof%je, roof%ks:roof%ke, source) = 0.
      endif

      select case (grav_bnd)
         case (bnd_periodic) ! probably also bnd_neumann
            call substract_average(level_max, source)
         case (bnd_dirichlet)
#ifdef JEANS_PROBLEM
            if (jeans_mode == 1) roof%mgvar(roof%is:roof%ie, roof%js:roof%je, roof%ks:roof%ke, source) = &
                 &         roof%mgvar(roof%is:roof%ie, roof%js:roof%je, roof%ks:roof%ke, source) - fpiG * jeans_d0 ! remove density bias
#endif /* JEANS_PROBLEM */
         case (bnd_givenval) ! convert potential into a layer of imaginary mass (subtract second derivative normal to computational domain boundary)
            !> \deprecated BEWARE: cylindrical factors go here
            if (is_external(XLO)) roof%mgvar(roof%is,         roof%js:roof%je, roof%ks:roof%ke, source) = &
                 &                roof%mgvar(roof%is,         roof%js:roof%je, roof%ks:roof%ke, source) - &
                 &                roof%bnd_x(                 roof%js:roof%je, roof%ks:roof%ke, LOW)  * 2. * roof%idx2 / fpiG
            if (is_external(XHI)) roof%mgvar(        roof%ie, roof%js:roof%je, roof%ks:roof%ke, source) = &
                 &                roof%mgvar(        roof%ie, roof%js:roof%je, roof%ks:roof%ke, source) - &
                 &                roof%bnd_x(                 roof%js:roof%je, roof%ks:roof%ke, HIGH) * 2. * roof%idx2 / fpiG
            if (is_external(YLO)) roof%mgvar(roof%is:roof%ie, roof%js,         roof%ks:roof%ke, source) = &
                 &                roof%mgvar(roof%is:roof%ie, roof%js,         roof%ks:roof%ke, source) - &
                 &                roof%bnd_y(roof%is:roof%ie,                  roof%ks:roof%ke, LOW)  * 2. * roof%idy2 / fpiG
            if (is_external(YHI)) roof%mgvar(roof%is:roof%ie,         roof%je, roof%ks:roof%ke, source) = &
                 &                roof%mgvar(roof%is:roof%ie,         roof%je, roof%ks:roof%ke, source) - &
                 &                roof%bnd_y(roof%is:roof%ie,                  roof%ks:roof%ke, HIGH) * 2. * roof%idy2 / fpiG
            if (is_external(ZLO)) roof%mgvar(roof%is:roof%ie, roof%js:roof%je, roof%ks,         source) = &
                 &                roof%mgvar(roof%is:roof%ie, roof%js:roof%je, roof%ks,         source) - &
                 &                roof%bnd_z(roof%is:roof%ie, roof%js:roof%je,                  LOW)  * 2. * roof%idz2 / fpiG
            if (is_external(ZHI)) roof%mgvar(roof%is:roof%ie, roof%js:roof%je,         roof%ke, source) = &
                 &                roof%mgvar(roof%is:roof%ie, roof%js:roof%je,         roof%ke, source) - &
                 &                roof%bnd_z(roof%is:roof%ie, roof%js:roof%je,                  HIGH) * 2. * roof%idz2 / fpiG
         case default
            call die("[multigrid_gravity:init_source] Unknown boundary type")
      end select

      call check_dirty(level_max, source, "init_src")

   end subroutine init_source

!!$ ============================================================================
!>
!! \brief This routine manages old copies of potential for recycling.
!<

   subroutine store_solution(history)

      use mpisetup,          only: t
      use multigridmpifuncs, only: mpi_multigrid_bnd
      use multigridvars,     only: roof, bnd_isolated, bnd_givenval, solution, level_max, mg_nb, extbnd_extrapolate, extbnd_mirror

      implicit none

      type(soln_history), intent(inout) :: history !< inner or outer potential history to store recent solution

      if (nold <= 0) return

      if (history%valid) then
         history%last = 1 + mod(history%last, nold)
      else
         history%old(:)%time = t ! prevents extrapolation too early
      endif

      history%old(history%last)%soln(:, :, :) = roof%mgvar(:, :, :, solution)
      history%old(history%last)%time = t
      history%valid = .true.

      ! Update guardcells of the solution before leaving. This can be done in higher-level routines that collect all the gravity contributions, but would be less safe.
      ! Extrapolate isolated boundaries, remember that grav_bnd is messed up by multigrid_solve_*
      if (grav_bnd == bnd_isolated .or. grav_bnd == bnd_givenval) then
         call mpi_multigrid_bnd(level_max, solution, mg_nb, extbnd_extrapolate)
      else
         call mpi_multigrid_bnd(level_max, solution, mg_nb, extbnd_mirror)
      endif

   end subroutine store_solution

!!$ ============================================================================
!>
!! \brief Multigrid gravity driver. This is the only multigrid routine intended to be called from the gravity module.
!! This routine is also responsible for communicating the solution to the rest of world via sgp array.
!<

   subroutine multigrid_solve_grav(dens)

      use timer,         only: set_timer
      use arrays,        only: sgp
      use grid,          only: cg
      use mpisetup,      only: xdim, ydim, zdim, has_dir
      use dataio_pub,    only: die
      use multipole,     only: multipole_solver
      use multigridvars, only: roof, solution, bnd_isolated, bnd_dirichlet, bnd_givenval, mg_nb, tot_ts, ts

      implicit none

      real, dimension(:,:,:), intent(in) :: dens !< input source field

      logical :: isolated
      integer :: isb, ieb, jsb, jeb, ksb, keb

      ts =  set_timer("multigrid", .true.)
      if ( (has_dir(xdim) .and. cg%is-mg_nb <= 0) .or. &
           (has_dir(ydim) .and. cg%js-mg_nb <= 0) .or. &
           (has_dir(zdim) .and. cg%ks-mg_nb <= 0) )    &
           call die("[multigrid_gravity:multigrid_solve_grav] Current implementation requires at least 2 guardcells in the hydro part")

      isolated = (grav_bnd == bnd_isolated) !> \deprecated BEWARE: not elegant; probably there should be two global grav_bnd variables

      if (isolated) then
         grav_bnd = bnd_dirichlet
         vstat%cprefix = "Gi-"
      else
#ifdef COSM_RAYS
         vstat%cprefix = "G-"
#else
         vstat%cprefix = ""
#endif /* COSM_RAYS */
      endif

      call init_source(dens)

      call vcycle_hg(inner)

      !> \todo move to multigridvars and init_multigrid
      if (has_dir(xdim)) then
         isb = cg%is-mg_nb
         ieb = cg%ie+mg_nb
      else
         isb = 1
         ieb = 1
      endif

      if (has_dir(ydim)) then
         jsb = cg%js-mg_nb
         jeb = cg%je+mg_nb
      else
         jsb = 1
         jeb = 1
      endif

      if (has_dir(zdim)) then
         ksb = cg%ks-mg_nb
         keb = cg%ke+mg_nb
      else
         ksb = 1
         keb = 1
      endif

      sgp(isb:ieb, jsb:jeb, ksb:keb) = roof%mgvar(:, :, :, solution)

      if (isolated) then
         grav_bnd = bnd_givenval

         vstat%cprefix = "Go-"
         call multipole_solver
         call init_source

         call vcycle_hg(outer)

         sgp(isb:ieb, jsb:jeb, ksb:keb) = sgp(isb:ieb, jsb:jeb, ksb:keb) + roof%mgvar(:, :, :, solution)

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

      use mpisetup,           only: master, nproc, cbuff_len
      use timer,              only: set_timer
      use multigridhelpers,   only: set_dirty, check_dirty, mg_write_log, brief_v_log, do_ascii_dump, numbered_ascii_dump
      use multigridbasefuncs, only: norm_sq, restrict_all, substract_average
      use dataio_pub,         only: msg, die, warn
      use multigridvars,      only: roof, base, source, solution, correction, defect, verbose_vcycle, &
           &                        bnd_periodic, level_min, level_max, stdout, tot_ts, ts

      implicit none

      type(soln_history), intent(inout) :: history !< inner or outer potential history used for initializing first guess

      real,    parameter :: suspicious_factor = 1.05 !> \deprecated If the norm decreases too slowly then dump diagnostic output (BEWARE: this option is for tests only)
      integer            :: l, v
      real               :: norm_rhs, norm_lhs, norm_old, norm_lowest
      logical            :: dump_every_step, dump_result
      logical, save      :: norm_was_zero = .false.
      integer, parameter       :: fmtlen = 32
      character(len=fmtlen)    :: fmt
      character(len=cbuff_len) :: dname

      inquire(file = "_dump_every_step_", EXIST=dump_every_step) ! use for debug only
      inquire(file = "_dump_result_", EXIST=dump_result)

      do_ascii_dump = do_ascii_dump .or. dump_every_step .or. dump_result

      ! On single CPU use FFT if possible because it is faster. Can be disabled by prefer_rbgs_relaxation = .true.
      if (nproc == 1 .and. roof%fft_type /= fft_none) then
         call set_dirty(solution)
         call fft_solve_roof
         if (trust_fft_solution) then
            write(msg, '(3a)')"[multigrid_gravity:vcycle_hg] FFT solution trusted, skipping ", trim(vstat%cprefix), "cycle."
            call mg_write_log(msg, stdout)
            return
         endif
      else
         call init_solution(history)
      endif

      call norm_sq(source, norm_rhs)
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

!      if (.not. history%valid .and. prefer_rbgs_relaxation) call approximate_solution(level_max, source, solution) ! not necessary when init_solution called FFT
! difficult statement: for approximate_solution_fft it requires to pass a flag to use guardcells instead of prolonging faces.
! how much does it improve? (make a benchmark at some point)

      ! iterations
      do v = 0, max_cycles

         call set_dirty(defect)
         call residual(level_max, source, solution, defect)
         if (grav_bnd == bnd_periodic) call substract_average(level_max, defect)
         call check_dirty(level_max, defect, "residual")

         call norm_sq(defect, norm_lhs)
         ts = set_timer("multigrid")
         tot_ts = tot_ts + ts
         if (master .and. verbose_vcycle) then
            if (norm_old/norm_lhs < 1.e5) then
               fmt='(3a,i3,a,f12.9,a,f8.2,a,f7.3)'
            else
               fmt='(3a,i3,a,f12.9,a,es8.2,a,f7.3)'
            endif
            write(msg, fmt)"[multigrid_gravity] ", trim(vstat%cprefix), "Cycle:", v, " norm/rhs= ", norm_lhs/norm_rhs, " reduction factor= ", norm_old/norm_lhs, "   dt_wall= ", ts
            call mg_write_log(msg, stdout)
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
                  if (.not. verbose_vcycle) call brief_v_log(vstat)
                  call die("[multigrid_gravity:vcycle_hg] Serious nonconvergence detected.")
                  !In such case one may increase nsmool, decrease refinement depth or use FFT
               endif
            endif
         endif

         norm_old = norm_lhs

         ! the Huang-Greengard V-cycle
         call restrict_all(defect)

         call set_dirty(correction)
         base%mgvar(:, :, :, correction) = 0.

         do l = level_min, level_max
            call approximate_solution(l, defect, correction)
            call check_dirty(l, correction, "Vup relax+")
         enddo
         roof%mgvar     (roof%is:roof%ie, roof%js:roof%je, roof%ks:roof%ke, solution) = &
              roof%mgvar(roof%is:roof%ie, roof%js:roof%je, roof%ks:roof%ke, solution) + &
              roof%mgvar(roof%is:roof%ie, roof%js:roof%je, roof%ks:roof%ke, correction)

      enddo

      if (v > max_cycles) then
         if (master .and. norm_lhs/norm_rhs > norm_tol) call warn("[multigrid_gravity:vcycle_hg] Not enough V-cycles to achieve convergence.")
         v = max_cycles
      endif

      vstat%norm_final = norm_lhs/norm_rhs
      if (.not. verbose_vcycle) call brief_v_log(vstat)

      call check_dirty(level_max, solution, "final_solution")

      call store_solution(history)

   end subroutine vcycle_hg

!!$ ============================================================================
!>
!! \brief Calculate the residuum for the Poisson equation.
!<

   subroutine residual(lev, src, soln, def)

      use dataio_pub,            only: die
      use multigridvars,         only: level_min, level_max, ngridvars

      implicit none

      integer, intent(in) :: lev  !< level for which approximate the solution
      integer, intent(in) :: src  !< index of source in lvl()%mgvar
      integer, intent(in) :: soln !< index of solution in lvl()%mgvar
      integer, intent(in) :: def  !< index of defect in lvl()%mgvar

      if (any( [ src, soln, def ] <= 0) .or. any( [ src, soln, def ] > ngridvars)) call die("[multigrid_gravity:residual] Invalid variable index")
      if (lev < level_min .or. lev > level_max) call die("[multigrid_gravity:residual] Invalid level number")

      select case (ord_laplacian)
      case (2)
         call residual2(lev, src, soln, def)
      case (4)
         call residual4(lev, src, soln, def)
      case default
         call die("[multigrid_gravity:residual] The parameter 'ord_laplacian' must be 2 or 4")
      end select

   end subroutine residual

!!$ ============================================================================
!>
!! \brief 2nd order Laplacian
!<

   subroutine residual2(lev, src, soln, def)

      use mpisetup,           only: xdim, ydim, zdim, has_dir, eff_dim
      use multigridhelpers,   only: multidim_code_3D
      use multigridmpifuncs,  only: mpi_multigrid_bnd
      use multigridvars,      only: lvl, NDIM, extbnd_antimirror

      implicit none

      integer, intent(in) :: lev  !< level for which approximate the solution
      integer, intent(in) :: src  !< index of source in lvl()%mgvar
      integer, intent(in) :: soln !< index of solution in lvl()%mgvar
      integer, intent(in) :: def  !< index of defect in lvl()%mgvar

      real                :: L0, Lx, Ly, Lz
      integer :: k

      call mpi_multigrid_bnd(lev, soln, 1, extbnd_antimirror) ! no corners required

      ! Coefficients for a simplest 3-point Laplacian operator: [ 1, -2, 1 ]
      ! for 2D and 1D setups appropriate elements of [ Lx, Ly, Lz ] should be == 0.
      Lx = lvl(lev)%idx2
      Ly = lvl(lev)%idy2
      Lz = lvl(lev)%idz2
      L0 = -2. * (Lx + Ly + Lz)

      ! Possible optimization candidate: reduce cache misses (secondary importance, cache-aware implementation required)
      ! Explicit loop over k gives here better performance than array operation due to less cache misses (at least on 32^3 and 64^3 arrays)
      !> \deprecated BEWARE: cylindrical factors go here
      if (eff_dim == NDIM .and. .not. multidim_code_3D) then
         do k = lvl(lev)%ks, lvl(lev)%ke
            lvl(       lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js  :lvl(lev)%je,   k,   def)        = &
                 & lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js  :lvl(lev)%je,   k,   src)        - &
                 ( lvl(lev)%mgvar(lvl(lev)%is-1:lvl(lev)%ie-1, lvl(lev)%js  :lvl(lev)%je,   k,   soln)       + &
                 & lvl(lev)%mgvar(lvl(lev)%is+1:lvl(lev)%ie+1, lvl(lev)%js  :lvl(lev)%je,   k,   soln)) * Lx - &
                 ( lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js-1:lvl(lev)%je-1, k,   soln)       + &
                 & lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js+1:lvl(lev)%je+1, k,   soln)) * Ly - &
                 ( lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js  :lvl(lev)%je,   k-1, soln)       + &
                 & lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js  :lvl(lev)%je,   k+1, soln)) * Lz - &
                 & lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js  :lvl(lev)%je,   k,   soln)  * L0
         enddo
      else
         ! In 3D this implementation can give a bit more cache misses, few times more writes and significantly more instructions executed than monolithic 3D above
         do k = lvl(lev)%ks, lvl(lev)%ke
            lvl(       lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js  :lvl(lev)%je,   k,   def)   = &
                 & lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js  :lvl(lev)%je,   k,   src)   - &
                 & lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js  :lvl(lev)%je,   k,   soln)  * L0
            if (has_dir(xdim)) &
                 & lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js  :lvl(lev)%je,   k,   def)   = &
                 & lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js  :lvl(lev)%je,   k,   def)   - &
                 ( lvl(lev)%mgvar(lvl(lev)%is-1:lvl(lev)%ie-1, lvl(lev)%js  :lvl(lev)%je,   k,   soln)  + &
                 & lvl(lev)%mgvar(lvl(lev)%is+1:lvl(lev)%ie+1, lvl(lev)%js  :lvl(lev)%je,   k,   soln)) * Lx
            if (has_dir(ydim)) &
                 & lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js  :lvl(lev)%je,   k,   def)   = &
                 & lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js  :lvl(lev)%je,   k,   def)   - &
                 ( lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js-1:lvl(lev)%je-1, k,   soln)  + &
                 & lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js+1:lvl(lev)%je+1, k,   soln)) * Ly
            if (has_dir(zdim)) &
                 & lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js  :lvl(lev)%je,   k,   def)   = &
                 & lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js  :lvl(lev)%je,   k,   def)   - &
                 ( lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js  :lvl(lev)%je,   k-1, soln)  + &
                 & lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js  :lvl(lev)%je,   k+1, soln)) * Lz
         enddo
      endif

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

   subroutine residual4(lev, src, soln, def)

      use dataio_pub,         only: die, warn
      use mpisetup,           only: master, eff_dim
      use multigridmpifuncs,  only: mpi_multigrid_bnd
      use multigridvars,      only: lvl, NDIM, bnd_givenval, extbnd_antimirror

      implicit none

      integer, intent(in) :: lev  !< level for which approximate the solution
      integer, intent(in) :: src  !< index of source in lvl()%mgvar
      integer, intent(in) :: soln !< index of solution in lvl()%mgvar
      integer, intent(in) :: def  !< index of defect in lvl()%mgvar

      real, parameter     :: L4_scaling = 1./12. ! with L4_strength = 1. this gives an L4 approximation for finite differences approach
      integer, parameter  :: L2w = 2             ! #layers of boundary cells for L2 operator

      real                :: c21, c41, c42 !, c20, c40
      real                :: L0, Lx1, Lx2, Ly1, Ly2, Lz1, Lz2, Lx, Ly, Lz

      logical, save       :: firstcall = .true.
      integer             :: i, j, k

      if (eff_dim<NDIM) call die("[multigrid_gravity:residual4] Only 3D is implemented")

      if (firstcall) then
         if (master) call warn("[multigrid_gravity:residual4] residual order 4 is experimental.")
         firstcall = .false.
      endif

      call mpi_multigrid_bnd(lev, soln, 2, extbnd_antimirror) ! no corners required

      c21 = 1.
      c42 = - L4_scaling * L4_strength
      c41 = c21 + 4. * L4_scaling * L4_strength
      !c20 = -2.
      !c40 = c20 - 6. * L4_strength

      Lx1 = c41 * lvl(lev)%idx2
      Ly1 = c41 * lvl(lev)%idy2
      Lz1 = c41 * lvl(lev)%idz2
      Lx2 = c42 * lvl(lev)%idx2
      Ly2 = c42 * lvl(lev)%idy2
      Lz2 = c42 * lvl(lev)%idz2
!      L0  = c40 * (lvl(lev)%idx2 + lvl(lev)%idy2 + lvl(lev)%idz2 )
      L0 = -2. * (Lx1 + Lx2 + Ly1 + Ly2 + Lz1 + Lz2)

      !> \deprecated BEWARE: cylindrical factors go here
      lvl(     lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js  :lvl(lev)%je,   lvl(lev)%ks  :lvl(lev)%ke,   def)        = &
           lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js  :lvl(lev)%je,   lvl(lev)%ks  :lvl(lev)%ke,   src)        - &
           lvl(lev)%mgvar(lvl(lev)%is-2:lvl(lev)%ie-2, lvl(lev)%js  :lvl(lev)%je,   lvl(lev)%ks  :lvl(lev)%ke,   soln) * Lx2 - &
           lvl(lev)%mgvar(lvl(lev)%is+2:lvl(lev)%ie+2, lvl(lev)%js  :lvl(lev)%je,   lvl(lev)%ks  :lvl(lev)%ke,   soln) * Lx2 - &
           lvl(lev)%mgvar(lvl(lev)%is-1:lvl(lev)%ie-1, lvl(lev)%js  :lvl(lev)%je,   lvl(lev)%ks  :lvl(lev)%ke,   soln) * Lx1 - &
           lvl(lev)%mgvar(lvl(lev)%is+1:lvl(lev)%ie+1, lvl(lev)%js  :lvl(lev)%je,   lvl(lev)%ks  :lvl(lev)%ke,   soln) * Lx1 - &
           lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js-2:lvl(lev)%je-2, lvl(lev)%ks  :lvl(lev)%ke,   soln) * Ly2 - &
           lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js+2:lvl(lev)%je+2, lvl(lev)%ks  :lvl(lev)%ke,   soln) * Ly2 - &
           lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js-1:lvl(lev)%je-1, lvl(lev)%ks  :lvl(lev)%ke,   soln) * Ly1 - &
           lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js+1:lvl(lev)%je+1, lvl(lev)%ks  :lvl(lev)%ke,   soln) * Ly1 - &
           lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js  :lvl(lev)%je,   lvl(lev)%ks-2:lvl(lev)%ke-2, soln) * Lz2 - &
           lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js  :lvl(lev)%je,   lvl(lev)%ks+2:lvl(lev)%ke+2, soln) * Lz2 - &
           lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js  :lvl(lev)%je,   lvl(lev)%ks-1:lvl(lev)%ke-1, soln) * Lz1 - &
           lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js  :lvl(lev)%je,   lvl(lev)%ks+1:lvl(lev)%ke+1, soln) * Lz1 - &
           lvl(lev)%mgvar(lvl(lev)%is  :lvl(lev)%ie,   lvl(lev)%js  :lvl(lev)%je,   lvl(lev)%ks  :lvl(lev)%ke,   soln) * L0

      ! WARNING: not optimized
      if (grav_bnd == bnd_givenval) then ! probably also in some other cases
         ! Use L2 Laplacian in two layers of cells next to the boundary because L4 seems to be incompatible with present image mass construction
         Lx = c21 * lvl(lev)%idx2
         Ly = c21 * lvl(lev)%idy2
         Lz = c21 * lvl(lev)%idz2
         L0 = -2. * (Lx + Ly + Lz)

         do k = lvl(lev)%ks, lvl(lev)%ke
            do j = lvl(lev)%js, lvl(lev)%je
               do i = lvl(lev)%is, lvl(lev)%ie
                  if ( i<lvl(lev)%is+L2w .or. i>lvl(lev)%ie-L2w .or. j<lvl(lev)%js+L2w .or. j>lvl(lev)%je-L2w .or. k<lvl(lev)%ks+L2w .or. k>lvl(lev)%ke-L2w) then
                     lvl(       lev)%mgvar(i,   j,   k,   def)   = lvl(lev)%mgvar(i,   j,   k,   src)        - &
                          ( lvl(lev)%mgvar(i-1, j,   k,   soln)  + lvl(lev)%mgvar(i+1, j,   k,   soln)) * Lx - &
                          ( lvl(lev)%mgvar(i,   j-1, k,   soln)  + lvl(lev)%mgvar(i,   j+1, k,   soln)) * Ly - &
                          ( lvl(lev)%mgvar(i,   j,   k-1, soln)  + lvl(lev)%mgvar(i,   j,   k+1, soln)) * Lz - &
                          & lvl(lev)%mgvar(i,   j,   k,   soln)  * L0
                  endif
               enddo
            enddo
         enddo
      endif

   end subroutine residual4

!!$ ============================================================================
!>
!! \brief This routine has to find an approximate solution for given source field and implemented differential operator
!<

   subroutine approximate_solution(lev, src, soln)

      use dataio_pub,         only: die
      use multigridhelpers,   only: check_dirty
      use multigridbasefuncs, only: prolong_level
      use multigridvars,      only: level_min, level_max, ngridvars, correction

      implicit none

      integer, intent(in) :: lev  !< level for which approximate the solution
      integer, intent(in) :: src  !< index of source in lvl()%mgvar
      integer, intent(in) :: soln !< index of solution in lvl()%mgvar

      if (any( [ src, soln ] <= 0) .or. any( [ src, soln ] > ngridvars)) call die("[multigrid_gravity:approximate_solution] Invalid variable index.")
      if (lev < level_min .or. lev > level_max) call die("[multigrid_gravity:approximate_solution] Invalid level number.")

      call check_dirty(lev, src, "approx_soln src-")

      if (lev == level_min .and. .not. gb_no_fft) then
         call gb_fft_solve(src, soln)
      else
         if (prefer_rbgs_relaxation) then
            call check_dirty(lev, soln, "approx_soln soln-")
            call approximate_solution_rbgs(lev, src, soln)
         else
            call approximate_solution_fft(lev, src, soln)
         endif
      endif

      if (prefer_rbgs_relaxation .and. soln == correction .and. lev <  level_max) call prolong_level(lev, correction)
      !> \deprecated BEWARE other implementations of the multigrid algorithm may be incompatible with prolongation called from here

      call check_dirty(lev, soln, "approx_soln soln+")

   end subroutine approximate_solution

!!$ ============================================================================
!>
!! \brief Red-Black Gauss-Seidel relaxation.
!!
!! \details This is the most costly routine in a serial run. Try to find optimal values for nsmool and nsmoob.
!! This routine also depends a lot on communication so it  may limit scalability of the multigrid.
!! \todo Implement convergence check on base level (not very important since we have a FFT solver for base level)
!<

   subroutine approximate_solution_rbgs(lev, src, soln)

      use mpisetup,           only: xdim, ydim, zdim, has_dir, eff_dim
      use multigridhelpers,   only: dirty_debug, check_dirty, multidim_code_3D, dirty_label
      use multigridmpifuncs,  only: mpi_multigrid_bnd
      use multigridvars,      only: lvl, level_min, NDIM, extbnd_antimirror

      implicit none

      integer, intent(in) :: lev  !< level for which approximate the solution
      integer, intent(in) :: src  !< index of source in lvl()%mgvar
      integer, intent(in) :: soln !< index of solution in lvl()%mgvar

      integer, parameter :: RED_BLACK = 2 !< the checkerboard requires two sweeps

      integer :: n, j, k, i1, j1, k1, id, jd, kd
      integer :: nsmoo

      if (lev == level_min) then
         nsmoo = nsmoob
      else
         nsmoo = nsmool
      endif

      do n = 1, RED_BLACK*nsmoo
         call mpi_multigrid_bnd(lev, soln, 1, extbnd_antimirror) ! no corners are required here

         if (dirty_debug) then
            write(dirty_label, '(a,i5)')"relax soln- smoo=", n
            call check_dirty(lev, soln, dirty_label)
         endif

         ! Possible optimization: this is the most costly part of the RBGS relaxation (instruction count, read and write data, L1 and L2 read cache miss)
         ! do n = 1, nsmoo
         !    call mpi_multigrid_bnd(lev, soln, 1, extbnd_antimirror)
         !    relax single layer of red cells at all faces
         !    call mpi_multigrid_bnd(lev, soln, 1, extbnd_antimirror)
         !    relax interior cells (except for single layer of cells at all faces), first red, then 1-cell behind black one.
         !    relax single layer of black cells at all faces
         ! enddo

         ! with explicit outer loops it is easier to describe a 3-D checkerboard :-)

         if (eff_dim==NDIM .and. .not. multidim_code_3D) then
            do k = lvl(lev)%ks, lvl(lev)%ke
               do j = lvl(lev)%js, lvl(lev)%je
                  i1 = lvl(lev)%is + mod(n+j+k, RED_BLACK)
                  lvl(          lev          )%mgvar(i1  :lvl(lev)%ie  :2, j,   k,   soln) = &
                       lvl(lev)%rx * (lvl(lev)%mgvar(i1-1:lvl(lev)%ie-1:2, j,   k,   soln) + lvl(lev)%mgvar(i1+1:lvl(lev)%ie+1:2, j,   k,   soln)) + &
                       lvl(lev)%ry * (lvl(lev)%mgvar(i1  :lvl(lev)%ie  :2, j-1, k,   soln) + lvl(lev)%mgvar(i1:  lvl(lev)%ie:  2, j+1, k,   soln)) + &
                       lvl(lev)%rz * (lvl(lev)%mgvar(i1  :lvl(lev)%ie  :2, j,   k-1, soln) + lvl(lev)%mgvar(i1:  lvl(lev)%ie:  2, j,   k+1, soln)) - &
                       lvl(lev)%r  *  lvl(lev)%mgvar(i1  :lvl(lev)%ie  :2, j,   k,   src)
               enddo
            enddo
         else
            ! In 3D this variant significantly increases instruction count and also some data read
            i1 = lvl(lev)%is; id = 1 ! mv to multigridvars, init_multigrid
            j1 = lvl(lev)%js; jd = 1
            k1 = lvl(lev)%ks; kd = 1
            if (has_dir(xdim)) then
               id = RED_BLACK
            else if (has_dir(ydim)) then
               jd = RED_BLACK
            else if (has_dir(zdim)) then
               kd = RED_BLACK
            endif

            if (kd == RED_BLACK) k1 = lvl(lev)%ks + mod(n, RED_BLACK)
            do k = k1, lvl(lev)%ke, kd
               if (jd == RED_BLACK) j1 = lvl(lev)%js + mod(n+k, RED_BLACK)
               do j = j1, lvl(lev)%je, jd
                  if (id == RED_BLACK) i1 = lvl(lev)%is + mod(n+j+k, RED_BLACK)
                  lvl(      lev)%mgvar(i1  :lvl(lev)%ie  :id, j,   k,   soln) = &
                       & (1. - Jacobi_damp)* lvl(lev)%mgvar(i1  :lvl(lev)%ie  :id, j,   k,   soln) &
                       &     - Jacobi_damp * lvl(lev)%mgvar(i1  :lvl(lev)%ie  :id, j,   k,   src)  * lvl(lev)%r
                  if (has_dir(xdim)) &
                       lvl (lev)%mgvar(i1  :lvl(lev)%ie  :id, j,   k,   soln) = lvl(lev)%mgvar(i1:  lvl(lev)%ie:  id, j,   k,   soln)  + &
                       &       Jacobi_damp *(lvl(lev)%mgvar(i1-1:lvl(lev)%ie-1:id, j,   k,   soln) + lvl(lev)%mgvar(i1+1:lvl(lev)%ie+1:id, j,   k,   soln)) * lvl(lev)%rx
                  if (has_dir(ydim)) &
                       lvl (lev)%mgvar(i1  :lvl(lev)%ie  :id, j,   k,   soln) = lvl(lev)%mgvar(i1:  lvl(lev)%ie:  id, j,   k,   soln)  + &
                       &       Jacobi_damp *(lvl(lev)%mgvar(i1  :lvl(lev)%ie  :id, j-1, k,   soln) + lvl(lev)%mgvar(i1:  lvl(lev)%ie:  id, j+1, k,   soln)) * lvl(lev)%ry
                  if (has_dir(zdim)) &
                       lvl (lev)%mgvar(i1  :lvl(lev)%ie  :id, j,   k,   soln) = lvl(lev)%mgvar(i1:  lvl(lev)%ie:  id, j,   k,   soln)  + &
                       &       Jacobi_damp *(lvl(lev)%mgvar(i1  :lvl(lev)%ie  :id, j,   k-1, soln) + lvl(lev)%mgvar(i1:  lvl(lev)%ie:  id, j,   k+1, soln)) * lvl(lev)%rz
               enddo
            enddo
         endif

         if (dirty_debug) then
            write(dirty_label, '(a,i5)')"relax soln+ smoo=", n
            call check_dirty(lev, soln, dirty_label)
         endif

      enddo

   end subroutine approximate_solution_rbgs

!!$ ============================================================================
!>
!! \brief FFT given-boundary Poisson solver applied to local domain. Should require less communication than RBGS implementation.
!!
!! \todo test a configuration with wider area being subjected to FFT (sizes would no longer be 2**n) to avoid the need of relaxation
!<

   subroutine approximate_solution_fft(lev, src, soln)

      use grid,               only: D_x, D_y, D_z
      use mpisetup,           only: xdim, ydim, zdim, has_dir, eff_dim
      use dataio_pub,         only: die, warn
      use multigridhelpers,   only: dirty_debug, check_dirty, dirtyL, multidim_code_3D
      use multigridmpifuncs,  only: mpi_multigrid_bnd
      use multigridvars,      only: lvl, LOW, HIGH, NDIM, extbnd_antimirror

      implicit none

      integer, intent(in) :: lev  !< level for which approximate the solution
      integer, intent(in) :: src  !< index of source in lvl()%mgvar
      integer, intent(in) :: soln !< index of solution in lvl()%mgvar

      integer :: nf, n

      do nf = 1, nsmoof
         lvl(lev)%src(:, :, :) = lvl(lev)%mgvar(lvl(lev)%is:lvl(lev)%ie, lvl(lev)%js:lvl(lev)%je, lvl(lev)%ks:lvl(lev)%ke, src)

         if (lvl(lev)%fft_type == fft_dst) then !correct boundaries on non-periodic local domain
            if (nf == 1) then
               call make_face_boundaries(lev, soln)
            else
               call mpi_multigrid_bnd(lev, soln, 1, extbnd_antimirror)
               if (has_dir(xdim)) then
                  lvl(lev)%bnd_x(:, :, LOW)  = 0.5* sum (lvl(lev)%mgvar(lvl(lev)%is-1:lvl(lev)%is, lvl(lev)%js:lvl(lev)%je, lvl(lev)%ks:lvl(lev)%ke, soln), 1)
                  lvl(lev)%bnd_x(:, :, HIGH) = 0.5* sum (lvl(lev)%mgvar(lvl(lev)%ie:lvl(lev)%ie+1, lvl(lev)%js:lvl(lev)%je, lvl(lev)%ks:lvl(lev)%ke, soln), 1)
               endif
               if (has_dir(ydim)) then
                  lvl(lev)%bnd_y(:, :, LOW)  = 0.5* sum (lvl(lev)%mgvar(lvl(lev)%is:lvl(lev)%ie, lvl(lev)%js-1:lvl(lev)%js, lvl(lev)%ks:lvl(lev)%ke, soln), 2)
                  lvl(lev)%bnd_y(:, :, HIGH) = 0.5* sum (lvl(lev)%mgvar(lvl(lev)%is:lvl(lev)%ie, lvl(lev)%je:lvl(lev)%je+1, lvl(lev)%ks:lvl(lev)%ke, soln), 2)
               endif
               if (has_dir(zdim)) then
                  lvl(lev)%bnd_z(:, :, LOW)  = 0.5* sum (lvl(lev)%mgvar(lvl(lev)%is:lvl(lev)%ie, lvl(lev)%js:lvl(lev)%je, lvl(lev)%ks-1:lvl(lev)%ks, soln), 3)
                  lvl(lev)%bnd_z(:, :, HIGH) = 0.5* sum (lvl(lev)%mgvar(lvl(lev)%is:lvl(lev)%ie, lvl(lev)%js:lvl(lev)%je, lvl(lev)%ke:lvl(lev)%ke+1, soln), 3)
               endif
            endif

            if (dirty_debug) then
               if (has_dir(xdim) .and. any(abs(lvl(lev)%bnd_x(:, :, :)) > dirtyL)) call warn("approximate_solution_fft dirty bnd_x")
               if (has_dir(ydim) .and. any(abs(lvl(lev)%bnd_y(:, :, :)) > dirtyL)) call warn("approximate_solution_fft dirty bnd_y")
               if (has_dir(zdim) .and. any(abs(lvl(lev)%bnd_z(:, :, :)) > dirtyL)) call warn("approximate_solution_fft dirty bnd_z")
            endif

            if (has_dir(xdim)) then
               lvl(lev)%src(1,            :, :) = lvl(lev)%src(1,            :, :) - lvl(lev)%bnd_x(:, :, LOW)  * 2. * lvl(lev)%idx2
               lvl(lev)%src(lvl(lev)%nxb, :, :) = lvl(lev)%src(lvl(lev)%nxb, :, :) - lvl(lev)%bnd_x(:, :, HIGH) * 2. * lvl(lev)%idx2
            endif
            if (has_dir(ydim)) then
               lvl(lev)%src(:, 1,            :) = lvl(lev)%src(:, 1,            :) - lvl(lev)%bnd_y(:, :, LOW)  * 2. * lvl(lev)%idy2
               lvl(lev)%src(:, lvl(lev)%nyb, :) = lvl(lev)%src(:, lvl(lev)%nyb, :) - lvl(lev)%bnd_y(:, :, HIGH) * 2. * lvl(lev)%idy2
            endif
            if (has_dir(zdim)) then
               lvl(lev)%src(:, :, 1           ) = lvl(lev)%src(:, :, 1           ) - lvl(lev)%bnd_z(:, :, LOW)  * 2. * lvl(lev)%idz2
               lvl(lev)%src(:, :, lvl(lev)%nzb) = lvl(lev)%src(:, :, lvl(lev)%nzb) - lvl(lev)%bnd_z(:, :, HIGH) * 2. * lvl(lev)%idz2
            endif
         endif

         call fft_convolve(lev)

         lvl(lev)%mgvar(lvl(lev)%is:lvl(lev)%ie, lvl(lev)%js:lvl(lev)%je, lvl(lev)%ks:lvl(lev)%ke, soln) = lvl(lev)%src(:, :, :)

         call check_dirty(lev, soln, "approx_soln fft+")

         !> \deprecated BEWARE use has_dir() here in a way that does not degrade performance

         !relax the boundaries
         do n = 1, nsmool
            call mpi_multigrid_bnd(lev, soln, 1, extbnd_antimirror)
            ! Possible optimization: This is a quite costly part of the local FFT solver
            if (fft_full_relax) then
               if (eff_dim == NDIM .and. .not. multidim_code_3D) then
                  lvl                    (lev)%mgvar(lvl(lev)%is:lvl(lev)%ie,     lvl(lev)%js:lvl(lev)%je,     lvl(lev)%ks:lvl(lev)%ke,     soln)  = &
                       lvl(lev)%rx * (lvl(lev)%mgvar(lvl(lev)%is-1:lvl(lev)%ie-1, lvl(lev)%js:lvl(lev)%je,     lvl(lev)%ks:lvl(lev)%ke,     soln)  + &
                       &              lvl(lev)%mgvar(lvl(lev)%is+1:lvl(lev)%ie+1, lvl(lev)%js:lvl(lev)%je,     lvl(lev)%ks:lvl(lev)%ke,     soln)) + &
                       lvl(lev)%ry * (lvl(lev)%mgvar(lvl(lev)%is:lvl(lev)%ie,     lvl(lev)%js-1:lvl(lev)%je-1, lvl(lev)%ks:lvl(lev)%ke,     soln)  + &
                       &              lvl(lev)%mgvar(lvl(lev)%is:lvl(lev)%ie,     lvl(lev)%js+1:lvl(lev)%je+1, lvl(lev)%ks:lvl(lev)%ke,     soln)) + &
                       lvl(lev)%rz * (lvl(lev)%mgvar(lvl(lev)%is:lvl(lev)%ie,     lvl(lev)%js:lvl(lev)%je,     lvl(lev)%ks-1:lvl(lev)%ke-1, soln)  + &
                       &              lvl(lev)%mgvar(lvl(lev)%is:lvl(lev)%ie,     lvl(lev)%js:lvl(lev)%je,     lvl(lev)%ks+1:lvl(lev)%ke+1, soln)) - &
                       lvl(lev)%r  *  lvl(lev)%mgvar(lvl(lev)%is:lvl(lev)%ie,     lvl(lev)%js:lvl(lev)%je,     lvl(lev)%ks:lvl(lev)%ke,     src)
               else

                  call die("[multigrid_gravity:approximate_solution_fft] fft_full_relax is allowed only for 3D at the moment")

                  ! An additional array (lvl(lev)%src would be good enough) is required here to assemble partial results or use red-black passes
                  lvl                    (lev)%mgvar(lvl(lev)%is:lvl(lev)%ie,     lvl(lev)%js:lvl(lev)%je,     lvl(lev)%ks:lvl(lev)%ke,     soln)  = &
                       - lvl(lev)%r * lvl(lev)%mgvar(lvl(lev)%is:lvl(lev)%ie,     lvl(lev)%js:lvl(lev)%je,     lvl(lev)%ks:lvl(lev)%ke,     src)
                  if (has_dir(xdim)) &
                       lvl               (lev)%mgvar(lvl(lev)%is:lvl(lev)%ie,     lvl(lev)%js:lvl(lev)%je,     lvl(lev)%ks:lvl(lev)%ke,     soln)  = &
                       lvl               (lev)%mgvar(lvl(lev)%is:lvl(lev)%ie,     lvl(lev)%js:lvl(lev)%je,     lvl(lev)%ks:lvl(lev)%ke,     soln)  + &
                       lvl(lev)%rx * (lvl(lev)%mgvar(lvl(lev)%is-1:lvl(lev)%ie-1, lvl(lev)%js:lvl(lev)%je,     lvl(lev)%ks:lvl(lev)%ke,     soln)  + &
                       &              lvl(lev)%mgvar(lvl(lev)%is+1:lvl(lev)%ie+1, lvl(lev)%js:lvl(lev)%je,     lvl(lev)%ks:lvl(lev)%ke,     soln))
                  if (has_dir(ydim)) &
                       lvl               (lev)%mgvar(lvl(lev)%is:lvl(lev)%ie,     lvl(lev)%js:lvl(lev)%je,     lvl(lev)%ks:lvl(lev)%ke,     soln)  = &
                       lvl               (lev)%mgvar(lvl(lev)%is:lvl(lev)%ie,     lvl(lev)%js:lvl(lev)%je,     lvl(lev)%ks:lvl(lev)%ke,     soln)  + &
                       lvl(lev)%ry * (lvl(lev)%mgvar(lvl(lev)%is:lvl(lev)%ie,     lvl(lev)%js-1:lvl(lev)%je-1, lvl(lev)%ks:lvl(lev)%ke,     soln)  + &
                       &              lvl(lev)%mgvar(lvl(lev)%is:lvl(lev)%ie,     lvl(lev)%js+1:lvl(lev)%je+1, lvl(lev)%ks:lvl(lev)%ke,     soln))
                  if (has_dir(zdim)) &
                       lvl               (lev)%mgvar(lvl(lev)%is:lvl(lev)%ie,     lvl(lev)%js:lvl(lev)%je,     lvl(lev)%ks:lvl(lev)%ke,     soln)  = &
                       lvl               (lev)%mgvar(lvl(lev)%is:lvl(lev)%ie,     lvl(lev)%js:lvl(lev)%je,     lvl(lev)%ks:lvl(lev)%ke,     soln)  + &
                       lvl(lev)%rz * (lvl(lev)%mgvar(lvl(lev)%is:lvl(lev)%ie,     lvl(lev)%js:lvl(lev)%je,     lvl(lev)%ks-1:lvl(lev)%ke-1, soln)  + &
                       &              lvl(lev)%mgvar(lvl(lev)%is:lvl(lev)%ie,     lvl(lev)%js:lvl(lev)%je,     lvl(lev)%ks+1:lvl(lev)%ke+1, soln))
               endif
            else
               ! relax only two layers of cells (1 is  significantly worse, 3 does not improve much)
               ! edges are relaxed twice, corners are relaxed three times which seems to be good

               if (has_dir(xdim)) then
                  lvl                    (lev)%mgvar(lvl(lev)%is    :lvl(lev)%is+D_x,   lvl(lev)%js    :lvl(lev)%je,     lvl(lev)%ks    :lvl(lev)%ke,     soln)  = & ! -X
                       lvl(lev)%rx * (lvl(lev)%mgvar(lvl(lev)%is-D_x:lvl(lev)%is,       lvl(lev)%js    :lvl(lev)%je,     lvl(lev)%ks    :lvl(lev)%ke,     soln)  + &
                       &              lvl(lev)%mgvar(lvl(lev)%is+D_x:lvl(lev)%is+2*D_x, lvl(lev)%js    :lvl(lev)%je,     lvl(lev)%ks    :lvl(lev)%ke,     soln)) + &
                       lvl(lev)%ry * (lvl(lev)%mgvar(lvl(lev)%is    :lvl(lev)%is+D_x,   lvl(lev)%js-D_y:lvl(lev)%je-D_y, lvl(lev)%ks    :lvl(lev)%ke,     soln)  + &
                       &              lvl(lev)%mgvar(lvl(lev)%is    :lvl(lev)%is+D_x,   lvl(lev)%js+D_y:lvl(lev)%je+D_y, lvl(lev)%ks    :lvl(lev)%ke,     soln)) + &
                       lvl(lev)%rz * (lvl(lev)%mgvar(lvl(lev)%is    :lvl(lev)%is+D_x,   lvl(lev)%js    :lvl(lev)%je,     lvl(lev)%ks-D_z:lvl(lev)%ke-D_z, soln)  + &
                       &              lvl(lev)%mgvar(lvl(lev)%is    :lvl(lev)%is+D_x,   lvl(lev)%js    :lvl(lev)%je,     lvl(lev)%ks+D_z:lvl(lev)%ke+D_z, soln)) - &
                       lvl(lev)%r  *  lvl(lev)%mgvar(lvl(lev)%is    :lvl(lev)%is+D_x,   lvl(lev)%js    :lvl(lev)%je,     lvl(lev)%ks    :lvl(lev)%ke,     src)

                  lvl                    (lev)%mgvar(lvl(lev)%ie-D_x  :lvl(lev)%ie,     lvl(lev)%js    :lvl(lev)%je,     lvl(lev)%ks    :lvl(lev)%ke,     soln)  = & ! +X
                       lvl(lev)%rx * (lvl(lev)%mgvar(lvl(lev)%ie-2*D_x:lvl(lev)%ie-D_x, lvl(lev)%js    :lvl(lev)%je,     lvl(lev)%ks    :lvl(lev)%ke,     soln)  + &
                       &              lvl(lev)%mgvar(lvl(lev)%ie      :lvl(lev)%ie+D_x, lvl(lev)%js    :lvl(lev)%je,     lvl(lev)%ks    :lvl(lev)%ke,     soln)) + &
                       lvl(lev)%ry * (lvl(lev)%mgvar(lvl(lev)%ie-D_x  :lvl(lev)%ie,     lvl(lev)%js-D_y:lvl(lev)%je-D_y, lvl(lev)%ks    :lvl(lev)%ke,     soln)  + &
                       &              lvl(lev)%mgvar(lvl(lev)%ie-D_x  :lvl(lev)%ie,     lvl(lev)%js+D_y:lvl(lev)%je+D_y, lvl(lev)%ks    :lvl(lev)%ke,     soln)) + &
                       lvl(lev)%rz * (lvl(lev)%mgvar(lvl(lev)%ie-D_x  :lvl(lev)%ie,     lvl(lev)%js    :lvl(lev)%je,     lvl(lev)%ks-D_z:lvl(lev)%ke-D_z, soln)  + &
                       &              lvl(lev)%mgvar(lvl(lev)%ie-D_x  :lvl(lev)%ie,     lvl(lev)%js    :lvl(lev)%je,     lvl(lev)%ks+D_z:lvl(lev)%ke+D_z, soln)) - &
                       lvl(lev)%r  *  lvl(lev)%mgvar(lvl(lev)%ie-D_x  :lvl(lev)%ie,     lvl(lev)%js    :lvl(lev)%je,     lvl(lev)%ks    :lvl(lev)%ke,     src)
               endif

               if (has_dir(ydim)) then
                  lvl                    (lev)%mgvar(lvl(lev)%is    :lvl(lev)%ie,     lvl(lev)%js    :lvl(lev)%js+D_y,   lvl(lev)%ks    :lvl(lev)%ke,     soln)  = & ! -Y
                       lvl(lev)%rx * (lvl(lev)%mgvar(lvl(lev)%is-D_x:lvl(lev)%ie-D_x, lvl(lev)%js    :lvl(lev)%js+D_y,   lvl(lev)%ks    :lvl(lev)%ke,     soln)  + &
                       &              lvl(lev)%mgvar(lvl(lev)%is+D_x:lvl(lev)%ie+D_x, lvl(lev)%js    :lvl(lev)%js+D_y,   lvl(lev)%ks    :lvl(lev)%ke,     soln)) + &
                       lvl(lev)%ry * (lvl(lev)%mgvar(lvl(lev)%is    :lvl(lev)%ie,     lvl(lev)%js-D_y:lvl(lev)%js,       lvl(lev)%ks    :lvl(lev)%ke,     soln)  + &
                       &              lvl(lev)%mgvar(lvl(lev)%is    :lvl(lev)%ie,     lvl(lev)%js+D_y:lvl(lev)%js+2*D_y, lvl(lev)%ks    :lvl(lev)%ke,     soln)) + &
                       lvl(lev)%rz * (lvl(lev)%mgvar(lvl(lev)%is    :lvl(lev)%ie,     lvl(lev)%js    :lvl(lev)%js+D_y,   lvl(lev)%ks-D_z:lvl(lev)%ke-D_z, soln)  + &
                       &              lvl(lev)%mgvar(lvl(lev)%is    :lvl(lev)%ie,     lvl(lev)%js    :lvl(lev)%js+D_y,   lvl(lev)%ks+D_z:lvl(lev)%ke+D_z, soln)) - &
                       lvl(lev)%r  *  lvl(lev)%mgvar(lvl(lev)%is    :lvl(lev)%ie,     lvl(lev)%js    :lvl(lev)%js+D_y,   lvl(lev)%ks    :lvl(lev)%ke,     src)

                  lvl                    (lev)%mgvar(lvl(lev)%is    :lvl(lev)%ie,     lvl(lev)%je-D_y  :lvl(lev)%je,     lvl(lev)%ks    :lvl(lev)%ke,     soln)  = & ! +Y
                       lvl(lev)%rx * (lvl(lev)%mgvar(lvl(lev)%is-D_x:lvl(lev)%ie-D_x, lvl(lev)%je-D_y  :lvl(lev)%je,     lvl(lev)%ks    :lvl(lev)%ke,     soln)  + &
                       &              lvl(lev)%mgvar(lvl(lev)%is+D_x:lvl(lev)%ie+D_x, lvl(lev)%je-D_y  :lvl(lev)%je,     lvl(lev)%ks    :lvl(lev)%ke,     soln)) + &
                       lvl(lev)%ry * (lvl(lev)%mgvar(lvl(lev)%is    :lvl(lev)%ie,     lvl(lev)%je-2*D_y:lvl(lev)%je-D_y, lvl(lev)%ks    :lvl(lev)%ke,     soln)  + &
                       &              lvl(lev)%mgvar(lvl(lev)%is    :lvl(lev)%ie,     lvl(lev)%je      :lvl(lev)%je+D_y, lvl(lev)%ks    :lvl(lev)%ke,     soln)) + &
                       lvl(lev)%rz * (lvl(lev)%mgvar(lvl(lev)%is    :lvl(lev)%ie,     lvl(lev)%je-D_y  :lvl(lev)%je,     lvl(lev)%ks-D_z:lvl(lev)%ke-D_z, soln)  + &
                       &              lvl(lev)%mgvar(lvl(lev)%is    :lvl(lev)%ie,     lvl(lev)%je-D_y  :lvl(lev)%je,     lvl(lev)%ks+D_z:lvl(lev)%ke+D_z, soln)) - &
                       lvl(lev)%r  *  lvl(lev)%mgvar(lvl(lev)%is    :lvl(lev)%ie,     lvl(lev)%je-D_y  :lvl(lev)%je,     lvl(lev)%ks    :lvl(lev)%ke,     src)
               endif

               if (has_dir(zdim)) then
                  lvl                    (lev)%mgvar(lvl(lev)%is    :lvl(lev)%ie,     lvl(lev)%js    :lvl(lev)%je,     lvl(lev)%ks    :lvl(lev)%ks+D_z,   soln)  = & ! -Z
                       lvl(lev)%rx * (lvl(lev)%mgvar(lvl(lev)%is-D_x:lvl(lev)%ie-D_x, lvl(lev)%js    :lvl(lev)%je,     lvl(lev)%ks    :lvl(lev)%ks+D_z,   soln)  + &
                       &              lvl(lev)%mgvar(lvl(lev)%is+D_x:lvl(lev)%ie+D_x, lvl(lev)%js    :lvl(lev)%je,     lvl(lev)%ks    :lvl(lev)%ks+D_z,   soln)) + &
                       lvl(lev)%ry * (lvl(lev)%mgvar(lvl(lev)%is    :lvl(lev)%ie,     lvl(lev)%js-D_y:lvl(lev)%je-D_y, lvl(lev)%ks    :lvl(lev)%ks+D_z,   soln)  + &
                       &              lvl(lev)%mgvar(lvl(lev)%is    :lvl(lev)%ie,     lvl(lev)%js+D_y:lvl(lev)%je+D_y, lvl(lev)%ks    :lvl(lev)%ks+D_z,   soln)) + &
                       lvl(lev)%rz * (lvl(lev)%mgvar(lvl(lev)%is    :lvl(lev)%ie,     lvl(lev)%js    :lvl(lev)%je,     lvl(lev)%ks-D_z:lvl(lev)%ks,       soln)  + &
                       &              lvl(lev)%mgvar(lvl(lev)%is    :lvl(lev)%ie,     lvl(lev)%js    :lvl(lev)%je,     lvl(lev)%ks+D_z:lvl(lev)%ks+2*D_z, soln)) - &
                       lvl(lev)%r  *  lvl(lev)%mgvar(lvl(lev)%is    :lvl(lev)%ie,     lvl(lev)%js    :lvl(lev)%je,     lvl(lev)%ks    :lvl(lev)%ks+D_z,   src)

                  lvl                    (lev)%mgvar(lvl(lev)%is    :lvl(lev)%ie,     lvl(lev)%js    :lvl(lev)%je,     lvl(lev)%ke-D_z  :lvl(lev)%ke ,    soln)  = & ! +Z
                       lvl(lev)%rx * (lvl(lev)%mgvar(lvl(lev)%is-D_x:lvl(lev)%ie-D_x, lvl(lev)%js    :lvl(lev)%je,     lvl(lev)%ke-D_z  :lvl(lev)%ke,     soln)  + &
                       &              lvl(lev)%mgvar(lvl(lev)%is+D_x:lvl(lev)%ie+D_x, lvl(lev)%js    :lvl(lev)%je,     lvl(lev)%ke-D_z  :lvl(lev)%ke,     soln)) + &
                       lvl(lev)%ry * (lvl(lev)%mgvar(lvl(lev)%is    :lvl(lev)%ie,     lvl(lev)%js-D_y:lvl(lev)%je-D_y, lvl(lev)%ke-D_z  :lvl(lev)%ke,     soln)  + &
                       &              lvl(lev)%mgvar(lvl(lev)%is    :lvl(lev)%ie,     lvl(lev)%js+D_y:lvl(lev)%je+D_y, lvl(lev)%ke-D_z  :lvl(lev)%ke,     soln)) + &
                       lvl(lev)%rz * (lvl(lev)%mgvar(lvl(lev)%is    :lvl(lev)%ie,     lvl(lev)%js    :lvl(lev)%je,     lvl(lev)%ke-2*D_z:lvl(lev)%ke-D_z, soln)  + &
                       &              lvl(lev)%mgvar(lvl(lev)%is    :lvl(lev)%ie,     lvl(lev)%js    :lvl(lev)%je,     lvl(lev)%ke      :lvl(lev)%ke+D_z, soln)) - &
                       lvl(lev)%r  *  lvl(lev)%mgvar(lvl(lev)%is    :lvl(lev)%ie,     lvl(lev)%js    :lvl(lev)%je,     lvl(lev)%ke-D_z  :lvl(lev)%ke,     src)
               endif

            endif
         enddo

         call check_dirty(lev, soln, "approx_soln relax+")

      enddo

   end subroutine approximate_solution_fft

!!$ ============================================================================
!>
!! \brief This routine prepares boundary values for local-FFT solver
!<

   subroutine make_face_boundaries(lev, soln)

      use mpisetup,           only: nproc
      use multigridbasefuncs, only: zero_boundaries, prolong_faces
      use dataio_pub,         only: warn
      use multigridvars,      only: bnd_givenval, bnd_periodic, level_min

      implicit none

      integer, intent(in) :: lev  !< level for which approximate the solution
      integer, intent(in) :: soln !< index of solution in lvl()%mgvar

      if (grav_bnd == bnd_periodic .and. nproc == 1) then
         call zero_boundaries(lev)
      else
         if (lev > level_min) then
            call prolong_faces(lev, soln)
         else
            if (grav_bnd /= bnd_givenval) call zero_boundaries(lev)
            call warn("m:mfb WTF?")
         endif
      endif

   end subroutine make_face_boundaries

!!$ ============================================================================
!>
!! \brief Obtain an exact solution on base level. (wrapper)
!<

   subroutine gb_fft_solve(src, soln)

      implicit none

      integer, intent(in) :: src  !< index of source in lvl()%mgvar
      integer, intent(in) :: soln !< index of solution in lvl()%mgvar

      if (gb_solve_gather) then
         call gb_fft_solve_gather(src, soln)
      else
         call gb_fft_solve_sendrecv(src, soln)
      endif

   end subroutine gb_fft_solve

!!$ ============================================================================
!>
!! \brief Obtain an exact solution on base level.
!! \details The source is gathered on the master PE (possible bottleneck), solution is obtained through FFT, then its parts are sent to all PEs.
!! Unfortunately non-blocking communication probably will not change much here.
!! Alternative implementation: set gb_src to 0, copy local part to it, then call MPI_Allreduce and solve with FFT on each PE (no need to communicate the solution)
!<

   subroutine gb_fft_solve_sendrecv(src, soln)

      use mpisetup,      only: nproc, proc, master, ierr, comm3d, status, xdim, ydim, zdim
      use mpi,           only: MPI_DOUBLE_PRECISION
      use multigridvars, only: gb, gb_cartmap, base

      implicit none

      integer, intent(in) :: src  !< index of source in lvl()%mgvar
      integer, intent(in) :: soln !< index of solution in lvl()%mgvar

      integer :: p

      if (master) then

         ! collect the source on the base level
         gb%src(gb_cartmap(0)%lo(xdim):gb_cartmap(0)%up(xdim), &
              & gb_cartmap(0)%lo(ydim):gb_cartmap(0)%up(ydim), &
              & gb_cartmap(0)%lo(zdim):gb_cartmap(0)%up(zdim)) = &
              base%mgvar(base%is:base%ie, base%js:base%je, base%ks:base%ke, src)
         do p = 1, nproc -1
            call MPI_Recv(gb%src(gb_cartmap(p)%lo(xdim):gb_cartmap(p)%up(xdim), &
                 &               gb_cartmap(p)%lo(ydim):gb_cartmap(p)%up(ydim), &
                 &               gb_cartmap(p)%lo(zdim):gb_cartmap(p)%up(zdim)), &
                 &        base%nxb*base%nyb*base%nzb, MPI_DOUBLE_PRECISION, p, p, comm3d, status, ierr)
         enddo

         call fft_convolve(gb%level)

         ! send all the pieces away
         base%mgvar(base%is:base%ie, base%js:base%je, base%ks:base%ke, soln) = &
              gb%src(gb_cartmap(0)%lo(xdim):gb_cartmap(0)%up(xdim), &
              &      gb_cartmap(0)%lo(ydim):gb_cartmap(0)%up(ydim), &
              &      gb_cartmap(0)%lo(zdim):gb_cartmap(0)%up(zdim))
         do p = 1, nproc -1
            call MPI_Send(gb%src(gb_cartmap(p)%lo(xdim):gb_cartmap(p)%up(xdim), &
                 &               gb_cartmap(p)%lo(ydim):gb_cartmap(p)%up(ydim), &
                 &               gb_cartmap(p)%lo(zdim):gb_cartmap(p)%up(zdim)), &
                 &        base%nxb*base%nyb*base%nzb, MPI_DOUBLE_PRECISION, p, p, comm3d, ierr)
         enddo

      else
         call MPI_Send(base%mgvar(base%is:base%ie, base%js:base%je, base%ks:base%ke, src),  base%nxb*base%nyb*base%nzb, MPI_DOUBLE_PRECISION, 0, proc, comm3d, ierr)
         call MPI_Recv(base%mgvar(base%is:base%ie, base%js:base%je, base%ks:base%ke, soln), base%nxb*base%nyb*base%nzb, MPI_DOUBLE_PRECISION, 0, proc, comm3d, status, ierr)
      endif

   end subroutine gb_fft_solve_sendrecv

!!$ ============================================================================
!>
!! \brief Alternative version using MPI_Gather
!!
!! \todo using 1d ffts and gather to reduce dimensions one by one could prevent possible bottleneck and could be (?) faster:
!!       1. N_x * N_y independent gathers reduce Z dim
!!       2. fft in Z dim on N_x*N_y proc
!!       3. N_x independent gathers of (2.) reduce X dim
!!       4. fft in Y dim on N_x
!!       5. gather (4.) on master proc
!!       6. fft in X dim on master
!!       7. Scatter
!<

   subroutine gb_fft_solve_gather(src, soln)

      use mpisetup,      only: nproc, master, ierr, comm3d, xdim, ydim, zdim
      use mpi,           only: MPI_DOUBLE_PRECISION
      use multigridvars, only: gb, gb_cartmap, base

      implicit none

      integer, intent(in) :: src  !< index of source in lvl()%mgvar
      integer, intent(in) :: soln !< index of solution in lvl()%mgvar

      integer :: p

      call MPI_Gather(base%mgvar(base%is:base%ie, base%js:base%je, base%ks:base%ke, src),  base%nxb*base%nyb*base%nzb, MPI_DOUBLE_PRECISION, &
           &          gb_src_temp, base%nxb*base%nyb*base%nzb, MPI_DOUBLE_PRECISION, 0, comm3d, ierr)

      if (master) then
         do p = 0, nproc-1
            gb%src(gb_cartmap(p)%lo(xdim):gb_cartmap(p)%up(xdim), &
                   gb_cartmap(p)%lo(ydim):gb_cartmap(p)%up(ydim), &
                   gb_cartmap(p)%lo(zdim):gb_cartmap(p)%up(zdim)) = gb_src_temp(:,:,:,p)
         enddo

         call fft_convolve(gb%level)

         do p = 0, nproc-1
            gb_src_temp(:,:,:,p) = &
               gb%src(gb_cartmap(p)%lo(xdim):gb_cartmap(p)%up(xdim), &
                      gb_cartmap(p)%lo(ydim):gb_cartmap(p)%up(ydim), &
                      gb_cartmap(p)%lo(zdim):gb_cartmap(p)%up(zdim))
         enddo
      endif

      call MPI_Scatter(gb_src_temp,  base%nxb*base%nyb*base%nzb, MPI_DOUBLE_PRECISION, &
                       base%mgvar(base%is:base%ie, base%js:base%je, base%ks:base%ke, soln), base%nxb*base%nyb*base%nzb, MPI_DOUBLE_PRECISION, &
                       0, comm3d, ierr)

   end subroutine gb_fft_solve_gather

!!$ ============================================================================
!>
!! \brief Solve finest level if allowed (typically on single CPU)
!<

   subroutine fft_solve_roof

      use multigridvars,    only: roof, source, solution

      implicit none

      if (roof%fft_type == fft_none) return

      roof%src(:, :, :) = roof%mgvar(roof%is:roof%ie, roof%js:roof%je, roof%ks:roof%ke, source)
      call fft_convolve(roof%level)
      roof%mgvar(roof%is:roof%ie, roof%js:roof%je, roof%ks:roof%ke, solution) = roof%src(:, :, :)

   end subroutine fft_solve_roof

!!$ ============================================================================
!>
!! \brief Do the FFT convolution
!<

   subroutine fft_convolve(level)

      use dataio_pub,    only: die
      use multigridvars, only: lvl

      implicit none

      integer, intent(in) :: level !< level at which make the convolution

      ! do the convolution in Fourier space; lvl(level)%src(:,:,:) -> lvl(level)%fft{r}(:,:,:)
      call dfftw_execute(lvl(level)%planf)

      select case (lvl(level)%fft_type)
         case (fft_rcr)
            lvl(level)%fft  = lvl(level)%fft  * lvl(level)%Green3D
         case (fft_dst)
            lvl(level)%fftr = lvl(level)%fftr * lvl(level)%Green3D
         case default
            call die("[multigrid_gravity:fft_convolve] Unknown FFT type.")
      end select

      call dfftw_execute(lvl(level)%plani) ! lvl(level)%fft{r}(:,:,:) -> lvl(level)%src(:,:,:)

   end subroutine fft_convolve

#else /* !(MULTIGRID && GRAV) */
#warning This should not happen. Probably the multigrid_gravity.F90 file is included in object directory by mistake.
#endif /* !(MULTIGRID && GRAV) */

end module multigrid_gravity
