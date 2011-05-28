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
   character(len=cbuff_len) :: grav_bnd_str                           !< Type of gravitational boundary conditions.

   ! boundaries
   integer                     :: grav_bnd                            !< boundary type for computational domain
   !integer                     :: grav_extbnd_mode                    !< external boundary mode

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
!! <tr><td>grav_bnd_str          </td><td>"periodic"/"dirichlet"</td><td>string of chars</td><td>\copydoc multigrid_gravity::grav_bnd_str          </td></tr>
!! </table>
!! \n \n
!<
   subroutine init_multigrid_grav

      use multigridvars,      only: bnd_periodic, bnd_dirichlet, bnd_isolated, bnd_invalid, correction, mg_nb, ngridvars, single_base
      use constants,          only: GEO_XYZ, GEO_RPZ, BND_PER
      use multipole,          only: use_point_monopole, lmax, mmax, ord_prolong_mpole, coarsen_multipole, interp_pt2mom, interp_mom2pot
      use mpisetup,           only: buffer_dim, comm, comm3d, ierr, master, slave, ibuff, cbuff, rbuff, lbuff, dom, has_dir, eff_dim, geometry_type, is_uneven
      use mpi,                only: MPI_CHARACTER, MPI_DOUBLE_PRECISION, MPI_INTEGER, MPI_LOGICAL, MPI_COMM_NULL
      use dataio_pub,         only: par_file, ierrh, namelist_errh, compare_namelist, cmdl_nml  ! QA_WARN required for diff_nml
      use dataio_pub,         only: msg, die, warn

      implicit none

      integer       :: periodic_bnd_cnt   !< counter of periodic boundaries in existing directions
      logical, save :: frun = .true.      !< First run flag

      namelist /MULTIGRID_GRAVITY/ norm_tol, vcycle_abort, max_cycles, nsmool, nsmoob, &
           &                       overrelax, overrelax_x, overrelax_y, overrelax_z, Jacobi_damp, L4_strength, nsmoof, ord_laplacian, ord_time_extrap, &
           &                       prefer_rbgs_relaxation, base_no_fft, fft_full_relax, fft_patient, trust_fft_solution, &
           &                       coarsen_multipole, lmax, mmax, ord_prolong_mpole, use_point_monopole, interp_pt2mom, interp_mom2pot, &
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
      base_no_fft            = .false.
      prefer_rbgs_relaxation = .false.
      fft_full_relax         = .false.
      fft_patient            = .false.
      interp_pt2mom          = .false.
      interp_mom2pot         = .false.

      periodic_bnd_cnt = count(dom%periodic(:) .and. has_dir(:))

      if (periodic_bnd_cnt == eff_dim) then
         grav_bnd_str = "periodic"
      else
         grav_bnd_str = "dirichlet"
      endif

      if (master) then

         diff_nml(MULTIGRID_GRAVITY)

         ! FIXME when ready
         if (geometry_type == GEO_RPZ) then
            call warn("[multigrid_gravity:init_multigrid_grav] cylindrical geometry support is under development.")
            ! switch off FFT-related bits
            base_no_fft = .true.
            prefer_rbgs_relaxation = .true.
            ord_laplacian = 2
            L4_strength = 0.
            ! ord_prolong_mpole = 0
         else if (geometry_type /= GEO_XYZ) then
            call die("[multigrid_gravity:init_multigrid_grav] non-cartesian geometry not implemented yet.")
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
         base_no_fft             = lbuff(3)
         prefer_rbgs_relaxation  = lbuff(4)
         fft_full_relax          = lbuff(5)
         fft_patient             = lbuff(6)
         interp_pt2mom           = lbuff(7)
         interp_mom2pot          = lbuff(8)

         grav_bnd_str   = cbuff(1)(1:len(grav_bnd_str))

      endif

      ngridvars = max(ngridvars, correction)

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

      if (periodic_bnd_cnt == eff_dim) then ! fully periodic domain
         if (grav_bnd /= bnd_periodic .and. master) call warn("[multigrid_gravity:init_multigrid_grav] Ignoring non-periodic boundary conditions for gravity on a fully periodic domain.")
         grav_bnd = bnd_periodic
      else if (periodic_bnd_cnt > 0 .and. periodic_bnd_cnt < eff_dim) then
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
      if (comm3d /= MPI_COMM_NULL .and. .not. base_no_fft) then
         base_no_fft = .true.
         if (master) call warn("[multigrid_gravity:init_multigrid_grav] comm3d disables use of FFT at base level")
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

      !! Sanity checks
      if (abs(ord_laplacian) > 2*mg_nb) call die("[multigrid_gravity:init_multigrid_grav] not enough guardcells for given laplacian operator order")

   end subroutine init_multigrid_grav

!!$ ============================================================================
!>
!! \brief Initialization - continued after allocation of everything interesting
!<

   subroutine init_multigrid_grav_post(mb_alloc)

      use multigridvars,    only: lvl, plvl, roof, base, bnd_periodic, bnd_dirichlet, bnd_isolated, vcycle_stats, is_mg_uneven, need_general_pf, single_base
      use mpisetup,         only: master, nproc, geometry_type, dom, proc, comm3d
      use multigridhelpers, only: vcycle_stats_init, dirty_debug, dirtyH
      use constants,        only: pi, dpi, GEO_XYZ
      use dataio_pub,       only: die, warn
      use multipole,        only: init_multipole, coarsen_multipole
      use mpi,              only: MPI_COMM_NULL
      use grid,             only: cg

      implicit none

      real, intent(inout)              :: mb_alloc               !< Allocation counter

      type(soln_history), pointer      :: os
      real, allocatable, dimension(:)  :: kx, ky, kz             !< FFT kernel directional components for convolution
      integer, dimension(6)            :: aerr                   !> \deprecated BEWARE: hardcoded magic integer. Update when you change number of simultaneous error checks
      integer :: i, j
      type(plvl), pointer :: curl

      need_general_pf = comm3d == MPI_COMM_NULL .or. single_base .or. is_mg_uneven

      if (need_general_pf .and. coarsen_multipole /= 0) then
         coarsen_multipole = 0
         if (master) call warn("[multigrid_gravity:init_multigrid_grav_post] multipole coarsening on uneven domains or with comm3d == MPI_COMM_NULL is not implemented yet.")
      endif

      call mpi_multigrid_prep_grav !supplement to mpi_multigrid_prep

      ! solution recycling
      ord_time_extrap = min(nold_max-1, max(-1, ord_time_extrap))
      nold = ord_time_extrap + 1
      if (nold > 0) then
         do j = 1, 2 ! inner and outer arrays
            nullify(os)
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
                  os%old(i)%time= -huge(1.0)
                  if (dirty_debug) os%old(i)%soln(:, :, :) = dirtyH
               enddo
               if (any(aerr(1:nold) /= 0)) call die("[multigrid_gravity:init_multigrid_grav_post] Allocation error: os%old(:)%soln")
               os%valid = .false.
               os%last  = 1
            endif
         enddo
      endif

      cg%sgp%arr(:,:,:) = 0. !Initialize all the guardcells, even those which does not impact the solution

      curl => base
      do while (associated(curl))
         ! this should work correctly also when eff_dim < 3
         if (curl%empty) then
            curl%fft_type = fft_none
            curl%r  = 0.
            curl%rx = 0.
            curl%ry = 0.
            curl%rz = 0.
         else
            curl%r  = overrelax   / 2.
            curl%rx = curl%dvol2 * curl%idx2
            curl%ry = curl%dvol2 * curl%idy2
            curl%rz = curl%dvol2 * curl%idz2
            curl%r  = curl%r  / (curl%rx + curl%ry + curl%rz)
            curl%rx = overrelax_x * curl%rx * curl%r
            curl%ry = overrelax_y * curl%ry * curl%r
            curl%rz = overrelax_z * curl%rz * curl%r
            curl%r  = curl%r  * curl%dvol2

            !>
            !! \deprecated BEWARE: some of the above invariants may be not optimally defined - the convergence ratio drops when dx /= dy or dy /= dz or dx /= dz
            !! and overrelaxation factors are required to get any convergence (often poor)
            !<

            if (prefer_rbgs_relaxation) then
               curl%fft_type = fft_none
            else if (grav_bnd == bnd_periodic .and. nproc == 1) then
               curl%fft_type = fft_rcr
            else if (grav_bnd == bnd_periodic .or. grav_bnd == bnd_dirichlet .or. grav_bnd == bnd_isolated) then
               curl%fft_type = fft_dst
            else
               curl%fft_type = fft_none
            endif
         endif
         curl => curl%finer
      enddo

      ! data related to local and global base-level FFT solver
      if (.not. base%empty) then
         if (base_no_fft) then
            base%fft_type = fft_none
         else
            select case (grav_bnd)
               case (bnd_periodic)
                  base%fft_type = fft_rcr
               case (bnd_dirichlet, bnd_isolated)
                  base%fft_type = fft_dst
               case default
                  base%fft_type = fft_none
                  if (master) call warn("[multigrid_gravity:init_multigrid_grav_post] base_no_fft set but no suitable boundary conditions found. Reverting to RBGS relaxation.")
            end select
         endif
      endif

      !special initialization of global base-level FFT-related data
      if (any(lvl(:)%fft_type /= fft_none) .and. geometry_type /= GEO_XYZ) &
           call die("[multigrid_gravity:init_multigrid_grav_post] FFT is not allowed in non-cartesian coordinates.")

      ! FFT solver storage and data
      curl => base
      do while (associated(curl))

         curl%planf = 0
         curl%plani = 0

         if (curl%fft_type /= fft_none .and. .not. curl%empty) then

            if (geometry_type /= GEO_XYZ) call die("[multigrid_gravity:init_multigrid_grav_post] FFT is not allowed in non-cartesian coordinates.")

            select case (curl%fft_type)
               case (fft_rcr)
                  curl%nxc = curl%nxb / 2 + 1
               case (fft_dst)
                  curl%nxc = curl%nxb
               case default
                  call die("[multigrid_gravity:init_multigrid_grav_post] Unknown FFT type.")
            end select

            if (allocated(curl%Green3D) .or. allocated(curl%src)) call die("[multigrid_gravity:init_multigrid_grav_post] Green3D or src arrays already allocated")
            allocate(curl%Green3D(curl%nxc, curl%nyb, curl%nzb), stat=aerr(1))
            allocate(curl%src    (curl%nxb, curl%nyb, curl%nzb), stat=aerr(2))
            if (any(aerr(1:2) /= 0)) call die("[multigrid_gravity:init_multigrid_grav_post] Allocation error: curl%Green3D or curl%src.")
            mb_alloc = mb_alloc + size(curl%Green3D) + size(curl%src)

            if (allocated(kx)) deallocate(kx)
            if (allocated(ky)) deallocate(ky)
            if (allocated(kz)) deallocate(kz)
            allocate(kx(curl%nxc), stat=aerr(1))
            allocate(ky(curl%nyb), stat=aerr(2))
            allocate(kz(curl%nzb), stat=aerr(3))
            if (any(aerr(1:3) /= 0)) call die("[multigrid_gravity:init_multigrid_grav_post] Allocation error: k[xyz]")

            select case (curl%fft_type)

              ! curl%fft_norm is set such that the following sequence gives identity:
              ! call dfftw_execute(curl%planf); curl%fftr(:, :, :) = curl%fftr(:, :, :) * curl%fft_norm ; call dfftw_execute(curl%plani)

               case (fft_rcr)
                  if (allocated(curl%fft)) call die("[multigrid_gravity:init_multigrid_grav_post] fft or Green3D array already allocated")
                  allocate(curl%fft(curl%nxc, curl%nyb, curl%nzb), stat=aerr(1))
                  if (aerr(1) /= 0) call die("[multigrid_gravity:init_multigrid_grav_post] Allocation error: fft.")
                  mb_alloc  = mb_alloc + 2*size(curl%fft)

                  curl%fft_norm = 1. / real( curl%nxb * curl%nyb * curl%nzb ) ! No 4 pi G factor here because the source was already multiplied by it

                  ! FFT local solver initialization for 2nd order (3-point) Laplacian
                  ! sin(k*x-d) - 2.*sin(k*x) + sin(k*x+d) = 2 * (cos(d)-1) * sin(k*x) = -4 * sin(d/2)**2 * sin(k*x)
                  ! For 4th order: a*sin(k*x) + b*(sin(k*x-d) + sin(k*x+d)) + c*(sin(k*x-2*d) + sin(k*x+2*d)), a+2*b+2*c == 0 it would be:
                  ! 4*(a+b+(a+2*b)*cos(d)) * sin(d/2)**2 * sin(k*x)
                  ! For 6th order: a*sin(k*x) + b*(sin(k*x-d) + sin(k*x+d)) + c*(sin(k*x-2*d) + sin(k*x+2*d)) + e*(sin(k*x-3*d) + sin(k*x+3*d)), a+2*b+2*c+2*e == 0 it would be:
                  ! 2*(3*a+4*b+2*c+4*(a+2*b+c)*cos(d)+2*(a+2*(b+c))*cos(2*d)) * sin(d/2)**2 * sin(k*x)
                  ! asymptotically: -d**2/2 for d<pi

                  kx(:) = curl%idx2 * (cos(dpi/curl%nxb*(/(j, j=0, curl%nxc-1)/)) - 1.)
                  ky(:) = curl%idy2 * (cos(dpi/curl%nyb*(/(j, j=0, curl%nyb-1)/)) - 1.)
                  kz(:) = curl%idz2 * (cos(dpi/curl%nzb*(/(j, j=0, curl%nzb-1)/)) - 1.)
                  call dfftw_plan_dft_r2c_3d(curl%planf, curl%nxb, curl%nyb, curl%nzb, curl%src, curl%fft, fftw_flags)
                  call dfftw_plan_dft_c2r_3d(curl%plani, curl%nxb, curl%nyb, curl%nzb, curl%fft, curl%src, fftw_flags)

               case (fft_dst)

                  if (allocated(curl%fftr)) call die("[multigrid_gravity:init_multigrid_grav_post] fftr array already allocated")
                  allocate(curl%fftr(curl%nxc, curl%nyb, curl%nzb), stat=aerr(1))
                  if (aerr(1) /= 0) call die("[multigrid_gravity:init_multigrid_grav_post] Allocation error: fftr.")
                  mb_alloc  = mb_alloc + size(curl%fftr)

                  curl%fft_norm = 1. / (8. * real( curl%nxb * curl%nyb * curl%nzb ))
                  kx(:) = curl%idx2 * (cos(pi/curl%nxb*(/(j, j=1, curl%nxc)/)) - 1.)
                  ky(:) = curl%idy2 * (cos(pi/curl%nyb*(/(j, j=1, curl%nyb)/)) - 1.)
                  kz(:) = curl%idz2 * (cos(pi/curl%nzb*(/(j, j=1, curl%nzb)/)) - 1.)
                  call dfftw_plan_r2r_3d(curl%planf, curl%nxb, curl%nyb, curl%nzb, curl%src,  curl%fftr, FFTW_RODFT10, FFTW_RODFT10, FFTW_RODFT10, fftw_flags)
                  call dfftw_plan_r2r_3d(curl%plani, curl%nxb, curl%nyb, curl%nzb, curl%fftr, curl%src,  FFTW_RODFT01, FFTW_RODFT01, FFTW_RODFT01, fftw_flags)

               case default
                  call die("[multigrid_gravity:init_multigrid_grav_post] Unknown FFT type.")
            end select

            ! compute Green's function for 7-point 3D discrete laplacian
            do i = 1, curl%nxc
               do j = 1, curl%nyb
                  where ( (kx(i) + ky(j) + kz(:)) /= 0 )
                     curl%Green3D(i,j,:) = 0.5 * curl%fft_norm / (kx(i) + ky(j) + kz(:))
                  elsewhere
                     curl%Green3D(i,j,:) = 0.0
                  endwhere
               enddo
            enddo

         endif
         curl => curl%finer
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
!! \brief Set up communication for face prolongation
!! \todo implement also prolongation of coarsened multipoles
!<

   subroutine mpi_multigrid_prep_grav

      use constants,     only: xdim, ydim, zdim, ndims, LO, HI, LONG
      use dataio_pub,    only: warn, die
      use mpisetup,      only: has_dir, proc, nproc, is_neigh, procmask, inflate_req, req
      use multigridvars, only: base, lvl, plvl, pr_segment, ord_prolong_face_norm, is_external, need_general_pf
#ifdef DEBUG
      use piernikdebug,  only: aux_R
#endif /* DEBUG */

      implicit none

      integer :: d, g, j, lh, hl, l, nl
      logical, dimension(xdim:zdim, LO:HI) :: neigh
      logical :: sharing, corner, face
      integer(kind=8), dimension(xdim:zdim) :: ijks, per
      logical, dimension(xdim:zdim) :: dmask
      integer(kind=8), dimension(xdim:zdim, LO:HI) :: coarsened, b_layer
      type(pr_segment), pointer :: seg
      type(plvl), pointer   :: curl                   !> current level (a pointer sliding along the linked list)
      logical :: is_internal_fine
      integer, parameter :: max_opfn = 2
      real, dimension(0:max_opfn) :: opfn_c_ff, opfn_c_cf

      if (.not. need_general_pf) return

      call inflate_req(size([LO, HI]) * 2 * nproc * ndims)

      if (ord_prolong_face_norm > max_opfn) ord_prolong_face_norm = max_opfn
      select case (ord_prolong_face_norm)
         case (0)
            opfn_c_ff(:) = [ 0.5, 0., 0. ]
            opfn_c_cf(:) = [ 1.0, 0., 0. ]
         case (1)
#ifdef DEBUG
            ! the maximum convergence is at aux_R(1) = 0.25 +/- 0.05
            opfn_c_ff(:) = 0.5 * [ 1+aux_R(1), -aux_R(1), 0. ]
            ! the maximum convergence is at aux_R(2) = -0.05 +/- 0.05 and 0
            opfn_c_cf(:) = [ 1.0 + 2*aux_R(2), -aux_R(2), 0. ]
#else /* !DEBUG */
            opfn_c_ff(:) = [  5., -1., 0. ] / 8.  ! adjusted experimentally
            opfn_c_cf(:) = [ 18.,  1., 0. ] / 20. ! adjusted experimentally
#endif /* !DEBUG */
         case (2)
#ifdef DEBUG
            ! the maximum convergence is at aux_R(1) = 0.35 +/- 0.1 and aux_R(3) = 0.20 +- 0.05
            opfn_c_ff(:) = 0.5 * [ 1+aux_R(1), -2.*aux_R(1)+aux_R(3), aux_R(1)-aux_R(3) ]
            ! the maximum convergence is at aux_R(2) = -0.01 +/- 0.02 and aux_R(4) = -0.125 +- 0.025
            opfn_c_cf(:) = [ 1.0 + 2*aux_R(2), -2.*aux_R(2)+aux_R(4), aux_R(2)-aux_R(4) ]
#else /* !DEBUG */
            opfn_c_ff(:) = [ 27., -10., 3. ] / 40. ! adjusted experimentally
            opfn_c_cf(:) = [  8.,  -1., 1. ] / 8.
#endif /* !DEBUG */
         case default
            call die("mg:mmpg opfn_c_[cf]f(:)")
      end select

      curl => base
      do while (associated(curl))

         ijks(:) = curl%ijkse(:, LO) - curl%off(:)  ! add this to convert absolute cell coordinates to local indices. (+mg_nb - off(:))
         per(:) = 0
         where (curl%dom%periodic(:)) per(:) = curl%dom%n_d(:)

         do d = xdim, zdim
            dmask(:) = .false.
            dmask(d) = .true.
            do lh = LO, HI
               if (has_dir(d) .and. .not. curl%empty) then
                  hl = LO+HI-lh

                  ! find coarse target for receiving data to be prolonged
                  if (associated(curl%coarser) .and. .not. is_external(d, lh)) then
                     procmask(:) = 0
                     ! two layers of cells are required for even locations (+ two layers per each interpolation order)
                     ! one layer of cells is required for odd locations (the local domain face is exactly at the centers of coarse cells)
                     coarsened(:, :) = curl%dom%se(proc, :, :)/2
                     coarsened(d, hl) = coarsened(d, lh)
                     select case (lh)
                        case (LO)
                           if (mod(curl%off(d),               2_LONG) == 0) coarsened(d, :) = coarsened(d, :) + [ -1-ord_prolong_face_norm,   ord_prolong_face_norm ]
                        case (HI)
                           if (mod(curl%off(d) + curl%n_b(d), 2_LONG) == 0) coarsened(d, :) = coarsened(d, :) + [   -ord_prolong_face_norm, 1+ord_prolong_face_norm ]
                     end select

                     do j = 0, nproc-1
                        call is_neigh(coarsened(:,:), curl%coarser%dom%se(j, :, :), neigh(:,:), sharing, corner, face, per)
                        if (sharing) procmask(j) = 1
                     enddo
                     allocate(curl%pfc_tgt(d, lh)%seg(count(procmask(:) /= 0)))

                     g = 0
                     do j = 0, nproc-1
                        if (procmask(j) /= 0) then
                           g = g + 1
                           if (.not. allocated(curl%pfc_tgt(d, lh)%seg) .or. g>ubound(curl%pfc_tgt(d, lh)%seg, dim=1)) call die("mg:mmpg pfc_tgt g>")
                           seg => curl%pfc_tgt(d, lh)%seg(g)
                           if (allocated(seg%buf)) then
                              call warn("mg:mmpg o seg%buf a a")
                              deallocate(seg%buf)
                           endif
                           seg%proc = j
                           ! find cross-section of own face segment with refined coarse segment
                           b_layer(:, :) = curl%dom%se(proc, :, :)
                           b_layer(d, hl) = b_layer(d, lh)
                           !b_layer(d, lh) = b_layer(d, lh) + 2*lh-LO-HI ! extend to two layers of buffer

                           where (.not. dmask(:)) ! find extents perpendicular to d
                              seg%se(:, LO) = max(b_layer(:, LO), curl%coarser%dom%se(j, :, LO)*2  )
                              seg%se(:, HI) = min(b_layer(:, HI), curl%coarser%dom%se(j, :, HI)*2+1)
                           endwhere
                           seg%se(d, :) = b_layer(d, :)
                           !if (j /= proc)
                           allocate(seg%buf(seg%se(xdim, HI)/2-seg%se(xdim, LO)/2 + 1, &
                                &           seg%se(ydim, HI)/2-seg%se(ydim, LO)/2 + 1, &
                                &           seg%se(zdim, HI)/2-seg%se(zdim, LO)/2 + 1))
                           ! not counted in mb_alloc

                           seg%se(:, LO) = seg%se(:, LO) - curl%off(:) !+ ijks(:)
                           seg%se(:, HI) = seg%se(:, HI) - curl%off(:) !+ ijks(:)
                        endif
                     enddo

                  endif

                  ! find fine target(s) for sending the data to be prolonged
                  if (associated(curl%finer)) then
                     procmask(:) = 0
                     do j = 0, nproc-1
                        is_internal_fine = curl%finer%dom%periodic(d)
                        coarsened(:, :) = curl%finer%dom%se(j, :, :)/2
                        coarsened(d, hl) = coarsened(d, lh)
                        select case (lh)
                           case (LO)
                              if (mod(curl%finer%dom%se(j, d, LO),     2_LONG) == 0) coarsened(d, :) = coarsened(d, :) + [ -1-ord_prolong_face_norm,   ord_prolong_face_norm ]
                              is_internal_fine = is_internal_fine .or. (curl%finer%dom%se(j, d, lh) /= 0)
                           case (HI)
                              if (mod(curl%finer%dom%se(j, d, HI) + 1, 2_LONG) == 0) coarsened(d, :) = coarsened(d, :) + [   -ord_prolong_face_norm, 1+ord_prolong_face_norm ]
                              is_internal_fine = is_internal_fine .or. (curl%finer%dom%se(j, d, lh) + 1 < curl%finer%dom%n_d(d))
                        end select
                        if (is_internal_fine) then
                           call is_neigh(coarsened(:, :), curl%dom%se(proc, :, :), neigh(:,:), sharing, corner, face, per)
                           if (sharing) procmask(j) = 1
                        endif
                     enddo
                     allocate(curl%pff_tgt(d, lh)%seg(count(procmask(:) /= 0)))

                     g = 0
                     do j = 0, nproc-1
                        if (procmask(j) /= 0) then
                           g = g + 1
                           if (.not. allocated(curl%pff_tgt(d, lh)%seg) .or. g>ubound(curl%pff_tgt(d, lh)%seg, dim=1)) call die("mg:mmpg pff_tgt g>")
                           seg => curl%pff_tgt(d, lh)%seg(g)
                           if (allocated(seg%buf)) then
                              call warn("mg:mmpg o seg%buf a a")
                              deallocate(seg%buf)
                           endif
                           seg%proc = j

                           ! find cross-section of own segment with coarsened fine face segment
                           coarsened(:, :) = curl%finer%dom%se(j, :, :)
                           coarsened(d, hl) = coarsened(d, lh)
                           coarsened(:, :) = coarsened(:, :)/2
                           where (.not. dmask(:))
                              seg%se(:, LO) = max(curl%dom%se(proc, :, LO), coarsened(:, LO))
                              seg%se(:, HI) = min(curl%dom%se(proc, :, HI), coarsened(:, HI))
                           endwhere
                           seg%se(d, :) = coarsened(d, :)
                           !if (j /= proc)
                           allocate(seg%buf(seg%se(xdim, HI)-seg%se(xdim, LO) + 1, &
                                &           seg%se(ydim, HI)-seg%se(ydim, LO) + 1, &
                                &           seg%se(zdim, HI)-seg%se(zdim, LO) + 1))

                           coarsened(:, :) = curl%finer%dom%se(j, :, :)
                           coarsened(d, hl) = coarsened(d, lh)
                           coarsened(d, lh) = coarsened(d, lh) + 2*lh-LO-HI ! extend to two layers of buffer
                           coarsened(:, :) = coarsened(:, :)/2
                           coarsened(d, :) = coarsened(d, :) + [ -ord_prolong_face_norm, ord_prolong_face_norm ]

                           seg%se(:, LO) = max(curl%dom%se(proc, :, LO), coarsened(:, LO)) + ijks(:)
                           seg%se(:, HI) = min(curl%dom%se(proc, :, HI), coarsened(:, HI)) + ijks(:)

                           coarsened(d, :) = coarsened(d, :) - [ -ord_prolong_face_norm, ord_prolong_face_norm ] ! revert broadening
                           allocate(seg%f_lay(seg%se(d, HI) - seg%se(d, LO) + 1))
                           do l = 1, size(seg%f_lay(:))
                              seg%f_lay(l)%layer = l + int(seg%se(d, LO), kind=4) - 1
                              nl = int(minval(abs(seg%f_lay(l)%layer - ijks(d) - coarsened(d, :))), kind=4)
                              if (mod(curl%finer%dom%se(j, d, lh) + lh - LO, 2_LONG) == 0) then ! fine face at coarse face
                                 seg%f_lay(l)%coeff = opfn_c_ff(nl)
                              else                                                              ! fine face at coarse center
                                 seg%f_lay(l)%coeff = opfn_c_cf(nl)
                              endif
                           enddo

                        endif
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
      use multipole,     only: cleanup_multipole
      use multigridvars, only: lvl, plvl, base, tgt_list

      implicit none

      integer :: i, g, ib
      integer, parameter :: nseg = 2*(HI-LO+1)*ndims
      type(tgt_list), dimension(nseg) :: io_tgt
      type(plvl), pointer :: curl

      call cleanup_multipole

      do i = 1, nold
         if (allocated(inner%old(i)%soln)) deallocate(inner%old(i)%soln)
         if (allocated(outer%old(i)%soln)) deallocate(outer%old(i)%soln)
      enddo

      if (allocated(vstat%factor)) deallocate(vstat%factor)
      if (allocated(vstat%time)) deallocate(vstat%time)

      if (allocated(lvl)) then
         curl => base
         do while (associated(curl))
            if (allocated(curl%fft))     deallocate(curl%fft)
            if (allocated(curl%fftr))    deallocate(curl%fftr)
            if (allocated(curl%src))     deallocate(curl%src)
            if (allocated(curl%Green3D)) deallocate(curl%Green3D)

            if (curl%planf /= 0) call dfftw_destroy_plan(curl%planf)
            if (curl%plani /= 0) call dfftw_destroy_plan(curl%plani)

            io_tgt(1:nseg) = [ curl%pfc_tgt, curl%pff_tgt ]
            do ib = 1, nseg
               if (allocated(io_tgt(ib)%seg)) then
                  do g = 1, ubound(io_tgt(ib)%seg, dim=1)
                     if (allocated(io_tgt(ib)%seg(g)%buf)) deallocate(io_tgt(ib)%seg(g)%buf)
                     if (allocated(io_tgt(ib)%seg(g)%f_lay)) deallocate(io_tgt(ib)%seg(g)%f_lay)
                  enddo
                  deallocate(io_tgt(ib)%seg)
               endif
            enddo

            curl => curl%finer
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
      use multigridvars,    only: lvl, plvl, base, roof, stdout, solution

      implicit none

      type(soln_history), intent(inout) :: history !< inner or outer potential history for recycling

      integer :: p0, p1, p2, ordt
      real, dimension(3)  :: dt_fac
      type(plvl), pointer :: curl

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
            curl => base
            do while (associated(curl))
               curl%mgvar(:, :, :, solution) = 0.
               curl => curl%finer
            enddo
            history%old(:)%time = -huge(1.0)
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

      call check_dirty(roof%level, solution, "init_soln")

   end subroutine init_solution

!!$ ============================================================================
!>
!! \brief Make a local copy of source (density) and multiply by 4 pi G
!<

   subroutine init_source(dens)

#ifdef JEANS_PROBLEM
      use problem_pub,        only: jeans_d0, jeans_mode ! hack for tests
#endif /* JEANS_PROBLEM */
      use units,              only: fpiG
      use grid,               only: cg
      use dataio_pub,         only: die
      use multigridhelpers,   only: set_dirty, check_dirty
      use multigridbasefuncs, only: substract_average
      use multigridvars,      only: roof, source, is_external, bnd_periodic, bnd_dirichlet, bnd_givenval
      use mpisetup,           only: geometry_type
      use constants,          only: GEO_RPZ, LO, HI, xdim, ydim, zdim

      implicit none

      real, optional, dimension(:,:,:), intent(in)  :: dens !< input source field or nothing for empty space
      real :: fac
      integer :: i

      call set_dirty(source)

      if (present(dens)) then
         roof%mgvar(roof%is:roof%ie, roof%js:roof%je, roof%ks:roof%ke, source) = fpiG * dens(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke)
      else
         if (grav_bnd /= bnd_givenval) call die("[multigrid_gravity:init_source] empty space allowed only for given value boundaries.")
         roof%mgvar(roof%is:roof%ie, roof%js:roof%je, roof%ks:roof%ke, source) = 0.
      endif

      select case (grav_bnd)
         case (bnd_periodic) ! probably also bnd_neumann
            call substract_average(roof%level, source)
         case (bnd_dirichlet)
#ifdef JEANS_PROBLEM
            if (jeans_mode == 1) roof%mgvar(roof%is:roof%ie, roof%js:roof%je, roof%ks:roof%ke, source) = &
                 &         roof%mgvar(roof%is:roof%ie, roof%js:roof%je, roof%ks:roof%ke, source) - fpiG * jeans_d0 ! remove density bias
#endif /* JEANS_PROBLEM */
         case (bnd_givenval) ! convert potential into a layer of imaginary mass (subtract second derivative normal to computational domain boundary)
            if (is_external(xdim, LO)) then
               fac = 2. * roof%idx2 / fpiG
               if (geometry_type == GEO_RPZ .and. roof%x(roof%is) /= 0.) fac = fac - 1./(roof%dx * roof%x(roof%is) * fpiG) !> BEWARE is it roof%x(ie), roof%x(ie+1) or something in the middle?
               roof%mgvar       (roof%is, roof%js:roof%je, roof%ks:roof%ke, source) = &
                    & roof%mgvar(roof%is, roof%js:roof%je, roof%ks:roof%ke, source) - &
                    & roof%bnd_x(         roof%js:roof%je, roof%ks:roof%ke, LO) * fac
            endif
            if (is_external(xdim, HI)) then
               fac = 2. * roof%idx2 / fpiG
               if (geometry_type == GEO_RPZ .and. roof%x(roof%ie) /= 0.) fac = fac - 1./(roof%dx * roof%x(roof%ie) * fpiG) !> BEWARE is it roof%x(ie), roof%x(ie+1) or something in the middle?
               roof%mgvar       (roof%ie, roof%js:roof%je, roof%ks:roof%ke, source) = &
                    & roof%mgvar(roof%ie, roof%js:roof%je, roof%ks:roof%ke, source) - &
                    & roof%bnd_x(         roof%js:roof%je, roof%ks:roof%ke, HI) * fac
            endif
            if (is_external(ydim, LO)) then
               if (geometry_type == GEO_RPZ) then
                  do i = roof%is, roof%ie
                     if (roof%x(i) /= 0.) then
                        roof%mgvar       (i, roof%js, roof%ks:roof%ke, source) = &
                             & roof%mgvar(i, roof%js, roof%ks:roof%ke, source) - &
                             & roof%bnd_y(i,          roof%ks:roof%ke, LO) * 2. * roof%idy2 / fpiG / roof%x(i)**2
                     endif
                  enddo
               else
                  roof%mgvar       (roof%is:roof%ie, roof%js, roof%ks:roof%ke, source) = &
                       & roof%mgvar(roof%is:roof%ie, roof%js, roof%ks:roof%ke, source) - &
                       & roof%bnd_y(roof%is:roof%ie,          roof%ks:roof%ke, LO) * 2. * roof%idy2 / fpiG
               endif
            endif
            if (is_external(ydim, HI)) then
               if (geometry_type == GEO_RPZ) then
                  do i = roof%is, roof%ie
                     if (roof%x(i) /= 0.) then
                        roof%mgvar       (i, roof%je, roof%ks:roof%ke, source) = &
                             & roof%mgvar(i, roof%je, roof%ks:roof%ke, source) - &
                             & roof%bnd_y(i,          roof%ks:roof%ke, HI) * 2. * roof%idy2 / fpiG / roof%x(i)**2
                     endif
                  enddo
               else
                  roof%mgvar       (roof%is:roof%ie, roof%je, roof%ks:roof%ke, source) = &
                       & roof%mgvar(roof%is:roof%ie, roof%je, roof%ks:roof%ke, source) - &
                       & roof%bnd_y(roof%is:roof%ie,          roof%ks:roof%ke, HI) * 2. * roof%idy2 / fpiG
               endif
            endif
            if (is_external(zdim, LO)) roof%mgvar(roof%is:roof%ie, roof%js:roof%je, roof%ks, source) = &
                 &                     roof%mgvar(roof%is:roof%ie, roof%js:roof%je, roof%ks, source) - &
                 &                     roof%bnd_z(roof%is:roof%ie, roof%js:roof%je,          LO) * 2. * roof%idz2 / fpiG
            if (is_external(zdim, HI)) roof%mgvar(roof%is:roof%ie, roof%js:roof%je, roof%ke, source) = &
                 &                     roof%mgvar(roof%is:roof%ie, roof%js:roof%je, roof%ke, source) - &
                 &                     roof%bnd_z(roof%is:roof%ie, roof%js:roof%je,          HI) * 2. * roof%idz2 / fpiG
            !> \todo compactify the above mess
         case default
            call die("[multigrid_gravity:init_source] Unknown boundary type")
      end select

      call check_dirty(roof%level, source, "init_src")

   end subroutine init_source

!!$ ============================================================================
!>
!! \brief This routine manages old copies of potential for recycling.
!<

   subroutine store_solution(history)

      use mpisetup,          only: t
      use multigridmpifuncs, only: mpi_multigrid_bnd
      use multigridvars,     only: roof, bnd_isolated, bnd_givenval, solution, mg_nb, extbnd_extrapolate, extbnd_mirror

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
         call mpi_multigrid_bnd(roof%level, solution, mg_nb, extbnd_extrapolate)
      else
         call mpi_multigrid_bnd(roof%level, solution, mg_nb, extbnd_mirror)
      endif

   end subroutine store_solution

!!$ ============================================================================
!>
!! \brief Multigrid gravity driver. This is the only multigrid routine intended to be called from the gravity module.
!! This routine is also responsible for communicating the solution to the rest of world via sgp array.
!<

   subroutine multigrid_solve_grav(dens)

      use constants,     only: xdim, ydim, zdim
      use timer,         only: set_timer
      use grid,          only: cg
      use mpisetup,      only: has_dir
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
#else /* !COSM_RAYS */
         vstat%cprefix = ""
#endif /* !COSM_RAYS */
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

      cg%sgp%arr(isb:ieb, jsb:jeb, ksb:keb) = roof%mgvar(:, :, :, solution)

      if (isolated) then
         grav_bnd = bnd_givenval

         vstat%cprefix = "Go-"
         call multipole_solver
         call init_source

         call vcycle_hg(outer)

         cg%sgp%arr(isb:ieb, jsb:jeb, ksb:keb) = cg%sgp%arr(isb:ieb, jsb:jeb, ksb:keb) + roof%mgvar(:, :, :, solution)

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

      use constants,          only: cbuff_len
      use mpisetup,           only: master, nproc
      use timer,              only: set_timer
      use multigridhelpers,   only: set_dirty, check_dirty, mg_write_log, brief_v_log, do_ascii_dump, numbered_ascii_dump
      use multigridbasefuncs, only: norm_sq, restrict_all, substract_average
      use dataio_pub,         only: msg, die, warn
      use multigridvars,      only: plvl, roof, base, source, solution, correction, defect, verbose_vcycle, bnd_periodic, stdout, tot_ts, ts

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
      type(plvl), pointer :: curl

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

!      if (.not. history%valid .and. prefer_rbgs_relaxation) call approximate_solution(roof%level, source, solution) ! not necessary when init_solution called FFT
! difficult statement: for approximate_solution_fft it requires to pass a flag to use guardcells instead of prolonging faces.
! how much does it improve? (make a benchmark at some point)

      ! iterations
      do v = 0, max_cycles

         call set_dirty(defect)
         call residual(roof%level, source, solution, defect)
         if (grav_bnd == bnd_periodic) call substract_average(roof%level, defect)
         call check_dirty(roof%level, defect, "residual")

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

         curl => base
         do while (associated(curl))
            call approximate_solution(curl%level, defect, correction)
            call check_dirty(curl%level, correction, "Vup relax+")
            curl => curl%finer
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

      call check_dirty(roof%level, solution, "final_solution")

      call store_solution(history)

   end subroutine vcycle_hg

!!$ ============================================================================
!>
!! \brief Calculate the residuum for the Poisson equation.
!<

   subroutine residual(lev, src, soln, def)

      use dataio_pub,            only: die
      use multigridvars,         only: ngridvars, base, roof

      implicit none

      integer, intent(in) :: lev  !< level for which approximate the solution
      integer, intent(in) :: src  !< index of source in lvl()%mgvar
      integer, intent(in) :: soln !< index of solution in lvl()%mgvar
      integer, intent(in) :: def  !< index of defect in lvl()%mgvar

      if (any( [ src, soln, def ] <= 0) .or. any( [ src, soln, def ] > ngridvars)) call die("[multigrid_gravity:residual] Invalid variable index")
      if (lev < base%level .or. lev > roof%level) call die("[multigrid_gravity:residual] Invalid level number")

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

      use constants,          only: xdim, ydim, zdim, ndims, GEO_XYZ, GEO_RPZ
      use dataio_pub,         only: die
      use mpisetup,           only: has_dir, eff_dim, geometry_type
      use multigridhelpers,   only: multidim_code_3D
      use multigridmpifuncs,  only: mpi_multigrid_bnd
      use multigridvars,      only: lvl, plvl, extbnd_antimirror

      implicit none

      integer, intent(in) :: lev  !< level for which approximate the solution
      integer, intent(in) :: src  !< index of source in lvl()%mgvar
      integer, intent(in) :: soln !< index of solution in lvl()%mgvar
      integer, intent(in) :: def  !< index of defect in lvl()%mgvar

      real    :: L0, Lx, Ly, Lz, Lx1
      integer :: i, j, k
      type(plvl), pointer :: curl

      curl => lvl(lev)

      call mpi_multigrid_bnd(lev, soln, 1, extbnd_antimirror) ! no corners required

      ! Coefficients for a simplest 3-point Laplacian operator: [ 1, -2, 1 ]
      ! for 2D and 1D setups appropriate elements of [ Lx, Ly, Lz ] should be == 0.
      Lx = curl%idx2
      Ly = curl%idy2
      Lz = curl%idz2
      L0 = -2. * (Lx + Ly + Lz)

      ! Possible optimization candidate: reduce cache misses (secondary importance, cache-aware implementation required)
      ! Explicit loop over k gives here better performance than array operation due to less cache misses (at least on 32^3 and 64^3 arrays)

      select case (geometry_type)
         case (GEO_XYZ)
            if (eff_dim == ndims .and. .not. multidim_code_3D) then
               do k = curl%ks, curl%ke
                  curl%mgvar       (curl%is  :curl%ie,   curl%js  :curl%je,   k,   def)        = &
                       & curl%mgvar(curl%is  :curl%ie,   curl%js  :curl%je,   k,   src)        - &
                       ( curl%mgvar(curl%is-1:curl%ie-1, curl%js  :curl%je,   k,   soln)       + &
                       & curl%mgvar(curl%is+1:curl%ie+1, curl%js  :curl%je,   k,   soln)) * Lx - &
                       ( curl%mgvar(curl%is  :curl%ie,   curl%js-1:curl%je-1, k,   soln)       + &
                       & curl%mgvar(curl%is  :curl%ie,   curl%js+1:curl%je+1, k,   soln)) * Ly - &
                       ( curl%mgvar(curl%is  :curl%ie,   curl%js  :curl%je,   k-1, soln)       + &
                       & curl%mgvar(curl%is  :curl%ie,   curl%js  :curl%je,   k+1, soln)) * Lz - &
                       & curl%mgvar(curl%is  :curl%ie,   curl%js  :curl%je,   k,   soln)  * L0
               enddo
            else
               ! In 3D this implementation can give a bit more cache misses, few times more writes and significantly more instructions executed than monolithic 3D above
               do k = curl%ks, curl%ke
                  curl%mgvar       (curl%is  :curl%ie,   curl%js  :curl%je,   k,   def)   = &
                       & curl%mgvar(curl%is  :curl%ie,   curl%js  :curl%je,   k,   src)   - &
                       & curl%mgvar(curl%is  :curl%ie,   curl%js  :curl%je,   k,   soln)  * L0
                  if (has_dir(xdim)) &
                       & curl%mgvar(curl%is  :curl%ie,   curl%js  :curl%je,   k,   def)   = &
                       & curl%mgvar(curl%is  :curl%ie,   curl%js  :curl%je,   k,   def)   - &
                       ( curl%mgvar(curl%is-1:curl%ie-1, curl%js  :curl%je,   k,   soln)  + &
                       & curl%mgvar(curl%is+1:curl%ie+1, curl%js  :curl%je,   k,   soln)) * Lx
                  if (has_dir(ydim)) &
                       & curl%mgvar(curl%is  :curl%ie,   curl%js  :curl%je,   k,   def)   = &
                       & curl%mgvar(curl%is  :curl%ie,   curl%js  :curl%je,   k,   def)   - &
                       ( curl%mgvar(curl%is  :curl%ie,   curl%js-1:curl%je-1, k,   soln)  + &
                       & curl%mgvar(curl%is  :curl%ie,   curl%js+1:curl%je+1, k,   soln)) * Ly
                  if (has_dir(zdim)) &
                       & curl%mgvar(curl%is  :curl%ie,   curl%js  :curl%je,   k,   def)   = &
                       & curl%mgvar(curl%is  :curl%ie,   curl%js  :curl%je,   k,   def)   - &
                       ( curl%mgvar(curl%is  :curl%ie,   curl%js  :curl%je,   k-1, soln)  + &
                       & curl%mgvar(curl%is  :curl%ie,   curl%js  :curl%je,   k+1, soln)) * Lz
               enddo
            endif
         case (GEO_RPZ)
            Lx = curl%idx2
            Lz = curl%idz2
            do k = curl%ks, curl%ke
               do j = curl%js, curl%je
                  do i = curl%is, curl%ie
                     if (curl%x(i) /= 0.) then !> \todo convert Ly, Lx1 and L0 into precomputed arrays
                        Ly = curl%idy2 / curl%x(i)**2 ! cylindrical factor
                        Lx1 = 0.5 / (curl%dx * curl%x(i))
                     else
                        Ly = 0.
                        Lx1 = 0.
                     endif
                     L0 = -2. * (Lx + Ly + Lz)
                     curl%mgvar(i, j, k, def) = curl%mgvar(i, j, k, src) - curl%mgvar(i, j, k, soln) * L0
                     if (has_dir(xdim)) curl%mgvar(i,   j,   k,   def)  = curl%mgvar(i,   j,   k,   def)        - &
                          &            (curl%mgvar(i+1, j,   k,   soln) + curl%mgvar(i-1, j,   k,   soln)) * Lx - &
                          &            (curl%mgvar(i+1, j,   k,   soln) - curl%mgvar(i-1, j,   k,   soln)) * Lx1    ! cylindrical term
                     if (has_dir(ydim)) curl%mgvar(i,   j,   k,   def)  = curl%mgvar(i,   j,   k,   def)        - &
                          &            (curl%mgvar(i,   j+1, k,   soln) + curl%mgvar(i,   j-1, k,   soln)) * Ly
                     if (has_dir(zdim)) curl%mgvar(i,   j,   k,   def)  = curl%mgvar(i,   j,   k,   def)        - &
                          &            (curl%mgvar(i,   j,   k+1, soln) + curl%mgvar(i,   j,   k-1, soln)) * Lz
                  enddo
               enddo
            enddo
         case default
            call die("[multigrid_gravity:residual2] Unsupported geometry.")
      end select

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
      use multigridvars,      only: lvl, plvl, bnd_givenval, extbnd_antimirror
      use constants,          only: ndims

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
      type(plvl), pointer :: curl

      curl => lvl(lev)

      if (eff_dim<ndims) call die("[multigrid_gravity:residual4] Only 3D is implemented")

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

      Lx1 = c41 * curl%idx2
      Ly1 = c41 * curl%idy2
      Lz1 = c41 * curl%idz2
      Lx2 = c42 * curl%idx2
      Ly2 = c42 * curl%idy2
      Lz2 = c42 * curl%idz2
!      L0  = c40 * (curl%idx2 + curl%idy2 + curl%idz2 )
      L0 = -2. * (Lx1 + Lx2 + Ly1 + Ly2 + Lz1 + Lz2)

      !> \deprecated BEWARE: cylindrical factors go here
      curl%mgvar     (curl%is  :curl%ie,   curl%js  :curl%je,   curl%ks  :curl%ke,   def)        = &
           curl%mgvar(curl%is  :curl%ie,   curl%js  :curl%je,   curl%ks  :curl%ke,   src)        - &
           curl%mgvar(curl%is-2:curl%ie-2, curl%js  :curl%je,   curl%ks  :curl%ke,   soln) * Lx2 - &
           curl%mgvar(curl%is+2:curl%ie+2, curl%js  :curl%je,   curl%ks  :curl%ke,   soln) * Lx2 - &
           curl%mgvar(curl%is-1:curl%ie-1, curl%js  :curl%je,   curl%ks  :curl%ke,   soln) * Lx1 - &
           curl%mgvar(curl%is+1:curl%ie+1, curl%js  :curl%je,   curl%ks  :curl%ke,   soln) * Lx1 - &
           curl%mgvar(curl%is  :curl%ie,   curl%js-2:curl%je-2, curl%ks  :curl%ke,   soln) * Ly2 - &
           curl%mgvar(curl%is  :curl%ie,   curl%js+2:curl%je+2, curl%ks  :curl%ke,   soln) * Ly2 - &
           curl%mgvar(curl%is  :curl%ie,   curl%js-1:curl%je-1, curl%ks  :curl%ke,   soln) * Ly1 - &
           curl%mgvar(curl%is  :curl%ie,   curl%js+1:curl%je+1, curl%ks  :curl%ke,   soln) * Ly1 - &
           curl%mgvar(curl%is  :curl%ie,   curl%js  :curl%je,   curl%ks-2:curl%ke-2, soln) * Lz2 - &
           curl%mgvar(curl%is  :curl%ie,   curl%js  :curl%je,   curl%ks+2:curl%ke+2, soln) * Lz2 - &
           curl%mgvar(curl%is  :curl%ie,   curl%js  :curl%je,   curl%ks-1:curl%ke-1, soln) * Lz1 - &
           curl%mgvar(curl%is  :curl%ie,   curl%js  :curl%je,   curl%ks+1:curl%ke+1, soln) * Lz1 - &
           curl%mgvar(curl%is  :curl%ie,   curl%js  :curl%je,   curl%ks  :curl%ke,   soln) * L0

      ! WARNING: not optimized
      if (grav_bnd == bnd_givenval) then ! probably also in some other cases
         ! Use L2 Laplacian in two layers of cells next to the boundary because L4 seems to be incompatible with present image mass construction
         Lx = c21 * curl%idx2
         Ly = c21 * curl%idy2
         Lz = c21 * curl%idz2
         L0 = -2. * (Lx + Ly + Lz)

         do k = curl%ks, curl%ke
            do j = curl%js, curl%je
               do i = curl%is, curl%ie
                  if ( i<curl%is+L2w .or. i>curl%ie-L2w .or. j<curl%js+L2w .or. j>curl%je-L2w .or. k<curl%ks+L2w .or. k>curl%ke-L2w) then
                     curl%mgvar       (i,   j,   k,   def)   = curl%mgvar(i,   j,   k,   src)        - &
                          ( curl%mgvar(i-1, j,   k,   soln)  + curl%mgvar(i+1, j,   k,   soln)) * Lx - &
                          ( curl%mgvar(i,   j-1, k,   soln)  + curl%mgvar(i,   j+1, k,   soln)) * Ly - &
                          ( curl%mgvar(i,   j,   k-1, soln)  + curl%mgvar(i,   j,   k+1, soln)) * Lz - &
                          & curl%mgvar(i,   j,   k,   soln)  * L0
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
      use multigridvars,      only: lvl, plvl, base, roof, ngridvars, correction

      implicit none

      integer, intent(in) :: lev  !< level for which approximate the solution
      integer, intent(in) :: src  !< index of source in lvl()%mgvar
      integer, intent(in) :: soln !< index of solution in lvl()%mgvar
      type(plvl), pointer :: curl

      if (any( [ src, soln ] <= 0) .or. any( [ src, soln ] > ngridvars)) call die("[multigrid_gravity:approximate_solution] Invalid variable index.")
      if (lev < base%level .or. lev > roof%level) call die("[multigrid_gravity:approximate_solution] Invalid level number.")

      curl => lvl(lev)

      call check_dirty(lev, src, "approx_soln src-")

      if (curl%fft_type /= fft_none) then
         call approximate_solution_fft(lev, src, soln)
      else
         call check_dirty(lev, soln, "approx_soln soln-")
         call approximate_solution_rbgs(lev, src, soln)
      endif

      if (prefer_rbgs_relaxation .and. soln == correction .and. .not. associated(curl, roof)) call prolong_level(lev, correction)
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

      use dataio_pub,         only: die
      use constants,          only: xdim, ydim, zdim, ndims, GEO_XYZ, GEO_RPZ
      use mpisetup,           only: has_dir, eff_dim, geometry_type
      use multigridhelpers,   only: dirty_debug, check_dirty, multidim_code_3D, dirty_label
      use multigridmpifuncs,  only: mpi_multigrid_bnd
      use multigridvars,      only: lvl, plvl, base, extbnd_antimirror

      implicit none

      integer, intent(in) :: lev  !< level for which approximate the solution
      integer, intent(in) :: src  !< index of source in lvl()%mgvar
      integer, intent(in) :: soln !< index of solution in lvl()%mgvar

      integer, parameter :: RED_BLACK = 2 !< the checkerboard requires two sweeps

      integer :: n, i, j, k, i1, j1, k1, id, jd, kd
      integer :: nsmoo
      real    :: crx, crx1, cry, crz, cr
      type(plvl), pointer :: curl

      curl => lvl(lev)

      if (curl%empty) return

      if (associated(curl, base)) then
         nsmoo = nsmoob
      else
         nsmoo = nsmool
      endif

      if (geometry_type == GEO_RPZ .and. .not. multidim_code_3D) call die("[multigrid_gravity:approximate_solution_rbgs] multidim_code_3D = .false. not implemented")

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

         if (eff_dim==ndims .and. .not. multidim_code_3D) then
            do k = curl%ks, curl%ke
               do j = curl%js, curl%je
                  i1 = curl%is + mod(n+j+k, RED_BLACK)
                  if (geometry_type == GEO_RPZ) then
!!$                     curl%mgvar(i1  :curl%ie  :2, j,   k,   soln) = &
!!$                          curl%rx * (curl%mgvar(i1-1:curl%ie-1:2, j,   k,   soln) + curl%mgvar(i1+1:curl%ie+1:2, j,   k,   soln)) + &
!!$                          curl%ry * (curl%mgvar(i1  :curl%ie  :2, j-1, k,   soln) + curl%mgvar(i1:  curl%ie:  2, j+1, k,   soln)) + &
!!$                          curl%rz * (curl%mgvar(i1  :curl%ie  :2, j,   k-1, soln) + curl%mgvar(i1:  curl%ie:  2, j,   k+1, soln)) - &
!!$                          curl%r  *  curl%mgvar(i1  :curl%ie  :2, j,   k,   src)  + &
!!$                          curl%rx * (curl%mgvar(i1+1:curl%ie+1:2, j,   k,   soln) - curl%mgvar(i1-1:curl%ie-1:2, j,   k,   soln)) * fac(i1:curl%ie:2)
                     call die("[multigrid_gravity:approximate_solution_rbgs] This variant of relaxation loop is not implemented for cylindrical coordinates.")
                  else
                     curl%mgvar(i1  :curl%ie  :2, j,   k,   soln) = &
                          curl%rx * (curl%mgvar(i1-1:curl%ie-1:2, j,   k,   soln) + curl%mgvar(i1+1:curl%ie+1:2, j,   k,   soln)) + &
                          curl%ry * (curl%mgvar(i1  :curl%ie  :2, j-1, k,   soln) + curl%mgvar(i1:  curl%ie:  2, j+1, k,   soln)) + &
                          curl%rz * (curl%mgvar(i1  :curl%ie  :2, j,   k-1, soln) + curl%mgvar(i1:  curl%ie:  2, j,   k+1, soln)) - &
                          curl%r  *  curl%mgvar(i1  :curl%ie  :2, j,   k,   src)
                  endif
               enddo
            enddo
         else
            ! In 3D this variant significantly increases instruction count and also some data read
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

            if (kd == RED_BLACK) k1 = curl%ks + mod(n, RED_BLACK)
            select case (geometry_type)
               case (GEO_XYZ)
                  do k = k1, curl%ke, kd
                     if (jd == RED_BLACK) j1 = curl%js + mod(n+k, RED_BLACK)
                     do j = j1, curl%je, jd
                        if (id == RED_BLACK) i1 = curl%is + mod(n+j+k, RED_BLACK)
                        curl%mgvar                           (i1:  curl%ie  :id, j,   k,   soln) = &
                             & (1. - Jacobi_damp)* curl%mgvar(i1  :curl%ie  :id, j,   k,   soln) - &
                             &       Jacobi_damp * curl%mgvar(i1  :curl%ie  :id, j,   k,   src)  * curl%r
                        if (has_dir(xdim))         curl%mgvar(i1  :curl%ie  :id, j,   k,   soln) = curl%mgvar(i1:  curl%ie:  id, j,   k,   soln)  + &
                             &       Jacobi_damp *(curl%mgvar(i1-1:curl%ie-1:id, j,   k,   soln) + curl%mgvar(i1+1:curl%ie+1:id, j,   k,   soln)) * curl%rx
                        if (has_dir(ydim))         curl%mgvar(i1  :curl%ie  :id, j,   k,   soln) = curl%mgvar(i1:  curl%ie:  id, j,   k,   soln)  + &
                             &       Jacobi_damp *(curl%mgvar(i1  :curl%ie  :id, j-1, k,   soln) + curl%mgvar(i1:  curl%ie:  id, j+1, k,   soln)) * curl%ry
                        if (has_dir(zdim))         curl%mgvar(i1  :curl%ie  :id, j,   k,   soln) = curl%mgvar(i1:  curl%ie:  id, j,   k,   soln)  + &
                             &       Jacobi_damp *(curl%mgvar(i1  :curl%ie  :id, j,   k-1, soln) + curl%mgvar(i1:  curl%ie:  id, j,   k+1, soln)) * curl%rz
                     enddo
                  enddo
               case (GEO_RPZ)
                  do k = k1, curl%ke, kd
                     if (jd == RED_BLACK) j1 = curl%js + mod(n+k, RED_BLACK)
                     do j = j1, curl%je, jd
                        if (id == RED_BLACK) i1 = curl%is + mod(n+j+k, RED_BLACK)
                        do i = i1, curl%ie, id
                           cr  = overrelax / 2.
                           crx = curl%dvol2 * curl%idx2 * curl%x(i)**2
                           cry = curl%dvol2 * curl%idy2
                           crz = curl%dvol2 * curl%idz2 * curl%x(i)**2
                           cr  = cr  / (crx + cry + crz)
                           crx = overrelax_x * crx * cr
                           cry = overrelax_y * cry * cr
                           crz = overrelax_z * crz * cr
                           cr  = cr * curl%dvol2 * curl%x(i)**2

                           crx1 = 2. * curl%x(i) / curl%dx
                           if (crx1 /= 0.) crx1 = 1./crx1
                           curl%mgvar                           (i,   j,   k,   soln) = &
                                & (1. - Jacobi_damp)* curl%mgvar(i,   j,   k,   soln) - &
                                &       Jacobi_damp * curl%mgvar(i,   j,   k,   src)  * cr
                           if (has_dir(xdim))         curl%mgvar(i,   j,   k,   soln) = curl%mgvar(i,   j,   k,   soln)  + &
                                &       Jacobi_damp *(curl%mgvar(i-1, j,   k,   soln) + curl%mgvar(i+1, j,   k,   soln)) * crx + &
                                &       Jacobi_damp *(curl%mgvar(i+1, j,   k,   soln) - curl%mgvar(i-1, j,   k,   soln)) * crx * crx1
                           if (has_dir(ydim))         curl%mgvar(i,   j,   k,   soln) = curl%mgvar(i,   j,   k,   soln)  + &
                                &       Jacobi_damp *(curl%mgvar(i,   j-1, k,   soln) + curl%mgvar(i,   j+1, k,   soln)) * cry
                           if (has_dir(zdim))         curl%mgvar(i,   j,   k,   soln) = curl%mgvar(i,   j,   k,   soln)  + &
                                &       Jacobi_damp *(curl%mgvar(i,   j,   k-1, soln) + curl%mgvar(i,   j,   k+1, soln)) * crz
                        enddo
                     enddo
                  enddo
               case default
                  call die("[multigrid_gravity:approximate_solution_rbgs] Unsupported geometry.")
            end select
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

      use constants,          only: LO, HI, ndims, xdim, ydim, zdim, GEO_XYZ
      use grid,               only: D_x, D_y, D_z
      use mpisetup,           only: has_dir, eff_dim, geometry_type
      use dataio_pub,         only: die, warn
      use multigridhelpers,   only: dirty_debug, check_dirty, dirtyL, multidim_code_3D
      use multigridmpifuncs,  only: mpi_multigrid_bnd
      use multigridvars,      only: lvl, plvl, base, extbnd_antimirror, single_base

      implicit none

      integer, intent(in) :: lev  !< level for which approximate the solution
      integer, intent(in) :: src  !< index of source in lvl()%mgvar
      integer, intent(in) :: soln !< index of solution in lvl()%mgvar

      integer :: nf, n, nsmoo
      type(plvl), pointer :: curl

      curl => lvl(lev)

      if (curl%empty) return

      if (curl%fft_type == fft_none) call die("[multigrid_gravity:approximate_solution_fft] unknown FFT type")

      if (geometry_type /= GEO_XYZ) call die("[multigrid_gravity:approximate_solution_fft] FFT is not allowed in non-cartesian coordinates.")

      do nf = 1, nsmoof
         curl%src(:, :, :) = curl%mgvar(curl%is:curl%ie, curl%js:curl%je, curl%ks:curl%ke, src)

         if (curl%fft_type == fft_dst) then !correct boundaries on non-periodic local domain
            if (nf == 1 .and. .not. associated(curl, base)) then
               call make_face_boundaries(lev, soln)
            else
               call mpi_multigrid_bnd(lev, soln, 1, extbnd_antimirror)
               if (has_dir(xdim)) then
                  curl%bnd_x(:, :, LO) = 0.5* sum (curl%mgvar(curl%is-1:curl%is, curl%js:curl%je, curl%ks:curl%ke, soln), 1)
                  curl%bnd_x(:, :, HI) = 0.5* sum (curl%mgvar(curl%ie:curl%ie+1, curl%js:curl%je, curl%ks:curl%ke, soln), 1)
               endif
               if (has_dir(ydim)) then
                  curl%bnd_y(:, :, LO) = 0.5* sum (curl%mgvar(curl%is:curl%ie, curl%js-1:curl%js, curl%ks:curl%ke, soln), 2)
                  curl%bnd_y(:, :, HI) = 0.5* sum (curl%mgvar(curl%is:curl%ie, curl%je:curl%je+1, curl%ks:curl%ke, soln), 2)
               endif
               if (has_dir(zdim)) then
                  curl%bnd_z(:, :, LO) = 0.5* sum (curl%mgvar(curl%is:curl%ie, curl%js:curl%je, curl%ks-1:curl%ks, soln), 3)
                  curl%bnd_z(:, :, HI) = 0.5* sum (curl%mgvar(curl%is:curl%ie, curl%js:curl%je, curl%ke:curl%ke+1, soln), 3)
               endif
            endif

            if (dirty_debug) then
               if (has_dir(xdim) .and. any(abs(curl%bnd_x(:, :, :)) > dirtyL)) call warn("approximate_solution_fft dirty bnd_x")
               if (has_dir(ydim) .and. any(abs(curl%bnd_y(:, :, :)) > dirtyL)) call warn("approximate_solution_fft dirty bnd_y")
               if (has_dir(zdim) .and. any(abs(curl%bnd_z(:, :, :)) > dirtyL)) call warn("approximate_solution_fft dirty bnd_z")
            endif

            if (has_dir(xdim)) then
               curl%src(1,        :, :) = curl%src(1,        :, :) - curl%bnd_x(:, :, LO) * 2. * curl%idx2
               curl%src(curl%nxb, :, :) = curl%src(curl%nxb, :, :) - curl%bnd_x(:, :, HI) * 2. * curl%idx2
            endif
            if (has_dir(ydim)) then
               curl%src(:, 1,        :) = curl%src(:, 1,        :) - curl%bnd_y(:, :, LO) * 2. * curl%idy2
               curl%src(:, curl%nyb, :) = curl%src(:, curl%nyb, :) - curl%bnd_y(:, :, HI) * 2. * curl%idy2
            endif
            if (has_dir(zdim)) then
               curl%src(:, :, 1       ) = curl%src(:, :, 1       ) - curl%bnd_z(:, :, LO) * 2. * curl%idz2
               curl%src(:, :, curl%nzb) = curl%src(:, :, curl%nzb) - curl%bnd_z(:, :, HI) * 2. * curl%idz2
            endif
         endif

         call fft_convolve(lev)

         curl%mgvar(curl%is:curl%ie, curl%js:curl%je, curl%ks:curl%ke, soln) = curl%src(:, :, :)

         call check_dirty(lev, soln, "approx_soln fft+")

         !> \deprecated BEWARE use has_dir() here in a way that does not degrade performance

         if (associated(curl, base) .and. single_base) then !> do not relax if it is a whole domain (BEWARE: simplified check)
            nsmoo = 0
         else
            nsmoo = nsmool
         endif

         !relax the boundaries
         do n = 1, nsmoo
            call mpi_multigrid_bnd(lev, soln, 1, extbnd_antimirror)
            ! Possible optimization: This is a quite costly part of the local FFT solver
            if (fft_full_relax) then
               if (eff_dim == ndims .and. .not. multidim_code_3D) then
                  curl%mgvar                (curl%is:curl%ie,     curl%js:curl%je,     curl%ks:curl%ke,     soln)  = &
                       curl%rx * (curl%mgvar(curl%is-1:curl%ie-1, curl%js:curl%je,     curl%ks:curl%ke,     soln)  + &
                       &          curl%mgvar(curl%is+1:curl%ie+1, curl%js:curl%je,     curl%ks:curl%ke,     soln)) + &
                       curl%ry * (curl%mgvar(curl%is:curl%ie,     curl%js-1:curl%je-1, curl%ks:curl%ke,     soln)  + &
                       &          curl%mgvar(curl%is:curl%ie,     curl%js+1:curl%je+1, curl%ks:curl%ke,     soln)) + &
                       curl%rz * (curl%mgvar(curl%is:curl%ie,     curl%js:curl%je,     curl%ks-1:curl%ke-1, soln)  + &
                       &          curl%mgvar(curl%is:curl%ie,     curl%js:curl%je,     curl%ks+1:curl%ke+1, soln)) - &
                       curl%r  *  curl%mgvar(curl%is:curl%ie,     curl%js:curl%je,     curl%ks:curl%ke,     src)
               else

                  call die("[multigrid_gravity:approximate_solution_fft] fft_full_relax is allowed only for 3D at the moment")

                  ! An additional array (curl%src would be good enough) is required here to assemble partial results or use red-black passes
                  curl%mgvar                (curl%is:curl%ie,     curl%js:curl%je,     curl%ks:curl%ke,     soln)  = &
                       - curl%r * curl%mgvar(curl%is:curl%ie,     curl%js:curl%je,     curl%ks:curl%ke,     src)
                  if (has_dir(xdim)) &
                       &          curl%mgvar(curl%is:curl%ie,     curl%js:curl%je,     curl%ks:curl%ke,     soln)  = &
                       &          curl%mgvar(curl%is:curl%ie,     curl%js:curl%je,     curl%ks:curl%ke,     soln)  + &
                       curl%rx * (curl%mgvar(curl%is-1:curl%ie-1, curl%js:curl%je,     curl%ks:curl%ke,     soln)  + &
                       &          curl%mgvar(curl%is+1:curl%ie+1, curl%js:curl%je,     curl%ks:curl%ke,     soln))
                  if (has_dir(ydim)) &
                       &          curl%mgvar(curl%is:curl%ie,     curl%js:curl%je,     curl%ks:curl%ke,     soln)  = &
                       &          curl%mgvar(curl%is:curl%ie,     curl%js:curl%je,     curl%ks:curl%ke,     soln)  + &
                       curl%ry * (curl%mgvar(curl%is:curl%ie,     curl%js-1:curl%je-1, curl%ks:curl%ke,     soln)  + &
                       &          curl%mgvar(curl%is:curl%ie,     curl%js+1:curl%je+1, curl%ks:curl%ke,     soln))
                  if (has_dir(zdim)) &
                       &          curl%mgvar(curl%is:curl%ie,     curl%js:curl%je,     curl%ks:curl%ke,     soln)  = &
                       &          curl%mgvar(curl%is:curl%ie,     curl%js:curl%je,     curl%ks:curl%ke,     soln)  + &
                       curl%rz * (curl%mgvar(curl%is:curl%ie,     curl%js:curl%je,     curl%ks-1:curl%ke-1, soln)  + &
                       &          curl%mgvar(curl%is:curl%ie,     curl%js:curl%je,     curl%ks+1:curl%ke+1, soln))
               endif
            else
               ! relax only two layers of cells (1 is  significantly worse, 3 does not improve much)
               ! edges are relaxed twice, corners are relaxed three times which seems to be good

               if (has_dir(xdim)) then
                  curl%mgvar                (curl%is    :curl%is+D_x,   curl%js    :curl%je,     curl%ks    :curl%ke,     soln)  = & ! -X
                       curl%rx * (curl%mgvar(curl%is-D_x:curl%is,       curl%js    :curl%je,     curl%ks    :curl%ke,     soln)  + &
                       &          curl%mgvar(curl%is+D_x:curl%is+2*D_x, curl%js    :curl%je,     curl%ks    :curl%ke,     soln)) + &
                       curl%ry * (curl%mgvar(curl%is    :curl%is+D_x,   curl%js-D_y:curl%je-D_y, curl%ks    :curl%ke,     soln)  + &
                       &          curl%mgvar(curl%is    :curl%is+D_x,   curl%js+D_y:curl%je+D_y, curl%ks    :curl%ke,     soln)) + &
                       curl%rz * (curl%mgvar(curl%is    :curl%is+D_x,   curl%js    :curl%je,     curl%ks-D_z:curl%ke-D_z, soln)  + &
                       &          curl%mgvar(curl%is    :curl%is+D_x,   curl%js    :curl%je,     curl%ks+D_z:curl%ke+D_z, soln)) - &
                       curl%r  *  curl%mgvar(curl%is    :curl%is+D_x,   curl%js    :curl%je,     curl%ks    :curl%ke,     src)

                  curl%mgvar                (curl%ie-D_x  :curl%ie,     curl%js    :curl%je,     curl%ks    :curl%ke,     soln)  = & ! +X
                       curl%rx * (curl%mgvar(curl%ie-2*D_x:curl%ie-D_x, curl%js    :curl%je,     curl%ks    :curl%ke,     soln)  + &
                       &          curl%mgvar(curl%ie      :curl%ie+D_x, curl%js    :curl%je,     curl%ks    :curl%ke,     soln)) + &
                       curl%ry * (curl%mgvar(curl%ie-D_x  :curl%ie,     curl%js-D_y:curl%je-D_y, curl%ks    :curl%ke,     soln)  + &
                       &          curl%mgvar(curl%ie-D_x  :curl%ie,     curl%js+D_y:curl%je+D_y, curl%ks    :curl%ke,     soln)) + &
                       curl%rz * (curl%mgvar(curl%ie-D_x  :curl%ie,     curl%js    :curl%je,     curl%ks-D_z:curl%ke-D_z, soln)  + &
                       &          curl%mgvar(curl%ie-D_x  :curl%ie,     curl%js    :curl%je,     curl%ks+D_z:curl%ke+D_z, soln)) - &
                       curl%r  *  curl%mgvar(curl%ie-D_x  :curl%ie,     curl%js    :curl%je,     curl%ks    :curl%ke,     src)
               endif

               if (has_dir(ydim)) then
                  curl%mgvar                (curl%is    :curl%ie,     curl%js    :curl%js+D_y,   curl%ks    :curl%ke,     soln)  = & ! -Y
                       curl%rx * (curl%mgvar(curl%is-D_x:curl%ie-D_x, curl%js    :curl%js+D_y,   curl%ks    :curl%ke,     soln)  + &
                       &          curl%mgvar(curl%is+D_x:curl%ie+D_x, curl%js    :curl%js+D_y,   curl%ks    :curl%ke,     soln)) + &
                       curl%ry * (curl%mgvar(curl%is    :curl%ie,     curl%js-D_y:curl%js,       curl%ks    :curl%ke,     soln)  + &
                       &          curl%mgvar(curl%is    :curl%ie,     curl%js+D_y:curl%js+2*D_y, curl%ks    :curl%ke,     soln)) + &
                       curl%rz * (curl%mgvar(curl%is    :curl%ie,     curl%js    :curl%js+D_y,   curl%ks-D_z:curl%ke-D_z, soln)  + &
                       &          curl%mgvar(curl%is    :curl%ie,     curl%js    :curl%js+D_y,   curl%ks+D_z:curl%ke+D_z, soln)) - &
                       curl%r  *  curl%mgvar(curl%is    :curl%ie,     curl%js    :curl%js+D_y,   curl%ks    :curl%ke,     src)

                  curl%mgvar                (curl%is    :curl%ie,     curl%je-D_y  :curl%je,     curl%ks    :curl%ke,     soln)  = & ! +Y
                       curl%rx * (curl%mgvar(curl%is-D_x:curl%ie-D_x, curl%je-D_y  :curl%je,     curl%ks    :curl%ke,     soln)  + &
                       &          curl%mgvar(curl%is+D_x:curl%ie+D_x, curl%je-D_y  :curl%je,     curl%ks    :curl%ke,     soln)) + &
                       curl%ry * (curl%mgvar(curl%is    :curl%ie,     curl%je-2*D_y:curl%je-D_y, curl%ks    :curl%ke,     soln)  + &
                       &          curl%mgvar(curl%is    :curl%ie,     curl%je      :curl%je+D_y, curl%ks    :curl%ke,     soln)) + &
                       curl%rz * (curl%mgvar(curl%is    :curl%ie,     curl%je-D_y  :curl%je,     curl%ks-D_z:curl%ke-D_z, soln)  + &
                       &          curl%mgvar(curl%is    :curl%ie,     curl%je-D_y  :curl%je,     curl%ks+D_z:curl%ke+D_z, soln)) - &
                       curl%r  *  curl%mgvar(curl%is    :curl%ie,     curl%je-D_y  :curl%je,     curl%ks    :curl%ke,     src)
               endif

               if (has_dir(zdim)) then
                  curl%mgvar                (curl%is    :curl%ie,     curl%js    :curl%je,     curl%ks    :curl%ks+D_z,   soln)  = & ! -Z
                       curl%rx * (curl%mgvar(curl%is-D_x:curl%ie-D_x, curl%js    :curl%je,     curl%ks    :curl%ks+D_z,   soln)  + &
                       &          curl%mgvar(curl%is+D_x:curl%ie+D_x, curl%js    :curl%je,     curl%ks    :curl%ks+D_z,   soln)) + &
                       curl%ry * (curl%mgvar(curl%is    :curl%ie,     curl%js-D_y:curl%je-D_y, curl%ks    :curl%ks+D_z,   soln)  + &
                       &          curl%mgvar(curl%is    :curl%ie,     curl%js+D_y:curl%je+D_y, curl%ks    :curl%ks+D_z,   soln)) + &
                       curl%rz * (curl%mgvar(curl%is    :curl%ie,     curl%js    :curl%je,     curl%ks-D_z:curl%ks,       soln)  + &
                       &          curl%mgvar(curl%is    :curl%ie,     curl%js    :curl%je,     curl%ks+D_z:curl%ks+2*D_z, soln)) - &
                       curl%r  *  curl%mgvar(curl%is    :curl%ie,     curl%js    :curl%je,     curl%ks    :curl%ks+D_z,   src)

                  curl%mgvar                (curl%is    :curl%ie,     curl%js    :curl%je,     curl%ke-D_z  :curl%ke ,    soln)  = & ! +Z
                       curl%rx * (curl%mgvar(curl%is-D_x:curl%ie-D_x, curl%js    :curl%je,     curl%ke-D_z  :curl%ke,     soln)  + &
                       &          curl%mgvar(curl%is+D_x:curl%ie+D_x, curl%js    :curl%je,     curl%ke-D_z  :curl%ke,     soln)) + &
                       curl%ry * (curl%mgvar(curl%is    :curl%ie,     curl%js-D_y:curl%je-D_y, curl%ke-D_z  :curl%ke,     soln)  + &
                       &          curl%mgvar(curl%is    :curl%ie,     curl%js+D_y:curl%je+D_y, curl%ke-D_z  :curl%ke,     soln)) + &
                       curl%rz * (curl%mgvar(curl%is    :curl%ie,     curl%js    :curl%je,     curl%ke-2*D_z:curl%ke-D_z, soln)  + &
                       &          curl%mgvar(curl%is    :curl%ie,     curl%js    :curl%je,     curl%ke      :curl%ke+D_z, soln)) - &
                       curl%r  *  curl%mgvar(curl%is    :curl%ie,     curl%js    :curl%je,     curl%ke-D_z  :curl%ke,     src)
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
      use multigridvars,      only: bnd_givenval, bnd_periodic, base, single_base

      implicit none

      integer, intent(in) :: lev  !< level for which approximate the solution
      integer, intent(in) :: soln !< index of solution in lvl()%mgvar

      if (grav_bnd == bnd_periodic .and. (nproc == 1 .or. (lev == base%level .and. single_base) ) ) then
         call zero_boundaries(lev)
      else
         if (lev > base%level) then
            call prolong_faces(lev, soln)
         else
            if (grav_bnd /= bnd_givenval) call zero_boundaries(lev)
            call warn("m:mfb WTF?")
         endif
      endif

   end subroutine make_face_boundaries

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

   subroutine fft_convolve(lev)

      use dataio_pub,    only: die
      use multigridvars, only: lvl, plvl

      implicit none

      integer, intent(in) :: lev !< level at which make the convolution

      type(plvl), pointer :: curl

      curl => lvl(lev)

      ! do the convolution in Fourier space; curl%src(:,:,:) -> curl%fft{r}(:,:,:)
      call dfftw_execute(curl%planf)

      select case (curl%fft_type)
         case (fft_rcr)
            curl%fft  = curl%fft  * curl%Green3D
         case (fft_dst)
            curl%fftr = curl%fftr * curl%Green3D
         case default
            call die("[multigrid_gravity:fft_convolve] Unknown FFT type.")
      end select

      call dfftw_execute(curl%plani) ! curl%fft{r}(:,:,:) -> curl%src(:,:,:)

   end subroutine fft_convolve

end module multigrid_gravity
