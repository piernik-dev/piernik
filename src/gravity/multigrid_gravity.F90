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

   use constants,          only: cbuff_len, ndims
   use multigrid_vstats,   only: vcycle_stats
   use multigrid_old_soln, only: soln_history

   implicit none

   private
   public :: multigrid_grav_par, init_multigrid_grav, cleanup_multigrid_grav, multigrid_solve_grav, init_multigrid_grav_ext, residual

   include "fftw3.f"
   ! constants from fftw3.f
   !   integer, parameter :: FFTW_MEASURE=0, FFTW_PATIENT=32, FFTW_ESTIMATE=64
   !   integer, parameter :: FFTW_RODFT01=8, FFTW_RODFT10=9

   ! namelist parameters
   real               :: norm_tol                                     !< stop V-cycle iterations when the ratio of norms ||residual||/||source|| is below this value
   real               :: overrelax                                    !< overrealaxation factor (if < 1. then works as underrelaxation), provided for tests
   real, dimension(ndims) :: overrelax_xyz                            !< directional overrealaxation factor for fine tuning convergence ratio when cell spacing is not equal in all 3 directions. Use with care, patience and lots of hope.
   real               :: Jacobi_damp                                  !< omega factor for damped Jacobi relaxation. Jacobi_damp == 1 gives undamped method. Try 0.5 in 1D.
   real               :: vcycle_abort                                 !< abort the V-cycle when lhs norm raises by this factor
   real               :: L4_strength                                  !< strength of the 4th order terms in the Laplace operator; 0.: 2nd, 1.: 4th direct, 0.5: 4th integral
   integer(kind=4)    :: max_cycles                                   !< Maximum allowed number of V-cycles
   integer(kind=4)    :: nsmoob                                       !< smoothing cycles on coarsest level when cannot use FFT. (a convergence check would be much better)
   integer(kind=4)    :: ord_laplacian                                !< Laplace operator order; allowed values are 2 (default) and 4 (experimental, not fully implemented)
   logical            :: trust_fft_solution                           !< Bypass the V-cycle, when doing FFT on whole domain, make sure first that FFT is properly set up.
   logical            :: base_no_fft                                  !< Deny solving the coarsest level with FFT. Can be very slow.
   logical            :: prefer_rbgs_relaxation                       !< Prefer relaxation over FFT local solver. Typically faster.
   !> \todo allow to perform one or more V-cycles with FFT method, the switch to the RBGS (may save one V-cycle in some cases)
   logical            :: fft_patient                                  !< Spend more time in init_multigrid to find faster fft plan
   character(len=cbuff_len) :: grav_bnd_str                           !< Type of gravitational boundary conditions.
   logical            :: require_FFT                                  !< .true. if we use FFT solver anywhere (and need face prolongation)

   integer            :: fftw_flags = FFTW_MEASURE                    !< or FFTW_PATIENT on request

   ! solution recycling
   type(soln_history), target :: inner, outer                         !< storage for recycling the inner and outer potentials

   ! miscellaneous
   type(vcycle_stats) :: vstat                                        !< V-cycle statistics

contains

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
!! <tr><td>nsmool                </td><td>4      </td><td>integer value  </td><td>\copydoc multigridvars::nsmool                    </td></tr>
!! <tr><td>nsmoob                </td><td>100    </td><td>integer value  </td><td>\copydoc multigrid_gravity::nsmoob                </td></tr>
!! <tr><td>overrelax             </td><td>1.     </td><td>real value     </td><td>\copydoc multigrid_gravity::overrelax             </td></tr>
!! <tr><td>overrelax_xxyz(ndims) </td><td>1.     </td><td>real value     </td><td>\copydoc multigrid_gravity::overrelax_xyz         </td></tr>
!! <tr><td>Jacobi_damp           </td><td>1.     </td><td>real value     </td><td>\copydoc multigrid_gravity::jacobi_damp           </td></tr>
!! <tr><td>L4_strength           </td><td>1.0    </td><td>real value     </td><td>\copydoc multigrid_gravity::l4_strength           </td></tr>
!! <tr><td>nsmoof                </td><td>1      </td><td>integer value  </td><td>\copydoc multigridvars::nsmoof                    </td></tr>
!! <tr><td>ord_laplacian         </td><td>2      </td><td>integer value  </td><td>\copydoc multigrid_gravity::ord_laplacian         </td></tr>
!! <tr><td>ord_time_extrap       </td><td>1      </td><td>integer value  </td><td>\copydoc multigrid_gravity::ord_time_extrap       </td></tr>
!! <tr><td>prefer_rbgs_relaxation</td><td>.true. </td><td>logical        </td><td>\copydoc multigrid_gravity::prefer_rbgs_relaxation</td></tr>
!! <tr><td>base_no_fft           </td><td>.false.</td><td>logical        </td><td>\copydoc multigrid_gravity::base_no_fft           </td></tr>
!! <tr><td>fft_full_relax        </td><td>.false.</td><td>logical        </td><td>\copydoc multigridvars::fft_full_relax            </td></tr>
!! <tr><td>fft_patient           </td><td>.false.</td><td>logical        </td><td>\copydoc multigrid_gravity::fft_patient           </td></tr>
!! <tr><td>trust_fft_solution    </td><td>.false.</td><td>logical        </td><td>\copydoc multigrid_gravity::trust_fft_solution    </td></tr>
!! <tr><td>coarsen_multipole     </td><td>1      </td><td>integer value  </td><td>\copydoc multipole::coarsen_multipole             </td></tr>
!! <tr><td>lmax                  </td><td>16     </td><td>integer value  </td><td>\copydoc multipole::lmax                          </td></tr>
!! <tr><td>mmax                  </td><td>-1     </td><td>integer value  </td><td>\copydoc multipole::mmax                          </td></tr>
!! <tr><td>ord_prolong_mpole     </td><td>-2     </td><td>integer value  </td><td>\copydoc multipole::ord_prolong_mpole             </td></tr>
!! <tr><td>use_point_monopole    </td><td>.false.</td><td>logical        </td><td>\copydoc multipole::use_point_monopole            </td></tr>
!! <tr><td>interp_pt2mom         </td><td>.false.</td><td>logical        </td><td>\copydoc multipole::interp_pt2mom                 </td></tr>
!! <tr><td>interp_mom2pot        </td><td>.false.</td><td>logical        </td><td>\copydoc multipole::interp_mom2pot                </td></tr>
!! <tr><td>multidim_code_3D      </td><td>.false.</td><td>logical        </td><td>\copydoc multigridvars::multidim_code_3d          </td></tr>
!! <tr><td>grav_bnd_str          </td><td>"periodic"/"dirichlet"</td><td>string of chars</td><td>\copydoc multigrid_gravity::grav_bnd_str          </td></tr>
!! </table>
!! The list is active while \b "GRAV" and \b "MULTIGRID" are defined.
!! \n \n
!<
   subroutine multigrid_grav_par

      use constants,     only: GEO_XYZ, GEO_RPZ, BND_PER, O_LIN, O_D2, O_I2, I_ONE
      use dataio_pub,    only: nh  ! QA_WARN required for diff_nml
      use dataio_pub,    only: msg, die, warn
      use domain,        only: dom, is_multicg !, is_uneven
      use mpisetup,      only: master, slave, ibuff, cbuff, rbuff, lbuff, nproc, piernik_MPI_Bcast
      use multigridvars, only: single_base, bnd_invalid, bnd_isolated, bnd_periodic, bnd_dirichlet, grav_bnd, fft_full_relax, multidim_code_3D, nsmool, nsmoof
      use multigrid_old_soln, only: nold_max, ord_time_extrap
      use multipole,     only: use_point_monopole, lmax, mmax, ord_prolong_mpole, coarsen_multipole, interp_pt2mom, interp_mom2pot

      implicit none

      integer       :: periodic_bnd_cnt   !< counter of periodic boundaries in existing directions
      logical, save :: frun = .true.      !< First run flag

      namelist /MULTIGRID_GRAVITY/ norm_tol, vcycle_abort, max_cycles, nsmool, nsmoob, &
           &                       overrelax, overrelax_xyz, Jacobi_damp, L4_strength, nsmoof, ord_laplacian, ord_time_extrap, &
           &                       prefer_rbgs_relaxation, base_no_fft, fft_full_relax, fft_patient, trust_fft_solution, &
           &                       coarsen_multipole, lmax, mmax, ord_prolong_mpole, use_point_monopole, interp_pt2mom, interp_mom2pot, multidim_code_3D, &
           &                       grav_bnd_str

      if (.not.frun) call die("[multigrid_gravity:multigrid_grav_par] Called more than once.")
      frun = .false.

      ! Default values for namelist variables
      norm_tol               = 1.e-6
      overrelax              = 1.
      overrelax_xyz(:)       = 1.
      Jacobi_damp            = 1.
      vcycle_abort           = 2.
      L4_strength            = 1.0

      coarsen_multipole      = 0
!      if (is_uneven) coarsen_multipole = 0
      lmax                   = 16
      mmax                   = -1 ! will be automatically set to lmax unless explicitly limited in problem.par
      max_cycles             = 20
      nsmool                 = 4
      nsmoob                 = 100
      nsmoof                 = 1
      ord_laplacian          = O_I2
      ord_prolong_mpole      = O_D2
      ord_time_extrap        = O_LIN

      use_point_monopole     = .false.
      trust_fft_solution     = .true.
      base_no_fft            = .true. !> \todo restore .false. when it will be possible
      prefer_rbgs_relaxation = .true. !> \warning it seems that FFT local solver is broken somewhere. \todo restore .false. as a default
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
            call warn("[multigrid_gravity:multigrid_grav_par] cylindrical geometry support is under development.")
            ! switch off FFT-related bits
            base_no_fft = .true.
            prefer_rbgs_relaxation = .true.
            ord_laplacian = O_I2
            L4_strength = 0.
            ! ord_prolong_mpole = O_INJ
         else if (dom%geometry_type /= GEO_XYZ) then
            call die("[multigrid_gravity:multigrid_grav_par] non-cartesian geometry not implemented yet.")
         endif

         if (is_multicg .and. .not. base_no_fft) then
            call warn("[multigrid_gravity:multigrid_grav_par] base_no_fft set to .true. for multicg configuration")
            base_no_fft = .true.
         endif

         if (.not. prefer_rbgs_relaxation .and. nproc > 1) then
            prefer_rbgs_relaxation = .true.
            call warn("[multigrid_gravity:multigrid_grav_par] FFT local solver disabled for multithreaded runs due to unresolved incompatibilities with mew grid features")
         endif

         rbuff(1) = norm_tol
         rbuff(2) = overrelax
         rbuff(3:5) = overrelax_xyz
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

      call piernik_MPI_Bcast(cbuff, cbuff_len)
      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(rbuff)
      call piernik_MPI_Bcast(lbuff)

      if (slave) then

         norm_tol       = rbuff(1)
         overrelax      = rbuff(2)
         overrelax_xyz  = rbuff(3:5)
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
                 call die("[multigrid_gravity:multigrid_grav_par] cannot enforce periodic boundaries for gravity on a not fully periodic domain")
            grav_bnd = bnd_periodic
         case ("dirichlet", "dir")
            grav_bnd = bnd_dirichlet
         case default
            call die("[multigrid_gravity:multigrid_grav_par] Non-recognized boundary description.")
      end select

      if (periodic_bnd_cnt == dom%eff_dim) then ! fully periodic domain
         if (grav_bnd /= bnd_periodic .and. master) call warn("[multigrid_gravity:multigrid_grav_par] Ignoring non-periodic boundary conditions for gravity on a fully periodic domain.")
         grav_bnd = bnd_periodic
      else if (periodic_bnd_cnt > 0 .and. periodic_bnd_cnt < dom%eff_dim) then
         if (master) call warn("[multigrid_gravity:multigrid_grav_par] Mixing periodic and non-periodic boundary conditions for gravity is experimental.")
         prefer_rbgs_relaxation = .true.
         base_no_fft = .true.
      endif
!!$      select case (grav_bnd)
!!$         case (bnd_periodic)
!!$            grav_extbnd_mode = BND_NONE
!!$         case (bnd_isolated, bnd_dirichlet, bnd_givenval)
!!$            grav_extbnd_mode = BND_NEGREF
!!$         case default
!!$            call die("[multigrid_gravity:multigrid_grav_par] Unsupported grav_bnd.")
!!$            !grav_extbnd_mode = BND_NONE
!!$      end select

      if (.not. (grav_bnd == bnd_periodic .or. grav_bnd == bnd_dirichlet .or. grav_bnd == bnd_isolated) .and. .not. base_no_fft) then
         base_no_fft = .true.
         if (master) call warn("[multigrid_gravity:multigrid_grav_par] Use of FFT not allowed by current boundary type/combination.")
      endif

      single_base = .not. base_no_fft

      if (.not. prefer_rbgs_relaxation .and. any([ overrelax, overrelax_xyz(:) ] /= 1.)) then
         if (master) call warn("[multigrid_gravity:multigrid_grav_par] Overrelaxation is disabled for FFT local solver.")
         overrelax = 1.
         overrelax_xyz(:) = 1.
      endif

      if (master) then
         if ((Jacobi_damp <= 0. .or. Jacobi_damp>1.)) then
            write(msg, '(a,g12.5,a)')"[multigrid_gravity:multigrid_grav_par] Jacobi_damp = ",Jacobi_damp," is outside (0, 1] interval."
            call warn(msg)
         endif
         if (overrelax /= 1. .or. any(overrelax_xyz(:) /= 1.)) then
            write(msg, '(a,f8.5,a,3f8.5,a)')"[multigrid_gravity:multigrid_grav_par] Overrelaxation factors: global = ", overrelax, ", directional = [", overrelax_xyz(:), "]"
            call warn(msg)
         endif
      endif

      if (fft_patient) fftw_flags = FFTW_PATIENT

      ! solution recycling
      ord_time_extrap = min(nold_max-I_ONE, max(-I_ONE, ord_time_extrap))
      associate (nold => ord_time_extrap + 1)
      if (nold > 0) then
         call inner%init_history(nold, "i")
         if (grav_bnd == bnd_isolated) call outer%init_history(nold, "o")
      endif
      end associate

      call vstat%init(max_cycles)

   end subroutine multigrid_grav_par

!> \brief Initialization - continued after allocation of everything interesting

   subroutine init_multigrid_grav

      use cg_leaves,           only: leaves
      use cg_level_coarsest,   only: coarsest
      use cg_level_connected,  only: cg_level_connected_T
      use cg_level_finest,     only: finest
      use constants,           only: GEO_XYZ, sgp_n, fft_none, fft_dst, fft_rcr, dsetnamelen
      use dataio_pub,          only: die, warn, printinfo, msg
      use domain,              only: dom
      use mpisetup,            only: master, nproc
      use multigrid_fftapprox, only: mpi_multigrid_prep_grav
      use multigridvars,       only: bnd_periodic, bnd_dirichlet, bnd_isolated, grav_bnd
      use multipole,           only: init_multipole, coarsen_multipole
      use named_array_list,    only: qna

      implicit none

      type(cg_level_connected_T), pointer :: curl
      character(len=dsetnamelen) :: FFTn

      if (coarsen_multipole /= 0) then
         coarsen_multipole = 0
         if (master) call warn("[multigrid_gravity:init_multigrid_grav] multipole coarsening temporarily disabled")
      endif

      call leaves%set_q_value(qna%ind(sgp_n), 0.) !Initialize all the guardcells, even those which does not impact the solution

      curl => coarsest%level
      do while (associated(curl))

         if (prefer_rbgs_relaxation) then
            curl%fft_type = fft_none
         else if (grav_bnd == bnd_periodic .and. nproc == 1) then
            curl%fft_type = fft_rcr
         else if (grav_bnd == bnd_periodic .or. grav_bnd == bnd_dirichlet .or. grav_bnd == bnd_isolated) then
            curl%fft_type = fft_dst
         else
            curl%fft_type = fft_none
         endif

         curl => curl%finer
      enddo

      base_no_fft = base_no_fft .or. (coarsest%level%tot_se /= 1)

      ! data related to local and global base-level FFT solver
      if (base_no_fft) then
         coarsest%level%fft_type = fft_none
      else
         select case (grav_bnd)
         case (bnd_periodic)
            coarsest%level%fft_type = fft_rcr
         case (bnd_dirichlet, bnd_isolated)
            coarsest%level%fft_type = fft_dst
         case default
            coarsest%level%fft_type = fft_none
            if (master) call warn("[multigrid_gravity:init_multigrid_grav] base_no_fft unset but no suitable boundary conditions found. Reverting to RBGS relaxation.")
         end select
      endif

      require_FFT = .false.

      ! FFT solver storage and data
      curl => coarsest%level
      do while (associated(curl))

         if (curl%fft_type /= fft_none) then
            require_FFT = .true.
            if (dom%geometry_type /= GEO_XYZ) call die("[multigrid_gravity:init_multigrid_grav] FFT is not allowed in non-cartesian coordinates.")
         endif

         if (master) then
            select case (curl%fft_type)
               case (fft_rcr)
                  FFTn="RCR"
               case (fft_dst)
                  FFTn="DST"
               case (fft_none)
                  FFTn="none"
               case default
                  call die("[multigrid_gravity:init_multigrid_grav] Unknown FFT type.")
            end select
            write(msg,'(a,i3,2a)')"[multigrid_gravity:init_multigrid_grav] Level ",curl%level_id,", FFT: ", trim(FFTn)
            call printinfo(msg)
         endif

         curl => curl%finer
      enddo

      if (require_FFT) call mpi_multigrid_prep_grav !supplement to mpi_multigrid_prep

      if (finest%level%fft_type == fft_none .and. trust_fft_solution) then
         if (master) call warn("[multigrid_gravity:init_multigrid_grav] cannot trust FFT solution on the finest level.")
         trust_fft_solution = .false.
      endif

      if (grav_bnd == bnd_isolated) call init_multipole

   end subroutine init_multigrid_grav

!> \brief Cleanup

   subroutine cleanup_multigrid_grav

      use multipole,      only: cleanup_multipole

      implicit none

      call cleanup_multipole
      call vstat%cleanup
      call dfftw_cleanup
      call inner%cleanup_history
      call outer%cleanup_history

   end subroutine cleanup_multigrid_grav

!> set up pointers for cg%mg initialization

   subroutine init_multigrid_grav_ext(after_label)

      use grid_container_ext, only: cg_ext, cg_extptrs

      implicit none

      character(len=*), intent(in) :: after_label

      procedure(cg_ext), pointer :: mgg_cg_init_p, mgg_cg_cleanup_p

      mgg_cg_init_p    => mgg_cg_init
      mgg_cg_cleanup_p => mgg_cg_cleanup
      call cg_extptrs%extend(mgg_cg_init_p, mgg_cg_cleanup_p, "multigrid_gravity", after_label)

   end subroutine init_multigrid_grav_ext

!> \brief Allocate some multigrid-specific arrays

   subroutine mgg_cg_init(cg)

      use cg_level_coarsest,  only: coarsest
      use cg_level_connected, only: cg_level_connected_T
      use constants,          only: xdim, ydim, zdim, fft_rcr, fft_dst, fft_none, pi, dpi, zero, half, one
      use dataio_pub,         only: die
      use domain,             only: dom
      use grid_cont,          only: grid_container

      implicit none

      type(grid_container), pointer,  intent(inout) :: cg
      type(cg_level_connected_T), pointer :: curl
      real, allocatable, dimension(:)  :: kx, ky, kz             !< FFT kernel directional components for convolution
      integer :: i, j

      ! this should work correctly also when dom%eff_dim < 3
      cg%mg%r  = overrelax / 2.
      cg%mg%rx = cg%dvol2 * cg%idx2
      cg%mg%ry = cg%dvol2 * cg%idy2
      cg%mg%rz = cg%dvol2 * cg%idz2
      cg%mg%r  = cg%mg%r / (cg%mg%rx + cg%mg%ry + cg%mg%rz)
      cg%mg%rx = overrelax_xyz(xdim) * cg%mg%rx * cg%mg%r
      cg%mg%ry = overrelax_xyz(ydim) * cg%mg%ry * cg%mg%r
      cg%mg%rz = overrelax_xyz(zdim) * cg%mg%rz * cg%mg%r
      cg%mg%r  = cg%mg%r * cg%dvol2
      !>
      !! \deprecated BEWARE: some of the above invariants may be not optimally defined - the convergence ratio drops when dx /= dy or dy /= dz or dx /= dz
      !! and overrelaxation factors are required to get any convergence (often poor)
      !<

      ! FFT solver storage and data
      curl => coarsest%level
      do while (associated(curl))
         if (cg%level_id == curl%level_id) exit
         curl => curl%finer
      enddo

      if (.not. associated(curl)) call die("[multigrid_gravity:mgg_cg_init] level not found")
      if (cg%level_id /= curl%level_id) call die("[multigrid_gravity:mgg_cg_init] wrong level found")

      cg%mg%planf = 0
      cg%mg%plani = 0

      if (curl%fft_type /= fft_none) then

         select case (curl%fft_type)
            case (fft_rcr)
               cg%mg%nxc = cg%nxb / 2 + 1
            case (fft_dst)
               cg%mg%nxc = cg%nxb
            case default
               call die("[multigrid_gravity:mgg_cg_init] Unknown FFT type.")
         end select

         if (allocated(cg%mg%Green3D) .or. allocated(cg%mg%src)) call die("[multigrid_gravity:mgg_cg_init] Green3D or src arrays already allocated")
         allocate(cg%mg%Green3D(cg%mg%nxc, cg%nyb, cg%nzb))
         allocate(cg%mg%src    (cg%nxb,    cg%nyb, cg%nzb))

         allocate(kx(cg%mg%nxc), ky(cg%nyb), kz(cg%nzb))

         select case (curl%fft_type)

            ! cg%mg%fft_norm is set such that the following sequence gives identity:
            ! call dfftw_execute(cg%mg%planf); cg%mg%fftr(:, :, :) = cg%mg%fftr(:, :, :) * cg%mg%fft_norm ; call dfftw_execute(cg%mg%plani)

            case (fft_rcr)
               if (allocated(cg%mg%fft)) call die("[multigrid_gravity:mgg_cg_init] fft or Green3D array already allocated")
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

               if (allocated(cg%mg%fftr)) call die("[multigrid_gravity:mgg_cg_init] fftr array already allocated")
               allocate(cg%mg%fftr(cg%mg%nxc, cg%nyb, cg%nzb))

               cg%mg%fft_norm = one / (8. * real( product(cg%n_b(:), mask=dom%has_dir(:)) ))
               kx(:) = cg%idx2 * (cos(pi/cg%nxb*[( j, j=1, cg%mg%nxc )]) - one)
               ky(:) = cg%idy2 * (cos(pi/cg%nyb*[( j, j=1, cg%nyb )]) - one)
               kz(:) = cg%idz2 * (cos(pi/cg%nzb*[( j, j=1, cg%nzb )]) - one)
               call dfftw_plan_r2r_3d(cg%mg%planf, cg%nxb, cg%nyb, cg%nzb, cg%mg%src,  cg%mg%fftr, FFTW_RODFT10, FFTW_RODFT10, FFTW_RODFT10, fftw_flags)
               call dfftw_plan_r2r_3d(cg%mg%plani, cg%nxb, cg%nyb, cg%nzb, cg%mg%fftr, cg%mg%src,  FFTW_RODFT01, FFTW_RODFT01, FFTW_RODFT01, fftw_flags)

            case default
               call die("[multigrid_gravity:mgg_cg_init] Unknown FFT type.")
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

   end subroutine mgg_cg_init

!> \brief Deallocate what was allocated in mg_cg_init

   subroutine mgg_cg_cleanup(cg)

      use grid_cont, only: grid_container
!!$      use constants,     only: LO, HI, ndims
!!$      use grid_cont,   only: tgt_list

      implicit none

      type(grid_container), pointer,  intent(inout) :: cg
!!$      integer :: g, ib
!!$      integer, parameter :: nseg = 2*(HI-LO+1)*ndims
!!$      type(tgt_list), dimension(nseg) :: io_tgt

      if (allocated(cg%mg%fft))     deallocate(cg%mg%fft)
      if (allocated(cg%mg%fftr))    deallocate(cg%mg%fftr)
      if (allocated(cg%mg%src))     deallocate(cg%mg%src)
      if (allocated(cg%mg%Green3D)) deallocate(cg%mg%Green3D)

      if (cg%mg%planf /= 0) call dfftw_destroy_plan(cg%mg%planf)
      if (cg%mg%plani /= 0) call dfftw_destroy_plan(cg%mg%plani)

!!$         io_tgt(1:nseg) = [ cgl%cg%mg%pfc_tgt, cgl%cg%mg%pff_tgt ]
!!$            do ib = 1, nseg
!!$               if (allocated(io_tgt(ib)%seg)) then
!!$                  do g = lbound(io_tgt(ib)%seg, dim=1), ubound(io_tgt(ib)%seg, dim=1)
!!$                     if (allocated(io_tgt(ib)%seg(g)%buf)) deallocate(io_tgt(ib)%seg(g)%buf)
!!$                     if (allocated(io_tgt(ib)%seg(g)%f_lay)) deallocate(io_tgt(ib)%seg(g)%f_lay)
!!$                  enddo
!!$                  deallocate(io_tgt(ib)%seg)
!!$               endif
!!$            enddo

   end subroutine mgg_cg_cleanup

!>
!! \brief Make a local copy of source (density) and multiply by 4 pi G
!!
!! \details Typically i_all_dens is a copy of fluidindex::iarr_all_sg.
!! Passing this as an argument allows for independent computation of the potential for several density fields if necessary.
!! \todo compactify the following more (if possible)
!<

   subroutine init_source(i_all_dens)

      use cg_list_global, only: all_cg
      use constants,      only: GEO_RPZ, LO, HI, xdim, ydim, zdim
      use dataio_pub,     only: die
      use domain,         only: dom
      use cg_list,        only: cg_list_element
      use cg_leaves,      only: leaves
      use grid_cont,      only: grid_container
      use multigridvars,  only: source, bnd_periodic, bnd_dirichlet, bnd_givenval, grav_bnd
      use units,          only: fpiG
      use particle_pub,   only: pset
#ifdef JEANS_PROBLEM
      use problem_pub,    only: jeans_d0, jeans_mode ! hack for tests
#endif /* JEANS_PROBLEM */

      implicit none

      integer(kind=4), dimension(:), intent(in) :: i_all_dens !< indices to selfgravitating fluids

      real                           :: fac
      integer                        :: i, side
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg

      call all_cg%set_dirty(source)

      if (size(i_all_dens) > 0) then
         cgl => leaves%first
         do while (associated(cgl))
            cg => cgl%cg
            cg%q(source)%arr(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) = fpiG * sum(cg%u(i_all_dens, cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke), dim=1)
            cgl => cgl%nxt
         enddo
         call pset%map(source, fpiG)
      else
         call leaves%set_q_value(source, 0.)
      endif

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
               do side = LO, HI
                  if (cg%ext_bnd(xdim, side)) then
                     fac = 2. * cg%idx2 / fpiG
                     if (dom%geometry_type == GEO_RPZ .and. cg%x(cg%ijkse(xdim,side)) /= 0.) fac = fac - 1./(cg%dx * cg%x(cg%ijkse(xdim,side)) * fpiG) !> BEWARE is it cg%x(ie), cg%x(ie+1) or something in the middle?
                     cg%q(source)%arr       (cg%ijkse(xdim,side), cg%js:cg%je, cg%ks:cg%ke) = &
                          & cg%q(source)%arr(cg%ijkse(xdim,side), cg%js:cg%je, cg%ks:cg%ke) - &
                          & cg%mg%bnd_x(                          cg%js:cg%je, cg%ks:cg%ke, side) * fac
                  endif
               enddo
               do side = LO, HI
                  if (cg%ext_bnd(ydim, side)) then
                     if (dom%geometry_type == GEO_RPZ) then
                        do i = cg%is, cg%ie
                           if (cg%x(i) /= 0.) cg%q(source)%arr(i, cg%ijkse(ydim,side), cg%ks:cg%ke) = &
                                &             cg%q(source)%arr(i, cg%ijkse(ydim,side), cg%ks:cg%ke) - &
                                &             cg%mg%bnd_y     (i,                      cg%ks:cg%ke, side) * 2. * cg%idy2 / fpiG / cg%x(i)**2
                        enddo
                     else
                        cg%q(source)%arr     (cg%is:cg%ie, cg%ijkse(ydim,side), cg%ks:cg%ke) = &
                           & cg%q(source)%arr(cg%is:cg%ie, cg%ijkse(ydim,side), cg%ks:cg%ke) - &
                           & cg%mg%bnd_y     (cg%is:cg%ie,                      cg%ks:cg%ke, side) * 2. * cg%idy2 / fpiG
                     endif
                  endif
               enddo
               do side = LO, HI
                  if (cg%ext_bnd(zdim, side)) cg%q(source)%arr(cg%is:cg%ie, cg%js:cg%je, cg%ijkse(zdim,side)) = &
                       &                      cg%q(source)%arr(cg%is:cg%ie, cg%js:cg%je, cg%ijkse(zdim,side)) - &
                       &                      cg%mg%bnd_z     (cg%is:cg%ie, cg%js:cg%je, side) * 2. * cg%idz2 / fpiG
               enddo
               cgl => cgl%nxt
            enddo

         case default
            call die("[multigrid_gravity:init_source] Unknown boundary type")
      end select

      call leaves%check_dirty(source, "init_src")

   end subroutine init_source

!>
!! \brief Multigrid gravity driver. This is the only multigrid routine intended to be called from the gravity module.
!! This routine is also responsible for communicating the solution to the rest of world via sgp array.
!<

   subroutine multigrid_solve_grav(i_all_dens)

      use cg_leaves,        only: leaves
      use constants,        only: sgp_n
      use multigridvars,    only: solution, tot_ts, ts, grav_bnd, bnd_dirichlet, bnd_givenval, bnd_isolated, all_dirty
      use multipole,        only: multipole_solver
      use named_array_list, only: qna
      use timer,            only: set_timer

      implicit none

      integer(kind=4), dimension(:), intent(in) :: i_all_dens !< indices to selfgravitating fluids

      integer :: grav_bnd_global
      integer(kind=4), dimension(0) :: empty_array !< trick to avoid compiler warnings on possibly uninitialized i_all_dens.0 in init_source

      ts =  set_timer("multigrid", .true.)

      call all_dirty

      grav_bnd_global = grav_bnd

      if (grav_bnd_global == bnd_isolated) then
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

      call leaves%q_copy(solution, qna%ind(sgp_n))

      if (grav_bnd_global == bnd_isolated) then
         grav_bnd = bnd_givenval

         vstat%cprefix = "Go-"
         call multipole_solver
         call init_source(empty_array)

         call vcycle_hg(outer)

         call leaves%q_add(solution, qna%ind(sgp_n)) ! add solution to sgp

      endif

      grav_bnd = grav_bnd_global
      ts = set_timer("multigrid")
      tot_ts = tot_ts + ts

   end subroutine multigrid_solve_grav

!>
!! \brief The solver. Here we choose an adaptation of the Huang-Greengard V-cycle.
!! For more difficult problems, like variable coefficient diffusion equation a more sophisticated V-cycle may be more effective.
!<

   subroutine vcycle_hg(history)

      use cg_leaves,           only: leaves
      use cg_list_global,      only: all_cg
      use cg_level_coarsest,   only: coarsest
      use cg_level_connected,  only: cg_level_connected_T
      use cg_level_finest,     only: finest
      use constants,           only: cbuff_len, fft_none
      use dataio_pub,          only: msg, die, warn, printinfo
      use global,              only: do_ascii_dump
      use mpisetup,            only: master, nproc
      use multigridvars,       only: source, solution, correction, defect, verbose_vcycle, stdout, tot_ts, ts, grav_bnd, bnd_periodic
      use multigrid_old_soln,  only: soln_history
      use timer,               only: set_timer

      implicit none

      type(soln_history), intent(inout) :: history !< inner or outer potential history used for initializing first guess

      integer :: v
      real    :: norm_rhs, norm_lhs, norm_old, norm_lowest
      logical :: dump_every_step, dump_result
      logical, save      :: norm_was_zero = .false.
      real,    parameter :: suspicious_factor = 1.05 !> \deprecated If the norm decreases too slowly then dump diagnostic output (BEWARE: this option is for tests only)
      integer, parameter :: fmtlen = 32
      character(len=fmtlen)    :: fmt
      character(len=cbuff_len) :: dname
      integer(kind=4), dimension(4)    :: mg_fields
      type(cg_level_connected_T), pointer :: curl

      inquire(file = "_dump_every_step_", EXIST=dump_every_step) ! use for debug only
      inquire(file = "_dump_result_", EXIST=dump_result)
      write(dname,'(2a)')trim(vstat%cprefix),"mdump"
      mg_fields = [ source, solution, defect, correction ]

      do_ascii_dump = do_ascii_dump .or. dump_every_step .or. dump_result

      ! On single CPU use FFT if possible because it is faster. Can be disabled by prefer_rbgs_relaxation = .true.
      if (nproc == 1 .and. finest%level%fft_type /= fft_none) then
         call all_cg%set_dirty(solution)
         call fft_solve_roof
         if (trust_fft_solution) then
            write(msg, '(3a)')"[multigrid_gravity:vcycle_hg] FFT solution trusted, skipping ", trim(vstat%cprefix), "cycle."
            call printinfo(msg, stdout)
            return
         endif
      else
         call history%init_solution(vstat%cprefix)
      endif

      norm_lhs = 0.
      norm_rhs = leaves%norm_sq(source)
      norm_old = norm_rhs
      norm_lowest = norm_rhs

      if (norm_rhs == 0.) then ! empty domain => potential == 0.
         if (master .and. .not. norm_was_zero) call warn("[multigrid_gravity:vcycle_hg] No gravitational potential for an empty space.")
         norm_was_zero = .true.
         call history%store_solution
         return
      else
         if (master .and. norm_was_zero) call warn("[multigrid_gravity:vcycle_hg] Spontaneous mass creation detected!")
         norm_was_zero = .false.
      endif

      ! if (.not. history%valid .and. prefer_rbgs_relaxation) call approximate_solution(finest, source, solution)
      ! not necessary when init_solution called FFT
      ! In single-process run, the above call may often save one V-cycle, but it always costs almost like one V-cycle.
      ! \todo test it in massively-parallel run (and most probably throw out for simplicity)
      ! difficult statement: for approximate_solution_fft it requires to pass a flag to use guardcells instead of prolonging faces.

      ! iterations
      do v = 0, max_cycles

         call all_cg%set_dirty(defect)
         call residual(leaves, source, solution, defect)
         call leaves%check_dirty(defect, "residual")
         if (grav_bnd == bnd_periodic) call leaves%subtract_average(defect)

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

         if (v>0 .and. norm_old/norm_lhs <= suspicious_factor) call all_cg%numbered_ascii_dump(mg_fields, dname, v)
         if (dump_result .and. norm_lhs/norm_rhs <= norm_tol) call all_cg%numbered_ascii_dump(mg_fields, dname)

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
         call finest%level%restrict_to_floor_q_1var(defect)

         call all_cg%set_dirty(correction)
         call coarsest%level%set_q_value(correction, 0.)

         curl => coarsest%level
         do while (associated(curl))
            call approximate_solution(curl, defect, correction)
            call curl%check_dirty(correction, "Vup relax+")
            curl => curl%finer
         enddo
         call leaves%q_add(correction, solution)

         if (dump_every_step) call all_cg%numbered_ascii_dump(mg_fields, dname, v)

      enddo

#ifdef DEBUG
      call residual(leaves, source, solution, defect)
#endif /* DEBUG */
      if (v > max_cycles) then
         if (master .and. norm_lhs/norm_rhs > norm_tol) call warn("[multigrid_gravity:vcycle_hg] Not enough V-cycles to achieve convergence.")
         v = max_cycles
      endif

      vstat%norm_final = norm_lhs/norm_rhs
      if (.not. verbose_vcycle) call vstat%brief_v_log

      call leaves%check_dirty(solution, "final_solution")

      call history%store_solution

   end subroutine vcycle_hg

!>
!! \brief Calculate the residuum for the Poisson equation.
!!
!! \details Different orders of the Laplace operator result in different quality of the approximation of potential
!! * Second order (O_I2) is the simplest operator (3-point in 1D, 5-point in 2D and 7-point in 3D). Its error is proportional to cg%dx**2
!! * Fourth order (O_I4) has width of 5 cells in each direction (5, 9 and 13-points in 1D, 2D and 3D, respectively). Its error is proportional to cg%dx**4.
!!   The value of its coefficients depends on interpretation of variables (e.g. point values, integral over cell)
!! * Fourth order Mehrstellen (-O_I4) is compact (takes 3, 9, and 27 solution values in 1D, 2D and 3D),
!!   but requires values of adjacent source cells as well (3, 5 and 7 points in 1D, 2D and 3D).
!!   Its error is proportional to cg%dx**4 and tends to be smaller than the error of simple fourth order operator by a factor of 3..4.
!!   \todo check how the coefficients depend on interpretation of cell values (center point values vs integral over cell)
!!
!! Note that each kind of Laplace operator requires its own approximate solver (relaxation scheme, Green function for FFT) for optimal convergence.
!! \warning Relaxation is implemented only for second order operator.
!<

   subroutine residual(cg_llst, src, soln, def)

      use cg_leaves,  only: cg_leaves_T
      use constants,  only: O_I2, O_I4
      use dataio_pub, only: die

      implicit none

      class(cg_leaves_T), intent(in) :: cg_llst !< pointer to a level for which we approximate the solution
      integer(kind=4),    intent(in) :: src     !< index of source in cg%q(:)
      integer(kind=4),    intent(in) :: soln    !< index of solution in cg%q(:)
      integer(kind=4),    intent(in) :: def     !< index of defect in cg%q(:)

      select case (ord_laplacian)
      case (O_I2)
         call residual2(cg_llst, src, soln, def)
      case (O_I4)
         call residual4(cg_llst, src, soln, def)
      case (-O_I4)
         call residual_Mehrstellen(cg_llst, src, soln, def)
      case default
         call die("[multigrid_gravity:residual] The parameter 'ord_laplacian' must be 2 or 4 or -4")
      end select

   end subroutine residual

!> \brief 2nd order Laplacian

   subroutine residual2(cg_llst, src, soln, def)

      use cg_list,       only: cg_list_element
      use cg_leaves,     only: cg_leaves_T
      use constants,     only: xdim, ydim, zdim, ndims, GEO_XYZ, GEO_RPZ, half, BND_NEGREF
      use dataio_pub,    only: die
      use domain,        only: dom
      use grid_cont,     only: grid_container
      use multigridvars, only: multidim_code_3D

      implicit none

      class(cg_leaves_T), intent(in) :: cg_llst !< pointer to a level for which we approximate the solution
      integer(kind=4),    intent(in) :: src     !< index of source in cg%q(:)
      integer(kind=4),    intent(in) :: soln    !< index of solution in cg%q(:)
      integer(kind=4),    intent(in) :: def     !< index of defect in cg%q(:)

      integer                         :: i, j, k
      real                            :: L0, Lx, Ly, Lz
      real, dimension(:), allocatable :: Lx1, Ly_a, L0_a
      type(cg_list_element), pointer  :: cgl
      type(grid_container),  pointer  :: cg

      call cg_llst%arr3d_boundaries(soln, bnd_type = BND_NEGREF)
      ! corners are required for non-cartesian decompositions because current implementation of arr3d_boundaries may use overlapping buffers at triple points

      ! Possible optimization candidate: reduce cache misses (secondary importance, cache-aware implementation required)
      ! Explicit loop over k gives here better performance than array operation due to less cache misses (at least on 32^3 and 64^3 arrays)
      cgl => cg_llst%first
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
!               Lx = cg%idx2 ! already set
!               Lz = cg%idz2 ! already set
               !> \todo convert Lx1, Ly_a and L0_a into precomputed arrays
               allocate(Lx1(size(cg%inv_x)), Ly_a(size(cg%inv_x)), L0_a(size(cg%inv_x)))
               Ly_a(cg%is:cg%ie) = cg%idy2 * cg%inv_x(cg%is:cg%ie)**2 ! cylindrical factor
               Lx1 (cg%is:cg%ie) = half * (cg%idx * cg%inv_x(cg%is:cg%ie))
               L0_a(cg%is:cg%ie) = -2. * (Lx + Ly_a(cg%is:cg%ie) + Lz)
               do k = cg%ks, cg%ke
                  do j = cg%js, cg%je
                     do i = cg%is, cg%ie
                        cg%q(def)%arr(i, j, k) = cg%q(src)%arr (i,   j,   k)   - cg%q(soln)%arr(i,   j,   k)    * L0_a(i)
                        if (dom%has_dir(xdim))   cg%q(def)%arr (i,   j,   k)   = cg%q(def) %arr(i,   j,   k)         - &
                             &                  (cg%q(soln)%arr(i+1, j,   k)   + cg%q(soln)%arr(i-1, j,   k))   * Lx - &
                             &                  (cg%q(soln)%arr(i+1, j,   k)   - cg%q(soln)%arr(i-1, j,   k))   * Lx1(i)    ! cylindrical term
                        if (dom%has_dir(ydim))   cg%q(def)%arr (i,   j,   k)   = cg%q(def) %arr(i,   j,   k)         - &
                             &                  (cg%q(soln)%arr(i,   j+1, k)   + cg%q(soln)%arr(i,   j-1, k))   * Ly_a(i)
                        if (dom%has_dir(zdim))   cg%q(def)%arr (i,   j,   k)   = cg%q(def) %arr(i,   j,   k)         - &
                             &                  (cg%q(soln)%arr(i,   j,   k+1) + cg%q(soln)%arr(i,   j,   k-1)) * Lz
                     enddo
                  enddo
               enddo
               deallocate(Lx1, Ly_a, L0_a)
            case default
               call die("[multigrid_gravity:residual2] Unsupported geometry.")
         end select
         cgl => cgl%nxt
      enddo

   end subroutine residual2

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

   subroutine residual4(cg_llst, src, soln, def)

      use cg_list,       only: cg_list_element
      use cg_leaves,     only: cg_leaves_T
      use constants,     only: ndims, idm2, xdim, ydim, zdim, BND_NEGREF
      use dataio_pub,    only: die, warn
      use domain,        only: dom
      use grid_cont,     only: grid_container
      use mpisetup,      only: master
      use multigridvars, only: grav_bnd, bnd_givenval
      use named_array,   only: p3

      implicit none

      class(cg_leaves_T), intent(in) :: cg_llst !< pointer to a level for which we approximate the solution
      integer(kind=4),    intent(in) :: src     !< index of source in cg%q(:)
      integer(kind=4),    intent(in) :: soln    !< index of solution in cg%q(:)
      integer(kind=4),    intent(in) :: def     !< index of defect in cg%q(:)

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

      call cg_llst%arr3d_boundaries(soln, bnd_type = BND_NEGREF)

      c21 = 1.
      c42 = - L4_scaling * L4_strength
      c41 = c21 + 4. * L4_scaling * L4_strength
      !c20 = -2.
      !c40 = c20 - 6. * L4_strength

      cgl => cg_llst%first
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

!>
!! \brief Compute the residual using compact, 4th order Mehrstellen formula
!!
!! \details 2nd order Laplace operator contains some 4th order terms on the potential size. To improve accuracy we may add a little laplacian of the source on the RHS
!! to match these 4th order terms. Careful calculation leads to such numerical operators:
!!
!!                     |   1   |                    | 1   4 1 |
!! 2D: residuum = 1/12 | 1 8 1 | source - 1/6 h**-2 | 4 -20 4 | solution
!!                     |   1   |                    | 1   4 1 |
!!
!!                     |       | |   1   | |       |                    |   1   | | 1   2 1 | |   1   |
!! 3D: residuum = 1/12 |   1   | | 1 6 1 | |   1   | source - 1/6 h**-2 | 1 2 1 | | 2 -24 2 | | 1 2 1 | solution
!!                     |       | |   1   | |       |                    |   1   | | 1   2 1 | |   1   |
!!
!! where h = cg%dx = cg%dy = cg%dz
!! On an unequally spaced grid (cg%dx /= cg%dy /= cg%dz), the stencil convolved with the solution is a bit more complicated:
!!
!!                 |  1  -2  1 |               |  1  10  1 |
!! 2D: 1/12 dx**-2 | 10 -20 10 | + 1/12 dy**-2 | -2 -20 -2 |
!!                 |  1  -2  1 |               |  1  10  1 |
!!
!!                 |        | | 1  -2 1 | |        |               |   1   | |  1   8  1 | |   1   |               |   1   | |     -2    | |   1   |
!! 3D: 1/12 dx**-2 | 1 -2 1 | | 8 -16 8 | | 1 -2 1 | + 1/12 dy**-2 |  -2   | | -2 -16 -2 | |  -2   | + 1/12 dz**-2 | 1 8 1 | | -2 -16 -2 | | 1 8 1 |
!!                 |        | | 1  -2 1 | |        |               |   1   | |  1   8  1 | |   1   |               |   1   | |     -2    | |   1   |
!!
!! For even higher order of accuracy look for HODIE schemes which are more general then the one described above.
!!
!! \warning The implementation is somehow incomplete - different coefficients are required to do efficient RBGS relaxation
!!
!! \warning This implementation does not support cylindrical coordinates
!<

   subroutine residual_Mehrstellen(cg_llst, src, soln, def)

      use cg_list,       only: cg_list_element
      use cg_leaves,     only: cg_leaves_T
      use constants,     only: ndims, xdim, ydim, zdim, BND_NEGREF, LO, HI
      use dataio_pub,    only: warn
      use domain,        only: dom
      use grid_cont,     only: grid_container
      use mpisetup,      only: master
      use multigridvars, only: grav_bnd, bnd_givenval
      use named_array,   only: p3

      implicit none

      class(cg_leaves_T), intent(in) :: cg_llst !< pointer to a level for which we approximate the solution
      integer(kind=4),    intent(in) :: src     !< index of source in cg%q(:)
      integer(kind=4),    intent(in) :: soln    !< index of solution in cg%q(:)
      integer(kind=4),    intent(in) :: def     !< index of defect in cg%q(:)

      real                :: L0, Lx, Ly, Lz, Lxy, Lxz, Lyz
      integer, parameter  :: L2w = 2             ! #layers of boundary cells for L2 operator

      logical, save       :: firstcall = .true.
      integer             :: i, j, k
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg
      integer(kind=4), dimension(ndims,ndims,LO:HI) :: idm

      if (firstcall) then
         if (master) call warn("[multigrid_gravity:residual_Mehrstellen] residual order 4 is experimental.")
         firstcall = .false.
      endif

      call cg_llst%arr3d_boundaries(soln, bnd_type = BND_NEGREF)
      call cg_llst%arr3d_boundaries(src,  bnd_type = BND_NEGREF)

      idm = 0
      do i = xdim, zdim
         if (dom%has_dir(i)) idm(i, i, :) = 1
      enddo

      cgl => cg_llst%first
      do while (associated(cgl))
         cg => cgl%cg

         Lxy = 0. ; if (dom%has_dir(xdim) .and. dom%has_dir(ydim)) Lxy = (cg%idx2 + cg%idy2) / 12.
         Lxz = 0. ; if (dom%has_dir(xdim) .and. dom%has_dir(zdim)) Lxz = (cg%idx2 + cg%idz2) / 12.
         Lyz = 0. ; if (dom%has_dir(ydim) .and. dom%has_dir(zdim)) Lyz = (cg%idy2 + cg%idz2) / 12.
         Lx  = 0. ; if (dom%has_dir(xdim)) Lx = cg%idx2 - 2. * (Lxy + Lxz)
         Ly  = 0. ; if (dom%has_dir(ydim)) Ly = cg%idy2 - 2. * (Lxy + Lyz)
         Lz  = 0. ; if (dom%has_dir(zdim)) Lz = cg%idz2 - 2. * (Lxz + Lyz)
         L0  = -2. * (Lx + Ly + Lz) - 4. * (Lxy + Lxz + Lyz)

         p3 => cg%q(def)%span(cg%ijkse)

         p3 = (12.-2*dom%eff_dim)/12.*cg%q(src)%span(cg%ijkse) - L0*cg%q(soln)%span(cg%ijkse)
         if (dom%has_dir(xdim)) p3 = p3 + &
              (cg%q(src)%span(cg%ijkse-idm(xdim,:,:)) + cg%q(src)%span(cg%ijkse+idm(xdim,:,:)))/12. &
              - Lx * (cg%q(soln)%span(cg%ijkse-idm(xdim,:,:)) + cg%q(soln)%span(cg%ijkse+idm(xdim,:,:)))
         if (dom%has_dir(ydim)) p3 = p3 + &
              (cg%q(src)%span(cg%ijkse-idm(ydim,:,:)) + cg%q(src)%span(cg%ijkse+idm(ydim,:,:)))/12. &
              - Ly * (cg%q(soln)%span(cg%ijkse-idm(ydim,:,:)) + cg%q(soln)%span(cg%ijkse+idm(ydim,:,:))) &
              - Lxy * (cg%q(soln)%span(cg%ijkse-idm(xdim,:,:)-idm(ydim,:,:)) + cg%q(soln)%span(cg%ijkse-idm(xdim,:,:)+idm(ydim,:,:)) + &
              &        cg%q(soln)%span(cg%ijkse+idm(xdim,:,:)-idm(ydim,:,:)) + cg%q(soln)%span(cg%ijkse+idm(xdim,:,:)+idm(ydim,:,:)) )
         if (dom%has_dir(zdim)) p3 = p3 + &
              (cg%q(src)%span(cg%ijkse-idm(zdim,:,:)) + cg%q(src)%span(cg%ijkse+idm(zdim,:,:)))/12. &
              - Lz * (cg%q(soln)%span(cg%ijkse-idm(zdim,:,:)) + cg%q(soln)%span(cg%ijkse+idm(zdim,:,:))) &
              - Lxz * (cg%q(soln)%span(cg%ijkse-idm(xdim,:,:)-idm(zdim,:,:)) + cg%q(soln)%span(cg%ijkse-idm(xdim,:,:)+idm(zdim,:,:)) + &
              &        cg%q(soln)%span(cg%ijkse+idm(xdim,:,:)-idm(zdim,:,:)) + cg%q(soln)%span(cg%ijkse+idm(xdim,:,:)+idm(zdim,:,:)) ) &
              - Lyz * (cg%q(soln)%span(cg%ijkse-idm(zdim,:,:)-idm(ydim,:,:)) + cg%q(soln)%span(cg%ijkse-idm(zdim,:,:)+idm(ydim,:,:)) + &
              &        cg%q(soln)%span(cg%ijkse+idm(zdim,:,:)-idm(ydim,:,:)) + cg%q(soln)%span(cg%ijkse+idm(zdim,:,:)+idm(ydim,:,:)) )
         ! OPT: Perhaps in 3D this would be more efficient without ifs.
         ! OPT: Try direct operations on cg%q(soln)%arr without calling span
         ! 2D and 1D rns do not need too much optimization as they will be rarely used in production runs with selfgravity

         ! WARNING: not optimized
         if (grav_bnd == bnd_givenval) then ! probably also in some other cases
            ! Use L2 Laplacian in two layers of cells next to the boundary because L4 seems to be incompatible with present image mass construction
            !> \todo try to fix this
            Lx  = 0. ; if (dom%has_dir(xdim)) Lx = cg%idx2
            Ly  = 0. ; if (dom%has_dir(ydim)) Ly = cg%idy2
            Lz  = 0. ; if (dom%has_dir(zdim)) Lz = cg%idz2
            L0 = -2. * (Lx + Ly + Lz)

            do k = cg%ks, cg%ke
               do j = cg%js, cg%je
                  do i = cg%is, cg%ie
                     if ( i<cg%is+L2w .or. i>cg%ie-L2w .or. j<cg%js+L2w .or. j>cg%je-L2w .or. k<cg%ks+L2w .or. k>cg%ke-L2w) then
                        cg%q(def)%arr        (i,         j,         k)         = cg%q(src)%arr (i,         j,         k      )         - &
                             ( cg%q(soln)%arr(i-dom%D_x, j,         k)         + cg%q(soln)%arr(i+dom%D_x, j,         k      ))   * Lx - &
                             ( cg%q(soln)%arr(i,         j-dom%D_x, k)         + cg%q(soln)%arr(i,         j+dom%D_x, k      ))   * Ly - &
                             ( cg%q(soln)%arr(i,         j,         k-dom%D_x) + cg%q(soln)%arr(i,         j,         k+dom%D_x)) * Lz - &
                             & cg%q(soln)%arr(i,         j,         k)                                                            * L0
                     endif
                  enddo
               enddo
            enddo
         endif
         cgl => cgl%nxt
      enddo

   end subroutine residual_Mehrstellen

!!$ ============================================================================
!>
!! \brief This routine has to find an approximate solution for given source field and implemented differential operator
!<

   subroutine approximate_solution(curl, src, soln)

      use cg_level_connected,  only: cg_level_connected_T
!!$      use constants,           only: fft_none
!!$      use multigrid_fftapprox, only: approximate_solution_fft

      implicit none

      type(cg_level_connected_T), pointer, intent(in) :: curl !< pointer to a level for which we approximate the solution
      integer(kind=4),                     intent(in) :: src  !< index of source in cg%q(:)
      integer(kind=4),                     intent(in) :: soln !< index of solution in cg%q(:)

      call curl%check_dirty(src, "approx_soln src-")

!!$      if (curl%fft_type /= fft_none) then
!!$         call approximate_solution_fft(curl, src, soln)
!!$      else
         call approximate_solution_rbgs(curl, src, soln)
!!$      endif

      call curl%check_dirty(soln, "approx_soln soln+")

   end subroutine approximate_solution

!>
!! \brief Red-Black Gauss-Seidel relaxation.
!!
!! \details This is the most costly routine in a serial run. Try to find optimal values for nsmool and nsmoob.
!! This routine also depends a lot on communication so it  may limit scalability of the multigrid.
!!
!! This implementation is suitable for 2nd order Laplace operator (routine residual2).
!! 4th order operators (residual4 and residual_Mehrstellen) require different relaxation operators or the convergence rate will be severily limited.
!<

   subroutine approximate_solution_rbgs(curl, src, soln)

      use cg_level_coarsest,  only: coarsest
      use cg_level_connected, only: cg_level_connected_T
      use cg_list,            only: cg_list_element
      use cg_list_dataop,     only: dirty_label
      use constants,          only: xdim, ydim, zdim, ndims, GEO_XYZ, GEO_RPZ, I_ONE, BND_NEGREF, O_I4
      use dataio_pub,         only: die
      use domain,             only: dom
      use global,             only: dirty_debug
      use grid_cont,          only: grid_container
      use multigridvars,      only: correction, multidim_code_3D, nsmool

      implicit none

      type(cg_level_connected_T), pointer, intent(in) :: curl !< pointer to a level for which we approximate the solution
      integer(kind=4),                     intent(in) :: src  !< index of source in cg%q(:)
      integer(kind=4),                     intent(in) :: soln !< index of solution in cg%q(:)

      integer, parameter :: RED_BLACK = 2 !< the checkerboard requires two sweeps

      integer :: n, i, j, k, i1, j1, k1, id, jd, kd
      integer :: nsmoo
      real    :: crx, crx1, cry, crz, cr, L0, Lx, Ly, Lz, Lxy, Lxz, Lyz
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg

      if (associated(curl, coarsest%level)) then
         nsmoo = nsmoob
         !> \todo Implement automatic convergence check on coarsest level (not very important when we have a FFT solver for coarsest level)
      else
         nsmoo = nsmool
         if (soln == correction) call curl%coarser%prolong_q_1var(soln) ! make sure that prolongation is called only in ascending (coarse -> fine) part of V-cycle.
         !> \warning this may be incompatible with V-cycles other than Huang - Greengard
      endif

      if (dom%geometry_type == GEO_RPZ .and. .not. multidim_code_3D) call die("[multigrid_gravity:approximate_solution_rbgs] multidim_code_3D = .false. not implemented")

!      if (ord_laplacian == -O_I4) call curl%arr3d_boundaries(src, bnd_type = BND_NEGREF) ! required when we use 7-point source term, not just 1-point
      ! Also required when we want to eliminate some communication of soln at a cost of expanding relaxated area into guardcells

      do n = 1, RED_BLACK*nsmoo
         call curl%arr3d_boundaries(soln, bnd_type = BND_NEGREF)

         if (dirty_debug) then
            write(dirty_label, '(a,i5)')"relax soln- smoo=", n
            call curl%check_dirty(soln, dirty_label, expand=I_ONE)
         endif
         cgl => curl%first
         do while (associated(cgl))
            cg => cgl%cg

            ! Possible optimization: this is the most costly part of the RBGS relaxation (instruction count, read and write data, L1 and L2 read cache miss)
            ! do n = 1, nsmoo
            !    call curl%arr3d_boundaries(soln, bnd_type = BND_NEGREF)
            !    relax single layer of red cells at all faces
            !    call curl%arr3d_boundaries(soln, bnd_type = BND_NEGREF)
            !    relax interior cells (except for single layer of cells at all faces), first red, then 1-cell behind black one.
            !    relax single layer of black cells at all faces
            ! enddo

            ! with explicit outer loops it is easier to describe a 3-D checkerboard :-)

            if (ord_laplacian == -O_I4) then
               if (dom%geometry_type == GEO_RPZ) call die("[multigrid_gravity:approximate_solution_rbgs] Relaxation for Mehrstellen not implemented for noncartesian grid")
               if (dom%eff_dim==ndims) then
                  Lxy = 0. ; if (dom%has_dir(xdim) .and. dom%has_dir(ydim)) Lxy = (cg%idx2 + cg%idy2) / 12.
                  Lxz = 0. ; if (dom%has_dir(xdim) .and. dom%has_dir(zdim)) Lxz = (cg%idx2 + cg%idz2) / 12.
                  Lyz = 0. ; if (dom%has_dir(ydim) .and. dom%has_dir(zdim)) Lyz = (cg%idy2 + cg%idz2) / 12.
                  Lx  = 0. ; if (dom%has_dir(xdim)) Lx = cg%idx2 - 2. * (Lxy + Lxz)
                  Ly  = 0. ; if (dom%has_dir(ydim)) Ly = cg%idy2 - 2. * (Lxy + Lyz)
                  Lz  = 0. ; if (dom%has_dir(zdim)) Lz = cg%idz2 - 2. * (Lxz + Lyz)
                  L0  = 2. * (Lx + Ly + Lz) + 4. * (Lxy + Lxz + Lyz)
                  do k = cg%ks, cg%ke
                     do j = cg%js, cg%je
                        i1 = cg%is + int(mod(n+cg%is+j+k, int(RED_BLACK)), kind=4)
                        cg%q(soln)%arr(i1  :cg%ie  :2, j,   k) = &
                             &   (cg%q(soln)%arr(i1-1:cg%ie-1:2, j,   k  ) + cg%q(soln)%arr(i1+1:cg%ie+1:2, j,   k  ))*Lx/L0   &
                             +   (cg%q(soln)%arr(i1  :cg%ie  :2, j-1, k  ) + cg%q(soln)%arr(i1  :cg%ie  :2, j+1, k  ))*Ly/L0   &
                             +   (cg%q(soln)%arr(i1  :cg%ie  :2, j,   k-1) + cg%q(soln)%arr(i1  :cg%ie  :2, j,   k+1))*Lz/L0   &
                             +  ((cg%q(soln)%arr(i1-1:cg%ie-1:2, j-1, k  ) + cg%q(soln)%arr(i1+1:cg%ie+1:2, j-1, k  ))         &
                             +   (cg%q(soln)%arr(i1-1:cg%ie-1:2, j+1, k  ) + cg%q(soln)%arr(i1+1:cg%ie+1:2, j+1, k  )))*Lxy/L0 &
                             +  ((cg%q(soln)%arr(i1  :cg%ie  :2, j-1, k-1) + cg%q(soln)%arr(i1  :cg%ie  :2, j+1, k-1))         &
                             +   (cg%q(soln)%arr(i1  :cg%ie  :2, j-1, k+1) + cg%q(soln)%arr(i1  :cg%ie  :2, j+1, k+1)))*Lyz/L0 &
                             +  ((cg%q(soln)%arr(i1-1:cg%ie-1:2, j,   k-1) + cg%q(soln)%arr(i1+1:cg%ie+1:2, j,   k-1))         &
                             +   (cg%q(soln)%arr(i1-1:cg%ie-1:2, j,   k+1) + cg%q(soln)%arr(i1+1:cg%ie+1:2, j,   k+1)))*Lxz/L0 &
                             -    cg%q(src)%arr (i1  :cg%ie  :2, j,   k  ) / L0
! For some weird reasons the formula that comes directly from the Mehrstellen operator worsens convergence
! Note that it requires enabling call curl%arr3d_boundaries(src, ...) above
!!$                             - ((12-2*dom%eff_dim)*cg%q(src)%arr (i1  :cg%ie  :2, j,   k  )                                                     &
!!$                             &  + cg%q(src)%arr (i1-1:cg%ie-1:2, j,   k  ) + cg%q(src)%arr (i1+1:cg%ie+1:2, j,   k  )          &
!!$                             &  + cg%q(src)%arr (i1  :cg%ie  :2, j-1, k  ) + cg%q(src)%arr (i1  :cg%ie  :2, j+1, k  )          &
!!$                             &  + cg%q(src)%arr (i1  :cg%ie  :2, j,   k-1) + cg%q(src)%arr (i1  :cg%ie  :2, j,   k+1))/(12. * L0)
                     enddo
                  enddo
               else
                  call die("[multigrid_gravity:approximate_solution_rbgs] Relaxation for Mehrstellen implemented only in 3D")
               endif
            else
               if (dom%eff_dim==ndims .and. .not. multidim_code_3D) then
                  do k = cg%ks, cg%ke
                     do j = cg%js, cg%je
                        i1 = cg%is + int(mod(n+cg%is+j+k, int(RED_BLACK)), kind=4)
                        if (dom%geometry_type == GEO_RPZ) then
!!$                  cg%q(soln)%arr(i1  :cg%ie  :2, j,   k) = &
!!$                       cg%mg%rx * (cg%q(soln)%arr(i1-1:cg%ie-1:2, j,   k  ) + cg%q(soln)%arr(i1+1:cg%ie+1:2, j,   k))   + &
!!$                       cg%mg%ry * (cg%q(soln)%arr(i1  :cg%ie  :2, j-1, k  ) + cg%q(soln)%arr(i1  :cg%ie  :2, j+1, k))   + &
!!$                       cg%mg%rz * (cg%q(soln)%arr(i1  :cg%ie  :2, j,   k-1) + cg%q(soln)%arr(i1  :cg%ie  :2, j,   k+1)) - &
!!$                       cg%mg%r  *  cg%q(src)%arr( i1  :cg%ie  :2, j,   k  )  + &
!!$                       cg%mg%rx * (cg%q(soln)%arr(i1+1:cg%ie+1:2, j,   k  ) - cg%q(soln)%arr(i1-1:cg%ie-1:2, j,   k)) * fac(i1:cg%ie:2)
                           call die("[multigrid_gravity:approximate_solution_rbgs] This variant of relaxation loop is not implemented for cylindrical coordinates.")
                        else
                           cg%q(soln)%arr(i1  :cg%ie  :2, j,   k) = &
                                cg%mg%rx * (cg%q(soln)%arr(i1-1:cg%ie-1:2, j,   k)   + cg%q(soln)%arr(i1+1:cg%ie+1:2, j,   k))   + &
                                cg%mg%ry * (cg%q(soln)%arr(i1  :cg%ie  :2, j-1, k)   + cg%q(soln)%arr(i1  :cg%ie  :2, j+1, k))   + &
                                cg%mg%rz * (cg%q(soln)%arr(i1  :cg%ie  :2, j,   k-1) + cg%q(soln)%arr(i1  :cg%ie  :2, j,   k+1)) - &
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

                  if (kd == RED_BLACK) k1 = cg%ks + int(mod(n+cg%ks, int(RED_BLACK)), kind=4)
                  select case (dom%geometry_type)
                     case (GEO_XYZ)
                        do k = k1, cg%ke, kd
                           if (jd == RED_BLACK) j1 = cg%js + int(mod(n+cg%js+k, int(RED_BLACK)), kind=4)
                           do j = j1, cg%je, jd
                              if (id == RED_BLACK) i1 = cg%is + int(mod(n+cg%is+j+k, int(RED_BLACK)), kind=4)
                              cg%q(soln)%arr                           (i1  :cg%ie  :id, j,   k)   = &
                                   & (1. - Jacobi_damp)* cg%q(soln)%arr(i1  :cg%ie  :id, j,   k)   - &
                                   &       Jacobi_damp * cg%q(src)%arr (i1  :cg%ie  :id, j,   k)   * cg%mg%r
                              if (dom%has_dir(xdim))     cg%q(soln)%arr(i1  :cg%ie  :id, j,   k)   = cg%q(soln)%arr(i1  :cg%ie  :id, j,   k)    + &
                                   &       Jacobi_damp *(cg%q(soln)%arr(i1-1:cg%ie-1:id, j,   k)   + cg%q(soln)%arr(i1+1:cg%ie+1:id, j,   k))   * cg%mg%rx
                              if (dom%has_dir(ydim))     cg%q(soln)%arr(i1  :cg%ie  :id, j,   k)   = cg%q(soln)%arr(i1  :cg%ie  :id, j,   k)    + &
                                   &       Jacobi_damp *(cg%q(soln)%arr(i1  :cg%ie  :id, j-1, k)   + cg%q(soln)%arr(i1  :cg%ie  :id, j+1, k))   * cg%mg%ry
                              if (dom%has_dir(zdim))     cg%q(soln)%arr(i1  :cg%ie  :id, j,   k)   = cg%q(soln)%arr(i1  :cg%ie  :id, j,   k)    + &
                                   &       Jacobi_damp *(cg%q(soln)%arr(i1  :cg%ie  :id, j,   k-1) + cg%q(soln)%arr(i1  :cg%ie  :id, j,   k+1)) * cg%mg%rz
                           enddo
                        enddo
                     case (GEO_RPZ)
                        do k = k1, cg%ke, kd
                           if (jd == RED_BLACK) j1 = cg%js + int(mod(n+cg%js+k, int(RED_BLACK)), kind=4)
                           do j = j1, cg%je, jd
                              if (id == RED_BLACK) i1 = cg%is + int(mod(n+cg%is+j+k, int(RED_BLACK)), kind=4)
                              do i = i1, cg%ie, id
                                 cr  = overrelax / 2.
                                 crx = cg%dvol2 * cg%idx2 * cg%x(i)**2
                                 cry = cg%dvol2 * cg%idy2
                                 crz = cg%dvol2 * cg%idz2 * cg%x(i)**2
                                 cr  = cr / (crx + cry + crz)
                                 crx = overrelax_xyz(xdim)* crx * cr
                                 cry = overrelax_xyz(ydim)* cry * cr
                                 crz = overrelax_xyz(zdim)* crz * cr
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
            endif
            cgl => cgl%nxt
         enddo

         if (dirty_debug) then
            write(dirty_label, '(a,i5)')"relax soln+ smoo=", n
            call curl%check_dirty(soln, dirty_label)
         endif

      enddo

   end subroutine approximate_solution_rbgs

!> \brief Solve finest level if allowed (single cg and single thread)

   subroutine fft_solve_roof

      use cg_level_finest,     only: finest
      use constants,           only: fft_none
      use dataio_pub,          only: die
      use grid_cont,           only: grid_container
      use multigrid_fftapprox, only: fft_convolve
      use multigridvars,       only: source, solution
      use named_array,         only: p3

      implicit none

      type(grid_container), pointer :: cg

      if (associated(finest%level%first)) then
         if (associated(finest%level%first%nxt)) call die("[multigrid_gravity:fft_solve_roof] multicg not possible")
      endif

      if (finest%level%fft_type == fft_none) return

      cg => finest%level%first%cg
      p3 => cg%q(source)%span(cg%ijkse)
      cg%mg%src(:, :, :) = p3
      call fft_convolve(finest%level)
      p3 => cg%q(solution)%span(cg%ijkse)
      p3 = cg%mg%src(:, :, :)

   end subroutine fft_solve_roof

end module multigrid_gravity
