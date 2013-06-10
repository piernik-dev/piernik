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

   use constants,          only: cbuff_len, dsetnamelen
   use multigrid_vstats,   only: vcycle_stats
   use multigrid_old_soln, only: soln_history

   implicit none

   private
   public :: multigrid_grav_par, init_multigrid_grav, cleanup_multigrid_grav, multigrid_solve_grav, init_multigrid_grav_ext, ordL, invalidate_history

   include "fftw3.f"
   ! constants from fftw3.f
   !   integer, parameter :: FFTW_MEASURE=0, FFTW_PATIENT=32, FFTW_ESTIMATE=64
   !   integer, parameter :: FFTW_RODFT01=8, FFTW_RODFT10=9

   ! namelist parameters
   real               :: norm_tol                                     !< stop V-cycle iterations when the ratio of norms ||residual||/||source|| is below this value
   real               :: vcycle_abort                                 !< abort the V-cycle when lhs norm raises by this factor
   real               :: vcycle_giveup                                !< exit the V-cycle when convergence ratio drops below that level
   integer(kind=4)    :: max_cycles                                   !< Maximum allowed number of V-cycles
   integer(kind=4)    :: nsmoob                                       !< smoothing cycles on coarsest level when cannot use FFT. (a convergence check would be much better)
   integer(kind=4)    :: ord_laplacian                                !< Laplace operator order; allowed values are 2, -4 (default) and 4 (not fully implemented)
   integer(kind=4)    :: ord_laplacian_outer                          !< Laplace operator order for isolated boundaries (useful as long as -4 is not fully implemented)
   logical            :: trust_fft_solution                           !< Bypass the V-cycle, when doing FFT on whole domain, make sure first that FFT is properly set up.
   logical            :: base_no_fft                                  !< Deny solving the coarsest level with FFT. Can be very slow.
   logical            :: prefer_rbgs_relaxation                       !< Prefer relaxation over FFT local solver. Typically faster.
   !> \todo allow to perform one or more V-cycles with FFT method, the switch to the RBGS (may save one V-cycle in some cases)
   logical            :: fft_patient                                  !< Spend more time in init_multigrid to find faster fft plan
   character(len=cbuff_len) :: grav_bnd_str                           !< Type of gravitational boundary conditions.
   logical            :: require_FFT                                  !< .true. if we use FFT solver anywhere (and need face prolongation)

   ! Conjugate gradients
   logical            :: use_CG                                       !< .true. if we want to use multigrid-preconditioned conjugate gradient iterations
   logical            :: use_CG_outer                                 !< .true. if we want to use multigrid-preconditioned conjugate gradient iterations for outer potential
   logical            :: use_MGpreconditioning                        !< .true. if we want to use multigrid preconditioner
   character(len=dsetnamelen), parameter :: cg_corr_n = "cg_correction" !< correction vector for CG
   integer(kind=4)    :: cg_corr                                      !< index of the cg-correction vector

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
!! <tr><td>vcycle_giveup         </td><td>1.5    </td><td>real value     </td><td>\copydoc multigrid_gravity::vcycle_giveup         </td></tr>
!! <tr><td>max_cycles            </td><td>20     </td><td>integer value  </td><td>\copydoc multigrid_gravity::max_cycles            </td></tr>
!! <tr><td>nsmool                </td><td>dom%nb </td><td>integer value  </td><td>\copydoc multigridvars::nsmool                    </td></tr>
!! <tr><td>nsmoob                </td><td>100    </td><td>integer value  </td><td>\copydoc multigrid_gravity::nsmoob                </td></tr>
!! <tr><td>overrelax             </td><td>1.     </td><td>real value     </td><td>\copydoc multigrid_gravity::overrelax             </td></tr>
!! <tr><td>overrelax_xxyz(ndims) </td><td>1.     </td><td>real value     </td><td>\copydoc multigrid_gravity::overrelax_xyz         </td></tr>
!! <tr><td>Jacobi_damp           </td><td>1.     </td><td>real value     </td><td>\copydoc multigrid_gravity::jacobi_damp           </td></tr>
!! <tr><td>L4_strength           </td><td>1.0    </td><td>real value     </td><td>\copydoc multigrid_gravity::l4_strength           </td></tr>
!! <tr><td>nsmoof                </td><td>1      </td><td>integer value  </td><td>\copydoc multigridvars::nsmoof                    </td></tr>
!! <tr><td>ord_laplacian         </td><td>-4     </td><td>integer value  </td><td>\copydoc multigrid_gravity::ord_laplacian         </td></tr>
!! <tr><td>ord_laplacian_outer   </td><td>2      </td><td>integer value  </td><td>\copydoc multigrid_gravity::ord_laplacian_outer   </td></tr>
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
!! <tr><td>use_CG                </td><td>.false.</td><td>logical        </td><td>\copydoc multigrid_gravity::use_CG                </td></tr>
!! <tr><td>use_CG_outer          </td><td>.false.</td><td>logical        </td><td>\copydoc multigrid_gravity::use_CG_outer          </td></tr>
!! <tr><td>use_MGpreconditioning </td><td>.true. </td><td>logical        </td><td>\copydoc multigrid_gravity::use_MGpreconditioning </td></tr>
!! <tr><td>grav_bnd_str          </td><td>"periodic"/"dirichlet"</td><td>string of chars</td><td>\copydoc multigrid_gravity::grav_bnd_str          </td></tr>
!! </table>
!! The list is active while \b "GRAV" and \b "MULTIGRID" are defined.
!! \n \n
!<
   subroutine multigrid_grav_par

      use cg_list_global,     only: all_cg
      use constants,          only: GEO_XYZ, GEO_RPZ, BND_PER, O_LIN, O_D2, O_I2, O_D4, I_ONE, INVALID
      use dataio_pub,         only: nh  ! QA_WARN required for diff_nml
      use dataio_pub,         only: msg, die, warn
      use domain,             only: dom, is_multicg !, is_uneven
      use mpisetup,           only: master, slave, ibuff, cbuff, rbuff, lbuff, nproc, piernik_MPI_Bcast
      use multigridvars,      only: single_base, bnd_invalid, bnd_isolated, bnd_periodic, bnd_dirichlet, grav_bnd, fft_full_relax, multidim_code_3D, nsmool, nsmoof, &
           &                        overrelax, overrelax_xyz, Jacobi_damp, L4_strength
      use multigrid_old_soln, only: nold_max, ord_time_extrap
      use multipole,          only: use_point_monopole, lmax, mmax, ord_prolong_mpole, coarsen_multipole, interp_pt2mom, interp_mom2pot
      use named_array_list,   only: qna

      implicit none

      integer       :: periodic_bnd_cnt   !< counter of periodic boundaries in existing directions
      logical, save :: frun = .true.      !< First run flag

      namelist /MULTIGRID_GRAVITY/ norm_tol, vcycle_abort, vcycle_giveup, max_cycles, nsmool, nsmoob, use_CG, use_CG_outer, use_MGpreconditioning, &
           &                       overrelax, overrelax_xyz, Jacobi_damp, L4_strength, nsmoof, ord_laplacian, ord_laplacian_outer, ord_time_extrap, &
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
      vcycle_giveup          = 1.5
      L4_strength            = 1.0

      coarsen_multipole      = 0
!      if (is_uneven) coarsen_multipole = 0
      lmax                   = 16
      mmax                   = -1 ! will be automatically set to lmax unless explicitly limited in problem.par
      max_cycles             = 20
      nsmool                 = -1  ! best to set it to dom%nb or its multiply
      nsmoob                 = 100
      nsmoof                 = 1
      select case (dom%geometry_type)
         case (GEO_XYZ)
            ord_laplacian    = O_D4
         case (GEO_RPZ)
            ord_laplacian    = O_I2
         case default
            ord_laplacian    = INVALID
      end select
      ord_laplacian_outer    = ord_laplacian
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
      use_CG                 = .false.
      use_CG_outer           = .false.
      use_MGpreconditioning  = .true.

      periodic_bnd_cnt = count(dom%periodic(:) .and. dom%has_dir(:))

      if (periodic_bnd_cnt == dom%eff_dim) then
         grav_bnd_str = "periodic"
      else
         grav_bnd_str = "dirichlet"
      endif

      if (master) then

         diff_nml(MULTIGRID_GRAVITY)

         if (nsmool < 0) nsmool = -nsmool * dom%nb

         ! FIXME when ready
         select case (dom%geometry_type)
            case (GEO_XYZ) ! do nothing
            case (GEO_RPZ)
               ! switch off FFT-related bits
               base_no_fft = .true.
               prefer_rbgs_relaxation = .true.
               if (any([ ord_laplacian, ord_laplacian_outer ] /= O_I2) .and. master) call warn("[multigrid_gravity:multigrid_grav_par] Laplacian order forced to 2]")
               ord_laplacian = O_I2
               ord_laplacian_outer = ord_laplacian
               L4_strength = 0.
               ! ord_prolong_mpole = O_INJ
            case default
               call die("[multigrid_gravity:multigrid_grav_par] Unsupported geometry.")
         end select

         if (is_multicg .and. .not. base_no_fft) then
            call warn("[multigrid_gravity:multigrid_grav_par] base_no_fft set to .true. for multicg configuration")
            base_no_fft = .true.
         endif

         if (.not. prefer_rbgs_relaxation .and. nproc > 1) then
            prefer_rbgs_relaxation = .true.
            call warn("[multigrid_gravity:multigrid_grav_par] FFT local solver disabled for multithreaded runs due to unresolved incompatibilities with mew grid features")
         endif

         if (ord_laplacian_outer /= ord_laplacian) call warn("[multigrid_gravity:multigrid_grav_par] ord_laplacian_outer /= ord_laplacian")

         rbuff(1) = norm_tol
         rbuff(2) = overrelax
         rbuff(3:5) = overrelax_xyz
         rbuff(6) = Jacobi_damp
         rbuff(7) = vcycle_abort
         rbuff(8) = vcycle_giveup
         rbuff(9) = L4_strength

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
         ibuff(11) = ord_laplacian_outer

         lbuff(1)  = use_point_monopole
         lbuff(2)  = trust_fft_solution
         lbuff(3)  = base_no_fft
         lbuff(4)  = prefer_rbgs_relaxation
         lbuff(5)  = fft_full_relax
         lbuff(6)  = fft_patient
         lbuff(7)  = interp_pt2mom
         lbuff(8)  = interp_mom2pot
         lbuff(9)  = multidim_code_3D
         lbuff(10) = use_CG
         lbuff(11) = use_CG_outer
         lbuff(12) = use_MGpreconditioning

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
         vcycle_giveup  = rbuff(8)
         L4_strength    = rbuff(9)

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
         ord_laplacian_outer = ibuff(11)

         use_point_monopole      = lbuff(1)
         trust_fft_solution      = lbuff(2)
         base_no_fft             = lbuff(3)
         prefer_rbgs_relaxation  = lbuff(4)
         fft_full_relax          = lbuff(5)
         fft_patient             = lbuff(6)
         interp_pt2mom           = lbuff(7)
         interp_mom2pot          = lbuff(8)
         multidim_code_3D        = lbuff(9)
         use_CG                  = lbuff(10)
         use_CG_outer            = lbuff(11)
         use_MGpreconditioning   = lbuff(12)

         grav_bnd_str   = cbuff(1)(1:len(grav_bnd_str))

      endif

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

      if (use_CG) then
         call all_cg%reg_var(cg_corr_n)
         cg_corr   = qna%ind(cg_corr_n)
      endif
      if (use_CG .or. use_CG_outer) then
         if (master) call warn("[multigrid_gravity:multigrid_grav_par] Multigrid-preconditioned conjugate gradient solver is experimental!")
      endif

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
      logical, save :: firstcall = .true.

      if (coarsen_multipole /= 0) then
         coarsen_multipole = 0
         if (master) call warn("[multigrid_gravity:init_multigrid_grav] multipole coarsening temporarily disabled")
      endif

      if (firstcall) call leaves%set_q_value(qna%ind(sgp_n), 0.) !Initialize all the guardcells, even those which does not impact the solution

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

      if (grav_bnd == bnd_isolated .and. firstcall) call init_multipole
      firstcall = .false.

      call invalidate_history

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
      use multigridvars,      only: overrelax, overrelax_xyz

      implicit none

      type(grid_container), pointer,  intent(inout) :: cg
      type(cg_level_connected_T), pointer :: curl
      real, allocatable, dimension(:)  :: kx, ky, kz             !< FFT kernel directional components for convolution
      integer :: i, j

      ! this should work correctly also when dom%eff_dim < 3
      cg%mg%r  = overrelax / 2.
      cg%mg%rx = cg%dvol**2 * cg%idx2
      cg%mg%ry = cg%dvol**2 * cg%idy2
      cg%mg%rz = cg%dvol**2 * cg%idz2
      cg%mg%r  = cg%mg%r / (cg%mg%rx + cg%mg%ry + cg%mg%rz)
      cg%mg%rx = overrelax_xyz(xdim) * cg%mg%rx * cg%mg%r
      cg%mg%ry = overrelax_xyz(ydim) * cg%mg%ry * cg%mg%r
      cg%mg%rz = overrelax_xyz(zdim) * cg%mg%rz * cg%mg%r
      cg%mg%r  = cg%mg%r * cg%dvol**2
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

!> \brief Mark the historical solutions as invalid

   subroutine invalidate_history

      implicit none

      inner%valid = .false.
      outer%valid = .false.

   end subroutine invalidate_history

!>
!! \brief Make a local copy of source (density) and multiply by 4 pi G
!!
!! \details Typically i_all_dens is a copy of fluidindex::iarr_all_sg.
!! Passing this as an argument allows for independent computation of the potential for several density fields if necessary.
!! \todo compactify the following more (if possible)
!<

   subroutine init_source(i_all_dens)

      use cg_list_global, only: all_cg
      use constants,      only: GEO_RPZ, LO, HI, xdim, ydim, zdim, O_I4
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
      logical                        :: apply_src_Mcorrection

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

               apply_src_Mcorrection = any(cg%ext_bnd(:,:)) .and. (ord_laplacian_outer == -O_I4) ! an improvement for Mehrstellen Laplace operator

               if (apply_src_Mcorrection) cg%wa(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) = 0.
               do side = LO, HI
                  if (cg%ext_bnd(xdim, side)) then
                     fac = 2. * cg%idx2 / fpiG
                     if (dom%geometry_type == GEO_RPZ .and. cg%x(cg%ijkse(xdim,side)) /= 0.) fac = fac - 1./(cg%dx * cg%x(cg%ijkse(xdim,side)) * fpiG) !> BEWARE is it cg%x(ie), cg%x(ie+1) or something in the middle?
                     cg%q(source)%arr       (cg%ijkse(xdim,side), cg%js:cg%je, cg%ks:cg%ke) = &
                          & cg%q(source)%arr(cg%ijkse(xdim,side), cg%js:cg%je, cg%ks:cg%ke) - &
                          & cg%mg%bnd_x(                          cg%js:cg%je, cg%ks:cg%ke, side) * fac
                     if (apply_src_Mcorrection) cg%wa(cg%ijkse(xdim,side), cg%js:cg%je, cg%ks:cg%ke) = cg%wa(cg%ijkse(xdim,side), cg%js:cg%je, cg%ks:cg%ke) + 1
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
                     if (apply_src_Mcorrection) cg%wa(cg%is:cg%ie, cg%ijkse(ydim,side), cg%ks:cg%ke) = cg%wa(cg%is:cg%ie, cg%ijkse(ydim,side), cg%ks:cg%ke) + 1
                  endif
               enddo
               do side = LO, HI
                  if (cg%ext_bnd(zdim, side)) then
                     cg%q(source)%arr       (cg%is:cg%ie, cg%js:cg%je, cg%ijkse(zdim,side)) = &
                          & cg%q(source)%arr(cg%is:cg%ie, cg%js:cg%je, cg%ijkse(zdim,side)) - &
                          & cg%mg%bnd_z     (cg%is:cg%ie, cg%js:cg%je, side) * 2. * cg%idz2 / fpiG
                     if (apply_src_Mcorrection) cg%wa(cg%is:cg%ie, cg%js:cg%je, cg%ijkse(zdim,side)) = cg%wa(cg%is:cg%ie, cg%js:cg%je, cg%ijkse(zdim,side)) + 1
                  endif
               enddo
               if (apply_src_Mcorrection) then
                  where (cg%wa(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) == 2.) &
                       cg%q(source)%arr(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) = cg%q(source)%arr(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) * 5./6.
                  where (cg%wa(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) == 3.) &
                       cg%q(source)%arr(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) = cg%q(source)%arr(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) * 2./3.
               endif
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

      call poisson_solver(inner)

      call leaves%q_copy(solution, qna%ind(sgp_n))

      if (grav_bnd_global == bnd_isolated) then
         grav_bnd = bnd_givenval

         vstat%cprefix = "Go-"
         call multipole_solver
         call init_source(empty_array)

         call poisson_solver(outer)

         call leaves%q_add(solution, qna%ind(sgp_n)) ! add solution to sgp

      endif

      grav_bnd = grav_bnd_global
      ts = set_timer("multigrid")
      tot_ts = tot_ts + ts

   end subroutine multigrid_solve_grav

!> \brief Chose the desired poisson solver

   subroutine poisson_solver(history)

      use multigrid_old_soln, only: soln_history
      use multigridvars,      only: grav_bnd, bnd_givenval

      implicit none

      type(soln_history), intent(inout) :: history !< inner or outer potential history used for initializing first guess

      if (grav_bnd == bnd_givenval) then
         if (use_CG_outer) then
            call mgpcg(history)
         else
            call vcycle_hg(history)
         endif
      else
         if (use_CG) then
            call mgpcg(history)
         else
            call vcycle_hg(history)
         endif
      endif

   end subroutine poisson_solver

!>
!! \brief Multigrid-preconditioned conjugate gradient solver
!!
!! \details This solver can converge on domains that require high cell anisotropy e.g cg%dx >> cg%dy
!!
!! Note that forcing alpha = 1 and beta = 0 reduces this algorithm to bare multigrid
!!
!! \todo Find out why it converges so poorly for outer potential
!<

   subroutine mgpcg(history)

      use cg_leaves,          only: leaves
      use cg_list_dataop,     only: ind_val
      use dataio_pub,         only: msg, printinfo
      use mpisetup,           only: master
      use multigrid_Laplace,  only: residual_order, vT_A_v_order
      use multigrid_old_soln, only: soln_history
      use multigridvars,      only: source, solution, defect, correction, tot_ts, ts
      use timer,              only: set_timer

      implicit none

      type(soln_history), intent(inout) :: history !< inner or outer potential history used for initializing first guess

      real    :: alpha, beta, dc
      integer :: it
      real    :: norm_rhs, norm_lhs, norm_old

      beta = 0.
      call history%init_solution(vstat%cprefix)
      norm_rhs = leaves%norm_sq(source)
      norm_old = norm_rhs

      call residual_order(ordL(), leaves, source, solution, defect) !    {r}_0 := {b} - {A x}_0
      call single_v_cycle(defect, correction) !     {z}_0 := {M}^{-1} {r}_0
      call leaves%q_copy(correction, cg_corr) !   {p}_0 := {z}_0
      do it = 0, max_cycles !    k := 0 \,     repeat

         dc = leaves%scalar_product(defect, correction)
         alpha = dc/vT_A_v_order(ordL(), cg_corr) !\alpha_k := {{r}_k^{T} {z}_k}/{{p}_k^{T} {A p}_k}
         call leaves%q_lin_comb( [ ind_val(solution, 1.), ind_val(cg_corr, alpha) ], solution) !{x}_{k+1} := {x}_k + \alpha_k {p}_k
         call residual_order(ordL(), leaves, source, solution, defect) ! {r}_{k+1} := {r}_k - \alpha_k {A p}_k
         norm_lhs = leaves%norm_sq(defect)
         ts = set_timer("multigrid")
         tot_ts = tot_ts + ts
         write(msg,'(a,i3,a,f12.8,a,f9.2,a,f11.7,g14.6,a,f8.3)')" MG-PCG: ", it, " lhs/rhs= ",norm_lhs/norm_rhs, " improvement= ",norm_old/norm_lhs, &
              &                                                 " alpha= ", alpha, beta, " time=", ts
         if (master) call printinfo(msg)
         if (norm_lhs/norm_rhs <= norm_tol) exit ! if rk+1 is sufficiently small then exit loop endif
         norm_old = norm_lhs
         call single_v_cycle(defect, correction) ! {z}_{k+1} := {M}^{-1} {r}_{k+1}
         beta = leaves%scalar_product(defect, correction)/dc ! \beta_k := {{z}_{k+1}^{T} {r}_{k+1}}/{{z}_k^{T} {r}_k}
         call leaves%q_lin_comb( [ ind_val(correction, 1.), ind_val(cg_corr, beta) ], cg_corr ) ! {p}_{k+1} := {z}_{k+1} + \beta_k {p}_k

      enddo !k := k + 1;    end repeat

      call history%store_solution       !    The result is xk+1

   end subroutine mgpcg

!> \brief A single Hueang-Greengard V-cycle

   subroutine single_v_cycle(def, corr)

      use cg_leaves,          only: leaves
      use cg_list_global,     only: all_cg
      use cg_level_coarsest,  only: coarsest
      use cg_level_connected, only: cg_level_connected_T
      use cg_level_finest,    only: finest

      implicit none

      integer(kind=4), intent(in) :: def  !< Defect (source, residual)
      integer(kind=4), intent(in) :: corr !< Approximate solution for the defect (correction)

      type(cg_level_connected_T), pointer :: curl


      if (use_MGpreconditioning) then
         ! the Huang-Greengard V-cycle
         call finest%level%restrict_to_floor_q_1var(def)

         call all_cg%set_dirty(corr)
         call coarsest%level%set_q_value(corr, 0.)

         curl => coarsest%level
         do while (associated(curl))
            call approximate_solution(curl, def, corr)
            call curl%check_dirty(corr, "Vup1 relax+")
            curl => curl%finer
         enddo
      else
         call leaves%q_copy(def, corr) ! non-preconditioned cg
      endif

   end subroutine single_v_cycle

!>
!! \brief The solver. Here we choose an adaptation of the Huang-Greengard V-cycle.
!! For more difficult problems, like variable coefficient diffusion equation a more sophisticated V-cycle may be more effective.
!<

   subroutine vcycle_hg(history)

      use cg_leaves,          only: leaves
      use cg_list_global,     only: all_cg
      use cg_level_coarsest,  only: coarsest
      use cg_level_connected, only: cg_level_connected_T
      use cg_level_finest,    only: finest
      use constants,          only: cbuff_len, fft_none
      use dataio_pub,         only: msg, die, warn, printinfo
      use global,             only: do_ascii_dump
      use mpisetup,           only: master, nproc
      use multigridvars,      only: source, solution, correction, defect, verbose_vcycle, stdout, tot_ts, ts, grav_bnd, bnd_periodic
      use multigrid_Laplace,  only: residual_order
      use multigrid_old_soln, only: soln_history
      use timer,              only: set_timer

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
         call residual_order(ordL(), leaves, source, solution, defect)
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

         if (v>0 .and. norm_old/norm_lhs <= vcycle_giveup) then
            if (master) then
               write(msg, '(a,f6.1)')"[multigrid_gravity:vcycle_hg] Poor convergence detected. Giving up. norm_tol missed by a factor of ",norm_lhs/norm_rhs/norm_tol
               call warn(msg)
            endif
            exit
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
      call residual_order(ordL(), leaves, source, solution, defect)
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

!> \brief Select appropriate order of laplacian, depending on which solve operation we're on

   function ordL()

      use multigridvars, only: grav_bnd, bnd_givenval

      implicit none

      integer(kind=4) :: ordL

      ordL = ord_laplacian
      if (grav_bnd == bnd_givenval) ordL = ord_laplacian_outer

   end function ordL

!> \brief This routine has to find an approximate solution for given source field and implemented differential operator

   subroutine approximate_solution(curl, src, soln)

      use cg_level_coarsest,  only: coarsest
      use cg_level_connected, only: cg_level_connected_T
      use constants,          only: BND_NEGREF
      use multigridvars,      only: correction, nsmool
      use multigrid_Laplace,  only: approximate_solution_order
!!$      use constants,           only: fft_none
!!$      use multigrid_fftapprox, only: approximate_solution_fft

      implicit none

      type(cg_level_connected_T), pointer, intent(in) :: curl !< pointer to a level for which we approximate the solution
      integer(kind=4),                     intent(in) :: src  !< index of source in cg%q(:)
      integer(kind=4),                     intent(in) :: soln !< index of solution in cg%q(:)

      integer :: nsmoo

      call curl%check_dirty(src, "approx_soln src-")
!!$      if (curl%fft_type /= fft_none) then
!!$         call approximate_solution_fft(curl, src, soln)
!!$      else
      if (associated(curl, coarsest%level)) then
         nsmoo = nsmoob
         !> \todo Implement automatic convergence check on coarsest level (not very important when we have a FFT solver for coarsest level)
      else
         nsmoo = nsmool
         if (soln == correction) call curl%coarser%prolong_q_1var(soln, bnd_type = BND_NEGREF) ! make sure that prolongation is called only in ascending (coarse -> fine) part of V-cycle.
         !> \warning this may be incompatible with V-cycles other than Huang - Greengard
      endif

      call approximate_solution_order(ordL(), curl, src, soln, nsmoo)
!!$      endif

      call curl%check_dirty(soln, "approx_soln soln+")

   end subroutine approximate_solution

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
