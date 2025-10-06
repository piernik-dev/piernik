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
!>
!! \brief This module provides global simulation variables such as t or nstep and some numerical parameters, like cfl
!!
!! In this module following namelists of parameters are specified:
!! \copydetails global::init_global
!<
module global

   use barrier,      only: extra_barriers
   use constants,    only: cbuff_len, xdim, zdim
   use mpi_wrappers, only: MPI_wrapper_stats

   implicit none

   private
   public :: cleanup_global, init_global, &
        &    cfl, cfl_max, cflcontrol, disallow_negatives, disallow_CRnegatives, dn_negative, ei_negative, cr_negative, tstep_attempt, &
        &    dt, dt_initial, dt_max_grow, dt_shrink, dt_cur_shrink, dt_min, dt_max, dt_old, dt_full, dtm, t, t_saved, nstep, nstep_saved, max_redostep_attempts, &
        &    repetitive_steps, integration_order, limiter, limiter_b, smalld, smallei, smallp, use_smalld, use_smallei, interpol_str, &
        &    relax_time, grace_period_passed, cfr_smooth, skip_sweep, geometry25D, &
        &    dirty_debug, do_ascii_dump, show_n_dirtys, no_dirty_checks, sweeps_mgu, use_fargo, print_divB, do_external_corners, prefer_merged_MPI, waitall_timeout, &
        &    divB_0_method, cc_mag, glm_alpha, use_eglm, cfl_glm, ch_grid, w_epsilon, psi_bnd, ord_mag_prolong, ord_fluid_prolong, which_solver, use_uhi, is_split


   logical         :: dn_negative = .false.
   logical         :: ei_negative = .false.
   logical         :: cr_negative = .false.
   logical         :: disallow_negatives, disallow_CRnegatives
   logical         :: repetitive_steps         !< repetitve fluid step if cfl condition is violated (significantly increases mem usage)
   logical         :: dirty_debug              !< Allow initializing arrays with some insane values and checking if these values can propagate
   integer(kind=4) :: show_n_dirtys            !< use to limit the amount of printed messages on dirty values found
   logical         :: do_ascii_dump            !< to dump, or not to dump: that is a question (ascii)
   logical         :: no_dirty_checks          !< Temporarily disable dirty checks
   integer(kind=4) :: nstep, nstep_saved
   real            :: t, dt, dt_old, dtm, t_saved
   real            :: dt_full                  !< timestep value which is a subject of shrinking while repeating step
   real            :: dt_cur_shrink            !< currently used dt and CFL number shrinking (used while redoing step)
   integer         :: divB_0_method            !< encoded method of making div(B) = 0 (currently DIVB_CT or DIVB_HDC)
   logical         :: cc_mag                   !< use cell-centered magnetic field
   integer(kind=4) :: psi_bnd                  !< BND_INVALID or enforce some other psi boundary
   integer         :: tstep_attempt            !< /= 0 when we retry timesteps
   integer         :: which_solver             !< one of RTVD_SPLIT, HLLC_SPLIT, RIEMANN_SPLIT or RIEMANN_UNSPLIT
   logical         :: is_split                 !< do we use directional splitting or not?
   logical         :: use_uhi = .false.        !< .true. ⇒ apply BCs to uhi
   real, parameter :: cfl_unsplit = 0.3        !< maximum recommended CFL factor for the unsplit solver

   ! Namelist variables

   real    :: dt_initial               !< if >0. : initial timestep; if 0. or < -1. : automatic timestep; reduced automatic timestep otherwise
   real    :: dt_max_grow              !< maximum timestep growth rate
   real    :: dt_shrink                !< dt shrink rate when timestep retry is used
   real    :: dt_min                   !< minimum allowed timestep
   real    :: dt_max                   !< maximum allowed timestep
   real    :: cfl                      !< desired Courant–Friedrichs–Lewy number
   real    :: cfl_max                  !< warning threshold for the effective CFL number achieved
   integer(kind=4) :: max_redostep_attempts  !< limitation for a number of redoing step attempts (Note: Something might be terribly wrong if a single step requires too many reductions)
   logical :: use_smalld               !< correct density when it gets lower than smalld
   logical :: use_smallei              !< correct internal energy density when it gets lower then smallei
   logical :: geometry25D              !< include source terms in reduced dimension for 2D simulations
   real    :: smallp                   !< artificial infimum for pressure
   real    :: smalld                   !< artificial infimum for density
   real    :: smallc                   !< artificial infimum for freezing speed
   real    :: smallei                  !< artificial infimum for internal energy density
   !>
   !! small number used to smooth freezing speed, especially handy in dust with random noise in velocity field.
   !! \f$c_{\textrm{fr}} = \sqrt{v^2 + \frac{1}{2}(\max{v} - \min{v})c_{\textrm{fr}}^{\textrm{smooth}}} + \ldots\f$
   !<
   real    :: cfr_smooth
   real    :: relax_time                              !< relaxation/grace time, additional physics will be turned off until global::t >= global::relax_time
   integer(kind=4), protected    :: integration_order !< Runge-Kutta time integration order (1 - 1st order (Euler), 2 - 2nd order (RK2))
   character(len=cbuff_len)      :: limiter           !< type of flux limiter
   character(len=cbuff_len)      :: limiter_b         !< type of flux limiter for magnetic field in the Riemann solver
   character(len=cbuff_len)      :: cflcontrol        !< type of cfl control just before/after each sweep (possibilities: 'none', 'warn', 'redo', 'flex', 'auto')
   character(len=cbuff_len)      :: interpol_str      !< type of interpolation
   character(len=cbuff_len)      :: divB_0            !< human-readable method of making div(B) = 0 (currently CT or HDC)
   character(len=cbuff_len)      :: psi_bnd_str       !< "default" for general boundaries or override ith something special
   logical, dimension(xdim:zdim) :: skip_sweep        !< allows to skip sweep in chosen direction
   logical                       :: sweeps_mgu        !< Mimimal Guardcell Update in sweeps
   logical                       :: use_fargo         !< use Fast Eulerian Transport for differentially rotating disks
   integer(kind=4)               :: print_divB        !< if >0 then print div(B) estimates each print_divB steps
   real                          :: glm_alpha         !< damping factor for the psi field
   logical                       :: use_eglm          !< use E-GLM?
   real                          :: cfl_glm           !< "CFL" for chspeed in divergence cleaning
   logical                       :: ch_grid           !< When true use grid properties to estimate ch (psi wave propagation speed). Use gas properties otherwise.
   real                          :: w_epsilon         !< small number for safe evaluation of weights in WENO interpolation
   integer(kind=4)               :: ord_mag_prolong   !< prolongation order for B and psi
   integer(kind=4)               :: ord_fluid_prolong !< prolongation order for u
   logical                       :: do_external_corners  !< when .true. then perform boundary exchanges inside external guardcells
   character(len=cbuff_len)      :: solver_str           !< allow to switch between RIEMANN and RTVD without recompilation

   namelist /NUMERICAL_SETUP/ cfl, cflcontrol, disallow_negatives, disallow_CRnegatives, cfl_max, use_smalld, use_smallei, smalld, smallei, smallc, smallp, dt_initial, dt_max_grow, dt_shrink, dt_min, dt_max, &
        &                     max_redostep_attempts, limiter, limiter_b, relax_time, integration_order, cfr_smooth, skip_sweep, geometry25D, sweeps_mgu, print_divB, &
        &                     use_fargo, divB_0, glm_alpha, use_eglm, cfl_glm, ch_grid, interpol_str, w_epsilon, psi_bnd_str, ord_mag_prolong, ord_fluid_prolong, do_external_corners, solver_str

   logical :: prefer_merged_MPI  !< prefer internal_boundaries_MPI_merged over internal_boundaries_MPI_1by1
   real :: waitall_timeout       !< when > 0. then replace MPI_Waitall with MPI_Test* calls and print some diagnostics it the timeout is reached

   namelist /PARALLEL_SETUP/ extra_barriers, prefer_merged_MPI, MPI_wrapper_stats, waitall_timeout

contains

!-----------------------------------------------------------------------------
!>
!! \brief Routine to set up global properties of the simulation
!!
!! \n \n
!! @b NUMERICAL_SETUP
!! \n \n
!! <table border="+1">
!!   <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!!   <tr><td>cfl                  </td><td>0.7      </td><td>real value between 0.0 and 1.0       </td><td>\copydoc global::cfl                  </td></tr>
!!   <tr><td>cfl_max              </td><td>0.9      </td><td>real value between cfl and 1.0       </td><td>\copydoc global::cfl_max              </td></tr>
!!   <tr><td>cflcontrol           </td><td>redo     </td><td>string                               </td><td>\copydoc global::cflcontrol           </td></tr>
!!   <tr><td>max_redostep_attempts</td><td>10       </td><td>integer                              </td><td>\copydoc global::max_redostep_attempts</td></tr>
!!   <tr><td>smallp               </td><td>1.e-10   </td><td>real value                           </td><td>\copydoc global::smallp               </td></tr>
!!   <tr><td>smalld               </td><td>1.e-10   </td><td>real value                           </td><td>\copydoc global::smalld               </td></tr>
!!   <tr><td>use_smalld           </td><td>.true.   </td><td>logical value                        </td><td>\copydoc global::use_smalld           </td></tr>
!!   <tr><td>smallei              </td><td>1.e-10   </td><td>real value                           </td><td>\copydoc global::smallei              </td></tr>
!!   <tr><td>smallc               </td><td>1.e-10   </td><td>real value                           </td><td>\copydoc global::smallc               </td></tr>
!!   <tr><td>integration_order    </td><td>2        </td><td>1 or 2                               </td><td>\copydoc global::integration_order    </td></tr>
!!   <tr><td>cfr_smooth           </td><td>0.0      </td><td>real value                           </td><td>\copydoc global::cfr_smooth           </td></tr>
!!   <tr><td>dt_initial           </td><td>-1.      </td><td>positive real value or -1. .. 0.     </td><td>\copydoc global::dt_initial           </td></tr>
!!   <tr><td>dt_max_grow          </td><td>2.       </td><td>real value, should be > 1.           </td><td>\copydoc global::dt_max_grow          </td></tr>
!!   <tr><td>dt_shrink            </td><td>0.5      </td><td>real value, should be < 1.           </td><td>\copydoc global::dt_shrink            </td></tr>
!!   <tr><td>dt_min               </td><td>0.       </td><td>positive real value                  </td><td>\copydoc global::dt_min               </td></tr>
!!   <tr><td>dt_max               </td><td>0.       </td><td>positive real value                  </td><td>\copydoc global::dt_max               </td></tr>
!!   <tr><td>limiter              </td><td>vanleer  </td><td>string                               </td><td>\copydoc global::limiter              </td></tr>
!!   <tr><td>limiter_b            </td><td>moncen   </td><td>string                               </td><td>\copydoc global::limiter_b            </td></tr>
!!   <tr><td>relax_time           </td><td>0.0      </td><td>real value                           </td><td>\copydoc global::relax_time           </td></tr>
!!   <tr><td>skip_sweep           </td><td>F, F, F  </td><td>logical array                        </td><td>\copydoc global::skip_sweep           </td></tr>
!!   <tr><td>geometry25D          </td><td>F        </td><td>logical value                        </td><td>\copydoc global::geometry25d          </td></tr>
!!   <tr><td>sweeps_mgu           </td><td>F        </td><td>logical value                        </td><td>\copydoc global::sweeps_mgu           </td></tr>
!!   <tr><td>divB_0               </td><td>CT       </td><td>string                               </td><td>\copydoc global::divB_0               </td></tr>
!!   <tr><td>glm_alpha            </td><td>0.1      </td><td>real value                           </td><td>\copydoc global::glm_alpha            </td></tr>
!!   <tr><td>use_eglm             </td><td>false    </td><td>logical value                        </td><td>\copydoc global::use_eglm             </td></tr>
!!   <tr><td>print_divB           </td><td>100      </td><td>integer value                        </td><td>\copydoc global::print_divB           </td></tr>
!!   <tr><td>ch_grid              </td><td>true     </td><td>logical value                        </td><td>\copydoc global::ch_grid              </td></tr>
!!   <tr><td>w_epsilon            </td><td>1e-10    </td><td>real                                 </td><td>\copydoc global::w_epsilon            </td></tr>
!!   <tr><td>psi_bnd_str          </td><td>"default"</td><td>string                               </td><td>\copydoc global::psi_bnd_str          </td></tr>
!!   <tr><td>ord_mag_prolong      </td><td>2        </td><td>integer                              </td><td>\copydoc global::ord_mag_prolong      </td></tr>
!!   <tr><td>ord_fluid_prolong    </td><td>0        </td><td>integer                              </td><td>\copydoc global::ord_fluid_prolong    </td></tr>
!!   <tr><td>do_external_corners  </td><td>.false.  </td><td>logical                              </td><td>\copydoc global::do_external_corners  </td></tr>
!! </table>
!! \n \n
!! \n \n
!! @b PARALLEL_SETUP
!! \n \n
!! <table border="+1">
!!   <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!!   <tr><td>prefer_merged_MPI </td><td>.true.  </td><td>logical </td><td>\copydoc global::prefer_merged_MPI      </td></tr>
!!   <tr><td>extra_barriers    </td><td>.false. </td><td>logical </td><td>\copydoc mpi_wrapper::extra_barriers    </td></tr>
!!   <tr><td>MPI_wrapper_stats </td><td>.false. </td><td>logical </td><td>\copydoc mpi_wrapper::MPI_wrapper_stats </td></tr>
!!   <tr><td>waitall_timeout   </td><td>0.      </td><td>real    </td><td>\copydoc mpi_wrapper::waitall_timeout   </td></tr>
!! </table>
!! \n \n

!<
   subroutine init_global

      use bcast,      only: piernik_MPI_Bcast
      use constants,  only: big_float, one, PIERNIK_INIT_DOMAIN, INVALID, DIVB_CT, DIVB_HDC, &
           &                BND_INVALID, BND_ZERO, BND_REF, BND_OUT, I_ZERO, O_INJ, O_LIN, O_I2, INVALID, &
           &                RTVD_SPLIT, HLLC_SPLIT, RIEMANN_SPLIT, RIEMANN_UNSPLIT, GEO_XYZ, V_INFO, V_DEBUG, V_ESSENTIAL
      use dataio_pub, only: die, msg, warn, code_progress, printinfo, nh
      use domain,     only: dom
      use mpisetup,   only: cbuff, ibuff, lbuff, rbuff, master, slave

      implicit none

      if (code_progress < PIERNIK_INIT_DOMAIN) call die("[global:init_global] MPI not initialized.")

      dirty_debug = .false.
      dt_old      = -1.
      t           = 0.

      ! Begin processing of namelist parameters

      which_solver      = RIEMANN_SPLIT
      divB_0            = "HDC"  ! This is the default for the Riemann solver, for RTVD it will be changed to "CT" anyway

      ! For RIEMANN_SPLIT 'moncen' and 'vanleer' seem to be best for emag conservation with GLM for b_limiter
      ! Leave RTVD defaults as they were before the implementation of the Riemann HLLD solver
      ! limiter_b   = 'moncen'
      ! limiter     = limiter_b
      limiter     = 'vanleer'
      limiter_b   = limiter

#ifdef NBODY
      cflcontrol  = 'warn'
#else /* !NBODY */
      cflcontrol  = 'redo'
#endif /* !NBODY */
      interpol_str = 'linear'

      geometry25D = .false.
      no_dirty_checks = .false.
#ifdef MAGNETIC
      sweeps_mgu  = .false.
      print_divB  = 100
#else /* !MAGNETIC */
      sweeps_mgu  = .true.
      print_divB  = 0
#endif /* !MAGNETIC */

      cfl         = 0.7
      cfl_max     = 0.9
      cfr_smooth  = 0.0
      smallp      = big_float
      smalld      = big_float
      disallow_negatives = .true.
      disallow_CRnegatives = .false.
      use_smalld  = .true.
      use_smallei = .true.
      smallc      = 1.e-10
      smallei     = 1.e-10
      dt_initial  = -1.              !< -1. indicates automatic choice of initial timestep, -0.5 would give half of that
      dt_max_grow = 2.               !< for sensitive setups consider setting this as low as 1.1
      dt_shrink   = 0.5
      dt_min      = tiny(1.)
      dt_max      = huge(1.)
      relax_time  = 0.
      use_fargo   = .false.
      glm_alpha   = 0.1
      skip_sweep  = .false.
      use_eglm    = .false.
      cfl_glm     = cfl
      ch_grid     = .true.
      w_epsilon   = 1e-10
      psi_bnd_str = "default"
      integration_order  = 2
      max_redostep_attempts = 10
      ord_mag_prolong = O_I2           !< it looks like most f/c artifacts are gone just with cubic prolongation of magnetic guardcells
      ord_fluid_prolong = O_INJ        !< O_INJ and O_LIN ensure monotonicity and nonnegative density and energy
      do_external_corners =.false.
      solver_str = ""

      prefer_merged_MPI = .true.  ! Non-merged MPI in internal_boundaries are implemented without buffers, which can be faster, especially for bsize(:) larger than 3*16, but in some non-periodic setups internal_boundaries_MPI_1by1 has tag collisions, so merged_MPI is currently safer.
      MPI_wrapper_stats = .false.
      waitall_timeout   = 0.

      if (master) then

         if (.not.nh%initialized) call nh%init()

         open(newunit=nh%lun, file=nh%tmp1, status="unknown")
         write(nh%lun,nml=NUMERICAL_SETUP)
         close(nh%lun)
         open(newunit=nh%lun, file=nh%par_file)
         nh%errstr=""
         read(unit=nh%lun, nml=NUMERICAL_SETUP, iostat=nh%ierrh, iomsg=nh%errstr)
         close(nh%lun)
         call nh%namelist_errh(nh%ierrh, "NUMERICAL_SETUP")
         read(nh%cmdl_nml,nml=NUMERICAL_SETUP, iostat=nh%ierrh)
         call nh%namelist_errh(nh%ierrh, "NUMERICAL_SETUP", .true.)
         open(newunit=nh%lun, file=nh%tmp2, status="unknown")
         write(nh%lun,nml=NUMERICAL_SETUP)
         close(nh%lun)
         call nh%compare_namelist()

         open(newunit=nh%lun, file=nh%tmp1, status="unknown")
         write(nh%lun,nml=PARALLEL_SETUP)
         close(nh%lun)
         open(newunit=nh%lun, file=nh%par_file)
         nh%errstr=""
         read(unit=nh%lun, nml=PARALLEL_SETUP, iostat=nh%ierrh, iomsg=nh%errstr)
         close(nh%lun)
         call nh%namelist_errh(nh%ierrh, "PARALLEL_SETUP")
         read(nh%cmdl_nml,nml=PARALLEL_SETUP, iostat=nh%ierrh)
         call nh%namelist_errh(nh%ierrh, "PARALLEL_SETUP", .true.)
         open(newunit=nh%lun, file=nh%tmp2, status="unknown")
         write(nh%lun,nml=PARALLEL_SETUP)
         close(nh%lun)
         call nh%compare_namelist()

         ! Sanitize input parameters, if possible
         if (cfl <= 0. .or. cfl >1.0) call die("[global:init_global] CFL value should be >0. and <=1.")
         cfl_max = min(max(cfl_max, min(cfl*1.1, cfl+0.05, (1.+cfl)/2.) ), 1.0) ! automatically sanitize cfl_max
         if (integration_order > 2) call die ('[global:init_global]: "ORIG" scheme integration_order must be 1 or 2')

         if (dt_max_grow <= 1.01) then
            write(msg,'(2(a,g10.3))')"[global:init_global] dt_max_grow = ", dt_max_grow, " is low. Recommended values are in 1.1 .. 2.0 range."
            call warn(msg)
         endif

         if (dt_shrink > 0.99 .or. dt_shrink < 0.1) then
            write(msg,'(2(a,g10.3))')"[global:init_global] dt_shrink = ", dt_shrink, " is strange. Recommended values are in 0.1 .. 0.9 range."
            call warn(msg)
         endif

         cbuff(1) = limiter
         cbuff(2) = limiter_b
         cbuff(3) = cflcontrol
         cbuff(5) = divB_0
         cbuff(6) = interpol_str
         cbuff(7) = psi_bnd_str
         cbuff(8) = solver_str

         ibuff(1) = integration_order
         ibuff(2) = print_divB
         ibuff(3) = ord_mag_prolong
         ibuff(4) = ord_fluid_prolong
         ibuff(5) = max_redostep_attempts

         rbuff( 1) = smalld
         rbuff( 2) = smallc
         rbuff( 3) = smallp
         rbuff( 4) = smallei
         rbuff( 5) = cfl
         rbuff( 6) = cfr_smooth
         rbuff( 7) = dt_initial
         rbuff( 8) = dt_max_grow
         rbuff( 9) = dt_min
         rbuff(10) = dt_max
         rbuff(11) = cfl_max
         rbuff(12) = relax_time
         rbuff(13) = glm_alpha
         rbuff(14) = cfl_glm
         rbuff(15) = w_epsilon
         rbuff(16) = dt_shrink
         rbuff(17) = waitall_timeout

         lbuff(1)   = use_smalld
         lbuff(2)   = use_smallei
         lbuff(3:5) = skip_sweep
         lbuff(6)   = geometry25D
         lbuff(7)   = sweeps_mgu
         lbuff(8)   = use_fargo
         lbuff(9)   = use_eglm
         lbuff(10)  = ch_grid
         lbuff(11)  = do_external_corners
         lbuff(12)  = disallow_negatives
         lbuff(13)  = disallow_CRnegatives
         lbuff(14)  = prefer_merged_MPI
         lbuff(15)  = extra_barriers
         lbuff(16)  = MPI_wrapper_stats

      endif

      call piernik_MPI_Bcast(cbuff, cbuff_len)
      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(rbuff)
      call piernik_MPI_Bcast(lbuff)

      if (slave) then

         use_smalld            = lbuff(1)
         use_smallei           = lbuff(2)
         skip_sweep            = lbuff(3:5)
         geometry25D           = lbuff(6)
         sweeps_mgu            = lbuff(7)
         use_fargo             = lbuff(8)
         use_eglm              = lbuff(9)
         ch_grid               = lbuff(10)
         do_external_corners   = lbuff(11)
         disallow_negatives    = lbuff(12)
         disallow_CRnegatives  = lbuff(13)
         prefer_merged_MPI     = lbuff(14)
         extra_barriers        = lbuff(15)
         MPI_wrapper_stats     = lbuff(16)

         smalld                = rbuff( 1)
         smallc                = rbuff( 2)
         smallp                = rbuff( 3)
         smallei               = rbuff( 4)
         cfl                   = rbuff( 5)
         cfr_smooth            = rbuff( 6)
         dt_initial            = rbuff( 7)
         dt_max_grow           = rbuff( 8)
         dt_min                = rbuff( 9)
         dt_max                = rbuff(10)
         cfl_max               = rbuff(11)
         relax_time            = rbuff(12)
         glm_alpha             = rbuff(13)
         cfl_glm               = rbuff(14)
         w_epsilon             = rbuff(15)
         dt_shrink             = rbuff(16)
         waitall_timeout       = rbuff(17)

         limiter               = cbuff(1)
         limiter_b             = cbuff(2)
         cflcontrol            = cbuff(3)
         divB_0                = cbuff(5)
         interpol_str          = cbuff(6)
         psi_bnd_str           = cbuff(7)
         solver_str            = cbuff(8)

         integration_order     = ibuff(1)
         print_divB            = ibuff(2)
         ord_mag_prolong       = ibuff(3)
         ord_fluid_prolong     = ibuff(4)
         max_redostep_attempts = ibuff(5)

      endif

      is_split = .true.
      select case (solver_str)
         case ("")  ! leave the default
         case ("rtvd", "RTVD")
            which_solver = RTVD_SPLIT
         case ("hllc", "HLLC")
            which_solver = HLLC_SPLIT
            call warn("[global:init_global] The HLLC solver is not maintained and may be less acurate or nonfunctional on some setups. Don't use it for production runs.")
         case ("riemann", "Riemann", "RIEMANN", "RIEMANN_SPLIT", "riemann_split")
            which_solver = RIEMANN_SPLIT
         case ("UNSPLIT", "unsplit", "van_leer", "Riemann_unsplit", "RIEMANN_UNSPLIT", "riemann_unsplit")  ! This is not the most clear way to choose the solver
            which_solver = RIEMANN_UNSPLIT
            is_split = .false.
         case default
            call die("[global:init_global] unrecognized solver: '" // trim(solver_str) // "'")
      end select

      select case (which_solver)
         case (RTVD_SPLIT)
            divB_0 = "CT"  ! no other option
         case (RIEMANN_SPLIT)
            if (dom%geometry_type /= GEO_XYZ) call die("[global:init_global] Riemann solver is implemented only for cartesian geometry")
         case (HLLC_SPLIT)
#ifdef MAGNETIC
            call die("[global:init_global] MAGNETIC not compatible with HLLC")
#endif /* MAGNETIC */
         case (RIEMANN_UNSPLIT)
            if (dom%geometry_type /= GEO_XYZ) call die("[global:init_global] Unsplit riemann solver is implemented only for cartesian geometry")
         case default
            call die("[global:init_global] no solvers defined")
      end select

#ifndef MAGNETIC
      if (print_divB > 0) call warn("[global:init_global] No magnetic field: printing div(B) will be ignored.")
#endif /* !MAGNETIC */

#ifdef CORIOLIS
      if (which_solver /= RTVD_SPLIT) call die("[global:init_global] CORIOLIS has been implemented only for RTVD so far.")
#endif /* CORIOLIS */

      divB_0_method = INVALID
      select case (divB_0)
         case ("CT", "ct", "constrained transport", "Constrained Transport")
            divB_0_method = DIVB_CT
         case ("HDC", "hdc", "GLM", "glm", "divergence cleaning", "divergence diffusion")
            divB_0_method = DIVB_HDC
            if (master .and. .false.) call warn("[global:init_global] In case of problems with stability connected with checkerboard pattern in the psi field consider reducing CFL parameter (or just CFL_GLM). This solver also doesn't like sudden changes of timestep length.")
            ! ToDo: create a way to add this to the crash message.
#ifdef RESISTIVE
            call die("[global:init_global] RESISTIVE not yet implemented for DIVB_HDC")
#endif /* RESISTIVE */
         case default
            call die("[global:init_global] unrecognized divergence cleaning description.")
      end select

      if ((which_solver == RTVD_SPLIT) .and. (divB_0_method /= DIVB_CT)) then
         if (master) call warn("[global:init_global] RTVD works only with Constrained Transport. Enforcing.")
         divB_0_method = DIVB_CT
      endif

      if ((which_solver == RTVD_SPLIT) .and. (divB_0_method /= DIVB_CT)) then
         if (master) call warn("[global:init_global] RTVD works only with Constrained Transport. Enforcing.")
         divB_0_method = DIVB_CT
      endif

      !> reshape_b should carefully check things here
      cc_mag = .false.
      select case (divB_0_method)
         case (DIVB_HDC)
            cc_mag = .true.
         case (DIVB_CT)
            cc_mag = .false.
         case default
            call die("[global:init_global] unrecognized divergence cleaning method.")
      end select

      select case (psi_bnd_str)
         case ('default')
            psi_bnd = BND_INVALID  ! special value,; means: do not override domain boundaries
         case ('zero')
            psi_bnd = BND_ZERO
         case ('ref', 'refl', 'reflecting')
            psi_bnd = BND_REF
         case ('out', 'free')
            psi_bnd = BND_OUT
         case default
            call die("[global:init_global] unrecognized psi boundaries")
      end select

      if (master) then
         select case (which_solver)
            case (RTVD_SPLIT)
               call printinfo("    (M)HD solver: RTVD.", V_INFO)
            case (HLLC_SPLIT)
               call printinfo("    HD solver: HLLC.", V_INFO)
            case (RIEMANN_SPLIT)
               call printinfo("    (M)HD solver:Split Riemann.", V_INFO)
            case (RIEMANN_UNSPLIT)
               call printinfo("    (M)HD solver:Unsplit Riemann.", V_INFO)
            case default
               call die("[global:init_global] unrecognized hydro solver")
         end select

         if (.not. is_split) then
            if (cfl > 0.5) call warn("[global:init_global] Unsplit MHD solver chosen. CFL > 0.5 may lead to unexpected result.")
         endif
#ifdef MAGNETIC
         if (.not. is_split) then
            if (cfl > cfl_unsplit) call warn("[global:init_global] Unsplit MHD solver may be unstable with CFL > 0.3")
            if (cfl_glm > cfl_unsplit) call warn("[global:init_global] Unsplit MHD solver chosen. Ideal CFL_GLM = 0.3. Anything else may lead to unexpected result.")
         endif
#endif /* MAGNETIC */
      endif

#ifdef MAGNETIC
      if (master) then
         select case (divB_0_method)
            case (DIVB_HDC)
               call printinfo("    The div(B) constraint is maintained by Hyperbolic Cleaning (GLM).", V_INFO)
            case (DIVB_CT)
               call printinfo("    The div(B) constraint is maintained by Constrained Transport (2nd order).", V_INFO)
            case default
               call die("    The div(B) constraint is maintained by Uknown Something.")
         end select

         if (cc_mag) then
            call printinfo("    Magnetic field is cell-centered.", V_INFO)
         else
            call printinfo("    Magnetic field is face-centered (staggered).", V_INFO)
         endif
      endif
#endif /* MAGNETIC */

      if (master) &
#  ifdef MPIF08
           call printinfo("    use mpi_f08 (modern interface)", V_DEBUG)
#  else /* !MPIF08 */
           call printinfo("    use mpi (old interface)", V_DEBUG)
#  endif /* !MPIF08 */

      if (all(ord_fluid_prolong /= [O_INJ, O_LIN])) then
         write(msg, '(a,i3,a)')"[global:init_global] Prolongation order ", ord_fluid_prolong, " is not positive-definite and thus not allowed for density and energy. Degrading to injection (0)"
         if (master) call warn(msg)
         ord_fluid_prolong = O_INJ
         ! ToDo implement higher order monotonized prolongation scheme
      endif

      if (master .and. ord_fluid_prolong /= O_INJ) call warn("[global:init_global] Linear prolongation of fluid in AMR is experimental.")

      tstep_attempt = I_ZERO
      dt_cur_shrink = one
      dt = 0.

      if (waitall_timeout > 0. .and. master) then
         write(msg, '(a,g0.2,a)')"[global:init_global] Timeout for MPI_Waitall is ", waitall_timeout, " s."
         call printinfo(msg, V_ESSENTIAL)
      endif

   end subroutine init_global

!-----------------------------------------------------------------------------

   subroutine cleanup_global

      implicit none

   end subroutine cleanup_global

!-----------------------------------------------------------------------------

   logical function grace_period_passed()
      implicit none
      grace_period_passed = (t >= relax_time)
   end function grace_period_passed

end module global
