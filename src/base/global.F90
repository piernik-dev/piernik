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

   use constants, only: cbuff_len, xdim, zdim

   implicit none

   private
   public :: cleanup_global, init_global, &
        &    cfl, cfl_max, cflcontrol, cfl_violated, &
        &    dt, dt_initial, dt_max_grow, dt_min, dt_old, dtm, t, t_saved, nstep, nstep_saved, &
        &    integration_order, limiter, limiter_b, smalld, smallei, smallp, use_smalld, h_solver, interpol_str, &
        &    relax_time, grace_period_passed, cfr_smooth, repeat_step, skip_sweep, geometry25D, &
        &    dirty_debug, do_ascii_dump, show_n_dirtys, no_dirty_checks, sweeps_mgu, use_fargo, print_divB, &
        &    divB_0_method, force_cc_mag, glm_alpha, use_eglm, cfl_glm, ch_grid, w_epsilon

   real, parameter :: dt_default_grow = 2.
   logical         :: cfl_violated             !< True when cfl condition is violated
   logical         :: dirty_debug              !< Allow initializing arrays with some insane values and checking if these values can propagate
   integer(kind=4) :: show_n_dirtys            !< use to limit the amount of printed messages on dirty values found
   logical         :: do_ascii_dump            !< to dump, or not to dump: that is a question (ascii)
   logical         :: no_dirty_checks          !< Temporarily disable dirty checks
   integer(kind=4) :: nstep, nstep_saved
   real            :: t, dt, dt_old, dtm, t_saved
   integer         :: divB_0_method            !< encoded method of making div(B) = 0 (currently DIVB_CT or DIVB_HDC)
   logical         :: force_cc_mag             !< treat magnetic field as cell-centered in the Riemann solver (temporary hack)

   ! Namelist variables

   real    :: dt_initial               !< initial timestep
   real    :: dt_max_grow              !< maximum timestep growth rate
   real    :: dt_min                   !< minimum allowed timestep
   real    :: cfl                      !< desired Courant–Friedrichs–Lewy number
   real    :: cfl_max                  !< warning threshold for the effective CFL number achieved
   logical :: use_smalld               !< correct density when it gets lower than smalld
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
   integer(kind=4), protected    :: integration_order !< Runge-Kutta time integration order (1 - 1st order, 2 - 2nd order)
   character(len=cbuff_len)      :: limiter           !< type of flux limiter
   character(len=cbuff_len)      :: limiter_b         !< type of flux limiter for magnetic field in the Riemann solver
   character(len=cbuff_len)      :: cflcontrol        !< type of cfl control just before each sweep (possibilities: 'none', 'main', 'user')
   character(len=cbuff_len)      :: h_solver          !< type of hydro solver
   character(len=cbuff_len)      :: interpol_str      !< type of interpolation
   character(len=cbuff_len)      :: divB_0            !< human-readable method of making div(B) = 0 (currently CT or HDC)
   logical                       :: repeat_step       !< repeat fluid step if cfl condition is violated (significantly increases mem usage)
   logical, dimension(xdim:zdim) :: skip_sweep        !< allows to skip sweep in chosen direction
   logical                       :: sweeps_mgu        !< Mimimal Guardcell Update in sweeps
   logical                       :: use_fargo         !< use Fast Eulerian Transport for differentially rotating disks
   integer(kind=4)               :: print_divB        !< if >0 then print div(B) estimates each print_divB steps
   real                          :: glm_alpha         !< damping factor for the psi field
   logical                       :: use_eglm          !< use E-GLM?
   real                          :: cfl_glm           !< "CFL" for chspeed in divergence cleaning
   logical                       :: ch_grid           !< When true use grid properties to estimate ch (psi wave propagation speed). Use gas properties otherwise.
   real                          :: w_epsilon         !< small number for safe evaluation of weights in WENO interpolation

   namelist /NUMERICAL_SETUP/ cfl, cflcontrol, cfl_max, use_smalld, smalld, smallei, smallc, smallp, dt_initial, dt_max_grow, dt_min, &
        &                     repeat_step, limiter, limiter_b, relax_time, integration_order, cfr_smooth, skip_sweep, geometry25D, sweeps_mgu, print_divB, &
        &                     use_fargo, h_solver, divB_0, glm_alpha, use_eglm, cfl_glm, ch_grid, interpol_str, w_epsilon

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
!!   <tr><td>cfl              </td><td>0.7    </td><td>real value between 0.0 and 1.0       </td><td>\copydoc global::cfl              </td></tr>
!!   <tr><td>cfl_max          </td><td>0.9    </td><td>real value between cfl and 1.0       </td><td>\copydoc global::cfl_max          </td></tr>
!!   <tr><td>cflcontrol       </td><td>warn   </td><td>string                               </td><td>\copydoc global::cflcontrol       </td></tr>
!!   <tr><td>repeat_step      </td><td>.true. </td><td>logical value                        </td><td>\copydoc global::use_smalld       </td></tr>
!!   <tr><td>smallp           </td><td>1.e-10 </td><td>real value                           </td><td>\copydoc global::smallp           </td></tr>
!!   <tr><td>smalld           </td><td>1.e-10 </td><td>real value                           </td><td>\copydoc global::smalld           </td></tr>
!!   <tr><td>use_smalld       </td><td>.true. </td><td>logical value                        </td><td>\copydoc global::use_smalld       </td></tr>
!!   <tr><td>smallei          </td><td>1.e-10 </td><td>real value                           </td><td>\copydoc global::smallei          </td></tr>
!!   <tr><td>smallc           </td><td>1.e-10 </td><td>real value                           </td><td>\copydoc global::smallc           </td></tr>
!!   <tr><td>integration_order</td><td>2      </td><td>1 or 2 (or 3 - currently unavailable)</td><td>\copydoc global::integration_order</td></tr>
!!   <tr><td>cfr_smooth       </td><td>0.0    </td><td>real value                           </td><td>\copydoc global::cfr_smooth       </td></tr>
!!   <tr><td>dt_initial       </td><td>-1.    </td><td>positive real value or -1.           </td><td>\copydoc global::dt_initial       </td></tr>
!!   <tr><td>dt_max_grow      </td><td>2.     </td><td>real value > 1.1                     </td><td>\copydoc global::dt_max_grow      </td></tr>
!!   <tr><td>dt_min           </td><td>0.     </td><td>positive real value                  </td><td>\copydoc global::dt_min           </td></tr>
!!   <tr><td>limiter          </td><td>vanleer</td><td>string                               </td><td>\copydoc global::limiter          </td></tr>
!!   <tr><td>limiter_b        </td><td>moncen </td><td>string                               </td><td>\copydoc global::limiter_b        </td></tr>
!!   <tr><td>relax_time       </td><td>0.0    </td><td>real value                           </td><td>\copydoc global::relax_time       </td></tr>
!!   <tr><td>skip_sweep       </td><td>F, F, F</td><td>logical array                        </td><td>\copydoc global::skip_sweep       </td></tr>
!!   <tr><td>geometry25D      </td><td>F      </td><td>logical value                        </td><td>\copydoc global::geometry25d      </td></tr>
!!   <tr><td>sweeps_mgu       </td><td>F      </td><td>logical value                        </td><td>\copydoc global::sweeps_mgu       </td></tr>
!!   <tr><td>h_solver         </td><td>"rk2"  </td><td>string                               </td><td>\copydoc global::h_solver         </td></tr>
!!   <tr><td>divB_0           </td><td>CT     </td><td>string                               </td><td>\copydoc global::divB_0           </td></tr>
!!   <tr><td>glm_alpha        </td><td>0.1    </td><td>real value                           </td><td>\copydoc global::glm_alpha        </td></tr>
!!   <tr><td>use_eglm         </td><td>false  </td><td>logical value                        </td><td>\copydoc global::use_eglm         </td></tr>
!!   <tr><td>print_divB       </td><td>0      </td><td>integer value                        </td><td>\copydoc global::print_divB       </td></tr>
!!   <tr><td>ch_grid          </td><td>false  </td><td>logical value                        </td><td>\copydoc global::ch_grid          </td></tr>
!!   <tr><td>w_epsilon        </td><td>1e-10  </td><td>real                                 </td><td>\copydoc global::w_epsilon        </td></tr>
!! </table>
!! \n \n
!<
   subroutine init_global

      use constants,  only: big_float, PIERNIK_INIT_MPI, INVALID, DIVB_CT, DIVB_HDC
      use dataio_pub, only: die, msg, warn, code_progress, printinfo
      use dataio_pub, only: nh  ! QA_WARN required for diff_nml
      use mpisetup,   only: cbuff, ibuff, lbuff, rbuff, master, slave, piernik_MPI_Bcast

      implicit none

      if (code_progress < PIERNIK_INIT_MPI) call die("[global:init_global] MPI not initialized.")

      dt_old = -1.
      t = 0.

      ! Begin processing of namelist parameters

#ifdef RIEMANN
      ! 'moncen' and 'vanleer' seem to be best for emag conservation with GLM for b_limiter
      limiter_b   = 'moncen'
      limiter     = limiter_b
      divB_0      = "HDC"
#else /* ! RIEMANN */
      limiter     = 'vanleer'
      limiter_b   = limiter
      divB_0      = "CT"
#endif /* RIEMANN */
      cflcontrol  = 'warn'
      h_solver    = 'rk2'
      interpol_str = 'linear'
      repeat_step = .true.
      geometry25D = .false.
      no_dirty_checks = .false.
#ifdef MAGNETIC
      sweeps_mgu  = .false.
#else /* !MAGNETIC */
      sweeps_mgu  = .true.
#endif /* !MAGNETIC */

      cfl         = 0.7
      cfl_max     = 0.9
      cfr_smooth  = 0.0
      smallp      = big_float
      smalld      = big_float
      use_smalld  = .true.
      smallc      = 1.e-10
      smallei     = 1.e-10
      dt_initial  = -1.              !< negative value indicates automatic choice of initial timestep
      dt_max_grow = dt_default_grow  !< for sensitive setups consider setting this as low as 1.1
      dt_min      = tiny(1.)
      relax_time  = 0.
      integration_order  = 2
      use_fargo   = .false.
      glm_alpha   = 0.1
      skip_sweep  = .false.
      use_eglm    = .false.
      print_divB  = 100
      cfl_glm     = cfl
      ch_grid     = .false.
      w_epsilon   = 1e-10

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

         ! Sanitize input parameters, if possible
         if (cfl <= 0. .or. cfl >1.0) call die("[global:init_global] CFL value should be >0. and <=1.")
         cfl_max = min(max(cfl_max, min(cfl*1.1, cfl+0.05, (1.+cfl)/2.) ), 1.0) ! automatically sanitize cfl_max
         if (integration_order > 2) call die ('[global:init_global]: "ORIG" scheme integration_order must be 1 or 2')

         if (dt_max_grow < 1.01) then
            write(msg,'(2(a,g10.3))')"[global:init_global] dt_max_grow = ",dt_max_grow," is way too low. Resetting to ",dt_default_grow
            call warn(msg)
            dt_max_grow = dt_default_grow
         endif

         if (h_solver /= "muscl" .and. h_solver /= "rk2" .and. h_solver /= "rk3" ) then
            write(msg, *)"[fluidupdate:sweep_dsplit:warn_experimental] The scheme '", trim(h_solver), "' is experimental. Use 'muscl' or 'rk(n)' for production runs."
            call warn(msg)
         endif

         cbuff(1) = limiter
         cbuff(2) = limiter_b
         cbuff(3) = cflcontrol
         cbuff(4) = h_solver
         cbuff(5) = divB_0
         cbuff(6) = interpol_str

         ibuff(1) = integration_order
         ibuff(2) = print_divB

         rbuff( 1) = smalld
         rbuff( 2) = smallc
         rbuff( 3) = smallp
         rbuff( 4) = smallei
         rbuff( 5) = cfl
         rbuff( 6) = cfr_smooth
         rbuff( 7) = dt_initial
         rbuff( 8) = dt_max_grow
         rbuff( 9) = dt_min
         rbuff(10) = cfl_max
         rbuff(11) = relax_time
         rbuff(12) = glm_alpha
         rbuff(13) = cfl_glm
         rbuff(14) = w_epsilon

         lbuff(1)   = use_smalld
         lbuff(2)   = repeat_step
         lbuff(3:5) = skip_sweep
         lbuff(6)   = geometry25D
         lbuff(7)   = sweeps_mgu
         lbuff(8)   = use_fargo
         lbuff(9)   = use_eglm
         lbuff(10)  = ch_grid

      endif

      call piernik_MPI_Bcast(cbuff, cbuff_len)
      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(rbuff)
      call piernik_MPI_Bcast(lbuff)

      if (slave) then

         use_smalld    = lbuff(1)
         repeat_step   = lbuff(2)
         skip_sweep    = lbuff(3:5)
         geometry25D   = lbuff(6)
         sweeps_mgu    = lbuff(7)
         use_fargo     = lbuff(8)
         use_eglm      = lbuff(9)
         ch_grid       = lbuff(10)

         smalld      = rbuff( 1)
         smallc      = rbuff( 2)
         smallp      = rbuff( 3)
         smallei     = rbuff( 4)
         cfl         = rbuff( 5)
         cfr_smooth  = rbuff( 6)
         dt_initial  = rbuff( 7)
         dt_max_grow = rbuff( 8)
         dt_min      = rbuff( 9)
         cfl_max     = rbuff(10)
         relax_time  = rbuff(11)
         glm_alpha   = rbuff(12)
         cfl_glm     = rbuff(13)
         w_epsilon   = rbuff(14)

         limiter    = cbuff(1)
         limiter_b  = cbuff(2)
         cflcontrol = cbuff(3)
         h_solver   = cbuff(4)
         divB_0     = cbuff(5)
         interpol_str = cbuff(6)

         integration_order = ibuff(1)
         print_divB        = ibuff(2)

      endif

      divB_0_method = INVALID
      select case (divB_0)
         case ("CT", "ct", "constrained transport", "Constrained Transport")
            divB_0_method = DIVB_CT
         case ("HDC", "hdc", "GLM", "glm", "divergence cleaning", "divergence diffusion")
            divB_0_method = DIVB_HDC
            if (master .and. cfl * cfl_glm**2 > 0.25) call warn("[global] To avoid unphysical instabilities with GLM divergence cleaning consider reducing CFL or CFL_GLM, so CFL * CFL_GLM**2 would be < 0.25.")
         case default
            call die("[global:init_global] unrecognized divergence cleaning description.")
      end select

      !> reshape_b should carefully check things here
      force_cc_mag = .false.
      select case (divB_0_method)
         case (DIVB_HDC)
            force_cc_mag = .true.
            if (ch_grid .and. master) call warn("[global] ch_grid = .true. is risky")
         case (DIVB_CT)
            force_cc_mag = .false.
         case default
            call die("[global:init_global] unrecognized divergence cleaning method.")
      end select

#ifdef RTVD
      if (master) call printinfo("    (M)HD solver: RTVD.")
#endif /* RTVD */
#ifdef HLLC
      if (master) call printinfo("    HD solver: HLLC.")
#endif /*HLLC */
#ifdef RIEMANN
      if (master) call printinfo("    (M)HD solver: Riemann.")
#endif /* RIEMANN */

#ifdef MAGNETIC
      if (master) then
         select case (divB_0_method)
             case (DIVB_HDC)
                call printinfo("    The div(B) constraint is maintaineded by Hyperbolic Cleaning (GLM).")
             case (DIVB_CT)
                call printinfo("    The div(B) constraint is maintaineded by Constrained Transport (2nd order).")
             case default
                call die("    The div(B) constraint is maintaineded by Uknown Something.")
         end select

         if (force_cc_mag) then
            call printinfo("    Magnetic field is cell-centered.")
         else
            call printinfo("    Magnetic field is face-centered (staggered).")
         endif
      endif
#endif /* MAGNETIC */

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

!-----------------------------------------------------------------------------

end module global
