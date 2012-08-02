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
        &    dt, dt_initial, dt_max_grow, dt_min, dt_old, dtm, t, nstep, &
        &    integration_order, limiter, smalld, smallei, smallp, use_smalld, magic_mass, local_magic_mass, recent_magic_mass, &
        &    relax_time, grace_period_passed, cfr_smooth, repeat_step, skip_sweep, geometry25D, dirty_debug, do_ascii_dump

   real, parameter :: dt_default_grow = 2.
   logical         :: cfl_violated             !< True when cfl condition is violated
   logical         :: dirty_debug              !< Allow initializing arrays with some insane values and checking if these values can propagate
   logical         :: do_ascii_dump                      !< to dump, or not to dump: that is a question (ascii)
   integer(kind=4) :: nstep
   real            :: t, dt, dt_old, dtm
   real, dimension(:),   allocatable       :: magic_mass
   real, dimension(:),   allocatable, save :: local_magic_mass
   real, dimension(:,:), allocatable, save :: recent_magic_mass

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
   character(len=cbuff_len)      :: cflcontrol        !< type of cfl control just before each sweep (possibilities: 'none', 'main', 'user')
   logical                       :: repeat_step       !< repeat fluid step if cfl condition is violated (significantly increases mem usage)
   logical, dimension(xdim:zdim) :: skip_sweep        !< allows to skip sweep in chosen direction

   namelist /NUMERICAL_SETUP/ cfl, cflcontrol, cfl_max, use_smalld, smalld, smallei, smallc, smallp, dt_initial, dt_max_grow, dt_min, &
        &                     repeat_step, limiter, relax_time, integration_order, cfr_smooth, skip_sweep, geometry25D

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
!!   <tr><td>relax_time       </td><td>0.0    </td><td>real value                           </td><td>\copydoc global::relax_time       </td></tr>
!!   <tr><td>skip_sweep       </td><td>F, F, F</td><td>logical array                        </td><td>\copydoc global::skip_sweep       </td></tr>
!!   <tr><td>geometry25D      </td><td>F      </td><td>logical value                        </td><td>\copydoc global::geometry25D      </td></tr>
!! </table>
!! \n \n
!<
   subroutine init_global

      use constants,  only: big_float, PIERNIK_INIT_MPI
      use dataio_pub, only: die, msg, warn, code_progress
      use dataio_pub, only: par_file, ierrh, namelist_errh, compare_namelist, cmdl_nml, lun  ! QA_WARN required for diff_nml
      use mpi,        only: MPI_CHARACTER, MPI_INTEGER, MPI_DOUBLE_PRECISION, MPI_LOGICAL
      use mpisetup,   only: buffer_dim, cbuff, ibuff, lbuff, rbuff, master, slave, comm, mpi_err, FIRST

      implicit none

      if (code_progress < PIERNIK_INIT_MPI) call die("[global:init_global] MPI not initialized.")

      dt_old = -1.

      ! Begin processing of namelist parameters

      limiter     = 'vanleer'
      cflcontrol  = 'warn'
      repeat_step = .true.
      geometry25D = .false.

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

      if (master) then
         diff_nml(NUMERICAL_SETUP)

         ! Sanitize input parameters, if possible
         if (cfl <= 0. .or. cfl >1.0) call die("[global:init_global] CFL value should be >0. and <=1.")
         cfl_max = min(max(cfl_max, min(cfl*1.1, cfl+0.05, (1.+cfl)/2.) ), 1.0) ! automatically sanitize cfl_max
         if (integration_order > 2) call die ('[global:init_global]: "ORIG" scheme integration_order must be 1 or 2')

         if (dt_max_grow < 1.01) then
            write(msg,'(2(a,g10.3))')"[global:init_global] dt_max_grow = ",dt_max_grow," is way too low. Resetting to ",dt_default_grow
            call warn(msg)
            dt_max_grow = dt_default_grow
         endif

         cbuff(1) = limiter
         cbuff(2) = cflcontrol

         ibuff(1) = integration_order

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

         lbuff(1)   = use_smalld
         lbuff(2)   = repeat_step
         lbuff(3:5) = skip_sweep
         lbuff(6)   = geometry25D

      endif

      call MPI_Bcast(cbuff, cbuff_len*buffer_dim, MPI_CHARACTER,        FIRST, comm, mpi_err)
      call MPI_Bcast(ibuff,           buffer_dim, MPI_INTEGER,          FIRST, comm, mpi_err)
      call MPI_Bcast(rbuff,           buffer_dim, MPI_DOUBLE_PRECISION, FIRST, comm, mpi_err)
      call MPI_Bcast(lbuff,           buffer_dim, MPI_LOGICAL,          FIRST, comm, mpi_err)

      if (slave) then

         use_smalld    = lbuff(1)
         repeat_step   = lbuff(2)
         skip_sweep    = lbuff(3:5)
         geometry25D   = lbuff(6)

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

         limiter    = cbuff(1)
         cflcontrol = cbuff(2)

         integration_order = ibuff(1)

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

!-----------------------------------------------------------------------------

end module global
