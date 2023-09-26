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
!! \brief Module containing all main subroutines and functions that govern %gravity force in the code
!!
!!
!! In this module a namelist of parameters is specified:
!! \copydetails gravity::init_grav
!<
module gravity
! pulled by GRAV

   use constants, only: cbuff_len, ndims

   implicit none

   private
   public :: init_grav, init_terms_grav, grav_accel, source_terms_grav, grav_src_exec, grav_pot_3d, grav_type, get_gprofs, compute_h_gpot, update_gp, need_update
   public :: r_gc, ptmass, ptm_x, ptm_y, ptm_z, r_smooth, nsub, tune_zeq, tune_zeq_bnd, r_grav, n_gravr, user_grav, gprofs_target, ptm2_x, variable_gp

   integer, parameter         :: gp_stat_len   = 9
   integer, parameter         :: gproft_len    = 5
   character(len=gp_stat_len) :: gp_status             !< variable set as 'undefined' in grav_pot_3d when grav_accel is supposed to use
   character(len=gproft_len)  :: gprofs_target         !< variable set pointing gravity routine in hydrostatic_zeq ('accel' or ready gp array 'extgp')
   character(len=cbuff_len)   :: external_gp           !< variable allowing to choose external gravitational potential
   real, dimension(ndims)     :: g_dir                 !< vector used by GRAV_UNIFORM and GRAV_LINEAR type of %gravity
   real                       :: r_gc                  !< galactocentric radius of the local simulation region used by local Galactic type of %gravity in grav_accel
   real                       :: ptmass                !< mass value of point %gravity source used by GRAV_PTMASS, GRAV_PTMASSSTIFF, GRAV_PTMASSPURE, GRAV_PTFLAT type of %gravity
   real                       :: ptm_x                 !< point mass position x-component
   real                       :: ptm_y                 !< point mass position y-component
   real                       :: ptm_z                 !< point mass position z-component
   real                       :: r_smooth              !< smoothing radius in point mass %types of %gravity
   integer(kind=4)            :: nsub                  !< number of subcells while additionally cell division in z-direction is present during establishment of hydrostatic equilibrium
   integer(kind=4)            :: n_gravr               !< index of hyperbolic-cosinusoidal cutting of gravitational potential used by GRAV_PTMASS, GRAV_PTFLAT type of %gravity
   integer(kind=4)            :: ord_pot2accel         !< Gradient operator order for gravitational potential
   real                       :: r_grav                !< radius of gravitational potential cut used by GRAV_PTMASS, GRAV_PTFLAT type of %gravity
   real                       :: tune_zeq              !< z-component of %gravity tuning factor used by hydrostatic_zeq
   real                       :: tune_zeq_bnd          !< z-component of %gravity tuning factor supposed to be used in boundaries
   real                       :: ptmass2               !< mass of the secondary for Roche potential
   real                       :: ptm2_x                !< x-position of the secondary for Roche potential (y and z positions are assumed to be 0)
   real                       :: cmass_x               !< center of mass for Roche potential
   real                       :: Omega                 !< corotational angular velocity for Roche potential

   logical                    :: user_grav             !< use user defined grav_pot_3d
   logical                    :: variable_gp           !< if .true. then cg%gp is evaluated at every step
   logical                    :: restart_gpot, restart_hgpot, restart_gp, restart_sgp, restart_sgpm !< if .true. then write this grav part to the restart files
   logical                    :: need_update           !< a flag to indicate that source_terms_grav needs to be called (e.g. because psolver skipped this)

   integer(kind=4)            :: ig_rk2_1, ig_rk2_2    !< indices to gravitational potential fields for different stages of the RK scheme
#ifdef SELF_GRAV
   integer(kind=4) :: i_sgp, i_sgpm
#endif /* SELF_GRAV */

   interface

      !< COMMENT ME
      subroutine user_grav_pot_3d
         implicit none
      end subroutine user_grav_pot_3d

      subroutine gprofs_default(iia, jja)
         implicit none
         integer, intent(in) :: iia                         !< COMMENT ME
         integer, intent(in) :: jja                         !< COMMENT ME
      end subroutine gprofs_default

      subroutine grav_types(gp, ax, lhn, flatten)
         use axes_M,    only: axes
         use constants, only: ndims, LO, HI
         implicit none
         real, dimension(:,:,:), pointer                     :: gp        !< gravitational potential array pointer
         type(axes),                              intent(in) :: ax        !< axes of cell centers positions
         integer(kind=4), dimension(ndims,LO:HI), intent(in) :: lhn       !< indices of gp array
         logical,                       optional, intent(in) :: flatten   !< optional parameter, commonly used to support 2D (alternative) potential version
      end subroutine grav_types

      subroutine user_grav_accel(sweep, i1,i2, xsw, n, grav)
         implicit none
         integer(kind=4),   intent(in)  :: sweep            !< COMMENT ME
         integer,           intent(in)  :: i1               !< COMMENT ME
         integer,           intent(in)  :: i2               !< COMMENT ME
         integer(kind=4),   intent(in)  :: n                !< COMMENT ME
         real, dimension(n),intent(in)  :: xsw              !< COMMENT ME
         real, dimension(n),intent(out) :: grav             !< COMMENT ME
      end subroutine user_grav_accel

      subroutine grav_pot2accel_T(sweep, i1, i2, n, grav, istep, cg)
         use grid_cont,          only: grid_container
         implicit none
         integer(kind=4),               intent(in)  :: sweep      !< string of characters that points out the current sweep direction
         integer,                       intent(in)  :: i1         !< number of column in the first direction after one pointed out by sweep
         integer,                       intent(in)  :: i2         !< number of column in the second direction after one pointed out by sweep
         integer(kind=4),               intent(in)  :: n          !< number of elements of returned array grav \todo OPT: would size(grav) be faster a bit?
         real, dimension(n),            intent(out) :: grav       !< 1D array of gravitational acceleration values computed for positions from %xsw and returned by the routine
         integer,                       intent(in)  :: istep      !< istep=RK2_1 for halfstep, istep=RK2_2 for fullstep in 2nd order Runge-Kutta method
         type(grid_container), pointer, intent(in)  :: cg         !< current grid_container
      end subroutine grav_pot2accel_T

   end interface

   procedure(user_grav_pot_3d), pointer :: grav_pot_3d => NULL()
   procedure(user_grav_accel),  pointer :: grav_accel  => NULL()
   procedure(gprofs_default),   pointer :: get_gprofs  => NULL()
   procedure(grav_types),       pointer :: grav_type   => NULL()
   procedure(grav_pot2accel_T), pointer :: grav_pot2accel => NULL()

contains

!>
!! \brief Routine that sets the initial values of %gravity parameters from namelist @c GRAVITY
!!
!! \n \n
!! @b GRAVITY
!! \n \n
!! <table border="+1">
!! <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!! <tr><td>g_dir        </td><td>0.0    </td><td>real             </td><td>\copydoc gravity::g_dir        </td></tr>
!! <tr><td>r_gc         </td><td>8500   </td><td>real             </td><td>\copydoc gravity::r_gc         </td></tr>
!! <tr><td>ptmass       </td><td>0.0    </td><td>non-negative real</td><td>\copydoc gravity::ptmass       </td></tr>
!! <tr><td>ptm_x        </td><td>0.0    </td><td>real             </td><td>\copydoc gravity::ptm_x        </td></tr>
!! <tr><td>ptm_y        </td><td>0.0    </td><td>real             </td><td>\copydoc gravity::ptm_y        </td></tr>
!! <tr><td>ptm_z        </td><td>0.0    </td><td>real             </td><td>\copydoc gravity::ptm_z        </td></tr>
!! <tr><td>r_smooth     </td><td>0.0    </td><td>real             </td><td>\copydoc gravity::r_smooth     </td></tr>
!! <tr><td>external_gp  </td><td>'null' </td><td>to be listed     </td><td>\copydoc gravity::external_gp  </td></tr>
!! <tr><td>nsub         </td><td>10     </td><td>integer > 0      </td><td>\copydoc gravity::nsub         </td></tr>
!! <tr><td>tune_zeq     </td><td>1.0    </td><td>real             </td><td>\copydoc gravity::tune_zeq     </td></tr>
!! <tr><td>tune_zeq_bnd </td><td>1.0    </td><td>real             </td><td>\copydoc gravity::tune_zeq_bnd </td></tr>
!! <tr><td>r_grav       </td><td>1.e6   </td><td>real             </td><td>\copydoc gravity::r_grav       </td></tr>
!! <tr><td>n_gravr      </td><td>0      </td><td>real             </td><td>\copydoc gravity::n_gravr      </td></tr>
!! <tr><td>user_grav    </td><td>.false.</td><td>logical          </td><td>\copydoc gravity::user_grav    </td></tr>
!! <tr><td>variable_gp  </td><td>.false.</td><td>logical          </td><td>\copydoc gravity::variable_gp  </td></tr>
!! <tr><td>restart_gp   </td><td>.false.</td><td>logical          </td><td>\copydoc gravity::restart_gp   </td></tr>
!! <tr><td>restart_gpot </td><td>.false.</td><td>logical          </td><td>\copydoc gravity::restart_gpot </td></tr>
!! <tr><td>restart_hgpot</td><td>.false.</td><td>logical          </td><td>\copydoc gravity::restart_hgpot</td></tr>
!! <tr><td>restart_sgp  </td><td>.false.</td><td>logical          </td><td>\copydoc gravity::restart_sgp  </td></tr>
!! <tr><td>restart_sgpm </td><td>.false.</td><td>logical          </td><td>\copydoc gravity::restart_sgpm </td></tr>
!! <tr><td>gprofs_target</td><td>'extgp'</td><td>string of chars  </td><td>\copydoc gravity::gprofs_target</td></tr>
!! </table>
!! The list is active while \b "GRAV" is defined.
!! \n \n
!<
   subroutine init_grav

      use cg_list_global,   only: all_cg
      use constants,        only: PIERNIK_INIT_MPI, gp_n, gpot_n, hgpot_n, O_I2, O_I4
      use dataio_pub,       only: printinfo, warn, die, code_progress, nh
      use mpisetup,         only: ibuff, rbuff, cbuff, master, slave, lbuff, piernik_MPI_Bcast
      use named_array_list, only: qna
      use units,            only: newtong
#ifdef SELF_GRAV
      use constants,        only: sgp_n, sgpm_n
#endif /* SELF_GRAV */
#ifdef NBODY
      use constants,        only: gp1b_n, nbdn_n, prth_n
#ifdef NBODY_GRIDDIRECT
      use constants,        only: nbgp_n
#endif /* NBODY_GRIDDIRECT */
#endif /* NBODY */
#ifdef CORIOLIS
      use coriolis,         only: set_omega
#endif /* CORIOLIS */

      implicit none

      namelist /GRAVITY/ g_dir, r_gc, ptmass, ptm_x, ptm_y, ptm_z, r_smooth, external_gp, ptmass2, ptm2_x, &
                         nsub, tune_zeq, tune_zeq_bnd, r_grav, n_gravr, user_grav, gprofs_target, variable_gp, ord_pot2accel, &
                         restart_gp, restart_gpot, restart_hgpot, restart_sgp, restart_sgpm

      if (code_progress < PIERNIK_INIT_MPI) call die("[gravity:init_grav] mpi not initialized.")

#ifdef VERBOSE
      if (master) call printinfo("[gravity:init_grav] Commencing gravity module initialization")
#endif /* VERBOSE */

      g_dir         = 0.0
      r_gc          = 8500
      ptmass        = 0.0
      ptm_x         = 0.0
      ptm_y         = 0.0
      ptm_z         = 0.0
      r_smooth      = 0.0
      nsub          = 10
      tune_zeq      = 1.0
      tune_zeq_bnd  = 1.0
      r_grav        = 1.e6
      ptmass2       = 0.0
      ptm2_x        = -1.0

      n_gravr       = 0
      ord_pot2accel = O_I2

      gprofs_target = 'extgp'
      external_gp   = 'null'

      user_grav     = .false.
      variable_gp   = .false.
      restart_gp    = .false.
      restart_gpot  = .false.
      restart_hgpot = .false.
      restart_sgp   = .false.
      restart_sgpm  = .false.

      if (master) then

         if (.not.nh%initialized) call nh%init()
         open(newunit=nh%lun, file=nh%tmp1, status="unknown")
         write(nh%lun,nml=GRAVITY)
         close(nh%lun)
         open(newunit=nh%lun, file=nh%par_file)
         nh%errstr=''
         read(unit=nh%lun, nml=GRAVITY, iostat=nh%ierrh, iomsg=nh%errstr)
         close(nh%lun)
         call nh%namelist_errh(nh%ierrh, "GRAVITY")
         read(nh%cmdl_nml,nml=GRAVITY, iostat=nh%ierrh)
         call nh%namelist_errh(nh%ierrh, "GRAVITY", .true.)
         open(newunit=nh%lun, file=nh%tmp2, status="unknown")
         write(nh%lun,nml=GRAVITY)
         close(nh%lun)
         call nh%compare_namelist()

         ibuff(1)   = nsub
         ibuff(2)   = n_gravr
         ibuff(3)   = ord_pot2accel

         rbuff(1:3) = g_dir
         rbuff(4)   = r_gc
         rbuff(5)   = ptmass
         rbuff(6)   = ptm_x
         rbuff(7)   = ptm_y
         rbuff(8)   = ptm_z
         rbuff(9)   = tune_zeq
         rbuff(10)  = tune_zeq_bnd
         rbuff(11)  = r_smooth
         rbuff(12)  = r_grav
         rbuff(13)  = ptmass2
         rbuff(14)  = ptm2_x

         lbuff(1)   = user_grav
         lbuff(2)   = variable_gp
         lbuff(3)   = restart_gp
         lbuff(4)   = restart_gpot
         lbuff(5)   = restart_hgpot
         lbuff(6)   = restart_sgp
         lbuff(7)   = restart_sgpm

         cbuff(1)   = gprofs_target
         cbuff(2)   = external_gp

      endif

      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(lbuff)
      call piernik_MPI_Bcast(rbuff)
      call piernik_MPI_Bcast(cbuff, cbuff_len)

      if (slave) then

         nsub          = int(ibuff(1), kind=4)
         n_gravr       = ibuff(2)
         ord_pot2accel = ibuff(3)

         g_dir         = rbuff(1:3)
         r_gc          = rbuff(4)
         ptmass        = rbuff(5)
         ptm_x         = rbuff(6)
         ptm_y         = rbuff(7)
         ptm_z         = rbuff(8)
         tune_zeq      = rbuff(9)
         tune_zeq_bnd  = rbuff(10)
         r_smooth      = rbuff(11)
         r_grav        = rbuff(12)
         ptmass2       = rbuff(13)
         ptm2_x        = rbuff(14)

         user_grav     = lbuff(1)
         variable_gp   = lbuff(2)
         restart_gp    = lbuff(3)
         restart_gpot  = lbuff(4)
         restart_hgpot = lbuff(5)
         restart_sgp   = lbuff(6)
         restart_sgpm  = lbuff(7)

         gprofs_target = cbuff(1)(1:gproft_len)
         external_gp   = cbuff(2)

      endif

      cmass_x = ptm_x
      Omega = 0.
      if (ptmass+ptmass2 > 0.) then
         cmass_x = (ptmass*ptm_x + ptmass2*ptm2_x)/(ptmass+ptmass2)
         if (abs(ptm_x-ptm2_x) > 0.) Omega = sqrt(newtong*(ptmass+ptmass2)/(abs(ptm_x-ptm2_x))**3)
#ifdef CORIOLIS
         call set_omega(Omega)
#endif /* CORIOLIS */
      endif

      ! Declare arrays for potential and make shortcuts
      ! All gravitational potential should be recalculated after refinement changes
      call all_cg%reg_var(gpot_n,  restart_mode = res_at(restart_gpot) )
      call all_cg%reg_var(hgpot_n, restart_mode = res_at(restart_hgpot))
      call all_cg%reg_var(gp_n,    restart_mode = res_at(restart_gp)   )
#ifdef SELF_GRAV
      call all_cg%reg_var(sgp_n,   restart_mode = res_at(restart_sgp)  )
      call all_cg%reg_var(sgpm_n,  restart_mode = res_at(restart_sgpm) )
      i_sgp  = qna%ind(sgp_n)
      i_sgpm = qna%ind(sgpm_n)
#endif /* SELF_GRAV */
#ifdef NBODY
      call all_cg%reg_var(prth_n)
      call all_cg%reg_var(nbdn_n)
      call all_cg%reg_var(gp1b_n)
#ifdef NBODY_GRIDDIRECT
      call all_cg%reg_var(nbgp_n)
#endif /* NBODY_GRIDDIRECT */
#endif /* NBODY */

      ig_rk2_1 = qna%ind(hgpot_n)
      ig_rk2_2 = qna%ind(gpot_n)

      if (.not.user_grav) then
         grav_pot_3d => default_grav_pot_3d
#ifdef VERBOSE
         if (master) call warn("[gravity:init_grav] user_grav is set to false. Using default grav_pot_3d.")
#endif /* VERBOSE */
      endif

      select case (ord_pot2accel)
         case (O_I2)
            grav_pot2accel => grav_pot2accel_ord2
         case (O_I4)
            grav_pot2accel => grav_pot2accel_ord4
         case default
            call die("[gravity:init_grav] Unknown gradient operator")
      end select

      call init_grav_ext

      need_update = .false.

   end subroutine init_grav

   integer(kind=4) function res_at(incl_gt)

      use constants, only: AT_IGNORE, AT_OUT_B

      implicit none

      logical, intent(in) :: incl_gt

      res_at = AT_IGNORE
      if (incl_gt) res_at = AT_OUT_B

   end function res_at

!> Register gravity-specific initialization of cg

   subroutine init_grav_ext

      use grid_container_ext, only: cg_ext, cg_extptrs

      implicit none

      procedure(cg_ext), pointer :: g_cg_init_p, g_cg_cleanup_p

      g_cg_init_p => g_cg_init
      g_cg_cleanup_p => null()
      call cg_extptrs%extend(g_cg_init_p, g_cg_cleanup_p, "gravity")

   end subroutine init_grav_ext

!> \brief Associate gravity-specific pointers

   subroutine g_cg_init(cg)

      use constants,        only: gp_n, gpot_n, hgpot_n, base_level_id
      use grid_cont,        only: grid_container
      use named_array_list, only: qna
#ifdef NBODY
      use constants,        only: gp1b_n, nbdn_n, prth_n
#ifdef NBODY_GRIDDIRECT
      use constants,        only: nbgp_n
#endif /* NBODY_GRIDDIRECT */
#endif /* NBODY */

      implicit none

      type(grid_container), pointer,  intent(inout) :: cg

      if (cg%l%id >= base_level_id) then

         cg%gpot  => cg%q(qna%ind( gpot_n))%arr
         cg%hgpot => cg%q(qna%ind(hgpot_n))%arr
         cg%gp    => cg%q(qna%ind(   gp_n))%arr
         cg%gp(:,:,:) = 0.0
         !> \todo move the following to multigrid?
#ifdef SELF_GRAV
         cg%sgp   => cg%q(i_sgp)%arr
         cg%sgpm  => cg%q(i_sgpm)%arr
#endif /* SELF_GRAV */
#ifdef NBODY
         cg%prth  => cg%q(qna%ind( prth_n))%arr
         cg%nbdn  => cg%q(qna%ind( nbdn_n))%arr
         cg%gp1b  => cg%q(qna%ind( gp1b_n))%arr
#ifdef NBODY_GRIDDIRECT
         cg%nbgp  => cg%q(qna%ind( nbgp_n))%arr
#endif /* NBODY_GRIDDIRECT */
#endif /* NBODY */

      endif

      !> \todo Place a call to initialize gp here, not in default_grav_pot_3d

   end subroutine g_cg_init

!>
!! \brief Collect gravitational terms depending on whether they are taken from restart file or not
!! \todo check if source_terms_grav should be called here for restarted simulation or new (non-restarted) simulation with particles.
!<
   subroutine init_terms_grav

      use dataio_pub, only: restarted_sim

      implicit none

      if (restarted_sim) then
         if (.not.restart_gp) call grav_pot_3d
         if (.not.restart_gpot) call compute_h_gpot
      else
         call update_gp
      endif

   end subroutine init_terms_grav

!> \brief Recover or update self-gravity potential, update other potential if needed

   subroutine source_terms_grav

      use constants,         only: PPP_GRAV
      use ppp,               only: ppp_main
#ifdef SELF_GRAV
      use cg_leaves,         only: leaves
      use cg_list_dataop,    only: expanded_domain
      use dataio_pub,        only: warn, die, restarted_sim
      use fluidindex,        only: iarr_all_sg
      use mpisetup,          only: master
      use multigrid_gravity, only: multigrid_solve_grav, recover_sgpm, recover_sgp
#endif /* SELF_GRAV */

      implicit none

      character(len=*), parameter :: grav_label = "source_terms_grav"
#ifdef SELF_GRAV
      logical, save :: frun = .true.
      logical :: sgpm_initialized
#endif /* SELF_GRAV */

      call ppp_main%start(grav_label, PPP_GRAV)

      need_update = .false.

#ifdef SELF_GRAV

      sgpm_initialized = .true.
      if (frun) then
         sgpm_initialized = recover_sgpm() ! try to recover sgpm from old soln
      else
         call leaves%q_copy(i_sgp, i_sgpm)
      endif

      if (frun .and. restarted_sim) then
         if (.not. recover_sgp()) call die("[gravity:source_terms_grav] cannot recover sgp")
         ! Reproducibility of restarts strongly depends on avalability of multigrid history in the restart file.
         ! ToDo: simplify the management of various histories of potential.
      else
         call multigrid_solve_grav(iarr_all_sg)
      endif
      frun = .false.

      call leaves%leaf_arr3d_boundaries(i_sgp) !, nocorners=.true.)
      ! No solvers should require corner values for the potential. Unfortunately some problems may relay on it indirectly (e.g. streaming_instability).
      !> \todo OPT: identify what relies on corner values of the potential and change it to work without corners. Then enable nocorners in the above call for some speedup.

      if (.not. sgpm_initialized) then
         call leaves%q_copy(i_sgp, i_sgpm) ! add fake history for selfgravitating potential: pretend that nothing was changing there until domain was created
         !> Restarted runs will be slightly affected as the previous selfgravitating potential was forgotten
         !> First step in highly dynamical setups will behave as the potential was frozen before first timestep
         !> Solution? Take one step backwards just for calculating old potential? Sounds complicated.
         !> Another solution: don't use extrapolation, exploit rich history instead and call multigrid more often.
         if (master) call warn("[gravity:source_terms_grav] assigned sgpm = sgp")
      endif

      call expanded_domain%q_copy(i_sgp, i_sgpm) ! add fake history for selfgravitating potential: pretend that nothing was changing there until domain expanded
#endif /* SELF_GRAV */
      if (variable_gp) call grav_pot_3d

      call ppp_main%stop(grav_label, PPP_GRAV)

   end subroutine source_terms_grav

!> \brief update static potential (gp field) in case of grid changes. Assume multigrid has been called

   subroutine update_gp

      implicit none

      if (associated(grav_pot_3d)) then
         call grav_pot_3d
         call compute_h_gpot
      endif

   end subroutine update_gp

!>
!! \brief Compute estimates of gravitational potential for t + dt/2 and t + dt
!! by extrapolation from current and previous potential. This is required by
!! RK2 integration scheme (RK4 can use it too).
!!
!! Please note that other high order integration schemes may need estimates for somewhat different time points.
!!
!! ToDo: In RK2 it is possible to obtain density estimate for t + dt  from the midpoint:
!!     u* = 2 u_h - u, which is equivalent to Euler's method
!! The gravitational potential from u* should offer better gpot and hgpot estimates than extrapolation.
!! The multigrid should be then called after execution of final RK2 step, but its cost should be
!! relatively small as we expect that less V-cycles would be required there.
!<

   subroutine compute_h_gpot

      use cg_leaves,        only: leaves
      use constants,        only: gp_n, gpot_n, hgpot_n
      use named_array_list, only: qna
#ifdef SELF_GRAV
      use constants,        only: one, half, zero
      use func,             only: operator(.notequals.)
      use global,           only: dt, dtm
#endif /* SELF_GRAV */
#if defined(SELF_GRAV) || defined(NBODY_GRIDDIRECT)
      use cg_list_dataop,   only: ind_val
#endif /* SELF_GRAV || NBODY_GRIDDIRECT */
#ifdef NBODY_GRIDDIRECT
      use constants,        only: nbgp_n
#endif /* NBODY_GRIDDIRECT */

      implicit none

#ifdef SELF_GRAV
      real :: h

      if (dtm .notequals. zero) then
         h = dt/dtm
      else
         h = 0.0
      endif

#ifdef NBODY_GRIDDIRECT
      !> \todo correct it
      call leaves%q_lin_comb([ ind_val(qna%ind(gp_n), 1.), ind_val(qna%ind(nbgp_n), 1.), ind_val(i_sgp, one+h),      ind_val(i_sgpm, -h)     ], qna%ind(gpot_n))
      call leaves%q_lin_comb([ ind_val(qna%ind(gp_n), 1.), ind_val(qna%ind(nbgp_n), 1.), ind_val(i_sgp, one+half*h), ind_val(i_sgpm, -half*h)], qna%ind(hgpot_n))
#else /* !NBODY_GRIDDIRECT */
      call leaves%q_lin_comb([ ind_val(qna%ind(gp_n), 1.), ind_val(i_sgp, one+h),      ind_val(i_sgpm, -h)     ], qna%ind(gpot_n))
      call leaves%q_lin_comb([ ind_val(qna%ind(gp_n), 1.), ind_val(i_sgp, one+half*h), ind_val(i_sgpm, -half*h)], qna%ind(hgpot_n))
#endif /* !NBODY_GRIDDIRECT */

#else /* !SELF_GRAV */
      !> \deprecated BEWARE: as long as grav_pot_3d is called only in init_piernik this assignment probably don't need to be repeated more than once
#ifdef NBODY_GRIDDIRECT
      !> \todo correct it
      call leaves%q_lin_comb([ ind_val(qna%ind(gp_n), 1.), ind_val(qna%ind(nbgp_n), 1.)], qna%ind(gpot_n))
      call leaves%q_lin_comb([ ind_val(qna%ind(gp_n), 1.), ind_val(qna%ind(nbgp_n), 1.)], qna%ind(hgpot_n))
#else /* !NBODY_GRIDDIRECT */
      call leaves%q_copy(qna%ind(gp_n), qna%ind(gpot_n))
      call leaves%q_copy(qna%ind(gp_n), qna%ind(hgpot_n))
#endif /* !NBODY_GRIDDIRECT */
#endif /* !SELF_GRAV */

   end subroutine compute_h_gpot

   subroutine grav_null(gp, ax, lhn, flatten)

      use axes_M,    only: axes
      use constants, only: ndims, LO, HI

      implicit none

      real, dimension(:,:,:), pointer                     :: gp
      type(axes),                              intent(in) :: ax
      integer(kind=4), dimension(ndims,LO:HI), intent(in) :: lhn
      logical,                       optional, intent(in) :: flatten

      gp = 0.0

      if (.false. .and. flatten) gp(1, 1, 1) = ax%x(1)*real(lhn(1,1)) ! suppress compiler warnings on unused arguments

   end subroutine grav_null

   subroutine grav_uniform(gp, ax, lhn, flatten)

      use axes_M,     only: axes
      use constants,  only: ndims, LO, HI, xdim, ydim, zdim, GEO_XYZ
      use dataio_pub, only: die
      use domain,     only: dom

      implicit none

      real, dimension(:,:,:), pointer                     :: gp
      type(axes),                              intent(in) :: ax
      integer(kind=4), dimension(ndims,LO:HI), intent(in) :: lhn
      logical,                       optional, intent(in) :: flatten
      integer                                             :: i, j, k

      if (dom%geometry_type /= GEO_XYZ) call die("[gravity:grav_uniform] Non-cartesian geometry is not implemented yet.")

      do i = lhn(xdim,LO), lhn(xdim,HI)
         do j = lhn(ydim,LO), lhn(ydim,HI)
            do k = lhn(zdim,LO), lhn(zdim,HI)
               gp(i,j,k) = -(g_dir(xdim)*ax%x(i) + g_dir(ydim)*ax%y(j) + g_dir(zdim)*ax%z(k))
            enddo
         enddo
      enddo

      if (.false. .and. flatten) i=0 ! suppress compiler warnings on unused arguments

   end subroutine grav_uniform

   subroutine grav_linear(gp, ax, lhn, flatten)

      use axes_M,     only: axes
      use constants,  only: ndims, LO, HI, xdim, ydim, zdim, half, GEO_XYZ
      use dataio_pub, only: die
      use domain,     only: dom

      implicit none

      real, dimension(:,:,:), pointer                     :: gp
      type(axes),                              intent(in) :: ax
      integer(kind=4), dimension(ndims,LO:HI), intent(in) :: lhn
      logical,                       optional, intent(in) :: flatten
      integer                                             :: i, j, k

      if (dom%geometry_type /= GEO_XYZ) call die("[gravity:grav_linear] Non-cartesian geometry is not implemented yet.")

      do i = lhn(xdim,LO), lhn(xdim,HI)
         do j = lhn(ydim,LO), lhn(ydim,HI)
            do k = lhn(zdim,LO), lhn(zdim,HI)
               gp(i,j,k) = -half*(g_dir(xdim)*ax%x(i)**2 + g_dir(ydim)*ax%y(j)**2 + g_dir(zdim)*ax%z(k)**2)
            enddo
         enddo
      enddo

      if (.false. .and. flatten) i=0 ! suppress compiler warnings on unused arguments

   end subroutine grav_linear

   subroutine grav_ptmass_pure(gp, ax, lhn, flatten)

      use axes_M,     only: axes
      use constants,  only: ndims, LO, HI, xdim, ydim, GEO_XYZ
      use dataio_pub, only: die
      use domain,     only: dom
      use units,     only: newtong

      implicit none

      real, dimension(:,:,:), pointer                     :: gp
      type(axes),                              intent(in) :: ax
      integer(kind=4), dimension(ndims,LO:HI), intent(in) :: lhn
      logical,                       optional, intent(in) :: flatten

      integer                                             :: i, j
      real                                                :: rc2, GM, x2
      logical                                             :: do_flatten

      if (dom%geometry_type /= GEO_XYZ) call die("[gravity:grav_ptmass_pure] Non-cartesian geometry is not implemented yet.")

      if (present(flatten)) then
         do_flatten = flatten
      else
         do_flatten = .false.
      endif

      GM = newtong*ptmass

      do i = lhn(xdim,LO), lhn(xdim,HI)
         x2 = (ax%x(i) - ptm_x)**2
         do j = lhn(ydim,LO), lhn(ydim,HI)
            rc2 = x2 + (ax%y(j) - ptm_y)**2

            if (do_flatten) then
               gp(i,j,:) = -GM / sqrt(rc2)
            else
               gp(i,j,:) = -GM / sqrt( (ax%z(:) - ptm_z)**2 + rc2 )
            endif

         enddo
      enddo

   end subroutine grav_ptmass_pure

   !>
   !! \todo this procedure is incompatible with cg%cs_iso2
   !<
   subroutine grav_ptmass_softened(gp, ax, lhn, flatten)

      use axes_M,     only: axes
      use constants,  only: ndims, LO, HI, xdim, ydim, GEO_XYZ
      use dataio_pub, only: die
      use domain,     only: dom
      use fluidindex, only: flind
      use global,     only: smalld
      use units,      only: newtong

      implicit none

      real, dimension(:,:,:), pointer                     :: gp
      type(axes),                              intent(in) :: ax
      integer(kind=4), dimension(ndims,LO:HI), intent(in) :: lhn
      logical,                       optional, intent(in) :: flatten
      integer                                             :: i, j, ifl
      real                                                :: rc2, r_smooth2, GM, fr, x2, cs_iso2
      logical                                             :: do_flatten

      if (dom%geometry_type /= GEO_XYZ) call die("[gravity:grav_ptmass_softened] Non-cartesian geometry is not implemented yet.")

      if (present(flatten)) then
         do_flatten = flatten
      else
         do_flatten = .false.
      endif

      !cs_iso2   = maxval(flind%all_fluids(:)%fl%cs2)
      cs_iso2 = 0.0
      do ifl = lbound(flind%all_fluids, dim=1), ubound(flind%all_fluids, dim=1)
         cs_iso2 = max(cs_iso2, flind%all_fluids(ifl)%fl%cs2)
      enddo

      r_smooth2 = r_smooth**2
      GM        = newtong*ptmass

      do i = lhn(xdim,LO), lhn(xdim,HI)
         x2 = (ax%x(i) - ptm_x)**2
         do j = lhn(ydim,LO), lhn(ydim,HI)
            rc2 = x2 + (ax%y(j) - ptm_y)**2
            fr  = min( (sqrt(rc2)/r_grav)**n_gravr , 100.0)    !> \deprecated BEWARE: hardcoded value
            fr  = max( 1./cosh(fr), smalld*1.e-2)              !> \deprecated BEWARE: hadrcoded value
            fr  = -cs_iso2 * log(fr)

            if (do_flatten) then
               gp(i,j,:) = -GM / sqrt( rc2 + r_smooth2 ) + fr
            else
               gp(i,j,:) = -GM / sqrt( (ax%z(:) - ptm_z)**2 + rc2 + r_smooth2 ) + fr
            endif

         enddo
      enddo

   end subroutine grav_ptmass_softened

!>
!! \brief Roche potential for two bodies on circular orbits. The coordinate system corotates with the bodies, so they stay on the X axis forever.
!<

   subroutine grav_roche(gp, ax, lhn, flatten)

      use axes_M,     only: axes
      use constants,  only: ndims, LO, HI, ydim, zdim, half, GEO_XYZ
      use dataio_pub, only: die
      use domain,     only: dom
      use units,      only: newtong

      implicit none

      real, dimension(:,:,:), pointer                     :: gp
      type(axes),                              intent(in) :: ax
      integer(kind=4), dimension(ndims,LO:HI), intent(in) :: lhn
      logical,                       optional, intent(in) :: flatten
      integer                                             :: j, k
      real                                                :: GM1, GM2, z2, yz2, r_smooth2

      if (dom%geometry_type /= GEO_XYZ) call die("[gravity:grav_roche] Non-cartesian geometry is not implemented yet.")

      r_smooth2 = r_smooth**2
      GM1 =  newtong * ptmass
      GM2 =  newtong * ptmass2

      do k = lhn(zdim,LO), lhn(zdim,HI)
         z2 = ax%z(k)**2
         do j = lhn(ydim,LO), lhn(ydim,HI)
            yz2 = ax%y(j)**2 + z2
            gp(:,j,k) =  - GM1 / sqrt((ax%x(:) - ptm_x)**2  + yz2 + r_smooth2) &
                 &       - GM2 / sqrt((ax%x(:) - ptm2_x)**2 + yz2 + r_smooth2) &
                 &       - half * Omega**2 * ((ax%x(:) - cmass_x)**2 + yz2)
         enddo
      enddo

      if (.false. .and. flatten) j=0 ! suppress compiler warnings on unused arguments

   end subroutine grav_roche

!>
!! \details promote stiff-body rotation inside smoothing length, don't affect the global potential outside
!<
   subroutine grav_ptmass_stiff(gp, ax, lhn, flatten)

      use axes_M,     only: axes
      use constants,  only: ndims, LO, HI, xdim, ydim, zdim, half, GEO_XYZ, GEO_RPZ
      use dataio_pub, only: die
      use domain,     only: dom
      use units,      only: newtong

      implicit none

      real, dimension(:,:,:), pointer                     :: gp
      type(axes),                              intent(in) :: ax
      integer(kind=4), dimension(ndims,LO:HI), intent(in) :: lhn
      logical,                       optional, intent(in) :: flatten

      integer                                             :: j, k
      real                                                :: r_smooth2, gmr, gm, z2, yz2, ptx_cos

      r_smooth2 = r_smooth**2
      gm =  - newtong * ptmass
      gmr = half * gm / r_smooth

      do k = lhn(zdim,LO), lhn(zdim,HI)
         z2 = (ax%z(k) - ptm_z)**2
         select case (dom%geometry_type)
            case (GEO_XYZ)
               do j = lhn(ydim,LO), lhn(ydim,HI)
                  yz2 = z2 + (ax%y(j) - ptm_y)**2
                  gp(lhn(xdim,LO):lhn(xdim,HI), j, k) = potstiff(yz2 + (ax%x(lhn(xdim,LO):lhn(xdim,HI)) - ptm_x)**2)
               enddo
            case (GEO_RPZ)
               do j = lhn(ydim,LO), lhn(ydim,HI)
                  ptx_cos = 2*ptm_x*cos(ax%y(j) - ptm_y)
                  gp(lhn(xdim,LO):lhn(xdim,HI), j, k) = potstiff(z2 + ax%x(lhn(xdim,LO):lhn(xdim,HI))**2 + ptm_x**2 - ax%x(lhn(xdim,LO):lhn(xdim,HI))*ptx_cos)
               enddo
            case default
               call die("[gravity:grav_ptmass_stiff] Non-cartesian geometry is not implemented yet.")
         end select
      enddo

      if (.false. .and. flatten) j=0 ! suppress compiler warnings on unused arguments

   contains

      real elemental function potstiff(r2)

         implicit none

         real, intent(in) :: r2

         if (r2 < r_smooth2) then
            potstiff = gmr * (3. - r2/r_smooth2)
         else
            potstiff = gm / sqrt(r2)
         endif

      end function potstiff

   end subroutine grav_ptmass_stiff

!--------------------------------------------------------------------------
!>
!! \brief Routine that compute values of gravitational potential filling in gp array and setting gp_status character string \n\n
!! The type of %gravity is governed by external_gp value: \n\n
!! \details
!! \b GRAV_NULL, \b grav_null, \b null - gravitational potential array is set to zero \n\n
!! \b GRAV_UNIFORM, \b grav_unif, \b uniform - uniform type of %gravity in z-direction \n
!! \f$\Phi\left(z\right)= - const \cdot z \f$\n
!! where \f$ const \f$ is set by parameter @c g_z \n\n
!! \b GRAV_LINEAR, \b grav_lin, \b linear - linear type of %gravity growing along z-direction \n
!! \f$\Phi\left(z\right) = -1/2 \cdot const \cdot z^2\f$ \n\n
!! \b GRAV_PTMASS, \b ptmass_soft, \b softened \b ptmass - softened point mass type of %gravity \n
!! \f$\Phi\left(x,y,z\right)= - GM/\sqrt{x^2+y^2+z^2+r_{soft}^2}\f$ \n
!! where \f$r_{soft}\f$ is a radius of softening\n\n
!! \b GRAV_PTMASSPURE, \b ptmass_pure, \b ptmass - unsoftened point mass type of %gravity \n
!! \f$\Phi\left(x,y,z\right)= - GM/\sqrt{x^2+y^2+z^2}\f$ \n\n
!! \b GRAV_PTMASSSTIFF, \b ptmass_stiff, \b stiff \b ptmass - softened point mass type of %gravity with stiff-body rotation inside softening radius\n
!! \f$\Phi\left(x,y,z\right)= - GM/\sqrt{x^2+y^2+z^2}\f$ for \f$r > r_{soft}\f$ and \f$ GM/r_{soft} \left( - 3/2 + 1/2 {x^2+y^2+z^2}/r_{soft}^2 \right)\f$ inside \f$r_{soft}\f$ \n\n
!! \b GRAV_PTFLAT, \b flat_ptmass_soft, \b flat \b softened \b ptmass - planar, softened point mass type of %gravity \n
!! \f$\Phi\left(x,y,z\right)= - GM/\sqrt{x^2+y^2+r_{soft}^2}\f$ \n
!! where \f$r_{soft}\f$ is a radius of softening\n\n
!! \b flat_ptmass, \b flat \b ptmass - planar, pure point mass type of %gravity \n
!! \f$\Phi\left(x,y\right)= - GM/\sqrt{x^2+y^2}\f$ \n\n
!! where \f$r_{soft}\f$ is a radius of softening\n\n
!! \b GRAV_ROCHE, \b grav_roche, \b roche - Roche potential for two bodies on circular orbits. The coordinate system corotates with the bodies, so they stay on the X axis forever. \n\n
!! \b GRAV_USER, \b grav_user, \b user - not a standard type of %gravity, implemented by user in the routine grav_pot_user from gravity_user module.\n\n
!! If none of them is specified grav_accel specified by user is called.
!<

   subroutine default_grav_pot_3d

      use axes_M,     only: axes
      use cg_leaves,  only: leaves
      use cg_list,    only: cg_list_element
      use constants,  only: GEO_XYZ
      use dataio_pub, only: die, warn
      use domain,     only: dom
      use grid_cont,  only: grid_container
      use mpisetup,   only: master

      implicit none

      type(axes)                     :: ax
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         if (.not. cg%is_old) then

            call ax%allocate_axes(cg%lhn)
            ax%x(:) = cg%x(:)
            ax%y(:) = cg%y(:)
            ax%z(:) = cg%z(:)

            gp_status = ''

            if (dom%geometry_type /= GEO_XYZ) then
               select case (external_gp)
                  case ("null", "grav_null", "GRAV_NULL")
                     ! No gravity - no problem, selfgravity has to check the geometry during initialization
                  case ("user", "grav_user", "GRAV_USER")
                     ! The User knows what he/she is doing ...
                  case default ! standard cases do not support cylindrical geometry yet
                     if (master) call warn("[gravity:default_grav_pot_3d] Non-cartesian geometry may or may not be implemented correctly.")
               end select
            endif

            select case (external_gp)
               case ("null", "grav_null", "GRAV_NULL")
                  call grav_null(cg%gp, ax, cg%lhn)                     ; grav_type => grav_null
               case ("linear", "grav_lin", "GRAV_LINEAR")
                  call grav_linear(cg%gp, ax, cg%lhn)                   ; grav_type => grav_linear
               case ("uniform", "grav_unif", "GRAV_UNIFORM")
                  call grav_uniform(cg%gp, ax, cg%lhn)                  ; grav_type => grav_uniform
               case ("softened ptmass", "ptmass_soft", "GRAV_PTMASS")
                  call grav_ptmass_softened(cg%gp, ax, cg%lhn, .false.) ; grav_type => grav_ptmass_softened
               case ("stiff ptmass", "ptmass_stiff", "GRAV_PTMASSSTIFF")
                  call grav_ptmass_stiff(cg%gp, ax, cg%lhn)             ; grav_type => grav_ptmass_stiff
               case ("ptmass", "ptmass_pure", "GRAV_PTMASSPURE")
                  call grav_ptmass_pure(cg%gp, ax, cg%lhn, .false.)     ; grav_type => grav_ptmass_pure
               case ("flat softened ptmass", "flat_ptmass_soft", "GRAV_PTFLAT")
                  call grav_ptmass_softened(cg%gp, ax, cg%lhn, .true.)  ; grav_type => grav_ptmass_softened
               case ("flat ptmass", "flat_ptmass")
                  call grav_ptmass_pure(cg%gp, ax, cg%lhn, .true.)      ; grav_type => grav_ptmass_pure
               case ("roche", "grav_roche", "GRAV_ROCHE")
#ifndef CORIOLIS
                  call die("[gravity:default_grav_pot_3d] define CORIOLIS in piernik.def for Roche potential")
#endif /* !CORIOLIS */
                  call grav_roche(cg%gp, ax, cg%lhn)                    ; grav_type => grav_roche
               case ("user", "grav_user", "GRAV_USER")
                  call die("[gravity:default_grav_pot_3d] user 'grav_pot_3d' should be defined in initprob!")
               case default
                  gp_status = 'undefined'
                  call die("[gravity:default_grav_pot_3d] unimplemented external_gp: '" // trim(external_gp) // "'")
            end select

            call ax%deallocate_axes

         endif
         cgl => cgl%nxt
      enddo

   end subroutine default_grav_pot_3d

!>
!! \brief Routine that compute values of gravitational acceleration using gravitational potential array gp
!! \todo offer high order gradient as an option in parameter file
!<
   subroutine grav_pot2accel_ord2(sweep, i1, i2, n, grav, istep, cg)

      use constants,        only: idm2, ndims, pdims, ydim, half, LO, HI, GEO_XYZ, GEO_RPZ, RK2_1, RK2_2, EULER, NORMAL, ORTHO1, ORTHO2, INVALID
      use dataio_pub,       only: die
      use domain,           only: dom
      use grid_cont,        only: grid_container

      implicit none

      integer(kind=4),               intent(in)  :: sweep      !< string of characters that points out the current sweep direction
      integer,                       intent(in)  :: i1         !< number of column in the first direction after one pointed out by sweep
      integer,                       intent(in)  :: i2         !< number of column in the second direction after one pointed out by sweep
      integer(kind=4),               intent(in)  :: n          !< number of elements of returned array grav \todo OPT: would size(grav) be faster a bit?
      real, dimension(n),            intent(out) :: grav       !< 1D array of gravitational acceleration values computed for positions from %xsw and returned by the routine
      integer,                       intent(in)  :: istep      !< istep=RK2_1 for halfstep, istep=RK2_2 for fullstep in 2nd order Runge-Kutta method
      type(grid_container), pointer, intent(in)  :: cg         !< current grid_container

      integer, dimension(ndims,LO:HI)            :: ispan
      integer(kind=4)                            :: ig

      ! Gravitational acceleration is computed on right cell boundaries

      ! For more general schemes (higher order than RK2, non-canonical choices of RK2) we will need timestep fraction provided by the solver, not istep
      select case (istep)
      case (RK2_1)
         ig = ig_rk2_1
      case (RK2_2, EULER)
         ig = ig_rk2_2
      case default
         call die("[gravity:grav_pot2accel_ord2] Unsupported substep")
         ig = INVALID  ! suppress compiler warning
      end select

      ispan(pdims(sweep,ORTHO1),:) = i1
      ispan(pdims(sweep,ORTHO2),:) = i2
      ispan(pdims(sweep,NORMAL),LO:HI) = cg%lhn(sweep,LO:HI) + [2,0]

      grav(2:n-1) = half*reshape(cg%q(ig)%span(ispan-2*idm2(sweep,:,:)) - cg%q(ig)%span(ispan),[n-2]) / cg%dl(sweep)
      grav(1) = grav(2); grav(n) = grav(n-1)

      select case (dom%geometry_type)
         case (GEO_XYZ) ! Do nothing
         case (GEO_RPZ)
            if (sweep == ydim) grav = grav * cg%inv_x(i2)
         case default
            call die("[gravity:grav_pot2accel_ord2] Unsupported geometry")
      end select

   end subroutine grav_pot2accel_ord2

   subroutine grav_pot2accel_ord4(sweep, i1, i2, n, grav, istep, cg)

      use constants,        only: idm2, ndims, pdims, ydim, LO, HI, GEO_XYZ, GEO_RPZ, RK2_1, RK2_2, EULER, NORMAL, ORTHO1, ORTHO2, INVALID
      use dataio_pub,       only: die
      use domain,           only: dom
      use grid_cont,        only: grid_container

      implicit none

      integer(kind=4),               intent(in)  :: sweep      !< string of characters that points out the current sweep direction
      integer,                       intent(in)  :: i1         !< number of column in the first direction after one pointed out by sweep
      integer,                       intent(in)  :: i2         !< number of column in the second direction after one pointed out by sweep
      integer(kind=4),               intent(in)  :: n          !< number of elements of returned array grav \todo OPT: would size(grav) be faster a bit?
      real, dimension(n),            intent(out) :: grav       !< 1D array of gravitational acceleration values computed for positions from %xsw and returned by the routine
      integer,                       intent(in)  :: istep      !< istep=RK2_1 for halfstep, istep=RK2_2 for fullstep in 2nd order Runge-Kutta method
      type(grid_container), pointer, intent(in)  :: cg         !< current grid_container

      integer, dimension(ndims,LO:HI)            :: ispan
      integer(kind=4)                            :: ig
      real, parameter                            :: onetw = 1./12.

      ! Gravitational acceleration is computed on right cell boundaries

      select case (istep)
      case (RK2_1)
         ig = ig_rk2_1
      case (RK2_2, EULER)
         ig = ig_rk2_2
      case default
         call die("[gravity:grav_pot2accel_ord4] Unsupported substep")
         ig = INVALID  ! suppress compiler warning
      end select

      ispan(pdims(sweep,ORTHO1),:) = i1
      ispan(pdims(sweep,ORTHO2),:) = i2
      ispan(pdims(sweep,NORMAL),LO:HI) = cg%lhn(sweep,LO:HI) + [2,-2]

      grav(3:n-2) = onetw*reshape(cg%q(ig)%span(ispan+2*idm2(sweep,:,:)) - 8.*cg%q(ig)%span(ispan+  idm2(sweep,:,:)) + &
                               8.*cg%q(ig)%span(ispan-  idm2(sweep,:,:)) -    cg%q(ig)%span(ispan-2*idm2(sweep,:,:)),[n-4]) / cg%dl(sweep)
      grav(2) = grav(3); grav(n-1) = grav(n-2)
      grav(1) = grav(2); grav(n) = grav(n-1)

      select case (dom%geometry_type)
         case (GEO_XYZ) ! Do nothing
         case (GEO_RPZ)
            if (sweep == ydim) grav = grav * cg%inv_x(i2)
         case default
            call die("[gravity:grav_pot2accel_ord4] Unsupported geometry")
      end select

   end subroutine grav_pot2accel_ord4

!--------------------------------------------------------------------------
!>
!! \brief Routine that is an interface between grav_pot2accel and internal_sources
!<
   subroutine grav_src_exec(n, u, cg, sweep, i1, i2, istep, gravsrc)

      use fluidindex, only: flind, iarr_all_dn, iarr_all_mx
#ifndef ISO
      use fluidindex, only: iarr_all_en
#endif /* !ISO */
      use grid_cont,  only: grid_container

      implicit none

      integer(kind=4),               intent(in)  :: n          !< array size
      real, dimension(n, flind%all), intent(in)  :: u          !< vector of conservative variables
      type(grid_container), pointer, intent(in)  :: cg         !< current grid_container
      integer(kind=4),               intent(in)  :: sweep      !< string of characters that points out the current sweep direction
      integer,                       intent(in)  :: i1         !< number of column in the first direction after one pointed out by sweep
      integer,                       intent(in)  :: i2         !< number of column in the second direction after one pointed out by sweep
      integer,                       intent(in)  :: istep      !< istep=RK2_1 for halfstep, istep=RK2_2 for fullstep in 2nd order Runge-Kutta method
      real, dimension(n, flind%all), intent(out) :: gravsrc    !< u array update from sources
      real, dimension(n)                         :: gravacc    !< acceleration caused by gravitation
      logical                                    :: full_dim

      full_dim = n > 1
      gravsrc = 0.0
      if (.not.full_dim) return

      call grav_pot2accel(sweep, i1, i2, n, gravacc, istep, cg)
      gravsrc(:, iarr_all_mx) = spread(gravacc,2,flind%fluids) * u(:, iarr_all_dn)
#ifndef ISO
      gravsrc(:, iarr_all_en) = spread(gravacc,2,flind%fluids) * u(:, iarr_all_mx)
#endif /* !ISO */

   end subroutine grav_src_exec

end module gravity
