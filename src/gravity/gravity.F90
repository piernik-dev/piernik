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
!! \brief (DW) Module containing all main subroutines and functions that govern %gravity force in the code
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
   public :: init_grav, grav_accel, source_terms_grav, grav_pot2accel, grav_pot_3d, grav_pot_3d_called, grav_type, get_gprofs, grav_accel2pot, sum_potential, grav_pot_3d_bnd
   public :: g_dir, r_gc, ptmass, ptm_x, ptm_y, ptm_z, r_smooth, nsub, tune_zeq, tune_zeq_bnd, h_grav, r_grav, n_gravr, n_gravh, user_grav, gp_status, gprofs_target, ptmass2, ptm2_x

   integer, parameter         :: gp_stat_len   = 9
   integer, parameter         :: gproft_len    = 5
   character(len=gp_stat_len) :: gp_status       !< variable set as 'undefined' in grav_pot_3d when grav_accel is supposed to use
   character(len=gproft_len)  :: gprofs_target   !< variable set pointing gravity routine in hydrostatic_zeq ('accel' or ready gp array 'gparr')
   character(len=cbuff_len)   :: external_gp     !< variable allowing to choose external gravitational potential
   real, dimension(ndims)     :: g_dir           !< vector used by GRAV_UNIFORM and GRAV_LINEAR type of %gravity
   real    :: r_gc                  !< galactocentric radius of the local simulation region used by local Galactic type of %gravity in grav_accel
   real    :: ptmass                !< mass value of point %gravity source used by GRAV_PTMASS, GRAV_PTMASSSTIFF, GRAV_PTMASSPURE, GRAV_PTFLAT type of %gravity
   real    :: ptm_x                 !< point mass position x-component
   real    :: ptm_y                 !< point mass position y-component
   real    :: ptm_z                 !< point mass position z-component
   real    :: r_smooth              !< smoothing radius in point mass %types of %gravity
   integer(kind=4) :: nsub          !< number of subcells while additionally cell division in z-direction is present during establishment of hydrostatic equilibrium
   real    :: h_grav                !< altitude of acceleration cut used when n_gravh is set to non-zero
   real    :: r_grav                !< radius of gravitational potential cut used by GRAV_PTMASS, GRAV_PTFLAT type of %gravity
   integer :: n_gravh               !< index of hyperbolic-cosinusoidal cutting of acceleration; used when set to non-zero
   integer :: n_gravr               !< index of hyperbolic-cosinusoidal cutting of gravitational potential used by GRAV_PTMASS, GRAV_PTFLAT type of %gravity
   real    :: tune_zeq              !< z-component of %gravity tuning factor used by hydrostatic_zeq
   real    :: tune_zeq_bnd          !< z-component of %gravity tuning factor supposed to be used in boundaries
   real    :: ptmass2               !< mass of the secondary for Roche potential
   real    :: ptm2_x                !< x-position of the secondary for Roche potential (y and z positions are assumed to be 0)
   real    :: cmass_x               !< center of mass for Roche potential
   real    :: Omega                 !< corotational angular velocity for Roche potential

   logical :: grav_pot_3d_called = .false.
   logical :: user_grav             !< use user defined grav_pot_3d
   logical :: variable_gp           !< .true. if arrays::gp must be evaluated at every step

   interface

      !< COMMENT ME
      subroutine user_grav_pot_3d
         implicit none
      end subroutine user_grav_pot_3d

      subroutine gprofs_default(iia, jja)
         implicit none
         integer, intent(in) :: iia                    !< COMMENT ME
         integer, intent(in) :: jja                    !< COMMENT ME
      end subroutine gprofs_default

      subroutine grav_types(gp,ax,flatten)
         use types, only: axes
         implicit none
         real, dimension(:,:,:), pointer       :: gp        !< COMMENT ME
         type(axes),             intent(in)    :: ax        !< COMMENT ME
         logical,      optional, intent(in)    :: flatten   !< COMMENT ME
      end subroutine grav_types

      subroutine user_grav_accel(sweep, i1,i2, xsw, n, grav)
         implicit none
         integer(kind=4), intent(in)    :: sweep            !< COMMENT ME
         integer, intent(in)            :: i1               !< COMMENT ME
         integer, intent(in)            :: i2               !< COMMENT ME
         integer(kind=4), intent(in)    :: n                !< COMMENT ME
         real, dimension(n),intent(in)  :: xsw              !< COMMENT ME
         real, dimension(n),intent(out) :: grav             !< COMMENT ME
      end subroutine user_grav_accel

   end interface

   procedure(user_grav_pot_3d), pointer :: grav_pot_3d => NULL()
   procedure(user_grav_accel),  pointer :: grav_accel  => NULL()
   procedure(gprofs_default),   pointer :: get_gprofs  => NULL()
   procedure(grav_types),       pointer :: grav_type   => NULL()

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
!! <tr><td>h_grav       </td><td>1.e6   </td><td>real             </td><td>\copydoc gravity::h_grav       </td></tr>
!! <tr><td>r_grav       </td><td>1.e6   </td><td>real             </td><td>\copydoc gravity::r_grav       </td></tr>
!! <tr><td>n_gravr      </td><td>0      </td><td>real             </td><td>\copydoc gravity::n_gravr      </td></tr>
!! <tr><td>n_gravh      </td><td>0      </td><td>real             </td><td>\copydoc gravity::n_gravh      </td></tr>
!! <tr><td>user_grav    </td><td>.false.</td><td>logical          </td><td>\copydoc gravity::user_grav    </td></tr>
!! <tr><td>variable_gp  </td><td>.false.</td><td>logical          </td><td>\copydoc gravity::variable_gp  </td></tr>
!! <tr><td>gprofs_target</td><td>'gparr'</td><td>string of chars  </td><td>\copydoc gravity::gprofs_target</td></tr>
!! </table>
!! \n \n
!<
   subroutine init_grav

      use dataio_pub,    only: ierrh, par_file, namelist_errh, compare_namelist, cmdl_nml, lun, getlun    ! QA_WARN required for diff_nml
      use dataio_pub,    only: printinfo, warn, die, code_progress
      use constants,     only: PIERNIK_INIT_GRID, AT_OUT_B, AT_IGNORE, gp_n, gpot_n, hgpot_n
      use mpisetup,      only: ibuff, rbuff, cbuff, comm, ierr, master, slave, lbuff, buffer_dim, FIRST
      use mpi,           only: MPI_DOUBLE_PRECISION, MPI_INTEGER, MPI_LOGICAL, MPI_CHARACTER
      use units,         only: newtong
      use grid,          only: all_cg
      use gc_list,       only: cg_list_element
#ifdef SELF_GRAV
      use constants,     only: sgp_n, sgpm_n
#endif /* SELF_GRAV */
#ifdef CORIOLIS
      use coriolis,      only: set_omega
#endif /* CORIOLIS */

      implicit none

      type(cg_list_element), pointer :: cgl

      namelist /GRAVITY/ g_dir, r_gc, ptmass, ptm_x, ptm_y, ptm_z, r_smooth, external_gp, ptmass2, ptm2_x, &
                nsub, tune_zeq, tune_zeq_bnd, h_grav, r_grav, n_gravr, n_gravh, user_grav, gprofs_target, variable_gp

      if (code_progress < PIERNIK_INIT_GRID) call die("[gravity:init_grav] units or arrays not initialized.")

#ifdef VERBOSE
      if (master) call printinfo("[gravity:init_grav] Commencing gravity module initialization")
#endif /* VERBOSE */

      g_dir   = 0.0
      r_gc    = 8500
      ptmass  = 0.0
      ptm_x   = 0.0
      ptm_y   = 0.0
      ptm_z   = 0.0
      r_smooth= 0.0
      nsub    = 10
      tune_zeq     = 1.0
      tune_zeq_bnd = 1.0
      h_grav = 1.e6
      r_grav = 1.e6
      ptmass2 = 0.0
      ptm2_x = -1.0

      n_gravr = 0
      n_gravh = 0

      gprofs_target = 'gparr'
      external_gp   = 'null'

      user_grav = .false.
      variable_gp = .false.

      if (master) then

         diff_nml(GRAVITY)

         ibuff(1)   = nsub
         ibuff(2)   = n_gravr
         ibuff(3)   = n_gravh

         rbuff(1:3) = g_dir
         rbuff(4)   = r_gc
         rbuff(5)   = ptmass
         rbuff(6)   = ptm_x
         rbuff(7)   = ptm_y
         rbuff(8)   = ptm_z
         rbuff(9)  = tune_zeq
         rbuff(10)  = tune_zeq_bnd
         rbuff(11)  = r_smooth
         rbuff(12)  = h_grav
         rbuff(13)  = r_grav
         rbuff(14)  = ptmass2
         rbuff(15)  = ptm2_x

         lbuff(1)   = user_grav
         lbuff(2)   = variable_gp

         cbuff(1)   = gprofs_target
         cbuff(2)   = external_gp

      endif

      call MPI_Bcast(ibuff,           buffer_dim, MPI_INTEGER,          FIRST, comm, ierr)
      call MPI_Bcast(lbuff,           buffer_dim, MPI_LOGICAL,          FIRST, comm, ierr)
      call MPI_Bcast(rbuff,           buffer_dim, MPI_DOUBLE_PRECISION, FIRST, comm, ierr)
      call MPI_Bcast(cbuff, cbuff_len*buffer_dim, MPI_CHARACTER,        FIRST, comm, ierr)

      if (slave) then

         nsub                = int(ibuff(1), kind=4)
         n_gravr             = ibuff(2)
         n_gravh             = ibuff(3)

         g_dir               = rbuff(1:3)
         r_gc                = rbuff(4)
         ptmass              = rbuff(5)
         ptm_x               = rbuff(6)
         ptm_y               = rbuff(7)
         ptm_z               = rbuff(8)
         tune_zeq            = rbuff(9)
         tune_zeq_bnd        = rbuff(10)
         r_smooth            = rbuff(11)
         h_grav              = rbuff(12)
         r_grav              = rbuff(13)
         ptmass2             = rbuff(14)
         ptm2_x              = rbuff(15)

         user_grav           = lbuff(1)
         variable_gp         = lbuff(2)

         gprofs_target       = cbuff(1)(1:gproft_len)
         external_gp         = cbuff(2)

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
      call all_cg%reg_var(gpot_n, AT_IGNORE)
      call all_cg%reg_var(hgpot_n, AT_IGNORE)
      call all_cg%reg_var(gp_n, AT_OUT_B)
#ifdef SELF_GRAV
      call all_cg%reg_var(sgp_n, AT_IGNORE)
      call all_cg%reg_var(sgpm_n, AT_IGNORE)
#endif /* SELF_GRAV */

      cgl => all_cg%first
      do while (associated(cgl))
         cgl%cg%gpot => cgl%cg%get_na_ptr(gpot_n)
         cgl%cg%gpot(:,:,:) = 0.0
         cgl%cg%hgpot => cgl%cg%get_na_ptr(hgpot_n)
         cgl%cg%gp => cgl%cg%get_na_ptr(gp_n)
#ifdef SELF_GRAV
         cgl%cg%sgp => cgl%cg%get_na_ptr(sgp_n)
         cgl%cg%sgpm => cgl%cg%get_na_ptr(sgpm_n)
#endif /* SELF_GRAV */
         cgl => cgl%nxt
      enddo

      if (.not.user_grav) then
         grav_pot_3d => default_grav_pot_3d
#ifdef VERBOSE
         if (master) call warn("[gravity:init_grav] user_grav is set to false. Using default grav_pot_3d.")
#endif /* VERBOSE */
      endif

   end subroutine init_grav

   subroutine source_terms_grav

#ifdef SELF_GRAV
      use constants,         only: sgp_n
      use dataio_pub,        only: die
      use domain,            only: is_multicg
      use fluidindex,        only: iarr_all_sg
      use gc_list,           only: cg_list_element
      use grid,              only: all_cg
      use grid_cont,         only: grid_container
      use external_bnd,      only: arr3d_boundaries
#ifdef POISSON_FFT
      use poissonsolver,     only: poisson_solve
#endif /* POISSON_FFT */
#ifdef MULTIGRID
      use multigrid_gravity, only: multigrid_solve_grav
#endif /* MULTIGRID */
#endif /* SELF_GRAV */

      implicit none

#ifdef SELF_GRAV
      logical, save :: frun = .true.
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg

      if (is_multicg) call die("[gravity:source_terms_grav] multiple grid pieces per procesor not implemented yet") !nontrivial all cg% must be solved at a time (nontrivial for multigrid, rarely possible dor FFT poisson solver)

      cgl => all_cg%first
      do while (associated(cgl))
         cg => cgl%cg
         cg%sgpm = cg%sgp

#ifdef POISSON_FFT
         call poisson_solve( sum(cg%u(iarr_all_sg,:,:,:),1) )
#endif /* POISSON_FFT */
#ifdef MULTIGRID
         if (size(iarr_all_sg) == 1) then
            call multigrid_solve_grav(cg%u(iarr_all_sg(1),:,:,:))
         else
            call multigrid_solve_grav( sum(cg%u(iarr_all_sg,:,:,:),1) )
            !>
         !! \deprecated BEWARE Here a lot of heap space is required and some compilers may generate code that do segfaults for big enough domains.
            !! It is the weakest point of this type in Maclaurin test. Next one (in fluidboundaries.F90) is 8 times less sensitive.
            !<
         endif
#endif /* MULTIGRID */
         cgl => cgl%nxt
      enddo

      ! communicate boundary values for sgp(:, :, :) because multigrid solver gives at most 2 guardcells, while for hydro solver typically 4 is required.

!> \warning An improper evaluation of guardcell potential may occur when the multigrid boundary conditions doesn't match /BOUNDARIES/ namelist (e.g. isolated on periodic domain).
      call arr3d_boundaries(all_cg%first%cg%get_na_ind(sgp_n))

      if (frun) then
         cgl => all_cg%first
         do while (associated(cgl))
            cgl%cg%sgpm = cgl%cg%sgp
            cgl => cgl%nxt
         enddo
         frun = .false.
      endif
#endif /* SELF_GRAV */
      if (variable_gp) then
         call grav_pot_3d
         call grav_pot_3d_bnd
      endif

      call sum_potential

   end subroutine source_terms_grav

   subroutine grav_pot_3d_bnd
      use constants,         only: xdim, ydim, zdim, LO, HI, BND_OUT, BND_OUTH
      use grid,              only: all_cg
      use gc_list,           only: cg_list_element
      use grid_cont,         only: grid_container
      use domain,            only: dom
      implicit none
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer  :: cg
      integer :: i

      cgl => all_cg%first
      do while (associated(cgl))
         cg => cgl%cg

         if (dom%has_dir(xdim)) then
            if (cg%bnd(xdim,LO) >= BND_OUT .and. cg%bnd(xdim,LO) <= BND_OUTH) then
               do i = 1, cg%nb+1
                  cg%gp(i,:,:)               = cg%gp(cg%nb+2,:,:)
               enddo
            endif

            if (cg%bnd(xdim,HI) >= BND_OUT .and. cg%bnd(xdim,HI) <= BND_OUTH) then
               do i = 1, cg%nb+1
                  cg%gp(cg%n_(xdim)-cg%nb-1+i,:,:) = cg%gp(cg%n_(xdim)-cg%nb-1,:,:)
               enddo
            endif
         endif

         if (dom%has_dir(ydim)) then
            if (cg%bnd(ydim,LO) >= BND_OUT .and. cg%bnd(ydim,LO) <= BND_OUTH) then
               do i = 1, cg%nb+1
                  cg%gp(:,i,:)               = cg%gp(:,cg%nb+2,:)
               enddo
            endif

            if (cg%bnd(ydim,HI) >= BND_OUT .and. cg%bnd(ydim,HI) <= BND_OUTH) then
               do i = 1, cg%nb+1
                  cg%gp(:,cg%n_(ydim)-cg%nb-1+i,:) = cg%gp(:,cg%n_(ydim)-cg%nb-1,:)
               enddo
            endif
         endif

         if (dom%has_dir(zdim)) then
            if (cg%bnd(zdim,LO) >= BND_OUT .and. cg%bnd(zdim,LO) <= BND_OUTH) then
               do i = 1, cg%nb+1
                  cg%gp(:,:,i)               = cg%gp(:,:,cg%nb+2)
               enddo
            endif

            if (cg%bnd(zdim,HI) >= BND_OUT .and. cg%bnd(zdim,HI) <= BND_OUTH) then
               do i = 1, cg%nb+1
                  cg%gp(:,:,cg%n_(zdim)-cg%nb-1+i) = cg%gp(:,:,cg%n_(zdim)-cg%nb-1)
               enddo
            endif
         endif

         cgl => cgl%nxt
      enddo

   end subroutine grav_pot_3d_bnd

   subroutine sum_potential

      use constants, only: one, half
      use global,    only: dt, dtm
      use grid,      only: all_cg
      use gc_list,   only: cg_list_element
      use grid_cont, only: grid_container

      implicit none

      real :: h
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg

      if (dtm /= 0) then
         h = dt/dtm
      else
         h = 0.0
      endif

      cgl => all_cg%first
      do while (associated(cgl))
         cg => cgl%cg
#ifdef SELF_GRAV
         cg%gpot  = cg%gp + (one+h)     *cg%sgp -      h*cg%sgpm
         cg%hgpot = cg%gp + (one+half*h)*cg%sgp - half*h*cg%sgpm
#else /* !SELF_GRAV */
         !> \deprecated BEWARE: as long as grav_pot_3d is called only in init_piernik this assignment probably don't need to be repeated more than once
         cg%gpot  = cg%gp
         cg%hgpot = cg%gp
#endif /* !SELF_GRAV */
         cgl => cgl%nxt
      enddo

   end subroutine sum_potential

   subroutine grav_null(gp, ax, flatten)

      use types, only: axes

      implicit none

      real, dimension(:,:,:), pointer       :: gp
      type(axes), intent(in)                :: ax
      logical, intent(in), optional         :: flatten

      gp = 0.0

      if (.false. .and. flatten) gp(1, 1, 1) = ax%x(1) ! suppress compiler warnings on unused arguments

   end subroutine grav_null

   subroutine grav_uniform(gp, ax, flatten)

      use types, only: axes

      implicit none

      real, dimension(:,:,:), pointer       :: gp
      type(axes), intent(in)                :: ax
      logical, intent(in), optional         :: flatten
      integer :: i, j, k

      do i = 1, ubound(gp,1)
         do j = 1, ubound(gp,2)
            do k = 1, ubound(gp,3)
               gp(i,j,k) = -(g_dir(1)*ax%x(i) + g_dir(2)*ax%y(j) + g_dir(3)*ax%z(k))
            enddo
         enddo
      enddo

      if (.false. .and. flatten) i=0 ! suppress compiler warnings on unused arguments

   end subroutine grav_uniform

   subroutine grav_linear(gp, ax, flatten)

      use constants, only: half
      use types,     only: axes

      implicit none

      real, dimension(:,:,:), pointer       :: gp
      type(axes), intent(in)                :: ax
      logical, intent(in), optional         :: flatten
      integer :: i, j, k

      do i = 1, ubound(gp,1)
         do j = 1, ubound(gp,2)
            do k = 1, ubound(gp,3)
               gp(i,j,k) = -half*(g_dir(1)*ax%x(i)**2 + g_dir(2)*ax%y(j)**2 + g_dir(3)*ax%z(k)**2)
            enddo
         enddo
      enddo

      if (.false. .and. flatten) i=0 ! suppress compiler warnings on unused arguments

   end subroutine grav_linear

   subroutine grav_ptmass_pure(gp, ax, flatten)

      use units, only: newtong
      use types, only: axes

      implicit none

      real, dimension(:,:,:), pointer       :: gp
      type(axes), intent(in)                :: ax
      logical, intent(in), optional         :: flatten

      integer             :: i, j
      real                :: rc2, GM, x2
      logical             :: do_flatten

      if (present(flatten)) then
         do_flatten = flatten
      else
         do_flatten = .false.
      endif

      GM        = newtong*ptmass

      do i = 1, ubound(gp,1)
         x2 = (ax%x(i) - ptm_x)**2
         do j = 1, ubound(gp,2)
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
   subroutine grav_ptmass_softened(gp, ax, flatten)

      use fluidindex, only: flind
      use global,     only: smalld
      use types,      only: axes
      use units,      only: newtong

      implicit none

      real, dimension(:,:,:), pointer       :: gp
      type(axes), intent(in)                :: ax
      logical, intent(in), optional         :: flatten

      integer             :: i, j
      real                :: rc2, r_smooth2, GM, fr, x2
      real                :: cs_iso2
      logical             :: do_flatten

      if (present(flatten)) then
         do_flatten = flatten
      else
         do_flatten = .false.
      endif

      cs_iso2   = maxval(flind%all_fluids(:)%cs2)
      r_smooth2 = r_smooth**2
      GM        = newtong*ptmass

      do i = 1, ubound(gp,1)
         x2 = (ax%x(i) - ptm_x)**2
         do j = 1, ubound(gp,2)
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
!! \brief Roche potential for two bodies oncircular orbits. The coordinate system corotates with the bodies, so they stay on the X axis forever.
!<

   subroutine grav_roche(gp, ax, flatten)

      use constants, only: ydim, zdim, half
      use types,     only: axes
      use units,     only: newtong

      implicit none

      real, dimension(:,:,:), pointer       :: gp
      type(axes), intent(in)                :: ax
      logical, intent(in), optional         :: flatten

      integer :: j, k
      real    :: r_smooth2
      real    :: GM1, GM2, z2, yz2

      r_smooth2 = r_smooth**2
      GM1 =  newtong * ptmass
      GM2 =  newtong * ptmass2

      do k = 1, ubound(gp,zdim)
         z2 = ax%z(k)**2
         do j = 1, ubound(gp,ydim)
            yz2 = ax%y(j)**2 + z2
            gp(:,j,k) =  - GM1 / sqrt((ax%x(:) - ptm_x)**2  + yz2 + r_smooth2) &
                 &       - GM2 / sqrt((ax%x(:) - ptm2_x)**2 + yz2 + r_smooth2) &
                 &       - half * Omega**2 * ((ax%x(:) - cmass_x)**2 + yz2)
         enddo
      enddo

      if (.false. .and. flatten) j=0 ! suppress compiler warnings on unused arguments

   end subroutine grav_roche

   subroutine grav_ptmass_stiff(gp, ax, flatten)

      use constants, only: half
      use types,     only: axes
      use units,     only: newtong

      implicit none

      real, dimension(:,:,:), pointer       :: gp
      type(axes), intent(in)                :: ax
      logical, intent(in), optional         :: flatten

      integer :: i, j, k
      real    :: r_smooth2, r2, gmr, gm, z2, yz2
      ! promote stiff-body rotation inside smoothing length, don't affect the global potential outside

      r_smooth2 = r_smooth**2
      gm =  - newtong * ptmass
      gmr = half * gm / r_smooth

      do k = 1, ubound(gp,3)
         z2 = (ax%z(k) - ptm_z)**2
         do j = 1, ubound(gp,2)
            yz2 = z2 + (ax%y(j) - ptm_y)**2
            do i = 1, ubound(gp,1)
               r2 = yz2 + (ax%x(i) - ptm_x)**2
               if (r2 < r_smooth2) then
                  gp(i,j,k) = gmr * (3. - r2/r_smooth2)
               else
                  gp(i,j,k) = gm / sqrt(r2)
               endif
            enddo
         enddo
      enddo

      if (.false. .and. flatten) i=0 ! suppress compiler warnings on unused arguments

   end subroutine grav_ptmass_stiff

!--------------------------------------------------------------------------
!>
!! \brief Routine that compute values of gravitational potential filling in gp array and setting gp_status character string \n\n
!! The type of %gravity is governed by preprocessor: \n\n
!! \deprecated BEWARE: This is no longer true: external_gp does the magic
!! \details
!! GRAV_NULL - gravitational potential array is set to zero \n\n
!! GRAV_UNIFORM - uniform type of %gravity in z-direction \n
!! \f$\Phi\left(z\right)= - const \cdot z \f$\n
!! where \f$ const \f$ is set by parameter @c g_z \n\n
!! GRAV_LINEAR - linear type of %gravity growing along z-direction \n
!! \f$\Phi\left(z\right) = -1/2 \cdot const \cdot z^2\f$ \n\n
!! GRAV_PTMASS - softened point mass type of %gravity \n
!! \f$\Phi\left(x,y,z\right)= - GM/\sqrt{x^2+y^2+z^2+r_{soft}^2}\f$ \n
!! where \f$r_{soft}\f$ is a radius of softening\n\n
!! GRAV_PTMASSPURE - unsoftened point mass type of %gravity \n
!! \f$\Phi\left(x,y,z\right)= - GM/\sqrt{x^2+y^2+z^2}\f$ \n\n
!! GRAV_PTMASSSTIFF - softened point mass type of %gravity with stiff-body rotation inside softening radius\n
!! \f$\Phi\left(x,y,z\right)= - GM/\sqrt{x^2+y^2+z^2}\f$ for \f$r > r_{soft}\f$ and \f$ GM/r_{soft} \left( - 3/2 + 1/2 {x^2+y^2+z^2}/r_{soft}^2 \right)\f$ inside \f$r_{soft}\f$ \n\n
!! GRAV_PTFLAT - planar, softened point mass type of %gravity \n
!! \f$\Phi\left(x,y,z\right)= - GM/\sqrt{x^2+y^2+r_{soft}^2}\f$ \n
!! where \f$r_{soft}\f$ is a radius of softening\n\n
!! GRAV_USER - not a standard type of %gravity, implemented by user in the routine grav_pot_user from gravity_user module.\n\n
!<

   subroutine default_grav_pot_3d

      use constants,  only: GEO_XYZ
      use dataio_pub, only: die, warn
      use domain,     only: dom
      use grid,       only: all_cg
      use gc_list,    only: cg_list_element
      use grid_cont,  only: grid_container
      use mpisetup,   only: master
      use types,      only: axes

      implicit none

      type(axes) :: ax
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg

      cgl => all_cg%first
      do while (associated(cgl))
         cg => cgl%cg
         if (.not.allocated(ax%x)) allocate(ax%x(size(cg%x)))
         if (.not.allocated(ax%y)) allocate(ax%y(size(cg%y)))
         if (.not.allocated(ax%z)) allocate(ax%z(size(cg%z)))
         ax%x = cg%x
         ax%y = cg%y
         ax%z = cg%z

         gp_status = ''

         if (dom%geometry_type /= GEO_XYZ) then
            select case (external_gp)
               case ("null", "grav_null", "GRAV_NULL")
                  ! No gravity - no problem, selfgravity has to check the geometry during initialization
               case ("user", "grav_user", "GRAV_USER")
                  ! The User knows what he/she is doing ...
               case default ! standard cases do not support cylindrical geometry yet
                  call die("[gravity:default_grav_pot_3d] Non-cartesian geometry is not implemented.")
            end select
         endif

         select case (external_gp)
            case ("null", "grav_null", "GRAV_NULL")
               call grav_null(cg%gp,ax)                    ; grav_type => grav_null
            case ("linear", "grav_lin", "GRAV_LINEAR")
               call grav_linear(cg%gp,ax)                  ; grav_type => grav_linear
            case ("uniform", "grav_unif", "GRAV_UNIFORM")
               call grav_uniform(cg%gp,ax)                 ; grav_type => grav_uniform
            case ("softened ptmass", "ptmass_soft", "GRAV_PTMASS")
               call grav_ptmass_softened(cg%gp,ax,.false.) ; grav_type => grav_ptmass_softened
            case ("stiff ptmass", "ptmass_stiff", "GRAV_PTMASSSTIFF")
               call grav_ptmass_stiff(cg%gp,ax)            ; grav_type => grav_ptmass_stiff
            case ("ptmass", "ptmass_pure", "GRAV_PTMASSPURE")
               call grav_ptmass_pure(cg%gp,ax,.false.)     ; grav_type => grav_ptmass_pure
            case ("flat softened ptmass", "flat_ptmass_soft", "GRAV_PTFLAT")
               call grav_ptmass_softened(cg%gp,ax,.true.)  ; grav_type => grav_ptmass_softened
            case ("flat ptmass", "flat_ptmass")
               call grav_ptmass_pure(cg%gp,ax,.true.)      ; grav_type => grav_ptmass_pure
            case ("roche", "grav_roche", "GRAV_ROCHE")
#ifndef CORIOLIS
               call die("[gravity:default_grav_pot_3d] define CORIOLIS in piernik.def for Roche potential")
#endif /* !CORIOLIS */
               call grav_roche(cg%gp,ax)                   ; grav_type => grav_roche
            case ("user", "grav_user", "GRAV_USER")
               call die("[gravity:default_grav_pot_3d] user 'grav_pot_3d' should be defined in initprob!")
            case default
               gp_status = 'undefined'
         end select

!-----------------------

         if (gp_status == 'undefined') then
            if (associated(grav_accel)) then
               if (master) call warn("[gravity:default_grav_pot_3d]: using 'grav_accel' defined by user")
               call grav_accel2pot
            else
               call die("[gravity:default_grav_pot_3d]: GRAV is defined, but 'gp' is not initialized")
            endif
         endif

         if (allocated(ax%x)) deallocate(ax%x)
         if (allocated(ax%y)) deallocate(ax%y)
         if (allocated(ax%z)) deallocate(ax%z)

         cgl => cgl%nxt
      enddo

   end subroutine default_grav_pot_3d

!>
!! \brief Routine that compute values of gravitational acceleration using gravitational potential array gp
!<
   subroutine grav_pot2accel(sweep, i1, i2, n, grav, istep, cg)

      use constants, only: xdim, ydim, zdim, half
      use grid_cont, only: grid_container

      implicit none

      integer(kind=4), intent(in)    :: sweep      !< string of characters that points out the current sweep direction
      integer, intent(in)            :: i1         !< number of column in the first direction after one pointed out by sweep
      integer, intent(in)            :: i2         !< number of column in the second direction after one pointed out by sweep
      integer(kind=4), intent(in)    :: n          !< number of elements of returned array grav
      real, dimension(n),intent(out) :: grav       !< 1D array of gravitational acceleration values computed for positions from %xsw and returned by the routine
      integer, intent(in)            :: istep      !< istep=1 for halfstep, istep=2 for fullstep
      type(grid_container), pointer, intent(in) :: cg

!> \todo offer high order gradient as an option in parameter file
!      real, parameter :: onetw = 1./12.

! Gravitational acceleration is computed on right cell boundaries

!      if (istep==1) then
!         select case (sweep)
!            case (xdim)
!               grav(3:n-2) = onetw*(cg%hgpot(5:n,i1,i2) - 8.*cg%hgpot(4:n-1,i1,i2) + 8.*cg%hgpot(2:n-3,i1,i2) - cg%hgpot(1:n-4,i1,i2) )/dl(xdim)
!            case (ydim)
!               grav(3:n-2) = onetw*(cg%hgpot(i2,5:n,i1) - 8.*cg%hgpot(i2,4:n-1,i1) + 8.*cg%hgpot(i2,2:n-3,i1) - cg%hgpot(i2,1:n-4,i1) )/dl(xdim)
!            case (zdim)
!               grav(3:n-2) = onetw*(cg%hgpot(i1,i2,5:n) - 8.*cg%hgpot(i1,i2,4:n-1) + 8.*cg%hgpot(i1,i2,2:n-3) - cg%hgpot(i1,i2,1:n-4) )/dl(xdim)
!         end select
!      else
!         select case (sweep)
!            case (xdim)
!               grav(3:n-2) = onetw*(cg%gpot(5:n,i1,i2) - 8.*cg%gpot(4:n-1,i1,i2) + 8.*cg%gpot(2:n-3,i1,i2) - cg%gpot(1:n-4,i1,i2) )/dl(xdim)
!            case (ydim)
!               grav(3:n-2) = onetw*(cg%gpot(i2,5:n,i1) - 8.*cg%gpot(i2,4:n-1,i1) + 8.*cg%gpot(i2,2:n-3,i1) - cg%gpot(i2,1:n-4,i1) )/dl(xdim)
!            case (zdim)
!               grav(3:n-2) = onetw*(cg%gpot(i1,i2,5:n) - 8.*cg%gpot(i1,i2,4:n-1) + 8.*cg%gpot(i1,i2,2:n-3) - cg%gpot(i1,i2,1:n-4) )/dl(xdim)
!         end select
!      endif
!      grav(2) = grav(3); grav(n-1) = grav(n-2)
!      grav(1) = grav(2); grav(n) = grav(n-1)
      if (istep==1) then
         select case (sweep)
            case (xdim)
               grav(2:n-1) = half*(cg%hgpot(1:n-2,i1,i2) - cg%hgpot(3:n,i1,i2))/cg%dl(xdim)
            case (ydim)
               grav(2:n-1) = half*(cg%hgpot(i2,1:n-2,i1) - cg%hgpot(i2,3:n,i1))/cg%dl(ydim)
            case (zdim)
               grav(2:n-1) = half*(cg%hgpot(i1,i2,1:n-2) - cg%hgpot(i1,i2,3:n))/cg%dl(zdim)
         end select

      else
         select case (sweep)
            case (xdim)
               grav(2:n-1) = half*(cg%gpot(1:n-2,i1,i2) - cg%gpot(3:n,i1,i2))/cg%dl(xdim)
            case (ydim)
               grav(2:n-1) = half*(cg%gpot(i2,1:n-2,i1) - cg%gpot(i2,3:n,i1))/cg%dl(ydim)
            case (zdim)
               grav(2:n-1) = half*(cg%gpot(i1,i2,1:n-2) - cg%gpot(i1,i2,3:n))/cg%dl(zdim)
         end select
      endif

      grav(1) = grav(2); grav(n) = grav(n-1)
   end subroutine grav_pot2accel

!--------------------------------------------------------------------------
!>
!! \brief Routine that uses %gravity acceleration given in grav_accel to compute values of gravitational potential filling in gp array
!!
!! \details Gravitational potential gp(i,j,k) is defined in cell centers
!! First we find gp(:,:,:) in each block separately.
!! Instead of gp, gpwork is used to not change already computed possibly existing other parts of gravitational potential.
!<
   subroutine grav_accel2pot

      use constants,  only: xdim, ydim, zdim, ndims, MAXL, I_ONE
      use dataio_pub, only: die
      use domain,     only: is_mpi_noncart, is_multicg, cdd, dom
      use gc_list,    only: get_extremum
      use grid,       only: all_cg
      use grid_cont,  only: grid_container !, cg_list_element
      use mpi,        only: MPI_DOUBLE_PRECISION, MPI_COMM_NULL
      use mpisetup,   only: master, nproc, FIRST, LAST, comm, ierr, have_mpi
      use types,      only: value

      implicit none

      integer                                                          :: i, j, k, ip, px, py, pz
      integer, dimension(3)                                            :: pc
      real, allocatable, dimension(:,:,:), target                      :: gpwork
      real, dimension(:,:,:), pointer                                  :: p
      real, dimension(:), allocatable                                  :: gravrx, gravry, gravrz
      real                                                             :: dgpx_proc, dgpy_proc, dgpz_proc, ddgph
      real, dimension(FIRST:LAST)                                      :: dgpx_all,  dgpy_all,  dgpz_all
      real, dimension(0:cdd%psize(xdim)-1,0:cdd%psize(ydim)-1,0:cdd%psize(zdim)-1) :: dgpx,      dgpy,      dgpz,     ddgp
      type(value)                                                      :: gp_max
      type(grid_container), pointer :: cg

      cg => all_cg%first%cg
      if (is_multicg) call die("[gravity:grav_accel2pot] multiple grid pieces per procesor not implemented yet") !nontrivial

      if (any([allocated(gravrx), allocated(gravry), allocated(gravrz)])) call die("[gravity:grav_accel2pot] gravr[xyz] already allocated")
      allocate(gravrx(cg%n_(xdim)), gravry(cg%n_(ydim)), gravrz(cg%n_(zdim)))

      if (have_mpi .and. is_mpi_noncart) call die("[gravity:grav_accel2pot] is_mpi_noncart is not implemented") ! MPI_Cart_coords, psize, pcoords
      if (cdd%comm3d == MPI_COMM_NULL) call die("[gravity:grav_accel2pot] cdd%comm3d == MPI_COMM_NULL")

      allocate(gpwork(cg%n_(xdim), cg%n_(ydim), cg%n_(zdim)))
      gpwork(1,1,1) = 0.0

      call grav_accel(xdim, 1, 1, cg%xr(:), cg%n_(xdim), gravrx)
      do i = 1, cg%n_(xdim)-1
         gpwork(i+1,1,1) = gpwork(i,1,1) - gravrx(i)*cg%dl(xdim)
      enddo

      do i=1, cg%n_(xdim)
         call grav_accel(ydim, 1, i, cg%yr(:), cg%n_(ydim), gravry)
         do j = 1, cg%n_(ydim)-1
            gpwork(i,j+1,1) = gpwork(i,j,1) - gravry(j)*cg%dl(ydim)
         enddo
      enddo

      do i=1, cg%n_(xdim)
         do j=1, cg%n_(ydim)
            call grav_accel(zdim, i, j, cg%zr(:), cg%n_(zdim), gravrz)
            do k = 1, cg%n_(zdim)-1
               gpwork(i,j,k+1) = gpwork(i,j,k) - gravrz(k)*cg%dl(zdim)
            enddo
         enddo
      enddo

      dgpx_proc = gpwork(cg%ie+dom%D_x, cg%js,         cg%ks        )-gpwork(cg%is,cg%js,cg%ks)
      dgpy_proc = gpwork(cg%is,         cg%je+dom%D_y, cg%ks        )-gpwork(cg%is,cg%js,cg%ks)
      dgpz_proc = gpwork(cg%is,         cg%js,         cg%ke+dom%D_z)-gpwork(cg%is,cg%js,cg%ks)

      call MPI_Gather ( dgpx_proc, I_ONE, MPI_DOUBLE_PRECISION, dgpx_all, I_ONE, MPI_DOUBLE_PRECISION, FIRST, cdd%comm3d, ierr )
      call MPI_Gather ( dgpy_proc, I_ONE, MPI_DOUBLE_PRECISION, dgpy_all, I_ONE, MPI_DOUBLE_PRECISION, FIRST, cdd%comm3d, ierr )
      call MPI_Gather ( dgpz_proc, I_ONE, MPI_DOUBLE_PRECISION, dgpz_all, I_ONE, MPI_DOUBLE_PRECISION, FIRST, cdd%comm3d, ierr )

      if (master) then

         do ip = FIRST, LAST
            call MPI_Cart_coords(cdd%comm3d, ip, ndims, pc, ierr)

            px = pc(1)
            py = pc(2)
            pz = pc(3)

            dgpx(px,py,pz) = dgpx_all(ip)
            dgpy(px,py,pz) = dgpy_all(ip)
            dgpz(px,py,pz) = dgpz_all(ip)

         enddo

         ddgp(0,0,0) = 0.0
         do i = 1, cdd%psize(xdim)-1
            ddgp(i,0,0) = ddgp(i-1,0,0) + dgpx(i-1,0,0)
         enddo

         do i=0, cdd%psize(xdim)-1
            do j = 1, cdd%psize(ydim)-1
               ddgp(i,j,0) = ddgp(i,j-1,0) + dgpy(i,j-1,0)
            enddo
         enddo

         do i=0, cdd%psize(xdim)-1
            do j=0, cdd%psize(ydim)-1
               do k = 1, cdd%psize(zdim)-1
                  ddgp(i,j,k)= ddgp(i,j,k-1) + dgpz(i,j,k-1)
               enddo
            enddo
         enddo

      endif

      call MPI_Bcast(ddgp, nproc, MPI_DOUBLE_PRECISION, FIRST, comm, ierr)

      px = cdd%pcoords(xdim)
      py = cdd%pcoords(ydim)
      pz = cdd%pcoords(zdim)

      ddgph  = gpwork(1,1,1)-gpwork(cg%is,cg%js,cg%ks)
      gpwork = gpwork + ddgp(px,py,pz) + ddgph
      p => gpwork(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke)
      call get_extremum(p, MAXL, gp_max, cg)

      call MPI_Bcast(gp_max%val, I_ONE, MPI_DOUBLE_PRECISION, gp_max%proc, comm, ierr)
      gpwork = gpwork - gp_max%val

      cg%gp = gpwork
      if (allocated(gpwork)) deallocate(gpwork)

      deallocate(gravrx)
      deallocate(gravry)
      deallocate(gravrz)

   end subroutine grav_accel2pot

end module gravity
