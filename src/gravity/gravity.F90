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

   use constants, only: cbuff_len

   implicit none

   private
   public :: init_grav, grav_accel, source_terms_grav, grav_pot2accel, grav_pot_3d, grav_pot_3d_called, grav_type, get_gprofs, grav_accel2pot, sum_potential
   public :: g_dir, r_gc, ptmass, ptm_x, ptm_y, ptm_z, r_smooth, nsub, tune_zeq, tune_zeq_bnd, h_grav, r_grav, n_gravr, n_gravh, user_grav, gp_status, gprofs_target, ptmass2, ptm2_x

   integer, parameter         :: gp_stat_len   = 9
   integer, parameter         :: gproft_len    = 5
   character(len=gp_stat_len) :: gp_status       !< variable set as 'undefined' in grav_pot_3d when grav_accel is supposed to use
   character(len=gproft_len)  :: gprofs_target   !< variable set pointing gravity routine in hydrostatic_zeq ('accel' or ready gp array 'gparr')
   character(len=cbuff_len)   :: external_gp     !< variable allowing to choose external gravitational potential
   real, dimension(3)         :: g_dir           !< vector used by GRAV_UNIFORM and GRAV_LINEAR type of %gravity
   real    :: r_gc                  !< galactocentric radius of the local simulation region used by local Galactic type of %gravity in grav_accel
   real    :: ptmass                !< mass value of point %gravity source used by GRAV_PTMASS, GRAV_PTMASSSTIFF, GRAV_PTMASSPURE, GRAV_PTFLAT type of %gravity
   real    :: ptm_x                 !< point mass position x-component
   real    :: ptm_y                 !< point mass position y-component
   real    :: ptm_z                 !< point mass position z-component
   real    :: r_smooth              !< smoothing radius in point mass %types of %gravity
   integer :: nsub                  !< number of subcells while additionally cell division in z-direction is present during establishment of hydrostatic equilibrium
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

      subroutine gprofs_default(iia,jja)
         implicit none
         integer, intent(in) :: iia                    !< COMMENT ME
         integer, intent(in) :: jja                    !< COMMENT ME
      end subroutine gprofs_default

      subroutine grav_types(gp,ax,flatten)
         use types, only: axes
         implicit none
         real, dimension(:,:,:), intent(inout) :: gp        !< COMMENT ME
         type(axes),             intent(in)    :: ax        !< COMMENT ME
         logical,      optional, intent(in)    :: flatten   !< COMMENT ME
      end subroutine grav_types

      subroutine user_grav_accel(sweep, i1,i2, xsw, n, grav)
         implicit none
         integer, intent(in)            :: sweep            !< COMMENT ME
         integer, intent(in)            :: i1               !< COMMENT ME
         integer, intent(in)            :: i2               !< COMMENT ME
         integer, intent(in)            :: n                !< COMMENT ME
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

      use arrays,        only: gpot
      use dataio_pub,    only: ierrh, par_file, namelist_errh, compare_namelist, cmdl_nml    ! QA_WARN required for diff_nml
      use dataio_pub,    only: warn, die, code_progress
      use constants,     only: PIERNIK_INIT_ARRAYS
      use mpisetup,      only: ibuff, rbuff, cbuff, comm, ierr, master, slave, lbuff, buffer_dim
      use mpi,           only: MPI_DOUBLE_PRECISION, MPI_INTEGER, MPI_LOGICAL, MPI_CHARACTER
      use units,     only: newtong
#ifdef CORIOLIS
      use coriolis,      only: set_omega
#endif /* CORIOLIS */

      implicit none

      namelist /GRAVITY/ g_dir, r_gc, ptmass, ptm_x, ptm_y, ptm_z, r_smooth, external_gp, ptmass2, ptm2_x, &
                nsub, tune_zeq, tune_zeq_bnd, h_grav, r_grav, n_gravr, n_gravh, user_grav, gprofs_target, variable_gp

      if (code_progress < PIERNIK_INIT_ARRAYS) call die("[gravity:init_grav] units or arrays not initialized.")

#ifdef VERBOSE
      if (master) call warn("[gravity:init_grav] Commencing gravity module initialization")
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

      call MPI_Bcast(ibuff,           buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_Bcast(lbuff,           buffer_dim, MPI_LOGICAL,          0, comm, ierr)
      call MPI_Bcast(rbuff,           buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)
      call MPI_Bcast(cbuff, cbuff_len*buffer_dim, MPI_CHARACTER,        0, comm, ierr)

      if (slave) then

         nsub                = ibuff(1)
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

      gpot(:,:,:) = 0.0

      if (.not.user_grav) then
         grav_pot_3d => default_grav_pot_3d
#ifdef VERBOSE
         if (master) call warn("[gravity:init_grav] user_grav is set to false. Using default grav_pot_3d.")
#endif /* VERBOSE */
      endif

   end subroutine init_grav

   subroutine source_terms_grav

#ifdef SELF_GRAV
      use arrays,            only: u, sgp, sgpm
      use fluidindex,        only: iarr_all_sg
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

      sgpm = sgp

#ifdef POISSON_FFT
      call poisson_solve( sum(u(iarr_all_sg,:,:,:),1) )
#endif /* POISSON_FFT */
#ifdef MULTIGRID
      if (size(iarr_all_sg) == 1) then
         call multigrid_solve_grav(u(iarr_all_sg(1),:,:,:))
      else
         call multigrid_solve_grav( sum(u(iarr_all_sg,:,:,:),1) )
         !>
         !! \deprecated BEWARE Here a lot of heap space is required and some compilers may generate code that do segfaults for big enough domains.
         !! It is the weakest point of this type in Maclaurin test. Next one (in fluidboundaries.F90) is 8 times less sensitive.
         !<
      endif
#endif /* MULTIGRID */

      ! communicate boundary values for sgp(:, :, :) because multigrid solver gives at most 2 guardcells, while for hydro solver typically 4 is required.
      call all_sgp_boundaries
      if (frun) then
         sgpm = sgp
         frun = .false.
      endif
#endif /* SELF_GRAV */
      if (variable_gp) call grav_pot_3d
      call sum_potential

   end subroutine source_terms_grav

   subroutine sum_potential

      use arrays,   only: gpot, gp, hgpot
      use mpisetup, only: dt, dtm
#ifdef SELF_GRAV
      use arrays,   only: sgp, sgpm
#endif /* SELF_GRAV */

      implicit none
      real :: h

      if (dtm /= 0) then
         h = dt/dtm
      else
         h = 0.0
      endif

#ifdef SELF_GRAV
      gpot  = gp + (1.+h)    *sgp -     h*sgpm
      hgpot = gp + (1.+0.5*h)*sgp - 0.5*h*sgpm
#else /* !SELF_GRAV */
      !> \deprecated BEWARE: as long as grav_pot_3d is called only in init_piernik this assignment probably don't need to be repeated more than once
      gpot  = gp
      hgpot = gp
#endif /* !SELF_GRAV */

   end subroutine sum_potential

#ifdef SELF_GRAV

!> \warning An improper evaluation of guardcell potential may occur when the multigrid boundary conditions doesn't match /BOUNDARIES/ namelist (e.g. isolated on periodic domain).

   subroutine all_sgp_boundaries

      use arrays,        only: sgp
      use dataio_pub,    only: die
      use grid,          only: cg
      use mpi,           only: MPI_STATUS_SIZE, MPI_REQUEST_NULL
      use mpisetup,      only: comm3d, ierr, procxl, procxr, procyl, procyr, proczl, proczr, psize, &
           &                   bnd_xl, bnd_xr, bnd_yl, bnd_yr, bnd_zl, bnd_zr, has_dir
      use constants,     only: xdim, ydim, zdim

      implicit none

      integer, parameter                        :: nreq = 3 * 4
      integer, dimension(nreq)                  :: req3d
      integer, dimension(MPI_STATUS_SIZE, nreq) :: status3d
      integer                                   :: i

      req3d(:) = MPI_REQUEST_NULL

      if (has_dir(xdim)) then

         select case (bnd_xl)
            case ('per')
               do i = 1, ceiling(cg%nb/real(cg%nxb)) ! Repeating is important for domains that are narrower than their guardcells (e.g. cg%nxb = 2)
                  sgp(1:cg%nb, :, :) = sgp(cg%ieb:cg%ie, :, :)
               enddo
            case ('mpi')
               if (psize(xdim) > 1) then
                  call MPI_Isend(sgp(1, 1, 1), 1, cg%ARR_YZ_LEFT_DOM,  procxl, 12, comm3d, req3d(1), ierr)
                  call MPI_Irecv(sgp(1, 1, 1), 1, cg%ARR_YZ_LEFT_BND,  procxl, 22, comm3d, req3d(2), ierr)
               else
                  call die("[gravity:all_grav_boundaries] bnd_xl == 'mpi' && psize(xdim) <= 1")
               endif
            case ('she') !> \todo move appropriate code from poissonsolver::poisson_solve or do nothing. Or die until someone really needs SHEAR.
                call die("[gravity:all_grav_boundaries] bnd_xl == 'she' not implemented")
            case default ! Set gradient == 0 on the boundaries
               do i = 1, cg%nb
                  sgp(i, :, :) = sgp(cg%is, :, :)
               enddo
         end select

         select case (bnd_xr)
            case ('per')
               do i = 1, ceiling(cg%nb/real(cg%nxb))
                  sgp(cg%ie+1:cg%nx, :, :) = sgp(cg%is:cg%isb, :, :)
               enddo
            case ('mpi')
               if (psize(xdim) > 1) then
                  call MPI_Isend(sgp(1, 1, 1), 1, cg%ARR_YZ_RIGHT_DOM, procxr, 22, comm3d, req3d(3), ierr)
                  call MPI_Irecv(sgp(1, 1, 1), 1, cg%ARR_YZ_RIGHT_BND, procxr, 12, comm3d, req3d(4), ierr)
               else
                  call die("[gravity:all_grav_boundaries] bnd_xr == 'mpi' && psize(xdim) <= 1")
               endif
            case ('she')
                call die("[gravity:all_grav_boundaries] bnd_xl == 'she' not implemented")
            case default
               do i = 1, cg%nb
                  sgp(cg%ie+i, :, :) = sgp(cg%ie, :, :)
               enddo
         end select

      endif

      if (has_dir(ydim)) then

         select case (bnd_yl)
            case ('per')
               do i = 1, ceiling(cg%nb/real(cg%nyb))
                  sgp(:, 1:cg%nb, :) = sgp(:, cg%jeb:cg%je, :)
               enddo
            case ('mpi')
               if (psize(ydim) > 1) then
                  call MPI_Isend(sgp(1, 1, 1), 1, cg%ARR_XZ_LEFT_DOM,  procyl, 32, comm3d, req3d(5), ierr)
                  call MPI_Irecv(sgp(1, 1, 1), 1, cg%ARR_XZ_LEFT_BND,  procyl, 42, comm3d, req3d(6), ierr)
               else
                  call die("[gravity:all_grav_boundaries] bnd_yl == 'mpi' && psize(ydim) <= 1")
               endif
            case default
               do i = 1, cg%nb
                  sgp(:, i, :) = sgp(:, cg%js, :)
               enddo
         end select

         select case (bnd_yr)
            case ('per')
               do i = 1, ceiling(cg%nb/real(cg%nyb))
                  sgp(:, cg%je+1:cg%ny, :) = sgp(:, cg%js:cg%jsb, :)
               enddo
            case ('mpi')
               if (psize(ydim) > 1) then
                  call MPI_Isend(sgp(1, 1, 1), 1, cg%ARR_XZ_RIGHT_DOM, procyr, 42, comm3d, req3d(7), ierr)
                  call MPI_Irecv(sgp(1, 1, 1), 1, cg%ARR_XZ_RIGHT_BND, procyr, 32, comm3d, req3d(8), ierr)
               else
                  call die("[gravity:all_grav_boundaries] bnd_yr == 'mpi' && psize(ydim) <= 1")
               endif
            case default
               do i = 1, cg%nb
                  sgp(:, cg%je+i, :) = sgp(:, cg%je, :)
               enddo
         end select

      endif

      if (has_dir(zdim)) then

         select case (bnd_zl)
            case ('per')
               do i = 1, ceiling(cg%nb/real(cg%nzb))
                  sgp(:, :, 1:cg%nb) = sgp(:, :, cg%keb:cg%ke)
               enddo
            case ('mpi')
               if (psize(zdim) > 1) then
                  call MPI_Isend(sgp(1, 1, 1), 1, cg%ARR_XY_LEFT_DOM,  proczl, 52, comm3d, req3d(9), ierr)
                  call MPI_Irecv(sgp(1, 1, 1), 1, cg%ARR_XY_LEFT_BND,  proczl, 62, comm3d, req3d(10), ierr)
               else
                  call die("[gravity:all_grav_boundaries] bnd_zl == 'mpi' && psize(zdim) <= 1")
               endif
            case default
               do i = 1, cg%nb
                  sgp(:, :, i) = sgp(:, :, cg%ks)
               enddo
         end select

         select case (bnd_zr)
            case ('per')
               do i = 1, ceiling(cg%nb/real(cg%nzb))
                  sgp(:, :, cg%ke+1:cg%nz) = sgp(:, :, cg%ks:cg%ksb)
               enddo
            case ('mpi')
               if (psize(zdim) > 1) then
                  call MPI_Isend(sgp(1, 1, 1), 1, cg%ARR_XY_RIGHT_DOM, proczr, 62, comm3d, req3d(11), ierr)
                  call MPI_Irecv(sgp(1, 1, 1), 1, cg%ARR_XY_RIGHT_BND, proczr, 52, comm3d, req3d(12), ierr)
               else
                  call die("[gravity:all_grav_boundaries] bnd_zr == 'mpi' && psize(zdim) <= 1")
               endif
            case default
               do i = 1, cg%nb
                  sgp(:, :, cg%ke+i) = sgp(:, :, cg%ke)
               enddo
         end select

      endif

      call MPI_Waitall(nreq, req3d(:), status3d(:,:), ierr)

   end subroutine all_sgp_boundaries

#endif /* SELF_GRAV */

   subroutine grav_null(gp, ax, flatten)

      use types, only: axes

      implicit none

      real, dimension(:,:,:), intent(inout) :: gp
      type(axes), intent(in)                :: ax
      logical, intent(in), optional         :: flatten

      gp = 0.0

      if (.false. .and. flatten) gp(1, 1, 1) = ax%x(1) ! suppress compiler warnings on unused arguments

   end subroutine grav_null

   subroutine grav_uniform(gp, ax, flatten)

      use types, only: axes

      implicit none

      real, dimension(:,:,:), intent(inout) :: gp
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

      use types, only: axes

      implicit none

      real, dimension(:,:,:), intent(inout) :: gp
      type(axes), intent(in)                :: ax
      logical, intent(in), optional         :: flatten
      integer :: i, j, k

      do i = 1, ubound(gp,1)
         do j = 1, ubound(gp,2)
            do k = 1, ubound(gp,3)
               gp(i,j,k) = -0.5*(g_dir(1)*ax%x(i)**2 + g_dir(2)*ax%y(j)**2 + g_dir(3)*ax%z(k)**2)
            enddo
         enddo
      enddo

      if (.false. .and. flatten) i=0 ! suppress compiler warnings on unused arguments

   end subroutine grav_linear

   subroutine grav_ptmass_pure(gp, ax, flatten)

      use units,  only: newtong
      use types,      only: axes

      implicit none

      real, dimension(:,:,:), intent(inout) :: gp
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

   subroutine grav_ptmass_softened(gp, ax, flatten)

      use units,  only: newtong
      use fluidindex, only: flind
      use mpisetup,   only: smalld
      use types,      only: axes

      implicit none

      real, dimension(:,:,:), intent(inout) :: gp
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

      use units,    only: newtong
      use types,        only: axes

      implicit none

      real, dimension(:,:,:), intent(inout) :: gp
      type(axes), intent(in)                :: ax
      logical, intent(in), optional         :: flatten

      integer :: j, k
      real    :: r_smooth2
      real    :: GM1, GM2, z2, yz2

      r_smooth2 = r_smooth**2
      GM1 =  newtong * ptmass
      GM2 =  newtong * ptmass2

      do k = 1, ubound(gp,1)
         z2 = ax%z(k)**2
         do j = 1, ubound(gp,2)
            yz2 = ax%y(j)**2 + z2
            gp(:,j,k) =  - GM1 / sqrt((ax%x(:) - ptm_x)**2  + yz2 + r_smooth2) &
                 &       - GM2 / sqrt((ax%x(:) - ptm2_x)**2 + yz2 + r_smooth2) &
                 &       - 0.5 * Omega**2 * ((ax%x(:) - cmass_x)**2 + yz2)
         enddo
      enddo

      if (.false. .and. flatten) j=0 ! suppress compiler warnings on unused arguments

   end subroutine grav_roche

   subroutine grav_ptmass_stiff(gp, ax, flatten)

      use units,  only: newtong
      use types,      only: axes

      implicit none

      real, dimension(:,:,:), intent(inout) :: gp
      type(axes), intent(in)                :: ax
      logical, intent(in), optional         :: flatten

      integer :: i, j, k
      real    :: r_smooth2, r2, gmr, gm, z2, yz2
      ! promote stiff-body rotation inside smoothing length, don't affect the global potential outside

      r_smooth2 = r_smooth**2
      gm =  - newtong * ptmass
      gmr = 0.5 * gm / r_smooth

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

      use arrays,       only: gp
      use dataio_pub,   only: die, warn
      use grid,         only: cg
      use mpisetup,     only: master, geometry_type
      use types,        only: axes
      use constants,    only: GEO_XYZ

      implicit none
      type(axes) :: ax

      if (.not.allocated(ax%x)) allocate(ax%x(size(cg%x)))
      if (.not.allocated(ax%y)) allocate(ax%y(size(cg%y)))
      if (.not.allocated(ax%z)) allocate(ax%z(size(cg%z)))
      ax%x = cg%x
      ax%y = cg%y
      ax%z = cg%z

      gp_status = ''

      if (geometry_type /= GEO_XYZ) then
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
            call grav_null(gp,ax)                    ; grav_type => grav_null
         case ("linear", "grav_lin", "GRAV_LINEAR")
            call grav_linear(gp,ax)                  ; grav_type => grav_linear
         case ("uniform", "grav_unif", "GRAV_UNIFORM")
            call grav_uniform(gp,ax)                 ; grav_type => grav_uniform
         case ("softened ptmass", "ptmass_soft", "GRAV_PTMASS")
            call grav_ptmass_softened(gp,ax,.false.) ; grav_type => grav_ptmass_softened
         case ("stiff ptmass", "ptmass_stiff", "GRAV_PTMASSSTIFF")
            call grav_ptmass_stiff(gp,ax)            ; grav_type => grav_ptmass_stiff
         case ("ptmass", "ptmass_pure", "GRAV_PTMASSPURE")
            call grav_ptmass_pure(gp,ax,.false.)     ; grav_type => grav_ptmass_pure
         case ("flat softened ptmass", "flat_ptmass_soft", "GRAV_PTFLAT")
            call grav_ptmass_softened(gp,ax,.true.)  ; grav_type => grav_ptmass_softened
         case ("flat ptmass", "flat_ptmass")
            call grav_ptmass_pure(gp,ax,.true.)      ; grav_type => grav_ptmass_pure
         case ("roche", "grav_roche", "GRAV_ROCHE")
#ifndef CORIOLIS
            call die("[gravity:default_grav_pot_3d] define CORIOLIS in piernik.def for Roche potential")
#endif /* !CORIOLIS */
            call grav_roche(gp,ax)                   ; grav_type => grav_roche
         case ("user", "grav_user", "GRAV_USER")
            call die("[gravity:default_grav_pot_3d] user 'grav_pot_3d' should be defined in initprob!")
         case default
            gp_status = 'undefined'
      end select

!-----------------------

      if (gp_status .eq. 'undefined') then
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

   end subroutine default_grav_pot_3d

!>
!! \brief Routine that compute values of gravitational acceleration using gravitational potential array gp
!<
   subroutine grav_pot2accel(sweep, i1,i2, n, grav,istep)

      use arrays,    only: gpot, hgpot
      use grid,      only: cg
      use constants, only: xdim, ydim, zdim

      implicit none

      integer, intent(in)            :: sweep      !< string of characters that points out the current sweep direction
      integer, intent(in)            :: i1         !< number of column in the first direction after one pointed out by sweep
      integer, intent(in)            :: i2         !< number of column in the second direction after one pointed out by sweep
      integer, intent(in)            :: n          !< number of elements of returned array grav
      real, dimension(n),intent(out) :: grav       !< 1D array of gravitational acceleration values computed for positions from %xsw and returned by the routine
      integer, intent(in)            :: istep      !< istep=1 for halfstep, istep=2 for fullstep
!> \todo offer high order gradient as an option in parameter file
!      real, parameter :: onetw = 1./12.

! Gravitational acceleration is computed on right cell boundaries

!      if (istep==1) then
!         select case (sweep)
!            case (xdim)
!               grav(3:n-2) = onetw*(hgpot(5:n,i1,i2) - 8.*hgpot(4:n-1,i1,i2) + 8.*hgpot(2:n-3,i1,i2) - hgpot(1:n-4,i1,i2) )/dl(xdim)
!            case (ydim)
!               grav(3:n-2) = onetw*(hgpot(i2,5:n,i1) - 8.*hgpot(i2,4:n-1,i1) + 8.*hgpot(i2,2:n-3,i1) - hgpot(i2,1:n-4,i1) )/dl(xdim)
!            case (zdim)
!               grav(3:n-2) = onetw*(hgpot(i1,i2,5:n) - 8.*hgpot(i1,i2,4:n-1) + 8.*hgpot(i1,i2,2:n-3) - hgpot(i1,i2,1:n-4) )/dl(xdim)
!         end select
!      else
!         select case (sweep)
!            case (xdim)
!               grav(3:n-2) = onetw*(gpot(5:n,i1,i2) - 8.*gpot(4:n-1,i1,i2) + 8.*gpot(2:n-3,i1,i2) - gpot(1:n-4,i1,i2) )/dl(xdim)
!            case (ydim)
!               grav(3:n-2) = onetw*(gpot(i2,5:n,i1) - 8.*gpot(i2,4:n-1,i1) + 8.*gpot(i2,2:n-3,i1) - gpot(i2,1:n-4,i1) )/dl(xdim)
!            case (zdim)
!               grav(3:n-2) = onetw*(gpot(i1,i2,5:n) - 8.*gpot(i1,i2,4:n-1) + 8.*gpot(i1,i2,2:n-3) - gpot(i1,i2,1:n-4) )/dl(xdim)
!         end select
!      endif
!      grav(2) = grav(3); grav(n-1) = grav(n-2)
!      grav(1) = grav(2); grav(n) = grav(n-1)
      if (istep==1) then
         select case (sweep)
            case (xdim)
               grav(2:n-1) = 0.5*(hgpot(1:n-2,i1,i2) - hgpot(3:n,i1,i2))/cg%dl(xdim)
            case (ydim)
               grav(2:n-1) = 0.5*(hgpot(i2,1:n-2,i1) - hgpot(i2,3:n,i1))/cg%dl(ydim)
            case (zdim)
               grav(2:n-1) = 0.5*(hgpot(i1,i2,1:n-2) - hgpot(i1,i2,3:n))/cg%dl(zdim)
         end select

      else
         select case (sweep)
            case (xdim)
               grav(2:n-1) = 0.5*(gpot(1:n-2,i1,i2) - gpot(3:n,i1,i2))/cg%dl(xdim)
            case (ydim)
               grav(2:n-1) = 0.5*(gpot(i2,1:n-2,i1) - gpot(i2,3:n,i1))/cg%dl(ydim)
            case (zdim)
               grav(2:n-1) = 0.5*(gpot(i1,i2,1:n-2) - gpot(i1,i2,3:n))/cg%dl(zdim)
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

      use arrays,    only: gp
      use grid,      only: cg
      use mpi,       only: MPI_DOUBLE_PRECISION
      use mpisetup,  only: psize, pcoords, master, nproc, comm, comm3d, err, ierr, mpifind
      use constants, only: xdim, ydim, zdim, ndims

      implicit none

      integer                             :: i, j, k, ip, pgpmax
      real, allocatable, dimension(:,:,:) :: gpwork
      real                                :: gravrx(cg%nx), gravry(cg%ny), gravrz(cg%nz)
      real                                :: gp_max
      integer, dimension(3)               :: loc_gp_max
      integer                             :: proc_gp_max
      integer                             :: px, py, pz, pc(3)
      real                                :: dgpx_proc, dgpx_all(0:nproc-1), &
                                             dgpy_proc, dgpy_all(0:nproc-1), &
                                             dgpz_proc, dgpz_all(0:nproc-1), &
                                             dgpx(0:psize(xdim)-1,0:psize(ydim)-1,0:psize(zdim)-1), &
                                             dgpy(0:psize(xdim)-1,0:psize(ydim)-1,0:psize(zdim)-1), &
                                             dgpz(0:psize(xdim)-1,0:psize(ydim)-1,0:psize(zdim)-1), &
                                             ddgp(0:psize(xdim)-1,0:psize(ydim)-1,0:psize(zdim)-1)

      allocate(gpwork(cg%nx, cg%ny, cg%nz))
      gpwork(1,1,1) = 0.0

      call grav_accel(xdim, 1, 1, cg%xr(:), cg%nx, gravrx)
      do i = 1, cg%nx-1
         gpwork(i+1,1,1) = gpwork(i,1,1) - gravrx(i)*cg%dl(xdim)
      enddo

      do i=1, cg%nx
         call grav_accel(ydim, 1, i, cg%yr(:), cg%ny, gravry)
         do j = 1, cg%ny-1
            gpwork(i,j+1,1) = gpwork(i,j,1) - gravry(j)*cg%dl(ydim)
         enddo
      enddo

      do i=1, cg%nx
         do j=1, cg%ny
            call grav_accel(zdim, i, j, cg%zr(:), cg%nz, gravrz)
            do k = 1, cg%nz-1
               gpwork(i,j,k+1) = gpwork(i,j,k) - gravrz(k)*cg%dl(zdim)
            enddo
         enddo
      enddo

      dgpx_proc = gpwork(cg%is, 1, 1)-gpwork(1,1,1)
      dgpy_proc = gpwork(1, cg%js, 1)-gpwork(1,1,1)
      dgpz_proc = gpwork(1, 1, cg%ks)-gpwork(1,1,1)

      call MPI_Gather ( dgpx_proc, 1, MPI_DOUBLE_PRECISION, dgpx_all, 1, MPI_DOUBLE_PRECISION, 0, comm3d, ierr )
      call MPI_Gather ( dgpy_proc, 1, MPI_DOUBLE_PRECISION, dgpy_all, 1, MPI_DOUBLE_PRECISION, 0, comm3d, ierr )
      call MPI_Gather ( dgpz_proc, 1, MPI_DOUBLE_PRECISION, dgpz_all, 1, MPI_DOUBLE_PRECISION, 0, comm3d, ierr )

      if (master) then

         do ip = 0, nproc-1
            call MPI_Cart_coords(comm3d, ip, ndims, pc, ierr)

            px = pc(1)
            py = pc(2)
            pz = pc(3)

            dgpx(px,py,pz) = dgpx_all(ip)
            dgpy(px,py,pz) = dgpy_all(ip)
            dgpz(px,py,pz) = dgpz_all(ip)

         enddo

         ddgp(0,0,0) = 0.0
         do i = 1, psize(xdim)-1
            ddgp(i,0,0) = ddgp(i-1,0,0) + dgpx(i-1,0,0)
         enddo

         do i=0, psize(xdim)-1
            do j = 1, psize(ydim)-1
               ddgp(i,j,0) = ddgp(i,j-1,0) + dgpy(i,j-1,0)
            enddo
         enddo

         do i=0, psize(xdim)-1
            do j=0, psize(ydim)-1
               do k = 1, psize(zdim)-1
                  ddgp(i,j,k)= ddgp(i,j,k-1) + dgpz(i,j,k-1)
               enddo
            enddo
         enddo

      endif

      call MPI_Bcast(ddgp, nproc, MPI_DOUBLE_PRECISION, 0, comm, ierr)

      px = pcoords(1)
      py = pcoords(2)
      pz = pcoords(3)

      gpwork = gpwork + ddgp(px,py,pz)

      gp_max      = maxval(gpwork(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke))
      loc_gp_max  = maxloc(gpwork(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke)) &
                  + [ cg%nb, cg%nb, cg%nb ]

      call mpifind(gp_max, 'max', loc_gp_max, proc_gp_max)
      pgpmax = proc_gp_max

      call MPI_Bcast(gp_max, 1, MPI_DOUBLE_PRECISION, pgpmax, comm, ierr)
      gpwork = gpwork - gp_max

      gp=gpwork
      if (allocated(gpwork)) deallocate(gpwork)

   end subroutine grav_accel2pot

end module gravity
