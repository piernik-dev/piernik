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
!! \brief [DW] Module containing all main subroutines and functions that govern %gravity force in the code
!! \todo to check importance and usefulness of such parameters as g_y and n_gravr2 and tune_zeq_bnd
!!
!!
!! In this module a namelist of parameters is specified:
!! \copydetails gravity::init_grav
!<
module gravity
! pulled by GRAV
   implicit none

   private
   public :: init_grav, grav_accel, source_terms_grav, grav_pot2accel, grav_pot_3d, get_gprofs, grav_accel2pot
   public :: g_dir, r_gc, ptmass, ptm_x, ptm_y, ptm_z, r_smooth, nsub, tune_zeq, tune_zeq_bnd, h_grav, r_grav, n_gravr, n_gravr2, n_gravh, user_grav, gp_status, gprofs_target

   integer, parameter         :: gp_stat_len = 9
   integer, parameter         :: gproft_len  = 5
   character(LEN=gp_stat_len) :: gp_status       !< variable set as 'undefined' in grav_pot_3d when grav_accel is supposed to use
   character(LEN=gproft_len)  :: gprofs_target   !< variable set pointing gravity routine in hydrostatic_zeq ('accel' or ready gp array 'gparr')
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
   integer :: n_gravr2              !< similar to n_gravr <b>(currently not used)</b>
   real    :: tune_zeq              !< z-component of %gravity tuning factor used by hydrostatic_zeq
   real    :: tune_zeq_bnd          !< z-component of %gravity tuning factor supposed to be used in boundaries <b>(currently not used)</b>

   logical :: user_grav             !< use user defined grav_pot_3d

   interface
      subroutine user_grav_pot_3d
         implicit none
      end subroutine user_grav_pot_3d
      subroutine gprofs_default(iia,jja)
         implicit none
         integer, intent(in) :: iia, jja
      end subroutine gprofs_default
      subroutine user_grav_accel(sweep, i1,i2, xsw, n, grav)
         implicit none
         character(len=*), intent(in)   :: sweep
         integer, intent(in)            :: i1, i2
         integer, intent(in)            :: n
         real, dimension(n),intent(in)  :: xsw
         real, dimension(n),intent(out) :: grav
      end subroutine user_grav_accel
   end interface

   procedure(user_grav_pot_3d), pointer :: grav_pot_3d => NULL()
   procedure(user_grav_accel),  pointer :: grav_accel  => NULL()
   procedure(gprofs_default),   pointer :: get_gprofs  => NULL()

   contains

!>
!! \brief Routine that sets the initial values of %gravity parameters from namelist @c GRAVITY
!!
!! \n \n
!! @b GRAVITY
!! \n \n
!! <table border="+1">
!! <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!! <tr><td>g_dir       </td><td>0.0 </td><td>real             </td><td>\copydoc gravity::g_dir        </td></tr>
!! <tr><td>r_gc        </td><td>8500</td><td>real             </td><td>\copydoc gravity::r_gc         </td></tr>
!! <tr><td>ptmass      </td><td>0.0 </td><td>non-negative real</td><td>\copydoc gravity::ptmass       </td></tr>
!! <tr><td>ptm_x       </td><td>0.0 </td><td>real             </td><td>\copydoc gravity::ptm_x        </td></tr>
!! <tr><td>ptm_y       </td><td>0.0 </td><td>real             </td><td>\copydoc gravity::ptm_y        </td></tr>
!! <tr><td>ptm_z       </td><td>0.0 </td><td>real             </td><td>\copydoc gravity::ptm_z        </td></tr>
!! <tr><td>r_smooth    </td><td>0.0 </td><td>real             </td><td>\copydoc gravity::r_smooth     </td></tr>
!! <tr><td>nsub        </td><td>10  </td><td>integer > 0      </td><td>\copydoc gravity::nsub         </td></tr>
!! <tr><td>tune_zeq    </td><td>1.0 </td><td>real             </td><td>\copydoc gravity::tune_zeq     </td></tr>
!! <tr><td>tune_zeq_bnd</td><td>1.0 </td><td>real             </td><td>\copydoc gravity::tune_zeq_bnd </td></tr>
!! <tr><td>h_grav      </td><td>1.e6</td><td>real             </td><td>\copydoc gravity::h_grav       </td></tr>
!! <tr><td>r_grav      </td><td>1.e6</td><td>real             </td><td>\copydoc gravity::r_grav       </td></tr>
!! <tr><td>n_gravr     </td><td>0   </td><td>real             </td><td>\copydoc gravity::n_gravr      </td></tr>
!! <tr><td>n_gravr2    </td><td>0   </td><td>real             </td><td>\copydoc gravity::n_gravr2     </td></tr>
!! <tr><td>n_gravh     </td><td>0   </td><td>real             </td><td>\copydoc gravity::n_gravh      </td></tr>
!! </table>
!! \n \n
!<
   subroutine init_grav

      use arrays,        only: gpot
      use dataio_pub,    only: ierrh, par_file, namelist_errh, compare_namelist    ! QA_WARN required for diff_nml
      use dataio_pub,    only: warn
      use mpisetup,      only: ibuff, rbuff, cbuff, cbuff_len, buffer_dim, comm, ierr, proc, &
           &                   MPI_DOUBLE_PRECISION, MPI_INTEGER, MPI_LOGICAL, lbuff, MPI_CHARACTER

      implicit none


      namelist /GRAVITY/ g_dir, r_gc, ptmass, ptm_x, ptm_y, ptm_z, r_smooth, &
                nsub, tune_zeq, tune_zeq_bnd, h_grav, r_grav, n_gravr, n_gravr2, n_gravh, user_grav, gprofs_target

#ifdef VERBOSE
      call warn("[gravity:init_grav] Commencing gravity module initialization")
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
      n_gravr = 0
      n_gravr2= 0
      n_gravh = 0
      gprofs_target = 'gparr'

      user_grav = .false.

      if (proc == 0) then

         diff_nml(GRAVITY)

         ibuff(1)   = nsub
         ibuff(2)   = n_gravr
         ibuff(3)   = n_gravr2
         ibuff(4)   = n_gravh

         rbuff(1:3) = g_dir
         rbuff(5)   = r_gc
         rbuff(6)   = ptmass
         rbuff(7)   = ptm_x
         rbuff(8)   = ptm_y
         rbuff(9)   = ptm_z
         rbuff(10)  = tune_zeq
         rbuff(11)  = tune_zeq_bnd
         rbuff(12)  = r_smooth
         rbuff(13)  = h_grav
         rbuff(14)  = r_grav

         lbuff(1)   = user_grav
         cbuff(1)   = gprofs_target

      endif

      call MPI_Bcast(ibuff, buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_Bcast(lbuff, buffer_dim, MPI_LOGICAL,          0, comm, ierr)
      call MPI_Bcast(rbuff, buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)
      call MPI_Bcast(cbuff, cbuff_len*buffer_dim, MPI_CHARACTER, 0, comm, ierr)

      if (proc /= 0) then

         nsub                = ibuff(1)
         n_gravr             = ibuff(2)
         n_gravr2            = ibuff(3)
         n_gravh             = ibuff(4)

         g_dir               = rbuff(1:3)
         r_gc                = rbuff(5)
         ptmass              = rbuff(6)
         ptm_x               = rbuff(7)
         ptm_y               = rbuff(8)
         ptm_z               = rbuff(9)
         tune_zeq            = rbuff(10)
         tune_zeq_bnd        = rbuff(11)
         r_smooth            = rbuff(12)
         h_grav              = rbuff(13)
         r_grav              = rbuff(14)

         user_grav           = lbuff(1)
         gprofs_target       = cbuff(1)(1:gproft_len)

      endif

      gpot(:,:,:) = 0.0

      if (.not.user_grav) then
         grav_pot_3d => default_grav_pot_3d
#ifdef VERBOSE
         call warn("[gravity:init_grav] user_grav is set to false. Using default grav_pot_3d.")
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
         ! BEWARE Here a lot of heap space is required and some compilers may generate code that do segfaults for big enough domains.
         ! It is the weakest point of this type in Maclaurin test. Next one (in fluidboundaries.F90) is 8 times less sensitive.
      endif
#endif /* MULTIGRID */

      ! communicate boundary values for sgp(:, :, :) because multigrid solver gives at most 2 guardcells, while for hydro solver typically 4 is required.
      call all_sgp_boundaries
      if (frun) then
         sgpm = sgp
         frun = .false.
      endif
#endif /* SELF_GRAV */

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
      !// BEWARE: as long as grav_pot_3d is called only in init_piernik this assignment probably don't need to be repeated more than once
      gpot  = gp
      hgpot = gp
#endif /* !SELF_GRAV */

   end subroutine sum_potential

#ifdef SELF_GRAV

!// An improper evaluation of guardcell potential may occur when the multigrid boundary conditions doesn't match /BOUNDARIES/ namelist (e.g. isolated on periodic domain).

   subroutine all_sgp_boundaries

      use arrays,        only: sgp
      use dataio_pub,    only: die
      use grid,          only: nb, nxd, nyd, nzd
      use mpisetup,      only: comm3d, ierr, MPI_STATUS_SIZE, MPI_REQUEST_NULL, &
           &                   procxl, procxr, procyl, procyr, proczl, proczr, proc, pxsize, pysize, pzsize, &
           &                   bnd_xl, bnd_xr, bnd_yl, bnd_yr, bnd_zl, bnd_zr, &
           &                   ARR_YZ_LEFT_BND, ARR_YZ_RIGHT_BND, ARR_YZ_LEFT_DOM, ARR_YZ_RIGHT_DOM, &
           &                   ARR_XZ_LEFT_BND, ARR_XZ_RIGHT_BND, ARR_XZ_LEFT_DOM, ARR_XZ_RIGHT_DOM, &
           &                   ARR_XY_LEFT_BND, ARR_XY_RIGHT_BND, ARR_XY_LEFT_DOM, ARR_XY_RIGHT_DOM

      implicit none

      integer, parameter                        :: nreq = 3 * 4
      integer, dimension(nreq)                  :: req3d
      integer, dimension(MPI_STATUS_SIZE, nreq) :: status3d
      integer                                   :: i

      req3d(:) = MPI_REQUEST_NULL

      if (nxd /= 1) then

         select case (bnd_xl)
            case ('per')
               do i = 1, ceiling(nb/real(nxd)) ! Repeating is important for domains that are narrower than their guardcells (e.g. nxd = 2)
                  sgp(1:nb, :, :) = sgp(nxd+1:nxd+nb, :, :)
               enddo
            case ('mpi')
               if (pxsize > 1) then
                  call MPI_Isend(sgp(1, 1, 1), 1, ARR_YZ_LEFT_DOM,  procxl, 12, comm3d, req3d(1), ierr)
                  call MPI_Irecv(sgp(1, 1, 1), 1, ARR_YZ_LEFT_BND,  procxl, 22, comm3d, req3d(2), ierr)
               else
                  call die("[gravity:all_grav_boundaries] bnd_xl == 'mpi' && pxsize <= 1")
               endif
            case ('she') ! move appropriate code from poissonsolver::poisson_solve or do nothing. Or die until someone really needs SHEAR.
                call die("[gravity:all_grav_boundaries] bnd_xl == 'she' not implemented")
            case default ! Set gradient == 0 on the boundaries
               do i = 1, nb
                  sgp(i, :, :) = sgp(nb+1, :, :)
               enddo
         end select

         select case (bnd_xr)
            case ('per')
               do i = 1, ceiling(nb/real(nxd))
                  sgp(nxd+nb+1:nxd+2*nb, :, :) = sgp(nb+1:2*nb, :, :)
               enddo
            case ('mpi')
               if (pxsize > 1) then
                  call MPI_Isend(sgp(1, 1, 1), 1, ARR_YZ_RIGHT_DOM, procxr, 22, comm3d, req3d(3), ierr)
                  call MPI_Irecv(sgp(1, 1, 1), 1, ARR_YZ_RIGHT_BND, procxr, 12, comm3d, req3d(4), ierr)
               else
                  call die("[gravity:all_grav_boundaries] bnd_xr == 'mpi' && pxsize <= 1")
               endif
            case ('she')
                call die("[gravity:all_grav_boundaries] bnd_xl == 'she' not implemented")
            case default
               do i = 1, nb
                  sgp(nb+nxd+i, :, :) = sgp(nb+nxd, :, :)
               enddo
         end select

      endif

      if (nyd /= 1) then

         select case (bnd_yl)
            case ('per')
               do i = 1, ceiling(nb/real(nyd))
                  sgp(:, 1:nb, :) = sgp(:, nyd+1:nyd+nb, :)
               enddo
            case ('mpi')
               if (pysize > 1) then
                  call MPI_Isend(sgp(1, 1, 1), 1, ARR_XZ_LEFT_DOM,  procyl, 32, comm3d, req3d(5), ierr)
                  call MPI_Irecv(sgp(1, 1, 1), 1, ARR_XZ_LEFT_BND,  procyl, 42, comm3d, req3d(6), ierr)
               else
                  call die("[gravity:all_grav_boundaries] bnd_yl == 'mpi' && pysize <= 1")
               endif
            case default
               do i = 1, nb
                  sgp(:, i, :) = sgp(:, nb+1, :)
               enddo
         end select

         select case (bnd_yr)
            case ('per')
               do i = 1, ceiling(nb/real(nyd))
                  sgp(:, nyd+nb+1:nyd+2*nb, :) = sgp(:, nb+1:2*nb, :)
               enddo
            case ('mpi')
               if (pysize > 1) then
                  call MPI_Isend(sgp(1, 1, 1), 1, ARR_XZ_RIGHT_DOM, procyr, 42, comm3d, req3d(7), ierr)
                  call MPI_Irecv(sgp(1, 1, 1), 1, ARR_XZ_RIGHT_BND, procyr, 32, comm3d, req3d(8), ierr)
               else
                  call die("[gravity:all_grav_boundaries] bnd_yr == 'mpi' && pysize <= 1")
               endif
            case default
               do i = 1, nb
                  sgp(:, nb+nyd+i, :) = sgp(:, nb+nyd, :)
               enddo
         end select

      endif

      if (nzd /= 1) then

         select case (bnd_zl)
            case ('per')
               do i = 1, ceiling(nb/real(nzd))
                  sgp(:, :, 1:nb) = sgp(:, :, nzd+1:nzd+nb)
               enddo
            case ('mpi')
               if (pzsize > 1) then
                  call MPI_Isend(sgp(1, 1, 1), 1, ARR_XY_LEFT_DOM,  proczl, 52, comm3d, req3d(9), ierr)
                  call MPI_Irecv(sgp(1, 1, 1), 1, ARR_XY_LEFT_BND,  proczl, 62, comm3d, req3d(10), ierr)
               else
                  call die("[gravity:all_grav_boundaries] bnd_zl == 'mpi' && pzsize <= 1")
               endif
            case default
               do i = 1, nb
                  sgp(:, :, i) = sgp(:, :, nb+1)
               enddo
         end select

         select case (bnd_zr)
            case ('per')
               do i = 1, ceiling(nb/real(nzd))
                  sgp(:, :, nzd+nb+1:nzd+2*nb) = sgp(:, :, nb+1:2*nb)
               enddo
            case ('mpi')
               if (pzsize > 1) then
                  call MPI_Isend(sgp(1, 1, 1), 1, ARR_XY_RIGHT_DOM, proczr, 62, comm3d, req3d(11), ierr)
                  call MPI_Irecv(sgp(1, 1, 1), 1, ARR_XY_RIGHT_BND, proczr, 52, comm3d, req3d(12), ierr)
               else
                  call die("[gravity:all_grav_boundaries] bnd_zr == 'mpi' && pzsize <= 1")
               endif
            case default
               do i = 1, nb
                  sgp(:, :, nb+nzd+i) = sgp(:, :, nb+nzd)
               enddo
         end select

      endif

      call MPI_Waitall(nreq, req3d(:), status3d(:,:), ierr)

   end subroutine all_sgp_boundaries

#endif /* SELF_GRAV */

   subroutine grav_null
      use arrays, only: gp
      implicit none

      gp = 0.0

   end subroutine grav_null

   subroutine grav_uniform
      use arrays, only: gp
      use grid,   only: nx, ny, nz, x, y, z
      implicit none
      integer :: i, j, k
      do i = 1, nx
         do j = 1, ny
            do k = 1, nz
               gp(:,:,i) = -(g_dir(1)*x(i) + g_dir(2)*y(j) + g_dir(3)*z(i))
            enddo
         enddo
      enddo
   end subroutine grav_uniform

   subroutine grav_linear
      use arrays, only: gp
      use grid,   only: nx, ny, nz, x, y, z
      implicit none
      integer :: i, j, k
      do i = 1, nx
         do j = 1, ny
            do k = 1, nz
               gp(:,:,i) = -0.5*(g_dir(1)*x(i)**2 + g_dir(2)*y(j)**2 + g_dir(3)*z(i)**2)
            enddo
         enddo
      enddo
   end subroutine grav_linear

   subroutine grav_ptmass_pure(flatten)
      use arrays,     only: gp
      use constants,  only: newtong
      use grid,       only: nx, ny, x, y, z
      implicit none
      logical, intent(in) :: flatten
      integer             :: i, j
      real                :: rc2, GM, x2

      r_smooth2 = r_smooth**2
      GM        = newtong*ptmass

       do i = 1, nx
          x2 = (x(i) - ptm_x)**2
          do j = 1, ny
             rc2 = x2 + (y(j) - ptm_y)**2

             if (flatten) then
                gp(i,j,:) = -GM / sqrt(rc2)
             else 
                gp(i,j,:) = -GM / sqrt( (z(:) - ptm_z)**2 + rc2 )
             endif

          enddo
       enddo
   end subroutine grav_ptmass_pure

   subroutine grav_ptmass_softened(flatten)
      use arrays,     only: gp
      use constants,  only: newtong
      use grid,       only: nx, ny, x, y, z
      use initfluids, only: cs_iso2
      use mpisetup,   only: smalld
      implicit none
      logical, intent(in) :: flatten
      integer             :: i, j
      real                :: rc2, r_smooth2, GM, fr, x2

      r_smooth2 = r_smooth**2
      GM        = newtong*ptmass

      do i = 1, nx
         x2 = (x(i) - ptm_x)**2
         do j = 1, ny
            rc2 = x2 + (y(j) - ptm_y)**2
            fr  = min( (sqrt(rc2)/r_grav)**n_gravr , 100.0)    ! BEWARE: hardcoded value
            fr  = max( 1./cosh(fr), smalld*1.e-2)              ! BEWARE: hadrcoded value
            fr  = -cs_iso2 * log(fr)

            if (flatten) then
               gp(i,j,:) = -GM / sqrt( rc2 + r_smooth2 ) + fr
            else
               gp(i,j,:) = -GM / sqrt( (z(:) - ptm_z)**2 + rc2 + r_smooth2 ) + fr
            endif

         enddo
      enddo
   end subroutine grav_ptmass_softened

   subroutine grav_ptmass_stiff
      use arrays,     only: gp
      use constants,  only: newtong
      use grid,       only: nx, ny, nz, x, y, z
      use initfluids, only: cs_iso2
      implicit none
      integer :: i, j, k
      real    :: r_smooth2, r2, gmr, gm, z2, yz2
      !// promote stiff-body rotation inside smoothing length, don't affect the global potential outside

      r_smooth2 = r_smooth**2 ! can be used also i other GRAV_PTMASS* clauses
      gm =  - newtong * ptmass
      gmr = 0.5 * gm / r_smooth

      do k = 1, nz
         z2 = (z(k) - ptm_z)**2
         do j = 1, ny
            yz2 = z2 + (y(j) - ptm_y)**2
            do i = 1, nx
               r2 = yz2 + (x(i) - ptm_x)**2
               if (r2 < r_smooth2) then
                  gp(i,j,k) = gmr * (3. - r2/r_smooth2)
               else
                  gp(i,j,k) = gm / dsqrt(r2)
               endif
            enddo
         enddo
      enddo
   end subroutine grav_ptmass_stiff
   
!--------------------------------------------------------------------------
!>
!! \brief Routine that compute values of gravitational potential filling in gp array and setting gp_status character string \n\n
!! The type of %gravity is governed by preprocessor: \n\n
!! \details
!! @b GRAV_NULL - gravitational potential array is set to zero \n\n
!! @b GRAV_UNIFORM - uniform type of %gravity in z-direction \n
!! \f$\Phi\left(z\right)= - const \cdot z \f$\n
!! where \f$ const \f$ is set by parameter @c g_z \n\n
!! @b GRAV_LINEAR - linear type of %gravity growing along z-direction \n
!! \f$\Phi\left(z\right) = -1/2 \cdot const \cdot z^2\f$ \n\n
!! @b GRAV_PTMASS - softened point mass type of %gravity \n
!! \f$\Phi\left(x,y,z\right)= - GM/\sqrt{x^2+y^2+z^2+r_{soft}^2}\f$ \n
!! where \f$r_{soft}\f$ is a radius of softening\n\n
!! @b GRAV_PTMASSPURE - unsoftened point mass type of %gravity \n
!! \f$\Phi\left(x,y,z\right)= - GM/\sqrt{x^2+y^2+z^2}\f$ \n\n
!! @b GRAV_PTMASSSTIFF - softened point mass type of %gravity with stiff-body rotation inside softening radius\n
!! \f$\Phi\left(x,y,z\right)= - GM/\sqrt{x^2+y^2+z^2}\f$ for \f$r > r_{soft}\f$ and \f$ GM/r_{soft} \left( - 3/2 + 1/2 {x^2+y^2+z^2}/r_{soft}^2 \right)\f$ inside \f$r_{soft}\f$ \n\n
!! @b GRAV_PTFLAT - planar, softened point mass type of %gravity \n
!! \f$\Phi\left(x,y,z\right)= - GM/\sqrt{x^2+y^2+r_{soft}^2}\f$ \n
!! where \f$r_{soft}\f$ is a radius of softening\n\n
!! @b GRAV_USER - not a standard type of %gravity, implemented by user in the routine grav_pot_user from gravity_user module.\n\n
!<

   subroutine default_grav_pot_3d

      use dataio_pub,   only: die, warn
      use arrays,       only: gp
      use grid,         only: nx, ny, nz, x, y, z
      use initfluids,   only: cs_iso2
      use mpisetup,     only: smalld
#ifdef GRAV_PTMTYPE
      use constants,    only: newtong
#endif /* GRAV_PTMTYPE */
#ifdef GRAV_USER
      use gravity_user, only: grav_pot_user
#endif /* GRAV_USER */

      implicit none

#ifdef GRAV_PTMTYPE
      integer :: i, j, k
      real    :: r_smooth2, gm, gmr, z2, yz2
#endif /* GRAV_PTMTYPE */
#if defined GRAV_UNIFORM || defined GRAV_LINEAR
      integer :: i
#endif /* GRAV_UNIFORM || GRAV_LINEAR */
#if defined (GRAV_PTMASSPURE) || defined (GRAV_PTMASS) || defined (GRAV_PTMASSSTIFF)
      real    :: r2
#endif /* GRAV_PTMASSPURE || GRAV_PTMASS || GRAV_PTMASSSTIFF */
#if defined (GRAV_PTFLAT) || defined (GRAV_PTMASS)
      real    :: rc, fr
#endif /* GRAV_PTFLAT || GRAV_PTMASS */

      gp_status = ''

#ifdef GRAV_NULL
      call grav_null
#elif defined (GRAV_UNIFORM)
      call grav_uniform
#elif defined (GRAV_LINEAR)
      call grav_linear
#elif defined (GRAV_PTMASS)
      call grav_ptmass_softened(.false.)
#elif defined (GRAV_PTMASSSTIFF)
      call grav_ptmass_stiff
#elif defined (GRAV_PTMASSPURE)
      call grav_ptmass_pure(.false.)
#elif defined (GRAV_PTFLAT)
      call grav_ptmass_softened(.true.)
#elif defined (GRAV_USER)
      call grav_pot_user
#else /* !GRAV_(SPECIFIED) */
      gp_status = 'undefined'
!#warning 'GRAV declared, but gravity model undefined in grav_pot'
! niektore modele grawitacji realizowane sa za pomoca przyspieszenia
! (np. 'galactic') z ktorego liczony jest potencjal
#endif /* !GRAV_(SPECIFIED) */
!-----------------------

      if (gp_status .eq. 'undefined') then
         if (associated(grav_accel)) then
            call warn("[gravity:default_grav_pot_3d]: using 'grav_accel' defined by user")
            call grav_accel2pot
         else
            call die("[gravity:default_grav_pot_3d]: GRAV is defined, but 'gp' is not initialized")
         endif
      endif

   end subroutine default_grav_pot_3d

!>
!! \brief Routine that compute values of gravitational acceleration using gravitational potential array gp
!! \param sweep string of characters that points out the current sweep direction
!! \param i1 integer, number of column in the first direction after one pointed out by sweep
!! \param i2 integer, number of column in the second direction after one pointed out by sweep
!! \param n number of elements of returned array grav
!! \param grav 1D array of gravitational acceleration values computed for positions from xsw and returned by the routine
!<
   subroutine grav_pot2accel(sweep, i1,i2, n, grav,istep)

      use arrays, only: gpot, hgpot
      use grid,   only: dl, xdim, ydim, zdim
      implicit none

      character(len=*), intent(in)   :: sweep
      integer, intent(in)            :: i1, i2
      integer, intent(in)            :: n
      real, dimension(n),intent(out) :: grav
      integer, intent(in)            :: istep
!\todo offer high order gradient as an option in parameter file
!      real, parameter :: onetw = 1./12.

! Gravitational acceleration is computed on right cell boundaries

!      if (istep==1) then
!         select case (sweep)
!            case ('xsweep')
!               grav(3:n-2) = onetw*(hgpot(5:n,i1,i2) - 8.*hgpot(4:n-1,i1,i2) + 8.*hgpot(2:n-3,i1,i2) - hgpot(1:n-4,i1,i2) )/dl(xdim)
!            case ('ysweep')
!               grav(3:n-2) = onetw*(hgpot(i2,5:n,i1) - 8.*hgpot(i2,4:n-1,i1) + 8.*hgpot(i2,2:n-3,i1) - hgpot(i2,1:n-4,i1) )/dl(xdim)
!            case ('zsweep')
!               grav(3:n-2) = onetw*(hgpot(i1,i2,5:n) - 8.*hgpot(i1,i2,4:n-1) + 8.*hgpot(i1,i2,2:n-3) - hgpot(i1,i2,1:n-4) )/dl(xdim)
!         end select
!      else
!         select case (sweep)
!            case ('xsweep')
!               grav(3:n-2) = onetw*(gpot(5:n,i1,i2) - 8.*gpot(4:n-1,i1,i2) + 8.*gpot(2:n-3,i1,i2) - gpot(1:n-4,i1,i2) )/dl(xdim)
!            case ('ysweep')
!               grav(3:n-2) = onetw*(gpot(i2,5:n,i1) - 8.*gpot(i2,4:n-1,i1) + 8.*gpot(i2,2:n-3,i1) - gpot(i2,1:n-4,i1) )/dl(xdim)
!            case ('zsweep')
!               grav(3:n-2) = onetw*(gpot(i1,i2,5:n) - 8.*gpot(i1,i2,4:n-1) + 8.*gpot(i1,i2,2:n-3) - gpot(i1,i2,1:n-4) )/dl(xdim)
!         end select
!      endif
!      grav(2) = grav(3); grav(n-1) = grav(n-2)
!      grav(1) = grav(2); grav(n) = grav(n-1)
      if (istep==1) then
         select case (sweep)
            case ('xsweep')
               grav(2:n-1) = 0.5*(hgpot(1:n-2,i1,i2) - hgpot(3:n,i1,i2))/dl(xdim)
            case ('ysweep')
               grav(2:n-1) = 0.5*(hgpot(i2,1:n-2,i1) - hgpot(i2,3:n,i1))/dl(ydim)
            case ('zsweep')
               grav(2:n-1) = 0.5*(hgpot(i1,i2,1:n-2) - hgpot(i1,i2,3:n))/dl(zdim)
         end select

      else
         select case (sweep)
            case ('xsweep')
               grav(2:n-1) = 0.5*(gpot(1:n-2,i1,i2) - gpot(3:n,i1,i2))/dl(xdim)
            case ('ysweep')
               grav(2:n-1) = 0.5*(gpot(i2,1:n-2,i1) - gpot(i2,3:n,i1))/dl(ydim)
            case ('zsweep')
               grav(2:n-1) = 0.5*(gpot(i1,i2,1:n-2) - gpot(i1,i2,3:n))/dl(zdim)
         end select
      endif

      grav(1) = grav(2); grav(n) = grav(n-1)
   end subroutine grav_pot2accel

!--------------------------------------------------------------------------
!>
!! \brief Routine that uses %gravity acceleration given in grav_accel to compute values of gravitational potential filling in gp array
!<
   subroutine grav_accel2pot

      use arrays,   only: gp
      use grid,     only: dl, xdim, ydim, zdim, is, ie, js, je, ks, ke, nb, nx, ny, nz, zr, yr, xr
      use mpisetup, only: pxsize, pysize, pzsize, pcoords, proc, nproc, ndims, &
           &              comm, comm3d, err, ierr, MPI_DOUBLE_PRECISION, mpifind

      implicit none
      integer               :: i, j, k, ip, pgpmax
      real, allocatable     :: gpwork(:,:,:)
      real                  :: gravrx(nx), gravry(ny), gravrz(nz)
      real                  :: gp_max
      integer, dimension(3) :: loc_gp_max
      integer               :: proc_gp_max
      integer               :: px, py, pz, pc(3)
      real                  :: dgpx_proc, dgpx_all(0:nproc-1), &
                               dgpy_proc, dgpy_all(0:nproc-1), &
                               dgpz_proc, dgpz_all(0:nproc-1), &
                               dgpx(0:pxsize-1,0:pysize-1,0:pzsize-1), &
                               dgpy(0:pxsize-1,0:pysize-1,0:pzsize-1), &
                               dgpz(0:pxsize-1,0:pysize-1,0:pzsize-1), &
                               ddgp(0:pxsize-1,0:pysize-1,0:pzsize-1)


! Gravitational potential gp(i,j,k) is defined in cell centers
! First we find gp(:,:,:) in each block separately assuming:
! Instead of gp, gpwork is used to not change already computed gp when
! disk, bulge and halo parts of galactic potential are computed.
! gravpart is used to distinguish which part has to be computed.

      allocate(gpwork(nx,ny,nz))
      gpwork(1,1,1) = 0.0

      call grav_accel('xsweep',1,1, xr(:), nx, gravrx)
      do i = 1, nx-1
         gpwork(i+1,1,1) = gpwork(i,1,1) - gravrx(i)*dl(xdim)
      enddo

      do i=1,nx
         call grav_accel('ysweep',1, i, yr(:), ny, gravry)
         do j = 1, ny-1
            gpwork(i,j+1,1) = gpwork(i,j,1) - gravry(j)*dl(ydim)
         enddo
      enddo

      do i=1,nx
         do j=1,ny
            call grav_accel('zsweep', i, j, zr(:), nz, gravrz)
            do k = 1, nz-1
               gpwork(i,j,k+1) = gpwork(i,j,k) - gravrz(k)*dl(zdim)
            enddo
         enddo
      enddo

      dgpx_proc = gpwork(is,1,1)-gpwork(1,1,1)
      dgpy_proc = gpwork(1,js,1)-gpwork(1,1,1)
      dgpz_proc = gpwork(1,1,ks)-gpwork(1,1,1)

      call MPI_Gather ( dgpx_proc, 1, MPI_DOUBLE_PRECISION, &
                        dgpx_all,  1, MPI_DOUBLE_PRECISION, &
                        0, comm3d, ierr )
      call MPI_Gather ( dgpy_proc, 1, MPI_DOUBLE_PRECISION, &
                        dgpy_all,  1, MPI_DOUBLE_PRECISION, &
                        0, comm3d, ierr )
      call MPI_Gather ( dgpz_proc, 1, MPI_DOUBLE_PRECISION, &
                        dgpz_all,  1, MPI_DOUBLE_PRECISION, &
                        0, comm3d, ierr )


      if (proc .eq. 0) then

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
         do i = 1, pxsize-1
            ddgp(i,0,0) = ddgp(i-1,0,0) + dgpx(i-1,0,0)
         enddo

         do i=0, pxsize-1
            do j = 1, pysize-1
               ddgp(i,j,0) = ddgp(i,j-1,0) + dgpy(i,j-1,0)
            enddo
         enddo

         do i=0, pxsize-1
            do j=0, pysize-1
               do k = 1, pzsize-1
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

      gp_max      = maxval(gpwork(is:ie,js:je,ks:ke))
      loc_gp_max  = maxloc(gpwork(is:ie,js:je,ks:ke)) &
                  + (/nb,nb,nb/)

      call mpifind(gp_max, 'max', loc_gp_max, proc_gp_max)
      pgpmax = proc_gp_max

      call MPI_Bcast(gp_max, 1, MPI_DOUBLE_PRECISION, pgpmax, comm, ierr)
      gpwork = gpwork - gp_max

      gp=gpwork
      if (allocated(gpwork)) deallocate(gpwork)

   end subroutine grav_accel2pot

end module gravity
