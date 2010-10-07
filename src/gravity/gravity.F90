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
!    Initial implemetation of PIERNIK code was based on TVD split MHD code by
!    Ue-Li Pen
!        see: Pen, Arras & Wong (2003) for algorithm and
!             http://www.cita.utoronto.ca/~pen/MHD
!             for original source code "mhd.f90"
!
!    For full list of developers see $PIERNIK_HOME/license/pdt.txt
!
#include "piernik.def"
#include "macros.h"

!>
!! \brief [DW] Module containing all main subroutines and functions that govern %gravity force in the code
!! \todo to check importance and usefullness of such parameters as g_y and n_gravr2 and tune_zeq_bnd
!!
!!
!! In this module a namelist of parameters is specified:
!! \copydetails gravity::init_grav
!<
module gravity

   character(LEN=9) :: gp_status    !< variable set as 'undefined' in grav_pot_3d when grav_accel is supposed to use
   real    :: g_z                   !< z-component used by GRAV_UNIFORM type of %gravity
   real    :: g_y                   !< y-component of GRAV_UNIFORM constant <b>(currently not used)</b>
   real    :: dg_dz                 !< constant used by GRAV_LINEAR type of %gravity
   real    :: r_gc                  !< galactocentric radius of the local simulation region used by local Galactic type of %gravity in grav_accel
   real    :: ptmass                !< mass value of point %gravity source used by GRAV_PTMASS, GRAV_PTMASSSTIFF, GRAV_PTMASSPURE, GRAV_PTFLAT type of %gravity
   real    :: ptm_x                 !< point mass position x-component
   real    :: ptm_y                 !< point mass position y-component
   real    :: ptm_z                 !< point mass position z-component
   real    :: r_smooth              !< smoothing radius in point mass %types of %gravity
   integer :: nsub                  !< number of subcells while additionally cell division in z-direction is present during estabilishment of hydrostatic equilibrium
   real    :: h_grav                !< altitude of acceleration cut used when n_gravh is set to non-zero
   real    :: r_grav                !< radius of gravitational potential cut used by GRAV_PTMASS, GRAV_PTFLAT type of %gravity
   integer :: n_gravh               !< index of hiperbolic-cosinusoidal cutting of acceleration; used when set to non-zero
   integer :: n_gravr               !< index of hiperbolic-cosinusoidal cutting of gravitational potential used by GRAV_PTMASS, GRAV_PTFLAT type of %gravity
   integer :: n_gravr2              !< similar to n_gravr <b>(currently not used)</b>
   real    :: tune_zeq              !< z-component of %gravity tunning factor used by hydrostatic_zeq
   real    :: tune_zeq_bnd          !< z-component of %gravity tunning factor supposed to be used in boundaries <b>(currently not used)</b>

   logical :: user_grav             !< use user defined grav_pot_3d

   interface
      subroutine user_grav_pot_3d
      end subroutine user_grav_pot_3d
   end interface

   procedure(user_grav_pot_3d), pointer :: grav_pot_3d => NULL()

   contains

!>
!! \brief Routine that sets the initial values of %gravity parameters from namelist @c GRAVITY
!!
!! \n \n
!! @b GRAVITY
!! \n \n
!! <table border="+1">
!! <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!! <tr><td>g_z         </td><td>0.0 </td><td>real             </td><td>\copydoc gravity::g_z          </td></tr>
!! <tr><td>g_y         </td><td>0.0 </td><td>real             </td><td>\copydoc gravity::g_y          </td></tr>
!! <tr><td>dg_dz       </td><td>0.0 </td><td>real             </td><td>\copydoc gravity::dg_dz        </td></tr>
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

      use errh,     only : namelist_errh, warn
      use mpisetup, only : ibuff, rbuff, buffer_dim, comm, ierr, proc, &
           &               MPI_DOUBLE_PRECISION, MPI_INTEGER, MPI_LOGICAL, lbuff
      use arrays,   only : gpot
      use func,     only : compare_namelist
      use dataio_public, only : ierrh, par_file

      implicit none


      namelist /GRAVITY/ g_z,g_y, dg_dz, r_gc, ptmass, ptm_x, ptm_y, ptm_z, r_smooth, &
                nsub, tune_zeq, tune_zeq_bnd, h_grav, r_grav, n_gravr, n_gravr2, n_gravh, user_grav

#ifdef VERBOSE
      call warn("[gravity:init_grav] Commencing gravity module initialization")
#endif /* VERBOSE */

      g_z     = 0.0
      g_y     = 0.0
      dg_dz   = 0.0
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

      user_grav = .false.

      if (proc == 0) then

         diff_nml(GRAVITY)

         ibuff(1)  = nsub
         ibuff(2)  = n_gravr
         ibuff(3)  = n_gravr2
         ibuff(4)  = n_gravh

         rbuff(1)  = g_z
         rbuff(2)  = g_y
         rbuff(3)  = dg_dz
         rbuff(4)  = r_gc
         rbuff(5)  = ptmass
         rbuff(6)  = ptm_x
         rbuff(7)  = ptm_y
         rbuff(8)  = ptm_z
         rbuff(9)  = tune_zeq
         rbuff(10) = tune_zeq_bnd
         rbuff(11) = r_smooth
         rbuff(12) = h_grav
         rbuff(13) = r_grav

         lbuff(1)  = user_grav

      end if

      call MPI_BCAST(ibuff, buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_BCAST(lbuff, buffer_dim, MPI_LOGICAL,          0, comm, ierr)
      call MPI_BCAST(rbuff, buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

      if (proc /= 0) then

         nsub                = ibuff(1)
         n_gravr             = ibuff(2)
         n_gravr2            = ibuff(3)
         n_gravh             = ibuff(4)

         g_z                 = rbuff(1)
         g_y                 = rbuff(2)
         dg_dz               = rbuff(3)
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

         user_grav           = lbuff(1)

      endif

      gpot(:,:,:) = 0.0

      if(.not.user_grav) then
         grav_pot_3d => default_grav_pot_3d
#ifdef VERBOSE
         call warn("[gravity:init_grav] user_grav is set to false. Using default grav_pot_3d.")
#endif /* VERBOSE */
      endif

   end subroutine init_grav

   subroutine source_terms_grav

#if defined(MULTIGRID) || defined(POISSON_FFT)
      use fluidindex,    only : iarr_all_sg
      use arrays,        only : u, sgp, sgpm
#ifdef POISSON_FFT
      use poissonsolver, only : poisson_solve
#endif /* POISSON_FFT */
#ifdef MULTIGRID
      use multigrid,     only : multigrid_solve
#endif /* MULTIGRID */
#endif /* defined(MULTIGRID) || defined(POISSON_FFT) */

      implicit none

#if defined(MULTIGRID) || defined(POISSON_FFT)
      logical, save :: frun = .true.

      sgpm = sgp

#ifdef POISSON_FFT
      call poisson_solve( sum(u(iarr_all_sg,:,:,:),1) )
#endif /* POISSON_FFT */
#ifdef MULTIGRID
      if (size(iarr_all_sg) == 1) then
         call multigrid_solve(u(iarr_all_sg(1),:,:,:))
      else
         call multigrid_solve( sum(u(iarr_all_sg,:,:,:),1) )
         ! BEWARE Here a lot of heap space is required and some compilers may generate code that do segfaults for big enough domains.
         ! It is the weakest point of this type in Maclaurin test. Next one (in fluidboundaries.F90) is 8 times less sensitive.
      end if
#endif /* MULTIGRID */

      ! communicate boundary values for sgp(:, :, :) because multtigrid solver gives at most 2 guardcells, while for hydro solver typically 4 is required.
      call all_sgp_boundaries
      if(frun) then
         sgpm = sgp
         frun = .false.
      endif
#endif /* defined(MULTIGRID) || defined(POISSON_FFT) */

      call sum_potential

   end subroutine source_terms_grav

   subroutine sum_potential

      use mpisetup, only : dt, dtm
      use arrays,   only : gpot, gp, hgpot
#if defined(MULTIGRID) || defined(POISSON_FFT)
      use arrays,   only : sgp, sgpm
#endif /* MULTIGRID || POISSON_FFT */

      implicit none
      real :: h

      if (dtm /= 0) then
         h = dt/dtm
      else
         h = 0.0
      endif

#if defined(MULTIGRID) || defined(POISSON_FFT)
      gpot  = gp + (1.+h)    *sgp -     h*sgpm
      hgpot = gp + (1.+0.5*h)*sgp - 0.5*h*sgpm
#else /* MULTIGRID || POISSON_FFT */
      !BEWARE: as long as grav_pot_3d is called only in init_piernik this assignment probably don't need to be repeated more than once
      gpot  = gp
      hgpot = gp
#endif /* MULTIGRID || POISSON_FFT */

   end subroutine sum_potential

#if defined(MULTIGRID) || defined(POISSON_FFT)

! An improper evaluation of guardcell potential may occur when the multigrid boundary conditions doesn't match /BOUNDARIES/ namelist (e.g. isolated on periodic domain).

   subroutine all_sgp_boundaries

      use mpisetup,   only : comm3d, ierr, MPI_STATUS_SIZE, MPI_REQUEST_NULL, &
           &                 procxl, procxr, procyl, procyr, proczl, proczr, proc, pxsize, pysize, pzsize, &
           &                 bnd_xl, bnd_xr, bnd_yl, bnd_yr, bnd_zl, bnd_zr, &
           &                 ARR_YZ_LEFT_BND, ARR_YZ_RIGHT_BND, ARR_YZ_LEFT_DOM, ARR_YZ_RIGHT_DOM, &
           &                 ARR_XZ_LEFT_BND, ARR_XZ_RIGHT_BND, ARR_XZ_LEFT_DOM, ARR_XZ_RIGHT_DOM, &
           &                 ARR_XY_LEFT_BND, ARR_XY_RIGHT_BND, ARR_XY_LEFT_DOM, ARR_XY_RIGHT_DOM
      use arrays,     only : sgp
      use grid,       only : nb, nxd, nyd, nzd
      use errh,       only : die

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
               end do
            case ('mpi')
               if (pxsize > 1) then
                  call MPI_Isend(sgp(1, 1, 1), 1, ARR_YZ_LEFT_DOM,  procxl, 12, comm3d, req3d(1), ierr)
                  call MPI_Irecv(sgp(1, 1, 1), 1, ARR_YZ_LEFT_BND,  procxl, 22, comm3d, req3d(2), ierr)
               else
                  call die("[gravity:all_grav_boundaries] bnd_xl == 'mpi' && pxsize <= 1")
               end if
            case ('she') ! move appropriate code from poissonsolver::poisson_solve or do nothing. Or die until someone really needs SHEAR.
                call die("[gravity:all_grav_boundaries] bnd_xl == 'she' not implemented")
            case default ! Set gradient == 0 on the boundaries
               do i = 1, nb
                  sgp(i, :, :) = sgp(nb+1, :, :)
               end do
         end select

         select case (bnd_xr)
            case ('per')
               do i = 1, ceiling(nb/real(nxd))
                  sgp(nxd+nb+1:nxd+2*nb, :, :) = sgp(nb+1:2*nb, :, :)
               end do
            case ('mpi')
               if (pxsize > 1) then
                  call MPI_Isend(sgp(1, 1, 1), 1, ARR_YZ_RIGHT_DOM, procxr, 22, comm3d, req3d(3), ierr)
                  call MPI_Irecv(sgp(1, 1, 1), 1, ARR_YZ_RIGHT_BND, procxr, 12, comm3d, req3d(4), ierr)
               else
                  call die("[gravity:all_grav_boundaries] bnd_xr == 'mpi' && pxsize <= 1")
               end if
            case ('she')
                call die("[gravity:all_grav_boundaries] bnd_xl == 'she' not implemented")
            case default
               do i = 1, nb
                  sgp(nb+nxd+i, :, :) = sgp(nb+nxd, :, :)
               end do
         end select

      end if

      if (nyd /= 1) then

         select case (bnd_yl)
            case ('per')
               do i = 1, ceiling(nb/real(nyd))
                  sgp(:, 1:nb, :) = sgp(:, nyd+1:nyd+nb, :)
               end do
            case ('mpi')
               if (pysize > 1) then
                  call MPI_Isend(sgp(1, 1, 1), 1, ARR_XZ_LEFT_DOM,  procyl, 32, comm3d, req3d(5), ierr)
                  call MPI_Irecv(sgp(1, 1, 1), 1, ARR_XZ_LEFT_BND,  procyl, 42, comm3d, req3d(6), ierr)
               else
                  call die("[gravity:all_grav_boundaries] bnd_yl == 'mpi' && pysize <= 1")
               end if
            case default
               do i = 1, nb
                  sgp(:, i, :) = sgp(:, nb+1, :)
               end do
         end select

         select case (bnd_yr)
            case ('per')
               do i = 1, ceiling(nb/real(nyd))
                  sgp(:, nyd+nb+1:nyd+2*nb, :) = sgp(:, nb+1:2*nb, :)
               end do
            case ('mpi')
               if (pysize > 1) then
                  call MPI_Isend(sgp(1, 1, 1), 1, ARR_XZ_RIGHT_DOM, procyr, 42, comm3d, req3d(7), ierr)
                  call MPI_Irecv(sgp(1, 1, 1), 1, ARR_XZ_RIGHT_BND, procyr, 32, comm3d, req3d(8), ierr)
               else
                  call die("[gravity:all_grav_boundaries] bnd_yr == 'mpi' && pysize <= 1")
               end if
            case default
               do i = 1, nb
                  sgp(:, nb+nyd+i, :) = sgp(:, nb+nyd, :)
               end do
         end select

      end if

      if (nzd /= 1) then

         select case (bnd_zl)
            case ('per')
               do i = 1, ceiling(nb/real(nzd))
                  sgp(:, :, 1:nb) = sgp(:, :, nzd+1:nzd+nb)
               end do
            case ('mpi')
               if (pzsize > 1) then
                  call MPI_Isend(sgp(1, 1, 1), 1, ARR_XY_LEFT_DOM,  proczl, 52, comm3d, req3d(9), ierr)
                  call MPI_Irecv(sgp(1, 1, 1), 1, ARR_XY_LEFT_BND,  proczl, 62, comm3d, req3d(10), ierr)
               else
                  call die("[gravity:all_grav_boundaries] bnd_zl == 'mpi' && pzsize <= 1")
               end if
            case default
               do i = 1, nb
                  sgp(:, :, i) = sgp(:, :, nb+1)
               end do
         end select

         select case (bnd_zr)
            case ('per')
               do i = 1, ceiling(nb/real(nzd))
                  sgp(:, :, nzd+nb+1:nzd+2*nb) = sgp(:, :, nb+1:2*nb)
               end do
            case ('mpi')
               if (pzsize > 1) then
                  call MPI_Isend(sgp(1, 1, 1), 1, ARR_XY_RIGHT_DOM, proczr, 62, comm3d, req3d(11), ierr)
                  call MPI_Irecv(sgp(1, 1, 1), 1, ARR_XY_RIGHT_BND, proczr, 52, comm3d, req3d(12), ierr)
               else
                  call die("[gravity:all_grav_boundaries] bnd_zr == 'mpi' && pzsize <= 1")
               end if
            case default
               do i = 1, nb
                  sgp(:, :, nb+nzd+i) = sgp(:, :, nb+nzd)
               end do
         end select

      end if

      call MPI_WAITALL(nreq, req3d(:), status3d(:,:), ierr)

   end subroutine all_sgp_boundaries

#endif /* defined(MULTIGRID) || defined(POISSON_FFT) */

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
!! where \f$r_{soft}\f$ is a radius of softenning\n\n
!! @b GRAV_PTMASSPURE - unsoftened point mass type of %gravity \n
!! \f$\Phi\left(x,y,z\right)= - GM/\sqrt{x^2+y^2+z^2}\f$ \n\n
!! @b GRAV_PTMASSSTIFF - softened point mass type of %gravity with stiff-body rotation inside softening radius\n
!! \f$\Phi\left(x,y,z\right)= - GM/\sqrt{x^2+y^2+z^2}\f$ for \f$r > r_{soft}\f$ and \f$ GM/r_{soft} \left( - 3/2 + 1/2 {x^2+y^2+z^2}/r_{soft}^2 \right)\f$ inside \f$r_{soft}\f$ \n\n
!! @b GRAV_PTFLAT - planar, softened point mass type of %gravity \n
!! \f$\Phi\left(x,y,z\right)= - GM/\sqrt{x^2+y^2+r_{soft}^2}\f$ \n
!! where \f$r_{soft}\f$ is a radius of softenning\n\n
!! @b GRAV_USER - not a standard type of %gravity, implemented by user in the routine grav_pot_user from gravity_user module.\n\n
!<

   subroutine default_grav_pot_3d
      use arrays, only     : gp
#if defined GRAV_PTMASS || defined GRAV_PTMASSPURE || defined GRAV_PTFLAT
      use constants, only  : newtong
#endif /* GRAV_PTMASS || GRAV_PTMASSPURE || GRAV_PTFLAT */
      use grid, only       : nx,ny,nz,x,y,z
      use mpisetup, only   : smalld
      use initfluids, only : cs_iso2
#ifdef GRAV_USER
      use gravity_user, only : grav_pot_user
#endif /* GRAV_USER */
#if defined (GRAV_PTMASSPURE) || defined (GRAV_PTMASS) || defined (GRAV_PTFLAT) || defined (GRAV_PTMASSSTIFF)
      use constants, only : newtong
#endif /* GRAV_PTMASSPURE || GRAV_PTMASS || GRAV_PTFLAT || GRAV_PTMASSSTIFF */
      implicit none

#if defined (GRAV_PTMASSPURE) || defined (GRAV_PTMASS) || defined (GRAV_PTFLAT) || defined (GRAV_PTMASSSTIFF)
      integer :: i, j, k
      real    :: r_smooth2, gm, gmr, z2, yz2
#endif /* GRAV_PTMASSPURE || GRAV_PTMASS || GRAV_PTFLAT || GRAV_PTMASSSTIFF */
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
      gp(:,:,:) = 0.0      ! do nothing

#elif defined (GRAV_UNIFORM)
      do i = 1, nz
         gp(:,:,i) = -g_z*z(i)
      enddo

#elif defined (GRAV_LINEAR)
      do i = 1, nz
         gp(:,:,i) = -0.5 * dg_dz * z(i)**2
      enddo

#elif defined (GRAV_PTMASS)
       do i = 1, nx
          do j = 1, ny
             do k = 1, nz
               rc = dsqrt(x(i)**2+y(j)**2)
               r2 = x(i)**2 + y(j)**2 + z(k)**2
               fr = min( (rc/r_grav)**n_gravr , 100.0)
               fr = max(1./cosh(fr),smalld/100.)
               gp(i,j,k) = -newtong*ptmass / dsqrt(r2 + r_smooth**2)
               gp(i,j,k) = gp(i,j,k) - cs_iso2 * dlog(fr) ! *d0
             enddo
          enddo
       enddo

#elif defined (GRAV_PTMASSSTIFF)

       ! promote stiff-body rotation inside smoothing length, don't affect the global potential outside

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
               end if
             enddo
          enddo
       enddo

#elif defined (GRAV_PTMASSPURE)
       do i = 1, nx
          do j = 1, ny
             do k = 1, nz
               r2 = (x(i) - ptm_x)**2 + (y(j) - ptm_y)**2 + (z(k) - ptm_z)**2
               gp(i,j,k) = -newtong*ptmass / dsqrt(r2 + r_smooth**2)
             enddo
          enddo
       enddo

#elif defined (GRAV_PTFLAT)
       do i = 1, nx
          do j = 1, ny
             rc = dsqrt(x(i)**2+y(j)**2)
             fr = min( (rc/r_grav)**n_gravr , 100.0)
             fr = max(1./cosh(fr),smalld/100.)
             gp(i,j,:) = -newtong*ptmass / dsqrt(rc**2 + r_smooth**2)
             gp(i,j,:) = gp(i,j,:) - cs_iso2 * dlog(fr) ! *d0
          enddo
       enddo

#elif defined (GRAV_USER)
       call grav_pot_user()

#else /* GRAV_(SPECIFIED) */
       gp_status = 'undefined'
!#warning 'GRAV declared, but gravity model undefined in grav_pot'
! niektore modele grawitacji realizowane sa za pomoca przyspieszenia
! (np. 'galactic') z ktorego liczony jest potencjal
#endif /* GRAV_(SPECIFIED) */
!-----------------------

      if(gp_status .eq. 'undefined') then
         call grav_accel2pot
      endif

   end subroutine default_grav_pot_3d

!--------------------------------------------------------------------------
!>
!! \brief Routine that compute values of gravitational acceleration
!! \param sweep string of characters that points out the current sweep direction
!! \param i1 integer, number of column in the first direction after one pointed out by sweep
!! \param i2 integer, number of column in the second direction after one pointed out by sweep
!! \param xsw 1D position array in the direction pointed out by sweep
!! \param n number of elements of xsw array
!! \param grav 1D array of gravitational acceleration values computed for positions from xsw and returned by the routine
!! \n\n
!! one type of %gravity is implemented here: \n\n
!! local Galactic %gravity only in z-direction (see <a href="http://cdsads.u-strasbg.fr/abs/1998ApJ...497..759F">Ferriere K., 1998, Astrophys. Journal, 497, 759</a>)\n
!! \f[
!! F_z = 3.23 \cdot 10^8 \cdot \left[\left(-4.4 \cdot 10^{-9} \cdot exp\left(-\frac{(r_{gc}-r_{gc_{}Sun})}{(4.9kpc)}\right) \cdot \frac{z}{\sqrt{(z^2+(0.2kpc)^2)}}\right)
!! -\left( 1.7 \cdot 10^{-9} \cdot \frac{(r_{gc_{}Sun}^2 + (2.2kpc)^2)}{(r_{gc}^2 + (2.2kpc)^2)} \cdot \frac{z}{1kpc}\right) \right]
!! \f]
!! where \f$r_{gc}\f$ is galactocentric radius and \f$r_{gcSun}\f$ is the galactocentric radius of Sun.
!<

   subroutine grav_accel(sweep, i1,i2, xsw, n, grav)
      use grid, only     : nb
      use arrays, only   : gp
      use grid, only     : x,y,z,dl,xdim,ydim,zdim,nx,ny,nz
      use grid, only     : is,ie,js,je,ks,ke, xr,yr,zr
#if defined GRAV_PTMASS || defined GRAV_PTFLAT
      use mpisetup, only : smalld
#endif /* GRAV_PTMASS || GRAV_PTFLAT */
#if defined GRAV_ACC_USER
      use gravity_user, only : grav_accel_user
#endif /* GRAV_ACC_USER  */
      implicit none
      character(len=6), intent(in)   :: sweep
      integer, intent(in)            :: i1, i2
      integer, intent(in)            :: n
      real, dimension(:)             :: xsw
      real, dimension(n),intent(out) :: grav
      real                           :: x1, x2
#if defined GRAV_PTMASS || defined GRAV_PTFLAT
      real, dimension(n)             :: r, r32
#endif /* GRAV_PTMASS || GRAV_PTFLAT */

      select case (sweep)
         case('xsweep')
            x1  = y(i1)
            x2  = z(i2)
         case('ysweep')
            x1  = z(i1)
            x2  = x(i2)
         case('zsweep')
            x1  = x(i1)
            x2  = y(i2)
         case default ! just for suppressing compiler warning
            x1 = xsw(1)
      end select

      if (sweep == 'zsweep') then
!         grav = 0.1*(tanh(abs((2.0*xsw)**10.0)-30.0)+1.0)*0.5
         grav = - 0.1
      else
         grav=0.0
      endif

   end subroutine grav_accel

!>
!! \brief Routine that compute values of gravitational acceleration using gravitational potential array gp
!! \param sweep string of characters that points out the current sweep direction
!! \param i1 integer, number of column in the first direction after one pointed out by sweep
!! \param i2 integer, number of column in the second direction after one pointed out by sweep
!! \param n number of elements of returned array grav
!! \param grav 1D array of gravitational acceleration values computed for positions from xsw and returned by the routine
!<
   subroutine grav_pot2accel(sweep, i1,i2, n, grav,istep)
      use arrays, only : gpot, hgpot
      use grid, only   : dl, xdim, ydim, zdim
      implicit none
      character(len=6), intent(in)   :: sweep
      integer, intent(in)            :: i1, i2
      integer, intent(in)            :: n
      real, dimension(n),intent(out) :: grav
      integer, intent(in)            :: istep
!\todo offer high order gradient as an option in parameter file
!      real, parameter :: onetw = 1./12.

! Gravitational accelleration is computed on right cell boundaries

!      if(istep==1) then
!         select case(sweep)
!            case('xsweep')
!               grav(3:n-2) = onetw*(hgpot(5:n,i1,i2) - 8.*hgpot(4:n-1,i1,i2) + 8.*hgpot(2:n-3,i1,i2) - hgpot(1:n-4,i1,i2) )/dl(xdim)
!            case('ysweep')
!               grav(3:n-2) = onetw*(hgpot(i2,5:n,i1) - 8.*hgpot(i2,4:n-1,i1) + 8.*hgpot(i2,2:n-3,i1) - hgpot(i2,1:n-4,i1) )/dl(xdim)
!            case('zsweep')
!               grav(3:n-2) = onetw*(hgpot(i1,i2,5:n) - 8.*hgpot(i1,i2,4:n-1) + 8.*hgpot(i1,i2,2:n-3) - hgpot(i1,i2,1:n-4) )/dl(xdim)
!         end select
!      else
!         select case(sweep)
!            case('xsweep')
!               grav(3:n-2) = onetw*(gpot(5:n,i1,i2) - 8.*gpot(4:n-1,i1,i2) + 8.*gpot(2:n-3,i1,i2) - gpot(1:n-4,i1,i2) )/dl(xdim)
!            case('ysweep')
!               grav(3:n-2) = onetw*(gpot(i2,5:n,i1) - 8.*gpot(i2,4:n-1,i1) + 8.*gpot(i2,2:n-3,i1) - gpot(i2,1:n-4,i1) )/dl(xdim)
!            case('zsweep')
!               grav(3:n-2) = onetw*(gpot(i1,i2,5:n) - 8.*gpot(i1,i2,4:n-1) + 8.*gpot(i1,i2,2:n-3) - gpot(i1,i2,1:n-4) )/dl(xdim)
!         end select
!      endif
!      grav(2) = grav(3); grav(n-1) = grav(n-2)
!      grav(1) = grav(2); grav(n) = grav(n-1)
      if(istep==1) then
         select case(sweep)
            case('xsweep')
               grav(2:n-1) = 0.5*(hgpot(1:n-2,i1,i2) - hgpot(3:n,i1,i2))/dl(xdim)
            case('ysweep')
               grav(2:n-1) = 0.5*(hgpot(i2,1:n-2,i1) - hgpot(i2,3:n,i1))/dl(ydim)
            case('zsweep')
               grav(2:n-1) = 0.5*(hgpot(i1,i2,1:n-2) - hgpot(i1,i2,3:n))/dl(zdim)
         end select

      else
         select case(sweep)
            case('xsweep')
               grav(2:n-1) = 0.5*(gpot(1:n-2,i1,i2) - gpot(3:n,i1,i2))/dl(xdim)
            case('ysweep')
               grav(2:n-1) = 0.5*(gpot(i2,1:n-2,i1) - gpot(i2,3:n,i1))/dl(ydim)
            case('zsweep')
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

      use mpisetup, only : pxsize, pysize, pzsize, pcoords, proc, nproc, ndims, &
           &               comm, comm3d, err, ierr, MPI_DOUBLE_PRECISION, mpifind
      use arrays,   only : gp
      use grid,     only : dl, xdim, ydim, zdim, is, ie, js, je, ks, ke, nb, nx, ny, nz, zr, yr, xr

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

      call MPI_GATHER ( dgpx_proc, 1, MPI_DOUBLE_PRECISION, &
                        dgpx_all,  1, MPI_DOUBLE_PRECISION, &
                        0, comm3d, ierr )
      call MPI_GATHER ( dgpy_proc, 1, MPI_DOUBLE_PRECISION, &
                        dgpy_all,  1, MPI_DOUBLE_PRECISION, &
                        0, comm3d, ierr )
      call MPI_GATHER ( dgpz_proc, 1, MPI_DOUBLE_PRECISION, &
                        dgpz_all,  1, MPI_DOUBLE_PRECISION, &
                        0, comm3d, ierr )


      if(proc .eq. 0) then

         do ip = 0, nproc-1
            call MPI_CART_COORDS(comm3d, ip, ndims, pc, ierr)

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

      call MPI_BCAST(ddgp, nproc, MPI_DOUBLE_PRECISION, 0, comm, ierr)

      px = pcoords(1)
      py = pcoords(2)
      pz = pcoords(3)

      gpwork = gpwork + ddgp(px,py,pz)

      gp_max      = maxval(gpwork(is:ie,js:je,ks:ke))
      loc_gp_max  = maxloc(gpwork(is:ie,js:je,ks:ke)) &
                  + (/nb,nb,nb/)

      call mpifind(gp_max, 'max', loc_gp_max, proc_gp_max)
      pgpmax = proc_gp_max

      call MPI_BCAST(gp_max, 1, MPI_DOUBLE_PRECISION, pgpmax, comm, ierr)
      gpwork = gpwork - gp_max

      gp=gpwork
      if(allocated(gpwork)) deallocate(gpwork)

   end subroutine grav_accel2pot


end module gravity
