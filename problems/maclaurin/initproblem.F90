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
!! \brief Implementation of selfgravitating distributions, for which an analytical solution of the gravitational potential is known.
!!
!! \details The analytical cases implemented so far:
!! * Maclaurin spheroid with a1 as equatorial radius and e as polar eccentricity of the spheroid (e>0 gives oblate object, e<0 gives prolate object)
!! * Plummersphere with a1 as the Plummer radius and e ignored
!<

module initproblem

   use constants, only: dsetnamelen, cbuff_len

   implicit none

   private
   public :: read_problem_par, problem_initial_conditions, problem_pointers

   ! namelist parameters
   real :: x0                         !< X-position of the center
   real :: y0                         !< Y-position of the center
   real :: z0                         !< Z-position of the center
   real :: d0                         !< central density
   real :: a1                         !< equatorial radius of the Maclaurin spheroid, Plummer radius, ...
   real :: e                          !< polar eccentricity of the spheroid (e>0 gives oblate object, e<0 gives prolate object), ignored for Plummer
   real :: ref_thr                    !< refinement threshold
   real :: ref_eps                    !< smoother filter
   integer(kind=4) :: nsub            !< subsampling on the grid
   logical :: analytical_ext_pot      !< if .true. then bypass multipole solver and use analytical potential for external boundaries (debugging/developing only)
   character(len=cbuff_len) :: model  !< one of ("Maclaurin", "Plummer")

   namelist /PROBLEM_CONTROL/ x0, y0, z0, d0, a1, e, ref_thr, ref_eps, nsub, analytical_ext_pot, model

   ! private data
   real :: d1 !< ambient density
   real :: p0 !< pressure
   real :: a3 !< length of polar radius of the spheroid
   character(len=dsetnamelen), parameter :: apot_n = "apot"   !< name of the analytical potential field
   character(len=dsetnamelen), parameter :: asrc_n = "asrc"   !< name of the source field used for "ares" calculation (auxiliary space)
   character(len=dsetnamelen), parameter :: ares_n = "ares"   !< name of the numerical residuum with respect to analytical potential field
   character(len=dsetnamelen), parameter :: mpole_n = "mpole" !< name of the potential recovered solely from multipole moments
   integer :: i_model

   enum, bind(C)
      enumerator :: MACLAURIN, PLUMMER
   end enum

contains

!> \brief Set up custom pointers to tweak the code execution according to our needs

   subroutine problem_pointers

      use user_hooks,  only: finalize_problem, problem_post_IC
#ifdef HDF5
      use dataio_user, only: user_vars_hdf5, user_attrs_wr
#endif /* HDF5 */

      implicit none

      problem_post_IC  => compute_mpole
      finalize_problem => finalize_problem_maclaurin
#ifdef HDF5
      user_attrs_wr    => problem_initial_conditions_attrs
      user_vars_hdf5   => maclaurin_error_vars
#endif /* HDF5 */

   end subroutine problem_pointers

!> \brief Read the runtime parameters specified in the namelist and set up some auxiliary values

   subroutine read_problem_par

      use cg_list_global,        only: all_cg
      use constants,             only: pi, GEO_XYZ, GEO_RPZ, xdim, ydim, LO, HI, cbuff_len, INVALID, V_VERBOSE, V_INFO
      use dataio_pub,            only: die, warn, msg, printinfo, nh
      use domain,                only: dom
      use fluidindex,            only: iarr_all_dn
      use global,                only: smalld
      use mpisetup,              only: rbuff, ibuff, lbuff, cbuff, master, slave, piernik_MPI_Bcast
      use multigridvars,         only: ord_prolong
      use named_array_list,      only: wna
      use unified_ref_crit_list, only: urc_list
      use user_hooks,            only: ext_bnd_potential
#ifdef NBODY
      use constants,             only: I_ONE
      use particle_utils,        only: add_part_in_proper_cg
#endif /* NBODY */

      implicit none

      integer, parameter :: maxsub = 10  !< upper limit for subsampling
      integer(kind=4) :: id

      d1 = smalld                  ! ambient density

      ! namelist default parameter values
      x0           = 0.0                 !< x-coordinate of the spheroid center
      y0           = 0.0                 !< y-coordinate of the spheroid center
      z0           = 0.0                 !< z-coordinate of the spheroid center
      d0           = 1.0                 !< Density inside the sphere
      a1           = 1.0                 !< Equatorial semimajor axis
      e            = 0.0                 !< Eccentricity; e>0 for flattened spheroids, e<0 for elongated spheroids
      nsub         = 3                   !< Subsampling factor

      ref_thr      = 0.3    !< Refine if density difference is greater than this value
      ref_eps      = 0.01   !< refinement smoothing factor

      analytical_ext_pot = .false.

      model = "Maclaurin"

      if (master) then

         if (.not.nh%initialized) call nh%init()
         open(newunit=nh%lun, file=nh%tmp1, status="unknown")
         write(nh%lun,nml=PROBLEM_CONTROL)
         close(nh%lun)
         open(newunit=nh%lun, file=nh%par_file)
         nh%errstr=""
         read(unit=nh%lun, nml=PROBLEM_CONTROL, iostat=nh%ierrh, iomsg=nh%errstr)
         close(nh%lun)
         call nh%namelist_errh(nh%ierrh, "PROBLEM_CONTROL")
         read(nh%cmdl_nml,nml=PROBLEM_CONTROL, iostat=nh%ierrh)
         call nh%namelist_errh(nh%ierrh, "PROBLEM_CONTROL", .true.)
         open(newunit=nh%lun, file=nh%tmp2, status="unknown")
         write(nh%lun,nml=PROBLEM_CONTROL)
         close(nh%lun)
         call nh%compare_namelist()

         rbuff(1) = x0
         rbuff(2) = y0
         rbuff(3) = z0
         rbuff(4) = d0
         rbuff(5) = a1
         rbuff(6) = e
         rbuff(7) = ref_thr
         rbuff(9) = ref_eps

         ibuff(1) = nsub

         lbuff(1) = analytical_ext_pot

         cbuff(1) = model

      endif

      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(rbuff)
      call piernik_MPI_Bcast(lbuff)
      call piernik_MPI_Bcast(cbuff, cbuff_len)

      if (slave) then

         x0           = rbuff(1)
         y0           = rbuff(2)
         z0           = rbuff(3)
         d0           = rbuff(4)
         a1           = rbuff(5)
         e            = rbuff(6)
         ref_thr      = rbuff(7)
         ref_eps      = rbuff(9)

         nsub         = ibuff(1)

         analytical_ext_pot = lbuff(1)

         model        = cbuff(1)

      endif

      select case (trim(model))
         case ("Maclaurin", "maclaurin", "ml")
            i_model = MACLAURIN
         case ("Plummer", "plummer")
            i_model = PLUMMER
         case default
            i_model = INVALID
      end select
      if (i_model == INVALID) call die("[initproblem:read_problem_par] unrecognized potential model '" // trim(model) // "'")

      if (analytical_ext_pot .and. i_model == MACLAURIN) ext_bnd_potential => maclaurin2bnd_potential
      ! ToDo: implement plummer2bnd_potential

      if (a1 <= 0.) then ! point-like source
         a1 = 0.
         e = 0.
      endif

      if (i_model == MACLAURIN) then
         if (abs(e) >= 1.) call die("[initproblem:read_problem_par] |e|>=1.")
         if (e >= 0.) then            ! vertical axis
            a3 = a1 * sqrt(1. - e**2) ! oblate Maclaurin spheroid
         else
            a3 = a1 / sqrt(1. - e**2) ! prolate spheroid
         endif
      else
         a3 = a1
      endif
      p0 = 100.*tiny(d1)           ! pressure

      if (d0 < 0.) call die("[initproblem:read_problem_par] Negative average density.")

      if (nsub < 1) then
         if (master) call warn("[initproblem:read_problem_par] subsampling disabled.")
         nsub = 1
      else if (nsub > maxsub) then
         if (master)call warn("[initproblem:read_problem_par] too much subsampling.")
         nsub = maxsub
      endif

#ifdef NBODY
      if (abs(a1) < tiny(1.)) call add_part_in_proper_cg(I_ONE, d0, [ x0, y0, z0 ], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 0.0)
#endif /* NBODY */

      if (master) then
         if (a1 > 0.) then
            select case (i_model)
               case (MACLAURIN)
                  write(msg, '(3(a,g12.5),a)')"[initproblem:problem_initial_conditions] Set up Maclaurin spheroid with a1 and a3 axes = ", a1, ", ", a3, " (eccentricity = ", e, ")"
               case (PLUMMER)
                  write(msg, '(a,g12.5)')"[initproblem:problem_initial_conditions] Set up Plummer sphere with a = ", a1
               case default
                  write(msg, '(a)')"[initproblem:problem_initial_conditions] No idea what's going on"
            end select
         else
            write(msg, '(a)')"[initproblem:problem_initial_conditions] Set up point-like mass"
         endif
         call printinfo(msg, V_INFO)
         if (x0-a1<dom%edge(xdim, LO) .or. x0+a1>dom%edge(xdim, HI)) call warn("[initproblem:problem_initial_conditions] Part of the spheroid is outside the domain in the X-direction.")
         if ( (dom%geometry_type == GEO_XYZ .and. (y0-a1<dom%edge(ydim, LO) .or. y0+a1>dom%edge(ydim, HI))) .or. &
              (dom%geometry_type == GEO_RPZ .and. (atan2(a1,x0) > minval([y0-dom%edge(ydim, LO), dom%edge(ydim, HI)-y0]))) ) & ! will fail when some one adds 2*k*pi to y0
              call warn("[initproblem:problem_initial_conditions] Part of the spheroid is outside the domain")
         if (a1 > 0.) then
            write(msg,'(2(a,g12.5))')   "[initproblem:problem_initial_conditions] Density = ", d0, " mass = ", 4./3.*pi * a1**2 * a3 * d0
         else
            write(msg,'(a,g12.5)')   "[initproblem:problem_initial_conditions] Mass = ", d0
         endif
         call printinfo(msg, V_VERBOSE)
      endif

      call all_cg%reg_var(apot_n, ord_prolong = ord_prolong)
      call all_cg%reg_var(ares_n)
      call all_cg%reg_var(asrc_n)
      call all_cg%reg_var(mpole_n)

      ! Set up automatic refinement criteria on densities
      do id = lbound(iarr_all_dn, dim=1, kind=4), ubound(iarr_all_dn, dim=1, kind=4)
         !> \warning only selfgravitating fluids should be added
         select case (i_model)
            case (MACLAURIN)
               call urc_list%add_user_urcv(wna%fi, id, ref_thr, ref_eps, "Loechner", .true.)
            case (PLUMMER)
               call urc_list%add_user_urcv(wna%fi, id, ref_thr*d0, 0., "grad", .true.)
         end select
      enddo

   end subroutine read_problem_par

!> \brief Set up the initial conditions. Note that this routine can be called multiple times during initial iterations of refinement structure

   subroutine problem_initial_conditions

      use cg_cost_data,      only: I_IC
      use cg_leaves,         only: leaves
      use cg_list,           only: cg_list_element
      use constants,         only: GEO_XYZ, GEO_RPZ, xdim, ydim, zdim, LO, HI, V_DEBUG
      use dataio_pub,        only: die, msg, printinfo
      use domain,            only: dom
      use fluidindex,        only: iarr_all_dn, iarr_all_mx, iarr_all_my, iarr_all_mz
      use global,            only: dirty_debug, no_dirty_checks
      use grid_cont,         only: grid_container
      use named_array_list,  only: qna
      use mpisetup,          only: master
      use multigrid_Laplace, only: residual
      use units,             only: fpiG

      implicit none

      integer                        :: i, j, k, ii, jj, kk
      real                           :: xx, yy, zz, rr, dm, n_res
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg

      call compute_analytical_potential

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         call cg%costs%start

         if (a1 > 0.) then
            do k = cg%ks, cg%ke
               do j = cg%js, cg%je
                  do i = cg%is, cg%ie

                     !< \todo use subsampling only near the surface of the spheroid
                     dm = 0.
                     do kk = -nsub+1, nsub-1, 2
                        zz = ((cg%z(k) + kk*cg%dz/(2.*nsub) - z0)/a3)**2
                        do jj = -nsub+1, nsub-1, 2
                           do ii = -nsub+1, nsub-1, 2

                              select case (dom%geometry_type)
                                 case (GEO_XYZ)
                                    yy = ((cg%y(j) + jj*cg%dy/(2.*nsub) - y0)/a1)**2
                                    xx = ((cg%x(i) + ii*cg%dx/(2.*nsub) - x0)/a1)**2
                                    rr = xx + yy + zz
                                 case (GEO_RPZ)
                                    yy = cg%y(j) + jj*cg%dy/(2.*nsub) - y0
                                    xx = cg%x(i) + ii*cg%dx/(2.*nsub)
                                    rr = (xx**2 + x0**2 - 2. * xx * x0 * cos(yy))/a1**2 + zz
                                 case default
                                    call die("[initproblem:problem_initial_conditions] Unsupported dom%geometry_type")
                                    rr = 0.
                              end select

                              select case (i_model)
                                 case (MACLAURIN)
                                    if (rr <= 1.) then
                                       dm = dm + d0
                                    else
                                       dm = dm + d1
                                    endif
                                 case (PLUMMER)
                                    dm = dm + d0 * (1. + rr) ** (-5./2.)
                              end select

                           enddo
                        enddo
                     enddo
                     cg%u(iarr_all_dn, i, j, k) = dm / nsub**3
                     !> \warning What if more than one fluid is in use?
                     !> \warning What if some of the fluids are not selfgravitating?

                  enddo
               enddo
            enddo
         else
            cg%u(iarr_all_dn, RNG) = 0.0
         endif

         cg%u(iarr_all_mx, RNG) = 0.0
         cg%u(iarr_all_my, RNG) = 0.0
         cg%u(iarr_all_mz, RNG) = 0.0

#ifdef MAGNETIC
         call cg%set_constant_b_field([0., 0., 0.])
#endif /* MAGNETIC */

         call cg%costs%stop(I_IC)
         cgl => cgl%nxt
      enddo

      ! Compute residuum for analytical solution and print its norm
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         if (dirty_debug .and. .not. no_dirty_checks) then
            cg%wa(RNG) = cg%q(qna%ind(apot_n))%arr(RNG)
            cg%q(qna%ind(apot_n))%arr = huge(1.)
            cg%q(qna%ind(apot_n))%arr(RNG) = cg%wa(RNG)
            cg%q(qna%ind(asrc_n))%arr = huge(1.)
         endif
         cg%q(qna%ind(asrc_n))%arr(RNG) = fpiG * sum(cg%u(iarr_all_dn, RNG), dim=1)
         cgl => cgl%nxt
      enddo

      call residual(leaves, qna%ind(asrc_n), qna%ind(apot_n), qna%ind(ares_n))
      call leaves%check_dirty(qna%ind(ares_n), "a-residual")

      ! Clear residual next to the external boundary as it is affected by the way the potential is extrapolated into the guardcells.
      ! When ord_prolong /= 0, there still will remain some artifacts at fine/coarse boundaries touching external boundary
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         do i = xdim, zdim
            do j = LO, HI
               if (cg%ext_bnd(i, j)) then
                  select case (j + HI*i)
                     case (LO + HI*xdim)
                        cg%q(qna%ind(ares_n))%arr(cg%is, :, :) = 0.
                     case (HI + HI*xdim)
                        cg%q(qna%ind(ares_n))%arr(cg%ie, :, :) = 0.
                     case (LO + HI*ydim)
                        cg%q(qna%ind(ares_n))%arr(:, cg%js, :) = 0.
                     case (HI + HI*ydim)
                        cg%q(qna%ind(ares_n))%arr(:, cg%je, :) = 0.
                     case (LO + HI*zdim)
                        cg%q(qna%ind(ares_n))%arr(:, :, cg%ks) = 0.
                     case (HI + HI*zdim)
                        cg%q(qna%ind(ares_n))%arr(:, :, cg%ke) = 0.
                     case default
                        call die("[initproblem:problem_initial_conditions] Non-existing side.")
                  end select
               endif
            enddo
         enddo
         cgl => cgl%nxt
      enddo

      dm = leaves%norm_sq(qna%ind(asrc_n))
      n_res = leaves%norm_sq(qna%ind(ares_n))
      if (abs(dm) > tiny(1.)) then  ! An FP overflow will occur when n_res > dm/tiny(1.)
         write(msg, '(a,f13.10)')"[initproblem:problem_initial_conditions] Analytical norm residual/source= ", n_res/dm
         if (master) call printinfo(msg, V_DEBUG)
      else
         write(msg, '(2(a,g14.6))')"[initproblem:problem_initial_conditions] Analytical norm residual= ", n_res, " point mass= ", d0
         ! Is n_res ~ sqrt(cg%dvol) ?
         if (master) call printinfo(msg, V_DEBUG)
      endif

   end subroutine problem_initial_conditions

!> \brief Provides parameters useful for python/maclaurin.py in .h5 files

#ifdef HDF5
   subroutine problem_initial_conditions_attrs(file_id)

      use hdf5,  only: HID_T, SIZE_T
      use h5lt,  only: h5ltset_attribute_double_f
      use units, only: fpiG

      implicit none

      integer(HID_T),intent(in)  :: file_id

      integer(SIZE_T) :: bufsize
      integer(kind=4) :: error

      bufsize = 1
      call h5ltset_attribute_double_f(file_id, "/", "rho0", [d0],   bufsize, error)
      call h5ltset_attribute_double_f(file_id, "/", "fpiG", [fpiG], bufsize, error)
      call h5ltset_attribute_double_f(file_id, "/", "a1",   [a1],   bufsize, error)
      call h5ltset_attribute_double_f(file_id, "/", "x0",   [x0],   bufsize, error)
      call h5ltset_attribute_double_f(file_id, "/", "y0",   [y0],   bufsize, error)
      call h5ltset_attribute_double_f(file_id, "/", "z0",   [z0],   bufsize, error)

   end subroutine problem_initial_conditions_attrs
#endif /* HDF5 */

!> \brief wrapper for analytical models

   subroutine compute_analytical_potential

      implicit none

      select case (i_model)
         case (MACLAURIN)
            call compute_maclaurin_potential
         case (PLUMMER)
            call compute_plummer_potential
      end select

   contains

      subroutine compute_plummer_potential

         use cg_leaves,        only: leaves
         use cg_list,          only: cg_list_element
         use constants,        only: pi, GEO_XYZ, GEO_RPZ, PPP_PROB
         use dataio_pub,       only: die
         use domain,           only: dom
         use grid_cont,        only: grid_container
         use named_array_list, only: qna
         use ppp,              only: ppp_main
         use units,            only: newtong

         implicit none

         integer                        :: i, j, k, apot_i
         real                           :: r2
         real                           :: x02, y02, z02, cdphi
         type(cg_list_element), pointer :: cgl
         type(grid_container),  pointer :: cg
         character(len=*), parameter    :: cpp_label = "compute_plummer_potential"

         call ppp_main%start(cpp_label, PPP_PROB)

         apot_i = qna%ind(apot_n)

         cgl => leaves%first
         do while (associated(cgl))
            cg => cgl%cg

            do k = cg%ks, cg%ke
               z02 = (cg%z(k)-z0)**2
               do j = cg%js, cg%je
                  do i = cg%is, cg%ie
                     select case (dom%geometry_type)
                        case (GEO_XYZ)
                           y02 = (cg%y(j)-y0)**2
                           x02 = (cg%x(i)-x0)**2
                           r2 = x02 + y02 + z02
                        case (GEO_RPZ)
                           cdphi = cos(cg%y(j)-y0)
                           r2 = (cg%x(i)**2 + x0**2 - 2. * cg%x(i) * x0 * cdphi) + z02
                           x02 = r2 * cdphi**2
                           y02 = r2 - x02
                        case default
                           call die("[initproblem:compute_plummer_potential] Invalid dom%geometry_type.")
                           r2 = 0 ; x02 = 0 ; y02 = 0 ! suppress compiler warnings
                     end select

                     cg%q(apot_i)%arr(i, j, k) = - 4./3. * pi * newtong * d0 * a1**3 / sqrt(r2 + a1**2)

                  enddo
               enddo
            enddo

            cgl => cgl%nxt
         enddo

         call ppp_main%stop(cpp_label, PPP_PROB)

      end subroutine compute_plummer_potential

   end subroutine compute_analytical_potential

!>
!! \brief Calculate analytical potential for given Maclaurin spheroid.
!!
!! \details We assume that the spheroid is fully contained in the computational domain here.
!! Oblate potential formula from Ricker_2008ApJS..176..293R
!! Prolate potential formula from http://scienceworld.wolfram.com/physics/ProlateSpheroidGravitationalPotential.html
!<
   subroutine compute_maclaurin_potential

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: pi, GEO_XYZ, GEO_RPZ, PPP_PROB
      use dataio_pub,       only: warn, die
      use domain,           only: dom
      use grid_cont,        only: grid_container
      use mpisetup,         only: master
      use named_array_list, only: qna
      use ppp,              only: ppp_main
      use units,            only: newtong

      implicit none

      integer                        :: i, j, k, apot_i
      real                           :: potential, r2, rr
      real                           :: AA1, AA3, a12, a32, x02, y02, z02, lam, h, cdphi
      real, parameter                :: small_e = 1e-3
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg
      character(len=*), parameter    :: cmp_label = "compute_maclaurin_potential"

      call ppp_main%start(cmp_label, PPP_PROB)

      AA1 = 2./3. ; AA3 = 2./3.
      if (e < 0. .and. master) call warn("[initproblem:compute_maclaurin_potential] e<0. not fully implemented yet!")

      if (e > small_e) then
         AA1 = ( sqrt(1. - e**2) * asin(e) - e * (1. - e**2) ) / e**3
         AA3 = 2. * ( e - sqrt(1. - e**2) * asin(e) ) / e**3
      else if (e>0.) then            ! Taylor expansions
         AA1 = 2./3. - 2./15. * e**2
         AA3 = 2./3. + 4./15. * e**2
      else if (e<0.) then            ! ToDo: find analytical expressions for -e > small_e
         AA1 = 2./3. + 2./15. * e**2
         AA3 = 2./3. - 4./15. * e**2
      endif

      a12 = a1**2
      a32 = a3**2

      apot_i = qna%ind(apot_n)

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         do k = cg%ks, cg%ke
            z02 = (cg%z(k)-z0)**2
            do j = cg%js, cg%je
               do i = cg%is, cg%ie
                  if (abs(a12) > tiny(1.)) then
                     select case (dom%geometry_type)
                        case (GEO_XYZ)
                           y02 = (cg%y(j)-y0)**2
                           x02 = (cg%x(i)-x0)**2
                           r2 = (x02+y02)/a12 + z02/a32
                        case (GEO_RPZ)
                           cdphi = cos(cg%y(j)-y0)
                           r2 = (cg%x(i)**2 + x0**2 - 2. * cg%x(i) * x0 * cdphi)/a12 + z02/a32
                           x02 = r2 * cdphi**2
                           y02 = r2 - x02
                        case default
                           call die("[initproblem:compute_maclaurin_potential] Invalid dom%geometry_type.")
                           r2 = 0 ; x02 = 0 ; y02 = 0 ! suppress compiler warnings
                     end select
                     rr = r2 * a12
                  else
                     e = 0.
                     r2 = huge(1.) ; x02 = 0. ; y02 = 0. ; rr = 0. ! suppress compiler warnings
                     select case (dom%geometry_type)
                        case (GEO_XYZ)
                           y02 = (cg%y(j)-y0)**2
                           x02 = (cg%x(i)-x0)**2
                           rr = x02 + y02 + z02
                        case (GEO_RPZ)
                           cdphi = cos(cg%y(j)-y0)
                           rr = (cg%x(i)**2 + x0**2 - 2. * cg%x(i) * x0 * cos(cg%y(j)-y0)) + z02
                        case default
                           call die("[initproblem:compute_maclaurin_potential] Invalid dom%geometry_type.")
                     end select
                  endif
                  if (r2 > 1.) then
                     if (e > small_e) then
                        lam = 0.5 * (a12 + a32 - (x02 + y02 + z02))
                        lam = -lam + sqrt(lam**2 + a12 * z02 + a32 * (x02 + y02 - a12))
                        h = a1 * e / sqrt(a32 + lam)
                        ! for e < small_e the expressions (atan(h) - h / (1. + h**2)) and (h - atan(h)) should be replaced with Taylor expansions
                        potential = - 2. * a1 * a3 / e * (atan(h) - ((x02 + y02) * (atan(h) - h / (1. + h**2)) + 2. * z02 * (h - atan(h)))/(2. * a12 * e**2))
                     else if (e < -small_e) then
                        lam = 0.5 * (a12 + a32 - (x02 + y02 + z02))
                        lam = - lam + sqrt(lam**2 + a12 * z02 + a32 * (x02 + y02 - a12))
                        h = sqrt((a32 - a12)/(a12 + lam))
                        potential = - a12 * a3 / (a32 - a12) * ( &
                             &      (2.*(a32 - a12) + x02 + y02 - 2.*z02)/sqrt(a32 - a12) * log(h + sqrt(1. + h**2)) - &
                             &      (x02 + y02) * sqrt(a32 + lam)/(a12 + lam) + 2.*z02 / sqrt(a32 + lam) )
                     else
                        if (abs(a1) < tiny(1.)) then
                           if (abs(rr) > sum(cg%dl**2)/8.) then  ! Warning: may need some small factor for better smoothness
                              potential = - 1. / (pi * sqrt(rr)) ! A bit dirty: we interpret d0 as a mass for point-like sources
                           else
                              potential = - 8. / (pi * sqrt(sum(cg%dl**2))) *.85
                              ! The factor of 0.85 was guessed via minimizing L2 error norm for point-like setups
                           endif
                        else
                           potential = - 4./3. * a1**3 / sqrt(rr)
                        endif
                     endif
                  else
                     if (abs(e) > small_e) then
                        potential = - (AA1*(2*a12 - x02 - y02) + AA3 * (a32 - z02))
                     else
                        potential = - 2./3. * (3*a12 - rr)
                     endif
                  endif
                  cg%q(apot_i)%arr(i, j, k) = potential * pi * newtong * d0
               enddo
            enddo
         enddo

         cgl => cgl%nxt
      enddo

      call ppp_main%stop(cmp_label, PPP_PROB)

   end subroutine compute_maclaurin_potential

!> \brief Compute the L2 error norm of the error of computed potential with respect to the analytical solution

   subroutine finalize_problem_maclaurin

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: GEO_RPZ, pSUM, pMIN, pMAX, idlen, V_ESSENTIAL
      use dataio_pub,       only: msg, printinfo, warn
      use domain,           only: dom
      use grid_cont,        only: grid_container
      use mpisetup,         only: master, piernik_MPI_Allreduce
      use named_array_list, only: qna

      implicit none

      integer                        :: i, j, k, apot_i
      real, dimension(2)             :: norm, dev
      real                           :: potential, fac
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg
      character(len=idlen)           :: ffmt

      fac = 1.
      norm(:) = 0.
      dev(1) = huge(1.0)
      dev(2) = -dev(1)
      apot_i = qna%ind(apot_n)

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         if (.not. qna%exists(apot_n)) then
            if (master) call warn("[initproblem:finalize_problem_maclaurin] Cannot compare results with the analytical potential.")
            return
         endif

         do k = cg%ks, cg%ke
            do j = cg%js, cg%je
               do i = cg%is, cg%ie
                  if (cg%leafmap(i, j, k)) then
                     potential = cg%q(apot_i)%arr(i, j, k)
                     if (dom%geometry_type == GEO_RPZ) fac = cg%x(i)
                     norm(1) = norm(1) + (potential - cg%sgp(i, j, k))**2 * fac
                     norm(2) = norm(2) + potential**2 * fac
                     dev(1) = min(dev(1), (potential - cg%sgp(i, j, k))/potential)
                     dev(2) = max(dev(2), (potential - cg%sgp(i, j, k))/potential)
                  endif
               enddo
            enddo
         enddo
         cgl => cgl%nxt
      enddo

      call piernik_MPI_Allreduce(norm,   pSUM)
      call piernik_MPI_Allreduce(dev(1), pMIN)
      call piernik_MPI_Allreduce(dev(2), pMAX)

      if (master) then
         ffmt = "f12"
         if (any(abs(dev) > 1e6) .or. any(abs(dev) < 1e-4)) ffmt = "g14"
         write(msg,'(a,' // ffmt // '.6,a,2' // ffmt // '.6)')"[initproblem:finalize_problem_maclaurin] L2 error norm = ", sqrt(norm(1)/norm(2)), ", min and max error = ", dev(1:2)
         call printinfo(msg, V_ESSENTIAL)
      endif

   end subroutine finalize_problem_maclaurin

!> \brief Compute multipole potential field

   subroutine compute_mpole

      use multipole,        only: compute_mpole_potential
      use named_array_list, only: qna

      implicit none

      call compute_mpole_potential(qna%ind(mpole_n))

   end subroutine compute_mpole

!>
!! \brief This routine provides the "errp", "errm", "relerr" and "relerrm" values to be dumped to the .h5 file
!!
!! \details
!! * "errp"    is the difference between analytical potential and computed potential
!! * "relerr"  is the relative difference between analytical potential and multigrid solution
!! * "errm"    is the difference between analytical potential and multipole solution
!! * "relerrm" is the relative difference between analytical potential and multipole solution
!!
!! For "errm" and "relerr" use '$MULTIGRID_GRAVITY mpole_solver = "3D" /'
!! for realistic 3D potential evaluation in whole computational domain.
!! The default mpole_solver = "img_mass" will give only the "outer potential" correction.
!!
!! The values are calculated for nonperiodic, isolated case.
!! Any other configuration of boundary conditions will show large inaccuracies.
!<

   subroutine maclaurin_error_vars(var, tab, ierrh, cg)

      use dataio_pub,       only: die
      use grid_cont,        only: grid_container
      use named_array_list, only: qna

      implicit none

      character(len=*),               intent(in)    :: var
      real, dimension(:,:,:),         intent(inout) :: tab
      integer,                        intent(inout) :: ierrh
      type(grid_container), pointer,  intent(in)    :: cg

      if (.not. qna%exists(apot_n)) call die("[initproblem:maclaurin_error_vars] Cannot find apot_n")

      ierrh = 0
      select case (trim(var))
         case ("errp")
            tab(:,:,:) = cg%q(qna%ind(apot_n))%span(cg%ijkse) - cg%sgp(RNG)
         case ("relerr")
            where (abs(cg%q(qna%ind(apot_n))%span(cg%ijkse)) > tiny(1.))  ! An FP overflow may occur when sgp > apot / tiny()
               tab(:,:,:) = cg%sgp(RNG)/cg%q(qna%ind(apot_n))%span(cg%ijkse) -1.
            elsewhere
               tab(:,:,:) = 0.
            endwhere
         case ("errm")
            tab(:,:,:) = cg%q(qna%ind(apot_n))%span(cg%ijkse) - cg%q(qna%ind(mpole_n))%span(cg%ijkse)
         case ("relerrm")
            where (abs(cg%q(qna%ind(apot_n))%span(cg%ijkse)) > tiny(1.))  ! An FP overflow may occur when mpole > apot / tiny()
               tab(:,:,:) = cg%q(qna%ind(mpole_n))%span(cg%ijkse)/cg%q(qna%ind(apot_n))%span(cg%ijkse) -1.
            elsewhere
               tab(:,:,:) = 0.
            endwhere
         case default
            ierrh = -1
      end select

   end subroutine maclaurin_error_vars

!< \brief Analytical, point-like potential outside of semi-major axis

   real function ap_potential(x, y, z) result(phi)

      use constants, only: pi
      use units,     only: newtong

      implicit none

      real, intent(in) :: x, y, z

      real :: f

      if (a1 > 0.) then
         f = 4./3. * a1**3 * pi
      else
         f = 1.
      endif
      phi = - f * newtong * d0 / sqrt((x-x0)**2 + (y-y0)**2 + (z-z0)**2)

   end function ap_potential

!>
!! \brief Set up analytical potential at external boundaries
!!
!! \details This routine can be used to bypass multipole solver.
!! It can be used for diagnosing inaccuracies that come from laplacian or multigrid without
!! being bothered by limitations of the multipole representation and its limits.
!<

   subroutine maclaurin2bnd_potential

      use cg_leaves,  only: leaves
      use cg_list,    only: cg_list_element
      use constants,  only: xdim, ydim, zdim, LO, HI, GEO_XYZ
      use dataio_pub, only: die
      use domain,     only: dom
      use grid_cont,  only: grid_container
      use units,      only: fpiG

      implicit none

      integer :: i, j, k
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg

      if (dom%geometry_type /= GEO_XYZ) call die("[initproblem:maclaurin2bnd_potential] only cartesian geometry implemented")

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         if (any(cg%ext_bnd(xdim, :))) then
            do j = cg%js, cg%je
               do k = cg%ks, cg%ke
                  if (cg%ext_bnd(xdim, LO)) cg%mg%bnd_x(j, k, LO) = ap_potential(cg%fbnd(xdim, LO), cg%y(j), cg%z(k)) * fpiG
                  if (cg%ext_bnd(xdim, HI)) cg%mg%bnd_x(j, k, HI) = ap_potential(cg%fbnd(xdim, HI), cg%y(j), cg%z(k)) * fpiG
               enddo
            enddo
         endif

         if (any(cg%ext_bnd(ydim, :))) then
            do i = cg%is, cg%ie
               do k = cg%ks, cg%ke
                  if (cg%ext_bnd(ydim, LO)) cg%mg%bnd_y(i, k, LO) = ap_potential(cg%x(i), cg%fbnd(ydim, LO), cg%z(k)) * fpiG
                  if (cg%ext_bnd(ydim, HI)) cg%mg%bnd_y(i, k, HI) = ap_potential(cg%x(i), cg%fbnd(ydim, HI), cg%z(k)) * fpiG
               enddo
            enddo
         endif

         if (any(cg%ext_bnd(zdim, :))) then
            do i = cg%is, cg%ie
               do j = cg%js, cg%je
                  if (cg%ext_bnd(zdim, LO)) cg%mg%bnd_z(i, j, LO) = ap_potential(cg%x(i), cg%y(j), cg%fbnd(zdim, LO)) * fpiG
                  if (cg%ext_bnd(zdim, HI)) cg%mg%bnd_z(i, j, HI) = ap_potential(cg%x(i), cg%y(j), cg%fbnd(zdim, HI)) * fpiG
               enddo
            enddo
         endif
         cgl => cgl%nxt
      enddo

   end subroutine maclaurin2bnd_potential

end module initproblem
