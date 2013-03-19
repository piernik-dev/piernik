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

!> \brief Implementation of selfgravitating spheroid, for which an analytical solution of the gravitational potential is known.

module initproblem

   use constants, only: dsetnamelen

   implicit none

   private
   public :: read_problem_par, problem_initial_conditions, problem_pointers

   ! namelist parameters
   real :: x0     !< X-position of the spheroid
   real :: y0     !< Y-position of the spheroid
   real :: z0     !< Z-position of the spherodi
   real :: d0     !< density of the spheroid
   real :: a1     !< equatorial radius of the spheroid
   real :: e      !< polar eccentricity of the spheroid; e>0 gives oblate object, e<0 gives prolate object
   real :: ref_thr !< refinement threshold
   real :: deref_thr !< derefinement threshold
   integer(kind=4) :: nsub !< subsampling on the grid

   namelist /PROBLEM_CONTROL/ x0, y0, z0, d0, a1, e, ref_thr, deref_thr, nsub

   ! private data
   real :: d1 !< ambient density
   real :: p0 !< pressure
   real :: a3 !< length of polar radius of the spheroid
   character(len=dsetnamelen), parameter :: apot_n = "apot" !< name of the analytical potential field
   character(len=dsetnamelen), parameter :: asrc_n = "asrc" !< name of the source fiels used for "ares" calculation (auxiliary space)
   character(len=dsetnamelen), parameter :: ares_n = "ares" !< name of the numerical residuum with respect to analytical potential field
#ifdef MACLAURIN_PROBLEM
   character(len=dsetnamelen), parameter :: apt_n  = "apt"  !< name of the potential as it was due to point-like source
#endif /* MACLAURIN_PROBLEM */

contains

!> \brief Set up custom pointers to tweak the code execution according to our needs

   subroutine problem_pointers

      use dataio_user, only: user_vars_hdf5, user_attrs_wr
      use user_hooks,  only: finalize_problem, problem_refine_derefine

      implicit none

      user_attrs_wr    => problem_initial_conditions_attrs
      finalize_problem => finalize_problem_maclaurin
      user_vars_hdf5   => maclaurin_error_vars
      problem_refine_derefine => mark_surface

   end subroutine problem_pointers

!> \brief Read the runtime parameters specified in the namelist and set up some auxiliary values

   subroutine read_problem_par

      use cg_list_global, only: all_cg
      use constants,      only: pi, GEO_XYZ, GEO_RPZ, xdim, ydim, LO, HI
      use dataio_pub,     only: nh      ! QA_WARN required for diff_nml
      use dataio_pub,     only: die, warn, msg, printinfo
      use domain,         only: dom
      use global,         only: smalld
      use mpisetup,       only: rbuff, ibuff, master, slave, piernik_MPI_Bcast
      use multigridvars,  only: ord_prolong
      use particle_pub,   only: pset

      implicit none

      integer, parameter :: maxsub = 10  !< upper limit for subsampling
      d1 = smalld                  ! ambient density

      ! namelist default parameter values
      x0           = 0.0                 !< x-coordinate of the spheroid center
      y0           = 0.0                 !< y-coordinate of the spheroid center
      z0           = 0.0                 !< z-coordinate of the spheroid center
      d0           = 1.0                 !< Density inside the sphere
      a1           = 1.0                 !< Equatorial semimajor axis
      e            = 0.0                 !< Eccentricity; e>0 for flattened spheroids, e<0 for elongated spheroids
      nsub         = 3                   !< Subsampling factor

      ref_thr      = max(1e-3, 2.*smalld/d0)    !< Refine if density difference is greater than this value
      deref_thr    = max(ref_thr**2, smalld/d0) !< Derefine if density difference is smaller than this value

      if (master) then

         diff_nml(PROBLEM_CONTROL)

         rbuff(1) = x0
         rbuff(2) = y0
         rbuff(3) = z0
         rbuff(4) = d0
         rbuff(5) = a1
         rbuff(6) = e
         rbuff(7) = ref_thr
         rbuff(8) = deref_thr

         ibuff(1) = nsub

      endif

      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(rbuff)

      if (slave) then

         x0           = rbuff(1)
         y0           = rbuff(2)
         z0           = rbuff(3)
         d0           = rbuff(4)
         a1           = rbuff(5)
         e            = rbuff(6)
         ref_thr      = rbuff(7)
         deref_thr    = rbuff(8)

         nsub         = ibuff(1)

      endif

      if (a1 <= 0.) then ! point-like source
         a1 = 0.
         e = 0.
      endif

      if (abs(e) >= 1.) call die("[initproblem:read_problem_par] |e|>=1.")
      if (e >= 0.) then            ! vertical axis
         a3 = a1 * sqrt(1. - e**2) ! oblate Maclaurin spheroid
      else
         a3 = a1 / sqrt(1. - e**2) ! prolate spheroid
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

      if (ref_thr <= deref_thr) call die("[initproblem:read_problem_par] ref_thr <= deref_thr")

      if (a1 == 0.) call pset%add(d0, [ x0, y0, z0 ], [0.0, 0.0, 0.0])

      if (master) then
         if (a1 > 0.) then
            write(msg, '(3(a,g12.5),a)')"[initproblem:problem_initial_conditions] Set up spheroid with a1 and a3 axes = ", a1, ", ", a3, " (eccentricity = ", e, ")"
         else
            write(msg, '(a)')"[initproblem:problem_initial_conditions] Set up point-like mass"
         endif
         call printinfo(msg, .true.)
         if (x0-a1<dom%edge(xdim, LO) .or. x0+a1>dom%edge(xdim, HI)) call warn("[initproblem:problem_initial_conditions] Part of the spheroid is outside the domain in the X-direction.")
         if ( (dom%geometry_type == GEO_XYZ .and. (y0-a1<dom%edge(ydim, LO) .or. y0+a1>dom%edge(ydim, HI))) .or. &
              (dom%geometry_type == GEO_RPZ .and. (atan2(a1,x0) > minval([y0-dom%edge(ydim, LO), dom%edge(ydim, HI)-y0]))) ) & ! will fail when some one adds 2*k*pi to y0
              call warn("[initproblem:problem_initial_conditions] Part of the spheroid is outside the domain")
         if (a1 > 0.) then
            write(msg,'(2(a,g12.5))')   "[initproblem:problem_initial_conditions] Density = ", d0, " mass = ", 4./3.*pi * a1**2 * a3 * d0
         else
            write(msg,'(a,g12.5)')   "[initproblem:problem_initial_conditions] Mass = ", d0
         endif
         call printinfo(msg, .true.)
      endif

      call all_cg%reg_var(apot_n, ord_prolong = ord_prolong)
      call all_cg%reg_var(ares_n)
      call all_cg%reg_var(asrc_n)
#ifdef MACLAURIN_PROBLEM
      call all_cg%reg_var(apt_n)
#endif /* MACLAURIN_PROBLEM */

   end subroutine read_problem_par

!> \brief Set up the initial conditions. Note that this routine can be called multiple times during initial iterations of refinement structure

   subroutine problem_initial_conditions

      use cg_leaves,         only: leaves
      use cg_list,           only: cg_list_element
      use constants,         only: GEO_XYZ, GEO_RPZ, xdim, ydim, zdim, LO, HI
      use dataio_pub,        only: die, msg, printinfo
      use domain,            only: dom
      use fluidindex,        only: iarr_all_dn, iarr_all_mx, iarr_all_my, iarr_all_mz
      use global,            only: dirty_debug, no_dirty_checks
      use grid_cont,         only: grid_container
      use named_array_list,  only: qna
      use mpisetup,          only: master
      use multigrid_gravity, only: residual
      use units,             only: fpiG

      implicit none

      integer                        :: i, j, k, ii, jj, kk
      real                           :: xx, yy, zz, rr, dm
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg

      call compute_maclaurin_potential

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

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

                              if (rr <= 1.) then
                                 dm = dm + d0
                              else
                                 dm = dm + d1
                              endif

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
            cg%u(iarr_all_dn, cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) = 0.0
         endif

         cg%u(iarr_all_mx, cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) = 0.0
         cg%u(iarr_all_my, cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) = 0.0
         cg%u(iarr_all_mz, cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) = 0.0

#ifdef MAGNETIC
         cg%b(:, cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) = 0.0
#endif /* MAGNETIC */

         cgl => cgl%nxt
      enddo

      ! Compute residuum for analytical solution and print its norm
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         if (dirty_debug .and. .not. no_dirty_checks) then
            cg%wa(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) = cg%q(qna%ind(apot_n))%arr(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke)
            cg%q(qna%ind(apot_n))%arr = huge(1.)
            cg%q(qna%ind(apot_n))%arr(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) = cg%wa(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke)
            cg%q(qna%ind(asrc_n))%arr = huge(1.)
         endif
         cg%q(qna%ind(asrc_n))%arr(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) = fpiG * sum(cg%u(iarr_all_dn, cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke), dim=1)
         cgl => cgl%nxt
      enddo

      call residual(leaves, qna%ind(asrc_n), qna%ind(apot_n), qna%ind(ares_n))
      call leaves%check_dirty(qna%ind(ares_n), "a-residual")

      ! clear residual next to the external boundary as it is affected by the way the potential is extrapolated into the guardcells.
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

      write(msg,'(a,f13.10)')"[initproblem:problem_initial_conditions] Analytical norm residual/source= ",leaves%norm_sq(qna%ind(ares_n))/leaves%norm_sq(qna%ind(asrc_n))
      if (master) call printinfo(msg)

   end subroutine problem_initial_conditions

!> \brief Provides parameters useful for python/maclaurin.py in .h5 files

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

!>
!! \brief Calculate analytical potential for given spheroid.
!!
!! \details We assume that the spheroid is fully contained in the computational domain here.
!! Oblate potential formula from Ricker_2008ApJS..176..293R
!! Prolate potential formula from http://scienceworld.wolfram.com/physics/ProlateSpheroidGravitationalPotential.html
!<
   subroutine compute_maclaurin_potential

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: pi, GEO_XYZ, GEO_RPZ
      use dataio_pub,       only: warn, die
      use domain,           only: dom
      use grid_cont,        only: grid_container
      use mpisetup,         only: master
      use named_array_list, only: qna
#ifdef MACLAURIN_PROBLEM
      use problem_pub,      only: xs, as, ap_potential
#endif /* MACLAURIN_PROBLEM */
      use units,            only: newtong

      implicit none

      integer                        :: i, j, k, apot_i
      real                           :: potential, r2, rr
      real                           :: AA1, AA3, a12, a32, x02, y02, z02, lam, h, cdphi
      real, parameter                :: small_e = 1e-3
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg

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
                  if (a12 /= 0.) then
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
                        if (a1 == 0.) then
                           if (rr /= 0) then
                              potential = - 1. / (pi * sqrt(rr)) ! A bit dirty: we interpret d0 as a mass for point-like sources
                           else
                              potential = - sqrt(sum(cg%idl2))*0.5
                              ! The factor of 0.5 was guessed. It gives small departures of potential in the cell containing the particle for 4th order Laplacian.
                              ! For 2nd order Laplacian (default) 0.58 looks a bit better because it compensates 4th order errors of the operator.
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
#ifdef MACLAURIN_PROBLEM

      xs = [ x0, y0, z0 ]
      as = - 4./3. * a1**3 * pi * newtong * d0 !> \todo add correction for e /= 0

      apot_i = qna%ind(apt_n)
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         do k = cg%ks, cg%ke
            do j = cg%js, cg%je
               do i = cg%is, cg%ie
                  cg%q(apot_i)%arr(i, j, k) = ap_potential(cg%x(i), cg%y(j), cg%z(k))
               enddo
            enddo
         enddo

         cgl => cgl%nxt
      enddo
#endif /* MACLAURIN_PROBLEM */

   end subroutine compute_maclaurin_potential

!> \brief Compute the L2 error norm of the error of computed potential with respect to the analytical solution

   subroutine finalize_problem_maclaurin

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: GEO_RPZ, pSUM, pMIN, pMAX
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
         write(msg,'(a,f12.6,a,2f12.6)')"[initproblem:finalize_problem_maclaurin] L2 error norm = ", sqrt(norm(1)/norm(2)), ", min and max error = ", dev(1:2)
         call printinfo(msg)
      endif

   end subroutine finalize_problem_maclaurin

!>
!! \brief This routine provides the "apot" and "errp" variablesvalues to be dumped to the .h5 file
!!
!! \details
!! * "apot" is the analytical potential solution for cell centers
!! * "errp" is the difference between analytical potential and computed potential
!<

   subroutine maclaurin_error_vars(var, tab, ierrh, cg)

      use dataio_pub,       only: die
      use grid_cont,        only: grid_container
      use named_array_list, only: qna

      implicit none

      character(len=*),               intent(in)    :: var
      real(kind=4), dimension(:,:,:), intent(inout) :: tab
      integer,                        intent(inout) :: ierrh
      type(grid_container), pointer,  intent(in)    :: cg

      if (.not. qna%exists(apot_n)) call die("[initproblem:maclaurin_error_vars] Cannot find apot_n")

      ierrh = 0
      select case (trim(var))
         case ("errp")
            tab(:,:,:) = real(cg%q(qna%ind(apot_n))%span(cg%ijkse) - cg%sgp(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke), 4)
         case ("relerr")
            where (cg%q(qna%ind(apot_n))%span(cg%ijkse) /= 0.)
               tab(:,:,:) = real(cg%sgp(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke)/cg%q(qna%ind(apot_n))%span(cg%ijkse) -1., 4)
            elsewhere
               tab(:,:,:) = 0.
            endwhere
#ifdef MACLAURIN_PROBLEM
         case ("a-pt")
            tab(:,:,:) = real(cg%q(qna%ind(apot_n))%span(cg%ijkse) - cg%q(qna%ind(apt_n))%span(cg%ijkse), 4)
#endif /* MACLAURIN_PROBLEM */
         case default
            ierrh = -1
      end select

   end subroutine maclaurin_error_vars

!> \brief Request refinement along the surface of the ellipsoid. Derefine inside and outside the ellipsoid if possible.

   subroutine mark_surface

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use fluidindex,       only: iarr_all_dn
!      use named_array_list, only: wna
      use refinement,       only: ref_flag

      implicit none

      type(cg_list_element), pointer :: cgl
      real :: delta_dens, dmin, dmax
      integer :: id

!      call leaves%internal_boundaries_4d(wna%fi) !< enable it as soon as c2f and f2c routines will work

      cgl => leaves%first
      do while (associated(cgl))
         if (any(cgl%cg%leafmap)) then
            dmax = -huge(1.)
            dmin =  huge(1.)
            do id = lbound(iarr_all_dn, dim=1), ubound(iarr_all_dn, dim=1)
               dmax = max(dmax, maxval(cgl%cg%u(id, cgl%cg%is:cgl%cg%ie, cgl%cg%js:cgl%cg%je, cgl%cg%ks:cgl%cg%ke), mask=cgl%cg%leafmap))
               dmin = min(dmin, minval(cgl%cg%u(id, cgl%cg%is:cgl%cg%ie, cgl%cg%js:cgl%cg%je, cgl%cg%ks:cgl%cg%ke), mask=cgl%cg%leafmap))
            enddo
            delta_dens = dmax - dmin
            !> \warning only selfgravitating fluids should be checked
            cgl%cg%refine_flags = ref_flag( delta_dens >= ref_thr*d0, delta_dens < deref_thr*d0 )
         endif
         cgl => cgl%nxt
      enddo

   end subroutine mark_surface

end module initproblem
