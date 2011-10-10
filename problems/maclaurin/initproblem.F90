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

module initproblem

   use constants, only: dsetnamelen

   implicit none

   private
   public :: read_problem_par, init_prob, problem_pointers

   ! namelist parameters
   real    :: x0, y0, z0, d0, a1, e
   integer :: nsub
   namelist /PROBLEM_CONTROL/ x0, y0, z0, d0, a1, e, nsub

   ! pivate data
   real :: d1, p0, a3
   character(len=dsetnamelen), parameter :: apot_n = "apot"

contains

!-----------------------------------------------------------------------------

   subroutine problem_pointers

      use dataio_user, only: user_vars_hdf5, additional_attrs
      use user_hooks,  only: finalize_problem

      implicit none

      additional_attrs => init_prob_attrs
      finalize_problem => finalize_problem_maclaurin
      user_vars_hdf5   => maclaurin_error_vars

   end subroutine problem_pointers

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use constants,   only: pi
      use dataio_pub,  only: ierrh, par_file, namelist_errh, compare_namelist, cmdl_nml, lun, getlun      ! QA_WARN required for diff_nml
      use dataio_pub,  only: die, warn
      use global,      only: smalld
      use mpisetup,    only: ierr, rbuff, ibuff, master, slave, buffer_dim, comm, FIRST
      use mpi,         only: MPI_DOUBLE_PRECISION, MPI_INTEGER

      implicit none

      integer, parameter :: maxsub = 10  !< upper limit for subsampling

      ! namelist default parameter values
      x0           = 0.0                 !< x-coordinate of the spheroid center
      y0           = 0.0                 !< y-coordinate of the spheroid center
      z0           = 0.0                 !< z-coordinate of the spheroid center
      d0           = 1.0                 !< Density inside the sphere
      a1           = 1.0                 !< Equatorial semimajor axis
      e            = 0.0                 !< Eccentricity; e>0 for flattened spheroids, e<0 for elongated spheroids
      nsub         = 3                   !< Subsampling factor

      if (master) then

         diff_nml(PROBLEM_CONTROL)

         rbuff(1) = x0
         rbuff(2) = y0
         rbuff(3) = z0
         rbuff(4) = d0
         rbuff(5) = a1
         rbuff(6) = e

         ibuff(1) = nsub

      endif

      call MPI_Bcast(ibuff, buffer_dim, MPI_INTEGER,          FIRST, comm, ierr)
      call MPI_Bcast(rbuff, buffer_dim, MPI_DOUBLE_PRECISION, FIRST, comm, ierr)

      if (slave) then

         x0           = rbuff(1)
         y0           = rbuff(2)
         z0           = rbuff(3)
         d0           = rbuff(4)
         a1           = rbuff(5)
         e            = rbuff(6)

         nsub         = ibuff(1)

      endif

      if (abs(e) >= 1.) call die("[initproblem:read_problem_par] |e|>=1.")
      if (e >= 0.) then            ! vertical axis
         a3 = a1 * sqrt(1. - e**2) ! oblate Maclaurin spheroid
      else
         a3 = a1 / sqrt(1. - e**2) ! prolate spheroid
      endif
      d1 = smalld                  ! ambient density
      p0 = 100.*tiny(d1)           ! pressure

      if (d0 < 0.) call die("[initproblem:read_problem_par] Negative average density.")

      if (nsub < 1) then
         if (master) call warn("[initproblem:read_problem_par] subsampling disabled.")
         nsub = 1
      else if (nsub > maxsub) then
         if (master)call warn("[initproblem:read_problem_par] too much subsampling.")
         nsub = maxsub
      endif

      call compute_maclaurin_potential

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine init_prob

      use constants,   only: pi, GEO_XYZ, GEO_RPZ, xdim, ydim, LO, HI
      use dataio_pub,  only: msg, printinfo, warn, die
      use domain,      only: dom, geometry_type
      use grid,        only: all_cg
      use gc_list,     only: cg_list_element
      use grid_cont,   only: grid_container
      use initionized, only: gamma_ion, idni, imxi, imzi, ieni
      use mpisetup,    only: master

      implicit none

      integer :: i, j, k, ii, jj, kk
      real    :: xx, yy, zz, rr, dm
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg

      cgl => all_cg%first
      do while (associated(cgl))
         cg => cgl%cg

         do k = cg%ks, cg%ke
            do j = cg%js, cg%je
               do i = cg%is, cg%ie

                  !< \todo use subsampling only near the surface of the spheroid
                  dm = 0.
                  do kk = -nsub+1, nsub-1, 2
                     zz = ((cg%z(k) + kk*cg%dz/(2.*nsub) - z0)/a3)**2
                     do jj = -nsub+1, nsub-1, 2
                        do ii = -nsub+1, nsub-1, 2

                           select case (geometry_type)
                              case (GEO_XYZ)
                                 yy = ((cg%y(j) + jj*cg%dy/(2.*nsub) - y0)/a1)**2
                                 xx = ((cg%x(i) + ii*cg%dx/(2.*nsub) - x0)/a1)**2
                                 rr = xx + yy + zz
                              case (GEO_RPZ)
                                 yy = cg%y(j) + jj*cg%dy/(2.*nsub) - y0
                                 xx = cg%x(i) + ii*cg%dx/(2.*nsub)
                                 rr = (xx**2 + x0**2 - 2. * xx * x0 * cos(yy))/a1**2 + zz
                              case default
                                 call die("[initproblem:init_prob] Unsupported geometry_type")
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
                  cg%u(idni, i, j, k) = dm / nsub**3

               enddo
            enddo
         enddo

         cg%u(imxi:imzi, cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) = 0.0

#ifndef ISO
#ifdef MAGNETIC
         cg%b(:, cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) = 0.0
#endif /* MAGNETIC */
         cg%u(ieni, cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) = p0/(gamma_ion - 1.0)
#endif /* !ISO */
         cgl => cgl%nxt
      enddo

      if (master) then
         write(msg, '(3(a,g12.5),a)')"[initproblem:init_prob] Set up spheroid with a1 and a3 axes = ", a1, ", ", a3, " (eccentricity = ", e, ")"
         call printinfo(msg, .true.)
         if (x0-a1<dom%edge(xdim, LO) .or. x0+a1>dom%edge(xdim, HI)) call warn("[initproblem:init_prob] Part of the spheroid is outside the domain in the X-direction.")
         if ( (geometry_type == GEO_XYZ .and. (y0-a1<dom%edge(ydim, LO) .or. y0+a1>dom%edge(ydim, HI))) .or. &
              (geometry_type == GEO_RPZ .and. (atan2(a1,x0) > minval([y0-dom%edge(ydim, LO), dom%edge(ydim, HI)-y0]))) ) & ! will fail when some one adds 2*k*pi to y0
              call warn("[initproblem:init_prob] Part of the spheroid is outside the domain")
         write(msg,'(2(a,g12.5))')   "[initproblem:init_prob] Density = ", d0, " mass = ", 4./3.*pi * a1**2 * a3 * d0
         call printinfo(msg, .true.)
      endif

   end subroutine init_prob

!-----------------------------------------------------------------------------
!
! This routine provides parameters useful for python/maclaurin.py
!

   subroutine init_prob_attrs(file_id)

      use units, only: fpiG
      use hdf5,  only: HID_T, SIZE_T
      use h5lt,  only: h5ltset_attribute_double_f

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

   end subroutine init_prob_attrs

!-----------------------------------------------------------------------------
!
! Oblate potential formula from Ricker_2008ApJS..176..293R
! Prolate potential formula from http://scienceworld.wolfram.com/physics/ProlateSpheroidGravitationalPotential.html
!
   subroutine compute_maclaurin_potential

      use constants,   only: pi, GEO_XYZ, GEO_RPZ, AT_IGNORE
      use dataio_pub,  only: warn, die
      use domain,      only: geometry_type
      use grid,        only: all_cg
      use gc_list,     only: cg_list_element
      use grid_cont,   only: grid_container
      use mpisetup,    only: master
      use units,       only: newtong

      implicit none

      integer            :: i, j, k
      real               :: potential, r2, rr
      real               :: AA1, AA3, a12, a32, x02, y02, z02, lam, h, cdphi
      real, parameter    :: small_e = 1e-3
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg
      real, dimension(:,:,:), pointer :: apot

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

      call all_cg%reg_var(apot_n, AT_IGNORE)

      cgl => all_cg%first
      do while (associated(cgl))
         cg => cgl%cg

         apot => cg%get_na_ptr(apot_n)

         do k = cg%ks, cg%ke
            z02 = (cg%z(k)-z0)**2
            do j = cg%js, cg%je
               do i = cg%is, cg%ie
                  select case (geometry_type)
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
                        call die("[initproblem:compute_maclaurin_potential] Invalid geometry_type.")
                        r2 = 0 ; x02 = 0 ; y02 = 0 ! suppress compiler warnings
                  end select
                  rr = r2 * a12
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
                        potential = - 4./3. * a1**3 / sqrt(rr)
                     endif
                  else
                     if (abs(e) > small_e) then
                        potential = - (AA1*(2*a12 - x02 - y02) + AA3 * (a32 - z02))
                     else
                        potential = - 2./3. * (3*a12 - rr)
                     endif
                  endif
                  apot(i, j, k) = potential * pi * newtong * d0
               enddo
            enddo
         enddo

         cgl => cgl%nxt
      enddo

   end subroutine compute_maclaurin_potential

!-----------------------------------------------------------------------------
!
! Here we compute the L2 error norm of the error of computed potential with respect to the analytical solution
!

   subroutine finalize_problem_maclaurin

      use constants,  only: GEO_RPZ, I_ONE, I_TWO
      use dataio_pub, only: msg, printinfo, warn
      use domain,     only: geometry_type
      use grid,       only: all_cg
      use gc_list,    only: cg_list_element
      use grid_cont,  only: grid_container
      use mpi,        only: MPI_DOUBLE_PRECISION, MPI_SUM, MPI_MIN, MPI_MAX, MPI_IN_PLACE
      use mpisetup,   only: master, comm, ierr

      implicit none

      integer            :: i, j, k
      real, dimension(2) :: norm, dev
      real               :: potential, fac
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg
      real, dimension(:,:,:), pointer :: apot

      fac = 1.
      norm(:) = 0.
      dev(1) = huge(1.0)
      dev(2) = -dev(1)

      cgl => all_cg%first
      do while (associated(cgl))
         cg => cgl%cg

         apot => cg%get_na_ptr(apot_n)
         if (.not. associated(apot))then
            if (master) call warn("[initproblem:finalize_problem_maclaurin] Cannot compare results with the analytical potential.")
            return
         endif

         do k = cg%ks, cg%ke
            do j = cg%js, cg%je
               do i = cg%is, cg%ie
                  potential =  apot(i, j, k)
                  if (geometry_type == GEO_RPZ) fac = cg%x(i)
                  norm(1) = norm(1) + (potential - cg%sgp(i, j, k))**2 * fac
                  norm(2) = norm(2) + potential**2 * fac
                  dev(1) = min(dev(1), (potential - cg%sgp(i, j, k))/potential)
                  dev(2) = max(dev(2), (potential - cg%sgp(i, j, k))/potential)
               enddo
            enddo
         enddo
         cgl => cgl%nxt
      enddo

      call MPI_Allreduce(MPI_IN_PLACE, norm,   I_TWO, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, dev(1), I_ONE, MPI_DOUBLE_PRECISION, MPI_MIN, comm, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, dev(2), I_ONE, MPI_DOUBLE_PRECISION, MPI_MAX, comm, ierr)

      if (master) then
         write(msg,'(a,f12.6,a,2f12.6)')"[initproblem:finalize_problem_maclaurin] L2 error norm = ", sqrt(norm(1)/norm(2)), ", min and max error = ", dev(1:2)
         call printinfo(msg)
      endif

   end subroutine finalize_problem_maclaurin

!-----------------------------------------------------------------------------
!
! This routine provides the "apot" and "errp" variablesvalues to be dumped to the .h5 file
! "apot" is the analytical potential solution for cell centers
! "errp" is the difference between analytical potential and computed potential
!

   subroutine maclaurin_error_vars(var, tab, ierrh, cg)

      use grid_cont, only: grid_container

      implicit none

      character(len=*), intent(in)                    :: var
      real(kind=4), dimension(:,:,:), intent(inout)   :: tab
      integer, intent(inout)                          :: ierrh
      type(grid_container), pointer, intent(in)       :: cg
      real, dimension(:,:,:), pointer :: apot

      ierrh = 0
      apot => cg%get_na_ptr(apot_n)
      select case (trim(var))
         case ("apot")
            tab(:,:,:) = real(apot(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke), 4)
         case ("errp")
            tab(:,:,:) = real(apot(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) - cg%sgp(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke), 4)
         case default
            ierrh = -1
      end select

   end subroutine maclaurin_error_vars

end module initproblem
