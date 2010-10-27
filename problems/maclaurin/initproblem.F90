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
#include "piernik.def"
#include "macros.h"

module initproblem

   use problem_pub, only: problem_name, run_id

   real              :: x0, y0, z0, d0, a1, e, d1, p0, a3
   integer           :: nsub

   namelist /PROBLEM_CONTROL/ problem_name, run_id, x0, y0, z0, d0, a1, e, nsub

contains

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use constants,     only: pi
      use dataio_pub,    only: skip_advection, ierrh, msg, par_file, die, warn, namelist_errh, compare_namelist
      use mpisetup,      only: ierr, rbuff, cbuff, ibuff, proc, buffer_dim, comm, smalld, cbuff_len, &
           &                   MPI_CHARACTER, MPI_DOUBLE_PRECISION, MPI_INTEGER
      use types,         only: idlen

      implicit none

      integer, parameter :: maxsub = 10  !< upper limit for subsampling


      skip_advection = .true. ! skip sweeps in fluidupdate

      ! namelist default parameter values
      problem_name = 'Maclaurin sphere'  !< The default problem name
      run_id       = 'sph'               !< Auxiliary run identifier
      x0           = 0.0                 !< x-coordinate of the spheroid center
      y0           = 0.0                 !< y-coordinate of the spheroid center
      z0           = 0.0                 !< z-coordinate of the spheroid center
      d0           = 1.0                 !< Density inside the sphere
      a1           = 1.0                 !< Equatorial semimajor axis
      e            = 0.0                 !< Eccentricity; e>0 for flattened spheroids, e<0 for elongated spheroids
      nsub         = 3                   !< Subsampling factor

      if (proc == 0) then

         diff_nml(PROBLEM_CONTROL)

         cbuff(1) =  problem_name
         cbuff(2) =  run_id

         rbuff(1) = x0
         rbuff(2) = y0
         rbuff(3) = z0
         rbuff(4) = d0
         rbuff(5) = a1
         rbuff(6) = e

         ibuff(1) = nsub

      endif

      call MPI_Bcast(cbuff, cbuff_len*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_Bcast(ibuff,           buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_Bcast(rbuff,           buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

      if (proc /= 0) then

         problem_name = cbuff(1)
         run_id       = cbuff(2)(1:idlen)

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
         a3 = a1 * sqrt(1. - e**2) ! oblate, Maclaurin spheroid
      else
         a3 = a1 / sqrt(1. - e**2) ! prolate spheroid
      endif
      d1 = smalld                  ! ambient density
      p0 = 100.*tiny(d1)           ! pressure

      if (d0 < 0.) call die("[initproblem:read_problem_par] Negative average density.")

      if (nsub < 1) then
         if (proc == 0) call warn("[initproblem:read_problem_par] subsampling disabled.")
         nsub = 1
      else if (nsub > maxsub) then
         if (proc == 0)call warn("[initproblem:read_problem_par] too much subsampling.")
         nsub = maxsub
      endif

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine init_prob

      use arrays,        only: u, b
      use constants,     only: pi
      use dataio_pub,    only: msg, printinfo, warn
      use grid,          only: x, y, z, dx, dy, dz, nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax
      use initionized,   only: gamma_ion, idni, imxi, imzi, ieni
      use list_hdf5,     only: additional_attrs
      use mpisetup,      only: proc
      use types,         only: finalize_problem

      implicit none

      integer :: i, j, k, ii, jj, kk
      real    :: xx, yy, zz, dm

      do k = 1, nz
         do j = 1, ny
            do i = 1, nx

               dm = 0.
               do kk = -nsub+1, nsub-1, 2
                  zz = ((z(k) + kk*dz/(2.*nsub) - z0)/a3)**2
                  do jj = -nsub+1, nsub-1, 2
                     yy = ((y(j) + jj*dy/(2.*nsub) - y0)/a1)**2
                     do ii = -nsub+1, nsub-1, 2
                        xx = ((x(i) + ii*dx/(2.*nsub) - x0)/a1)**2

                        if (xx+yy+zz <= 1.) then
                           dm = dm + d0
                        else
                           dm = dm + d1
                        endif

                     enddo
                  enddo
               enddo
               u(idni, i, j, k) = dm / nsub**3

            enddo
         enddo
      enddo

      u(imxi:imzi, 1:nx, 1:ny, 1:nz) = 0.0

#ifndef ISO
      b(:, 1:nx, 1:ny, 1:nz) = 0.0
      u(ieni, 1:nx, 1:ny, 1:nz) = p0/(gamma_ion - 1.0)
#endif /* !ISO */

      if (proc == 0) then
         write(msg, '(3(a,g12.5),a)')"[initproblem:init_prob] Set up spheroid with a1 and a3 axes = ", a1, ", ", a3, " (eccentricity = ", e, ")"
         call printinfo(msg, .true.)
         if (x0-a1<xmin .or. x0+a1>xmax .or. y0-a1<ymin .or. y0+a1>ymax .or. z0-a3<zmin .or. z0+a3>zmax) &
              call warn("[initproblem:init_prob] Part of the spheroid is outside the domain")
         write(msg,'(2(a,g12.5))')   "[initproblem:init_prob] Density = ", d0, " mass = ", 4./3.*pi * a1**2 * a3 * d0
         call printinfo(msg, .true.)
      endif

      additional_attrs => init_prob_attrs
      finalize_problem => finalize_problem_maclaurin

   end subroutine init_prob

!-----------------------------------------------------------------------------

   subroutine init_prob_attrs(file_id)

      use constants, only: fpiG
      use hdf5,      only: HID_T, SIZE_T
      use h5lt,      only: h5ltset_attribute_double_f

      implicit none

      integer(HID_T),intent(in)  :: file_id

      integer(SIZE_T) :: bufsize = 1
      integer         :: error

      call h5ltset_attribute_double_f(file_id, "/", "rho0", [d0],   bufsize,error)
      call h5ltset_attribute_double_f(file_id, "/", "fpiG", [fpiG], bufsize,error)
      call h5ltset_attribute_double_f(file_id, "/", "a1",   [a1],   bufsize,error)
      call h5ltset_attribute_double_f(file_id, "/", "e",    [e],    bufsize,error)
      call h5ltset_attribute_double_f(file_id, "/", "x0",   [x0],   bufsize,error)
      call h5ltset_attribute_double_f(file_id, "/", "y0",   [y0],   bufsize,error)
      call h5ltset_attribute_double_f(file_id, "/", "z0",   [z0],   bufsize,error)

   end subroutine init_prob_attrs

!-----------------------------------------------------------------------------
!
! Oblate potential formula from Ricker_2008ApJS..176..293R
! Prolate potential formula from http://scienceworld.wolfram.com/physics/ProlateSpheroidGravitationalPotential.html
!
! ToDo: save the analytical potential in the restart file
!
   subroutine finalize_problem_maclaurin

      use arrays,        only: sgp
      use constants,     only: pi, newtong
      use dataio_pub,    only: msg, printinfo, warn
      use grid,          only: x, y, z, is, ie, js, je, ks, ke
      use mpisetup,      only: proc, comm3d, ierr, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_MIN, MPI_MAX, MPI_IN_PLACE

      implicit none

      integer            :: i, j, k
      real               :: potential, r2, rr
      real, dimension(2) :: norm, dev
      real               :: AA1 = 2./3., AA3 = 2./3., a12, a32, x02, y02, z02, lam, h
      real, parameter    :: small_e = 1e-3
      norm(:) = 0.
      dev(1) = huge(1.0)
      dev(2) = -dev(1)

      if (e < 0. .and. proc == 0) call warn("[initproblem:finalize_problem] e<0. not fully implemented yet!")

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
      do k = ks, ke
         z02 = (z(k)-z0)**2
         do j = js, je
            y02 = (y(j)-y0)**2
            do i = is, ie
               x02 = (x(i)-x0)**2
               r2 = (x02+y02)/a12 + z02/a32
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
                          (2.*(a32 - a12) + x02 + y02 - 2.*z02)/sqrt(a32 - a12) * log(h + sqrt(1. + h**2)) - &
                          (x02 + y02) * sqrt(a32 + lam)/(a12 + lam) + 2.*z02 / sqrt(a32 + lam) )
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
               potential = potential * pi * newtong * d0
               norm(1) = norm(1) + (potential - sgp(i, j, k))**2
               norm(2) = norm(2) + potential**2
               dev(1) = min(dev(1), (potential - sgp(i, j, k))/potential)
               dev(2) = max(dev(2), (potential - sgp(i, j, k))/potential)
            enddo
         enddo
      enddo

      call MPI_Allreduce(MPI_IN_PLACE, norm,   2, MPI_DOUBLE_PRECISION, MPI_SUM, comm3d, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, dev(1), 1, MPI_DOUBLE_PRECISION, MPI_MIN, comm3d, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, dev(2), 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm3d, ierr)

      if (proc == 0) then
         write(msg,'(a,f12.6,a,2f12.6)')"[initproblem:finalize_problem] L2 error norm = ", sqrt(norm(1)/norm(2)), ", min and max error = ", dev(1:2)
         call printinfo(msg)
      endif

   end subroutine finalize_problem_maclaurin

end module initproblem
