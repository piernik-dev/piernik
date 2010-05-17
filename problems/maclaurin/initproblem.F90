! $Id: initproblem.F90 1448 2009-12-15 03:45:01Z gawrysz $
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

module initproblem

   use problem_pub, only: problem_name, run_id

   real              :: x0, y0, z0, d0, a1, e, d1, p0, a3
   integer           :: nsub

   namelist /PROBLEM_CONTROL/  problem_name, run_id, x0, y0, z0, d0, a1, e, nsub

contains

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use grid, only : xmin, xmax, ymin, ymax, zmin, zmax
      use errh, only : namelist_errh, die
      use mpisetup, only : cwd, ierr, rbuff, cbuff, ibuff, proc, &
         MPI_CHARACTER, MPI_DOUBLE_PRECISION, MPI_INTEGER, &
         buffer_dim, comm, smalld
      use constants, only: pi

      implicit none

      integer, parameter :: maxsub = 10  !< upper limit for subsampling

      integer :: ierrh
      character(LEN=100) :: par_file, tmp_log_file

      ! namelist default parameter values
      problem_name = 'Maclaurin sphere'  !< The default problem name
      run_id       = 'sph'               !< Auxiliary run identifier
      x0           = 0.0                 !< x-coordinate of the spheroid center
      y0           = 0.0                 !< y-coordinate of the spheroid center
      z0           = 0.0                 !< z-coordinate of the spheroid center
      d0           = 1.0                 !< Density inside the sphere
      a1           = 1.0                 !< Equatorial semimajor axis
      e            = 0.0                 !< Eccentricity
      nsub         = 3                   !< Subsampling factor

      if(proc == 0) then
         par_file = trim(cwd)//'/problem.par'
         tmp_log_file = trim(cwd)//'/tmp.log'

         open(1,file=par_file)
            read(unit=1,nml=PROBLEM_CONTROL,iostat=ierrh)
            call namelist_errh(ierrh,'PROBLEM_CONTROL')
         close(1)
         open(3, file=tmp_log_file, position='append')
            write(3,nml=PROBLEM_CONTROL)
            write(3,*)
         close(3)
      endif

      if(proc == 0) then

         cbuff(1) =  problem_name
         cbuff(2) =  run_id

         rbuff(1) = x0
         rbuff(2) = y0
         rbuff(3) = z0
         rbuff(4) = d0
         rbuff(5) = a1
         rbuff(6) = e

         ibuff(1) = nsub

         call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
         call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
         call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

      else

         call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
         call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
         call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

         problem_name = cbuff(1)
         run_id       = cbuff(2)(1:3)

         x0           = rbuff(1)
         y0           = rbuff(2)
         z0           = rbuff(3)
         d0           = rbuff(4)
         a1           = rbuff(5)
         e            = rbuff(6)

         nsub         = ibuff(1) ! unused yet

      endif

      a3 = a1 * sqrt(1. - e**2) ! vertical axis
      d1 = smalld               ! ambient density
      p0 = 100.*tiny(d1)        ! pressure

      if (d0 < 0.) call die("[initproblem:read_problem_par] Negative average density.")
      if (x0-a1<xmin .or. x0+a1>xmax .or. y0-a1<ymin .or. y0+a1>ymax .or. z0-a3<zmin .or. z0+a3>zmax) &
           write(*,'(a)')"[initproblem:read_problem_par] Warning: part of the spheroid is outside the domain"

      if (nsub < 1) then
         if (proc == 0) write(*,'(a)')"[initproblem:read_problem_par] subsampling disabled"
         nsub = 1
      else if (nsub > maxsub) then
         if (proc == 0) write(*,'(a)')"[initproblem:read_problem_par] too much subsampling"
         nsub = maxsub
      end if

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine init_prob

      use mpisetup,     only: proc
      use arrays,       only: u, b
      use constants,    only: pi
      use grid,         only: x, y, z, dx, dy, dz, nx, ny, nz, xmin, ymin, zmin
      use initionized,  only: gamma_ion, idni, imxi, imzi, ieni
      use list_hdf5,    only: additional_attrs
      use types,        only: finalize_problem

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
                        end if

                     end do
                  end do
               end do
               u(idni, i, j, k) = dm / nsub**3

            end do
         end do
      end do

      u(imxi:imzi, 1:nx, 1:ny, 1:nz) = 0.0

#ifndef ISO
      b(:, 1:nx, 1:ny, 1:nz) = 0.0
      u(ieni, 1:nx, 1:ny, 1:nz) = p0/(gamma_ion - 1.0)
#endif /* ISO */

      if (proc == 0) then
         write(*,'(3(a,g12.5),a)')"[initproblem:init_prob] Set up spheroid with a1 and a3 axes = ", a1, ", ", a3, " (eccentricity = ", e, ")"
         write(*,'(2(a,g12.5))')  "[initproblem:init_prob] Density = ", d0, " mass = ", 4./3.*pi*a1**2*a3*d0
      endif

      additional_attrs => init_prob_attrs

      finalize_problem => finalize_problem_maclaurin
      return
   end subroutine init_prob

!-----------------------------------------------------------------------------

   subroutine init_prob_attrs(file_id)
      use hdf5, only : HID_T, SIZE_T
      use h5lt, only : h5ltset_attribute_double_f
      use constants, only : fpiG
      implicit none
      integer(HID_T),intent(in)  :: file_id
      integer(SIZE_T) :: bufsize = 1
      integer :: error

      call h5ltset_attribute_double_f(file_id, "/", "rho0", [d0],   bufsize,error)
      call h5ltset_attribute_double_f(file_id, "/", "fpiG", [fpiG], bufsize,error)
      call h5ltset_attribute_double_f(file_id, "/", "a1",   [a1],   bufsize,error)
      call h5ltset_attribute_double_f(file_id, "/", "x0",   [x0],   bufsize,error)
      call h5ltset_attribute_double_f(file_id, "/", "y0",   [y0],   bufsize,error)
      call h5ltset_attribute_double_f(file_id, "/", "z0",   [z0],   bufsize,error)

   end subroutine init_prob_attrs

!-----------------------------------------------------------------------------

   subroutine finalize_problem_maclaurin

      use mpisetup,    only: proc, comm3d, ierr, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_MIN, MPI_MAX, MPI_IN_PLACE
      use constants,   only: pi, newtong
      use grid,        only: x, y, z, is, ie, js, je, ks, ke
      use arrays,      only: mgp

      implicit none

      integer :: i, j, k
      real    :: potential, r2, rr, norm(2), dev(2)

      norm(:) = 0.
      dev(1) = huge(1.0)
      dev(2) = -dev(1)

      if (e> tiny(0.)) then
         if (proc == 0) write (*,'(a)')"[initproblem:finalize_problem] e>0 not implemented yet"
         return
      end if
      do k = ks, ke
         do j = js, je
            do i = is, ie
               !if (e == 0)
               r2 = ((x(i)-x0)/a1)**2 + ((y(j)-y0)/a1)**2 + ((z(k)-z0)/a3)**2
               rr = r2 * a1**2
               if (r2 > 1.) then
                  potential = - 4./3. * pi * newtong * d0 * a1**3 / sqrt(rr)
               else
                  potential = - 2./3. * pi * newtong * d0 * (3*a1**2 - rr)
               end if
               norm(1) = norm(1) + (potential - mgp(i, j, k))**2
               norm(2) = norm(2) + potential**2
               dev(1) = min(dev(1), (potential - mgp(i, j, k))/potential)
               dev(2) = max(dev(2), (potential - mgp(i, j, k))/potential)
            end do
         end do
      end do

      call MPI_Allreduce(MPI_IN_PLACE, norm,   2, MPI_DOUBLE_PRECISION, MPI_SUM, comm3d, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, dev(1), 1, MPI_DOUBLE_PRECISION, MPI_MIN, comm3d, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, dev(2), 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm3d, ierr)

      if (proc == 0) write(*,'(a,f12.5,a,2f12.5)')"[initproblem:finalize_problem] L2 error norm = ", sqrt(norm(1)/norm(2)), " min and max error = ", dev(1:2)

   end subroutine finalize_problem_maclaurin

end module initproblem
