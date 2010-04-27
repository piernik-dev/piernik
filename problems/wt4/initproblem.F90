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

module initproblem
   use mpisetup, only: cbuff_len

   real              :: gamma_loc
   integer           :: nsub
   character(len=cbuff_len) :: problem_name
   character(len=3)         :: run_id
   character(len=cbuff_len) :: input_file

   namelist /PROBLEM_CONTROL/  problem_name, run_id, input_file, gamma_loc

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
      problem_name = 'wengen4'                     !< The default problem name
      input_file   = '/dev/shm/wt4/test4-512.alt'  !< File with initial conditions
      run_id       = 'tst'                         !< Auxiliary run identifier
      gamma_loc    = 1.4                           !< gamma used for calculating initial T distribution

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
         cbuff(2) =  input_file

         rbuff(1) = gamma_loc

         call MPI_BCAST(cbuff, cbuff_len*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
         call MPI_BCAST(ibuff,           buffer_dim, MPI_INTEGER,          0, comm, ierr)
         call MPI_BCAST(rbuff,           buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

      else

         call MPI_BCAST(cbuff, cbuff_len*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
         call MPI_BCAST(ibuff,           buffer_dim, MPI_INTEGER,          0, comm, ierr)
         call MPI_BCAST(rbuff,           buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

         problem_name = trim(cbuff(1))
         run_id       = cbuff(2)(1:3)
         input_file   = trim(cbuff(3))

         gamma_loc    = rbuff(1)

      endif

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine init_prob

      use mpisetup,     only: proc
      use arrays,       only: u, b, cs_iso2_arr
      use grid,         only: nx, ny, nz, nb
      use initionized,  only: idni, imxi, imyi, imzi
      use list_hdf5,    only: additional_attrs

      implicit none

      integer :: i, j, k, ii, jj, kk
      real    :: xx, yy, zz, dm, val

      open(1,file=input_file, status='old')
      do k = nb+1, nz-nb
         do j = nb+1, ny-nb
            do i = nb+1, nx-nb
               read(1,*) val
               u(idni, i, j, k) = val
            end do
         end do
      end do

      do k = nb+1, nz-nb
         do j = nb+1, ny-nb
            do i = nb+1, nx-nb
               read(1,*) val
               u(imxi, i, j, k) = val * u(idni, i, j, k)
            end do
         end do
      end do

      do k = nb+1, nz-nb
         do j = nb+1, ny-nb
            do i = nb+1, nx-nb
               read(1,*) val
               u(imyi, i, j, k) = val * u(idni, i, j, k)
            end do
         end do
      end do

      do k = nb+1, nz-nb
         do j = nb+1, ny-nb
            do i = nb+1, nx-nb
               read(1,*) val
               u(imzi, i, j, k) = val * u(idni, i, j, k)
            end do
         end do
      end do
      ! U = ( kB * T ) / (mean_mol_weight * (gamma - 1))
      ! cs2 = (gamma) * kB * T / mean_mol_weight.
      !   => cs2 = U * (gamma)(gamma - 1)
      do k = nb+1, nz-nb
         do j = nb+1, ny-nb
            do i = nb+1, nx-nb
               read(1,*) val
               cs_iso2_arr(i, j, k) = val * (gamma_loc)*(gamma_loc - 1.0)
            end do
         end do
      end do


      do i = 1,nb
         u(:,i,:,:)               = u(:,nb+1,:,:)
         u(:,nx-nb+i,:,:)         = u(:,nx-nb,:,:)
         cs_iso2_arr(i,:,:)       = cs_iso2_arr(nb+1,:,:)
         cs_iso2_arr(nx-nb+i,:,:) = cs_iso2_arr(nx-nb,:,:)

         u(:,:,i,:)               = u(:,:,nb+1,:)
         u(:,:,ny-nb+i,:)         = u(:,:,ny-nb,:)
         cs_iso2_arr(:,i,:)       = cs_iso2_arr(:,nb+1,:)
         cs_iso2_arr(:,ny-nb+i,:) = cs_iso2_arr(:,ny-nb,:)

         u(:,:,:,i)               = u(:,:,:,nb+1)
         u(:,:,:,nz-nb+i)         = u(:,:,:,nz-nb)
         cs_iso2_arr(:,:,i)       = cs_iso2_arr(:,:,nb+1)
         cs_iso2_arr(:,:,nz-nb+i) = cs_iso2_arr(:,:,nz-nb)
      enddo

      write(*,*) '[initproblem:init_problem]: minval(dens)    = ', minval(u(idni,:,:,:)),      ' maxval(dens) = ',    maxval(u(idni,:,:,:))
      write(*,*) '[initproblem:init_problem]: minval(cs_iso2) = ', minval(cs_iso2_arr(:,:,:)), ' maxval(cs_iso2) = ', maxval(cs_iso2_arr(:,:,:))

      b(:, 1:nx, 1:ny, 1:nz) = 0.0
#if __GNUC__ >=4 && __GNUC_MINOR__ >=4
      additional_attrs => init_prob_attrs
#endif

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

      call h5ltset_attribute_double_f(file_id, "/", "fpiG", [fpiG], bufsize,error)

   end subroutine init_prob_attrs

end module initproblem
