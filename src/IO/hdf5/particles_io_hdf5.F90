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

module particles_io_hdf5
! pulled by NBODY && HDF5

   implicit none

   private
   public  :: write_nbody_hdf5, read_nbody_hdf5

   contains

   subroutine write_nbody_hdf5

      use constants, only: fnamelen, idlen, ndims
      use global,    only: t
      use hdf5,      only: h5open_f, h5close_f, h5fclose_f, h5dcreate_f, h5dclose_f, h5dwrite_f
      use hdf5,      only: h5screate_simple_f, h5sclose_f, h5fcreate_f
      use hdf5,      only: HID_T, HSIZE_T, SIZE_T, H5F_ACC_TRUNC_F, H5T_NATIVE_DOUBLE
      use h5lt,      only: h5ltset_attribute_int_f, h5ltset_attribute_double_f
      use particle_types, only: pset

      implicit none

      integer                           :: i, n_part
      real, dimension(:,:), allocatable :: pos_table, vel_table
      character(len=fnamelen)           :: hdf_name
      integer(kind=4)                   :: error, rank
      integer(HID_T)                    :: file_id, dataspace_id, dataset_id
      integer(SIZE_T)                   :: bufsize
      integer(HSIZE_T), dimension(2)    :: dimm
      character(len=idlen)              :: pvars = 'pos', vvars = 'vel'

      n_part = size(pset%p, dim=1)
      allocate(pos_table(n_part,ndims))
      allocate(vel_table(n_part,ndims))
      do i = 1, n_part
         pos_table(i,:) = pset%p(i)%pos(:)
         vel_table(i,:) = pset%p(i)%vel(:)
      enddo

      rank    = 2
      dimm    = [n_part,ndims]
      bufsize = 1

      write(hdf_name,'(a)') 'test_01.h5'

      call h5open_f(error)
      call h5fcreate_f(hdf_name, H5F_ACC_TRUNC_F, file_id, error)

      call h5ltset_attribute_int_f   (file_id, "/", "npart", [n_part], bufsize, error)
      call h5ltset_attribute_double_f(file_id, "/", "time",  [t],      bufsize, error)

      call h5screate_simple_f(rank, dimm, dataspace_id, error)
      call h5dcreate_f(file_id, pvars, H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, error)
      call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, pos_table, dimm, error)
      call h5dclose_f(dataset_id, error)
      call h5sclose_f(dataspace_id, error)

      call h5screate_simple_f(rank, dimm, dataspace_id, error)
      call h5dcreate_f(file_id, vvars, H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, error)
      call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, vel_table, dimm, error)
      call h5dclose_f(dataset_id, error)
      call h5sclose_f(dataspace_id, error)

      call h5fclose_f(file_id, error)
      call h5close_f(error)

      deallocate(pos_table, vel_table)

   end subroutine write_nbody_hdf5

   subroutine read_nbody_hdf5(table, n)

      use constants, only: fnamelen
      use hdf5,      only: h5open_f, h5close_f, h5fopen_f, h5fclose_f, h5dopen_f, h5dclose_f, h5dread_f
      use hdf5,      only: HID_T, HSIZE_T, H5F_ACC_RDONLY_F, H5T_NATIVE_DOUBLE
      use h5lt,      only: h5ltget_attribute_int_f

      implicit none

      integer,              intent(in)  :: n
      real, dimension(n,3), intent(out) :: table
      character(len=fnamelen)           :: hdf_name = 'test_01.h5'
      integer(kind=4)                   :: error, time
      integer(HID_T)                    :: file_id, dataset_id
      integer(HSIZE_T), dimension(2)    :: dimm
      character(len=9)                  :: vars='positions'
      integer,dimension(1)              :: rbuf

      dimm(1) = n
      dimm(2) = 3

      call h5open_f(error)
      call h5fopen_f(hdf_name, H5F_ACC_RDONLY_F, file_id, error)

      call h5ltget_attribute_int_f  (file_id,  "/", "time", rbuf, error)
      time = rbuf(1)

      call h5dopen_f(file_id, vars, dataset_id, error)
      call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, table, dimm, error)
      call h5dclose_f(dataset_id, error)

      call h5fclose_f(file_id, error)
      call h5close_f(error)

   end subroutine read_nbody_hdf5

end module particles_io_hdf5
