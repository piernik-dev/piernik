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

   subroutine write_nbody_hdf5(fname)

      use constants,      only: cwdlen, idlen, ndims, tmr_hdf
      use dataio_pub,     only: msg, printinfo, printio, thdf
      use global,         only: t
      use hdf5,           only: h5open_f, h5close_f, h5fclose_f, h5dcreate_f, h5dclose_f, h5dwrite_f, h5screate_simple_f, h5sclose_f, h5fcreate_f
      use hdf5,           only: HID_T, HSIZE_T, SIZE_T, H5F_ACC_TRUNC_F, H5T_NATIVE_DOUBLE
      use h5lt,           only: h5ltset_attribute_int_f, h5ltset_attribute_double_f
      use mpisetup,       only: master
      use particle_types, only: pset
      use timer,          only: set_timer

      implicit none

      character(len=*),      intent(in) :: fname
      character(len=cwdlen)             :: filename
      character(len=idlen)              :: mvars = 'mas', pvars = 'pos', vvars = 'vel'
      integer                           :: flen, i, n_part
      integer(kind=4)                   :: error, rank1 = 1, rank2 = 2
      integer(HID_T)                    :: file_id, dataspace_id, dataset_id
      integer(SIZE_T)                   :: bufsize
      integer(HSIZE_T), dimension(1)    :: dimm
      integer(HSIZE_T), dimension(2)    :: dimv
      real, dimension(:),   allocatable :: mas_h
      real, dimension(:,:), allocatable :: pos_h, vel_h

      thdf = set_timer(tmr_hdf,.true.)
      flen = len(trim(fname))
      write(filename,'(3a)') fname(1:flen-7), 'p', fname(flen-6:flen)
      if (master) then
         write(msg,'(a,1x,2a)') 'Writing datafile ', trim(filename), " ... "
         call printio(msg, .true.)
      endif

      n_part = size(pset%p, dim=1)
      allocate(pos_h(n_part,ndims), vel_h(n_part,ndims), mas_h(n_part))
      do i = 1, n_part
         mas_h(i)   = pset%p(i)%mass
         pos_h(i,:) = pset%p(i)%pos(:)
         vel_h(i,:) = pset%p(i)%vel(:)
      enddo

      dimm    = [n_part]
      dimv    = [n_part,ndims]
      bufsize = 1

      call h5open_f(error)
      call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

      call h5ltset_attribute_int_f   (file_id, "/", "npart", [n_part], bufsize, error)
      call h5ltset_attribute_double_f(file_id, "/", "time",  [t],      bufsize, error)

      call h5screate_simple_f(rank1, dimm, dataspace_id, error)
      call h5dcreate_f(file_id, mvars, H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, error)
      call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, mas_h, dimm, error)
      call h5dclose_f(dataset_id, error)
      call h5sclose_f(dataspace_id, error)

      call h5screate_simple_f(rank2, dimv, dataspace_id, error)
      call h5dcreate_f(file_id, pvars, H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, error)
      call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, pos_h, dimv, error)
      call h5dclose_f(dataset_id, error)
      call h5sclose_f(dataspace_id, error)

      call h5screate_simple_f(rank2, dimv, dataspace_id, error)
      call h5dcreate_f(file_id, vvars, H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, error)
      call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, vel_h, dimv, error)
      call h5dclose_f(dataset_id, error)
      call h5sclose_f(dataspace_id, error)

      call h5fclose_f(file_id, error)
      call h5close_f(error)

      deallocate(pos_h, vel_h, mas_h)

      thdf = set_timer(tmr_hdf)
      if (master) then
         write(msg,'(a6,f10.2,a2)') ' done ', thdf, ' s'
         call printinfo(msg, .true.)
      endif

   end subroutine write_nbody_hdf5

   subroutine read_nbody_hdf5(fname, table, n)

      use constants, only: idlen, ndims
      use hdf5,      only: h5open_f, h5close_f, h5fopen_f, h5fclose_f, h5dopen_f, h5dclose_f, h5dread_f
      use hdf5,      only: HID_T, HSIZE_T, H5F_ACC_RDONLY_F, H5T_NATIVE_DOUBLE
      use h5lt,      only: h5ltget_attribute_double_f, h5ltget_attribute_int_f

      implicit none

      character(len=*),           intent(in)  :: fname
      integer,                    intent(in)  :: n
      real, dimension(2,n,ndims), intent(out) :: table
      character(len=idlen)                    :: pvars = 'pos', vvars = 'vel'
      integer                                 :: n_part
      integer(kind=4)                         :: error
      integer(HID_T)                          :: file_id, dataset_id
      integer(HSIZE_T), dimension(2)          :: dimm
      integer, dimension(1)                   :: ibuf
      real                                    :: time
      real, dimension(1)                      :: rbuf

      call h5open_f(error)
      call h5fopen_f(trim(fname), H5F_ACC_RDONLY_F, file_id, error)

      call h5ltget_attribute_double_f(file_id, "/", "time",  rbuf, error) ; time   = rbuf(1)
      call h5ltget_attribute_int_f   (file_id, "/", "npart", ibuf, error) ; n_part = ibuf(1)

      dimm = [min(n_part,n), ndims]

      call h5dopen_f(file_id, pvars, dataset_id, error)
      call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, table(1,1:dimm(1),:), dimm, error)
      call h5dclose_f(dataset_id, error)

      call h5dopen_f(file_id, vvars, dataset_id, error)
      call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, table(2,1:dimm(1),:), dimm, error)
      call h5dclose_f(dataset_id, error)

      call h5fclose_f(file_id, error)
      call h5close_f(error)

   end subroutine read_nbody_hdf5

end module particles_io_hdf5
