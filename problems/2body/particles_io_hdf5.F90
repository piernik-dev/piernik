#include "piernik.h"
module particles_io_hdf5

   implicit none
   private
   public  :: write_hdf5, read_hdf5
   contains

#ifdef HDF5
   subroutine write_hdf5(table, n)

      use hdf5, only: h5open_f, h5close_f, h5fclose_f, h5dcreate_f, h5dclose_f, h5dwrite_f
      use hdf5, only: h5screate_simple_f, h5sclose_f, h5fcreate_f
      use hdf5, only: HID_T, HSIZE_T, SIZE_T, H5F_ACC_TRUNC_F, H5T_NATIVE_DOUBLE
      use h5lt, only: h5ltset_attribute_double_f, h5ltset_attribute_int_f
      implicit none
      
      character(len=128)                 :: hdf_name

      integer(kind=4)                    :: error, rank, iv
      integer(HID_T)                     :: file_id, dataspace_id, dataset_id
      integer(SIZE_T)                    :: bufsize
      integer(HSIZE_T), dimension(2)     :: dimm
      integer(kind=4) :: time=2

      !double precision, dimension(npart) :: temp
      !character(len=5), dimension(3)     :: vars
      character(len=9) :: vars='positions'
      character(len=7) :: dumpfile = 'test_01'
      
      
      
      integer :: n
      real, dimension(n, 3) :: table
      
      rank    = 2
      dimm(1) = n
      dimm(2) = 3
      bufsize = 1

      write(hdf_name,'(a,a3)') trim(dumpfile),'.h5'

      call h5open_f(error)
      call h5fcreate_f(hdf_name, H5F_ACC_TRUNC_F, file_id, error)

      call h5ltset_attribute_int_f   (file_id, "/", "time", [time], bufsize, error)
      !call h5ltset_attribute_double_f(file_id, "/", "time",     [t],  bufsize, error)

         call h5screate_simple_f(rank, dimm, dataspace_id, error)
         call h5dcreate_f(file_id, vars, H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, error)
         call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, table, dimm, error)
         call h5dclose_f(dataset_id, error)
         call h5sclose_f(dataspace_id, error)


      call h5fclose_f(file_id, error)
      call h5close_f(error)

   end subroutine write_hdf5


   subroutine read_hdf5(table,n)

      use hdf5,  only: h5open_f, h5close_f, h5fopen_f, h5fclose_f, h5dopen_f, h5dclose_f, h5dread_f, h5dget_space_f
      use hdf5,  only: HID_T, HSIZE_T, H5F_ACC_RDONLY_F, H5T_NATIVE_DOUBLE
      use h5lt,  only: h5ltget_attribute_double_f, h5ltget_attribute_int_f
      !use start, only: nstep
      implicit none


      character(len=128)                       :: hdf_name='test_01.h5'

      integer(kind=4)                          :: error, iv,time
      integer(HID_T)                           :: file_id, dataset_id
      integer(HSIZE_T), dimension(2)           :: dimm
      !character(len=5), dimension(4)           :: vars
      character(len=9)::  vars='positions'
      !
      integer,dimension(1) :: rbuf
      integer :: n
      real, dimension(n, 3),intent(out) :: table
      dimm(1) = n
      dimm(2) = 3

      call h5open_f(error)
      call h5fopen_f(hdf_name, H5F_ACC_RDONLY_F, file_id, error)

      call h5ltget_attribute_int_f  (file_id,  "/", "time", rbuf, error)
      time = rbuf(1)

      write(*,*) "time=", time


         call h5dopen_f(file_id, vars, dataset_id, error)
         call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, table, dimm, error)
         call h5dclose_f(dataset_id, error)

      !enddo
      !deallocate(temp)

      call h5fclose_f(file_id, error)
      call h5close_f(error)

   end subroutine read_hdf5
#endif /* HDF5 */

end module particles_io_hdf5
