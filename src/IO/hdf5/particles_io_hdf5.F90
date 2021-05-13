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

      use common_hdf5,    only: set_common_attributes
      use constants,      only: cwdlen, tmr_hdf, idlen
      use cg_particles_io, only: dump_cg_particles
      use dataio_pub,     only: msg, printinfo, printio, thdf
      use hdf5,           only: h5open_f, h5close_f, h5fcreate_f, h5fclose_f, HID_T, H5F_ACC_TRUNC_F
      use mpisetup,       only: master, proc
      use timer,          only: set_timer

      implicit none

      character(len=*), intent(in) :: fname
      character(len=cwdlen)        :: filename
      character(len=idlen)         :: proc_c
      integer                      :: flen
      integer(kind=4)              :: error
      integer(HID_T)               :: file_id

      thdf = set_timer(tmr_hdf,.true.)
      flen = len(trim(fname))
      write(proc_c,'(i3.3)') proc
      write(filename,'(6a)') fname(1:flen-7), 'p', fname(flen-6:flen-3),'_',proc_c,'.h5'
      write(msg,'(a,1x,2a)') 'Writing datafile ', trim(filename), " ... "
      call printio(msg, .true.)
      call set_common_attributes(filename)
      call h5open_f(error)
      call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

      call set_nbody_attributes(file_id)
      call dump_cg_particles(file_id)

      call h5fclose_f(file_id, error)
      call h5close_f(error)

      thdf = set_timer(tmr_hdf)
      if (master) then
         write(msg,'(a6,f10.2,a2)') ' done ', thdf, ' s'
         call printinfo(msg, .true.)
      endif

   end subroutine write_nbody_hdf5


   subroutine set_nbody_attributes(file_id)

      use dataio_pub,         only: require_problem_IC, piernik_hdf5_version2, problem_name, run_id, last_hdf_time, &
         &                          last_res_time, last_log_time, last_tsl_time, nres, nhdf, domain_dump
      use global,             only: t, dt, nstep
      use hdf5,               only: HID_T
      use particle_utils,     only: count_all_particles
      use set_get_attributes, only: set_attr

      implicit none

      integer(HID_T),  intent(in) :: file_id       !< File identifier
      integer(kind=4)             :: n_part

      n_part = count_all_particles()
      ! real attributes
      call set_attr(file_id, "time",          [t                     ]) !rr2
      call set_attr(file_id, "timestep",      [dt                    ]) !rr2
      call set_attr(file_id, "piernik",       [piernik_hdf5_version2 ]) !rr1, rr2
      call set_attr(file_id, "last_log_time", [last_log_time         ]) !rr2
      call set_attr(file_id, "last_tsl_time", [last_tsl_time         ]) !rr2
      call set_attr(file_id, "last_hdf_time", [last_hdf_time         ]) !rr2
      call set_attr(file_id, "last_res_time", [last_res_time         ]) !rr2

      ! integer attributes
      call set_attr(file_id, "npart",              [n_part            ]) !rr2
      call set_attr(file_id, "nstep",              [nstep             ]) !rr2
      call set_attr(file_id, "nres",               [nres              ]) !rr2
      call set_attr(file_id, "nhdf",               [nhdf              ]) !rr2
      call set_attr(file_id, "require_problem_IC", [require_problem_IC]) !rr2
      !> \todo  add number of pieces in the restart point/data dump

      ! string attributes
      call set_attr(file_id, "problem_name", [trim(problem_name)]) !rr2
      call set_attr(file_id, "domain",       [trim(domain_dump) ]) !rr2
      call set_attr(file_id, "run_id",       [trim(run_id)      ]) !rr2

   end subroutine set_nbody_attributes

   subroutine read_nbody_hdf5(fname, table, n)

      use constants, only: idlen, ndims
      use hdf5,      only: h5open_f, h5close_f, h5fopen_f, h5fclose_f, h5dopen_f, h5dclose_f, h5dread_f
      use hdf5,      only: HID_T, HSIZE_T, H5F_ACC_RDONLY_F, H5T_NATIVE_DOUBLE
      use h5lt,      only: h5ltget_attribute_double_f, h5ltget_attribute_int_f

      implicit none

      character(len=*),           intent(in)  :: fname
      integer(kind=4),            intent(in)  :: n
      real, dimension(2,n,ndims), intent(out) :: table
      character(len=idlen)                    :: pvars = 'pos', vvars = 'vel'
      integer(kind=4)                         :: n_part
      integer(kind=4)                         :: error
      integer(HID_T)                          :: file_id, dataset_id
      integer(HSIZE_T), dimension(2)          :: dimm
      integer(kind=4), dimension(1)           :: ibuf
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
