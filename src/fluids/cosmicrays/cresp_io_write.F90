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
!! \brief This module contains I/O instructions for writing data with CRESP.
!<

module cresp_io_write
! pulled by COSM_RAY_ELECTRONS

   private
   public   :: create_cresp_smap_fields, save_cresp_smap_h5

   contains

!  TODO include creation / writing to a new file when basic functionalities are working

   subroutine create_cresp_smap_fields(file_id)

      use constants,       only: LO, HI
      use cresp_io_common, only: n_g_cresp, n_g_smaps, n_a_dims, n_a_esmall, n_a_max_p_r, &
         &  n_a_clight, n_a_qbig, n_a_qbig, n_a_amin, n_a_amax, n_a_nmin, n_a_nmax, hdr_io
      use hdf5,            only: HID_T, h5gcreate_f, h5gclose_f
      use helpers_hdf5,    only: create_attribute

      implicit none

      integer(HID_T), intent(in) :: file_id
      integer(HID_T)             :: g_id, g_id_cresp
      integer                    :: error, i

      call h5gcreate_f(file_id, n_g_cresp, g_id_cresp, error)
         do i = LO, HI
            call h5gcreate_f(file_id, n_g_smaps(i), g_id, error)
               call create_attribute(g_id, n_a_dims,    [ hdr_io(i)%s_dim1, hdr_io(i)%s_dim2 ] )
               call create_attribute(g_id, n_a_esmall,  [ hdr_io(i)%s_es   ])
               call create_attribute(g_id, n_a_max_p_r, [ hdr_io(i)%s_pr   ])
               call create_attribute(g_id, n_a_qbig,    [ hdr_io(i)%s_qbig ])
               call create_attribute(g_id, n_a_clight,  [ hdr_io(i)%s_c    ])
               call create_attribute(g_id, n_a_amin,    [ hdr_io(i)%s_amin ])
               call create_attribute(g_id, n_a_amax,    [ hdr_io(i)%s_amax ])
               call create_attribute(g_id, n_a_nmin,    [ hdr_io(i)%s_nmin ])
               call create_attribute(g_id, n_a_nmax,    [ hdr_io(i)%s_nmax ])
            call h5gclose_f(g_id, error)
         enddo
      call h5gclose_f(g_id_cresp, error)

   end subroutine create_cresp_smap_fields
!---------------------------------------------------------------------------------------------------
   subroutine save_cresp_smap_h5(smap_data, bound, dsname, filename)

      use cresp_io_common, only: n_g_smaps, create_dataset_real8_dim2, check_file_group, file_has_group, &
         &  file_has_dataset
      use dataio_pub,      only: msg, printinfo
      use hdf5,            only: HID_T, h5close_f, h5fclose_f, h5fcreate_f, h5fopen_f, h5gclose_f,       &
         &  h5gopen_f, h5open_f, H5F_ACC_RDWR_F, H5F_ACC_TRUNC_F
      use helpers_hdf5,    only: create_attribute

      implicit none

      integer,                      intent(in) :: bound
      character(len=*),             intent(in) :: dsname
      logical                                  :: file_exist, has_group, has_dset
      integer(HID_T)                           :: file_id, g_id
      character(len=*),             intent(in) :: filename
      integer                                     :: error
      real, dimension(:,:), target, intent(inout) :: smap_data
      real(kind=8), pointer, dimension(:,:)       :: p_smap

      call check_file_group(filename, "/cresp", file_exist, has_group) ! Assumption: if H5 file has "/cresp" group, other fields should be present

      if (.not. file_exist) then
         write(msg,"(A,A,A)") "[cresp_io_write:save_cresp_smap_h5] File '", trim(filename),"' does not exist; must be created and initialized."
         call printinfo(msg)

         call h5open_f(error)
         call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, file_id, error)
         call h5fclose_f(file_id, error)
         call h5close_f(error)
      endif

      if (.not. has_group) then
         call h5open_f(error)
         call h5fopen_f(trim(filename), H5F_ACC_RDWR_F, file_id, error)
         call create_cresp_smap_fields(file_id)   ! Creates "/cresp" group and associated attrs
         call h5fclose_f(file_id, error)
         call h5close_f(error)
      endif

      has_dset = file_has_dataset(filename, "/"//n_g_smaps(bound)//"/"//dsname)

      if (.not. has_dset) then
         call h5open_f(error)
         call h5fopen_f(trim(filename), H5F_ACC_RDWR_F, file_id, error)
         call h5gopen_f(file_id, n_g_smaps(bound), g_id, error)
         p_smap => smap_data
         call create_dataset_real8_dim2(g_id, dsname, p_smap)
         call h5gclose_f(g_id, error)
         call h5fclose_f(file_id, error)
         p_smap => null()
      else
         write(msg,"(A,A,A)") "[cresp_io_write:save_cresp_smap_h5] Dataset ", "/"//n_g_smaps(bound)//"/"//dsname, " already present in file; omitting."
         call printinfo(msg)
      endif

   end subroutine save_cresp_smap_h5

end module cresp_io_write
