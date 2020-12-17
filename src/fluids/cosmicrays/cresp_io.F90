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
!! \brief This module contains I/O instructions and handlers for CRESP (NR) data.
!<
module cresp_io

! pulled by COSM_RAY_ELECTRONS

   use constants,    only: cbuff_len, LO, HI

   implicit none

   private
   public   :: check_file_group, file_has_group, file_has_dataset, save_smap_to_open, save_cresp_smap_h5, &
            &  read_cresp_smap_fields, create_cresp_smap_fields


   contains

!> \brief Create fields containing parameters of cresp solution maps, which are saved later in these groups.
!
      subroutine create_cresp_smap_fields(file_id)

      use constants,       only: LO, HI
      use cresp_helpers,   only: n_g_cresp, n_g_smaps, n_a_dims, n_a_esmall, n_a_max_p_r, &
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
!> \brief Create an external file to save ONLY cresp solution maps, optional, WIP TODO, FIXME
!
   subroutine save_cresp_smap_h5(smap_data, bound, dsname, filename)

      use cresp_helpers, only: n_g_smaps
      use dataio_pub,      only: msg, printinfo
      use hdf5,            only: HID_T, h5close_f, h5fclose_f, h5fcreate_f, h5fopen_f, h5gclose_f,       &
         &  h5gopen_f, h5open_f, H5F_ACC_RDWR_F, H5F_ACC_TRUNC_F
      use helpers_hdf5,    only: create_attribute

      implicit none

      integer,                      intent(in) :: bound
      character(len=*),             intent(in) :: dsname
      logical                                  :: file_exist, has_group, has_dset
      integer(HID_T)                           :: file_id
      character(len=*),             intent(in) :: filename
      integer                                     :: error
      real, dimension(:,:), target, intent(inout) :: smap_data

      call check_file_group(filename, "/cresp", file_exist, has_group) ! Assumption: if H5 file has "/cresp" group, other fields should be present

      if (.not. file_exist) then
         write(msg,"(A,A,A)") "[cresp_io:save_cresp_smap_h5] File '", trim(filename),"' does not exist; must be created and initialized."
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
         call save_smap_to_open(file_id, n_g_smaps(bound), dsname, smap_data)
      else
         write(msg,"(A,A,A)") "[cresp_io:save_cresp_smap_h5] Dataset ", "/"//n_g_smaps(bound)//"/"//dsname, " already present in file; omitting."
         call printinfo(msg)
      endif

   end subroutine save_cresp_smap_h5
!---------------------------------------------------------------------------------------------------
!> \brief Prepare and save rank-2 real array dataset in given file_id
!
   subroutine save_smap_to_open(file_id, path, dsname, dsdata)

      use helpers_hdf5,    only: create_dataset
      use hdf5,            only: HID_T, h5gopen_f, h5gclose_f

      implicit none

      real, dimension(:,:), target, intent(in) :: dsdata
      character(len=*),             intent(in) :: dsname, path
      integer(HID_T),               intent(in) :: file_id
      integer(HID_T)                           :: group_id
      integer                                  :: error
      real, dimension(:,:), pointer            :: pdsdata

      call h5gopen_f(file_id, path, group_id, error)
      pdsdata => dsdata
      call create_dataset(group_id, dsname, pdsdata)
      call h5gclose_f(group_id, error)
      pdsdata => null()


   end subroutine save_smap_to_open
!---------------------------------------------------------------------------------------------------
!> \brief Prepare and read cresp solution map fields from an external file.
!
   subroutine read_cresp_smap_fields(read_error, filename_opt)

      use constants,          only: cwdlen
      use cresp_helpers,      only: n_g_smaps, n_a_dims, n_a_esmall, n_a_max_p_r, n_a_clight,   &
        &  n_a_qbig, n_a_amin, n_a_amax, n_a_nmin, n_a_nmax, hdr_io, int_attrs, real_attrs
      use dataio_pub,         only: die, msg, printinfo
      use hdf5,               only: HID_T, H5F_ACC_RDONLY_F, h5close_f, h5fclose_f, h5open_f, h5fopen_f
      use set_get_attributes, only: get_attr

      implicit none

      integer                                    :: error, i, ia
      integer(HID_T)                             :: file_id
      character(len=cwdlen)                      :: filename
      character(len=*), optional,     intent(in) :: filename_opt
      logical                                    :: file_exist
      logical,                       intent(out) :: read_error
      integer(kind=4), dimension(:), allocatable :: ibuf
      real,            dimension(:), allocatable :: rbuf

      if (present(filename_opt)) then
         filename = filename_opt
      else
         filename = "CRESP_smaps.h5"   ! FIXME - I should not be here! TODO arrange the code to read smaps and header in common_hdf5
         write (msg,"(A,A)") "[cresp_io:read_cresp_smap_fields] Name of file to open: ", trim(filename)
         call printinfo(msg)
      endif

      inquire(file = trim(filename), exist = file_exist)

      if (.not. file_exist) then
         write(msg, "(A,A,A)") "[cresp_io:read_cresp_smap_fields] file '", trim(filename), "' does not exist."
         call printinfo(msg)
         read_error = .true.
         return
      else
         write(msg, "(A,A,A)") "[cresp_io:read_cresp_smap_fields] Proceeding with '", trim(filename), "' file."
         call printinfo(msg)
         read_error = .false.
      endif

      call h5open_f(error)
      call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, error)

! WARNING format version omitted - CRESP is not included in older versions of piernik

      do i = LO, HI              ! use statement for LO, HI in the upper level

         do ia = lbound(int_attrs, dim=1), ubound(int_attrs, dim=1)
            call get_attr(file_id, trim(int_attrs(ia)),  ibuf, n_g_smaps(i))
            if ((int_attrs(ia)) .eq. n_a_dims ) then
               hdr_io(i)%s_dim1    = ibuf(1)
               hdr_io(i)%s_dim2    = ibuf(2)
            else
               call die("[cresp_io:read_cresp_smap_fields] Non-recognized area_type.")
            endif
         enddo
         do ia = lbound(real_attrs, dim=1), ubound(real_attrs, dim=1)
            call get_attr(file_id, trim(real_attrs(ia)), rbuf, n_g_smaps(i))
            select case (real_attrs(ia))
               case (n_a_esmall)
                  hdr_io(i)%s_es      = rbuf(1)
               case (n_a_max_p_r)
                  hdr_io(i)%s_pr      = rbuf(1)
               case (n_a_qbig)
                  hdr_io(i)%s_qbig    = rbuf(1)
               case (n_a_clight)
                  hdr_io(i)%s_c       = rbuf(1)
               case (n_a_amin)
                  hdr_io(i)%s_amin    = rbuf(1)
               case (n_a_amax)
                  hdr_io(i)%s_amax    = rbuf(1)
               case (n_a_nmin)
                  hdr_io(i)%s_nmin    = rbuf(1)
               case (n_a_nmax)
                  hdr_io(i)%s_nmax    = rbuf(1)
               case default
                  call die("[cresp_io:read_cresp_smap_fields] Non-recognized area_type.")
            end select
         enddo
      enddo

      call h5fclose_f(file_id, error)
      call h5close_f(error)

   end subroutine read_cresp_smap_fields
!---------------------------------------------------------------------------------------------------
!> \brief Check if file "filename" exists and has group "groupname".
!
   subroutine check_file_group(filename, groupname, file_exist, has_group)

      implicit none

      character(len=*)             :: filename
      logical,       intent(inout) :: file_exist, has_group
      character(len=*), intent(in) :: groupname

      file_exist = .false.
      has_group  = .false.

      inquire(file = trim(filename), exist = file_exist)    ! check if file "filename" exists


      if (file_exist) then
         has_group = file_has_group(filename, groupname)
      endif

   end subroutine check_file_group
!---------------------------------------------------------------------------------------------------
!
!> \brief Check if HDF5 file "filename" has group "groupname" by trying to open it
!
   logical function file_has_group(filename, groupname) result (has_group)   ! WARNING may be slow!

      use constants,    only: I_ZERO, I_ONE
      use hdf5,         only: HID_T, H5F_ACC_RDONLY_F, h5close_f, h5fclose_f, h5fopen_f, &
                        &     h5gclose_f, h5gopen_f, h5eset_auto_f, h5open_f
      implicit none

      integer                       :: error
      integer(HID_T)                :: file_id, group_id
      character(len=*), intent(in)  :: filename, groupname

      has_group = .false.

      call h5open_f(error)
      call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, error)
      call h5eset_auto_f(I_ZERO, error)   ! Turn off printing messages
      call h5gopen_f(file_id, groupname, group_id, error)
      if (error .eq. 0) has_group = .true.
      call h5eset_auto_f(I_ONE, error)    ! Turn on  printing messages
      call h5gclose_f(group_id, error)
      call h5fclose_f(file_id, error)
      call h5close_f(error)

   end function file_has_group
!---------------------------------------------------------------------------------------------------
!
!> \brief Check if HDF5 file "filename" has dataset "dsname" by trying to open it
!
   logical function file_has_dataset(filename, dsname) result (has_dset)   ! WARNING may be slow!

      use constants,    only: I_ZERO, I_ONE
      use hdf5,         only: HID_T, H5F_ACC_RDONLY_F, h5close_f, h5dopen_f, h5fclose_f, &
                        &  h5fopen_f, h5dclose_f, h5eset_auto_f, h5open_f
      implicit none

      integer                       :: error
      integer(HID_T)                :: file_id, dset_id
      character(len=*), intent(in)  :: filename, dsname

      has_dset = .false.

      call h5open_f(error)
      call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, error)
      call h5eset_auto_f(I_ZERO, error)   ! Turn off printing messages
      call h5dopen_f(file_id, dsname, dset_id, error)
      if (error .eq. 0) has_dset = .true.
      call h5dclose_f(dset_id, error)
      call h5eset_auto_f(I_ONE, error)    ! Turn on  printing messages
      call h5fclose_f(file_id, error)
      call h5close_f(error)

   end function file_has_dataset
!---------------------------------------------------------------------------------------------------


end module cresp_io
