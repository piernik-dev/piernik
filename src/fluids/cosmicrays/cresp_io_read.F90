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
!! \brief This module contains I/O instructions for reading data with CRESP.
!<

module   cresp_io_read
! pulled by COSM_RAY_ELECTRONS

   use constants,    only: LO, HI

   implicit none

   private
   public   :: read_cresp_smap_fields, read_cresp_smap_h5

   contains

   subroutine read_cresp_smap_fields(read_error)

      use common_hdf5,        only: output_fname
      use constants,          only: cwdlen, RD
      use cresp_io_common,    only: n_g_smaps, n_a_dims, n_a_esmall, n_a_max_p_r, n_a_clight,   &
        &  n_a_qbig, n_a_amin, n_a_amax, n_a_nmin, n_a_nmax, hdr_io, int_attrs, real_attrs
      use dataio_pub,         only: die, msg, nres, printinfo
      use hdf5,               only: HID_T, H5F_ACC_RDONLY_F, h5close_f, h5fclose_f, h5open_f, h5fopen_f
      use set_get_attributes, only: get_attr

      implicit none

      integer                         :: error, i, ia
      integer(HID_T)                  :: file_id
      character(len=cwdlen)           :: filename
      logical                         :: file_exist
      logical, intent(inout)          :: read_error
      integer(kind=4), dimension(:), allocatable :: ibuf
      real,            dimension(:), allocatable :: rbuf

      filename = output_fname(RD, '.res', nres+1, bcast=.true.) ! TRY opening res with number > 0
      write (msg,"(A,A)") "[cresp_io_read:read_cresp_smap_fields] Name of file to open: ", trim(filename)
      call printinfo(msg)

      inquire(file = trim(filename), exist = file_exist)

      if (.not. file_exist) then
         write(msg, "(A,A,A)") "[cresp_io_common:read_cresp_smap_fields] file '", trim(filename), "' does not exist."
         call printinfo(msg)
         read_error = .true.
         return
      else
         write(msg, "(A,A,A)") "[cresp_io_common:read_cresp_smap_fields] Proceeding with '", trim(filename), "' file."
         call printinfo(msg)
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
               call die("[cresp_io_read:read_cresp_smap_fields] Non-recognized area_type.")
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
                  call die("[cresp_io_read:read_cresp_smap_fields] Non-recognized area_type.")
            end select
         enddo
      enddo

      call h5fclose_f(file_id, error)
      call h5close_f(error)

   end subroutine read_cresp_smap_fields
!---------------------------------------------------------------------------------------------------
!
!> \brief Opens h5 file, checks if dset is present and reads it
!  TODO optimize me: many overlapping open / close instructions.
!
   subroutine read_cresp_smap_h5(filename, dsetname, dset, read_error)

      use cresp_io_common, only: file_has_dataset
      use dataio_pub,      only: msg, warn
      use hdf5,            only: HID_T, h5fclose_f, h5close_f, h5fopen_f, &
                        &     h5open_f, H5F_ACC_RDONLY_F
      implicit none

      character(len=*),                  intent(in) :: dsetname
      real,   dimension(:,:), target, intent(inout) :: dset
      integer                                       :: error
      integer(HID_T)                                :: file_id
      character(len=*),                  intent(in) :: filename
      logical                                       :: has_dataset
      real, dimension(:,:), pointer                 :: pdset
      logical,                          intent(out) :: read_error

      pdset => null()
      has_dataset = file_has_dataset(filename, dsetname)
      read_error  = .not. has_dataset

      if (has_dataset) then
         call h5open_f(error)
         call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, error)

         pdset => dset(:,:)
         call read_real_arr2d_dset(file_id, dsetname, pdset)
         pdset => null()

         call h5fclose_f(file_id, error)
         call h5close_f(error)
      else
         write(msg, "(5a)") "[cresp_io_read:read_cresp_smap_h5] Problem reading dataset '", dsetname, "' from file '", trim(filename),"'."
         call warn(msg)
      endif

   end subroutine read_cresp_smap_h5
!---------------------------------------------------------------------------------------------------
!
!> \brief  Read 2D double precision dataset 'dsdata' of provided name 'dsetname' with 'file_id'
!  Requires file to be open, may be slow. WARNING Does not contain fail-safe instructions!
!
   subroutine read_real_arr2d_dset(file_id, dsetname, pa2d)

      use hdf5,      only: HID_T, h5dopen_f, h5dread_f, h5dopen_f, &
                     &  h5dclose_f, H5T_NATIVE_DOUBLE, HSIZE_T
      implicit none

      integer(HID_T)                 :: dset_id
      integer(HID_T),     intent(in) :: file_id
      integer                        :: error
      character(len=*),   intent(in) :: dsetname
      real, dimension(:,:),  pointer :: pa2d
      integer(HSIZE_T), dimension(2) :: pa2ddims

      if (associated(pa2d)) then
         pa2ddims = shape(pa2d)

         call h5dopen_f(file_id, dsetname, dset_id, error)
         call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, pa2d, pa2ddims, error)
         call h5dclose_f(dset_id, error)
      endif

   end subroutine read_real_arr2d_dset

end module cresp_io_read
