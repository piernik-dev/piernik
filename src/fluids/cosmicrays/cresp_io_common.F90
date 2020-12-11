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
!! \brief This module contains I/O instructions for writing / reading data for CRESP.
!<

module cresp_io_common
! pulled by COSM_RAY_ELECTRONS

   use constants,    only: cbuff_len, LO, HI

   implicit none

   private
   public   :: hdr_io, map_header, n_g_cresp, n_g_smaps, n_a_dims, n_a_esmall, n_a_max_p_r, n_a_clight,     &
      &  n_a_qbig, n_a_amin, n_a_amax, n_a_nmin, n_a_nmax, real_attrs, int_attrs, create_dataset_real8_dim2,&
      &  bound_name, dset_attrs

   character(len=*), parameter, dimension(LO:HI)      ::   n_g_smaps = (/ "cresp/smaps_LO", "cresp/smaps_UP" /)
   character(len=*), parameter :: n_g_cresp = "cresp", &
      &  n_a_dims = "dims", n_a_esmall = "e_small", n_a_max_p_r = "max_p_ratio", n_a_clight = "used_clight", &
      &  n_a_qbig = "q_big",  n_a_amin = "a_min", n_a_amax = "a_max", n_a_nmin = "n_min", n_a_nmax = "n_max"
   character(len=cbuff_len), parameter, dimension(8)  :: real_attrs = [ "e_small    ",  &
                                                         &              "max_p_ratio",  &
                                                         &              "used_clight",  &
                                                         &              "q_big      ",  &
                                                         &              "a_min      ",  &
                                                         &              "a_max      ",  &
                                                         &              "n_min      ",  &
                                                         &              "n_max      "   ]

   character(len=cbuff_len), parameter, dimension(2)  :: dset_attrs = [ "p_ratios",     &
                                                         &              "f_ratios"      ]
   character(len=cbuff_len), parameter, dimension(1)  :: int_attrs =  [ "dims       "   ]

   integer, parameter                                 :: blen = 2

   character(len=blen), dimension(LO:HI), parameter   :: bound_name = ['lo', 'up']

   type     map_header
      integer  :: s_dim1, s_dim2
      real     :: s_es
      real     :: s_pr
      real     :: s_qbig
      real     :: s_c
      real     :: s_amin, s_amax, s_nmin, s_nmax
   end type map_header

   type(map_header), dimension(LO:HI)  :: hdr_io

   contains

!
!> \brief Create double precision real dataset (rank-2 array) in the given place.
!
   subroutine create_dataset_real8_dim2(place, dname, ddata)   ! TODO might be moved to helpers_hdf5:create_dataset interface

      use constants,     only: I_TWO
      use hdf5,          only: h5dcreate_f, h5dclose_f, h5kind_to_type, h5sclose_f, h5screate_simple_f, &
          &                    h5dwrite_f, HID_T, HSIZE_T, H5T_NATIVE_DOUBLE
      use iso_c_binding, only: c_ptr, c_loc

      implicit none

      integer(HID_T),                intent(in) :: place !< object id where dataset will be created
      character(len=*),              intent(in) :: dname !< name of dataset
      real, pointer, dimension(:,:), intent(in) :: ddata !< data used to create dataset

      integer(HID_T)                            :: dset, space!, mem_type
      integer(kind=4)                           :: hdferr
      integer(HSIZE_T), dimension(2)            :: dims
      type(c_ptr)                               :: f_ptr

      dims = shape(ddata)
      call h5screate_simple_f(I_TWO, dims, space, hdferr)
      call h5dcreate_f(place, dname, H5T_NATIVE_DOUBLE, space, dset, hdferr)
      f_ptr = c_loc(ddata(1,1))
      call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, f_ptr, hdferr)
      call h5dclose_f(dset,  hdferr)
      call h5sclose_f(space, hdferr)

   end subroutine create_dataset_real8_dim2

end module cresp_io_common
