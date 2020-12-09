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
   public   :: gdf_create_cresp_smap_fields

   contains

!  TODO include creation / writing to a new file when basic functionalities are working

   subroutine gdf_create_cresp_smap_fields(file_id)

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

   end subroutine gdf_create_cresp_smap_fields

end module cresp_io_write
