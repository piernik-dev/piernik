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
!! \brief This module contains names and types useful in writing / reading data for CRESP.
!<

module cresp_helpers
! pulled by COSM_RAY_ELECTRONS

   use constants,    only: cbuff_len, LO, HI

   implicit none

   private
   public   :: hdr_io, map_header, n_g_cresp, n_g_smaps, n_a_dims, n_a_esmall, n_a_max_p_r, n_a_clight,     &
      &  n_a_qbig, n_a_amin, n_a_amax, n_a_nmin, n_a_nmax, real_attrs, int_attrs, extension, flen,          &
      &  bound_name, dset_attrs, enden_CMB

   character(len=*), parameter, dimension(LO:HI)      ::  n_g_smaps = [ "cresp/smaps_LO", "cresp/smaps_UP" ]
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

   integer, parameter                               :: extlen = 4, flen = 15
   character(len=extlen), parameter                 :: extension =  ".dat"

   type     map_header
      integer(kind=4) :: s_dim1, s_dim2
      real     :: s_es
      real     :: s_pr
      real     :: s_qbig
      real     :: s_c
      real     :: s_amin, s_amax, s_nmin, s_nmax
   end type map_header

   type(map_header), dimension(LO:HI)  :: hdr_io

   contains

!> \brief Calculate energy density of Cosmic Microwave Background at given epoch

   elemental real function enden_CMB(z)
      use units,     only: u_CMB
      use constants, only: one, four
      implicit none

      real, intent(in)  :: z

      enden_CMB = u_CMB * (one + z)**four

   end function enden_CMB


end module cresp_helpers
