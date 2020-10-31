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

!> \brief This module provides a function to find a pointer to a level based on level ID

module find_lev

   implicit none

   private
   public :: find_level

contains

!> \brief Find pointer to a level given by level_id

   function find_level(level_id) result(lev_p)

      use cg_level_base,      only: base
      use cg_level_connected, only: cg_level_connected_t
      use dataio_pub,         only: die

      implicit none

      integer(kind=4), intent(in) :: level_id         !< level number (relative to base level)

      type(cg_level_connected_t), pointer :: lev_p, curl

      nullify(lev_p)
      curl => base%level
      if (level_id >= base%level%l%id) then
         do while (associated(curl))
            if (curl%l%id == level_id) then
               lev_p => curl
               exit
            endif
            curl => curl%finer
         enddo
      else
         do while (associated(curl))
            if (curl%l%id == level_id) then
               lev_p => curl
               exit
            endif
            curl => curl%coarser
         enddo
      endif

      if (.not. associated(lev_p)) call die("[find_level:find_level] Cannot find level")

   end function find_level

end module
