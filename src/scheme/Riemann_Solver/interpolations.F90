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
!--------------------------------------------------------------------------------------------------------------

#include "piernik.def"

module interpolations
! pulled by RIEMANN

  implicit none
  private
  public :: set_interpolations, interpol

  interface
     subroutine interpolation(n,f,fl,fr)
       
       implicit none

       integer(kind=4),      intent(in)  :: n
       real, dimension(:,:), intent(in)  :: f
       real, dimension(:,:), intent(out) :: fl
       real, dimension(:,:), intent(out) ::fr

    end interface

    procedure(interpolation), pointer :: interpol => null()

  contains


    subroutine set_interpolations(interpol_str)

      use dataio_pub, only: die

      implicit none

      character(len=*), intent(in) :: interpol_str

      if (associated(interp)) call die("[interpolations:set_interpolations] interpol already associated")
      interpol => set_interpolation(interpol_str)
      
    end subroutine set_interpolations

    function set_interpolation(interp_str) result(interp)

      use dataio_pub, only: msg, die

      implicit none

      character(len=*), intent(in) :: interp_str

      procedure(interpolation), pointer :: interp

      select case(interp_str)
         case('linear', 'LINEAR')
            inetrp => linear_interpolation
         case default
            write(msg,'(2a)') "[interpolations:set_interpolation] unknown interpolation ", interp_str 
            call die(msg)
            inter => null()
         end select

    end function set_interpolation
      
    
end module interpolations
