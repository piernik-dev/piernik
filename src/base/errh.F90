! $Id$
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
!    Initial implemetation of PIERNIK code was based on TVD split MHD code by
!    Ue-Li Pen 
!        see: Pen, Arras & Wong (2003) for algorithm and
!             http://www.cita.utoronto.ca/~pen/MHD 
!             for original source code "mhd.f90" 
!   
!    For full list of developers see $PIERNIK_HOME/license/pdt.txt
!
#include "piernik.def"
!>
!! \brief (KK) Module responsible for error handling 
!!
!<
module errh

   implicit none
   contains

   subroutine namelist_errh(ierrh,nm)
      implicit none
      integer, intent(in) :: ierrh
      character(len=*), intent(in) :: nm
      
      select case (ierrh)
         case (19)
            write(*,*) "severe (19): Invalid reference to variable in ",trim(nm), " namelist"
            write(*,*) "One of the following conditions occurred: "
            write(*,*) "    * The variable was not a member of the namelist group."
            write(*,*) "    * An attempt was made to subscript a scalar variable."
            write(*,*) "    * A subscript of the array variable was out-of-bounds."
            write(*,*) "    * An array variable was specified with too many or too few subscripts for the variable."
            write(*,*) "    * An attempt was made to specify a substring of a noncharacter variable or array name."
            write(*,*) "    * A substring specifier of the character variable was out-of-bounds."
            write(*,*) "    * A subscript or substring specifier of the variable was not an integer constant."
            write(*,*) "    * An attempt was made to specify a substring by using an unsubscripted array variable."
            stop
         case (-1)
            write(*,*) "Namelist: ",trim(nm)," not found in problem.par"
            stop
         case (5010)
            write(*,*) "One of the variables found in problem.par doesn't belong to ",trim(nm), " namelist"
            stop
         case (239)
            write(*,*) "One of the variables found in problem.par doesn't belong to ",trim(nm), " namelist"
            stop
         case (0)
         case default
            write(*,*) 'Unknown error (', ierrh,') in namelist ',trim(nm)
            stop
      endselect

   end subroutine namelist_errh

!-----------------------------------------------------------------------------
end module errh
