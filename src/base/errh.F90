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

   private

   public :: die, warn, printinfo, namelist_errh

   include 'mpif.h'

   contains

!-----------------------------------------------------------------------------

   subroutine printinfo(nm, to_stdout)

      use dataio_public, only : colormessage, T_PLAIN, T_SILENT, T_INFO

      implicit none

      character(len=*), intent(in) :: nm
      logical, optional, intent(in) :: to_stdout

      if (present(to_stdout)) then
         if (to_stdout) then
            call colormessage(nm, T_PLAIN)
         else
            call colormessage(nm, T_SILENT)
         end if
      else
         call colormessage(nm, T_INFO)
      end if

   end subroutine printinfo

!-----------------------------------------------------------------------------

   subroutine warn(nm)

      use dataio_public, only : colormessage, T_WARN

      implicit none

      character(len=*), intent(in) :: nm

      call colormessage(nm, T_WARN)

   end subroutine warn

!-----------------------------------------------------------------------------
   !! BEWARE: routine is not finished, it should kill PIERNIK gracefully

   subroutine die(nm, allprocs)

      use dataio_public, only : colormessage, T_ERR

      implicit none

      character(len=*), intent(in)  :: nm
      integer, optional, intent(in) :: allprocs

      integer :: ierr

      call colormessage(nm, T_ERR)

      if(present(allprocs)) then
         if(allprocs /= 0) then
            call MPI_Barrier(MPI_COMM_WORLD, ierr)
            call MPI_Finalize(ierr)
         endif
      endif

      stop

   end subroutine die

!-----------------------------------------------------------------------------

   subroutine namelist_errh(ierrh,nm)

      use dataio_public, only : msg

      implicit none

      integer, intent(in) :: ierrh
      character(len=*), intent(in) :: nm

      select case (ierrh)
         case (19)
            call warn("One of the following conditions occurred: ")
            call warn("    * The variable was not a member of the namelist group.")
            call warn("    * An attempt was made to subscript a scalar variable.")
            call warn("    * A subscript of the array variable was out-of-bounds.")
            call warn("    * An array variable was specified with too many or too few subscripts for the variable.")
            call warn("    * An attempt was made to specify a substring of a noncharacter variable or array name.")
            call warn("    * A substring specifier of the character variable was out-of-bounds.")
            call warn("    * A subscript or substring specifier of the variable was not an integer constant.")
            call warn("    * An attempt was made to specify a substring by using an unsubscripted array variable.")
            call die("severe (19): Invalid reference to variable in the "//trim(nm)//" namelist")
         case (-1)
            call die("Namelist: "//trim(nm)//" not found in problem.par")
         case (239, 5010)
            call die("One of the variables found in problem.par doesn't belong to the "//trim(nm)//" namelist")
         case (0)
         case default
            write(msg, *)'Unknown error (', ierrh,') in namelist ',trim(nm)
            call die(msg)
      endselect

   end subroutine namelist_errh

!-----------------------------------------------------------------------------
end module errh
