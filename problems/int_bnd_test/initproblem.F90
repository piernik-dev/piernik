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
!    Initial implementation of PIERNIK code was based on TVD split MHD code by
!    Ue-Li Pen
!        see: Pen, Arras & Wong (2003) for algorithm and
!             http://www.cita.utoronto.ca/~pen/MHD
!             for original source code "mhd.f90"
!
!    For full list of developers see $PIERNIK_HOME/license/pdt.txt
!
#include "piernik.h"
#include "macros.h"

module initproblem

   use constants, only: dsetnamelen

   implicit none

   private
   public :: read_problem_par, problem_initial_conditions, problem_pointers

   ! namelist parameters
!   namelist /PROBLEM_CONTROL/

   ! other private data
   character(len=dsetnamelen), parameter :: fld_n = "fld"

contains

!-----------------------------------------------------------------------------

   subroutine problem_pointers

      use dataio_user, only: user_vars_hdf5, user_reg_var_restart

      implicit none

      user_vars_hdf5       => fld_var
      user_reg_var_restart => fld_reg

   end subroutine problem_pointers

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      implicit none

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine problem_initial_conditions

      use cg_list_global,   only: all_cg
      use constants,        only: xdim, zdim
      use data_hdf5,        only: write_hdf5
      use dataio_pub,       only: printinfo, msg, die
      use mpisetup,         only: master
      use named_array_list, only: qna, wna

      implicit none

      integer(kind=4) :: dir

      call fld_reg

      call set_up_ub
      call diffuse_idn
      if (master) call printinfo("ip:ip ub+diff")
      call write_hdf5

      call all_cg%internal_boundaries_4d(wna%fi)
      call all_cg%internal_boundaries_3d(qna%ind(fld_n))
      if (master) call printinfo("ip:ip ub+diff+ib")
      call write_hdf5

      call diffuse_idn
      if (master) call printinfo("ip:ip ub+ib+diff")
      call write_hdf5

      call set_up_ub
      do dir = xdim, zdim
         call all_cg%internal_boundaries_4d(wna%fi, dir=dir)
         write(msg, '(a,i1)')"ip:ip ub+diff+ib:",dir
         if (master) call printinfo(msg)
      enddo
      call diffuse_idn
      if (master) call printinfo("ip:ip ub+ibXYZ+diff")
      call write_hdf5

!!$      call set_up_ub
!!$      do dir = xdim, zdim
!!$         call all_cg%internal_boundaries_4d(wna%fi, dim=dir)
!!$         call all_cg%internal_boundaries_3d(qna%ind(fld_n), dim=dir)
!!$         write(msg, '(a,i1)')"ip:ip ub+ib:",dir
!!$         if (master) call printinfo(msg)
!!$         call write_hdf5
!!$         call diffuse_idn
!!$         write(msg, '(a,i1,a)')"ip:ip ub+ib:",dir,"+diff"
!!$         if (master) call printinfo(msg)
!!$         call write_hdf5
!!$         call clear_fld
!!$      enddo

      call die("[initproblem:problem_initial_conditions] End of test")

   end subroutine problem_initial_conditions

   subroutine set_up_ub

      use cg_list,        only: cg_list_element
      use cg_list_global, only: all_cg
      use constants,      only: xdim, ydim, zdim
      use domain,         only: dom
      use fluidindex,     only: flind
      use grid_cont,      only: grid_container

      implicit none

      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg

      cgl => all_cg%first
      do while (associated(cgl))
         cg => cgl%cg

         cg%b(:, :, :, :) = 0.
         cg%u(:, :, :, :) = 0.

         if (dom%has_dir(xdim)) then
            cg%u(flind%neu%idn, :cg%is-dom%D_x,  :, :) = cg%u(flind%neu%idn, :cg%is-dom%D_x,  :, :) + 1.
            cg%u(flind%neu%idn,  cg%ie+dom%D_x:, :, :) = cg%u(flind%neu%idn,  cg%ie+dom%D_x:, :, :) + 10.
         endif

         if (dom%has_dir(ydim)) then
            cg%u(flind%neu%idn, :, :cg%js-dom%D_y,  :) = cg%u(flind%neu%idn, :, :cg%js-dom%D_y,  :) + 100.
            cg%u(flind%neu%idn, :,  cg%je+dom%D_y:, :) = cg%u(flind%neu%idn, :,  cg%je+dom%D_y:, :) + 1000.
         endif

         if (dom%has_dir(zdim)) then
            cg%u(flind%neu%idn, :, :, :cg%ks-dom%D_z ) = cg%u(flind%neu%idn, :, :, :cg%ks-dom%D_z ) + 10000.
            cg%u(flind%neu%idn, :, :,  cg%ke+dom%D_z:) = cg%u(flind%neu%idn, :, :,  cg%ke+dom%D_z:) + 100000.
         endif

         cgl => cgl%nxt
      enddo

   end subroutine set_up_ub

   subroutine diffuse_idn

      use cg_list,          only: cg_list_element
      use cg_list_global,   only: all_cg
      use domain,           only: dom
      use fluidindex,       only: flind
      use grid_cont,        only: grid_container
      use named_array_list, only: qna

      implicit none

      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg

      integer :: i, j, k

      call clear_fld
      cgl => all_cg%first
      do while (associated(cgl))
         cg => cgl%cg
         do k = cg%ks, cg%ke
            do j = cg%js, cg%je
               do i = cg%is, cg%ie
                  cg%q(qna%ind(fld_n))%arr(i, j, k) = sum(cg%u(flind%neu%idn, i-dom%D_x:i+dom%D_x, j-dom%D_y:j+dom%D_y, k-dom%D_z:k+dom%D_z))
               enddo
            enddo
         enddo
         cgl => cgl%nxt
      enddo

   end subroutine diffuse_idn

   subroutine clear_fld

      use cg_list_global,   only: all_cg
      use cg_list,          only: cg_list_element
      use named_array_list, only: qna

      implicit none

      type(cg_list_element), pointer :: cgl

      cgl => all_cg%first
      do while (associated(cgl))
         cgl%cg%q(qna%ind(fld_n))%arr(:,:,:) = 0.
         cgl => cgl%nxt
      enddo

   end subroutine clear_fld

   subroutine fld_var(var, tab, ierrh, cg)

      use grid_cont,   only: grid_container
      use named_array_list, only: qna

      implicit none

      character(len=*), intent(in)                    :: var
      real(kind=4), dimension(:,:,:), intent(inout)   :: tab
      integer, intent(inout)                          :: ierrh
      type(grid_container), pointer, intent(in)       :: cg

      ierrh = 0
      select case (trim(var))
         case ("fld")
            tab(:,:,:) = real(cg%q(qna%ind(fld_n))%span(cg%ijkse), 4)
         case default
            ierrh = -1
      end select

   end subroutine fld_var

   subroutine fld_reg

      use cg_list_global, only: all_cg
      use constants,      only: AT_OUT_B

      implicit none

      call all_cg%reg_var(fld_n, restart_mode = AT_OUT_B)

   end subroutine fld_reg

end module initproblem
