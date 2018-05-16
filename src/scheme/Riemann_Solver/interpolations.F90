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
  public :: set_interpolations, interpol, linear

  interface
     subroutine interpolation(prim_var,prim_var_l,prim_var_r)
       
       implicit none
       
       real, dimension(:,:), intent(in)  :: prim_var
       real, dimension(:,:), intent(out) :: prim_var_l
       real, dimension(:,:), intent(out) :: prim_var_r

     end subroutine interpolation

  end interface

  procedure(interpolation), pointer :: interpol => null()

contains


  subroutine set_interpolations(interpol_str)

    use dataio_pub, only: die

    implicit none

    character(len=*), intent(in) :: interpol_str

    if (associated(interpol)) call die("[interpolations:set_interpolations] interpol already associated")
    interpol => set_interpolation(interpol_str)
      
  end subroutine set_interpolations

  function set_interpolation(interp_str) result(interp)

    use dataio_pub, only: msg, die

    implicit none

    character(len=*), intent(in) :: interp_str
    
    procedure(interpolation), pointer :: interp

    select case(interp_str)
    case('linear', 'LINEAR')
       interp => linear
    case default
       write(msg,'(2a)') "[interpolations:set_interpolation] unknown interpolation ", interp_str 
       call die(msg)
       interp => null()
    end select
    
  end function set_interpolation

  subroutine linear(q,ql,qr)
    
    use fluxlimiters, only: set_limiters,flimiter,blimiter
    use domain,       only: dom
    use constants,    only: half, GEO_XYZ
    use dataio_pub,   only: die
      
    implicit none

    real, dimension(:,:), intent(in)     :: q
    real, dimension(:,:), intent(out)    :: ql
    real, dimension(:,:), intent(out)    :: qr

    real, dimension(size(q,1),size(q,2)) :: dq_lim, dq_interp

    integer                              :: im1
    integer                              :: n 
    n = size(q,2)
      
    dq_lim = flimiter(q)

    if (dom%geometry_type /= GEO_XYZ) call die("[interpolations:linear] non-cartesian geometry not implemented yet.")

    dq_interp = half*dq_lim
    
    ql = q + dq_interp
    qr = q - dq_interp

    im1 = max(1,n-1) ! neighbouring indices
    !associate(im1 => i - Dom%D_x)
    ! shfit right state
    qr(:,1:im1) = qr(:,2:n)


    ! interpolation for the first and last points
    
    ql(:,1) = q(:,1)
    qr(:,n) = q(:,n)  
      
  end subroutine linear
    
end module interpolations
