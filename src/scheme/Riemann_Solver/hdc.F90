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
!    Hyperbolic divergence cleaning
!    A. Dedner et al., Journal of Computational Physics 175, 645-673 (2002)
!    Dr. Varadarajan Parthasarathy, Dr. Artur Gawryszczak, CAMK (Warszawa)
!    November 2017
!
#include "piernik.h"

#ifdef GLM

module hdc
! pulled by RIEMANN
  implicit none

  real, protected :: chspeed

  private
  public :: chspeed, update_chspeed, init_psi, glm_mhd

contains

!>
!! \brief recalcualte the speed of propagation of psi waves
!!
!! This can be perhaps evaluated somewhere in timestep and reused here.
!<

   subroutine update_chspeed()

      use cg_leaves,  only: leaves
      use cg_list,    only: cg_list_element
      use constants,  only: GEO_XYZ, small
      use dataio_pub, only: warn, die
      use domain,     only: dom
      use fluidindex, only: flind
      use fluids_pub, only: has_ion
      use fluidtypes, only: component_fluid
      use global,     only: use_fargo

      type(cg_list_element), pointer  :: cgl
      class(component_fluid), pointer :: fl
      integer                         :: i, j, k

      chspeed = huge(1.)
      if (has_ion) then
         if (use_fargo) call die("[hdc:update_chspeed] FARGO is not implemented here yet.")
         if (dom%geometry_type /= GEO_XYZ) call die("[hdc:update_chspeed] non-cartesian geometry not implemented yet.")
         chspeed = small
         fl => flind%ion
         cgl => leaves%first
         do while (associated(cgl))
            do k = cgl%cg%ks, cgl%cg%ke
               do j = cgl%cg%js, cgl%cg%je
                  do i = cgl%cg%is, cgl%cg%ie
                     chspeed = max(chspeed, maxval(abs(cgl%cg%u(fl%imx:fl%imz, i, j, k) / cgl%cg%u(fl%idn, i, j, k)) + fl%get_cs(i, j, k, cgl%cg%u, cgl%cg%b, cgl%cg%cs_iso2)))
                  enddo
               enddo
            enddo
            cgl => cgl%nxt
         enddo
      endif

    return
    
  end subroutine update_chspeed
!--------------------------------------------------------------------------------------------------------------
  subroutine init_psi

     use cg_list,          only: cg_list_element
     use cg_leaves,        only: leaves
     use constants,        only: psi_n
     use global,           only: psi_0
     use named_array_list, only: qna

     type(cg_list_element), pointer :: cgl

     if (qna%exists(psi_n)) then
        cgl => leaves%first
        do while (associated(cgl))
           cgl%cg%q(qna%ind(psi_n))%arr = psi_0
           cgl => cgl%nxt
        enddo
     endif

  end subroutine init_psi
  !---------------------------------------------------------------------------------------------------------------
  subroutine glm_mhd(n, psi_l, psi_r, b_ccl, b_ccr, b_cc, psif)

    ! external procedures
    
    use constants,  only: half, one, xdim

    ! arguments
    
    implicit none
    
    integer,              intent(in)    :: n
    real, dimension(:,:), intent(out) :: psif
    real, dimension(:,:), intent(in) :: psi_l
    real, dimension(:,:), intent(in) :: psi_r
    real, dimension(:,:), intent(out) :: b_cc
    real, dimension(:,:), intent(in) :: b_ccl
    real, dimension(:,:), intent(in) :: b_ccr

    ! local declarations

    integer                           :: i
    real, dimension(xdim)            :: glm_b
    real, dimension(size(psif,1),size(psif,2)) :: glm_psi

    do i = 1, n
       glm_b(xdim) = half*((b_ccl(xdim,i)+b_ccr(xdim,i)) - (one/chspeed)*(psi_r(1,i)-psi_l(1,i)))
       glm_psi(1,i)    = half*((psi_r(1,i)+psi_l(1,i)) - chspeed*(b_ccr(xdim,i)-b_ccl(xdim,i)))
   
       b_cc(xdim,i) = glm_psi(1,i)
       psif(1,i) = chspeed**2.0*glm_b(xdim)
    end do
    
  end subroutine glm_mhd
end module hdc

#endif
