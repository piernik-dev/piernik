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
  public :: chspeed, update_chspeed, init_psi, glm_mhd, glmdamping

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

!>
!! \brief calculate vector of fluxes for B_x and psi
!!
!! \details This follows Dedner et al. (2002), Eq. 42.
!! The result is already simplified: flux of B_x = psi_m and flux of psi = c_h^2 * B_x
!<

  subroutine glm_mhd(psi_l, psi_r, b_ccl, b_ccr, b_cc, psif)

    ! external procedures
    
    use constants,  only: half

    ! arguments
    
    implicit none
    
    real, dimension(:), intent(out) :: psif
    real, dimension(:), intent(in) :: psi_l
    real, dimension(:), intent(in) :: psi_r
    real, dimension(:), intent(out) :: b_cc
    real, dimension(:), intent(in) :: b_ccl
    real, dimension(:), intent(in) :: b_ccr


    b_cc = half * ((psi_r + psi_l) - chspeed * (b_ccr - b_ccl))
    psif = half * chspeed * (chspeed * (b_ccl + b_ccr) - (psi_r - psi_l))
    
  end subroutine glm_mhd

!-----------------------------------------------------------------------------------------------------------

!>
  !! Parabolic damping to psi
!<
  subroutine glmdamping(n, psi)
    
    use global,    only: cfl

    implicit none

    integer,                       intent(in)    :: n
    real, dimension(:),            intent(inout) :: psi

    real                                         :: glm_alpha ! Should be move to globals.F90 and possible over-writing from problem.par enabled.
    integer                                      :: i

    glm_alpha = 0.1

    do i= 1,n
       psi = psi*exp(-glm_alpha*cfl)
    end do
    
  end subroutine glmdamping

end module hdc

#endif /* GLM */
!--------------------------------------------------------------------------------------------
!subroutine glmpdamping(n, cg, dtodx, glm_cp, psi_dtn, psi)
!use constants, only: xdim
!use grid_cont, only: grid_container
!type(grid_container), pointer, intent(in)    :: cg
!real,                          intent(in)    :: dtodx
!real, dimension(:),            intent(out)   :: glm_cp
!real, dimension(:),            intent(out)   :: parb_damping
!glm_cp    = sqrt(cg%dl(xdim)*chspeed/glm_alpha)
!parb_damping   = exp(-dtodx*glm_alpha*chspeed) 
