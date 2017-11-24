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
! pulled by MAGNETIC && RIEMANN
  implicit none

  real, protected :: chspeed

  private
  public :: chspeed, update_chspeed, init_psi

contains

  subroutine update_chspeed() 

    use fluidindex, only: flind
    use fluids_pub, only: has_ion
    use dataio_pub, only: warn

    if(has_ion) then 
       chspeed = flind%ion%snap%cs_max%val
    else
       call warn('Kenobi is in exile!')
       ! set up some safe value or die
    end if
    
    return
    
  end subroutine update_chspeed
!--------------------------------------------------------------------------------------------------------------
  subroutine init_psi

     use cg_list,          only: cg_list_element
     use cg_leaves,        only: leaves
     use constants,        only: phi_n
     use global,           only: psi_0
     use named_array_list, only: qna

     type(cg_list_element), pointer :: cgl

     if (qna%exists(phi_n)) then
        cgl => leaves%first
        do while (associated(cgl))
           cgl%cg%q(qna%ind(phi_n))%arr = psi_0
           cgl => cgl%nxt
        enddo
     endif

  end subroutine init_psi
!---------------------------------------------------------------------------------------------------------------
end module hdc

#endif
