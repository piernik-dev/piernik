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
#define RNG 2:n-1

!>
!! \brief Computation of Cosmic Ray sources
!!
!!
!<

module sourcecosmicrays
! pulled by COSM_RAYS_SOURCES
   implicit none

   private
   public :: src_crn

contains

!==========================================================================================

   subroutine src_crn(uu, n, decrn)

      use cr_data,       only: eCRSP, cr_table, cr_tau, cr_sigma, icr_C12, icr_Be9, icr_Be10, icr_Li7, icr_N14, icr_O16
      use fluids_pub,    only: has_ion, has_neu
      use fluidindex,    only: flind

      implicit none

      integer(kind=4),                  intent(in)  :: n
      real, dimension(flind%all,n),     intent(in)  :: uu
      real, dimension(flind%crn%all,n), intent(out) :: decrn

! locals
      real, allocatable       :: dgas(:)
      real,         parameter :: gamma_lor = 10.
      real(kind=8), parameter :: speed_of_light = 3e10*1e6*365.*24.*60.*60. !< cm/Myr \deprecated BEWARE: this line breaks unit consistency, move it to units.F90 and use scaling
      real,         parameter :: ndim = 2.

      allocate(dgas(n))

      dgas(:) = 0.0

      if (has_ion) dgas(:) = dgas(:) + uu(flind%ion%idn,:)
      if (has_neu) dgas(:) = dgas(:) + uu(flind%neu%idn,:)

      decrn(:,:) = 0.0

      if (eCRSP(icr_Be10)) decrn(cr_table(icr_Be10),:) = decrn(cr_table(icr_Be10),:) - 1./ndim*uu(flind%crn%beg-1+cr_table(icr_Be10),:)/gamma_lor/cr_tau(cr_table(icr_Be10))

      if (eCRSP(icr_C12)) then
         if (eCRSP(icr_Li7)) then
            decrn(cr_table(icr_C12 ),:) = decrn(cr_table(icr_C12 ),:) - 1./ndim*dgas(:)*speed_of_light*cr_sigma(cr_table(icr_C12),cr_table(icr_Li7 ))*uu(flind%crn%beg-1+cr_table(icr_C12),:)
            decrn(cr_table(icr_Li7 ),:) = decrn(cr_table(icr_Li7 ),:) + 1./ndim*dgas(:)*speed_of_light*cr_sigma(cr_table(icr_C12),cr_table(icr_Li7 ))*uu(flind%crn%beg-1+cr_table(icr_C12),:)
         endif
         if (eCRSP(icr_Be9)) then
            decrn(cr_table(icr_C12 ),:) = decrn(cr_table(icr_C12 ),:) - 1./ndim*dgas(:)*speed_of_light*cr_sigma(cr_table(icr_C12),cr_table(icr_Be9 ))*uu(flind%crn%beg-1+cr_table(icr_C12),:)
            decrn(cr_table(icr_Be9 ),:) = decrn(cr_table(icr_Be9 ),:) + 1./ndim*dgas(:)*speed_of_light*cr_sigma(cr_table(icr_C12),cr_table(icr_Be9 ))*uu(flind%crn%beg-1+cr_table(icr_C12),:)
         endif
         if (eCRSP(icr_Be10)) then
            decrn(cr_table(icr_C12 ),:) = decrn(cr_table(icr_C12 ),:) - 1./ndim*dgas(:)*speed_of_light*cr_sigma(cr_table(icr_C12),cr_table(icr_Be10))*uu(flind%crn%beg-1+cr_table(icr_C12),:)
            decrn(cr_table(icr_Be10),:) = decrn(cr_table(icr_Be10),:) + 1./ndim*dgas(:)*speed_of_light*cr_sigma(cr_table(icr_C12),cr_table(icr_Be10))*uu(flind%crn%beg-1+cr_table(icr_C12),:)
         endif
      endif

      if (eCRSP(icr_N14)) then
         if (eCRSP(icr_Li7)) then
            decrn(cr_table(icr_N14 ),:) = decrn(cr_table(icr_N14 ),:) - 1./ndim*dgas(:)*speed_of_light*cr_sigma(cr_table(icr_N14),cr_table(icr_Li7 ))*uu(flind%crn%beg-1+cr_table(icr_N14),:)
            decrn(cr_table(icr_Li7 ),:) = decrn(cr_table(icr_Li7 ),:) + 1./ndim*dgas(:)*speed_of_light*cr_sigma(cr_table(icr_N14),cr_table(icr_Li7 ))*uu(flind%crn%beg-1+cr_table(icr_N14),:)
         endif
         if (eCRSP(icr_Be9)) then
            decrn(cr_table(icr_N14 ),:) = decrn(cr_table(icr_N14 ),:) - 1./ndim*dgas(:)*speed_of_light*cr_sigma(cr_table(icr_N14),cr_table(icr_Be9 ))*uu(flind%crn%beg-1+cr_table(icr_N14),:)
            decrn(cr_table(icr_Be9 ),:) = decrn(cr_table(icr_Be9 ),:) + 1./ndim*dgas(:)*speed_of_light*cr_sigma(cr_table(icr_N14),cr_table(icr_Be9 ))*uu(flind%crn%beg-1+cr_table(icr_N14),:)
         endif
         if (eCRSP(icr_Be10)) then
            decrn(cr_table(icr_N14 ),:) = decrn(cr_table(icr_N14 ),:) - 1./ndim*dgas(:)*speed_of_light*cr_sigma(cr_table(icr_N14),cr_table(icr_Be10))*uu(flind%crn%beg-1+cr_table(icr_N14),:)
            decrn(cr_table(icr_Be10),:) = decrn(cr_table(icr_Be10),:) + 1./ndim*dgas(:)*speed_of_light*cr_sigma(cr_table(icr_N14),cr_table(icr_Be10))*uu(flind%crn%beg-1+cr_table(icr_N14),:)
         endif
      endif

      if (eCRSP(icr_O16)) then
         if (eCRSP(icr_Li7)) then
            decrn(cr_table(icr_O16 ),:) = decrn(cr_table(icr_O16 ),:) - 1./ndim*dgas(:)*speed_of_light*cr_sigma(cr_table(icr_O16),cr_table(icr_Li7 ))*uu(flind%crn%beg-1+cr_table(icr_O16),:)
            decrn(cr_table(icr_Li7 ),:) = decrn(cr_table(icr_Li7 ),:) + 1./ndim*dgas(:)*speed_of_light*cr_sigma(cr_table(icr_O16),cr_table(icr_Li7 ))*uu(flind%crn%beg-1+cr_table(icr_O16),:)
         endif
         if (eCRSP(icr_Be9)) then
            decrn(cr_table(icr_O16 ),:) = decrn(cr_table(icr_O16 ),:) - 1./ndim*dgas(:)*speed_of_light*cr_sigma(cr_table(icr_O16),cr_table(icr_Be9 ))*uu(flind%crn%beg-1+cr_table(icr_O16),:)
            decrn(cr_table(icr_Be9 ),:) = decrn(cr_table(icr_Be9 ),:) + 1./ndim*dgas(:)*speed_of_light*cr_sigma(cr_table(icr_O16),cr_table(icr_Be9 ))*uu(flind%crn%beg-1+cr_table(icr_O16),:)
         endif
         if (eCRSP(icr_Be10)) then
            decrn(cr_table(icr_O16 ),:) = decrn(cr_table(icr_O16 ),:) - 1./ndim*dgas(:)*speed_of_light*cr_sigma(cr_table(icr_O16),cr_table(icr_Be10))*uu(flind%crn%beg-1+cr_table(icr_O16),:)
            decrn(cr_table(icr_Be10),:) = decrn(cr_table(icr_Be10),:) + 1./ndim*dgas(:)*speed_of_light*cr_sigma(cr_table(icr_O16),cr_table(icr_Be10))*uu(flind%crn%beg-1+cr_table(icr_O16),:)
         endif
      endif




   end subroutine src_crn

end module sourcecosmicrays
