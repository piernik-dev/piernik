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
!>
!! \brief (MH)
!<
module crcomposition
! pulled by COSM_RAYS_SOURCES
   implicit none

   public               ! QA_WARN no secrets are kept here

! this type looks useful but is unused.
!!$   integer, parameter :: isoname_len = 8
!!$   type cr_component
!!$      character(len=isoname_len) :: isotope     !< isotope name, eg. Be10
!!$      integer          :: index=      !< relative index (with respect to crn_beg)
!!$      real             :: abund=      !< initial abundance relative to H
!!$   end type cr_component

   integer, parameter :: icr_H1   = 1
   integer, parameter :: icr_C12  = 2
   integer, parameter :: icr_Be9  = 3
   integer, parameter :: icr_Be10 = 4
   integer, parameter :: icr_N14  = 5         !< from decay of Be7 with tau of 0.3 years
   integer, parameter :: icr_O16  = 6
   integer, parameter :: icr_Li7  = 7    !!! BEWARE: ncrn should be set up gerater than maximum isotope numeber,
                                  !!! which should be smaller than ncr_max (<10 currently)

!<====Cross sections for spallation from Garcia-Munoz 1987 (see also Longair)====>

   real, parameter :: mbarn=1e-27 !cm2   !!! BEWARE: this line breaks unit consistency, move it to constants.F90 and use scaling

   real, parameter :: sigma_C12_Li7  = 10*mbarn
   real, parameter :: sigma_C12_Be9  = 6*mbarn
   real, parameter :: sigma_C12_Be10 = 3.5*mbarn

   real, parameter :: sigma_N14_Li7  = 9.5*mbarn

   real, parameter :: sigma_O16_Li7  = 9.5*mbarn
   real, parameter :: sigma_O16_Be9  = 4.5*mbarn
   real, parameter :: sigma_O16_Be10 = 2*mbarn

!<====Decay half live times from Garcia-Munoz 1987====>

   real,parameter :: Myear=1d6*365*24*60*60 !s !!! BEWARE: this line breaks unit consistency, move it to constants.F90 and use scaling

   real, parameter   :: tau_Be10 = 1.6 !Myr !!! BEWARE: this line breaks unit consistency, move it to constants.F90 and use scaling

!<Initial source abundances (in numer density) relative to hydrogen (compare e.g. Longair)>

   real, parameter :: primary_C12  =  4.5e-3
   real, parameter :: primary_N14  =  1e-3
   real, parameter :: primary_O16  =  4e-3

end module crcomposition
