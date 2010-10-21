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
#include "piernik.def"
!>
!! \brief (KK)
!<
module fluidtypes
   type :: component
      integer :: all = 0              !< number of all variables in fluid/component
      integer :: beg = 0              !< beginning number of variables in fluid/component
      integer :: end = 0              !< end number of variables in fluid/component
      integer :: pos = 0              !< index denoting position of the fluid in the row of fluids
   end type component

   type :: component_fluid
      integer :: all = 0              !< number of all variables in fluid/component
      integer :: beg = 0              !< beginning number of variables in fluid/component
      integer :: end = 0              !< end number of variables in fluid/component
      integer :: pos = 0              !< index denoting position of the fluid in the row of fluids
      integer :: idn = -1
      integer :: imx = -1
      integer :: imy = -1
      integer :: imz = -1
      integer :: ien = -1
   end type component_fluid

   type :: var_numbers
      integer :: all         = 0     !< total number of fluid variables = the size of array \a u(:,:,;,:) in the first index
      integer :: fluids      = 0     !< number of fluids (ionized gas, neutral gas, dust)
      integer :: adiab       = 0     !< number of adiabatic fluids (indicating the presence of energy density in the vector of conservative variables)
      integer :: components  = 0     !< number of components, such as CRs, tracers, magnetic helicity (in future), whose formal description does not involve [???]
      integer :: fluids_sg   = 0     !< number of selfgravitating fluids (ionized gas, neutral gas, dust)

#ifdef IONIZED
      type(component_fluid) :: ion         !< numbers of variables for the ionized fluid
      integer, allocatable, dimension(:)  :: iarr_ion
      integer, allocatable, dimension(:)  :: iarr_ion_swpx, iarr_ion_swpy, iarr_ion_swpz
#endif /* IONIZED */
#ifdef NEUTRAL
      type(component_fluid) :: neu         !< numbers of variables for the neutral fluid
      integer, allocatable, dimension(:)  :: iarr_neu
      integer, allocatable, dimension(:)  :: iarr_neu_swpx, iarr_neu_swpy, iarr_neu_swpz
#endif /* NEUTRAL */
#ifdef DUST
      type(component_fluid) :: dst         !< numbers of variables for the dust fluid
      integer, allocatable, dimension(:)  :: iarr_dst
      integer, allocatable, dimension(:)  :: iarr_dst_swpx, iarr_dst_swpy, iarr_dst_swpz
#endif /* DUST */
#ifdef COSM_RAYS
      type(component) :: crs         !< numbers of variables in all cosmic ray components
      type(component) :: crn         !< numbers of variables in cosmic ray nuclear components
      type(component) :: cre         !< numbers of variables in cosmic ray electron components
#endif /* COSM_RAYS */
   end type var_numbers

end module fluidtypes
