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

module initfluids

#ifdef IONIZED
  use initionized, only : init_ionized, gamma_ion, cs_iso_ion, cs_iso_ion2
  use fluidindex,  only : i_ion
#endif /* IONIZED */

#ifdef NEUTRAL
  use initneutral, only : init_neutral, gamma_neu, cs_iso_neu, cs_iso_neu2
  use fluidindex,  only : i_neu
#endif /* NEUTRAL */

#ifdef DUST
  use initdust, only : init_dust
#endif /* DUST */

#ifdef COSM_RAYS
  use initcosmicrays, only : init_cosmicrays
#endif /* COSM_RAYS */

  use fluidindex,   only   : nfluid, fluid_index

  implicit none 

  real, allocatable :: gamma(:) 
  real              :: cs_iso, cs_iso2 

  contains

  subroutine init_fluids


#ifdef IONIZED
  call init_ionized
#endif /* IONIZED */

#ifdef NEUTRAL
  call init_neutral
#endif /* NEUTRAL */

#ifdef DUST
  call init_dust
#endif /* DUST */

#ifdef COSM_RAYS
  call init_cosmicrays
#endif /* COSM_RAYS */

  call fluid_index

  allocate(gamma(nfluid)) 
     
#if defined NEUTRAL && defined IONIZED 
    if(cs_iso_neu /= cs_iso_ion) &
        write(*,*) "WARNING: 'cs_iso_neu' and 'cs_iso_ion' should be equal" 
#endif /* defined NEUTRAL && defined IONIZED  */ 
   
#ifdef IONIZED 
    gamma(i_ion) = gamma_ion 
    cs_iso   = cs_iso_ion 
    cs_iso2  = cs_iso_ion2 
#endif /* IONIZED */ 
#ifdef NEUTRAL   
    gamma(i_neu) = gamma_neu 
    cs_iso  = cs_iso_neu  
    cs_iso2 = cs_iso_neu2 
#endif /* NEUTRAL  */  

  end subroutine init_fluids

end module initfluids
