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
!! \brief (MH) Module to initialise all fluids, fluid components, tracers 
!! and other dependent variables which are relevant for the current problem (doxy comments ready).
!!
!! The general purpose of the multi-fluid framework is to simplify all code 
!! modifications when a new fluid, component or variable is added.
!! We want to avoid, for example, touching boundary conditions routines, 
!! when new fluid is added. We formulate boundary conditions in a general manner 
!! and apply them to all densities in one instant, through the array indexes: 
!! iarr_all_dn, iarr_all_mx, etc ... (see comments to fluidindex module).
!!
!! \par DEFINITIONS
!!
!! \n \b Fluid: ingredient characterised by mass density, momenta and optionally 
!! energy density. Examples: ionised fluid, neutral fluid, dust fluid.
!! Variable nfluid (in fluidindex) counts fluids. 
!!
!! \n \b Non-isothermal \b fluid: the fluid engaging energy equation. 
!! Variable "nadiab" (defined in fluidindex) counts independent energy equations 
!! used for fluids description. 
!! \todo Change variable name "nadiab" to "nenerg". Reason: we are not limited 
!!       to isothermal and  adiabatic equos. Energy equation will be used 
!!       for non-adiabatic fluids in presence of cooling and heating. 
!!
!! \n \b Variable: Single quantity, such as gas density, momentum component, 
!!    energy density, CR energy density, etc ...
!!
!! \n \b Component: An ingredient which, contrary to fluids, does not't involve momenta. 
!!               Examples: CR energy density, or a set variables describing 
!!               several CR energy bins, or CR species.
!!
!! \n 
!! \n All these ingredients are organised in the module fluidindex into 
!! the array of conservative variables \a u(iflv,i,j,k), where \a iflv is the 
!! index of fluid variable. 
!!
!! \par The module initfluids is organised as follows:
!! \n (1)  All fluids defined in "piernik.def" are initialised subsequently.
!! \n (2)  The routine fluidindex is invoked to construct arrays of indexes for 
!!         each fluid, e.g.. ionised_index, neutral_index, etc. 
!!         See fluidindex for more details.
!! \n (3)  Physical parameters common for all fluids are computed if necessary.
!!     
!! \todo Change typo in: "ionized_index" to "ionised_index" everywhere in the code.
!! \todo Subdivide different fluids into species
!! \warning check if cs_iso and cs_neu are correctly defined (end of init_fluids 
!!  subroutine) for your purposes (if used).
!<

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

  real, allocatable :: gamma(:)          !< array containing adiabatic indices of all fluids, indexed by   ifluid = i_ion, i_neu, etc.
  real              :: cs_iso            !< isothermal sound speed for a mixture of fluids
  real              :: cs_iso2           !< square of isothermal sound speed for a mixture of fluids  

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
