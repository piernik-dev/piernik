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
!! \brief (MH) [R] In this module fluid variables of individual fluids are indexed to make use of the single array
!! \a u(:,:,:,:) containing all fluid variables.
!!
!! \par The purpose of this module is to compute:
!!
!! \n (1) positions of all fluid variables referenced through the first index of array u(:,:,:,:),
!! \n (2) %arrays of indexes to make easy reference to all gas densities, x,y,z-components of momenta,
!! energy densities, CR energy densities, transposed components of x,y,z-momenta in directional %sweeps, etc ...
!!
!! \n\b Note: Basic definitions and more information in initfluids module.
!!
!! \todo All stuff related specifically to magnetic field should be embraced
!! (everywhere in the code) with the precompiler definition "MAGNETIC" instead of "IONIZED"
!<

module fluidindex
! pulled by ANY
   use types, only: var_numbers

   implicit none

   public ! QA_WARN no secrets are kept here

   type(var_numbers),save :: nvar     !< COMMENT ME

   integer, parameter  :: nmag = 3    !< number of magnetic field components

   integer, parameter  :: ibx = 1     !< index of x-component of magnetic field
   integer, parameter  :: iby = 2     !< index of y-component of magnetic field
   integer, parameter  :: ibz = 3     !< index of z-component of magnetic field
   integer, parameter  :: idn = 1     !< position of density in the vector of conserv. variables for single fluid
   integer, parameter  :: imx = 2     !< position of x-mom. in the vector of conserv. variables for single fluid
   integer, parameter  :: imy = 3     !< position of y-mom. in the vector of conserv. variables for single fluid
   integer, parameter  :: imz = 4     !< position of z-mom. in the vector of conserv. variables for single fluid
   integer, parameter  :: ien = 5     !< position of energy density in the vector of conserv. variables for single fluid (only for adiabatic fluids)

   integer, parameter  :: icx = 1     !< index of x-component of current density
   integer, parameter  :: icy = 2     !< index of y-component of current density
   integer, parameter  :: icz = 3     !< index of z-component of current density

   integer, allocatable, dimension(:) :: iarr_all_dn   !< array of indexes pointing to mass densities of all fluids
   integer, allocatable, dimension(:) :: iarr_all_sg   !< array of indexes pointing to mass densities of all selfgravitating fluids
   integer, allocatable, dimension(:) :: iarr_all_mx   !< array of indexes pointing to mom. densities of all fluids
   integer, allocatable, dimension(:) :: iarr_all_my   !< array of indexes pointing to mom. densities of all fluids
   integer, allocatable, dimension(:) :: iarr_all_mz   !< array of indexes pointing to mom. densities of all fluids
   integer, allocatable, dimension(:) :: iarr_all_en   !< array of indexes pointing to ener. densities of all fluids
   integer, allocatable, dimension(:) :: iarr_all_crn   !< array of indexes pointing to ener. densities of all nuclear CR-components
   integer, allocatable, dimension(:) :: iarr_all_cre   !< array of indexes pointing to ener. densities of all electron CR-components
   integer, allocatable, dimension(:), target :: iarr_all_crs   !< array of indexes pointing to ener. densities of all CR-components
   integer, allocatable, dimension(:) :: iarr_all_swpx !< array (size = nvar) of all fluid indexes in the original order
   integer, allocatable, dimension(:) :: iarr_all_swpy !< array (size = nvar) of all fluid indexes with \a x and \a y components of mom. interchanged
   integer, allocatable, dimension(:) :: iarr_all_swpz !< array (size = nvar) of all fluid indexes with \a x and \a z components of mom. interchanged

   integer, allocatable, dimension(:) :: iarr_all_mag  !< array (size = nmag) of all magnetic field components
   integer, allocatable, dimension(:) :: iarr_mag_swpx !< array (size = nmag) of all mag. field indexes in the original order (same as iarr_all_mag)
   integer, allocatable, dimension(:) :: iarr_mag_swpy !< array (size = nmag) of all mag. field indexes \a x and \a z components interchanged
   integer, allocatable, dimension(:) :: iarr_mag_swpz !< array (size = nmag) of all mag. field indexes \a x and \a z components interchanged

   integer :: i_sg                                     !< index denoting position of the selfgravitating fluid in the row of fluids - should be an iarr_sg !

contains

   subroutine set_fluidindex_arrays(fl, have_ener)

      use types, only: component_fluid

      implicit none

      type(component_fluid), intent(inout) :: fl
      logical, intent(in)                  :: have_ener

#ifdef ISO
      logical :: fnord
#endif /* ISO */

      iarr_all_swpx(fl%beg:fl%end) = fl%iarr_swpx
      iarr_all_swpy(fl%beg:fl%end) = fl%iarr_swpy
      iarr_all_swpz(fl%beg:fl%end) = fl%iarr_swpz

      if (fl%is_selfgrav) then
         i_sg = i_sg + 1
         iarr_all_sg(i_sg) = fl%idn
      endif

      iarr_all_dn(fl%pos) = fl%idn
      iarr_all_mx(fl%pos) = fl%imx
      iarr_all_my(fl%pos) = fl%imy
      iarr_all_mz(fl%pos) = fl%imz

#ifdef ISO
      if (.false.) fnord = have_ener ! suppress compiler warning on unused argument
#else
      if (have_ener) iarr_all_en(fl%pos) = fl%ien
#endif /* ISO */

   end subroutine set_fluidindex_arrays

!>
!! \brief Subroutine fluid_index constructing all multi-fluid indexes used in other parts
!! of PIERNIK code
!<
    subroutine fluid_index

!      use diagnostics,    only: my_allocate
#ifdef IONIZED
      use initionized,    only: ionized_index
#endif /* IONIZED */
#ifdef NEUTRAL
      use initneutral,    only: neutral_index
#endif /* NEUTRAL */
#ifdef DUST
      use initdust,       only: dust_index
#endif /* DUST */

#ifdef COSM_RAYS
      use initcosmicrays, only: iarr_crn, iarr_cre, iarr_crs, cosmicray_index
#endif /* COSM_RAYS */

      implicit none

      integer                          :: i

      i_sg        = 0

#ifdef IONIZED
!  Compute indexes for the ionized fluid and update counters
      allocate(nvar%ion) ! BEWARE: not deallocated
      call ionized_index(nvar)
#endif /* IONIZED */

#ifdef NEUTRAL
!  Compute indexes for the neutral fluid and update counters
      allocate(nvar%neu)
      call neutral_index(nvar)
#endif /* NEUTRAL */

#ifdef DUST
!  Compute indexes for the dust fluid and update counters
      allocate(nvar%dst)
      call dust_index(nvar)
#endif /* DUST */

#ifdef COSM_RAYS
!  Compute indexes for the CR component and update counters
      call cosmicray_index(nvar)
#endif /* !COSM_RAYS */

! Allocate index arrays
#ifdef IONIZED
      allocate(iarr_mag_swpx(nmag),iarr_mag_swpy(nmag),iarr_mag_swpz(nmag),iarr_all_mag(nmag))
#endif /* IONIZED */
      allocate(iarr_all_swpx(nvar%all),iarr_all_swpy(nvar%all),iarr_all_swpz(nvar%all))
      allocate(iarr_all_dn(nvar%fluids),iarr_all_mx(nvar%fluids),iarr_all_my(nvar%fluids),iarr_all_mz(nvar%fluids))
      allocate(iarr_all_sg(nvar%fluids_sg))
#ifdef ISO
      allocate(iarr_all_en(0))
#else /* !ISO */
      allocate(iarr_all_en(nvar%adiab))
#endif /* !ISO */

#ifdef COSM_RAYS
      allocate(iarr_all_crn(nvar%crn%all))
      allocate(iarr_all_cre(nvar%cre%all))
      allocate(iarr_all_crs(nvar%crs%all))
#else /* !COSM_RAYS */
      allocate(iarr_all_crn(0))
      allocate(iarr_all_cre(0))
      allocate(iarr_all_crs(0))
#endif /* !COSM_RAYS */

#ifdef IONIZED
! Compute index arrays for magnetic field
      iarr_mag_swpx = [ibx,iby,ibz]
      iarr_mag_swpy = [iby,ibx,ibz]
      iarr_mag_swpz = [ibz,iby,ibx]
      iarr_all_mag  = [ibx,iby,ibz]
#endif /* IONIZED */

#ifdef IONIZED
! Compute index arrays for the ionized fluid
      call set_fluidindex_arrays(nvar%ion,.true.)
#endif /* IONIZED */

#ifdef NEUTRAL
! Compute index arrays for the neutral fluid
      call set_fluidindex_arrays(nvar%neu,.true.)
#endif /* NEUTRAL */

#ifdef DUST
! Compute index arrays for the dust fluid
      call set_fluidindex_arrays(nvar%dst,.false.)
#endif /* DUST */

#ifdef COSM_RAYS
! Compute index arrays for the CR components
      iarr_all_swpx(nvar%crs%beg:nvar%crs%end) = iarr_crs
      iarr_all_swpy(nvar%crs%beg:nvar%crs%end) = iarr_crs
      iarr_all_swpz(nvar%crs%beg:nvar%crs%end) = iarr_crs

      iarr_all_crn(1:nvar%crn%all) = iarr_crn
      iarr_all_cre(1:nvar%cre%all) = iarr_cre
      iarr_all_crs(1:nvar%crs%all) = iarr_crs
#endif /* COSM_RAYS */

      allocate(nvar%all_fluids(nvar%fluids)) ! BEWARE: not deallocated
      i = 1
#ifdef IONIZED
      nvar%all_fluids(i) = nvar%ion ; i = i + 1
#endif /* IONIZED */
#ifdef NEUTRAL
      nvar%all_fluids(i) = nvar%neu ; i = i + 1
#endif /* NEUTRAL */
#ifdef DUST
      nvar%all_fluids(i) = nvar%dst ; i = i + 1
#endif /* DUST */
   end subroutine fluid_index

   subroutine cleanup_fluidindex

      use diagnostics, only: my_deallocate

      implicit none

      call my_deallocate(iarr_mag_swpx)
      call my_deallocate(iarr_mag_swpy)
      call my_deallocate(iarr_mag_swpz)
      call my_deallocate(iarr_all_swpx)
      call my_deallocate(iarr_all_swpy)
      call my_deallocate(iarr_all_swpz)
      call my_deallocate(iarr_all_mag)
      call my_deallocate(iarr_all_dn)
      call my_deallocate(iarr_all_mx)
      call my_deallocate(iarr_all_my)
      call my_deallocate(iarr_all_mz)
      call my_deallocate(iarr_all_sg)
      call my_deallocate(iarr_all_en)

      call my_deallocate(iarr_all_crn)
      call my_deallocate(iarr_all_cre)
      call my_deallocate(iarr_all_crs)

   end subroutine cleanup_fluidindex

end module fluidindex
