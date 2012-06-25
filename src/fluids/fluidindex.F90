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
   use constants,  only: ndims, xdim, ydim, zdim
   use fluidtypes, only: var_numbers

   implicit none

   private :: ndims, xdim, ydim, zdim, var_numbers ! QA_WARN: prevent reexporting
   public ! QA_WARN no secrets are kept here

   type(var_numbers),save :: flind     !< COMMENT ME

   integer(kind=4), parameter  :: nmag = ndims     !< number of magnetic field components
   integer(kind=4), parameter  :: idn = 1          !< position of density in the vector of conserv. variables for single fluid
   integer(kind=4), parameter  :: imx = idn + xdim !< position of x-mom. in the vector of conserv. variables for single fluid
   integer(kind=4), parameter  :: imy = idn + ydim !< position of y-mom. in the vector of conserv. variables for single fluid
   integer(kind=4), parameter  :: imz = idn + zdim !< position of z-mom. in the vector of conserv. variables for single fluid
   integer(kind=4), parameter  :: ien = imz + 1    !< position of energy density in the vector of conserv. variables for single fluid (only for adiabatic fluids)
   integer(kind=4), parameter  :: icx = xdim       !< index of x-component of current density
   integer(kind=4), parameter  :: icy = ydim       !< index of y-component of current density
   integer(kind=4), parameter  :: icz = zdim       !< index of z-component of current density

   integer(kind=4), allocatable, dimension(:) :: iarr_all_dn   !< array of indexes pointing to mass densities of all fluids
   integer(kind=4), allocatable, dimension(:) :: iarr_all_sg   !< array of indexes pointing to mass densities of all selfgravitating fluids
   integer(kind=4), allocatable, dimension(:) :: iarr_all_mx   !< array of indexes pointing to mom. densities of all fluids
   integer(kind=4), allocatable, dimension(:) :: iarr_all_my   !< array of indexes pointing to mom. densities of all fluids
   integer(kind=4), allocatable, dimension(:) :: iarr_all_mz   !< array of indexes pointing to mom. densities of all fluids
   integer(kind=4), allocatable, dimension(:) :: iarr_all_en   !< array of indexes pointing to ener. densities of all fluids
   integer(kind=4), allocatable, dimension(:) :: iarr_all_crn   !< array of indexes pointing to ener. densities of all nuclear CR-components
   integer(kind=4), allocatable, dimension(:) :: iarr_all_cre   !< array of indexes pointing to ener. densities of all electron CR-components
   integer(kind=4), allocatable, dimension(:), target :: iarr_all_crs   !< array of indexes pointing to ener. densities of all CR-components
   integer(kind=4), allocatable, dimension(:), target :: iarr_all_trc   !< array of indexes pointing to tracers
   integer(kind=4), allocatable, dimension(:,:) :: iarr_all_swp !< array (size = flind) of all fluid indexes in the order depending on sweeps direction

   integer(kind=4), allocatable, dimension(:) :: iarr_all_mag  !< array (size = nmag) of all magnetic field components
   integer(kind=4), allocatable, dimension(:,:) :: iarr_mag_swp !< array (size = nmag) of all mag. field indexes in the order depending on sweeps direction

   integer(kind=4) :: i_sg                                     !< index denoting position of the selfgravitating fluid in the row of fluids - should be an iarr_sg !

contains

   subroutine set_fluidindex_arrays(fl, have_ener)

      use constants, only: I_ONE
      use fluidtypes, only: component_fluid

      implicit none

      class(component_fluid), intent(inout) :: fl
      logical, intent(in)                  :: have_ener

#ifdef ISO
      logical :: fnord
#endif /* ISO */

      iarr_all_swp(:,fl%beg:fl%end) = fl%iarr_swp(:,:)

      if (fl%is_selfgrav) then
         i_sg = i_sg + I_ONE
         iarr_all_sg(i_sg) = fl%idn
      endif

      iarr_all_dn(fl%pos) = fl%idn
      iarr_all_mx(fl%pos) = fl%imx
      iarr_all_my(fl%pos) = fl%imy
      iarr_all_mz(fl%pos) = fl%imz

#ifdef ISO
      if (.false.) fnord = have_ener ! suppress compiler warning on unused argument
#else /* !ISO */
      if (have_ener) iarr_all_en(fl%pos) = fl%ien
#endif /* !ISO */

   end subroutine set_fluidindex_arrays

!>
!! \brief Subroutine fluid_index constructing all multi-fluid indexes used in other parts
!! of PIERNIK code
!<
   subroutine fluid_index

!      use diagnostics,    only: my_allocate
      use constants,      only: xdim, zdim
      use fluids_pub,     only: has_dst, has_ion, has_neu
#if defined IONIZED || defined COSM_RAYS || defined TRACER
      use constants,      only: ydim
#endif /* IONIZED || COSM_RAYS || TRACER */
      use constants,      only: ndims
      use initionized,    only: ionized_index, ion_fluid
      use initneutral,    only: neutral_index, neutral_fluid
      use initdust,       only: dust_index, dust_fluid

#ifdef COSM_RAYS
      use initcosmicrays, only: iarr_crn, iarr_cre, iarr_crs, cosmicray_index
#endif /* COSM_RAYS */
#ifdef TRACER
      use inittracer,     only: tracer_index, iarr_trc
#endif /* TRACER */

      implicit none

      integer                          :: i

      i_sg        = 0

      if (has_ion) then
         !  Compute indexes for the ionized fluid and update counters
         allocate(ion_fluid::flind%ion)
         call ionized_index(flind)
      endif

      if (has_neu) then
         !  Compute indexes for the neutral fluid and update counters
         allocate(neutral_fluid::flind%neu)
         call neutral_index(flind)
      endif

      if (has_dst) then
      !  Compute indexes for the dust fluid and update counters
         allocate(dust_fluid::flind%dst)
         call dust_index(flind)
      endif

#ifdef COSM_RAYS
!  Compute indexes for the CR component and update counters
      call cosmicray_index(flind)
#endif /* !COSM_RAYS */

#ifdef TRACER
      call tracer_index(flind)
#endif /* TRACER */

! Allocate index arrays
      if (has_ion) allocate(iarr_mag_swp(ndims,nmag),iarr_all_mag(nmag))
      allocate(iarr_all_swp(xdim:zdim, flind%all))
      allocate(iarr_all_dn(flind%fluids),iarr_all_mx(flind%fluids),iarr_all_my(flind%fluids),iarr_all_mz(flind%fluids))
      allocate(iarr_all_sg(flind%fluids_sg))
#ifdef ISO
      allocate(iarr_all_en(0))
#else /* !ISO */
      allocate(iarr_all_en(flind%energ))
#endif /* !ISO */

#ifdef COSM_RAYS
      allocate(iarr_all_crn(flind%crn%all))
      allocate(iarr_all_cre(flind%cre%all))
      allocate(iarr_all_crs(flind%crs%all))
#else /* !COSM_RAYS */
      allocate(iarr_all_crn(0))
      allocate(iarr_all_cre(0))
      allocate(iarr_all_crs(0))
#endif /* !COSM_RAYS */

#ifdef TRACER
      allocate(iarr_all_trc(flind%trc%all))
#else /* !TRACER */
      allocate(iarr_all_trc(0))
#endif /* !TRACER */

      if (has_ion) then
         ! Compute index arrays for magnetic field
         iarr_mag_swp(xdim,:) = [xdim,ydim,zdim]
         iarr_mag_swp(ydim,:) = [ydim,xdim,zdim]
         iarr_mag_swp(zdim,:) = [zdim,ydim,xdim]
         iarr_all_mag         = [xdim,ydim,zdim]

         ! Compute index arrays for the ionized fluid
         call set_fluidindex_arrays(flind%ion,.true.)
      endif

      ! Compute index arrays for the neutral fluid
      if (has_neu) call set_fluidindex_arrays(flind%neu,.true.)

      ! Compute index arrays for the dust fluid
      if (has_dst) call set_fluidindex_arrays(flind%dst,.false.)

#ifdef COSM_RAYS
! Compute index arrays for the CR components
      iarr_all_swp(xdim,flind%crs%beg:flind%crs%end) = iarr_crs
      iarr_all_swp(ydim,flind%crs%beg:flind%crs%end) = iarr_crs
      iarr_all_swp(zdim,flind%crs%beg:flind%crs%end) = iarr_crs

      iarr_all_crn(1:flind%crn%all) = iarr_crn
      iarr_all_cre(1:flind%cre%all) = iarr_cre
      iarr_all_crs(1:flind%crs%all) = iarr_crs
#endif /* COSM_RAYS */

#ifdef TRACER
      iarr_all_swp(xdim,flind%trc%beg:flind%trc%end) = iarr_trc
      iarr_all_swp(ydim,flind%trc%beg:flind%trc%end) = iarr_trc
      iarr_all_swp(zdim,flind%trc%beg:flind%trc%end) = iarr_trc

      iarr_all_trc(1:flind%trc%all) = iarr_trc
#endif /* TRACER */

      allocate(flind%all_fluids(flind%fluids))

      i = 1
      if (has_ion) then
         flind%all_fluids(i)%fl => flind%ion
         i = i + 1
      endif

      if (has_neu) then
         flind%all_fluids(i)%fl => flind%neu
         i = i + 1
      endif

      if (has_dst) then
!         allocate(dust_fluid::flind%all_fluids(i)%fl) ! = flind%dst
         flind%all_fluids(i)%fl => flind%dst
         i = i + 1
      endif
   end subroutine fluid_index

   subroutine cleanup_fluidindex

      use diagnostics, only: my_deallocate
      use fluids_pub,  only: has_ion, has_neu, has_dst

      implicit none

      integer :: i

      call my_deallocate(iarr_mag_swp)
      call my_deallocate(iarr_all_swp)
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

      call my_deallocate(iarr_all_trc)

      do i = lbound(flind%all_fluids, dim=1), ubound(flind%all_fluids, dim=1)
         deallocate(flind%all_fluids(i)%fl%iarr)
         deallocate(flind%all_fluids(i)%fl%iarr_swp)
      enddo
      deallocate(flind%all_fluids)
      if (has_ion) deallocate(flind%ion) ! cannot check if allocated, because it is not allocatable
      if (has_neu) deallocate(flind%neu)
      if (has_dst) deallocate(flind%dst)

   end subroutine cleanup_fluidindex

end module fluidindex
