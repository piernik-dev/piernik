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
!! \brief (MH) In this module fluid variables of individual fluids are indexed to make use of the single array
!! \a u(:,:,:,:) containing all fluid variables. (doxy comments ready)
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

   use fluidtypes, only: var_numbers
   use types,      only: indx
   implicit none

   type(var_numbers),save :: nvar

   type(indx),save     :: ind         !< derived type variable storing all fluid-, component- and magnetic-field-indexes

   integer, parameter  :: nmag = 3    !< number of magnetic field components

   integer, parameter  :: ibx = 1     !< index of x-component of magnetic field
   integer, parameter  :: iby = 2     !< index of y-component of magnetic field
   integer, parameter  :: ibz = 3     !< index of z-component of magnetic field
   integer, parameter  :: idn = 1     !< position of density in the vector of conserv. variables for single fluid
   integer, parameter  :: imx = 2     !< position of x-mom. in the vector of conserv. variables for single fluid
   integer, parameter  :: imy = 3     !< position of y-mom. in the vector of conserv. variables for single fluid
   integer, parameter  :: imz = 4     !< position of z-mom. in the vector of conserv. variables for single fluid
   integer, parameter  :: ien = 5     !< position of energy density in the vector of conserv. variables for single fluid (only for adiabatic fluids)

#ifdef IONIZED
#ifdef RESISTIVE
   integer, parameter  :: icx = 1     !< index of x-component of current density
   integer, parameter  :: icy = 2     !< index of y-component of current density
   integer, parameter  :: icz = 3     !< index of z-component of current density
#endif /* RESISTIVE */
#endif /* IONIZED */

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

#ifdef IONIZED
   integer, allocatable, dimension(:) :: iarr_all_mag  !< array (size = nmag) of all magnetic field components
   integer, allocatable, dimension(:) :: iarr_mag_swpx !< array (size = nmag) of all mag. field indexes in the original order (same as iarr_all_mag)
   integer, allocatable, dimension(:) :: iarr_mag_swpy !< array (size = nmag) of all mag. field indexes \a x and \a z components interchanged
   integer, allocatable, dimension(:) :: iarr_mag_swpz !< array (size = nmag) of all mag. field indexes \a x and \a z components interchanged
#endif /* IONIZED */

   integer :: i_sg                                     !< index denoting position of the selfgravitating fluid in the row of fluids - should be an iarr_sg !

  contains

!>
!! \brief Subroutine fluid_index constructing all multi-fluid indexes used in other parts
!! of PIERNIK code
!! \param none - all arguments are global variables
!<
    subroutine fluid_index

#ifdef IONIZED
      use initionized,    only: idni, imxi, imyi, imzi, ionized_index, selfgrav_ion
      use initionized,    only: iarr_ion_swpx, iarr_ion_swpy, iarr_ion_swpz
#ifndef ISO
      use initionized,    only: ieni
#endif /* !ISO */
#endif /* IONIZED */

#ifdef NEUTRAL
      use initneutral,    only: idnn, imxn, imyn, imzn, neutral_index, selfgrav_neu
      use initneutral,    only: iarr_neu_swpx, iarr_neu_swpy, iarr_neu_swpz
#ifndef ISO
      use initneutral,    only: ienn
#endif /* !ISO */
#endif /* NEUTRAL */

#ifdef DUST
      use initdust,       only: idnd, imxd, imyd, imzd, dust_index, selfgrav_dst
      use initdust,       only: iarr_dst_swpx, iarr_dst_swpy, iarr_dst_swpz
#endif /* DUST */

#ifdef COSM_RAYS
      use initcosmicrays, only: iarr_crn, iarr_cre, iarr_crs, cosmicray_index
#endif /* COSM_RAYS */


      implicit none

      i_sg        = 0

#ifdef IONIZED
!  Compute indexes for the ionized fluid and update counters
      nvar%ion%beg    = nvar%all + 1
      call ionized_index(nvar%all,nvar%ion%all)
      nvar%ion%end    = nvar%all
      nvar%components = nvar%components + 1
      nvar%fluids     = nvar%fluids + 1
      nvar%ion%pos    = nvar%components
      if (selfgrav_ion)  nvar%fluids_sg = nvar%fluids_sg + 1
#ifndef ISO
      nvar%adiab = nvar%adiab + 1
#endif /* !ISO */
#endif /* IONIZED */

#ifdef NEUTRAL
!  Compute indexes for the neutral fluid and update counters
      nvar%neu%beg    = nvar%all + 1
      call neutral_index(nvar%all,nvar%neu%all)
      nvar%neu%end    = nvar%all
      nvar%components = nvar%components + 1
      nvar%fluids     = nvar%fluids + 1
      nvar%neu%pos    = nvar%components
      if (selfgrav_neu)  nvar%fluids_sg = nvar%fluids_sg + 1
#ifndef ISO
      nvar%adiab = nvar%adiab + 1
#endif /* !ISO */
#endif /* NEUTRAL */

#ifdef DUST
!  Compute indexes for the dust fluid and update counters
      nvar%dst%beg    = nvar%all + 1
      call dust_index(nvar%all,nvar%dst%all)
      nvar%dst%end    = nvar%all
      nvar%components = nvar%components + 1
      nvar%fluids     = nvar%fluids + 1
      nvar%dst%pos    = nvar%components
      if (selfgrav_dst)  nvar%fluids_sg = nvar%fluids_sg + 1
#endif /* DUST */

#ifdef COSM_RAYS
!  Compute indexes for the CR component and update counters
      nvar%crn%beg    = nvar%all + 1
      nvar%crs%beg    = nvar%crn%beg
      call cosmicray_index(nvar%all, nvar%crn%all, nvar%cre%all, nvar%crs%all)
      nvar%crn%end    = nvar%crn%beg + nvar%crn%all - 1
      nvar%cre%beg    = nvar%crn%end + 1
      nvar%cre%end    = nvar%all
      nvar%crs%end    = nvar%cre%end
      if (nvar%crn%all  /= 0) nvar%components = nvar%components + 1
      nvar%crn%pos = nvar%components
      if (nvar%cre%all  /= 0) nvar%components = nvar%components + 1
      nvar%cre%pos = nvar%components
#endif /* COSM_RAYS */

! Allocate index arrays
#ifdef IONIZED
      allocate(iarr_mag_swpx(nmag),iarr_mag_swpy(nmag),iarr_mag_swpz(nmag),iarr_all_mag(nmag))
#endif /* IONIZED */
      allocate(iarr_all_swpx(nvar%all),iarr_all_swpy(nvar%all),iarr_all_swpz(nvar%all))
      allocate(iarr_all_dn(nvar%fluids),iarr_all_mx(nvar%fluids),iarr_all_my(nvar%fluids),iarr_all_mz(nvar%fluids))
      allocate(iarr_all_sg(nvar%fluids_sg))
#ifndef ISO
      allocate(iarr_all_en(nvar%adiab))
#else /* !ISO */
      allocate(iarr_all_en(0))
#endif /* !ISO */

#ifdef COSM_RAYS
      allocate(iarr_all_crn(nvar%crn%all))
      allocate(iarr_all_cre(nvar%cre%all))
      allocate(iarr_all_crs(nvar%crs%all))
#else /* COSM_RAYS */
      allocate(iarr_all_crn(0))
      allocate(iarr_all_cre(0))
      allocate(iarr_all_crs(0))
#endif /* COSM_RAYS */

#ifdef IONIZED
! Compute index arrays for magnetic field
      iarr_mag_swpx = [ibx,iby,ibz]
      iarr_mag_swpy = [iby,ibx,ibz]
      iarr_mag_swpz = [ibz,iby,ibx]
      iarr_all_mag  = [ibx,iby,ibz]
      ind%bx = ibx; ind%by = iby; ind%bz = ibz
#endif /* IONIZED */

#ifdef IONIZED
! Compute index arrays for the ionized fluid
      iarr_all_swpx(nvar%ion%beg:nvar%ion%end) = iarr_ion_swpx
      iarr_all_swpy(nvar%ion%beg:nvar%ion%end) = iarr_ion_swpy
      iarr_all_swpz(nvar%ion%beg:nvar%ion%end) = iarr_ion_swpz

      if (selfgrav_ion) then
         i_sg = i_sg + 1
         iarr_all_sg(i_sg) = idni
      endif
      iarr_all_dn(nvar%ion%pos)      = idni ; ind%dni = idni
      iarr_all_mx(nvar%ion%pos)      = imxi ; ind%mxi = imxi
      iarr_all_my(nvar%ion%pos)      = imyi ; ind%myi = imyi
      iarr_all_mz(nvar%ion%pos)      = imzi ; ind%mzi = imzi
#ifndef ISO
      iarr_all_en(nvar%ion%pos)      = ieni ; ind%eni = ieni
#endif /* !ISO */
#endif /* IONIZED */

#ifdef NEUTRAL
! Compute index arrays for the neutral fluid
      iarr_all_swpx(nvar%neu%beg:nvar%neu%end) = iarr_neu_swpx
      iarr_all_swpy(nvar%neu%beg:nvar%neu%end) = iarr_neu_swpy
      iarr_all_swpz(nvar%neu%beg:nvar%neu%end) = iarr_neu_swpz

      if (selfgrav_neu) then
         i_sg = i_sg + 1
         iarr_all_sg(i_sg) = idnn
      endif
      iarr_all_dn(nvar%neu%pos)      = idnn ; ind%dnn = idnn
      iarr_all_mx(nvar%neu%pos)      = imxn ; ind%mxn = imxn
      iarr_all_my(nvar%neu%pos)      = imyn ; ind%myn = imyn
      iarr_all_mz(nvar%neu%pos)      = imzn ; ind%mzn = imzn
#ifndef ISO
      iarr_all_en(nvar%neu%pos)      = ienn ; ind%enn = ienn
#endif /* !ISO */
#endif /* NEUTRAL */

#ifdef DUST
! Compute index arrays for the dust fluid
      iarr_all_swpx(nvar%dst%beg:nvar%dst%end) = iarr_dst_swpx
      iarr_all_swpy(nvar%dst%beg:nvar%dst%end) = iarr_dst_swpy
      iarr_all_swpz(nvar%dst%beg:nvar%dst%end) = iarr_dst_swpz

      if (selfgrav_dst) then
         i_sg = i_sg + 1
         iarr_all_sg(i_sg) = idnd
      endif
      iarr_all_dn(nvar%dst%pos)      = idnd ; ind%dnd = idnd
      iarr_all_mx(nvar%dst%pos)      = imxd ; ind%mxd = imxd
      iarr_all_my(nvar%dst%pos)      = imyd ; ind%myd = imyd
      iarr_all_mz(nvar%dst%pos)      = imzd ; ind%mzd = imzd
#endif /* DUST */

#ifdef COSM_RAYS
! Compute index arrays for the CR components


      iarr_all_swpx(nvar%crs%beg:nvar%crs%end) = iarr_crs
      iarr_all_swpy(nvar%crs%beg:nvar%crs%end) = iarr_crs
      iarr_all_swpz(nvar%crs%beg:nvar%crs%end) = iarr_crs

      iarr_all_crn(1:nvar%crn%all) = iarr_crn
      iarr_all_cre(1:nvar%cre%all) = iarr_cre
      iarr_all_crs(1:nvar%crs%all) = iarr_crs

      ind%arr_crs => iarr_all_crs

#endif /* COSM_RAYS */

   end subroutine fluid_index

   subroutine cleanup_fluidindex

      implicit none

#ifdef IONIZED
      deallocate(iarr_mag_swpx, iarr_mag_swpy, iarr_mag_swpz, iarr_all_mag)
#endif /* IONIZED */
      deallocate(iarr_all_swpx, iarr_all_swpy, iarr_all_swpz)
      deallocate(iarr_all_dn, iarr_all_mx, iarr_all_my, iarr_all_mz)
      deallocate(iarr_all_sg)
#ifndef ISO
      deallocate(iarr_all_en)
#else /* !ISO */
      deallocate(iarr_all_en)
#endif /* !ISO */

#ifdef COSM_RAYS
      deallocate(iarr_all_crn)
      deallocate(iarr_all_cre)
      deallocate(iarr_all_crs)
#else /* COSM_RAYS */
      deallocate(iarr_all_crn)
      deallocate(iarr_all_cre)
      deallocate(iarr_all_crs)
#endif /* COSM_RAYS */

   end subroutine cleanup_fluidindex

end module fluidindex

