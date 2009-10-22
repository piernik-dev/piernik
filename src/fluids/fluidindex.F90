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

!>
!! \brief The purpose of this module is to compute: 
!! \n (1) positions of all fluid variables to be referenced through the first index of array u(:,:,:,:), 
!! \n (2) arrays of indexes to make easy reference to all gas densities, x,y,z-components of momenta, energy densities, CR energy densities, transposed components of x,y,z-momenta in directional sweeps, etc ...
!! 
!<

module fluidindex

   use types   
   implicit none

   type(indx),save     :: ind         !< derived type variable storing all fluid-, component- and magnetic-field-indexes
   
   integer, parameter  :: nmag = 3    !< number of magnetic field components                  
   integer             :: nvar 	      !< total number of fluid variables = the size of array \a u(:,:,;,:) in the first index
   integer             :: nfluid      !< number of fluids (ionized gas, neutral gas, dust)
   integer             :: nadiab      !< number of adiabatic fluids (indicating the presence of energy density in the vector of conservative variables for adiabatic fluids)
   integer             :: ncomponents !< number of components, such as CRs, tracers, magnetic helicity (in future),  whose formal description does not involve momenta and mass density

   integer, parameter  :: ibx = 1     !< index of x-component of magnetic field
   integer, parameter  :: iby = 2     !< index of y-component of magnetic field 
   integer, parameter  :: ibz = 3     !< index of z-component of magnetic field
   integer, parameter  :: idn = 1     !< position of density in the vector of conserv. variables for single fluid
   integer, parameter  :: imx = 2     !< position of x-mom. in the vector of conserv. variables for single fluid
   integer, parameter  :: imy = 3     !< position of y-mom. in the vector of conserv. variables for single fluid
   integer, parameter  :: imz = 4     !< position of z-mom. in the vector of conserv. variables for single fluid
   integer, parameter  :: ien = 5     !< position of energy density in the vector of conserv. variables for single fluid (only for adiabatic fluids)
   integer, parameter  :: icr = 1     !< position of CR-energy density in the vector of cons. variables for CR multiple components (may change if future, when multiple CR species will be taken into account)

#ifdef IONIZED
#ifdef RESISTIVE
   integer, parameter  :: icx = 1     !< index of x-component of current density
   integer, parameter  :: icy = 2     !< index of y-component of current density
   integer, parameter  :: icz = 3     !< index of z-component of current density
#endif /* RESISTIVE */
#endif /* IONIZED */

   
   integer, allocatable, dimension(:) :: iarr_all_dn   !< array of indexes pointing to mass densities of all fluids 
   integer, allocatable, dimension(:) :: iarr_all_mx   !< array of indexes pointing to mom. densities of all fluids
   integer, allocatable, dimension(:) :: iarr_all_my   !< array of indexes pointing to mom. densities of all fluids
   integer, allocatable, dimension(:) :: iarr_all_mz   !< array of indexes pointing to mom. densities of all fluids
   integer, allocatable, dimension(:) :: iarr_all_en   !< array of indexes pointing to ener. densities of all fluids
   integer, allocatable, dimension(:) :: iarr_all_cr   !< array of indexes pointing to ener. densities of all CR-components
   integer, allocatable, dimension(:) :: iarr_all_swpx !< array (size = nvar) of all fluid indexes in the original order
   integer, allocatable, dimension(:) :: iarr_all_swpy !< array (size = nvar) of all fluid indexes with \a x and \a y components of mom. interchanged 
   integer, allocatable, dimension(:) :: iarr_all_swpz !< array (size = nvar) of all fluid indexes with \a x and \a z components of mom. interchanged 

#ifdef IONIZED
   integer, allocatable, dimension(:) :: iarr_all_mag  !< array (size = nmag) of all magnetic field components
   integer, allocatable, dimension(:) :: iarr_mag_swpx !< array (size = nmag) of all mag. field indexes in the original order (same as iarr_all_mag)
   integer, allocatable, dimension(:) :: iarr_mag_swpy !< array (size = nmag) of all mag. field indexes \a x and \a z components interchanged
   integer, allocatable, dimension(:) :: iarr_mag_swpz !< array (size = nmag) of all mag. field indexes \a x and \a z components interchanged
#endif /* IONIZED */

#ifdef IONIZED
   integer :: nvar_ion                                 !< number of variables for the ion fluid (4 when iso, 5 when adiab)
   integer :: beg_ion                                  !< index of the first ion component (density) in the full vector of cons. variables
   integer :: end_ion                                  !< index of the last ion component (z-momentum or en. density) in the full vector of cons. variables
   integer :: i_ion                                    !< index denoting position of the ion fluid in the row of fluids 
#endif /* IONIZED */

#ifdef NEUTRAL
   integer :: nvar_neu                                 !< number of variables for the neutral fluid (4 when iso, 5 when adiab)
   integer :: beg_neu                                  !< index of the first neutral component (density) in the full vector of cons.
   integer :: end_neu                                  !< index of the last neutral component (z-momentum or en. density) in the full vector of cons. variables
   integer :: i_neu                                    !< index denoting position of the neutral fluid in the row of fluids
#endif /* NEUTRAL */

#ifdef DUST
   integer :: nvar_dst                                 !< number of variables for the dust fluid (4)
   integer :: beg_dst                                  !< index of the first dust component (density) in the full vector of cons.
   integer :: end_dst                                  !< index of the last dust component (z-momentum or en. density) in the full vector of cons. variables
   integer :: i_dst                                    !< index denoting position of the dust fluid in the row of fluids
#endif /* DUST */
   
#ifdef COSM_RAYS
   integer :: nvar_crs                                 !< number of variables for the CR component 
   integer :: beg_crs                                  !< index of the first CR component (density) in the full vector of cons.
   integer :: end_crs                                  !< index of the last CR component (z-momentum or en. density) in the full vector of cons. variables
   integer :: i_crs                                    !< index denoting position of CRs in the row of fluids and components
#endif /* COSM_RAYS */


  contains
   
    subroutine fluid_index
#ifdef IONIZED
      use initionized,    only : ionized_index
      use initionized,    only : iarr_ion_swpx, iarr_ion_swpy, iarr_ion_swpz
      use initionized,    only : idni,imxi,imyi,imzi
#ifndef ISO
      use initionized,    only : ieni
#endif /* ISO */ 
#endif /* IONIZED */

#ifdef NEUTRAL
      use initneutral,    only : neutral_index
      use initneutral,    only : iarr_neu_swpx, iarr_neu_swpy, iarr_neu_swpz
      use initneutral,    only : idnn,imxn,imyn,imzn
#ifndef ISO 
      use initneutral,    only : ienn
#endif /* ISO */      
#endif /* NEUTRAL */

#ifdef DUST
      use initdust,       only : dust_index
      use initdust,       only : iarr_dst_swpx, iarr_dst_swpy, iarr_dst_swpz
      use initdust,       only : idnd,imxd,imyd,imzd      
#endif /* DUST */

#ifdef COSM_RAYS
      use initcosmicrays, only : cosmicray_index
      use initcosmicrays, only : iarr_crs_swpx, iarr_crs_swpy, iarr_crs_swpz
      use initcosmicrays, only : iecr           
#endif /* COSM_RAYS */

      implicit none


      nvar   = 0
      ncomponents  = 0
      nfluid = 0
      nadiab = 0
      
     
#ifdef IONIZED
      nvar_ion  = 0
      beg_ion   = nvar + 1
      call ionized_index(nvar,nvar_ion) 
      end_ion   = nvar  
      ncomponents  = ncomponents + 1
      nfluid = nfluid + 1  
      i_ion = ncomponents 
#ifndef ISO 
      nadiab = nadiab + 1
#endif /* ISO */      
#endif /* IONIZED */

#ifdef NEUTRAL
      nvar_neu = 0
      beg_neu = nvar + 1
      call neutral_index(nvar,nvar_neu) 
      end_neu = nvar    
      ncomponents  = ncomponents + 1
      nfluid = nfluid + 1   
      i_neu = ncomponents 
#ifndef ISO 
      nadiab = nadiab + 1
#endif /* ISO */      
#endif /* NEUTRAL */

#ifdef DUST
      nvar_dst = 0
      beg_dst = nvar + 1
      call dust_index(nvar,nvar_dst) 
      end_dst = nvar    
      ncomponents  = ncomponents + 1
      nfluid = nfluid + 1   
      i_dst = ncomponents 
#endif /* DUST */

#ifdef COSM_RAYS
      nvar_crs   = 0
      beg_crs = nvar + 1
      call cosmicray_index(nvar,nvar_crs)  
      end_crs = nvar   
      ncomponents  = ncomponents + 1   
      i_crs = ncomponents                
#endif /* COSM_RAYS */     

#ifdef IONIZED
      allocate(iarr_mag_swpx(nmag),iarr_mag_swpy(nmag),iarr_mag_swpz(nmag),iarr_all_mag(nmag))
#endif /* IONIZED */
      allocate(iarr_all_swpx(nvar),iarr_all_swpy(nvar),iarr_all_swpz(nvar))
      allocate(iarr_all_dn(nfluid),iarr_all_mx(nfluid),iarr_all_my(nfluid),iarr_all_mz(nfluid))
#ifndef ISO      
      allocate(iarr_all_en(nadiab))  
#else
      allocate(iarr_all_en(0))
#endif /* ISO */
#ifdef COSM_RAYS      
      allocate(iarr_all_cr(nvar_crs))  
#else
      allocate(iarr_all_cr(0))
#endif /* COSM_RAYS */

#ifdef IONIZED
      iarr_mag_swpx = [ibx,iby,ibz]
      iarr_mag_swpy = [iby,ibx,ibz]
      iarr_mag_swpz = [ibz,iby,ibx]
      iarr_all_mag  = [ibx,iby,ibz]
      ind%bx = ibx; ind%by = iby; ind%bz = ibz
#endif /* IONIZED */


#ifdef IONIZED
      iarr_all_swpx(beg_ion:end_ion) = iarr_ion_swpx
      iarr_all_swpy(beg_ion:end_ion) = iarr_ion_swpy
      iarr_all_swpz(beg_ion:end_ion) = iarr_ion_swpz    

      iarr_all_dn(i_ion)      = idni ; ind%dni = idni
      iarr_all_mx(i_ion)      = imxi ; ind%mxi = imxi
      iarr_all_my(i_ion)      = imyi ; ind%myi = imyi
      iarr_all_mz(i_ion)      = imzi ; ind%mzi = imzi
#ifndef ISO
      iarr_all_en(i_ion)      = ieni ; ind%eni = ieni
#endif /* ISO */
#endif /* IONIZED */
      
#ifdef NEUTRAL
      iarr_all_swpx(beg_neu:end_neu) = iarr_neu_swpx
      iarr_all_swpy(beg_neu:end_neu) = iarr_neu_swpy
      iarr_all_swpz(beg_neu:end_neu) = iarr_neu_swpz    

      iarr_all_dn(i_neu)      = idnn ; ind%dnn = idnn
      iarr_all_mx(i_neu)      = imxn ; ind%mxn = imxn
      iarr_all_my(i_neu)      = imyn ; ind%myn = imyn
      iarr_all_mz(i_neu)      = imzn ; ind%mzn = imzn
#ifndef ISO
      iarr_all_en(i_neu)      = ienn ; ind%enn = ienn
#endif /* ISO */
#endif /* NEUTRAL */
      
#ifdef DUST
      iarr_all_swpx(beg_dst:end_dst) = iarr_dst_swpx
      iarr_all_swpy(beg_dst:end_dst) = iarr_dst_swpy
      iarr_all_swpz(beg_dst:end_dst) = iarr_dst_swpz    

      iarr_all_dn(i_dst)      = idnd ; ind%dnd = idnd
      iarr_all_mx(i_dst)      = imxd ; ind%mxd = imxd
      iarr_all_my(i_dst)      = imyd ; ind%myd = imyd
      iarr_all_mz(i_dst)      = imzd ; ind%mzd = imzd
#endif /* DUST */
      
#ifdef COSM_RAYS      
      iarr_all_swpx(beg_crs:end_crs) = iarr_crs_swpx
      iarr_all_swpy(beg_crs:end_crs) = iarr_crs_swpy
      iarr_all_swpz(beg_crs:end_crs) = iarr_crs_swpz    

      iarr_all_cr(1:nvar_crs) = iecr  ; ind%ecr = iecr
#endif /* COSM_RAYS */
   
!      write(*,*) 'fluid_index', iarr_ion_swpx
!      write(*,*) 'fluid_index', iarr_all_swpx

   end subroutine fluid_index


end module fluidindex

