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
!! \brief Definitions of compound types for fluids
!<
module fluidtypes
! pulled by ANY
   use types,     only: value

   implicit none

   private
   public :: component, component_fluid, phys_prop, var_numbers

   type :: phys_prop
      type(value) :: dens_min
      type(value) :: dens_max
      type(value) :: velx_max
      type(value) :: vely_max
      type(value) :: velz_max
      type(value) :: pres_max
      type(value) :: pres_min
      type(value) :: temp_max
      type(value) :: temp_min
      type(value) :: cs_max
   end type phys_prop

   type :: component
      integer(kind=4) :: all = 0   !< number of all variables in fluid/component
      integer(kind=4) :: beg = 0   !< beginning number of variables in fluid/component
      integer(kind=4) :: end = 0   !< end number of variables in fluid/component
      integer(kind=4) :: pos = 0   !< index denoting position of the fluid in the row of fluids
   end type component

   type, extends(component) :: component_fluid
      integer(kind=4) :: idn = -1      !< index denoting position of the fluid density in array arrays::u
      integer(kind=4) :: imx = -1      !< index denoting position of the fluid x-momentum in array arrays::u
      integer(kind=4) :: imy = -1      !< index denoting position of the fluid y-momentum in array arrays::u
      integer(kind=4) :: imz = -1      !< index denoting position of the fluid z-momentum in array arrays::u
      integer(kind=4) :: ien = -1      !< index denoting position of the fluid energy in array arrays::u

      real    :: cs    = 0.0   !< fluid's isothermal sound speed
      real    :: cs2   = 0.0   !< fluid's isothermal sound speed squared
      real    :: gam   = -1.0  !< fluid's adiabatic index
      real    :: gam_1 = -1.0  !< fluid's adiabatic index minus one

      logical :: is_selfgrav   = .false. !< True if fluid is selfgravitating
      logical :: is_magnetized = .false. !< True if fluid is magnetized
      logical :: has_energy    = .false. !< True if fluid has additional energy array

      integer :: tag = -1 !< Human readable tag describing fluid

      integer(kind=4), allocatable, dimension(:)   :: iarr
      integer(kind=4), allocatable, dimension(:,:) :: iarr_swp

      type(phys_prop) :: snap

      real :: c !< COMMENT ME (this quantity was previously a member of phys_prop, but is used in completely different way than other phys_prop% members

      contains
         procedure :: set_cs  => update_sound_speed
         procedure :: set_gam => update_adiabatic_index
         procedure :: info    => printinfo_component_fluid
         procedure :: get_tag
   end type component_fluid

   type :: var_numbers
      integer(kind=4) :: all         = 0      !< total number of fluid variables = the size of array \a u(:,:,:,:) in the first index
      integer(kind=4) :: fluids      = 0      !< number of fluids (ionized gas, neutral gas, dust)
      integer(kind=4) :: energ       = 0      !< number of non-isothermal fluids (indicating the presence of energy density in the vector of conservative variables)
      integer(kind=4) :: components  = 0      !< number of components, such as CRs, tracers, magnetic helicity (in future), whose formal description does not involve [???]
      integer (kind=4):: fluids_sg   = 0      !< number of selfgravitating fluids (ionized gas, neutral gas, dust)

      type(component_fluid), dimension(:), pointer :: all_fluids

      type(component_fluid), pointer :: ion         !< numbers of variables for the ionized fluid
      type(component_fluid), pointer :: neu         !< numbers of variables for the neutral fluid
      type(component_fluid), pointer :: dst         !< numbers of variables for the dust fluid

      !> \todo those vars should be converted to pointers
      type(component) :: crs         !< numbers of variables in all cosmic ray components
      type(component) :: crn         !< numbers of variables in cosmic ray nuclear components
      type(component) :: cre         !< numbers of variables in cosmic ray electron components
   end type var_numbers

   contains

      subroutine update_adiabatic_index(this,new_gamma)
         use dataio_pub, only: warn
         implicit none
         class(component_fluid) :: this
         real, intent(in)       :: new_gamma

         if (.not.this%has_energy) then
            call warn("Fluid does not have energy component")
            call warn("Updating gamma does not make much sense o.O")
         endif
         this%gam   = new_gamma
         this%gam_1 = new_gamma-1.0
      end subroutine update_adiabatic_index

      subroutine update_sound_speed(this,new_cs)
         implicit none
         class(component_fluid) :: this
         real, intent(in)       :: new_cs

         this%cs  = new_cs
         this%cs2 = new_cs**2
      end subroutine update_sound_speed

      function get_tag(this) result(tag)
         use constants, only: ION, NEU, DST, idlen
         implicit none
         class(component_fluid) :: this
         character(len=idlen)   :: tag

         select case (this%tag)
            case (ION)
               tag = "ION"
            case (NEU)
               tag = "NEU"
            case (DST)
               tag = "DST"
            case default
               tag = "---"
         end select
      end function get_tag

      subroutine printinfo_component_fluid(this)
         use dataio_pub,  only: msg, printinfo
         implicit none
         class(component_fluid) :: this

         write(msg,*) "idn   = ", this%idn;     call printinfo(msg)
         write(msg,*) "imx   = ", this%imx;     call printinfo(msg)
         write(msg,*) "imy   = ", this%imy;     call printinfo(msg)
         write(msg,*) "imz   = ", this%imz;     call printinfo(msg)
         write(msg,*) "ien   = ", this%ien;     call printinfo(msg)

         write(msg,*) "cs    = ", this%cs;      call printinfo(msg)
         write(msg,*) "cs2   = ", this%cs2;     call printinfo(msg)
         write(msg,*) "gam   = ", this%gam;     call printinfo(msg)
         write(msg,*) "gam_1 = ", this%gam_1;   call printinfo(msg)

         if (this%is_selfgrav) then
            write(msg,*) "Fluid is selfgravitating"; call printinfo(msg)
         endif
         if (this%is_magnetized) then
            write(msg,*) "Fluid is magnetized";      call printinfo(msg)
         endif
         if (this%has_energy) then
            write(msg,*) "Fluid has energy";         call printinfo(msg)
         endif
         write(msg,*) "TAG   = ", this%tag;     call printinfo(msg)

      end subroutine printinfo_component_fluid

end module fluidtypes
