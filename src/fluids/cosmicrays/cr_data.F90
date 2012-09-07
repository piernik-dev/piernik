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
#include "macros.h"
!>
!! \brief A place for global constants related to cosmic rays
!! \todo make units in this module consistent with units module
!! Cross sections for spallation from Garcia-Munoz 1987 (see also Longair)
!! Decay half live times from Garcia-Munoz 1987
!! Initial source abundances (in numer density) relative to hydrogen (compare e.g. Longair)
!<
module cr_data

! pulled by COSM_RAYS_SOURCES

   implicit none

   public               ! QA_WARN no secrets are kept here

   !! Composition
   ! Isotope list
   !> \deprecated BEWARE: changing the following order should not provide any differences, yet it does!
   enum, bind(C)
      enumerator :: icr_H1 = 1
      enumerator :: icr_C12
      enumerator :: icr_Be9
      enumerator :: icr_Be10
      enumerator :: icr_N14      !< from decay of Be7 with tau of 0.3 years
      enumerator :: icr_O16
      enumerator :: icr_Li7      !< \deprecated BEWARE: ncrn should be set up gerater than maximum isotope numeber, which should be smaller than ncr_max (<10 currently)
      enumerator :: icr_LAST     !< should be used nowhere despite with nicr
   end enum

   integer, parameter                                      :: nicr = icr_LAST - 1

   logical                                                 :: eH1, eLi7, eBe9, eBe10, eC12, eN14, eO16
   logical,                                dimension(nicr) :: eCRSP
   integer, parameter                                      :: specieslen = 6
   character(len=specieslen), allocatable, dimension(:)    :: cr_names
   integer,                   allocatable, dimension(:)    :: cr_table
   real,                      allocatable, dimension(:,:)  :: cr_sigma
   real,                      allocatable, dimension(:)    :: cr_tau
   real,                      allocatable, dimension(:)    :: cr_primary

!<====Cross sections for spallation from Garcia-Munoz 1987 (see also Longair)====>

   real, parameter :: Myear=1d6*365*24*60*60 !< s     \deprecated BEWARE: this line breaks unit consistency, move it to units.F90 and use scaling
   real, parameter :: mbarn=1e-27            !< cm2   \deprecated BEWARE: this line breaks unit consistency, move it to units.F90 and use scaling

   real, parameter :: sigma_C12_Li7  = 10   * mbarn
   real, parameter :: sigma_C12_Be9  =  6   * mbarn
   real, parameter :: sigma_C12_Be10 =  3.5 * mbarn

   real, parameter :: sigma_N14_Li7  =  9.5 * mbarn

   real, parameter :: sigma_O16_Li7  =  9.5 * mbarn
   real, parameter :: sigma_O16_Be9  =  4.5 * mbarn
   real, parameter :: sigma_O16_Be10 =  2   * mbarn

!<====Decay half live times from Garcia-Munoz 1987====>

   real, parameter :: tau_Be10 = 1.6 !< Myr \deprecated BEWARE: this line breaks unit consistency, move it to units.F90 and use scaling

!<Initial source abundances (in numer density) relative to hydrogen (compare e.g. Longair)>

   real, parameter :: primary_C12  =  4.5e-3
   real, parameter :: primary_N14  =  1.0e-3
   real, parameter :: primary_O16  =  4.0e-3

contains

!>
!! \todo add user species possibility
!<
   subroutine init_crsources(ncrn)

      use dataio_pub,      only: par_file, ierrh, namelist_errh, compare_namelist, cmdl_nml, lun   ! QA_WARN required for diff_nml
      use dataio_pub,      only: msg, printinfo
      use mpisetup,        only: lbuff, master, slave, piernik_MPI_Bcast

      implicit none

      integer, intent(in)                        :: ncrn
      integer                                    :: icr, i
      character(len=specieslen), dimension(nicr) :: eCRSP_names

      namelist /CR_SPECIES/ eH1, eLi7, eBe9, eBe10, eC12, eN14, eO16

      eH1   = .true.
      eLi7  = .false.
      eBe9  = .true.
      eBe10 = .true.
      eC12  = .true.
      eN14  = .false.
      eO16  = .false.

      if (master) then

         diff_nml(CR_SPECIES) ! Do not use one-line if here!

         lbuff(icr_H1)   = eH1
         lbuff(icr_C12)  = eC12
         lbuff(icr_Be9)  = eBe9
         lbuff(icr_Be10) = eBe10
         lbuff(icr_N14)  = eN14
         lbuff(icr_O16)  = eO16
         lbuff(icr_Li7)  = eLi7

      endif

      call piernik_MPI_Bcast(lbuff)

      if (slave) then

         eH1   = lbuff(icr_H1)
         eC12  = lbuff(icr_C12)
         eBe9  = lbuff(icr_Be9)
         eBe10 = lbuff(icr_Be10)
         eN14  = lbuff(icr_N14)
         eO16  = lbuff(icr_O16)
         eLi7  = lbuff(icr_Li7)

      endif

      eCRSP(1:7)       = [eH1, eC12, eBe9, eBe10, eN14, eO16, eLi7]
      eCRSP_names(1:7) = ['H1  ','C12 ','Be9 ','Be10','N14 ','O16 ','Li7 ']
      allocate(cr_names(ncrn), cr_table(nicr), cr_sigma(ncrn,ncrn), cr_tau(ncrn), cr_primary(ncrn))
      cr_names(:)   = ''
      cr_table(:)   = 0
      cr_sigma(:,:) = 0.0
      cr_tau(:)     = 1.0
      cr_primary(:) = 0.0

      icr = 0
      do i = 1, size(eCRSP)
         if (eCRSP(i)) then
            icr = icr + 1
            cr_table(i)   = icr
            cr_names(icr) = eCRSP_names(i)
            write(msg,'(a,a)') eCRSP_names(i), 'CR species is present'
            call printinfo(msg)
         endif
      enddo

      if (eCRSP(icr_C12)) then
         cr_primary(cr_table(icr_C12)) = primary_C12
         if (eCRSP(icr_Li7 )) cr_sigma(cr_table(icr_C12), cr_table(icr_Li7 )) = sigma_C12_Li7
         if (eCRSP(icr_Be9 )) cr_sigma(cr_table(icr_C12), cr_table(icr_Be9 )) = sigma_C12_Be9
         if (eCRSP(icr_Be10)) cr_sigma(cr_table(icr_C12), cr_table(icr_Be10)) = sigma_C12_Be10
      endif
      if (eCRSP(icr_N14)) then
         cr_primary(cr_table(icr_N14)) = primary_N14
         if (eCRSP(icr_Li7 )) cr_sigma(cr_table(icr_N14), cr_table(icr_Li7 )) = sigma_N14_Li7
      endif
      if (eCRSP(icr_O16)) then
         cr_primary(cr_table(icr_O16)) = primary_O16
         if (eCRSP(icr_Li7 )) cr_sigma(cr_table(icr_O16), cr_table(icr_Li7 )) = sigma_O16_Li7
         if (eCRSP(icr_Be9 )) cr_sigma(cr_table(icr_O16), cr_table(icr_Be9 )) = sigma_O16_Be9
         if (eCRSP(icr_Be10)) cr_sigma(cr_table(icr_O16), cr_table(icr_Be10)) = sigma_O16_Be10
      endif
      if (eCRSP(icr_Be10)) cr_tau(cr_table(icr_Be10)) = tau_Be10


   end subroutine init_crsources

end module cr_data

! this type looks useful but is unused.
!!$   integer, parameter :: isoname_len = 8
!!$   type :: cr_component
!!$      character(len=isoname_len) :: isotope     !< isotope name, eg. Be10
!!$      integer          :: index=      !< relative index (with respect to crn_beg)
!!$      real             :: abund=      !< initial abundance relative to H
!!$   end type cr_component
