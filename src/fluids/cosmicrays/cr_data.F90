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
!! \brief A place for global constants and parameters related to cosmic rays species
!!
!! Cross sections for spallation from Garcia-Munoz 1987 (see also Longair)
!!
!! Decay half live times from Garcia-Munoz 1987
!!
!! Initial source abundances (in numer density) relative to hydrogen (compare e.g. Longair)
!!
!! \todo make units in this module consistent with units module
!! \todo add user species possibility
!<
module cr_data

! pulled by COSM_RAYS

   implicit none

   public               ! QA_WARN no secrets are kept here

   !! Composition
   ! Isotope list
   !> \deprecated BEWARE: changing the following order should not provide any differences, yet it does!
   enum, bind(C)
      enumerator :: icr_E = 1
      enumerator :: icr_H1
      enumerator :: icr_C12
      enumerator :: icr_N14      !< from decay of Be7 with tau of 0.3 years
      enumerator :: icr_O16
      enumerator :: icr_Li7
      enumerator :: icr_Be9
      enumerator :: icr_Be10

      enumerator :: icr_LAST     !< should be used nowhere despite with nicr
   end enum
   enum, bind(C)
      enumerator :: PRES = 1     !< index for presence of the isotope/component
      enumerator :: ESS          !< index for grad_pcr essentiality of the isotope/component
      enumerator :: SPEC         !< index for energy spectrum treatment of the isotope/component
   end enum

   integer, parameter                                      :: nicr = icr_LAST - 1
   integer                                                 :: ncrsp_auto

   logical, dimension(PRES:SPEC)                           :: eE                 !< presence and grad_pcr essentiality of electrons
   logical, dimension(PRES:SPEC)                           :: eH1                !< presence and grad_pcr essentiality of H1 isotope
   logical, dimension(PRES:SPEC)                           :: eLi7               !< presence and grad_pcr essentiality of Li7 isotope
   logical, dimension(PRES:SPEC)                           :: eBe9               !< presence and grad_pcr essentiality of Be9 isotope
   logical, dimension(PRES:SPEC)                           :: eBe10              !< presence and grad_pcr essentiality of Be10 isotope
   logical, dimension(PRES:SPEC)                           :: eC12               !< presence and grad_pcr essentiality of C12 isotope
   logical, dimension(PRES:SPEC)                           :: eN14               !< presence and grad_pcr essentiality of N14 isotope
   logical, dimension(PRES:SPEC)                           :: eO16               !< presence and grad_pcr essentiality of O16 isotope
   logical,                                dimension(nicr) :: eCRSP              !< table of all isotopes presences
   integer, parameter                                      :: specieslen = 6     !< length of species names
   character(len=specieslen), allocatable, dimension(:)    :: cr_names           !< table of species names
   integer,                   allocatable, dimension(:)    :: cr_table           !< table of cr_data indices for CR species
   integer,                   allocatable, dimension(:)    :: cr_index           !< table of flind indices for CR species
   real,                      allocatable, dimension(:)    :: cr_mass            !< table of mass numbers for CR species
   real,                      allocatable, dimension(:,:)  :: cr_sigma           !< table of cross sections for spallation
   real,                      allocatable, dimension(:)    :: cr_sigma_N         !< table of cross sections for spallation
   real,                      allocatable, dimension(:)    :: cr_tau             !< table of decay half live times
   real,                      allocatable, dimension(:)    :: cr_primary         !< table of initial source abundances
   real,                      allocatable, dimension(:)    :: cr_Z               !< table of atomic numbers
   logical,                   allocatable, dimension(:)    :: cr_spectral        !< table of logicals about energy spectral treatment
   logical,                   allocatable, dimension(:)    :: cr_gpess           !< table of essentiality for grad_pcr calculation
   integer,                   allocatable, dimension(:)    :: icr_spc            !< table of cr_data indices for spectrally resolved CR species
   integer,                   allocatable, dimension(:)    :: iarr_spc           !< table of indices allowing to find spectrally resolved via icr_* (helpful for iarr_crspc)
!<====Mass number and atomic number of nuclei species====>

   real, parameter :: m_H1   = 1.
   real, parameter :: m_Li7  = 7.
   real, parameter :: m_Be9  = 9.
   real, parameter :: m_Be10 = 10.
   real, parameter :: m_C12  = 12.
   real, parameter :: m_N14  = 14.
   real, parameter :: m_O16  = 16.

   real, parameter :: Z_H1   = 1.
   real, parameter :: Z_Li7  = 3.
   real, parameter :: Z_Be9  = 4.
   real, parameter :: Z_Be10 = 4.
   real, parameter :: Z_C12  = 6.
   real, parameter :: Z_N14  = 7.
   real, parameter :: Z_O16  = 8.

!<====Cross sections for spallation from Garcia-Munoz 1987 (see also Longair)====>

   real, parameter :: sigma_C12_Li7  = 10   !< mbarn
   real, parameter :: sigma_C12_Be9  =  6   !< mbarn
   real, parameter :: sigma_C12_Be10 =  3.5 !< mbarn

   real, parameter :: sigma_N14_Li7  =  9.5 !< mbarn

   real, parameter :: sigma_O16_Li7  =  9.5 !< mbarn
   real, parameter :: sigma_O16_Be9  =  4.5 !< mbarn
   real, parameter :: sigma_O16_Be10 =  2   !< mbarn

   real, parameter :: sigma_P_P      =  25. !< mbarn

!<====Decay half live times from Garcia-Munoz 1987====>

   real, parameter :: tau_Be10 = 1.6 !< Myr

!<Initial source abundances (in numer density) relative to hydrogen (compare e.g. Longair)>

   real, parameter :: primary_C12  =  4.5e-3
   real, parameter :: primary_N14  =  1.0e-3
   real, parameter :: primary_O16  =  4.0e-3

   integer, dimension(3), parameter :: icr_prim = [icr_C12, icr_N14, icr_O16], icr_sec = [icr_Li7, icr_Be9, icr_Be10]

contains

!>
!! \brief Routine to set parameters values from namelist CR_SPECIES and auxiliary CR species arrays
!!
!! \n \n
!! @b CR_SPECIES
!! \n \n
!! <table border="+1">
!! <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!! <tr><td>eH1  </td><td>.true.           </td><td>2 logical values</td><td>\copydoc cr_data::eh1  </td></tr>
!! <tr><td>eLi7 </td><td>.false.          </td><td>2 logical values</td><td>\copydoc cr_data::eli7 </td></tr>
!! <tr><td>eBe9 </td><td>[.true., .false.]</td><td>2 logical values</td><td>\copydoc cr_data::ebe9 </td></tr>
!! <tr><td>eBe10</td><td>[.true., .false.]</td><td>2 logical values</td><td>\copydoc cr_data::ebe10</td></tr>
!! <tr><td>eC12 </td><td>[.true., .false.]</td><td>2 logical values</td><td>\copydoc cr_data::ec12 </td></tr>
!! <tr><td>eN14 </td><td>.false.          </td><td>2 logical values</td><td>\copydoc cr_data::en14 </td></tr>
!! <tr><td>eO16 </td><td>.false.          </td><td>2 logical values</td><td>\copydoc cr_data::eo16 </td></tr>
!! </table>
!! The list is active while \b "COSM_RAYS" is defined.
!! \n \n
!<
   subroutine init_cr_species

      use dataio_pub, only: nh
      use mpisetup,   only: lbuff, master, slave, piernik_MPI_Bcast

      implicit none

      namelist /CR_SPECIES/ eE, eH1, eLi7, eBe9, eBe10, eC12, eN14, eO16

      ! Only protons (p+) are dynamically important, we can neglect grad_pcr from heavier nuclei
      ! because of their lower abundancies: n(alpha) ~ 0.1 n(p+), other elements less abundant by orders of magnitude
#ifdef CRESP
      !eE    = [.true., .false., .true.]
#else /* !CRESP */
      eE    = .false.
#endif /* !CRESP */
      eH1   = [.true., .true., .false.]
      !eLi7  = .false.
      !eBe9  = .false.
      eBe10 = .false.
      eC12  = .false.
      !eN14  = .false.
      !eO16  = .false.

#define VS *3-2:3*

      if (master) then

         if (.not.nh%initialized) call nh%init()
         open(newunit=nh%lun, file=nh%tmp1, status="unknown")
         write(nh%lun,nml=CR_SPECIES)
         close(nh%lun)
         open(newunit=nh%lun, file=nh%par_file)
         nh%errstr=""
         read(unit=nh%lun, nml=CR_SPECIES, iostat=nh%ierrh, iomsg=nh%errstr)
         close(nh%lun)
         call nh%namelist_errh(nh%ierrh, "CR_SPECIES")
         read(nh%cmdl_nml,nml=CR_SPECIES, iostat=nh%ierrh)
         call nh%namelist_errh(nh%ierrh, "CR_SPECIES", .true.)
         open(newunit=nh%lun, file=nh%tmp2, status="unknown")
         write(nh%lun,nml=CR_SPECIES)
         close(nh%lun)
         call nh%compare_namelist() ! Do not use one-line if here!

         lbuff(icr_E    VS icr_E   ) = eE
         lbuff(icr_H1   VS icr_H1  ) = eH1
         lbuff(icr_C12  VS icr_C12 ) = eC12
         lbuff(icr_Be9  VS icr_Be9 ) = eBe9
         lbuff(icr_Be10 VS icr_Be10) = eBe10
         lbuff(icr_N14  VS icr_N14 ) = eN14
         lbuff(icr_O16  VS icr_O16 ) = eO16
         lbuff(icr_Li7  VS icr_Li7 ) = eLi7

      endif

      call piernik_MPI_Bcast(lbuff)

      if (slave) then

         eE    = lbuff(icr_E    VS icr_E   )
         eH1   = lbuff(icr_H1   VS icr_H1  )
         eC12  = lbuff(icr_C12  VS icr_C12 )
         eBe9  = lbuff(icr_Be9  VS icr_Be9 )
         eBe10 = lbuff(icr_Be10 VS icr_Be10)
         eN14  = lbuff(icr_N14  VS icr_N14 )
         eO16  = lbuff(icr_O16  VS icr_O16 )
         eLi7  = lbuff(icr_Li7  VS icr_Li7 )

      endif

#undef VS

      eCRSP(1:nicr) = [eE(PRES), eH1(PRES), eC12(PRES), eBe9(PRES), eBe10(PRES), eN14(PRES), eO16(PRES), eLi7(PRES)]
      ncrsp_auto = count(eCRSP)

   end subroutine init_cr_species

   subroutine cr_species_tables(ncrsp, crness)

      use dataio_pub, only: msg, printinfo, warn
      use mpisetup,   only: master
      use units,      only: me, mp, myr, mbarn, sigma_T

      implicit none

      integer(kind=4),       intent(in)    :: ncrsp
      logical, dimension(:), intent(inout) :: crness

      integer                                    :: i, icr, jcr, kcr
      character(len=specieslen), dimension(nicr) :: eCRSP_names
      logical,                   dimension(nicr) :: eCRSP_ess, eCRSP_spec
      real,                      dimension(nicr) :: eCRSP_mass, eCRSP_Z

      eCRSP_names(1:nicr) = ['e-  ', 'p+  ', 'C12 ', 'Be9 ', 'Be10', 'N14 ', 'O16 ', 'Li7 ']
      eCRSP_mass (1:nicr) = [me/mp,  m_H1,   m_C12,   m_Be9, m_Be10, m_N14,  m_O16,  m_Li7 ]
      eCRSP_Z    (1:nicr) = [   0.,  Z_H1,   Z_C12,   Z_Be9, Z_Be10, Z_N14,  Z_O16,  Z_Li7 ]
      eCRSP_ess  (1:nicr) = [eE(ESS) , eH1(ESS) , eC12(ESS) , eBe9(ESS) , eBe10(ESS) , eN14(ESS) , eO16(ESS) , eLi7(ESS) ]
      eCRSP_spec (1:nicr) = [eE(SPEC), eH1(SPEC), eC12(SPEC), eBe9(SPEC), eBe10(SPEC), eN14(SPEC), eO16(SPEC), eLi7(SPEC)]

      allocate(cr_names(ncrsp), cr_table(nicr), cr_index(nicr), cr_sigma(ncrsp,ncrsp), cr_tau(ncrsp), cr_primary(ncrsp), cr_mass(ncrsp), cr_Z(nicr), cr_spectral(ncrsp), cr_gpess(ncrsp),cr_sigma_N(nicr), icr_spc(count(eCRSP_spec .and. eCRSP)), iarr_spc(nicr))
      cr_names(:)    = ''
      cr_table(:)    = 0
      cr_index(:)    = 0
      icr_spc(:)     = 0
      iarr_spc(:)    = 0
      cr_sigma(:,:)  = 0.0
      cr_tau(:)      = 1.0
      cr_primary(:)  = 0.0
      cr_spectral(:) = .false.
      cr_gpess(:)    = .false.

      icr = 0 ; jcr = 0; kcr = 0
      do i = icr_E, size(eCRSP)
         if (eCRSP(i)) then
            icr = icr + 1
            cr_table(i)      = icr
            cr_names(icr)    = eCRSP_names(i)
            cr_mass(icr)     = eCRSP_mass(i)
            cr_Z(icr)        = eCRSP_Z(i)
            cr_spectral(icr) = eCRSP_spec(i)
            if (eCRSP_spec(i)) then
               kcr = kcr + 1
               icr_spc(kcr) = cr_table(i)
               iarr_spc(i)  = kcr
            endif
            cr_sigma_N(icr)  = (cr_Z(icr))**4/(cr_mass(icr)**2)*(me/mp)**2*sigma_T ! Schlickeiser, Cosmic ray astrophysics (2002), formula p.105
            if (icr == icr_E)  cr_sigma_N(icr) = sigma_T
            if (eCRSP_spec(i)) then
               if (i /= icr_E) then
                  write(msg, '(3a)') "[cr_data:init_cr_species] Energy spectral treatment for ", eCRSP_names(i), " under development, results will not be reliable."
                  call warn(msg)
               endif
            else
               jcr = jcr + 1
               cr_index(i) = jcr
            endif
            cr_gpess(icr)    = eCRSP_ess(i)
            if (master) then
               write(msg,'(a,a,a,l2)') spectral_or_not(cr_spectral(icr)), eCRSP_names(i), 'CR species is present; taken into account for grad_pcr: ', eCRSP_ess(i)
               call printinfo(msg)
            endif
         endif
      enddo

      if (ncrsp_auto < ncrsp) then
         cr_gpess(ncrsp_auto+1:ncrsp) = crness
         do i = ncrsp_auto+1, ncrsp
            write(msg,'(a,a,i2,a,l2)') spectral_or_not(cr_spectral(i)), 'user CR species no: ', i,' is present; taken into account for grad_pcr: ', cr_gpess(i)
            if (master) call printinfo(msg)
         enddo
      endif

      if (master) then
         write(msg, '(a,i3)') 'Total amount of CR species: ', ncrsp
         if (count(cr_spectral) > 0) write(msg, '(a,a,i3,a)') trim(msg), ' | ', count(cr_spectral), ' spectral component(s)'
         if (ncrsp - ncrsp_auto > 0) write(msg, '(a,a,i3,a)') trim(msg), ' | ', ncrsp - ncrsp_auto, ' user component(s).'
         call printinfo(msg)
      endif

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
      cr_sigma = cr_sigma * mbarn
      if (eCRSP(icr_Be10)) cr_tau(cr_table(icr_Be10)) = tau_Be10 * myr

   end subroutine cr_species_tables

   function spectral_or_not(sp) result(wr)

      implicit none

      logical, intent(in) :: sp
      integer, parameter  :: spl = 13
      character(len=spl)  :: wr

      if (sp) then
         wr = '    spectral '
      else
         wr = 'non-spectral '
      endif

   end function spectral_or_not

!> \brief cleanup routine

   subroutine cleanup_cr_species

      implicit none

      deallocate(cr_names, cr_table, cr_sigma, cr_tau, cr_primary, cr_mass)

   end subroutine cleanup_cr_species

end module cr_data

! this type looks useful but is unused.
!!$   integer, parameter :: isoname_len = 8
!!$   type :: cr_component
!!$      character(len=isoname_len) :: isotope     !< isotope name, eg. Be10
!!$      integer          :: index=      !< relative index (with respect to crn_beg)
!!$      real             :: abund=      !< initial abundance relative to H
!!$   end type cr_component
