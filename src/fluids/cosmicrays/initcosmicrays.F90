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
!! \brief Initialization of Cosmic Ray component
!!
!!
!!
!! In this module following namelist of parameters is specified:
!! \copydetails initcosmicrays::init_cosmicrays
!<
module initcosmicrays
! pulled by COSM_RAYS
   use constants, only: cbuff_len
   implicit none

   public ! QA_WARN no secrets are kept here
   private :: cbuff_len ! QA_WARN prevent reexport

   integer, parameter                  :: ncr_max = 102  !< maximum number of CR nuclear and electron components (\warning higher ncr_max limit would require changes in names of components in common_hdf5)
   ! namelist parameters
   integer(kind=4)                     :: ncrsp        !< number of CR components \deprecated BEWARE: ncrtot (sum of ncrsp and ncr2b) should not be higher than ncr_max = 102
   integer(kind=4)                     :: ncr_user     !< number of CR components supplementary in respect of listed in CR_SPECIES
   integer(kind=4)                     :: ncrn         !< number of CR non-spectral components
   integer(kind=4)                     :: nspc         !< number of CR spectral components
   integer(kind=4)                     :: ncrb         !< number of bins for CRESP
   integer(kind=4)                     :: ncr2b        !< 2*ncrb for CRESP
   integer(kind=4)                     :: ncrtot       !< number of all CR components \deprecated BEWARE: ncrtot (sum of ncrsp and ncr2b) should not be higher than ncr_max = 102
   real                                :: cfl_cr       !< CFL number for diffusive CR transport
   real                                :: smallecr     !< floor value for CR energy density
   real                                :: cr_active    !< parameter specifying whether CR pressure gradient is (when =1.) or isn't (when =0.) included in the gas equation of motion
   real                                :: cr_eff       !< conversion rate of SN explosion energy to CR energy (default = 0.1)
   real                                :: gamma_cr     !< adiabatic index of all CR non-spectral components
   real                                :: gamma_cr_1   !< gamma_cr - 1
   logical                             :: use_CRdiff   !< switch for diffusion of cosmic rays
   logical                             :: use_CRdecay  !< switch for spallation and decay of cosmic rays
   logical                             :: use_smallecr !< correct CR energy density when it gets lower than smallecr
   character(len=cbuff_len)            :: divv_scheme  !< scheme used to calculate div(v), see crhelpers for more details
   real, dimension(ncr_max)            :: K_cr_paral   !< array containing parallel diffusion coefficients of all CR nuclear components or maximal parallel diffusion coefficient value for CRESP
   real, dimension(ncr_max)            :: K_cr_perp    !< array containing perpendicular diffusion coefficients of all CR nuclear components or maximal perpendicular diffusion coefficient value for CRESP
   logical, dimension(ncr_max)         :: gpcr_ess_user !< if user CR species is essential for grad_pcr calculation
   integer(kind=4), allocatable, dimension(:) :: gpcr_ess_noncresp !< indexes of essentials for grad_pcr calculation for non-CRESP components
   ! public component data
   integer(kind=4), allocatable, dimension(:) :: iarr_crn !< array of indexes pointing to all CR nuclear components
   integer(kind=4), allocatable, dimension(:) :: iarr_crspc !< array of indexes pointing to all CR electron components
   integer(kind=4), allocatable, dimension(:) :: iarr_crs !< array of indexes pointing to all CR components
#ifdef CRESP
   integer(kind=4), allocatable, dimension(:) :: iarr_crspc_e !< array of indexes pointing to all CR electron energy components
   integer(kind=4), allocatable, dimension(:) :: iarr_crspc_n !< array of indexes pointing to all CR electron number density components
   integer(kind=4), allocatable, dimension(:,:) :: iarr_crspc2_e !< 2D iterable array of indexes pointing to spectrally resolved CR components' energy density
   integer(kind=4), allocatable, dimension(:,:) :: iarr_crspc2_n !< 2D iterable array of indexes pointing to spectrally resolved CR components' number density
#endif /* CRESP */

   real,    allocatable, dimension(:)  :: K_crs_paral  !< array containing parallel diffusion coefficients of all CR components
   real,    allocatable, dimension(:)  :: K_crs_perp   !< array containing perpendicular diffusion coefficients of all CR components
   !> \deprecated BEWARE Possible confusion: *_perp coefficients are not "perpendicular" but rather isotropic
   real                                :: def_dtcrs    !< default dt limitation due to diffusion
   logical                             :: K_crs_valid  !< condition to use dt_crs

contains

!>
!! \brief Routine to set parameters values from namelist COSMIC_RAYS
!!
!! \n \n
!! @b COSMIC_RAYS
!! \n \n
!! <table border="+1">
!! <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!! <tr><td>cfl_cr      </td><td>0.9    </td><td>real value</td><td>\copydoc initcosmicrays::cfl_cr     </td></tr>
!! <tr><td>smallecr    </td><td>0.0    </td><td>real value</td><td>\copydoc initcosmicrays::smallecr   </td></tr>
!! <tr><td>cr_active   </td><td>1.0    </td><td>real value</td><td>\copydoc initcosmicrays::cr_active  </td></tr>
!! <tr><td>cr_eff      </td><td>0.1    </td><td>real value</td><td>\copydoc initcosmicrays::cr_eff     </td></tr>
!! <tr><td>use_CRdiff  </td><td>.true. </td><td>logical   </td><td>\copydoc initcosmicrays::use_CRdiff </td></tr>
!! <tr><td>use_CRdecay </td><td>.false.</td><td>logical   </td><td>\copydoc initcosmicrays::use_CRdecay</td></tr>
!! <tr><td>ncr_user    </td><td>0      </td><td>integer   </td><td>\copydoc initcosmicrays::ncr_user   </td></tr>
!! <tr><td>ncrb        </td><td>0      </td><td>integer   </td><td>\copydoc initcosmicrays::ncrb       </td></tr>
!! <tr><td>gamma_cr    </td><td>4./3.  </td><td>real array</td><td>\copydoc initcosmicrays::gamma_cr   </td></tr>
!! <tr><td>K_cr_paral  </td><td>0      </td><td>real array</td><td>\copydoc initcosmicrays::k_cr_paral </td></tr>
!! <tr><td>K_cr_perp   </td><td>0      </td><td>real array</td><td>\copydoc initcosmicrays::k_cr_perp  </td></tr>
!! <tr><td>divv_scheme </td><td>''     </td><td>string    </td><td>\copydoc initcosmicrays::divv_scheme</td></tr>
!! <tr><td>gpcr_ess_user</td><td>.false.</td><td>logical array</td><td>\copydoc initcosmicrays::gpcr_ess_user</td></tr>
!! </table>
!! The list is active while \b "COSM_RAYS" is defined.
!! \n \n
!<
   subroutine init_cosmicrays

      use constants,       only: cbuff_len, I_ONE, I_TWO, half, big
      use cr_data,         only: init_cr_species, cr_species_tables, cr_gpess, cr_spectral, ncrsp_auto
      use diagnostics,     only: ma1d, ma2d, my_allocate
      use dataio_pub,      only: die, warn, nh
      use func,            only: operator(.notequals.)
      use mpisetup,        only: ibuff, rbuff, lbuff, cbuff, master, slave, piernik_MPI_Bcast

      implicit none

      integer(kind=4) :: nl, nn, icr
      real            :: maxKcrs

      namelist /COSMIC_RAYS/ cfl_cr, use_smallecr, smallecr, cr_active, cr_eff, use_CRdiff, use_CRdecay, divv_scheme, &
           &                 gamma_cr, K_cr_paral, K_cr_perp, ncr_user, ncrb, gpcr_ess_user

      call init_cr_species

      cfl_cr         = 0.9
      smallecr       = 0.0
      cr_active      = 1.0
      cr_eff         = 0.1       !  canonical conversion rate of SN en.-> CR (e_sn=10**51 erg)
      ncrsp          = ncrsp_auto
      ncr_user       = 0
      ncrb           = 0

      use_CRdiff     = .true.
      use_CRdecay    = .false.
      use_smallecr   = .true.

      gamma_cr       = 4./3.
      K_cr_paral(:)  = 0.0
      K_cr_perp(:)   = 0.0

      gpcr_ess_user  = .false.

      divv_scheme    = ''

      if (master) then

         if (.not.nh%initialized) call nh%init()
         open(newunit=nh%lun, file=nh%tmp1, status="unknown")
         write(nh%lun,nml=COSMIC_RAYS)
         close(nh%lun)
         open(newunit=nh%lun, file=nh%par_file)
         nh%errstr=""
         read(unit=nh%lun, nml=COSMIC_RAYS, iostat=nh%ierrh, iomsg=nh%errstr)
         close(nh%lun)
         call nh%namelist_errh(nh%ierrh, "COSMIC_RAYS")
         read(nh%cmdl_nml,nml=COSMIC_RAYS, iostat=nh%ierrh)
         call nh%namelist_errh(nh%ierrh, "COSMIC_RAYS", .true.)
         open(newunit=nh%lun, file=nh%tmp2, status="unknown")
         write(nh%lun,nml=COSMIC_RAYS)
         close(nh%lun)
         call nh%compare_namelist()
      endif

      rbuff(:) = huge(1.)                         ! mark unused entries to allow automatic determination of nn

      if (master) then

         cbuff(1) = divv_scheme

         ibuff(1) = ncr_user
         ibuff(2) = ncrb

         rbuff(1) = cfl_cr
         rbuff(2) = smallecr
         rbuff(3) = cr_active
         rbuff(4) = cr_eff
         rbuff(5) = gamma_cr

         lbuff(1) = use_CRdiff
         lbuff(2) = use_CRdecay
         lbuff(3) = use_smallecr

         ncrsp    = ncrsp + ncr_user
         nl       = 3                                     ! this must match the last lbuff() index above
         nn       = count(rbuff(:) < huge(1.), kind=4)    ! this must match the last rbuff() index above
         ibuff(ubound(ibuff, 1)    ) = nn
         ibuff(ubound(ibuff, 1) - 1) = nl

         if (nn + 2 * ncrsp > ubound(rbuff, 1)) call die("[initcosmicrays:init_cosmicrays] rbuff size exceeded.")
         if (nl + ncr_user  > ubound(lbuff, 1)) call die("[initcosmicrays:init_cosmicrays] lbuff size exceeded.")

         if (ncrsp > 0) then
            rbuff(nn+1      :nn+  ncrsp) = K_cr_paral(1:ncrsp)
            rbuff(nn+1+ncrsp:nn+2*ncrsp) = K_cr_perp (1:ncrsp)

            lbuff(nl+1:nl+ncr_user) = gpcr_ess_user(1:ncr_user)
         endif


      endif

      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(rbuff)
      call piernik_MPI_Bcast(lbuff)
      call piernik_MPI_Bcast(cbuff, cbuff_len)

      if (slave) then

         divv_scheme  = cbuff(1)

         ncr_user     = int(ibuff(1), kind=4)
         ncrb         = int(ibuff(2), kind=4)

         cfl_cr       = rbuff(1)
         smallecr     = rbuff(2)
         cr_active    = rbuff(3)
         cr_eff       = rbuff(4)
         gamma_cr     = rbuff(5)

         use_CRdiff   = lbuff(1)
         use_CRdecay  = lbuff(2)
         use_smallecr = lbuff(3)

         ncrsp        = ncrsp + ncr_user
         nn           = ibuff(ubound(ibuff, 1)    )    ! this must match the last rbuff() index above
         nl           = ibuff(ubound(ibuff, 1) - 1)    ! this must match the last lbuff() index above

         if (ncrsp > 0) then
            K_cr_paral(1:ncrsp) = rbuff(nn+1      :nn+  ncrsp)
            K_cr_perp (1:ncrsp) = rbuff(nn+1+ncrsp:nn+2*ncrsp)

            gpcr_ess_user(1:ncr_user) = lbuff(nl+1:nl+ncr_user)
         endif

      endif

      gamma_cr_1 = gamma_cr - 1.0

      call cr_species_tables(ncrsp, gpcr_ess_user(1:ncr_user))

      nspc = count(cr_spectral, kind=4)
      ncrn = ncrsp - nspc

      ncr2b  = I_TWO * ncrb
      ncrtot = ncr2b * nspc + ncrn

      if (any([ncrsp, ncrb] > ncr_max) .or. any([ncrsp, ncrb] < 0)) call die("[initcosmicrays:init_cosmicrays] ncr[nes] > ncr_max or ncr[nes] < 0")
      if (ncrtot == 0) call warn("[initcosmicrays:init_cosmicrays] ncrtot == 0; no cr components specified")

      ma1d = [ncrtot]
      call my_allocate(K_crs_paral, ma1d)
      call my_allocate(K_crs_perp,  ma1d)

      K_crs_paral(:) = 0.0
      K_crs_perp (:) = 0.0

      if (ncrsp > 0) then
         K_crs_paral(1:ncrn) = pack(K_cr_paral(1:ncrsp), .not.cr_spectral)
         K_crs_perp (1:ncrn) = pack(K_cr_perp (1:ncrsp), .not.cr_spectral)
      endif

      ma1d = [ncrn]
      call my_allocate(iarr_crn, ma1d)

      if (ncrb <= 0) then
         ma1d = 0
      else
         ma1d = [ncr2b * nspc]
      endif
      call my_allocate(iarr_crspc, ma1d) ! < iarr_crspc will point: (1:ncrb) - cre number per bin, (ncrb+1:2*ncrb) - cre energy per bin

#ifdef CRESP
      ma1d = [ncrb * nspc]
      call my_allocate(iarr_crspc_e, ma1d)
      call my_allocate(iarr_crspc_n, ma1d)
      ma2d = [nspc, ncrb]
      call my_allocate(iarr_crspc2_e, ma2d)
      call my_allocate(iarr_crspc2_n, ma2d)
#endif /* CRESP */
      ma1d = [ncrtot]
      call my_allocate(iarr_crs, ma1d)

      ma1d = [ int(count(cr_gpess .and. .not.cr_spectral), kind=4) ]
      call my_allocate(gpcr_ess_noncresp, ma1d)
      gpcr_ess_noncresp = pack([(icr, icr = I_ONE, count(.not.cr_spectral))], mask=(pack(cr_gpess, mask=(.not.cr_spectral))))

      def_dtcrs = big
      maxKcrs = maxval(K_cr_paral(1:ncrsp) + K_cr_perp(1:ncrsp), mask=.not.cr_spectral)
      K_crs_valid = (maxKcrs > 0)
      if (maxKcrs .notequals. 0.) def_dtcrs = cfl_cr * half / maxKcrs

   end subroutine init_cosmicrays

   subroutine cosmicray_index(flind)

      use constants,  only: I_ONE
      use fluidtypes, only: var_numbers

      implicit none

      type(var_numbers), intent(inout) :: flind
      integer(kind=4)                  :: icr, jnb

      flind%crn%beg = flind%all + I_ONE
      flind%crs%beg = flind%crn%beg

      flind%crn%all   = ncrn
      flind%crspc%all = ncr2b * nspc

      flind%crs%all = flind%crn%all + flind%crspc%all
      do icr = 1, ncrn
         iarr_crn(icr) = flind%all + icr
         iarr_crs(icr) = flind%all + icr
      enddo
      flind%all = flind%all + flind%crn%all

      do icr = I_ONE, ncr2b * nspc
         iarr_crspc(icr)      = flind%all + icr
         iarr_crs(ncrn + icr) = flind%all + icr
      enddo

      flind%all = flind%all + flind%crspc%all

      flind%crn%end   = flind%crn%beg + flind%crn%all - I_ONE
      flind%crspc%beg = flind%crn%end + I_ONE
      flind%crspc%end = flind%all
      flind%crs%end = flind%crspc%end
      if (flind%crn%all  /= 0) flind%components = flind%components + I_ONE
      flind%crn%pos = flind%components
      if (flind%crspc%all  /= 0) flind%components = flind%components + I_ONE
      flind%crspc%pos = flind%components

#ifdef CRESP
      flind%crspc%nbeg = flind%crn%end + I_ONE
      flind%crspc%nend = flind%crn%end + ncrb * nspc
      flind%crspc%ebeg = flind%crspc%nend + I_ONE
      flind%crspc%eend = flind%crspc%nend + ncrb * nspc

      do icr = 1, ncrb
         iarr_crspc_n(icr) = flind%crspc%nbeg - I_ONE + icr
         iarr_crspc_e(icr) = flind%crspc%ebeg - I_ONE + icr
      enddo

      do icr = 1, nspc        !< Arrange iterable indexes for each spectral species separately; first indexes for n, then e.
         do jnb = 1, ncrb
            iarr_crspc2_n(icr, jnb) = flind%crn%end + I_ONE + (icr - I_ONE) * ncrb + jnb
            iarr_crspc2_e(icr, jnb) = flind%crn%end + I_ONE + nspc * ncrb + (icr - I_ONE) * ncrb + jnb
         enddo
      enddo
#endif /* CRESP */

   end subroutine cosmicray_index

   subroutine cleanup_cosmicrays

      use diagnostics, only: my_deallocate

      implicit none

      call my_deallocate(iarr_crn)
      call my_deallocate(iarr_crspc)
      call my_deallocate(iarr_crs)
      call my_deallocate(K_crs_paral)
      call my_deallocate(K_crs_perp)
      call my_deallocate(gpcr_ess_noncresp)
#ifdef CRESP
      call my_deallocate(iarr_crspc_e)
      call my_deallocate(iarr_crspc_n)
#endif /* CRESP */

   end subroutine cleanup_cosmicrays

end module initcosmicrays
