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
#ifdef COSM_RAY_ELECTRONS
   use initcrspectrum, only: ncre, K_cre_paral_1, K_cre_perp_1, K_cre_pow, p_mid_fix, cre_gpcr_ess
#endif /* COSM_RAY_ELECTRONS */
   implicit none

   public ! QA_WARN no secrets are kept here
   private :: cbuff_len ! QA_WARN prevent reexport

   integer, parameter                  :: ncr_max = 102  !< maximum number of CR nuclear and electron components (\warning higher ncr_max limit would require changes in names of components in common_hdf5)
   ! namelist parameters
   integer(kind=4)                     :: ncrn         !< number of CR nuclear  components \deprecated BEWARE: ncrs (sum of ncrn and ncre) should not be higher than ncr_max = 9
#ifndef COSM_RAY_ELECTRONS
   integer(kind=4)                     :: ncre         !< number of CR electron components \deprecated BEWARE: ncrs (sum of ncrn and ncre) should not be higher than ncr_max = 9
#endif /* !COSM_RAY_ELECTRONS */
   integer(kind=4)                     :: ncrs         !< number of all CR components \deprecated BEWARE: ncrs (sum of ncrn and ncre) should not be higher than ncr_max = 9
   real                                :: cfl_cr       !< CFL number for diffusive CR transport
   real                                :: smallecr     !< floor value for CR energy density
   real                                :: cr_active    !< parameter specifying whether CR pressure gradient is (when =1.) or isn't (when =0.) included in the gas equation of motion
   real                                :: cr_eff       !< conversion rate of SN explosion energy to CR energy (default = 0.1)
   logical                             :: use_split    !< apply all diffusion operators at once (.false.) or use directional splittiong (.true.)
   real, dimension(ncr_max)            :: gamma_crn    !< array containing adiabatic indexes of all CR nuclear components
   real, dimension(ncr_max)            :: K_crn_paral  !< array containing parallel diffusion coefficients of all CR nuclear components
   real, dimension(ncr_max)            :: K_crn_perp   !< array containing perpendicular diffusion coefficients of all CR nuclear components
#ifndef COSM_RAY_ELECTRONS
   real, dimension(ncr_max)            :: gamma_cre    !< array containing adiabatic indexes of all CR electron components
#endif /* !COSM_RAY_ELECTRONS */
   real, dimension(ncr_max)            :: K_cre_paral  !< array containing parallel diffusion coefficients of all CR electron components (number density and energy density)
   real, dimension(ncr_max)            :: K_cre_perp   !< array containing perpendicular diffusion coefficients of all CR electron components (number density and energy density)
   character(len=cbuff_len)            :: divv_scheme  !< scheme used to calculate div(v), see crhelpers for more details
   logical, dimension(ncr_max)         :: crn_gpcr_ess !< if CRn species/energy-bin is essential for grad_pcr calculation
#ifndef COSM_RAY_ELECTRONS
   logical, dimension(ncr_max)         :: cre_gpcr_ess !< if CRe species/energy-bin is essential for grad_pcr calculation
#endif /* !COSM_RAY_ELECTRONS */
   integer(kind=4), allocatable, dimension(:) :: gpcr_essential !< crs indexes of essentials for grad_pcr calculation
   ! public component data
   integer(kind=4), allocatable, dimension(:) :: iarr_crn !< array of indexes pointing to all CR nuclear components
   integer(kind=4), allocatable, dimension(:) :: iarr_cre !< array of indexes pointing to all CR electron components
   integer(kind=4), allocatable, dimension(:) :: iarr_crs !< array of indexes pointing to all CR components
#ifdef COSM_RAY_ELECTRONS
   integer(kind=4), allocatable, dimension(:) :: iarr_cre_e !< array of indexes pointing to all CR electron energy components
   integer(kind=4), allocatable, dimension(:) :: iarr_cre_n !< array of indexes pointing to all CR electron number density components
#endif /* COSM_RAY_ELECTRONS */
   integer(kind=4), allocatable, dimension(:) :: iarr_crs_diff !< array of indexes pointing to all conventionally diffusing CR components

   real,    allocatable, dimension(:)  :: gamma_crs    ! < array containing adiabatic indexes of all CR components

   real,    allocatable, dimension(:)  :: K_crs_paral  !< array containing parallel diffusion coefficients of all CR components
   real,    allocatable, dimension(:)  :: K_crs_perp   !< array containing perpendicular diffusion coefficients of all CR components
   !> \deprecated BEWARE Possible confusion: *_perp coefficients are not "perpendicular" but rather isotropic

contains

!>
!! \brief Routine to set parameters values from namelist COSMIC_RAYS
!!
!! \n \n
!! @b COSMIC_RAYS
!! \n \n
!! <table border="+1">
!! <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!! <tr><td>cfl_cr      </td><td>0.9   </td><td>real value</td><td>\copydoc initcosmicrays::cfl_cr     </td></tr>
!! <tr><td>smallecr    </td><td>0.0   </td><td>real value</td><td>\copydoc initcosmicrays::smallecr   </td></tr>
!! <tr><td>cr_active   </td><td>1.0   </td><td>real value</td><td>\copydoc initcosmicrays::cr_active  </td></tr>
!! <tr><td>cr_eff      </td><td>0.1   </td><td>real value</td><td>\copydoc initcosmicrays::cr_eff     </td></tr>
!! <tr><td>use_split   </td><td>.true.</td><td>logical   </td><td>\copydoc initcosmicrays::use_split  </td></tr>
!! <tr><td>ncrn        </td><td>0     </td><td>integer   </td><td>\copydoc initcosmicrays::ncrn       </td></tr>
!! <tr><td>ncre        </td><td>0     </td><td>integer   </td><td>\copydoc initcosmicrays::ncre       </td></tr>
!! <tr><td>gamma_crn   </td><td>4./3. </td><td>real array</td><td>\copydoc initcosmicrays::gamma_crn  </td></tr>
!! <tr><td>gamma_cre   </td><td>4./3. </td><td>real array</td><td>\copydoc initcosmicrays::gamma_cre  </td></tr>
!! <tr><td>K_crn_paral </td><td>0     </td><td>real array</td><td>\copydoc initcosmicrays::k_crn_paral</td></tr>
!! <tr><td>K_crn_perp  </td><td>0     </td><td>real array</td><td>\copydoc initcosmicrays::k_crn_perp </td></tr>
!! <tr><td>K_cre_paral </td><td>0     </td><td>real array</td><td>\copydoc initcosmicrays::k_cre_paral</td></tr>
!! <tr><td>K_cre_perp  </td><td>0     </td><td>real array</td><td>\copydoc initcosmicrays::k_cre_perp </td></tr>
!! <tr><td>divv_scheme </td><td>''    </td><td>string    </td><td>\copydoc initcosmicrays::divv_scheme</td></tr>
!! <tr><td>crn_gpcr_ess</td><td>(1): .true.; (>2):.false.</td><td>logical</td><td>\copydoc initcosmicrays::crn_gpcr_ess</td></tr>
!! <tr><td>cre_gpcr_ess</td><td>.false.                  </td><td>logical</td><td>\copydoc initcosmicrays::cre_gpcr_ess</td></tr>
!! </table>
!! The list is active while \b "COSM_RAYS" is defined.
!! \n \n
!<
   subroutine init_cosmicrays

      use constants,       only: cbuff_len, I_ONE
      use diagnostics,     only: ma1d, my_allocate
      use dataio_pub,      only: nh   ! QA_WARN required for diff_nml
      use dataio_pub,      only: die, warn
      use mpisetup,        only: ibuff, rbuff, lbuff, cbuff, master, slave, piernik_MPI_Bcast
#ifdef COSM_RAYS_SOURCES
      use cr_data,         only: init_crsources
#endif /* COSM_RAYS_SOURCES */
      implicit none

      integer(kind=4) :: nn, icr, jcr
      integer         :: ne
#ifdef COSM_RAY_ELECTRONS
      namelist /COSMIC_RAYS/ cfl_cr, smallecr, cr_active, cr_eff, use_split, &
           &                 ncrn, gamma_crn, K_crn_paral, K_crn_perp, &
           &                 divv_scheme, crn_gpcr_ess
#else /* !COSM_RAY_ELECTRONS */
      namelist /COSMIC_RAYS/ cfl_cr, smallecr, cr_active, cr_eff, use_split, &
           &                 ncrn, gamma_crn, K_crn_paral, K_crn_perp, &
           &                 ncre, gamma_cre, K_cre_paral, K_cre_perp, &
           &                 divv_scheme, crn_gpcr_ess, cre_gpcr_ess
#endif /* COSM_RAY_ELECTRONS */

      cfl_cr     = 0.9
      smallecr   = 0.0
      cr_active  = 1.0
      cr_eff     = 0.1       !  canonical conversion rate of SN en.-> CR
      !  we fix E_SN=10**51 erg
      ncrn       = 0
#ifndef COSM_RAY_ELECTRONS
      ncre       = 0
      gamma_cre(:)   = 4./3.
      K_cre_paral(:) = 0.0
      K_cre_perp(:)  = 0.0

      cre_gpcr_ess(:) = .false.
#endif /* !COSM_RAY_ELECTRONS */
      use_split  = .true.

      gamma_crn(:)   = 4./3.
      K_crn_paral(:) = 0.0
      K_crn_perp(:)  = 0.0
      K_cre_paral(:) = 0.0
      K_cre_perp(:)  = 0.0

      crn_gpcr_ess(:) = .false.
      crn_gpcr_ess(1) = .true.       ! in most cases protons are the first ingredient of CRs and they are essential

      divv_scheme = ''

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

#ifndef MULTIGRID
      if (.not. use_split) call warn("[initcosmicrays:init_cosmicrays] No multigrid solver compiled in: use_split reset to .true.")
      use_split  = .true.
#endif /* !MULTIGRID */

      rbuff(:)   = huge(1.)                         ! mark unused entries to allow automatic determination of nn

      if (master) then

         ibuff(1)   = ncrn
#ifndef COSM_RAY_ELECTRONS
         ibuff(2)   = ncre
#endif /* COSM_RAY_ELECTRONS */

         rbuff(1)   = cfl_cr
         rbuff(2)   = smallecr
         rbuff(3)   = cr_active
         rbuff(4)   = cr_eff

         lbuff(1)   = use_split

         cbuff(1)   = divv_scheme

         nn         = count(rbuff(:) < huge(1.), kind=4)    ! this must match the last rbuff() index above
         ibuff(ubound(ibuff, 1)) = nn
         ne         = nn + 3 * ncrn
         if (ne + 3 * ncre > ubound(rbuff, 1)) call die("[initcosmicrays:init_cosmicrays] rbuff size exceeded.")

         if (ncrn > 0) then
            rbuff(nn+1       :nn+  ncrn) = gamma_crn  (1:ncrn)
            rbuff(nn+1+  ncrn:nn+2*ncrn) = K_crn_paral(1:ncrn)
            rbuff(nn+1+2*ncrn:nn+3*ncrn) = K_crn_perp (1:ncrn)

            lbuff(2:ncrn+1) = crn_gpcr_ess(1:ncrn)
         endif

         if (ncre > 0) then
#ifndef COSM_RAY_ELECTRONS
            rbuff(ne+1       :ne+ncre) = K_cre_paral(1:ncre)  !< K_cre_paral explicitly defined & broadcasted if CRESP not in use
            rbuff(ne+1+ncre:ne+2*ncre) = K_cre_perp (1:ncre)  !< K_cre_perp explicitly defined & broadcasted if CRESP not in use
            rbuff(ne+2*ncre+1:ne+3*ncre) = gamma_cre(1:ncre)  ! gamma_cre used only if CRESP module not used
            lbuff(ncrn+2:ncrn+1+ncre)    = cre_gpcr_ess(1:ncre)
#else
            lbuff(ncrn+2)                = cre_gpcr_ess
#endif /* !COSM_RAY_ELECTRONS */

         endif

      endif

      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(rbuff)
      call piernik_MPI_Bcast(lbuff)
      call piernik_MPI_Bcast(cbuff, cbuff_len)

      if (slave) then

         ncrn       = int(ibuff(1), kind=4)
#ifndef COSM_RAY_ELECTRONS
         ncre       = int(ibuff(2), kind=4)
#endif /* COSM_RAY_ELECTRONS */

         cfl_cr     = rbuff(1)
         smallecr   = rbuff(2)
         cr_active  = rbuff(3)
         cr_eff     = rbuff(4)

         use_split  = lbuff(1)

         nn         = ibuff(ubound(ibuff, 1))    ! this must match the last rbuff() index above
         ne         = nn + 3 * ncrn

         divv_scheme = cbuff(1)

         if (ncrn > 0) then
            gamma_crn  (1:ncrn) = rbuff(nn+1       :nn+  ncrn)
            K_crn_paral(1:ncrn) = rbuff(nn+1+  ncrn:nn+2*ncrn)
            K_crn_perp (1:ncrn) = rbuff(nn+1+2*ncrn:nn+3*ncrn)
            crn_gpcr_ess(1:ncrn) = lbuff(2:ncrn+1)
         endif

         if (ncre > 0) then
#ifndef COSM_RAY_ELECTRONS
            K_cre_paral(1:ncre) = rbuff(ne+1       :ne+ncre)    !< K_cre_paral explicitly defined & broadcasted if CRESP not in use
            K_cre_perp (1:ncre) = rbuff(ne+1+ncre:ne+2*ncre)    !< K_cre_perp explicitly defined & broadcasted if CRESP not in use
            gamma_cre(1:ncre)    = rbuff(ne+2*ncre+1:ne+3*ncre) ! gamma_cre used only if CRESP module not used
            cre_gpcr_ess(1:ncre) = lbuff(ncrn+2:ncrn+2+ncre)
#else
            cre_gpcr_ess         = lbuff(ncrn+2)
#endif /* !COSM_RAY_ELECTRONS */
         endif

      endif

      ncrs = ncrn + ncre
#ifdef COSM_RAY_ELECTRONS
      ncrs = ncrs + ncre ! ncrn + 2 * ncre overall
#endif /* COSM_RAY_ELECTRONS */

      if (any([ncrn, ncre] > ncr_max) .or. any([ncrn, ncre] < 0)) call die("[initcosmicrays:init_cosmicrays] ncr[nes] > ncr_max or ncr[nes] < 0")
      if (ncrs ==0) call warn("[initcosmicrays:init_cosmicrays] ncrs == 0; no cr components specified")

      ma1d = [ncrs]
      call my_allocate(gamma_crs,   ma1d)
      call my_allocate(K_crs_paral, ma1d)
      call my_allocate(K_crs_perp,  ma1d)

      gamma_crs  (:) = 0.0
      K_crs_paral(:) = 0.0
      K_crs_perp (:) = 0.0

      if (ncrn > 0) then
        gamma_crs  (1:ncrn) = gamma_crn  (1:ncrn)
        K_crs_paral(1:ncrn) = K_crn_paral(1:ncrn)
        K_crs_perp (1:ncrn) = K_crn_perp (1:ncrn)
      endif

#ifdef COSM_RAY_ELECTRONS
      if (ncre > 0) then
         K_crs_paral(ncrn+1:ncrn+ncre) = K_cre_paral_1  * (p_mid_fix**K_cre_pow)/(maxval(p_mid_fix)**K_cre_pow)          !< CRESP number density K
         K_crs_paral(ncrn+1+ncre:ncrn+2*ncre) = K_cre_paral_1 * (p_mid_fix**K_cre_pow)/(maxval(p_mid_fix)**K_cre_pow)    !< CRESP energy density K
         K_crs_perp(ncrn+1:ncrn+ncre)  = K_cre_perp_1 * (p_mid_fix**K_cre_pow)/(maxval(p_mid_fix)**K_cre_pow)             !< CRESP number density K
         K_crs_perp(ncrn+ncre+1:ncrn+2*ncre) = K_cre_perp_1 * (p_mid_fix**K_cre_pow)/(maxval(p_mid_fix)**K_cre_pow)      !< CRESP energy density K

         K_cre_paral(1:ncre)        =  K_crs_paral(ncrn+1:ncrn+ncre)               !< CRESP number density K
         K_cre_paral(ncre+1:2*ncre) =  K_crs_paral(ncrn+1+ncre:ncrn+2*ncre)        !< CRESP energy density K
         K_cre_perp(1:ncre)         =  K_crs_perp(ncrn+1:ncrn+ncre)                !< CRESP number density K
         K_cre_perp(ncre+1:2*ncre)  =  K_crs_perp(ncrn+1+ncre:ncrn+2*ncre)         !< CRESP energy density K
      endif
#endif /* COSM_RAY_ELECTRONS */
      ma1d = [ncrn]
      call my_allocate(iarr_crn, ma1d)

      if (ncre .le. 0) then
         ma1d = 0
      else
#ifdef COSM_RAY_ELECTRONS
         ma1d = [2*ncre]
#else
         ma1d = ncre
#endif /*COSM_RAY_ELECTRONS */
      endif
      call my_allocate(iarr_cre, ma1d) ! < iarr_cre will point: (1:ncre) - cre number per bin, (ncre+1:2*ncre) - cre energy per bin

#ifdef COSM_RAY_ELECTRONS
      ma1d = [ncre]
      call my_allocate(iarr_cre_e, ma1d)
      call my_allocate(iarr_cre_n, ma1d)
#endif /* COSM_RAY_ELECTRONS */
      ma1d = [ncrs]
      call my_allocate(iarr_crs, ma1d)

      ma1d = [ncrs]
      call my_allocate(iarr_crs_diff, ma1d)

#ifdef COSM_RAYS_SOURCES
      call init_crsources(ncrn, crn_gpcr_ess)
#endif /* COSM_RAYS_SOURCES */

      ma1d = [ int(count(crn_gpcr_ess), kind=4) ]
#ifndef COSM_RAY_ELECTRONS
      ma1d = [ int(count(crn_gpcr_ess) + int(count(cre_gpcr_ess)), kind=4) ]
#endif /* !COSM_RAY_ELECTRONS */

      call my_allocate(gpcr_essential, ma1d)
      jcr = 0
      do icr = 1, ncrn
         if (crn_gpcr_ess(icr)) then
            jcr = jcr + I_ONE
            gpcr_essential(jcr) = icr
         endif
      enddo
#ifndef COSM_RAY_ELECTRONS
      do icr = 1, ncre
        if (cre_gpcr_ess(icr)) then
            jcr = jcr + I_ONE
            gpcr_essential(jcr) = ncrn + 1
        endif
      enddo
#endif /* !COSM_RAY_ELECTRONS */
   end subroutine init_cosmicrays

   subroutine cosmicray_index(flind)

      use constants,    only: I_ONE
      use fluidtypes,   only: var_numbers

      implicit none

      type(var_numbers), intent(inout) :: flind
      integer(kind=4) :: icr, ncr

      flind%crn%beg    = flind%all + I_ONE
      flind%crs%beg    = flind%crn%beg

      flind%crn%all  = ncrn

      if (ncre .le. 0) then
            flind%cre%all  = 0
      else
#ifdef COSM_RAY_ELECTRONS
            flind%cre%all  = 2*ncre
#else
            flind%cre%all  = ncre
#endif /* COSM_RAY_ELECTRONS */
      endif

      flind%crs%all  = flind%crn%all + flind%cre%all
      do icr = 1, ncrn
         iarr_crn(icr)      = flind%all + icr
         iarr_crs(icr)      = flind%all + icr
      enddo
      flind%all = flind%all + flind%crn%all

      if (ncre.gt.0) then
       ncr = 0
#ifdef COSM_RAY_ELECTRONS
       ncr = 2 * ncre
#else
       ncr = ncre
#endif /* COSM_RAY_ELECTRONS */
        do icr = 1, ncr
          iarr_cre(icr)        = flind%all + icr
          iarr_crs(ncrn + icr) = flind%all + icr
       enddo
      endif

      flind%all = flind%all + flind%cre%all
      flind%crn%end = flind%crn%beg + flind%crn%all - I_ONE
      flind%cre%beg = flind%crn%end + I_ONE
      flind%cre%end = flind%all
      flind%crs%end = flind%cre%end
      if (flind%crn%all  /= 0) flind%components = flind%components + I_ONE
      flind%crn%pos = flind%components
      if (flind%cre%all  /= 0) flind%components = flind%components + I_ONE
      flind%cre%pos = flind%components

#ifdef COSM_RAY_ELECTRONS
     flind%cre%nbeg = flind%crn%end + I_ONE
     flind%cre%nend = flind%crn%end + ncre
     flind%cre%ebeg = flind%cre%nend + I_ONE
     flind%cre%eend = flind%cre%nend + ncre

     do icr = 1, ncre
        iarr_cre_n(icr) = flind%cre%nbeg - I_ONE + icr
        iarr_cre_e(icr) = flind%cre%ebeg - I_ONE + icr
     enddo
#endif /* COSM_RAY_ELECTRONS */
!      if ( ncre.eq.0) then
         iarr_crs_diff = iarr_crs
!      endif
   end subroutine cosmicray_index

   subroutine cleanup_cosmicrays
      use diagnostics, only: my_deallocate
      implicit none

      call my_deallocate(iarr_crn)
      call my_deallocate(iarr_cre)
      call my_deallocate(iarr_crs)
      call my_deallocate(gamma_crs)
      call my_deallocate(K_crs_paral)
      call my_deallocate(K_crs_perp)
      call my_deallocate(gpcr_essential)

   end subroutine cleanup_cosmicrays

end module initcosmicrays
