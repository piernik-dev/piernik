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
! pulled by COSM_RAY_ELECTRONS
   use cresp_vars
   use cresp_spectrum, only: cresp_deallocate_active_arrays, cresp_init_state, cresp_update_bin_index
   use constants, only: cbuff_len
   use units,     only: clight
   implicit none

   public ! QA_WARN no secrets are kept here
   private :: cbuff_len ! QA_WARN prevent reexport

   integer, parameter                  :: ncr_max = 9  !< maximum number of CR nuclear and electron components (\warning higher ncr_max limit would require changes in names of components in common_hdf5)

   ! namelist parameters
   integer(kind=4)                     :: ncrn         !< number of CR nuclear  components \deprecated BEWARE: ncrs (sum of ncrn and ncre) should not be higher than ncr_max = 9
   integer(kind=4)                     :: ncre         !< number of CR electron components \deprecated BEWARE: ncrs (sum of ncrn and ncre) should not be higher than ncr_max = 9
   integer(kind=4)                     :: ncrs         !< number of all CR components \deprecated BEWARE: ncrs (sum of ncrn and ncre) should not be higher than ncr_max = 9
   real                                :: cfl_cr       !< CFL number for diffusive CR transport
   real                                :: smallecr     !< floor value for CR energy density
   real                                :: cr_active    !< parameter specifying whether CR pressure gradient is (when =1.) or isn't (when =0.) included in the gas equation of motion
   real                                :: cr_eff       !< conversion rate of SN explosion energy to CR energy (default = 0.1)
   logical                             :: use_split    !< apply all diffusion operators at once (.false.) or use directional splittiong (.true.)
   real, dimension(ncr_max)            :: gamma_crn    !< array containing adiabatic indexes of all CR nuclear components
   real, dimension(ncr_max)            :: K_crn_paral  !< array containing parallel diffusion coefficients of all CR nuclear components
   real, dimension(ncr_max)            :: K_crn_perp   !< array containing perpendicular diffusion coefficients of all CR nuclear components
   real, dimension(ncr_max)            :: gamma_cre    !< array containing adiabatic indexes of all CR nuclear components
   real, dimension(ncr_max)            :: K_cre_paral  !< array containing parallel diffusion coefficients of all CR electron components
   real, dimension(ncr_max)            :: K_cre_perp   !< array containing perpendicular diffusion coefficients of all CR electron components
   character(len=cbuff_len)            :: divv_scheme  !< scheme used to calculate div(v), see crhelpers for more details
   logical, dimension(ncr_max)         :: crn_gpcr_ess !< if CRn species/energy-bin is essential for grad_pcr calculation
   logical, dimension(ncr_max)         :: cre_gpcr_ess !< if CRe species/energy-bin is essential for grad_pcr calculation
   integer(kind=4), allocatable, dimension(:) :: gpcr_essential !< crs indexes of essentials for grad_pcr calculation
   
   ! public component data
   integer(kind=4), allocatable, dimension(:) :: iarr_crn !< array of indexes pointing to all CR nuclear components
   integer(kind=4), allocatable, dimension(:) :: iarr_cre !< array of indexes pointing to all CR electron components
   integer(kind=4), allocatable, dimension(:) :: iarr_crs !< array of indexes pointing to all CR components

   real,    allocatable, dimension(:)  :: gamma_crs    !< array containing adiabatic indexes of all CR components
!    real,    allocatable, dimension(:)  :: K_crs_paral  !< array containing parallel diffusion coefficients of all CR components
!    real,    allocatable, dimension(:)  :: K_crs_perp   !< array containing perpendicular diffusion coefficients of all CR components
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

      namelist /COSMIC_RAYS/ cfl_cr, smallecr, cr_active, cr_eff, use_split, &
           &                 ncrn, gamma_crn, K_crn_paral, K_crn_perp, &
           &                 ncre, gamma_cre, K_cre_paral, K_cre_perp, &
           &                 divv_scheme, crn_gpcr_ess, cre_gpcr_ess

      cfl_cr     = 0.9
      smallecr   = 0.0
      cr_active  = 1.0
      cr_eff     = 0.1       !  canonical conversion rate of SN en.-> CR
      !  we fix E_SN=10**51 erg
      ncrn       = 0
      ncre       = 0

      use_split  = .true.

      gamma_crn(:)   = 4./3.
      K_crn_paral(:) = 0.0
      K_crn_perp(:)  = 0.0
      gamma_cre(:)   = 4./3.
      K_cre_paral(:) = 0.0
      K_cre_perp(:)  = 0.0

      crn_gpcr_ess(:) = .false.
      crn_gpcr_ess(1) = .true.       ! in most cases protons are the first ingredient of CRs and they are essential
      cre_gpcr_ess(:) = .false.

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
         call nh%compare_namelist() ! Do not use one-line if here!

      endif

#ifndef MULTIGRID
      if (.not. use_split) call warn("[initcosmicrays:init_cosmicrays] No multigrid solver compiled in: use_split reset to .true.")
      use_split  = .true.
#endif /* !MULTIGRID */

      rbuff(:)   = huge(1.)                         ! mark unused entries to allow automatic determination of nn

      if (master) then

         ibuff(1)   = ncrn
         ibuff(2)   = ncre

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
            rbuff(ne+1       :ne+  ncre) = gamma_cre  (1:ncre)
            rbuff(ne+1+  ncre:ne+2*ncre) = K_cre_paral(1:ncre)
            rbuff(ne+1+2*ncre:ne+3*ncre) = K_cre_perp (1:ncre)
            lbuff(ncrn+2:ncrn+ncre+1) = cre_gpcr_ess(1:ncre)
         endif

      endif

      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(rbuff)
      call piernik_MPI_Bcast(lbuff)
      call piernik_MPI_Bcast(cbuff, cbuff_len)

      if (slave) then

         ncrn       = int(ibuff(1), kind=4)
         ncre       = int(ibuff(2), kind=4)

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
            gamma_cre  (1:ncre) = rbuff(ne+1       :ne+  ncre)
            K_cre_paral(1:ncre) = rbuff(ne+1+  ncre:ne+2*ncre)
            K_cre_perp (1:ncre) = rbuff(ne+1+2*ncre:ne+3*ncre)
            cre_gpcr_ess(1:ncre) = lbuff(ncrn+2:ncrn+ncre+1)
         endif

      endif
      
          
!      ncrs = ncre + ncrn
            ncrs = (2*ncre+2) + ncrn   !!!!!

      if (any([ncrn, ncre] > ncr_max) .or. any([ncrn, ncre] < 0)) call die("[initcosmicrays:init_cosmicrays] ncr[nes] > ncr_max or ncr[nes] < 0")
      if (ncrs ==0) call warn("[initcosmicrays:init_cosmicrays] ncrs == 0; no cr components specified")

!       ma1d = [ncrs]
!       call my_allocate(gamma_crs,   ma1d)
!       call my_allocate(K_crs_paral, ma1d)
!       call my_allocate(K_crs_perp,  ma1d)

!        if (ncrn > 0) then
!           gamma_crs  (1:ncrn) = gamma_crn  (1:ncrn)
!           K_crs_paral(1:ncrn) = K_crn_paral(1:ncrn)
!           K_crs_perp (1:ncrn) = K_crn_perp (1:ncrn)
!        endif

       if (ncre > 0) then
            gamma_cre  (1:ncre) = 0 !rbuff(ne+1       :ne+  ncre)
            K_cre_paral(1:ncre) = 0 !rbuff(ne+1+  ncre:ne+2*ncre)
            K_cre_perp (1:ncre) = 0 !rbuff(ne+1+2*ncre:ne+3*ncre)
       endif



!       if (ncre > 0) then
!          gamma_crs  (ncrn+1:ncrs) = 0 !gamma_cre  (1:ncre) !<- cre gamma and K is supposed to be 0
!          K_crs_paral(ncrn+1:ncrs) = 0!K_cre_paral(1:ncre)
!          K_crs_perp (ncrn+1:ncrs) = 0!K_cre_perp (1:ncre)
!          allocate(cres_n(ncre))   !!!
!          allocate(cres_en(ncre))  !!!
!       endif
     
      ma1d = [ncrn]
      call my_allocate(iarr_crn, ma1d)
      ma1d = [2*ncre+2] !!!
      call my_allocate(iarr_cre, ma1d)
      ma1d = [ncrs]
      call my_allocate(iarr_crs, ma1d)

#ifdef COSM_RAYS_SOURCES
      call init_crsources(ncrn, crn_gpcr_ess)
#endif /* COSM_RAYS_SOURCES */

      ma1d = [ int(count(crn_gpcr_ess) + count(cre_gpcr_ess), kind=4) ]
      call my_allocate(gpcr_essential, ma1d)
      jcr = 0
      do icr = 1, ncrn
         if (crn_gpcr_ess(icr)) then
            jcr = jcr + I_ONE
            gpcr_essential(jcr) = icr
         endif
      enddo
      do icr = 1, ncre
         if (cre_gpcr_ess(icr)) then
            jcr = jcr + I_ONE
            gpcr_essential(jcr) = icr + ncrn
         endif
      enddo

   end subroutine init_cosmicrays
   
      subroutine cresp_init_state(dt)
      implicit none
      real(kind = 8)    :: dt
      integer           :: i, k
      real(kind=8)      ::  c ! width of bin

      all_edges = (/ (i,i=0,ncre) /)
      all_bins = (/ (i,i=1,ncre) /)

      crel%q = q_init

      w  = (log10(p_max_fix/p_min_fix))/dble(ncre-2)
      
      p_fix(1:ncre-1)  =  p_min_fix*ten**(w*dble(all_edges(1:ncre-1)-1))
      p_fix(0)    = zero     ! initial array of p used to identify fixed edges
      p_fix(ncre) = zero     ! extreme two edges are never fixed edges

      crel%p        = p_fix       ! actual array of p including free edges
      crel%p(0)     = p_lo
      crel%p(ncre)  = p_up
      
! Sorting bin edges - arbitrary chosen p_lo and p_up may need to be sorted to appear in growing order
      do k = ncre, 1, -1
         do i = 0, k-1
            if (crel%p(i)>crel%p(i+1)) then 
               c = crel%p(i)
               crel%p(i) = crel%p(i+1)
               crel%p(i+1) = c
            endif
         enddo
      enddo
      
      crel%f = f_init * (crel%p/p_min_fix)**(-q_init)

      crel%i_lo = 0
      crel%i_up = ncre
      
!       dt = 0.0d0
!       u_d = 0.0d0
      
!      call timestep(dt) <- will be called in the driver instead

       call cresp_update_bin_index(dt, p_lo, p_up, p_lo_next, p_up_next)

      
      
!      print*, 'In init_state'
!      print *, num_active_edges, size(active_edges), lbound(active_edges), ubound(active_edges)
!      print *, num_active_bins,  size(active_bins), lbound(active_bins), ubound(active_bins)
!      print *, 'p      =', active_edges(1:num_active_edges-1), & 
!                                                           crel%p(active_edges(1:num_active_edges-1))
!      print *, 'p      =', active_edges(2:num_active_edges),   &
!                                                           crel%p(active_edges(2:num_active_edges))
!      print *, 'init f =', active_edges(1:num_active_edges-1), & 
!                                                           crel%f(active_edges(1:num_active_edges-1))
!      print *, 'init q =', active_bins(1:num_active_bins),     & 
!                                                           crel%q(active_bins(1:num_active_bins))
      
      
       e = fq_to_e(crel%p(0:ncre-1), crel%p(1:ncre), crel%f(0:ncre-1), crel%q(1:ncre), active_bins)
!       
       n = fq_to_n(crel%p(0:ncre-1), crel%p(1:ncre), crel%f(0:ncre-1), crel%q(1:ncre), active_bins)


!       print *, 'init n =', n
!       print *, 'init e =', e 
! 
!        n_tot0 = sum(n)
!        e_tot0 = sum(e)
! 
!        print *, 'n_tot0 =', n_tot0
!        print *, 'e_tot0 =', e_tot0

   call cresp_deallocate_active_arrays
   
!   stop
   
    end subroutine cresp_init_state

   subroutine cosmicray_index(flind)

      use constants,    only: I_ONE
      use fluidtypes,   only: var_numbers

      implicit none

      type(var_numbers), intent(inout) :: flind
      integer(kind=4) :: icr

      flind%crn%beg    = flind%all + I_ONE
      flind%crs%beg    = flind%crn%beg

      
      flind%crn%all  = ncrn
      flind%cre%all  = 2*ncre+2 !!!
      flind%crs%all  = ncrs ! \crs deprecated, but we most likely need flind%crs%all to allocate part of cg%u responsible for all cosmic rays, so this remains unchanged
      

      do icr = 1, ncrn
         iarr_crn(icr)      = flind%all + icr
         iarr_crs(icr)      = flind%all + icr
      enddo
      flind%all = flind%all + flind%crn%all

      do icr = 1, (2*ncre+2)
         iarr_cre(icr)        = flind%all + icr
         iarr_crs(ncrn + icr) = flind%all + icr
      enddo
      flind%all = flind%all + flind%cre%all

      flind%crn%end = flind%crn%beg + flind%crn%all - I_ONE
      flind%cre%beg = flind%crn%end + I_ONE
      flind%cre%end = flind%all
      flind%crs%end = flind%cre%end
      if (flind%crn%all  /= 0) flind%components = flind%components + I_ONE
      flind%crn%pos = flind%components
      if (flind%cre%all  /= 0) flind%components = flind%components + I_ONE
      flind%cre%pos = flind%components

   end subroutine cosmicray_index
   

   subroutine cleanup_cosmicrays

      use diagnostics, only: my_deallocate

      implicit none
 
      call my_deallocate(iarr_crn)
      call my_deallocate(iarr_cre)
      call my_deallocate(iarr_crs)
      call my_deallocate(gamma_crs)
!       call my_deallocate(K_crs_paral)
!       call my_deallocate(K_crs_perp)

   end subroutine cleanup_cosmicrays

end module initcosmicrays