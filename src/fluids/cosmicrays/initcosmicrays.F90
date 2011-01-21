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
!! \brief (MH) Initialization of Cosmic Ray component
!!
!!
!!
!! In this module following namelist of parameters is specified:
!! \copydetails initcosmicrays::init_cosmicrays
!<
module initcosmicrays
! pulled by COSM_RAYS
   implicit none

   public ! QA_WARN no secrets are kept here

   integer, parameter                  :: ncr_max = 9  !< maximum number of CR nuclear and electron components

   ! namelist parameters
   integer                             :: ncrn         !< number of CR nuclear  components
   integer                             :: ncre         !< number of CR electron components
   integer                             :: ncrs         !< number of all CR components
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

   ! public component data
   integer, allocatable, dimension(:)  :: iarr_crn     !< array of indexes pointing to all CR nuclear components
   integer, allocatable, dimension(:)  :: iarr_cre     !< array of indexes pointing to all CR electron components
   integer, allocatable, dimension(:)  :: iarr_crs     !< array of indexes pointing to all CR components

   real,    allocatable, dimension(:)  :: gamma_crs    !< array containing adiabatic indexes of all CR components
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
   !! <tr><td>K_crn_paral </td><td>0     </td><td>real array</td><td>\copydoc initcosmicrays::K_crn_paral</td></tr>
   !! <tr><td>K_crn_perp  </td><td>0     </td><td>real array</td><td>\copydoc initcosmicrays::K_crn_perp </td></tr>
   !! <tr><td>K_cre_paral </td><td>0     </td><td>real array</td><td>\copydoc initcosmicrays::K_cre_paral</td></tr>
   !! <tr><td>K_cre_perp  </td><td>0     </td><td>real array</td><td>\copydoc initcosmicrays::K_cre_perp </td></tr>
   !! </table>
   !! \n \n
   !<
   subroutine init_cosmicrays

      use diagnostics,     only: ma1d, my_allocate
      use dataio_pub,      only: par_file, ierrh, namelist_errh, compare_namelist, cmdl_nml   ! QA_WARN required for diff_nml
      use dataio_pub,      only: die, warn
      use mpisetup,        only: master, slave, ibuff, rbuff, lbuff, comm, ierr, buffer_dim
      use mpi,             only: MPI_DOUBLE_PRECISION, MPI_INTEGER, MPI_LOGICAL

      implicit none

      integer                          :: nn
      integer                          :: ne

      namelist /COSMIC_RAYS/ cfl_cr, smallecr, cr_active, cr_eff, use_split, &
           &                 ncrn, gamma_crn, K_crn_paral, K_crn_perp, &
           &                 ncre, gamma_cre, K_cre_paral, K_cre_perp

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

      if (master) then
         diff_nml(COSMIC_RAYS) ! Do not use one-line if here!
      endif
#ifndef MULTIGRID
      if (.not. use_split) call warn("[initcosmicrays:init_cosmicrays] No multigrid solver compiled in: use_split reset to .true.")
      use_split  = .true.
#endif /* !MULTIGRID */
      rbuff(:)   = HUGE(1.)                         ! mark unused entries to allow automatic determination of nn

      if (master) then

         ibuff(1)   = ncrn
         ibuff(2)   = ncre

         rbuff(1)   = cfl_cr
         rbuff(2)   = smallecr
         rbuff(3)   = cr_active
         rbuff(4)   = cr_eff

         lbuff(1)   = use_split

         nn         = count(rbuff(:) < HUGE(1.))    ! this must match the last rbuff() index above
         ibuff(ubound(ibuff, 1)) = nn
         ne         = nn + 3 * ncrn
         if (ne + 3 * ncre > ubound(rbuff, 1)) call die("[initcosmicrays:init_cosmicrays] rbuff size exceeded.")

         if (ncrn > 0) then
            rbuff(nn+1       :nn+  ncrn) = gamma_crn  (1:ncrn)
            rbuff(nn+1+  ncrn:nn+2*ncrn) = K_crn_paral(1:ncrn)
            rbuff(nn+1+2*ncrn:nn+3*ncrn) = K_crn_perp (1:ncrn)
         endif

         if (ncre > 0) then
            rbuff(ne+1       :ne+  ncre) = gamma_cre  (1:ncre)
            rbuff(ne+1+  ncre:ne+2*ncre) = K_cre_paral(1:ncre)
            rbuff(ne+1+2*ncre:ne+3*ncre) = K_cre_perp (1:ncre)
         endif

      endif

      call MPI_Bcast(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_Bcast(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)
      call MPI_BCAST(lbuff,    buffer_dim, MPI_LOGICAL,          0, comm, ierr)

      if (slave) then

         ncrn       = ibuff(1)
         ncre       = ibuff(2)

         cfl_cr     = rbuff(1)
         smallecr   = rbuff(2)
         cr_active  = rbuff(3)
         cr_eff     = rbuff(4)

         use_split  = lbuff(1)

         nn         = ibuff(ubound(ibuff, 1))    ! this must match the last rbuff() index above
         ne         = nn + 3 * ncrn

         if (ncrn > 0) then
            gamma_crn  (1:ncrn) = rbuff(nn+1       :nn+  ncrn)
            K_crn_paral(1:ncrn) = rbuff(nn+1+  ncrn:nn+2*ncrn)
            K_crn_perp (1:ncrn) = rbuff(nn+1+2*ncrn:nn+3*ncrn)
         endif

         if (ncre > 0) then
            gamma_cre  (1:ncre) = rbuff(ne+1       :ne+  ncre)
            K_cre_paral(1:ncre) = rbuff(ne+1+  ncre:ne+2*ncre)
            K_cre_perp (1:ncre) = rbuff(ne+1+2*ncre:ne+3*ncre)
         endif

      endif

      ncrs = ncre + ncrn

      if (ncrs > ncr_max) call die("[initcosmicrays:init_cosmicrays] ncrs > ncr_max") !> \warning higher ncr_max limit would require changes in names of components in dataio_hdf5

      ma1d = [ncrs]
      call my_allocate(gamma_crs,   ma1d, "gamma_crs")
      call my_allocate(K_crs_paral, ma1d, "K_crs_paral")
      call my_allocate(K_crs_perp,  ma1d, "K_crs_perp")

      if (ncrn > 0) then
         gamma_crs  (1:ncrn) = gamma_crn  (1:ncrn)
         K_crs_paral(1:ncrn) = K_crn_paral(1:ncrn)
         K_crs_perp (1:ncrn) = K_crn_perp (1:ncrn)
      endif

      if (ncre > 0) then
         gamma_crs  (ncrn+1:ncrs) = gamma_cre  (1:ncre)
         K_crs_paral(ncrn+1:ncrs) = K_cre_paral(1:ncre)
         K_crs_perp (ncrn+1:ncrs) = K_cre_perp (1:ncre)
      endif

      ma1d = [ncrn]
      call my_allocate(iarr_crn, ma1d, "iarr_crn")
      ma1d = [ncre]
      call my_allocate(iarr_cre, ma1d, "iarr_cre")
      ma1d = [ncrs]
      call my_allocate(iarr_crs, ma1d, "iarr_crs")

   end subroutine init_cosmicrays

!

   subroutine cosmicray_index(flind)

      use grid,            only: cg
      use types,           only: var_numbers

      implicit none

      type(var_numbers), intent(inout) :: flind
      integer :: icr

      flind%crn%beg    = flind%all + 1
      flind%crs%beg    = flind%crn%beg

      flind%crn%all  = ncrn
      flind%cre%all  = ncre
      flind%crs%all  = ncrs

      do icr = 1, ncrn
         iarr_crn(icr)      =flind%all+icr
         iarr_crs(icr)      =flind%all+icr
      enddo
      flind%all = flind%all + flind%crn%all

      do icr = 1, ncre
         iarr_cre(icr)      =flind%all+icr
         iarr_crs(ncrn+icr) =flind%all+icr
      enddo
      flind%all = flind%all + flind%cre%all

      flind%crn%end    = flind%crn%beg + flind%crn%all - 1
      flind%cre%beg    = flind%crn%end + 1
      flind%cre%end    = flind%all
      flind%crs%end    = flind%cre%end
      if (flind%crn%all  /= 0) flind%components = flind%components + 1
      flind%crn%pos = flind%components
      if (flind%cre%all  /= 0) flind%components = flind%components + 1
      flind%cre%pos = flind%components

#ifdef NEW_HDF5
      call cr_add_hdf5(ncrs)
#else /* !NEW_HDF5 */
      if (.false.) icr = 0 * cg%is !suppress compiler warnings on unused arguments
#endif /* !NEW_HDF5 */

   end subroutine cosmicray_index

!

   subroutine cleanup_cosmicrays

      use diagnostics, only: my_deallocate

      implicit none

      call my_deallocate(iarr_crn)
      call my_deallocate(iarr_cre)
      call my_deallocate(iarr_crs)
      call my_deallocate(gamma_crs)
      call my_deallocate(K_crs_paral)
      call my_deallocate(K_crs_perp)

   end subroutine cleanup_cosmicrays

#ifdef NEW_HDF5
   subroutine cr_add_hdf5(flind_crs)

      use grid,      only: cg
      use arrays,    only: u
      use list_hdf5, only: add_lhdf5, lhdf5_info

      implicit none

      integer, intent(in)              :: flind_crs
      type(lhdf5_info) :: item
      integer          :: i

      item%p    => get_cr
      if (.not.allocated(item%ivec)) allocate(item%ivec(10))
      if (.not.allocated(item%rvec)) allocate(item%rvec(0))
      item%ivec  = [cg%nxb, cg%nyb, cg%nzb, cg%is, cg%ie, cg%js, cg%je, cg%ks, cg%ke]

      do i = 1, flind_crs

         write(item%key,'(A,I1)')  "ecr",i
         item%ivec(10) = iarr_crs(i)
         call add_lhdf5(item)

      enddo
   end subroutine cr_add_hdf5

   subroutine get_cr(ivec,rvec,outtab)
      use arrays,       only: u
      use dataio_pub,   only: die
      implicit none

      integer, dimension(:), intent(in)  :: ivec
      real,    dimension(:), intent(in)  :: rvec
      real, dimension(:,:,:), allocatable, intent(out) :: outtab

      if (allocated(outtab)) call die("[initcosmicrays:get_cr]: outtab already allocated")
      allocate(outtab(ivec(1),ivec(2),ivec(3)))
      outtab(:,:,:) = u(ivec(10),ivec(4):ivec(5),ivec(6):ivec(7),ivec(8):ivec(9))
      return
      if (.false.) write(0,*) rvec
   end subroutine get_cr

#endif /* NEW_HDF5 */

end module initcosmicrays
