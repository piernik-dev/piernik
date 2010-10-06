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
!#include "piernik.def"
#include "defines.c"

!>
!! \brief (MH) Initialization of Cosmic Ray component
!!
!!
!!
!! In this module following namelist of parameters is specified:
!! \copydetails initcosmicrays::init_cosmicrays
!<
module initcosmicrays

   implicit none

   integer, parameter                  :: ncr_max = 9  !< maximum number of CR nuclear and electron components
   integer                             :: ncrn         !< number of CR nuclear  components
   integer                             :: ncre         !< number of CR electron components

   real                                :: cfl_cr       !< CFL number for diffusive CR transport
   real                                :: smallecr     !< floor value for CR energy density
   real                                :: cr_active    !< parameter specifying whether CR pressure gradient is (when =1.) or isn't (when =0.) included in the gas equation of motion
   real                                :: cr_eff       !< conversion rate of SN explosion energy to CR energy (default = 0.1)

   real, dimension(ncr_max)            :: gamma_crn    !< array containing adiabatic indexes of all CR nuclear components
   real, dimension(ncr_max)            :: K_crn_paral  !< array containing parallel diffusion coefficients of all CR nuclear components
   real, dimension(ncr_max)            :: K_crn_perp   !< array containing perpendicular diffusion coefficients of all CR nuclear components
   real, dimension(ncr_max)            :: gamma_cre    !< array containing adiabatic indexes of all CR nuclear components
   real, dimension(ncr_max)            :: K_cre_paral  !< array containing parallel diffusion coefficients of all CR electron components
   real, dimension(ncr_max)            :: K_cre_perp   !< array containing perpendicular diffusion coefficients of all CR electron components

   integer, allocatable, dimension(:)  :: iarr_crn     !< array of indexes pointing to all CR nuclear components
   integer, allocatable, dimension(:)  :: iarr_cre     !< array of indexes pointing to all CR electron components
   integer, allocatable, dimension(:)  :: iarr_crs     !< array of indexes pointing to all CR components

   real,    allocatable, dimension(:)  :: gamma_crs    !< array containing adiabatic indexes of all CR components
   real,    allocatable, dimension(:)  :: K_crs_paral  !< array containing parallel diffusion coefficients of all CR components
   real,    allocatable, dimension(:)  :: K_crs_perp   !< array containing perpendicular diffusion coefficients of all CR components

   !BEWARE Possible confusion: *_perp coefficients are not "perpendicular" but rather isotropic

contains

   !>
   !! \brief Routine to set parameters values from namelist COSMIC_RAYS
   !!
   !! \n \n
   !! @b COSMIC_RAYS
   !! \n \n
   !! <table border="+1">
   !! <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
   !! <tr><td>cfl_cr    </td><td>0.9  </td><td>real value</td><td>\copydoc initcosmicrays::cfl_cr</td></tr>
   !! <tr><td>smallecr  </td><td>0.0  </td><td>real value</td><td>\copydoc initcosmicrays::smallecr</td></tr>
   !! <tr><td>cr_active </td><td>1.0  </td><td>real value</td><td>\copydoc initcosmicrays::cr_active</td></tr>
   !! <tr><td>gamma_cr  </td><td>4/3  </td><td>real value</td><td>\copydoc initcosmicrays::gamma_cr</td></tr>
   !! <tr><td>cr_eff    </td><td>0.1  </td><td>real value</td><td>\copydoc initcosmicrays::cr_eff</td></tr>
   !! <tr><td>K_cr_paral</td><td>0.0  </td><td>real value</td><td>\copydoc initcosmicrays::K_cr_paral</td></tr>
   !! <tr><td>K_cr_perp </td><td>0.0  </td><td>real value</td><td>\copydoc initcosmicrays::K_cr_perp</td></tr>
   !! </table>
   !! \n \n
   !<
   subroutine init_cosmicrays

      use errh,     only: namelist_errh, die
      use mpisetup, only: proc, ibuff, rbuff, comm, ierr, MPI_DOUBLE_PRECISION, MPI_INTEGER, buffer_dim
      use dataio_public, only: par_file, cwd
      use func,        only : compare_namelist

      implicit none

      integer            :: ierrh
      integer            :: nn
      integer            :: ne

      namelist /COSMIC_RAYS/ cfl_cr, smallecr, cr_active, cr_eff, &
           &                 ncrn, gamma_crn, K_crn_paral, K_crn_perp, &
           &                 ncre, gamma_cre, K_cre_paral, K_cre_perp

      cfl_cr     = 0.9
      smallecr   = 0.0
      cr_active  = 1.0
      cr_eff     = 0.1       !  canonical conversion rate of SN en.-> CR
      !  we fix E_SN=10**51 erg
      ncrn       = 0
      ncre       = 0

      gamma_crn(:)   = 4./3.
      K_crn_paral(:) = 0.0
      K_crn_perp(:)  = 0.0
      gamma_cre(:)   = 4./3.
      K_cre_paral(:) = 0.0
      K_cre_perp(:)  = 0.0

      if (proc == 0) then

         diff_nml(COSMIC_RAYS)

         ibuff(1)   = ncrn
         ibuff(2)   = ncre

         rbuff(1)   = cfl_cr
         rbuff(2)   = smallecr
         rbuff(3)   = cr_active
         rbuff(4)   = cr_eff

         nn         = 4         ! WARNING: this must match the last rbuff() index above
         ne         = nn + 3 * ncrn
         if (ne + 3 * ncre > ubound(rbuff, 1)) call die("[initcosmicrays:init_cosmicrays] rbuff size exceeded.")

         if(ncrn > 0) then
            rbuff(nn+1       :nn+  ncrn) = gamma_crn  (1:ncrn)
            rbuff(nn+1+  ncrn:nn+2*ncrn) = K_crn_paral(1:ncrn)
            rbuff(nn+1+2*ncrn:nn+3*ncrn) = K_crn_perp (1:ncrn)
         endif

         if(ncre > 0) then
            rbuff(ne+1       :ne+  ncre) = gamma_cre  (1:ncre)
            rbuff(ne+1+  ncre:ne+2*ncre) = K_cre_paral(1:ncre)
            rbuff(ne+1+2*ncre:ne+3*ncre) = K_cre_perp (1:ncre)
         endif

      end if

      call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

      if (proc /= 0) then

         ncrn       = ibuff(1)
         ncre       = ibuff(2)

         cfl_cr     = rbuff(1)
         smallecr   = rbuff(2)
         cr_active  = rbuff(3)
         cr_eff     = rbuff(4)

         nn         = 4         ! WARNING: this must match the last rbuff() index above
         ne         = nn + 3 * ncrn

         if(ncrn > 0) then
            gamma_crn  (1:ncrn) = rbuff(nn+1       :nn+  ncrn)
            K_crn_paral(1:ncrn) = rbuff(nn+1+  ncrn:nn+2*ncrn)
            K_crn_perp (1:ncrn) = rbuff(nn+1+2*ncrn:nn+3*ncrn)
         endif

         if(ncre > 0) then
            gamma_cre  (1:ncre) = rbuff(ne+1       :ne+  ncre)
            K_cre_paral(1:ncre) = rbuff(ne+1+  ncre:ne+2*ncre)
            K_cre_perp (1:ncre) = rbuff(ne+1+2*ncre:ne+3*ncre)
         endif

      endif

      if (ncrn > ncr_max) call die("[initcosmicrays:init_cosmicrays] ncrn > ncr_max")
      if (ncre > ncr_max) call die("[initcosmicrays:init_cosmicrays] ncre > ncr_max")

      allocate(gamma_crs(ncrn+ncre),K_crs_paral(ncrn+ncre),K_crs_perp(ncrn+ncre))

      if(ncrn > 0) then
         gamma_crs  (1:ncrn) = gamma_crn  (1:ncrn)
         K_crs_paral(1:ncrn) = K_crn_paral(1:ncrn)
         K_crs_perp (1:ncrn) = K_crn_perp (1:ncrn)
      endif

      if(ncre > 0) then
         gamma_crs  (ncrn+1:ncrn+ncre) = gamma_cre  (1:ncre)
         K_crs_paral(ncrn+1:ncrn+ncre) = K_cre_paral(1:ncre)
         K_crs_perp (ncrn+1:ncrn+ncre) = K_cre_perp (1:ncre)
      endif

   end subroutine init_cosmicrays

   subroutine cosmicray_species

   end subroutine cosmicray_species

   subroutine cosmicray_index(nvar, nvar_crn, nvar_cre, nvar_crs)

      implicit none

      integer :: nvar
      integer :: nvar_crn
      integer :: nvar_cre
      integer :: nvar_crs
      integer :: icr

      nvar_crn      = ncrn
      nvar_cre      = ncre
      nvar_crs      = ncrn + ncre

      allocate(iarr_crn(ncrn))
      allocate(iarr_cre(ncre))
      allocate(iarr_crs(nvar_crs))

      do icr = 1, ncrn
         iarr_crn(icr)      =nvar+icr
         iarr_crs(icr)      =nvar+icr
      enddo
      nvar = nvar + nvar_crn

      do icr = 1, ncre
         iarr_cre(icr)      =nvar+icr
         iarr_crs(ncrn+icr) =nvar+icr
      enddo
      nvar = nvar + nvar_cre

#ifdef NEW_HDF5
      call cr_add_hdf5(ncrn+ncre)
#endif /* NEW_HDF5 */

   end subroutine cosmicray_index

   subroutine cleanup_cosmicrays

      implicit none

      if (allocated(iarr_crn)) deallocate(iarr_crn)
      if (allocated(iarr_cre)) deallocate(iarr_cre)
      if (allocated(iarr_crs)) deallocate(iarr_crs)
      if (allocated(gamma_crs)) deallocate(gamma_crs, K_crs_paral, K_crs_perp) !BEWARE: simplified allocated() check

   end subroutine cleanup_cosmicrays

#ifdef NEW_HDF5
   subroutine cr_add_hdf5(nvar_crs)

      use arrays,    only : u
      use grid,      only : nxb, nyb, nzb, is, ie, js, je, ks, ke
      use list_hdf5, only : add_lhdf5, lhdf5_info

      implicit none

      integer, intent(in) :: nvar_crs

      type(lhdf5_info) :: item
      integer          :: i

      item%p    => get_cr
      item%sz   = (/nxb, nyb, nzb/)
      item%ind  = (/is, ie, js, je, ks, ke/)

      do i = 1, nvar_crs

         write(item%key,'(A,I1)')  "ecr",i
         item%opti = (/iarr_crs(i),0,0,0,0/)       !< optional integer passed to func
         call add_lhdf5(item)

      enddo

   contains

      function get_cr(i,sz,opti) result (outtab)

         implicit none

         integer, dimension(6), intent(in)    :: i
         integer, dimension(3), intent(in)    :: sz
         integer, dimension(5), intent(in)    :: opti
         real, dimension(sz(1),sz(2),sz(3))   :: outtab

         outtab(:,:,:) = u(opti(1),i(1):i(2),i(3):i(4),i(5):i(6))

      end function get_cr

   end subroutine cr_add_hdf5
#endif /* NEW_HDF5 */

end module initcosmicrays
