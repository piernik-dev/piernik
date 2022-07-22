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
!! \brief Module containing subroutines and functions that govern supernovae insert
!<
module snsources
! pulled by SN_SRC

   use constants, only: ndims, LO, HI

   implicit none

   private
   public :: random_sn, init_snsources, r_sn, nsn
#ifdef COSM_RAYS
   public :: cr_sn, amp_ecr_sn
#endif /* COSM_RAYS */
#ifdef HDF5
   public :: read_snsources_from_restart, write_snsources_to_restart
#endif /* HDF5 */
#ifdef SHEAR
   public :: sn_shear
#endif /* SHEAR */

   integer(kind=4)                         :: nsn, nsn_last
   real                                    :: e_sn                !< energy per supernova [erg]
   real                                    :: f_sn                !< frequency of SN
   real                                    :: f_sn_kpc2           !< frequency of SN per kpc^2
   real                                    :: h_sn                !< galactic height in SN gaussian distribution ?
   real                                    :: r_sn                !< radius of SN
   integer(kind=4), dimension(ndims,LO:HI) :: auxper
   real                                    :: gnorm               !< gauss distribution normalization factor
#ifdef COSM_RAYS
   real                                    :: amp_ecr_sn          !< cosmic ray explosion amplitude
#endif /* COSM_RAYS */

   namelist /SN_SOURCES/ e_sn, h_sn, r_sn, f_sn_kpc2

contains
!>
!! \brief Routine to set parameter values from namelist SN_SOURCES
!!
!! \n \n
!! @b SN_SOURCES
!! \n \n
!! <table border="+1">
!! <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!! <tr><td>e_sn     </td><td>4.96e6</td><td>real value</td><td>\copydoc snsources::e_sn     </td></tr>
!! <tr><td>h_sn     </td><td>0.0   </td><td>real value</td><td>\copydoc snsources::h_sn     </td></tr>
!! <tr><td>r_sn     </td><td>0.0   </td><td>real value</td><td>\copydoc snsources::r_sn     </td></tr>
!! <tr><td>f_sn_kpc2</td><td>0.0   </td><td>real value</td><td>\copydoc snsources::f_sn_kpc2</td></tr>
!! </table>
!! The list is active while \b "SN_SRC" is defined.
!! \n \n
!<
   subroutine init_snsources

      use constants,      only: PIERNIK_INIT_GRID, xdim, ydim, zdim, pi, I_ONE
      use dataio_pub,     only: die, code_progress, nh
      use domain,         only: dom
      use mpisetup,       only: rbuff, master, slave, piernik_MPI_Bcast
      use units,          only: erg, kpc, Myr
#ifdef COSM_RAYS
      use initcosmicrays, only: cr_eff
#endif /* COSM_RAYS */

      implicit none

      integer(kind=4) :: id

      if (code_progress < PIERNIK_INIT_GRID) call die("[snsources:init_snsources] grid or fluids/cosmicrays not initialized.")

      e_sn      = 1.e51
      h_sn      = 0.0
      r_sn      = 0.0
      f_sn_kpc2 = 0.0

      if (master) then
         if (.not.nh%initialized) call nh%init()
         open(newunit=nh%lun, file=nh%tmp1, status="unknown")
         write(nh%lun,nml=SN_SOURCES)
         close(nh%lun)
         open(newunit=nh%lun, file=nh%par_file)
         nh%errstr=""
         read(unit=nh%lun, nml=SN_SOURCES, iostat=nh%ierrh, iomsg=nh%errstr)
         close(nh%lun)
         call nh%namelist_errh(nh%ierrh, "SN_SOURCES")
         read(nh%cmdl_nml,nml=SN_SOURCES, iostat=nh%ierrh)
         call nh%namelist_errh(nh%ierrh, "SN_SOURCES", .true.)
         open(newunit=nh%lun, file=nh%tmp2, status="unknown")
         write(nh%lun,nml=SN_SOURCES)
         close(nh%lun)
         call nh%compare_namelist()

         rbuff(1) = e_sn
         rbuff(2) = h_sn
         rbuff(3) = r_sn
         rbuff(4) = f_sn_kpc2
      endif

      call piernik_MPI_Bcast(rbuff)

      if (slave) then
         e_sn      = rbuff(1)
         h_sn      = rbuff(2)
         r_sn      = rbuff(3)
         f_sn_kpc2 = rbuff(4)
      endif

      if (all(dom%has_dir)) then
         gnorm = 1./(16./9.*pi*r_sn**3)
      else
         gnorm = 1./(pi*r_sn**2)
      endif

#ifdef COSM_RAYS
      amp_ecr_sn = cr_eff * e_sn * erg * gnorm
#endif /* COSM_RAYS */

      f_sn = f_sn_kpc2 / Myr * dom%L_(xdim) * dom%L_(ydim) / kpc**2

      nsn_last = 0
      nsn      = 0

      auxper = 0
      do id = xdim, zdim
         if (dom%periodic(id) .and. dom%has_dir(id)) auxper(id,:) = [-I_ONE, I_ONE]
      enddo

   end subroutine init_snsources

!>
!! \brief Main routine to insert one supernova event
!<
   subroutine random_sn

      use global,     only: t
      use repeatstep, only: repeat_step

      implicit none

      real, dimension(ndims) :: snpos
      integer                :: isn, nsn_per_timestep

      if (.not. repeat_step()) nsn_last = nsn

      nsn = int(t * f_sn, kind=4)
      nsn_per_timestep = nsn - nsn_last

      do isn = 1, nsn_per_timestep

         call rand_coords(snpos)

#ifdef COSM_RAYS
         call cr_sn(snpos, amp_ecr_sn)
#endif /* COSM_RAYS */

      enddo

   end subroutine random_sn

!--------------------------------------------------------------------------
#ifdef COSM_RAYS
!>
!! \brief Routine that inserts an amount of cosmic ray energy around the position of supernova
!! \param pos real, dimension(3), array of supernova position components
!<
   subroutine cr_sn(pos, ampl)

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: xdim, ydim, zdim
      use domain,           only: dom
      use grid_cont,        only: grid_container
      use cr_data,          only: cr_index, cr_table, cr_mass, cr_primary, eCRSP, icr_H1, icr_C12, icr_N14, icr_O16, rel_abound
      use initcosmicrays,   only: iarr_crn
#ifdef CRESP
      use cresp_crspectrum, only: cresp_get_scaled_init_spectrum
      !use initcrspectrum,   only: cresp, cre_eff, smallcree, use_cresp, nspc
      !use initcosmicrays,   only: iarr_crspc_n, iarr_crspc_e
      use initcosmicrays,   only: iarr_crspc2_e, iarr_crspc2_n, nspc, iarr_crspc
      use initcrspectrum,   only: expan_order, smallcree, cresp, cre_eff, use_cresp
#endif /* CRESP */

      implicit none

      real, dimension(ndims), intent(in) :: pos
      real,                   intent(in) :: ampl
      integer                            :: i, j, k, ipm, jpm, icr
      real                               :: decr, ysna
      real, dimension(ndims)             :: posr
      type(cg_list_element), pointer     :: cgl
      type(grid_container),  pointer     :: cg
#ifdef SHEAR
      real, dimension(3)                 :: ysnoi
#endif /* SHEAR */
#ifdef CRESP
      real                               :: e_tot_sn
#endif /* CRESP */

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

#ifdef SHEAR
         ysnoi(2) = pos(ydim)
         call sn_shear(cg, ysnoi)
#else /* !SHEAR */
         ysna = pos(ydim)
#endif /* !SHEAR */

         do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
            posr(zdim) = ((cg%z(k)-pos(zdim))/r_sn)**2
            do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
               do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)

                  decr = 0.0
                  do ipm = auxper(xdim,LO), auxper(xdim,HI)
                     posr(xdim) = ((cg%x(i)-pos(xdim) + real(ipm)*dom%L_(xdim))/r_sn)**2
#ifdef SHEAR
                     ysna = ysnoi(ipm+2)
#endif /* SHEAR */
                     do jpm = auxper(ydim,LO), auxper(ydim,HI)
                        posr(ydim) = ((cg%y(j)-ysna + real(jpm)*dom%L_(ydim))/r_sn)**2
                        ! BEWARE:  for num < -744.6 the exp(num) is the underflow
                        decr = decr + exp(-sum(posr, mask=dom%has_dir))
                     enddo
                  enddo
                  decr = decr * ampl


#ifdef CRESP
                  if (use_cresp) then

!                     e_tot_sn = decr * cre_eff
!                     if (e_tot_sn > smallcree) then!
!                        cresp%n =  0.0;  cresp%e = 0.0
!                        call cresp_get_scaled_init_spectrum(cresp%n, cresp%e, e_tot_sn) !< injecting source spectrum scaled with e_tot_sn
!                        cg%u(iarr_crspc_n,i,j,k) = cg%u(iarr_crspc_n,i,j,k) + cresp%n   !< update, TODO need to talk to the team if this should be inside if-clause
!                        cg%u(iarr_crspc_e,i,j,k) = cg%u(iarr_crspc_e,i,j,k) + cresp%e   !< if outside, cresp%n and cresp%e needs to be zeroed
!                     endif
!                  endif

                    do icr = 1, nspc
                        if (icr==icr_H1) e_tot_sn = cre_eff(icr) * decr
                        if (e_tot_sn > smallcree .and. use_cresp) then
!                         cresp%n = 1.e-4 ;  cresp%e = 1.e-2
                            call cresp_get_scaled_init_spectrum(cresp%n, cresp%e, e_tot_sn, icr)
                        !stop

                            cg%u(iarr_crspc2_n(icr,:),i,j,k) = cg%u(iarr_crspc2_n(icr,:),i,j,k) + rel_abound(icr)*cresp%n
                            cg%u(iarr_crspc2_e(icr,:),i,j,k) = cg%u(iarr_crspc2_e(icr,:),i,j,k) + rel_abound(icr)*cresp%e
                        !print *, 'cresp n ' , cresp%n,  'cresp e ' , cresp%e
                        !print *, 'i : ', i, ' j :', j
                        !if (i == 29 .and. j == 29) then
                           !print *, icr
                           !print *, ' cresp%n ', rel_abound(icr)*cresp%n
                           !print *, 'cresp%e ', rel_abound(icr)*cresp%e
                           !endif
                         endif
                    enddo
                  endif
#else /* CRESP */
               if (eCRSP(icr_H1 )) cg%u(iarr_crn(cr_index(icr_H1 )),i,j,k) = cg%u(iarr_crn(cr_index(icr_H1 )),i,j,k) + decr
               if (eCRSP(icr_C12)) cg%u(iarr_crn(cr_index(icr_C12)),i,j,k) = cg%u(iarr_crn(cr_index(icr_C12)),i,j,k) + cr_primary(cr_table(icr_C12)) * cr_mass(cr_table(icr_C12)) * decr
               if (eCRSP(icr_N14)) cg%u(iarr_crn(cr_index(icr_N14)),i,j,k) = cg%u(iarr_crn(cr_index(icr_N14)),i,j,k) + cr_primary(cr_table(icr_N14)) * cr_mass(cr_table(icr_N14)) * decr
               if (eCRSP(icr_O16)) cg%u(iarr_crn(cr_index(icr_O16)),i,j,k) = cg%u(iarr_crn(cr_index(icr_O16)),i,j,k) + cr_primary(cr_table(icr_O16)) * cr_mass(cr_table(icr_O16)) * decr
#endif /* CRESP */
               enddo
            enddo
         enddo

         cgl => cgl%nxt
      enddo

   end subroutine cr_sn
#endif /* COSM_RAYS */
!--------------------------------------------------------------------------
!>
!! \brief Routine that determines the position of next supernova
!! \return pos @e real,  @e dimension(3), array of supernova position components
!<
   subroutine rand_coords(pos)

      use constants,   only: xdim, ydim, zdim
      use domain,      only: dom

      implicit none

      real, dimension(ndims), intent(out) :: pos
      real, dimension(4)                  :: rand

      call random_number(rand)
      pos(xdim:ydim) = dom%edge(xdim:ydim, LO)+ dom%L_(xdim:ydim)*rand(xdim:ydim)

      if (dom%has_dir(zdim)) then
         pos(zdim) = h_sn * gasdev(rand(3),rand(4))
      else
         pos(zdim) = 0.0
      endif

   end subroutine rand_coords

!-----------------------------------------------------------------------
#ifdef SHEAR
   subroutine sn_shear(cg, ysnoi)

      use constants,   only: ydim
      use dataio_pub,  only: die
      use domain,      only: is_multicg, dom
      use grid_cont,   only: grid_container
      use shear,       only: delj, eps

      implicit none

      type(grid_container), pointer     :: cg
      real, dimension(3), intent(inout) :: ysnoi
      integer                           :: jsn, jremap
      real                              :: dysn, epsi, epso

      if (is_multicg) call die("[snsources:sn_shear] multiple grid pieces per processor not implemented yet") !nontrivial SHEAR

      jsn  = cg%js+int((ysnoi(2)-dom%edge(ydim, LO))*cg%idy)
      dysn  = mod(ysnoi(2), cg%dy)

      epsi   = eps*cg%dy
      epso   = -epsi

!  outer boundary
      jremap = jsn - delj
      jremap = mod(mod(jremap, cg%nyb)+cg%nyb, cg%nyb)
      if (jremap <= (cg%lh1(ydim,LO))) jremap = jremap + cg%nyb

      ysnoi(1) = cg%y(jremap) + epso + dysn

!  inner boundary
      jremap = jsn + delj
      jremap = mod(jremap, cg%nyb)+cg%nyb
      if (jremap >= (cg%lh1(ydim,HI))) jremap = jremap - cg%nyb

      ysnoi(3) = cg%y(jremap) + epsi + dysn

   end subroutine sn_shear
#endif /* SHEAR */
!-----------------------------------------------------------------------

!>
!! \brief Function that generates values of normal distribution (from Numerical Recipes)
!! \return x random value of uniform distribution
!! \return y random value of uniform distribution
!! \return @e real, random value of normal distribution
!! \todo Change function @c gasdev into another one
!<
      function gasdev(x,y)

         implicit none

         real               :: x, y, x1, y1, r, fac, gasdev
         real,    save      :: gset
         integer, save      :: iset = 0
         real, dimension(2) :: rand

         r = 2.0
         if (iset == 0) then
            do
               x1 = 2.0*x - 1.0
               y1 = 2.0*y - 1.0
               r  = x1**2 + y1**2
               if (r >= 1.0) then
                  call random_number(rand)
                  x = rand(1)
                  y = rand(2)
               else
                  exit
               endif
            enddo
            fac=sqrt(-2.*log(r)/r)
            gset=x1*fac
            gasdev=y1*fac
            iset=1
         else
            gasdev=gset
            iset=0
         endif

      end function gasdev

#ifdef HDF5
   subroutine write_snsources_to_restart(file_id)

      use hdf5, only: HID_T, SIZE_T
      use h5lt, only: h5ltset_attribute_int_f

      implicit none

      integer(HID_T), intent(in) :: file_id
      integer(SIZE_T)            :: bufsize
      integer(kind=4)            :: error
      integer(kind=4), dimension(1) :: lnsnbuf

      bufsize = 1
      lnsnbuf(bufsize) = nsn
      call h5ltset_attribute_int_f(file_id, "/", "nsn", lnsnbuf, bufsize, error)

   end subroutine write_snsources_to_restart

!-----------------------------------------------------------------------------

   subroutine read_snsources_from_restart(file_id)

      use cmp_1D_mpi, only: compare_array1D
      use hdf5,       only: HID_T
      use h5lt,       only: h5ltget_attribute_int_f

      implicit none

      integer(HID_T), intent(in)         :: file_id
      integer(kind=4)                    :: error
      integer(kind=4), dimension(:), allocatable :: lnsnbuf

      if (.not.allocated(lnsnbuf)) allocate(lnsnbuf(1))
      call h5ltget_attribute_int_f(file_id, "/", "nsn", lnsnbuf, error)
      call compare_array1D(lnsnbuf)
      nsn = lnsnbuf(1)
      deallocate(lnsnbuf)

   end subroutine read_snsources_from_restart
#endif /* HDF5 */

end module snsources
