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
!! \brief Module containing subroutines and functions that govern supernovae insert
!<
module snsources
! pulled by SN_SRC
   implicit none
   private
   public ::  random_sn, init_snsources, r_sn
#ifdef SHEAR
   public :: sn_shear

   real, dimension(3) :: ysnoi
#endif /* SHEAR */
   integer, save      :: nsn, nsn_last
   real,    save      :: dt_sn_prev, ecr_supl, decr_supl
   real,    save      :: gset
   integer, save      :: irand, iset

   real, parameter    :: ethu = 7.0**2/(5.0/3.0-1.0) * 1.0    !< thermal energy unit=0.76eV/cm**3 for c_si= 7km/s, n=1/cm^3 gamma=5/3

   real               :: amp_ecr_sn          !< cosmic ray explosion amplitude in units: e_0 = 1/(5/3-1)*rho_0*c_s0**2  rho_0=1.67e-24g/cm**3, c_s0 = 7km/s
   real               :: f_sn                !< frequency of SN
   real               :: f_sn_kpc2           !< frequency of SN per kpc^2
   real               :: h_sn                !< galactic height in SN gaussian distribution ?
   real               :: r_sn                !< radius of SN

!   namelist /SN_SOURCES/ amp_ecr_sn, f_sn, h_sn, r_sn, f_sn_kpc2
   namelist /SN_SOURCES/ h_sn, r_sn, f_sn_kpc2

contains
!>
!! \brief Routine to set parameter values from namelist SN_SOURCES
!!
!! \n \n
!! @b SN_SOURCES
!! \n \n
!! <table border="+1">
!! <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!! <tr><td>h_sn     </td><td>0.0  </td><td>real value</td><td>\copydoc snsources::h_sn     </td></tr>
!! <tr><td>r_sn     </td><td>0.0  </td><td>real value</td><td>\copydoc snsources::r_sn     </td></tr>
!! <tr><td>f_sn_kpc2</td><td>0.0  </td><td>real value</td><td>\copydoc snsources::f_sn_kpc2</td></tr>
!! </table>
!! The list is active while \b "SN_SRC" is defined.
!! \n \n
!<
   subroutine init_snsources

      use constants,      only: PIERNIK_INIT_GRID, xdim, ydim
      use dataio_pub,     only: ierrh, par_file, namelist_errh, compare_namelist, cmdl_nml, lun                  ! QA_WARN required for diff_nml
      use dataio_pub,     only: die, code_progress
      use domain,         only: dom
      use initcosmicrays, only: cr_eff
      use mpisetup,       only: rbuff, master, slave, piernik_MPI_Bcast

      implicit none

      if (code_progress < PIERNIK_INIT_GRID) call die("[snsources:init_snsources] grid or fluids/cosmicrays not initialized.")

!      amp_ecr_sn = 0.0    !> \todo set sane default values
      f_sn       = 0.0    !
      h_sn       = 0.0
      r_sn       = 0.0

      if (master) then
         diff_nml(SN_SOURCES)
!         rbuff(1)   = amp_ecr_sn
!         rbuff(2)   = f_sn
         rbuff(3)   = h_sn
         rbuff(4)   = r_sn
         rbuff(5)   = f_sn_kpc2
      endif

      call piernik_MPI_Bcast(rbuff)

      if (slave) then
!        amp_ecr_sn  = rbuff(1)
!        f_sn        = rbuff(2)
         h_sn        = rbuff(3)
         r_sn        = rbuff(4)
         f_sn_kpc2   = rbuff(5)
      endif

      amp_ecr_sn = 4.96e6*cr_eff/r_sn**3

      if (dom%has_dir(xdim)) then
         f_sn = f_sn_kpc2 * dom%L_(xdim)/1000.0 !\deprecated magic numbers
      else
         f_sn = f_sn_kpc2 * 2.0*r_sn/1000.0
      endif

      if (dom%has_dir(ydim)) then
         f_sn = f_sn * dom%L_(ydim)/1000.0
      else
         f_sn = f_sn * 2.0*r_sn/1000.0
      endif

   end subroutine init_snsources
!>
!! \brief Main routine to insert one supernova event
!<
   subroutine random_sn

      use constants, only: small
      use global,    only: t

      implicit none
      real :: dt_sn
!      real, dimension(2) :: orient
      real, dimension(3) :: snpos
      integer :: isn, nsn_per_timestep

      dt_sn = 1./(f_sn+small)

      nsn = int(t/dt_sn)
      nsn_per_timestep = nsn - nsn_last
      nsn_last = nsn

      do isn = 1, nsn_per_timestep

         call rand_coords(snpos)

#ifdef COSM_RAYS
         call cr_sn(snpos)
#endif /* COSM_RAYS */

      enddo ! isn
      return
   end subroutine random_sn

!--------------------------------------------------------------------------
!>
!! \brief Routine that inserts an amount of cosmic ray energy around the position of supernova
!! \param pos real, dimension(3), array of supernova position components
!<
   subroutine cr_sn(pos)

      use cg_leaves,      only: leaves
      use cg_list,        only: cg_list_element
      use constants,      only: xdim, ydim, zdim
      use domain,         only: dom
      use grid_cont,      only: grid_container
#ifdef COSM_RAYS_SOURCES
      use cr_data,        only: icr_H1, icr_C12, icr_N14, icr_O16, cr_table, cr_primary, eCRSP
      use initcosmicrays, only: iarr_crn
#endif /* COSM_RAYS_SOURCES */

      implicit none

      real, dimension(3), intent(in) :: pos
      integer                        :: i, j, k, ipm, jpm
      real                           :: decr, xsn, ysn, zsn, ysna, zr
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg

      xsn = pos(1)
      ysn = pos(2)
      zsn = pos(3)

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         do k=1, cg%n_(zdim)
            zr = (cg%z(k)-zsn)**2
            do j=1, cg%n_(ydim)
               do i=1, cg%n_(xdim)

                  decr = 0.0
                  do ipm=-1,1
#ifdef SHEAR
                     ysna = ysnoi(ipm+2)
#else /* !SHEAR */
                     ysna = ysn
#endif /* !SHEAR */
                     do jpm=-1,1
                        decr = decr + exp(-((cg%x(i)-xsn +real(ipm)*dom%L_(xdim))**2  &
                             +              (cg%y(j)-ysna+real(jpm)*dom%L_(ydim))**2 + zr)/r_sn**2)
                     enddo
                  enddo
                  decr = decr * amp_ecr_sn *ethu

#ifdef COSM_RAYS_SOURCES
                  if (eCRSP(icr_H1 )) cg%u(iarr_crn(cr_table(icr_H1 )),i,j,k) = cg%u(iarr_crn(cr_table(icr_H1 )),i,j,k) + decr
                  if (eCRSP(icr_C12)) cg%u(iarr_crn(cr_table(icr_C12)),i,j,k) = cg%u(iarr_crn(cr_table(icr_C12)),i,j,k) + cr_primary(cr_table(icr_C12))*12*decr
                  if (eCRSP(icr_N14)) cg%u(iarr_crn(cr_table(icr_N14)),i,j,k) = cg%u(iarr_crn(cr_table(icr_N14)),i,j,k) + cr_primary(cr_table(icr_N14))*14*decr
                  if (eCRSP(icr_O16)) cg%u(iarr_crn(cr_table(icr_O16)),i,j,k) = cg%u(iarr_crn(cr_table(icr_O16)),i,j,k) + cr_primary(cr_table(icr_O16))*16*decr
#endif /* COSM_RAYS_SOURCES */

               enddo
            enddo
         enddo

         cgl => cgl%nxt
      enddo

   end subroutine cr_sn

!--------------------------------------------------------------------------
!>
!! \brief Routine that determines the position of next supernova
!! \return pos @e real,  @e dimension(3), array of supernova position components
!<
   subroutine rand_coords(pos)

      use constants,   only: xdim, ydim, zdim, LO
      use domain,      only: dom
#ifdef SHEAR
      use cg_leaves,   only: leaves
      use grid_cont,   only: grid_container
#endif /* SHEAR */

      implicit none

#ifdef SHEAR
      type(grid_container), pointer   :: cg
#endif /* SHEAR */
      real, dimension(3), intent(out) :: pos
      real, dimension(4)              :: rand
      real                            :: xsn, ysn, zsn, znorm

      call random_number(rand)
      xsn = dom%edge(xdim, LO)+ dom%L_(xdim)*rand(1)
      ysn = dom%edge(ydim, LO)+ dom%L_(ydim)*rand(2)

      if (dom%has_dir(zdim)) then
         irand = irand+4
         znorm = gasdev(rand(3),rand(4))
         zsn = h_sn*znorm
      else
         zsn = 0.0
      endif

#ifdef SHEAR
      cg => leaves%first%cg
      ysnoi = 0.0 ; ysnoi(2) = ysn
      call sn_shear(cg, ysnoi)
#endif /* SHEAR */

      pos(1) = xsn
      pos(2) = ysn
      pos(3) = zsn
      return

   end subroutine rand_coords

!-----------------------------------------------------------------------
#ifdef SHEAR
   subroutine sn_shear(cg, ysnoi)

      use constants,   only: ydim, LO
      use dataio_pub,  only: die
      use domain,      only: is_multicg, dom
      use grid_cont,   only: grid_container
      use shear,       only: delj, eps

      implicit none

      type(grid_container), pointer     :: cg
      real, dimension(3), intent(inout) :: ysnoi
      integer                           :: jsn, jremap
      real                              :: dysn, epsi, epso

      if (is_multicg) call die("[snsources:sn_shear] multiple grid pieces per procesor not implemented yet") !nontrivial SHEAR

      jsn  = cg%js+int((ysnoi(2)-dom%edge(ydim, LO))*cg%idy)
      dysn  = mod(ysnoi(2), cg%dy)

      epsi   = eps*cg%dy
      epso   = -epsi

!  outer boundary
      jremap = jsn - delj
      jremap = mod(mod(jremap, cg%nyb)+cg%nyb, cg%nyb)
      if (jremap <= (cg%js-1)) jremap = jremap + cg%nyb

      ysnoi(1) = cg%y(jremap) + epso + dysn

!  inner boundary
      jremap = jsn + delj
      jremap = mod(jremap, cg%nyb)+cg%nyb
      if (jremap >= (cg%je+1)) jremap = jremap - cg%nyb

      ysnoi(3) = cg%y(jremap) + epsi + dysn

   end subroutine sn_shear
#endif /* SHEAR */
!-----------------------------------------------------------------------

!>
!! \brief Function that generates values of normal distribution (from Numerical Recipies)
!! \return x random value of uniform distribution
!! \return y random value of uniform distribution
!! \return @e real, random value of normal distribution
!! \todo Change function @c gasdev into another one
!<
      function gasdev(x,y)

         implicit none

         real               :: x, y, x1, y1, r, fac, gasdev
         real,    save      :: gset
         integer, save      :: iset, irand
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
                  irand = irand+2
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

end module snsources
