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

module initproblem

   implicit none

   private
   public :: read_problem_par, problem_initial_conditions, problem_pointers

   real :: d0, p0, bx0, by0, bz0, x0, y0, z0, r0, beta_cr, amp_cr1, amp_cr2

   namelist /PROBLEM_CONTROL/ d0, p0, bx0, by0, bz0, x0, y0, z0, r0, beta_cr, amp_cr1, amp_cr2

contains

!-----------------------------------------------------------------------------

   subroutine problem_pointers

      implicit none

   end subroutine problem_pointers

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use bcast,      only: piernik_MPI_Bcast
      use dataio_pub, only: die, nh
      use domain,     only: dom
      use func,       only: operator(.equals.)
      use mpisetup,   only: rbuff, master, slave

      implicit none

      d0             = 1.0e5       !< density
      p0             = 1.0         !< pressure
      bx0            =   0.        !< Magnetic field component x
      by0            =   0.        !< Magnetic field component y
      bz0            =   0.        !< Magnetic field component z
      x0             = 0.0         !< x-position of the blob
      y0             = 0.0         !< y-position of the blob
      z0             = 0.0         !< z-position of the blob
      r0             = merge(0., 5.* minval(dom%L_(:)/dom%n_d(:), mask=dom%has_dir(:)), dom%eff_dim == 0)  !< radius of the blob

      beta_cr        = 0.0         !< ambient level
      amp_cr1        = 1.0         !< amplitude of the blob
      amp_cr2        = 0.1*amp_cr1 !< amplitude for the second species

      if (master) then

         if (.not.nh%initialized) call nh%init()
         open(newunit=nh%lun, file=nh%tmp1, status="unknown")
         write(nh%lun,nml=PROBLEM_CONTROL)
         close(nh%lun)
         open(newunit=nh%lun, file=nh%par_file)
         nh%errstr=""
         read(unit=nh%lun, nml=PROBLEM_CONTROL, iostat=nh%ierrh, iomsg=nh%errstr)
         close(nh%lun)
         call nh%namelist_errh(nh%ierrh, "PROBLEM_CONTROL")
         read(nh%cmdl_nml,nml=PROBLEM_CONTROL, iostat=nh%ierrh)
         call nh%namelist_errh(nh%ierrh, "PROBLEM_CONTROL", .true.)
         open(newunit=nh%lun, file=nh%tmp2, status="unknown")
         write(nh%lun,nml=PROBLEM_CONTROL)
         close(nh%lun)
         call nh%compare_namelist()

         rbuff(1)  = d0
         rbuff(2)  = p0
         rbuff(3)  = bx0
         rbuff(4)  = by0
         rbuff(5)  = bz0
         rbuff(6)  = x0
         rbuff(7)  = y0
         rbuff(8)  = z0
         rbuff(9)  = r0
         rbuff(10) = beta_cr
         rbuff(11) = amp_cr1
         rbuff(12) = amp_cr2

      endif

      call piernik_MPI_Bcast(rbuff)

      if (slave) then

         d0        = rbuff(1)
         p0        = rbuff(2)
         bx0       = rbuff(3)
         by0       = rbuff(4)
         bz0       = rbuff(5)
         x0        = rbuff(6)
         y0        = rbuff(7)
         z0        = rbuff(8)
         r0        = rbuff(9)
         beta_cr   = rbuff(10)
         amp_cr1   = rbuff(11)
         amp_cr2   = rbuff(12)

      endif

      if (r0 .equals. 0.) call die("[initproblem:read_problem_par] r0 == 0")

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine problem_initial_conditions

      use cg_leaves,      only: leaves
      use cg_list,        only: cg_list_element
      use constants,      only: xdim, ydim, zdim
      use domain,         only: dom
      use fluidindex,     only: flind
      use fluidtypes,     only: component_fluid
      use func,           only: ekin, emag
      use grid_cont,      only: grid_container
#ifdef COSM_RAYS
      use allreduce,      only: piernik_MPI_Allreduce
      use constants,      only: ndims, LO, HI, pMAX, BND_PER
      use cr_data,        only: eCRSP, icr_H1, icr_C12, cr_index
      use dataio_pub,     only: msg, warn, printinfo
      use func,           only: operator(.equals.), operator(.notequals.)
      use initcosmicrays, only: iarr_crn, iarr_crs, gamma_cr_1, K_cr_paral, K_cr_perp
      use mpisetup,       only: master
#endif /* COSM_RAYS */

      implicit none

      class(component_fluid), pointer :: fl
#ifdef COSM_RAYS
      integer                         :: i, j, k, icr, ipm, jpm, kpm
      integer, dimension(ndims,LO:HI) :: mantle
      real                            :: decr, r2, maxv
#endif /* COSM_RAYS */
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg

      fl => flind%ion

! Uniform equilibrium state

      !cs_iso = sqrt(p0/d0)

      if (.not.dom%has_dir(xdim)) bx0 = 0. ! ignore B field in nonexistent direction
      if (.not.dom%has_dir(ydim)) by0 = 0.
      if (.not.dom%has_dir(zdim)) bz0 = 0.

#ifdef COSM_RAYS
      if ((bx0**2 + by0**2 + bz0**2 .equals. 0.) .and. (any(K_cr_paral(:) .notequals. 0.) .or. any(K_cr_perp(:) .notequals. 0.))) then
         call warn("[initproblem:problem_initial_conditions] No magnetic field is set, K_cr_* also have to be 0.")
         K_cr_paral(:) = 0.
         K_cr_perp(:)  = 0.
      endif

      mantle = 0
      do i = xdim, zdim
         if (any(dom%bnd(i,:) == BND_PER)) mantle(i,:) = [-1,1] !> for periodic boundary conditions
      enddo
#endif /* COSM_RAYS */

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         call cg%set_constant_b_field([bx0, by0, bz0])
         cg%u(fl%idn,RNG) = d0
         cg%u(fl%imx:fl%imz,RNG) = 0.0

#ifndef ISO
         cg%u(fl%ien,RNG) = p0/fl%gam_1 + ekin(cg%u(fl%imx,RNG), cg%u(fl%imy,RNG), cg%u(fl%imz,RNG), cg%u(fl%idn,RNG)) + &
              &             emag(cg%b(xdim,RNG), cg%b(ydim,RNG), cg%b(zdim,RNG))
#endif /* !ISO */

#ifdef COSM_RAYS
         cg%u(iarr_crn, RNG) = 0.0
         if (eCRSP(icr_H1 )) cg%u(iarr_crn(cr_index(icr_H1 )), RNG) = beta_cr * fl%cs2 * cg%u(fl%idn, RNG) / gamma_cr_1
         if (eCRSP(icr_C12)) cg%u(iarr_crn(cr_index(icr_C12)), RNG) = beta_cr * fl%cs2 * cg%u(fl%idn, RNG) / gamma_cr_1

! Explosions
         do k = cg%ks, cg%ke
            do j = cg%js, cg%je
               do i = cg%is, cg%ie

                  decr = 0.0
                  do ipm = mantle(xdim,LO), mantle(xdim,HI)
                     do jpm = mantle(xdim,LO), mantle(xdim,HI)
                        do kpm = mantle(xdim,LO), mantle(xdim,HI)
                           r2 = (cg%x(i)-x0+real(ipm)*dom%L_(xdim))**2+(cg%y(j)-y0+real(jpm)*dom%L_(ydim))**2+(cg%z(k)-z0+real(kpm)*dom%L_(zdim))**2
                           if (r2/r0**2 < 0.9999*log(huge(1.))) &  ! preventing numerical underflow
                                decr = decr + exp(-r2/r0**2)
                        enddo
                     enddo
                  enddo
                  if (eCRSP(icr_H1 )) cg%u(iarr_crn(cr_index(icr_H1 )), i, j, k) = cg%u(iarr_crn(cr_index(icr_H1 )), i, j, k) + amp_cr1*decr
                  if (eCRSP(icr_C12)) cg%u(iarr_crn(cr_index(icr_C12)), i, j, k) = cg%u(iarr_crn(cr_index(icr_C12)), i, j, k) + amp_cr2*decr
               enddo
            enddo
         enddo
#endif /* COSM_RAYS */

         cgl => cgl%nxt
      enddo

#ifdef COSM_RAYS
      do icr = 1, flind%crs%all

         maxv = - huge(1.)
         cgl => leaves%first
         do while (associated(cgl))
            associate (cg => cgl%cg)
               maxv = max(maxv, maxval(cg%u(iarr_crs(icr), RNG)))
            end associate
            cgl => cgl%nxt
         enddo

         call piernik_MPI_Allreduce(maxv, pMAX)
         if (master) then
            write(msg,*) '[initproblem:problem_initial_conditions] icr=', icr, ' maxecr =', maxv
            call printinfo(msg)
         endif

      enddo
#endif /* COSM_RAYS */

   end subroutine problem_initial_conditions

end module initproblem
