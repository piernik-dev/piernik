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

module initproblem

   implicit none

   private
   public :: read_problem_par, init_prob

   real               :: t_sn
   integer            :: n_sn
   integer            :: norm_step
   real               :: d0, p0, bx0, by0, bz0, x0, y0, z0, r0, beta_cr, amp_cr

   namelist /PROBLEM_CONTROL/  d0, p0, bx0, by0, bz0, &
                               x0, y0, z0, r0, &
                               beta_cr, amp_cr, &
                               norm_step

contains

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use dataio_pub,    only: ierrh, par_file, namelist_errh, compare_namelist, cmdl_nml ! QA_WARN required for diff_nml
      use dataio_pub,    only: die
      use grid,          only: cg
      use mpi,           only: MPI_INTEGER, MPI_DOUBLE_PRECISION
      use mpisetup,      only: ibuff, rbuff, buffer_dim, comm, ierr, master, slave

      implicit none

      t_sn = 0.0

      d0           = 1.0e5     !< density
      p0           = 1.0       !< pressure
      bx0          =   0.      !< Magnetic field component x
      by0          =   0.      !< Magnetic field component y
      bz0          =   0.      !< Magnetic field component z
      x0           = 0.0       !< x-position of the blob
      y0           = 0.0       !< y-position of the blob
      z0           = 0.0       !< z-position of the blob
      r0           = 5.*cg%dxmn   !< radius of the blob

      beta_cr      = 0.0       !< ambient level
      amp_cr       = 1.0       !< amplitude of the blob

      norm_step    = 10        !< how often to compute the norm (in steps)

      if (master) then

         diff_nml(PROBLEM_CONTROL)

         rbuff(1) = d0
         rbuff(2) = p0
         rbuff(3) = bx0
         rbuff(4) = by0
         rbuff(5) = bz0
         rbuff(6) = x0
         rbuff(7) = y0
         rbuff(8) = z0
         rbuff(9) = r0
         rbuff(10)= beta_cr
         rbuff(11)= amp_cr

         ibuff(1) = norm_step

      endif

      call MPI_Bcast(ibuff, buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_Bcast(rbuff, buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

      if (slave) then

         d0           = rbuff(1)
         p0           = rbuff(2)
         bx0          = rbuff(3)
         by0          = rbuff(4)
         bz0          = rbuff(5)
         x0           = rbuff(6)
         y0           = rbuff(7)
         z0           = rbuff(8)
         r0           = rbuff(9)
         beta_cr      = rbuff(10)
         amp_cr       = rbuff(11)

         norm_step    = ibuff(1)

      endif

      if (r0 == 0.) call die("[initproblem:read_problem_par] r0 == 0")

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine init_prob

      use arrays,         only: b, u
      use dataio_pub,     only: msg, warn, printinfo
      use fluidindex,     only: ibx, iby, ibz, flind
      use grid,           only: cg
      use initcosmicrays, only: iarr_crn, iarr_crs, gamma_crn, K_crn_paral, K_crn_perp
      use initionized,    only: idni, imxi, imzi, ieni, gamma_ion
      use mpisetup,       only: comm3d, ierr, master, xdim, ydim, zdim, has_dir, dom
      use mpi,            only: MPI_IN_PLACE, MPI_INTEGER, MPI_MAX
#ifdef COSM_RAYS_SOURCES
      use crcomposition,  only: icr_H1, icr_C12
#endif /* COSM_RAYS_SOURCES */

      implicit none

      integer :: i, j, k
      integer :: icr
      real    :: cs_iso
      real    :: xsn, ysn, zsn
      real    :: r2, maxv
      integer :: ipm, jpm, kpm

#ifndef COSM_RAYS_SOURCES
      integer, parameter :: icr_H1 = 1, icr_C12 = 2
#endif /* !COSM_RAYS_SOURCES */

      ! BEWARE: temporary fix
      xsn = x0
      ysn = y0
      zsn = z0

! Uniform equilibrium state

      cs_iso = sqrt(p0/d0)

      if (.not.has_dir(xdim)) bx0 = 0. ! ignore B field in nonexistent direction
      if (.not.has_dir(ydim)) by0 = 0.
      if (.not.has_dir(zdim)) bz0 = 0.

      if ((bx0**2 + by0**2 + bz0**2 == 0.) .and. (any(K_crn_paral(:) /= 0.) .or. any(K_crn_perp(:) /= 0.))) then
         call warn("[initproblem:init_prob] No magnetic field is set, K_crn_* also have to be 0.")
         K_crn_paral(:) = 0.
         K_crn_perp(:)  = 0.
      endif

      b(ibx, :, :, :) = bx0
      b(iby, :, :, :) = by0
      b(ibz, :, :, :) = bz0
      u(idni, :, :, :) = d0
      u(imxi:imzi, :, :, :) = 0.0

#ifndef ISO
      do k = 1, cg%nz
         do j = 1, cg%ny
            do i = 1, cg%nx
               u(ieni,i,j,k) = p0/(gamma_ion-1.0) + &
                    &          0.5*sum(u(imxi:imzi,i,j,k)**2,1)/u(idni,i,j,k) + &
                    &          0.5*sum(b(:,i,j,k)**2,1)
            enddo
         enddo
      enddo
#endif /* !ISO */

#ifdef COSM_RAYS
      do icr = 1, flind%crs%all
         u(iarr_crs(icr), :, :, :) =  beta_cr*cs_iso**2 * u(idni, :, :, :)/(gamma_crn(icr)-1.0)
      enddo

! Explosions
      do icr = 1, flind%crn%all
         do k = cg%ks, cg%ke
            do j = cg%js, cg%je
               do i = cg%is, cg%ie

                  do ipm=-1,1
                     do jpm=-1,1
                        do kpm=-1,1

                           r2 = (cg%x(i)-xsn+real(ipm)*dom%Lx)**2+(cg%y(j)-ysn+real(jpm)*dom%Ly)**2+(cg%z(k)-zsn+real(kpm)*dom%Lz)**2
                           if (icr == icr_H1) then
                              u(iarr_crn(icr), i, j, k) = u(iarr_crn(icr), i, j, k) + amp_cr*exp(-r2/r0**2)
                           elseif (icr == icr_C12) then
                              u(iarr_crn(icr), i, j, k) = u(iarr_crn(icr), i, j, k) + amp_cr*0.1*exp(-r2/r0**2) ! BEWARE: magic number
                           else
                              u(iarr_crn(icr), i, j, k) = 0.0
                           endif

                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo

      do icr = 1, flind%crs%all
         maxv = maxval(u(iarr_crs(icr),:,:,:))
         call MPI_Allreduce(MPI_IN_PLACE, maxv, 1, MPI_INTEGER, MPI_MAX, comm3d, ierr)
         if (master) then
            write(msg,*) '[initproblem:init_prob] icr=',icr,' maxecr =',maxv
            call printinfo(msg)
         endif
      enddo
#endif /* COSM_RAYS */

   end subroutine init_prob

end module initproblem
