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
#include "piernik.def"
#include "macros.h"

module initproblem

! Initial condition for the cosmic ray driven dynamo
! Based on Parker instability setup
! Written by: M. Hanasz, February 2006
! Modified by M.Hanasz for CR-driven dynamo

   use problem_pub, only: problem_name, run_id

   real :: d0, bxn,byn,bzn, h_sn, f_sn_kpc2, amp_cr, beta_cr, r_sn
   real :: ethu, f_sn, amp_ecr_sn, alpha, x0, y0, z0

   namelist /PROBLEM_CONTROL/  problem_name, run_id, &
                               d0, &
                               bxn,byn,bzn, &
                               x0, y0, z0, alpha, &
                               r_sn, h_sn, f_sn_kpc2, amp_cr, beta_cr
   contains

!-----------------------------------------------------------------------------

   subroutine read_problem_par
      use dataio_public, only: ierrh, msg, par_file
      use errh,          only: namelist_errh
      use func,          only: compare_namelist
      use mpisetup,      only: cbuff_len, cbuff, rbuff, buffer_dim, comm, ierr, proc, &
                               MPI_CHARACTER, MPI_DOUBLE_PRECISION

      implicit none

      problem_name = 'xxx'
      run_id  = 'aaa'
      d0     = 1.0
      bxn    = 0.0
      byn    = 1.0
      bzn    = 0.0
      x0     = 0.0
      y0     = 0.0
      z0     = 0.0

      if(proc == 0) then

         diff_nml(PROBLEM_CONTROL)

         cbuff(1) =  problem_name
         cbuff(2) =  run_id

         rbuff(1)  = d0
         rbuff(2)  = bxn
         rbuff(3)  = byn
         rbuff(4)  = bzn
         rbuff(5)  = x0
         rbuff(6)  = y0
         rbuff(7)  = z0
         rbuff(8)  = amp_cr
         rbuff(9)  = beta_cr
         rbuff(10) = r_sn
         rbuff(11) = h_sn
         rbuff(12) = f_sn_kpc2
         rbuff(13) = alpha

      endif

      call MPI_Bcast(cbuff, cbuff_len*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_Bcast(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

      if (proc /= 0) then

         problem_name = cbuff(1)
         run_id       = cbuff(2)(1:3)

         d0           = rbuff(1)
         bxn          = rbuff(2)
         byn          = rbuff(3)
         bzn          = rbuff(4)
         x0           = rbuff(5)
         y0           = rbuff(6)
         z0           = rbuff(7)
         amp_cr       = rbuff(8)
         beta_cr      = rbuff(9)
         r_sn         = rbuff(10)
         h_sn         = rbuff(11)
         f_sn_kpc2    = rbuff(12)
         alpha        = rbuff(13)

      endif

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine init_prob
      use arrays,         only: u, b, dprof
      use fluidindex,     only: ibx,iby,ibz
      use grid,           only: x, y, z, nx, ny, nz, nxd, nyd, xmin, xmax, ymin, ymax, lx, ly
      use hydrostatic,    only: hydrostatic_zeq
      use initcosmicrays, only: gamma_crs, iarr_crs, cr_eff
      use initfluids,     only: cs_iso2
      use initionized,    only: idni, imxi, imyi, imzi
      use mpisetup,       only: proc, smalld
#ifdef SHEAR
      use shear,          only: qshear, omega
#endif /* SHEAR */

      implicit none

      integer :: i,j,k
      real :: b0, csim2

      ethu = 7.0**2/(5.0/3.0-1.0) * 1.0    ! thermal energy unit=0.76eV/cm**3
                                           ! for c_si= 7km/s, n=1/cm^3
                                           ! gamma=5/3

      amp_ecr_sn = 4.96e6*cr_eff/r_sn**3   ! cosmic ray explosion amplitude
                                           ! in units:
                                           ! e_0 = 1/(5/3-1)*rho_0*c_s0**2
                                           ! rho_0=1.67e-24g/cm**3,
                                           ! c_s0 = 7km/s

      f_sn = f_sn_kpc2                     ! SN freq

      if(nxd /=1) then
         f_sn = f_sn * (xmax-xmin)/1000.0
      else
         f_sn = f_sn * 2.0*r_sn/1000.0
      endif

      if(nyd /=1) then
         f_sn = f_sn * (ymax-ymin)/1000.0
      else
         f_sn = f_sn * 2.0*r_sn/1000.0
      endif

      write(*,*) proc, csim2
!   Secondary parameters

      b0 = sqrt(2.*alpha*d0*cs_iso2)

      csim2 = cs_iso2*(1.0+alpha)

      call hydrostatic_zeq(1, 1, d0, csim2, dprof)

      do k = 1,nz
         do j = 1,ny
            do i = 1,nx
               u(idni,i,j,k)   = max(smalld,dprof(k))

               u(imxi,i,j,k) = 0.0
               u(imyi,i,j,k) = 0.0
               u(imzi,i,j,k) = 0.0
#ifdef SHEAR
               u(imyi,i,j,k) = -qshear*omega*x(i)*u(idni,i,j,k)
#endif /* SHEAR */

#ifndef ISO
               u(ieni,i,j,k)   = cs_iso2/(gamma_ion-1.0) * u(idni,i,j,k) &
                               + 0.5*(u(imxi,i,j,k)**2 + u(imyi,i,j,k)**2 + &
                                      u(imzi,i,j,k)**2 ) / u(idni,i,j,k)
#endif /* !ISO */
#ifdef COSM_RAYS
               u(iarr_crs,i,j,k)   =  beta_cr*cs_iso2 * u(idni,i,j,k)/( gamma_crs - 1.0 )
#ifdef GALAXY
! Single SN explosion in x0,y0,z0 at t = 0 if amp_cr /= 0
               u(iarr_crs,i,j,k)= u(iarr_crs,i,j,k) &
                     + amp_cr*exp(-((x(i)- x0    )**2 + (y(j)- y0    )**2 + (z(k)-z0)**2)/r_sn**2) &
                     + amp_cr*exp(-((x(i)-(x0+Lx))**2 + (y(j)- y0    )**2 + (z(k)-z0)**2)/r_sn**2) &
                     + amp_cr*exp(-((x(i)- x0    )**2 + (y(j)-(y0+Ly))**2 + (z(k)-z0)**2)/r_sn**2) &
                     + amp_cr*exp(-((x(i)-(x0+Lx))**2 + (y(j)-(y0+Ly))**2 + (z(k)-z0)**2)/r_sn**2)
#endif /* GALAXY */
#endif /* COSM_RAYS */
            enddo
         enddo
      enddo

      do k = 1,nz
         do j = 1,ny
            do i = 1,nx
               b(ibx,i,j,k)   = b0*sqrt(u(idni,i,j,k)/d0)* bxn/sqrt(bxn**2+byn**2+bzn**2)
               b(iby,i,j,k)   = b0*sqrt(u(idni,i,j,k)/d0)* byn/sqrt(bxn**2+byn**2+bzn**2)
               b(ibz,i,j,k)   = b0*sqrt(u(idni,i,j,k)/d0)* bzn/sqrt(bxn**2+byn**2+bzn**2)
#ifndef ISO
               u(ieni,i,j,k)   = u(ieni,i,j,k) +0.5*sum(b(:,i,j,k)**2,1)
#endif /* !ISO */
            enddo
         enddo
      enddo

      return
   end subroutine init_prob

end module initproblem

