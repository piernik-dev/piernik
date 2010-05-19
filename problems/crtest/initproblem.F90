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

module initproblem

   use problem_pub, only: problem_name, run_id

   real    :: t_sn
   integer :: n_sn, ierrh
   real    :: d0, p0, bx0, by0, bz0, x0, y0, z0, r0, beta_cr, amp_cr

   namelist /PROBLEM_CONTROL/  problem_name, run_id,      &
                               d0, p0, bx0, by0, bz0, &
                               x0, y0, z0, r0, &
                               beta_cr, amp_cr

   contains

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use mpisetup, only: MPI_CHARACTER, MPI_INTEGER, MPI_DOUBLE_PRECISION, &
           &              cbuff, ibuff, rbuff, buffer_dim, comm, ierr, proc
      use grid,     only: dxmn
      use errh,     only: namelist_errh

      implicit none

      t_sn = 0.0

      problem_name = 'aaa'
      run_id  = 'aaa'
      d0      = 1.0
      p0      = 1.0
      bx0     =   0.
      by0     =   0.
      bz0     =   0.
      x0      = 0.0
      y0      = 0.0
      z0      = 0.0
      r0      = dxmn/2.

      beta_cr    = 0.0
      amp_cr     = 1.0

      if(proc .eq. 0) then
         open(1,file='problem.par')
         read(unit=1,nml=PROBLEM_CONTROL,iostat=ierrh)
         call namelist_errh(ierrh,'PROBLEM_CONTROL')
         write(*,nml=PROBLEM_CONTROL)
         close(1)

         open(3, file='tmp.log', position='append')
         write(3,nml=PROBLEM_CONTROL)
         write(3,*)
         close(3)
      endif


      if (proc == 0) then

         cbuff(1) =  problem_name
         cbuff(2) =  run_id

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

      end if

      call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

      if (proc /= 0) then

         problem_name = cbuff(1)
         run_id       = cbuff(2)(1:3)

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

      endif

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine init_prob

      use fluidindex,     only : ibx,iby,ibz
      use initionized,    only : idni,imxi,imyi,imzi,ieni
      use initcosmicrays, only : gamma_crs, iarr_crs, ncrn, ncre
      use initionized,    only : gamma_ion
      use arrays,         only : b, u
      use grid,           only : nx, ny, nz, nb, ks, ke, x, y, z
      use errh,           only : die

      implicit none

      integer :: i, j, k
      real    :: cs_iso, r2
      integer :: iecr = -1
      integer, parameter :: icr = 1 !< Only first CR component

      if (ncrn+ncre >= icr) then
         iecr = iarr_crs(icr)
      else
         call die("[initproblem:init_prob] No CR components defined.")
      end if

! Uniform equilibrium state

      cs_iso = sqrt(p0/d0)

      do k = 1,nz
         do j = 1,ny
            do i = 1,nx

               b(ibx,i,j,k)       = bx0
               b(iby,i,j,k)       = by0
               b(ibz,i,j,k)       = bz0

               u(idni,i,j,k)      = d0
               u(imxi:imzi,i,j,k) = 0.0
#ifndef ISO
               u(ieni,i,j,k)      = p0/(gamma_ion-1.0)
               u(ieni,i,j,k)      = u(ieni,i,j,k) &
                           + 0.5*sum(u(imxi:imzi,i,j,k)**2,1)/u(idni,i,j,k)
               u(ieni,i,j,k)      = u(ieni,i,j,k) + 0.5*sum(b(:,i,j,k)**2,1)
#endif /* ISO */

#ifdef COSM_RAYS
               u(iecr,i,j,k)      =  beta_cr*cs_iso**2 * u(idni,i,j,k)/(gamma_crs(icr)-1.0)
#endif /* COSM_RAYS */

            enddo
         enddo
      enddo

! Explosions

#ifdef COSM_RAYS
      do k = ks,ke
         do j = nb+1,ny-nb
            do i = nb+1,nx-nb
               r2 = (x(i)-x0)**2+(y(j)-y0)**2+(z(k)-z0)**2
               u(iecr,i,j,k)= u(iecr,i,j,k) + amp_cr*exp(-r2/r0**2)
            enddo
         enddo
      enddo
#endif /* COSM_RAYS */

      write(*,*) 'maxecr =',maxval(u(iecr,:,:,:))
      write(*,*) amp_cr

      return
   end subroutine init_prob

end module initproblem
