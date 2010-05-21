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

   integer            :: norm_step
   real               :: d0, p0, bx0, by0, bz0, x0, y0, z0, r0, beta_cr, amp_cr

   namelist /PROBLEM_CONTROL/  problem_name, run_id,      &
                               d0, p0, bx0, by0, bz0, &
                               x0, y0, z0, r0, &
                               beta_cr, amp_cr, &
                               norm_step

   contains

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use mpisetup, only: MPI_CHARACTER, MPI_INTEGER, MPI_DOUBLE_PRECISION, &
           &              cbuff, ibuff, rbuff, buffer_dim, comm, ierr, proc, cwd
      use grid,     only: dxmn
      use errh,     only: namelist_errh

      implicit none

      integer            :: ierrh
      character(LEN=100) :: par_file, tmp_file

      par_file = trim(cwd)//'/problem.par'
      tmp_file = trim(cwd)//'/tmp.log'

      problem_name = 'aaa'
      run_id       = 'aaa'
      d0           = 1.0       !< density
      p0           = 1.0       !< pressure
      bx0          =   0.      !< Magnetic field component x
      by0          =   0.      !< Magnetic field component y
      bz0          =   0.      !< Magnetic field component z
      x0           = 0.0       !< x-position of the blob
      y0           = 0.0       !< y-position of the blob
      z0           = 0.0       !< z-position of the blob
      r0           = dxmn/2.   !< radius of the blob

      beta_cr      = 0.0       !< ambient level
      amp_cr       = 1.0       !< amplitude of the blob

      norm_step    = 10        !< how often to compute the norm (in steps)

      if(proc .eq. 0) then
         open(1,file=par_file)
         read(unit=1,nml=PROBLEM_CONTROL,iostat=ierrh)
         call namelist_errh(ierrh,'PROBLEM_CONTROL')
         write(*,nml=PROBLEM_CONTROL)
         close(1)

         open(3, file=tmp_file, position='append')
         write(3,nml=PROBLEM_CONTROL)
         write(3,*)
         close(3)
      endif

      if (proc == 0) then

         cbuff(1) = problem_name
         cbuff(2) = run_id

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

         norm_step    = ibuff(1)

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
      use types,          only : problem_customize_solution, finalize_problem

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


      b(ibx, 1:nx, 1:ny, 1:nz) = bx0
      b(iby, 1:nx, 1:ny, 1:nz) = by0
      b(ibz, 1:nx, 1:ny, 1:nz) = bz0
      u(idni, 1:nx, 1:ny, 1:nz) = d0
      u(imxi:imzi, 1:nx, 1:ny, 1:nz) = 0.0

#ifndef ISO
      do k = 1,nz
         do j = 1,ny
            do i = 1,nx
               u(ieni,i,j,k)      = p0/(gamma_ion-1.0) + &
                    &               0.5*sum(u(imxi:imzi,i,j,k)**2,1)/u(idni,i,j,k) + &
                    &               0.5*sum(b(:,i,j,k)**2,1)
            enddo
         enddo
      enddo
#endif /* ISO */

! Explosions

#ifdef COSM_RAYS
      u(iecr, 1:nx, 1:ny, 1:nz)      =  beta_cr*cs_iso**2 * u(idni, 1:nx, 1:ny, 1:nz)/(gamma_crs(icr)-1.0)
      do k = ks,ke
         do j = nb+1,ny-nb
            do i = nb+1,nx-nb
               r2 = (x(i)-x0)**2+(y(j)-y0)**2+(z(k)-z0)**2
               u(iecr,i,j,k)= u(iecr,i,j,k) + amp_cr*exp(-r2/r0**2)
            enddo
         enddo
      enddo

      write(*,*) 'maxecr =',maxval(u(iecr,:,:,:)), amp_cr

#endif /* COSM_RAYS */

      problem_customize_solution => check_norm
      finalize_problem           => check_norm

      return
   end subroutine init_prob

!-----------------------------------------------------------------------------

   subroutine check_norm

      use mpisetup,       only : proc, comm3d, ierr, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_MIN, MPI_MAX, MPI_IN_PLACE, t, nstep
      use grid,           only : x, y, z, is, ie, js, je, ks, ke
      use arrays,         only : b, u
      use initcosmicrays, only : iarr_crs, ncrn, ncre, K_crn_paral, K_crn_perp
      use errh,           only : die
      use dataio_public,  only : code_progress, PIERNIK_FINISHED, halfstep

      implicit none

      integer :: i, j, k
      real    :: r_par2, r_perp2, delx, dely, delz, magb2, ampt, r0_par2, r0_perp2, crt
      integer :: iecr = -1
      integer, parameter :: icr = 1 !< Only first CR component
      real, dimension(2) :: norm, dev

      integer, save :: nn = 0

      if (code_progress < PIERNIK_FINISHED .and. (mod(nstep, norm_step) /=0 .or. halfstep)) return

      if (ncrn+ncre >= icr) then
         iecr = iarr_crs(icr)
      else
         call die("[initproblem:init_prob] No CR components defined.")
      end if

      norm(:) = 0.
      dev(1) = huge(1.0)
      dev(2) = -dev(1)

      magb2 = bx0**2 + by0**2 + bz0**2

      r0_par2  = r0**2 + 4 * (K_crn_paral(icr) + K_crn_perp(icr)) * t
      r0_perp2 = r0**2 + 4 * K_crn_perp(icr) * t
      ampt     = amp_cr * r0**2 / sqrt(r0_par2 * r0_perp2)

      !write(*,'(a,5g15.6)')"[i:cn] tBrrA:",t,magb2,r0_par2,r0_perp2,ampt

      do k = ks, ke
         delz = z(k) - z0
         do j = js, je
            dely = y(j) - y0
            do i = is, ie
               delx = x(i) - x0

               r_par2 = (bx0*delx + by0*dely + bz0*delz)**2/magb2 ! square of the distance form the center of the bump in direction parallel to the magnetic field
               r_perp2 = delx**2 + dely**2 + delz**2 - r_par2
               crt = ampt * exp( - r_par2/r0_par2 - r_perp2/r0_perp2)

               !if (nn==400) write(*,'(a,2i4,a,2f10.3,a,2f10.3,a,2f17.10)')"icn ",i,j," Dxy= ",delx,dely," r= ",sqrt(r_par2),sqrt(r_perp2)," crt= ",crt, u(iecr, i, j, k)

               norm(1) = norm(1) + (crt - u(iecr, i, j, k))**2
               norm(2) = norm(2) + crt**2
               dev(1) = min(dev(1), (crt - u(iecr, i, j, k)))
               dev(2) = max(dev(2), (crt - u(iecr, i, j, k)))

            end do
         end do
      end do

      call MPI_Allreduce(MPI_IN_PLACE, norm,   2, MPI_DOUBLE_PRECISION, MPI_SUM, comm3d, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, dev(1), 1, MPI_DOUBLE_PRECISION, MPI_MIN, comm3d, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, dev(2), 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm3d, ierr)

      if (proc == 0) write(*,'(a,f12.5,a,2f12.5)')"[initproblem:check_norm] L2 error norm = ", sqrt(norm(1)/norm(2)), " min and max error = ", dev(1:2)

      nn = nn + 1

   end subroutine check_norm

end module initproblem
