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

! Initial condition for Keplerian disk
! Written by: M. Hanasz, March 2006

   use constants,    only: cbuff_len

   implicit none

   private
   public  :: read_problem_par, init_prob, problem_pointers

   real    :: sigma0, Rin, R0, HtoR, eps, amp
   character(len=cbuff_len) :: sigma_model

   namelist /PROBLEM_CONTROL/  sigma0, amp, Rin, R0, HtoR, sigma_model, eps

contains

!-----------------------------------------------------------------------------

   subroutine problem_pointers

      implicit none

   end subroutine problem_pointers

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use dataio_pub,    only: par_file, ierrh, namelist_errh, compare_namelist, cmdl_nml         ! QA_WARN required for diff_nml
      use constants,     only: cbuff_len
      use mpisetup,      only: cbuff, rbuff, master, slave, comm, ierr, buffer_dim
      use mpi,           only: MPI_CHARACTER, MPI_DOUBLE_PRECISION

      implicit none

      Sigma0  = 1.0
      Rin     = 1.0e-4
      R0      = 1.0
      HtoR    = 1.0
      sigma_model = 'hayashi'
      eps     = 1.0
      amp     = 1.e-5

      if (master) then

         diff_nml(PROBLEM_CONTROL)

         cbuff(1) =  sigma_model

         rbuff(1) = sigma0
         rbuff(2) = Rin
         rbuff(3) = R0
         rbuff(4) = HtoR
         rbuff(5) = eps
         rbuff(6) = amp

      endif

      call MPI_Bcast(cbuff, cbuff_len*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_Bcast(rbuff,           buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

      if (slave) then

         sigma_model  = cbuff(1)

         sigma0       = rbuff(1)
         Rin          = rbuff(2)
         R0           = rbuff(3)
         HtoR         = rbuff(4)
         eps          = rbuff(5)
         amp          = rbuff(6)

      endif

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   real function dens_Rdistr(R,Rin,n)

      implicit none

      real, intent(in) :: R,Rin,n
      real :: ninv

      ninv = 1./n
      dens_Rdistr = (R - Rin)**ninv / R**(2.+ninv)

   end function dens_Rdistr

   subroutine init_prob

      use constants,   only: pi, dpi
      use domain,      only: dom
      use global,      only: smalld
      use gravity,     only: ptmass
      use grid,        only: cga
      use grid_cont,   only: cg_list_element, grid_container
      use initneutral, only: idnn, imxn, imyn, imzn, cs_iso_neu, cs_iso_neu2
      use initdust,    only: idnd, imxd, imyd, imzd
      use units,       only: newtong
#ifndef ISO
      use initneutral, only: ienn, gamma_neu
      use global,      only: smallei
#endif /* !ISO */

      implicit none

      integer :: i,j,k
      real :: xi,yj,zk, rc, H0, sqr_gm, rho0
      real :: n,norm, H, ninv
      real :: gradP, iOmega, ilook, gradgp
      real, dimension(:, :), allocatable :: noise
      real, dimension(:), allocatable :: omega,omegad
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg
#ifndef ISO
      real :: vx, vy, vz
#endif /* !ISO */

      cgl => cga%cg_leafs%cg_l(1)
      do while (associated(cgl))
         cg => cgl%cg

         allocate(omega(cg%nx),omegad(cg%nx), noise(3, cg%nz))
!   Secondary parameters

         sqr_gm = sqrt(newtong*ptmass)

         n = 0.5*Rin / (R0 - Rin)
         ninv = 1./n
         if (sigma_model == 'hayashi') then
            sigma0 = 0.2* R0**(-1.5)
         endif
         H0 = R0 * HtoR
         cs_iso_neu2 = H0 * (pi*newtong) * sigma0
         cs_iso_neu = sqrt(cs_iso_neu2)

         rho0 = sigma0 / (sqrt(dpi)*H0)

         norm = 1. / dens_Rdistr(R0,Rin,n)

         do k = 1, cg%nz
            zk = cg%z(k)
            do j = 1, cg%ny
               yj = cg%y(j)
               do i = 1, cg%nx
                  xi = cg%x(i)
                  rc = sqrt(xi**2+yj**2)
                  H = HtoR * rc
                  cg%u%arr(idnn,i,j,k) = max(rho0 * norm * dens_RdistR(rc,Rin,n) * exp(- 0.25 * zk**2 / H**2 ), smalld)
                  cg%u%arr(idnd,i,j,k) = eps*cg%u%arr(idnn,i,j,k)
               enddo
            enddo
         enddo

         do i = 2, cg%nx-1   ! 2d
            rc= cg%x(i)*sqrt(2.0)
            gradgp=  0.5*(cg%gp%arr(i+1,i+1,max(cg%nz/2,1))-cg%gp%arr(i-1,i-1,max(cg%nz/2,1)))/cg%dx/sqrt(2.)
            gradp = -0.5*(cg%u%arr(idnn,i+1,i+1,max(cg%nz/2,1))-cg%u%arr(idnn,i-1,i-1,max(cg%nz/2,1)))/cg%dx /sqrt(2.)*cs_iso_neu2
            omega(i)  = sqrt( abs( (gradgp-gradp)/rc ) )
            omegad(i) = sqrt( abs(    gradgp/rc      ) )
         enddo
         omega(1)  = omega(2);  omega(cg%nx)  = omega(cg%nx-1)
         omegad(1) = omegad(2); omegad(cg%nx) = omegad(cg%nx-1)

         call random_seed()

         do j = 1, cg%ny
            yj = cg%y(j)
            do i = 1, cg%nx
               xi = cg%x(i)
               rc = sqrt(xi**2+yj**2)
               call random_number(noise)

               ilook = (rc-dom%xmin)/cg%dx/sqrt(2.) + 0.5 + cg%nb
               iOmega = omega(int(ilook))+(rc-cg%x(int(ilook))*sqrt(2.))*(omega(int(ilook)+1)-omega(int(ilook))) &
                    &   / (cg%x(int(ilook)+1)-cg%x(int(ilook)))/sqrt(2.)
!
!
               cg%u%arr(imxn,i,j,:) = -yj*iOmega*cg%u%arr(idnn,i,j,:)
               cg%u%arr(imyn,i,j,:) =  xi*iOmega*cg%u%arr(idnn,i,j,:)
               cg%u%arr(imzn,i,j,:) = 0.0
#ifndef ISO
               cg%u%arr(ienn,i,j,:) = cs_iso_neu2/(gamma_neu-1.0)*cg%u%arr(idnn,i,j,:)
               cg%u%arr(ienn,i,j,:) = max(cg%u%arr(ienn,i,j,:), smallei)
               cg%u%arr(ienn,i,j,:) = cg%u%arr(ienn,i,j,:) +0.5*(vx**2+vy**2+vz**2)*cg%u%arr(idnn,i,j,:)
#endif /* !ISO */

               iOmega = omegad(int(ilook))+(rc-cg%x(int(ilook))*sqrt(2.))*(omegad(int(ilook)+1)-omegad(int(ilook))) &
                    &   /(cg%x(int(ilook)+1)-cg%x(int(ilook)))/sqrt(2.)
               cg%u%arr(imxd,i,j,:) = -yj*iOmega*cg%u%arr(idnd,i,j,:) + amp*(noise(1,:)-0.5)
               cg%u%arr(imyd,i,j,:) =  xi*iOmega*cg%u%arr(idnd,i,j,:) + amp*(noise(2,:)-0.5)
               cg%u%arr(imzd,i,j,:) = 0.0 + amp*(noise(3,:)-0.5)

            enddo
         enddo

         deallocate(omega, omegad, noise)

         cgl => cgl%nxt
      enddo

   end subroutine init_prob

end module initproblem

