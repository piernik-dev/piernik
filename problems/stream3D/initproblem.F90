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

      use dataio_pub,    only: par_file, ierrh, namelist_errh, compare_namelist, cmdl_nml, lun         ! QA_WARN required for diff_nml
      use constants,     only: cbuff_len
      use mpisetup,      only: cbuff, rbuff, master, slave, comm, mpi_err, buffer_dim, FIRST
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

      call MPI_Bcast(cbuff, cbuff_len*buffer_dim, MPI_CHARACTER,        FIRST, comm, mpi_err)
      call MPI_Bcast(rbuff,           buffer_dim, MPI_DOUBLE_PRECISION, FIRST, comm, mpi_err)

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
      dens_Rdistr = max((R - Rin),Rin)**ninv / R**(2.+ninv)

   end function dens_Rdistr

   subroutine init_prob

      use constants,   only: pi, dpi, xdim, ydim, zdim, LO
      use domain,      only: dom
      use global,      only: smalld
      use gravity,     only: ptmass
      use grid,        only: leaves
      use gc_list,     only: cg_list_element
      use grid_cont,   only: grid_container
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

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         allocate(omega(cg%n_(xdim)),omegad(cg%n_(xdim)), noise(3, cg%n_(zdim)))
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

         do k = 1, cg%n_(zdim)
            zk = cg%z(k)
            do j = 1, cg%n_(ydim)
               yj = cg%y(j)
               do i = 1, cg%n_(xdim)
                  xi = cg%x(i)
                  rc = sqrt(xi**2+yj**2)
                  H = HtoR * rc
                  cg%u(idnn,i,j,k) = max(rho0 * norm * dens_RdistR(rc,Rin,n) * exp(- 0.25 * zk**2 / H**2 ), smalld)
                  cg%u(idnd,i,j,k) = eps*cg%u(idnn,i,j,k)
               enddo
            enddo
         enddo

         do i = 2, cg%n_(xdim)-1   ! 2d
            rc= cg%x(i)*sqrt(2.0)
            gradgp=  0.5*(cg%gp(i+1,i+1,max(cg%n_(zdim)/2,1))-cg%gp(i-1,i-1,max(cg%n_(zdim)/2,1)))/cg%dx/sqrt(2.)
            gradp = -0.5*(cg%u(idnn,i+1,i+1,max(cg%n_(zdim)/2,1))-cg%u(idnn,i-1,i-1,max(cg%n_(zdim)/2,1)))/cg%dx /sqrt(2.)*cs_iso_neu2
            omega(i)  = sqrt( abs( (gradgp-gradp)/rc ) )
            omegad(i) = sqrt( abs(    gradgp/rc      ) )
         enddo
         omega(1)  = omega(2);  omega(cg%n_(xdim))  = omega(cg%n_(xdim)-1)
         omegad(1) = omegad(2); omegad(cg%n_(xdim)) = omegad(cg%n_(xdim)-1)

         call random_seed()

         do j = 1, cg%n_(ydim)
            yj = cg%y(j)
            do i = 1, cg%n_(xdim)
               xi = cg%x(i)
               rc = sqrt(xi**2+yj**2)
               call random_number(noise)

               ilook = (rc-dom%edge(xdim, LO))/cg%dx/sqrt(2.) + 0.5 + dom%nb
               iOmega = omega(int(ilook))+(rc-cg%x(int(ilook))*sqrt(2.))*(omega(int(ilook)+1)-omega(int(ilook))) &
                    &   / (cg%x(int(ilook)+1)-cg%x(int(ilook)))/sqrt(2.)
!
!
               cg%u(imxn,i,j,:) = -yj*iOmega*cg%u(idnn,i,j,:)
               cg%u(imyn,i,j,:) =  xi*iOmega*cg%u(idnn,i,j,:)
               cg%u(imzn,i,j,:) = 0.0
#ifndef ISO
               cg%u(ienn,i,j,:) = cs_iso_neu2/(gamma_neu-1.0)*cg%u(idnn,i,j,:)
               cg%u(ienn,i,j,:) = max(cg%u(ienn,i,j,:), smallei)
               cg%u(ienn,i,j,:) = cg%u(ienn,i,j,:) +0.5*(vx**2+vy**2+vz**2)*cg%u(idnn,i,j,:)
#endif /* !ISO */

               iOmega = omegad(int(ilook))+(rc-cg%x(int(ilook))*sqrt(2.))*(omegad(int(ilook)+1)-omegad(int(ilook))) &
                    &   /(cg%x(int(ilook)+1)-cg%x(int(ilook)))/sqrt(2.)
               cg%u(imxd,i,j,:) = -yj*iOmega*cg%u(idnd,i,j,:) + amp*(noise(1,:)-0.5)
               cg%u(imyd,i,j,:) =  xi*iOmega*cg%u(idnd,i,j,:) + amp*(noise(2,:)-0.5)
               cg%u(imzd,i,j,:) = 0.0 + amp*(noise(3,:)-0.5)

            enddo
         enddo

         deallocate(omega, omegad, noise)

         cgl => cgl%nxt
      enddo

   end subroutine init_prob

end module initproblem

